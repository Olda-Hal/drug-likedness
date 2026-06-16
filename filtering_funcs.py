import pathlib
import gemmi


DEFAULT_EXCLUDED = {
    # water / solvents
    "HOH", "WAT", "DOD", "H2O", "SOL",

    # common crystallization additives
    "GOL", "EDO", "DMS", "PEG", "ACT", "EDT", "MPD", "FMT",

    # ions
    "SO4", "PO4", "CL", "NA", "K", "MG", "CA", "ZN", "FE", "CU"
}


def heavy_atom_count(residue: gemmi.Residue) -> int:
    return sum(1 for atom in residue if not atom.is_hydrogen())


def is_organic_like(residue: gemmi.Residue) -> bool:
    """Returns True if residue contains at least one carbon atom."""
    return any(atom.element.name == "C" for atom in residue)


def ligand_score(residue: gemmi.Residue) -> float:
    """
    Heuristic scoring function for ligand selection.

    Higher score = more likely to be a drug-like small molecule.
    """

    ha = heavy_atom_count(residue)

    # filter out noise and large polymers
    if ha < 5:
        return -1.0
    if ha > 120:
        return -1.0

    # organic molecules are preferred
    organic_bonus = 1.5 if is_organic_like(residue) else 0.0

    # prefer drug-like size range (~10–60 heavy atoms)
    size_penalty = -abs(ha - 30) / 30.0

    return organic_bonus + size_penalty


def extract_best_ligand(pdb_id: str) -> tuple[str, gemmi.Residue] | None:
    pdb_id = pdb_id.lower()
    path = pathlib.Path(f"downloads/{pdb_id}.cif")

    if not path.is_file():
        print(f"Missing file: {path}")
        return None

    structure = gemmi.read_structure(str(path))

    best_residue = None
    best_score = float("-inf")

    for model in structure:
        for chain in model:
            for residue in chain:

                residue_name = residue.name.strip().upper()

                # 1. exclude known solvent / ion molecules
                if residue_name in DEFAULT_EXCLUDED:
                    continue

                # 2. skip protein residues
                if residue.het_flag == "A":
                    continue

                # 3. quick noise filter
                ha = heavy_atom_count(residue)
                if ha < 3:
                    continue

                # 4. score candidate
                score = ligand_score(residue)

                if score > best_score:
                    best_score = score
                    best_residue = residue

    if best_residue is None:
        print(f"No ligand found in {pdb_id}")
        return None

    return pdb_id, best_residue