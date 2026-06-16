from rdkit import Chem
import requests
import gemmi
from rdkit.Chem import Crippen
from rdkit.Chem import rdmolops
from rdkit import RDLogger

RDLogger.DisableLog('rdApp.*')


def get_molecular_weight(ligand: gemmi.Residue) -> float:
    return sum(atom.element.weight for atom in ligand)


def get_hb_donor_count(ligand: gemmi.Residue) -> int:
    result = 0
    no_atoms = [  # nitrogen and oxygen atoms
                atom for atom in ligand
                if atom.element.name in ['N', 'O']
                ]
    hydrogen_atoms = [  # hydrogen atoms
                      atom for atom in ligand
                      if atom.is_hydrogen()
                     ]

    for atom in no_atoms:
        for h_atom in hydrogen_atoms:
            # a standard bond is ~1.0 Å, 1.2 Å provides a safer ceiling
            if atom.pos.dist(h_atom.pos) <= 1.2:
                result += 1
                break

    return result


def get_hb_acceptor_count(ligand: gemmi.Residue) -> int:
    result = 0
    for atom in ligand:
        if atom.element.name in ['N', 'O']:
            result += 1
    return result


def get_lipophilicity(ligand: gemmi.Residue) -> float | None:
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{ligand.name}/property/XLogP/json"

    response = requests.get(url)
    if response.status_code == 200:
        content = response.json()
        return content["PropertyTable"]["Properties"][0].get("XLogP")
    return None


def get_hb_donor_count_pubchem(ligand: gemmi.Residue) -> int | None:
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{ligand.name}/property/HBondDonorCount/json"

    response = requests.get(url)
    if response.status_code == 200:
        content = response.json()
        return content["PropertyTable"]["Properties"][0].get("HBondDonorCount")
    return None


def get_hb_acceptor_count_pubchem(ligand: gemmi.Residue) -> int | None:
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{ligand.name}/property/HBondAcceptorCount/json"

    response = requests.get(url)
    if response.status_code == 200:
        content = response.json()
        return content["PropertyTable"]["Properties"][0].get("HBondAcceptorCount")
    return None

def residue_to_rdkit_mol(residue: gemmi.Residue, cutoff=1.9):
    mol = Chem.RWMol()

    atoms = list(residue)

    atom_map = {}

    for i, atom in enumerate(atoms):
        rd_atom = Chem.Atom(atom.element.name)
        idx = mol.AddAtom(rd_atom)
        atom_map[i] = idx

    conf = Chem.Conformer(len(atoms))

    for i, atom in enumerate(atoms):
        pos = atom.pos
        conf.SetAtomPosition(i, (pos.x, pos.y, pos.z))

    mol.AddConformer(conf)

    n = len(atoms)

    for i in range(n):
        for j in range(i + 1, n):

            dist = atoms[i].pos.dist(atoms[j].pos)

            if dist < cutoff:
                if mol.GetBondBetweenAtoms(i, j) is None:
                    mol.AddBond(i, j, Chem.BondType.SINGLE)

    mol = mol.GetMol()

    try:
        Chem.SanitizeMol(mol)
    except Exception:
        try:
            Chem.SanitizeMol(
                mol,
                sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES
            )
        except Exception:
            pass

    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

    return mol


def get_lipophilicity_local(ligand: gemmi.Residue) -> float:
    mol = residue_to_rdkit_mol(ligand)
    return Crippen.MolLogP(mol)

def get_tolopogical_polar_surface_area(ligand: gemmi.Residue) -> float:
    mol = residue_to_rdkit_mol(ligand)
    return Chem.rdMolDescriptors.CalcTPSA(mol)

def get_rotatable_bond_count(ligand: gemmi.Residue) -> int:
    mol = residue_to_rdkit_mol(ligand)
    return Chem.rdMolDescriptors.CalcNumRotatableBonds(mol)

def get_ring_count(ligand: gemmi.Residue) -> int:
    mol = residue_to_rdkit_mol(ligand)
    return Chem.rdMolDescriptors.CalcNumRings(mol)

def get_fsp3(ligand: gemmi.Residue) -> float:
    mol = residue_to_rdkit_mol(ligand)
    return Chem.rdMolDescriptors.CalcFractionCSP3(mol)

def get_aromatic_ring_count(ligand: gemmi.Residue) -> int:
    mol = residue_to_rdkit_mol(ligand)
    return Chem.rdMolDescriptors.CalcNumAromaticRings(mol)

def get_formal_charge(ligand: gemmi.Residue) -> int:
    mol = residue_to_rdkit_mol(ligand)
    charge = rdmolops.GetFormalCharge(mol)
    return charge

def get_heavy_atom_count(ligand: gemmi.Residue) -> int:
    mol = residue_to_rdkit_mol(ligand)
    return mol.GetNumHeavyAtoms()

def get_atom_count(ligand: gemmi.Residue) -> int:
    mol = residue_to_rdkit_mol(ligand)
    return mol.GetNumAtoms()
