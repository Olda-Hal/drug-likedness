import requests
import gemmi


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
