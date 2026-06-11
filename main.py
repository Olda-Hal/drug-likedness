import argparse
from turtle import title
import gemmi
import json
import pathlib
import re
import requests
from datetime import datetime
from property_funcs import *
import numpy as np
from visualisation import *

# mol weight, hb_donor_count, hb_acceptor_count, hb_donor_count_pubchem, hb_acceptor_count_pubchem, lipophilicity, satisfies_rules
LipinskiData = tuple[float, int, int, int | None, int | None, float | None, bool]
TITLES = ["Molecular Weight",
          "Hydrogen Bond Donor Count",
          "Hydrogen Bond Acceptor Count",
          "Lipophilicity", 
          "Satisfies Lipinski's Rule of Five"]

def parse_input(file_path: str, verbose: bool) -> list[str]:
    if pathlib.Path(file_path).is_file():
        pdb_ids = []
        with open(file_path) as file:
            for line in file.readlines():
                pdb_ids.extend(line.split())
        if verbose:
            print("Successfully read PDB IDs: ", " ".join(pdb_ids), "\n")
        return pdb_ids

    raise ValueError("Invalid input file path provided.")


def check_pdb_ids(pdb_ids: list[str], error: bool,
                  verbose: bool) -> list[str]:
    new_pdb_ids = []
    for pdb_id in pdb_ids:
        pdb_id = pdb_id.strip().lower()  # sanitize input
        if len(pdb_id) == 4:  # convert old ID format to new format
            pdb_id = "pdb_0000" + pdb_id
        if bool(re.match(r"^pdb_[0-9a-z]{8}$", pdb_id)):
            new_pdb_ids.append(pdb_id)
        elif error:
            raise ValueError(f"An invalid PDB ID was provided: {pdb_id}.")
        elif verbose:
            print(f"An invalid PDB ID was provided: {pdb_id}. Skipping.")

    return new_pdb_ids


def download_structures(pdb_ids: list[str], db: str,
                        error: bool, verbose: bool) -> list[str]:
    valid_pdb_ids = []
    # create a directory for the downloaded structure files if it doesn't exist
    pathlib.Path("downloads").mkdir(exist_ok=True)

    for pdb_id in pdb_ids:
        if db == "RCSB_PDB":
            url = f"https://files.rcsb.org/download/{pdb_id}.cif"
        else:
            url = f"https://www.ebi.ac.uk/pdbe/entry-files/download/{pdb_id}.cif"
        response = requests.get(url)

        if response.status_code == 200:
            with open(f"downloads/{pdb_id}.cif", 'wb') as file:
                file.write(response.content)
            if verbose:
                print(f"Successfully downloaded structure for PDB ID: {pdb_id}")
            valid_pdb_ids.append(pdb_id)
        elif error:
            raise ValueError(f"Unable to download PDB with ID: {pdb_id}. Maybe it doesn't exist?")
        elif verbose:
            print(f"Unable to download PDB with ID: {pdb_id}. Maybe it doesn't exist? Skipping.")

    if verbose:
        print()
    return valid_pdb_ids


def filter_ligands(pdb_ids: list[str]) -> list[tuple[str, gemmi.Residue]]:
    result_ligands = []
    excluded_molecules = {
        "HOH", "WAT", "DOD",  # water
        "GOL", "EDO", "DMS", "PEG", "ACT", "EDT",  # common solvents
        "SO4", "PO4", "CL", "NA", "K", "MG", "CA", "ZN"  # common ions
    }

    for pdb_id in pdb_ids:
        structure = gemmi.read_structure("downloads/" + pdb_id + ".cif")
        most_significant_ligand = None
        best_ha_count = 0  # highest heavy atom count

        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.het_flag == 'H' and residue.name \
                            not in excluded_molecules:  # it is a ligand
                        curr_ha_count = 0
                        for atom in residue:
                            if not atom.is_hydrogen():
                                curr_ha_count += 1
                        if curr_ha_count > best_ha_count:
                            best_ha_count = curr_ha_count
                            most_significant_ligand = residue

        if most_significant_ligand is not None:
            result_ligands.append((pdb_id, most_significant_ligand))
        else:
            print(f"No significant ligand found in the PDB file {pdb_id}.cif.")

    return result_ligands


def get_lipinski_data(ligands: list[gemmi.Residue], compare: bool,
                      verbose: bool) -> list[LipinskiData]:
    lipinski_data: list[LipinskiData] = []

    for ligand in ligands:
        molecular_weight = get_molecular_weight(ligand)
        hb_donor_count = get_hb_donor_count(ligand)
        hb_acceptor_count = get_hb_acceptor_count(ligand)
        lipophilicity = get_lipophilicity(ligand)

        violations = 0
        if molecular_weight >= 500: 
            violations += 1
        if lipophilicity is None or lipophilicity > 5:
            violations += 1
        if hb_donor_count > 5:
            violations += 1
        if hb_acceptor_count > 10:
            violations += 1

        satisfies_rules = violations <= 1

        if lipophilicity is None:
            print(f"Could not load the lipophilicity for ligand {ligand.name}, will assume lipophilicity > 5.")

        if verbose:
            if satisfies_rules:
                print(f"Ligand {ligand.name} does satisfy the Lipinski's rule of five.")
            else:
                print(f"Ligand {ligand.name} does not satisfy the Lipinski's rule of five.")

        if compare:
            hb_donor_count_pubchem = get_hb_donor_count_pubchem(ligand)
            hb_acceptor_count_pubchem = get_hb_acceptor_count_pubchem(ligand)
            if hb_donor_count_pubchem is None or hb_acceptor_count_pubchem is None:
                print("Could not load the PubChem data for comparison.")
                lipinski_data.append((molecular_weight, hb_donor_count,
                                      hb_acceptor_count, None, None, lipophilicity, satisfies_rules))
            else:
                donor_diff = abs(hb_donor_count - hb_donor_count_pubchem)
                acceptor_diff = abs(hb_acceptor_count - hb_acceptor_count_pubchem)
                print("Hydrogen Bond Donor Count")
                print(f"Calculated values: {hb_donor_count}, PubChem values: {hb_donor_count_pubchem}, difference: {donor_diff}")
                print("Hydrogen Bond Acceptor Count")
                print(f"Calculated values: {hb_acceptor_count}, PubChem values: {hb_acceptor_count_pubchem}, difference: {acceptor_diff}\n")
                lipinski_data.append((molecular_weight, hb_donor_count,
                                      hb_acceptor_count, hb_donor_count_pubchem,
                                      hb_acceptor_count_pubchem, lipophilicity, satisfies_rules))
        else:
            lipinski_data.append((molecular_weight, hb_donor_count,
                                  hb_acceptor_count, None, None, lipophilicity, satisfies_rules))

        if verbose:
            print()

    return lipinski_data


def export_results(data: list[tuple[str, gemmi.Residue, LipinskiData]],
                   output_filename: str | None) -> None:
    results = []

    for pdb_id, ligand, lipinski_data in data:
        M, D, A, D_pc, A_pc, L, satisfies_rules = lipinski_data
        lipinski_export = {
            "molecular_weight": M,
            "hydrogen_bond_donor_count": D,
            "hydrogen_bond_acceptor_count": A,
            "lipophilicity": L
        }
        curr_data = {
            "pdb_id": pdb_id,
            "ligand_name": ligand.name,
            "satisfies_lipinski": satisfies_rules,
            "data": lipinski_export
        }

        if D_pc is not None:
            pubchem_data = {
                "hydrogen_bond_donor_count": D_pc,
                "hydrogen_bond_acceptor_count": A_pc
            }
            curr_data["pubchem_data"] = pubchem_data

        results.append(curr_data)

    pathlib.Path("results").mkdir(parents=True, exist_ok=True)
    if output_filename is None:
        output_path = "results/" + datetime.now().strftime("%Y-%m-%d_%H-%M-%S") + ".json"
    else:
        output_path = "results/" + output_filename + ".json"

    with open(output_path, 'w') as file:
        json.dump(results, file, indent=4)


def visualize_results(data: list[tuple[str, gemmi.Residue, LipinskiData]], compare: bool) -> None:
    if not data:
        print("No data to visualize.")
        return

    # 1. Data preparation
    names = [item[0] for item in data]
    
    mol_weights = [item[2][0] for item in data]
    hb_donors_local = [item[2][1] for item in data]
    hb_acceptors_local = [item[2][2] for item in data]
    
    # Replace None values with np.nan for safe visualization (will remain completely empty)
    hb_donors_pubchem = [item[2][3] if item[2][3] is not None else np.nan for item in data]
    hb_acceptors_pubchem = [item[2][4] if item[2][4] is not None else np.nan for item in data]
    lipophilicities = [item[2][5] if item[2][5] is not None else np.nan for item in data]
    lipinsky_results = [item[2][6] for item in data]

    graph_generator(names, mol_weights, hb_donors_local, hb_acceptors_local,
                    hb_donors_pubchem, hb_acceptors_pubchem, lipophilicities,
                    lipinsky_results, compare)


# creates a file describing a query from the user. 
# this file will be deleted after final export. 
# If any file is found it will be able to continue the previous query in case of an error or a crash.
def query_job():
    pass

if __name__ == "__main__":
    # -------------------------------- [ Parsing arguments ] --------------------------------

    parser = argparse.ArgumentParser(description="Drug Likeness Analyzer")
    parser.add_argument("-f", "--input-path",
                        help="Specifies a path of an input file with PDB IDs.")
    parser.add_argument("-i", "--ids", nargs="+",
                        help="Specifies input PDB IDs. You can enter any amount.")
    parser.add_argument("-o", "--output-name",
                        help="Specifies a name for the output file. (Current date and time by default.)")
    parser.add_argument("-d", "--db", choices=["RCSB_PDB", "PDBe"],
                        help="Specifies a source database for ligand lookup.", default="RCSB_PDB")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Enables verbose mode. (Disabled by default.)")
    parser.add_argument("-e", "--error", action="store_true",
                        help="Toggles throwing an exception if an invalid PDB ID is provided. (Default: False, skips invalid IDs.)")
    parser.add_argument("-x", "--disable-visualisation", action="store_true",
                        help="Disables visualisation. (Enabled by default.)")
    parser.add_argument("-c", "--compare", action="store_true",
                        help="Enables comparing PubChem registry values with calculated values. (Disabled by default.)")
    args = parser.parse_args()

    # ---------------------------------------------------------------------------------------

    # ------------------------------ [ Validating arguments ] -------------------------------

    if args.ids is None and args.input_path is None:
        parser.error("No PDB IDs provided. Please specify either --ids or --input-path.")
    if args.ids and args.input_path:
        parser.error("Both --ids and --input-path were provided. Please specify only one.")

    pdb_ids = args.ids
    if args.input_path:
        pdb_ids = parse_input(args.input_path, args.verbose)
    # filter PDB IDs by their format validity
    pdb_ids = check_pdb_ids(pdb_ids, args.error, args.verbose)
    # keep only PDB IDs we found and downloaded the PDB file for
    pdb_ids = download_structures(pdb_ids, args.db, args.error, args.verbose)

    # ---------------------------------------------------------------------------------------

    # ---------------------------------- [ Main pipeline ] ----------------------------------

    # keep only PDB IDs for which we found a significant ligand
    filtered_data = filter_ligands(pdb_ids)
    if len(filtered_data) == 0:
        raise RuntimeError("No significant ligands were found in the provided input PDBs. Exiting.")
    pdb_ids_temp, ligands_temp = zip(*filtered_data)  # unpack new filtered data
    pdb_ids = list(pdb_ids_temp)
    ligands = list(ligands_temp)
    lipinski_data = get_lipinski_data(ligands, args.compare, args.verbose)
    result_data = list(zip(pdb_ids, ligands, lipinski_data))
    export_results(result_data, args.output_name)
    if not args.disable_visualisation:
        visualize_results(result_data, compare=args.compare)

    # ---------------------------------------------------------------------------------------
