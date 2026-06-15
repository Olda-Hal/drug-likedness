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
import tempfile

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


def download_structure(pdb_id: str, db: str,
                        error: bool, verbose: bool) -> str | None:
    # create a directory for the downloaded structure files if it doesn't exist
    pathlib.Path("downloads").mkdir(exist_ok=True)

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
        return pdb_id
    elif error:
        raise ValueError(f"Unable to download PDB with ID: {pdb_id}. Maybe it doesn't exist?")
    elif verbose:
        print(f"Unable to download PDB with ID: {pdb_id}. Maybe it doesn't exist? Skipping.")

    if verbose:
        print()


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


def export_results(data: list[tuple[str, gemmi.Residue | argparse.Namespace, LipinskiData]],
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


def visualize_results(data: list[tuple[str, gemmi.Residue | argparse.Namespace, LipinskiData]], compare: bool) -> None:
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
def retry_job() -> dict | None:
    query_path = pathlib.Path(tempfile.gettempdir()) / "drug_likeness_query.json"
    if query_path.exists():
        with open(query_path) as file:
            query_data = json.load(file)
        print("A previous query was found. Do you want to continue it? (y/N)")
        answer = input().strip().lower()
        if answer == "y":
            return query_data
        else:
            query_path.unlink()  # delete the file if the user doesn't want to continue

    return None

# deletes the old tempfile, and overwrites it with the new query data.
def start_job(args: argparse.Namespace) -> None:
    query_path = pathlib.Path(tempfile.gettempdir()) / "drug_likeness_query.json"
    if query_path.exists():
        query_path.unlink()  # delete old file

    query_data = vars(args).copy()
    query_data["results"] = []
    query_data["processed_pdb_ids"] = []

    with open(query_path, 'w') as file:
        json.dump(query_data, file, indent=4)
    
# exports the results to the output file and deletes the tempfile with the query data, so that it doesn't offer to continue the previous query on the next run.
def finish_job() -> None:
    query_path = pathlib.Path(tempfile.gettempdir()) / "drug_likeness_query.json"
    if query_path.exists():
        with open(query_path) as file:
            query_data = json.load(file)
        results = get_completed_result_data(query_data)
        if results:
            export_results(results, query_data.get("output_name"))

        query_path.unlink()

def update_query_results(new_results: list[tuple[str, gemmi.Residue | argparse.Namespace, LipinskiData]],
                        processed_pdb_ids: list[str]) -> None:
    query_path = pathlib.Path(tempfile.gettempdir()) / "drug_likeness_query.json"
    if query_path.exists():
        with open(query_path) as file:
            query_data = json.load(file)
        serialized_results = []
        for pdb_id, ligand, lipinski_data in new_results:
            molecular_weight, hb_donor_count, hb_acceptor_count, hb_donor_count_pubchem, hb_acceptor_count_pubchem, lipophilicity, satisfies_rules = lipinski_data

            curr_data = {
                "pdb_id": pdb_id,
                "ligand_name": ligand.name,
                "satisfies_lipinski": satisfies_rules,
                "data": {
                    "molecular_weight": molecular_weight,
                    "hydrogen_bond_donor_count": hb_donor_count,
                    "hydrogen_bond_acceptor_count": hb_acceptor_count,
                    "lipophilicity": lipophilicity
                }
            }

            if hb_donor_count_pubchem is not None and hb_acceptor_count_pubchem is not None:
                curr_data["pubchem_data"] = {
                    "hydrogen_bond_donor_count": hb_donor_count_pubchem,
                    "hydrogen_bond_acceptor_count": hb_acceptor_count_pubchem
                }

            serialized_results.append(curr_data)

        query_data["results"] = serialized_results
        query_data["processed_pdb_ids"] = processed_pdb_ids
        with open(query_path, 'w') as file:
            json.dump(query_data, file, indent=4)


def get_completed_result_data(query_data: dict | None) -> list[tuple[str, gemmi.Residue | argparse.Namespace, LipinskiData]]:
    if not query_data:
        return []

    completed_raw = query_data.get("results")
    if not isinstance(completed_raw, list):
        return []

    completed_data: list[tuple[str, gemmi.Residue | argparse.Namespace, LipinskiData]] = []
    for item in completed_raw:
        if not isinstance(item, dict):
            continue

        pdb_id = item.get("pdb_id")
        ligand_name = item.get("ligand_name")
        data = item.get("data")
        if not isinstance(pdb_id, str) or not isinstance(ligand_name, str) or not isinstance(data, dict):
            continue

        pubchem_data = item.get("pubchem_data")
        if not isinstance(pubchem_data, dict):
            pubchem_data = {}

        try:
            molecular_weight = float(data["molecular_weight"])
            hb_donor_count = int(data["hydrogen_bond_donor_count"])
            hb_acceptor_count = int(data["hydrogen_bond_acceptor_count"])
            lipophilicity_raw = data.get("lipophilicity")
            lipophilicity = None if lipophilicity_raw is None else float(lipophilicity_raw)
            satisfies_rules = bool(item["satisfies_lipinski"])
        except (KeyError, TypeError, ValueError):
            continue

        d_pc_raw = pubchem_data.get("hydrogen_bond_donor_count")
        a_pc_raw = pubchem_data.get("hydrogen_bond_acceptor_count")
        d_pc = None if d_pc_raw is None else int(d_pc_raw)
        a_pc = None if a_pc_raw is None else int(a_pc_raw)

        lipinski_data: LipinskiData = (
            molecular_weight,
            hb_donor_count,
            hb_acceptor_count,
            d_pc,
            a_pc,
            lipophilicity,
            satisfies_rules
        )

        # For resumed results we only need a name-like object for export formatting.
        completed_data.append((pdb_id, argparse.Namespace(name=ligand_name), lipinski_data))

    return completed_data

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
    parser.add_argument("-r", "--retry", action="store_true",
                        help="Enables retrying a previous query if it was interrupted. (Disabled by default.)")
    args = parser.parse_args()
    query_data: dict | None = None
    completed_result_data: list[tuple[str, gemmi.Residue | argparse.Namespace, LipinskiData]] = []
    completed_pdb_ids: set[str] = set()

    # ---------------------------------------------------------------------------------------

    # ------------------------------ [ Validating arguments ] -------------------------------
    if args.retry:
        query_data = retry_job()
        if query_data is not None:
            resumed_args = {k: v for k, v in query_data.items() if k not in {"results", "processed_pdb_ids"}}
            args = argparse.Namespace(**resumed_args)
            completed_result_data = get_completed_result_data(query_data)
            stored_processed_pdb_ids = query_data.get("processed_pdb_ids", [])
            if isinstance(stored_processed_pdb_ids, list):
                completed_pdb_ids = set(stored_processed_pdb_ids)
            else:
                completed_pdb_ids = set()
            if not completed_pdb_ids:
                completed_pdb_ids = {item[0] for item in completed_result_data}
    if args.ids is None and args.input_path is None:
        parser.error("No PDB IDs provided. Please specify either --ids or --input-path.")
    if args.ids and args.input_path:
        parser.error("Both --ids and --input-path were provided. Please specify only one.")

    pdb_ids = args.ids
    if args.input_path:
        pdb_ids = parse_input(args.input_path, args.verbose)
    # filter PDB IDs by their format validity
    pdb_ids = check_pdb_ids(pdb_ids, args.error, args.verbose)

    # In retry mode process only structures that are not already present in temp results.
    if completed_pdb_ids:
        pending_pdb_ids = [pdb_id for pdb_id in pdb_ids if pdb_id not in completed_pdb_ids]
        if args.verbose:
            skipped = len(pdb_ids) - len(pending_pdb_ids)
            if skipped > 0:
                print(f"Retry mode: skipping {skipped} already processed PDB IDs.")
        pdb_ids = pending_pdb_ids

    # create or continue the temp job state before processing structures one by one
    if not args.retry or query_data is None:
        start_job(args)

    # ---------------------------------------------------------------------------------------

    # ---------------------------------- [ Main pipeline ] ----------------------------------

    combined_result_data = completed_result_data.copy()
    processed_pdb_ids = set(completed_pdb_ids)

    for pdb_id in pdb_ids:
        downloaded_pdb_id = download_structure(pdb_id, args.db, args.error, args.verbose)
        processed_pdb_ids.add(pdb_id)
        if downloaded_pdb_id is None:
            update_query_results(combined_result_data, sorted(processed_pdb_ids))
            continue

        filtered_data = filter_ligands([downloaded_pdb_id])
        if len(filtered_data) == 0:
            if args.verbose:
                print(f"No significant ligand found in the PDB file {downloaded_pdb_id}.cif.")
            update_query_results(combined_result_data, sorted(processed_pdb_ids))
            continue

        _, ligand = filtered_data[0]
        lipinski_data = get_lipinski_data([ligand], args.compare, args.verbose)[0]
        combined_result_data.append((downloaded_pdb_id, ligand, lipinski_data))
        update_query_results(combined_result_data, sorted(processed_pdb_ids))

    if len(combined_result_data) == 0:
        raise RuntimeError("No significant ligands were found in the provided input PDBs. Exiting.")

    finish_job()

    if not args.disable_visualisation:
        visualize_results(combined_result_data, compare=args.compare)

    # ---------------------------------------------------------------------------------------
