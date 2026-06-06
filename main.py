import argparse
from typing import List
import requests
import pathlib


def check_pdb_ids(pdbs: str, error_on_invalid_id: bool = False, is_file: bool = False) -> List[str]:
    # TODO: pdb id sanity check
    pass


def pdb_exists(pdb_id: str) -> bool:
    pdb_id = pdb_id.lower()

    url = f"https://files.rcsb.org/download/{pdb_id}.cif"

    r = requests.head(url)

    return r.status_code == 200


def download_structures(pdb_ids: List[str], db: str, error_on_invalid_id: bool = False) -> List[pathlib.Path]:
    downloaded_files: List[pathlib.Path] = []
    for pdb_id in pdb_ids:
        if not pdb_exists(pdb_id):
            if error_on_invalid_id:
                raise ValueError(f"Invalid PDB ID: {pdb_id}")
            else:
                print(f"Warning: PDB ID {pdb_id} does not exist. Skipping.")
                continue

        # creating a folder for the downloaded structures if it doesn't exist
        pathlib.Path("structures").mkdir(parents=True, exist_ok=True)
        match db:
            case "RSCB PDB":
                url = f"https://files.rcsb.org/download/{pdb_id}.cif"
            case 
        response = requests.get(url)
        if response.status_code == 200:
            with open(f"structures/{pdb_id}.cif", "wb") as f:
                f.write(response.content)
            print(f"Downloaded structure for PDB ID: {pdb_id}")
            downloaded_files.append(pathlib.Path(f"structures/{pdb_id}.cif"))
        else:
            print(f"Failed to download structure for PDB ID: {pdb_id}")

    return downloaded_files


def filter_ligands(structure_files: List[pathlib.Path]):
    # TODO: Implement ligand filtering logic
    pass


def check_rules(ligands: List[str]):
    # TODO: Implement rule checking logic
    pass


def export_results(results: List[str]):
    # TODO: Implement result export logic
    pass


def visualize_results(results: List[str]):
    # TODO: Implement result visualization logic;
    pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Drug Likeness Analyzer")
    parser.add_argument("-f", "--ids_file_path",
                        help="path for file with PDB IDs to analyse", type=str)
    parser.add_argument("-I", "--ids", nargs="*",
                        help="PDB Ids to analyse", type=str)
    parser.add_argument(
        "-d", "--db", help="specify a DB for ligand lookup", default="RSCB PDB", type=str)
    parser.add_argument("-r", "--error_on_invalid_id",
                        help="throws an exception if invalid ID is provided", action="store_false")
    args = parser.parse_args()
    print(args)
    if not args.ids and not args.ids_file_path:
        parser.error(
            "No PDB IDs provided. Please specify either --ids or --ids_file_path.")
    if args.ids and args.ids_file_path:
        parser.error(
            "Both --ids and --ids_file_path provided. Please specify only one.")

    if args.ids:
        ids = check_pdb_ids(args.ids, args.error_on_invalid_id, False)
    else:
        ids = check_pdb_ids(args.ids_file_path, args.error_on_invalid_id, True)

    structure_files = download_structures(
        args.ids, args.db, args.error_on_invalid_id)
    ligands = filter_ligands(structure_files)
    result = check_rules(ligands)
    export_results(result)
    visualize_results(result)
