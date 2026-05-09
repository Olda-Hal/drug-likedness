import argparse
from typing import List
from rcsbapi.data import DataQuery

def check_pdb_ids(pdbs: str, error_on_invalid_id: bool = False, is_file: bool = False):
    # TODO: pdb id sanity check
    pass

def download_structures(pdb_ids: List[str], db: str, error_on_invalid_id: bool = False):
    pass
    

def filter_ligands(structure_file):
    # TODO: Implement ligand filtering logic
    pass

def check_rules(ligands):
    # TODO: Implement rule checking logic
    pass

def export_results(results):
    # TODO: Implement result export logic
    pass

def visualize_results(results):
    # TODO: Implement result visualization logic;
    pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Drug Likeness Analyzer")
    parser.add_argument("-f", "--ids_file_path", help="path for file with PDB IDs to analyse", type=str)
    parser.add_argument("-I", "--ids", nargs="*", help="PDB Ids to analyse", type=str)
    parser.add_argument("-d", "--db", help="specify a DB for ligand lookup", default="RSCB PDB", type=str)
    parser.add_argument("-r", "--error_on_invalid_id", help="throws an exception if invalid ID is provided", action="store_false")
    args = parser.parse_args()
    if not args.ids and not args.ids_file_path:
        parser.error("No PDB IDs provided. Please specify either --ids or --ids_file_path.")
    if args.ids and args.ids_file_path:
        parser.error("Both --ids and --ids_file_path provided. Please specify only one.")
    print(args)
    
    if args.ids:
        ids = check_pdb_ids(args.ids, args.error_on_invalid_id, False)
    else:
        ids = check_pdb_ids(args.ids_file_path, args.error_on_invalid_id, True)
        
    structure_file = download_structures(["pdb_00001goj"], args.db)
    ligands = filter_ligands(structure_file)
    result = check_rules(ligands)
    export_results(result)
    visualize_results(result)