# drug-likeness

This repo contains a tool for filtering molecules based on their drug-likeness properties. The tool uses various criteria to evaluate the potential of a molecule to be a viable drug candidate.

# installation

To install the tool, follow these steps:
```
git clone https://github.com/Olda-hal/drug-likedness.git
python -m venv venv
source /venv/bin/activate
pip install -r requirements.txt

# or if you prefere to support corporations like OpenAI, use UV:
uv sync
```

# usage

to run the command line tool, use the following command:
```
./main.py <ARGS>
```

following are the available arguments:
| Option                          | Description                                                                                                                                                         |
| ------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `-f`, `--input-path`            | Path to a file containing PDB IDs to process.                                                                                                                       |
| `-i`, `--ids`                   | One or more PDB IDs provided directly via the command line. Multiple IDs can be specified.                                                                          |
| `-o`, `--output-name`           | Name of the output file. If not specified, a timestamp-based name is used.                                                                                          |
| `-d`, `--db`                    | Source database used for ligand lookup. Supported values: `RCSB_PDB` (default) and `PDBe`.                                                                          |
| `-v`, `--verbose`               | Enables verbose logging and additional diagnostic output.                                                                                                           |
| `-e`, `--error`                 | Raises an exception when an invalid PDB ID is encountered. By default, invalid IDs are skipped.                                                                     |
| `-x`, `--disable-visualisation` | Disables generation of visualizations. Visualizations are enabled by default.                                                                                       |
| `-c`, `--compare`               | Compares calculated molecular properties with values obtained from the PubChem registry.                                                                            |
| `-r`, `--retry`                 | Attempts to resume and continue a previously interrupted query.                                                                                                     |
| `-q`, `--quick`                 | Runs the tool in quick mode for faster execution. Disables several non-essential checks and features, including progress and resume state saving. **Experimental.** |

example usage:
```
# Process a single PDB entry
./main.py --ids 1CRN

# Process multiple PDB entries
./main.py --ids 1CRN 4HHB 2PTC

# Read PDB IDs from a file
./main.py --input-path pdb_ids.txt

# Use the PDBe database with verbose output
./main.py --ids 1CRN --db PDBe --verbose

# Run in quick mode
./main.py --ids 1CRN --quick

# Resume a previously interrupted run
./main.py --retry
```

# results
you will find the results in the `results` directory. The output will be JSON files containing the calculated properties for each processed PDB entry. If visualizations are enabled, you will be also shown the graph.