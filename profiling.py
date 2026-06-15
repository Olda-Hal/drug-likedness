import cProfile
import sys
from main import pipeline


if __name__ == "__main__":
    # forward command-line arguments (e.g. ligand IDs) to pipeline
    args = sys.argv[1:]
    prof = cProfile.Profile()
    original_argv = sys.argv
    sys.argv = [original_argv[0], *args]
    try:
        prof.runcall(pipeline)
    finally:
        sys.argv = original_argv
    prof.print_stats(sort="cumulative")