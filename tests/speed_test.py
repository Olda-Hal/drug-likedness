from __future__ import annotations

import pathlib
import subprocess
import sys
import time
import unittest
from typing import Iterable


PDB_IDS = [
	"1TUP", "1HHO", "1BNA", "1CRN", "1AKE", "1L2Y", "1CAG", "1D66", "1E8L", "1FAT",
	"1GFL", "1HIV", "1I1B", "1J8K", "1K5N", "1LMB", "1MBS", "1N0R", "1O3T", "1PGA",
	"1Q2W", "1R6J", "1S0Q", "1TND", "1U9S", "1V54", "1W9R", "1X8B", "1YIT", "1ZIH",
	"2A0B", "2B6A", "2C7D", "2D3U", "2E4E", "2F4K", "2G5U", "2H5J", "2I2A", "2J3Q",
	"2K39", "2L6O", "2MNR", "2N8D", "2O4I", "2P4A", "2Q5S", "2R1R", "2SIC", "2TAA",
	"3A0U", "3B5V", "3C7D", "3D3Z", "3E8Y", "3F6M", "3G5U", "3H0H", "3I4V", "3J3Y",
	"3K0N", "3L1C", "3M8O", "3NIR", "3O5Q", "3P0G", "3Q4A", "3R2X", "3S5M", "3T4L",
	"4A0I", "4B7Q", "4C3N", "4D5S", "4E2Z", "4F1A", "4G6K", "4H2D", "4I9P", "4J8R",
	"4K3L", "4L7M", "4M2C", "4N6E", "4O4F", "4P2J", "4Q9V", "4R1B", "4S7H", "4T5Y",
	"5A0M", "5B3T", "5C2K", "5D7N", "5E1P", "5F4R", "5G2W", "5H6A", "5I8D", "5J1F",
]



def _run_speed_test(pdb_ids: Iterable[str]) -> float:
	root_dir = pathlib.Path(__file__).resolve().parents[1]
	command = [
		sys.executable,
		str(root_dir / "main.py"),
		"-i",
		*list(pdb_ids),
		"-x",
	]

	start = time.perf_counter()
	result = subprocess.run(
		command,
		cwd=root_dir,
		capture_output=True,
		text=True,
		check=False,
	)
	duration = time.perf_counter() - start

	if result.returncode != 0:
		raise subprocess.CalledProcessError(
			result.returncode,
			command,
			output=result.stdout,
			stderr=result.stderr,
		)

	return duration


class SpeedTest(unittest.TestCase):
	def test_speed_pdb_ids(self) -> None:
		duration = _run_speed_test(PDB_IDS)
		per_id = duration / len(PDB_IDS)

		print(f"Speed test over {len(PDB_IDS)} PDB IDs")
		print(f"Total time: {duration:.4f} s")
		print(f"Average per 1 ID: {per_id:.6f} s")
		print(f"Estimated for 100 IDs: {per_id * 100:.4f} s")

