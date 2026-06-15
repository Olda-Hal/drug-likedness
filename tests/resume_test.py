from __future__ import annotations

import json
import os
import pathlib
import shutil
import subprocess
import sys
import tempfile
import time
import unittest

# this test uses requests lib shimming to simulate network interactions and control the test environment
# (please dont ask, this is some dark magic...)

class ResumeTest(unittest.TestCase):
	def setUp(self) -> None:
		self.root_dir = pathlib.Path(__file__).resolve().parents[1]
		self.query_path = pathlib.Path(tempfile.gettempdir()) / "drug_likeness_query.json"
		self.output_name = "resume_test_output"
		self.output_path = self.root_dir / "results" / f"{self.output_name}.json"
		self.shim_dir = pathlib.Path(tempfile.mkdtemp(prefix="drug_likeness_requests_"))
		self._write_requests_shim()
		self._cleanup_artifacts()

	def tearDown(self) -> None:
		self._cleanup_artifacts()
		shutil.rmtree(self.shim_dir, ignore_errors=True)

	def _cleanup_artifacts(self) -> None:
		if self.query_path.exists():
			self.query_path.unlink()

		if self.output_path.exists():
			self.output_path.unlink()

	def _get_test_pdb_ids(self) -> list[str]:
		pdb_files = sorted(self.root_dir.joinpath("downloads").glob("pdb_*.cif"))
		pdb_ids = [file_path.stem for file_path in pdb_files]
		if len(pdb_ids) < 100:
			self.fail(f"Expected at least 100 local PDB files, found {len(pdb_ids)}.")
		return pdb_ids[:100]

	def _write_requests_shim(self) -> None:
		shim_code = f'''from __future__ import annotations

import pathlib
import time

ROOT_DIR = pathlib.Path({str(self.root_dir)!r})
DOWNLOAD_DIR = ROOT_DIR / "downloads"


class Response:
	def __init__(self, status_code: int, content: bytes = b"", payload: dict | None = None) -> None:
		self.status_code = status_code
		self.content = content
		self._payload = payload or {{}}

	def json(self) -> dict:
		return self._payload


def _pubchem_payload(url: str) -> dict:
	if "/XLogP/" in url:
		key = "XLogP"
		value: float | int = 2.5
	elif "/HBondDonorCount/" in url:
		key = "HBondDonorCount"
		value = 1
	else:
		key = "HBondAcceptorCount"
		value = 2

	return {{"PropertyTable": {{"Properties": [{{key: value}}]}}}}


def get(url: str) -> Response:
	time.sleep(0.01)
	if "files.rcsb.org/download/" in url or "www.ebi.ac.uk/pdbe/entry-files/download/" in url:
		file_name = url.rsplit("/", 1)[-1]
		content = (DOWNLOAD_DIR / file_name).read_bytes()
		return Response(200, content=content)

	if "pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/" in url:
		return Response(200, payload=_pubchem_payload(url))

	return Response(404)
'''
		(self.shim_dir / "requests.py").write_text(shim_code)

	def _run_main(self, *, resume: bool = False, input_text: str | None = None) -> subprocess.CompletedProcess[str]:
		env = os.environ.copy()
		pythonpath_parts = [str(self.shim_dir)]
		if env.get("PYTHONPATH"):
			pythonpath_parts.append(env["PYTHONPATH"])
		env["PYTHONPATH"] = os.pathsep.join(pythonpath_parts)

		command = [
			sys.executable,
			str(self.root_dir / "main.py"),
			"-i",
			*self._get_test_pdb_ids(),
			"-o",
			self.output_name,
			"-x",
		]
		if resume:
			command.append("-r")

		return subprocess.run(
			command,
			cwd=self.root_dir,
			env=env,
			input=input_text,
			text=True,
			capture_output=True,
			check=False,
		)

	def _run_and_kill_after_partial_progress(self, threshold: int = 20) -> subprocess.CompletedProcess[str]:
		env = os.environ.copy()
		pythonpath_parts = [str(self.shim_dir)]
		if env.get("PYTHONPATH"):
			pythonpath_parts.append(env["PYTHONPATH"])
		env["PYTHONPATH"] = os.pathsep.join(pythonpath_parts)

		command = [
			sys.executable,
			str(self.root_dir / "main.py"),
			"-i",
			*self._get_test_pdb_ids(),
			"-o",
			self.output_name,
			"-x",
		]

		process = subprocess.Popen(
			command,
			cwd=self.root_dir,
			env=env,
			stdin=subprocess.PIPE,
			stdout=subprocess.PIPE,
			stderr=subprocess.PIPE,
			text=True,
		)

		try:
			deadline = time.monotonic() + 30
			while time.monotonic() < deadline:
				if process.poll() is not None:
					stdout, stderr = process.communicate()
					self.fail(
						"The first run finished before it could be interrupted. "
						f"stdout={stdout!r} stderr={stderr!r}"
					)

				if self.query_path.exists():
					query_data = json.loads(self.query_path.read_text())
					processed = query_data.get("processed_pdb_ids", [])
					if isinstance(processed, list) and len(processed) >= threshold:
						process.kill()
						stdout, stderr = process.communicate(timeout=10)
						return subprocess.CompletedProcess(command, process.returncode, stdout, stderr)

				time.sleep(0.05)

			process.kill()
			stdout, stderr = process.communicate(timeout=10)
			return subprocess.CompletedProcess(command, process.returncode, stdout, stderr)
		finally:
			if process.poll() is None:
				process.kill()
				process.communicate(timeout=10)

	def test_resume_after_forced_interrupt(self) -> None:
		first_run = self._run_and_kill_after_partial_progress()
		self.assertNotEqual(first_run.returncode, 0, first_run.stderr)
		self.assertTrue(self.query_path.exists(), "Interrupted run should leave the resume state file behind.")

		query_data = json.loads(self.query_path.read_text())
		self.assertGreater(len(query_data["processed_pdb_ids"]), 0)
		self.assertLess(len(query_data["processed_pdb_ids"]), 100)
		self.assertGreater(len(query_data["results"]), 0)
		self.assertLessEqual(len(query_data["results"]), len(query_data["processed_pdb_ids"]))

		second_run = self._run_main(resume=True, input_text="y\n")
		self.assertEqual(second_run.returncode, 0, second_run.stderr)
		self.assertIn("A previous query was found.", second_run.stdout)
		self.assertFalse(self.query_path.exists(), "Successful resume should clean up the resume state file.")
		self.assertTrue(self.output_path.exists(), "Resume run should export the final results file.")

		results = json.loads(self.output_path.read_text())
		self.assertGreaterEqual(len(results), len(query_data["results"]))


if __name__ == "__main__":
	unittest.main()
