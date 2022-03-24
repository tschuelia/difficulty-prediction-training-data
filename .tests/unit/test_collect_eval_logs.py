import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_collect_eval_logs():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/collect_eval_logs/data")
        expected_path = PurePosixPath(".tests/unit/collect_eval_logs/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)
        shutil.copy(".tests/12540_3.phy", workdir)
        shutil.copy(".tests/config.yaml", workdir)

        # dbg
        print("results/12540_3.phy/output_files/raxmlng/evaluation/AllEvalLogs.log", file=sys.stderr)

        # Run the test job.
        sp.check_output([
            "python",
            "-m",
            "snakemake", 
            "results/12540_3.phy/output_files/raxmlng/evaluation/AllEvalLogs.log",
            "-F", 
            "-j1",
            "--keep-target-files",
    
            "--directory",
            workdir,
        ])

        common.RAxMLNGLogChecker(data_path, expected_path, workdir, run_mode="eval").check()
