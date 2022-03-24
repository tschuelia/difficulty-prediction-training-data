import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_raxmlng_rfdistance_search_trees():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/raxmlng_rfdistance_search_trees/data")
        expected_path = PurePosixPath(".tests/unit/raxmlng_rfdistance_search_trees/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)
        shutil.copy(".tests/12540_3.phy", workdir)
        shutil.copy(".tests/config.yaml", workdir)

        # dbg
        print("results/12540_3.phy/output_files/raxmlng/inference/inference.raxml.rfDistances results/12540_3.phy/output_files/raxmlng/inference/inference.raxml.rfDistances.log", file=sys.stderr)

        # Run the test job.
        sp.check_output([
            "python",
            "-m",
            "snakemake", 
            "results/12540_3.phy/output_files/raxmlng/inference/inference.raxml.rfDistances",
            "results/12540_3.phy/output_files/raxmlng/inference/inference.raxml.rfDistances.log",
            "-F", 
            "-j1",
            "--keep-target-files",
    
            "--directory",
            workdir,
        ])

        common.RFDistanceChecker(data_path, expected_path, workdir).check()
