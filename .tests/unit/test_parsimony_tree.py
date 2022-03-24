import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_parsimony_tree():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/parsimony_tree/data")
        expected_path = PurePosixPath(".tests/unit/parsimony_tree/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)
        shutil.copy(".tests/12540_3.phy", workdir)
        shutil.copy(".tests/config.yaml", workdir)

        # dbg
        print("results/12540_3.phy/output_files/parsimony/seed_1.raxml.startTree results/12540_3.phy/output_files/parsimony/seed_1.raxml.log", file=sys.stderr)

        # Run the test job.
        sp.check_output([
            "python",
            "-m",
            "snakemake", 
            "results/12540_3.phy/output_files/parsimony/seed_1.raxml.startTree",
            "results/12540_3.phy/output_files/parsimony/seed_1.raxml.log",
            "-F", 
            "-j1",
            "--keep-target-files",
    
            "--directory",
            workdir,
        ])

        common.RAxMLNGChecker(data_path, expected_path, workdir, run_mode="start").check()