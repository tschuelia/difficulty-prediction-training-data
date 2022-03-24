import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_iqtree_significance_tests_on_eval_trees():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/iqtree_significance_tests_on_eval_trees/data")
        expected_path = PurePosixPath(".tests/unit/iqtree_significance_tests_on_eval_trees/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)
        shutil.copy(".tests/12540_3.phy", workdir)
        shutil.copy(".tests/config.yaml", workdir)

        # dbg
        print("results/12540_3.phy/output_files/iqtree/significance.iqtree results/12540_3.phy/output_files/iqtree/significance.iqtree.log", file=sys.stderr)

        # Run the test job.
        sp.check_output([
            "python",
            "-m",
            "snakemake", 
            "results/12540_3.phy/output_files/iqtree/significance.iqtree",
            "results/12540_3.phy/output_files/iqtree/significance.iqtree.log",
            "-F", 
            "-j1",
            "--keep-target-files",
    
            "--directory",
            workdir,
        ])

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file), 
        # also see common.py.
        common.OutputChecker(data_path, expected_path, workdir).check()

        # TODO: change the behaviour: instead of comparing it byte by byte we want to check
        # - if .iqtree contains the same test results
        # - if .iqtree.log contains the correct function call
