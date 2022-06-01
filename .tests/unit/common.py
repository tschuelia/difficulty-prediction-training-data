"""
Common code for unit testing of rules generated with Snakemake 6.13.1.
"""
import pytest
from Bio import Phylo
from io import StringIO
from pathlib import Path
import subprocess as sp
import os
import ete3
import warnings
import yaml
import pandas as pd
import pickle
import sqlite3
import json

import sys

sys.path.append(os.path.join(os.getcwd(), "rules/scripts"))

from raxmlng_parser import get_all_raxmlng_llhs, get_all_parsimony_scores


class OutputChecker:
    def __init__(self, data_path, expected_path, workdir):
        self.data_path = data_path
        self.expected_path = expected_path
        self.workdir = workdir

    def check(self):
        input_files = set(
            (Path(path) / f).relative_to(self.data_path)
            for path, subdirs, files in os.walk(self.data_path)
            for f in files
        )
        expected_files = set(
            (Path(path) / f).relative_to(self.expected_path)
            for path, subdirs, files in os.walk(self.expected_path)
            for f in files
        )
        unexpected_files = set()
        for path, subdirs, files in os.walk(self.workdir):
            for f in files:
                f = (Path(path) / f).relative_to(self.workdir)
                if (
                    str(f).startswith(".snakemake")
                    or str(f) == "config.yaml"
                    or str(f) == "12540_3.phy"
                ):
                    continue
                if f in expected_files:
                    self.compare_files(self.workdir / f, self.expected_path / f)
                elif f in input_files:
                    # ignore input files
                    pass
                else:
                    unexpected_files.add(f)
        if unexpected_files:
            warnings.warn(
                message="Unexpected files:\n{}".format(
                    "\n".join(sorted(map(str, unexpected_files))),
                ),
                category=UserWarning,
            )

    def compare_files(self, generated_file, expected_file):
        sp.check_output(["cmp", generated_file, expected_file])


class NewickTreeChecker(OutputChecker):
    def compare_tree_files(self, generated_tree_string, expected_tree_string):
        # first: check if the tree topologies are equal
        generated_tree = ete3.Tree(generated_tree_string)
        expected_tree = ete3.Tree(expected_tree_string)

        rfdist = generated_tree.robinson_foulds(expected_tree, unrooted_trees=True)[0]
        assert rfdist == pytest.approx(0)

        # second: compare if the branch lengths are approximately equal
        generated_tree = Phylo.read(StringIO(generated_tree_string), "newick")
        expected_tree = Phylo.read(StringIO(expected_tree_string), "newick")

        generated_branch_lenths = [
            node.branch_length
            for node in generated_tree.find_clades(branch_length=True)
        ]
        expected_branch_lenths = [
            node.branch_length for node in expected_tree.find_clades(branch_length=True)
        ]

        # sort the branch lengths before comparing them
        # this should cancel out effects of different display options for equal tree topologies
        generated_branch_lenths.sort()
        expected_branch_lenths.sort()

        assert all([gen == pytest.approx(exp, abs=0.01) for gen, exp in zip(generated_branch_lenths, expected_branch_lenths)])

    def compare_files(self, generated_file, expected_file):
        generated_trees = [l.strip() for l in open(generated_file).readlines() if l]
        exprected_trees = [l.strip() for l in open(expected_file).readlines() if l]

        for generated_tree, expected_tree in zip(generated_trees, exprected_trees):
            self.compare_tree_files(generated_tree, expected_tree)


class RAxMLNGLogChecker(OutputChecker):
    def __init__(self, data_path, expected_path, workdir, run_mode="inference"):
        super().__init__(data_path, expected_path, workdir)
        self.config = self.workdir.joinpath("config.yaml")
        with open(self.config) as f:
            self.raxmlng_command = yaml.safe_load(f)["software"]["raxml-ng"]["command"]

        self.run_mode = run_mode

    def _get_all_values(self, log_file, string_identifier):
        log = [l.strip() for l in open(log_file).readlines() if l]

        values = []

        for line in log:
            if string_identifier in line:
                values.append(line.rsplit(":", 1)[1].strip())

        return values

    def compare_log_files(self, generated_file, expected_file):
        if self.run_mode == "inference" or self.run_mode == "eval":
            generated_llhs = get_all_raxmlng_llhs(generated_file)
            expected_llhs = get_all_raxmlng_llhs(expected_file)

            assert all([gen == pytest.approx(exp, abs=0.1) for gen, exp in zip(generated_llhs, expected_llhs)])

        generated_pars_score = get_all_parsimony_scores(generated_file)
        expected_pars_score = get_all_parsimony_scores(expected_file)
        assert generated_pars_score == expected_pars_score

        generated_seeds = self._get_all_values(generated_file, "random seed:")
        expected_seeds = self._get_all_values(expected_file, "random seed:")

        assert generated_seeds == expected_seeds

        generated_run_modes = self._get_all_values(generated_file, "run mode:")
        expected_run_modes = self._get_all_values(expected_file, "run mode:")

        assert generated_run_modes == expected_run_modes

        generated_start_trees = self._get_all_values(generated_file, "start tree(s):")
        expected_start_trees = self._get_all_values(expected_file, "start tree(s):")

        assert generated_start_trees == expected_start_trees

    def compare_files(self, generated_file, expected_file):
        self.compare_log_files(generated_file, expected_file)


class RAxMLNGChecker(OutputChecker):
    def __init__(self, data_path, expected_path, workdir, run_mode="inference"):
        super().__init__(data_path, expected_path, workdir)
        self.config = self.workdir.joinpath("config.yaml")
        with open(self.config) as f:
            self.raxmlng_command = yaml.safe_load(f)["software"]["raxml-ng"]["command"]

        self.run_mode = run_mode

    def _compare_model_files(self, generated_file, expected_file):
        generated_model = open(generated_file).readline()
        expected_model = open(expected_file).readline()

        # TODO: this will most likely not work if I run this on another machine
        # the flaoting points might not be exactly identical
        # but for now this is good enough

        assert generated_model == expected_model

    def compare_files(self, generated_file, expected_file):
        if generated_file.suffix in [".bestTree", ".startTree"]:
            newick_checker = NewickTreeChecker(
                self.data_path, self.expected_path, self.workdir
            )
            newick_checker.compare_tree_files(
                open(generated_file).readline(), open(expected_file).readline()
            )
        elif generated_file.suffix == ".bestModel":
            self._compare_model_files(generated_file, expected_file)
        elif generated_file.suffix == ".log":
            log_checker = RAxMLNGLogChecker(self.data_path, self.expected_path, self.workdir, self.run_mode)
            log_checker.compare_log_files(generated_file, expected_file)


class RFDistanceChecker(OutputChecker):
    def _compare_log_files(self, generated_file, expected_file):
        generated_log = [l.strip() for l in open(generated_file).readlines() if l]
        expected_log = [l.strip() for l in open(expected_file).readlines() if l]

        for gen, exp in zip(generated_log, expected_log):
            identifiers = [
                "Average absolute RF distance in this tree set:",
                "Average relative RF distance in this tree set:",
                "Number of unique topologies in this tree set:",
            ]
            for identifier in identifiers:
                if identifier in gen:
                    assert identifier in exp
                    gen_rfdist = float(gen.rsplit(":", 1)[1])
                    exp_rfdist = float(exp.rsplit(":", 1)[1])

                    assert gen_rfdist == exp_rfdist

    def compare_files(self, generated_file, expected_file):
        if generated_file.suffix == ".log":
            self._compare_log_files(generated_file, expected_file)
        elif generated_file.suffix == ".rfDistances":
            # for the pairwise rfDistances use the parent's bytewise comparison
            super().compare_files(generated_file, expected_file)


class PandasDataFrameChecker(OutputChecker):
    def compare_files(self, generated_file, expected_file):
        # since the UUID is different each run, ignore this column when comparing the dataframes
        generated_df = pd.read_parquet(generated_file).drop("uuid", axis=1)
        expected_df = pd.read_parquet(expected_file).drop("uuid", axis=1)

        assert set(generated_df.columns) == set(expected_df.columns)

        columns_to_compare = [col for col in generated_df.columns if col != "difficult"]

        pd.testing.assert_frame_equal(generated_df[columns_to_compare], expected_df[columns_to_compare])


class IQTreeSignificanceTestChecker(OutputChecker):
    def _compare_result_files(self, generated_file, expected_file):
        generated_res = [l.strip() for l in open(generated_file).readlines() if l]
        expected_res = [l.strip() for l in open(expected_file).readlines() if l]

        result_start = -1
        result_end = -1

        for i, (gen, exp) in enumerate(zip(generated_res, expected_res)):
            if gen.startswith("Tree") and "logL" in gen:
                assert exp.startswith("Tree") and "logL" in exp
                result_start = i + 2
            elif "logL difference from the maximal logl in the set" in gen:
                assert "logL difference from the maximal logl in the set" in exp
                result_end = i

        for gen, exp in zip(generated_res[result_start:result_end], expected_res[result_start:result_end]):
            gen_res = gen.split()
            exp_res = exp.split()

            for g, e in zip(gen_res, exp_res):
                assert float(g) == float(e)

    def _compare_log_files(self, generated_file, expected_file):
        generated_log = [l.strip() for l in open(generated_file).readlines() if l]
        expected_log = [l.strip() for l in open(expected_file).readlines() if l]

        for gen, exp in zip(generated_log, expected_log):
            if gen.startswith("Seed:"):
                assert exp.startswith("Seed:")
                gen_seed = int(gen.split(":", 1)[1].split("(")[0].strip())
                exp_seed = int(exp.split(":", 1)[1].split("(")[0].strip())

                assert gen_seed == exp_seed

    def compare_files(self, generated_file, expected_file):
        if generated_file.suffix == ".iqtree":
            self._compare_result_files(generated_file, expected_file)
        elif generated_file.suffix == ".log":
            self._compare_log_files(generated_file, expected_file)


class FilteredClusterChecker(OutputChecker):
    def compare_files(self, generated_file, expected_file):
        if generated_file.suffix == ".pkl":
            generated_clusters = pickle.load(open(generated_file, "rb"))
            expected_clusters = pickle.load(open(expected_file, "rb"))

            for g in generated_clusters:
                assert g in expected_clusters

            for e in expected_clusters:
                assert e in generated_clusters
        elif generated_file.suffix == ".trees":
            newick_checker = NewickTreeChecker(
                self.data_path, self.expected_path, self.workdir
            )
            newick_checker.compare_files(generated_file, expected_file)


class DatabaseChecker(OutputChecker):
    def compare_files(self, generated_file, expected_file):
        if generated_file.suffix == ".sqlite3":
            gen_con = sqlite3.connect(generated_file)
            exp_con = sqlite3.connect(expected_file)

            for table in ["Dataset", "RaxmlNGTree", "ParsimonyTree"]:
                generated = pd.read_sql_query(f"SELECT * FROM {table}", gen_con)
                expected = pd.read_sql_query(f"SELECT * FROM {table}", exp_con)

                assert generated.shape == expected.shape
                assert generated.columns.tolist().sort() == expected.columns.tolist().sort()


class MSAFeatureChecker(OutputChecker):
    def compare_files(self, generated_file, expected_file):
        if generated_file.suffix == ".json":
            gen_json = json.load(open(generated_file))
            exp_json = json.load(open(expected_file))

            assert gen_json.keys() == exp_json.keys()

            for k, gen_v in gen_json.items():
                exp_v = exp_json.get(k)
                assert gen_v == exp_v