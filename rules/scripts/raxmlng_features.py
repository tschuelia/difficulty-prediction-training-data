from tempfile import TemporaryDirectory

from custom_types import *
from utils import run_cmd


class RAxMLNG:
    def __init__(self, exe_path: ExePath):
        self.exe_path = exe_path

    def __base_cmd(
        self, msa_file: FilePath, model: Model, prefix: FilePath, **kwargs
    ) -> Command:
        additional_settings = []
        for key, value in kwargs.items():
            if value is None:
                additional_settings += [f"--{key}"]
            else:
                additional_settings += [f"--{key}", str(value)]

        return [
            self.exe_path,
            "--msa",
            msa_file,
            "--model",
            model,
            "--prefix",
            prefix,
            *additional_settings,
        ]

    def run_alignment_parse(
        self, msa_file: FilePath, model: Model, prefix: str, **kwargs
    ) -> None:
        cmd = self.__base_cmd(msa_file, model, prefix, parse=None, **kwargs)
        run_cmd(cmd)

    def infer_parsimony_trees(
        self,
        msa_file: FilePath,
        model: Model,
        prefix: str,
        n_trees: int = 100,
        **kwargs,
    ) -> None:
        cmd = self.__base_cmd(
            msa_file, model, prefix, start=None, tree=f"pars{{{n_trees}}}", **kwargs
        )
        run_cmd(cmd)

    def parse_rfdist_log(self, log_file):
        for line in open(log_file).readlines():
            line = line.strip()
            if "Average relative RF distance in this tree set:" not in line:
                continue
            _, rfdist = line.rsplit(" ", 1)
            return float(rfdist)

        raise ValueError("Relative RF Distance not found in file", log_file)

    def get_rel_rfdistance(self, trees_file: FilePath, **kwargs) -> float:
        additional_settings = []
        for key, value in kwargs.items():
            if value is None:
                additional_settings += [f"--{key}"]
            else:
                additional_settings += [f"--{key}", str(value)]

        cmd = [
            self.exe_path,
            "--rfdist",
            trees_file,
            *additional_settings,
        ]
        run_cmd(cmd)

        return self.parse_rfdist_log(trees_file + ".raxml.log")

    def parse_log(self, log_file: FilePath) -> Tuple[int, float, float]:
        patterns = None
        gaps = None
        invariant = None
        for line in open(log_file).readlines():
            if line.startswith("Alignment sites"):
                # number of alignment patterns
                # Alignment sites / patterns: 1940 / 933
                _, numbers = line.split(":")
                _, patterns = [int(el) for el in numbers.split("/")]
            elif line.startswith("Gaps"):
                # proportion of gaps
                _, number = line.split(":")
                percentage, _ = number.strip().split(" ")
                gaps = float(percentage) / 100.0
            elif line.startswith("Invariant sites"):
                # proportion invariant sites
                _, number = line.split(":")
                percentage, _ = number.strip().split(" ")
                invariant = float(percentage) / 100.0

        if patterns is None or gaps is None or invariant is None:
            raise ValueError("Error parsing raxml-ng log")

        return patterns, gaps, invariant

    def get_patterns_gaps_invariant(
        self, msa_file: FilePath, model: Model
    ) -> Tuple[int, float, float]:
        with TemporaryDirectory() as tmpdir:
            prefix = tmpdir + "/parse"
            self.run_alignment_parse(msa_file, model, prefix)
            return self.parse_log(f"{prefix}.raxml.log")
