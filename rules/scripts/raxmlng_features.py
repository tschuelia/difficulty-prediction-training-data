from tempfile import TemporaryDirectory

from utils import run_cmd


class RAxMLNG:
    def __init__(self, exe_path):
        self.exe_path = exe_path

    def __base_cmd(self, msa_file, model, prefix, **kwargs):
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

    def run_alignment_parse(self, msa_file, model, prefix, **kwargs):
        cmd = self.__base_cmd(msa_file, model, prefix, parse=None, **kwargs)
        run_cmd(cmd)

    def infer_parsimony_trees(self, msa_file, model, prefix, n_trees=100, **kwargs):
        cmd = self.__base_cmd(msa_file, model, prefix, start=None, tree=f"pars{{{n_trees}}}", **kwargs)
        run_cmd(cmd)

    def get_rel_rfdistance(self, trees_file, **kwargs):
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

        outlog = trees_file + ".raxml.log"

        for line in open(outlog).readlines():
            line = line.strip()
            if "Average relative RF distance in this tree set:" not in line:
                continue
            _, rfdist = line.rsplit(" ", 1)
            return float(rfdist)

        raise ValueError("Relative RF Distance not found in file", outlog)

    def get_patterns_gaps_invariant(self, msa_file, model):
        with TemporaryDirectory() as tmpdir:
            prefix = tmpdir + "/parse"
            self.run_alignment_parse(msa_file, model, prefix)
            patterns = None
            gaps = None
            invariant = None
            for line in open(f"{prefix}.raxml.log").readlines():
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