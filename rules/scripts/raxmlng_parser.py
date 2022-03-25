import regex
from tempfile import TemporaryDirectory
import subprocess

from custom_types import *
from utils import (
    get_single_value_from_file,
    get_multiple_values_from_file,
    read_file_contents,
)


def get_raxmlng_abs_rf_distance(log_file: FilePath) -> float:
    STR = "Average absolute RF distance in this tree set:"
    return get_single_value_from_file(log_file, STR)


def get_raxmlng_rel_rf_distance(log_file: FilePath) -> float:
    STR = "Average relative RF distance in this tree set:"
    return get_single_value_from_file(log_file, STR)


def get_raxmlng_num_unique_topos(log_file: FilePath) -> float:
    STR = "Number of unique topologies in this tree set:"
    return get_single_value_from_file(log_file, STR)


def get_cleaned_rf_dist(raw_line: str) -> Tuple[TreeIndex, TreeIndex, float, float]:
    line_regex = regex.compile(r"(\d+)\s+(\d+)\s+(\d+)\s+(\d+\.\d+)\s*")
    tree_idx1, tree_idx2, plain_dist, normalized_dist = regex.search(
        line_regex, raw_line
    ).groups()
    return int(tree_idx1), int(tree_idx2), float(plain_dist), float(normalized_dist)


def read_rfdistances(
    rfdistances_file_path: FilePath,
) -> Tuple[TreeTreeIndexed, TreeTreeIndexed]:
    with open(rfdistances_file_path) as f:
        rfdistances = f.readlines()

    abs_res = {}
    rel_res = {}

    for line in rfdistances:
        idx1, idx2, plain, norm = get_cleaned_rf_dist(line)
        abs_res[(idx1, idx2)] = plain
        rel_res[(idx1, idx2)] = norm

    return abs_res, rel_res


def get_raxmlng_llh(raxmlng_file: FilePath) -> float:
    STR = "Final LogLikelihood:"
    return get_single_value_from_file(raxmlng_file, STR)


def get_raxmlng_starting_llh(raxmlng_file: FilePath) -> float:
    content = read_file_contents(raxmlng_file)
    for line in content:
        if "Initial branch length optimization" not in line:
            continue
        # [00:00:00 -8735.928562] Initial branch length optimization
        numbers, _ = line.split("]")
        _, llh = numbers.split(" ")
        return float(llh)

    raise ValueError("The given file does not contain the starting llh")


def get_all_raxmlng_llhs(raxmlng_file: FilePath) -> List[float]:
    STR = "Final LogLikelihood:"
    return get_multiple_values_from_file(raxmlng_file, STR)


def get_best_raxmlng_llh(raxmlng_file: FilePath) -> float:
    all_llhs = get_all_raxmlng_llhs(raxmlng_file)
    return max(all_llhs)


def get_raxmlng_time_from_line(line: str) -> float:
    # two cases now:
    # either the run was cancelled an rescheduled
    if "restarts" in line:
        # line looks like this: "Elapsed time: 5562.869 seconds (this run) / 91413.668 seconds (total with restarts)"
        _, right = line.split("/")
        value = right.split(" ")[1]
        return float(value)

    # ...or the run ran in one sitting...
    else:
        # line looks like this: "Elapsed time: 63514.086 seconds"
        value = line.split(" ")[2]
        return float(value)


def get_raxmlng_elapsed_time(log_file: FilePath) -> float:
    content = read_file_contents(log_file)

    for line in content:
        if "Elapsed time:" not in line:
            continue
        else:
            return get_raxmlng_time_from_line(line)

    raise ValueError(
        f"The given input file {log_file} does not contain the elapsed time."
    )


def get_raxmlng_runtimes(log_file: FilePath) -> List[float]:
    content = read_file_contents(log_file)

    all_times = []

    for line in content:
        if "Elapsed time:" not in line:
            continue
        else:
            all_times.append(get_raxmlng_time_from_line(line))

    if not all_times:
        raise ValueError(
            f"The given input file {log_file} does not contain the elapsed time."
        )

    return all_times


def get_raxmlng_num_spr_rounds(log_file: FilePath) -> Tuple[int, int]:
    slow_regex = regex.compile("SLOW\s+spr\s+round\s+(\d+)")
    fast_regex = regex.compile("FAST\s+spr\s+round\s+(\d+)")

    content = read_file_contents(log_file)

    slow = 0
    fast = 0

    for line in content:
        if "SLOW spr round" in line:
            m = regex.search(slow_regex, line)
            if m:
                slow_round = int(m.groups()[0])
                slow = max(slow, slow_round)
        if "FAST spr round" in line:
            m = regex.search(fast_regex, line)
            if m:
                fast_round = int(m.groups()[0])
                fast = max(fast, fast_round)

    return slow, fast


def rel_rfdistance_starting_final(
    newick_starting: Newick,
    newick_final: Newick,
    raxmlng_executable: Command = "raxml-ng",
) -> float:
    with TemporaryDirectory() as tmpdir:
        trees = tmpdir + ".trees"
        with open(trees, "w") as f:
            f.write(newick_starting.strip() + "\n" + newick_final.strip())

        cmd = [raxmlng_executable, "--rfdist", trees, "--prefix", tmpdir]

        subprocess.check_output(cmd)

        log_file = tmpdir + ".raxml.log"

        return get_raxmlng_rel_rf_distance(log_file)


def get_model_parameter_estimates(raxmlng_file: FilePath) -> Tuple[str, str, str]:
    """
    For now just store everyting as string, different models result in different strings
    and I don't want to commit to parsing just yet
    TODO: adapt this for multiple partitions, in this case return a dict instead of a string with the partition name as key and the corresponding string as value
    """
    content = read_file_contents(raxmlng_file)

    rate_het = None
    base_freq = None
    subst_rates = None

    for line in content:
        if line.startswith("Rate heterogeneity"):
            _, res = line.split(":", 1)
            rate_het = res.strip()
        if line.startswith("Base frequencies"):
            _, res = line.split(":", 1)
            base_freq = res.strip()
        if line.startswith("Substitution rates"):
            _, res = line.split(":", 1)
            subst_rates = res.strip()

    return rate_het, base_freq, subst_rates


def get_all_parsimony_scores(log_file: FilePath) -> List[float]:
    content = read_file_contents(log_file)

    scores = []

    for line in content:
        if "Parsimony score" not in line:
            continue

        _, score = line.split(":")
        score = int(score.strip())
        scores.append(score)

    return scores


def get_patterns_gaps_invariant(log_file: FilePath) -> Tuple[int, float, float]:
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
