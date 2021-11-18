import regex

from utils import get_single_value_from_file, get_multiple_values_from_file, read_file_contents


def get_raxml_abs_rf_distance(log_file):
    STR = "Average absolute RF distance in this tree set:"
    return get_single_value_from_file(log_file, STR)


def get_raxml_rel_rf_distance(log_file):
    STR = "Average relative RF distance in this tree set:"
    return get_single_value_from_file(log_file, STR)


def get_raxml_num_unique_topos(log_file):
    STR = "Number of unique topologies in this tree set:"
    return get_single_value_from_file(log_file, STR)


def get_cleaned_rf_dist(raw_line):
    line_regex = regex.compile(r"(\d+)\s+(\d+)\s+(\d+)\s+(\d+\.\d+)\s*")
    tree_idx1, tree_idx2, plain_dist, normalized_dist = regex.search(line_regex, raw_line).groups()
    return int(tree_idx1), int(tree_idx2), float(plain_dist), float(normalized_dist)


def read_rfdistances(rfdistances_file_path):
    with open(rfdistances_file_path) as f:
        rfdistances = f.readlines()

    abs_res = {}
    rel_res = {}

    for line in rfdistances:
        idx1, idx2, plain, norm = get_cleaned_rf_dist(line)
        abs_res[(idx1, idx2)] = plain
        rel_res[(idx1, idx2)] = norm

    return abs_res, rel_res


def get_all_raxml_llhs(raxml_file):
    STR = "Final LogLikelihood:"
    return get_multiple_values_from_file(raxml_file, STR)


def get_best_raxml_llh(raxml_file):
    all_llhs = get_all_raxml_llhs(raxml_file)
    return max(all_llhs)


def get_raxml_elapsed_time(log_file):
    content = read_file_contents(log_file)

    all_times = []

    for line in content:
        if "Elapsed time:" not in line:
            continue
        # two cases now:
        # either the run was cancelled an rescheduled
        if "restarts" in line:
            # line looks like this: "Elapsed time: 5562.869 seconds (this run) / 91413.668 seconds (total with restarts)"
            _, right = line.split("/")
            value = right.split(" ")[1]
            all_times.append(float(value))

        # ...or the run ran in one sitting...
        else:
            # line looks like this: "Elapsed time: 63514.086 seconds"
            value = line.split(" ")[2]
            all_times.append(float(value))

    if not all_times:
        raise ValueError(
            f"The given input file {log_file} does not contain the elapsed time."
        )

    return all_times
