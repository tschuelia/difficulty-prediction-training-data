from utils import read_file_contents
import regex


def get_all_parsimonator_parsimony_scores(log_file):
    content = read_file_contents(log_file)
    score_regex = regex.compile(r"Parsimony tree \[0\] with length (\d+) .*")

    scores = []

    for line in content:
        if "with length" not in line:
            continue
        line = line.strip()
        # Parsimony tree [0] with length 472 computed in 0.003233 seconds written to file:
        m = regex.match(score_regex, line)
        if m:
            score = m.groups()[0]
            scores.append(int(score))

    return scores


def get_all_parsimonator_runtimes(log_file):
    content = read_file_contents(log_file)
    runtime_regex = regex.compile(r"Parsimony tree \[0\] with length \d+ computed in (\d+(?:[\.]?\d+)?(?:[e][-+]?\d+)?) seconds")

    runtimes = []

    for line in content:
        if "with length" not in line:
            continue
        line = line.strip()
        # Parsimony tree [0] with length 472 computed in 0.003233 seconds written to file:
        m = regex.match(runtime_regex, line)
        if m:
            runtime = m.groups()[0]
            runtimes.append(float(runtime))

    return runtimes