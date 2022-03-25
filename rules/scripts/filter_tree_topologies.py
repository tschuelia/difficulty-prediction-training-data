from custom_types import *
from utils import read_file_contents
from raxmlng_parser import get_raxmlng_num_unique_topos, get_raxmlng_abs_rf_distance

import pickle


def get_rfdist_clusters(log_path, all_trees):
    content = read_file_contents(log_path)
    clusters = []

    for line in content:
        line = line.strip()
        if not line.startswith("["):
            continue
        # the line is a string representation of a list of ints
        # like this: [1, 2, 3, 4, ]
        tree_ids = eval(line)
        cluster = set()

        for id in tree_ids:
            cluster.add(all_trees[id])

        clusters.append(cluster)

    return clusters



def filter_tree_topologies(
        eval_trees: List[Newick],
        log_path: FilePath,
):
    num_trees = len(eval_trees)

    if num_trees > 1:
        num_topos = get_raxmlng_num_unique_topos(log_path)

        if num_topos > 1:
            clusters = get_rfdist_clusters(log_path, eval_trees)

            assert len(clusters) == num_topos
        else:
            clusters = [set(eval_trees)]

    else:
        clusters = [set(eval_trees)]

    # for each cluster: keep only one tree as representative of the cluster
    unique_trees = [next(iter(cluster)) for cluster in clusters]

    # sanity checks
    assert sum([len(s) for s in clusters]) <= num_trees

    return unique_trees, clusters


if __name__ == "__main__":
    eval_trees = [l.strip() for l in open(snakemake.input.all_eval_trees).readlines()]
    log_file = snakemake.input.eval_trees_rfdistances_log

    unique_trees, clusters = filter_tree_topologies(
            eval_trees=eval_trees,
            log_path=log_file,
    )

    with open(snakemake.output.filtered_trees, "w") as f:
        f.write("\n".join(unique_trees))

    with open(snakemake.output.clusters, "wb") as f:
        pickle.dump(clusters, f)