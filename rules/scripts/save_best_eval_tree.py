from raxmlng_parser import get_all_raxmlng_llhs
from utils import read_file_contents


def get_best_tree_and_llh(eval_trees, eval_logs):
    all_trees = read_file_contents(eval_trees)
    all_llhs = get_all_raxmlng_llhs(eval_logs)

    # get the tree with the highest likelihood
    return max(zip(all_llhs, all_trees), key=lambda x: x[0])


if __name__ == "__main__":
    _, best_tree = get_best_tree_and_llh(
        eval_trees=snakemake.input.all_eval_trees,
        eval_logs=snakemake.input.all_eval_logs
    )

    open(snakemake.output.best_eval_tree, "w").write(best_tree)