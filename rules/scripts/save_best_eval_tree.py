from raxmlng_parser import get_all_raxmlng_llhs
from utils import read_file_contents

all_trees = read_file_contents(snakemake.input.all_eval_trees)
all_llhs = get_all_raxmlng_llhs(snakemake.input.all_eval_logs)

# get the tree with the highest likelihood
_, best_tree = max(zip(all_llhs, all_trees), key=lambda x: x[0])

open(snakemake.output.best_eval_tree, "w").write(best_tree)