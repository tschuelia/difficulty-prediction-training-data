from raxml_parser import get_all_raxml_llhs
from utils import read_file_contents

all_trees = read_file_contents(snakemake.input.all_eval_trees)
all_llhs = get_all_raxml_llhs(snakemake.input.all_eval_logs)

# get the tree with the highest likelihood
best_tree = sorted(zip(all_llhs, all_trees), key=lambda x: x[0])[-1][1]

open(snakemake.output.best_eval_tree, "w").write(best_tree)