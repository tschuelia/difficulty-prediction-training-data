import pickle
from iqtree_statstest_parser import get_iqtree_results, get_iqtree_results_for_eval_tree_str


iqtree_results = get_iqtree_results(snakemake.input.iqtree_results)
with open(snakemake.input.clusters, "rb") as f:
    clusters = pickle.load(f)

eval_trees = open(snakemake.input.eval_trees).readlines()
plausible_trees = []

for newick_eval in eval_trees:
    if not newick_eval:
        continue

    statstest_results, _ = get_iqtree_results_for_eval_tree_str(iqtree_results, newick_eval, clusters)
    plausible = statstest_results["plausible"]

    if plausible:
        plausible_trees.append(newick_eval.strip())

with open(snakemake.output.all_plausible_trees, "w") as f:
    f.write("\n".join(plausible_trees))
