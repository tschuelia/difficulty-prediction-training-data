from iqtree_statstest_parser import get_iqtree_results

iqtree_results = get_iqtree_results(snakemake.input.iqtree_results)
eval_trees = open(snakemake.input.eval_trees).readlines()

if len(iqtree_results) != len(eval_trees):
    raise ValueError(f"The numer of IQ-Tree test results ({len(iqtree_results)}) "
                     f"does not match the number of eval trees ({len(eval_trees)})."
                     f"Check the IQ-Tree log file ({snakemake.input.iqtree_results}) and the parser.")

plausible_trees = []

for newick_eval, test_results in zip(eval_trees, iqtree_results):
    if test_results["plausible"]:
        plausible_trees.append(newick_eval.strip())

with open(snakemake.output.all_plausible_trees, "w") as f:
    f.write("\n".join(plausible_trees))
