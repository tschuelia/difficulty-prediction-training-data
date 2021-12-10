import os

configfile: "config.yaml"

raxmlng_command = config["software"]["raxml-ng"]["command"]
iqtree_command = config["software"]["iqtree"]["command"]
parsimonator_command = config["software"]["parsimonator"]["command"]

num_pars_trees = config["_debug"]["_num_pars_trees"]
num_rand_trees  = config["_debug"]["_num_rand_trees"]
num_parsimony_trees = config["_debug"]["_num_parsimony_trees"]

pars_seeds = range(num_pars_trees)
rand_seeds = range(num_pars_trees, num_pars_trees + num_rand_trees)
# Parsimonator requires seeds greater than 1
parsimony_seeds = range(1, num_parsimony_trees + 1)

# TODO: resolve duplicate names
msa_paths = config["msa_paths"]

# This assumes, that each msa
msa_names = [os.path.split(pth)[1] for pth in msa_paths]

msas = dict(zip(msa_names, msa_paths))

outdir = config["outdir"]
db_path = outdir + "{msa}/"
output_files_dir = outdir + "{msa}/output_files/"

# File paths for RAxML-NG files
output_files_raxmlng_dir = output_files_dir + "raxmlng/"
# tree inference
raxmlng_tree_inference_dir = output_files_raxmlng_dir + "inference/"
raxmlng_tree_inference_prefix_pars = raxmlng_tree_inference_dir + "pars_{seed}"
raxmlng_tree_inference_prefix_rand = raxmlng_tree_inference_dir + "rand_{seed}"
# tree evaluation
raxmlng_tree_eval_dir = output_files_raxmlng_dir + "evaluation/"
raxmlng_tree_eval_prefix_pars = raxmlng_tree_eval_dir + "pars_{seed}"
raxmlng_tree_eval_prefix_rand = raxmlng_tree_eval_dir + "rand_{seed}"

# File paths for IQ-Tree files
output_files_iqtree_dir = output_files_dir + "iqtree/"

# File paths for parsimony trees
output_files_parsimony_trees = output_files_dir + "parsimonator/"
parsimonator_tree_file_name = output_files_parsimony_trees + "RAxML_parsimonyTree.seed_{seed}.tree"
parsimonator_log_file_name = output_files_parsimony_trees + "RAxML_parsimonyTree.seed_{seed}.log"

rule all:
    input:
        expand(f"{db_path}data.sqlite3", msa=msa_names)


include: "rules/raxmlng_tree_inference.smk"
include: "rules/raxmlng_tree_evaluation.smk"
include: "rules/collect_data.smk"
include: "rules/raxmlng_rfdistance.smk"
include: "rules/iqtree_significance_tests.smk"
include: "rules/msa_features.smk"
include: "rules/parsimonator.smk"
include: "rules/save_data.smk"