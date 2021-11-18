import os

configfile: "config.yaml"

raxmlng_command = config["software"]["raxml-ng"]["command"]

num_pars_trees = config["_debug"]["_num_pars_trees"]
num_rand_trees  = config["_debug"]["_num_rand_trees"]

pars_seeds = range(num_pars_trees)
rand_seeds = range(num_pars_trees, num_pars_trees + num_rand_trees)

# TODO: resolve duplicate names
msa_paths = config["msa_paths"]

# This assumes, that each msa
msa_names = [os.path.split(pth)[1] for pth in msa_paths]

msas = dict(zip(msa_names, msa_paths))

outdir = config["outdir"]
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

rule all:
    input:
        expand(f"{raxmlng_tree_inference_dir}inference.raxml.rfDistances", msa=msa_names),
        expand(f"{raxmlng_tree_eval_dir}eval.raxml.rfDistances", msa=msa_names)


include: "rules/raxmlng_tree_inference.smk"
include: "rules/raxmlng_tree_evaluation.smk"
include: "rules/collect_data.smk"
include: "rules/raxmlng_rfdistance.smk"