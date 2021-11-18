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
print(msa_names)

msas = dict(zip(msa_names, msa_paths))
print(msas)

outdir = config["outdir"]
output_files_dir = outdir + "{msa}/output_files/"
raxmlng_tree_inference_dir = output_files_dir + "raxmlng/inference/"
raxmlng_tree_inference_prefix_pars = raxmlng_tree_inference_dir + "pars_{seed}"
raxmlng_tree_inference_prefix_rand = raxmlng_tree_inference_dir + "rand_{seed}"

raxmlng_tree_eval_dir = output_files_dir + "raxmlng/evaluation/"
raxmlng_tree_eval_prefix_pars = raxmlng_tree_eval_dir + "pars_{seed}"
raxmlng_tree_eval_prefix_rand = raxmlng_tree_eval_dir + "rand_{seed}"

rule all:
    input:
        raxmlng_pars_search_trees = expand(raxmlng_tree_inference_prefix_pars + ".raxml.bestTree", msa=msa_names, seed=pars_seeds),
        raxmlng_rand_search_trees = expand(raxmlng_tree_inference_prefix_rand + ".raxml.bestTree", msa=msa_names, seed=pars_seeds)


include: "rules/raxmlng_tree_inference.smk"