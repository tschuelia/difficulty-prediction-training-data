import os
import sys
import pathlib

sys.path.append("rules/scripts")
from pypythia.msa import MSA

configfile: "config.yaml"

raxmlng_command = config["software"]["raxml-ng"]["command"]
iqtree_command = config["software"]["iqtree"]["command"]

num_pars_trees = config["_debug"]["_num_pars_trees"]
num_rand_trees  = config["_debug"]["_num_rand_trees"]
num_parsimony_trees = config["_debug"]["_num_parsimony_trees"]

pars_seeds = range(num_pars_trees)
rand_seeds = range(num_pars_trees, num_pars_trees + num_rand_trees)
parsimony_seeds = range(1, num_parsimony_trees + 1)

# TODO: resolve duplicate names
msa_paths = config["msa_paths"]
part_paths_raxmlng = []
part_paths_iqtree = []
partitioned = False

if isinstance(msa_paths[0], list):
    # in this case the MSAs are partitioned
    msa_paths, part_paths_raxmlng, part_paths_iqtree = zip(*msa_paths)
    partitioned = True

# This assumes, that each msa
msa_names = [os.path.split(pth)[1] for pth in msa_paths]
msas = dict(zip(msa_names, msa_paths))

data_types = {}
for msa, name in zip(msa_paths, msa_names):
    msa = MSA(msa)
    data_types[name] = msa.data_type

if partitioned:
    raxmlng_models = dict(list(zip(msa_names, part_paths_raxmlng)))
    iqtree_models = dict(list(zip(msa_names, part_paths_iqtree)))
else:
    # infer the data type for each MSA
    raxmlng_models = []
    iqtree_models = []
    for name, msa in msas.items():
        msa = MSA(msa)
        raxmlng_model = msa.get_raxmlng_model()
        raxmlng_models.append((name, raxmlng_model))

        if msa.data_type == "MORPH":
            iqtree_models.append((name, "MK"))
        else:
            iqtree_models.append((name, f"{raxmlng_model}4+FO"))

    raxmlng_models = dict(raxmlng_models)
    iqtree_models = dict(iqtree_models)



outdir = pathlib.Path(config["outdir"])
db_path = outdir / "{msa}/"
output_files_dir = outdir / "{msa}" / "output_files"

# File paths for RAxML-NG files
output_files_raxmlng_dir = output_files_dir / "raxmlng"
# tree inference
raxmlng_tree_inference_dir = output_files_raxmlng_dir / "inference"
raxmlng_tree_inference_prefix_pars = raxmlng_tree_inference_dir / "pars_{seed}"
raxmlng_tree_inference_prefix_rand = raxmlng_tree_inference_dir / "rand_{seed}"
# tree evaluation
raxmlng_tree_eval_dir = output_files_raxmlng_dir / "evaluation"
raxmlng_tree_eval_prefix_pars = raxmlng_tree_eval_dir / "pars_{seed}"
raxmlng_tree_eval_prefix_rand = raxmlng_tree_eval_dir / "rand_{seed}"
# bootstrapping
raxmlng_bootstrap_prefix =  output_files_raxmlng_dir / "bootstrap"

# File paths for IQ-Tree files
output_files_iqtree_dir = output_files_dir / "iqtree"

# File paths for parsimony trees
output_files_parsimony_trees = output_files_dir / "parsimony"
parsimony_tree_file_name = output_files_parsimony_trees / "seed_{seed}.raxml.startTree"
parsimony_log_file_name = output_files_parsimony_trees / "seed_{seed}.raxml.log"


rule all:
    input:
        expand(db_path / "training_data.parquet", msa=msa_names),
        expand(db_path / "raxmlng_tree_data.parquet", msa=msa_names)


include: "rules/raxmlng_tree_inference.smk"
include: "rules/raxmlng_tree_evaluation.smk"
include: "rules/collect_data.smk"
include: "rules/raxmlng_rfdistance.smk"
include: "rules/iqtree_significance_tests.smk"
include: "rules/msa_features.smk"
include: "rules/parsimony.smk"
include: "rules/save_data.smk"
include: "rules/tree_metrics.smk"