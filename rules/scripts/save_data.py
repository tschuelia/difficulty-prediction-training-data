import json
import numpy as np
import uuid

from database import *
from iqtree_statstest_parser import get_iqtree_results
from raxmlng_parser import (
    get_all_raxmlng_llhs,
    get_raxmlng_llh,
    get_raxmlng_elapsed_time,
    get_raxmlng_starting_llh,
    get_raxmlng_num_spr_rounds,
    rel_rfdistance_starting_final,
    get_model_parameter_estimates,
)
from utils import read_file_contents

from pypythia.raxmlng_parser import get_raxmlng_rfdist_results
from pypythia.msa import MSA

db.init(snakemake.output.database)
db.connect()
db.create_tables(
    [
        Dataset,
        RaxmlNGTree,
        ParsimonyTree
    ]
)
dataset_name = snakemake.wildcards.msa
raxmlng_command = snakemake.params.raxmlng_command

# tree search
pars_search_trees = snakemake.input.pars_search_trees
pars_starting_trees = snakemake.input.pars_starting_trees
pars_search_logs = snakemake.input.pars_search_logs
rand_search_trees = snakemake.input.rand_search_trees
rand_search_logs = snakemake.input.rand_search_logs
search_logs_collected = snakemake.input.search_logs_collected
search_rfdistance = snakemake.input.search_rfdistance

# eval
pars_eval_trees = snakemake.input.pars_eval_trees
pars_eval_logs = snakemake.input.pars_eval_logs
rand_eval_trees = snakemake.input.rand_eval_trees
rand_eval_logs = snakemake.input.rand_eval_logs
eval_logs_collected = snakemake.input.eval_logs_collected
eval_rfdistance = snakemake.input.eval_rfdistance
best_eval_tree = open(snakemake.input.best_eval_tree).readline().strip()

# plausible
plausible_rfdistance = snakemake.input.plausible_rfdistance
plausible_trees_collected = read_file_contents(snakemake.input.plausible_trees_collected)
iqtree_results = get_iqtree_results(snakemake.input.iqtree_results)

# msa features
with open(snakemake.input.msa_features) as f:
    msa_features = json.load(f)

# parsimony trees
parsimony_trees = snakemake.input.parsimony_trees
parsimony_logs = snakemake.input.parsimony_logs
parsimony_rfdistance = snakemake.input.parsimony_rfdistance

llhs_search = get_all_raxmlng_llhs(search_logs_collected)
llhs_eval = get_all_raxmlng_llhs(eval_logs_collected)

num_searches = len(pars_search_trees) + len(rand_search_trees)
data_type = MSA(snakemake.params.msa).data_type

# for the starting tree features, we simply take the first parsimony tree inference
single_tree = pars_search_trees[0]
single_tree_log = pars_search_logs[0]
single_tree_starting = pars_starting_trees[0]

slow_spr, fast_spr = get_raxmlng_num_spr_rounds(single_tree_log)
starting_llh = get_raxmlng_starting_llh(single_tree_log)
final_llh = get_raxmlng_llh(single_tree_log)
newick_starting = open(single_tree_starting).readline()
newick_final = open(single_tree).readline()
rate_het, base_freq, subst_rates = get_model_parameter_estimates(single_tree_log)

num_topos_search, avg_rfdist_search, _ = get_raxmlng_rfdist_results(search_rfdistance)
num_topos_eval, avg_rfdist_eval, _ = get_raxmlng_rfdist_results(eval_rfdistance)
num_topos_plausible, avg_rfdist_plausible, _ = get_raxmlng_rfdist_results(plausible_rfdistance)
num_topos_parsimony, avg_rfdist_parsimony, _ = get_raxmlng_rfdist_results(parsimony_rfdistance)

# fmt: off
dataset_dbobj = Dataset.create(
    uuid        = uuid.uuid4().hex,
    verbose_name= dataset_name,
    data_type = data_type,

    # Label features
    num_searches=num_searches,

    avg_rfdist_search   = avg_rfdist_search,
    num_topos_search    = num_topos_search,
    mean_llh_search     = np.mean(llhs_search),
    std_llh_search      = np.std(llhs_search),

    avg_rfdist_eval = avg_rfdist_search,
    num_topos_eval  = num_topos_eval,
    mean_llh_eval   = np.mean(llhs_eval),
    std_llh_eval    = np.std(llhs_eval),

    avg_rfdist_plausible    = avg_rfdist_plausible,
    num_topos_plausible     = num_topos_plausible,
    # we will update this information after inserting the trees to the database
    mean_llh_plausible  = None,
    std_llh_plausible   = None,
    num_trees_plausible = None,
    proportion_plausible= None,

    # Single inference features
    num_slow_spr_rounds             = slow_spr,
    num_fast_spr_rounds             = fast_spr,
    llh_starting_tree               = starting_llh,
    llh_final_tree                  = final_llh,
    rfdistance_starting_final       = rel_rfdistance_starting_final(newick_starting, newick_final, raxmlng_command),
    llh_difference_starting_final   = final_llh - starting_llh,
    rate_heterogeneity_final        = rate_het,
    eq_frequencies_final            = base_freq,
    substitution_rates_final        = subst_rates,
    newick_starting                 = newick_starting,
    newick_final                    = newick_final,

    # MSA Features
    num_taxa                = msa_features["taxa"],
    num_sites               = msa_features["sites"],
    num_patterns            = msa_features["patterns"],
    proportion_gaps         = msa_features["gaps"],
    proportion_invariant    = msa_features["invariant"],
    entropy                 = msa_features["entropy"],
    column_entropies        = msa_features["column_entropies"],
    bollback                = msa_features["bollback"],
    treelikeness            = msa_features["treelikeness"],

    # Parsimony Trees Features
    avg_rfdist_parsimony    = avg_rfdist_parsimony,
    num_topos_parsimony     = num_topos_parsimony,
)


# fmt: on
def save_raxmlng_tree(search_trees, search_logs, eval_trees, eval_logs, starting_type):
    plausible_llhs = []

    for (search_tree, search_log, eval_tree, eval_log, statstest_results) in zip(search_trees, search_logs, eval_trees, eval_logs, iqtree_results):
        newick_eval = open(eval_tree).readline().strip()
        tests = statstest_results["tests"]

        RaxmlNGTree.create(
            dataset=dataset_dbobj,
            dataset_uuid=dataset_dbobj.uuid,
            uuid=uuid.uuid4().hex,

            # Search trees
            starting_type=starting_type,
            newick_search=open(search_tree).readline(),
            llh_search=get_raxmlng_llh(search_log),
            compute_time_search=get_raxmlng_elapsed_time(search_log),

            # Eval trees
            newick_eval=newick_eval,
            llh_eval=get_raxmlng_llh(eval_log),
            compute_time_eval=get_raxmlng_elapsed_time(eval_log),
            is_best=newick_eval == best_eval_tree,

            # Plausible trees
            plausible=statstest_results["plausible"],
            bpRell=tests["bp-RELL"]["score"],
            bpRell_significant=tests["bp-RELL"]["significant"],
            pKH=tests["p-KH"]["score"],
            pKH_significant=tests["p-KH"]["significant"],
            pSH=tests["p-SH"]["score"],
            pSH_significant=tests["p-SH"]["significant"],
            pWKH=tests["p-WKH"]["score"],
            pWKH_significant=tests["p-WKH"]["significant"],
            pWSH=tests["p-WSH"]["score"],
            pWSH_significant=tests["p-WSH"]["significant"],
            cELW=tests["c-ELW"]["score"],
            cELW_significant=tests["c-ELW"]["significant"],
            pAU=tests["p-AU"]["score"],
            pAU_significant=tests["p-AU"]["significant"],
        )

        if statstest_results["plausible"]:
            plausible_llhs.append(get_raxmlng_llh(eval_log))

    return plausible_llhs

# store the parsimony and random raxml-ng trees in the database
plausible_llhs_pars = save_raxmlng_tree(pars_search_trees, pars_search_logs, pars_eval_trees, pars_eval_logs, "parsimony")
plausible_llhs_rand = save_raxmlng_tree(rand_search_trees, rand_search_logs, rand_eval_trees, rand_eval_logs, "random")

plausible_llhs = plausible_llhs_pars + plausible_llhs_rand
dataset_dbobj.update(
    {
        "mean_llh_plausible": np.mean(plausible_llhs),
        "std_llh_plausible": np.std(plausible_llhs),
        "num_trees_plausible": len(plausible_trees_collected),
        "proportion_plausible": len(plausible_trees_collected) / num_searches,
    }
).execute()

# store the parsimonator parsimony trees in the database
parsimony_trees = open(parsimony_trees).readlines()
parsimony_trees = [tree.strip() for tree in parsimony_trees if tree]


for tree in parsimony_trees:
    ParsimonyTree.create(
        uuid            = uuid.uuid4(),
        dataset         = dataset_dbobj,
        dataset_uuid    = dataset_dbobj.uuid,
        newick_tree     = tree
    )
