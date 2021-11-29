import json
import numpy as np
import pickle
import uuid

from database import *
from iqtree_statstest_parser import get_iqtree_results, get_iqtree_results_for_eval_tree_str
from raxml_parser import (
    get_raxmlng_abs_rf_distance,
    get_raxmlng_num_unique_topos,
    get_all_raxmlng_llhs,
    get_raxmlng_llh,
    get_raxmlng_elapsed_time
)

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

# tree search
pars_search_trees = snakemake.input.pars_search_trees
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

# plausible
iqtree_results = get_iqtree_results(snakemake.input.iqtree_results)
with open(snakemake.input.clusters, "rb") as f:
    clusters = pickle.load(f)

# msa features
with open(snakemake.input.msa_features) as f:
    msa_features = json.load(f)

llhs_search = get_all_raxmlng_llhs(search_logs_collected)
llhs_eval = get_all_raxmlng_llhs(eval_logs_collected)

dataset_dbobj = Dataset.create(
    uuid=uuid.uuid4().hex,
    verbose_name=dataset_name,

    # Label features
    num_searches=len(pars_search_trees) + len(rand_search_trees),

    avg_rfdist_search=get_raxmlng_abs_rf_distance(search_rfdistance),
    num_topos_search=get_raxmlng_num_unique_topos(search_rfdistance),
    mean_llh_search=np.mean(llhs_search),
    std_llh_search=np.std(llhs_search),

    avg_rfdist_eval=get_raxmlng_abs_rf_distance(eval_rfdistance),
    num_topos_eval=get_raxmlng_num_unique_topos(eval_rfdistance),
    mean_llh_eval=np.mean(llhs_eval),
    std_llh_eval=np.std(llhs_eval),

    avg_rfdist_plausible=None,
    num_topos_plausible=None,
    mean_llh_plausible=None,
    std_llh_plausible=None,
    num_trees_plausible=None,
    proportion_plausible=None,

    # Single inference features
    num_slow_spr_rounds=None,
    num_fast_spr_rounds=None,
    parsimony_score_starting_tree=None,
    parsimony_score_final_tree=None,
    llh_starting_tree=None,
    llh_final_tree=None,
    rfdistance_starting_final=None,
    llh_difference_starting_final=None,
    eq_frequencies_final=None,
    substitution_rates_final=None,
    average_branch_length_final=None,
    std_branch_length_final=None,
    total_branch_length_final=None,
    minimum_branch_length_final=None,
    maximum_branch_length_final=None,
    newick_starting=None,
    newick_final=None,

    # MSA Features
    num_taxa=msa_features["taxa"],
    num_sites=msa_features["sites"],
    num_patters=msa_features["patterns"],
    proportion_gaps=msa_features["gaps"],
    proportion_invariant=msa_features["invariant"],
    entropy=msa_features["entropy"],
    column_entropies=msa_features["column_entropies"],
    bollback=msa_features["bollback"],
    treelikeness=msa_features["treelikeness"],
    char_frequencies=msa_features["char_frequencies"],
)


def save_raxmlng_tree(search_trees, search_logs, eval_trees, eval_logs, starting_type):
    for (search_tree, search_log, eval_tree, eval_log) in zip(search_trees, search_logs, eval_trees, eval_logs):
        newick_eval = open(eval_tree).readline()
        statstest_results, cluster_id = get_iqtree_results_for_eval_tree_str(iqtree_results, newick_eval, clusters)
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

            # Plausible trees
            plausible=statstest_results["plausible"],
            cluster_id=cluster_id,

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


save_raxmlng_tree(pars_search_trees, pars_search_logs, pars_eval_trees, pars_eval_logs, "parsimony")
save_raxmlng_tree(rand_search_trees, rand_search_logs, rand_eval_trees, rand_eval_logs, "random")