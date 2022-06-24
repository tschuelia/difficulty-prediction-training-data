"""
This file contains lists and explanations for the features in the training_data.parquet file
"""
LABEL = "difficult"

"""
Subset of features we present in the main Pythia paper that we trained the final predictor on
"""

FINAL_FEATURES = [
    "num_patterns/num_taxa",
    "num_sites/num_taxa",
    "proportion_gaps",
    "proportion_invariant",
    "entropy",
    "bollback",
    "avg_rfdist_parsimony",
    "proportion_unique_topos_parsimony",
]

"""
All features we present in the Pythia supplementary information  we used for our prediction experiments 
"""
ALL_FEATURES = [
    "num_slow_spr_rounds",
    "num_fast_spr_rounds",
    "rfdistance_starting_final",
    "average_branch_length_final",
    "std_branch_length_final",
    "total_branch_length_final",
    "minimum_branch_length_final",
    "maximum_branch_length_final",
    "num_taxa",
    "num_sites",
    "num_patterns",
    "proportion_gaps",
    "proportion_invariant",
    "entropy",
    "column_entropies",
    "bollback",
    "treelikeness",
    "avg_rfdist_parsimony",
    "num_topos_parsimony",
    "num_topos_parsimony/num_trees_parsimony",
    "num_patterns/num_taxa",
    "num_sites/num_taxa",
    "proportion_unique_topos_parsimony",
]

"""
Additional features we compute using the presented methods, but decided against experimenting with them
"""

ADDITIONAL_FEATURES = [
    "llh_starting_tree",
    "llh_final_tree",
    "llh_difference_starting_final",
    "rate_heterogeneity_final",
    "eq_frequencies_final",
    "substitution_rates_final",
    "newick_starting",
    "newick_final",
    "char_frequencies",
    "freq_a",
    "freq_t",
    "freq_c",
    "freq_g",
    "mean_parsimony_score",
    "std_parsimony_score",
]

"""
Features used for label generation.
CAREFUL: if you use the LABEL column, you should not use these features to train a predictor. 
These features were used to quantify the difficulty and are therefore directly correlated to the difficulty.
"""

LABEL_GENERATION_FEATURES = [
    "avg_rfdist_eval",
    "num_topos_eval",
    "mean_llh_eval",
    "std_llh_eval",
    "num_topos_eval/num_trees_eval",
    "avg_rfdist_search",
    "num_topos_search",
    "mean_llh_search",
    "std_llh_search",
    "num_topos_search/num_trees_search",
    "avg_rfdist_plausible",
    "num_topos_plausible",
    "mean_llh_plausible",
    "std_llh_plausible",
    "num_trees_plausible",
    "proportion_plausible",
    "num_topos_plausible/num_trees_plausible",
]
