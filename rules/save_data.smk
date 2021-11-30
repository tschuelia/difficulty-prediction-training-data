rule save_data:
    input:
        # Tree seach tree files and logs
        pars_search_trees   = expand(raxmlng_tree_inference_prefix_pars + ".raxml.bestTree",seed=pars_seeds,allow_missing=True),
        pars_starting_trees = expand(raxmlng_tree_inference_prefix_pars + ".raxml.startTree",seed=pars_seeds,allow_missing=True),
        pars_search_logs    = expand(raxmlng_tree_inference_prefix_pars + ".raxml.inference.log",seed=pars_seeds,allow_missing=True),
        rand_search_trees   = expand(raxmlng_tree_inference_prefix_rand + ".raxml.bestTree", seed=rand_seeds, allow_missing=True),
        rand_search_logs    = expand(raxmlng_tree_inference_prefix_rand + ".raxml.inference.log",seed=rand_seeds,allow_missing=True),
        search_logs_collected = f"{raxmlng_tree_inference_dir}AllSearchLogs.log",

        # Tree search tree RFDistance logs
        search_rfdistance = f"{raxmlng_tree_inference_dir}inference.raxml.rfDistances.log",

        # Eval tree files and logs
        pars_eval_trees = expand(raxmlng_tree_eval_prefix_pars + ".raxml.bestTree", seed=pars_seeds, allow_missing=True),
        pars_eval_logs  = expand(raxmlng_tree_eval_prefix_pars + ".raxml.eval.log",seed=pars_seeds,allow_missing=True),
        rand_eval_trees = expand(raxmlng_tree_eval_prefix_rand + ".raxml.bestTree", seed=rand_seeds, allow_missing=True),
        rand_eval_logs  = expand(raxmlng_tree_eval_prefix_rand + ".raxml.eval.log",seed=rand_seeds,allow_missing=True),
        eval_logs_collected = f"{raxmlng_tree_eval_dir}AllEvalLogs.log",

        # Eval tree RFDistance logs
        eval_rfdistance = f"{raxmlng_tree_eval_dir}eval.raxml.rfDistances.log",

        # Plausible tree RFDistance logs
        plausible_rfdistance = f"{raxmlng_tree_eval_dir}plausible.raxml.rfDistances.log",
        plausible_trees_collected = f"{raxmlng_tree_eval_dir}AllPlausibleTrees.trees",

        # IQ-Tree significance test results and clusters
        iqtree_results  = f"{output_files_iqtree_dir}significance.iqtree",
        clusters        = f"{output_files_iqtree_dir}filteredEvalTrees.clusters.pkl",

        # MSA Features
        msa_features = f"{output_files_dir}msa_features.json",
    output:
        database = f"{db_path}data.sqlite3"

    params:
        raxmlng_command = raxmlng_command
    script:
        "scripts/save_data.py"