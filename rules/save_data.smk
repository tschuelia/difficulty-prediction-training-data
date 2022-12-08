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
        best_eval_tree      = f"{raxmlng_tree_eval_dir}BestEvalTree.tree",
        eval_logs_collected = f"{raxmlng_tree_eval_dir}AllEvalLogs.log",

        # Eval tree RFDistance logs
        eval_rfdistance = f"{raxmlng_tree_eval_dir}eval.raxml.rfDistances.log",

        # Plausible tree RFDistance logs
        plausible_rfdistance = f"{raxmlng_tree_eval_dir}plausible.raxml.rfDistances.log",
        plausible_trees_collected = f"{raxmlng_tree_eval_dir}AllPlausibleTrees.trees",

        # IQ-Tree significance test results
        iqtree_results  = f"{output_files_iqtree_dir}significance.iqtree",

        # MSA Features
        msa_features = f"{output_files_dir}msa_features.json",

        # Parsimony Trees and logs
        parsimony_trees = f"{output_files_parsimony_trees}AllParsimonyTrees.trees",
        parsimony_logs = f"{output_files_parsimony_trees}AllParsimonyLogs.log",
        parsimony_rfdistance = f"{output_files_parsimony_trees}parsimony.raxml.rfDistances.log",
    output:
        database = "{msa}_data.sqlite3"
    params:
        raxmlng_command = raxmlng_command,
        msa             = lambda wildcards: msas[wildcards.msa],
    script:
        "scripts/save_data.py"


rule move_db:
    # due to an issue with our lab webservers, I cannot directly create the database on the mounted fs
    # therefore I creat it in the current workdir and then move it to the mounted fs
    input:
        "{msa}_data.sqlite3"
    output:
        database = f"{db_path}data.sqlite3"
    shell:
        "mv {input} {output}"


rule database_to_training_dataframe:
    input:
        database = rules.move_db.output.database,
    output:
        training_data = f"{db_path}training_data.parquet",
        raxmlng_tree_data = f"{db_path}raxmlng_tree_data.parquet",
    params:
        num_pars_trees = num_pars_trees,
        num_rand_trees = num_rand_trees,
        num_parsimony_trees = num_parsimony_trees
    script:
        "scripts/database_to_dataframe.py"