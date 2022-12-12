rule collect_search_trees:
    """
    Rule that collects all search trees for one dataset in one file.
    """
    input:
        raxmlng_pars_search_trees = expand(raxmlng_tree_inference_prefix_pars.with_suffix(".raxml.bestTree"), seed=pars_seeds, allow_missing=True),
        raxmlng_rand_search_trees = expand(raxmlng_tree_inference_prefix_rand.with_suffix(".raxml.bestTree"), seed=rand_seeds, allow_missing=True)
    output:
        all_search_trees = raxmlng_tree_inference_dir / "AllSearchTrees.trees"
    shell:
        "cat {input.raxmlng_pars_search_trees} {input.raxmlng_rand_search_trees} > {output.all_search_trees}"


rule collect_search_logs:
    """
    Rule that collects all search logs for one dataset in one file.
    """
    input:
        raxmlng_pars_search_logs = expand(raxmlng_tree_inference_prefix_pars.with_suffix(".raxml.inference.log"), seed=pars_seeds, allow_missing=True),
        raxmlng_rand_search_logs = expand(raxmlng_tree_inference_prefix_rand.with_suffix(".raxml.inference.log"), seed=rand_seeds, allow_missing=True)
    output:
        all_search_logs = raxmlng_tree_inference_dir / "AllSearchLogs.log"
    shell:
        "cat {input.raxmlng_pars_search_logs} {input.raxmlng_rand_search_logs} > {output.all_search_logs}"



rule collect_eval_trees:
    """
    Rule that collects all eval trees for one dataset in one file.
    """
    input:
        raxmlng_pars_eval_trees = expand(raxmlng_tree_eval_prefix_pars.with_suffix(".raxml.bestTree"), seed=pars_seeds, allow_missing=True),
        raxmlng_rand_eval_trees = expand(raxmlng_tree_eval_prefix_rand.with_suffix(".raxml.bestTree"), seed=rand_seeds, allow_missing=True)
    output:
        all_eval_trees = raxmlng_tree_eval_dir / "AllEvalTrees.trees"
    shell:
        "cat {input.raxmlng_pars_eval_trees} {input.raxmlng_rand_eval_trees} > {output.all_eval_trees}"


rule collect_eval_logs:
    """
    Rule that collects all eval logs for one dataset in one file.
    """
    input:
        raxmlng_pars_eval_logs = expand(raxmlng_tree_eval_prefix_pars.with_suffix(".raxml.eval.log"), seed=pars_seeds, allow_missing=True),
        raxmlng_rand_eval_logs = expand(raxmlng_tree_eval_prefix_rand.with_suffix(".raxml.eval.log"), seed=rand_seeds, allow_missing=True)
    output:
        all_eval_logs = raxmlng_tree_eval_dir / "AllEvalLogs.log"
    shell:
        "cat {input.raxmlng_pars_eval_logs} {input.raxmlng_rand_eval_logs} > {output.all_eval_logs}"


rule save_best_eval_tree:
    """
    Rule that saves the best eval tree for the dataset in one file.
    The best tree is the eval tree with the highest log-likelihood score.
    """
    input:
        all_eval_trees = rules.collect_eval_trees.output.all_eval_trees,
        all_eval_logs= rules.collect_eval_logs.output.all_eval_logs,
    output:
        best_eval_tree = raxmlng_tree_eval_dir / "BestEvalTree.tree"
    script:
        "scripts/save_best_eval_tree.py"


rule collect_plausible_trees:
    """
    Rule that collects all plausible trees for one dataset in one file.
    """
    input:
        iqtree_results = output_files_iqtree_dir / "significance.iqtree",
        eval_trees = rules.collect_eval_trees.output.all_eval_trees
    output:
        all_plausible_trees = raxmlng_tree_eval_dir / "AllPlausibleTrees.trees",
    script:
        "scripts/collect_plausible_trees.py"


rule collect_parsimony_trees:
    """
    Rule that collects all parsimony trees inferred with RAxML-NG in one file. 
    """
    input:
        parsimony_trees = expand(parsimony_tree_file_name, seed=parsimony_seeds, allow_missing=True),
    output:
        all_trees = output_files_parsimony_trees / "AllParsimonyTrees.trees"
    shell:
        "cat {input.parsimony_trees} > {output.all_trees}"


rule collect_parsimony_logs:
    """
    Rule that collects all parsimony trees in one file. 
    """
    input:
        parsimony_logs = expand(parsimony_log_file_name, seed=parsimony_seeds, allow_missing=True),
    output:
        all_logs = output_files_parsimony_trees / "AllParsimonyLogs.log"
    shell:
        "cat {input.parsimony_logs} > {output.all_logs}"
