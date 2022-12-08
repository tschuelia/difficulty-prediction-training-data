rule compute_tree_metrics:
    input:
        eval_trees = rules.collect_eval_trees.output.all_eval_trees
    output:
        tree_metrics =  f"{output_files_dir}tree_metrics.json"
    script:
        "scripts/collect_tree_metrics.py"
