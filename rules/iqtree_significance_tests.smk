rule iqtree_significance_tests_on_eval_trees:
    """
    Perfoms all significance tests as implemented in IQ-Tree on the set of eval trees.
    As reference tree for estimating the model parameters, we pass the best tree 
    (i.e. with the highest log-likelihood) of the dataset
    """
    input:
        all_eval_trees  = rules.collect_eval_trees.output.all_eval_trees,
        best_tree       = rules.save_best_eval_tree.output.best_eval_tree
    output:
        summary     = output_files_iqtree_dir / "significance.iqtree",
        iqtree_log  = output_files_iqtree_dir / "significance.iqtree.log",
    params:
        msa         = lambda wildcards: msas[wildcards.msa],
        data_type   = lambda wildcards: data_types[wildcards.msa],
        prefix      = str(output_files_iqtree_dir / "significance"),
        model       = lambda wildcards: iqtree_models[wildcards.msa],
        model_str   = "-p" if partitioned else "-m",
        threads     = config["software"]["iqtree"]["threads"]
    log:
        output_files_iqtree_dir / "significance.iqtree.snakelog",
    run:
        morph = "-st MORPH " if params.data_type == "MORPH" else ""
        shell("{iqtree_command} "
        "-s {params.msa} "
        "{morph} "
        "{params.model_str} {params.model} "
        "-pre {params.prefix} "
        "-z {input.all_eval_trees} "
        "-te {input.best_tree} "
        "-n 0 "
        "-zb 10000 "
        "-zw "
        "-au "
        "-nt {params.threads} "
        "-treediff "
        "-seed 0 "
        "> {log} ")
