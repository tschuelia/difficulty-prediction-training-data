rule raxmlng_rfdistance_search_trees:
    """
    Rule that computes the RF-Distances between all search trees using RAxML-NG.
    """
    input:
        all_search_trees = rules.collect_search_trees.output.all_search_trees
    output:
        rfDist      = f"{raxmlng_tree_inference_dir}inference.raxml.rfDistances",
        rfDist_log  = f"{raxmlng_tree_inference_dir}inference.raxml.rfDistances.log",
    params:
        prefix = f"{raxmlng_tree_inference_dir}inference"
    log:
        f"{raxmlng_tree_inference_dir}inference.raxml.rfDistances.snakelog",
    shell:
        "{raxmlng_command} "
        "--rfdist "
        "--tree {input.all_search_trees} "
        "--prefix {params.prefix} "
        ">> {output.rfDist_log} "


rule raxmlng_rfdistance_eval_trees:
    """
    Rule that computes the RF-Distances between all eval trees using RAxML-NG.
    """
    input:
        all_eval_trees = rules.collect_eval_trees.output.all_eval_trees
    output:
        rfDist      = f"{raxmlng_tree_eval_dir}eval.raxml.rfDistances",
        rfDist_log  = f"{raxmlng_tree_eval_dir}eval.raxml.rfDistances.log",
    params:
        prefix = f"{raxmlng_tree_eval_dir}eval"
    log:
        f"{raxmlng_tree_eval_dir}eval.raxml.rfDistances.snakelog",
    shell:
        "{raxmlng_command} "
        "--rfdist "
        "--tree {input.all_eval_trees} "
        "--prefix {params.prefix} "
        ">> {output.rfDist_log} "

