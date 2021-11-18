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
        f"{raxmlng_tree_inference_dir}.raxml.rfDistances.snakelog",
    shell:
        "{raxmlng_command} "
        "--rfdist "
        "--tree {input.all_search_trees} "
        "--prefix {params.prefix} "
        ">> {output.rfDist_log} "
