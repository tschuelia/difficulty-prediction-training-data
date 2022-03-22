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


rule raxmlng_rfdistance_plausible_trees:
    """
    Rule that computes the RF-Distances between all plausible trees using RAxML-NG.
    """
    input:
        all_plausible_trees = rules.collect_plausible_trees.output.all_plausible_trees
    output:
        rfDist      = f"{raxmlng_tree_eval_dir}plausible.raxml.rfDistances",
        rfDist_log  = f"{raxmlng_tree_eval_dir}plausible.raxml.rfDistances.log",
    params:
        prefix = f"{raxmlng_tree_eval_dir}plausible"
    log:
        f"{raxmlng_tree_eval_dir}plausible.raxml.rfDistances.snakelog",
    run:
        num_plausible = len(open(input.all_plausible_trees).readlines())
        # we need this distinction because RAxML-NG requires more than one tree in the input file
        # in order to compute the RF Distance
        # but there might be no plausible trees for a given dataset
        if num_plausible <= 1:
            # write 0.0 as RF-Distance in a dummy log
            with open(output.rfDist_log, "w") as f:
                f.write("""
                Number of unique topologies in this tree set: 1
                Average absolute RF distance in this tree set: 0.0
                Average relative RF distance in this tree set: 0.0
                """)

            with open(output.rfDist, "w") as f:
                f.write("0 1 0.0 0.0")
        else:
            shell("{raxmlng_command} --rfdist --tree {input.all_plausible_trees} --prefix {params.prefix} >> {output.rfDist_log}")


rule raxmlng_rfdistance_parsimony_trees:
    """
    Rule that computes the RF-Distances between all parsimony trees inferred with Parsimonator using RAxML-NG.
    """
    input:
        all_parsimony_trees = f"{output_files_parsimony_trees}AllParsimonyTrees.trees",
    output:
        rfDist      = f"{output_files_parsimony_trees}parsimony.raxml.rfDistances",
        rfDist_log  = f"{output_files_parsimony_trees}parsimony.raxml.rfDistances.log",
    params:
        prefix = f"{output_files_parsimony_trees}parsimony"
    log:
        f"{output_files_parsimony_trees}parsimony.raxml.rfDistances.snakelog",
    shell:
        "{raxmlng_command} "
        "--rfdist "
        "--tree {input.all_parsimony_trees} "
        "--prefix {params.prefix} "
        ">> {output.rfDist_log} "
