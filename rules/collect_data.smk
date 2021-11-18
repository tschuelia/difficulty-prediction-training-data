rule collect_search_trees:
    """
    Rule that collects all search trees for one dataset in one file.
    """
    input:
        raxmlng_pars_search_trees = expand(raxmlng_tree_inference_prefix_pars + ".raxml.bestTree", seed=pars_seeds, allow_missing=True),
        raxmlng_rand_search_trees = expand(raxmlng_tree_inference_prefix_rand + ".raxml.bestTree", seed=pars_seeds, allow_missing=True)
    output:
        all_search_trees = f"{output_files_raxmlng_dir}AllSearchTrees.trees"
    shell:
        "cat {input.raxmlng_pars_search_trees} {input.raxmlng_rand_search_trees} > {output.all_search_trees}"
