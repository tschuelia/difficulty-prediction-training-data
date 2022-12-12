rule reevaluate_raxml_pars_tree:
    """
    Rule that re-evaluates the given parsimony search tree.
    """
    input:
        best_tree_of_run    = raxmlng_tree_inference_prefix_pars.with_suffix(".raxml.bestTree")
    output:
        log         = raxmlng_tree_eval_prefix_pars.with_suffix(".raxml.log"),
        best_tree   = raxmlng_tree_eval_prefix_pars.with_suffix(".raxml.bestTree"),
        eval_log    = raxmlng_tree_eval_prefix_pars.with_suffix(".raxml.eval.log")
    params:
        prefix  = raxmlng_tree_eval_prefix_pars,
        msa     = lambda wildcards: msas[wildcards.msa],
        model   = lambda wildcards: raxmlng_models[wildcards.msa],
        threads = config["software"]["raxml-ng"]["threads"]
    log:
        raxmlng_tree_eval_prefix_pars.with_suffix(".snakelog")
    shell:
        "{raxmlng_command} "
        "--eval "
        "--tree {input.best_tree_of_run} "
        "--msa {params.msa} "
        "--model {params.model} "
        "--prefix {params.prefix} "
        "--threads {params.threads} "
        "--seed 0 "
        "> {output.eval_log} "


rule reevaluate_raxml_rand_tree:
    """
    Rule that re-evaluates the given random search tree.
    """
    input:
        best_tree_of_run    = raxmlng_tree_inference_prefix_rand.with_suffix(".raxml.bestTree")
    output:
        log         = raxmlng_tree_eval_prefix_rand.with_suffix(".raxml.log"),
        best_tree   = raxmlng_tree_eval_prefix_rand.with_suffix(".raxml.bestTree"),
        eval_log    = raxmlng_tree_eval_prefix_rand.with_suffix(".raxml.eval.log")
    params:
        prefix  = raxmlng_tree_eval_prefix_rand,
        msa     = lambda wildcards: msas[wildcards.msa],
        model   = lambda wildcards: raxmlng_models[wildcards.msa],
        threads = config["software"]["raxml-ng"]["threads"]
    log:
        raxmlng_tree_eval_prefix_rand.with_suffix(".snakelog")
    shell:
        "{raxmlng_command} "
        "--eval "
        "--tree {input.best_tree_of_run} "
        "--msa {params.msa} "
        "--model {params.model} "
        "--prefix {params.prefix} "
        "--threads {params.threads} "
        "--seed 0 "
        "> {output.eval_log} "