rule reevaluate_raxml_pars_tree:
    """
    Rule that re-evaluates the given parsimony search tree.
    """
    input:
        best_tree_of_run    = f"{raxmlng_tree_inference_prefix_pars}.raxml.bestTree"
    output:
        log         = f"{raxmlng_tree_eval_prefix_pars}.raxml.log",
        best_tree   = f"{raxmlng_tree_eval_prefix_pars}.raxml.bestTree",
        eval_log    = f"{raxmlng_tree_eval_prefix_pars}.raxml.eval.log",
    params:
        prefix  = raxmlng_tree_eval_prefix_pars,
        msa     = lambda wildcards: msas[wildcards.msa],
        model   = lambda wildcards: raxmlng_model[wildcards.msa] if partitioned else raxmlng_model
    log:
        f"{raxmlng_tree_eval_prefix_pars}.snakelog"
    shell:
        "{raxmlng_command} "
        "--eval "
        "--tree {input.best_tree_of_run} "
        "--msa {params.msa} "
        "--model {params.model} "
        "--prefix {params.prefix} "
        "--seed 0 "
        "> {output.eval_log} "


rule reevaluate_raxml_rand_tree:
    """
    Rule that re-evaluates the given random search tree.
    """
    input:
        best_tree_of_run    = f"{raxmlng_tree_inference_prefix_rand}.raxml.bestTree"
    output:
        log         = f"{raxmlng_tree_eval_prefix_rand}.raxml.log",
        best_tree   = f"{raxmlng_tree_eval_prefix_rand}.raxml.bestTree",
        eval_log    = f"{raxmlng_tree_eval_prefix_rand}.raxml.eval.log",
    params:
        prefix  = raxmlng_tree_eval_prefix_rand,
        msa     = lambda wildcards: msas[wildcards.msa],
        model   = lambda wildcards: raxmlng_model[wildcards.msa] if partitioned else raxmlng_model
    log:
        f"{raxmlng_tree_eval_prefix_rand}.snakelog"
    shell:
        "{raxmlng_command} "
        "--eval "
        "--tree {input.best_tree_of_run} "
        "--msa {params.msa} "
        "--model {params.model} "
        "--prefix {params.prefix} "
        "--seed 0 "
        "> {output.eval_log} "