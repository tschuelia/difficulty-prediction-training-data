rule raxmlng_pars_tree:
    """
    Rule that infers a single tree based on a parsimony starting tree using RAxML-NG.
    """
    output:
        raxml_best_tree     = f"{raxmlng_tree_inference_prefix_pars}.raxml.bestTree",
        raxml_starting_tree = f"{raxmlng_tree_inference_prefix_pars}.raxml.startTree",
        raxml_best_model    = f"{raxmlng_tree_inference_prefix_pars}.raxml.bestModel",
        raxml_log           = f"{raxmlng_tree_inference_prefix_pars}.raxml.inference.log",
    params:
        prefix  = raxmlng_tree_inference_prefix_pars,
        msa     = lambda wildcards: msas[wildcards.msa],
        model   = lambda wildcards: raxmlng_models[wildcards.msa],
        threads = config["software"]["raxml-ng"]["threads"]
    log:
        f"{raxmlng_tree_inference_prefix_pars}.snakelog",
    shell:
        "{raxmlng_command} "
        "--msa {params.msa} "
        "--model {params.model} "
        "--prefix {params.prefix} "
        "--seed {wildcards.seed} "
        "--threads {params.threads} "
        "--tree pars{{1}} "
        "> {output.raxml_log} "


rule raxmlng_rand_tree:
    """
    Rule that infers a single tree based on a random starting tree using RAxML-NG.
    """
    output:
        raxml_best_tree     = f"{raxmlng_tree_inference_prefix_rand}.raxml.bestTree",
        raxml_best_model    = f"{raxmlng_tree_inference_prefix_rand}.raxml.bestModel",
        raxml_log           = f"{raxmlng_tree_inference_prefix_rand}.raxml.inference.log",
    params:
        prefix  = raxmlng_tree_inference_prefix_rand,
        msa     = lambda wildcards: msas[wildcards.msa],
        model   = lambda wildcards: raxmlng_models[wildcards.msa],
        threads = config["software"]["raxml-ng"]["threads"]
    log:
        f"{raxmlng_tree_inference_prefix_rand}.snakelog",
    shell:
        "{raxmlng_command} "
        "--msa {params.msa} "
        "--model {params.model} "
        "--prefix {params.prefix} "
        "--seed {wildcards.seed} "
        "--threads {params.threads} "
        "--tree rand{{1}} "
        "> {output.raxml_log} "
