rule raxmlng_pars_tree:
    """
    Rule that infers a single tree based on a parsimony starting tree using RAxML-NG.
    """
    output:
        raxml_best_tree     = raxmlng_tree_inference_prefix_pars.with_suffix(".raxml.bestTree"),
        raxml_starting_tree = raxmlng_tree_inference_prefix_pars.with_suffix(".raxml.startTree"),
        raxml_best_model    = raxmlng_tree_inference_prefix_pars.with_suffix(".raxml.bestModel"),
        raxml_log           = raxmlng_tree_inference_prefix_pars.with_suffix(".raxml.inference.log")
    params:
        prefix  = raxmlng_tree_inference_prefix_pars,
        msa     = lambda wildcards: msas[wildcards.msa],
        model   = lambda wildcards: raxmlng_models[wildcards.msa],
        threads = config["software"]["raxml-ng"]["threads"]
    log:
        raxmlng_tree_inference_prefix_pars.with_suffix(".snakelog"),
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
        raxml_best_tree     = raxmlng_tree_inference_prefix_rand.with_suffix(".raxml.bestTree"),
        raxml_best_model    = raxmlng_tree_inference_prefix_rand.with_suffix(".raxml.bestModel"),
        raxml_log           = raxmlng_tree_inference_prefix_rand.with_suffix(".raxml.inference.log")
    params:
        prefix  = raxmlng_tree_inference_prefix_rand,
        msa     = lambda wildcards: msas[wildcards.msa],
        model   = lambda wildcards: raxmlng_models[wildcards.msa],
        threads = config["software"]["raxml-ng"]["threads"]
    log:
        raxmlng_tree_inference_prefix_rand.with_suffix(".snakelog")
    shell:
        "{raxmlng_command} "
        "--msa {params.msa} "
        "--model {params.model} "
        "--prefix {params.prefix} "
        "--seed {wildcards.seed} "
        "--threads {params.threads} "
        "--tree rand{{1}} "
        "> {output.raxml_log} "
