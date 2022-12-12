def tree_metrics_input():
    input_files = {
        "database": rules.move_db.output.database,
        "eval_trees": rules.collect_eval_trees.output.all_eval_trees,
    }

    if config["bootstrap_based_metrics"]:
        input_files["bootstraps"] = rules.raxmlng_bootstrap.output.bootstraps

    return input_files


rule raxmlng_bootstrap:
    output:
        bootstraps = f"{raxmlng_bootstrap_prefix}.raxml.bootstraps",
        bootstrap_log = f"{raxmlng_bootstrap_prefix}.raxml.log"
    params:
        msa = lambda wildcards: msas[wildcards.msa],
        model = lambda wildcards: raxmlng_models[wildcards.msa],
        prefix = raxmlng_bootstrap_prefix,
        threads = config["software"]["raxml-ng"]["threads"]
    log:
        raxmlng_log = f"{raxmlng_bootstrap_prefix}.log"
    shell:
        "{raxmlng_command} " +
        " --bootstrap " +
        "--msa {params.msa} "
        "--model {params.model} "
        "--bs-trees autoMRE "
        "--prefix {params.prefix} "
        "--threads {params.threads} "
        "--seed 0 "
        "> {log.raxmlng_log}"


rule compute_and_collect_tree_metrics:
    # input:
    #     database    = rules.move_db.output.database,
    #     bootstraps  = rules.raxmlng_bootstrap.output.bootstraps,
    #     eval_trees  = rules.collect_eval_trees.output.all_eval_trees,
    input:
        tree_metrics_input()
    output:
        raxmlng_tree_data = f"{db_path}raxmlng_tree_data.parquet",
    params:
        bootstrap_based_metrics = bootstrap_based_metrics,
        raxmlng_command = raxmlng_command
    script:
        "scripts/compute_and_collect_tree_metrics.py"
