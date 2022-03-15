rule parsimony_tree:
    output:
        parsimony_tree  = parsimony_tree_file_name,
        log             = parsimony_log_file_name,
    params:
        msa     = lambda wildcards: msas[wildcards.msa],
        prefix  = output_files_parsimony_trees + "seed_{seed}",
        model   = lambda wildcards: raxmlng_model[wildcards.msa] if partitioned else raxmlng_model

    run:
        # Use RAxML-NG
        cmd = [
            "{raxmlng_command} ",
            "--start "
            "--msa {params.msa} ",
            "--tree pars{{1}} ",
            "--prefix {params.prefix} ",
            "--model {params.model} ",
            "> {output.log} "
        ]
        shell("".join(cmd))
