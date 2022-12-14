rule parsimony_tree:
    output:
        parsimony_tree  = output_files_parsimony_trees / "parsimony.raxml.startTree",
        log             = output_files_parsimony_trees / "parsimony.raxml.log"
    params:
        msa     = lambda wildcards: msas[wildcards.msa],
        prefix  = str(output_files_parsimony_trees / "parsimony"),
        model   = lambda wildcards: raxmlng_models[wildcards.msa]
    run:
        # Use RAxML-NG
        cmd = [
            "{raxmlng_command} ",
            "--start "
            "--msa {params.msa} ",
            "--tree pars{{100}} ",
            "--prefix {params.prefix} ",
            "--model {params.model} ",
            "--seed 0 "
            "> {output.log} "
        ]
        shell("".join(cmd))
