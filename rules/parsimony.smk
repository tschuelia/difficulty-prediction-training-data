# Depending on the data type we use either Parsimonator or RAxML-NG
data_type = config["data_type"]
model = config["software"]["raxml-ng"]["model"]
threads = config["software"]["raxml-ng"]["threads"]

rule parsimony_tree:
    output:
        parsimony_tree=parsimony_tree_file_name,
        log=parsimony_log_file_name,
    params:
        msa     =lambda wildcards: msas[wildcards.msa],
        prefix  = output_files_parsimony_trees + "seed_{seed}"
    run:
        if data_type == "DNA":
            # Use Parsimonator
            cmd = [
                # first run parsimonator to generate one parsimony tree
                "{parsimonator_command} ",
                "-s {params.msa} ",
                "-N 1 ",
                "-p {wildcards.seed} ",
                "-n seed_{wildcards.seed} ",
                "> {output.log} ",
                # copy the output files to the desired place
                " && mv RAxML_parsimonyTree.seed_{wildcards.seed}.0 {output.parsimony_tree} ",
                # finally remove the output files from the current workdir
                " && rm RAxML_info.seed_{wildcards.seed} ",
            ]
            shell("".join(cmd))
        else:
            # Use RAxML-NG
            cmd = [
                "{raxmlng_command} ",
                "--start "
                "--msa {params.msa} ",
                "--tree pars{{1}} ",
                "--prefix {params.prefix} ",
                "--model {model} ",
                "--threads {threads} "
                "> {output.log} "
            ]
            shell("".join(cmd))
