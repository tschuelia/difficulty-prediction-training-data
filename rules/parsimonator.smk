rule parsimonator_parsimony_tree:
    output:
        parsimony_tree = parsimonator_tree_file_name,
        log = parsimonator_log_file_name,
    params:
        msa = lambda wildcards: msas[wildcards.msa]
    shell:
        # first run parsimonator to generate one parsimony tree
        "{parsimonator_command} " 
        "-s {params.msa} "
        "-N 1 "
        "-p {wildcards.seed} "
        "-n seed_{wildcards.seed} "
        "> {output.log} "
        # copy the output files to the desired place
        " && mv RAxML_parsimonyTree.seed_{wildcards.seed}.0 {output.parsimony_tree} "
        # finally remove the output files from the current workdir
        " && rm RAxML_info.seed_{wildcards.seed} "

