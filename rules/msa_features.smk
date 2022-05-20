rule compute_msa_features:
    output:
        msa_features =  f"{output_files_dir}msa_features.json"
    params:
        msa                 = lambda wildcards: msas[wildcards.msa],
        model               = lambda wildcards: raxmlng_models[wildcards.msa],
        raxmlng_command     = raxmlng_command
    script:
        "scripts/collect_msa_features.py"
