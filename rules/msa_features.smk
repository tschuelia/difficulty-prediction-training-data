model = config["software"]["raxml-ng"]["model"]
data_type = config["data_type"]

rule compute_msa_features:
    output:
        msa_features =  f"{output_files_dir}msa_features.json"
    params:
        msa                 = lambda wildcards: msas[wildcards.msa],
        raxmlng_command     = raxmlng_command,
        model               = model,
        data_type           = data_type
    script:
        "scripts/collect_msa_features.py"
