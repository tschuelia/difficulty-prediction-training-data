import json

from msa_features import *

msa_file = snakemake.params.msa
raxmlng_command = snakemake.params.raxmlng_command
data_type = snakemake.params.data_type
model = snakemake.params.model
msa = read_alignment(msa_file, data_type=data_type)

msa_features = {
    "taxa": get_number_of_taxa(msa),
    "sites": get_number_of_sites(msa),
    "patterns": get_number_of_patterns(msa_file, raxmlng_executable=raxmlng_command, model=model),
    "gaps": get_percentage_of_gaps(msa_file, raxmlng_executable=raxmlng_command, model=model),
    "invariant": get_percentage_of_invariant_sites(msa_file, raxmlng_executable=raxmlng_command, model=model),
    "entropy": get_msa_avg_entropy(msa),
    "column_entropies": get_msa_column_entropies(msa),
    "bollback": bollback_multinomial(msa),
    "treelikeness": treelikeness_score(msa, data_type),
    "char_frequencies": get_character_frequencies(msa)
}

with open(snakemake.output.msa_features, "w") as f:
    json.dump(msa_features, f)
