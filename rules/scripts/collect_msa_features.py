import json

from msa_features import MSA
from raxmlng_features import RAxMLNG

msa_file = snakemake.params.msa
model = snakemake.params.model

msa = MSA(msa_file)
raxmlng = RAxMLNG(snakemake.params.raxmlng_command)

patterns, gaps, invariant = raxmlng.get_patterns_gaps_invariant(msa_file, model)

msa_features = {
    "taxa": msa.get_number_of_taxa(),
    "sites": msa.get_number_of_sites(),
    "patterns": patterns,
    "gaps": gaps,
    "invariant": invariant,
    "entropy": msa.get_avg_entropy(),
    "column_entropies": msa.get_column_entropies(),
    "bollback": msa.bollback_multinomial(),
    "treelikeness": msa.treelikeness_score(),
    "char_frequencies": msa.get_character_frequencies()
}

with open(snakemake.output.msa_features, "w") as f:
    json.dump(msa_features, f)
