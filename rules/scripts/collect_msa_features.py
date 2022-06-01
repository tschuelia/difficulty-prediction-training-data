import json

from pyphypred.msa import MSA
from pyphypred.raxmlng import RAxMLNG

msa_file = snakemake.params.msa
model = snakemake.params.model

msa = MSA(msa_file)
raxmlng = RAxMLNG(snakemake.params.raxmlng_command)

patterns, gaps, invariant = raxmlng.get_patterns_gaps_invariant(msa_file, model)

msa_features = {
    "taxa": msa.number_of_taxa(),
    "sites": msa.number_of_sites(),
    "patterns": patterns,
    "gaps": gaps,
    "invariant": invariant,
    "entropy": msa.entropy(),
    "column_entropies": msa.column_entropies(),
    "bollback": msa.bollback_multinomial(),
    "treelikeness": msa.treelikeness_score(),
}

with open(snakemake.output.msa_features, "w") as f:
    json.dump(msa_features, f)
