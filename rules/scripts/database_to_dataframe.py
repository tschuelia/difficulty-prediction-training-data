import sqlite3
import pandas as pd
import numpy as np

from custom_types import *


def get_difficulty_labels(df: pd.DataFrame) -> List[float]:
    """
    difficult if:
    - avg_rfdist_plausible is close to 1.0 -> + val
    - num_topos_plausible/num_trees_plausible is close to 1.0 -> + val
    - proportion_plausible is closer to 0.0 -> + (1-val)
    """
    labels = []

    for idx, row in df.iterrows():
        diff_proba = 0
        ct = 0

        for col in [
            "avg_rfdist_eval",
            "avg_rfdist_plausible",
            "num_topos_eval/num_trees_eval",
            "num_topos_plausible/num_trees_plausible",
        ]:
            if (row[col] > -np.inf) and (row[col] < np.inf) and (row[col] is not None):
                diff_proba += row[col]
                ct += 1

        for col in [
            "proportion_plausible"
        ]:
            if (
                    (row[col] > -np.inf)
                    and (row[col] < np.inf)
                    and (row[col] is not None)
            ):
                diff_proba += 1 - row[col]
                ct += 1
        labels.append(diff_proba / ct)

    return labels


if __name__ == "__main__":
    db_path = snakemake.input.database
    parquet_path = snakemake.output.dataframe
    num_pars_trees = snakemake.params.num_pars_trees
    num_rand_trees = snakemake.params.num_rand_trees
    num_trees = num_rand_trees + num_pars_trees
    num_parsimony_trees = snakemake.params.num_parsimony_trees

    con = sqlite3.connect(db_path)

    df = pd.read_sql_query("SELECT * FROM dataset", con)

    # fmt: off
    df["num_topos_plausible/num_trees_plausible"]   = df["num_topos_plausible"] / df["num_trees_plausible"]
    df["num_topos_parsimony/num_trees_parsimony"]   = df["num_topos_parsimony"] / num_parsimony_trees
    df["num_topos_search/num_trees_search"]         = df["num_topos_search"] / num_trees
    df["num_topos_eval/num_trees_eval"]             = df["num_topos_eval"] / num_trees
    df["num_patterns/num_taxa"]                     = df["num_patterns"] / df["num_taxa"]
    df["num_sites/num_taxa"]                        = df["num_sites"] / df["num_taxa"]
    df["difficult"] = get_difficulty_labels(df)
    # fmt: on

    df.to_parquet(parquet_path)

    con.close()