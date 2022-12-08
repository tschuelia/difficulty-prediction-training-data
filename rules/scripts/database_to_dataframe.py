import sqlite3
import pandas as pd
import numpy as np

from custom_types import *
from tree_metrics import get_tree_characteristics


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


def save_training_data():
    num_pars_trees = snakemake.params.num_pars_trees
    num_rand_trees = snakemake.params.num_rand_trees
    num_trees = num_rand_trees + num_pars_trees

    num_parsimony_trees = snakemake.params.num_parsimony_trees

    con = sqlite3.connect(snakemake.input.database)

    df = pd.read_sql_query("SELECT * FROM dataset", con)

    # fmt: off
    df["num_topos_plausible/num_trees_plausible"] = df["num_topos_plausible"] / df["num_trees_plausible"]
    df["num_topos_parsimony/num_trees_parsimony"] = df["num_topos_parsimony"] / num_parsimony_trees
    df["num_topos_search/num_trees_search"] = df["num_topos_search"] / num_trees
    df["num_topos_eval/num_trees_eval"] = df["num_topos_eval"] / num_trees
    df["num_patterns/num_taxa"] = df["num_patterns"] / df["num_taxa"]
    df["num_sites/num_taxa"] = df["num_sites"] / df["num_taxa"]
    df["difficult"] = get_difficulty_labels(df)
    # fmt: on

    df.to_parquet(snakemake.output.training_data)

    con.close()


def save_raxmlng_tree_data_and_add_tree_characteristics():
    con = sqlite3.connect(snakemake.input.database)
    df = pd.read_sql_query("SELECT * FROM RaxmlNGTree", con)
    con.close()

    # for each tree: get the tree characteristics and append the data to the dataframe
    tree_metrics = []
    for idx, row in df.iterrows():
        newick_tree = row.newick_eval
        metrics = get_tree_characteristics(newick_tree)
        metrics = pd.DataFrame(metrics, index=row.index)
        tree_metrics.append(metrics)

    df.to_parquet(snakemake.output.raxmlng_tree_data)



if __name__ == "__main__":
    save_training_data()
    save_raxmlng_tree_data_and_add_tree_characteristics()
