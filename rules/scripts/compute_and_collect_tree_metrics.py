import pandas as pd
import sqlite3

from bootstrap_metrics import get_bootstrap_metrics
from tree_metrics import get_tree_characteristics


def bootstrap_metrics_all_trees(raxmlng, df):
    bootstrap_metrics = []

    for idx, row in df.iterrows():
        newick_tree = row.newick_eval
        bootstrap = get_bootstrap_metrics(raxmlng, newick_tree, snakemake.input.bootstraps)
        bootstrap = pd.DataFrame(bootstrap, index=[idx])
        bootstrap_metrics.append(bootstrap)

    return pd.concat(bootstrap_metrics)


def save_raxmlng_tree_data_and_add_tree_characteristics():
    con = sqlite3.connect(snakemake.input.database)
    df = pd.read_sql_query("SELECT * FROM RaxmlNGTree", con)
    con.close()

    # for each tree: get the tree characteristics and append the data to the dataframe
    tree_metrics = []

    for idx, row in df.iterrows():
        newick_tree = row.newick_eval
        metrics = get_tree_characteristics(newick_tree)
        metrics = pd.DataFrame(metrics, index=[idx])
        tree_metrics.append(metrics)

    if do_bootstrapping:
        bootstrap_metrics = bootstrap_metrics_all_trees(snakemake.params.raxmlng_command, df)
    else:
        bootstrap_metrics = pd.DataFrame()

    tree_metrics = pd.concat(tree_metrics)
    df = pd.concat([df, tree_metrics, bootstrap_metrics], axis=1)
    df.to_parquet(snakemake.output.raxmlng_tree_data)


if __name__ == "__main__":
    do_bootstrapping = snakemake.params.bootstrap_based_metrics
    save_raxmlng_tree_data_and_add_tree_characteristics()