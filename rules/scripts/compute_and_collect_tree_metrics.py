import pandas as pd
import sqlite3

from bootstrap_metrics import get_bootstrap_metrics
from tree_metrics import get_tree_characteristics


def save_raxmlng_tree_data_and_add_tree_characteristics():
    con = sqlite3.connect(snakemake.input.database)
    df = pd.read_sql_query("SELECT * FROM RaxmlNGTree", con)
    con.close()

    raxmlng = snakemake.params.raxmlng_command
    bootstrap_replicates = snakemake.input.bootstraps

    # for each tree: get the tree characteristics and append the data to the dataframe
    tree_metrics = []
    bootstrap_metrics = []

    for idx, row in df.iterrows():
        newick_tree = row.newick_eval
        metrics = get_tree_characteristics(newick_tree)
        metrics = pd.DataFrame(metrics, index=[idx])
        tree_metrics.append(metrics)

        bootstrap = get_bootstrap_metrics(raxmlng, newick_tree, bootstrap_replicates)
        bootstrap = pd.DataFrame(bootstrap, index=[idx])
        bootstrap_metrics.append(bootstrap)

    tree_metrics = pd.concat(tree_metrics)
    bootstrap_metrics = pd.concat(bootstrap_metrics)
    df = pd.concat([df, tree_metrics, bootstrap_metrics], axis=1)
    df.to_parquet(snakemake.output.raxmlng_tree_data)


if __name__ == "__main__":
    save_raxmlng_tree_data_and_add_tree_characteristics()