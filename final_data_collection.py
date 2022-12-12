import pandas as pd
import pathlib
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--dir",
        type=str,
        required=True,
        help="Results directory to look for parquet files.",
    )

    args = parser.parse_args()

    directory = pathlib.Path(args.dir)

    dfs = []

    for msa_dir in directory.iterdir():
        if not msa_dir.is_dir():
            continue
        dataset_parquet = msa_dir / "training_data.parquet"
        raxmlng_parquet = msa_dir / "raxmlng_tree_data.parquet"

        if not dataset_parquet.exists() or not raxmlng_parquet.exists():
            print(f"Skipping {msa_dir}: not all parquet files present.")

        df = pd.read_parquet(dataset_parquet)

        rax_df = pd.read_parquet(raxmlng_parquet)
        rax_df.is_best = rax_df.is_best.astype(bool)

        rax_df = rax_df.loc[rax_df.is_best]

        assert rax_df.shape[0] >= 1, f"The data for MSA {msa_dir} does not seem to contain a best_tree, please check the RAxML-NG logs..."

        rax_df = rax_df.drop(["id", "uuid", "dataset_id", "dataset_uuid"], axis=1)
        rax_df = rax_df.head(1)

        df = pd.concat([df, rax_df], axis=1)
        dfs.append(df)

    final_df = pd.concat(dfs, ignore_index=True)
    final_df.to_parquet(directory / "all_data.parquet")






