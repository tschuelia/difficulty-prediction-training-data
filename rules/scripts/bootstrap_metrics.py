import os
import statistics
from ete3 import Tree

from snakemake import shell
from tempfile import TemporaryDirectory


def _run_raxmlng_support(raxmlng_command, tree_file, bootstraps, prefix):
    shell(f"{raxmlng_command} --support --tree {tree_file} --bs-trees {bootstraps} --prefix {prefix} --threads 2")


def _extract_support_values(newick):
    tree = Tree(newick)
    return [n.support for n in tree.traverse() if not n.is_leaf() and not n.is_root()]


def get_bootstrap_support_values_for_tree(raxmlng_command, newick_tree, bootstraps):
    with TemporaryDirectory() as tmpdir:
        prefix = os.path.join(tmpdir, "bootstrap")
        tree_file = os.path.join(tmpdir, "newick.tree")

        with open(tree_file, "w") as f:
            f.write(newick_tree)

        _run_raxmlng_support(raxmlng_command, tree_file, bootstraps, prefix)

        tree_with_support_values = f"{prefix}.raxml.support"
        return _extract_support_values(open(tree_with_support_values).readline())


def get_bootstrap_metrics(raxmlng_command, newick_tree, bootstraps):
    support_values = get_bootstrap_support_values_for_tree(raxmlng_command, newick_tree, bootstraps)
    high_support_values = [s for s in support_values if s > 70.0]
    proportion_high = len(high_support_values) / len(support_values)

    data = {
        "average": statistics.mean(support_values),
        "standard_dev": statistics.stdev(support_values),
        "median": statistics.median(support_values),
        "minimum": min(support_values),
        "maximum": max(support_values),
        "total": sum(support_values),
        "proportion_greater_70": proportion_high,
        "raw_values": [support_values],
    }

    return data

