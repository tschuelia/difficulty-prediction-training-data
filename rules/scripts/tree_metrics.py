import statistics
from ete3 import Tree
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr

R_SPECTRAL = importr("RPANDA")
R_APE = importr("ape")


def _get_stats_for_batch(brlen_batch, suffix):
    return {
        f"mean_{suffix}": statistics.mean(brlen_batch),
        f"median_{suffix}": statistics.median(brlen_batch),
        f"stdev_{suffix}": statistics.stdev(brlen_batch),
        f"total_{suffix}": sum(brlen_batch),
        f"min_{suffix}": min(brlen_batch),
        f"max_{suffix}": max(brlen_batch),
    }


def get_internal_external_brlens(newick_tree):
    tree = Tree(newick_tree)
    external_brlens = []
    internal_brlens = []

    internal_nodes = []

    for node in tree.traverse():
        if node.is_leaf():
            external_brlens.append(node.dist)
        else:
            internal_brlens.append(node.dist)
            internal_nodes.append(node)

    return internal_brlens, external_brlens


def compute_brlen_statistics(internal_brlens, external_brlens):
    all_brlens = external_brlens + internal_brlens

    data = {
        **_get_stats_for_batch(external_brlens, "external_brlens"),
        **_get_stats_for_batch(internal_brlens, "internal_brlens"),
        **_get_stats_for_batch(all_brlens, "all_brlens")
    }

    return data


def get_treeness(brlen_stats):
    sum_internal = brlen_stats["total_internal_brlens"]
    sum_total = brlen_stats["total_all_brlens"]

    return {"treeness": sum_internal / sum_total}


def get_percentage_near_zero_brlens(brlens):
    near_zero_brlens = [b for b in brlens if b < 0.0001]
    return len(near_zero_brlens) / len(brlens)


def spectral_analysis(newick_tree):
    phylo = f"read.tree(text='{newick_tree}')"
    try:
        spectral = R_SPECTRAL.spectR(robjects.r(phylo))
        spectral = {key: spectral.rx2(key)[0] for key in spectral.names}
        return spectral
    except:
        return None


def tree_diameter(newick_tree):
    tree = Tree(newick_tree)

    diameter = 0.0

    for node in tree.traverse():
        _, max_dist = node.get_farthest_node()
        diameter = max(diameter, max_dist)

    return diameter


def get_tree_characteristics(newick_tree):
    internal_brlens, external_brlens = get_internal_external_brlens(newick_tree)
    brlen_stats = compute_brlen_statistics(internal_brlens, external_brlens)

    data = {
        "internal_brlens": [internal_brlens],
        "external_brlens": [external_brlens],
        **brlen_stats,
        "near_zero_percentage_internal": get_percentage_near_zero_brlens(internal_brlens),
        "near_zero_percentage_external": get_percentage_near_zero_brlens(external_brlens),
        **spectral_analysis(newick_tree),
        "diameter": tree_diameter(newick_tree)
    }

    return data