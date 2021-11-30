from Bio import Phylo
import numpy as np


def get_tree_object(newick_str):
    trees = list(Phylo.NewickIO.Parser.from_string(newick_str).parse())
    return trees[0]


def get_all_branch_lengths_for_tree(newick_str):
    tree = get_tree_object(newick_str)
    return [node.branch_length for node in tree.find_clades(branch_length=True)]


def get_total_branch_length_for_tree(newick_str):
    all_brlens = get_all_branch_lengths_for_tree(newick_str)
    return sum(all_brlens)


def get_min_branch_length_for_tree(newick_str):
    all_brlens = get_all_branch_lengths_for_tree(newick_str)
    return min(all_brlens)


def get_max_branch_length_for_tree(newick_str):
    all_brlens = get_all_branch_lengths_for_tree(newick_str)
    return max(all_brlens)


def get_avg_branch_lenghts_for_tree(newick_str):
    all_brlens = get_all_branch_lengths_for_tree(newick_str)
    return np.mean(all_brlens)


def get_std_branch_lenghts_for_tree(newick_str):
    all_brlens = get_all_branch_lengths_for_tree(newick_str)
    return np.std(all_brlens)