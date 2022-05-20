from fixtures import *

from filter_tree_topologies import *


def test_get_rfdist_clusters(raxmlng_rfdistance_log_many_trees, list_of_many_newick_trees):
    clusters = get_rfdist_clusters(raxmlng_rfdistance_log_many_trees, list_of_many_newick_trees)

    assert isinstance(clusters, list)
    assert len(clusters) == 6
    assert all([isinstance(c, set) for c in clusters])

    expected_cluster_lengths = [7, 52, 13, 2, 8, 4]

    assert all([len(c) == exp for c, exp in zip(clusters, expected_cluster_lengths)])


def test_filter_tree_topologies_many_trees(raxmlng_rfdistance_log_many_trees, list_of_many_newick_trees):
    unique_trees, clusters = filter_tree_topologies(list_of_many_newick_trees, raxmlng_rfdistance_log_many_trees)

    assert isinstance(unique_trees, list)
    assert len(unique_trees) == 6
    assert all([isinstance(t, str) for t in unique_trees])

    assert isinstance(clusters, list)
    assert len(clusters) == 6
    assert all([isinstance(c, set) for c in clusters])

    expected_cluster_lengths = [7, 52, 13, 2, 8, 4]

    assert all([len(c) == exp for c, exp in zip(clusters, expected_cluster_lengths)])


def test_filter_tree_topologies_single_tree(newick_tree1):
    unique_trees, clusters = filter_tree_topologies([newick_tree1], "")

    assert isinstance(unique_trees, list)
    assert len(unique_trees) == 1
    assert all([isinstance(t, str) for t in unique_trees])

    assert isinstance(clusters, list)
    assert len(clusters) == 1
    assert all([isinstance(c, set) for c in clusters])

    assert len(clusters[0]) == 1


def test_filter_tree_topologies_two_trees(newick_tree1, raxmlng_rfdistance_log_identical_trees):
    unique_trees, clusters = filter_tree_topologies([newick_tree1, newick_tree1], raxmlng_rfdistance_log_identical_trees)

    assert isinstance(unique_trees, list)
    assert len(unique_trees) == 1
    assert all([isinstance(t, str) for t in unique_trees])

    assert isinstance(clusters, list)
    assert len(clusters) == 1
    assert all([isinstance(c, set) for c in clusters])

    assert len(clusters[0]) == 1