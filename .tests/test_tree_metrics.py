from fixtures import *

from tree_metrics import *

from Bio import Phylo


def test_get_tree_object(newick_tree1):
    tree = get_tree_object(newick_tree1)

    assert isinstance(tree, Phylo.Newick.Tree)


def test_get_all_branch_lengths_for_tree(newick_tree1):
    branch_lengths = get_all_branch_lengths_for_tree(newick_tree1)
    branch_lengths.sort()

    expected = [
        1e-06,
        1e-06,
        0.001127,
        0.001152,
        0.001155,
        0.002272,
        0.011648,
        0.060709,
        0.071413,
    ]

    assert isinstance(branch_lengths, list)
    assert len(branch_lengths) == 9
    assert all([isinstance(bl, float) for bl in branch_lengths])
    assert all(
        [
            bl == pytest.approx(exp, abs=1e-6)
            for bl, exp in zip(branch_lengths, expected)
        ]
    )


def test_get_total_branch_length_for_tree(newick_tree1):
    total_brlen = get_total_branch_length_for_tree(newick_tree1)
    expected = 0.149478

    assert isinstance(total_brlen, float)
    assert total_brlen == pytest.approx(expected, abs=1e-6)


def test_get_min_branch_length_for_tree(newick_tree1):
    min_brlen = get_min_branch_length_for_tree(newick_tree1)
    expected = 1e-6

    assert isinstance(min_brlen, float)
    assert min_brlen == pytest.approx(expected, abs=1e-6)


def test_get_max_branch_length_for_tree(newick_tree1):
    max_brlen = get_max_branch_length_for_tree(newick_tree1)
    expected = 0.071413

    assert isinstance(max_brlen, float)
    assert max_brlen == pytest.approx(expected, abs=1e-6)


def test_get_avg_branch_length_for_tree(newick_tree1):
    avg_brlen = get_avg_branch_lengths_for_tree(newick_tree1)
    expected = 0.01660866

    assert isinstance(avg_brlen, float)
    assert avg_brlen == pytest.approx(expected, abs=1e-6)


def test_get_std_branch_length_for_tree(newick_tree1):
    std_brlen = get_std_branch_lengths_for_tree(newick_tree1)
    expected = 0.0267655

    assert isinstance(std_brlen, float)
    assert std_brlen == pytest.approx(expected, abs=1e-6)
