from fixtures import *

from save_best_eval_tree import *


def test_get_best_tree(multiple_trees_path, raxmlng_multiple_logs):
    best_llh, best_tree = get_best_tree_and_llh(multiple_trees_path, raxmlng_multiple_logs)

    assert isinstance(best_tree, str)
    assert isinstance(best_llh, float)

    assert best_llh == pytest.approx(-1662.268, abs=0.1)

    tree_strings = open(multiple_trees_path).readlines()
    assert best_tree == tree_strings[1].strip()
