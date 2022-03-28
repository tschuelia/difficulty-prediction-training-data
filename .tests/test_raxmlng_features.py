from tempfile import TemporaryDirectory

from fixtures import *


def test_get_patterns_gaps_invariant(raxmlng, example_msa_path):
    patterns, gaps, invariant = raxmlng.get_patterns_gaps_invariant(example_msa_path, "GTR+G")

    assert isinstance(patterns, int)
    assert isinstance(gaps, float)
    assert isinstance(invariant, float)

    assert patterns == 241
    assert gaps == pytest.approx(0.0789, abs=0.01)
    assert invariant == pytest.approx(0.5979, abs=0.01)


def test_get_rfdistance_results(raxmlng, newick_tree1, newick_tree2):
    with TemporaryDirectory() as tmpdir:
        trees = os.path.join(tmpdir, "trees.trees")

        with open(trees, "w") as f:
            f.write(newick_tree1 + "\n" + newick_tree2)

        num_topos, rel_rfdist, abs_rfdist = raxmlng.get_rfdistance_results(trees)

        assert isinstance(num_topos, float)
        assert num_topos == pytest.approx(1.0, abs=1e-6)

        assert isinstance(rel_rfdist, float)
        assert rel_rfdist == pytest.approx(0.0, abs=0.01)

        assert isinstance(abs_rfdist, float)
        assert abs_rfdist == pytest.approx(0.0, abs=0.01)


def test_infer_parsimony_trees(raxmlng, example_msa_path):
    n_trees = 2
    with TemporaryDirectory() as tmpdir:
        raxmlng.infer_parsimony_trees(
            msa_file=example_msa_path,
            model="GTR+G",
            prefix=os.path.join(tmpdir, "parsimony"),
            n_trees=n_trees
        )

        generated_files = os.listdir(tmpdir)

        assert "parsimony.raxml.startTree" in generated_files
        assert "parsimony.raxml.log" in generated_files

        generated_trees = [t for t in open(os.path.join(tmpdir, "parsimony.raxml.startTree")).readlines() if t]
        assert len(generated_trees) == n_trees
