import pytest

from fixtures import *

from raxmlng_parser import *


def test_get_raxmlng_llh(raxmlng_inference_log):
    llh = get_raxmlng_llh(raxmlng_inference_log)

    assert isinstance(llh, float)
    assert llh == pytest.approx(-1662.27, abs=0.1)


def test_get_raxmlng_llh_raises_value_error(raxmlng_rfdistance_log):
    with pytest.raises(ValueError):
        get_raxmlng_llh(raxmlng_rfdistance_log)


def test_get_raxmlng_starting_llh(raxmlng_inference_log):
    llh = get_raxmlng_starting_llh(raxmlng_inference_log)

    assert isinstance(llh, float)
    assert llh == pytest.approx(-2260.177, abs=0.1)


def test_get_raxmlng_starting_llh_raises_value_error(raxmlng_rfdistance_log):
    with pytest.warns(Warning):
        get_raxmlng_starting_llh(raxmlng_rfdistance_log)


def test_get_all_raxmlng_llhs(raxmlng_multiple_logs):
    llhs = get_all_raxmlng_llhs(raxmlng_multiple_logs)

    assert isinstance(llhs, list)
    assert len(llhs) == 2
    assert all([isinstance(llh, float) for llh in llhs])
    assert llhs[0] == pytest.approx(-1666.269, abs=0.1)
    assert llhs[1] == pytest.approx(-1662.269, abs=0.1)


def test_get_best_raxmlng_llh(raxmlng_multiple_logs):
    llh = get_best_raxmlng_llh(raxmlng_multiple_logs)

    assert isinstance(llh, float)
    assert llh == pytest.approx(-1662.269, abs=0.1)


def test_get_raxmlng_elapsed_time(raxmlng_inference_log):
    time = get_raxmlng_elapsed_time(raxmlng_inference_log)

    assert isinstance(time, float)
    assert time == pytest.approx(0.025, abs=0.01)


def test_get_raxmlng_elapsed_time_with_restarts(raxmlng_inference_restarted_log):
    time = get_raxmlng_elapsed_time(raxmlng_inference_restarted_log)

    assert isinstance(time, float)
    assert time == pytest.approx(111.111, abs=0.01)


def test_get_raxmlng_elapsed_time_raises_value_error(raxmlng_rfdistances):
    with pytest.raises(ValueError):
        get_raxmlng_elapsed_time(raxmlng_rfdistances)


def test_get_raxmlng_runtimes(raxmlng_multiple_logs):
    times = get_raxmlng_runtimes(raxmlng_multiple_logs)

    assert isinstance(times, list)
    assert len(times) == 2
    assert all([isinstance(t, float) for t in times])
    assert times[0] == pytest.approx(0.025, abs=0.01)
    assert times[1] == pytest.approx(0.027, abs=0.01)


def test_get_raxmlng_runtimes_raises_value_error(raxmlng_rfdistances):
    with pytest.raises(ValueError):
        get_raxmlng_runtimes(raxmlng_rfdistances)


def test_get_raxmlng_num_spr_rounds(raxmlng_inference_log):
    slow, fast = get_raxmlng_num_spr_rounds(raxmlng_inference_log)

    assert isinstance(slow, int)
    assert isinstance(fast, int)

    assert slow == 1
    assert fast == 1


def test_get_raxmlng_num_spr_rounds_no_spr_rounds_logged(raxmlng_rfdistance_log):
    slow, fast = get_raxmlng_num_spr_rounds(raxmlng_rfdistance_log)

    assert isinstance(slow, int)
    assert isinstance(fast, int)

    assert slow == 0
    assert fast == 0


def test_rel_rfdistance_starting_final(newick_tree1, newick_tree2, raxmlng_command):
    rfdist = rel_rfdistance_starting_final(newick_tree1, newick_tree2, raxmlng_command)

    assert isinstance(rfdist, float)
    assert rfdist == pytest.approx(0.0, abs=0.01)


def test_get_all_parsimony_scores(raxmlng_multiple_logs):
    scores = get_all_parsimony_scores(raxmlng_multiple_logs)

    assert isinstance(scores, list)
    assert len(scores) == 2
    assert all([isinstance(score, int) for score in scores])
    assert scores[0] == 119
    assert scores[1] == 123


def test_get_patterns_gaps_invariant(raxmlng, raxmlng_inference_log):
    patterns, gaps, invariant = get_patterns_gaps_invariant(raxmlng_inference_log)

    assert isinstance(patterns, int)
    assert isinstance(gaps, float)
    assert isinstance(invariant, float)

    assert patterns == 50
    assert gaps == pytest.approx(0.0163, abs=0.01)
    assert invariant == pytest.approx(0.8970, abs=0.01)


def test_get_patterns_gaps_invariant_raises_value_error(raxmlng_rfdistance_log):
    with pytest.raises(ValueError):
        get_patterns_gaps_invariant(raxmlng_rfdistance_log)
