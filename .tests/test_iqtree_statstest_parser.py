from fixtures import *
from iqtree_statstest_parser import *


def test_get_relevant_section(iqtree_siginficance_log):
    section = get_relevant_section(iqtree_siginficance_log)

    print(section[0])

    assert isinstance(section, list)
    assert len(section) > 1
    assert all([isinstance(line, str) for line in section])
    assert section[0].startswith(START_STRING)


def test_get_relevant_section_raises_value_error_for_incorrect_files(
        raxmlng_inference_log,
):
    with pytest.raises(ValueError):
        get_relevant_section(raxmlng_inference_log)


def test_get_cleaned_table_entries(iqtree_siginficance_log):
    section = get_relevant_section(iqtree_siginficance_log)
    table_entries = get_cleaned_table_entries(section)

    assert isinstance(table_entries, list)
    assert len(table_entries) == 7
    assert all([isinstance(tree_id, int) for tree_id, _, _, _ in table_entries])
    assert all([isinstance(llh, float) for _, llh, _, _ in table_entries])
    assert all([isinstance(deltaL, float) for _, _, deltaL, _ in table_entries])
    assert all([isinstance(results, list) for _, _, _, results in table_entries])

    assert all([len(results) == 7 for _, _, _, results in table_entries])
    assert all([isinstance(r, str) for _, _, _, results in table_entries for r in results])


def test_get_cleaned_table_entries_raises_value_error_for_incorrect_section(raxmlng_inference_log):
    section = open(raxmlng_inference_log).readlines()
    with pytest.raises(ValueError):
        get_cleaned_table_entries(section)


def test_get_names_of_performed_tests(iqtree_siginficance_log):
    section = get_relevant_section(iqtree_siginficance_log)
    test_names = get_names_of_performed_tests(section)

    assert isinstance(test_names, list)
    assert len(test_names) == 7
    assert test_names == ["bp-RELL", "p-KH", "p-SH", "p-WKH", "p-WSH", "c-ELW", "p-AU"]


def test_get_names_of_performed_tests_raises_value_error_for_incorrect_section(raxmlng_inference_log):
    section = open(raxmlng_inference_log).readlines()

    with pytest.raises(ValueError):
        get_names_of_performed_tests(section)


def test_get_iqtree_results(iqtree_siginficance_log):
    results = get_iqtree_results(iqtree_siginficance_log)

    assert isinstance(results, list)
    assert all([isinstance(res, dict) for res in results])
    assert all([len(res.keys()) == 4 for res in results])
    assert all([len(res["tests"].keys()) == 7 for res in results])

    section = get_relevant_section(iqtree_siginficance_log)
    test_names = get_names_of_performed_tests(section)

    assert all([list(res["tests"].keys()) == test_names for res in results])

    res0 = results[0]
    assert res0["logL"] == pytest.approx(-64759.701, abs=0.1)
    assert res0["deltaL"] == pytest.approx(25.028)
    assert res0["tests"]["bp-RELL"]["score"] == pytest.approx(0.0503)
    assert res0["tests"]["p-KH"]["score"] == pytest.approx(0.185)
    assert res0["tests"]["p-SH"]["score"] == pytest.approx(0.268)
    assert res0["tests"]["p-WKH"]["score"] == pytest.approx(0.185)
    assert res0["tests"]["p-WSH"]["score"] == pytest.approx(0.516)
    assert res0["tests"]["c-ELW"]["score"] == pytest.approx(0.0502)
    assert res0["tests"]["p-AU"]["score"] == pytest.approx(0.262)
    assert res0["plausible"]

    res4 = results[4]
    assert not res4["plausible"]


def test_get_iqtree_results_with_single_result(iqtree_siginficance_default_case_log):
    results = get_iqtree_results(iqtree_siginficance_default_case_log)

    assert len(results) == 1
    assert results[0] == {
                "deltaL": 0,
                "plausible": 1,
                "tests": {
                    'bp-RELL': {
                        'score': 1,
                        'significant': True
                    },
                    'p-KH': {
                        'score': 1,
                        'significant': True
                    },
                    'p-SH': {
                        'score': 1,
                        'significant': True
                    },
                    'p-WKH': {
                        'score': 1,
                        'significant': True
                    },
                    'p-WSH': {
                        'score': 1,
                        'significant': True
                    },
                    'c-ELW': {
                        'score': 1,
                        'significant': True
                    },
                    'p-AU': {
                        'score': 1,
                        'significant': True
                    }
                }
            }


def test_get_iqtree_results_for_tree_str(iqtree_siginficance_log, filtered_trees_cluster, list_of_trees_for_iqtree_statstests):
    results = get_iqtree_results(iqtree_siginficance_log)

    expected_cluster_ids = [0, 1, 2, 1, 1]

    for i, tree in enumerate(list_of_trees_for_iqtree_statstests):
        res, cluster_id = get_iqtree_results_for_eval_tree_str(results, tree, filtered_trees_cluster)

        assert cluster_id == expected_cluster_ids[i]
        assert res["plausible"]


def test_get_iqtree_results_for_tree_str_raises_value_error_for_unknown_tree(iqtree_siginficance_log, filtered_trees_cluster, newick_tree1):
    results = get_iqtree_results(iqtree_siginficance_log)

    with pytest.raises(ValueError):
        get_iqtree_results_for_eval_tree_str(results, newick_tree1, filtered_trees_cluster)