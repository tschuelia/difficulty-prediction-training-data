import pytest

from fixtures import *
from utils import *


def test_get_value_from_line():
    line = "Test number: 100.0"
    value = get_value_from_line(line, "number")

    assert isinstance(value, float)
    assert value == pytest.approx(100)


def test_get_value_from_line_raises_value_error_if_string_not_in_line():
    line = "Test number: 100.0"

    with pytest.raises(ValueError):
        get_value_from_line(line, "bananans")


def test_get_single_value_from_file(raxmlng_inference_log):
    value = get_single_value_from_file(raxmlng_inference_log, "random seed:")

    assert isinstance(value, float)
    assert value == pytest.approx(0)


def test_get_single_value_from_file_raises_value_error_if_string_not_in_file(raxmlng_inference_log):
    with pytest.raises(ValueError):
        get_single_value_from_file(raxmlng_inference_log, "bananans")


def test_get_multiple_values_from_file(raxmlng_multiple_logs):
    values = get_multiple_values_from_file(raxmlng_multiple_logs, "random seed:")

    assert isinstance(values, list)
    assert all([isinstance(value, float) for value in values])
    assert values[0] == pytest.approx(0)
    assert values[1] == pytest.approx(1)


def test_get_multiple_values_from_file_raises_value_error_if_string_not_in_file(raxmlng_multiple_logs):
    with pytest.raises(ValueError):
        get_multiple_values_from_file(raxmlng_multiple_logs, "bananans")

