import subprocess

from custom_types import *


def read_file_contents(file_path: FilePath) -> list[str]:
    with open(file_path) as f:
        content = f.readlines()

    return [l.strip() for l in content]


def get_value_from_line(line: str, search_string: str) -> float:
    line = line.strip()
    if search_string in line:
        _, value = line.rsplit(" ", 1)
        return float(value)

    raise ValueError(
        f"The given line '{line}' does not contain the search string '{search_string}'."
    )


def get_single_value_from_file(input_file: FilePath, search_string: str) -> float:
    with open(input_file) as f:
        lines = f.readlines()

    for l in lines:
        if search_string in l:
            return get_value_from_line(l, search_string)

    raise ValueError(
        f"The given input file {input_file} does not contain the search string '{search_string}'."
    )


def get_multiple_values_from_file(
    input_file: FilePath, search_string: str
) -> list[float]:
    with open(input_file) as f:
        lines = f.readlines()

    values = []
    for l in lines:
        if search_string in l:
            values.append(get_value_from_line(l, search_string))

    if len(values) == 0:
        raise ValueError(
            f"The given input file {input_file} does not contain the search string '{search_string}'."
        )

    return values


def run_cmd(cmd: Command) -> None:
    try:
        subprocess.check_output(cmd)
    except Exception as e:
        print(f"Error running command \"{' '.join(cmd)}\"")
        raise e
