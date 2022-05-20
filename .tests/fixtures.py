import pytest
import os
import yaml
import pickle

# unfortunately this is the way to go here since with snakemake I can't use relative imports (I think...)
import sys

sys.path.append(os.path.join(os.getcwd(), "rules/scripts"))

from msa_features import MSA
from raxmlng_features import RAxMLNG

@pytest.fixture
def example_msa_path():
    cwd = os.getcwd()
    return f"{cwd}/.tests/data/DNA/0.phy"


@pytest.fixture
def dna_phylip_msa(example_msa_path):
    return MSA(example_msa_path)


@pytest.fixture
def dna_fasta_msa():
    cwd = os.getcwd()
    fasta_file = f"{cwd}/.tests/data/DNA/0.fasta"
    return MSA(fasta_file)


@pytest.fixture
def small_msa():
    cwd = os.getcwd()
    phylip_file = f"{cwd}/.tests/data/DNA/small.fasta"
    return MSA(phylip_file)


@pytest.fixture
def raxmlng_inference_log():
    cwd = os.getcwd()
    return f"{cwd}/.tests/data/logs/raxml.inference.log"


@pytest.fixture
def raxmlng_inference_restarted_log():
    cwd = os.getcwd()
    return f"{cwd}/.tests/data/logs/raxml.inference.restart.log"


@pytest.fixture
def iqtree_siginficance_log():
    cwd = os.getcwd()
    return f"{cwd}/.tests/data/logs/significance.iqtree"


@pytest.fixture
def iqtree_siginficance_default_case_log():
    cwd = os.getcwd()
    return f"{cwd}/.tests/data/logs/significance.defaultCase.iqtree"



@pytest.fixture
def filtered_trees_cluster():
    cwd = os.getcwd()
    return pickle.load(open(f"{cwd}/.tests/data/logs/clusters.pkl", "rb"))


@pytest.fixture
def raxmlng_rfdistance_log():
    cwd = os.getcwd()
    return f"{cwd}/.tests/data/logs/raxml.rfdistance.log"


@pytest.fixture
def raxmlng_rfdistance_log_many_trees():
    cwd = os.getcwd()
    return f"{cwd}/.tests/data/logs/raxml.rfdistance.many.log"


@pytest.fixture
def raxmlng_rfdistance_log_identical_trees():
    cwd = os.getcwd()
    return f"{cwd}/.tests/data/logs/raxml.rfdistance.identical.log"


@pytest.fixture
def raxmlng_rfdistances():
    cwd = os.getcwd()
    return f"{cwd}/.tests/data/logs/raxml.rfdistances"


@pytest.fixture
def raxmlng_multiple_logs():
    cwd = os.getcwd()
    return f"{cwd}/.tests/data/logs/raxml.inference.multiple.log"


@pytest.fixture
def multiple_trees_path():
    cwd = os.getcwd()
    return f"{cwd}/.tests/data/trees/multiple.trees"


@pytest.fixture
def newick_tree1():
    cwd = os.getcwd()
    tf = f"{cwd}/.tests/data/trees/start.tree"
    return open(tf).readline().strip()


@pytest.fixture
def newick_tree2():
    cwd = os.getcwd()
    tf = f"{cwd}/.tests/data/trees/final.tree"
    return open(tf).readline().strip()


@pytest.fixture
def list_of_many_newick_trees():
    cwd = os.getcwd()
    tf = f"{cwd}/.tests/data/trees/many.trees"
    return [t.strip() for t in open(tf).readlines() if t]


@pytest.fixture
def list_of_trees_for_iqtree_statstests():
    cwd = os.getcwd()
    tf = f"{cwd}/.tests/data/trees/iqtree_statstest.trees"
    return [t.strip() for t in open(tf).readlines() if t]


@pytest.fixture
def raxmlng_command():
    cwd = os.getcwd()
    config = f"{cwd}/.tests/config.yaml"
    with open(config) as f:
        return yaml.safe_load(f)["software"]["raxml-ng"]["command"]


@pytest.fixture
def raxmlng(raxmlng_command):
    return RAxMLNG(raxmlng_command)