import os
import pytest

# unfortunately this is the way to go here since with snakemake I can't use relative imports (I think...)
import sys

sys.path.append(os.path.join(os.getcwd(), "rules/scripts"))

from msa_features import MSA


@pytest.fixture
def dna_phylip_msa():
    cwd = os.getcwd()
    phylip_file = f"{cwd}/tests/data/DNA/0.phy"
    return MSA(phylip_file)


@pytest.fixture
def dna_fasta_msa():
    cwd = os.getcwd()
    fasta_file = f"{cwd}/tests/data/DNA/0.fasta"
    return MSA(fasta_file)


@pytest.fixture
def small_msa():
    cwd = os.getcwd()
    phylip_file = f"{cwd}/tests/data/DNA/small.fasta"
    return MSA(phylip_file)


class TestMSAFeatures:
    def test_get_msa_file_format(self, dna_phylip_msa, dna_fasta_msa):
        ff = dna_phylip_msa._get_file_format()
        assert ff == "phylip-relaxed"

        ff = dna_fasta_msa._get_file_format()
        assert ff == "fasta"

    def test_guess_msa_file_data_type(self):
        cwd = os.getcwd()
        for true_type in ["DNA", "AA"]:
            base_dir = f"{cwd}/tests/data/{true_type}/"
            for msa_file in os.listdir(base_dir):
                msa_file = base_dir + msa_file
                msa = MSA(msa_file)
                guessed_type = msa.guess_data_type()
                assert guessed_type == true_type

    def test_number_of_taxa(self, dna_phylip_msa):
        assert dna_phylip_msa.get_number_of_taxa() == 68

    def test_number_of_sites(self, dna_phylip_msa):
        assert dna_phylip_msa.get_number_of_sites() == 766

    def test_get_avg_entropy(self, small_msa):
        # TODO
        pass

    def test_bollback_multinomial(self, small_msa):
        # TODO
        pass

    def test_treelikeness_score(self, small_msa):
        # TODO
        pass

    def test_get_character_frequencies(self, small_msa):
        char_frequencies = small_msa.get_character_frequencies()
        assert char_frequencies["A"] == 280
        assert char_frequencies["C"] == 160
        assert char_frequencies["T"] == 260
        assert char_frequencies["G"] == 210
        assert char_frequencies["N"] == 3640
        assert char_frequencies["-"] == 670

        assert (
            sum([v for v in char_frequencies.values()])
            == small_msa.get_number_of_taxa() * small_msa.get_number_of_sites()
        )
