from fixtures import *


class TestMSAFeatures:
    def test_get_msa_file_format(self, dna_phylip_msa, dna_fasta_msa):
        ff = dna_phylip_msa._get_file_format()
        assert ff == "phylip-relaxed"

        ff = dna_fasta_msa._get_file_format()
        assert ff == "fasta"

    def test_guess_msa_file_data_type(self):
        cwd = os.getcwd()
        for true_type in ["DNA", "AA"]:
            base_dir = f"{cwd}/.tests/data/{true_type}/"
            for msa_file in os.listdir(base_dir):
                msa_file = base_dir + msa_file
                msa = MSA(msa_file)
                guessed_type = msa.guess_data_type()
                assert guessed_type == true_type

    def test_number_of_taxa(self, dna_phylip_msa):
        assert dna_phylip_msa.number_of_taxa() == 68

    def test_number_of_sites(self, dna_phylip_msa):
        assert dna_phylip_msa.number_of_sites() == 766

    def test_get_avg_entropy(self, small_msa):
        # TODO
        pass

    def test_bollback_multinomial(self, small_msa):
        # TODO
        pass

    def test_treelikeness_score(self, small_msa):
        # TODO
        pass
