import os.path

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceMatrix
from collections import Counter
from itertools import product
import math
import numpy as np
import random
from tempfile import NamedTemporaryFile

from custom_types import *

STATE_CHARS = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ!\"#$%&'()*+,/:;<=>@[\\]^_{|}~"
DNA_CHARS = "ATUCGMRWSYKVHDBN"
GAP_CHARS = "-?."


class MSA:
    def __init__(self, msa_file: FilePath, msa_name: str = None):
        if msa_name:
            self.msa_name = msa_name
        else:
            self.msa_name = os.path.split(msa_file)[1]

        self.msa_file = msa_file
        self.data_type = self.guess_data_type()
        with NamedTemporaryFile(mode="w") as tmpfile:
            if self.data_type == "DNA":
                self._convert_dna_msa_to_biopython_format(tmpfile)
            else:
                self._convert_aa_msa_to_biopython_format(tmpfile)

            tmpfile.flush()
            self.msa = AlignIO.read(tmpfile.name, format=self._get_file_format())

    def _convert_dna_msa_to_biopython_format(self, tmpfile: NamedTemporaryFile) -> None:
        """
        The unknonwn char in DNA MSA files for Biopython to work
        has to be "-" instead of "X" or "N" -> replace all occurences
        All "?" are gaps -> convert to "-"
        Also all "U" need to be "T"
        """
        with open(self.msa_file) as f:
            repl = f.read().replace("X", "-")
            repl = repl.replace("N", "-")
            repl = repl.replace("?", "-")
            repl = repl.replace("U", "T")
            repl = repl.upper()

        tmpfile.write(repl)

    def _convert_aa_msa_to_biopython_format(self, tmpfile: NamedTemporaryFile) -> None:
        """
        The unknonwn char in AA MSA files for Biopython to work
        All "?" are gaps -> convert to "-"
        """
        with open(self.msa_file) as f:
            repl = f.read().replace("?", "-")
            repl = repl.upper()

        tmpfile.write(repl)

    def _get_file_format(self) -> FileFormat:
        first_line = open(self.msa_file).readline().strip()

        try:
            # phylip file contains two integer numbers in the first line separated by whitespace
            _num1, _num2, *_ = first_line.split()
            _num1 = int(_num1)
            _num2 = int(_num2)
            # in case these conversions worked, the file is (most likely) in phylip format
            return "phylip-relaxed"
        except:
            # if the MSA is in fasta format, the first line should start with a '>' character
            if first_line.startswith(">"):
                return "fasta"

        raise ValueError(
            f"The file type of this MSA could not be autodetected, please check file."
        )

    def guess_data_type(self) -> DataType:
        format = self._get_file_format()
        msa_content = open(self.msa_file).readlines()

        sequence_chars = set()

        if format == "phylip-relaxed":
            # the sequences start in the second line and the schema is "taxon_name [SEQUENCE]"
            for line in msa_content[1:]:
                line = line.strip()
                if not line:
                    continue
                try:
                    # ...with some files the sequences are split into multiple lines
                    # in this case the taxon name is skipped and splitting will raise a ValueError
                    # if this happens, we simply add the rest of the block to the set
                    _, sequence = line.split(None, 1)
                except ValueError:
                    # read the line as is
                    sequence = line
                # remove whitespace and add the characters to the set of unique characters of this MSA
                sequence = sequence.strip().replace(" ", "")
                for char in set(sequence):
                    sequence_chars.add(char)
        elif format == "fasta":
            seen_taxon_name = False
            # taxon names start with a ">" followed by the sequence
            for line in msa_content:
                line = line.strip()
                if line.startswith(">"):
                    seen_taxon_name = True
                    continue
                else:
                    if seen_taxon_name:
                        for char in set(line):
                            sequence_chars.add(char)
                        seen_taxon_name = False
        else:
            raise ValueError(
                f"Unsupported MSA file format {format}. Supported formats are phylip and fasta."
            )

        # now check whether the sequence_chars contain only DNA and GAP chars or not
        if all([(c in DNA_CHARS) or (c in GAP_CHARS) for c in sequence_chars]):
            return "DNA"
        else:
            return "AA"

    def get_number_of_taxa(self) -> int:
        return len(self.msa)

    def get_number_of_sites(self) -> int:
        return self.msa.get_alignment_length()

    def get_column_entropies(self) -> List[float]:
        def _remove_gaps_from_sequence(seq):
            for char in GAP_CHARS:
                seq = seq.replace(char, "")
            return seq

        entropies = []
        for i in range(self.msa.get_alignment_length()):
            column = _remove_gaps_from_sequence(self.msa[:, i]).upper()
            entropy = 0

            for char in STATE_CHARS:
                count = str.count(column, char)
                if count == 0:
                    entropy_x = 0
                else:
                    prob = count / len(column)
                    entropy_x = prob * math.log2(prob)

                entropy += entropy_x

            entropy = -entropy

            assert (
                entropy >= 0
            ), f"Entropy negative, check computation. Entropy is {entropy}"

            entropies.append(entropy)
        return entropies

    def get_avg_entropy(self) -> float:
        return np.mean(self.get_column_entropies())

    def bollback_multinomial(self) -> float:
        """
        Compute the bollback multinomial statistic on the msa file
        According to Bollback, JP: Bayesian model adequacy and choice in phylogenetics (2002)
        """
        msa_length = self.get_number_of_sites()

        sites = []
        for i in range(msa_length):
            sites.append(self.msa[:, i])

        site_counts = Counter(sites)
        mult = 0
        for i in site_counts:
            N_i = site_counts[i]
            mult += N_i * math.log(N_i)

        mult = mult - msa_length * math.log(msa_length)
        return mult

    def _get_distance_matrix(self, num_samples: int) -> DistanceMatrix:
        """
        For large MSAs (i.e. more than num_samples taxa), computing the distance matrix
        is computationally very expensive.
        So for large MSAs, we rather compute the distance matrix on a subsample of at most num_samples sequences
        """
        if self.get_number_of_taxa() > num_samples:
            sample_population = range(self.get_number_of_taxa())
            selection = sorted(random.sample(sample_population, num_samples))
            _msa = MultipleSeqAlignment([self.msa[el] for el in selection])
        else:
            _msa = self.msa

        model = "blastn" if self.data_type == "DNA" else "blosum62"
        calculator = DistanceCalculator(model=model)
        return calculator.get_distance(_msa)

    def treelikeness_score(self, n_samples: int = 100) -> float:
        """
        Compute the treelikeness score according to
        δ Plots: A Tool for Analyzing Phylogenetic Distance Data, Holland, Huber, Dress and Moulton (2002)
        https://doi.org/10.1093/oxfordjournals.molbev.a004030
        """
        num_samples = min(self.get_number_of_taxa(), n_samples)
        dm = self._get_distance_matrix(num_samples)

        options = list(range(len(dm)))

        frac = num_samples // 4
        X = options[:frac]
        Y = options[frac : 2 * frac]
        U = options[2 * frac : 3 * frac]
        V = options[3 * frac :]

        res = product(X, Y, U, V)
        deltas = []
        for x, y, u, v in res:
            dxv = abs(dm[x, v])
            dyu = abs(dm[y, u])
            dxu = abs(dm[x, u])
            dyv = abs(dm[y, v])
            dxy = abs(dm[x, y])
            duv = abs(dm[u, v])

            dxv_yu = dxv + dyu
            dxu_yv = dxu + dyv
            dxy_uv = dxy + duv

            vals = sorted([dxv_yu, dxu_yv, dxy_uv])
            smallest = vals[0]
            intermediate = vals[1]
            largest = vals[2]

            numerator = largest - intermediate
            denominator = largest - smallest

            if denominator == 0:
                delta = 0
            else:
                delta = numerator / denominator
            assert delta >= 0
            assert delta <= 1
            deltas.append(delta)

        return np.mean(deltas)

    def get_character_frequencies(self) -> Dict[str, int]:
        counts = {}

        for char in STATE_CHARS:
            counts[char] = 0

        for char in GAP_CHARS:
            counts[char] = 0

        for sequence in self.msa:
            sequence = sequence.seq

            for char in STATE_CHARS:
                counts[char] += sequence.count(char)

            for char in GAP_CHARS:
                counts[char] += sequence.count(char)

        return counts