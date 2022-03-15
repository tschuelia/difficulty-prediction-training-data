from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator
from collections import Counter
from itertools import product
import math
import numpy as np
import random
import subprocess
from tempfile import TemporaryDirectory, NamedTemporaryFile

from custom_types import FilePath

STATE_CHARS = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ!\"#$%&'()*+,/:;<=>@[\\]^_{|}~"
DNA_CHARS = "ATUCGMRWSYKVHDBN"
GAP_CHARS = "-?."


def _convert_dna_msa_to_biopython_format(msa_file: FilePath, tmpfile: FilePath) -> None:
    """
    The unknonwn char in DNA MSA files for Biopython to work
    has to be "N" instead of "X" -> replace all occurences
    All "?" are gaps -> convert to "-"
    Also all "U" need to be "T"
    """
    with open(msa_file) as f:
        repl = f.read().replace("X", "N")
        repl = repl.replace("?", "-")
        repl = repl.replace("U", "T")
        repl = repl.upper()

    tmpfile.write(repl)


def _convert_aa_msa_to_biopython_format(msa_file, tmpfile):
    """
    The unknonwn char in AA MSA files for Biopython to work
    All "?" are gaps -> convert to "-"
    """
    with open(msa_file) as f:
        repl = f.read().replace("?", "-")
        repl = repl.upper()

    tmpfile.write(repl)


def _get_msa_file_format(msa_file):
    first_line = open(msa_file).readline().strip()

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

    raise ValueError(f"The file type of this MSA could not be autodetected, please check file.")


def guess_msa_file_data_type(msa_file):
    format = _get_msa_file_format(msa_file)
    msa_content = open(msa_file).readlines()

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
        raise ValueError(f"Unsupported MSA file format {format}. Supported formats are phylip and fasta.")

    # now check whether the sequence_chars contain only DNA and GAP chars or not
    if all([(c in DNA_CHARS) or (c in GAP_CHARS) for c in sequence_chars]):
        return "DNA"
    else:
        return "AA"


def read_alignment(msa_file):
    data_type = guess_msa_file_data_type(msa_file)
    with NamedTemporaryFile(mode="w") as tmpfile:
        if data_type == "DNA":
            _convert_dna_msa_to_biopython_format(msa_file, tmpfile)
        else:
            _convert_aa_msa_to_biopython_format(msa_file, tmpfile)

        tmpfile.flush()
        return AlignIO.read(tmpfile.name, format=_get_msa_file_format(msa_file))


def get_number_of_taxa(msa):
    return len(msa)


def get_number_of_sites(msa):
    return msa.get_alignment_length()


def remove_gaps_from_sequence(seq):
    for char in GAP_CHARS:
        seq = seq.replace(char, "")
    return seq


def get_column_entropy(column):
    # column_entropy = - sum(for every nucleotide x) {count(x)*log2(Prob(nuc x in col i))}
    column = remove_gaps_from_sequence(column).upper()
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

    assert entropy >= 0, f"Entropy negative, check computation. Entropy is {entropy}"

    return entropy


def get_msa_column_entropies(msa):
    return [get_column_entropy(msa[:, i]) for i in range(msa.get_alignment_length())]


def get_msa_avg_entropy(msa):
    return np.mean(get_msa_column_entropies(msa))


def bollback_multinomial(msa):
    """
    Compute the bollback multinomial statistic on the msa file
    According to Bollback, JP: Bayesian model adequacy and choice in phylogenetics (2002)
    """
    msa_length = msa.get_alignment_length()

    sites = []
    for i in range(msa_length):
        sites.append(msa[:, i])

    site_counts = Counter(sites)
    mult = 0
    for i in site_counts:
        N_i = site_counts[i]
        mult += N_i * math.log(N_i)

    mult = mult - msa_length * math.log(msa_length)
    return mult


def _get_distance_matrix(msa, num_samples, data_type):
    """
    For large MSAs (i.e. more than num_samples taxa), computing the distance matrix
    is computationally very expensive.
    So for large MSAs, we rather compute the distance matrix on a subsample of at most num_samples sequences
    """
    if get_number_of_taxa(msa) > num_samples:
        sample_population = range(get_number_of_taxa(msa))
        selection = sorted(random.sample(sample_population, num_samples))
        msa = MultipleSeqAlignment([msa[el] for el in selection])

    model = "blastn" if data_type == "DNA" else "blosum62"
    calculator = DistanceCalculator(model=model)
    dm = calculator.get_distance(msa)
    return dm


def treelikeness_score(msa, data_type):
    """
    Compute the treelikeness score according to
    δ Plots: A Tool for Analyzing Phylogenetic Distance Data, Holland, Huber, Dress and Moulton (2002)
    https://doi.org/10.1093/oxfordjournals.molbev.a004030
    """
    num_samples = min(get_number_of_taxa(msa), 100)
    dm = _get_distance_matrix(msa, num_samples, data_type)

    options = list(range(len(dm)))

    frac = num_samples // 4
    X = options[:frac]
    Y = options[frac:2 * frac]
    U = options[2 * frac:3 * frac]
    V = options[3 * frac:]

    res = product(X, Y, U, V)
    deltas = []
    for x, y, u, v in res:
        dxv = np.abs(dm[x, v])
        dyu = np.abs(dm[y, u])
        dxu = np.abs(dm[x, u])
        dyv = np.abs(dm[y, v])
        dxy = np.abs(dm[x, y])
        duv = np.abs(dm[u, v])

        # assert dm[x, v] >= 0, (x, v, dm[x, v])
        # assert dm[y, u] >= 0, (y, u, dm[y, u])
        # assert dm[x, u] >= 0, (x, u, dm[x, u])
        # assert dm[y, v] >= 0, (y, v, dm[y, v])
        # assert dm[x, y] >= 0, (x, y, dm[x, y])
        # assert dm[u, v] >= 0, (u, v, dm[u, v])

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


def _run_raxmlng_alignment_parse(msa_file, raxmlng_executable, model, tmpdir):
    cmd = [
        raxmlng_executable,
        "--parse",
        "--msa",
        msa_file,
        "--model",
        model,
        "--prefix",
        tmpdir
    ]

    subprocess.check_output(cmd)


def get_number_of_patterns(msa_file, raxmlng_executable="raxml-ng", model="GTR+G"):
    """
    The easiest way to get the number of patterns for the msa_file is to run
    raxml-ng with the --parse option and parse the results.
    """
    with TemporaryDirectory() as tmpdir:
        _run_raxmlng_alignment_parse(msa_file, raxmlng_executable, model, tmpdir)

        log_file = tmpdir + ".raxml.log"

        for line in open(log_file).readlines():
            line = line.strip()
            if not line.startswith("Alignment sites"):
                continue
            # Alignment sites / patterns: 1940 / 933
            _, numbers = line.split(":")
            sites, patterns = numbers.split("/")
            return int(patterns)


def get_percentage_of_invariant_sites(msa_file, raxmlng_executable="raxml-ng", model="GTR+G"):
    """
    The easiest way to get the number of patterns for the msa_file is to run
    raxml-ng with the --parse option and parse the results.
    """
    with TemporaryDirectory() as tmpdir:
        _run_raxmlng_alignment_parse(msa_file, raxmlng_executable, model, tmpdir)

        log_file = tmpdir + ".raxml.log"

        for line in open(log_file).readlines():
            line = line.strip()
            if not line.startswith("Invariant sites"):
                continue
            # Invariant sites: 80.77 %
            _, number = line.split(":")
            percentage, _ = number.strip().split(" ")
            return float(percentage) / 100.0


def get_percentage_of_gaps(msa_file, raxmlng_executable="raxml-ng", model="GTR+G"):
    """
    The easiest way to get the number of patterns for the msa_file is to run
    raxml-ng with the --parse option and parse the results.
    """
    with TemporaryDirectory() as tmpdir:
        _run_raxmlng_alignment_parse(msa_file, raxmlng_executable, model, tmpdir)

        log_file = tmpdir + ".raxml.log"

        for line in open(log_file).readlines():
            line = line.strip()
            if not line.startswith("Gaps"):
                continue
            # Gaps: 20.05 %
            _, number = line.split(":")
            percentage, _ = number.strip().split(" ")
            return float(percentage) / 100.0


def get_character_frequencies(msa):
    counts = {}

    for char in STATE_CHARS:
        counts[char] = 0

    for sequence in msa:
        sequence = sequence.seq

        for char in STATE_CHARS:
            counts[char] += sequence.count(char)

    return counts

