from collections import defaultdict


def reverse_complement(seq):
    seq_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return "".join([seq_dict[base] for base in reversed(seq)])


def kmerize(ip_string, kmer_size):
    """
    function that kmerizes an input_string and returns a list of kmers
    """
    return [ip_string[i:i + kmer_size] for i in range(0, len(ip_string) - kmer_size + 1, 1)]


def make_set_from_kmer_abundance(ip_file, kmer_size):
    """
    function that settifies a large file with kmers and count
    file looks like : ATGCT 101
    """
    kmers = []
    with open(ip_file, 'r') as kmers_fp:
        for line in kmers_fp:
            kmers.append(line[:kmer_size])
            kmers.append(reverse_complement(line[:kmer_size]))
    return set(kmers)


def make_dict_from_kmer_abundance (ip_file, kmer_size) :
    kmer_dicts = defaultdict(int)
    with open(ip_file,'r') as file_handle1 :
        for line in file_handle1 :
            current_abundance = int(line.split(' ')[1])
            kmer_dicts[line[:kmer_size]] = current_abundance
            # kmer_dicts[reverse_complement(line[:kmer_size])] = current_abundance
    return kmer_dicts
