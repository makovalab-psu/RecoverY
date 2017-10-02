from collections import defaultdict
import subprocess
import os

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


def fastq_chopper(num_pieces, ip_file, tmp_dir):
    '''

    :param num_pieces:
    :param ip_file:
    :param ip_dir:
    :return: list of output files after chopping
    '''

    # if tmp directory doesn't exist, create the directory
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    # if tmp directory is non-empty, then empty it out
    file_list = [f for f in os.listdir(tmp_dir)]
    for f in file_list:
        os.remove(tmp_dir + '/' + f)

    # now we split input

    # first, check number of lines in input file
    try:
        wc_out = subprocess.check_output(["wc", "-l", str(ip_file)])
    except subprocess.CalledProcessError, e:
        print "wc stdout output:\n", e.output

    num_lines_in_ip_file = int(wc_out.strip().split(' ')[0])
    print "Before : Number of lines in ip file : ", num_lines_in_ip_file
    if num_lines_in_ip_file % num_pieces != 0:
        num_lines_in_ip_file = num_lines_in_ip_file + num_pieces - (num_lines_in_ip_file % num_pieces)
    print "After: Number of lines in ip file : ", num_lines_in_ip_file

    num_lines_in_tmp_pieces = num_lines_in_ip_file / num_pieces
    print "Before : Number of lines in each tmp piece : ", num_lines_in_tmp_pieces
    if num_lines_in_tmp_pieces % 4 != 0:
        num_lines_in_tmp_pieces = num_lines_in_tmp_pieces + 4 - (num_lines_in_tmp_pieces % 4)
        # num_lines_in_tmp_pieces = num_lines_in_tmp_pieces - (num_lines_in_tmp_pieces%4)
    print "After: Number of lines in each tmp piece : ", num_lines_in_tmp_pieces

    # num_lines_str = str(num_lines)
    cmd = ["split", "-l", str(num_lines_in_tmp_pieces), str(ip_file), tmp_dir + "/" + "piece_r1_"]
    subprocess.call(cmd)

    file_list = [f for f in os.listdir(tmp_dir)]
    return file_list
