from __future__ import division
import kmers
from collections import defaultdict
import numpy as np
import os
import sys


def make_new_kmer_table(old_kmer_table_fp, threshold, new_kmer_table_fp):
    for line in old_kmer_table_fp:
        abundance = str(line).split(' ')[1]
        if int(abundance) > threshold:
            new_kmer_table_fp.write(line)


def kmerPaint(trusted_kmers, reads_kmers, ip_dir, abundance_threshold, kmer_size=25):
    """
    :param trusted_kmers:
    :param reads_kmers:
    :param ip_dir:
    :param abundance_threshold:
    :param kmer_size:
    :return:
    """
    reads_Ymers = ip_dir + "/Ymer_table"
    op_file = ip_dir + "/trusted_DSK_counts_acc_to_fsY"

    # if file already exists, remove it
    if os.path.isfile(op_file):
        os.remove(op_file)

    # if abundance threshold has not been changed from default, then calculate it
    if abundance_threshold == 0:
        print "Testing validity of trusted_kmers file and reading into memory... "
        trusted_kmers_set = kmers.make_set_from_kmer_abundance(trusted_kmers, kmer_size)
        print "Trusted_kmers file is valid"

        print "Testing validity of reads_kmers file and reading into memory... "
        fsY_kmers_dict = kmers.make_dict_from_kmer_abundance(reads_kmers, kmer_size)
        print "Kmers_from_reads file is valid"

        print "Done, now going to find abundance of trusted kmers acc. to all reads_kmers"
        trusted_dict_from_fsY_dict = defaultdict(int)
        for trusted in trusted_kmers_set:
            trusted_dict_from_fsY_dict[trusted] = fsY_kmers_dict[trusted]

        all_abundances = [v for _, v in trusted_dict_from_fsY_dict.iteritems() if v > 0]
        if not all_abundances:
            print "Error : None of the trusted_kmers are found in kmers_from_reads"
            print "Please check if trusted_kmers file is from an identical or related species, " \
                  "or if read_kmers_file is too small"
            kmers.exit_gracefully()
        all_abundances = np.array(all_abundances)
        abundance_threshold = np.percentile(all_abundances, 5)
        print "The 5% threshold calculated using trusted kmers is :", abundance_threshold
    else:
        print "The abundance threshold set by user is :", abundance_threshold

    # OPTIONAL : create an output file with new kmerCounts
    '''
    with open(op_file, "w") as output_handle:
        for k, v in trusted_dict_from_fsY_dict.iteritems():
            to_write = str(k) + ' ' + str(v) + '\n'
            output_handle.write(to_write)
    '''

    print "Making the Ymer_table"
    with open(reads_kmers) as reads_kmers_fp, open(reads_Ymers, "w") as reads_Ymers_fp:
        make_new_kmer_table(reads_kmers_fp, abundance_threshold, reads_Ymers_fp)


# kmerPaint()
