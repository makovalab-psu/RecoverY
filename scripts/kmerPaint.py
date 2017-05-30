from __future__ import division
import kmers
from collections import defaultdict
import numpy as np
import os



def make_new_kmer_table(old_kmer_table_fp, threshold, new_kmer_table_fp):
    for line in old_kmer_table_fp:
        abundance = str(line).split(' ')[1]
        if int(abundance) > threshold:
            new_kmer_table_fp.write(line)


def kmerPaint(kmer_size=25):
    """
    # This program paints histograms
    # It needs as input two files
    # file1 : fsY_kmer_table : the DSK output of fsY reads
    # file2 : trusted_kmers : DSK output of trusted single copy genes
    output : creates a new_kmer_table
    """
    ip_dir = "data"
    trusted_kmers = ip_dir + "/trusted_kmers"
    reads_kmers = ip_dir + "/kmers_from_reads"
    reads_Ymers = ip_dir + "/Ymer_table"

    op_file = ip_dir + "/trusted_DSK_counts_acc_to_fsY"

    # if file already exists, remove it
    if os.path.isfile(op_file):
        os.remove(op_file)

    print "Started kmerPaint"

    print "Creating a set from trusted kmers"
    trusted_kmers_set = kmers.make_set_from_kmer_abundance(trusted_kmers, kmer_size)

    print "Creating a dict from reads_kmers"
    fsY_kmers_dict = kmers.make_dict_from_kmer_abundance(reads_kmers, kmer_size)

    print "Done, now going to find abundance of trusted kmers acc. to all reads_kmers"
    trusted_dict_from_fsY_dict = defaultdict(int)
    for trusted in trusted_kmers_set:
        trusted_dict_from_fsY_dict[trusted] = fsY_kmers_dict[trusted]
    print "Found abundances for all trusted kmers"

    print "Finding the 5% threshold"
    all_abundances = [v for _, v in trusted_dict_from_fsY_dict.iteritems() if v > 0]
    all_abundances = np.array(all_abundances)
    threshold = np.percentile(all_abundances, 5)
    print "The 5% threshold is :", threshold

    print "Creating an output file with new kmer_counts"
    # OPTIONAL : create an output file with new kmerCounts
    with open(op_file, "w") as output_handle:
        for k, v in trusted_dict_from_fsY_dict.iteritems():
            to_write = str(k) + ' ' + str(v) + '\n'
            output_handle.write(to_write)

    print "Making the Ymer_table"
    # make the post-kmerPaint Ymer_table

    with open(reads_kmers) as reads_kmers_fp, open(reads_Ymers, "w") as reads_Ymers_fp:
        make_new_kmer_table(reads_kmers_fp, threshold, reads_Ymers_fp)


# kmerPaint()
