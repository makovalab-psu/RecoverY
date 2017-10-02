import kmers
import os
from Bio import SeqIO


def classify_as_Y_chr(ip_file, kmer_size=25, strictness=20):
    """
    :param kmer_size: 25
    :param strictness: 20
    :return: 1 and create file R1
    """
    print "classify has been called on file : "
    Ymer_table = "data/Ymer_table"
    ip_tmp_dir = "tmp_r1_pieces/"
    op_tmp_dir = "op_tmp_r1_pieces_after_classify/"
    op_file = op_tmp_dir + ip_file + "_op_r1"
    ip_file = ip_tmp_dir + ip_file
    print ip_file


    print "Running classify_reads() to shortlist fwd Y-reads from ip dataset"

    print "First making Ymer set ... this might take a while"
    Ymer_set = kmers.make_set_from_kmer_abundance(Ymer_table, kmer_size)
    print "Ymer set ready, time to classify"

    # now that you have a fresh op_file ...
    with open(ip_file, "r") as ip_reads, open(op_file, "a") as op_reads:
        for seq_record in SeqIO.parse(ip_reads, "fastq"):
            curr_seq = str(seq_record.seq)

            kmers_from_seq = kmers.kmerize(curr_seq, kmer_size)
            matches = 0
            for kmer in kmers_from_seq:
                if kmer in Ymer_set:
                    matches += 1

            # print "# of matches is :", matches

            if matches > strictness:
                # Y_seq_records.append(seq_record)
                SeqIO.write(seq_record, op_reads, "fastq")
    print "Classification done"


# classify_as_Y_chr()
