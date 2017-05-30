import kmers
import os
from Bio import SeqIO


def classify_as_Y_chr(kmer_size=25, strictness=20):
    """
    :param kmer_size: 25
    :param strictness: 20
    :return: 1 and create file R1
    """

    ip_dir = "data"
    ip_file = ip_dir + "/r1.fastq"
    Ymer_table = ip_dir + "/Ymer_table"
    op_dir = "output"
    op_file = op_dir + "/op_r1.fastq"

    print "Running classify_reads() to shortlist fwd Y-reads from ip dataset"

    print "First making Ymer set ... this might take a while"
    Ymer_set = kmers.make_set_from_kmer_abundance(Ymer_table, kmer_size)
    print "Ymer set ready, time to classify"

    # if folder "output" doesn't exist, create the folder
    if not os.path.exists(op_dir):
        os.makedirs(op_dir)

    # if file op_r1.fastq exists, remove it
    if os.path.isfile(op_file):
        os.remove(op_file)

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
