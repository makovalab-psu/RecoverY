from scripts import kmerPaint, classify_as_Y_chr, find_mates, kmers
import argparse
import multiprocessing as mp
import os
import shutil
import glob
from functools import partial
import sys

def main():
    """
    The following input files are required in ./data folder
    r1.fastq, r2.fastq, kmers_from_reads, trusted_kmers
    """
    
    parser = argparse.ArgumentParser(description='RecoverY selects Y-specific reads from an enriched data set.')
    parser.add_argument('--read_length', help='Set read length (defaults to 150)', required=False)
    parser.add_argument('--kmer_size', help='Set kmer size (defaults to 25)', required=False)
    parser.add_argument('--Ymer_match_threshold', help='Set Y-mer match threshold (default is calculated by formula '
                                                       ': 0.4(l-k+1-2kl/100))', required=False)
    parser.add_argument('--abundance_threshold', help='Set abundance threshold for classifying kmers as Y-mers '
                                                      '(default is calculated by using Trusted k-mers)', required=False)
    parser.add_argument('--threads', help='Set number of threads for RecoverY (defaults to 2)', required=False)
    parser.add_argument('--plots', help='Use this to generate k-mer abundance plots if you have '
                                        'matplotlib & seaborn (defaults to False)', action='store_true', required=False)
    args = vars(parser.parse_args())

    # set read_len from argument or using default here
    if not args['read_length']:
        read_len = 150
    else:
        try:
            read_len = int(args['read_length'])
        except ValueError:
            print "Error : read_length provided is not an integer"
            kmers.exit_gracefully()

    # set kmer_size from argument or using default here
    if not args['kmer_size']:
        k_size = 25
    else:
        try:
            # check if k-mer size provided by user is an integer
            k_size = int(args['kmer_size'])
            try:
                # check if kmers_from_reads file exists
                test_open = open("data/kmers_from_reads")
                first_line = test_open.readline()
                dsk_kmer_size = len(first_line.strip().split(' ')[0])
                if k_size != dsk_kmer_size :
                    print "Error : kmer_size provided is not the same as DSK kmer_size"
                    kmers.exit_gracefully()
                test_open.close()
            except IOError:
                # kmers_from_reads file does not exist
                print "Unable to locate reads_from_kmers file. " \
                      "Please check /data folder and uncompress tar file if provided."
                kmers.exit_gracefully()
        except ValueError:
            # k-mer size provided by user is not an integer!
            print "Error : kmer_size provided is not an integer"
            kmers.exit_gracefully()

    # check if kmer_size is greater than read_length
    if read_len-k_size < 0 :
        print "Error : kmer_size provided is larger than read_length"
        kmers.exit_gracefully()

    # set match_threshold from argument or using default here
    if not args['Ymer_match_threshold']:
        strictness = int(0.4 * (read_len - k_size + 1 - (2*k_size*read_len/100)))
    else:
        try :
            strictness = int(args['Ymer_match_threshold'])
        except ValueError:
            print "Error : Ymer_match_threshold provided is not an integer"
            kmers.exit_gracefully()

    # set num_threads from argument or using default here
    if not args['threads']:
        num_threads = 2
    else:
        try:
            num_threads = int(args['threads'])
        except ValueError:
            print "Error : threads provided is not an integer"
            kmers.exit_gracefully()

    print "RecoverY starting with : "
    print "number of processors : ", num_threads
    print "read length : ", read_len
    print "kmer-size : ", k_size
    print "Y-mer match threshold : ", strictness

    # manage input files
    ip_dir = "data"
    ip_r1_file_name = "r1.fastq"
    ip_file_r1 = ip_dir + "/" + ip_r1_file_name
    ip_r2_file_name = "r2.fastq"
    ip_file_r2 = ip_dir + "/" + ip_r2_file_name
    ip_kmers_from_reads_file_name = "kmers_from_reads"
    ip_file_kmers_from_reads = ip_dir + "/" + ip_kmers_from_reads_file_name
    ip_trusted_kmers_file_name = "trusted_kmers"
    ip_file_trusted_kmers = ip_dir + "/" + ip_trusted_kmers_file_name
    all_ip_files = [ip_r1_file_name, ip_r2_file_name, ip_kmers_from_reads_file_name, ip_trusted_kmers_file_name]

    # if folder "data" doesn't exist, throw an error
    if not os.path.exists(ip_dir):
        print "Input directory does not exist"
        kmers.exit_gracefully()

    # if any of the input files don't exist, throw an error
    files_provided_in_ip_folder = [f for f in os.listdir(ip_dir)]
    for file_name in all_ip_files :
        if file_name not in files_provided_in_ip_folder :
            print "Missing file : ", file_name, " please provide."
            kmers.exit_gracefully()

    # manage output files
    op_dir = "output"
    op_r1_file_name = "op_r1.fastq"
    op_file_r1 = op_dir + "/" + op_r1_file_name
    op_r2_file_name = "op_r2.fastq"
    op_file_r2 = op_dir + "/" + op_r2_file_name
    
    # if folder "output" doesn't exist, create the folder
    if not os.path.exists(op_dir):
        os.makedirs(op_dir)
    # if file output/op_r1.fastq exists, remove it
    if os.path.isfile(op_file_r1):
        os.remove(op_file_r1)
    if os.path.isfile(op_file_r2):
        os.remove(op_file_r2)

    # manage tmp files
    op_tmp_dir = "op_tmp_r1_pieces_after_classify"

    # if folder "tmp_output" doesn't exist, create the folder
    if not os.path.exists(op_tmp_dir):
        os.makedirs(op_tmp_dir)

    # if tmp directory is non-empty, then empty it out
    file_list = [f for f in os.listdir(op_tmp_dir)]
    for f in file_list:
        os.remove(op_tmp_dir + '/' + f)

    op_tmp_dir = "op_tmp_r2_pieces_after_mate_finding"
    # if folder "tmp_output" doesn't exist, create the folder
    if not os.path.exists(op_tmp_dir):
        os.makedirs(op_tmp_dir)

    # if tmp directory is non-empty, then empty it out
    file_list = [f for f in os.listdir(op_tmp_dir)]
    for f in file_list:
        os.remove(op_tmp_dir + '/' + f)

    # Ready to start RecoverY
    print "Started RecoverY"
    if not args['abundance_threshold']:
        abundance = 0
        kmerPaint.kmerPaint(ip_file_trusted_kmers, ip_file_kmers_from_reads, ip_dir, abundance, k_size)
    else:
        try:
            # check if abundance threshold provided by user is an integer
            abundance = int(args['abundance_threshold'])
            kmerPaint.kmerPaint(ip_file_trusted_kmers, ip_file_kmers_from_reads, ip_dir, abundance, k_size)
        except ValueError:
            print "Error : Abundance threshold provided is not an integer"
            kmers.exit_gracefully()

    # plot if needed
    if args['plots']:
        from scripts import plot_kmers
        print "Generating kmer plot"
        plot_kmers.plot_kmers()

    print "Chopping input R1 reads into smaller files..."
    list_of_ip_files_r1 = kmers.fastq_chopper(num_threads, ip_file_r1, "tmp_r1_pieces")

    print "Chopping input R2 reads into smaller files..."
    list_of_ip_files_r2 = kmers.fastq_chopper(num_threads, ip_file_r2, "tmp_r2_pieces")

    print "Classifying reads in parallel..."
    pool = mp.Pool(processes=num_threads)
    pool.map(partial(classify_as_Y_chr.classify_as_Y_chr, k=k_size, strict=strictness),
             [file_name for file_name in list_of_ip_files_r1])

    print "Finding mates in parallel..."
    pool = mp.Pool(processes=num_threads)
    pool.map(find_mates.find_mates, [file_name for file_name in list_of_ip_files_r2])

    print "Concatenating R1 reads..."
    with open(op_file_r1, 'wb') as outfile:
        for filename in sorted(glob.glob('./op_tmp_r1_pieces_after_classify/*')):
            with open(filename, 'rb') as readfile:
                shutil.copyfileobj(readfile, outfile)

    print "Concatenating R2 reads..."
    with open(op_file_r2, 'wb') as outfile:
        for filename in sorted(glob.glob('./op_tmp_r2_pieces_after_mate_finding/*')):
            with open(filename, 'rb') as readfile:
                shutil.copyfileobj(readfile, outfile)

    print "RecoverY completed successfully"

if __name__ == "__main__":
    main()

