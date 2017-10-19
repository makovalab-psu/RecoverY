from scripts import kmerPaint, classify_as_Y_chr, find_mates, kmers
import argparse
# un-comment line below if you have installed matplotlib
# from scripts import plot_kmers
import multiprocessing as mp
import os
import shutil
import glob

def main():
    """
    The following input files are required in ./data folder
    r1.fastq, r2.fastq, kmers_from_reads, trusted_kmers
    """
    
    parser = argparse.ArgumentParser(description='RecoverY selects Y-specific reads from an enriched data set')
    parser.add_argument('--threads', help='Set number of threads for RecoverY (defaults to 2)', required=False)
    args = vars(parser.parse_args())

    # set num_threads from argument or using default here
    if not args['threads']:
        num_threads = 2
    else:
        num_threads = int(args['threads'])

    print "RecoverY starting with number of processors : ", num_threads 
    
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

    op_tmp_dir = "op_tmp_r1_pieces_after_classify"
    # if folder "output" doesn't exist, create the folder
    if not os.path.exists(op_tmp_dir):
        os.makedirs(op_tmp_dir)

    # if tmp directory is non-empty, then empty it out
    file_list = [f for f in os.listdir(op_tmp_dir)]
    for f in file_list:
        os.remove(op_tmp_dir + '/' + f)


    op_tmp_dir = "op_tmp_r2_pieces_after_mate_finding"
    # if folder "output" doesn't exist, create the folder
    if not os.path.exists(op_tmp_dir):
        os.makedirs(op_tmp_dir)

    # if tmp directory is non-empty, then empty it out
    file_list = [f for f in os.listdir(op_tmp_dir)]
    for f in file_list:
        os.remove(op_tmp_dir + '/' + f)


    print "Started RecoverY"
    #print "Please specify kmer-size, strictness and data directory"

    # declare defaults
    print "Using default of k=25, strictness=20, and input folder='data'"

    #print "Running kmerPaint"
    kmerPaint.kmerPaint()

    # Un-comment the lines below if you have matplotlib and seaborn
    #print "Generating kmer plot"
    #plot_kmers.plot_kmers()

    print "Chopping input R1 reads into smaller files..."
    list_of_ip_files_r1 = kmers.fastq_chopper(num_threads, "data/r1.fastq", "tmp_r1_pieces")
    for file_name in list_of_ip_files_r1 :
        print file_name

    print "Chopping input R2 reads into smaller files..."
    list_of_ip_files_r2 = kmers.fastq_chopper(num_threads, "data/r2.fastq", "tmp_r2_pieces")
    for file_name in list_of_ip_files_r2 :
        print file_name

    print "Classifying reads in parallel..."
    pool = mp.Pool(processes=num_threads)
    pool.map(classify_as_Y_chr.classify_as_Y_chr, [file_name for file_name in list_of_ip_files_r1])

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

