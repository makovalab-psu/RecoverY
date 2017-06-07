from scripts import kmerPaint, classify_as_Y_chr, find_mates, plot_kmers

def main():
    """
    The following input files are required in ./data folder
    r1.fastq, r2.fastq, kmers_from_reads, trusted_kmers
    """
    print "Started RecoverY"
    #print "Please specify kmer-size, strictness and data directory"

    # declare defaults
    print "Using default of k=25, strictness=20, and input folder='data'"

    print "Running kmerPaint"
    kmerPaint.kmerPaint()

    # Un-comment the lines below if you have matplotlib and seaborn
    # print "Generating kmer plot"
    # plot_kmers.plot_kmers()

    print "Shortlisting Y-reads"
    classify_as_Y_chr.classify_as_Y_chr()

    print "Finding mates"
    find_mates.find_mates()

    print "RecoverY completed successfully"

if __name__ == "__main__":
    main()

