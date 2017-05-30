from scripts import kmerPaint, classify_as_Y_chr, find_mates

def main():
    """
    The following input files are required in ./data folder
    r1.fastq, r2.fastq, kmers_from_reads, trusted_kmers
    """
    print "Going to run RecoverY..."
    print "Please specify kmer-size, strictness and data directory"

    # wait for user input

    # user does not input, so use default
    print "Using default of k=25, strictness=20, and folder='data'"

    print "Going to run kmerPaint"
    kmerPaint.kmerPaint()

    print "Going to shortlist Y-reads"
    classify_as_Y_chr.classify_as_Y_chr()

    print "Going to find mates"
    find_mates.find_mates()

    print "RecoverY completed successfully"

if __name__ == "__main__":
    main()

