import os
from Bio import SeqIO

def clean_record_id(seq_record):
    """
    cleans up the record ID if there is a "/"
    :return: a str with clean id
    """
    if "/" in seq_record.id:
        clean_id = str(seq_record.id)
        clean_id = clean_id.split('/')[0]
    else:
        clean_id = seq_record.id
    return clean_id


def find_mates():
    """
    find mates for all the post RecoverY R1 reads
    """
    ref_ip_dir = "data"
    ref_ip_file = ref_ip_dir + "/r2.fastq"
    op_dir = "output"
    r1_file = op_dir + "/op_r1.fastq"
    r2_file = op_dir + "/op_r2.fastq"

    print "Running find_mates () to find mates for fwd Y-reads"

    # if file op_r2.fastq exists, remove it
    if os.path.isfile(r2_file):
        os.remove(r2_file)

    with open(ref_ip_file) as ref_reads, open(r1_file) as r1_reads, open(r2_file, "a") as r2_reads:
        ref_record_iterator = SeqIO.parse(ref_reads, "fastq")
        for r1_record in SeqIO.parse(r1_reads, "fastq"):
            clean_r1_id = clean_record_id(r1_record)
            curr_ref_record = next(ref_record_iterator)
            clean_ref_id = clean_record_id(curr_ref_record)

            # compare R1 record and ref R2 record
            # if they are the same, then write it to output file, else keep looping
            while clean_r1_id != clean_ref_id:
                curr_ref_record = next(ref_record_iterator)
                clean_ref_id = clean_record_id(curr_ref_record)

            # complement successfully found
            SeqIO.write(curr_ref_record, r2_reads, "fastq")

# find_mates()
