#!/bin/sh

#uses dsk2.2.0
#run as ./run_dsk.sh FASTQ_file_to_be_kmerized

if [ $# -ne 1 ]
then
	echo "Usage: $0 <FASTQ_file_to_be_kmerized>"
	exit 1
fi


R1_fsY_reads=$1

echo $R1_fsY_reads

time dsk-v2.2.0-bin-Linux/bin/dsk -file $R1_fsY_reads -abundance-min 0 -kmer-size 25 -out R1_dsk -verbose 0

time dsk-v2.2.0-bin-Linux/bin/h5dump -y -d dsk/histogram R1_dsk.h5 | grep "^\ *[0-9]" | tr -d " " | tr -d "," | paste - - > R1_dsk.histo

time dsk-v2.2.0-bin-Linux/bin/dsk2ascii -file R1_dsk -out pre_Threshold_reference_table  -verbose 0


rm -rf dsk_output
mkdir dsk_output
mv R1* dsk_output/

rm -rf tables
mkdir tables
mv pre_Threshold_reference_table tables/


