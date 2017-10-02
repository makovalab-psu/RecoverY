# RecoverY

RecoverY is a tool for shortlisting enriched reads from a sequencing dataset, based on k-mer abundance. Specifically, it can be used for isolating Y-specific reads from a Y flow-sorted dataset.

### Usage  

    python recoverY.py
  	
Important parameters for the user to choose are : 


**kmer-size** : 
- the size of k used while iterating through every read 
- this must be the same as DSK's kmer-size
- usually optimal in the range [25, 31] for Illumina 150x150 bp reads


**strictness** : 
- the # of successful matches to the Ymer table required per read, before classifying a read as Y-specific 
- usually optimal in the range [20, 50] for Illumina 150x150 bp reads



### Installation 

	git clone https://github.com/makovalab-psu/RecoverY
	cd RecoverY


### Dependencies 

The latest DSK binary (v2.2.0 for Linux) is provided in the dependency folder. 
See https://gatb.inria.fr/software/dsk/ for alternate versions. 
    
Numpy and Biopython can be installed as follows :

    pip install numpy
    pip install biopython
    

### Input

The following input files are required in ./data folder. 
Note that currently RecoverY expects the folder to be named "data".
    	
	
	r1.fastq : Enriched raw reads (first in pair) 
	r2.fastq : Enriched raw reads (second in pair) 
	kmers_from_reads : kmer counts from DSK for r1.fastq
	trusted_kmers : kmer counts from DSK for human Y single copy genes

The input folder and file names can be changed by the user within the program. 

### Generating k-mer counts with DSK

The ./dependency folder contains a DSK binary and a script that help generate k-mer counts required for RecoverY. Usage is as follows :

    cd dependency
    ./run_dsk.sh <FASTQ_FILE>

In this case, FASTQ_FILE is r1.fastq. 
The kmer_counts table will be generated in 

    ./dependency/dsk_output/kmer_counts_from_dsk

### Output 

The ./output folder contains :

 	op_r1.fastq
	op_r2.fastq

These are the Y-reads files produced by RecoverY.  


### Example

The data folder contains an example reads dataset and kmer tables. 
It can be used to test if RecoverY runs to completion. 

Before running recoverY.py, please navigate to the data folder and un-compress the tar.xz file : 

	cd data/
	tar xf kmers_from_reads.tar.xz


### Generating k-mer plots 

Matplotlib and Seaborn are required to generate k-mer plots. 

           pip install matplotlib
           pip install seaborn

After installation, please un-comment the following line from recoverY.py :

	print "Generating kmer plot"
	plot_kmers.plot_kmers()


### Scripts 

The following scripts are included with this distribution of RecoverY, and are automatically run by recovery.py as part of the pipeline. Users may consider them separately for custom needs if required. 

	
**kmers.py** 
	
	a set of general purpose functions to work with kmers

**kmerPaint.py**
	
	input : trusted_kmers and reads_from_kmers 
	output : Ymer_table with new abundance threshold

**classify_as_Y_chr.py**
	
	input : all raw reads (first in pair) and Ymer table
	output : Y-specific reads according to RecoverY algorithm (first in pair)

**find_mates.py** 

	input : all raw reads (second in pair) and Y-specific reads accoding to RecoverY algorithm (first in pair)
	output : Y-specific reads according to RecoverY algorithm (second in pair)


### License
This program is released under the MIT License. Please see LICENSE.md for details


### Citation
Please cite this Github repository if you use this tool in your research. Thanks !
https://github.com/makovalab-psu/RecoverY
