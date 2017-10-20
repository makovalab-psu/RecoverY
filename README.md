# RecoverY

RecoverY is a tool for shortlisting enriched reads from a sequencing dataset, based on k-mer abundance. Specifically, it can be used for isolating Y-specific reads from a Y flow-sorted dataset.

## Usage 

Before running RecoverY, the following input files are required in ./data folder.

	r1.fastq : Enriched raw reads (first in pair) 
	r2.fastq : Enriched raw reads (second in pair) 
	kmers_from_reads : kmer counts for r1.fastq
	trusted_kmers : kmer counts for human Y single copy genes

See below for a description of the file formats. Note that currently the names of the "data" folder and of the files are hardcoded into RecoverY.

To run RecoverY, 

	cd RecoverY
	python recoverY.py [--help] [--read_length READ_LENGTH] [--kmer_size KMER_SIZE]
                   [--Ymer_match_threshold YMER_MATCH_THRESHOLD]
                   [--threads THREADS] [--plots]

Important parameters for the user to choose are : 

- --kmer\_size (default: 25): kmer size used for classifying reads. This must be the same as DSK's k-mer-size. We recommend a value between 25 and 31 for Illumina 150x150 bp reads. This value is used for calculating Ymer\_match\_threshold below

- --read\_length (default:150) : The length of longest un-trimmed read input to RecoverY. This value is used for calculating Ymer\_match\_threshold below
	
- --Ymer\_match\_threshold : the number of k-mers a read must match to the Ymer table in order to be classified as coming from the Y. The default value is calculated by the formula : 0.4 * (read_len - kmer_size + 1 - (2\*kmer\_size\*read\_len/100)). User may change this, but we recommend a value between 20 and 50 for Illumina 150x150 bp reads.
	
- --threads (default: 2): The number of threads that RecoverY should use. 

- --help: print usage information.

- --plots: generates a k-mer abundance plot in ```output/kmerplot.png``` (requires Matplotlib and Seaborn python packages to be installed, see below). This plot visualizes the abundance threshold selected by RecoverY, by plotting the abundance of raw read k-mers as well as trusted k-mers. See [here](/img/kmerplot.png) for an example. 

For example, you can run RecoverY as follows. 

	python recoverY.py --read_length 250 --k_size 31 --Ymer_match_threshold 50 --threads 8
	
The output of RecoverY are two files: ```output/op_r1.fastq``` and ```output/op_r2.fastq```. 
This is a subset of the read pairs in ```data/r1.fastq``` and ```data/r2.fastq``` that are deemed to have originated on the Y chromosome. 


### Trusted kmers 

The trusted_kmers file consists of (kmer_sequence,count) pairs separated by a whitespace. Below is an example of the lexicographically smallest 5 trusted k-mers :
	
	AAAAAAAAAAAAAAAAGAAAAACAA 1
	AAAAAAAAAAAAAAACAAGCTGAAT 1
	AAAAAAAAAAAAAAAGAAAAACAAA 1
	AAAAAAAAAAAAAACAAGCTGAATG 1
	AAAAAAAAAAAAAAGAAAAACAAAA 1


The trusted_kmers file consists of all the single copy k-mers obtained by k-merizing known single copy regions on the Y. 
These trusted k-mers are used as a proxy to determine the abundance threshold for Y-mers. 
This file is obtained by running DSK on X-degenerate gene sequences from human Y, and from these, extracting only k-mers with a count of 1. 

### k-mers from reads file
TODO

## Installation 

To download, 

	git clone https://github.com/makovalab-psu/RecoverY
	
RecoverY also requires the numpy and biopython python packages in order to run.
Additionally, if you use the ```--plot``` option, RecoverY needs the matplotlib and seaborn packages for python.
However, RecoverY can be run without the matplotlib or seaborn packages, as long as the ```--plot``` option is not used.
These packages can be installed on many  systems as follows:

    pip install numpy
    pip install biopython
    pip install matplotlib
    pip install seaborn

RecoverY also uses the k-mer counter DSK. The latest DSK binaries (v2.2.0 for Linux 64 bit and v2.2.0 for Mac OSX) are provided in the dependency folder. Thus, if you are using either of these operating systems, DSK need not be installed, and you may use the binaries as provided. For other operating systems, or if alternate versions or functionality of DSK is desired, see https://gatb.inria.fr/software/dsk/.

## Example

The data folder contains an example reads dataset and kmer tables. 
It can be used to test if RecoverY runs to completion. 

Before running recoverY.py, please navigate to the data folder and un-compress the tar.xz file : 

	cd data/
	tar xf kmers_from_reads.tar.xz

Subsequently, RecoverY can be run as : 

	cd ../
	python recoverY.py
	
**Results :**

The data/r1.fastq and data/r2.fastq were generated from hg38 using wg-sim.
Thus, each FASTQ record header has the chromosome of origin for a given read. 

Using grep and wc commands, one can check if RecoverY has correctly retrieved most of the Y-reads. 

	grep "@chrY" data/r1.fastq | wc -l
	grep "@chrY" output/op_r1.fastq | wc -l


### Generating k-mer counts with DSK

The ./dependency folder contains DSK binaries and a script that helps generate k-mer counts required for RecoverY. There are separate binaries for Linux 64 bit and Mac OSX. Usage is as follows (example shown below is for a Linux system) :

    cd dependency
    ./run_dsk_Linux.sh <FASTQ_file> <kmer_size>


If the k-mer counts file for raw reads (r1.fastq) is not already provided, the user may need to generate k-mer counts manually using DSK. To generate k-mer counts with DSK, the following steps are needed : 

    cd dependency 
    ln -s ../data/r1.fastq   # make sure the correct reads file is provided to DSK
    ./run_dsk_Linux.sh r1.fastq 25  


The kmer\_counts table will be generated in :

    dependency/dsk_output/kmers_from_reads


This file can be copied or linked to the data folder so that RecoverY can use it : 

    cd ../data
    ln -s ../dependency/kmers_from_reads 



## Scripts 

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
If you use RecoverY in your research, please cite 

[RecoverY : K-mer based read classification for Y-chromosome specific sequencing and assembly](https://doi.org/10.1101/148114), 

Samarth Rangavittal, Robert S. Harris, Monika Cechova, Marta Tomaszkiewicz, Rayan Chikhi, Kateryna Makova, Paul Medvedev

bioRxiv 2017.
