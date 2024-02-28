# compbio_pipeline_project
This pipeline is for the computational biology (COMP483) Spring 2024

## Description 
Human herpesvirus 5 is also known as Human cytomegalovirus and is typically abbreviated as HCMV. Cheng et al. 2017 (https://www.ncbi.nlm.nih.gov/pubmed/29158406) sequences the transcriptomes of HCMV post infection. This pipeline uses bioinformatics tools (Bowtie2, SPAdes, and command-line BLAST) to find the most similar HCMV strain to the two patient donor HCMV transcriptomes

This code is meant for and test in Python3.8


## Before you start

### Required Software


### Required Python packages
- os
- argparse
- glob
- re
- BioPython

###  Configure Environment 

- Command to clone this repo to your machine 
	`git clone https://github.com/haleyatkins18/compbio_pipeline_project.git`

- Then cd into pipeline_project directory 
	`cd pipeline_project`
	
### Running the pipeline

- There are two sets of data that this pipeline can run: sample set of trimmed fastq files (faster) and the full dataset of untrimmed fastq files. The results of either dataset ouput to the "PipelineProject_HA.
- To run the full dataset set run this command from your command line 
python loggy.py
- To run the sample dataset run this commmand from your command line
python loggy.py --sample

### Output files 
- All files and results are found in the PipelineProject_HA directory
- blast
--Blast output results of blasting the longest contig from the SPAdes assembly against the Betaherpesvirinae subfamily

- bowtie_index 
--indexes created by bowtie

- bowtie2_output 
--output from bowtie2

- PipelineProejct.log

-spades
--assembly files 
