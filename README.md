# compbio_pipeline_project
This pipeline is for computational biology (COMP483) Spring 2024

## Description 
Human herpesvirus 5 is also known as Human cytomegalovirus and is typically abbreviated as HCMV. Cheng et al. 2017 (https://www.ncbi.nlm.nih.gov/pubmed/29158406) sequences the transcriptomes of HCMV post infection. This pipeline uses bioinformatics tools (Bowtie2, SPAdes, and command-line BLAST) to find the most similar HCMV strain to the two patient donor HCMV transcriptomes

This code is meant for and tested in Python3.8


## Before you start

### Required Software

[Bowtie2](https://github.com/BenLangmead/bowtie2)
[SPAdes](https://github.com/ablab/spades)
[SRA](https://github.com/ncbi/sra-tools)

### Required Python packages
- os
- argparse
- glob
- re
- BioPython

###  Configure Environment 

- Command to clone this repo to your machine 
	```git clone https://github.com/haleyatkins18/compbio_pipeline_project.git```

- Then cd into pipeline_project directory 
	```cd pipeline_project```
	
### Running the pipeline

- There are two sets of data that this pipeline can run: sample set of trimmed fastq files (faster) and the full dataset of untrimmed fastq files. The results of either dataset ouput to the "PipelineProject_HA.
- To run the full dataset set run this command from your command line 
	```python gen_assemble.py```
- To run the sample dataset run this commmand from your command line
	```python loggy.py --sample```

### Output files 
- All files and results are found in the PipelineProject_HA directory

- Each of the commands generate more files, but this only includes the files that are pertinent to running the pipeline. More information on the other files can be found at the software links.
 
- blast
--Blast output results of blasting the longest contig from the SPAdes assembly against the Betaherpesvirinae subfamily

- bowtie_index 
--indexes created by bowtie

- bowtie2_output 
--output from bowtie2

- PipelineProejct.log

-spades
-- contigs.fasta contains the assembled contigs  
