# compbio_pipeline_project
This pipeline is for computational biology (COMP483) Spring 2024

## Description 
Human herpesvirus 5 is also known as Human cytomegalovirus and is typically abbreviated as HCMV. Cheng et al. 2017 (https://www.ncbi.nlm.nih.gov/pubmed/29158406) sequences the transcriptomes of HCMV post infection. The goal of this project is to compare HCMV transcriptomes 2- and 6-days post-infection (dpi). This pipeline follows the genome assembly approach using Bowtie2, SPAdes, and command-line BLAST to find the most similar HCMV strain to the two patient donor HCMV transcriptomes

This code is meant for and tested in Python 3.8


## Before you start

### Required Software

- [Bowtie2](https://github.com/BenLangmead/bowtie2)
- [SPAdes](https://github.com/ablab/spades)
- [SRA](https://github.com/ncbi/sra-tools)

### Required Python packages
- os
- argparse
- glob
- re
- BioPython

###  Configure Environment 

- Command to clone this repo to your machine 
	-  ```git clone https://github.com/haleyatkins18/compbio_pipeline_project.git```

- Then cd into pipeline_project directory 
	- ```cd pipeline_project```
	
### Running the pipeline

- Now that this github is cloned and your environment is configured correctly, you can actually run this pipeline. 
- There are two sets of data that this pipeline can run: sample set of trimmed fastq files and the full dataset of untrimmed fastq files. The results of either dataset ouput to the "PipelineProject_HA.
- To run the full dataset set run this command from your command line 
	- ```python gen_assemble.py```
- To run the sample dataset run this commmand from your command line
	- ```python loggy.py --sample```
- The sample dataset was created using the comand:
  	- ``` head -n 40000 data.fastq > sampledata.fastq``` and then creating the new, sample fastq files. This is usedul beacuase it takes less time to run these trimmed samples rather than the full dataset.

- Looking at the log file
  	- ```cat PipelineProject_HA.log```

### Output files 
- All files and results are found in the PipelineProject_HA directory

- Each of the commands generate more files, but this only includes the files that are pertinent to running the pipeline. More information on the other files can be found at the software links.
 
- blast
	-  Blast output results of blasting the longest contig from the SPAdes assembly against the Betaherpesvirinae subfamily

- bowtie_index 
	- Index.fasta contains all the sequences you want to index
   	- *.bt2 files contain the genome indeces, epresenting different components of the index
- bowtie2_output 
	- Mapped fastq files for reach paired read
   	- Sam files are the human-readable text files that contain alignment information for each read
 
- PipelineProject.log
  	- Contains all of the desired output
  	- Number of reads in each transcriptome before and after the Bowtie2 mapping
  	-  SPAdes command used
  	-  The number of contigs with a length > 1000
  	-  The length of the assembly (the total number of bp in all of the contigs > 1000 bp in length)
  	-  Subject accession,
	- Percent identity, Alignment length, Start of alignment in query, End of alignment in query, Start of alignment in subject, End of alignment in subject, Bit score, E-value, and Subject Title of the top ten BLAST hits

- spades
	- Contigs.fasta contains the assembled contigs  
