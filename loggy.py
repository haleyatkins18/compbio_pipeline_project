#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 07:48:49 2024

@author: haleyatkins
"""

import os
import argparse
import glob
import re
from Bio import Entrez
from Bio import SeqIO
from os import popen

def convert_sra_to_fastq(main):
                                                #create directory for untrimmed_fastqs
    untrimmed_fastq_dir = os.path.join(main, "untrimmed_fastqs")
    os.makedirs(untrimmed_fastq_dir, exist_ok=True)
                                                #change directory to untrimmed_fastqs
    os.chdir(untrimmed_fastq_dir)
                                                #command for fastq-dump for each SRA file 
    fastq_dump_cmd = "fastq-dump -I --split-files"
    for sra_file in glob.glob(os.path.join(main, 'SRAfiles', '*')):
        os.system(f"{fastq_dump_cmd} {sra_file}")

def sra_files(main):                        #function to get the SRA files
                                                #create SRAfiles directory
    sra_files_dir = os.path.join(main, 'SRAfiles')
    os.makedirs(sra_files_dir, exist_ok=True)

                                                #change directory to SRAfiles
    os.chdir(sra_files_dir)

                                                #download the SRA files
    donor_urls = [
        "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030",
        "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033",
        "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044",
        "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045"
    ]
    for url in donor_urls:                      #for each of the donors in the list get the downloads from the wab
        os.system(f"wget {url}")

def get_donor(srr_no):                          #assigning each SRA to the matching donor name 
    donor_mapping = {
        "SRR5660030": "Donor 1 (2dpi)",
        "SRR5660033": "Donor 1 (6dpi)",
        "SRR5660044": "Donor 3 (2dpi)",
        "SRR5660045": "Donor 3 (6dpi)"
    }
    return donor_mapping.get(srr_no, "Unknown Donor")       #return the donor name for each SRA but if it's not the right name then return unkown 

def build_index(main):                                      #function to build the bowtie index       
                                                            #create the bowtie index directory
    bowtie_index_dir = os.path.join(main, 'bowtie_index')
    os.makedirs(bowtie_index_dir, exist_ok=True)            #make the bowtie directory 

                                                            #change directory to bowtie_index
    os.chdir(bowtie_index_dir)

                                                            #fetch HCMV sequence and build index
    handle = Entrez.efetch(db="nucleotide", id='NC_006273.2', rettype="fasta")
    records = list(SeqIO.parse(handle, "fasta"))                #parse through the fasta files 
    SeqIO.write(records, os.path.join(bowtie_index_dir, 'index.fasta'), "fasta")    #write the results to the index fasta files

                                                                                    #builds the bowtie index
    bowtie_cmd = f"bowtie2-build {os.path.join(bowtie_index_dir, 'index.fasta')} HCMV"
    os.system(bowtie_cmd)

def bowtie(main, input_reads_path):         
                                                                        #create the bowtie2_output directory
    bowtie_output_dir = os.path.join(main, 'bowtie2_output')
    os.makedirs(bowtie_output_dir, exist_ok=True)

                                                                        #open log file
    log_file = open(os.path.join(main, 'PipelineProject.log'), "a")         #append to the log file

                                                                    #change directory to input_reads_path
    os.chdir(input_reads_path)

                                                                    #iterate over fastq files in pairs
    for fastq1, fastq2 in zip(sorted(glob.glob(os.path.join(input_reads_path, "*_1.fastq"))),
                               sorted(glob.glob(os.path.join(input_reads_path, "*_2.fastq")))):
        match = re.search(r"SRR\d+", fastq1)                        #search for SRR pattern in the fastq file
        SRR_no = match.group()                                      #match the substring found above

                                                            #count the number of reads before bowtie filtering
        start_reads = popen(f"grep -c '@SRR' {fastq1}").read().strip()

        # Run Bowtie2
        index = os.path.join(main, 'bowtie_index', 'HCMV')
        bowtie_cmd = f"bowtie2 -x {index} -1 {fastq1} -2 {fastq2} -S {os.path.join(bowtie_output_dir, f'{SRR_no}_map.sam')} --al-conc {os.path.join(bowtie_output_dir, f'{SRR_no}_mapped_%.fq')}"
        os.system(bowtie_cmd)

                                                                    #count the reads after bowtie filtering
        end_reads = popen(f"grep -c '@SRR' {os.path.join(bowtie_output_dir, f'{SRR_no}_mapped_1.fq')}").read().strip()

                                                                    #write the read counts information to the log file
        log_file.write(f"{get_donor(SRR_no)} had {start_reads} read pairs before Bowtie2 filtering and {end_reads} read pairs after\n")

    log_file.write("\n\n\n\n")
    log_file.close()

def spades(main, log_file):                                         #spades function
    os.chdir(main)
    
                                                                    #makes a new directory for spades output
    os.makedirs(main + "/spades")
    
                                                #using glob to go through all files that end with fastq adnd add to the list
    b_list = []
    for name in glob.glob("**/*.fq", recursive=True):
        b_list.append(name)

                                            #sort the list of fastqs
    b_list = sorted(b_list)
    
                                                #paired-end spades assembly with only asssembler, outputs to spades/assembly
    spades_cmd = "spades.py " + "-k 77, 99, 127 -t 8 --only-assembler " + \
                 "--pe-1 1 " + b_list[0] + " --pe-2 1 " + b_list[1] + \
                 " --pe-1 2 " + b_list[2] + " --pe-2 2 " + b_list[3] + \
                 " --pe-1 3 " + b_list[4] + " --pe-2 3 " + b_list[5] + \
                 " --pe-1 4 " + b_list[6] + " --pe-2 4 " + b_list[7] + \
                 " -o " + "spades/" + "assembly/"
                                                                     #gets each pair from the list and runs through the spades command
                                                                     #writes the spades command to the log file
    log_file.write("# Spades Assembly\n")
    log_file.write(spades_cmd + "\n\n\n\n")

    os.system(spades_cmd)
    
def calc_contigs(main, log_file):
    os.chdir(main + "/spades/assembly/")
    
                                                                #parsing the contigs fasta from spades assembly
    records = SeqIO.parse(main + "/spades/assembly/contigs.fasta", format = "fasta")
    
                                            #find the longest contig and how many contigs greater than 1000 bp from the spades assembly
    num_contigs = 0
    contig_length = 0 
    longest_contig = ""
    for record in records:                  #if the length of the current seq is greater than the longest then that becomes the longest one
        if len(record.seq) > len(longest_contig):
            longest_contig = record
        if len(record.seq) > 1000:              #if teh length of the seq is greater than 1000 then add to the list 
            num_contigs += 1
            contig_length += len(record.seq)
            
                                                #write calc_contig info to the log file
    log_file.write("# Summary stats from SPAdes assembly output" + "\n")
    log_file.write("There are " + str(num_contigs) + " contigs > 1000 bp in the assembly." + "\n")
    log_file.write("There are " + str(contig_length) + " bp in the assembly." + "\n")
    log_file.write("\n" + "\n" + "\n" + "\n")
        
    return(longest_contig)                                  #return the longest contig

def blast(main, log_file):
                                                            #create blast directory
    blast_dir = os.path.join(main, 'blast')
    os.makedirs(blast_dir, exist_ok=True)

                                                        #change directory to blast directory
    os.chdir(blast_dir)

                                                        #create the HCMV database
    with Entrez.efetch(db="nucleotide", id='NC_006273.2', rettype="fasta") as handle:
        records = SeqIO.parse(handle, "fasta")
        SeqIO.write(records, "HCMV.fasta", "fasta")

                                                        #make BLAST database
    os.system("makeblastdb -in HCMV.fasta -dbtype nucl -out HCMV_db")

                                                        #get the longest contig
    longest_contig = calc_contigs(main, log_file)

                                                            #write longest contig to query.fasta
    with open("query.fasta", "w") as query_file:
        query_file.write(f">{longest_contig.id}\n{longest_contig.seq}\n")

                                                            #run blastn
    os.system("blastn -query query.fasta -db HCMV_db -num_threads 4 -max_hsps 1 -max_target_seqs 10 -out HCMV.txt -outfmt '6 sacc pident length qstart qend sstart send bitscore evalue stitle'")

                                                            #write the top 10 hit of BLAST output to log file
    log_file.write("# BLAST output\nsacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle\n")
    with open("HCMV.txt", "r") as blast_output:
        log_file.write(blast_output.read())
        
def run_pipeline(main, sample=False):                           #actually running the pipeline
    input_reads_path = os.path.join(os.getcwd(), "trimmed_fastq_files") if sample else os.path.join(main, "untrimmed_fastqs")

    log_file_path = os.path.join(main, "PipelineProject.log")
    with open(log_file_path, "a") as log_file:
        if not sample:
            sra_files(main)
            convert_sra_to_fastq(main)
        build_index(main)
        bowtie(main, input_reads_path)
        spades(main, log_file)
        blast(main, log_file)                                               #pass the log_file argument to blast function

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", help="run the pipeline on the trimmed fastq dataset", action="store_true")
    args = parser.parse_args()

    main_dir = os.path.join(os.getcwd(), "PipelineProject_HA")
    os.makedirs(main_dir, exist_ok=True)

    with open(os.path.join(main_dir, "PipelineProject.log"), "w") as log_file:
        log_file.write("## Pipeline Log ##\n\n")

    run_pipeline(main_dir, args.sample)