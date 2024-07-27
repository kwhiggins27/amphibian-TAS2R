#!/bin/bash
#SBATCH --job-name=combotrim  # Job name
#SBATCH --mail-type=END,FAIL     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=youremailaddress@yourinstitute    # Where to send mail
#SBATCH --mem=400gb          # Job memory request, down from 200 and closer to the 64gb I think you're using per instance
#SBATCH --nodes=1           # ensure cores are on one node
#SBATCH --ntasks=1          # run a single task
#SBATCH --cpus-per-task=5       # number of cores/threads requested, up from 4 you're asking for now - see too that I changed --runThreadN below to match
#SBATCH --partition=20
#SBATCH --output=combotrim.log # Standard output and error log
#SBATCH --error=combotrim.err

path_in="/lab/solexa_public/Weng/230404_WIGTC-NOVASEQ1A_AHWMGTDSX5/FASTQ/L319_"
forward="_R1_001.fastq.gz"
reverse="_R2_001.fastq.gz"

#echo "$path_in""$1""$forward"

java \
-jar /lab/solexa_weng/testtube/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
-threads 5 \
-phred33 "$path_in""$1""$forward" "$path_in""$1""$reverse" \
"$1"_trim_b_R1.fastq.gz "$1"_trim_b_R1U.fastq.gz \
"$1"_trim_b_R2.fastq.gz "$1"_trim_b_R2U.fastq.gz \
ILLUMINACLIP:/archive/weng/2022.05.23-30905/solexa_weng/Seq_data/Projects/Kate_Higgins2/trimmed_reads/Roberto_adapters.fa:3:30:10 LEADING:15 TRAILING:15 SLIDINGWINDOW:4:20 MINLEN:36


# java \
# -jar /lab/solexa_weng/testtube/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
# -threads 3 \
# -phred33 "$path_in""$1""$forward" "$path_in""$1""$reverse" \
# "$1"_trim_R1.fastq.gz "$1"_trim_R1U.fastq.gz \
# "$1"_trim_R2.fastq.gz "$1"_trim_R2U.fastq.gz \
# ILLUMINACLIP:/archive/weng/2022.05.23-30905/solexa_weng/Seq_data/Projects/Kate_Higgins2/trimmed_reads/Roberto_adapters.fa:3:30:10 LEADING:15 TRAILING:15 SLIDINGWINDOW:4:20 MINLEN:36
