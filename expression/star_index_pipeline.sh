#!/bin/bash
#SBATCH --job-name=bull_index  # Job name
#SBATCH --mail-type=END,FAIL     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=youremailaddress@yourinstitute   # Where to send mail
#SBATCH --mem=200gb          # Job memory request, down from 200 and closer to the 64gb I think you're using per instance
#SBATCH --nodes=1           # ensure cores are on one node
#SBATCH --ntasks=1          # run a single task
#SBATCH --cpus-per-task=6       # number of cores/threads requested, up from 4 you're asking for now - see too that I changed --runThreadN below to match
#SBATCH --partition=20
#SBATCH --output=logs/bull_index_%j.log # Standard output and error log
#SBATCH --error=logs/bull_index_%j.err


STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir $1 \
--genomeFastaFiles $2 \
--sjdbGTFfile $3 \
--genomeChrBinNbits 12 \
--limitGenomeGenerateRAM 200000000000 \
--sjdbOverhang 99 \
--genomeSAsparseD 2
