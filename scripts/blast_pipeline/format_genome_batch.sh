#!/bin/bash
#SBATCH --job-name=gen  # Job name
#SBATCH --mail-type=END,FAIL     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=youremailaddress@yourinstitute   # Where to send mail
#SBATCH --mem=4gb          # Job memory request, down from 200 and closer to the 64gb I think you\'re using per instance
#SBATCH --nodes=1           # ensure cores are on one node
#SBATCH --ntasks=1          # run a single task
#SBATCH --cpus-per-task=1       # number of cores/threads requested, up from 4 you\'re asking for now - see too that I changed --runThreadN below to match
#SBATCH --partition=20
#SBATCH --output=logs/gen_%j.log # Standard output and error log
#SBATCH --error=logs/gen_%j.err

sbatch format_genome.sh	GCA_000225785.1_LatCha1_genomic.fna	GCA_000225785.1  Latimeria_chalumnae  coelacanth 7897 fish
