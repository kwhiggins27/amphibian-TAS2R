#!/bin/bash
#SBATCH --job-name=gen  # Job name
#SBATCH --mem=4gb          # Job memory request, down from 200 and closer to the 64gb I think you\'re using per instance
#SBATCH --nodes=1           # ensure cores are on one node
#SBATCH --ntasks=1          # run a single task
#SBATCH --cpus-per-task=1       # number of cores/threads requested, up from 4 you\'re asking for now - see too that I changed --runThreadN below to match
#SBATCH --partition=20
#SBATCH --output=logs/gen_%j.log # Standard output and error log
#SBATCH --error=logs/gen_%j.err

#sbatch format_genome.sh	GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna	GCA_016699485.1 Gallus_gallus  chicken 9031 bird
sbatch format_genome.sh	GCA_009819775.1_bPhoRub2.pri_genomic.fna	GCA_009819775.1 Gallus_gallus American_flamingo 9031 bird
sbatch format_genome.sh GCA_027887145.1_bAmmCau1.pri_genomic.fna GCA_027887145.1 Gallus_gallus birds 9031 bird
