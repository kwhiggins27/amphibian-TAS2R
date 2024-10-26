#!/bin/bash
#SBATCH --job-name=blastR  # Job name
#SBATCH --mem=10gb          # Job memory request, down from 200 and closer to the 64gb I think you're using per instance
#SBATCH --nodes=1           # ensure cores are on one node
#SBATCH --ntasks=1          # run a single task
#SBATCH --cpus-per-task=10       # number of cores/threads requested, up from 4 you're asking for now - see too that I changed --runThreadN below to match
#SBATCH --partition=20
#SBATCH --output=blastR.log # Standard output and error log
#SBATCH --error=blastR.err

shopt -s extglob

cd $1

source config_unix_species.sh

query=$1"/for_blastR.fasta"
out=$1"/reciprocal_blast.tsv"

makeblastdb \
-in $reference_genome \
-dbtype nucl

tblastn \
-query $query \
-db $reference_genome \
-outfmt "6 qacc sacc evalue sstart send pident stitle sseq qstart qend" \
-num_threads 10 \
-out $out
