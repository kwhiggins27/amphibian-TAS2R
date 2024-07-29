#!/bin/bash
#SBATCH --job-name=blastF  # Job name
#SBATCH --mem=50gb          # Job memory request, down from 200 and closer to the 64gb I think you're using per instance
#SBATCH --nodes=1           # ensure cores are on one node
#SBATCH --ntasks=1          # run a single task
#SBATCH --cpus-per-task=10       # number of cores/threads requested, up from 4 you're asking for now - see too that I changed --runThreadN below to match
#SBATCH --partition=20
#SBATCH --output=blastF.log # Standard output and error log
#SBATCH --error=blastF.err

shopt -s extglob

cd $1

source config_unix_species.sh

echo $directory
echo $genome_location

cd $directory

query="../../scripts/blast_pipeline/bitter_reference_2023.fasta"
out=$directory"/forward_blast.tsv"

makeblastdb \
-in $genome_location \
-dbtype nucl

tblastn \
-query $query \
-db $genome_location \
-outfmt "6 qacc sacc evalue sstart send pident stitle sseq qstart qend" \
-num_threads 10 \
-out $out
