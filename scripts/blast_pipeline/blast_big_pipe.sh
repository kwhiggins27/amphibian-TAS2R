#!/bin/bash
#SBATCH --job-name=blastBig  # Job name
#SBATCH --mem=10gb          # Job memory request, down from 200 and closer to the 64gb I think you're using per instance
#SBATCH --nodes=1           # ensure cores are on one node
#SBATCH --ntasks=1          # run a single task
#SBATCH --cpus-per-task=10       # number of cores/threads requested, up from 4 you're asking for now - see too that I changed --runThreadN below to match
#SBATCH --partition=20
#SBATCH --output=blastBig.log # Standard output and error log
#SBATCH --error=blastBig.err

shopt -s extglob

cd $1

source config_unix_species.sh
accession=$3

look_for_genome=$(dirname $genome_location)
subdivided="${look_for_genome}/${accession}_subdivided"

fasta=$subdivided"/$2.fasta"

cd $subdivided

query="../../scripts/blast_pipeline/bitter_reference_2023.fasta"
makeblastdb \
-in $fasta \
-dbtype nucl


tblastn \
-query $query  \
-db $fasta \
-outfmt "6 qacc sacc evalue sstart send pident stitle sseq qstart qend" \
-num_threads 10 \
-out $directory/blastF/"$2"_local_BLAST.tsv
