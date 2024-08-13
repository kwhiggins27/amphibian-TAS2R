#!/bin/bash
#SBATCH --job-name=busco  # Job name
#SBATCH --mem=50gb          # Job memory request, down from 200 and closer to the 64gb I think you're using per instance
#SBATCH --nodes=1           # ensure cores are on one node
#SBATCH --ntasks=1          # run a single task
#SBATCH --cpus-per-task=10       # number of cores/threads requested, up from 4 you're asking for now - see too that I changed --runThreadN below to match
#SBATCH --partition=20
#SBATCH --output=logs/busco_%j.log # Standard output and error log
#SBATCH --error=logs/busco_%j.err

remaining="../../results/busco/busco_needed.txt"
 #"/lab/wengpj01/vertebrate_pipeline/all_accessions_remaining.txt"
complete="../../../../results/busco/list_complete.txt"
attempted="../../../../results/busco/list_attempted.txt"
failed="../../../../results/busco/list_failed.txt"

accession="$1"
sed -i "/$accession/d" "$remaining"

touch $attempted
echo $accession >> $attempted
#
# set -e
set -o pipefail

export BUSCO_CONFIG_FILE="../../../../scripts/phylogenetics/config.ini"
export AUGUSTUS_CONFIG_PATH="../../../../scripts/phylogenetics/config/"

mkdir -p ../../results/busco/busco_results/$accession

cd ../../results/busco/busco_results/$accession

busco \
-m genome \
-i ../../../../genomes/$1*.fna \
-o busco_output \
-f \
-c 10 \
-l vertebrata_odb10 \
--augustus

sed -i "/$accession/d" "$attempted"
echo $accession >> $complete
