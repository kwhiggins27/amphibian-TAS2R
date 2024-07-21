#!/bin/bash
#SBATCH --job-name=100K_KH  # Job name
#SBATCH --mail-type=NONE     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=youremailaddress@yourinstitute    # Where to send mail
#SBATCH --mem=50gb          # Job memory request, down from 200 and closer to the 64gb I think you're using per instance
#SBATCH --nodes=1           # ensure cores are on one node
#SBATCH --ntasks=1          # run a single task
#SBATCH --cpus-per-task=10       # number of cores/threads requested, up from 4 you're asking for now - see too that I changed --runThreadN below to match
#SBATCH --partition=20
#SBATCH --output=logs/100K_%j.log # Standard output and error log
#SBATCH --error=logs/100K_%j.err


#Use this script for finding nearst neighbor
max_skip=1000000

# Read the accession numbers from the file and process them one by one
while IFS= read -r accession; do
  # Run the find_clusters_gtf.py script for each accession
  echo "Processing accession: $accession"
  if ./cluster_pipeline_nearest_neighbor.sh "$accession" $max_skip; then
    echo "Successfully processed accession: $accession"
  else
    echo $accession >> /lab/wengpj01/vertebrate_pipeline/cluster_fail_pt3.txt
  fi
done < "/lab/wengpj01/vertebrate_pipeline/accessions_to_keep.txt"
