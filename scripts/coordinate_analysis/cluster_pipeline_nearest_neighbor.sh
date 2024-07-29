#!/bin/bash
#SBATCH --job-name=clusters  # Job name
#SBATCH --mail-type=FAIL,END     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=youremailaddress@yourinstitute     # Where to send mail
#SBATCH --mem=200gb          # Job memory request, down from 200 and closer to the 64gb I think you're using per instance
#SBATCH --nodes=1           # ensure cores are on one node
#SBATCH --ntasks=1          # run a single task
#SBATCH --cpus-per-task=20       # number of cores/threads requested, up from 4 you're asking for now - see too that I changed --runThreadN below to match
#SBATCH --partition=20
#SBATCH --output=logs/clusters_%j.log # Standard output and error log
#SBATCH --error=logs/clusters_%j.err

accession=$1

max_skip=$2

cd ../../subdirs/$accession

#Need to run conda activate /lab/wengpj01/phyl_conda before starting

bitter_list=*_bitter_protein_confirmed2xv1.fasta

if [ ! -f mafft.fas ]; then
  mafft --maxiterate 1000 \
  --genafpair \
  --thread 20 \
  $bitter_list \
  > mafft.fas
fi

if [ ! -f distmat_1.dist ]; then
  distmat \
  -sequence mafft.fas \
  -protmethod 1 \
  distmat_1.dist
fi

if [ ! -f distmat_2.dist ]; then
  awk -F'\t' 'NR == 7 {
      # Process the header row separately
      gsub(/^[[:space:]]+|[[:space:]]+$/, "");  # Remove leading and trailing spaces
      for (i = 1; i <= NF; i++) {
          printf("%s%s", $i, (i == NF) ? "\n" : "\t");
      }
  }
  NR > 7 {
      # Process data rows, excluding the first 7 lines
      for (i = 2; i <= NF - 2; i++) {
          gsub(/^[[:space:]]+|[[:space:]]+$/, "", $i);  # Remove leading and trailing spaces
          printf("%s%s", $i, (i == NF - 2) ? "\n" : "\t");
      }

  }'  distmat_1.dist > distmat_2.dist
fi

gtf=*.gtf

clusters=clusters_found_1005_$max_skip.csv

touch ../../results/coordinate_analysis/summary_clusters_withnearest_$max_skip.csv

# # Read in the distance matrices as pandas DataFrames
# evo = pd.read_csv('/lab/wengpj01/vertebrate_pipeline/trees/mini_wood_frog.dist', sep='\t')

python_script=../../scripts/coordinate_analysis/cluster_pipeline_nearest_neighbor.py
$python_script $accession $gtf $clusters distmat_2.dist $max_skip
