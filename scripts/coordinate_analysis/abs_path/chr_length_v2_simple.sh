#!/bin/bash
#SBATCH --job-name=gen_size  # Job name
#SBATCH --mail-type=FAIL,END     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=youremailaddress@yourinstitute   # Where to send mail
#SBATCH --mem=5gb          # Job memory request, down from 200 and closer to the 64gb I think you're using per instance
#SBATCH --nodes=1           # ensure cores are on one node
#SBATCH --ntasks=1          # run a single task
#SBATCH --cpus-per-task=1       # number of cores/threads requested, up from 4 you're asking for now - see too that I changed --runThreadN below to match
#SBATCH --partition=20
#SBATCH --output=logs/gen_size_%j.log # Standard output and error log
#SBATCH --error=logs/gen_size_%j.err

while IFS= read -r accession; do
  # Specify the path to the directory you want to change to
  cd "/lab/wengpj01/vertebrate_pipeline/subdirs/$accession"
    # if [ -e "chromsizes_1007.csv" ]; then
    #     rm "chromsizes_1007.csv"
    #     fi
  # echo $accession
  # Use a for loop to iterate over matching files
  for genome_location in /lab/wengpj01/genomes/ncbi-genomes-2023-05-24/"$accession"*.fna; do
  #for genome_location in /lab/wengpj01/genomes/ncbi-genomes-2023-05-24/"$accession"_KH.fasta; do
    if [ -e "$genome_location" ]; then
      samtools faidx $genome_location
      # Process the genome file and add a column to chromsizes.csv
      faidx_output=$(faidx -i chromsizes "$genome_location")

      # Save the result to chromsizes.csv
      echo "$faidx_output" > chromsizes.csv

      echo "Successfully processed accession: $accession"
    else
      echo "Genome file not found for accession: $accession"
    fi
  done
done < "/lab/wengpj01/vertebrate_pipeline/accessions_to_keep.txt"
#done < "/lab/wengpj01/vertebrate_pipeline/big_genomes.txt"
