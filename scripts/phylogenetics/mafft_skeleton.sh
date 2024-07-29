#!/bin/bash
#SBATCH --job-name=mafft  # Job name
#SBATCH --mem=200gb          # Job memory request, down from 200 and closer to the 64gb I think you're using per instance
#SBATCH --nodes=1           # ensure cores are on one node
#SBATCH --ntasks=1          # run a single task
#SBATCH --cpus-per-task=10       # number of cores/threads requested, up from 4 you're asking for now - see too that I changed --runThreadN below to match
#SBATCH --partition=20
#SBATCH --output=logs/mafft_%j.log # Standard output and error log
#SBATCH --error=logs/mafft_%j.err

# Make sure you've created and activated the conda environment described by requirements.txt

#Used for adding seqs to skeleton tree below
#mafft --add new_sequences --reorder existing_alignment > output
mafft --add ../../results/pipeline/all_pro_chromosomal_20231006.fasta \
--thread 10 \
--maxiterate 1000 \
--reorder ../../results/phylogenetics/skeleton_nopseu_reference_genafpair \
> ../../results/phylogenetics/skel_nopseu_reference_withadded_maxit1K_1006


# ##Used this for creating tree with all non-pseu in reference. Sequences from Behrens et al, 2021 and Li and Zhang, 2014 as described in text.
# mafft --maxiterate 1000 --genafpair \
# bitter_reference_2023_nopseu.fasta \
# > ../../results/phylogenetics/skeleton_nopseu_reference_genafpair
#
#
