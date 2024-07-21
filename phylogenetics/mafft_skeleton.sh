#!/bin/bash
#SBATCH --job-name=mafft  # Job name
#SBATCH --mail-type=FAIL,END     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=youremailaddress@yourinstitute    # Where to send mail
#SBATCH --mem=200gb          # Job memory request, down from 200 and closer to the 64gb I think you're using per instance
#SBATCH --nodes=1           # ensure cores are on one node
#SBATCH --ntasks=1          # run a single task
#SBATCH --cpus-per-task=10       # number of cores/threads requested, up from 4 you're asking for now - see too that I changed --runThreadN below to match
#SBATCH --partition=20
#SBATCH --output=logs/mafft_%j.log # Standard output and error log
#SBATCH --error=logs/mafft_%j.err

#Used for adding seqs to skeleton tree below
#mafft --add new_sequences --reorder existing_alignment > output
mafft --add /lab/wengpj01/vertebrate_pipeline/all_pro_chromosomal_20231006.fasta \
--thread 10 \
--maxiterate 1000 \
--reorder /lab/wengpj01/vertebrate_pipeline/trees/skeleton_nopseu_reference_genafpair \
> /lab/wengpj01/vertebrate_pipeline/trees/skel_nopseu_reference_withadded_maxit1K_1006


# ##Used this for creating tree with all non-pseu in reference.
# mafft --maxiterate 1000 --genafpair \
# /lab/solexa_weng/playground/Kate_Higgins/other_vertebrates/bitter_reference_2023_nopseu.fasta \
# > /lab/wengpj01/vertebrate_pipeline/trees/skeleton_nopseu_reference_genafpair
#
#
# #Used this to attempt tree with all genes.
# mafft --maxiterate 1000 \
# --genafpair \
# --thread 20 \
# /lab/wengpj01/vertebrate_pipeline/all_pro_chromosomal.txt \
# > /lab/wengpj01/vertebrate_pipeline/trees/full_test_genafpair
