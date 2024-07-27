#!/bin/bash
#SBATCH --job-name=iqtree  # Job name
#SBATCH --mail-type=FAIL,END     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=youremailaddress@yourinstitute     # Where to send mail
#SBATCH --mem=200gb          # Job memory request, down from 200 and closer to the 64gb I think you're using per instance
#SBATCH --nodes=1           # ensure cores are on one node
#SBATCH --ntasks=1          # run a single task
#SBATCH --cpus-per-task=10       # number of cores/threads requested, up from 4 you're asking for now - see too that I changed --runThreadN below to match
#SBATCH --partition=20
#SBATCH --output=logs/iqtree_%j.log # Standard output and error log
#SBATCH --error=logs/iqtree_%j.err

# iqtree -s /lab/wengpj01/vertebrate_pipeline/trees/skel_nopseu_reference_withadded_nocolon.fas -m MFP -T AUTO --abayes

iqtree -s /lab/wengpj01/vertebrate_pipeline/trees/to_add_1006.fas -m JTT+F+R8 --abayes -T AUTO

# iqtree -s /lab/wengpj01/vertebrate_pipeline/trees/skeleton_nopseu_reference_genafpair.fas -m MFP -T AUTO --abayes
