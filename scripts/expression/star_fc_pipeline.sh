#!/bin/bash
#SBATCH --job-name=x_combo   # Job name
#SBATCH --mem=100gb                     # Job memory request, down from 200 and closer to the 64gb I think you're using per instance
#SBATCH --nodes=1                      # ensure cores are on one node
#SBATCH --ntasks=1                     # run a single task
#SBATCH --cpus-per-task=10              # number of cores/threads requested, up from 4 you're asking for now - see too that I changed --runThreadN below to match
#SBATCH --partition=20
#SBATCH --output=logs/x_combo_%j.log   # Standard output and error log
#SBATCH --error=logs/x_combo_%j.err

#"1" F_READS
#2 R_READS
#3 star index
#4 star out --> "$STAR_out"
#5 "$GTF"
#6 featureoutput file --> "$FC_out"
#7 star out with variable and ending
#8 var

#"$F_READS" "$R_READS" "$STAR_index" "$STAR_out""$var" "$GTF" "$FC_out" "$temp" "$var"


mkdir -p "$4"
mkdir -p "$6"

ulimit -n 20000

#Run STAR
STAR \
--runMode alignReads \
--genomeDir "$3" \
--readFilesIn "$1" "$2" \
--readFilesCommand zcat \
--runThreadN 10 \
--outFileNamePrefix "$4" \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--quantMode TranscriptomeSAM \
--genomeSAsparseD 2

echo "Done with STAR"
date +%Y%m%d%H%M%S

echo "$6/""$8""_fc.tsv"

#featureCounts
##Added -M on 5/2/23
featureCounts \
-a "$5" \
-o "$6/""$8""_fc.tsv" \
"$7"


sleep 90
