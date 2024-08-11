#!/bin/bash
#SBATCH --job-name=repeatM   # Job name
#SBATCH --mem=50gb                     # Job memory request, down from 200 and closer to the 64gb I think you're using per instance
#SBATCH --nodes=1                      # ensure cores are on one node
#SBATCH --ntasks=1                     # run a single task
#SBATCH --cpus-per-task=20              # number of cores/threads requested, up from 4 you're asking for now - see too that I changed --runThreadN below to match
#SBATCH --partition=20
#SBATCH --output=logs/repeatM_%j.log   # Standard output and error log
#SBATCH --error=logs/repeatM_%j.err

#RepeatMasker -lib <fasta file of repeats to mask>  <seqfiles(s) in fasta format>

accession=$1

cd ../../subdirs/$accession

common=$(cat "for_py2.py" | grep -o 'common="[^"]*' | cut -d'"' -f2)
echo $common
big_genome=$(grep -o 'big_genome="[^"]*' "config_unix_species.sh" | cut -d'"' -f2)
echo $big_genome
if [ "$big_genome" = "yes" ]; then
    genome_location="../../genomes/${accession}_KH.fasta"
else
    genome_location=$(grep -o 'genome_location="[^"]*' "for_py2.py" | cut -d'"' -f2)
fi
echo $genome_location

awk 'BEGIN{OFS="\t"} {
    n=split($1, a, "_");
    if (length(a[n]) > 5) {
        sub(/.*_/, "", $1);
    } else {
        $1 = shortname "_" $1;
    }
    gsub(/\.0/, "", $1);  # Remove all instances of .0
    print;
}' shortname=$common clusters_margins.bed | sed 's/\.0//g' > clusters_margins2.bed

awk 'BEGIN{OFS="\t"} {
    n=split($1, a, "_");
    if (length(a[n]) > 5) {
        sub(/.*_/, "", $1);
    } else {
        $1 = shortname "_" $1;
    }
    gsub(/\.0/, "", $1);  # Remove all instances of .0
    print;
}' shortname=$common singletons_margins.bed | sed 's/\.0//g' > singletons_margins2.bed

awk 'BEGIN{OFS="\t"} {
    n=split($1, a, "_");
    if (length(a[n]) > 5) {
        sub(/.*_/, "", $1);
    } else {
        $1 = shortname "_" $1;
    }
    gsub(/\.0/, "", $1);  # Remove all instances of .0
    print;
}' shortname=$common random_windows.bed | sed 's/\.0//g' > random_windows2.bed

# awk 'BEGIN{OFS="\t"} {n=split($1, a, "_"); if (length(a[n]) > 5) sub(/.*_/, "", $1); print}' singletons_margins.bed > singletons_margins2.bed

singleton="singletons_margins2.bed"
clusters="clusters_margins2.bed"
random_regions="random_windows2.bed"

library="../../results/coordinate_analysis/repeat/references/${common}-families.fa"

fasta_singleton="../../results/coordinate_analysis/repeat/mini_run/${common}_singleton.fasta"
fasta_cluster="../../results/coordinate_analysis/repeat/mini_run/${common}_cluster.fasta"
fasta_random="../../results/coordinate_analysis/repeat/mini_run/${common}_random.fasta"

repeat_singleton="../../results/coordinate_analysis/repeat/mini_run/${common}_singleton_repeat.fasta"
repeat_cluster="../../results/coordinate_analysis/repeat/mini_run/${common}_cluster_repeat.fasta"
repeat_cluster="../../results/coordinate_analysis/repeat/mini_run/${common}_random_repeat.fasta"

bedtools getfasta -fo $fasta_singleton -fi $genome_location -bed $singleton
bedtools getfasta -fo $fasta_cluster -fi $genome_location -bed $clusters
bedtools getfasta -fo $fasta_random -fi $genome_location -bed $random_regions

rm -f "$repeat_singleton"
rm -f "$repeat_cluster"
rm -f "$repeat_random"

RepeatMasker \
-lib $library \
-dir ../../results/coordinate_analysis/repeat/mini_run \
$fasta_singleton

RepeatMasker \
-lib $library \
-dir ../../results/coordinate_analysis/repeat/mini_run \
$fasta_cluster

RepeatMasker \
-lib $library \
-dir ../../results/coordinate_analysis/repeat/mini_run \
$fasta_random
