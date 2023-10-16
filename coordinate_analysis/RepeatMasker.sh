#!/bin/bash
#SBATCH --job-name=repeatM   # Job name
#SBATCH --mail-type=NONE          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kwh1@wi.mit.edu     # Where to send mail
#SBATCH --mem=50gb                     # Job memory request, down from 200 and closer to the 64gb I think you're using per instance
#SBATCH --nodes=1                      # ensure cores are on one node
#SBATCH --ntasks=1                     # run a single task
#SBATCH --cpus-per-task=20              # number of cores/threads requested, up from 4 you're asking for now - see too that I changed --runThreadN below to match
#SBATCH --partition=20
#SBATCH --output=logs/repeatM_%j.log   # Standard output and error log
#SBATCH --error=logs/repeatM_%j.err

#RepeatMasker -lib <fasta file of repeats to mask>  <seqfiles(s) in fasta format>

accession=$1

cd /lab/wengpj01/vertebrate_pipeline/subdirs/$accession

common=$(cat "/lab/wengpj01/vertebrate_pipeline/subdirs/$accession/for_py2.py" | grep -o 'common="[^"]*' | cut -d'"' -f2)
echo $common
big_genome=$(grep -o 'big_genome="[^"]*' "/lab/wengpj01/vertebrate_pipeline/subdirs/$accession/config_unix_species.sh" | cut -d'"' -f2)
if [ "$big_genome" = "yes" ]; then
    genome_location="/lab/wengpj01/genomes/ncbi-genomes-2023-05-24/${accession}_KH.fasta"
else
    genome_location=$(grep -o 'genome_location="/[^"]*' "/lab/wengpj01/vertebrate_pipeline/subdirs/$accession/for_py2.py" | cut -d'"' -f2)
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
singleton=singletons_margins2.bed
clusters=clusters_margins2.bed
random_regions=random_windows2.bed

library="/lab/wengpj01/repeat/references/${common}-families.fa"

fasta_singleton="/lab/wengpj01/repeat/amphibian_near_TAS2R/${common}_singleton.fasta"
fasta_cluster="/lab/wengpj01/repeat/amphibian_near_TAS2R/${common}_cluster.fasta"
fasta_random="/lab/wengpj01/repeat/amphibian_near_TAS2R/${common}_random.fasta"

repeat_singleton="/lab/wengpj01/repeat/amphibian_near_TAS2R/${common}_singleton_repeat.fasta"
repeat_cluster="/lab/wengpj01/repeat/amphibian_near_TAS2R/${common}_cluster_repeat.fasta"
repeat_cluster="/lab/wengpj01/repeat/amphibian_near_TAS2R/${common}_random_repeat.fasta"

bedtools getfasta -fo $fasta_singleton -fi $genome_location -bed $singleton
bedtools getfasta -fo $fasta_cluster -fi $genome_location -bed $clusters
bedtools getfasta -fo $fasta_random -fi $genome_location -bed $random_regions

cd /lab/wengpj01/repeat

rm -f "$repeat_singleton"
rm -f "$repeat_cluster"
rm -f "$repeat_random"

RepeatMasker \
-lib $library \
-dir /lab/wengpj01/repeat/amphibian_near_TAS2R \
$fasta_singleton

RepeatMasker \
-lib $library \
-dir /lab/wengpj01/repeat/amphibian_near_TAS2R \
$fasta_cluster

RepeatMasker \
-lib $library \
-dir /lab/wengpj01/repeat/amphibian_near_TAS2R \
$fasta_random

# RepeatMasker \
# -lib /lab/wengpj01/repeat/axolotl-families.fa \
# -dir custom_library_near_TAS2R \
# /lab/wengpj01/repeat/near_TAS2R/axolotl_tas2r.fasta

# RepeatMasker \
# -lib /lab/wengpj01/repeat/bullfrog-families.fa \
# -dir custom_library_near_TAS2R \
# /lab/wengpj01/repeat/near_TAS2R/bullfrog_tas2r.fasta

# RepeatMasker \
# -lib /lab/wengpj01/repeat/cane-families.fa \
# -dir custom_library_near_TAS2R \
# /lab/wengpj01/repeat/near_TAS2R/cane_tas2r.fasta
#
# RepeatMasker \
# -lib /lab/wengpj01/repeat/xenopus-families.fa \
# -dir custom_library_near_TAS2R \
# /lab/wengpj01/repeat/near_TAS2R/xenopus_tas2r.fasta
#
# RepeatMasker \
# -lib /lab/wengpj01/repeat/terribilis-families.fa \
# -dir custom_library_near_TAS2R \
# /lab/wengpj01/repeat/near_TAS2R/terribilis_tas2r.fasta

# RepeatMasker \
# -lib /lab/wengpj01/repeat/axolotl-families.fa \
# -dir custom_library_sampled_control \
# /lab/wengpj01/repeat/sampled_gtf/axolotl_sampled.fasta
#
# RepeatMasker \
# -lib /lab/wengpj01/repeat/bullfrog-families.fa \
# -dir custom_library_sampled_control \
# /lab/wengpj01/repeat/sampled_gtf/bullfrog_sampled.fasta
#
# RepeatMasker \
# -lib /lab/wengpj01/repeat/cane-families.fa \
# -dir custom_library_sampled_control \
# /lab/wengpj01/repeat/sampled_gtf/cane_sampled.fasta

# RepeatMasker \
# -lib /lab/wengpj01/repeat/xenopus-families.fa \
# -dir custom_library_sampled_control \
# /lab/wengpj01/repeat/sampled_gtf/xenopus_sampled.fasta


#RepeatMasker -dir default_repeat_genome /lab/solexa_weng/playground/Kate_Higgins/other_vertebrates/terribilis/P.terribilis.gapclosed.fasta

#RepeatMasker -dir default_repeat  /lab/wengpj01/repeat/xenopus_near_tas2r.fasta
