#!/bin/bash
#SBATCH --job-name=format  # Job name
#SBATCH --mail-type=NONE     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=youremailaddress@yourinstitute   # Where to send mail
#SBATCH --mem=50gb          # Job memory request, down from 200 and closer to the 64gb I think you're using per instance
#SBATCH --nodes=1           # ensure cores are on one node
#SBATCH --ntasks=1          # run a single task
#SBATCH --cpus-per-task=1       # number of cores/threads requested, up from 4 you're asking for now - see too that I changed --runThreadN below to match
#SBATCH --partition=20
#SBATCH --output=format.log # Standard output and error log
#SBATCH --error=format.err

#1: genome without folder
#2: accessation
#3: latin
#4: common
#5: taxa
#6: reference

genome="/lab/wengpj01/genomes/ncbi-genomes-2023-06-07/"$1


shortname=$2

common=$4

echo $shortname

taxa=$6

if [ "$taxa" = "mammal" ]; then
  reference_sp="human"
  reference_genome="/archive/weng/2022.09.13-31823/lab/solexa_weng/playground/Kate_Higgins/human/hg38.fa"
  coordinates="/lab/solexa_weng/playground/Kate_Higgins/other_vertebrates/human/human_tas2r_reference.csv"
elif [ "$taxa" = "amphibian" ]; then
  reference_sp="xenopus"
  reference_genome="/lab/wengpj01/xenopus/Xenopus_Tropicalis_2022.fasta"
  coordinates="/lab/solexa_weng/playground/Kate_Higgins/other_vertebrates/xenopus/xenopus_TAS2R_reference_20230412.csv"
elif [ "$taxa" = "reptile" ]; then
  reference_sp="anole"
  reference_genome="/lab/solexa_weng/playground/Kate_Higgins/other_vertebrates/anole/anole_genome.fasta"
  coordinates="/lab/solexa_weng/playground/Kate_Higgins/other_vertebrates/anole/anole_tas2r_reference.csv"
elif [ "$taxa" = "bird" ]; then
  reference_sp="chicken"
  reference_genome="/lab/solexa_weng/playground/Kate_Higgins/other_vertebrates/chicken/chicken_genome.fasta"
  coordinates="/lab/solexa_weng/playground/Kate_Higgins/other_vertebrates/chicken/chicken_tas2r_reference.csv"
elif [ "$taxa" = "fish" ]; then
  reference_sp="coelacanth"
  reference_genome="/lab/solexa_weng/playground/Kate_Higgins/other_vertebrates/coelacanth/coelacanth_genome.fasta"
  coordinates="/lab/solexa_weng/playground/Kate_Higgins/other_vertebrates/coelacanth/coelacanth_tas2r_reference.csv"
else
  echo "Unknown taxa: $species"
  exit 1
fi

echo $reference_sp

#underscore=$(awk '/^>/ { sub(">", ""); split($1, arr, " "); if (index(arr[1], "_") > 0) { print "yes"; exit } } END { print "no" }' "$genome_unzipped")
underscore=$(awk 'NR==1 && /^>.* / { if (index($1, "_") > 0) { print "yes" } else { print "no" }; exit }' "$genome")


echo $underscore

now=$(date +%Y%m%d%H)

config_full="/lab/wengpj01/vertebrate_pipeline/configs/config_full_${shortname}.txt"
#for_py="/lab/wengpj01/vertebrate_pipeline/configs/for_py2_${shortname}.py"

mkdir -p "/lab/wengpj01/vertebrate_pipeline/subdirs/${shortname}"

echo $config_full

touch $config_full
#touch $for_py


echo "BLAST_species=\"$common\"" > $config_full
echo "reference_sp=\"$reference_sp\"" >> $config_full
echo "num_threads_forward=5" >> $config_full
echo "num_threads_reverse=5" >> $config_full
echo "tasted=\"bitter\"" >> $config_full
echo "run=\"v1\"" >> $config_full
echo "type=\"genome\"" >> $config_full
echo "genome_location=\"$genome\"" >> $config_full
echo "reference_genome=\"$reference_genome\"" >> $config_full
echo "big_genome=\"no\"" >> $config_full
echo "min_aa=200" >> $config_full
echo "max_aa=500" >> $config_full
echo "go_to=0" >> $config_full
echo "existing_directory=\"\"" >> $config_full

echo "#do not add anything after this point" >> $config_full
echo "directory=\"/lab/wengpj01/vertebrate_pipeline/subdirs/${shortname}\"" >> $config_full

echo "##Variables for python" >> $config_full
echo "BLAST_species=\"$common\"" >> $config_full
echo "reference_species=\"$reference_sp\"" >> $config_full
echo "min_aa=200" >> $config_full
echo "max_aa=500" >> $config_full
echo "tasted=\"bitter\"" >> $config_full
echo "run=\"v1\"" >> $config_full
echo "type=\"genome\"" >> $config_full
echo "chrom_name_have_underscore=\"$underscore\"" >> $config_full
echo "genome_location=\"$genome\"" >> $config_full
echo "reference_genome=\"reference_genome\"" >> $config_full
echo "borders=\"new\"">> $config_full
echo "full_pull=BLAST_species + \"_full_pull.sh\"">> $config_full
echo "coordinates=\"$coordinates\"">> $config_full
echo "accession=\"$2\"">> $config_full
echo "latin=\"$3\"">> $config_full
echo "common=\"$4\"">> $config_full
echo "taxa=\"$5\"">> $config_full
echo "folder=\"/lab/wengpj01/vertebrate_pipeline/subdirs/${shortname}\"" >> $config_full
echo "logfile=\"/lab/wengpj01/vertebrate_pipeline/20230604_summaryruns.csv\"" >> $config_full

# #!/usr/bin/env python
# BLAST_species="ornatum"
# reference_species="xenopus"
# min_aa=200
# max_aa=500
# tasted="bitter"
# run="v1"
# type="genome"
# chrom_name_have_underscore="yes"
# genome_location="/archive/weng/2022.09.13-31823/lab/solexa_weng/playground/Kate_Higgins/other_vertebrates/ornatum/Platyplectrum_ornatum_frog.fna"
# reference_genome="/lab/wengpj01/xenopus/Xenopus_Tropicalis_2022.fasta"
# borders="new" #wide or narrow, where narrow shown to work better in terribilis
# #BLAST_file="/lab/solexa_weng/playground/Kate_Higgins/other_vertebrates/xenopus/xenopus_bitter_230403_blast.tsv"
# #expanded_file=BLAST_species + "_expanded2023_seq.fasta"
# #expanded_file="/lab/solexa_weng/playground/Kate_Higgins/other_vertebrates/xenopus/xenopus_expanded230403_fortest.fasta"
# full_pull=BLAST_species + "_full_pull.sh"
#
# #blast_reciprocal="/lab/solexa_weng/playground/Kate_Higgins/other_vertebrates/xenopus/xenopus_reciprocal230403_Thurs_BLAST.tsv"
#
# coordinates="/lab/solexa_weng/playground/Kate_Higgins/other_vertebrates/xenopus/xenopus_TAS2R_reference_20230412.csv"
# folder="/lab/solexa_weng/playground/Kate_Higgins/other_vertebrates/ornatum/20230520195012/"


# ##for unix
# BLAST_species="lungfish"
# reference_sp="coelacanth"
# now=$(date +%Y%m%d%H)
# num_threads_forward=5
# num_threads_reverse=5
# tasted="bitter"
# run="v1"
# type="genome"
# genome_location="/lab/wengpj01/non-amphibians/lungfish/subdivided"
# reference_genome="/lab/solexa_weng/playground/Kate_Higgins/other_vertebrates/coelacanth/coelacanth_genome.fasta"
# big_genome="yes"
# min_aa=200
# max_aa=500
#
# go_to=0 #0 for beginning with 1st blast, 1 for py1, 2 for fullpull, 3 for py2, 4 for blastR, 5 for py3
# existing_directory="/lab/solexa_weng/playground/Kate_Higgins/other_vertebrates/terribilis/20230514183109"
#
#
# #do not add anything after this point
# directory="/lab/wengpj01/non-amphibians/lungfish/20230525100749"/
