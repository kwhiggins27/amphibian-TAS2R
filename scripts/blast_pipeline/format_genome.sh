#!/bin/bash
#SBATCH --job-name=format  # Job name
#SBATCH --mail-type=NONE     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=youremailaddress@yourinstitute   # Where to send mail
#SBATCH --mem=50gb          # Job memory request, down from 200 and closer to the 64gb I think you're using per instance
#SBATCH --nodes=1           # ensure cores are on one node
#SBATCH --ntasks=1          # run a single task
#SBATCH --cpus-per-task=1       # number of cores/threads requested, up from 4 you're asking for now - see too that I changed --runThreadN below to match
#SBATCH --partition=20
#SBATCH --output=logs/format.log # Standard output and error log
#SBATCH --error=logs/format.err

#1: genome without folder
#2: accessation
#3: latin
#4: common
#5: taxa
#6: reference

genome="../../genomes/$1"

shortname=$2

common=$4

echo $shortname

taxa=$6

if [ "$taxa" = "mammal" ]; then
  reference_sp="human"
  reference_genome="../../references/reference_genomes/hg38.fa"
  coordinates="../../references/reference_coordinates/human_tas2r_reference.csv"
elif [ "$taxa" = "amphibian" ]; then
  reference_sp="xenopus"
  reference_genome="../../references/reference_genomes/Xenopus_Tropicalis_2022.fasta"
  coordinates="../../references/reference_coordinates/xenopus_TAS2R_reference_20230412.csv"
elif [ "$taxa" = "reptile" ]; then
  reference_sp="anole"
  reference_genome="../../references/reference_genomes/anole_genome.fasta"
  coordinates="../../references/reference_coordinates/anole_tas2r_reference.csv"
elif [ "$taxa" = "bird" ]; then
  reference_sp="chicken"
  reference_genome="../../references/reference_genomes/chicken_genome.fasta"
  coordinates="../../references/reference_coordinates/chicken_tas2r_reference.csv"
elif [ "$taxa" = "fish" ]; then
  reference_sp="coelacanth"
  reference_genome="../../references/reference_genomes/coelacanth_genome.fasta"
  coordinates="../../references/reference_coordinates/coelacanth_tas2r_reference.csv"
else
  echo "Unknown taxa: $species"
  exit 1
fi

echo $reference_sp

#underscore=$(awk '/^>/ { sub(">", ""); split($1, arr, " "); if (index(arr[1], "_") > 0) { print "yes"; exit } } END { print "no" }' "$genome_unzipped")
underscore=$(awk 'NR==1 && /^>.* / { if (index($1, "_") > 0) { print "yes" } else { print "no" }; exit }' "$genome")


echo $underscore

now=$(date +%Y%m%d%H)

config_full="../../configs/config_full_${shortname}.txt"

mkdir -p "../../subdirs/${shortname}"

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
echo "directory=\"../../subdirs/${shortname}\"" >> $config_full

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
echo "folder=\"../../subdirs/${shortname}\"" >> $config_full
echo "logfile=\"../../logs/summaryruns.csv\"" >> $config_full
