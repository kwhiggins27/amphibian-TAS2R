#!/bin/bash
#SBATCH --job-name=exp_pipe  # Job name
#SBATCH --mail-type=END,FAIL     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=youremailaddress@yourinstitute  # Where to send mail
#SBATCH --mem=4gb          # Job memory request, down from 200 and closer to the 64gb I think you're using per instance
#SBATCH --nodes=1           # ensure cores are on one node
#SBATCH --ntasks=1          # run a single task
#SBATCH --cpus-per-task=1       # number of cores/threads requested, up from 4 you're asking for now - see too that I changed --runThreadN below to match
#SBATCH --partition=20
#SBATCH --output=logs/exp_pipe_%j.log # Standard output and error log
#SBATCH --error=logs/exp_pipe_%j.err

## Define variables
# species="axolotl"
# genome="/lab/wengpj01/axolotl/axolotl_mygenome.fasta"
# GTF="/lab/wengpj01/expression_pipeline/results_20231004/axolotl_from_pipeline.gtf"
# Bitter_TS="/lab/wengpj01/vertebrate_pipeline/subdirs/GCA_002915635.3/final_genes_7TM.csv"
# summary_table="/lab/wengpj01/expression_pipeline/results_20231004/axolotl_old.csv"
# go_to=0 #0 starts from beginning (star index), 1 starts with star/fc, 2 starts with python script
# jobname="ax_stfc"
#
# STAR_index="/lab/wengpj01/star/axolotl/axolotl_star_index_20231004/"
# STAR_out="/lab/wengpj01/axolotl/star_expression_20231004/"
# FC_out="/lab/wengpj01/axolotl/feature_counts_20231004/"

# # #
species="bullfrog"
genome="/lab/wengpj01/bullfrog/Lithobates_catesbeianus_bullfrog_2022.fasta"
GTF="/lab/wengpj01/expression_pipeline/results_20231004/bullfrog_from_pipeline.gtf"
Bitter_TS="/lab/solexa_weng/playground/Kate_Higgins/other_vertebrates/bullfrog/20231004072453/final_genes_7TM.csv"
summary_table="/lab/wengpj01/expression_pipeline/results_20231004/bullfrog_old.csv"
go_to=0 #0 starts from beginning (star index), 1 starts with star/fc, 2 starts with python script
jobname="b_stfc"

STAR_index="/lab/wengpj01/star/bullfrog/bullfrog_star_index_20231004/"
STAR_out="/lab/wengpj01/bullfrog/star_expression_20231004/"
FC_out="/lab/wengpj01/bullfrog/feature_counts_20231004/"

# # #

# species="cane"
# genome="/lab/wengpj01/cane/canetoad.v2.2.fasta"
# GTF="/lab/wengpj01/expression_pipeline/results_20231004/cane_from_pipeline.gtf"
# Bitter_TS="/lab/solexa_weng/playground/Kate_Higgins/other_vertebrates/cane/20231003083531/final_genes_7TM.csv"
# summary_table="/lab/wengpj01/expression_pipeline/results_20231004/cane_old.csv"
# go_to=1 #0 starts from beginning (star index), 1 starts with star/fc, 2 starts with python script
# jobname="c_stfc"
#
# STAR_index="/lab/wengpj01/star/cane/cane_star_index_20231004/"
# STAR_out="/lab/wengpj01/cane/star_expression_20231004/"
# FC_out="/lab/wengpj01/cane/feature_counts_20231004/"
#
# # #
#
# species="terribilis"
# genome="/lab/solexa_weng/playground/Kate_Higgins/other_vertebrates/terribilis/P.terribilis.gapclosed.fasta"
# GTF="/lab/wengpj01/expression_pipeline/results_20231004/terribilis_from_pipeline.gtf"
# Bitter_TS="/lab/solexa_weng/playground/Kate_Higgins/other_vertebrates/terribilis/20231003083547/final_genes_7TM.csv"
# summary_table="/lab/wengpj01/expression_pipeline/results_20231004/terribilis_old.csv"
# go_to=1 #0 starts from beginning (star index), 2 starts with star/fc, 3 starts with python script
# jobname="t_stfc"
#
# STAR_index="/lab/wengpj01/star/terribilis/terribilis_star_index_20231004/"
# STAR_out="/lab/wengpj01/terribilis/star_expression_20231004/"
# FC_out="/lab/wengpj01/terribilis/feature_counts_20231004/"
#
# # #
#
# species="xenopus"
# genome="/lab/wengpj01/xenopus/Xenopus_Tropicalis_2022.fasta"
# GTF="/lab/wengpj01/expression_pipeline/results_20231004/xenopus_from_pipeline.gtf"
# Bitter_TS="/lab/solexa_weng/playground/Kate_Higgins/other_vertebrates/xenopus/20231003083557/final_genes_7TM.csv"
# summary_table="/lab/wengpj01/expression_pipeline/results_20231004/xenopus_old.csv"
# go_to=1 #0 starts from beginning (star index), 1 starts with star/fc, 2 starts with python script
# jobname="x_stfc"
#
# STAR_index="/lab/wengpj01/star/xenopus/xenopus_star_index_20231004/"
# STAR_out="/lab/wengpj01/xenopus/star_expression_20231004/"
# FC_out="/lab/wengpj01/xenopus/feature_counts_20231004/"
#
# #
path_in="/lab/wengpj01/trimmed/"
ending="Aligned.sortedByCoord.out.bam"
forward="_trim_R1.fastq.gz"
reverse="_trim_R2.fastq.gz"


##Remember that some axolotl trimmed files have a different name!

#Create directories, if they don't already exist
mkdir -p $STAR_index
mkdir -p $STAR_out
mkdir -p $FC_out

#Create star index
if [ $go_to -gt 0 ]
then
  echo "STAR index already made"
else
  checkme="$STAR_index""SAindex"
  echo $checkme
  if [ -f "$checkme" ]
  then
      echo "STAR index already made (backup)"
  else
    sbatch -W ./star_index_pipeline.sh $STAR_index $genome $GTF
  fi
fi

#Run script that has STAR and featureCounts
#This shouldn't continue until all jobs are complete
if [ "$species" = "axolotl" ]; then
  var=(511_S10_L004 512_S11_L004 513_S12_L004 514_S13_L004 515_S14_L004 516_S15_L004 517_S16_L004 518_S17_L004 519_S18_L004 520_S19_L004 521_S20_L004 522_S21_L004 523_S22_L004 524_S23_L004 525_S24_L004 526_S25_L004 527_S26_L004 528_S27_L004 529_S28_L004 530_S29_L004 531_S30_L004)
elif [ "$species" = "terribilis" ]; then
  var=(409_S1_L001 410_S2_L001 411_S3_L001 412_S4_L001 413_S5_L001 414_S6_L001 415_S7_L001 416_S8_L001 417_S9_L001 418_S10_L001 419_S11_L001 420_S12_L001 421_S13_L001 422_S14_L001 423_S15_L001 424_S16_L001 425_S17_L001 426_S18_L001 427_S19_L001 428_S20_L001 429_S21_L001)
elif [ "$species" = "cane" ]; then
  var=(451_S43_L003 452_S44_L003 453_S45_L003 454_S46_L003 455_S47_L003 456_S48_L003 457_S49_L003 458_S50_L003 459_S51_L003 460_S52_L003 461_S53_L003 462_S54_L003 463_S55_L003 464_S56_L003 465_S57_L003 466_S58_L003 467_S59_L003 468_S60_L003 469_S61_L003 470_S62_L003 471_S63_L003)
elif [ "$species" = "bullfrog" ]; then
  var=(472_S64_L004 473_S65_L004 474_S66_L004 475_S67_L004 476_S68_L004 477_S69_L004 478_S70_L004 479_S71_L004 480_S72_L004 481_S73_L004 482_S74_L004 483_S75_L004 484_S76_L004 485_S77_L004 486_S78_L004 487_S79_L004 488_S80_L004 489_S81_L004 490_S82_L004 491_S83_L004 492_S84_L004)
elif [ "$species" = "xenopus" ]; then
  var=(430_S22_L002 431_S23_L002 432_S24_L002 433_S25_L002 434_S26_L002 435_S27_L002 436_S28_L002 437_S29_L002 438_S30_L002 439_S31_L002 440_S32_L002 441_S33_L002 442_S34_L002 443_S35_L002 444_S36_L002 445_S37_L002 446_S38_L002 447_S39_L002 448_S40_L002 449_S41_L002 450_S42_L002)
else
  echo "Unknown species: $species"
  exit 1
fi


if [ $go_to -gt 1 ]
then
  echo "STAR already complete"
else
  for var in "${var[@]}"
  do
      temp="$STAR_out""$var""$ending"
      F_READS="$path_in""$var""$forward"
      R_READS="$path_in""$var""$reverse"

      sbatch -J $jobname star_fc_pipeline.sh "$F_READS" "$R_READS" "$STAR_index" "$STAR_out""$var" "$GTF" "$FC_out" "$temp" "$var" &
  done
fi

#Make sure completely done
  # Wait for all jobs with job name "myjobname" to complete
while [ "$(squeue -u kwh1 -h -n $jobname | wc -l)" -gt 0 ]; do
  #while [ $(squeue -h -o "%j" -u $kwh1 | grep "starfc" | wc -l) -gt 0 ]; do
      sleep 60
done

#Run python script that combines featureCounts data into one big data table
if [ $go_to -gt 2 ]
then
  echo "Aren't you done already?"
else
  echo $species $FC_out $Bitter_TS $summary_table "$GTF"
  sbatch -W tabulate_expression2.py $species $FC_out $Bitter_TS $summary_table $GTF
fi
