#!/bin/bash
#SBATCH --job-name=bull # Job name
#SBATCH --mail-type=END,FAIL     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kwh1@wi.mit.edu   # Where to send mail
#SBATCH --mem=4gb          # Job memory request, down from 200 and closer to the 64gb I think you're using per instance
#SBATCH --nodes=1           # ensure cores are on one node
#SBATCH --ntasks=1          # run a single task
#SBATCH --cpus-per-task=1       # number of cores/threads requested, up from 4 you're asking for now - see too that I changed --runThreadN below to match
#SBATCH --partition=20
#SBATCH --output=logs/pipe_%j.log # Standard output and error log
#SBATCH --error=logs/pipe_%j.err

shopt -s extglob

##read config file.  Create unix and python configs, and populate unix config
config_full="/lab/solexa_weng/playground/Kate_Higgins/other_vertebrates/pipeline/config_pipeline_bullfrog.sh"
config_unix="/lab/solexa_weng/playground/Kate_Higgins/other_vertebrates/pipeline/config_unix.sh"
config_unix2="/lab/solexa_weng/playground/Kate_Higgins/other_vertebrates/pipeline/config_unix2.sh"
for_py="/lab/solexa_weng/playground/Kate_Higgins/other_vertebrates/pipeline/for_py2.py"

sed '/##Variables for python/Q' $config_full > $config_unix

##Load unix variables
source $config_unix

sed '/#do not add anything after this point/q' $config_unix > $config_unix2
echo "directory=\"$directory\"/" >> $config_unix2
z
source $config_unix2

set -e
set -o pipefail

#echo "$directory"
##generate python config file
echo "#!/usr/bin/env python" > $for_py
sed '1,/##Variables for python/d' $config_full >> $for_py
echo "folder=\"$directory\"" >> $for_py

##create directory and navigate there
if [ $go_to -gt 0 ] ##whenever this loop appears, it's allowing alternate entry points (in case of restart)
then
  cd $existing_directory
  source $existing_directory/config_unix_species.sh
else
  mkdir $directory
  cd $directory
  cp "$config_unix2" "$directory""/config_unix_species.sh"
fi

##generate log file
big_log="log_$BLAST_species"_"$now.log"
echo "##Logfile for Kate's blast pipeline" > $big_log


##move config files into directory
if [ $go_to -gt 0 ]
then
  echo "Already dealt with config files previously"
else
  cd $directory
  cp $config_full $directory/log_config.txt
  mv $for_py $directory
  #cp $big_log $directory
  mkdir TMHMM
fi



#if big genome, check to see if split version.  If not, make it
subdivided="/lab/wengpj01/"$BLAST_species"/subdivided"
echo "subdivided=\"$subdivided\"" >> "config_unix_species.sh"

echo $subdivided
echo $big_genome" Is this a big genome?"

subdivided_test=$subdivided"/checkdone.txt"
splitter="/lab/solexa_weng/playground/Kate_Higgins/other_vertebrates/pipeline/split.sh"

if [ "$big_genome" == "yes" ]
then
  #mkdir -p $subdivided
    if [ -f "$subdivided_test" ]
      then
        echo "Directory already made"
    else
      mkdir -p "$subdivided"
      job_id0=`sbatch -W  $splitter $directory |sed 's/Submitted batch job //g'`
    fi
else
  echo "Not a big genome"
fi


#run forward blast in unix, write to log.  If big genome, use alternate version
echo "Starting blast" > $big_log
blastF="/lab/solexa_weng/playground/Kate_Higgins/other_vertebrates/pipeline/blastF.sh"
blast_big="/lab/solexa_weng/playground/Kate_Higgins/other_vertebrates/pipeline/blast_big.sh"

if [ $go_to -gt 0 ]
then
  echo "already have blast"
else
  if [ "$big_genome" == "yes" ]
  then
    names=$subdivided"/names.txt"
    # cat $names | parallel -pipe sbatch -W  $blast_big #$directory
    #sbatch -W  $blast_big $directory $names_line

    echo "#!/bin/bash" > torunparallel.sh
    echo "#SBATCH --partition=20" >> torunparallel.sh
    echo "#SBATCH --mail-type=NONE" >> torunparallel.sh
    echo "#SBATCH --mail-user=kwh1@wi.mit.edu" >> torunparallel.sh
    while read n || [[ -n $p ]]
    do
      echo "sbatch -p 20 $blast_big $directory $n" >> torunparallel.sh
    done < $names

    mkdir blastF

    sbatch -W  torunparallel.sh

    while [ "$(squeue -u kwh1 -h -n blastBig | wc -l)" -gt 0 ]; do
      #while [ $(squeue -h -o "%j" -u $kwh1 | grep "starfc" | wc -l) -gt 0 ]; do
          sleep 60
    done


  else
    sbatch -W  $blastF $directory
    echo "Forward blast complete" >> $big_log
  fi
fi



# ##run first python script, write to log.  If big genome, use alternate version
## In brief, py1 makes dataframe out of blast hits and defines region around each hit.  It then writes the script full_pull.
echo $go_to

if [ $go_to -gt 1 ]
then
  echo "already have first py"
  if [ "$big_genome" == "no" ]
  then
    chmod +rwx $BLAST_species"_full_pull.sh"
  fi
else
  echo "need to run py1"
  py1="/lab/solexa_weng/playground/Kate_Higgins/other_vertebrates/pipeline/bitter_pt1.py"
  py1_big="/lab/solexa_weng/playground/Kate_Higgins/other_vertebrates/pipeline/bitter_pt1_Big.py"
  if [ "$big_genome" == "yes" ]
  then
    if [ $go_to -gt 0 ]
    then
      echo "run py1_big with existing directory" >> $big_log
      mkdir -p pull
      sbatch -W  $py1_big $existing_directory
      echo "py1_big complete" >> $big_log
    else
      echo "run py1_big in new directory" >> $big_log
      mkdir -p pull
      sbatch -W  $py1_big $directory
      echo "py1_big complete" >> $big_log
    fi
  else
    if [ $go_to -gt 0 ]
    then
      echo "run py1 with existing directory" >> $big_log
      sbatch -W  $py1 $existing_directory
      # chmod +rwx $BLAST_species"_full_pull.sh"
      echo "py1 complete" >> $big_log
    else
      echo "run py1 in new directory" >> $big_log
      sbatch -W  $py1 $directory
      chmod +rwx $BLAST_species"_full_pull.sh"
      echo "py1 complete" >> $big_log
    fi
  fi
fi

##unix: full_pull.  Get bases immediately up and downstream from each blast hit
if [ $go_to -gt 2 ]
then
  echo "already have full pull" >> $big_log
else
  if [ "$big_genome" == "yes" ]
  then

    names=$subdivided"/names.txt"
    # cat $names | parallel -pipe sbatch -W  $blast_big #$directory
    #sbatch -W  $blast_big $directory $names_line

    echo "#!/bin/bash" > big_full_pull.sh
    echo "#SBATCH --partition=20" >> big_full_pull.sh
    echo "#SBATCH --mail-type=NONE" >> big_full_pull.sh
    echo "#SBATCH --mail-user=kwh1@wi.mit.edu" >> big_full_pull.sh
    while read n || [[ -n $p ]]
    do
      part_aa="pull/""$n""aa_full_pull.sh"
      part_ab="pull/""$n""ab_full_pull.sh"
      chmod +rwx "$part_aa"
      chmod +rwx "$part_ab"
      echo "sbatch -W -p 20 $part_aa" >> big_full_pull.sh
      echo "sbatch -W -p 20 $part_ab" >> big_full_pull.sh
    done < $names

    sbatch -W  big_full_pull.sh

    # parallel sbatch -W pull/ ::: *.sh
  else
    pull=$BLAST_species"_full_pull.sh"
    sbatch -W  ./$pull
    echo "full_pull complete" >> $big_log
  fi
fi

#run second python script, write to log
#for each expanded region, find the longest ORF and add it (plus new coordinates) to table from py1.  Apply filters and remove duplicates.
if [ $go_to -gt 3 ]
then
  echo "already have py2" >> $big_log
else
  py2="/lab/solexa_weng/playground/Kate_Higgins/other_vertebrates/pipeline/bitter_pt2.py"
  if [ $go_to -gt 0 ]
  then
    sbatch -W  $py2 $existing_directory
    echo "py2 complete" >> $big_log
  else
    sbatch -W  $py2 $directory
    echo "py2 complete" >> $big_log
  fi
fi

##run reciprocal blast in unix, write to log
if [ $go_to -gt 4 ]
then
  echo "already have blastR" >> $big_log
else
  blastR="/lab/solexa_weng/playground/Kate_Higgins/other_vertebrates/pipeline/blastR.sh"
  if [ $go_to -gt 0 ]
  then
    sbatch -W  $blastR $existing_directory
    echo "blastR complete" >> $big_log
  else
    sbatch -W  $blastR $directory
    echo "blastR complete" >> $big_log
  fi
fi


##run third python script, write to log
##identify reciprocal blast hits that fall within coordinates and match them back to original forward blast hits.
##identify ORFs that may be artifically long (in frame M before true start M) and shorten to region of blast hit
##use tmhmm to validate that these hits have 7 TM regions
##create final fasta of confirmed genes and gtf to use for expression analysis
py3="/lab/solexa_weng/playground/Kate_Higgins/other_vertebrates/pipeline/bitter_pt3_TMbed_forv6.py"
if [ $go_to -gt 0 ]
then
  sbatch -W  $py3 $existing_directory
  echo "py3 complete" >> $big_log
else
  sbatch -W  $py3 $directory
  echo "py3 complete" >> $big_log
fi


#all done!
echo "pipeline done" >> $big_log
