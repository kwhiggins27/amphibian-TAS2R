#!/bin/bash
#SBATCH --job-name=pipe  # Job name
#SBATCH --mail-type=NONE     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=youremailaddress@yourinstitute   # Where to send mail
#SBATCH --mem=4gb          # Job memory request, down from 200 and closer to the 64gb I think you're using per instance
#SBATCH --nodes=1           # ensure cores are on one node
#SBATCH --ntasks=1          # run a single task
#SBATCH --cpus-per-task=1       # number of cores/threads requested, up from 4 you're asking for now - see too that I changed --runThreadN below to match
#SBATCH --partition=20
#SBATCH --output=/lab/wengpj01/vertebrate_pipeline/logs/pipe_%j.log # Standard output and error log
#SBATCH --error=/lab/wengpj01/vertebrate_pipeline/logs/pipe_%j.err

shopt -s extglob

accession="$1"

##read config file.  Create unix and python configs, and populate unix config
config_full="../../configs/config_full_$accession.txt"
config_unix="../../configs/config_unix_$accession.txt"
config_unix2="../../configs/config_unix_${accession}_2.txt"
for_py="../../configs/for_py2_$accession.py"

remaining="../../progress/all_accessions_remaining_rerun.txt"
complete="../../progress/list_complete_rerun.txt"
attempted="../../progress/list_attempted_rerun.txt"
failed="../../progress/list_failed_rerun.txt"

now2=$(date +%Y%m%d%H%M%S)

sed -i "/$accession/d" "$remaining"

sed '/##Variables for python/Q' $config_full > $config_unix

##Load unix variables
source $config_unix

sed '/#do not add anything after this point/q' $config_unix > $config_unix2
echo "directory=\"$directory\"/" >> $config_unix2

source $config_unix2

touch $attempted
echo $accession >> $attempted
#
# set -e
set -o pipefail

# Custom function to accession to list_failed if pipelinefails
fail_condition() {
    echo "$accession" >> "$failed"
    sed -i "/$accession/d" "$attempted"
    exit 1
}

#echo "$directory"
##generate python config file
echo "#!/usr/bin/env python" > $for_py
sed '1,/##Variables for python/d' $config_full >> $for_py
# echo "folder=\"$directory\"" >> $for_py

##create directory and navigate there
if [ $go_to -gt 0 ] ##whenever this loop appears, it's allowing alternate entry points (in case of restart)
then
  cd $existing_directory
  source $existing_directory/config_unix_species.sh
else
  mkdir -p $directory || fail_condition
  cd $directory
  cp "$config_unix2" "$directory""/config_unix_species.sh"
fi

##generate log file
big_log="log_$BLAST_species"_"$now2.log"
echo "##Logfile for blast pipeline" > $big_log


##move config files into directory
if [ $go_to -gt 0 ]
then
  echo "Already dealt with config files previously"
else
  cd $directory || fail_condition
  cp $config_full $directory/log_config.txt
  cp $for_py $directory/for_py2.py
  #cp $big_log $directory
  mkdir -p TMHMM
fi



#if big genome, check to see if split version.  If not, make it
subdivided="${genome_location%/*}/${accession}_subdivided"
echo "subdivided=\"$subdivided\"" >> "config_unix_species.sh"

echo $subdivided
echo $big_genome" Is this a big genome?"

subdivided_test=$subdivided"/checkdone.txt"
splitter="split_try3_pipe.sh"

if [ "$big_genome" == "yes" ]
then
  #mkdir -p $subdivided
    if [ -f "$subdivided_test" ]
      then
        echo "Directory already made"
    else
      mkdir -p "$subdivided"
      job_id0=`sbatch -W  $splitter $directory $accession |sed 's/Submitted batch job //g'`
    fi
else
  echo "Not a big genome"
fi


#run forward blast in unix, write to log.  If big genome, use alternate version
echo "Starting blast" > $big_log
blastF="blastF.sh"
blast_big="blast_big_pipe.sh"

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
    echo "#SBATCH --mail-user=youremailaddress@yourinstitute" >> torunparallel.sh
    while read n || [[ -n $p ]]
    do
      echo "sbatch -p 20 $blast_big $directory $n $accession" >> torunparallel.sh
    done < $names

    mkdir -p blastF

    sbatch -W  torunparallel.sh || fail_condition

    while [ "$(squeue -u kwh1 -h -n blastBig | wc -l)" -gt 0 ]; do
      #while [ $(squeue -h -o "%j" -u $kwh1 | grep "starfc" | wc -l) -gt 0 ]; do
          sleep 60
    done


  else
    sbatch -W  $blastF $directory || fail_condition
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
  py1="bitter_pt1.py"
  py1_big="bitter_pt1_Big.py"
  if [ "$big_genome" == "yes" ]
  then
    if [ $go_to -gt 0 ]
    then
      echo "run py1_big with existing directory" >> $big_log
      mkdir -p pull
      sbatch -W  $py1_big $existing_directory || fail_condition
      echo "py1_big complete" >> $big_log
    else
      echo "run py1_big in new directory" >> $big_log
      mkdir -p pull
      sbatch -W  $py1_big $directory || fail_condition
      echo "py1_big complete" >> $big_log
    fi
  else
    if [ $go_to -gt 0 ]
    then
      echo "run py1 with existing directory" >> $big_log
      sbatch -W  $py1 $existing_directory || fail_condition
      # chmod +rwx $BLAST_species"_full_pull.sh"
      echo "py1 complete" >> $big_log
    else
      echo "run py1 in new directory" >> $big_log
      rm -f $pull
      rm -f "expanded_file.fasta"
      sbatch -W  $py1 $directory || fail_condition
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
    sorted_names=$subdivided"/names_sorted.txt"
    sort $names > $sorted_names
    # cat $names | parallel -pipe sbatch -W  $blast_big #$directory
    #sbatch -W  $blast_big $directory $names_line
    rm -f "expanded_file.fasta"

    echo "#!/bin/bash" > big_full_pull.sh
    echo "#SBATCH --partition=20" >> big_full_pull.sh
    echo "#SBATCH --mail-type=NONE" >> big_full_pull.sh
    echo "#SBATCH --mail-user=youremailaddress@yourinstitute" >> big_full_pull.sh
    while read n || [[ -n $p ]]
    do
      part_aa="pull/""$n""_full_pull.sh"
      chmod +rwx "$part_aa"
      echo "sbatch -W -p 20 $part_aa" >> big_full_pull.sh
    done < $sorted_names

    sbatch -W  big_full_pull.sh || fail_condition

    # parallel sbatch -W pull/ ::: *.sh
  else
    rm -f "expanded_file.fasta"
    pull=$BLAST_species"_full_pull.sh"
    sbatch -W  ./$pull || fail_condition
    echo "full_pull complete" >> $big_log
  fi
fi

#run second python script, write to log
#for each expanded region, find the longest ORF and add it (plus new coordinates) to table from py1.  Apply filters and remove duplicates.
if [ $go_to -gt 3 ]
then
  echo "already have py2" >> $big_log
else
  py2="bitter_pt2.py"
  if [ $go_to -gt 0 ]
  then
    sbatch -W  $py2 $existing_directory || fail_condition
    echo "py2 complete" >> $big_log
  else
    sbatch -W  $py2 $directory || fail_condition
    echo "py2 complete" >> $big_log
  fi
fi

##run reciprocal blast in unix, write to log
if [ $go_to -gt 4 ]
then
  echo "already have blastR" >> $big_log
else
  blastR="blastR.sh"
  if [ $go_to -gt 0 ]
  then
    sbatch -W  $blastR $existing_directory || fail_condition
    echo "blastR complete" >> $big_log
  else
    sbatch -W  $blastR $directory || fail_condition
    echo "blastR complete" >> $big_log
  fi
fi


##run third python script, write to log
##identify reciprocal blast hits that fall within coordinates and match them back to original forward blast hits.
##identify ORFs that may be artifically long (in frame M before true start M) and shorten to region of blast hit
##use tmhmm to validate that these hits have 7 TM regions
##create final fasta of confirmed genes and gtf to use for expression analysis
py3="bitter_pt3_TMbed.py"
if [ $go_to -gt 0 ]
then
  while true; do
    running=$(squeue -u kwh1 -n py3 -h | wc -l)
    if [ $running -gt 50 ]; then
      sleep 60
    else
      sbatch -W $py3 $existing_directory || fail_condition
      if [ $? -eq 0 ]; then
        break
      else
        exit 1  # Exit the script with a non-zero status
      fi
    fi
  done
else
  while true; do
    running=$(squeue -u kwh1 -n py3 -h | wc -l)
    if [ $running -gt 50 ]; then
      sleep 30
    else
      sbatch -W $py3 $directory || fail_condition
      if [ $? -eq 0 ]; then
        break
      else
        exit 1  # Exit the script with a non-zero status
      fi
    fi
  done
  echo "py3 complete" >> $big_log
fi

#all done!
echo "pipeline done" >> $big_log
touch $complete
#echo $accession >> $failed
echo $accession >> $complete


sed -i "/$accession/d" "$attempted"
