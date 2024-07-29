#!/bin/bash
#SBATCH --job-name=checksubmit  # Job name
#SBATCH --mem=2gb          # Job memory request, down from 200 and closer to the 64gb I think you're using per instance
#SBATCH --nodes=1           # ensure cores are on one node
#SBATCH --ntasks=1          # run a single task
#SBATCH --cpus-per-task=1       # number of cores/threads requested, up from 4 you're asking for now - see too that I changed --runThreadN below to match
#SBATCH --partition=20
#SBATCH --output=logs/checksubmit.log # Standard output and error log
#SBATCH --error=logs/checksubmit.err

#Update your username in line 25
#Update minimum number of desired jobs to run in parallel in line 17

remaining="../../progress/all_accessions_remaining.txt"
less_than=25
to_add=5

counter=0

while [ $counter -lt 600 ]; do
    chmod +rwx $remaining
    # Count the number of running jobs for the specified user and job name
    job_count=$(squeue -u kwh1 -n pipe -h | wc -l)
    # splitC=$(squeue -u kwh1 -n split3 -h | wc -l)
    # bigB=$(squeue -u kwh1 -n blastBig -h | wc -l)
    # py1=$(squeue -u kwh1 -n py1B -h | wc -l)
    #job_count=$(($splitC + $bigB + $py1))

    # Check if the job count is less than 5
    if [ $job_count -lt $less_than ]; then
        # Count the number of lines with text in the file
        line_count=$(grep -c '[^[:space:]]' $remaining)
        (( counter++ ))
        # Check if the line count is 0
        if [ $line_count -eq 0 ]; then
            break
        fi

        # Adjust the parameter for sbatch command based on the line count
        if [ $line_count -lt $less_than ]; then
            ./add_new_jobs.sh $line_count $remaining
            sleep 3
            python ./progress.py
            #sbatch ./new_batch_notifier.sh
        else
            ./add_new_jobs.sh $to_add $remaining
            sleep 3
            python ./progress.py
            #sbatch ./new_batch_notifier.sh
        fi
    fi

    # Sleep for 2 minutes
    sleep 30

    ./progress.py
done
