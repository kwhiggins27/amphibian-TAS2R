#!/bin/bash
#SBATCH --job-name=checksubmit  # Job name
#SBATCH --mem=2gb          # Job memory request, down from 200 and closer to the 64gb I think you're using per instance
#SBATCH --nodes=1           # ensure cores are on one node
#SBATCH --ntasks=1          # run a single task
#SBATCH --cpus-per-task=1       # number of cores/threads requested, up from 4 you're asking for now - see too that I changed --runThreadN below to match
#SBATCH --partition=20
#SBATCH --output=logs/checksubmit.log # Standard output and error log
#SBATCH --error=logs/checksubmit.err

remaining="../../results/busco/busco_needed.txt"
less_than=35
to_add=5

counter=0

while [ $counter -lt 600 ]; do
    chmod +rwx $remaining
    # Count the number of running jobs for the specified user and job name
    job_count=$(squeue -u kwh1 -n busco -h | wc -l)

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
            ./add_busco.sh $line_count $remaining
            sleep 3
        else
            ./add_busco.sh $to_add $remaining
            sleep 3
        fi
    fi

    # Sleep for 2 minutes
    sleep 30

done
