#!/usr/bin/env python
#SBATCH --job-name=randomgen   # Job name
#SBATCH --mem=100gb                     # Job memory request, down from 200 and closer to the 64gb I think you're using per instance
#SBATCH --nodes=1                      # ensure cores are on one node
#SBATCH --ntasks=1                     # run a single task
#SBATCH --cpus-per-task=1              # number of cores/threads requested, up from 4 you're asking for now - see too that I changed --runThreadN below to match
#SBATCH --partition=20
#SBATCH --output=logs/randomgen_%j.log   # Standard output and error log
#SBATCH --error=logs/randomgen_%j.err

import os
import pandas as pd
import random

# Function to generate random windows for a single accession
def generate_random_windows(accession):
    # Define the path to the chromsizes.csv file
    chromsizes_path = f"../../subdirs/{accession}/chromsizes.csv"

    if not os.path.exists(chromsizes_path):
        print(f"Chromsizes file not found for accession: {accession}. Skipping.")
        return

    # Specify column names
    column_names = ["chromosome", "length", "ignore"]

    # Read the chromsizes.csv file into a DataFrame with column names
    chromsizes = pd.read_csv(chromsizes_path, sep="\t", names=column_names)

 # Check if the resulting DataFrame contains only False values
    if chromsizes.all().all() == False:
        print(f"DataFrame for {chromsizes_path} contains only False values. Skipping to the next item.")
        return

    # Initialize variables
    running_total_length = 0
    long_enough = []
    largest_start = 0
    for_bed = pd.DataFrame(columns=["chromosome", "start", "stop"])

    # Fill "before" and "after" columns with empty values
    chromsizes["before"] = ""
    chromsizes["after"] = ""

    # Iterate through each chromosome
    for index, row in chromsizes.iterrows():
        chromosome = row["chromosome"]
        length = row["length"]

        if length > 100000:
            long_enough.append(chromosome)
            before = running_total_length
            running_total_length += length
            after = running_total_length

            # Update largest_start if needed
            if largest_start < running_total_length - 100000:
                largest_start = running_total_length - 100000

            chromsizes.at[index, "before"] = before
            chromsizes.at[index, "after"] = after

    chromsizes["before"] = pd.to_numeric(chromsizes["before"])
    chromsizes["after"] = pd.to_numeric(chromsizes["after"])
    print(running_total_length)
    # Generate 10 random numbers
    random_starts = []
    random_loc = []
    # for _ in range(10):
    #     while True:
    #         i = random.randint(0, largest_start)
    #         j = chromsizes[(chromsizes["before"] < i) & (chromsizes["after"] >= i)]
    #         if not j.empty and (j["length"].values[0] - 100000 >= (i - j["after"]).values[0]):
    #             random_starts.append(i)
    #             break
    #
    for _ in range(10):
        while True:
            i = random.randint(0, running_total_length-100000)
            j = chromsizes[(chromsizes["before"] < i) & (chromsizes["after"] >= i)]
            # Initialize k
            k = None

            # Find the row before j
            if not j.empty:
                j_index = j.index[0]
                if j_index > 0:  # Check if j_index is greater than 0
                    k = chromsizes.loc[j_index - 1]
                else:
                    k = None  # Set k to None when there's no previous row

                if (j["length"].values[0] - 100000) >= (i - j["before"].values[0]):
                    random_starts.append(i)
                    random_loc.append(i - j["before"].values[0])
                    break



    # Create the for_bed DataFrame
    for index, start in enumerate(random_starts):
        chromosome = chromsizes.loc[(chromsizes["before"] <= start) & (chromsizes["after"] > start), "chromosome"].values[0]
        start = random_loc[index]
        stop = start + 100000
        for_bed.loc[index] = [chromosome, start, stop]

    # Save for_bed as a bed file
    for_bed.to_csv(f"../../subdirs/{accession}/random_windows.bed", sep='\t', header=False, index=False)

    print("yay!")
# Read accession numbers from accessions_to_keep.txt
# accession_file = "/lab/wengpj01/vertebrate_pipeline/big_genomes.txt"
accession_file = "../../results/all_accessions_used.txt"
with open(accession_file, 'r') as f:
    accessions = [line.strip() for line in f]

# Process each accession
for accession in accessions:
    generate_random_windows(accession)
    print(f"Completed:{accession}")
