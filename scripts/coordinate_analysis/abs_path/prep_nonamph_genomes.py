#!/usr/bin/env python
import pandas as pd
import os
import random
# import importlib
import sys
import re

# Read the dataframe from match_accession_class.csv
df = pd.read_csv("/lab/wengpj01/vertebrate_pipeline/20231006_run/match_accession_class.csv")

selected_accessions = []

for class_name in df['class'].unique():
    if class_name != 'Amphibia':
        # Get Accessions for the current class
        accessions = df[df['class'] == class_name]['Accession'].tolist()

        # Randomly select up to four values, or all if fewer than four
        num_to_select = min(4, len(accessions))
        selected_accessions.extend(random.sample(accessions, num_to_select))
        # print(selected_accessions)

# print(selected_accessions)
# Specify the output file path
output_file_path = "/lab/wengpj01/repeat/non_amph_pt1.txt"
# Create a set to store unique entries
unique_entries = set()

for accession in selected_accessions:
    # Change directory to /lab/wengpj01/vertebrate_pipeline/subdirs/Accession
    directory_path = f"/lab/wengpj01/vertebrate_pipeline/subdirs/{accession}"
    os.chdir(directory_path)
    sys.path.insert(0,directory_path)

    # Initialize variables to store common and genome_location values
    common = None
    genome_location = None

    # Open the file for_py2.py and extract the values of "common" and "genome_location"
    with open("for_py2.py", "r") as for_py2_file:
        for line in for_py2_file:
            common_match = re.match(r'\s*common\s*=\s*"([^"]+)"', line)
            genome_location_match = re.match(r'\s*genome_location\s*=\s*"([^"]+)"', line)

            if common_match:
                common = common_match.group(1)
            elif genome_location_match:
                genome_location = genome_location_match.group(1)

            # If both values are found, exit the loop
            if common and genome_location:
                break

    print(common)

    # # Source the file for_py2.py to access the variables
    # importlib.import_module('for_py2')
    # from for_py2 import *

    # Build the required command
    command = f"BuildDatabase -name references/{common} {genome_location}"

    # Add the command to the set of unique entries
    unique_entries.add(command)

output_file_path1 = "/lab/wengpj01/repeat/non_amph_pt2.txt"
# After processing all accessions, write unique entries to the output file
with open(output_file_path1, "a") as output_file:
    # Append the nohup command for each unique entry
    for entry in unique_entries:
        nohup_script = f"nohup RepeatModeler -database {entry} -threads 10 -LTRStruct >& run.out &\n"
        output_file.write(nohup_script)
