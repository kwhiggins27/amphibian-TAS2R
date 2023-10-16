#!/usr/bin/env python
import os
import csv

# Define the directory path
directory_path = '/lab/wengpj01/repeat/amphibian_near_TAS2R'

# Initialize dictionaries to store percentages for each category
percentages_singleton = {}
percentages_cluster = {}
percentages_random = {}

# Loop through the files in the directory
for filename in os.listdir(directory_path):
    if filename.endswith('singleton.fasta.tbl'):
        category = 'singleton'
    elif filename.endswith('cluster.fasta.tbl'):
        category = 'cluster'
    elif filename.endswith('random.fasta.tbl'):
        category = 'random'
    else:
        continue  # Skip files that don't match the expected pattern

    common_name = filename.split('_')[0]
    with open(os.path.join(directory_path, filename), 'r', encoding='latin-1', errors='ignore') as file:
        for line in file:
            # if "bases masked" in line:
            #     percentage = float(line.split('(')[1].split('%')[0].strip())
            #     if category == 'singleton':
            #         percentages_singleton[common_name] = percentage
            #     elif category == 'cluster':
            #         percentages_cluster[common_name] = percentage
            #     elif category == 'random':
            #         percentages_random[common_name] = percentage
            if line.startswith("LTR elements:"):
                percentage = float(line.split()[-2])
                if category == 'singleton':
                    percentages_singleton[common_name] = percentage
                elif category == 'cluster':
                    percentages_cluster[common_name] = percentage
                elif category == 'random':
                    percentages_random[common_name] = percentage

# # Print the table header
# print("common_name percentage_singleton percentage_cluster percentage_random")
#
# # Loop through the common names and print the percentages
# for common_name in sorted(set(percentages_singleton.keys()) | set(percentages_cluster.keys()) | set(percentages_random.keys())):
#     print(f"{common_name} {percentages_singleton.get(common_name, 0)} {percentages_cluster.get(common_name, 0)} {percentages_random.get(common_name, 0)}")
# Define the CSV file name
csv_file = 'amphibians_LTR.csv'

# Write the table to a CSV file
with open(csv_file, 'w', newline='') as csvfile:
    fieldnames = ["common_name", "percentage_singleton", "percentage_cluster", "percentage_random"]
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    writer.writeheader()
    for common_name in sorted(set(percentages_singleton.keys()) | set(percentages_cluster.keys()) | set(percentages_random.keys())):
        writer.writerow({
            "common_name": common_name,
            "percentage_singleton": percentages_singleton.get(common_name, 0),
            "percentage_cluster": percentages_cluster.get(common_name, 0),
            "percentage_random": percentages_random.get(common_name, 0)
        })
