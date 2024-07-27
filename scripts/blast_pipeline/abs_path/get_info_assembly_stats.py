#!/usr/bin/env python
#SBATCH --job-name=getinfo  # Job name
#SBATCH --mail-type=END,FAIL     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=youremailaddress@yourinstitute   # Where to send mail
#SBATCH --mem=4gb          # Job memory request, down from 200 and closer to the 64gb I think you're using per instance
#SBATCH --nodes=1           # ensure cores are on one node
#SBATCH --ntasks=1          # run a single task
#SBATCH --cpus-per-task=1       # number of cores/threads requested, up from 4 you're asking for now - see too that I changed --runThreadN below to match
#SBATCH --partition=20
#SBATCH --output=getinfo_%j.log # Standard output and error log
#SBATCH --error=getinfo_%j.err


import csv
import os
import pandas as pd
import subprocess
from Bio import Entrez

# Set your email address (required by NCBI)
Entrez.email = 'UPDATE ME'
Entrez.api_key= 'UPDATE ME'

csv_file = '/lab/wengpj01/vertebrate_pipeline/toy_accessions.csv'  # Replace with your CSV file path
output = '/lab/wengpj01/genomes/toy_accession_stats.csv'


def fetch_stats(accession):
    command = (
        f"esearch -db assembly -query '{accession}' | "
        "esummary | "
        "xtract -pattern DocumentSummary -block Stat "
        "-if General@category -equals assembly_type "
        "-or Stat@category -equals contig_count "
        "-or Stat@category -equals contig_l50 "
        "-or Stat@category -equals contig_n50 "
        "-or Stat@category -equals total_length "
        "-element Stat"
    )
    output = subprocess.check_output(command, shell=True, text=True)
    stats_string = output.strip().split('\n')[0]  # Assuming only one line of output
    stats = stats_string.split('\t')

    print(stats)

    coverage_command = (
        f"esearch -db assembly -query '{accession}' | "
        "esummary | "
        "xtract -pattern DocumentSummary -element Coverage"
    )


    coverage = subprocess.check_output(coverage_command, shell=True, text=True).strip()

    stats.insert(0, accession)  # Insert accession at the beginning

    stats.append(coverage)

    return stats

def main():
    data = []
    with open(csv_file, 'r') as file:
        for line in file:
            accession = line.strip()
            stats = fetch_stats(accession)
            data.append(stats)

    #columns = ['accession','contig_count', 'contig_l50', 'contig_n50', 'total_length', 'coverage']
    columns = ['accession','assembly_type','contig_count', 'contig_l50', 'contig_n50', 'total_length', 'coverage']
    df = pd.DataFrame(data, columns=columns)
    df.to_csv(output, index=False)

if __name__ == '__main__':
    main()







# def append_stats_to_csv(output, csv_file):
#     stats = output.strip().split('\t')
#     with open(csv_file, 'a', newline='') as file:
#         writer = csv.writer(file)
#         writer.writerow(stats)
#
# # Create an empty list to store results
# results = []
#
# # Read the CSV file
# with open(csv_file, 'r') as file:
#     reader = csv.reader(file)
#     for row in reader:
#         if row:  # Skip empty lines
#             accession = row[0].strip()
#             print(accession)
#
#             # Construct the command
#             command = (
#                 f"esearch -db assembly -query '{accession}' | "
#                 "esummary | "
#                 "xtract -pattern DocumentSummary -block Stat \
#         -if Stat@category -equals contig_count -or Stat@category -equals contig_l50 \
#         -or Stat@category -equals contig_n50 -or Stat@category -equals total_length \
#         -sep \":\" -def \"NA\" -element Stat"
#             )
#
#             # Run the command and capture the output
#             species = os.popen(command).read().strip()
#
#             # Extract latin and common variables from species
#             if '(' in species and ')' in species:
#                 start_index = species.index('(')
#                 end_index = species.index(')')
#                 latin = species[:start_index].rstrip().replace(' ', '_')
#                 common = species[start_index + 1:end_index].strip().replace(' ', '_')
#             else:
#                 latin = species.replace(' ', '_')
#                 common = ''
#             # Add the results to the list
#             results.append({
#                 'accession': accession,
#                 'latin': latin,
#                 'common': common,
#             })
#
# # Create a dataframe from the results
# df = pd.DataFrame(results)
#
# # Save the dataframe as a new CSV file
# output_csv = output  # Replace with your desired output file path
# df.to_csv(output_csv, index=False)
#
# print("CSV file saved successfully!")
