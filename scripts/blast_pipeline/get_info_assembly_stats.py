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

# Need to update entrez email and api key in lines 23 and 24


import csv
import os
import pandas as pd
import subprocess
from Bio import Entrez

# Set your email address (required by NCBI)
Entrez.email = 'UPDATE ME'
Entrez.api_key= 'UPDATE ME'

csv_file = '../../progress/toy_accessions.csv'  # Replace with your CSV file path
output = '../../progress/toy_accession_stats.csv'


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
