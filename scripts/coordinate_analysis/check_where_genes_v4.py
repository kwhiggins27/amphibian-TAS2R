#!/usr/bin/env python
#SBATCH --job-name=localize  # Job name
#SBATCH --mail-type=NONE     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=youremailaddress@yourinstitute # Where to send mail
#SBATCH --mem=50gb          # Job memory request, down from 200 and closer to the 64gb I think you're using per instance
#SBATCH --nodes=1           # ensure cores are on one node
#SBATCH --ntasks=1          # run a single task
#SBATCH --cpus-per-task=10       # number of cores/threads requested, up from 4 you're asking for now - see too that I changed --runThreadN below to match
#SBATCH --partition=20
#SBATCH --output=logs/localize_%j.log # Standard output and error log
#SBATCH --error=logs/localize_%j.err

import pandas as pd
import glob
import sys
import os

threshold_short = 1000000
threshold_ends = 0.1

output_file = "../../results/coordinate_analysis/position_within_chromosome_1021.csv"
output_file1 = "../../results/coordinate_analysis/position_within_chromosome_all_genes_1021.csv"

# Read the list of accessions from all_accessions.txt
with open('../../results/all_accessions_used.txt', 'r') as accession_file:
    accessions = accession_file.read().splitlines()

# Iterate through each accession
for accession in accessions:
    new_directory = f"../../subdirs/{accession}"
    os.chdir(new_directory)

    # Read chromsizes.csv into a pandas dataframe
    chromsizes = pd.read_csv('chromsizes_1007.csv', sep='\t', header=None, names=['chr', 'length','mini'])
    print(chromsizes)

    # Find and read all *_from_pipeline.gtf files into a list of dataframes
    gtf_files = glob.glob('*_from_pipeline.gtf')
    gene_dfs = []

    # print(gtf_files)

    # Skip processing if no GTF files are found
    if not gtf_files:
        print(f"No GTF files found for accession {accession}. Skipping.")
        continue

    for gtf_file in gtf_files:
        df = pd.read_csv(gtf_file, sep='\t', header=None, names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])

        # Strip text before the last underscore in the 'seqname' column
        df['seqname'] = df['seqname'].str.replace(r'^.*_', '', regex=True)

        gene_dfs.append(df)

    # Concatenate all dataframes into a single dataframe
    genes = pd.concat(gene_dfs, ignore_index=True)

    # Iterate through rows of the 'genes' dataframe and find corresponding chromosome length
    genes['chr_length'] = None  # Create an empty column to store chromosome lengths
    genes['chr_mini'] = None
    genes['ends'] = None
    genes['from_end'] = None

    for i, gene_row in genes.iterrows():
        seqname = gene_row['seqname']

        # Find the corresponding row in 'chromsizes'
        chrom_row = chromsizes[chromsizes['chr'] == seqname]

        if not chrom_row.empty:
            chr_length = chrom_row.iloc[0]['length']
            chr_mini = chrom_row.iloc[0]['mini']
            genes.at[i, 'chr_length'] = chr_length
            genes.at[i, 'chr_mini'] = chr_mini
            i_mini = chr_mini

            if not chrom_row.empty:
                chr_length = chrom_row.iloc[0]['length']
                chr_mini = chrom_row.iloc[0]['mini']
                genes.at[i, 'chr_length'] = chr_length
                genes.at[i, 'chr_mini'] = chr_mini

                try:
                    to_end = (chr_length - gene_row["start"]) / chr_length
                    to_start = (gene_row["start"]) / chr_length

                    if to_start < to_end:
                        genes.at[i, 'from_end'] = to_start
                        i_ends = to_start
                    else:
                        genes.at[i, 'from_end'] = to_end
                        i_ends = to_end

                    if (to_start < threshold_ends) or (to_end < threshold_ends):
                        genes.at[i, 'ends'] = "True"
                    else:
                        genes.at[i, 'ends'] = "False"

                    i_chromosome = gene_row["seqname"]
                    i_start = gene_row["start"]
                    i_stop = gene_row["end"]
                    new_row1 = f"{accession},{i_chromosome}, {i_start}, {i_stop},{i_mini},{i_ends}\n"

                    with open(output_file1, 'a') as file:
                        file.write(new_row1)
                except ZeroDivisionError:
                    print(f"ZeroDivisionError for accession {accession}. Skipping this entry.")
                    continue

    # Calculate the fraction of rows where 'chr_mini' is True

    length_genes=len(genes)

    if length_genes > 0:
        fraction_chr_mini_true = len(genes[genes['chr_mini'] == "TRUE"]) / len(genes)

        # Calculate the fraction of rows where 'ends' is True
        fraction_near_ends = len(genes[genes['ends'] == "TRUE"]) / len(genes)

        # Calculate the fraction of rows where 'ends' is True and 'chr_mini' is False
        denominator = len(genes[genes['chr_mini'] == "FALSE"])

        if denominator > 0:
            fraction_near_ends_long = len(genes[(genes['ends'] == "TRUE") & (genes['chr_mini'] == "FALSE")]) / denominator
        else:
            fraction_near_ends_long = 0.0
    else:
        continue
    # Check if the file position_within_chromosome.csv exists


    if not os.path.isfile(output_file):
        # If it doesn't exist, create the file with header
        with open(output_file, 'w') as file:
            file.write("accession,fraction_chr_mini_true,fraction_near_ends,fraction_near_ends_long\n")

    # Append the information as a new row
    new_row = f"{accession},{fraction_chr_mini_true},{fraction_near_ends},{fraction_near_ends_long}\n"

    with open(output_file, 'a') as file:
        file.write(new_row)
