#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import glob
import sys
import os

accession="GCA_000001405.29"
threshold_short=1000000
threshold_ends=0.1

new_directory = f"/../../subdirs/{accession}"
os.chdir(new_directory)


# Read chromsizes.csv into a pandas dataframe
chromsizes = pd.read_csv('chromsizes.csv', sep='\t', header=None, names=['chr', 'length'])

# Find and read all *_from_pipeline.gtf files into a list of dataframes
gtf_files = glob.glob('*_from_pipeline.gtf')
gene_dfs = []

for gtf_file in gtf_files:
    df = pd.read_csv(gtf_file, sep='\t', header=None, names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])
    gene_dfs.append(df)

    # Strip text before the last underscore in the 'seqname' column
    df['seqname'] = df['seqname'].str.replace(r'^.*_', '', regex=True)

# Concatenate all dataframes into a single dataframe
genes = pd.concat(gene_dfs, ignore_index=True)

# Iterate through rows of the 'genes' dataframe and find corresponding chromosome length
genes['chr_length'] = None  # Create an empty column to store chromosome lengths
genes['chr_mini'] = None
genes['ends'] = None

for i, gene_row in genes.iterrows():
    seqname = gene_row['seqname']

    # Find the corresponding row in 'chromsizes'
    chrom_row = chromsizes[chromsizes['chr'] == seqname]

    if not chrom_row.empty:
        chr_length = chrom_row.iloc[0]['length']
        genes.at[i, 'chr_length'] = chr_length

        to_end = (chr_length - gene_row["start"])/chr_length
        to_start = (gene_row["start"])/chr_length

        if (to_start < threshold_ends) or (to_end < threshold_ends):
            genes.at[i, 'ends'] = "True"
        else:
            genes.at[i, 'ends'] = "False"

        if chr_length < threshold_short:
            genes.at[i, 'chr_mini'] = "True"
        else:
            genes.at[i, 'chr_mini'] = "False"


# Calculate the fractions
fraction_chr_mini_true = len(genes[genes['chr_mini'] == "True"]) / len(genes)
fraction_near_ends = len(genes[genes['ends'] == "True"]) / len(genes)
fraction_near_ends_long = len(genes[(genes['ends'] == "True") & (genes['chr_mini'] == "False")]) / len(genes[genes['chr_mini'] == "False"])

# Save the results to a file
output_path = "../../results/coordinate_analysis/near_end.txt"
with open(output_path, 'w') as f:
    f.write(f"Fraction of genes on short chromosomes (chr_mini=True): {fraction_chr_mini_true}\n")
    f.write(f"Fraction of genes near chromosome ends (ends=True): {fraction_near_ends}\n")
    f.write(f"Fraction of genes near ends on long chromosomes (ends=True and chr_mini=False): {fraction_near_ends_long}\n")
