#!/usr/bin/env python
import pandas as pd
import os

# Read the list of accessions from accessions_to_keep.txt
accession_file = '../../results/accessions_mini_run.txt'
with open(accession_file, 'r') as file:
    accessions = [line.strip() for line in file]

# Iterate through each accession and create the expanded BED file
for accession in accessions:
    # Create the input BED file path
    bed_file = f'../../subdirs/{accession}/singletons.bed'

        # Check if the file exists before proceeding
    if not os.path.exists(bed_file):
        print(f"File not found for accession: {accession}. Skipping to the next accession.")
        continue

    # Define the column names
    columns = ['chromosome', 'start', 'stop']

    # Read the BED file into a Pandas DataFrame
    df = pd.read_csv(bed_file, sep='\t', names=columns)

    # Create a new DataFrame to store the expanded entries
    expanded_df = pd.DataFrame(columns=columns)

    # Iterate through each line of the original DataFrame
    for index, line in df.iterrows():
        # Create a new row with start = max(0, line["start"] - 100000) and stop = line["start"]
        new_entry1 = {
            'chromosome': line['chromosome'],
            'start': max(0, line['start'] - 1000000),
            'stop': line['start']
        }
        new_entry1_df = pd.DataFrame([new_entry1])
        expanded_df = pd.concat([expanded_df, new_entry1_df], ignore_index=True)

        # Create a new row with start = line["start"] and stop = line["stop"] + 100000
        new_entry2 = {
            'chromosome': line['chromosome'],
            'start': line['stop'],
            'stop': line['stop'] + 1000000
        }
        new_entry2_df = pd.DataFrame([new_entry2])
        expanded_df = pd.concat([expanded_df, new_entry2_df], ignore_index=True)

    # Define the output BED file path
    output_file = f'../../subdirs/{accession}/singletons_margins_1M.bed'

    # Save the expanded DataFrame as a BED file
    expanded_df.to_csv(output_file, sep='\t', header=False, index=False)

    # Create the input BED file path
    bed_file = f'../../subdirs/{accession}/clusters.bed'

    # Define the column names
    columns = ['chromosome', 'start', 'stop']

    # Read the BED file into a Pandas DataFrame
    df = pd.read_csv(bed_file, sep='\t', names=columns)

    # Create a new DataFrame to store the expanded entries
    expanded_df = pd.DataFrame(columns=columns)

    # Iterate through each line of the original DataFrame
    for index, line in df.iterrows():
        # Create a new row with start = max(0, line["start"] - 100000) and stop = line["start"]
        new_entry1 = {
            'chromosome': line['chromosome'],
            'start': max(0, line['start'] - 1000000),
            'stop': line['start']
        }
        expanded_df = pd.concat([expanded_df, new_entry1_df], ignore_index=True)

        # Create a new row with start = line["start"] and stop = line["stop"] + 100000
        new_entry2 = {
            'chromosome': line['chromosome'],
            'start': line['stop'],
            'stop': line['stop'] + 1000000
        }
        expanded_df = pd.concat([expanded_df, new_entry2_df], ignore_index=True)


    # Define the output BED file path
    output_file = f'../../subdirs/{accession}/clusters_margins_1M.bed'
    # Save the expanded DataFrame as a BED file
    expanded_df.to_csv(output_file, sep='\t', header=False, index=False)
