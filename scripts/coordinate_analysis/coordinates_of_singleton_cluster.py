#!/usr/bin/env python
import pandas as pd
import os
import glob

# Define a function to process a single accession
def process_accession(accession):
    try:
        # Set the working directory
        data_directory = f"../../subdirs/{accession}"
        os.chdir(data_directory)

        # Define the file paths
        clusters_file_path = "clusters_found_1000000.csv"
        gtf_files = glob.glob(f"*_from_pipeline.gtf")

        # Read clusters data into a pandas DataFrame
        clusters = pd.read_csv(clusters_file_path)
        # Remove all text before (and including) the second to last underscore in the "chromosome" column of the clusters DataFrame
        clusters["chromosome"] = clusters["chromosome"].apply(lambda x: "_".join(x.split('_')[-2:]))
        clusters['start_cl'], clusters['stop_cl'] = clusters[['start_cl', 'stop_cl']].values.T.tolist()
        mask = clusters['stop_cl'] < clusters['start_cl']
        clusters.loc[mask, ['start_cl', 'stop_cl']] = clusters.loc[mask, ['stop_cl', 'start_cl']]

        print(clusters["chromosome"])

        # Create an empty DataFrame to store both clusters and singleton genes
        bed_data = []

        # Read GTF file as a TSV
        with open(gtf_files[0], 'r') as gtf_file:
            counter = 0
            for line in gtf_file:
                if line.startswith('#'):
                    continue  # Skip comment lines

                fields = line.strip().split('\t')
                chromosome_full = fields[0]
                parts = chromosome_full.split('_')
                chromosome = "_".join(parts[-2:])
                start = int(fields[3])
                stop = int(fields[4])

                # Check if the gene is contained in one of the clusters
                matching_clusters = clusters[(clusters['chromosome'] == chromosome) & (clusters['start_cl'] <= start) & (clusters['stop_cl'] >= stop)]

                # If no matching cluster is found, add it as a singleton gene
                if matching_clusters.empty:
                    bed_data.append((chromosome, start, stop))
                    counter += 1

        # Convert the bed_data list to a DataFrame
        bed_df = pd.DataFrame(bed_data, columns=["chromosome", "start", "stop"])

        # Save the clusters and singleton genes to a single BED file
        bed_df.to_csv(f"singletons.bed", sep='\t', header=False, index=False)
        print(counter)

        clusters.to_csv("clusters.bed", sep='\t', header=False, index=False, columns=["chromosome", "start_cl", "stop_cl"])
    except Exception as e:
        # Handle the exception (you can log the error, print a message, etc.)
        print(f"Error processing accession {accession}: {str(e)}")
        pass
# Read accession numbers from a file
with open("../../results/accessions_mini_run.txt", "r") as accession_file:
    for line in accession_file:
        accession = line.strip()
        process_accession(accession)
