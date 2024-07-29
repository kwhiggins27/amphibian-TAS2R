#!/usr/bin/env python

import pandas as pd


# Initialize an empty dictionary to store clusters
clusters = {}

# Read the file
file_path = '../../results/phylogenetics/clusters_output_with_names_v5_q2_withBH.txt'

with open(file_path, 'r') as file:
    data = file.read()

# Split data into clusters
cluster_texts = data.split('\n\n')

# Process clusters
for cluster_text in cluster_texts:
    if cluster_text:
        lines = cluster_text.split('\n')
        cluster_name = lines[0].strip()
        rows_text = [line.strip(" '[]") for line in lines if line.startswith('Rows')]
        if rows_text:
            rows = rows_text[0].split(': ')[1].replace('"', '').replace("'", '').replace("[", "").replace("]", "")
            rows = [row.strip() for row in rows.split(',')]
            truncated_rows = []
            for row in rows:
                second_underscore_index = row.find('_', row.find('_') + 1)
                if second_underscore_index != -1:
                    truncated_row = row[:second_underscore_index]
                    truncated_rows.append(truncated_row)
                else:
                    truncated_rows.append(row)
            clusters[cluster_name] = truncated_rows


# Read the CSV file into a DataFrame
converter_file_path = '../../results/coordinate_analysis/genes_clusters_etc3.csv'
converter = pd.read_csv(converter_file_path)

# Create an empty dictionary to store clusterB lists
clusterB_dict = {}

# Iterate through each cluster in clusters
for cluster_name, cluster_accessions in clusters.items():
    clusterB = []
    for accession in cluster_accessions:
        # Find matching rows in the converter DataFrame
        matching_row = converter[converter['Accession'] == accession]
        if not matching_row.empty:
            # Append the corresponding PlottingClade to the clusterB list
            clusterB.extend(matching_row['order'].tolist())
            #clusterB.extend(matching_row['latin.x'].tolist())

    # Add the clusterB list to the dictionary
    clusterB_dict[f'{cluster_name}_B'] = clusterB

# Print the generated clusterB lists
for cluster_name, clusterB_list in clusterB_dict.items():
    print(f"{cluster_name}: {clusterB_list}")

# Save the generated clusterB lists to a text file
output_file = '../../results/phylogenetics/clusters_output_with_names_v5_q2_withBH_order.txt'
with open(output_file, 'w') as file:
    for cluster_name, clusterB_list in clusterB_dict.items():
        file.write(f"{cluster_name}: {clusterB_list}\n")
