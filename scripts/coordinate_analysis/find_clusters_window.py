#!/usr/bin/env python
import pandas as pd
import glob
import sys
import csv

accession=sys.argv[1]

# Use glob to match files with the specified wildcard
gtf_files = glob.glob(f'../../subdirs/{accession}/*_from_pipeline.gtf')

# Define the sorting function as shown in the previous response
def gtf_sort_key(line):
    fields = line.strip().split('\t')
    seqname = fields[0]
    start = int(fields[3])
    return (seqname, start)

for gtf_file in gtf_files:
    # Create the sorted output file name
    output_file = gtf_file.replace('_from_pipeline.gtf', '_sorted.gtf')

    # Read the GTF file into a list of lines
    with open(gtf_file, 'r') as gtf_file:
        gtf_lines = gtf_file.readlines()

    # Sort the GTF lines using the custom sorting key
    sorted_gtf_lines = sorted(gtf_lines[1:], key=gtf_sort_key)

    # Write the sorted GTF lines to the output file
    with open(output_file, 'w') as sorted_file:
        # Write the header line first
        sorted_file.write(gtf_lines[0])

        # Write the sorted GTF lines
        for line in sorted_gtf_lines:
            sorted_file.write(line)

    # Create a DataFrame from the sorted GTF lines
    gtf_df = pd.read_csv(output_file, sep='\t', header=None)

    # # Append the DataFrame to the list
    # gtf_files1.append(gtf_df)

# Initialize an empty DataFrame to store cluster information
clusters_found = pd.DataFrame(columns=['cluster number', 'start_cl', 'stop_cl', 'chromosome', 'num_cl'])

#Determine cluster size
max_skip=int(sys.argv[2])
max_skip_str=str(sys.argv[2])


# Function to check if two genes are part of the same cluster
def is_same_cluster(row1, row2):
    return row1['seqname'] == row2['seqname'] and abs(row1['start'] - row2['start']) <= max_skip

# Use glob to match files with the specified wildcard
gtf_files = glob.glob(f'../../subdirs/{accession}/*_sorted.gtf')


# Iterate through the GTF files
for gtf_file in gtf_files:
    # Read the GTF file with tab and semi-colon delimiters
    input = pd.read_csv(gtf_file, sep='\t', header=None, engine='python')
    input.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']

    # input = gtf_file
    # Initialize variables to track cluster information for this file
    current_cluster = None
    start_cl = None
    stop_cl = None
    cluster_num = 0
    num_cl = 0

# Iterate through the rows of the GTF DataFrame
    for i in range(len(input)):
        row = input.iloc[i]

        # If no cluster is started yet, check if this row can start a new cluster
        if current_cluster is None:
            for j in range(i + 1, len(input)):
                next_row = input.iloc[j]
                if is_same_cluster(row, next_row):
                    current_cluster = cluster_num
                    start_cl = row['start']
                    num_cl = 0
                    break

            # If a cluster has already started but is incomplete, check if the next gene belongs to the same cluster
        if current_cluster is not None:
            if i == len(input) - 1 or not is_same_cluster(row, input.iloc[i + 1]):
                stop_cl = row['end']
                current_cluster = None
                cluster_num += 1
                num_cl += 1
                # Get the 'seqname' of the first entry in the cluster for the 'chromosome' column
                chromosome = input.iloc[i]['seqname']
                clusters_found = pd.concat([clusters_found, pd.DataFrame([{'cluster number': cluster_num, 'start_cl': start_cl, 'stop_cl': stop_cl, 'num_cl': num_cl, 'chromosome': chromosome}])], ignore_index=True)
            elif is_same_cluster(row, input.iloc[i + 1]):
                num_cl += 1

clusters_found['size'] = abs(clusters_found['stop_cl'] - clusters_found['start_cl'])
clusters_found['spacing'] = (clusters_found['size'] / clusters_found['num_cl']).apply(lambda x: round(x, 2))

# Save the clusters_found DataFrame to a CSV file
# clusters_found.to_csv('/lab/wengpj01/vertebrate_pipeline/subdirs/GCA_028390025.1/clusters_found.csv', index=False)
clusters_found.to_csv(f'../../subdirs/{accession}/clusters_found_{max_skip}.csv', index=False)

#fraction of genes in clusters of 2 or more
total_genes = len(input)  # Total number of genes in the GTF file
clusters_with_2_or_more_genes = clusters_found[clusters_found['num_cl'] >= 2]
genes_in_clusters_of_2_or_more = clusters_with_2_or_more_genes['num_cl'].sum()

fraction_genes_in_clusters_of_2_or_more = round(genes_in_clusters_of_2_or_more / total_genes, 2) if total_genes != 0 else 'NA'

# Initialize variables to track the cluster with the largest num_cl
largest_cluster_num = None
largest_num_cl = 0

# Iterate through the clusters_found DataFrame
for index, row in clusters_found.iterrows():
    if row['num_cl'] > largest_num_cl:
        largest_num_cl = row['num_cl']
        largest_cluster_num = row['cluster number']

# Find the spacing of the largest cluster
if largest_cluster_num is not None:
    largest_cluster = clusters_found[clusters_found['cluster number'] == largest_cluster_num]
    genes_in_largest_cluster = int(largest_cluster['num_cl'])
    spacing_big = round(largest_cluster['spacing'].values[0] / 1000, 2) if largest_cluster['num_cl'].values[0] != 0 else 'NA'
else:
    spacing_big = None
    genes_in_largest_cluster = 0


# Calculate total_len (sum of all values of size)
total_len = round(clusters_found['size'].sum()/1000,2)

# Calculate total_num (sum of all values of num_cl)
total_num = clusters_found['num_cl'].sum()

# Calculate total_space (total_len divided by total_num)
total_space = round(total_len / total_num, 2) if total_num != 0 else 'NA'

# Calculate the fraction of genes in the largest cluster
fraction_genes_in_largest_cluster = round(genes_in_largest_cluster / total_genes, 2) if total_genes != 0 else 'NA'



# Open the CSV file for writing
# csv_filename = '/lab/wengpj01/vertebrate_pipeline/summary_clusters.csv'  # Change the filename as needed
csv_filename = glob.glob(f'../../results/coordinate_analysis/summary_clusters_{max_skip_str}.csv')
print(max_skip_str)
print(csv_filename)
with open(csv_filename[0], 'a', newline='') as csv_file:
    writer = csv.writer(csv_file)

    # Write the header row if the file is empty
    if csv_file.tell() == 0:
        writer.writerow(['Accession', 'Number of Genes', 'Number of Clusters', 'Fraction Clustered','Fraction in Biggest Cluster', 'Spacing of Biggest Cluster(kb)', 'Total Len(kb)', 'Total Space (kb/receptor)'])

    # Write the data row
    writer.writerow([accession, total_genes, len(clusters_found), fraction_genes_in_clusters_of_2_or_more,fraction_genes_in_largest_cluster, spacing_big, total_len, total_space])

# Print a message to confirm the write operation
print(f'Data written to {csv_filename}')
