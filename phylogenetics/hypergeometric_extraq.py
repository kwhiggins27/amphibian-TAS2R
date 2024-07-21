#!/usr/bin/env python
#SBATCH --job-name=hyperG # Job name
#SBATCH --mail-type=FAIL,END    # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=youremailaddress@yourinstitute    # Where to send mail
#SBATCH --mem=20gb          # Job memory request, down from 200 and closer to the 64gb I think you're using per instance
#SBATCH --nodes=1           # ensure cores are on one node
#SBATCH --ntasks=1          # run a single task
#SBATCH --cpus-per-task=1       # number of cores/threads requested, up from 4 you're asking for now - see too that I changed --runThreadN below to match
#SBATCH --partition=20
#SBATCH --output=logs/hyperG_%j.log # Standard output and error log
#SBATCH --error=logs/hyperG_%j.err

import pandas as pd
from scipy.cluster import hierarchy
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import hypergeom

#by accession
# Read the CSV file without specifying dtype (default is object)
calc_m = pd.read_csv('/lab/wengpj01/vertebrate_pipeline/search_neighborhood/calc_m7.csv', index_col=0)
calc_m = calc_m.apply(pd.to_numeric, errors='coerce')
calc_m = calc_m.fillna(0)

calc_n = pd.read_csv('/lab/wengpj01/vertebrate_pipeline/search_neighborhood/calc_n7.csv',index_col=0)
calc_n = calc_n.apply(pd.to_numeric, errors='coerce')
calc_n = calc_n.fillna(0)
#by cluster
calc_q = pd.read_csv('/lab/wengpj01/vertebrate_pipeline/search_neighborhood/calc_q7.csv',index_col=0)
calc_q = calc_q.apply(pd.to_numeric, errors='coerce')
calc_q = calc_q.fillna(0)
calc_k = pd.read_csv('/lab/wengpj01/vertebrate_pipeline/search_neighborhood/calc_k7.csv', index_col=0)
calc_k = calc_k.apply(pd.to_numeric, errors='coerce')
calc_k = calc_k.fillna(0)



accepted_range = pd.read_csv('/lab/wengpj01/vertebrate_pipeline/search_neighborhood/maximal_range_all.tsv', sep='\t')  #This df is notable for columns chromosome, start, stop, accession
accepted_range['id'] = (
    accepted_range['accession'] + "_" +
    accepted_range['chromosome'].astype(str) + ":" +
    accepted_range['start'].astype(str) + "_" +
    accepted_range['stop'].astype(str)
)

hg_p = pd.DataFrame(index=accepted_range["id"], columns=accepted_range["id"])



# Iterating through pairs of accepted_range_A and accepted_range_B
for index_A, row_A in accepted_range.iterrows():
    chromosome_A = str(row_A['chromosome'])
    start_A = row_A['start']
    stop_A = row_A['stop']
    accession_A = row_A['accession'].split(':')[0]  # Extracting accession for species A

    accepted_range2 = accepted_range[index_A:]

    for index_B, row_B in accepted_range2.iterrows():
        chromosome_B = str(row_B['chromosome'])
        start_B = row_B['start']
        stop_B = row_B['stop']
        accession_B = row_B['accession'].split(':')[0]  # Extracting accession for species B
        # calc_m.reset_index(drop=True, inplace=True)
        # calc_n.reset_index(drop=True, inplace=True)
        calc_m.reset_index(drop=True, inplace=True)
        calc_n.reset_index(drop=True, inplace=True)
        calc_q.reset_index(drop=True, inplace=True)
        calc_k.reset_index(drop=True, inplace=True)

        try:
            m_value = calc_m.iloc[index_A, index_B]
            n_value = calc_n.iloc[index_A, index_B]
            M_value = m_value + n_value
            q_value = calc_q.iloc[index_A, index_B]
            k_value = calc_k.iloc[index_A, index_B]
            print("start here")
            if M_value == 0:
                p_value = 1
                hg_p.iloc[index_A, index_B] = p_value
                print("M_value is 0")
            elif q_value < 2:
                p_value = 1
                hg_p.iloc[index_A, index_B] = p_value
                print("q_value < 2")
            else:
                p_value = 1 - hypergeom.cdf(q_value - 1, M_value, m_value, k_value)
                hg_p.iloc[index_A, index_B] = p_value

        except KeyError:
            print("problem!")
            # Skip this iteration if the accession number is not present in calc_m
            # print(f"Accession number {accession_A} or {accession_B} not found in calc_m, skipping this pair.")
            continue

hg_p.to_csv('/lab/wengpj01/vertebrate_pipeline/search_neighborhood/hg_p_value_v7_q2.csv')

# # Convert hg_p to a numeric DataFrame (it may contain strings like 'NaN' due to missing values)
# hg_p_numeric = hg_p.apply(pd.to_numeric, errors='coerce')
#
# # Perform hierarchical clustering on rows and columns
# row_linkage = hierarchy.linkage(hg_p_numeric.values, method='average')
# col_linkage = hierarchy.linkage(hg_p_numeric.values.T, method='average')
#
# # Obtain the reordered indices from row dendrogram
# row_dendrogram = hierarchy.dendrogram(row_linkage, no_plot=True)
# reordered_indices = row_dendrogram['leaves']
#
# # Reorder the rows of hg_p_numeric based on hierarchical clustering
# hg_p_reordered = hg_p_numeric.iloc[reordered_indices]
#
# # Save hg_p_reordered to a new file (replace 'your_path_here' with your desired file path)
# hg_p_reordered.to_csv('/lab/wengpj01/vertebrate_pipeline/search_neighborhood/hg_p_reordered.csv', index=True, header=True)
