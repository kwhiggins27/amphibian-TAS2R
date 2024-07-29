#!/usr/bin/env python
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
plt.rcParams['pdf.fonttype'] = 42

# Read the first dataframe from 20231006_basics.csv
my_work = pd.read_csv('../../results/coordinate_analysis/20231006_basics.csv', header=0)

# Read the second dataframe from genome_size_data_011123_11_59_51.csv
genome_database = pd.read_csv('../../results/coordinate_analysis/genome_size_data_011123_11_59_51.csv', header=0)

# Replace underscores with spaces in the "latin" column of my_work
my_work['latin'] = my_work['latin'].str.replace('_', ' ')

# Create a dictionary to store the average C-values for each unique "latin" value
average_c_values = {}

# Iterate through each row in genome_database and calculate the mean C-value for each unique "Species"
for _, row in genome_database.iterrows():
    species = row['Species']
    c_value = row['C-value']

    # Check if the C-value is a valid number (float or int)
    if pd.notna(c_value) and str(c_value).replace('.', '', 1).isdigit():
        c_value = float(c_value)
        if species in average_c_values:
            average_c_values[species].append(c_value)
        else:
            average_c_values[species] = [c_value]

# Calculate the average C-value for each unique "latin" value and store it in a new list
avg_c_values = []
for _, row in my_work.iterrows():
    species = row['latin']
    if species in average_c_values:
        c_values = average_c_values[species]
        avg_c_value = sum(c_values) / len(c_values)
        avg_c_values.append(avg_c_value)
    else:
        avg_c_values.append(None)  # Set to None if no matching species found

# Add the "genome_size" column to my_work
my_work['genome_size'] = avg_c_values

# Manually specify the entries
manual_entries = {
    "Bombina bombina": 11.385,
    "Hymenochirus boettgeri": 4.013,
    "Xenopus borealis": 3.56,
    "Pyxicephalus adspersus": 1.4
}

# Update the "genome_size" column in my_work with the manual entries
my_work.loc[my_work['latin'].isin(manual_entries.keys()), 'genome_size'] = my_work['latin'].map(manual_entries)

# Create a new dataframe "combo" with rows where "genome_size" is not 0 or NaN
combo = my_work[my_work['genome_size'].notna() & (my_work['genome_size'] != 0)]

amphibian_accessions = [
    "GCA_000004195.4", "GCA_002915635.3", "GCA_004786255.1", "GCA_009667805.1",
    "GCA_014858855.1", "GCA_017654675.1", "GCA_018994145.1", "GCA_019447015.1",
    "GCA_019512145.1", "GCA_019857665.1", "GCA_024363595.1", "GCA_026652325.1",
    "GCA_027358695.2", "GCA_027579735.1", "GCA_027789765.1", "GCA_027917425.1",
    "GCA_028564925.1", "GCA_029206835.1", "GCA_029499605.1", "GCA_901001135.2",
    "GCA_901765095.2", "GCA_902459505.2", "GCA_905171765.1", "GCA_905171775.1"
]

# Add the "Amphibian" column based on the list of accession values
combo['Amphibian'] = combo['Accession'].isin(amphibian_accessions)

# Convert boolean values to "TRUE" and "FALSE"
combo['Amphibian'] = combo['Amphibian'].map({True: 'TRUE', False: 'FALSE'})


# Create a scatter plot with different colors for "Amphibian" values
plt.scatter(
    combo[combo['Amphibian'] == 'TRUE']['genome_size'],
    combo[combo['Amphibian'] == 'TRUE']['Number of Genes'],
    label='Amphibian', c='green'
)
plt.scatter(
    combo[combo['Amphibian'] == 'FALSE']['genome_size'],
    combo[combo['Amphibian'] == 'FALSE']['Number of Genes'],
    label='Non-Amphibian', c='blue'
)

plt.xlabel('Genome Size (C-value)')
plt.ylabel('Number of TAS2Rs')
plt.title('Genome Size vs Number of Genes')

# Add a legend
plt.legend()

# Save the plot as a PDF
plt.savefig('../../results/coordinate_analysis/genome_vs_TAS2Rcount2.pdf')


# Calculate the correlation coefficient and p-value for all rows where Amphibian is TRUE
true_amphibian = combo[combo['Amphibian'] == 'TRUE']
corr_true, p_value_true = pearsonr(true_amphibian['genome_size'], true_amphibian['Number of Genes'])

# Calculate the correlation coefficient and p-value for all rows where Amphibian is FALSE
false_amphibian = combo[combo['Amphibian'] == 'FALSE']
corr_false, p_value_false = pearsonr(false_amphibian['genome_size'], false_amphibian['Number of Genes'])

# Calculate the correlation coefficient and p-value for all rows
corr_all, p_value_all = pearsonr(combo['genome_size'], combo['Number of Genes'])

# Print the results
print("Correlation coefficient (Amphibian TRUE):", corr_true)
print("P-value (Amphibian TRUE):", p_value_true)

print("Correlation coefficient (Amphibian FALSE):", corr_false)
print("P-value (Amphibian FALSE):", p_value_false)

print("Correlation coefficient (All rows):", corr_all)
print("P-value (All rows):", p_value_all)
