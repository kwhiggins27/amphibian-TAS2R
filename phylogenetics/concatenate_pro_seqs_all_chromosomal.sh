#!/bin/bash

# Set the path to the directory containing subdirectories
directory="/lab/wengpj01/vertebrate_pipeline/subdirs"

# Create or overwrite the output files
output_file="/lab/wengpj01/vertebrate_pipeline/all_pro_chromosomal_20231005.txt"
problem_file="/lab/wengpj01/vertebrate_pipeline/pro_problems_20231005.txt"
touch "$output_file" "$problem_file"

# Path to the file containing the list of accessions to keep
accessions_to_keep_file="/lab/wengpj01/vertebrate_pipeline/accessions_to_keep.txt"

# Loop through the list of accessions to keep
while IFS= read -r accession; do
    subdir="$directory/$accession"

    if [ -d "$subdir" ]; then
        # Look for files matching the pattern '*_bitter_protein_confirmed2xv1.fasta' in the subdirectory
        confirmed_na_file=$(find "$subdir" -type f -name '*_bitter_protein_confirmed2xv1.fasta' -print -quit)

        if [ -n "$confirmed_na_file" ]; then
            # Get the subdirectory shortname
            subdir_shortname=$(basename "$subdir")

            # Prepend the shortname to the header lines and append the contents to all_chromosomal.txt
            sed "s/^>/>${subdir_shortname}_/" "$confirmed_na_file" >> "$output_file"
        else
            # If the file is not found, append the subdirectory name to problems.txt
            echo "$(basename "$subdir")" >> "$problem_file"
        fi
    fi
done < "$accessions_to_keep_file"


# # Set the path to the directory containing subdirectories
# directory="/lab/wengpj01/vertebrate_pipeline/subdirs"
#
# # Create or overwrite the output files
# touch /lab/wengpj01/vertebrate_pipeline/all_pro_chromosomal_20231005.txt /lab/wengpj01/vertebrate_pipeline/pro_problems_20231005.txt
#
# # Loop through subdirectories
# for subdir in "$directory"/*; do
#     if [ -d "$subdir" ]; then
#         # Look for files matching the pattern '*_bitter_protein_confirmed2xv1.fasta' in the subdirectory
#         confirmed_na_file=$(find "$subdir" -type f -name '*_bitter_protein_confirmed2xv1.fasta' -print -quit)
#
#         if [ -n "$confirmed_na_file" ]; then
#             # Get the subdirectory shortname
#             subdir_shortname=$(basename "$subdir")
#
#             # Prepend the shortname to the header lines and append the contents to all_chromosomal.txt
#             sed "s/^>/>${subdir_shortname}_/" "$confirmed_na_file" >>  /lab/wengpj01/vertebrate_pipeline/all_pro_chromosomal_20231005.txt
#         else
#             # If the file is not found, append the subdirectory name to problems.txt
#             echo "$(basename "$subdir")" >> /lab/wengpj01/vertebrate_pipeline/pro_problems_20231005.txt
#         fi
#     fi
# done
