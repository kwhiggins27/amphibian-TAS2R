#!/bin/bash
#SBATCH --job-name=split3  # Job name
#SBATCH --mail-type=END,FAIL     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=youremailaddress@yourinstitute   # Where to send mail
#SBATCH --mem=200gb          # Job memory request, down from 200 and closer to the 64gb I think you're using per instance
#SBATCH --nodes=1           # ensure cores are on one node
#SBATCH --ntasks=1          # run a single task
#SBATCH --cpus-per-task=1       # number of cores/threads requested, up from 4 you're asking for now - see too that I changed --runThreadN below to match
#SBATCH --partition=20
#SBATCH --output=split3_%j.log # Standard output and error log
#SBATCH --error=split3_%j.err

shopt -s extglob

cd $1

source config_unix_species.sh

yes="yes"

# Create a directory for subdivided FASTA files
accession=$2
species=$BLAST_species
subdivided="${genome_location%/*}/${accession}_subdivided"
mkdir -p "$subdivided"
subdivided_test=$subdivided"/checkdone.txt"

if [ "$big_genome" == "$yes" ]
then
    if [ -f "$subdivided_test" ]
    # if (( "$subdivided" > 500000000)); then
      then
        echo "Directory already made"
    fi
else
  echo "Not a big genome"
fi

cd $subdivided
save_originals=$subdivided/originals

mkdir -p $save_originals


# Set the input file path
#input_file="/lab/wengpj01/non-amphibians/lungfish/lungfish_genome.fasta"
input_file=$genome_location

# Set the output directory path
#output_dir="/lab/wengpj01/non-amphibians/lungfish/subdivided2"
output_dir=$subdivided


# Split the large FASTA file into separate entries
csplit -s -f "${output_dir}/entry_" -b '%02d.fasta' "${input_file}" '/^>/' '{*}'

# Rename the output files based on the entry's header row
for file in ${output_dir}/*; do
  header=$(head -n 1 "$file" | sed 's/^>//;s/ .*//')
  if [[ "$file" != "${output_dir}/${header}.fasta" ]]; then
    mv "$file" "${output_dir}/${header}.fasta"
fi
done

# Check the size of each file in the subdivided directory
file_list=(${output_dir}/*.fasta)
file_count=${#file_list[@]}

# Function to combine neighboring files into a hybrid file
combine_files() {
  local combined_file="${output_dir}/${species}_$1.fasta"
  local temp_file="${output_dir}/temp.fasta"
  local file_size=0
  local i=$2

  # Combine neighboring files until the total size exceeds 100MB
  while [ $file_size -lt 100000000 ] && [ $i -lt $file_count ]; do
    file="${file_list[$i]}"
    file_size=$(stat -c %s "$file")
    cat "$file" >> "$temp_file"
    rm "$file"
    i=$((i + 1))
  done

  # Only perform the mv command if the combined file is different
  if ! cmp -s "$temp_file" "$combined_file"; then
    mv "$temp_file" "$combined_file"
  else
    rm "$temp_file"
  fi
}

# Process the files in the subdivided directory
file_index=0
species_index=1

while [ $file_index -lt $file_count ]; do
  file="${file_list[$file_index]}"
  file_size=$(stat -c %s "$file")

  if [ $file_size -gt 100000000 ]; then
    # If the file size is >100MB, leave it as is
    species_file="${output_dir}/${species}_${species_index}.fasta"
    mv "$file" "$species_file"
    species_index=$((species_index + 1))
    file_index=$((file_index + 1))
  else
    # If the file size is <100MB, combine with neighboring files
    combine_files "$species_index" "$file_index"
    species_index=$((species_index + 1))
    file_index=$((file_index + 1))
  fi
done



# # Function to check file size and subdivide if necessary
# subdivide_files() {
#   local file_index=$1
#
#   while [ $file_index -lt $file_count ]; do
#     file="${file_list[$file_index]}"
#     file_size=$(stat -c %s "$file" 2>/dev/null || echo 0)
#
#     if [ $file_size -gt 500000000 ]; then
#       # If file size is greater than 500MB, split into two files
#       name=`basename -s .fasta $file`
#       name2=$name\ab.fasta
#       split $name.fasta -n 2 --additional-suffix=.fasta $name  #"$file" "${file}."
#       sed -i "1s/^/\>$name\n/" "$name2"
#       rm "$file"
#       file_index=$((file_index + 1))
#     else
#       # File size is below 500MB, move to the next file
#       break
#     fi
#   done
# }
#
# # Subdivide files until all are below 500MB
# start_index=0
#
# while [ $start_index -lt $file_count ]; do
#   subdivide_files "$start_index"
#   start_index=$((start_index + 1))
# done

# Iterate over the fasta files in the directory
for file in "$subdivided"/*.fasta; do
  # Get the file size in bytes
  file_size=$(stat -c %s "$file")

  # Extract the file name without extension
  filename=$(basename "$file" .fasta)

  # Split the file into multiple parts
  split -b 500M "$file" "${filename}_"

  # Rename the output files and add a new header row
  suffixes=({a..z})
  suffix_index=0

  for output_file in "${filename}_"*; do
    # Get the suffix from the output file name
    suffix="${output_file#${filename}_}"

    # Create a temporary file to store the modified contents
    temp_file="${output_file}.tmp"

    # Check if the output file has a header row in the first line
    first_line=$(head -n 1 "$output_file")
    if [[ "$first_line" =~ ^\>.+ ]]; then
      # Delete the existing header row
      sed -i '1d' "$output_file"
    fi

    # Generate the new header row
    new_header=">${filename}_${suffix}"

    # Add the new header row at the start of the output file
    echo "$new_header" | cat - "$output_file" > "$temp_file" && mv "$temp_file" "$output_file"

    suffix_index=$((suffix_index + 1))
    if [ "$suffix_index" -ge "${#suffixes[@]}" ]; then
      suffix_index=0
    fi

    mv ${filename}_${suffix} ${filename}_${suffix}.fasta

    seqtk seq -l 50 "$file" > "${file}.tmp"
    mv "${file}.tmp" "$file"
  done
      echo "Splitting and reheadering complete."
      mv $file $save_originals
  done

# # Adjust line length to 50 nucleic acids (except the last line)
# for file in "${file_list[@]}"; do
#   seqtk seq -l 50 "$file" > "${file}.tmp"
#   mv "${file}.tmp" "$file"
# done

# Delete files of size 0 in the subdivided directory
#find "${output_dir}" -type f -size 0 -delete
find "${output_dir}" -type f -size -50c -delete

for file in *.fasta; do
    seqtk seq -l 50 "$file" > "${file%.fasta}.tmp"
    mv "${file%.fasta}.tmp" "$file"
done

echo "done with split" > checkdone.txt

#Save the contig names to a file named "names.txt"
rm -f names.txt

find . -maxdepth 1 -type f -name "*.fasta" -exec basename {} .fasta \; > names.txt
