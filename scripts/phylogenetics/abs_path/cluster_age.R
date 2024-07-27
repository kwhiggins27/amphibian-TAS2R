library(ape)
library(phytools)
library(paleotree)
library(stringr)

setwd("/Volumes/wengpj01/vertebrate_pipeline/search_neighborhood")

# Read the file
file_path <- "/Volumes/wengpj01/vertebrate_pipeline/search_neighborhood/clusters_output_with_names_v7_q2_withBH.txt"
file_content <- readLines(file_path)

# Initialize variables
clusters <- list()
current_cluster <- NULL

# Function to extract accession number up to the second underscore
extract_accession <- function(string) {
  unlist(regmatches(string, regexpr("^([[:alnum:]_]+\\.[[:alnum:]]+)", string)))  # Extract the pattern
}

# Process each line in the file
for (line in file_content) {
  if (grepl("Cluster", line)) {
    # Extract cluster name
    current_cluster <- gsub("[: ]", "", line)
    current_accessions <- list()
    clusters[[current_cluster]] <- current_accessions
  } else if (grepl("Rows:", line)) {
    # Extract accession numbers from the 'Rows' section
    accession_numbers <- unlist(strsplit(gsub("Rows: \\['|'\\]", "", line), "', '"))
    accession_numbers <- lapply(accession_numbers, extract_accession)
    clusters[[current_cluster]] <- accession_numbers
  }
}

#Read in converter file
file_path <- "/Volumes/wengpj01/vertebrate_pipeline/20231021_run/genes_clusters_etc3.csv"  # Update this path if needed
converter <- read.csv(file_path)

# Initialize clusters_latin
clusters_latin <- list()

# Function to find Latin name from converter dataset
find_latin_name <- function(accession) {
  row <- converter[converter$Accession == accession, ]
  if (nrow(row) > 0) {
    return(row$latin.x)
  } else {
    return(NA)  # If no corresponding value is found
  }
}
find_PC <- function(accession) {
  row <- converter[converter$Accession == accession, ]
  if (nrow(row) > 0) {
    return(row$PlottingClade)
  } else {
    return(NA)  # If no corresponding value is found
  }
}

# Process each cluster
for (cluster_name in names(clusters)) {
  current_cluster <- clusters[[cluster_name]]
  accession_list <- list()
  
  # Process each accession in the cluster
  for (accession in current_cluster) {
    latin_name <- find_latin_name(accession)
    accession_list[[accession]] <- latin_name
  }
  
  clusters_latin[[cluster_name]] <- accession_list
}

clusters_PC <- list()

# Process each cluster
for (cluster_name in names(clusters)) {
  current_cluster <- clusters[[cluster_name]]
  accession_list <- list()
  
  # Process each accession in the cluster
  for (accession in current_cluster) {
    latin_name <- find_PC(accession)
    accession_list[[accession]] <- latin_name
  }
  
  clusters_PC[[cluster_name]] <- accession_list
}
# Read the file containing the list of Latin names
not_in_tree_file <- "/Volumes/wengpj01/phylogenetics/species_tree/not_in_tree.txt"
not_in_tree <- readLines(not_in_tree_file)
not_in_tree <- c(not_in_tree, "Coregonus_sp._'balchen'")

# Identify and remove Latin names present in not_in_tree from clusters_latin
final_latin <- clusters_latin

# Function to remove duplicates from a list while maintaining order
remove_duplicates <- function(lst) {
  unique_items <- unique(lst)
  unique_items[match(unique_items, lst)]
}

# Remove duplicates from each cluster in final_latin
final_latin_unique <- lapply(final_latin, remove_duplicates)

for (cluster_name in names(final_latin)) {
  taxa <- final_latin[[cluster_name]]
  taxa_list <- unlist(taxa)
  
  # Remove Latin names present in not_in_tree
  taxa_list <- taxa_list[!taxa_list %in% not_in_tree]
  
  # Update final_latin with modified taxa_list
  final_latin[[cluster_name]] <- list(taxa_list)
}


tr=read.tree("/Volumes/wengpj01/phylogenetics/species_tree/Timetree_ultrametric.tre")
# Initialize cluster_age list


cluster_age <- list()  # Initialize the list to store cluster ages

#for below, was final_latin_unique but changed to final_latin
for (cluster_name in names(final_latin)) {
  # Get the MRCA for the current cluster
  taxa <- final_latin[[cluster_name]]
  taxa_list <- unlist(taxa)
  
  tryCatch({
    mrca=getMRCA(tr,taxa_list)
    tmrca=dateNodes(tr)[mrca]
    cluster_age[[cluster_name]] <- tmrca
  }, error = function(e) {
    if (grepl("missing value where TRUE/FALSE needed", e$message)) {
      cluster_age[[cluster_name]] <- NA
      print(taxa_list)
    } else {
      stop("Unexpected error occurred:", e$message)
    }
  })
}
cluster_age
#getMRCA(tr,c("Agelaius_phoeniceus", "Ammodramus_caudacutus", "Diglossa_brunneiventris", "Junco_hyemalis", "Camarhynchus_parvulus"))

# Extract ages from cluster_age
cluster_ages <- unlist(cluster_age, use.names = FALSE)

# Create a histogram
pdf("scratch2/cluster_age2.pdf")
hist(cluster_ages, main = "Age Distribution of Clusters", xlab = "Minimum Age (Million Years)", ylab = "Frequency", col = "#166A53", border = "black", breaks = seq(0, max(cluster_ages) + 50, by = 50))
dev.off()

clusters_with_Anura <- names(Filter(function(x) any("Anura" %in% unlist(x)), clusters_PC))

# Filter cluster_age for clusters containing "Anura"
cluster_age_Anura <- cluster_age[clusters_with_Anura]

# Extract non-NA values for the histogram
values_to_plot <- unlist(cluster_age_Anura, use.names = FALSE)
values_to_plot <- values_to_plot[!is.na(values_to_plot)]

pdf("cluster_age_anura.pdf")
hist(values_to_plot, breaks = seq(0, max(values_to_plot) + 50, by = 50),
     main = "Histogram of Cluster Ages with 'Anura'", xlab = "Cluster Age")
dev.off()


# Initialize empty vectors to store data
cluster_names <- character()
ages <- numeric()
item_counts <- numeric()

# # Loop through cluster_age and num_items to match data
# for (cluster_name in names(cluster_age)) {
#   # Check if cluster_age has a non-zero value and num_items is available for the cluster
#   if (!is.null(cluster_age[[cluster_name]]) && length(cluster_age[[cluster_name]]) > 0 &&
#       !is.na(num_items[cluster_name])) {
#     # Append the cluster name, age, and number of items to the respective vectors
#     cluster_names <- c(cluster_names, rep(cluster_name, length(cluster_age[[cluster_name]])))
#     ages <- c(ages, unlist(cluster_age[[cluster_name]]))
#     item_counts <- c(item_counts, rep(num_items[cluster_name], length(cluster_age[[cluster_name]])))
#   }
# }
# 
# # Create a data frame from the vectors
# data <- data.frame(Cluster = cluster_names, Age = ages, Num_Items = item_counts)
# 
# # Remove rows where Age is NA or 0 and Num_Items is 0
# data <- data[complete.cases(data$Age, data$Num_Items) & data$Num_Items != 0 & data$Age != 0, ]



# # Plot
# pdf("cluster_age_size.pdf")
# plot(data$Age, data$Num_Items, xlab = "Cluster Age", ylab = "Number of Items in the Original List", 
#      main = "Cluster Age vs Number of Items in Original List")
# dev.off()

age_filter <- (cluster_age)
filtered_cluster_names <- names(cluster_age)[!is.na(age_filter) & age_filter >= 50 & age_filter <= 100]

# Exclude numeric(0) values from cluster_age
filtered_cluster_age <- cluster_age[filtered_cluster_names]
filtered_cluster_age <- filtered_cluster_age[lengths(filtered_cluster_age) > 0]

# Extract values from clusters_PC for the filtered clusters
filtered_clusters_PC <- clusters_PC[filtered_cluster_names]
# Exclude null values
filtered_clusters_PC <- filtered_clusters_PC[sapply(filtered_clusters_PC, function(x) !all(is.na((x))))]
######

# Generating all possible pairs of tip labels
all_pairs <- combn(tr$tip.label, 2, simplify = TRUE)

# Randomly select 1000 pairs from all_pairs
set.seed(123)  # Setting seed for reproducibility
sampled_indices <- sample(ncol(all_pairs), 1000)
sampled_pairs <- all_pairs[, sampled_indices]
# Initialize a list to store the results
sim_node_age <- list()

# Loop through each pair and compute MRCA and TMRCA
for (i in 1:ncol(sampled_pairs)) {
  pair <- sampled_pairs[, i]
  mrca <- getMRCA(tr, pair)
  tmrca <- dateNodes(tr)[mrca]
  sim_node_age[[i]] <- tmrca
}


# Plotting the histogram of node ages
pdf("scratch2/sim_node_ages.pdf")
hist(unlist(sim_node_age), main = "Simulation: pairwise node ages", xlab = "Node Ages", ylab = "Frequency")
dev.off()
#######

# Assuming 'tr' is your phylogenetic tree
node_ages <- branching.times(tr)
max_branch <-max(node_ages)

pdf("scratch2/node_ages.pdf")
hist(node_ages, main = "Node ages", xlab = "Node Ages", ylab = "Frequency")
dev.off()

# Initialize an empty data frame
cluster_table <- data.frame(Cluster = character(), Full_Names = character(), Age = numeric(), stringsAsFactors = FALSE)

# Process each line in the file
for (line in file_content) {
  if (grepl("Cluster", line)) {
    # Extract cluster name
    current_cluster <- gsub("[: ]", "", line)
    current_rows <- NULL
  } else if (grepl("^Rows:", line)) {
    # Extract the entire string from the line starting with 'Rows:'
    current_rows <- gsub("^Rows: ", "", line)
    
    # Check if 'Rows' section is not empty
    if (length(trimws(current_rows)) > 0) {
      # Count the number of occurrences of "GCA" in Full_Names
      num_occurrences <- str_count(current_rows, "GCA")
      
      # Get Latin names directly from final_latin
      latin_names <- final_latin[[current_cluster]]
      
      # Create a data frame for the current cluster
      cluster_data <- data.frame(Cluster = current_cluster, Full_Names = current_rows, Latin_Names = paste(latin_names, collapse = ", "), Age = ifelse(current_cluster %in% names(cluster_age), cluster_age[[current_cluster]], ""), Num_occurrences = num_occurrences, stringsAsFactors = FALSE)
      
      # Append the data frame to the cluster_table
      cluster_table <- rbind(cluster_table, cluster_data)
    }
  }
}

# Save cluster_table as a CSV file
write.csv(cluster_table, file = "/Users/kwhiggins27/Desktop/OneDrive - Harvard University/Weng Lab/Amphibian_paper/Figures/cluster_table.csv", row.names = FALSE)

# Plot
pdf("cluster_age.pdf")
plot(cluster_table$Age, data$Num_Items, xlab = "Cluster Age", ylab = "Number of Items in the Original List",
     main = "Cluster Age vs Number of Items in Original List")
dev.off()

