library(MASS)
library(ggplot2)
library(taxize)
library(readr)
library(dplyr)

setwd("/Volumes/wengpj01/vertebrate_pipeline/coordinates/")

# Set the file path to your CSV file
csv_file <- "position_within_chromosome_all_genes_added_big_1026.csv"

genes_accessions_input <- "/Volumes/wengpj01/vertebrate_pipeline/match_accession_order_and_class.csv"
# genes_accessions_input <- "/Volumes/wengpj01/vertebrate_pipeline/20231005_run/20231005_basics.csv"
genes_accessions <- read_csv(genes_accessions_input)


df_location <- read_csv(csv_file)
colnames(df_location) <- c("accession", "chromosome","start", "mini", "ends")

combo_df <- merge(genes_accessions, df_location, by.x = "Accession", by.y = "accession", all.x = TRUE)
write.csv(combo_df, file = "/Volumes/wengpj01/vertebrate_pipeline/20240115/end_summary2.csv", row.names = FALSE)

clusters <- read.csv("/Volumes/wengpj01/vertebrate_pipeline/search_neighborhood/maximal_range_all.csv")
# Sort clusters dataframe by accession, chromosome, and start
clusters <- clusters[order(clusters$accession, clusters$chromosome, clusters$start), ]

# Create a column 'cluster_name' in clusters dataframe using the counter for each accession
clusters$cluster_name <- with(clusters, ave(accession, accession, FUN = function(x) paste0(x, "_", seq_along(x))))

# Merge clusters with df_ordered_combo based on Accession and chromosome
merged_data <- merge(combo_df, clusters, by.x = c("Accession", "chromosome"), by.y = c("accession", "chromosome"), all.x = TRUE)

# Convert columns to integers
merged_data$start.x <- as.integer(merged_data$start.x)
merged_data$start.y <- as.integer(merged_data$start.y)
merged_data$stop.y <- as.integer(merged_data$stop)

# Filter and find matches based on conditions
filtered_data <- merged_data %>%
  filter(start.y < start.x, start.x < stop.y)
filtered_data <- unique(filtered_data)

duplicates_start_x <- filtered_data$start.x[duplicated(filtered_data$start.x)]

filtered_data <- filtered_data %>%
  group_by(cluster_name) %>%
  mutate(type = ifelse(n() == 1, "singleton", "cluster")) %>%
  ungroup()

write.csv(filtered_data, file = "/Volumes/wengpj01/vertebrate_pipeline/20240115/end_summary.csv", row.names = FALSE)


# Subset rows where type is singleton
singleton_data <- subset(filtered_data, type == "singleton")

# Calculate average and median values for ends in singleton data
singleton_avg <- mean(singleton_data$ends)
singleton_median <- median(singleton_data$ends)

# Subset rows where type is cluster
cluster_data <- subset(filtered_data, type == "cluster")

# Calculate average and median values for ends in cluster data
cluster_avg <- mean(cluster_data$ends)
cluster_median <- median(cluster_data$ends)

# Print the results
cat("Singleton Data:\n")
cat("Average ends:", singleton_avg, "\n")
cat("Median ends:", singleton_median, "\n\n")

cat("Cluster Data:\n")
cat("Average ends:", cluster_avg, "\n")
cat("Median ends:", cluster_median, "\n")

t_test <-t.test(singleton_data$ends, cluster_data$ends)

nrow(cluster_data)/(nrow(singleton_data) + nrow(cluster_data))

#unique_accessions <- unique(df_location$accession)
# 
# unique_latin <- c()
# 
# unique_latin[unique_latin == "Crotalus_viridis_viridis"] <- "Crotalus_viridis"

# Iterate through unique_accessions
# for (acc in unique_accessions) {
#   # Find the matching row in genes_accessions
#   match_row <- df_location[df_location$accession == acc, ]
#   
#   # Check if a matching row was found
#   if (nrow(match_row) > 0) {
#     # Add the latin value to the unique_latin list
#     unique_latin <- c(unique_latin, match_row$latin)
#   } else {
#     # If no matching row was found, add NA to unique_latin
#     unique_latin <- c(unique_latin, NA)
#   }
# }
# 
# get_taxonomic_info <- function(latin_name) {
#   # Get the taxonomic information for the Latin name from NCBI
#   tax_info <- tax_name(latin_name, get = c("class", "order", "family"), db = "ncbi")
#   
#   # Extract the order and family names (assuming they are single names)
#   class_name <- ifelse(length(tax_info$class) > 0, tax_info$class[[1]], NA)
#   order_name <- ifelse(length(tax_info$order) > 0, tax_info$order[[1]], NA)
#   family_name <- ifelse(length(tax_info$family) > 0, tax_info$family[[1]], NA)
#   
#   return(c(class = class_name, order = order_name, family = family_name))
# }
# 
# # Apply the function to the dataframe
# # taxonomic_info <- sapply(unique_latin, get_taxonomic_info)
# # taxonomic_info <- t(taxonomic_info)
# taxonomic_info <- read.csv("tax_info_unique.csv")
# 
# 
# write.csv(taxonomic_info, file = "tax_info_unique.csv", row.names = FALSE)
# taxonomic_info <- cbind(taxonomic_info, Accessions = unique_accessions)
# taxonomic_info <- cbind(taxonomic_info, Latin = unique_latin)

# taxonomic_info_list <- list()
# 
# # Iterate through unique_latin values
# for (latin_name in unique_latin) {
#   tryCatch({
#     # Get the taxonomic information for the Latin name from NCBI
#     tax_info <- tax_name(latin_name, get = c("class", "order", "family"), db = "ncbi")
#     
#     # Extract the order and family names (assuming they are single names)
#     class_name <- ifelse(length(tax_info$class) > 0, tax_info$class[[1]], NA)
#     order_name <- ifelse(length(tax_info$order) > 0, tax_info$order[[1]], NA)
#     family_name <- ifelse(length(tax_info$family) > 0, tax_info$family[[1]], NA)
#     
#     # Create a list of taxonomic information for the current Latin name
#     tax_info_list <- list(class = class_name, order = order_name, family = family_name)
#     
#     # Add the taxonomic information to the taxonomic_info_list
#     taxonomic_info_list <- c(taxonomic_info_list, list(tax_info_list))
#   }, error = function(e) {
#     # If an error occurs, add an empty element to the taxonomic_info_list
#     taxonomic_info_list <- c(taxonomic_info_list, list(list()))
#   })
# }
# taxonomic_info_df <- data.frame(do.call(rbind, taxonomic_info_list))
# subset_taxonomic_info_df <- taxonomic_info_df[1:637, ]
# rows_with_na <- which(rowSums(is.na(taxonomic_info_df)) > 0)
# na_rows_list <- as.list(na_rows)
# 
# latin_with_na <- unique_latin[na_rows]
# 
# taxonomic_info_list_pt2 <- list()
# 
# # Iterate through unique_latin values
# for (latin_name in latin_with_na) {
#   tryCatch({
#     # Get the taxonomic information for the Latin name from NCBI
#     tax_info <- tax_name(latin_name, get = c("class", "order", "family"), db = "ncbi")
#     
#     # Extract the order and family names (assuming they are single names)
#     class_name <- ifelse(length(tax_info$class) > 0, tax_info$class[[1]], NA)
#     order_name <- ifelse(length(tax_info$order) > 0, tax_info$order[[1]], NA)
#     family_name <- ifelse(length(tax_info$family) > 0, tax_info$family[[1]], NA)
#     
#     # Create a list of taxonomic information for the current Latin name
#     tax_info_list <- list(class = class_name, order = order_name, family = family_name)
#     
#     # Add the taxonomic information to the taxonomic_info_list
#     taxonomic_info_list <- c(taxonomic_info_list_pt2, list(tax_info_list))
#   }, error = function(e) {
#     # If an error occurs, add an empty element to the taxonomic_info_list
#     taxonomic_info_list_pt2 <- c(taxonomic_info_list_pt2, list(list()))
#   })
# }
# taxonomic_info_df_pt2 <- data.frame(do.call(rbind, taxonomic_info_list_pt2))
# 
# # Merge df_location and taxonomic_info based on the "accessions" column
# merged_df <- merge(df_location, taxonomic_info, by.x = "accession", by.y = "Accessions", all.x = TRUE)
# 
# # Rename the merged column to "class"
# merged_df <- rename(merged_df, class = class.y)






# Create a subset table where i_mini is TRUE
subset_df <- combo_df[combo_df$mini == TRUE, ]
data <- subset_df



# Create a histogram of i_ends
pdf("fraction_near_ends_all_genes_1021.pdf")
ggplot(subset_df, aes(x = ends)) +
  geom_histogram(binwidth = 0.01, fill = "blue", color = "black") +
  labs(title = "Histogram of distance from end",
       x = "ends",
       y = "Frequency") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()

unique_accessions <- unique(subset_df$Accession)
num_unique_accessions <- length(unique_accessions)

amphibia_subset <- subset_df[subset_df$PlottingClade %in% c("Anura", "Gymnophiona", "Caudata"), ]
unique_accessions <- unique(amphibia_subset$Accession)
num_unique_accessions <- length(unique_accessions)
frogs_subset <- subset_df[subset_df$PlottingClade %in% c("Anura"), ]
salamanders_subset <- subset_df[subset_df$PlottingClade %in% c("Caudata"), ]
amphibia_subset <- amphibia_subset[complete.cases(amphibia_subset), ]
not_amphibia_subset <- subset_df[!(subset_df$PlottingClade %in% c("Anura", "Gymnophiona", "Caudata")), ]
not_amphibia_subset <- not_amphibia_subset[complete.cases(not_amphibia_subset), ]

reptile <- subset_df[subset_df$PlottingClade %in% c("Squamata"),]
reptile_high_family <- c("Cordylidae", "Dactyloidae", "Eublepharidae", "Lacertidae", "Phrynosomatidae", "Rhineuridae", "Sphaerodactylidae")
high_reptile <- reptile[reptile$family %in% reptile_high_family]

t_test_ends <-t.test(amphibia_subset$ends, not_amphibia_subset$ends,alternative = "less")
sd(amphibia_subset$ends)
sd(not_amphibia_subset$ends)
ks_result <- ks.test(amphibia_subset$ends, not_amphibia_subset$ends)


pdf("fraction_near_ends_amphibians_1025.pdf")
ggplot(amphibia_subset, aes(x = ends)) +
  geom_histogram(binwidth = 0.01, fill = "blue", color = "black") +
  labs(title = "Histogram of distance from end for amphibians",
       x = "ends",
       y = "Frequency") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()

pdf("fraction_near_ends_frogs_1025.pdf")
ggplot(frogs_subset, aes(x = ends)) +
  geom_histogram(binwidth = 0.01, fill = "blue", color = "black") +
  labs(title = "Histogram of distance from end for amphibians",
       x = "ends",
       y = "Frequency") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()

pdf("fraction_near_ends_salamanders_1025.pdf")
ggplot(salamanders_subset, aes(x = ends)) +
  geom_histogram(binwidth = 0.01, fill = "blue", color = "black") +
  labs(title = "Histogram of distance from end for amphibians",
       x = "ends",
       y = "Frequency") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()

pdf("fraction_near_ends_notamphibians_1025.pdf")
ggplot(not_amphibia_subset, aes(x = ends)) +
  geom_histogram(binwidth = 0.01, fill = "blue", color = "black") +
  labs(title = "Histogram of distance from end for non-amphibians",
       x = "ends",
       y = "Frequency") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()

pdf("fraction_near_ends_highreptiles_1025.pdf")
ggplot(high_reptile, aes(x = ends)) +
  geom_histogram(binwidth = 0.01, fill = "blue", color = "black") +
  labs(title = "Histogram of distance from end for lizards and snakes",
       x = "ends",
       y = "Frequency") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()

median(not_amphibia_subset$ends)
mean(high_reptile$ends)

accession_counts <- subset_df %>%
  group_by(accession) %>%
  summarise(entry_count = n())

# Calculate the average value of 'ends' per accession
average_ends <- subset_df %>%
  group_by(accession) %>%
  summarise(average_ends = mean(ends))

# Perform a correlation test
correlation_test <- cor.test(accession_counts$entry_count, average_ends$average_ends)

p_value <- correlation_test$p.value

# Create a scatter plot
pdf("ends_vs_genes.pdf")
ggplot(data = merge(accession_counts, average_ends, by = "accession"), aes(y = entry_count, x = average_ends)) +
  geom_point() +
  labs(x = "Average Distance from Ends", y = "Average Number Genes") +
  ggtitle("Scatter Plot of Average Ends vs. Average Number of Genes")
dev.off()










# Create a histogram of i_ends
pdf("fraction_near_ends_all_genes_includingshort.pdf")
ggplot(df_location, aes(x = ends)) +
  geom_histogram(binwidth = 0.01, fill = "blue", color = "black") +
  labs(title = "Histogram of distance from end",
       x = "ends",
       y = "Frequency") +
  theme_minimal()
dev.off()

# Create bins (adjust the number and width of bins as needed)
num_bins <- 20
bins <- seq(min(data$ends), max(data$ends), length.out = num_bins + 1)

# Cut the data into bins and count frequencies
bin_counts <- cut(data$ends, breaks = bins, include.lowest = TRUE, right = TRUE)
observed_freq <- table(bin_counts)

# Calculate the expected frequency assuming a uniform distribution
total_data_points <- length(data$ends)
expected_freq <- rep(total_data_points / num_bins, num_bins)

# Calculate Chi-squared statistic
chi_squared <- sum((observed_freq - expected_freq)^2 / expected_freq)

# Calculate degrees of freedom
df <- num_bins - 1

# Calculate p-value
p_value <- 1 - pchisq(chi_squared, df)

# Set significance level
alpha <- 0.05

# Perform hypothesis test
if (p_value < alpha) {
  cat("Reject the null hypothesis: Data differs from a uniform distribution.\n")
} else {
  cat("Fail to reject the null hypothesis: Data does not significantly differ from a uniform distribution.\n")
}

# Print the Chi-squared statistic and p-value
cat("Chi-squared:", chi_squared, "\n")
cat("Degrees of Freedom:", df, "\n")
cat("P-value (Scientific Notation):", format(p_value, scientific = TRUE), "\n")


##########

# Calculate lambda (mean rate of gene occurrence)
##lambda <- length(gene_coordinates) / chromosome_length
#lambda <- mean(subset_df$ends)/0.5

## Generate a sequence of values for the Poisson distribution
#x <- seq(0, 0.5, length.out = 50)

## Calculate the Poisson probability mass function
#poisson_prob <- dpois(x, lambda)

## Create a histogram to visualize your data

histogram <- hist(subset_df$ends, breaks = 20, main = "Gene Distribution (Poisson Assumption)",
     xlab = "Normalized Gene Coordinates", ylab = "Frequency")

# Calculate the density within each bin
#bin_width <- diff(histogram$breaks)
#density_in_bins <- histogram$counts / (length(subset_df$ends) * bin_width)

# Extract the counts and midpoints from the histogram
counts <- histogram$counts
midpoints <- histogram$mids

# Define a function to calculate the negative log-likelihood for a Poisson distribution
poisson_log_likelihood <- function(lambda) {
  -sum(dpois(midpoints, lambda, log = TRUE) * counts)
}


# Initial guess for the Poisson parameter (you can change this if needed)
initial_lambda <- 0.1

# Fit the Poisson distribution using maximum likelihood estimation (MLE)
poisson_fit <- fitdistr(midpoints, poisson_log_likelihood, start = list(lambda = initial_lambda))

# Get the estimated Poisson parameter
lambda_estimate <- poisson_fit$estimate
## Find the maximum density
#max_density <- max(density_in_bins)

pdf("fraction_near_ends_all_genes_poisson.pdf")
hist(subset_df$ends, breaks = 20, main = "Gene Distribution (Poisson Assumption)",
     xlab = "Normalized Gene Coordinates", ylab = "Frequency")
# Overlay the Poisson distribution
lines(x, poisson_prob * max_density, col = "red")
dev.off()
