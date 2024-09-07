library(MASS)
library(ggplot2)
library(taxize)
library(readr)
library(dplyr)

## Absolute path to working directory, where figures and files will be produced
setwd("/Volumes/wengpj01/vertebrate_pipeline/20240905/")

# Path to CSV file containing locations of all genes
csv_file <- "/Volumes/wengpj01/amphibian-TAS2R/results/coordinate_analysis/plot_generation/position_within_chromosome_all_genes_added_big_1109.csv"

## File with phylogenetics info about each assembly
genes_accessions_input <- "/Volumes/wengpj01/amphibian-TAS2R/results/coordinate_analysis/plot_generation/match_accession_order_and_class.csv"

## File defining the range surrounding each singleton or cluster
clusters <- read.csv("/Volumes/wengpj01/amphibian-TAS2R/results/coordinate_analysis/plot_generation/maximal_range_all.csv")

## No changes necessary beyond this point
#########

genes_accessions <- read_csv(genes_accessions_input)


df_location <- read_csv(csv_file)
colnames(df_location) <- c("accession", "chromosome","start", "mini", "ends")

combo_df <- merge(genes_accessions, df_location, by.x = "Accession", by.y = "accession", all.x = TRUE)
write.csv(combo_df, file = "end_summary2.csv", row.names = FALSE)

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

write.csv(filtered_data, file = "end_summary.csv", row.names = FALSE)


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


median(not_amphibia_subset$ends)

accession_counts <- subset_df %>%
  group_by(Accession) %>%
  summarise(entry_count = n())

# Calculate the average value of 'ends' per accession
average_ends <- subset_df %>%
  group_by(Accession) %>%
  summarise(average_ends = mean(ends))

# Perform a correlation test
correlation_test <- cor.test(accession_counts$entry_count, average_ends$average_ends)

p_value <- correlation_test$p.value

# Create a scatter plot
pdf("ends_vs_genes.pdf")
ggplot(data = merge(accession_counts, average_ends, by = "Accession"), aes(y = entry_count, x = average_ends)) +
  geom_point() +
  labs(x = "Average Distance from Ends", y = "Average Number Genes") +
  ggtitle("Scatter Plot of Average Ends vs. Average Number of Genes")
dev.off()


# This removes rows where all elements are NA
data <- data[!apply(data, 1, function(row) all(is.na(row))), ]


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


