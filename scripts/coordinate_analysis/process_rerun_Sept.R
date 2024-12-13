library(dplyr)
library(readr)
library(taxize)
library(phytools)
library(ape)
library(ggplot2)

## Set your working directory.  This is where plots will be saved.
setwd("/Volumes/wengpj01/vertebrate_pipeline/20240905")

##Key input files -- update with absolute path of each.  Examples included.
genes_accessions_input <- "/Volumes/wengpj01/amphibian-TAS2R/results/coordinate_analysis/plot_generation/20231006_basics_addedplottingclade.csv"
clusters_accessions_input <- "/Volumes/wengpj01/amphibian-TAS2R/results/coordinate_analysis/plot_generation/summary_clusters_withnearest_1000000.csv"
genome_info_input <- "/Volumes/wengpj01/amphibian-TAS2R/results/coordinate_analysis/plot_generation/accessions_pt1_assembly_combo.csv"
position_coordinates_input <- "/Volumes/wengpj01/amphibian-TAS2R/results/coordinate_analysis/plot_generation/position_within_chromosome_1005.csv"
genome_size_input <- "/Volumes/wengpj01/amphibian-TAS2R/results/coordinate_analysis/plot_generation/genome_size.csv"

## Set your API key.  This is only necessary for the slow step in 72-75.  I recommend running this once and then commenting out
api_key <- Sys.getenv("CHANGE THIS")

## Intermediate file.  Can save and reload this to restart after the slow step in 72-75
int_file <- "/Volumes/wengpj01/amphibian-TAS2R/results/coordinate_analysis/plot_generation/tax_info_RM_batch.csv"

##No changes necessary after this point
###########

#Read in all dataframes
genes_accessions <- read_csv(genes_accessions_input)
clusters_accessions <- read_csv(clusters_accessions_input)
genome_info <- read_csv(genome_info_input)
genome_info <- genome_info %>%
  select(1:6)
position_coordinates <- read_csv(position_coordinates_input)
genome_size <- read_csv(genome_size_input, col_names = c("common", "C_value", "method", "accession", "TAS2R_count", "order"))

# Merge genes_accessions with clusters_accessions by "Accessions" using a left join
merged_temp <- merge(genes_accessions, clusters_accessions, by = "Accession", all.x = TRUE)

# Merge the result with genome_info by "Accessions" using a left join
merged_main <- merge(merged_temp, genome_info, by.x = "Accession", by.y = "accession", all.x = TRUE)
merged_main2 <- merge(merged_main, position_coordinates, by.x = "Accession", by.y = "accession", all.x = TRUE)

# Remove the "Number of Genes.y" column
merged_main <- merged_main %>%
  subset(select = -`Number of Genes.y`)

# Rename "Number of Genes.x" to "genes"
merged_main <- merged_main %>%
  rename(genes = `Number of Genes.x`)

# Rename "Number of Clusters" to "clusters"
merged_main <- merged_main %>%
  rename(clusters = `Number of Clusters`)

merged_main_no_duplicates <- merged_main[!duplicated(merged_main), ]



##Position on chromosome code (from is_gene_near_end.R)

get_taxonomic_info <- function(latin_name) {
  # Get the taxonomic information for the Latin name from NCBI
  tax_info <- tax_name(latin_name, get = c("class", "order", "family"), db = "ncbi")

  # Extract the order and family names (assuming they are single names)
  class_name <- ifelse(length(tax_info$class) > 0, tax_info$class[[1]], NA)
  order_name <- ifelse(length(tax_info$order) > 0, tax_info$order[[1]], NA)
  family_name <- ifelse(length(tax_info$family) > 0, tax_info$family[[1]], NA)

  # Pause for 0.25 seconds
  Sys.sleep(0.25)

  return(c(class = class_name, order = order_name, family = family_name))
}

## Apply the function to the dataframe
# taxonomic_info <- sapply(unique_latin, get_taxonomic_info)
# taxonomic_info <- t(taxonomic_info)
# int_file <- sapply(merged_main_no_duplicates$latin, get_taxonomic_info)
taxonomic_info_2 <- read_csv(int_file)

# write.csv(taxonomic_info_2, file = "tax_info_RM_batch.csv", row.names = FALSE)
#Transpose taxonomic_info_2
transposed_tax_info <- t(taxonomic_info_2)
colnames(transposed_tax_info) <- c("class", "order", "family")

# Check if they have the same number of rows
if (nrow(merged_main_no_duplicates) == nrow(transposed_tax_info)) {
  # Bind the two dataframes together horizontally
  merged_with_taxa <- cbind(merged_main_no_duplicates, transposed_tax_info)
} else {
  print("The number of rows is not the same after transposition.")
}

merged_with_taxa$genes <- as.numeric(as.character(merged_with_taxa$genes))
merged_with_taxa$clusters <- as.numeric(as.character(merged_with_taxa$clusters))
merged_with_taxa$FractionClustered <- as.numeric(as.character(merged_with_taxa$FractionClustered))


merged_with_taxa$genes_per_cluster <- (merged_with_taxa$genes / merged_with_taxa$clusters) * merged_with_taxa$FractionClustered

merged_with_taxa$class <- ifelse(
  is.na(merged_with_taxa$class) & (merged_with_taxa$order == "Testudines" | merged_with_taxa$order == "Crocodylia"),
  "Lepidosauria",
  ifelse(
    is.na(merged_with_taxa$class) & merged_with_taxa$order == "Ceratodontiformes",
    "Ceratodontoidei",
    merged_with_taxa$class
  )
)

# Subset the DataFrame to include only "Accession" and "Order" columns
subset_df1 <- merged_with_taxa[, c("Accession", "class")]

# Specify the file path where you want to save the subset DataFrame
output_file <- "match_accession_class2.csv"

# Export the subset DataFrame to a CSV file
write.csv(subset_df1, file = output_file, row.names = FALSE)

clade_colors_df <- data.frame(
  PlottingClade = c(
    "Mammalia", "Aves", "Squamata", "Testudines", "Crocodylia", "Anura","Caudata", 
    "Gymnophiona", "Actinopteri", "Chondrichthyes",
    "Hyperoartia", "Ceratodontoidei", "Cladistia"
  ),
  HexCode = c(
    "#D55E00", "#5B066B", "#AC71F3",
    "#9245E3", "#EC7EEF",  "#166A53", "#AEF2B5", "#36AB59","#063781", "#56B4E9", "#0072B2", "#ACF1FF",
    "#9EBADC"
  )
)

# Merge the color mapping dataframe with your data
merged_with_taxa <- merge(merged_with_taxa, clade_colors_df, by = "PlottingClade", all.x = TRUE)


amphibia_subset <- merged_with_taxa[merged_with_taxa$class == "Amphibia", ]
amphibia_subset <- amphibia_subset[complete.cases(amphibia_subset), ]
reptile_high_family <- c("Cordylidae", "Dactyloidae", "Eublepharidae", "Lacertidae", "Phrynosomatidae", "Rhineuridae", "Sphaerodactylidae")
not_amphibia_subset <- subset(merged_with_taxa, class != "Amphibia")
high_reptile <- merged_with_taxa[merged_with_taxa$family %in% reptile_high_family, ]
fr_sal <- subset(merged_with_taxa, order %in% c("Anura", "Caudata"))
has_clusters <- merged_with_taxa[!(is.na(merged_with_taxa$clusters) | merged_with_taxa$clusters == ""), ]


#Fraction in clusters
matrix_data <- matrix(merged_with_taxa$'genes', nrow = nrow(merged_with_taxa), ncol = 1, dimnames = list(merged_with_taxa$latin, "genes"))
graphics.off()



# Create a box and whisker plot for genes_by_class
pdf("genes_by_class_boxplot6.pdf")
ggplot(merged_with_taxa, aes(
  x = factor(PlottingClade, levels = clade_colors_df$PlottingClade), 
  y = `Number of Genes.x`, 
  fill = factor(PlottingClade, levels = clade_colors_df$PlottingClade)
)) +
  geom_boxplot() +
  scale_fill_manual(values = clade_colors_df$HexCode) +
  labs(x = "PlottingClade", y = "Number of Genes", title = "Boxplots by PlottingClade") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    axis.line = element_line(),  # Show x and y axes lines
    axis.ticks = element_line()  # Show ticks on the axes
  ) +
  theme(legend.position = "none") # Remove the legend
dev.off()


pdf("clusters_by_class_boxplot6.pdf")
ggplot(has_clusters, aes(x = factor(PlottingClade, levels = clade_colors_df$PlottingClade), y = `clusters`, fill = factor(PlottingClade, levels = clade_colors_df$PlottingClade))) +
  geom_boxplot() +
  scale_fill_manual(values = clade_colors_df$HexCode) +
  labs(x = "PlottingClade", y = "Number of Clusters", title = "Boxplots by PlottingClade") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) +  # Remove the grid 
  theme(legend.position = "none") # Print the plot to the PDF
dev.off()

pdf("genes_per_cluster6.pdf")
ggplot(has_clusters, aes(x = factor(PlottingClade, levels = clade_colors_df$PlottingClade), y = `genes_per_cluster`, fill = factor(PlottingClade, levels = clade_colors_df$PlottingClade))) +
  geom_boxplot() +
  scale_fill_manual(values = clade_colors_df$HexCode) +
  labs(x = "PlottingClade", y = "Genes per cluster", title = "Boxplots by PlottingClade") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) +  # Remove the grid 
  theme(legend.position = "none") # Print the plot to the PDF
dev.off()

pdf("nearest6.pdf")
ggplot(has_clusters, aes(x = factor(PlottingClade, levels = clade_colors_df$PlottingClade), y = `Nearest_in_cluster`, fill = factor(PlottingClade, levels = clade_colors_df$PlottingClade))) +
  geom_boxplot() +
  scale_fill_manual(values = clade_colors_df$HexCode) +
  labs(x = "PlottingClade", y = "Number of Genes", title = "Boxplots by PlottingClade") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) +  # Remove the grid 
  theme(legend.position = "none") # Print the plot to the PDF
dev.off()

merged_with_taxa$FractionClustered <- ifelse(merged_with_taxa$`Number of Genes.x` > 0 & is.na(merged_with_taxa$`Fraction Clustered`), 0, merged_with_taxa$`Fraction Clustered`)
subset_merged <- merged_with_taxa[merged_with_taxa$`Number of Genes.x` > 0, ]

pdf("fraction_clustered.pdf")
ggplot(subset_merged, aes(x = factor(PlottingClade, levels = clade_colors_df$PlottingClade), y = `Fraction Clustered`, fill = factor(PlottingClade, levels = clade_colors_df$PlottingClade))) +
  geom_boxplot() +
  scale_fill_manual(values = clade_colors_df$HexCode) +
  labs(x = "PlottingClade", y = "Fraction Clustered", title = "Boxplots by PlottingClade") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) +  # Remove the grid 
  theme(legend.position = "none") # Print the plot to the PDF
dev.off()

pdf("cluster_spacing.pdf")
ggplot(has_clusters, aes(x = factor(PlottingClade, levels = clade_colors_df$PlottingClade), y = `Total Space (kb/receptor)`, fill = factor(PlottingClade, levels = clade_colors_df$PlottingClade))) +
  geom_boxplot() +
  scale_fill_manual(values = clade_colors_df$HexCode) +
  labs(x = "PlottingClade", y = "Total Space (kb/receptor)", title = "Boxplots by PlottingClade") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) +  # Remove the grid 
  theme(legend.position = "none") # Print the plot to the PDF
dev.off()



cor.test(merged_with_taxa$genes, merged_with_taxa$contig_n50)

# Remove everything after the space and convert to float
merged_with_taxa$coverage <- as.integer(sub(" .*", "", merged_with_taxa$coverage))


# Filter rows where both columns are numeric
merged_with_taxa2 <- merged_with_taxa[!is.na(merged_with_taxa$total_length) & !is.na(merged_with_taxa$genes), ]
merged_with_taxa2_coverage <- merged_with_taxa2[!is.na(merged_with_taxa2$coverage),]

cor.test(merged_with_taxa2_coverage$genes, merged_with_taxa2_coverage$coverage)

####Genome size
subset_amphibian <- genome_size[genome_size$order == "amphibian", ]
subset_not_amphibian <- genome_size[genome_size$order != "amphibian", ]

cor_test_not_amphibian <- cor.test(subset_not_amphibian$C_value, subset_not_amphibian$TAS2R_count)
cor_test_amphibian <- cor.test(subset_amphibian$C_value, subset_amphibian$TAS2R_count)
cor_test_all <- cor.test(genome_size$C_value, genome_size$TAS2R_count)

