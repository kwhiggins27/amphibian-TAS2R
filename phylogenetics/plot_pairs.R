# Load required library
library(ggplot2)
library(dplyr)

# Read the CSV file and assign column names
data <- read.csv("/Volumes/wengpj01/phylogenetics/orth_tree/summary_manypairs_forplot.csv")
colnames(data) <- c("species", "onetoone", "onetomany", "manytomany", "onetozero", "manytozero", "total")

# Function to create stacked bar plot and save as PDF
create_stacked_bar_plot_and_save <- function(species1, species2, filename) {
  # Filter data for the specified species
  species_data <- subset(data, species %in% c(species1, species2))
  
  # Subset the data for the selected species
  plot_data <- species_data[species_data$species %in% c(species1, species2), ]
  plot_data <- plot_data[, c("species", "manytozero", "onetozero", "manytomany", "onetomany", "onetoone")]
  
  # Reshape data for ggplot
  plot_data <- reshape2::melt(plot_data, id.vars = "species")
  # Define colors
  colors <- c("#97ADDA", "#0773B3", "#B6DCAF", "#699E8D", "#CF966F")
  # Create a stacked bar plot using ggplot2
  plot <- ggplot(plot_data, aes(fill = variable, y = value, x = species)) +
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_manual(values = colors) +
    labs(x = "Species", y = "Values", fill = "Category", title = "Comparison of Species Pairs") +
    theme_minimal() +
    theme(legend.position = "top") +
    theme(legend.position = "top",
          axis.line = element_line(color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  
  # Save the plot as a PDF file
  ggsave(filename, plot, width = 8, height = 6)
}

# Example: Generate stacked bar plot for species "human" and "marmoset" and save as PDF
Species_A <- "human"
Species_B <- "marmoset"
file_path <- "/Volumes/wengpj01/phylogenetics/orth_tree/pairwise/"
file_name <- paste0(Species_A, "_", Species_B, ".pdf")
output_file <- paste0(file_path, file_name)

create_stacked_bar_plot_and_save(Species_A, Species_B, output_file)

