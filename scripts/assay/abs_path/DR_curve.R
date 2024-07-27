# Load required R packages
library(tidyverse)
library(drc)
library(cowplot)
library(grid)
library(gridExtra)
# Set the working directory
setwd("/Volumes/wengpj01/assay/20231127")

# Load the raw data
raw_file <- "/Volumes/wengpj01/assay/data_DRcurve_ver231125.csv"
raw <- read.csv(raw_file, header = TRUE, row.names = NULL)

# Load the converter data
converter_file <- "/Volumes/wengpj01/assay/Converter_file.csv"
converter <- read.csv(converter_file)

# Merge the dataframes based on the matching condition
raw <- merge(raw, converter, by.x = "gene", by.y = "Aki", all.x = TRUE)

# Rename the columns
names(raw)[names(raw) == "gene"] <- "old_gene"
names(raw)[names(raw) == "Paper"] <- "gene"

raw["AUCe4"] = raw["AUC"]/10000
raw["concentration1000"]= raw["concentration"] * 10000
# names(raw)[names(raw) == "Seq"] <- NULL

# Create a new column for log2 concentration
raw$log10dose <- ifelse(raw$concentration1000 > 0, log10(raw$concentration1000), 0)

# Function to plot concentration vs. AUC
plot_concentration_vs_AUC <- function(dataframe, compound_name, gene_name, output_name = "output.pdf") {
  # Filter the dataframe to select rows that match the given criteria
  subset_match <- subset(dataframe, compound == compound_name & gene == gene_name)
  
  # Check if there are any matching rows
  if (nrow(subset_match) == 0) {
    cat(paste("No matching data found for compound '", compound_name, "' and gene '", gene_name, "'.\n", sep = ""))
  } else {
    # Clean the data by removing rows with missing or zero values in concentration or AUC
    subset_match <- subset_match %>% filter(concentration1000 > 0, !is.na(AUC))
    
    # Check if there are any valid data points after cleaning
    if (nrow(subset_match) == 0) {
      cat("No valid data points to plot after data cleaning.\n")
    } else {
      # Create the plot
      p <- ggplot(subset_match, aes(x = concentration1000, y = AUC)) +
        geom_point() +
        labs(title = paste(compound_name, " Concentration vs AUC for ", gene_name),
             x = "Concentration (log10 scale)",
             y = "AUC") +
        scale_x_continuous(trans = "log10", breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", scales::math_format(10^.x))) +
        theme_minimal()
      
      # Check for valid data points
      if (sum(!is.na(subset_match$AUCe4)) > 0) {
        # Save the plot to the specified output_name
        pdf(output_name)
        plot(p)
        dev.off()
      } else {
        cat("No valid data points for linear regression.\n")
      }
    }
  }
}


# Call the function to generate the plot
plot_concentration_vs_AUC(raw, "Marinobufagenin", "cane_56", "test_DR.pdf")

afl_cane_56 <- subset(raw, gene == "cane_56" & compound == "AlfatoxinB1")
afl_cane_56_mdl <- drm(AUC ~ concentration, data =MBG_cane_56, fct = LL.4())

pdf("afl_cane_56_mdl.pdf")
plot(afl_cane_56_mdl,broken=TRUE,type="all")
dev.off()

#######
# Create data frames for each compound
aflatoxin_data <- subset(raw, compound == "AflatoxinB1")
cinobufagin_data <- subset(raw, compound == "Cinobufagin")
marinobufagenin_data <- subset(raw, compound == "Marinobufagenin")
heliotrine_data <-subset(raw, compound == "Heliotrine")

# Create a new data frame that combines the three data frames
desired_order <- c("AflatoxinB1", "Cinobufagin", "Marinobufagenin","Heliotrine")

combined_data <- rbind(aflatoxin_data, cinobufagin_data, marinobufagenin_data)
combined_data$compound <- factor(combined_data$compound, levels = desired_order)


# Get unique combinations of gene and compound
unique_combinations <- unique(combined_data[, c("gene", "compound")])

# Set the number of columns for the grid
num_columns <- 3

# Set the desired height in inches
pdf_height <- 40  # Adjust this value as needed

# Create a function to save plots to a PDF file
save_plots_to_pdf <- function(plots, pdf_filename) {
  pdf(pdf_filename, height = pdf_height)  # Set the height here
  print(plots)
  dev.off()
}

# Initialize an empty list to store grobs
mdl_list <- list()

# Iterate through unique combinations and create and save grobs
for (i in 1:nrow(unique_combinations)) {
  gene_name <- unique_combinations[i, "gene"]
  compound_name <- unique_combinations[i, "compound"]
  
  # Create a subset for the current combination
  subset_data <- subset(raw, gene == gene_name & compound == compound_name)
  
  # Check if there are data points for the combination
  if (nrow(subset_data) > 0) {
    # Create a plot for the current combination
    mdl <- drm(AUCe4 ~ concentration1000, data = subset_data, fct = LL.4())
    
    # Append mdl to mdl_list
    mdl_list <- append(mdl_list, list(mdl))
  }
}

# Initialize an empty list to store plots
plots_list <- list()

# Iterate through unique combinations and generate plots
for (i in 1:nrow(unique_combinations)) {
  gene_name <- unique_combinations[i, "gene"]
  compound_name <- unique_combinations[i, "compound"]
  
  # Create a subset for the current combination
  subset_data <- subset(raw, gene == gene_name & compound == compound_name)
  
  # Check if there are data points for the combination
  if (nrow(subset_data) > 0) {
    # Calculate mean and standard error
    summary_data <- subset_data %>%
      group_by(concentration1000) %>%
      summarise(mean_AUC = mean(AUCe4), mean_se = sd(AUCe4) / sqrt(n()))  # Calculating standard error
    
    # Create a ggplot2 plot for the current combination with mean and standard error
    p <- ggplot(data = summary_data, aes(x = concentration1000, y = mean_AUC)) +
      geom_errorbar(aes(ymin = mean_AUC - mean_se, ymax = mean_AUC + mean_se), width = 0.1) +
      geom_point() +
      labs(title = paste(gene_name, " vs. ", compound_name),
           x = "Concentration (log10 scale, uM)",
           y = "AUC x 10^-4") +
      scale_x_continuous(trans = "log10", breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      theme_minimal()
    
    # Extract the mdl from mdl_list for this combination
    mdl <- mdl_list[[i]]
    
    # Extract the coefficients b, c, d, and e from the model
    b_coeff <- coef(mdl)["b:(Intercept)"]
    c_coeff <- coef(mdl)["c:(Intercept)"]
    d_coeff <- coef(mdl)["d:(Intercept)"]
    e_coeff <- coef(mdl)["e:(Intercept)"]
    
    # Define the curve using the coefficients
    curve_data <- data.frame(concentration1000 = summary_data$concentration1000)
    curve_data$AUCe4 <- predict(mdl, newdata = curve_data)
    
    # Calculate the range of concentration values
    min_concentration <- min(curve_data$concentration1000)
    max_concentration <- max(curve_data$concentration1000)
    
    # Create a new data frame with more points within the same range
    new_curve_data <- data.frame(concentration1000 = seq(min_concentration, max_concentration, length.out = 10 * nrow(curve_data)))
    
    # Calculate predicted AUC values for the new data
    new_curve_data$AUCe4 <- predict(mdl, newdata = new_curve_data)
    
    # Add the curve to the plot
    p <- p + geom_line(data = new_curve_data, aes(x = concentration1000, y = AUCe4), color = "red")
    
    # Add the ggplot2 plot to the list
    plots_list[[length(plots_list) + 1]] <- p
  }
}

# Arrange the ggplot2 plots into a grid
grid_plots <- cowplot::plot_grid(plotlist = plots_list, ncol = num_columns)



# Save the grid of plots to a PDF file
pdf("plots_grid9.pdf")
print(grid_plots)
dev.off()

#####
library(dplyr)
# Subset data for each compound
aflatoxin_data <- subset(raw, compound == "AflatoxinB1")
cinobufagin_data <- subset(raw, compound == "Cinobufagin")
marinobufagenin_data <- subset(raw, compound == "Marinobufagenin")
heliotrine_data <-subset(raw, compound == "Heliotrine")

# Create a new data frame that combines the three data frames
desired_order <- c("AflatoxinB1", "Cinobufagin","Heliotrine", "Marinobufagenin")

combined_data <- bind_rows(aflatoxin_data, cinobufagin_data, heliotrine_data, marinobufagenin_data,)
combined_data$compound <- factor(combined_data$compound, levels = desired_order)

# Get unique combinations of gene and compound
unique_combinations <- unique(combined_data[, c("gene", "compound")])

# Initialize an empty list to store models
mdl_list <- list()

# Initialize an empty list to store final plots
final_plots <- list()
unique_genes <- unique(combined_data$gene)

# Iterate through unique genes to generate combined plots for each gene
for (gene_name in unique_genes) {
  # Create a subset for the current gene
  gene_data <- subset(combined_data, gene == gene_name)
  
  # Check if there are data points for the gene
  if (nrow(gene_data) > 0) {
    # Create an empty list to store compound-specific plots for the gene
    gene_plots <- list()
    
    # Iterate through unique compounds within the gene
    unique_compounds <- unique(gene_data$compound)
    
    for (compound_name in unique_compounds) {
      # Create a subset for the current compound within the gene
      compound_data <- subset(gene_data, compound == compound_name)
      
      # Check if there are data points for the compound
      if (nrow(compound_data) > 0) {
        # Fit model for the compound
        mdl <- drm(AUCe4 ~ concentration1000, data = compound_data, fct = LL.4())
        
        # Predicted values using the model
        curve_data <- data.frame(concentration1000 = compound_data$concentration1000)
        curve_data$AUCe4 <- predict(mdl, newdata = curve_data)
        
        # Store the compound-specific plot
        gene_plots[[compound_name]] <- list(data = compound_data, curve = curve_data)
      }
    }
    
    # Store the gene-specific plot containing data and curves for all compounds
    final_plots[[gene_name]] <- gene_plots
  }
}
# Define a custom color palette for compounds
compound_palette <- c("AflatoxinB1" = "#D55E00", 
                      "Cinobufagin" = "#166A53", 
                      "Heliotrine" = "#063781",
                      "Marinobufagenin" = "#AC71F3") 

shape_palette <- c("AflatoxinB1" ="circle",#"19",  # Filled circle
                     "Cinobufagin" = "square",#5,  # Filled diamond
                      "Marinobufagenin" ="triangle", #2,  # Filled triangle point-up
                      "Heliotrine" ="diamond")


# Iterate through final_plots and generate plots for each gene
plot_list <- list()

for (gene_name in unique_genes) {
  gene_plots <- final_plots[[gene_name]]
  
  # Create a new ggplot for the gene
  p <- ggplot() +
    labs(title = paste("Combined Plot for", gene_name),
         x = "Concentration (log10 scale, uM)",
         y = "AUC x 10^-4") +
    scale_x_continuous(trans = "log10", breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    theme_minimal()
  
  # Add mean, standard deviation error bars, and model curves for each gene/compound pair to the plot
  for (compound_name in names(gene_plots)) {
    compound_data <- gene_plots[[compound_name]]$data
    compound_curve <- gene_plots[[compound_name]]$curve
    
    # Calculate mean and standard deviation for the compound data
    summary_data <- compound_data %>%
      group_by(concentration1000) %>%
      summarise(mean_AUC = mean(AUCe4), sd_AUC = sd(AUCe4), mean_se = sd(AUCe4) / sqrt(n()))  # Calculating standard deviation and standard error
    
    # Assigning shapes to different compounds
    compound_shape <- ifelse(compound_name == "AflatoxinB1", "circle",#"19",  # Filled circle
                             ifelse(compound_name == "Cinobufagin", "square",#5,  # Filled diamond
                                    ifelse(compound_name == "Marinobufenagin", "triangle", #2,  # Filled triangle point-up
                                           ifelse(compound_name == "Heliotrine", "diamond", "square"))))#"3, 16"))))  # Filled triangle point-down, a different shape as default
    # Add mean, standard deviation error bars, and model curves to the plot, mapping shape to the compound
    # p <- p + geom_errorbar(data = summary_data, aes(x = concentration1000, ymin = mean_AUC - mean_se, ymax = mean_AUC + mean_se), width = 0.1, color = compound_palette[compound_name]) +
    #   geom_point(data = summary_data, aes(x = concentration1000, y = mean_AUC, shape = compound_shape), color = unname(compound_palette[compound_name]), show.legend = TRUE) +
    #   geom_line(data = compound_curve, aes(x = concentration1000, y = AUCe4), color = unname(compound_palette[compound_name]), show.legend = FALSE) +
    #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    # Add mean, standard deviation error bars to the plot, mapping shape to the compound
    p <- p + geom_errorbar(data = summary_data, aes(x = concentration1000, ymin = mean_AUC - mean_se, ymax = mean_AUC + mean_se), width = 0.1, color = compound_palette[compound_name])
    
    # Add points with shapes for each compound, omitting shapes for the lines
    p <- p + geom_point(data = summary_data, aes(x = concentration1000, y = mean_AUC), shape = unname(shape_palette[compound_name]), color = unname(compound_palette[compound_name]), size = 3, show.legend =FALSE)
    
    # Add lines for model curves
    p <- p + geom_line(data = compound_curve, aes(x = concentration1000, y = AUCe4), color = unname(compound_palette[compound_name]), show.legend = FALSE)

  # Remove legend for the lines (optional)
    p <- p + guides(color = FALSE)
    
    p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  }
  # Add the plot to the plot_list
  plot_list[[gene_name]] <- p
}

# Arrange plots in a grid layout
pdf("combined_gene_plots_grid4.pdf", width = 10, height = 8)  # Adjust width and height as needed
multiplot(plotlist = plot_list, cols = 3)
dev.off()


