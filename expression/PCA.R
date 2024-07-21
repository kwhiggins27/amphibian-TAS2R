library(readr)
library(dplyr)
library(ggfortify)
setwd("/Volumes/wengpj01/expression_pipeline/results_20231004/")

axolotl_file <- "axolotl_old.csv"
axolotl_raw <- read_csv(axolotl_file)


# Skip the first column, set the second column as row names, and select columns 3 to the end
axolotl_processed <- axolotl_raw[-1, -1, drop = FALSE]

# Set rownames to the values in the second column
axolotl_processed <- as.data.frame(axolotl_processed)
rownames(axolotl_processed) <- axolotl_processed$id
axolotl_processed <- axolotl_processed[,-1]

# Transpose the axolotl_processed data frame (or tibble)
transposed_axolotl <- t(axolotl_processed)
transposed_axolotl <- as.data.frame(transposed_axolotl)
tissue_values <- c(
  "1tongue", "2brain", "3stomach", "4intestines", "5vskin", "6dskin", "7liver",
  "1tongue", "2brain", "3stomach", "4intestines", "5vskin", "6dskin", "7liver",
  "1tongue", "2brain", "3stomach", "4intestines", "5vskin", "6dskin", "7liver"
)
bullalt <-c(
    "fresh", "fresh", "later", "fresh", "fresh", "fresh", "later",
    "fresh", "fresh", "later", "later", "later", "fresh", "later",
    "fresh", "fresh", "later", "later", "fresh", "fresh", "later"
  )

bullRIN <-c(
  "good", "good", "poor", "good", "good", "good", "good",
  "good", "good", "poor", "good", "good", "good", "good",
  "good", "good", "good", "good", "good", "good", "poor"
)

axRIN <-c(
  "good", "good", "good", "good", "good", "good", "good",
  "good", "good", "good", "good", "good", "good", "poor",
  "good", "good", "good", "good", "good", "good", "good"
)

# Add the 'tissue' column to transposed_axolotl
# transposed_axolotl['tissue'] <- tissue_values
transposed_axolotl <- cbind(tissue_values, axRIN,transposed_axolotl)

# Find columns with all zeros
zero_columns <- which(colSums(transposed_axolotl == 0) == nrow(transposed_axolotl))

# Remove the identified columns
transposed_axolotl <- transposed_axolotl[,-zero_columns] #skip for xenopus
df <- transposed_axolotl[,-1]
df <- df[,-1]

pca_res <- prcomp(df, scale. = TRUE)
pdf("Taxolotl_PCA_RIN.pdf")
autoplot(pca_res, data=transposed_axolotl, colour="tissue_values", shape="axRIN")
dev.off()

