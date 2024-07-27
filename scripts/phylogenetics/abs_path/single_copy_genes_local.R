singlibrary(paleotree)
library(dplyr)
library(ape)
library(phytools)
library(tidyr)
library(ggplot2)
library(readr)

setwd("/Users/kwhiggins27/Desktop/OneDrive - Harvard University/Weng Lab/Amphibian_paper/NZ_orth")

# Create an empty dataframe to store overall results
overall_results <- data.frame(
  gene_A = character(),
  gene_B = character(),
  taxon = character(),
  result_message = character(),
  stringsAsFactors = FALSE
)
# Read in the different trees for each taxon
amphibian_tree <- read.tree("amphibians2.newick")
bird_tree <- read.tree("birds2.newick")
mammal_tree <- read.tree("mammals2.newick")
reptile_tree <- read.tree("reptile.newick")
fish_tree <- read.tree("rayfinned.newick")

# Define ordered_tips for each taxon
ordered_tips_amphibian <- rev(amphibian_tree$tip.label)
ordered_tips_bird <- rev(bird_tree$tip.label)
ordered_tips_mammal <- rev(mammal_tree$tip.label)
ordered_tips_reptile <- rev(reptile_tree$tip.label)
ordered_tips_fish <- rev(fish_tree$tip.label)

# Create empty dataframes for each taxon
df_ordered_tips4_amphibian <- data.frame(gene = character(), scg = integer())
df_ordered_tips4_bird <- data.frame(gene = character(), scg = integer())
df_ordered_tips4_mammal <- data.frame(gene = character(), scg = integer())
df_ordered_tips4_reptile <- data.frame(gene = character(), scg = integer())
df_ordered_tips4_fish <- data.frame(gene = character(), scg = integer())

species_tree <- read.tree("Timetree_ultrametric.tre")

genes_accessions_input <- "match_accession_order_and_class.csv"
genes_accessions <- read_csv(genes_accessions_input)


# Create counters to track gene family names for each taxon
amphibian_counter <- 1
bird_counter <- 1
mammal_counter <- 1
reptile_counter <- 1
fish_counter <- 1

# Function to assign gene family names based on taxon
assign_gene_family_name <- function(taxon) {
  switch(
    taxon,
    amphibian = {
      name <- paste("amphibian", amphibian_counter, sep = "_")
      amphibian_counter <<- amphibian_counter + 1
      return(name)
    },
    bird = {
      name <- paste("bird", bird_counter, sep = "_")
      bird_counter <<- bird_counter + 1
      return(name)
    },
    mammal = {
      name <- paste("mammal", mammal_counter, sep = "_")
      mammal_counter <<- mammal_counter + 1
      return(name)
    },
    reptile = {
      name <- paste("reptile", reptile_counter, sep = "_")
      reptile_counter <<- reptile_counter + 1
      return(name)
    },
    fish = {
      name <- paste("fish", fish_counter, sep = "_")
      fish_counter <<- fish_counter + 1
      return(name)
    },
    "unknown"
  )
}


analyze_genes <- function(gene_A, gene_B, taxon) {
  process_taxon <- function(taxon) {
    if (taxon == "amphibian") {
      tr <- read.tree("amphibians2.newick")
    } else if (taxon == "bird") {
      tr <- read.tree("birds2.newick")
    } else if (taxon == "mammal") {
      tr <- read.tree("mammals2.newick")
    } else if (taxon == "reptile") {
      tr <- read.tree("reptile.newick")
    } else if (taxon == "fish") {
      tr <- read.tree("rayfinned.newick")
    } else {
      print("Error: wrong taxa")
      return(NULL)
    }
    return(tr)
  }
  print(taxon)
  tr <- process_taxon(taxon)
  if (!is.null(tr)) {
    print("Tree successfully loaded.")
  }
  
  
  # find_latin_value <- function(accession) {
  #   file_path <- paste("/Volumes/wengpj01/vertebrate_pipeline/subdirs/", accession, "/for_py2.py", sep = "")
  #   found_line <- NULL
  #   con <- file(file_path, "r")
  # 
  #   while (TRUE) {
  #     line <- readLines(con, n = 1)
  #     if (length(line) == 0) {
  #       break
  #     }
  #     match <- regmatches(line, regexec('latin="([^"]+)"', line))
  #     if (length(match[[1]]) > 1) {
  #       found_line <- match[[1]][2]
  #       break
  #     }
  #   }
  #   close(con)
# 
#     if (!is.null(found_line)) {
#       return(found_line)
#     } else {
#       return("Value for 'latin' not found or does not match the expected pattern in the file.")
#     }
#   }
  
  find_latin_value2 <- function(ordered_tips4, genes_accessions) {
    latin_values <- character(length(ordered_tips4))  # Initialize an empty vector to store latin values
    
    for (i in seq_along(ordered_tips4)) {
      var <- ordered_tips4[i]
      # Find row indices where Accession matches var
      matching_rows <- which(genes_accessions$Accession == var)
      
      if (length(matching_rows) > 0) {
        # Extract corresponding latin value from matching row in genes_accessions dataframe
        latin_values[i] <- genes_accessions$latin[matching_rows[1]]  # Consider the first match if multiple matches exist
      } else {
        latin_values[i] <- NA  # If no match is found, assign NA
      }
    }
    
    return(latin_values)
  }
  
  get_latin_from_gene <- function(gene_name) {
    gene_split <- unlist(strsplit(gene_name, "_"))
    accession <- paste(gene_split[1:3], collapse = "_")
    accession <- gsub("^(.+?_.+?)_", "\\1.", accession)
    latin_value <- find_latin_value2(accession, genes_accessions)
    
    return(latin_value)
  }
  
  sp_a <- which(species_tree$tip.label == get_latin_from_gene(gene_A))
  sp_b <- which(species_tree$tip.label == get_latin_from_gene(gene_B))
  
  
  if (get_latin_from_gene(gene_B) == "Heteronetta_atricapilla") {
    sp_b <- which(species_tree$tip.label == "Sylvia_atricapilla")
  } else {
    print("sp_b is NULL or does not match the condition.")
  }
  mrca_sp <- getMRCA(species_tree, c(sp_a, sp_b))
  tr_sp <- extract.clade(species_tree, mrca_sp)
  tmrca <- dateNodes(species_tree)[mrca_sp]
  node_age <- unname(tmrca)
  
  ordered_tips_sp <- rev(tr_sp$tip.label)
  
  
  check_all_in <- function(values, reference) {
    all_present <- all(values %in% reference)
    return(all_present)
  }
  
  ordered_tips <- rev(tr$tip.label)
  a <- which(tr$tip.label == gene_A)
  b <- which(tr$tip.label == gene_B)
  
  mrca <- getMRCA(tr, c(a, b))
  tr4 <- extract.clade(tr, mrca)
  
  print("Done with mrca")
  
  ordered_tips4 <- rev(tr4$tip.label)
  
  print("about to start latin converter")
  latin4 <- lapply(ordered_tips4, get_latin_from_gene)
  print("finished latin converter")
  all_in_ordered_tips_sp <- check_all_in(unlist(latin4), ordered_tips_sp)
  
  if (!all_in_ordered_tips_sp) {
    print("Not the most distant genes!")
    
    # Combine all elements of latin4 and ordered_tips_sp into vectors
    all_latin_values <- unique(unlist(latin4))
    all_ordered_tips_sp_values <- unique(unlist(ordered_tips_sp))
    
    extract_info <- function(value) {
      split_value <- unlist(strsplit(value, "_"))
      
      # Check if there are more than one underscore
      if (length(split_value) > 2) {
        modified_value <- paste(split_value[1:2], collapse = "_")  # Keep info before the second underscore
      } else {
        modified_value <- value  # Keep the value as is
      }
      
      return(modified_value)
    }
    
    # Apply the function to all values in all_latin_values and all_ordered_tips_sp_values
    all_latin_values <- lapply(all_latin_values, extract_info)
    all_ordered_tips_sp_values <- lapply(all_ordered_tips_sp_values, extract_info)
    
    
    # Find values in latin4 not present in ordered_tips_sp
    values_not_in_ordered_tips_sp <- setdiff(all_latin_values, all_ordered_tips_sp_values)
    
    print(values_not_in_ordered_tips_sp)

    
    # Function to check if a match is found
    check_match <- function(placeholder) {
      which(species_tree$tip.label == placeholder)
    }
    
    # Iterate through values_not_in_ordered_tips_sp to check for matches
    matching_indices <- lapply(values_not_in_ordered_tips_sp, check_match)
    
    # Identify indices that do not have a match
    missing_indices <- which(sapply(matching_indices, length) == 0)
    missing_values <- unlist(values_not_in_ordered_tips_sp[missing_indices])
    print(missing_values)
    # Drop values without a match from all_ordered_tips_sp_values and latin4
    if (length(missing_values) > 0) {
      all_ordered_tips_sp_values5 <- all_ordered_tips_sp_values[!all_ordered_tips_sp_values %in% missing_values]
      all_latin_values5 <- all_latin_values[!all_latin_values %in% missing_values]
      latin5 <- latin4[!latin4 %in% missing_values]
    } else {
      all_latin_values5 <-all_latin_values
      all_ordered_tips_sp_values5 <- all_ordered_tips_sp_values
      latin5 <- latin4
      
    }
    
    # Initialize an empty list to store the results
    indices_list2 <- list()
    
    # Iterate through each value in values_not_in_ordered_tips_sp
    indices_list2 <- lapply(values_not_in_ordered_tips_sp, function(placeholder) {
      which(species_tree$tip.label == placeholder)
    })
    print(indices_list2)
    
    
    
    combined_values <- unique(c(sp_a, sp_b, unlist(indices_list2)))
    
    
    print(combined_values)
    mrca_sp <- getMRCA(species_tree, combined_values)
    tr_sp <- extract.clade(species_tree, mrca_sp)
    tmrca <- dateNodes(species_tree)[mrca_sp]
    print(node_age)
    node_age <- unname(tmrca)
    print(node_age)
    
    ordered_tips_sp <- rev(tr_sp$tip.label)
    all_ordered_tips_sp_values <- unique(unlist(ordered_tips_sp))  #this is new
    all_ordered_tips_sp_values <- lapply(all_ordered_tips_sp_values, extract_info)
    print("here")
    
    all_in_ordered_tips_sp <- check_all_in(unlist(all_latin_values5), all_ordered_tips_sp_values)
    
    values_not_in_ordered_tips_sp2 <- setdiff(unlist(all_latin_values5), all_ordered_tips_sp_values)
    print("still missing")
    print(values_not_in_ordered_tips_sp2)
    
    
    
    if (all_in_ordered_tips_sp) {
      print("These weren't the most distant species, but the problem was fixed")
    } else {
      print("These weren't the most distant species, and still aren't")
    }
  } else {
    all_latin_values <- unique(unlist(latin4))
    all_ordered_tips_sp_values <- unique(unlist(ordered_tips_sp))
    latin5 <- latin4
    all_ordered_tips_sp_values5 <- all_ordered_tips_sp_values
    all_latin_values5 <- all_latin_values
    # Proceed with no comment
    # Your further processing or code goes here
  }
  
  if (exists("indices_list2")) {
    all_equal <- identical(sp_a, sp_b)
    if (all_equal) {
      for (indices in indices_list2) {
        all_equal <- all_equal && identical(sp_a, indices)
        if (!all_equal) break
      }
    }
  } else {
    all_equal <- identical(sp_a, sp_b)
  }
  print("All equal")
  print(all_equal)
  
  if (all_in_ordered_tips_sp) {
    count_occurrences <- sapply(all_ordered_tips_sp_values5, function(x) sum(latin5 == x))  #changed from all_latin_values
    result_df <- data.frame(accession = all_ordered_tips_sp_values5, count = count_occurrences)
    count_exact_1 <- sum(result_df$count == 1)
    count_more_2 <- sum(result_df$count >= 2)
    prop_exact_1 <- count_exact_1 / nrow(result_df)
    prop_more_than_1 <- count_more_2 / nrow(result_df)
    
    print(count_exact_1)
    print(nrow(result_df))
    print(prop_exact_1)
    print(count_more_2)
    print(prop_more_than_1)
    print(node_age)
    print(count_exact_1)
    print("all equal")
    print(all_equal)
    
    if (all_equal) {
      result_message <- "This is just one gene."
      scg_value <- "nope"  # Set scg to 1 for single copy gene
    } else if (prop_exact_1 >= 0.9 && prop_more_than_1 <= 0.05 && tmrca >= 50 && count_exact_1 >= 2) {
      result_message <- "This is an obligate single copy gene, scg 1."
      scg_value <- 1  # Set scg to 1 for single copy gene
    } else if (prop_exact_1 >= 0.75 && prop_exact_1 < 0.9 && prop_more_than_1 <= 0.05 && tmrca >= 50 && count_exact_1 >= 2) {
      result_message <- "This is an obligate single copy gene by a relaxed definition scg 2."
      scg_value <- 2  # Set scg to 2 for single copy gene
    } else if (prop_exact_1 >= 0.5 && prop_exact_1 < 0.75 && prop_more_than_1 <= 0.05 && tmrca >= 50 && count_exact_1 >= 2) {
      result_message <- "This is an obligate single copy gene by a relaxed definition scg 3."
      scg_value <- 3  # Set scg to 2 for single copy gene
    } else if (prop_exact_1 < 0.4 && prop_more_than_1 >= 0.5 && tmrca >= 50 && count_exact_1 >= 2) {
      result_message <- "This is a multi-copy gene in many lineages."
      scg_value <- 10  # Set scg to 2 for single copy gene
    } else if (prop_exact_1 < 0.4 && prop_more_than_1 >= 0.3 && tmrca >= 50 && count_exact_1 >= 2) {
      result_message <- "This is a multi-copy gene in many lineages."
      scg_value <- 11  # Set scg to 2 for single copy gene
    } else if (tmrca >= 50 && nrow(result_df) >= 2) {
      result_message <- "This is NOT an obligate single copy gene."
      scg_value <- 0  # Set scg to 2 for single copy gene
    } else {
      result_message <- "Not the most distant genes!"
      scg_value <- "nope"  # Set scg to 0 for not a single copy gene
    }
  }

    
  # Determine the gene family name based on taxon
  print("about to assign gene family name")
  gene_family_name <- assign_gene_family_name(taxon)
  
  # Create a dataframe for ordered_tips4 with scg column
  df_ordered_tips4 <- data.frame(gene = ordered_tips4, scg = scg_value, gene_family = gene_family_name, tmrca = rep(tmrca, length(ordered_tips4)))
  
  
  # Create a dataframe for ordered_tips4 with scg and gene family name columns
  print("about to append to df")
  # Update df_ordered_tips4 based on taxon
  if (taxon == "amphibian") {
    df_ordered_tips4_amphibian <<- rbind(df_ordered_tips4_amphibian, df_ordered_tips4)
  } else if (taxon == "bird") {
    df_ordered_tips4_bird <<- rbind(df_ordered_tips4_bird, df_ordered_tips4)
  } else if (taxon == "mammal") {
    df_ordered_tips4_mammal <<- rbind(df_ordered_tips4_mammal, df_ordered_tips4)
  } else if (taxon == "reptile") {
    df_ordered_tips4_reptile <<- rbind(df_ordered_tips4_reptile, df_ordered_tips4)
  } else if (taxon == "fish") {
    df_ordered_tips4_fish <<- rbind(df_ordered_tips4_fish, df_ordered_tips4)
  }
  
  
  
  return(result_message)
}


# # Example usage:
# gene_A <- "GCA_014858855_1_Asiatic_toad_toad_CM026465_1_501747465_501748433"
# gene_B <- "GCA_004786255_1_frogs_and_toads_and_toads_CM016421_1_37234200_37235180"
# taxon <- "amphibian"
# 
# result_message <- analyze_genes(gene_A, gene_B, taxon)
# print(result_message)
# print(df_ordered_tips4_amphibian)
# ###
# 
# gene_A <- "GCA_009769605_1_Abyssinian_ground_hornbill_ground_hornbill_CM020024_1_98380618_98381565"
# gene_B <- "GCA_000146605_4_turkey_CM000962_2_185213563_185214498"
# taxon <- "bird"
# 
# result_message <- analyze_genes(gene_A, gene_B, taxon)
# print(result_message)
# 
# ###
# 
# gene_A <- "GCA_014905855_1_human_AP023465_1_10241571_10242470"
# gene_B <- "GCA_016433145_1_Agile_Gracile_Mouse_Opossum_Gracile_Mouse_Opossum_CM028235_1_439259536_439260480"
# taxon <- "mammal"
# 
# result_message <- analyze_genes(gene_A, gene_B, taxon)
# print(result_message)
# 
# ###
# gene_A <- "GCA_014905855_1_human_AP023472_1_11032450_11033406"
# gene_B <- "GCA_009760805_1_Mountain_hare_hare_VUAX01000167_1_759270_760196"
# taxon <- "mammal"
# result_message <- analyze_genes(gene_A, gene_B, taxon)
# print(result_message)


###AMPHIBIAN SECTION
taxon <- "amphibian"
###
gene_A <- "GCA_002915635_3_axolotl_axolotl_5_aa_148491355_148492293"
gene_B <- "GCA_901765095_2_caecilians_LR594645_1_34359436_34360374"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_014858855_1_Asiatic_toad_toad_JABFFT010000389_1_65011_65925"
gene_B <- "GCA_027579735_1_fire_bellied_toad_toad_fire_bellied_toad_11_aa_11707780_11708682"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#problem out of place GCA_027579735_1_fire_bellied_toad_toad_fire_bellied_toad_11_aa_11758776_11759729

gene_A <- "GCA_019512145_1_Tungara_frog_frog_WNYA01013758_1_4155_5105"
gene_B <- "GCA_004786255_1_frogs_and_toads_and_toads_CM016422_1_1136010_1136927"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_029499605_1_Sardinian_treefrog_treefrog_CM056042_1_230491015_230491950"
gene_B <- "GCA_027358695_2_plains_spadefoot_toad_spadefoot_toad_CM049776_1_46916759_46917670"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_905171765_1_common_toad_toad_LR991673_1_5924757_5925635"
gene_B <- "GCA_905171775_1_common_frog_frog_LR991685_1_2720396_2721289"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_019512145_1_Tungara_frog_frog_WNYA01000075_1_181447_182346"
gene_B <- "GCA_028390025_1_corroboree_frog_frog_corroboree_frog_7_aa_478409232_478410380"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_014858855_1_Asiatic_toad_toad_JABFFT010000493_1_246450_247337"
gene_B <- "GCA_029574335_1_plateau_brown_frog_brown_frog_CM056054_1_4639019_4639915"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_028564925_1_wood_frog_frog_CM053031_1_288072376_288073260"
gene_B <- "GCA_019857665_1_Puerto_Rican_coqui_Rican_coqui_WNTK01003886_1_30253_31146"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_027789765_1_hourglass_treefrog_treefrog_CM051064_1_119024895_119025962"
gene_B <- "GCA_905171765_1_common_toad_toad_LR991673_1_10119440_10120324"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_027789765_1_hourglass_treefrog_treefrog_CM051064_1_119125250_119126254"
gene_B <- "GCA_028564925_1_wood_frog_frog_CM053031_1_288004220_288005101"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_014858855_1_Asiatic_toad_toad_CM026472_1_200464445_200465344"
gene_B <- "GCA_029499605_1_Sardinian_treefrog_treefrog_CM056042_1_234558251_234559153"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

##problem:"GCA_009667805_1_Leishan_spiny_toad_spiny_toad_CM019069_1_157468637_157469539"

gene_A <- "GCA_000004195_4_tropical_clawed_frog_clawed_frog_CM004451_2_12114029_12114970"
gene_B <- "GCA_019447015_1_Congo_dwarf_clawed_frog_dwarf_clawed_frog_JAACNH010000905_1_95401_96450"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_000004195_4_tropical_clawed_frog_clawed_frog_CM004451_2_12071867_12072754"
gene_B <- "GCA_017654675_1_African_clawed_frog_clawed_frog_CM030357_1_109598598_109599515"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_000004195_4_tropical_clawed_frog_clawed_frog_CM004451_2_12033009_12033881"
gene_B <- "GCA_019447015_1_Congo_dwarf_clawed_frog_dwarf_clawed_frog_CM033477_1_11261228_11262127"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

##problem GCA_019447015_1_Congo_dwarf_clawed_frog_dwarf_clawed_frog_CM033477_1_3291890_3292795  ##THIS ONE ISN'T WORKING
gene_A <- "GCA_014858855_1_Asiatic_toad_toad_CM026472_1_199891061_199891975"
gene_B <- "GCA_905171775_1_common_frog_frog_LR991685_1_3497372_3498307"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_028564925_1_wood_frog_frog_CM053031_1_293387975_293388916"
gene_B <- "GCA_029206835_1_mountain_yellow_legged_frog_yellow_legged_frog_mountain_yellow_legged_frog_12_aa_183455658_183456587"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_014858855_1_Asiatic_toad_toad_CM026472_1_199973405_199974313"
gene_B <- "GCA_905171775_1_common_frog_frog_LR991685_1_5062476_5063402"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_019512145_1_Tungara_frog_frog_WNYA01001295_1_41610_42593"
gene_B <- "GCA_004786255_1_frogs_and_toads_and_toads_CM016422_1_1164568_1165545"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_028564925_1_wood_frog_frog_CM053031_1_288277407_288278321"
gene_B <- "GCA_004786255_1_frogs_and_toads_and_toads_CM016422_1_1273359_1274270"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_014858855_1_Asiatic_toad_toad_JABFFT010000493_1_224131_225045"
gene_B <- "GCA_004786255_1_frogs_and_toads_and_toads_CM016422_1_3115181_3116104"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#potential problem -- only 2
gene_A <- "GCA_004786255_1_frogs_and_toads_and_toads_CM016422_1_2987530_2988444"
gene_B <- "GCA_027789765_1_hourglass_treefrog_treefrog_CM051064_1_119070365_119071303"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)


gene_A <- "GCA_009667805_1_Leishan_spiny_toad_spiny_toad_CM019069_1_157406932_157407786"
gene_B <- "GCA_004786255_1_frogs_and_toads_and_toads_CM016422_1_1259287_1260297"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_014858855_1_Asiatic_toad_toad_JABFFT010000493_1_201761_202732"
gene_B <- "GCA_019512145_1_Tungara_frog_frog_WNYA01000455_1_91317_92201"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

# gene_A <- "GCA_014858855_1_Asiatic_toad_toad_JABFFT010000493_1_147241_148140"
# gene_B <- "GCA_004786255_1_frogs_and_toads_and_toads_CM016422_1_1499590_1500507"
# result_message <- analyze_genes(gene_A, gene_B, taxon)
# print(result_message)

#Possible problem -- expanded to include several pairs of singletons
gene_A <- "GCA_014858855_1_Asiatic_toad_toad_JABFFT010000493_1_147241_148140"
gene_B <- "GCA_019512145_1_Tungara_frog_frog_WNYA01000455_1_96365_97297"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_027358695_2_plains_spadefoot_toad_spadefoot_toad_CM049776_1_46912765_46913673"
gene_B <- "GCA_027579735_1_fire_bellied_toad_toad_fire_bellied_toad_11_aa_14128821_14129765"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_019512145_1_Tungara_frog_frog_WNYA01001295_1_57325_58245"
gene_B <- "GCA_028390025_1_corroboree_frog_frog_corroboree_frog_7_aa_479704408_479705328"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_004786255_1_frogs_and_toads_and_toads_CM016422_1_1503205_1504110"
gene_B <- "GCA_009667805_1_Leishan_spiny_toad_spiny_toad_CM019069_1_159683860_159684813"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#Also very big, but I think it makes sense
gene_A <- "GCA_014858855_1_Asiatic_toad_toad_JABFFT010000493_1_77970_78896"
gene_B <- "GCA_019447015_1_Congo_dwarf_clawed_frog_dwarf_clawed_frog_CM033477_1_3262881_3263831"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_905171765_1_common_toad_toad_LR991673_1_2955112_2956014"
gene_B <- "GCA_029574335_1_plateau_brown_frog_brown_frog_CM056054_1_2669045_2669953"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_905171775_1_common_frog_frog_LR991685_1_3625707_3626639"
gene_B <- "GCA_019857665_1_Puerto_Rican_coqui_Rican_coqui_WNTK01006178_1_23745_24662"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

##FLAG
gene_A <- "GCA_014858855_1_Asiatic_toad_toad_CM026472_1_200138266_200139198"
gene_B <- "GCA_029206835_1_mountain_yellow_legged_frog_yellow_legged_frog_mountain_yellow_legged_frog_12_aa_183094486_183095424" #GCA_028390025_1_corroboree_frog_frog_corroboree_frog_7_aa_483032947_483033879"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_028564925_1_wood_frog_frog_CM053031_1_297339281_297340210"
gene_B <- "GCA_027789765_1_hourglass_treefrog_treefrog_CM051064_1_118943985_118944896"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_028564925_1_wood_frog_frog_CM053031_1_297380450_297381469"
gene_B <- "GCA_004786255_1_frogs_and_toads_and_toads_CM016422_1_3137570_3138487"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#f
gene_A <- "GCA_019512145_1_Tungara_frog_frog_WNYA01000075_1_258816_259739"
gene_B <- "GCA_028390025_1_corroboree_frog_frog_corroboree_frog_7_aa_482660567_482661517"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_028564925_1_wood_frog_frog_CM053031_1_288707572_288708495"
gene_B <- "GCA_905171775_1_common_frog_frog_LR991685_1_5070276_5071205"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#Includes one weird gene that may actually be outside
gene_A <- "GCA_014858855_1_Asiatic_toad_toad_CM026472_1_199884165_199885133"
gene_B <- "GCA_028564925_1_wood_frog_frog_CM053031_1_292543538_292544449"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)


gene_A <- "GCA_017654675_1_African_clawed_frog_clawed_frog_CM030357_1_109642933_109643850"
gene_B <- "GCA_018994145_1_Yunnan_mustache_toad_mustache_toad_CM032244_1_6421194_6422105"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_027358695_2_plains_spadefoot_toad_spadefoot_toad_CM049776_1_46852584_46853462"
gene_B <- "GCA_027358695_2_plains_spadefoot_toad_spadefoot_toad_CM049776_1_46829978_46830883"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_014858855_1_Asiatic_toad_toad_JABFFT010000493_1_65534_66478"
gene_B <- "GCA_004786255_1_frogs_and_toads_and_toads_CM016422_1_1531378_1532340"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_027579735_1_fire_bellied_toad_toad_fire_bellied_toad_11_aa_14211401_14212345"
gene_B <- "GCA_027358695_2_plains_spadefoot_toad_spadefoot_toad_CM049776_1_46816763_46817680"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#Just two -- potential problem
gene_A <- "GCA_009667805_1_Leishan_spiny_toad_spiny_toad_CM019069_1_159519186_159520223"
gene_B <- "GCA_018994145_1_Yunnan_mustache_toad_mustache_toad_CM032244_1_6427951_6428988"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#Salamanders only -- potential problem
gene_A <- "GCA_002915635_3_axolotl_axolotl_12_ac_28158667_28159563"
gene_B <- "GCA_002915635_3_axolotl_axolotl_5_aa_130217360_130218265"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#Salamanders only -- potential problem
gene_A <- "GCA_002915635_3_axolotl_axolotl_5_aa_123204605_123205531"
gene_B <- "GCA_026652325_1_Iberian_ribbed_newt_ribbed_newt_Iberian_ribbed_newt_11_aa_74937123_74938034"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#Salamanders only -- potential problem
gene_A <- "GCA_002915635_3_axolotl_axolotl_5_aa_130194722_130195657"
gene_B <- "GCA_026652325_1_Iberian_ribbed_newt_ribbed_newt_Iberian_ribbed_newt_11_aa_75100509_75101444"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#Salamanders only -- potential problem
gene_A <- "GCA_002915635_3_axolotl_axolotl_7_ac_134675603_134676664"
gene_B <- "GCA_026652325_1_Iberian_ribbed_newt_ribbed_newt_Iberian_ribbed_newt_9_ac_2505310_2506368"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#Salamanders only -- potential problem
gene_A <- "GCA_026652325_1_Iberian_ribbed_newt_ribbed_newt_Iberian_ribbed_newt_11_ab_73185902_73186831"
gene_B <- "GCA_002915635_3_axolotl_axolotl_5_ab_181376179_181377093"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#Salamanders only -- potential problem
gene_A <- "GCA_002915635_3_axolotl_axolotl_5_ab_236676538_236677443"
gene_B <- "GCA_026652325_1_Iberian_ribbed_newt_ribbed_newt_Iberian_ribbed_newt_11_ab_76991863_76992831"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_014858855_1_Asiatic_toad_toad_CM026465_1_501774788_501775762"
gene_B <- "GCA_027579735_1_fire_bellied_toad_toad_fire_bellied_toad_2_aa_506862050_506863039"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_014858855_1_Asiatic_toad_toad_CM026465_1_501705824_501706864"
gene_B <- "GCA_905171775_1_common_frog_frog_LR991680_1_249689937_249690923"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_019512145_1_Tungara_frog_frog_CM033641_1_218024251_218025219"
gene_B <- "GCA_004786255_1_frogs_and_toads_and_toads_CM016421_1_37239288_37240241"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_014858855_1_Asiatic_toad_toad_CM026465_1_501676332_501677324"
gene_B <- "GCA_027789765_1_hourglass_treefrog_treefrog_CM051058_1_41350507_41351487"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_019512145_1_Tungara_frog_frog_CM033641_1_217994430_217995401"
gene_B <- "GCA_029206835_1_mountain_yellow_legged_frog_yellow_legged_frog_mountain_yellow_legged_frog_1_ab_300503055_300504068"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_027358695_2_plains_spadefoot_toad_spadefoot_toad_CM049770_1_110107707_110108693"
gene_B <- "GCA_018994145_1_Yunnan_mustache_toad_mustache_toad_CM032245_1_210280893_210281870"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_014858855_1_Asiatic_toad_toad_CM026465_1_501747465_501748433"
gene_B <- "GCA_029499605_1_Sardinian_treefrog_treefrog_CM056035_1_403500750_403501700"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_027917425_1_eastern_narrow_mouthed_toad_narrow_mouthed_toad_CM051234_1_428353369_428354346"
gene_B <- "GCA_004786255_1_frogs_and_toads_and_toads_CM016421_1_37234200_37235180"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

##problem, not in right location GCA_027579735_1_fire_bellied_toad_toad_fire_bellied_toad_2_aa_506906501_506907481

##caec only -- potential problem
gene_A <- "GCA_901001135_2_two_lined_caecilian_caecilian_CAAJIB020000526_1_14492_15472"
gene_B <- "GCA_901001135_2_two_lined_caecilian_caecilian_LR584402_1_9669013_9670002"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

##caec only -- potential problem
gene_A <- "GCA_901001135_2_two_lined_caecilian_caecilian_LR584402_1_9646915_9647859"
gene_B <- "GCA_901765095_2_caecilians_LR594645_1_26945060_26946016"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_014858855_1_Asiatic_toad_toad_CM026465_1_501844240_501845199"
gene_B <- "GCA_018994145_1_Yunnan_mustache_toad_mustache_toad_CM032245_1_210123304_210124386"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_009667805_1_Leishan_spiny_toad_spiny_toad_CM019063_1_198437920_198439026"
gene_B <- "GCA_019512145_1_Tungara_frog_frog_CM033641_1_218084856_218085914"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)


#single species -- potential problem
gene_A <- "GCA_027579735_1_fire_bellied_toad_toad_fire_bellied_toad_2_aa_506743741_506744730"
gene_B <- "GCA_027579735_1_fire_bellied_toad_toad_fire_bellied_toad_2_aa_506803443_506804459"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#problem GCA_027579735_1_fire_bellied_toad_toad_fire_bellied_toad_2_aa_506851924_506852922

gene_A <- "GCA_905171765_1_common_toad_toad_LR991668_1_296633933_296634853"
gene_B <- "GCA_004786255_1_frogs_and_toads_and_toads_CM016421_1_37192424_37193335"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)


#only salamanders -- potential problem
gene_A <- "GCA_002915635_3_axolotl_axolotl_25_aa_14872627_14873721"
gene_B <- "GCA_026652325_1_Iberian_ribbed_newt_ribbed_newt_Iberian_ribbed_newt_10_ab_198446393_198447358"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#only caec -- potential problem
gene_A <- "GCA_902459505_2_caecilians_LR699150_1_110059228_110060145"
gene_B <- "GCA_901001135_2_two_lined_caecilian_caecilian_LR584392_1_250834798_250835712"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)


gene_A <- "GCA_002915635_3_axolotl_axolotl_9_aa_176215663_176216577"
gene_B <- "GCA_026652325_1_Iberian_ribbed_newt_ribbed_newt_Iberian_ribbed_newt_3_aa_141457559_141458461"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_026652325_1_Iberian_ribbed_newt_ribbed_newt_Iberian_ribbed_newt_11_ab_84168128_84169057"
gene_B <- "GCA_002915635_3_axolotl_axolotl_28_aa_382392241_382393167"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#single species -- potential problem
gene_A <- "GCA_002915635_3_axolotl_axolotl_10_ab_165971549_165972451"
gene_B <- "GCA_002915635_3_axolotl_axolotl_10_ab_168048163_168049302"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_002915635_3_axolotl_axolotl_19_ab_184081965_184082867"
gene_B <- "GCA_026652325_1_Iberian_ribbed_newt_ribbed_newt_Iberian_ribbed_newt_15_aa_49305218_49306231"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#single species -- potential problem
gene_A <- "GCA_026652325_1_Iberian_ribbed_newt_ribbed_newt_Iberian_ribbed_newt_15_aa_49826458_49827360"
gene_B <- "GCA_026652325_1_Iberian_ribbed_newt_ribbed_newt_Iberian_ribbed_newt_15_aa_53396192_53397094"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)


#problem GCA_026652325_1_Iberian_ribbed_newt_ribbed_newt_Iberian_ribbed_newt_15_aa_49454709_49455611
#problem GCA_002915635_3_axolotl_axolotl_19_ab_183909210_183910250

gene_A <- "GCA_026652325_1_Iberian_ribbed_newt_ribbed_newt_Iberian_ribbed_newt_15_aa_49952580_49953842"
gene_B <- "GCA_002915635_3_axolotl_axolotl_19_ab_187740161_187740937"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)


#single species -- potential problem
gene_A <- "GCA_026652325_1_Iberian_ribbed_newt_ribbed_newt_Iberian_ribbed_newt_15_aa_49755686_49756666"
gene_B <- "GCA_026652325_1_Iberian_ribbed_newt_ribbed_newt_Iberian_ribbed_newt_15_aa_49563759_49564646"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#problem GCA_026652325_1_Iberian_ribbed_newt_ribbed_newt_Iberian_ribbed_newt_15_aa_49252583_49253467
#problem GCA_002915635_3_axolotl_axolotl_5_ab_185023144_185024040

#single species -- potential problem
gene_A <- "GCA_026652325_1_Iberian_ribbed_newt_ribbed_newt_Iberian_ribbed_newt_11_ab_44098055_44098960"
gene_B <- "GCA_026652325_1_Iberian_ribbed_newt_ribbed_newt_Iberian_ribbed_newt_11_ab_44213784_44214689"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#problem GCA_026652325_1_Iberian_ribbed_newt_ribbed_newt_Iberian_ribbed_newt_11_ab_57111428_57112330
#problem GCA_026652325_1_Iberian_ribbed_newt_ribbed_newt_Iberian_ribbed_newt_11_ab_56347319_56348236

gene_A <- "GCA_002915635_3_axolotl_axolotl_8_aa_105409569_105410516"
gene_B <- "GCA_026652325_1_Iberian_ribbed_newt_ribbed_newt_Iberian_ribbed_newt_9_ab_241220391_241221335"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#problem GCA_026652325_1_Iberian_ribbed_newt_ribbed_newt_Iberian_ribbed_newt_9_ab_240925263_240926210

#single species -- potential problem
gene_A <- "GCA_026652325_1_Iberian_ribbed_newt_ribbed_newt_Iberian_ribbed_newt_9_ab_241184397_241185341"
gene_B <- "GCA_026652325_1_Iberian_ribbed_newt_ribbed_newt_Iberian_ribbed_newt_9_ab_241021279_241022208"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#single species -- potential problem
gene_A <- "GCA_901765095_2_caecilians_LR594645_1_32907710_32908705"
gene_B <- "GCA_901765095_2_caecilians_LR594645_1_32911825_32912844"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_002915635_3_axolotl_axolotl_13_ab_43114917_43115888"
gene_B <- "GCA_026652325_1_Iberian_ribbed_newt_ribbed_newt_Iberian_ribbed_newt_12_ac_109137125_109138093"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#single species -- potential problem
gene_A <- "GCA_002915635_3_axolotl_axolotl_4_aa_10641732_10642688"
gene_B <- "GCA_002915635_3_axolotl_axolotl_4_aa_11749262_11750161"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_019857665_1_Puerto_Rican_coqui_Rican_coqui_WNTK01009322_1_16912_17889"
gene_B <- "GCA_004786255_1_frogs_and_toads_and_toads_CM016419_1_4652847_4653755"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_019857665_1_Puerto_Rican_coqui_Rican_coqui_WNTK01001922_1_30157_31122"
gene_B <- "GCA_004786255_1_frogs_and_toads_and_toads_CM016419_1_4648961_4649926"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_009667805_1_Leishan_spiny_toad_spiny_toad_CM019066_1_42607281_42608225"
gene_B <- "GCA_027358695_2_plains_spadefoot_toad_spadefoot_toad_CM049772_1_118948146_118949108"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_000004195_4_tropical_clawed_frog_clawed_frog_CM004447_2_159990732_159991658"
gene_B <- "GCA_019447015_1_Congo_dwarf_clawed_frog_dwarf_clawed_frog_CM033473_1_294655621_294656550"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_000004195_4_tropical_clawed_frog_clawed_frog_CM004447_2_159998140_159999093"
gene_B <- "GCA_024363595_1_Kenyan_clawed_frog_clawed_frog_CM044439_1_171783744_171784697"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_017654675_1_African_clawed_frog_clawed_frog_CM030348_1_168704134_168705111"
gene_B <- "GCA_024363595_1_Kenyan_clawed_frog_clawed_frog_CM044439_1_171778260_171779234"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_000004195_4_tropical_clawed_frog_clawed_frog_CM004447_2_160001627_160002565"
gene_B <- "GCA_017654675_1_African_clawed_frog_clawed_frog_CM030348_1_168717532_168718428"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_019512145_1_Tungara_frog_frog_WNYA01004441_1_7429_8337"
gene_B <- "GCA_004786255_1_frogs_and_toads_and_toads_CM016419_1_4666122_4667105"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_014858855_1_Asiatic_toad_toad_CM026468_1_50831813_50832778"
gene_B <- "GCA_004786255_1_frogs_and_toads_and_toads_CM016419_1_4590084_4591019"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_014858855_1_Asiatic_toad_toad_CM026468_1_48907419_48908474"
gene_B <- "GCA_027789765_1_hourglass_treefrog_treefrog_CM051061_1_117974784_117975698"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

# #including this one in below
# gene_A <- "GCA_014858855_1_Asiatic_toad_toad_CM026468_1_49370677_49371612"
# gene_B <- "GCA_905171765_1_common_toad_toad_LR991670_1_38633175_38634110"
# result_message <- analyze_genes(gene_A, gene_B, taxon)
# print(result_message)

#problem, might be part of above GCA_014858855_1_Asiatic_toad_toad_CM026468_1_49451068_49452033


##MASSIVE parallel expansion
gene_A <- "GCA_014858855_1_Asiatic_toad_toad_CM026465_1_411493683_411494777"
gene_B <- "GCA_019512145_1_Tungara_frog_frog_WNYA01000358_1_177003_177905"
result_message <- analyze_genes(gene_A, gene_B, taxon)

#problem GCA_019512145_1_Tungara_frog_frog_WNYA01000451_1_108650_109807
print(result_message)

#single species -- potential problem
gene_A <- "GCA_019857665_1_Puerto_Rican_coqui_Rican_coqui_WNTK01000826_1_115150_116055"
gene_B <- "GCA_019857665_1_Puerto_Rican_coqui_Rican_coqui_WNTK01003773_1_41112_42035"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_014858855_1_Asiatic_toad_toad_CM026468_1_50175202_50176017"
gene_B <- "GCA_029499605_1_Sardinian_treefrog_treefrog_CM056037_1_9726159_9727109"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_019512145_1_Tungara_frog_frog_WNYA01000451_1_26835_27887"
gene_B <- "GCA_027789765_1_hourglass_treefrog_treefrog_CM051061_1_118347338_118348279"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_019512145_1_Tungara_frog_frog_WNYA01005051_1_13423_14373"
gene_B <- "GCA_027789765_1_hourglass_treefrog_treefrog_CM051061_1_117263910_117264968"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#problem GCA_019512145_1_Tungara_frog_frog_WNYA01033205_1_678_1643

gene_A <- "GCA_027917425_1_eastern_narrow_mouthed_toad_narrow_mouthed_toad_CM051237_1_56520477_56521424"
gene_B <- "GCA_004786255_1_frogs_and_toads_and_toads_CM016419_1_5255801_5256964"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_028564925_1_wood_frog_frog_CM053028_1_518433921_518434958"
gene_B <- "GCA_004786255_1_frogs_and_toads_and_toads_CM016419_1_5836403_5837362"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_004786255_1_frogs_and_toads_and_toads_CM016419_1_5891078_5891944"
gene_B <- "GCA_027917425_1_eastern_narrow_mouthed_toad_narrow_mouthed_toad_CM051236_1_4136977_4137924"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#single species -- possible problem
gene_A <- "GCA_028390025_1_corroboree_frog_frog_corroboree_frog_4_aa_45999808_46000740"
gene_B <- "GCA_028390025_1_corroboree_frog_frog_corroboree_frog_4_aa_46601863_46602774"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_028564925_1_wood_frog_frog_CM053028_1_533139520_533140401"
gene_B <- "GCA_029206835_1_mountain_yellow_legged_frog_yellow_legged_frog_JARGYO010000094_1_62793198_62794097"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_028564925_1_wood_frog_frog_CM053028_1_533465719_533466600"
gene_B <- "GCA_004786255_1_frogs_and_toads_and_toads_CM016419_1_4641894_4642823"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_019857665_1_Puerto_Rican_coqui_Rican_coqui_WNTK01004753_1_17642_18583"
gene_B <- "GCA_027789765_1_hourglass_treefrog_treefrog_CM051061_1_118341205_118342257"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_028564925_1_wood_frog_frog_CM053028_1_533385137_533386042"
gene_B <- "GCA_004786255_1_frogs_and_toads_and_toads_CM016419_1_4608273_4609232"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_028564925_1_wood_frog_frog_CM053028_1_533451421_533452341"
gene_B <- "GCA_905171765_1_common_toad_toad_LR991670_1_38824396_38825346"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_029499605_1_Sardinian_treefrog_treefrog_CM056037_1_10654276_10655190"
gene_B <- "GCA_028390025_1_corroboree_frog_frog_corroboree_frog_4_aa_46790330_46791229"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_027358695_2_plains_spadefoot_toad_spadefoot_toad_CM049772_1_119024290_119025216"
gene_B <- "GCA_018994145_1_Yunnan_mustache_toad_mustache_toad_CM032237_1_52475439_52476362"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_000004195_4_tropical_clawed_frog_clawed_frog_CM004447_2_159983992_159984951"
gene_B <- "GCA_027579735_1_fire_bellied_toad_toad_fire_bellied_toad_4_aa_69315844_69316761"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#problem GCA_901765095_2_caecilians_LR594645_1_32902235_32903170
#problem GCA_002915635_3_axolotl_axolotl_7_aa_130525418_130526398

gene_A <- "GCA_014858855_1_Asiatic_toad_toad_CM026470_1_290788996_290790012"
gene_B <- "GCA_027579735_1_fire_bellied_toad_toad_fire_bellied_toad_9_aa_266394810_266395808"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_009667805_1_Leishan_spiny_toad_spiny_toad_CM019068_1_52313721_52314719"
gene_B <- "GCA_027358695_2_plains_spadefoot_toad_spadefoot_toad_CM049780_1_21057744_21058685"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_009667805_1_Leishan_spiny_toad_spiny_toad_CM019068_1_3909329_3910390"
gene_B <- "GCA_018994145_1_Yunnan_mustache_toad_mustache_toad_CM032241_1_185994387_185995448"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_000004195_4_tropical_clawed_frog_clawed_frog_CM004449_2_18861902_18862864"
gene_B <- "GCA_019447015_1_Congo_dwarf_clawed_frog_dwarf_clawed_frog_CM033475_1_18282151_18283119"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#problem GCA_019447015_1_Congo_dwarf_clawed_frog_dwarf_clawed_frog_CM033475_1_18341697_18342650

gene_A <- "GCA_000004195_4_tropical_clawed_frog_clawed_frog_CM004449_2_18855735_18856757"
gene_B <- "GCA_019447015_1_Congo_dwarf_clawed_frog_dwarf_clawed_frog_CM033475_1_18051542_18052519"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_002915635_3_axolotl_PGSH02019610_1_43444_44397"
gene_B <- "GCA_026652325_1_Iberian_ribbed_newt_ribbed_newt_Iberian_ribbed_newt_11_aa_12861746_12862699"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#problem GCA_026652325_1_Iberian_ribbed_newt_ribbed_newt_Iberian_ribbed_newt_11_aa_12777811_12778770
#problem GCA_901765095_2_caecilians_LR594645_1_33379835_33380842

gene_A <- "GCA_014858855_1_Asiatic_toad_toad_CM026468_1_48668754_48669653"
gene_B <- "GCA_029499605_1_Sardinian_treefrog_treefrog_CM056037_1_12081763_12082731"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_009667805_1_Leishan_spiny_toad_spiny_toad_CM019066_1_42699650_42700687"
gene_B <- "GCA_027579735_1_fire_bellied_toad_toad_fire_bellied_toad_4_aa_68353144_68354082"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_027579735_1_fire_bellied_toad_toad_fire_bellied_toad_4_aa_67674339_67675271"
gene_B <- "GCA_027579735_1_fire_bellied_toad_toad_fire_bellied_toad_4_aa_67739371_67740303"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_027579735_1_fire_bellied_toad_toad_fire_bellied_toad_4_aa_68246215_68247129"
gene_B <- "GCA_027579735_1_fire_bellied_toad_toad_fire_bellied_toad_4_aa_68293917_68294825"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_027579735_1_fire_bellied_toad_toad_fire_bellied_toad_4_aa_69038483_69039415"
gene_B <- "GCA_027579735_1_fire_bellied_toad_toad_fire_bellied_toad_4_aa_69139141_69140106"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_029499605_1_Sardinian_treefrog_treefrog_CM056037_1_12065283_12066188"
gene_B <- "GCA_019857665_1_Puerto_Rican_coqui_Rican_coqui_WNTK01000954_1_96238_97170"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_014858855_1_Asiatic_toad_toad_CM026468_1_48761674_48762642"
gene_B <- "GCA_027789765_1_hourglass_treefrog_treefrog_CM051061_1_118526344_118527393"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_014858855_1_Asiatic_toad_toad_CM026468_1_48741373_48742293"
gene_B <- "GCA_019512145_1_Tungara_frog_frog_WNYA01000090_1_218377_219297"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_019857665_1_Puerto_Rican_coqui_Rican_coqui_WNTK01001005_1_101195_102121"
gene_B <- "GCA_027789765_1_hourglass_treefrog_treefrog_CM051061_1_118389357_118390475"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_014858855_1_Asiatic_toad_toad_CM026468_1_48715105_48716013"
gene_B <- "GCA_027789765_1_hourglass_treefrog_treefrog_CM051061_1_118806406_118807566"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_029499605_1_Sardinian_treefrog_treefrog_CM056037_1_10417645_10418556"
gene_B <- "GCA_019857665_1_Puerto_Rican_coqui_Rican_coqui_WNTK01000954_1_42872_43900"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_028564925_1_wood_frog_frog_CM053028_1_533666591_533667514"
gene_B <- "GCA_029206835_1_mountain_yellow_legged_frog_yellow_legged_frog_JARGYO010000094_1_61800463_61801395"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#problem: GCA_027917425_1_eastern_narrow_mouthed_toad_narrow_mouthed_toad_CM051237_1_56560849_56561766

gene_A <- "GCA_028564925_1_wood_frog_frog_CM053028_1_533646387_533647310"
gene_B <- "GCA_029206835_1_mountain_yellow_legged_frog_yellow_legged_frog_JARGYO010000094_1_61905725_61906744"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_028564925_1_wood_frog_frog_CM053028_1_533806870_533807799"
gene_B <- "GCA_029206835_1_mountain_yellow_legged_frog_yellow_legged_frog_JARGYO010000094_1_61611963_61613027"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_905171775_1_common_frog_frog_LR991683_1_27656291_27657214"
gene_B <- "GCA_004786255_1_frogs_and_toads_and_toads_CM016419_1_1075222_1076127"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#single species -- potential problem
gene_A <- "GCA_027579735_1_fire_bellied_toad_toad_fire_bellied_toad_4_aa_68897684_68898631"
gene_B <- "GCA_027579735_1_fire_bellied_toad_toad_fire_bellied_toad_4_aa_68971485_68972432"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#single species -- potential problem
gene_A <- "GCA_027579735_1_fire_bellied_toad_toad_fire_bellied_toad_4_aa_67866381_67867442"
gene_B <- "GCA_027579735_1_fire_bellied_toad_toad_fire_bellied_toad_4_aa_68148228_68149151"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_000004195_4_tropical_clawed_frog_clawed_frog_CM004447_2_160040586_160041506"
gene_B <- "GCA_024363595_1_Kenyan_clawed_frog_clawed_frog_CM044440_1_136020781_136021701"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_000004195_4_tropical_clawed_frog_clawed_frog_CM004447_2_160048579_160049499"
gene_B <- "GCA_024363595_1_Kenyan_clawed_frog_clawed_frog_CM044439_1_171890563_171891477"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_000004195_4_tropical_clawed_frog_clawed_frog_CM004447_2_160037873_160038811"
gene_B <- "GCA_024363595_1_Kenyan_clawed_frog_clawed_frog_CM044439_1_171861105_171861947"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_019447015_1_Congo_dwarf_clawed_frog_dwarf_clawed_frog_CM033473_1_293731646_293732545"
gene_B <- "GCA_019447015_1_Congo_dwarf_clawed_frog_dwarf_clawed_frog_CM033473_1_294366810_294367841"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_000004195_4_tropical_clawed_frog_clawed_frog_CM004447_2_160031157_160032143"
gene_B <- "GCA_024363595_1_Kenyan_clawed_frog_clawed_frog_CM044440_1_136006536_136007471"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_028564925_1_wood_frog_frog_CM053028_1_532128171_532129058"
gene_B <- "GCA_004786255_1_frogs_and_toads_and_toads_CM016419_1_1100455_1101363"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_019857665_1_Puerto_Rican_coqui_Rican_coqui_WNTK01002590_1_27144_28073"
gene_B <- "GCA_027789765_1_hourglass_treefrog_treefrog_CM051061_1_118159956_118160861"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_905171765_1_common_toad_toad_LR991670_1_38777756_38778661"
gene_B <- "GCA_019857665_1_Puerto_Rican_coqui_Rican_coqui_WNTK01002590_1_71085_72053"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_019857665_1_Puerto_Rican_coqui_Rican_coqui_WNTK01002464_1_56079_57044"
gene_B <- "GCA_027789765_1_hourglass_treefrog_treefrog_CM051061_1_118051561_118052511"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#could have made this part of above
gene_A <- "GCA_027789765_1_hourglass_treefrog_treefrog_CM051061_1_118072427_118073338"
gene_B <- "GCA_027789765_1_hourglass_treefrog_treefrog_CM051061_1_118093907_118095121"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_028564925_1_wood_frog_frog_CM053028_1_533583160_533584062"
gene_B <- "GCA_004786255_1_frogs_and_toads_and_toads_CM016419_1_4671831_4672721"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#single species -- potential problem
gene_A <- "GCA_028390025_1_corroboree_frog_frog_corroboree_frog_4_aa_46931631_46932563"
gene_B <- "GCA_028390025_1_corroboree_frog_frog_corroboree_frog_4_aa_47679119_47679961"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_009667805_1_Leishan_spiny_toad_spiny_toad_CM019066_1_42665094_42666020"
gene_B <- "GCA_027358695_2_plains_spadefoot_toad_spadefoot_toad_CM049772_1_118925963_118926865"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#single species -- potential problem
gene_A <- "GCA_027579735_1_fire_bellied_toad_toad_fire_bellied_toad_4_aa_66286378_66287310"
gene_B <- "GCA_027579735_1_fire_bellied_toad_toad_fire_bellied_toad_4_aa_66463336_66464247"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_000004195_4_tropical_clawed_frog_clawed_frog_CM004447_2_160025846_160026766"
gene_B <- "GCA_019447015_1_Congo_dwarf_clawed_frog_dwarf_clawed_frog_CM033473_1_294472669_294473541"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_000004195_4_tropical_clawed_frog_clawed_frog_CM004447_2_160012672_160013574"
gene_B <- "GCA_019447015_1_Congo_dwarf_clawed_frog_dwarf_clawed_frog_CM033473_1_294625503_294626369"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)


gene_A <- "GCA_902459505_2_caecilians_LR699146_1_442155494_442156504"
gene_B <- "GCA_901001135_2_two_lined_caecilian_caecilian_LR584387_1_496819530_496820549"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#single species -- potential problem
gene_A <- "GCA_902459505_2_caecilians_LR699146_1_448014052_448015044"
gene_B <- "GCA_901765095_2_caecilians_LR594633_1_282897510_282898532"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_901001135_2_two_lined_caecilian_caecilian_LR584387_1_496481113_496482150"
gene_B <- "GCA_901765095_2_caecilians_LR594633_1_282884006_282885061"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_026652325_1_Iberian_ribbed_newt_ribbed_newt_Iberian_ribbed_newt_9_ad_244232043_244232969"
gene_B <- "GCA_002915635_3_axolotl_axolotl_7_aa_150835345_150836511"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_026652325_1_Iberian_ribbed_newt_ribbed_newt_Iberian_ribbed_newt_11_aa_74562572_74563492"
gene_B <- "GCA_002915635_3_axolotl_axolotl_16_ab_324466481_324467404"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_002915635_3_axolotl_axolotl_5_aa_122623391_122624509"
gene_B <- "GCA_026652325_1_Iberian_ribbed_newt_ribbed_newt_Iberian_ribbed_newt_11_aa_70239458_70240411"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_027789765_1_hourglass_treefrog_treefrog_CM051063_1_42518061_42518990"
gene_B <- "GCA_019512145_1_Tungara_frog_frog_WNYA01000589_1_34770_35717"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_014858855_1_Asiatic_toad_toad_CM026470_1_291456606_291457550"
gene_B <- "GCA_029499605_1_Sardinian_treefrog_treefrog_CM056041_1_79639357_79640301"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_014858855_1_Asiatic_toad_toad_CM026470_1_291458904_291459848"
gene_B <- "GCA_027789765_1_hourglass_treefrog_treefrog_CM051063_1_42503256_42504200"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#single species -- potential problem
gene_A <- "GCA_028390025_1_corroboree_frog_frog_corroboree_frog_3_ab_140310928_140311839"
gene_B <- "GCA_028390025_1_corroboree_frog_frog_corroboree_frog_3_ab_140410252_140411184"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_028564925_1_wood_frog_frog_CM053036_1_88317254_88318198"
gene_B <- "GCA_004786255_1_frogs_and_toads_and_toads_CM016425_1_26923478_26924425"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_028564925_1_wood_frog_frog_CM053036_1_88305452_88306399"
gene_B <- "GCA_027917425_1_eastern_narrow_mouthed_toad_narrow_mouthed_toad_CM051235_1_596636938_596637885"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_027917425_1_eastern_narrow_mouthed_toad_narrow_mouthed_toad_CM051235_1_596644011_596644949"
gene_B <- "GCA_027789765_1_hourglass_treefrog_treefrog_CM051063_1_42474524_42475471"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_019512145_1_Tungara_frog_frog_WNYA01000589_1_87391_88335"
gene_B <- "GCA_004786255_1_frogs_and_toads_and_toads_CM016425_1_26937965_26938909"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_014858855_1_Asiatic_toad_toad_CM026470_1_291604352_291605299"
gene_B <- "GCA_019512145_1_Tungara_frog_frog_WNYA01012115_1_458_1405"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_014858855_1_Asiatic_toad_toad_CM026470_1_291631939_291632865"
gene_B <- "GCA_028390025_1_corroboree_frog_frog_corroboree_frog_3_ab_140360927_140361868"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_028564925_1_wood_frog_frog_CM053036_1_88370349_88371287"
gene_B <- "GCA_004786255_1_frogs_and_toads_and_toads_CM016425_1_26932958_26933902"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_027917425_1_eastern_narrow_mouthed_toad_narrow_mouthed_toad_CM051235_1_596548803_596549750"
gene_B <- "GCA_004786255_1_frogs_and_toads_and_toads_CM016425_1_26928686_26929633"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_009667805_1_Leishan_spiny_toad_spiny_toad_CM019068_1_52659980_52660933"
gene_B <- "GCA_027358695_2_plains_spadefoot_toad_spadefoot_toad_CM049780_1_21146214_21147158"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#problem GCA_027579735_1_fire_bellied_toad_toad_fire_bellied_toad_9_aa_265312167_265313102
#problem GCA_027579735_1_fire_bellied_toad_toad_fire_bellied_toad_9_aa_265375768_265376742

#very old!  Consider splitting into two
gene_A <- "GCA_000004195_4_tropical_clawed_frog_clawed_frog_CM004449_2_18691391_18692386"
gene_B <- "GCA_027358695_2_plains_spadefoot_toad_spadefoot_toad_CM049780_1_21155392_21156324"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#single species -- potential problem
gene_A <- "GCA_002915635_3_axolotl_axolotl_17_aa_215545311_215546222"
gene_B <- "GCA_002915635_3_axolotl_axolotl_17_aa_215871481_215872389"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#problem GCA_901765095_2_caecilians_LR594645_1_33203963_33204889

gene_A <- "GCA_002915635_3_axolotl_axolotl_5_aa_59110188_59111207"
gene_B <- "GCA_026652325_1_Iberian_ribbed_newt_ribbed_newt_Iberian_ribbed_newt_11_aa_33378386_33379357"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_002915635_3_axolotl_axolotl_5_aa_50920824_50921777"
gene_B <- "GCA_026652325_1_Iberian_ribbed_newt_ribbed_newt_Iberian_ribbed_newt_11_aa_29126305_29127264"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_026652325_1_Iberian_ribbed_newt_ribbed_newt_Iberian_ribbed_newt_11_aa_25570890_25571855"
gene_B <- "GCA_002915635_3_axolotl_axolotl_5_aa_48313896_48314858"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#single species -- potential problem
gene_A <- "GCA_002915635_3_axolotl_axolotl_5_aa_34847698_34848648"
gene_B <- "GCA_002915635_3_axolotl_axolotl_5_aa_34880522_34881391"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#single species -- potential problem
gene_A <- "GCA_026652325_1_Iberian_ribbed_newt_ribbed_newt_Iberian_ribbed_newt_11_aa_29186694_29187647"
gene_B <- "GCA_026652325_1_Iberian_ribbed_newt_ribbed_newt_Iberian_ribbed_newt_11_aa_25534020_25534979"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_002915635_3_axolotl_axolotl_5_aa_20861205_20862152"
gene_B <- "GCA_026652325_1_Iberian_ribbed_newt_ribbed_newt_Iberian_ribbed_newt_11_aa_17328976_17329938"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_026652325_1_Iberian_ribbed_newt_ribbed_newt_Iberian_ribbed_newt_5_ab_499568518_499569465"
gene_B <- "GCA_002915635_3_axolotl_axolotl_17_aa_226525863_226526771"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_901001135_2_two_lined_caecilian_caecilian_LR584388_1_347605380_347606291"
gene_B <- "GCA_002915635_3_axolotl_axolotl_14_ab_481866378_481867322"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_017654675_1_African_clawed_frog_clawed_frog_CM030352_1_13200406_13201296"
gene_B <- "GCA_019447015_1_Congo_dwarf_clawed_frog_dwarf_clawed_frog_CM033475_1_18872516_18873589"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_009667805_1_Leishan_spiny_toad_spiny_toad_CM019068_1_52745717_52746628"
gene_B <- "GCA_027358695_2_plains_spadefoot_toad_spadefoot_toad_CM049780_1_21170566_21171477"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_009667805_1_Leishan_spiny_toad_spiny_toad_CM019068_1_53052293_53053189"
gene_B <- "GCA_018994145_1_Yunnan_mustache_toad_mustache_toad_CM032241_1_135943033_135943932"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_009667805_1_Leishan_spiny_toad_spiny_toad_CM019068_1_53077212_53078153"
gene_B <- "GCA_027358695_2_plains_spadefoot_toad_spadefoot_toad_CM049780_1_21178852_21179757"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#single species -- potential problem
gene_A <- "GCA_027579735_1_fire_bellied_toad_toad_fire_bellied_toad_9_aa_264998439_264999353"
gene_B <- "GCA_027579735_1_fire_bellied_toad_toad_fire_bellied_toad_9_aa_265014979_265015893"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#single species -- potential problem
gene_A <- "GCA_901765095_2_caecilians_LR594645_1_33300853_33301776"
gene_B <- "GCA_901765095_2_caecilians_LR594645_1_33317670_33318635"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_905171765_1_common_toad_toad_LR991672_1_101326515_101327444"
gene_B <- "GCA_018994145_1_Yunnan_mustache_toad_mustache_toad_CM032241_1_136294727_136295662"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_019512145_1_Tungara_frog_frog_WNYA01000269_1_168355_169248"
gene_B <- "GCA_019857665_1_Puerto_Rican_coqui_Rican_coqui_WNTK01004906_1_487_1383"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_014858855_1_Asiatic_toad_toad_CM026475_1_77199350_77200276"
gene_B <- "GCA_019512145_1_Tungara_frog_frog_WNYA01001056_1_47426_48325"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_019512145_1_Tungara_frog_frog_CM033647_1_11195931_11196848"
gene_B <- "GCA_028390025_1_corroboree_frog_frog_corroboree_frog_12_aa_150298624_150299529"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_009667805_1_Leishan_spiny_toad_spiny_toad_CM019075_1_16515076_16515951"
gene_B <- "GCA_018994145_1_Yunnan_mustache_toad_mustache_toad_CM032240_1_75329931_75330806"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_000004195_4_tropical_clawed_frog_clawed_frog_CM004450_2_121656610_121657698"
gene_B <- "GCA_027358695_2_plains_spadefoot_toad_spadefoot_toad_CM049778_1_20154286_20155182"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)


gene_A <- "GCA_009667805_1_Leishan_spiny_toad_spiny_toad_CM019075_1_16310754_16311641"
gene_B <- "GCA_009667805_1_Leishan_spiny_toad_spiny_toad_CM019075_1_16426907_16427887"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)


#problem GCA_018994145_1_Yunnan_mustache_toad_mustache_toad_CM032240_1_75434058_75434966

gene_A <- "GCA_009667805_1_Leishan_spiny_toad_spiny_toad_CM019075_1_16113975_16114850"#"GCA_009667805_1_Leishan_spiny_toad_spiny_toad_CM019075_1_16281301_16282212"
gene_B <- "GCA_027358695_2_plains_spadefoot_toad_spadefoot_toad_CM049778_1_20090972_20091904"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_009667805_1_Leishan_spiny_toad_spiny_toad_CM019075_1_16257013_16257930"#"GCA_009667805_1_Leishan_spiny_toad_spiny_toad_CM019075_1_16281301_16282212"
gene_B <- "GCA_018994145_1_Yunnan_mustache_toad_mustache_toad_CM032240_1_75612864_75613784"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#single species -- potential problem
gene_A <- "GCA_027358695_2_plains_spadefoot_toad_spadefoot_toad_CM049778_1_20100162_20101094"
gene_B <- "GCA_027358695_2_plains_spadefoot_toad_spadefoot_toad_CM049778_1_20102462_20103352"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_027358695_2_plains_spadefoot_toad_spadefoot_toad_CM049778_1_20112296_20113222"
gene_B <- "GCA_018994145_1_Yunnan_mustache_toad_mustache_toad_CM032240_1_75316538_75317479"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#problem GCA_027358695_2_plains_spadefoot_toad_spadefoot_toad_CM049778_1_20141418_20142293

gene_A <- "GCA_009667805_1_Leishan_spiny_toad_spiny_toad_RXON01001720_1_44015_44953"
gene_B <- "GCA_027358695_2_plains_spadefoot_toad_spadefoot_toad_CM049778_1_20147137_20148033"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_027358695_2_plains_spadefoot_toad_spadefoot_toad_CM049778_1_20138499_20139437"
gene_B <- "GCA_027358695_2_plains_spadefoot_toad_spadefoot_toad_CM049778_1_20151156_20152067"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_009667805_1_Leishan_spiny_toad_spiny_toad_CM019075_1_1883591_1884538"
gene_B <- "GCA_027358695_2_plains_spadefoot_toad_spadefoot_toad_CM049778_1_20109560_20110471"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_019512145_1_Tungara_frog_frog_CM033641_1_144070289_144071185"
gene_B <- "GCA_027789765_1_hourglass_treefrog_treefrog_CM051068_1_20138901_20139998"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_014858855_1_Asiatic_toad_toad_CM026468_1_65996695_65997612"
gene_B <- "GCA_019512145_1_Tungara_frog_frog_WNYA01001056_1_72396_73274"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#problem GCA_019857665_1_Puerto_Rican_coqui_Rican_coqui_WNTK01015483_1_9520_10398

gene_A <- "GCA_014858855_1_Asiatic_toad_toad_CM026475_1_77255350_77256237"
gene_B <- "GCA_027789765_1_hourglass_treefrog_treefrog_CM051068_1_20005168_20006022"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_014858855_1_Asiatic_toad_toad_CM026475_1_83072603_83073481"
gene_B <- "GCA_019512145_1_Tungara_frog_frog_WNYA01016017_1_2657_3535"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#single species -- potential problem
gene_A <- "GCA_019857665_1_Puerto_Rican_coqui_Rican_coqui_WNTK01006831_1_3106_3996"
gene_B <- "GCA_019857665_1_Puerto_Rican_coqui_Rican_coqui_WNTK01007498_1_18374_19267"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#problem GCA_027789765_1_hourglass_treefrog_treefrog_CM051068_1_18771735_18772610
#problem GCA_029499605_1_Sardinian_treefrog_treefrog_CM056045_1_92675711_92676592
#problem GCA_019512145_1_Tungara_frog_frog_WNYA01038552_1_177_1061
#problem GCA_019857665_1_Puerto_Rican_coqui_Rican_coqui_WNTK01003228_1_5779_6666
#problem GCA_019512145_1_Tungara_frog_frog_WNYA01016845_1_3405_4328
#problem GCA_029206835_1_mountain_yellow_legged_frog_yellow_legged_frog_mountain_yellow_legged_frog_11_aa_78556017_78556928

gene_A <- "GCA_905171765_1_common_toad_toad_LR991677_1_83958592_83959518"
gene_B <- "GCA_027789765_1_hourglass_treefrog_treefrog_CM051068_1_18988691_18989602"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_029574335_1_plateau_brown_frog_brown_frog_CM056049_1_728705692_728706597"
gene_B <- "GCA_029499605_1_Sardinian_treefrog_treefrog_CM056045_1_92700812_92701711"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_014858855_1_Asiatic_toad_toad_CM026475_1_77189313_77190266"
gene_B <- "GCA_019512145_1_Tungara_frog_frog_CM033647_1_11205026_11205979"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_019512145_1_Tungara_frog_frog_WNYA01000400_1_132858_133721"
gene_B <- "GCA_019857665_1_Puerto_Rican_coqui_Rican_coqui_WNTK01024859_1_4561_5424"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_014858855_1_Asiatic_toad_toad_CM026475_1_72552146_72553003"
gene_B <- "GCA_027789765_1_hourglass_treefrog_treefrog_CM051068_1_20159432_20160469"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_014858855_1_Asiatic_toad_toad_CM026475_1_77278876_77279763"
gene_B <- "GCA_905171765_1_common_toad_toad_LR991677_1_84205195_84206082"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_019512145_1_Tungara_frog_frog_WNYA01000080_1_19853_20761"
gene_B <- "GCA_027789765_1_hourglass_treefrog_treefrog_CM051068_1_20096217_20097116"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_029499605_1_Sardinian_treefrog_treefrog_CM056045_1_92579295_92580203"
gene_B <- "GCA_019857665_1_Puerto_Rican_coqui_Rican_coqui_WNTK01004536_1_31430_32335"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#problem GCA_019857665_1_Puerto_Rican_coqui_Rican_coqui_WNTK01006652_1_8672_9574
#problem GCA_019512145_1_Tungara_frog_frog_WNYA01000719_1_100327_101295

gene_A <- "GCA_014858855_1_Asiatic_toad_toad_CM026475_1_83049899_83050888"
gene_B <- "GCA_019512145_1_Tungara_frog_frog_CM033647_1_9692480_9693406"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#very old!
gene_A <- "GCA_014858855_1_Asiatic_toad_toad_CM026475_1_78545873_78546781"
gene_B <- "GCA_024363595_1_Kenyan_clawed_frog_clawed_frog_CM044445_1_107329706_107330611"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_028390025_1_corroboree_frog_frog_corroboree_frog_12_aa_149778458_149779450"
gene_B <- "GCA_027358695_2_plains_spadefoot_toad_spadefoot_toad_CM049778_1_20144945_20145934"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_027917425_1_eastern_narrow_mouthed_toad_narrow_mouthed_toad_CM051235_1_401858718_401859569"
gene_B <- "GCA_028390025_1_corroboree_frog_frog_corroboree_frog_12_aa_149561362_149562264"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#problem GCA_009667805_1_Leishan_spiny_toad_spiny_toad_RXON01001720_1_1859_2812"
#problem GCA_019512145_1_Tungara_frog_frog_WNYA01000719_1_97243_98247

gene_A <- "GCA_009667805_1_Leishan_spiny_toad_spiny_toad_CM019068_1_29098136_29099107"
gene_B <- "GCA_027579735_1_fire_bellied_toad_toad_fire_bellied_toad_9_aa_286117346_286118269"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_014858855_1_Asiatic_toad_toad_CM026470_1_362529797_362530768"
gene_B <- "GCA_019857665_1_Puerto_Rican_coqui_Rican_coqui_WNTK01006898_1_16386_17363"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_014858855_1_Asiatic_toad_toad_CM026470_1_362541119_362542075"
gene_B <- "GCA_027789765_1_hourglass_treefrog_treefrog_CM051063_1_15652565_15653491"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#single species -- potential problem
gene_A <- "GCA_028390025_1_corroboree_frog_frog_corroboree_frog_3_ab_251130567_251131520"
gene_B <- "GCA_028390025_1_corroboree_frog_frog_corroboree_frog_3_ab_251252580_251253533"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_028564925_1_wood_frog_frog_CM053036_1_3182453_3183448"
gene_B <- "GCA_905171775_1_common_frog_frog_LR991687_1_1959073_1959960"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_028564925_1_wood_frog_frog_CM053036_1_3223888_3224859"
gene_B <- "GCA_029206835_1_mountain_yellow_legged_frog_yellow_legged_frog_mountain_yellow_legged_frog_8_aa_320538354_320539436"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#problem GCA_004786255_1_frogs_and_toads_and_toads_CM016425_1_9147610_9148563

#out of place  #single species -- potential problem
gene_A <- "GCA_027579735_1_fire_bellied_toad_toad_fire_bellied_toad_11_aa_11839376_11840290"
gene_B <- "GCA_027579735_1_fire_bellied_toad_toad_fire_bellied_toad_11_aa_9255718_9256635"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#out of place
gene_A <- "GCA_009667805_1_Leishan_spiny_toad_spiny_toad_CM019069_1_157878488_157879396"
gene_B <- "GCA_009667805_1_Leishan_spiny_toad_spiny_toad_CM019069_1_157610420_157611370"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#out of place
gene_A <- "GCA_014858855_1_Asiatic_toad_toad_CM026468_1_48797876_48798838"
gene_B <- "GCA_004786255_1_frogs_and_toads_and_toads_CM016419_1_4592676_4593545"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#problem GCA_026652325_1_Iberian_ribbed_newt_ribbed_newt_Iberian_ribbed_newt_9_ad_257110961_257111875

###MAMMAL
taxon <- "mammal"

#1
gene_A <- "GCA_014905855_1_human_AP023467_1_121879574_121880449"
gene_B <- "GCA_004115215_4_platypus_CM014212_1_37671841_37672944"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#2
gene_A <- "GCA_014905855_1_human_AP023467_1_142420492_142421415"
gene_B <- "GCA_027887165_1_gray_short_tailed_opossum_short_tailed_opossum_CM051227_1_218843627_218844574"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#3
gene_A <- "GCA_014905855_1_human_AP023467_1_142386062_142387018"
gene_B <- "GCA_027887165_1_gray_short_tailed_opossum_short_tailed_opossum_CM051227_1_218687429_218688346"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#4
gene_A <- "GCA_029281585_1_western_lowland_gorilla_lowland_gorilla_CM055453_1_150905550_150906449"
gene_B <- "GCA_027887165_1_gray_short_tailed_opossum_short_tailed_opossum_CM051227_1_216741171_216742088"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#5
gene_A <- "GCA_014905855_1_human_AP023467_1_142126152_142127168"
gene_B <- "GCA_016433145_1_Agile_Gracile_Mouse_Opossum_Gracile_Mouse_Opossum_CM028239_1_106471300_106472268"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#6
gene_A <- "GCA_014905855_1_human_AP023467_1_142164799_142165770"
gene_B <- "GCA_027887165_1_gray_short_tailed_opossum_short_tailed_opossum_CM051227_1_218504863_218505825"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#7
gene_A <- "GCA_014905855_1_human_AP023465_1_10241571_10242470"
gene_B <- "GCA_027887165_1_gray_short_tailed_opossum_short_tailed_opossum_CM051225_1_383671900_383672844"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#8
gene_A <- "GCA_015852505_1_Australian_echidna_echidna_CM027596_1_176237397_176238362"
gene_B <- "GCA_004115215_4_platypus_CM014212_1_37537979_37538908"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#9
gene_A <- "GCA_000001635_9_house_mouse_mouse_CM000999_3_13263694_13264605"
gene_B <- "GCA_007646695_3_Tasmanian_wolf_wolf_CM040578_1_292606627_292607538"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#10
gene_A <- "GCA_014905855_1_human_AP023467_1_140853616_140854617"
gene_B <- "GCA_019393635_1_monito_del_monte_del_monte_CM033369_1_103412949_103413917"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#11
gene_A <- "GCA_014905855_1_human_AP023467_1_140645401_140646351"
gene_B <- "GCA_019393635_1_monito_del_monte_del_monte_CM033369_1_103611103_103612011"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#12
gene_A <- "GCA_014905855_1_human_AP023472_1_10812703_10813659"
gene_B <- "GCA_023159225_1_Pacific_pocket_mouse_pocket_mouse_CM041290_1_5156605_5157546"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#13
gene_A <- "GCA_011100635_1_common_brushtail_brushtail_CM021908_1_376121526_376122458"
gene_B <- "GCA_019393635_1_monito_del_monte_del_monte_CM033367_1_162704840_162705787"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#14
gene_A <- "GCA_014905855_1_human_AP023472_1_10817140_10818069"
gene_B <- "GCA_028571685_1_yellow_spotted_hyrax_hyrax_CM053142_1_73207293_73208222"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#15
gene_A <- "GCA_014905855_1_human_AP023472_1_10820225_10821163"
gene_B <- "GCA_028533215_1_North_Sulawesi_babirusa_Sulawesi_babirusa_CM052690_1_40321382_40322317"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#16
gene_A <- "GCA_011100635_1_common_brushtail_brushtail_CM021911_1_3093162_3094073"
gene_B <- "GCA_027887165_1_gray_short_tailed_opossum_short_tailed_opossum_CM051224_1_536859416_536860423"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#17
gene_A <- "GCA_000003025_6_pig_CM000816_5_61227816_61228733"
gene_B <- "GCA_023851605_1_southern_tamandua_tamandua_CM043333_1_73534496_73535434"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#18
#This is very large
gene_A <- "GCA_014905855_1_human_AP023472_1_10919443_10920354"
gene_B <- "GCA_028646465_1_greater_Indian_rhinoceros_Indian_rhinoceros_CM053529_1_11422840_11423802"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#19
gene_A <- "GCA_014905855_1_human_AP023472_1_11032450_11033406"
gene_B <- "GCA_028571685_1_yellow_spotted_hyrax_hyrax_CM053142_1_73335384_73336319"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#20
gene_A <- "GCA_009806435_2_rabbit_VIYN02003076_1_404441_405394"
gene_B <- "GCA_014633375_1_American_pika_pika_CM025752_1_19889374_19890330"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#21
gene_A <- "GCA_000165445_3_gray_mouse_lemur_mouse_lemur_CM007667_1_99851021_99851950"
gene_B <- "GCA_015220235_1_southern_two_toed_sloth_two_toed_sloth_CM026697_1_122115653_122116603"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#22
gene_A <- "GCA_014905855_1_human_AP023467_1_140671313_140672473"
gene_B <- "GCA_019393635_1_monito_del_monte_del_monte_CM033369_1_103546344_103547231"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#problem GCA_030035585_1_dugong_CM057468_1_128773908_128774828

###BIRDS

taxon <- "bird"

#1
gene_A <- "GCA_000146605_4_turkey_CM000963_2_106292103_106293122"
gene_B <- "GCA_011075105_1_birds_CM021733_1_113217225_113218178"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#2
gene_A <- "GCA_000247815_2_Collared_flycatcher_flycatcher_CM001991_1_110976344_110977297"
gene_B <- "GCA_013398505_2_white_breasted_antbird_antbird_CM031047_1_111141286_111142242"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#3
gene_A <- "GCA_000247815_2_Collared_flycatcher_flycatcher_CM001991_1_111000448_111001410"
gene_B <- "GCA_013398505_2_white_breasted_antbird_antbird_CM031047_1_111205638_111206600"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#4
gene_A <- "GCA_000247815_2_Collared_flycatcher_flycatcher_CM001991_1_110982449_110983402"
gene_B <- "GCA_009829145_1_lance_tailed_manakin_manakin_CM020535_1_111174821_111175768"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#5
gene_A <- "GCA_000738735_6_hooded_crow_crow_CM022407_2_109446664_109447587"
gene_B <- "GCA_009741485_1_Superb_fairywren_fairywren_CM019216_1_17339535_17340482"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#6
gene_A <- "GCA_000738735_6_hooded_crow_crow_CM022407_2_109464900_109465853"
gene_B <- "GCA_009741485_1_Superb_fairywren_fairywren_CM019216_1_17313818_17314777"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#7
gene_A <- "GCA_013398505_2_white_breasted_antbird_antbird_CM031047_1_111172660_111173607"
gene_B <- "GCA_009829145_1_lance_tailed_manakin_manakin_CM020535_1_111182041_111182982"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)



#8
gene_A <- "GCA_000738735_6_hooded_crow_crow_CM022407_2_109451665_109452618"
gene_B <- "GCA_013398505_2_white_breasted_antbird_antbird_CM031047_1_111186954_111187907"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#problem GCA_009819775_1_American_flamingo_flamingo_CM020217_1_121792898_121793851
#problem GCA_027574665_1_kagu_CM050650_1_119208351_119209304

#9
gene_A <- "GCA_009769605_1_Abyssinian_ground_hornbill_ground_hornbill_CM020024_1_51790451_51791404"
gene_B <- "GCA_020746105_1_Suruca_trogon_trogon_CM036619_1_98842435_98843388"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#10
gene_A <- "GCA_017639485_1_birds_CM030009_1_123269331_123270284"
gene_B <- "GCA_028858725_1_speckled_mousebird_mousebird_CM054346_1_41347296_41348249"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#11
gene_A <- "GCA_017639555_1_birds_CM030195_1_7015122_7016078"
gene_B <- "GCA_015220805_1_Red_fronted_tinkerbird_tinkerbird_CM026782_1_5682748_5683683"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#12
gene_A <- "GCA_028858705_1_Whooping_crane_crane_CM054307_1_6829726_6830679"
gene_B <- "GCA_020800305_1_South_Island_takahe_Island_takahe_CM036907_1_120101833_120102786"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#13
gene_A <- "GCA_949628215_1_European_shag_shag_OX451224_1_6936896_6937849"
gene_B <- "GCA_013368605_1_birds_CM023732_1_126739531_126740484"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#14
gene_A <- "GCA_003957555_2_Annas_hummingbird_hummingbird_CM012116_1_109234763_109235716"
gene_B <- "GCA_009769465_1_red_crested_turaco_turaco_CM019967_1_66225931_66226884"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#15
gene_A <- "GCA_004027225_2_Kakapo_CM013767_2_57151595_57152545"
gene_B <- "GCA_028858755_1_blue_and_yellow_macaw_macaw_CM054379_1_19963552_19964502"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#16
gene_A <- "GCA_026413225_1_great_bustard_bustard_CM048554_1_120147490_120148467"
gene_B <- "GCA_009769525_1_yellow_throated_sandgrouse_sandgrouse_CM020152_1_6081981_6082931"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)


##
#17
gene_A <- "GCA_028389875_1_lesser_rhea_rhea_CM051763_1_123214431_123215384"
gene_B <- "GCA_016128335_2_emu_CM027947_2_7426861_7427913"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#18
gene_A <- "GCA_000146605_4_turkey_CM000962_2_185213563_185214498"
gene_B <- "GCA_009769605_1_Abyssinian_ground_hornbill_ground_hornbill_CM020024_1_98380618_98381565"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#19
gene_A <- "GCA_000247815_2_Collared_flycatcher_flycatcher_CM001988_1_91783798_91784718"
gene_B <- "GCA_016128335_2_emu_CM027945_2_130798482_130799417"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

##REPTILE
#1
taxon <- "reptile"
gene_A <- "GCA_009819535_1_lizards_and_snakes_and_snakes_CM020431_1_24564949_24565932"
gene_B <- "GCA_027244095_1_graceful_crag_lizard_crag_lizard_CM050091_1_155932168_155933142"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#2
taxon <- "reptile"
gene_A <- "GCA_000090745_2_green_anole_anole_CM000938_1_147312068_147312976"
gene_B <- "GCA_011800845_1_common_lizard_lizard_CM022340_1_67048130_67049047"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#3
gene_A <- "GCA_000090745_2_green_anole_anole_CM000938_1_149660566_149661498"
gene_B <- "GCA_020142125_1_Desert_horned_lizard_horned_lizard_CM034703_1_259060556_259061485"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#problem GCA_000090745_2_green_anole_anole_CM000938_1_149645374_149646306
#problem GCA_000090745_2_green_anole_anole_CM000938_1_149712533_149713438
#4
gene_A <- "GCA_000090745_2_green_anole_anole_CM000938_1_149716271_149717173"
gene_B <- "GCA_020142125_1_Desert_horned_lizard_horned_lizard_CM034703_1_259130515_259131417"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#5
gene_A <- "GCA_016801065_1_plateau_fence_lizard_fence_lizard_CM028773_1_57618877_57619779"
gene_B <- "GCA_019175285_1_fence_lizard_lizard_CM032820_1_76042732_76043634"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#6
gene_A <- "GCA_000090745_2_green_anole_anole_CM000938_1_152394250_152395140"
gene_B <- "GCA_019175285_1_fence_lizard_lizard_CM032820_1_76095729_76096625"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#7
gene_A <- "GCA_000090745_2_green_anole_anole_CM000938_1_149683227_149684126"
gene_B <- "GCA_019175285_1_fence_lizard_lizard_CM032820_1_76081609_76082448"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#8
gene_A <- "GCA_016801065_1_plateau_fence_lizard_fence_lizard_CM028773_1_57658235_57659149"
gene_B <- "GCA_020142125_1_Desert_horned_lizard_horned_lizard_CM034703_1_259073946_259074872"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#9
gene_A <- "GCA_004329235_1_Common_wall_lizard_wall_lizard_CM014759_1_9195902_9196918"
gene_B <- "GCA_029931775_1_tarantolino_CM057264_1_14101260_14102378"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#10
gene_A <- "GCA_000241765_5_western_painted_turtle_painted_turtle_CM002657_2_92571340_92572296"
gene_B <- "GCA_029931775_1_tarantolino_CM057270_1_11043852_11044871"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#11
gene_A <- "GCA_004329235_1_Common_wall_lizard_wall_lizard_CM014745_1_1180936_1182117"
gene_B <- "GCA_020142125_1_Desert_horned_lizard_horned_lizard_CM034716_1_8413506_8414570" 
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#12
gene_A <- "GCA_027244095_1_graceful_crag_lizard_crag_lizard_CM050090_1_244740526_244741491"
gene_B <- "GCA_029931775_1_tarantolino_CM057270_1_11102660_11103379" 
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#13
gene_A <- "GCA_007399415_1_Goodes_thornscrub_tortoise_thornscrub_tortoise_CM017296_1_14263159_14264088"
gene_B <- "GCA_019425775_1_Swinhoes_soft_shelled_turtle_soft_shelled_turtle_CM033437_1_23687376_23688296"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#14
gene_A <- "GCA_007399415_1_Goodes_thornscrub_tortoise_thornscrub_tortoise_CM017296_1_217480381_217481319"
gene_B <- "GCA_026122505_1_Aldabra_giant_tortoise_giant_tortoise_CM047479_1_155550766_155551704" 
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#15
gene_A <- "GCA_000090745_2_green_anole_anole_CM000939_1_202851747_202852676"
gene_B <- "GCA_028640845_1_snakes_CM054098_1_8082001_8082927"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#17
gene_A <- "GCA_000090745_2_green_anole_anole_CM000939_1_202846349_202847317"
gene_B <- "GCA_029931775_1_tarantolino_CM057272_1_7459677_7460630"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#18
gene_A <- "GCA_004329235_1_Common_wall_lizard_wall_lizard_CM014759_1_4210965_4211912"
gene_B <- "GCA_027172205_1_Aeolian_wall_lizard_wall_lizard_CM049766_1_4134139_4135050" 
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#19
gene_A <- "GCA_004329235_1_Common_wall_lizard_wall_lizard_CM014759_1_4285121_4286032"
gene_B <- "GCA_009819535_1_lizards_and_snakes_and_snakes_CM020431_1_24614472_24615380"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#20
gene_A <- "GCA_004329235_1_Common_wall_lizard_wall_lizard_CM014759_1_4316967_4317908"
gene_B <- "GCA_011800845_1_common_lizard_lizard_CM022354_1_24946898_24947809"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#21
gene_A <- "GCA_004329235_1_Common_wall_lizard_wall_lizard_CM014759_1_4219979_4220890"
gene_B <- "GCA_009819535_1_lizards_and_snakes_and_snakes_CM020431_1_24622595_24623518"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_004329235_1_Common_wall_lizard_wall_lizard_CM014759_1_4306044_4306967"
gene_B <- "GCA_009819535_1_lizards_and_snakes_and_snakes_CM020431_1_24581228_24582151"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#problem GCA_011800845_1_common_lizard_lizard_CM022354_1_24957791_24958723
#problem GCA_011800845_1_common_lizard_lizard_CM022354_1_25011268_25012197


#single species --potentential problem
gene_A <- "GCA_030035675_1_Florida_worm_lizard_worm_lizard_CM057605_1_151746577_151747497"
gene_B <- "GCA_030035675_1_Florida_worm_lizard_worm_lizard_CM057605_1_151770839_151771756"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#single species --potentential problem
gene_A <- "GCA_029931775_1_tarantolino_CM057264_1_12044558_12045511"
gene_B <- "GCA_029931775_1_tarantolino_CM057264_1_12318714_12319658"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)


gene_A <- "GCA_000241765_5_western_painted_turtle_painted_turtle_ML621257_1_2728447_2729379"
gene_B <- "GCA_019425775_1_Swinhoes_soft_shelled_turtle_soft_shelled_turtle_CM033437_1_212300196_212301128"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_000090745_2_green_anole_anole_AAWZ02040413_1_3957_4862"
gene_B <- "GCA_020142125_1_Desert_horned_lizard_horned_lizard_CM034708_1_6618981_6619913"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_000090745_2_green_anole_anole_GL343381_1_415858_416832"
gene_B <- "GCA_019175285_1_fence_lizard_lizard_JAGXEY010004659_1_886_1794"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_000090745_2_green_anole_anole_GL343381_1_253683_254636"
gene_B <- "GCA_019175285_1_fence_lizard_lizard_CM032820_1_223432835_223433746"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_000090745_2_green_anole_anole_GL343381_1_272859_273767"
gene_B <- "GCA_020142125_1_Desert_horned_lizard_horned_lizard_CM034703_1_130957356_130958285"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_016801065_1_plateau_fence_lizard_fence_lizard_CM028773_1_198077124_198078044"
gene_B <- "GCA_019175285_1_fence_lizard_lizard_CM032820_1_223613560_223614600"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_000090745_2_green_anole_anole_GL343381_1_198025_198939"
gene_B <- "GCA_020142125_1_Desert_horned_lizard_horned_lizard_CM034703_1_130703626_130704546"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_004329235_1_Common_wall_lizard_wall_lizard_CM014759_1_4382219_4383145"
gene_B <- "GCA_030035675_1_Florida_worm_lizard_worm_lizard_CM057605_1_151638303_151639238"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#single species -- potential problem
gene_A <- "GCA_027244095_1_graceful_crag_lizard_crag_lizard_CM050091_1_155943173_155944099"
gene_B <- "GCA_027244095_1_graceful_crag_lizard_crag_lizard_CM050091_1_155966503_155967429"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_027244095_1_graceful_crag_lizard_crag_lizard_CM050091_1_155985271_155986206"
gene_B <- "GCA_029931775_1_tarantolino_CM057264_1_12520574_12521515"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_028583425_1_Leopard_gecko_gecko_CM053294_1_1478169_1479113"
gene_B <- "GCA_029931775_1_tarantolino_CM057264_1_12517856_12518794"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_000241765_5_western_painted_turtle_painted_turtle_CM002657_2_63297351_63298277"
gene_B <- "GCA_028017835_1_European_pond_turtle_pond_turtle_CM051432_1_142017176_142018105"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_016161935_1_Reevess_turtle_turtle_CM027978_1_73708883_73709806"
gene_B <- "GCA_007399415_1_Goodes_thornscrub_tortoise_thornscrub_tortoise_CM017298_1_150473356_150474210"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_015237465_2_Green_sea_turtle_sea_turtle_CM026900_1_134691895_134692833"
gene_B <- "GCA_030012505_1_hawksbill_sea_turtle_sea_turtle_CM057285_1_139076868_139077779"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_000241765_5_western_painted_turtle_painted_turtle_CM002657_2_63339507_63340427"
gene_B <- "GCA_027887155_1_diamondback_terrapin_terrapin_CM051160_1_72033659_72034579"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#problem GCA_019425775_1_Swinhoes_soft_shelled_turtle_soft_shelled_turtle_CM033440_1_132567821_132568726
#problem GCA_019425775_1_Swinhoes_soft_shelled_turtle_soft_shelled_turtle_CM033440_1_132545728_132546651

gene_A <- "GCA_007399415_1_Goodes_thornscrub_tortoise_thornscrub_tortoise_CM017298_1_151188506_151189441"
gene_B <- "GCA_019425775_1_Swinhoes_soft_shelled_turtle_soft_shelled_turtle_CM033440_1_133019563_133020498"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)


gene_A <- "GCA_000241765_5_western_painted_turtle_painted_turtle_CM002657_2_113843147_113844076"
gene_B <- "GCA_019425775_1_Swinhoes_soft_shelled_turtle_soft_shelled_turtle_CM033440_1_36352658_36353602"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#problem "GCA_028583425_1_Leopard_gecko_gecko_CM053277_1_56434222_56435160"


##FISH

#1
taxon <- "fish"
gene_A <- "GCA_020184715_1_zebrafish_CM035062_1_8785507_8786481"
gene_B <- "GCA_022985175_1_rohu_CM040943_1_35928581_35929564"#"GCA_902713425_2_sterlet_OV754645_1_3742484_3743443"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#2
gene_A <- "GCA_008369825_1_bony_fishes_fishes_CM017904_1_20348123_20349088"
gene_B <- "GCA_027580225_1_oriental_weatherfish_weatherfish_CM050767_1_4876998_4877969"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#problem not in tree GCA_015846415_1_bony_fishes_fishes_CM027556_1_2252996_2253967
#3
gene_A <- "GCA_029783895_1_bony_fishes_fishes_CM056468_1_20252008_20252979"
gene_B <- "GCA_024868665_1_bony_fishes_fishes_CM045789_1_21180072_21181043"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#problem GCA_019155185_1_bony_fishes_fishes_CM032806_1_1453571_1454542
#problem GCA_028031935_1_bony_fishes_fishes_CM051498_1_27977238_27978203
#4
gene_A <- "GCA_019703515_2_Chinese_sucker_sucker_CM033865_2_5525094_5526068"
gene_B <- "GCA_025860055_1_razorback_sucker_sucker_CM047335_1_35555767_35556741"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#5
gene_A <- "GCA_001660625_3_channel_catfish_catfish_CM004418_2_14099251_14100225"
gene_B <- "GCA_023375975_1_Mexican_tetra_tetra_CM041920_1_13227842_13228816"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#problem GCA_902362185_1_milkfish_LR697106_1_53198809_53199810
#6
gene_A <- "GCA_000238955_5_zebra_mbuna_mbuna_CM009191_2_9273080_9274018"
gene_B <- "GCA_023634145_1_grayling_CM042401_1_20770197_20771174"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#7
gene_A <- "GCA_013347855_1_European_eel_eel_CM023574_1_12805285_12806259"
gene_B <- "GCA_019176425_1_tarpon_CM032878_1_10029719_10030699"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#8
gene_A <- "GCA_023856365_1_bony_fishes_fishes_CM043636_1_9878681_9879655"
gene_B <- "GCA_900964775_1_Asian_bonytongue_bonytongue_LR584082_1_26428452_26429426"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#9
gene_A <- "GCA_000242695_1_spotted_gar_gar_CM001404_1_13234574_13235548"
gene_B <- "GCA_017591415_1_bowfin_CM030127_1_32825655_32826629"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#10
gene_A <- "GCA_017654505_1_Mississippi_paddlefish_paddlefish_CM030280_1_2918076_2919035"
gene_B <- "GCA_902713425_2_sterlet_OV754645_1_3742484_3743443"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#11
gene_A <- "GCA_000242695_1_spotted_gar_gar_CM001407_1_36789853_36790839"
gene_B <- "GCA_902713425_2_sterlet_OV754670_1_104812105_104813088"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#12
gene_A <- "GCA_020184715_1_zebrafish_CM035062_1_33754007_33755002"
gene_B <- "GCA_900700375_2_denticle_herring_herring_LR535823_1_10834401_10835375"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#13
gene_A <- "GCA_000972845_2_large_yellow_croaker_yellow_croaker_CM011638_1_18923207_18924154"
gene_B <- "GCA_023856365_1_bony_fishes_fishes_CM043636_1_4522064_4523020"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#14
gene_A <- "GCA_013347855_1_European_eel_eel_CM023574_1_16060331_16061260"
gene_B <- "GCA_022829145_1_bony_fishes_fishes_CM040744_1_28660001_28661041"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#problem
#15
gene_A <- "GCA_902713425_2_sterlet_OV754629_1_6972647_6973609"
gene_B <- "GCA_017654505_1_Mississippi_paddlefish_paddlefish_CM030271_1_24327344_24328294"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#problem

#16 -- flag
gene_A <- "GCA_017654505_1_Mississippi_paddlefish_paddlefish_CM030261_1_11279909_11280838"
gene_B <- "GCA_902713425_2_sterlet_OV754630_1_2968848_2969795"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#17
gene_A <- "GCA_020184715_1_zebrafish_CM035063_1_40277028_40277969"
gene_B <- "GCA_027580225_1_oriental_weatherfish_weatherfish_CM050764_1_32171144_32172175"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#18
gene_A <- "GCA_020184715_1_zebrafish_CM035063_1_24813294_24814211"
gene_B <- "GCA_021613375_1_bony_fishes_fishes_CM038709_1_37191286_37192248"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#19
gene_A <- "GCA_001660625_3_channel_catfish_catfish_CM004419_2_12510333_12511271"
gene_B <- "GCA_023375975_1_Mexican_tetra_tetra_CM041919_1_16927882_16928847"#"GCA_019176425_1_tarpon_CM032886_1_35582396_35583358"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#20
gene_A <- "GCA_015220715_1_red_bellied_piranha_piranha_CM026726_1_39354591_39355517"
gene_B <- "GCA_030014385_1_bony_fishes_fishes_CM057363_1_27932151_27933092"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#21
gene_A <- "GCA_006149115_2_sockeye_salmon_salmon_CM016829_1_39403030_39403977"
gene_B <- "GCA_021917145_1_delta_smelt_smelt_CM038994_1_10863922_10864872"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#23
gene_A <- "GCA_002021735_2_coho_salmon_salmon_CM007743_2_41867083_41868030"
gene_B <- "GCA_023658055_1_eulachon_CM042867_1_10751482_10752429"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#24
gene_A <- "GCA_017589495_2_allis_shad_shad_CM030048_1_34659688_34660620"
gene_B <- "GCA_021917145_1_delta_smelt_smelt_CM038994_1_10858100_10859026"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#25
gene_A <- "GCA_018136845_1_African_bonytongue_bonytongue_CM030897_1_3941228_3942121"
gene_B <- "GCA_900964775_1_Asian_bonytongue_bonytongue_LR584077_1_27619487_27620413"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#26
gene_A <- "GCA_018136845_1_African_bonytongue_bonytongue_CM030897_1_3944183_3945091"
gene_B <- "GCA_900964775_1_Asian_bonytongue_bonytongue_LR584077_1_27622081_27622986"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#27
gene_A <- "GCA_018136845_1_African_bonytongue_bonytongue_CM030897_1_3955201_3956103"
gene_B <- "GCA_900964775_1_Asian_bonytongue_bonytongue_LR584077_1_27695870_27696781"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#28
#single species -- potential problem
gene_A <- "GCA_900964775_1_Asian_bonytongue_bonytongue_LR584077_1_27646939_27647859"
gene_B <- "GCA_900964775_1_Asian_bonytongue_bonytongue_LR584077_1_27691761_27692687"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#problem GCA_900964775_1_Asian_bonytongue_bonytongue_LR584077_1_27704834_27705739
#problem GCA_018136845_1_African_bonytongue_bonytongue_CM030897_1_3947349_3948401
#29
gene_A <- "GCA_013347855_1_European_eel_eel_CM023579_1_20789561_20790514"
gene_B <- "GCA_019176425_1_tarpon_CM032886_1_35582396_35583358"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#30
gene_A <- "GCA_000633615_2_guppy_CM002707_1_35291117_35292031"
gene_B <- "GCA_027744805_2_ringed_pipefish_pipefish_CM050986_1_10397183_10398082"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#31
gene_A <- "GCA_000972845_2_large_yellow_croaker_yellow_croaker_CM011632_1_9817665_9818654"
gene_B <- "GCA_949987615_1_pollack_OX465126_1_40899026_40899976"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#32

gene_A <- "GCA_949987555_1_large_eye_snaggletooth_snaggletooth_OX465220_1_21178075_21179022"
gene_B <- "GCA_017639675_1_bony_fishes_fishes_CM030177_1_3244837_3245784"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#33
gene_A <- "GCA_002021735_2_coho_salmon_salmon_CM007719_2_54677541_54678479"
gene_B <- "GCA_000721915_3_northern_pike_pike_CM002841_3_19747318_19748247"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#34
gene_A <- "GCA_002021735_2_coho_salmon_salmon_ML762283_1_8340692_8341642"
gene_B <- "GCA_000721915_3_northern_pike_pike_CM002836_3_33186684_33187664"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#35
gene_A <- "GCA_002021735_2_coho_salmon_salmon_CM007731_2_55256996_55257940"
gene_B <- "GCA_905237065_2_Atlantic_salmon_salmon_HG993286_1_9649700_9650629"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#36
gene_A <- "GCA_029465975_1_Alaska_blackfish_blackfish_CM055738_1_2778864_2779820"
gene_B <- "GCA_000721915_3_northern_pike_pike_CM002845_3_37354061_37354993"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#37
gene_A <- "GCA_017589495_2_allis_shad_shad_CM030048_1_15593197_15594123"
gene_B <- "GCA_902362185_1_milkfish_LR697115_1_35657502_35658626"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#38
gene_A <- "GCA_900700375_2_denticle_herring_herring_LR535832_1_16318375_16319394"
gene_B <- "GCA_902362185_1_milkfish_LR697115_1_4705556_4706554"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#39
gene_A <- "GCA_001660625_3_channel_catfish_catfish_CM004419_2_28398276_28399211"
gene_B <- "GCA_014805685_1_bony_fishes_fishes_CM025887_1_22374220_22375233"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#40
gene_A <- "GCA_014805685_1_bony_fishes_fishes_CM025887_1_22360145_22361068"
gene_B <- "GCA_030014385_1_bony_fishes_fishes_CM057363_1_3110999_3111934"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#41
gene_A <- "GCA_001660625_3_channel_catfish_catfish_CM004419_2_28423997_28424989"
gene_B <- "GCA_030014155_1_bony_fishes_fishes_CM057318_1_19588512_19589474"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#42
gene_A <- "GCA_019097595_1_bony_fishes_fishes_CM032599_1_12718323_12719255"
gene_B <- "GCA_030014155_1_bony_fishes_fishes_CM057318_1_19592170_19593105"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#43
gene_A <- "GCA_015220715_1_red_bellied_piranha_piranha_CM026726_1_17225209_17226207"
gene_B <- "GCA_017165825_1_African_pike_characin_pike_characin_CM029369_1_9767376_9768365"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)
#44
gene_A <- "GCA_015220715_1_red_bellied_piranha_piranha_CM026726_1_17149306_17150307"
gene_B <- "GCA_029633875_1_trahira_CM056254_1_45013913_45014902"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

#problem GCA_023375975_1_Mexican_tetra_tetra_CM041919_1_37805588_37806589

gene_A <- "GCA_001660625_3_channel_catfish_catfish_CM004419_2_28433098_28434048"
gene_B <- "GCA_027579695_1_lesser_salmon_catfish_salmon_catfish_CM050628_1_72467173_72468126"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_015220715_1_red_bellied_piranha_piranha_CM026726_1_1166414_1167262"
gene_B <- "GCA_023375975_1_Mexican_tetra_tetra_CM041919_1_37754441_37755415"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_015220715_1_red_bellied_piranha_piranha_CM026726_1_17191947_17192927"
gene_B <- "GCA_021613375_1_bony_fishes_fishes_CM038709_1_16834703_16835695"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_017165825_1_African_pike_characin_pike_characin_CM029369_1_9741562_9742557"
gene_B <- "GCA_023375975_1_Mexican_tetra_tetra_CM041919_1_37790333_37791298"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_015220715_1_red_bellied_piranha_piranha_CM026726_1_17200766_17201683"
gene_B <- "GCA_021613375_1_bony_fishes_fishes_CM038709_1_16862800_16863741"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

gene_A <- "GCA_017654505_1_Mississippi_paddlefish_paddlefish_CM030286_1_2575642_2576640"
gene_B <- "GCA_902713425_2_sterlet_OV754649_1_3633004_3634005"
result_message <- analyze_genes(gene_A, gene_B, taxon)
print(result_message)

# Find the subset of ordered_tips_amphibian that is not in df_ordered_tips4_amphibian$gene
amphibian_remaining <- subset(ordered_tips_amphibian, !(ordered_tips_amphibian %in% df_ordered_tips4_amphibian$gene))
bird_remaining <- subset(ordered_tips_bird, !(ordered_tips_bird %in% df_ordered_tips4_bird$gene))
mammal_remaining <- subset(ordered_tips_mammal, !(ordered_tips_mammal %in% df_ordered_tips4_mammal$gene))
reptile_remaining <- subset(ordered_tips_reptile, !(ordered_tips_reptile %in% df_ordered_tips4_reptile$gene))
fish_remaining <- subset(ordered_tips_fish, !(ordered_tips_fish %in% df_ordered_tips4_fish$gene))

# Function to calculate summary for a specific taxon dataframe
calculate_summary <- function(df_ordered_tips4_taxon) {
  summary_table <- df_ordered_tips4_taxon %>%
    group_by(gene_family) %>%
    summarise(
      num_scg_1 = sum(scg == 1),
      num_scg_0 = sum(scg == 0),
      num_scg_2 = sum(scg == 2),
      num_scg_3 = sum(scg == 3),
      num_scg_10 = sum(scg ==10),
      num_scg_11 = sum(scg ==11)
    ) %>%
    summarise(
      total_gene_families = n(),
      total_scg_1 = sum(num_scg_1),
      total_scg_0 = sum(num_scg_0),
      total_scg_2 = sum(num_scg_2),
      total_scg_3 = sum(num_scg_3),
      total_scg_10 = sum(num_scg_10),
      total_scg_11 = sum(num_scg_11),
      percent_scg_1 = total_scg_1 / (total_scg_1 + total_scg_0 + total_scg_2 + total_scg_3 + total_scg_10+ total_scg_11),
      percent_scg_2 = total_scg_2 / (total_scg_1 + total_scg_0 + total_scg_2 + total_scg_3 + total_scg_10 + total_scg_11),
      percent_scg_3= total_scg_3 / (total_scg_1 + total_scg_0 + total_scg_2 + total_scg_3 + total_scg_10 + total_scg_11),
      percent_scg_10= total_scg_10 / (total_scg_1 + total_scg_0 + total_scg_2 + total_scg_3 + total_scg_10 + total_scg_11),
      percent_scg_11= total_scg_11 / (total_scg_1 + total_scg_0 + total_scg_2 + total_scg_3 + total_scg_10 + total_scg_11)
      # percent_scg_1 = total_scg_1 / (total_scg_1 + total_scg_0 + total_scg_2 + num_scg_3 + num_scg_10+ num_scg_11),
      # percent_scg_2 = total_scg_2 / (total_scg_1 + total_scg_0 + total_scg_2 + num_scg_3 + num_scg_10 + num_scg_11),
      # percent_scg_3= total_scg_3 / (total_scg_1 + total_scg_0 + total_scg_2 + num_scg_3 + num_scg_10 + num_scg_11),
      # percent_scg_10= total_scg_10 / (total_scg_1 + total_scg_0 + total_scg_2 + num_scg_3 + num_scg_10 + num_scg_11),
      # percent_scg_11= total_scg_11 / (total_scg_1 + total_scg_0 + total_scg_2 + num_scg_3 + num_scg_10 + num_scg_11)
    )
  
  return(summary_table)
}
# Calculate summary for amphibian dataframe
summary_amphibian <- calculate_summary(df_ordered_tips4_amphibian)
summary_mammal <- calculate_summary(df_ordered_tips4_mammal)
summary_bird <- calculate_summary(df_ordered_tips4_bird)
summary_reptile <- calculate_summary(df_ordered_tips4_reptile)
summary_fish <- calculate_summary(df_ordered_tips4_fish)

# Printing summaries
print("Summary for Amphibian")
print(summary_amphibian)

print("Summary for Mammal")
print(summary_mammal)

print("Summary for Bird")
print(summary_bird)

# Create a vector of average scg values
average_scg <- c(mean(summary_amphibian$percent_scg_1), mean(summary_bird$percent_scg_1), mean(summary_mammal$percent_scg_1),mean(summary_reptile$percent_scg_1),mean(summary_fish$percent_scg_1) )
average_scg2 <- c(mean(summary_amphibian$percent_scg_2), mean(summary_bird$percent_scg_2), mean(summary_mammal$percent_scg_2),mean(summary_reptile$percent_scg_2),mean(summary_fish$percent_scg_2) )
average_scg3 <- c(mean(summary_amphibian$percent_scg_3), mean(summary_bird$percent_scg_3), mean(summary_mammal$percent_scg_3),mean(summary_reptile$percent_scg_3),mean(summary_fish$percent_scg_3) )
average_scg10 <- c(mean(summary_amphibian$percent_scg_10), mean(summary_bird$percent_scg_10), mean(summary_mammal$percent_scg_10),mean(summary_reptile$percent_scg_10),mean(summary_fish$percent_scg_10) )
average_scg11 <- c(mean(summary_amphibian$percent_scg_11), mean(summary_bird$percent_scg_11), mean(summary_mammal$percent_scg_11),mean(summary_reptile$percent_scg_11),mean(summary_fish$percent_scg_11) )

average_scg_total <- c((mean(summary_amphibian$percent_scg_1) + mean(summary_amphibian$percent_scg_2) + mean(summary_amphibian$percent_scg_3)), 
                        mean(summary_bird$percent_scg_1+ mean(summary_bird$percent_scg_2) + mean(summary_bird$percent_scg_3)), 
                        mean(summary_mammal$percent_scg_1) + mean(summary_mammal$percent_scg_2) + mean(summary_mammal$percent_scg_3),
                        mean(summary_reptile$percent_scg_1) + mean(summary_reptile$percent_scg_2) + mean(summary_reptile$percent_scg_3), 
                        mean(summary_fish$percent_scg_1) + mean(summary_fish$percent_scg_2) + mean(summary_fish$percent_scg_3) )

# Create a vector of taxa names
taxa <- c("Amphibian", "Bird", "Mammal", "Reptile", "Fish")

# # Create a bar plot
# pdf("Fraction_of_all_genes_that_are_single_copy_genes_min3_50percent.pdf")
# barplot(average_scg, names.arg = taxa, col = "skyblue",
#         main = "Fraction of all genes that are single copy",
#         xlab = "Taxa", ylab = "Average SCG Value")
# dev.off()

# Combine average_scg and average_scg2 into a matrix
values_matrix <- rbind(average_scg, average_scg2, average_scg3)

# Create the stacked bar plot
pdf("Stacked_barplot_average_scg_and_average_scg2.pdf")
barplot(values_matrix, col = c("#166A53", "#AC71F3", "#063781"),
        main = "Fraction of all genes that are single copy",
        xlab = "Taxa", ylab = "Proportion single copy",
        legend.text = c("1 copy in more than 90%","75 to 90%", "50 to 75%"),
        args.legend = list(x = "topleft"),
  names.arg = taxa)
# Add text labels for each bar

pdf("Stacked_barplot_average_scg_and_average_total.pdf")
barplot(average_scg_total, col = c("#166A53"),
        main = "Fraction of all genes that are copy number constrained",
        xlab = "Taxa", ylab = "Proportion copy number constrained",
        names.arg = taxa)
# Add text labels for each bar

dev.off()

values_matrix <- rbind(average_scg10, average_scg11)

# Create the stacked bar plot
pdf("Stacked_barplot_average_scg_and_average_scg10.pdf")
barplot(values_matrix, col = c("#166A53", "#AC71F3"),
        main = "Fraction of all genes that are multiple copy",
        xlab = "Taxa", ylab = "Proportion multiple copy",
        legend.text = c("2+ copies in more than 50%", "2+ copies in more than 30%"),
        args.legend = list(x = "topleft"),
        names.arg = taxa)
# Add text labels for each bar

dev.off()


# Function to get scg value for unique gene_family
get_scg_for_unique_gene_family <- function(df_ordered_tips4_taxon) {
  unique_gene_families <- unique(df_ordered_tips4_taxon$gene_family)
  
  gene_family_scg <- sapply(unique_gene_families, function(gf) {
    row <- df_ordered_tips4_taxon %>%
      filter(gene_family == gf) %>%
      slice(1)
    
    if (nrow(row) > 0) {
      return(row$scg)
    } else {
      return(NA)
    }
  })
  
  result <- data.frame(gene_family = unique_gene_families, scg = gene_family_scg)
  return(result)
}

# Get scg values for unique gene_family in df_ordered_tips4_amphibian
scg_for_unique_gene_family_amphibian <- get_scg_for_unique_gene_family(df_ordered_tips4_amphibian)
scg_for_unique_gene_family_mammal <- get_scg_for_unique_gene_family(df_ordered_tips4_mammal)
scg_for_unique_gene_family_bird <- get_scg_for_unique_gene_family(df_ordered_tips4_bird)
scg_for_unique_gene_family_reptile <- get_scg_for_unique_gene_family(df_ordered_tips4_reptile)
scg_for_unique_gene_family_fish <- get_scg_for_unique_gene_family(df_ordered_tips4_fish)

# Calculate average scg for unique gene_family in df_ordered_tips4_amphibian
average_scg_amphibian <- mean(scg_for_unique_gene_family_amphibian$scg, na.rm = TRUE)
proportion_scg_1_amphibian <- sum(scg_for_unique_gene_family_amphibian$scg == 1, na.rm = TRUE) / length(scg_for_unique_gene_family_amphibian$scg)
proportion_scg_2_amphibian <- sum(scg_for_unique_gene_family_amphibian$scg == 2, na.rm = TRUE) / length(scg_for_unique_gene_family_amphibian$scg)
proportion_scg_3_amphibian <- sum(scg_for_unique_gene_family_amphibian$scg == 3, na.rm = TRUE) / length(scg_for_unique_gene_family_amphibian$scg)
proportion_scg_10_amphibian <- sum(scg_for_unique_gene_family_amphibian$scg == 10, na.rm = TRUE) / length(scg_for_unique_gene_family_amphibian$scg)
proportion_scg_11_amphibian <- sum(scg_for_unique_gene_family_amphibian$scg == 11, na.rm = TRUE) / length(scg_for_unique_gene_family_amphibian$scg)
print(paste("Average scg for unique gene_family in df_ordered_tips4_amphibian:", average_scg_amphibian))

# Calculate average scg for unique gene_family in df_ordered_tips4_bird
average_scg_bird <- mean(scg_for_unique_gene_family_bird$scg, na.rm = TRUE)
proportion_scg_1_bird <- sum(scg_for_unique_gene_family_bird$scg == 1, na.rm = TRUE) / length(scg_for_unique_gene_family_bird$scg)
proportion_scg_2_bird <- sum(scg_for_unique_gene_family_bird$scg == 2, na.rm = TRUE) / length(scg_for_unique_gene_family_bird$scg)
proportion_scg_3_bird <- sum(scg_for_unique_gene_family_bird$scg == 3, na.rm = TRUE) / length(scg_for_unique_gene_family_bird$scg)
proportion_scg_10_bird <- sum(scg_for_unique_gene_family_bird$scg == 10, na.rm = TRUE) / length(scg_for_unique_gene_family_bird$scg)
proportion_scg_11_bird <- sum(scg_for_unique_gene_family_bird$scg == 11, na.rm = TRUE) / length(scg_for_unique_gene_family_bird$scg)
print(paste("Average scg for unique gene_family in df_ordered_tips4_bird", average_scg_bird))

# Calculate average scg for unique gene_family in df_ordered_tips4_mammal
average_scg_mammal <- mean(scg_for_unique_gene_family_mammal$scg, na.rm = TRUE)
proportion_scg_1_mammal <- sum(scg_for_unique_gene_family_mammal$scg == 1, na.rm = TRUE) / length(scg_for_unique_gene_family_mammal$scg)
proportion_scg_2_mammal <- sum(scg_for_unique_gene_family_mammal$scg == 2, na.rm = TRUE) / length(scg_for_unique_gene_family_mammal$scg)
proportion_scg_3_mammal <- sum(scg_for_unique_gene_family_mammal$scg == 3, na.rm = TRUE) / length(scg_for_unique_gene_family_mammal$scg)
proportion_scg_10_mammal <- sum(scg_for_unique_gene_family_mammal$scg == 10, na.rm = TRUE) / length(scg_for_unique_gene_family_mammal$scg)
proportion_scg_11_mammal <- sum(scg_for_unique_gene_family_mammal$scg == 11, na.rm = TRUE) / length(scg_for_unique_gene_family_mammal$scg)
print(paste("Average scg for unique gene_family in df_ordered_tips4_mammal:", average_scg_mammal))

# Calculate average scg for unique gene_family in df_ordered_tips4_reptile
average_scg_reptile <- mean(scg_for_unique_gene_family_reptile$scg, na.rm = TRUE)
proportion_scg_1_reptile <- sum(scg_for_unique_gene_family_reptile$scg == 1, na.rm = TRUE) / length(scg_for_unique_gene_family_reptile$scg)
proportion_scg_2_reptile <- sum(scg_for_unique_gene_family_reptile$scg == 2, na.rm = TRUE) / length(scg_for_unique_gene_family_reptile$scg)
proportion_scg_3_reptile <- sum(scg_for_unique_gene_family_reptile$scg == 3, na.rm = TRUE) / length(scg_for_unique_gene_family_reptile$scg)
proportion_scg_10_reptile <- sum(scg_for_unique_gene_family_reptile$scg == 10, na.rm = TRUE) / length(scg_for_unique_gene_family_reptile$scg)
proportion_scg_11_reptile <- sum(scg_for_unique_gene_family_reptile$scg == 11, na.rm = TRUE) / length(scg_for_unique_gene_family_reptile$scg)
print(paste("Average scg for unique gene_family in df_ordered_tips4_reptile:", average_scg_reptile))

# Calculate average scg for unique gene_family in df_ordered_tips4_fish
average_scg_fish <- mean(scg_for_unique_gene_family_fish$scg, na.rm = TRUE)
proportion_scg_1_fish <- sum(scg_for_unique_gene_family_fish$scg == 1, na.rm = TRUE) / length(scg_for_unique_gene_family_fish$scg)
proportion_scg_2_fish <- sum(scg_for_unique_gene_family_fish$scg == 2, na.rm = TRUE) / length(scg_for_unique_gene_family_fish$scg)
proportion_scg_3_fish <- sum(scg_for_unique_gene_family_fish$scg == 3, na.rm = TRUE) / length(scg_for_unique_gene_family_fish$scg)
proportion_scg_10_fish <- sum(scg_for_unique_gene_family_fish$scg == 10, na.rm = TRUE) / length(scg_for_unique_gene_family_fish$scg)
proportion_scg_11_fish <- sum(scg_for_unique_gene_family_fish$scg == 11, na.rm = TRUE) / length(scg_for_unique_gene_family_fish$scg)
print(paste("Average scg for unique gene_family in df_ordered_tips4_fish:", average_scg_fish))


# Create a vector of average scg values
# average_scg <- c(average_scg_amphibian, average_scg_bird, average_scg_mammal)
proportion_scg_1 <- c(proportion_scg_1_amphibian, proportion_scg_1_bird, proportion_scg_1_mammal, proportion_scg_1_reptile, proportion_scg_1_fish )
proportion_scg_2 <- c(proportion_scg_2_amphibian, proportion_scg_2_bird, proportion_scg_2_mammal, proportion_scg_2_reptile, proportion_scg_2_fish )
proportion_scg_3 <- c(proportion_scg_3_amphibian, proportion_scg_3_bird, proportion_scg_3_mammal, proportion_scg_3_reptile, proportion_scg_3_fish )
proportion_scg_10 <- c(proportion_scg_10_amphibian, proportion_scg_10_bird, proportion_scg_10_mammal, proportion_scg_10_reptile, proportion_scg_10_fish )
proportion_scg_11 <- c(proportion_scg_11_amphibian, proportion_scg_11_bird, proportion_scg_11_mammal, proportion_scg_11_reptile, proportion_scg_11_fish )

# Create a vector of taxa names
taxa <- c("Amphibian", "Bird", "Mammal", "Reptile", "Fish")

# # Create a bar plot
# pdf("Fraction_of_genes_that_are_multi_copy_over30.pdf")
# barplot(proportion_scg_10, names.arg = taxa, col = "skyblue",
#         main = "Fraction of all genes that are multi-copy (>30%)",
#         xlab = "Taxa", ylab = "Average SCG Value")
# dev.off()

# Combine average_scg and average_scg2 into a matrix
values_matrix <- rbind(proportion_scg_10, proportion_scg_11)

# Create the stacked bar plot
pdf("Multi_copy_gene_family.pdf")
barplot(values_matrix, col = c("#166A53", "#AC71F3"),
        main = "Fraction of all genes families that are multi copy",
        xlab = "Taxa", ylab = "Proportion multi-copy",
        legend.text = c("2+ copies in more than 50%", "2+ copies in more than 30%"),
        args.legend = list(x = "topleft"),
        names.arg = taxa)

# Combine average_scg and average_scg2 into a matrix
# values_matrix <- rbind(proportion_scg_1, proportion_scg_2, proportion_scg_3)
# pdf("Multi_copy_gene_family.pdf")
# barplot(average_scg10, col = c("#166A53"),
#         main = "Fraction of all genes families that are multi copy",
#         xlab = "Taxa", ylab = "Proportion multi-copy",
#         args.legend = list(x = "topleft"),
#         names.arg = taxa)

# Add text labels for each bar

dev.off()

values_matrix <- rbind(proportion_scg_1, proportion_scg_2,proportion_scg_3 )
pdf("Single_copy_gene_family.pdf")
barplot(values_matrix, col = c("#166A53", "#AC71F3",  "#063781"),
        main = "Fraction of all genes families that are single copy",
        xlab = "Taxa", ylab = "Proportion single-copy",
        legend.text = c("1 copy in more than 90%","75 to 90%", "50 to 75%"),
        args.legend = list(x = "topleft"),
        names.arg = taxa)

dev.off()

  #Function to filter unique gene_family values based on scg value
# Function to filter unique gene_family values based on scg value
get_unique_scg_values <- function(df_ordered_tips4_taxon) {
  # Group by gene_family and filter rows with unique gene_family values
  unique_gene_families <- df_ordered_tips4_taxon %>%
    group_by(gene_family) %>%
    ungroup()
  
  
  # Filter rows with scg=0 and scg=1 respectively and combine them into a unified list
  unique_scg <- unique_gene_families %>%
    group_by(scg, gene_family) %>%
    summarise(across(everything(), first)) %>%
    ungroup() %>%
    arrange(scg)
  
  return(as.data.frame(unique_scg))
  #return (unique_gene_families)
}

# Apply the function to the dataframes for each taxon
unique_scg_amphibian <- get_unique_scg_values(df_ordered_tips4_amphibian)
unique_scg_bird <- get_unique_scg_values(df_ordered_tips4_bird)
unique_scg_mammal <- get_unique_scg_values(df_ordered_tips4_mammal)
# 
# mean_tmrca_scg0_amphibian <- mean(unique_scg_amphibian$tmrca[unique_scg_amphibian$scg == 0], na.rm = TRUE)
# mean_tmrca_scg0_amphibian
# mean_tmrca_scg1_amphibian <- mean(unique_scg_amphibian$tmrca[unique_scg_amphibian$scg == 1], na.rm = TRUE)
# mean_tmrca_scg1_amphibian
# mean_tmrca_scg0_bird <- mean(unique_scg_bird$tmrca[unique_scg_bird$scg == 0], na.rm = TRUE)
# mean_tmrca_scg0_bird
# mean_tmrca_scg1_bird <- mean(unique_scg_bird$tmrca[unique_scg_bird$scg == 1], na.rm = TRUE)
# mean_tmrca_scg1_bird
# mean_tmrca_scg0_mammal <- mean(unique_scg_mammal$tmrca[unique_scg_mammal$scg == 0], na.rm = TRUE)
# mean_tmrca_scg0_mammal
# mean_tmrca_scg1_mammal <- mean(unique_scg_mammal$tmrca[unique_scg_mammal$scg == 1], na.rm = TRUE)
# mean_tmrca_scg1_mammal

genes_accessions$Accession <- gsub("\\.", "_", genes_accessions$Accession)
# Combining data frames by column binding (cbind)
df_ordered_combo <- rbind(df_ordered_tips4_mammal, df_ordered_tips4_amphibian, df_ordered_tips4_bird, df_ordered_tips4_fish, df_ordered_tips4_reptile)

# Assuming 'gene' column exists in df_ordered_combo
df_ordered_combo$accession <- sub("^([^_]+_[^_]+_[^_]+).*", "\\1", df_ordered_combo$gene)

# Initialize an empty vector to store scg123 counts
scg123_counts <- numeric()

# Loop through each row in genes_accessions
for (i in 1:nrow(genes_accessions)) {
  # Extract Accession value for each row in genes_accessions
  current_accession <- genes_accessions$Accession[i]
  
  # Subset df_ordered_combo based on current Accession value
  subset_df <- subset(df_ordered_combo, accession == current_accession)
  
  # Count the number of rows where scg equals 1, 2, or 3 in the subset
  count_scg123 <- sum(subset_df$scg %in% 1:3)
  
  # Append the count to scg123_counts vector
  scg123_counts <- c(scg123_counts, count_scg123)
}

# Add scg123 column to genes_accessions with the counts
genes_accessions$scg123 <- scg123_counts

pdf("NumGenes_vs_NumSCO.pdf")
plot(genes_accessions$`Number of Genes`, genes_accessions$scg123,
     xlab = "Number of Genes", ylab = "Number of SCO",
     main = "Scatter Plot of Number of Genes vs Numer of SCOs")
dev.off()

# Calculate correlation coefficient and p-value
correlation_result <- cor.test(genes_accessions$`Number of Genes`, genes_accessions$scg123)

# Access the correlation coefficient and p-value
correlation_coefficient <- correlation_result$estimate
p_value <- correlation_result$p.value

# Print the correlation coefficient and p-value
print(paste("Correlation Coefficient:", correlation_coefficient))
print(paste("P-value:", p_value))

genes_accessions$prop_scg <- genes_accessions$scg123/genes_accessions$`Number of Genes`

# Calculate correlation coefficient and p-value
correlation_result <- cor.test(genes_accessions$`Number of Genes`, genes_accessions$prop_scg)

# Access the correlation coefficient and p-value
correlation_coefficient <- correlation_result$estimate
p_value <- correlation_result$p.value

# Print the correlation coefficient and p-value
print(paste("Correlation Coefficient:", correlation_coefficient))
print(paste("P-value:", p_value))

# Subset the data for gene_accessions$class == "Amphibia"
subset_data <- subset(genes_accessions, class == "Amphibia")

# Calculate correlation coefficient and p-value for the subset
correlation_result <- cor.test(subset_data$`Number of Genes`, subset_data$prop_scg)

# Access the correlation coefficient and p-value
correlation_coefficient <- correlation_result$estimate
p_value <- correlation_result$p.value

# Print the correlation coefficient and p-value
print(paste("Correlation Coefficient (Amphibia subset):", correlation_coefficient))
print(paste("P-value (Amphibia subset):", p_value))

# Subset the data for gene_accessions$class == "Amphibia"
subset_data <- subset(genes_accessions, class == "Aves")

# Calculate correlation coefficient and p-value for the subset
correlation_result <- cor.test(subset_data$`Number of Genes`, subset_data$prop_scg)

# Access the correlation coefficient and p-value
correlation_coefficient <- correlation_result$estimate
p_value <- correlation_result$p.value

# Print the correlation coefficient and p-value
print(paste("Correlation Coefficient ( subset):", correlation_coefficient))
print(paste("P-value ( subset):", p_value))

df_ordered_combo$Accession <- sub("^(.*?)_(.*?)_(.*)", "\\1_\\2.\\3", df_ordered_combo$accession)

# Split the 'gene' column by underscores
gene_split <- strsplit(df_ordered_combo$gene, "_")

# Extract the last element from the split strings and assign it to 'stop'
df_ordered_combo$stop <- sapply(gene_split, function(x) tail(x, 1))

# Extract the second to last element from the split strings and assign it to 'start'
df_ordered_combo$start <- sapply(gene_split, function(x) if(length(x) > 1) x[length(x) - 1] else "")

# chromosome_part <- lapply(gene_split, function(x) {
#   if (length(x) > 3) {
#     print(length(x[length(x) - 2]))
#     chromosome <- paste("axolotl_",x[(length(x) - 3)], x[(length(x) - 2)], collapse = "_")
#     return(chromosome)
#   } else if (length(x[length(x) - 2]) == 2 && df_ordered_combo$accession == "GCA_008782695_1") {
#     chromosome <- paste("muntjak_",x[(length(x) - 3)], x[(length(x) - 2)], collapse = "_")
#     return(chromosome)
#   } else if (length(x[length(x) - 2]) == 2 && df_ordered_combo$accession == "GCA_019279795_1") {
#     chromosome <- paste("West_African_lungfish_",x[(length(x) - 3)], x[(length(x) - 2)], collapse = "_")
#     return(chromosome)
#   } else if (length(x[length(x) - 2]) == 2 && df_ordered_combo$accession == "GCA_026652325_1") {
#     chromosome <- paste("Iberian_ribbed_newt_",x[(length(x) - 3)], x[(length(x) - 2)], collapse = "_")
#     return(chromosome)
#   } else if (length(x[length(x) - 2]) == 2 && df_ordered_combo$accession == "GCA_027579735_1") {
#     chromosome <- paste("fire-bellied_toad_",x[(length(x) - 3)], x[(length(x) - 2)], collapse = "_")
#     return(chromosome)
#   } else if (length(x[length(x) - 2]) == 2 && df_ordered_combo$accession == "GCA_028390025_1") {
#     chromosome <- paste("corroboree_frog_",x[(length(x) - 3)], x[(length(x) - 2)], collapse = "_")
#     return(chromosome)
#   } else if (length(x[length(x) - 2]) == 2 && df_ordered_combo$accession == "GCA_029206835_1") {
#     chromosome <- paste("mountain_yellow-legged_frog_",x[(length(x) - 3)], x[(length(x) - 2)], collapse = "_")
#     return(chromosome)
#   } else if (length(x) > 3) {
#     chromosome <- paste(x[(length(x) - 3):(length(x) - 2)], collapse = "_")
#     return(chromosome)
#   } else {
#     return(NA)  # If the pattern doesn't match, assign NA or adjust according to your needs
#   }
# })
chromosome_part <- lapply(seq_along(gene_split), function(i) {
  x <- gene_split[[i]]
  accession <- df_ordered_combo$accession[i]
  # print(x)
  # print(length(x)-2)
  # print(nchar(x[length(x) - 2]))
  if (length(x) > 3 && nchar(x[length(x) - 2]) == 2) {
    if (accession == "GCA_002915635_3") {
      chromosome <- paste("axolotl", x[(length(x) - 3)], x[(length(x) - 2)], sep = "_")
      return(chromosome)
    } else if (accession == "GCA_008782695_1") {
      chromosome <- paste("muntjak", x[(length(x) - 3)], x[(length(x) - 2)], sep = "_")
      return(chromosome)
    } else if (accession == "GCA_019279795_1") {
      chromosome <- paste("West_African_lungfish", x[(length(x) - 3)], x[(length(x) - 2)], sep = "_")
      return(chromosome)
    } else if (accession == "GCA_026652325_1") {
      chromosome <- paste("Iberian_ribbed_newt", x[(length(x) - 3)], x[(length(x) - 2)], sep = "_")
      return(chromosome)
    } else if (accession == "GCA_027579735_1") {
      chromosome <- paste("fire-bellied_toad", x[(length(x) - 3)], x[(length(x) - 2)], sep = "_")
      return(chromosome)
    } else if (accession == "GCA_028390025_1") {
      chromosome <- paste("corroboree_frog", x[(length(x) - 3)], x[(length(x) - 2)], sep = "_")
      return(chromosome)
    } else if (accession == "GCA_029206835_1") {
      chromosome <- paste("mountain_yellow-legged_frog", x[(length(x) - 3)], x[(length(x) - 2)], sep = "_")
      return(chromosome)
    } else {
      chromosome <- paste(x[(length(x) - 3):(length(x) - 2)], collapse = "_")
      return(chromosome)
    }
  } else {
    chromosome <- paste(x[(length(x) - 3):(length(x) - 2)], collapse = "_")
    return(chromosome)
    return(NA)  # If the pattern doesn't match, assign NA or adjust according to your needs
  }
})
# Assign the extracted chromosome information to the 'chromosome' column
df_ordered_combo$chromosome <- unlist(chromosome_part)
#(length(x) > 3 && nchar(x[length(x) - 2]) == 2)
for (i in 1:length(df_ordered_combo$chromosome)) {
  if (length(gene_split[[i]]) > 3 && nchar(gene_split[[i]][length(gene_split[[i]]) - 2]) == 2 &&
      df_ordered_combo$accession[i] == "GCA_002915635_3") {
    # Condition met: Leave the third to last underscore as it is
    # Do nothing; no replacement needed
  } else if (length(gene_split[[i]]) > 3 && nchar(gene_split[[i]][length(gene_split[[i]]) - 2]) == 2 &&
              df_ordered_combo$accession[i] == "GCA_008782695_1"){
  } else if (length(gene_split[[i]]) > 3 && nchar(gene_split[[i]][length(gene_split[[i]]) - 2]) == 2 &&
             df_ordered_combo$accession[i] == "GCA_019279795_1"){
  } else if (length(gene_split[[i]]) > 3 && nchar(gene_split[[i]][length(gene_split[[i]]) - 2]) == 2 &&
             df_ordered_combo$accession[i] == "GCA_026652325_1"){
  } else if (length(gene_split[[i]]) > 3 && nchar(gene_split[[i]][length(gene_split[[i]]) - 2]) == 2 &&
             df_ordered_combo$accession[i] == "GCA_028390025_1"){
  } else if (length(gene_split[[i]]) > 3 && nchar(gene_split[[i]][length(gene_split[[i]]) - 2]) == 2 &&
             df_ordered_combo$accession[i] == "GCA_029206835_1"){
  } else if (length(gene_split[[i]]) > 3 && nchar(gene_split[[i]][length(gene_split[[i]]) - 2]) == 2 &&
             df_ordered_combo$accession[i] == "GCA_027579735_1"){
  }
 else {
    # Replace the third to last underscore with a period
    df_ordered_combo$chromosome[i] <- sub("^(.*?)_(.*?)", "\\1.\\2", df_ordered_combo$chromosome[i])
  }
}

# Read in maximal_range_all.csv as clusters dataframe
clusters <- read.csv("maximal_range_all.csv")
genes_accessions <- read_csv(genes_accessions_input)

# Sort clusters dataframe by accession, chromosome, and start
clusters <- clusters[order(clusters$accession, clusters$chromosome, clusters$start), ]

# Create a column 'cluster_name' in clusters dataframe using the counter for each accession
clusters$cluster_name <- with(clusters, ave(accession, accession, FUN = function(x) paste0(x, "_", seq_along(x))))


# Sort df_ordered_combo dataframe by Accession, chromosome, and start
df_ordered_combo <- df_ordered_combo %>%
  arrange(Accession, chromosome, start)

# Merge clusters with df_ordered_combo based on Accession and chromosome
merged_data <- merge(df_ordered_combo, clusters, by.x = c("Accession", "chromosome"), by.y = c("accession", "chromosome"))

# Convert columns to integers
merged_data$start.x <- as.integer(merged_data$start.x)
merged_data$start.y <- as.integer(merged_data$start.y)
merged_data$stop.x <- as.integer(merged_data$stop.x)
merged_data$stop.y <- as.integer(merged_data$stop.y)

# Filter and find matches based on conditions
filtered_data <- merged_data %>%
  filter(start.y < start.x, start.x < stop.y, start.y < stop.x, stop.x < stop.y)
filtered_data <- unique(filtered_data)

duplicates_start_x <- filtered_data$start.x[duplicated(filtered_data$start.x)]



# Assuming the columns with varying values are col1, col2, and col3
sorted_data <- filtered_data %>%
  arrange(Accession, chromosome, start.x, scg) %>%
  distinct(Accession, chromosome, start.x, .keep_all = TRUE)

# Remove duplicates where all values are the same except possibly scg
final_filtered_data <- sorted_data %>%
  group_by(Accession, chromosome, start.x) %>%
  filter(scg == min(scg)) %>%
  ungroup()

# Group final_filtered_data by cluster_name and scg values, then summarize to count occurrences
summary_clusters <- final_filtered_data %>%
  group_by(cluster_name, scg) %>%
  summarize(count = n()) %>%
  pivot_wider(names_from = scg, values_from = count, values_fill = list(count = 0))

# Create columns scg_123, scg_0, scg_1011 by summing appropriate scg values
summary_clusters <- summary_clusters %>%
  mutate(
    scg_123 = `1` + `2` + `3`,
    scg_0 = `0`,
    scg_1011 = `10` + `11`,
    scg_none = `nope`
  ) %>%
  select(cluster_name, scg_123, scg_0, scg_1011, scg_none)

summary_clusters$total = summary_clusters$scg_123 +  summary_clusters$scg_1011 +  summary_clusters$scg_0+ summary_clusters$scg_none 
summary_clusters$prop = summary_clusters$scg_123/summary_clusters$total
summary_clusters <- summary_clusters %>%
  mutate(accession = sub("_[^_]+$", "", cluster_name))

summary_clusters_merged <- merge(summary_clusters, genes_accessions, by.x = "accession", by.y = "Accession", all.x = TRUE, all.y=FALSE)
summary_clusters_merged$type <- ifelse(summary_clusters_merged$total == 1, 1,
                                       ifelse(summary_clusters_merged$total > 1, 0, ""))

# Define singleton2 as a subset of summary_clusters_merged for type == singleton
singleton2 <- summary_clusters_merged[summary_clusters_merged$type == 1, ]

# Define cluster2 as a subset of summary_clusters_merged for type == cluster
cluster2 <- summary_clusters_merged[summary_clusters_merged$type == 0, ]

pdf("Prop_cluster_by_class.pdf")
ggplot(cluster2, aes(x = class, y = prop)) +
  geom_boxplot() +
  labs(x = "Class", y = "Proportion") +
  ggtitle("Boxplot of Proportion Cluster by Class")
dev.off()

pdf("Prop_singleton_by_class.pdf")
ggplot(singleton2, aes(x = class, y = prop)) +
  geom_boxplot() +
  labs(x = "Class", y = "Proportion") +
  ggtitle("Boxplot of Proportion Singleton by Class")
dev.off()

overall_summary <- summary_clusters %>%
  group_by(accession) %>%
  filter(total > 1) %>%
  summarise(
    proportion_cluster = mean(prop > 0.3)
  )

overall_summary2 <- summary_clusters %>%
  group_by(accession) %>%
  filter(total == 1) %>%
  summarise(
    proportion_singleton = mean(prop > 0.3)
  )

# Group summary_clusters by 'accession' and count the rows where total >= 2
cluster_counts <- summary_clusters %>%
  group_by(accession) %>%
  summarise(clusters = sum(total >= 2, na.rm = TRUE))

singleton_counts <- summary_clusters %>%
  group_by(accession) %>%
  summarise(singletons = sum(total ==1, na.rm = TRUE))

# Merge the cluster counts into overall_summary based on 'accession'
overall_summary <- merge(overall_summary, overall_summary2, by.x = "accession", by.y = "accession", all.x = TRUE, all.y=TRUE)
overall_summary <- merge(overall_summary, cluster_counts, by.x = "accession", by.y = "accession", all.x = TRUE)
overall_summary <- merge(overall_summary, singleton_counts, by.x = "accession", by.y = "accession", all.x = TRUE)


summary_merged <- merge(overall_summary, genes_accessions, by.x = "accession", by.y = "Accession", all.x = TRUE, all.y=FALSE)  #all.y=TRUE

pdf("Prop_cluster_by_class_by_accession.pdf")
ggplot(summary_merged, aes(x = class, y = proportion_cluster)) +
  geom_boxplot() +
  labs(x = "Plotting Clade", y = "Class") +
  ggtitle("Boxplot of Proportion Cluster by Class")
dev.off()

pdf("Prop_singleton_by_class_by_accession.pdf")
ggplot(summary_merged, aes(x = class, y = proportion_singleton)) +
  geom_boxplot() +
  labs(x = "Class", y = "Proportion") +
  ggtitle("Boxplot of Proportion Singleton by Class")
dev.off()

complete_summary <- summary_merged[!is.na(summary_merged$proportion_singleton), ]

cluster_summary <- c(sum((cluster2[cluster2$class == "Mammalia", "scg_123"]))/sum((cluster2[cluster2$class == "Mammalia", "total"])),
                     sum((cluster2[cluster2$class == "Aves", "scg_123"]))/sum((cluster2[cluster2$class == "Aves", "total"])),
                     sum((cluster2[cluster2$class == "Amphibia", "scg_123"]))/sum((cluster2[cluster2$class == "Amphibia", "total"])),
                     sum((cluster2[cluster2$class == "Lepidosauria", "scg_123"]))/sum((cluster2[cluster2$class == "Lepidosauria", "total"])),
                     sum((cluster2[cluster2$class == "Actinopteri", "scg_123"]))/sum((cluster2[cluster2$class == "Actinopteri", "total"])))

singleton_summary <- c(sum((singleton2[singleton2$class == "Mammalia", "scg_123"]))/sum((singleton2[singleton2$class == "Mammalia", "total"])),
                     sum((singleton2[singleton2$class == "Aves", "scg_123"]))/sum((singleton2[singleton2$class == "Aves", "total"])),
                     sum((singleton2[singleton2$class == "Amphibia", "scg_123"]))/sum((singleton2[singleton2$class == "Amphibia", "total"])),
                     sum((singleton2[singleton2$class == "Lepidosauria", "scg_123"]))/sum((singleton2[singleton2$class == "Lepidosauria", "total"])),
                     sum((singleton2[singleton2$class == "Actinopteri", "scg_123"]))/sum((singleton2[singleton2$class == "Actinopteri", "total"])))

taxa <- c("Mammal", "Bird", "Amphibian", "Reptile", "Fish")

# Combine cluster and singleton summaries into a matrix
summary_matrix <- rbind(cluster_summary, singleton_summary)

# Create bar plot
pdf("Summary_cluster_SCO_status.pdf")
barplot(summary_matrix, beside = TRUE, col = c("#D55E00", "#166A53"),
        legend.text = c("Cluster", "Singleton"), xlab = "Taxa", ylab = "Proportion of SCO / total genes",
        names.arg = taxa, main = "Proportion of SCO by Taxa")
dev.off()


duplicated_genes <- df_ordered_tips4_amphibian$gene[duplicated(df_ordered_tips4_amphibian$gene)]
unique(duplicated_genes)



for_gene_family <- merge(summary_clusters_merged, final_filtered_data, by.x = "cluster_name", by.y = "cluster_name", all.x = TRUE, all.y=FALSE)  #all.y=TRUE
for_gene_family <- mutate(for_gene_family, type = as.numeric(type))


gene_family_summary <- for_gene_family %>%
  group_by(gene_family) %>%
  summarise(typ_avg = mean(type, na.rm = TRUE))

pdf("gene_family_summary.pdf")
hist(gene_family_summary$typ_avg, main = "Histogram of typ_avg", xlab = "typ_avg")
dev.off()

# Assuming gene_family_summary contains your summarized data
gene_family_summary$class <- sub("^(.*?)_.*", "\\1", gene_family_summary$gene_family)

# Create boxplots
pdf("gene_family_summary_class.pdf")
boxplot(typ_avg ~ class, data = gene_family_summary, main = "Boxplot of typ_avg by gene_family class", xlab = "Gene Family Class", ylab = "typ_avg")
dev.off()


gene_family_summary <- gene_family_summary %>%
  rowwise() %>%
  mutate(
    scg = for_gene_family$scg[for_gene_family$gene_family == gene_family][1],
    median = median(for_gene_family$prop[for_gene_family$gene_family == gene_family], na.rm = TRUE),
    min = min(for_gene_family$prop[for_gene_family$gene_family == gene_family], na.rm = TRUE)
  )

categorize_scg <- function(value) {
  if (value %in% c(1, 2, 3)) {
    return("1-3")
  } else if (value %in% c(10, 11)) {
    return("10-11")
  } else {
    return("0-nope")
  }
}

# Apply the function to create a new column 'scg_category' in gene_family_summary
gene_family_summary <- gene_family_summary %>%
  mutate(scg_category = sapply(scg, categorize_scg))

# Create a boxplot of median grouped by scg_category
pdf("prop_scg_by_gene_family.pdf")
boxplot(gene_family_summary$median ~ gene_family_summary$scg_category, 
        xlab = "SCG Category", ylab = "Median", 
        main = "Median Prop for Gene Family by SCG Category")
dev.off()


# Create a boxplot for the subsets
pdf("hist_cluster_SCO.pdf")
hist(cluster2$prop, main = "Histogram of prop values", xlab = "prop values", col = "skyblue")
dev.off()

pdf("hist_cluster_SCO_amphibian.pdf")
hist(cluster2$prop[cluster2$class == "Amphibia"], main = "Histogram of prop values for amphibians", xlab = "prop values", col = "lightgreen")
dev.off()

pdf("hist_cluster_SCO_mammal.pdf")
hist(cluster2$prop[cluster2$class == "Mammalia"], main = "Histogram of prop values for mammals", xlab = "prop values", col = "purple")
dev.off()

pdf("hist_cluster_size.pdf")
hist(cluster2$total, main = "Histogram of cluster size", xlab = "# genes", col = "skyblue")
dev.off()

pdf("hist_cluster_size_amphibian.pdf")
hist(cluster2$total[cluster2$class == "Amphibia"], main = "Histogram of cluster size for amphibians", xlab = "# genes", col = "lightgreen")
dev.off()

pdf("hist_cluster_size_mammal.pdf")
hist(cluster2$total[cluster2$class == "Mammalia"], main = "Histogram of cluster size for mammals", xlab = "# genes", col = "purple")
dev.off()
######

taxon <-"Mammalia"

# Filtering 'for_gene_family' dataframe where class == "Mammalia"
subset_gene_family <- for_gene_family[for_gene_family$class == taxon, ]

# Generating cluster sizes from 2 to 150
cluster_size <- 2:max(subset_gene_family$total)

if (taxon %in% c("Mammalia", "Mammal")) {
  prop <- cluster_summary[1]
} else if (taxon %in% c("Aves", "Bird")) {
  prop <- cluster_summary[2]
} else if (taxon %in% c("Amphibia", "Amphibian")) {
  prop <- cluster_summary[3]
} else if (taxon %in% c("Lepidosauria", "Reptile")) {
  prop <- cluster_summary[4]
} else if (taxon %in% c("Actinopteri", "Fish")) {
  prop <- cluster_summary[5]
} else {
  print("Not recognized taxon")
}


# Calculating expected_SCO (20% of cluster size)
expected_SCO <- cluster_size * prop


# Filter 'subset_gene_family' to exclude rows where 'total' is 1
subset_filtered <- subset_gene_family[subset_gene_family$total != 1, ]

# Extract unique cluster sizes from the filtered data
unique_cluster_sizes <- unique(subset_filtered$total)

# Count occurrences of each unique cluster size in the filtered 'total' column
cluster_sizes_frequency <- table(subset_filtered$total)

# Create 'expected2' dataframe with unique cluster sizes and their frequencies
expected2 <- data.frame(cluster_size = as.integer(unique_cluster_sizes),
                        expected_SCO = double(length(unique_cluster_sizes)),
                        freq = integer(length(unique_cluster_sizes)))

# Assign frequency counts to the 'freq' column in 'expected2'
expected2$freq <- as.integer(cluster_sizes_frequency[as.character(expected2$cluster_size)]/expected2$cluster_size)
expected2$freq[is.na(expected2$freq)] <- 0

expected2$zero <- dbinom(0, size = expected2$cluster_size, prob = prop) 
expected2$one <- dbinom(1, size = expected2$cluster_size, prob = prop)

# expected2$exactlyone <- expected2$cluster_size * prop * (1 - prop)^(expected2$cluster_size - 1)

# Calculate the probability of zero or one success
expected2$oneorzero <- (expected2$one + expected2$zero)*expected2$freq
expected2$exactlyone <- expected2$one * expected2$freq
#expected2$twoplus <-  1-expected2$zero - expected2$exactlyone
expected2$twoplus <- (1 - expected2$zero - expected2$one) * expected2$freq
expected2$twoplus_given_one <- expected2$twoplus / (expected2$exactlyone + expected2$twoplus) * expected2$freq

prop_twoplus_given_one <- sum(expected2$twoplus_given_one, na.rm = TRUE)/sum(expected2$freq, na.rm = TRUE)


# Initialize empty columns 'min1' and 'over1' in 'expected2'
expected2$min1 <- integer(nrow(expected2))
expected2$over1 <- integer(nrow(expected2))

# Iterate through each row of 'expected2'
for (i in 1:nrow(expected2)) {
  # Subset 'subset_gene_family' based on 'total == cluster_size' for each row
  subset_data <- subset_gene_family[subset_gene_family$total == expected2$cluster_size[i], ]
  
  # Calculate the count where 'scg_123' >= 1 and save it as 'min1' for each row
  expected2$min1[i] <- length(subset_data$scg_123[subset_data$scg_123 >= 1])/expected2$cluster_size[i]
  
  # Calculate the count where 'scg_123' > 1 and save it as 'over1' for each row
  expected2$over1[i] <- length(subset_data$scg_123[subset_data$scg_123 > 1])/expected2$cluster_size[i]
}

# 
# # Calculate the proportion of clusters with additional Type A genes among clusters with at least one Type A gene
proportion_additional_A_expected <- sum(expected2$over1) / sum(expected2$min1)

# Print the proportion
cat("Actual proportion of clusters with additional SCO among clusters with at least one SCO:", proportion_additional_A_expected, "\n")
cat("Theoretical proportion of clusters with additional SCO among clusters with at least one SCO:", prop_twoplus_given_one, "\n")

# Creating a plot of expected2$cluster_size against twoplus and over1
pdf("actual_vs_expected")
ggplot(expected2, aes(x = cluster_size)) +
  geom_line(aes(y = twoplus, color = "Assuming no bias"), linewidth = 1) +
  geom_line(aes(y = over1, color = "Actual data"), linewidth = 1) +
  labs(title = "Probability of Two or More Successes vs. Cluster Size",
       x = "Cluster Size",
       y = "Probability") +
  scale_color_manual(name = "Lines",
                     values = c("Assuming no bias" = "green", "Actual data" = "blue")) +
  theme_minimal()

dev.off()

# Observed proportion
observed_proportion <- proportion_additional_A_expected

# Expected proportion
expected_proportion <- prop_twoplus_given_one

# Total number of clusters with at least one Type A gene
total_clusters <- sum(expected2$min1)

# Total number of clusters with additional Type A genes among those with at least one Type A gene
observed_additional_A <- sum(expected2$over1)

# Calculate pooled proportion
pooled_proportion <- (observed_additional_A + expected_proportion * total_clusters) / (total_clusters + total_clusters)

# Calculate standard error
SE <- sqrt(pooled_proportion * (1 - pooled_proportion) * (1 / total_clusters + 1 / total_clusters))

# Calculate z-score
z_score <- (observed_proportion - expected_proportion) / SE

# Calculate p-value
p_value <- 2 * pnorm(-abs(z_score))

# Print z-score and p-value
cat("Z-score:", z_score, "\n")
cat("P-value:", p_value, "\n")

position_coordinates_input <- "/Volumes/wengpj01/vertebrate_pipeline/coordinates/position_within_chromosome_all_genes_added_big_1109.csv"
position_coordinates <- read_csv(position_coordinates_input)

for_coordinate <- merge(final_filtered_data, position_coordinates, by.x = c("chromosome", "start.x"), by.y = c("chromosome", "start"))
for_coordinate <- merge(for_coordinate, summary_clusters_merged, by.x = c("cluster_name"), by.y = c("cluster_name"))
for_coordinate_cluster <- for_coordinate %>%
  distinct(cluster_name, .keep_all = TRUE)

cluster3 <- for_coordinate_cluster %>%
  filter(type == 0)

singleton3 <- for_coordinate_cluster %>%
  filter(type == 1)

cor.test(cluster3$ends, cluster3$prop, method = "pearson")
#if ends larger, prop smaller

for_coordinate_cluster %>%
  filter(type == 1) %>%
  summarise(mean_ends = mean(ends, na.rm = TRUE))

for_coordinate_cluster %>%
  filter(type == 1, scg_123 == 1) %>%
  summarise(mean_ends = mean(ends, na.rm = TRUE))

subset_type1_scg0 <- for_coordinate_cluster %>% filter(type == 1, scg_123==0)
subset_type1_scg1 <- for_coordinate_cluster %>% filter(type == 1, scg_123 == 1)

# Performing t-test
t.test(subset_type1_scg0$ends, subset_type1_scg1$ends)

cluster_macro <- for_coordinate_cluster %>% filter(mini == TRUE, type == 0)
cluster_short <- for_coordinate_cluster %>% filter(mini == FALSE, type == 0)
singleton_macro <- for_coordinate_cluster %>% filter(mini == TRUE, type == 1)
singleton_short <- for_coordinate_cluster %>% filter(mini == FALSE, type == 1)

(nrow(cluster_macro)/ (nrow(cluster_macro) + nrow(singleton_macro)))
(nrow(cluster_short)/ (nrow(cluster_short) + nrow(singleton_short)))

amp_singleton_notSCO <-for_coordinate_cluster %>% filter(type == 1, scg_123==0, class=="Amphibia")
mam_singleton_notSCO <-for_coordinate_cluster %>% filter(type == 1, scg_123==0, class=="Mammalia")
bird_singleton_notSCO <-for_coordinate_cluster %>% filter(type == 1, scg_123==0, class=="Aves")
amp_singleton_SCO <-for_coordinate_cluster %>% filter(type == 1, scg_123==1, class=="Amphibia")
mam_singleton_SCO <-for_coordinate_cluster %>% filter(type == 1, scg_123==1, class=="Mammalia")
bird_singleton_SCO <-for_coordinate_cluster %>% filter(type == 1, scg_123==1, class=="Aves")

amp_cluster_notSCO <-for_coordinate_cluster %>% filter(type == 0, scg_123==0, class=="Amphibia")
mam_cluster_notSCO <-for_coordinate_cluster %>% filter(type == 0, scg_123==0, class=="Mammalia")
bird_cluster_notSCO <-for_coordinate_cluster %>% filter(type == 0, scg_123==0, class=="Aves")
amp_cluster_SCO <-for_coordinate_cluster %>% filter(type == 0, scg_123>=1, class=="Amphibia")
mam_cluster_SCO <-for_coordinate_cluster %>% filter(type == 0, scg_123>=1, class=="Mammalia")
bird_cluster_SCO <-for_coordinate_cluster %>% filter(type == 0, scg_123>=1, class=="Aves")

amp_singleton <-for_coordinate_cluster %>% filter(type == 1, class=="Amphibia")
amp_cluster <-for_coordinate_cluster %>% filter(type == 0, class=="Amphibia")
mam_singleton <-for_coordinate_cluster %>% filter(type == 1, class=="Mammalia")
mam_cluster <-for_coordinate_cluster %>% filter(type == 0, class=="Mammalia")
bird_singleton <-for_coordinate_cluster %>% filter(type == 1, class=="Aves")
bird_cluster <-for_coordinate_cluster %>% filter(type == 0, class=="Aves")

singleton_SCO <- for_coordinate_cluster %>% filter(type == 1, scg_123==1)
singleton_notSCO <- for_coordinate_cluster %>% filter(type == 1, scg_123==0)
cluster_SCO <- for_coordinate_cluster %>% filter(type == 0, scg_123>=1)
cluster_notSCO <- for_coordinate_cluster %>% filter(type == 0, scg_123==0)



mean(amp_singleton_SCO$ends)
mean(mam_singleton_SCO$ends)
mean(amp_singleton_notSCO$ends)
mean(mam_singleton_notSCO$ends)


t.test(amp_singleton_SCO$ends, amp_singleton_notSCO$ends)


t.test(mam_singleton_SCO$ends, mam_singleton_notSCO$ends)

t.test(bird_singleton_SCO$ends, bird_singleton_notSCO$ends)

t.test(amp_cluster_SCO$ends, amp_cluster_notSCO$ends)
t.test(bird_cluster_SCO$ends, bird_cluster_notSCO$ends)
t.test(mam_cluster_SCO$ends, mam_cluster_notSCO$ends)

common_wood_file <-"/Users/kwhiggins27/Desktop/OneDrive - Harvard University/Weng Lab/Amphibian_paper/NZ_orth/common_wood_annotated_genes_to_plot.csv"
common_wood <- read_csv(common_wood_file)
colnames(common_wood)[4] = "chromosome"
colnames(common_wood)[5] = "stop"

common_wood_merged <- merge(common_wood, final_filtered_data, by.x = c("chromosome", "stop"), by.y = c("chromosome", "stop.x"))

nrow(common_wood_merged %>%
  filter(scg %in% c(1, 2, 3), `1:1` == 1))

nrow(common_wood_merged %>%
       filter(scg %in% c(1, 2, 3), `1:many` == 1))

nrow(common_wood_merged %>%
       filter(scg %in% c(1, 2, 3), `many:many` == 1))

nrow(common_wood_merged %>%
       filter(scg %in% c(1, 2, 3), `1:0` == 1))

nrow(common_wood_merged %>%
       filter(scg %in% c(1, 2, 3), `many:0` == 1))
###
nrow(common_wood_merged %>%
       filter(scg %in% c(0), `1:1` == 1))

nrow(common_wood_merged %>%
       filter(scg %in% c(0), `1:many` == 1))

nrow(common_wood_merged %>%
       filter(scg %in% c(0), `many:many` == 1))

nrow(common_wood_merged %>%
       filter(scg %in% c(0), `1:0` == 1))

nrow(common_wood_merged %>%
       filter(scg %in% c(0), `many:0` == 1))

####

human_file <- "/Users/kwhiggins27/Desktop/OneDrive - Harvard University/Weng Lab/Amphibian_paper/NZ_orth/human_marmoset_annotated_genes_to_plot.csv"
human <- read.csv(human_file)

human_merged <- merge(human, final_filtered_data, by.x = c("stop"), by.y = c("stop.x"))


# Assuming genes_accessions is your data frame
subset_genes <- subset(genes_accessions, `Number of Genes` > 1)

# Calculate the proportion
proportion <- length(intersect(unique(subset_genes$Accession), unique(cluster2$accession))) / length(unique(subset_genes$Accession))
print(proportion)

sum(cluster2$total, na.rm = TRUE) / (sum(cluster2$total, na.rm = TRUE) + sum(singleton2$total, na.rm = TRUE))

sum(singleton2$scg_123)/ (sum(singleton2$scg_123) + sum(cluster2$scg_123))
