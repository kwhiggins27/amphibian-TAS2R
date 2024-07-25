#Activate libraries
library(tidyverse)

#Read dataset
data_sc <- readr::read_csv("data_BitterAssay_ver240513.csv")

name_convert <- c("Axo_1" = "axolotl_54", "Axo_2" = "axolotl_10", "Axo_3" = "axolotl_42",
                  "Bul_1" = "bullfrog_51", "Bul_2" = "bullfrog_12", "Bul_3" = "bullfrog_17",
                  "Bul_4" = "bullfrog_61", "Can_1" = "cane_58", "Can_2" = "cane_17",
                  "Can_3" = "cane_56", "Can_4" = "cane_54", "Ter_1" = "dart_18",
                  "Ter_3" = "dart_15", "Xen_1" = "clawed_31", "Xen_2" = "clawed_20",
                  "Xen_3" = "clawed_02", "Xen_5" = "clawed_18", "Xen_6" = "clawed_23")
data_sc$gene <-  stringr::str_replace_all(data_sc$gene, name_convert)

#lists of ligands and genes
ligand_list <- unique(data_sc$compound)
order_compound <- c("AflatoxinB1", "Batrachotoxin", "Cinobufagin", "Heliotrine", "Marinobufagenin", 
                    "Swainsonine", "a-thujone", "Amarogentin", "Arbutin", "Aristolochic_acid", 
                    "Camphor", "Chloramphenicol", "Chloroquine", "Colchicine", "Coumarin", 
                    "Denatonium_benzoate", "Diphenylthiourea", "Genistein", "Helicin", "Papaverine", 
                    "Picrotoxin", "PROP", "PTC", "Quinine", "Salicin", 
                    "Strychnine", "Xanthotoxin", "Yohimbine", "Buffer")

gene_list <- unique(data_sc$gene) %>%
  sort()

#Set minimum and maximum values for y axis 
plot_limit_y <- tibble(
  well = "empty",
  compound = "Buffer",
  concentration = rep(0, length(gene_list)),
  gene = gene_list,
  group = "empty",
  AUC = c(10e+4, 1.3e+4, 10e+4, 
          5e+4, 3e+4, 2e+4, 12e+4, 
          5e+4, 3e+4, 8e+4, 3e+4, 
          5e+4, 2e+4, 3e+4, 3e+4, 2e+4, 
          12e+4, 2e+4),
  file = "empty"
  )

#Merge dataset with limits for y axis
data_sc_plot <- bind_rows(
  mutate(data_sc, alpha = 1),
  mutate(plot_limit_y, alpha = 0)
  ) %>%
  transform(., compound = factor(compound, levels = order_compound)) %>%
  mutate(., AUCx4 = AUC/10000)

#Plot barplots
bar_screening <- ggplot(data_sc_plot[data_sc_plot$alpha == "1",], aes(x = compound, y = AUCx4, fill = compound)) +
    stat_summary(fun = mean, geom = "bar", alpha = 0.8, colour = "black", width = 0.7) +
    geom_point(size = 1, fill = "#888888", colour = "#888888", shape = 21) +
    stat_summary(fun.data = mean_se, geom = "errorbar", colour = "black", width = 0.2) +
    geom_point(data = data_sc_plot[data_sc_plot$alpha == "0",], size = 1.2, fill = "black", colour = "black", shape = 21, alpha = 0) +
    geom_hline(yintercept = 0) +
    scale_y_continuous(expand = c(0,0), limits = c(0, NA)) +
    scale_fill_manual(values = c(rep("#ff8888", 6), rep("#99ccff", 22), "#ffffff")) +
    ylab("Receptor activity (AUC x10^4)") +
    theme_classic() +
    facet_wrap(facets = "gene", ncol = 3, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.line.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = "none")

plot(bar_screening)