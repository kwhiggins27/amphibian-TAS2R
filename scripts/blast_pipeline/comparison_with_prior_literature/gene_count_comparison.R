library(tidyverse)

# Read dataset
## accession of assembly which crashed the gene identification pipeline.
special_characters <- read.delim("accession_special_character.txt", header = F) %>%
  as.matrix()

## read gene count data 
data_kat <- read_csv("gene_count_Kate2024.csv") %>%
  dplyr::filter(!(Accession %in% special_characters))
data_pol <- read_csv("gene_count_Policarpo2024.csv")

## read busco analysis results
data_busco <- read_csv("../../../results/coordinate_analysis/plot_generation/busco_full_data.csv")  %>%
  dplyr::filter(!(Accession %in% special_characters))


# Make dataset
data_common <- inner_join(data_kat, data_pol[-(1:3)], by="Accession") %>%
  mutate(diff_vs_complete7tm = Complete7tm - genes_Kate,
         diff_vs_completeTotal = Complete_Total - genes_Kate,
         diff_category_vs_complete7tm = ifelse(Complete7tm > genes_Kate, "pol>kat",
                                       ifelse(Complete7tm == genes_Kate, "pol=kat", "pol<kat")),
         diff_category_vs_completeTotal = ifelse(Complete_Total > genes_Kate, "pol>kat",
                                                 ifelse(Complete_Total == genes_Kate, "pol=kat", "pol<kat"))) %>%
  inner_join(., data_busco[-(2:3)], by="Accession") %>%
    mutate(busco_category = ifelse(used_busco >= 0.9, "BUSCO90",
                                  ifelse(used_busco >= 0.8, "BUSCO80",
                                          ifelse(used_busco >= 0.7, "BUSCO70", "under70"))))


# Summary
table(data_common$diff_category_vs_complete7tm)
table(data_common$diff_category_vs_completeTotal)


# Regression
LM_vs_complete7tm <- lm(formula = Complete7tm ~ genes_Kate, data = data_common)
summary(LM_vs_complete7tm)

LM_vs_complete_Total <- lm(formula = Complete_Total ~ genes_Kate, data = data_common)
summary(LM_vs_complete_Total)


# Scatter plot
col_clade <- c("#223D7C", "#116B54", "#5A236A", "#B4DBAE",
               "#B4E3F5", "#aaaaaa", "#9EBADD", "#30AB58",
               "#333333", "#D46027", "#AC82DF", "#7D25D9")

data_common %>%
  ggplot(aes(x = genes_Kate,
             y = Complete7tm, 
             color = PlottingClade)) +
  geom_point(size = 2, mapping = aes(shape = busco_category)) +
  xlim(0, 220) +
  ylim(0, 220) +
  labs(x = "TAS2R number identified in this study",
       y = "TAS2R number identified in Policarpo et al. 2024",
       colour = "Clade") +
  scale_color_manual(values = col_clade) +
  geom_smooth(method = "lm", 
              formula = 'y~x', 
              se = FALSE,
              color = "black", 
              linewidth = 0.7) +
  theme_classic()

ggsave("comparison_vs_policarpo2024_complete7TM_cor.pdf", width = 12, height = 6.75)


# Scatter plot to present differences of gene counts between studies
data_common %>%
  ggplot(aes(x = latin,
             y = diff_vs_complete7tm,
             color = PlottingClade)) +
  geom_point(size = 2, mapping = aes(shape = busco_category)) +
  labs(x = "Species (alphabetical order)",
       y = "Differences in TAS2R counts between two studies",
       colour = "Clade") +
  scale_color_manual(values = col_clade) +
  scale_y_continuous(limits = c(-20, 20), 
                     breaks = seq(-20, 30, 10)) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggsave("comparison_vs_policarpo2024_complete7TM_diff.pdf", width = 12, height = 6)
