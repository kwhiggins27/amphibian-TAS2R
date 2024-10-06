library(tidyverse)

#Read dataset
data_busco <- read_tsv("../../results/coordinate_analysis/plot_generation/busco5_5_0_human_zebrafish.tsv")
data_position <- read_csv("../../results/coordinate_analysis/plot_generation/variation_bitter_receptors_within_species.csv") %>%
  left_join(., data_busco[-2], by = "accession")


#Gene position among human genome
assembly_order_level_human <- c("GCA_009914755.4",
                               "GCA_024586135.1", "GCA_011064465.2", "GCA_018873775.2", "GCA_021951015.1", "GCA_021950905.1", "GCA_014905855.1", "GCA_016695395.2", "GCA_016700455.2", "GCA_020497085.1", "GCA_020497115.1", "GCA_015476435.1",
                               "GCA_001292825.2", "GCA_000306695.2", "GCA_000002115.2", "GCA_000002125.2", "GCA_000212995.1", "GCA_001712695.1")

color_chr_human <- c("#D46027", "#30AB58", "#AC82DF") 

data_position %>%
  dplyr::filter(species == "human") %>%
  select(-strain) %>%
  na.omit() %>%
  transform(., accession = factor(accession, levels = assembly_order_level_human)) %>%
  transform(., Chr = factor(Chr, levels = c("Chr5", "Chr7", "Chr12"))) %>%
  ggplot(aes(x = accession, y = position, color = Chr, group = gene_name)) +
  geom_point(size = 2, aes(shape = assembly_level)) +
  geom_line() +
  ylim(0, 0.5) +
  coord_flip() +
  scale_color_manual(values = color_chr_human) +
  theme_classic()

ggsave("position_human.pdf", width = 12, height = 6.75)


data_position %>%
  dplyr::filter(species == "human") %>%
  select(-strain) %>%
  na.omit() %>%
  transform(., accession = factor(accession, levels = assembly_order_level_human)) %>%
  transform(., Chr = factor(Chr, levels = c("Chr5", "Chr7", "Chr12"))) %>%
  ggplot(aes(x = contig_N50, y = position, color = accession)) +
  geom_smooth(method = "lm", color = "#555555", se = FALSE) +
  geom_point(size = 1.5, aes(shape = assembly_level)) +
  ylim(0, 0.5) +
  scale_x_log10() +
  facet_wrap(Chr~gene_name, ncol = 6) +
  theme_bw() +
  theme(strip.background = element_blank())

ggsave("position_human_vs_contigN50.pdf", width = 15, height = 10)


data_position %>%
  dplyr::filter(species == "human") %>%
  select(-strain) %>%
  na.omit() %>%
  transform(., accession = factor(accession, levels = assembly_order_level_human)) %>%
  transform(., Chr = factor(Chr, levels = c("Chr5", "Chr7", "Chr12"))) %>%
  ggplot(aes(x = Complete, y = position, color = accession)) +
  geom_smooth(method = "lm", color = "#555555", se = FALSE) +
  geom_point(size = 1.5, aes(shape = assembly_level)) +
  ylim(0, 0.5) +
  xlim(94, 98) +
  facet_wrap(Chr~gene_name, ncol = 6) +
  theme_bw() +
  theme(strip.background = element_blank())

ggsave("position_human_vs_busco.pdf", width = 15, height = 10)


#Gene position among zebrafish genome
assembly_order_level_zf <- c("GCA_033170195.2", "GCA_018400075.1", "GCA_020184715.1", "GCA_944039275.1",
                             "GCA_008692375.1", "GCA_903684855.2", "GCA_903684865.1", "GCA_000002035.4",
                             "GCA_903798185.1", "GCA_903798175.1", "GCA_903798165.1")

color_chr_zf <- c("#D46027", "#30AB58") 

data_position %>%
  dplyr::filter(species == "zebrafish") %>%
  transform(., accession = factor(accession, levels = assembly_order_level_zf)) %>%
  ggplot(aes(x = accession, y = position, color = Chr, group = gene_name)) +
  geom_point(size = 2, aes(shape = strain)) +
  geom_line() +
  ylim(0, 0.5) +
  coord_flip() +
  scale_color_manual(values = color_chr_zf) +
  theme_classic()

ggsave("position_zebrafish.pdf", width = 12, height = 6.75)

data_position %>%
  dplyr::filter(species == "zebrafish") %>%
  transform(., accession = factor(accession, levels = assembly_order_level_zf)) %>%
  ggplot(aes(x = contig_N50, y = position, color = accession)) +
  geom_smooth(method = "lm", color = "#555555", se = FALSE) +
  geom_point(size = 1.5, aes(shape = assembly_level)) +
  ylim(0, 0.5) +
  facet_wrap(Chr~gene_name, ncol = 3) +
  theme_bw() +
  theme(strip.background = element_blank())

ggsave("position_zebrafish_vs_contigN50.pdf", width = 8, height = 5.9)


data_position %>%
  dplyr::filter(species == "zebrafish") %>%
  transform(., accession = factor(accession, levels = assembly_order_level_zf)) %>%
  ggplot(aes(x = Complete, y = position, color = accession)) +
  geom_smooth(method = "lm", color = "#555555", se = FALSE) +
  geom_point(size = 1.5, aes(shape = assembly_level)) +
  xlim(75, 100) +
  ylim(0, 0.5) +
  facet_wrap(Chr~gene_name, ncol = 3) +
  theme_bw() +
  theme(strip.background = element_blank())

ggsave("position_zebrafish_vs_busco.pdf", width = 8, height = 5.9)


#Liner regression analysis
data_position_lm <- data_position %>%
  group_by(gene_name) %>%
  nest() %>%
  mutate(lm_vs_busco = map(data, ~lm(position~Complete, data = .)),
         lm_vs_contigN50 = map(data, ~lm(position~contig_N50, data = .))
  )
