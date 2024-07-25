#Activate libraries
library(tidyverse)
library(scales)

#Read dataset (csv format)
data_dr <- readr::read_csv("data_DRcurve_ver240516.csv")

name_convert <- c("Axo_1" = "axolotl_54", "Bul_1" = "bullfrog_51","Bul_4" = "bullfrog_61", 
                  "Can_1" = "cane_58", "Can_3" = "cane_56", "Can_4" = "cane_54", 
                  "Xen_2" = "clawed_20", "Xen_6" = "clawed_23")
data_dr$gene <-  stringr::str_replace_all(data_dr$gene, name_convert)

#Set minimum and maximum values for y axis 
plot_limit_y <- tibble(
  well = "empty",
  compound = c("AflatoxinB1",
               "Cinobufagin",
               "AflatoxinB1",
               "AflatoxinB1",
               "AflatoxinB1",
               "AflatoxinB1",
               "AflatoxinB1",
               "AflatoxinB1"),
  concentration = c(rep(0.1, 8)),
  gene = c("axolotl_54",
           "bullfrog_51",
           "bullfrog_61",
           "cane_58",
           "cane_56",
           "cane_54",
           "clawed_20",
           "clawed_23"),
  group = "DRcurve",
  AUC = 0,
  file = "empty",
  AUCc = c(6, 0.9, 7, 1.5, 4, 0.9, 0.9, 2.5)
  )

#Merge dataset with limits for y axis
data_dr_plot <- bind_rows(
  mutate(data_dr, alpha = 1),
  mutate(plot_limit_y, alpha = 0)
  )

#plot lineplots
line_DR <- data_dr_plot %>%
  ggplot(aes(x = concentration, y = AUCc, colour = compound, fill = compound, shape = compound, linetype = compound)) +
    stat_summary(data = data_dr_plot[data_dr_plot$alpha == 1,], fun.data = mean_se, geom = "errorbar",width = 0.05, linetype = "solid") +
    stat_summary(data = data_dr_plot[data_dr_plot$alpha == 1,], fun = mean, geom = "line", size = 0.7) +
    stat_summary(data = data_dr_plot[data_dr_plot$alpha == 1,], fun = mean, geom = "point", size = 2.5) +
    stat_summary(data = data_dr_plot[data_dr_plot$alpha == 0,], fun = mean, geom = "point", size = 3, alpha = 0) +
    scale_x_log10(breaks = 10^(-3:1), labels = trans_format("log10", math_format(10^.x))) +
    scale_color_manual(values = c("#D56127", "#136A54", "#9374B4")) +
    scale_fill_manual(values = c("#D56127", "#136A54", "#9374B4")) +
    scale_shape_manual(values = c(21, 23, 24, 25, 22)) +
    scale_linetype_manual(values = c(1, 1, 1, 1, 1)) +
    xlab("concentration (mM)") +
    ylab("Receptor activity (AUC x10^4)") +
    ylim(0, NA) +
    facet_wrap(facets = "gene", ncol = 4, scales = "free") +
    theme_classic() +
    theme(legend.position = "none")

plot(line_DR)