##################################################################################################
#
# Figure 1. Performance - total mass and growth rate
#
# Growth rate figure is adapted from:
# https://datadrumstick.wordpress.com/2016/09/03/chick-weight-vs-diet-a-case-for-one-way-anova/
# [accessed 20191212]
#
##################################################################################################

# packages
library(tidyverse)  # version 1.3.0
library(ggpubr)     # version 0.2.4
library(cowplot)    # version 1.0.0
library(ggsci)      # version 2.9



#read in data
mass <- read_csv("data/zootechnical/bird_mass.csv") %>%
                mutate(diet = if_else(diet == "gos", "GOS", diet))

# repeated measure birds for growth rate calculation
rm_birds <- c("e_102", "e_104", "e_109", "e_112", "e_113", "e_117", "e_119", 
              "e_124", "e_176", "e_173", "e_72", "e_81", "e_82", "e_92", "e_95", 
              "e_96", "e_223", "e_209", "e_199", "e_210")

# plot growth rate
grate.p <- mass %>% 
            filter(sample_title %in% rm_birds | trial == "sup") %>% 
            group_by(diet, age) %>% 
            summarise(median = median(total_mass), 
                      sem = sd(total_mass)/sqrt(length(total_mass))) %>% 
            ggplot(mapping = aes(x = as.numeric(age), 
                                 y = median, shape = diet, color = diet)) + 
                geom_line(aes(linetype = diet, color = diet), size=1) +
                geom_errorbar(aes(ymin = median - sem, ymax = median + sem)) +
                geom_point() +
                scale_linetype_manual(values = c(ctl = "solid", 
                                                 GOS = "twodash", 
                                                 target = "dashed")) +
                scale_colour_manual(values = c("#000000", "#E69F00", "grey59")) +
                theme_bw() +
                theme(plot.margin = unit(c(0, 0, 0, 0), "mm"),
                      legend.position = c(0.875, 0.095),
                      axis.line = element_blank(),
                      panel.border = element_rect(fill = NA,
                                    colour = "black", linetype = 1, size = 0.7),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      axis.text = element_text(size = 10, colour = "black"),
                      axis.title = element_text(size = 10),
                      legend.title = element_blank(),
                      legend.background = element_rect(
                                          linetype = 1, size = 0.5, colour = 1),
                      legend.text = element_text(size = 10), 
                      legend.margin = margin(c(0.7,2,1,1), unit = "mm"),
                      legend.key = element_blank(),
                      plot.title = element_text(size = 12, face = "bold", hjust = 0)) +
            labs(title= "A Trial 1", x = "Age (days)", y = "Bird mass (g)") +
            scale_y_continuous(limits = c(0,3000), expand = c(0,0)) +
            scale_x_continuous(breaks = seq(0, 35, 5), limits = c(0,35)) 



final.tot.mass <- mass %>% filter(age == 35 & trial != "sup")
  
p.T1 <- final.tot.mass %>% 
              filter(trial == "T1") %>% 
              ggboxplot(. , y = "total_mass", x = "diet", 
                        fill = "diet", palette = c("#ffffff", "#E69F00"), 
                        add = "jitter", outlier.shape = NA) +
                geom_hline(yintercept = 2144, linetype = "dashed") +
                theme(plot.margin = unit(c(0, 0, 0, 0), "mm"),
                      legend.position = "none",
                      panel.border = element_rect(linetype = 1, 
                                                size = 0.7, fill = NA),
                      axis.line = element_blank(), 
                      axis.text = element_text(size = 10),
                      axis.title = element_text(size = 10),
                      plot.title = element_text(size = 12, face="bold")) +
                labs(title="B Trial 1", x = "", y = "Bird mass (g)") +
                stat_compare_means(method = "t.test", size = 3.5, 
                                   label = "p.format", label.x = 1.3) +
                scale_y_continuous(breaks = seq(1500, 3500, 1000), 
                                   limits = c(1400,3400)) + 
                annotate("text",  x = 2.07, y = 1510, label = "ctl = 1.478", size = 3.5) +
                annotate("text", x = 1.945, y = 1400, label = "GOS = 1.455", size = 3.5)



p.T2 <- final.tot.mass %>% 
              filter(trial == "T2") %>% 
              arrange(diet) %>% 
              ggboxplot(. , y = "total_mass", x = "diet", 
                        fill = "diet", palette = c("#ffffff", "#E69F00"), 
                        add = "jitter", outlier.shape = NA) +
                geom_hline(yintercept = 2144, linetype = "dashed") +
                theme(plot.margin = unit(c(0, 0, 0, 0), "mm"),
                      legend.position = "none",
                      panel.border = element_rect(linetype = 1, 
                                                    size = 0.7, fill = NA),
                      axis.line = element_blank(), 
                      axis.text.y = element_text(size = 10),
                      axis.title.y = element_text(size = 10),
                      axis.text.x = element_blank(), 
                      plot.title = element_text(size = 12, face="bold")) +
                labs(title="C Trial 2", x = "", y = "Bird mass (g)") +
                stat_compare_means(method = "t.test", size = 3.5,
                                   label = "p.format", label.x = 1.2, label.y = 3450) +
                scale_y_continuous(breaks = seq(1500, 3500, 1000), limits = c(1400,3600)) + 
                annotate("text",  x = 2.05, y = 1700, label = "ctl = 1.536", size = 3.5) +
                annotate("text", x = 1.93, y = 1470, label = "GOS = 1.379", size = 3.5)




p.T3 <- final.tot.mass %>% 
              filter(trial == "T3") %>% 
              ggboxplot(. , y = "total_mass", x = "diet", add = "jitter",
                      fill = "diet", palette = c("#ffffff", "#E69F00"), 
                      outlier.shape = NA  ) +
                  geom_hline(yintercept = 2144, linetype = "dashed") +
                  theme(plot.margin = unit(c(0, 0, 0, 0), "mm"),
                        legend.position = "none",
                        panel.border = element_rect(linetype = 1, 
                                                      size = 0.7, fill = NA),
                        axis.line = element_blank(), 
                        axis.text.y = element_text(size = 10),
                        axis.text.x = element_text(size = 10),
                        axis.title.y = element_text(size = 10),
                        plot.title = element_text(size = 12, face = "bold")) +
                        labs(title="D Trial 3", x = "", y = "Bird mass (g)") +
                  stat_compare_means(method = "t.test", size = 3.5, 
                                     label = "p.format", label.y = 3450, 
                                     label.x = 1.2) +
                  scale_y_continuous(breaks = seq(1500, 3500, 1000), 
                                                        limits = c(1400,3600)) + 
                  annotate("text",  x = 2.05, y = 1700, label = "ctl = 1.379", size = 3.5) +
                  annotate("text", x = 1.93, y = 1470, label = "GOS = 1.341", size = 3.5)




# print .tiff
tiff(filename="results/figures/fig1_performance.tiff", 
                       width=200, height=100, units="mm", res=300)

ggdraw() +
  draw_plot(grate.p, x = 0, y = 0, width = 0.45, height = 1) +
  draw_plot(p.T1, x = 0.46, y = 0.025, width = 0.25, height = 0.97) +
  draw_plot(p.T2, x = 0.72, y = 0.526, width = 0.25, height = 0.47) +
  draw_plot(p.T3, x = 0.72, y = 0.02, width = 0.25, height = 0.52)
  
dev.off()




