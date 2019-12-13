##############################################################################
#
# Figure 2. Differential OTUs
#
# including:
# OTU-level stack barcharts
# LEfSE outputs
# absolute abundance
#
##############################################################################

# packages
library(tidyverse)  # version 1.3.0
library(cowplot)    # version 1.0.0
library(ggpubr)     # version 0.2.3


# read in OTUs
OTU <- read_tsv("data/mothur/broiler_gos.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.shared") %>% 
        select(-c(label, numOtus))  %>% 
        mutate(row_sum = rowSums(select(., starts_with("Otu")))) %>%  
        mutate_if(is.numeric, funs(100*(. / row_sum))) %>% 
        select (-row_sum) %>% 
        gather(OTU, value, -Group) %>% 
        spread(Group, value)

# calculate n
#names(OTU) %>% 
#            as_tibble() %>% 
#            filter(value != "OTU") %>% 
#            rename("sample" = value) %>% 
#            mutate(diet = word(sample, 1, sep = fixed('_')),
#                   age = as.numeric(word(sample, 2, sep = fixed('_')))) %>% 
#            group_by(diet, age) %>% 
#            summarise(seq_n = n())


# read in taxonomy
tax <- read_tsv("data/mothur/broiler_gos.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.cons.taxonomy")
  
# merge datasets
# drop controls
df <- inner_join(tax, OTU, by = "OTU") %>% 
                      select(-c("gos_8_1", "ctl_35_8"))
    
# LEFSE
lefse_08 <- "data/mothur/broiler_gos.0.03.day08.0.03.filter.0.03.subsample.0.03.lefse_summary"
lefse_15 <- "data/mothur/broiler_gos.0.03.day15.0.03.filter.0.03.subsample.0.03.lefse_summary"
lefse_22 <- "data/mothur/broiler_gos.0.03.day22.0.03.filter.0.03.subsample.0.03.lefse_summary"
lefse_35 <- "data/mothur/broiler_gos.0.03.day35.0.03.filter.0.03.subsample.0.03.lefse_summary"

# absolute abundance
abund <- read_tsv("data/zootechnical/abs_abundance.txt") %>% 
                filter(trial == "T1" & seq_rep != "gos_8_1") %>% 
                mutate(L_crispatus_copy_no = if_else(is.na(L_crispatus_copy_no),
                                       2.36, L_crispatus_copy_no))



# n
#abund %>% group_by(diet, age) %>% summarise(seq_n = n())

# 'tidy' format abundance data
abund_tidy <- abund %>% 
                select(age, diet, L_johnsonii_copy_no, L_crispatus_copy_no) %>% 
                group_by(diet, age) %>% 
                gather(key = "target", value = "count", -c(diet, age)) %>% 
                group_by(diet, age, target)



# 'tidy' format otu data
df_otu <- df %>% 
            gather("sample", "RA", -c(OTU, Size, Taxonomy)) %>% 
            separate(sample, c("diet", "age", "replicate"), "_", remove = FALSE) %>%
            mutate(genus = word(Taxonomy, 6, sep = fixed(';'))) %>%
            select(c(sample, diet, age, RA, OTU, genus))


## OTUs

# list top 10 most abundant OTUs
top_10 <- df_otu %>% 
            group_by(OTU) %>% 
            summarise(mean_RA = mean(RA)) %>% 
            arrange(desc(mean_RA)) %>% 
            head(10) %>% 
            pull(OTU)


# subset top10 OTUs
df_otu_top10 <- df_otu %>%
                      subset(OTU %in% top_10) %>% 
                      mutate(Genus = paste0(str_remove_all(genus, "[()01569]"), " (", OTU, ")")) %>% 
                      select(-c(OTU, genus))

  
# compile non-top 10 OTUs as 'other'
df_otu_top10_other <- df_otu %>% 
                          subset(OTU %in% top_10) %>% 
                          group_by(sample, diet, age) %>% 
                          summarize(other = 100 - sum(RA)) %>% 
                          gather(other, RA, -c(sample, diet, age)) %>% 
                          rename(Genus = "other") %>% 
                          bind_rows(df_otu_top10)

# OTU sanity check!
# df_otu_top10_other %>% group_by(sample) %>% summarise(sum = sum(RA))

# summarize labels
ids <- c("ctl_8_1"="ctl_1", "ctl_8_2"="ctl_2", "ctl_8_3"="ctl_3", 
         "ctl_8_4"="ctl_4", "ctl_8_5"="ctl_5", "ctl_8_6"="ctl_6", 
         "ctl_8_7"="ctl_7", # n = 7
         "gos_8_2"="GOS_2", "gos_8_3"="GOS_3", "gos_8_4"="GOS_4", 
         "gos_8_5"="GOS_5", "gos_8_6"="GOS_6", "gos_8_7"="GOS_7", # n = 6
         "ctl_15_1"="ctl_1", "ctl_15_2"="ctl_2", "ctl_15_3"="ctl_3", 
         "ctl_15_4"="ctl_4", "ctl_15_5"="ctl_5", "ctl_15_6"="ctl_6", 
         "ctl_15_7"="ctl_7", # n = 7
         "gos_15_1"="GOS_1", "gos_15_2"="GOS_2", "gos_15_3"="GOS_3", 
         "gos_15_4"="GOS_4", "gos_15_5"="GOS_5", "gos_15_6"="GOS_6", # n = 6
         "ctl_22_1"="ctl_1", "ctl_22_2"="ctl_2", "ctl_22_3"="ctl_3", 
         "ctl_22_4"="ctl_4", "ctl_22_5"="ctl_5", "ctl_22_6"="ctl_6", 
         "ctl_22_7"="ctl_7", # n = 7
         "gos_22_1"="GOS_1", "gos_22_2"="GOS_2", "gos_22_3"="GOS_3", 
         "gos_22_4"="GOS_4", "gos_22_5"="GOS_5", "gos_22_6"="GOS_6", 
         "gos_22_7"="GOS_7", # n = 7
         "ctl_35_1"="ctl_1", "ctl_35_2"="ctl_2", "ctl_35_3"="ctl_3", 
         "ctl_35_4"="ctl_4", "ctl_35_5"="ctl_5", "ctl_35_6"="ctl_6", 
         "ctl_35_7"="ctl_7", # n =7 
         "gos_35_1"="GOS_1", "gos_35_2"="GOS_2", "gos_35_3"="GOS_3", 
         "gos_35_4"="GOS_4", "gos_35_5"="GOS_5", "gos_35_6"="GOS_6") # n = 6


# list OTU ID order
otu_list <- c("other", "Lactobacillus (Otu0010)", "Butyricicoccus (Otu0009)",
              "Clostridium_IV (Otu0008)","Lachnospiraceae_unclassified (Otu0007)",
              "Lactobacillus (Otu0006)","Lachnospiraceae_unclassified (Otu0005)",
              "Lachnospiraceae_unclassified (Otu0004)", 
              "Lachnospiraceae_unclassified (Otu0003)", 
              "Enterobacteriaceae_unclassified (Otu0002)",  "Faecalibacterium (Otu0001)" )


# list OTU colours
otu_colours <- c("#F0A3FF", "#add8e6", "#9DCC00", "#993F00", "#8F7C00","#C20088",
                 "#94FFB5","#808080","#FFCC99","#003380", "#2BCE48")


# plot OTU barchart
otu.p <- df_otu_top10_other %>%
              mutate(Genus = factor(Genus, levels = otu_list)) %>%
              ggplot(aes(sample, RA, fill = Genus)) +
                geom_bar(stat = "identity") + 
                facet_grid(. ~ as.numeric(age), drop = TRUE, scale = "free", space = "free_x") +
                scale_fill_manual(values = otu_colours, breaks = otu_list) +
                ylab("Relative abundance (%)") + 
                scale_y_continuous(expand = c(0,0)) + 
                scale_x_discrete(labels = ids) +
                theme_bw() + 
                theme(strip.background = element_rect(fill = "white"), 
                      panel.spacing = unit(0.3, "lines"), 
                      text = element_text(size = 7.5),
                      axis.text.x = element_text(angle = 90, hjust = 1, 
                                      vjust = 0.5, colour = "black"), 
                      axis.text.y = element_text(colour = "black") ) +
                guides(fill = guide_legend(reverse = TRUE))


## families

# summarise at family level
# tidy format family data
df_fam <- df %>% 
            gather("sample", "RA", -c(OTU, Size, Taxonomy)) %>% 
            separate(sample, c("diet", "age", "replicate"), "_", remove = FALSE) %>%
            mutate(Family = word(Taxonomy, 5, sep = fixed(';'))) %>%
            select(c(sample, diet, age, RA, Family)) %>% 
            group_by(sample, diet, age, Family) %>% 
            summarise(sum_RA = sum(RA))


# sanity check!
#df_fam %>% group_by(sample) %>% summarise(sum = sum(sum_RA))


# subset to only families >1% RA
# rename families (strip bootstrap values from ID)
df_fam_sub <- df_fam  %>% 
                  filter(sum_RA > 1) %>% 
                  mutate(Family = paste0(str_remove_all(Family, "[()019]")))


# compile non-top 10 families as 'other'
df_fam_sub_other <- df_fam_sub %>% 
                            group_by(sample, diet, age) %>% 
                            summarize(other = 100 - sum(sum_RA)) %>% 
                            gather(other, sum_RA, -c(sample, diet, age)) %>% 
                            rename("Family" = other) %>% 
                            bind_rows(df_fam_sub)


# list family order
fam_list <- c("other", "Ruminococcaceae", "Peptostreptococcaceae", 
              "Lactobacillaceae", "Lachnospiraceae", "Firmicutes_unclassified", 
              "Erysipelotrichaceae", "Enterococcaceae", "Enterobacteriaceae", 
              "Coriobacteriaceae", "Clostridiales_unclassified", 
              "Bacteria_unclassified", "Bacillales_unclassified")
                 

# list family colours  
fam_colours <- c("#F0A3FF", "#2BCE48", "#993F00", "#8F7C00","#add8e6", "#94FFB5", 
                 "#4f4e4e", "#C20088", "#003380", "#ff7400", "#9DCC00", "#FFCC99", 
                 "#808080")


# plot family barchart
fam.p <- df_fam_sub_other %>% 
            mutate(Family = factor(Family, levels = fam_list)) %>%
            ggplot(aes(sample, sum_RA, fill = Family)) +
              geom_bar(stat = "identity") + 
              facet_grid(. ~ as.numeric(age), drop = TRUE, scale = "free", space = "free_x") +
              scale_fill_manual(values = fam_colours, breaks = fam_list) +
              ylab("Relative abundance (%)") +
              scale_y_continuous(expand = c(0,0)) + 
              scale_x_discrete(labels = ids) +
              theme_bw() + 
              theme(strip.background = element_rect(fill = "white"), 
                    panel.spacing = unit(0.3, "lines"), 
                    text = element_text(size = 7.5),
                    axis.text.x = element_text(angle = 90, hjust = 1, 
                                         vjust = 0.5, colour = "black"), 
                    axis.text.y = element_text(colour = "black")  ) +
              guides(fill = guide_legend(reverse = TRUE, override.aes = list(size = 0.3)))




### LEfSE ###

# bind lefse data together
df_lefse <- bind_rows(read_tsv(lefse_08) %>% mutate(age = 8),
                         read_tsv(lefse_15) %>% mutate(age = 15),
                         read_tsv(lefse_22) %>% mutate(age = 22),
                         read_tsv(lefse_35) %>% mutate(age = 35)) %>% 
            drop_na(LDA)


# format LEfSE data
# add OTU IDs
df_lefse_format <- df_lefse %>% 
                      inner_join(tax, by = "OTU") %>% 
                      mutate(Genus = paste0(str_remove_all(
                                      word(Taxonomy, 6, sep = fixed(';')), 
                                      "[()01569]"), 
                                      " (", OTU, ")"),
                             Class = if_else(str_sub(Class, 1, 3) == "gos", "GOS",
                                             str_sub(Class, 1, 3)))
  

# plot LEfSE data
lefse.p <- df_lefse_format %>% 
                  ggplot(aes(x = Genus, y = LDA, fill = Class)) +
                  geom_bar(stat = "identity", position = "dodge", colour = "black") + 
                  scale_fill_manual(values = c("white", "#E69F00")) +
                  theme_bw() +
                  theme(text = element_text(size = 7, colour="#000000"),
                        axis.text.y = element_text(size = 5),
                        strip.background = element_rect(fill = "#ffffff", color = "#000000"),
                        panel.grid.major.y = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.grid.major.x = element_line(size = 0.2, colour = "grey50") ) +
                  scale_y_continuous(breaks=seq(0, 5, 1), limits=c(0,5), expand=c(0,0)) +
                  ylab(expression("LDA score" ~ (log[10]))) + 
                  facet_grid(age ~ .) +
                  coord_flip()




### abs abundance ###

# stats...
lc_stat <- abund_tidy %>% 
              summarise(count = list(count)) %>% 
              filter(target == "L_crispatus_copy_no") %>% 
              spread(diet, count) %>% 
              group_by(age, target) %>% 
              mutate(wilcox = wilcox.test(unlist(ctl),unlist(gos))$p.value)

lj_stat <- abund_tidy %>% 
              summarise(count = list(count)) %>% 
              filter(target == "L_johnsonii_copy_no") %>% 
              spread(diet, count) %>% 
              group_by(age, target) %>% 
              mutate(wilcox = wilcox.test(unlist(ctl),unlist(gos))$p.value)


# plot L. crispatus abundance
cris.p <- ggboxplot(abund_tidy %>% filter(target == "L_crispatus_copy_no"), 
                    x = "age", y = "count", add = "jitter", outlier.shape = NA,
                    fill = "diet", palette = c("#ffffff", "#E69F00")) + 
                    theme(legend.position="right", 
                          panel.background = element_blank(), 
                          panel.grid.major = element_blank(), 
                          panel.grid.minor = element_blank(),
                          legend.title = element_blank(),
                          text = element_text(size = 7),
                          legend.background = element_rect(linetype = 1, 
                                                           size = 0.3, colour = 1)) +
                     ylab("Otu0006 genome" ~ (log[10])) + xlab("Age (days)") +
                     annotate(geom="text", x=c(1,2,3,4), y=c(9.75,9.75,9.75,9.75), 
                              label=c(paste("p ==",round(lc_stat$wilcox[1],3)), 
                                      paste("p ==",round(lc_stat$wilcox[2],3)),
                                      paste("p ==",round(lc_stat$wilcox[3],3)),
                                      paste("p ==",round(lc_stat$wilcox[4],3))), 
                          size = 2.1, parse=T) +
                    scale_y_continuous(breaks = seq(2, 10, 2), limits = c(2, 10), 
                                       expand = c(0,0))
                              

# plot L. johnsonii abundance
Ljohn.p <- ggboxplot(abund_tidy %>% filter(target == "L_johnsonii_copy_no"), 
                     x = "age", y = "count", add = "jitter", outlier.shape = NA,
                     fill = "diet", palette=c("#ffffff", "#E69F00")) + 
                     theme(legend.position="right", 
                           panel.background = element_blank(), 
                           panel.grid.major = element_blank(), 
                           panel.grid.minor = element_blank(),
                           legend.title = element_blank(),
                           text = element_text(size = 7),
                           legend.background = element_rect(linetype = 1, 
                                                           size = 0.3, colour = 1)  ) +
                      ylab("Otu0010 genome" ~ (log[10])) + xlab("Age (days)") +
                      scale_y_continuous(breaks = seq(2, 10, 2), 
                                         limits = c(2, 10), expand = c(0,0)) +
                      annotate(geom="text", x=c(1,2,3,4), y=c(9.75,9.75,9.75,9.75), 
                               label=c(paste("p ==",round(lj_stat$wilcox[1],3)), 
                                       paste("p ==",round(lj_stat$wilcox[2],3)),
                                       paste("p ==",round(lj_stat$wilcox[3],3)),
                                       paste("p ==",round(lj_stat$wilcox[4],3))), 
                               size = 2.1, parse=T)



abund.p <- plot_grid(cris.p + theme(legend.position="none"),
                     Ljohn.p + theme(legend.position="none"),
                     align = 'h',
                     ncol = 2)



# get legend
legend <- get_legend(lefse.p)


# print plot
tiff(filename="results/figures/fig2_differential_otu.tiff", 
     width = 315, height = 170, units = "mm", res = 300)

ggdraw() +
  draw_plot(otu.p, x = 0.02, y = 0.5, width = 0.605, height = 0.47) +
  draw_plot(fam.p, x = 0.02, y = 0, width = 0.553, height = 0.47) +
  draw_plot(lefse.p + theme(legend.position = "none"), 
            x = 0.63, y = 0.5, width = 0.32, height = 0.48) +
  draw_plot(abund.p, x = 0.63, y = 0, width = 0.37, height = 0.48) +
  draw_plot_label(c("A", "B", "C", "D"), c(0, 0, 0.61, 0.61), c(1, 0.5, 1, 0.48), size = 12) +
  draw_grob(legend, x = 0.875, y = 0.51, width = 0.2, height = 0.5)


dev.off()











