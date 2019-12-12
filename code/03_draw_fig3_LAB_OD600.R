########################################################################
#
# Figure 3 In vitro growth of lactobacilli on galacto-oligosaccharides
#
########################################################################

# packages 
library(tidyverse)  # version 1.3.1
library(ggpubr)     # version 0.2.4

# read in data
df <- read_csv("data/zootechnical/od600.csv")

# format
od_mean <- df %>% 
          group_by(organism, replicate, medium) %>% 
          summarise(mean = mean(OD600)) 


od_blank <- od_mean %>% 
                spread(replicate, mean) %>% ungroup() %>% 
                mutate_at(.vars = vars(a, b, c), .funs = funs(. - blank)) %>% 
                select(-blank)


od_nc <- od_blank %>% 
          gather("replicate", "od", -organism, -medium) %>% 
          spread(medium, -organism) %>% 
          mutate_at(.vars = vars(`GOS 0.5%`, `glucose 0.5%`), .funs = funs(. - `No additional carbon (NC)`)) %>% 
          select(-`No additional carbon (NC)`) %>% 
          gather("diet", "od", -replicate, -organism)

od.p <- ggboxplot(od_nc, x = "organism", y = "od",
                  color = "diet", palette =c("#000000", "#E69F00"),
                  add = "jitter", shape = "diet",
                  ylab = "OD (600)", ylim = c(0, 0.5)) +
                  theme(legend.position="right",
                        panel.border = element_blank()) +
                  rotate_x_text(45)


tiff("results/figures/fig3_invitro_lab.tiff", width = 140, height = 100, units = "mm", res = 300)

od.p

dev.off()

 

