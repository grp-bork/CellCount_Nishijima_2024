library(tidyverse)
library(ggthemes)
library(ggrastr)
library(data.table)
library(egg)
library(scales)
library(ggpubr)

pal <- tableau_color_pal("Superfishel Stone")
col <- pal(10)

d1 <- read_rds("out/rds/external_validation.GALAXY.mOTUs.rds")
d2 <- read_rds("out/rds/external_validation.MetaCardis.mOTUs.rds")

md1 <- fread("data/GALAXY_load.tsv")
md2 <- fread("data/MetaCardis_load.tsv")

d2$cohort <- md1$cohort
d1$cohort <- md2$cohort

## GALAXY/MicrobLiver
p1 <- ggplot(d2, aes(x = count, y = pred)) +
  theme_bw() +
  geom_abline(linetype = 2, col = "gray") +  
  geom_point_rast(stroke = 0, col = col[4], size = 0.8, alpha = 0.8) +
  stat_cor(size = 3) +
  geom_smooth(method = "lm") +
  xlab("Microbial load") +
  ylab("Predicted value") +
  ggtitle("GALAXY/MicrobLiver") +
  facet_wrap(~cohort, ncol = 6)


## MetaCardis
p2 <- ggplot(d1, aes(x = count, y = pred)) +
  theme_bw() +
  geom_abline(linetype = 2, col = "gray") +  
  geom_point_rast(stroke = 0, col = col[6], size = 0.8, alpha = 0.8) +
  stat_cor(size = 3) +
  geom_smooth(method = "lm") +
  xlab("Microbial load") +
  ylab("Predicted value") +
  ggtitle("MetaCardis") +
  facet_wrap(~cohort, ncol = 6)



## merge
p <- ggarrange(p1, p2, heights = c(2, 1.2), ncol = 1)
ggsave(p, filename = "out/supplementary_figure7.pdf", width = 12, height = 7)
