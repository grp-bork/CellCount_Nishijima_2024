library(tidyverse)
library(ggthemes)
library(ggrastr)
library(egg)
library(scales)
library(ggpubr)

pal <- tableau_color_pal("Superfishel Stone")
col <- pal(10)


d1.int <- read_rds("out/rds/internal_validation.GALAXY.KO.rds")
d1.ext <- read_rds("out/rds/external_validation.GALAXY.KO.rds")
d2.int <- read_rds("out/rds/internal_validation.MetaCardis.KO.rds")
d2.ext <- read_rds("out/rds/external_validation.MetaCardis.KO.rds")
d3.int <- read_rds("out/rds/internal_validation.Vandeputte.rds")
d3.ext <- read_rds("out/rds/external_validation.Vandeputte.rds")


## Internal validation
## GALAXY/MicrobLiver (internal)
p1.int <- ggplot(d1.int, aes(x = count, y = pred)) +
  theme_classic() +
  ggtitle("GALAXY/MicrobLiver") +
  xlab("Microbial load") +
  ylab("Predicted value") +
  stat_cor() +
  geom_abline(linetype = 2, col = "gray") +
  geom_point_rast(size = 0.1, col = col[4]) +
  geom_smooth(method = "lm", se = T) 
p1.int


## MetaCardis (internal)
p2.int <- ggplot(d2.int, aes(x = count, y = pred)) +
  theme_classic() +
  ggtitle("MetaCardis") +
  xlab("Microbial load") +
  ylab("Predicted value") +    
  stat_cor() +  
  geom_abline(linetype = 2, col = "gray") +
  geom_point_rast(size = 0.1, col = col[6]) +
  geom_smooth(method = "lm", se = T) 
p2.int


## External validation
## GALAXY/MicrobLiver (external)
p1.ext <- ggplot(d1.ext, aes(x = count, y = pred)) +
  theme_classic() +
  ggtitle("GALAXY/MicrobLiver") +
  xlab("Microbial load") +
  ylab("Predicted value") +
  stat_cor() +
  geom_abline(linetype = 2, col = "gray") +
  geom_point_rast(size = 0.1, col = col[4]) +
  geom_smooth(method = "lm", se = T) 
p1.ext


## MetaCardis (internal)
p2.ext <- ggplot(d2.ext, aes(x = count, y = pred)) +
  theme_classic() +
  ggtitle("MetaCardis") +
  xlab("Microbial load") +
  ylab("Predicted value") +    
  stat_cor() +  
  geom_abline(linetype = 2, col = "gray") +
  geom_point_rast(size = 0.1, col = col[6]) +
  geom_smooth(method = "lm", se = T) 
p2.ext



## Vandeputte 16S
p3.int <- ggplot(d3.int, aes(x = count, y = pred)) +
  theme_classic() +
  ggtitle("Vandeputte 2021") +
  xlab("Microbial load") +
  ylab("Predicted value") +
  stat_cor() +
  geom_abline(linetype = 2, col = "gray") +
  geom_point_rast(size = 0.1, col = col[2]) +
  geom_smooth(method = "lm", se = T) 
p3.int

p3.ext <- ggplot(d3.ext, aes(x = count, y = pred)) +
  theme_classic() +
  ggtitle("Vandeputte 2021") +
  xlab("Microbial load") +
  ylab("Predicted value") +    
  stat_cor() +  
  geom_abline(linetype = 2, col = "gray") +
  geom_point_rast(size = 0.1, col = col[2]) +
  geom_smooth(method = "lm", se = T) 
p3.ext


n <- ggplot() + theme_minimal()
p <- egg::ggarrange(p1.int, p2.int, p1.ext, p2.ext, p3.int, p3.ext, ncol = 4)
ggsave(p, filename = "out/supplementary_figure6.pdf", width = 12, height = 6)
