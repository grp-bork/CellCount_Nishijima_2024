library(tidyverse)
library(ggthemes)
library(ggrastr)
library(egg)
pal <- tableau_color_pal("Superfishel Stone")
col <- pal(10)

## species
cor1 <- read_rds("out/rds/cor1.rds")
cor2 <- read_rds("out/rds/cor2.rds")

df.all <- data.frame(x = cor1$r, y = cor2$r)

p.sp <- ggplot(df.all, aes(x = x, y = y)) +
  theme_classic() +
  xlab("Microbial load-species correlation\n(GALAXY/MicrobLiver)") +
  ylab("Microbial load-species correlation\n(MetaCardis)") +  
  ggpubr::stat_cor() +
  geom_point_rast(alpha = 0.5, size = 1, stroke = 0, col = col[4]) +
  geom_smooth(method = "lm", se = T) 
p.sp



## KO

## function for filtering
modify_format <- function(d, thre = 1E-5, pseud = 1E-4){
  shannon <- vegan::diversity(d)
  simpson <- vegan::diversity(d, index = "simpson")
  richness <- apply(d > 0, 1, sum)
  
  keep <- apply(d, 2, mean) > thre & apply(d > 0, 2, mean) > 0.1
  d <- d[, keep]
  d <- d %>% select(-contains("-1"))
  
  d <- log10(d + pseud)
  
  d$`Shannon diversity` <- shannon
  d$`Simpson diversity` <- simpson
  d$`Species richness` <- richness
  
  d %>% dim() %>% print()
  return(d)
}

## read data
d1 <- fread("data/GALAXY_KO.tsv") %>% data.frame(check.names = F) %>% column_to_rownames("V1")
d2 <- fread("data/MetaCardis_KO.tsv") %>% data.frame(check.names = F) %>% column_to_rownames("V1")

d1 <- modify_format(d1, 1E-6, 1E-6)
d2 <- modify_format(d2, 1E-6, 1E-6)

d1 <- d1 %>% select(-contains("diversity")) %>% select(-contains("richness"))
d2 <- d2 %>% select(-contains("diversity")) %>% select(-contains("richness"))

md1 <- read.delim("data/GALAXY_load.tsv")
md2 <- read.delim("data/MetaCardis_load.tsv")



## ext shared KO
keep1 <- colnames(d1) %in% colnames(d2)
keep2 <- colnames(d2) %in% colnames(d1)

d1 <- d1[, keep1] %>% data.frame(check.names = F)
d2 <- d2[, keep2] %>% data.frame(check.names = F)
##


## calc correlation
cor1 <- psych::corr.test(d1, log10(md1$count), method = "pearson")
cor2 <- psych::corr.test(d2, log10(md2$count), method = "pearson")

df.all <- data.frame(x = cor1$r, y = cor2$r)

p.ko <- ggplot(df.all, aes(x = x, y = y)) +
  theme_classic() +
  xlab("Microbial load-KO correlation\n(GALAXY/MicrobLiver)") +
  ylab("Microbial load-KO correlation\n(MetaCardis)") +  
  ggpubr::stat_cor() +
  geom_point_rast(alpha = 0.5, size = 1, stroke = 0, col = col[4]) +
  geom_smooth(method = "lm", se = T) 
p.ko


## merge plots
p <- ggarrange(p.sp, p.ko, ncol = 2)
ggsave(p, filename = "out/supplementary_figure3.pdf", width = 7, height = 3.5)


## save intermediate results
df1 <- data.frame(num = cor1$r, ko = str_replace(rownames(cor1$r), "\\:.*", ""), cohort = "GALAXY/MicrobLiver")
df2 <- data.frame(num = cor2$r, ko = str_replace(rownames(cor2$r), "\\:.*", ""), cohort = "MetaCardis")
saveRDS(df1, file = "out/rds/cor1.ko.rds")
saveRDS(df2, file = "out/rds/cor2.ko.rds")

