library(tidyverse)
library(ggthemes)
library(data.table)
library(ggrastr)
library(vegan)
library(egg)

pal <- tableau_color_pal("Tableau 20")
col <- pal(20)

## define function
modify_format <- function(d, thre = 1E-5, pseud = 1E-4){
  keep <- apply(d, 2, mean) > thre & apply(d > 0, 2, mean) > 0.1
  d <- d[, keep]
  d %>% dim() %>% print()
  return(d)
}

plot_envfit <- function(d, count, e, title){
  set.seed(1)
  
  ## prepare metadata
  md <- data.frame(count = count)
  md$bacteroides <- ifelse(e$enterotype == "Bact_Phoc", 1, 0)
  md$prevotella <- ifelse(e$enterotype == "Prevotella", 1, 0)
  md$firmicutes <- ifelse(e$enterotype == "Firmicutes", 1, 0)  
  ##
  
  ## mds
  d.min <- min(d[d != 0])/2
  d.dist <- vegdist(log10(d + d.min), method = "euclidean")
  d.mds <- metaMDS(d.dist)
  
  ## envfit
  df <- data.frame(x = d.mds$points[, 1], y = d.mds$points[, 2], count = count, enterotype = e$enterotype)
  d.envfit <- envfit(d.mds, md)
  
  ## prepare labels
  df.cont <- data.frame(scores(d.envfit, "vectors")) * ordiArrowMul(d.envfit) 
  df.cont <- df.cont %>% rownames_to_column()
  df.cont$rowname <- df.cont$rowname %>% 
    str_replace("enterotype", "") %>% 
    str_replace("count", "Microbial load") %>% 
    str_replace("bacteroides", "Bacteroides type") %>% 
    str_replace("prevotella", "Prevotella type") %>% 
    str_replace("firmicutes", "Firmicutes type") %>%     
    str_replace(" ", "\n")
  df.cont
  df.cont <- df.cont[2:4, ]
  
  
  ## plot result    
  fold <- (max(d.mds$points[, 2]) / max(df.cont$NMDS2)) / 2
  
  col <- colorRampPalette(c(col[1], col[2], "cornsilk", col[11]))(100)
  p <- ggplot(df, aes(x = x, y = y)) +
    theme_classic() +
    ggtitle(title) +
    xlab("MDS1") +
    ylab("MDS2") +
    geom_vline(xintercept = 0, linetype = 2, col = "gray", size = 0.5) +
    geom_hline(yintercept = 0, linetype = 2, col = "gray", size = 0.5) +
    geom_point_rast(data = df, aes(col = count), size = 0.3) +    
    scale_color_gradientn(colours = col) +
    theme(legend.position = "right") +
    theme(legend.text = element_text(size = 8), legend.key.size = unit(0.2, "cm")) +
    theme(legend.title = element_blank())
  p
  
  p2 <- p + geom_segment(data = df.cont, aes(x = 0, y = 0, xend = NMDS1*fold, yend = NMDS2*fold), col = "tomato", arrow = arrow(length = unit(0.2, "cm"))) +
    geom_text(data = df.cont, aes(x = NMDS1*fold, y = NMDS2*fold, label = rowname), size = 2)
  
  return(p2)
}

## read data
d1 <- fread("data/GALAXY_mOTUs_v25.tsv") %>% data.frame(check.names = F) %>% column_to_rownames("V1")
d2 <- fread("data/MetaCardis_mOTUs_v25.tsv") %>% data.frame(check.names = F) %>% column_to_rownames("V1")

d1 <- modify_format(d1)
d2 <- modify_format(d2)

md1 <- read.delim("data/GALAXY_load.tsv")
md2 <- read.delim("data/MetaCardis_load.tsv")

e1 <- read.delim("data/Galaxy_enterotypes.tsv", header = T)
e2 <- read.delim("data/MetaCardis_enterotypes.tsv", header = T)

d <- d1
e <- e1
count <- log10(md1$count)
md <- md1
title <- "GALAXY"

## mds + envfit
p1 <- plot_envfit(d1, log10(md1$count), e1, "GALAXY/MicrobLiver")
p2 <- plot_envfit(d2, log10(md2$count), e2, "MetaCardis")

## save plots
p <- ggarrange(p1, p2, ncol = 2)
ggsave(p, filename = "out/figure_1A.pdf", width = 7, height = 2.5)


