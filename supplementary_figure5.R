library(tidyverse)
library(clusterProfiler)
library(DOSE)
library(ggridges)
library(colorspace)

## define function
enrichment_analysis <- function(df, title){
  set.seed(1)
  res <- list()
  
  ## read pathway file
  md <- read.delim("data/pathway.list", header = F, stringsAsFactors = F)
  colnames(md) <- c("ko", "map", "path")
  md$path <- as.factor(md$path)
  
  ## extract only KOs in df
  keep <- which(md$ko %in% df$ko)
  
  ## proportion of detected KOs in each pathway
  prop <- table(md$path[keep]) / table(md$path)
  path.list <- names(table(md$path))[prop > 0.5]
  
  ## use only pathway with 50% KOs
  keep <- which(md$path %in% path.list)
  md <- md[keep, ]
  keep <- which(df$ko %in% md$ko)
  df <- df[keep, ]
  
  ## prepare input files for GSEA
  cor <- as.numeric(df$num)
  names(cor) <- df$ko
  
  cor <- rev(sort(cor))
  disease2gene <- md[, c(2, 1)]
  disease2name <- md[, c(2, 3)]
  
  ## GSEA
  y <- GSEA(cor, TERM2GENE = disease2gene, TERM2NAME = disease2name, minGSSize = 5, pvalueCutoff = 0.2)
  df <- y@result[, 1:8]
  
  p <- ggplot(df, aes(x = enrichmentScore, y = reorder(Description, enrichmentScore, mean), size = -log10(p.adjust), col = -log10(p.adjust))) +
    theme_bw() +
    geom_point() +
    ggtitle(title) +
    theme(axis.title.y = element_blank()) +
    labs(x = "Enrichment score") +
    geom_vline(xintercept = 0, linetype = 2, col = "gray") +
    scale_color_continuous_sequential(palette = "Reds", l1 = 20, c2 = 70, p1 = 1)
  p
  return(p)
}

df1 <- read_rds("out/rds/cor1.ko.rds")
df2 <- read_rds("out/rds/cor2.ko.rds")
p1 <- enrichment_analysis(df1, "GALAXY/MicrobLiver")
p2 <- enrichment_analysis(df2, "MetaCardis")

## merge plots
p <- egg::ggarrange(p1, p2, ncol = 2)
ggsave(p, filename = "out/supplementary_figure5.pdf", width = 12, height = 5)


