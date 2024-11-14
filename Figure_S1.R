

library(tidyverse)
library(ggthemes)
library(data.table)
library(ggrastr)
library(clusterProfiler)
library(colorspace)

pal <- tableau_color_pal("Superfishel Stone")
col <- pal(10)

d1 <- read.delim("data/GALAXY_load.tsv")
d2 <- read.delim("data/MetaCardis_load.tsv")

## compare histogram between GALAXY and MetaCardis (Figure S1B)
df1 <- data.frame(count = d1$count, data = "GALAXY/MicrobLiver")
df2 <- data.frame(count = d2$count, data = "MetaCardis")
df <- rbind(df1, df2)

p_hist <- ggplot(df, aes(x = log10(count), y = ..density.., fill = data, color = data)) +
  theme_bw() +
  geom_histogram(alpha = 0, position = "identity", bins = 50) +
  scale_x_log10() +
  xlab("Microbial load (log10)") +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values = col[c(4, 6)]) +
  scale_color_manual(values = col[c(4, 6)]) +  
  theme(legend.position = "right", legend.title = element_blank())
p_hist
ggsave(p_hist, filename = "out/Figure_S1B.pdf", width = 4.5, height = 2.5)



## compare microbial load associations (Figure S1C and D)
## speceis
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
p <- egg::ggarrange(p.sp, p.ko, ncol = 2)
ggsave(p, filename = "out/Figure_S1CD.pdf", width = 7, height = 3.5)


## save intermediate results
df1 <- data.frame(num = cor1$r, ko = str_replace(rownames(cor1$r), "\\:.*", ""), cohort = "GALAXY/MicrobLiver")
df2 <- data.frame(num = cor2$r, ko = str_replace(rownames(cor2$r), "\\:.*", ""), cohort = "MetaCardis")
saveRDS(df1, file = "out/rds/cor1.ko.rds")
saveRDS(df2, file = "out/rds/cor2.ko.rds")



## microbial load-species correlations at the genus level (Figure S1E)
df <- read_rds("out/rds/cor_df.rds")
df <- df %>% filter(name != "Shannon diversity") %>% filter(name != "Simpson diversity") %>% filter(name != "Species richness")
df$id <- df$name %>% str_replace(".*\\[", "\\") %>% str_replace("\\]", "")

## get GTDB taxonomy
tax <- read.delim("data/motus2GTDB.txt", header = T, stringsAsFactors = F)

gtdb.genus <- c()
for(i in 1:nrow(df)){
  print(i)
  i.id <- df$id[i] %>% as.numeric()
  temp <- tax %>% filter(mOTU == i.id)  
  
  if(nrow(temp) == 1){
    
    ## add taxonomic names to unclassified genus 
    if(grepl("\\d", temp$Genus) & !grepl("\\d", temp$Family)){
      temp$Genus <- paste0(temp$Genus, " [", temp$Family, "]")
    }else{
      if(grepl("\\d", temp$Genus) & !grepl("\\d", temp$Order)){
        temp$Genus <- paste0(temp$Genus, " [", temp$Order, "]")
      }else{
        if(grepl("\\d", temp$Genus) & !grepl("\\d", temp$Class)){
          temp$Genus <- paste0(temp$Genus, " [", temp$Class, "]")
        }
      }
    }
    
    gtdb.genus[i] <- temp$Genus
    print(temp$Genus)
  }else{
    gtdb.genus[i] <- NA
  }
}
df$genus <- gtdb.genus

## use only genus with more than 4
minimum_num <- 4 
keep.genus <- names(table(df$genus))[table(df$genus) > minimum_num]
keep <- df$genus %in% keep.genus
df <- df[keep, ]

p <- ggplot(df, aes(y = reorder(genus, -cor, median), x = cor, fill = cohort)) +
  theme_classic() +
  ggforestplot::geom_stripes(odd = "#80808010", even = "#00000000") +  
  geom_vline(xintercept = 0, linetype = 2) +    
  geom_boxplot(outlier.size = 0.5, size = 0.3) +
  xlab("Pearson correlation") +
  ylab("") +
  scale_fill_manual(values = col[c(4, 6)]) +
  theme(strip.text.x = element_blank()) +  
  theme(legend.position = "top", legend.title = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face="italic")) +
  theme(legend.position = c(0.92, 0.85), legend.key.size = unit(0.7, "line")) +
  coord_flip()
p
ggsave(p, filename = "out/Figure_S1E.pdf", width = 8, height = 3, useDingbats=FALSE)




## KEGG KO enrichment analysis (Figure S1F)
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
ggsave(p, filename = "out/Figure_S1F.pdf", width = 12, height = 5)