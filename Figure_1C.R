library(tidyverse)
library(ggthemes)
library(data.table)
library(psych)
library(ggforestplot)

pal <- tableau_color_pal("Superfishel Stone")
col <- pal(10)

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
d1 <- fread("data/GALAXY_mOTUs_v25.tsv") %>% data.frame(check.names = F) %>% column_to_rownames("V1")
d2 <- fread("data/MetaCardis_mOTUs_v25.tsv") %>% data.frame(check.names = F) %>% column_to_rownames("V1")

d1 <- modify_format(d1)
d2 <- modify_format(d2)

md1 <- read.delim("data/GALAXY_load.tsv")
md2 <- read.delim("data/MetaCardis_load.tsv")



## ext shared species
keep1 <- colnames(d1) %in% colnames(d2)
keep2 <- colnames(d2) %in% colnames(d1)

d1 <- d1[, keep1] %>% data.frame(check.names = F)
d2 <- d2[, keep2] %>% data.frame(check.names = F)
##




## calc correlation
cor1 <- psych::corr.test(d1, log10(md1$count), method = "pearson")
cor2 <- psych::corr.test(d2, log10(md2$count), method = "pearson")
df1 <- data.frame(cor = cor1$r, cohort = "GALAXY/MicrobLiver")
df2 <- data.frame(cor = cor2$r, cohort = "MetaCardis")
df1$name <- rownames(df1)
df2$name <- rownames(df2)

df <- rbind(df1, df2)

## modify species name
df$name <- df$name %>% 
  str_replace("ref_mOTU_v25_", "") %>% 
  str_replace("meta_mOTU_v25_", "") %>% 
  str_replace("species incertae sedis", "sp.")
rownames(df) <- c()


## boxplot for top and bottom 20 species
keep <- !df$name %>% str_detect("Shannon|Simpson|richness")
df.div <- df[!keep, ]
div.names <- df$name[!keep] %>% unique()
df <- df[keep, ]

top20 <- tapply(df$cor, df$name, mean) %>% sort() %>% names() %>% head(n = 20)
tail20 <- tapply(df$cor, df$name, mean) %>% sort() %>% names() %>% tail(n = 20)

keep <- df$name %in% c(top20, tail20)
df2 <- rbind(df[keep, ], df.div)
df2$cate <- ifelse(str_detect(df2$name, "Shannon|Simpson|richness"), "Diversity", "Speceis")

## plot
p <- ggplot(df2, aes(y = reorder(name, -cor, median), x = cor)) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3) +
  xlab("Pearson correlation") +
  geom_stripes(odd = "#80808010", even = "#00000000") +
  ylab("") +
  scale_color_manual(values = col[c(4, 6)]) +
  theme(legend.position = "right", legend.title = element_blank()) +
  geom_point(aes(col = cohort), size = 1) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  facet_grid(~cate, space = "free", scales = "free") +
  theme(legend.position = c(0.9, 0.8), legend.key.size = unit(0.7, "line")) +
  
  theme(strip.text.x = element_blank()) +
  coord_flip()
p
ggsave(p, filename = "out/figure_1C.pdf", width = 10, height = 4)

## save intermediate results
saveRDS(df, file = "out/rds/cor_df.rds")
saveRDS(cor1, file = "out/rds/cor1.rds")
saveRDS(cor2, file = "out/rds/cor2.rds")
