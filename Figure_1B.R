library(tidyverse)
library(ggthemes)
library(data.table)
library(ggrastr)
library(ggpubr)

col <- c("#6d2f20", "#b75347", "#df7e66", "#e09351", "#edc775", "#94b594", "#224b5e")

## read data
md1 <- read.delim("data/GALAXY_load.tsv")
md2 <- read.delim("data/MetaCardis_load.tsv")

e1 <- read.delim("data/Galaxy_enterotypes.tsv", header = T)
e2 <- read.delim("data/MetaCardis_enterotypes.tsv", header = T)

## GALAXY/MicrobLiver
df1 <- data.frame(id = md1$ID, count = md1$count, enterotype = factor(e1$enterotype))
df1$enterotype <- df1$enterotype %>% str_replace("Bact_Phoc", "Bacteroides") %>% paste0(. , " type") %>% factor()

pair <- combn(levels(df1$enterotype), 2, simplify=F)
p1 <- ggplot(df1, aes(x = enterotype, y = log10(count), fill = enterotype)) +
  theme_classic() +
  ylab("Microbial load") +
  ggtitle("GALAXY/MicrobLiver") +  
  scale_y_continuous(limits = c(9.5, 11.9)) +
  stat_compare_means(comparisons = pair, label = "p.signif", tip.length = 0, label.y = c(11.3, 11.5, 11.7)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, face = "italic")) +
  theme(legend.position = "none", axis.title.x = element_blank()) +
  geom_jitter_rast(size = 0.05, width = 0.1, alpha = 0.1) +
  scale_fill_manual(values = col[c(5, 6, 3)]) +
  geom_boxplot(outlier.size = -1)
p1

## MetaCardis
df2 <- data.frame(id = md2$ID, count = md2$count, enterotype = factor(e2$enterotype))
df2$enterotype <- df2$enterotype %>% str_replace("Bact_Phoc", "Bacteroides") %>% paste0(. , " type") %>% factor()

pair <- combn(levels(df2$enterotype), 2, simplify=F)
p2 <- ggplot(df2, aes(x = enterotype, y = log10(count), fill = enterotype)) +
  theme_classic() +
  ylab("Microbial load") +
  ggtitle("MetaCardis") +  
  scale_y_continuous(limits = c(9.8, 12.1)) +
  stat_compare_means(comparisons = pair, label = "p.signif", tip.length = 0, label.y = c(11.6, 11.8, 12)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, face = "italic")) +
  theme(legend.position = "none", axis.title.x = element_blank()) +
  geom_jitter_rast(size = 0.05, width = 0.1, alpha = 0.1) +
  scale_fill_manual(values = col[c(5, 6, 3)]) +
  geom_boxplot(outlier.size = -1)
p2


## save plots
p <- ggarrange(p1, p2, ncol = 2)
ggsave(p, filename = "out/figure_1B.pdf", width = 4, height = 3)


