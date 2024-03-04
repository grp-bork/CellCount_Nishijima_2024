library(tidyverse)
library(ggthemes)
pal <- tableau_color_pal("Superfishel Stone")
col <- pal(10)

## genus level analysis
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
  theme_mine() +
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
ggsave(p, filename = "out/supplementary_figure4.pdf", width = 8, height = 3, useDingbats=FALSE)

