library(tidyverse)
library(ggthemes)
pal <- tableau_color_pal("Superfishel Stone")
col <- pal(10)

d1 <- read.delim("data/GALAXY_load.tsv")
d2 <- read.delim("data/MetaCardis_load.tsv")

## compare histogram between GALAXY and MetaCardis
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
ggsave(p_hist, filename = "out/supplementary_figure2.pdf", width = 4.5, height = 2.5)


