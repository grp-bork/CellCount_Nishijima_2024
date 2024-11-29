library(tidyverse)
library(ggthemes)
library(ggrastr)
library(egg)
library(scales)
library(ggpubr)

pal <- tableau_color_pal("Superfishel Stone")
ss <- pal(10)

d <- list()
d[[1]] <- read_rds("out/model/model_galaxy.motus25.rds.in_validation.rds")
d[[2]] <- read_rds("out/model/model_metacardis.motus25.rds.in_validation.rds")
d[[3]] <- read_rds("out/model/model_galaxy.motus25.rds.ex_validation.rds")
d[[4]] <- read_rds("out/model/model_metacardis.motus25.rds.ex_validation.rds")

d[[5]] <- read_rds("out/model/model_galaxy.ko.rds.in_validation.rds")
d[[6]] <- read_rds("out/model/model_metacardis.ko.rds.in_validation.rds")
d[[7]] <- read_rds("out/model/model_galaxy.ko.rds.ex_validation.rds")
d[[8]] <- read_rds("out/model/model_metacardis.ko.rds.ex_validation.rds")

d[[9]] <- read_rds("out/model/model_Vandeputte.rds.in_validation.rds")
d[[10]] <- read_rds("out/model/model_Vandeputte.rds.ex_validation.rds")

plot <- list()
i <- 1
for(i in 1:10){
  print(i)
  res <- d[[i]]
  
  pcc <- cor.test(res$count, res$pred, use = "pairwise.complete.obs", method = "pearson")
  pcc.r <- round(pcc$estimate, digits = 2)
  pcc.p <-  scientific(pcc$p.value, digits = 2) 
  
  lab1 <- sprintf(" italic(R): %s", pcc.r) %>% parse(text = .)
  lab2 <- sprintf(" italic(P): %s", pcc.p) %>% parse(text = .)
  
  xmin <- min(res$count)
  ymax <- max(res$pred)    
  
  if(i == 1 | i == 3){
    ## GALAXY
    col <- ss[4]
  }else{
    ## MetaCardis
    col <- ss[6] 
  }
  
  res %>% head()
  
  p <- ggplot(res, aes(x = count, y = pred)) +
    theme_mine() +
    xlab("Microbial load (log10)") +
    ylab("Predicted value") +
    
    theme(plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), "lines")) +
    annotate(geom = "text", x = xmin, y = ymax, label = lab1, hjust = 0, vjust = 0.1, size = 3) +
    annotate(geom = "text", x = xmin, y = ymax, label = lab2, hjust = 0, vjust = 1.6, size = 3) +
    geom_abline(linetype = 2, col = "gray") +
    theme(plot.title = element_text(size = 9)) +
    #geom_point_rast(alpha = 0.4, size = 0.6, col = col, stroke = 0) +
    geom_point_rast(alpha = 0.4, size = 0.6, col = "gray50", stroke = 0) +
    
    geom_smooth(method = "lm", se = T, col = t20[2]) 
  p
  
  plot[[i]] <- p
}

fig2 <- egg::ggarrange(plot[[1]], plot[[4]], plot[[5]], plot[[8]], plot[[9]], plot[[3]], plot[[2]], plot[[7]], plot[[6]], plot[[10]], nrow = 2)
ggsave(fig2, filename = "out/figures/figure_2.pdf", width = 12, height = 4)