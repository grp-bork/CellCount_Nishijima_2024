library(tidyverse)
library(ggthemes)
library(ggrastr)
library(egg)
library(scales)
pal <- tableau_color_pal("Tableau 20")
col <- pal(20)

## external validation
external_validation_function <- function(model, d, count, title){
  table <- list()
  temp <- predict(model, d)
  
  res <- data.frame(count = count, pred = temp)
  scc <- cor.test(res$count, res$pred, use = "pairwise.complete.obs", method = "spearman")
  scc.v <- round(scc$estimate, digits = 2)
  scc.p <-  scientific(scc$p.value, digits = 2)  
  
  cor.test(res$count, res$pred, use = "pairwise.complete.obs", method = "pearson")
  pcc <- cor.test(res$count, res$pred, use = "pairwise.complete.obs", method = "pearson")
  
  pcc.r <- round(pcc$estimate, digits = 2)
  pcc.p <-  scientific(pcc$p.value, digits = 2)  
  
  lab1 <- sprintf(" italic(R): %s", pcc.r) %>% parse(text = .)
  lab2 <- sprintf(" italic(P): %s", pcc.p) %>% parse(text = .)
  lab1 %>% mode()
  
  xmin <- min(res$count)
  ymax <- max(res$pred)    
  
  xmin
  cor_pred <- ggplot(res, aes(x = count, y = pred)) +
    theme_classic() +
    ggtitle(paste0(title, " (ext. validation)")) +
    xlab("Microbial load") +
    ylab("Predicted value") +    
    annotate(geom = "text", x = xmin, y = ymax, label = lab1, hjust = 0, vjust = 0.1) +
    annotate(geom = "text", x = xmin, y = ymax, label = lab2, hjust = 0, vjust = 1.6) +
    geom_abline(linetype = 2, col = "gray") +
    #geom_point_rast(alpha = 0.4, size = 0.1) +
    geom_point_rast(size = 0.1) +
    geom_smooth(method = "lm", se = F, col = col[1]) 
  cor_pred
  
  table[[1]] <- cor_pred
  table[[2]] <- res
  return(table)
}


filter_minor_species_and_add_shannon <- function(d, c, thre = 1E-4, pseud = 1E-4){
  shannon <- vegan::diversity(d)
  d$`Shannon diversity` <- shannon
  
  keep <- apply(d, 2, mean) > thre & apply(d > 0, 2, mean) > 0.1
  d <- d[, keep]
  return(d)
}


## prepare same df for prediction
prepare_same_label <- function(d.tr, d.ev){
  d.tr <- d.tr %>% select(-contains("outcome"))
  d.tr %>% dim() %>% print()
  d.ev %>% dim() %>% print()  
  
  keep <- colnames(d.ev) %in% colnames(d.tr)
  d.ev <- d.ev[, keep]
  
  keep <- !colnames(d.tr) %in% colnames(d.ev)
  table(keep)
  
  d.ev <- d.ev %>% data.frame(check.names = F)
  #d.ev[, colnames(d.tr)[keep]]
  
  if(sum(keep) > 0){
    #missing <- colnames(d.tr)[keep][1]
    d.ev[, colnames(d.tr)[keep]] <- 0
    keep <- colnames(d.ev) %in% colnames(d.tr)
    d.ev <- d.ev[, keep]
  }
  dim(d.ev)
  
  d.ev <- d.ev[, order(match(colnames(d.ev), colnames(d.tr)))] 
  dim(d.ev)
  return(d.ev)
}

files <- read.delim("data/validation_files.txt", header = T, check.names = F)
files
files$model <- files$model %>% paste0("out/model/", .)
files$training_data <- files$training_data %>% paste0("data/", .)
files$validation_data <- files$validation_data %>% paste0("data/", .)
files

i <- 3
res <- c()
for(i in 1:nrow(files)){
    
  ## read model
  model <- read_rds(files$model[i])
  
  ## read data for validation  
  d.validation <- files$validation_data[i]
  c.validation <- files$validation_data[i] %>% 
    str_replace("mOTU.*\\.", "load.") %>% 
    str_replace("KO.*\\.", "load.") %>%     
    str_replace("_16S", "_load")
  c.validation
  
  d.ev <- fread(d.validation) %>% data.frame(check.names = F) %>% column_to_rownames("V1")
  c.ev <- read.delim(c.validation, row.names = 1, header = T, check.names = F)
  
  if(!grepl("Vande", files$model[i])){
    d.ev <- filter_minor_species_and_add_shannon(d.ev, c.ev) ## filter some samples and add Shannon diversity
    c.ev <- c.ev$count %>% log10()
  }else{
    d.ev <- d.ev
    c.ev <- c.ev$Cell_count %>% log10()
  }
  
  ## prepare same labels with the same order with that of training data
  d.tr <- model$trainingData
  d.tr <- d.tr[, -ncol(d.tr)]
  d.ev <- prepare_same_label(d.tr, d.ev)
  
  # sanity check
  identical(colnames(d.ev), colnames(d.tr)) %>% print()
  
  ## scaling
  d.min <- min(d.ev[d.ev != 0])/2
  d.ev <- log10(d.ev + d.min) %>% scale()
  
  ## evaluate the model
  pred <- predict(model, d.ev)
  df <- data.frame(pred = pred, count = c.ev, data = data, validation = "external")
  
  ## save result
  out <- paste0("out/rds/internal_validation.", data, ".rds")  
  saveRDS(df, file = out)
}

