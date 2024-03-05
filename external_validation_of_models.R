library(tidyverse)
library(ggthemes)
library(ggrastr)
library(data.table)
library(egg)
library(scales)
library(vegan)

pal <- tableau_color_pal("Tableau 20")
col <- pal(20)


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

## read files
files <- read.delim("data/validation_files.txt", header = T, check.names = F)
files$model <- files$model %>% paste0("out/model/", .)
files$training_data <- files$training_data %>% paste0("data/", .)
files$validation_data <- files$validation_data %>% paste0("data/", .)

## validate each file
for(i in 1:nrow(files)){
  print(i)

  data <- files$training_data[i] %>% str_remove("data/") %>% str_remove(".tsv") %>% 
  str_remove("_v25") %>% 
  str_remove(".2021.*") %>%     
  str_replace("_", ".")
  
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
  
  d.ev$`Shannon diversity` <- diversity(d.ev)
  c.ev <- c.ev$count %>% log10()
  
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
  
  ## save result
  df <- data.frame(id = rownames(d.ev), pred = pred, count = c.ev, data = data, validation = "external")
  df %>% head()
  
  out <- paste0("out/rds/external_validation.", data, ".rds")  
  saveRDS(df, file = out)
}

