# R script to construct prediction models for fecal microbial load

# load library
library(tidyverse)
library(data.table)
library(ggpubr)
library(vegan)
library(doMC)
library(MLmetrics)
library(caret)
library(xgboost)
registerDoMC(cores = 4)

## define function
construct_model <- function(d, count, cutoff = 1E-4, pseud = 1E-4, data){
  set.seed(0)
  
  ## add Shannon diversity
  d$`Shannon diversity` <- diversity(d)
  
  ## filtering minor species
  keep <- apply(d, 2, mean) > cutoff & apply(d > 0, 2, mean) > 0.1
  d <- d[, keep]
  
  ## scaling
  d <- scale(log10(d + pseud))
  
  #train_control <- trainControl(method = "repeatedcv", savePredictions = T, repeats = 5)
  train_control <- trainControl(method = "repeatedcv", savePredictions = T, repeats = 1)
  tune_grid <- expand.grid(
    nrounds=c(100, 500, 1000),
    max_depth = c(3:5),
    eta = c(0.01, 0.05, 0.1),
    gamma = c(0.01),
    colsample_bytree = c(0.75),
    subsample = 1,
    min_child_weight = 1
  )
  
  model <- caret::train(d, count, trControl = train_control, method = "xgbTree", metric = "Rsquared", tuneGrid = tune_grid)
  best <- model$bestTune
  pred <- model$pred %>% filter(nrounds == best$nrounds & max_depth == best$max_depth & eta == best$eta & gamma == best$gamma & colsample_bytree == best$colsample_bytree & min_child_weight == best$min_child_weight & subsample == best$subsample)
  
  count.pred <- tapply(pred$pred, pred$rowIndex, mean)
  df <- data.frame(pred = count.pred, count = count, data = data, validation = "internal")
  
  ## save result
  out <- paste0("out/model/model.", data, ".rds")
  saveRDS(model, file = out)
  
  #out <- paste0("out/internal_validation.", data, ".pdf")
  #ggsave(p, filename = out, width = 3, height = 3)
  out <- paste0("out/rds/internal_validation.", data, ".rds")  
  saveRDS(df, file = out)
}

## read data
## GALAXY/MicrobLiver
d <- fread("data/GALAXY_mOTUs_v25.tsv") %>% data.frame(check.names = F) %>% column_to_rownames("V1")
c <- read.delim("data/GALAXY_load.tsv", header = T, check.names = F)
count <- c$count %>% log10()
construct_model(d, count, 1E-4, 1E-4, "GALAXY.mOTUs")

## MetaCardis
d <- fread("data/MetaCardis_mOTUs_v25.tsv") %>% data.frame(check.names = F) %>% column_to_rownames("V1")
c <- read.delim("data/MetaCardis_load.tsv", header = T, check.names = F)
count <- c$count %>% log10()
construct_model(d, count, 1E-4, 1E-4, "MetaCardis.mOTUs")



## GALAXY/MicrobLiver (KO)
d <- fread("data/GALAXY_KO.tsv") %>% data.frame(check.names = F) %>% column_to_rownames("V1")
c <- read.delim("data/GALAXY_load.tsv", header = T, check.names = F)
count <- c$count %>% log10()
construct_model(d, count, 1E-6, 1E-7, "GALAXY.KO")


## MetaCardis (KO)
d <- fread("data/MetaCardis_KO.tsv") %>% data.frame(check.names = F) %>% column_to_rownames("V1")
c <- read.delim("data/MetaCardis_load.tsv", header = T, check.names = F)
count <- c$count %>% log10()
construct_model(d, count, 1E-6, 1E-7, "MetaCardis.KO")
