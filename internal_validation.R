library(tidyverse)

## get path for each model
files <- list.files("out/model/", pattern = "rds") %>% paste0("out/model/", .)

## validate each model
for(i in files){
  print(i)
  data <- i %>% str_remove(".*model/model.") %>% str_remove(".rds")
  
  ## read model
  model <- read_rds(i)
  best <- model$bestTune
  pred <- model$pred %>% filter(nrounds == best$nrounds & max_depth == best$max_depth & eta == best$eta & gamma == best$gamma & colsample_bytree == best$colsample_bytree & min_child_weight == best$min_child_weight & subsample == best$subsample)
  
  ## prepare result
  count <- tapply(pred$obs, pred$rowIndex, mean)
  count.pred <- tapply(pred$pred, pred$rowIndex, mean)
  df <- data.frame(pred = count.pred, count = count, data = data, validation = "internal")
  
  ## save result
  out <- paste0("out/rds/internal_validation.", data, ".rds")  
  saveRDS(df, file = out)
}
