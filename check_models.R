

model1 <- read_rds("out/model/model.GALAXY.mOTUs.rds")
model2 <- read_rds("out/model/model.MetaCardis.mOTUs.rds")

model3 <- read_rds("out/model/model.GALAXY.KO.rds")
model4 <- read_rds("out/model/model.MetaCardis.KO.rds")


model1$trainingData %>% dim()
model2$trainingData %>% dim()

model1$trainingData[1:5, 1:5]
model2$trainingData[1:5, 1:5]

model2$pred %>% dim()

getwd()
model <- read_rds("../../R/out/model/model_metacardis.motus25.rds")
model <- model[[1]]
model$trainingData %>% dim()
model$trainingData[1:5, 1:5]
model$pred %>% dim()
