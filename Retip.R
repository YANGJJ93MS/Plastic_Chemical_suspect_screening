reticulate::install_python(version = "3.10")
install.packages('rJava', repos='http://cran.rstudio.com/')
install.packages('devtools', version='2.4.5', repos='http://cran.rstudio.com/')
install.packages('caret', version='6.0-94', repos='http://cran.rstudio.com/')
install.packages('ggplot2', version='3.5.1', repos='http://cran.rstudio.com/')
install.packages('rcdk', version='3.8.1', repos='http://cran.rstudio.com/')
install.packages('rcdklibs', version='2.9', repos='http://cran.rstudio.com/')
install.packages('doParallel', version='1.0.17', repos='http://cran.rstudio.com/')
install.packages('stringi', version='1.8.4', repos='http://cran.rstudio.com/')
install.packages('lattice', version='0.22-5', repos='http://cran.rstudio.com/')
install.packages('randomForest', version='4.7-1.1', repos='http://cran.rstudio.com/')
install.packages('xgboost', version='1.7.7.1', repos='http://cran.rstudio.com/')
install.packages('brnn', version='0.9.3', repos='http://cran.rstudio.com/')
install.packages('lightgbm', version='4.3.0', repos='http://cran.rstudio.com/')
install.packages('h2o', version='3.44.0.3')
install.packages('gtable', version='0.3.5', repos='http://cran.rstudio.com/')
install.packages('grid', version='4.4.0', repos='http://cran.rstudio.com/')
install.packages('gridExtra', version='2.3', repos='http://cran.rstudio.com/')
install.packages('reticulate', version='1.37', repos='http://cran.rstudio.com/')
devtools::install_github('olobion/Retiplib')
devtools::install_github('olobion/Retip')

library(Retip)
library(readxl)
library(reticulate)

setwd(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))
keras_installed <- TRUE
if (!keras_installed) {
  keras::install_keras()
}

prep.wizard()
# setwd("D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/RT_prediction")
rp2 <- readxl::read_excel("D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/RT_prediction/RT_library_DTSC_SMILES_fromCAS_retip.xlsx",
                          col_types = c("text", "text", "text", "numeric"))
descs <-getCD(rp2)
#clean dataset from NA and low variance value
db_rt <- proc.data(descs)
# chem.space(db_rt, target = "T3DB")

#> Split in training and testing using caret::createDataPartition
set.seed(101)
in_training <- caret::createDataPartition(db_rt$XLogP, p = .8, list = FALSE)
training <- db_rt[in_training, ]
testing <- db_rt[-in_training, ]

build_models <- TRUE

if (build_models) {
  rf  <- fit.rf(training)
  saveRDS(rf, "rf_model.rds")
  
  brnn <- fit.brnn(training)
  saveRDS(brnn, "brnn_model.rds")
  
  keras <- fit.keras(training, testing)
  save_model_hdf5(keras, filepath = "keras_model.h5")
  
  lightgbm <- fit.lightgbm(training, testing)
  saveRDS(lightgbm, "lightgbm_model.rds")
  
  xgb <- fit.xgboost(training)
  saveRDS(xgb, "xgb_model.rds")
  
  aml <- fit.automl.h2o(training)
  # It saves by default the best model
  h2o::h2o.saveModel(aml@leader, "automl_h2o_model")
} else {
  xgb <- readRDS("D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/RT_prediction/xgb_model.rds")
  rf <- readRDS("D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/RT_prediction/rf_model.rds")
  brnn <- readRDS("D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/RT_prediction/brnn_model.rds")
  keras <- keras::load_model_hdf5("D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/RT_prediction/keras_model.h5")
  lightgbm <- readRDS("D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/RT_prediction/lightgbm_model.rds")
  h2o::h2o.init(nthreads = -1, strict_version_check=FALSE)
  aml <- h2o.loadModel("automl_h2o_model/###") # replace ### with the name of your saved model
}

#> first you have to put the testing dataframe and then the name of the models
#> you have computed
stat <- get.score(testing, xgb, rf, brnn, keras, lightgbm, aml)
stat

#> Last value is the title of the graphics.
p.model(testing, m = xgb, title = "XGBoost - PlasticChemicals")

#> Last value is the title of the graphics.
p.model(testing, m = rf, title = "Random forest - PlasticChemicals")


#> Last value is the title of the graphics.
p.model(testing, m = brnn, title = "BRNN - PlasticChemicals")

#> Last value is the title of the graphics.
p.model(testing, m = lightgbm, title = "LightGBM - PlasticChemicals")

#> Last value is the title of the graphics.
p.model(testing, m = keras, title = "Keras - PlasticChemicals")

#> Last value is the title of the graphics.
p.model(testing, m = aml, title = "H2O autoML - PlasticChemicals")

p.model.features(xgb, mdl_type = "xgb")

##xgb generate best MAE, select XGB
# install.packages("webchem")
library(data.table)
library(webchem)
library(rcdk)
xgb <- readRDS("D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/RT_prediction/xgb_model.rds")
plasticpath <- "D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/RT_prediction/PlasticChem_parent.xlsx"
rp_ext <- readxl::read_excel(plasticpath,col_types = c('text','text','text'))
rp_ext <- as.data.table(rp_ext)

# rp_ext$validsmi <- unlist(lapply(rp_ext$SMILES,is.smiles))
# rp_ext_invalid <- rp_ext[which(rp_ext$validsmi=='FALSE')]

rp_ext1 <- rp_ext[1:2000]
rp_ext2 <- rp_ext[2001:2325]
rp_ext3 <- rp_ext[2326:2327]
rp_ext4 <- rp_ext[2329:dim(rp_ext)[1]]
View(rp_ext3)

rp_ext_desc1 <- getCD(rp_ext1)
rp_ext_desc2 <- getCD(rp_ext2)
rp_ext_desc3 <- getCD(rp_ext3)
rp_ext_desc4 <- getCD(rp_ext4)

rp_ext_pred_xgb2 <- Retip::RT.spell(training, rp_ext_desc2, model = xgb)
rp_ext_pred_xgb3 <- Retip::RT.spell(training, rp_ext_desc3, model = xgb)
rp_ext_pred_xgb4 <- Retip::RT.spell(training, rp_ex_desc4, model = xgb)
# write.csv(rp_ext_pred_xgb,file = "D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/RT_prediction/PlasticChem_parent_retip_2000.csv")
write.csv(rp_ext_pred_xgb2,file = "D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/RT_prediction/PlasticChem_parent_retip_2001-2325.csv")
write.csv(rp_ext_pred_xgb3,file = "D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/RT_prediction/PlasticChem_parent_retip_2326-2327.csv")
write.csv(rp_ext_pred_xgb4,file = "D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/RT_prediction/PlasticChem_parent_retip_2329-5144.csv")

rp_ext_pred_xgb1 <- fread("D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/RT_prediction/PlasticChem_parent_retip_2000.csv")
rp_ext_pred_xgb2 <- fread("D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/RT_prediction/PlasticChem_parent_retip_2001-2325.csv")
rp_ext_pred_xgb3 <- fread("D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/RT_prediction/PlasticChem_parent_retip_2326-2327.csv")
rp_ext_pred_xgb4 <- fread("D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/RT_prediction/PlasticChem_parent_retip_2329-5144.csv")

newpred = rbind(rp_ext_pred_xgb1,rp_ext_pred_xgb2,
                rp_ext_pred_xgb3,rp_ext_pred_xgb4)
#left join 
rp_ext[newpred,on='NAME', predRT := RTP]
rp_ext[,predRT_lower:=predRT-3.86][,predRT_higher:=predRT+3.86]

path="D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/RT_prediction/"
filelst <- list.files(path="D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/RT_prediction/",
           pattern = ".csv")
#output predicted RT
origpredRT <- fread(paste0(path,filelst[2]))
# Define the prefixes to exclude
prefixes <- c('ring', 'bond', 'chain', 'group', 'atom')
# Create a pattern that matches any of the prefixes at the start of the column names
pattern <- paste0("^", prefixes, collapse = "|")
# Deselect columns starting with the specified prefixes using data.table syntax
origpredRT_subset <- origpredRT[, !grepl(pattern, names(origpredRT)), with = FALSE]
origpredRT_subset[['newPred_rtip']] = rp_ext[['predRT']]
origpredRT_subset[, minevsrtip:= abs(newPred_rtip - `Predicted RT`)]
origpredRT_subset[['predRT_higher']] = rp_ext[['predRT_higher']]
origpredRT_subset[['predRT_lower']] = rp_ext[['predRT_lower']]
write.csv(origpredRT_subset, file =paste0(path,'PlasticChem_predRT.csv'))

#plot chemspace of my dataset and the new dataset
  ##filter not calculable smiles
  rp_ext_desc1 <- getCD(rp_ext1)
  rp_ext_desc2 <- getCD(rp_ext2)
  rp_ext_desc3 <- getCD(rp_ext3)
  rp_ext_desc4 <- getCD(rp_ext4)
  plastic <- rbind(rp_ext_desc1,
                   rp_ext_desc2,
                   rp_ext_desc3,
                   rp_ext_desc4)
  plastic <- proc.data(plastic)
  
  #
  db_rt_new <- db_rt[,2:dim(db_rt)[2]]
  chem.space(db_rt_new, plastic, title ='RTlibrary_PlasticChem')
  
