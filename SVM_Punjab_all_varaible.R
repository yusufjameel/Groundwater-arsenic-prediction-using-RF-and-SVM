setwd("/Users/yusufjameel/Dropbox/Arsenic_MIT/Arsenic_prediction_indus_valley")

rm(list = ls(all=T))

library(rsample)      
library(randomForest) 
library(ranger)       
library(caret)        
library(h2o)         
library(tidyverse)
library(readxl)
library(data.table)
library(leaflet)
library(ggridges)
library(ggExtra)
library(corrplot)
library(raster)
library(tibble)
library(sf)
library(rnaturalearthdata)
library(tidyverse)
library(rnaturalearthhires)
library(kernlab)
library(e1071)


###############################

As_Indus = read.csv('Indus_Arsenic_kit_vars.csv', as.is =T)
unique(As_Indus$As__ppb)
As_Indus$As__ppb = floor(As_Indus$As__ppb)
unique(As_Indus$Recalibrated_As)
As_Indus$Recalibrated_As[As_Indus$As__ppb == 0] = 6
As_Indus$Recalibrated_As[As_Indus$As__ppb == 10] = 36
As_Indus$Recalibrated_As[As_Indus$As__ppb == 25] = 28
As_Indus$Recalibrated_As[As_Indus$As__ppb == 50] = 46
As_Indus$Recalibrated_As[As_Indus$As__ppb == 100] = 50
As_Indus$Recalibrated_As[As_Indus$As__ppb == 200] = 81
As_Indus$Recalibrated_As[As_Indus$As__ppb == 300] = 116
As_Indus$Recalibrated_As[As_Indus$As__ppb == 500] = 258
As_Indus$Recalibrated_As[As_Indus$As__ppb == 1000] = 719

str(As_Indus) ### check the class of each variable

###########convert factors  from numeric
As_Indus$sol_grtgroup_usda.soiltax_c_250m_s0..0cm_1950..2017_v0.2_punjab_250m_FACT = as.factor(As_Indus$sol_grtgroup_usda.soiltax_c_250m_s0..0cm_1950..2017_v0.2_punjab_250m_FACT)
As_Indus$sol_texture.class_usda.tt_m_250m_b0..0cm_1950..2017_v0.2_punjab_250m_FACT  = as.factor(As_Indus$sol_texture.class_usda.tt_m_250m_b0..0cm_1950..2017_v0.2_punjab_250m_FACT )
As_Indus$sol_texture.class_usda.tt_m_250m_b200..200cm_1950..2017_v0.2_punjab_250m_FACT = as.factor(As_Indus$sol_texture.class_usda.tt_m_250m_b200..200cm_1950..2017_v0.2_punjab_250m_FACT)
As_Indus$dtm_geology_usgs.ecotapestry_c_250m_s0..0cm_2014_v1.0_punjab_250m_FACT   = as.factor(As_Indus$dtm_geology_usgs.ecotapestry_c_250m_s0..0cm_2014_v1.0_punjab_250m_FACT)
As_Indus$dtm_landform_usgs.ecotapestry_c_250m_s0..0cm_2014_v1.0_punjab_250m_FACT  = as.factor(As_Indus$dtm_landform_usgs.ecotapestry_c_250m_s0..0cm_2014_v1.0_punjab_250m_FACT )
As_Indus$HYSOGs250m_punjab_250m_FACT    = as.factor(As_Indus$HYSOGs250m_punjab_250m_FACT)
As_Indus$MCD12Q1_IGBP2_punjab_250m_FACT   = as.factor(As_Indus$MCD12Q1_IGBP2_punjab_250m_FACT)

### reclassify arsenic categories
As_Indus$As_cat[As_Indus$Recalibrated_As < 10] = 0
As_Indus$As_cat[As_Indus$Recalibrated_As >= 10] = 1
As_Indus$As_cat = as.factor(As_Indus$As_cat)

As_Indus= As_Indus[complete.cases(As_Indus), ]

colnames(As_Indus)
all_data_select = As_Indus[,c(68, 6:8, 12:67)]

for (i in 1:56){
   colnames(all_data_select)[i+4] = paste0("var", i)
}

##################Import raster - they will be needed for predictions
rastlist <- list.files(path = "./Punjab", pattern='.tif$', all.files=TRUE, 
                       full.names=TRUE)
rast = lapply(rastlist, raster)

for (i in 1:56){
   names(rast[[i]]) = paste0("var",i)
}

all_stack <- raster::stack( rast)

dataframe_raster <- raster::as.data.frame(all_stack, centroids = T, xy = T)
dataframe_raster <- dataframe_raster %>% 
   drop_na()


colnames(dataframe_raster)[1] = "Long"
colnames(dataframe_raster)[2] = "Lat"

##########################################assign projection
wgs.84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #### projection

######################################reclassify rasters incorrectly identified as numeric into factos
dataframe_raster$Depth__m = 100
dataframe_raster = dataframe_raster[,c(59,2,1,3:58)]
dataframe_raster$var16 = as.factor(dataframe_raster$var16)
dataframe_raster$var18 = as.factor(dataframe_raster$var18)
dataframe_raster$var28 = as.factor(dataframe_raster$var28)
dataframe_raster$var29 = as.factor(dataframe_raster$var29)
dataframe_raster$var43 = as.factor(dataframe_raster$var43)
dataframe_raster$var44 = as.factor(dataframe_raster$var44)

######add all the levels from raster to the test and training data: 
##Training data have only few levels

all_data_select$var41 = as.numeric(all_data_select$var41)
levels(all_data_select$var16) <- levels(dataframe_raster$var16)
levels(all_data_select$var18) <- levels(dataframe_raster$var18)
levels(all_data_select$var28) <- levels(dataframe_raster$var28)
levels(all_data_select$var29) <- levels(dataframe_raster$var29)
levels(all_data_select$var43) <- levels(dataframe_raster$var43)
levels(all_data_select$var44) <- levels(dataframe_raster$var44)

######split the data in test and training
set.seed(123)
trainIndex <- createDataPartition(all_data_select$As_cat, p = .7, 
                                  list = FALSE, 
                                  times = 1)

iso_train <- all_data_select[trainIndex,]
iso_test <- all_data_select[-trainIndex,]

#iso_train = iso_train[sample(nrow(iso_train), 8000), ]

#########run the svm model - try linear, radial and polnomial kernal (using Caret package)
##polynomial kernel will take few hours to run
train_control <- trainControl(method="repeatedcv", number=5, repeats=5,
                              summaryFunction=twoClassSummary,sampling = "down",classProbs = TRUE)


svm1 <- train(As_cat~., data = iso_train, method = "svmLinear", trControl = train_control,  
              preProcess = c("center","scale"), tuneGrid = expand.grid(C = seq(0, 4, length = 40)))

summary(svm1) ### check output

svm2 <- train(As_cat~., data = iso_train, method = "svmRadial", trControl = train_control,  
              preProcess = c("center","scale"), tuneLength = 10, metric = "ROC")
summary(svm2) ### check output

svm3 <- train(As_cat~., data = iso_train, method = "svmPoly", trControl = train_control,  
             preProcess = c("center","scale"), tuneLength = 10, metric = "ROC")
summary(svm3) ### check output

####################run the radial kernel using e1071 package
svm_e1071 <- tune(svm, As_cat~., data = iso_train, kernel = "radial",  
  ranges = list(cost = seq(0, 100, length = 20), gamma = seq(0, 10, length = 20)))



### THIS IS THE BEST MODEL BASED UPON THE DIFFERENT TUNING PARAMETERS I TRIED
## we will use this model for our predictions
# kernel = radial
#gamma = 0.01
#cost = 10
svm_final = svm(As_cat~., data = iso_train, kernel = "radial",  
           cost =10, gamma = 0.01, probability =TRUE )
pred_svm <- predict(svm_final, iso_test,probability =TRUE )
head(pred_svm)


 ###############################get confusion matrix
model.df = data.frame(iso_test$As_cat, attr(pred_svm, "probabilities")[,2])
colnames(model.df) <- c("true", "predicted")
par(pty = "s")
pROC::roc(model.df$true, model.df$predicted, plot=TRUE,
          legacy.axes=TRUE, percent=TRUE, 
          xlab="False Positive Percentage", ylab="True Postive Percentage", 
          col="#377eb8", lwd=4, print.auc=TRUE)

##calcuate the optimale cutoff point
cp <- cutpointr(model.df, predicted, true, 
                method = maximize_metric, metric = sum_sens_spec, pos_class = 0)

plot(cp)

cutoff = cp$optimal_cutpoint
model.df$As_model_cat =  ifelse(model.df$predicted >= cutoff , 1, 0)  

model.df$As_model_cat = as.factor(model.df$As_model_cat)

confusionMatrix(model.df$As_model_cat, model.df$true)

  
  
  

##########################################predict Arsenic based on SVM 
rf_output <- predict(svm_final, dataframe_raster, type = "prob")

#########get the probability of arsenic >10 ppb for punjab

rf_model_Punjab <- tibble(longitude = dataframe_raster$Long,
                          latitude = dataframe_raster$Lat,
                          arsenic = rf_output[,2])

rf_model_raster <- raster::rasterFromXYZ(rf_model_Punjab)
proj4string(rf_model_raster) = wgs.84


## Plot the prediction map 
#######add some cites for references
long = c(74.30,71.50,73.2,74.8,73.1,76.8)
lat = c(31.45,30.16,31.6,31.6,33.5,30.7)
city = c(1,1,1,1,1,1)
cities = data.frame(long,lat,city)
cities_points <- cities %>% 
  st_as_sf(coords = c('long', 'lat'), crs = wgs.84)

##########################get state and country boundry
states_ind = rnaturalearth::ne_states(country = "india")
Punjab_ind = subset(states_ind, states_ind$name %in% c("Punjab"))
states_pak = rnaturalearth::ne_states(country = "pakistan")
Punjab_pak = subset(states_pak, states_pak$name %in% c("Punjab"))
ab <- aggregate(rbind(Punjab_pak,Punjab_ind))
ab ### only the state of punjab in India and Pakistan

#########################get rivers in the basin
river_punjab<- st_read(dsn ="./world_rivers_dSe/world_rivers_dSe.shp")
indus =   subset(river_punjab , river_punjab$Name1 == "Indus")
ravi =   subset(river_punjab , river_punjab$Name1 == "Ravi")
chenab =   subset(river_punjab, river_punjab$Name1 == "Chenab")
sutlej =   subset(river_punjab, river_punjab$Name1 == "Sutlej")
jhelum =   subset(river_punjab, river_punjab$Name1 == "Jhelum")


punjab_crop <- crop(rf_model_raster, ab)
punjab_crop <- mask(punjab_crop, ab)
st_crs(river_punjab)==st_crs(punjab_crop)


## use tmap for plotting
tm_shape(punjab_crop) + 
  tm_raster(n=8, palette="-RdBu", midpoint = 0.5) +
  tm_shape(cities_points) + tm_dots(size = 0.25, shape=20) +
  tm_shape( ab) + tm_polygons(col = "gray99",border.col = "black", lwd = 0.75,alpha = 0.1,legend.show = FALSE)+
  tm_layout(main.title  = "Random forest prediction for As > 10 ppb", legend.bg.color = "white", legend.bg.alpha=.5,
            legend.frame=F, legend.outside = F, bg.color = "gray80") + 
  tm_shape(indus) + tm_lines(col = "darkgreen",  alpha = 0.7, lwd = 1, lty = "solid") +
  tm_shape(ravi) + tm_lines(col = "darkgreen", alpha = 0.7, lwd = 1, lty = "solid") +
  tm_shape(chenab) + tm_lines(col = "darkgreen",  alpha = 0.7, lwd = 1, lty = "solid") +
  tm_shape(sutlej) + tm_lines(col = "darkgreen",  alpha = 0.7, lwd = 1, lty = "solid") +
  tm_shape(jhelum) + tm_lines(col = "darkgreen",  alpha = 0.7, lwd = 1, lty = "solid") +
  tm_grid(lwd = 0.25, labels.size = 1, n.x = 4, n.y =4) 

  