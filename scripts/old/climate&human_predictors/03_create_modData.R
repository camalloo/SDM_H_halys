###################################################
############ Create the modelling data ############
###################################################

library(terra)
library(fuzzySim)
library(tidyverse)

folder <- 'climate_human' # name of the mother folder

######### Load data ----------
world <- geodata::world(path = 'data/geodata/', resolution = 1)
world <- world[!world$GID_0 == 'ATA'] # world map no Antarctica
# plot(world)

preds <- rast(file.path('output', folder, 'modellingData/predictors.tif'))
plot(preds[[1:ceiling(nlyr(preds)/2)]], legend = F)
plot(preds[[(ceiling(nlyr(preds)/2)+1):nlyr(preds)]], legend = F)

## load data from only_climate folder
mult <-readRDS(file.path('output', 'only_climate', 'modellingData/ModelRegion/region_mult.rds'))
occs <- read.csv(file.path('output', 'only_climate', 'modellingData/clean_occs.csv'))
occs <- cleanCoords(occs, coord.cols = c('long', 'lat'))
occs <- vect(occs, geom = c('long', 'lat'), crs = crs(world))
######### 


######### Modelling region ----------
## compute the buffer
b <- buffer(occs, width = mult)
## get the region
region <- aggregate(b)
region <- crop(region, world)
plot(region, col = 'gold', border = NA, main = 'Model training region')
plot(world, add = T)
points(occs, cex = 0.1)

## crop predictors according to model region
predsCropped <- crop(preds, region, mask = T, touches = F)
names(predsCropped)
# plot(predsCropped[[1:ceiling(nlyr(predsCropped)/2)]], legend = F)
# plot(predsCropped[[(ceiling(nlyr(predsCropped)/2)+1):nlyr(predsCropped)]], legend = F)
## plot one
plot(predsCropped[[1]])
plot(world, lw = 0.5, add = T)
## save to disk
writeRaster(predsCropped, filename = file.path('output', folder, 'modellingData/predsCropped.tif'),
            filetype = 'GTiff', gdal = c('COMPRESS=DEFLATE'), overwrite = T)
rm(preds, region)
gc()
######### 


######### Data gridding ----------
## presence VS. absence of record
gridData <- gridRecords(rst = predsCropped, pres.coords = occs)
nrow(gridData)
## clean from NAs in predictors if any
nas <- gridData |> 
  filter(if_any(everything(), is.na))
if (nrow(nas) > 0) {
  gridData <- drop_na(gridData)
  
  # nasSPAT <- vect(nas, geom = c('x', 'y'), crs = crs(world))
  # plet(predsCropped[['blt']]) |>
  #   points(occs, col = 'red') |>
  #   points(nasSPAT)
}

## save the raw gridded data to csv so I have them ready to go
write.csv(data.frame(gridData), file.path('output', folder, 'modellingData/GriddedData.csv'), row.names = F)
######### 


######## Create the modelling data --------
# gridData <- read.csv(file.path('output', folder, 'modellingData/GriddedData.csv'))
nrow(gridData)
table(gridData$presence)  # number of presences and absences

## thin the absences to be x the presences
multiplier = 1 # multiplication factor of presences to create pseudo-absences
data <- selectAbsences(gridData, 
                       sp.cols = 'presence', coord.cols = c('x', 'y'), CRS = crs(occs),
                       mult.p = multiplier, seed = 1312, df = T, bias = F)
plot(world, add = T, main = 'Modelling data')
table(data$presence)

## save test and train dfs to disk
write.csv(data.frame(data), file.path('output', folder, 'modellingData/modData.csv'), row.names = F)
######### 


######## Predictors analysis (VIF + pairwse correlation) on only added vars --------
data <- read.csv(file.path('output', folder, 'modellingData/modData.csv'))
humanLayers <- c('air', 'blt', 'prt', 'rds', 'rdt')
vars <- collinear::collinear(data, responses = 'presence', 
                             predictors = names(data)[names(data) %in% humanLayers],
                             max_cor = 0.6, max_vif = 5)
merged <- c(read.csv(file.path('output', 'only_climate', 'modellingData/selectedVariables.csv'))[, 2],
            vars$presence$selection)
write.csv(data.frame(preds = merged), file.path('output', folder, 'modellingData/selectedVariables.csv'))

names(predsCropped)
merged
######### 


######## Clean memory --------
rm(list = ls())
gc()

