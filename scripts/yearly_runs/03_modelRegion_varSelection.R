###############################################################
############ Create the model region and apply VIF ############
###############################################################

library(terra)
library(fuzzySim)
library(tidyverse)

folder <- 'baseline_models_yearly_runs' # name of the mother folder

######### Load data ----------
world <- geodata::world(path = 'data/geodata/', resolution = 1)
world <- world[!world$GID_0 == 'ATA'] # world map no Antarctica

occs <- read.csv(file.path('output', folder, 'modellingData/clean_occs.csv'))
occs <- cleanCoords(occs, coord.cols = c('long', 'lat'), uncert.col = NULL, abs.col = NULL)
occs <- vect(occs, geom = c('long', 'lat'), crs = crs(world))

plot(world, lwd = 0, main = 'GBIF presence occurrences')
plot(occs, col = 'red', cex = 0.3, add = T)
plot(world, lwd = 0.5, add = T)

preds <- rast(file.path('output', folder, 'modellingData/predictors.tif'))
plot(preds[[1:ceiling(nlyr(preds)/2)]], legend = F)
plot(preds[[(ceiling(nlyr(preds)/2)+1):nlyr(preds)]], legend = F)
######### 


######### Modelling region ----------
# region <- getRegion(occs, type = 'width', width_mult = 0.05)

## instead of using FuzzySim::getRegion() I compute it manually so I can use the same buffer later for the yearly runs
## compute the width of the occurrences
w <- width(aggregate(occs))
scaling <- 0.05 # scaling factor for the width
wd <- w * scaling ## compute the diameter around each point
## save the mult object for later use
saveRDS(wd, file.path('output', folder, 'modellingData/ModelRegion/region_mult.rds'))
## compute the buffer
b <- buffer(occs, width = wd)
## get the region
region <- aggregate(b)
region <- crop(region, world)
plot(world, lwd = 0, main = 'Model training region', font.main = 1)
plot(region, col = 'gold', border = NA, add = T)
plot(world, add = T)
points(occs, cex = 0.1)
## save to disk
writeVector(region, file.path('output', folder, 'modellingData/ModelRegion/modellingRegionPoly.shp'), overwrite = T)

## crop predictors according to model region
predsCropped <- crop(preds, region, mask = T, touches = F)
# plot(predsCropped[[1:ceiling(nlyr(predsCropped)/2)]], legend = F)
# plot(predsCropped[[(ceiling(nlyr(predsCropped)/2)+1):nlyr(predsCropped)]], legend = F)
## plot one
plot(predsCropped[[1]])
plot(world, lw = 0.5, add = T)
## save to disk
writeRaster(predsCropped, filename = file.path('output', folder, 'modellingData/predsCropped.tif'),
            filetype = 'GTiff', datatype = 'FLT4S', gdal = c('COMPRESS=DEFLATE'), overwrite = T)
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
  nasSPAT <- vect(nas, geom = c('x', 'y'), crs = crs(world))
  gridData <- drop_na(gridData)
  
  plet(predsCropped[['dem']]) |>
    points(occs, col = 'red') |>
    points(nasSPAT)
}
table(gridData$presence)
## plot result
# plot(world)
# points(subset(gridData, gridData$presence == 0, select = c('x', 'y')), 
#        col = 'red', cex = 0.08)
# points(subset(gridData, gridData$presence == 1, select = c('x', 'y')), 
#        col = 'blue', cex = 0.08)
# legend("top", inset = c(0, 0.1), bg = "white", cex = 0.8, horiz = T, y.intersp = 0.01, xpd = T,
#        legend = c(paste('N total =', nrow(gridData)),
#                   paste('N pres =', sum(gridData$presence, na.rm = TRUE)),
#                   paste('N background =', nrow(gridData[gridData$presence == 0,]))), 
#        pch = c(NA, 21, 21), col = 'black', pt.bg = c(NA, 'blue', 'red'))
## interactive plot
# plet(vect(subset(gridData, presence == 1, select = c('x', 'y')), 
#           crs = crs(occs)), col = 'blue') |> 
#   points(vect(subset(gridData, presence == 0, select = c('x', 'y')), 
#               crs = crs(occs)), col = 'red')
## save the raw gridded data to csv so I have them ready to go
write.csv(data.frame(gridData), file.path('output', folder, 'modellingData/GriddedData.csv'), row.names = F)
######### 


######## Create the modelling data (select absences) --------
# gridData <- read.csv(file.path('output', folder, 'modellingData/GriddedData.csv'))
nrow(gridData)
table(gridData$presence)  # number of presences and absences
bias <- geodata::footprint(year = 2009, path = 'data/geodata/')
bias <- project(bias, preds)
bias <- crop(bias, region, mask = T)
plot(bias)
# ## first select x5 absences biased towards builtup areas (decrease GBIF sampling bias)
# biasedAbs <- selectAbsences(gridData, 
#                             sp.cols = "presence", coord.cols = c("x", "y"), CRS = crs(occs),
#                             mult.p = 5, seed = 1312, df = T, bias = bias)
# ## then remove these absences from the pool of total absences
# biasedAbs <- biasedAbs |> 
#   filter(presence == 0)
# gridData <- setdiff(gridData, biasedAbs)
# nrow(gridData)
# ## after that subset presence*(multiplier/2) randomly selected absences from the total absences - biased absences
# multiplier = 5 # multiplication factor of presences to create pseudo-absences
# data <- selectAbsences(gridData, 
#                        sp.cols = "presence", coord.cols = c("x", "y"), CRS = crs(occs),
#                        mult.p = multiplier, seed = 1312, df = T)
# ## finally, merge the biased absences to the random absences 
# data <- rbind(data, biasedAbs)
# message('The final presence:absence ratio is 1 + multilpier:\n1:', 5 + multiplier)
data <- selectAbsences(gridData,
                       sp.cols = "presence", coord.cols = c("x", "y"), CRS = crs(occs),
                       mult.p = 10, seed = 1312, df = T, bias = bias)
plot(world, add = T, main = 'Modelling data')
nrow(data)
table(data$presence)
## save test and train dfs to disk
write.csv(data.frame(data), file.path('output', folder, 'modellingData/modData.csv'), row.names = F)
######### 


######## Predictors analysis (VIF + pairwse correlation) --------
data <- read.csv(file.path('output', folder, 'modellingData/modData.csv'))
vars <- collinear::collinear(data, responses = 'presence', predictors = names(data)[-c(1:4)],
                             max_cor = 0.75, max_vif = 5)
write.csv(data.frame(preds = vars$presence$selection), file.path('output', folder, 'modellingData/selectedVariables.csv'))

names(predsCropped)
vars$presence$selection
######### 


######## Clean memory --------
rm(list = ls())
gc()

