library("terra")
library("fuzzySim")


################# Load Data_________________________________
folder <- 'AllPredictors' # name of the mother folder
world <- geodata::world(path = 'data/geodata/', resolution = 1)
world <- world[!world$GID_0 == 'ATA'] # world map no Antarctica
plot(world)
bio <- rast('data/CHELSAbioclim/historicCropped/bioAll.tif')





################# Predictors porcessing_________________________________
## plot them
plot(bio[[1:11]])
plot(bio[[12:22]])
plot(bio[[23:33]])





################# Process occurrences_________________________________
world <- world[!world$GID_0 %in% c('ATA', 'CHN', 'TWN', 'JPN', 'PRK', 'KOR', 'MAC', 'HKG')] # remove Antarctica & native countries (China, Japan, Taiwan, Koreas, Macao, Hong Kong)
# europe <- dataClean |> 
#   filter(continent == "EUROPE" | countryCode == 'TR') |> 
#   CoordinateCleaner::cc_dupl(lon = "decimalLongitude", lat = "decimalLatitude",value = "clean")

## spatialize the data and assign crs
occs <- vect(dataClean, geom = c("decimalLongitude", "decimalLatitude"), crs = 'epsg:4326', keepgeom = T)
plot(world, main = 'Occurrences used in training')
points(occs, col = 'red', cex = 0.3)

## additional check on occurrences
rmvOccs <- occs[extract(world, occs)$GID_0 %in% c('ATA', 'CHN', 'TWN', 'JPN', 'PRK', 'KOR', 'MAC', 'HKG')]
if (nrow(rmvOccs) == 0) {
  ## if 0, remove the object
  print('no additional occs removed')
  rm(rmvOccs)
} else {
  print(paste('removing'), nrow(rmvOccs), 'additional occs')
  occs <- occs - rmvOccs
  plot(world)
  points(occs, col = 'red', cex = 0.3)
  rm(rmvOccs)
}
# plot(occs, col = 'red', cex = 0.3)
# plot(world, add = T)





################# PCA_________________________________
######## process predictors --------
## aggregate the variables to ~10 km
res(bio)
bio <- aggregate(bio, fact = 10, fun = 'mean', na.rm = T)
res(bio)
## clean the names of the bioclim variables
names(bio) <- gsub('CHELSA_|_1981-2010_V.2.1', '', names(bio))
plot(bio[[1:11]])
plot(bio[[12:22]])
plot(bio[[23:33]])

## rasters contain NAs within the world maps
## I need to transform them to 0 otherwise they mess up the pca
bio <- subst(bio, NA, 0)
plot(bio[[1:11]])
plot(bio[[12:22]])
plot(bio[[23:33]])
## crop again to keep only land
#bio <- crop(bio, world, mask = T)


######## pca on the cropped rasters --------
pca <- prcomp(bio, scale. = T)
summary(pca)
plot(pca)
biplot(pca, xlabs = rep("", nrow(pca$x)))
rotations <- as.data.frame(pca$rotation)
rotationsABS <- abs(rotations)
## transform back to a raster stack with values the PC scores of the PCs that explain > 95% of variance [PC8, so a stack of 8 layers, one for each PC included]
PCscores <- predict(bio, pca, index = 1:8)
## keep only land
PCscores <- crop(PCscores, world, mask = T, snap = 'out', touches = T)
## plot it
plot(PCscores)
## plot one
plot(PCscores[[1]])
plot(world, lw = 0.5, add = T)


plot(region, border = 'magenta', col = NA, add = T)


## save to disk
writeRaster(PCscores, paste0('output/', folder, '/modellingData/ModelRegion/PCscoresWorld.tif'), filetype = 'GTiff', overwrite = T)






################# Modelling region_________________________________
region <- getRegion(occs, type = 'width', width_mult = 0.05)
plot(world, add = T)
## keep only land
region <- crop(region, world)
regionEXT <- crop(world, ext(region))
plot(region, col = 'yellow', main = 'Model training region')
plot(regionEXT, add = T)
points(occs, col = 'red', cex = 0.3)

## save to disk
writeVector(region, paste0('output/', folder, '/modellingData/ModelRegion/modellingRegionPoly.shp'), overwrite = T)
writeVector(regionEXT, paste0('output/', folder, '/modellingData/ModelRegion/modellingRegionEXT.shp'), overwrite = T)

## crop predictors according to model region
bioCropped <- crop(PCscores, region, mask = T, snap = 'out', touches = T)
plot(bioCropped[[1:11]])
plot(bioCropped[[12:22]])
plot(bioCropped[[23:33]])
## plot one
plot(bioCropped[[1]])
plot(world, lw = 0.5, add = T)
plot(region, border = 'magenta', col = NA, add = T)
## save to disk
writeRaster(bioCropped, paste0('output/', folder, '/modellingData/ModelRegion/predsCropped.tif'), filetype = 'GTiff', overwrite = T)





################# Data gridding_________________________________
world <- geodata::world(path = 'data/geodata/', resolution = 1)[!world$GID_0 == 'ATA']
gridData <- gridRecords(rst = bioCropped, pres.coords = occs)
plot(world)
points(subset(gridData, gridData$presence == 0, select = c('x', 'y')), col = 'red', cex = 0.08)
points(subset(gridData, gridData$presence == 1, select = c('x', 'y')), col = 'blue', cex = 0.08)
legend("top", inset = c(0, 0), bg = "white", cex = 0.8, horiz = TRUE, y.intersp = 0.01, xpd = T,
       legend = c(paste('N total =', nrow(gridData)),
                  paste('N pres =', sum(gridData$presence, na.rm = TRUE)),
                  paste('N background =', nrow(gridData[gridData$presence == 0,]))), 
       pch = c(NA, 21, 21), col = 'black', pt.bg = c(NA, 'blue', 'red'))

## save the modelling data to csv so I have them ready to go
write.csv(data.frame(gridData), paste0('output/', folder, '/modellingData/GridDataAll.csv'), row.names = F)


######## Create the modelling data --------
nrow(gridData)
table(gridData$presence)  # number of presences and absences

## thin the absences to be x the presences and add a bias layer to prioritize choice of absence close to human infrastructures
#---------- NOT SURE ABOUT USING THE BIAS LAYER if used change fuzzySim::selectAbsences(..., bias = biasLayer)
# biasLayer <- geodata::footprint(year = 2009, path = "data/")
# biasLayer <- crop(biasLayer, regionEXT)
# plot(biasLayer)
multiplier = 1 # multiplication factor of presences to create pseudo-absences
data <- fuzzySim::selectAbsences(gridData, 
                                 sp.cols = "presence", coord.cols = c("x", "y"), 
                                 mult.p = multiplier, seed = 1312, df = T, bias = F)
plot(world, add = T, main = 'Modelling data')
table(data$presence)

## plot training data to see the differences from all gridded data
par(mfrow = c(2, 1))
plot(world, main = paste('unfiltered data:', nrow(gridData)))
points(subset(gridData, gridData$presence == 0, select = c('x', 'y')), col = 'red', cex = 0.08)
points(subset(gridData, gridData$presence == 1, select = c('x', 'y')), col = 'blue', cex = 0.08)
plot(world, main = paste('thinned data:', nrow(data)))
points(subset(data, data$presence == 0, select = c('x', 'y')), col = 'pink3', cex = 0.08)
points(subset(data, data$presence == 1, select = c('x', 'y')), col = 'blue', cex = 0.08)
# points(subset(testData, testData$presence == 1, select = c('x', 'y')), col = 'green', cex = 0.03) # add test points
par(mfrow = c(1, 1))
legend('right', xpd = T, cex = .8,
       legend = c(paste('Presence:                 ', sum(data$presence)), 
                  paste('Background all:        ', nrow(gridData[gridData$presence == 0, ])),
                  paste('Background thinned:', sum(data$presence == 0, na.rm = T))),  
       pch = 21, col = 'black', pt.bg = c('blue', 'red', 'pink3'))

## save test and train dfs to disk
write.csv(data.frame(data), paste0('output/', folder, '/modellingData/modData.csv'), row.names = F)





# PCA PERFORMED SO VIF NOT USEFUL BECAUSE PCSCORES ARE UNCORRELATED
################# Variables selection_________________________________
# vars <- collinear::collinear(data, response = 'presence', predictors = names(data)[-c(1:4)], 
#                              max_cor = 0.7, max_vif = 10)
# write.csv(as.data.frame(vars), paste0('output/', folder, '/modellingData/selectedVariables.csv'), row.names = F)
vars <- names(PCscores)
write.csv(as.data.frame(vars), paste0('output/', folder, '/modellingData/selectedVariables.csv'), row.names = F)





################# Clean memory_________________________________
rm(list = ls())
gc()
