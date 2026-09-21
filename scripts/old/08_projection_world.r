library("terra")
library("maxnet")
library("mgcv")
library("embarcadero")
library("fuzzySim")


################# Load data_________________________________
folder <- 'world' # name of the mother folder

data <- read.csv(paste0('output/', folder, '/modellingData/modData.csv')) # train dataset
vars <- read.csv(paste0('output/', folder, '/modellingData/selectedVariables.csv'))[ , 1] # variables after collinearity analysis
bio <- rast('data/CHELSAbioclim/historicCropped/bioAll.tif')
world <- geodata::world(path = 'data/geodata/', resolution = 1)
world <- world[!world$GID_0 == 'ATA'] # world map no Antarctica
region <- vect(paste0('output/', folder, '/modellingData/ModelRegion/modellingRegionPoly.shp')) # modelling region polygon

## check which models to load based on AUC
AUCs <- read.csv(paste0('output/', folder, '/models/evaluation/AUCall.csv'))
AUCs

## load the models (RF overfitting so not used)
mxtMod <- readRDS(paste0('output/', folder, '/models/mxt_mod.rds'))
gamMod <- readRDS(paste0('output/', folder, '/models/gam_mod.rds'))
## BART MODEL cannot be saved so always re-fit it if not loaded on env (KEEP THE SEED)
bartMod <- dbarts::bart(x.train = data[, vars], y.train = data[, 'presence'],
                        ntree = 100, keeptrees = T, verbose = F, seed = 1312)




################# Project models in space_________________________________
######## check bioclim vars --------
plot(bio[vars])
res(bio)


######## run projection --------
## maxent
mxtProj <- predict(bio, mxtMod, type = 'cloglog', na.rm = T, cores = 4, cpkgs = 'maxnet')
names(mxtProj) <- 'MXT.pred'
## plot
plot(mxtProj, main = 'MXT', range = c(0, 1))

## gam
gamProj <- predict(bio, gamMod, type = 'response')
names(gamProj) <- 'GAM.pred'
gamProjFav <- Fav(pred = gamProj, sample.preval = prevalence(model = gamMod))
names(gamProjFav) <- 'GAM.fav'
## plot
par(mfrow = c(2, 1))
plot(gamProj, main = 'GAM', range = c(0, 1))
plot(gamProjFav, main = 'GAM Fav', range = c(0, 1))

## bart
bartProj <- predict2.bart(bartMod, raster::stack(bio), quantiles = c(.05, .95), splitby = 5, quiet = T)
names(bartProj) <- c('BART.pred', 'BART.lower.int', 'BART.upper.int')
bartProj <- rast(bartProj)
crs(bartProj) <- 'epsg:4326'
## convert to prevalence-independent favourability
bartProjFav <- Fav(pred = bartProj, sample.preval = prevalence(model = bartMod))
names(bartProjFav) <- c('BART.fav', 'BART.fav.lower.int', 'BART.fav.upper.int')
## add the uncertainty range
bartProj$uncertainty.range <- bartProj[[3]] - bartProj[[2]] # to pred
bartProjFav$uncertainty.range <- bartProjFav[[3]] - bartProjFav[[2]] # to fav
## plot
par(mfrow = c(1, 1))
plot(bartProj, range = c(0, 1))
title(main = 'BART', outer = T, line = -3)
plot(bartProjFav, range = c(0, 1))
title(main = 'BART Fav', outer = T, line = -3)

## save all to disk
writeRaster(mxtProj, filename = paste0('output/', folder, '/projections/rast/mxt_proj.tif'), gdal = c('COMPRESS=DEFLATE'), overwrite = TRUE)
writeRaster(gamProj, filename = paste0('output/', folder, '/projections/rast/gam_proj.tif'), gdal = c('COMPRESS=DEFLATE'), overwrite = TRUE)
writeRaster(gamProjFav, filename = paste0('output/', folder, '/projections/rast/gam_projFav.tif'), gdal = c('COMPRESS=DEFLATE'), overwrite = TRUE)
writeRaster(bartProj, filename = paste0('output/', folder, '/projections/rast/bart_proj.tif'), gdal = c('COMPRESS=DEFLATE'), overwrite = TRUE)
writeRaster(bartProjFav, filename = paste0('output/', folder, '/projections/rast/bart_projFav.tif'), gdal = c('COMPRESS=DEFLATE'), overwrite = TRUE)


######## combine projections in a single raster --------
## check crs
crs(mxtProj, describe = T)
crs(gamProjFav, describe = T)
crs(bartProjFav, describe = T)

## combine in a single object
## I USE FAVOURABILITY FOR BART AND GAM SO IT IS MORE COMPARABLE TO MAXENT
allProj <- c(mxtProj, gamProjFav, bartProjFav[[1]])
allProj
names(allProj) <- sub('.pred|.fav', '', names(allProj))
plot(allProj, range = c(0, 1))
title(main = 'Model projections', outer = T, line = -3)

################# Compute ensemble_________________________________
## prepare the AUCs df
mods <- names(allProj)
modsAUC <- unlist(AUCs[grep(paste0(tolower(mods), collapse = '|'), colnames(AUCs))])
names(modsAUC) <- mods
## convert to df and save only used AUCs
modsAUCdf <- matrix(modsAUC, nrow = 1)
colnames(modsAUCdf) <- names(modsAUC)
write.csv(modsAUCdf, paste0('output/', folder, '/models/evaluation/AUCused.csv'), row.names = F)

## run the weighted ensemble mean based on AUC score
ensProj <- weighted.mean(allProj, w = modsAUC)
names(ensProj) <- 'w.mean'
## compute the variance
ensVar <- app(allProj, var)

## save all to disk
writeRaster(ensProj, filename = paste0('output/', folder, '/projections/rast/ensemble_proj.tif'), gdal = c('COMPRESS=DEFLATE'), overwrite = TRUE)
writeRaster(ensVar, filename = paste0('output/', folder, '/projections/rast/ensemble_variance.tif'), gdal = c('COMPRESS=DEFLATE'), overwrite = TRUE)

## plot them
par(mfrow = c(2, 1))
plot(ensProj, main = ' Ensemble mean weighted on AUC', range = c(0, 1))
plot(ensVar, main = ' Ensemble variance', range = c(0, 1))





################# Clean memory_________________________________
rm(list = ls())
gc()
