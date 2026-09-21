#############################################################################################
############ Run the models with all observation and project on the entire world ############
#############################################################################################

library(terra)
library(mgcv)
library(fuzzySim)
library(dbarts)
library(embarcadero)
library(randomForest)

######### Load data ----------
folder <- 'climate_human' # name of the mother folder

preds <- rast(file.path('output', folder, 'modellingData/predictors.tif'))
data <- read.csv(file.path('output', folder, 'modellingData/modData.csv')) # train dataset
dataSPAT <- vect(data, geom = c('x', 'y'), crs = 'EPSG:4326') # spatialize it

world <- geodata::world(path = 'data/geodata/', resolution = 1)
world <- world[!world$GID_0 == 'ATA'] # world map no Antactica

## load data from only_climate folder
region <- vect(file.path('output', 'only_climate', 'modellingData/ModelRegion/modellingRegionPoly.shp')) # modelling region polygon

selectedMods <- read.csv(file.path('output', 'only_climate', 'models/crossValidation/selectedModels.csv'))

vars <- read.csv(file.path('output', 'only_climate', 'modellingData/selectedVariables.csv'))[, 2] # variables after collinearity analysis
## add the human layers to vars
vars <- c(vars, 'air', 'blt', 'prt', 'rds', 'rdt')

## check the data
names(preds)
vars
######### 


######### Run complete models on world ----------
## selected models
selectedMods$mod
## time tracking
(start <- Sys.time())

######## GAM - PRESENCE/ABSENCE (non-linear fitting)
## formula
gamForm <- as.formula(paste0('presence ~ ', paste0('s(', vars, ', bs = \'cr\')', collapse = '+')))
## model
gamMod <- gam(gamForm, family = binomial, data = data, model = T)
## prediction
gamPred <- predict(preds, gamMod, type = 'response', na.rm = T,
                   filename = file.path('output', folder, 'models/rast/gam_predWorld.tif'), overwrite = T, 
                   wopt = list(gdal = c('COMPRESS=DEFLATE')))
## plot
plot(gamPred, main = 'GAM')
plot(world, lwd = .2, add = T)
plot(region, add = T, border = 'magenta', col = NA)
## convert to prevalence-independent favourability
gamFav <- Fav(pred = gamPred, sample.preval = prevalence(model = gamMod))
## plot
plot(gamFav, main = 'GAM favourability', range = c(0, 1))
plot(world, lwd = .2, add = T)
plot(region, add = T, border = 'magenta', col = NA)
## save to disk
writeRaster(gamFav, filename = file.path('output', folder, 'models/rast/gam_favWorld.tif'), 
            gdal = c('COMPRESS=DEFLATE'), overwrite = T)
## clean unused memory
gc()

######## BART - PRESENCE/ABSENCE (bayesian sum of trees)
## model
bartMod <- bart(x.train = data[, vars], y.train = data[, 'presence'],
                ntree = 100, keeptrees = T, verbose = F, seed = 1312)
## prediction
bartPred <- predict2.bart(bartMod, raster::stack(preds), quantiles = c(.05, .95), 
                          splitby = 10, quiet = T)
summary(bartPred)
names(bartPred) <- c('BART.pred', 'BART.lower.interval', 'BART.upper.interval')
bartPred <- rast(bartPred)
## convert to prevalence-independent favourability
bartFav <- Fav(pred = bartPred, sample.preval = prevalence(model = bartMod))
names(bartFav) <- c('BART.fav', 'BART.fav.lower.interval', 'BART.fav.upper.interval')
## add the uncertainty range
bartPred$BART.pred.uncertainty.range <- bartPred[[3]] - bartPred[[2]] # to pred
bartFav$BART.fav.uncertainty.range <- bartFav[[3]] - bartFav[[2]]    # to fav
## plot
plot(bartPred)
plot(bartFav)
plot(bartFav[[1]], main = 'BART mean posterior favourability')
plot(world, lwd = .2, add = T)
plot(region, add = T, border = 'magenta', col = NA)
## BART MODEL CANNOT be saved to rds because it looses information
writeRaster(bartPred, filename = file.path('output', folder, 'models/rast/bart_predWorld.tif'), 
            gdal = c('COMPRESS=DEFLATE'), overwrite = T)
writeRaster(bartFav, filename = file.path('output', folder, 'models/rast/bart_favWorld.tif'), 
            gdal = c('COMPRESS=DEFLATE'), overwrite = T)
## clean unused memory
gc()

######## Random Forest - PRESENCE/ABSENCE (average of trees)
data4rf <- data
data4rf$presence <- as.factor(data4rf$presence)
## formula
rfForm <- formula(paste0('presence ~ ', paste(vars, collapse = '+')))
## model
rfMod <- randomForest(rfForm, data4rf, na.action = na.exclude)
## prediction
rfPred <- predict(preds, rfMod, type = 'prob', na.rm = T, index = 2, cores = 10, cpkgs = 'randomForest',
                  filename = file.path('output', folder, 'models/rast/rf_predWorld.tif'), overwrite = T,
                  wopt = list(gdal = c('COMPRESS=DEFLATE')))
## plot
plot(rfPred, main = 'RF')
plot(world, lwd = .2, add = T)
plot(region, add = T, border = 'magenta', col = NA)
## convert to prevalence-independent favourability
rfFav <- Fav(pred = rfPred, sample.preval = prevalence(model = rfMod))
## plot
plot(rfFav, main = 'RF favourability')
plot(world, lwd = .2, add = T)
plot(region, add = T, border = 'magenta', col = NA)
## save to disk
writeRaster(rfFav, filename = file.path('output', folder, 'models/rast/rf_favWorld.tif'), 
            gdal = c('COMPRESS=DEFLATE'), overwrite = T)
## clean unused memory
gc()

## time tracking
(end <- Sys.time())
(length <- difftime(end, start, units = 'mins'))
######### 

## ensemble
ens <- mean(bartFav[[1]], gamFav, rfFav)
writeRaster(ens, filename = file.path('output', folder, 'models/rast/ensemble.tif'), 
            gdal = c('COMPRESS=DEFLATE'), overwrite = T)
plot(ens, main = 'ens ~ climate + human')

######## Clean memory --------
rm(list = setdiff(ls(), c('mxtMod', 'gamMod', 'bartMod', 'rfMod')))
gc()
