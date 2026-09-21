library("terra")
library("embarcadero")
library("randomForest")
library("gam")
library("maxnet")
library("modEvA")
library("blockCV")
library("predicts")
library("fuzzySim")


################# Load Data_________________________________
folder <- 'AllPredictors' # name of the mother folder

bioCropped <- rast(paste0('output/', folder, '/modellingData/ModelRegion/predsCropped.tif'))
data <- read.csv(paste0('output/', folder, '/modellingData/modData.csv')) # train dataset
dataSPAT <- vect(data, geom = c('x', 'y'), crs = 'EPSG:4326') # spatialize it
region <- vect(paste0('output/', folder, '/modellingData/ModelRegion/modellingRegionPoly.shp')) # modelling region polygon
regionEXT <- vect(paste0('output/', folder, '/modellingData/ModelRegion/modellingRegionEXT.shp')) # modelling region extent
world <- geodata::world(path = 'data/geodata/', resolution = 1)[!world$GID_0 == 'ATA'] # world map no Antactica
vars <- read.csv(paste0('output/', folder, '/modellingData/selectedVariables.csv'), # variables after collinearity analysis
                 stringsAsFactors = F, header = T)[, 1]

## check data
plot(bioCropped[[1]])
nrow(data)
table(data$presence)
plot(world, main = 'Presence and background points')
points(subset(data, data$presence == 0, select = c('x', 'y')), col = 'pink3', cex = 0.08)
points(subset(data, data$presence == 1, select = c('x', 'y')), col = 'blue', cex = 0.08)
plot(world, lw = 0.3, add = T)
legend('bottomleft', horiz = T, inset = c(0, 0.06), xpd = T,
       bty = 'n', bg = 'white', cex = 0.8,
       legend = c('presence', 'background'), pch = 21, col = 'black', pt.bg = c('blue', 'pink3'))
vars





################# Blocks cross validation_________________________________
######## create blocks and define the number of folds for cross-validation --------
## METHOD 1: use spatial autocorrelation of points or predictor variables to select blocks size quantitatively
## METHOD 2: empirically set it through iterations and compute the size of the blocks as 
##           the sqrt of the area divided by a factor i (from iteration)
method = 2      # select which method to use for block size computation
fact = 130       # the factor for method 2
if (method == 1){
  dataSPATsf <- sf::st_as_sf(data, coords = c('x', 'y'), crs = 'EPSG:4326')     # transform to sf object
  spatAutocorrPoints <- cv_spatial_autocor(x = dataSPATsf, column = 'presence', 
                                     plot = T, progress = T)
  summary(spatAutocorrPoints)
  blockSize <- spatAutocorrPoints$range
  rm(dataSPATsf)     # remove sf object because not needed
} else if (method == 2) {
  pol <- as.polygons(bioCropped*0)      # transform the whole region to polygon
  plot(pol, col = 'orange', border = NA)
  blockSize <- sqrt(expanse(pol) / fact)
  rm(pol)
}
blockSize/1000 # block area in km

## clean unused memory
gc()


######## define the k folds for cv and divide the region in blocks of k classes --------
k = 5
## create the blocks and divide them in five folds
blocks <- cv_spatial(dataSPAT,
                     column = 'presence', r = bioCropped, size = blockSize, 
                     k = k, selection = 'random', hexagon = T,
                     seed = 1312)
## assign the belonging fold to each observation in the dataset (also the spatial one for plotting)
data$fold <- blocks$folds_ids
dataSPAT$fold <- blocks$folds_ids

## check the number of test data for each fold
(n <- table(data[, c('fold', 'presence')]))
(nt <- rowSums(n))  # n total rows per fold
(np <- n[ , 2])  # n presence rows per fold

## plot the folds and presence points
plot(dataSPAT, 'fold', cex = 0.4, col = hcl.colors(k, 'Peach'), legend = F)
plot(regionEXT, add = T, lw = .4)
points(subset(dataSPAT, dataSPAT$presence == 1), pch = 4, cex = 0.2)
plot(vect(blocks$blocks), border = "blue", add = TRUE, lw = .7)
legend("bottom", inset = c(0, -0.05), xpd = T, bty = 'n',
       legend = paste("block size:", round(blockSize / 1000), "km"))
legend("top", inset = c(0, 0), xpd = T, cex = 0.7, bg = 'white', horiz = T, 
       legend = paste("fold", 1:k, ": N =", np, "/", nt), 
       fill = hcl.colors(k, 'Peach'))

## save the blocked data to csv
write.csv(data.frame(data), paste0('output/', folder, '/modellingData/blockedData.csv'), row.names = F)

## clean unused memory
gc()





################# Train the cv models_________________________________
folds <- sort(unique(data$fold))
## create a df to store results
predsCV <- data[ , c(1:4, grep("fold", names(data)))]
## create a df with as.factor(data$presence) for the RF model
data4rf <- data
data4rf$presence <- as.factor(data4rf$presence)

## iterate through the folds and train the models
source('https://raw.githubusercontent.com/AMBarbosa/unpackaged/master/predict_bart_df')
for (f in folds) {
  startTime <- proc.time() # start time of the fold
  ## train the model on all data but a fold f
  ## but do the prediction using all data
  foldTrainingData <- subset(data, data$fold != f)
  
  message(paste('fold', f, ' - Bioclim'))
  bcMod_fold <- envelope(x = subset(foldTrainingData, presence == 1, select = vars))
  bcPred_fold <- predict(bcMod_fold, data)
  predsCV[ , paste0('BC_fold_', f, '_pred')] <- bcPred_fold
  predsCV[ , paste0('BC_fold_', f, '_pred_reclass')] <- quantReclass(bcPred_fold)
  
  message(paste('fold', f, ' - Maxent'))
  mxtMod_fold <- maxnet(p = foldTrainingData[, 'presence'], data = foldTrainingData[, vars],
                        f = maxnet.formula(foldTrainingData[, 'presence'], foldTrainingData[, vars]))
  mxtPred_fold <- predict(mxtMod_fold, data, type = 'cloglog')
  predsCV[ , paste0('MXT_fold_', f, '_pred')] <- mxtPred_fold
  
  message(paste('fold', f, ' - GLM'))
  glmFormula <- reformulate(termlabels = vars, response = 'presence')
  glmMod_fold <- glm(formula = glmFormula, family = binomial, data = foldTrainingData)
  glmPred_fold <- predict(glmMod_fold, data, type = 'response')
  glmFav_fold <- Fav(pred = glmPred_fold, sample.preval = prevalence(model = glmMod_fold))
  predsCV[ , paste0('GLM_fold_', f, '_pred')] <- glmPred_fold
  predsCV[ , paste0('GLM_fold_', f, '_fav')] <- glmFav_fold
  
  message(paste('fold', f, ' - GAM'))
  gamFormula <- as.formula(paste0('presence ~ ', paste0('s(', vars, ',4)', collapse = '+')))
  gamMod_fold <- gam(gamFormula, family = binomial, data = foldTrainingData)
  gamPred_fold <- predict(gamMod_fold, data, type = 'response')
  gamFav_fold <- Fav(pred = gamPred_fold, sample.preval = prevalence(model = gamMod_fold))
  predsCV[ , paste0('GAM_fold_', f, '_pred')] <- gamPred_fold
  predsCV[ , paste0('GAM_fold_', f, '_fav')] <- gamFav_fold
  
  message(paste('fold', f, ' - BART'))
  bartMod_fold <- bart(x.train = foldTrainingData[, vars], y.train = foldTrainingData[, 'presence'],
                       keeptrees = TRUE, verbose = FALSE, seed = 1312)
  bartPred_fold <- predict_bart_df(bartMod_fold, data)
  bartFav_fold <- Fav(pred = bartPred_fold, sample.preval = prevalence(model = bartMod_fold))
  predsCV[ , paste0('BART_fold_', f, '_pred')] <- bartPred_fold
  predsCV[ , paste0('BART_fold_', f, '_fav')] <- bartFav_fold
  
  message(paste('fold', f, ' - RF'))
  foldTrainingData$presence <- as.factor(foldTrainingData$presence)
  rfFormula <- formula(paste0('presence ~ ', paste(vars, collapse = '+')))
  rfMod_fold <- randomForest(rfFormula, foldTrainingData, na.action = na.exclude)
  rfPred_fold <- predict(rfMod_fold, data4rf, type = 'prob', index = 2)[,2]
  rfFav_fold <- Fav(pred = rfPred_fold, sample.preval = prevalence(model = rfMod_fold))
  predsCV[ , paste0('RF_fold_', f, '_pred')] <- rfPred_fold
  predsCV[ , paste0('RF_fold_', f, '_fav')] <- rfFav_fold

  endTime <- proc.time() # end time of the fold
  message('Elapsed fold ', f, ' time: ', (endTime-startTime)[3], ' sec') # elapsed fold time
  message('')
  gc()
}

## save outputs to disk
write.csv(predsCV, paste0('output/', folder, '/models/crossValidation/crossPreds.csv'), row.names = FALSE)





################# Clean memory_________________________________
rm(list = ls())
gc()
