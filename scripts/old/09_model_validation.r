library("tidyverse")
library("terra")


################# Load data_________________________________
folder <- 'world' # name of the mother folder

dataRaw <- read.csv("data/GBIF/occurrence.csv")
world <- geodata::world(path = 'data/geodata/', resolution = 1)
world <- world[!world$GID_0 == 'ATA'] # world map no Antarctica
region <- vect(paste0('output/', folder, '/modellingData/ModelRegion/modellingRegionPoly.shp')) # modelling region polygon
## data needed for cutoff threshold computation
preds <- read.csv(paste0('output/', folder, '/models/allModelsPredictions.csv')) # prediction in the study area

## single models
mxtProj <- rast(paste0('output/', folder, '/projections/rast/mxt_proj.tif'))
gamProj <- rast(paste0('output/', folder, '/projections/rast/gam_proj.tif'))
gamProjFav <- rast(paste0('output/', folder, '/projections/rast/gam_projFav.tif'))
bartProj <- rast(paste0('output/', folder, '/projections/rast/bart_proj.tif'))
bartProjFav <- rast(paste0('output/', folder, '/projections/rast/bart_projFav.tif'))
par(mfrow = c(2,2))
plot(mxtProj,  main = 'MXT')
plot(gamProjFav,  main = 'GAM')
plot(bartProjFav[[1]],  main = 'BART')
## ensemble
ensProj <- rast(paste0('output/', folder, '/projections/rast/ensemble_proj.tif'))
ensVar <- rast(paste0('output/', folder, '/projections/rast/ensemble_variance.tif'))
par(mfrow = c(2,1))
plot(ensProj, main = ' Ensemble mean weighted on AUC', range = c(0, 1))
plot(region, add = T,  border = 'magenta', col = NA)
plot(ensVar, main = ' Ensemble variance', range = c(0, 1))
plot(region, add = T,  border = 'magenta', col = NA)
par(mfrow = c(1,1))
legend('left', legend = 'Training area', lty = 1,lwd = 1.5, border = F, col = 'magenta', 
       inset = 0.1, cex = 1, bty = 'n', xpd = T)




################# Get occurrences from native areas_________________________________
## clean the dataframe
dataClean <- dataRaw |> 
  select(gbifID, occurrenceID, basisOfRecord, occurrenceStatus, species, continent, eventDate, countryCode,
         decimalLatitude, decimalLongitude, verbatimEventDate, year, month, day, startDayOfYear, endDayOfYear) |> 
  mutate(eventDate = as.Date(eventDate, "%d/%m/%Y")) |> 
  ## make the date from day, month, year when available
  mutate(eventDate = coalesce(eventDate, as.Date(paste(day, month, year, sep = "/"), "%d/%m/%Y"))) |> 
  ## make the date from doy and year when available 
  mutate(eventDate = coalesce(eventDate, as.Date(startDayOfYear-1, origin = paste0(year, "/01/01")))) |>
  ## convert verbatim dates
  mutate(eventDate = coalesce(eventDate, as.Date(verbatimEventDate, format = "%d/%m/%Y %H:%M"))) |> 
  ## drop no coords or date
  drop_na(decimalLatitude, decimalLongitude, eventDate)

## remove duplicates, check for errors, and get native occurrences
nativeOccs <- dataClean |> 
  filter(countryCode %in% c("CN", "TW", "JP", "KP", "KR", "HK", "MO")) |>
  CoordinateCleaner::cc_dupl(lon = "decimalLongitude", lat = "decimalLatitude", addition = "eventDate", value = "clean") |> 
  CoordinateCleaner::clean_coordinates(lon = "decimalLongitude", lat = "decimalLatitude", countries = "countryCode", 
                                       tests = c("centroids", "equal", "gbif", "outliers", "seas", "zeros"), value = "clean") |> 
  select(gbifID, occurrenceID, basisOfRecord, occurrenceStatus, species, continent, eventDate, countryCode,
         decimalLatitude, decimalLongitude)
write.csv(nativeOccs, paste0('output/', folder, '/modellingData/nativeOccs.csv'), row.names = F)

## plot native occs
nativeOccsSP <- terra::vect(nativeOccs, geom = c('decimalLongitude', 'decimalLatitude'), crs = 'epsg:4326')

plot(world)
points(nativeOccsSP, col = 'blue', cex = .5, pch = 4)
## zoom
plot(ensProj, ext = ext(nativeOccsSP)+10, main = 'ensemble focus on native regions')
plot(nativeOccsSP, col = 'red', cex = .4, ext = ext(nativeOccsSP)+10, add = T)
plot(world, add = T)
inset(ensProj, lwd = 0.1, scale = 0.3, loc = "bottomright",
      box = ext(nativeOccsSP)+10, pbox = list(col = "red", lwd = 2), background = "white")





################# Model projection evaluation_________________________________
######## custom function to calculate cutoff maximising TSS -------- 
maxTSScutoff <- function(presence, prediction, valueStart = 0, valueEnd = 1, nsplits = 101) {
  ## create a df to store results
  output <- data.frame(matrix(nrow = nsplits, ncol = 2))
  colnames(output) <- c('thresh', 'TSS')
  ## divide in splits
  splits <- seq(valueStart, valueEnd, length.out = nsplits)
  ## for each cutoff thresh
  for (ls in 1:length(splits)) {
    ## get the value of the split
    split = splits[ls]
    ## get the prediction results
    predPresence <- as.numeric(prediction >= split)
    
    ## compute confusion matrix
    TP = sum(predPresence == 1 & presence == 1) ## true positive
    TN = sum(predPresence == 0 & presence == 0) ## true negative
    FP = sum(predPresence == 1 & presence == 0) ## false positive
    FN = sum(predPresence == 0 & presence == 1) ## false negative
    
    ## compute sensitivity and specificity
    sensitivity = ifelse((TP + FN) == 0, NA, TP / (TP + FN)) # true positive rate
    specificity = ifelse((FP + TN) == 0, NA, TN / (FP + TN)) # true negative rate
    ## compute TSS
    TSS = sensitivity + specificity - 1
    
    # append result to df
    output$thresh[ls] <- split
    output$TSS[ls] <- TSS
  }
  ## extract the cutoff thresh at which TSS is highest
  cutoffThresh = output$thresh[output$TSS == max(output$TSS)]
  cat(paste('Maximum TSS:', max(output$TSS), '\nCutoff threshold:', cutoffThresh))
  ## return it
  return(cutoffThresh)
}


######## calculate cutoff threshold based on the predictions in the study area --------
## maxent
mxtCut <- maxTSScutoff(presence = preds[, 'presence'], prediction = preds[, 'mxtPred'], nsplits = 1001)
## maxent
gamCut <- maxTSScutoff(presence = preds[, 'presence'], prediction = preds[, 'gamFav'], nsplits = 1001)
## maxent
bartCut <- maxTSScutoff(presence = preds[, 'presence'], prediction = preds[, 'bartFav'], nsplits = 1001)

## select the minimum cutoff
cutoff <- min(c(mxtCut, gamCut, bartCut))
cutoff


######## apply cutoff threshold to the validation dataset to calculate metrics --------
######## my validation dataset has only presences so I compute sensitivity (true positive rate) which doesn't require absences

## extract predicted values at presence locations
gridded <- terra::extract(ensProj, nativeOccsSP, ID = F, xy = T)
gridded <- gridded |> 
  select(x, y, w.mean) |> 
  rename(predictedValue = w.mean) |> 
  mutate(presence = 1,
         cutoffPres = as.numeric(predictedValue >= cutoff))

## plot the binary map
binaryMap <- ifel(ensProj >= cutoff, 1, 0)
levels(binaryMap) <- data.frame(id = c(1, 0), value = c('presence', 'absence'))
plot(binaryMap, col = c('white', 'gold3'), main = 'Predicted presence/absence')
plot(world, add = T, lwd = 0.4)

## compute true positive
truePos <- sum(gridded[, 'cutoffPres'] == 1 & gridded[, 'presence'] == 1)
## compute false negative
falseNeg <- sum(gridded[, 'cutoffPres'] == 0 & gridded[, 'presence'] == 1)
## compute sensitivity
sensitivity <- truePos / (truePos + falseNeg)
mtext(paste0('With a cutoff threshold of ', round(cutoff, 3), ', the sensitivity of the enselble is ', round(sensitivity, 3)),
     side = 1, line = 3)





################# Clean memory_________________________________
rm(list = ls())
gc()

