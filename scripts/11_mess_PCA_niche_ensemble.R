############################################################
############ MESS, niche analysis, and ensemble ############
############################################################

library(terra)
library(tidyverse)
library(tidyterra)
library(predicts)
library(FactoMineR)
library(factoextra)
library(dismo)
library(pheatmap)
library(patchwork)

######### Load data ----------
folder <- 'baseline_models_yearly_runs' # name of the mother folder

world <- geodata::world(path = 'data/geodata/', resolution = 1)
world <- world[!world$GID_0 == 'ATA'] # world map no Antactica

native <- read.csv(file.path('output', folder, 'modellingData/native_occs.csv'))
invasive <- read.csv(file.path('output', folder, 'modellingData/invasive_occs.csv'))

## models predictions
avgAllObs <- rast(file.path('output', folder, 'models/rast/all_obs_mods_avg_fav.tif'))
avgNative <- rast(file.path('output', folder, 'models/rast/native_mods_avg_fav.tif'))
avgInvasive <- rast(file.path('output', folder, 'models/rast/invasive_mods_avg_fav.tif'))
avgYearly <- rast(file.path('output', folder, 'models/rast/yearly_no_weights_mods_avg_fav.tif'))
## models evaluation
allEvalsAgg <- read.csv(file.path('output', folder, 'models/evaluation/final_evaluation_all_models.csv'))

## predictors
bios <- rast(file.path('output', folder,'modellingData/predictors.tif')) # not cropped
names(bios)[names(bios) == 'hursmean'] <- 'hurs'
## training regions for mess analysis
allObsTrain <-  read.csv(file.path('output', folder, 'modellingData/modData.csv')) # all observations training area
nativeTrain <-  read.csv(file.path('output', folder, 'modellingData/native_modData.csv')) # native training area
invasiveTrain <-  read.csv(file.path('output', folder, 'modellingData/invasive_modData.csv')) # invasive training area
recursiveTrain <- readRDS(file.path('output', folder, 'models/yearly_dataStorage.Rdata')) # recursive training areas
######### 


######### MESS analysis ----------
## all observations
messAllObs <- predicts::mess(x = bios, v = allObsTrain[c('bio01', 'bio04', 'bio07', 'bio12', 'bio15', 'hursmean')])
messAllObs <- ifel(messAllObs >= 0, NA, messAllObs)
plot(messAllObs, col = heat.colors(10), main = 'mess - all observations', font.main = 1)
plot(world, add = T)
names(messAllObs) <- 'messAllObs'

## native
messNative <- predicts::mess(x = bios, v = nativeTrain[c('bio01', 'bio04', 'bio07', 'bio12', 'bio15', 'hursmean')])
messNative <- ifel(messNative >= 0, NA, messNative)
plot(messNative, col = heat.colors(10), main = 'mess - native', font.main = 1)
plot(world, add = T)
names(messNative) <- 'messNative'

## invasive
messInvasive <- predicts::mess(x = bios, v = invasiveTrain[c('bio01', 'bio04', 'bio07', 'bio12', 'bio15', 'hursmean')])
messInvasive <- ifel(messInvasive >= 0, NA, messInvasive)
plot(messInvasive, col = heat.colors(10), main = 'mess - invasive', font.main = 1)
plot(world, add = T)
names(messInvasive) <- 'messInvasive'

## recursive models
years <- names(recursiveTrain)
messRecursive <- NULL
for (y in years) {
  cat('year', y, '\n')
  thisTrain <- recursiveTrain[[y]]
  
  thisMess <- predicts::mess(x = bios, v = thisTrain[c('bio01', 'bio04', 'bio07', 'bio12', 'bio15', 'hursmean')])
  thisMess <- ifel(thisMess >= 0, NA, thisMess)
  plot(thisMess, col = heat.colors(10), main = paste0('mess - ', y), font.main = 1)
  plot(world, add = T)
  names(thisMess) <- as.character(y)
  
  messRecursive <- c(messRecursive, thisMess)
}
messRecursive <- rast(messRecursive)

## save all mess to file
writeRaster(messAllObs, filename = file.path('output', folder, 'mess/all_obs_mess.tif'), 
            gdal = c('COMPRESS=DEFLATE'), overwrite = T)
writeRaster(messNative, filename = file.path('output', folder, 'mess/native_mess.tif'), 
            gdal = c('COMPRESS=DEFLATE'), overwrite = T)
writeRaster(messInvasive, filename = file.path('output', folder, 'mess/invasive_mess.tif'), 
            gdal = c('COMPRESS=DEFLATE'), overwrite = T)
writeRaster(messRecursive, filename = file.path('output', folder, 'mess/recursive_mess.tif'), 
            gdal = c('COMPRESS=DEFLATE'), overwrite = T)
######### 


######### Consensus of the different predictions ----------
stack <- c(avgAllObs[['mods_mean_fav']], 
           avgNative[['mods_mean_fav']], 
           avgInvasive[['mods_mean_fav']])
names(stack) <- c('allObs', 'native', 'invasive')
## add the yearly runs
rec <- avgYearly[[grep('mods_mean_fav_', names(avgYearly), value = T)]]
names(rec) <- gsub('mods_mean_fav_', '', names(rec))
stack <- c(stack, rec)
predictionsDF <- as.data.frame(stack, xy = TRUE, na.rm = TRUE)
predictionsDF <- predictionsDF |> 
  rename(AO = allObs, IV = invasive, NT = native )
rm(rec)

## correlation and clustering of the different models
correlation <- cor(predictionsDF[, -c(1:2)])
corrplot::corrplot(correlation, method = 'ellipse', type = 'upper', addCoefasPercent = F,
                   tl.col = 'black', diag = T)
clust <- hclust(as.dist(1 - correlation))
clustPLOT <- fviz_dend(clust, k = 4, repel = T, type = 'rectangle', rect = T, 
                       lower_rect = -0.15, cex = 1.3) +
  labs(title = '(C) Models Clusters') +
  theme(text = element_text(size = 20), legend.position = 'none')
clustPLOT
######### 


######### PCA ----------
## transform predictor variables to df
biosDF <- as.data.frame(bios)
sum(is.na(biosDF))
## remove bio07 because it is not used in modelling
# biosDF <- biosDF |> 
#   dplyr::select(-bio07)
## divide the observations in the four clusters of correlation analysis
cl1 <- rbind(native, invasive |> filter(year(eventDate) <= 2008)) |> 
  mutate(tag = 'Cl1') # tag each df
cl2 <- rbind(invasive |> filter(year(eventDate) > 2008 & year(eventDate) <= 2015)) |> 
  mutate(tag = 'Cl2') # tag each df
cl3 <- rbind(invasive |> filter(year(eventDate) > 2015 & year(eventDate) <= 2018)) |> 
  mutate(tag = 'Cl3') # tag each df
cl4 <- rbind(invasive |> filter(year(eventDate) > 2018)) |> 
  mutate(tag = 'Cl4') # tag each df
## subset the occurrences so that we have equal numbers of observations in each clusters
# nRows <- min(nrow(cl1), nrow(cl2), nrow(cl3), nrow(cl4))
# set.seed(1312) # maybe repeat multiple times
# clusteredOccurrences <- rbind(slice_sample(cl1, n = nRows), slice_sample(cl2, n = nRows),
#                            slice_sample(cl3, n = nRows), slice_sample(cl4, n = nRows))
clusteredOccurrences <- rbind(cl1, cl2, cl3, cl4)
## extract the values of the predictors at each location
extractedVals <- terra::extract(bios, clusteredOccurrences[, c('long', 'lat')], ID = F)
colSums(is.na(extractedVals))
## merge the points to the extracted values
clusteredOccurrences <- clusteredOccurrences |> bind_cols(extractedVals)
colSums(is.na(clusteredOccurrences))
clusteredOccurrences <- drop_na(clusteredOccurrences)
colSums(is.na(clusteredOccurrences))
## bind the values at the occurrence points to the predictor variables df
nrow(biosDF)
biosDF_occs <- rbind(biosDF, clusteredOccurrences[, names(biosDF)])
## run PCA and check it## run PCA and chbio07eck it
pca <- PCA(biosDF_occs, scale.unit = T,
           ind.sup = (nrow(biosDF)+1):nrow(biosDF_occs), # extracted values at location are used as supplementary individuals (dont influence PCA but are projected onto PCA space)
           graph = F)
fviz_pca_var(pca)
pca$eig
fviz_eig(pca)
(vars <- as.data.frame(pca$var$coord))
fviz_contrib(pca, choice = 'var', axes = 1) # TEMPERATURE ON PC1
fviz_contrib(pca, choice = 'var', axes = 2) # PRECIPITATION ON PC2
## get the coordinates of the supplementary individuals
suppCoords <- as.data.frame(pca$ind.sup$coord) |>
  mutate(tag = clusteredOccurrences$tag)
## calculate the centroid of each cluster
centroids <- suppCoords |>
  group_by(tag) |>
  summarise(Dim.1 = mean(Dim.1, na.rm = TRUE),
            Dim.2 = mean(Dim.2, na.rm = TRUE),
            .groups = 'drop')
## get a percentile of points in each cluster
percentile <- suppCoords %>%
  group_by(tag) |>
  filter(Dim.1 >= quantile(Dim.1, 0.25) & Dim.1 <= quantile(Dim.1, 0.75) &
      Dim.2 >= quantile(Dim.2, 0.25) & Dim.2 <= quantile(Dim.2, 0.75)) |> 
  slice(chull(Dim.1, Dim.2))

pcaPLOT <- fviz_pca_var(pca, col.circle = NA, labelsize = 6)
pcaPLOT <- pcaPLOT +
  # geom_polygon(data = percentile, 
  #              aes(x = Dim.1, y = Dim.2, colour = tag), fill = NA) +
  geom_point(data = suppCoords |> 
               filter(Dim.2 >= -1 & Dim.2 <= 1 & Dim.1 >= -1) |> 
               sample_n(500),
             aes(x = Dim.1, y = Dim.2, colour = tag), alpha = 0.6, shape = 16) +
  geom_point(data = centroids, 
             aes (x = Dim.1, y = Dim.2, fill = tag), size = 7, shape = 23) +
  scale_color_manual(values = c('Cl4' = '#F8766D',
                                'Cl3' = '#7CAE00',
                                'Cl2' = '#00BFC4',
                                'Cl1' = '#C77CFF'),
                     # labels = c('cl4' ='2019-2025',
                     #            'cl3' ='2016-2018',
                     #            'cl2' ='2009-2015',
                     #            'cl1' ='2002-2008'),
                     aesthetics = c('colour', 'fill')) +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  labs(fill = 'Cluster\nCentroids',
       colour = 'Sampled\nOccurrences\n(n = 500)',
       title = '(A) PCA',
       x = paste0('PC1 (', round(pca$eig[1, 2], 2), '%)'),
       y = paste0('PC2 (', round(pca$eig[2, 2], 2), '%)')) +
  theme(text = element_text(size = 20),
        legend.key.size = unit(0.8, 'cm'))
pcaPLOT

## reconstruct the path of the clusters 
## by looking at the median value of the biolcims in each cluster of observations
## I do it only on hurs-bio15 and bio1-bio4 because they have a good gradient along the first two PCs
clustTraj <- clusteredOccurrences|> 
  group_by(tag) |> 
  summarise(bio01_median = median(bio01, na.rm = TRUE),
            bio01_sd = sd(bio01, na.rm = TRUE),
            bio04_median = median(bio04, na.rm = TRUE),
            bio04_sd = sd(bio04, na.rm = TRUE),
            bio07_median = median(bio07, na.rm = TRUE),
            bio07_sd = sd(bio07, na.rm = TRUE),
            bio12_median = median(bio12, na.rm = TRUE),
            bio12_sd = sd(bio12, na.rm = TRUE),
            bio15_median = median(bio15, na.rm = TRUE),
            bio15_sd = sd(bio15, na.rm = TRUE),
            hurs_median = median(hurs, na.rm = TRUE),
            hurs_sd  = sd(hurs, na.rm = TRUE))

long <- clustTraj |> 
  mutate(years = case_when(tag == 'Cl1' ~ '2002-2008',
                           tag == 'Cl2' ~ '2009-2015',
                           tag == 'Cl3' ~ '2016-2018',
                           tag == 'Cl4' ~ '2019-2025')) |> 
  pivot_longer(cols = -c(years, tag),
               names_to = c("variable", "stat"),
               names_pattern = "^(.*)_(median|sd)$",
               values_to = "value") |>
  pivot_wider(names_from = stat, values_from = value)
traj_data <- long |>
  group_by(variable) |>
  arrange(match(tag, c('2002-2008', '2009-2015', '2016-2018', '2019-2025')), .by_group = TRUE) |>
  mutate(xend = lead(tag),
         yend = lead(median)) |>
  filter(!is.na(xend)) |>
  ungroup()

clustTrajPLOT <- ggplot(data = long, aes(x = tag)) +
  geom_errorbar(aes(ymin = median-sd/2, ymax = median+sd/2), colour = 'darkgrey') +
  geom_point(aes(y = median, colour = years), size = 5) +
  scale_colour_manual(labels = c('2019-2025' = 'Cl4',
                                 '2016-2018' = 'Cl3',
                                 '2009-2015' = 'Cl2',
                                 '2002-2008' = 'Cl1'),
                      values = c('2019-2025' = '#F8766D',
                                 '2016-2018' = '#7CAE00',
                                 '2009-2015' = '#00BFC4',
                                 '2002-2008' = '#C77CFF')) +
  geom_segment(data = traj_data,
               aes(x = tag, xend = xend, y = median, yend = yend, group = 1, linewidth = 'Trajectory'),
               arrow = arrow(length = unit(0.25, 'cm'), type = 'open'),
               color = 'black') +
  scale_linewidth_manual(values = c('Trajectory' = 0.5)) +
  facet_wrap(~variable, scales = 'free') +
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  labs(colour = 'Clusters',
       linewidth = NULL,
       title = '(B) Predictors Trajectories',
       x = NULL, y = NULL) +
  theme(strip.background = element_blank(),strip.text = element_text(hjust = 0),
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
clustTrajPLOT

######### 


######### Niche overlap ----------
similarity <- matrix(NA, nrow = nlyr(stack), ncol = nlyr(stack), 
                       dimnames = list(names(stack), names(stack)))
for (i in 1:(nlyr(stack)-1)) {
  for (j in (i+1):nlyr(stack)) {
    ## horizontal is I similarity (right of NA)
    sim_I <- nicheOverlap(raster::raster(stack[[i]]), raster::raster(stack[[j]]), 
                        stat = 'I', mask = F)
    similarity[i, j] <- sim_I
    ## vertical is D similarity (left of NA)
    sim_D <- nicheOverlap(raster::raster(stack[[i]]), raster::raster(stack[[j]]), 
                        stat = 'D', mask = F)
    similarity[j, i] <- sim_D
  }
}
write.csv(similarity, file.path('output', folder, 'niche/similarity.csv'), row.names = F)
similarity <- read.csv(file.path('output', folder, 'niche/similarity.csv'), check.names = F)
new_labels <- c('AO', 'NT', 'IV')
# Rename first 3 columns
colnames(similarity)[1:3] <- new_labels
# Rename first 3 rows
rownames(similarity)[1:3] <- new_labels
customOrder <- c('NT', 2002, 2003, 2007:2025, 'IV', 'AO')
rownames(similarity) <- colnames(similarity)
similarityOrdered <- similarity[customOrder, customOrder]

pheatmap(similarityOrdered,
  cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE,
  number_color = 'black', fontsize_number = 8, color = colorRampPalette(c('white', '#FDD835', '#E53935'))(100),
  main = '(B) Niche Similarity Index (lower D, upper I)')
## heatmap created with ggplot
row_order <- rownames(similarityOrdered)
col_order <- colnames(similarityOrdered)
m <- as.matrix(similarityOrdered)
sim <- data.frame(Var1 = rep(rownames(m), times = ncol(m)),
                  Var2 = rep(colnames(m), each = nrow(m)),
                  value = as.vector(m))
sim$Var1 <- factor(sim$Var1, levels = rev(row_order))  # rev so first row is on top
sim$Var2 <- factor(sim$Var2, levels = col_order)

ggplot(sim, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = 'white') +
  geom_text(aes(label = ifelse(is.na(value), 'NA', round(value, 2))),
            size = 2.5) +
  scale_fill_gradient2(low = 'orange', mid = 'white', high = 'red2', 
                      midpoint = 0.7, limits = c(0.4, 1), na.value = 'grey',
                      name = 'Similarity') +
  guides(fill = guide_colorbar(barheight = unit(6, 'cm'))) +
  theme_minimal() +
  theme(text = element_text(size = 20),
        plot.subtitle = element_text(size = 15),
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
        axis.title = element_blank()) +
  coord_fixed() +
  labs(title = 'Niche Similarity Index',
       subtitle = 'Lower Schoener’s D, upper Warren’s I')
######### 


######### Ensemble of all models ----------
## load mess rasters
messAllObs <- rast(file.path('output', folder, 'mess/all_obs_mess.tif'))
messNative <- rast(file.path('output', folder, 'mess/native_mess.tif'))
messInvasive <- rast(file.path('output', folder, 'mess/invasive_mess.tif'))
messRecursive <- rast(file.path('output', folder, 'mess/recursive_mess.tif'))
## get the average mess
allMess <- c(messAllObs, messInvasive, messNative, messRecursive)
avgMess <- app(allMess, mean)
plot(avgMess)
plot(world, add = T)
## get all predictions in a single stack
allPredictions <- c(avgAllObs[[grep('mean', names(avgAllObs))]],
                    avgInvasive[[grep('mean', names(avgInvasive))]],
                    avgNative[[grep('mean', names(avgNative))]],
                    avgYearly[[grep('mean', names(avgYearly))]])
names(allPredictions[[1:3]]) <- c('allObs', 'invasive', 'native')
names(allPredictions)
names(allMess)
## mask all predictions with its relative mess
allPredictionsMasked <- mask(allPredictions, allMess, inverse = T)
## plot all simulations
defaultPlotOptions <- par(no.readonly = T) # store original plotting options
par(oma = c(0, 0, 0, 2))
plot(allPredictionsMasked[[1:12]], font.main = 1, legend = F)
plot(allPredictionsMasked[[1]],
     legend.only = TRUE,
     axes = F,
     plg = list(x = 194, y = 180, size = c(1.5, 0.2)))
par(oma = c(0, 0, 0, 2))
plot(allPredictionsMasked[[13:24]], font.main = 1, legend = F)
plot(allPredictionsMasked[[1]],
     legend.only = TRUE,
     axes = F,
     plg = list(x = 194, y = 180, size = c(1.5, 0.2)))
par(defaultPlotOptions) # restore plotting window to original

weights <- allEvalsAgg |> 
  dplyr::select(year, auc_mean, boyce_mean) |> 
  mutate(weight = (auc_mean+boyce_mean) / 2)
## compute ensemble weighted on AUC
ensemble <- weighted.mean(allPredictions, w = weights$weight)
plot(ensemble)
## normal ensemble
ensembleNoWeights <- app(allPredictions, mean)
plot(ensembleNoWeights)
## ensemble variance
ensembleVar <- app(allPredictions, var)

## mask ensemble with mess
avgMessBinary <- ifel(avgMess <= 0, 1, 0) ## convert to bianry
plot(avgMessBinary)
maskedEns <- mask(ensemble, avgMessBinary, inverse = T)
maskedEnsVar <- mask(ensembleVar, avgMessBinary, inverse = T)
plot(maskedEns)

ens <- c(maskedEns, maskedEnsVar)
names(ens) <- c('w.ensemble', 'variance')
## plot them nicely
defaultPlotOptions <- par(no.readonly = T) # store original plotting options
par(mfrow = c(2, 1))
plot(ens[[1]], main = 'w.ensemble', font.main = 1, plg = list(size = c(1, 0.3)))
plot(world, lwd = 0.2, border = 'black', add = T)
plot(ens[[2]], main = 'variance', font.main = 1, col = hcl.colors(n = 100, 'Green-Orange'), 
     plg = list(size = c(1, 0.3)))
plot(world, lwd = 0.5, border = 'black', add = T)
par(defaultPlotOptions) # restore plotting window to original
## save to file
writeRaster(ens, filename = file.path('output', folder, 'ensemble/ens.tif'), 
            gdal = c('COMPRESS=DEFLATE'), overwrite = T)

## niche similarity ens VS. allObs
similarity_I <- nicheOverlap(raster::raster(ens[[1]]), raster::raster(avgAllObs[[1]]), 
                             stat = 'I', mask = F)
similarity_D <- nicheOverlap(raster::raster(ens[[1]]), raster::raster(avgAllObs[[1]]), 
                             stat = 'D', mask = F)
## plot ens VS. allObs
defaultPlotOptions <- par(no.readonly = T) # store original plotting options
par(mfrow = c(2, 1))
plot(world, lwd = 0, main = 'w.ensemble', font.main = 1)
plot(ens[[1]], plg = list(size = c(1, 0.5)), add = T)
plot(world, lwd = 0.5, border = 'grey', add = T)
plot(world, lwd = 0, main = 'all observations', font.main = 1)
plot(mask(avgAllObs[[grep('mean', names(avgAllObs))]], messAllObs, inverse = T),
     plg = list(size = c(1, 0.5)), add = T)
plot(world, lwd = 0.5, border = 'grey', add = T)
par(defaultPlotOptions) # restore plotting window to original

## four different ensembles according to the clusters
evasion <- c(allPredictions[['native']], 
             allPredictions[[grep(paste(2002:2008, collapse = '|'), names(allPredictions))]])
evasionMess <- app(c(messNative, 
                     messRecursive[[grep(paste(2002:2008, collapse = '|'), names(messRecursive))]]), mean)
evasionMessBinary <- ifel(evasionMess <= 0, 1, 0)
evasionW <- allEvalsAgg |> 
  filter(year == 0 | year >= 2002 & year <= 2008) |>
  dplyr::select(year, auc_mean, boyce_mean) |>  
  mutate(weight = (auc_mean+boyce_mean) / 2)
evasionEns <- weighted.mean(evasion, w = evasionW$weight)
names(evasionEns) <- 'Evasion'
evasionEns <- mask(evasionEns, evasionMessBinary, inverse = T)

expansion <- allPredictions[[grep(paste(2009:2015, collapse = '|'), names(allPredictions))]]
expansionMess <- app(messRecursive[[grep(paste(2009:2015, collapse = '|'), names(messRecursive))]], mean)
expansionMessBinary <- ifel(expansionMess <= 0, 1, 0)
expansionW <- allEvalsAgg |> 
  filter(year > 2008 & year <= 2015) |> 
  dplyr::select(year, auc_mean, boyce_mean) |> 
  mutate(weight = (auc_mean+boyce_mean) / 2)
expansionEns <- weighted.mean(expansion, w = expansionW$weight)
names(expansionEns) <- 'Expansion'
expansionEns <- mask(expansionEns, expansionMessBinary, inverse = T)

consolidation <- allPredictions[[grep(paste(2016:2018, collapse = '|'), names(allPredictions))]]
consolidationMess <- app(messRecursive[[grep(paste(2016:2018, collapse = '|'), names(messRecursive))]], mean)
consolidationMessBinary <- ifel(consolidationMess <= 0, 1, 0)
consolidationW <- allEvalsAgg |> 
  filter(year > 2015 & year <= 2018) |> 
  dplyr::select(year, auc_mean, boyce_mean) |> 
  mutate(weight = (auc_mean+boyce_mean) / 2)
consolidationEns <- weighted.mean(consolidation, w = consolidationW$weight)
names(consolidationEns) <- 'Consolidation'
consolidationEns <- mask(consolidationEns, consolidationMessBinary, inverse = T)

saturation <- c(allPredictions[['allObs']], allPredictions[['invasive']], 
                allPredictions[[grep(paste(2019:2025, collapse = '|'), names(allPredictions))]])
saturationMess <- app(c(messAllObs, messInvasive,
                        messRecursive[[grep(paste(2019:2025, collapse = '|'), names(messRecursive))]]), mean)
saturationMessBinary <- ifel(saturationMess <= 0, 1, 0)
saturationW <- allEvalsAgg |> 
  filter(year > 2018 & year <= 2025 | year == 1 | year == 2) |> 
  dplyr::select(year, auc_mean, boyce_mean) |> 
  mutate(weight = (auc_mean+boyce_mean) / 2)
saturationEns <- weighted.mean(saturation, w = saturationW$weight)
names(saturationEns) <- 'Saturation'
saturationEns <- mask(saturationEns, saturationMessBinary, inverse = T)

## nice plot with all clusters and a single legend
all_ens <- c(evasionEns, expansionEns,
             consolidationEns, saturationEns)
oldNames <- names(all_ens)
names(all_ens)  <- c('Cl1', 'Cl2', 'Cl3', 'Cl4')
wd <- ggplot(data = world) +
  geom_spatraster(data = all_ens) +
  scale_x_continuous(n.breaks = 4) +
  scale_y_continuous(n.breaks = 4) +
  scale_fill_viridis_c(na.value = 'transparent') +
  geom_spatvector(fill = NA) +
  facet_wrap(~lyr, nrow = 4) +
  labs(title = '(A) Global') +
  theme_minimal() +
  theme(strip.text = element_text(hjust = 0),
        text = element_text(size = 20),
        axis.text = element_text(size = 8))
wd

## create strips with zooms on specific areas
## North America
NorthA <- crop(world, ext(-180,-25,8,90))
plot(NorthA)
NorthAEns <- rast()
for (i in 1:nlyr(all_ens)) {
  thisE <- all_ens[[i]]

  p <- crop(thisE, NorthA)
  # plot(p, main = names(p))
  # plot(world, add = T)
  NorthAEns <- c(NorthAEns, p)
}
na <- ggplot(data = NorthA) +
  geom_spatraster(data = NorthAEns) +
  scale_x_continuous(n.breaks = 4) +
  scale_y_continuous(n.breaks = 4) +
  scale_fill_viridis_c(na.value = 'transparent') +
  geom_spatvector(fill = NA) +
  facet_wrap(~lyr, nrow = 1) +
  labs(title = '(B) North America') +
  theme_minimal() +
  theme(strip.text = element_text(hjust = 0),
        text = element_text(size = 20),
        axis.text = element_text(size = 8))
na

## Europe
Europe <- crop(world, ext(-25,50,30,72))
plot(Europe)
EuropeEns <- rast()
for (i in 1:nlyr(all_ens)) {
  thisE <- all_ens[[i]]
  
  p <- crop(thisE, Europe)
  # plot(p, main = names(p))
  # plot(world, add = T)
  EuropeEns <- c(EuropeEns, p)
}
eu <- ggplot(data = Europe) +
  geom_spatraster(data = EuropeEns) +
  scale_x_continuous(n.breaks = 4) +
  scale_y_continuous(n.breaks = 4) +
  scale_fill_viridis_c(na.value = 'transparent') +
  geom_spatvector(fill = NA) +
  facet_wrap(~lyr, nrow = 1) +
  labs(title = '(C) Europe') +
  theme_minimal() +
  theme(strip.text = element_text(hjust = 0),
        text = element_text(size = 20),
        axis.text = element_text(size = 8))
eu

## Native areas (South East Asia)
SEAsia <- crop(world, ext(65,170,-12,60))
plot(SEAsia)
SEAsiaEns <- rast()
for (i in 1:nlyr(all_ens)) {
  thisE <- all_ens[[i]]
  
  p <- crop(thisE, SEAsia)
  # plot(p, main = names(p))
  # plot(world, add = T)
  SEAsiaEns <- c(SEAsiaEns, p)
}
sa <- ggplot(data = SEAsia) +
  geom_spatraster(data = SEAsiaEns) +
  scale_x_continuous(n.breaks = 4) +
  scale_y_continuous(n.breaks = 4) +
  scale_fill_viridis_c(na.value = 'transparent') +
  geom_spatvector(fill = NA) +
  facet_wrap(~lyr, nrow = 1) +
  labs(title = '(D) Southeast Asia') +
  theme_minimal() +
  theme(strip.text = element_text(hjust = 0),
        text = element_text(size = 20),
        axis.text = element_text(size = 8))
sa

## add all plots together with patchwork
wd/na/eu/sa +
  plot_layout(guides = 'collect', widths = c(1, 2)) &
  labs(fill = NULL)
######### 


######### Overlap with Koppen-geiger ----------
kg <- rast('data/otherLayers/beck_color_coded.tif')
plot(kg)
plot(world, add = T)
climateBeck <- c('Af', 'Am', 'Aw', 'BWh', 'BWk', 'BSh', 'BSk', 'Csa', 'Csb', 'Csc', 'Cwa', 'Cwb', 'Cwc', 'Cfa', 'Cfb',
                 'Cfc', 'Dsa', 'Dsb', 'Dsc', 'Dsd', 'Dwa', 'Dwb', 'Dwc', 'Dwd', 'Dfa', 'Dfb', 'Dfc', 'Dfd', 'ET', 'EF')
classBeck <- data.frame(ID = 1:30, climate = climateBeck)
kg <- as.factor(kg)
levels(kg) <- classBeck
## native
kgNative <- terra::extract(kg, native[, c('long', 'lat')])
climates <- unique(c(kgNative$climate))
target_IDs <- classBeck$ID[classBeck$climate %in% climates]
kg_filtered <- match(kg, target_IDs) 
kg_filtered <- mask(kg, kg_filtered)
plot(kg_filtered, main = 'native')
## invasive
kgInvasive <- terra::extract(kg, invasive[, c('long', 'lat')])
climates <- unique(c(kgInvasive$climate))
target_IDs <- classBeck$ID[classBeck$climate %in% climates]
kg_filtered <- match(kg, target_IDs) 
kg_filtered <- mask(kg, kg_filtered)
plot(kg_filtered, main = 'invasive')
######### 

