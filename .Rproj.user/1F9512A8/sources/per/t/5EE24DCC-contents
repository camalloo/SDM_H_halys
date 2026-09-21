#####################################################
############ Explore and analyse results ############
#####################################################

library(tidyverse)
library(terra)
library(fuzzySim)
library(pROC)
library(modEvA)
library(patchwork)

######### Load data ----------
folder <- 'baseline_models_yearly_runs' # name of the mother folder

world <- geodata::world(path = 'data/geodata/', resolution = 1)
world <- world[!world$GID_0 == 'ATA'] # world map no Antactica

native <- read.csv(file.path('output', folder, 'modellingData/native_occs.csv')) |> 
  mutate(year = year(eventDate),
         status = 'n')
invasive <- read.csv(file.path('output', folder, 'modellingData/invasive_occs.csv')) |> 
  mutate(year = year(eventDate),
         status = 'i') |> 
  ## filter obs before 1996 (first invasive discovery)
  filter(year >= 1996)
plot(world, main = 'observations')

GBIF2026 <- read.csv('data/observations/GBIF_validation/occurrence.csv')

dataStorage <- readRDS(file.path('output', folder, 'models/yearly_dataStorage.Rdata'))
region <- vect(file.path('output', folder, 'modellingData/ModelRegion/modellingRegionPoly.shp'))

## all observations models
avgAllObs <- rast(file.path('output', folder, 'models/rast/all_obs_mods_avg_fav.tif'))
plot(avgAllObs[[1]], main = 'all observations model (AO)', font.main = 1,
     plg = list(size = c(1, 0.5)))
allObsVarimp <- read.csv(file.path('output', folder, 'models/evaluation/all_obs_varimp.csv'))
allObsRc <- read.csv(file.path('output', folder, 'models/evaluation/all_obs_response_curves.csv'))
## native models
avgNative <- rast(file.path('output', folder, 'models/rast/native_mods_avg_fav.tif'))
plot(avgNative[[1]], main = 'native model (NT)', font.main = 1,
     plg = list(size = c(1, 0.5)))
nativeVarimp <- read.csv(file.path('output', folder, 'models/evaluation/native_varimp.csv'))
nativeRc <- read.csv(file.path('output', folder, 'models/evaluation/native_response_curves.csv'))
## invasive models
avgInvasive <- rast(file.path('output', folder, 'models/rast/invasive_mods_avg_fav.tif'))
plot(avgInvasive[[1]], main = 'invasive model (IV)', font.main = 1,
     plg = list(size = c(1, 0.5)))
invasiveVarimp <- read.csv(file.path('output', folder, 'models/evaluation/invasive_varimp.csv'))
invasiveRc <- read.csv(file.path('output', folder, 'models/evaluation/invasive_response_curves.csv'))
## yearly models
avgYearly <- rast(file.path('output', folder, 'models/rast/yearly_no_weights_mods_avg_fav.tif'))
plot(avgYearly[[grep('2025', names(avgYearly))]])
yearlyVarimp <- read.csv(file.path('output', folder, 'models/evaluation/yearly_no_weights_varImp.csv'))
yearlyVarimp$Variable <- gsub('hursmean', 'hurs', yearlyVarimp$Variable)
yearlyRc <- read.csv(file.path('output', folder, 'models/evaluation/yearly_no_weights_Rc.csv'))
yearlyRc$var <- gsub('hursmean', 'hurs', yearlyRc$var)
######### 


######### Clean 2026 GBIF observations to validate models ----------
dataClean <- GBIF2026 |> 
  dplyr::select(gbifID, occurrenceID, basisOfRecord, occurrenceStatus, species, continent, eventDate, countryCode,
         decimalLatitude, decimalLongitude, verbatimEventDate, year, month, day, startDayOfYear, endDayOfYear) |> 
  mutate(eventDate = as.Date(eventDate, '%d/%m/%Y')) |> 
  ## make the date from day, month, year when available
  mutate(eventDate = coalesce(eventDate, as.Date(paste(day, month, year, sep = '/'), '%d/%m/%Y'))) |> 
  ## make the date from doy and year when available 
  mutate(eventDate = coalesce(eventDate, as.Date(startDayOfYear-1, origin = paste0(year, '/01/01')))) |>
  ## convert verbatim dates
  mutate(eventDate = coalesce(eventDate, as.Date(verbatimEventDate, format = '%d/%m/%Y %H:%M'))) |> 
  ## drop no coords or date
  drop_na(decimalLatitude, decimalLongitude, eventDate)
## remove duplicated records and check for other error in coordinates --> FINAL DF
dataClean <- dataClean |> 
  CoordinateCleaner::cc_dupl(lon = 'decimalLongitude', lat = 'decimalLatitude', addition = 'eventDate', value = 'clean') |> 
  CoordinateCleaner::clean_coordinates(lon = 'decimalLongitude', lat = 'decimalLatitude', countries = 'countryCode', 
                                       tests = c('centroids', 'equal', 'gbif', 'outliers', 'seas', 'zeros'),
                                       value = 'clean') |> 
  dplyr::select(gbifID, occurrenceID, basisOfRecord, occurrenceStatus, species, continent, eventDate, year, month, day, countryCode,
         lat = decimalLatitude, long = decimalLongitude)
## additional cleaning with FuzzySim
dataClean <- fuzzySim::cleanCoords(dataClean, coord.cols = c('long', 'lat'), uncert.col = NULL)
#########  


######### Create a set of background points for validation ----------
## load the region and sample background points
spatValidation <- vect(dataClean, geom = c('long', 'lat'), crs = crs(world))
## plot
plot(world, main = '2026 validation data', font.main = 1)
points(spatValidation, pch = 20, cex = 0.3, col = 'green')

## create a grid of presence and absence points on the world
stackedEns <- c(avgAllObs[['mods_mean_fav']], avgYearly[['mods_mean_fav_2025']],
                avgNative[['mods_mean_fav']], avgInvasive[['mods_mean_fav']]) ## stack the ensembles
bg <- gridRecords(rst = stackedEns, pres.coords = spatValidation)
table(bg$presence)
## subset absence points
# bg <- selectAbsences(bg, sp.cols = 'presence', coord.cols = c('x', 'y'), CRS = crs(spatValidation),
#                      n = 200000, seed = 1, plot = F, df = T, verbosity = 0)
bg <- bg |>
  filter(presence == 0) |> 
  dplyr::select(presence, x, y)
table(bg$presence)
bgSPAT <- vect(bg, geom = c('x', 'y'), crs = crs(world))
plot(world, main = 'validation dataset')
points(bgSPAT, col = 'blue', cex = 0.5)
points(spatValidation, col = 'red', cex = 0.5)
#########  


######### Validate iterative recursive models ----------
## prepare observation dataframes for validation
valInvasive <- native |>  # invasive are validated with native
  mutate(presence = 1) |> 
  dplyr::select(presence, x = long, y = lat)
valNative <- invasive |> # native are validated with invasive
  mutate(presence = 1) |> 
  dplyr::select(presence, x = long, y = lat)
val2026 <- dataClean |>
  mutate(presence = 1) |> # 2026 data used in all validations
  dplyr::select(presence, x = long, y = lat)

(years <- as.integer(names(dataStorage)))
## run validation script
iter <- 10                            # define the number of iteration for the validation
allEvals <- data.frame()              # store the results
sampleSize <- min(nrow(valInvasive),  # use the same number of presences for all models in validation
                  nrow(valNative), 
                  nrow(val2026))

set.seed(1312) # seed for random sampling of validation presences
for (idx in 1:iter) {
  cat('iteration', idx, '\n')
  for (yr in 1:length(years)) {
    y <- years[yr]
    cat('\t\tEvaluating model year', y, '\n')
    
    ## get the validation set
    if (y == max(years)) { # if last year evaluate with 2026 GBIF records
      val <- val2026 |>
        slice_sample(n = sampleSize)
      val <- val2026
    } else {
      val <- invasive |> 
        filter(year > y) |> 
        mutate(presence = 1) |> 
        dplyr::select(presence, x = long, y = lat)
      val <- rbind(val, val2026) |>
        slice_sample(n = sampleSize)
      val <- rbind(val, val2026)
    }
    valSPAT <- vect(val, geom = c('x', 'y'), crs = crs(world))
    
    ## 1. handle absences maintaining a 1:1 prevalence ratio
    # bgrnd <- bg |> slice_sample(n = sampleSize)       # OR
    ## 2. handle absences maintaining a 10:1 prevalence ratio
    bgrnd <- bg |> slice_sample(n = nrow(val)*10)
    bgrndSPAT <- vect(bgrnd, geom = c('x', 'y'), crs = crs(world))
    
    ## get the ensemble prediction for this year to be evaluated
    toEval <- avgYearly[[paste0('mods_mean_fav_', y)]]
    
    ## extract prediction at validation set
    prediction <- gridRecords(rst = toEval, pres.coords = valSPAT, abs.coords = bgrndSPAT)
    names(prediction) <- c('presence', 'x', 'y', 'cells', 'pred')
    
    ## AUC
    roc <- roc(response = prediction$presence, 
               predictor = prediction$pred, 
               quiet = TRUE)
    auc <- as.numeric(auc(roc))
    auc_ci <- as.numeric(ci.auc(roc))  # lower, mean, upper
    
    ## Boyce index
    boyce <- Boyce(obs = prediction$presence, pred = prediction$pred, plot = F)$B
    
    ## Miller calibration
    miller <- MillerCalib(obs = prediction$presence, pred = prediction$pred, plot = F)$slope
    
    ## store results
    allEvals <- rbind(allEvals,
                      data.frame(iteration = idx,
                                 tag = 'yearly',
                                 year = y,
                                 auc = auc,
                                 auc_lower = auc_ci[1],
                                 auc_upper = auc_ci[3],
                                 boyce = boyce,
                                 miller = miller))
  }
}
#########  


######### Validate all observations model ----------
set.seed(1312)
for (idx in 1:iter) {
  cat('iteration', idx, '\n')
  
  ## validation presences
  val <- dataClean |>  ## use 2026 GBIF records to validate it
    mutate(presence = 1) |> 
    dplyr::select(presence, x = long, y = lat) |> 
    slice_sample(n = sampleSize)
  valSPAT <- vect(val, geom = c('x', 'y'), crs = crs(world))
  
  ## 1. handle absences maintaining a 1:1 prevalence ratio
  # bgrnd <- bg |> slice_sample(n = sampleSize)       # OR
  ## 2. handle absences maintaining a 10:1 prevalence ratio
  bgrnd <- bg |> slice_sample(n = nrow(val)*10)
  bgrndSPAT <- vect(bgrnd, geom = c('x', 'y'), crs = crs(world))
  
  ## extract prediction at validation presences
  prediction <- gridRecords(rst = avgAllObs[['mods_mean_fav']], pres.coords = valSPAT, abs.coords = bgrndSPAT)
  names(prediction) <- c('presence', 'x', 'y', 'cells', 'pred')
  ## AUC
  roc_all <- roc(response = prediction$presence, 
                 predictor = prediction$pred, 
                 quiet = TRUE)
  auc_all <- as.numeric(auc(roc_all))
  auc_all_ci <- as.numeric(ci.auc(roc_all))  # lower, mean, upper
  ## Boyce
  boyce <- Boyce(obs = prediction$presence, pred = prediction$pred, plot = F)$B
  ## append
  allEvals <- rbind(allEvals,
                    data.frame(iteration = idx,
                               tag = 'all_obs_model',
                               year = 2,
                               auc = auc_all,
                               auc_lower = auc_all_ci[1],
                               auc_upper = auc_all_ci[3],
                               boyce = boyce,
                               miller = miller))
}
#########


######### Validate native model ----------
set.seed(1312)
for (idx in 1:iter) {
  cat('iteration', idx, '\n')
  
  ## validation presences
  val <- rbind(valNative, val2026) |> 
    slice_sample(n = sampleSize)
  valSPAT <- vect(val, geom = c('x', 'y'), crs = crs(world))
  
  ## 1. handle absences maintaining a 1:1 prevalence ratio
  # bgrnd <- bg |> slice_sample(n = sampleSize)       # OR
  ## 2. handle absences maintaining a 10:1 prevalence ratio
  bgrnd <- bg |> slice_sample(n = nrow(val)*10)
  bgrndSPAT <- vect(bgrnd, geom = c('x', 'y'), crs = crs(world))
  
  ## extract prediction at validation presences
  prediction <- gridRecords(rst = avgNative[['mods_mean_fav']], pres.coords = valSPAT, abs.coords = bgrndSPAT)
  names(prediction) <- c('presence', 'x', 'y', 'cells', 'pred')
  ## AUC
  roc_all <- roc(response = prediction$presence, 
                 predictor = prediction$pred, 
                 quiet = TRUE)
  auc_all <- as.numeric(auc(roc_all))
  auc_all_ci <- as.numeric(ci.auc(roc_all))  # lower, mean, upper
  ## Boyce
  boyce <- Boyce(obs = prediction$presence, pred = prediction$pred, plot = F)$B
  ## append
  allEvals <- rbind(allEvals,
                    data.frame(iteration = idx,
                               tag = 'native_model',
                               year = 0,
                               auc = auc_all,
                               auc_lower = auc_all_ci[1],
                               auc_upper = auc_all_ci[3],
                               boyce = boyce,
                               miller = miller))
}
#########


######### Validate invasive model ----------
set.seed(1312)
for (idx in 1:iter) {
  cat('iteration', idx, '\n')
  
  ## validation presences
  val <- rbind(valInvasive, val2026) |> 
    slice_sample(n = sampleSize)
  valSPAT <- vect(val, geom = c('x', 'y'), crs = crs(world))
  
  ## 1. handle absences maintaining a 1:1 prevalence ratio
  # bgrnd <- bg |> slice_sample(n = sampleSize)       # OR
  ## 2. handle absences maintaining a 10:1 prevalence ratio
  bgrnd <- bg |> slice_sample(n = nrow(val)*10)
  bgrndSPAT <- vect(bgrnd, geom = c('x', 'y'), crs = crs(world))
  
  ## extract prediction at validation presences
  prediction <- gridRecords(rst = avgInvasive[['mods_mean_fav']], pres.coords = valSPAT, abs.coords = bgrndSPAT)
  names(prediction) <- c('presence', 'x', 'y', 'cells', 'pred')
  ## AUC
  roc_all <- roc(response = prediction$presence, 
                 predictor = prediction$pred, 
                 quiet = TRUE)
  auc_all <- as.numeric(auc(roc_all))
  auc_all_ci <- as.numeric(ci.auc(roc_all))  # lower, mean, upper
  ## Boyce
  boyce <- Boyce(obs = prediction$presence, pred = prediction$pred, plot = F)$B
  ## append
  allEvals <- rbind(allEvals,
                    data.frame(iteration = idx,
                               tag = 'invasive_model',
                               year = 1,
                               auc = auc_all,
                               auc_lower = auc_all_ci[1],
                               auc_upper = auc_all_ci[3],
                               boyce = boyce,
                               miller = miller))
}

## aggregate all evals
allEvalsAgg <- allEvals |> 
  group_by(tag, year) |> 
  summarise(auc_mean = mean(auc),
            auc_sd = sd(auc),
            boyce_mean = mean(boyce),
            boyce_sd = sd(boyce),
            mill_mean = mean(miller),
            mill_sd = sd(miller))
write.csv(allEvalsAgg, 
          file.path('output', folder, 
                    'models/evaluation/final_evaluation_all_models.csv'),
          row.names = F)
#########  


######### Plot performance metrics ----------
allEvalsAgg <- read.csv(file.path('output', folder, 'models/evaluation/final_evaluation_all_models.csv'))
years <- unique(allEvalsAgg$year)

## CBI
cbiPLOT <- ggplot(allEvalsAgg |> dplyr::filter(tag == 'yearly' & !year %in% c(0, 1, 2)),
                  aes(x = year)) +
  geom_ribbon(aes(y = boyce_mean, ymin = boyce_mean-boyce_sd, ymax = boyce_mean+boyce_sd), 
              colour = 'blue', fill = 'blue', linewidth = 0.1, alpha = 0.4) +
  geom_line(aes(y = boyce_mean), colour = 'blue') +
  # geom_hline(aes(yintercept = 1), linetype = 'solid', colour = 'black') +
  geom_hline(aes(yintercept = allEvalsAgg$boyce_mean[allEvalsAgg$year == 2], linetype = 'AO')) +
  geom_hline(aes(yintercept = allEvalsAgg$boyce_mean[allEvalsAgg$year == 0], linetype = 'NT')) +
  geom_hline(aes(yintercept = allEvalsAgg$boyce_mean[allEvalsAgg$year == 1], linetype = 'IV')) +
  scale_x_continuous(breaks = years) +
  scale_y_continuous(limits = c(0.3, 1)) +
  scale_linetype_manual(values = c('AO' = 'dotdash',
                                   'NT' = 'dashed',
                                   'IV' = 'dotted'),
                        breaks = c('NT', 'IV', 'AO'),
                        labels = c('AO' = 'All Obs. (AO)',
                                   'NT' = 'Native (NT)',
                                   'IV' = 'Invasive (IV)')) +
  labs(x = 'Training Year',
       y = 'Metric Score',
       title = '(A) CBI',
       linetype = 'Static SDMs') +
  theme_bw() +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))

## AUC
aucPLOT <- ggplot(allEvalsAgg |> dplyr::filter(tag == 'yearly' & !year %in% c(0, 1, 2)),
                  aes(x = year)) +
  geom_ribbon(aes(y = auc_mean, ymin = auc_mean-auc_sd, ymax = auc_mean+auc_sd), 
              colour = 'red', fill = 'red', linewidth = 0.1, alpha = 0.4) +
  geom_line(aes(y = auc_mean), colour = 'red') +
  geom_hline(aes(yintercept = allEvalsAgg$auc_mean[allEvalsAgg$year == 2], linetype = 'AO')) +
  geom_hline(aes(yintercept = allEvalsAgg$auc_mean[allEvalsAgg$year == 0], linetype = 'NT')) +
  geom_hline(aes(yintercept = allEvalsAgg$auc_mean[allEvalsAgg$year == 1], linetype = 'IV')) +
  scale_x_continuous(breaks = years) +
  scale_y_continuous(limits = c(0.3, 1)) +
  scale_linetype_manual(values = c('AO' = 'dotdash',
                                   'NT' = 'dashed',
                                   'IV' = 'dotted'),
                        breaks = c('NT', 'IV', 'AO'),
                        labels = c('AO' = 'All Obs. (AO)',
                                   'NT' = 'Native (NT)',
                                   'IV' = 'Invasive (IV)')) +
  labs(x = 'Training Year',
       y = 'Metric Score',
       title = '(B) AUC',
       linetype = 'Static SDMs') +
  theme_bw() +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))

evaluationPLOTS <- cbiPLOT+aucPLOT + plot_layout(guides = 'collect', axes = 'collect')
evaluationPLOTS

# "The recursive SDM required a minimum training dataset of approximately 400–1000 georeferenced invasive occurrences (reached around 2015–2016) 
# before achieving reliable predictive discrimination of future invasion locations. 
# Below this threshold, evaluation metrics were unstable and near-random, reflecting fundamental data limitation rather than model failure. 
# After this threshold, the recursive model converged rapidly to the predictive performance of a model trained on the complete 24-year invasion record, 
# demonstrating that the species' invaded climatic niche was effectively characterised by 2019."


######### Analyse the variable importance and response curve evolution ----------
## response curves of yearly models
ggplot(yearlyRc, aes(x = x, y = y)) +
  geom_line(aes(color = factor(year))) +
  facet_grid(mod~var, scale = 'free') + 
  scale_x_continuous(n.breaks = 4) +
  labs(title = '(A) Partial Dependence - Dynamic Models',
       x = 'Variable',
       y = 'Prediction',
       color = 'Year') +
  guides(colour = guide_legend(ncol = 1)) +
  theme_bw() +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 70, vjust = 1, hjust = 1))
  
## variable importance yearly models
# ggplot(yearlyVarimp, aes(x = Importance, y = reorder(Variable, Importance))) +
#   geom_col() +
#   facet_grid(year~mod) + ggtitle('variable importance bars')
ggplot(yearlyVarimp, aes(x = year, y = Importance)) +
  geom_line(aes(colour = Variable)) +
  facet_wrap(~mod) + 
  scale_x_continuous(n.breaks = 4) +
  labs(title = '(B) Variable Importance - Dynamic Models',
       x = 'Year',
       y = 'Importance',
       color = 'Variable') +
  theme_bw() +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 70, vjust = 1, hjust = 1))

## same plot structure as baseline of varimp plot
ggplot(yearlyVarimp, aes(x = Variable, y = Importance, group = factor(year))) +
  geom_line(aes(colour = factor(year))) +
  geom_point(aes(colour = factor(year))) +
  facet_wrap(~mod) + 
  labs(title = '(B) Variable Importance - Dynamic Models',
       x = 'Variable',
       y = 'Importance',
       color = 'Year') +
  theme_bw() +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1))

## variable importance of baseline models
baselineVarimp <- rbind(allObsVarimp |> mutate(tag = 'AO'), 
                        nativeVarimp |> mutate(tag = 'NT'), 
                        invasiveVarimp |> mutate(tag = 'IV')) |> 
  dplyr::select(mod,tag, Variable, Importance) |> 
  mutate(Variable = gsub('hursmean', 'hurs', Variable))

ggplot(baselineVarimp, aes(x = Variable, y = Importance, group = tag, colour = tag)) +
  geom_line(aes(colour = tag)) +
  geom_point()+
  facet_wrap(~mod) + 
  labs(title = '(B) Variable Importance - Static Models',
       x = 'Variable',
       y = 'Importance',
       color = 'Model') +
  theme_bw() +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 70, vjust = 1, hjust = 1))

## same varimp as dynamic models
ggplot(baselineVarimp, aes(x = tag, y = Importance, group = Variable)) +
  geom_line(aes(colour = Variable)) +
  facet_wrap(~mod) + 
  labs(title = '(B) Variable Importance - Static Models',
       x = 'Model',
       y = 'Importance',
       colour = 'Variable') +
  theme_bw() +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1))

## response curves baseline models
baselineRc <- rbind(allObsRc |> dplyr::select(-X) |> mutate(year = 'AO'), 
                    nativeRc |> dplyr::select(-X) |> mutate(year = 'NT'), 
                    invasiveRc |> dplyr::select(-X) |> mutate(year = 'IV')) |> 
  mutate(var = gsub('hursmean', 'hurs', var))
ggplot(baselineRc, aes(x = x, y = y)) +
  geom_line(aes(color = tag)) +
  facet_grid(mod~var, scales = 'free') +
  labs(title = '(A) Partial Dependence - Static Models',
       x = 'Variable',
       y = 'Prediction',
       color = 'Model') +
  guides(colour = guide_legend(ncol = 1)) +
  theme_bw() +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1))

## combine response curves and variable importance of all models in a single plot
## response curves
ggplot() +
  geom_line(data = yearlyRc, aes(x = x, y = y, color = factor(year))) +
  geom_line(data = baselineRc, aes(x = x, y = y, linetype = tag), color = 'black') +
  scale_linetype_manual(
    name = 'Static SDMs',
    values = c( 'AO' = 'solid', 'IV' = 'dotted', 'NT' = 'dashed')) +
  facet_grid(mod~var, scales = 'free') +
  labs(title = '(A) Partial Dependence Plots',
       x = 'Variable',
       y = 'Prediction',
       color = 'Dynamic SDMs') +
  guides(colour = guide_legend(ncol = 1)) +
  theme_bw() +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

## variable importance
years <- sort(unique(na.omit(yearlyVarimp$year))) |> as.character()
combinedVarimp <- bind_rows(baselineVarimp |> 
                              mutate(progression = as.factor(tag)), 
                            yearlyVarimp |> 
                              mutate(progression = as.factor(year)))
combinedVarimp <- combinedVarimp |> 
  mutate(progression = as.character(progression)) |> 
  mutate(progression = factor(progression,
                              levels = c('NT', years, 'IV', 'AO'))) |>
  arrange(progression)

ggplot(data = combinedVarimp, aes(x = progression, y = Importance)) +
  geom_line(aes(group = Variable, color = Variable)) +
  facet_wrap(~mod) +
  labs(title = '(B) Variable Importance Plot',
       x = 'Year',
       y = 'Importance',
       color = 'Variable') +
  theme_bw() +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


## use clusters from cluster analysis (script 11_) to create pdps
# 1. Create a uniform x-grid for each variable
allRcCombined <- rbind(yearlyRc, baselineRc)
x_grid_df <- allRcCombined |> 
  reframe(
    x = seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = 100),
    .by = var
  )

# 2. Interpolate y values onto the standardized x-grid
yearlyRc_interpolated <- yearlyRc |> 
  mutate(tag = case_when(
    year <= 2008 ~ 'cl1',
    year <= 2015 ~ 'cl2',
    year <= 2018 ~ 'cl3',
    TRUE         ~ 'cl4'
  )) |> 
  nest(data = -c(mod, year, tag, var)) |> 
  inner_join(x_grid_df, by = "var", relationship = "many-to-many") |> 
  mutate(
    y_interp = map2_dbl(data, x, ~ approx(.x$x, .x$y, xout = .y)$y)
  ) |> 
  filter(!is.na(y_interp))

# 3. Summarize across shared x_std points
yearlyRcCLUST <- yearlyRc_interpolated |> 
  summarise(
    n_obs  = n(),
    y_mean = mean(y_interp, na.rm = TRUE),
    y_sd   = sd(y_interp, na.rm = TRUE),
    .by = c(tag, var, x)
  )

ggplot(data = yearlyRcCLUST) +
  geom_ribbon(aes(x = x, 
                  ymin = y_mean-y_sd/2, ymax = y_mean+y_sd/2, 
                  fill = tag), alpha = 0.3) +
  geom_line(aes(x = x, y = y_mean, colour = tag)) +
  scale_fill_manual(values = c('cl4' = '#F8766D',
                               'cl3' = '#7CAE00',
                               'cl2' = '#00BFC4',
                               'cl1' = '#C77CFF'),
                    # labels = c('cl4' = '2019-2025',
                    #            'cl3' = '2016-2018',
                    #            'cl2' = '2009-2015',
                    #            'cl1' = '2002-2008')
                    ) +
  scale_colour_manual(values = c('cl4' = '#F8766D',
                                 'cl3' = '#7CAE00',
                                 'cl2' = '#00BFC4',
                                 'cl1' = '#C77CFF'),
                      # labels = c('cl4' = '2019-2025',
                      #            'cl3' = '2016-2018',
                      #            'cl2' = '2009-2015',
                      #            'cl1' = '2002-2008')
                      ) +
  facet_wrap(~var, scale = 'free') +
  labs(title = 'Partial Dependence Plots',
       x = 'Variable',
       y = 'Prediction',
       colour = 'Clusters',
       fill = 'Clusters') +
  theme_bw() +
  theme(strip.background = element_blank(),strip.text = element_text(hjust = 0),
        text = element_text(size = 20),
        legend.position = c(1-0.3/2,0.5/2))
