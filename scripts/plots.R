###################################################
############ Plot and generate figures ############
###################################################

library(terra)
library(tidyverse)

folder <- 'baseline_models_yearly_runs'

######### Load data ----------
world <- geodata::world(path = 'data/geodata/', resolution = 1)
world <- world[!world$GID_0 == 'ATA'] # world map no Antactica
region <- vect(file.path('output', folder, 'modellingData/ModelRegion/modellingRegionPoly.shp')) # modelling region polygon

data <- read.csv(file.path('output', folder, 'modellingData/modData.csv'))
dataSPAT <- vect(data, geom = c('x', 'y'), crs = 'EPSG:4326')
## number of presence pixels
nrow(filter(data, data$presence == 1))

## folded data
folds <- read.csv(file.path('output', folder, 'modellingData/FoldedData.csv'))
foldsSPAT <- vect(folds, geom = c('x', 'y'), crs = 'EPSG:4326')
## check the number of test data for each fold
(n <- table(folds[, c('fold', 'presence')]))
(nt <- rowSums(n))  # n total rows per fold
(np <- n[ , 2])  # n presence rows per fold

## discrimination metrics
discrimination <- c('AUC', 'TSS', 'Kappa', 'Spearman')
## calibration metrics
calibration <- c('Miller', 'Brier', 'Boyce')
## folds evaluation
evalResults <- read.csv(file.path('output', folder, 'models/crossValidation/foldsEval_ALL.csv'))
## prepare df
ev <- evalResults |> 
  pivot_longer(cols = -fold, names_to = 'mod', values_to = 'score') |> 
  arrange(mod) |> 
  separate(mod, into = c('mod', 'metric'), sep = '_') |> 
  mutate(thresh = thresh[metric],
         normScore = pmax( # force the lower limit to -1
           case_when(
             ## handle Miller who is an intorno of 1
             metric == 'Miller' ~ 1 - abs(score-1)/thresh,
             ## handle Brier where 0 is best result, so pivot it as 1 is best
             metric == 'Brier' & score <= thresh ~ ((1-score) - (1-thresh)) / (1 - (1-thresh)),
             metric == 'Brier' & score > thresh ~ ((1-score) - (1-thresh)) / (1-thresh),
             ## other cases
             score >= thresh ~ (score - thresh) / (1 - thresh),
             score < thresh ~ (score - thresh) / thresh),
           -1))
#########


######### Make plots ----------
## presences and study area
plot(region, col = 'gold', border = NA, main = 'Study area')
plot(world, add = T)
points(dataSPAT[dataSPAT$presence == 1], col = 'red', cex = 0.3)
legend('bottomleft', cex = 1.4, xpd = T, bty = 'n', inset = c(0, 0.15), x.intersp = 0.3,
       legend = 'presence', pch = 16, col = 'red')


## blocks and cross validation
plot(foldsSPAT, 'fold', cex = 0.4, col = hcl.colors(k, 'viridis'),main = 'Block CV folds', 
     legend = T, plg = list(x = 170, y = 90, xpd = T, 
                            x.intersp = 0.1, y.intersp = 0.3, pt.cex = 2, cex = 2))
plot(world, add = T, lw = .4)
plot(vect(blocks$blocks), border = "brown", add = TRUE, lw = .7)


## discrimination metrics
ggplot(data = ev |> 
         filter(metric %in% discrimination),
       aes(x = metric, y = normScore, col = mod)) +
  geom_point(position = position_dodge(width = 1), size = 2) +
  ylim(-1, 1.1) + 
  stat_summary(aes(group = mod), fun = mean, 
               geom = 'point', size = 2, shape = 4, color = 'black',
               position = position_dodge(width = 1)) +
  geom_vline(xintercept = seq(0.5, length(discrimination) + 0.5, by = 1),
             color = "black", linewidth = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = 'Discrimination metrics') +
  theme(axis.text.x = element_text(size = 10), axis.title.x = element_blank())
## calibration metrics
ggplot(data = ev |> 
         filter(metric %in% calibration),
       aes(x = metric, y = normScore, col = mod)) +
  geom_point(position = position_dodge(width = 1), size = 2) +
  ylim(-1, 1.1) + 
  stat_summary(aes(group = mod), fun = mean, 
               geom = 'point', size = 2, shape = 4, color = 'black',
               position = position_dodge(width = 1)) +
  geom_vline(xintercept = seq(0.5, length(calibration) + 0.5, by = 1),
             color = "black", linewidth = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = 'Calibration metrics') +
  theme(axis.text.x = element_text(size = 10), axis.title.x = element_blank())

