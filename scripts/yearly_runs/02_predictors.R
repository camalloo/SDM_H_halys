########################################################
############ Download and select predictors ############
########################################################

library(tidyverse)
library(terra)

folder <- 'baseline_models_yearly_runs' # name of the mother folder

######### World map ----------
world <- geodata::world(path = 'data/geodata/', resolution = 1)
world <- world[!world$GID_0 == 'ATA'] # no Antarctica 

######### Bioclim+ ----------
# dir.create('data/CHELSAbioclim/raw')
# ## download vars
# links <- read.table('data/CHELSAbioclim/bioclim_chelsa.txt', header = FALSE)[,1]
# for (link in links) {
#   name <- basename(link)
#   if (grepl('pet|cmi|bio10|bio11|bio16|bio17', name)) next # unsure about these
#   download.file(link, destfile = file.path('data/CHELSAbioclim/raw', name))
# }
# 
# ## offload memory
# rm(link, links, name)
# gc()

## load vars
bio <- rast(list.files(paste0('data/CHELSAbioclim/raw'), full.names = T))
nlyr(bio) # check if all loaded
## handle names and assign value of 0 to NAs in gdd10
names(bio) <- gsub('CHELSA_|_1981.*', '', names(bio))
names(bio)
## aggregate to 10km
bio <- aggregate(bio, fact = 10, fun = 'mean', na.rm = T, cores = 10)
res(bio)
## mask the world
bio <- crop(bio, world, snap = 'out', mask = T)

## plot them
plot(bio[[1:ceiling(nlyr(bio)/2)]], legend = F)
plot(bio[[(ceiling(nlyr(bio)/2)+1):nlyr(bio)]], legend = F)
######### 


######### Topology ----------
# dem <- geodata::elevation_global(path = 'data/geodata/', res = 5)
# names(dem) <- 'dem'
# dem <- crop(dem, world, snap = 'out', mask = T)
# plot(dem)
######### 


######### stack all vars together ----------
# preds <- c(bio, dem)
preds <- bio
nlyr(preds)
# ## scale values to z-score
# stats <- global(preds, c('mean', 'sd'), na.rm = T)
# preds <- (preds-stats$mean) / stats$sd
## plot all predictors
plot(preds[[1:ceiling(nlyr(preds)/2)]], legend = F)
plot(preds[[(ceiling(nlyr(preds)/2)+1):nlyr(preds)]], legend = F)
## offload memory
rm(bio, dem)
gc()

## save to disk
writeRaster(preds, filename = file.path('output', folder, 'modellingData/predictors.tif'),
            filetype = 'GTiff', gdal = c('COMPRESS=DEFLATE'), overwrite = T)
######### 

