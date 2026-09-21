########################################################
############ Download and select predictors ############
########################################################

library(tidyverse)
library(terra)

folder <- 'climate_human' # name of the mother folder

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
## remove some rasters
bio <- bio[[!names(bio) %in% c('gdd10', 'gsl', 'sfcWindmean', 'rsdsmean')]]
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
dem <- geodata::elevation_global(path = 'data/geodata/', res = 5)
names(dem) <- 'dem'
dem <- crop(dem, world, snap = 'out', mask = T)
plot(dem)
######### 


######### Human ----------
## Calculate distances from human features in Mollweide crs for faster distance calculations
moll <- '+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'
## template raster
wr <- project(world, moll)
r <- rast(ext(wr), resolution = 10000, crs = moll) # 10000m resolution
plot(wr)
## rasterize the world and invert to calculate only on land
# wr <- project(world, moll)
# plot(wr)
# wr <- rasterize(wr, r, field = 1)
# wr[is.na(wr)] <- 0
# wr <- mask(wr, project(world, moll), inverse = T)
## plot to check
# plot(wr)

## Airports distance
airports <- read.csv('data/otherLayers/airport_volume_airport_locations.csv')
airports <- airports |> 
  select(orig = Orig, name = Name, country = Country.Name, lat = Airport1Latitude, long = Airport1Longitude) |> 
  filter(!is.na(long), !is.na(lat))
## spat vector
airportsSPAT <- vect(airports, geom = c('long', 'lat'), keepgeom = T, crs = 'EPSG:4326')
airportsSPAT <- project(airportsSPAT, moll)
plot(wr)
points(airportsSPAT, col = 'red')
airportsSPAT <- rasterize(airportsSPAT, r)
plot(airportsSPAT, col = 'red')
## compute distance
air <- distance(airportsSPAT) #, filename = 'data/otherLayers/airports_calculated_Distance.tif', overwrite = T)
names(air) <- 'air'
plot(air)

## Ports distance
ports <- vect('data/otherLayers/attributed_ports.geojson')
ports <- project(ports, moll)
plot(wr)
points(ports, col = 'red')
ports <- rasterize(ports, r)
plot(ports, col = 'red')
## compute distance
prt <- distance(ports) #, filename = 'data/otherLayers/ports_calculated_Distance.tif', overwrite = T)
names(prt) <- 'prt'
plot(prt)

## Roads distance
roads <- vect('data/otherLayers/gROADS_v1.gdb')
roads <- roads[roads$FCLASS <= 2] # keep highway and primary roads
roads <- project(roads, moll)
# plot(wr)
# lines(roads, col = 'red')
roads <- rasterize(roads, r, field = 1, touches = T)
plot(roads)
## compute distance
rds <- distance(roads) #, filename = 'data/otherLayers/gROADS_calculated_Distance.tif', overwrite = T)
names(rds) <- 'rds'
plot(rds)

## offload memory
rm(wr, r, moll,
   airports, airportsSPAT, 
   roads, 
   ports)
gc()

## BuiltUp cover fraction
blt <- rast('data/otherLayers/BuiltUp_CoverFraction_layer.tif')
names(blt) <- 'blt'
plot(blt)

## Roads density
rdt <- rast('data/otherLayers/grip4_total_dens_m_km2.asc')
names(rdt) <- 'rdt'
rdt <- subst(rdt, NA, 0) # transform NAs to 0
plot(rdt, col = c('black', 'darkgreen', 'green', 'yellow', 'orange', 'red'), breaks = c(-Inf, 0, 10, 100, 250, 1000, Inf))
######### 


######### stack all vars together ----------
air_proj <- project(air, bio, threads = T)
air_proj
blt_proj <- project(blt, bio, threads = T)
blt_proj
prt_proj <- project(prt, bio, threads = T)
prt_proj
rds_proj <- project(rds, bio, threads = T)
rds_proj
rdt_proj <- project(rdt, bio, threads = T)
rdt_proj

## stack all layers together to have the predictors
preds <- c(bio, dem, air_proj, blt_proj, prt_proj, rds_proj, rdt_proj)
nlyr(preds)
preds <- crop(preds, world, snap = 'out', mask = T)

## plot all predictors
plot(preds[[1:ceiling(nlyr(preds)/2)]], legend = F)
plot(preds[[(ceiling(nlyr(preds)/2)+1):nlyr(preds)]], legend = F)
## offload memory
rm(stats, air, bio, blt, dem, prt, rds, rdt,
   air_proj, blt_proj, prt_proj, rds_proj, rdt_proj)
gc()

## save to disk
writeRaster(preds, filename = file.path('output', folder, 'modellingData/predictors.tif'),
            filetype = 'GTiff', gdal = c('COMPRESS=DEFLATE'), overwrite = T)
######### 

