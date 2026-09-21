library("terra")


################# Predictors_________________________________
######## Download chelsa bioclims+ --------
# folder <- "data/CHELSAbioclim"
# dir.create(file.path(folder, 'raw'))
# 
# origin <- 'https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio/CHELSA_'
# tail <- '_1981-2010_V.2.1.tif'
# 
# variables <- c('cmi_range', 'cmi_mean', 'fcf', 'gdd10', 'gddlgd10', 'gdgfgd10',
#                'hurs_mean', 'hurs_range', 'npp', 'pet_penman_mean', 'pet_penman_range',
#                'sfcWind_mean', 'sfcWind_range', 'vpd_mean', 'vpd_range', 'swb')
# 
# for (v in variables) {
#   thisVar <- paste0(origin, v, tail)
#   download.file(thisVar, destfile = file.path(folder, 'raw', basename(thisVar)))
# }
# ## download radiation (not harmonized names)
# thisVar <- "https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio/CHELSA_rsds_1981-2010_mean_V.2.1.tif"
# download.file(thisVar, destfile = file.path("data/CHELSAbioclim/raw/CHELSA_rsds_mean_1981-2010_V.2.1.tif"))
# thisVar <- "https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio/CHELSA_rsds_1981-2010_range_V.2.1.tif"
# download.file(thisVar, destfile = file.path("data/CHELSAbioclim/raw/CHELSA_rsds_range_1981-2010_V.2.1.tif"))

# ## download 19 bioclimatic variables from chelsa
# links <- read.table("data/CHELSAbioclim/CHELSA_bioclim_links_present.txt", header = FALSE)[,1]
# for (link in links) {
#   name <- basename(link)
#   download.file(link, destfile = file.path(folder, 'raw', name))
# }

######## Load CHELSA bioclims 1981-2010 (1km) --------
bio <- rast(list.files(paste0("data/CHELSAbioclim/", folder, "/"), full.names = T))
nlyr(bio) # check if all were loaded

plot(bio[[1:16]]) # all
plot(bio[[17:32]]) # all
plot(bio[[33:37]]) # all

plot(bio[[1]]) # one and add world map
plot(world, lw = 0.5, add = T)

## crop them to include only land
world <- world[!world$GID_0 == 'ATA'] # world map no Antarctica
plot(world)
bio <- crop(bio, world, mask = T)
## check with plots if all worked
plot(bio[[1:16]])
plot(bio[[17:31]])

## dem from CHELSA
dem <- rast('data/otherLayers/dem_latlong.nc')
dem <- crop(dem, world, mask = T)
plot(dem)

## land cover/use from copernicus land monitoring service
builtup <- rast('data/otherLayers/Builtup_cover.tif')
builtup <- aggregate(builtup, fact = 10, fun = 'mean', na.rm = T)
builtup <- crop(builtup, world , mask = T)
builtup <- project(builtup, dem, method = 'bilinear', threads = T)
plot(builtup)

## group all to a single raster
bio <- c(bio, dem, builtup)
## save all the bios in the folder
writeRaster(bio, filename = paste0('data/CHELSAbioclim/historicCropped/bioAll.tif'), gdal = c('COMPRESS=DEFLATE'), filetype = 'GTiff', overwrite = T)





################# Create output folders_________________________________
folder <- 'AllPredictors' # name of the mother folder

if (!dir.exists(paste0('output/', folder, '/'))) {
  dir.create(paste0('output/', folder, '/'))
  dir.create(paste0('output/', folder, '/models/'))
  dir.create(paste0('output/', folder, '/models/crossValidation/'))
  dir.create(paste0('output/', folder, '/models/evaluation/'))
  dir.create(paste0('output/', folder, '/models/rast/'))
  dir.create(paste0('output/', folder, '/modellingData/'))
  dir.create(paste0('output/', folder, '/modellingData/ModelRegion/'))
  dir.create(paste0('output/', folder, '/projections/'))
  dir.create(paste0('output/', folder, '/projections/rast'))
} else {
  print('output dirs already exist')
}


