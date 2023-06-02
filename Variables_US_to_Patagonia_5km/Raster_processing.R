require(terra)
require(dplyr)
setwd("D:/Variables")

"D:/Variables/Variables_US_to_Patagonia_5km"

# read extent
p <- terra::vect("./ExtentAmerica2.gpkg")
plot(p)
# CHELSA 5 km
l <- list.files("./CHELSA_5km/1981_2010/", full.names = TRUE)
i=1
for(i in 1:length(l)){
  r <- terra::rast(l[i])
  r <- r %>% terra::crop(., p) %>% terra::mask(., p)
  terra::writeRaster(r, file.path("D:/Variables/Variables_US_to_Patagonia_5km/CHELSA-1981_2010", basename(l[i])))
}

# CHELSA 5 km
l <- list.files("./CHELSA_5km/1981_2010/", full.names = TRUE)
i=1
for(i in 1:length(l)){
  r <- terra::rast(l[i])
  r <- r %>% terra::crop(., p) %>% terra::mask(., p)
  terra::writeRaster(r, file.path("D:/Variables/Variables_US_to_Patagonia_5km/CHELSA-1981_2010", basename(l[i])))
}

# Elevation
dd <- "D:/Variables/Variables_US_to_Patagonia_5km/Elevation"
dir.create(dd)
r <- "./Elevation/wc2.1_30s_elev1.tif" %>% terra::rast()
base_r <- "./Variables_US_to_Patagonia_5km/CHELSA-1981_2010/ai.tif" %>% terra::rast()
r <- crop(r, base_r)
r <- resample(r, base_r)
r <- mask(r, base_r)
names(r) <- "elevation"
terra::writeRaster(r, file.path(dd, "elevation.tif"))

# Geomorpho
dd <- "D:/Variables/Variables_US_to_Patagonia_5km/Geomorpho"
base_r <- "./Variables_US_to_Patagonia_5km/CHELSA-1981_2010/ai.tif" %>% terra::rast()
l <- list.files("./Geomorpho", full.names = TRUE, pattern = ".tif$")
i=1
for(i in 1:length(l)){
  r <- terra::rast(l[i])
  r <- crop(r, base_r)
  r <- resample(r, base_r)
  r <- mask(r, base_r)
  names(r) <- strsplit(basename(l[i]), "_")[[1]][2]
  terra::writeRaster(r, file.path(dd, paste0(names(r), ".tif")))
}
dir.create(dd)
list.files()

# SoilGrids V2
dd <- "D:/Variables/Variables_US_to_Patagonia_5km/SoilGrids_v2"
dir.create(dd)
base_r <- "./Variables_US_to_Patagonia_5km/CHELSA-1981_2010/ai.tif" %>% terra::rast()
l <- list.files("./SoilGridsV2/1km/", full.names = TRUE, pattern = ".tif$")
l <- l[-1]
i=1
# Variables with problem? bdod (bulk density)
rm(r2)
for(i in 1:length(l)){
  r <- terra::rast(l[i])
  r <- terra::project(r, crs(base_r))
  r <- crop(r, base_r)
  r <- aggregate(r, fact=4, fun="mean", na.rm=T)
  r <- resample(r, base_r)
  r <- mask(r, base_r)
  names(r) <- gsub(".tif", "", gsub("_mean_1000", "", basename(l[i])))
  terra::writeRaster(r, file.path(dd, paste0(names(r), ".tif")))
}
dir.create(dd)
list.files()


# Crop ecoregions
base_r <- "./Variables_US_to_Patagonia_5km/CHELSA-1981_2010/ai.tif" %>% terra::rast()
p <- "Terrestrial Ecoregions of the World/Ecoregions2017.shp" %>% terra::vect()
p <- terra::crop(p, base_r)
plot(p)

dd <- "./Variables_US_to_Patagonia_5km/Ecoregions"
dir.create(dd)
p <- p[,names(p)%in%c("ECO_NAME", "ECO_ID")]
terra::writeVector(p, file.path(dd, "Ecoregions.gpkg"))
