## %######################################################%##
#                                                          #
####  Estimate environmental distance for species with  ####
####                      <5 occurrences                ####
#                                                          #
## %######################################################%##
# devtools::install_github("sjevelazco/flexsdm")
{
  require(dplyr)
  require(terra)
  require(flexsdm)
  require(here)
  require(progress)
  require(raster)
  require(dismo)
  require(ggplot2)
}


## %######################################################%##
#                                                          #
####             Read occurrences databases             ####
#                                                          #
## %######################################################%##
occ <- data.table::fread("occurrences_cleaned_final_FILTERED_menores5.txt") %>% tibble()

# Environmental variables - Current conditions
env <-
  "./Variablesusadas_filtradas/" %>%
  list.files(., full.names = TRUE) %>%
  terra::rast()
env %>% names()
names(env) <- gsub("_0-5cm", "", names(env))
env_names <- names(env)

# Extract env conditions
occ <- sdm_extract(occ, x = "x", y = "y", env_layer = env, filter_na = TRUE)

# Count number of presences
n_occ <- occ %>%
  pull(species) %>%
  table() %>%
  sort()

# Colombia boundaries
# col <- rnaturalearth::ne_countries(country = "colombia", scale = "large")
# col <- terra::vect(col)
# terra::writeVector(col, "Colombia.gpkg")
col <- terra::vect("Colombia.gpkg")

# Recorte para el area de proyeccion
env_projection <- terra::crop(env, terra::vect("./area_projection.gpkg"))
env_projection <- env_projection %>%
  terra::crop(., col) %>%
  terra::mask(., col)
env_projection <- raster::stack(env_projection)


## %######################################################%##
#                                                          #
####        Loop for modeling species with >15          ####
#                                                          #
## %######################################################%##
sp <- names(n_occ[n_occ > 1])


# Extraccion de los nombres de las variables
env_names <- names(env)
names(env_names)[env_names %in% c("bio1", "bio13", "bio15", "bio2", "elevation")] <- "climatica"
names(env_names)[env_names %in% c("clay", "sand")] <- "suelo"
names(env_names)[env_names %in% "cti"] <- "hidrologica"


db_env_used <- readxl::read_excel("./Lista limpiaFinal.xlsx", sheet = 1)
db_env_used <- db_env_used %>% dplyr::select(species = accepted_name_bin, climatica, suelo, hidrologica)

# Loop para procesar las especies
i <- 21
for (i in 1:length(sp)) {
  message(paste("Estimating sp", i, sp[i]))
  pa <- occ %>% dplyr::filter(species == sp[i])

  env_used <- db_env_used %>%
    dplyr::filter(species == sp[i]) %>%
    dplyr::select(-species) %>%
    data.frame() %>%
    unlist() %>%
    na.omit() %>%
    names()

  env_used <- env_names[names(env_names) %in% env_used]
  names(env_used) <- NULL

  #### domain - Gower distance ####

  prd <- dismo::domain(pa %>% dplyr::select({
    env_used
  }) %>% as.matrix()) 
  prd <- prd %>% 
    raster::predict(., env_projection[[names(prd@presence)]]) %>%
    terra::rast()
  buf <- terra::buffer(
    terra::vect(pa[c("x", "y")],
      geom = c("x", "y"), crs = crs(prd)
    ),
    width = 50000
  )
  buf <- (terra::rasterizeGeom(buf, prd) > 0) %>% terra::mask(., prd)
  prd[!buf] <- 0
  prd %>% plot()
  prd <- c(prd, prd)
  terra::writeRaster(prd,
    here(
      "1-SDM/2_Outputs/1_Current/FinalModels",
      paste0(sp[i], ".tif")
    ),
    overwrite = TRUE
  )
}


