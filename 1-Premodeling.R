################################
# Install                      #
################################
# Codes for installing package from CRAN:


install.packages("devtools")
install.packages("tidyverse")
install.packages("readr")
install.packages("countrycode")
install.packages("rangeBuilder")
install.packages("sf")
install.packages("terra")
install.packages("rworldmap")
install.packages("maps")
install.packages("ridigbio")
install.packages("rgbif")
install.packages("BIEN")
install.packages("rinat")

# Codes for installing package from GitHub:
devtools::install_github("idiv-biodiversity/LCVP")
devtools::install_github("idiv-biodiversity/lcvplants")
install.packages("rnaturalearthhires", repos = "http://packages.ropensci.org", type = "source")
devtools::install_github("ropensci/rnaturalearthdata")
devtools::install_github("ropensci/rgnparser")
rgnparser::install_gnparser()
# In case of trouble you can install gnparser see the help at  https://github.com/gnames/gnparser#install-with-homebrew

devtools::install_github("brunobrr/bdc")
devtools::install_github("liibre/rocc")
devtools::install_github("sjevelazco/flexsdm")
devtools::install_github("andrefaa/ENMTML")




## %######################################################%##
#                                                          #
####                   3-Pre-modeling                   ####
#                                                          #
## %######################################################%##
{
  require(dplyr)
  require(terra)
  require(flexsdm)
  require(here)
  require(progress)
  require(ape)
}
getwd()

## %######################################################%##
#                                                          #
####             Create directory structure             ####
#                                                          #
## %######################################################%##
dir <- flexsdm::sdm_directory(
  main_dir = file.path(getwd(), "1-SDM"),
  calibration_area = TRUE,
  algorithm = c("gam", "gau", "gbm", "glm", "max", "net", "raf", "svm", "mean"), #mean es el ensamble for bees
  ensemble = NULL,
  threshold = FALSE,
  return_vector = TRUE
)
# dir[1] %>% fs::dir_tree(., recurse = TRUE)
dir %>% head()



# Occ database
unfilt_occ <- data.table::fread(here("occ_final_unfiltered.txt")) %>% 
  dplyr::tibble()
names(unfilt_occ)
unfilt_occ <- unfilt_occ %>% dplyr::select(id, species, x, y)


## %######################################################%##
#                                                          #
####                  Calibration area                  ####
#                                                          #
## %######################################################%##
here()

# Ecoregion
eco <- terra::vect("./Variables_US_to_Patagonia_5km/Ecoregions/Ecoregions.gpkg")

# Process species
sp <- unfilt_occ$species %>%
  unique() %>%
  sort()
i=1
for (i in 1:length(sp)) {
  x2 <-
    flexsdm::calib_area(
      data = unfilt_occ[unfilt_occ$species == sp[i], ],
      x = "x",
      y = "y",
      method = c("mask", eco, "ECO_ID")
    )
  x2$ECO_NAME <- NULL
  terra::writeVector(x2, 
                     file.path(dir[6], paste0(sp[i], ".gpkg")), overwrite = TRUE)
}
rm(unfilt_occ)


##%######################################################%##
#                                                          #
####               Evaluar la correlacion               ####
####          entre las variables ambientales           ####
#                                                          #
##%######################################################%##
require(ggplot2)
library(corrplot)

env_variables <- list.files("./Variav/", pattern = ".tif", full.names = TRUE) %>% 
  terra::rast() # una variable ambiental
names(env_variables)
corr <- correct_colinvar(env_layer = env_variables, method = c('pearson', th='0.7'))
mtrx <- corr$cor_table 
mtrx[mtrx<0.7] <- 0
corrplot(mtrx, type="upper", order="hclust")
  # ggsave(filename = "./Variablesusadas/correlationplot.png")

# filt escrever o vector com os nomes das variaveis ambientais selecionadas
filt <- c("silt_0-5cm", "hurs_mean", "bio17", "bio18",  "cmi_mean", "bio10", "bio13", "bio14", "bio8", "bio6", "bio11", "bio4", "bio3", "bio7", "bio19", "bio16", "bio1")
filt_m <- !colnames(mtrx)%in%filt
corrplot(mtrx[filt_m,filt_m], type="upper", order="hclust")
env_variables <- env_variables[[!names(env_variables)%in%filt]]
names(env_variables)

dir.create("./Variaveisusadas_filtradas")
# terra::writeRaster(env_variables, file.path("./Variaveisusadas_filtradas", paste0(names(env_variables), ".tif")))


##%######################################################%##
#                                                          #
####            Filtering species occurrence            ####
#                                                          #
##%######################################################%##
# Registros de especies no filtradas
db0 <- data.table::fread(here("occ_final_unfiltered.txt")) %>% 
  dplyr::tibble()
db0 <- db0 %>% dplyr::select(id = id, species = species, x, y)

nocc <- db0 %>% dplyr::group_by(species) %>% 
  dplyr::count() %>% 
  dplyr::arrange(species)
nocc
nocc <- nocc %>% mutate(need_bias_corr = FALSE)
nocc[nocc$n > 50, 3] <- TRUE

sp <- nocc %>% dplyr::filter(need_bias_corr==T) %>% pull(species)

env <- list.files("./Variaveisusadas_filtradas/", 
                  pattern = ".tif", full.names = TRUE) %>% 
  terra::rast() # una variable ambientales
# Eliminamos aquellas especificas de especies fosoriales y aquaticas
env <- env[[!names(env)%in%c("cti", "clay_0-5cm", "sand_0-5cm")]]
names(env)

list_filt_occ <- list()
for(i in 1:length(sp)){
  message(i/length(sp))
  bin_4 <- 
    flexsdm::occfilt_env(
      data = db0[db0$species==sp[i],],
      x = "x",
      y = "y",
      id = "id",
      env_layer = env,
      nbins = 4
    )
  bin_6 <- 
    flexsdm::occfilt_env(
      data = db0[db0$species==sp[i],],
      x = "x",
      y = "y",
      id = "id",
      env_layer = env,
      nbins = 6
    )
  bin_8 <- 
    flexsdm::occfilt_env(
      data = db0[db0$species==sp[i],],
      x = "x",
      y = "y",
      id = "id",
      env_layer = env,
      nbins = 8
    )
  bin_10 <- 
    flexsdm::occfilt_env(
      data = db0[db0$species==sp[i],],
      x = "x",
      y = "y",
      id = "id",
      env_layer = env,
      nbins = 10
    )
  
  # Calculo de autocorrelacion espacial
  list_bin <- list(bin_4, bin_6, bin_8, bin_10)
  names(list_bin) <- c('bin_4', 'bin_6', 'bin_8', 'bin_10')
  imoran <- list()
  for(ii in 1:4){
    coord <- list_bin[[ii]] %>% dplyr::select(x, y)
    data <- data.frame(terra::extract(env, coord))[-1]
    distm <- dist(coord)
    distm <- as.matrix(distm)
    distm <- 1/distm
    diag(distm) <- 0
    try(imoran[[ii]] <-
          apply(data, 2, function(x)
            ape::Moran.I(x, distm, na.rm = T)[c(1, 4)] %>% unlist) %>% data.frame() %>% as_tibble()
    )
    try(imoran[[ii]] <- imoran[[ii]][1,])
    try(imoran[[ii]]$mean_bin <- apply(imoran[[ii]], 1, mean))
    imoran[[ii]]$nrecords <- nrow(coord)
  }
  names(imoran) <- c('bin_4', 'bin_6', 'bin_8', 'bin_10')
  imoran <- bind_rows(imoran, .id="nbins")
  # Seleccion 
  bin_selected <- imoran %>%
    filter(mean_bin <= mean(mean_bin)) %>%
    filter(nrecords == max(nrecords)) %>%
    pull("nbins")
  bin_selected <- bin_selected[1]
  # imoran <- data.frame(species=sp[i], imoran)
  list_filt_occ[[i]] <- list_bin[[bin_selected]] %>% mutate(species=sp[i])
}

# dind data.frames
list_filt_occ_final <- bind_rows(list_filt_occ)
# merge unfiltered species with filtered one
db1 <- db0 %>% filter(!species%in%sp)
# remove occ in a same cell
filtr <- env$elevation
filtr[!is.na(filtr)] <- as.data.frame(filtr, cells = T)[,"cell"]
names(filtr) <- "cell"
flexsdm::sdm_extract(db1)

db1$cell <- terra::extract(filtr, db1 %>% select(x, y))[,"cell"]
db1 <- db1 %>% 
  group_by(species) %>% 
  filter(!is.na(cell)&!duplicated(cell)) %>% 
  group_by() %>% 
  select(-cell)

filt_occ <- bind_rows(db1, list_filt_occ_final)
# readr::write_tsv(filt_occ, "occurrences_cleaned_final_FILTERED.txt") RESULTADO FINAL, txt con la lista completa de <50 no filtrada y >50 filtrada. 

# seleccionar especies con >=5 registros 
occ <- readr::read_tsv(here("occurrences_cleaned_final_FILTERED.txt"))
filt <- occ %>% group_by(species) %>% count
filt <- filt %>% filter(n>=5) %>% pull(species)
occ_mayor5 <- occ %>% filter(species%in%filt) #filtra mayor que 5 
occ_menor5 <- occ %>% filter(!species%in%filt) #filtra menor que 5 
# readr::write_tsv(occ_mayor5, here("occurrences_cleaned_final_FILTERED_mayore5.txt"))
# readr::write_tsv(occ_menor5, here("occurrences_cleaned_final_FILTERED_menores5.txt"))
#Generamos dos listas, porque las sp que tienen menos de 5 registros u ocurrencias van a ser modeladas diferente. 

## %######################################################%##
#                                                          #
####      Psd-abs, background points and partition      ####
####       for species without spatial partition        ####
#                                                          #
## %######################################################%##
occ <- readr::read_tsv(here("occurrences_cleaned_final_FILTERED_mayore5.txt"))
# occ$species %>% table %>% sort

occ$pr_ab <- 1 

# Lista con las areas de calibracion 
clibarea <- list.files(
  "./1-SDM/1_Inputs/3_Calibration_area/",
  pattern = "gpkg$",
  full.names = TRUE
)
names(clibarea) <- clibarea %>% basename() %>% gsub(".gpkg$", "", .)
clibarea["Xenodon severus"]
grep("Urotheca", clibarea, value = TRUE) # En esta parte indicamos la búsqueda en la base de datos las coincidencias con el nombre que damos ahí...

pred <- "./Variablesusadas_filtradas/" %>% list.files(full.names = T)
pred <- rast(pred)
pred <- pred[[c("bio1", "elevation")]]
plot(pred)
pred <- homogenize_na(pred) #Aquí homogenizamos dando NA en la parte de los raster de TODAS las variables donde no posiblemente no hay ningun dato
alayer <- pred[[1]]

db_bg <- db_pa <- list() # object to store backgroud points (db_bg) and pseudo-absences (db_pa). Aquí haremos la partición de los datos para pseudoausencias y presencias, esto va determinado por el número de ocurrencias...
i=1
sp <- occ$species %>% unique() %>% sort()
for (i in 1:length(sp)) {
  message("species ", i)
  coord <- occ %>%
    dplyr::filter(species == sp[i]) %>%
    dplyr::select(x, y, pr_ab)
  v <- terra::vect(clibarea[[sp[i]]])
  r <- alayer %>%
    crop(., v) %>%
    mask(., v)
  
  db_pa[[i]] <-
    set.seed(15) %>% #función set.seed deja un valor fijo en esta selección de datos
    flexsdm::sample_pseudoabs(
      data = coord,
      x = "x",
      y = "y",
      n = nrow(coord) * 2, # Doble de ps-ausencias que presencias (para ANN, SVM, RF)
      method = c("geo_const", width = "50000"), # restriccion espacial 
      rlayer = r
    ) %>% 
    dplyr::bind_rows(coord, .) 
  #Ahora para los datos que son mayores o iguales a 15 kfold "sencillo" 
      if(nrow(coord)>=15){
        db_pa[[i]] <- db_pa[[i]] %>% 
        flexsdm::part_random( # particion k-fold de los datos
          data = .,
          pr_ab = "pr_ab",
          method = c(
            method = "kfold",
            folds = 5
          )
        )
      } else {
        db_pa[[i]] <- db_pa[[i]] %>%   #Pero aquí se realiza el kfold con 5 réplicas para las sp que tienen <15 ya que si no es así, al haber tan pocos datos puede haber sobreajuste del modelo 
        flexsdm::part_random( # particion k-fold de los datos
          data = .,
          pr_ab = "pr_ab",
          method = c(method = "rep_kfold",
                     folds = 5,
                     replicates = 5)
        )
      }
    
  # table(db_pa[[i]]$pr_ab, db_pa[[i]]$.part)

  db_bg[[i]] <-
    flexsdm::sample_background(
      data = coord,
      x = "x",
      y = "y",
      n = 10000,
      method = "random",
      rlayer = r,
    ) 
  if(nrow(coord)>=15){
    db_bg[[i]] <- db_bg[[i]] %>% 
      flexsdm::part_random( # particion k-fold de los datos
        data = .,
        pr_ab = "pr_ab",
        method = c(
          method = "kfold",
          folds = 5
        )
      )
  } else {
    db_bg[[i]] <- db_bg[[i]] %>% 
      flexsdm::part_random( # particion k-fold de los datos
        data = .,
        pr_ab = "pr_ab",
        method = c(method = "rep_kfold",
                   folds = 5,
                   replicates = 5)
      )
  }
}

names(db_pa) <- names(db_bg) <- sp
db_pa <- bind_rows(db_pa, .id = "species")
db_bg <- bind_rows(db_bg, .id = "species")
db_pa <- db_pa %>% arrange(species, desc(pr_ab))

data.table::fwrite(db_pa, file.path("./1-SDM/1_Inputs/1_Occurrences", "1_occ_presabs_randompart.gz")) # gz version comprimida de txt 
data.table::fwrite(db_bg, file.path("1-SDM/1_Inputs/1_Occurrences", "1_occ_bkground_randompart.gz"))

#Finalizamos con dos megadataframes (comprimidos de txt)en formato gz, una de las presencias y pseudoausencias donde se usó el método normal de kfold y la otra de background points 