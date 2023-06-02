## %######################################################%##
#                                                          #
####                     4-Modeling                     ####
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
  require(ggplot2)
}
memory.limit(1000000)


##%######################################################%##
#                                                          #
####             Read occurrences databases             ####
#                                                          #
##%######################################################%##
occ <- data.table::fread(here("1-SDM/1_Inputs/1_Occurrences", "1_occ_presabs_randompart.gz")) %>% tibble()
bkg <- data.table::fread(here("1-SDM/1_Inputs/1_Occurrences", "1_occ_bkground_randompart.gz")) %>% tibble()

# Environmental variables - Current conditions
env <-
  "./Variablesusadas_filtradas/" %>% 
  list.files(., full.names = TRUE) %>%
  terra::rast()
env %>% names
names(env) <- gsub("_0-5cm", "", names(env))
env_names <- names(env)

# Extract env conditions
occ <- sdm_extract(occ, x = "x", y = "y", env_layer = env, filter_na = TRUE)
bkg <- sdm_extract(bkg, x = "x", y = "y", env_layer = env, filter_na = TRUE)

# Count number of presences
n_occ <- occ %>% dplyr::filter(pr_ab==1) %>% pull(species) %>% table() %>% sort()


##%######################################################%##
#                                                          #
####        Loop for modeling species with >15          ####
#                                                          #
##%######################################################%##
perf_dir <- here("1-SDM/2_Outputs/0_Model_performance")
sp <- names(n_occ[n_occ<15])

# Recorte para el area de proyeccion 
env_projection <- terra::crop(env,  terra::vect("./area_projection.gpkg"))

# Extraccion de los nombres de las variables 
env_names <- names(env)
names(env_names)[env_names%in%c("bio1", "bio13", "bio15", "bio2", "elevation")] <- "climatica"
names(env_names)[env_names%in%c("clay", "sand")] <- "suelo"
names(env_names)[env_names%in%"cti"] <- "hidrologica"


db_env_used <- readxl::read_excel("./Lista limpiaFinal.xlsx", sheet=1)
db_env_used <- db_env_used %>% dplyr::select(species=accepted_name_bin, climatica, suelo, hidrologica )


# Funcion para reemplazar las ausencias por background points
replace_absences <- function(x, bg){
  result <- dplyr::bind_rows(x %>% dplyr::filter(pr_ab==1), bg)
  return(result)
}

# Loop para procesar las especies
i=1
for (i in 1:length(sp)) {
  message(paste("Modeling sp", i, sp[i]))
  pa <- occ %>% dplyr::filter(species == sp[i])
  b <- bkg %>% dplyr::filter(species == sp[i])
  pa <- pa %>% dplyr::select(-".part")
  b <- b %>% dplyr::select(-".part")
  
  env_used <- db_env_used %>% 
    dplyr::filter(species== sp[i]) %>% 
    dplyr::select(-species) %>% 
    data.frame() %>% 
    unlist() %>% 
    na.omit() %>% 
    names()
  
  env_not_used <- env_names[!names(env_names)%in%env_used]
  pa <- pa %>% dplyr::select(-{env_not_used})
  b <- b %>% dplyr::select(-{env_not_used})
  
  env_used <- env_names[names(env_names)%in%env_used]
  names(env_used) <- NULL
  
  #### Boosted regression trees ####
  m_gbm <- flexsdm::esm_gbm(
    data = pa,
    response = "pr_ab",
    predictors = env_used,
    partition = ".part",
    thr = "max_sorensen"
  )
  
  #### Maximum entropy ####
  try(m_max <- flexsdm::esm_max(
    data = pa,
    response = "pr_ab",
    predictors = env_used,
    background = b,
    partition = ".part",
    thr = "max_sorensen"
  ))
  
  
  #### Neural Network ####
  m_net <- esm_net(
    data = pa,
    response = "pr_ab",
    predictors = env_used,
    partition = ".part",
    thr = "max_sorensen"
  )
  
  #### Random forest ####
  m_raf <- tune_raf(
    data = pa,
    response = "pr_ab",
    predictors = env_used,
    partition = ".part",
    grid = expand.grid(mtry = seq(1, length(env_names), 1)),
    thr = "max_sorensen",
    metric = "SORENSEN",
    n_cores = 5
  )
  
  if(length(m_raf)>1){
    h <- m_raf$hyper_performance
    h %>% ggplot(aes(x = mtry, y = SORENSEN_mean)) +
      geom_line() +
      theme_classic()
    ggsave(filename = here(perf_dir, paste0(sp[i],' hyp_raf_sorensen',  '.png')), dpi=200)
    h %>% ggplot(aes(x = mtry, y = AUC_mean)) +
      geom_line() +
      theme_classic()
    ggsave(filename = here(perf_dir, paste0(sp[i],' hyp_raf_auc',  '.png')), dpi=200)
    h %>% ggplot(aes(x = mtry, y = TSS_mean)) +
      geom_line() +
      theme_classic()
    ggsave(filename = here(perf_dir, paste0(sp[i],' hyp_raf_tss',  '.png')), dpi=200)
    readr::write_tsv(x = h, file = here(perf_dir, paste0(sp[i], " hyp_raf.txt")))
  }
  
  
  #### Support Vector Machine ####
  m_svm <- esm_svm(
    data = pa,
    response = "pr_ab",
    predictors = env_used,
    partition = ".part",
    thr = "max_sorensen"
  )
  
  #### Generalized Additive Model ####
  n_t <- flexsdm:::n_training(data = pa, partition = ".part")
  
  candidate_k <- 20
  while(any(n_t < flexsdm:::n_coefficients(data = replace_absences(pa, b), 
                                           predictors = env_used, 
                                           k = candidate_k))){
    candidate_k <- candidate_k-3
  }
  
  m_gam <- esm_gam(
    data = replace_absences(pa, b),
    response = "pr_ab",
    predictors = env_used,
    partition = ".part",
    thr = "max_sorensen",
    k = candidate_k
  )
  
  #### Generalized Linear Models ####
  if(sum(pa$pr_ab==1) >= 2){
    m_glm <- esm_glm(
      data = replace_absences(pa, b),
      response = "pr_ab",
      predictors = env_used,
      partition = ".part",
      thr = "max_sorensen",
      poly = 2
    )
  }
  
  
  #### Gaussian Process ####
  m_gau <- esm_gau(
    data = pa,
    response = "pr_ab",
    predictors = env_used,
    partition = ".part",
    # background = b,
    thr = "max_sorensen"
  )
  
  
  models <- grep("m_", ls(), value = TRUE)
  filt <- sapply(models, function(x) {
    length(get(x))
  })
  models <- models[filt > 0]
  
  # Filter by performance
  filt <- flexsdm::sdm_summarize(lapply(models, get))
  filt_perf <- 
    filt$SORENSEN_mean >= 0.7 &
    filt$thr_value != 0 &
    filt$thr_value != 1
  models_perf <- models[filt_perf]
  
  ##### Ensemble ####
  # if(length(models_perf)>1){
  #   m_ensemble <-
  #     flexsdm::fit_ensemble(
  #       lapply(models_perf, get),
  #       ens_method = c("mean"),
  #       thr_model = "max_sorensen",
  #       metric = "SORENSEN",
  #       thr = "max_sorensen"
  #     )  
  # } else if(length(models_perf)==1) {
  #   m_ensemble <-
  #     flexsdm::fit_ensemble(
  #       lapply(c(models_perf, models_perf), get),
  #       ens_method = c("mean"),
  #       thr_model = "max_sorensen",
  #       metric = "SORENSEN",
  #       thr = "max_sorensen"
  #     )
  # }
  
  models <- grep("m_", ls(), value = TRUE)
  filt <- sapply(models, function(x) {
    length(get(x))
  })
  models <- models[filt > 0]
  
  # Model performance
  performance <- flexsdm::sdm_summarize(lapply(models, function(x) {
    if (length(get(x)) > 0) {
      get(x)
    }
  }))
  
  readr::write_tsv(x = performance, file = here(perf_dir, paste0(sp[i], "_models_performance.txt")))
  
  if(length(models_perf)>=1){
    ##%######################################################%##
    #                                                          #
    ####             Predict individual models              ####
    #                                                          #
    ##%######################################################%##
    # models <- models[models!="m_ensemble"]
    models <- models[models!="m_ensemble"]
    
    message("Predicting models for species ", sp[i], " ", i)
    
    
    models_object <- lapply(models[filt_perf], function(x) {
      get(x)
    })
    models_object[[1]]$esm_model
    
    prd <-
      flexsdm::sdm_predict( # ver error in sdm_predict 
        models = models_object,
        pred = env_projection,
        thr = c("max_sorensen"),
        con_thr = TRUE,
        clamp = TRUE,
        pred_type = 'cloglog',
        predict_area = NULL
      )
    
    for(mm in 1:length(prd)){
      terra::writeRaster(prd[[mm]],
                         # here('1-SDM/2_Outputs/1_Current',
                         file.path('D:/SerpientesColombia/Tesis/Modelos/1-SDM/2_Outputs/1_Current',
                              'Algorithm',
                              names(prd[mm]),
                              paste0(sp[i], '.tif'))
                         , overwrite=TRUE)
    }
    rm(prd)
    # Calculat ensamble only using rasters based on mean approach
    
    # m_ensemble 
    prd <- 
    file.path('D:/SerpientesColombia/Tesis/Modelos/1-SDM/2_Outputs/1_Current',
              'Algorithm', names(prd[mm]), paste0(sp[i], '.tif'))
    b
    
    terra::writeRaster(prd[[1]],
                       here('1-SDM/2_Outputs/1_Current',
                            'Algorithm',
                            "mean",
                            paste0(sp[i], '.tif'))
                       , overwrite=TRUE)
    rm(list=grep("m_", ls(), value = TRUE))
  }
}






##%######################################################%##
#                                                          #
####            Save final performance table             ####
#                                                          #
##%######################################################%##

thr <- "./1-SDM/2_Outputs/0_Model_performance/" %>% list.files(pattern = "_models_performance_2.txt", full.names = TRUE)
thr2 <- lapply(thr, readr::read_tsv)
names(thr2) <- gsub("_models_performance_2.txt", "", basename(thr))
thr2 <- bind_rows(thr2, .id="species")
# readr::write_tsv(thr2,  "./1-SDM/2_Outputs/0_Model_performance/00_model_performance.txt")


