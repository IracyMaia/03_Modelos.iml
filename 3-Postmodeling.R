## %######################################################%##
#                                                          #
####               Models post-processing                ####
#                                                          #
## %######################################################%##
{
  require(terra)
  require(dplyr)
  require(flexsdm)
  require(here)
  require(ggplot2)
}

## %######################################################%##
#                                                          #
####               Correct overprediction               ####
#                                                          #
## %######################################################%##
bmcp <- function(records, x, y, buffer, cont_suit) {
  data_pl <- data.frame(records[, c(x, y)])
  data_pl <- data_pl[grDevices::chull(data_pl), ]
  data_pl <- data.frame(
    object = 1, part = 1, data_pl,
    hole = 0
  )
  data_pl <- terra::vect(as.matrix(data_pl), type = "polygons")
  terra::crs(data_pl) <- terra::crs(cont_suit)
  data_pl <- terra::buffer(data_pl, width = buffer)
  hull <- terra::rasterize(data_pl, cont_suit)
  hull[is.na(hull)] <- 0
  result <- cont_suit * hull
  return(result)
}

# Colombia boundaries
col <- terra::vect("Colombia.gpkg")

occ <- readr::read_tsv("occurrences_cleaned_final_final.txt") # registros no filtrados!!!
occ <- occ %>% dplyr::select(species = accepted_name_bin, x, y)

# occ$species %>% unique %>% sort %>% grep("tif", ., value = T)
d <- "1-SDM/2_Outputs/1_Current/Algorithm/mean/" %>% list.files(pattern = ".tif", full.names = TRUE)
names(d) <- gsub(".tif$", "", basename(d))
occ <- occ %>% dplyr::filter(species %in% names(d))

sp <- names(d)
i <- 1
for (i in 1:length(sp)) {
  p <- occ %>%
    dplyr::filter(species == sp[i])
  r <- terra::rast(d[sp[i]])
  r <- bmcp(p, x = "x", y = "y", buffer = 100000, cont_suit = r) %>%
    terra::crop(., col) %>%
    terra::mask(., col)
  plot(r[[1]])
  points(p[-1])
  terra::writeRaster(r,
    here(
      "1-SDM/2_Outputs/1_Current/FinalModels",
      paste0(sp[i], ".tif")
    ),
    overwrite = TRUE
  )
}





## %######################################################%##
#                                                          #
####               Richness maps based on               ####
####       continuous suitability above threshold       ####
#                                                          #
## %######################################################%##

d <- list.files("./1-SDM/2_Outputs/1_Current/FinalModels/", pattern = ".tif$", full.names = TRUE)
names(d) <- basename(d) %>% gsub(".tif$", "", .)

check_colom_distr <- data.frame(species = names(d), zero_suitability = NA)
r <- terra::rast(d[1])[[2]]

plot(r)
for (i in 2:length(d)) {
  message("sp ", i)
  r0 <- terra::rast(d[i])[[2]]
  check_colom_distr[i, 2] <- terra::global(r0, max, na.rm = T)
  r <- r + r0
}
plot(r, col = pals::viridis(20))
hist(r)
# terra::writeRaster(r, "1-SDM/2_Outputs/1_Current/richnes_map.tif")

# Especies sin valores de "suitability" en Colombia
check_colom_distr %>% filter(zero_suitability == 0)
check_colom_distr %>% dplyr::arrange(zero_suitability)

# readr::write_tsv(check_colom_distr, "1-SDM/2_Outputs/1_Current/species_with_problem.txt")


## %######################################################%##
#                                                          #
####     relacion entre areas protegidas y riqueza      ####
#                                                          #
## %######################################################%##
# Poligono de colombia
colombia <- terra::vect("./Colombia.gpkg")

# areas protegidas
pas <- terra::vect("./areas_protegidas/WDPA_WDOECM_Mar2023_Public_COL_shp-polygons.shp")
pas <- terra::crop(pas, colombia)
plot(pas)

# mapa de riqueza
rich_col <- terra::rast("1-SDM/2_Outputs/1_Current/richnes_map.tif")
plot(rich_col)
plot(rich_col > 50)
plot(rich_col > 40)

# Rasterizar las AP
pas_r <- rasterize(pas, rich_col, cover = TRUE)
plot(pas_r)
pas_r <- pas_r > 0.5
pas_r[!is.na(rich_col) & is.na(pas_r)] <- 0
plot(pas_r)

# Relacion entre AP y riqueza
rich_pa <- c(rich_col, pas_r) %>%
  as.data.frame() %>%
  as_tibble() %>%
  na.omit()
names(rich_pa) <- c("Riqueza", "Protegido")

require(ggplot2)

p1 <- ggplot(rich_pa, aes(x = Riqueza)) +
  geom_histogram(aes(fill = Protegido),
    alpha = 0.5,
    position = "identity", bins = 30
  ) +
  labs(y = "N° de celdas") +
  theme_light() +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = pals::viridis(10)[c(1, 10)])
p1
class(p1)

# ggsave(plot = p1, "./ResultadosCongreso/Figura 2.png", dpi = 300,
#        units = "cm", width = 18, height = 15)


## %######################################################%##
#                                                          #
#### Representacion dentro de unidades de conservacion  ####
#                                                          #
## %######################################################%##
# Modelos de las especies
d <- list.files("./1-SDM/2_Outputs/1_Current/FinalModels/", pattern = ".tif$", full.names = TRUE)
names(d) <- basename(d) %>% gsub(".tif$", "", .)

df <- data.frame(species = names(d), total_area = NA, area_within_pas = NA)

i <- 45
for (i in 1:length(d)) {
  message("species n: ", i)
  r <- terra::rast(d[i])
  r <- (r[[2]] > 0)
  df[i, "total_area"] <- terra::global(r, sum, na.rm = TRUE)
  df[i, "area_within_pas"] <- terra::mask(r, pas_r, maskvalue = 0) %>% global(., sum, na.rm = TRUE)
  # plot(terra::mask(r, pas_r, maskvalue=0), col=pals::parula(2))
  # plot(r, col=pals::parula(2))
}

# Filtrar las especies que no tienen predicciones dentro de Colombia
df2 <- df %>% dplyr::filter(total_area > 0)

# Calculo de proporciones
df2 <- df2 %>% dplyr::mutate(proportion = area_within_pas / total_area)

# readr::write_tsv(df2, "./ResultadosCongreso/representacion dentro de AP.txt")

# Grafico
summary(df2$proportion)

# boxplot
df2 <- readr::read_tsv("./ResultadosCongreso/representacion dentro de AP.txt")

p3 <- ggplot(df2, aes(x = "", y = proportion)) +
  geom_violin() +
  geom_boxplot(fill = "blue", alpha = 0.5) +
  theme_light() +
  labs(y = "Representatividad dentro de AP")
p3
# trasformar el area total a km2
df2$total_area <- df2$total_area * 5 * 5
range(df2$total_area)

# 900000 # km2
# seq(0, 1150000, by=50000)
# cut(df2$total_area, breaks=c(0,1000), include.lowest=TRUE)
p4 <- ggplot(df2, aes(x = total_area, y = proportion)) +
  geom_point() +
  theme_light() +
  labs(
    x = "Area de distribucion (km2)",
    y = "Representatividad dentro de AP"
  )

require(patchwork)
p34 <- p3 + p4 + 
  plot_layout(widths = c(1, 3))

ggsave(plot = p34, "./ResultadosCongreso/Figura 3.png", 
       dpi = 300, units = "cm", 
       width = 18, height = 15)

