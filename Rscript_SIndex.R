# -----------------------------------------------------------------
# -----------------------------------------------------------------
# Exercise: Agroclimatic zoning 
#
# Goal: generate a suitability map for Cassava in Angola
# 
# Climate database: monthly data of avg temperature and rainfall (https://www.worldclim.org/data/bioclim.html#google_vignette)
# Soil database: https://files.isric.org/soilgrids/former/2017-03-10/aggregated/10km/
# Crop requirements: according to FAO-EcoCrop database
#
# Leonardo A. Monteiro, PhD
# 15th of May, 2021
# -----------------------------------------------------------------
# -----------------------------------------------------------------

# ====================================================
# --- setting working directory
# ====================================================
setwd("C:/Users/lomo0003/FAO_task/")

# ====================================================
# --- packages
# ====================================================
pkgs <- c("raster", "data.table", "ggplot2", "maptools", "rgeos")
lapply(pkgs, require, character.only = T)

# ====================================================
# --- Angola elevation and borders
# ====================================================
elevation <- getData("alt", country = "Angola", level = 0)
country <- getData("GADM", country = "Angola", level = 0)

writeRaster(elevation, "./outputs/elev_Angola.tif")
writeSpatialShape(country, "./shp/shp_Angola_borders1.shp")

shp_angola <- readShapeSpatial("./shp/shp_Angola_borders1.shp")

# ====================================================
# --- Input files
# ====================================================

#--- climate (monthly averages 1970-2000)
tair <- stack("./climate/processed/Tair.tif")
prec <- sum(stack("./climate/processed/Precipitation.tif"))

# --- defining thresholds (absmin, optmin, optmax and absmax) 

# air temperature
Tabsmin <- 10
Toptmin <- 20
Toptmax <- 29
Tabsmax <- 35
  
# precipitation
Pabsmin <- 500
Poptmin <- 1000
Poptmax <- 1500
Pabsmax <- 5000

# pH
pHabsmin <- 4
pHoptmin <- 5.5
pHoptmax <- 8
pHabsmax <- 9

# soil texture
TextSand1 <- 52
TextSand2 <- 65
TextClay1 <- 15
TextClay2 <- 27

# ====================================================
# --- Index for temperature (Itair) 
# ====================================================
TS <- raster()

for(j in 1:nlayers(tair)){
  t1 <- tair[[j]]
  t2 <- as.data.frame(t1, xy = T)
  colnames(t2)[3] <- "tair"
  
  t2$index <- ifelse(t2$tair < Tabsmin, 0, ifelse(
    t2$tair < Toptmin & t2$tair > Tabsmin,  100 * (t2$tair - Tabsmin) / (Toptmin - Tabsmin), ifelse(
      t2$tair < Toptmax & t2$tair > Toptmin,  100, ifelse(
        t2$tair < Tabsmax & t2$tair > Toptmax, 100 * (1 - (t2$tair - Toptmax) / (Tabsmax - Toptmax)), 0
        ))))
  t3 <- t2[, c("x", "y", "index")]
  r_t3 <- rasterFromXYZ(t3)
  names(r_t3) <- j
  r_t3[r_t3 == 0] <- NA
  TS <- stack(TS, r_t3)
}  

Itair <- TS[[6]] 

  # -- saving the output raster
writeRaster(Itair, "./outputs/Itair.tif", overwrite = T)
 
# ====================================================
# --- Index for precipitation (Iprec) 
# ====================================================
p1 <- prec
p2 <- as.data.frame(p1, xy = T)
colnames(p2)[3] <- "pp"

p2$index <- ifelse(p2$pp < Pabsmin, 0, ifelse(
  p2$pp < Poptmin & p2$pp > Pabsmin,  100 * (p2$pp - Pabsmin) / (Poptmin - Pabsmin), ifelse(
    p2$pp < Poptmax & p2$pp > Poptmin,  100, ifelse(
      p2$pp < Pabsmax & p2$pp > Poptmax, 100 * (1 - (p2$pp - Poptmax) / (Pabsmax - Poptmax)), 0
))))

p3 <- p2[, c("x", "y", "index")]
r_p3 <- rasterFromXYZ(p3)
names(r_p3) <- p
r_p3[r_p3 == 0] <- NA
Iprec <- r_p3  

  # -- saving the output raster
writeRaster(Iprec, "./outputs/Iprec.tif", overwrite = T)  

# ====================================================
# --- Index for pH (IpH)
# ====================================================
pH1 <- raster("./soil/pH_AGGREGATED_0_100cm.tif")
pH2 <- as.data.frame(pH1, xy =T)
colnames(pH2)[3] <- "pH"

pH2$index <- ifelse(pH2$pH < pHabsmin, 0, ifelse(
  pH2$pH < pHoptmin & pH2$pH > pHabsmin,  100 * (pH2$pH - pHabsmin) / (pHoptmin - pHabsmin), ifelse(
    pH2$pH < pHoptmax & pH2$pH > pHoptmin,  100, ifelse(
      pH2$pH < pHabsmax & pH2$pH > pHoptmax, 100 * (1 - (pH2$pH - pHoptmax) / (pHabsmax - pHoptmax)), 0
    ))))

pH3 <- pH2[, c("x", "y", "index")]
r_pH3 <- rasterFromXYZ(pH3)
names(r_pH3) <- 'pH'
r_pH3[r_pH3 == 0] <- NA
IpH <- r_pH3

  # -- resampling IpH to the same res of climate data
IpH_res <- resample(IpH, prec)
IpH_res[IpH_res > 100] <- 100
  
  # -- saving the output raster
writeRaster(IpH_res, "./outputs/IpH.tif", overwrite = T)

# ====================================================
# --- Index for soil texture (Isoiltext)
# ====================================================
sand <- raster("./soil/Sand_AGGREGATED_0_100cm.tif") # by mean of all depths
clay <- raster("./soil/Clay_AGGREGATED_0_100cm.tif") # by mean of all depths

sc <- brick(sand, clay)

df_sc <- na.omit(as.data.frame(sc, xy = T))
colnames(df_sc) <- c("x", "y", "sand", "clay")

df_sc$index <- ifelse(df_sc$sand < TextSand2, 100, ifelse(
  df_sc$sand < TextSand1 | df_sc$clay < TextClay2, 100, ifelse(
    df_sc$clay < TextClay1, 100, 25)))
df_sc1 <- df_sc[, c("x", "y", "index")]

r_texture <- rasterFromXYZ(df_sc1)

  # -- resampling to the same res of climate data
Isoiltext <- resample(r_texture, prec)
  
  # -- saving the output raster
writeRaster(Isoiltext, "./outputs/Isoiltext.tif", overwrite = T)


# =============================================================================================
# =============================================================================================
# CALCULATING THE INDEX OF LAND SUITABILITY FOR CASSAVA
# =============================================================================================
# =============================================================================================

# --- input files
tair <- raster("./outputs/Itair.tif")
prec <- raster("./outputs/Iprec.tif")
pH   <- raster("./outputs/IpH.tif")
text <- raster("./outputs/Isoiltext.tif")
elev <- raster("./outputs/elev_Angola.tif")

final_index <- 0.5 * (tair + (0.6 * pH + 0.2 * text + 0.2 * 100))

# --- resampling elevation 
elev_res <- resample(elev, final_index)

# --- filtering out pixels with elevation greater than 1500 m
final_index[elev_res > 1500] <- 0

# --- saving the output file 
writeRaster(final_index, "./outputs/SI_Cassava_Angola.tif", overwrite = T)

























