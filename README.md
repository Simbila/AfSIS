# AfSIS
Data Science
#' OCP-Tanzania cropping system productivity indicators
#' W. Simbila, August 2017

require(downloader)
require(rgdal)
require(raster)
require(quantreg)

# Data setup----------------------------------------------------------------
# Create a data folder in your current working directory
dir.create("OCP_npp", showWarnings=F)
setwd("./OCP_npp")

# Download
# MobileSurvey data
download("https://www.dropbox.com/s/m1ilf42a5xybhgo/TZ_geos_082317.csv?raw=1", mode="wb")
crps <- read.table("TZ_geos_082317.csv?raw=1", header=T, sep=",")

#Productivity grids
download("https://www.dropbox.com/s/m1ilf42a5xybhgo/TZ_geos_082317.csv?raw=1", mode="wb")
glist <- list.files(pattern="tif", full.names=T)
grids <- stack(glist)

# Overlay with grided covariates ---------------------------------------------
# Project survey coords to grid CRS
crps.proj <- as.data.frame(project(cbindCcrps$lon,crps$lat), "+proj=laea +ellps=WGS84 +long_0=35 +lat_0=6 +units=m +no_defs"))
colnames(crps.proj) <- c("x","y") ## laea coordinates
crps <- cbind(crps,crps.proj)
coordinates(crps) <- ~x+y
projection(crps) <- projection(grid)

#Extract gridded variables at MobileSurvey locations
crpsgrid <- extract(grids,crps)
crps <- as.data.frame(crps)
crps <- cbind.data.frame(crps,crpsgrid)
crps <- unique(na.omit(crps)) ## includes only unique & complete records

# Quantile regressions ---------------------------------------------------------
# Long-term Average Net Primary Productivity (NPP, kg/h yr)
# 2000-2014 MODIS MOD17A3 data (ftp://africagrids.net/500m/MOD17A3H)
NPPa <- rq(I(NPPa*10000)~MZP+SGP+LGP+RCP+RCP+OCP+LVS,tau=c(0.10,0.5,0.9), data=crps)
print(NPPa)

# Long-term interannual NPP standard deviation
# 2000-2014 MODIS MOD17A3 data (ftp://africagrids.net/500m/MOD17A3H)
NPPs <- rq(I(NPPs*10000)~MZP+SGP+LGP+RCP+OCP+LVS, tau=c(0.10,0.5,0.9), =crps)
print(NPPs)

# Mean Annual Precipitation (MAP,mm/yr)
# 2000-2014 CHIRPS data (ftp://africagrids.net/5000m/CHIRPS/Annual/sum)
crps$MAP <- ifelse(crps$MAP==0,NA,crps$MAP)
MAP <- rq(MAP~MZP+SGP+LGP+RCP+OCP+LVS,tau=c(0.10,0.5,0.9), data=crps)
print(MAP)

# Rain Use Efficiency (NPPa/MAP)
crps$RUE <- (crps$NPPa*1000)/crps$MAP
RUE <- rq(RUE~MZP+SGP+LGP+RCP+OCP+LVS, tau=c(0.10,0.5,0.9), data=crps)
print(RUE)

# NPP residual (following median regression against MAP)
NPM <- rq(I(NPPa*10000)~MAP, tau=0.5, data=crps)
crps$NPPr <- crps$NPPa*10000 - predict(NPM, crps)
NPPr <- rq(NPPr~MZP+SGP+LGP+RCP+OCP+LVS, tau=c(0.1,0.5,0.9), data=crps)
print(NPPr)
