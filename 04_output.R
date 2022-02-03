# Copyright 2021 Province of British Columbia
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
# http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and limitations under the License.

# https://cran.r-project.org/web/packages/bayesplot/vignettes/plotting-mcmc-draws.html

# For Friday 5 Nov 2021
# In OneNote or a word document put together data wrangling output and analysis results
# Create / save plots showing capture and effort data (to make the case for just using 1996/97 and 1997/98 for retrospective data)
# Create graphs with simulated output (nmix and SCR)
# Create violin plots with actual data output (nmix and SCR)
# Clean code and write (brief) methods so know what I did when going back later

#####################################################################################
# 04_output.R
# script to produce output (figures / tables) from trap data and results
# written by Joanna Burgar (Joanna.Burgar@gov.bc.ca) - 02-Feb-2022
#####################################################################################

version$major
version$minor
R_version <- paste0("R-",version$major,".",version$minor)

.libPaths(paste0("C:/Program Files/R/",R_version,"/library")) # to ensure reading/writing libraries from C drive
tz = Sys.timezone() # specify timezone in BC
 

# Load Packages
list.of.packages <- c("tidyverse", "lubridate","chron","bcdata", "bcmaps","sf", "rgdal", "Cairo","OpenStreetMap", "ggmap",
                      "leaflet", "lunar", "zoo", "colortools", "RColorBrewer", "viridis","osmdata", "ggspatial") # "gridExtra", "grid"


# Check you have them and load them
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)
rm(list.of.packages, new.packages) # for housekeeping
#####################################################################################

# ###--- function to retrieve geodata from BCGW
# 
# retrieve_geodata_aoi <- function (ID=ID){
#   aoi.geodata <- bcdc_query_geodata(ID) %>%
#     filter(BBOX(st_bbox(aoi))) %>%
#     collect()
#   aoi.geodata <- aoi.geodata %>% st_intersection(aoi)
#   aoi.geodata$Area_km2 <- st_area(aoi.geodata)*1e-6
#   aoi.geodata <- drop_units(aoi.geodata)
#   return(aoi.geodata)
# }

#################################################################################
load("./data/01_load.RData")
glimpse(lttrap) # marten live trap trap data
colnames(lttrap[[1]])[2:3] <- c("x","y")
glimpse(ltdat)
ltdat %>% group_by(`Trapping session`) %>% count(Status)

# `Trapping session` Status     n
# 1 96-97              MAAM     280
# 2 97-98              MAAM     273
# 3 98-99              MAAM      80
# 4 99-00              MAAM     160

m.obs.96 <- ltdat %>% filter(`Trapping session`=="96-97") %>% group_by(Trap_ID) %>% count(Status)
colnames(m.obs.96)[3] <- c("Martens Trapped")
m.obs.97 <- ltdat %>% filter(`Trapping session`=="97-98") %>% group_by(Trap_ID) %>% count(Status)
m.obs.98 <- ltdat %>% filter(`Trapping session`=="98-99") %>% group_by(Trap_ID) %>% count(Status)
m.obs.99 <- ltdat %>% filter(`Trapping session`=="99-00") %>% group_by(Trap_ID) %>% count(Status)

###--- Load formatted data
load("./out/MartenData_1996.Rda")
glimpse(marten.data)
lt96.traps <- marten.data$traps
lt96.traps$x <- lt96.traps$x*1000
lt96.traps$y <- lt96.traps$y*1000
lt96.traps <- left_join(lt96.traps, lttrap[[1]] %>% select(Trap_ID, x, y), by=c("x", "y"))
lt96.traps <- left_join(lt96.traps, m.obs.96 %>% select(-Status))
lt96.traps[is.na(lt96.traps)] <- 0

###--- create sf object of trap data
lt96.traps.sf <- st_as_sf(lt96.traps, coords=c("x","y"), crs=26910)

###--- view OSM data and download appropriate section for study area
lt96.traps.latlon <- st_transform(lt96.traps.sf, crs=4326)
lt96_bbox <- st_bbox(lt96.traps.latlon)

# use latlon for entire study area

LAT1 = lt96_bbox[2] ; LAT2 = lt96_bbox[4]
LON1 = lt96_bbox[3] ; LON2 = lt96_bbox[1]

#our background map
# library("OpenStreetMap")
# library("rJava")
map <- openmap(c(LAT2+0.05,LON1+0.05), c(LAT1-0.05,LON2-0.05), 
               zoom = NULL,
               type = c("osm", "stamen-toner", "stamen-terrain","stamen-watercolor", "esri","esri-topo")[6],
               mergeTiles = TRUE)

## OSM CRS :: "+proj=merc +a=6378137 +b=6378137 +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +k=1 +units=m +nadgrids=@null +wktext +no_defs"
map.latlon <- openproj(map, projection = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

lt96.traps.latlon$Longitude <- st_coordinates(lt96.traps.latlon)[,1]
lt96.traps.latlon$Latitude <- st_coordinates(lt96.traps.latlon)[,2]


lt96_plot <- OpenStreetMap::autoplot.OpenStreetMap(map.latlon)  +
  labs(title = "Live Trap Locations\nWilliston Basin", subtitle = "1996-1997", x = "Longitude", y="Latitude")+
  geom_sf(data=lt96.traps.latlon[lt96.traps.latlon$`Martens Trapped`!=0,],
          aes(x=Longitude, y=Latitude, size=`Martens Trapped`), col="darkblue")+
  geom_sf(data=lt96.traps.latlon[lt96.traps.latlon$`Martens Trapped`==0,],
          aes(x=Longitude, y=Latitude), col="cadetblue")


Cairo(file="out/lt96_plot.PNG",type="png",width=2200,height=2000,pointsize=12,bg="white",dpi=300)
lt96_plot
dev.off()
