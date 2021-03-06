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
list.of.packages <- c("tidyverse", "lubridate","chron","bcdata", "bcmaps","sf", "rgdal", "nngeo","Cairo","OpenStreetMap", "ggmap",
                      "leaflet", "lunar", "zoo", "colortools", "RColorBrewer", "viridis","osmdata", "ggspatial", "gridExtra","cowplot")


# Check you have them and load them
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)
rm(list.of.packages, new.packages) # for housekeeping
#####################################################################################

load("./data/01_load.RData")
glimpse(lttrap) # marten live trap trap data
glimpse(ltdat)
ltdat %>% group_by(`Trapping session`) %>% count(Status)

# `Trapping session` Status     n
# 1 96-97              MAAM     280
# 2 97-98              MAAM     273
# 3 98-99              MAAM      80
# 4 99-00              MAAM     160


###--- Load formatted data
create_trap_map <- function(lttrap=lttrap, ltdat=ltdat, year=year){
  # lttrap=lttrap[[2]]
  # ltdat=ltdat
  # year="97"
   
  load(paste0("./out/MartenData_19",year,".RData"))
  colnames(lttrap)[2:3] <- c("x","y")
  
  traps <- marten.data$traps
  traps$x <- traps$x*1000
  traps$y <- traps$y*1000
  traps <- left_join(traps, lttrap %>% select(Trap_ID, x, y), by=c("x", "y"))
  
  
  m.obs <- ltdat %>% filter(grepl(year, `Trapping session`)) %>% group_by(Trap_ID) %>% count(Status)
  colnames(m.obs)[3] <- c("Martens per Trap")
  
  traps <- left_join(traps, m.obs %>% select(-Status))
  traps[is.na(traps)] <- 0
  
  ###--- create sf object of trap data
  traps.sf <- st_as_sf(traps, coords=c("x","y"), crs=26910)
  
  ###--- view OSM data and download appropriate section for study area
  traps.latlon <- st_transform(traps.sf, crs=4326)
  bbox <- st_bbox(traps.latlon)
  
  # use latlon for entire study area
  
  LAT1 = bbox[2] ; LAT2 = bbox[4]
  LON1 = bbox[3] ; LON2 = bbox[1]
  
  #our background map
  map <- openmap(c(LAT2+0.05,LON1+0.05), c(LAT1-0.05,LON2-0.05), 
                 zoom = NULL,
                 type = c("osm", "stamen-toner", "stamen-terrain","stamen-watercolor", "esri","esri-topo")[6],
                 mergeTiles = TRUE)
  
  ## OSM CRS :: "+proj=merc +a=6378137 +b=6378137 +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +k=1 +units=m +nadgrids=@null +wktext +no_defs"
  map.latlon <- openproj(map, projection = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  
  traps.latlon$Longitude <- st_coordinates(traps.latlon)[,1]
  traps.latlon$Latitude <- st_coordinates(traps.latlon)[,2]
  
  trapyear <- as.numeric(paste0("19",year))
  plot.subtitle <- paste0(trapyear," - ",trapyear+1)
  
  traps.plot <- OpenStreetMap::autoplot.OpenStreetMap(map.latlon)  +
    labs(subtitle = plot.subtitle)+
    geom_sf(data=traps.latlon[traps.latlon$`Martens per Trap`!=0,],
            aes(x=Longitude, y=Latitude, size=`Martens per Trap`), col="darkblue")+
    geom_sf(data=traps.latlon[traps.latlon$`Martens per Trap`==0,],
            aes(x=Longitude, y=Latitude), col="cadetblue")
  
   return(list(traps.sf=traps.sf, traps.plot=traps.plot))
}

###--- create trap plots for the 4 live traps years
years <- c("96", "97", "98", "99")

traps.out = list()
for(i in 2:length(years)){
  traps.out[[i]] <- create_trap_map(lttrap=lttrap[[i]], ltdat=ltdat, year=years[i])
}

Cairo(file="out/lt99_plot.PNG",type="png",width=2200,height=2000,pointsize=12,bg="white",dpi=300)
traps.out[[4]]$traps.plot
dev.off()

grid.arrange(traps.out[[1]]$traps.plot,traps.out[[2]]$traps.plot,
             traps.out[[3]]$traps.plot,traps.out[[4]]$traps.plot, 
             nrow=2,
             top="Williston Basin Live Trap Locations")

#####################################################################################

###--- determine distance of traps and group those that are too close
# need to convert to utm for m distance
# espg 26910
# to visualise trap distances
# traps.utm <- st_transform(traps.out[[2]]$traps.sf, crs=26910)
# st_bbox(traps.utm)
# 
# nrow(traps.utm)
# traps.utm$Trap_ID
# traps.utm.3km <- st_buffer(traps.utm, dist=3000)
# traps.utm.5km <- st_buffer(traps.utm, dist=5000)
# traps.utm.10km <- st_buffer(traps.utm, dist=10000)
# 
# ggplot(traps.utm) +
#   # geom_sf(aes(col = TrapGroup)) +
#   geom_sf_label(aes(label = Trap_ID))
# 
# ggplot()+
#   geom_sf(data=traps.utm)+
#   geom_sf(data=traps.utm.3km, fill=NA, col="blue")+
#   geom_sf(data=traps.utm.5km, fill=NA, col="red")+
#   geom_sf(data=traps.utm.10km, fill=NA, col="green")
#   
# #- distances between traps
# traps.dist <- st_nn(traps.utm, traps.utm, k=nrow(traps.utm), maxdist=5000, returnDist = T, sparse=TRUE)
# 
# not sure going by distance is the best method as lots overlap in distance...
# changing tactics and setting down a 'grid' and grouping traps that fall within arbitrary grid
# some traps will be closer to others in different grid but not sure how else to group

################################################################################
###--- function to create grid
# note that this function uses hexagon grids to group points (traps)

find_grid <- function (input=input, cellsize=cellsize){
  
  aoi_utm <- st_transform(input, crs=26910) # to have in metres for specifying grid cell size
  aoi_grid <- st_make_grid(st_bbox(aoi_utm), cellsize=cellsize, square=FALSE) #  grid for entire AOI (rectangle)
  
  tg.dist <- st_nn(aoi_utm, aoi_grid, k=1, returnDist = T)
  aoi_utm$tg_value <- tg.dist$nn
  
  tmp <- as.data.frame(aoi_utm %>% count(tg_value) %>% st_drop_geometry())
  tmp$Trap_Grp <- rownames(tmp)
  aoi_utm$Trap_Grp <- tmp$Trap_Grp[match(aoi_utm$tg_value, tmp$tg_value)]
  aoi_utm <- aoi_utm %>% dplyr::select(-c(tg_value))
  
  return(list(aoi_utm=aoi_utm, aoi_grid=aoi_grid))
}

# run the find_grid function for all trap years
traps.grid <- list()
for(i in 1:4){
  traps.grid[[i]] <- find_grid(input=traps.out[[i]]$traps.sf, cellsize=5000)
}

t96 <- traps.grid[[1]]$aoi_utm
t97 <- traps.grid[[2]]$aoi_utm
t98 <- traps.grid[[3]]$aoi_utm
t99 <- traps.grid[[4]]$aoi_utm

ggplot()+
  geom_sf(data=t99, aes(col=Trap_Grp))

ggplot(t99)+
  geom_sf(aes(col=Trap_Grp))+
  geom_sf_label(aes(label = Trap_Grp))

hist(as.numeric(t96$Trap_Grp), breaks = length(unique(t96$Trap_Grp)))
hist(as.numeric(t97$Trap_Grp), breaks = length(unique(t96$Trap_Grp)))
hist(as.numeric(t98$Trap_Grp), breaks = length(unique(t96$Trap_Grp)))
hist(as.numeric(t99$Trap_Grp), breaks = length(unique(t96$Trap_Grp)))

save(traps.grid, file = paste0("./out/retro_traps_grp.RData"))


##################################################################################
###--- Create map for recent trap data
load("out/rec.data.out.RData")
# rec.data.out <- list(rec_effort=rec_effort, rec_ydata=rec_ydata)
rec.data.out$rec_grid_output$aoi_utm
traps.centroid <- as.data.frame(rec.data.out$grid_centroid_utm)
traps.centroid$Martens <- rowSums(rec.data.out$rec_ydata)
traps.centroid$MartenOcc <- ifelse(traps.centroid$Martens>1,1,0)

traps.sf <- st_as_sf(traps.centroid, coords=c("X","Y"), crs=26910)

grid <- rec.data.out$rec_grid_output$fishnet_grid_sf
grid <- st_as_sf(grid)
grid.latlon <- st_transform(grid, crs=4326)

grid$focus <- rec.data.out$rec_grid_output$aoi_utm$Grid_Focus[match(grid$grid_id,rec.data.out$rec_grid_output$aoi_utm$Trap_Grp)]

hs2019_gridtraps_plot <- ggplot()+
  geom_sf(data=grid, fill=NA)+
  geom_sf(data=rec.data.out$rec_grid_output$aoi_utm, aes(colour=Grid_Focus))+
  scale_colour_manual(values = c("black","red"))+
  theme(legend.position = "none")


Cairo(file="out/hs2019_gridtraps_plot.PNG",type="png",width=2200,height=2000,pointsize=12,bg="white",dpi=300)
hs2019_gridtraps_plot
dev.off()

# output for Table 2
Marten19 <- as.data.frame(rowSums(rec.data.out$rec_ydata))
colnames(Marten19) <-"Marten19"
Marten19$grid_id <- rownames(Marten19)
Marten19$n <- rec.data.out$rec_effort[,1]
Marten19$Effort19 <- rowSums(rec.data.out$rec_effort)
Marten19$DtnRate <- Marten19$Marten19/Marten19$Effort19
# sum(Marten96$Marten96); sum(Marten97$Marten97); sum(Marten19$Marten19)

sum(Marten19$n); nrow(Marten19); mean(Marten19$n); standard_error(Marten19$n) # 184 traps; 122 grid cells; 1.5 mean, 0.1 SE
nrow(Marten19[Marten19$Marten19>0,])/nrow(Marten19) # 0.31 naive occupancy
mean(Marten19$Effort19); standard_error(Marten19$Effort19) # 6.0, 0.4
mean(Marten19$DtnRate)*100; standard_error(Marten19$DtnRate)*100 # 5.9, 0.9


###--- view OSM data and download appropriate section for study area
traps.latlon <- st_transform(traps.sf, crs=4326)
bbox <- st_bbox(grid.latlon)

# use latlon for entire study area

LAT1 = bbox[2] ; LAT2 = bbox[4]
LON1 = bbox[3] ; LON2 = bbox[1]

#our background map
map <- openmap(c(LAT2+0.05,LON1+0.05), c(LAT1-0.05,LON2-0.05), 
               zoom = NULL,
               type = c("osm", "stamen-toner", "stamen-terrain","stamen-watercolor", "esri","esri-topo")[6],
               mergeTiles = TRUE)

## OSM CRS :: "+proj=merc +a=6378137 +b=6378137 +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +k=1 +units=m +nadgrids=@null +wktext +no_defs"
map.latlon <- openproj(map, projection = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

traps.latlon$Longitude <- st_coordinates(traps.latlon)[,1]
traps.latlon$Latitude <- st_coordinates(traps.latlon)[,2]

traps.plot.2019 <- OpenStreetMap::autoplot.OpenStreetMap(map.latlon)  +
  labs(subtitle = "2019-2020")+
  geom_sf(data=traps.latlon[traps.latlon$Martens!=0,],
          aes(x=Longitude, y=Latitude, size=Martens), col="darkblue")+
  geom_sf(data=traps.latlon[traps.latlon$Martens==0,],
          aes(x=Longitude, y=Latitude), col="cadetblue")+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  theme(legend.position = "bottom")

Cairo(file="out/hs2019_cnt_plot.PNG",type="png",width=2200,height=2000,pointsize=12,bg="white",dpi=300)
traps.plot.2019
dev.off()

traps.plot.Occ.2019 <- OpenStreetMap::autoplot.OpenStreetMap(map.latlon)  +
  geom_sf(data=traps.latlon[traps.latlon$MartenOcc==1,],
          aes(x=Longitude, y=Latitude), col="darkblue")+
  geom_sf(data=traps.latlon[traps.latlon$MartenOcc==0,],
          aes(x=Longitude, y=Latitude), col="cadetblue")+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank()) 

Cairo(file="out/hs2019_occ_plot.PNG",type="png",width=2200,height=2000,pointsize=12,bg="white",dpi=300)
traps.plot.Occ.2019
dev.off()

###################################################
###--- Create map for 1996 and 1997 trap data
# output for Table 2 in report

# standard error function
standard_error <- function(x, na.rm=T) sd(x, na.rm=T) / sqrt(length(x)) 

load("out/retro.data.out9697.RData")
load("out/aoi_grid.RData")
names(retro.data.out[[2]])
aoi_grid

grid_centroid_utm <- st_coordinates(st_centroid(aoi_grid))
traps.centroid <- as.data.frame(grid_centroid_utm)
traps.centroid <- cbind(traps.centroid, aoi_grid[,4:7] %>% st_drop_geometry())

traps.sf <- st_as_sf(traps.centroid, coords=c("X","Y"), crs=26910)
Marten96 <- retro.data.out[[1]]$traps.sf %>% count(grid_id) %>% st_drop_geometry()
Marten96$Marten96 <- rowSums(retro.data.out[[1]]$y_21day)
Marten96$Effort96 <- rowSums(retro.data.out[[1]]$effort.21days)
Marten96$DtnRate <- Marten96$Marten96/Marten96$Effort96

Marten97 <- retro.data.out[[2]]$traps.sf %>% count(grid_id) %>% st_drop_geometry()
Marten97$Marten97 <- rowSums(retro.data.out[[2]]$y_21day)
Marten97$Effort97 <- rowSums(retro.data.out[[2]]$effort.21days)
Marten97$DtnRate <- Marten97$Marten97/Marten97$Effort97

# output for Table 2
sum(Marten96$n); nrow(Marten96); mean(Marten96$n); standard_error(Marten96$n) # 77 traps; 47 grid cells; 1.6 mean, 0.12 SE
nrow(Marten96[Marten96$Marten96>0,])/nrow(Marten96) # 0.85 naive occupancy
mean(Marten96$Effort96); standard_error(Marten96$Effort96) # 70.0, 7.9
mean(Marten96$DtnRate)*100; standard_error(Marten96$DtnRate)*100 # 14.9, 3.1

sum(Marten97$n); nrow(Marten97); mean(Marten97$n); standard_error(Marten97$n) # 77 traps; 47 grid cells; 1.6 mean, 0.12 SE
nrow(Marten97[Marten97$Marten97>0,])/nrow(Marten97) # 0.85 naive occupancy
mean(Marten97$Effort97); standard_error(Marten97$Effort97) # 51.2, 5.6
mean(Marten97$DtnRate, na.rm=T)*100; standard_error(Marten97$DtnRate)*100 # 7.4, 1.0

traps.sf$Marten96 <- Marten96$n[match(rownames(traps.sf),Marten96$grid_id)]
traps.sf$Marten97 <- Marten96$n[match(rownames(traps.sf),Marten97$grid_id)]

grid.latlon <- st_transform(aoi_grid, crs=4326)

###--- view OSM data and download appropriate section for study area
traps.latlon <- st_transform(traps.sf, crs=4326)
bbox <- st_bbox(grid.latlon)

# use latlon for entire study area

LAT1 = bbox[2] ; LAT2 = bbox[4]
LON1 = bbox[3] ; LON2 = bbox[1]

#our background map
map <- openmap(c(LAT2+0.05,LON1+0.05), c(LAT1-0.05,LON2-0.05), 
               zoom = NULL,
               type = c("osm", "stamen-toner", "stamen-terrain","stamen-watercolor", "esri","esri-topo")[6],
               mergeTiles = TRUE)

## OSM CRS :: "+proj=merc +a=6378137 +b=6378137 +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +k=1 +units=m +nadgrids=@null +wktext +no_defs"
map.latlon <- openproj(map, projection = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

traps.latlon$Longitude <- st_coordinates(traps.latlon)[,1]
traps.latlon$Latitude <- st_coordinates(traps.latlon)[,2]
traps.latlon[is.na(traps.latlon)] <-0


traps.plot <- OpenStreetMap::autoplot.OpenStreetMap(map.latlon)  +
  labs(subtitle = "1997-1998")+
  geom_sf(data=traps.latlon[traps.latlon$Marten97!=0,],
          aes(x=Longitude, y=Latitude, size=Marten97), col="darkblue")+
  geom_sf(data=traps.latlon[traps.latlon$Marten96!=0,],
          aes(x=Longitude, y=Latitude, size=Marten96), col="cadetblue")+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  theme(legend.position = "bottom")


Cairo(file="out/lt1997_cnt_plot.PNG",type="png",width=2200,height=2000,pointsize=12,bg="white",dpi=300)
traps.plot.1997
dev.off()

traps.plot.1996 <- OpenStreetMap::autoplot.OpenStreetMap(map.latlon)  +
  labs(subtitle = "1996-1997")+
  geom_sf(data=traps.latlon[traps.latlon$Martens!=0,],
          aes(x=Longitude, y=Latitude, size=Martens), col="darkblue")+
  geom_sf(data=traps.latlon[traps.latlon$Martens==0,],
          aes(x=Longitude, y=Latitude), col="cadetblue")+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  theme(legend.position = "bottom")

Cairo(file="out/lt1996_cnt_plot.PNG",type="png",width=2200,height=2000,pointsize=12,bg="white",dpi=300)
traps.plot.1996
dev.off()

traps.plot.Occ.1997 <- OpenStreetMap::autoplot.OpenStreetMap(map.latlon)  +
  geom_sf(data=traps.latlon[traps.latlon$MartenOcc==1,],
          aes(x=Longitude, y=Latitude), col="darkblue")+
  geom_sf(data=traps.latlon[traps.latlon$MartenOcc==0,],
          aes(x=Longitude, y=Latitude), col="cadetblue")+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank()) 

Cairo(file="out/lt1997_occ_plot.PNG",type="png",width=2200,height=2000,pointsize=12,bg="white",dpi=300)
traps.plot.Occ.1997
dev.off()

traps.plot.Occ.1996 <- OpenStreetMap::autoplot.OpenStreetMap(map.latlon)  +
  geom_sf(data=traps.latlon[traps.latlon$MartenOcc==1,],
          aes(x=Longitude, y=Latitude), col="darkblue")+
  geom_sf(data=traps.latlon[traps.latlon$MartenOcc==0,],
          aes(x=Longitude, y=Latitude), col="cadetblue")+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank()) 

Cairo(file="out/lt1996_occ_plot.PNG",type="png",width=2200,height=2000,pointsize=12,bg="white",dpi=300)
traps.plot.Occ.1996
dev.off()

###################################################

Cairo(file="out/WB_marten969719.PNG",type="png",width=3800,height=2200,pointsize=14,bg="white",dpi=300)
plot_grid(traps.plot.1996,traps.plot.1997,traps.plot.2019,nrow=1,
          rel_heights = c(1,1,1.5), rel_widths = c(0.95,1,1.15))
dev.off()

Cairo(file="out/WB_traps969719.PNG",type="png",width=3800,height=2200,pointsize=14,bg="white",dpi=300)
plot_grid(lt1996_gridtraps_plot,lt1997_gridtraps_plot,hs2019_gridtraps_plot,nrow=1,
          rel_heights = c(1,1,1.5), rel_widths = c(0.95,1,1.15))
dev.off()
