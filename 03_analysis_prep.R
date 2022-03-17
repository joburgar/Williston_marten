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

#####################################################################################
# 03_analysis_prep.R
# script to prep output for modelling
# written by Joanna Burgar (Joanna.Burgar@gov.bc.ca) - 13-Oct-2021
#####################################################################################
version$major
version$minor
R_version <- paste0("R-",version$major,".",version$minor)

.libPaths(paste0("C:/Program Files/R/",R_version,"/library")) # to ensure reading/writing libraries from C drive
tz = Sys.timezone() # specify timezone in BC

# Load Packages
list.of.packages <- c("bcdata","bcmaps","tidyverse", "lubridate","chron","sf","Cairo", 
                      "sf", "nngeo", "units","OpenStreetMap", "ggmap","rgdal","PNWColors")

# Check you have them and load them
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)
rm(list.of.packages, new.packages) # for housekeeping

################################################################################
#- function to group traps based on grid cell
find_grid <- function (input=input, cellsize=cellsize){
  # input=traps.sf
  # cellsize=5000
  aoi_utm <- st_transform(input, crs=26910) # to have in metres for specifying grid cell size
  aoi_utm$TrapNum <- rownames(aoi_utm)
  aoi_grid <- st_make_grid(st_bbox(aoi_utm), what="polygons", cellsize=cellsize, square=FALSE) #  grid for entire AOI (rectangle)
  
  # subset grid to just cells of interest
  aoi_grid <- aoi_grid[aoi_utm] 
  # To sf and add grid ID
  fishnet_grid_sf = st_sf(aoi_grid) %>%
    # add grid ID
    mutate(grid_id = 1:length(lengths(aoi_grid)))
  
  fishnet_grid_sf$Area_km2 <- st_area(fishnet_grid_sf)*1e-6
  fishnet_grid_sf <- drop_units(fishnet_grid_sf)
  fishnet_grid_sf %>% summarise(sum(Area_km2)) %>% st_drop_geometry() #714.5 km2 study area each cell is 21.65 km2
  
  tg.dist <- st_nn(aoi_utm, fishnet_grid_sf, k=1, returnDist = T)
  aoi_utm$Trap_Grp <- unlist(tg.dist$nn)
  
  Trap_Grp <- aoi_utm %>% st_drop_geometry()
  
  return(list(aoi_utm=aoi_utm, fishnet_grid_sf=fishnet_grid_sf, Trap_Grp=Trap_Grp))
}

#####################################################################################

###--- function to retrieve geodata from BCGW
droplevels.sfc = function(x, except, exclude, ...) x

retrieve_geodata_aoi <- function (ID=ID, aoi=aoi){
  aoi.geodata <- bcdc_query_geodata(ID) %>%
    filter(BBOX(st_bbox(aoi))) %>%
    collect()
  aoi.geodata <- aoi.geodata %>% st_intersection(aoi)
  aoi.geodata$Area_km2 <- st_area(aoi.geodata)*1e-6
  aoi.geodata <- drop_units(aoi.geodata)
  aoi.geodata <- droplevels.sfc(aoi.geodata, except = geometry)
  return(aoi.geodata)
}


#####################################################################################
###--- function to retrieve data from downloaded gdb or shapefile
# Read the feature class

retrieve_gdb_shp_aoi <- function (dsn=dsn, layer=layer){
  # if a gdb then the dsn should be the path all the way to the '.gdb'
  aoi.geodata <- st_read(dsn=dsn,layer=layer) %>%
    st_transform(crs=3005) %>% st_intersection(aoi)
  aoi.geodata$Area_km2 <- st_area(aoi.geodata)*1e-6
  aoi.geodata <- drop_units(aoi.geodata)
  aoi.geodata <- droplevels.sfc(aoi.geodata, except = geometry)
  return(aoi.geodata)
}


#####################################################################################
####################################################################################
###--- grab covariate data based on grid layout for each survey and also area of interest for 1997 and current

out.files <- list.files("./out/", pattern="*.RData")
out.files <- out.files[grepl("MartenData", out.files)]
out.files <- out.files[!grepl("2020", out.files)]

retro.traps.out <- list()

for(r in 1:length(out.files)){
  load(paste0("./out/",out.files[r]))
  traps <- marten.data$traps*1000
  traps.sf <- st_as_sf(traps, coords=c("x","y"), crs=26910)
  
  
  retro.traps.out[[r]] <- list(traps.sf=traps.sf)
}


load("out/recent_occ_data.RData")
recent_traps.sf <- st_as_sf(recent_occ_data[[2]], coords=c("Easting","Northing"), crs=26910)

all.traps <- rbind(retro.traps.out[[1]]$traps.sf %>% select(geometry),
                   retro.traps.out[[2]]$traps.sf %>% select(geometry),
                   recent_traps.sf %>% select(geometry))

# now create a 5 km grid buffer following the angle of the actual trap locations
aoi <- st_buffer(all.traps %>% st_transform(crs=26910), dist=5000)
aoi <- find_grid(input=aoi, cellsize = 3000)
aoi_grid <- aoi$fishnet_grid_sf

ggplot()+
  geom_sf(data=aoi_grid)+
  geom_sf(data=retro.traps.out[[1]]$traps.sf, col="black")+
  geom_sf(data=retro.traps.out[[2]]$traps.sf, col="red")+
  geom_sf(data=recent_traps.sf, col="blue")

retro_96 <- st_join(retro.traps.out[[1]]$traps.sf, aoi_grid) %>% 
  dplyr::select(grid_id) %>%
  st_drop_geometry()

retro_97 <- st_join(retro.traps.out[[2]]$traps.sf, aoi_grid) %>% 
  dplyr::select(grid_id) %>%
  st_drop_geometry()

recent_19 <- st_join(recent_traps.sf, aoi_grid) %>% 
  dplyr::select(grid_id) %>%
  st_drop_geometry()

# Takes too long - clipped in ArcCatalog instead
# ogrListLayers(paste0(getwd(),"/data/VRI2002_VEG_COMP_LYR_R1_POLY_FINAL.gdb"))
# retrieve_gdb_shp_aoi(dsn=paste0(getwd(),"/data/VRI2002_VEG_COMP_LYR_R1_POLY_FINAL.gdb"),
#                      layer="VEG_COMP_LYR_R1_POLY_FINALV4")

# save.image("03_analysis.prep.RData")
###--- create covariate dataframe to put covariate info
aoi_grid_id <- aoi_grid$grid_id
# save(aoi_grid, file = paste0("./out/aoi_grid.RData"))
# load("data/aoi_grid.RData") # for the complete aoi_grid at the end of this script

cov.df <- as.data.frame(array(seq_len(length(aoi_grid_id)),c(nrow(aoi_grid),1))) # for all cov data
colnames(cov.df) <- c("grid_id")

# load covariates from bcdata
# using the bc data warehouse option to clip to aoi
aoi <- aoi_grid %>% st_transform(3005) # for regional cov data
aoi_centroid <- st_centroid(aoi) %>% dplyr::select(grid_id) %>% st_transform(crs=26910) # to have in m

# all traps in SBS so this might not be a useful measure
# biogeoclimatic zones
# bcdc_search("Biogeoclimatic zone", res_format = "wms")
# 3: BEC Map (other, wms, kml)
# ID: f358a53b-ffde-4830-a325-a5a03ff672c3
# Name: bec-map
aoi.BEC <- retrieve_geodata_aoi(ID = "f358a53b-ffde-4830-a325-a5a03ff672c3", aoi=aoi)
aoi.BEC %>% group_by(MAP_LABEL) %>% summarise(Area_km2=sum(Area_km2)) %>% st_drop_geometry
aoi.BEC %>% group_by(ZONE) %>% summarise(Area_km2=sum(Area_km2)) %>% st_drop_geometry
aoi.BEC %>% summarise(Area_km2=sum(Area_km2)) %>% st_drop_geometry # 823 km2 for 1996 study area, at grid cell size # 1863 km2 for 2020
# proportion of SBS in cell? proportion of ESSF # 139 km2 of ESSF and 685 km2 of SBS overall (also <1 km2 of BAFA)

Cairo(file="out/ALL_BEC.PNG",type="png",width=2200,height=2000,pointsize=12,bg="white",dpi=300)
ggplot()+
  geom_sf(data=aoi.BEC, aes(fill=ZONE))+
  geom_sf(data=retro.traps.out[[1]]$traps.sf, col="black")+
  geom_sf(data=retro.traps.out[[2]]$traps.sf, col="red")+
  geom_sf(data=recent_traps.sf, col="blue")
dev.off()

cov.df$SBS_prop <- NA
for(i in seq_len(nrow(cov.df))){
  tmp <- aoi %>% st_transform(crs=26910) %>% filter(aoi$grid_id==i)
  tmp2 <- st_intersection(tmp, aoi.BEC %>% st_transform(crs=26910) %>% filter(ZONE=="SBS"))
  cov.area <- sum(drop_units(tmp2 %>% st_area()*1e-6))
  cov.prop <- cov.area/tmp$Area_km2
  cov.df$SBS_prop[i] <- cov.prop
}


# watercourses layer
# bcdc_search("NTS BC River", res_format = "wms")
aoi.RLW <- retrieve_geodata_aoi(ID = "414be2d6-f4d9-4f32-b960-caa074c6d36b", 
                                aoi= aoi %>% st_transform(crs=3005)) # to download enough area for distance measure
#distance of centroid to watercourse layer might be an option
aoi.RLW %>% dplyr::count(DESCRIPTION) %>% st_drop_geometry()
aoi.RLW <- aoi.RLW %>% filter(DESCRIPTION!="Land")
ggplot()+
  geom_sf(data=aoi)+
  geom_sf(data=aoi.RLW, aes(fill=DESCRIPTION))

rlw.dist <- st_nn(aoi_centroid, aoi.RLW %>% st_transform(crs=26910), k=1, returnDist = T)
cov.df$RLW_dist <- unlist(rlw.dist$dist)
cov.df$RLW_type <- unlist(rlw.dist$nn)
cov.df$RLW_type <- aoi.RLW$DESCRIPTION[match(cov.df$RLW_type,rownames(aoi.RLW))]

# transportation layer (Digital Road Atlas)
# bcdc_search("road", res_format = "wms")
aoi.DRA <- retrieve_geodata_aoi(ID = "bb060417-b6e6-4548-b837-f9060d94743e", aoi=aoi %>% st_transform(crs=3005))
# think about density per cell
aoi.DRA$Length_m <- st_length(aoi.DRA)
aoi.DRA %>% summarise(sum(Length_m)) %>% st_drop_geometry
aoi.DRA %>% count(FEATURE_TYPE) %>% st_drop_geometry
aoi.DRA$Year <- year(aoi.DRA$DATA_CAPTURE_DATE)
#filter to the appropriate year - only roads built the year of or before study
aoi.DRA.97 <- aoi.DRA %>% filter(Year<1998)

cov.df$RD_density97 <- NA
for(i in seq_len(nrow(cov.df))){
  tmp <- aoi %>% st_transform(crs=26910) %>% filter(aoi$grid_id==i)
  tmp2 <- st_intersection(tmp, aoi.DRA.97 %>% st_transform(crs=26910))
  cov.length <- sum(drop_units(tmp2 %>% st_length()*1e-3))
  cov.prop <- cov.length/tmp$Area_km2
  cov.df$RD_density97[i] <- cov.prop
}


aoi.DRA.96 <- aoi.DRA %>% filter(Year<1997)

cov.df$RD_density96 <- NA
for(i in seq_len(nrow(cov.df))){
  tmp <- aoi %>% st_transform(crs=26910) %>% filter(aoi$grid_id==i)
  tmp2 <- st_intersection(tmp, aoi.DRA.96 %>% st_transform(crs=26910))
  cov.length <- sum(drop_units(tmp2 %>% st_length()*1e-3))
  cov.prop <- cov.length/tmp$Area_km2
  cov.df$RD_density96[i] <- cov.prop
}

aoi.DRA.19 <- aoi.DRA %>% filter(Year<2020)

cov.df$RD_density19 <- NA
for(i in seq_len(nrow(cov.df))){
  tmp <- aoi %>% st_transform(crs=26910) %>% filter(aoi$grid_id==i)
  tmp2 <- st_intersection(tmp, aoi.DRA.19 %>% st_transform(crs=26910))
  cov.length <- sum(drop_units(tmp2 %>% st_length()*1e-3))
  cov.prop <- cov.length/tmp$Area_km2
  cov.df$RD_density19[i] <- cov.prop
}

# vegetation data (VRI)
# bcdc_search("VRI", res_format = "wms")
aoi.VRI <- retrieve_geodata_aoi(ID = "2ebb35d8-c82f-4a17-9c96-612ac3532d55", aoi=aoi %>% st_transform(crs=3005))
# for retrospective data, use the 2002 VRI, already clipped to Williston Basin
aoi.VRIh <- aoi.VRI %>% filter(!is.na(PROJ_HEIGHT_1))

# proportion of VRI with projected height >= 20 m
cov.df$TREE20_prop_recent <- NA
for(i in seq_len(nrow(cov.df))){
  tmp <- aoi %>% st_transform(crs=26910) %>% filter(aoi$grid_id==i)
  tmp2 <- st_intersection(tmp, aoi.VRIh %>% filter(PROJ_HEIGHT_1>=20) %>% st_transform(crs=26910))
  cov.area <- sum(drop_units(tmp2 %>% st_area()*1e-6))
  cov.prop <- cov.area/tmp$Area_km2
  cov.df$TREE20_prop_recent[i] <- cov.prop
}

# proportion of VRI with canopy >45%
aoi.VRI$CROWN_CLOSURE_CLASS_CD <- as.numeric(aoi.VRI$CROWN_CLOSURE_CLASS_CD) # for
aoi.VRIc <- aoi.VRI %>% filter(!is.na(CROWN_CLOSURE_CLASS_CD))

# summary(aoi.VRI$CROWN_CLOSURE_CLASS_CD, na.rm=T)
cov.df$CANOPY_prop_recent <- NA
for(i in seq_len(nrow(cov.df))){
  tmp <- aoi %>% st_transform(crs=26910) %>% filter(aoi$grid_id==i)
  tmp2 <- st_intersection(tmp, aoi.VRIc %>% filter(CROWN_CLOSURE>=45) %>% st_transform(crs=26910))
  cov.area <- sum(drop_units(tmp2 %>% st_area()*1e-6))
  cov.prop <- cov.area/tmp$Area_km2
  cov.df$CANOPY_prop_recent[i] <- cov.prop
}

aoi.VRIh$edge <- ifelse(aoi.VRIh$PROJ_HEIGHT_1<3, "EdgeL3",
                        ifelse(aoi.VRIh$PROJ_HEIGHT_1>=3, "EdgeU3", NA))

aoi.VRI.edge <-aoi.VRIh %>% filter(!is.na(edge)) %>% group_by(edge) %>%
  summarise(across(geometry, ~ st_union(.)), .groups = "keep") %>%
  summarise(across(geometry, ~ st_combine(.)))

# aoi.VRI.edge %>% filter() st_cast(aoi.VRI.edge$geometry, "MULTILINESTRING")

# length of forest edge
# first created polygons of VRI heights >=3 m and < 3 m
cov.df$EDGE_density_recent <- NA
for(i in seq_len(nrow(cov.df))){
  tmp <- aoi %>% st_transform(crs=26910) %>% filter(aoi$grid_id==i)
  tmp2 <- st_intersection(tmp, aoi.VRI.edge %>% filter(edge=="EdgeL3") %>% st_transform(crs=26910))
  tmp2 <- st_cast(tmp2, "MULTILINESTRING")
  cov.length <- sum(drop_units(tmp2 %>% st_length()*1e-3))
  cov.prop <- cov.length/tmp$Area_km2
  cov.df$EDGE_density_recent[i] <- cov.prop
}


### for retro data
aoi.VRI <- st_read(dsn="./data", layer="VRI2002_WBaoi")
aoi.VRI$PROJ_HEIGHT_1 <- aoi.VRI$PROJ_HEIGH

summary(aoi.VRI$PROJ_HEIGHT_1)
ggplot()+
  geom_sf(data=aoi)+
  geom_sf(data=aoi.VRI, aes(fill=as.numeric(PROJ_HEIGHT_1)))

aoi.VRI$PROJ_HEIGHT_1 <- as.numeric(aoi.VRI$PROJ_HEIGHT_1)
aoi.VRI %>% filter(!is.na(PROJ_HEIGHT_1)) %>%  
  summarise(mean = mean(PROJ_HEIGHT_1), min = min(PROJ_HEIGHT_1), max=max(PROJ_HEIGHT_1), sd = sd(PROJ_HEIGHT_1)) %>% st_drop_geometry()
#mean   min   max    sd
#19.6   0.2  42.1  9.38
# from literature (Breault et al. 2021) going with
aoi.VRIh <- aoi.VRI %>% filter(!is.na(PROJ_HEIGHT_1))

# proportion of VRI with projected height >= 20 m
cov.df$TREE20_prop_retro <- NA
for(i in seq_len(nrow(cov.df))){
  tmp <- aoi %>% st_transform(crs=26910) %>% filter(aoi$grid_id==i)
  tmp2 <- st_intersection(tmp, aoi.VRIh %>% filter(PROJ_HEIGHT_1>=20) %>% st_transform(crs=26910))
  cov.area <- sum(drop_units(tmp2 %>% st_area()*1e-6))
  cov.prop <- cov.area/tmp$Area_km2
  cov.df$TREE20_prop_retro[i] <- cov.prop
}

# proportion of VRI with canopy >45%
aoi.VRI$CROWN_CLOSURE <- as.numeric(aoi.VRI$CROWN_CLOS) # for retrospective data
aoi.VRIc <- aoi.VRI %>% filter(!is.na(CROWN_CLOSURE))

# summary(aoi.VRI$CROWN_CLOSURE_CLASS_CD, na.rm=T)
cov.df$CANOPY_prop_retro <- NA
for(i in seq_len(nrow(cov.df))){
  tmp <- aoi %>% st_transform(crs=26910) %>% filter(aoi$grid_id==i)
  tmp2 <- st_intersection(tmp, aoi.VRIc %>% filter(CROWN_CLOSURE>=45) %>% st_transform(crs=26910))
  cov.area <- sum(drop_units(tmp2 %>% st_area()*1e-6))
  cov.prop <- cov.area/tmp$Area_km2
  cov.df$CANOPY_prop_retro[i] <- cov.prop
}

aoi.VRIh$edge <- ifelse(aoi.VRIh$PROJ_HEIGHT_1<3, "EdgeL3",
                        ifelse(aoi.VRIh$PROJ_HEIGHT_1>=3, "EdgeU3", NA))

aoi.VRI.edge <-aoi.VRIh %>% filter(!is.na(edge)) %>% group_by(edge) %>%
  summarise(across(geometry, ~ st_union(.)), .groups = "keep") %>%
  summarise(across(geometry, ~ st_combine(.)))

# aoi.VRI.edge %>% filter() st_cast(aoi.VRI.edge$geometry, "MULTILINESTRING")

# length of forest edge
# first created polygons of VRI heights >=3 m and < 3 m
cov.df$EDGE_density_retro <- NA
for(i in seq_len(nrow(cov.df))){
  tmp <- aoi %>% st_transform(crs=26910) %>% filter(aoi$grid_id==i)
  tmp2 <- st_intersection(tmp, aoi.VRI.edge %>% filter(edge=="EdgeL3") %>% st_transform(crs=26910))
  tmp2 <- st_cast(tmp2, "MULTILINESTRING")
  cov.length <- sum(drop_units(tmp2 %>% st_length()*1e-3))
  cov.prop <- cov.length/tmp$Area_km2
  cov.df$EDGE_density_retro[i] <- cov.prop
}

# cutbock (Consolidated Cutblocks)
# bcdc_search("cutblock", res_format = "wms")
aoi.CUT <- retrieve_geodata_aoi(ID = "b1b647a6-f271-42e0-9cd0-89ec24bce9f7", aoi=aoi %>% st_transform(crs=3005))
# as.data.frame(aoi.CUT %>% group_by(HARVEST_YEAR) %>% summarise(Area_km2=sum(Area_km2)) %>% st_drop_geometry())
# 198/823 # 24% cut in 1996
# 343/823 # 41% cut in 2022

# ggplot()+
#   geom_sf(data=aoi)+
#   geom_sf(data=aoi.CUT %>% filter(HARVEST_YEAR < 1997), aes(fill=HARVEST_YEAR))
#proportion of cell harvested

aoi.CUT96 <- aoi.CUT %>% filter(HARVEST_YEAR < 1997)

cov.df$HARVEST_prop96 <- NA
for(i in seq_len(nrow(cov.df))){
  tmp <- aoi %>% st_transform(crs=26910) %>% filter(aoi$grid_id==i)
  tmp2 <- st_intersection(tmp, aoi.CUT96 %>% st_transform(crs=26910))
  cov.area <- sum(drop_units(tmp2 %>% st_area()*1e-6))
  cov.prop <- cov.area/tmp$Area_km2
  cov.df$HARVEST_prop96[i] <- cov.prop
}

aoi.CUT97 <- aoi.CUT %>% filter(HARVEST_YEAR < 1998)

cov.df$HARVEST_prop97 <- NA
for(i in seq_len(nrow(cov.df))){
  tmp <- aoi %>% st_transform(crs=26910) %>% filter(aoi$grid_id==i)
  tmp2 <- st_intersection(tmp, aoi.CUT97 %>% st_transform(crs=26910))
  cov.area <- sum(drop_units(tmp2 %>% st_area()*1e-6))
  cov.prop <- cov.area/tmp$Area_km2
  cov.df$HARVEST_prop97[i] <- cov.prop
}

aoi.CUT19 <- aoi.CUT %>% filter(HARVEST_YEAR < 2020)

cov.df$HARVEST_prop19 <- NA
for(i in seq_len(nrow(cov.df))){
  tmp <- aoi %>% st_transform(crs=26910) %>% filter(aoi$grid_id==i)
  tmp2 <- st_intersection(tmp, aoi.CUT19 %>% st_transform(crs=26910))
  cov.area <- sum(drop_units(tmp2 %>% st_area()*1e-6))
  cov.prop <- cov.area/tmp$Area_km2
  cov.df$HARVEST_prop19[i] <- cov.prop
}

### traplines
aoi.TRP <- st_read(dsn="./data", layer="Traplines_in_study_area")
aoi.TRP <- aoi.TRP %>% arrange(TRAPLINE_1)
aoi.TRP$TRP_DNSTY_9195 <- c(2.00, 0.82, 16.89, 78.77)
aoi.TRP$TRP_DNSTY_9296 <- c(1.58, 0.70, 25.89, 79.16)
aoi.TRP$TRP_DNSTY_1418 <- c(0.59, 2.09, 0.00, 0.00)

ggplot()+
  # geom_sf(data=aoi.TRP, aes(fill=as.factor(TRP_DNSTY_9195)))+
  geom_sf(data=aoi.TRP, aes(fill=as.factor(TRP_DNSTY_9296)))+
  # geom_sf(data=aoi.TRP, aes(fill=as.factor(TRP_DNSTY_1418)))+
  geom_sf(data=aoi_grid, fill=NA)

tmp <- st_join(aoi_grid %>% st_transform(crs=26910), 
               aoi.TRP %>% st_transform(crs=26910) %>% select(TRP_DNSTY_9195, TRP_DNSTY_9296, TRP_DNSTY_1418), largest=TRUE)


cov.df <- left_join(cov.df, tmp %>% select(grid_id, TRP_DNSTY_9195, TRP_DNSTY_9296, TRP_DNSTY_1418) %>% st_drop_geometry())


ggplot()+
  # geom_sf(data=tmp, aes(fill=TRP_DNSTY_1418))+
  geom_sf(data=tmp, aes(fill=TRP_DNSTY_9296))+
  # geom_sf(data=tmp, aes(fill=TRP_DNSTY_9195))+
  geom_sf(data=retro.traps.out[[1]]$traps.sf, col="black")+
  geom_sf(data=retro.traps.out[[2]]$traps.sf, col="red")+
  geom_sf(data=recent_traps.sf, col="blue")

# wildfire (Fire Perimeters)
# bcdc_search("fire", res_format = "wms")
# no wildfires

################################################################################
# write.csv(cov.df,"data/ALL.covdata.csv")
cov.df <- read.csv("data/ALL.covdata.csv", row.names = 1)
head(cov.df)
cov.df$retro_96 <- cov.df$retro_97 <- cov.df$recent_19 <- NA
cov.df$retro_96 <- case_when(cov.df$grid_id %in% unique(retro_96$grid_id) ~ 1)
cov.df$retro_97 <- case_when(cov.df$grid_id %in% unique(retro_97$grid_id) ~ 1)
cov.df$recent_19 <- case_when(cov.df$grid_id %in% unique(recent_19$grid_id) ~ 1)

summary(cov.df)
cov.df %>% filter(RLW_dist==0) %>% count(RLW_type)
cov.df$use <- 1
cov.df$use <- case_when(cov.df$RLW_dist==0 & cov.df$RLW_type!="Man-made waterbody" ~ 1,
                        cov.df$RLW_dist==0 & cov.df$RLW_type=="Man-made waterbody" ~ 0,
                        TRUE ~ as.numeric(cov.df$use))

sum(cov.df$recent_19, na.rm=T)
sum(cov.df$retro_96, na.rm=T)
sum(cov.df$retro_97, na.rm=T)

aoi_grid$use <- cov.df$use
aoi_grid$recent_19 <- cov.df$recent_19
aoi_grid$retro_96 <- cov.df$retro_96
aoi_grid$retro_97 <- cov.df$retro_97

aoi_grid <- aoi_grid %>% rowwise() %>% mutate(Yrs_surveyed = sum(recent_19,retro_96,retro_97, na.rm = T))

pal = pnw_palette(name="Cascades",n=3,type="discrete")

Cairo(file="out/marten_YrsSurveyed_map_969719.PNG",type="png",width=3400,height=2400,pointsize=14,bg="white",dpi=300)
ggplot()+
  geom_sf(data=aoi_grid %>% filter(use==1 | Yrs_surveyed>0), fill=NA)+
  geom_sf(data=aoi_grid %>% filter(use==1 | Yrs_surveyed>0) %>% filter(Yrs_surveyed >0), aes(fill=as.factor(Yrs_surveyed)))+
  theme(legend.position="bottom")+
  theme(legend.title=element_blank())+
  scale_fill_manual(values=pal)
dev.off()

grid_touse <- aoi_grid %>% filter(use==1 | Yrs_surveyed>0) %>% dplyr::select(grid_id) %>% st_drop_geometry()
cov.df$grid_touse <- case_when(cov.df$grid_id %in% grid_touse$grid_id ~ 1)
# cov.df[is.na(cov.df)] <- 0
summary(cov.df)

write.csv(cov.df,"data/covdata.csv")
save(aoi_grid, file=paste0("data/aoi_grid.RData"))


############################--- RETROSPECTIVE DATA ---##########################
# load data if not running concurrently (use run_all.R to source)

out.files

# new grouping of traps to keep consistent grid cells with recent grid cell
trap_grp <- list(retro_96, retro_97)

retro.data.out <- list()

for(r in 1:length(trap_grp)){
  load(paste0("./out/",out.files[r]))
  # marten.data$observations
  # marten.data$trap.oper
  trap.oper <- marten.data$trap.oper
  traps.open <- rowSums(trap.oper)
  
  # add week, and 21 day, and respective occasions to the daylookup
  daylookup <- marten.data$daylookup
  # glimpse(daylookup)
  daylookup$Week <- week(daylookup$Date)
  num.weeks <- round(nrow(daylookup)/7+1)
  week.occ <- rep(1:num.weeks, each = 7)
  daylookup$Occ_week <- week.occ[1:nrow(daylookup)]
  num.21days <- round(nrow(daylookup)/21+1)
  days21.occ <- rep(1:num.21days, each = 21)
  daylookup$Occ_21day <- days21.occ[1:nrow(daylookup)]

  # create a covariate and dataset for weekly occasions
  # covariate is number of days trap open per week (0-7)
  
  # make sure to only use full weeks in occasion week (i.e., remove last week if not containing 7 days)
  weeks.to.use <- daylookup %>% count(Occ_week)
  weeks.to.use <- weeks.to.use[weeks.to.use$n==7,]$Occ_week
  
  days21.to.use <- daylookup %>% count(Occ_21day)
  days21.to.use <- days21.to.use[days21.to.use$n==21,]$Occ_21day
  
  # observation covariates need to be in dim(M,J) where M = number of sites, J = number of sampling occasions
  
  # first need to group traps
  traps <- marten.data$traps*1000
  traps.sf <- st_as_sf(traps, coords=c("x","y"), crs=26910)
  traps.grouped <- trap_grp[[r]]
  colnames(traps.grouped)[1] <- "Trap_Grp"
  traps.grouped$TrapNum <- rownames(traps.grouped)
  num.trp.grps <- length(unique(traps.grouped$Trap_Grp))
  traps.sf$grid_id <- traps.grouped$Trap_Grp
  
  # now back to effort and observations
  week.effort <- as.data.frame(array(NA, dim=c(num.trp.grps, length(weeks.to.use))))
  colnames(week.effort) <- c(weeks.to.use)
  rownames(week.effort) <- unique(sort(traps.grouped$Trap_Grp))
  
  # create trap.oper based on trap groups
  tg_trap.oper <- as.data.frame(trap.oper)
  tg_trap.oper$Trap_Grp <- traps.grouped$Trap_Grp[match(rownames(tg_trap.oper),traps.grouped$TrapNum)]
  tg_trap.oper <- as.data.frame(tg_trap.oper %>% group_by(Trap_Grp) %>% summarise(across(everything(), sum)))
  rownames(tg_trap.oper) <- tg_trap.oper[,1]
  tg_trap.oper <- as.matrix(tg_trap.oper[2:ncol(tg_trap.oper)])

  eff.col <- 1
  for(i in 1:length(week.effort)){
    # i=8
    tmp1 <- as.data.frame(t(tg_trap.oper))
    tmp1$Occ_week <- daylookup$Occ_week[match(rownames(tmp1), as.character(daylookup$Date))]
    tmp1 <- tmp1 %>% filter(Occ_week %in% weeks.to.use)
    tmp1 %>% count(Occ_week)
    tmp2 <- tmp1 %>% filter(Occ_week==i) %>% colSums()
    
    week.effort[,eff.col] <- tmp2[1:nrow(tg_trap.oper)]
    
    eff.col <- eff.col + 1
  }
  
  # trying to find equal 21 day periods for effort - omit the last few days/weeks not in a 21 day period
  tmp <- floor(length(week.effort)/3)
  tmp2 <- floor(tmp*3/3)
  tmp3 <- tmp2*3
  week.effort.trunc <- week.effort[1:tmp3]
  effort.21days <- sapply(seq(1,tmp3,by=3),function(i) rowSums(week.effort.trunc[,i:(i+2)]))
  
  tot.week.effort <- rowSums(week.effort) # total effort per trap 
  tot.21day.effort <- rowSums(effort.21days) # total effort per trap 
  
  # tot.week.effort == tot.21day.effort # not always true, must be due to removing the last few weeks
  
  # have a covariate (obs.cov) of trap effort per week
  # bin occasions from day to week to emulate Poisson data for nmix and better work in occupancy
  
  # nmix.create.y.obs.cov <- function(num.days = num.days){
  # i=1
  traps.open <- rowSums(trap.oper)
  # observations %>% filter(TrapNumber %in% tmp2) 
  
  # need to find the observation data that corresponds to those dates
  # add in the 1 when a marten was in a trap and a 0 if trap open but not in trap
  # create a matrix appending rows as the loop cycles...not sure about this bit
  
 
  observations <- marten.data$observations
  observations$Trap_Grp <- traps.grouped$Trap_Grp[match(observations$TrapNumber, traps.grouped$TrapNum)]
  observations$Count <- 1
  observations <- observations %>% group_by(Trap_Grp, Date_obs) %>% summarise(Count = sum(Count))
  
  observations <- observations %>% arrange(Trap_Grp, Date_obs)
  obs.wide <- pivot_wider(observations, names_from = Trap_Grp, values_from = Count, values_fill = 0)
  obs.wide <- obs.wide %>%  arrange(Date_obs) %>% rename(Date="Date_obs")
  dim(obs.wide) # 
  duplicated(obs.wide$Date) # all should be false
  # need to add in dates and sites
  # add in dates, transpose to add in sites and then transpose back to become the y matrix
  
  # create y for all sites and all occasions 
  # y = R (number of sites/traps) x J (number of sampling periods) with y as repeated counts
  y_all <- left_join(daylookup, obs.wide)
  y_all[is.na(y_all)] <- 0
  y_allT <- as.data.frame(t(y_all %>% dplyr::select(-c(YDay, Date, Occ, Week, Occ_week, Occ_21day))))
  dim(y_allT); sum(y_allT[1:dim(y_all)[1]])
  
  y_allT$Trap_Grp <- colnames(obs.wide[2:ncol(obs.wide)])
  
  tmp1 <- as.data.frame(array(NA, dim=c(num.trp.grps, 0)))
  tmp1$Trap_Grp <- as.character(unique(sort(traps.grouped$Trap_Grp)))
  tmp2 <- left_join(tmp1, y_allT)
  y_day <- tmp2[,2:dim(y_all)[1]]
  rownames(y_day) <- tmp2$Trap_Grp
  colnames(y_day) <- seq_len(ncol(y_day))
  y_day[is.na(y_day)] <- 0
  sum(y_day)
  
  # create y for all sites, weekly, 21 day occasions 
  tmp3 <- as.data.frame(t(y_day))
  rownames(tmp3) <- seq_len(dim(y_day)[2])
  tmp3$Occ_week <- daylookup$Occ_week[match(rownames(tmp3), daylookup$Occ)]
  y_week <- tmp3 %>% group_by(Occ_week) %>% summarise_at(1:num.trp.grps, sum)
  y_week <- as.matrix(t(y_week %>% dplyr::select(-Occ_week)))
  y_week <- y_week[,1:length(weeks.to.use)]
  sum(y_week)
  
  tmp3$Occ_21day <- daylookup$Occ_21day[match(rownames(tmp3), daylookup$Occ)]
  y_21day <- tmp3 %>% group_by(Occ_21day) %>% summarise_at(1:num.trp.grps, sum)
  y_21day <- as.matrix(t(y_21day %>% dplyr::select(-Occ_21day)))
  y_21day <- y_21day[,1:length(days21.to.use)]
  sum(y_day)  # 273 observations 
  
  sum(y_week) == sum(y_day) # same number of observations in weekly and full dataset
  # may not be the same because of omitted days
  
  retro.data.out[[r]] <- list(y_week=y_week, y_day=y_day, y_21day=y_21day, week.effort=week.effort, weeks.to.use=weeks.to.use, 
                              daylookup=daylookup, effort.21days=effort.21days, days21.to.use=days21.to.use,traps.sf=traps.sf)
  
}

save(retro.data.out, file = paste0("./out/retro.data.out9697.RData"))

# Look at detections per week and per trap
# rowSums(retro.data.out[[1]]$y_21day);colSums(retro.data.out[[1]]$y_21day)
# rowSums(retro.data.out[[2]]$y_21day);colSums(retro.data.out[[1]]$y_21day)

# retro_year <- c("1996", "1997", "1998", "1999")
# for(i in 1:4){
# Cairo(file=paste0("out/weekly_det_",retro_year[i],".PNG"),type="png",width=2000,height=2000,pointsize=15,bg="white",dpi=300)
# plot(colSums(retro.data.out[[i]]$y_week), ylab="Weeky marten capture count", xlab="Week", 
#      main=paste0("Live Capture Data - ",retro_year[i],"/",as.numeric(retro_year[i])-1899))
# dev.off()
# }

# 1999/00 detections high for first 3 weeks and then drop off - unlikely for models to converge
# different sampling effort in 1999/00 - concentrating only on fisher and moving traps to target recapturing fisher
# 1998/99 might also be worth ignoring as trapping effort became much more focused on fisher

#####################################################################################
############################--- RETROSPECTIVE DATA ---##########################
# load data if not running concurrently (use run_all.R to source)
# load("out/retro.data.out.RData")
# load("out/recent_occ_data.RData")
###--- DECISION 2022-Mar-03 go with grid cells for recent data as naming convention is fisher row/column
###--- UPDATE - 2022-Mar-04 changed mind to go with hexagons, same as retro data for comparison and to consider marten (not fisher) home range size

hsdat <- recent_occ_data[[1]]
traps.df <- recent_occ_data[[2]]

ltraps.sf <- st_as_sf(traps.df, coords=c("Easting","Northing"), crs=26910)
# ggplot()+
#   geom_sf(data=ltraps.sf, aes(fill=Grid_Focus, col=Grid_Focus))
# 
# ggplot()+
#   geom_sf(data=ltraps.sf, aes(fill=Grid, col=Grid))
# 
# ggplot()+
#   geom_sf(data=aoi_grid)+
#   geom_sf(data=ltraps.sf, aes(fill=Grid, col=Grid))
  
ltraps.sf <- st_join(ltraps.sf, aoi_grid %>% dplyr::select(grid_id), left=TRUE, largest=TRUE)
ggplot()+
  geom_sf(data=ltraps.sf, aes(fill=grid_id, col=grid_id))

hsdat$grid_id <- ltraps.sf$grid_id[match(hsdat$`Sample Station Label`, ltraps.sf$Station)]
hsdat %>% count(`Sampling Session`)
hsdat.mart <- hsdat
hsdat.mart$Occ <- as.factor(str_sub(hsdat.mart$`Sampling Session`,-1))
hsdat.mart <- hsdat.mart %>% dplyr::select(-`Study Area Name`, -Species, -`Sampling Session`) %>% rename("Station"=1, "Animal_ID"=3)
hsdat.mart$Date <- ymd(hsdat.mart$Date)
as.data.frame(hsdat.mart %>% group_by(grid_id) %>% count(Occ))
count(hsdat.mart, Station) # 44 stations with detections
count(hsdat.mart, grid_id) #32 marten and fisher focused grid cells with detections (occassions grouped)

hsdat.mart %>% arrange(grid_id) %>% count(grid_id, Occ)
hsdat.mart <- hsdat.mart %>% arrange(grid_id, Date, Occ)
hsdat.mart$Count <- 1
observations <- hsdat.mart %>% group_by(Station, Occ) %>% summarise(Count=sum(Count))
obs.wide <- pivot_wider(observations, names_from = Occ, values_from = Count, values_fill = 0)
obs.wide <- obs.wide %>%  arrange(Station)
duplicated(obs.wide$Station) # all should be false

traps.obs <- left_join(ltraps.sf %>% dplyr::select(Station, grid_id) %>% st_drop_geometry(), obs.wide)
traps.obs[is.na(traps.obs)] <- 0
colnames(traps.obs)[3:6] <- c("Occ3","Occ2","Occ4","Occ1") 

rec_ydata <- traps.obs %>% dplyr::select(-Station) %>% group_by(grid_id) %>% summarise_at(vars("Occ3":"Occ1"), sum, na.rm=TRUE)
rec_ydata <- as.data.frame(rec_ydata)
row.names(rec_ydata) <- rec_ydata$grid_id
rec_ydata <- rec_ydata %>% dplyr::select(-grid_id)
rec_ydata <- rec_ydata[c("Occ1","Occ2","Occ3","Occ4")]
rec_ydata <- as.matrix(rec_ydata)
dim(rec_ydata); sum(rec_ydata)

traps.obs$effCount <- 1
rec_effort <- as.data.frame(traps.obs %>% group_by(grid_id) %>% dplyr::count(effCount) %>% dplyr::select(-effCount))
rec_effort$Occ4 <- rec_effort$Occ3 <- rec_effort$Occ2 <- rec_effort$Occ1 <- rec_effort$n
row.names(rec_effort) <- rec_effort$grid_id
rec_effort$n <- rec_effort$grid_id <- NULL
rec_effort <- as.matrix(rec_effort)

###--- SINCE GROUPING COV DATA WITH HEXAGONS, NEED TO GROUP OBS DATA WITH HEX TOO ---###
rec.data.out <- list(rec_effort=rec_effort, rec_ydata=rec_ydata)
save(rec.data.out, file = paste0("./out/rec.data.out.RData"))

aoi_grid %>% count(recent_19, na.rm=T) %>% st_drop_geometry(); nrow(rec_ydata)
aoi_grid %>% count(retro_96, na.rm=T) %>% st_drop_geometry(); length(unique(retro_96$grid_id))
aoi_grid %>% count(retro_97, na.rm=T) %>% st_drop_geometry(); length(unique(retro_97$grid_id))

################################################################################
save.image("data/03_analysis_prep.RData")
# load("data/03_analysis_prep.RData")
################################################################################