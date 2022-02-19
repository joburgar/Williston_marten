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
list.of.packages <- c("tidyverse", "lubridate","chron","sf","Cairo", "sf", "nngeo", "units","OpenStreetMap", "ggmap")

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
#- function to retrieve geodata from BCGW
retrieve_geodata_aoi <- function (ID=ID){
  aoi.geodata <- bcdc_query_geodata(ID) %>%
    filter(BBOX(st_bbox(aoi))) %>%
    collect()
  aoi.geodata <- aoi.geodata %>% st_intersection(aoi)
  aoi.geodata$Area_km2 <- st_area(aoi.geodata)*1e-6
  aoi.geodata <- drop_units(aoi.geodata)
  return(aoi.geodata)
}

#####################################################################################
############################--- RETROSPECTIVE DATA ---##########################
# load data if not running concurrently (use run_all.R to source)
# load("out/MartenData_1996.RData")
# nm.area <- diff(marten.data$xlim)*diff(marten.data$ylim)/100	# Density reported per 100 sq km
# nm.area # 51.21623 100 km2 or 5122 km2 total area (is this true?)

load("out/retro_traps_grp.RData")
str(traps.grid)

out.files <- list.files("./out/", pattern="*.RData")
out.files <- out.files[grepl("MartenData", out.files)]
out.files <- out.files[!grepl("2020", out.files)]

# retro.data.out <- vector('list', length(out.files))
# str(retro.data.out)

retro.data.out <- list()

for(r in 1:length(out.files)){
  load(paste0("./out/",out.files[r]))
  marten.data$observations
  marten.data$trap.oper
  trap.oper <- marten.data$trap.oper
  traps.open <- rowSums(trap.oper)
  traps.open[order(traps.open)]
  
  # add week, and 21 day, and respective occasions to the daylookup
  daylookup <- marten.data$daylookup
  glimpse(daylookup)
  daylookup$Week <- week(daylookup$Date)
  num.weeks <- round(nrow(daylookup)/7+1)
  week.occ <- rep(1:num.weeks, each = 7)
  daylookup$Occ_week <- week.occ[1:nrow(daylookup)]
  num.21days <- round(nrow(daylookup)/21+1)
  days21.occ <- rep(1:num.21days, each = 21)
  length(days21.occ)
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
  grid_output <- find_grid(input=traps.sf, cellsize=5000)
  traps.grouped <- grid_output$aoi_utm %>% st_drop_geometry()
  traps.grouped$TrapNum <- rownames(traps.grouped)
  num.trp.grps <- length(unique(traps.grouped$Trap_Grp))
  
  # now back to effort and observations
  week.effort <- as.data.frame(array(NA, dim=c(num.trp.grps, length(weeks.to.use))))
  colnames(week.effort) <- weeks.to.use
  rownames(week.effort) <- rownames(unique(traps.grouped$Trap_Grp))
  
  # create trap.oper based on trap groups
  tg_trap.oper <- as.data.frame(trap.oper)
  tg_trap.oper$Trap_Grp <- traps.grouped$Trap_Grp[match(rownames(tg_trap.oper),traps.grouped$TrapNum)]
  tg_trap.oper <- as.data.frame(tg_trap.oper %>% group_by(Trap_Grp) %>% summarise(across(everything(), sum)))
  rownames(tg_trap.oper) <- tg_trap.oper[,1]
  tg_trap.oper <- as.matrix(tg_trap.oper[2:ncol(tg_trap.oper)])
  dim(tg_trap.oper)
  
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
  dim(y_allT)
  
  y_allT$Trap_Grp <- colnames(obs.wide[2:ncol(obs.wide)])
  
  tmp1 <- as.data.frame(array(NA, dim=c(num.trp.grps, 0)))
  tmp1$Trap_Grp <- rownames(tmp1)
  tmp2 <- left_join(tmp1, y_allT)
  y_day <- as.matrix(tmp2[,2:ncol(tmp2)])
  colnames(y_day) <- seq_len(ncol(trap.oper))
  rownames(y_day) <- seq_len(num.trp.grps)
  y_day[is.na(y_day)] <- 0
  
  # create y for all sites, weekly, 21 day occasions 
  tmp3 <- as.data.frame(t(y_day))
  tmp3$Occ_week <- daylookup$Occ_week[match(rownames(tmp3), daylookup$Occ)]
  y_week <- tmp3 %>% group_by(Occ_week) %>% summarise_at(1:num.trp.grps, sum)
  y_week <- as.matrix(t(y_week %>% dplyr::select(-Occ_week)))
  y_week <- y_week[,1:length(weeks.to.use)]
  
  tmp3$Occ_21day <- daylookup$Occ_21day[match(rownames(tmp3), daylookup$Occ)]
  y_21day <- tmp3 %>% group_by(Occ_21day) %>% summarise_at(1:num.trp.grps, sum)
  y_21day <- as.matrix(t(y_21day %>% dplyr::select(-Occ_21day)))
  y_21day <- y_21day[,1:length(days21.to.use)]
  
  
  sum(y_21day)# 234 observations
  sum(y_week) # 273 observations
  sum(y_day)  # 273 observations 
  # omitted 3 days from the week occasion because last week not a full week
  
  sum(y_week) == sum(y_day) # same number of observations in weekly and full dataset
  # may not be the same because of omitted days
  
  retro.data.out[[r]] <- list(y_week=y_week, y_day=y_day, y_21day=y_21day, week.effort=week.effort, weeks.to.use=weeks.to.use, 
                              daylookup=daylookup, effort.21days=effort.21days, days21.to.use=days21.to.use, grid_output=grid_output)
  
}


# Look at detections per week and per trap
# glimpse(retro.data.out[[1]])
# rowSums(retro.data.out[[1]]$y_week);colSums(retro.data.out[[1]]$y_week)
# rowSums(retro.data.out[[1]]$y_day);colSums(retro.data.out[[1]]$y_day)
# rowSums(retro.data.out[[2]][[1]]);colSums(retro.data.out[[2]][[1]])
# rowSums(retro.data.out[[3]][[1]]);colSums(retro.data.out[[3]][[1]])
# rowSums(retro.data.out[[4]][[1]]);colSums(retro.data.out[[4]][[1]])

retro_year <- c("1996", "1997", "1998", "1999")
for(i in 1:4){
Cairo(file=paste0("out/weekly_det_",retro_year[i],".PNG"),type="png",width=2000,height=2000,pointsize=15,bg="white",dpi=300)
plot(colSums(retro.data.out[[i]]$y_week), ylab="Weeky marten capture count", xlab="Week", 
     main=paste0("Live Capture Data - ",retro_year[i],"/",as.numeric(retro_year[i])-1899))
dev.off()
}

# 1999/00 detections high for first 3 weeks and then drop off - unlikely for models to converge
# different sampling effort in 1999/00 - concentrating only on fisher and moving traps to target recapturing fisher
# 1998/99 might also be worth ignoring as trapping effort became much more focused on fisher


####################################################################################
### START HERE AGAIN, NOW THAT FISHNET GRID MAKES SENSE
# NEED TO REALLY THINK ABOUT COVARIATES TO BRING IN
###--- grab covariate data
glimpse(retro.data.out[[1]]$grid_output)
aoi <- retro.data.out[[1]]$grid_output$fishnet_grid_sf
traps.sf <- retro.data.out[[1]]$grid_output$aoi_utm

ggplot()+
  geom_sf(data=aoi)+
  geom_sf(data=traps.sf)

#- Elevation
# aoi_raster <- cded_raster(aoi) 
# plot(aoi_raster)
# plot(traps.sf[1] %>% st_transform(crs = 4326), add= TRUE) # as a check
# 
# aoi.cded <- rasterToPoints(aoi_raster) # convert to points for join
# aoi.cded.sf <- st_as_sf(as.data.frame(aoi.cded), coords = c("x","y"), crs = 4326) # create spatial layer
# aoi.cded.utm <- st_transform(aoi.cded.sf, crs = 26910) # convert to utm for join distance
# 
# sa.elev.dist <- st_nn(aoi, aoi.cded.utm, k=1, returnDist = T)
# tmp$elev.dist <- unlist(sa.elev.dist$nn)
# tmp$elev <- aoi.cded.utm$elevation[match(tmp$elev,rownames(aoi.cded.utm))]
# summary(tmp$elev)



# load covariates from bcdata
# using the bc data warehouse option to clip to aoi
aoi <- aoi %>% st_transform(3005)

# biogeoclimatic zones
# bcdc_search("Biogeoclimatic zone", res_format = "wms")
# 3: BEC Map (other, wms, kml)
# ID: f358a53b-ffde-4830-a325-a5a03ff672c3
# Name: bec-map
aoi.BEC <- retrieve_geodata_aoi(ID = "f358a53b-ffde-4830-a325-a5a03ff672c3")
aoi.BEC %>% group_by(MAP_LABEL) %>% summarise(Area_km2=sum(Area_km2)) %>% st_drop_geometry
aoi.BEC %>% group_by(ZONE) %>% summarise(Area_km2=sum(Area_km2)) %>% st_drop_geometry
aoi.BEC %>% summarise(Area_km2=sum(Area_km2)) %>% st_drop_geometry # 823 km2 for study area, at grid cell size
# proportion of SBS in cell? proportion of ESSF # 139 km2 of ESSF and 685 km2 of SBS overall (also <1 km2 of BAFA)
ggplot()+
  geom_sf(data=aoi)+
  geom_sf(data=aoi.BEC, aes(fill=ZONE))

# watercourses layer
# bcdc_search("NTS BC River", res_format = "wms")
aoi.RLW <- retrieve_geodata_aoi(ID = "414be2d6-f4d9-4f32-b960-caa074c6d36b") 
#distance of centroid to watercourse layer might be an option
ggplot()+
  geom_sf(data=aoi)+
  geom_sf(data=aoi.RLW, aes(fill=DESCRIPTION))

# transportation layer (Digital Road Atlas)
# bcdc_search("road", res_format = "wms")
aoi.DRA <- retrieve_geodata_aoi(ID = "bb060417-b6e6-4548-b837-f9060d94743e")
# think about density per cell
aoi.DRA$Length_m <- st_length(aoi.DRA)
aoi.DRA %>% summarise(sum(Length_m)) %>% st_drop_geometry
aoi.DRA %>% count(FEATURE_TYPE) %>% st_drop_geometry
ggplot()+
  geom_sf(data=aoi)+
  geom_sf(data=aoi.DRA, aes(fill=FEATURE_TYPE))

aoi.DRA %>% st_intersection(aoi_grid[,1])
st_geometry(aoi_grid)

# vegetation data (VRI)
# bcdc_search("VRI", res_format = "wms")
aoi.VRI <- retrieve_geodata_aoi(ID = "2ebb35d8-c82f-4a17-9c96-612ac3532d55")
summary(aoi.VRI$PROJ_HEIGHT_1)
ggplot()+
  geom_sf(data=aoi)+
  geom_sf(data=aoi.VRI, aes(fill=PROJ_HEIGHT_1))

aoi.VRI %>% filter(!is.na(PROJ_HEIGHT_1)) %>%  
  summarise(mean = mean(PROJ_HEIGHT_1), min = min(PROJ_HEIGHT_1), max=max(PROJ_HEIGHT_1), sd = sd(PROJ_HEIGHT_1)) %>% st_drop_geometry()
#mean   min   max    sd
#19.6   0.2  42.1  9.38
aoi.VRI$PROJ_HEIGHT_1_cat <- as.factor(ifelse(aoi.VRI$PROJ_HEIGHT_1 < 10, "H0-10",
                                              ifelse(aoi.VRI$PROJ_HEIGHT_1 < 20, "H10-20",
                                                     ifelse(aoi.VRI$PROJ_HEIGHT_1 < 30, "H20-30",
                                                            ifelse(aoi.VRI$PROJ_HEIGHT_1 < 40, "H30-40",
                                                                   ifelse(aoi.VRI$PROJ_HEIGHT_1 < 50, "H40+"))))))# remove NAs
aoi.VRI <- aoi.VRI[complete.cases(aoi.VRI$PROJ_HEIGHT_1),]
ggplot()+
  geom_sf(data = aoi.VRI, aes(fill=PROJ_HEIGHT_1_cat, col=NA)) +
  scale_fill_brewer(palette="Greens") +
  scale_color_brewer(palette="Greens") +
  geom_sf(data = aoi , lwd=1, col="red", fill=NA) +
  theme(legend.title=element_blank())

# cutbock (Consolidated Cutblocks)
# bcdc_search("cutblock", res_format = "wms")
aoi.CUT <- retrieve_geodata_aoi(ID = "b1b647a6-f271-42e0-9cd0-89ec24bce9f7")
as.data.frame(aoi.CUT %>% group_by(HARVEST_YEAR) %>% summarise(Area_km2=sum(Area_km2)) %>% st_drop_geometry())
aoi.CUT %>% filter(HARVEST_YEAR < 1997) %>% summarise(Area_km2=sum(Area_km2)) %>% st_drop_geometry()
aoi.CUT %>% summarise(Area_km2=sum(Area_km2)) %>% st_drop_geometry()
# 198/823 # 24% cut in 1996
# 343/823 # 41% cut in 2022
ggplot()+
  geom_sf(data=aoi)+
  geom_sf(data=aoi.CUT %>% filter(HARVEST_YEAR < 1997), aes(fill=HARVEST_YEAR))
#proportion of cell harvested


# wildfire (Fire Perimeters)
# bcdc_search("fire", res_format = "wms")
# no wildfires