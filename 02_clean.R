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
# 02_clean.R
# script to clean marten trap and detection data from 2000 and 2020
# written by Joanna Burgar (Joanna.Burgar@gov.bc.ca) - 10-Oct-2021
#####################################################################################
version$major
version$minor
R_version <- paste0("R-",version$major,".",version$minor)

.libPaths(paste0("C:/Program Files/R/",R_version,"/library")) # to ensure reading/writing libraries from C drive
tz = Sys.timezone() # specify timezone in BC

# Load Packages
list.of.packages <- c("tidyverse", "lubridate","chron","sf","sp","raster","rgeos","rgdal", "concaveman","Cairo","OpenStreetMap", "ggmap")

# Check you have them and load them
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)
#####################################################################################
###--- load data if not running concurrently
# load("./data/01_load.RData")

############################--- RETROSPECTIVE DATA ---###########################
glimpse(lttrap)
glimpse(ltdat)

### create daylookup for each occasion

# first check how many occasions for each year
# this isn't accurate as 
for(i in 1:length(lttrap)){
  # i=1
num.days <- as.Date(colnames(lttrap[[i]])[ncol(lttrap[[i]])]) - as.Date(colnames(lttrap[[i]][4]))
print(num.days)
}

# Time difference of 177 days
# Time difference of 195 days
# Time difference of 158 days
# Time difference of 132 days

# note that not all traps were open for the entire sampling period
# on some dates no traps were open
# need to create day look up table and fill in missing dates with 0
# create Rdata object for each year of data, with all necessary model input data
for(i in 1:length(lttrap)){
  # i=2
  # rm(i)
  survey.open <- as.Date(colnames(lttrap[[i]][4]))
  survey.close <- as.Date(colnames(lttrap[[i]])[ncol(lttrap[[i]])])
  
  survey.open.year <- year(survey.open)
  
  survey.openJ <- yday(survey.open-ifelse(survey.open.year=="1996",1,0)) # leap year in 1996
  survey.closeJ <- yday(survey.close)
  survey.length <- (365-survey.openJ) + survey.closeJ

  daylookup <- as.data.frame(matrix(0, nrow = survey.length+1, ncol = 2))
  colnames(daylookup) <- c("YDay","Date")
  daylookup$YDay <- c(survey.openJ:365,1:survey.closeJ)
  daylookup$Date <- seq(survey.open, survey.close, by="days")
  daylookup$Occ <- rownames(daylookup)


# check the number of captures for each session
  # ltdat %>% count(`Trapping session`)
  # glimpse(lttrap[i])

###--- wrangle trap data
  traps.df <- as.data.frame(lttrap[i])
  traps.df <- traps.df %>% dplyr::select(Trap_ID:NAD83.northing) %>% rename("Easting"=2, "Northing"=3)
  traps <- data.frame(Site = traps.df$Trap_ID, x = traps.df$Easting, y = traps.df$Northing)
  traps$TrapNumber <- 1:nrow(traps.df) # Index the traps by order 

# create a trap operability matrix
  tmp1 <- as.data.frame(lttrap[[i]])
  tmp1$TrapNumber <- traps$TrapNumber[match(tmp1$Trap_ID,traps$Site)]
  tmp1 <- tmp1[order(tmp1$TrapNumber),] # to ensure order by TrapNumber
  
  tmp2 <- t(tmp1[4:ncol(tmp1)])
  colnames(tmp2) <- tmp1$TrapNumber
  
  trap.oper <- matrix(NA, nrow=nrow(daylookup),ncol=0)
  rownames(trap.oper) <- as.character(daylookup$Date)
  sum(tmp1[,4:ncol(lttrap[[i]])]) # 3362 occasions traps were open
  
  trap.oper <- merge(trap.oper, tmp2, by="row.names", all.x=TRUE)
  trap.oper[is.na(trap.oper)] <- 0
  rownames(trap.oper) <- trap.oper[,1]
  trap.oper <- trap.oper[2:ncol(trap.oper)]
  
  # dim(trap.oper)
  # colSums(trap.oper)
  # rowSums(trap.oper)
  # sum(trap.oper)

  trap.oper <- t(trap.oper) # trap operability matrix with rows as traps (row name = TrapNumber) and columns as Date (occasion)


###--- wrangle detection data
  ltdat %>% count(`Trapping session`)
  ltdat$`Trapping session` <- as.factor(ltdat$`Trapping session`)
  TS.levels <- levels(ltdat$`Trapping session`)
  dat.mart <- ltdat %>% filter(`Trapping session`==TS.levels[i]) %>% 
    dplyr::select(-`Trapping session`, -Status) %>% rename("Date_obs"=`Date checked`)
  head(dat.mart)

# Match trap numbering and dat.mart
  dat.mart$TrapNumber <- traps$TrapNumber[match(dat.mart$Trap_ID, traps$Site)]

# Scale to km
  coord.scale <- 1000
  traps.scale <- traps
  traps.scale <- traps[,c("x", "y")]/coord.scale

# Create a buffer around the traps by 15km.
  pts <- SpatialPoints(traps.scale[,c("x", "y")])
  b.r <- buffer(pts, width = 15, dissolve = TRUE)
  plot(b.r)
  points(pts, col = "red", pch = 4)
  bb <- bbox(b.r)

marten.data <- list(xlim = bb[1,], ylim = bb[2,], trap.oper = trap.oper,
                    traps = traps.scale, observations = dat.mart, daylookup = daylookup)
save(marten.data, file = paste0("./out/MartenData_",survey.open.year,".Rdata"))
}

# load("out/MartenData_1999.Rdata")
# summary(marten.data$observations) # checked all Rda objects = appear correct
glimpse(marten.data)
############################--- CURRENT DATA ---###########################

###--- wrangle trap data
traps.df <- as.data.frame(hstrap)
names(traps.df)

traps.df <- traps.df %>% dplyr::select(-`UTM Zone Sample Station`) %>% rename("Station"=1, "Easting"=2, "Northing"=3, "Grid_Cell"=4)
traps.df %>% count(Grid_Cell) # 106 grid cells
traps <- data.frame(Grid_Cell = traps.df$Grid_Cell, x = traps.df$Easting, y = traps.df$Northing)
traps <- traps[!duplicated(traps$Grid_Cell),]
traps <- traps %>% arrange(Grid_Cell)

traps$Grid_Num <- 1:nrow(traps) # Index the traps by order 

# Add covariate to determine if the trap was in a fisher or marten cell (i.e., cell size varies)
traps.df$tmp <- str_count(traps.df$Grid_Cell, "-")
traps.df$Grid_Focus <- ifelse(traps.df$tmp==3, "marten","fisher")
traps.df$tmp <- NULL
traps.df %>% filter(Grid_Focus=="fisher") %>% count(Grid_Cell)# 66 fisher focused cells
traps.df %>% filter(Grid_Focus=="marten") %>% count(Grid_Cell)# 40 marten focused cells
# traps.df %>% count(Grid_Focus)

traps.sf <- st_as_sf(traps.df, coords=c("Easting","Northing"), crs=26910)
ggplot()+
  geom_sf(data=traps.sf, aes(fill=Grid_Focus, col=Grid_Focus))

# subset data to just the marten focused Grid Cells
traps.marten <- data.frame(Grid_Cell = traps.df[traps.df$Grid_Focus=="marten",]$Grid_Cell, 
                           x = traps.df[traps.df$Grid_Focus=="marten",]$Easting, y = traps.df[traps.df$Grid_Focus=="marten",]$Northing)
traps.marten <- traps.marten[!duplicated(traps.marten$Grid_Cell),]
traps.marten <- traps.marten %>% arrange(Grid_Cell) 
# tmp <- gsub(".*\\-", "", traps.sf$Grid_Cell)
# traps.sf$Grid <- as.factor(substr(traps.sf$Grid_Cell, 1, nchar(traps.sf$Grid_Cell)-nchar(tmp)-1))
# traps.sf %>% count(Grid) %>% st_drop_geometry() # only 


###--- wrangle detection data
hsdat <- hsdat %>% filter(Comments=="Sample genotyped")
# comment out above code if trying with random thinning
# try SCR with thinning with random thinning for unidentified encounters
# https://onlinelibrary.wiley.com/doi/10.1002/ece3.7091
# includes the unidentified but genotyped samples

recaps <- hsdat %>% count(`Animal ID`)
nrow(hsdat) # 47 detections
recaps %>% summarise(mean(n), min(n), max(n), se=sd(n)/sqrt(nrow(recaps)))
nrow(recaps) # 37 individuals detected
# `mean(n)` `min(n)` `max(n)`     se
#   1      1.27        1        3 0.0999
recaps %>% filter(n==1) # 30 animals only caught once
recaps %>% filter(n==2) # 4 animals caught twice
recaps %>% filter(n==3) # 3 animals caught three times

hsdat %>% count(`Sampling Session`)
hsdat.mart <- hsdat
hsdat.mart$Occ <- as.factor(str_sub(hsdat.mart$`Sampling Session`,-1))
hsdat.mart <- hsdat.mart %>% dplyr::select(-`Study Area Name`, -Species, -`Sampling Session`) %>% rename("Station"=1, "Animal_ID"=3)
hsdat.mart$Date <- ymd(hsdat.mart$Date)
hsdat.mart$Grid_Cell <- traps.df$Grid_Cell[match(hsdat.mart$Station, traps.df$Station)]
as.data.frame(hsdat.mart %>% arrange(Grid_Cell, Station))
as.data.frame(hsdat.mart %>% group_by(Grid_Cell, Station) %>% count(Occ))
count(hsdat.mart, Station) # 44 stations with detections
count(hsdat.mart, Grid_Cell) # 42 grid cells with detections

# Create numeric values for animal ID and grid cells
# animal ID
edf <- hsdat.mart %>% dplyr::select(Animal_ID, Occ, Grid_Cell, Sex)
edf <- edf %>% arrange(Animal_ID, Occ, Grid_Cell, Sex)  

animal <- edf[!duplicated(edf$Animal_ID),]
animal <- animal %>% arrange(Animal_ID)
animal$Animal_Num <- row.names(animal)

sex <- animal$Sex
sex <- ifelse(sex=="M",0,1) # 0 is male, 1 is female
# Now put it on the animal:

edf$Animal_Num <- as.numeric(animal$Animal_Num[match(edf$Animal_ID, animal$Animal_ID)])
edf$Sex <- as.numeric(ifelse(edf$Sex=="M",0,1))
edf$Grid_Num <- as.numeric(traps$Grid_Num[match(edf$Grid_Cell, traps$Grid_Cell)])

edf.marten <- edf %>% filter(Grid_Cell %in% traps.marten$Grid_Cell)
unique(edf.marten$Grid_Cell) # 17 Grid Cells with detections (out of a possible 40)
edf.marten %>% count(Animal_ID) # 11 animals detected
edf.marten %>% count(Sex) # 11 male and 8 females

# create detection histories for each marten, by Grid cell rather than station
# keeping code in case need to copy in other scripts, but going with edf code above for this project
# first create larger detection history table
# hsdat.mart$Count <- 1
# observations <- as.data.frame(hsdat.mart  %>% dplyr::select(Animal_ID, Occ, Count, Grid_Cell) %>% 
#   pivot_wider(names_from=Occ, values_from=Count, values_fill = 0))
# observations <- as.data.frame(observations[order(observations$Grid_Cell),])
# observations <- observations[c("Grid_Cell","Animal_ID","1","2","3","4")]
# sum(observations[c("1","2","3","4")]) # 47 detections
# 
# unique.marten <- unique(observations$Animal_ID)
# 
# marten.dh <- array(0, c(37,106,4))
# marten.dh
# dim(marten.dh)
# i=1
# for(i in 1:length(unique.marten)){
#   indi.obs <- matrix(NA, nrow=length(unique(traps.df$Grid_Cell)),ncol=0)
#   rownames(indi.obs) <- as.character(unique(traps.df$Grid_Cell))
#   tmp <- observations %>% filter(Animal_ID==unique.marten[i]) %>% dplyr::select(-Animal_ID)
#   rownames(tmp) <- as.character(tmp$Grid_Cell)
#   tmp$Grid_Cell <- NULL
#   indi.obs <- merge(indi.obs, tmp, by="row.names", all.x=TRUE)
#   indi.obs[is.na(indi.obs)] <- 0
#   indi.obs <- indi.obs[2:ncol(indi.obs)]
#   marten.dh[i,,] <- indi.obs
# }


# Scale to km
coord.scale <- 1000
traps.scale <- traps
traps.scale <- traps[,c("x", "y")]/coord.scale

buffer <- 5 #5 km unit buffer

traps.sc <- as.data.frame(cbind(traps.scale$x-min(traps.scale$x-buffer), traps.scale$y-min(traps.scale$y-buffer)))
colnames(traps.sc) <- c("x","y")

xlim = range(traps.sc[,1])+c(-buffer,buffer)
ylim = range(traps.sc[,2])+c(-buffer,buffer)
area <- diff(xlim)*diff(ylim)/100	# Density reported per 100 sq km

area <- diff(marten.data$xlim)*diff(marten.data$ylim)/100	# Density reported per 100 sq km

747/103.7 # 7.2, 10.3, 15.7
564/103.7 # 5.4
759/103.7 # 7.3



# other data to save
marten.hsdata <- list(J = nrow(traps),
                      area = area,
                      xlim = xlim,
                      ylim = ylim,
                      traps = traps.sc,
                      edf = edf,
                      sex = sex)

save(marten.hsdata, file = paste0("./out/MartenData_2020.Rda"))

# For the marten only trap data
coord.scale <- 1000
traps.scale <- traps.marten[,c("x", "y")]/coord.scale
row.names(traps.scale) <- traps.marten$Grid_Cell

buffer <- 5 #5 km unit buffer

traps.scale.C1 <- traps.scale[1:20,]
traps.scale.C2 <- traps.scale[21:40,]

traps.sc.C1 <- as.data.frame(cbind(traps.scale.C1$x-min(traps.scale.C1$x-buffer), traps.scale.C1$y-min(traps.scale.C1$y-buffer)))
traps.sc.C2 <- as.data.frame(cbind(traps.scale.C2$x-min(traps.scale.C2$x-buffer), traps.scale.C2$y-min(traps.scale.C2$y-buffer)))
colnames(traps.sc.C1) <- c("x","y")
rownames(traps.sc.C1) <- rownames(traps.scale.C1)
colnames(traps.sc.C2) <- c("x","y")
rownames(traps.sc.C2) <- rownames(traps.scale.C2)

xlim.C1 = range(traps.sc.C1[,1])+c(-buffer,buffer)
ylim.C1 = range(traps.sc.C1[,2])+c(-buffer,buffer)
area.C1 <- diff(xlim.C1)*diff(ylim.C1)/100	# Density reported per 100 sq km
area.C1 # 2.92

xlim.C2 = range(traps.sc.C2[,1])+c(-buffer,buffer)
ylim.C2 = range(traps.sc.C2[,2])+c(-buffer,buffer)
area.C2 <- diff(xlim.C2)*diff(ylim.C2)/100	# Density reported per 100 sq km
area.C2 # 2.94

edf.marten <- edf.marten %>% arrange(Animal_ID, Occ, Grid_Cell, Sex)  

animal <- edf.marten[!duplicated(edf.marten$Animal_ID),]
animal <- animal %>% arrange(Animal_ID)
animal$Animal_Num <- row.names(animal)


grid_cell <- as.data.frame(array(NA, dim=c(nrow(traps.scale),0)))
grid_cell$Grid_Cell <- rownames(traps.scale)
grid_cell$Grid_Num <- rownames(grid_cell)

# Now put it on the animal:
edf.marten$Animal_Num <- as.numeric(animal$Animal_Num[match(edf.marten$Animal_ID, animal$Animal_ID)])
edf.marten$Grid_Num <- as.numeric(grid_cell$Grid_Num[match(edf.marten$Grid_Cell, grid_cell$Grid_Cell)])

traps.sc.C1$Grid_Num <- grid_cell$Grid_Num[match(rownames(traps.sc.C1), grid_cell$Grid_Cell)]
traps.sc.C2$Grid_Num <- grid_cell$Grid_Num[match(rownames(traps.sc.C2), grid_cell$Grid_Cell)]

Sex.C1 <- edf.marten %>% filter(Grid_Num<21) %>% dplyr::select(Sex)
Sex.C2 <- edf.marten %>% filter(Grid_Num>20) %>% dplyr::select(Sex)

# Check that traps are in same order as 
martenGrid.hsdata <- list(J = nrow(traps.sc.C1), # 20 traps in each cluster
                      area = list(area.C1, area.C2),
                      xlim = list(xlim.C1, xlim.C2),
                      ylim = list(ylim.C1, ylim.C2),
                      traps = list(traps.sc.C1,traps.sc.C2),
                      edf = list(edf.marten %>% filter(Grid_Num<21),edf.marten %>% filter(Grid_Num>20)),
                      sex = list(Sex.C1$Sex, Sex.C2$Sex))

save(martenGrid.hsdata, file = paste0("./out/MartenGridData_2020.Rda"))


###--- for visualization of recaps

# # Create a buffer around the traps (15 km)
# For visulaization of recaps
# pts <- SpatialPoints(traps.scale[,c("x", "y")])
# b.r <- buffer(pts, width = 15, dissolve = TRUE)
# plot(b.r)
# points(pts, col = "red", pch = 4)
# bb <- bbox(b.r)

# plot captures and recaptures
# hsdat.mart$Easting <- traps.df$Easting[match(hsdat.mart$Station, traps.df$Station)]
# hsdat.mart$Northing <- traps.df$Northing[match(hsdat.mart$Station, traps.df$Station)]
# 
# cap1.coords <- hsdat.mart %>% filter(Animal_ID %in% recaps[recaps$n==1,]$`Animal ID`) %>% dplyr::select("Easting", "Northing")
# cap1.coords.scaled <- cap1.coords/coord.scale
# cap1 <- SpatialPoints(cap1.coords.scaled[,c("Easting", "Northing")])
# 
# cap2.coords <- hsdat.mart %>% filter(Animal_ID %in% recaps[recaps$n==2,]$`Animal ID`) %>% dplyr::select("Easting", "Northing")
# cap2.coords.scaled <- cap2.coords/coord.scale
# cap2 <- SpatialPoints(cap2.coords.scaled[,c("Easting", "Northing")])
# 
# cap3.coords <- hsdat.mart %>% filter(Animal_ID %in% recaps[recaps$n==3,]$`Animal ID`) %>% dplyr::select("Easting", "Northing")
# cap3.coords.scaled <- cap3.coords/coord.scale
# cap3 <- SpatialPoints(cap3.coords.scaled[,c("Easting", "Northing")])
# 
# b.r <- buffer(pts, width = 15, dissolve = TRUE)
# plot(b.r)
# points(pts, col="black", pch=19, cex=0.5)
# points(cap1.coords.scaled, col = "red", pch = 4,  cex=1.5)
# points(cap2.coords.scaled, col = "blue", pch = 4, cex=1.5)
# points(cap3.coords.scaled, col = "green", pch = 4, cex=1.5)

