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

.libPaths("C:/Program Files/R/R-4.0.5/library") # to ensure reading/writing libraries from C drive

# Load Packages
list.of.packages <- c("tidyverse", "lubridate","chron","sf","sp","raster","rgdal", "Cairo","OpenStreetMap", "ggmap")

# Check you have them and load them
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)
#####################################################################################
###--- load data if not running concurrently
getwd()
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
# create Rda object for each year of data, with all necessary SC input data
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
                    traps = traps.scale, observations = dat.mart)
save(marten.data, file = paste0("./out/MartenData_",survey.open.year,".Rda"))
}

# load("out/MartenData_1999.Rda")
# summary(marten.data$observations) # checked all Rda objects = appear correct
############################--- CURRENT DATA ---###########################
glimpse(hstrap)
summary(hstrap)

# hsdat <- hsdat %>% filter(Comments=="Sample genotyped")
# try SCR with thinning with random thinning for unidentified encounters
# https://onlinelibrary.wiley.com/doi/10.1002/ece3.7091
# includes the unidentified but genotyped samples

recaps <- hsdat %>% count(`Animal ID`)
recaps %>% summarise(mean(n), min(n), max(n), se=sd(n)/sqrt(nrow(recaps)))
recaps %>% filter(n==1) # 30 animals only caught once
recaps %>% filter(n==2) # 4 animals caught twice
recaps %>% filter(n==3) # 3 animals caught three times

###--- wrangle trap data
traps.df <- as.data.frame(hstrap)
names(traps.df)
traps.df <- traps.df %>% dplyr::select(-`UTM Zone Sample Station`) %>% rename("Station"=1, "Easting"=2, "Northing"=3, "Grid_Cell"=4)
traps.df %>% count(Grid_Cell) # 106 grid cells
traps <- data.frame(Site = traps.df$Station, x = traps.df$Easting, y = traps.df$Northing)
traps$TrapNumber <- 1:nrow(traps.df) # Index the traps by order 

# 4 sampling occasions 
hstrap %>% count(`Grid Cell`)


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
hsdat %>% count(`Sampling Session`)
hsdat$Occ <- as.factor(str_sub(hsdat$`Sampling Session`,-1))
hsdat <- hsdat %>% dplyr::select(`Sample Station Label`, Date, `Animal ID`, Sex, Occ) %>% rename("Station"=1, "Animal_ID"=3)

hsdat$Grid_Cell <- traps.df$Grid_Cell[match(hsdat$Station, traps.df$Station)]
hsdat %>% count(Grid_Cell)
as.data.frame(hsdat %>% group_by(Grid_Cell, Station) %>% count(Occ))

hsdat.mart <- hsdat %>% 
dat.mart <- ltdat %>% filter(`Trapping session`==ltdat$`Trapping session`[i]) %>% 
  dplyr::select(-`Trapping session`, -Status) %>% rename("Date_obs"=`Date checked`)
head(dat.mart)

# Match trap numbering and dat.mart
hsdat.mart$TrapNumber <- traps$TrapNumber[match(hsdat.mart$Site, traps$Site)]

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

# plot captures and recaptures
hsdat$Easting <- traps.df$Easting[match(hsdat$`Sample Station Label`, traps.df$Station)]
hsdat$Northing <- traps.df$Northing[match(hsdat$`Sample Station Label`, traps.df$Station)]

cap1.coords <- hsdat %>% filter(`Animal ID` %in% recaps[recaps$n==1,]$`Animal ID`) %>% dplyr::select("Easting", "Northing")
cap1.coords.scaled <- cap1.coords/coord.scale
cap1 <- SpatialPoints(cap1.coords.scaled[,c("Easting", "Northing")])

cap2.coords <- hsdat %>% filter(`Animal ID` %in% recaps[recaps$n==2,]$`Animal ID`) %>% dplyr::select("Easting", "Northing")
cap2.coords.scaled <- cap2.coords/coord.scale
cap2 <- SpatialPoints(cap2.coords.scaled[,c("Easting", "Northing")])

cap3.coords <- hsdat %>% filter(`Animal ID` %in% recaps[recaps$n==3,]$`Animal ID`) %>% dplyr::select("Easting", "Northing")
cap3.coords.scaled <- cap3.coords/coord.scale
cap3 <- SpatialPoints(cap3.coords.scaled[,c("Easting", "Northing")])

b.r <- buffer(pts, width = 15, dissolve = TRUE)
plot(b.r)
points(pts, col="black", pch=19, cex=0.5)
points(cap1.coords.scaled, col = "red", pch = 4,  cex=1.5)
points(cap2.coords.scaled, col = "blue", pch = 4, cex=1.5)
points(cap3.coords.scaled, col = "green", pch = 4, cex=1.5)


marten.hsdata <- list(xlim = bb[1,], ylim = bb[2,], trap.oper = trap.oper,
                    traps = traps.scale, observations = hsdat.mart)
save(marten.hsdata, file = paste0("./out/MartenData_2020.Rda"))