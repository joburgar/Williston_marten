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
# 01_load.R
# script to load elk collar and EPU data
# written by Joanna Burgar (Joanna.Burgar@gov.bc.ca) - 17-Feb-2021
#####################################################################################

.libPaths("C:/Program Files/R/R-4.0.5/library") # to ensure reading/writing libraries from C drive
tz = Sys.timezone() # specify timezone in BC

# Load Packages
list.of.packages <- c("tidyverse", "lubridate","chron","bcdata", "bcmaps","sf", "rgdal", "readxl", "Cairo", "coda",
                      "OpenStreetMap", "ggmap", "truncnorm", "doParallel", "nimble", "scrbook", "xtable")

# Check you have them and load them
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)
#####################################################################################

###--- function to retrieve geodata from BCGW

retrieve_geodata_aoi <- function (ID=ID){
  aoi.geodata <- bcdc_query_geodata(ID) %>%
    filter(BBOX(st_bbox(aoi))) %>%
    collect()
  aoi.geodata <- aoi.geodata %>% st_intersection(aoi)
  aoi.geodata$Area_km2 <- st_area(aoi.geodata)*1e-6
  aoi.geodata <- drop_units(aoi.geodata)
  return(aoi.geodata)
}

#################################################################################

############################--- RETROSPECTIVE DATA ---###########################
###--- load marten live trap trap data
lttrap96 <- read_excel("data/WFI_marten_trapping_captures_211006.xlsx",
                  sheet = "TDF 96-97", na="NULL", col_types="text") %>%  type_convert()
as.Date(colnames(lttrap96[,ncol(lttrap96)])) - as.Date(colnames(lttrap96[,4])) # 177 days


lttrap97 <- read_excel("data/WFI_marten_trapping_captures_211006.xlsx",
                       sheet = "TDF 97-98", na="NULL", col_types="text") %>%  type_convert()
as.Date(colnames(lttrap97[,ncol(lttrap97)])) - as.Date(colnames(lttrap97[,4])) # 195 days


lttrap98 <- read_excel("data/WFI_marten_trapping_captures_211006.xlsx",
                       sheet = "TDF 98-99", na="NULL", col_types="text") %>%  type_convert()
as.Date(colnames(lttrap98[,ncol(lttrap98)])) - as.Date(colnames(lttrap98[,4])) # 158 days

lttrap99 <- read_excel("data/WFI_marten_trapping_captures_211006.xlsx",
                       sheet = "TDF 99-00", na="NULL", col_types="text") %>%  type_convert()
as.Date(colnames(lttrap99[,ncol(lttrap99)])) - as.Date(colnames(lttrap99[,4])) # 132 days


###--- load marten live trap detection data (2000)
ltdat <- read_excel("data/WFI_marten_trapping_captures_211006.xlsx",
                   sheet = "Marten captures", na="NULL", col_types="text") %>%  type_convert()

###--- load marten hair snag trap data
hstrap <- read_excel("data/Williston_Fisher_Marten_SPI_submission_210929.xlsm",
                       sheet = "Sample Station Information", na="NULL", col_types="text") %>%  
  type_convert() %>% select(`Sample Station Label`, `UTM Zone Sample Station`, `Easting Sample Station`, `Northing Sample Station`, `Grid Cell`)
glimpse(hstrap)
summary(hstrap)

############################--- CURRENT DATA ---###########################
###--- load marten hair snag data (2020)
hsdat <- read_excel("data/Williston_Fisher_Marten_SPI_submission_210929.xlsm",
                    sheet = "Biological Sample Collection", na="NULL", col_types="text") %>%
  type_convert() %>% filter(Species=="M-MAAM") %>% select(`Study Area Name`, `Sample Station Label`, Date, `Species`,`Animal ID`,`Comments`, `Sex`, `Sampling Session`)
glimpse(hsdat)

hsdat <- hsdat %>% filter(Comments=="Sample genotyped")
recaps <- hsdat %>% count(`Animal ID`)
recaps %>% summarise(mean(n), min(n), max(n), se=sd(n)/sqrt(nrow(recaps)))
recaps %>% filter(n==1) # 30 animals only caught once
recaps %>% filter(n==2) # 4 animals caught twice
recaps %>% filter(n==3) # 3 animals caught three times
# 