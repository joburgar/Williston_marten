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
# script to load marten trap and detection data from 2000 and 2020
# written by Joanna Burgar (Joanna.Burgar@gov.bc.ca) - 08-Oct-2021
#####################################################################################
version$major
version$minor
R_version <- paste0("R-",version$major,".",version$minor)

.libPaths(paste0("C:/Program Files/R/",R_version,"/library")) # to ensure reading/writing libraries from C drive
tz = Sys.timezone() # specify timezone in BC

# Load Packages
list.of.packages <- c("tidyverse", "lubridate","chron","bcdata", "bcmaps","sf", "rgdal", "readxl", "Cairo","OpenStreetMap", "ggmap", "truncnorm", "xtable")

# Check you have them and load them
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)
rm(list.of.packages, new.packages) # for housekeeping
#################################################################################

############################--- RETROSPECTIVE DATA ---###########################
###--- load marten live trap trap data
lttrap96 <- read_excel("data/WFI_marten_trapping_captures_211006.xlsx",
                  sheet = "TDF 96-97", na="NULL", col_types="text") %>%  type_convert()

lttrap97 <- read_excel("data/WFI_marten_trapping_captures_211006.xlsx",
                       sheet = "TDF 97-98", na="NULL", col_types="text") %>%  type_convert()

lttrap98 <- read_excel("data/WFI_marten_trapping_captures_211006.xlsx",
                       sheet = "TDF 98-99", na="NULL", col_types="text") %>%  type_convert()

lttrap99 <- read_excel("data/WFI_marten_trapping_captures_211006.xlsx",
                       sheet = "TDF 99-00", na="NULL", col_types="text") %>%  type_convert()

lttrap <- list(lttrap96, lttrap97, lttrap98, lttrap99)
rm(lttrap96, lttrap97, lttrap98, lttrap99)

###--- load marten live trap detection data (2000)
ltdat <- read_excel("data/WFI_marten_trapping_captures_211006.xlsx",
                   sheet = "Marten captures", na="NULL", col_types="text") %>%  type_convert()
# glimpse(ltdat)

############################--- CURRENT DATA ---###########################
###--- load marten hair snag trap data
hstrap <- read_excel("data/Williston_Fisher_Marten_SPI_submission_210929.xlsm",
                     sheet = "Sample Station Information", na="NULL", col_types="text") %>%  
  type_convert() %>% dplyr::select(`Sample Station Label`, `UTM Zone Sample Station`, `Easting Sample Station`, `Northing Sample Station`, `Grid Cell`)
# glimpse(hstrap)

###--- load marten hair snag data (2020)
hsdat <- read_excel("data/Williston_Fisher_Marten_SPI_submission_210929.xlsm",
                    sheet = "Biological Sample Collection", na="NULL") %>%
  filter(Species=="M-MAAM") %>% dplyr::select(`Study Area Name`, `Sample Station Label`, Date, `Species`,`Animal ID`,`Comments`, `Sex`, `Sampling Session`)
# glimpse(hsdat)

################################################################################
save.image("data/01_load.RData")
################################################################################
