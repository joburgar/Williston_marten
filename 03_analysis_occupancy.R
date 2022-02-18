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
# 03_analysis_occupancy.R
# script to run occupancy models fit to live trap (2000) and hair snag (2020) data
# written by Joanna Burgar (Joanna.Burgar@gov.bc.ca) - 17-Feb-2021
#####################################################################################
version$major
version$minor
R_version <- paste0("R-",version$major,".",version$minor)

.libPaths(paste0("C:/Program Files/R/",R_version,"/library")) # to ensure reading/writing libraries from C drive
tz = Sys.timezone() # specify timezone in BC

# Load Packages
list.of.packages <- c("tidyverse","unmarked","Cairo","tictoc")

# Check you have them and load them
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)
rm(list.of.packages, new.packages) # for housekeeping


############################--- SET-UP DATA ---##########################
# Run occupancy models for unmarked and marked live trap data
# consider that male home ranges ~5.25 km2 and female home ranges ~3.16 km2 (Eric Lofroth's MSc thesis)
# larger in fire disturbed landscapes in BC/WA during winter, ~13 km2 for males and ~8.5 km2 for females  (Volkmann & Hodges 2021)

# source("03_analysis_prep.R") # if not working concurrently

# 1999/00 detections high for first 3 weeks and then drop off - unlikely for models to converge
# different sampling effort in 1999/00 - concentrating only on fisher and moving traps to target recapturing fisher
# 1998/99 might also be worth ignoring as trapping effort became much more focused on fisher
#####################################################################################
# borrowed code from https://jamesepaterson.github.io/jamespatersonblog/2020-09-01_occupancyintroduction.html

# Load detection history (100 sites with 10 visits each)
detection_history <- retro.data.out[[1]]$y_21day
detection_history[detection_history > 0] <- 1 # change to 0 and 1 only

# Examine data
# head(detection_history)

# Create unmarkedFrameOccu that holds the data
sample.unmarkedFrame_simple <- unmarkedFrameOccu( # y is a matrix with observed detection history 
  # (0's and 1's, one row per site, one column per survey)
  y = as.matrix(detection_history)) 

# S4 class for occupancy model data
summary(sample.unmarkedFrame_simple)


# Build basic single-season occupancy model with intercepts only (one estimate for detection, one for occupancy)
occu.m1 <- occu(formula = ~1 # detection formula first
                ~1, # occupancy formula second, 
                data = sample.unmarkedFrame_simple)

summary(occu.m1) # Show AIC, estimates (on logit scale), SE, z-scores

# To get real estimate of occupancy (with 95% CI)
predict(occu.m1, 
        newdata = data.frame(site = 1),
        type = "state")
##   Predicted         SE     lower     upper
## 1 0.6209272 0.04861255 0.5221585 0.7105956

# To get real estimate of detection (with 95% CI)
predict(occu.m1, 
        newdata = data.frame(site = 1),
        type = "det")
##   Predicted         SE     lower    upper
## 1  0.478317 0.02018546 0.4389719 0.517933
# Equivalent to inverse logit
boot::inv.logit(coef(occu.m1)[1]) # Real estimate of occupancy
##  psi(Int) 
## 0.6209272
boot::inv.logit(coef(occu.m1)[2]) # Real estimate of detection
##   p(Int) 
## 0.478317

# Load covariate data
effort <- retro.data.out[[1]]$effort.21days

# Build a new unmarkedFramOccu
sample.unmarkedFrame_cov <- unmarkedFrameOccu( # y is a matrix with observed detection history 
  # (0's and 1's, one row per site, one column per survey)
  y = as.matrix(detection_history),
  # obsCovs = observation covariates in a list, 
  # each variable has site rows x survey columns
  obsCovs = list(effort = effort))

# siteCovs = dataframe with site rows x column variables
# siteCovs = site_cov) 

# S4 class for occupancy model data
summary(sample.unmarkedFrame_cov)

occu.m2 <- occu(formula = ~effort # detection formula first
                ~1,
                # ~forest + agri, # occupancy formula second,
                data = sample.unmarkedFrame_cov)

# Summarize
summary(occu.m1)
summary(occu.m2)
