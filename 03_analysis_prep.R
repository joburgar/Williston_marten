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
list.of.packages <- c("tidyverse","parallel","unmarked", "nimble","nimbleSCR","MCMCvis","coda","Cairo","basicMCMCplots","tictoc")

# Check you have them and load them
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)
rm(list.of.packages, new.packages) # for housekeeping

################################################################################
# if issues running nimble, run the following line
writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")

################################################################################

# Some of the NIMBLE examples generate html pages with comparisons
# of different MCMCs.  These will be placed in the path contained in
# the variable 'outputDirectory'.  If you do not provide 'outputDirectory',
# the default location will be a subdirectory 'your_comparison_pages'
# of your current working directory.  Your current working directory
# can be obtained by getwd().
if(!exists('outputDirectory')) {
  outputDirectory <- file.path(getwd(), "your_comparison_pages")
  dir.create(outputDirectory, showWarnings = FALSE)
}

############################--- RETROSPECTIVE DATA ---##########################
# Run N-mixutre model for unmarked live trap data
# https://cran.r-project.org/web/packages/nimbleEcology/vignettes/Introduction_to_nimbleEcology.html
# Adapting from https://rdrr.io/cran/nimbleEcology/man/dNmixture.html
# Also need to review Applied Hierarchical Modeling code (converted to NIMBLE and available on github)
# https://github.com/nimble-training/AHMnimble/blob/master/Chapter_6/Section_6p11_setup.R
# unmarked pcount: https://cran.r-project.org/web/packages/unmarked/unmarked.pdf
# https://github.com/nimble-training/AHMnimble/tree/master/Chapter_6

# load data if not running concurrently
# load("out/MartenData_1996.Rdata")
# nm.area <- diff(marten.data$xlim)*diff(marten.data$ylim)/100	# Density reported per 100 sq km
# nm.area # 51.21623 100 km2 or 5122 km2 total area (is this true?)

out.files <- list.files("./out/", pattern="*.Rdata")
out.files <- out.files[!grepl("2020", out.files)]

retro.data.out <- vector('list', length(out.files))
str(retro.data.out)

for(r in 1:length(out.files)){
  
  load(paste0("./out/",out.files[r]))
  # https://github.com/nimble-training/AHMnimble/blob/master/Chapter_6/Section_6p11_setup.R
  marten.data$observations
  marten.data$observations %>% count(Trap_ID)
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
  week.effort <- as.data.frame(array(NA, dim=c(nrow(trap.oper), length(weeks.to.use))))
  colnames(week.effort) <- weeks.to.use
  rownames(week.effort) <- rownames(trap.oper)
  
  eff.col <- 1
  for(i in 1:length(week.effort)){
    # i=8
    tmp1 <- as.data.frame(t(trap.oper))
    tmp1$Occ_week <- daylookup$Occ_week[match(rownames(tmp1), as.character(daylookup$Date))]
    tmp1 <- tmp1 %>% filter(Occ_week %in% weeks.to.use)
    tmp1 %>% count(Occ_week)
    tmp2 <- tmp1 %>% filter(Occ_week==i) %>% colSums()
    
    week.effort[,eff.col] <- tmp2[1:nrow(trap.oper)]
    
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
  
  # for nmixture model, have a covariate (obs.cov) of trap effort per week
  # change occasions from day to week to emulate Poisson data
  
  # nmix.create.y.obs.cov <- function(num.days = num.days){
  # i=1
  traps.open <- rowSums(trap.oper)
  observations <- marten.data$observations
  # observations %>% filter(TrapNumber %in% tmp2) 
  
  # need to find the observation data that corresponds to those dates
  # add in the 1 when a marten was in a trap and a 0 if trap open but not in trap
  # create a matrix appending rows as the loop cycles...not sure about this bit
  
  y <- as.data.frame(array(NA, dim=c(nrow(trap.oper), length(weeks.to.use)))) # to create the y object
  dim(y) # 77 traps by 25 occasions
  
  observations$count <- 1
  observations <- observations %>% arrange(TrapNumber, Date_obs) %>% dplyr::select(-Trap_ID)
  obs.wide <- pivot_wider(observations, names_from = TrapNumber, values_from = count, values_fill = 0)
  obs.wide <- obs.wide %>%  arrange(Date_obs) %>% rename(Date="Date_obs")
  dim(obs.wide) # 96 x 54
  duplicated(obs.wide$Date) # all should be false
  
  # need to add in dates and sites
  # add in dates, transpose to add in sites and then transpose back to become the y matrix
  
  # create y for all sites and all occasions 
  # y = R (number of sites/traps) x J (number of sampling periods) with y as repeated counts
  y_all <- left_join(daylookup, obs.wide)
  y_all[is.na(y_all)] <- 0
  y_allT <- as.data.frame(t(y_all %>% dplyr::select(-c(YDay, Date, Occ, Week, Occ_week, Occ_21day))))
  dim(y_allT)
  
  y_allT$TrapNum <- colnames(obs.wide[2:ncol(obs.wide)])
  
  tmp1 <- as.data.frame(array(NA, dim=c(nrow(trap.oper), 0)))
  tmp1$TrapNum <- rownames(tmp1)
  tmp2 <- left_join(tmp1, y_allT)
  y_day <- as.matrix(tmp2[,2:ncol(tmp2)])
  colnames(y_day) <- seq_len(ncol(trap.oper))
  rownames(y_day) <- seq_len(nrow(trap.oper))
  y_day[is.na(y_day)] <- 0
  
  # create y for all sites, weekly, 21 day occasions 
  tmp3 <- as.data.frame(t(y_day))
  tmp3$Occ_week <- daylookup$Occ_week[match(rownames(tmp3), daylookup$Occ)]
  y_week <- tmp3 %>% group_by(Occ_week) %>% summarise_at(1:nrow(trap.oper), sum)
  y_week <- as.matrix(t(y_week %>% dplyr::select(-Occ_week)))
  y_week <- y_week[,1:length(weeks.to.use)]
  
  tmp3$Occ_21day <- daylookup$Occ_21day[match(rownames(tmp3), daylookup$Occ)]
  y_21day <- tmp3 %>% group_by(Occ_21day) %>% summarise_at(1:nrow(trap.oper), sum)
  y_21day <- as.matrix(t(y_21day %>% dplyr::select(-Occ_21day)))
  y_21day <- y_21day[,1:length(days21.to.use)]
  
  
  
  sum(y_21day)# 277 observations
  sum(y_week) # 279 observations
  sum(y_day)  # 280 observations 
  # omitted 3 days from the week occasion because last week not a full week
  
  sum(y_week) == sum(y_day) # same number of observations in weekly and full dataset
  # not the same because of omitted days
  
  retro.data.out[[r]] <- list(y_week=y_week, y_day=y_day, y_21day=y_21day, week.effort=week.effort, weeks.to.use=weeks.to.use, 
                              daylookup=daylookup, effort.21days=effort.21days, days21.to.use=days21.to.use)
  
}


# Look at detections per week and per trap
glimpse(retro.data.out[[1]])
rowSums(retro.data.out[[1]]$y_week);colSums(retro.data.out[[1]]$y_week)
rowSums(retro.data.out[[1]]$y_day);colSums(retro.data.out[[1]]$y_day)
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
