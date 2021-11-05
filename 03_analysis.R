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
# 03_analysis.R
# script to run SCR models fit to live trap (2000) and hair snag (2020) data
# written by Joanna Burgar (Joanna.Burgar@gov.bc.ca) - 13-Oct-2021
#####################################################################################

.libPaths("C:/Program Files/R/R-4.0.5/library") # to ensure reading/writing libraries from C drive

# Load Packages
list.of.packages <- c("tidyverse","parallel","unmarked", "nimble","scrbook","nimbleSCR","MCMCvis","coda","Cairo","basicMCMCplots","tictoc")

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
# load("out/MartenData_1996.Rda")
# nm.area <- diff(marten.data$xlim)*diff(marten.data$ylim)/100	# Density reported per 100 sq km
# nm.area # 51.21623 100 km2 or 5122 km2 total area (is this true?)

out.files <- list.files("./out/", pattern="*.Rda")
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
  
  # add week and week occasion to the daylookup
  daylookup <- marten.data$daylookup
  glimpse(daylookup)
  daylookup$Week <- week(daylookup$Date)
  num.weeks <- round(nrow(daylookup)/7+1)
  week.occ <- rep(1:num.weeks, each = 7)
  daylookup$Occ_week <- week.occ[1:nrow(daylookup)]
  
  # create a covariate and dataset for weekly occasions
  # covariate is number of days trap open per week (0-7)
  
  # make sure to only use full weeks in occasion week (i.e., remove last week if not containing 7 days)
  weeks.to.use <- daylookup %>% count(Occ_week)
  weeks.to.use <- weeks.to.use[weeks.to.use$n==7,]$Occ_week
  
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
  
  tot.effort <- rowSums(week.effort) # total effort per trap
  
  
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
  y_allT <- as.data.frame(t(y_all %>% dplyr::select(-c(YDay, Date, Occ, Week, Occ_week))))
  dim(y_allT)
  
  y_allT$TrapNum <- colnames(obs.wide[2:ncol(obs.wide)])
  
  tmp1 <- as.data.frame(array(NA, dim=c(nrow(trap.oper), 0)))
  tmp1$TrapNum <- rownames(tmp1)
  tmp2 <- left_join(tmp1, y_allT)
  y_day <- as.matrix(tmp2[,2:ncol(tmp2)])
  colnames(y_day) <- seq_len(ncol(trap.oper))
  rownames(y_day) <- seq_len(nrow(trap.oper))
  y_day[is.na(y_day)] <- 0
  
  # create y for all sites and weekly occasions 
  tmp3 <- as.data.frame(t(y_day))
  tmp3$Occ_week <- daylookup$Occ_week[match(rownames(tmp3), daylookup$Occ)]
  y_week <- tmp3 %>% group_by(Occ_week) %>% summarise_at(1:nrow(trap.oper), sum)
  y_week <- as.matrix(t(y_week %>% dplyr::select(-Occ_week)))
  y_week <- y_week[,1:length(weeks.to.use)]
  
  sum(y_week) # 279 observations
  sum(y_day)  # 280 observations 
  # omitted 3 days from the week occasion because last week not a full week
  
  sum(y_week) == sum(y_day) # same number of observations in weekly and full dataset
  # not the same because of omitted days
  
  retro.data.out[[r]] <- list(y_week, y_day, week.effort, weeks.to.use, daylookup)
  
}


# Look at detections per week and per trap
rowSums(retro.data.out[[1]][[1]]);colSums(retro.data.out[[1]][[1]])
rowSums(retro.data.out[[2]][[1]]);colSums(retro.data.out[[2]][[1]])
rowSums(retro.data.out[[3]][[1]]);colSums(retro.data.out[[3]][[1]])
rowSums(retro.data.out[[4]][[1]]);colSums(retro.data.out[[4]][[1]])

Cairo(file="out/weekly_det_1996.PNG",type="png",width=2000,height=2000,pointsize=15,bg="white",dpi=300)
plot(colSums(retro.data.out[[1]][[1]]), ylab="Weeky marten capture count", xlab="Week", main="Live Capture Data - 1996/97")
dev.off()

Cairo(file="out/weekly_det_1997.PNG",type="png",width=2000,height=2000,pointsize=15,bg="white",dpi=300)
plot(colSums(retro.data.out[[2]][[1]]), ylab="Weeky marten capture count", xlab="Week", main="Live Capture Data - 1997/98")
dev.off()

Cairo(file="out/weekly_det_1998.PNG",type="png",width=2000,height=2000,pointsize=15,bg="white",dpi=300)
plot(colSums(retro.data.out[[3]][[1]]), ylab="Weeky marten capture count", xlab="Week", main="Live Capture Data - 1998/99")
dev.off()

Cairo(file="out/weekly_det_1999.PNG",type="png",width=2000,height=2000,pointsize=15,bg="white",dpi=300)
plot(colSums(retro.data.out[[4]][[1]]), ylab="Weeky marten capture count", xlab="Week", main="Live Capture Data - 1999/00")
dev.off()
# 1999/00 detections high for first 3 weeks and then drop off - unlikley for models to converge

#####################################################################################

#---(begin AHMnimble header)---
# This file contains code adapted from the file R_BUGS_code_AHM_Vol_1_20170519.R
# available at https://www.mbr-pwrc.usgs.gov/pubanalysis/keryroylebook/, with
# permission of the authors.  It may have been modified so that examples originally
# written to work with JAGS or WinBUGS will work with NIMBLE.  More information
# about NIMBLE can be found at https://r-nimble.org, https://github.com/nimble-dev/nimble,
# and https://cran.r-project.org/web/packages/nimble/index.html.
#
# The file  R_BUGS_code_AHM_Vol_1_20170519.R contains the following header:
# -----(begin R_BUGS_code_AHM_Vol_1_20170519.R header)-----
# =========================================================================
#
#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#
#   Marc Kéry & J. Andy Royle
#
#   *** This is the text file with all R and BUGS code from the book ***
#
#   Created 2 Dec 2015 based on draft from 21 Oct 2015
#
# =========================================================================
### Last change: 19 May 2017 by Mike Meredith
# Incorporated errata up to 19 May 2017
# Code updated to implement new names and changes to functions in AHMbook package 0.1.4
# Use built-in data sets instead of .csv files
# In Chapter 11, replaced 'Y', 'Ysum' and 'Yaug' with lower case 'y', 'ysum' and 'yaug'
#  to match the code in the printed book.
#
# -----(end R_BUGS_code_AHM_Vol_1_20170519.R header)-----

# This file was created by:
#
# Jacob Levine and Perry de Valpine
#
#---(end AHMnimble header)---

#####################################################################################
# Specify model in BUGS language
# This corresponds to "model2.txt" in original AHM code.
# without effort
Section6p3_code <- nimbleCode( {
  # Priors
  # lambda ~ dgamma(0.001, 0.001) # came with generic code
  lambda ~ dunif(0,100) # for an uninformative prior
  p ~ dunif(0, 1)
  # Likelihood
  for (i in 1:M) {
    N[i] ~ dpois(lambda)      # State model
    for (j in 1:J) {
      C[i,j] ~ dbin(p, N[i]) # Observation model
    }
  }
})


# MCMC settings
ni <- 25000   ;   nt <- 20   ;   nb <- 5000   ;   nc <- 3

# Parameters monitored
params_noeffort <- c("lambda", "p")

# Bundle data without effort

# list(y_week, y_day, week.effort, weeks.to.use, daylookup)
retro_weekly_noeffort <- vector("list", length(retro.data.out))
for(i in 1:length(retro.data.out)){
  y_week <- retro.data.out[[i]][[1]]
  # y_week <- as.matrix(y_week)
  ndata <- list(C = y_week, M = nrow(y_week), J = ncol(y_week))
  str(ndata)
  
  # Specify initial values
  Nst <- apply(y_week, 1, max)       # Avoid data/model/inits conflict
  inits <- function(){list(N = Nst)}
  
  out <- nimbleMCMC(code = Section6p3_code, 
                    constants = ndata, 
                    inits = inits,
                    monitors = params_noeffort,
                    nburnin = nb, 
                    niter = ni,
                    nchains = nc,
                    samplesAsCodaMCMC = TRUE)
  retro_weekly_noeffort[[i]] <- out
  
}

save(retro_weekly_noeffort, file = paste0("./out/retro_weekly_noeffort_mcmcoutput.Rda"))
# load("out/retro_weekly_noeffort_mcmcoutput.Rda")

out96 <- MCMCsummary(retro_weekly_noeffort[[1]])
out97 <- MCMCsummary(retro_weekly_noeffort[[2]])
out98 <- MCMCsummary(retro_weekly_noeffort[[3]])
out99 <- MCMCsummary(retro_weekly_noeffort[[4]])

(retro.data.out[[4]][1])
retro.out <- rbind(out96, out97, out98, out99)
retro.out$Year <- rep(c("1996/97","1997/98","1998/99","1999/20"), each=2, times=1)
retro.out$param <- rep(c("lambda","p"), each=1, times=4)
retro.out # n.eff low and Gelman high fo 1999/00 = did not converge and should remove from output

#- lambda
nmix.retro.plot.lambda <- ggplot(data = retro.out[retro.out$param=="lambda" & retro.out$Year!="1999/20",]) +
  theme_bw() + theme(strip.background = element_rect(fill = "white", colour = "white")) +
  theme(panel.grid = element_blank())+
  geom_point(aes(x = Year, y = mean*100), size=2) +
  geom_linerange(aes(x = Year, y = mean*100, ymin=`2.5%`*100, ymax= `97.5%`*100)) +
  xlab("Year") +
  ylab("Total Estimated Abundance (i.e., lambda * 100)")+
  ggtitle("Williston Basin marten density estimates;\nlive trap data fit to n-mixture models")

Cairo(file="out/nmix.retro.plot.lambda.PNG",
      type="png",
      width=3000,
      height=2200,
      pointsize=15,
      bg="white",
      dpi=300)
nmix.retro.plot.lambda
dev.off()

#- p
nmix.retro.plot.pall <- ggplot(data = retro.out[retro.out$param=="p",]) +
  theme_bw() + theme(strip.background = element_rect(fill = "white", colour = "white")) +
  theme(panel.grid = element_blank())+
  geom_point(aes(x = Year, y = mean), size=2) +
  geom_linerange(aes(x = Year, y = mean, ymin=`2.5%`, ymax= `97.5%`)) +
  xlab("Year") +
  ylab("Probabilty of Detection")+
  ggtitle("Williston Basin marten detection probabilities;\nlive trap data fit to n-mixture models")

Cairo(file="out/nmix.retro.plot.pall.PNG",
      type="png",
      width=3000,
      height=2200,
      pointsize=15,
      bg="white",
      dpi=300)
nmix.retro.plot.pall 
dev.off()

# excluding 1999/2000
nmix.retro.plot.p <- ggplot(data = retro.out[retro.out$param=="p" & retro.out$Year!="1999/20",]) +
  theme_bw() + theme(strip.background = element_rect(fill = "white", colour = "white")) +
  theme(panel.grid = element_blank())+
  geom_point(aes(x = Year, y = mean), size=2) +
  geom_linerange(aes(x = Year, y = mean, ymin=`2.5%`, ymax= `97.5%`)) +
  xlab("Year") +
  ylab("Probabilty of Detection")+
  ggtitle("Williston Basin marten detection probabilities;\nlive trap data fit to n-mixture models")

Cairo(file="out/nmix.retro.plot.p.PNG",
      type="png",
      width=3000,
      height=2200,
      pointsize=15,
      bg="white",
      dpi=300)
nmix.retro.plot.p 
dev.off()


#####################################################################################

# Section6p4_code <- nimbleCode( {
#   # Priors
#   for(k in 1:3) {                # Loop over 3/4 years
#     alpha0[k] ~ dunif(-10, 10) # Detection intercepts
#     alpha1[k] ~ dunif(-10, 10) # Detection slopes
#     beta0[k] ~ dunif(-10, 10)  # Abundance intercepts
#     beta1[k] ~ dunif(-10, 10)  # Abundance slopes
#   }
#   
#   # Likelihood
#   # Ecological model for true abundance
#   for (i in 1:M){
#     N[i] ~ dpois(lambda[i])
#     log(lambda[i]) <- beta0[year[i]] + beta1[year[i]]
#     # Some intermediate derived quantities
#     # critical[i] <- step(2-N[i])# yields 1 whenever N is 2 or less
#     # z[i] <- step(N[i]-0.5)     # Indicator for occupied site
#     # Observation model for replicated counts
#     for (j in 1:J){
#       C[i,j] ~ dbin(p[i,j], N[i])
#       logit(p[i,j]) <- alpha0[j] + alpha1[j] * effort[i,j]
#     }
#   }
#   
#   # Derived quantities; unnecessary when running for inference purpose
#   Nocc <- sum(z[1:M])         # Number of occupied sites among sample of M
#   Ntotal <- sum(N[1:M])       # Total population size at M sites combined
#   Nyear[1] <- sum(N[1:77])  # Total abundance for sites in year 1 (1996-1997) # 78+105
#   Nyear[2] <- sum(N[78:183]) # Total abundance for sites in year 2 (1997-1998) # 105
#   Nyear[3] <- sum(N[184:283])# Total abundance for sites in year 3 (1998-1999) # 101
#   # Nyear[4] <- sum(N[284:362])# Total abundance for sites in year 4 (1999-2000)   # 79+286
#   # for(k in 1:M){         # Predictions of lambda and p ...
#   #   for(level in 1:3){    #    ... for each level of year and effort factors
#   #     lam.pred[k, level] <- exp(beta0[level] + beta1[level])
#   #     logit(p.pred[k, level]) <- alpha0[level] + alpha1[level] * Xeffort[k]
#   #   }
#   # }
#   # N.critical <- sum(critical[1:M]) # Number of populations with critical size
# })

# M = traps
# J = occasions
# K = years


nmix_effort <- nimbleCode( {
  # Priors - separate for each year to be able to monitor output
  for(k in 1:K) {                # Loop over 3 years
    alpha0[k] ~ dunif(-10, 10) # Detection intercepts
    alpha1[k] ~ dunif(-10, 10) # Detection slopes
    beta0[k] ~ dunif(-10, 10)  # Abundance intercepts
    beta1[k] ~ dunif(-10, 10)  # Abundance slopes
    
    # Likelihood
    # Ecological model for true abundance
    for (i in 1:length(M)[1]){
      N[i,k] ~ dpois(lambda[i,k])
      log(lambda[i,k]) <- beta0[year[i,k]] + beta1[year[i,k]]
      
      # Observation model for replicated counts 
      for (j in 1:J){
        C[i,j,k] ~ dbin(p[i,j,k], N[i,k])
        logit(p[i,j,k]) <- alpha0[j,k] + alpha1[j,k] * effort[i,j,k]
      }
    }
  }
  
  Nyear[k] <- sum(z[1:M,k])         
  D[k] <- Nyear[k]/area
}
)

# compareMCMCs for comparing jags/BUGS and nimble - calls MCMCsuite
# needs modelcode to be exactly the same

# To run as multi-year model, treat each year as a factor, and create arrays with slice for each year
y_week_9697 <- as.matrix(retro.data.out[[1]][[1]])
y_week_9798 <- as.matrix(retro.data.out[[2]][[1]])
y_week_9899 <- as.matrix(retro.data.out[[3]][[1]])
y_week_9900 <- as.matrix(retro.data.out[[4]][[1]])
# dimensions are traps by occasions

# M = traps
# J = occasions
# K = years
# M <- c(nrow(y_week_9697),nrow(y_week_9798),nrow(y_week_9899),nrow(y_week_9900)) # number of traps
# J <- c(ncol(y_week_9697),ncol(y_week_9798),ncol(y_week_9899),ncol(y_week_9900)) # number of occasions
M <- c(nrow(y_week_9697),nrow(y_week_9798),nrow(y_week_9899)) # number of traps
J <- c(ncol(y_week_9697),ncol(y_week_9798),ncol(y_week_9899)) # number of occasions
K <- 3 # try first with 3 years to see if it will run 

y_week_all <- array(c(y_week_9697, y_week_9798, y_week_9899), dim=c(max(M),max(J),K))
dim(y_week_all)
sum(y_week_all)

# create year covariate
year <- array(0, dim=c(max(M),K))
dim(year)
year[,1] <- c(rep(1,M[1]),rep(NA,max(M)-M[1]))
year[,2] <- c(rep(2,M[2]),rep(NA,max(M)-M[2]))
year[,3] <- c(rep(3,M[3]),rep(NA,max(M)-M[3]))


# create effort covariate
effort_9697 <- as.matrix(retro.data.out[[1]][[3]])
effort_9798 <- as.matrix(retro.data.out[[2]][[3]])
effort_9899 <- as.matrix(retro.data.out[[3]][[3]])
effort_9900 <- as.matrix(retro.data.out[[4]][[3]])

# effort_all <- array(c(effort_9697, effort_9798, effort_9899), dim=c(max(M),max(J),K))
effort_all <- array(c(effort_9697, effort_9798), dim=c(max(M),max(J),K))
dim(effort_all)
sum(effort_all)


# Bundle data
# test data
y_week_all.test <- y_week_all[1:min(M),1:min(J),1:3]
effort_all.test <- effort_all[1:min(M),1:min(J),1:3]
year.test <- year[1:77,]
ydata_all.test <- list(C=y_week_all.test, effort=effort_all.test, year=year.test)

constants.test <- list(
  J = min(J),
  area = nm.area,
  M = min(M),
  K = K)

# actual data
ydata_all <- list(C=y_week_all, effort=effort_all, year=year)


constants <- list( # not sure how this works as need to have constants as constants??
  J = J,
  area = nm.area,
  M = M,
  K = K)


N.init = apply(y_week_all, c(1,3), sum)
N.init = ifelse(N.init >=1, 1, 0)

inits = list(N = N.init, 
             alpha0 = rnorm(K), 
             alpha1 = rnorm(K), 
             beta0 = rnorm(K), 
             beta1 = rnorm(K))

# Parameters monitored
params <- c("alpha0", "alpha1", "beta0", "beta1", "Nyear", "D") 

# MCMC settings to test
ni <- 250   ;   nt <- 1  ;   nb <- 50   ;   nc <- 1
# MCMC settings
# nc <- 3   ;   ni <- 22000   ;   nb <- 2000   ;   nt <- 10

nmix.effort.input <- list(constants, ydata_all)
save(nmix.effort.input, file = paste0("./out/nmix.effort.input.Rda"))

# Test with no latent N
nmixR <- nimbleModel(code = nmix_effort,
                     data=ydata_all,
                     constants = constants,
                     inits = inits)

nmixR$calculate()
nmixR$initializeInfo()

# compile model to C++#
nmixC <- compileNimble(nmixR, showCompilerOutput = F)
# MCMC sampler configurations
mcmcspec<-configureMCMC(nmixR, monitors=params)
# build the MCMC specifications
nmixMCMC <- buildMCMC(mcmcspec)
# complile the code in S+
CnmixMCMC <- compileNimble(nmixMCMC, project = nmixR, resetFunctions = TRUE)
# run MCMC
tic()
nmix.results1 <- runMCMC(CnmixMCMC, niter = ni, nburnin=nb,thin=nt,nchains=nc, setSeed = 500)
toc()



save("nmix.results1",file=paste0("out/retro_weekly_effort_mcmcoutput.RData"))
str(nmix.results1)

MCMCsummary(nmix.results1[,1:3], round = 4)
# saved traceplots
chainsPlot(nmix.results1,
           var = c("N", "D"))


require(mcmcplots)
colnames(out)
mcmcplot(out[,1:16])
MCMCsummary(out[[1]])

as.mcmc.list(out)
nrow(y_week)
out_all <- MCMCsummary(out[,1:16])
str(out)
77+105+101+(4*22)


# Calculate ESS effective sample size
# adjusted sample size = effective sample size
# if a "long" Markov chain has only generated a short effective sample size, consider a longer run
apply(samples, 2, effectiveSize)


#####################################################################################
##################------ CURRENT DATA ------##################
#####################################################################################

# run for 2 clusters of marten grid cells
# data wrangled and code fixed by Daniel Eacker (through nimble user group listserv)
load("out/MartenGridData_2020.Rda")

# just marten cells
J <- martenGrid.hsdata$J  # 2 clusters of 20 traps each
area <- c(martenGrid.hsdata$area[[1]],martenGrid.hsdata$area[[2]])
xlim <- array(0,c(2,2))
xlim[1,] <- martenGrid.hsdata$xlim[[1]]
xlim[2,] <- martenGrid.hsdata$xlim[[2]]

ylim <- array(0,c(2,2))
ylim[1,] <- martenGrid.hsdata$ylim[[1]]
ylim[2,] <- martenGrid.hsdata$ylim[[2]]

traps <- array(0,c(J,2,2))
traps.C1 <- as.matrix(martenGrid.hsdata$traps[[1]][,c("x","y")])
traps[,,1] <- traps.C1
traps.C2 <- as.matrix(martenGrid.hsdata$traps[[2]][,c("x","y")])
traps[,,2] <- traps.C2


edf.marten <- martenGrid.hsdata$edf
trials <- cbind(rep(4,20),rep(4,20))
sex <- martenGrid.hsdata$sex
G <- 2 # 2 clusters
M <- 200 # augmented population

# traps are in 2 clusters: 1-20, 21-40
plot(traps[,,1])
plot(traps[,,2])

# y for Cluster 1
y.C1 <- array(0, dim = c(M, 20, 4)) # augmented pop = 200, 20 traps and 4 occasions
# Add the captures as 1s with the good old cbind trick.
edf.marten.C1 <- edf.marten[[1]] %>% filter(Grid_Num < 21)
y.C1[cbind(edf.marten.C1$Animal_Num,edf.marten.C1$Grid_Num, edf.marten.C1$Occ)] <- 1
sum(y.C1) == nrow(edf.marten.C1)     # 6 detections = 3 animals, 1 with 3 detections, 1 with 2 and 1 with 1

# get rid of zeros so observed animals come first
y.C1 = y.C1[which(apply(y.C1,1,sum)>0),,]
dim(y.C1) # 3 rows = 1 for each observed animal

# y for Cluster 2
y.C2 <- array(0, dim = c(M, 20, 4)) # augmented pop = 200, 20 traps and 4 occasions
# Add the captures as 1s with the good old cbind trick.
edf.marten.C2 <- edf.marten[[2]] %>% filter(Grid_Num > 20)
y.C2[cbind(edf.marten.C2$Animal_Num,edf.marten.C2$Grid_Num-20, edf.marten.C2$Occ)] <- 1 # changed the Grid_Num to reflect rownum in traps df
sum(y.C2) == nrow(edf.marten.C2)     # 13 detections = 8 animals, 2 with 3 detections, 1 with 2, and 5 with 1
edf.marten.C2 %>% count(Animal_ID)

# get rid of zeros so observed animals come first
y.C2 = y.C2[which(apply(y.C2,1,sum)>0),,]
dim(y.C2)

# Now let's speed it up by summing over all 4 occasions for a binomial dist.
n0 = c(length(which(apply(y.C1,1,sum)>0)),length(which(apply(y.C2,1,sum)>0)))
y_all.C1 <- apply(y.C1, c(1,2), sum)
dim(y_all.C1)
y_all.C2 <- apply(y.C2, c(1,2), sum)
y_all <- array(0,c(M,J,2)) # needs to be an array, not a list
y_all[1:n0[1],,1] <- y_all.C1
y_all[1:n0[2],,2] <- y_all.C2

z_all <- array(1,c(M,J,2)) # needs to be an array, not a list

SCR_bern <- nimbleCode({
  sigma ~ dunif(0,100) # uninformative prior
  psi ~ dbeta(1,1)
  p0 ~ dunif(0,1)
  
  for(g in 1:G){
    for(i in 1:M){
      z[i,g] ~ dbern(psi)
      s[i,1,g]~dunif(0,20) # traps are centered and scaled, start at 0 and end <20 for both xlim and ylim
      s[i,2,g]~dunif(0,20)
      
      for(j in 1:J){
        d2[i,j,g]<- sqrt((s[i,1,g]-traps[j,1,g])^2 + (s[i,2,g]-traps[j,2,g])^2)
        p[i,j,g]<- z[i,g]*p0*exp(-d2[i,j,g]^2/(sigma*sigma^2))
      }
      
      y[i,1:J,g] ~ dbinom_vector(size = trials[1:J], prob = p[i,1:J,g])
      
    }
    
    N[g]<-sum(z[1:M,g])
    D[g]<-N[g]/area[g]
  }
}
)

constants<- list(
  J = J,
  area = area,
  M = M,
  G = G)

# get average capture locations for detected individuals at starting activity center locations
st=array(NA, c(M,2,G))
for(g in 1:G){
  for(i in 1:n0[g]){ # augmented
    if(sum(y_all[i,,g])==1){
      st[i,1:2,g] = traps[y_all[i,,g],,g]
    }else {
      st[i,1:2,g] = apply(traps[y_all[i,,g],,g], 2, mean)
    }
  }
  for(i in (n0[g]+1):M){
    st[i,1:2,g] = runif(2, 0, 20)
  }}  

data <- list(
  y = y_all, traps = traps, trials=rep(4,J))

z.init = apply(y_all, c(1,3), sum)
z.init = ifelse(z.init >=1, 1, 0)

inits = list(z = z.init, p0 = runif(1,0.05, 1), psi = mean(z.init), sigma = runif(1, 2, 5), s = st)

dim(y_all)

params <- c('sigma', 'p0', 'psi', 'N', 'D')

# MCMC settings to test
# ni <- 250   ;   nt <- 1  ;   nb <- 50   ;   nc <- 1
# MCMC settings for actual run
ni <- 100000   ;   nt <- 20   ;   nb <- 5000   ;   nc <- 3

# Test with no latent N
scrR <- nimbleModel(code = SCR_bern,
                    data=data,
                    constants = constants,
                    inits = inits)

scrR$calculate()
scrR$initializeInfo()

# compile model to C++#
scrC <- compileNimble(scrR, showCompilerOutput = F)
# MCMC sampler configurations
mcmcspec<-configureMCMC(scrR, monitors=params)
# build the MCMC specifications
scrMCMC <- buildMCMC(mcmcspec)
# complile the code in S+
CscrMCMC <- compileNimble(scrMCMC, project = scrR, resetFunctions = TRUE)
# run MCMC
tic()
results3 <- runMCMC(CscrMCMC, niter = ni, nburnin=nb,thin=nt,nchains=nc, setSeed = 500)
toc()
# results1 = ni = 25000 # 482.63/60 # 8 min # poor mixing, don't go with these few ni
# results2 = ni = 50000 # 930.06/60 # 15 min # better mixing
# results3 = ni = 100000 # 2145.6/60 # 36 min # similar to results 2 
# go with ni=50000

save("results2",file=paste0("out/SCRbern2020.RData"))
str(results2)
load("out/SCRbern2020.RData")

MCMCsummary(results3, round = 4)
MCMCsummary(results2, round = 4)

# saved traceplots
chainsPlot(results2,
           var = c("N", "D"))
chainsPlot(results2,
           var = c("sigma","p0", "psi"))