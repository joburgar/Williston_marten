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
# 03a_simulations.R
# script to run N-mixxture models fit to live trap (1996 - 2000) simulated data
# script to run SCR models fit to genetic hair snag (2020) simulated data
# adapted by Joanna Burgar (Joanna.Burgar@gov.bc.ca) - 27-Oct-2021 & 2-Nov-2021
#####################################################################################

# if issues running nimble, run the following line
writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
Sys.which("make")
## "C:\\rtools40\\usr\\bin\\make.exe"

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

# 6.3. Simulation and analysis of the simplest possible N-mixture model
# ------------------------------------------------------------------------

R_version <- paste0("R-",version$major,".",version$minor)

.libPaths(paste0("C:/Program Files/R/",R_version,"/library")) # to ensure reading/writing libraries from C drive
tz = Sys.timezone() # specify timezone in BC

# Load Packages
list.of.packages <- c("tidyverse","nimble","nimbleSCR","mcmcplots","MCMCvis","coda","Cairo","doParallel")
# Check you have them and load them
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)
rm(list.of.packages, new.packages) # for housekeeping

#####################################################################################
load("out/retro.data.out.RData")
load("out/rec.data.out.RData")

nmix.sim.data.function <- function(ydata=ydata,effdata=effdata,lambda=lambda, p=p){
  # ydata <- retro.data.out[[1]]$y_21day
  # effdata <- retro.data.out[[1]]$effort.21days
  # lambda <- 3
  # p <- 0.3
  
  print(dim(ydata))
  print(dim(ydata))
  
  M <- nrow(ydata)                     # Number of trap groups
  J <- ncol(ydata)                      # Number of occasions
  C <- sim.effort <- matrix(NA, nrow = M, ncol = J) # to contain the obs. data and effort data
  
  for(j in 1:J){
    sim.effort[,j] <- round(rnorm(n=M, mean=mean(effdata), sd=sd(effdata)))
  }
  print(sum(sim.effort));print(sum(effdata))
  
  # Generate local abundance data (the truth)
  N <- rpois(n = M, lambda = lambda)
  print(sum(N)) #118
  
  # Conduct repeated measurements (generate replicated counts)
  for(j in 1:J){
    C[,j] <- rbinom(n = M, size = N, prob = p)
  }
  
  # # Look at data
  # # The truth ....
  # table(N)                    # True abundance distribution
  # sum(N)                      # True total population size at M sites
  # sum(N>0)                    # True number of occupied sites
  # mean(N)                     # True mean abundance (estimate of lambda)
  # 
  # # ... and the observations
  # table(apply(C, 1, max))     # Observed abundance distribution (max count)
  # sum(apply(C, 1, max))       # Observed total population size at M sites
  # sum(apply(C, 1, max)>0)     # Observed number of occupied sites
  # mean(apply(C, 1, max))      # Observed mean "relative abundance"
  
  print(sum(C));print(sum(ydata))
  
 return(list(M=M,J=J,N=N,C=C,sim.effort=sim.effort)) 
}

#################################################################################
# Determine sample sizes and simulate observed data array C
# Modified for the Williston Basin, based off of 1996, 1997 and 2019 data
# 21 day trapping sessions in 8,9,4 occasions

# Parameter values 
# based on output from run models, have lambda = 3 for 1996, 2.5 for 1997 and 2.0 for 2019
# lambda <- 3 #1996               # Expected abundance
# based on output from run models, have p = 0.02 for 1996, 0.002 for 1997, and 0.0002 for 2019
# p <- 0.02                    # Probability of detection (per individual)
nmix.sim.data.96 <- nmix.sim.data.function(ydata=retro.data.out[[1]]$y_21day,
                                           effdata=retro.data.out[[1]]$effort.21days,
                                           lambda=3, p=0.3)

nmix.sim.data.97 <- nmix.sim.data.function(ydata=retro.data.out[[2]]$y_21day,
                                           effdata=retro.data.out[[2]]$effort.21days,
                                           lambda=2.5, p=0.2)

nmix.sim.data.19 <- nmix.sim.data.function(ydata=rec.data.out$rec_ydata,
                                           effdata=rec.data.out$rec_effort,
                                           lambda=2, p=0.1)

# nmix.sim.data <- list(nmix.sim.data.96,nmix.sim.data.97,nmix.sim.data.19)
# save(nmix.sim.data, file = paste0("out/nmix_sim_969719_data.RData"))
# load("out/nmix_sim_969719_data.RData")
# nmix.sim.data.96 <- nmix.sim.data[[1]]
# nmix.sim.data.97 <- nmix.sim.data[[2]]
# nmix.sim.data.19 <- nmix.sim.data[[3]]

M <- c(nrow(nmix.sim.data.96$C),nrow(nmix.sim.data.97$C),nrow(nmix.sim.data.19$C)) # number of traps
J <- c(ncol(nmix.sim.data.96$C),ncol(nmix.sim.data.97$C),ncol(nmix.sim.data.19$C)) # number of occasions
K <- length(M)

y_all <- array(c(nmix.sim.data.96$C, nmix.sim.data.97$C, nmix.sim.data.19$C),
               dim=c(max(M),max(J),K))
dim(y_all)
sum(y_all)

effort_all <- array(c(nmix.sim.data.96$sim.effort,nmix.sim.data.97$sim.effort,nmix.sim.data.19$sim.effort),
                    dim=c(max(M),max(J),K))
dim(effort_all)
sum(effort_all)

# simple function to standardize variables
std2=function(x){
  (x - mean(x,na.rm=TRUE))/(2*sd(x,na.rm=TRUE))
}
effort_all = std2(effort_all)


#####################################################################################
# Model code and constants
# Bayesian analysis of the model using NIMBLE:

# from Daniel Eacker (google group listserv initially)
nmix_effort <- nimbleCode( {
  # Priors - separate for each year to be able to monitor output
  alpha0 ~ dnorm(0, 0.5) # intercept on detection probability (a less informative prior perhaps), try plotting: hist(plogis(rnorm(10000,0,sd=sqrt(2))))
  beta0 ~ dnorm(0, 0.1) # intercept on abundance
  alpha2 ~ dnorm(0, 0.001)
  alpha1[1] <- 0 # corner-point constraint for year effects on detection probability
  beta1[1] <- 0 # corner-point constraint for year effects on latent abundance
  
  for(k in 2:K) {                # Loop over 2 years (this is equivalent to coding separate indicator variables)
    beta1[k] ~ dunif(-8,8) # year effects on latent abundance
    alpha1[k] ~ dunif(-8,8) # year effects on detection probability
  }
  
  # Likelihood
  # Ecological model for true abundance
  for(k in 1:K) {
    for (i in 1:M[k]){ # length(M)[1] - you needed to use a nested indexing here since there are a different number of sites in different years
      N[i,k] ~ dpois(lambda[i,k])
      log(lambda[i,k]) <- beta0 + beta1[k]
      # Observation model for replicated counts
      for (j in 1:J[k]){
        y[i,j,k] ~ dbin(p[i,j,k], N[i,k]) # I tend to call the data "y", but this isn't a big deal
        logit(p[i,j,k]) <- alpha0 + alpha1[k] * alpha2*effort[i,j,k] # I use a corner-point constraint to model the effects of year
      }
    }
  }
})

constants <- list(J = J, 
                  M = M,
                  K = K)

data <- list(y = y_all, 
             effort = effort_all)

N.init = apply(y_all, c(1,3), function(x) max(x, na.rm=TRUE)+1)
N.init[is.na(N.init)] <- 1 # to deal with NA's

inits = list(N = N.init,
             alpha0 = rnorm(1,0,0.5),
             alpha1 = c(NA,rnorm(2,0,0.5)),
             alpha2 = rnorm(1,0,0.5),
             beta0 = rnorm(1,0,0.5),
             beta1 = c(NA,rnorm(2,0,0.5)))

# Parameters monitored
params <- c("alpha0", "alpha1", "beta0", "beta1")

# MCMC settings
nc <- 3   ;   ni <- 50000   ;   nb <- 5000

nmix.sim <- list()
for(i in 1:30){
nmix.sim[[i]] <- nimbleMCMC(code = nmix_effort, 
                  data=data,
                  constants = constants, 
                  inits = inits,
                  monitors = params,
                  nburnin = nb, 
                  niter = ni,
                  nchains = nc,
                  samplesAsCodaMCMC = TRUE)
}

save(nmix.sim, file = paste0("out/nmix_sim_969719_mcmcoutput.RData"))
# load("out/nmix_sim_969719_mcmcoutput.RData")
# str(nmix.sim)

nmix.sim.mcmcout <- list()
for(i in 1:length(nmix.sim)){
  nmix.sim.mcmcout[[i]] <- MCMCsummary(nmix.sim[[i]],round = 4)
}


nmix.sim2.df <- as.data.frame(unlist(nmix.sim2.mcmcout))
colnames(nmix.sim2.df) <- "mcmcoutput"
nmix.sim.df <- as.data.frame(unlist(nmix.sim.mcmcout))
colnames(nmix.sim.df) <- "mcmcoutput"
nmix.sim.df <- rbind(nmix.sim2.df, nmix.sim.df)

nrow(nmix.sim.df)
head(nmix.sim.df)
nmix.sim.mcmcout[[1]]
nmix.sim.df[1:8,]
colnames(nmix.sim.df)[1] <- c("value")
rownames(nmix.sim.mcmcout[[1]])
nmix.sim.df$estimate <- rep(c("mean","sd","CI_2.5","CI_50","CI_97.5","Rhat","n.eff"), each=8, times=30)
nmix.sim.df$Sim <- rep(paste0("Sim",seq_len(30)),each=8*7)
nmix.sim.df$param <- rep(rep(rownames(nmix.sim.mcmcout[[1]]),each=1,times=7), times=30)

nmix.sim.df.alpha <- nmix.sim.df %>% filter(grepl("alpha",param))%>% filter(estimate %in% c("CI_2.5","CI_50","CI_97.5"))
nmix.sim.df.alpha <- nmix.sim.df.alpha %>% group_by(Sim,estimate) %>% mutate(new = ifelse(param != 'alpha0', value + value[param == 'alpha0'], value) )

nmix.sim.df.beta <- nmix.sim.df %>% filter(grepl("beta",param)) %>% filter(estimate %in% c("CI_2.5","CI_50","CI_97.5"))
nmix.sim.df.beta <- nmix.sim.df.beta %>% group_by(Sim,estimate) %>% mutate(new = ifelse(param != 'beta0', value + value[param == 'beta0'], value) )

nmix.sim.wide.alpha <- nmix.sim.df.alpha %>% dplyr::select(-value) %>% filter(param!="alpha0") %>% pivot_wider(names_from = estimate, values_from = new)
nmix.sim.wide.beta <- nmix.sim.df.beta %>% dplyr::select(-value) %>% filter(param!="beta0") %>% pivot_wider(names_from = estimate, values_from = new)


#- probability of detection - NEED TO BACK TRANSFORM CONSIDERING EFFORT VARIABLE WAS SCALED
# logit(p[i,j,k]) <- alpha0 + alpha1[k] * alpha2*effort[i,j,k]
# not sure I have it totally correct, but regardless it shows how imprecise probability of detection gets with fewer occassions
# nmix.sim.wide.alpha <- nmix.sim.wide.alpha %>% ungroup()


nmix.sim.wide.alpha <- nmix.sim.wide.alpha %>% 
  mutate(CI_2.5_sd = case_when(param=="alpha1[1]" ~ boot::inv.logit(CI_2.5/(2*sd(nmix.sim.data.96$sim.effort))),
                               param=="alpha1[2]" ~ boot::inv.logit(CI_2.5/(2*sd(nmix.sim.data.97$sim.effort))),
                               param=="alpha1[3]" ~ boot::inv.logit(CI_2.5/(2*sd(nmix.sim.data.19$sim.effort)))))


nmix.sim.wide.alpha <- nmix.sim.wide.alpha %>% mutate(CI_50_sd = case_when(param=="alpha1[1]" ~ boot::inv.logit(CI_50/(2*sd(nmix.sim.data.96$sim.effort))),
                                                                            param=="alpha1[2]" ~ boot::inv.logit(CI_50/(2*sd(nmix.sim.data.97$sim.effort))),
                                                                            param=="alpha1[3]" ~ boot::inv.logit(CI_50/(2*sd(nmix.sim.data.19$sim.effort)))))


nmix.sim.wide.alpha <- nmix.sim.wide.alpha %>% mutate(CI_97.5_sd = case_when(param=="alpha1[1]" ~ boot::inv.logit(CI_97.5/(2*sd(nmix.sim.data.96$sim.effort))),
                                                                            param=="alpha1[2]" ~ boot::inv.logit(CI_97.5/(2*sd(nmix.sim.data.97$sim.effort))),
                                                                            param=="alpha1[3]" ~ boot::inv.logit(CI_97.5/(2*sd(nmix.sim.data.19$sim.effort)))))

nmix.sim.wide.alpha$param <- as.character(rep(c("1996-1997","1997-1998","2019-2020"),each=1, times=30))

# for horizontal line by facet wrap (param)
Pgmean <- nmix.sim.wide.alpha %>%
  group_by(param) %>%
  summarise(Pg = mean(CI_50))
Pgmean$Pg <- c(0.3, 0.2, 0.1)


nmix.sim.plot.dp <- ggplot(data = nmix.sim.wide.alpha) +
  theme_bw() + theme(strip.background = element_rect(fill = "white", colour = "white")) +
  theme(panel.grid = element_blank())+
  geom_point(aes(x = Sim, y = CI_50_sd), size=2) +
  geom_hline(data= Pgmean, aes(yintercept = Pg), col="grey") +
  geom_linerange(aes(x = Sim, y = CI_50_sd, ymin=CI_2.5_sd, ymax=CI_97.5_sd)) +
  theme(axis.text.x = element_blank()) +
  xlab("Simulation Runs") +
  ylab("Estimated Probability of Detection")+
  # ggtitle("Simulations of n-mixture models with 100 sampling occasions,\nlambda=c(2,1,0.5), and p=c(0.5,0.4,0.3)")+
  facet_wrap(~ param)

nmix.sim.plot.dp

Cairo(file="out/nmix.sim.plot.dp.PNG",type="png",width=3000,height=2200,pointsize=15,bg="white",dpi=300)
nmix.sim.plot.dp
dev.off()

#- latent abundance
Ng96 <- sum(nmix.sim.data.96$N)/nmix.sim.data.96$M # put in as 3 (3.0)
Ng97 <- sum(nmix.sim.data.97$N)/nmix.sim.data.97$M # put in as 2.5 (2.5)
Ng19 <- sum(nmix.sim.data.19$N)/nmix.sim.data.19$M # put in as 2 (2.0)

# to rename the facets by survey year
nmix.sim.wide.beta$param <- as.character(rep(c("1996-1997","1997-1998","2019-2020"),each=1, times=30))

# for horizontal line by facet wrap (param)
Ngmean <- nmix.sim.wide.beta %>%
  group_by(param) %>%
  summarise(Ng = mean(CI_50))
Ngmean$Ngactual <- c(3.0, 2.6, 1.9)

Ngmean$Nactual <- c(sum(nmix.sim.data.96$N),sum(nmix.sim.data.97$N),sum(nmix.sim.data.19$N))

nmix.sim.wide.beta$J <- (rep(c(nmix.sim.data.96$M,nmix.sim.data.97$M,nmix.sim.data.19$M),each=1, times=30))
nmix.sim.wide.beta$CI_2.5N <- nmix.sim.wide.beta$CI_2.5*nmix.sim.wide.beta$J
nmix.sim.wide.beta$CI_50N <- nmix.sim.wide.beta$CI_50*nmix.sim.wide.beta$J
nmix.sim.wide.beta$CI_97.5N <- nmix.sim.wide.beta$CI_97.5*nmix.sim.wide.beta$J

nmix.sim.plot.la <- ggplot(data = nmix.sim.wide.beta) +
  theme_bw() + theme(strip.background = element_rect(fill = "white", colour = "white")) +
  theme(panel.grid = element_blank())+
  geom_hline(data=Ngmean, aes(yintercept=Ngactual), col="grey") +
  geom_point(aes(x = Sim, y = CI_50), size=2) +
  geom_linerange(aes(x = Sim, y = CI_50, ymin=CI_2.5, ymax=CI_97.5)) +
  theme(axis.text.x = element_blank()) +
  xlab("Simulation Runs") +
  ylab("Estimated Median Abundance per Grid Cell")+
  ggtitle("N-mixture Model Simulations")+
  facet_wrap(~ param)

Cairo(file="out/nmix.sim.plot.la.PNG",type="png",width=3000,height=2200,pointsize=15,bg="white",dpi=300)
nmix.sim.plot.la
dev.off()

nmix.sim.plot.N <- ggplot(data = nmix.sim.wide.beta) +
  theme_bw() + theme(strip.background = element_rect(fill = "white", colour = "white")) +
  theme(panel.grid = element_blank())+
  geom_hline(data=Ngmean, aes(yintercept=Nactual), col="grey") +
  geom_point(aes(x = Sim, y = CI_50N), size=2) +
  geom_linerange(aes(x = Sim, y = CI_50N, ymin=CI_2.5N, ymax=CI_97.5N)) +
  theme(axis.text.x = element_blank()) +
  xlab("Simulation Runs") +
  ylab("Estimated Median Abundance per Study Area")+
  ggtitle("N-mixture Model Simulations")+
  facet_wrap(~ param)

Cairo(file="out/nmix.sim.plot.N.PNG",type="png",width=3000,height=2200,pointsize=15,bg="white",dpi=300)
nmix.sim.plot.N
dev.off()

#####################################################################################
##################------ CURRENT DATA ------##################
#####################################################################################
# e2dist from scrbook
e2dist <- function (x, y) {
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

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


###--- simulate the data - this is for one cluster only

#####################################################################################
# Model code and constants
# Bayesian analysis of the model using NIMBLE:
# no clusters
SCR_bern_nocluster <- nimbleCode({
  sigma ~ dunif(0,100) # uninformative prior
  psi ~ dbeta(1,1)
  p0 ~ dunif(0,1)
  
  for(i in 1:M){
    z[i] ~ dbern(psi)
    s[i,1]~dunif(0,20) # traps are centered and scaled, start at 0 and end <20 for both xlim and ylim
    s[i,2]~dunif(0,20)
    
    for(j in 1:J){
      d2[i,j]<- sqrt((s[i,1]-traps[j,1])^2 + (s[i,2]-traps[j,2])^2)
      p[i,j]<- z[i]*p0*exp(-d2[i,j]^2/(sigma*sigma^2))
    }
    
    y[i,1:J] ~ dbinom_vector(size = trials[1:J], prob = p[i,1:J])
    
  }
  
  N<-sum(z[1:M])
  D<-N/area
}
)

# keep in mind that WGS84 lat/long espg = 4326; BC Albers espg = 3005; NAD83 / UTM zone 10N espg = 26910 
### Simulate for Kara's marten traps
# first upload and wrangle Kara's trap data
library(sf)
om_marten <- read.csv("data/marten_deployment_pts_v3.csv", row.names=1)
om_marten_sf <- st_as_sf(om_marten, coords=c("Long","Lat"), crs = 4326)
om_marten_utm <- st_transform(om_marten_sf, crs=26910)
st_coordinates(om_marten_utm)

ggplot()+
  geom_sf(data = om_marten_utm)

# For the marten only trap data
coord.scale <- 1000
buffer <- 5 #5 km unit buffer

traps.scale <- as.data.frame(st_coordinates(om_marten_utm)/coord.scale)

traps.sc <- as.data.frame(cbind(traps.scale$X-min(traps.scale$X-buffer), traps.scale$Y-min(traps.scale$Y-buffer)))
colnames(traps.sc) <- c("x","y")
rownames(traps.sc) <- rownames(traps.scale)
plot(traps.sc)
xlim = range(traps.sc[,1])+c(-buffer,buffer)
ylim = range(traps.sc[,2])+c(-buffer,buffer)
area <- diff(xlim)*diff(ylim)/100	# Density reported per 100 sq km
area # 3.72 or 372 km2

round(area*15) # 56 population based on expected density at 15 marten per 100 sq km
round(area*25) # 93 population based on expected density at 25 per 100 sq km
round(area*35) # 130 population based on expected density at 35 per 100 sq km

###--- for all models
K <- 4 # number of occasions
N = 56
M = 200
J = 42
p0 <- 0.5   # define parameters of encounter probability
sigma <- 1.5 # scale parameter of half-normal

params <- c('sigma', 'p0', 'psi', 'N', 'D')
ni <- 50000  ;   nb <- 5000   ;   nc <- 3

###--- now run 10 simulations per parameter above
# alter N as per 15, 25 and 35 marten per 100 sq km density
# will also need to change output save at end of loop

for(t in 10:10){
  # simulate activity centres 
  sx <- runif(N, xlim[1],xlim[2])
  sy <- runif(N, ylim[1],ylim[2])
  S <- cbind(sx, sy)
  
  # compute distance matrix
  D <- e2dist(S, traps.sc) # distance of each individual from each trap
  
  # Parameter values
  p0 <- p0   # define parameters of encounter probability
  sigma <- sigma # scale parameter of half-normal
  alpha1 <- 1/(2*sigma*sigma) # convert to coefficient on distance
  
  # Compute probability of encounter
  probcap1 <- plogis(-2.5)*exp(-alpha1*D*D)
  
  # Generate the encounters of every individual in every trap
  ntraps <- nrow(traps.sc)
  
  Y1 <- matrix(NA, nrow=N, ncol=ntraps)
  for(i in 1:nrow(Y1)){
    Y1[i,] <- rbinom(ntraps,K,probcap1[i,])
    }
  Y1sum = Y1[which(apply(Y1,1,sum)>0),]
  
  y_sim <- array(0,c(M,J)) # needs to be an array, not a list
  
  n0 = c(length(which(apply(Y1,1,sum)>0)))
  y_sim[1:n0,] <- Y1sum
  
  # get average capture locations for detected individuals at starting activity center locations
  traps = traps.sc
  st=array(NA, c(M,2))
  for(i in 1:M){
    st[i,1:2] = runif(2, 0, 20)
    }  
  
  constants<- list(J = J, area = area, M = M)
  
  data <- list(y = y_sim, traps = traps, trials=rep(4,J))
  
  z.init = apply(y_sim, c(1,2), sum)
  z.init = ifelse(rowSums(z.init)>=1, 1, 0)
  
  inits = list(z = z.init, p0 = runif(1,0.05, 1), psi = mean(z.init), sigma = runif(1, 2, 5), s = st)
  
  SCR.sim.out <- nimbleMCMC(code = SCR_bern_nocluster, 
                          data = data,
                          constants = constants,
                          inits = inits,
                          monitors = params,
                          niter = ni, 
                          nburnin = nb,
                          nchains = nc,
                          samplesAsCodaMCMC = TRUE)
  
  SCR.sim <- MCMCsummary(SCR.sim.out, round = 4)
  save("SCR.sim", file=paste0("out/scrsim_Kara/SCR_SimN56","_",t,".RData"))
  }


###--- for all models
K <- 4 # number of occasions
N = 93
M = 200
J = 42
p0 <- 0.5   # define parameters of encounter probability
sigma <- 1.5 # scale parameter of half-normal

params <- c('sigma', 'p0', 'psi', 'N', 'D')
ni <- 50000  ;   nb <- 5000   ;   nc <- 3

###--- now run 10 simulations per parameter above
# alter N as per 15, 25 and 35 marten per 100 sq km density
# will also need to change output save at end of loop

for(t in 1:10){
  # simulate activity centres 
  sx <- runif(N, xlim[1],xlim[2])
  sy <- runif(N, ylim[1],ylim[2])
  S <- cbind(sx, sy)
  
  # compute distance matrix
  D <- e2dist(S, traps.sc) # distance of each individual from each trap
  
  # Parameter values
  p0 <- p0   # define parameters of encounter probability
  sigma <- sigma # scale parameter of half-normal
  alpha1 <- 1/(2*sigma*sigma) # convert to coefficient on distance
  
  # Compute probability of encounter
  probcap1 <- plogis(-2.5)*exp(-alpha1*D*D)
  
  # Generate the encounters of every individual in every trap
  ntraps <- nrow(traps.sc)
  
  Y1 <- matrix(NA, nrow=N, ncol=ntraps)
  for(i in 1:nrow(Y1)){
    Y1[i,] <- rbinom(ntraps,K,probcap1[i,])
  }
  Y1sum = Y1[which(apply(Y1,1,sum)>0),]
  
  y_sim <- array(0,c(M,J)) # needs to be an array, not a list
  
  n0 = c(length(which(apply(Y1,1,sum)>0)))
  y_sim[1:n0,] <- Y1sum
  
  # get average capture locations for detected individuals at starting activity center locations
  traps = traps.sc
  st=array(NA, c(M,2))
  for(i in 1:M){
    st[i,1:2] = runif(2, 0, 20)
  }  
  
  constants<- list(J = J, area = area, M = M)
  
  data <- list(y = y_sim, traps = traps, trials=rep(4,J))
  
  z.init = apply(y_sim, c(1,2), sum)
  z.init = ifelse(rowSums(z.init)>=1, 1, 0)
  
  inits = list(z = z.init, p0 = runif(1,0.05, 1), psi = mean(z.init), sigma = runif(1, 2, 5), s = st)
  
  SCR.sim.out <- nimbleMCMC(code = SCR_bern_nocluster, 
                            data = data,
                            constants = constants,
                            inits = inits,
                            monitors = params,
                            niter = ni, 
                            nburnin = nb,
                            nchains = nc,
                            samplesAsCodaMCMC = TRUE)
  
  SCR.sim <- MCMCsummary(SCR.sim.out, round = 4)
  save("SCR.sim", file=paste0("out/scrsim_Kara/SCR_SimN93","_",t,".RData"))
}

###---
# Load simulated runs
# load("out/scrsim/SCR.sim.out_Sim01_1.RData")
list.sims <- list.files("out/scrsim_Kara/")


scrsim01.out <- vector('list', length(list.sims))
for(i in 1:length(list.sims)){
  load(paste0("out/scrsim_Kara/",list.sims[i]))
  scrsim01.out[[i]] <- SCR.sim
}

#3 sims: N130 = d35, N56 = d15, N93 = d25

scrsim01.df <- as.data.frame(unlist(scrsim01.out))
nrow(scrsim01.df)
glimpse(scrsim01.df)
scrsim01.df[1:14,]
colnames(scrsim01.df)[1] <- c("value")
scrsim01.df$estimate <- rep(c("mean","sd","CI_2.5","CI_50","CI_97.5","Rhat","n.eff"), each=5, time=length(list.sims))
scrsim01.df$param <- rep(c("D","N","p0","psi","sigma"),each=1, time=7*length(list.sims))
scrsim01.df$Sim <- rep(c("d35","d15","d25"),each=5*7*10)
scrsim01.df$Run <- rep(paste0("Run", seq_len(10)),each=5*7)
glimpse(scrsim01.df)
scrsim01.wide <- pivot_wider(scrsim01.df, names_from = estimate, values_from = value)

#- Density

hline_dat = data.frame(Sim=c("d15", "d25","d35"),threshold=c(15,25,35))

scr.sim.plot.density <- ggplot(data = scrsim01.wide[grepl("D",scrsim01.wide$param),]) +
  theme_bw() + theme(strip.background = element_rect(fill = "white", colour = "white")) +
  theme(panel.grid = element_blank())+
  geom_point(aes(x = Run, y = mean), size=2) +
  geom_hline(data=hline_dat, aes(yintercept=threshold), colour="grey")+
  geom_linerange(aes(x = Run, y = mean, ymin=CI_2.5, ymax= CI_97.5)) +
  theme(axis.text.x = element_blank()) +
  xlab("Simulation Runs") +
  ylab("Mean Density (marten per 100 sq km")+
  ggtitle("Simulations of SCR models in the Omineca:\n42 traps, 4 occassions, p=0.5, sigma=1.5")+
  facet_wrap(~ Sim)

Cairo(file="out/scrsim_Kara_plot.density.PNG",
      type="png",
      width=3000,
      height=2200,
      pointsize=15,
      bg="white",
      dpi=300)
scr.sim.plot.density
dev.off()


#####################################################################################
# Start sets of simulations
# with clusters
# no clusters
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

# MCMC settings for actual run
ni <- 50000   ;   nt <- 20   ;   nb <- 5000   ;   nc <- 3

# create function to run through simulations
scr.function <- function(simname=simname, xlims=xlim, ylims=ylim, traps1=traps.C1, traps2=traps.C2,
                         N=50, J=20, G=2, p0=0.5, sigma=1.5, K=4, M=200){
  
  # SCR.sim <- vector('list', numsim)
  # names(SCR.sim) <- paste0('SCR.sim', seq_along(SCR.sim))
  # for(s in 1:numsim){
    
    # simulate activity centres
    sx1 <- runif(N, xlims[1,1],xlims[1,2])
    sy1 <- runif(N, ylims[1,1],ylims[1,2])
    S1 <- cbind(sx1, sy1)
    
    sx2 <- runif(N, xlims[2,1],xlims[2,2])
    sy2 <- runif(N, ylims[2,1],ylims[2,2])
    S2 <- cbind(sx2, sy2)
    
    # compute distance matrix
    D1 <- e2dist(S1, traps1) # distance of each individual from each trap
    D2 <- e2dist(S2, traps2) # distance of each individual from each trap
    
    # Parameter values
    p0 <- p0   # define parameters of encounter probability
    sigma <- sigma # scale parameter of half-normal
    alpha1 <- 1/(2*sigma*sigma) # convert to coefficient on distance
    
    # Compute probability of encounter
    probcap1 <- plogis(-2.5)*exp(-alpha1*D1*D1)
    probcap2 <- plogis(-2.5)*exp(-alpha1*D2*D2)
    
    # Generate the encounters of every individual in every trap
    ntraps <- J
    
    Y1 <- matrix(NA, nrow=N, ncol=ntraps)
    for(y in 1:nrow(Y1)){
      Y1[y,] <- rbinom(ntraps,K,probcap1[y,])
    }
    
    Y2 <- matrix(NA, nrow=N, ncol=ntraps)
    for(y in 1:nrow(Y2)){
      Y2[y,] <- rbinom(ntraps,K,probcap2[y,])
    }
    
    Y0 <-  matrix(0, 1, ncol=ntraps)
    
    Y1sum = Y1[which(apply(Y1,1,sum)>0),]
    # dim(Y1sum) # 1 row for each observed animal
    Y2sum = Y2[which(apply(Y2,1,sum)>0),]
    # dim(Y2sum) # 1 row for each observed animal

    Y1use <- if (sum(Y1) > 0) Y1sum else Y0
    Y2use <- if (sum(Y2) > 0) Y2sum else Y0
    
    y_sim <- array(0,c(M,J,2)) # needs to be an array, not a list
    
    n0 = c(length(which(apply(Y1use,1,sum)>0)),length(which(apply(Y2use,1,sum)>0)))
    
    y_sim[1:n0[1],,1] <- Y1use
    y_sim[1:n0[2],,2] <- Y2use
    # dim(y_sim)
    # sum(y_sim)
    # class(traps)
    constants<- list(
      J = J,
      area = area,
      M = M,
      G = G)
    
    # get average capture locations for detected individuals at starting activity center locations
    st=array(NA, c(M,2,G))
    for(g in 1:G){
      for(i in 1:n0[g]){ # augmented
        if(sum(y_sim[i,,g])==1){
          st[i,1:2,g] = traps[y_sim[i,,g],,g]
        }else {
          st[i,1:2,g] = apply(as.matrix(traps[y_sim[i,,g],,g]), 2, mean)
        }
      }
      for(i in (n0[g]+1):M){
        st[i,1:2,g] = runif(2, 0, 20)
      }}  
    
    
    data <- list(y = y_sim, traps = traps, trials=rep(4,J))
    
    z.init = apply(y_sim, c(1,3), sum)
    z.init = ifelse(z.init >=1, 1, 0)
    
    inits = list(z = z.init, p0 = runif(1,0.05, 1), psi = mean(z.init), sigma = runif(1, 2, 5), s = st)
    
    params <- c('sigma', 'p0', 'psi', 'N', 'D')
    
    SCR.sim.out <- nimbleMCMC(code = SCR_bern, 
                              data = data,
                              constants = constants,
                              inits = inits,
                              monitors = params,
                              niter = ni, 
                              nburnin = nb,
                              nchains = nc,
                              # thin=nt,
                              samplesAsCodaMCMC = TRUE)
    SCR.sim <- MCMCsummary(SCR.sim.out, round = 4)
    return(SCR.sim)
    
    }

# Sim01 = N = 50, p0 = 0.5, sigma = 2
# Sim02 = N = 50, p0 = 0.4, sigma = 2
# Sim03 = N = 50, p0 = 0.3, sigma = 2
# Sim04 = N = 50, p0 = 0.5, sigma = 1.5
# Sim05 = N = 50, p0 = 0.4, sigma = 1.5
# Sim06 = N = 50, p0 = 0.3, sigma = 1.5
# Sim07 = N = 50, p0 = 0.5, sigma = 1
# Sim08 = N = 50, p0 = 0.4, sigma = 1
# Sim09 = N = 50, p0 = 0.3, sigma = 1
# Sim10 = N = 100, p0 = 0.5, sigma = 2


# running simulations with 2 clusters, each with 20 traps, open for 4 occasions
# real population at each trap = 44, varying p0 and sigma
for(i in 1:10){
  sim.out <- scr.function(simname=c("Sim05"), N=50, J=20, G=2, K=4, p0=0.4, sigma=1.5)
  save("sim.out", file=paste0("out/scrsim/SCR.sim.out_Sim05","_",i,".RData"))
}

for(i in 1:10){
  sim.out <- scr.function(simname=c("Sim06"), N=50, J=20, G=2, K=4, p0=0.3, sigma=1.5)
  save("sim.out", file=paste0("out/scrsim/SCR.sim.out_Sim06","_",i,".RData"))
}

for(i in 1:10){
  sim.out <- scr.function(simname=c("Sim07"), N=50, J=20, G=2, K=4, p0=0.5, sigma=1)
  save("sim.out", file=paste0("out/scrsim/SCR.sim.out_Sim07","_",i,".RData"))
}

for(i in 1:10){
  sim.out <- scr.function(simname=c("Sim08"), N=50, J=20, G=2, K=4, p0=0.4, sigma=1)
  save("sim.out", file=paste0("out/scrsim/SCR.sim.out_Sim08","_",i,".RData"))
}

for(i in 1:10){
  sim.out <- scr.function(simname=c("Sim09"), N=50, J=20, G=2, K=4, p0=0.3, sigma=1)
  save("sim.out", file=paste0("out/scrsim/SCR.sim.out_Sim09","_",i,".RData"))
}

scr.Sim01 <- scr.function(simname=c("Sim01"), N=50, J=20, G=2, K=4, p0=0.5, sigma=2, numsim=10)
scr.Sim02 <- scr.function(simname=c("Sim02"), N=50, J=20, G=2, K=4, p0=0.4, sigma=2, numsim=10)
scr.Sim03 <- scr.function(simname=c("Sim03"), N=50, J=20, G=2, K=4, p0=0.3, sigma=2, numsim=10)
scr.Sim04 <- scr.function(simname=c("Sim04"), N=50, J=20, G=2, K=4, p0=0.5, sigma=1.5, numsim=10)
scr.Sim05 <- scr.function(simname=c("Sim05"), N=50, J=20, G=2, K=4, p0=0.4, sigma=1.5, numsim=10)
scr.Sim06 <- scr.function(simname=c("Sim06"), N=50, J=20, G=2, K=4, p0=0.3, sigma=1.5, numsim=10)
scr.Sim07 <- scr.function(simname=c("Sim07"), N=50, J=20, G=2, K=4, p0=0.5, sigma=1, numsim=10)
scr.Sim08 <- scr.function(simname=c("Sim08"), N=50, J=20, G=2, K=4, p0=0.4, sigma=1, numsim=10)
scr.Sim09 <- scr.function(simname=c("Sim09"), N=50, J=20, G=2, K=4, p0=0.3, sigma=1, numsim=10)


###---
# Load simulated runs
# load("out/scrsim/SCR.sim.out_Sim01_1.RData")
list.sims <- list.files("out/scrsim/")

list.sims01 <- list.sims[grepl("Sim01",list.sims)]


scrsim01.out <- vector('list', length(list.sims01))
for(i in 1:length(list.sims01)){
  load(paste0("out/scrsim/",list.sims01[i]))
  scrsim01.out[[i]] <- sim.out
}

scrsim01.df <- as.data.frame(unlist(scrsim01.out))
nrow(scrsim01.df)
head(scrsim01.df)
scrsim01.df[1:14,]
colnames(scrsim01.df)[1] <- c("value")
scrsim01.df$estimate <- rep(c("mean","sd","CI_2.5","CI_50","CI_97.5","Rhat","n.eff"), each=7, time=length(list.sims01))
scrsim01.df$param <- rep(c("D1","D2","N1","N2","p0","psi","sigma"),each=1, time=7*length(list.sims01))
scrsim01.df$Sim <- rep(paste0("Sim0",seq_len(length(list.sims01))),each=7*7)
scrsim01.df$Run <- rep(paste0("Run", seq_len(10)),each=7*7)
scrsim01.wide <- pivot_wider(scrsim01.df, names_from = estimate, values_from = value)

glimpse(scrsim01.wide)

#- Density
50/area

scr.sim.plot.density <- ggplot(data = scrsim01.wide[grepl("D",scrsim01.wide$param),]) +
  theme_bw() + theme(strip.background = element_rect(fill = "white", colour = "white")) +
  theme(panel.grid = element_blank())+
  geom_point(aes(x = Run, y = mean), size=2) +
  geom_hline(yintercept = c(50/area), col="grey") +
  geom_linerange(aes(x = Run, y = mean, ymin=CI_2.5, ymax= CI_97.5)) +
  theme(axis.text.x = element_blank()) +
  xlab("Simulation Runs") +
  ylab("Mean Density (marten per 100 sq km")+
  ggtitle("Simulations of SCR models in the Williston Basin (2020):\n2 clusters, 4 occassions, N=50, p0=0.5, sigma=2")+
  facet_wrap(~ param)

Cairo(file="out/scr.sim.plot.density.PNG",
      type="png",
      width=3000,
      height=2200,
      pointsize=15,
      bg="white",
      dpi=300)
scr.sim.plot.density
dev.off()

#- p
scr.sim.plot.p0sigma <- ggplot(data = scrsim01.wide[grepl("p0|sigma",scrsim01.wide$param),]) +
  theme_bw() + theme(strip.background = element_rect(fill = "white", colour = "white")) +
  theme(panel.grid = element_blank())+
  geom_point(aes(x = Run, y = mean), size=2) +
  geom_linerange(aes(x = Run, y = mean, ymin=CI_2.5, ymax= CI_97.5)) +
  theme(axis.text.x = element_blank()) +
  xlab("Simulation Runs") +
  ylab("Mean of the Estimated Parameter")+
  ggtitle("Simulations of SCR models in the Williston Basin (2020):\n2 clusters, 4 occassions, N=50, p0=0.5, sigma=2")+
  facet_wrap(~ param, scales="free_y")

Cairo(file="out/scr.sim.plot.p0sigma.PNG",
      type="png",
      width=3000,
      height=2200,
      pointsize=15,
      bg="white",
      dpi=300)
scr.sim.plot.p0sigma
dev.off()


nmix.sim.df[1:25,] 
mean.lambda <- 0.9385
M*mean.lambda # estimated N
sum(N) # compared to true N

