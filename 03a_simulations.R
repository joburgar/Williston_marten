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

.libPaths("C:/Program Files/R/R-4.1.1/library") # to ensure reading/writing libraries from C drive

# Load Packages
list.of.packages <- c("tidyverse","nimble","mcmcplots","MCMCvis","coda","Cairo","nimbleSCR","tictoc","basicMCMCplots")
# Check you have them and load them
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)
rm(list.of.packages, new.packages) # for housekeeping

#####################################################################################
load("out/MartenData_1996.Rda")
glimpse(marten.data)
trap.oper <- marten.data$trap.oper
daylookup <- marten.data$daylookup

# add week and week occasion to the daylookup
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

count = 1 
for(i in 1:length(week.effort)){
  # i=8
  tmp1 <- as.data.frame(t(trap.oper))
  tmp1$Occ_week <- daylookup$Occ_week[match(rownames(tmp1), as.character(daylookup$Date))]
  tmp1 <- tmp1 %>% filter(Occ_week %in% weeks.to.use)
  tmp1 %>% count(Occ_week)
  tmp2 <- tmp1 %>% filter(Occ_week==i) %>% colSums()
  
  week.effort[,count] <- tmp2[1:nrow(trap.oper)]
  
  count <- count + 1
}

tot.effort <- rowSums(week.effort) # total effort per trap
summary(tot.effort)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3.00   23.00   39.00   43.31   65.00   98.00
# try sample sessions (J) = 25, 50, 75, 100 to somewhat mimic actual effort

# Choose sample sizes and prepare observed data array C
# Modified for the Williston Basin, based off of 1996 data
# set.seed(24)                # So we all get same data set
M <- nrow(trap.oper)                     # Number of sites (williston Basin in 1996) = 77
J <- ncol(trap.oper)                      # Number of abu. measurements per site (rep. counts) = 178
C <- matrix(NA, nrow = M, ncol = J) # to contain the obs. data

M <- 77                     # Number of sites (williston Basin in 1996) = 77
J <- 178                      # Number of abu. measurements per site (rep. counts) = 178
C <- matrix(NA, nrow = M, ncol = J) # to contain the obs. data


# Parameter values
lambda <- 1               # Expected abundance
p <- 0.4                    # Probability of detection (per individual)

# Generate local abundance data (the truth)
N <- rpois(n = M, lambda = lambda)
# sum(N)

# Conduct repeated measurements (generate replicated counts)
for(j in 1:J){
  C[,j] <- rbinom(n = M, size = N, prob = p)
}

# Look at data
# The truth ....
table(N)                    # True abundance distribution
sum(N)                      # True total population size at M sites
sum(N>0)                    # True number of occupied sites
mean(N)                     # True mean abundance (estimate of lambda)

# ... and the observations
table(apply(C, 1, max))     # Observed abundance distribution (max count)
sum(apply(C, 1, max))       # Observed total population size at M sites
sum(apply(C, 1, max)>0)     # Observed number of occupied sites
mean(apply(C, 1, max))      # Observed mean "relative abundance"

head(cbind(N=N, count1=C[,1], count2=C[,2])) # First 6 sites

cor(C)[1,2]

#####################################################################################
# Model code and constants
# Bayesian analysis of the model using NIMBLE:

# Specify model in BUGS language:
# This code corresponds to "model1.txt" in the AHM code
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


# Parameters monitored
params <- c("lambda", "p")

# MCMC settings
ni <- 25000   ;   nt <- 20   ;   nb <- 5000   ;   nc <- 3

#####################################################################################
# Start sets of simulations

# create function to run through simulations
nmix.function <- function(simname=simname, M=77, J=J, lambda=lambda, p=p, numsim=25){
  
  nmix.sim <- vector('list', numsim)
  names(nmix.sim) <- paste0('nmix.sim', seq_along(nmix.sim))
  for(i in seq_along(nmix.sim)){
    
    C <- matrix(NA, nrow = M, ncol = J) # to contain the obs. data
    
    # Generate local abundance data (the truth)
    N <- rpois(n = M, lambda = lambda)

    # Conduct repeated measurements (generate replicated counts)
    for(j in 1:J){
      C[,j] <- rbinom(n = M, size = N, prob = p)
    }
    
    # Bundle and summarize data set
    win.data <- list(C = C, M = nrow(C), J = ncol(C))
    
    # Specify initial values
    Nst <- apply(C, 1, max)       # Avoid data/model/inits conflict
    inits <- function(){list(N = Nst)}
    
    nmix.sim.out <- nimbleMCMC(code = Section6p3_code, 
                                 constants = win.data, 
                                 inits = inits,
                                 monitors = params,
                                 niter = ni, 
                                 nburnin = nb,
                                 nchains = nc,
                                 samplesAsCodaMCMC = TRUE)
    nmix.sim[[i]] <- MCMCsummary(nmix.sim.out, round = 4)
  }
  
  save("nmix.sim",file=paste0("out/nmixsim/out.nmix.",simname,".RData"))
  return(nmix.sim)
  
}


# Run nimble from R with output as Coda MCMC object
# run 25 models for each simulated data set
# J = 100, 50
# lambda = 2, 1, 0.5
# p = 0.5, 0.4, 0.3

# Sim01 = J = 100, lambda = 2, p = 0.5
# Sim02 = J = 100, lambda = 2, p = 0.4
# Sim03 = J = 100, lambda = 2, p = 0.3
# Sim04 = J = 100, lambda = 1, p = 0.5
# Sim05 = J = 100, lambda = 1, p = 0.4
# Sim06 = J = 100, lambda = 1, p = 0.3
# Sim07 = J = 100, lambda = 0.5, p = 0.5
# Sim08 = J = 100, lambda = 0.5, p = 0.4
# Sim09 = J = 100, lambda = 0.5, p = 0.3
# Sim11 = J = 50, lambda = 2, p = 0.5
# Sim12 = J = 50, lambda = 2, p = 0.4
# Sim13 = J = 50, lambda = 2, p = 0.3
# Sim14 = J = 50, lambda = 1, p = 0.5
# Sim15 = J = 50, lambda = 1, p = 0.4
# Sim16 = J = 50, lambda = 1, p = 0.3
# Sim17 = J = 50, lambda = 0.5, p = 0.5
# Sim18 = J = 50, lambda = 0.5, p = 0.4
# Sim19 = J = 50, lambda = 0.5, p = 0.3


###--- Simulation runs
# for(i in 1:10){
#   (sum(rpois(100,2)))
# }
# 
# mean(204, 207, 211, 184, 232, 210, 191, 227, 180, 183)

nmix.Sim01 <- nmix.function(simname=c("Sim01"), J=100, lambda=2, p=0.5)
nmix.Sim02 <- nmix.function(simname=c("Sim02"), J=100, lambda=2, p=0.4)
nmix.Sim03 <- nmix.function(simname=c("Sim03"), J=100, lambda=2, p=0.3)
nmix.Sim04 <- nmix.function(simname=c("Sim04"), J=100, lambda=1, p=0.5)
nmix.Sim05 <- nmix.function(simname=c("Sim05"), J=100, lambda=1, p=0.4)
nmix.Sim06 <- nmix.function(simname=c("Sim06"), J=100, lambda=1, p=0.3)
nmix.Sim07 <- nmix.function(simname=c("Sim07"), J=100, lambda=0.5, p=0.5)
nmix.Sim08 <- nmix.function(simname=c("Sim08"), J=100, lambda=0.5, p=0.4)
nmix.Sim09 <- nmix.function(simname=c("Sim09"), J=100, lambda=0.5, p=0.3)


# for some reason the J=50 sims didn't run - not sure why but letting it go for now
# nmix.Sim11 <- nmix.function(simname=c("Sim11"), J=50, lambda=2, p=0.5)
# nmix.Sim12 <- nmix.function(simname=c("Sim12"), J=50, lambda=2, p=0.4)
# nmix.Sim13 <- nmix.function(simname=c("Sim13"), J=50, lambda=2, p=0.3)
# nmix.Sim14 <- nmix.function(simname=c("Sim14"), J=50, lambda=1, p=0.5)
# nmix.Sim15 <- nmix.function(simname=c("Sim15"), J=50, lambda=1, p=0.4)
# nmix.Sim16 <- nmix.function(simname=c("Sim16"), J=50, lambda=1, p=0.3)
# nmix.Sim17 <- nmix.function(simname=c("Sim17"), J=50, lambda=0.5, p=0.5)
# nmix.Sim18 <- nmix.function(simname=c("Sim18"), J=50, lambda=0.5, p=0.4)
# nmix.Sim19 <- nmix.function(simname=c("Sim19"), J=50, lambda=0.5, p=0.3)


###---
# Load simulated runs
# forgot to save actual N - should be ~200
load("out/nmixsim/out.nmix.Sim01.RData")
list.sims <- list.files("out/nmixsim/")

nmix.sim.out <- vector('list', length(list.sims))
for(i in 1:length(list.sims)){
  load(paste0("out/nmixsim/",list.sims[i]))
  nmix.sim.out[[i]] <- nmix.sim
  }

nmix.sim.df <- as.data.frame(unlist(nmix.sim.out))
nrow(nmix.sim.df)
head(nmix.sim.df)
nmix.sim.df[1:14,]
colnames(nmix.sim.df)[1] <- c("value")
nmix.sim.df$estimate <- rep(c("mean","sd","CI_2.5","CI_50","CI_97.5","Rhat","n.eff"), each=2, time=25*length(list.sims))
nmix.sim.df$param <- rep(c("lambda","p"),each=1, time=7*25*length(list.sims))
nmix.sim.df$Sim <- rep(paste0("Sim0",seq_len(length(list.sims))),each=14*25)
nmix.sim.df$Run <- rep(paste0("Run", seq_len(25)),each=14, time=length(list.sims))
nmix.sim.wide <- pivot_wider(nmix.sim.df, names_from = estimate, values_from = value)

# plot simulations
# keep in mind
# Sim01 = J = 100, lambda = 2, p = 0.5
# Sim02 = J = 100, lambda = 2, p = 0.4
# Sim03 = J = 100, lambda = 2, p = 0.3
# Sim04 = J = 100, lambda = 1, p = 0.5
# Sim05 = J = 100, lambda = 1, p = 0.4
# Sim06 = J = 100, lambda = 1, p = 0.3
# Sim07 = J = 100, lambda = 0.5, p = 0.5
# Sim08 = J = 100, lambda = 0.5, p = 0.4
# Sim09 = J = 100, lambda = 0.5, p = 0.3
glimpse(nmix.sim.wide)

#- lambda
nmix.sim.plot.lambda <- ggplot(data = nmix.sim.wide[nmix.sim.wide$param=="lambda",]) +
  theme_bw() + theme(strip.background = element_rect(fill = "white", colour = "white")) +
  theme(panel.grid = element_blank())+
  geom_point(aes(x = Run, y = mean*100), size=2) +
  geom_hline(yintercept = c(200, 100, 50), col="grey") +
  geom_linerange(aes(x = Run, y = mean*100, ymin=CI_2.5*100, ymax= CI_97.5*100)) +
  theme(axis.text.x = element_blank()) +
  xlab("Simulation Runs") +
  ylab("Total Estimated Abundance (i.e., lambda * 100)")+
  ggtitle("Simulations of n-mixture models with 100 sampling occasions,\nlambda=c(2,1,0.5), and p=c(0.5,0.4,0.3)")+
  facet_wrap(~ Sim)

Cairo(file="out/nmix.sim.plot.lambda.PNG",
      type="png",
      width=3000,
      height=2200,
      pointsize=15,
      bg="white",
      dpi=300)
nmix.sim.plot.lambda
dev.off()

#- p
nmix.sim.plot.p <- ggplot(data = nmix.sim.wide[nmix.sim.wide$param=="p",]) +
  theme_bw() + theme(strip.background = element_rect(fill = "white", colour = "white")) +
  theme(panel.grid = element_blank())+
  geom_point(aes(x = Run, y = mean), size=2) +
  geom_hline(yintercept = c(0.5, 0.4, 0.3), col="grey") +
  geom_linerange(aes(x = Run, y = mean, ymin=CI_2.5, ymax= CI_97.5)) +
  theme(axis.text.x = element_blank()) +
  xlab("Simulation Runs") +
  ylab("Estimated Probabilty of Detection")+
  ggtitle("Simulations of n-mixture models with 100 sampling occasions,\nlambda=c(2,1,0.5), and p=c(0.5,0.4,0.3)")+
  facet_wrap(~ Sim)

Cairo(file="out/nmix.sim.plot.p.PNG",
      type="png",
      width=3000,
      height=2200,
      pointsize=15,
      bg="white",
      dpi=300)
nmix.sim.plot.p
dev.off()


nmix.sim.df[1:25,] 
mean.lambda <- 0.9385
M*mean.lambda # estimated N
sum(N) # compared to true N

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

# e2dist from scrbook
e2dist <- function (x, y) {
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}


###--- simulate the data
N <- round(area[1]*15) # population based on expected density
N # 44

K <- 4 # number of occasions

# simulate activity centres
sx1 <- runif(N, xlim[1,1],xlim[1,2])
sy1 <- runif(N, ylim[1,1],ylim[1,2])
S1 <- cbind(sx1, sy1)

sx2 <- runif(N, xlim[2,1],xlim[2,2])
sy2 <- runif(N, ylim[2,1],ylim[2,2])
S2 <- cbind(sx2, sy2)

# compute distance matrix
D1 <- e2dist(S, traps.C1) # distance of each individual from each trap
D2 <- e2dist(S, traps.C2) # distance of each individual from each trap

# Parameter values
p0 <- 0.1   # define parameters of encounter probability
sigma <- 1.5 # scale parameter of half-normal
alpha1 <- 1/(2*sigma*sigma) # convert to coefficient on distance

# Compute probability of encounter
probcap1 <- plogis(-2.5)*exp(-alpha1*D1*D1)
probcap2 <- plogis(-2.5)*exp(-alpha1*D2*D2)

# Generate the encounters of every individual in every trap
ntraps <- 20


Y1 <- matrix(NA, nrow=N, ncol=ntraps)
for(i in 1:nrow(Y1)){
  Y1[i,] <- rbinom(ntraps,K,probcap1[i,])
}

Y2 <- matrix(NA, nrow=N, ncol=ntraps)
for(i in 1:nrow(Y2)){
  Y2[i,] <- rbinom(ntraps,K,probcap2[i,])
}

sum(Y1)
sum(Y2)
sum(N)

# View(simSCR0)
# function to simulate data
# function (N = 100, K = 20, alpha0 = -2.5, sigma = 0.5, discard0 = TRUE, 
#           array3d = FALSE, rnd = NULL) 
# {
#   if (!is.null(rnd)) 
#     set.seed(rnd)
#   traplocs <- cbind(sort(rep(1:5, 5)), rep(1:5, 5))
#   Dmat <- e2dist(traplocs, traplocs)
#   ntraps <- nrow(traplocs)
#   plot(traplocs)
#   buffer <- 2
#   Xl <- min(traplocs[, 1] - buffer)
#   Xu <- max(traplocs[, 1] + buffer)
#   Yl <- min(traplocs[, 2] - buffer)
#   Yu <- max(traplocs[, 2] + buffer)
#   sx <- runif(N, Xl, Xu)
#   sy <- runif(N, Yl, Yu)
#   S <- cbind(sx, sy)
#   D <- e2dist(S, traplocs)
#   alpha1 <- 1/(2 * sigma * sigma)
#   probcap <- plogis(alpha0) * exp(-alpha1 * D * D)
#   Y <- matrix(NA, nrow = N, ncol = ntraps)
#   for (i in 1:nrow(Y)) {
#     Y[i, ] <- rbinom(ntraps, K, probcap[i, ])
#   }
#   if (discard0) {
#     totalcaps <- apply(Y, 1, sum)
#     Y <- Y[totalcaps > 0, ]
#   }
#   dimnames(Y) <- list(1:nrow(Y), paste("trap", 1:ncol(Y), 
#                                        sep = ""))
#   if (array3d) {
#     Y <- array(NA, dim = c(N, ntraps, K))
#     for (i in 1:nrow(Y)) {
#       for (j in 1:ntraps) {
#         Y[i, j, 1:K] <- rbinom(K, 1, probcap[i, j])
#       }
#     }
#     if (discard0) {
#       Y2d <- apply(Y, c(1, 2), sum)
#       ncaps <- apply(Y2d, 1, sum)
#       Y <- Y[ncaps > 0, , ]
#     }
#   }
#   list(Y = Y, traplocs = traplocs, xlim = c(Xl, Xu), ylim = c(Yl,Yu), N = N, alpha0 = alpha0, alpha1 = alpha1, sigma = sigma, 
#        K = K)
# }



# get rid of zeros so observed animals come first
Y1 = Y1[which(apply(Y1,1,sum)>0),]
dim(Y1) # 1 row for each observed animal

Y2 = Y2[which(apply(Y2,1,sum)>0),]
dim(Y2) # 1 row for each observed animal


y_sim <- array(0,c(M,J,2)) # needs to be an array, not a list

n0 = c(length(which(apply(Y1,1,sum)>0)),length(which(apply(Y2,1,sum)>0)))

y_sim[1:n0[1],,1] <- Y1
y_sim[1:n0[2],,2] <- Y2
dim(y_sim)
sum(y_sim)

#####################################################################################
# Model code and constants
# Bayesian analysis of the model using NIMBLE:

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
    if(sum(y_sim[i,,g])==1){
      st[i,1:2,g] = traps[y_sim[i,,g],,g]
    }else {
      st[i,1:2,g] = apply(traps[y_sim[i,,g],,g], 2, mean)
    }
  }
  for(i in (n0[g]+1):M){
    st[i,1:2,g] = runif(2, 0, 20)
  }}  

data <- list(
  y = y_sim, traps = traps, trials=rep(4,J))

z.init = apply(y_sim, c(1,3), sum)
z.init = ifelse(z.init >=1, 1, 0)

inits = list(z = z.init, p0 = runif(1,0.05, 1), psi = mean(z.init), sigma = runif(1, 2, 5), s = st)

params <- c('sigma', 'p0', 'psi', 'N', 'D')

# MCMC settings to test
# ni <- 250   ;   nt <- 1  ;   nb <- 50   ;   nc <- 1
# MCMC settings for actual run
ni <- 50000   ;   nt <- 20   ;   nb <- 5000   ;   nc <- 3

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
sim1 <- runMCMC(CscrMCMC, niter = ni, nburnin=nb,thin=nt,nchains=nc)
toc()
# sim1 = ni = 50000 # 952.56/60 # 16 min
str(sim1)

MCMCsummary(sim1, round = 4)

chainsPlot(sim1,
           var = c("N", "D", "sigma"))

# mean      sd    2.5%     50%   97.5% Rhat n.eff
# D[1]  13.0977  5.5887  4.1129 12.3388 25.7058 1.00  1606
# D[2]  14.3724  5.3685  5.6884 13.6045 26.5288 1.00  1706
# N[1]  38.2142 16.3059 12.0000 36.0000 75.0000 1.00  1606
# N[2]  42.2576 15.7845 16.7250 40.0000 78.0000 1.00  1706
# p0     0.1084  0.0423  0.0427  0.1029  0.2046 1.00  2428
# psi    0.2025  0.0800  0.0702  0.1941  0.3808 1.00  1631
# sigma  1.7776  0.4545  1.3205  1.6801  2.8438 1.02   687


#####################################################################################
# Start sets of simulations

# MCMC settings for actual run
ni <- 50000   ;   nt <- 20   ;   nb <- 5000   ;   nc <- 3

# create function to run through simulations
scr.function <- function(simname=simname, xlims=xlim, ylims=ylim, traps1=traps.C1, traps2=traps.C2,
                         N=44, J=20, G=2, p0=0.5, sigma=1.5, K=4, M=200, numsim=25){
  
  SCR.sim <- vector('list', numsim)
  names(SCR.sim) <- paste0('SCR.sim', seq_along(SCR.sim))
  for(s in seq_along(SCR.sim)){
    
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
    
    Y1 = Y1[which(apply(Y1,1,sum)>0),]
    # dim(Y1) # 1 row for each observed animal
    
    Y2 = Y2[which(apply(Y2,1,sum)>0),]
    # dim(Y2) # 1 row for each observed animal
    
    y_sim <- array(0,c(M,J,2)) # needs to be an array, not a list
    
    n0 = c(length(which(apply(Y1,1,sum)>0)),length(which(apply(Y2,1,sum)>0)))
    
    y_sim[1:n0[1],,1] <- Y1
    y_sim[1:n0[2],,2] <- Y2
    # dim(y_sim)
    # sum(y_sim)
    # 
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
          st[i,1:2,g] = apply(traps[y_sim[i,,g],,g], 2, mean)
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
                              thin=nt,
                              samplesAsCodaMCMC = TRUE)
    SCR.sim[[s]] <- MCMCsummary(SCR.sim.out, round = 4)
  }
  
  save("SCR.sim",file=paste0("out/scrsim/out.SCR.",simname,".RData"))
  return(SCR.sim)
  
}


# Sim01 = N = 44, p0 = 0.5, sigma = 2
# Sim02 = N = 44, p0 = 0.4, sigma = 2
# Sim03 = N = 44, p0 = 0.3, sigma = 2
# Sim04 = N = 44, p0 = 0.5, sigma = 1.5
# Sim05 = N = 44, p0 = 0.4, sigma = 1.5
# Sim06 = N = 44, p0 = 0.3, sigma = 1.5
# Sim07 = N = 44, p0 = 0.5, sigma = 1
# Sim08 = N = 44, p0 = 0.4, sigma = 1
# Sim09 = N = 44, p0 = 0.3, sigma = 1

# running simulations with 2 clusters, each with 20 traps, open for 4 occasions
# real population at each trap = 44, varying p0 and sigma
scr.Sim01 <- scr.function(simname=c("Sim01"), N=44, J=20, G=2, K=4, p0=0.5, sigma=2, numsim=25)
scr.Sim02 <- scr.function(simname=c("Sim02"), N=44, J=20, G=2, K=4, p0=0.4, sigma=2, numsim=25)
scr.Sim03 <- scr.function(simname=c("Sim03"), N=44, J=20, G=2, K=4, p0=0.3, sigma=2, numsim=25)
scr.Sim04 <- scr.function(simname=c("Sim04"), N=44, J=20, G=2, K=4, p0=0.5, sigma=1.5, numsim=25)
scr.Sim05 <- scr.function(simname=c("Sim05"), N=44, J=20, G=2, K=4, p0=0.4, sigma=1.5, numsim=25)
scr.Sim06 <- scr.function(simname=c("Sim06"), N=44, J=20, G=2, K=4, p0=0.3, sigma=1.5, numsim=25)
scr.Sim07 <- scr.function(simname=c("Sim07"), N=44, J=20, G=2, K=4, p0=0.5, sigma=1, numsim=25)
scr.Sim08 <- scr.function(simname=c("Sim08"), N=44, J=20, G=2, K=4, p0=0.4, sigma=1, numsim=25)
scr.Sim08 <- scr.function(simname=c("Sim09"), N=44, J=20, G=2, K=4, p0=0.3, sigma=1, numsim=25)

