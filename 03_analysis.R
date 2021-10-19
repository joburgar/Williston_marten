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
list.of.packages <- c("tidyverse","parallel", "nimbleEcology","nimble", "scrbook","nimbleSCR","basicMCMCplots","coda","Cairo")

# Check you have them and load them
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)
rm(list.of.packages, new.packages) # for housekeeping

################################################################################
# if issues running nimble, run the following line
writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")


############################--- RETROSPECTIVE DATA ---##########################
# Run N-mixutre model for unmarked live trap data
# https://cran.r-project.org/web/packages/nimbleEcology/vignettes/Introduction_to_nimbleEcology.html
# Adapting from https://rdrr.io/cran/nimbleEcology/man/dNmixture.html
# Also need to review Applied Hierarchical Modeling code (converted to NIMBLE and available on github)
# https://github.com/nimble-training/AHMnimble/blob/master/Chapter_6/Section_6p11_setup.R


# load data if not running concurrently
# load("out/MartenData_1996.Rda")
glimpse(marten.data)


# Set up constants and initial values for defining the model
len <- 5 # length of dataset
dat <- c(1,2,0,1,5) # A vector of observations
lambda <- 10 # mean abundance

# load
prob <- c(0.2, 0.3, 0.2, 0.1, 0.4) # A vector of detection probabilities

# Define code for a nimbleModel
nc <- nimbleCode({
  x[1:5] ~ dNmixture_v(lambda, prob = prob[1:5],
                       Nmin = -1, Nmax = -1, len = 5)
  
  lambda ~ dunif(0, 1000)
  
  for (i in 1:5) {
    prob[i] ~ dunif(0, 1)
  }
})

# Build the model
nmix <- nimbleModel(nc,
                    data = list(x = dat),
                    inits = list(lambda = lambda,
                                 prob = prob))
# Calculate log probability of data from the model
nmix$calculate()

############################--- CURRENT DATA ---###########################
############################--- NIMBLE CODE ---###########################

###--- SCR nimble code from Paul van Dam-Bates
SCR_bern_sex <- nimbleCode({
  sigma ~ dunif(0,1000) # uninformative prior
  # sigma ~ dgamma(20,20) # uninformative prior
  psi ~ dbeta(1,1)
  psex ~ dbeta(1,1)
  lambda ~ dunif(0,10)
  
  for(i in 1:M){
    sex[i] ~ dbern(psex)
    z[i] ~ dbern(psi)
    X[i,1]~dunif(xlim[1],xlim[2])
    X[i,2]~dunif(ylim[1],ylim[2])
    d2[i,1:J]<- (X[i,1]-traps[1:J,1])^2 + (X[i,2]-traps[1:J,2])^2
    # Because detectability didn't vary per trap night we aggregated!
    # For comparison to Poisson we will use a hazard half normal detection function.
    # That way we are using the same lambda, encounter rate.
    # See Augustine paper as he does this as well.
    p[i,1:J]<- z[i]*(1-exp(-lambda*exp(-d2[i, 1:J]/(2*sigma*sigma) )) )
    #From Daniel Turek in nimbleSCR package. Fast binomial! Avoids loopin.
    y[i,1:J] ~ dbinom_vector(size = trials[1:J], prob = p[i,1:J])
  }
  N <- sum(z[1:M])
  D <- N/area
})

SCR_bern <- nimbleCode({
  sigma ~ dunif(0,1000) # uninformative prior
  # sigma ~ dgamma(20,20) # uninformative prior
  psi ~ dbeta(1,1)
  lambda ~ dunif(0,10)
  
  for(i in 1:M){
    z[i] ~ dbern(psi)
    X[i,1]~dunif(xlim[1],xlim[2])
    X[i,2]~dunif(ylim[1],ylim[2])
    d2[i,1:J]<- (X[i,1]-traps[1:J,1])^2 + (X[i,2]-traps[1:J,2])^2
    p[i,1:J]<- z[i]*(1-exp(-lambda*exp(-d2[i, 1:J]/(2*sigma*sigma) )) )
    y[i,1:J] ~ dbinom_vector(size = trials[1:J], prob = p[i,1:J])
  }
  N <- sum(z[1:M])
  D <- N/area
})

############################--- COMPILING DATA ---###########################
###--- actual data
load("out/MartenData_2020.Rda")

J <- marten.hsdata$J
area <- marten.hsdata$area
xlim <- marten.hsdata$xlim
ylim <- marten.hsdata$ylim
traps <- marten.hsdata$traps
edf <- marten.hsdata$edf
trials <- rep(4,nrow(traps))
sex <- marten.hsdata$sex

M <- 1000
y <- array(0, dim = c(M, 106, 4))
# Add the captures as 1s with the good old cbind trick.
glimpse(edf)
y[cbind(edf$Animal_Num,edf$Grid_Num, edf$Occ)] <- 1
sum(y) == nrow(edf) 	# It worked right?

# Now let's speed it up by summing over all 4 nights for a binomial dist.
y_all <- apply(y, 1:2, sum)

constants <- list(
  J = nrow(traps),
  area = area,
  xlim = xlim,
  ylim = ylim,
  traps = traps,
  M = M,
  trials = rep(4,nrow(traps))
)


data.sex <- list(y = y_all, 
             z =  c(rep(1, max(edf$Animal_Num)), rep(NA, M-max(edf$Animal_Num))),
             sex = c(sex, rep(NA, M-max(edf$Animal_Num))))

data <- list(y = y_all, 
                 z =  c(rep(1, max(edf$Animal_Num)), rep(NA, M-max(edf$Animal_Num))))

###--- Simulated data
# simulation to see what results should look like
set.seed(10)

area # Density reported per 100 sq km

mu <- 5 # density per 100 km2
N <- rpois(1, mu*area) # generate population
N # 554

s <- data.frame(s.x = runif(N, xlim[1], xlim[2]),
                s.y = runif(N, ylim[1], ylim[2]))

sigma <- 1
lambda0 <- 1
J <- nrow(traps) # nb of traps
K <- length(unique(edf$Occ)) # nb capture occasions

M <- 1000

# Now let's speed it up by summing over all 4 nights for a binomial dist.
yy <- array(NA, c(N, J, K))
for(j in 1:J) {
  dist <- sqrt((traps.sc$x[j] - s$s.x)^2 + (traps.sc$y[j] - s$s.y)^2)
  lambda <- lambda0 * exp(-dist^2 / (2 * sigma^2))
  for(k in 1:K) {
    yy[,j,k] <- rpois(N, lambda)
  }
}
n_all <- apply(yy, 1:2, sum)
sum(n_all) # 75 animals in the detection matrix
dim(n_all)

M_rest <- array(0, c(M-nrow(n_all), ncol(n_all)))
dim(M_rest)

sim_all <- rbind(n_all, M_rest)
sum(sim_all) == sum(n_all) #TRUE

sim.data.sex <- list(y = sim_all, 
                 z =  c(rep(1, max(nrow(n_all))), rep(NA, M-max(nrow(n_all)))),
                 sex = rbinom(M, 1, 0.4))
sim.data <- list(y = sim_all, 
                     z =  c(rep(1, max(nrow(n_all))), rep(NA, M-max(nrow(n_all)))))

############################--- RUN MODEL ---###########################
# run for both simulated data and actual data
Rmodel <- nimbleModel(SCR_bern_sex, constants, sim.data.sex)
Rmodel <- nimbleModel(SCR_bern, constants, sim.data)

conf <- configureMCMC(Rmodel)
# conf$setMonitors(c('sigma', 'lambda', 'psi', 'N', 'D', 'psex'))
conf$setMonitors(c('sigma', 'lambda', 'psi', 'N', 'D'))

# Use a block update on locations. Saves time.
conf$removeSamplers('X')
for(i in 1:M) conf$addSampler(target = paste0('X[', i, ', 1:2]'), 
                              type = 'RW_block', silent = TRUE)

Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run compiled model. Takes ~11 minutes on my machine (M=300, niter=10000) or ~40 minutes (M=500, niter=20000).
print(Sys.time())
samplesList <- runMCMC(Cmcmc,
                       niter = 10000,
                       nburnin = 5000,
                       nchains = 3)
print(Sys.time())

samples <- rbind(samplesList[[1]],
                 samplesList[[2]],
                 samplesList[[2]])

str(samples)

# Calculate ESS effective sample size
# adjusted sample size = effective sample size
# if a "long" Markov chain has only generated a short effective sample size, consider a longer run
apply(samples, 2, effectiveSize)

# Produce trace an density plots.
chainsPlot(samplesList,
           var = c("N", "sigma", "lambda"))

chainsPlot(samplesList,
           var = c("D", "psi","psex"))
# Display summary stats. Compare to the values used to simulate data, in particular N=161, σ=0.5 and λ0=2.
summary(samples)





Rmodel <- nimbleModel(SCR_bern, constants, sim.data)

conf <- configureMCMC(Rmodel)
conf$setMonitors(c('sigma', 'lambda', 'psi', 'N', 'D', 'psex'))

# Use a block update on locations. Saves time.
conf$removeSamplers('X')
for(i in 1:M) conf$addSampler(target = paste0('X[', i, ', 1:2]'), 
                              type = 'RW_block', silent = TRUE)

Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run compiled model. Takes ~11 minutes on my machine (M=300, niter=10000) or ~40 minutes (M=500, niter=20000).
print(Sys.time())
samplesList <- runMCMC(Cmcmc,
                       niter = 10000,
                       nburnin = 5000,
                       nchains = 3)
print(Sys.time())

samples <- rbind(samplesList[[1]],
                 samplesList[[2]],
                 samplesList[[2]])

str(samples)

# Calculate ESS effective sample size
# adjusted sample size = effective sample size
# if a "long" Markov chain has only generated a short effective sample size, consider a longer run
apply(samples, 2, effectiveSize)

# Produce trace an density plots.
chainsPlot(samplesList,
           var = c("N", "sigma", "lambda"))

chainsPlot(samplesList,
           var = c("D", "psi","psex"))
# Display summary stats. Compare to the values used to simulate data, in particular N=161, σ=0.5 and λ0=2.
summary(samples)




# # https://esajournals.onlinelibrary.wiley.com/doi/10.1002/ecs2.3385
# # https://cran.r-project.org/web/packages/nimbleSCR/vignettes/Simulate_and_fit_SCR_models_with_dbinomLocal_normal.html
# # Follow similar example as above - create discretized habitat grid represented by grid cell centers
# # mask out areas not suitable (i.e., water)
# # consider relatively small cell size of 'marten cell' as per trap grid set up
# # think it's the grid cell that either has or doesn't have a detection (need to confirm)
# # above would work considering the hair snag traps did move within a grid between occasions
# 
# # Define the model
# code <- nimbleCode({
#   # sigma ~ dgamma(20, 20) # centred on ~ 500 ha home range (00_SC_Informed_Sigma.R)
#   sigma ~ dunif(0, 10) # vague prior
#   lam0 ~ dunif(0, 10)
#   psi ~ dbeta(1, 1)
#   for(i in 1:M) {
#     z[i] ~ dbern(psi)
#     s[i,1] ~ dunif(xlim[1], xlim[2])
#     s[i,2] ~ dunif(ylim[1], ylim[2])
#     dist[i,1:J] <- (s[i,1] - X[1:J,1])^2 + (s[i,2] - X[1:J,2])^2
#     lam[i,1:J] <- exp(-dist[i,1:J] / (2 * sigma^2)) * z[i]
#   }
#   for(j in 1:J){
#     bigLambda[j] <- lam0 * sum(lam[1:M,j])
#     for(k in 1:K) {
#       n[j,k] ~ dpois(bigLambda[j])
#     }
#   }
#   N <- sum(z[1:M])
# })
# 








SCR_bern <- nimbleCode({
  sigma ~ dunif(0,1000) # uninformative prior
  # sigma ~ dgamma(20,20) # uninformative prior
  psi ~ dbeta(1,1)
  lambda ~ dunif(0,10)
  
  for(i in 1:M){
    z[i] ~ dbern(psi)
    X[i,1]~dunif(xlim[1],xlim[2])
    X[i,2]~dunif(ylim[1],ylim[2])
    d2[i,1:J]<- (X[i,1]-traps[1:J,1])^2 + (X[i,2]-traps[1:J,2])^2
    
    p[i,1:J]<- z[i]*(1-exp(-lambda*exp(-d2[i, 1:J]/(2*sigma*sigma) )) )
    y[i,1:J] ~ dbinom_vector(size = trials[1:J], prob = p[i,1:J])
  }
  N <- sum(z[1:M])
  D <- N/area
})


# 
# 
# # Define constants, data and inits
# glimpse(marten.hsdata)
# 
# # run for simulated data
# M <- 200 # data augmented population, aim for double expected population
# # J <- nrow(marten.hsdata$traps) # nb of traps
# J <- nrow(traps.sc)
# # K <- ncol(marten.hsdata$observations) # nb capture occasions
# K <- ncol(n)
# constants <- list(M = M, 
#                   K = K, 
#                   J = J)
# 
# 
# data <- list(n = n, 
#              X = traps.sc, 
#              xlim = xlim, 
#              ylim = ylim)
# 
# 
# s <- cbind(runif(M, xlim[1], xlim[2]), 
#            runif(M, ylim[1], ylim[2]))
# 
# z <- rep(1, M)
# 
# inits <- list(sigma = 0.5, 
#               lam0 = 1, 
#               s = s, 
#               z = z,
#               psi = 0.5)
# 
# # Build R model (not compiled yet)
# Rmodel <- nimbleModel(code = code, 
#                       constants = constants, 
#                       data = data, 
#                       inits = inits)
# 
# # Check whether the model is fully initialized
# # If you failed at providing initial values for some parameters (e.g. ψ), you’ll get NAs.
# Rmodel$calculate()
# 
# # Now compile the model in C++.
# # if issues check this: https://cran.r-project.org/bin/windows/Rtools/#putting-rtools-on-the-path
# Cmodel <- compileNimble(Rmodel)
# 
# # The R and C models are exactly the same versions of the model.
# calculate(Cmodel)
# 
# # You can simulate from prior
# Cmodel$simulate('lam0')
# calculate(Cmodel)
# 
# # Specify MCMC.
# conf <- configureMCMC(Rmodel,
#                       monitors = c("N", "sigma","lam0", "psi"))
# 
# # Build an executable MCMC.
# Rmcmc <- buildMCMC(conf)
# 
# # Compile in C++.
# Cmcmc <- compileNimble(Rmcmc, project = Cmodel)
# 
# # Run compiled model (do not run the uncompiled model with runMCMC(Rmcmc,100)).
# samples <- runMCMC(Cmcmc,100)
# 
# # Explore
# dim(samples)
# colnames(samples)
# samplesSummary(samples)
# 
# # Run compiled model. Takes ~15-30 minutes on my machine (M=200 vs M=300).
# print(Sys.time())
# 
# samplesList <- runMCMC(Cmcmc,
#                        niter = 10000,
#                        nburnin = 5000,
#                        nchains = 3)
# 
# print(Sys.time())
# 
# samples <- rbind(samplesList[[1]],
#                  samplesList[[2]],
#                  samplesList[[2]])
# 
# str(samples)
# 
# # Calculate ESS effective sample size
# # adjusted sample size = effective sample size
# # if a "long" Markov chain has only generated a short effective sample size, consider a longer run
# apply(samples, 2, effectiveSize)
# 
# # Produce trace an density plots.
# chainsPlot(samplesList,
#            var = c("D", "sigma", "lambda"))
# 
# # Display summary stats. Compare to the values used to simulate data, in particular N=161, σ=0.5 and λ0=2.
# summary(samples)
# # N              lam0            psi             sigma       
# # Min.   : 46.0   Min.   :1.944   Min.   :0.1464   Min.   :0.3247  
# # 1st Qu.:169.0   1st Qu.:4.297   1st Qu.:0.5611   1st Qu.:0.4093  
# # Median :213.0   Median :5.317   Median :0.7100   Median :0.4354  
# # Mean   :208.1   Mean   :5.572   Mean   :0.6927   Mean   :0.4380  
# # 3rd Qu.:251.0   3rd Qu.:6.631   3rd Qu.:0.8380   3rd Qu.:0.4641  
# # Max.   :300.0   Max.   :9.994   Max.   :0.9999   Max.   :0.5783
# 


# library(parallel)
# library(coda)
# nc <- 3 # number of chains
# cl<-makeCluster(nc,timeout=5184000)
# inits <- function() {list(sigma = 0.5, 
#                           lam0 = 0.1, 
#                           s = s, 
#                           z = z,
#                           psi = 0.5)}
# myCode <- code
# nimbledata <- data
# nimbleconstants <- constants
# params <- c("N", "psi", "sigma", "lam0")
# clusterExport(cl, c("myCode", "inits", "nimbledata", "nimbleconstants", "params"))
# for (j in seq_along(cl)) {
#   set.seed(j)
#   init <- inits()
#   clusterExport(cl[j], "init")
# }
# out <- clusterEvalQ(cl, {
#   library(nimble)
#   library(coda)
#   model <- nimbleModel(code = myCode, name = "myCode",
#                        constants = nimbleconstants, data = nimbledata,
#                        inits = init)
#   Cmodel <- compileNimble(model)
#   modelConf <- configureMCMC(model)
#   modelConf$addMonitors(params)
#   modelMCMC <- buildMCMC(modelConf)
#   CmodelMCMC <- compileNimble(modelMCMC, project = myCode)
#   out1 <- runMCMC(CmodelMCMC, niter = 100) #10000
#   return(as.mcmc(out1))
# })
# 
# out.mcmc <- as.mcmc(out)
# traceplot(out.mcmc[, "N"]
#           
#           ## If has not converged, continue sampling
#           start <- Sys.time()
#           out2 <- clusterEvalQ(cl, {
#             out1 <- runMCMC(CmodelMCMC, niter = 20000)
#             return(as.mcmc(out1))
#           })
#           
#           out.mcmc.update1 <- as.mcmc(out2)
#           
#           out.mcmc.bind <- mcmc.list()
#           for (i in seq_len(nc)) {
#             out.mcmc.bind[[i]] <- mcmc(rbind(out.mcmc[[i]], out.mcmc.update1[[i]]))
#           }
#           traceplot(out.mcmc.bind[, "N"]
