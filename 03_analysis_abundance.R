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
# 03_analysis_abundance.R
# script to run SCR models fit to live trap (2000) and hair snag (2020) data
# written by Joanna Burgar (Joanna.Burgar@gov.bc.ca) - 13-Oct-2021
#####################################################################################
version$major
version$minor
R_version <- paste0("R-",version$major,".",version$minor)

.libPaths(paste0("C:/Program Files/R/",R_version,"/library")) # to ensure reading/writing libraries from C drive
tz = Sys.timezone() # specify timezone in BC

# Load Packages
list.of.packages <- c("tidyverse","parallel","unmarked", "nimble","nimbleSCR","MCMCvis","coda","Cairo","basicMCMCplots","tictoc","bayesplot")

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

# source("03_analysis_prep.R") # if not working concurrently

# nm.area <- diff(marten.data$xlim)*diff(marten.data$ylim)/100	# Density reported per 100 sq km
# nm.area # 51.21623 100 km2 or 5122 km2 total area (is this true?)

# 1999/00 detections high for first 3 weeks and then drop off - unlikely for models to converge
# different sampling effort in 1999/00 - concentrating only on fisher and moving traps to target recapturing fisher
# 1998/99 might also be worth ignoring as trapping effort became much more focused on fisher
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

# #####################################################################################
# # Specify model in BUGS language
# # This corresponds to "model2.txt" in original AHM code.
# # without effort
# Section6p3_code <- nimbleCode( {
#   # Priors
#   # lambda ~ dgamma(0.001, 0.001) # came with generic code
#   lambda ~ dunif(0,100) # for an uninformative prior
#   p ~ dunif(0, 1)
#   # Likelihood
#   for (i in 1:M) {
#     N[i] ~ dpois(lambda)      # State model
#     for (j in 1:J) {
#       C[i,j] ~ dbin(p, N[i]) # Observation model
#     }
#   }
# })
# 
# 
# # MCMC settings
# ni <- 25000   ;   nt <- 20   ;   nb <- 5000   ;   nc <- 3
# 
# # Parameters monitored
# params_noeffort <- c("lambda", "p")
# 
# # Bundle data without effort
# 
# # list(y_week, y_day, week.effort, weeks.to.use, daylookup)
# retro_weekly_noeffort <- vector("list", length(retro.data.out))
# for(i in 1:length(retro.data.out)){
#   y_week <- retro.data.out[[i]][[1]]
#   # y_week <- as.matrix(y_week)
#   ndata <- list(C = y_week, M = nrow(y_week), J = ncol(y_week))
#   str(ndata)
#   
#   # Specify initial values
#   Nst <- apply(y_week, 1, max)       # Avoid data/model/inits conflict
#   inits <- function(){list(N = Nst)}
#   
#   out <- nimbleMCMC(code = Section6p3_code, 
#                     constants = ndata, 
#                     inits = inits,
#                     monitors = params_noeffort,
#                     nburnin = nb, 
#                     niter = ni,
#                     nchains = nc,
#                     samplesAsCodaMCMC = TRUE)
#   retro_weekly_noeffort[[i]] <- out
#   
# }
# 
# save(retro_weekly_noeffort, file = paste0("./out/retro_weekly_noeffort_mcmcoutput.Rda"))
# # load("out/retro_weekly_noeffort_mcmcoutput.Rda")
# 
# out96 <- MCMCsummary(retro_weekly_noeffort[[1]])
# out97 <- MCMCsummary(retro_weekly_noeffort[[2]])
# out98 <- MCMCsummary(retro_weekly_noeffort[[3]])
# out99 <- MCMCsummary(retro_weekly_noeffort[[4]])
# 
# (retro.data.out[[4]][1])
# retro.out <- rbind(out96, out97, out98, out99)
# retro.out$Year <- rep(c("1996/97","1997/98","1998/99","1999/20"), each=2, times=1)
# retro.out$param <- rep(c("lambda","p"), each=1, times=4)
# retro.out # n.eff low and Gelman high fo 1999/00 = did not converge and should remove from output
# 
# #- lambda
# nmix.retro.plot.lambda <- ggplot(data = retro.out[retro.out$param=="lambda" & retro.out$Year!="1999/20",]) +
#   theme_bw() + theme(strip.background = element_rect(fill = "white", colour = "white")) +
#   theme(panel.grid = element_blank())+
#   geom_point(aes(x = Year, y = mean*100), size=2) +
#   geom_linerange(aes(x = Year, y = mean*100, ymin=`2.5%`*100, ymax= `97.5%`*100)) +
#   xlab("Year") +
#   ylab("Total Estimated Abundance (i.e., lambda * 100)")+
#   ggtitle("Williston Basin marten abundance estimates;\nlive trap data fit to n-mixture models")
# 
# Cairo(file="out/nmix.retro.plot.lambda.PNG",
#       type="png",
#       width=3000,
#       height=2200,
#       pointsize=15,
#       bg="white",
#       dpi=300)
# nmix.retro.plot.lambda
# dev.off()
# 
# #- p
# nmix.retro.plot.pall <- ggplot(data = retro.out[retro.out$param=="p",]) +
#   theme_bw() + theme(strip.background = element_rect(fill = "white", colour = "white")) +
#   theme(panel.grid = element_blank())+
#   geom_point(aes(x = Year, y = mean), size=2) +
#   geom_linerange(aes(x = Year, y = mean, ymin=`2.5%`, ymax= `97.5%`)) +
#   xlab("Year") +
#   ylab("Probabilty of Detection")+
#   ggtitle("Williston Basin marten detection probabilities;\nlive trap data fit to n-mixture models")
# 
# Cairo(file="out/nmix.retro.plot.pall.PNG",
#       type="png",
#       width=3000,
#       height=2200,
#       pointsize=15,
#       bg="white",
#       dpi=300)
# nmix.retro.plot.pall 
# dev.off()
# 
# # excluding 1999/2000
# nmix.retro.plot.p <- ggplot(data = retro.out[retro.out$param=="p" & retro.out$Year!="1999/20",]) +
#   theme_bw() + theme(strip.background = element_rect(fill = "white", colour = "white")) +
#   theme(panel.grid = element_blank())+
#   geom_point(aes(x = Year, y = mean), size=2) +
#   geom_linerange(aes(x = Year, y = mean, ymin=`2.5%`, ymax= `97.5%`)) +
#   xlab("Year") +
#   ylab("Probabilty of Detection")+
#   ggtitle("Williston Basin marten detection probabilities;\nlive trap data fit to n-mixture models")
# 
# Cairo(file="out/nmix.retro.plot.p.PNG",
#       type="png",
#       width=3000,
#       height=2200,
#       pointsize=15,
#       bg="white",
#       dpi=300)
# nmix.retro.plot.p 
# dev.off()


#####################################################################################
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


###--- wrangle data
# try nmix for all traps, using the same 21 day sampling bins as for occupancy and recent data
load("out/rec.data.out.RData")
load("out/retro.data.out.RData")

# To run as multi-year model, treat each year as a factor, and create arrays with slice for each year
# Only do so for 1996, 1997, and 2019
y_21day_9697 <- as.matrix(retro.data.out[[1]]$y_21day)
y_21day_9798 <- as.matrix(retro.data.out[[2]]$y_21day)
y_21day_2019 <- as.matrix(rec.data.out$rec_ydata)
# dimensions are traps by occasions

# M = traps
# J = occasions
# K = years
M <- c(nrow(y_21day_9697),nrow(y_21day_9798),nrow(y_21day_2019)) # number of traps
J <- c(ncol(y_21day_9697),ncol(y_21day_9798),ncol(y_21day_2019)) # number of occasions
K <- 3 # try first with 3 years as omitting 1998/99 and 1999/00 but including 2019

y_all <- array(c(y_21day_9697, y_21day_9798, y_21day_2019), dim=c(max(M),max(J),K))
dim(y_all)
sum(y_all)

# create year covariate
# year <- array(0, dim=c(max(M),K))
# dim(year)
# year[,1] <- c(rep(1,M[1]),rep(NA,max(M)-M[1]))
# year[,2] <- c(rep(2,M[2]),rep(NA,max(M)-M[2]))
# year[,3] <- c(rep(3,M[3]),rep(NA,max(M)-M[3]))

# create effort covariate
effort_9697 <- as.matrix(retro.data.out[[1]]$effort.21days)
effort_9798 <- as.matrix(retro.data.out[[2]]$effort.21days)
effort_2019 <- as.matrix(rec.data.out$rec_effort)

# effort_all <- array(c(effort_9697, effort_9798, effort_9899), dim=c(max(M),max(J),K))
effort_all <- array(c(effort_9697,effort_9798, effort_2019), dim=c(max(M),max(J),K))
dim(effort_all)
sum(effort_all)

# I use the scale by 2 SDs here following Gelman (2006):
# simple function to standardize variables
std2=function(x){
    (x - mean(x,na.rm=TRUE))/(2*sd(x,na.rm=TRUE))
  }
effort_all = std2(effort_all)

# ### note that I just used your ".Rda" - note that R data files are now saved with the extention (".Rdata"), but the old .Rda still works
# load("out/nmix.effort.input.Rda")
# str(nmix.effort.input)
# constants = nmix.effort.input[[1]]
# ydata_all = nmix.effort.input[[2]][1:2] # I didn't need the year variable
# names(ydata_all)[1] = "y"

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

# MCMC settings to test
# ni <- 250   ;   nt <- 1  ;   nb <- 50   ;   nc <- 1
# MCMC settings
nc <- 3   ;   ni <- 50000   ;   nb <- 5000   #nt <- 10

out <- nimbleMCMC(code = nmix_effort, 
                  data=data,
                  constants = constants, 
                  inits = inits,
                  monitors = params,
                  nburnin = nb, 
                  niter = ni,
                  nchains = nc,
                  samplesAsCodaMCMC = TRUE)

save(out, file = paste0("out/nmix_21day_969719_effort_mcmcoutput.RData"))
# save(out, file = paste0("out/nmix_21day_9719_effort_mcmcoutput.RData"))
# load("out/nmix_21day_969719_effort_mcmcoutput.RData")

chainsPlot(out) # traceplots saved - dp is better!

str(out)
head(out)

MCMCsummary(out,round = 4)
nmix.mcmc.summary <- MCMCsummary(out,round = 4)
str(nmix.mcmc.summary)
nmix.mcmc.summary$param <- row.names(nmix.mcmc.summary)
write.csv(nmix.mcmc.summary,"out/nmix_output_969719.csv")

nmix.df.alpha <- nmix.mcmc.summary %>% filter(grepl("alpha",param))%>% select("param","2.5%","50%","97.5%")
nmix.df.alpha <- nmix.df.alpha %>% pivot_longer(!param, names_to="estimate",values_to="value" )

nmix.df.alpha <- nmix.df.alpha %>% group_by(estimate) %>% mutate(new = ifelse(param != 'alpha0', value + value[param == 'alpha0'], value) )
nmix.df.alpha <- nmix.df.alpha %>% dplyr::select(-value) %>% filter(param!="alpha0") %>% pivot_wider(names_from = estimate, values_from = new)

# back calculate for detection probability
nmix.df.alpha$effort_2sd <- c(2*sd(effort_all[,,1]),2*sd(effort_all[,,2]),2*sd(effort_all[,,3]))
nmix.df.alpha <- nmix.df.alpha %>% mutate(across(`2.5%`:`97.5%`, ~ boot::inv.logit(.x*effort_2sd)))
write.csv(nmix.df.alpha,"out/nmix_alpha_969719.csv")


varnames(out) <- c("alpha0","dp 1996-1997","dp 1997-1998","dp 2019-2020",
                   "beta0","1996-1997","1997-1998", "2019-2020" )
str(out)
posterior <- as.matrix(out)
dimnames(posterior)
color_scheme_set("teal")
str(posterior)

dimnames(posterior)
# get annual coefficients and backtransform from scaled effort variable
posterior[,2] <- posterior[,1]+posterior[,2]
posterior[,3] <- posterior[,1]+posterior[,3]
posterior[,4] <- posterior[,1]+posterior[,4]

posterior[,2] <- boot::inv.logit(posterior[,2]*(2*sd(effort_all[,,1])))
posterior[,3] <- boot::inv.logit(posterior[,3]*(2*sd(effort_all[,,2])))
posterior[,4] <- boot::inv.logit(posterior[,4]*(2*sd(effort_all[,,3])))

Cairo(file="out/nmix_dp969719.PNG",type="png",width=2800,height=2200,pointsize=14,bg="white",dpi=300)
# mcmc_intervals(posterior, pars = c("dp 1996-1997","dp 1997-1998","dp 2019-2020"))+ # year effects on latent abundance
#   labs(
#     title = "Annual Detection Probability (dp) Estimates",
#     subtitle = "N-mixture models fit to Williston Basin live trapping (retrospective) and hair snag (current) data"
#   )
mcmc_areas(posterior, pars = c("dp 1996-1997","dp 1997-1998","dp 2019-2020"))+ # year effects on latent abundance
  labs(
    title = "Estimated Detection Probability (dp)",
    subtitle = "N-mixture models fit to Williston Basin live trapping (retrospective) and hair snag (current) data"
  )
dev.off()

dimnames(posterior)
# get annual coefficients
posterior[,6] <- posterior[,5]+posterior[,6]
posterior[,7] <- posterior[,5]+posterior[,7]
posterior[,8] <- posterior[,5]+posterior[,8]

Cairo(file="out/nmix_Ng969719.PNG",type="png",width=2800,height=2200,pointsize=14,bg="white",dpi=300)
mcmc_areas(posterior, pars = c("1996-1997","1997-1998","2019-2020"), area_method = "equal area")+ # year effects on latent abundance
  labs(
    title = "Estimated Median Abundance per Grid Cell",
    subtitle = "N-mixture models fit to Williston Basin live trapping (retrospective) and hair snag (recent) data"
  )
dev.off()


# Calculate ESS effective sample size
# adjusted sample size = effective sample size
# if a "long" Markov chain has only generated a short effective sample size, consider a longer run
# apply(posterior, 2, effectiveSize)


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
ni <- 50000  ;   nb <- 5000   ;   nc <- 3 #  nt <- 20 

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
results3 <- runMCMC(CscrMCMC, niter = ni, nburnin=nb,nchains=nc, setSeed = 500)
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