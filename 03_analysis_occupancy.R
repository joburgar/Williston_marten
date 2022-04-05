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
list.of.packages <- c("tidyverse","unmarked","Cairo","tictoc","AICcmodavg", "MuMIn","PNWColors","sf")

# Check you have them and load them
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)
rm(list.of.packages, new.packages) # for housekeeping

############################--- LOAD DATA ---##########################

load("out/retro.data.out9697.RData")
load("out/rec.data.out.RData")

str(retro.data.out)

cov.df <- read.csv("data/covdata.csv",row.names=1)
head(cov.df)
# cov.data.files <- c("data/aoi.covdata.1996.csv","data/retro.covdata.1996.csv",
#                     "data/aoi.covdata.1997.csv","data/retro.covdata.1997.csv",
#                     "data/aoi.covdata.2020.csv","data/rec.covdata.2020.csv")
# 
# ALLcov.df <-do.call(rbind,lapply(cov.data.files, read.csv))
# ALLcov.df$X <- ALLcov.df$grid_id <- NULL
# cov.df.nrows <- sapply( cov.data.files, function(f) nrow(read.csv(f)) )
# 
# ALLcov.df$Year <- c(rep(1996,cov.df.nrows[1]+cov.df.nrows[2]),
#                     rep(1997,cov.df.nrows[3]+cov.df.nrows[4]), 
#                     rep(2019,cov.df.nrows[5]+cov.df.nrows[6]))
# ALLcov.df$Area <- c(rep("Area",cov.df.nrows[1]), rep("Traps",cov.df.nrows[2]),
#                     rep("Area",cov.df.nrows[3]), rep("Traps",cov.df.nrows[4]),
#                     rep("Area",cov.df.nrows[5]), rep("Traps",cov.df.nrows[6]))
# ALLcov.df %>% group_by(Year) %>% count(Area)
# ALLcov.df$Area_Year <- paste(ALLcov.df$Area, ALLcov.df$Year)


cov.df <- cov.df %>% filter(grid_touse==1) %>% select(-use, -grid_touse,-RLW_type)
head(cov.df)
cov.df <- cov.df %>% rename("TRP_DNSTY_96"=TRP_DNSTY_9195, "TRP_DNSTY_97"=TRP_DNSTY_9296, "TRP_DNSTY_19"=TRP_DNSTY_1418)

names.cov.df <- names(cov.df)
col96 <- names.cov.df[grepl("retro|96",names.cov.df)]
col96 <- col96[!grepl("97",col96)]
cov_retro_96 <- cov.df %>% dplyr::select("grid_id","SBS_prop","RLW_dist",all_of(col96))
cov_retro_96$Year <- "1996-1997"

col97 <- names.cov.df[grepl("retro|97",names.cov.df)]
col97 <- col97[!grepl("96",col97)]
cov_retro_97 <- cov.df %>% dplyr::select("grid_id","SBS_prop","RLW_dist",all_of(col97))
cov_retro_97$Year <- "1997-1998"

col19 <- names.cov.df[grepl("recent|19",names.cov.df)]
cov_recent_19 <- cov.df %>% dplyr::select("grid_id","SBS_prop","RLW_dist",all_of(col19))
cov_recent_19$Year <- "2019-2020"

names(cov_retro_96);names(cov_retro_97);names(cov_recent_19)
new.names.cov.df <- c("grid_id","SBS_prop","RLW_dist","RD_density","TREE20_prop","CANOPY_prop",
                      "EDGE_density","HARVEST_prop","Area","TRP_DNSTY","Year")
names(cov_retro_96) <- names(cov_retro_97) <- names(cov_recent_19) <- new.names.cov.df

ALLcov.df <- rbind(cov_retro_96, cov_retro_96 %>% filter(Area==1), 
                   cov_retro_97, cov_retro_97 %>% filter(Area==1),
                   cov_recent_19, cov_recent_19 %>% filter(Area==1))

nrow(cov_retro_96 %>% filter(Area==1)); nrow(cov_retro_97%>% filter(Area==1)); nrow(cov_recent_19 %>% filter(Area==1))
nrow(cov_retro_96); nrow(cov_retro_97); nrow(cov_recent_19)

ALLcov.df$Area_Year <- c(rep("Area 1996",times=nrow(cov_retro_96)),
                         rep("Traps 1996", times=nrow(cov_retro_96 %>% filter(Area==1))),
                         rep("Area 1997",times=nrow(cov_retro_97)),
                         rep("Traps 1997", times=nrow(cov_retro_97 %>% filter(Area==1))),
                         rep("Area 2019",times=nrow(cov_recent_19)),
                         rep("Traps 2019", times=nrow(cov_recent_19 %>% filter(Area==1))))

ALLcov.df %>% group_by(Area_Year) %>% count(Area,Year)

ALLcov.df_longer <- ALLcov.df %>% dplyr::select(-Area,-Year,-grid_id) %>%
  pivot_longer(!(Area_Year), names_to = "Covariate", values_to = "Values")

names(ALLcov.df_longer)
unique(ALLcov.df_longer$Covariate)
ALLcov.df_longer$Covariate <- as.factor(recode(ALLcov.df_longer$Covariate, SBS_prop = "Proportion SBS", RLW_dist = "Distance to Water",
       RD_density = "Road Density", TREE20_prop = "Proportion Trees > 20 m", CANOPY_prop = "Proportion Canopy > 45%",
       EDGE_density = "Forest Edge Density", HARVEST_prop = "Proportion Harvested", TRP_DNSTY = "Marten Trap Harvest Density"))
levels(ALLcov.df_longer$Covariate)

ALLcov.df_longer$Covariate <- fct_relevel(ALLcov.df_longer$Covariate, "Proportion SBS", "Proportion Trees > 20 m","Proportion Canopy > 45%", "Distance to Water","Proportion Harvested",
            "Forest Edge Density", "Road Density", "Marten Trap Harvest Density")

ALLcov.df_longer$Area_Year <- fct_relevel(ALLcov.df_longer$Area_Year,  "Traps 1996","Area 1996", "Traps 1997","Area 1997","Traps 2019","Area 2019")

pal = pnw_palette(name="Cascades",n=6,type="discrete")

cov.plot <- ggplot(ALLcov.df_longer %>% filter(Covariate!="Marten Trap Harvest Density"),
                   aes(x=Covariate, y=Values, fill=as.factor(Area_Year))) +
  geom_boxplot(notch=FALSE) +
  facet_wrap(~Covariate, scale="free",nrow=2)+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank()) +
  theme(legend.position="bottom")+
  theme(legend.title=element_blank())+
  scale_fill_manual(values=pal)

Cairo(file="out/marten_ALLcov_plot_969719.PNG",type="png",width=3400,height=2400,pointsize=14,bg="white",dpi=300)
cov.plot
dev.off()

# Run occupancy models for unmarked and marked live trap data
# consider that male home ranges ~5.25 km2 and female home ranges ~3.16 km2 (Eric Lofroth's MSc thesis)
# larger in fire disturbed landscapes in BC/WA during winter, ~13 km2 for males and ~8.5 km2 for females  (Volkmann & Hodges 2021)

# 1999/00 detections high for first 3 weeks and then drop off - unlikely for models to converge
# different sampling effort in 1999/00 - concentrating only on fisher and moving traps to target recapturing fisher
# ignoring 1998/99 as well, as trapping effort became much more focused on fisher
#####################################################################################
# Following from Kery and Chandler for 'Dynamic occupancy models in unmarked'
# https://cran.r-project.org/web/packages/unmarked/vignettes/colext.pdf
# site = i; survey = j; season=t; yijt = 1 if detected
# replicate surveys at a site during a single season are indpendent (or dependency must be modeled)
# occurrence state does not change over replicate surveys at site i during season t
# there are no false positive errors

head(ALLcov.df)

ALLcov.df_grid_to_use <- ALLcov.df %>% filter(grepl("Traps", Area_Year)) %>% select(grid_id)
ALLcov.df_grid_to_use <- unique(ALLcov.df_grid_to_use$grid_id)
length(ALLcov.df_grid_to_use)

M <- length(ALLcov.df_grid_to_use) # 150; number of sites
J <- max(length(retro.data.out[[1]]$days21.to.use),  # number of sampling occasions
         length(retro.data.out[[2]]$days21.to.use),
         ncol(rec.data.out$rec_effort))
T <- 3 # number of years / seasons


# create site covariates
# 4 covariates that vary by site only
# geographic range (utm x and utm y), proportion of SBS, distance to water
load("data/aoi_grid.RData")
# ggplot()+
#   geom_sf(data=aoi_grid %>% filter(grid_id %in% ALLcov.df_grid_to_use), aes(fill=Yrs_surveyed))

grid_centroid_utm <- as.data.frame(st_coordinates(st_centroid(aoi_grid)))
# grid_centroid_utm$grid_id <- rownames(grid_centroid_utm)
# grid_centroid_to_use <- grid_centroid_utm %>% filter(grid_id %in% ALLcov.df_grid_to_use)

ALLcov.df$utm_x <- grid_centroid_utm$X[match(ALLcov.df$grid_id, rownames(grid_centroid_utm))]
ALLcov.df$utm_y <- grid_centroid_utm$Y[match(ALLcov.df$grid_id, rownames(grid_centroid_utm))]

colext.df <- ALLcov.df %>% filter(grid_id %in% ALLcov.df_grid_to_use) %>% filter(grepl("Area",Area_Year))
colext.df %>% count(Year) # 150 rows for each year
colext.df <- colext.df %>% arrange(Year, grid_id) # arrange by Year then by grid_id for consistency

# the 4 site covariates don't change by year so can be subset to just the 1996-1997 year values

colext.site.cov <- colext.df %>% filter(Year=="1996-1997") %>% dplyr::select(utm_x, utm_y, RLW_dist,SBS_prop)
summary(colext.site.cov)
# scale covariates for convergence - keep SBS as proportion 0-1 other ones with mean 0 and 1 SD scale
colext.site.covS <- cbind(as.data.frame(scale(colext.site.cov[,1:3])),colext.site.cov[4])
summary(colext.site.covS)


# annual covariates
# yearly site covariates (M rows and T columns)
year <- c("1996-1997","1997-1998","2019-2020")
year <- matrix(year, nrow(yy),3,byrow=TRUE)

names(colext.df) # recall that colext.df is arranged by grid_id and then by year
colext.annual.cov <- colext.df %>% dplyr::select(RD_density, EDGE_density, TRP_DNSTY, TREE20_prop, CANOPY_prop, HARVEST_prop)
colext.annual.cov <- colext.df %>% dplyr::select(RD_density, EDGE_density, TRP_DNSTY, TREE20_prop, CANOPY_prop, HARVEST_prop)
colext.annual.covS <- cbind(as.data.frame(scale(colext.annual.cov[,1:3])),colext.annual.cov[,4:6])
colext.annual.covS[is.na(colext.annual.covS)] <- 0 # change the trapping density to 0 (mean) for unknown pixels

y1 <- 1:nrow(yy)
y2 <- (nrow(yy)+1):(2*nrow(yy))
y3 <- (2*nrow(yy)+1):(3*nrow(yy))


create_annual_cov_matrix <- function(M=150, T=3, cov.df=cov.df, cov.name=cov.name){
  cov_array <- array(NA, dim=c(M, T))
  cov_array[,1] <- unlist(cov.df[y1,][c(cov.name)])
  cov_array[,2] <- unlist(cov.df[y2,][c(cov.name)])
  cov_array[,3] <- unlist(cov.df[y3,][c(cov.name)])
  return(cov_array)
}

names(colext.annual.covS)

RD_dnsty <- create_annual_cov_matrix(cov.df=colext.annual.covS, cov.name="RD_density")
EDGE_dnsty <- create_annual_cov_matrix(cov.df=colext.annual.covS, cov.name="EDGE_density")
TRP_dnsty <- create_annual_cov_matrix(cov.df=colext.annual.covS, cov.name="TRP_DNSTY")
TREE20_prop <- create_annual_cov_matrix(cov.df=colext.annual.covS, cov.name="TREE20_prop")
CANOPY_prop <- create_annual_cov_matrix(cov.df=colext.annual.covS, cov.name="CANOPY_prop")
HARVEST_prop <- create_annual_cov_matrix(cov.df=colext.annual.covS, cov.name="HARVEST_prop")


################################################################################
# detection data
y_all <- array(NA,dim=c(length(ALLcov.df_grid_to_use),0))
rownames(y_all) <- as.character(ALLcov.df_grid_to_use)

y96.df <- retro.data.out[[1]]$y_21day
y96 <- merge(y_all, y96.df, by="row.names", all.x=TRUE)
y96 <- arrange(y96, Row.names)
rownames(y96) <- y96$Row.names
y96$Row.names <- NULL
y96 <- as.matrix(y96)
y96[y96 > 0] <- 1 # change to 0 and 1 only


y97.df <- retro.data.out[[2]]$y_21day
y97 <- merge(y_all, y97.df, by="row.names", all.x=TRUE)
y97 <- arrange(y97, Row.names)
rownames(y97) <- y97$Row.names
y97$Row.names <- NULL
y97 <- as.matrix(y97)
y97[y97 > 0] <- 1 # change to 0 and 1 only


y19.df <- rec.data.out$rec_ydata
y19 <- merge(y_all, y19.df, by="row.names", all.x=TRUE)
y19 <- arrange(y19, Row.names)
rownames(y19) <- y19$Row.names
y19$Row.names <- NULL
y19 <- as.matrix(y19)
y19[y19 > 0] <- 1 # change to 0 and 1 only

y <- array(NA, dim=c(M,J,T))
y[,,1] <- cbind(y96, array(NA, dim=c(nrow(y96),J-ncol(y96))))
y[,,2] <- cbind(y97, array(NA, dim=c(nrow(y97),J-ncol(y97))))
y[,,3] <- cbind(y19, array(NA, dim=c(nrow(y19),J-ncol(y19))))

sum(y, na.rm=T); sum(sum(y96, na.rm=T),sum(y97, na.rm = T),sum(y19, na.rm = T))

yy <- matrix(y, M, J*T)

################################################################################
# effort data
e_all <- array(NA,dim=c(length(ALLcov.df_grid_to_use),0))
rownames(e_all) <- as.character(ALLcov.df_grid_to_use)

e96.df <- retro.data.out[[1]]$effort.21days
e96 <- merge(e_all, e96.df, by="row.names", all.x=TRUE)
e96 <- arrange(e96, Row.names)
rownames(e96) <- e96$Row.names
e96$Row.names <- NULL
e96 <- as.matrix(e96)

e97.df <- retro.data.out[[2]]$effort.21days
e97 <- merge(e_all, e97.df, by="row.names", all.x=TRUE)
e97 <- arrange(e97, Row.names)
rownames(e97) <- e97$Row.names
e97$Row.names <- NULL
e97 <- as.matrix(e97)

e19.df <- rec.data.out$rec_effort
e19 <- merge(e_all, e19.df, by="row.names", all.x=TRUE)
e19 <- arrange(e19, Row.names)
rownames(e19) <- e19$Row.names
e19$Row.names <- NULL
e19 <- as.matrix(e19)

e <- array(NA, dim=c(M,J,T))
e[,,1] <- cbind(e96, array(NA, dim=c(nrow(e96),J-ncol(e96))))
e[,,2] <- cbind(e97, array(NA, dim=c(nrow(e97),J-ncol(e97))))
e[,,3] <- cbind(e19, array(NA, dim=c(nrow(e19),J-ncol(e19))))

sum(e, na.rm=T); sum(sum(e96, na.rm=T),sum(e97, na.rm = T),sum(e19, na.rm = T))
ee <- matrix(e, M, J*T)

# same scaling as other covariates
# sd.ee <- sd(c(ee), na.rm=TRUE)
# mean.ee <- mean(ee, na.rm=TRUE)
# eeS <- (ee - mean.ee) / sd.ee

# use the scale by 2 SDs here following Gelman (2006) if want the same as n-mixture
# simple function to standardize variables
std2=function(x){
  (x - mean(x,na.rm=TRUE))/(2*sd(x,na.rm=TRUE))
}

eeS = std2(ee)
summary(eeS, na.rm=T)
################################################################################

###--- combine data for simple data frame - only year as yearly site covariates
sUMF <- unmarkedMultFrame(
  y = yy,
  yearlySiteCovs = list(year = year),
  numPrimary=T)
summary(sUMF)

# Model with all constant parameters
m0 <- colext(psiformula= ~1, gammaformula = ~ 1, epsilonformula = ~ 1,
                 pformula = ~ 1, data = sUMF, method="BFGS")
summary(m0)

names(m0)
# [1] "psi" "col" "ext" "det"
# occupancy, colonization, extinction, detection

# back-transform to original scale using the inverse-logit function (plogis)
# all parameters were estimated on the logit scale

plogis(coef(m0))
# psi(Int)     col(Int)     ext(Int)       p(Int) 
# 0.9999999326 0.0001749503 0.2651643237 0.2167176622

# can also use backTransform function to back transform estimate, se and confidence intervals
backTransform(m0, type="psi")
confint(backTransform(m0, type="psi"))

###--- dynamic occupancy model with full year-dependence in the parameters describing occupancy dynamics and also detection
# as year is a factor, this analysis is parameterized in terms of an intercept and effects representing differences
# this means that the parameter for the first year is the intercept and the effects denote differences between the parameter values in all other years,
# relative to the parameter value in the first year, which serves as a reference level
# a means parameterization may be more practical for simple presentation
# means parameterization can be specified by adding a -1 to the formula for the time-dependent parameters
m1 <- colext(psiformula = ~1, # First-year occupancy
               gammaformula = ~ year-1, # Colonization
               epsilonformula = ~ year-1, # Extinction
               pformula = ~ year-1, # Detection
               data = sUMF)
m1
backTransform(m1, type="psi")
confint(backTransform(m1, type="psi"))

nd <- data.frame(year=c("1996-1997","1997-1998"))
E.ext <- predict(m1, type='ext', newdata=nd)
E.col <- predict(m1, type='col', newdata=nd)
nd <- data.frame(year=c("1996-1997","1997-1998","2019-2020"))
E.det <- predict(m1, type='det', newdata=nd)
# plogis(coef(m1)) # checking to see same back transformation with predict and plogis = TRUE

op <- par(mfrow=c(3,1))
with(E.ext, { # Plot for extinction probability
  plot(1:2, Predicted, pch=1, xaxt='n', xlab='Year',
       ylab=expression(paste('Extinction probability ( ', epsilon, ' )')),
       ylim=c(0,1))
  axis(1, at=1:2, labels=nd$year[1:2])
  arrows(1:2, lower, 1:2, upper, code=3, angle=90, length=0.03, col=4)
  points((1:2), 1-phi, col=1, lwd = 1, pch=16)
  legend(7, 1, c('Parameter', 'Estimate'), col=c(1,4), pch=c(16, 1),
         cex=0.8)
})

with(E.col, { # Plot for colonization probability
  plot(1:2, Predicted, pch=1, xaxt='n', xlab='Year',
       ylab=expression(paste('Colonization probability ( ', gamma, ' )')),
       ylim=c(0,1), col=4)
  axis(1, at=1:2, labels=nd$year[1:2])
  arrows(1:2, lower, 1:2, upper, code=3, angle=90, length=0.03, col=4)
  points((1:2), gamma, col=1, lwd = 1, pch=16)
  legend(7, 1, c('Parameter', 'Estimate'), col=c(1,4), pch=c(16, 1),
         cex=0.8)
})
with(E.det, { # Plot for detection probability: note 3 years
  plot(1:3, Predicted, pch=1, xaxt='n', xlab='Year',
       ylab=expression(paste('Detection probability ( ', p, ' )')),
       ylim=c(0,1), col=4)
  axis(1, at=1:3, labels=nd$year)
  arrows(1:3, lower, 1:3, upper, code=3, angle=90, length=0.03, col=4)
  points((1:3), p, col=1, lwd = 1, pch=16)
  legend(7.5, 1, c('Parameter','Estimate'), col=c(1,4), pch=c(16, 1),
         cex=0.8)
})
par(op)


m1 <- nonparboot(m1, B = 10) # B should be 1000 for more complex models
cbind(smoothed=smoothed(m1)[2,], SE=m1@smoothed.mean.bsse[2,])

# # turnover (not sure this is relevant, just following along colext example)
# turnover <- function(fm) {
#   psi.hat <- plogis(coef(fm, type="psi"))
#   if(length(psi.hat) > 1)
#     stop("this function only works if psi is scalar")
#   T <- getData(fm)@numPrimary
#   tau.hat <- numeric(T-1)
#   gamma.hat <- plogis(coef(fm, type="col"))
#   phi.hat <- 1 - plogis(coef(fm, type="ext"))
#   if(length(gamma.hat) != T-1 | length(phi.hat) != T-1)
#     stop("this function only works if gamma and phi T-1 vectors")
#   for(t in 2:T) {
#     psi.hat[t] <- psi.hat[t-1]*phi.hat[t-1] +
#       (1-psi.hat[t-1])*gamma.hat[t-1]
#     tau.hat[t-1] <- gamma.hat[t-1]*(1-psi.hat[t-1]) / psi.hat[t]
#   }
#   return(tau.hat)
# }

# parametric bootstrap for turnover function
# again, not sure this is relevant
# pb <- parboot(m1, statistic=turnover, nsim=2)
# pb <- nonparboot(m1, statistic=turnover, nsim=2)
# turnCI <- cbind(pb@t0,
#                   t(apply(pb@t.star, 2, quantile, probs=c(0.025, 0.975))))
# colnames(turnCI) <- c("tau", "lower", "upper")
# turnCI

###--- Goodness of fit - experimental with dynamic occupancy models
# doesn't work with missing values....best to simulate?
# chisq <- function(fm) {
#   umf <- getData(fm)
#   y <- getY(umf)
#   sr <- fm@sitesRemoved
#   if(length(sr)>0)
#     y <- y[-sr,,drop=FALSE]
#   fv <- fitted(fm, na.rm=TRUE)
#   y[is.na(fv)] <- NA
#   sum((y-fv)^2/(fv*(1-fv)))
# }
# pb.gof <- parboot(m0, statistic=chisq, nsim=100)

################################################################################
# find the best supported model within a covariate category and then add complexity
# 1: simply test the effect of year
# 2: test the effect of sampling effort
# 3: test the site covariates
# 4: move on to the more comple yearly site covariates

names(colext.annual.covS)
# now try with site covariates as well as year
umf <- unmarkedMultFrame(
  y = yy,
  siteCovs = colext.site.covS,
  yearlySiteCovs = list(year = year, 
                        road=RD_dnsty, edge=EDGE_dnsty, trap=TRP_dnsty, tree20=TREE20_prop, canopy=CANOPY_prop, harvest=HARVEST_prop),
  obsCovs = list(effort=eeS), # use the scaled values for eeS (same values as in N-mixture)
  numPrimary=T)
summary(umf)

# Model with all constant parameters
# A model with constant parameters
fm0 <- colext(~1, ~1, ~1, ~1, umf)

# Like fm0, but with year-dependent detection
fmy1 <- colext(~1, ~1, ~1, ~year, umf)
# Like fm0, but with year-dependent colonization and extinction
fmy2 <- colext(~1, ~year-1, ~year-1, ~1, umf)
# A fully time-dependent model
fmy3 <- colext(~1, ~year-1, ~year-1, ~year, umf)

Yearmodels <- fitList('psi(.)gam(.)eps(.)p(.)' = fm0,
                      'psi(.)gam(.)eps(.)p(Y)' = fmy1, # best model
                      'psi(.)gam(Y)eps(Y)p(.)' = fmy2,
                      'psi(.)gam(Y)eps(Y)p(Y)' = fmy3)

Yms <- modSel(Yearmodels)
Yms

# add effort in to the detection portion and test with and without year
fme1 <- colext(~1, ~1, ~1, ~effort, umf)
fme2 <- colext(~1, ~1, ~1, ~effort+I(effort^2), umf)
fmye1 <- colext(~1, ~1, ~1, ~year+effort, umf)
fmye2 <- colext(~1, ~1, ~1, ~year+effort+I(effort^2), umf)

Effortmodels <- fitList('psi(.)gam(.)eps(.)p(.)' = fm0,
                      'psi(.)gam(.)eps(.)p(Y)' = fmy1,
                      'psi(.)gam(.)eps(.)p(E)' = fme1,
                      'psi(.)gam(.)eps(.)p(YE)' = fmye1,
                      'psi(.)gam(.)eps(.)p(E2)' = fme2,
                      'psi(.)gam(.)eps(.)p(YE2)' = fmye2) # clearly the best model

Ems <- modSel(Effortmodels)
Ems 

summary(fme2)
plogis(coef(fmye2)) # checking to see same back transformation with predict and plogis = TRUE

################################################################################
# Now consider the covariates in their groupings
# Sampling bias and year have been taken into account
# Forest Structure: SBS_prop, TREE20_prop, CANOPY_prop
# Landscape: RLW_dist, I(RLW_dist^2), utm_x, utm_y
# Human Disturbance: HARVEST_prop, EDGE_density, I(EDGE_density^2), RD_density, I(RD_density^2), TRP_DNSTY

###--- Forest Structure (on occupancy)
# Like fmye2 with variations on forest structure
fmFS1 <- colext(~SBS_prop, ~1, ~1, ~year+effort+I(effort^2), umf, control=list(maxit=1000, trace=TRUE, REPORT=1))
fmFS2 <- colext(~TREE20_prop, ~1, ~1, ~year+effort+I(effort^2),umf,control=list(maxit=1000, trace=TRUE, REPORT=1))
fmFS3 <- colext(~CANOPY_prop, ~1, ~1, ~year+effort+I(effort^2),umf,control=list(maxit=1000, trace=TRUE, REPORT=1))
fmFS4 <- colext(~SBS_prop + TREE20_prop, ~1, ~1, ~year+effort+I(effort^2), umf,control=list(maxit=1000, trace=TRUE, REPORT=1))
fmFS5 <- colext(~TREE20_prop + CANOPY_prop, ~1, ~1, ~year+effort+I(effort^2),umf, control=list(maxit=1000, trace=TRUE, REPORT=1))
fmFS6 <- colext(~CANOPY_prop + SBS_prop, ~1, ~1, ~year+effort+I(effort^2),umf, control=list(maxit=1000, trace=TRUE, REPORT=1))
fmFS7 <- colext(~CANOPY_prop + SBS_prop + TREE20_prop, ~1, ~1, ~year+effort+I(effort^2),umf,control=list(maxit=1000, trace=TRUE, REPORT=1))

ForestStructuremodels <- fitList('psi(.)gam(.)eps(.)p(YE2)' = fmye2,
                      'psi(SBS)gam(.)eps(.)p(YE2)' = fmFS1,
                      'psi(TREE)gam(.)eps(.)p(YE2)' = fmFS2,
                      'psi(CANOPY)gam(.)eps(.)p(YE2)' = fmFS3,
                      'psi(SBS_TREE)gam(.)eps(.)p(YE2)' = fmFS4,
                      'psi(TREE_CANOPY)gam(.)eps(.)p(YE2)' = fmFS5,
                      'psi(CANOPY_SBS)gam(.)eps(.)p(YE2)' = fmFS6,
                      'psi(CANOPY_SBS_TREE)gam(.)eps(.)p(YE2)' = fmFS7)

FSms <- modSel(ForestStructuremodels)
FSms # no real difference with or without SBS (came in as 0.26 difference from "top" model of fmye2)

summary(fmFS1)

###--- Landscape (on occupancy): RLW_dist, I(RLW_dist^2), utm_x, utm_y
# Like fmye2 & fmFS1 with variations on forest structure
fmLS1 <- colext(~SBS_prop + RLW_dist, ~1, ~1, ~year+effort+I(effort^2), umf, control=list(maxit=1000, trace=TRUE, REPORT=1))
fmLS2 <- colext(~SBS_prop + RLW_dist + I(RLW_dist^2), ~1, ~1, ~year+effort+I(effort^2), umf, control=list(maxit=1000, trace=TRUE, REPORT=1))
fmLS3 <- colext(~SBS_prop + utm_x, ~1, ~1, ~year+effort+I(effort^2), umf, control=list(maxit=1000, trace=TRUE, REPORT=1))
fmLS4 <- colext(~SBS_prop + utm_y, ~1, ~1, ~year+effort+I(effort^2), umf, control=list(maxit=1000, trace=TRUE, REPORT=1))
fmLS5 <- colext(~SBS_prop + RLW_dist + I(RLW_dist^2)+utm_x, ~1, ~1, ~year+effort+I(effort^2), umf, control=list(maxit=1000, trace=TRUE, REPORT=1))
fmLS6 <- colext(~SBS_prop + RLW_dist + I(RLW_dist^2)+utm_y, ~1, ~1, ~year+effort+I(effort^2), umf, control=list(maxit=1000, trace=TRUE, REPORT=1))
fmLS7 <- colext(~SBS_prop +utm_x +utm_y, ~1, ~1, ~year+effort+I(effort^2), umf, control=list(maxit=1000, trace=TRUE, REPORT=1))
fmLS8 <- colext(~SBS_prop + RLW_dist + I(RLW_dist^2)+utm_x+utm_y, ~1, ~1, ~year+effort+I(effort^2), umf, control=list(maxit=1000, trace=TRUE, REPORT=1))
fmLS9 <- colext(~RLW_dist + I(RLW_dist^2)+utm_x+utm_y, ~1, ~1, ~year+effort+I(effort^2), umf, control=list(maxit=1000, trace=TRUE, REPORT=1))
fmLS10 <- colext(~RLW_dist + I(RLW_dist^2)+utm_y, ~1, ~1, ~year+effort+I(effort^2), umf, control=list(maxit=1000, trace=TRUE, REPORT=1))

Landscapemodels <- fitList('psi(.)gam(.)eps(.)p(YE2)' = fmye2,
                           'psi(SBS_W)gam(.)eps(.)p(YE2)' = fmFS1,
                           'psi(SBS_W2)gam(.)eps(.)p(YE2)' = fmLS1,
                           'psi(SBS_utmx)gam(.)eps(.)p(YE2)' = fmLS3,
                           'psi(SBS_utmy)gam(.)eps(.)p(YE2)' = fmLS4,
                           'psi(SBS_W2_utmx)gam(.)eps(.)p(YE2)' = fmLS5,
                           'psi(SBS_W2_utmy)gam(.)eps(.)p(YE2)' = fmLS6,
                           'psi(SBS_utmx_utmy)gam(.)eps(.)p(YE2)' = fmLS7,
                           'psi(SBS_W2_utmx_utmy)gam(.)eps(.)p(YE2)' = fmLS8,
                           # 'psi(W2_utmx_utmy)gam(.)eps(.)p(YE2)' = fmLS9,
                           'psi(W2_utmy)gam(.)eps(.)p(YE2)' = fmLS10)

LSms <- modSel(Landscapemodels)
LSms # slightly better without SBS, but within 2 delta - no clear winners other than including W2 and possibly utmy

summary(fmLS10)

###--- Human Disturbance (on occupancy): HARVEST_prop, EDGE_density, I(EDGE_density^2), RD_density, I(RD_density^2), TRP_DNSTY
# Like fmLS10 with variations on human disturbance
fmHD1 <- colext(~RLW_dist + I(RLW_dist^2)+utm_y + HARVEST_prop,
                ~1, ~1, ~year+effort+I(effort^2), umf, control=list(maxit=1000, trace=TRUE, REPORT=1))
fmHD2 <- colext(~RLW_dist + I(RLW_dist^2)+utm_y + EDGE_dnsty + I(EDGE_dnsty^2),
                ~1, ~1, ~year+effort+I(effort^2), umf, control=list(maxit=1000, trace=TRUE, REPORT=1))
fmHD3 <- colext(~RLW_dist + I(RLW_dist^2)+utm_y + RD_dnsty + I(RD_dnsty^2),
                ~1, ~1, ~year+effort+I(effort^2), umf, control=list(maxit=1000, trace=TRUE, REPORT=1))
fmHD4 <- colext(~RLW_dist + I(RLW_dist^2)+utm_y + TRP_dnsty,
                ~1, ~1, ~year+effort+I(effort^2), umf, control=list(maxit=1000, trace=TRUE, REPORT=1))
fmHD5 <- colext(~RLW_dist + I(RLW_dist^2)+utm_y + HARVEST_prop + EDGE_dnsty + I(EDGE_dnsty^2),
                ~1, ~1, ~year+effort+I(effort^2), umf, control=list(maxit=1000, trace=TRUE, REPORT=1))
fmHD6 <- colext(~RLW_dist + I(RLW_dist^2)+utm_y + HARVEST_prop +RD_dnsty + I(RD_dnsty^2) ,
                ~1, ~1, ~year+effort+I(effort^2), umf, control=list(maxit=1000, trace=TRUE, REPORT=1))
fmHD7 <- colext(~RLW_dist + I(RLW_dist^2)+utm_y + HARVEST_prop + TRP_dnsty,
                ~1, ~1, ~year+effort+I(effort^2), umf, control=list(maxit=1000, trace=TRUE, REPORT=1))
fmHD8 <- colext(~RLW_dist + I(RLW_dist^2)+utm_y + EDGE_dnsty + I(EDGE_dnsty^2) + RD_dnsty + I(RD_dnsty^2) ,
                ~1, ~1, ~year+effort+I(effort^2), umf, control=list(maxit=1000, trace=TRUE, REPORT=1))
fmHD9 <- colext(~RLW_dist + I(RLW_dist^2)+utm_y + EDGE_dnsty + I(EDGE_dnsty^2) + TRP_dnsty,
                ~1, ~1, ~year+effort+I(effort^2), umf, control=list(maxit=1000, trace=TRUE, REPORT=1))
fmHD10 <- colext(~RLW_dist + I(RLW_dist^2)+utm_y + RD_dnsty + I(RD_dnsty^2) + TRP_dnsty,
                ~1, ~1, ~year+effort+I(effort^2), umf, control=list(maxit=1000, trace=TRUE, REPORT=1))
fmHD11 <- colext(~RLW_dist + I(RLW_dist^2)+utm_y + HARVEST_prop + EDGE_dnsty + I(EDGE_dnsty^2) + RD_dnsty + I(RD_dnsty^2) + TRP_dnsty,
                 ~1, ~1, ~year+effort+I(effort^2), umf, control=list(maxit=1000, trace=TRUE, REPORT=1))

Disturbancemodels <- fitList('psi(.)gam(.)eps(.)p(YE2)' = fmye2,
                              'psi(W2_utmy)gam(.)eps(.)p(YE2)' = fmLS10,
                              'psi(W2_utmy_H)gam(.)eps(.)p(YE2)' = fmHD1,
                              'psi(W2_utmy_E2)gam(.)eps(.)p(YE2)' = fmHD2, # top by 3.33
                              'psi(W2_utmy_R2)gam(.)eps(.)p(YE2)' = fmHD3,
                              'psi(W2_utmy_T)gam(.)eps(.)p(YE2)' = fmHD4,
                              'psi(W2_utmy_HE2)gam(.)eps(.)p(YE2)' = fmHD5,
                              'psi(W2_utmy_HR2)gam(.)eps(.)p(YE2)' = fmHD6,
                              'psi(W2_utmy_HT)gam(.)eps(.)p(YE2)' = fmHD7,
                              'psi(W2_utmy_E2R2)gam(.)eps(.)p(YE2)' = fmHD8,
                              'psi(W2_utmy_E2T)gam(.)eps(.)p(YE2)' = fmHD9,
                              'psi(W2_utmy_R2T)gam(.)eps(.)p(YE2)' = fmHD10,
                              'psi(W2_utmy_HE2R2T)gam(.)eps(.)p(YE2)' = fmHD11)

HDms <- modSel(Disturbancemodels)
HDms # slightly better without SBS, but within 2 delta - no clear winners other than including W2 and possibly utmy


fmHD2 <- nonparboot(fmHD2, B = 1000) # B should be 1000 for more complex models
save(fmHD2, file = paste0("./out/colext_top_model_fmHD2.RData"))

cbind(smoothed=smoothed(fmHD2)[2,], SE=fmHD2@smoothed.mean.bsse[2,])


TopModels <- fitList('psi(.)gam(.)eps(.)p(Y)' = fmy1,
                     'psi(.)gam(.)eps(.)p(YE2)' = fmye2,
                     'psi(SBS)gam(.)eps(.)p(YE2)' = fmFS1,
                     'psi(W2_utmy)gam(.)eps(.)p(YE2)' = fmLS10,
                     'psi(W2_utmy_E2)gam(.)eps(.)p(YE2)' = fmHD2)

top_models_ms <- modSel(TopModels)
top_models_ms
coef(top_models_ms)
SE(top_models_ms)
toExport <- as(top_models_ms, "data.frame")
write.csv(toExport, "out/colext_top_model_ms.csv")

# can also use backTransform function to back transform estimate, se and confidence intervals
toExport2 <- plogis(coef(fmHD2))
print(toExport2)
write.csv(toExport2, "out/colext_top_model_backtransform.csv")

summary(fmHD2)
# only things significant are year and effort on detection probability
# no need to graph the occupancy covariates, instead go with figures for detection probability

# op <- par(mfrow=c(2,1), mai=c(0.8,0.8,0.1,0.1))
# nd <- data.frame(utm_y=seq(-2,2, length=50),
#                  RLW_dist=seq(-1,5.5, length=50),
#                  EDGE_dnsty=seq(-1.5,3.5, length=50))
# # do this for W2, EDGE and utmY on occupancy
# E.psi <- predict(fmLS10, type="psi", newdata=nd, appendData=TRUE)
# E.psi$utmyOrig <- E.psi$utm_y*(sd(colext.site.cov$utm_y, na.rm=T)) + mean(colext.site.cov$utm_y, na.rm = T)
# E.psi$RLWOrig <- E.psi$RLW_dist*(sd(colext.site.cov$RLW_dist, na.rm=T)) + mean(colext.site.cov$RLW_dist, na.rm = T)
# 
# 
# with(E.psi, {
#   plot(utmyOrig, Predicted, ylim=c(0,1), type="l",
#        xlab="South-North Gradient",
#        ylab=expression(hat(psi)), cex.lab=0.8, cex.axis=0.8)
#   lines(utmyOrig, Predicted+1.96*SE, col=gray(0.7))
#   lines(utmyOrig, Predicted-1.96*SE, col=gray(0.7))
# })
# with(E.psi, {
#   plot(RLWOrig, Predicted, ylim=c(0,1), type="l",
#        xlab="Distance to Watercourse",
#        ylab=expression(hat(psi)), cex.lab=0.8, cex.axis=0.8)
#   lines(RLWOrig, Predicted+1.96*SE, col=gray(0.7))
#   lines(RLWOrig, Predicted-1.96*SE, col=gray(0.7))
# })
# do this for effort and year on detection

# need to do this for year as well
# nd <- data.frame(utm_y=seq(-2,2, length=50),
#                  RLW_dist=seq(-1,5.5, length=50),
#                  EDGE_dnsty=seq(-1.5,3.5, length=50),
#                  year=factor("1996-1997", levels=c(unique(year))),
#                  effort=seq(-0.3,3.5,length=50))

Cairo(file="out/colext_top_model_pdetfig.PNG",type="png",width=3000,height=2600,pointsize=13,bg="white",dpi=300)

nd <- data.frame(utm_y=0,
                 RLW_dist=0,
                 EDGE_dnsty=0,
                 year=factor("1996-1997", levels=c(unique(year))),
                 effort=seq(-0.3,3.5,length=50))

E.p <- predict(fmHD2, type="det", newdata=nd, appendData=TRUE)
# back transform effort to determine best detection probability per effort
# looks like peak of detection between 20-40 trap days of effort (# traps and days open per grid cell)
op <- par(mfrow=c(2,1), mai=c(0.9,0.9,0.2,0.2))
E.p$effOrig <- E.p$effort*(2*sd(ee, na.rm=T)) + mean(ee, na.rm = T)
with(E.p, {
  plot(effOrig, Predicted, ylim=c(0,1), type="l",col=4,
       xlab="Trap Effort", ylab=expression( italic(p) ),
       cex.lab=0.8, cex.axis=0.8)
  lines(effOrig, Predicted+1.96*SE, col=gray(0.7))
  lines(effOrig, Predicted-1.96*SE, col=gray(0.7))
})

nd <- data.frame(year=c("1996-1997","1997-1998","2019-2020"),
                 effort=0)
E.p <- predict(fmHD2, type='det', newdata=nd,appendData=TRUE)
# plogis(coef(m1)) # checking to see same back transformation with predict and plogis = TRUE
with(E.p, { # Plot for detection probability: note 3 years
  plot(1:3, Predicted, pch=1, xaxt='n', xlab='Year',
       ylab=expression( italic(p)),
       cex.lab=0.8, cex.axis=0.8,
       ylim=c(0,1), col=4)
  axis(1, at=1:3, labels=nd$year)
  arrows(1:3, lower, 1:3, upper, code=3, angle=90, length=0.03, col=4)
})
par(op)

dev.off()
