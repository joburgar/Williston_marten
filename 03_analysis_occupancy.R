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

M <- length(ALLcov.df_grid_to_use) # 150; number of sites
J <- max(length(retro.data.out[[1]]$days21.to.use),  # number of sampling occassions
         length(retro.data.out[[2]]$days21.to.use),
         ncol(rec.data.out$rec_effort))
T <- 3 # number of years / sesasons


# data(crossbill)
# head(crossbill); nrow(crossbill)
# summary(crossbill)
head(ALLcov.df)

ALLcov.df_grid_to_use <- ALLcov.df %>% filter(grepl("Traps", Area_Year)) %>% select(grid_id)
ALLcov.df_grid_to_use <- unique(ALLcov.df_grid_to_use$grid_id)
length(ALLcov.df_grid_to_use)

# create site covariates
colext.df <- ALLcov.df %>% filter(grid_id %in% ALLcov.df_grid_to_use) %>% filter(grepl("Area",Area_Year))
colext.df$RLW_dist_sq <- colext.df$RLW_dist*colext.df$RLW_dist
colext.df$RD_density_sq <- colext.df$RD_density*colext.df$RD_density
colext.df$EDGE_density_sq <- colext.df$EDGE_density*colext.df$EDGE_density
colext.site <- colext.df %>% select(-grid_id, -Year, -Area_Year, -Area)
head(colext.site)

# scale the site covariates
colext.site <- as.data.frame(scale(colext.site))
colext.site[is.na(colext.site)] <- 0 # do this to make the NAs for trap density the same as the mean value (could also just delete but then deletes from model)
summary(colext.site)
colext.site <- as.matrix(colext.site)

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

# yearly covariates (think it's just one row for each site...)
year <- c("1996-1997","1997-1998","2019-2020")
year <- matrix(year, nrow(yy),3,byrow=TRUE)

simUMF <- unmarkedMultFrame(
  y = yy,
  yearlySiteCovs = list(year = year),
  numPrimary=T)
summary(simUMF)

# Model with all constant parameters
m0 <- colext(psiformula= ~1, gammaformula = ~ 1, epsilonformula = ~ 1,
                 pformula = ~ 1, data = simUMF, method="BFGS")
summary(m0)

names(m0)
# [1] "psi" "col" "ext" "det"
# occupancy, colonization, extinction, detection

# backTransform(m0, type="psi")
# backTransform(m0, type="col")
# backTransform(m0, type="ext")
# backTransform(m0, type="det")

confint(backTransform(m0, type="psi"))


m1 <- colext(psiformula = ~1, # First-year occupancy
               gammaformula = ~ year-1, # Colonization
               epsilonformula = ~ year-1, # Extinction
               pformula = ~ year-1, # Detection
               data = simUMF)
m1


nd <- data.frame(year=c("1996-1997","1997-1998"))
E.ext <- predict(m1, type='ext', newdata=nd)
E.col <- predict(m1, type='col', newdata=nd)
nd <- data.frame(year=c("1996-1997","1997-1998","2019-2020"))
E.det <- predict(m1, type='det', newdata=nd)

op <- par(mfrow=c(3,1), mai=c(0.6, 0.6, 0.1, 0.1))
with(E.ext, { # Plot for extinction probability
  plot(1:2, Predicted, pch=1, xaxt='n', xlab='Year',
       ylab=expression(paste('Extinction probability ( ', epsilon, ' )')),
       ylim=c(0,1), col=4)
  axis(1, at=1:2, labels=nd$year[1:2])
  arrows(1:2, lower, 1:2, upper, code=3, angle=90, length=0.03, col=4)
  points((1:2)-0.1, 1-phi, col=1, lwd = 1, pch=16)
  legend(7, 1, c('Parameter', 'Estimate'), col=c(1,4), pch=c(16, 1),
         cex=0.8)
})

with(E.col, { # Plot for colonization probability
  plot(1:2, Predicted, pch=1, xaxt='n', xlab='Year',
       ylab=expression(paste('Colonization probability ( ', gamma, ' )')),
       ylim=c(0,1), col=4)
  axis(1, at=1:2, labels=nd$year[1:2])
  arrows(1:2, lower, 1:2, upper, code=3, angle=90, length=0.03, col=4)
  points((1:2)-0.1, gamma, col=1, lwd = 1, pch=16)
  legend(7, 1, c('Parameter', 'Estimate'), col=c(1,4), pch=c(16, 1),
         cex=0.8)
})
with(E.det, { # Plot for detection probability: note 3 years
  plot(1:3, Predicted, pch=1, xaxt='n', xlab='Year',
       ylab=expression(paste('Detection probability ( ', p, ' )')),
       ylim=c(0,1), col=4)
  axis(1, at=1:3, labels=nd$year)
  arrows(1:3, lower, 1:3, upper, code=3, angle=90, length=0.03, col=4)
  points((1:3)-0.1, p, col=1, lwd = 1, pch=16)
  legend(7.5, 1, c('Parameter','Estimate'), col=c(1,4), pch=c(16, 1),
         cex=0.8)
})
par(op)


# detection_history <- retro.data.out[[1]]$y_21day
# detection_history <- retro.data.out[[2]]$y_21day
detection_history <- rec.data.out$rec_ydata
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

# To get real estimate of detection (with 95% CI)
predict(occu.m1, 
        newdata = data.frame(site = 1),
        type = "det")

boot::inv.logit(coef(occu.m1)[1]) # Real estimate of occupancy
boot::inv.logit(coef(occu.m1)[2]) # Real estimate of detection

# Load covariate data
head(ALLcov.df)
grid_centroid_utm <- as.data.frame(st_coordinates(st_centroid(aoi_grid)))

ALLcov.df$utm_x <- grid_centroid_utm$X[match(ALLcov.df$grid_id, rownames(grid_centroid_utm))]
ALLcov.df$utm_y <- grid_centroid_utm$Y[match(ALLcov.df$grid_id, rownames(grid_centroid_utm))]

# cov.df <- ALLcov.df %>% filter(Area_Year=="Traps 1996") %>% select(-grid_id, -Area, -Year, -Area_Year)
# cov.df <- ALLcov.df %>% filter(Area_Year=="Traps 1997") %>% select(-grid_id, -Area, -Year, -Area_Year)
cov.df <- ALLcov.df %>% filter(Area_Year=="Traps 2019") %>% select(-grid_id, -Area, -Year, -Area_Year)
summary(cov.df)

cov.df$RLW_dist_sq <- cov.df$RLW_dist*cov.df$RLW_dist
cov.df$RD_density_sq <- cov.df$RD_density*cov.df$RD_density
cov.df$EDGE_density_sq <- cov.df$EDGE_density*cov.df$EDGE_density
dim(cov.df)
summary(cov.df)

cov.df_scaled <- as.data.frame(scale(cov.df))
summary(cov.df_scaled)
cov.df_scaled[is.na(cov.df_scaled)] <- 0

# scale(x,center=min(x),scale=diff(range(x)))
# https://stackoverflow.com/questions/5468280/scale-a-series-between-two-points

# effort <- retro.data.out[[1]]$effort.21days
# effort <- retro.data.out[[2]]$effort.21days
effort <- rec.data.out$rec_effort

# Build a new unmarkedFramOccu
# sample.unmarkedFrame_cov <- unmarkedFrameOccu( # y is a matrix with observed detection history 
#   # (0's and 1's, one row per site, one column per survey)
#   y = as.matrix(detection_history),
#   # obsCovs = observation covariates in a list, 
#   # each variable has site rows x survey columns
#   obsCovs = list(effort = effort),
#   # siteCovs = dataframe with site rows x column variables
#   siteCovs = cov.df)

sample.unmarkedFrame_cov_scaled <- unmarkedFrameOccu( # y is a matrix with observed detection history 
  # (0's and 1's, one row per site, one column per survey)
  y = as.matrix(detection_history),
  # obsCovs = observation covariates in a list, 
  # each variable has site rows x survey columns
  obsCovs = list(effort = effort),
  # siteCovs = dataframe with site rows x column variables
  siteCovs = cov.df_scaled)
# S4 class for occupancy model data
summary(sample.unmarkedFrame_cov_scaled)

occu.m2 <- occu(formula = ~effort # detection formula first
                ~1,
                data = sample.unmarkedFrame_cov_scaled)


occu.full <- occu(formula = ~effort # detection formula first
                ~ TRP_DNSTY + utm_x + utm_y + SBS_prop + RLW_dist + RLW_dist_sq + RD_density + RD_density_sq + TREE20_prop + 
                  CANOPY_prop + EDGE_density + EDGE_density_sq + HARVEST_prop, # occupancy formula second,
                data = sample.unmarkedFrame_cov_scaled)


# Summarize
summary(occu.m1)
summary(occu.m2)
summary(occu.full)

# dredge all possible combinations of the occupancy covariates
occ_dredge <- dredge(occu.full)

# model comparison to explore the results for occupancy
mc <- as.data.frame(occ_dredge) %>% 
  select(starts_with("psi"), df, AICc, delta, weight)
# # shorten names for printing
names(mc) <- names(mc) %>%
  str_remove("psi") %>% 
  coalesce(names(mc))

# take a quick peak at the model selection table
mst <- mutate_all(mc, ~ round(., 3)) %>% 
  head(20) %>% 
  knitr::kable()

# for the 2019 data, occu.m3 is the best model (next top model is > delta 2 away; still not great for aic weight)
# occu.m3 <- occu(formula = ~effort # detection formula first
#                 ~utm_x + SBS_prop + RLW_dist + RLW_dist_sq + RD_density_sq +
#                   EDGE_density + HARVEST_prop, # occupancy formula second,
#                 data = sample.unmarkedFrame_cov_scaled)
# summary(occu.m3)
# 
# occ_gof <- mb.gof.test(occu.m3, nsim = 1000, plot.hist = FALSE) # up to nsim=1000 when actually checking final model
# # hide the chisq table to give simpler output
# occ_gof$chisq.table <- NULL
# occ_gof
# 
# boot::inv.logit(coef(occu.m3)[1]) #occupancy
# boot::inv.logit(coef(occu.m3)[9]) #detection (per sampling session)

# for the 1997 data, use model averaging to get info for the top models within delta aicc
# select models with the most support for model averaging (< 2 delta aicc)
occ_dredge_delta <- get.models(occ_dredge, subset = delta <= 2)

# average models based on model weights 
occ_avg <- model.avg(occ_dredge_delta, fit = TRUE)
coef(occ_avg)


# recent_occ_output <- list(occu.m2=occu.m2, occu.m3=occu.m3, mc=mc,occ_dredge=occ_dredge,occ_gof=occ_gof)
# save(recent_occ_output, file = paste0("./out/recent_occ_output.RData"))
# # load("out/recent_occ_output.RData")
# recent_occ_output$mc
# write.csv(recent_occ_output$mc,"out/recent_occ_mc.csv")
# write.csv(as.data.frame(summary(recent_occ_output$occu.m3)),"out/recent_occu.m3.csv")

retro_occ_output <- list(occu.m2=occu.m2, mc=mc,occ_dredge=occ_dredge,occ_avg=occ_avg)
save(retro_occ_output, file = paste0("./out/retro97_occ_output.RData"))
load("out/retro97_occ_output.RData")
write.csv(retro_occ_output$mc,"out/retro_occ_mc.csv")
write.csv(as.data.frame(coef(retro_occ_output$occ_avg)),"out/retro_occ_avg.csv")


# 1996 models did not converge - only the null model converged with psi(Int) 0.9161621 and p(Int) 0.3087824 (but not a good fit)

# boot::inv.logit(coef(recent_occ_output$occu.m2)) # Real estimate of occupancy / detection
# boot::inv.logit(coef(retro_occ_output$occu.m2)) # Real estimate of occupancy / detection

boot::inv.logit(coef(recent_occ_output$occu.m3)[1]);boot::inv.logit(coef(recent_occ_output$occu.m3)[9])
boot::inv.logit(coef(retro_occ_output$occ_avg)[1]);boot::inv.logit(coef(retro_occ_output$occ_avg)[5])

