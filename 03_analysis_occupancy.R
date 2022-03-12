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

load("out/retro.data.out.RData")
# load("out/rec.data.out.RData")
cov.df <- read.csv("data/retro.covdata.1997.csv",row.names=1)
# cov.df <- read.csv("data/rec.covdata.2020.csv",row.names=1)


# ALLcov.df <- as.data.frame(rbind(read.csv("data/retro.covdata.1996.csv",row.names=1),
#                                  read.csv("data/retro.covdata.1997.csv",row.names=1),
#                                  read.csv("data/rec.covdata.2020.csv",row.names=1)))
# ALLcov.df$Year <- c(rep(1996,38),rep(1997,38), rep(2020,86))

cov.data.files <- c("data/aoi.covdata.1996.csv","data/retro.covdata.1996.csv",
                    "data/aoi.covdata.1997.csv","data/retro.covdata.1997.csv",
                    "data/aoi.covdata.2020.csv","data/rec.covdata.2020.csv")

ALLcov.df <-do.call(rbind,lapply(cov.data.files, read.csv))
ALLcov.df$X <- ALLcov.df$grid_id <- NULL
cov.df.nrows <- sapply( cov.data.files, function(f) nrow(read.csv(f)) )


ALLcov.df$Year <- c(rep(1996,cov.df.nrows[1]+cov.df.nrows[2]),
                    rep(1997,cov.df.nrows[3]+cov.df.nrows[4]), 
                    rep(2019,cov.df.nrows[5]+cov.df.nrows[6]))
ALLcov.df$Area <- c(rep("Area",cov.df.nrows[1]), rep("Traps",cov.df.nrows[2]),
                    rep("Area",cov.df.nrows[3]), rep("Traps",cov.df.nrows[4]),
                    rep("Area",cov.df.nrows[5]), rep("Traps",cov.df.nrows[6]))
ALLcov.df %>% group_by(Year) %>% count(Area)
ALLcov.df$Area_Year <- paste(ALLcov.df$Area, ALLcov.df$Year)

ALLcov.df_longer <- ALLcov.df %>% dplyr::select(-RLW_type,-Area,-Year) %>%
  pivot_longer(!(Area_Year), names_to = "Covariate", values_to = "Values")

names(ALLcov.df_longer)
unique(ALLcov.df_longer$Covariate)
ALLcov.df_longer$Covariate <- as.factor(recode(ALLcov.df_longer$Covariate, SBS_prop = "Proportion SBS", RLW_dist = "Distance to Water",
       RD_density = "Road Density", TREE20_prop = "Proportion Trees > 20 m", CANOPY_prop = "Proportion Canopy > 45%",
       EDGE_density = "Forest Edge Density", HARVEST_prop = "Proportion Harvested"))
levels(ALLcov.df_longer$Covariate)

ALLcov.df_longer$Covariate <- fct_relevel(ALLcov.df_longer$Covariate, "Proportion Canopy > 45%", "Proportion Harvested","Proportion SBS", "Proportion Trees > 20 m",
            "Distance to Water","Forest Edge Density", "Road Density")

ALLcov.df_longer$Area_Year <- fct_relevel(ALLcov.df_longer$Area_Year,  "Traps 1996","Area 1996", "Traps 1997","Area 1997","Traps 2019","Area 2019")

pal = pnw_palette(name="Cascades",n=6,type="discrete")

cov.plot <- ggplot(ALLcov.df_longer, aes(x=Covariate, y=Values, fill=as.factor(Area_Year))) +
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

# source("03_analysis_prep.R") # if not working concurrently

# 1999/00 detections high for first 3 weeks and then drop off - unlikely for models to converge
# different sampling effort in 1999/00 - concentrating only on fisher and moving traps to target recapturing fisher
# 1998/99 might also be worth ignoring as trapping effort became much more focused on fisher
#####################################################################################
# borrowed code from https://jamesepaterson.github.io/jamespatersonblog/2020-09-01_occupancyintroduction.html

# Load detection history (100 sites with 10 visits each)
detection_history <- retro.data.out[[2]]$y_21day
# detection_history <- rec.data.out$rec_ydata
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
grid_centroid_utm <- st_coordinates(st_centroid(retro.data.out[[2]]$grid_output$fishnet_grid_sf))
cov.df$utm_x <- grid_centroid_utm[,1]
cov.df$utm_y <- grid_centroid_utm[,2]

cov.df$RLW_type <- cov.df$grid_id <- NULL

# cov.df$utm_x <- rec.data.out$grid_centroid_utm[,1]
# cov.df$utm_y <- rec.data.out$grid_centroid_utm[,2]

summary(cov.df)

cov.df$RLW_dist_sq <- cov.df$RLW_dist*cov.df$RLW_dist
cov.df$RD_density_sq <- cov.df$RD_density*cov.df$RD_density
cov.df$EDGE_density_sq <- cov.df$EDGE_density*cov.df$EDGE_density
dim(cov.df)
summary(cov.df)

cov.df_scaled <- as.data.frame(scale(cov.df))
summary(cov.df_scaled)


# scale(x,center=min(x),scale=diff(range(x)))
# https://stackoverflow.com/questions/5468280/scale-a-series-between-two-points

effort <- retro.data.out[[2]]$effort.21days
# effort <- rec.data.out$rec_effort

# Build a new unmarkedFramOccu
sample.unmarkedFrame_cov <- unmarkedFrameOccu( # y is a matrix with observed detection history 
  # (0's and 1's, one row per site, one column per survey)
  y = as.matrix(detection_history),
  # obsCovs = observation covariates in a list, 
  # each variable has site rows x survey columns
  obsCovs = list(effort = effort),
  # siteCovs = dataframe with site rows x column variables
  siteCovs = cov.df)

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
                ~ utm_x + utm_y + SBS_prop + RLW_dist + RLW_dist_sq + RD_density + RD_density_sq + TREE20_prop + 
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

