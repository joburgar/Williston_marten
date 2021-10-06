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

#############################################################
# 00_SC_Informed Sigma.R
# adapated from Chandler & Royle 2013 (appendix)
# created by Joanna Burgar, 08-May-2018
# updated by Joanna Burgar, 04-Aug-2021 for Leroy Gonzales (streamlined)
# script for estimating calculating vaguely informed priors
#############################################################

#########################################################################
# Informative sigma priors (from Chandler & Royle)

###--- home range and call attentuation to sigma prior
# home range size in ha
# detection range in metres

# for martens in sub-boreal spruce biogeoclimatic zone
# Lofroth's 1993 MSc thesis, telemetry data suggested
# adult male marten home range mean size = 5.52 km2, range 2.95-8.53 km2
# adult female marten home range mean size = 4.55 km2m range 1.25-10.05 km2

sigma_prior <- function (hrsize = hrszie, detect_min = detect_min, detect_max = detect_max){
  hrlength <- sqrt(hrsize/pi*10000)

  hr_detect_min_area <- (hrlength+detect_min)^2*pi / 10000 # 65.62 ha, and
  hr_detect_max_area <- (hrlength+detect_max)^2*pi / 10000 # 135.61 ha

  # Following Royle et. al (2011), and assuming a
  # chi-squared distribution with 2 degrees of freedom,
  # the range of sigma is given by

  sigma_lower <- sqrt(hr_detect_min_area*10000/pi)/sqrt(5.99)   # 187 m
  sigma_upper <- sqrt(hr_detect_max_area*10000/pi)/sqrt(5.99)  # 268 m

  # Assuming a grid spacing of 1 unit = 100 m, aim to have a prior with most
  # of the density between:

  pd_lower <- sigma_lower/100
  pd_upper <- sigma_upper/100

  return(list(pd_lower, pd_upper))
}


###---
# now use the output to determine the range of dgamma by adjusting the last shape and scale
# one example, for a home range of 40 ha with a call detection range between 100-300 m
# hrsize <- 10
# detect_min <- 100
# detect_max <- 300

# two example, for a home range of 10 ha no detection range
# hrsize <- 10
# detect_min <- 0
# detect_max <- 0

sigma_prior(hrsize = 10, detect_min = 0, detect_max = 0)
# [1] 0.7289734
# so in this example, want the density of the simga prior at 0.7289734 so a peak ~0.7
# for a SC or SCR model, we want 2*sigma spacing of cameras or 0.7*2 = 1.46
# assuming a grid spacing of 100 m we want cameras spaced 1.46 units or ~ 146 m

qgamma(c(0.001,0.5,0.999),60,27) #  1.439910 2.209889 3.215138 - tight home ranges centred on 40
curve(dgamma(x,60,27), col='black',xlim=c(0,5), ylim=c(0,2))

