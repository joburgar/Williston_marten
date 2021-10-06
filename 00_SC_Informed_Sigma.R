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
# script for estimating / calculating vaguely informed priors
#############################################################

#########################################################################
# Informative sigma priors (from Chandler & Royle)

###--- home range and call attentuation to sigma prior
# home range size in ha
# detection range in metres

sigma_prior <- function (hrsize = hrszie, detect_min = detect_min, detect_max = detect_max){
  hrlength <- sqrt(hrsize/pi*10000)

  hr_detect_min_area <- (hrlength+detect_min)^2*pi / 10000
  hr_detect_max_area <- (hrlength+detect_max)^2*pi / 10000

  # Following Royle et. al (2011), and assuming a
  # chi-squared distribution with 2 degrees of freedom,
  # the range of sigma is given by

  sigma_lower <- sqrt(hr_detect_min_area*10000/pi)/sqrt(5.99)
  sigma_upper <- sqrt(hr_detect_max_area*10000/pi)/sqrt(5.99)

  # Assuming a grid spacing of 1 unit = 1000 m or 1 km, aim to have a prior with most
  # of the density between:

  pd_lower <- sigma_lower/1000
  pd_upper <- sigma_upper/1000

  return(list(pd_lower, pd_upper))
}


###---
# now use the output to determine the range of dgamma by adjusting the last shape and scale

# for martens in sub-boreal spruce biogeoclimatic zone
# Lofroth's 1993 MSc thesis, telemetry data suggested
# adult male marten home range mean size = 5.52 km2, range 2.95-8.53 km2
# adult female marten home range mean size = 4.55 km2m range 1.25-10.05 km2

# to find out the range go with 125 ha to 1005 ha, with the bulk centred on 500 ha
# hrsize <- 125 # lower
# hrsize <- 1005 # upper
# hrsize <- 500 # mid
# detect_min <- 0
# detect_max <- 0

sigma_prior(hrsize = 500, detect_min = 0, detect_max = 0)
# [1] 0.258
# so in this example, want the density of the simga prior to range between 0.258 and 0.731
# with a peak at 0.515
 
# for a SC or SCR model, we want 2*sigma spacing of cameras ~ 0.515*2 = 1.03 km (~0.5 - 1.5 km)
# no smaller than 0.258*2 = 0.516 km and no greater than 0.731*2 = 1.462 km
# assuming a grid spacing of 1 km we want traps spaced ~ 1 km apart

qgamma(c(0.001,0.5,0.999),20,20) #  1.439910 2.209889 3.215138 - tight home ranges centred on 40
curve(dgamma(x,20,20), col='black',xlim=c(0,3), ylim=c(0,2))

