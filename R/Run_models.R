#======================================================================================================
#
# Plant community change in mountain hay meadows
#
# Run the models using rstanarm
#
#======================================================================================================

rm(list=ls(all=TRUE))

#------------------------------------------------------------------------------------------------------
# Settings
#------------------------------------------------------------------------------------------------------
# Libraries
library(tidyverse)
library(rstanarm)

# Load data
load("RData/surv.RData")
load("RData/sites.RData")

# Settings for rstanarm
ncores = 4

#------------------------------------------------------------------------------------------------------
# Temporal trend of community measures
#------------------------------------------------------------------------------------------------------
# Temporal trend of species richness
mod <- stan_glmer(SR ~ yr + (yr|aID_STAO), data = surv, family = poisson, cores = ncores)
save(mod, file = "Modelfit/SR-trend.RData")

# Temporal trend of Ellenberg temperature value
mod <- stan_lmer(T ~ yr + (yr|aID_STAO), data = surv, cores = ncores)
save(mod, file = "Modelfit/T-trend.RData")

# Temporal trend of Ellenberg humidity value
mod <- stan_lmer(F ~ yr + (yr|aID_STAO), data = surv, cores = ncores)
save(mod, file = "Modelfit/F-trend.RData")

# Temporal trend of Ellenberg nutrient value
mod <- stan_lmer(N ~ yr + (yr|aID_STAO), data = surv, cores = ncores)
save(mod, file = "Modelfit/N-trend.RData")

# Temporal trend of Ellenberg light value
mod <- stan_lmer(L ~ yr + (yr|aID_STAO), data = surv, cores = ncores)
save(mod, file = "Modelfit/L-trend.RData")

#------------------------------------------------------------------------------------------------------
# Difference in tunover between first/second and second/third survey period
#------------------------------------------------------------------------------------------------------
tt <- rbind(
  sites %>% transmute(
    change = change1, 
    nochange = nochange1,
    period = 1),
  sites %>% transmute(
    change = change2, 
    nochange = nochange2,
    period = 2))
mod <- stan_glm(cbind(change, nochange) ~ period, data = tt, family = binomial, cores = ncores)
save(mod, file = "Modelfit/turnover-trend.RData")
