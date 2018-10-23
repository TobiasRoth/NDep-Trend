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
# Source overall settings
source("R/Settings.R")

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
    aID_STAO = aID_STAO,
    change = change1, 
    nochange = nochange1,
    Temperature = Temperature,
    Precipitation = Precipitation,
    NTOT = NTOT,
    Inclination = Inclination, 
    period = 1),
  sites %>% transmute(
    aID_STAO = aID_STAO,
    change = change2, 
    nochange = nochange2,
    Temperature = Temperature,
    Precipitation = Precipitation,
    NTOT = NTOT,
    Inclination = Inclination, 
    period = 2))

tt$SR <- ((tt$change + tt$nochange) - 50) / 10
mod <- stan_glmer(cbind(change, nochange) ~ period + SR + Temperature + Precipitation + NTOT + Inclination + (1|aID_STAO), 
                data = tt, family = binomial, cores = ncores)
save(mod, file = "Modelfit/turnover-trend.RData")

#------------------------------------------------------------------------------------------------------
# Model for the colonization and local survival probabilities
#------------------------------------------------------------------------------------------------------
coldat <- coldat %>% left_join(sites) %>% dplyr::filter(!is.na(N))
mod <- glmer(Occ ~ NTOT * N + (1|aID_SP) + (1|aID_STAO), data = coldat, family = binomial)
save(mod, file = "Modelfit/colonization-N.RData")

survdat <- survdat %>% left_join(sites) %>% dplyr::filter(!is.na(N))
mod <- glmer(Occ ~ NTOT * N + (1|aID_SP) + (1|aID_STAO), data = survdat, family = binomial)
save(mod, file = "Modelfit/localsurvival-N.RData")


