#======================================================================================================
#
# Plant community change in mountain hay meadows
#
# Calc temporal trend for all species
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
# Calculate trends for all species
#------------------------------------------------------------------------------------------------------
# spec <- pl %>% group_by(aID_SP, T, F, N, L) %>%
#   dplyr::summarise(Nsites = n_distinct(aID_STAO),
#                    Nrecords = n())
# save(spec, file = "RData/spec.RData")
load("RData/spec.RData")

for(s in 1:nrow(spec)) {
  if(spec$Nsites[s] >= 5) {try({
    tt <- surv %>% left_join(pl %>% filter(aID_SP == spec$aID_SP[s]) %>% dplyr::select(aID_KD, Occ)) %>% 		
      replace_na(list(Occ = 0))
    mod <- stan_glmer(Occ ~ yr + (1|aID_STAO), data = tt, family = binomial, 
                      prior_intercept = normal(location = 0, scale = 10), prior = normal(location = 0, scale = 2.5), 
                      cores = ncores)
    trend <- as.data.frame(mod, pars = "yr")$yr
    save(trend, file = paste0("SpeciesTrends/", spec$aID_SP[s], ".RData"))
    spec[s, "Filename"] <- paste0("SpeciesTrends/", spec$aID_SP[s], ".RData")
  })}
}

save(spec, file = "RData/spec.RData")
