#======================================================================================================
#
# Plant community change in mountain hay meadows
#
# Compare gains and losses
#
#======================================================================================================

rm(list=ls(all=TRUE))

#------------------------------------------------------------------------------------------------------
# Settings and load data
#------------------------------------------------------------------------------------------------------
# Libraries
library(tidyverse)
library(arm)

# Load data
load("RData/pl.RData")
load("RData/sites.RData")

#------------------------------------------------------------------------------------------------------
# Compile absence data for all species x site combinations
#------------------------------------------------------------------------------------------------------
d <- expand.grid(aID_STAO = sites$aID_STAO, aID_SP = unique(pl$aID_SP)) %>% as.tibble() %>% 
  left_join(pl %>% group_by(aID_SP) %>% summarise(T = min(T), F = min(F), N = min(N), L = min(L))) %>% 
  left_join(pl %>% filter(yearP >= 2003 & yearP<= 2007) %>% 
              transmute(aID_STAO = aID_STAO, aID_SP = aID_SP, Occ1 = Occ)) %>% 
  left_join(pl %>% filter(yearP >= 2008 & yearP <= 2012) %>% 
              transmute(aID_STAO = aID_STAO, aID_SP = aID_SP, Occ2 = Occ)) %>% 
  left_join(pl %>% filter(yearP >= 2013 & yearP <= 2017) %>% 
              transmute(aID_STAO = aID_STAO, aID_SP = aID_SP, Occ3 = Occ)) %>% 
  replace_na(replace = list(Occ1 = 0, Occ2 = 0, Occ3 = 0)) %>% 
  left_join(sites) 

#------------------------------------------------------------------------------------------------------
# Analyse colonization
#------------------------------------------------------------------------------------------------------
summary(glmer(Occ2 ~ T + F + N + L + (1|aID_SP) + (1|aID_STAO), data = d %>% filter(Occ1 == 0), family = binomial))
summary(glmer(Occ3 ~ T + F + N + L + (1|aID_SP) + (1|aID_STAO), data = d %>% filter(Occ2 == 0), family = binomial))

#------------------------------------------------------------------------------------------------------
# Analyse local survival
#------------------------------------------------------------------------------------------------------
summary(glmer(Occ2 ~ T + F + N + L + (1|aID_SP) + (1|aID_STAO), data = d %>% filter(Occ1 == 1), family = binomial))
summary(glmer(Occ3 ~ T + F + N + L + (1|aID_SP) + (1|aID_STAO), data = d %>% filter(Occ2 == 1), family = binomial))

summary(glmer(Occ2 ~ T + F + NTOT2007 * N + L + (1|aID_SP) + (1|aID_STAO), data = d %>% filter(Occ1 == 1), family = binomial))
summary(glmer(Occ3 ~ T + F + NTOT2007 * N + L + (1|aID_SP) + (1|aID_STAO), data = d %>% filter(Occ2 == 1), family = binomial))








