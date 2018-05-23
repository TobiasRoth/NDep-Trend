#======================================================================================================
#
# Plant community change in mountain hay meadows
#
# Settings for manuscript
#
#======================================================================================================

rm(list=ls(all=TRUE))

#------------------------------------------------------------------------------------------------------
# Libraries
#------------------------------------------------------------------------------------------------------
library(knitr)
library(kableExtra)
library(xtable)
library(tidyverse)
library(broom)
library(rstanarm)
library(simba)
library(Rmisc)

#------------------------------------------------------------------------------------------------------
# Load data
#------------------------------------------------------------------------------------------------------
load("RData/surv.RData")
load("RData/pl.RData")
load("RData/sites.RData")
load("RData/coldat.RData")
load("RData/survdat.RData")

#------------------------------------------------------------------------------------------------------
# Settings 
#------------------------------------------------------------------------------------------------------
# Knitr-Settings
opts_chunk$set(echo = FALSE, hide = TRUE, cache = TRUE, warning = FALSE, message = FALSE,
               fig.asp = .4, fig.width = 8, out.width = "100%")

# Stting for withing text juncs
inline_hook <- function(x) {
  if (is.numeric(x)) {
    x <- format(x, nsmall = 2, digits = 2)
  }
  x
}
knit_hooks$set(inline = inline_hook)

# ggplot settings
theme_set(theme_classic() + 
            theme(legend.position="bottom") + 
            theme(legend.title = element_blank()))
pd <- position_dodge(0.2) 
tcol <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", 
          "#D55E00", "#CC79A7")

# Settings for xtable()
options(xtable.comment = FALSE)
ndigits = 2

#------------------------------------------------------------------------------------------------------
# Functions to sample random communities and compare with species gained or lost 
#------------------------------------------------------------------------------------------------------
# Number of simulations
nsim <- 100

# Function to sample community and callculate CM
getsamples <- function(s, nspec, d, EL) {
  sample(d %>% pull(paste(EL)), nspec, replace = FALSE) %>% mean(na.rm = TRUE)
}

# Prepare data from one site and calculate mean and CrI for sample and mean for real community
getsim <- function(x, EL, gain) {
  tt <- pl %>% filter(aID_STAO == x) %>% 
    group_by(aID_SP, T, F, N, L) %>% 
    dplyr::summarise(Occ1 = sum(Visit == 1),
                     Occ2 = sum(Visit == 2),
                     Occ3 = sum(Visit == 3))
  if(gain)  ausw <- (tt$Occ1 == 0 & tt$Occ2 == 1) | (tt$Occ2 == 0 & tt$Occ3 == 1)
  if(!gain) ausw <- (tt$Occ1 == 1 & tt$Occ2 == 0) | (tt$Occ2 == 1 & tt$Occ3 == 0)
  tmp <- map_dbl(1:nsim, getsamples, nspec = sum(ausw), d = tt, EL = EL)
  list(
    nchange = sum(ausw),
    meanreal = mean(tt[ausw, paste(EL)] %>% unlist, na.rm = TRUE),
    mean = mean(tmp, na.rm = TRUE),
    lower = quantile(tmp, probs = 0.05, na.rm = TRUE),
    upper = quantile(tmp, probs = 0.95, na.rm = TRUE))
}

# Function to prepare the res file for plotting
makesim <- function(EL, gain, sitem) {
  res <- tibble(aID_STAO = sites$aID_STAO, Sitemeasure = pull(sites[, paste(sitem)]))
  tmp <- map(res$aID_STAO, getsim, EL = EL, gain = gain)
  res$meanreal <- tmp %>% map_dbl("meanreal")
  res$mean <- tmp %>% map_dbl("mean")
  res$lower <- tmp %>% map_dbl("lower")
  res$upper <- tmp %>% map_dbl("upper")
  res
}


