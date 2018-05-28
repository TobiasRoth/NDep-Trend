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
# library(xtable)
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




