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
library(arm)
library(nlme)
library(simba)

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

# Darstellung von R-Berechnungen innerhalb von Text
inline_hook <- function(x) {
  if (is.numeric(x)) {
    x <- format(x, nsmall = 2, digits = 2)
  }
  x
}
knit_hooks$set(inline = inline_hook)

# Plot settings
theme_set(theme_classic() + 
            theme(legend.position="bottom") + 
            theme(legend.title = element_blank()))
pd <- position_dodge(0.2) 
tcol <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", 
          "#D55E00", "#CC79A7")

# Remove comment for xtable
options(xtable.comment = FALSE)

# Number of posterior samples
nsim <- 100

#------------------------------------------------------------------------------------------------------
# Hilfsfunktionen
#------------------------------------------------------------------------------------------------------
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots == 1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


