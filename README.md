# Plant community change in Swiss mountain hay meadows 2003-2017
We infered different drivers of community change (i.e. Nitrogen deposition, climate warming, land-use change) in Swiss mountain hay meadows. The data were obtained from the Swiss biodiversity monitoring ([BDM](http://www.biodiversitymonitoring.ch/en/home.html)). 

This repository contains the data and R-code for the analyses. The `Manuscript.Rmd` file is an [Rmarkdown](https://rmarkdown.rstudio.com) that can be used to fully reporduce our manuscript. Wenn rendered it is transformed into the formated `Manuscript.pdf`. 

Additionally the repository contains a folder that contain additional R-scripts and a folder with the data:

## Additional R-scripts
The folder `R` contains the following R-scripts:

- `Data_processing.R`: This file compiles the data from the BDM SQLite database. Note, that to run this script one needs acces to this data base. Furthermore it loads an excel file that contains the coordinates of the BDM-sites that are in mountain hay meadows. This file is protected since the locations of the BDM sites should not be published. The file saves the data from these locations as different .RData files. See below for a description of the .RData files.
- `Settings.R`: This file is sourced by `Manuscript.Rmd`. It loads the libraries, the RData and defines some formating settings for the figures.

## Data files
The following .RData files are produced by the `Data_processing.R` file and stored in the `RData` folder:

- `sites.RData`:
- `surv.RData`:
- `survdat.RData`:
- `pl.RData`:
- `coldat.RData`:


