# Plant community change in Swiss mountain hay meadows 2003-2017
We infered different drivers of community change (i.e. Nitrogen deposition, climate warming, land-use change) in Swiss mountain hay meadows. The data were obtained from the Swiss biodiversity monitoring ([BDM](http://www.biodiversitymonitoring.ch/en/home.html)). 

This repository contains the data and R-code for the analyses. The `Manuscript.Rmd` file is an [Rmarkdown](https://rmarkdown.rstudio.com) that can be used to fully reporduce the manuscript. Wenn rendered it is transformed into the formated `Manuscript.pdf`. 

Likewise The `Appendix_A.Rmd` file is an [Rmarkdown](https://rmarkdown.rstudio.com) that can be used to fully reproduce the comparision of colonizing or disappearing species with random communities. Wenn rendered it is doing the simulations, stores the result of the simulation in the folder `Communitysim` and produces the formated `Appendx.pdf`. 

The repository contains also the folder `R` that contain additional R-scripts and the folder `RData` with the data stored as `.RData` files:

## Additional R-scripts
The folder `R` contains the following R-scripts:

- `Data_processing.R`: This file compiles the data from the BDM SQLite database. Note, that to run this script one needs acces to this data base. Furthermore it loads an excel file that contains the coordinates of the BDM-sites that are in mountain hay meadows. This file is protected since the locations of the BDM sites should not be published. The file saves the data from these locations as different .RData files. See below for a description of the .RData files.
- `Settings.R`: This file is sourced by `Manuscript.Rmd`. It loads the libraries, the RData and defines some formating settings for the figures.
- `Run_models.R`: This R scripts run the different models for the manuscript and stores the result in the folder `Modelfit`. We run the model in a separate script because some of them take several minuts to run.

## Data files
The following .RData files are produced by the `Data_processing.R` file and stored in the `RData` folder:

- `sites.RData`: It contains site specific data such as the values for the four gradients as well as the Nitrogen deposition estimated for the years 1980, 1990, 2000, 2010, 1015. It also contains the number of species that remained constant between the first and second survey (nochange1), that differed between first and second (change1), the number of constant species between second and third survey (nochange2) and the numer of species that differed (change2).
- `surv.RData`: This table contains one line per survey and columns for the survey year (yearPl), the year when the survey was planed (yearP), and the measures for community structure such as the number of species (SR), the mean Ellenberg values for temperature (T), humidity (F), nutrients (N) and Light (L).
- `pl.RData`: Contains the single plant observations with the survey year (yearPl), in which of the three surveys the observation was made (Visits), and the species Ellenberg values for temperature (T), humidity (F), nutrients (N) and light (L).
- `survdat.RData`: This table contain for all species that were observed during the first survey, if they were observed (then Occ = 1) or not observed (then Occ = 0) during the second survey. All these data are coded with Period = 1. The same was also done for all species that were observed during the second survey whether or not they were observed during the third survey. All these data are coded with Period = 2. Additionally the table contains the species Ellenberg values for temperature (T), humidity (F), nutrients (N) and light (L).
- `coldat.RData`: This table contain for all species that were not observed during the first survey, if they were observed (then Occ = 1) or not observed (then Occ = 0) during the second survey. All these data are coded with Period = 1. The same was also done for all species that were not observed during the second survey whether or not they were observed during the third survey. All these data are coded with Period = 2. Additionally the table contains the species Ellenberg values for temperature (T), humidity (F), nutrients (N) and light (L).







