# Plant community change in Swiss mountain hay meadows 2003-2017
This repository contains all the materials needed to reproduce the analyses in 

Roth T, Kohli L, Bühler C, Rihm B, Meuli RG, Meier R, Amrhein V. (accepted for publication). Species turnover reveals hidden effects of decreasing nitrogen deposition in mountain hay meadows. [PeerJ 7:e6347 DOI 10.7717/peerj.6347](https://peerj.com/preprints/27230/).

## Introduction
We inferred different drivers of community change (i.e. Nitrogen deposition, climate warming, land-use change) in Swiss mountain hay meadows. The data were obtained from the Swiss biodiversity monitoring ([BDM](http://www.biodiversitymonitoring.ch/en/home.html)). The results are currently published as preeprint in [PeerJ](https://peerj.com/preprints/27230/).

This repository contains the data and R-code for the analyses. The `Manuscript.Rmd` file is an [Rmarkdown](https://rmarkdown.rstudio.com) file that can be used to fully reproduce the manuscript. When rendered it is transformed into the formatted `Manuscript.pdf`.

Likewise The `Appendix_A.Rmd` is an [Rmarkdown](https://rmarkdown.rstudio.com) file that can be used to fully reproduce the comparison of colonizing or disappearing species with random communities. When rendered it is doing the simulations, stores the result of the simulation in the folder `Communitysim` and produces the formatted `Appendx.pdf`. 

The repository contains also the folder `Geodata` that contain the spatial data to produce a figure of the locations of the surveyed plots, `R` that contain additional R-scripts (see next section for more details), the folder `RData` with the data stored as `.RData` files (see below for more details) and `Settings` that contain a list of the R-packages (including the version number) that were used to produce the manuscript.

## Additional R-scripts
The folder `R` contains the following R-scripts:

- `Data_processing.R`: This file compiles the data from the BDM SQLite database. Note, that to run this script one needs access to this data base. Furthermore it loads an excel file that contains the coordinates of the BDM-sites that are in mountain hay meadows. This file is protected since the locations of the BDM sites should not be published. The file saves the data from these locations as different .RData files. See below for a description of the .RData files.
- `Settings.R`: This file is sourced by `Manuscript.Rmd`. It loads the libraries, the RData and defines some formatting settings for the figures.
- `Run_models.R`: This R scripts run the different models for the manuscript and stores the result in the folder `Modelfit`. We run the model in a separate script because some of them take several minutes to run.

## Data files
The following .RData files are produced by the `Data_processing.R` file and stored in the `RData` folder:

- `coldat.RData`: This table contain for all species that were not observed during the first survey, if they were observed (then Occ = 1) or not observed (then Occ = 0) during the second survey. All these data are coded with Period = 1. The same was also done for all species that were not observed during the second survey whether or not they were observed during the third survey. All these data are coded with Period = 2. Additionally the table contains the species Ellenberg values for temperature (T), humidity (F), nutrients (N) and light (L).
- `pl.RData`: Contains the single plant observations with the survey year (yearPl), in which of the three surveys the observation was made (Visits), and the species Ellenberg values for temperature (T), humidity (F), nutrients (N) and light (L).
- `sites.RData`: It contains site specific data such as the values for the four gradients as well as the Nitrogen deposition estimated for the years 1980, 1990, 2000, 2010, 1015. It also contains the number of species that remained constant between the first and second survey (nochange1), that differed between first and second (change1), the number of constant species between second and third survey (nochange2) and the numer of species that differed (change2).
- `surv.RData`: This table contains one line per survey and columns for the survey year (yearPl), the year when the survey was planed (yearP), and the measures for community structure such as the number of species (SR), the mean Ellenberg values for temperature (T), humidity (F), nutrients (N) and Light (L).
- `survdat.RData`: This table contain for all species that were observed during the first survey, if they were observed (then Occ = 1) or not observed (then Occ = 0) during the second survey. All these data are coded with Period = 1. The same was also done for all species that were observed during the second survey whether or not they were observed during the third survey. All these data are coded with Period = 2. Additionally the table contains the species Ellenberg values for temperature (T), humidity (F), nutrients (N) and light (L).
- `turnover_double.RData`: Data on pseudo-turnover of species of independent replicate surveys that were conducted by two different surveyors during the same year on the same site.
