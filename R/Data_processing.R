#======================================================================================================
#
# Plant community change in mountain hay meadows
#
# DATA PREPARATION
#
#======================================================================================================

rm(list=ls(all=TRUE))

#------------------------------------------------------------------------------------------------------
# Settings
#------------------------------------------------------------------------------------------------------
# Libraries
library(tidyverse)
library(simba)

# Connection to data base
db <- src_sqlite(path = "~/Documents/Dropbox/DB_BDM.db", create = FALSE)

#------------------------------------------------------------------------------------------------------
# Selection of sites
#------------------------------------------------------------------------------------------------------
# Import list with sites
sel <- read_delim("Data_protected/1366 K_Standort.csv", delim = ";") %>% 
  filter(E23_1366 == "ja") %>% 
  transmute(aID_STAO = aIdStaoZ9)

# Select only sites where three surveys are available since 2003
tmp <- tbl(db, "KD_Z9") %>% 
  as.tibble() %>% 
  filter(!is.na(match(aID_STAO, sel$aID_STAO)) & !is.na(yearPl) & yearP >= 2003 & yearP <= 2017) %>% 
  dplyr::select(aID_KD, aID_STAO, yearPl) %>% 
  left_join(tbl(db, "PL") %>% group_by(aID_KD) %>% dplyr::summarise(SR = n()), copy = TRUE) %>% 
  filter(SR >= 0)
tmp <- table(tmp$aID_STAO)
sum(tmp != 3)
sel <- names(tmp[tmp == 3])

#------------------------------------------------------------------------------------------------------
# Prepare site data 
#------------------------------------------------------------------------------------------------------
# Add Turnover to 'surv'
sites <- tibble(aID_STAO = as.integer(sel)) %>% 
  left_join(tbl(db, "RAUMDATEN_Z9") %>% 
              transmute(aID_STAO = aID_STAO, Elevation = Hoehe, Inclination = Neig, Aspect = Expos, CACO3 = CACO3, 
                        Temperature = bio1, Precipitation = bio12), copy = TRUE) %>% 
  left_join(tbl(db, "RAUMDATEN_Z9_NABO") %>% 
              transmute(aID_STAO = aID_STAO, pH = pH), by = c("aID_STAO" = "aID_STAO"), copy = TRUE) %>% 
  left_join(tbl(db, "RAUMDATEN_Z9_NDEP") %>% 
              transmute(aID_STAO = aID_STAO, NTOT_1980 = NTOT1980_Wiesen, NTOT_1990 = NTOT1990_Wiesen, 
                        NTOT_2000 = NTOT2000_Wiesen, NTOT_2010 = NTOT2010_Wiesen, NTOT_2015 = NTOT2015_Wiesen), copy = TRUE)

# Standardize covariates
sites$Elevation <- (sites$Elevation - 1000) / 500
sites$Temperature <- (sites$Temperature - 50) / 10
sites$Precipitation <- (sites$Precipitation - 1000) / 200
sites$NTOT <- (sites$NTOT_2010 - 10) / 10
sites$Inclination <- (sites$Inclination - 10) / 10
sites$pH <- sites$pH - 6

#------------------------------------------------------------------------------------------------------
# Prepare survey data 
#------------------------------------------------------------------------------------------------------
# Select surveys
surv <- tbl(db, "KD_Z9") %>% 
  as.tibble() %>% 
  filter(!is.na(match(aID_STAO, sel)) & !is.na(yearPl) & yearP >= 2003 & yearP <= 2017) %>% 
  transmute(aID_KD = aID_KD, aID_STAO = aID_STAO, yearP = yearP, yearPl = yearPl, Visit = as.integer(floor((yearP - 1998) / 5)))

# Add community measures to 'surv'
surv <- surv %>% left_join(
  tbl(db, "PL") %>% 
    filter(!Z7) %>% 
    left_join(tbl(db, "Traits_PL")) %>%
    group_by(aID_KD) %>% 
    dplyr::summarise(
      SR = n(),
      T = mean(T, na.rm = TRUE), 
      F = mean(F, na.rm = TRUE),
      N = mean(N, na.rm = TRUE),
      L = mean(L, na.rm = TRUE)), copy = TRUE)

# Standardize years
surv$yr <- (surv$yearPl - 2010) / 10

#------------------------------------------------------------------------------------------------------
# Prepare plant data
#------------------------------------------------------------------------------------------------------
pl <- tbl(db, "PL") %>% 
  filter(Z7 == 0 & !is.na(aID_SP)) %>% 
  as.tibble() %>% 
  filter(!is.na(match(aID_KD, surv$aID_KD))) %>% 
  left_join(tbl(db, "KD_Z9"), copy = TRUE) %>% 
  transmute(aID_KD = aID_KD, aID_STAO = aID_STAO, aID_SP = aID_SP, yearP = yearP, yearPl = yearPl, 
            Visit = as.integer(floor((yearP - 1998) / 5)),
            Occ = as.integer(1)) %>% 
  left_join(tbl(db, "Traits_PL") %>% dplyr::select(aID_SP, T, F, N, L), copy = TRUE) 

#------------------------------------------------------------------------------------------------------
# Add data for turnover calculation to 'sites'
#------------------------------------------------------------------------------------------------------
for(i in 1:nrow(sites)) {
  tt <- sim(pl %>% filter(aID_STAO == sites$aID_STAO[i] & Visit <= 2) %>% dplyr::select(aID_KD, aID_SP, Occ), method = "cocogaston", listin = TRUE, listout = TRUE)
  sites[i, "nochange1"] <- tt$a
  sites[i, "change1"] <- tt$b + tt$c
  tt <- sim(pl %>% filter(aID_STAO == sites$aID_STAO[i] & Visit >= 2) %>% dplyr::select(aID_KD, aID_SP, Occ), method = "cocogaston", listin = TRUE, listout = TRUE)
  sites[i, "nochange2"] <- tt$a
  sites[i, "change2"] <- tt$b + tt$c
}

#------------------------------------------------------------------------------------------------------
# Turnover between two survey from the same year
#------------------------------------------------------------------------------------------------------
# Select sites with two surveys per year
turnover_double <- tbl(db, "KD_Z9") %>% 
  filter(Aufnahmetyp == "Doppelaufnahme_Z9_Pflanzen") %>% 
  filter(yearPl >= 2003 & yearPl <= 2017) %>% 
  as.tibble() %>% 
  filter(!is.na(match(aID_STAO, sites$aID_STAO))) %>% 
  transmute(aID_KD_K = aID_KD,  aID_STAO = aID_STAO, yearPl = yearPl) %>% 
  left_join(surv %>% transmute(aID_KD_R = aID_KD,  aID_STAO = aID_STAO, yearPl = yearP)) %>% 
  filter(!is.na(aID_KD_R))

getturnover <- function(x) {
  tbl(db, "PL") %>% 
    filter(aID_KD == turnover_double$aID_KD_K[x] | 
             aID_KD == turnover_double$aID_KD_R[x]) %>% 
    filter(!is.na(aID_SP)) %>% 
    transmute(aID_KD = aID_KD, aID_SP = aID_SP, Occ = 1) %>% 
    sim(method = "cocogaston", listin = TRUE, listout = TRUE) %>% 
    pull(cocogaston)
}

turnover_double$turnover <- map_dbl(1:nrow(turnover_double), getturnover)

#------------------------------------------------------------------------------------------------------
# Compile data for colonization and local survival calculations
#------------------------------------------------------------------------------------------------------
d <- expand.grid(aID_STAO = sites$aID_STAO, aID_SP = unique(pl$aID_SP)) %>% as.tibble() %>% 
  left_join(pl %>% group_by(aID_SP) %>% dplyr::summarise(T = min(T), F = min(F), N = min(N), L = min(L))) %>% 
  left_join(pl %>% filter(yearP >= 2003 & yearP<= 2007) %>% 
              transmute(aID_STAO = aID_STAO, aID_SP = aID_SP, Occ1 = Occ)) %>% 
  left_join(pl %>% filter(yearP >= 2008 & yearP <= 2012) %>% 
              transmute(aID_STAO = aID_STAO, aID_SP = aID_SP, Occ2 = Occ)) %>% 
  left_join(pl %>% filter(yearP >= 2013 & yearP <= 2017) %>% 
              transmute(aID_STAO = aID_STAO, aID_SP = aID_SP, Occ3 = Occ)) %>% 
  replace_na(replace = list(Occ1 = 0, Occ2 = 0, Occ3 = 0))

# Prepare colonization data
coldat <- rbind(
  d %>% filter(Occ1 == 0) %>% 
    transmute(aID_STAO = aID_STAO, aID_SP =aID_SP, 
              T = T, F = F, N = N, L = L, Occ = Occ2, Period = 1),
  d %>% filter(Occ2 == 0) %>% 
    transmute(aID_STAO = aID_STAO, aID_SP =aID_SP,
              T = T, F = F, N = N, L = L, Occ = Occ3, Period = 2))

# Remove species that never colonized a site
coldat <- coldat %>% left_join(coldat %>% group_by(aID_SP) %>% dplyr::summarise(ausw = sum(Occ)>0)) %>%
  filter(ausw)

# Prepare survival dat
survdat <- rbind(
  d %>% filter(Occ1 == 1) %>% 
    transmute(aID_STAO = aID_STAO, aID_SP =aID_SP, 
              T = T, F = F, N = N, L = L, Occ = Occ2, Period = 1),
  d %>% filter(Occ2 == 1) %>% 
    transmute(aID_STAO = aID_STAO, aID_SP =aID_SP,
              T = T, F = F, N = N, L = L, Occ = Occ3, Period = 2))

#------------------------------------------------------------------------------------------------------
# Change siteID and save data
#------------------------------------------------------------------------------------------------------
siteIDs <- data.frame(
  aID_STAO = sites$aID_STAO,
  siteID = as.integer(factor(sites$aID_STAO))
  )

surv$aID_STAO <- siteIDs$siteID[match(surv$aID_STAO, siteIDs$aID_STAO)]
save(surv, file = "RData/surv.RData")

pl$aID_STAO <- siteIDs$siteID[match(pl$aID_STAO, siteIDs$aID_STAO)]
save(pl, file = "RData/pl.RData")

sites$aID_STAO <- siteIDs$siteID[match(sites$aID_STAO, siteIDs$aID_STAO)]
save(sites, file = "RData/sites.RData")

coldat$aID_STAO <- siteIDs$siteID[match(coldat$aID_STAO, siteIDs$aID_STAO)]
save(coldat, file = "RData/coldat.RData")

survdat$aID_STAO <- siteIDs$siteID[match(survdat$aID_STAO, siteIDs$aID_STAO)]
save(survdat, file = "RData/survdat.RData")

turnover_double$aID_STAO <- siteIDs$siteID[match(turnover_double$aID_STAO, siteIDs$aID_STAO)]
save(turnover_double, file = "RData/turnover_double.RData")
