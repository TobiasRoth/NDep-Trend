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

# Select surveys and add site data
surv <- tbl(db, "KD_Z9") %>% 
  as.tibble() %>% 
  filter(!is.na(match(aID_STAO, sel)) & !is.na(yearPl) & yearP >= 2003 & yearP <= 2017) %>% 
  transmute(aID_KD = aID_KD, aID_STAO = aID_STAO, yearP = yearP, yearPl = yearPl, Visit = as.integer(floor((yearP - 1998) / 5))) %>% 
  left_join(tbl(db, "RAUMDATEN_Z9") %>% 
              transmute(aID_STAO = aID_STAO, Elevation = Hoehe, Inclination = Neig, 
                       Temperature = bio1, Precipitation = bio12), copy = TRUE) %>% 
  left_join(tbl(db, "RAUMDATEN_Z9_NDEP") %>% 
              transmute(aID_STAO = aID_STAO, NTOT = NTOT2010_Wiesen), copy = TRUE)

# Add community measures to 'surv'
surv <- surv %>% left_join(
  tbl(db, "PL") %>% 
    filter(!Z7) %>% 
    left_join(tbl(db, "Traits_PL")) %>%
    group_by(aID_KD) %>% 
    dplyr::summarise(
      SR = n(),
      SR_oligo = sum(N <= 2, na.rm = TRUE),
      SR_eutro = sum(N >= 4, na.rm = TRUE),
      SLA = mean(log(SLA), na.rm = TRUE),
      CH = mean(log(CH), na.rm = TRUE),
      SM = mean(log(SM), na.rm = TRUE),
      SLA_sd = sd(log(SLA), na.rm = TRUE),
      CH_sd = sd(log(CH), na.rm = TRUE),
      SM_sd = sd(log(SM), na.rm = TRUE),
      T = mean(T, na.rm = TRUE), 
      F = mean(F, na.rm = TRUE),
      N = mean(N, na.rm = TRUE),
      L = mean(L, na.rm = TRUE),
      MV = mean(MV, na.rm = TRUE)), copy = TRUE)

#------------------------------------------------------------------------------------------------------
# Selection plants data
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
# Sites data 
#------------------------------------------------------------------------------------------------------
# Add Turnover to 'surv'
sites <- data.frame(aID_STAO = as.integer(sel)) %>% 
  as.tibble() %>% 
  left_join(tbl(db, "RAUMDATEN_Z9") %>% 
              transmute(aID_STAO = aID_STAO, Elevation = Hoehe, Inclination = Neig, 
                            Temperature = bio1, Precipitation = bio12), copy = TRUE) %>% 
  left_join(tbl(db, "RAUMDATEN_Z9_NDEP") %>% 
              transmute(aID_STAO = aID_STAO, NTOT = NTOT2010_Wiesen), copy = TRUE)

for(i in 1:nrow(sites)) {
  tt <- sim(pl %>% filter(aID_STAO == sites$aID_STAO[i] & Visit <= 2) %>% dplyr::select(aID_KD, aID_SP, Occ), method = "cocogaston", listin = TRUE, listout = TRUE)
  sites[i, "nochange1"] <- tt$a
  sites[i, "change1"] <- tt$b + tt$c
  tt <- sim(pl %>% filter(aID_STAO == sites$aID_STAO[i] & Visit >= 2) %>% dplyr::select(aID_KD, aID_SP, Occ), method = "cocogaston", listin = TRUE, listout = TRUE)
  sites[i, "nochange2"] <- tt$a
  sites[i, "change2"] <- tt$b + tt$c
}

#------------------------------------------------------------------------------------------------------
# Data transformation
#------------------------------------------------------------------------------------------------------
surv$yr <- (surv$yearPl - 2010) / 10
surv$Elevation <- (surv$Elevation - 1000) / 500
surv$Temperature <- (surv$Temperature - 50) / 10
surv$Precipitation <- (surv$Precipitation - 1000) / 200
surv$NTOT <- (surv$NTOT - 10) / 10
surv$Inclination <- (surv$Inclination - 10) / 10

sites$Elevation <- (sites$Elevation - 1000) / 500
sites$Temperature <- (sites$Temperature - 50) / 10
sites$Precipitation <- (sites$Precipitation - 1000) / 200
sites$NTOT <- (sites$NTOT - 10) / 10
sites$Inclination <- (sites$Inclination - 10) / 10

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
  replace_na(replace = list(Occ1 = 0, Occ2 = 0, Occ3 = 0)) %>% 
  left_join(sites) 

# Prepare colonization data
coldat <- rbind(
  d %>% filter(Occ1 == 0) %>% 
    transmute(aID_STAO = aID_STAO, aID_SP =aID_SP, 
              Temperature = Temperature, Precipitation = Precipitation,
              NTOT = NTOT, Inclination = Inclination,
              T = T, F = F, N = N, L = L, Occ = Occ2, Period = 1),
  d %>% filter(Occ2 == 0) %>% 
    transmute(aID_STAO = aID_STAO, aID_SP =aID_SP, 
              Temperature = Temperature, Precipitation = Precipitation,
              NTOT = NTOT, Inclination = Inclination,
              T = T, F = F, N = N, L = L, Occ = Occ3, Period = 2))%>% 
  filter(!is.na(T) & !is.na(F) & !is.na(N) & !is.na(L)) 
coldat <- coldat %>% left_join(coldat %>% group_by(aID_SP) %>% dplyr::summarise(ausw = sum(Occ)>0)) %>% 
  filter(ausw)

# Prepare survival dat
survdat <- rbind(
  d %>% filter(Occ1 == 1) %>% 
    transmute(aID_STAO = aID_STAO, aID_SP =aID_SP, 
              Temperature = Temperature, Precipitation = Precipitation,
              NTOT = NTOT, Inclination = Inclination,
              T = T, F = F, N = N, L = L, Occ = Occ2, Period = 1),
  d %>% filter(Occ2 == 1) %>% 
    transmute(aID_STAO = aID_STAO, aID_SP =aID_SP, 
              Temperature = Temperature, Precipitation = Precipitation,
              NTOT = NTOT, Inclination = Inclination,
              T = T, F = F, N = N, L = L, Occ = Occ3, Period = 2)) %>% 
  filter(!is.na(T) & !is.na(F) & !is.na(N) & !is.na(L))

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

