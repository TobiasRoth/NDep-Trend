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
sel <- read_delim("Data/1366 K_Standort.csv", delim = ";") %>% 
  filter(E23_1366 == "ja") %>% 
  transmute(aID_STAO = aIdStaoZ9)

# Select only sites where three surveys are available since 2003
tmp <- tbl(db, "KD_Z9") %>% 
  as.tibble() %>% 
  filter(!is.na(match(aID_STAO, sel$aID_STAO)) & !is.na(yearPl) & yearP >= 2003 & yearP <= 2017) %>% 
  dplyr::select(aID_KD, aID_STAO, yearPl) %>% 
  left_join(tbl(db, "PL") %>% group_by(aID_KD) %>% summarise(SR = n()), copy = TRUE) %>% 
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
              dplyr::select(aID_STAO, Hoehe, Neig, Expos, AS09_72, AS09_17, CACO3, Vernass, BGR_6, 
                            bio1, bio12, bio13, bio14), copy = TRUE) %>% 
  left_join(tbl(db, "RAUMDATEN_Z9_NDEP") %>% 
              dplyr::select(aID_STAO, NTOT2007), copy = TRUE)
  
# Add community measures to 'surv'
surv <- surv %>% left_join(
  tbl(db, "PL") %>% 
    filter(!Z7) %>% 
    left_join(tbl(db, "Traits_PL")) %>%
    group_by(aID_KD) %>% 
    summarise(SR = n(),
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
              L = mean(L, na.rm = TRUE)), copy = TRUE)

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
              dplyr::select(aID_STAO, Hoehe, Neig, Expos, AS09_72, AS09_17, CACO3, Vernass, BGR_6, 
                            bio1, bio12, bio13, bio14), copy = TRUE) %>% 
  left_join(tbl(db, "RAUMDATEN_Z9_NDEP") %>% 
              dplyr::select(aID_STAO, NTOT2007), copy = TRUE)

for(i in 1:nrow(sites)) {
  sites[i, "Turnover"] <- sim(pl %>% filter(aID_STAO == sites$aID_STAO[i]) %>% dplyr::select(aID_KD, aID_SP, Occ), method = "simpson", 
      listin = TRUE, listout = TRUE)$simpson %>% mean
}
  
#------------------------------------------------------------------------------------------------------
# Save data
#------------------------------------------------------------------------------------------------------
save(surv, file = "RData/surv.RData")
save(pl, file = "RData/pl.RData")
save(sites, file = "RData/sites.RData")

