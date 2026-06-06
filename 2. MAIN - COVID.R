
library(dplyr)
library(abind)
library(reshape2)
library(ggplot2)
library(writexl)
library(openxlsx)
library(tidyr)
library(data.table)

# Script om resultaten runnen om sectie 3.2 van deelopdracht 
# "Opdracht Werkgroep AG2024 - COVID19 - nieuwere modellen v20240301" tot 2023 
# te runnen. 

# Dependencies: Alvorens deze script gerund wordt, dient "1. MAIN - AG2022.R" 
# tot r. 151 (Projecting mortality in future) ook gedraaid te worden voor 
# de Mu_xt forecasts.

run_name <- "4. AG2026 incl. startpuntmiddeling, COVID 2022-2023"
run_folder <- here("Runs", "Deelopdracht 7", runname)

input_folder <- here("Input data", "Deelopdracht 6")

path_death <- here(input_folder, "Deaths_NL - weekly.txt")
path_expo  <- here(input_folder, "Exposures_NL - weekly_2023_opnieuw_v2.txt")
path_phi  <- here(input_folder, "phi.xlsx")

dfDeath_week <- read.table(path_death,header=TRUE,dec=".",sep="")
dfExposure_week  <- read.table(path_expo,header=TRUE,dec=".",sep="")
dfPhi <- read.xlsx(path_phi, sheet = "Phi")

phi_long <- dfPhi %>%
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "Phi") %>%  
  mutate_if(is.character, as.numeric)  

dfDagen_per_week <- read.xlsx(path_phi, sheet = "Dagen_per_week")
rownames(dfDagen_per_week) <- dfDagen_per_week$Week
dfDagen_per_week$Week <- NULL
dfDagen_per_week_relatief <- dfDagen_per_week/colSums(dfDagen_per_week, na.rm=TRUE) 


vanaf_jaar <- 2022
tot_jaar <- 2025
iAge_treshold_beta <- 0
aantal_jaar <- length(vanaf_jaar:tot_jaar)
max_leeftijd <- 90

# Frakbeta cluster from age 0 until iAge_indiv_from - 1
bClusterAges <- TRUE
iAge_indiv_from <- 41

data <- c()
for(gender in c("Male", "Female")) {
  
  mQ_tx_m <- read.xlsx(here(run_folder, paste0(runname, "_BE_proj_all.xlsx")), sheet = paste("NL",gender, sep=""))
  mQ_tx_m <- apply(mQ_tx_m, 2, as.numeric)
  mMu_tx_m <- -log(1 - t(mQ_tx_m))
  #mMu_tx_m <- -log(1 - objAG_BE_all[["lQg_NL_tx"]][[gender]])
  
  mMu_tx_m_pre_covid <- t(mMu_tx_m[rownames(mMu_tx_m) %in% c(vanaf_jaar:tot_jaar),])
  rownames(mMu_tx_m_pre_covid) <- 0:120
  
  dfExposure_adj <- 
    dfExposure_week %>% 
    right_join(reshape2::melt(mMu_tx_m_pre_covid, varnames = c("Age", "Year")), by=c("Year","Age")) %>% 
    right_join(phi_long, by=c("Week", "Year")) %>%
    mutate_if(is.character, as.numeric)
  dfExposure_adj <- na.omit(dfExposure_adj)
  
  E_xtw <- dfExposure_adj[[gender]] * dfExposure_adj$value * dfExposure_adj$Phi
  
  dfDeath_week<- dfDeath_week %>%
    filter(Year >= vanaf_jaar & Year <= tot_jaar)
  
  D_xtw <- as.numeric(as.character(dfDeath_week[[gender]]))
  
  Year <- as.numeric(dfDeath_week$Year)
  Age <- as.numeric(as.character(dfDeath_week$Age))
  Week <- as.numeric(as.character(dfDeath_week$Week))
  
  data <- rbind(data, data.frame(year=Year, age=Age, week=Week, dxt=D_xtw,ext=E_xtw,gender=gender))
}

Data_weekly_male <- data %>% 
  filter(gender=="Male") %>% 
  filter(age <= max_leeftijd) %>%
  filter(year >= vanaf_jaar) %>%
  filter(year <= tot_jaar)

Data_weekly_female <- data %>% 
  filter(gender=="Female") %>% 
  filter(age <= max_leeftijd) %>%
  filter(year >= vanaf_jaar) %>%
  filter(year <= tot_jaar)

source("./Scripts/2. Fit AG-model 2022 - weekmodel.R")
lRes_kalibratie_m <- fnFitWeekmodel_withClustering(
  Data_weekly_male,
  iAge_treshold_beta,
  dRelTol = 1e-18,
  iStep = 2,
  bClusterAges = bClusterAges,
  iAge_indiv_from = iAge_indiv_from
)

lRes_kalibratie_v <- fnFitWeekmodel_withClustering(
  Data_weekly_female,
  iAge_treshold_beta,
  dRelTol = 1e-18,
  iStep = 2,
  bClusterAges = bClusterAges,
  iAge_indiv_from = iAge_indiv_from
)

# lRes_kalibratie_m <- fnFitWeekmodel(Data_weekly_male, iAge_treshold_beta, dRelTol=1e-18, iStep=2)
# lRes_kalibratie_v <- fnFitWeekmodel(Data_weekly_female, iAge_treshold_beta, dRelTol=1e-18, iStep=2)

#pas normalisatie toe
const.beta_m <- sum(lRes_kalibratie_m$frakbeta[56:(max_leeftijd+1)])
const.beta_v <- sum(lRes_kalibratie_v$frakbeta[56:(max_leeftijd+1)])

#const.beta_m <- sqrt(sum(lRes_kalibratie_m$frakbeta[56:(max_leeftijd+1)]^2))
#const.beta_v <- sqrt(sum(lRes_kalibratie_v$frakbeta[56:(max_leeftijd+1)]^2))

vFrakbeta_m <- lRes_kalibratie_m$frakbeta/const.beta_m
vFrakbeta_v <- lRes_kalibratie_v$frakbeta/const.beta_v

# FrakKappa matrix (54 x #years) contains non-active elements for years with 52 vs 53 months 
mFrakKappa_m <- lRes_kalibratie_m$frakkappa * const.beta_m
mFrakKappa_v <- lRes_kalibratie_v$frakkappa * const.beta_v

dtData_CBS_sel <- readRDS("dtData_CBS_sel_zp.rds")
dtData_CBS_sel <- as.data.table(dtData_CBS_sel)
dtData_CBS_sel <- dtData_CBS_sel %>% 
  filter(Year >= vanaf_jaar) %>%
  filter(Year <= tot_jaar) %>%
  filter(Age >= iAge_treshold_beta) %>%
  filter(Age <= 90)

#write.xlsx(dtData_CBS_sel, file="test_b.xlsx")

# Initialize output objects
dtKappa_week_M <- dtBeta_M <- dtKappa_week_F <- dtBeta_F <- NULL
dtFrakX <- dtFrakB <- dtFrakX_pre <- dtFrakB_pre <- NULL
vG <- dtData_CBS_sel[ , unique(Gender)]
dtYearWeek <- data.table(YearWeek = sort(dtData_CBS_sel[,unique(YearWeek)]))
dtData_CBS_sel[ , Bx := NA_real_]
dtData_CBS_sel[ , Kt := NA_real_]
dtData_CBS_sel[ , FrakBx := NA_real_]
dtData_CBS_sel[ , FrakXt := NA_real_]

# Step 1: Estimate age and period effects on NL using weekly data from the years 2022 and 2023
vY_temp <- unique(dtData_CBS_sel$YearWeek)
vX_temp <- unique(dtData_CBS_sel$Age)

dtKappa_week_M <- rbind(dtKappa_week_M, data.table(Gender = "Male", YearWeek = vY_temp, Kappa = reshape2::melt(mFrakKappa_m)$value))
dtKappa_week_F <- rbind(dtKappa_week_F, data.table(Gender = "Female", YearWeek = vY_temp, Kappa = reshape2::melt(mFrakKappa_v)$value))
dtKappa_week_NL <- rbind(dtKappa_week_M, dtKappa_week_F)

dtBeta_M <- rbind(dtBeta_M, data.table(Gender = "Male", Age = vX_temp, Beta = vFrakbeta_m[(iAge_treshold_beta+1):(max_leeftijd+1)]))
dtBeta_F <- rbind(dtBeta_F, data.table(Gender = "Female", Age = vX_temp, Beta = vFrakbeta_v[(iAge_treshold_beta+1):(max_leeftijd+1)]))
dtBeta_NL <- rbind(dtBeta_M, dtBeta_F)

for(gender in vG) {
  # gender = vG[1]
  for(age in vX_temp) dtData_CBS_sel[Gender==gender & Age==age, Bx := dtBeta_NL[Gender==gender & Age==age, Beta]]
  for(yearweek in vY_temp) dtData_CBS_sel[Gender==gender & YearWeek==yearweek, Kt := dtKappa_week_NL[Gender==gender & YearWeek==yearweek, Kappa]]
}

fnFindBx <- function(dtData_temp) { uniroot(f = fnBx_LHS, interval=c(-1,1), tol = 1E-18, dtData_temp = dtData_temp)$root}
fnBx_LHS <- function(dtData_temp, dFrakBx) {
  dtData_inner <- data.table::copy(dtData_temp)
  dtData_inner <- dtData_inner[ , j = list(INNER = sum(Nfrac * (exp(dFrakBx * FrakXt) - Phi_fit * exp(Bx * Kt)))),
                                by = c("Year", "Mu")]
  return(dtData_inner[ , sum(Mu * INNER)])
}

# Step 2: Aggregate week effects into annual effects  
for(gender in vG) {
  # gender = vG[1]
  # Compute mathfrakX_t
  for(year in vanaf_jaar:tot_jaar) {
    # year = 2020
    dtData_temp <- data.table::copy(dtData_CBS_sel[Gender==gender & Year==year, ])
    dtData_temp[ , INNER := Nfrac * Phi_fit * exp(Bx * Kt) ]
    dtData_aggr2 <- dtData_temp[ , j = list(INNER_SUM = sum(INNER)), by = Age ]
    dtData_aggr <- dtData_temp[ , j = list(INNER = sum( Nfrac  * Phi_fit * exp(Bx * Kt))), by=Age]
    dtData_aggr2[ , sum(log(INNER_SUM))]
    dtFrakX_pre <- rbind(dtFrakX_pre, data.table(Year = year,
                                                 Gender = gender,
                                                 FrakXt = dtData_aggr[ , sum(log(INNER))]))
    dtFrakX <- rbind(dtFrakX, data.table(Year = year,
                                         Gender = gender,
                                         FrakXt = dtData_aggr[ , sum(log(INNER))]))
    dtData_CBS_sel[Gender==gender & Year==year, FrakXt := dtData_aggr[ , sum(log(INNER))] ]
  }
  
  # Compute mathfrakB_t
  for(age in vX_temp) {
    # age = vX_temp[1]
    dtData_temp <- data.table::copy(dtData_CBS_sel[Gender==gender & Age==age,])
    dtFrakB_pre <- rbind(dtFrakB_pre, data.table(Age = age,
                                                 Gender = gender,
                                                 FrakBx =  (dtData_temp)))
    dtFrakB <- rbind(dtFrakB, data.table(Age = age,
                                         Gender = gender,
                                         FrakBx = fnFindBx(dtData_temp)))
  }
}
# Apply constraints
for(gender in vG) {
  # gender = vG[1]
  vFrakB <- dtFrakB[Gender==gender, FrakBx]
  dtFrakB[Gender==gender, FrakBx := FrakBx / sum(vFrakB)]
  dtFrakX[Gender==gender, FrakXt := FrakXt * sum(vFrakB)]
}
dtData_CBS_sel[ , FrakBx := NULL][ , FrakXt := NULL]
dtData_CBS_sel <- merge(x=dtData_CBS_sel, y=dtFrakB, by=c("Gender", "Age"))
dtData_CBS_sel <- merge(x=dtData_CBS_sel, y=dtFrakX, by=c("Gender", "Year"))

dfFrakKappa_M <- data.frame(
  1:length(mFrakKappa_m),
  reshape2::melt(mFrakKappa_m)$value
)
colnames(dfFrakKappa_M) <- c("Tijd", "FrakKappa")

dfFrakKappa_F <- data.frame(
  1:length(mFrakKappa_v),
  reshape2::melt(mFrakKappa_v)$value
)
colnames(dfFrakKappa_F) <- c("Tijd", "FrakKappa")

dfFrakBeta_M <- data.frame(
  iAge_treshold_beta:max_leeftijd,
  vFrakbeta_m[(iAge_treshold_beta+1):(max_leeftijd+1)]
)
colnames(dfFrakBeta_M) <- c("Leeftijd", "FrakBeta")

dfFrakBeta_F <- data.frame(
  iAge_treshold_beta:max_leeftijd,
  vFrakbeta_v[(iAge_treshold_beta+1):(max_leeftijd+1)]
)
colnames(dfFrakBeta_F) <- c("Leeftijd", "FrakBeta")

# plot(dfFrakKappa_M, type='l')
# grid(nx = 4, # X-axis divided in two sections
#      ny = 6, # Y-axis divided in three sections
#      lty = 2, col = "gray", lwd = 2)
# plot(dfFrakKappa_F, type='l')
# plot(dfFrakBeta_M, type='l')
# plot(dfFrakBeta_F, type='l')

mFrakX_t <- dtData_CBS_sel[(Week==1 & Age==55), .(Gender, Year, FrakXt)]
FrakBeta_long <- dtData_CBS_sel[(Week==1 & Year==tot_jaar), .(Age, Gender, FrakBx)]
print(mFrakX_t)


# Afronden op 15 signifante cijfers
dfFrakKappa_M <- dfFrakKappa_M %>% mutate(across(where(is.numeric), ~ signif(., digits = 15)))
dfFrakKappa_F <- dfFrakKappa_F %>% mutate(across(where(is.numeric), ~ signif(., digits = 15)))
dfFrakBeta_M <- dfFrakBeta_M %>% mutate(across(where(is.numeric), ~ signif(., digits = 15)))
dfFrakBeta_F <- dfFrakBeta_F %>% mutate(across(where(is.numeric), ~ signif(., digits = 15)))
mFrakX_t <- mFrakX_t %>% mutate(across(where(is.numeric), ~ signif(., digits = 15)))
FrakBeta_long <- FrakBeta_long %>% mutate(across(where(is.numeric), ~ signif(., digits = 15)))

output_covid_to_excel <- list("FrakKappa M" = dfFrakKappa_M, 
                              "FrakKappa F" = dfFrakKappa_F,
                              "FrakBeta M" = dfFrakBeta_M,
                              "FrakBeta F" = dfFrakBeta_F, 
                              "FrakX" = mFrakX_t, 
                              "FrakBx" = FrakBeta_long)


browser()
write.xlsx(output_covid_to_excel, file=here(run_folder, paste0(runname, "_COVID-19 parametersschattingen.xlsx")))





