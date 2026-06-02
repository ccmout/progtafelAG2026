# remotes::install_github("opensdmx/rsdmx")
library(rsdmx)
library(dplyr)
library(tidyverse)
library(stringr)

rm(list = ls())

land_list <- c('AT', 'BE', 'DK', 'FI', 'FR', 'DE', 'IS', 'IE', 'LU', 'NL', 'NO', 'SE', 'CH', 'GB')
dir_output <- file.path("C:/Users/cmout/OneDrive - Goudse/Documenten/AG-AI/Prognosetafel/AG2024/Data/EU21NL22/2024-05 C43/Output_EuroStat")
year_EUR_ub <- 2019
# myUrl <- "https://stats.oecd.org/restsdmx/sdmx.ashx/GetData/MIG/TOT../OECD?startTime=2000&endTime=2011"
# myUrl <- "https://ec.europa.eu/eurostat/api/dissemination/sdmx/2.1/data/NAMA_10_GDP/A.CP_MEUR.B1GQ.BE+LU+NL"
# myUrl <- "https://ec.europa.eu/eurostat/api/dissemination/sdmx/2.1/data/NAMA_10_GDP/A.CP_MEUR.B1GQ.NL"
# myUrl <- "https://ec.europa.eu/eurostat/api/dissemination/sdmx/2.1/dataflow/ESTAT/DEMO_PJAN?references=children"
# myUrl <- "https://ec.europa.eu/eurostat/api/dissemination/sdmx/2.1/codelist/ESTAT/DEMO_PJAN"

myUrl <- "https://ec.europa.eu/eurostat/api/dissemination/sdmx/2.1/data/DEMO_MAGER?startPeriod=2020&endPeriod=2023"
dataset <- readSDMX(myUrl)
stats_C <- as.data.frame(dataset)
stats_C <- stats_C %>% filter(geo %in% land_list)

myUrl <- "https://ec.europa.eu/eurostat/api/dissemination/sdmx/2.1/data/DEMO_MAGEC?startPeriod=2020&endPeriod=2023"
dataset <- readSDMX(myUrl)
stats_D <- as.data.frame(dataset)
stats_D <- stats_D %>% filter(geo %in% land_list)

myUrl <- "https://ec.europa.eu/eurostat/api/dissemination/sdmx/2.1/data/DEMO_PJAN?startPeriod=2020&endPeriod=2023"
dataset <- readSDMX(myUrl)
stats_P <- as.data.frame(dataset)
stats_P <- stats_P %>% filter(geo %in% land_list)


myUrl_FRFX <- "https://ec.europa.eu/eurostat/api/dissemination/sdmx/2.1/data/DEMO_MAGER?startPeriod=2012&endPeriod=2012"
dataset_FRFX <- readSDMX(myUrl_FRFX)
stats_C_FRFX <- as.data.frame(dataset_FRFX)
stats_C_FRFX <- stats_C_FRFX %>% filter(geo %in% c('FR', 'FX'))
corr_C_FRFX <- as.data.frame((stats_C_FRFX %>% filter(geo == 'FX') %>% select(obsValue)) - (stats_C_FRFX %>% filter(geo == 'FR') %>% select(obsValue)))
corr_C_FRFX$geo <- 'FR'
corr_C_FRFX$sex <- stats_C_FRFX$sex[stats_C_FRFX$geo == 'FR']
corr_C_FRFX$age <- stats_C_FRFX$age[stats_C_FRFX$geo == 'FR']
corr_C_FRFX <- corr_C_FRFX %>% relocate(geo, sex, age, .before = obsValue)
stats_C <- stats_C %>% left_join(corr_C_FRFX, join_by(geo == geo, sex == sex, age == age))
stats_C$obsValue.y <- ifelse(is.na(stats_C$obsValue.y), 0, stats_C$obsValue.y)
stats_C$obsValue <- stats_C$obsValue.x + stats_C$obsValue.y
stats_C$obsValue.x <- NULL
stats_C$obsValue.y <- NULL

stats_C$unit <- NULL
stats_C$freq <- NULL
stats_C <- stats_C %>% filter(!(age %in% c("TOTAL", "UNK")))
stats_C$age <- ifelse(stats_C$age == "Y_LT1", 0, ifelse(stats_C$age == "Y_OPEN", 100, as.integer(str_sub(stats_C$age, start = 2))))
stats_C <- pivot_wider(stats_C, names_from = sex, values_from = obsValue)
stats_C <- rename(stats_C, Age = age, Year = obsTime, Female = F, Male = M, Total = T)
stats_C <- relocate(stats_C, Year, .before = Age)
stats_C <- stats_C %>% arrange(geo, Year, Age)
stats_C <- stats_C %>% group_by(geo, Year) %>%
  mutate(F_exposure = (Female - lead(Female, 1))/2,
         M_exposure = (Male - lead(Male, 1))/2,
         T_exposure = (Total - lead(Total, 1))/2)
stats_C$F_exposure[stats_C$Age == 0] <- stats_C$Female[stats_C$Age == 0] - stats_C$Female[stats_C$Age == 1]/2
stats_C$M_exposure[stats_C$Age == 0] <- stats_C$Male[stats_C$Age == 0] - stats_C$Male[stats_C$Age == 1]/2
stats_C$T_exposure[stats_C$Age == 0] <- stats_C$Total[stats_C$Age == 0] - stats_C$Total[stats_C$Age == 1]/2


myUrl_FRFX <- "https://ec.europa.eu/eurostat/api/dissemination/sdmx/2.1/data/DEMO_MAGEC?startPeriod=2012&endPeriod=2012"
dataset_FRFX <- readSDMX(myUrl_FRFX)
stats_D_FRFX <- as.data.frame(dataset_FRFX)
stats_D_FRFX <- stats_D_FRFX %>% filter(geo %in% c('FR', 'FX'))
corr_D_FRFX <- as.data.frame((stats_D_FRFX %>% filter(geo == 'FX') %>% select(obsValue)) - (stats_D_FRFX %>% filter(geo == 'FR') %>% select(obsValue)))
corr_D_FRFX$geo <- 'FR'
corr_D_FRFX$sex <- stats_D_FRFX$sex[stats_D_FRFX$geo == 'FR']
corr_D_FRFX$age <- stats_D_FRFX$age[stats_D_FRFX$geo == 'FR']
corr_D_FRFX <- corr_D_FRFX %>% relocate(geo, sex, age, .before = obsValue)
stats_D <- stats_D %>% left_join(corr_D_FRFX, join_by(geo == geo, sex == sex, age == age))
stats_D$obsValue.y <- ifelse(is.na(stats_D$obsValue.y), 0, stats_D$obsValue.y)
stats_D$obsValue <- stats_D$obsValue.x + stats_D$obsValue.y
stats_D$obsValue.x <- NULL
stats_D$obsValue.y <- NULL

stats_D$unit <- NULL
stats_D$freq <- NULL
stats_D <- stats_D %>% filter(!(age %in% c("TOTAL", "UNK")))
stats_D$age <- ifelse(stats_D$age == "Y_LT1", 0, ifelse(stats_D$age == "Y_OPEN", 100, as.integer(str_sub(stats_D$age, start = 2))))
stats_D <- pivot_wider(stats_D, names_from = sex, values_from = obsValue)
stats_D <- rename(stats_D, Age = age, Year = obsTime, Female = F, Male = M, Total = T)
stats_D <- relocate(stats_D, Year, .before = Age)
stats_D <- stats_D %>% arrange(geo, Year, Age)


myUrl_FRFX <- "https://ec.europa.eu/eurostat/api/dissemination/sdmx/2.1/data/DEMO_PJAN?startPeriod=2012&endPeriod=2012"
dataset_FRFX <- readSDMX(myUrl_FRFX)
stats_P_FRFX <- as.data.frame(dataset_FRFX)
stats_P_FRFX <- stats_P_FRFX %>% filter(geo %in% c('FR', 'FX'))
corr_P_FRFX <- as.data.frame((stats_P_FRFX %>% filter(geo == 'FX') %>% select(obsValue)) - (stats_P_FRFX %>% filter(geo == 'FR') %>% select(obsValue)))
corr_P_FRFX$geo <- 'FR'
corr_P_FRFX$sex <- stats_P_FRFX$sex[stats_P_FRFX$geo == 'FR']
corr_P_FRFX$age <- stats_P_FRFX$age[stats_P_FRFX$geo == 'FR']
corr_P_FRFX <- corr_P_FRFX %>% relocate(geo, sex, age, .before = obsValue)
stats_P <- stats_P %>% left_join(corr_P_FRFX, join_by(geo == geo, sex == sex, age == age))
stats_P$obsValue.y <- ifelse(is.na(stats_P$obsValue.y), 0, stats_P$obsValue.y)
stats_P$obsValue <- stats_P$obsValue.x + stats_P$obsValue.y
stats_P$obsValue.x <- NULL
stats_P$obsValue.y <- NULL

stats_P$unit <- NULL
stats_P$freq <- NULL
stats_P$OBS_FLAG <- NULL
stats_P <- stats_P %>% filter(!(age %in% c("TOTAL", "UNK")))
stats_P$age <- ifelse(stats_P$age == "Y_LT1", 0, ifelse(stats_P$age == "Y_OPEN", 100, as.integer(str_sub(stats_P$age, start = 2))))
stats_P <- pivot_wider(stats_P, names_from = sex, values_from = obsValue)
stats_P <- rename(stats_P, Age = age, Year = obsTime, Female = F, Male = M, Total = T)
stats_P <- relocate(stats_P, Year, .before = Age)
stats_P <- stats_P %>% arrange(geo, Age, Year)
stats_P <- stats_P %>% group_by(geo, Age) %>%
  mutate(F_exposure = (Female + lead(Female, 1))/2, 
         M_exposure = (Male + lead(Male, 1))/2, 
         T_exposure = (Total + lead(Total, 1))/2)

stats_P <- stats_P %>% inner_join(stats_C, by = c("geo", "Year", "Age"), suffix = c("_P", "_C"))
stats_P$F_exposure <- stats_P$F_exposure_P + stats_P$F_exposure_C/6
stats_P$M_exposure <- stats_P$M_exposure_P + stats_P$M_exposure_C/6
stats_P$T_exposure <- stats_P$T_exposure_P + stats_P$T_exposure_C/6
stats_P <- stats_P[!(is.na(stats_P$T_exposure)),]
stats_P$Female_P <- NULL
stats_P$Male_P <- NULL
stats_P$Total_P <- NULL
stats_P$Female_C <- NULL
stats_P$Male_C <- NULL
stats_P$Total_C <- NULL
stats_P$F_exposure_P <- NULL
stats_P$M_exposure_P <- NULL
stats_P$T_exposure_P <- NULL
stats_P$F_exposure_C <- NULL
stats_P$M_exposure_C <- NULL
stats_P$T_exposure_C <- NULL
stats_P$OBS_FLAG <- NULL
colnames(stats_P)[4:6] <- c("Female", "Male", "Total")
stats_P$Female <- as.numeric(format(round(replace_na(stats_P$Female, 0), 2), nsmall = 2))
stats_P$Male <- as.numeric(format(round(replace_na(stats_P$Male, 0), 2), nsmall = 2))
stats_P$Total <- as.numeric(format(round(replace_na(stats_P$Total, 0), 2), nsmall = 2))
stats_P <- stats_P %>% arrange(geo, Year, Age)


for (i in 1:length(land_list)) {
  filename_out <- file.path(dir_output, paste0("Deaths", "_"  , land_list[i], "_EuroStat", ".txt"))
  df_write <- subset(stats_D %>% filter(geo == land_list[i] & Year <= year_EUR_ub), select = -c(geo))
  write_delim(df_write, filename_out, delim = "\t")
  filename_out <- file.path(dir_output, paste0("Exposures", "_"  , land_list[i], "_EuroStat", ".txt"))
  df_write <- subset(stats_P %>% filter(geo == land_list[i] & Year <= year_EUR_ub), select = -c(geo))
  write_delim(df_write, filename_out, delim = "\t")
}


