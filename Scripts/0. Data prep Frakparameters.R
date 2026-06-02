
path <- paste0("Input data/")
sMainDir <- "C:/Users/zeusp/OneDrive/Documents/Werkgroep Prognosetafels/Sterftetafel AG2024/5. Rapportage"
setwd(sMainDir)

path_death <- paste0(path, "/Deaths_NL - weekly.txt")
path_expo  <- paste0(path, "/Exposures_NL - weekly.txt")
path_phi  <- paste0(path, "/phi.txt")

dfDeath_week <- read.table(path_death,header=TRUE,dec=".",sep="")
dfExposure_week  <- read.table(path_expo,header=TRUE,dec=".",sep="")
dfPhi <- read.xlsx("Input data/Phi.xlsx", sheet = "Phi")

deaths_long <- dfDeath_week %>%
  pivot_longer(cols = c("Male", "Female"), names_to = "Gender", values_to = "D")
exposure_long <- dfExposure_week %>%
  pivot_longer(cols = c("Male", "Female"), names_to = "Gender", values_to = "E")  

phi_long <- dfPhi %>%
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "Phi_fit") %>%  
  mutate_if(is.character, as.numeric)  

dfDagen_per_week <- read.xlsx("Input data/Phi.xlsx", sheet = "Dagen_per_week")
dfDagen_per_week_relatief <- dfDagen_per_week/colSums(dfDagen_per_week, na.rm=TRUE) 
dfDagen_per_week_relatief$Week <- 0:53
dfDagen_per_week_relatief_long <- dfDagen_per_week_relatief%>%
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "Nfrac") %>%  
  mutate_if(is.character, as.numeric) 

dfDagen_per_week_relatief <- read.xlsx("Input data/Phi.xlsx", sheet = "Nfrac")
dfDagen_per_week_relatief$Week <- 0:53
dfDagen_per_week_relatief_long <- dfDagen_per_week_relatief%>%
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "Nfrac") %>%  
  mutate_if(is.character, as.numeric) 




mQ_tx_m <- read.xlsx("Output/AG2024(EU21NL22) preCOVID/AG2024_BE_proj_all.xlsx", sheet = "NLMale")
mQ_tx_m <- apply(mQ_tx_m, 2, as.numeric)
mMu_tx_m <- -log(1 - t(mQ_tx_m))

mMu_tx_m_pre_covid <- t(mMu_tx_m[rownames(mMu_tx_m) %in% c(2020, 2021, 2022, 2023),])
rownames(mMu_tx_m_pre_covid) <- 0:120
mMu_tx_m_pre_covid <- melt(mMu_tx_m_pre_covid)
names(mMu_tx_m_pre_covid) <- c("Age", "Year", "Mu") 
mMu_tx_m_pre_covid$Gender <- "Male"

mQ_tx_v <- read.xlsx("Output/AG2024(EU21NL22) preCOVID/AG2024_BE_proj_all.xlsx", sheet = "NLFemale")
mQ_tx_v <- apply(mQ_tx_v, 2, as.numeric)
mMu_tx_v <- -log(1 - t(mQ_tx_v))

mMu_tx_v_pre_covid <- t(mMu_tx_v[rownames(mMu_tx_v) %in% c(2020, 2021, 2022, 2023),])
rownames(mMu_tx_v_pre_covid) <- 0:120
mMu_tx_v_pre_covid <- melt(mMu_tx_v_pre_covid)
names(mMu_tx_v_pre_covid) <- c("Age", "Year", "Mu") 
mMu_tx_v_pre_covid$Gender <- "Female"

mMu_tx_pre_covid <- rbind(mMu_tx_m_pre_covid, mMu_tx_v_pre_covid)

Data <- deaths_long[, c("Age", "Week", "Year", "Gender", "D")] %>%
          right_join(exposure_long[, c("Age", "Week", "Year", "Gender", "E")], by=c("Age", "Gender", "Week", "Year")) %>%
          right_join(phi_long, by=c("Week", "Year")) %>%
          right_join(mMu_tx_pre_covid, by=c("Age", "Year", "Gender")) %>%
          right_join((dfDagen_per_week_relatief_long), by=c("Week", "Year"))

Data$YearWeek <- as.numeric(paste0(Data$Year, sprintf("%02d", Data$Week)))



saveRDS(na.omit(Data), "dtData_CBS_sel_zp.rds")





