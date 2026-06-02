library(openxlsx)
library(readr)
library(stringr)
library(dplyr)
library(tidyr)
rm(list = ls())

# deze parameters moeten geinitialiseerd worden vlak voordat ze gebruikt worden in de sourcecode
# dit is nodig omdat een source() aanroep de instellingen van deze parameters verwijdert
# switch_load_CBS <- TRUE
# switch_merge_CBS_HMD <- TRUE
# switch_load_EuroStat <- TRUE
# switch_merge_EuroStat_HMD <- TRUE

dir_data <- file.path("C:/Users/cmout/OneDrive - Goudse/Documenten/AG-AI/Prognosetafel/AG2024/Data/EU21NL22/2024-05 C43/Input_HMD")
dir_output_HMD <- file.path("C:/Users/cmout/OneDrive - Goudse/Documenten/AG-AI/Prognosetafel/AG2024/Data/EU21NL22/2024-05 C43/Output_HMD")
file_list_in <-            c("AUT", "BEL", "DNK", "FIN", "FRACNP", "DEUTNP", "DEUTE", "DEUTW", "ISL", "IRL", "LUX", "NLD", "NOR", "SWE", "CHE", "GBR_NP", "GBRTENW") #land-id van bestanden in HMD

# file_list_out en EuroStat_land_list eisen herhaling van initialisatie nadat een source() is aangeroepen
# daarom is het nodig deze initialisatie te herhalen na een source() aanroep
# file_list_out_zonder_GE <- c("AUS", "BE" , "DNK", "FIN", "FR"    ,       "GEE"  , "GEW"  , "ICE", "IRE", "LUX", "NL" , "NOR", "SWE", "SWI", "UK"    , "EW") #land-id voor bestanden uit HMD, GEE en GEW nog niet samengevoegd
file_list_out <-           c("AUS", "BE" , "DNK", "FIN", "FR"    , "GEP"   , "GEE"  , "GEW"  , "ICE", "IRE", "LUX", "NL" , "NOR", "SWE", "SWI", "UK"    , "EW") #land-id voor bestanden uit HMD 
EuroStat_land_list <-      c('AT' , 'BE' , 'DK' , 'FI' , 'FR'    , 'DE'    , 'DEE'  , 'DEW'  , 'IS' , 'IE' , 'LU' , 'NL' , 'NO' , 'SE' , 'CH' , 'GB'    , 'EW') #2-letterig land-id in EuroStat, 'EW' = Engeland + Wales is geen 2 letterige landcode
lft_aggr_lb <- 105
year_lb <- 1970
year_ub <- 2019
year_DEUTE_DEUTW_merged <- 1990

for (datafile in c("Deaths", "Exposures")) {
  for (i in 1:length(file_list_in)) {
    file <- file.path(dir_data, paste0(file_list_in[i], ".", datafile,  "_1x1.txt"))
    temp <- readr::read_table(file, col_names = TRUE, skip = 2) %>% filter(Year >= year_lb & Year <= year_ub)
    temp2 <- 
      temp %>% filter(temp$Age %in% c(c(lft_aggr_lb:109), "110+")) %>% group_by(Year) %>% summarise(Female = sum(Female), Male = sum(Male), Total = sum(Total))
    temp2$Age <- lft_aggr_lb
    temp2 <- relocate(temp2, Age, .after = Year)
    temp3 <- rbind(temp2, temp %>% filter(temp$Age %in% c(0:lft_aggr_lb - 1)))
    temp3$Age <- as.integer(temp3$Age)
    temp3$Female <- format(round(temp3$Female, 2), nsmall = 2)
    temp3$Male <-format(round(temp3$Male, 2), nsmall = 2)
    temp3$Total <-format(round(temp3$Total, 2), nsmall = 2)
    temp3 <- temp3 %>% arrange(Year, Age)
    filename_out <- file.path(dir_output_HMD, paste0(datafile, "_"  , file_list_out[i], "_HMD", ".txt"))
    write_delim(temp3, filename_out, delim = "\t")
  }
}

# voeg West-Duitsland (= GEW) toe aan Duitsland (= GEP) voor de jaren voor 1990
for (datafile in c("Deaths", "Exposures")) {
  file_P <- file.path(dir_output_HMD, paste0(datafile, "_GEP", "_HMD", ".txt"))
  temp_P <- readr::read_table(file_P, col_names = TRUE)
  file_W <- file.path(dir_output_HMD, paste0(datafile, "_GEW", "_HMD", ".txt"))
  temp_W <- readr::read_table(file_W, col_names = TRUE)
  temp <- rbind(temp_W[temp_W$Year < year_DEUTE_DEUTW_merged,], temp_P)
  # temp <- temp_E %>% inner_join(temp_W, by = c("Year", "Age"), suffix = c("_E", "_W")) %>% 
  # mutate(Female = Female_E + Female_W, Male = Male_E + Male_W, Total = Total_E + Total_W) %>% 
  # subset(select = c(Year, Age, Female, Male, Total)) %>%
  # arrange(Year, Age)
  temp$Age <- as.integer(temp$Age)
  temp$Female <- format(round(temp$Female, 2), nsmall = 2)
  temp$Male <-format(round(temp$Male, 2), nsmall = 2)
  temp$Total <-format(round(temp$Total, 2), nsmall = 2)
  filename_out <- file.path(dir_output_HMD, paste0(datafile, "_GEP", "_HMD", ".txt"))
  write_delim(temp, filename_out, delim = "\t")
}  


# haal de CBS data op om de nog ontbrekende jaren in de HMD aan te vullen met de jaren afkomstig van het CBS 
sourcefile_CBS <- "C:/Users/cmout/OneDrive - Goudse/Documenten/AG-AI/Prognosetafel/AG2024/Data/scripts/20240501_extract_Statline_data_v01.R"
switch_load_CBS <- TRUE
if(switch_load_CBS) source(sourcefile_CBS)

switch_merge_CBS_HMD <- TRUE
if(switch_merge_CBS_HMD) {
  dir_output <- file.path("C:/Users/cmout/OneDrive - Goudse/Documenten/AG-AI/Prognosetafel/AG2024/Data/EU21NL22/2024-05 C43/Output")
  dir_output_HMD <- file.path("C:/Users/cmout/OneDrive - Goudse/Documenten/AG-AI/Prognosetafel/AG2024/Data/EU21NL22/2024-05 C43/Output_HMD")
  dir_data_CBS <- "C:/Users/cmout/OneDrive - Goudse/Documenten/AG-AI/Prognosetafel/AG2024/Data/EU21NL22/2024-05 C43/Output_CBS"
  
  for (datafile in c("Deaths", "Exposures")) {
    file_HMD <- file.path(dir_output_HMD, paste0(datafile, "_NL_HMD", ".txt"))
    temp_HMD <- readr::read_table(file_HMD, col_names = TRUE)
    file_CBS <- file.path(dir_data_CBS, paste0(datafile, "_NL_CBS", ".txt"))
    temp_CBS <- readr::read_table(file_CBS, col_names = TRUE)  
    temp_CBS <- temp_CBS %>% anti_join(temp_HMD, join_by(Year == Year, Age == Age))
    temp_HMD <- rbind(temp_HMD, temp_CBS)
    temp_HMD$Female <- format(round(temp_HMD$Female, 2), nsmall = 2)
    temp_HMD$Male <- format(round(temp_HMD$Male, 2), nsmall = 2)
    temp_HMD$Total <- format(round(temp_HMD$Total, 2), nsmall = 2)
    filename_out <- file.path(dir_output, paste0(datafile, "_NL", ".txt"))
    write_delim(temp_HMD, filename_out, delim = "\t")
    filename_out <- file.path(dir_output, paste0(datafile, "_NL"  , ".xlsx"))
    write.xlsx(list("NL" = temp_HMD), filename_out, col.Names = TRUE, overwrite = TRUE)
  }
}

# haal de EuroStat data op en voeg die toe aan de HMD bestanden
# dit toevoegen hoeft niet te gebeuren voor NL, en voor UK en EW
# mag wel, maar zal niks toevoegen aan de betreffende bestanden
sourcefile_EuroStat <- "C:/Users/cmout/OneDrive - Goudse/Documenten/AG-AI/Prognosetafel/AG2024/Data/scripts/20240501_extract_EuroStat_allCountries_v01.R"
switch_load_EuroStat <- TRUE
if(switch_load_EuroStat) source(sourcefile_EuroStat)

switch_merge_EuroStat_HMD <- TRUE
if(switch_merge_EuroStat_HMD) {
  dir_output <- file.path("C:/Users/cmout/OneDrive - Goudse/Documenten/AG-AI/Prognosetafel/AG2024/Data/EU21NL22/2024-05 C43/Output")
  dir_output_HMD <- file.path("C:/Users/cmout/OneDrive - Goudse/Documenten/AG-AI/Prognosetafel/AG2024/Data/EU21NL22/2024-05 C43/Output_HMD")
  dir_data_EuroStat <- "C:/Users/cmout/OneDrive - Goudse/Documenten/AG-AI/Prognosetafel/AG2024/Data/EU21NL22/2024-05 C43/Output_EuroStat"
  file_list_out <- c("AUS", "BE" , "DNK", "FIN", "FR"    , "GEP"   , "GEE"  , "GEW"  , "ICE", "IRE", "LUX", "NL" , "NOR", "SWE", "SWI", "UK"    , "EW") #land-id voor bestanden uit HMD 
  EuroStat_land_list <-      c('AT' , 'BE' , 'DK' , 'FI' , 'FR'    , 'DE'    , 'DEE'  , 'DEW'  , 'IS' , 'IE' , 'LU' , 'NL' , 'NO' , 'SE' , 'CH' , 'GB'    , 'EW') #2-letterig land-id in EuroStat, 'EW' = Engeland + Wales is geen 2 letterige landcode
  
  for (datafile in c("Deaths", "Exposures")) {
    for (i in 1:length(file_list_out)) {
      file_HMD <- file.path(dir_output_HMD, paste0(datafile, "_", file_list_out[i], "_HMD", ".txt"))
      temp_HMD <- readr::read_table(file_HMD, col_names = TRUE)
      if(!(EuroStat_land_list[i] %in% c('DEE', 'DEW', 'EW'))) {
        file_EuroStat <- file.path(dir_data_EuroStat, paste0(datafile, "_", EuroStat_land_list[i], "_EuroStat", ".txt"))
        temp_EuroStat <- readr::read_table(file_EuroStat, col_names = TRUE)
        temp_EuroStat$Age <- as.integer(temp_EuroStat$Age)
        temp_EuroStat$Year <- as.integer(temp_EuroStat$Year)
        temp_EuroStat <- temp_EuroStat %>% anti_join(temp_HMD, join_by(Year == Year, Age == Age))
        temp_HMD <- rbind(temp_HMD, temp_EuroStat)
      }
      temp_HMD$Female <- format(round(temp_HMD$Female, 2), nsmall = 2)
      temp_HMD$Male <- format(round(temp_HMD$Male, 2), nsmall = 2)
      temp_HMD$Total <- format(round(temp_HMD$Total, 2), nsmall = 2)
      filename_out <- file.path(dir_output, paste0(datafile, "_", file_list_out[i], ".txt"))
      if(EuroStat_land_list[i] != 'NL') write_delim(temp_HMD, filename_out, delim = "\t")
      filename_out <- file.path(dir_output, paste0(datafile, "_", file_list_out[i], ".xlsx"))
      if(EuroStat_land_list[i] != 'NL') write.xlsx(list("Blad1" = temp_HMD), filename_out, col.Names = TRUE, overwrite = TRUE)
    }
  }
}
