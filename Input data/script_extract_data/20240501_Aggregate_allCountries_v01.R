library(readr)
library(tidyr)
library(openxlsx)
library(dplyr)

rm(list = ls())

dir_output <- file.path("C:/Users/cmout/OneDrive - Goudse/Documenten/AG-AI/Prognosetafel/AG2024/Data/EU21NL22/2024-05 C43/Output")
file_list_out <- c("AUS", "BE", "DNK", "FIN", "FR", "GEP", "ICE", "IRE", "LUX", "NL", "NOR", "SWE", "SWI", "UK") #landen om TOTAAL-EU te maken

for (datafile in c("Deaths", "Exposures")) {
  for (i in 1:length(file_list_out)) {
    file <- file.path(dir_output, paste0(datafile, "_", file_list_out[i], ".txt"))
    temp <- readr::read_table(file, col_names = TRUE)
    temp$landcode <- file_list_out[i]
    temp <- temp %>% relocate(landcode, .before = Year)
    if(i > 1) temp2 <- rbind(temp, temp2) else temp2 <- temp
  }
  temp2 <- temp2 %>% group_by(Year, Age) %>% summarise(Female = sum(Female), Male = sum(Male), Total = sum(Total)) %>% arrange(Year, Age)
  temp2$Age <- as.integer(temp2$Age)
  temp2$Female <- format(round(temp2$Female, 2), nsmall = 2)
  temp2$Male <-format(round(temp2$Male, 2), nsmall = 2)
  temp2$Total <-format(round(temp2$Total, 2), nsmall = 2)
  filename_out <- file.path(dir_output, paste0(datafile, "_TOT"  , ".txt"))
  write_delim(temp2, filename_out, delim = "\t")
  filename_out <- file.path(dir_output, paste0(datafile, "_TOT"  , ".xlsx"))
  write.xlsx(list("TOT" = temp2), filename_out, col.Names = TRUE, overwrite = TRUE)
}

