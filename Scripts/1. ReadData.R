Read_CountryGenderData_List = function(path, country.total.strings, genders=c("F", "M"), ages=NULL, years=NULL)
{
  countriesdata <- list()
  
  for(country in country.total.strings)
  {
    for(gender in genders)
    {
      countriesdata[[paste0(country, "-", gender)]] <- Read_CountryData(path=path, country.string=country, genders=gender, ages=ages, years=years)
    }
  }
  countriesdata
}

Read_CountryData = function(path, country.string, genders=c("F", "M"), ages=NULL, years=NULL)
{
  path_death <- paste0(path, "/Deaths_", country.string, ".txt")
  path_expo  <- paste0(path, "/Exposures_", country.string, ".txt")
  death <- read.table(path_death,header=TRUE,dec=".",sep="")
  expo  <- read.table(path_expo,header=TRUE,dec=".",sep="")
  
  data <- c()
  for(G in genders)
  {
    if(G=="F")
    {
      Dxt <- as.numeric(as.character(death$Female))
      Ext <- as.numeric(as.character(expo$Female))
    }
    else
    {
      Dxt <- as.numeric(as.character(death$Male))
      Ext <- as.numeric(as.character(expo$Male))
    }
    
    Mxt    <- Dxt/Ext
    Year   <- as.numeric(death$Year)
    Age    <- as.numeric(as.character(death$Age))
    Age[is.na(Age)] <- 110  
    cntry  <- rep(country.string, dim(death)[1])
    gnder  <- rep(G, dim(death)[1])
    Cg     <- paste(cntry, gnder, sep=".")
    Cohort <- Year-Age 
    
    data <- rbind(data, data.frame(dxt=Dxt,ext=Ext,mxt=Mxt,year=Year,age=Age,country=cntry,gender=gnder,cg=Cg,cohort=Cohort))
  }

  if(!is.null(ages))
    data <- subset(data, data$age %in% ages)
  if(!is.null(years))
    data <- subset(data, data$year %in% years)
  
  data
}


fnReadData <- function(sPath, Countries, Gender, vX, vY) {

	objTemp <- Read_CountryGenderData_List(path=sPath, country.total.string=Countries, genders=Gender, ages=vX, years=vY)

	#create one large dataframe of all countries data
	countrydf <- do.call(rbind, objTemp)
	
	#extract ages, years
	ages  <- sort(unique(countrydf$age))
	years <- sort(unique(countrydf$year))

	#prepare data for fitting common trend (step 1)
	age.factor  <- as.factor(countrydf$age)
	year.factor <- as.factor(countrydf$year)
	Ext <- as.matrix.xtabs(xtabs(countrydf$ext ~ age.factor + year.factor))
	Ext[Ext==0] <- 1e-10
	Dxt <- as.matrix.xtabs(xtabs(countrydf$dxt ~ age.factor + year.factor))
	Dxt[Dxt==0] <- 1e-10
	dimnames(Dxt) <- dimnames(Ext) <- list(Age=vX, Year=vY)
	list(Ext=Ext, Dxt=Dxt)

}		

fnReadData_output <- function(sPath, Countries, Gender, vX, vY) {
  
  objTemp <- Read_CountryGenderData_List(path=sPath, country.total.string=Countries, genders=Gender, ages=vX, years=vY)
  
  #create one large dataframe of all countries data
  countrydf <- do.call(rbind, objTemp)
  
  #extract ages, years
  ages  <- sort(unique(countrydf$age))
  years <- sort(unique(countrydf$year))
  
  #prepare data for fitting common trend (step 1)
  age.factor  <- as.factor(countrydf$age)
  year.factor <- as.factor(countrydf$year)
  Ext <- as.matrix.xtabs(xtabs(countrydf$ext ~ age.factor + year.factor))
  #Ext[Ext==0] <- 1e-10
  Dxt <- as.matrix.xtabs(xtabs(countrydf$dxt ~ age.factor + year.factor))
  #Dxt[Dxt==0] <- 1e-10
  dimnames(Dxt) <- dimnames(Ext) <- list(Age=vX, Year=vY)
  list(Ext=Ext, Dxt=Dxt)
  
}	


	