  rm(list=ls(all=TRUE))
  
  # The objective of this script is to estimate the AG2016 model
  # using different calibration periods. For every calibration
  # period, we shall save the model as an .rds-file, such that
  # it can be loaded at a later stage.
  
  # Further, we construct best estimate mortality forecasts that
  # can be used to compute remaining life expectancies (RLEs).
  # If we plot the RLEs from sequential calibrations, we can
  # entangle the impact of different starting level and different
  # parameter estimates.
  source("./Scripts/1. ReadData.R")
  source("./Scripts/2. Fit AG-model 2022.R")
  source("./Scripts/3. new startingpoint 2026.R")
  source("./Scripts/4. CalculationLifeExpectancy.R")
  
  # Load packages
  library(systemfit)
  library(MASS)
  library(DescTools)
  #library(xlsx)
  library(mvtnorm)
  # ## toevoeging vanwege xlsx package
  library(writexl)
  library(openxlsx2)
  library(readxl)
  library(here)
  
  ###############################################################################
  ### Define table to run
  ###############################################################################
  # MorTable <- 'EU19NL19(AG2024)'
  # MorTable <- 'EU19NL19(AG2026)'
  MorTable <- 'EU23NL24_zonderUK'
  MorTable_inp <-  MorTable #'EU23NL24_zonderUK' # added for purpose of using other input data, default is equal to MorTable
  #MorTable_inp <-  'EU20NL20 (COVID-SC v10)'
  
  # Choices: 
  # 1)  AG2016: EU14NL15
  # 2)  AG2018: EU16NL17, etc. 
  
    
  if(MorTable == "EU16NL17"){
    eyNL <- 2017; eyEU <- 2016
  } else if(MorTable == "EU17NL17"){
    eyNL <- 2017; eyEU <- 2017
  } else if(MorTable == "EU16NL17"){
    eyNL <- 2017; eyEU <- 2016
  } else if(MorTable == "EU17NL18"){
    eyNL <- 2018; eyEU <- 2017
  } else if(MorTable == "EU18NL19"){
    eyNL <- 2019; eyEU <- 2018
  } else if(MorTable == "EU19NL19"){
    eyNL <- 2019; eyEU <- 2019
  } else if(MorTable == "EU19NL20"){
    eyNL <- 2020; eyEU <- 2019
  } else if(MorTable == "EU20NL20"){
    eyNL <- 2020; eyEU <- 2020
  } else if(MorTable == "EU21NL22"){
    eyNL <- 2019; eyEU <- 2019
  } else if(MorTable == "EU19NL19(AG2024)"){
    eyNL <- 2019; eyEU <- 2019
  } else if(MorTable == "EU19NL19(AG2026)"){
    eyNL <- 2019; eyEU <- 2019
  } else if(MorTable == "EU23NL24_zonderUK"){
    eyNL <- 2024; eyEU <- 2023
  }
    
  ###############################################################################
  ### Initialize variables
  ###############################################################################
  
  ### Ages
  vX	<- 0:90
  vY_EU <- 1970:eyEU
  vY_NL_s <- 1983   # adapted for new model change
  vY_NL <- vY_NL_s:eyNL   #not a variable, but still included for completeness
  # vCountries <- c("NL","GEP","FR","DNK","SWI","LUX","IRE","UK","NOR","BE","FIN","SWE","AUS","ICE") # c("NL","DUW","FR","DNK","SWI","LUX","IRE","UK","NOR","BE","FIN","SWE","AUS","ICE")
  vCountries <- c("NL","GEP","FR","DNK","SWI","LUX","IRE","NOR","BE","FIN","SWE","AUS","ICE")
  sCountry <- c("NL")
  
  # dynamic variables 
  fixed_Sims <- FALSE
  iSims <- 100 # 100          # PUT TO 0 TO USE PREVIOUS SIMULATIONS
  iT_fut <- 231 ## 250
  quantiles <- c(0.50, 0.975, 0.025)  # gewenste kwantielen van levensverwachting
  YoI <- c(2020, 2022)  # gewenste jaren voor levensverwachting o.b.v. leeftijden voor simulatie
  AoI <- c(1,66) #select ages of interest for calculation of life expetancy = age + 1 for purpose of column. Thus 1 = age 0
  vCountries_LE <- c("NL","EU") # select countries for which life expectancy is calculated; always have NL as first
  
  # runname <- "3. incl startpunt, excl Covid"
  # run_folder <- here("Runs", "Raymond validatie", runname)
  run_folder <- here("Runs", "deelopdracht_1c. excl UK, excl startpunt, excl Covid")
  runname <- "EU23NL24_zonderUK"
  if(!dir.exists(run_folder)) dir.create(run_folder, recursive=TRUE)
  
  # Kannisto
  kannisto <- FALSE # TRUE, sinds AG2022 moet deze op FALSE staan
  
  # Startpunt middeling
  # apply_new_startingpoint <- TRUE
  apply_new_startingpoint <- FALSE # voor deze deelopdracht moeten de COVID correcties niet worden gebruikt
  if (apply_new_startingpoint==TRUE){
    s <- 5
    c <- 25
    startpunt_convergentie <- 2019
  }
  
  # Covid
  applyCOVID <- FALSE
  if (applyCOVID==TRUE){
    COVIDjaren   <- 2022:2023
    COVID_parameterbestand <- here(run_folder, paste0(runname, "_COVID-19 parameters.xlsx"))
    eta_discount <- 0.75
    eta_limit <- 0
  }  

  ###############################################################################
  ### Read in Dxt and Ext
  ###############################################################################
  
  # Load mortality data
  lData_EU <- list()
  lData_NL <- list()
  lData_Country <- list()
  sInputFolder <- here("Input data", MorTable_inp)
  
  for(gender in c("M", "F")) {
    # Read data from files	(to check: how many levels upwards in the folder structure?)
    lData_EU[[gender]] <- fnReadData(sPath=sInputFolder, Countries=vCountries, Gender=gender, vX=vX, vY=vY_EU)
    lData_NL[[gender]] <- fnReadData(sPath=sInputFolder, Countries=sCountry,   Gender=gender, vX=vX, vY=vY_NL)
  }
  
  sTableName <- paste0("AG", tail(vY_NL,1)+1)
 
  ### Create the output folder if it does not already exist
  
 
  #########################################
  ### Calibration of the AG-model
  #########################################
  # create the weight vector, allowing to exclude specific years from the calibration (like the COVID period)
  w_EU <- rep_len(1, length(vY_EU))
  names(w_EU) <- vY_EU
  w_EU["2020"] <- w_EU["2021"] <- w_EU["2022"] <- 0
  w_NL <- rep_len(1, length(vY_NL))
  names(w_NL) <- vY_NL
  w_NL["2020"] <- w_NL["2021"] <- w_NL["2022"] <- 0
  
  objAG_est <- fnFitAG2020(lData_EU, lData_NL, w_EU, w_NL, vX, vY_EU, vY_NL, parRest=2)
  #saveRDS(object=objAG_TS_est, file=here(run_folder, paste(runname, "_TS_est.rds", sep="")))

  #########################################
  ### Projecting mortality in future
  #########################################
  if(kannisto==TRUE){
    # objAG_BE <- fnProjectAG_BE_kannisto(objAG=objAG_est, iT_fut=iT_fut)
    objAG_BE_all <- fnProjectAG_BE_all_kannisto(objAG=objAG_est, iT_fut=iT_fut)
    
  } else {
    # objAG_BE <- fnProjectAG_BE_alternatiefB_HellingEnConstante(objAG=objAG_est, iT_fut=iT_fut)
    if(sum(ifelse(w_EU == 0, 1, 0)) > 0) {
      objAG_filt <- objAG_est
      for(gender in c("M", "F")){
        objAG_filt[[gender]]$mExt_EU <- objAG_filt[[gender]]$mExt_EU[, w_EU != 0]
        objAG_filt[[gender]]$mDxt_EU <- objAG_filt[[gender]]$mDxt_EU[, w_EU != 0]
        objAG_filt[[gender]]$vK <- objAG_filt[[gender]]$vK[w_EU != 0]
        objAG_filt[[gender]]$vK_old <- objAG_filt[[gender]]$vK_old[w_EU != 0]
        objAG_filt[[gender]]$mMxt_EU <- objAG_filt[[gender]]$mMxt_EU[, w_EU != 0]
        objAG_filt[[gender]]$mMxt_EU_old <- objAG_filt[[gender]]$mMxt_EU_old[, w_EU != 0]
        objAG_filt[[gender]]$Dxt_EU_hat <- objAG_filt[[gender]]$Dxt_EU_hat[, w_EU != 0]
        objAG_filt[[gender]]$vY_EU <- objAG_filt[[gender]]$vY_EU[w_EU != 0]
        
        objAG_filt[[gender]]$mExt_NL <- objAG_filt[[gender]]$mExt_NL[, w_NL != 0]
        objAG_filt[[gender]]$mDxt_NL <- objAG_filt[[gender]]$mDxt_NL[, w_NL != 0]
        objAG_filt[[gender]]$vKappa <- objAG_filt[[gender]]$vKappa[w_NL != 0]
        objAG_filt[[gender]]$mMxt_NL <- objAG_filt[[gender]]$mMxt_NL[, w_NL != 0]
        objAG_filt[[gender]]$mDxt_NL_hat <- objAG_filt[[gender]]$mDxt_NL_hat[, w_NL != 0]
        objAG_filt[[gender]]$vY_NL <- objAG_filt[[gender]]$vY_NL[w_NL != 0]
      }

      # Estimate the TS model and extract regression results
      obj_filt_TS <- fnEstimateTS_AG2020(objAG=objAG_filt, dRelTol=1E-18)  
      saveRDS(obj_filt_TS, here(run_folder, paste0(runname, "_filt_TS.rds")))
      mCov_filt  <- obj_filt_TS$Sigma
      mRho_filt  <- obj_filt_TS$Rho
      vDrift_filt <- obj_filt_TS$Drift
      vMR_filt    <- c(1, 1, obj_filt_TS$MR)
      A <- diag(vMR_filt)
      A_cumsum <- matrix(data = 0, nrow = 4, ncol = 4)
      A_cumsum <- diag(4)
      start_jaar <- as.character(min(as.integer(names(w_EU)[w_EU == 0]))-1)
      temp_start <- c(objAG_filt[[1]]$vK[start_jaar], objAG_filt[[2]]$vK[start_jaar], objAG_filt[[1]]$vKappa[start_jaar], objAG_filt[[2]]$vKappa[start_jaar])
      names(temp_start) <- c("K_tM", "K_tV", "Kappa_tM", "Kappa_tV")
      for(i in c(names(w_EU)[w_EU == 0])){
        temp <- drop(t(vDrift_filt %*% A_cumsum) + (A %^% (which(i == c(names(w_EU)[w_EU == 0])))) %*% temp_start)
        objAG_filt[[1]]$vK[i] <- temp[1]
        objAG_filt[[2]]$vK[i] <- temp[2]
        objAG_filt[[1]]$vK_old[i] <- objAG_filt[[1]]$vK[i]
        objAG_filt[[2]]$vK_old[i] <- objAG_filt[[2]]$vK[i]
        objAG_filt[[1]]$vKappa[i] <- temp[3]
        objAG_filt[[2]]$vKappa[i] <- temp[4]
        A_cumsum <- A_cumsum + (A %^% (which(i == c(names(w_EU)[w_EU == 0]))))
      }
    }
    objAG_est$M$vK <- objAG_filt$M$vK[sort(names(objAG_filt$M$vK))]
    objAG_est$F$vK <- objAG_filt$F$vK[sort(names(objAG_filt$F$vK))]
    objAG_est$M$vK_old <- objAG_filt$M$vK_old[sort(names(objAG_filt$M$vK_old))]
    objAG_est$F$vK_old <- objAG_filt$F$vK_old[sort(names(objAG_filt$F$vK_old))]
    objAG_est$M$vKappa <- objAG_filt$M$vKappa[sort(names(objAG_filt$M$vKappa))]
    objAG_est$F$vKappa <- objAG_filt$F$vKappa[sort(names(objAG_filt$F$vKappa))]

    for( gender in c("M", "F") ) {
      objAG_est[[gender]]$mMxt_EU <-
        fnMxtFromPars(alpha=objAG_est[[gender]]$vA, beta=objAG_est[[gender]]$vB, kappa=objAG_est[[gender]]$vK, ages=objAG_est[[gender]]$vX,
                      years=objAG_est[[gender]]$vY_EU)
      objAG_est[[gender]]$mMxt_EU_old <- 
        fnMxtFromPars(alpha=objAG_est[[gender]]$vA, beta=objAG_est[[gender]]$vB, kappa=objAG_est[[gender]]$vK_old, 
                      ages=objAG_est[[gender]]$vX, years=objAG_est[[gender]]$vY_EU[objAG_est[[gender]]$vY_EU %in% vY_EU])
      objAG_est[[gender]]$mMxt_NL <- 
        objAG_est[[gender]]$mMxt_EU[ , as.character(vY_NL)] * 
        fnMxtFromPars(alpha=objAG_est[[gender]]$vAlpha, beta=objAG_est[[gender]]$vBeta, kappa=objAG_est[[gender]]$vKappa, ages=vX, years=vY_NL)
      objAG_est[[gender]]$Dxt_EU_hat <- objAG_est[[gender]]$mExt_EU * objAG_est[[gender]]$mMxt_EU_old
      objAG_est[[gender]]$mDxt_NL_hat <- objAG_est[[gender]]$mExt_NL * objAG_est[[gender]]$mMxt_NL
      objAG_est[[gender]]$vY_EU <- objAG_est[[gender]]$vY_EU[objAG_est[[gender]]$vY_EU %in% vY_EU]
      objAG_est[[gender]]$mMxt_EU <- objAG_est[[gender]]$mMxt_EU_old
    }
  }
  objAG_BE_all <- fnProjectAG_BE_all_alternatiefB_HellingEnConstante(objAG=objAG_est, iT_fut=iT_fut)
  # objAG_est <- objAG_BE_all$objAG
 
  ### Change startingpoint projected mortaility and middelen
  if (apply_new_startingpoint==TRUE){
    
    lData_NL_incl_age_91 <- list()
    for(gender in c("M", "F")) {
      lData_NL_incl_age_91[[gender]] <- fnReadData(sPath=sInputFolder, Countries=sCountry, Gender=gender, vX=0:91, vY=vY_NL) ### Notice we also read age 91, because it is used for averaging for age 90.
    }
    
    # objAG_BE <- fnApplyNewStartingpoints(s,
    #                                      c,
    #                                      startpunt_convergentie,
    #                                      objAG_BE,
    #                                      lData_NL_incl_age_91,
    #                                      eyNL,
    #                                      vY_NL)
    objAG_BE_all <- fnApplyNewStartingpoints_all(s,
                                                 c,
                                                 startpunt_convergentie,
                                                 objAG_BE_all,
                                                 lData_NL_incl_age_91,
                                                 eyNL,
                                                 vY_NL)
  }

  ### Apply COVID parameters over the Pre-Covid results
  if (applyCOVID==TRUE) {
    frak_beta_m  <- openxlsx::read.xlsx(COVID_parameterbestand, sheet = "Parameterwaarden", colNames=FALSE, cols = 6 , rows = 10:130)
    frak_kappa_m <- openxlsx::read.xlsx(COVID_parameterbestand, sheet = "Parameterwaarden", colNames=FALSE, cols = 11, rows = 60:(59+length(COVIDjaren)))

    #eta_m        <- openxlsx::read.xlsx(COVID_parameterbestand, sheet = "Parameterwaarden", colNames=FALSE, cols = 9 , rows = 68    )
    frak_beta_f  <- openxlsx::read.xlsx(COVID_parameterbestand, sheet = "Parameterwaarden", colNames=FALSE, cols = 18, rows = 10:130)
    frak_kappa_f <- openxlsx::read.xlsx(COVID_parameterbestand, sheet = "Parameterwaarden", colNames=FALSE, cols = 23, rows = 60:(59+length(COVIDjaren)))
    #eta_f        <- openxlsx::read.xlsx(COVID_parameterbestand, sheet = "Parameterwaarden", colNames=FALSE, cols = 21, rows = 68    )
    
    objCOVIDparameters <- list()
    objCOVIDparameters[["M"]] <- list(COVIDjaren = COVIDjaren, frak_beta = frak_beta_m, frak_kappa = frak_kappa_m)
    objCOVIDparameters[["F"]] <- list(COVIDjaren = COVIDjaren, frak_beta = frak_beta_f, frak_kappa = frak_kappa_f)
    
    # objAG_BE <- fnApplyCOVIDparameters(objAG_BE, 
    #                                    objCOVIDparameters,
    #                                    eta_discount,
    #                                    eta_limit)
    objAG_BE_all <- fnApplyCOVIDparameters_all(objAG_BE_all, 
                                               objCOVIDparameters,
                                               eta_discount,
                                               eta_limit)
  }

  fnWriteBE_all(objAG_BE=objAG_BE_all, sFile=here(run_folder, paste0(runname, "_BE_proj_all.xlsx")))
  saveRDS(object=objAG_BE_all, file=here(run_folder, paste0(runname, "_BE.rds")))
  fnWritePars(objAG=objAG_est, sFile=here(run_folder, paste0(runname, "_pars.xlsx")))
  saveRDS(object=objAG_est, file=here(run_folder, paste0(runname, "_est.rds")))

  # Calculate projected kappa
  Kappa_proj <- Project_kappa(objAG=objAG_est, iT_fut=iT_fut)
  saveRDS(object=Kappa_proj, file=here(run_folder, paste0(runname, "_kappa.rds")))


  # # Stochastic simulations
  # if (iSims == 0){
  #   objAG_sim <- readRDS(file=paste0(sOutputFolder, sTableName, "_sim.rds", sep=""))
  # } else{
  #   if(kannisto==TRUE){
  #     objAG_sim <- fnProjectAG_sim_kannisto(objAG=objAG_est, iT_fut=iT_fut, iSims=iSims)
  #   }else{
  #     objAG_sim <- fnProjectAG_sim_alternatiefB(objAG=objAG_est, iT_fut=iT_fut, iSims=iSims)
  #   }
  #   saveRDS(object=objAG_sim, file=paste0(sOutputFolder, sTableName, "_sim.rds", sep=""))
  # }
  # 
  # 
  ###############################################################################
  ### Write pdf mortality fit vs actual and prediction for all ages and both genders
  ###############################################################################

#   pdf(file.path(sOutputFolder, "Mortality_fit_vs_prediction.pdf"), height=7, width=15)
#   vG <- c("Males", "Females")
#   vCols <- c( "blue", "magenta")
#   vT_plot <- 1983:2030
#   for(age in 1:length(vX)) {
#     if(age%%4==1) par(mfcol=c(2,4), las=1, mar=c(3,4,3,2))
#     for(gender in 1:2) {
#       # Observed mortality
#       vQ_obs <- 1 - exp( - objAG_est[[gender]]$mDxt_NL[O(vX[age]), ] / objAG_est[[gender]]$mExt_NL[as.character(vX[age]), ] )
#       vT_obs <- as.integer(names(vQ_obs))
#       # Fitted and projected mortality
#       vQ_model <- objAG_BE[[gender]][ , as.character(vX[age])]
#       vT_model <- as.integer(names(vQ_model))
# 
#       # Plot results:
#       plot(x=NULL, y=NULL, xlim=range(vT_plot), ylim=range(vQ_obs[vT_obs %in% vT_plot], vQ_model[vT_model %in% vT_plot]), xlab=NA, ylab=NA, main=NULL)
#       grid(); box()
#       points(x=vT_obs[vT_obs %in% vT_plot], y=vQ_obs[vT_obs %in% vT_plot], pch=16, col=vCols[gender])
#       lines(x=vT_model[vT_model %in% vT_plot], y=vQ_model[vT_model %in% vT_plot], lty=1, col="black", lwd=1.5)
#       legend("topright", legend=paste0(vG[gender], " aged ", vX[age]), adj=0.1, bg="white")
#     }
#   }
#   dev.off(
#   # 
#   # ###############################################################################
#   # ### Output data - Dxt and Ext per country for all years incl. years not used for estimation purposes
#   # ###############################################################################
#   # 
#   # #Set output years
#   # vY_out    <- 1970:eyEU #individual EU countries
#   # #vY_out    <- 1970:2022 #individual EU countries
#   # 
#   # lData_output_aggr <- list()
#   # for(gender in c("M", "F")) {
#   #   # Output EU
#   #   lData_output_aggr[[paste0("Ext_EU_", gender)]] <- fnReadData_output(sPath=sInputFolder, Countries=vCountries, Gender=gender, vX=vX, vY=vY_out)[[1]]
#   #   lData_output_aggr[[paste0("Dxt_EU_", gender)]] <- fnReadData_output(sPath=sInputFolder, Countries=vCountries, Gender=gender, vX=vX, vY=vY_out)[[2]]
#   #   # Output NL
#   #   lData_output_aggr[[paste0("Ext_NL_", gender)]] <- fnReadData_output(sPath=sInputFolder, Countries=sCountry, Gender=gender, vX=vX, vY=c(1970:eyNL))[[1]]
#   #   lData_output_aggr[[paste0("Dxt_NL_", gender)]] <- fnReadData_output(sPath=sInputFolder, Countries=sCountry, Gender=gender, vX=vX, vY=c(1970:eyNL))[[2]]
#   #   for(country in setdiff(vCountries,"NL")) {
#   #     # Output individual EU countries
#   #     lData_output_aggr[[paste0("Ext_", country, "_", gender)]] <- fnReadData_output(sPath=sInputFolder, Countries=country,   Gender=gender, vX=vX, vY=vY_out)[[1]]
#   #     lData_output_aggr[[paste0("Dxt_", country, "_", gender)]] <- fnReadData_output(sPath=sInputFolder, Countries=country,   Gender=gender, vX=vX, vY=vY_out)[[2]]
#   #   }
#   # }
#   # openxlsx::write.xlsx(x=lData_output_aggr, file=paste0(sOutputFolder, sTableName, " - DxtExt_dataset.xlsx"), rowNames=TRUE, colNames=TRUE)
#   # 
#   # ###############################################################################
#   # ### Calculate historical Period Life Expectancy - multiple countries based on kannisto
#   # ###############################################################################
#   # 
#   # nCountries <- length(vCountries_LE)  # define number of selected countries for LE
#   # 
#   # #Initialize output
#   # LE_P_hist <- list(); lData <- list(); lQtx_hist <- list();
#   # 
#   # for(country in vCountries_LE){
#   #   lData <- list()
#   #   for(gender in c("M", "F")) {
#   #     # Get historical qtx
#   #     if(country=="EU"){ sCountry <- vCountries }else{sCountry=country}
#   #     if(country=="NL"){ vY=c(1970:eyEU) }else{ vY = vY_out}
#   #     lData[[paste0(country,"_",gender)]] <- fnReadData(sPath=sInputFolder, Countries=sCountry, Gender=gender, vX=vX, vY=vY)
#   #     lData[[paste0(country,"_",gender)]]$mxt <- lData[[paste0(country,"_",gender)]]$Dxt/lData[[paste0(country,"_",gender)]]$Ext
#   #     lData[[paste0(country,"_",gender)]]$qxt <- 1-exp(-lData[[paste0(country,"_",gender)]]$mxt)
#   # 
#   #     # Close mortality table
#   # 
#   #     qtx_oldAges <- 1 - exp( -(ApplyKannisto_mtx(mtx=t(lData[[paste0(country,"_",gender)]]$mxt), ages=vX, kannisto_nages=30, kannisto_nobs=11)))
#   #     lQtx_hist[[paste0(country,"_",gender)]] <- cbind(t(lData[[paste0(country,"_",gender)]]$qxt), qtx_oldAges)
#   #     dimnames(lQtx_hist[[paste0(country,"_",gender)]])[[2]] <- 0:120
#   #   }
#   # }
#   # #Calculate period Life expectancy historical qxt
#   # for(age in AoI){
#   #   for(i in 1:length(lQtx_hist)){
#   #     name <- names(lQtx_hist)[[i]]
#   #     LE_P_hist[[paste0(age-1)]][[name]] <- fnLifeExpPer(lQtx_hist[[name]],age, age_ult=120)
#   #   }
#   # }
#   # 
#   # lOutput_LE_P_hist <- list();
#   # for(age in AoI){
#   #   LE_P_hist_temp <- data.frame(LE_P_hist[[paste0(age-1)]])
#   #   rownames(LE_P_hist_temp) <- vY_EU
#   #   lOutput_LE_P_hist[[paste0("LE",age-1,"P")]] <- LE_P_hist_temp #data.frame(output_LE_P[[paste0(age-1)]])
#   # }
#   # 
#   # openxlsx::write.xlsx(x=lOutput_LE_P_hist, file=paste0(sOutputFolder, sTableName, "_Period_Life_expectancy_hist.xlsx"), rowNames=TRUE, colNames=TRUE)
#   # 
#   ###############################################################################
#   ### Calculate Life Expectancy - NL & EU
#   ###############################################################################
# 
#   #Calculate period and cohort LE for future years - only for EU and NL (based on parameters AG table )
ages_LE <- 0:99

start_year <- vY_NL[1] #first year 
start_yearEU <- vY_EU[1] 
proj_year <- vY_NL[length(vY_NL)]+1 #first projection year 
proj_yearEU <- vY_EU[length(vY_EU)]+1 #first projection year for EU, equal to NL or -1 depending on the data update 
vY_NL2 <- rownames(objAG_BE_all$lQg_NL_tx$Male) # all years NL 
vY_EU2 <- rownames(objAG_BE_all$lQg_EU_tx$Male) #all years EU

LE_P_male <- t(fnLifeExpPer(objAG_BE_all$lQg_NL_tx$Male[vY_NL2 >=start_year,], age_ult=120))[1:100, 1:118]
rownames(LE_P_male) <- ages_LE
colnames(LE_P_male) <- start_year:(start_year + 117)

LE_C_male <- t(fnLifeExpCoh(objAG_BE_all$lQg_NL_tx$Male[vY_NL2 >=start_year,], age_ult=120))[1:100, 1:118]
rownames(LE_C_male) <- ages_LE
colnames(LE_C_male) <- start_year:(start_year + 117)

LE_P_female <- t(fnLifeExpPer(objAG_BE_all$lQg_NL_tx$Female[vY_NL2 >=start_year,], age_ult=120))[1:100, 1:118]
rownames(LE_P_female) <- ages_LE
colnames(LE_P_female) <- start_year:(start_year + 117)

LE_C_female <- t(fnLifeExpCoh(objAG_BE_all$lQg_NL_tx$Female[vY_NL2 >=start_year,], age_ult=120))[1:100, 1:118] 
rownames(LE_C_female) <- ages_LE
colnames(LE_C_female) <- start_year:(start_year + 117)

write_xlsx(list(Male_P = LE_P_male, 
                          Male_C = LE_C_male, 
                          Female_P = LE_P_female, 
                          Female_C = LE_C_female), 
                     file=here(run_folder, paste0(runname, "_LE.xlsx")), rowNames = TRUE)

# vG <- c("Males", "Females")
vCols <- c( "lightblue", "red", "darkblue", "magenta")
vLty <- c(1,1,2,2)
vT_plot <- 1983:2100
vLabels <- c("Periodic Male", "Periodic Female", "Cohort Male", "Cohort Female")
function_name <- "function_temp" 
for(age in AoI) {
  pdf(file=here(run_folder, paste0("periode_ en cohort_LE_", vX[age], ".pdf")), height=7, width=15)
  # #  if(age%%4==1) par(mfcol=c(2,4), las=1, mar=c(3,4,3,2))
  for(gender in c("Male", "Female")) {
    for(LE_type in c("Periodic", "Cohort")) {
      real_function_name <- paste0("fnLifeExp", substr(LE_type, start = 1, stop = 3))
      arguments <- list(aQtxg = objAG_BE_all$lQg_NL_tx[[gender]][vY_NL2 >=start_year,], age_ult = 120)
      calculated_LE <- t(do.call(get(real_function_name), arguments))[age, 1:118]
      assign(paste0("LE_", LE_type, "_", gender), calculated_LE)
    }
  }
  # Plot results:
  plot_data <- cbind(get("LE_Periodic_Male"), get("LE_Periodic_Female"), get("LE_Cohort_Male"), get("LE_Cohort_Female")) 
  y_range <- range(plot_data)
  par(mar=c(8,4,4,2), xpd = TRUE)
  plot(x=NULL, y=NULL, xlim=range(vT_plot), ylim=y_range, xlab=NA, ylab=NA, main=NULL)
  grid(); box()
  #points(x=vT_obs[vT_obs %in% vT_plot], y=vQ_obs[vT_obs %in% vT_plot], pch=16, col=vCols[gender])
  for(i in 1:4){lines(x=vT_plot, y=plot_data[,i], lty=vLty[i], col=vCols[i], lwd=1.5)}
  legend("bottom", inset=c(0, -0.32),legend=paste0(vLabels, " aged ", vX[age]), col=vCols, lty=vLty, lwd=1.5, adj=0.1, bg="white", xpd = TRUE)
  dev.off()
}

# vect <-
# rownames(objAG_BE_all$lQg_NL_tx$Male)[vY_NL2 >= start_year] # vect_EU <-
# rownames(objAG_BE_all$lQg_EU_tx$Male)[vY_EU2 >= start_year] #




# LE_P[[paste0(age-1)]][["NL_M"]] <- c(rep(NA, length(vY_EU2) -
# length(LE_P_NL_M)), LE_P_NL_M) # LE[[paste0(age-1)]][["NL_M"]] <- c(rep(NA,
# length(vY_EU2) - length(LE_NL_M)), LE_NL_M) # LE_P[[paste0(age-1)]][["NL_F"]]
# <- c(rep(NA, length(vY_EU2) - length(LE_P_NL_M)), LE_P_NL_F) #
# LE[[paste0(age-1)]][["NL_F"]] <- c(rep(NA, length(vY_EU2) - length(LE_NL_M)),
# LE_NL_F)
# 
# # LE_P[[paste0(age-1)]][["EU_M"]] <-
# fnLifeExpPer(objAG_BE_all$lQg_EU_tx$Male[vY_EU2 >= start_yearEU,],age,age_ult=120) # LE[[paste0(age-1)]][["EU_M"]] <-
# fnLifeExpCoh(objAG_BE_all$lQg_EU_tx$Male[vY_EU2 >= start_yearEU,],age,age_ult=120) # LE_P[[paste0(age-1)]][["EU_F"]] <-
# fnLifeExpPer(objAG_BE_all$lQg_EU_tx$Female[vY_EU2 >= start_yearEU,],age,age_ult=120) # LE[[paste0(age-1)]][["EU_F"]]
# <-fnLifeExpCoh(objAG_BE_all$lQg_EU_tx$Female[vY_EU2 >= start_yearEU,],age,age_ult=120)
# 
# 
#   for(age in AoI){
#     LE_P_temp <- data.frame(LE_P[[paste0(age-1)]])
#     rownames(LE_P_temp) <- vY_EU2
#     lOutput_LE[[paste0("LE",age-1,"P")]] <- LE_P_temp #data.frame(output_LE_P[[paste0(age-1)]])
#     LE_temp <- data.frame(LE[[paste0(age-1)]])
#     rownames(LE_temp) <- vY_EU2
#     lOutput_LE[[paste0("LE",age-1)]] <- LE_temp #data.frame(output_LE[[paste0(age-1)]])
#   }
# openxlsx::write.xlsx(x=lOutput_LE, file=paste0(sOutputFolder, sTableName,
# "_Life_expectancy.xlsx"), rowNames=TRUE, colNames=TRUE) #
# saveRDS(object=lOutput_LE, file=paste0(sOutputFolder, sTableName, "_LE.rds",
# sep="")) # # # ############################## # # y_o_i <- 2021:2031      #
# years on interest for summary table # # Create summary table for Covid impact
# (only forecasts NL) # output_sum <- matrix(NA,length(y_o_i),length(AoI)*2*2)
# # #omvang tabel = aantal leeftijden * aantal geslacht & aantal soorten LE #
# #   # for (i in 1:2){   # male/female indicator; 1 = male, 2 = female
# #   #   for (age_i in 1:length(AoI)){
# #   #     output_sum[,age_i] <- LE[[age_i]]$NL_M[(y_o_i-proj_year+1)]
# #   #     output_sum[,age_i+2*length(AoI)] <- LE_P[[age_i]]$NL_M[(y_o_i-proj_year+1)]
# #   #     output_sum[,i+age_i] <- LE[[age_i]]$NL_F[(y_o_i-proj_year+1)]
# #   #     output_sum[,i+age_i+2*length(AoI)] <- LE_P[[age_i]]$NL_F[(y_o_i-proj_year+1)]
#   #
#   #   }
#   # }
# # dimnames(output_sum) <- list(y_o_i,c(paste0("CLE.",AoI[1]-1,".M"),
# paste0("CLE.",AoI[2]-1,".M"),paste0("CLE.",AoI[1]-1,".F"),
# paste0("CLE.",AoI[2]-1,".F"),paste0("PLE.",AoI[1]-1,".M"),
# paste0("PLE.",AoI[2]-1,".M"),paste0("PLE.",AoI[1]-1,".F"),
# paste0("PLE.",AoI[2]-1,".F"))) # # lOutput_LE[[paste0("Summary_NL")]] <-
# data.frame(output_sum) # # openxlsx::write.xlsx(x=lOutput_LE,
# file=paste0(sOutputFolder, sTableName, "_Life_expectancy.xlsx"),
# rowNames=TRUE, colNames=TRUE) # saveRDS(object=lOutput_LE,
# file=paste0(sOutputFolder, sTableName, "_LE.rds", sep=""))
#   # 
#   # 
#   # ###############################################################################
#   # ### Simulations: Calculate Life Expectancy over the years for given ages
#   # ###############################################################################
#   # 
#   # 
#   # # Stochastic simulations
#   # if (iSims == 0){
#   #   LE_P_sim <- readRDS(file=paste0(sOutputFolder, "LE_P_years_sim.rds", sep=""))
#   #   LE_sim <- readRDS(file=paste0(sOutputFolder, "LE_years_sim.rds", sep=""))
#   #   iSimsA <- dim(LE_P_sim)[4]
#   #   iT_futA <- dim(LE_P_sim)[3]
#   # } else{
#   #   #Calculate period and cohort LE for future years - only for NL
#   #   #vY_NL2 <- rownames(objAG_BE_all$lQg_NL_tx$Male) #all years NL
#   #   vY_NL2 <- rownames(objAG_sim[[1]])
#   # 
#   #   LE_P_sim <- array(data=0, dim=c(2, length(AoI), sum(vY_NL2 >= proj_year), iSims))  # geslachten x ages of interest x prognosejaren x simulaties
#   #   LE_sim   <- array(data=0, dim=c(2, length(AoI), sum(vY_NL2 >= proj_year), iSims))  # geslachten x ages of interest x prognosejaren x simulaties
#   # 
#   #   for(i in 1:length(AoI)){
#   #     for(s in 1:iSims){
#   #       print(paste0("Levensverwachting ",s+iSims*(i-1)," van ",iSims*length(AoI)))
#   #       age <- AoI[i]
#   # 
#   #       # Males
#   #       LE_P_sim[1,i,,s] <- fnLifeExpPer(objAG_sim$Males[vY_NL2 >= proj_year,,s],age, age_ult = 120)
#   #       LE_sim[1,i,,s] <- fnLifeExpCoh(objAG_sim$Males[vY_NL2 >= proj_year,,s],age, age_ult = 120)
#   # 
#   #       # Females
#   #       LE_P_sim[2,i,,s] <- fnLifeExpPer(objAG_sim$Females[vY_NL2 >= proj_year,,s],age, age_ult = 120)
#   #       LE_sim[2,i,,s] <-fnLifeExpCoh(objAG_sim$Females[vY_NL2 >= proj_year,,s],age, age_ult = 120)
#   #       iSimsA <- iSims
#   #       iT_futA <- iT_fut
#   #     }
#   #   }
#   # 
#   # 
#   #   # Opslaan van levensverwachtingen obv simulaties
#   #   saveRDS(object=LE_P_sim, file=paste0(sOutputFolder, "LE_P_years_sim.rds", sep=""))
#   #   saveRDS(object=LE_sim, file=paste0(sOutputFolder, "LE_years_sim.rds", sep=""))
#   # 
#   # }
#   # 
#   # # Bepaal kwantielen van de levensverwachtingen
#   # 
#   # lOutput_LE_sim <- list();lOutput_LE_sim_ages <- list();
#   # for (i in 1:length(quantiles)) {
#   #   LE_P_sim_quantile <- apply(LE_P_sim, c(1,2,3), quantile, quantiles[i])
#   #   LE_sim_quantile <- apply(LE_sim, c(1,2,3), quantile, quantiles[i])
#   # 
#   #   # Maak een output dataframe
#   #   for (g in 1:2) {
#   #     for (j in 1:length(AoI)) {
#   #       gender <- ifelse(g==1,"Male","Female")
#   #       age <- AoI[j]
#   # 
#   #       if (g==1 & j==1) {
#   #         df_output_LE_P_sim <- as.data.frame(LE_P_sim_quantile[g,j,])
#   #         colnames(df_output_LE_P_sim) <- paste0("LE_period_",gender,", age_",age-1)
#   # 
#   #         df_output_LE_sim <- as.data.frame(LE_sim_quantile[g,j,])
#   #         colnames(df_output_LE_sim) <- paste0("LE_cohort_",gender,", age_",age-1)
#   # 
#   #       } else {
#   #         df_output_LE_P_sim <- cbind(df_output_LE_P_sim, as.data.frame(LE_P_sim_quantile[g,j,]))
#   #         colnames(df_output_LE_P_sim)[length(colnames(df_output_LE_P_sim))] <- paste0("LE period, ",gender,", age=",age-1)
#   # 
#   #         df_output_LE_sim <- cbind(df_output_LE_sim, as.data.frame(LE_sim_quantile[g,j,]))
#   #         colnames(df_output_LE_sim)[length(colnames(df_output_LE_sim))] <- paste0("LE cohort, ",gender,", age=",age-1)
#   # 
#   #       }
#   #     }
#   #   }
#   # 
#   #   # Creeren output set levensverwachtingen obv simulaties
#   #   lOutput_LE_sim[[paste0("LE_period_quant_",quantiles[i])]] <- df_output_LE_P_sim
#   #   lOutput_LE_sim[[paste0("LE_cohort_quant_",quantiles[i])]] <- df_output_LE_sim
#   # 
#   #   if(i==2){
#   #     df_output_LE_sim_q2 <- df_output_LE_sim
#   #   } else if (i==1){
#   #     df_output_LE_sim_q1 <- df_output_LE_sim
#   #   }
#   # }
#   # 
#   # openxlsx::write.xlsx(x=lOutput_LE_sim, file=paste0(sOutputFolder, "Life_expectancy_years_simulations.xlsx"), rowNames=TRUE, colNames=TRUE)
#   # 
#   # ###############################################################################
#   # ### GRAPH: Best estimate & confidence level
#   # ###############################################################################
#   # sY <- 1+eyNL
#   # eY <- sY + iT_futA -1
#   # png("Life Expectancy.png", height=600, width=600)
#   # type <- c(1,1)
#   # title <- paste0("Cohort life expectancy age 0 incl. quantiles: ",quantiles[2]," - ", quantiles[3], ". n=", iSimsA)
#   # plot(x=NULL, y=NULL,xlim = range(1,iT_futA),ylim = range(85,105), xlab="Year", ylab="Life Expectancy", main =title, xaxt = 'n')
#   # lines(df_output_LE_sim_q1[1],lty = 1, col = "blue")
#   # lines(df_output_LE_sim[1],lty = 2, col = "blue")
#   # lines(df_output_LE_sim_q2[1],lty = 2, col = "blue")
#   # lines(df_output_LE_sim_q1[3],lty = 1, col = "magenta")
#   # lines(df_output_LE_sim[3],lty = 2, col = "magenta")
#   # lines(df_output_LE_sim_q2[3],lty = 2, col = "magenta")
#   # axis(1, at=1:iT_futA, labels=sY:eY)
#   # legend("bottomright", legend = c("Male", "Female"), col = c("blue", "magenta"), lty = type )
#   # dev.off()
#   # 
#   # ###############################################################################
#   # ### Simulations: Calculate Life Expectancy over the ages for a given year
#   # ###############################################################################
#   # # Stochastic simulations
#   # if (iSims == 0){
#   #   LE_sim <- readRDS(file=paste0(sOutputFolder, "LE_ages_sim.rds", sep=""))
#   # } else{
#   #   #Calculate period and cohort LE for given years - only for NL
#   #   vY_NL2 <- as.numeric(rownames(objAG_BE_all$lQg_NL_tx$Male)) #all years NL
#   #   ages <- as.numeric(colnames(objAG_BE_all$lQg_NL_tx$Male))+1 #all ages, + 1 for purpose of column. Thus 1 = age 0
#   # 
#   #   LE_P_sim <- array(data=0, dim=c(2, length(ages), length(YoI), iSims))  # geslachten x ages of interest x prognosejaren x simulaties
#   #   LE_sim   <- array(data=0, dim=c(2, length(ages), length(YoI), iSims))  # geslachten x ages of interest x prognosejaren x simulaties
#   # 
#   #   for(i in 1:length(YoI)){
#   #     for(s in 1:iSims){
#   #       print(paste0("Levensverwachting ",s+iSims*(i-1)," van ",iSims*length(YoI)))
#   #       year <- YoI[i]
#   #       LE_sim[1,,i,s] <- fnLifeExpCoh(as.matrix(t(objAG_sim$Males[vY_NL2 >= YoI[i],,s])), 1, age_ult = 120)    # Males
#   #       LE_sim[2,,i,s] <- fnLifeExpCoh(as.matrix(t(objAG_sim$Females[vY_NL2 >= YoI[i],,s])), 1, age_ult = 120)  # Females
#   #     }
#   #   }
#   # 
#   #   # Opslaan van levensverwachtingen obv simulaties
#   #   saveRDS(object=LE_sim, file=paste0(sOutputFolder, "LE_ages_sim.rds", sep=""))
#   # }
#   # 
#   # # Bepaal afwijkende kwantielen van de levensverwachtingen (indien gewenst)
#   # 
#   # for (i in 1:length(quantiles)) {
#   #   LE_sim_quantile <- apply(LE_sim, c(1,2,3), quantile, quantiles[i])
#   #   # Maak een output dataframe
#   #   for (g in 1:2) {
#   #     for (j in 1:length(YoI)) {
#   #       gender <- ifelse(g==1,"Male","Female")
#   #       year <- YoI[j]
#   #       if (g==1 & j==1) {
#   #         df_output_LE_sim <- as.data.frame(LE_sim_quantile[g,,j])
#   #         colnames(df_output_LE_sim) <- paste0("LE cohort, ",gender,", year=",year)
#   #       } else {
#   #         df_output_LE_sim <- cbind(df_output_LE_sim, as.data.frame(LE_sim_quantile[g,,j]))
#   #         colnames(df_output_LE_sim)[length(colnames(df_output_LE_sim))] <- paste0("LE cohort, ",gender,", year=",year)
#   # 
#   #       }
#   #     }
#   #   }
#   # 
#   #   # Creeren output set levensverwachtingen obv simulaties
#   #   lOutput_LE_sim_ages[[paste0("LE_cohort_quant_",quantiles[i])]] <- df_output_LE_sim
#   # }
#   # 
#   # openxlsx::write.xlsx(x=lOutput_LE_sim_ages, file=paste0(sOutputFolder, "Life_expectancy_ages_simulations.xlsx"), rowNames=TRUE, colNames=TRUE)
#   # 
#   # 
#   # 
#   # 
