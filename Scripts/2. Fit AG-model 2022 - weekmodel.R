library(dplyr)

# The following code is written by Frank van Berkum, based on earlier code
# used within the working group Prognosetafels. This code is extended with 
# fnFitWeekmodel3.2 to model monthly seasonality. 

resize_data <- function(dfData, result_matrix){
  
  is_value_present <- ("2020" %in% colnames(dfData) | 
                       "2024" %in% colnames(dfData) |
                       "2025" %in% colnames(dfData))

  if (dim(result_matrix)[1]>1){
    if (is_value_present) {
      # first_column <- result_matrix[, 1]
      # new_first_column <- c(first_column[-1], NA)
      # result_matrix[, 1] <- new_first_column
      
      first_column <- result_matrix[, "2024"]
      new_first_column <- c(first_column[-1], NA)
      result_matrix[, "2024"] <- new_first_column
      
      first_column <- result_matrix[, "2025"]
      new_first_column <- c(first_column[-1], NA)
      result_matrix[, "2025"] <- new_first_column
      
      result_matrix <- result_matrix[-(54 * (1:ceiling(nrow(result_matrix) / 54))), ]
      
    }
  }
  return(result_matrix)
}

slice_weekly_data <- function(iAge_teshold_beta, dfData){

  result_matrix <- as.matrix(
    dfData %>%
      filter(age>=iAge_teshold_beta) %>% 
      dplyr::select(3:ncol(.))
  )
  result_matrix <- resize_data(dfData, result_matrix)
  #print(result_matrix)
  
  #print(result_matrix)
  return(result_matrix)
}

generate_kappa_matrix <- function(iYears, iWeeks){
  mKappa <- matrix(rep(2, iYears * iWeeks), ncol=iYears)     # vector iYears - average iyears  as starting point (creating for instance 3,2,1,0,-1,-2,-3)

  return(mKappa)
}


optim_frakkappa_multi_year <- function(mKappa, vBeta, dfExposures, dfDeaths){
  
  iYears <- dim(mKappa)[2]
  for (t in 1:(iYears)) {
    
    sYear <- colnames(dfExposures%>% dplyr::select(t+2))
    if (sYear %in% c("2020","2024", "2025")) iWeekCorrection<-0 else iWeekCorrection<-1
    for (w in 1:53){
      k1 <- mKappa[w,t];  
      k0 <- k1 - 1 # Starting condition
      while(abs(k1-k0) > 0.1) {  
        
        mExtw_tilde <- dfExposures %>% 
          filter(week==w-iWeekCorrection) %>% 
          dplyr::select(t+2)

        mDxtw_tilde <- dfDeaths %>% 			    
          filter(week==w-iWeekCorrection) %>% 
          dplyr::select(t+2)
        
        k0 <- k1;    
        k1 <- k0 - (sum( (exp(vBeta * k0) * vBeta * mExtw_tilde)) - sum(mDxtw_tilde * vBeta)) / 
          sum( (exp(vBeta * k0) * vBeta^2 * mExtw_tilde) ) 
        
        #print(mExtw_tilde)
        # cat(t, w, w-iWeekCorrection, sum( (mExtw_tilde) ), '\n')
        #cat(t, w, k1, k0, abs(k1-k0), abs(k1-k0) > 0.1, '\n')
        
      }
      #cat(w,t,k1, '\n')
      mKappa[w,t] <- k1
    }
  }
  return(mKappa)
}


optim_frakkappa_single_year <- function(mKappa, vBeta, dfExposures, dfDeaths){
  
  num_weeks <- max(dfExposures$week)
  if (num_weeks==53) iWeekCorrection<-0 else iWeekCorrection<-1
  
  for (w in (1:length(mKappa))){
    #cat(w, w-iWeekCorrection, '\n')
    k1 <- mKappa[w];  
    k0 <- k1 - 1 # Starting condition
    while(abs(k1-k0) > 0.1) {  
      mExtw_tilde <- dfExposures %>% 
        filter(week==w-iWeekCorrection) %>% 
        dplyr::select(3)
      
      mDxtw_tilde <- dfDeaths %>% 			    
        filter(week==w-iWeekCorrection) %>% 
        dplyr::select(3)
      k0 <- k1;    
      k1 <- k0 - (sum( (exp(vBeta * k0) * vBeta * mExtw_tilde)) - sum(mDxtw_tilde * vBeta)) / 
        sum( (exp(vBeta * k0) * vBeta^2 * mExtw_tilde) ) 
      
    }
    mKappa[w] <- k1
  }
  return(mKappa)
}


# A = Fit Lee Carter Process
fnFitWeekmodel <- function(dfData, iAge_teshold_beta, dRelTol=1E-6, iStep=2) {

	# iStep=1 : kappa[t1] = 0
	# iStep=2 : sum(kappa[t]) = 0
	iAges <- length(unique(dfData$age))	# number of ages
	iYears <- length(unique(dfData$year))	# number of years
	iWeeks <- length(unique(dfData[dfData$year == min(dfData$year), ]$week))	# number of weeks
	
	dfDeaths <- reshape2::dcast(dfData, age + week ~ year, value.var="dxt",fill = 0)
	dfExposures <- reshape2::dcast(dfData, age + week ~ year, value.var="ext",fill = 0)
	
	# Give initial estimates for vAlpha(x), vBeta(x) and vKappa(t) 
	vBeta  <- c(rep(0, iAge_teshold_beta), rep(0.04, iAges-iAge_teshold_beta))         # 1/number of ages is starting point (rep --> starting point for every age)
	const.beta	<- sum(vBeta)
	vBeta	<- vBeta / const.beta	
	
	mKappa <- generate_kappa_matrix(iYears, iWeeks)
	
	# Initialise iteration of maximum likelihood
	dL0 <- -1000000;  dL1 <- -999999;  iteration <- 0
 
	# We continue to iterate if the improvement in log-likelihood exceeds 10^(-18)
	while( abs(dL1-dL0) > dRelTol) {

		iteration <- iteration + 1;   dL0 <- dL1
    # Stage 1 optimise over vBeta(x)
		for (x in iAge_teshold_beta:(iAges-1)) {
			b1 <- vBeta[x+1];  
			b0 <- b1 - 1 # Starting condition
			while(abs(b1-b0) > 0.1) {
			  
			  mExtw_tilde <- dfExposures %>% 
			    filter(age==x) %>% 
			    dplyr::select(3:ncol(.))
			  mExtw_tilde <- resize_data(dfExposures, mExtw_tilde)
			  
			  mDxtw_tilde <- dfDeaths %>% 			    
			    filter(age==x) %>% 
			    dplyr::select(3:ncol(.))
			  mDxtw_tilde <- resize_data(dfDeaths, mDxtw_tilde)
			  
			  b0 <- b1;    
				b1 <- b0 - (sum( (exp(b0 * mKappa) * mKappa * mExtw_tilde) ) - sum(mDxtw_tilde*mKappa)) / 
				sum( (exp(b0 * mKappa) * mKappa^2 * mExtw_tilde) ) 
			}
			vBeta[x+1] <- b1
			
		}
    
		# Stage 2 optimise over vKappa(w,t)
		if (iYears==1) {
		  mKappa <- optim_frakkappa_single_year(mKappa, vBeta, dfExposures, dfDeaths)
		} else{
		  mKappa <- optim_frakkappa_multi_year(mKappa, vBeta, dfExposures, dfDeaths)
	
		}
	  
		# Now apply the constraints
		#mKappa[,1] <- c(mKappa[-1, 1], NA)
		#mKappa <- mKappa[-nrow(mKappa), ]
		#const.kappa	<- ifelse(iStep==1, mKappa[1], mean(mKappa_hat))
		# const.beta	<- sqrt(sum(vBeta^2))
		#const.beta	<- sqrt(sum(vBeta^2))
		const.beta	<- sum(vBeta)  
		const.kappa <- sign(mean(mKappa))
		#const.kappa <- 1
		
		mKappa	<- mKappa * const.kappa * const.beta      # Rescales fit
		vBeta	<- (vBeta/const.beta) / const.kappa 		                      # Ensures that betas  sum to 1
    #print(vBeta)
		
		mDxtw <- slice_weekly_data(iAge_teshold_beta, dfDeaths)
		mExtw <- slice_weekly_data(iAge_teshold_beta, dfExposures)
		#print(mExtw)
		
		mMu_hat <- c()
		for (x in (iAge_teshold_beta):(iAges-1)){
		  mMu_tilde <- exp(vBeta[x+1] * mKappa)
		  #print(vBeta[x+1])
		  mMu_hat <- rbind(mMu_hat, mMu_tilde)
		}
		
		log_ExMu <- log(mExtw * mMu_hat)
		#log_ExMu[!is.finite(log_ExMu)] <- 0
		
		ExMu <- mExtw * mMu_hat
		#ExMu[!is.finite(ExMu)] <- 0
		
		dL1 <- sum( mDxtw*log_ExMu - ExMu - lgamma(mDxtw+1) )
		#print(mKappa)
		#print(vBeta)
		cat(dL1, dL1-dL0, '\n')
		
	}		 # end while loop
	#dimnames(mMxt) <- list(Age=vX, Year=vY)
	colnames(mKappa) <- c(unique(dfData$year))
	list(frakbeta=vBeta, frakkappa=mKappa, mxt=mMu_hat) #output 
}
