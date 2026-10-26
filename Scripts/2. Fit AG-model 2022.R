

# The following code is written by Frank van Berkum, based on earlier code
# used within the working group Prognosetafels

set.seed(42)

# function to create forces of mortality from a, b, k
fnMxtFromPars <- function(alpha, beta, kappa, ages, years) {
  mxt <- exp( alpha %o% rep(1, length(kappa)) + beta %o% kappa)
	dimnames(mxt) <- list(Age=ages, Year=years)
	# mxt[,kappa == 0] <- 999
	return(mxt)
}

# A = Fit Lee Carter Process
# fnFitLC <- function(mDxt, mExt, vX, vY, dRelTol=1E-18, iStep=2) {
fnFitLC <- function(mDxt, mExt, w_C, vX, vY, dRelTol=1E-18, iStep=2) {
    
	# iStep=1 : kappa[t1] = 0
	# iStep=2 : sum(kappa[t]) = 0
	mMxt <- mDxt/mExt		# matrix of (raw) central death rates
	iAges <- length(vX)	# number of ages
	iYears <- length(vY)	# number of years

	# Give initial estimates for vAlpha(x), vBeta(x) and vKappa(t) 
	vAlpha <- rowMeans(log(mMxt))         # Row average of log(mMxt) as starting point
	vBeta  <- rep(1/iAges, iAges)         # 1/number of ages is starting point (rep --> starting point for every age)
	vKappa <- (iYears:1)-(iYears+1)/2     # vector iYears - average iyears  as starting point (creating for instance 3,2,1,0,-1,-2,-3)

	# Initialise iteration of maximum likelihood
	dL0 <- -1000000;  dL1 <- -999999;  iteration <- 0
	# We continue to iterate if the improvement in log-likelihood exceeds 10^(-18)
	while( abs(dL1-dL0) > dRelTol) {

		iteration <- iteration + 1;   dL0 <- dL1

		# Stage 1 optimise over vBeta(x)
		for (j in 1:iAges) {
			b1 <- vBeta[j];  
			b0 <- b1 - 1 # Starting condition
			while(abs(b1-b0) > 0.1) {
				b0 <- b1;  
				# b1 <- b0 - (sum( (exp(vAlpha[j] + b0 * vKappa) * vKappa * mExt[j, ]) ) - sum(mDxt[j, ]*vKappa)) / 
				# sum( (exp(vAlpha[j] + b0 * vKappa) * vKappa^2 * mExt[j, ]) ) 
				b1 <- b0 - (sum( (exp(vAlpha[j] + b0 * vKappa) * (vKappa*w_C) * mExt[j, ]) ) - sum(mDxt[j, ]*(vKappa*w_C))) / 
				  sum( (exp(vAlpha[j] + b0 * vKappa) * (vKappa^2) * w_C * mExt[j, ]) ) 
			}
			vBeta[j] <- b1
		}

		# Stage 2 optimise over vKappa(t)
		for (i in 1:iYears) {
			ifelse(w_C[i] > 0, {
		    k1 <- vKappa[i];  
			  k0 <- k1 - 1 # Starting condition
			  while(abs(k1-k0) > 0.1) {  
				  k0 <- k1; 
				  # k1 <- k0 - (sum( (exp(vAlpha + vBeta * k0) * vBeta * mExt[, i])) - sum(mDxt[, i] * vBeta)) / 
				  # sum( (exp(vAlpha + vBeta * k0) * vBeta^2 * mExt[, i]) ) 
				  k1 <- k0 - (sum( (exp(vAlpha + vBeta * k0) * vBeta * mExt[, i]*w_C[i])) - sum(mDxt[, i]*w_C[i] * vBeta)) / 
				    sum( (exp(vAlpha + vBeta * k0) * vBeta^2 * mExt[, i]*w_C[i]) )
			  }
			  vKappa[i] <- k1
			}, 
			{vKappa[i] <- 0})
		}

		# Now apply the constraints
    # const.kappa	<- ifelse(iStep==1, vKappa[1], mean(vKappa))
		const.kappa	<- ifelse(iStep == 1, vKappa[1], sum(vKappa) / sum(w_C > 0))
		# const.beta	<- sqrt(sum(vBeta^2))
		const.beta	<- sum(vBeta)
		vKappa	<- const.beta*(ifelse(w_C > 0, vKappa-const.kappa, 0))        # Ensures that kappas sum to 0
		# vKappa	<- const.beta*(vKappa-const.kappa)        # Ensures that kappas sum to 0
		vBeta	<- vBeta/const.beta		                      # Ensures that betas  sum to 1
		vAlpha	<- vAlpha+vBeta*const.beta*const.kappa    # Alpha needs to be adjusted to compensate

		for(j in 1:iAges) vAlpha[j] <- log(sum(mDxt[j, ]*w_C)) - log(sum(mExt[j, ]*exp(vBeta[j]*w_C*vKappa)))
		# for(j in 1:iAges) vAlpha[j] <- log(sum(mDxt[j, ])) - log(sum(mExt[j, ]*exp(vBeta[j]*vKappa)))
		for(i in 1:iYears) mMxt[, i] <- exp(vAlpha + vBeta * vKappa[i] * w_C[i])
		# for(i in 1:iYears) mMxt[, i] <- exp(vAlpha + vBeta * vKappa[i])

		# dL1 <- sum( mDxt*log(mExt*mMxt) - mExt*mMxt - lgamma(mDxt+1) )
		dL1 <- sum( t(t(mDxt)*w_C)*log(mExt*mMxt) - t(t(mExt)*w_C)*mMxt - lgamma(t(t(mDxt)*w_C)+1) )
	}		 # end while loop
	dimnames(mMxt) <- list(Age=vX, Year=vY)
	list(alpha=as.numeric(vAlpha), beta=vBeta, kappa=vKappa, ages=vX, years=vY, mxt=mMxt) #output 
}


# Houd deze functie op alleen schatten! Dat is handiger voor plotten van parameters!
# Output: parameters en gefitte intensities (EU, deltaNL en NL)
# fnFitAG2020 <- function(lData_EU, lData_NL, vX, vY_EU, vY_NL, parRest=2, dRelTol=1E-20) {
fnFitAG2020 <- function(lData_EU, lData_NL, w_EU, w_NL, vX, vY_EU, vY_NL, parRest=2, dRelTol=1E-20) {
	# Fit common trend for both genders
	objM1_EU <- list()
	for( gender in c("M", "F") ) {
		# objM1_EU[[gender]] <- fnFitLC(mDxt=lData_EU[[gender]]$Dxt, mExt=lData_EU[[gender]]$Ext, vX=vX, vY=vY_EU, dRelTol=dRelTol, iStep=parRest)
		objM1_EU[[gender]] <- fnFitLC(mDxt=lData_EU[[gender]]$Dxt, mExt=lData_EU[[gender]]$Ext, w_EU, vX=vX, vY=vY_EU, dRelTol=dRelTol, iStep=parRest)
	}
	# Extend kappa with ML kappa's if necessary
	for( gender in c("M", "F") ) {	
		if( max(vY_NL) > max(vY_EU) ) { 
			objM1_EU[[gender]]$kappa_old  <- objM1_EU[[gender]]$kappa
			dDrift <- mean(diff(objM1_EU[[gender]]$kappa[!(objM1_EU[[gender]]$kappa==0)]))
			# dDrift <- mean( diff( objM1_EU[[gender]]$kappa ) )
			objM1_EU[[gender]]$kappa <- c(objM1_EU[[gender]]$kappa, tail(objM1_EU[[gender]]$kappa,1) + dDrift * 1:(max(vY_NL)-max(vY_EU)) )
			names(objM1_EU[[gender]]$kappa) <- min(vY_EU):max(vY_NL)
			objM1_EU[[gender]]$years <- min(vY_EU):max(vY_NL)
		} else {
			objM1_EU[[gender]]$kappa_old <- objM1_EU[[gender]]$kappa
		}
	}
	vY_EU_old <- vY_EU
	vY_EU <- min(vY_EU):max(vY_NL)
	#fit country specific deviation
	objM1M1_NL <- list()
	# build to fit length EU and NL
	sY <- min(vY_NL) - min(vY_EU)+1
	sE <- length(objM1_EU[[gender]]$years)  # last element EY
	for( gender in c("M", "F") ) {
	  mxt <- fnMxtFromPars(alpha=objM1_EU[[gender]]$alpha, beta=objM1_EU[[gender]]$beta, kappa=objM1_EU[[gender]]$kappa, ages=objM1_EU[[gender]]$ages, years=objM1_EU[[gender]]$years)
		# mExt_adj <- lData_NL[[gender]]$Ext * mxt[ , as.character(vY_NL)]
	  mExt_adj <- lData_NL[[gender]]$Ext * mxt[, as.character(vY_NL)]
		# objM1M1_NL[[gender]] <- fnFitLC(mDxt=lData_NL[[gender]]$Dxt, mExt=mExt_adj, vX=vX, vY=vY_NL, dRelTol=dRelTol, iStep=2)
		objM1M1_NL[[gender]] <- fnFitLC(mDxt=lData_NL[[gender]]$Dxt, mExt=mExt_adj, w_NL, vX=vX, vY=vY_NL, dRelTol=dRelTol, iStep=2)
	}

	# Define output variables
	objAG <- list()
	for( gender in c("M", "F") ) {
		mDxt_EU = lData_EU[[gender]]$Dxt
		mExt_EU = lData_EU[[gender]]$Ext
		mDxt_NL = lData_NL[[gender]]$Dxt
		mExt_NL = lData_NL[[gender]]$Ext
		vX = vX
		vY_EU = vY_EU
		vY_NL = vY_NL
		vA = objM1_EU[[gender]]$alpha
		vB = objM1_EU[[gender]]$beta
		vK = objM1_EU[[gender]]$kappa
		vK_old = objM1_EU[[gender]]$kappa_old
		vAlpha = objM1M1_NL[[gender]]$alpha
		vBeta = objM1M1_NL[[gender]]$beta
		vKappa = objM1M1_NL[[gender]]$kappa
		mMxt_EU = fnMxtFromPars(alpha=objM1_EU[[gender]]$alpha, beta=objM1_EU[[gender]]$beta, kappa=objM1_EU[[gender]]$kappa, ages=vX, years=vY_EU)
		mMxt_EU_old = fnMxtFromPars(alpha=objM1_EU[[gender]]$alpha, beta=objM1_EU[[gender]]$beta, kappa=objM1_EU[[gender]]$kappa_old, ages=vX, years=vY_EU_old)
		mMxt_NL = mMxt_EU[ , as.character(vY_NL)] * fnMxtFromPars(alpha=objM1M1_NL[[gender]]$alpha, beta=objM1M1_NL[[gender]]$beta, kappa=objM1M1_NL[[gender]]$kappa, ages=vX, years=vY_NL)
    Dxt_EU_hat <- mExt_EU * mMxt_EU_old
		mDxt_NL_hat <- mExt_NL * mMxt_NL

		objAG[[gender]] <- list(
		  mDxt_EU = mDxt_EU, 
		  mExt_EU = mExt_EU,
		  mDxt_NL = mDxt_NL, 
		  mExt_NL = mExt_NL,
			vX = vX, 
			vY_EU = vY_EU, 
			vY_NL = vY_NL,
			vA = signif(vA, digits = 15),
			vB = signif(vB, digits = 15),
			vK = signif(vK, digits = 15),
			vK_old = vK_old,
			vAlpha = signif(vAlpha, digits = 15),
			vBeta = signif(vBeta, digits = 15),
			vKappa = signif(vKappa, digits = 15),
      mMxt_EU = mMxt_EU, 
			mMxt_EU_old = mMxt_EU_old, 
			mMxt_NL = mMxt_NL, 
			Dxt_EU_hat = Dxt_EU_hat, 
			mDxt_NL_hat = mDxt_NL_hat)
		}

	return(objAG)
}



# Houd deze functie op alleen schatten! Dat is handiger voor plotten van parameters!
# Output: parameters en gefitte intensities (EU, deltaNL en NL)
fnFitAG2020tbvDO <- function(lData_EU, lData_NL, vX, vY_EU, vY_NL, parRest=2, objCOVIDparameters, dRelTol=1E-20) {
  
  # Fit common trend for both genders
  objM1_EU <- list()
  for( gender in c("M", "F") ) {
    objM1_EU[[gender]] <- fnFitLC(mDxt=lData_EU[[gender]]$Dxt, mExt=lData_EU[[gender]]$Ext, vX=vX, vY=vY_EU, dRelTol=dRelTol, iStep=parRest)
  }
  
  # Extend kappa with ML kappa's if necessary
  for( gender in c("M", "F") ) {	
    # Lees parameters uit lijst zodat de notatie verderop korter wordt

    if( max(vY_NL) > max(vY_EU) ) { 
      objM1_EU[[gender]]$kappa_old  <- objM1_EU[[gender]]$kappa
      dDrift <- mean( diff( objM1_EU[[gender]]$kappa ) )
      objM1_EU[[gender]]$kappa <- c(objM1_EU[[gender]]$kappa, tail(objM1_EU[[gender]]$kappa,1) + dDrift * 1:(max(vY_NL)-max(vY_EU)) )
      objM1_EU[[gender]]$years <- min(vY_EU):max(vY_NL)
      
    } else {
      objM1_EU[[gender]]$kappa_old <- objM1_EU[[gender]]$kappa
    }
  }
  vY_EU_old <- vY_EU
  vY_EU <- min(vY_EU):max(vY_NL)
  
  #fit country specific deviation
  objM1M1_NL <- list()
  # build to fit length EU and NL
  sY <- min(vY_NL) - min(vY_EU)+1
  sE <- length(objM1_EU[[gender]]$years)  # last element EY
  for( gender in c("M", "F") ) {
    mxt <- fnMxtFromPars(alpha=objM1_EU[[gender]]$alpha, beta=objM1_EU[[gender]]$beta, kappa=objM1_EU[[gender]]$kappa, ages=objM1_EU[[gender]]$ages, years=objM1_EU[[gender]]$years)
    mExt_adj <- lData_NL[[gender]]$Ext * mxt[ , as.character(vY_NL)]
    objM1M1_NL[[gender]] <- fnFitLC(mDxt=lData_NL[[gender]]$Dxt, mExt=mExt_adj, vX=vX, vY=vY_NL, dRelTol=dRelTol, iStep=2)
  }
  
  # Define output variables
  objAG <- list()
  for( gender in c("M", "F") ) {
    COVIDjaren <- objCOVIDparameters[[gender]]$COVIDjaren
    frak_beta  <- objCOVIDparameters[[gender]]$frak_beta$X1
    frak_kappa <- objCOVIDparameters[[gender]]$frak_kappa$X1
    eta        <- objCOVIDparameters[[gender]]$eta$X1
    mDxt_EU = lData_EU[[gender]]$Dxt
    mExt_EU = lData_EU[[gender]]$Ext
    mDxt_NL = lData_NL[[gender]]$Dxt
    mExt_NL = lData_NL[[gender]]$Ext
    vX = vX
    vY_EU = vY_EU
    vY_NL = vY_NL
    vA = objM1_EU[[gender]]$alpha
    vB = objM1_EU[[gender]]$beta
    vK = objM1_EU[[gender]]$kappa
    vK_old = objM1_EU[[gender]]$kappa_old
    vAlpha = objM1M1_NL[[gender]]$alpha
    vBeta = objM1M1_NL[[gender]]$beta
    vKappa = objM1M1_NL[[gender]]$kappa
    mMxt_EU = fnMxtFromPars(alpha=objM1_EU[[gender]]$alpha, beta=objM1_EU[[gender]]$beta, kappa=objM1_EU[[gender]]$kappa, ages=vX, years=vY_EU)
    mMxt_EU_old = fnMxtFromPars(alpha=objM1_EU[[gender]]$alpha, beta=objM1_EU[[gender]]$beta, kappa=objM1_EU[[gender]]$kappa_old, ages=vX, years=vY_EU_old)
    mMxt_NL = mMxt_EU[ , as.character(vY_NL)] * fnMxtFromPars(alpha=objM1M1_NL[[gender]]$alpha, beta=objM1M1_NL[[gender]]$beta, kappa=objM1M1_NL[[gender]]$kappa, ages=vX, years=vY_NL)
    mDxt_EU_hat <- mExt_EU * mMxt_EU_old
    mDxt_NL_hat <- mExt_NL * mMxt_NL
   

      
    objAG[[gender]] <- list(mDxt_EU = mDxt_EU, mExt_EU = mExt_EU,
                            mDxt_NL = mDxt_NL, mExt_NL = mExt_NL,
                            vX = vX, vY_EU = vY_EU, vY_NL = vY_NL,
                            
                            vA = vA, vB = vB, vK = vK, vK_old = vK_old,
                            vAlpha = vAlpha, vBeta = vBeta, vKappa = vKappa,
                            
                            mMxt_EU = mMxt_EU, mMxt_EU_old = mMxt_EU_old, mMxt_NL = mMxt_NL, 
                            mDxt_EU_hat = mDxt_EU_hat, mDxt_NL_hat = mDxt_NL_hat )
  }
  
  return(objAG)
}



fnLogitCorr <- function(dRho) {
	# The below transformation is to ensure -1 < dRho < 1, 
	# and the output of this function is then -\infty to \infty
	dP <- (dRho+1) / 2
	return( log(dP / (1-dP)) )
}
fnInvLogitCorr <- function(dInv) {
	dP <- exp(dInv) / (exp(dInv) + 1)
	dRho <- 2*dP - 1
	return(dRho)
}

# Define the log likelihood function of the original AG2020 TS specification
# This function will later be used to optimize, so we use transformations of original parameters
fnLL_TS_AG2020_trans <- function(objAG, vPars) {

	vK_M <- objAG[["M"]]$vK
	vK_F <- objAG[["F"]]$vK
	vKappa_M <- objAG[["M"]]$vKappa
	vKappa_F <- objAG[["F"]]$vKappa

	# vPars consists of:
	# (theta[M], theta[F], a[M], a[F], c[M], c[F], 
	#  sigma[1], sigma[2], sigma[3], sigma[3],	
	#  rho[1,2], rho[1,3], rho[1,4], 
	#  rho[2,3], rho[2,4], rho[3,4] )

	dTheta_M <- vPars[1]
	dTheta_F <- vPars[2]
	dA_M <- fnInvLogitCorr(dInv=vPars[3])	
	dA_F <- fnInvLogitCorr(dInv=vPars[4])	
	dC_M <- vPars[5]	
	dC_F <- vPars[6]
	vSigma <- exp(vPars[7:10])
	mRho <- matrix(data=1, nrow=4, ncol=4)
	mRho[1,2] <- mRho[2,1] <- fnInvLogitCorr(dInv=vPars[11])
	mRho[1,3] <- mRho[3,1] <- fnInvLogitCorr(vPars[12])
	mRho[1,4] <- mRho[4,1] <- fnInvLogitCorr(vPars[13])
	mRho[2,3] <- mRho[3,2] <- fnInvLogitCorr(vPars[14])
	mRho[2,4] <- mRho[4,2] <- fnInvLogitCorr(vPars[15])
	mRho[3,4] <- mRho[4,3] <- fnInvLogitCorr(vPars[16])
	mSigma <- (vSigma %o% vSigma) * mRho
	if(sum(w_EU == 0) > 0) {
  	A2 <- diag(2)
  	A2_cumsum <- diag(2)
  	A4 <- diag(c(1, 1, dA_M, dA_F))
  	A4_cumsum <- diag(4)
  	A2_cumsum <- vSigma[1:2, 1:2]
  	A4_cumsum <- vSigma
  	for(i in 1:(count(w_EU[w_EU == 0])-1)) {
  	  A2_cumsum <- A2_cumsum + (A2 %^% i) %*% vSigma %*% t(A2 %^% i)
  	  A4_cumsum <- A4_cumsum + (A4 %^% i) %*% vSigma %*% t(A4 %^% i)
  	}
	}
browser()
	# Construct a matrix with all error terms
	# Note that the matrix has NA values at some places!
	vY_EU <- objAG[[1]]$vY_EU
	vY_NL <- objAG[[1]]$vY_NL

	vY_diff <- tail(vY_EU, -1)
	vY_p1 <- setdiff( tail(vY_EU, -1), tail(vY_NL, -1) )
	vY_p2 <- tail(vY_NL, -1) 

	mError <- matrix(data=NA, nrow=length(vY_EU)-1, ncol=4, dimnames=list(Year=tail(vY_EU,-1), TS=c("EU-M", "EU-F", "NL-M", "NL-F")))
	mError[as.character(vY_diff), "EU-M"] <- diff(vK_M) - dTheta_M
	mError[as.character(vY_diff), "EU-F"] <- diff(vK_F) - dTheta_F
	mError[as.character(vY_p2), "NL-M"] <- tail(vKappa_M, -1) - dA_M * head(vKappa_M, -1) - dC_M
	mError[as.character(vY_p2), "NL-F"] <- tail(vKappa_F, -1) - dA_F * head(vKappa_F, -1) - dC_F

  if(fixed_Sims) set.seed(5379)
	dLogLik <- sum( dmvnorm(x=mError[as.character(vY_p1), c("EU-M", "EU-F")], mean=rep(0,2), sigma=mSigma[1:2,1:2], log=TRUE)) + 
				sum( dmvnorm(x=mError[as.character(vY_p2), ], mean=rep(0,4), sigma=mSigma, log=TRUE))

	# Remember to return minus the loglikelihood: 
	# we want to maximize, but optim minimizes
	return( - dLogLik )
}
# Define the log likelihood function of the original AG2020 TS specification using the actual parameters
fnLL_TS_AG2020_orig <- function(objAG, objPars) {

	vK_M <- objAG[["M"]]$vK
	vK_F <- objAG[["F"]]$vK
	vKappa_M <- objAG[["M"]]$vKappa
	vKappa_F <- objAG[["F"]]$vKappa

	dTheta_M <- objPars$Drift[1]
	dTheta_F <- objPars$Drift[2]
	dA_M <- objPars$MR[1]
	dA_F <- objPars$MR[2]	

	mSigma <- objPars$Sigma

	# Construct a matrix with all error terms
	# Note that the matrix has NA values at some places!
	vY_EU <- objAG[[1]]$vY_EU
	vY_NL <- objAG[[1]]$vY_NL

	vY_diff <- tail(vY_EU, -1)
	vY_p1 <- setdiff( tail(vY_EU, -1), tail(vY_NL, -1) )
	vY_p2 <- tail(vY_NL, -1) 

	mError <- matrix(data=NA, nrow=length(vY_EU)-1, ncol=4, dimnames=list(Year=tail(vY_EU,-1), TS=c("EU-M", "EU-F", "NL-M", "NL-F")))
	mError[as.character(vY_diff), "EU-M"] <- diff(vK_M) - dTheta_M
	mError[as.character(vY_diff), "EU-F"] <- diff(vK_F) - dTheta_F
	mError[as.character(vY_p2), "NL-M"] <- tail(vKappa_M, -1) - dA_M * head(vKappa_M, -1) - dC_M
	mError[as.character(vY_p2), "NL-F"] <- tail(vKappa_F, -1) - dA_F * head(vKappa_F, -1) - dC_F


	dLogLik <- sum( dmvnorm(x=mError[as.character(vY_p1), c("EU-M", "EU-F")], mean=rep(0,2), sigma=mSigma[1:2,1:2], log=TRUE)) + 
				sum( dmvnorm(x=mError[as.character(vY_p2), ], mean=rep(0,4), sigma=mSigma, log=TRUE))

	cat("higher is better!\n"); flush.console()
	# Remember to return minus the loglikelihood: 
	# we want to maximize, but optim minimizes
	return( dLogLik )
}
# The parameters in the loglikehood function are often transformed versions
# Therefore, we define a function to obtain the parameters at the original level
# by applying the inverse functions of the transformations
# For correlation parameters: adjusted version of logit transform, see also functions above
# For variance parameters: log(sigma) is used in the function, so exp(par) results in sigma
fnReturnEstimatedPars_AG2020 <- function(vEstimatedPars) {
	mRho <- matrix(data=1, ncol=4, nrow=4)
	mRho[1,2] <- mRho[2,1] <- fnInvLogitCorr(vEstimatedPars[11])
	mRho[1,3] <- mRho[3,1] <- fnInvLogitCorr(vEstimatedPars[12])
	mRho[1,4] <- mRho[4,1] <- fnInvLogitCorr(vEstimatedPars[13])
	mRho[2,3] <- mRho[3,2] <- fnInvLogitCorr(vEstimatedPars[14])
	mRho[2,4] <- mRho[4,2] <- fnInvLogitCorr(vEstimatedPars[15])
	mRho[3,4] <- mRho[4,3] <- fnInvLogitCorr(vEstimatedPars[16])
	vSigma <- exp(vEstimatedPars[7:10])
	mSigma <- vSigma %o% vSigma * mRho
	list( 
	  Drift = signif(vEstimatedPars[c(1,2,5,6)], digits=15),
	  MR = signif(fnInvLogitCorr(vEstimatedPars[3:4]), digits=15),
	  Sigma = signif(mSigma, digits = 15),
	  Rho = signif(mRho, digits = 15)
	)
}
fnEstimateTS_AG2020 <- function(objAG, dRelTol=1E-18) {
	
	# vPars consists of (transformed versions of):
	# (theta[M], theta[F], a[M], a[F], c[M], c[F],
	#  sigma[1], sigma[2], sigma[3], sigma[4],	
	#  rho[1,2], rho[1,3], rho[1,4], 
	#  rho[2,3], rho[2,4], rho[3,4] )
  
  ### 
  if (fixed_Sims) set.seed(42)

	vStartPar <- c(rnorm(n=2,mean=-1,sd=0.25), fnLogitCorr(dRho=0.50), fnLogitCorr(dRho=0.5), rnorm(n=2, mean=0, sd=0.25),
				log(rnorm(n=4, mean=2, sd=0.25)),
				fnLogitCorr(0.15), fnLogitCorr(0.25), fnLogitCorr(0.35), 
				fnLogitCorr(0.45), fnLogitCorr(0.55), fnLogitCorr(0.65))
	objOptim <- optim(par=vStartPar, fn=fnLL_TS_AG2020_trans, objAG=objAG, method="BFGS", 
				control=list(reltol=dRelTol, fnscale=100, maxit=1E5)); 
	if(objOptim$convergence!=0) {
		stop("Optimization did not converge succesfully!"); flush.console()
	}
	lEstimatedPars <- fnReturnEstimatedPars_AG2020(as.numeric(objOptim$par))
	
	return(lEstimatedPars)
}

fnWritePars <- function(objAG, sFile) {
  
  # Number of ages and years
  iX <- length(objAG[[1]]$vA)
  iY <- length(objAG[[1]]$vY_NL)
  iYEU <- length(objAG[[1]]$vY_EU)
  # iYEU <- length(objAG[[1]]$vY_EU)
  
  # Age effects
  mAge_effects <- matrix(data=NA, nrow=iX, ncol=4*2)
  rownames(mAge_effects) <- 0:(iX-1)
  colnames(mAge_effects) <- c("A_M", "B_M", "alpha_M", "beta_M", "A_F", "B_F", "alpha_F", "beta_F")
  for(gender in 1:2) {
    mAge_effects[,1 + (gender-1)*4] <- objAG[[gender]]$vA
    mAge_effects[,2 + (gender-1)*4] <- objAG[[gender]]$vB
    mAge_effects[,3 + (gender-1)*4] <- objAG[[gender]]$vAlpha
    mAge_effects[,4 + (gender-1)*4] <- objAG[[gender]]$vBeta
  }
  # write.xlsx(x=mAge_effects, file=sFile, sheetName="Age effects", col.names=TRUE, row.names=TRUE)

  # Period effects
  mPeriod_effects <- matrix(data=NA, nrow=iYEU+max(0, max(objAG[[1]]$vY_NL)-max(objAG[[1]]$vY_EU)), ncol=4)
  # mPeriod_effects <- matrix(data=NA, nrow=iYEU
  rownames(mPeriod_effects) <- seq(from=min(objAG[[1]]$vY_EU), to=min(objAG[[1]]$vY_EU)+iYEU, by=1)
  colnames(mPeriod_effects) <- c("K_M", "kappa_M", "K_F", "kappa_F")
  for(gender in 1:2) {
    mPeriod_effects[,1 + (gender-1)*2] <- objAG[[gender]]$vK
    mPeriod_effects[(iYEU-iY+1):iYEU,2 + (gender-1)*2] <- objAG[[gender]]$vKappa
  }
  # write.xlsx(x=mPeriod_effects, file=sFile, sheetName="Period effects", col.names=TRUE, row.names=TRUE, append=TRUE)

  # Drift, constant mean reversion parameters
  mTS_pars <- matrix(data=NA, nrow=6, ncol=1)
  rownames(mTS_pars) <- c("theta_M", "theta_F","c_M", "c_F", "a_M", "a_F")
  objTS <- readRDS(here(run_folder, paste0(runname, "_filt_TS.rds")))
  # objTS <- fnEstimateTS_AG2020(objAG, dRelTol=1E-18)
  
  ### Check if parameters equal the parameters used in calculations.
  # objTS_to_check <- readRDS(here(run_folder, paste0(runname, "_TS.rds")))
  # print(identical(objTS, objTS_to_check))
  mTS_pars[1:4,] <- objTS$Drift
  mTS_pars[5:6,] <- objTS$MR
  # write.xlsx(x=mTS_pars, file=sFile, sheetName="TS pars", col.names=FALSE, row.names=TRUE, append=TRUE)
  
  # Covariance matrix
  mSigma <- objTS$Sigma
  mSigma_CSO <- mSigma
  mSigma_CSO[1,1] <- mSigma[1,1]
  mSigma_CSO[2,2] <- mSigma[3,3]
  mSigma_CSO[3,3] <- mSigma[2,2]
  mSigma_CSO[4,4] <- mSigma[4,4]
  mSigma_CSO[1,2] <- mSigma_CSO[2,1] <- mSigma[1,3]
  mSigma_CSO[1,3] <- mSigma_CSO[3,1] <- mSigma[1,2]
  mSigma_CSO[1,4] <- mSigma_CSO[4,1] <- mSigma[1,4]
  
  mSigma_CSO[2,3] <- mSigma_CSO[3,2] <- mSigma[2,3]
  mSigma_CSO[2,4] <- mSigma_CSO[4,2] <- mSigma[3,4]
  mSigma_CSO[3,4] <- mSigma_CSO[4,3] <- mSigma[2,4]
  
  colnames(mSigma) <- rownames(mSigma) <- c("eps_M", "eps_F", "delta_M", "delta_F")
  colnames(mSigma_CSO) <- rownames(mSigma_CSO) <- c("eps_M", "delta_M", "eps_F", "delta_F")
  # write.xlsx(x=mSigma, file=sFile, sheetName="Sigma", col.names=TRUE, row.names=TRUE, append=TRUE)
  # write.xlsx(x=mSigma_CSO, file=sFile, sheetName="Sigma_CSO", col.names=TRUE, row.names=TRUE, append=TRUE)
  
  # Cholesky decomposition of covariance matrix
  mChol <- chol(mSigma_CSO)
  colnames(mChol) <- rownames(mChol) <- c("eps_M", "delta_M", "eps_F", "delta_F")
  # write.xlsx(x=mChol, file=sFile, sheetName="Cholesky_CSO", col.names=TRUE, row.names=TRUE, append=TRUE)
  
  list_output <- list(mAge_effects, mPeriod_effects, mTS_pars,mSigma_CSO,mChol)
  names(list_output) <- cbind("Age effects","Period effects","TS pars","Sigma_CSO","Cholesky_CSO")
  
  saveRDS(object=objTS, file=here(run_folder, paste0(runname, "filt_TS.rds")))
  openxlsx::write.xlsx(x=list_output, file=sFile, rowNames=TRUE, colNames=TRUE)
  
}

fnWriteBE <- function(objAG_BE, sFile) {
	
	# Write best estimate projection to spreadsheet (including fit on historical data!)

  # Write best estimate projection to spreadsheet (including fit on historical data!)
  list_output <- list()
  for(gender in c("Male", "Female")) {
    # Output NL
    list_output[[paste0("NL", gender)]] <- t(objAG_BE$lQg_NL_tx[[gender]])
  }
  write_xlsx(x=list_output, file=sFile, rowNames=TRUE, colNames=TRUE)
  
  
  
}
fnWriteBE_all <- function(objAG_BE, sFile) {
	
	# Write best estimate projection to spreadsheet (including fit on historical data!)
  list_output <- list()
  for(gender in c("Male", "Female")) {
    # Output NL
    list_output[[paste0("NL", gender)]] <- t(objAG_BE$lQg_NL_tx[[gender]])
    # Output EU
    list_output[[paste0("EU", gender)]] <- t(objAG_BE$lQg_EU_tx[[gender]])
  }
  write_xlsx(x=list_output, file=sFile, rowNames=TRUE, colNames=TRUE)
 
}

# Generate best estimate mortality forecast
fnProjectAG_BE_kannisto <- function(objAG, iT_fut=20) {
  
  # Estimate the TS model and extract regression results
  objTS <- fnEstimateTS_AG2020(objAG=objAG, dRelTol=1E-18)  
  mCov  <- objTS$Sigma
  mRho  <- objTS$Rho
  vDrift <- objTS$Drift
  vMR    <- c(1, 1, objTS$MR)
  
  # Simulate the time series
  vY_NL <- objAG[[1]]$vY_NL
  iT_his <- length(vY_NL)		# length of the historical time series
  mTS_BE <- matrix(data=0, nrow=4, ncol=iT_his + iT_fut, dimnames=list(TS=c("EU-M", "EU-F", "NL-M", "NL-F"), Year=seq(from=vY_NL[1], by=1, length.out=iT_his+iT_fut)))
  
  # First populate historical kappas
  mTS_BE["EU-M", as.character(vY_NL)] <- tail(objAG[["M"]]$vK, iT_his)
  mTS_BE["EU-F", as.character(vY_NL)] <- tail(objAG[["F"]]$vK, iT_his)
  mTS_BE["NL-M", as.character(vY_NL)] <- tail(objAG[["M"]]$vKappa, iT_his)
  mTS_BE["NL-F", as.character(vY_NL)] <- tail(objAG[["F"]]$vKappa, iT_his)
  
  # Now create forecasts, but loop over the years instead of over simualations
  for(i in 1:iT_fut) {
    mTS_BE[ , iT_his + i] <- vDrift + vMR * mTS_BE[ , iT_his + i - 1] 
  }
  
  # Reconstruct mortality rates
  iX_compl <- length(0:120) # length(objAG[[1]]$vX)	
  iX_cal <- length(objAG[[1]]$vX)	
  iX_kan <- iX_compl - iX_cal
  vX <- objAG[[1]]$vX
  vOnes <- rep(1, iT_his+iT_fut)
  lMg_tx <- list()
  lQg_tx <- list()
  vG <- c("Male", "Female")
  for(g in 1:2) {
    # Extract parameter estimates to enable shorter notation below
    vA <- objAG[[g]]$vA
    vB <- objAG[[g]]$vB
    vAlpha <- objAG[[g]]$vAlpha
    vBeta <- objAG[[g]]$vBeta
    # Define output object
    lMg_tx[[vG[g]]] <- matrix(data=0, nrow=iT_his+iT_fut, ncol=iX_compl)
    # Loop over ages and years
    for(i in 1:(iT_his+iT_fut)) {
      for(j in 1:iX_cal) {
        lMg_tx[[vG[g]]][i, j] <- exp( vA[j] + vB[j] * mTS_BE[g, i] + vAlpha[j] + vBeta[j] * mTS_BE[2+g, i] )
      }
    }
    # Close mortality table		
    lMg_tx[[vG[g]]][ , iX_cal + 1:iX_kan] <- ApplyKannisto_mtx(mtx=lMg_tx[[vG[g]]][ , 1:iX_cal], ages=vX, kannisto_nages=iX_kan, kannisto_nobs=11)
    lQg_tx[[vG[g]]] <- 1 - exp( - lMg_tx[[vG[g]]] )
    dimnames(lQg_tx[[vG[g]]]) <- list(Year=seq(from=vY_NL[1], by=1, length.out=iT_his+iT_fut), Age=0:120)
  }
  
  return(lQg_tx)
  
}

fnProjectAG_BE_all_kannisto <- function(objAG, iT_fut=20) {
  
  # Estimate the TS model and extract regression results
  objTS <- fnEstimateTS_AG2020(objAG=objAG, dRelTol=1E-18)  
  mCov  <- objTS$Sigma
  mRho  <- objTS$Rho
  vDrift <- objTS$Drift
  vMR    <- c(1, 1, objTS$MR)
  
  # Simulate the time series
  vY_NL <- objAG[[1]]$vY_NL
  vY_EU <- objAG[[1]]$vY_EU
  iT_his <- length(vY_EU)		# length of the historical time series
  iT_his_diff <- length(vY_EU)-length(vY_NL)+1
  mTS_BE <- matrix(data=0, nrow=4, ncol=iT_his + iT_fut, dimnames=list(TS=c("EU-M", "EU-F", "NL-M", "NL-F"), Year=seq(from=vY_EU[1], by=1, length.out=iT_his+iT_fut)))
  
  # First populate historical kappas
  mTS_BE["EU-M", as.character(vY_EU)] <- tail(objAG[["M"]]$vK, iT_his)
  mTS_BE["EU-F", as.character(vY_EU)] <- tail(objAG[["F"]]$vK, iT_his)
  mTS_BE["NL-M", as.character(vY_NL)] <- tail(objAG[["M"]]$vKappa, iT_his)
  mTS_BE["NL-F", as.character(vY_NL)] <- tail(objAG[["F"]]$vKappa, iT_his)
  
  # Now create forecasts, but loop over the years instead of over simulations
  for(i in 1:iT_fut) {
    mTS_BE[ , iT_his + i] <- vDrift + vMR * mTS_BE[ , iT_his + i - 1] 
  }
  
  # Reconstruct mortality rates
  iX_compl <- length(0:120) # length(objAG[[1]]$vX)	
  iX_cal <- length(objAG[[1]]$vX)	
  iX_kan <- iX_compl - iX_cal
  vX <- objAG[[1]]$vX
  vOnes <- rep(1, iT_his+iT_fut)
  lQg_EU_tx <- lQg_NL_tx <- list()
  vG <- c("Male", "Female")
  
  for(g in 1:2) {
    # Extract parameter estimates to enable shorter notation below
    vA <- objAG[[g]]$vA
    vB <- objAG[[g]]$vB
    vAlpha <- objAG[[g]]$vAlpha
    vBeta <- objAG[[g]]$vBeta
    # Define output object
    mMu_EU_tx <- mMu_NL_tx <- matrix(data=0, nrow=iT_his+iT_fut, ncol=iX_compl)
    # Loop over ages and years
    for(i in 1:(iT_his+iT_fut)) {
      for(j in 1:iX_cal) {
        mMu_EU_tx[i, j] <- exp( vA[j] + vB[j] * mTS_BE[g, i] )
        if(i >= iT_his_diff){
          mMu_NL_tx[i, j] <- exp( vA[j] + vB[j] * mTS_BE[g, i] + vAlpha[j] + vBeta[j] * mTS_BE[2+g, i] )
        }
      }
    }
    
    # Close mortality table		
    mMu_EU_tx[ , iX_cal + 1:iX_kan] <- ApplyKannisto_mtx(mtx=mMu_EU_tx[ , 1:iX_cal], ages=vX, kannisto_nages=iX_kan, kannisto_nobs=11)
    mMu_NL_tx <- mMu_NL_tx[iT_his_diff :nrow(mMu_NL_tx),]
    mMu_NL_tx[ , iX_cal + 1:iX_kan] <- ApplyKannisto_mtx(mtx=mMu_NL_tx[ , 1:iX_cal], ages=vX, kannisto_nages=iX_kan, kannisto_nobs=11)
    lQg_EU_tx[[vG[g]]] <- 1 - exp( - mMu_EU_tx )
    lQg_NL_tx[[vG[g]]] <- 1 - exp( - mMu_NL_tx )
    dimnames(lQg_EU_tx[[vG[g]]]) <- list(Year=seq(from=vY_EU[1], by=1, length.out=iT_his+iT_fut), Age=0:120)
    dimnames(lQg_NL_tx[[vG[g]]]) <- list(Year=seq(from=vY_NL[1], by=1, length.out=iT_his+iT_fut-iT_his_diff+1), Age=0:120)
  }
  
  list(lQg_NL_tx=lQg_NL_tx, lQg_EU_tx=lQg_EU_tx)
}


# fnProjectAG_BE_alternatiefB_HellingEnConstantetbvDO <- function(objAG, iT_fut=20,output,gender) {
#   
#   # Estimate the TS model and extract regression results
#   objTS <- fnEstimateTS_AG2020(objAG=objAG, dRelTol=1E-18)  
#   mCov  <- objTS$Sigma
#   mRho  <- objTS$Rho
#   vDrift <- objTS$Drift
#   vMR    <- c(1, 1, objTS$MR)
#   
#   # Simulate the time series
#   vY_NL <- objAG[[1]]$vY_NL
#   iT_his <- length(vY_NL)		# length of the historical time series
#   iT_his_diff <- length(vY_EU)-length(vY_NL)+1
#   mTS_BE <- matrix(data=0, nrow=4, ncol=iT_his + iT_fut, dimnames=list(TS=c("EU-M", "EU-F", "NL-M", "NL-F"), Year=seq(from=vY_NL[1], by=1, length.out=iT_his+iT_fut)))
#   
#   # First populate historical kappas
#   mTS_BE["EU-M", as.character(vY_NL)] <- tail(objAG[["M"]]$vK, iT_his)
#   mTS_BE["EU-F", as.character(vY_NL)] <- tail(objAG[["F"]]$vK, iT_his)
#   mTS_BE["NL-M", as.character(vY_NL)] <- tail(objAG[["M"]]$vKappa, iT_his)
#   mTS_BE["NL-F", as.character(vY_NL)] <- tail(objAG[["F"]]$vKappa, iT_his)
#   
#   # Now create forecasts, but loop over the years instead of over simualations
#   for(i in 1:iT_fut) {
#     mTS_BE[ , iT_his + i] <- vDrift + vMR * mTS_BE[ , iT_his + i - 1] 
#   }
#   
#   # Reconstruct mortality rates
#   iX_compl <- length(0:120) # length(objAG[[1]]$vX)	
#   iX_cal <- length(objAG[[1]]$vX)	
#   iX_kan <- iX_compl - iX_cal
#   vX <- objAG[[1]]$vX
#   vOnes <- rep(1, iT_his+iT_fut)
#   lQg_NL_tx <- list()
#   vG <- c("Male", "Female")
#   if(gender == "M"){
#     g <- 1
#   } else g<-2
#   
#     # Extract parameter estimates to enable shorter notation below
#     vA <- objAG[[g]]$vA
#     vB <- objAG[[g]]$vB
#     vAlpha <- objAG[[g]]$vAlpha
#     vBeta <- objAG[[g]]$vBeta
#     
#     # Bepaal sterfte intensiteiten t/m leeftijd 90 door geschatte paramters toe te passen voor EU en NL
#     # Define output object
#     mMu_EU_tx <- mMu_NL_tx <- matrix(data=0, nrow=iT_his+iT_fut, ncol=iX_compl)
#     # Loop over ages and years
#     for(i in 1:(iT_his+iT_fut)) {
#       for(j in 1:iX_cal) {
#         mMu_EU_tx[i, j] <- exp( vA[j] + vB[j] * mTS_BE[g, i] )
#         if(i >= iT_his_diff){
#           mMu_NL_tx[i, j] <- exp( vA[j] + vB[j] * mTS_BE[g, i] + vAlpha[j] + vBeta[j] * mTS_BE[2+g, i] )
#         }
#       }
#   }
# 
# if(output == "EU"){
#   return(mMu_EU_tx)
# } else return(mMu_NL_tx)
#   
# }

# Generate best estimate mortality forecast
fnProjectAG_BE_alternatiefB_HellingEnConstante <- function(objAG, iT_fut=20) {
  
  # Estimate the TS model and extract regression results
  objTS <- fnEstimateTS_AG2020(objAG=objAG, dRelTol=1E-18)  
  mCov  <- objTS$Sigma
  mRho  <- objTS$Rho
  vDrift <- objTS$Drift
  vMR    <- c(1, 1, objTS$MR)
  
  # Simulate the time series
  vY_NL <- objAG[[1]]$vY_NL
  iT_his <- length(vY_NL)		# length of the historical time series
  iT_his_diff <- length(vY_EU)-length(vY_NL)+1
  mTS_BE <- matrix(data=0, nrow=4, ncol=iT_his + iT_fut, dimnames=list(TS=c("EU-M", "EU-F", "NL-M", "NL-F"), Year=seq(from=vY_NL[1], by=1, length.out=iT_his+iT_fut)))
  
  # First populate historical kappas
  mTS_BE["EU-M", as.character(vY_NL)] <- tail(objAG[["M"]]$vK, iT_his)
  mTS_BE["EU-F", as.character(vY_NL)] <- tail(objAG[["F"]]$vK, iT_his)
  mTS_BE["NL-M", as.character(vY_NL)] <- tail(objAG[["M"]]$vKappa, iT_his)
  mTS_BE["NL-F", as.character(vY_NL)] <- tail(objAG[["F"]]$vKappa, iT_his)
  
  # Now create forecasts, but loop over the years instead of over simualations
  for(i in 1:iT_fut) {
    mTS_BE[ , iT_his + i] <- vDrift + vMR * mTS_BE[ , iT_his + i - 1] 
  }
  
  # Reconstruct mortality rates
  iX_compl <- length(0:120) # length(objAG[[1]]$vX)	
  iX_cal <- length(objAG[[1]]$vX)	
  iX_kan <- iX_compl - iX_cal
  vX <- objAG[[1]]$vX
  vOnes <- rep(1, iT_his+iT_fut)
  lQg_NL_tx <- list()
  vG <- c("Male", "Female")
  for(g in 1:2) {
    # Extract parameter estimates to enable shorter notation below
    vA <- objAG[[g]]$vA
    vB <- objAG[[g]]$vB
    vAlpha <- objAG[[g]]$vAlpha
    vBeta <- objAG[[g]]$vBeta
    
    # Bepaal sterfte intensiteiten t/m leeftijd 90 door geschatte paramters toe te passen voor EU en NL
    # Define output object
    mMu_EU_tx <- mMu_NL_tx <- matrix(data=0, nrow=iT_his+iT_fut, ncol=iX_compl)
    # Loop over ages and years
    for(i in 1:(iT_his+iT_fut)) {
      for(j in 1:iX_cal) {
        mMu_EU_tx[i, j] <- exp( vA[j] + vB[j] * mTS_BE[g, i] )
        if(i >= iT_his_diff){
          mMu_NL_tx[i, j] <- exp( vA[j] + vB[j] * mTS_BE[g, i] + vAlpha[j] + vBeta[j] * mTS_BE[2+g, i] )
        }
      }
    }
    
    
    ####################################################################
    # Alternatief B: Extrapoleer A, B, alpha en beta (lineair of m.b.v. Kannisto uitkomsten)
    kannisto_nages <- 30
    kannisto_nobs <- 11
    x_kannisto <- (iX_cal):(tail(vX,1)+kannisto_nages)
    
    ####################################################################
    # Europa
    # Extrapoleer B door parameters vanuit lineaire regressie over de kannisto observaties (lft 80 t/m 90) toe te passen
    vB_extrap <- log(vB)
    B_lm <- lm(log(tail(vB,kannisto_nobs)) ~ tail(vX,kannisto_nobs))
    for (i in x_kannisto) { vB_extrap[i+1] <- B_lm$coefficients[1] + B_lm$coefficients[2]*i }
    vB_extrap <- exp(vB_extrap)
    
    # Pas eerst "gewoon" Kannisto toe over de geschatte Europese sterfte intensiteiten
    mMu_EU_tx[ , iX_cal + 1:iX_kan] <- ApplyKannisto_mtx(mtx=mMu_EU_tx[ , 1:iX_cal], ages=vX, kannisto_nages=iX_kan, kannisto_nobs=11)
    # Stel A voor leeftijden 91 t/m 120 gelijk aan de log sterfte intensiteit van de laatste steekproefperiode
    vA_extrap <- rep(NA,iX_cal+iX_kan)
    vA_extrap[vX+1] <- vA
    vA_extrap[x_kannisto+1] <- log(mMu_EU_tx[length(vY_EU), x_kannisto+1]) - vB_extrap[x_kannisto+1] * mTS_BE[g,length(vY_EU)]
    vA_extrap <- pmin(0,vA_extrap) #alleen nodig voor leeftijden > 120
    
    ####################################################################
    # Nederland
    # Extrapoleer alpha lineair naar 0
    vAlpha_extrap <- vAlpha
    for (i in x_kannisto) { vAlpha_extrap[i+1] <- vAlpha_extrap[i] - (tail(vAlpha,1) / (kannisto_nages)) }

    # Pas eerst "gewoon" Kannisto toe over de geschatte Europese sterfte intensiteiten
    mMu_NL_tx[(iT_his_diff:nrow(mMu_NL_tx)), (iX_cal+1:iX_kan)] <- ApplyKannisto_mtx(mtx=mMu_NL_tx[(iT_his_diff:nrow(mMu_NL_tx)), 1:iX_cal], ages=vX, kannisto_nages=iX_kan, kannisto_nobs=11)
    # Bepaal het verschil tussen de NL en EU sterfte intensiteiten
    mMu_NL_tx_diff <- mMu_NL_tx / mMu_EU_tx[(dim(mMu_EU_tx)[1]-dim(mMu_NL_tx)[1]+1):(dim(mMu_EU_tx)[1]),]
    # Bepaal beta door voor de laatste steekproefperiode de geextrapoleerde alpha af te halen van log(Mu) en te delen door kappa
    vBeta_extrap <- (log(mMu_NL_tx_diff[length(vY_EU),]) - vAlpha_extrap)/mTS_BE[2+g,length(vY_EU)]
    
    ####################################################################
    # Bepaal de sterfte intensiteiten voor leeftijden 91 t/m 120 o.b.v. de doorgetrokken A, B, alpha en beta (_extrap) voor alle jaren 
    # Loop over ages and years
    for(i in 1:(iT_his+iT_fut)) {
      for(j in x_kannisto+1) {
        mMu_EU_tx[i, j] <- exp( vA_extrap[j] + vB_extrap[j] * mTS_BE[g, i] )
        if(i >= iT_his_diff){
          mMu_NL_tx[i, j] <- exp( vA_extrap[j] + vB_extrap[j] * mTS_BE[g, i] + vAlpha_extrap[j] + vBeta_extrap[j] * mTS_BE[2+g, i] )
        }
      }
    }
    
    # Kort de Nederlandse sterfte intensiteiten in (starten bij 1983)
    mMu_NL_tx <- mMu_NL_tx[iT_his_diff:nrow(mMu_NL_tx),]
    
    ####################################################################
    
    #Translate mtx to qtx
    lQg_NL_tx[[vG[g]]] <- signif(1 - exp( - mMu_NL_tx ), digits = 15)
    dimnames(lQg_NL_tx[[vG[g]]]) <- list(Year=seq(from=vY_NL[1], by=1, length.out=iT_his+iT_fut-iT_his_diff+1), Age=0:120)
    
    
  }
  
  return(lQg_NL_tx)
  
}

fnProjectAG_BE_all_alternatiefB_HellingEnConstante <- function(objAG, iT_fut=20) {

  # Estimate the TS model and extract regression results
  # objTS <- fnEstimateTS_AG2020(objAG=objAG, dRelTol=1E-18)  
  # saveRDS(objTS, here(run_folder, paste0(runname, "_TS.rds")))
  objTS <- readRDS(here(run_folder, paste0(runname, "_filt_TS.rds")))
  mCov  <- objTS$Sigma
  mRho  <- objTS$Rho
  vDrift <- objTS$Drift
  vMR    <- c(1, 1, objTS$MR)
  
  # Simulate the time series
  vY_NL <- objAG[[1]]$vY_NL
  vY_EU <- objAG[[1]]$vY_EU
  iT_his <- length(vY_EU)		# length of the historical time series
  iT_his_diff <- length(vY_EU)-length(vY_NL)+1
  mTS_BE <- matrix(data=0, nrow=4, ncol=iT_his + iT_fut, dimnames=list(TS=c("EU-M", "EU-F", "NL-M", "NL-F"), Year=seq(from=vY_EU[1], by=1, length.out=iT_his+iT_fut)))
  
  # First populate historical kappas
  mTS_BE["EU-M", as.character(vY_EU)] <- tail(objAG[["M"]]$vK, iT_his)
  mTS_BE["EU-F", as.character(vY_EU)] <- tail(objAG[["F"]]$vK, iT_his)
  mTS_BE["NL-M", as.character(vY_NL)] <- tail(objAG[["M"]]$vKappa, iT_his)
  mTS_BE["NL-F", as.character(vY_NL)] <- tail(objAG[["F"]]$vKappa, iT_his)
  
  # Now create forecasts, but loop over the years instead of over simulations
  for(i in 1:iT_fut) {
    mTS_BE[ , iT_his + i] <- vDrift + vMR * mTS_BE[ , iT_his + i - 1] 
  }
  
  # Reconstruct mortality rates
  iX_compl <- length(0:120) # length(objAG[[1]]$vX)	
  iX_cal <- length(objAG[[1]]$vX)	
  iX_kan <- iX_compl - iX_cal
  vX <- objAG[[1]]$vX
  vOnes <- rep(1, iT_his+iT_fut)
  lQg_EU_tx <- lQg_NL_tx <- list()
  vG <- c("Male", "Female")
  
  for(g in 1:2) {
    # Extract parameter estimates to enable shorter notation below
    vA <- objAG[[g]]$vA
    vB <- objAG[[g]]$vB
    vAlpha <- objAG[[g]]$vAlpha
    vBeta <- objAG[[g]]$vBeta
    
    # Bepaal sterfte intensiteiten t/m leeftijd 90 door geschatte paramters toe te passen voor EU en NL
    # Define output object
    mMu_EU_tx <- mMu_NL_tx <- matrix(data=0, nrow=iT_his+iT_fut, ncol=iX_compl)
    # Loop over ages and years
    for(i in 1:(iT_his+iT_fut)) {
      for(j in 1:iX_cal) {
        mMu_EU_tx[i, j] <- exp( vA[j] + vB[j] * mTS_BE[g, i] )
        if(i >= iT_his_diff){
          mMu_NL_tx[i, j] <- exp( vA[j] + vB[j] * mTS_BE[g, i] + vAlpha[j] + vBeta[j] * mTS_BE[2+g, i] )
        }
      }
    }
    
    ####################################################################
    # Alternatief B: Extrapoleer A, B, alpha en beta (lineair of m.b.v. Kannisto uitkomsten)
    kannisto_nages <- 30
    kannisto_nobs <- 11
    x_kannisto <- (iX_cal):(tail(vX,1)+kannisto_nages)
    
    ####################################################################
    # Europa
    # Extrapoleer B door parameters vanuit lineaire regressie over de kannisto observaties (lft 80 t/m 90) toe te passen
    vB_extrap <- log(vB)
    B_lm <- lm(log(tail(vB,kannisto_nobs)) ~ tail(vX,kannisto_nobs))
    for (i in x_kannisto) { vB_extrap[i+1] <- B_lm$coefficients[1] + B_lm$coefficients[2]*i }
    vB_extrap <- exp(vB_extrap)
	  vB_extrap[1:91] <- vB #niet per se nodig. heb gecheckt of er niets veranderd op 15 significante cijfers bij log en exp transformatie
    
    # Pas eerst "gewoon" Kannisto toe over de geschatte Europese sterfte intensiteiten
    mMu_EU_tx[ , iX_cal + 1:iX_kan] <- ApplyKannisto_mtx(mtx=mMu_EU_tx[ , 1:iX_cal], ages=vX, kannisto_nages=iX_kan, kannisto_nobs=11)
    # Stel A voor leeftijden 91 t/m 120 gelijk aan de log sterfte intensiteit van de laatste steekproefperiode
    vA_extrap <- rep(NA,iX_cal+iX_kan)
    vA_extrap[vX+1] <- vA
    vA_extrap[x_kannisto+1] <- log(mMu_EU_tx[length(vY_EU), x_kannisto+1]) - vB_extrap[x_kannisto+1] * mTS_BE[g,length(vY_EU)]
    vA_extrap <- pmin(0,vA_extrap) #alleen nodig voor leeftijden > 120
    
    ####################################################################
    # Nederland
    # Extrapoleer alpha lineair naar 0
    vAlpha_extrap <- vAlpha
    for (i in x_kannisto) { vAlpha_extrap[i+1] <- vAlpha_extrap[i] - (tail(vAlpha,1) / (kannisto_nages)) }

    # Pas eerst "gewoon" Kannisto toe over de geschatte Europese sterfte intensiteiten
    mMu_NL_tx[(iT_his_diff:nrow(mMu_NL_tx)), (iX_cal+1:iX_kan)] <- ApplyKannisto_mtx(mtx=mMu_NL_tx[(iT_his_diff:nrow(mMu_NL_tx)), 1:iX_cal], ages=vX, kannisto_nages=iX_kan, kannisto_nobs=11)
    # Bepaal het verschil tussen de NL en EU sterfte intensiteiten
    mMu_NL_tx_diff <- mMu_NL_tx / mMu_EU_tx[(dim(mMu_EU_tx)[1]-dim(mMu_NL_tx)[1]+1):(dim(mMu_EU_tx)[1]),]
    # Bepaal beta door voor de laatste steekproefperiode de geextrapoleerde alpha af te halen van log(Mu) en te delen door kappa

	  vBeta_extrap <- vBeta
    
    vBeta_extrap[x_kannisto + 1] <-
      (log(mMu_NL_tx_diff[length(vY_EU), x_kannisto + 1]) -
         vAlpha_extrap[x_kannisto + 1]) /
      mTS_BE[2 + g, length(vY_EU)]
    
	  #vBeta_extrap <- (log(mMu_NL_tx_diff[length(vY_EU),]) - vAlpha_extrap)/mTS_BE[2+g,length(vY_EU)]

	  vA_extrap <- signif(vA_extrap,15)
      vB_extrap <- signif(vB_extrap,15)
      vAlpha_extrap <- signif(vAlpha_extrap,15)
      vBeta_extrap <- signif(vBeta_extrap,15)
	
    ####################################################################
    # Bepaal de sterfte intensiteiten voor leeftijden 91 t/m 120 o.b.v. de doorgetrokken A, B, alpha en beta (_extrap) voor alle jaren 
    # Loop over ages and years
    for(i in 1:(iT_his+iT_fut)) {
      for(j in x_kannisto+1) {
        mMu_EU_tx[i, j] <- exp( vA_extrap[j] + vB_extrap[j] * mTS_BE[g, i] )
        if(i >= iT_his_diff){
          mMu_NL_tx[i, j] <- exp( vA_extrap[j] + vB_extrap[j] * mTS_BE[g, i] + vAlpha_extrap[j] + vBeta_extrap[j] * mTS_BE[2+g, i] )
        }
      }
    }
    
    # Kort de Nederlandse sterfte intensiteiten in (starten bij 1983)
    mMu_NL_tx <- mMu_NL_tx[iT_his_diff:nrow(mMu_NL_tx),]
    
    ####################################################################
    
    #Translate mtx to qtx
    lQg_EU_tx[[vG[g]]] <- signif(1 - exp( - mMu_EU_tx ), digits = 15)
    lQg_NL_tx[[vG[g]]] <- signif(1 - exp( - mMu_NL_tx ), digits = 15)
    dimnames(lQg_EU_tx[[vG[g]]]) <- list(Year=seq(from=vY_EU[1], by=1, length.out=iT_his+iT_fut), Age=0:120)
    dimnames(lQg_NL_tx[[vG[g]]]) <- list(Year=seq(from=vY_NL[1], by=1, length.out=iT_his+iT_fut-iT_his_diff+1), Age=0:120)
    
    #Replace alpha and beta parameters with extrapolated ones
    objAG[[g]]$vA <- vA_extrap
    objAG[[g]]$vB <- vB_extrap
    objAG[[g]]$vAlpha <- vAlpha_extrap
    objAG[[g]]$vBeta <- vBeta_extrap
  }
  
  return(list(lQg_NL_tx=lQg_NL_tx, 
       lQg_EU_tx=lQg_EU_tx, 
       objAG=objAG))
}



# Project Kappa 
Project_kappa <- function(objAG, iT_fut=20) {
 
  # Estimate the TS model and extract regression results
  # objTS <- fnEstimateTS_AG2020(objAG=objAG, dRelTol=1E-18)
  objTS <- readRDS(here(run_folder, paste0(runname, "_filt_TS.rds")))
  mCov  <- objTS$Sigma
  mRho  <- objTS$Rho
  vDrift <- objTS$Drift
  vMR    <- c(1, 1, objTS$MR)
  
  # Simulate the time series
  vY_NL <- objAG[[1]]$vY_NL
  iT_his <- length(vY_NL)		# length of the historical time series
  mTS_BE <- matrix(data=0, nrow=4, ncol=iT_his + iT_fut, dimnames=list(TS=c("EU-M", "EU-F", "NL-M", "NL-F"), Year=seq(from=vY_NL[1], by=1, length.out=iT_his+iT_fut)))
  
  # First populate historical kappas
  mTS_BE["EU-M", as.character(vY_NL)] <- tail(objAG[["M"]]$vK, iT_his)
  mTS_BE["EU-F", as.character(vY_NL)] <- tail(objAG[["F"]]$vK, iT_his)
  mTS_BE["NL-M", as.character(vY_NL)] <- tail(objAG[["M"]]$vKappa, iT_his)
  mTS_BE["NL-F", as.character(vY_NL)] <- tail(objAG[["F"]]$vKappa, iT_his)
  
  # Now create forecasts, but loop over the years instead of over simulations
  for(i in 1:iT_fut) {
    mTS_BE[ , iT_his + i] <- vDrift + vMR * mTS_BE[ , iT_his + i - 1] 
  }
  return(mTS_BE)
}



# Generate stochastic simulations for mortality forecasts
fnProjectAG_sim_kannisto <- function(objAG, iT_fut=120, iSims=10) {
  
  # Estimate the TS model and extract regression results
  objTS <- fnEstimateTS_AG2020(objAG=objAG, dRelTol=1E-18)  
  
  mCov  <- objTS$Sigma
  mRho  <- objTS$Rho
  vDrift <- objTS$Drift
  vMR    <- c(1, 1, objTS$MR)
  
  # Simulate the time series
  vY_NL <- objAG[[1]]$vY_NL
  iT_his <- length(vY_NL)		# length of the historical time series
  mTS_BE <- matrix(data=0, nrow=4, ncol=iT_his + iT_fut, dimnames=list(TS=c("EU-M", "EU-F", "NL-M", "NL-F"), Year=seq(from=vY_NL[1], by=1, length.out=iT_his+iT_fut)))  
  aTS_sim <- array(data=0, dim=c(4, iT_his + iT_fut, iSims), dimnames=list(TS=c("EU-M", "EU-F", "NL-M", "NL-F"), Year=seq(from=vY_NL[1], by=1, length.out=iT_his+iT_fut), Sim=1:iSims))
  # First populate historical kappas
  aTS_sim["EU-M", as.character(vY_NL), ] <- matrix(data=tail(objAG[["M"]]$vK, iT_his), nrow=iT_his, ncol=iSims)
  aTS_sim["EU-F", as.character(vY_NL), ] <- matrix(data=tail(objAG[["F"]]$vK, iT_his), nrow=iT_his, ncol=iSims)
  aTS_sim["NL-M", as.character(vY_NL), ] <- matrix(data=tail(objAG[["M"]]$vKappa, iT_his), nrow=iT_his, ncol=iSims)
  aTS_sim["NL-F", as.character(vY_NL), ] <- matrix(data=tail(objAG[["F"]]$vKappa, iT_his), nrow=iT_his, ncol=iSims)
  # Now create forecasts, but loop over the years instead of over simulations
  if (fixed_Sims) set.seed(654987)
  for(i in 1:iT_fut) {
    mEps <- t(mvrnorm(n=iSims, mu=rep(0, 4), Sigma=mCov))
    for(j in 1:4) {
      aTS_sim[j, iT_his + i, ] <- vMR[j] * aTS_sim[j, iT_his + i-1, ] + vDrift[j] + mEps[j,]
    }
  }
  
  # Reconstruct mortality rates
  vX_compl <- 0:120	
  iX_compl <- length(vX_compl)
  iX_cal <- length(objAG[[1]]$vX)	
  iX_kan <- iX_compl - iX_cal
  vX <- objAG[[1]]$vX
  lMg_txi <- list()
  lQg_txi <- list()
  vG <- c("Males", "Females")
  
  for(g in 1:2) {
    # Extract parameter estimates to enable shorter notation below
    vA <- objAG[[g]]$vA
    vB <- objAG[[g]]$vB
    vAlpha <- objAG[[g]]$vAlpha
    vBeta <- objAG[[g]]$vBeta
    # Define output object
    # lMg_txi[[g]] <- array(data=0, dim=c(iT_his+iT_fut, iX_compl, iSims))
    lMg_txi[[vG[g]]] <- array(data=0, dim=c(iT_his+iT_fut, iX_compl, iSims), dimnames=list(Year=seq(from=min(objAG[[1]]$vY_NL), by=1, length=iT_his+iT_fut), Age=vX_compl, Sim=1:iSims))	
    # Loop over ages and years
    for(i in 1:(iT_his+iT_fut)) {
      for(j in 1:iX_cal) {
        lMg_txi[[vG[g]]][i, j, ] <- exp( vA[j] + vB[j] * aTS_sim[g, i, ] + vAlpha[j] + vBeta[j] * aTS_sim[2+g, i, ] )
        
      }
    }
    # Close mortality table		
    lMg_txi[[vG[g]]][ , iX_cal + 1:iX_kan, ] <- ApplyKannisto_mtxi(mtxi=lMg_txi[[vG[g]]][,1:iX_cal,], ages=vX, kannisto_nages=iX_kan, kannisto_nobs=11)
    lQg_txi[[vG[g]]] <- 1 - exp(-lMg_txi[[vG[g]]])
  }
  
  return(lQg_txi)
  
}



# Generate stochastic simulations for mortality forecasts
fnProjectAG_sim_alternatiefB <- function(objAG, iT_fut=120, iSims=10) {
  
  # Estimate the TS model and extract regression results
  objTS <- fnEstimateTS_AG2020(objAG=objAG, dRelTol=1E-18)  
  mCov  <- objTS$Sigma
  mRho  <- objTS$Rho
  vDrift <- objTS$Drift
  vMR    <- c(1, 1, objTS$MR)
  
  # Simulate the time series
  vY_NL <- objAG[[1]]$vY_NL
  iT_his <- length(vY_NL)		# length of the historical time series
  mTS_BE <- matrix(data=0, nrow=4, ncol=iT_his + iT_fut, dimnames=list(TS=c("EU-M", "EU-F", "NL-M", "NL-F"), Year=seq(from=vY_NL[1], by=1, length.out=iT_his+iT_fut)))  
  aTS_sim <- array(data=0, dim=c(4, iT_his + iT_fut, iSims), dimnames=list(TS=c("EU-M", "EU-F", "NL-M", "NL-F"), Year=seq(from=vY_NL[1], by=1, length.out=iT_his+iT_fut), Sim=1:iSims))
  # First populate historical kappas
  aTS_sim["EU-M", as.character(vY_NL), ] <- matrix(data=tail(objAG[["M"]]$vK, iT_his), nrow=iT_his, ncol=iSims)
  aTS_sim["EU-F", as.character(vY_NL), ] <- matrix(data=tail(objAG[["F"]]$vK, iT_his), nrow=iT_his, ncol=iSims)
  aTS_sim["NL-M", as.character(vY_NL), ] <- matrix(data=tail(objAG[["M"]]$vKappa, iT_his), nrow=iT_his, ncol=iSims)
  aTS_sim["NL-F", as.character(vY_NL), ] <- matrix(data=tail(objAG[["F"]]$vKappa, iT_his), nrow=iT_his, ncol=iSims)
  # Now create forecasts, but loop over the years instead of over simulations
  if (fixed_Sims) set.seed(654987)
  for(i in 1:iT_fut) {
    mEps <- t(mvrnorm(n=iSims, mu=rep(0, 4), Sigma=mCov))
    for(j in 1:4) {
      aTS_sim[j, iT_his + i, ] <- vMR[j] * aTS_sim[j, iT_his + i-1, ] + vDrift[j] + mEps[j,]
    }
  }
  remove(mEps)
  
  # Reconstruct mortality rates
  vX_compl <- 0:120	
  iX_compl <- length(vX_compl)
  iX_cal <- length(objAG[[1]]$vX)	
  iX_kan <- iX_compl - iX_cal
  vX <- objAG[[1]]$vX
  lQg_txi <- list()
  vG <- c("Males", "Females")
  
  for(g in 1:2) {
    # Extract parameter estimates to enable shorter notation below
    vA <- objAG[[g]]$vA
    vB <- objAG[[g]]$vB
    vAlpha <- objAG[[g]]$vAlpha
    vBeta <- objAG[[g]]$vBeta
    # Define output object
    aMu_txi <- array(data=0, dim=c(iT_his+iT_fut, iX_compl, iSims), dimnames=list(Year=seq(from=min(objAG[[1]]$vY_NL), by=1, length=iT_his+iT_fut), Age=vX_compl, Sim=1:iSims))	
    # Loop over ages and years
    for(i in 1:(iT_his+iT_fut)) {
      for(j in 1:iX_compl) {
        aMu_txi[i, j, ] <- exp( vA[j] + vB[j] * aTS_sim[g, i, ] + vAlpha[j] + vBeta[j] * aTS_sim[2+g, i, ] )
        
      }
    }
    
    # # Close mortality table		
    lQg_txi[[vG[g]]] <- 1 - exp(-aMu_txi)
  }
  return(lQg_txi)
  
}


# Closing function for a single mortality forecast
ApplyKannisto_mtx = function(mtx, ages, kannisto_nages=30, kannisto_nobs=15)
{
  if(is.null(mtx))
    return(NULL)

  # om de vele warnings dat log(...) resulteert in NaN te voorkomen
  mtx[(mtx>= 1) | (mtx <= 0)] <- NA

  #extract dimensions of mtxi
  n.ahead <- dim(mtx)[1]
  n.age   <- dim(mtx)[2]
    
  ###########################################################
  #extend using kannisto for each year 
  kannisto_obsages         <- tail(ages, kannisto_nobs)
  kannisto_extrapolateages <- tail(ages, 1) + 1:kannisto_nages
  kannisto_mtx            <- matrix(data=0, nrow=n.ahead, ncol=kannisto_nages)

  for(t in 1:(n.ahead)) {
      #obtain mx used for extrapolation
      kannisto_mx_obs <- tail(mtx[t,], kannisto_nobs)
    
      #apply transformation
      kannisto_mx_obs <- log(kannisto_mx_obs/(1-kannisto_mx_obs))
      
      #apply the regression and obtain phi1 and phi2
      kannisto_obsages[is.na(kannisto_mx_obs)] <-  NA
      kannisto_lm <- lm(kannisto_mx_obs ~ kannisto_obsages, na.action = na.omit)
      kannisto_phi1 <- exp(kannisto_lm$coefficient[1])
      kannisto_phi2 <- kannisto_lm$coefficient[2]
      
      #calculate qx and put into array
      kannisto_mtx[t, 1:kannisto_nages] <- kannisto_phi1*exp(kannisto_phi2*kannisto_extrapolateages)/(1+kannisto_phi1*exp(kannisto_phi2*kannisto_extrapolateages))	
 
  }
  kannisto_mtx
}
# Closing function for mortality forecast simulations
ApplyKannisto_mtxi = function(mtxi, ages, kannisto_nages=30, kannisto_nobs=15)
{
  if(is.null(mtxi))
    return(NULL)
  
  #extract dimensions of mtxi
  n.ahead <- dim(mtxi)[1]
  n.age   <- dim(mtxi)[2]
  n.sim   <- dim(mtxi)[3]
    
  ###########################################################
  #extend using kannisto for each year and simulation
  kannisto_obsages         <- tail(ages, kannisto_nobs)
  kannisto_extrapolateages <- tail(ages, 1) + 1:kannisto_nages
  kannisto_mtxi            <- array(dim=c(n.ahead, kannisto_nages, n.sim))
  
  for(t in 1:(n.ahead)) {
    
    for(i in 1:n.sim) {
      #obtain mx used for extrapolation
      kannisto_mx_obs <- tail(mtxi[t,,i], kannisto_nobs)

      #apply transformation
      kannisto_mx_obs <- log(kannisto_mx_obs/(1-kannisto_mx_obs))
      
      #apply the regression and obtain phi1 and phi2
      kannisto_lm <- lm(kannisto_mx_obs ~ kannisto_obsages)
      kannisto_phi1 <- exp(kannisto_lm$coefficient[1])
      kannisto_phi2 <- kannisto_lm$coefficient[2]
      
      #calculate qx and put into array
      kannisto_mtxi[t, 1:kannisto_nages, i] <- kannisto_phi1*exp(kannisto_phi2*kannisto_extrapolateages)/(1+kannisto_phi1*exp(kannisto_phi2*kannisto_extrapolateages))
	
    }
  }
  
  kannisto_mtxi 
}


fnApplyCOVIDparameters = function(objAG_BE, 
                                  objCOVIDparameters,
                                  eta_discount = 0.95,
                                  eta_limit = 0.25) {
  
  # Loop over the geslachten
  for (g in 1:2) {
    
    # Lees parameters uit lijst zodat de notatie verderop korter wordt
    COVIDjaren <- objCOVIDparameters[[g]]$COVIDjaren
    frak_beta  <- objCOVIDparameters[[g]]$frak_beta$X1
    frak_kappa <- objCOVIDparameters[[g]]$frak_kappa$X1
    
    # Transformeer matrix met sterftekansen naar een matrix met sterfte-intensiteiten  
    lMg_tx <- - log(1 - objAG_BE[[g]])
    
    # Pas de COVID parameters toe, per jaar
    jaren <- as.numeric(rownames(lMg_tx))
    
    # Loop over de jaren in de tafel
    for (t in jaren) {
      idx <- match(t, rownames(lMg_tx))  # index van betreffende covidjaar
      
      # Er zijn 3 situaties:
      #   het is een jaar v??r corona -> doe niks
      #   het is een jaar met corona weekdata -> pas ingelezen parameters toe
      #   of het is een prognosejaar -> trek frak_kappa door met eta en pas parameters daarna toe
      
      if (t %in% COVIDjaren) {
        # pas de ingelezen COVID paramters toe
        idx_covid <- match(t, COVIDjaren)  # index van betreffende covidjaar
        lMg_tx[idx,] <- lMg_tx[idx,] * exp(frak_beta * frak_kappa[idx_covid])
        
      } else if (t > max(COVIDjaren)) {
        
        # trek frak_kappa door
        idx_covid <- match(max(COVIDjaren), COVIDjaren)  # index van het laatste beschikbare covidjaar
        h <- t - max(COVIDjaren)
        
        frak_kappa_start <- frak_kappa[idx_covid]
        eta_star <- eta_limit + ((1 - eta_limit) * eta_discount^h)
        frak_kappa_t <- frak_kappa_start * eta_star
        
        # pas de COVID paramters toe
        lMg_tx[idx,] <- lMg_tx[idx,] * exp(frak_beta * frak_kappa_t)
      }
    }
    
    # Transformeer matrix met sterfte-intensiteiten terug naar een matrix met sterftekansen
    objAG_BE[[g]] <- 1 - exp(-lMg_tx)
    
    # Rond sterftekansen af op 15 significante cijfers
    objAG_BE[[g]] <- signif(objAG_BE[[g]], digits = 15)
  }
  
  return(objAG_BE)
}


fnApplyCOVIDparameters_all = function(objAG_BE_all, 
                                      objCOVIDparameters,
                                      eta_discount = 0.95,
                                      eta_limit = 0.25
                                      ) {
  
  # Loop over the geslachten
  for (g in 1:2) {
    
    # Lees parameters uit lijst zodat de notatie verderop korter wordt
    COVIDjaren <- objCOVIDparameters[[g]]$COVIDjaren
    frak_beta  <- objCOVIDparameters[[g]]$frak_beta$X1
    frak_kappa <- objCOVIDparameters[[g]]$frak_kappa$X1
    
    # Transformeer matrix met sterftekansen naar een matrix met sterfte-intensiteiten  
    lMg_tx <- - log(1 - objAG_BE_all$lQg_NL_tx[[g]])
    
    # Pas de COVID parameters toe, per jaar
    jaren <- as.numeric(rownames(lMg_tx))
    
    # Loop over de jaren in de tafel
    for (t in jaren) {
      idx <- match(t, rownames(lMg_tx))  # index van betreffende covidjaar
      
      # Er zijn 3 situaties:
      #   het is een jaar v??r corona -> doe niks
      #   het is een jaar met corona weekdata -> pas ingelezen parameters toe
      #   of het is een prognosejaar -> trek frak_kappa door met eta en pas parameters daarna toe
      
      if (t %in% COVIDjaren) {
        # pas de ingelezen COVID paramters toe
        idx_covid <- match(t, COVIDjaren)  # index van betreffende covidjaar
        lMg_tx[idx,] <- lMg_tx[idx,] * exp(frak_beta * frak_kappa[idx_covid])
        
      } else if (t > max(COVIDjaren)) {
        
        # trek frak_kappa door
        idx_covid <- match(max(COVIDjaren), COVIDjaren)  # index van het laatste beschikbare covidjaar
        h <- t - max(COVIDjaren)
        
        frak_kappa_start <- frak_kappa[idx_covid]
        eta_star <- eta_limit + ((1 - eta_limit) * eta_discount^h)
        frak_kappa_t <- frak_kappa_start * eta_star
        
        # pas de COVID paramters toe
        lMg_tx[idx,] <- lMg_tx[idx,] * exp(frak_beta * frak_kappa_t)
      }
    }
    
    # Transformeer matrix met sterfte-intensiteiten terug naar een matrix met sterftekansen
    objAG_BE_all$lQg_NL_tx[[g]] <- 1 - exp(-lMg_tx)
    
    objAG_BE_all$lQg_NL_tx[[g]] <- signif(objAG_BE_all$lQg_NL_tx[[g]], digits = 15)
  }
  
  return(objAG_BE_all)
}

