## Code to determine cohort and period life expectancy 

fnLifeExpCoh<- function(aQtxg, age_ult) {
  # Assumption is that the array does not contain any historical years that should not be used for determination of life expectancy!
  
  years <- dim(aQtxg)[1]
  ages <- dim(aQtxg)[2]
  
  P <- 1 - aQtxg
  
  # create empty matrices cohort life expectancy
  ec <- matrix(nrow = years,ncol = ages); tp <- matrix(0,ncol = ages, nrow = ages)
  for (t in 1:years){
    x <- 1
    tp[x,] <- P[t,]; 
    for(x in 2:(age_ult+1)){
      if(x+t-1>years){
        tp[x,1:((age_ult+1)-x+1)] <- tp[x-1,1:((age_ult+1)-x+1)]*P[years,x:(age_ult+1)]
      } else {
        tp[x,1:((age_ult+1)-x+1)] <- tp[x-1,1:((age_ult+1)-x+1)]*P[x+t-1,x:(age_ult+1)]
      }
    }
    "Calculate life expectancy for all starting ages "
    ec[t,] <- +0.5 + apply(tp,2,sum)
  }
  
  return(ec)
}

fnLifeExpPer<- function(aQtxg, age_ult) {
  # Assumption is that the array does not contain any historical years that should not be used for determination of life expectancy!
  
  years <- dim(aQtxg)[1]
  ages <- dim(aQtxg)[2]
  
  P <- 1 - aQtxg
  
  # create empty matrices period life expectancy
  ep <- matrix(nrow = years,ncol = ages); tp <- matrix(0,ncol = ages, nrow = ages)
  for (t in 1:years){
    x <- 1
    tp[x,] <- P[t,]; 
    for(x in 2:(age_ult+1)){
      tp[x,1:((age_ult+1)-x+1)] <- tp[x-1,1:((age_ult+1)-x+1)]*P[t,x:(age_ult+1)]
    }
    "Calculate life expectancy for all starting ages "
    ep[t,] <- +0.5 + apply(tp,2,sum)
  }
  
  return(ep)
}

