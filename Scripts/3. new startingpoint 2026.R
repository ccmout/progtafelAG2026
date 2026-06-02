
fnApplyNewStartingpoints <- function(s, c, startpunt_convergentie, objAG_BE, lData_NL, t_max, VY_NL_s) {

  for (geslacht in c("Male", "Female")) {
  
    gender <- substring(geslacht,1,1)
    
    ### Jaren waarover gemiddeld wordt
    years <- seq(t_max, t_max - s + 1, by = -1)
    t_ref <- mean(years)
    
   
    years <- seq(t_max, t_max - s + 1, by = -1)
    t_ref <- mean(years)
    
    Dxt_years_pre <- lData_NL[[gender]][["Dxt"]][,(years-vY_NL_s+1), drop = FALSE]
    Dxt_years <- Dxt_years_pre
    Ext_years_pre <- lData_NL[[gender]][["Ext"]][,(years-vY_NL_s+1), drop = FALSE]
    Ext_years <- Ext_years_pre
    
    Mxt_years_pre <- log(Dxt_years/Ext_years)
    Mxt_years <- Mxt_years_pre[1:91, ] ### Remove Age 91
    
    ### Voor leeftijden 2 t/m 90, middel over lft - 1, lft, lfd + 1
    for (lft in 3:91) {
      Mxt_years[lft,] <- (Mxt_years_pre[lft - 1,] + Mxt_years_pre[lft,] + Mxt_years_pre[lft + 1,])/3
    }
    
    
    
    ### Middel log(mu)
    mM <- rowSums(Mxt_years)/s
    
    ### Leeftijden Kannisto
    iX_compl <- length(0:120)
    iX_cal <- length(vX)	
    iX_kan <- iX_compl - iX_cal
    
    ### Geschatte mu (take exponent)
    mM_exp = rbind(exp(mM), exp(mM)) #doubled just so formula works
    
    ### Kanisto methodiek, startpunten (log) mu
    Mxt_years_kan <- log(ApplyKannisto_mtx(mtx=mM_exp[ , 1:iX_cal], ages=vX, kannisto_nages=iX_kan, kannisto_nobs=11)[1,])
    mM_tot <- exp(c(mM, Mxt_years_kan))
    names(mM_tot) <- 0:120
    
    ### Mu projected and mortality impr
    mu_proj <- -log(1 - objAG_BE_all$lQg_NL_tx[[geslacht]])
    mort_impr <- -log( mu_proj[-1, ] / mu_proj[-nrow(mu_proj), ])
    rownames(mort_impr) <- rownames(mu_proj)[-1]
    
    for (i in 1:121) {
      
      future_years <- ceiling(t_ref):rownames(mort_impr)[nrow(mort_impr)]
      mu_adj <- matrix(,nrow = length(future_years), ncol = 121)
      rownames(mu_adj) <- future_years
      colnames(mu_adj) <- 0:120
      
      
      idx <- which(rownames(mort_impr) == ceiling(t_ref))
     
      count = 1
      mu_adj[count,] <- mM_tot
      count = count + 1
      
      for (j in (idx+1):nrow(mort_impr)) {
        mu_adj[count,] <- mu_adj[(count-1),]*exp(-mort_impr[j,])
        count = count+1
      }
      
      
      mu_orig = mu_proj[(nrow(mu_proj)-nrow(mu_adj)+1):nrow(mu_proj),]
      
      q_adj <-  1 - exp(-mu_adj)
      q_orig <- 1 - exp(-mu_orig)
     
      ### Neem oorspronkelijke geschatte sterfkansen
      q_c = q_orig
      
      ### Voor convergentie, nemen we de adjusted sterftekansen
      for (year in (ceiling(t_ref):startpunt_convergentie)){
        year_index <- as.character(year)
        q_c[year_index,] <- q_adj[year_index,]
        objAG_BE[[geslacht]][year_index, ] = q_c[year_index, ]
      }
      
      ### Convergeren vanaf startpunt (2019)
      for (h in (0:c)){
        year_index = as.character(startpunt_convergentie + h)
        q_c[year_index, ] = (
          (1 - min(h, c) / c) * q_adj[year_index, ] + min(h, c) / c * q_orig[year_index,]
        )
        objAG_BE[[geslacht]][year_index, ] = q_c[year_index, ]
      }
    }
  }
  objAG_BE
}
################################################################################

fnApplyNewStartingpoints_all <- function(s, c, startpunt_convergentie, objAG_BE_all, lData_NL, t_max, VY_NL_s) {
  
  for (geslacht in c("Male", "Female")) {
    
    gender <- substring(geslacht,1,1)
    ### Jaren waarover gemiddeld wordt
    years <- seq(t_max, t_max - s + 1, by = -1)
    t_ref <- mean(years)
    
    ### Deaths and Exposures
    Dxt_years <- lData_NL[[gender]][["Dxt"]][,(years-vY_NL_s+1), drop = FALSE]
    Ext_years <- lData_NL[[gender]][["Ext"]][,(years-vY_NL_s+1), drop = FALSE]
    
    Dxt_years_pre <- lData_NL[[gender]][["Dxt"]][,(years-vY_NL_s+1), drop = FALSE]
    Dxt_years <- Dxt_years_pre
    Ext_years_pre <- lData_NL[[gender]][["Ext"]][,(years-vY_NL_s+1), drop = FALSE]
    Ext_years <- Ext_years_pre
    
    Mxt_years_pre <- log(Dxt_years/Ext_years)
    Mxt_years <- Mxt_years_pre[1:91,] ### Remove age 91
    
    for (lft in 3:91) {
      Mxt_years[lft,] <- (Mxt_years_pre[lft - 1,] + Mxt_years_pre[lft,] + Mxt_years_pre[lft + 1,])/3
    }
    
    ### Middel log(mu)
    mM <- rowSums(Mxt_years)/s
    
    ### Leeftijden Kannisto
    iX_compl <- length(0:120)
    iX_cal <- length(vX)	
    iX_kan <- iX_compl - iX_cal
    
    ### Geschatte mu (take exponent)
    mM_exp = rbind(exp(mM), exp(mM)) #doubled just so formula works
    
    ### Kanisto methodiek, startpunten (log) mu
    Mxt_years_kan <- log(ApplyKannisto_mtx(mtx=mM_exp[ , 1:iX_cal], ages=vX, kannisto_nages=iX_kan, kannisto_nobs=11)[1,])
    mM_tot <- exp(c(mM, Mxt_years_kan))
    names(mM_tot) <- 0:120
    
    ### Mu projected and mortality impr
    mu_proj <- -log(1 - objAG_BE_all$lQg_NL_tx[[geslacht]])
    mort_impr <- -log( mu_proj[-1, ] / mu_proj[-nrow(mu_proj), ])
    rownames(mort_impr) <- rownames(mu_proj)[-1]
    
    for (i in 1:121) {
      
      future_years <- ceiling(t_ref):rownames(mort_impr)[nrow(mort_impr)]
      mu_adj <- matrix(,nrow = length(future_years), ncol = 121)
      rownames(mu_adj) <- future_years
      colnames(mu_adj) <- 0:120
      
      
      idx <- which(rownames(mort_impr) == ceiling(t_ref))
      
      count = 1
 
      mu_adj[count,] <- mM_tot
      count = count + 1
      
      for (j in (idx+1):nrow(mort_impr)) {
        mu_adj[count,] <- mu_adj[(count-1),]*exp(-mort_impr[j,])
        count = count+1
      }
      
      mu_orig = mu_proj[(nrow(mu_proj)-nrow(mu_adj)+1):nrow(mu_proj),]
      
      q_adj <-  1 - exp(-mu_adj)
      q_orig <- 1 - exp(-mu_orig)
      
      ### Neem oorspronkelijke geschatte sterfkansen
      q_c = q_orig
      
      ### Voor convergentie, nemen we de adjusted sterftekansen
      for (year in (ceiling(t_ref):startpunt_convergentie)){
        year_index <- as.character(year)
        q_c[year_index,] <- q_adj[year_index,]
        objAG_BE_all$lQg_NL_tx[[geslacht]][year_index, ] = q_c[year_index, ]
      }
      
      ### Convergeren vanaf startpunt (2019)
      for (h in (0:c)){
        year_index = as.character(startpunt_convergentie + h)
        q_c[year_index, ] = (
          (1 - min(h, c) / c) * q_adj[year_index, ] + min(h, c) / c * q_orig[year_index,]
        )
        objAG_BE_all$lQg_NL_tx[[geslacht]][year_index, ] = q_c[year_index, ]
      }
    }
  }
  return(objAG_BE_all)
}