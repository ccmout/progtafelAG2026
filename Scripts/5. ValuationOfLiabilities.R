library(readxl)
library(openxlsx)
library(dplyr)
library(here)

### VOORdat je dit runt, doe ff ctrl + shift + f10

runname <- "4. AG2026 incl. startpuntmiddeling, COVID 2022-2023"
run_folder <- here("Runs", "Deelopdracht 7", runname)

Q2024M <- read_excel(here(run_folder, paste0(runname, "_BE_proj_all.xlsx")), sheet = "NLMale")
Q2024F <- read_excel(here(run_folder, paste0(runname, "_BE_proj_all.xlsx")), sheet = "NLFemale")

waarderingsjaar <- 2027
Q2024M <- Q2024M %>% dplyr::select(as.character(waarderingsjaar):last_col())
Q2024F <- Q2024F %>% dplyr::select(as.character(waarderingsjaar):last_col())

mortality_table <- array(c(as.matrix(Q2024M), as.matrix(Q2024F)), dim = c(nrow(as.matrix(Q2024M)), ncol(as.matrix(Q2024M)), 2))
# Qlist <- list(mortaility_table)

objPF <- readRDS("Input data/objPF.rds")

# fnPV_factors <- function(aQtxg, dtPF, dR=0.03, iXR=65) {  ##ONlY LIABILITIES
#   # Assumption is that the array does not contain any historical years that should not be used for valuation!
#   vX <- 0:120
#   iN <- 14
#   mOutput <- matrix(data=0, nrow=iN, ncol=3)
#   for(participant in 1:iN) {
#     # participant = 1
#     sGeslacht <- dtPF[participant, 2]
#     geslacht <- ifelse(sGeslacht=="M", 1, 2)
#     leeftijd <- dtPF[participant, 1]
#     
#     vAges <- seq(from=leeftijd, to=120, by=1)
#     
#     # Make sure at least enough rows are available
#     mPx <- matrix(data=0, ncol=3, nrow=120)
#     vDF_p <- (1+dR)^( -(0:119) )		
#     vDF_u <- (1+dR)^( -(1:120) )
#     mCF <- matrix(data=seq(from=leeftijd, by=1, length.out=120), ncol=3, nrow=120, byrow=FALSE)
#     mCF[,1] <- 1 * (mCF[,1] >= iXR)
#     mCF[,c(2,3)] <- 1
#     
#     # OP - Hoofdverzekerde
#     if(dtPF[participant, 3] > 0) {	
#       vAge_ind <- (vX >= leeftijd)
#       vQx <- diag(aQtxg[vAge_ind, , geslacht])
#       vPx <- 1-vQx	
#       vPx_cum <- c(1, cumprod(vPx))
#       vkPx_p <- head(vPx_cum, -1)
#       vkPx_u <- tail(vPx_cum, -1)
#       iProj <- length(vkPx_u)
#       mOutput[participant, 1] <- sum(mCF[1:iProj,1] * (vkPx_p * vDF_p[1:iProj] + vkPx_u * vDF_u[1:iProj]) / 2)
#     }
#     
#     # Ingegaan NP - Nabestaande
#     if(dtPF[participant, 5] > 0) {	
#       vAge_ind <- (vX >= leeftijd)
#       vQx <- diag(aQtxg[vAge_ind, , geslacht])	# gender does not need to be switched, it is already switched in the data table
#       vPx <- 1-vQx	
#       vPx_cum <- c(1, cumprod(vPx))
#       vkPx_p <- head(vPx_cum, -1)
#       vkPx_u <- tail(vPx_cum, -1)
#       iProj <- length(vkPx_u)
#       mOutput[participant, 3] <- sum(mCF[1:iProj,3] * (vkPx_p * vDF_p[1:iProj] + vkPx_u * vDF_u[1:iProj]) / 2)
#     }
#     
#     # Latent NP - Medeverzekerde
#     if(dtPF[participant, 4] > 0) {	
#       # Obtain unadjusted mortality rates
#       leeftijd_spouse <- ifelse(geslacht==1, leeftijd - 3, leeftijd + 3)
#       vAge_ind_m <- (vX >= leeftijd)
#       vAge_ind_s <- (vX >= leeftijd_spouse)
#       vQx_m <- diag(aQtxg[vAge_ind_m, , geslacht])
#       vQx_s <- diag(aQtxg[vAge_ind_s, , 3 - geslacht])
#       vQx_s_star <- vQx_s
#       vAge_sm <- vX[vX >= leeftijd_spouse] + ifelse(geslacht==1, +3, -3)
#       vQx_s_star[vAge_sm < iXR] <- 0	# hx = 1
#       if(geslacht==1) vQx_m <- c(vQx_m, 1, 1, 1) else { vQx_s <- c(vQx_s, 1, 1, 1);   vQx_s_star <- c(vQx_s_star, 1, 1, 1) }
#       vX_m <- seq(from = leeftijd, by=1, length.out=length(vQx_m))
#       
#       vPx_m <- 1-vQx_m
#       vPx_s <- 1-vQx_s
#       vPx_s_star <- 1-vQx_s_star
#       vkPx_m <- c(1, head(cumprod(vPx_m), -1))
#       vkPx_s <- c(1, head(cumprod(vPx_s), -1))
#       vkPx_s_star <- c(1, head(cumprod(vPx_s_star), -1))
#       
#       # Adjust mortality rates of spouse to 1 if main is not yet retired
#       vtPxy <- numeric(length(vQx_s))
#       for(j in 2:length(vQx_s)) {
#         vtPxy[j] <- vtPxy[j-1] * vPx_s[j-1] + vQx_m[j-1] * vkPx_m[j-1] * vkPx_s_star[j] * ifelse(vX_m[j]>iXR, 1, sqrt(vPx_s[j-1]))
#       }
#       mOutput[participant, 2] <- sum( vtPxy * vDF_p[1:length(vtPxy)])
#     }
#   }
#   
#   return(mOutput)
#   
# }

fnValueLiabilities <- function(aQtxg, dtPF, dR=0.03, iXR=65, liabilities=2) {
  # liabilities = 1 indicates calculation of PV's for liabilities, value of 2 equals calculation of premium factors 
  # Assumption is that the array does not contain any historical years that should not be used for valuation!
  iN <- nrow(dtPF)
  iC <- dim(dtPF)[2] - 2  # minus two for columns containing age & gender
  if(liabilities ==1 ){
    mOutput <- matrix(data=0, nrow=iN, ncol= iC)  
  } else{
    mOutput <- matrix(data=0, nrow=iN, ncol= iC+1)    # more extensive output than input due to distinction PP risk & capital
  }
  #vX <- as.integer(dimnames(aQtxg)[[2]])
  vX <- 0:120
  for(participant in 1:iN) {
    #participant = 1
    sGeslacht <- dtPF[participant, 2]
    geslacht <- ifelse(sGeslacht=="M", 1, 2)
    leeftijd <- dtPF[participant, 1]
    
    vAges <- seq(from=leeftijd, to=120, by=1)
    
    # Make sure at least enough rows are available
    mPx <- matrix(data=0, ncol=3, nrow=120)
    vDF_p <- (1+dR)^( -(0:119) )		
    vDF_u <- (1+dR)^( -(1:120) )
    mCF <- matrix(data=seq(from=leeftijd, by=1, length.out=120), ncol=iC, nrow=120, byrow=FALSE)
    mCF[,1] <- dtPF[participant, 3] * (mCF[,1] >= iXR)
    mCF[,2] <- dtPF[participant, 4] 
    if(liabilities==1){
      mCF[,3] <- dtPF[participant, 5]   
    }else{
      mCF[,3] <- 0
    }
    
    # Hoofdverzekerde
    if(dtPF[participant, 3] > 0) {	
      vAge_ind <- (vX >= leeftijd)  
      vQx <- diag(aQtxg[vAge_ind, , geslacht])
      vPx <- 1-vQx	
      vPx_cum <- c(1, cumprod(vPx))
      vkPx_p <- head(vPx_cum, -1)
      vkPx_u <- tail(vPx_cum, -1)
      iProj <- length(vkPx_u)
      mOutput[participant, 1] <- sum(mCF[1:iProj,1] * (vkPx_p * vDF_p[1:iProj] + vkPx_u * vDF_u[1:iProj]) / 2) 
    }
    
    # Medeverzekerde
    if(dtPF[participant, 4] > 0) {	
      # Obtain unadjusted mortality rates
      leeftijd_spouse <- ifelse(geslacht==1, leeftijd - 3, leeftijd + 3)
      vAge_ind_m <- (vX >= leeftijd)
      vAge_ind_s <- (vX >= leeftijd_spouse)
      vQx_m <- diag(aQtxg[vAge_ind_m, , geslacht])
      vQx_s <- diag(aQtxg[vAge_ind_s, , 3 - geslacht])
      print(leeftijd_spouse)
      vQx_s_star <- vQx_s
      vAge_sm <- vX[vX >= leeftijd_spouse] + ifelse(geslacht==1, +3, -3)
      if(liabilities == 1){
        vQx_s_star[vAge_sm < iXR] <- 0	# hx = 1		
        if(geslacht==1) vQx_m <- c(vQx_m, 1, 1, 1) else { vQx_s <- c(vQx_s, 1, 1, 1);   vQx_s_star <- c(vQx_s_star, 1, 1, 1) }
        vX_m <- seq(from = leeftijd, by=1, length.out=length(vQx_m))
        
        vPx_m <- 1-vQx_m
        vPx_s <- 1-vQx_s
        vPx_s_star <- 1-vQx_s_star
        vkPx_m <- c(1, head(cumprod(vPx_m), -1))
        vkPx_s <- c(1, head(cumprod(vPx_s), -1))
        vkPx_s_star <- c(1, head(cumprod(vPx_s_star), -1))
        
        # Adjust mortality rates of spouse to 1 if main is not yet retired
        vtPxy <- numeric(length(vQx_s))
        
        #print(aQtxg[, vAge_ind_s, 2])
        for(j in 2:length(vQx_s)) {
          vtPxy[j] <- vtPxy[j-1] * vPx_s[j-1] + vQx_m[j-1] * vkPx_m[j-1] * vkPx_s_star[j] * ifelse(vX_m[j]>iXR, 1, sqrt(vPx_s[j-1]))
        }
        mOutput[participant, 2] <- sum( mCF[1:length(vtPxy),2] * vtPxy * vDF_p[1:length(vtPxy)])
      } else{ #only for premium impact; difference in calculation 
        ## PP opbouw
        hx <- vQx_s; hx[] <- 1  # create array with equal length to vQx for marriage frequency
        LatPP <- hx[] - 1
        for(age in 1:length(vQx_m)){
          if(age <= length(vQx_s)){
            if(vAge_sm[age]<=iXR){
              vQx_s_star[age] <- 1
            } else{
              vQx_s_star[age] <- (1-vQx_s[age-1])*vQx_s_star[age-1]
              hx[age-1] <- vQx_s_star[age]^0.5*vQx_s_star[age-1]^0.5
            }
            if(age>1){
              LatPP[age] <- (1-vQx_s[age-1])*LatPP[age-1]+vPx_cum[age-1]*vQx_m[age-1]*hx[age-1]*sqrt(1-vQx_s[age-1])
            }
          }
        }
        PP_opbouw <- sum(mCF[1:length(LatPP),2] * LatPP * vDF_p[1:length(LatPP)])
        ## PP risicobasis (incl toekenning)
        dienstjaren <- 40; PP_perc <- 0.7
        PP_risico_bedrag <- dtPF[participant, 3]*dienstjaren*PP_perc  
        PP_risico_toekenning_bedrag<- dtPF[participant, 3]*(iXR-1-leeftijd)*PP_perc  
        tpy_LAtPP <- LatPP;tpy_LAtPP[2] <- 1    # Both arrays start with 0 observation & equally long
        for(i in 3:length(tpy_LAtPP)){
          tpy_LAtPP[i] <- tpy_LAtPP[i-1] * (1-vQx_s[i-1])
        }
        PP_risico <- PP_risico_bedrag*sum(tpy_LAtPP* vDF_p[1:length(tpy_LAtPP)])*vQx_m[1]*hx[1]*(1+vQx_s[1])^(-0.5)
        PP_risico_toekenning <- PP_risico_toekenning_bedrag*sum(tpy_LAtPP* vDF_p[1:length(tpy_LAtPP)])*vQx_m[1]*hx[1]*(1+vQx_s[1])^(-0.5)
        
        mOutput[participant, 2] <- PP_opbouw + PP_risico_toekenning
        mOutput[participant, 3] <- PP_risico
      }
    }
    
    # Nabestaande, enkel voor verplichting niet voor premie
    if(liabilities ==1){
      if(dtPF[participant, 5] > 0) {	
        vAge_ind <- (vX >= leeftijd)
        vQx <- diag(aQtxg[vAge_ind, , geslacht])	# gender does not need to be switched, it is already switched in the data table
        vPx <- 1-vQx	
        vPx_cum <- c(1, cumprod(vPx))
        vkPx_p <- head(vPx_cum, -1)
        vkPx_u <- tail(vPx_cum, -1)
        iProj <- length(vkPx_u)
        mOutput[participant, 3] <- sum(mCF[1:iProj,3] * (vkPx_p * vDF_p[1:iProj] + vkPx_u * vDF_u[1:iProj]) / 2)
      }
    }
    
    
  }
  
  return(mOutput)
}

fnOPNPTP <- function(mOutput,liabilities=1) {
  if(liabilities ==1 ){
    vOutput <- c(sum(mOutput[,1]), sum(mOutput[,2]), sum(mOutput[,3]), sum(mOutput))
  } else{
    vOutput <- c(sum(mOutput[,1]), sum(mOutput[,2]), sum(mOutput[,3]), sum(mOutput[,1],mOutput[,2]), sum(mOutput[,1],mOutput[,3]))
  }
  return(vOutput)
}

fnOPNPTPgeslacht <- function(mOutput, dtPF, liabilities=1) {
  idx_man <- which(dtPF$Geslacht == "M")
  idx_vrouw <- which(dtPF$Geslacht == "V")
  if(liabilities ==1 ){
    vOutput <- c(sum(mOutput[idx_man ,1]), 
                sum(mOutput[idx_man ,2]), 
                sum(mOutput[idx_man ,3]), 
                sum(mOutput[idx_man, ]), 
                sum(mOutput[idx_vrouw ,1]), 
                sum(mOutput[idx_vrouw ,2]), 
                sum(mOutput[idx_vrouw ,3]),
                sum(mOutput[idx_vrouw, ]),
                sum(mOutput[, 1]),
                sum(mOutput[, 2]),
                sum(mOutput[, 3]),
                sum(mOutput))
  } else{
    vOutput <- c(sum(mOutput[idx_man, 1]), 
                 sum(mOutput[idx_man, 2]), 
                 sum(mOutput[idx_man, 3]), 
                 sum(mOutput[idx_man,1],mOutput[idx_man,2]), 
                 sum(mOutput[idx_man,1],mOutput[idx_man,3]),
                 sum(mOutput[idx_vrouw,1]), 
                 sum(mOutput[idx_vrouw,2]), 
                 sum(mOutput[idx_vrouw,3]), 
                 sum(mOutput[idx_vrouw,1],mOutput[idx_vrouw,2]), 
                 sum(mOutput[idx_vrouw,1],mOutput[idx_vrouw,3]),
                 sum(mOutput[,1]), 
                 sum(mOutput[,2]), 
                 sum(mOutput[,3]), 
                 sum(mOutput[,1],mOutput[,2]), 
                 sum(mOutput[,1],mOutput[,3]))
  }
  return(vOutput)
}

fnValueLiabilities_endResults <- function(aQtxg, dtPF, dR=0.03, iXR=65, liabilities=1) {
  mOutput <- fnValueLiabilities(aQtxg=aQtxg, dtPF=dtPF, dR=dR, iXR=iXR, liabilities) 
  vOutput <- fnOPNPTP(mOutput, liabilities)
  geslacht_output <- fnOPNPTPgeslacht(mOutput, dtPF, liabilities)
  return(list(vOutput = vOutput,
              vOutputInclGeslacht = geslacht_output))
}


################### Zes Portefeuilles ###############################

Output <- matrix(0, nrow = 6, ncol = 6)
colnames(Output) <- c("Tafel", "Fonds", "OP", "Latent_NP", "Ingegaan_NP", "Totaal")
Q_names <- c("Q")
Q_col <- c(Q_names, Q_names, Q_names, Q_names, Q_names, Q_names)
Output[,1] <- Q_col

Output2 <- matrix(0, nrow = 6, ncol = 14)
colnames(Output2) <- c("Tafel", "Fonds", 
                      "OP man", "Latent_NP man", "Ingegaan_NP man", "Totaal man",
                      "OP vrouw", "Latent_NP vrouw", "Ingegaan_NP vrouw", "Totaal vrouw",
                      "OP", "Latent_NP", "Ingegaan_NP", "Totaal")
Q_names <- c("Q")
Q_col <- c(Q_names, Q_names, Q_names, Q_names, Q_names, Q_names)
Output2[,1] <- Q_col

n <- 1
for (p in 1:6){
  print(p)
  Output[n,2] <- p
  Output2[n,2] <- p
  results <- fnValueLiabilities_endResults(mortality_table, objPF[[p]], liabilities = 1)
  Output[n,3:6] <- results$vOutput
  Output2[n, 3:14] <- results$vOutputInclGeslacht
  n <- n + 1
  write.xlsx(list(Output = Output,
                  output2 = Output2),
             here(run_folder, paste0(runname, " portefeuille_resultaten.xlsx", sep="")))
  }

################## Factoren ############################

factoren_portefeuille <- objPF[[1]]
factoren_portefeuille[, (ncol(factoren_portefeuille) - 2): ncol(factoren_portefeuille)] <- 1
factoren_portefeuille[, (ncol(factoren_portefeuille) - 2): ncol(factoren_portefeuille)] <- fnValueLiabilities(aQtxg=mortality_table, dtPF=factoren_portefeuille, dR=0.03, iXR=65, liabilities=1)
write.xlsx(factoren_portefeuille, 
           here(run_folder, paste0(runname, " factoren.xlsx")))

########################## Premie ######################################
# 
objPF <- readRDS("Input data/objPF.rds")
objPremie <- matrix(0, nrow = 8, ncol = 4)
objPremie[,1] <- rbind(30, 40, 50, 60, 30, 40, 50, 60)
objPremie[,2] <- rbind("M", "M", "M", "M", "V", "V", "V", "V")
objPremie[,3] <- rbind(600, 750, 800, 600, 400, 500, 550, 400)
objPremie[,4] <- rbind(420, 525, 560, 420, 280, 350, 385, 280)
objPremie <- as.data.frame(objPremie)
cols_to_fix <- c(1,3,4)
objPremie[cols_to_fix] <- lapply(objPremie[cols_to_fix], as.integer)
colnames(objPremie) <- c("Leeftijd", "Geslacht", "OP", "PP (lat.)")
objPremie$`NP (ing.)` <- 0

Output_base <- matrix(0, nrow = 1, ncol = 7)
colnames(Output_base) <- c("Tafel", "Fonds", "OP", "PP_opbouw + PP_risico_toekenning", "PP_risico",	"OP + PP_opbouw + PP_risico_toekenning", "OP + PP_risico")
Q_names <- c("Q")
Q_col <- c(Q_names)
Output_base[,1] <- Q_col

Output_base2 <- matrix(0, nrow = 1, ncol = 17)
colnames(Output_base2) <-c("Tafel", "Fonds",
                          "OP man", "PP_opbouw + PP_risico_toekenning man", "PP_risico man",	"OP + PP_opbouw + PP_risico_toekenning man", "OP + PP_risico man",
                          "OP vrouw", "PP_opbouw + PP_risico_toekenning vrouw", "PP_risico vrouw",	"OP + PP_opbouw + PP_risico_toekenning vrouw", "OP + PP_risico vrouw",
                          "OP", "PP_opbouw + PP_risico_toekenning", "PP_risico",	"OP + PP_opbouw + PP_risico_toekenning", "OP + PP_risico")
Q_names <- c("Q")
Q_col <- c(Q_names)
Output_base2[,1] <- Q_col

n <- 1
for (p in 1:1){
  Output_base[n,2] <- p
  Output_base2[n, 2] <- p
  results <- fnValueLiabilities_endResults(mortality_table, objPremie, iXR=68, liabilities = 2)
  Output_base[n,3:7] <- results$vOutput
  Output_base2[n,3:17] <- results$vOutputInclGeslacht
  n <- n + 1
}

write.xlsx(list(output = Output_base,
                output2 = Output_base2),
           here(run_folder, paste0(runname, " Liabilities_premies.xlsx")))


