model{
  ## ## ## ## ## ## ## ## ## ## ## ## ## ## #
  ####      Simplified L2SWBM model      ####
  ###          JAGS/MODEL SCRIPT          ###
  ###      Dani Cohn - May 17, 2021       ###
  ###     Jennani Jayaram - July 2021     ###
  ###  Based on the Lake Chilwa WB Model  ###
  ## ## ## ## ## ## ## ## ## ## ## ## ## ## #   
  
  ## This is a long script (300+ lines)! ##
  ## If using RStudio, you can use the section finder at the bottom-left of the code window to navigate. ##
  ## Quick location guide:    Priors loop ~line 20    Conversion cms to mm ~line 50    Likelihood function loop ~line 60 ##
  ## (SFD components ~line 135)    Water balance equation (WBE) rolling window ~line 140    Biases ~line 185 ##
  ## Precisions ~line 225     Posterior Predictive WBE ~line 260   Posterior Predictive by component ~line 310 ##
  
  #### For each month ####
  for (j in PosteriorStartMonth:PosteriorEndMonth){
    ##### Priors #####
    ###### Superior ######
    SupPrecip[j]      ~dgamma(PrecipShapeSup[m[j]],PrecipRateSup[m[j]])         # For P, E, R, precision calculated empirically
    SupEvap[j]        ~dnorm(EvapMeanSup[m[j]],EvapPrecisionSup[m[j]])          # For Q_out, precision is probabilistic!
    SupRunoff[j]      =exp(SupLogRunoff[j])
    SupLogRunoff[j]   ~dnorm(LogRunoffMeanSup[m[j]],LogRunoffPrecisionSup[m[j]])
    SupOutflow[j]     ~dnorm(OutflowMeanSup[m[j]],OutflowPrecisionSup[m[j]])
    ####### Add in St. Marys River SFD code here #######
    
    ###### Michigan-Huron ######
    MHPrecip[j]      ~dgamma(PrecipShapeMH[m[j]],PrecipRateMH[m[j]])            # For P, E, R, precision calculated empirically
    MHEvap[j]        ~dnorm(EvapMeanMH[m[j]],EvapPrecisionMH[m[j]])             # For Q_out, precision is probabilistic!
    MHRunoff[j]      =exp(MHLogRunoff[j])
    MHLogRunoff[j]   ~dnorm(LogRunoffMeanMH[m[j]],LogRunoffPrecisionMH[m[j]])
    MHOutflow[j]     ~dnorm(OutflowMeanMH[m[j]],OutflowPrecisionMH[m[j]])
    ####### Add in St. Clair River SFD code here #######
   
    ###### LAKE St. Clair ######
    ## Uses NBS instead of P, E, R - easier calculations
    LkStClNBS[j]        ~dnorm(NBSMeanLkStCl[m[j]],NBSPrecisionLkStCl[m[j]])
    LkStClOutflow[j]     ~dnorm(OutflowMeanLkStCl[m[j]],OutflowPrecisionLkStCl[m[j]])
    ####### Add in Detroit River SFD code here #######
    
    ###### Erie ######
    EriPrecip[j]      ~dgamma(PrecipShapeEri[m[j]],PrecipRateEri[m[j]])         # For P, E, R, precision calculated empirically
    EriEvap[j]        ~dnorm(EvapMeanEri[m[j]],EvapPrecisionEri[m[j]])          # For Q_out, precision is probabilistic!
    EriRunoff[j]      =exp(EriLogRunoff[j])
    EriLogRunoff[j]   ~dnorm(LogRunoffMeanEri[m[j]],LogRunoffPrecisionEri[m[j]])
    EriOutflow[j]     ~dnorm(OutflowMeanEri[m[j]],OutflowPrecisionEri[m[j]])
    ####### Add in Niagara River SFD code here #######
    
    ##### cms to mm conversion ######
    #SupInflow_mm[j]    = 1000*((SupOutflow[j]   *secondsInADay*dayVector[j])/(supArea))
    MHInflow_mm[j]     = 1000*((SupOutflow[j]    *secondsInADay*dayVector[j])/(mhgArea))
    LkStClInflow_mm[j] = 1000*((MHOutflow[j]     *secondsInADay*dayVector[j])/(lkstclArea))
    EriInflow_mm[j]    = 1000*((LkStClOutflow[j] *secondsInADay*dayVector[j])/(eriArea))
    
    SupOutflow_mm[j]    = 1000*((SupOutflow[j]   *secondsInADay*dayVector[j])/(supArea))
    MHOutflow_mm[j]     = 1000*((MHOutflow[j]    *secondsInADay*dayVector[j])/(mhgArea))
    LkStClOutflow_mm[j] = 1000*((LkStClOutflow[j]*secondsInADay*dayVector[j])/(lkstclArea))
    EriOutflow_mm[j]    = 1000*((EriOutflow[j]   *secondsInADay*dayVector[j])/(eriArea))
    
    ##### Likelihood function #####
    ## Adding in additional datasets: use similar format as below, just increase the number ##
    ###### Superior ######
    ySupPrecip1[j]      ~dnorm(ySupPrecip1Mean[j],ySupPrecip1Precision) ## NOAA GLERL GLM HMD ##
    ySupPrecip1Mean[j]  =SupPrecip[j]+ySupPrecip1Bias[m[j]]
    #ySupPrecip2[j]      ~dnorm(ySupPrecip2Mean[j],ySupPrecip2Precision) ## USACE AHPS ##
    #ySupPrecip2Mean[j]  =SupPrecip[j]+ySupPrecip2Bias[m[j]]
    ySupPrecip3[j]      ~dnorm(ySupPrecip3Mean[j],ySupPrecip3Precision) ## USACE AHPS ##
    ySupPrecip3Mean[j]  =SupPrecip[j]+ySupPrecip3Bias[m[j]]
    ySupPrecip4[j]      ~dnorm(ySupPrecip4Mean[j],ySupPrecip4Precision) ## USACE AHPS ##
    ySupPrecip4Mean[j]  =SupPrecip[j]+ySupPrecip4Bias[m[j]]
    ySupPrecip5[j]      ~dnorm(ySupPrecip5Mean[j],ySupPrecip5Precision) ## USACE AHPS ##
    ySupPrecip5Mean[j]  =SupPrecip[j]+ySupPrecip5Bias[m[j]]
    ySupPrecip6[j]      ~dnorm(ySupPrecip6Mean[j],ySupPrecip6Precision) ## USACE AHPS ##
    ySupPrecip6Mean[j]  =SupPrecip[j]+ySupPrecip6Bias[m[j]]
    #ySupPrecip7[j]      ~dnorm(ySupPrecip7Mean[j],ySupPrecip7Precision) ## USACE AHPS ##
    #ySupPrecip7Mean[j]  =SupPrecip[j]+ySupPrecip7Bias[m[j]]
    ySupPrecip8[j]      ~dnorm(ySupPrecip8Mean[j],ySupPrecip8Precision) ## USACE AHPS ##
    ySupPrecip8Mean[j]  =SupPrecip[j]+ySupPrecip8Bias[m[j]]
    ySupPrecip9[j]      ~dnorm(ySupPrecip9Mean[j],ySupPrecip9Precision) ## USACE AHPS ##
    ySupPrecip9Mean[j]  =SupPrecip[j]+ySupPrecip9Bias[m[j]]
    #
    ySupEvap1[j]        ~dnorm(ySupEvap1Mean[j],ySupEvap1Precision) ## NOAA GLERL GLM HMD ##
    ySupEvap1Mean[j]    =SupEvap[j]+ySupEvap1Bias[m[j]]
    #ySupEvap2[j]        ~dnorm(ySupEvap2Mean[j],ySupEvap2Precision) ## USACE AHPS ##
    #ySupEvap2Mean[j]    =SupEvap[j]+ySupEvap2Bias[m[j]]
    ySupEvap3[j]        ~dnorm(ySupEvap3Mean[j],ySupEvap3Precision) ## USACE AHPS ##
    ySupEvap3Mean[j]    =SupEvap[j]+ySupEvap3Bias[m[j]]
    ySupEvap4[j]        ~dnorm(ySupEvap4Mean[j],ySupEvap4Precision) ## USACE AHPS ##
    ySupEvap4Mean[j]    =SupEvap[j]+ySupEvap4Bias[m[j]]
    ySupEvap5[j]        ~dnorm(ySupEvap5Mean[j],ySupEvap5Precision) ## USACE AHPS ##
    ySupEvap5Mean[j]    =SupEvap[j]+ySupEvap5Bias[m[j]]
    #
    ySupRunoff1[j]      ~dnorm(ySupRunoff1Mean[j],ySupRunoff1Precision) ## NOAA GLERL GLM HMD ##
    ySupRunoff1Mean[j]  =SupRunoff[j]+ySupRunoff1Bias[m[j]]
    #ySupRunoff2[j]      ~dnorm(ySupRunoff2Mean[j],ySupRunoff2Precision) ## USACE AHPS ##
    #ySupRunoff2Mean[j]  =SupRunoff[j]+ySupRunoff2Bias[m[j]]
    ySupRunoff3[j]      ~dnorm(ySupRunoff3Mean[j],ySupRunoff3Precision) ## USACE AHPS ##
    ySupRunoff3Mean[j]  =SupRunoff[j]+ySupRunoff3Bias[m[j]]
    ySupRunoff4[j]      ~dnorm(ySupRunoff4Mean[j],ySupRunoff4Precision) ## USACE AHPS ##
    ySupRunoff4Mean[j]  =SupRunoff[j]+ySupRunoff4Bias[m[j]]
    ySupRunoff5[j]      ~dnorm(ySupRunoff5Mean[j],ySupRunoff5Precision) ## USACE AHPS ##
    ySupRunoff5Mean[j]  =SupRunoff[j]+ySupRunoff5Bias[m[j]]
    #
    ySupOutflow1[j]     ~dnorm(ySupOutflow1Mean[j],ySupOutflow1Precision) ## St Marys IGS ##
    ySupOutflow1Mean[j] =SupOutflow[j]+ySupOutflow1Bias[m[j]]
    ySupOutflow2[j]     ~dnorm(ySupOutflow2Mean[j],ySupOutflow2Precision) ## St Marys Flow Accounting ##
    ySupOutflow2Mean[j] =SupOutflow[j]+ySupOutflow2Bias[m[j]]
    #ySupOutflow3[j]     ~dnorm(ySupOutflow3Mean[j],ySupOutflow3Precision) ## St Marys Flow Accounting ##
    #ySupOutflow3Mean[j] =SupOutflow[j]+ySupOutflow3Bias[m[j]]
    
    ###### Michigan-Huron ######
    yMHPrecip1[j]      ~dnorm(yMHPrecip1Mean[j],yMHPrecip1Precision) ## NOAA GLERL GLM HMD ##
    yMHPrecip1Mean[j]  =MHPrecip[j]+yMHPrecip1Bias[m[j]]
    yMHPrecip2[j]      ~dnorm(yMHPrecip2Mean[j],yMHPrecip2Precision) ## USACE AHPS ##
    yMHPrecip2Mean[j]  =MHPrecip[j]+yMHPrecip2Bias[m[j]]
    yMHPrecip3[j]      ~dnorm(yMHPrecip3Mean[j],yMHPrecip3Precision) ## USACE AHPS ##
    yMHPrecip3Mean[j]  =MHPrecip[j]+yMHPrecip3Bias[m[j]]
    yMHPrecip4[j]      ~dnorm(yMHPrecip4Mean[j],yMHPrecip4Precision) ## USACE AHPS ##
    yMHPrecip4Mean[j]  =MHPrecip[j]+yMHPrecip4Bias[m[j]]
    yMHPrecip5[j]      ~dnorm(yMHPrecip5Mean[j],yMHPrecip5Precision) ## USACE AHPS ##
    yMHPrecip5Mean[j]  =MHPrecip[j]+yMHPrecip5Bias[m[j]]
    yMHPrecip6[j]      ~dnorm(yMHPrecip6Mean[j],yMHPrecip6Precision) ## USACE AHPS ##
    yMHPrecip6Mean[j]  =MHPrecip[j]+yMHPrecip6Bias[m[j]]
    yMHPrecip7[j]      ~dnorm(yMHPrecip7Mean[j],yMHPrecip7Precision) ## USACE AHPS ##
    yMHPrecip7Mean[j]  =MHPrecip[j]+yMHPrecip7Bias[m[j]]
    yMHPrecip8[j]      ~dnorm(yMHPrecip8Mean[j],yMHPrecip8Precision) ## USACE AHPS ##
    yMHPrecip8Mean[j]  =MHPrecip[j]+yMHPrecip8Bias[m[j]]
    yMHPrecip9[j]      ~dnorm(yMHPrecip9Mean[j],yMHPrecip9Precision) ## USACE AHPS ##
    yMHPrecip9Mean[j]  =MHPrecip[j]+yMHPrecip9Bias[m[j]]
    #
    yMHEvap1[j]        ~dnorm(yMHEvap1Mean[j],yMHEvap1Precision) ## NOAA GLERL GLM HMD ##
    yMHEvap1Mean[j]    =MHEvap[j]+yMHEvap1Bias[m[j]]
    #yMHEvap2[j]        ~dnorm(yMHEvap2Mean[j],yMHEvap2Precision) ## USACE AHPS ##
    #yMHEvap2Mean[j]    =MHEvap[j]+yMHEvap2Bias[m[j]]
    yMHEvap3[j]        ~dnorm(yMHEvap3Mean[j],yMHEvap3Precision) ## USACE AHPS ##
    yMHEvap3Mean[j]    =MHEvap[j]+yMHEvap3Bias[m[j]]
    yMHEvap4[j]        ~dnorm(yMHEvap4Mean[j],yMHEvap4Precision) ## USACE AHPS ##
    yMHEvap4Mean[j]    =MHEvap[j]+yMHEvap4Bias[m[j]]
    yMHEvap5[j]        ~dnorm(yMHEvap5Mean[j],yMHEvap5Precision) ## USACE AHPS ##
    yMHEvap5Mean[j]    =MHEvap[j]+yMHEvap5Bias[m[j]]
    #
    yMHRunoff1[j]      ~dnorm(yMHRunoff1Mean[j],yMHRunoff1Precision) ## NOAA GLERL GLM HMD ##
    yMHRunoff1Mean[j]  =MHRunoff[j]+yMHRunoff1Bias[m[j]]
    #yMHRunoff2[j]      ~dnorm(yMHRunoff2Mean[j],yMHRunoff2Precision) ## USACE AHPS ##
    #yMHRunoff2Mean[j]  =MHRunoff[j]+yMHRunoff2Bias[m[j]]
    yMHRunoff3[j]      ~dnorm(yMHRunoff3Mean[j],yMHRunoff3Precision) ## USACE AHPS ##
    yMHRunoff3Mean[j]  =MHRunoff[j]+yMHRunoff3Bias[m[j]]
    yMHRunoff4[j]      ~dnorm(yMHRunoff4Mean[j],yMHRunoff4Precision) ## USACE AHPS ##
    yMHRunoff4Mean[j]  =MHRunoff[j]+yMHRunoff4Bias[m[j]]
    yMHRunoff5[j]      ~dnorm(yMHRunoff5Mean[j],yMHRunoff5Precision) ## USACE AHPS ##
    yMHRunoff5Mean[j]  =MHRunoff[j]+yMHRunoff5Bias[m[j]]
    #
    yMHOutflow1[j]     ~dnorm(yMHOutflow1Mean[j],yMHOutflow1Precision) ## St Clair IGS ##
    yMHOutflow1Mean[j] =MHOutflow[j]+yMHOutflow1Bias[m[j]]
    yMHOutflow2[j]     ~dnorm(yMHOutflow2Mean[j],yMHOutflow2Precision) ## St Clair SFD.ADVM ##
    yMHOutflow2Mean[j] =MHOutflow[j]+yMHOutflow2Bias[m[j]]
    #yMHOutflow3[j]     ~dnorm(yMHOutflow3Mean[j],yMHOutflow3Precision) ## St Clair SFD.ADVM ##
    #yMHOutflow3Mean[j] =MHOutflow[j]+yMHOutflow3Bias[m[j]]
    
    ###### LAKE St. Clair ######
    ## Uses NBS instead of P, E, R - easier calculations
    yLkStClNBS1[j]      ~dnorm(yLkStClNBS1Mean[j],yLkStClNBS1Precision) ## NOAA GLERL GLM HMD ##
    yLkStClNBS1Mean[j]  =(LkStClNBS[j]+yLkStClNBS1Bias[m[j]]) #*Mult1Bias[m[j]]
    #yLkStClNBS2[j]      ~dnorm(yLkStClNBS2Mean[j],yLkStClNBS2Precision) ## USACE AHPS ##
    #yLkStClNBS2Mean[j]  =(LkStClNBS[j]+yLkStClNBS2Bias[m[j]]) #*Mult2Bias[m[j]]
    yLkStClNBS3[j]      ~dnorm(yLkStClNBS3Mean[j],yLkStClNBS3Precision) ## USACE AHPS ##
    yLkStClNBS3Mean[j]  =(LkStClNBS[j]+yLkStClNBS3Bias[m[j]]) #*Mult2Bias[m[j]]
    yLkStClNBS4[j]      ~dnorm(yLkStClNBS4Mean[j],yLkStClNBS4Precision) ## USACE AHPS ##
    yLkStClNBS4Mean[j]  =(LkStClNBS[j]+yLkStClNBS4Bias[m[j]]) #*Mult2Bias[m[j]]
    yLkStClNBS5[j]      ~dnorm(yLkStClNBS5Mean[j],yLkStClNBS5Precision) ## USACE AHPS ##
    yLkStClNBS5Mean[j]  =(LkStClNBS[j]+yLkStClNBS5Bias[m[j]]) #*Mult2Bias[m[j]]
    yLkStClNBS6[j]      ~dnorm(yLkStClNBS6Mean[j],yLkStClNBS6Precision) ## USACE AHPS ##
    yLkStClNBS6Mean[j]  =(LkStClNBS[j]+yLkStClNBS6Bias[m[j]]) #*Mult2Bias[m[j]]
    #
#sk     yLkStClOutflow1[j]     ~dnorm(yLkStClOutflow1Mean[j],yLkStClOutflow1Precision) ## Detroit IGS ##
#sk     yLkStClOutflow1Mean[j] =LkStClOutflow[j]+yLkStClOutflow1Bias[m[j]]
#sk     yLkStClOutflow2[j]     ~dnorm(yLkStClOutflow2Mean[j],yLkStClOutflow2Precision) ## Detroit SFD.ADVM ##
#sk    yLkStClOutflow2Mean[j] =LkStClOutflow[j]+yLkStClOutflow2Bias[m[j]]
    #yLkStClOutflow3[j]     ~dnorm(yLkStClOutflow3Mean[j],yLkStClOutflow3Precision) ## Detroit SFD.ADVM ##
    #yLkStClOutflow3Mean[j] =LkStClOutflow[j]+yLkStClOutflow3Bias[m[j]]
    
    ###### Erie ######
    yEriPrecip1[j]      ~dnorm(yEriPrecip1Mean[j],yEriPrecip1Precision) ## NOAA GLERL GLM HMD ##
    yEriPrecip1Mean[j]  =EriPrecip[j]+yEriPrecip1Bias[m[j]]
    #yEriPrecip2[j]      ~dnorm(yEriPrecip2Mean[j],yEriPrecip2Precision) ## USACE AHPS ##
    #yEriPrecip2Mean[j]  =EriPrecip[j]+yEriPrecip2Bias[m[j]]
    yEriPrecip3[j]      ~dnorm(yEriPrecip3Mean[j],yEriPrecip3Precision) 
    yEriPrecip3Mean[j]  =EriPrecip[j]+yEriPrecip3Bias[m[j]]
    yEriPrecip4[j]      ~dnorm(yEriPrecip4Mean[j],yEriPrecip4Precision) 
    yEriPrecip4Mean[j]  =EriPrecip[j]+yEriPrecip4Bias[m[j]]
    yEriPrecip5[j]      ~dnorm(yEriPrecip5Mean[j],yEriPrecip5Precision) 
    yEriPrecip5Mean[j]  =EriPrecip[j]+yEriPrecip5Bias[m[j]]
    yEriPrecip6[j]      ~dnorm(yEriPrecip6Mean[j],yEriPrecip6Precision) 
    yEriPrecip6Mean[j]  =EriPrecip[j]+yEriPrecip6Bias[m[j]]
    yEriPrecip7[j]      ~dnorm(yEriPrecip7Mean[j],yEriPrecip7Precision) 
    yEriPrecip7Mean[j]  =EriPrecip[j]+yEriPrecip7Bias[m[j]]
    yEriPrecip8[j]      ~dnorm(yEriPrecip8Mean[j],yEriPrecip8Precision) 
    yEriPrecip8Mean[j]  =EriPrecip[j]+yEriPrecip8Bias[m[j]]
    yEriPrecip9[j]      ~dnorm(yEriPrecip9Mean[j],yEriPrecip9Precision)
    yEriPrecip9Mean[j]  =EriPrecip[j]+yEriPrecip9Bias[m[j]]
    #
    yEriEvap1[j]        ~dnorm(yEriEvap1Mean[j],yEriEvap1Precision) ## NOAA GLERL GLM HMD ##
    yEriEvap1Mean[j]    =EriEvap[j]+yEriEvap1Bias[m[j]]
    #yEriEvap2[j]        ~dnorm(yEriEvap2Mean[j],yEriEvap2Precision) ## USACE AHPS ##
    #yEriEvap2Mean[j]    =EriEvap[j]+yEriEvap2Bias[m[j]]
    yEriEvap3[j]        ~dnorm(yEriEvap3Mean[j],yEriEvap3Precision) 
    yEriEvap3Mean[j]    =EriEvap[j]+yEriEvap3Bias[m[j]]
    yEriEvap4[j]        ~dnorm(yEriEvap4Mean[j],yEriEvap4Precision) 
    yEriEvap4Mean[j]    =EriEvap[j]+yEriEvap4Bias[m[j]]
    yEriEvap5[j]        ~dnorm(yEriEvap5Mean[j],yEriEvap5Precision) 
    yEriEvap5Mean[j]    =EriEvap[j]+yEriEvap5Bias[m[j]]
    #
    yEriRunoff1[j]      ~dnorm(yEriRunoff1Mean[j],yEriRunoff1Precision) ## NOAA GLERL GLM HMD ##
    yEriRunoff1Mean[j]  =EriRunoff[j]+yEriRunoff1Bias[m[j]]
    #yEriRunoff2[j]      ~dnorm(yEriRunoff2Mean[j],yEriRunoff2Precision) ## USACE AHPS ##
    #yEriRunoff2Mean[j]  =EriRunoff[j]+yEriRunoff2Bias[m[j]]
    yEriRunoff3[j]      ~dnorm(yEriRunoff3Mean[j],yEriRunoff3Precision) 
    yEriRunoff3Mean[j]  =EriRunoff[j]+yEriRunoff3Bias[m[j]]
    yEriRunoff4[j]      ~dnorm(yEriRunoff4Mean[j],yEriRunoff4Precision) 
    yEriRunoff4Mean[j]  =EriRunoff[j]+yEriRunoff4Bias[m[j]]
    yEriRunoff5[j]      ~dnorm(yEriRunoff5Mean[j],yEriRunoff5Precision) 
    yEriRunoff5Mean[j]  =EriRunoff[j]+yEriRunoff5Bias[m[j]]
    #
    yEriOutflow1[j]     ~dnorm(yEriOutflow1Mean[j],yEriOutflow1Precision) ## Niagara SFD.ADVM ##
    yEriOutflow1Mean[j] =EriOutflow[j]+yEriOutflow1Bias[m[j]]
    #yEriOutflow2[j]     ~dnorm(yEriOutflow2Mean[j],yEriOutflow2Precision) ## Niagara Coordinated (exact same values as SFD.ADVM) ##
    #yEriOutflow2Mean[j] =EriOutflow[j]+yEriOutflow2Bias[m[j]]
    #yEriOutflow3[j]     ~dnorm(yEriOutflow3Mean[j],yEriOutflow3Precision) ## Niagara Coordinated (exact same values as SFD.ADVM) ##
    #yEriOutflow3Mean[j] =EriOutflow[j]+yEriOutflow3Bias[m[j]]
  }
  
  #### For rolling window and Water balance Equation ####
  for (k in RollPeriod:PosteriorEndMonth){ 
    ##### Superior #####
    ## WB Equation ##
    ySupRStore[k]~dnorm(SupRStore[k],ySupRStorePreci)
    SupRStore[k]=(sum(SupPrecip[(k-RollPeriod+1):k])
                 -sum(SupEvap[(k-RollPeriod+1):k])
                 +sum(SupRunoff[(k-RollPeriod+1):k])
                 +sum(SupError[m[(k-RollPeriod+1):k]])
                 -sum(SupOutflow_mm[(k-RollPeriod+1):k]))
    ## Added parentheses here to definitively contain SupRStore equation ##
    ##### Michigan-Huron #####
    ## WB Equation ##
    yMHRStore[k]~dnorm(MHRStore[k],yMHRStorePreci)
    MHRStore[k]=(sum(MHPrecip[(k-RollPeriod+1):k])
      -sum(MHEvap[(k-RollPeriod+1):k])
      +sum(MHRunoff[(k-RollPeriod+1):k])
      +sum(MHError[m[(k-RollPeriod+1):k]])
      -sum(MHOutflow_mm[(k-RollPeriod+1):k])
      #+sum(SupOutflow_mm[(k-RollPeriod+1):k])#JLC
      +sum(MHInflow_mm[(k-RollPeriod+1):k])) ## MH Inflow = Sup Outflow ##
    ## Added parentheses here to definitively contain MHRStore equation ##
    ##### LAKE St. Clair #####
    ## WB Equation ##
    yLkStClRStore[k]~dnorm(LkStClRStore[k],yLkStClRStorePreci)
    LkStClRStore[k]=(sum(LkStClNBS[(k-RollPeriod+1):k])
                  +sum(LkStClError[m[(k-RollPeriod+1):k]])
                  -sum(LkStClOutflow_mm[(k-RollPeriod+1):k])
                  #+sum(MHOutflow_mm[(k-RollPeriod+1):k])#JLC
                  +sum(LkStClInflow_mm[(k-RollPeriod+1):k])) ## LkStCl Inflow = MH Outflow ##
    ## Added parentheses here to definitively contain LkStClRStore equation ##
    ##### Erie #####
    ## WB Equation ##
    yEriRStore[k]~dnorm(EriRStore[k],yEriRStorePreci)
    EriRStore[k]=(sum(EriPrecip[(k-RollPeriod+1):k])
                 -sum(EriEvap[(k-RollPeriod+1):k])
                 +sum(EriRunoff[(k-RollPeriod+1):k])
                 +sum(EriError[m[(k-RollPeriod+1):k]])
                 -sum(EriOutflow_mm[(k-RollPeriod+1):k])
                 #+sum(LkStClOutflow_mm[(k-RollPeriod+1):k])#JLC
                 +sum(EriInflow_mm[(k-RollPeriod+1):k])) ## Eri Inflow = LkStCl Outflow ##
    ## Added parentheses here to definitively contain EriRStore equation ##
  }
  
  #### Bias ####
  for (jp in 1:12){
    #Mult1Bias[jp] ~ dnorm(1,0.1)
    #Mult2Bias[jp] ~ dnorm(1,0.1)
    ## Adding in additional datasets: use similar format as below, just increase the number ##
    ##### Superior #####
    ySupPrecip1Bias[jp]      ~dnorm(0,0.01)        
#    ySupPrecip2Bias[jp]      ~dnorm(0,0.01)
    ySupPrecip3Bias[jp]      ~dnorm(0,0.01)
    ySupPrecip4Bias[jp]      ~dnorm(0,0.01)
    ySupPrecip5Bias[jp]      ~dnorm(0,0.01)
    ySupPrecip6Bias[jp]      ~dnorm(0,0.01)
#    ySupPrecip7Bias[jp]      ~dnorm(0,0.01)
    ySupPrecip8Bias[jp]      ~dnorm(0,0.01)
    ySupPrecip9Bias[jp]      ~dnorm(0,0.01)
    ySupEvap1Bias[jp]        ~dnorm(0,0.01)
#    ySupEvap2Bias[jp]        ~dnorm(0,0.01)
    ySupEvap3Bias[jp]        ~dnorm(0,0.01)
    ySupEvap4Bias[jp]        ~dnorm(0,0.01)
    ySupEvap5Bias[jp]        ~dnorm(0,0.01)
    ySupRunoff1Bias[jp]      ~dnorm(0,0.01)
#    ySupRunoff2Bias[jp]      ~dnorm(0,0.01)
    ySupRunoff3Bias[jp]      ~dnorm(0,0.01)
    ySupRunoff4Bias[jp]      ~dnorm(0,0.01)
    ySupRunoff5Bias[jp]      ~dnorm(0,0.01)
    ySupOutflow1Bias[jp]     ~dnorm(0,0.01)
    ySupOutflow2Bias[jp]     ~dnorm(0,0.01)
#    ySupOutflow3Bias[jp]     ~dnorm(0,0.01)
    SupError[jp]             ~dnorm(0,0.01)
    ##### Michigan-Huron #####
    yMHPrecip1Bias[jp]      ~dnorm(0,0.01)        
    yMHPrecip2Bias[jp]      ~dnorm(0,0.01)
    yMHPrecip3Bias[jp]      ~dnorm(0,0.01)
    yMHPrecip4Bias[jp]      ~dnorm(0,0.01)
    yMHPrecip5Bias[jp]      ~dnorm(0,0.01)
    yMHPrecip6Bias[jp]      ~dnorm(0,0.01)
    yMHPrecip7Bias[jp]      ~dnorm(0,0.01)
    yMHPrecip8Bias[jp]      ~dnorm(0,0.01)
    yMHPrecip9Bias[jp]      ~dnorm(0,0.01)
    yMHEvap1Bias[jp]        ~dnorm(0,0.01)
#    yMHEvap2Bias[jp]        ~dnorm(0,0.01)
    yMHEvap3Bias[jp]        ~dnorm(0,0.01)
    yMHEvap4Bias[jp]        ~dnorm(0,0.01)
    yMHEvap5Bias[jp]        ~dnorm(0,0.01)
    yMHRunoff1Bias[jp]      ~dnorm(0,0.01)
#    yMHRunoff2Bias[jp]      ~dnorm(0,0.01)
    yMHRunoff3Bias[jp]      ~dnorm(0,0.01)
    yMHRunoff4Bias[jp]      ~dnorm(0,0.01)
    yMHRunoff5Bias[jp]      ~dnorm(0,0.01)
    yMHOutflow1Bias[jp]     ~dnorm(0,0.01)
    yMHOutflow2Bias[jp]     ~dnorm(0,0.01)
#    yMHOutflow3Bias[jp]     ~dnorm(0,0.01)
    MHError[jp]             ~dnorm(0,0.01)
    ##### LAKE St. Clair #####
    yLkStClNBS1Bias[jp]     ~dnorm(0,0.01)       
#    yLkStClNBS2Bias[jp]     ~dnorm(0,0.01) 
    yLkStClNBS3Bias[jp]     ~dnorm(0,0.01) 
    yLkStClNBS4Bias[jp]     ~dnorm(0,0.01) 
    yLkStClNBS5Bias[jp]     ~dnorm(0,0.01) 
    yLkStClNBS6Bias[jp]     ~dnorm(0,0.01) 
#sk     yLkStClOutflow1Bias[jp] ~dnorm(0,0.01)
#sk    yLkStClOutflow2Bias[jp] ~dnorm(0,0.01)
#    yLkStClOutflow3Bias[jp] ~dnorm(0,0.01)
    LkStClError[jp]         ~dnorm(0,0.01)
    ##### Erie #####
    yEriPrecip1Bias[jp]     ~dnorm(0,0.01)        
#    yEriPrecip2Bias[jp]     ~dnorm(0,0.01)
    yEriPrecip3Bias[jp]     ~dnorm(0,0.01)
    yEriPrecip4Bias[jp]     ~dnorm(0,0.01)
    yEriPrecip5Bias[jp]     ~dnorm(0,0.01)
    yEriPrecip6Bias[jp]     ~dnorm(0,0.01)
    yEriPrecip7Bias[jp]     ~dnorm(0,0.01)
    yEriPrecip8Bias[jp]     ~dnorm(0,0.01)
    yEriPrecip9Bias[jp]     ~dnorm(0,0.01)
    yEriEvap1Bias[jp]       ~dnorm(0,0.01)
#    yEriEvap2Bias[jp]       ~dnorm(0,0.01)
    yEriEvap3Bias[jp]       ~dnorm(0,0.01)
    yEriEvap4Bias[jp]       ~dnorm(0,0.01)
    yEriEvap5Bias[jp]       ~dnorm(0,0.01)
    yEriRunoff1Bias[jp]     ~dnorm(0,0.01)
#    yEriRunoff2Bias[jp]     ~dnorm(0,0.01)
    yEriRunoff3Bias[jp]     ~dnorm(0,0.01)
    yEriRunoff4Bias[jp]     ~dnorm(0,0.01)
    yEriRunoff5Bias[jp]     ~dnorm(0,0.01)
    yEriOutflow1Bias[jp]    ~dnorm(0,0.01)
#    yEriOutflow2Bias[jp]    ~dnorm(0,0.01)
#    yEriOutflow3Bias[jp]    ~dnorm(0,0.01)
    EriError[jp]            ~dnorm(0,0.01)
  }
  
  #### Precision for Observation ####
  ## Adding in additional datasets: use similar format as below, just increase the number ##
  ##### Superior #####
  ySupRStorePreci            ~dlnorm(-log(100),1)#dgamma(0.01,0.01)
  ySupPrecip1Precision       ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
#  ySupPrecip2Precision       ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  ySupPrecip3Precision       ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  ySupPrecip4Precision       ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  ySupPrecip5Precision       ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  ySupPrecip6Precision       ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
#  ySupPrecip7Precision       ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  ySupPrecip8Precision       ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  ySupPrecip9Precision       ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  ySupEvap1Precision         ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
#  ySupEvap2Precision         ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  ySupEvap3Precision         ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  ySupEvap4Precision         ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  ySupEvap5Precision         ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  ySupRunoff1Precision       ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
#  ySupRunoff2Precision       ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  ySupRunoff3Precision       ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  ySupRunoff4Precision       ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  ySupRunoff5Precision       ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  ySupOutflow1Precision      ~dlnorm(-log(10000),1)#dgamma(0.1,0.1)
  ySupOutflow2Precision      ~dlnorm(-log(10000),1)#dgamma(0.1,0.1)
#  ySupOutflow3Precision      ~dlnorm(-log(10000),1)#dgamma(0.1,0.1)
  ##### Michigan-Huron #####
  yMHRStorePreci            ~dlnorm(-log(100),1)#dgamma(0.01,0.01)
  yMHPrecip1Precision       ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  yMHPrecip2Precision       ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  yMHPrecip3Precision       ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  yMHPrecip4Precision       ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  yMHPrecip5Precision       ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  yMHPrecip6Precision       ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  yMHPrecip7Precision       ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  yMHPrecip8Precision       ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  yMHPrecip9Precision       ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  yMHEvap1Precision         ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
#  yMHEvap2Precision         ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  yMHEvap3Precision         ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  yMHEvap4Precision         ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  yMHEvap5Precision         ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  yMHRunoff1Precision       ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
#  yMHRunoff2Precision       ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  yMHRunoff3Precision       ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  yMHRunoff4Precision       ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  yMHRunoff5Precision       ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  yMHOutflow1Precision      ~dlnorm(-log(10000),1)#dgamma(0.1,0.1)
   yMHOutflow2Precision      ~dlnorm(-log(10000),1)#dgamma(0.1,0.1)
#  yMHOutflow3Precision      ~dlnorm(-log(10000),1)#dgamma(0.1,0.1)
  ##### LAKE St. Clair #####
  yLkStClRStorePreci        ~dlnorm(-log(100),1)#dgamma(0.01,0.01)
  yLkStClNBS1Precision      ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
#  yLkStClNBS2Precision      ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  yLkStClNBS3Precision      ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  yLkStClNBS4Precision      ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  yLkStClNBS5Precision      ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  yLkStClNBS6Precision      ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
#sk   yLkStClOutflow1Precision  ~dlnorm(-log(10000),1)#dgamma(0.1,0.1)
#sk  yLkStClOutflow2Precision  ~dlnorm(-log(10000),1)#dgamma(0.1,0.1)
#  yLkStClOutflow3Precision  ~dlnorm(-log(10000),1)#dgamma(0.1,0.1)
  ##### Erie #####
  yEriRStorePreci           ~dlnorm(-log(100),1)#dgamma(0.01,0.01)
  yEriPrecip1Precision      ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
#  yEriPrecip2Precision      ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  yEriPrecip3Precision      ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  yEriPrecip4Precision      ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  yEriPrecip5Precision      ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  yEriPrecip6Precision      ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  yEriPrecip7Precision      ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  yEriPrecip8Precision      ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  yEriPrecip9Precision      ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  yEriEvap1Precision        ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
#  yEriEvap2Precision        ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  yEriEvap3Precision        ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  yEriEvap4Precision        ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  yEriEvap5Precision        ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  yEriRunoff1Precision      ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
#  yEriRunoff2Precision      ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  yEriRunoff3Precision      ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  yEriRunoff4Precision      ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  yEriRunoff5Precision      ~dlnorm(-log(100),1)#dgamma(0.1,0.1)
  yEriOutflow1Precision     ~dlnorm(-log(10000),1)#dgamma(0.1,0.1)
#  yEriOutflow2Precision     ~dlnorm(-log(10000),1)#dgamma(0.1,0.1)
#  yEriOutflow3Precision     ~dlnorm(-log(10000),1)#dgamma(0.1,0.1)
  
  #### Posterior Predictive (PP) Distribution ####
  for (jp in PosteriorStartMonth:PosteriorEndMonth){
    ##### Lake level in each month #####
    ###### Superior ######
    ## WB Equation ##
    ySupDStorePP[jp]~dnorm(SupDStore[jp],ySupRStorePreci)
    SupDStore[jp]=(
      SupPrecip[jp]
      -SupEvap[jp]
      +SupRunoff[jp]
      -SupOutflowPP_mm[jp]
      +SupError[m[jp]]
    )
    ###### Michigan-Huron ######
    ## WB Equation ##
    yMHDStorePP[jp]~dnorm(MHDStore[jp],yMHRStorePreci)
    MHDStore[jp]=(
      MHPrecip[jp]
      -MHEvap[jp]
      +MHRunoff[jp]
      -MHOutflowPP_mm[jp]
      +MHInflowPP_mm[jp]#JLC
      #+SupOutflowPP_mm[jp] ## MH Inflow = Sup Outflow ##
      +MHError[m[jp]]
    )
    ###### LAKE St. Clair ######
    ## WB Equation ##
    yLkStClDStorePP[jp]~dnorm(LkStClDStore[jp],yLkStClRStorePreci)
    LkStClDStore[jp]=(
      LkStClNBS[jp]
      -LkStClOutflowPP_mm[jp]
      +LkStClInflowPP_mm[jp]#JLC
      #+MHOutflowPP_mm[jp] ## LkStCl Inflow = MH Outflow ##
      +LkStClError[m[jp]]
    )
    ###### Erie ######
    ## WB Equation ##
    yEriDStorePP[jp]~dnorm(EriDStore[jp],yEriRStorePreci)
    EriDStore[jp]=(
      EriPrecip[jp]
      -EriEvap[jp]
      +EriRunoff[jp]
      -EriOutflowPP_mm[jp]
      +EriInflowPP_mm[jp]#JLC
      #+LkStClOutflowPP_mm[jp]#+MHOutflowPP_mm[jp] ## Eri Inflow = MH Outflow ##  --JLC error detected?
      +EriError[m[jp]]
    )
    #SupInflowPP_mm[jp]    = 1000*((SupOutflow[jp]   *secondsInADay*dayVector[jp])/(supArea))
    MHInflowPP_mm[jp]     = 1000*((SupOutflow[jp]    *secondsInADay*dayVector[jp])/mhgArea)
    LkStClInflowPP_mm[jp] = 1000*((MHOutflow[jp]     *secondsInADay*dayVector[jp])/lkstclArea)
    EriInflowPP_mm[jp]    = 1000*((LkStClOutflow[jp] *secondsInADay*dayVector[jp])/eriArea)
    
    ## Added by JLC
    SupOutflowPP_mm[jp]    = 1000*((SupOutflow[jp]   *secondsInADay*dayVector[jp])/(supArea))
    MHOutflowPP_mm[jp]     = 1000*((MHOutflow[jp]    *secondsInADay*dayVector[jp])/mhgArea)
    LkStClOutflowPP_mm[jp] = 1000*((LkStClOutflow[jp]*secondsInADay*dayVector[jp])/lkstclArea)
    EriOutflowPP_mm[jp]    = 1000*((EriOutflow[jp]   *secondsInADay*dayVector[jp])/eriArea)
  }
    
    ##### Setting PP distributions by component #####
    ## Adding in additional datasets: use similar format as below, just increase the number ##
    for (jp in PosteriorStartMonth:PosteriorEndMonth){
      ###### Superior ######
      ySupPrecip1PP[jp]    ~dnorm(ySupPrecip1Mean[jp],ySupPrecip1Precision)
#      ySupPrecip2PP[jp]    ~dnorm(ySupPrecip2Mean[jp],ySupPrecip2Precision)
      ySupPrecip3PP[jp]    ~dnorm(ySupPrecip3Mean[jp],ySupPrecip3Precision)
      ySupPrecip4PP[jp]    ~dnorm(ySupPrecip4Mean[jp],ySupPrecip4Precision)
      ySupPrecip5PP[jp]    ~dnorm(ySupPrecip5Mean[jp],ySupPrecip5Precision)
      ySupPrecip6PP[jp]    ~dnorm(ySupPrecip6Mean[jp],ySupPrecip6Precision)
#      ySupPrecip7PP[jp]    ~dnorm(ySupPrecip7Mean[jp],ySupPrecip7Precision)
      ySupPrecip8PP[jp]    ~dnorm(ySupPrecip8Mean[jp],ySupPrecip8Precision)
      ySupPrecip9PP[jp]    ~dnorm(ySupPrecip9Mean[jp],ySupPrecip9Precision)
      ySupEvap1PP[jp]      ~dnorm(ySupEvap1Mean[jp],ySupEvap1Precision)
#      ySupEvap2PP[jp]      ~dnorm(ySupEvap2Mean[jp],ySupEvap2Precision)
      ySupEvap3PP[jp]      ~dnorm(ySupEvap3Mean[jp],ySupEvap3Precision)
      ySupEvap4PP[jp]      ~dnorm(ySupEvap4Mean[jp],ySupEvap4Precision)
      ySupEvap5PP[jp]      ~dnorm(ySupEvap5Mean[jp],ySupEvap5Precision)
      ySupRunoff1PP[jp]    ~dnorm(ySupRunoff1Mean[jp],ySupRunoff1Precision)
#      ySupRunoff2PP[jp]    ~dnorm(ySupRunoff2Mean[jp],ySupRunoff2Precision)
      ySupRunoff3PP[jp]    ~dnorm(ySupRunoff3Mean[jp],ySupRunoff3Precision)
      ySupRunoff4PP[jp]    ~dnorm(ySupRunoff4Mean[jp],ySupRunoff4Precision)
      ySupRunoff5PP[jp]    ~dnorm(ySupRunoff5Mean[jp],ySupRunoff5Precision)
      ySupOutflow1PP[jp]   ~dnorm(ySupOutflow1Mean[jp],ySupOutflow1Precision)
      ySupOutflow2PP[jp]   ~dnorm(ySupOutflow2Mean[jp],ySupOutflow2Precision)
#      ySupOutflow3PP[jp]   ~dnorm(ySupOutflow3Mean[jp],ySupOutflow3Precision)
      ###### Michigan-Huron ######
      yMHPrecip1PP[jp]     ~dnorm(yMHPrecip1Mean[jp],yMHPrecip1Precision)
      yMHPrecip2PP[jp]     ~dnorm(yMHPrecip2Mean[jp],yMHPrecip2Precision)
      yMHPrecip3PP[jp]     ~dnorm(yMHPrecip3Mean[jp],yMHPrecip3Precision)
      yMHPrecip4PP[jp]     ~dnorm(yMHPrecip4Mean[jp],yMHPrecip4Precision)
      yMHPrecip5PP[jp]     ~dnorm(yMHPrecip5Mean[jp],yMHPrecip5Precision)
      yMHPrecip6PP[jp]     ~dnorm(yMHPrecip6Mean[jp],yMHPrecip6Precision)
      yMHPrecip7PP[jp]     ~dnorm(yMHPrecip7Mean[jp],yMHPrecip7Precision)
      yMHPrecip8PP[jp]     ~dnorm(yMHPrecip8Mean[jp],yMHPrecip8Precision)
      yMHPrecip9PP[jp]     ~dnorm(yMHPrecip9Mean[jp],yMHPrecip9Precision)
      yMHEvap1PP[jp]       ~dnorm(yMHEvap1Mean[jp],yMHEvap1Precision)
#      yMHEvap2PP[jp]       ~dnorm(yMHEvap2Mean[jp],yMHEvap2Precision)
      yMHEvap3PP[jp]       ~dnorm(yMHEvap3Mean[jp],yMHEvap3Precision)
      yMHEvap4PP[jp]       ~dnorm(yMHEvap4Mean[jp],yMHEvap4Precision)
      yMHEvap5PP[jp]       ~dnorm(yMHEvap5Mean[jp],yMHEvap5Precision)
      yMHRunoff1PP[jp]     ~dnorm(yMHRunoff1Mean[jp],yMHRunoff1Precision)
#      yMHRunoff2PP[jp]     ~dnorm(yMHRunoff2Mean[jp],yMHRunoff2Precision)
      yMHRunoff3PP[jp]     ~dnorm(yMHRunoff3Mean[jp],yMHRunoff3Precision)
      yMHRunoff4PP[jp]     ~dnorm(yMHRunoff4Mean[jp],yMHRunoff4Precision)
      yMHRunoff5PP[jp]     ~dnorm(yMHRunoff5Mean[jp],yMHRunoff5Precision)
      yMHOutflow1PP[jp]    ~dnorm(yMHOutflow1Mean[jp],yMHOutflow1Precision)
      yMHOutflow2PP[jp]    ~dnorm(yMHOutflow2Mean[jp],yMHOutflow2Precision)
#      yMHOutflow3PP[jp]    ~dnorm(yMHOutflow3Mean[jp],yMHOutflow3Precision)
      ###### LAKE St. Clair ######
      yLkStClNBS1PP[jp] ~dnorm(yLkStClNBS1Mean[jp],yLkStClNBS1Precision)
#      yLkStClNBS2PP[jp] ~dnorm(yLkStClNBS2Mean[jp],yLkStClNBS2Precision)
      yLkStClNBS3PP[jp] ~dnorm(yLkStClNBS3Mean[jp],yLkStClNBS3Precision)
      yLkStClNBS4PP[jp] ~dnorm(yLkStClNBS4Mean[jp],yLkStClNBS4Precision)
      yLkStClNBS5PP[jp] ~dnorm(yLkStClNBS5Mean[jp],yLkStClNBS5Precision)
      yLkStClNBS6PP[jp] ~dnorm(yLkStClNBS6Mean[jp],yLkStClNBS6Precision)
#sk       yLkStClOutflow1PP[jp]~dnorm(yLkStClOutflow1Mean[jp],yLkStClOutflow1Precision)
#sk      yLkStClOutflow2PP[jp]~dnorm(yLkStClOutflow2Mean[jp],yLkStClOutflow2Precision)
#      yLkStClOutflow3PP[jp]~dnorm(yLkStClOutflow3Mean[jp],yLkStClOutflow3Precision)
      ###### Erie ######
      yEriPrecip1PP[jp]    ~dnorm(yEriPrecip1Mean[jp],yEriPrecip1Precision)
#      yEriPrecip2PP[jp]    ~dnorm(yEriPrecip2Mean[jp],yEriPrecip2Precision)
      yEriPrecip3PP[jp]    ~dnorm(yEriPrecip3Mean[jp],yEriPrecip3Precision)
      yEriPrecip4PP[jp]    ~dnorm(yEriPrecip4Mean[jp],yEriPrecip4Precision)
      yEriPrecip5PP[jp]    ~dnorm(yEriPrecip5Mean[jp],yEriPrecip5Precision)
      yEriPrecip6PP[jp]    ~dnorm(yEriPrecip6Mean[jp],yEriPrecip6Precision)
      yEriPrecip7PP[jp]    ~dnorm(yEriPrecip7Mean[jp],yEriPrecip7Precision)
      yEriPrecip8PP[jp]    ~dnorm(yEriPrecip8Mean[jp],yEriPrecip8Precision)
      yEriPrecip9PP[jp]    ~dnorm(yEriPrecip9Mean[jp],yEriPrecip9Precision)
      yEriEvap1PP[jp]      ~dnorm(yEriEvap1Mean[jp],yEriEvap1Precision)
#      yEriEvap2PP[jp]      ~dnorm(yEriEvap2Mean[jp],yEriEvap2Precision)
      yEriEvap3PP[jp]      ~dnorm(yEriEvap3Mean[jp],yEriEvap3Precision)
      yEriEvap4PP[jp]      ~dnorm(yEriEvap4Mean[jp],yEriEvap4Precision)
      yEriEvap5PP[jp]      ~dnorm(yEriEvap5Mean[jp],yEriEvap5Precision)
      yEriRunoff1PP[jp]    ~dnorm(yEriRunoff1Mean[jp],yEriRunoff1Precision)
#      yEriRunoff2PP[jp]    ~dnorm(yEriRunoff2Mean[jp],yEriRunoff2Precision)
      yEriRunoff3PP[jp]    ~dnorm(yEriRunoff3Mean[jp],yEriRunoff3Precision)
      yEriRunoff4PP[jp]    ~dnorm(yEriRunoff4Mean[jp],yEriRunoff4Precision)
      yEriRunoff5PP[jp]    ~dnorm(yEriRunoff5Mean[jp],yEriRunoff5Precision)
      yEriOutflow1PP[jp]   ~dnorm(yEriOutflow1Mean[jp],yEriOutflow1Precision)
#      yEriOutflow2PP[jp]   ~dnorm(yEriOutflow2Mean[jp],yEriOutflow2Precision)
#      yEriOutflow3PP[jp]   ~dnorm(yEriOutflow3Mean[jp],yEriOutflow3Precision)
    }
    
}