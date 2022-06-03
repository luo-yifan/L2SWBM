## ## ## ## ## ## ## ## ## ## ## ## ## ## #
####      Simplified L2SWBM model      ####
###            PLOTTER SCRIPT           ###
###      Dani Cohn - May 17, 2021       ###
###  Based on the Lake Chilwa WB Model  ###
## ## ## ## ## ## ## ## ## ## ## ## ## ## #

## This is a long script (700+ lines)!    
## If using RStudio, you can use the section finder at the bottom-left of the code window to navigate.

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## #
#### 1. set up working directory of the data ####
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## #

# Plotting.on.iteration.time=20000
# Plotting.on.iteration.time=iters
print("opening plotter")
library(rjags)
setwd("~/Desktop/CGLGP_CIA/L2SWBM/Sarahs_L2S")
load(paste(outputname,startAnalysisYear,endAnalysisYear,'_',iters/1000,'k_roll1.RData',sep=''))
fileNameBase=paste(outputname,iters/1000,'k',sep='')

eriOut = jSumStats_MSD[445:660,1]
eriOut = t(matrix(eriOut, nrow = 12, ncol = 18))
eriBOM = erieBOM[1201:1416,3]
eriBOM = t(matrix(eriBOM, nrow = 12, ncol = 18))

stClOut = jSumStats_MSD[1969:2184,1]
stClOut = t(matrix(stClOut, nrow = 12, ncol = 18))
stClBOM = stClairBOM[1201:1416,3]
stClBOM = t(matrix(stClBOM, nrow = 12, ncol = 18))

miHuOut = jSumStats_MSD[3061:3276,1]
miHuOut = t(matrix(miHuOut, nrow = 12, ncol = 18))
miHuBOM = miHuronBOM[1201:1416,3]
miHuBOM = t(matrix(miHuBOM, nrow = 12, ncol = 18))

supeOut = jSumStats_MSD[4585:4800,1]
supeOut = t(matrix(supeOut, nrow = 12, ncol = 18))
supeBOM = superiorBOM[1201:1416,3]
supeBOM = t(matrix(supeBOM, nrow = 12, ncol = 18))

diffSM = supeBOM - miHuBOM
diffMC = miHuBOM - stClBOM
diffCE = stClBOM - eriBOM

## Generate list of vectors to hold regression data
# 1 - Lake Superior Intercept      4 - Lake Michigan & Huron Intercept      7 - Lake St Claire Intercept    10 - Lake Erie Intercept
# 2 - Lake Superior Coefficient    5 - Lake Michigan & Huron Coefficient    8 - Lake St Claire Coefficient  11 = Lake Erie Coefficient
# 3 - Lake Superior Diff           6 - Lake Michigan & Huron Diff           9 - Lake St Claire Diff
regressData <- vector(mode = "list", length = 11)
for (i in 1:11) regressData[[i]]<- vector(length = 12)

predSup = matrix(nrow = 18, ncol = 12)
predMhn = matrix(nrow = 18, ncol = 12)
predStc = matrix(nrow = 18, ncol = 12)
predEri = matrix(nrow = 18, ncol = 12)

for (i in 1:12) {
   model_superior <- lm(log(supeOut[,i]) ~ log(supeBOM[,i]) + log(diffSM[,i]))
   regressData[[1]][i] <- model_superior$coefficients[1]
   regressData[[2]][i] <- model_superior$coefficients[2]
   regressData[[3]][i] <- model_superior$coefficients[3]
   predSup[,i] = exp(model_superior$fitted.values)
   
   model_mhn <- lm(log(miHuOut[,i]) ~ log(miHuBOM[,i]) + log(diffMC[,i]))
   regressData[[4]][i] <- model_mhn$coefficients[1]
   regressData[[5]][i] <- model_mhn$coefficients[2]
   regressData[[6]][i] <- model_mhn$coefficients[3]
   predMhn[,i] = exp(model_mhn$fitted.values)
   
   model_stc <- lm(log(stClOut[,i]) ~ log(stClBOM[,i]) + log(diffCE[,i]))
   regressData[[7]][i] <- model_stc$coefficients[1]
   regressData[[8]][i] <- model_stc$coefficients[2]
   regressData[[9]][i] <- model_stc$coefficients[3]
   predStc[,i] = exp(model_stc$fitted.values)
   
   model_eri <- lm(log(eriOut[,i]) ~ log(eriBOM[,i]))
   regressData[[10]][i] <- model_eri$coefficients[1]
   regressData[[11]][i] <- model_eri$coefficients[2]
   predEri[,i] = exp(model_eri$fitted.values)
}

# for (i in 1:12) {
#    predSup[,i] = 1000*(predSup[,i]*secondsInADay*daysInMonthsWithLeap[i])/supArea
#    predMhn[,i] = 1000*(predMhn[,i]*secondsInADay*daysInMonthsWithLeap[i])/mhgArea
#    predStc[,i] = 1000*(predStc[,i]*secondsInADay*daysInMonthsWithLeap[i])/lkstclArea
#    predEri[,i] = 1000*(predEri[,i]*secondsInADay*daysInMonthsWithLeap[i])/eriArea
# }

predSup = matrix(t(predSup), nrow = 216, ncol = 1)

print(predSup)

predMhn = matrix(t(predMhn), nrow = 216, ncol = 1)
predStc = matrix(t(predStc), nrow = 216, ncol = 1)
predEri = matrix(t(predEri), nrow = 216, ncol = 1)


## ## ## ## ## ## ## ## ## ## ## ## ##
#### 2. send summary table to csv ####
## ## ## ## ## ## ## ## ## ## ## ## ##

# write.table(jSum,paste(fileNameBase,'_stats.csv',sep=''),sep=',',quote=FALSE )
## Moved above line to end of main wrapper script. ##
## Then if there's a plotter issue, the stats csv is still available. ##


## ## ## ## ## ## ## ## ## ## ## #
#### 3. yMHDStore & yMHRStore ####
## ## ## ## ## ## ## ## ## ## ## #

## This section deleted by Dani b/c DStore and PP (1Y, 2Y,...) not used yet ##

## ## ## ## ## ## ## ## ## ## ##
#### 3.5 c, d, f value Plots ####
## ## ## ## ## ## ## ## ## ## ##

#### Makes PDFs of c, d, f values
##### SFD c #####
# pdf(paste("SFDc_ValuePlot_SupAddIn.pdf"),width=11,height=7.5)
# # plot(jSample[1][,2162],type="l") ## c ##
# hist(jSample[1][,paste("SFDc")][[1]], freq = FALSE, main=NULL,xlab=NULL, xlim = c(100,400))
# input.c = seq(100,400,.1)
# lines(x = input.c,
#       y = dnorm(input.c, mean = 250, sd = 50),
#       col = "red")
# abline(v = 201.03, col = "blue", lwd = 2)
# title(main="SFD equation parameter 'c' from L2SWBM \nfor Fort Gratiot - Dry Dock reach",cex.main=1.75,xlab="Value of 'c'")
# legend("topright",lty=c(1,1),bty="n",legend = c("Prior","Value used by USACE","L2SWBM output"),ncol=1,cex=0.9,col = c("red","blue","grey80"))
# # legend("topright",fill="grey80",bty="n",legend = c("L2SWBM output (posterior)"),ncol=1,yjust=1,cex=0.9,border="grey80") ## ncol=4 here to shift this legend entry to the left ##
# dev.off()
# 
# ##### SFD d #####
# pdf(paste("SFDd_ValuePlot_SupAddIn.pdf"),width=11,height=7.5)
# # plot(jSample[1][,2163],type="l") ## d ##
# hist(jSample[1][,paste("SFDd")][[1]], freq = FALSE, main=NULL,xlab=NULL, xlim = c(-1,4))
# input.d = seq(-3,6,.01)
# lines(x = input.d,
#       y = dnorm(input.d, mean = 1.6, sd = 2),
#       col = "red")
# abline(v = 1.732, col = "blue", lwd = 2)
# title(main="SFD equation parameter 'd' from L2SWBM \nfor Fort Gratiot - Dry Dock reach",cex.main=1.75,xlab="Value of 'd'")
# legend("topright",lty=c(1,1),bty="n",legend = c("Prior","Value used by USACE","L2SWBM output"),ncol=1,cex=0.9,col = c("red","blue","grey80"))
# dev.off()
# 
# ##### SFD f #####
# pdf(paste("SFDf_ValuePlot_SupAddIn.pdf"),width=11,height=7.5)
# # plot(jSample[1][,2164],type="l") ## f ##
# hist(jSample[1][,paste("SFDf")][[1]], freq = FALSE, main=NULL,xlab=NULL, xlim = c(-2,3))
# input.f = seq(-3,3,.1)
# lines(x = input.f,
#       y = dnorm(input.f, mean = 0.5, sd = 1),
#       col = "red")
# abline(v = 0.428, col = "blue", lwd = 2)
# title(main="SFD equation parameter 'f' from L2SWBM \nfor Fort Gratiot - Dry Dock reach",cex.main=1.75,xlab="Value of 'f'")
# legend("topright",lty=c(1,1),bty="n",legend = c("Prior","Value used by USACE","L2SWBM output"),ncol=1,cex=0.9,col = c("red","blue","grey80"))
# dev.off()


#### Makes PDF of "OutflowSFDPrecisionMH" values
##### OutflowSFDPrecisionMH #####
# pdf(paste("OutflowSFDPrecision_ValuePlot_SupAddIn.pdf"),width=11,height=7.5)
# hist(jSample[1][,paste("OutflowSFDPrecisionMH")][[1]], freq = FALSE, main=NULL,xlab=NULL)
# # input.c = seq(100,400,.1)
# # lines(x = input.c,
# #       y = dnorm(input.c, mean = 250, sd = 50),
# #       col = "red")
# # abline(v = 201.03, col = "blue", lwd = 2)
# title(main="SFD-Outflow Precision from L2SWBM \nfor Fort Gratiot - Dry Dock reach",cex.main=1.75)
# legend("topright",lty=c(1,1),bty="n",legend = c("L2SWBM output"),ncol=1,cex=0.9,col = c("grey80"))
# # legend("topright",fill="grey80",bty="n",legend = c("L2SWBM output (posterior)"),ncol=1,yjust=1,cex=0.9,border="grey80") ## ncol=4 here to shift this legend entry to the left ##
# dev.off()


## ## ## ## ## ## ## ## ## ##
#### 4. Time Series Plots ####
## ## ## ## ## ## ## ## ## ##

yearRange  =startAnalysisYear:endAnalysisYear
nMonths    =12*length(yearRange)
# compLimits=c(-50,400)
# diversionLimits=c(-25,100)
# compLimitsLabs=c(0,100,200,300,400,500)
# compLimitsTicks=seq(-50,450,50)
# diversionLimitsLabs=c(0,25,50,75,100)
# diversionLimitsTicks=seq(-25,100,25)
# storeLimits=c(-650,650)
# clairFlowLimits=c(7000,20000)


## ## ## ## ## ## ## ## ## #
#### 5. Time Series M-H ####
## ## ## ## ## ## ## ## ## #

pdf(paste(fileNameBase, '_TimeSeries_MHnoSFD.pdf', sep=''), width = 10, height = 7.5);
par(mfrow=c(6,1)) # HD adjusted to make all panels fall inside one page 
## COMMENT ON 17 FEB: PLEASE NOTE THAT WE PLOT OUT THE TOTAL NUMBER OF SIX VARIABLES
##                    THUS YOU SHOULD SET mfrow = c(6,1)
##                    IF DIVERSION IS ALSO PLOTTED, SET mfrow=c(7,1) AND SO ON
par(mar = c(0,0,0,0))
par(oma = c(4,4,4,4))

## ADDED ON FEB 16, 2021 ##
## to separate time series plots into several 5-year blocks ##
nYears = length(yearRange)
intv   = 5
nPlots = ceiling(nYears/intv)

for (ii in 1:nPlots) {
   xlim1 = (ii-1)*intv*12+1
   xlim2 = ii*intv*12
   ##### Inflow #####
   yILimit = range(c(yMHInfl1,yMHInfl2),
#                     yMHInfl3),
                   na.rm=T)
   yIDiff = diff(yILimit)
   plot(c(0), c(0), type = "n", col = 8, lwd = 3, axes = FALSE, xlim = c(xlim1,xlim2), ylim = c(yILimit[1]-0.05*yIDiff,yILimit[2]+0.05*yIDiff));#c(yILimit[1]*0.9,yILimit[2]*1.1));
   box()
   axis(2, cex.axis = 0.9);
   # axis(2, at=seq(0,700,100), labels=FALSE, cex.axis = 0.9);
   axis(4, labels=FALSE, cex.axis = 0.9);
   mtext(paste("Inflow (cms)"), side = 2, line = 2.5, cex=0.8);
   axis(1, at=seq(0,nMonths,12), labels=FALSE);
   #
   xpol = c(1:nMonths,nMonths:1)
   ypol = c(jSum[paste('SupOutflow[',1:nMonths,']', sep=''),3],
            jSum[paste('SupOutflow[',nMonths:1,']', sep=''),7])
   polygon(x = xpol,y = ypol,col="grey80",border="grey80")
   #
   lines(yMHInfl1, col = "red",  lwd = 1)
   lines(yMHInfl2, col = "blue",  lwd = 1,lty=2)
   #lines(yMHInfl3, col = "green4",  lwd = 1,lty=3)
   legend("topleft",lty=c(1,2,3),bty="n",legend = c("IGS","Flow.Accounting")
 #         ,"JJ21")
      ,ncol=3,cex=0.9,col = c("red","blue")) #,"green4"
   legend("top",fill="grey80",bty="n",legend = c("L2SWBM.95CI"),ncol=4,cex=0.9,border="grey80") ## ncol=4 here to shift this legend entry to the left ##
   #
   title(main='Lake Michigan-Huron Time Series', outer=TRUE,cex.main=1.75);
   
   #### Precip ####
   yPLimit = range(c(yMHPrecip1,yMHPrecip2,yMHPrecip3,yMHPrecip4,
                     yMHPrecip5,yMHPrecip6,yMHPrecip7,yMHPrecip8,
                     yMHPrecip9),na.rm=T)
   plot(c(0), c(0), type = "n", col = 8, lwd = 3, axes = FALSE, xlim = c(xlim1,xlim2), ylim = yPLimit*1.1); 
   box()
   # axis(2, cex.axis = 0.9); 
   axis(4, labels=T, cex.axis = 0.9); 
   axis(2, labels=FALSE, cex.axis = 0.9); 
   mtext(paste("Precip (mm)"), side = 2, line = 2.5, cex=0.8); 
   axis(1, at=seq(0,nMonths,12), labels=FALSE);
   #
   ypol = c(jSum[paste('MHPrecip[',1:nMonths,']', sep=''),3],
            jSum[paste('MHPrecip[',nMonths:1,']', sep=''),7])
   polygon(x = xpol,y = ypol,col="grey80",border="grey80")
   #
   lines(yMHPrecip1, col = "red",  lwd = 1)
   lines(yMHPrecip2, col = "blue",  lwd = 1,lty=2)
   lines(yMHPrecip3, col = "green4", lwd = 1,lty=4)
   lines(yMHPrecip4, col = "navy",  lwd = 1,lty=5)
   lines(yMHPrecip5, col = "magenta",  lwd = 1)
   lines(yMHPrecip6, col = "steelblue",  lwd = 1,lty=2)
   lines(yMHPrecip7, col = "deeppink",  lwd = 1,lty=4)
   lines(yMHPrecip8, col = "purple",  lwd = 1,lty=5)
   lines(yMHPrecip9, col = "dodgerblue",  lwd = 1)
   legend("topleft",lty=c(1,2,4,5,1,2,4,5,1),bty="n",
          legend = c("NOAA GLERL GLM HMD","GLERL AHPS PROVISIONAL","USACE AHPS","ECCC WCPS",
                     "ECCC CaPA","NWS MPE","HISTORICAL COORDINATED","USACE Thiessen",
                     "MERGED MPE CaPA"),ncol=4,cex=0.9,col = c("red","blue","green4","navy",
                                                               "magenta","steelblue","deeppink",
                                                               "purple","dodgerblue"))
   legend("top",fill="grey80",bty="n",legend = c("L2SWBM.95CI"),ncol=4,cex=0.9,border="grey80") ## ncol=4 here to shift this legend entry to the left ##
   
   #### Evap ####
   yELimit = range(c(yMHEvap1,
                     #yMHEvap2,
                     yMHEvap3,yMHEvap4,yMHEvap5),na.rm=T)
   plot(c(0), c(0), type = "n", col = 8, lwd = 3, axes = FALSE, xlim = c(xlim1,xlim2), ylim = yELimit*1.1); 
   box()
   # axis(2, cex.axis = 0.9); 
   axis(2, labels=T, cex.axis = 0.9); 
   axis(4, labels=FALSE, cex.axis = 0.9); 
   mtext(paste("Evap (mm)"), side = 2, line = 2.5, cex=0.8); 
   axis(1, at=seq(0,nMonths,12), labels=FALSE);
   #
   ypol = c(jSum[paste('MHEvap[',1:nMonths,']', sep=''),3],
            jSum[paste('MHEvap[',nMonths:1,']', sep=''),7])
   polygon(x = xpol,y = ypol,col="grey80",border="grey80")
   #
   lines(yMHEvap1, col = "red",  lwd = 1)
   #lines(yMHEvap2, col = "blue",  lwd = 1,lty=2)
   lines(yMHEvap3, col = "green4", lwd = 1,lty=4)
   lines(yMHEvap4, col = "navy",  lwd = 1,lty=5)
   lines(yMHEvap5, col = "magenta",  lwd = 1)
   legend("topleft",lty=c(1,4,5,1),
          #lty=c(1,2,4,5,1),
          bty="n",legend = c("NOAA GLERL GLM HMD",
                                                        #"GLM HMD PROVISIONAL",
                                                      "USACE AHPS","ECCC WCPS","GLERL FVCOM"),
          ncol=3,cex=0.9,col = c("red",
                                 #"blue",
                                 "green4","navy","magenta"))
   legend("top",fill="grey80",bty="n",legend = c("L2SWBM.95CI"),ncol=4,cex=0.9,border="grey80") ## ncol=4 here to shift this legend entry to the left ##
   
   #### Runoff ####
   yRLimit = range(c(yMHRunoff1,
                     #yMHRunoff2,
                     yMHRunoff3,yMHRunoff4,yMHRunoff5),na.rm=T)
   plot(c(0), c(0), type = "n", col = 8, lwd = 3, axes = FALSE, xlim = c(xlim1,xlim2), ylim = yRLimit*1.1); 
   box()
   axis(4, cex.axis = 0.9); 
   # axis(2, at=seq(0,700,100), labels=FALSE, cex.axis = 0.9); 
   axis(2, labels=FALSE, cex.axis = 0.9); 
   mtext(paste("Runoff (mm)"), side = 2, line = 2.5, cex=0.8); 
   axis(1, at=seq(0,nMonths,12), labels=FALSE);
   #
   ypol = c(jSum[paste('MHRunoff[',1:nMonths,']', sep=''),3],
            jSum[paste('MHRunoff[',nMonths:1,']', sep=''),7])
   polygon(x = xpol,y = ypol,col="grey80",border="grey80")
   #
   lines(yMHRunoff1, col = "red",  lwd = 1)
   #lines(yMHRunoff2, col = "blue",  lwd = 1,lty=2)
   lines(yMHRunoff3, col = "green4", lwd = 1,lty=4)
   lines(yMHRunoff4, col = "navy",  lwd = 1,lty=5)
   lines(yMHRunoff5, col = "magenta",  lwd = 1)
   legend("topleft",lty=c(1,4,5,1),
          #lty=c(1,2,4,5,1),
          bty="n",legend = c("NOAA GLERL GLM HMD",
                                                        #"GLM HMD PROVISIONAL",
                                                        "USACE AHPS","ECCC WCPS","ECCC WATFLOOD"),
          ncol=4,cex=0.9,col = c("red",
                                 #"blue",
                                 "green4","navy","magenta"))
   legend("top",fill="grey80",bty="n",legend = c("L2SWBM.95CI"),ncol=4,cex=0.9,border="grey80") ## ncol=4 here to shift this legend entry to the left ##
   
   #### Outflow ####
   # yMHSFDOutflow = jSum[paste('SFDc'),1]*((((yMHOutFtGr[217:444]+yMHOutDryD[217:444])/2)-167)^jSum[paste('SFDd'),1])*((yMHOutFtGr[217:444] - yMHOutDryD[217:444])^jSum[paste('SFDf'),1])
   yQLimit = range(c(yMHOut1,yMHOut2
                     #yMHOut3
                     ),na.rm=T)
   yQDiff = diff(yQLimit)
   # yQLimit = range(c(yMHSFDOut1),na.rm=T)
   # yQLimit = c(0,8000)
   plot(c(0), c(0), type = "n", col = 8, lwd = 3, axes = FALSE, xlim = c(xlim1,xlim2), ylim = c(yQLimit[1]-0.05*yQDiff,yQLimit[2]+0.05*yQDiff));#c(yQLimit[1]*0.9,yQLimit[2]*1.1));
   box()
   axis(2, cex.axis = 0.9);
   # axis(2, at=seq(0,700,100), labels=FALSE, cex.axis = 0.9);
   axis(4, labels=FALSE, cex.axis = 0.9);
   mtext(paste("Outflow (cms)"), side = 2, line = 2.5, cex=0.8);
   axis(1, at=seq(0,nMonths,12), labels=FALSE);
   #
   ypol = c(jSum[paste('MHOutflow[',1:nMonths,']', sep=''),3], ## 2.5% value ##
            jSum[paste('MHOutflow[',nMonths:1,']', sep=''),7]) ## 97.5% value ##
   polygon(x = xpol,y = ypol,col="grey80",border="grey80")
   #
   # ypol2 = c(jSum[paste('OutflowSFDMeanMH[',1:nMonths,']', sep=''),3], ## 2.5% value ##
   #          jSum[paste('OutflowSFDMeanMH[',nMonths:1,']', sep=''),7]) ## 97.5% value ##
   # polygon(x = xpol,y = ypol2,border="mediumpurple2")
   #
   # lines(yMHOutflow, col = "red",  lwd = 1)
   lines(yMHOut1, col = "red",  lwd = 1,lty=1)
   lines(yMHOut2, col = "blue",  lwd = 1,lty=2)
 #  lines(yMHOut3, col = "green4",  lwd = 1,lty=4)
   lines(predMhn, col = "navy", lwd = 1, lty=5)
   # lines(yMHSFDOutflow,col = "seagreen4",lwd = 1, lty=2) ## Paired with legend code on next line! If this is uncommented, uncomment legend too! ##
   # legend("topleft",lty=c(1,2,2,1),bty="n",legend = c("IGS","SFD.ADVM","SFD from model"),ncol=3,cex=0.9,col = c("red","blue","seagreen4"))
   legend("topleft",lty=c(1,2,5),
          #lty=c(1,2,4,5),
          bty="n",legend = c("IGS","SFD.ADVM (*EXCLUDED FROM MODEL*)",
                                                      #"JJ21",
                                                      "LM Fit"),ncol=4,cex=0.9,col = c("red","blue",
                                                                                       #"green4",
                                                                                       "navy")) ## Comment out if above legend code is active ##
   legend("top",fill="grey80",bty="n",legend = c("L2SWBM.95CI.MHOut"),ncol=2,cex=0.9,border="grey80") ## ncol=4 here to shift this legend entry to the left ##
   # legend("topright",fill="mediumpurple2",bty="n",legend = c("L2SWBM.95CI.OutSFD"),ncol=3,cex=0.9,border="mediumpurple2") ## ncol=4 here to shift this legend entry to the left ##
   
   #### Dstore ####
   yDSLimit = range(yMHDStore,na.rm=T)
   plot(yMHDStore, type = "n", col = "black", lwd = 2, axes = FALSE, ylim = c(yDSLimit[1],yDSLimit[2]*1.1), xlim = c(xlim1,xlim2)); 
   box()
   abline(h=0, col = 8)
   axis(4,cex.axis = 0.9);
   axis(2, labels=FALSE, cex.axis = 0.9); 
   axis(1, at=seq(0,nMonths,12), labels=FALSE);
   axis(1, at=seq(0,nMonths-1,12)+6, labels=yearRange, tick=FALSE);
   mtext('DStore (mm)', side = 2, line = 2.5, cex = 0.8)
   #
   ypol = c(jSum[paste('MHDStore[',1:nMonths,']', sep=''),3], ## Reads in the 2.5% column. ##
            jSum[paste('MHDStore[',nMonths:1,']', sep=''),7]) ## Reads in the 97.5% column. Between these is 95%!! ##
   polygon(x = xpol,y = ypol,col="grey80",border="grey80")
   #
   lines(yMHDStore, col = "red",  lwd = 1,lty=1)
   legend("topleft",lty=1,bty="n",legend = c("Observed DStore"),ncol=2,cex=0.9,col = c("red","blue"))
   legend("top",fill="grey80",bty="n",legend = c("Inferred DStore"),ncol=7,cex=0.9,border="grey80") ## ncol=7 here to shift this legend entry to the left ##
   
}

dev.off();




## ## ## ## ## ## ## ## ## ## ##
#### 6. Time Series Superior ####
## ## ## ## ## ## ## ## ## ## ##

pdf(paste(fileNameBase, '_TimeSeriesSuperior.pdf', sep=''), width = 10, height = 7.5);
par(mfrow=c(6,1)) # HD adjusted to make all panels fall inside one page 
## COMMENT ON 17 FEB: PLEASE NOTE THAT WE PLOT OUT THE TOTAL NUMBER OF SIX VARIABLES
##                    THUS YOU SHOULD SET mfrow = c(6,1)
##                    IF DIVERSION IS ALSO PLOTTED, SET mfrow=c(7,1) AND SO ON
par(mar = c(0,0,0,0))
par(oma = c(4,4,4,4))

## ADDED ON FEB 16, 2021 ## ADJUSTED ON APR 6, 2021 ##
## to separate time series plots into several 5-year blocks ##
nYears = length(yearRange)
intv   = 5
nPlots = ceiling(nYears/intv)

for (ii in 1:nPlots) {
   xlim1 = (ii-1)*intv*12+1
   xlim2 = ii*intv*12
   
   ###### Inflow ######
   plot(c(0), c(0), type = "n", col = 8, lwd = 3, axes = FALSE);
   box()
   axis(2, labels=FALSE, cex.axis = 0.9);
   # axis(2, at=seq(0,700,100), labels=FALSE, cex.axis = 0.9);
   axis(4, labels=FALSE, cex.axis = 0.9);
   mtext(paste("Inflow (cms)"), side = 2, line = 2.5, cex=0.8);
   axis(1, at=seq(0,nMonths,12), labels=FALSE);
   #
   legend("topleft",lty=1,lwd=0,bty="n",legend = c("Note: Lake Superior does not have an inflow source."),ncol=1,cex=0.9)
   #
   title(main='Lake Superior Time Series', outer=TRUE,cex.main=1.75);
   
   #### Precip ####
   yPLimit = range(c(ySupPrecip1,
                     #ySupPrecip2,
                     ySupPrecip3,
                     ySupPrecip4,ySupPrecip5,ySupPrecip6,
                     #ySupPrecip7,
                     ySupPrecip8,ySupPrecip9),na.rm=T)
   plot(c(0), c(0), type = "n", col = 8, lwd = 3, axes = FALSE, xlim = c(xlim1,xlim2), ylim = yPLimit*1.1); 
   box()
   # axis(2, cex.axis = 0.9); 
   axis(4, labels=T, cex.axis = 0.9); 
   axis(2, labels=FALSE, cex.axis = 0.9); 
   mtext(paste("Precip (mm)"), side = 2, line = 2.5, cex=0.8); 
   axis(1, at=seq(0,nMonths,12), labels=FALSE);
   #
   ypol = c(jSum[paste('SupPrecip[',1:nMonths,']', sep=''),3],
            jSum[paste('SupPrecip[',nMonths:1,']', sep=''),7])
   polygon(x = xpol,y = ypol,col="grey80",border="grey80")
   #
   lines(ySupPrecip1, col = "red",  lwd = 1)
   #lines(ySupPrecip2, col = "blue",  lwd = 1,lty=2)
   lines(ySupPrecip3, col = "green4", lwd = 1,lty=4)
   lines(ySupPrecip4, col = "navy",  lwd = 1,lty=5)
   lines(ySupPrecip5, col = "magenta",  lwd = 1)
   lines(ySupPrecip6, col = "steelblue",  lwd = 1,lty=2)
   #lines(ySupPrecip7, col = "deeppink",  lwd = 1,lty=4)
   lines(ySupPrecip8, col = "purple",  lwd = 1,lty=5)
   lines(ySupPrecip9, col = "dodgerblue",  lwd = 1)
   legend("topleft",lty=c(1,4,5,1,2,5,1),
          #lty=c(1,2,4,5,1,2,4,5,1),
          bty="n",
          legend = c("NOAA GLERL GLM HMD",
                     #"GLERL AHPS PROVISIONAL",
                     "USACE AHPS","ECCC WCPS",
                     "ECCC CaPA","NWS MPE",
                     #"HISTORICAL COORDINATED",
                     "USACE Thiessen",
                     "MERGED MPE CaPA"),ncol=4,cex=0.9,col = c("red",
                                                               #"blue",
                                                               "green4","navy","magenta","steelblue",
                                                               #"deeppink",
                                                               "purple","dodgerblue"))
   legend("top",fill="grey80",bty="n",legend = c("L2SWBM.95CI"),ncol=4,cex=0.9,border="grey80") ## ncol=4 here to shift this legend entry to the left ##
   
   #### Evap ####
   yELimit = range(c(ySupEvap1,#ySupEvap2,
                     ySupEvap3,ySupEvap4,ySupEvap5),na.rm=T)
   plot(c(0), c(0), type = "n", col = 8, lwd = 3, axes = FALSE, xlim = c(xlim1,xlim2), ylim = yELimit*1.1); 
   box()
   # axis(2, cex.axis = 0.9); 
   axis(2, labels=T, cex.axis = 0.9); 
   axis(4, labels=FALSE, cex.axis = 0.9); 
   mtext(paste("Evap (mm)"), side = 2, line = 2.5, cex=0.8); 
   axis(1, at=seq(0,nMonths,12), labels=FALSE);
   #
   ypol = c(jSum[paste('SupEvap[',1:nMonths,']', sep=''),3],
            jSum[paste('SupEvap[',nMonths:1,']', sep=''),7])
   polygon(x = xpol,y = ypol,col="grey80",border="grey80")
   #
   lines(ySupEvap1, col = "red",  lwd = 1)
   #lines(ySupEvap2, col = "blue",  lwd = 1,lty=2)
   lines(ySupEvap3, col = "green4", lwd = 1,lty=4)
   lines(ySupEvap4, col = "navy",  lwd = 1,lty=5)
   lines(ySupEvap5, col = "magenta",  lwd = 1)
   legend("topleft",lty=c(1,4,5,1), #lty=c(1,2,4,5,1),
          bty="n",legend = c("NOAA GLERL GLM HMD",
                                                        #"GLM HMD PROVISIONAL",
                                                        "USACE AHPS","ECCC WCPS","GLERL FVCOM"),
          ncol=3,cex=0.9,col = c("red",
                                 #"blue",
                                 "green4","navy","magenta"))
   legend("top",fill="grey80",bty="n",legend = c("L2SWBM.95CI"),ncol=4,cex=0.9,border="grey80") ## ncol=4 here to shift this legend entry to the left ##
   
   #### Runoff ####
   yRLimit = range(c(ySupRunoff1,
                     #ySupRunoff2,
                     ySupRunoff3,ySupRunoff4,ySupRunoff5),na.rm=T)
   plot(c(0), c(0), type = "n", col = 8, lwd = 3, axes = FALSE, xlim = c(xlim1,xlim2), ylim = yRLimit*1.1); 
   box()
   axis(4, cex.axis = 0.9); 
   # axis(2, at=seq(0,700,100), labels=FALSE, cex.axis = 0.9); 
   axis(2, labels=FALSE, cex.axis = 0.9); 
   mtext(paste("Runoff (mm)"), side = 2, line = 2.5, cex=0.8); 
   axis(1, at=seq(0,nMonths,12), labels=FALSE);
   #
   ypol = c(jSum[paste('SupRunoff[',1:nMonths,']', sep=''),3],
            jSum[paste('SupRunoff[',nMonths:1,']', sep=''),7])
   polygon(x = xpol,y = ypol,col="grey80",border="grey80")
   #
   lines(ySupRunoff1, col = "red",  lwd = 1)
   #lines(ySupRunoff2, col = "blue",  lwd = 1,lty=2)
   lines(ySupRunoff3, col = "green4", lwd = 1,lty=4)
   lines(ySupRunoff4, col = "navy",  lwd = 1,lty=5)
   lines(ySupRunoff5, col = "magenta",  lwd = 1)
   legend("topleft",lty=c(1,4,5,1), #lty=c(1,2,4,5,1),
          bty="n",legend = c("NOAA GLERL GLM HMD",
                        #"GLM HMD PROVISIONAL",
                        "USACE AHPS","ECCC WCPS","ECCC WATFLOOD"),
          ncol=4,cex=0.9,col = c("red",
                                 #"blue",
                                 "green4","navy","magenta"))
   legend("top",fill="grey80",bty="n",legend = c("L2SWBM.95CI"),ncol=4,cex=0.9,border="grey80") ## ncol=4 here to shift this legend entry to the left ##
   
   #### Outflow ####
   yQLimit = range(c(ySupOut1,ySupOut2
                     #,ySupOut3
                     ),na.rm=T)
   yQDiff = diff(yQLimit)
   # yQLimit = c(0,8000)
   plot(c(0), c(0), type = "n", col = 8, lwd = 3, axes = FALSE, xlim = c(xlim1,xlim2), ylim = c(yQLimit[1]-0.05*yQDiff,yQLimit[2]+0.05*yQDiff));#c(yQLimit[1]*0.9,yQLimit[2]*1.1));
   box()
   axis(2, cex.axis = 0.9);
   # axis(2, at=seq(0,700,100), labels=FALSE, cex.axis = 0.9);
   axis(4, labels=FALSE, cex.axis = 0.9);
   mtext(paste("Outflow (cms)"), side = 2, line = 2.5, cex=0.8);
   axis(1, at=seq(0,nMonths,12), labels=FALSE);
   #
   ypol = c(jSum[paste('SupOutflow[',1:nMonths,']', sep=''),3], ## 2.5% value ##
            jSum[paste('SupOutflow[',nMonths:1,']', sep=''),7]) ## 97.5% value ##
   polygon(x = xpol,y = ypol,col="grey80",border="grey80")
   #
   lines(ySupOut1, col = "red",  lwd = 1,lty=1)
   lines(ySupOut2, col = "blue",  lwd = 1,lty=2)
   #lines(ySupOut3, col = "green4",  lwd = 1,lty=4)
   lines(predSup, col = "navy", lwd = 1, lty=5)
   legend("topleft",lty=c(1,2,5), #lty=c(1,2,4,5),
          bty="n",legend = c("IGS","Flow Accounting",
                                                      #"JJ21",
                                                      "LM Fit"),ncol=4,cex=0.9,col = c("red","blue",#"green4",
                                                                                       "navy"))
   legend("top",fill="grey80",bty="n",legend = c("L2SWBM.95CI"),ncol=4,cex=0.9,border="grey80") ## ncol=4 here to shift this legend entry to the left ##
   
   #### Dstore ####
   yDSLimit = range(ySupDStore,na.rm=T)
   plot(ySupDStore, type = "n", col = "black", lwd = 2, axes = FALSE, ylim = c(yDSLimit[1],yDSLimit[2]*1.1), xlim = c(xlim1,xlim2)); 
   box()
   abline(h=0, col = 8)
   axis(4,cex.axis = 0.9);
   axis(2, labels=FALSE, cex.axis = 0.9); 
   axis(1, at=seq(0,nMonths,12), labels=FALSE);
   axis(1, at=seq(0,nMonths-1,12)+6, labels=yearRange, tick=FALSE);
   mtext('DStore (mm)', side = 2, line = 2.5, cex = 0.8)
   #
   ypol = c(jSum[paste('SupDStore[',1:nMonths,']', sep=''),3], ## Reads in the 2.5% column. ##
            jSum[paste('SupDStore[',nMonths:1,']', sep=''),7]) ## Reads in the 97.5% column. Between these is 95%!! ##
   polygon(x = xpol,y = ypol,col="grey80",border="grey80")
   #
   lines(ySupDStore, col = "red",  lwd = 1,lty=1)
   legend("topleft",lty=1,bty="n",legend = c("Observed DStore"),ncol=2,cex=0.9,col = c("red","blue"))
   legend("top",fill="grey80",bty="n",legend = c("Inferred DStore"),ncol=7,cex=0.9,border="grey80") ## ncol=7 here to shift this legend entry to the left ##
   
}

dev.off();





## ## ## ## ## ## ## ## ## ## ## ## ##
#### 7. Time Series Lake St Clair ####
## ## ## ## ## ## ## ## ## ## ## ## ##

pdf(paste(fileNameBase, '_TimeSeries_LakeStClair.pdf', sep=''), width = 10, height = 7.5);
par(mfrow=c(6,1)) # HD adjusted to make all panels fall inside one page 
## COMMENT ON 17 FEB: PLEASE NOTE THAT WE PLOT OUT THE TOTAL NUMBER OF SIX VARIABLES
##                    THUS YOU SHOULD SET mfrow = c(6,1)
##                    IF DIVERSION IS ALSO PLOTTED, SET mfrow=c(7,1) AND SO ON
par(mar = c(0,0,0,0))
par(oma = c(4,4,4,4))

## ADDED ON FEB 16, 2021 ##
## to separate time series plots into several 5-year blocks ##
nYears = length(yearRange)
intv   = 5
nPlots = ceiling(nYears/intv)

for (ii in 1:nPlots) {
   xlim1 = (ii-1)*intv*12+1
   xlim2 = ii*intv*12
   ##### Inflow #####
   yILimit = range(c(yLkStClInfl1,yLkStClInfl2
                     #,yLkStClInfl3
                     ),na.rm=T)
   yIDiff = diff(yILimit)
   plot(c(0), c(0), type = "n", col = 8, lwd = 3, axes = FALSE, xlim = c(xlim1,xlim2), ylim = c(yILimit[1]-0.05*yIDiff,yILimit[2]+0.05*yIDiff));#c(yILimit[1]*0.9,yILimit[2]*1.1));
   box()
   axis(2, cex.axis = 0.9);
   # axis(2, at=seq(0,700,100), labels=FALSE, cex.axis = 0.9);
   axis(4, labels=FALSE, cex.axis = 0.9);
   mtext(paste("Inflow (cms)"), side = 2, line = 2.5, cex=0.8);
   axis(1, at=seq(0,nMonths,12), labels=FALSE);
   #
   xpol = c(1:nMonths,nMonths:1)
   ypol = c(jSum[paste('MHOutflow[',1:nMonths,']', sep=''),3],
            jSum[paste('MHOutflow[',nMonths:1,']', sep=''),7])
   polygon(x = xpol,y = ypol,col="grey80",border="grey80")
   #
   lines(yLkStClInfl1, col = "red",  lwd = 1)
   lines(yLkStClInfl2, col = "blue",  lwd = 1,lty=2)
   #lines(yLkStClInfl3, col = "blue",  lwd = 1,lty=3)
   legend("topleft",lty=c(1,2), #lty=c(1,2,3),
          bty="n",legend = c("IGS","Flow.Accounting (*EXCLUDED FROM MODEL*)"
                                                    #,"JJ21"
                                                    ),ncol=2,cex=0.9,col = c("red","blue")) #,"green4"
   legend("top",fill="grey80",bty="n",legend = c("L2SWBM.95CI"),ncol=4,cex=0.9,border="grey80") ## ncol=4 here to shift this legend entry to the left ##
   #
   title(main='Lake St. Clair Time Series', outer=TRUE,cex.main=1.75);
   
   #### NBS ####
   yPLimit = range(c(yLkStClNBS1,
                     #yLkStClNBS2,
                     yLkStClNBS3,yLkStClNBS4,yLkStClNBS5,yLkStClNBS6),na.rm=T)
   plot(c(0), c(0), type = "n", col = 8, lwd = 3, axes = FALSE, xlim = c(xlim1,xlim2), ylim = yPLimit*1.1); 
   box()
   # axis(2, cex.axis = 0.9); 
   axis(4, labels=T, cex.axis = 0.9); 
   axis(2, labels=FALSE, cex.axis = 0.9); 
   mtext(paste("NBS (mm)"), side = 2, line = 2.5, cex=0.8); 
   axis(1, at=seq(0,nMonths,12), labels=FALSE);
   #
   ypol = c(jSum[paste('LkStClNBS[',1:nMonths,']', sep=''),3],
            jSum[paste('LkStClNBS[',nMonths:1,']', sep=''),7])
   polygon(x = xpol,y = ypol,col="grey80",border="grey80")
   #
   lines(yLkStClNBS1, col = "red",  lwd = 1)
   #lines(yLkStClNBS2, col = "blue",  lwd = 1,lty=2)
   lines(yLkStClNBS3, col = "green4",  lwd = 1,lty=4)
   lines(yLkStClNBS4, col = "navy",  lwd = 1,lty=5)
   lines(yLkStClNBS5, col = "magenta",  lwd = 1)
   lines(yLkStClNBS6, col = "steelblue",  lwd = 1,lty=2)
   legend("topleft",lty=c(1,4,5,1,2), #lty=c(1,2,4,5,1,2),
          bty="n",legend = c("NOAA GLERL GLM HMD",
                                                          #"GLERL AHPS PROVISIONAL",
                                                        "USACE AHPS","ECCC CaPA","ECCC WATFLOOD"),
          ncol=4,cex=0.9,col = c("red",
                                 #"blue",
                                 "green4","navy","magenta","steelblue"))   
   legend("top",fill="grey80",bty="n",legend = c("L2SWBM.95CI"),ncol=2,cex=0.9,border="grey80") ## ncol=4 here to shift this legend entry to the left ##
   
   #### Precip, Evap, Runoff ####
   plot(c(0), c(0), type = "n", col = 8, lwd = 3, axes = FALSE);
   box()
   axis(2, labels=FALSE, cex.axis = 0.9);
   # axis(2, at=seq(0,700,100), labels=FALSE, cex.axis = 0.9);
   axis(4, labels=FALSE, cex.axis = 0.9);
   mtext(paste("Precip, Evap, Runoff"), side = 2, line = 2.5, cex=0.8);
   axis(1, at=seq(0,nMonths,12), labels=FALSE);
   #
   legend("topleft",lty=1,lwd=0,bty="n",legend = c("Note: Lake St. Clair does not use Precip, Evap, Runoff. NBS is used instead."),ncol=1,cex=0.9)
   
   #### Precip, Evap, Runoff ####
   plot(c(0), c(0), type = "n", col = 8, lwd = 3, axes = FALSE);
   box()
   axis(2, labels=FALSE, cex.axis = 0.9);
   # axis(2, at=seq(0,700,100), labels=FALSE, cex.axis = 0.9);
   axis(4, labels=FALSE, cex.axis = 0.9);
   mtext(paste(""), side = 2, line = 2.5, cex=0.8);
   axis(1, at=seq(0,nMonths,12), labels=FALSE);
   #
   legend("topleft",lty=1,lwd=0,bty="n",legend = c("Note: Lake St. Clair does not use Precip, Evap, Runoff. NBS is used instead."),ncol=1,cex=0.9)

   #### Outflow ####
   # yLkStClSFDOutflow = jSum[paste('SFDc'),1]*((((yLkStClOutFtGr[217:444]+yLkStClOutDryD[217:444])/2)-167)^jSum[paste('SFDd'),1])*((yLkStClOutFtGr[217:444] - yLkStClOutDryD[217:444])^jSum[paste('SFDf'),1])
   yQLimit = range(c(3500,6500
      #yLkStClOut1,yLkStClOut2
                     #,yLkStClOut3
                     ),na.rm=T)
   yQDiff = diff(yQLimit)
   # yQLimit = range(c(yLkStClSFDOut1),na.rm=T)
   # yQLimit = c(0,8000)
   plot(c(0), c(0), type = "n", col = 8, lwd = 3, axes = FALSE, xlim = c(xlim1,xlim2), ylim = c(yQLimit[1]-0.05*yQDiff,yQLimit[2]+0.05*yQDiff));#c(yQLimit[1]*0.9,yQLimit[2]*1.1));
   box()
   axis(2, cex.axis = 0.9);
   # axis(2, at=seq(0,700,100), labels=FALSE, cex.axis = 0.9);
   axis(4, labels=FALSE, cex.axis = 0.9);
   mtext(paste("Outflow (cms)"), side = 2, line = 2.5, cex=0.8);
   axis(1, at=seq(0,nMonths,12), labels=FALSE);
   #
   ypol = c(jSum[paste('LkStClOutflow[',1:nMonths,']', sep=''),3], ## 2.5% value ##
            jSum[paste('LkStClOutflow[',nMonths:1,']', sep=''),7]) ## 97.5% value ##
   polygon(x = xpol,y = ypol,col="grey80",border="grey80")
   #
   # ypol2 = c(jSum[paste('OutflowSFDMeanLkStCl[',1:nMonths,']', sep=''),3], ## 2.5% value ##
   #          jSum[paste('OutflowSFDMeanLkStCl[',nMonths:1,']', sep=''),7]) ## 97.5% value ##
   # polygon(x = xpol,y = ypol2,border="mediumpurple2")
   #
   # lines(yLkStClOutflow, col = "red",  lwd = 1)
   lines(yLkStClOut1, col = "red",  lwd = 1,lty=1)
   lines(yLkStClOut2, col = "blue",  lwd = 1,lty=2)
   #lines(yLkStClOut3, col = "green4",  lwd = 1,lty=4)
   lines(predStc, col = "navy", lwd = 1,lty=5)
   # lines(yLkStClSFDOutflow,col = "seagreen4",lwd = 1, lty=2) ## Paired with legend code on next line! If this is uncommented, uncomment legend too! ##
   # legend("topleft",lty=c(1,2,2,1),bty="n",legend = c("IGS","SFD.ADVM","SFD from model"),ncol=3,cex=0.9,col = c("red","blue","seagreen4"))
   legend("topleft",lty=c(1,2,5), #lty=c(1,2,4,5),
          bty="n",legend = c("IGS","SFD.ADVM",#"JJ21",
                                                      "LM Fit
                                                      "),ncol=2,cex=0.9,col = c("red","blue",
                                                                                #"green4",
                                                                                "navy")) ## Comment out if above legend code is active ##
   legend("top",fill="grey80",bty="n",legend = c("L2SWBM.95CI.LkStClOut"),ncol=2,cex=0.9,border="grey80") ## ncol=4 here to shift this legend entry to the left ##
   # legend("topright",fill="mediumpurple2",bty="n",legend = c("L2SWBM.95CI.OutSFD"),ncol=3,cex=0.9,border="mediumpurple2") ## ncol=4 here to shift this legend entry to the left ##
   
   #### Dstore ####
   yDSLimit = range(yLkStClDStore,na.rm=T)
   plot(yLkStClDStore, type = "n", col = "black", lwd = 2, axes = FALSE, ylim = c(yDSLimit[1],yDSLimit[2]*1.1), xlim = c(xlim1,xlim2)); 
   box()
   abline(h=0, col = 8)
   axis(4,cex.axis = 0.9);
   axis(2, labels=FALSE, cex.axis = 0.9); 
   axis(1, at=seq(0,nMonths,12), labels=FALSE);
   axis(1, at=seq(0,nMonths-1,12)+6, labels=yearRange, tick=FALSE);
   mtext('DStore (mm)', side = 2, line = 2.5, cex = 0.8)
   #
   ypol = c(jSum[paste('LkStClDStore[',1:nMonths,']', sep=''),3], ## Reads in the 2.5% column. ##
            jSum[paste('LkStClDStore[',nMonths:1,']', sep=''),7]) ## Reads in the 97.5% column. Between these is 95%!! ##
   polygon(x = xpol,y = ypol,col="grey80",border="grey80")
   #
   lines(yLkStClDStore, col = "red",  lwd = 1,lty=1)
   legend("topleft",lty=1,bty="n",legend = c("Observed DStore"),ncol=2,cex=0.9,col = c("red","blue"))
   legend("top",fill="grey80",bty="n",legend = c("Inferred DStore"),ncol=7,cex=0.9,border="grey80") ## ncol=7 here to shift this legend entry to the left ##
   
}

dev.off();





## ## ## ## ## ## ## ## ## ##
#### 7. Time Series Erie ####
## ## ## ## ## ## ## ## ## ##

pdf(paste(fileNameBase, '_TimeSeries_Erie.pdf', sep=''), width = 10, height = 7.5);
par(mfrow=c(6,1)) # HD adjusted to make all panels fall inside one page 
## COMMENT ON 17 FEB: PLEASE NOTE THAT WE PLOT OUT THE TOTAL NUMBER OF SIX VARIABLES
##                    THUS YOU SHOULD SET mfrow = c(6,1)
##                    IF DIVERSION IS ALSO PLOTTED, SET mfrow=c(7,1) AND SO ON
par(mar = c(0,0,0,0))
par(oma = c(4,4,4,4))

## ADDED ON FEB 16, 2021 ##
## to separate time series plots into several 5-year blocks ##
nYears = length(yearRange)
intv   = 5
nPlots = ceiling(nYears/intv)

for (ii in 1:nPlots) {
   xlim1 = (ii-1)*intv*12+1
   xlim2 = ii*intv*12
   ##### Inflow #####
   yILimit = range(c(yEriInfl1,yEriInfl2
                     #,yEriInfl3
                     ),na.rm=T)
   yIDiff = diff(yILimit)
   plot(c(0), c(0), type = "n", col = 8, lwd = 3, axes = FALSE, xlim = c(xlim1,xlim2), ylim = c(yILimit[1]-0.05*yIDiff,yILimit[2]+0.05*yIDiff));#c(yILimit[1]*0.9,yILimit[2]*1.1));
   box()
   axis(2, cex.axis = 0.9);
   # axis(2, at=seq(0,700,100), labels=FALSE, cex.axis = 0.9);
   axis(4, labels=FALSE, cex.axis = 0.9);
   mtext(paste("Inflow (cms)"), side = 2, line = 2.5, cex=0.8);
   axis(1, at=seq(0,nMonths,12), labels=FALSE);
   #
   xpol = c(1:nMonths,nMonths:1)
   ypol = c(jSum[paste('LkStClOutflow[',1:nMonths,']', sep=''),3],
            jSum[paste('LkStClOutflow[',nMonths:1,']', sep=''),7])
   polygon(x = xpol,y = ypol,col="grey80",border="grey80")
   #
   lines(yEriInfl1, col = "red",  lwd = 1)
   lines(yEriInfl2, col = "blue",  lwd = 1,lty=2)
   #lines(yEriInfl3, col = "green4",  lwd = 1,lty=3)
   legend("topleft",lty=c(1,2,3),bty="n",legend = c("IGS","Flow.Accounting"
                                                    #,"JJ21"
                                                    ),ncol=2,cex=0.9,col = c("red","blue")) #,"green4"
   legend("top",fill="grey80",bty="n",legend = c("L2SWBM.95CI"),ncol=4,cex=0.9,border="grey80") ## ncol=4 here to shift this legend entry to the left ##
   #
   title(main='Lake Erie Time Series', outer=TRUE,cex.main=1.75);
   
   #### Precip ####
   yPLimit = range(c(yEriPrecip1,
                     #yEriPrecip2,
                     yEriPrecip3,
                     yEriPrecip4,yEriPrecip5,yEriPrecip6,
                     yEriPrecip7,yEriPrecip8,yEriPrecip9),na.rm=T)
   plot(c(0), c(0), type = "n", col = 8, lwd = 3, axes = FALSE, xlim = c(xlim1,xlim2), ylim = yPLimit*1.1); 
   box()
   # axis(2, cex.axis = 0.9); 
   axis(4, labels=T, cex.axis = 0.9); 
   axis(2, labels=FALSE, cex.axis = 0.9); 
   mtext(paste("Precip (mm)"), side = 2, line = 2.5, cex=0.8); 
   axis(1, at=seq(0,nMonths,12), labels=FALSE);
   #
   ypol = c(jSum[paste('EriPrecip[',1:nMonths,']', sep=''),3],
            jSum[paste('EriPrecip[',nMonths:1,']', sep=''),7])
   polygon(x = xpol,y = ypol,col="grey80",border="grey80")
   #
   lines(yEriPrecip1, col = "red",  lwd = 1)
   #lines(yEriPrecip2, col = "blue",  lwd = 1,lty=2)
   lines(yEriPrecip3, col = "green4", lwd = 1,lty=4)
   lines(yEriPrecip4, col = "navy",  lwd = 1,lty=5)
   lines(yEriPrecip5, col = "magenta",  lwd = 1)
   lines(yEriPrecip6, col = "steelblue",  lwd = 1,lty=2)
   lines(yEriPrecip7, col = "deeppink",  lwd = 1,lty=4)
   lines(yEriPrecip8, col = "purple",  lwd = 1,lty=5)
   lines(yEriPrecip9, col = "dodgerblue",  lwd = 1)
   legend("topleft",lty=c(1,4,5,1,2,4,5,1), #lty=c(1,2,4,5,1,2,4,5,1),
          bty="n",
          legend = c("NOAA GLERL GLM HMD",
                     #"GLERL AHPS PROVISIONAL",
                     "USACE AHPS","ECCC WCPS",
                     "ECCC CaPA","NWS MPE","HISTORICAL COORDINATED","USACE Thiessen",
                     "MERGED MPE CaPA"),ncol=4,cex=0.9,col = c("red",
                                                               #"blue",
                                                               "green4","navy",
                                                               "magenta","steelblue","deeppink",
                                                               "purple","dodgerblue"))
   legend("top",fill="grey80",bty="n",legend = c("L2SWBM.95CI"),ncol=4,cex=0.9,border="grey80") ## ncol=4 here to shift this legend entry to the left ##
   
   #### Evap ####
   yELimit = range(c(yEriEvap1,#yEriEvap2,
                     yEriEvap3,yEriEvap4,yEriEvap5),na.rm=T)
   plot(c(0), c(0), type = "n", col = 8, lwd = 3, axes = FALSE, xlim = c(xlim1,xlim2), ylim = yELimit*1.1); 
   box()
   # axis(2, cex.axis = 0.9); 
   axis(2, labels=T, cex.axis = 0.9); 
   axis(4, labels=FALSE, cex.axis = 0.9); 
   mtext(paste("Evap (mm)"), side = 2, line = 2.5, cex=0.8); 
   axis(1, at=seq(0,nMonths,12), labels=FALSE);
   #
   ypol = c(jSum[paste('EriEvap[',1:nMonths,']', sep=''),3],
            jSum[paste('EriEvap[',nMonths:1,']', sep=''),7])
   polygon(x = xpol,y = ypol,col="grey80",border="grey80")
   #
   lines(yEriEvap1, col = "red",  lwd = 1)
   #lines(yEriEvap2, col = "blue",  lwd = 1,lty=2)
   lines(yEriEvap3, col = "green4", lwd = 1,lty=4)
   lines(yEriEvap4, col = "navy",  lwd = 1,lty=5)
   lines(yEriEvap5, col = "magenta",  lwd = 1)
   legend("topleft",lty=c(1,4,5,1), #lty=c(1,2,4,5,1),
          bty="n",legend = c("NOAA GLERL GLM HMD",
                                                        #"GLM HMD PROVISIONAL",
                                                        "USACE AHPS","ECCC WCPS","GLERL FVCOM"),
          ncol=3,cex=0.9,col = c("red",#"blue",
                                 "green4","navy","magenta"))
   legend("top",fill="grey80",bty="n",legend = c("L2SWBM.95CI"),ncol=4,cex=0.9,border="grey80") ## ncol=4 here to shift this legend entry to the left ##
   
   #### Runoff ####
   yRLimit = range(c(yEriRunoff1,#yEriRunoff2,
                     yEriRunoff3,yEriRunoff4,yEriRunoff5),na.rm=T)
   plot(c(0), c(0), type = "n", col = 8, lwd = 3, axes = FALSE, xlim = c(xlim1,xlim2), ylim = yRLimit*1.1); 
   box()
   axis(4, cex.axis = 0.9); 
   # axis(2, at=seq(0,700,100), labels=FALSE, cex.axis = 0.9); 
   axis(2, labels=FALSE, cex.axis = 0.9); 
   mtext(paste("Runoff (mm)"), side = 2, line = 2.5, cex=0.8); 
   axis(1, at=seq(0,nMonths,12), labels=FALSE);
   #
   ypol = c(jSum[paste('EriRunoff[',1:nMonths,']', sep=''),3],
            jSum[paste('EriRunoff[',nMonths:1,']', sep=''),7])
   polygon(x = xpol,y = ypol,col="grey80",border="grey80")
   #
   lines(yEriRunoff1, col = "red",  lwd = 1)
   #lines(yEriRunoff2, col = "blue",  lwd = 1,lty=2)
   lines(ySupRunoff3, col = "green4", lwd = 1,lty=4)
   lines(ySupRunoff4, col = "navy",  lwd = 1,lty=5)
   lines(ySupRunoff5, col = "magenta",  lwd = 1)
   legend("topleft",lty=c(1,4,5,1), #lty=c(1,2,4,5,1),
          bty="n",legend = c("NOAA GLERL GLM HMD",#"GLM HMD PROVISIONAL",
                                                        "USACE AHPS","ECCC WCPS","ECCC WATFLOOD"),
          ncol=4,cex=0.9,col = c("red",#"blue",
                                 "green4","navy","magenta"))
   legend("top",fill="grey80",bty="n",legend = c("L2SWBM.95CI"),ncol=4,cex=0.9,border="grey80") ## ncol=4 here to shift this legend entry to the left ##
   
   #### Outflow ####
   # yEriSFDOutflow = jSum[paste('SFDc'),1]*((((yEriOutFtGr[217:444]+yEriOutDryD[217:444])/2)-167)^jSum[paste('SFDd'),1])*((yEriOutFtGr[217:444] - yEriOutDryD[217:444])^jSum[paste('SFDf'),1])
   yQLimit = range(c(yEriOut1),na.rm=T) #,yEriOut2,yEriOut3
   yQDiff = diff(yQLimit)
   # yQLimit = range(c(yEriSFDOut1),na.rm=T)
   # yQLimit = c(0,8000)
   plot(c(0), c(0), type = "n", col = 8, lwd = 3, axes = FALSE, xlim = c(xlim1,xlim2), ylim = c(yQLimit[1]-0.05*yQDiff,yQLimit[2]+0.05*yQDiff));#c(yQLimit[1]*0.9,yQLimit[2]*1.1));
   box()
   axis(2, cex.axis = 0.9);
   # axis(2, at=seq(0,700,100), labels=FALSE, cex.axis = 0.9);
   axis(4, labels=FALSE, cex.axis = 0.9);
   mtext(paste("Outflow (cms)"), side = 2, line = 2.5, cex=0.8);
   axis(1, at=seq(0,nMonths,12), labels=FALSE);
   #
   ypol = c(jSum[paste('EriOutflow[',1:nMonths,']', sep=''),3], ## 2.5% value ##
            jSum[paste('EriOutflow[',nMonths:1,']', sep=''),7]) ## 97.5% value ##
   polygon(x = xpol,y = ypol,col="grey80",border="grey80")
   #
   # ypol2 = c(jSum[paste('OutflowSFDMeanEri[',1:nMonths,']', sep=''),3], ## 2.5% value ##
   #          jSum[paste('OutflowSFDMeanEri[',nMonths:1,']', sep=''),7]) ## 97.5% value ##
   # polygon(x = xpol,y = ypol2,border="mediumpurple2")
   #
   # lines(yEriOutflow, col = "red",  lwd = 1)
   lines(yEriOut1, col = "red",  lwd = 1,lty=1)
   #lines(yEriOut2, col = "blue",  lwd = 1,lty=2) ## Exact same data values as yEriOut2!! ##
   #lines(yEriOut3, col = "green4",  lwd = 1,lty=4)
   lines(predEri, col = "navy", lwd = 1, lty = 5)
   # lines(yEriSFDOutflow,col = "seagreen4",lwd = 1, lty=2) ## Paired with legend code on next line! If this is uncommented, uncomment legend too! ##
   # legend("topleft",lty=c(1,2,2,1),bty="n",legend = c("IGS","SFD.ADVM","SFD from model"),ncol=3,cex=0.9,col = c("red","blue","seagreen4"))
   legend("topleft",lty=c(1,5), #lty=c(1,2,4,5),
          bty="n",legend = c("SFD.ADVM",
                                                      #"Coordinated","JJ21",
                                                      "LM Fit"),ncol=4,cex=0.9,col = c("red",
                                                                                       #"blue","green4",
                                                                                       "navy")) ## Comment out if above legend code is active ##
   legend("top",fill="grey80",bty="n",legend = c("L2SWBM.95CI.EriOut"),ncol=2,cex=0.9,border="grey80") ## ncol=4 here to shift this legend entry to the left ##
   # legend("topright",fill="mediumpurple2",bty="n",legend = c("L2SWBM.95CI.OutSFD"),ncol=3,cex=0.9,border="mediumpurple2") ## ncol=4 here to shift this legend entry to the left ##
   
   #### Dstore ####
   yDSLimit = range(yEriDStore,na.rm=T)
   plot(yEriDStore, type = "n", col = "black", lwd = 2, axes = FALSE, ylim = c(yDSLimit[1],yDSLimit[2]*1.1), xlim = c(xlim1,xlim2)); 
   box()
   abline(h=0, col = 8)
   axis(4,cex.axis = 0.9);
   axis(2, labels=FALSE, cex.axis = 0.9); 
   axis(1, at=seq(0,nMonths,12), labels=FALSE);
   axis(1, at=seq(0,nMonths-1,12)+6, labels=yearRange, tick=FALSE);
   mtext('DStore (mm)', side = 2, line = 2.5, cex = 0.8)
   #
   ypol = c(jSum[paste('EriDStore[',1:nMonths,']', sep=''),3], ## Reads in the 2.5% column. ##
            jSum[paste('EriDStore[',nMonths:1,']', sep=''),7]) ## Reads in the 97.5% column. Between these is 95%!! ##
   polygon(x = xpol,y = ypol,col="grey80",border="grey80")
   #
   lines(yEriDStore, col = "red",  lwd = 1,lty=1)
   legend("topleft",lty=1,bty="n",legend = c("Observed DStore"),ncol=2,cex=0.9,col = c("red","blue"))
   legend("top",fill="grey80",bty="n",legend = c("Inferred DStore"),ncol=7,cex=0.9,border="grey80") ## ncol=7 here to shift this legend entry to the left ##
   
}

dev.off()



