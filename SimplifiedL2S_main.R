## ## ## ## ## ## ## ## ## ## ## ## ## ## #
####      Simplified L2SWBM model      ####
###         MAIN/WRAPPER SCRIPT         ###
###     Dani Cohn - May 17, 2021        ###
###    Jennani Jayaram - July 2021      ###
## ## ## ## ## ## ## ## ## ## ## ## ## ## #

## This is a long script (1000+ lines)! ##
## If using RStudio, you can use the section finder at the bottom-left of the code window to navigate. ##
## Quick location guide:    Set up & data read-in ~line 20    Prior ranges ~line 305      Channel flow config ~line 520 ##
## Prior vectors ~line 550     Loop for prior vectors ~line 610     DStore/RStore likelihood functions ~line 730 ##
## Other likelihood functions ~line 780     R list pass to jags (variables needed in jags script) ~line 830 ##
## Parameters to monitor (export to stats.csv) ~line 905      Running the jags model ~line 1030 ##
## Save work, stats, plots ~line 1075 ##

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## #
#### 1. Set up working directory and data ####
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## #
# library packages & setwd
library(rjags);
setwd("~/Desktop/CGLGP_CIA/L2SWBM/Sarahs_L2S/")

## Add-in by Dani ##
##### Config setup (workaround) ####
iters = as.integer(as.vector(200000))
halfIters = iters/2;
## End add-in ##
##### Set output file(s) name #####
outputname = 'June2_trial_200k'


##### Read in data and clean up #####
## To add in an additional dataset, you need to add a column(s) to the notated lines ##
##### Superior #####
superiorBOM = read.csv('input/SUP_BOM_MM.csv', TRUE, skip=6);
superiorBOM$date <- as.Date(paste(superiorBOM$Year, superiorBOM$Month, "01", sep="-"), format = "%Y-%m-%d");
superiorDS = cbind(
  superiorBOM[-(nrow(superiorBOM)),c(1:2,5)],
  (superiorBOM[-1,3] - superiorBOM[-(nrow(superiorBOM)),3])*1000
);

superiorOutflow = read.csv('input/StMarysMonthlyMeanFlows.csv', TRUE, skip=8);
superiorOutflow = superiorOutflow[,c(1,2,3,4,5)] ## DATASETS ## Let's do IGS (3) and Flow Accounting (4) datasets. 6 is new date format ##
superiorOutflow$date = as.Date(paste(superiorOutflow$Year, superiorOutflow$Month, "01", sep="-"), format = "%Y-%m-%d")
superiorOutflowSrc = c("***MODEL RAN WITHOUT OBSERVATIONS***")
if(ncol(superiorOutflow) >= 3){
  colnames(superiorOutflow) = c('Year', 'Month', colnames(superiorOutflow)[3:length(colnames(superiorOutflow))])
  superiorOutflowSrc = colnames(superiorOutflow)[3:length(colnames(superiorOutflow))]
}

superiorDiversion = read.csv('input/LongLacOgokiMonthlyMeanFlows.csv', TRUE, skip=6);
superiorDiversion$date <- as.Date(paste(superiorDiversion$Year, superiorDiversion$Month, "01", sep="-"), format = "%Y-%m-%d")
superiorDiversion[is.na(superiorDiversion[,3]),3] = 91;
superiorDiversion_Prior = superiorDiversion;
superiorDiversion = superiorDiversion[,c(1,2,3,4)] ## DATASETS ## 
superiorDiversionSrc = c("***MODEL RAN WITHOUT OBSERVATIONS***")
if(ncol(superiorDiversion) >= 3){
  colnames(superiorDiversion) = c('Year', 'Month', colnames(superiorDiversion)[3:length(colnames(superiorDiversion))])
  superiorDiversionSrc = colnames(superiorDiversion)[3:length(colnames(superiorDiversion))]
}
superiorPrecip = read.csv('input/SUP_lake_Prec.csv', TRUE, na.strings='-9999.9', skip=4);
superiorPrecip[superiorPrecip < 0] = NA;
superiorPrecip = superiorPrecip[,c(1,2,3,4,5,6,7,8,9,10,11)] ## DATASETS ## Let's do NOAA GLERL GLM HMD (3) and USACE AHPS (5) datasets. 13 is new date format
superiorPrecip$date <- as.Date(paste(superiorPrecip$Year, superiorPrecip$Month, "01", sep="-"), format = "%Y-%m-%d")
superiorPrecip_Prior = superiorPrecip
superiorPrecipSrc = c("***MODEL RAN WITHOUT OBSERVATIONS***")
if(ncol(superiorPrecip) >= 3){
  colnames(superiorPrecip) = c('Year', 'Month', colnames(superiorPrecip)[3:length(colnames(superiorPrecip))])
  superiorPrecipSrc = colnames(superiorPrecip)[3:length(colnames(superiorPrecip))]
}
superiorEvap = read.csv('input/SUP_lake_Evap.csv', TRUE, na.strings='-9999.9', skip=4);
superiorEvap[superiorEvap < -5000] = NA
superiorEvap = superiorEvap[,c(1,2,3,4,5,6,7)] ## DATASETS ## Let's do NOAA GLERL GLM HMD (3) and USACE AHPS (5) datasets. 10 is new date format
superiorEvap$date <- as.Date(paste(superiorEvap$Year, superiorEvap$Month, "01", sep="-"), format = "%Y-%m-%d")
superiorEvap_Prior = superiorEvap
superiorEvapSrc = c("***MODEL RAN WITHOUT OBSERVATIONS***")
if(ncol(superiorEvap) >= 3){
  colnames(superiorEvap) = c('Year', 'Month', colnames(superiorEvap)[3:length(colnames(superiorEvap))])
  superiorEvapSrc = colnames(superiorEvap)[3:length(colnames(superiorEvap))]
}
superiorRunoff = read.csv('input/SUP_lake_Runoff.csv', TRUE, na.strings='None', skip=4);
superiorRunoff[superiorRunoff < 0] = NA
superiorRunoff = superiorRunoff[,c(1,2,3,4,5,6,7)] ## DATASETS ## Let's do NOAA GLERL GLM HMD (3) and USACE AHPS (5) datasets. 10 is new date format
superiorRunoff$date <- as.Date(paste(superiorRunoff$Year, superiorRunoff$Month, "01", sep="-"), format = "%Y-%m-%d")
superiorRunoff_Prior = superiorRunoff
superiorRunoffSrc = c("***MODEL RAN WITHOUT OBSERVATIONS***")
if(ncol(superiorRunoff) >= 3){
  colnames(superiorRunoff) = c('Year', 'Month', colnames(superiorRunoff)[3:length(colnames(superiorRunoff))])
  superiorRunoffSrc = colnames(superiorRunoff)[3:length(colnames(superiorRunoff))]
}

##### Michigan-Huron #####
miHuronBOM = read.csv('input/MHG_BOM_MM.csv', TRUE, skip=6);
miHuronBOM$date <- as.Date(paste(miHuronBOM$Year, miHuronBOM$Month, "01", sep="-"), format = "%Y-%m-%d");
miHuronDS = cbind(
  miHuronBOM[-(nrow(miHuronBOM)),c(1:2,5)],
  (miHuronBOM[-1,3] - miHuronBOM[-(nrow(miHuronBOM)),3])*1000
);

miHuronOutflow = read.csv('input/StClairMonthlyMeanFlows.csv', TRUE, skip=8);
miHuronOutflow = miHuronOutflow[,c(1,2,3,4,5)] ## DATASETS ## Let's do IGS (3) and SFD.ADVM (4) datasets. 6 is new date format
miHuronOutflow$date <- as.Date(paste(miHuronOutflow$Year, miHuronOutflow$Month, "01", sep="-"), format = "%Y-%m-%d")
miHuronOutflowSrc = c("***MODEL RAN WITHOUT OBSERVATIONS***")
if(ncol(miHuronOutflow) >= 3){
  colnames(miHuronOutflow) = c('Year', 'Month', colnames(miHuronOutflow)[3:length(colnames(miHuronOutflow))])
  miHuronOutflowSrc = colnames(miHuronOutflow)[3:length(colnames(miHuronOutflow))]
}

miHuronDiversion = read.csv('input/ChicagoMonthlyMeanFlows.csv', TRUE, skip=6);
miHuronDiversion$date <- as.Date(paste(miHuronDiversion$Year, miHuronDiversion$Month, "01", sep="-"), format = "%Y-%m-%d")
miHuronDiversion[is.na(miHuronDiversion[,3]),3] = 91;
miHuronDiversion_Prior = miHuronDiversion;
miHuronDiversion = miHuronDiversion[,c(1,2,3,4)] ## DATASETS ##
miHuronDiversionSrc = c("***MODEL RAN WITHOUT OBSERVATIONS***")
if(ncol(miHuronDiversion) >= 3){
  colnames(miHuronDiversion) = c('Year', 'Month', colnames(miHuronDiversion)[3:length(colnames(miHuronDiversion))])
  miHuronDiversionSrc = colnames(miHuronDiversion)[3:length(colnames(miHuronDiversion))]
}
miHuronPrecip = read.csv('input/MHG_lake_Prec.csv', TRUE, na.strings='-9999.9', skip=4);
miHuronPrecip[miHuronPrecip < 0] = NA;
miHuronPrecip = miHuronPrecip[,c(1,2,3,4,5,6,7,8,9,10,11)] ## DATASETS ## Let's do NOAA GLERL GLM HMD (3) and USACE AHPS (5) datasets. 13 is new date format
miHuronPrecip$date <- as.Date(paste(miHuronPrecip$Year, miHuronPrecip$Month, "01", sep="-"), format = "%Y-%m-%d")
miHuronPrecip_Prior = miHuronPrecip
miHuronPrecipSrc = c("***MODEL RAN WITHOUT OBSERVATIONS***")
if(ncol(miHuronPrecip) >= 3){
  colnames(miHuronPrecip) = c('Year', 'Month', colnames(miHuronPrecip)[3:length(colnames(miHuronPrecip))])
  miHuronPrecipSrc = colnames(miHuronPrecip)[3:length(colnames(miHuronPrecip))]
}
miHuronEvap = read.csv('input/MHG_lake_Evap.csv', TRUE, na.strings='-9999.9', skip=4);
miHuronEvap[miHuronEvap < -5000] = NA
miHuronEvap = miHuronEvap[,c(1,2,3,4,5,6,7)] ## DATASETS ## Let's do NOAA GLERL GLM HMD (3) and USACE AHPS (5) datasets. 10 is new date format
miHuronEvap$date <- as.Date(paste(miHuronEvap$Year, miHuronEvap$Month, "01", sep="-"), format = "%Y-%m-%d")
miHuronEvap_Prior = miHuronEvap
miHuronEvapSrc = c("***MODEL RAN WITHOUT OBSERVATIONS***")
if(ncol(miHuronEvap) >= 3){
  colnames(miHuronEvap) = c('Year', 'Month', colnames(miHuronEvap)[3:length(colnames(miHuronEvap))])
  miHuronEvapSrc = colnames(miHuronEvap)[3:length(colnames(miHuronEvap))]
}
miHuronRunoff = read.csv('input/MHG_lake_Runoff.csv', TRUE, na.strings='None', skip=4);
miHuronRunoff[miHuronRunoff < 0] = NA
miHuronRunoff = miHuronRunoff[,c(1,2,3,4,5,6,7)] ## DATASETS ## Let's do NOAA GLERL GLM HMD (3) and USACE AHPS (5) datasets. 10 is new date format
miHuronRunoff$date <- as.Date(paste(miHuronRunoff$Year, miHuronRunoff$Month, "01", sep="-"), format = "%Y-%m-%d")
miHuronRunoff_Prior = miHuronRunoff
miHuronRunoffSrc = c("***MODEL RAN WITHOUT OBSERVATIONS***")
if(ncol(miHuronRunoff) >= 3){
  colnames(miHuronRunoff) = c('Year', 'Month', colnames(miHuronRunoff)[3:length(colnames(miHuronRunoff))])
  miHuronRunoffSrc = colnames(miHuronRunoff)[3:length(colnames(miHuronRunoff))]
}

miHuronInflow = read.csv('input/StMarysMonthlyMeanFlows.csv', TRUE, skip=8);
miHuronInflow = miHuronInflow[,c(1,2,3,4,5)] ## DATASETS ## Let's do IGS (3) and Flow Accounting (4) datasets. 6 is new date format
miHuronInflow$date <- as.Date(paste(miHuronInflow$Year, miHuronInflow$Month, "01", sep="-"), format = "%Y-%m-%d")
miHuronInflow_Prior = miHuronInflow; 
miHuronInflowSrc = c("***MODEL RAN WITHOUT OBSERVATIONS***")
if(ncol(miHuronInflow) >= 3){
  colnames(miHuronInflow) = c('Year', 'Month', colnames(miHuronInflow)[3:length(colnames(miHuronInflow))])
  miHuronInflowSrc = colnames(miHuronInflow)[3:length(colnames(miHuronInflow))]
}

##### LAKE St. Clair #####
stClairBOM = read.csv('input/STC_BOM_MM.csv', TRUE, skip=6);
stClairBOM$date <- as.Date(paste(stClairBOM$Year, stClairBOM$Month, "01", sep="-"), format = "%Y-%m-%d");
stClairDS = cbind(
  stClairBOM[-(nrow(stClairBOM)),c(1:2,5)],
  (stClairBOM[-1,3] - stClairBOM[-(nrow(stClairBOM)),3])*1000
);

stClairInflow = read.csv('input/StClairMonthlyMeanFlows.csv', TRUE, skip=8);
stClairInflow = stClairInflow[,c(1,2,3,4,5)] ## DATASETS ## Let's do IGS (3) and Flow Accounting (4) datasets. 6 is new date format
stClairInflow$date <- as.Date(paste(stClairInflow$Year, stClairInflow$Month, "01", sep="-"), format = "%Y-%m-%d")
stClairInflow_Prior = stClairInflow;
stClairInflowSrc = c("***MODEL RAN WITHOUT OBSERVATIONS***")
if(ncol(stClairInflow) >= 3){
  colnames(stClairInflow) = c('Year', 'Month', colnames(stClairInflow)[3:length(colnames(stClairInflow))])
  stClairInflowSrc = colnames(stClairInflow)[3:length(colnames(stClairInflow))]
}
## Net Basin Supply is best to use for Lake St. Clair here ##
stClairNBS = read.csv('input/STC_lake_NBS.csv', TRUE, na.strings='-9999.9', skip=4);
stClairNBS[stClairNBS < -5000] = NA
stClairNBS = stClairNBS[,c(1,2,3,4,5,6,7,8)] ## DATASETS ## Let's do NOAA GLERL GLM HMD (3) and USACE AHPS (5) datasets. 12 is new date format
stClairNBS$date <- as.Date(paste(stClairNBS$Year, stClairNBS$Month, "01", sep="-"), format = "%Y-%m-%d")
stClairNBS_Prior = stClairNBS;
stClairNBSSrc = c("***MODEL RAN WITHOUT OBSERVATIONS***")
if(ncol(stClairNBS) >= 3){
  colnames(stClairNBS) = c('Year', 'Month', colnames(stClairNBS)[3:length(colnames(stClairNBS))])
  stClairNBSSrc = colnames(stClairNBS)[3:length(colnames(stClairNBS))]
}

stClairOutflow = read.csv('input/DetroitMonthlyMeanFlows.csv', TRUE, skip=8);
stClairOutflow = stClairOutflow[,c(1,2,3,4,5)] ## DATASETS ## Let's do IGS (3) and SFD.ADVM (4) datasets. 6 is new date format
stClairOutflow$date <- as.Date(paste(stClairOutflow$Year, stClairOutflow$Month, "01", sep="-"), format = "%Y-%m-%d")
stClairOutflowSrc = c("***MODEL RAN WITHOUT OBSERVATIONS***")
if(ncol(stClairOutflow) >= 3){
  colnames(stClairOutflow) = c('Year', 'Month', colnames(stClairOutflow)[3:length(colnames(stClairOutflow))])
  stClairOutflowSrc = colnames(stClairOutflow)[3:length(colnames(stClairOutflow))]
}

##### Erie #####
erieBOM = read.csv('input/ERI_BOM_MM.csv', TRUE, skip=6);
erieBOM$date <- as.Date(paste(erieBOM$Year, erieBOM$Month, "01", sep="-"), format = "%Y-%m-%d");
erieDS = cbind(
  erieBOM[-(nrow(erieBOM)),c(1:2,5)],
  (erieBOM[-1,3] - erieBOM[-(nrow(erieBOM)),3])*1000
);

erieOutflow = read.csv('input/NiagaraWellandMonthlyMeanFlows.csv', TRUE, skip=7);
erieOutflow = erieOutflow[,c(1,2,3,4)] ## DATASETS ## Let's do the SFD.ADVM (3) and Coordinated (4, exact same values as SFD.ADVM!) datasets. 5 is new date format
erieOutflow$date <- as.Date(paste(erieOutflow$Year, erieOutflow$Month, "01", sep="-"), format = "%Y-%m-%d")
erieOutflowSrc = c("***MODEL RAN WITHOUT OBSERVATIONS***")
if(ncol(erieOutflow) >= 3){
  colnames(erieOutflow) = c('Year', 'Month', colnames(erieOutflow)[3:length(colnames(erieOutflow))])
  erieOutflowSrc = colnames(erieOutflow)[3:length(colnames(erieOutflow))]
}

eriePrecip = read.csv('input/ERI_lake_Prec.csv', TRUE, na.strings='-9999.9', skip=4);
eriePrecip[eriePrecip < 0] = NA;
eriePrecip = eriePrecip[,c(1,2,3,4,5,6,7,8,9,10,11)] ## DATASETS ## Let's do NOAA GLERL GLM HMD (3) and USACE AHPS (5) datasets. 13 is new date format
eriePrecip$date <- as.Date(paste(eriePrecip$Year, eriePrecip$Month, "01", sep="-"), format = "%Y-%m-%d")
eriePrecip_Prior = eriePrecip
eriePrecipSrc = c("***MODEL RAN WITHOUT OBSERVATIONS***")
if(ncol(eriePrecip) >= 3){
  colnames(eriePrecip) = c('Year', 'Month', colnames(eriePrecip)[3:length(colnames(eriePrecip))])
  eriePrecipSrc = colnames(eriePrecip)[3:length(colnames(eriePrecip))]
}
erieEvap = read.csv('input/ERI_lake_Evap.csv', TRUE, na.strings='-9999.9', skip=4);
erieEvap[erieEvap < -5000] = NA
erieEvap = erieEvap[,c(1,2,3,4,5,6,7)] ## DATASETS ## Let's do NOAA GLERL GLM HMD (3) and USACE AHPS (5) datasets. 10 is new date format
erieEvap$date <- as.Date(paste(erieEvap$Year, erieEvap$Month, "01", sep="-"), format = "%Y-%m-%d")
erieEvap_Prior = erieEvap
erieEvapSrc = c("***MODEL RAN WITHOUT OBSERVATIONS***")
if(ncol(erieEvap) >= 3){
  colnames(erieEvap) = c('Year', 'Month', colnames(erieEvap)[3:length(colnames(erieEvap))])
  erieEvapSrc = colnames(erieEvap)[3:length(colnames(erieEvap))]
}
erieRunoff = read.csv('input/ERI_lake_Runoff.csv', TRUE, na.strings='None', skip=4);
erieRunoff[erieRunoff < 0] = NA
erieRunoff = erieRunoff[,c(1,2,3,4,5,6,7)] ## DATASETS ## Let's do NOAA GLERL GLM HMD (3) and USACE AHPS (5) datasets. 10 is new date format
erieRunoff$date <- as.Date(paste(erieRunoff$Year, erieRunoff$Month, "01", sep="-"), format = "%Y-%m-%d")
erieRunoff_Prior = erieRunoff
erieRunoffSrc = c("***MODEL RAN WITHOUT OBSERVATIONS***")
if(ncol(erieRunoff) >= 3){
  colnames(erieRunoff) = c('Year', 'Month', colnames(erieRunoff)[3:length(colnames(erieRunoff))])
  erieRunoffSrc = colnames(erieRunoff)[3:length(colnames(erieRunoff))]
}
erieInflow = read.csv('input/DetroitMonthlyMeanFlows.csv', TRUE, skip=8);
erieInflow = erieInflow[,c(1,2,3,4,5)] ## DATASETS ## Let's do IGS (3) and Flow Accounting (4) datasets. 6 is new date format
erieInflow$date <- as.Date(paste(erieInflow$Year, erieInflow$Month, "01", sep="-"), format = "%Y-%m-%d")
erieInflow_Prior = erieInflow;
erieInflowSrc = c("***MODEL RAN WITHOUT OBSERVATIONS***")
if(ncol(erieInflow) >= 3){
  colnames(erieInflow) = c('Year', 'Month', colnames(erieInflow)[3:length(colnames(erieInflow))])
  erieInflowSrc = colnames(erieInflow)[3:length(colnames(erieInflow))]
}


##### Prior Ranges #####
# set up prior range
startPriorYear = 1982;
startPriorMonth = 1;
endPriorYear = 1999;
endPriorMonth = 12;
dateAll = as.Date(miHuronBOM[,5]) # just using miHuronBOM b/c all csvs are same date range
PriorStart = which(as.numeric(format(dateAll,"%Y"))== startPriorYear & as.numeric(format(dateAll,"%m"))== startPriorMonth)
PriorEnd = which(as.numeric(format(dateAll,"%Y"))== endPriorYear & as.numeric(format(dateAll,"%m"))== endPriorMonth)
PriorRange=PriorStart:PriorEnd


###### Superior ######
SupBOMPriorRange0=which(as.numeric(format(as.Date(superiorBOM[,5]),"%Y"))==startPriorYear 
                       & as.numeric(format(as.Date(superiorBOM[,5]),"%m"))==startPriorMonth)
SupBOMPriorRange=which(as.numeric(format(as.Date(superiorBOM[,5]),"%Y"))==endPriorYear
                      & as.numeric(format(as.Date(superiorBOM[,5]),"%m"))==endPriorMonth)
PriorSupBOMData=superiorBOM[SupBOMPriorRange0:SupBOMPriorRange,]
# PriorSupBOMData[PriorSupBOMData==0]=NA
SupDivPriorRange0=which(as.numeric(format(as.Date(superiorDiversion[,4]),"%Y"))==startPriorYear 
                       & as.numeric(format(as.Date(superiorDiversion[,4]),"%m"))==startPriorMonth)
SupDivPriorRange=which(as.numeric(format(as.Date(superiorDiversion[,4]),"%Y"))==endPriorYear
                      & as.numeric(format(as.Date(superiorDiversion[,4]),"%m"))==endPriorMonth)
PriorSupDivData=superiorDiversion[SupDivPriorRange0:SupDivPriorRange,]
# PriorSupDivData[PriorSupDivData==0]=NA
SupDSPriorRange0=which(as.numeric(format(as.Date(superiorDS[,3]),"%Y"))==startPriorYear 
                      & as.numeric(format(as.Date(superiorDS[,3]),"%m"))==startPriorMonth)
SupDSPriorRange=which(as.numeric(format(as.Date(superiorDS[,3]),"%Y"))==endPriorYear
                     & as.numeric(format(as.Date(superiorDS[,3]),"%m"))==endPriorMonth)
PriorSupDSData=superiorDS[SupDSPriorRange0:SupDSPriorRange,]
# PriorSupDSData[PriorSupDSData==0]=NA
SupOutPriorRange0=which(as.numeric(format(as.Date(superiorOutflow[,6]),"%Y"))==startPriorYear
                        & as.numeric(format(as.Date(superiorOutflow[,6]),"%m"))==startPriorMonth)
SupOutPriorRange=which(as.numeric(format(as.Date(superiorOutflow[,6]),"%Y"))==endPriorYear
                       & as.numeric(format(as.Date(superiorOutflow[,6]),"%m"))==endPriorMonth)
PriorSupOutData=superiorOutflow[SupOutPriorRange0:SupOutPriorRange,]
# PriorSupOutData[PriorSupOutData==0]=NA
SupPrecPriorRange0=which(as.numeric(format(as.Date(superiorPrecip[,12]),"%Y"))==startPriorYear 
                         & as.numeric(format(as.Date(superiorPrecip[,12]),"%m"))==startPriorMonth)
SupPrecPriorRange=which(as.numeric(format(as.Date(superiorPrecip[,12]),"%Y"))==endPriorYear
                        & as.numeric(format(as.Date(superiorPrecip[,12]),"%m"))==endPriorMonth)
PriorSupPrecData=superiorPrecip[SupPrecPriorRange0:SupPrecPriorRange,]
# PriorSupPrecData[PriorSupPrecData[,3]==0,3]=NA
# PriorSupPrecData[PriorSupPrecData[,4]==0,4]=NA
SupEvapPriorRange0=which(as.numeric(format(as.Date(superiorEvap[,8]),"%Y"))==startPriorYear 
                        & as.numeric(format(as.Date(superiorEvap[,8]),"%m"))==startPriorMonth)
SupEvapPriorRange=which(as.numeric(format(as.Date(superiorEvap[,8]),"%Y"))==endPriorYear
                       & as.numeric(format(as.Date(superiorEvap[,8]),"%m"))==endPriorMonth)
PriorSupEvapData=superiorEvap[SupEvapPriorRange0:SupEvapPriorRange,]
# PriorSupEvapData[PriorSupEvapData==0]=NA
SupRunPriorRange0=which(as.numeric(format(as.Date(superiorRunoff[,8]),"%Y"))==startPriorYear 
                       & as.numeric(format(as.Date(superiorRunoff[,8]),"%m"))==startPriorMonth)
SupRunPriorRange=which(as.numeric(format(as.Date(superiorRunoff[,8]),"%Y"))==endPriorYear
                      & as.numeric(format(as.Date(superiorRunoff[,8]),"%m"))==endPriorMonth)
PriorSupRunData=superiorRunoff[SupRunPriorRange0:SupRunPriorRange,]
# PriorSupRunData[PriorSupRunData==0]=NA

###### Michigan-Huron ######
### New GLs code (based on original code) ###
MHBOMPriorRange0=which(as.numeric(format(as.Date(miHuronBOM[,5]),"%Y"))==startPriorYear 
                  & as.numeric(format(as.Date(miHuronBOM[,5]),"%m"))==startPriorMonth)
MHBOMPriorRange=which(as.numeric(format(as.Date(miHuronBOM[,5]),"%Y"))==endPriorYear
                 & as.numeric(format(as.Date(miHuronBOM[,5]),"%m"))==endPriorMonth)
PriorMHBOMData=miHuronBOM[MHBOMPriorRange0:MHBOMPriorRange,]
# PriorMHBOMData[PriorMHBOMData==0]=NA
MHDivPriorRange0=which(as.numeric(format(as.Date(miHuronDiversion[,4]),"%Y"))==startPriorYear 
                       & as.numeric(format(as.Date(miHuronDiversion[,4]),"%m"))==startPriorMonth)
MHDivPriorRange=which(as.numeric(format(as.Date(miHuronDiversion[,4]),"%Y"))==endPriorYear
                      & as.numeric(format(as.Date(miHuronDiversion[,4]),"%m"))==endPriorMonth)
PriorMHDivData=miHuronDiversion[MHDivPriorRange0:MHDivPriorRange,]
# PriorMHDivData[PriorMHDivData==0]=NA
MHDSPriorRange0=which(as.numeric(format(as.Date(miHuronDS[,3]),"%Y"))==startPriorYear 
                       & as.numeric(format(as.Date(miHuronDS[,3]),"%m"))==startPriorMonth)
MHDSPriorRange=which(as.numeric(format(as.Date(miHuronDS[,3]),"%Y"))==endPriorYear
                      & as.numeric(format(as.Date(miHuronDS[,3]),"%m"))==endPriorMonth)
PriorMHDSData=miHuronDS[MHDSPriorRange0:MHDSPriorRange,]
# PriorMHDSData[PriorMHDSData==0]=NA
MHInflPriorRange0=which(as.numeric(format(as.Date(miHuronInflow[,6]),"%Y"))==startPriorYear 
                        & as.numeric(format(as.Date(miHuronInflow[,6]),"%m"))==startPriorMonth)
MHInflPriorRange=which(as.numeric(format(as.Date(miHuronInflow[,6]),"%Y"))==endPriorYear
                       & as.numeric(format(as.Date(miHuronInflow[,6]),"%m"))==endPriorMonth)
PriorMHInflData=miHuronInflow[MHInflPriorRange0:MHInflPriorRange,]
MHOutPriorRange0=which(as.numeric(format(as.Date(miHuronOutflow[,6]),"%Y"))==1950 
                       & as.numeric(format(as.Date(miHuronOutflow[,6]),"%m"))==startPriorMonth)
MHOutPriorRange=which(as.numeric(format(as.Date(miHuronOutflow[,6]),"%Y"))==1985
                      & as.numeric(format(as.Date(miHuronOutflow[,6]),"%m"))==endPriorMonth)
PriorMHOutData=miHuronOutflow[MHOutPriorRange0:MHOutPriorRange,]
# PriorMHOutData[PriorMHOutData==0]=NA
MHPrecPriorRange0=which(as.numeric(format(as.Date(miHuronPrecip[,12]),"%Y"))==startPriorYear 
                        & as.numeric(format(as.Date(miHuronPrecip[,12]),"%m"))==startPriorMonth)
MHPrecPriorRange=which(as.numeric(format(as.Date(miHuronPrecip[,12]),"%Y"))==endPriorYear
                       & as.numeric(format(as.Date(miHuronPrecip[,12]),"%m"))==endPriorMonth)
PriorMHPrecData=miHuronPrecip[MHPrecPriorRange0:MHPrecPriorRange,]
# PriorMHPrecData[PriorMHPrecData[,3]==0,3]=NA
# PriorMHPrecData[PriorMHPrecData[,4]==0,4]=NA
MHEvapPriorRange0=which(as.numeric(format(as.Date(miHuronEvap[,8]),"%Y"))==startPriorYear 
                       & as.numeric(format(as.Date(miHuronEvap[,8]),"%m"))==startPriorMonth)
MHEvapPriorRange=which(as.numeric(format(as.Date(miHuronEvap[,8]),"%Y"))==endPriorYear
                      & as.numeric(format(as.Date(miHuronEvap[,8]),"%m"))==endPriorMonth)
PriorMHEvapData=miHuronEvap[MHEvapPriorRange0:MHEvapPriorRange,]
# PriorMHEvapData[PriorMHEvapData==0]=NA
MHRunPriorRange0=which(as.numeric(format(as.Date(miHuronRunoff[,8]),"%Y"))==startPriorYear 
                        & as.numeric(format(as.Date(miHuronRunoff[,8]),"%m"))==startPriorMonth)
MHRunPriorRange=which(as.numeric(format(as.Date(miHuronRunoff[,8]),"%Y"))==endPriorYear
                       & as.numeric(format(as.Date(miHuronRunoff[,8]),"%m"))==endPriorMonth)
PriorMHRunData=miHuronRunoff[MHRunPriorRange0:MHRunPriorRange,]
# PriorMHRunData[PriorMHRunData==0]=NA


###### LAKE St Clair ######
### New GLs code (based on original code) ###
LkStClBOMPriorRange0=which(as.numeric(format(as.Date(stClairBOM[,5]),"%Y"))==startPriorYear 
                       & as.numeric(format(as.Date(stClairBOM[,5]),"%m"))==startPriorMonth)
LkStClBOMPriorRange=which(as.numeric(format(as.Date(stClairBOM[,5]),"%Y"))==endPriorYear
                      & as.numeric(format(as.Date(stClairBOM[,5]),"%m"))==endPriorMonth)
PriorLkStClBOMData=stClairBOM[LkStClBOMPriorRange0:LkStClBOMPriorRange,]
# PriorLkStClBOMData[PriorLkStClBOMData==0]=NA
LkStClDSPriorRange0=which(as.numeric(format(as.Date(stClairDS[,3]),"%Y"))==startPriorYear 
                      & as.numeric(format(as.Date(stClairDS[,3]),"%m"))==startPriorMonth)
LkStClDSPriorRange=which(as.numeric(format(as.Date(stClairDS[,3]),"%Y"))==endPriorYear
                     & as.numeric(format(as.Date(stClairDS[,3]),"%m"))==endPriorMonth)
PriorLkStClDSData=stClairDS[LkStClDSPriorRange0:LkStClDSPriorRange,]
# PriorLkStClDSData[PriorLkStClDSData==0]=NA
LkStClOutPriorRange0=which(as.numeric(format(as.Date(stClairOutflow[,6]),"%Y"))==1950 
                       & as.numeric(format(as.Date(stClairOutflow[,6]),"%m"))==startPriorMonth)
LkStClOutPriorRange=which(as.numeric(format(as.Date(stClairOutflow[,6]),"%Y"))==1985
                      & as.numeric(format(as.Date(stClairOutflow[,6]),"%m"))==endPriorMonth)
PriorLkStClOutData=stClairOutflow[LkStClOutPriorRange0:LkStClOutPriorRange,]
# PriorLkStClOutData[PriorLkStClOutData==0]=NA
LkStClInflPriorRange0=which(as.numeric(format(as.Date(stClairInflow[,6]),"%Y"))==startPriorYear 
                        & as.numeric(format(as.Date(stClairInflow[,6]),"%m"))==startPriorMonth)
LkStClInflPriorRange=which(as.numeric(format(as.Date(stClairInflow[,6]),"%Y"))==endPriorYear
                       & as.numeric(format(as.Date(stClairInflow[,6]),"%m"))==endPriorMonth)
PriorLkStClInflData=stClairInflow[LkStClInflPriorRange0:LkStClInflPriorRange,]
LkStClNBSPriorRange0=which(as.numeric(format(as.Date(stClairNBS[,9]),"%Y"))==startPriorYear 
                           & as.numeric(format(as.Date(stClairNBS[,9]),"%m"))==startPriorMonth)
LkStClNBSPriorRange=which(as.numeric(format(as.Date(stClairNBS[,9]),"%Y"))==endPriorYear
                          & as.numeric(format(as.Date(stClairNBS[,9]),"%m"))==endPriorMonth)
PriorLkStClNBSData=stClairNBS[LkStClNBSPriorRange0:LkStClNBSPriorRange,]
# PriorLkStClNBSData[PriorLkStClNBSData==0]=NA

###### Erie ######
### New GLs code (based on original code) ###
EriBOMPriorRange0=which(as.numeric(format(as.Date(erieBOM[,5]),"%Y"))==startPriorYear 
                       & as.numeric(format(as.Date(erieBOM[,5]),"%m"))==startPriorMonth)
EriBOMPriorRange=which(as.numeric(format(as.Date(erieBOM[,5]),"%Y"))==endPriorYear
                      & as.numeric(format(as.Date(erieBOM[,5]),"%m"))==endPriorMonth)
PriorEriBOMData=erieBOM[EriBOMPriorRange0:EriBOMPriorRange,]
# PriorEriBOMData[PriorEriBOMData==0]=NA
EriDSPriorRange0=which(as.numeric(format(as.Date(erieDS[,3]),"%Y"))==startPriorYear 
                      & as.numeric(format(as.Date(erieDS[,3]),"%m"))==startPriorMonth)
EriDSPriorRange=which(as.numeric(format(as.Date(erieDS[,3]),"%Y"))==endPriorYear
                     & as.numeric(format(as.Date(erieDS[,3]),"%m"))==endPriorMonth)
PriorEriDSData=erieDS[EriDSPriorRange0:EriDSPriorRange,]
# PriorEriDSData[PriorEriDSData==0]=NA
EriInflPriorRange0=which(as.numeric(format(as.Date(erieInflow[,6]),"%Y"))==startPriorYear 
                         & as.numeric(format(as.Date(erieInflow[,6]),"%m"))==startPriorMonth)
EriInflPriorRange=which(as.numeric(format(as.Date(erieInflow[,6]),"%Y"))==endPriorYear
                        & as.numeric(format(as.Date(erieInflow[,6]),"%m"))==endPriorMonth)
PriorEriInflData=erieInflow[EriInflPriorRange0:EriInflPriorRange,]
EriOutPriorRange0=which(as.numeric(format(as.Date(erieOutflow[,5]),"%Y"))==1950 ## Change to [,4] if only one dataset used
                        & as.numeric(format(as.Date(erieOutflow[,5]),"%m"))==startPriorMonth) ## Change to [,4] if only one dataset used
EriOutPriorRange=which(as.numeric(format(as.Date(erieOutflow[,5]),"%Y"))==1985 ## Change to [,4] if only one dataset used
                       & as.numeric(format(as.Date(erieOutflow[,5]),"%m"))==endPriorMonth) ## Change to [,4] if only one dataset used
PriorEriOutData=erieOutflow[EriOutPriorRange0:EriOutPriorRange,]
# PriorEriOutData[PriorEriOutData==0]=NA
EriPrecPriorRange0=which(as.numeric(format(as.Date(eriePrecip[,12]),"%Y"))==startPriorYear 
                         & as.numeric(format(as.Date(eriePrecip[,12]),"%m"))==startPriorMonth)
EriPrecPriorRange=which(as.numeric(format(as.Date(eriePrecip[,12]),"%Y"))==endPriorYear
                        & as.numeric(format(as.Date(eriePrecip[,12]),"%m"))==endPriorMonth)
PriorEriPrecData=eriePrecip[EriPrecPriorRange0:EriPrecPriorRange,]
# PriorEriPrecData[PriorEriPrecData[,3]==0,3]=NA
# PriorEriPrecData[PriorEriPrecData[,4]==0,4]=NA
EriEvapPriorRange0=which(as.numeric(format(as.Date(erieEvap[,8]),"%Y"))==startPriorYear 
                        & as.numeric(format(as.Date(erieEvap[,8]),"%m"))==startPriorMonth)
EriEvapPriorRange=which(as.numeric(format(as.Date(erieEvap[,8]),"%Y"))==endPriorYear
                       & as.numeric(format(as.Date(erieEvap[,8]),"%m"))==endPriorMonth)
PriorEriEvapData=erieEvap[EriEvapPriorRange0:EriEvapPriorRange,]
# PriorEriEvapData[PriorEriEvapData==0]=NA
EriRunPriorRange0=which(as.numeric(format(as.Date(erieRunoff[,8]),"%Y"))==startPriorYear 
                       & as.numeric(format(as.Date(erieRunoff[,8]),"%m"))==startPriorMonth)
EriRunPriorRange=which(as.numeric(format(as.Date(erieRunoff[,8]),"%Y"))==endPriorYear
                      & as.numeric(format(as.Date(erieRunoff[,8]),"%m"))==endPriorMonth)
PriorEriRunData=erieRunoff[EriRunPriorRange0:EriRunPriorRange,]
# PriorEriRunData[PriorEriRunData==0]=NA
### End new GLs code block ###

# set up posterior range
startAnalysisYear = 2000;
startAnalysisMonth = 1;
endAnalysisYear = 2017;
endAnalysisMonth = 12;
dateAll = as.Date(miHuronBOM[,5]) # just using miHuronBOM b/c all csvs are same date range
idStart = which(as.numeric(format(dateAll,"%Y"))== startAnalysisYear & as.numeric(format(dateAll,"%m"))== startAnalysisMonth)
idEnd = which(as.numeric(format(dateAll,"%Y"))== endAnalysisYear & as.numeric(format(dateAll,"%m"))== endAnalysisMonth)
PosteriorRange=idStart:idEnd

##### Channel flow configurations (for mm to cms conversion) - added Feb 16 & May 24, 2021 #####
supArea    = 8.1925e10
mhgArea    = 1.1685e11
lkstclArea = 1.1090e9
eriArea    = 2.5404e10
daysInMonths = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31);
daysInMonthsWithLeap = c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31);
secondsInADay = 24*60*60;
dayVector = NULL;
yearVector = NULL;
m = NULL;
for(i in startAnalysisYear:endAnalysisYear){
  for(mo in 1:12){
    if(i == startAnalysisYear & mo < startAnalysisMonth){
      next;
    }
    if(i == endAnalysisYear & mo > endAnalysisMonth){
      next;
    }
    yearVector = c(yearVector, i);
    m = c(m, mo)
    if(i %% 4 == 0){
      dayVector = c(dayVector, daysInMonthsWithLeap[mo])
    }
    else{
      dayVector = c(dayVector, daysInMonths[mo])
    }
  }
}

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
#### 2. Parameter for Prior Distribution ####
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

##### Set prior vectors #####
###### Superior ######
PrecipShapeSup          =NULL
PrecipRateSup           =NULL
PrecipPrecisionSup      =NULL
EvapMeanSup             =NULL
EvapPrecisionSup        =NULL
LogRunoffMeanSup        =NULL
LogRunoffPrecisionSup   =NULL
OutflowShapeSup         =NULL
OutflowMeanSup          =NULL
OutflowShapeSup         =NULL
OutflowMeanSup          =NULL
OutflowPrecisionSup     =NULL
InflowShapeSup          =NULL
InflowMeanSup           =NULL
InflowPrecisionSup      =NULL
###### Michigan-Huron ######
PrecipShapeMH           =NULL
PrecipRateMH            =NULL
PrecipPrecisionMH       =NULL
EvapMeanMH              =NULL
EvapPrecisionMH         =NULL
LogRunoffMeanMH         =NULL
LogRunoffPrecisionMH    =NULL
OutflowShapeMH          =NULL
OutflowMeanMH           =NULL
OutflowPrecisionMH      =NULL
InflowShapeMH           =NULL
InflowMeanMH            =NULL
InflowPrecisionMH       =NULL
###### LAKE St. Clair ######
NBSMeanLkStCl           =NULL
NBSPrecisionLkStCl      =NULL
OutflowShapeLkStCl      =NULL
OutflowMeanLkStCl       =NULL
OutflowPrecisionLkStCl  =NULL
InflowShapeLkStCl       =NULL
InflowMeanLkStCl        =NULL
InflowPrecisionLkStCl   =NULL
###### Erie ######
PrecipShapeEri          =NULL
PrecipRateEri           =NULL
PrecipPrecisionEri      =NULL
EvapMeanEri             =NULL
EvapPrecisionEri        =NULL
LogRunoffMeanEri        =NULL
LogRunoffPrecisionEri   =NULL
OutflowShapeEri         =NULL
OutflowMeanEri          =NULL
OutflowPrecisionEri     =NULL
InflowShapeEri          =NULL
InflowMeanEri           =NULL
InflowPrecisionEri      =NULL

##### Assigning variable columns #####
Precipitation.Column.Number=Runoff.Column.Number=Evap.Column.Number=Outflow.Column.Number=NBS.Column.Number=Inflow.Column.Number=c(3,4)
Diversion.Column.Number = c(3)

##### Prior distribution loop #####
for (i in 1:12){
  ###### Precipitation ######
  ## Superior ##
  datPrecipSup      =PriorSupPrecData[as.numeric(format(as.Date(PriorSupPrecData[,12]),"%m"))==i,Precipitation.Column.Number]
  datPrecipSup      =as.numeric(do.call(c,datPrecipSup))
  x.bar_s          =mean(datPrecipSup,na.rm = TRUE)
  theta_s          =log(x.bar_s)-mean(log(datPrecipSup),na.rm = TRUE)
  PrecipShapeSup[i] =1/(4*theta_s)*(1+sqrt(1+4*theta_s/3))
  PrecipRateSup[i]  =PrecipShapeSup[i]/x.bar_s
  ## Michigan-Huron ##
  datPrecipMH      =PriorMHPrecData[as.numeric(format(as.Date(PriorMHPrecData[,12]),"%m"))==i,Precipitation.Column.Number]
  datPrecipMH      =as.numeric(do.call(c,datPrecipMH))
  x.bar_s          =mean(datPrecipMH,na.rm = TRUE)
  theta_s          =log(x.bar_s)-mean(log(datPrecipMH),na.rm = TRUE)
  PrecipShapeMH[i] =1/(4*theta_s)*(1+sqrt(1+4*theta_s/3))
  PrecipRateMH[i]  =PrecipShapeMH[i]/x.bar_s
  ## Erie ##
  datPrecipEri     =PriorEriPrecData[as.numeric(format(as.Date(PriorEriPrecData[,12]),"%m"))==i,Precipitation.Column.Number]
  datPrecipEri     =as.numeric(do.call(c,datPrecipEri))
  x.bar_s          =mean(datPrecipEri,na.rm = TRUE)
  theta_s          =log(x.bar_s)-mean(log(datPrecipEri),na.rm = TRUE)
  PrecipShapeEri[i]=1/(4*theta_s)*(1+sqrt(1+4*theta_s/3))
  PrecipRateEri[i] =PrecipShapeEri[i]/x.bar_s
  ###### Evaporation ######
  ## Superior ##
  EvapDatSup         =PriorSupEvapData[as.numeric(format(as.Date(PriorSupEvapData[,8]),"%m"))==i,Evap.Column.Number]
  EvapDatSup         =as.numeric(do.call(c,EvapDatSup))
  EvapMeanSup[i]     =mean(EvapDatSup,na.rm = TRUE)
  EvapPrecisionSup[i]=(1/var(EvapDatSup,na.rm = TRUE))
  ## Michigan-Huron ##
  EvapDatMH         =PriorMHEvapData[as.numeric(format(as.Date(PriorMHEvapData[,8]),"%m"))==i,Evap.Column.Number]
  EvapDatMH         =as.numeric(do.call(c,EvapDatMH))
  EvapMeanMH[i]     =mean(EvapDatMH,na.rm = TRUE)
  EvapPrecisionMH[i]=(1/var(EvapDatMH,na.rm = TRUE))
  ## Erie ##
  EvapDatEri         =PriorEriEvapData[as.numeric(format(as.Date(PriorEriEvapData[,8]),"%m"))==i,Evap.Column.Number]
  EvapDatEri         =as.numeric(do.call(c,EvapDatEri))
  EvapMeanEri[i]     =mean(EvapDatEri,na.rm = TRUE)
  EvapPrecisionEri[i]=(1/var(EvapDatEri,na.rm = TRUE))
  ###### Runoff ######
  ## Superior ##
  logRunoffDatSup         =log(PriorSupRunData[as.numeric(format(as.Date(PriorSupRunData[,8]),"%m"))==i,Runoff.Column.Number])
  logRunoffDatSup         =as.numeric(do.call(c,logRunoffDatSup))
  LogRunoffMeanSup[i]     =mean(logRunoffDatSup,na.rm = TRUE)
  LogRunoffPrecisionSup[i]=(1/var(logRunoffDatSup,na.rm = TRUE))
  ## Michigan-Huron ##
  logRunoffDatMH         =log(PriorMHRunData[as.numeric(format(as.Date(PriorMHRunData[,8]),"%m"))==i,Runoff.Column.Number])
  logRunoffDatMH         =as.numeric(do.call(c,logRunoffDatMH))
  LogRunoffMeanMH[i]     =mean(logRunoffDatMH,na.rm = TRUE)
  LogRunoffPrecisionMH[i]=(1/var(logRunoffDatMH,na.rm = TRUE))
  ## Erie ##
  logRunoffDatEri         =log(PriorEriRunData[as.numeric(format(as.Date(PriorEriRunData[,8]),"%m"))==i,Runoff.Column.Number])
  logRunoffDatEri         =as.numeric(do.call(c,logRunoffDatEri))
  LogRunoffMeanEri[i]     =mean(logRunoffDatEri,na.rm = TRUE)
  LogRunoffPrecisionEri[i]=(1/var(logRunoffDatEri,na.rm = TRUE))
  ###### Net Basin Supply (NBS) ######
  ## LAKE St. Clair ##
  NBSDatLkStCl           =PriorLkStClNBSData[as.numeric(format(as.Date(PriorLkStClNBSData[,9]),"%m"))==i,NBS.Column.Number]
  NBSDatLkStCl           =as.numeric(do.call(c,NBSDatLkStCl))
  NBSMeanLkStCl[i]       =mean(NBSDatLkStCl,na.rm = TRUE)
  NBSPrecisionLkStCl[i]  =(1/var(NBSDatLkStCl,na.rm = TRUE))
  # ###### Diversions (added in by Dani) ######
  ## Superior ##
  # DivDatSup = PriorSupDivData[as.numeric(format(as.Date(PriorSupDivData[,4]),"%m"))==i,Diversion.Column.Number]
  # DivDatSup = as.numeric(do.call(c,DivDatSup))
  # DivMeanSup[i]=mean(DivDatSup,na.rm = TRUE)
  # DivPrecisionSup[i]=(1/var(DivDatSup,na.rm = TRUE))
  # ## Michigan-Huron ##
  # DivDatMH = PriorMHDivData[as.numeric(format(as.Date(PriorMHDivData[,4]),"%m"))==i,Diversion.Column.Number]
  # DivDatMH = as.numeric(do.call(c,DivDatMH))
  # DivMeanMH[i]=mean(DivDatMH,na.rm = TRUE)
  # DivPrecisionMH[i]=(1/var(DivDatMH,na.rm = TRUE))
  ###### Outflow (added in by Dani) ######
  ## Superior ##
  OutflowDatSup         =PriorSupOutData[as.numeric(format(as.Date(PriorSupOutData[,6]),"%m"))==i,Outflow.Column.Number]
  OutflowDatSup         =as.numeric(do.call(c,OutflowDatSup))
  OutflowMeanSup[i]     =mean(OutflowDatSup,na.rm = TRUE)
  OutflowPrecisionSup[i]=(1/var(OutflowDatSup,na.rm = TRUE))
  ## Michigan-Huron ##
  OutflowDatMH          =PriorMHOutData[as.numeric(format(as.Date(PriorMHOutData[,6]),"%m"))==i,Outflow.Column.Number]
  OutflowDatMH          =as.numeric(do.call(c,OutflowDatMH))
  OutflowMeanMH[i]      =mean(OutflowDatMH,na.rm = TRUE)
  OutflowPrecisionMH[i] =(1/var(OutflowDatMH,na.rm = TRUE))
  ## LAKE St. Clair ##
  OutflowDatLkStCl          =PriorLkStClOutData[as.numeric(format(as.Date(PriorLkStClOutData[,6]),"%m"))==i,Outflow.Column.Number]
  OutflowDatLkStCl          =as.numeric(do.call(c,OutflowDatLkStCl))
  OutflowMeanLkStCl[i]      =mean(OutflowDatLkStCl,na.rm = TRUE)
  OutflowPrecisionLkStCl[i] =(1/var(OutflowDatLkStCl,na.rm = TRUE))
  ## Erie ##
  OutflowDatEri         =PriorEriOutData[as.numeric(format(as.Date(PriorEriOutData[,5]),"%m"))==i,Outflow.Column.Number] ## Change to [,4] if only one dataset used
  OutflowDatEri         =as.numeric(do.call(c,OutflowDatEri))
  OutflowMeanEri[i]     =mean(OutflowDatEri,na.rm = TRUE)
  OutflowPrecisionEri[i]=(1/var(OutflowDatEri,na.rm = TRUE))
  ###### Inflow (added in by Dani) ######
  ## Superior has no Inflow source ##
  ## Michigan-Huron ##                                                ## Need to set inflow here as Sup out?? ##
  InflowDatMH         =PriorMHInflData[as.numeric(format(as.Date(PriorMHInflData[,6]),"%m"))==i,Inflow.Column.Number]
  InflowDatMH         =as.numeric(do.call(c,InflowDatMH))
  InflowMeanMH[i]     =mean(InflowDatMH,na.rm = TRUE)
  InflowPrecisionMH[i]=(1/var(InflowDatMH,na.rm = TRUE))
  ## LAKE St. Clair ##
  InflowDatLkStCl         =PriorLkStClInflData[as.numeric(format(as.Date(PriorLkStClInflData[,6]),"%m"))==i,Inflow.Column.Number]
  InflowDatLkStCl         =as.numeric(do.call(c,InflowDatLkStCl))
  InflowMeanLkStCl[i]     =mean(InflowDatLkStCl,na.rm = TRUE)
  InflowPrecisionLkStCl[i]=(1/var(InflowDatLkStCl,na.rm = TRUE))
  ## Erie ##
  InflowDatEri         =PriorEriInflData[as.numeric(format(as.Date(PriorEriInflData[,6]),"%m"))==i,Inflow.Column.Number]
  InflowDatEri         =as.numeric(do.call(c,InflowDatEri))
  InflowMeanEri[i]     =mean(InflowDatEri,na.rm = TRUE)
  InflowPrecisionEri[i]=(1/var(InflowDatEri,na.rm = TRUE))
}

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
#### 3. Observation for likelihood function ####
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
m=as.numeric(format(as.Date(miHuronBOM[,5]),"%m"))[PosteriorRange]
RollPeriod=1

##### Fill in observation (likelihood function y) #####
# waterlevel 1
## Superior ##
superiorwaterLevel = as.numeric(superiorBOM[,3])*1000
ySupDStore=c(superiorwaterLevel[2:length(superiorwaterLevel)]-superiorwaterLevel[1:length(superiorwaterLevel)-1])
ySupRStore=rep(NA,length(ySupDStore)) # RStore = rollperiod cumulative storage
for (rp in RollPeriod:length(ySupDStore)){
  # ySupRStore[rp]=superiorwaterLevel[rp+1] - superiorwaterLevel[(rp-RollPeriod+1)]
  ySupRStore[rp]=sum(ySupDStore[(rp-RollPeriod+1):rp])
}
ySupDStore=ySupDStore[PosteriorRange]
ySupDStore[ySupDStore=='NaN']=NA
ySupRStore=ySupRStore[PosteriorRange]
ySupRStore[ySupRStore=='NaN']=NA
## Michigan-Huron ##
miHuronwaterLevel = as.numeric(miHuronBOM[,3])*1000
yMHDStore=c(miHuronwaterLevel[2:length(miHuronwaterLevel)]-miHuronwaterLevel[1:length(miHuronwaterLevel)-1])
yMHRStore=rep(NA,length(yMHDStore)) # RStore = rollperiod cumulative storage
for (rp in RollPeriod:length(yMHDStore)){
  # yMHRStore[rp]=miHuronwaterLevel[rp+1] - miHuronwaterLevel[(rp-RollPeriod+1)]
  yMHRStore[rp]=sum(yMHDStore[(rp-RollPeriod+1):rp])
}
yMHDStore=yMHDStore[PosteriorRange]
yMHDStore[yMHDStore=='NaN']=NA
yMHRStore=yMHRStore[PosteriorRange]
yMHRStore[yMHRStore=='NaN']=NA
## LAKE St. Clair ##
stClairwaterLevel = as.numeric(stClairBOM[,3])*1000
yLkStClDStore=c(stClairwaterLevel[2:length(stClairwaterLevel)]-stClairwaterLevel[1:length(stClairwaterLevel)-1])
yLkStClRStore=rep(NA,length(yLkStClDStore)) # RStore = rollperiod cumulative storage
for (rp in RollPeriod:length(yLkStClDStore)){
  # yLkStClRStore[rp]=stClairwaterLevel[rp+1] - stClairwaterLevel[(rp-RollPeriod+1)]
  yLkStClRStore[rp]=sum(yLkStClDStore[(rp-RollPeriod+1):rp])
}
yLkStClDStore=yLkStClDStore[PosteriorRange]
yLkStClDStore[yLkStClDStore=='NaN']=NA
yLkStClRStore=yLkStClRStore[PosteriorRange]
yLkStClRStore[yLkStClRStore=='NaN']=NA
## Erie ##
eriewaterLevel = as.numeric(erieBOM[,3])*1000
yEriDStore=c(eriewaterLevel[2:length(eriewaterLevel)]-eriewaterLevel[1:length(eriewaterLevel)-1])
yEriRStore=rep(NA,length(yEriDStore)) # RStore = rollperiod cumulative storage
for (rp in RollPeriod:length(yEriDStore)){
  # yEriRStore[rp]=eriewaterLevel[rp+1] - eriewaterLevel[(rp-RollPeriod+1)]
  yEriRStore[rp]=sum(yEriDStore[(rp-RollPeriod+1):rp])
}
yEriDStore=yEriDStore[PosteriorRange]
yEriDStore[yEriDStore=='NaN']=NA
yEriRStore=yEriRStore[PosteriorRange]
yEriRStore[yEriRStore=='NaN']=NA

##### Other likelihood functions #####
###### Superior ######
ySupPrecip1    =as.numeric(superiorPrecip[PosteriorRange,3])
#ySupPrecip2    =as.numeric(superiorPrecip[PosteriorRange,4])
ySupPrecip3    =as.numeric(superiorPrecip[PosteriorRange,5])
ySupPrecip4    =as.numeric(superiorPrecip[PosteriorRange,6])
ySupPrecip5    =as.numeric(superiorPrecip[PosteriorRange,7])
ySupPrecip6    =as.numeric(superiorPrecip[PosteriorRange,8])
#ySupPrecip7    =as.numeric(superiorPrecip[PosteriorRange,9])
ySupPrecip8    =as.numeric(superiorPrecip[PosteriorRange,10])
ySupPrecip9    =as.numeric(superiorPrecip[PosteriorRange,11])
ySupEvap1      =as.numeric(superiorEvap[PosteriorRange,3])
#ySupEvap2      =as.numeric(superiorEvap[PosteriorRange,4])
ySupEvap3      =as.numeric(superiorEvap[PosteriorRange,5])
ySupEvap4      =as.numeric(superiorEvap[PosteriorRange,6])
ySupEvap5      =as.numeric(superiorEvap[PosteriorRange,7])
ySupRunoff1    =as.numeric(superiorRunoff[PosteriorRange,3])
#ySupRunoff2    =as.numeric(superiorRunoff[PosteriorRange,4])
ySupRunoff3    =as.numeric(superiorRunoff[PosteriorRange,5])
ySupRunoff4    =as.numeric(superiorRunoff[PosteriorRange,6])
ySupRunoff5    =as.numeric(superiorRunoff[PosteriorRange,7])
ySupOut1       =as.numeric(superiorOutflow[PosteriorRange,3]) ## St Marys IGS ##
ySupOut2       =as.numeric(superiorOutflow[PosteriorRange,4]) ## St Marys Flow Accounting ##
#ySupOut3       =as.numeric(superiorOutflow[PosteriorRange,5]) ## JJ 2021 ##
###### Michigan-Huron ######
yMHPrecip1    =as.numeric(miHuronPrecip[PosteriorRange,3])
yMHPrecip2    =as.numeric(miHuronPrecip[PosteriorRange,4])
yMHPrecip3    =as.numeric(miHuronPrecip[PosteriorRange,5])
yMHPrecip4    =as.numeric(miHuronPrecip[PosteriorRange,6])
yMHPrecip5    =as.numeric(miHuronPrecip[PosteriorRange,7])
yMHPrecip6    =as.numeric(miHuronPrecip[PosteriorRange,8])
yMHPrecip7    =as.numeric(miHuronPrecip[PosteriorRange,9])
yMHPrecip8    =as.numeric(miHuronPrecip[PosteriorRange,10])
yMHPrecip9    =as.numeric(miHuronPrecip[PosteriorRange,11])
yMHEvap1      =as.numeric(miHuronEvap[PosteriorRange,3])
#yMHEvap2      =as.numeric(miHuronEvap[PosteriorRange,4])
yMHEvap3      =as.numeric(miHuronEvap[PosteriorRange,5])
yMHEvap4      =as.numeric(miHuronEvap[PosteriorRange,6])
yMHEvap5      =as.numeric(miHuronEvap[PosteriorRange,7])
yMHRunoff1    =as.numeric(miHuronRunoff[PosteriorRange,3])
#yMHRunoff2    =as.numeric(miHuronRunoff[PosteriorRange,4])
yMHRunoff3    =as.numeric(miHuronRunoff[PosteriorRange,5])
yMHRunoff4    =as.numeric(miHuronRunoff[PosteriorRange,6])
yMHRunoff5    =as.numeric(miHuronRunoff[PosteriorRange,7])
yMHOut1       =as.numeric(miHuronOutflow[PosteriorRange,3]) ## St Clair IGS ##
yMHOut2       =as.numeric(miHuronOutflow[PosteriorRange,4]) ## St Clair SFD.ADVM ##
#yMHOut3       =as.numeric(miHuronOutflow[PosteriorRange,5]) ## JJ 2021 ##
yMHInfl1      =as.numeric(miHuronInflow[PosteriorRange,3])
yMHInfl2      =as.numeric(miHuronInflow[PosteriorRange,4])
#yMHInfl3      =as.numeric(miHuronInflow[PosteriorRange,5]) ## JJ 2021 ##
###### LAKE St. Clair ######
yLkStClNBS1   =as.numeric(stClairNBS[PosteriorRange,3])#*c(rep(1,12),rep(NA,length(PosteriorRange)-12))
#yLkStClNBS2   =as.numeric(stClairNBS[PosteriorRange,4])#*c(rep(1,12),rep(NA,length(PosteriorRange)-12))
yLkStClNBS3   =as.numeric(stClairNBS[PosteriorRange,5])#*c(rep(1,12),rep(NA,length(PosteriorRange)-12))
yLkStClNBS4   =as.numeric(stClairNBS[PosteriorRange,6])#*c(rep(1,12),rep(NA,length(PosteriorRange)-12))
yLkStClNBS5   =as.numeric(stClairNBS[PosteriorRange,7])#*c(rep(1,12),rep(NA,length(PosteriorRange)-12))
yLkStClNBS6   =as.numeric(stClairNBS[PosteriorRange,8])#*c(rep(1,12),rep(NA,length(PosteriorRange)-12))
#sk   yLkStClOut1   =as.numeric(stClairOutflow[PosteriorRange,3]) ## Detroit IGS ##
#sk   yLkStClOut2   =as.numeric(stClairOutflow[PosteriorRange,4]) ## Detroit SFD.ADVM ##
#yLkStClOut3   =as.numeric(stClairOutflow[PosteriorRange,5]) ## JJ 2021 ##
yLkStClInfl1  =as.numeric(stClairInflow[PosteriorRange,3])
yLkStClInfl2  =as.numeric(stClairInflow[PosteriorRange,4])
#yLkStClInfl3  =as.numeric(stClairInflow[PosteriorRange,5]) ## JJ 2021 ##
###### Erie ######
#Noise = rnorm(length(PosteriorRange),mean=0,sd=1000)
yEriPrecip1   =as.numeric(eriePrecip[PosteriorRange,3])
#yEriPrecip2   =as.numeric(eriePrecip[PosteriorRange,4])
yEriPrecip3   =as.numeric(eriePrecip[PosteriorRange,5])
yEriPrecip4   =as.numeric(eriePrecip[PosteriorRange,6])
yEriPrecip5   =as.numeric(eriePrecip[PosteriorRange,7])
yEriPrecip6   =as.numeric(eriePrecip[PosteriorRange,8])
yEriPrecip7   =as.numeric(eriePrecip[PosteriorRange,9])
yEriPrecip8   =as.numeric(eriePrecip[PosteriorRange,10])
yEriPrecip9   =as.numeric(eriePrecip[PosteriorRange,11])
yEriEvap1     =as.numeric(erieEvap[PosteriorRange,3])
#yEriEvap2     =as.numeric(erieEvap[PosteriorRange,4])
yEriEvap3     =as.numeric(erieEvap[PosteriorRange,5])
yEriEvap4     =as.numeric(erieEvap[PosteriorRange,6])
yEriEvap5     =as.numeric(erieEvap[PosteriorRange,7])
yEriRunoff1   =as.numeric(erieRunoff[PosteriorRange,3])
#yEriRunoff2   =as.numeric(erieRunoff[PosteriorRange,4])
yEriRunoff3   =as.numeric(erieRunoff[PosteriorRange,5])
yEriRunoff4   =as.numeric(erieRunoff[PosteriorRange,6])
yEriRunoff5   =as.numeric(erieRunoff[PosteriorRange,7])
yEriOut1      =as.numeric(erieOutflow[PosteriorRange,3]) ## Niagara SFD.ADVM ##
#yEriOut2      =as.numeric(erieOutflow[PosteriorRange,4])+NA ## Niagara Coordinated - EXACT same values as SFD.ADVM data! ##
#yEriOut3      =as.numeric(erieOutflow[PosteriorRange,5]) ## JJ 2021 ##
yEriInfl1     =as.numeric(erieInflow[PosteriorRange,3])
yEriInfl2     =as.numeric(erieInflow[PosteriorRange,4])
#yEriInfl3     =as.numeric(erieInflow[PosteriorRange,5])

## ## ## ## ## ## ## ## 
#### 4. Call to JAGS ####
## ## ## ## ## ## ## ## 

PosteriorStartMonth=idStart - (idStart-1);
PosteriorEndMonth=idEnd - (idStart-1);

##### R list pass to JAGS #####
inputDataCoreJAGS=list(
  "PosteriorStartMonth" =PosteriorStartMonth,
  "PosteriorEndMonth"   =PosteriorEndMonth,
  "RollPeriod"          =RollPeriod,
  "m"                   =m,
  "dayVector"           =dayVector,
  "secondsInADay"       =secondsInADay,
  "supArea"             =supArea,
  "mhgArea"             =mhgArea,
  "lkstclArea"          =lkstclArea,
  "eriArea"             =eriArea,
  "ySupRStore"           =ySupRStore,
  "PrecipShapeSup"       =PrecipShapeSup,
  "PrecipRateSup"        =PrecipRateSup,
  "EvapMeanSup"          =EvapMeanSup,
  "EvapPrecisionSup"     =EvapPrecisionSup,
  "LogRunoffMeanSup"     =LogRunoffMeanSup,
  "LogRunoffPrecisionSup"=LogRunoffPrecisionSup,
  "OutflowMeanSup"       =OutflowMeanSup,
  "OutflowPrecisionSup"  =OutflowPrecisionSup,
  "ySupPrecip1"          =ySupPrecip1,
#  "ySupPrecip2"          =ySupPrecip2,
  "ySupPrecip3"          =ySupPrecip3,
  "ySupPrecip4"          =ySupPrecip4,
  "ySupPrecip5"          =ySupPrecip5,
  "ySupPrecip6"          =ySupPrecip6,
#  "ySupPrecip7"          =ySupPrecip7,
  "ySupPrecip8"          =ySupPrecip8,
  "ySupPrecip9"          =ySupPrecip9,
  "ySupEvap1"            =ySupEvap1,
#  "ySupEvap2"            =ySupEvap2,
  "ySupEvap3"            =ySupEvap3,
  "ySupEvap4"            =ySupEvap4,
  "ySupEvap5"            =ySupEvap5,
  "ySupRunoff1"          =ySupRunoff1,
#  "ySupRunoff2"          =ySupRunoff2,
  "ySupRunoff3"          =ySupRunoff3,
  "ySupRunoff4"          =ySupRunoff4,
  "ySupRunoff5"          =ySupRunoff5,
  "ySupOutflow1"         =ySupOut1,
  "ySupOutflow2"         =ySupOut2,
#  "ySupOutflow3"         =ySupOut3,
  "yMHRStore"           =yMHRStore,
  "PrecipShapeMH"       =PrecipShapeMH,
  "PrecipRateMH"        =PrecipRateMH,
  "EvapMeanMH"          =EvapMeanMH,
  "EvapPrecisionMH"     =EvapPrecisionMH,
  "LogRunoffMeanMH"     =LogRunoffMeanMH,
  "LogRunoffPrecisionMH"=LogRunoffPrecisionMH,
  "OutflowMeanMH"       =OutflowMeanMH,
  "OutflowPrecisionMH"  =OutflowPrecisionMH,
  "yMHPrecip1"          =yMHPrecip1,
  "yMHPrecip2"          =yMHPrecip2,
  "yMHPrecip3"          =yMHPrecip3,
  "yMHPrecip4"          =yMHPrecip4,
  "yMHPrecip5"          =yMHPrecip5,
  "yMHPrecip6"          =yMHPrecip6,
  "yMHPrecip7"          =yMHPrecip7,
  "yMHPrecip8"          =yMHPrecip8,
  "yMHPrecip9"          =yMHPrecip9,
  "yMHEvap1"            =yMHEvap1,
#  "yMHEvap2"            =yMHEvap2,
  "yMHEvap3"            =yMHEvap3,
  "yMHEvap4"            =yMHEvap4,
  "yMHEvap5"            =yMHEvap5,
  "yMHRunoff1"          =yMHRunoff1,
#  "yMHRunoff2"          =yMHRunoff2,
  "yMHRunoff3"          =yMHRunoff3,
  "yMHRunoff4"          =yMHRunoff4,
  "yMHRunoff5"          =yMHRunoff5,
  "yMHOutflow1"         =yMHOut1,
  "yMHOutflow2"         =yMHOut2,
#  "yMHOutflow3"         =yMHOut3,
  "yLkStClRStore"         =yLkStClRStore,
  "NBSMeanLkStCl"         =NBSMeanLkStCl,
  "NBSPrecisionLkStCl"    =NBSPrecisionLkStCl,
  "OutflowMeanLkStCl"     =OutflowMeanLkStCl,
  "OutflowPrecisionLkStCl"=OutflowPrecisionLkStCl,
  "yLkStClNBS1"           =yLkStClNBS1,
#  "yLkStClNBS2"           =yLkStClNBS2,
  "yLkStClNBS3"           =yLkStClNBS3,
  "yLkStClNBS4"           =yLkStClNBS4,
  "yLkStClNBS5"           =yLkStClNBS5,
  "yLkStClNBS6"           =yLkStClNBS6,
#sk     "yLkStClOutflow1"       =yLkStClOut1,
#sk     "yLkStClOutflow2"       =yLkStClOut2,
#  "yLkStClOutflow3"       =yLkStClOut3,
  "yEriRStore"           =yEriRStore,
  "PrecipShapeEri"       =PrecipShapeEri,
  "PrecipRateEri"        =PrecipRateEri,
  "EvapMeanEri"          =EvapMeanEri,
  "EvapPrecisionEri"     =EvapPrecisionEri,
  "LogRunoffMeanEri"     =LogRunoffMeanEri,
  "LogRunoffPrecisionEri"=LogRunoffPrecisionEri,
  "OutflowMeanEri"       =OutflowMeanEri,
  "OutflowPrecisionEri"  =OutflowPrecisionEri,
  "yEriPrecip1"          =yEriPrecip1,
#  "yEriPrecip2"          =yEriPrecip2,
  "yEriPrecip3"          =yEriPrecip3,
  "yEriPrecip4"          =yEriPrecip4,
  "yEriPrecip5"          =yEriPrecip5,
  "yEriPrecip6"          =yEriPrecip6,
  "yEriPrecip7"          =yEriPrecip7,
  "yEriPrecip8"          =yEriPrecip8,
  "yEriPrecip9"          =yEriPrecip9,
  "yEriEvap1"            =yEriEvap1,
#  "yEriEvap2"            =yEriEvap2,
  "yEriEvap3"            =yEriEvap3,
  "yEriEvap4"            =yEriEvap4,
  "yEriEvap5"            =yEriEvap5,
  "yEriRunoff1"          =yEriRunoff1,
#  "yEriRunoff2"          =yEriRunoff2,
  "yEriRunoff3"          =yEriRunoff3,
  "yEriRunoff4"          =yEriRunoff4,
  "yEriRunoff5"          =yEriRunoff5,
  "yEriOutflow1"         =yEriOut1
#  "yEriOutflow2"         =yEriOut2,
#  "yEriOutflow3"         =yEriOut3
  )

##### Variables that are of interest to us #####
ParamsToMonitor=c(
  "SupPrecip",
  "SupEvap",
  "SupRunoff",
  "SupOutflow",
  "SupOutflow_mm",
  "SupDStore",
  "SupRStore",
  "SupError",
  "ySupPrecip1Bias",
#  "ySupPrecip2Bias",
  "ySupPrecip3Bias",
  "ySupPrecip4Bias",
  "ySupPrecip5Bias",
  "ySupPrecip6Bias",
#  "ySupPrecip7Bias",
  "ySupPrecip8Bias",
  "ySupPrecip9Bias",
  "ySupEvap1Bias",
#  "ySupEvap2Bias",
  "ySupEvap3Bias",
  "ySupEvap4Bias",
  "ySupEvap5Bias",
  "ySupRunoff1Bias",
#  "ySupRunoff2Bias",
  "ySupRunoff3Bias",
  "ySupRunoff4Bias",
  "ySupRunoff5Bias",
  "ySupOutflow1Bias",
  "ySupOutflow2Bias",
#  "ySupOutflow3Bias",
  "ySupRStorePreci",
  "ySupPrecip1Precision",
#  "ySupPrecip2Precision",
  "ySupPrecip3Precision",
  "ySupPrecip4Precision",
  "ySupPrecip5Precision",
  "ySupPrecip6Precision",
#  "ySupPrecip7Precision",
  "ySupPrecip8Precision",
  "ySupPrecip9Precision",
  "ySupEvap1Precision",
#  "ySupEvap2Precision",
  "ySupEvap3Precision",
  "ySupEvap4Precision",
  "ySupEvap5Precision",
  "ySupRunoff1Precision",
#  "ySupRunoff2Precision",
  "ySupRunoff3Precision",
  "ySupRunoff4Precision",
  "ySupRunoff5Precision",
  "ySupOutflow1Precision",
  "ySupOutflow2Precision",
#  "ySupOutflow3Precision",
  "ySupPrecip1PP",
#  "ySupPrecip2PP",
  "ySupPrecip3PP",
  "ySupPrecip4PP",
  "ySupPrecip5PP",
  "ySupPrecip6PP",
#  "ySupPrecip7PP",
  "ySupPrecip8PP",
  "ySupPrecip9PP",
  "ySupEvap1PP",
#  "ySupEvap2PP",
  "ySupEvap3PP",
  "ySupEvap4PP",
  "ySupEvap5PP",
  "ySupRunoff1PP",
#  "ySupRunoff2PP",
  "ySupRunoff3PP",
  "ySupRunoff4PP",
  "ySupRunoff5PP",
  "ySupOutflow1PP",
  "ySupOutflow2PP",
#  "ySupOutflow3PP",
  "ySupDStorePP", ## End Superior monitored parameters ##
  "MHPrecip",
  "MHEvap",
  "MHRunoff",
  "MHOutflow",
  "MHOutflow_mm",
  "MHDStore",
  "MHRStore",
  "MHError", ## When incorporating SFD: probably want to add in SFDc, SFDd, SFDf to track too ##
  "yMHPrecip1Bias",
  "yMHPrecip2Bias",
  "yMHPrecip3Bias",
  "yMHPrecip4Bias",
  "yMHPrecip5Bias",
  "yMHPrecip6Bias",
  "yMHPrecip7Bias",
  "yMHPrecip8Bias",
  "yMHPrecip9Bias",
  "yMHEvap1Bias",
#  "yMHEvap2Bias",
  "yMHEvap3Bias",
  "yMHEvap4Bias",
  "yMHEvap5Bias",
  "yMHRunoff1Bias",
#  "yMHRunoff2Bias",
  "yMHRunoff3Bias",
  "yMHRunoff4Bias",
  "yMHRunoff5Bias",
  "yMHOutflow1Bias",
  "yMHOutflow2Bias",
#  "yMHOutflow3Bias",
  "yMHRStorePreci",
  "yMHPrecip1Precision",
  "yMHPrecip2Precision",
  "yMHPrecip3Precision",
  "yMHPrecip4Precision",
  "yMHPrecip5Precision",
  "yMHPrecip6Precision",
  "yMHPrecip7Precision",
  "yMHPrecip8Precision",
  "yMHPrecip9Precision",
  "yMHEvap1Precision",
#  "yMHEvap2Precision",
  "yMHEvap3Precision",
  "yMHEvap4Precision",
  "yMHEvap5Precision",
  "yMHRunoff1Precision",
#  "yMHRunoff2Precision",
  "yMHRunoff3Precision",
  "yMHRunoff4Precision",
  "yMHRunoff5Precision",
  "yMHOutflow1Precision",
  "yMHOutflow2Precision",
#  "yMHOutflow3Precision",
  "yMHPrecip1PP",
  "yMHPrecip2PP",
  "yMHPrecip3PP",
  "yMHPrecip4PP",
  "yMHPrecip5PP",
  "yMHPrecip6PP",
  "yMHPrecip7PP",
  "yMHPrecip8PP",
  "yMHPrecip9PP",
  "yMHEvap1PP",
#  "yMHEvap2PP",
  "yMHEvap3PP",
  "yMHEvap4PP",
  "yMHEvap5PP",
  "yMHRunoff1PP",
#  "yMHRunoff2PP",
  "yMHRunoff3PP",
  "yMHRunoff4PP",
  "yMHRunoff5PP",
"yMHOutflow1PP",
   "yMHOutflow2PP",
#  "yMHOutflow3PP",
  "yMHDStorePP", ## H.D. ADDED ON 12 FEB 2021 ## End Michigan-Huron monitored parameters ##
  "LkStClNBS",
  "LkStClOutflow",
  "LkStClOutflow_mm",
  "LkStClDStore",
  "LkStClRStore",
  "LkStClError",
  "yLkStClNBS1Bias",
#  "yLkStClNBS2Bias",
  "yLkStClNBS3Bias",
  "yLkStClNBS4Bias",
  "yLkStClNBS5Bias",
  "yLkStClNBS6Bias",
#sk   "yLkStClOutflow1Bias",
#sk   "yLkStClOutflow2Bias",
#  "yLkStClOutflow3Bias",
  "yLkStClRStorePreci",
  "yLkStClNBS1Precision",
#  "yLkStClNBS2Precision",
  "yLkStClNBS3Precision",
  "yLkStClNBS4Precision",
  "yLkStClNBS5Precision",
  "yLkStClNBS6Precision",
#sk   "yLkStClOutflow1Precision",
#sk  "yLkStClOutflow2Precision",
#  "yLkStClOutflow3Precision",
  "yLkStClNBS1PP",
#  "yLkStClNBS2PP",
  "yLkStClNBS3PP",
  "yLkStClNBS4PP",
  "yLkStClNBS5PP",
  "yLkStClNBS6PP",
#sk   "yLkStClOutflow1PP",
#sk  "yLkStClOutflow2PP",
#  "yLkStClOutflow3PP",
  "yLkStClDStorePP", ## End LAKE St. Clair monitored parameters ##
  #"Mult1Bias",
  #"Mult2Bias",
  "EriPrecip",
  "EriEvap",
  "EriRunoff",
  "EriOutflow",
  "EriOutflow_mm",
  "EriDStore",
  "EriRStore",
  "EriError",
  "yEriOutflow1",#--JLC
#  "yEriOutflow2",#--JLC
#  "yEriOutflow3",#--JLC
  "yEriOutflow1Mean",#--JLC
#  "yEriOutflow2Mean",#--JLC
#  "yEriOutflow3Mean",#--JLC
  "yEriPrecip1Bias",
#  "yEriPrecip2Bias",
  "yEriPrecip3Bias",
  "yEriPrecip4Bias",
  "yEriPrecip5Bias",
  "yEriPrecip6Bias",
  "yEriPrecip7Bias",
  "yEriPrecip8Bias",
  "yEriPrecip9Bias",
  "yEriEvap1Bias",
#  "yEriEvap2Bias",
  "yEriEvap3Bias",
  "yEriEvap4Bias",
  "yEriEvap5Bias",
  "yEriRunoff1Bias",
#  "yEriRunoff2Bias",
  "yEriRunoff3Bias",
  "yEriRunoff4Bias",
  "yEriRunoff5Bias",
  "yEriOutflow1Bias",
#  "yEriOutflow2Bias",
#  "yEriOutflow3Bias",
  "yEriRStorePreci",
  "yEriPrecip1Precision",
#  "yEriPrecip2Precision",
  "yEriPrecip3Precision",
  "yEriPrecip4Precision",
  "yEriPrecip5Precision",
  "yEriPrecip6Precision",
  "yEriPrecip7Precision",
  "yEriPrecip8Precision",
  "yEriPrecip9Precision",
  "yEriEvap1Precision",
#  "yEriEvap2Precision",
  "yEriEvap3Precision",
  "yEriEvap4Precision",
  "yEriEvap5Precision",
  "yEriRunoff1Precision",
#  "yEriRunoff2Precision",
  "yEriRunoff3Precision",
  "yEriRunoff4Precision",
  "yEriRunoff5Precision",
  "yEriOutflow1Precision",
#  "yEriOutflow2Precision",
#  "yEriOutflow3Precision",
  "yEriPrecip1PP",
#  "yEriPrecip2PP",
  "yEriPrecip3PP",
  "yEriPrecip4PP",
  "yEriPrecip5PP",
  "yEriPrecip6PP",
  "yEriPrecip7PP",
  "yEriPrecip8PP",
  "yEriPrecip9PP",
  "yEriEvap1PP",
#  "yEriEvap2PP",
  "yEriEvap3PP",
  "yEriEvap4PP",
  "yEriEvap5PP",
  "yEriRunoff1PP",
#  "yEriRunoff2PP",
  "yEriRunoff3PP",
  "yEriRunoff4PP",
  "yEriRunoff5PP",
  "yEriOutflow1PP",
#  "yEriOutflow2PP",
#  "yEriOutflow3PP",
  "yEriDStorePP" ## End Erie monitored parameters ##
)

##### Call to JAGS #####
##### Console output #####
cat(date(),paste(RollPeriod,' ROLL, ERA5 DATA, 3CHAINS, ',halfIters*2/1000,'K ITERATIONS\n\n', sep=''))
##### Progress output #####
cat(date(),'ADAPTING SAMPLER TO MODEL...')
##### Initializes and adapts the model #####
jMod = jags.model(file = "SimplifiedL2S_jagsModel.jags.r", data = inputDataCoreJAGS, n.chains = 2) 

# Progress output for burning
cat(date(),paste0('UPDATE STEP (',halfIters/1000,'K BURNIN)...\n'))
update(jMod, halfIters) 

# Sampling
cat(date(),'SAMPLING... (',halfIters/1000,'K, THINNING DOWN TO 3K)\n')
jSample = coda.samples(jMod, ParamsToMonitor, halfIters, ceiling(halfIters/1000), na.rm=TRUE)

cat(date(),'COMPUTING STATS...\n')
jSumStats     = summary(jSample)
jSumStats_MSD = jSumStats$statistics[,1:2]
jSumStats_Q   = jSumStats$quantiles
jSumEff       = effectiveSize(jSample)

print(jSumStats$quantiles[grepl(pattern="LkStClError",rownames(jSumStats$statistics)),])
print(jSumStats$quantiles[grepl(pattern="yLkStClRStorePreci",rownames(jSumStats$statistics)),]^-0.5)
print(jSumStats$quantiles[grepl(pattern="StClNBS[0-9]Bias",rownames(jSumStats$statistics)),])
print(jSumStats$quantiles[grepl(pattern="yLkStClNBS[0-9]Precision",rownames(jSumStats$statistics)),]^-0.5)




# # Get R-Hats - NOT USED!
# cat('GETTING GELMAN-RUBIN STAT (COMPUTATIONALLY EXPENSIVE)...\n')
# jRHat = gelman.diag(jSample, multivariate=FALSE)
# jRHatEsts = jRHat$psrf

# Print the output
rn     = rownames(jSumStats_MSD)
ro_MSD = match(rn, rownames(jSumStats_MSD))
ro_Q   = match(rn, rownames(jSumStats_Q))
ro_eff = match(rn, names(jSumEff))
# ro_rhat = match(rn, rownames(jRHatEsts))

jSum = cbind(
  jSumStats_MSD,
  jSumStats_Q[ro_Q,],
  # jRHatEsts[ro_rhat,],
  jSumEff[ro_eff]
)

colnames(jSum)[8] = 'n.eff'

print(jSum[grepl(pattern="LkStClError",rownames(jSum)),])
print(jSum[grepl(pattern="yLkStClRStorePreci",rownames(jSum)),]^-0.5)
print(jSum[grepl(pattern="StClNBS[0-9]Bias",rownames(jSum)),])
print(jSum[grepl(pattern="yLkStClOutflow[0-9]Precision",rownames(jSum)),]^-0.5)
print(jSum[grepl(pattern="yEriOutflow[0-9]Precision",rownames(jSum)),]^-0.5)

plot(jSample[[1]][,"yLkStClRStorePreci"]^-0.5)
#plot(jSample[[1]][,"yLkStClNBS1Precision"]^-0.5)
plot(jSample[[1]][,"yLkStClNBS2Precision"]^-0.5)
plot(jSample[[1]][,"yLkStClNBS3Precision"]^-0.5)
plot(jSample[[1]][,"yLkStClNBS4Precision"]^-0.5)
plot(jSample[[1]][,"yLkStClNBS5Precision"]^-0.5)
plot(jSample[[1]][,"yLkStClNBS6Precision"]^-0.5)

plot(jSample[[1]][,"ySupOutflow1Precision"]^-0.5)
plot(jSample[[1]][,"ySupOutflow2Precision"]^-0.5)
plot(jSample[[1]][,"yMHOutflow1Precision"]^-0.5)
plot(jSample[[1]][,"yMHOutflow2Precision"]^-0.5)
#plot(jSample[[1]][,"yMHOutflow3Precision"]^-0.5)
plot(jSample[[1]][,"yLkStClOutflow1Precision"]^-0.5)
plot(jSample[[1]][,"yLkStClOutflow2Precision"]^-0.5)
plot(jSample[[1]][,"yEriRStorePreci"]^-0.5)
plot(jSample[[1]][,"yEriPrecip1Precision"]^-0.5)
#plot(jSample[[1]][,"yEriPrecip2Precision"]^-0.5)
plot(jSample[[1]][,"yEriPrecip3Precision"]^-0.5)
plot(jSample[[1]][,"yEriPrecip4Precision"]^-0.5)
plot(jSample[[1]][,"yEriPrecip5Precision"]^-0.5)
plot(jSample[[1]][,"yEriPrecip6Precision"]^-0.5)
plot(jSample[[1]][,"yEriPrecip7Precision"]^-0.5)
plot(jSample[[1]][,"yEriPrecip8Precision"]^-0.5)
plot(jSample[[1]][,"yEriPrecip9Precision"]^-0.5)
plot(jSample[[1]][,"yEriEvap1Precision"]^-0.5)
#plot(jSample[[1]][,"yEriEvap2Precision"]^-0.5)
plot(jSample[[1]][,"yEriEvap3Precision"]^-0.5)
plot(jSample[[1]][,"yEriEvap4Precision"]^-0.5)
plot(jSample[[1]][,"yEriEvap5Precision"]^-0.5)
plot(jSample[[1]][,"yEriRunoff1Precision"]^-0.5)
#plot(jSample[[1]][,"yEriRunoff2Precision"]^-0.5)
plot(jSample[[1]][,"yEriRunoff3Precision"]^-0.5)
plot(jSample[[1]][,"yEriRunoff4Precision"]^-0.5)
plot(jSample[[1]][,"yEriRunoff5Precision"]^-0.5)
plot(jSample[[1]][,"yEriOutflow1Precision"]^-0.5)
#plot(jSample[[1]][,"yEriOutflow2Precision"]^-0.5)


plot(jSample[[1]][,"yEriOutflow1Mean[1]"])
p#lot(jSample[[1]][,"yEriOutflow2Mean[1]"])

print(jSum[grepl(pattern="EriOutflow_mm",rownames(jSum)),])


cat(date(),'SAVING WORK...\n')
workspacename = paste0(outputname,startAnalysisYear,endAnalysisYear,'_',(halfIters*2)/1000,'k_roll',RollPeriod,'.RData')
save.image(workspacename)
cat(date(),'DONE!\n')

cat(date(),'CREATING STATS CSV...\n')
Plotting.on.iteration.time=3000
fileNameBase=paste(outputname,Plotting.on.iteration.time/1000,'k',sep='')
write.table(jSum,paste(fileNameBase,'_stats.csv',sep=''),sep=',',quote=FALSE )
cat(date(),'DONE!\n')

cat(date(),'CREATING PLOTS BY 5 YEAR INTERVAL...\n')
source('SimplifiedL2S_plotter_5y.R')
cat(date(),'DONE!\n')
