##############################################################################################################################
##############################################################################################################################
##R CODE FOR AN EXPERT-BASED STOCHASTIC POPULATION FORECAST, USING AUTOREGRESSIVE MODELS WITH RANDOM COEFFICIENTS, APPLIED TO ALASKA
##
##EDDIE HUNSINGER, SEPTEMBER 2010 (UPDATED JUNE 2020)
##http://www.demog.berkeley.edu/~eddieh/
##edyhsgr@gmail.com
##
##TO SEE THE RELATED PAPER (WITH REFERENCES), GO TO http://www.demog.berkeley.edu/~eddieh/documents/ExpertForecastPaper.pdf
##
##IF YOU WOULD LIKE TO USE, SHARE OR REPRODUCE ANY INFORMATION OR IDEAS FROM THIS WORK, BE SURE TO CITE THE SOURCE
##
##MUCH OF THIS PROJECTION CODE IS BASED ON CODE THAT I DEVELOPED WHILE WORKING FOR THE ALASKA DEPARTMENT OF LABOR AND WORKFORCE DEVELOPMENT
##I ALSO BENEFITED FROM FACULTY, STAFF AND FELLOW STUDENTS AT THE UC BERKELEY DEPARTMENT OF DEMOGRAPHY
##THERE IS NO WARRANTY FOR THIS CODE
##############################################################################################################################
##############################################################################################################################

##############################################################################################################################
##############################################################################################################################
##SELECT THE INPUT DATA
##############################################################################################################################
##############################################################################################################################

##DIMENSIONS
##
##NUMBER OF ITERATIONS
##SET TO 10,000 
ITER<-10000
##SIZE OF PROJECTION MATRIX
SIZE<-21
##NUMBER OF PROJECTION STEPS
STEPS<-6
BASEANDSTEPS<-STEPS+1

##SURVIVAL PARAMETERS
##
##YEAR 2000 US LIFE TABLE lx CURVES FROM THE NATIONAL CENTER FOR HEALTH STATISTICS
Survival<-read.table(file="https://github.com/AppliedDemogToolbox/Hunsinger_StochasticForecast/raw/master/lx2000_US_NCHS.csv",header=TRUE,sep=",")
lxF<-Survival$F_2000
lxM<-Survival$M_2000
##"BA" IS THE BRASS RELATIONAL LOGIT MODEL ALPHA. IT IS CALIBRATED FOR THE JUMP-OFF PERIOD AND TO FOLLOW A FIXED RANDOM WALK WITH DRIFT MODEL,
##WITH THE DRIFT TERM SET TO .03 FOR EACH FIVE-YEAR STEP, AND THE RANDOM TERM WITH MEAN 0 AND STANDARD ERROR .0182
BA_baseF<-.03
BA_baseM<-.08
BA_consF<-.03
BA_consM<-.03
BA_e<-.0182

##FERTILITY PARAMETERS
##
##YEAR 2005 US FERTILITY RATES FROM THE NATIONAL CENTER FOR HEALTH STATISTICS, SUMMED TO 1
Fertility<-read.table(file="https://github.com/AppliedDemogToolbox/Hunsinger_StochasticForecast/raw/master/Fx2005_US_NCHS.csv",header=TRUE,sep=",")
PropFx<-c(Fertility$PropFx)
##FRACTION FEMALE AT BIRTH
ffab<-.4886
##"TFR" IS THE TOTAL FERTILITY RATE. IT IS SET TO AVERAGE 2.3 FOR 2/3 OF THE FORECAST PERIODS, 2.1 FOR 1/3 OF THE FORECAST PERIODS, AND TO FOLLOW A FIRST-ORDER AUTOREGRESSIVE MODEL,
##WITH UNIFORM RANDOM VALUES BETWEEN .5 AND 1 USED FOR THE AUTOREGRESSIVE TERM, AND THE RANDOM TERM WITH MEAN 0 AND STANDARD ERROR .243
TFR_base<-2.37
TFR_cons<-c(2.3,2.1)
TFR_consprob<-c(2/3,1/3)
TFR_ar<-array(runif(ITER,.5,1))
TFR_e<-.243

##NETMIGRATION PARAMETERS
##
##1995-2000 MIGRATION PROFILES FROM THE 2000 CENSUS, SUMMED TO 1
Migration<-read.table(file="https://github.com/AppliedDemogToolbox/Hunsinger_StochasticForecast/raw/master/MigProf2000_AKCACO_Census2000.csv",header=TRUE,sep=",")
PropInM<-c(Migration$AKIn)
PropInF<-c(Migration$AKIn)
PropOutM<-c(Migration$AKOut)
PropOutF<-c(Migration$AKOut)
##THE OUT MIGRATION RATIO, WHICH IS MULTIPLIED BY THE OUT MIGRATION PROFILE AND POPULATION FOR EACH SEX
OutRate<-.2
##"TIR" IS THE NET MIGRATION RATIO. IT IS SET TO AVERAGE 0 FOR ALL FORECAST PERIODS, AND TO FOLLOW A FIRST-ORDER AUTOREGRESSIVE MODEL, 
##WITH UNIFORM RANDOM VALUES BETWEEN -.5 AND 1 USED FOR THE AUTOREGRESSIVE TERM, AND THE RANDOM TERM WITH MEAN 0 AND STANDARD ERROR .0303
TIR_base<-.0025
TIR_cons<-.00
TIR_ar<-array(runif(ITER,-.5,1))
TIR_e<-.0243
##PIPELINE EVENT. PROBABILITY OF SUCH AN EVENT IN A SINGLE FORECAST PERIOD IS SET TO .1 FOR POSITIVE, AND .05 FOR NEGATIVE; LEVEL IS .1 POSITIVE OR .1 NEGATIVE.
PipeProb<-c(.05,.20)
PipeLev<-c(-.1,.1,0)

##BASE POPULATION
##
##VINTAGE 2009 POPULATION ESTIMATES FOR 2005, FROM THE ALASKA DEPARTMENT OF LABOR AND WORKFORCE DEVELOPMENT (http://almis.labor.state.ak.us/?PAGEID=67&SUBID=115)
K05<-read.table(file="https://github.com/AppliedDemogToolbox/Hunsinger_StochasticForecast/raw/master/AgeSex2005_AK_AKDOLv2009.csv",header=TRUE,sep=",")
KF05<-K05$AK_F_2005
KM05<-K05$AK_M_2005

##############################################################################################################################
##############################################################################################################################
##RUN THE FORECAST CODE
##############################################################################################################################
##############################################################################################################################

##SURVIVAL
BA_error<-array(0,c(BASEANDSTEPS,ITER))
for(i in 2:BASEANDSTEPS){BA_error[i,]<-rnorm(ITER,0,BA_e)}
BAF<-array(0,c(BASEANDSTEPS,ITER))
BAF[1,]<-BA_baseF
for(i in 2:BASEANDSTEPS){BAF[i,]<-BAF[i-1,]+BA_consF+BA_error[i,]}
BAF<-t(BAF)
BrassF00<-data.frame(Alpha=BAF[,1],Beta=1)
BrassF05<-data.frame(Alpha=BAF[,2],Beta=1)
BrassF10<-data.frame(Alpha=BAF[,3],Beta=1)
BrassF15<-data.frame(Alpha=BAF[,4],Beta=1)
BrassF20<-data.frame(Alpha=BAF[,5],Beta=1)
BrassF25<-data.frame(Alpha=BAF[,6],Beta=1)
BrassF30<-data.frame(Alpha=BAF[,7],Beta=1)
BAM<-array(0,c(BASEANDSTEPS,ITER))
BAM[1,]<-BA_baseM
for(i in 2:BASEANDSTEPS){BAM[i,]<-BAM[i-1,]+BA_consM+BA_error[i,]}
BAM<-t(BAM)
BrassM00<-data.frame(Alpha=BAM[,1],Beta=1)
BrassM05<-data.frame(Alpha=BAM[,2],Beta=1)
BrassM10<-data.frame(Alpha=BAM[,3],Beta=1)
BrassM15<-data.frame(Alpha=BAM[,4],Beta=1)
BrassM20<-data.frame(Alpha=BAM[,5],Beta=1)
BrassM25<-data.frame(Alpha=BAM[,6],Beta=1)
BrassM30<-data.frame(Alpha=BAM[,7],Beta=1)

##FERTILITY
fmab <- 1-ffab
Fx<-array(0,c(SIZE,ITER))
Fx[1:SIZE,]<-PropFx
TFRFORE<-array(0,c(BASEANDSTEPS,ITER))
TFRFORE[1,]<-TFR_base
TFRCONSPROB<-sample(TFR_cons,ITER,replace=TRUE,prob=TFR_consprob)
for(i in 2:BASEANDSTEPS){TFRFORE[i,]<-TFR_ar[]*TFRFORE[i-1,]+TFRCONSPROB[]*(1-TFR_ar[])+rnorm(ITER,0,TFR_e)}
TFRFORE<-t(TFRFORE)
TFRFORE[TFRFORE<0]<-0

##NETMIGRATION
TIRFORE<-array(0,c(BASEANDSTEPS,ITER))
TIRFORE[1,]<-TIR_base
for(i in 2:BASEANDSTEPS){TIRFORE[i,]<-TIR_ar[]*TIRFORE[i-1,]+TIR_cons*(1-TIR_ar[])+rnorm(ITER,0,TIR_e)}
TIRFORE<-t(TIRFORE)
PipeLevProb<-array(0,c(BASEANDSTEPS,ITER))
for(i in 2:BASEANDSTEPS){PipeLevProb[i,]<-sample(PipeLev,ITER,replace=TRUE,prob=c(PipeProb[1]/STEPS,PipeProb[2]/STEPS,1-(PipeProb[1]/STEPS+PipeProb[2]/STEPS)))}
PipeLevProb<-t(PipeLevProb)
TIRFORE<-TIRFORE+PipeLevProb

##CALCULATE THE Yx FOR THE lx'S
YxM<-YxF<-NULL
for (i in 1:length(lxF)){YxF[i]<-.5*log(lxF[i]/(1-lxF[i]))}
for (i in 1:length(lxM)){YxM[i]<-.5*log(lxM[i]/(1-lxM[i]))}

##Improve Survival and make lx's for each period
lxF30<-lxF25<-lxF20<-lxF15<-lxF10<-lxF05<-lxF00<-array(0,c(SIZE+1,ITER))
lxM30<-lxM25<-lxM20<-lxM15<-lxM10<-lxM05<-lxM00<-array(0,c(SIZE+1,ITER))
for (i in 1:length(lxF)){lxF00[i,]<-1/(1+exp(-2*BrassF00$Alpha-2*BrassF00$Beta*YxF[i]))}
for (i in 1:length(lxM)){lxM00[i,]<-1/(1+exp(-2*BrassM00$Alpha-2*BrassM00$Beta*YxM[i]))}
for (i in 1:length(lxF)){lxF05[i,]<-1/(1+exp(-2*BrassF05$Alpha-2*BrassF05$Beta*YxF[i]))}
for (i in 1:length(lxM)){lxM05[i,]<-1/(1+exp(-2*BrassM05$Alpha-2*BrassM05$Beta*YxM[i]))}
for (i in 1:length(lxF)){lxF10[i,]<-1/(1+exp(-2*BrassF10$Alpha-2*BrassF10$Beta*YxF[i]))}
for (i in 1:length(lxM)){lxM10[i,]<-1/(1+exp(-2*BrassM10$Alpha-2*BrassM10$Beta*YxM[i]))}
for (i in 1:length(lxF)){lxF15[i,]<-1/(1+exp(-2*BrassF15$Alpha-2*BrassF15$Beta*YxF[i]))}
for (i in 1:length(lxM)){lxM15[i,]<-1/(1+exp(-2*BrassM15$Alpha-2*BrassM15$Beta*YxM[i]))}
for (i in 1:length(lxF)){lxF20[i,]<-1/(1+exp(-2*BrassF20$Alpha-2*BrassF20$Beta*YxF[i]))}
for (i in 1:length(lxM)){lxM20[i,]<-1/(1+exp(-2*BrassM20$Alpha-2*BrassM20$Beta*YxM[i]))}
for (i in 1:length(lxF)){lxF25[i,]<-1/(1+exp(-2*BrassF25$Alpha-2*BrassF25$Beta*YxF[i]))}
for (i in 1:length(lxM)){lxM25[i,]<-1/(1+exp(-2*BrassM25$Alpha-2*BrassM25$Beta*YxM[i]))}
for (i in 1:length(lxF)){lxF30[i,]<-1/(1+exp(-2*BrassF30$Alpha-2*BrassF30$Beta*YxF[i]))}
for (i in 1:length(lxM)){lxM30[i,]<-1/(1+exp(-2*BrassM30$Alpha-2*BrassM30$Beta*YxM[i]))}

##MAKE nLx's FOR EACH PERIOD
LxF30<-LxF25<-LxF20<-LxF25<-LxF15<-LxF10<-LxF05<-array(0,c(SIZE,ITER))
LxM30<-LxM25<-LxM20<-LxM15<-LxM15<-LxM10<-LxM05<-array(0,c(SIZE,ITER))
##**THIS IS A LITTLE OFF FOR THE FIRST AGE GROUP**
for (i in 1:SIZE){LxF05[i,]<-2.5*(lxF05[i,]+lxF05[i+1,])}
for (i in 1:SIZE){LxM05[i,]<-2.5*(lxM05[i,]+lxM05[i+1,])}
for (i in 1:SIZE){LxF10[i,]<-2.5*(lxF10[i,]+lxF10[i+1,])}
for (i in 1:SIZE){LxM10[i,]<-2.5*(lxM10[i,]+lxM10[i+1,])}
for (i in 1:SIZE){LxF15[i,]<-2.5*(lxF15[i,]+lxF15[i+1,])}
for (i in 1:SIZE){LxM15[i,]<-2.5*(lxM15[i,]+lxM15[i+1,])}
for (i in 1:SIZE){LxF20[i,]<-2.5*(lxF20[i,]+lxF20[i+1,])}
for (i in 1:SIZE){LxM20[i,]<-2.5*(lxM20[i,]+lxM20[i+1,])}
for (i in 1:SIZE){LxF25[i,]<-2.5*(lxF25[i,]+lxF25[i+1,])}
for (i in 1:SIZE){LxM25[i,]<-2.5*(lxM25[i,]+lxM25[i+1,])}
for (i in 1:SIZE){LxF30[i,]<-2.5*(lxF30[i,]+lxF30[i+1,])}
for (i in 1:SIZE){LxM30[i,]<-2.5*(lxM30[i,]+lxM30[i+1,])}

##TABLE e0
e0MFORE<-array(0,c(BASEANDSTEPS,ITER))
for (i in 1:ITER){e0MFORE[2,i]<-sum(LxM05[,i])}
for (i in 1:ITER){e0MFORE[3,i]<-sum(LxM10[,i])}
for (i in 1:ITER){e0MFORE[4,i]<-sum(LxM15[,i])}
for (i in 1:ITER){e0MFORE[5,i]<-sum(LxM20[,i])}
for (i in 1:ITER){e0MFORE[6,i]<-sum(LxM25[,i])}
for (i in 1:ITER){e0MFORE[7,i]<-sum(LxM30[,i])}
e0FFORE<-array(0,c(BASEANDSTEPS,ITER))
for (i in 1:ITER){e0FFORE[2,i]<-sum(LxF05[,i])}
for (i in 1:ITER){e0FFORE[3,i]<-sum(LxF10[,i])}
for (i in 1:ITER){e0FFORE[4,i]<-sum(LxF15[,i])}
for (i in 1:ITER){e0FFORE[5,i]<-sum(LxF20[,i])}
for (i in 1:ITER){e0FFORE[6,i]<-sum(LxF25[,i])}
for (i in 1:ITER){e0FFORE[7,i]<-sum(LxF30[,i])}

##MAKE nSx's FOR EACH PERIOD
SxF30<-SxF25<-SxF20<-SxF25<-SxF15<-SxF10<-SxF05<-array(0,c(SIZE-1,ITER))
SxM30<-SxM25<-SxM20<-SxM15<-SxM15<-SxM10<-SxM05<-array(0,c(SIZE-1,ITER))
for (i in 1:SIZE-1){SxF05[i,]<-(LxF05[i+1,]/LxF05[i,])}
for (i in 1:SIZE-1){SxM05[i,]<-(LxM05[i+1,]/LxM05[i,])}
for (i in 1:SIZE-1){SxF10[i,]<-(LxF10[i+1,]/LxF10[i,])}
for (i in 1:SIZE-1){SxM10[i,]<-(LxM10[i+1,]/LxM10[i,])}
for (i in 1:SIZE-1){SxF15[i,]<-(LxF15[i+1,]/LxF15[i,])}
for (i in 1:SIZE-1){SxM15[i,]<-(LxM15[i+1,]/LxM15[i,])}
for (i in 1:SIZE-1){SxF20[i,]<-(LxF20[i+1,]/LxF20[i,])}
for (i in 1:SIZE-1){SxM20[i,]<-(LxM20[i+1,]/LxM20[i,])}
for (i in 1:SIZE-1){SxF25[i,]<-(LxF25[i+1,]/LxF25[i,])}
for (i in 1:SIZE-1){SxM25[i,]<-(LxM25[i+1,]/LxM25[i,])}
for (i in 1:SIZE-1){SxF30[i,]<-(LxF30[i+1,]/LxF30[i,])}
for (i in 1:SIZE-1){SxM30[i,]<-(LxM30[i+1,]/LxM30[i,])}

##PUT THE Sx DATA INTO THE SUBDIAGONAL OF WHAT WILL BE THE LESLIE MATRICES
SF30<-SF25<-SF20<-SF15<-SF10<-SF05<-array(0,c(SIZE,SIZE,ITER))
SM30<-SM25<-SM20<-SM15<-SM10<-SM05<-array(0,c(SIZE,SIZE,ITER))
for (i in 1:ITER){SF05[,,i]<-rbind(0,cbind(diag(SxF05[,i]),0))}
for (i in 1:ITER){SF10[,,i]<-rbind(0,cbind(diag(SxF10[,i]),0))}
for (i in 1:ITER){SF15[,,i]<-rbind(0,cbind(diag(SxF15[,i]),0))}
for (i in 1:ITER){SF20[,,i]<-rbind(0,cbind(diag(SxF20[,i]),0))}
for (i in 1:ITER){SF25[,,i]<-rbind(0,cbind(diag(SxF25[,i]),0))}
for (i in 1:ITER){SF30[,,i]<-rbind(0,cbind(diag(SxF30[,i]),0))}
for (i in 1:ITER){SM05[,,i]<-rbind(0,cbind(diag(SxM05[,i]),0))}
for (i in 1:ITER){SM10[,,i]<-rbind(0,cbind(diag(SxM10[,i]),0))}
for (i in 1:ITER){SM15[,,i]<-rbind(0,cbind(diag(SxM15[,i]),0))}
for (i in 1:ITER){SM20[,,i]<-rbind(0,cbind(diag(SxM20[,i]),0))}
for (i in 1:ITER){SM25[,,i]<-rbind(0,cbind(diag(SxM25[,i]),0))}
for (i in 1:ITER){SM30[,,i]<-rbind(0,cbind(diag(SxM30[,i]),0))}

##SPECIAL CALCULATION FOR OPEN-ENDED AGE GROUP OF LESLIE MATRICES
for (i in 1:ITER){SF05[SIZE,SIZE,i]<-SF05[SIZE,SIZE-1,i]<-(LxF05[SIZE,i]/(LxF05[SIZE,i]+LxF05[SIZE-1,i]))}
for (i in 1:ITER){SF10[SIZE,SIZE,i]<-SF10[SIZE,SIZE-1,i]<-(LxF10[SIZE,i]/(LxF10[SIZE,i]+LxF10[SIZE-1,i]))}
for (i in 1:ITER){SF15[SIZE,SIZE,i]<-SF15[SIZE,SIZE-1,i]<-(LxF15[SIZE,i]/(LxF15[SIZE,i]+LxF15[SIZE-1,i]))}
for (i in 1:ITER){SF20[SIZE,SIZE,i]<-SF20[SIZE,SIZE-1,i]<-(LxF20[SIZE,i]/(LxF20[SIZE,i]+LxF20[SIZE-1,i]))}
for (i in 1:ITER){SF25[SIZE,SIZE,i]<-SF25[SIZE,SIZE-1,i]<-(LxF25[SIZE,i]/(LxF25[SIZE,i]+LxF25[SIZE-1,i]))}
for (i in 1:ITER){SF30[SIZE,SIZE,i]<-SF30[SIZE,SIZE-1,i]<-(LxF30[SIZE,i]/(LxF30[SIZE,i]+LxF30[SIZE-1,i]))}
for (i in 1:ITER){SM05[SIZE,SIZE,i]<-SM05[SIZE,SIZE-1,i]<-(LxM05[SIZE,i]/(LxM05[SIZE,i]+LxM05[SIZE-1,i]))}
for (i in 1:ITER){SM10[SIZE,SIZE,i]<-SM10[SIZE,SIZE-1,i]<-(LxM10[SIZE,i]/(LxM10[SIZE,i]+LxM10[SIZE-1,i]))}
for (i in 1:ITER){SM15[SIZE,SIZE,i]<-SM15[SIZE,SIZE-1,i]<-(LxM15[SIZE,i]/(LxM15[SIZE,i]+LxM15[SIZE-1,i]))}
for (i in 1:ITER){SM20[SIZE,SIZE,i]<-SM20[SIZE,SIZE-1,i]<-(LxM20[SIZE,i]/(LxM20[SIZE,i]+LxM20[SIZE-1,i]))}
for (i in 1:ITER){SM25[SIZE,SIZE,i]<-SM25[SIZE,SIZE-1,i]<-(LxM25[SIZE,i]/(LxM25[SIZE,i]+LxM25[SIZE-1,i]))}
for (i in 1:ITER){SM30[SIZE,SIZE,i]<-SM30[SIZE,SIZE-1,i]<-(LxM30[SIZE,i]/(LxM30[SIZE,i]+LxM30[SIZE-1,i]))}

##PUT FERTILITY INTO AGE PROFILES
TFRs2005<-TFRFORE[,2]
TFRs2005<-as.matrix(TFRs2005)
TFRs2010<-TFRFORE[,3]
TFRs2010<-as.matrix(TFRs2010)
TFRs2015<-TFRFORE[,4]
TFRs2015<-as.matrix(TFRs2015)
TFRs2020<-TFRFORE[,5]
TFRs2020<-as.matrix(TFRs2020)
TFRs2025<-TFRFORE[,6]
TFRs2025<-as.matrix(TFRs2025)
TFRs2030<-TFRFORE[,7]
TFRs2030<-as.matrix(TFRs2030)

TFR2005<-array(0,c(ITER,SIZE))
TFR2005[1:ITER,]<-TFRs2005
TFR2010<-array(0,c(ITER,SIZE))
TFR2010[1:ITER,]<-TFRs2010
TFR2015<-array(0,c(ITER,SIZE))
TFR2015[1:ITER,]<-TFRs2015
TFR2020<-array(0,c(ITER,SIZE))
TFR2020[1:ITER,]<-TFRs2020
TFR2025<-array(0,c(ITER,SIZE))
TFR2025[1:ITER,]<-TFRs2025
TFR2030<-array(0,c(ITER,SIZE))
TFR2030[1:ITER,]<-TFRs2030

TFR2005<-t(TFR2005)
TFR2010<-t(TFR2010)
TFR2015<-t(TFR2015)
TFR2020<-t(TFR2020)
TFR2025<-t(TFR2025)
TFR2030<-t(TFR2030)

Fert2005<-TFR2005*Fx
Fert2010<-TFR2010*Fx
Fert2015<-TFR2015*Fx
Fert2020<-TFR2020*Fx
Fert2025<-TFR2025*Fx
Fert2030<-TFR2030*Fx

##MAKE MIGRATION AGE PROFILES
IxM<-array(0,c(SIZE,ITER))
IxM[1:SIZE,]<-PropInM
IxF<-array(0,c(SIZE,ITER))
IxF[1:SIZE,]<-PropInF

OxM<-array(0,c(SIZE,ITER))
OxM[1:SIZE,]<-PropOutM
OxF<-array(0,c(SIZE,ITER))
OxF[1:SIZE,]<-PropOutF

##MAKE THE LESLIE MATRICES FOR FEMALES
BF30<-BF25<-BF20<-BF15<-BF10<-BF05<-0*SF05
  ##NOTE: THE TOP-ROW CALCULATION IS A LITTLE OFF (SEE DIVIDING BY 10 AT THE FRONT - I BELIEVE TO HELP MANAGE SEPARATE MALE AND FEMALE MATRICES)
for(j in 1:SIZE-1)
  {BF05[1,j,1:ITER]<-(LxF05[1]/10)*(Fert2005[j,]+Fert2005[j+1,]*(SxF05[j]))*ffab}
AF05 = SF05 + BF05
for(j in 1:SIZE-1)
  {BF10[1,j,1:ITER]<-(LxF10[1]/10)*(Fert2010[j,]+Fert2010[j+1,]*(SxF10[j]))*ffab}
AF10 = SF10 + BF10
for(j in 1:SIZE-1)
  {BF15[1,j,1:ITER]<-(LxF15[1]/10)*(Fert2015[j,]+Fert2015[j+1,]*(SxF15[j]))*ffab}
AF15 = SF15 + BF15
for(j in 1:SIZE-1)
  {BF20[1,j,1:ITER]<-(LxF20[1]/10)*(Fert2020[j,]+Fert2020[j+1,]*(SxF20[j]))*ffab}
AF20 = SF20 + BF20
for(j in 1:SIZE-1)
  {BF25[1,j,1:ITER]<-(LxF25[1]/10)*(Fert2025[j,]+Fert2025[j+1,]*(SxF25[j]))*ffab}
AF25 = SF25 + BF25
for(j in 1:SIZE-1)
  {BF30[1,j,1:ITER]<-(LxF30[1]/10)*(Fert2030[j,]+Fert2030[j+1,]*(SxF30[j]))*ffab}
AF30 = SF30 + BF30

##MAKE ARRAYS TO HOLD THE DATA
KF05<-array(KF05,c(SIZE,1,ITER))
KF10<-array(0,c(SIZE,1,ITER))
KF35<-KF30<-KF25<-KF20<-KF15<-KF10

##PROJECT THE FEMALE POPULATION (NATURAL INCREASE, LESS OUT MIGRATION, PLUS IN MIGRATION 
##(IN MIGRATION IS OUT MIGRATION SUM PLUS NET MIGRATION SUM)
Out2005<-array(0,c(SIZE,ITER))
for (i in 1:ITER){Out2005[,i]<-(OxF[,i]*OutRate*sum(KF05[,,i]))}
In2005<-array(0,c(SIZE,ITER))
for (i in 1:ITER){In2005[,i]<-(sum(Out2005[,i])+(sum(KF05[,,i])*TIRFORE[i,2]))*IxF[,i]}
for (i in 1:ITER){KF10[,,i]<- (AF05[,,i]%*%KF05[,,i])+(In2005[,i])-(Out2005[,i])}

Out2010<-array(0,c(SIZE,ITER))
for (i in 1:ITER){Out2010[,i]<-(OxF[,i]*OutRate*sum(KF10[,,i]))}
In2010<-array(0,c(SIZE,ITER))
for (i in 1:ITER){In2010[,i]<-(sum(Out2010[,i])+(sum(KF10[,,i])*TIRFORE[i,3]))*IxF[,i]}
for (i in 1:ITER){KF15[,,i]<- (AF10[,,i]%*%KF10[,,i])+(In2010[,i])-(Out2010[,i])}

Out2015<-array(0,c(SIZE,ITER))
for (i in 1:ITER){Out2015[,i]<-(OxF[,i]*OutRate*sum(KF15[,,i]))}
In2015<-array(0,c(SIZE,ITER))
for (i in 1:ITER){In2015[,i]<-(sum(Out2015[,i])+(sum(KF15[,,i])*TIRFORE[i,4]))*IxF[,i]}
for (i in 1:ITER){KF20[,,i]<- (AF15[,,i]%*%KF15[,,i])+(In2015[,i])-(Out2015[,i])}

Out2020<-array(0,c(SIZE,ITER))
for (i in 1:ITER){Out2020[,i]<-(OxF[,i]*OutRate*sum(KF20[,,i]))}
In2020<-array(0,c(SIZE,ITER))
for (i in 1:ITER){In2020[,i]<-(sum(Out2020[,i])+(sum(KF20[,,i])*TIRFORE[i,5]))*IxF[,i]}
for (i in 1:ITER){KF25[,,i]<- (AF20[,,i]%*%KF20[,,i])+(In2020[,i])-(Out2020[,i])}

Out2025<-array(0,c(SIZE,ITER))
for (i in 1:ITER){Out2025[,i]<-(OxF[,i]*OutRate*sum(KF25[,,i]))}
In2025<-array(0,c(SIZE,ITER))
for (i in 1:ITER){In2025[,i]<-(sum(Out2025[,i])+(sum(KF25[,,i])*TIRFORE[i,6]))*IxF[,i]}
for (i in 1:ITER){KF30[,,i]<- (AF25[,,i]%*%KF25[,,i])+(In2025[,i])-(Out2025[,i])}

Out2030<-array(0,c(SIZE,ITER))
for (i in 1:ITER){Out2030[,i]<-(OxF[,i]*OutRate*sum(KF30[,,i]))}
In2030<-array(0,c(SIZE,ITER))
for (i in 1:ITER){In2030[,i]<-(sum(Out2030[,i])+(sum(KF30[,,i])*TIRFORE[i,7]))*IxF[,i]}
for (i in 1:ITER){KF35[,,i]<- (AF30[,,i]%*%KF30[,,i])+(In2030[,i])-(Out2030[,i])}

##MAKE THE LESLIE MATRICES FOR MALES
BM30<-BM25<-BM20<-BM15<-BM10<-BM05<-0*SF05
  ##NOTE: THE TOP-ROW CALCULATION IS A LITTLE OFF (SEE DIVIDING BY 10 AT THE FRONT - I BELIEVE TO HELP MANAGE SEPARATE MALE AND FEMALE MATRICES)
for(j in 1:SIZE-1)
  {BM05[1,j,1:ITER] <- (LxM05[1]/10)*(Fert2005[j,]+Fert2005[j+1,]*(SxF05[j]))*fmab}
AM05 = SM05 + BM05
for(j in 1:SIZE-1)
  {BM10[1,j,1:ITER] <- (LxM10[1]/10)*(Fert2010[j,]+Fert2010[j+1,]*(SxF10[j]))*fmab}
AM10 = SM10 + BM10
for(j in 1:SIZE-1)
  {BM15[1,j,1:ITER] <- (LxM15[1]/10)*(Fert2015[j,]+Fert2015[j+1,]*(SxF15[j]))*fmab}
AM15 = SM15 + BM15
for(j in 1:SIZE-1)
  {BM20[1,j,1:ITER] <- (LxM20[1]/10)*(Fert2020[j,]+Fert2020[j+1,]*(SxF20[j]))*fmab}
AM20 = SM20 + BM20
for(j in 1:SIZE-1)
  {BM25[1,j,1:ITER] <- (LxM25[1]/10)*(Fert2025[j,]+Fert2025[j+1,]*(SxF25[j]))*fmab}
AM25 = SM25 + BM25
for(j in 1:SIZE-1)
  {BM30[1,j,1:ITER] <- (LxM30[1]/10)*(Fert2030[j,]+Fert2030[j+1,]*(SxF30[j]))*fmab}
AM30 = SM30 + BM30

##MAKE ARRAYS TO HOLD THE DATA
KM05<-array(KM05,c(SIZE,1,ITER))
KM10<-array(0,c(SIZE,1,ITER))
KBM10<-array(0,c(SIZE,1,ITER))
KSM10<-array(0,c(SIZE,1,ITER))
KM35<-KM30<-KM25<-KM20<-KM15<-KM10
KBM35<-KBM30<-KBM25<-KBM20<-KBM15<-KBM10
KSM35<-KSM30<-KSM25<-KSM20<-KSM15<-KSM10

##PROJECT THE MALE POPULATION (NATURAL INCREASE, LESS OUT MIGRATION, PLUS IN MIGRATION 
##(IN MIGRATION IS OUT MIGRATION SUM PLUS NET MIGRATION SUM)
OutM2005<-array(0,c(SIZE,ITER))
for (i in 1:ITER){OutM2005[,i]<-(OxM[,i]*OutRate*sum(KM05[,,i]))}
InM2005<-array(0,c(SIZE,ITER))
for (i in 1:ITER){InM2005[,i]<-(sum(OutM2005[,i])+(sum(KM05[,,i])*TIRFORE[i,2]))*IxM[,i]}
for (i in 1:ITER){KBM10[,,i]<- (BM05[,,i]%*%KF05[,,i])}
for (i in 1:ITER){KSM10[,,i]<- (SM05[,,i]%*%KM05[,,i])+(InM2005[,i])-(OutM2005[,i])}
for (i in 1:ITER){KM10[,,i]<-(KBM10[,,i]+KSM10[,,i])}

OutM2010<-array(0,c(SIZE,ITER))
for (i in 1:ITER){OutM2010[,i]<-(OxM[,i]*OutRate*sum(KM10[,,i]))}
InM2010<-array(0,c(SIZE,ITER))
for (i in 1:ITER){InM2010[,i]<-(sum(OutM2010[,i])+(sum(KM10[,,i])*TIRFORE[i,3]))*IxM[,i]}
for (i in 1:ITER){KBM15[,,i]<- (BM10[,,i]%*%KF10[,,i])}
for (i in 1:ITER){KSM15[,,i]<- (SM10[,,i]%*%KM10[,,i])+(InM2010[,i])-(OutM2010[,i])}
for (i in 1:ITER){KM15[,,i]<-(KBM15[,,i]+KSM15[,,i])}

OutM2015<-array(0,c(SIZE,ITER))
for (i in 1:ITER){OutM2015[,i]<-(OxM[,i]*OutRate*sum(KM15[,,i]))}
InM2015<-array(0,c(SIZE,ITER))
for (i in 1:ITER){InM2015[,i]<-(sum(OutM2015[,i])+(sum(KM15[,,i])*TIRFORE[i,4]))*IxM[,i]}
for (i in 1:ITER){KBM20[,,i]<- (BM15[,,i]%*%KF15[,,i])}
for (i in 1:ITER){KSM20[,,i]<- (SM15[,,i]%*%KM15[,,i])+(InM2015[,i])-(OutM2015[,i])}
for (i in 1:ITER){KM20[,,i]<-(KBM20[,,i]+KSM20[,,i])}

OutM2020<-array(0,c(SIZE,ITER))
for (i in 1:ITER){OutM2020[,i]<-(OxM[,i]*OutRate*sum(KM20[,,i]))}
InM2020<-array(0,c(SIZE,ITER))
for (i in 1:ITER){InM2020[,i]<-(sum(OutM2020[,i])+(sum(KM20[,,i])*TIRFORE[i,5]))*IxM[,i]}
for (i in 1:ITER){KBM25[,,i]<- (BM20[,,i]%*%KF20[,,i])}
for (i in 1:ITER){KSM25[,,i]<- (SM20[,,i]%*%KM20[,,i])+(InM2020[,i])-(OutM2020[,i])}
for (i in 1:ITER){KM25[,,i]<-(KBM25[,,i]+KSM25[,,i])}

OutM2025<-array(0,c(SIZE,ITER))
for (i in 1:ITER){OutM2025[,i]<-(OxM[,i]*OutRate*sum(KM25[,,i]))}
InM2025<-array(0,c(SIZE,ITER))
for (i in 1:ITER){InM2025[,i]<-(sum(OutM2025[,i])+(sum(KM25[,,i])*TIRFORE[i,6]))*IxM[,i]}
for (i in 1:ITER){KBM30[,,i]<- (BM25[,,i]%*%KF25[,,i])}
for (i in 1:ITER){KSM30[,,i]<- (SM25[,,i]%*%KM25[,,i])+(InM2025[,i])-(OutM2025[,i])}
for (i in 1:ITER){KM30[,,i]<-(KBM30[,,i]+KSM30[,,i])}

OutM2030<-array(0,c(SIZE,ITER))
for (i in 1:ITER){OutM2030[,i]<-(OxM[,i]*OutRate*sum(KM30[,,i]))}
InM2030<-array(0,c(SIZE,ITER))
for (i in 1:ITER){InM2030[,i]<-(sum(OutM2030[,i])+(sum(KM30[,,i])*TIRFORE[i,7]))*IxM[,i]}
for (i in 1:ITER){KBM35[,,i]<- (BM30[,,i]%*%KF30[,,i])}
for (i in 1:ITER){KSM35[,,i]<- (SM30[,,i]%*%KM30[,,i])+(InM2030[,i])-(OutM2030[,i])}
for (i in 1:ITER){KM35[,,i]<-(KBM35[,,i]+KSM35[,,i])}

##MAKE TABLES OF DATA OF INTEREST
Netx2005<-In2005-Out2005+InM2005-OutM2005
Netx2010<-In2010-Out2010+InM2010-OutM2010
Netx2015<-In2015-Out2015+InM2015-OutM2015
Netx2020<-In2020-Out2020+InM2020-OutM2020
Netx2025<-In2025-Out2025+InM2025-OutM2025
Netx2030<-In2030-Out2030+InM2030-OutM2030

Net2005<-array(0,ITER)
for (i in 1:ITER){Net2005[i]<-(sum(Netx2005[,i]))}
Net2010<-array(0,ITER)
for (i in 1:ITER){Net2010[i]<-(sum(Netx2010[,i]))}
Net2015<-array(0,ITER)
for (i in 1:ITER){Net2015[i]<-(sum(Netx2015[,i]))}
Net2020<-array(0,ITER)
for (i in 1:ITER){Net2020[i]<-(sum(Netx2020[,i]))}
Net2025<-array(0,ITER)
for (i in 1:ITER){Net2025[i]<-(sum(Netx2025[,i]))}
Net2030<-array(0,ITER)
for (i in 1:ITER){Net2030[i]<-(sum(Netx2030[,i]))}
Net<-c(Net2005,Net2010,Net2015,Net2020,Net2025,Net2030)

KT35<-KT30<-KT25<-KT20<-KT15<-KT10<-KT05<-array(0,c(ITER))
for (i in 1:ITER){KT05[i]<-sum(KF05[,,i])+sum(KM05[,,i])}
for (i in 1:ITER){KT10[i]<-sum(KF10[,,i])+sum(KM10[,,i])}
for (i in 1:ITER){KT15[i]<-sum(KF15[,,i])+sum(KM15[,,i])}
for (i in 1:ITER){KT20[i]<-sum(KF20[,,i])+sum(KM20[,,i])}
for (i in 1:ITER){KT25[i]<-sum(KF25[,,i])+sum(KM25[,,i])}
for (i in 1:ITER){KT30[i]<-sum(KF30[,,i])+sum(KM30[,,i])}
for (i in 1:ITER){KT35[i]<-sum(KF35[,,i])+sum(KM35[,,i])}
KT<-array(c(KT05,KT10,KT15,KT20,KT25,KT30,KT35),c(ITER,7))

KY35<-KY30<-KY25<-KY20<-KY15<-KY10<-KY05<-array(0,c(ITER))
for (i in 1:ITER){KY05[i]<-sum(KF05[0:4,,i])+sum(KM05[0:4,,i])}
for (i in 1:ITER){KY10[i]<-sum(KF10[0:4,,i])+sum(KM10[0:4,,i])}
for (i in 1:ITER){KY15[i]<-sum(KF15[0:4,,i])+sum(KM15[0:4,,i])}
for (i in 1:ITER){KY20[i]<-sum(KF20[0:4,,i])+sum(KM20[0:4,,i])}
for (i in 1:ITER){KY25[i]<-sum(KF25[0:4,,i])+sum(KM25[0:4,,i])}
for (i in 1:ITER){KY30[i]<-sum(KF30[0:4,,i])+sum(KM30[0:4,,i])}
for (i in 1:ITER){KY35[i]<-sum(KF35[0:4,,i])+sum(KM35[0:4,,i])}
KY<-array(c(KY05,KY10,KY15,KY20,KY25,KY30,KY35),c(ITER,7))

KW35<-KW30<-KW25<-KW20<-KW15<-KW10<-KW05<-array(0,c(ITER))
for (i in 1:ITER){KW05[i]<-sum(KF05[5:13,,i])+sum(KM05[5:13,,i])}
for (i in 1:ITER){KW10[i]<-sum(KF10[5:13,,i])+sum(KM10[5:13,,i])}
for (i in 1:ITER){KW15[i]<-sum(KF15[5:13,,i])+sum(KM15[5:13,,i])}
for (i in 1:ITER){KW20[i]<-sum(KF20[5:13,,i])+sum(KM20[5:13,,i])}
for (i in 1:ITER){KW25[i]<-sum(KF25[5:13,,i])+sum(KM25[5:13,,i])}
for (i in 1:ITER){KW30[i]<-sum(KF30[5:13,,i])+sum(KM30[5:13,,i])}
for (i in 1:ITER){KW35[i]<-sum(KF35[5:13,,i])+sum(KM35[5:13,,i])}
KW<-array(c(KW05,KW10,KW15,KW20,KW25,KW30,KW35),c(ITER,7))

KO35<-KO30<-KO25<-KO20<-KO15<-KO10<-KO05<-array(0,c(ITER))
for (i in 1:ITER){KO05[i]<-sum(KF05[14:21,,i])+sum(KM05[14:21,,i])}
for (i in 1:ITER){KO10[i]<-sum(KF10[14:21,,i])+sum(KM10[14:21,,i])}
for (i in 1:ITER){KO15[i]<-sum(KF15[14:21,,i])+sum(KM15[14:21,,i])}
for (i in 1:ITER){KO20[i]<-sum(KF20[14:21,,i])+sum(KM20[14:21,,i])}
for (i in 1:ITER){KO25[i]<-sum(KF25[14:21,,i])+sum(KM25[14:21,,i])}
for (i in 1:ITER){KO30[i]<-sum(KF30[14:21,,i])+sum(KM30[14:21,,i])}
for (i in 1:ITER){KO35[i]<-sum(KF35[14:21,,i])+sum(KM35[14:21,,i])}
KO<-array(c(KO05,KO10,KO15,KO20,KO25,KO30,KO35),c(ITER,7))

OA35<-OA30<-OA25<-OA20<-OA15<-OA10<-OA05<-array(0,c(ITER))
for (i in 1:ITER){OA05[i]<-(sum(KF05[14:21,,i])+sum(KM05[14:21,,i]))/(sum(KF05[5:13,,i])+sum(KM05[5:13,,i]))}
for (i in 1:ITER){OA10[i]<-(sum(KF10[14:21,,i])+sum(KM10[14:21,,i]))/(sum(KF10[5:13,,i])+sum(KM10[5:13,,i]))}
for (i in 1:ITER){OA15[i]<-(sum(KF15[14:21,,i])+sum(KM15[14:21,,i]))/(sum(KF15[5:13,,i])+sum(KM15[5:13,,i]))}
for (i in 1:ITER){OA20[i]<-(sum(KF20[14:21,,i])+sum(KM20[14:21,,i]))/(sum(KF20[5:13,,i])+sum(KM20[5:13,,i]))}
for (i in 1:ITER){OA25[i]<-(sum(KF25[14:21,,i])+sum(KM25[14:21,,i]))/(sum(KF25[5:13,,i])+sum(KM25[5:13,,i]))}
for (i in 1:ITER){OA30[i]<-(sum(KF30[14:21,,i])+sum(KM30[14:21,,i]))/(sum(KF30[5:13,,i])+sum(KM30[5:13,,i]))}
for (i in 1:ITER){OA35[i]<-(sum(KF35[14:21,,i])+sum(KM35[14:21,,i]))/(sum(KF35[5:13,,i])+sum(KM35[5:13,,i]))}
OA<-array(c(OA05,OA10,OA15,OA20,OA25,OA30,OA35),c(ITER,7))

YA35<-YA30<-YA25<-YA20<-YA15<-YA10<-YA05<-array(0,c(ITER))
for (i in 1:ITER){YA05[i]<-(sum(KF05[1:4,,i])+sum(KM05[1:4,,i]))/(sum(KF05[5:13,,i])+sum(KM05[5:13,,i]))}
for (i in 1:ITER){YA10[i]<-(sum(KF10[1:4,,i])+sum(KM10[1:4,,i]))/(sum(KF10[5:13,,i])+sum(KM10[5:13,,i]))}
for (i in 1:ITER){YA15[i]<-(sum(KF15[1:4,,i])+sum(KM15[1:4,,i]))/(sum(KF15[5:13,,i])+sum(KM15[5:13,,i]))}
for (i in 1:ITER){YA20[i]<-(sum(KF20[1:4,,i])+sum(KM20[1:4,,i]))/(sum(KF20[5:13,,i])+sum(KM20[5:13,,i]))}
for (i in 1:ITER){YA25[i]<-(sum(KF25[1:4,,i])+sum(KM25[1:4,,i]))/(sum(KF25[5:13,,i])+sum(KM25[5:13,,i]))}
for (i in 1:ITER){YA30[i]<-(sum(KF30[1:4,,i])+sum(KM30[1:4,,i]))/(sum(KF30[5:13,,i])+sum(KM30[5:13,,i]))}
for (i in 1:ITER){YA35[i]<-(sum(KF35[1:4,,i])+sum(KM35[1:4,,i]))/(sum(KF35[5:13,,i])+sum(KM35[5:13,,i]))}
YA<-array(c(YA05,YA10,YA15,YA20,YA25,YA30,YA35),c(ITER,7))

TA35<-TA30<-TA25<-TA20<-TA15<-TA10<-TA05<-array(0,c(ITER))
for (i in 1:ITER){TA05[i]<-(sum(KF05[1:4,,i])+sum(KM05[1:4,,i])+sum(KF05[14:21,,i])+sum(KM05[14:21,,i]))/(sum(KF05[5:13,,i])+sum(KM05[5:13,,i]))}
for (i in 1:ITER){TA10[i]<-(sum(KF10[1:4,,i])+sum(KM10[1:4,,i])+sum(KF10[14:21,,i])+sum(KM10[14:21,,i]))/(sum(KF10[5:13,,i])+sum(KM10[5:13,,i]))}
for (i in 1:ITER){TA15[i]<-(sum(KF15[1:4,,i])+sum(KM15[1:4,,i])+sum(KF15[14:21,,i])+sum(KM15[14:21,,i]))/(sum(KF15[5:13,,i])+sum(KM15[5:13,,i]))}
for (i in 1:ITER){TA20[i]<-(sum(KF20[1:4,,i])+sum(KM20[1:4,,i])+sum(KF20[14:21,,i])+sum(KM20[14:21,,i]))/(sum(KF20[5:13,,i])+sum(KM20[5:13,,i]))}
for (i in 1:ITER){TA25[i]<-(sum(KF25[1:4,,i])+sum(KM25[1:4,,i])+sum(KF25[14:21,,i])+sum(KM25[14:21,,i]))/(sum(KF25[5:13,,i])+sum(KM25[5:13,,i]))}
for (i in 1:ITER){TA30[i]<-(sum(KF30[1:4,,i])+sum(KM30[1:4,,i])+sum(KF30[14:21,,i])+sum(KM30[14:21,,i]))/(sum(KF30[5:13,,i])+sum(KM30[5:13,,i]))}
for (i in 1:ITER){TA35[i]<-(sum(KF35[1:4,,i])+sum(KM35[1:4,,i])+sum(KF35[14:21,,i])+sum(KM35[14:21,,i]))/(sum(KF35[5:13,,i])+sum(KM35[5:13,,i]))}
TA<-array(c(TA05,TA10,TA15,TA20,TA25,TA30,TA35),c(ITER,7))

GR0535<-array(0,c(ITER))
for (i in 1:ITER){GR0535[i]<-(1/(STEPS*5))*log(KT35[i]/KT05[i])}

##############################################################################################################################
##############################################################################################################################
##FIGURES AND TABLES
##############################################################################################################################
##############################################################################################################################

#FIGURE 1
TFR<-array(0,c(ITER,13))
for(i in 1:ITER){TFR[i,]<-c(2.58,2.16,2.40,2.35,2.58,2.48,TFRFORE[i,])}
for(i in 1:50){plot(TFR[i,],type="l",ylim=c(.5,4.5),xlim=c(0,13),col="green",xlab="Time Period",ylab="Total Fertility Rate",axes=F,lwd=1)}
axis(side=1,at=1:13,labels=c("1970-75","1975-80","1980-85","1985-90","1990-95","1995-00","2000-05","2005-10","2010-15","2015-20","2020-25","2025-30","2030-35"),cex.axis=0.8)
axis(side=2,cex.axis=0.8)
for(i in 1:50){points(TFR[i,],type="l",col=sample(6),lwd=1)}
ksdl<-array(0,13)
for (i in 1:13){ksdl[i]<-(quantile(TFR[,i],.05))}
ksdu<-array(0,13)
for (i in 1:13){ksdu[i]<-(quantile(TFR[,i],.95))}
ksdl2<-array(0,13)
for (i in 1:13){ksdl2[i]<-(quantile(TFR[,i],.05))}
ksdu2<-array(0,13)
for (i in 1:13){ksdu2[i]<-(quantile(TFR[,i],.95))}
points(ksdl,type="l",lwd=8)
points(ksdu,type="l",lwd=8)
points(ksdl2,type="l",lwd=8)
points(ksdu2,type="l",lwd=8)
title("TOTAL FERTILITY RATE: HISTORICAL AND 50 FORECASTS",cex.main=1)
mtext(side=3,line=.75,text="WITH 90% CONFIDENCE INTERVALS (BOLD BLACK)",font=2,cex=1)
dev.copy2eps(file="C:/Users/Eddie/Desktop/Forecast_AK/Figure1.eps") 
Sys.sleep(5)

#FIGURE 2
TNRFORE<-array(0,c(ITER,13))
for(i in 1:ITER){TNRFORE[i,]<-c(0.15387,0.00229,0.18100,-0.07252,0.00713,-0.01774,TIRFORE[i,])}
for(i in 1:50){plot(TNRFORE[i,],type="l",ylim=c(-.2,.2),xlim=c(0,13),col="green",xlab="Time Period",ylab="Net Migration Ratio",axes=F)}
axis(side=1,at=1:13,labels=c("1970-75","1975-80","1980-85","1985-90","1990-95","1995-00","2000-05","2005-10","2010-15","2015-20","2020-25","2025-30","2030-35"),cex.axis=0.8)
axis(side=2,cex.axis=0.8)
for(i in 1:50){points(TNRFORE[i,],type="l",col=sample(6),lwd=1)}
ksdl<-array(0,13)
for (i in 1:13){ksdl[i]<-(quantile(TNRFORE[,i],.05))}
ksdu<-array(0,13)
for (i in 1:13){ksdu[i]<-(quantile(TNRFORE[,i],.95))}
ksdl2<-array(0,13)
for (i in 1:13){ksdl2[i]<-(quantile(TNRFORE[,i],.05))}
ksdu2<-array(0,13)
for (i in 1:13){ksdu2[i]<-(quantile(TNRFORE[,i],.95))}
points(ksdl,type="l",lwd=8)
points(ksdu,type="l",lwd=8)
points(ksdl2,type="l",lwd=8)
points(ksdu2,type="l",lwd=8)
title("NET MIGRATION RATIO: HISTORICAL AND 50 FORECASTS",cex.main=1)
mtext(side=3,line=.75,text="WITH 90% CONFIDENCE INTERVALS (BOLD BLACK)",font=2,cex=1)
dev.copy2eps(file="C:/Users/Eddie/Desktop/Forecast_AK/Figure2.eps") 
Sys.sleep(5)

#FIGURE 3
Net<-array(0,c(ITER,13))
for(i in 1:ITER){Net[i,]<-c(47469,881,75984,-39444,3942,-10674,1601,Net2005[i],Net2010[i],Net2015[i],Net2020[i],Net2025[i],Net2030[i])}
for(i in 1:50){plot(Net[i,],type="l",ylim=c(-100000,100000),xlim=c(0,13),col="green",xlab="Time Period",ylab="Net Migration",axes=F,lwd=1)}
axis(side=1,at=1:13,labels=c("1970-75","1975-80","1980-85","1985-90","1990-95","1995-00","2000-05","2005-10","2010-15","2015-20","2020-25","2025-30","2030-35"),cex.axis=0.8)
axis(side=2,cex.axis=0.8)
for(i in 1:50){points(Net[i,],type="l",col=sample(6),lwd=1)}
ksdl<-array(0,13)
for (i in 1:13){ksdl[i]<-(quantile(Net[,i],.05))}
ksdu<-array(0,13)
for (i in 1:13){ksdu[i]<-(quantile(Net[,i],.95))}
ksdl2<-array(0,13)
for (i in 1:13){ksdl2[i]<-(quantile(Net[,i],.05))}
ksdu2<-array(0,13)
for (i in 1:13){ksdu2[i]<-(quantile(Net[,i],.95))}
points(ksdl,type="l",lwd=8)
points(ksdu,type="l",lwd=8)
points(ksdl2,type="l",lwd=8)
points(ksdu2,type="l",lwd=8)
title("NET MIGRATION: HISTORICAL AND 50 FORECASTS",cex.main=1)
mtext(side=3,line=.75,text="WITH 90% CONFIDENCE INTERVALS (BOLD BLACK)",font=2,cex=1)
dev.copy2eps(file="C:/Users/Eddie/Desktop/Forecast_AK/Figure3.eps") 
Sys.sleep(5)

#FIGURE 4
e0F<-array(0,c(13,ITER))
for(i in 1:ITER){e0F[,i]<-c(74.6,75.8,77.0,78.1,78.9,79.4,79.9,e0FFORE[2:7,i])}
for(i in 1:50){plot(e0F[,i],type="l",ylim=c(65,90),xlim=c(0,13),col="blue",xlab="Time Period",ylab="Life Expectancy at Birth",axes=F,lwd=1)}
axis(side=1,at=1:13,labels=c("1970-75","1975-80","1980-85","1985-90","1990-95","1995-00","2000-05","2005-10","2010-15","2015-20","2020-25","2025-30","2030-35"),cex.axis=0.8)
axis(side=2,cex.axis=0.8)
for(i in 1:50){points(e0F[,i],type="l",ylim=c(70,85),xlim=c(0,13),col=sample(6),lwd=1)}
ksdl<-array(0,13)
for (i in 1:13){ksdl[i]<-(quantile(e0F[i,],.05))}
ksdu<-array(0,13)
for (i in 1:13){ksdu[i]<-(quantile(e0F[i,],.95))}
ksdl2<-array(0,13)
for (i in 1:13){ksdl2[i]<-(quantile(e0F[i,],.05))}
ksdu2<-array(0,13)
for (i in 1:13){ksdu2[i]<-(quantile(e0F[i,],.95))}
points(ksdl,type="l",lwd=8)
points(ksdu,type="l",lwd=8)
points(ksdl2,type="l",lwd=8)
points(ksdu2,type="l",lwd=8)
e0M<-array(0,c(13,ITER))
for(i in 1:ITER){e0M[,i]<-c(66.8,68.1,69.5,70.9,72.4,74.0,75.2,e0MFORE[2:7,i])}
for(i in 1:50){points(e0M[,i],type="l",ylim=c(65,80),xlim=c(0,13),col=sample(6),lwd=1)}
ksdl<-array(0,13)
for (i in 1:13){ksdl[i]<-(quantile(e0M[i,],.05))}
ksdu<-array(0,13)
for (i in 1:13){ksdu[i]<-(quantile(e0M[i,],.95))}
ksdl2<-array(0,13)
for (i in 1:13){ksdl2[i]<-(quantile(e0M[i,],.05))}
ksdu2<-array(0,13)
for (i in 1:13){ksdu2[i]<-(quantile(e0M[i,],.95))}
points(ksdl,type="l",lwd=8)
points(ksdu,type="l",lwd=8)
points(ksdl2,type="l",lwd=8)
points(ksdu2,type="l",lwd=8)
title("LIFE EXPECTANCY AT BIRTH: HISTORICAL AND 50 FORECASTS",cex.main=1)
mtext(side=3,line=.75,text="WITH 90% CONFIDENCE INTERVALS (BOLD BLACK)",font=2,cex=1)
mtext(side=3,line=-.5,text="(Historical estimates are interpolated by the author from estimates for 1970, 1980, 1990 and 2000.)",font=1,cex=.75)
mtext(side=1,line=-18,text="Female",font=2,cex=1)
mtext(side=1,line=-10,adj=.6,text="Male",font=2,cex=1)
dev.copy2eps(file="C:/Users/Eddie/Desktop/Forecast_AK/Figure4.eps") 
Sys.sleep(5)

#FIGURE 5
for(i in 1:50){plot(Fert2030[,i],type="l",ylim=c(0,1.5),xlim=c(0,SIZE),col="purple",xlab="Age",ylab="Fertility Rate",axes=F,lwd=1)}
axis(side=1,at=1:21,labels=seq(0,100,5),cex.axis=0.8)
axis(side=2,cex.axis=0.8)
for(i in 1:50){points(Fert2030[,i],type="l",col=sample(6),lwd=1)}
ksdl<-array(0,21)
for (i in 1:21){ksdl[i]<-(quantile(Fert2030[i,],.05))}
ksdu<-array(0,21)
for (i in 1:21){ksdu[i]<-(quantile(Fert2030[i,],.95))}
ksdl2<-array(0,21)
for (i in 1:21){ksdl2[i]<-(quantile(Fert2030[i,],.05))}
ksdu2<-array(0,21)
for (i in 1:21){ksdu2[i]<-(quantile(Fert2030[i,],.95))}
points(ksdl,type="l",lwd=8)
points(ksdu,type="l",lwd=8)
points(ksdl2,type="l",lwd=8)
points(ksdu2,type="l",lwd=8)
title("50 AGE SPECIFIC FERTILITY RATE FORECASTS FOR 2030-2035",cex.main=1)
mtext(side=3,line=.75,text="WITH 90% CONFIDENCE INTERVALS (BOLD BLACK)",font=2,cex=1)
dev.copy2eps(file="C:/Users/Eddie/Desktop/Forecast_AK/Figure5.eps") 
Sys.sleep(5)

#FIGURE 6
for(i in 1:50){plot(Netx2030[,i],type="l",ylim=c(-15000,20000),xlim=c(0,SIZE),col="purple",xlab="Age",ylab="Net Migration",axes=F,lwd=1)}
axis(side=1,at=1:21,labels=seq(0,100,5),cex.axis=0.8)
axis(side=2,cex.axis=0.8)
for(i in 1:50){points(Netx2030[,i],type="l",col=sample(6),lwd=1)}
ksdl<-array(0,21)
for (i in 1:21){ksdl[i]<-(quantile(Netx2030[i,],.05))}
ksdu<-array(0,21)
for (i in 1:21){ksdu[i]<-(quantile(Netx2030[i,],.95))}
ksdl2<-array(0,21)
for (i in 1:21){ksdl2[i]<-(quantile(Netx2030[i,],.05))}
ksdu2<-array(0,21)
for (i in 1:21){ksdu2[i]<-(quantile(Netx2030[i,],.95))}
points(ksdl,type="l",lwd=8)
points(ksdu,type="l",lwd=8)
points(ksdl2,type="l",lwd=8)
points(ksdu2,type="l",lwd=8)
title("50 NET MIGRATION BY AGE FORECASTS FOR 2030-2035",cex.main=1)
mtext(side=3,line=.75,text="WITH 90% CONFIDENCE INTERVALS (BOLD BLACK)",font=2,cex=1)
dev.copy2eps(file="C:/Users/Eddie/Desktop/Forecast_AK/Figure6.eps") 
Sys.sleep(5)

#FIGURE 7
for(i in 1:50){plot(lxF30[,i],type="l",ylim=c(0,1),xlim=c(0,SIZE+1),col="purple",xlab="Age",ylab="Survivorship",axes=F,lwd=1)}
axis(side=1,at=1:22,labels=seq(0,105,5),cex.axis=0.8)
axis(side=2,cex.axis=0.8)
for(i in 1:50){points(lxF30[,i],type="l",col=sample(6),lwd=1)}
ksdl<-array(0,22)
for (i in 1:22){ksdl[i]<-(quantile(lxF30[i,],.05))}
ksdu<-array(0,22)
for (i in 1:22){ksdu[i]<-(quantile(lxF30[i,],.95))}
ksdl2<-array(0,22)
for (i in 1:22){ksdl2[i]<-(quantile(lxF30[i,],.05))}
ksdu2<-array(0,22)
for (i in 1:22){ksdu2[i]<-(quantile(lxF30[i,],.95))}
points(ksdl,type="l",lwd=5)
points(ksdu,type="l",lwd=5)
points(ksdl2,type="l",lwd=5)
points(ksdu2,type="l",lwd=5)
for(i in 1:50){points(lxM30[,i],type="l",col=sample(6),lwd=1)}
ksdl<-array(0,22)
for (i in 1:22){ksdl[i]<-(quantile(lxM30[i,],.05))}
ksdu<-array(0,22)
for (i in 1:22){ksdu[i]<-(quantile(lxM30[i,],.95))}
ksdl2<-array(0,22)
for (i in 1:22){ksdl2[i]<-(quantile(lxM30[i,],.05))}
ksdu2<-array(0,22)
for (i in 1:22){ksdu2[i]<-(quantile(lxM30[i,],.95))}
points(ksdl,type="l",lwd=4)
points(ksdu,type="l",lwd=4)
points(ksdl2,type="l",lwd=4)
points(ksdu2,type="l",lwd=4)
title("50 PERIOD LIFE TABLE SURVIVORSHIP FORECASTS FOR 2030-2035",cex.main=1)
mtext(side=3,line=.75,text="WITH 90% CONFIDENCE INTERVALS (BOLD BLACK)",font=2,cex=1)
mtext(side=1,line=-22.5,adj=.79,text="Female",font=2,cex=1)
mtext(side=1,line=-16.5,adj=.65,text="Male",font=2,cex=1)
dev.copy2eps(file="C:/Users/Eddie/Desktop/Forecast_AK/Figure7.eps") 
Sys.sleep(5)

#FIGURE 8
KTH<-array(0,c(ITER,14))
for(i in 1:ITER){KTH[i,]<-c(308500,384100,419800,543900,553171,601581,627533,664334,KT[i,2:7])}
for(i in 1:50){plot(KTH[i,],type="l",ylim=c(0,mean(KT35)*1.5),xlim=c(0,14),col="black",xlab="Year",ylab="Population",axes=F,lwd=1)}
axis(side=1,at=1:14,labels=c("1970","1975","1980","1985","1990","1995","2000","2005","2010","2015","2020","2025","2030","2035"),cex.axis=0.8)
axis(side=2,cex.axis=0.8)
for(i in 1:50){points(KTH[i,],type="l",xlim=c(0,14),col=sample(6),lwd=1)}
ksdl<-array(0,14)
for (i in 1:14){ksdl[i]<-(quantile(KTH[,i],.05))}
ksdu<-array(0,14)
for (i in 1:14){ksdu[i]<-(quantile(KTH[,i],.95))}
ksdl2<-array(0,14)
for (i in 1:14){ksdl2[i]<-(quantile(KTH[,i],.05))}
ksdu2<-array(0,14)
for (i in 1:14){ksdu2[i]<-(quantile(KTH[,i],.95))}
points(ksdl,type="l",lwd=8)
points(ksdu,type="l",lwd=8)
points(ksdl2,type="l",lwd=8)
points(ksdu2,type="l",lwd=8)
title("TOTAL POPULATION: HISTORICAL AND 50 FORECASTS",cex.main=1)
mtext(side=3,line=.75,text="WITH 90% CONFIDENCE INTERVALS (BOLD BLACK)",font=2,cex=1)
dev.copy2eps(file="C:/Users/Eddie/Desktop/Forecast_AK/Figure8.eps") 
Sys.sleep(5)

#FIGURE 9
hist(GR0535*100,50,xlim=c(-4,4),col="gold",xlab="2005 to 2035 Growth Rate",ylab="Frequency",main="TOTAL POPULATION GROWTH RATE: HISTOGRAM",cex.main=1,lwd=1,axes=F)
axis(side=1,at=-4:4,labels=c("-4%","-3%","-2%","-1%","0%","1%","2%","3%","4%"),cex.axis=0.8)
axis(side=2,cex.axis=0.8)
mtext(side=3,line=.75,text="OF THE FORECASTS FOR 2005 THROUGH 2035",font=2,cex=1)
dev.copy2eps(file="C:/Users/Eddie/Desktop/Forecast_AK/Figure9.eps") 
Sys.sleep(5)

#FIGURE 10
YAH<-array(0,c(ITER,7))
for(i in 1:ITER){YAH[i,]<-c(YA[i,1:7])}
for(i in 1:50){plot(YAH[i,],type="l",ylim=c(0,1),xlim=c(0,7),col="black",xlab="Year",ylab="Young-Age Dependency Ratio",axes=F,lwd=1)}
axis(side=1,at=1:7,labels=c("2005","2010","2015","2020","2025","2030","2035"),cex.axis=0.8)
axis(side=2,cex.axis=0.8)
for(i in 1:50){points(YAH[i,],type="l",xlim=c(0,7),col=sample(6),lwd=1)}
ksdl<-array(0,7)
for (i in 1:7){ksdl[i]<-(quantile(YAH[,i],.05))}
ksdu<-array(0,7)
for (i in 1:7){ksdu[i]<-(quantile(YAH[,i],.95))}
ksdl2<-array(0,7)
for (i in 1:7){ksdl2[i]<-(quantile(YAH[,i],.05))}
ksdu2<-array(0,7)
for (i in 1:7){ksdu2[i]<-(quantile(YAH[,i],.95))}
points(ksdl,type="l",lwd=8)
points(ksdu,type="l",lwd=8)
points(ksdl2,type="l",lwd=8)
points(ksdu2,type="l",lwd=8)
title("YOUNG-AGE DEPENDENCY RATIO: 50 FORECASTS",cex.main=1)
mtext(side=3,line=.75,text="WITH 90% CONFIDENCE INTERVALS (BOLD BLACK)",font=2,cex=1)
dev.copy2eps(file="C:/Users/Eddie/Desktop/Forecast_AK/Figure10.eps") 
Sys.sleep(5)

#FIGURE 11
OAH<-array(0,c(ITER,7))
for(i in 1:ITER){OAH[i,]<-c(OA[i,1:7])}
for(i in 1:50){plot(OAH[i,],type="l",ylim=c(0,1),xlim=c(0,7),col="black",xlab="Year",ylab="Old-Age Dependency Ratio",axes=F,lwd=1)}
axis(side=1,at=1:7,labels=c("2005","2010","2015","2020","2025","2030","2035"),cex.axis=0.8)
axis(side=2,cex.axis=0.8)
for(i in 1:50){points(OAH[i,],type="l",xlim=c(0,7),col=sample(6),lwd=1)}
ksdl<-array(0,7)
for (i in 1:7){ksdl[i]<-(quantile(OAH[,i],.05))}
ksdu<-array(0,7)
for (i in 1:7){ksdu[i]<-(quantile(OAH[,i],.95))}
ksdl2<-array(0,7)
for (i in 1:7){ksdl2[i]<-(quantile(OAH[,i],.05))}
ksdu2<-array(0,7)
for (i in 1:7){ksdu2[i]<-(quantile(OAH[,i],.95))}
points(ksdl,type="l",lwd=8)
points(ksdu,type="l",lwd=8)
points(ksdl2,type="l",lwd=8)
points(ksdu2,type="l",lwd=8)
title("OLD-AGE DEPENDENCY RATIO: 50 FORECASTS",cex.main=1)
mtext(side=3,line=.75,text="WITH 90% CONFIDENCE INTERVALS (BOLD BLACK)",font=2,cex=1)
dev.copy2eps(file="C:/Users/Eddie/Desktop/Forecast_AK/Figure11.eps") 
Sys.sleep(5)

#FIGURE *
TAH<-array(0,c(ITER,7))
for(i in 1:ITER){TAH[i,]<-c(TA[i,1:7])}
for(i in 1:50){plot(TAH[i,],type="l",ylim=c(0,1),xlim=c(0,7),col="black",xlab="Year",ylab="Total Dependency Ratio",axes=F,lwd=1)}
axis(side=1,at=1:7,labels=c("2005","2010","2015","2020","2025","2030","2035"),cex.axis=0.8)
axis(side=2,cex.axis=0.8)
for(i in 1:50){points(TAH[i,],type="l",xlim=c(0,7),col=sample(6),lwd=1)}
ksdl<-array(0,7)
for (i in 1:7){ksdl[i]<-(quantile(TAH[,i],.05))}
ksdu<-array(0,7)
for (i in 1:7){ksdu[i]<-(quantile(TAH[,i],.95))}
ksdl2<-array(0,7)
for (i in 1:7){ksdl2[i]<-(quantile(TAH[,i],.05))}
ksdu2<-array(0,7)
for (i in 1:7){ksdu2[i]<-(quantile(TAH[,i],.95))}
points(ksdl,type="l",lwd=8)
points(ksdu,type="l",lwd=8)
points(ksdl2,type="l",lwd=8)
points(ksdu2,type="l",lwd=8)
title("TOTAL DEPENDENCY RATIO: 50 FORECASTS",cex.main=1)
mtext(side=3,line=.75,text="WITH 90% CONFIDENCE INTERVALS (BOLD BLACK)",font=2,cex=1)
#dev.copy2eps(file="C:/Users/Eddie/Desktop/Forecast_AK/Figure*.eps") 
Sys.sleep(5)

#FIGURE 12 (CARL MASON'S GREAT PYRAMID FUNCTION, WITH CI'S ADDED-- NEED GPLOTS PACKAGE)
library(gplots)
poppyr3<-function(male,female,malesdl,femalesdl,malesdu,femalesdu,cat){
split.screen(figs=rbind(c(0,.58,0,1),c(.43,1,0,1)))
screen(1)
barplot2(male,horiz=T,names=cat,plot.ci = TRUE,ci.l=malesdl,ci.u=malesdu,ci.color ="black",ci.lwd=1,space=0,
xlim=c(50000,0),col="dodgerblue")
title("Male",line=-3,cex.main=1)
screen(2)
barplot2(female,horiz=T,axisnames=F,plot.ci=TRUE,ci.l=femalesdl,ci.u=femalesdu,ci.color="black",ci.lwd=1,space=0,
xlim=c(0,50000),col="gold")
title("Female",line=-3,cex.main=1)
close.screen(all=T)}
male<-array(0,c(SIZE))
for (i in 1:SIZE){male[i]<-quantile(KM05[i,,],.5)}
malesdl<-array(0,c(SIZE))
for (i in 1:SIZE){malesdl[i]<-(quantile(KM05[i,,],.05))}
malesdu<-array(0,c(SIZE))
for (i in 1:SIZE){malesdu[i]<-(quantile(KM05[i,,],.95))}
female<-array(0,c(SIZE))
for (i in 1:SIZE){female[i]<-quantile(KF05[i,,],.5)}
femalesdl<-array(0,c(SIZE))
for (i in 1:SIZE){femalesdl[i]<-quantile(KF05[i,,],.05)}
femalesdu<-array(0,c(SIZE))
for (i in 1:SIZE){femalesdu[i]<-quantile(KF05[i,,],.95)}
cat<-seq(0,100,5)
poppyr3(male,female,malesdl,femalesdl,malesdu,femalesdu,cat)
mtext(side=3,line=1,text="ALASKA 2005 POPULATION ESTIMATES BY AGE AND SEX",font=2,cex=1)
mtext(side=3,line=0,text="",font=2,cex=1)
mtext(side=2,line=3,text="Age",font=1,cex=1)
mtext(side=1,line=3,text="Population",font=1,cex=1)
dev.copy2eps(file="C:/Users/Eddie/Desktop/Forecast_AK/Figure12.eps") 
Sys.sleep(5)

#FIGURE 13 (CARL MASON'S GREAT PYRAMID FUNCTION, WITH CI'S ADDED-- NEED GPLOTS PACKAGE)
library(gplots)
poppyr3<-function(male,female,malesdl,femalesdl,malesdu,femalesdu,cat){
split.screen(figs=rbind(c(0,.58,0,1),c(.43,1,0,1)))
screen(1)
barplot2(male,horiz=T,names=cat,plot.ci = TRUE,ci.l=malesdl,ci.u=malesdu,ci.color ="black",ci.lwd=1,space=0,
xlim=c(50000,0),col="dodgerblue")
title("Male",line=-3,cex.main=1)
screen(2)
barplot2(female,horiz=T,axisnames=F,plot.ci=TRUE,ci.l=femalesdl,ci.u=femalesdu,ci.color="black",ci.lwd=1,space=0,
xlim=c(0,50000),col="gold")
title("Female",line=-3,cex.main=1)
close.screen(all=T)}
male<-array(0,c(SIZE))
for (i in 1:SIZE){male[i]<-quantile(KM35[i,,],.5)}
malesdl<-array(0,c(SIZE))
for (i in 1:SIZE){malesdl[i]<-(quantile(KM35[i,,],.05))}
malesdu<-array(0,c(SIZE))
for (i in 1:SIZE){malesdu[i]<-(quantile(KM35[i,,],.95))}
female<-array(0,c(SIZE))
for (i in 1:SIZE){female[i]<-quantile(KF35[i,,],.5)}
femalesdl<-array(0,c(SIZE))
for (i in 1:SIZE){femalesdl[i]<-quantile(KF35[i,,],.05)}
femalesdu<-array(0,c(SIZE))
for (i in 1:SIZE){femalesdu[i]<-quantile(KF35[i,,],.95)}
cat<-seq(0,100,5)
poppyr3(male,female,malesdl,femalesdl,malesdu,femalesdu,cat)
mtext(side=3,line=1,text="ALASKA 2035 POPULATION FORECAST BY AGE AND SEX:",font=2,cex=1)
mtext(side=3,line=0,text="MEDIANS WITH 90% CONFIDENCE INTERVALS",font=2,cex=1)
mtext(side=2,line=3,text="Age",font=1,cex=1)
mtext(side=1,line=3,text="Population",font=1,cex=1)
dev.copy2eps(file="C:/Users/Eddie/Desktop/Forecast_AK/Figure13.eps")

#TABLE 1
quantile(KT05,c(.005,.05,.1,.25,.5,.75,.90,.95,.995))
quantile(KT15,c(.005,.05,.1,.25,.5,.75,.90,.95,.995))
quantile(KT25,c(.005,.05,.1,.25,.5,.75,.90,.95,.995))
quantile(KT35,c(.005,.05,.1,.25,.5,.75,.90,.95,.995))

#TABLE 2
quantile(GR0535*100,c(.005,.05,.1,.25,.5,.75,.90,.95,.995))

#TABLE 3
quantile(YA05,c(.005,.05,.1,.25,.5,.75,.90,.95,.995))
quantile(YA15,c(.005,.05,.1,.25,.5,.75,.90,.95,.995))
quantile(YA25,c(.005,.05,.1,.25,.5,.75,.90,.95,.995))
quantile(YA35,c(.005,.05,.1,.25,.5,.75,.90,.95,.995))

#TABLE 4
quantile(OA05,c(.005,.05,.1,.25,.5,.75,.90,.95,.995))
quantile(OA15,c(.005,.05,.1,.25,.5,.75,.90,.95,.995))
quantile(OA25,c(.005,.05,.1,.25,.5,.75,.90,.95,.995))
quantile(OA35,c(.005,.05,.1,.25,.5,.75,.90,.95,.995))

#TABLE *
quantile(TA05,c(.005,.05,.1,.25,.5,.75,.90,.95,.995))
quantile(TA15,c(.005,.05,.1,.25,.5,.75,.90,.95,.995))
quantile(TA25,c(.005,.05,.1,.25,.5,.75,.90,.95,.995))
quantile(TA35,c(.005,.05,.1,.25,.5,.75,.90,.95,.995))

#write.table(###, file="G:/###/###.csv", sep=",")

