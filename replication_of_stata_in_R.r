## https://github.com/EconometricsBySimulation/RStata/wiki/Dictionary:-Stata-to-R
## http://www.princeton.edu/~otorres/RStata.pdf
library(foreign)
library(data.table)
library(DataCombine)
library(Matching)
library(MatchIt)
library(AER)
library(quantreg)
library(rgenoud)
orig_data <- read.dta("/Users/gili/drive/minerva/2/cs112/assignments/proposal research/foriegn_aid/PerniciousEffectUNSC/USaid_UNSC_selectorate_data_replication.dta")

### DATA PRE PROCESS AND FEATURE EXTRACTION ###
# lagged aid
orig_data <- orig_data[order(orig_data$ccode, orig_data$year),]

## lagging func
L <- function(var_name,  n) {
  Data <- data.frame(orig_data$year, orig_data$ccode, orig_data[,var_name])
  colnames(Data) <- c("year", "ccode", var_name)
  Data <- Data[order(Data$ccode, Data$year),]
  
  DataSlid2 <- slide(Data, Var = var_name, GroupVar = "ccode",
                     slideBy = n, keepInvalid=TRUE )
  
  return(as.numeric(unlist(DataSlid2[4])))
}

#####################
#### first file #####
#####################
# HowMuchAidAndUNSCreplication.smcl

orig_data$Lln_totaid = L("ln_totaid96", -1)
orig_data$Ltotaid96 = L("totaid96", -1)

orig_data$year4SC = ifelse(orig_data$T0==0&orig_data$T1==0&orig_data$T2==0&orig_data$T3==0&orig_data$T4==0&orig_data$unmem==1, 0, NA)
orig_data$year4SC = ifelse( orig_data$T0==1&orig_data$unmem==1, 1, orig_data$year4SC)

orig_data$F2.T0 <- L("T0",2)
orig_data$F1.T0 <- L("T0",1)

orig_data$year2SC = ifelse(orig_data$T0==0&orig_data$T1==0&orig_data$T2==0&orig_data$T3==0&orig_data$T4==0&orig_data$unmem==1&orig_data$F1.T0==0&orig_data$F2.T0==0, 0, NA)
orig_data$year2SC = ifelse( orig_data$T0==1&orig_data$unmem==1, 1, orig_data$year2SC)

## T0 corresponds to year of election to SC 
##  compare average aid in years t-1 and t-2 with years t1 and t2  to find average increase in aid

##any aid as the inverse of no aid
orig_data$anyaid = 1 - orig_data$noaid

orig_data$L1anyaid=L("anyaid", -1) 

orig_data$anyaidWB= ifelse(!is.na(orig_data$AidGDP), 0, NA)
orig_data$anyaidWB= ifelse(orig_data$AidGDP>0, 1, orig_data$AidGDP)

orig_data$FanyaidWB=L("anyaidWB",1)
orig_data$LanyaidWB=L("anyaidWB",-1)

orig_data$LAidGDP=L("AidGDP",-1)

orig_data$noaidtminus2=L("noaid",-2)
orig_data$noaidtminus1=L("noaid",-1)
orig_data$noaidtplus1=L("noaid",1)
orig_data$noaidtplus2=L("noaid",2)

orig_data$noaidbefore= 0
orig_data$noaidbefore = ifelse(orig_data$noaidtminus1==1|orig_data$noaidtminus2==1|orig_data$noaid==1, 1, orig_data$noaidbefore)

orig_data$noaidduring= 0
orig_data$noaidduring = ifelse(orig_data$noaidtplus1==1|orig_data$noaidtplus2==1, 1, orig_data$noaidduring)

orig_data$aidbefore=(L("totaid96",-1)+L("totaid96",-2)) /2
orig_data$aidduring=(L("totaid96",1)+ L("totaid96",2) )/2
orig_data$aidafter=(L("totaid96",3)+L("totaid96",4)) /2                     
orig_data$F12lnTotAid96=log(orig_data$aidduring)
orig_data$lnaidduring=log(orig_data$aidduring)

orig_data$UNSC=orig_data$scmem
orig_data$UNSC= ifelse(orig_data$T0==1,1,orig_data$UNSC)
orig_data$UNSCimptyr = orig_data$imptyr * orig_data$UNSC

orig_data$endUNSC <- ifelse(orig_data$unmem==1, 0, NA)
orig_data$endUNSC <- ifelse(orig_data$T3==1|orig_data$T4==1, 1, orig_data$unmem)


## Does UNSC increase the likelihood of getting aid ##
table(orig_data[orig_data$hic==0,]$UNSC, orig_data[orig_data$hic==0,]$anyaid)
#        0    1
# 0 1951 4098
# 1   46  458

table(orig_data[(orig_data$hic==0& orig_data$L1anyaid==1),]$UNSC, orig_data[(orig_data$hic==0& orig_data$L1anyaid==1),]$anyaid)
#        0    1
# 0  135 3842
# 1    9  439

summary(glm( anyaid~ L1anyaid+demaut+scmem+lpop+lrgdpch+logdstab , data=orig_data[orig_data$hic ==0,],family=binomial(link = "probit")))

'''
Call:
glm(formula = anyaid ~ L1anyaid + demaut + scmem + lpop + lrgdpch + 
logdstab, family = binomial(link = "probit"), data = orig_data[orig_data$hic == 
0, ])

Deviance Residuals: 
Min       1Q   Median       3Q      Max  
-3.2059   0.1416   0.1811   0.2113   2.1962  

Coefficients:
Estimate Std. Error z value Pr(>|z|)    
(Intercept)  1.56764    0.81544   1.922 0.054549 .  
L1anyaid     2.92672    0.08041  36.396  < 2e-16 ***
demaut       0.44444    0.12303   3.613 0.000303 ***
scmem        0.11849    0.17930   0.661 0.508727    
lpop        -0.01144    0.02310  -0.495 0.620518    
lrgdpch     -0.16984    0.04637  -3.663 0.000250 ***
logdstab    -0.13843    0.07166  -1.932 0.053401 .  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for binomial family taken to be 1)

Null deviance: 3116.6  on 3956  degrees of freedom
Residual deviance: 1163.2  on 3950  degrees of freedom
(4403 observations deleted due to missingness)
AIC: 1177.2

Number of Fisher Scoring iterations: 7
'''
summary(glm( anyaid~ L1anyaid+demaut+UNSC+lpop+lrgdpch+logdstab , data=orig_data[orig_data$hic ==0,],family=binomial(link = "probit")))
'''
Call:
glm(formula = anyaid ~ L1anyaid + demaut + UNSC + lpop + lrgdpch + 
logdstab, family = binomial(link = "probit"), data = orig_data[orig_data$hic == 
0, ])

Deviance Residuals: 
Min       1Q   Median       3Q      Max  
-3.3197   0.1389   0.1816   0.2126   2.2010  

Coefficients:
Estimate Std. Error z value Pr(>|z|)    
(Intercept)  1.57735    0.81506   1.935 0.052959 .  
L1anyaid     2.92056    0.08031  36.365  < 2e-16 ***
demaut       0.44053    0.12310   3.579 0.000345 ***
UNSC         0.27690    0.15818   1.751 0.080017 .  
lpop        -0.01498    0.02319  -0.646 0.518312    
lrgdpch     -0.17217    0.04653  -3.700 0.000215 ***
logdstab    -0.13500    0.07149  -1.888 0.058973 .  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for binomial family taken to be 1)

Null deviance: 3116.6  on 3956  degrees of freedom
Residual deviance: 1160.0  on 3950  degrees of freedom
(4403 observations deleted due to missingness)
AIC: 1174

Number of Fisher Scoring iterations: 7
'''
summary(glm( anyaid~ L1anyaid+demaut+T0+T1+T2+T3+T4+lpop+lrgdpch+logdstab , data=orig_data[orig_data$hic ==0,],family=binomial(link = "probit")))
'''
Call:
glm(formula = anyaid ~ L1anyaid + demaut + T0 + T1 + T2 + T3 + 
T4 + lpop + lrgdpch + logdstab, family = binomial(link = "probit"), 
data = orig_data[orig_data$hic == 0, ])

Deviance Residuals: 
Min       1Q   Median       3Q      Max  
-3.3172   0.1326   0.1804   0.2151   2.2111  

Coefficients:
Estimate Std. Error z value Pr(>|z|)    
(Intercept)  1.546437   0.820199   1.885 0.059370 .  
L1anyaid     2.945800   0.082011  35.920  < 2e-16 ***
demaut       0.450218   0.123632   3.642 0.000271 ***
T0           0.625840   0.304937   2.052 0.040135 *  
T1          -0.176258   0.210415  -0.838 0.402217    
T2           0.573829   0.313727   1.829 0.067389 .  
T3          -0.184555   0.210681  -0.876 0.381033    
T4           0.003113   0.228069   0.014 0.989110    
lpop        -0.013951   0.023373  -0.597 0.550597    
lrgdpch     -0.171618   0.046762  -3.670 0.000243 ***
logdstab    -0.134288   0.072000  -1.865 0.062167 .  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for binomial family taken to be 1)

Null deviance: 3116.6  on 3956  degrees of freedom
Residual deviance: 1152.2  on 3946  degrees of freedom
(4403 observations deleted due to missingness)
AIC: 1174.2

Number of Fisher Scoring iterations: 7
'''

na.data <- cbind(
  as.numeric(orig_data$year), 
  as.numeric(orig_data$ccode), #country code
  as.numeric(orig_data$scmem), #SC member
  as.numeric(orig_data$UNSC), 
  as.numeric(orig_data$democ), #democracy score (1-10)
  as.numeric(orig_data$autoc), #autocracy score (1-10)
  as.numeric(orig_data$lpop), # log population 
  as.numeric(orig_data$GDPpcWB), # GDP per capita
  as.numeric(orig_data$AID), # total aid in $
  as.numeric(orig_data$demaut), #polity score betwwen 0-1
  as.numeric(orig_data$anyaid), #aid
  as.numeric(orig_data$lrgdpch),
  as.numeric(orig_data$logdstab),
  as.numeric(orig_data$L1anyaid),
  as.numeric(orig_data$hic),
  as.numeric(orig_data$T3),
  as.numeric(orig_data$endUNSC),
  as.numeric(orig_data$lpopWB),
  as.numeric(orig_data$lGDPpcWB)

)

na.data<- as.data.frame(na.data)

colnames(na.data) <- c("year", "ccode", "scmem", "UNSC","democ", "autoc", 
                       "lpop", "GDPpcWB", "AID",
                       "demaut","anyaid",
                       "lrgdpch","logdstab","L1anyaid","hic",
                       "T3", "endUNSC", "lpopWB", "lGDPpcWB")

# drop missing treatment vars
completeFunc <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}
na.data<-completeFunc(na.data, "UNSC")
na.data<-completeFunc(na.data, "T3")
na.data<-completeFunc(na.data, "endUNSC")
na.data<-completeFunc(na.data, "L1anyaid")


# Fill na with the mean for the rest 
# (linear regression shoudl not be impacted by the mean)
for(i in 1:ncol(na.data)){
  na.data[is.na(na.data[,i]), i] <- mean(na.data[,i], na.rm = TRUE)
}


## NN matching line 512
Tr=na.data[na.data$hic==0,]$UNSC
Y=na.data[na.data$hic==0,]$anyaid
X=cbind(na.data[na.data$hic==0,]$L1anyaid,na.data[na.data$hic==0,]$demaut,
        na.data[na.data$hic==0,]$lpop,na.data[na.data$hic==0,]$lrgdpch,
        na.data[na.data$hic==0,]$logdstab)


match.out <- Match(Y=Y, Tr=Tr, X=X, M=4, estimand='ATT',version="standard",replace=TRUE)

summary(match.out)
'''
Estimate...  0.012525 
AI SE......  0.011099 
T-stat.....  1.1284 
p.val......  0.25914 

Original number of observations..............  6550 
Original number of treated obs...............  507 
Matched number of observations...............  507 
Matched number of observations  (unweighted).  2039 

'''
## nn match line 526
Tr=na.data[(na.data$hic==0 & na.data$L1anyaid ==1),]$UNSC
Y=na.data[(na.data$hic==0 & na.data$L1anyaid ==1),]$anyaid
X=cbind(na.data[(na.data$hic==0 & na.data$L1anyaid ==1),]$lpop,na.data[(na.data$hic==0 & na.data$L1anyaid ==1),]$demaut,
        na.data[(na.data$hic==0 & na.data$L1anyaid ==1),]$lrgdpch,na.data[(na.data$hic==0 & na.data$L1anyaid ==1),]$logdstab)


match.out <- Match(Y=Y, Tr=Tr, X=X, M=2, estimand='ATT',version="standard",replace=TRUE)

summary(match.out)
'''


Estimate...  0.0062685 
AI SE......  0.010458 
T-stat.....  0.5994 
p.val......  0.54891 

Original number of observations..............  4554 
Original number of treated obs...............  457 
Matched number of observations...............  457 
Matched number of observations  (unweighted).  920 


'''

## nn match line 541

Tr=na.data[(na.data$hic==0 & na.data$L1anyaid ==1),]$T3
Y=na.data[(na.data$hic==0 & na.data$L1anyaid ==1),]$anyaid
X=cbind(na.data[(na.data$hic==0 & na.data$L1anyaid ==1),]$lpop,na.data[(na.data$hic==0 & na.data$L1anyaid ==1),]$demaut,
        na.data[(na.data$hic==0 & na.data$L1anyaid ==1),]$lrgdpch,na.data[(na.data$hic==0 & na.data$L1anyaid ==1),]$logdstab)


match.out <- Match(Y=Y, Tr=Tr, X=X, M=2, estimand='ATT',version="standard",replace=TRUE)

summary(match.out)
'''

Estimate...  -0.025475 
AI SE......  0.017568 
T-stat.....  -1.4501 
p.val......  0.14704 

Original number of observations..............  4554 
Original number of treated obs...............  155 
Matched number of observations...............  155 
Matched number of observations  (unweighted).  316 
'''

## nn match line 564 X

Tr=na.data[(na.data$hic==0 & na.data$L1anyaid ==1),]$endUNSC
Y=na.data[(na.data$hic==0 & na.data$L1anyaid ==1),]$anyaid
X=cbind(na.data[(na.data$hic==0 & na.data$L1anyaid ==1),]$lpop,na.data[(na.data$hic==0 & na.data$L1anyaid ==1),]$demaut,
        na.data[(na.data$hic==0 & na.data$L1anyaid ==1),]$lrgdpch,na.data[(na.data$hic==0 & na.data$L1anyaid ==1),]$logdstab)

match.out <- Match(Y=Y, Tr=Tr,M=4, X=X,estimand='ATT', replace= TRUE)

summary(match.out)
'''

Estimate...  0.16525 
AI SE......  0.042963 
T-stat.....  3.8463 
p.val......  0.00011991 

Original number of observations..............  4554 
Original number of treated obs...............  4410 
Matched number of observations...............  4410 
Matched number of observations  (unweighted).  17856 

'''

## NN match line 581 X

Tr=na.data[(na.data$hic==0),]$endUNSC
Y=na.data[(na.data$hic==0),]$anyaid
X=cbind(na.data[(na.data$hic==0),]$lpop,na.data[(na.data$hic==0),]$demaut,
        na.data[(na.data$hic==0),]$lrgdpch,na.data[(na.data$hic==0),]$logdstab,
        na.data[(na.data$hic==0),]$L1anyaid)

C <- cov(X)
match.out <- Match(Y=Y, Tr=Tr, X=X,Weight.matrix=C,  M=4, estimand='ATT', replace= TRUE)
summary(match.out)
'''
Estimate...  0.18793 
AI SE......  0.034126 
T-stat.....  5.5069 
p.val......  3.6529e-08 

Original number of observations..............  6550 
Original number of treated obs...............  5384 
Matched number of observations...............  5384 
Matched number of observations  (unweighted).  21621  
'''
############################

#####################
#### second file ####
#####################
#USaidUNSCreplication.smcl

orig_data$L4demaut=L("demaut",-4)
orig_data$delta4demaut=L("demaut",4)-orig_data$demaut
orig_data$L4delta4demaut=L("delta4demaut",4)

### TABLE 2
## line 20
summary(glm( T0~ demaut+lpopWB+lGDPpcWB+tau_lead, 
             data=orig_data[(orig_data$unmem==1 & orig_data$P5 !=1 & !(orig_data$T1==1|orig_data$T2==1|orig_data$T3==1|orig_data$T4==1)),],
             family=binomial(link = "probit")))
'''
Call:
glm(formula = T0 ~ demaut + lpopWB + lGDPpcWB + tau_lead, family = binomial(link = "probit"), 
data = orig_data[(orig_data$unmem == 1 & orig_data$P5 != 
1 & !(orig_data$T1 == 1 | orig_data$T2 == 1 | orig_data$T3 == 
1 | orig_data$T4 == 1)), ])

Deviance Residuals: 
Min       1Q   Median       3Q      Max  
-0.6329  -0.3515  -0.2745  -0.2195   2.9577  

Coefficients:
Estimate Std. Error z value Pr(>|z|)    
(Intercept) -5.190181   0.505660 -10.264  < 2e-16 ***
demaut       0.034031   0.114880   0.296 0.767053    
lpopWB       0.169346   0.027281   6.208 5.38e-10 ***
lGDPpcWB     0.107236   0.028696   3.737 0.000186 ***
tau_lead     0.006408   0.117932   0.054 0.956670    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for binomial family taken to be 1)

Null deviance: 1356.4  on 3585  degrees of freedom
Residual deviance: 1295.1  on 3581  degrees of freedom
(3728 observations deleted due to missingness)
AIC: 1305.1

Number of Fisher Scoring iterations: 6
'''

## line 42
summary(glm( T0~ demaut+lpopWB+lGDPpcWB+tau_lead+AidGDP, 
             data=orig_data[(orig_data$unmem==1 & orig_data$P5 !=1 & !(orig_data$T1==1|orig_data$T2==1|orig_data$T3==1|orig_data$T4==1)),],
             family=binomial(link = "probit")))
'''
Call:
glm(formula = T0 ~ demaut + lpopWB + lGDPpcWB + tau_lead + AidGDP, 
family = binomial(link = "probit"), data = orig_data[(orig_data$unmem == 
1 & orig_data$P5 != 1 & !(orig_data$T1 == 1 | orig_data$T2 == 
1 | orig_data$T3 == 1 | orig_data$T4 == 1)), ])

Deviance Residuals: 
Min       1Q   Median       3Q      Max  
-0.5567  -0.3149  -0.2668  -0.2212   3.1794  

Coefficients:
Estimate Std. Error z value Pr(>|z|)    
(Intercept) -4.345768   0.702120  -6.189 6.04e-10 ***
demaut      -0.022428   0.129492  -0.173    0.862    
lpopWB       0.142271   0.031999   4.446 8.74e-06 ***
lGDPpcWB     0.058319   0.043576   1.338    0.181    
tau_lead     0.028365   0.130801   0.217    0.828    
AidGDP      -0.008990   0.008208  -1.095    0.273    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for binomial family taken to be 1)

Null deviance: 1019.59  on 2990  degrees of freedom
Residual deviance:  984.94  on 2985  degrees of freedom
(4323 observations deleted due to missingness)
AIC: 996.94

Number of Fisher Scoring iterations: 6
'''

## line 88
##### MISSING TAU

orig_data$ln_totaid = log1p(orig_data$totaid)
orig_data$Lln_totaid = L("ln_totaid",-1)

## line 110
summary(glm( T0~ demaut+lpop+lrgdpch+tau_lead+Lln_totaid, 
             data=orig_data[(orig_data$unmem==1 & orig_data$P5 !=1 & !(orig_data$T1==1|orig_data$T2==1|orig_data$T3==1|orig_data$T4==1)),],
             family=binomial(link = "probit")))
'''
Call:
glm(formula = T0 ~ demaut + lpop + lrgdpch + tau_lead + Lln_totaid, 
family = binomial(link = "probit"), data = orig_data[(orig_data$unmem == 
1 & orig_data$P5 != 1 & !(orig_data$T1 == 1 | orig_data$T2 == 
1 | orig_data$T3 == 1 | orig_data$T4 == 1)), ])

Deviance Residuals: 
Min       1Q   Median       3Q      Max  
-0.6463  -0.3449  -0.2799  -0.2268   3.0710  

Coefficients:
Estimate Std. Error z value Pr(>|z|)    
(Intercept) -3.91639    0.40294  -9.720  < 2e-16 ***
demaut       0.20395    0.10966   1.860   0.0629 .  
lpop         0.18467    0.02665   6.929 4.23e-12 ***
lrgdpch      0.06688    0.04157   1.609   0.1077    
tau_lead     0.08987    0.11085   0.811   0.4175    
Lln_totaid  -0.05445    0.02256  -2.413   0.0158 *  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for binomial family taken to be 1)

Null deviance: 1553.7  on 4010  degrees of freedom
Residual deviance: 1480.8  on 4005  degrees of freedom
(3303 observations deleted due to missingness)
AIC: 1492.8

Number of Fisher Scoring iterations: 6
'''

####### Effects of Membership of Security Council ######

orig_data$F1.T0<- L("T0", 1)
orig_data$F2.T0<- L("T0", 2)
orig_data$F3.T0<- L("T0", 3)

orig_data$year4SC = NA
orig_data$year4SC = ifelse(orig_data$T0==0&orig_data$T1==0&orig_data$T2==0&orig_data$T3==0&orig_data$T4==0&
                             orig_data$F1.T0==0&orig_data$F2.T0==0&orig_data$F3.T0==0&orig_data$unmem==1,0,orig_data$year4SC)

orig_data$year4SC = ifelse(orig_data$T0==1&orig_data$unmem==1 ,1, orig_data$year4SC )


orig_data$year2SC = NA
orig_data$year2SC= ifelse(orig_data$T0==0&orig_data$T1==0&orig_data$T2==0&orig_data$T3==0&orig_data$T4==0&orig_data$F1.T0==0&orig_data$F2.T0==0&orig_data$unmem==1, 0, orig_data$year2SC)
orig_data$year2SC= ifelse(orig_data$T0==1&orig_data$unmem==1, 1,orig_data$year2SC)

orig_data$L.totaid96 <- L("totaid96",-1)

orig_data$increaseAID = ifelse (orig_data$unmem==1, 0, NA)
orig_data$increaseAID = ifelse ((((orig_data$totaid96-orig_data$L.totaid96)/orig_data$L.totaid96 >.5) & orig_data$Lanyaid==1) ,1 ,orig_data$increaseAID )

orig_data$F1.totaid96 <- L("totaid96", 1)
orig_data$F2.totaid96 <- L("totaid96", 2)
orig_data$F.anyaid <- L("anyaid", 1)
orig_data$F2.anyaid<- L("anyaid", 2)

#/* 50% increase in aid in year */
orig_data$increaseAID= ifelse(((orig_data$F1.totaid96-orig_data$L.totaid96)/orig_data$L.totaid96>.5)& orig_data$Lanyaid==1, 1, orig_data$increaseAID)

#/* 50% increase in aid in next year */
orig_data$increaseAID= ifelse(((orig_data$F2.totaid96-orig_data$L.totaid96)/orig_data$L.totaid96>.5)& orig_data$Lanyaid==1, 1, orig_data$increaseAID)

# /* 50% increase in aid in year t+2 */
orig_data$increaseAID=ifelse(orig_data$anyaid==1&  orig_data$Lanyaid==0, 1, orig_data$increaseAID)

# /*start aid in year t*/ 
orig_data$increaseAID=ifelse(orig_data$F.anyaid==1&  orig_data$Lanyaid==0, 1, orig_data$increaseAID) 

# /* start aid in t+1 */
orig_data$increaseAID=ifelse( orig_data$F2.anyaid==1&  orig_data$Lanyaid==0, 1,orig_data$increaseAID) 

# /* start aid in t+2 */ 
orig_data$UNSCmore= ifelse(orig_data$unmem==1, 0, NA )
orig_data$UNSCmore=ifelse(orig_data$increaseAID==1 & orig_data$year4SC==1, 1, orig_data$UNSCmore )

###/******** Growth in Aid ********/###
orig_data$F3.totaid96<- L("totaid96", 3)
orig_data$F4.totaid96<- L("totaid96", 4)

orig_data$USaidgrowth=(orig_data$F2.totaid96+orig_data$F1.totaid96)/(2*orig_data$totaid96)
orig_data$USaidgrowth34=(orig_data$F3.totaid96+orig_data$F4.totaid96)/(2*orig_data$totaid96)
orig_data$LnUSaid12=log1p((orig_data$F2.totaid96+orig_data$F1.totaid96)/2)

orig_data$F1.AID<- L("AID", 1)
orig_data$F2.AID<- L("AID", 2)
orig_data$F3.AID<- L("AID", 3)
orig_data$F4.AID<- L("AID", 4)
orig_data$F1.AidGDP<- L("AidGDP", 1)
orig_data$F2.AidGDP<- L("AidGDP", 2)
orig_data$F1.OIL<- L("OIL", 1)
orig_data$F2.OIL<- L("OIL", 2)

orig_data$AID12=log1p((orig_data$F2.AID+orig_data$F1.AID)/2)
orig_data$AIDGDP12=(orig_data$F2.AidGDP+orig_data$F1.AidGDP)/2
orig_data$OILGDP12=(orig_data$F2.OIL+orig_data$F1.OIL)/2

orig_data$AIDgrowth=(orig_data$F2.AID+orig_data$F1.AID)/(2*orig_data$AID)
orig_data$AIDgrowth34=(orig_data$F3.AID+orig_data$F4.AID)/(2*orig_data$AID)

orig_data$demautUSaidgrowth=orig_data$demaut*orig_data$USaidgrowth
orig_data$demautUSaidgrowth34=orig_data$demaut*orig_data$USaidgrowth34
orig_data$demautLnUSaid12=orig_data$demaut*orig_data$LnUSaid12
orig_data$demautAID12=orig_data$demaut*orig_data$AID12
orig_data$demautAIDGDP12=orig_data$demaut*orig_data$AIDGDP12
orig_data$demautOILGDP12=orig_data$demaut*orig_data$OILGDP12
orig_data$demautAIDgrowth=orig_data$demaut*orig_data$AIDgrowth
orig_data$demautAIDgrowth34=orig_data$demaut*orig_data$AIDgrowth34
orig_data$demautOIL=orig_data$demaut*orig_data$OIL
orig_data$demautAidGDP=orig_data$demaut*orig_data$AidGDP

### /********** military expenditure ********/ ###
orig_data$lnmilex=log1p(orig_data$milex)
orig_data$lnmilper=log1p(orig_data$milper)

orig_data$Llnmilper=L("lnmilper",-1)
orig_data$L2lnmilper=L("lnmilper",-2)

orig_data$Llnmilex=L("lnmilex",-1)
orig_data$L2lnmilex=L("lnmilex",-2)

orig_data$F4.milex=L("lnmilex",4)
orig_data$F2.milex = L("lnmilex",2)

orig_data$F4.milper=L("lnmilper",4)
orig_data$F2.milper = L("lnmilper",2)

orig_data$delta4milex=100*(orig_data$F4.milex-orig_data$milex)/orig_data$milex
orig_data$delta2milex=100*(orig_data$F2.milex-orig_data$milex)/orig_data$milex

orig_data$delta4milper=100*(orig_data$F4.milper-orig_data$milper)/orig_data$milper
orig_data$delta2milper=100*(orig_data$F2.milper-orig_data$milper)/orig_data$milper

### /******** leader survival **************/ ###
orig_data$leadchange=ifelse(orig_data$unmem==1, orig_data$LEADERCHANGE, NA)

orig_data$F1.LEADERCHANGE=L("LEADERCHANGE",1)
orig_data$F2.LEADERCHANGE=L("LEADERCHANGE",2)
orig_data$F3.LEADERCHANGE=L("LEADERCHANGE",3)
orig_data$F4.LEADERCHANGE=L("LEADERCHANGE",4)

orig_data$leadchange=ifelse(orig_data$unmem==1 & orig_data$T0==1 , orig_data$F1.LEADERCHANGE+orig_data$F2.LEADERCHANGE/2 , orig_data$leadchange)
orig_data$leadchange1234=ifelse(orig_data$unmem==1, orig_data$LEADERCHANGE, NA)
orig_data$leadchange1234=ifelse(orig_data$unmem==1 & orig_data$T0==1, (orig_data$F1.LEADERCHANGE+orig_data$F2.LEADERCHANGE+orig_data$F3.LEADERCHANGE+orig_data$F4.LEADERCHANGE)/4, orig_data$leadchange1234)

### /*************** GDP *******************/ ###

orig_data$F4.rgdpl <- L("rgdpl",4)
orig_data$F2.rgdpl <- L("rgdpl",2)
orig_data$F4.GDPpcWB <- L("GDPpcWB",4)
orig_data$F2.GDPpcWB <- L("GDPpcWB",2)
orig_data$F4.rgdpch <- L("rgdpch",4)
orig_data$F2.rgdpch <- L("rgdpch",2)

orig_data$F4.demaut<-  L("demaut",4)
  
orig_data$delta4rgdpl=100*(orig_data$F4.rgdpl-orig_data$rgdpl)/orig_data$rgdpl
orig_data$delta4GDPpcWB=100*(orig_data$F4.GDPpcWB- orig_data$GDPpcWB)/ orig_data$GDPpcWB
orig_data$delta4rgdpch=100*(orig_data$F4.rgdpch-orig_data$rgdpch)/orig_data$rgdpch

orig_data$delta2rgdpl=100*(orig_data$F2.rgdpl-orig_data$rgdpl)/orig_data$rgdpl
orig_data$delta2GDPpcWB=100*(orig_data$F2.GDPpcWB- orig_data$GDPpcWB)/ orig_data$GDPpcWB
orig_data$delta2rgdpch=100*(orig_data$F2.rgdpch-orig_data$rgdpch)/orig_data$rgdpch

orig_data$delta4demaut=orig_data$F4.demaut-orig_data$demaut
orig_data$L4delta4demaut= L("delta4demaut",-4)
orig_data$delta2demaut=L("demaut", 2)-orig_data$demaut
orig_data$F4demaut=L("demaut",4)
orig_data$F2demaut=L("demaut",2)
orig_data$Ldemaut=L("demaut",-1)

orig_data$demautyear4SC=orig_data$demaut*orig_data$year4SC
orig_data$demautyear2SC=orig_data$demaut*orig_data$year2SC

orig_data$delta4kg=L("kg",4)-orig_data$kg
orig_data$delta2kg=L("kg",2)-orig_data$kg

orig_data$delta4tau=L("tau_lead",4)-orig_data$tau_lead
orig_data$delta2tau=L("tau_lead",2)-orig_data$tau_lead

orig_data$F4tau=L("tau_lead",4)
orig_data$F2tau=L("tau_lead",2)
orig_data$Ltau=L("tau_lead",-1)

orig_data$delta4aid=100*(orig_data$F4.totaid96-orig_data$totaid96)/orig_data$totaid96
orig_data$delta2aid=100*(orig_data$F2.totaid96-orig_data$totaid96)/orig_data$totaid96

orig_data$delta4score=L("score",4)-orig_data$score
orig_data$delta2score=L("score",2)-orig_data$score

orig_data$F4score=L("score",4)
orig_data$F2score=L("score",2)
orig_data$Lscore=L("score",-1)

orig_data$delta2mass=L("mass", 2)-orig_data$mass
orig_data$delta4mass=L("mass", 4)-orig_data$mass

orig_data$L4lpopWB= L("lpopWB",-4)
orig_data$L4lGDPpcWB=L("lGDPpcWB",-4)

### /********************** GDP ***********/ ###

orig_data$L4delta4GDPpcWB=L("delta4GDPpcWB",-4)
orig_data$L4delta4rgdpch=L("delta4rgdpch",-4)



## TABLE 3 ##
########### delta 2 years ############
na.data <- cbind(
  as.numeric(orig_data$year), 
  as.numeric(orig_data$ccode), #country code
  as.numeric(orig_data$scmem), #SC member
  as.numeric(orig_data$UNSC), 
  as.numeric(orig_data$democ), #democracy score (1-10)
  as.numeric(orig_data$autoc), #autocracy score (1-10)
  as.numeric(orig_data$lpop), # log population 
  as.numeric(orig_data$GDPpcWB), # GDP per capita
  as.numeric(orig_data$AID), # total aid in $
  as.numeric(orig_data$demaut), #polity score betwwen 0-1
  as.numeric(orig_data$anyaid), #aid
  as.numeric(orig_data$lrgdpch),
  as.numeric(orig_data$logdstab),
  as.numeric(orig_data$L1anyaid),
  as.numeric(orig_data$hic),
  as.numeric(orig_data$T3),
  as.numeric(orig_data$endUNSC),
  as.numeric(orig_data$lpopWB),
  as.numeric(orig_data$lGDPpcWB),
  as.numeric(orig_data$year4SC),
  as.numeric(orig_data$delta4GDPpcWB),
  as.numeric(orig_data$unmem),
  as.numeric(orig_data$delta4demaut),
  as.numeric(orig_data$delta4score),
  as.numeric(orig_data$delta4tau),
  as.numeric(orig_data$delta2GDPpcWB),
  as.numeric(orig_data$delta2demaut),
  as.numeric(orig_data$delta2score),
  as.numeric(orig_data$delta2tau),
  as.numeric(orig_data$year2SC)

)

na.data<- as.data.frame(na.data)

colnames(na.data) <- c("year", "ccode", "scmem", "UNSC","democ", "autoc", 
                       "lpop", "GDPpcWB", "AID",
                       "demaut","anyaid",
                       "lrgdpch","logdstab","L1anyaid","hic",
                       "T3", "endUNSC", "lpopWB", "lGDPpcWB",
                       "year4SC","delta4GDPpcWB","unmem",
                       "delta4demaut", "delta4score", "delta4tau",
                       "delta2GDPpcWB",
                       "delta2demaut", "delta2score", "delta2tau",
                       "year2SC")


# drop missing 

# drop missing treatment vars
completeFunc <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}
na.data<-completeFunc(na.data, "year4SC")
#na.data<-completeFunc(na.data, "year2SC")


# Fill na with the mean for the rest 
# (linear regression shoudl not be impacted by the mean)
for(i in 1:ncol(na.data)){
  na.data[is.na(na.data[,i]), i] <- mean(na.data[,i], na.rm = TRUE)
}


## NN matching line 1847 (DELTA GDP)

Tr=na.data[na.data$unmem==1,]$year4SC
Y=na.data[na.data$unmem==1,]$delta4GDPpcWB

X=cbind(na.data[na.data$unmem==1,]$demaut,
        na.data[na.data$unmem==1,]$lpopWB,
        na.data[na.data$unmem==1,]$lGDPpcWB)

C <- cov(X)

GDP4.match.out <- Match(Y=Y, Tr=Tr, X=X,Weight.matrix=C, M=45, estimand='ATT', replace = TRUE)
summary(GDP4.match.out)

'''
Estimate...  -2.3803 
AI SE......  0.88478 
T-stat.....  -2.6903 
p.val......  0.0071393 

Original number of observations..............  6728 
Original number of treated obs...............  253 
Matched number of observations...............  253 
Matched number of observations  (unweighted).  13575 
'''
GDP4.mb <- MatchBalance(year4SC ~ lGDPpcWB+lpop+demaut,data=na.data, match.out=GDP4.match.out, nboots=500)

'''
Before Matching Minimum p.value: < 2.22e-16 
Variable Name(s): lpop demaut  Number(s): 2 3 

After Matching Minimum p.value: < 2.22e-16 
Variable Name(s): lpop demaut  Number(s): 2 3 
'''
## NN matching ? (DELTA DEMOC)
Tr=na.data[na.data$unmem==1,]$year4SC
Y=na.data[na.data$unmem==1,]$delta4demaut

X=cbind(
        na.data[na.data$unmem==1,]$lpopWB,
        na.data[na.data$unmem==1,]$lGDPpcWB)

C <- cov(X)

DEMOC4.match.out <- Match(Y=Y, Tr=Tr, X=X,Weight.matrix=C, M=50, estimand='ATT', replace = TRUE)
summary(DEMOC4.match.out)
'''

Estimate...  -0.014535 
AI SE......  0.010966 
T-stat.....  -1.3255 
p.val......  0.18502 

Original number of observations..............  6728 
Original number of treated obs...............  253 
Matched number of observations...............  253 
Matched number of observations  (unweighted).  43569 
'''
DEMOC4.mb <- MatchBalance(year4SC ~ lGDPpcWB+lpop+demaut,data=na.data, match.out=DEMOC4.match.out, nboots=500)
'''
Before Matching Minimum p.value: < 2.22e-16 
Variable Name(s): lpop demaut  Number(s): 2 3 

After Matching Minimum p.value: < 2.22e-16 
Variable Name(s): lpop demaut  Number(s): 2 3
'''

## NN matching  (DELTA PRESS) X
Tr=na.data[na.data$unmem==1,]$year4SC
Y=na.data[na.data$unmem==1,]$delta4score

X=cbind(
  na.data[na.data$unmem==1,]$demaut,
  na.data[na.data$unmem==1,]$lpopWB,
  na.data[na.data$unmem==1,]$lGDPpcWB)

C <- cov(X)

PRESS4.match.out <- Match(Y=Y, Tr=Tr, X=X,Weight.matrix=C, M=5, estimand='ATT', replace = TRUE)
summary(PRESS4.match.out)
'''
Estimate...  0.22588 
AI SE......  0.18992 
T-stat.....  1.1893 
p.val......  0.23431 

Original number of observations..............  6728 
Original number of treated obs...............  253 
Matched number of observations...............  253 
Matched number of observations  (unweighted).  4212 
'''

PRESS4.mb <- MatchBalance(year4SC ~ lGDPpcWB+lpop+demaut,data=na.data, match.out=PRESS4.match.out, nboots=500)
'''
Before Matching Minimum p.value: < 2.22e-16 
Variable Name(s): lpop demaut  Number(s): 2 3 

After Matching Minimum p.value: < 2.22e-16 
Variable Name(s): lpop  Number(s): 2 
'''

## NN matching ? (DELTA US ALIGN)
Tr=na.data[na.data$unmem==1,]$year4SC
Y=na.data[na.data$unmem==1,]$delta4tau

X=cbind(
  na.data[na.data$unmem==1,]$demaut,
  na.data[na.data$unmem==1,]$lpop,
  na.data[na.data$unmem==1,]$lGDPpcWB)

C <- cov(X)

USALIGN4.match.out <- Match(Y=Y, Tr=Tr, X=X,Weight.matrix=C, M=6, estimand='ATT', replace = TRUE)
summary(USALIGN4.match.out)
'''
Estimate...  0.0066259 
AI SE......  0.0071656 
T-stat.....  0.92468 
p.val......  0.35513 

Original number of observations..............  6728 
Original number of treated obs...............  253 
Matched number of observations...............  253 
Matched number of observations  (unweighted).  1924 
'''
USALIGN4.mb <- MatchBalance(year4SC ~ lGDPpcWB+lpop+demaut,data=na.data, match.out=USALIGN4.match.out, nboots=500)

'''
Before Matching Minimum p.value: < 2.22e-16 
Variable Name(s): lpop demaut  Number(s): 2 3 

After Matching Minimum p.value: 0.00086155 
Variable Name(s): lpop  Number(s): 2 

'''

########### delta 2 years ############

na.data <- cbind(
  as.numeric(orig_data$year), 
  as.numeric(orig_data$ccode), #country code
  as.numeric(orig_data$scmem), #SC member
  as.numeric(orig_data$UNSC), 
  as.numeric(orig_data$democ), #democracy score (1-10)
  as.numeric(orig_data$autoc), #autocracy score (1-10)
  as.numeric(orig_data$lpop), # log population 
  as.numeric(orig_data$GDPpcWB), # GDP per capita
  as.numeric(orig_data$AID), # total aid in $
  as.numeric(orig_data$demaut), #polity score betwwen 0-1
  as.numeric(orig_data$anyaid), #aid
  as.numeric(orig_data$lrgdpch),
  as.numeric(orig_data$logdstab),
  as.numeric(orig_data$L1anyaid),
  as.numeric(orig_data$hic),
  as.numeric(orig_data$T3),
  as.numeric(orig_data$endUNSC),
  as.numeric(orig_data$lpopWB),
  as.numeric(orig_data$lGDPpcWB),
  as.numeric(orig_data$year4SC),
  as.numeric(orig_data$delta4GDPpcWB),
  as.numeric(orig_data$unmem),
  as.numeric(orig_data$delta4demaut),
  as.numeric(orig_data$delta4score),
  as.numeric(orig_data$delta4tau),
  as.numeric(orig_data$delta2GDPpcWB),
  as.numeric(orig_data$delta2demaut),
  as.numeric(orig_data$delta2score),
  as.numeric(orig_data$delta2tau),
  as.numeric(orig_data$year2SC)
  
)

na.data<- as.data.frame(na.data)

colnames(na.data) <- c("year", "ccode", "scmem", "UNSC","democ", "autoc", 
                       "lpop", "GDPpcWB", "AID",
                       "demaut","anyaid",
                       "lrgdpch","logdstab","L1anyaid","hic",
                       "T3", "endUNSC", "lpopWB", "lGDPpcWB",
                       "year4SC","delta4GDPpcWB","unmem",
                       "delta4demaut", "delta4score", "delta4tau",
                       "delta2GDPpcWB",
                       "delta2demaut", "delta2score", "delta2tau",
                       "year2SC")


# drop missing 

# drop missing treatment vars
completeFunc <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}
#na.data<-completeFunc(na.data, "year4SC")
na.data<-completeFunc(na.data, "year2SC")


# Fill na with the mean for the rest 
# (linear regression shoudl not be impacted by the mean)
for(i in 1:ncol(na.data)){
  na.data[is.na(na.data[,i]), i] <- mean(na.data[,i], na.rm = TRUE)
}


## NN matching  (DELTA GDP)

Tr=na.data[na.data$unmem==1,]$year2SC
Y=na.data[na.data$unmem==1,]$delta2GDPpcWB

X=cbind(na.data[na.data$unmem==1,]$demaut,
        na.data[na.data$unmem==1,]$lpopWB,
        na.data[na.data$unmem==1,]$lGDPpcWB)

C <- cov(X)

GDP2.match.out <- Match(Y=Y, Tr=Tr, X=X,Weight.matrix=C, M=3, estimand='ATT', replace = TRUE)
summary(GDP2.match.out)
'''
Estimate...  -0.55012 
AI SE......  0.69907 
T-stat.....  -0.78693 
p.val......  0.43132 

Original number of observations..............  5774 
Original number of treated obs...............  253 
Matched number of observations...............  253 
Matched number of observations  (unweighted).  3215
'''
GDP2.mb <- MatchBalance(year2SC ~ lGDPpcWB+lpop+demaut,data=na.data, match.out=GDP2.match.out, nboots=500)
'''
Before Matching Minimum p.value: < 2.22e-16 
Variable Name(s): lpop demaut  Number(s): 2 3 

After Matching Minimum p.value: < 2.22e-16 
Variable Name(s): lpop  Number(s): 2   
'''
## NN matching ? (DELTA DEMOC)
Tr=na.data[na.data$unmem==1,]$year2SC
Y=na.data[na.data$unmem==1,]$delta2demaut

X=cbind(
  na.data[na.data$unmem==1,]$lpopWB,
  na.data[na.data$unmem==1,]$lGDPpcWB)

C <- cov(X)

DEMOC2.match.out <- Match(Y=Y, Tr=Tr, X=X,Weight.matrix=C, M=1, estimand='ATT', replace = TRUE)
summary(DEMOC2.match.out)
'''
Estimate...  -0.010419 
AI SE......  0.01057 
T-stat.....  -0.98572 
p.val......  0.32427 

Original number of observations..............  5774 
Original number of treated obs...............  253 
Matched number of observations...............  253 
Matched number of observations  (unweighted).  28640 
'''

DEMOC2.mb <- MatchBalance(year2SC ~ lGDPpcWB+lpop+demaut,data=na.data, match.out=DEMOC2.match.out, nboots=500)
'''
Before Matching Minimum p.value: < 2.22e-16 
Variable Name(s): lpop demaut  Number(s): 2 3 

After Matching Minimum p.value: < 2.22e-16 
Variable Name(s): lpop demaut  Number(s): 2 3
'''

## NN matching ? (DELTA PRESS) X
Tr=na.data[na.data$unmem==1,]$year2SC
Y=na.data[na.data$unmem==1,]$delta2score

X=cbind(
  na.data[na.data$unmem==1,]$demaut,
  na.data[na.data$unmem==1,]$lpopWB,
  na.data[na.data$unmem==1,]$lGDPpcWB)

C <- cov(X)

PRESS2.match.out <- Match(Y=Y, Tr=Tr, X=X,Weight.matrix=C, M=1, estimand='ATT', replace = TRUE)
summary(PRESS2.match.out)
'''
Estimate...  0.41079 
AI SE......  0.28217 
T-stat.....  1.4558 
p.val......  0.14544 

Original number of observations..............  5774 
Original number of treated obs...............  253 
Matched number of observations...............  253 
Matched number of observations  (unweighted).  2798  
'''

PRESS2.mb <- MatchBalance(year2SC ~ lGDPpcWB+lpop+demaut,data=na.data, match.out=PRESS2.match.out, nboots=500)
'''
Before Matching Minimum p.value: < 2.22e-16 
Variable Name(s): lpop demaut  Number(s): 2 3 

After Matching Minimum p.value: < 2.22e-16 
Variable Name(s): lpop  Number(s): 2
'''

## NN matching ? (DELTA US ALIGN)
Tr=na.data[na.data$unmem==1,]$year2SC
Y=na.data[na.data$unmem==1,]$delta2tau

X=cbind(
  na.data[na.data$unmem==1,]$demaut,
  na.data[na.data$unmem==1,]$lpop,
  na.data[na.data$unmem==1,]$lGDPpcWB)

C <- cov(X)

USALIGN2.match.out <- Match(Y=Y, Tr=Tr, X=X,Weight.matrix=C, M=5, estimand='ATT', replace = TRUE)
summary(USALIGN2.match.out)
'''
Estimate...  0.008815 
AI SE......  0.008091 
T-stat.....  1.0895 
p.val......  0.27594 

Original number of observations..............  5774 
Original number of treated obs...............  253 
Matched number of observations...............  253 
Matched number of observations  (unweighted).  1560 
'''
USALIGN2.mb <- MatchBalance(year2SC ~ lGDPpcWB+lpop+demaut,data=na.data, match.out=USALIGN2.match.out, nboots=500)

'''
Before Matching Minimum p.value: < 2.22e-16 
Variable Name(s): lpop demaut  Number(s): 2 3 

After Matching Minimum p.value: 0.0038154 
Variable Name(s): lGDPpcWB  Number(s): 1
'''

########################################
########### PERSOANL CODE ##############
########################################

'''
PLAN
1.GenMatch
2. IV reg
3. QR
'''

### GenMatch ###

########### delta 4 years ############

na.data <- cbind(
  as.numeric(orig_data$year), 
  as.numeric(orig_data$ccode), #country code
  as.numeric(orig_data$scmem), #SC member
  as.numeric(orig_data$UNSC), 
  as.numeric(orig_data$democ), #democracy score (1-10)
  as.numeric(orig_data$autoc), #autocracy score (1-10)
  as.numeric(orig_data$lpop), # log population 
  as.numeric(orig_data$GDPpcWB), # GDP per capita
  as.numeric(orig_data$AID), # total aid in $
  as.numeric(orig_data$demaut), #polity score betwwen 0-1
  as.numeric(orig_data$anyaid), #aid
  as.numeric(orig_data$lrgdpch),
  as.numeric(orig_data$logdstab),
  as.numeric(orig_data$L1anyaid),
  as.numeric(orig_data$hic),
  as.numeric(orig_data$T3),
  as.numeric(orig_data$endUNSC),
  as.numeric(orig_data$lpopWB),
  as.numeric(orig_data$lGDPpcWB),
  as.numeric(orig_data$year4SC),
  as.numeric(orig_data$delta4GDPpcWB),
  as.numeric(orig_data$unmem),
  as.numeric(orig_data$delta4demaut),
  as.numeric(orig_data$delta4score),
  as.numeric(orig_data$delta4tau),
  as.numeric(orig_data$delta2GDPpcWB),
  as.numeric(orig_data$delta2demaut),
  as.numeric(orig_data$delta2score),
  as.numeric(orig_data$delta2tau),
  as.numeric(orig_data$year2SC)
  
)

na.data<- as.data.frame(na.data)

colnames(na.data) <- c("year", "ccode", "scmem", "UNSC","democ", "autoc", 
                       "lpop", "GDPpcWB", "AID",
                       "demaut","anyaid",
                       "lrgdpch","logdstab","L1anyaid","hic",
                       "T3", "endUNSC", "lpopWB", "lGDPpcWB",
                       "year4SC","delta4GDPpcWB","unmem",
                       "delta4demaut", "delta4score", "delta4tau",
                       "delta2GDPpcWB",
                       "delta2demaut", "delta2score", "delta2tau",
                       "year2SC")


# drop missing 

# drop missing treatment vars
completeFunc <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}
na.data<-completeFunc(na.data, "year4SC")
#na.data<-completeFunc(na.data, "year2SC")


# Fill na with the mean for the rest 
# (linear regression shoudl not be impacted by the mean)
for(i in 1:ncol(na.data)){
  na.data[is.na(na.data[,i]), i] <- mean(na.data[,i], na.rm = TRUE)
}


## Gen matching  (DELTA GDP)

Tr=na.data[na.data$unmem==1,]$year4SC
Y=na.data[na.data$unmem==1,]$delta4GDPpcWB

X=cbind(na.data[na.data$unmem==1,]$demaut,
        na.data[na.data$unmem==1,]$lpopWB,
        na.data[na.data$unmem==1,]$lGDPpcWB)

C <- cov(X)

GDP4.genout <- GenMatch(Tr=Tr, X=X, estimand="ATT", M=1, pop.size=16, max.generations=10, wait.generations=1)
GDP4.mout <- Match(Y=Y, Tr=Tr, X=X, estimand="ATT", Weight.matrix=GDP4.genout, M=1)

summary(GDP4.mout)

'''
Estimate...  -0.94757 
AI SE......  0.83625 
T-stat.....  -1.1331 
p.val......  0.25716 

Original number of observations..............  6728 
Original number of treated obs...............  253 
Matched number of observations...............  253 
Matched number of observations  (unweighted).  3188 
'''
GDP4.G.mb <- MatchBalance(year4SC ~ lGDPpcWB+lpop+demaut,data=na.data, match.out=GDP4.mout, nboots=500)
'''
Before Matching Minimum p.value: < 2.22e-16 
Variable Name(s): lpop demaut  Number(s): 2 3 

After Matching Minimum p.value: < 2.22e-16 
Variable Name(s): lpop  Number(s): 2  
'''

# multivairate

X=na.data[,!names(na.data) %in%
            c("year4SC","delta4GDPpcWB")]

C <- cov(X)

GDP4.genout <- GenMatch(Tr=Tr, X=X, estimand="ATT", M=1, pop.size=16, max.generations=10, wait.generations=1)
GDP4.mout <- Match(Y=Y, Tr=Tr, X=X, estimand="ATT", Weight.matrix=GDP4.genout, M=1)

summary(GDP4.mout)
'''
Estimate...  0.63876 
AI SE......  0.502 
T-stat.....  1.2724 
p.val......  0.20322 

Original number of observations..............  6728 
Original number of treated obs...............  253 
Matched number of observations...............  253 
Matched number of observations  (unweighted).  253 
'''
GDP4.G.mb<- MatchBalance(year4SC ~ lGDPpcWB+lpop+demaut,data=na.data, match.out=GDP4.mout, nboots=500)
'''
Before Matching Minimum p.value: < 2.22e-16 
Variable Name(s): lpop demaut  Number(s): 2 3 

After Matching Minimum p.value: 0.12634 
Variable Name(s): lpop  Number(s): 2 
'''

## NN matching  (DELTA DEMOC)
Tr=na.data[na.data$unmem==1,]$year4SC
Y=na.data[na.data$unmem==1,]$delta4demaut

X=cbind(
  na.data[na.data$unmem==1,]$lpopWB,
  na.data[na.data$unmem==1,]$lGDPpcWB)

C <- cov(X)

DEMOC4.genout <- GenMatch(Tr=Tr, X=X, estimand="ATT", M=1, pop.size=16, max.generations=10, wait.generations=1)
DEMOC4.mout <- Match(Y=Y, Tr=Tr, X=X, estimand="ATT", Weight.matrix=DEMOC4.genout, M=1)

summary(DEMOC4.mout)
'''
Estimate...  -0.00010363 
AI SE......  0.011426 
T-stat.....  -0.0090704 
p.val......  0.99276 

Original number of observations..............  6728 
Original number of treated obs...............  253 
Matched number of observations...............  253 
Matched number of observations  (unweighted).  32398 
'''

mb <- MatchBalance(year4SC ~ lGDPpcWB+lpop+demaut,data=na.data, match.out=DEMOC4.mout, nboots=500)
'''
Before Matching Minimum p.value: < 2.22e-16 
Variable Name(s): lpop demaut  Number(s): 2 3 

After Matching Minimum p.value: < 2.22e-16 
Variable Name(s): lpop demaut  Number(s): 2 3 
'''

# multivairate

X=na.data[,!names(na.data) %in%
            c("year4SC","delta4demaut")]

C <- cov(X)

DEMOC4.genout <- GenMatch(Tr=Tr, X=X, estimand="ATT", M=1, pop.size=16, max.generations=10, wait.generations=1)
DEMOC4.mout <- Match(Y=Y, Tr=Tr, X=X, estimand="ATT", Weight.matrix=DEMOC4.genout, M=1)

summary(DEMOC4.mout)
'''
Estimate...  0.0084387 
AI SE......  0.0067026 
T-stat.....  1.259 
p.val......  0.20803 

Original number of observations..............  6728 
Original number of treated obs...............  253 
Matched number of observations...............  253 
Matched number of observations  (unweighted).  253 
'''
DEMOC4.G.mb <- MatchBalance(year4SC ~ lGDPpcWB+lpop+demaut,data=na.data, match.out=DEMOC4.mout, nboots=500)
'''
Before Matching Minimum p.value: < 2.22e-16 
Variable Name(s): lpop demaut  Number(s): 2 3 

After Matching Minimum p.value: 0.23454 
Variable Name(s): lpop  Number(s): 2 
'''

## NN matching (DELTA PRESS) X
Tr=na.data[na.data$unmem==1,]$year4SC
Y=na.data[na.data$unmem==1,]$delta4score

X=cbind(
  na.data[na.data$unmem==1,]$demaut,
  na.data[na.data$unmem==1,]$lpopWB,
  na.data[na.data$unmem==1,]$lGDPpcWB)

C <- cov(X)

PRESS4.genout <- GenMatch(Tr=Tr, X=X, estimand="ATT", M=1, pop.size=16, max.generations=10, wait.generations=1)
PRESS4.mout <- Match(Y=Y, Tr=Tr, X=X, estimand="ATT", Weight.matrix=PRESS4.genout, M=1)

summary(PRESS4.mout)
'''
Estimate...  0.34899 
AI SE......  0.28771 
T-stat.....  1.213 
p.val......  0.22514 

Original number of observations..............  6728 
Original number of treated obs...............  253 
Matched number of observations...............  253 
Matched number of observations  (unweighted).  3188  
'''

PRESS4.G.mb <- MatchBalance(year4SC ~ lGDPpcWB+lpop+demaut,data=na.data, match.out=PRESS4.mout, nboots=500)
'''
Before Matching Minimum p.value: < 2.22e-16 
Variable Name(s): lpop demaut  Number(s): 2 3 

After Matching Minimum p.value: < 2.22e-16 
Variable Name(s): lpop  Number(s): 2 
'''
##mutlivariate 
# multivairate

X=na.data[,!names(na.data) %in%
            c("year4SC","delta4score")]

C <- cov(X)

PRESS4.genout <- GenMatch(Tr=Tr, X=X, estimand="ATT", M=1, pop.size=16, max.generations=10, wait.generations=1)
PRESS4.mout <- Match(Y=Y, Tr=Tr, X=X, estimand="ATT", Weight.matrix=PRESS4.genout, M=1)

summary(PRESS4.mout)
'''
Estimate...  -0.1578 
AI SE......  0.18001 
T-stat.....  -0.87664 
p.val......  0.38068 

Original number of observations..............  6728 
Original number of treated obs...............  253 
Matched number of observations...............  253 
Matched number of observations  (unweighted).  253
'''
PRESS4.G.mb <- MatchBalance(year4SC ~ lGDPpcWB+lpop+demaut,data=na.data, match.out=PRESS4.mout, nboots=500)
'''
Before Matching Minimum p.value: < 2.22e-16 
Variable Name(s): lpop demaut  Number(s): 2 3 

After Matching Minimum p.value: 0.13471  
Variable Name(s): lpop  Number(s): 2
'''

## NN matching  (DELTA US ALIGN)
Tr=na.data[na.data$unmem==1,]$year4SC
Y=na.data[na.data$unmem==1,]$delta4tau

X=cbind(
  na.data[na.data$unmem==1,]$demaut,
  na.data[na.data$unmem==1,]$lpop,
  na.data[na.data$unmem==1,]$lGDPpcWB)

C <- cov(X)


USALIGN4.genout <- GenMatch(Tr=Tr, X=X, estimand="ATT", M=1, pop.size=16, max.generations=10, wait.generations=1)
USALIGN4.mout <- Match(Y=Y, Tr=Tr, X=X, estimand="ATT", Weight.matrix=USALIGN4.genout, M=1)

summary(USALIGN4.mout)
'''
Estimate...  0.0013652 
AI SE......  0.0090894 
T-stat.....  0.15019 
p.val......  0.88061 

Original number of observations..............  6728 
Original number of treated obs...............  253 
Matched number of observations...............  253 
Matched number of observations  (unweighted).  714 
'''
USALIGN4.G.mb <- MatchBalance(year4SC ~ lGDPpcWB+lpop+demaut,data=na.data, match.out=USALIGN4.mout, nboots=500)
'''
Before Matching Minimum p.value: < 2.22e-16 
Variable Name(s): lpop demaut  Number(s): 2 3 

After Matching Minimum p.value: < 2.22e-16 
Variable Name(s): lpop  Number(s): 2 
'''


# multivairate

X=na.data[,!names(na.data) %in%
            c("year4SC","delta4tau")]

C <- cov(X)

USALIGN4.genout <- GenMatch(Tr=Tr, X=X, estimand="ATT", M=1, pop.size=16, max.generations=10, wait.generations=1)
USALIGN4.mout <- Match(Y=Y, Tr=Tr, X=X, estimand="ATT", Weight.matrix=USALIGN4.genout, M=1)

summary(USALIGN4.mout)
'''
Estimate...  -0.0074856 
AI SE......  0.0059227 
T-stat.....  -1.2639 
p.val......  0.20627 

Original number of observations..............  6728 
Original number of treated obs...............  253 
Matched number of observations...............  253 
Matched number of observations  (unweighted).  253
'''
USALIGN4.G.mb <- MatchBalance(year4SC ~ lGDPpcWB+lpop+demaut,data=na.data, match.out=USALIGN4.mout, nboots=500)
'''
Before Matching Minimum p.value: < 2.22e-16 
Variable Name(s): lpop demaut  Number(s): 2 3 

After Matching Minimum p.value: 0.26514 
Variable Name(s): lpop  Number(s): 2 

'''

## 4 years conclusion ##
match_table <- matrix(c( paste(round(GDP4.match.out$est[1], digits = 4), "(",round(GDP4.match.out$se[1], digits = 4), ")",sep=""),
GDP4.mb$AMsmallest.p.value,
paste(round(GDP4.mout$est[1], digits = 4), "(",round(GDP4.mout$se[1], digits = 4), ")",sep=""),
GDP4.G.mb$AMsmallest.p.value,

paste(round(DEMOC4.match.out$est[1], digits = 4), "(",round(DEMOC4.match.out$se[1], digits = 4), ")",sep=""),
DEMOC4.mb$AMsmallest.p.value,
paste(round(DEMOC4.mout$est[1], digits = 4), "(",round(DEMOC4.mout$se[1], digits = 4), ")",sep=""),
DEMOC4.G.mb$AMsmallest.p.value,

paste(round(PRESS4.match.out$est[1], digits = 4), "(",round(PRESS4.match.out$se[1], digits = 4), ")",sep=""),
PRESS4.mb$AMsmallest.p.value,
paste(round(PRESS4.mout$est[1], digits = 4), "(",round(PRESS4.mout$se[1], digits = 4), ")",sep=""),
PRESS4.G.mb$AMsmallest.p.value,

paste(round(USALIGN4.match.out$est[1], digits = 4), "(",round(USALIGN4.match.out$se[1], digits = 4), ")",sep=""),
USALIGN4.mb$AMsmallest.p.value,
paste(round(USALIGN4.mout$est[1], digits = 4), "(",round(USALIGN4.mout$se[1], digits = 4), ")",sep=""),
USALIGN4.G.mb$AMsmallest.p.value),

ncol=4,byrow=TRUE)

colnames(match_table) <- c("Original result","Original balance", "GenMatch result","GenMatch balance")
rownames(match_table) <- c("%∆GDPpc","∆Democracy","∆Press Freedoms", "∆U.S. Alignment")
match_table <- as.table(match_table)


print ("4 years")
match_table

##
# 2. IV reg

#y ~ ex + en | ex + in

# matched dataset
matched_for_gdp <- na.data

treat <- na.data[GDP4.genout$matches[,1],]
ctrl <- na.data[GDP4.genout$matches[,2],]

matched_for_gdp <- rbind(treat, ctrl)

m <- ivreg(delta4GDPpcWB ~  demaut +lpopWB + AID |   demaut +lpopWB +year4SC , data=matched_for_gdp)
summary(m)
'''
Call:
ivreg(formula = delta4GDPpcWB ~ demaut + lpopWB + AID | demaut + 
lpopWB + year4SC, data = matched_for_gdp)

Residuals:
Min      1Q  Median      3Q     Max 
-52.814  -6.260   1.758   5.610  55.817 

Coefficients:
Estimate Std. Error t value Pr(>|t|)
(Intercept) -24.705740 155.464237  -0.159    0.874
demaut        7.353657   9.207911   0.799    0.425
lpopWB        1.817670  11.402255   0.159    0.873
AID          -0.001939   0.072980  -0.027    0.979

Residual standard error: 13.11 on 502 degrees of freedom
Multiple R-Squared: 0.07228,	Adjusted R-squared: 0.06674 
Wald test: 15.45 on 3 and 502 DF,  p-value: 1.247e-09 
'''

matched_for_democ <- na.data

treat <- na.data[DEMOC4.genout$matches[,1],]
ctrl <- na.data[DEMOC4.genout$matches[,2],]

matched_for_democ <- rbind(treat, ctrl)

m <- ivreg(delta4demaut ~   GDPpcWB+lpopWB + AID |   GDPpcWB +lpopWB +year4SC , data=matched_for_democ)
summary(m)
'''
Call:
ivreg(formula = delta4demaut ~ GDPpcWB + lpopWB + AID | GDPpcWB + 
lpopWB + year4SC, data = matched_for_democ)

Residuals:
Min        1Q    Median        3Q       Max 
-0.855923 -0.014442  0.006923  0.044156  0.827646 

Coefficients:
Estimate Std. Error t value Pr(>|t|)
(Intercept)  3.756e-01  1.536e+00   0.244    0.807
GDPpcWB      2.509e-06  8.608e-06   0.292    0.771
lpopWB      -2.822e-02  1.126e-01  -0.251    0.802
AID          2.064e-04  7.334e-04   0.281    0.778

Residual standard error: 0.174 on 502 degrees of freedom
Multiple R-Squared: -0.2308,	Adjusted R-squared: -0.2381 
Wald test: 0.1679 on 3 and 502 DF,  p-value: 0.918 
'''

matched_for_press <- na.data

treat <- na.data[PRESS4.genout$matches[,1],]
ctrl <- na.data[PRESS4.genout$matches[,2],]

matched_for_press <- rbind(treat, ctrl)

m <- ivreg(delta4score ~  demaut+ GDPpcWB+lpopWB + AID |   demaut+GDPpcWB +lpopWB +year4SC , data=matched_for_press)
summary(m)
'''
Call:
ivreg(formula = delta4score ~ demaut + GDPpcWB + lpopWB + AID | 
demaut + GDPpcWB + lpopWB + year4SC, data = matched_for_press)

Residuals:
Min      1Q  Median      3Q     Max 
-45.426  -8.485  -0.982   2.407 125.026 

Coefficients:
Estimate Std. Error t value Pr(>|t|)
(Intercept) -8.789e+01  8.182e+02  -0.107    0.915
demaut      -1.362e+00  1.604e+01  -0.085    0.932
GDPpcWB     -5.345e-04  4.631e-03  -0.115    0.908
lpopWB       6.528e+00  6.057e+01   0.108    0.914
AID         -4.166e-02  3.795e-01  -0.110    0.913

Residual standard error: 16 on 501 degrees of freedom
Multiple R-Squared: -30.68,	Adjusted R-squared: -30.93 
Wald test: 0.03116 on 4 and 501 DF,  p-value: 0.9981 
'''

matched_for_usalign <- na.data

treat <- na.data[USALIGN4.genout$matches[,1],]
ctrl <- na.data[USALIGN4.genout$matches[,2],]

matched_for_usalign <- rbind(treat, ctrl)

m <- ivreg(delta4tau ~ demaut+GDPpcWB +lpopWB + AID |   demaut+GDPpcWB +lpopWB +year4SC , data=matched_for_usalign)
summary(m)

'''
Call:
ivreg(formula = delta4tau ~ demaut + GDPpcWB + lpopWB + AID | 
demaut + GDPpcWB + lpopWB + year4SC, data = matched_for_usalign)

Residuals:
Min      1Q  Median      3Q     Max 
-5.7323 -0.9985 -0.1296  0.3564 16.2992 

Coefficients:
Estimate Std. Error t value Pr(>|t|)
(Intercept) -1.123e+01  2.661e+02  -0.042    0.966
demaut      -3.620e-02  1.431e+00  -0.025    0.980
GDPpcWB     -6.468e-05  1.526e-03  -0.042    0.966
lpopWB       8.268e-01  1.957e+01   0.042    0.966
AID         -5.210e-03  1.231e-01  -0.042    0.966

Residual standard error: 1.957 on 501 degrees of freedom
Multiple R-Squared: -242.8,	Adjusted R-squared: -244.8 
Wald test: 0.002678 on 4 and 501 DF,  p-value: 1 
'''

##
# 3. QR
## same original 3 vars :demaut, lpopWB, lGDPpcWB


#GDP


Y=matched_for_gdp$delta4GDPpcWB
X=cbind(pop=matched_for_gdp$lpopWB , demaut=matched_for_gdp$demaut)


QR_M=rq(Y ~X, tau=seq(0.1, .9, by=0.1))
plot(summary(QR_M))

# DEMOC


Y=matched_for_democ$delta4demaut

X=cbind(pop=matched_for_democ$lpopWB , demaut=matched_for_democ$lGDPpcWB)


QR_M=rq(Y ~X, tau=seq(0.1, .9, by=0.1))
plot(summary(QR_M))

matched_for_usalign <- na.data


# PRESS



Y=matched_for_press$delta4score

X=cbind(pop=matched_for_press$lpopWB , demaut=matched_for_press$demaut, demaut=matched_for_press$lGDPpcWB)


QR_M=rq(Y ~X, tau=seq(0.1, .9, by=0.1))
plot(summary(QR_M))



## US ALIGNMENT


Y=matched_for_usalign$delta4tau
X=cbind(pop=matched_for_usalign$lpopWB , demaut=matched_for_usalign$demaut ,GDP = matched_for_usalign$GDPpcWB)


QR_M=rq(Y ~X, tau=seq(0.1, .9, by=0.1))
plot(summary(QR_M))

