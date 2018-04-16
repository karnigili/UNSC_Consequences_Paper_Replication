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
library(rbounds)
library(xtable)
library(mfx)
library(spelling)

###################
### IMPORT DATA ###
###################
orig_data <- read.dta("/Users/gili/drive/minerva/2/cs112/assignments/proposal research/foriegn_aid/PerniciousEffectUNSC/USaid_UNSC_selectorate_data_replication.dta")

###############################################
### DATA PRE PROCESS AND FEATURE EXTRACTION ###
###############################################

# panel the data - order bt years and ccode
orig_data <- orig_data[order(orig_data$ccode, orig_data$year),]

## lagging func
## using slide to mimic L and F from STATA
L <- function(var_name,  n) {
  Data <- data.frame(orig_data$year, orig_data$ccode, orig_data[,var_name])
  colnames(Data) <- c("year", "ccode", var_name)
  Data <- Data[order(Data$ccode, Data$year),]
  
  DataSlid2 <- slide(Data, Var = var_name, GroupVar = "ccode",
                     slideBy = n, keepInvalid=TRUE )
  
  return(as.numeric(unlist(DataSlid2[4])))
}


## HowMuchAidAndUNSCreplication.smcl ##

orig_data$Ltotaid96 = L("totaid96", -1)

orig_data$F2.T0 <- L("T0",2)
orig_data$F1.T0 <- L("T0",1)
orig_data$F3.T0 <- L("T0",3)
# setting the treaatment variables:

# 1. year4SC
# signifies a country which was elected to UNSC this year and not already a current member (or was in the past 4 years)
orig_data$year4SC = ifelse(orig_data$T0==0&orig_data$T1==0&orig_data$T2==0&orig_data$T3==0&orig_data$T4==0&orig_data$F1.T0==0&orig_data$F2.T0==0&orig_data$F3.T0==0&orig_data$unmem==1 , 0, NA)
orig_data$year4SC = ifelse(orig_data$T0==1&orig_data$unmem==1, 1, orig_data$year4SC)

# 1. year2SC
# signifies a country which was elected to UNSC this year and not already a current member (or was in the past 2 years)
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

## paired t test to comapre aid before and during
t.test(orig_data$aidduring,orig_data$aidbefore, paired=TRUE)

'''
Paired t-test

data:  orig_data$aidduring and orig_data$aidbefore
t = 7.154, df = 7349, p-value = 9.241e-13
alternative hypothesis: true difference in means is not equal to 0

95 percent confidence interval:
2.900485 5.089964
sample estimates:
mean of the differences 
3.995225 
'''


#UNSC includes the year of election as memebership
orig_data$UNSC=orig_data$scmem
orig_data$UNSC= ifelse(orig_data$T0==1,1,orig_data$UNSC)
orig_data$UNSCimptyr = orig_data$imptyr * orig_data$UNSC

orig_data$endUNSC <- ifelse(orig_data$unmem==1, 0, NA)
orig_data$endUNSC <- ifelse(orig_data$T3==1|orig_data$T4==1, 1, orig_data$unmem)


## Does UNSC increase the likelihood of getting aid ##
t_unsc <- orig_data[orig_data$hic==0,]$UNSC
t_anyaid <- orig_data[orig_data$hic==0,]$anyaid

table(t_unsc, t_anyaid)
'''
t_anyaid
t_unsc    0    1
0 1951 4098
1   46  458
'''
t_unsc <- orig_data[(orig_data$hic==0& orig_data$L1anyaid==1),]$UNSC
t_anyaid <- orig_data[(orig_data$hic==0& orig_data$L1anyaid==1),]$anyaid
table(t_unsc, t_anyaid)
'''
t_anyaid
t_unsc    0    1
0  135 3842
1    9  439
'''

#### The impact of SC membership on aid using probit regression
## * Probit regression * describes the classification probability. The coefficient values are
## given on a normal distributed CDF, thus represet the z-value that corresponds to a specific cumulative probability
## we can interpret them as the standardtaized normalized diffrence associated with each covariates
##
## *The intercept* is simply the expected mean value of Y where X=0.
## Dummy coded variables have values of 0 for the reference group and 1 for the comparison group. 
## Since the intercept is the expected mean value when X=0, 
## it is the mean value only for the reference group (when all other X=0).

# These probit regressions using similar covriates to measure the impact of UNSC membership on aid, each one
# uses different binary covariate.

#1. scmem, simply whether a country was a member at a given year
# results: the intercept = 1.56 -> 56% increase in aid corresponding to scmem covariate
probit <- glm( anyaid~ L1anyaid+demaut+scmem+lpop+lrgdpch+logdstab , data=orig_data[orig_data$hic == 0,],family=binomial(link = "probit"))
summary(probit)

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

#2. UNSC, including the year of election 
# results: the intercept = 1.57 -> 57% increase in aid corresponding to UNSC covariate
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

#3. T0+T1+T2+T3+T4, election year, membership years, and the two following years 
# results: the intercept = 1.54 -> 54% increase in aid corresponding to Ts covariate
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

# conclusion: there is an increase of ~56% in aid to UNSC members

#### TABLE2: correaltion between UNSC membership to the 4 outcomes : gdp, democracy, press, US alignment ##
# Table 2 contains probit analyses that assess which factors influence election to the UNSC

## line 20
table2_1 <- summary(glm( T0~ demaut+lpopWB+lGDPpcWB+tau_lead, 
                         data=orig_data[(orig_data$unmem==1 & orig_data$P5 !=1 & !(orig_data$T1==1|orig_data$T2==1|orig_data$T3==1|orig_data$T4==1)),],
                         family=binomial(link = "probit")))

table2_1
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
table2_2 <- summary(glm( T0~ demaut+lpopWB+lGDPpcWB+tau_lead+AidGDP, 
                         data=orig_data[(orig_data$unmem==1 & orig_data$P5 !=1 & !(orig_data$T1==1|orig_data$T2==1|orig_data$T3==1|orig_data$T4==1)),],
                         family=binomial(link = "probit")))
table2_2
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

table_2 <- matrix(c(
  paste(round(table2_1$coefficients[2], 4),"(", round(table2_1$coefficients[7], 4),")"), paste(round(table2_2$coefficients[2], 4),"(", round(table2_2$coefficients[7], 4),")"),
  paste(round(table2_1$coefficients[3], 4),"(", round(table2_1$coefficients[8], 4),")"), paste(round(table2_2$coefficients[3], 4),"(", round(table2_2$coefficients[8], 4),")"),
  paste(round(table2_1$coefficients[4], 4),"(", round(table2_1$coefficients[9], 4),")"), paste(round(table2_2$coefficients[4], 4),"(", round(table2_2$coefficients[9], 4),")"),
  paste(round(table2_1$coefficients[5], 4),"(", round(table2_1$coefficients[10], 4),")"), paste(round(table2_2$coefficients[5], 4),"(", round(table2_2$coefficients[10], 4),")"),
  "", paste(round(table2_2$coefficients[6], 4),"(", round(table2_2$coefficients[11], 4),")"),
  "","",
  paste(round(table2_1$coefficients[1], 4),"(", round(table2_1$coefficients[6], 4),")"), paste(round(table2_2$coefficients[1], 4),"(", round(table2_2$coefficients[6], 4),")"),
  table2_1$df[1]+table2_1$df[2], table2_2$df[1]+table2_2$df[2]),ncol=2,byrow=TRUE)

colnames(table_2) <- c("Probit Analysis of Election to","UNSC for all Non-Permanent UN Members")
rownames(table_2) <- c("Democracy","Log(Population)","Log(GDPpc)", "U.S. Security Alignment","Aid/GDP(%)","", "Constant", "Observations")
table_2 <- as.table(table_2)

table2_latex <- xtable(table_2, , digits = 0)

'''
\begin{table}[ht]
\centering
\begin{tabular}{rll}
\hline
& Probit Analysis of Election to & UNSC for all Non-Permanent UN Members \\ 
\hline
Democracy & 0.034 ( 0.1149 ) & -0.0224 ( 0.7021 ) \\ 
Log(Population) & 0.1693 ( 0.0273 ) & 0.1423 ( 0.1295 ) \\ 
Log(GDPpc) & 0.1072 ( 0.0287 ) & 0.0583 ( 0.032 ) \\ 
U.S. Security Alignment & 0.0064 ( 0.1179 ) & 0.0284 ( 0.0436 ) \\ 
Aid/GDP(\%) &  & -0.009 ( 0.1308 ) \\ 
&  &  \\ 
Constant & -5.1902 ( 0.5057 ) & -4.3458 ( -0.009 ) \\ 
Observations & 3586 & 2991 \\ 
\hline
\end{tabular}
\end{table}
'''

### the other 3 cols of the table requires a var names TAU which I could not replicate

####### Effects of Membership of Security Council ######

## feature exctartion 

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



## TABLE 3: shows the average treatment effect of UNSC membership using a matching design ##
## Using NN matching and linear regression
## four prospective outcomes: GDP, democracy score, press freedom, and US alignment
## each is regressed on Log(Population), Log(GDPpc), and Democracy score.
## We measure the effect on the span of 4 years and 2 years seperately.



## 4 years ##

# create a table containing smaller group of covaraites, to preform regression we need an table without Null values, reducing the
# number of covariates enables maintaing a large numbers of observations
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


# Drop missing values

# 1. drop missing values for the treatment variables
completeFunc <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}
na.data<-completeFunc(na.data, "year4SC")


# Fill the null values with the mean for the rest 
# (linear regression should not be impacted by observations with the values of the mean of the population)
for(i in 1:ncol(na.data)){
  na.data[is.na(na.data[,i]), i] <- mean(na.data[,i], na.rm = TRUE)
}


## NN matching line 1847 (DELTA-GDP)

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


## 2 years ##
## equivilant process
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


# drop Null values
# treatment covariates
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


table_3 <- matrix(c( paste(round(GDP4.match.out$est[1], digits = 4), "(",round(GDP4.match.out$se[1], digits = 4), ")",sep=""),
                     paste(round(GDP2.match.out$est[1], digits = 4), "(",round(GDP2.match.out$se[1], digits = 4), ")",sep=""),
      
                     paste(round(DEMOC4.match.out$est[1], digits = 4), "(",round(DEMOC4.match.out$se[1], digits = 4), ")",sep=""),
                     paste(round(DEMOC2.match.out$est[1], digits = 4), "(",round(DEMOC2.match.out$se[1], digits = 4), ")",sep=""),
                     
                     paste(round(PRESS4.match.out$est[1], digits = 4), "(",round(PRESS4.match.out$se[1], digits = 4), ")",sep=""),
                     paste(round(PRESS2.match.out$est[1], digits = 4), "(",round(PRESS2.match.out$se[1], digits = 4), ")",sep=""),
                    
                      paste(round(USALIGN4.match.out$est[1], digits = 4), "(",round(USALIGN4.match.out$se[1], digits = 4), ")",sep=""),
                      paste(round(USALIGN2.match.out$est[1], digits = 4), "(",round(USALIGN2.match.out$se[1], digits = 4), ")",sep="")),
                     
                      ncol=2,byrow=TRUE)

colnames(table_3) <- c("Average Treatment Effect of UNSC Four Years","Average Treatment Effect of UNSC Two Years")
rownames(table_3) <- c("%∆GDPpc","∆Democracy","∆Press Freedoms", "∆U.S. Alignment")
table_3 <- as.table(table_3)
latex_table_3 <- xtable(table_3)
'''
\begin{table}[ht]
\centering
\begin{tabular}{rll}
\hline
& Average Treatment Effect of UNSC Four Years & Average Treatment Effect of UNSC Two Years \\ 
\hline
\%∆GDPpc & -2.8789(0.9562) & -0.5501(0.6991) \\ 
∆Democracy & -0.0177(0.0112) & -0.0104(0.0106) \\ 
∆Press Freedoms & 0.2717(0.2519) & 0.4108(0.2822) \\ 
∆U.S. Alignment & 0.0063(0.0081) & 0.0088(0.0081) \\ 
\hline
\end{tabular}
\end{table}
'''


########################################
########### PERSOANL CODE ##############
########################################

'''
PLAN - 4 years impact only

1. balance tests
2. GenMatch + balance
- sensitivy test
3. IV reg
4. QR
'''

### prep data ###

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


# Drop missing values

# Drop missing treatment vars for treat var
na.data<-completeFunc(na.data, "year4SC")

# Fill na with the mean value for the rest 
# (linear regression shoudl not be impacted by the mean)
for(i in 1:ncol(na.data)){
  na.data[is.na(na.data[,i]), i] <- mean(na.data[,i], na.rm = TRUE)
}

## balance ##

GDP4.mb <- MatchBalance(year4SC ~ lGDPpcWB+lpop+demaut,data=na.data, match.out=GDP4.match.out, nboots=500)

'''
Before Matching Minimum p.value: < 2.22e-16 
Variable Name(s): lGDPpcWB lpop demaut  Number(s): 1 2 3 

After Matching Minimum p.value: < 2.22e-16 
Variable Name(s): lGDPpcWB lpop demaut  Number(s): 1 2 3 
'''

DEMOC4.mb <- MatchBalance(year4SC ~ lGDPpcWB+lpop+demaut,data=na.data, match.out=DEMOC4.match.out, nboots=500)
'''
Before Matching Minimum p.value: < 2.22e-16 
Variable Name(s): lpop demaut  Number(s): 2 3 

After Matching Minimum p.value: < 2.22e-16 
Variable Name(s): lpop demaut  Number(s): 2 3
'''
PRESS4.mb <- MatchBalance(year4SC ~ lGDPpcWB+lpop+demaut,data=na.data, match.out=PRESS4.match.out, nboots=500)
'''
Before Matching Minimum p.value: < 2.22e-16 
Variable Name(s): lpop demaut  Number(s): 2 3 

After Matching Minimum p.value: < 2.22e-16 
Variable Name(s): lpop  Number(s): 2 
'''

USALIGN4.mb <- MatchBalance(year4SC ~ lGDPpcWB+lpop+demaut,data=na.data, match.out=USALIGN4.match.out, nboots=500)

'''
Before Matching Minimum p.value: < 2.22e-16 
Variable Name(s): lpop demaut  Number(s): 2 3 

After Matching Minimum p.value: 0.00086155 
Variable Name(s): lpop  Number(s): 2 

'''

### GenMatch ###

# for each one of the four outcomes I ran 2 types of GenMatch, first; using the
# 3 given covariates, then a multivariate matching using all covariates excluding the
# tretment and outcome ones

# 1.a. GDP / 3 covs
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

# 1.b. GDP / multivariate

X=na.data[,!names(na.data) %in%
            c("year4SC","delta4GDPpcWB")]

C <- cov(X)

GDP4.genout <- GenMatch(Tr=Tr, X=X, estimand="ATT", M=5, pop.size=20, max.generations=10, wait.generations=1)
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

After Matching Minimum p.value: 0.0040535 
Variable Name(s): lpop  Number(s): 2  
'''

# 2.a. Democracy / 3 covs
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

# 2.b. Democracy/ multivairate

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

After Matching Minimum p.value: 0.00036485 
Variable Name(s): lpop  Number(s): 2  
'''

# 3.a. Press / 3 covs
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

# 3.b. Press / 3 multivairiate
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

# 4.a US alignmnet / 3 cov
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

# 4.b US alignmnet / multivairate

X=na.data[,!names(na.data) %in%
            c("year4SC","delta4tau")]

C <- cov(X)

USALIGN4.genout <- GenMatch(Tr=Tr, X=X, estimand="ATT", M=1, pop.size=16, max.generations=10, wait.generations=1)
USALIGN4.mout <- Match(Y=Y, Tr=Tr, X=X, estimand="ATT", Weight.matrix=USALIGN4.genout, M=1)

summary(USALIGN4.mout)
'''
Estimate...  -0.00086388 
AI SE......  0.005446 
T-stat.....  -0.15863 
p.val......  0.87396 

Original number of observations..............  5413 
Original number of treated obs...............  253 
Matched number of observations...............  253 
Matched number of observations  (unweighted).  253 
'''
USALIGN4.G.mb <- MatchBalance(year4SC ~ lGDPpcWB+lpop+demaut,data=na.data, match.out=USALIGN4.mout, nboots=500)
'''
Before Matching Minimum p.value: < 2.22e-16 
Variable Name(s): lpop demaut  Number(s): 2 3 

After Matching Minimum p.value: 0.0049466 
Variable Name(s): lGDPpcWB  Number(s): 1 

'''

## 4 years conclusion ##
match_table <- matrix(c( paste(round(GDP4.match.out$est[1], digits = 4), "(",round(GDP4.match.out$se[1], digits = 4), ")",sep=""),
                         round(GDP4.mb$AMsmallest.p.value, digits = 4),
                         "0.007",
                         paste(round(GDP4.mout$est[1], digits = 4), "(",round(GDP4.mout$se[1], digits = 4), ")",sep=""),
                         round(GDP4.G.mb$AMsmallest.p.value, digits = 4),
                         "0.83",
                         
                         paste(round(DEMOC4.match.out$est[1], digits = 4), "(",round(DEMOC4.match.out$se[1], digits = 4), ")",sep=""),
                         round(DEMOC4.mb$AMsmallest.p.value, digits = 4),
                         "0.185", 
                         paste(round(DEMOC4.mout$est[1], digits = 4), "(",round(DEMOC4.mout$se[1], digits = 4), ")",sep=""),
                         round(DEMOC4.G.mb$AMsmallest.p.value, digits = 4),
                         "0.86",
                         
                         paste(round(PRESS4.match.out$est[1], digits = 4), "(",round(PRESS4.match.out$se[1], digits = 4), ")",sep=""),
                         round(PRESS4.mb$AMsmallest.p.value, digits = 4),
                         "0.23",
                         paste(round(PRESS4.mout$est[1], digits = 4), "(",round(PRESS4.mout$se[1], digits = 4), ")",sep=""),
                         round(PRESS4.G.mb$AMsmallest.p.value, digits = 4),
                         "0.38",
                         
                         paste(round(USALIGN4.match.out$est[1], digits = 4), "(",round(USALIGN4.match.out$se[1], digits = 4), ")",sep=""),
                         round(USALIGN4.mb$AMsmallest.p.value, digits = 4),
                         "0.35",
                         paste(round(USALIGN4.mout$est[1], digits = 4), "(",round(USALIGN4.mout$se[1], digits = 4), ")",sep=""),
                         round(USALIGN4.G.mb$AMsmallest.p.value, digits = 4),
                         "0.87"),
                      
                      ncol=6,byrow=TRUE)

colnames(match_table) <- c("Original result","original p_value","Original balance", "GenMatch result","GenMatch balance", "GenMatch p_value")
rownames(match_table) <- c("%∆GDPpc","∆Democracy","∆Press Freedoms", "∆U.S. Alignment")
match_table <- as.table(match_table)


print ("4 years")
match_table
'''
                Original result original p_value Original balance GenMatch result GenMatch balance GenMatch p_value
%∆GDPpc         -2.8789(0.9562) 0                0.007            -0.1415(0.6662) 0.225           0.83            
∆Democracy      -0.0177(0.0112) 0                0.185            0.0018(0.0103)  0.1125          0.86            
∆Press Freedoms 0.2717(0.2519)  0                0.23             0.2692(0.2661)  0.1347         0.38            
∆U.S. Alignment 0.0063(0.0081)  0.0201           0.35             -9e-04(0.0054)  0.2651          0.87            
> 
'''
latex_match<- xtable(match_table)
'''
% latex table generated in R 3.4.3 by xtable 1.8-2 package
% Mon Apr 16 16:53:34 2018
\begin{table}[ht]
\centering
\begin{tabular}{rllllll}
\hline
& Original result & original p\_value & Original balance & GenMatch result & GenMatch balance & GenMatch p\_value \\ 
\hline
\%∆GDPpc & -2.8789(0.9562) & 0 & 0.007 & -0.1415(0.6662) & 0.225 & 0.83 \\ 
∆Democracy & -0.0177(0.0112) & 0 & 0.185 & 0.0018(0.0103) & 0.1125 & 0.86 \\ 
∆Press Freedoms & 0.2717(0.2519) & 0 & 0.23 & 0.2692(0.2661) & 0.1347 & 0.38 \\ 
∆U.S. Alignment & 0.0063(0.0081) & 0.0201 & 0.35 & -9e-04(0.0054) & 0.2651 & 0.87 \\ 
\hline
\end{tabular}
\end{table}
'''

## conclusions:
# The banalce given by the authors is basically 0. Meaning the bias in their result is very high
# by increasing the balnce we have decreased the bias, yet due to decrease in the degrees of freeedom (the number of
# observations), we have lost the statistical power and the vairance of the results increased. This can be showen in the
# shift from relatively low p-values to higher p-values.

## sensitivit test ##
# 1.a. original GDP
psens(GDP4.match.out, Gamma=2, GammaInc=.2)
'''
Rosenbaum Sensitivity Test for Wilcoxon Signed Rank P-Value 
 
Unconfounded estimate ....  0 

 Gamma Lower bound Upper bound
   1.0           0      0.0000
   1.2           0      0.0000
   1.4           0      0.9005
   1.6           0      1.0000
   1.8           0      1.0000
   2.0           0      1.0000

 Note: Gamma is Odds of Differential Assignment To
 Treatment Due to Unobserved Factors 
'''
# 2.a. original democracy
psens(DEMOC4.match.out, Gamma=2, GammaInc=.2)
'''
Rosenbaum Sensitivity Test for Wilcoxon Signed Rank P-Value 
 
Unconfounded estimate ....  0 

Gamma Lower bound Upper bound
1.0           0      0.0000
1.2           0      0.9997
1.4           0      1.0000
1.6           0      1.0000
1.8           0      1.0000
2.0           0      1.0000

Note: Gamma is Odds of Differential Assignment To
Treatment Due to Unobserved Factors 
'''
# 3.a. original press freedom
psens(PRESS4.match.out, Gamma=2, GammaInc=.2)
'''
Rosenbaum Sensitivity Test for Wilcoxon Signed Rank P-Value 
 
Unconfounded estimate ....  0.0045 

 Gamma Lower bound Upper bound
   1.0      0.0045      0.0045
   1.2      0.0000      0.1300
   1.4      0.0000      0.5474
   1.6      0.0000      0.8846
   1.8      0.0000      0.9844
   2.0      0.0000      0.9987

 Note: Gamma is Odds of Differential Assignment To
 Treatment Due to Unobserved Factors 
'''
# 4.a. original US alignment
psens(USALIGN4.match.out, Gamma=2, GammaInc=.2)
'''
 Rosenbaum Sensitivity Test for Wilcoxon Signed Rank P-Value 
 
Unconfounded estimate ....  0.391 

Gamma Lower bound Upper bound
1.0      0.3910      0.3910
1.2      0.0002      0.9985
1.4      0.0000      1.0000
1.6      0.0000      1.0000
1.8      0.0000      1.0000
2.0      0.0000      1.0000

Note: Gamma is Odds of Differential Assignment To
Treatment Due to Unobserved Factors 
'''

## Conlcusion of preforming Rosenbaum Sensitivity Test for Wilcoxon Signed Rank P-Value:

# The original results for the first three outcomes (GDP, democracy, and press freedom) showed statistical 
# significance considering the impact UNCS membership has on them. These results are highly sensitive, as shwon 
# in the above tables.

# The p-value at \gamma=1 ia assumed to represent the true p-value for assuming there is no
# hidden bias due to an unobserved confounder.
# We observe that for the p-value to cross the 0.05 significance threshold, the \gamma needs to increase to
# 1.2-1.4. Meaning that the odds of a country being in the UNSC are only 1.2-1.4 times
# higher beacuse of diffrent values in an unobserved covariate (regardless og seemingly identical matched covariates).
# This imlies that even a small unobserved difference in a covariate would change our inference.

# 1.b. original GDP
hlsens(GDP4.match.out, Gamma=2, GammaInc=.1)
'''
Rosenbaum Sensitivity Test for Hodges-Lehmann Point Estimate 
 
Unconfounded estimate ....  -2.4051 

 Gamma Lower bound Upper bound
   1.0     -2.4051   -2.405100
   1.1     -2.5051    0.094855
   1.2     -2.5051    0.094855
   1.3     -2.5051    0.094855
   1.4     -2.5051    0.094855
   1.5     -2.9051    0.094855
   1.6     -3.3051    0.394860
   1.7     -3.8051    0.794860
   1.8     -4.2051    1.094900
   1.9     -4.6051    1.394900
   2.0     -5.0051    1.694900

 Note: Gamma is Odds of Differential Assignment To
 Treatment Due to Unobserved Factors 
'''

# 2.b.
hlsens(DEMOC4.match.out, Gamma=2, GammaInc=.1)
'''
Rosenbaum Sensitivity Test for Hodges-Lehmann Point Estimate 

Unconfounded estimate ....  0 

Gamma Lower bound Upper bound
1.0 -1.9247e-05 -1.9247e-05
1.1 -1.0002e-01  9.9981e-02
1.2 -1.0002e-01  9.9981e-02
1.3 -1.0002e-01  9.9981e-02
1.4 -1.0002e-01  9.9981e-02
1.5 -1.0002e-01  9.9981e-02
1.6 -1.0002e-01  9.9981e-02
1.7 -1.0002e-01  9.9981e-02
1.8 -1.0002e-01  9.9981e-02
1.9 -1.0002e-01  9.9981e-02
2.0 -1.0002e-01  9.9981e-02

Note: Gamma is Odds of Differential Assignment To
Treatment Due to Unobserved Factors 
'''
#3.b. original press freedom
hlsens(PRESS4.match.out, Gamma=2, GammaInc=.1)
'''
Rosenbaum Sensitivity Test for Hodges-Lehmann Point Estimate 
 
Unconfounded estimate ....  1 

Gamma Lower bound Upper bound
1.0    1.000000         1.0
1.1   -0.099976         1.1
1.2   -0.099976         1.1
1.3   -0.099976         1.1
1.4   -0.099976         1.1
1.5   -0.099976         1.1
1.6   -0.099976         1.1
1.7   -0.099976         1.1
1.8   -0.099976         1.1
1.9   -0.099976         1.1
2.0   -0.099976         1.1

Note: Gamma is Odds of Differential Assignment To
Treatment Due to Unobserved Factors 
'''

# 4.b. original US alignment
hlsens(USALIGN4.match.out, Gamma=2, GammaInc=.1)
'''
Rosenbaum Sensitivity Test for Hodges-Lehmann Point Estimate 
 
Unconfounded estimate ....  3e-04 

Gamma Lower bound Upper bound
1.0  0.00026426  0.00026426
1.1 -0.09973600  0.10026000
1.2 -0.09973600  0.10026000
1.3 -0.09973600  0.10026000
1.4 -0.09973600  0.10026000
1.5 -0.09973600  0.10026000
1.6 -0.09973600  0.10026000
1.7 -0.09973600  0.10026000
1.8 -0.09973600  0.10026000
1.9 -0.09973600  0.10026000
2.0 -0.09973600  0.10026000

Note: Gamma is Odds of Differential Assignment To
Treatment Due to Unobserved Factors 
'''
## Conlcusion of preforming Rosenbaum Sensitivity Test for Hodges-Lehmann Point Estimate:

# This test provides the bounds for the additive effect due to treatment. Signifying the bias introduced in the results median.
# The change in sign implies a different direction of the treatment and will be the threshold
# we use to determine the senstivity.

# This test shows evenhigher sensitivity since all outcomes flip the sign of treatment
# at \gamma =1.1. Generally, we conclude that even if it seems that the UNSC has a negative impact
# on the mentioned outcomes, the findings are sensitive to possible hidden bias due to an unobserved confounder.



### IV regression ###

# Encorgement design enables using
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


QR_M=rq(Y~X, tau=seq(0.1, .9, by=0.1))
plot(summary(QR_M))

# DEMOC
Y=matched_for_democ$delta4demaut
X=cbind(pop=matched_for_democ$lpopWB , demaut=matched_for_democ$lGDPpcWB)


QR_M=rq(Y ~X, tau=seq(0.1, .9, by=0.1))
plot(summary(QR_M))




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




