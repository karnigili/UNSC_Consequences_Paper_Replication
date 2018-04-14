/* Replication batch file for "the pernicious effects of the UNSC-- JCR " */
/* Alastair Smith Novermber 20th 2008 */
/* I want to test the predictions of selectorate aid-for-policy favors model */
/* Use Kuzuemko and Werker (JPR 2006) UN security council data */ 
/* data taken straight from JPE website -- data organized by World bank Codes */ 

/* This file runs in stata 10 and requires the data file USaid_UNSC_selectorate_data_replication.dta */

#delimit;
clear; set mem 80m;  
set more off; 


use USaid_UNSC_selectorate_data_replication.dta; 


sort ccode year;
drop if (ccode ==ccode[_n-1]& year==year[_n-1]);

tsset ccode year;  
gen Lln_totaid=L.ln_totaid96;

/***************************** Election to UN Security Council **************************/


/*** Examine regions ***/ 
gen latin =0 if ccode~=.;  replace latin=1 if ccode>20 & ccode<200;
gen WestEur =0 if ccode~=.; replace WestEur=1 if (ccode >=200 & ccode < 355)
			& ccode~=265 & ccode~=290&ccode~=310& ccode~=315 &ccode~=316; 
replace WestEur=1 if (ccode >=375 & ccode < 400); /* scandinavia */
replace WestEur=1 if ccode==20|ccode==900|ccode==920; /*Can, Aus and NZ grouped as other in W Europe category*/ 

gen EastEur=0 if ccode~=.; replace EastEur=1 if ccode==265 | ccode==290|ccode==310| ccode==315 |ccode==316
									|(ccode>=355 &ccode<=373); 

gen Africa=0 if ccode~=.; replace Africa=1 if ccode>=400 & ccode<=625;
gen Asia=0 if ccode~=.; replace Asia=1 if ccode>=630 &ccode~=900& ccode~=920& ccode~=666; 
 
capture log close; log using USaidUNSCreplication,replace; 

#delimit;
tsset ccode year; capture gen L4demaut=L4.demaut; capture gen delta4demaut=F4.demaut-demaut; 
capture gen L4delta4demaut=L4.delta4demaut;

probit T0 demaut lpopWB  lGDPpcWB tau_lead if unmem==1 & P5~=1 & ~(T1==1|T2==1|T3==1|T4==1); 
probit T0 demaut lpopWB  lGDPpcWB tau_lead AidGDP if unmem==1 & P5~=1 & ~(T1==1|T2==1|T3==1|T4==1); 
 
probit T0 demaut lpop  lrgdpch  tau_lead if unmem==1 & P5~=1 & ~(T1==1|T2==1|T3==1|T4==1); 
probit T0 demaut lpop  lrgdpch  tau if unmem==1 & P5~=1 & ~(T1==1|T2==1|T3==1|T4==1); 
probit T0 demaut lpop  lrgdpch  tau_lead  Lln_totaid if unmem==1 & P5~=1 & ~(T1==1|T2==1|T3==1|T4==1); 

probit T0 demaut L4delta4demaut lpopWB  lGDPpcWB tau_lead if unmem==1 & P5~=1 & ~(T1==1|T2==1|T3==1|T4==1); 


/* Eastern Europe nations only elected in odd years */
gen even=(year/2==ceil(year/2)); 
gen EastEven=EastEur*even;  /* nation can not be elected to UNSC if this variable is a 1 */
#delimit;
probit T0 demaut lpopWB  lGDPpcWB tau_lead if unmem==1 & P5~=1 & ~(T1==1|T2==1|T3==1|T4==1|EastEven==1); 
probit T0 demaut lpopWB  lGDPpcWB tau_lead AidGDP if unmem==1 & P5~=1 & ~(T1==1|T2==1|T3==1|T4==1|EastEven==1); 
 
probit T0 demaut lpop  lrgdpch  tau_lead if unmem==1 & P5~=1 & ~(T1==1|T2==1|T3==1|T4==1|EastEven==1); 
probit T0 demaut lpop  lrgdpch  tau if unmem==1 & P5~=1 & ~(T1==1|T2==1|T3==1|T4==1|EastEven==1); 
 probit T0 demaut lpop  lrgdpch  tau_lead  Lln_totaid if unmem==1 & P5~=1 & ~(T1==1|T2==1|T3==1|T4==1|EastEven==1); 
probit T0 demaut L4delta4demaut lpopWB  lGDPpcWB tau_lead if unmem==1 & P5~=1 & ~(T1==1|T2==1|T3==1|T4==1|EastEven==1);




/* Latin America */
preserve; keep if latin==1;
probit T0 demaut lpopWB  lGDPpcWB tau_lead if unmem==1 & P5~=1 & ~(T1==1|T2==1|T3==1|T4==1|EastEven==1); 
probit T0 demaut lpopWB  lGDPpcWB tau_lead AidGDP if unmem==1 & P5~=1 & ~(T1==1|T2==1|T3==1|T4==1|EastEven==1); 
 
probit T0 demaut lpop  lrgdpch  tau_lead if unmem==1 & P5~=1 & ~(T1==1|T2==1|T3==1|T4==1|EastEven==1); 
probit T0 demaut lpop  lrgdpch  tau if unmem==1 & P5~=1 & ~(T1==1|T2==1|T3==1|T4==1|EastEven==1); 
 probit T0 demaut lpop  lrgdpch  tau_lead  Lln_totaid if unmem==1 & P5~=1 & ~(T1==1|T2==1|T3==1|T4==1|EastEven==1); 
probit T0 demaut L4delta4demaut lpopWB  lGDPpcWB tau_lead if unmem==1 & P5~=1 & ~(T1==1|T2==1|T3==1|T4==1|EastEven==1); 
restore; 

/*West Europe*/
#delimit;
preserve; keep if WestEur==1;
probit T0 demaut lpopWB  lGDPpcWB tau_lead if unmem==1 & P5~=1 & ~(T1==1|T2==1|T3==1|T4==1); 

 probit T0 demaut lpop  lrgdpch  tau_lead if unmem==1 & P5~=1 & ~(T1==1|T2==1|T3==1|T4==1);
probit T0 demaut lpop  lrgdpch  tau if unmem==1 & P5~=1 & ~(T1==1|T2==1|T3==1|T4==1);  
probit T0 demaut lpop  lrgdpch  tau_lead  Lln_totaid if unmem==1 & P5~=1 & ~(T1==1|T2==1|T3==1|T4==1); 
probit T0 demaut L4delta4demaut lpopWB  lGDPpcWB tau_lead if unmem==1 & P5~=1 & ~(T1==1|T2==1|T3==1|T4==1); 
restore; 
/* East Europe */
preserve; keep if EastEur==1;
probit T0 demaut lpopWB  lGDPpcWB tau_lead if unmem==1 & P5~=1 & ~(T1==1|T2==1|T3==1|T4==1|EastEven==1); 
probit T0 demaut lpop  lrgdpch  tau_lead if unmem==1 & P5~=1 & ~(T1==1|T2==1|T3==1|T4==1|EastEven==1);
probit T0 demaut lpop  lrgdpch  tau if unmem==1 & P5~=1 & ~(T1==1|T2==1|T3==1|T4==1|EastEven==1);  
probit T0 demaut lpop  lrgdpch  tau_lead  Lln_totaid if unmem==1 & P5~=1 & ~(T1==1|T2==1|T3==1|T4==1|EastEven==1); 
probit T0 demaut L4delta4demaut lpopWB  lGDPpcWB tau_lead if unmem==1 & P5~=1 & ~(T1==1|T2==1|T3==1|T4==1|EastEven==1); 

restore; 
/* Africa */
preserve; keep if Africa==1;
probit T0 demaut lpopWB  lGDPpcWB tau_lead if unmem==1 & P5~=1 & ~(T1==1|T2==1|T3==1|T4==1); 
probit T0 demaut lpopWB  lGDPpcWB tau_lead AidGDP if unmem==1 & P5~=1 & ~(T1==1|T2==1|T3==1|T4==1); 
 
probit T0 demaut lpop  lrgdpch  tau_lead if unmem==1 & P5~=1 & ~(T1==1|T2==1|T3==1|T4==1);
probit T0 demaut lpop  lrgdpch  tau if unmem==1 & P5~=1 & ~(T1==1|T2==1|T3==1|T4==1);  
probit T0 demaut lpop  lrgdpch  tau_lead  Lln_totaid if unmem==1 & P5~=1 & ~(T1==1|T2==1|T3==1|T4==1); 
probit T0 demaut L4delta4demaut lpopWB  lGDPpcWB tau_lead if unmem==1 & P5~=1 & ~(T1==1|T2==1|T3==1|T4==1); 
restore; 
/* Asia */
preserve; keep if Asia==1;
probit T0 demaut lpopWB  lGDPpcWB tau_lead if unmem==1 & P5~=1 & ~(T1==1|T2==1|T3==1|T4==1); 
probit T0 demaut lpopWB  lGDPpcWB tau_lead AidGDP if unmem==1 & P5~=1 & ~(T1==1|T2==1|T3==1|T4==1); 
 
probit T0 demaut lpop  lrgdpch  tau_lead if unmem==1 & P5~=1 & ~(T1==1|T2==1|T3==1|T4==1); 
probit T0 demaut lpop  lrgdpch  tau if unmem==1 & P5~=1 & ~(T1==1|T2==1|T3==1|T4==1); 
probit T0 demaut lpop  lrgdpch  tau_lead  Lln_totaid if unmem==1 & P5~=1 & ~(T1==1|T2==1|T3==1|T4==1); 
probit T0 demaut L4delta4demaut lpopWB  lGDPpcWB tau_lead if unmem==1 & P5~=1 & ~(T1==1|T2==1|T3==1|T4==1); 
restore; 
 

/******************* Effects of Membership of Security Council **********************/
#delimit;
sort ccode year;
tsset ccode year;
gen year4SC=.;
replace year4SC=0 if T0==0&T1==0&T2==0&T3==0&T4==0&F1.T0==0&F2.T0==0&F2.T0==0&F3.T0==0&unmem==1 ;
replace year4SC=1 if T0==1&unmem==1 ;
gen year2SC=.;
replace year2SC=0 if T0==0&T1==0&T2==0&T3==0&T4==0&F1.T0==0&F2.T0==0&unmem==1 ;
replace year2SC=1 if T0==1&unmem==1;
capture gen anyaid=1-noaid;
capture gen Lanyaid=L.anyaid;

gen increaseAID=0 if unmem==1 ;
#delimit;
replace increaseAID=1 if ((totaid96-L.totaid96)/L.totaid96>.5)& Lanyaid==1; /* 50% increase in aid in year */
replace increaseAID=1 if ((F1.totaid96-L.totaid96)/L.totaid96>.5)& Lanyaid==1;  /* 50% increase in aid in next year */
replace increaseAID=1 if ((F2.totaid96-L.totaid96)/L.totaid96>.5)& Lanyaid==1; /* 50% increase in aid in year t+2 */
replace increaseAID=1 if (anyaid==1&  Lanyaid==0);  /* start aid in year t*/ 
replace increaseAID=1 if (F.anyaid==1&  Lanyaid==0); /* start aid in t+1 */
replace increaseAID=1 if (F2.anyaid==1&  Lanyaid==0); /* start aid in t+2 */ 

gen UNSCmore=0 if unmem==1;
replace UNSCmore=1 if (increaseAID==1 & year4SC==1); 


#delimit;
/******** Growth in Aid ********/
sort ccode year; tsset ccode year; 
gen USaidgrowth=(F2.totaid96+F1.totaid96)/(2*totaid96);
gen USaidgrowth34=(F3.totaid96+F4.totaid96)/(2*totaid96);
gen LnUSaid12=ln((F2.totaid96+F1.totaid96)/2); 
gen AID12=ln((F2.AID+F1.AID)/2);
gen AIDGDP12=(F2.AidGDP+F1.AidGDP)/2; 
gen OILGDP12=(F2.OIL+F1.OIL)/2; 
gen AIDgrowth=(F2.AID+F1.AID)/(2*AID);
gen AIDgrowth34=(F3.AID+F4.AID)/(2*AID);

gen demautUSaidgrowth=demaut*USaidgrowth;
gen demautUSaidgrowth34=demaut*USaidgrowth34;
gen demautLnUSaid12=demaut*LnUSaid12;
gen demautAID12=demaut*AID12;
gen demautAIDGDP12=demaut*AIDGDP12;
gen demautOILGDP12=demaut*OILGDP12;
gen demautAIDgrowth=demaut*AIDgrowth;
gen demautAIDgrowth34=demaut*AIDgrowth34;

gen demautOIL=demaut*OIL;
gen demautAidGDP=demaut*AidGDP;

/********** military expenditure ********/
tsset ccode year;
gen lnmilex=ln(milex);
 gen lnmilper=ln(milper);
gen Llnmilper=L.lnmilper;
gen L2lnmilper=L2.lnmilper;
gen Llnmilex=L.lnmilex;
gen L2lnmilex=L2.lnmilex;
gen delta4milex=100*(F4.milex-milex)/milex;
gen delta2milex=100*(F2.milex-milex)/milex;
gen delta4milper=100*(F4.milper-milper)/milper;
gen delta2milper=100*(F2.milper-milper)/milper;


/******** leader survival **************/
gen leadchange=LEADERCHANGE if unmem==1;
replace leadchange=(F1.LEADERCHANGE+F2.LEADERCHANGE)/2 if unmem==1 & T0==1;
gen leadchange1234=LEADERCHANGE if unmem==1;
replace leadchange1234=(F1.LEADERCHANGE+F2.LEADERCHANGE+F3.LEADERCHANGE+F4.LEADERCHANGE)/4 if unmem==1 & T0==1;

/*************** GDP *******************/

gen delta4rgdpl=100*(F4.rgdpl-rgdpl)/rgdpl; 
gen delta4GDPpcWB=100*(F4.GDPpcWB- GDPpcWB)/ GDPpcWB; 
gen delta4rgdpch=100*(F4.rgdpch-rgdpch)/rgdpch; 

gen delta2rgdpl=100*(F2.rgdpl-rgdpl)/rgdpl; 
gen delta2GDPpcWB=100*(F2.GDPpcWB- GDPpcWB)/ GDPpcWB; 
gen delta2rgdpch=100*(F2.rgdpch-rgdpch)/rgdpch;

capture gen delta4demaut=F4.demaut-demaut; 
capture gen L4delta4demaut=L4.delta4demaut;
gen delta2demaut=F2.demaut-demaut; 
gen F4demaut=F4.demaut;gen F2demaut=F2.demaut;gen Ldemaut=L1.demaut;
gen demautyear4SC=demaut*year4SC;
gen demautyear2SC=demaut*year2SC;


gen delta4kg=F4.kg-kg; 
gen delta2kg=F2.kg-kg; 

gen delta4tau=F4.tau_lead-tau_lead;
gen delta2tau=F2.tau_lead-tau_lead;
gen F4tau=F4.tau_lead;gen F2tau=F2.tau_lead;gen Ltau=L1.tau_lead;



gen delta4aid=100*(F4.totaid96-totaid96)/totaid96; 
gen delta2aid=100*(F2.totaid96-totaid96)/totaid96;

gen delta4score=F4.score-score; 
gen delta2score=F2.score-score;
gen F4score=F4.score;gen F2score=F2.score;gen Lscore=L1.score;


gen delta2mass=F2.mass-mass; 
gen delta4mass=F4.mass-mass; 

tsset ccode year;
gen L4lpopWB= L4.lpopWB; gen L4lGDPpcWB=L4.lGDPpcWB; 


/* Does Aid increase? */
#delimit;
ttest USaidgrowth, by(year4SC);
ttest USaidgrowth34, by(year4SC);
ttest AIDgrowth, by(year4SC);
ttest AIDgrowth34, by(year4SC);



#delimit;
ttest delta4aid ,by(year4SC);
ttest delta4rgdpl ,by(year4SC);
ttest delta4rgdpch ,by(year4SC);
ttest delta4GDPpcWB ,by(year4SC);
ttest delta4demaut ,by(year4SC);
ttest delta4kg ,by(year4SC);
ttest delta4tau ,by(year4SC);
ttest delta4score ,by(year4SC);
ttest delta4mass ,by(year4SC);
ttest leadchange ,by(year4SC);
ttest leadchange1234 ,by(year4SC);
ttest delta2aid ,by(year2SC);
ttest delta2rgdpl ,by(year2SC);
ttest delta4rgdpch ,by(year2SC);
ttest delta2GDPpcWB ,by(year2SC);
ttest delta2demaut ,by(year2SC);
ttest delta2kg ,by(year2SC);
ttest delta2tau ,by(year2SC);
ttest delta2score ,by(year2SC);
ttest delta2score ,by(year2SC);
ttest leadchange ,by(year2SC); 





/********************** GDP ***********/ 
#delimit;
ttest delta4rgdpch ,by(year4SC);
ttest delta4GDPpcWB ,by(year4SC);
tsset ccode year;
gen L4delta4GDPpcWB=L4.delta4GDPpcWB;
gen L4delta4rgdpch=L4.delta4rgdpch; 
ttest L4delta4rgdpch ,by(year4SC);
ttest L4delta4GDPpcWB ,by(year4SC);
nnmatch delta4GDPpcWB year4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, tc(att)  m(4);; 


xtreg delta4rgdpch year4SC demautyear4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(regyr);
xtreg delta2rgdpch year2SC demautyear2SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(regyr);
xtreg delta4GDPpcWB year4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(regyr);
xtreg delta4GDPpcWB year4SC demautyear4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(regyr);
xtreg delta2GDPpcWB year2SC demautyear2SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(regyr);
#delimit;
xtreg delta4GDPpcWB year4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(ccode);
xtreg delta4GDPpcWB year4SC demautyear4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(ccode);
xtreg delta2GDPpcWB year2SC demautyear2SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(ccode);


/*** prior difference ***/
xtreg L4delta4GDPpcWB year4SC L4demaut L4lpopWB  L4lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(regyr);
xtreg L4delta4GDPpcWB year4SC L4demaut L4lpopWB  L4lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(ccode);



/**** moderating effect of institutions ****/
ttest delta4GDPpcWB if demaut>.7,by(year4SC);
ttest delta4GDPpcWB if demaut<=.7,by(year4SC);
ttest delta2GDPpcWB if demaut>.7,by(year4SC);
ttest delta2GDPpcWB if demaut<=.7,by(year4SC);
nnmatch delta4GDPpcWB year4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1 &demaut>.7, tc(att)  m(4);
nnmatch delta4GDPpcWB year4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1 &demaut<=.7, tc(att)  m(4);

xtreg delta4GDPpcWB year4SC demautyear4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(regyr);




/***************** Democracy *************************/ 
#delimit;
ttest delta4demaut ,by(year4SC);
ttest delta2demaut ,by(year4SC);
tsset ccode year;
capture gen L2delta2demaut=L2.delta2demaut;
capture gen L4delta4demaut=L4.delta4demaut; capture gen F4demaut=F4.demaut; 
ttest L2delta2demaut ,by(year4SC);
nnmatch delta4demaut  year4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1 , tc(att)  m(4);
nnmatch delta2demaut  year4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1 , tc(att)  m(4);

nnmatch delta4demaut  year4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1 &demaut>.7, tc(att)  m(4);
nnmatch delta4demaut  year4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1 &demaut<=.7, tc(att)  m(4);

nnmatch delta2demaut  year4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1 &demaut>.7, tc(att)  m(4);
nnmatch delta2demaut  year4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1 &demaut<=.7, tc(att)  m(4);

#delimit;
xtreg delta4demaut year4SC demautyear4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(regyr);
xtreg delta2demaut year2SC demautyear2SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(regyr);
xtreg delta4demaut year4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(regyr);

xtreg delta4demaut year4SC demautyear4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(ccode);
xtreg delta2demaut year2SC demautyear2SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(ccode);
xtreg delta4demaut year4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(ccode);

#delimit;
xtreg F4demaut demaut year4SC demautyear4SC  lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(regyr);
xtreg F4demaut year4SC  lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(regyr);

xtreg F4demaut demaut year4SC demautyear4SC  lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(ccode);
xtreg F4demaut demaut year4SC  lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(ccode);

xtreg F2demaut demaut year4SC demautyear4SC  lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(ccode);
xtreg F2demaut demaut year4SC  lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(ccode);



#delimit;
/*** prior difference ***/
xtreg L4delta4demaut year4SC  L4demaut L4lpopWB  L4lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(regyr);
xtreg L4delta4demaut year4SC  L4demaut L4lpopWB  L4lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(ccode);
xtreg demaut year4SC  L4demaut L4lpopWB  L4lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(regyr);
xtreg demaut year4SC  L4demaut L4lpopWB  L4lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(ccode);




/****** pressfreedom *************************/ 

#delimit;
ttest delta4score ,by(year4SC);
ttest delta2score ,by(year4SC);
tsset ccode year;
gen L2delta2score=L2.delta2score;
ttest L2delta2score ,by(year4SC);

nnmatch delta4score  year4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1 , tc(att)  m(4);
nnmatch delta2score  year4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1 , tc(att)  m(4);

nnmatch delta4score  year4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1 &demaut>.7, tc(att)  m(4);
nnmatch delta4score  year4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1 &demaut<=.7, tc(att)  m(4);

nnmatch delta2score  year4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1 &demaut>.7, tc(att)  m(4);
nnmatch delta2score  year4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1 &demaut<=.7, tc(att)  m(4);



#delimit;
xtreg delta4score year4SC demautyear4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(regyr);
xtreg delta2score year2SC demautyear2SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(regyr);
xtreg delta4score year4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(regyr);

xtreg delta4score year4SC demautyear4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(ccode);
xtreg delta2score year2SC demautyear2SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(ccode);
xtreg delta4score year4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(ccode);

xtreg F4score score year4SC demautyear4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(regyr);
xtreg F2score score year2SC demautyear2SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(regyr);
xtreg F4score score year4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(regyr);

xtreg F4score score  year4SC demautyear4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(ccode);
xtreg F2score score  year2SC demautyear2SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(ccode);
xtreg F4score score  year4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(ccode);


/*** prior difference ***/
#delimit;
tsset ccode year; capture gen L4delta4score=L4.delta4score;
xtreg L4delta4score year4SC   L4.demaut L4.lpopWB  L4.lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(ccode);

/******* US security alignments **************/ 
#delimit;
ttest delta4tau ,by(year4SC);
ttest delta2tau ,by(year4SC);
tsset ccode year;
capture gen L2delta2tau=L2.delta2tau;capture gen L4delta4tau=L4.delta4tau;
ttest L2delta2tau ,by(year4SC);
ttest L4delta4tau ,by(year4SC);

nnmatch delta4tau  year4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1 , tc(att)  m(4);
nnmatch delta2tau  year4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1 , tc(att)  m(4);

nnmatch delta4tau  year4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1 &demaut>.7, tc(att)  m(4);
nnmatch delta4tau  year4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1 &demaut<=.7, tc(att)  m(4);

nnmatch delta2tau  year4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1 &demaut>.7, tc(att)  m(4);
nnmatch delta2tau  year4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1 &demaut<=.7, tc(att)  m(4);

#delimit;
xtreg delta4tau year4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(regyr);
xtreg delta2tau year2SC   demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(regyr);
xtreg delta4tau year4SC   demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(ccode);
xtreg delta2tau year2SC   demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(ccode);

xtreg F4tau tau year4SC demautyear4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(regyr);
xtreg F2tau tau year2SC demautyear2SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(regyr);
xtreg F4tau tau year4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(regyr);

xtreg F4tau tau  year4SC demautyear4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(ccode);
xtreg F2tau tau  year2SC demautyear2SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(ccode);
xtreg F4tau tau  year4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(ccode);




			/*** prior difference ***/
#delimit;
tsset ccode year; capture gen L4delta4tau=L4.delta4tau;
xtreg L4delta4tau year4SC   L4.demaut L4.lpopWB  L4.lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(ccode);
xtreg L4delta4tau year4SC   demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(ccode);


/******* UN GA voting affinity **************/ 
#delimit;
tsset ccode year;
capture gen delta4s3uni=F4.s3uni-s3uni;
capture gen delta2s3uni=F2.s3uni-s3uni;
capture gen F4s3uni=F4.s3uni;
capture gen F2s3uni=F2.s3uni;
ttest delta4s3uni ,by(year4SC);
ttest delta2s3uni ,by(year4SC);

nnmatch delta4s3uni  year4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1 , tc(att)  m(4);
nnmatch delta2s3uni  year4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1 , tc(att)  m(4);

nnmatch delta4s3uni  year4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1 &demaut>.7, tc(att)  m(4);
nnmatch delta4s3uni  year4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1 &demaut<=.7, tc(att)  m(4);

nnmatch delta2s3uni  year4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1 &demaut>.7, tc(att)  m(4);
nnmatch delta2s3uni  year4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1 &demaut<=.7, tc(att)  m(4);



#delimit;
xtreg delta4s3uni year4SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(regyr);
xtreg delta2s3uni year2SC   demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(regyr);
xtreg delta4s3uni year4SC   demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(ccode);
xtreg delta2s3uni year2SC   demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(ccode);
xtreg F4s3uni s3uni year4SC  demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(regyr);
xtreg F4s3uni s3uni  year2SC demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(regyr);
xtreg F4s3uni s3uni  year4SC  demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(ccode);
xtreg F4s3uni s3uni  year2SC  demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(ccode);




			/*** prior difference ***/
#delimit;
tsset ccode year; capture gen L4delta4s3uni=L4.delta4s3uni;
capture gen L2delta2s3uni=L2.delta2s3uni;
xtreg L4delta4s3uni year4SC   L4demaut L4lpopWB  L4lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(regyr);
xtreg L4delta4s3uni year4SC   L4demaut L4lpopWB  L4lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(ccode);


/************** Transparency International CPI **********/
#delimit;
tsset ccode year; capture gen delta2TI=F2.TI-TI;capture gen delta4TI=F4.TI-TI; capture gen F4TI=F4.TI;
xtreg delta4TI year4SC demautyear4SC demaut lpopWB  lGDPpcWB TI if (year4SC==0|year4SC==1)&unmem==1, fe i(regyr);
xtreg delta4TI year4SC demautyear4SC demaut lpopWB  lGDPpcWB TI if (year4SC==0|year4SC==1)&unmem==1, fe i(ccode);
xtreg delta2TI year4SC demautyear4SC demaut lpopWB  lGDPpcWB TI if (year4SC==0|year4SC==1)&unmem==1, fe i(regyr);
xtreg delta2TI year4SC demautyear4SC demaut lpopWB  lGDPpcWB TI if (year4SC==0|year4SC==1)&unmem==1, fe i(ccode);

#delimit;
tab year, gen(ydum);
areg delta4TI year4SC demautyear4SC demaut lpopWB  lGDPpcWB TI ydum* cotrendreg* egytrend 
	if (year4SC==0|year4SC==1)&unmem==1,  absorb(ccode)  cluster(ccode);
#delimit;
areg F4TI year4SC demautyear4SC demaut lpopWB  lGDPpcWB TI ydum* cotrendreg* egytrend 
	if (year4SC==0|year4SC==1)&unmem==1,  absorb(ccode)  cluster(ccode);

nnmatch delta4TI year4SC demaut lpopWB  lGDPpcWB TI if (year4SC==0|year4SC==1)&unmem==1;

/*********************** tables ****************************/


/*********** table with ccode fixed effects ***************************/
#delimit;

xtreg delta4GDPpcWB year4SC demautyear4SC  demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(ccode);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 delta4GDPpcWB year4SC demautyear4SC  demaut lpopWB  lGDPpcWB 
	using C:\al2009\UN_securitycouncil\writeup\effect.out,replace  addstat(R2, e(r2), Log Lik, e(ll));
#delimit;
xtreg delta2GDPpcWB year4SC demautyear4SC  demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(ccode);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 delta2GDPpcWB year4SC demautyear4SC  demaut lpopWB  lGDPpcWB 
	using C:\al2009\UN_securitycouncil\writeup\effect.out,append  addstat(R2, e(r2), Log Lik, e(ll));

#delimit;
xtreg F4demaut year4SC demautyear4SC  demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(ccode);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effect.out,append  addstat(R2, e(r2), Log Lik, e(ll));
#delimit;
xtreg F2demaut year4SC demautyear4SC  demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(ccode);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effect.out,append  addstat(R2, e(r2), Log Lik, e(ll));

#delimit;
xtreg F4score year4SC demautyear4SC  demaut lpopWB  lGDPpcWB score if (year4SC==0|year4SC==1)&unmem==1, fe i(ccode);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effect.out,append  addstat(R2, e(r2), Log Lik, e(ll));
#delimit;
xtreg F2score year4SC demautyear4SC  demaut lpopWB  lGDPpcWB score if (year4SC==0|year4SC==1)&unmem==1, fe i(ccode);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effect.out,append  addstat(R2, e(r2), Log Lik, e(ll));

#delimit;
xtreg F4tau tau year4SC demautyear4SC  demaut lpopWB  lGDPpcWB  if (year4SC==0|year4SC==1)&unmem==1 & ccode~=316 , fe i(ccode);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effect.out,append  addstat(R2, e(r2), Log Lik, e(ll));
#delimit;
xtreg F2tau tau year4SC demautyear4SC  demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(ccode);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effect.out,append  addstat(R2, e(r2), Log Lik, e(ll))  ;
#delimit;
xtreg F4s3uni s3uni year4SC demautyear4SC  demaut lpopWB  lGDPpcWB  if (year4SC==0|year4SC==1)&unmem==1, fe i(ccode);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effect.out,append  addstat(R2, e(r2), Log Lik, e(ll));
#delimit;
xtreg F2s3uni s3uni year4SC demautyear4SC  demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(ccode);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effect.out,append  addstat(R2, e(r2), Log Lik, e(ll)) word;


/*********** table with regyr fixed effects ***************************/
#delimit;
xtreg delta4GDPpcWB year4SC demautyear4SC  demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(regyr);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 delta4GDPpcWB year4SC demautyear4SC  demaut lpopWB  lGDPpcWB 
	using C:\al2009\UN_securitycouncil\writeup\effect_regyr.out,replace  addstat(R2, e(r2), Log Lik, e(ll));
#delimit;
xtreg delta2GDPpcWB year4SC demautyear4SC  demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(regyr);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 delta2GDPpcWB year4SC demautyear4SC  demaut lpopWB  lGDPpcWB 
	using C:\al2009\UN_securitycouncil\writeup\effect_regyr.out,append  addstat(R2, e(r2), Log Lik, e(ll));

#delimit;
xtreg F4demaut year4SC demautyear4SC  demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(regyr);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effect_regyr.out,append  addstat(R2, e(r2), Log Lik, e(ll));
#delimit;
xtreg F2demaut year4SC demautyear4SC  demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(regyr);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effect_regyr.out,append  addstat(R2, e(r2), Log Lik, e(ll));

#delimit;
xtreg F4score year4SC demautyear4SC  demaut lpopWB  lGDPpcWB score if (year4SC==0|year4SC==1)&unmem==1, fe i(regyr);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effect_regyr.out,append  addstat(R2, e(r2), Log Lik, e(ll));
#delimit;
xtreg F2score year4SC demautyear4SC  demaut lpopWB  lGDPpcWB score if (year4SC==0|year4SC==1)&unmem==1, fe i(regyr);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effect_regyr.out,append  addstat(R2, e(r2), Log Lik, e(ll));

#delimit;
xtreg F4tau tau year4SC demautyear4SC  demaut lpopWB  lGDPpcWB  if (year4SC==0|year4SC==1)&unmem==1, fe i(regyr);

test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effect_regyr.out,append  addstat(R2, e(r2), Log Lik, e(ll));
#delimit;
xtreg F2tau tau year4SC demautyear4SC  demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(regyr);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effect_regyr.out,append  addstat(R2, e(r2), Log Lik, e(ll))  ;
#delimit;
xtreg F4s3uni s3uni year4SC demautyear4SC  demaut lpopWB  lGDPpcWB  if (year4SC==0|year4SC==1)&unmem==1, fe i(regyr);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effect_regyr.out,append  addstat(R2, e(r2), Log Lik, e(ll));
#delimit;
xtreg F2s3uni s3uni year4SC demautyear4SC  demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1, fe i(regyr);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effect_regyr.out,append  addstat(R2, e(r2), Log Lik, e(ll))  word;



/********* ccode fixed effect with region year quadtics **********/

#delimit;

areg delta4GDPpcWB year4SC demautyear4SC demaut lpopWB  lGDPpcWB score tau s3uni ydum* cotrendreg* egytrend 
	if (year4SC==0|year4SC==1)&unmem==1,  absorb(ccode)  cluster(ccode);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effect_all_fe.out,replace  addstat(R2, e(r2), Log Lik, e(ll));
#delimit;
areg delta4GDPpcWB year4SC demautyear4SC demaut lpopWB  lGDPpcWB ydum* cotrendreg* egytrend 
	if (year4SC==0|year4SC==1)&unmem==1,  absorb(ccode)  cluster(ccode);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effect_all_fe.out,append  addstat(R2, e(r2), Log Lik, e(ll));
#delimit;
areg delta2GDPpcWB year4SC demautyear4SC demaut lpopWB  lGDPpcWB ydum* cotrendreg* egytrend 
	if (year4SC==0|year4SC==1)&unmem==1,  absorb(ccode)  cluster(ccode);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effect_all_fe.out,append  addstat(R2, e(r2), Log Lik, e(ll));
#delimit;
areg F4demaut year4SC demautyear4SC demaut lpopWB  lGDPpcWB ydum* cotrendreg* egytrend 
	if (year4SC==0|year4SC==1)&unmem==1,  absorb(ccode)  cluster(ccode);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effect_all_fe.out,append  addstat(R2, e(r2), Log Lik, e(ll));
#delimit;
areg F2demaut year4SC demautyear4SC demaut lpopWB  lGDPpcWB ydum* cotrendreg* egytrend 
	if (year4SC==0|year4SC==1)&unmem==1,  absorb(ccode)  cluster(ccode);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effect_all_fe.out,append  addstat(R2, e(r2), Log Lik, e(ll));

#delimit;
areg F4score year4SC demautyear4SC demaut lpopWB  lGDPpcWB score ydum* cotrendreg* egytrend 
	if (year4SC==0|year4SC==1)&unmem==1,  absorb(ccode)  cluster(ccode);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effect_all_fe.out,append  addstat(R2, e(r2), Log Lik, e(ll));

#delimit;
areg F2score year4SC demautyear4SC demaut lpopWB  lGDPpcWB score ydum* cotrendreg* egytrend 
	if (year4SC==0|year4SC==1)&unmem==1,  absorb(ccode)  cluster(ccode);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effect_all_fe.out,append  addstat(R2, e(r2), Log Lik, e(ll));


#delimit;
areg F4tau year4SC demautyear4SC demaut lpopWB  lGDPpcWB tau ydum* cotrendreg* egytrend 
	if (year4SC==0|year4SC==1)&unmem==1,  absorb(ccode)  cluster(ccode);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effect_all_fe.out,append  addstat(R2, e(r2), Log Lik, e(ll));
#delimit;
areg F2tau year4SC demautyear4SC demaut lpopWB  lGDPpcWB tau ydum* cotrendreg* egytrend 
	if (year4SC==0|year4SC==1)&unmem==1,  absorb(ccode)  cluster(ccode);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effect_all_fe.out,append  addstat(R2, e(r2), Log Lik, e(ll)) word;

#delimit;
areg F4s3uni year4SC demautyear4SC demaut lpopWB  lGDPpcWB s3uni ydum* cotrendreg* egytrend 
	if (year4SC==0|year4SC==1)&unmem==1,  absorb(ccode)  cluster(ccode);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effect_all_fe.out,append  addstat(R2, e(r2), Log Lik, e(ll)) word;



/* robustness test that results not due to African cases with decline in democracy prior to entry to UNSC */ 
#delimit;
capture  gen precoup=( L4delta4demaut<-0.1 & T0==1);
/* repeat tables excluding these cases */ 

/*********** table with ccode fixed effects ***************************/
#delimit;

xtreg delta4GDPpcWB year4SC demautyear4SC  demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1&precoup==0, fe i(ccode);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 delta4GDPpcWB year4SC demautyear4SC  demaut lpopWB  lGDPpcWB 
	using C:\al2009\UN_securitycouncil\writeup\effectR.out,replace  addstat(R2, e(r2), Log Lik, e(ll));
#delimit;
xtreg delta2GDPpcWB year4SC demautyear4SC  demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1&precoup==0, fe i(ccode);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 delta2GDPpcWB year4SC demautyear4SC  demaut lpopWB  lGDPpcWB 
	using C:\al2009\UN_securitycouncil\writeup\effectR.out,append  addstat(R2, e(r2), Log Lik, e(ll));

#delimit;
xtreg F4demaut year4SC demautyear4SC  demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1&precoup==0, fe i(ccode);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effectR.out,append  addstat(R2, e(r2), Log Lik, e(ll));
#delimit;
xtreg F2demaut year4SC demautyear4SC  demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1&precoup==0, fe i(ccode);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effectR.out,append  addstat(R2, e(r2), Log Lik, e(ll));

#delimit;
xtreg F4score year4SC demautyear4SC  demaut lpopWB  lGDPpcWB score if (year4SC==0|year4SC==1)&unmem==1&precoup==0, fe i(ccode);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effectR.out,append  addstat(R2, e(r2), Log Lik, e(ll));
#delimit;
xtreg F2score year4SC demautyear4SC  demaut lpopWB  lGDPpcWB score if (year4SC==0|year4SC==1)&unmem==1&precoup==0, fe i(ccode);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effectR.out,append  addstat(R2, e(r2), Log Lik, e(ll));

#delimit;
xtreg F4tau tau year4SC demautyear4SC  demaut lpopWB  lGDPpcWB  if (year4SC==0|year4SC==1)&unmem==1 & ccode~=316 &precoup==0, fe i(ccode);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effectR.out,append  addstat(R2, e(r2), Log Lik, e(ll));
#delimit;
xtreg F2tau tau year4SC demautyear4SC  demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1 &precoup==0, fe i(ccode);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effectR.out,append  addstat(R2, e(r2), Log Lik, e(ll))  ;
#delimit;
xtreg F4s3uni s3uni year4SC demautyear4SC  demaut lpopWB  lGDPpcWB  if (year4SC==0|year4SC==1)&unmem==1&precoup==0, fe i(ccode);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effectR.out,append  addstat(R2, e(r2), Log Lik, e(ll));
#delimit;
xtreg F2s3uni s3uni year4SC demautyear4SC  demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1&precoup==0, fe i(ccode);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effectR.out,append  addstat(R2, e(r2), Log Lik, e(ll)) word;


/*********** table with regyr fixed effects ***************************/
#delimit;
xtreg delta4GDPpcWB year4SC demautyear4SC  demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1&precoup==0, fe i(regyr);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 delta4GDPpcWB year4SC demautyear4SC  demaut lpopWB  lGDPpcWB 
	using C:\al2009\UN_securitycouncil\writeup\effect_regyrR.out,replace  addstat(R2, e(r2), Log Lik, e(ll));
#delimit;
xtreg delta2GDPpcWB year4SC demautyear4SC  demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1&precoup==0, fe i(regyr);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 delta2GDPpcWB year4SC demautyear4SC  demaut lpopWB  lGDPpcWB 
	using C:\al2009\UN_securitycouncil\writeup\effect_regyrR.out,append  addstat(R2, e(r2), Log Lik, e(ll));

#delimit;
xtreg F4demaut year4SC demautyear4SC  demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1&precoup==0, fe i(regyr);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effect_regyrR.out,append  addstat(R2, e(r2), Log Lik, e(ll));
#delimit;
xtreg F2demaut year4SC demautyear4SC  demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1&precoup==0, fe i(regyr);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effect_regyrR.out,append  addstat(R2, e(r2), Log Lik, e(ll));

#delimit;
xtreg F4score year4SC demautyear4SC  demaut lpopWB  lGDPpcWB score if (year4SC==0|year4SC==1)&unmem==1&precoup==0, fe i(regyr);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effect_regyrR.out,append  addstat(R2, e(r2), Log Lik, e(ll));
#delimit;
xtreg F2score year4SC demautyear4SC  demaut lpopWB  lGDPpcWB score if (year4SC==0|year4SC==1)&unmem==1&precoup==0, fe i(regyr);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effect_regyrR.out,append  addstat(R2, e(r2), Log Lik, e(ll));

#delimit;
xtreg F4tau tau year4SC demautyear4SC  demaut lpopWB  lGDPpcWB  if (year4SC==0|year4SC==1)&unmem==1&precoup==0, fe i(regyr);

test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effect_regyrR.out,append  addstat(R2, e(r2), Log Lik, e(ll));
#delimit;
xtreg F2tau tau year4SC demautyear4SC  demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1&precoup==0, fe i(regyr);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effect_regyrR.out,append  addstat(R2, e(r2), Log Lik, e(ll))  ;
#delimit;
xtreg F4s3uni s3uni year4SC demautyear4SC  demaut lpopWB  lGDPpcWB  if (year4SC==0|year4SC==1)&unmem==1&precoup==0, fe i(regyr);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effect_regyrR.out,append  addstat(R2, e(r2), Log Lik, e(ll));
#delimit;
xtreg F2s3uni s3uni year4SC demautyear4SC  demaut lpopWB  lGDPpcWB if (year4SC==0|year4SC==1)&unmem==1&precoup==0, fe i(regyr);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effect_regyrR.out,append  addstat(R2, e(r2), Log Lik, e(ll))  word;



/********* ccode fixed effect with region year quadtics **********/

#delimit;
areg delta4GDPpcWB year4SC demautyear4SC demaut lpopWB  lGDPpcWB score tau s3uni ydum* cotrendreg* egytrend 
	if (year4SC==0|year4SC==1)&unmem==1&precoup==0,  absorb(ccode)  cluster(ccode);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effect_all_feR.out,replace  addstat(R2, e(r2), Log Lik, e(ll));
#delimit;
areg delta4GDPpcWB year4SC demautyear4SC demaut lpopWB  lGDPpcWB ydum* cotrendreg* egytrend 
	if (year4SC==0|year4SC==1)&unmem==1&precoup==0,  absorb(ccode)  cluster(ccode);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effect_all_feR.out,append  addstat(R2, e(r2), Log Lik, e(ll));
#delimit;
areg delta2GDPpcWB year4SC demautyear4SC demaut lpopWB  lGDPpcWB ydum* cotrendreg* egytrend 
	if (year4SC==0|year4SC==1)&unmem==1&precoup==0,  absorb(ccode)  cluster(ccode);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effect_all_feR.out,append  addstat(R2, e(r2), Log Lik, e(ll));
#delimit;
areg F4demaut year4SC demautyear4SC demaut lpopWB  lGDPpcWB ydum* cotrendreg* egytrend 
	if (year4SC==0|year4SC==1)&unmem==1&precoup==0,  absorb(ccode)  cluster(ccode);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effect_all_feR.out,append  addstat(R2, e(r2), Log Lik, e(ll));
#delimit;
areg F2demaut year4SC demautyear4SC demaut lpopWB  lGDPpcWB ydum* cotrendreg* egytrend 
	if (year4SC==0|year4SC==1)&unmem==1&precoup==0,  absorb(ccode)  cluster(ccode);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effect_all_feR.out,append  addstat(R2, e(r2), Log Lik, e(ll));

#delimit;
areg F4score year4SC demautyear4SC demaut lpopWB  lGDPpcWB score ydum* cotrendreg* egytrend 
	if (year4SC==0|year4SC==1)&unmem==1&precoup==0,  absorb(ccode)  cluster(ccode);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effect_all_feR.out,append  addstat(R2, e(r2), Log Lik, e(ll));

#delimit;
areg F2score year4SC demautyear4SC demaut lpopWB  lGDPpcWB score ydum* cotrendreg* egytrend 
	if (year4SC==0|year4SC==1)&unmem==1&precoup==0,  absorb(ccode)  cluster(ccode);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effect_all_feR.out,append  addstat(R2, e(r2), Log Lik, e(ll));


#delimit;
areg F4tau year4SC demautyear4SC demaut lpopWB  lGDPpcWB tau ydum* cotrendreg* egytrend 
	if (year4SC==0|year4SC==1)&unmem==1&precoup==0,  absorb(ccode)  cluster(ccode);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effect_all_feR.out,append  addstat(R2, e(r2), Log Lik, e(ll));
#delimit;
areg F2tau year4SC demautyear4SC demaut lpopWB  lGDPpcWB tau ydum* cotrendreg* egytrend 
	if (year4SC==0|year4SC==1)&unmem==1&precoup==0,  absorb(ccode)  cluster(ccode);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effect_all_feR.out,append  addstat(R2, e(r2), Log Lik, e(ll)) word;

#delimit;
areg F4s3uni year4SC demautyear4SC demaut lpopWB  lGDPpcWB s3uni ydum* cotrendreg* egytrend 
	if (year4SC==0|year4SC==1)&unmem==1&precoup==0,  absorb(ccode)  cluster(ccode);
test year4SC demautyear4SC  ; test year4SC+ demautyear4SC  ==0;
outreg2 
	using C:\al2009\UN_securitycouncil\writeup\effect_all_feR.out,append  addstat(R2, e(r2), Log Lik, e(ll)) word;

