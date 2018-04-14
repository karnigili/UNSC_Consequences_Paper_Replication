/* Alastair Smith January 2009 --- how does UN SC affect amount of aid */ 
/* use data generated in C:\al2009\UN_securitycouncil\Werker_UNseatDATA\Alastair_organize_UN_aid_data.do */ 

/*************************************** Describe how level aid changes in UNSC *****************/ 

#delimit; 
clear; set mem 100m;
use C:\al2009\UN_securitycouncil\Werker_UNseatDATA\USaid_UNSC_selectorate_data_temp.dta; 
use USaid_UNSC_selectorate_data_replication.dta; 
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


capture log close;
log using HowMuchAidAndUNSCreplication,replace; 
sort ccode year;

drop if (ccode==ccode[_n-1]& year==year[_n-1]); 
tsset ccode year; 
gen Lln_totaid=L.ln_totaid96;
gen Ltotaid96=L.totaid96;
gen lntotaid96=ln(totaid96);

#delimit;
sort ccode year;
tsset ccode year;
gen year4SC=.;
replace year4SC=0 if T0==0&T1==0&T2==0&T3==0&T4==0&F1.T0==0&F2.T0==0&F2.T0==0&F3.T0==0&unmem==1 ;
replace year4SC=1 if T0==1&unmem==1 ;

/**** T0 corresponds to year of election to SC ****/
/** compare average aid in years t-1 and t-2 with years t1 and t2  to find average increase in aid **/ 
tsset ccode year; 
gen anyaid=1-noaid;
gen F1anyaid=F1.anyaid; 
gen L1anyaid=L1.anyaid; 
gen anyaidWB=0 if AidGDP~=.;
replace anyaidWB=1 if AidGDP>0;
gen FanyaidWB=F1.anyaidWB;
gen LanyaidWB=L1.anyaidWB;
gen LAidGDP=L.AidGDP;
gen Ltau=L.tau; 


gen noaidtminus2=L2.noaid; gen noaidtminus1=L1.noaid; gen noaidtplus1=F1.noaid; gen noaidtplus2=F2.noaid;
gen noaidbefore=0; replace noaidbefore=1 if noaidtminus1==1|noaidtminus2==1|noaid==1; 
gen noaidduring=0; replace noaidduring=1 if noaidtplus1==1|noaidtplus2==1; 
gen aidbefore=(L1.totaid96+L2.totaid96) /2; 
gen aidduring=(F1.totaid96+F2.totaid96) /2; 
gen aidafter=(F3.totaid96+F4.totaid96) /2;

gen F12lnTotAid96=ln((F1.totaid96+F2.totaid96) /2); 
gen lnaidduring=ln(aidduring);
gen UNSC=scmem; replace UNSC=1 if T0==1;
gen  UNSCimptyr= UNSC*imptyr;
#delimit;
/* Does UNSC increase the likelihood of getting aid */ 
tab  UNSC anyaid if hic==0, all row;
tab  UNSC anyaid if  L1anyaid==0&hic==0, all row;
tab  UNSC anyaid if  L1anyaid==1&hic==0, all row; 



probit  anyaid L1anyaid  demaut scmem lpop  lrgdpch  logdstab if hic==0;
probit  anyaid L1anyaid  demaut UNSC  lpop  lrgdpch  logdstab if hic==0;
probit  anyaid L1anyaid  demaut UNSC imptyr UNSCimptyr lpop  lrgdpch  logdstab if hic==0;
probit  anyaid L1anyaid  demaut T0 T1 T2 T3 T4 lpop  lrgdpch  logdstab if hic==0;
test T0 T1 T2; test T0+T1+T2==0;
probit  anyaid L1anyaid  demaut T0 T1 T2 T3 T4 lpop  lrgdpch  logdstab tau Ltau bigwar if hic==0 & L1anyaid==0;
test T0 T1 T2; test T0+T1+T2==0;
probit  anyaid L1anyaid  demaut T0 T1 T2 T3 T4 lpop  lrgdpch  logdstab tau Ltau bigwar if hic==0 & L1anyaid==1;
test T0 T1 T2; test T0+T1+T2==0;


probit  anyaid L1anyaid  demaut lpop  lrgdpch  logdstab if hic==0& UNSC==1; /* looking just at the nations entering UNSC we see that the factors that matter in general -- distance wealth size etc  have no effect */ 
#delimit;
li ccode name year Ltotaid96 T0 T1 T2 if UNSC ==1 &  L1anyaid  ==1 & anyaid ==0 &hic==0; 
li ccode name year totaid96 T0 T1 T2 if UNSC ==1 &  L1anyaid  ==0 & anyaid ==1 &hic==0; 
nnmatch anyaid UNSC  demaut lpop  lrgdpch  logdstab if hic==0 & L1anyaid ==0, tc(att)  m(4);
nnmatch anyaid UNSC  L1anyaid demaut lpop  lrgdpch  logdstab if hic==0 ,tc(att)  m(4);
nnmatch anyaid UNSC  demaut lpop  lrgdpch  logdstab if hic==0 & L1anyaid ==1, tc(att)  ;
nnmatch anyaid T3  demaut lpop  lrgdpch  logdstab if hic==0 & L1anyaid ==1, tc(att)  ;
gen endUNSC=0 if unmem==1; replace endUNSC=1 if T3==1|T4==1; 
#delimit;
nnmatch anyaid endUNSC  demaut lpop  lrgdpch  logdstab if hic==0 & L1anyaid ==1, tc(att) m(4) ;
#delimit;
nnmatch anyaid endUNSC  L1anyaid demaut lpop  lrgdpch  logdstab if hic==0 ,tc(att) m(4) ;

/* repeat with WB WDI measure of variables */ 
#delimit;
probit  anyaid L1anyaid  demaut scmem lpopWB  lGDPpcWB  logdstab if hic==0;
probit  anyaid L1anyaid  demaut UNSC  lpopWB  lGDPpcWB  logdstab if hic==0;
probit  anyaid L1anyaid  demaut UNSC imptyr UNSCimptyr lpopWB  lGDPpcWB  logdstab if hic==0;
probit  anyaid L1anyaid  demaut T0 T1 T2 T3 T4 lpopWB  lGDPpcWB  logdstab if hic==0;
test T0 T1 T2; test T0+T1+T2==0;
probit  anyaid L1anyaid  demaut T0 T1 T2 T3 T4 lpopWB  lGDPpcWB  logdstab tau Ltau bigwar if hic==0 & L1anyaid==0;
test T0 T1 T2; test T0+T1+T2==0;
probit  anyaid L1anyaid  demaut T0 T1 T2 T3 T4 lpopWB  lGDPpcWB  logdstab tau Ltau bigwar if hic==0 & L1anyaid==1;



probit  anyaid L1anyaid  demaut lpopWB  lGDPpcWB  logdstab if hic==0& UNSC==1; /* looking just at the nations entering UNSC we see that the factors that matter in general -- distance wealth size etc  have no effect */ 


nnmatch anyaid UNSC  demaut lpopWB  lGDPpcWB  logdstab if hic==0 & L1anyaid ==0, tc(att)  m(4);
nnmatch anyaid UNSC  L1anyaid demaut lpopWB  lGDPpcWB  logdstab if hic==0 ,tc(att)  m(4);
nnmatch anyaid UNSC  demaut lpopWB  lGDPpcWB  logdstab if hic==0 & L1anyaid ==1, tc(att)  ;
nnmatch anyaid T3  demaut lpopWB  lGDPpcWB  logdstab if hic==0 & L1anyaid ==1, tc(att)  ;
capture gen endUNSC=0 if unmem==1; replace endUNSC=1 if T3==1|T4==1; 
#delimit;
nnmatch anyaid endUNSC  demaut lpopWB  lGDPpcWB  logdstab if hic==0 & L1anyaid ==1, tc(att) m(4) ;
#delimit;
nnmatch anyaid endUNSC  L1anyaid demaut lpopWB  lGDPpcWB  logdstab if hic==0 ,tc(att) m(4) ;

/* table */ 
#delimit;

ttest anyaid if L1anyaid   ==0 &hic==0, by(UNSC);



probit  anyaid L1anyaid UNSC scmem T0 endUNSC demaut lpopWB  lGDPpcWB  logdstab  if hic==0;


probit  anyaid L1anyaid   UNSC demaut lpopWB  lGDPpcWB  logdstab if hic==0;

probit  anyaid L1anyaid  scmem T0 endUNSC demaut lpopWB  lGDPpcWB  logdstab if hic==0;
probit  anyaid L1anyaid   UNSC lpopWB demaut lGDPpcWB  logdstab if hic==0& L1anyaid==0;
probit  anyaid L1anyaid   UNSC lpopWB demaut lGDPpcWB  logdstab if hic==0& L1anyaid==1;
probit  anyaid L1anyaid UNSC  endUNSC demaut lpopWB  lGDPpcWB  logdstab if hic==0& L1anyaid==1;



#delimit;
ttest anyaid if L1anyaid   ==0 &hic==0, by(UNSC);
capture gen endUNSC=0 if unmem==1; replace endUNSC=1 if T3==1|T4==1;
nnmatch anyaid scmem  demaut lpopWB  lGDPpcWB    if hic==0 ,tc(att) m(4) ;
nnmatch anyaid scmem  demaut lpopWB  lGDPpcWB    if hic==0 & L1anyaid ==0, tc(att) m(4) ;
nnmatch anyaid endUNSC  demaut lpopWB  lGDPpcWB    if hic==0 & L1anyaid ==1, tc(att) m(4) ;




/* Amount of USAID */ 
#delimit; 
xtreg lntotaid96 Lln_totaid demaut  T0 T1 T2 T3 T4 lpop  lrgdpch  logdstab tau bigwar 
						if hic==0 , fe i(regyr);
test T0 T1 T2; test T0+T1+T2==0;
#delimit; 
xtreg lntotaid96 Lln_totaid demaut  T0 T1 T2 T3 T4 lpop  lrgdpch  logdstab tau bigwar 
						if hic==0 & L1anyaid==1 &anyaid==1 , fe i(regyr);
test T0 T1 T2; test T0+T1+T2==0;
#delimit; 
xtreg lntotaid96 Lln_totaid demaut UNSC  imptyr UNSCimptyr lpop  lrgdpch  logdstab tau bigwar 
						if hic==0 & L1anyaid==1 &anyaid==1 , fe i(regyr);

/***** overall level of aid -- not just US aid */ 
#delimit; 
xtreg AidGDP  LAidGDP demaut UNSC  lpop  lrgdpch  logdstab tau bigwar 
						if hic==0 & L1anyaid==1 &anyaid==1 , fe i(regyr);
/* repeat with WB WDI measures */
#delimit; 
xtreg lntotaid96 Lln_totaid demaut  T0 T1 T2 T3 T4 lpopWB  lGDPpcWB    logdstab tau bigwar 
						if hic==0 , fe i(regyr);
test T0 T1 T2; test T0+T1+T2==0;
#delimit; 
xtreg lntotaid96 Lln_totaid demaut  T0 T1 T2 T3 T4 lpopWB  lGDPpcWB    logdstab tau bigwar 
						if hic==0 & L1anyaid==1 &anyaid==1 , fe i(regyr);
test T0 T1 T2; test T0+T1+T2==0;
#delimit; 
xtreg lntotaid96 Lln_totaid demaut UNSC  imptyr UNSCimptyr lpopWB  lGDPpcWB    logdstab tau bigwar 
						if hic==0 & L1anyaid==1 &anyaid==1 , fe i(regyr);

/***** overall level of aid -- not just US aid */ 
#delimit; 
xtreg AidGDP  LAidGDP demaut UNSC  lpopWB  lGDPpcWB    logdstab tau bigwar 
						if hic==0 & L1anyaid==1 &anyaid==1 , fe i(regyr);


/*************************************************************/
/* Now look at analyses in the spirit of Kuziemko and Werker */ 
/*************************************************************/
#delimit;
capture gen endUNSC=0 if unmem==1; replace endUNSC=1 if T3==1|T4==1;
reg ln_totaid96 scmem if hic==0  ;
reg ln_totaid96 scmem if hic==0  & L1anyaid==1 &anyaid==1  ;
reg ln_totaid96 UNSC if hic==0  ;
reg ln_totaid96 UNSC if hic==0  & L1anyaid==1 &anyaid==1  ;
areg ln_totaid96 scmem ydum* cotrendreg* egytrend if hic==0  , absorb(ccode)  cluster(ccode);
areg ln_totaid96 scmem ydum* cotrendreg* egytrend if hic==0 & L1anyaid==1 &anyaid==1  , absorb(ccode)  cluster(ccode);
areg ln_totaid96 UNSC endUNSC ydum* cotrendreg* egytrend if hic==0 , absorb(ccode)  cluster(ccode);
areg ln_totaid96 UNSC endUNSC ydum* cotrendreg* egytrend if hic==0 & L1anyaid==1 &anyaid==1  , absorb(ccode)  cluster(ccode);
areg  ln_totaid96 scnyt1 scnyt2 scnyt3 ydum* cotrendreg* egytrend if insample & unmem==1 , absorb(ccode)  cluster(ccode) ;
areg  ln_totaid96 scnyt1 scnyt2 scnyt3 ydum* cotrendreg* egytrend if insample & unmem==1 & L1anyaid==1 &anyaid==1  , absorb(ccode)  cluster(ccode) ;
areg ln_totaid96 UNSC endUNSC UNSC bigwar  demaut lGDPpcWB    logdstab tau bigwar ydum* cotrendreg* egytrend if hic==0 , absorb(ccode)  cluster(ccode);
areg ln_totaid96 UNSC endUNSC UNSC bigwar  demaut lGDPpcWB    logdstab tau bigwar ydum* cotrendreg* egytrend if hic==0 & L1anyaid==1 &anyaid==1  , absorb(ccode)  cluster(ccode);

/*** now control for prior level of aid */ 
reg ln_totaid96 Lln_totaid  scmem if hic==0  ;
reg ln_totaid96 Lln_totaid  scmem if hic==0  & L1anyaid==1 &anyaid==1  ;
reg ln_totaid96 Lln_totaid  UNSC if hic==0  ;
reg ln_totaid96 Lln_totaid  UNSC if hic==0  & L1anyaid==1 &anyaid==1  ;
areg ln_totaid96 Lln_totaid  scmem ydum* cotrendreg* egytrend if hic==0  , absorb(ccode)  cluster(ccode);
areg ln_totaid96 Lln_totaid  scmem ydum* cotrendreg* egytrend if hic==0 & L1anyaid==1 &anyaid==1  , absorb(ccode)  cluster(ccode);
areg ln_totaid96 Lln_totaid UNSC endUNSC ydum* cotrendreg* egytrend if hic==0 , absorb(ccode)  cluster(ccode);
areg ln_totaid96 Lln_totaid UNSC endUNSC ydum* cotrendreg* egytrend if hic==0 & L1anyaid==1 &anyaid==1  , absorb(ccode)  cluster(ccode);
areg  ln_totaid96 Lln_totaid scnyt1 scnyt2 scnyt3 ydum* cotrendreg* egytrend if insample & unmem==1 , absorb(ccode)  cluster(ccode) ;
areg  ln_totaid96 Lln_totaid scnyt1 scnyt2 scnyt3 ydum* cotrendreg* egytrend if insample & unmem==1 & L1anyaid==1 &anyaid==1  , absorb(ccode)  cluster(ccode) ;
areg ln_totaid96 Lln_totaid UNSC endUNSC UNSC bigwar  demaut lGDPpcWB    logdstab tau bigwar ydum* cotrendreg* egytrend if hic==0 , absorb(ccode)  cluster(ccode);
areg ln_totaid96 Lln_totaid UNSC endUNSC UNSC bigwar  demaut lGDPpcWB    logdstab tau bigwar ydum* cotrendreg* egytrend if hic==0 & L1anyaid==1 &anyaid==1  , absorb(ccode)  cluster(ccode);









