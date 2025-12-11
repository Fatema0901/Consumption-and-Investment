
  ___  ____  ____  ____  ____ ®
 /__    /   ____/   /   ____/      StataNow 19.5
___/   /   /___/   /   /___/       BE—Basic Edition

 Statistics and Data Science       Copyright 1985-2025 StataCorp LLC
                                   StataCorp
                                   4905 Lakeway Drive
                                   College Station, Texas 77845 USA
                                   800-782-8272        https://www.stata.com
                                   979-696-4600        service@stata.com

Stata license: Single-user , expiring 15 Mar 2026
Serial number: 301909014077
  Licensed to: Fatema
               TTU

Notes:
      1. Unicode is supported; see help unicode_advice.

. clear

. cd "C:\Users\fatem\OneDrive\Desktop\Fall 2025\Consumption\Term Paper"
C:\Users\fatem\OneDrive\Desktop\Fall 2025\Consumption\Term Paper

. set fredkey 38f24e70221b0a608d3e539bd2a46833

. *********REAL PERSONAL CONSUMPTION EXPENDITURES*********

. 
. import fred PCECC96, clear

Summary
-----------------------------------------------------------------------------------------------------------------------
Series ID                    Nobs    Date range                Frequency
-----------------------------------------------------------------------------------------------------------------------
PCECC96                      314     1947-01-01 to 2025-04-01  Quarterly
-----------------------------------------------------------------------------------------------------------------------
# of series imported: 1
   highest frequency: Quarterly
    lowest frequency: Quarterly

. 
. gen daten = qofd(daten)
.
. format daten %tq

. 
. save pce_clean.dta, replace
(file pce_clean.dta not found)
file pce_clean.dta saved

. 
. 
. 
. *********PERSONAL SAVING RATE*********

. 
. import fred PSAVERT, clear

Summary
-----------------------------------------------------------------------------------------------------------------------
Series ID                    Nobs    Date range                Frequency
-----------------------------------------------------------------------------------------------------------------------
PSAVERT                      801     1959-01-01 to 2025-09-01  Monthly
-----------------------------------------------------------------------------------------------------------------------
# of series imported: 1
   highest frequency: Monthly
    lowest frequency: Monthly

. 
. gen daten = mofd(daten)
.
. format daten %tm

. 
. save psavert_clean.dta, replace
(file psavert_clean.dta not found)
file psavert_clean.dta saved

. 
. 
. 
. *********ECONOMIC POLICY UNCERTAINTY INDEX*********

. 
. import fred USEPUINDXD, clear

Summary
-----------------------------------------------------------------------------------------------------------------------
Series ID                    Nobs    Date range                Frequency
-----------------------------------------------------------------------------------------------------------------------
USEPUINDXD                   14948   1985-01-01 to 2025-12-04  Daily, 7-Day
-----------------------------------------------------------------------------------------------------------------------
# of series imported: 1
   highest frequency: Daily
    lowest frequency: Daily

. 
. gen daten = mofd(daten)
.
. format daten %tm

. 
. collapse (mean) USEPUINDXD, by(daten)

. 
. save usepu_clean.dta, replace
(file usepu_clean.dta not found)
file usepu_clean.dta saved


. *******Coverting Personal COnsumption Expenditure into Monthly Frequency******

. use pce_clean.dta, clear


. gen mdate = ym(year(daten), 3 * quarter(daten) - 2)

. 
. format mdate %tm

. keep mdate PCECC96

. 
. save pce_monthly.dta, replace
(file pce_monthly.dta not found)
file pce_monthly.dta saved

. save "C:\Users\fatem\OneDrive\Desktop\Fall 2025\Consumption\Term Paper\Rough.dta"
file C:\Users\fatem\OneDrive\Desktop\Fall 2025\Consumption\Term Paper\Rough.dta saved

. use pce_monthly.dta, clear

. rename mdate daten

. save pce_monthly.dta, replace
file pce_monthly.dta saved

. *** Now Merging all the Dataset***

. use pce_monthly.dta, clear

. merge 1:1 daten using psavert_clean.dta

    Result                      Number of obs
    -----------------------------------------
    Not matched                         1,095
        from master                       304  (_merge==1)
        from using                        791  (_merge==2)

    Matched                                10  (_merge==3)
    -----------------------------------------

. drop _merge

. 
. merge 1:1 daten using uncertainty_clean.dta

    Result                      Number of obs
    -----------------------------------------
    Not matched                           319
        from master                       312  (_merge==1)
        from using                          7  (_merge==2)

    Matched                               484  (_merge==3)
    -----------------------------------------

. keep if _merge == 3
(319 observations deleted)

. 
. drop _merge

. 
. save merged_final.dta, replace
file merged_final.dta saved

. describe

Contains data from merged_final.dta
 Observations:           484                  
    Variables:             6                  8 Dec 2025 20:52
-------------------------------------------------------------------------------------------------------------------------------
Variable      Storage   Display    Value
    name         type    format    label      Variable label
-------------------------------------------------------------------------------------------------------------------------------
datestr         str10   %-10s                 observation date
daten           int     %tm                   numeric (daily) date
PSAVERT         float   %9.0g                 Personal Saving Rate
daten_num       float   %9.0g                 
PCECC96_mon     double  %10.0g                Interpolation of PCECC96_interp on daten
USEPUINDXM      float   %9.0g                 Economic Policy Uncertainty Index for United States
-------------------------------------------------------------------------------------------------------------------------------
Sorted by: 
     Note: Dataset has changed since last saved.

. ****Stationarity Check + Transformation****

. gen d_pcecc96 = PCECC96_mon - PCECC96_mon[_n-1]
(1 missing value generated)

. 
. gen d_psavert = PSAVERT - PSAVERT[_n-1]
(1 missing value generated)

. format daten %tm

. 
. tsset daten

Time variable: daten, 2721m1 to 3946m1, but with gaps
        Delta: 1 month

. *** Running ADF Test ***
. summarize d_pcecc96 d_psavert

    Variable |        Obs        Mean    Std. dev.       Min        Max
-------------+---------------------------------------------------------
   d_pcecc96 |        483    23.03215    48.71605    -409.38   384.3955
   d_psavert |        483   -.0093168    1.620402      -13.8       19.4

. count if missing(d_pcecc96)
  1

. 
. count if missing(d_psavert)
  1

. list daten d_pcecc96 d_psavert if missing(d_pcecc96) | missing(d_psavert)

     +------------------------------+
     |  daten   d_pce~96   d_psav~t |
     |------------------------------|
  1. | 2721m1          .          . |
     +------------------------------+

. drop if missing(d_pcecc96) | missing(d_psavert)
(1 observation deleted)
. browse

. count
  483

. describe d_pcecc96 d_psavert

Variable      Storage   Display    Value
    name         type    format    label      Variable label
-------------------------------------------------------------------------------------------------------------------------------
d_pcecc96       float   %9.0g                 
d_psavert       float   %9.0g

. list daten d_pcecc96 d_psavert if missing(d_pcecc96) | missing(d_psavert)

. tsset

Time variable: daten, 2723m8 to 3946m1, but with gaps
        Delta: 1 month

. tsset daten, monthly

Time variable: daten, 2723m8 to 3946m1, but with gaps
        Delta: 1 month

. gen monthly_date = mofd(daten)

. 
. format monthly_date %tm

. 
. tsset monthly_date

Time variable: monthly_date, 1985m2 to 2025m4
        Delta: 1 month

. dfuller d_pcecc96, lags(12)

Augmented Dickey–Fuller test for unit root

Variable: d_pcecc96                       Number of obs  = 470
                                          Number of lags =  12

H0: Random walk without drift, d = 0

                                       Dickey–Fuller
                   Test      -------- critical value ---------
              statistic           1%           5%          10%
--------------------------------------------------------------
 Z(t)            -5.948       -3.442       -2.871       -2.570
--------------------------------------------------------------
MacKinnon approximate p-value for Z(t) = 0.0000.

. 
. dfuller d_psavert, lags(12)

Augmented Dickey–Fuller test for unit root

Variable: d_psavert                       Number of obs  = 470
                                          Number of lags =  12

H0: Random walk without drift, d = 0

                                       Dickey–Fuller
                   Test      -------- critical value ---------
              statistic           1%           5%          10%
--------------------------------------------------------------
 Z(t)            -6.518       -3.442       -2.871       -2.570
--------------------------------------------------------------
MacKinnon approximate p-value for Z(t) = 0.0000.

. *** Plotting the Differenced Series ***

. tsline d_pcecc96 d_psavert

. graph export "C:\Users\fatem\OneDrive\Desktop\Fall 2025\Consumption\Term Paper\Plot for differenced series.jpg", as(jpg) name
> ("Graph") quality(100)
file C:\Users\fatem\OneDrive\Desktop\Fall 2025\Consumption\Term Paper\Plot for differenced series.jpg written in JPEG format

. *** VAR Analysis ***

. varsoc d_pcecc96 d_psavert

Lag-order selection criteria

   Sample: 1985m6 thru 2025m4                              Number of obs = 479
  +---------------------------------------------------------------------------+
  | Lag |    LL      LR      df    p     FPE       AIC      HQIC      SBIC    |
  |-----+---------------------------------------------------------------------|
  |   0 | -3399.15                     5040.54    14.201   14.2079   14.2184  |
  |   1 | -3200.99  396.32    4  0.000 2240.79   13.3903   13.4109   13.4426  |
  |   2 |  -3165.4  71.179    4  0.000  1963.9   13.2584   13.2927   13.3455  |
  |   3 | -3133.73  63.328    4  0.000 1749.67   13.1429   13.1909   13.2649  |
  |   4 | -3090.37   86.72*   4  0.000 1484.51*  12.9786*  13.0402*  13.1354* |
  +---------------------------------------------------------------------------+
   * optimal lag
   Endogenous: d_pcecc96 d_psavert
    Exogenous: _cons

. var d_pcecc96 d_psavert, lags(1/4)

Vector autoregression

Sample: 1985m6 thru 2025m4                      Number of obs     =        479
Log likelihood =  -3090.372                     AIC               =   12.97859
FPE            =   1484.515                     HQIC              =   13.04021
Det(Sigma_ml)  =   1377.021                     SBIC              =   13.13535

Equation           Parms      RMSE     R-sq      chi2     P>chi2
----------------------------------------------------------------
d_pcecc96             9     28.7751   0.6597   928.6654   0.0000
d_psavert             9     1.32811   0.3425   249.5393   0.0000
----------------------------------------------------------------

------------------------------------------------------------------------------
             | Coefficient  Std. err.      z    P>|z|     [95% conf. interval]
-------------+----------------------------------------------------------------
d_pcecc96    |
   d_pcecc96 |
         L1. |   .8589068   .0425184    20.20   0.000     .7755723    .9422413
         L2. |   .1524603   .0510492     2.99   0.003     .0524057    .2525149
         L3. |  -.5081946   .0519127    -9.79   0.000    -.6099416   -.4064477
         L4. |   .3961771   .0443254     8.94   0.000     .3093009    .4830532
             |
   d_psavert |
         L1. |    13.1351   .9960486    13.19   0.000     11.18288    15.08732
         L2. |   8.170301   1.247558     6.55   0.000     5.725131    10.61547
         L3. |   2.150238   1.291209     1.67   0.096    -.3804851    4.680961
         L4. |   2.367334   1.127691     2.10   0.036     .1571012    4.577567
             |
       _cons |   2.544517    1.72948     1.47   0.141    -.8452011    5.934234
-------------+----------------------------------------------------------------
d_psavert    |
   d_pcecc96 |
         L1. |  -.0133559   .0019624    -6.81   0.000    -.0172022   -.0095096
         L2. |  -.0084106   .0023562    -3.57   0.000    -.0130286   -.0037926
         L3. |   .0045572    .002396     1.90   0.057    -.0001389    .0092533
         L4. |  -.0057107   .0020458    -2.79   0.005    -.0097204   -.0017009
             |
   d_psavert |
         L1. |  -.6117218   .0459725   -13.31   0.000    -.7018262   -.5216174
         L2. |  -.2631223   .0575809    -4.57   0.000    -.3759787   -.1502658
         L3. |  -.1098188   .0595956    -1.84   0.065     -.226624    .0069864
         L4. |  -.0979306   .0520484    -1.88   0.060    -.1999436    .0040824
             |
       _cons |   .5070781   .0798239     6.35   0.000     .3506262    .6635301
------------------------------------------------------------------------------

. varstable

   Eigenvalue stability condition
  +----------------------------------------+
  |        Eigenvalue        |   Modulus   |
  |--------------------------+-------------|
  |  -.8264396               |    .82644   |
  |   .4333398 +  .6526279i  |   .783394   |
  |   .4333398 -  .6526279i  |   .783394   |
  |   .6768016               |   .676802   |
  |  -.4440253 +  .3089674i  |   .540943   |
  |  -.4440253 -  .3089674i  |   .540943   |
  |    .209097 +  .4560064i  |   .501661   |
  |    .209097 -  .4560064i  |   .501661   |
  +----------------------------------------+
   All the eigenvalues lie inside the unit circle.
   VAR satisfies stability condition.


. *** Set up for IRF analysis ***

. irf set irf_results, replace
(file irf_results.irf created)
(file irf_results.irf now active)

. irf create model_irf, step(12)
(file irf_results.irf updated)

. irf graph oirf

. graph export "C:\Users\fatem\OneDrive\Desktop\Fall 2025\Consumption\Term Paper\IRF graph.jpg", as(jpg) name("Graph") quality(
> 100)
file C:\Users\fatem\OneDrive\Desktop\Fall 2025\Consumption\Term Paper\IRF graph.jpg written in JPEG format

. irf table fevd

Results from model_irf

-------------------------------------------------------------------------------------------------------------------
         |      (1)         (1)         (1)         (2)         (2)         (2)         (3)         (3)         (3)  
    Step |     fevd       Lower       Upper        fevd       Lower       Upper        fevd       Lower       Upper  
---------+---------------------------------------------------------------------------------------------------------
       0 |        0           0           0           0           0           0           0           0           0
       1 |        1           1           1     .020707    -.004532     .045947           0           0           0
       2 |  .815932      .76209     .869773      .04353     .019556     .067505     .184068     .130227      .23791
       3 |  .760255     .687288     .833223     .082283     .052187     .112379     .239745     .166777     .312712
       4 |  .723224     .636391     .810057     .080729     .051398      .11006     .276776     .189943     .363609
       5 |  .724178     .637245     .811111     .081333     .051882     .110785     .275822     .188889     .362755
       6 |   .72346     .636814     .810107      .08524     .052117     .118363      .27654     .189893     .363186
       7 |  .722287     .634903     .809671      .08577     .052063     .119478     .277713     .190329     .365097
       8 |  .723064     .635916     .810212     .085766     .052059     .119473     .276936     .189788     .364084
       9 |  .723213      .63608     .810346     .087725     .052414     .123036     .276787     .189654      .36392
      10 |  .721221     .633723      .80872      .08768     .052405     .122956     .278779      .19128     .366277
      11 |  .721209     .633678     .808739     .087766     .052425     .123107     .278791     .191261     .366322
      12 |  .721488     .633983     .808993     .088103     .052339     .123866     .278512     .191007     .366017
-------------------------------------------------------------------------------------------------------------------

---------------------------------------------
         |      (4)         (4)         (4)  
    Step |     fevd       Lower       Upper  
---------+-----------------------------------
       0 |        0           0           0
       1 |  .979293     .954053     1.00453
       2 |   .95647     .932495     .980444
       3 |  .917717     .887621     .947813
       4 |  .919271      .88994     .948602
       5 |  .918667     .889215     .948118
       6 |   .91476     .881637     .947883
       7 |   .91423     .880522     .947937
       8 |  .914234     .880527     .947941
       9 |  .912275     .876964     .947586
      10 |   .91232     .877044     .947595
      11 |  .912234     .876893     .947575
      12 |  .911897     .876134     .947661
---------------------------------------------
95% lower and upper bounds reported.
(1) irfname = model_irf, impulse = d_pcecc96, and response = d_pcecc96.
(2) irfname = model_irf, impulse = d_pcecc96, and response = d_psavert.
(3) irfname = model_irf, impulse = d_psavert, and response = d_pcecc96.
(4) irfname = model_irf, impulse = d_psavert, and response = d_psavert.

. irf graph fevd

. graph export "C:\Users\fatem\OneDrive\Desktop\Fall 2025\Consumption\Term Paper\IRF Graph FEVD.jpg", as(jpg) name("Graph") qua
> lity(100)
file C:\Users\fatem\OneDrive\Desktop\Fall 2025\Consumption\Term Paper\IRF Graph FEVD.jpg written in JPEG format

*** Robustness Check ***
 *** Granger Casuality Test ***

. vargranger

   Granger causality Wald tests
  +------------------------------------------------------------------+
  |          Equation           Excluded |   chi2     df Prob > chi2 |
  |--------------------------------------+---------------------------|
  |         d_pcecc96          d_psavert |   180.2     4    0.000    |
  |         d_pcecc96                ALL |   180.2     4    0.000    |
  |--------------------------------------+---------------------------|
  |         d_psavert          d_pcecc96 |  158.53     4    0.000    |
  |         d_psavert                ALL |  158.53     4    0.000    |
  +------------------------------------------------------------------+

*** Cholesky Ordering ***
. irf create model_cholesky, replace
(irfname model_cholesky not found in irf_results.irf)
(file irf_results.irf updated)

. 
. irf graph oirf

. irf graph oirf, impulse(d_pcecc96) response(d_psavert) irf(model_cholesky)

. graph export "C:\Users\fatem\OneDrive\Desktop\Fall 2025\Consumption\Term Paper\How savings respond to uncertainty shocks.jpg"
> , as(jpg) name("Graph") quality(100)
file C:\Users\fatem\OneDrive\Desktop\Fall 2025\Consumption\Term Paper\How savings respond to uncertainty shocks.jpg written in 
> JPEG format

. irf table oirf, impulse(d_pcecc96) response(d_psavert) irf(model_cholesky)

Results from model_cholesky

---------------------------------------------
         |      (1)         (1)         (1)  
    Step |     oirf       Lower       Upper  
---------+-----------------------------------
       0 | -.189312    -.306515     -.07211
       1 | -.264884    -.395361    -.134407
       2 |  -.32165    -.440087    -.203213
       3 | -.011081    -.131181      .10902
       4 | -.046134    -.152612     .060345
       5 |   .10726     .026404     .188117
       6 | -.039445     -.12026     .041369
       7 | -.000625    -.064661     .063412
       8 | -.075284    -.125847    -.024722
---------------------------------------------
95% lower and upper bounds reported.
(1) irfname = model_cholesky, impulse = d_pcecc96, and response = d_psavert.

. *** Alternative Cholesky Ordering ***

. var d_psavert d_pcecc96, lags(1/4)

Vector autoregression

Sample: 1985m6 thru 2025m4                      Number of obs     =        479
Log likelihood =  -3090.372                     AIC               =   12.97859
FPE            =   1484.515                     HQIC              =   13.04021
Det(Sigma_ml)  =   1377.021                     SBIC              =   13.13535

Equation           Parms      RMSE     R-sq      chi2     P>chi2
----------------------------------------------------------------
d_psavert             9     1.32811   0.3425   249.5393   0.0000
d_pcecc96             9     28.7751   0.6597   928.6654   0.0000
----------------------------------------------------------------

------------------------------------------------------------------------------
             | Coefficient  Std. err.      z    P>|z|     [95% conf. interval]
-------------+----------------------------------------------------------------
d_psavert    |
   d_psavert |
         L1. |  -.6117218   .0459725   -13.31   0.000    -.7018262   -.5216174
         L2. |  -.2631223   .0575809    -4.57   0.000    -.3759787   -.1502658
         L3. |  -.1098188   .0595956    -1.84   0.065     -.226624    .0069864
         L4. |  -.0979306   .0520484    -1.88   0.060    -.1999436    .0040824
             |
   d_pcecc96 |
         L1. |  -.0133559   .0019624    -6.81   0.000    -.0172022   -.0095096
         L2. |  -.0084106   .0023562    -3.57   0.000    -.0130286   -.0037926
         L3. |   .0045572    .002396     1.90   0.057    -.0001389    .0092533
         L4. |  -.0057107   .0020458    -2.79   0.005    -.0097204   -.0017009
             |
       _cons |   .5070781   .0798239     6.35   0.000     .3506262    .6635301
-------------+----------------------------------------------------------------
d_pcecc96    |
   d_psavert |
         L1. |    13.1351   .9960486    13.19   0.000     11.18288    15.08732
         L2. |   8.170301   1.247558     6.55   0.000     5.725131    10.61547
         L3. |   2.150238   1.291209     1.67   0.096    -.3804851    4.680961
         L4. |   2.367334   1.127691     2.10   0.036     .1571012    4.577567
             |
   d_pcecc96 |
         L1. |   .8589068   .0425184    20.20   0.000     .7755723    .9422413
         L2. |   .1524603   .0510492     2.99   0.003     .0524057    .2525149
         L3. |  -.5081946   .0519127    -9.79   0.000    -.6099416   -.4064477
         L4. |   .3961771   .0443254     8.94   0.000     .3093009    .4830532
             |
       _cons |   2.544517    1.72948     1.47   0.141    -.8452011    5.934234
------------------------------------------------------------------------------
. irf create cholesky_rev, replace
(irfname cholesky_rev not found in irf_results.irf)
(file irf_results.irf updated)

. irf graph oirf, impulse(d_pcecc96) response(d_psavert) irf(cholesky_rev)

. graph save "Graph" "C:\Users\fatem\OneDrive\Desktop\Fall 2025\Consumption\Term Paper\Cholesky Rev.gph"
file C:\Users\fatem\OneDrive\Desktop\Fall 2025\Consumption\Term Paper\Cholesky Rev.gph saved

. *** Different Lag Lengths ***

. * Lag 2

. 
. var d_pcecc96 d_psavert, lags(1/2)

Vector autoregression

Sample: 1985m4 thru 2025m4                      Number of obs     =        481
Log likelihood =  -3177.733                     AIC               =   13.25461
FPE            =   1956.387                     HQIC              =   13.28873
Det(Sigma_ml)  =   1876.706                     SBIC              =   13.34142

Equation           Parms      RMSE     R-sq      chi2     P>chi2
----------------------------------------------------------------
d_pcecc96             5     33.2891   0.5388   562.0076   0.0000
d_psavert             5     1.34231   0.3214   227.7943   0.0000
----------------------------------------------------------------

------------------------------------------------------------------------------
             | Coefficient  Std. err.      z    P>|z|     [95% conf. interval]
-------------+----------------------------------------------------------------
d_pcecc96    |
   d_pcecc96 |
         L1. |   .8444329   .0440112    19.19   0.000     .7581726    .9306933
         L2. |   -.059542   .0494602    -1.20   0.229    -.1564822    .0373981
             |
   d_psavert |
         L1. |   13.79329    1.13767    12.12   0.000      11.5635    16.02308
         L2. |   6.714427   1.252655     5.36   0.000     4.259268    9.169586
             |
       _cons |   5.169049   1.784733     2.90   0.004     1.671037    8.667062
-------------+----------------------------------------------------------------
d_psavert    |
   d_pcecc96 |
         L1. |  -.0138498   .0017747    -7.80   0.000     -.017328   -.0103715
         L2. |  -.0063829   .0019944    -3.20   0.001    -.0102918    -.002474
             |
   d_psavert |
         L1. |  -.6091432   .0458742   -13.28   0.000    -.6990549   -.5192315
         L2. |  -.2104155   .0505107    -4.17   0.000    -.3094147   -.1114164
             |
       _cons |   .4525673   .0719656     6.29   0.000     .3115173    .5936173
------------------------------------------------------------------------------

. 
. irf create model_lag2, replace
(irfname model_lag2 not found in irf_results.irf)
(file irf_results.irf updated)

. 
. irf graph oirf, impulse(d_pcecc96) response(d_psavert) irf(model_lag2)

. 
. 
. 
. * Lag 6

. 
. var d_pcecc96 d_psavert, lags(1/6)

Vector autoregression

Sample: 1985m8 thru 2025m4                      Number of obs     =        477
Log likelihood =  -3029.927                     AIC               =   12.81311
FPE            =   1258.132                     HQIC              =   12.90243
Det(Sigma_ml)  =   1128.158                     SBIC              =   13.04027

Equation           Parms      RMSE     R-sq      chi2     P>chi2
----------------------------------------------------------------
d_pcecc96            13     28.1997   0.6773   1001.132   0.0000
d_psavert            13     1.23178   0.4407   375.8429   0.0000
----------------------------------------------------------------

------------------------------------------------------------------------------
             | Coefficient  Std. err.      z    P>|z|     [95% conf. interval]
-------------+----------------------------------------------------------------
d_pcecc96    |
   d_pcecc96 |
         L1. |   .8812219   .0455226    19.36   0.000     .7919993    .9704445
         L2. |   .1574601   .0609284     2.58   0.010     .0380426    .2768777
         L3. |  -.6221538   .0579923   -10.73   0.000    -.7358166    -.508491
         L4. |   .4669292     .05586     8.36   0.000     .3574457    .5764127
         L5. |   .0454926    .059757     0.76   0.446     -.071629    .1626142
         L6. |   -.120432   .0476969    -2.52   0.012    -.2139163   -.0269477
             |
   d_psavert |
         L1. |   12.94221   .9930591    13.03   0.000     10.99585    14.88857
         L2. |   7.855203   1.257139     6.25   0.000     5.391256    10.31915
         L3. |   1.635752   1.315509     1.24   0.214    -.9425994    4.214103
         L4. |   3.271349   1.298118     2.52   0.012     .7270847    5.815614
         L5. |   2.025652   1.284645     1.58   0.115    -.4922055    4.543509
         L6. |   1.524036   1.121679     1.36   0.174    -.6744137    3.722487
             |
       _cons |   4.699791   1.879165     2.50   0.012     1.016696    8.382887
-------------+----------------------------------------------------------------
d_psavert    |
   d_pcecc96 |
         L1. |  -.0091504   .0019885    -4.60   0.000    -.0130477   -.0052531
         L2. |  -.0184171   .0026614    -6.92   0.000    -.0236333   -.0132009
         L3. |   .0148182   .0025331     5.85   0.000     .0098533     .019783
         L4. |  -.0010974     .00244    -0.45   0.653    -.0058797    .0036849
         L5. |   -.020423   .0026102    -7.82   0.000     -.025539   -.0153071
         L6. |   .0162265   .0020834     7.79   0.000     .0121431      .02031
             |
   d_psavert |
         L1. |  -.5760964   .0433774   -13.28   0.000    -.6611146   -.4910782
         L2. |  -.2832166   .0549126    -5.16   0.000    -.3908433   -.1755898
         L3. |  -.0692102   .0574623    -1.20   0.228    -.1818342    .0434138
         L4. |  -.1421091   .0567026    -2.51   0.012    -.2532441    -.030974
         L5. |  -.1701518   .0561141    -3.03   0.002    -.2801334   -.0601702
         L6. |   .1059063   .0489956     2.16   0.031     .0098766     .201936
             |
       _cons |   .3967772   .0820831     4.83   0.000     .2358973    .5576571
------------------------------------------------------------------------------

. 
. irf create model_lag6, replace
(irfname model_lag6 not found in irf_results.irf)
(file irf_results.irf updated)

. 
. irf graph oirf, impulse(d_pcecc96) response(d_psavert) irf(model_lag6)

. graph export "C:\Users\fatem\OneDrive\Desktop\Fall 2025\Consumption\Term Paper\Lag 6.jpg", as(jpg) name("Graph") quality(100)
file C:\Users\fatem\OneDrive\Desktop\Fall 2025\Consumption\Term Paper\Lag 6.jpg written in JPEG format

. irf create model_lag2, replace
(file irf_results.irf updated)

. irf graph oirf, impulse(d_pcecc96) response(d_psavert) irf(model_lag2)

. graph export "C:\Users\fatem\OneDrive\Desktop\Fall 2025\Consumption\Term Paper\Lag 2.jpg", as(jpg) name("Graph") quality(100)
file C:\Users\fatem\OneDrive\Desktop\Fall 2025\Consumption\Term Paper\Lag 2.jpg written in JPEG format
.
. *** Using sVAR with short-run restrictions ***
. matrix A = (1, 0 \ ., 1)

. svar d_pcecc96 d_psavert, aeq(A) lags(1/4)
Estimating short-run parameters

Iteration 0:  Log likelihood = -395682.98  
Iteration 1:  Log likelihood = -264674.86  (backed up)
Iteration 2:  Log likelihood = -219561.83  (backed up)
Iteration 3:  Log likelihood = -204027.03  (backed up)
Iteration 4:  Log likelihood = -198677.58  (backed up)
Iteration 5:  Log likelihood = -196835.48  (backed up)
Iteration 6:  Log likelihood = -196201.14  (backed up)
Iteration 7:  Log likelihood = -195982.71  (backed up)
Iteration 8:  Log likelihood = -195907.49  (backed up)
Iteration 9:  Log likelihood = -195881.59  (backed up)
Iteration 10: Log likelihood = -195872.67  (backed up)
Iteration 11: Log likelihood =  -195869.6  (backed up)
Iteration 12: Log likelihood = -195868.54  (backed up)
Iteration 13: Log likelihood = -195868.18  (backed up)
Iteration 14: Log likelihood = -195868.05  (backed up)
Iteration 15: Log likelihood = -195868.01  (backed up)
Iteration 16: Log likelihood = -195867.99  (backed up)
Iteration 17: Log likelihood = -195867.99  (backed up)
Iteration 18: Log likelihood = -195867.99  (backed up)
Iteration 19: Log likelihood = -195867.99  (backed up)
Iteration 20: Log likelihood = -195867.99  (backed up)
Iteration 21: Log likelihood = -195867.99  (backed up)
Iteration 22: Log likelihood = -195867.99  (backed up)
Iteration 23: Log likelihood = -195867.99  (backed up)
Iteration 24: Log likelihood = -195867.99  (backed up)
Iteration 25: Log likelihood = -195867.99  (backed up)
Iteration 26: Log likelihood = -195867.99  (backed up)
Iteration 27: Log likelihood = -195867.99  (backed up)
Iteration 28: Log likelihood = -195867.99  (backed up)
Iteration 29: Log likelihood = -195867.99  (backed up)
Iteration 30: Log likelihood = -195867.99  (backed up)

Structural vector autoregression

 ( 1)  [/A]1_1 = 1
 ( 2)  [/A]1_2 = 0
 ( 3)  [/A]2_2 = 1
 ( 4)  [/B]1_1 = 1
 ( 5)  [/B]1_2 = 0
 ( 6)  [/B]2_1 = 0
 ( 7)  [/B]2_2 = 1

Sample: 1985m6 thru 2025m4                      Number of obs     =        479
Overidentified model                            Log likelihood    =    -195868

------------------------------------------------------------------------------
             | Coefficient  Std. err.      z    P>|z|     [95% conf. interval]
-------------+----------------------------------------------------------------
/A           |
         1_1 |          1  (constrained)
         2_1 |   .0066418   .0456912     0.15   0.884    -.0829112    .0961949
         1_2 |          0  (constrained)
         2_2 |          1  (constrained)
-------------+----------------------------------------------------------------
/B           |
         1_1 |          1  (constrained)
         2_1 |          0  (constrained)
         1_2 |          0  (constrained)
         2_2 |          1  (constrained)
------------------------------------------------------------------------------
LR test of identifying restrictions: chi2(2) = 3.9e+05     Prob > chi2 = 0.000

. irf create model_svar, replace
(irfname model_svar not found in irf_results.irf)
(file irf_results.irf updated)

. irf graph oirf, impulse(d_pcecc96) response(d_psavert) irf(model_svar)

. graph save "Graph" "C:\Users\fatem\OneDrive\Desktop\Fall 2025\Consumption\Term Paper\sVAR.gph"
file C:\Users\fatem\OneDrive\Desktop\Fall 2025\Consumption\Term Paper\sVAR.gph saved
.
. *** Forecast Error Variance Decomposition Robustness ***

. irf create model_fevdcheck, replace
(irfname model_fevdcheck not found in irf_results.irf)
(file irf_results.irf updated)

. 
. irf table fevd, irf(model_fevdcheck)

Results from model_fevdcheck

-------------------------------------------------------------------------------------------------------------------
         |      (1)         (1)         (1)         (2)         (2)         (2)         (3)         (3)         (3)  
    Step |     fevd       Lower       Upper        fevd       Lower       Upper        fevd       Lower       Upper  
---------+---------------------------------------------------------------------------------------------------------
       0 |        0           0           0           0           0           0           0           0           0
       1 |        1           1           1     .020707    -.004532     .045947           0           0           0
       2 |  .815932      .76209     .869773      .04353     .019556     .067505     .184068     .130227      .23791
       3 |  .760255     .687288     .833223     .082283     .052187     .112379     .239745     .166777     .312712
       4 |  .723224     .636391     .810057     .080729     .051398      .11006     .276776     .189943     .363609
       5 |  .724178     .637245     .811111     .081333     .051882     .110785     .275822     .188889     .362755
       6 |   .72346     .636814     .810107      .08524     .052117     .118363      .27654     .189893     .363186
       7 |  .722287     .634903     .809671      .08577     .052063     .119478     .277713     .190329     .365097
       8 |  .723064     .635916     .810212     .085766     .052059     .119473     .276936     .189788     .364084
-------------------------------------------------------------------------------------------------------------------

---------------------------------------------
         |      (4)         (4)         (4)  
    Step |     fevd       Lower       Upper  
---------+-----------------------------------
       0 |        0           0           0
       1 |  .979293     .954053     1.00453
       2 |   .95647     .932495     .980444
       3 |  .917717     .887621     .947813
       4 |  .919271      .88994     .948602
       5 |  .918667     .889215     .948118
       6 |   .91476     .881637     .947883
       7 |   .91423     .880522     .947937
       8 |  .914234     .880527     .947941
---------------------------------------------
95% lower and upper bounds reported.
(1) irfname = model_fevdcheck, impulse = d_pcecc96, and response = d_pcecc96.
(2) irfname = model_fevdcheck, impulse = d_pcecc96, and response = d_psavert.
(3) irfname = model_fevdcheck, impulse = d_psavert, and response = d_pcecc96.
(4) irfname = model_fevdcheck, impulse = d_psavert, and response = d_psavert.

. 


