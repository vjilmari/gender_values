---
title: "Examinations of VBMT from equal number of country-time samples"
output: 
  html_document: 
    toc: true
    keep_md: true
---



# Packages

Two packages are not in CRAN and need different installation
* devtools::install_github("vjilmari/vjihelpers")
* install.packages("ggflags", repos = c("https://jimjam-slam.r-universe.dev","https://cloud.r-project.org")) 


``` r
library(multid)
library(rio)
library(dplyr)
library(lme4)
library(emmeans)
library(vjihelpers)
library(MetBrewer)
library(ggplot2)
library(finalfit)
library(ggflags) 
library(lmerTest)
library(metafor)
library(ggpubr)
library(Hmisc)
library(r2mlm)
library(tidyr)
library(stringr)
library(apaTables)
library(tibble)
```

# Custom functions


``` r
# toster function for equivalence testing

tost_z <- function(est, se, low, high, alpha = 0.10) {
  # z statistics
  z_low  <- (est - low)  / se     # H0: theta <= low  vs H1: theta > low
  z_high <- (est - high) / se     # H0: theta >= high vs H1: theta < high
  
  # one-sided p-values
  p_low  <- 1 - pnorm(z_low)
  p_high <- pnorm(z_high)
  
  # CI corresponding to TOST (1 - 2*alpha)
  z_crit <- qnorm(1 - alpha)
  ci_low  <- est - z_crit * se
  ci_high <- est + z_crit * se
  
  # equivalence decision
  equivalent <- (p_low < alpha) && (p_high < alpha)
  
  list(
    estimate     = est,
    se           = se,
    low_bound    = low,
    high_bound   = high,
    alpha        = alpha,
    z_low        = z_low,
    p_low        = p_low,
    z_high       = z_high,
    p_high       = p_high,
    ci_level     = 1 - 2*alpha,
    ci_lower     = ci_low,
    ci_upper     = ci_high,
    equivalent   = equivalent
  )
}


# TOST function for equivalence testing with t distribution

tost_t <- function(est, se, low, high, df, alpha = 0.10) {
  # t statistics
  t_low  <- (est - low)  / se   # H0: theta <= low  vs H1: theta > low
  t_high <- (est - high) / se   # H0: theta >= high vs H1: theta < high
  
  # one-sided p-values
  p_low  <- 1 - pt(t_low,  df = df)
  p_high <-     pt(t_high, df = df)
  
  # CI corresponding to TOST (1 - 2*alpha)
  t_crit <- qt(1 - alpha, df = df)
  ci_low  <- est - t_crit * se
  ci_high <- est + t_crit * se
  
  # equivalence decision
  equivalent <- (p_low < alpha) && (p_high < alpha)
  
  list(
    estimate   = est,
    se         = se,
    df         = df,
    low_bound  = low,
    high_bound = high,
    alpha      = alpha,
    t_low      = t_low,
    p_low      = p_low,
    t_high     = t_high,
    p_high     = p_high,
    ci_level   = 1 - 2 * alpha,
    ci_lower   = ci_low,
    ci_upper   = ci_high,
    equivalent = equivalent
  )
}
```


# Data

* Loads the preprocessed European Social Survey data file made within "Preparations_ESS" code


``` r
load("../data/fdat.rdata")
```

* Include ISO data for country-names


``` r
# read country names
ISO<-read.csv2("../data/ISO.csv")
```


# Exclusions

* exclude participants with missing values or gender


``` r
value.vars<-
  c("con","tra",
  "ben","uni",
  "sdi","sti",
  "hed","ach",
  "pow","sec")

fdat$miss_values<-
  rowSums(is.na(fdat[,value.vars]))
table(fdat$miss_values)
```

```
## 
##      0      1      2      3      4      5      6      7      8      9     10 
## 500044   2791    835    442    267    204    142    143    156    177  35470
```

``` r
fdat<-fdat %>%
  filter(miss_values==0 & !is.na(gndr.bin))
```

# Construct value-based gender-typicality measure

* Use 100 men and women from each country x time data fold for training set


``` r
# set seed number for reproducibility in training/testing split
set.seed(13032023)
# check the cross-validation fold sizes
table(fdat$cntry_time)
```

```
## 
##  AL_6  AT_1 AT_11  AT_2  AT_3  AT_7  AT_8  AT_9  BE_1 BE_10 BE_11  BE_2  BE_3  BE_4  BE_5  BE_6  BE_7 
##  1117  2254  2314  2198  2348  1795  1993  2477  1830  1334  1577  1771  1796  1754  1699  1862  1767 
##  BE_8  BE_9 BG_10 BG_11  BG_3  BG_4  BG_5  BG_6  BG_9  CH_1 CH_10 CH_11  CH_2  CH_3  CH_4  CH_5  CH_6 
##  1759  1756  2697  2218  1295  2144  2371  2179  1926  2024  1505  1368  2110  1780  1753  1491  1483 
##  CH_7  CH_8  CH_9 CY_11  CY_3  CY_4  CY_5  CY_6  CY_9  CZ_1 CZ_10  CZ_2  CZ_4  CZ_5  CZ_6  CZ_7  CZ_8 
##  1521  1504  1517   667   978  1210  1053  1110   773  1208  2369  2557  1986  2335  1973  1862  2252 
##  CZ_9  DE_1 DE_11  DE_2  DE_3  DE_4  DE_5  DE_6  DE_7  DE_8  DE_9  DK_1  DK_2  DK_3  DK_4  DK_5  DK_6 
##  2343  2819  2381  2840  2884  2732  3007  2935  3006  2821  2328  1470  1458  1461  1581  1564  1621 
##  DK_7  DK_9 EE_10 EE_11  EE_2  EE_3  EE_4  EE_5  EE_6  EE_7  EE_8  EE_9  ES_1 ES_11  ES_2  ES_3  ES_4 
##  1483  1554  1538  1282  1948  1466  1646  1793  2345  2036  2007  1899  1712  1833  1623  1847  2562 
##  ES_5  ES_6  ES_7  ES_8  ES_9  FI_1 FI_10 FI_11  FI_2  FI_3  FI_4  FI_5  FI_6  FI_7  FI_8  FI_9  FR_1 
##  1881  1871  1907  1929  1619  1763  1561  1524  1701  1649  1901  1649  2158  2050  1903  1735  1355 
## FR_10 FR_11  FR_2  FR_3  FR_4  FR_5  FR_6  FR_7  FR_8  FR_9  GB_1 GB_10 GB_11  GB_2  GB_3  GB_4  GB_5 
##  1951  1745  1699  1983  2067  1723  1960  1902  2057  1982  1798  1131  1529  1864  2353  2311  2374 
##  GB_6  GB_7  GB_8  GB_9  GR_1 GR_10 GR_11  GR_2  GR_4  GR_5 HR_10 HR_11  HR_4  HR_5  HR_9  HU_1 HU_10 
##  2261  2231  1942  2183  2551  2768  2745  2399  2063  2669  1564  1548  1430  1601  1781  1634  1816 
## HU_11  HU_2  HU_3  HU_4  HU_5  HU_6  HU_7  HU_8  HU_9  IE_1 IE_10 IE_11  IE_2  IE_3  IE_4  IE_5  IE_6 
##  2117  1460  1462  1430  1473  1968  1520  1458  1643  1916  1751  1985  1187  1589  1757  2400  2616 
##  IE_7  IE_8  IE_9  IL_1 IL_11  IL_4  IL_5  IL_6  IL_7  IL_8 IS_10 IS_11  IS_2  IS_6  IS_8  IS_9 IT_10 
##  2380  2746  2189  2279   893  2382  2212  2378  2351  2366   886   825   524   739   841   844  2573 
## IT_11  IT_6  IT_8  IT_9 LT_10 LT_11  LT_5  LT_6  LT_7  LT_8  LT_9  LU_2 LV_11  LV_4  LV_9 ME_10 ME_11 
##  2783   909  2531  2660  1606  1337  1632  2108  2241  2079  1677  1614  1225  1970   891  1248  1581 
##  ME_9 MK_10  NL_1 NL_10 NL_11  NL_2  NL_3  NL_4  NL_5  NL_6  NL_7  NL_8  NL_9  NO_1 NO_10 NO_11  NO_2 
##  1188  1400  2337  1466  1677  1858  1860  1724  1801  1828  1823  1669  1657  1819  1408  1318  1575 
##  NO_3  NO_4  NO_5  NO_6  NO_7  NO_8  NO_9  PL_1 PL_11  PL_2  PL_3  PL_4  PL_5  PL_6  PL_7  PL_8  PL_9 
##  1550  1391  1530  1610  1423  1530  1396  2065  1423  1683  1685  1596  1719  1866  1594  1675  1443 
##  PT_1 PT_10 PT_11  PT_2  PT_3  PT_4  PT_5  PT_6  PT_7  PT_8  PT_9  RO_4 RS_11  RS_9  RU_3  RU_4  RU_5 
##  1482  1827  1366  2024  2182  2337  2139  2138  1242  1254  1045  2104  1511  1969  2339  2446  2557 
##  RU_6  RU_8  SE_1 SE_11  SE_2  SE_3  SE_4  SE_5  SE_6  SE_7  SE_8  SE_9  SI_1 SI_10 SI_11  SI_2  SI_3 
##  2429  2374  1682  1204  1678  1604  1556  1463  1838  1761  1526  1510  1488  1232  1227  1384  1465 
##  SI_4  SI_5  SI_6  SI_7  SI_8  SI_9 SK_10 SK_11  SK_2  SK_3  SK_4  SK_5  SK_6  SK_9  TR_2  TR_4 UA_11 
##  1257  1369  1244  1189  1295  1307  1395  1414  1425  1711  1789  1803  1827  1061  1790  2305  2589 
##  UA_2  UA_3  UA_4  UA_5  UA_6  XK_6 
##  1896  1885  1766  1779  2064  1244
```

``` r
# run the analysis
value_typ<-
  D_regularized(data=fdat,mv.vars = value.vars,
                group.var = "gndr.bin",group.values = c(1,0),
                out = T,fold = T,fold.var = "cntry_time",size = 100,
                pcc = T,auc=T,pred.prob = T,append.data=T)
```
## Gender typicality results and description


### Training sample size


``` r
# n folds in the training phase
length(unique(fdat$cntry_time))
```

```
## [1] 278
```

``` r
length(unique(fdat$cntry_time))*200
```

```
## [1] 55600
```

### Coefficients for value variables in training set


``` r
round(coefficients(value_typ$cv.mod,s = "lambda.min"),2)
```

```
## 11 x 1 sparse Matrix of class "dgCMatrix"
##             lambda.min
## (Intercept)       0.70
## con               0.09
## tra              -0.03
## ben              -0.27
## uni              -0.08
## sdi               0.10
## sti               0.06
## hed               0.06
## ach               0.07
## pow               0.14
## sec              -0.17
```

``` r
plot(value_typ$cv.mod)
```

![](VBMT_balanced_sampling_of_country_time_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

### Description of gender differences/prediction in testing dataset


``` r
# save to summary tab
sum_tab<-value_typ$D
# print the country X time -fold results
round(sum_tab,2)
```

```
##        n.1  n.0   m.1   m.0 sd.1 sd.0 pooled.sd diff    D pcc.1 pcc.0 pcc.total  auc pooled.sd.1
## AL_6   419  498  0.05 -0.01 0.35 0.36      0.35 0.06 0.17  0.54  0.51      0.52 0.54        0.38
## AT_1   941 1113  0.16 -0.03 0.43 0.45      0.44 0.19 0.44  0.65  0.53      0.59 0.63        0.38
## AT_11  875 1239  0.03 -0.12 0.34 0.34      0.34 0.15 0.46  0.54  0.65      0.60 0.63        0.38
## AT_2   910 1088  0.21 -0.01 0.40 0.43      0.42 0.22 0.52  0.70  0.52      0.60 0.65        0.38
## AT_3   994 1154  0.19 -0.02 0.43 0.40      0.42 0.21 0.51  0.70  0.51      0.60 0.65        0.38
## AT_7   753  842  0.07 -0.10 0.36 0.37      0.36 0.18 0.48  0.58  0.62      0.60 0.63        0.38
## AT_8   792 1001  0.12 -0.04 0.40 0.42      0.41 0.16 0.40  0.62  0.55      0.58 0.62        0.38
## AT_9  1043 1234  0.07 -0.12 0.36 0.39      0.38 0.19 0.49  0.59  0.62      0.60 0.64        0.38
## BE_1   845  785  0.07 -0.06 0.37 0.39      0.38 0.13 0.35  0.61  0.58      0.59 0.60        0.38
## BE_10  570  564  0.05 -0.11 0.32 0.34      0.33 0.17 0.51  0.58  0.63      0.61 0.64        0.38
## BE_11  702  675  0.02 -0.09 0.33 0.36      0.34 0.11 0.33  0.55  0.58      0.56 0.60        0.38
## BE_2   770  801  0.09 -0.09 0.35 0.38      0.36 0.19 0.51  0.62  0.59      0.60 0.64        0.38
## BE_3   740  856  0.07 -0.07 0.36 0.36      0.36 0.14 0.39  0.60  0.56      0.58 0.61        0.38
## BE_4   760  794  0.08 -0.07 0.34 0.35      0.34 0.15 0.43  0.62  0.55      0.58 0.62        0.38
## BE_5   719  780  0.06 -0.08 0.33 0.36      0.35 0.14 0.41  0.59  0.59      0.59 0.62        0.38
## BE_6   811  851  0.07 -0.08 0.32 0.35      0.33 0.15 0.44  0.60  0.56      0.58 0.62        0.38
## BE_7   795  772  0.05 -0.08 0.33 0.33      0.33 0.13 0.40  0.58  0.56      0.57 0.61        0.38
## BE_8   784  775  0.08 -0.08 0.33 0.34      0.34 0.17 0.49  0.60  0.56      0.58 0.63        0.38
## BE_9   765  791  0.08 -0.06 0.32 0.36      0.34 0.14 0.40  0.59  0.55      0.57 0.61        0.38
## BG_10 1171 1326  0.18  0.03 0.38 0.42      0.40 0.16 0.39  0.71  0.45      0.58 0.61        0.38
## BG_11  957 1061  0.25  0.05 0.41 0.45      0.43 0.21 0.48  0.73  0.46      0.59 0.64        0.38
## BG_3   409  686 -0.02 -0.19 0.42 0.45      0.44 0.17 0.39  0.50  0.67      0.61 0.61        0.38
## BG_4   848 1096 -0.06 -0.25 0.43 0.46      0.44 0.19 0.42  0.44  0.70      0.58 0.62        0.38
## BG_5   942 1229  0.06 -0.12 0.37 0.42      0.40 0.17 0.44  0.59  0.59      0.59 0.62        0.38
## BG_6   830 1149  0.05 -0.13 0.38 0.41      0.40 0.18 0.45  0.57  0.61      0.60 0.63        0.38
## BG_9   763  963  0.08 -0.07 0.39 0.44      0.42 0.15 0.37  0.59  0.56      0.58 0.60        0.38
## CH_1   874  950  0.07 -0.09 0.34 0.34      0.34 0.16 0.47  0.59  0.60      0.59 0.63        0.38
## CH_10  672  633  0.06 -0.11 0.33 0.33      0.33 0.17 0.51  0.60  0.61      0.60 0.65        0.38
## CH_11  590  578  0.02 -0.12 0.31 0.35      0.33 0.14 0.42  0.56  0.64      0.60 0.62        0.38
## CH_2   834 1076  0.06 -0.10 0.35 0.35      0.35 0.17 0.48  0.58  0.61      0.60 0.63        0.38
## CH_3   702  878  0.04 -0.12 0.34 0.32      0.33 0.16 0.48  0.54  0.65      0.61 0.63        0.38
## CH_4   697  856  0.11 -0.14 0.35 0.34      0.34 0.25 0.71  0.63  0.67      0.65 0.70        0.38
## CH_5   662  629  0.05 -0.11 0.33 0.33      0.33 0.16 0.49  0.57  0.63      0.60 0.64        0.38
## CH_6   643  640  0.06 -0.10 0.31 0.34      0.32 0.16 0.48  0.59  0.60      0.60 0.64        0.38
## CH_7   659  662  0.08 -0.11 0.31 0.33      0.32 0.19 0.61  0.62  0.64      0.63 0.67        0.38
## CH_8   681  623  0.07 -0.10 0.31 0.32      0.32 0.17 0.53  0.60  0.64      0.62 0.65        0.38
## CH_9   667  650  0.05 -0.11 0.31 0.33      0.32 0.16 0.50  0.56  0.63      0.60 0.64        0.38
## CY_11  203  264 -0.03 -0.26 0.42 0.43      0.42 0.22 0.52  0.45  0.73      0.61 0.65        0.38
## CY_3   368  410 -0.03 -0.12 0.36 0.41      0.39 0.09 0.23  0.51  0.59      0.55 0.57        0.38
## CY_4   512  498 -0.01 -0.13 0.37 0.39      0.38 0.12 0.33  0.51  0.59      0.55 0.58        0.38
## CY_5   375  478 -0.09 -0.21 0.37 0.46      0.42 0.12 0.28  0.43  0.67      0.56 0.58        0.38
## CY_6   391  519 -0.11 -0.27 0.39 0.43      0.41 0.16 0.38  0.40  0.72      0.58 0.60        0.38
## CY_9   263  310 -0.19 -0.27 0.41 0.40      0.41 0.08 0.19  0.31  0.73      0.54 0.55        0.38
## CZ_1   472  536  0.10 -0.18 0.46 0.44      0.45 0.27 0.61  0.59  0.68      0.64 0.67        0.38
## CZ_10  928 1241  0.32  0.15 0.40 0.45      0.43 0.17 0.40  0.79  0.35      0.54 0.60        0.38
## CZ_2  1089 1268  0.24 -0.04 0.42 0.48      0.45 0.29 0.64  0.73  0.52      0.62 0.68        0.38
## CZ_4   863  923  0.25  0.02 0.42 0.47      0.45 0.23 0.52  0.74  0.48      0.61 0.64        0.38
## CZ_5  1069 1066  0.27  0.04 0.41 0.45      0.43 0.23 0.54  0.76  0.44      0.60 0.64        0.38
## CZ_6   875  898  0.28  0.07 0.41 0.48      0.45 0.20 0.46  0.77  0.42      0.59 0.62        0.38
## CZ_7   704  958  0.30  0.08 0.40 0.43      0.42 0.22 0.52  0.79  0.40      0.56 0.64        0.38
## CZ_8   989 1063  0.29  0.16 0.38 0.41      0.39 0.14 0.35  0.79  0.33      0.55 0.59        0.38
## CZ_9   921 1222  0.26  0.07 0.40 0.44      0.42 0.19 0.44  0.75  0.43      0.56 0.62        0.38
## DE_1  1249 1370  0.05 -0.13 0.41 0.42      0.42 0.18 0.44  0.56  0.64      0.60 0.63        0.38
## DE_11 1092 1089 -0.07 -0.25 0.35 0.36      0.36 0.17 0.49  0.42  0.75      0.58 0.64        0.38
## DE_2  1267 1373  0.10 -0.11 0.40 0.39      0.39 0.21 0.53  0.59  0.62      0.61 0.64        0.38
## DE_3  1319 1365  0.10 -0.10 0.38 0.39      0.39 0.19 0.50  0.60  0.60      0.60 0.64        0.38
## DE_4  1342 1190  0.05 -0.15 0.38 0.38      0.38 0.20 0.53  0.55  0.66      0.60 0.65        0.38
## DE_5  1444 1363 -0.01 -0.17 0.36 0.36      0.36 0.16 0.44  0.49  0.69      0.58 0.62        0.38
##       pooled.sd.0 pooled.sd.total d.sd.total
## AL_6          0.4            0.39       0.16
## AT_1          0.4            0.39       0.49
## AT_11         0.4            0.39       0.39
## AT_2          0.4            0.39       0.56
## AT_3          0.4            0.39       0.54
## AT_7          0.4            0.39       0.45
## AT_8          0.4            0.39       0.42
## AT_9          0.4            0.39       0.48
## BE_1          0.4            0.39       0.34
## BE_10         0.4            0.39       0.43
## BE_11         0.4            0.39       0.29
## BE_2          0.4            0.39       0.48
## BE_3          0.4            0.39       0.36
## BE_4          0.4            0.39       0.38
## BE_5          0.4            0.39       0.37
## BE_6          0.4            0.39       0.38
## BE_7          0.4            0.39       0.34
## BE_8          0.4            0.39       0.43
## BE_9          0.4            0.39       0.35
## BG_10         0.4            0.39       0.40
## BG_11         0.4            0.39       0.53
## BG_3          0.4            0.39       0.44
## BG_4          0.4            0.39       0.48
## BG_5          0.4            0.39       0.45
## BG_6          0.4            0.39       0.46
## BG_9          0.4            0.39       0.40
## CH_1          0.4            0.39       0.41
## CH_10         0.4            0.39       0.43
## CH_11         0.4            0.39       0.35
## CH_2          0.4            0.39       0.43
## CH_3          0.4            0.39       0.41
## CH_4          0.4            0.39       0.63
## CH_5          0.4            0.39       0.41
## CH_6          0.4            0.39       0.40
## CH_7          0.4            0.39       0.50
## CH_8          0.4            0.39       0.43
## CH_9          0.4            0.39       0.42
## CY_11         0.4            0.39       0.57
## CY_3          0.4            0.39       0.23
## CY_4          0.4            0.39       0.32
## CY_5          0.4            0.39       0.31
## CY_6          0.4            0.39       0.40
## CY_9          0.4            0.39       0.20
## CZ_1          0.4            0.39       0.70
## CZ_10         0.4            0.39       0.44
## CZ_2          0.4            0.39       0.74
## CZ_4          0.4            0.39       0.59
## CZ_5          0.4            0.39       0.59
## CZ_6          0.4            0.39       0.53
## CZ_7          0.4            0.39       0.56
## CZ_8          0.4            0.39       0.35
## CZ_9          0.4            0.39       0.48
## DE_1          0.4            0.39       0.47
## DE_11         0.4            0.39       0.45
## DE_2          0.4            0.39       0.54
## DE_3          0.4            0.39       0.50
## DE_4          0.4            0.39       0.52
## DE_5          0.4            0.39       0.41
##  [ reached 'max' / getOption("max.print") -- omitted 220 rows ]
```

``` r
# range in gender differences across folds with fold-specific SDs
range(sum_tab$D)
```

```
## [1] 0.04560813 0.71413463
```

``` r
# range in gender differences across folds with SD pooled across all folds
range(sum_tab$d.sd.total)
```

```
## [1] 0.04898065 0.73592517
```

``` r
# mean gender difference across folds with SD pooled across all folds
mean(sum_tab$d.sd.total)
```

```
## [1] 0.412898
```

``` r
# average probability of correct classification (pcc)
mean(sum_tab$pcc.total)
```

```
## [1] 0.5802061
```

``` r
# pcc for men
mean(sum_tab$pcc.1)
```

```
## [1] 0.5954792
```

``` r
# pcc for women
mean(sum_tab$pcc.0)
```

```
## [1] 0.5718254
```

``` r
# average area under the curve
mean(sum_tab$auc)
```

```
## [1] 0.6178933
```

``` r
# print smallest and largest gender differences
sum_tab[sum_tab$d.sd.total==min(sum_tab$d.sd.total),]
```

```
##       n.1 n.0         m.1         m.0      sd.1      sd.0 pooled.sd       diff          D pcc.1     pcc.0
## MK_10 538 662 0.003747672 -0.01534145 0.4044009 0.4296955 0.4185463 0.01908912 0.04560813   0.5 0.5196375
##       pcc.total       auc pooled.sd.1 pooled.sd.0 pooled.sd.total d.sd.total
## MK_10 0.5108333 0.5125535   0.3773205   0.3998527       0.3897277 0.04898065
```

``` r
sum_tab[sum_tab$d.sd.total==max(sum_tab$d.sd.total),]
```

```
##       n.1  n.0       m.1         m.0      sd.1      sd.0 pooled.sd      diff         D    pcc.1     pcc.0
## CZ_2 1089 1268 0.2423504 -0.04446001 0.4194121 0.4770804 0.4513545 0.2868104 0.6354438 0.728191 0.5205047
##      pcc.total       auc pooled.sd.1 pooled.sd.0 pooled.sd.total d.sd.total
## CZ_2 0.6164616 0.6751419   0.3773205   0.3998527       0.3897277  0.7359252
```

``` r
# correlation between men and women in male-typicality across the folds
cor(sum_tab$m.1,sum_tab$m.0)
```

```
## [1] 0.9204389
```

``` r
# correlation between men and women in male-typicality deviations across the folds
cor(sum_tab$sd.1,sum_tab$sd.0)
```

```
## [1] 0.8574715
```



# Obtain VBMT with another training sample

Limit this training sample so that there is equal number of people from all participating countries. Choose 200 men and 200 women from each country. In this case, the training data is also a lot smaller.



``` r
# rerun the analysis
value_typ2<-
  D_regularized(data=fdat,mv.vars = value.vars,
                group.var = "gndr.bin",group.values = c(1,0),
                out = T,fold = T,fold.var = "cntry",size = 200,
                pcc = T,auc=T,pred.prob = T,append.data=T)
```


## Country-year-individual variance accounted

Combine the separate frames


``` r
d1<-value_typ$preds
d2<-value_typ2$preds

d3<-left_join(x=d1,y=d2[,c("cntry","essround","idno","pred")],
              by=c("cntry","essround","idno"))
```


Separate by gender and calculate orthogonal predictors at country, country-round, and individual level by centering.


``` r
d3_men<-d3[d3$gndr.c==0.5,]
d3_women<-d3[d3$gndr.c==-0.5,]

d3_men_cntry<-
  d3_men %>%
  group_by(cntry)%>%
  summarise(y_cntry_mean=mean(pred.y,na.rm=T),
            x_cntry_mean=mean(pred.x,na.rm=T))

d3_men_cntry_essround<-
  d3_men %>%
  group_by(cntry,essround) %>%
  summarise(y_cntry_essround_mean=mean(pred.y,na.rm=T),
            x_cntry_essround_mean=mean(pred.x,na.rm=T))
```

```
## `summarise()` has regrouped the output.
## ℹ Summaries were computed grouped by cntry and essround.
## ℹ Output is grouped by cntry.
## ℹ Use `summarise(.groups = "drop_last")` to silence this message.
## ℹ Use `summarise(.by = c(cntry, essround))` for per-operation grouping (`?dplyr::dplyr_by`) instead.
```

``` r
d3_men<-
  left_join(
    x=d3_men,
    y=d3_men_cntry,
    by="cntry"
  )

d3_men<-
  left_join(
    x=d3_men,
    y=d3_men_cntry_essround,
    by=c("cntry","essround")
  )

# center individual scores around cntry-essround means
d3_men$y.crc<-d3_men$pred.y-d3_men$y_cntry_essround_mean
d3_men$x.crc<-d3_men$pred.x-d3_men$x_cntry_essround_mean

# center cntry-essround means around country means

d3_men$y_cntry_essround_mean.c<-
  d3_men$y_cntry_essround_mean-d3_men$y_cntry_mean
d3_men$x_cntry_essround_mean.c<-
  d3_men$x_cntry_essround_mean-d3_men$x_cntry_mean


d3_women_cntry<-
  d3_women %>%
  group_by(cntry)%>%
  summarise(y_cntry_mean=mean(pred.y,na.rm=T),
            x_cntry_mean=mean(pred.x,na.rm=T))

d3_women_cntry_essround<-
  d3_women %>%
  group_by(cntry,essround) %>%
  summarise(y_cntry_essround_mean=mean(pred.y,na.rm=T),
            x_cntry_essround_mean=mean(pred.x,na.rm=T))
```

```
## `summarise()` has regrouped the output.
## ℹ Summaries were computed grouped by cntry and essround.
## ℹ Output is grouped by cntry.
## ℹ Use `summarise(.groups = "drop_last")` to silence this message.
## ℹ Use `summarise(.by = c(cntry, essround))` for per-operation grouping (`?dplyr::dplyr_by`) instead.
```

``` r
d3_women<-
  left_join(
    x=d3_women,
    y=d3_women_cntry,
    by="cntry"
  )

d3_women<-
  left_join(
    x=d3_women,
    y=d3_women_cntry_essround,
    by=c("cntry","essround")
  )

# center individual scores around cntry-essround means
d3_women$y.crc<-d3_women$pred.y-d3_women$y_cntry_essround_mean
d3_women$x.crc<-d3_women$pred.x-d3_women$x_cntry_essround_mean

# center cntry-essround means around country means

d3_women$y_cntry_essround_mean.c<-d3_women$y_cntry_essround_mean-d3_women$y_cntry_mean
d3_women$x_cntry_essround_mean.c<-d3_women$x_cntry_essround_mean-d3_women$x_cntry_mean
```


### Variance accounted for men


``` r
## empty model

mod0_men<-lmer(pred.x~(1|cntry/essround),data=d3_men,REML=F,
               weights = pspwght,subset=!is.na(d3_men$y.crc),
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

summary(mod0_men)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: pred.x ~ (1 | cntry/essround)
##    Data: d3_men
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
##  Subset: !is.na(d3_men$y.crc)
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##  200484.8  200525.5 -100238.4  200476.8    195985 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -7.6845 -0.6125 -0.0016  0.5835  8.3568 
## 
## Random effects:
##  Groups         Name        Variance Std.Dev.
##  essround:cntry (Intercept) 0.003099 0.05567 
##  cntry          (Intercept) 0.008181 0.09045 
##  Residual                   0.148640 0.38554 
## Number of obs: 195989, groups:  essround:cntry, 278; cntry, 39
## 
## Fixed effects:
##             Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)  0.10130    0.01518 37.92766   6.673 6.89e-08 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

``` r
getVC(mod0_men,round=8)
```

```
##              grp        var1 var2      sdcor       vcov
## 1 essround:cntry (Intercept) <NA> 0.05567122 0.00309929
## 2          cntry (Intercept) <NA> 0.09044840 0.00818091
## 3       Residual        <NA> <NA> 0.38553844 0.14863989
```

``` r
## predictors model

mod1_men<-lmer(pred.x~y_cntry_mean+y_cntry_essround_mean.c+y.crc+(1|cntry/essround),
               data=d3_men,REML=F,weights = pspwght,subset=!is.na(d3_men$y.crc),
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

summary(mod1_men)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: pred.x ~ y_cntry_mean + y_cntry_essround_mean.c + y.crc + (1 |      cntry/essround)
##    Data: d3_men
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
##  Subset: !is.na(d3_men$y.crc)
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## -596494.2 -596422.9  298254.1 -596508.2    195982 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.2305 -0.5826  0.0272  0.6046  8.9047 
## 
## Random effects:
##  Groups         Name        Variance  Std.Dev.
##  essround:cntry (Intercept) 1.667e-05 0.004083
##  cntry          (Intercept) 7.941e-05 0.008911
##  Residual                   2.551e-03 0.050506
## Number of obs: 195989, groups:  essround:cntry, 278; cntry, 39
## 
## Fixed effects:
##                          Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)             2.085e-02  1.833e-03 3.817e+01   11.38 7.89e-14 ***
## y_cntry_mean            1.154e+00  1.716e-02 3.893e+01   67.25  < 2e-16 ***
## y_cntry_essround_mean.c 1.154e+00  5.810e-03 2.487e+02  198.68  < 2e-16 ***
## y.crc                   1.153e+00  3.444e-04 1.958e+05 3348.20  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) y_cnt_ y_c__.
## y_cntry_men -0.594              
## y_cntry_s_.  0.002  0.011       
## y.crc       -0.002  0.000  0.002
```

``` r
getVC(mod1_men,round=8)
```

```
##              grp        var1 var2      sdcor       vcov
## 1 essround:cntry (Intercept) <NA> 0.00408333 0.00001667
## 2          cntry (Intercept) <NA> 0.00891139 0.00007941
## 3       Residual        <NA> <NA> 0.05050573 0.00255083
```

``` r
# variance accounted

(data.frame(getVC(mod0_men,round=8))[,"vcov"]-data.frame(getVC(mod1_men,round=8))[,"vcov"])/
  data.frame(getVC(mod0_men,round=8))[,"vcov"]
```

```
## [1] 0.9946213 0.9902933 0.9828389
```


### Variance accounted for women


``` r
## empty model

mod0_women<-lmer(pred.x~(1|cntry/essround),data=d3_women,REML=F,weights = pspwght,
                 subset=!is.na(d3_women$y.crc),
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

summary(mod0_women)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: pred.x ~ (1 | cntry/essround)
##    Data: d3_women
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
##  Subset: !is.na(d3_women$y.crc)
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##  267512.7  267554.2 -133752.4  267504.7    234489 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.7229 -0.6151  0.0089  0.6088  7.5508 
## 
## Random effects:
##  Groups         Name        Variance Std.Dev.
##  essround:cntry (Intercept) 0.003422 0.05849 
##  cntry          (Intercept) 0.009832 0.09915 
##  Residual                   0.153819 0.39220 
## Number of obs: 234493, groups:  essround:cntry, 278; cntry, 39
## 
## Fixed effects:
##             Estimate Std. Error       df t value Pr(>|t|)   
## (Intercept) -0.05339    0.01658 37.00440  -3.221  0.00266 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

``` r
getVC(mod0_women,round=8)
```

```
##              grp        var1 var2      sdcor       vcov
## 1 essround:cntry (Intercept) <NA> 0.05849437 0.00342159
## 2          cntry (Intercept) <NA> 0.09915411 0.00983154
## 3       Residual        <NA> <NA> 0.39219720 0.15381865
```

``` r
## predictors model

mod1_women<-lmer(pred.x~y_cntry_mean+y_cntry_essround_mean.c+y.crc+(1|cntry/essround),
                 data=d3_women,REML=F,weights = pspwght,subset=!is.na(d3_women$y.crc),
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

summary(mod1_women)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: pred.x ~ y_cntry_mean + y_cntry_essround_mean.c + y.crc + (1 |      cntry/essround)
##    Data: d3_women
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
##  Subset: !is.na(d3_women$y.crc)
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## -709046.2 -708973.7  354530.1 -709060.2    234486 
## 
## Scaled residuals: 
##    Min     1Q Median     3Q    Max 
## -8.653 -0.590  0.018  0.599  8.485 
## 
## Random effects:
##  Groups         Name        Variance  Std.Dev.
##  essround:cntry (Intercept) 1.957e-05 0.004423
##  cntry          (Intercept) 9.105e-05 0.009542
##  Residual                   2.392e-03 0.048910
## Number of obs: 234493, groups:  essround:cntry, 278; cntry, 39
## 
## Fixed effects:
##                          Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)             2.027e-02  1.951e-03 4.121e+01   10.39 4.43e-13 ***
## y_cntry_mean            1.178e+00  1.727e-02 4.099e+01   68.25  < 2e-16 ***
## y_cntry_essround_mean.c 1.170e+00  5.922e-03 2.457e+02  197.53  < 2e-16 ***
## y.crc                   1.159e+00  3.011e-04 2.343e+05 3850.75  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) y_cnt_ y_c__.
## y_cntry_men  0.589              
## y_cntry_s_.  0.011  0.007       
## y.crc       -0.001  0.000  0.003
```

``` r
getVC(mod1_women,round=8)
```

```
##              grp        var1 var2      sdcor       vcov
## 1 essround:cntry (Intercept) <NA> 0.00442330 0.00001957
## 2          cntry (Intercept) <NA> 0.00954197 0.00009105
## 3       Residual        <NA> <NA> 0.04891001 0.00239219
```

``` r
# variance accounted

(data.frame(getVC(mod0_women,round=8))[,"vcov"]-data.frame(getVC(mod1_women,round=8))[,"vcov"])/
  data.frame(getVC(mod0_women,round=8))[,"vcov"]
```

```
## [1] 0.9942804 0.9907390 0.9844480
```


# Session information


``` r
s<-sessionInfo()
print(s,locale=FALSE)
```

```
## R version 4.6.0 (2026-04-24 ucrt)
## Platform: x86_64-w64-mingw32/x64
## Running under: Windows 11 x64 (build 22631)
## 
## Matrix products: default
##   LAPACK version 3.12.1
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] tibble_3.3.1          apaTables_2.0.8       stringr_1.6.0         tidyr_1.3.2          
##  [5] r2mlm_0.3.8           nlme_3.1-169          Hmisc_5.2-5           ggpubr_0.6.3         
##  [9] metafor_5.0-1         numDeriv_2016.8-1.1   metadat_1.6-0         lmerTest_3.2-1       
## [13] ggflags_0.0.4         finalfit_1.1.0        ggplot2_4.0.3         MetBrewer_0.2.0      
## [17] vjihelpers_0.0.0.9000 emmeans_2.0.3         lme4_2.0-1            Matrix_1.7-5         
## [21] dplyr_1.2.1           rio_1.3.0             multid_1.0.2          knitr_1.51           
## [25] rmarkdown_2.31       
## 
## loaded via a namespace (and not attached):
##  [1] Rdpack_2.6.6       pROC_1.19.0.1      gridExtra_2.3      rlang_1.2.0        magrittr_2.0.5    
##  [6] otel_0.2.0         rockchalk_1.8.164  compiler_4.6.0     vctrs_0.7.3        pkgconfig_2.0.3   
## [11] shape_1.4.6.1      fastmap_1.2.0      backports_1.5.1    nloptr_2.2.1       purrr_1.2.2       
## [16] xfun_0.57          glmnet_5.0         jomo_2.7-6         cachem_1.1.0       kutils_1.73       
## [21] jsonlite_2.0.0     pan_1.9            broom_1.0.13       cluster_2.1.8.2    R6_2.6.1          
## [26] bslib_0.10.0       stringi_1.8.7      RColorBrewer_1.1-3 car_3.1-5          boot_1.3-32       
## [31] rpart_4.1.27       jquerylib_0.1.4    estimability_1.5.1 Rcpp_1.1.1-1.1     iterators_1.0.14  
## [36] base64enc_0.1-6    R.utils_2.13.0     splines_4.6.0      nnet_7.3-20        tidyselect_1.2.1  
## [41] rstudioapi_0.18.0  abind_1.4-8        yaml_2.3.12        codetools_0.2-20   lattice_0.22-9    
## [46] plyr_1.8.9         withr_3.0.2        S7_0.2.2           coda_0.19-4.1      evaluate_1.0.5    
## [51] foreign_0.8-91     survival_3.8-6     zip_2.3.3          pillar_1.11.1      carData_3.0-6     
## [56] mice_3.19.0        checkmate_2.3.4    foreach_1.5.2      reformulas_0.4.4   generics_0.1.4    
## [61] mathjaxr_2.0-0     scales_1.4.0       minqa_1.2.8        xtable_1.8-8       glue_1.8.1        
## [66] tools_4.6.0        data.table_1.18.4  openxlsx_4.2.8.1   ggsignif_0.6.4     forcats_1.0.1     
## [71] mvtnorm_1.3-7      grid_4.6.0         rbibutils_2.4.1    colorspace_2.1-2   htmlTable_2.5.0   
## [76] Formula_1.2-5      cli_3.6.6          gtable_0.3.6       R.methodsS3_1.8.2  rstatix_0.7.3     
## [81] sass_0.4.10        digest_0.6.39      htmlwidgets_1.6.4  farver_2.1.2       htmltools_0.5.9   
## [86] R.oo_1.27.1        lifecycle_1.0.5    mitml_0.4-5        MASS_7.3-65
```

