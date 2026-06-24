---
title: "Analysis for Examining the Gender Equality Paradox in Values Using a New Measure of Value-Based Gender Typicality"
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
library(multilevel)
library(misty)
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

![](Analysis_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

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

# Data for the analysis

* use the testing partition of the data set (diff_dat)
* basic descriptives


``` r
diff_dat<-
  value_typ$preds

#diff_dat$age.c<-as.numeric(diff_dat$age.c)
#diff_dat$rlgdgr.c<-as.numeric(diff_dat$rlgdgr.c)
#diff_dat$eduyrs.c<-as.numeric(diff_dat$eduyrs.c)

# data exclusions (include countries with at least two rounds of data)
n_rounds<-diff_dat %>% group_by(cntry) %>%
  summarise(n_unique_essround = n_distinct(essround))

# join number of rounds to data
diff_dat<-
  left_join(x=diff_dat,
            y=n_rounds,
            by="cntry")

# filter only countries with more than 1 round and non-missing gender variable
diff_dat<-diff_dat %>%
  dplyr::filter(n_unique_essround>1 &
                  !is.na(gndr))

# cross-tabulate sample sizes and ESS rounds by country
table(diff_dat$essround,diff_dat$cntry)
```

```
##     
##        AT   BE   BG   CH   CY   CZ   DE   DK   EE   ES   FI   FR   GB   GR   HR   HU   IE   IL   IS   IT
##   1  2054 1630    0 1824    0 1008 2619 1270    0 1512 1563 1155 1598 2351    0 1434 1716 2079    0    0
##   2  1998 1571    0 1910    0 2357 2640 1258 1748 1423 1501 1499 1664 2199    0 1260  987    0  324    0
##   3  2148 1596 1095 1580  778    0 2684 1261 1266 1647 1449 1783 2153    0    0 1262 1389    0    0    0
##   4     0 1554 1944 1553 1010 1786 2532 1381 1446 2362 1701 1867 2111 1863 1230 1230 1557 2182    0    0
##   5     0 1499 2171 1291  853 2135 2807 1364 1593 1681 1449 1523 2174 2469 1401 1273 2200 2012    0    0
##   6     0 1662 1979 1283  910 1773 2735 1421 2145 1671 1958 1760 2061    0    0 1768 2416 2178  539  709
##   7  1595 1567    0 1321    0 1662 2806 1283 1836 1707 1850 1702 2031    0    0 1320 2180 2151    0    0
##   8  1793 1559    0 1304    0 2052 2621    0 1807 1729 1703 1857 1742    0    0 1258 2546 2166  641 2331
##   9  2277 1556 1726 1317  573 2143 2128 1354 1699 1419 1535 1782 1983    0 1581 1443 1989    0  644 2460
##   10    0 1134 2497 1305    0 2169    0    0 1338    0 1361 1751  931 2568 1364 1616 1551    0  686 2373
##   11 2114 1377 2018 1168  467    0 2181    0 1082 1633 1324 1545 1329 2545 1348 1917 1785  693  625 2583
##     
##        LT   LV   ME   NL   NO   PL   PT   RS   RU   SE   SI   SK   TR   UA
##   1     0    0    0 2137 1619 1865 1282    0    0 1482 1288    0    0    0
##   2     0    0    0 1658 1375 1483 1824    0    0 1478 1184 1225 1590 1696
##   3     0    0    0 1660 1350 1485 1982    0 2139 1404 1265 1511    0 1685
##   4     0 1770    0 1524 1191 1396 2137    0 2246 1356 1057 1589 2105 1566
##   5  1432    0    0 1601 1330 1519 1939    0 2357 1263 1169 1603    0 1579
##   6  1908    0    0 1628 1410 1666 1938    0 2229 1638 1044 1627    0 1864
##   7  2041    0    0 1623 1223 1394 1042    0    0 1561  989    0    0    0
##   8  1879    0    0 1469 1330 1475 1054    0 2174 1326 1095    0    0    0
##   9  1477  691  988 1457 1196 1243  845 1769    0 1310 1107  861    0    0
##   10 1406    0 1048 1266 1208    0 1627    0    0    0 1032 1195    0    0
##   11 1137 1025 1381 1477 1118 1223 1166 1311    0 1004 1027 1214    0 2389
```

``` r
# range of sample sizes
range(table(diff_dat$cntry))
```

```
## [1]  3080 25753
```

``` r
# value-based gender-typicality histogram
hist(diff_dat$pred)
```

![](Analysis_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

``` r
# scale/standardize with SD pooled across all country-time-gender folds
FM_pooled_sd<-value_typ$D[1,"pooled.sd.total"]
FM_pooled_sd
```

```
## [1] 0.3897277
```

``` r
# standardized
diff_dat$FM.z<-diff_dat$pred/FM_pooled_sd
hist(diff_dat$FM.z)
```

![](Analysis_files/figure-html/unnamed-chunk-10-2.png)<!-- -->

``` r
# recode time to start from 2002=0

diff_dat$year<-
  case_when(
    diff_dat$essround==1~2002,  
    diff_dat$essround==2~2004,
    diff_dat$essround==3~2006,
    diff_dat$essround==4~2008,
    diff_dat$essround==5~2010,
    diff_dat$essround==6~2012,
    diff_dat$essround==7~2014,
    diff_dat$essround==8~2016,
    diff_dat$essround==9~2018,
    diff_dat$essround==10~2020,
    diff_dat$essround==11~2023
  )

diff_dat$year.c<-diff_dat$year-2002
```

## Add GEI-variable

* GEI = Reverse-coded Gender Inequality Index (GII)



``` r
GII<-read.csv("../data/GII.csv")

# Obtain GII only for ESS countries
GII_in_ESS_d<-
  GII %>%
  filter(ISO2 %in% unique(diff_dat$cntry))

# Make long format data of GII, so that country x year has own rows

long_GII_in_ESS_d <- GII_in_ESS_d %>%
  tidyr::pivot_longer(
    cols = starts_with("gii_") & !matches("avg"),  # <-- exclude the "_avg" columns
    names_to = "year",
    values_to = "gii",
    names_prefix = "gii_"
  ) %>%
  mutate(year = as.integer(year)) %>%
  filter(year >= 2002 & year <= 2023) %>%
  mutate(gei = 1 - gii) %>%
  dplyr::select(ISO2, iso3, country, year, gii, gei, gii_2002_2023_avg, gei_2002_2023_avg)


# describe gei average
psych::describe(GII_in_ESS_d$gei_2002_2023_avg)
```

```
##    vars  n mean   sd median trimmed  mad  min  max range  skew kurtosis   se
## X1    1 33 0.87 0.07   0.87    0.87 0.06 0.63 0.96  0.34 -1.06     1.37 0.01
```

``` r
# one is missing, see which one
GII_in_ESS_d[is.na(GII_in_ESS_d$gei_2002_2023_avg),"country"]
```

```
## [1] "Ukraine"
```

``` r
# get means and SDs for standardizing
gei_mean<-mean(GII_in_ESS_d$gei_2002_2023_avg,na.rm=T)
gei_sd<-sd(GII_in_ESS_d$gei_2002_2023_avg,na.rm=T)

# standardize 
# year specific scores
long_GII_in_ESS_d$gei.z<-(long_GII_in_ESS_d$gei-gei_mean)/gei_sd

# get the country means to a separate frame

GII_country_means<-
  long_GII_in_ESS_d %>% group_by(ISO2) %>%
  summarise(gei.cm=mean(gei,na.rm=T),
            gei.z.cm=mean(gei.z,na.rm=T))


long_GII_in_ESS_d<-left_join(
  x=long_GII_in_ESS_d,
  y=GII_country_means,
  by="ISO2"
)


# center both the raw score and standardized scores
long_GII_in_ESS_d$gei.cmc<-(long_GII_in_ESS_d$gei-long_GII_in_ESS_d$gei.cm)
long_GII_in_ESS_d$gei.z.cmc<-(long_GII_in_ESS_d$gei.z-long_GII_in_ESS_d$gei.z.cm)

#View(long_GII_in_ESS_d)

# averages across 2002-2023
#long_GII_in_ESS_d$gei_2002_2023_avg.z<-(long_GII_in_ESS_d$gei_2002_2023_avg-gei_mean)/gei_sd

#long_GII_in_ESS_d$gei.cmc<-long_GII_in_ESS_d$gei-long_GII_in_ESS_d$gei_2002_2023_avg

# add year to ESS data-frame


# link datasets by country and year

diff_dat<-
  left_join(
    x=diff_dat,
    y=long_GII_in_ESS_d[,c("ISO2","year","gei","gei.cm","gei.cmc",
                           "gei.z","gei.z.cm","gei.z.cmc")],
    by=c("cntry"="ISO2","year")
)
```


## Add GDI-variable

* GDI = Gender Development Index


``` r
GDI<-read.csv("../data/GDI.csv")

# Obtain GDI only for ESS countries
GDI_in_ESS_d<-
  GDI %>%
  filter(ISO2 %in% unique(diff_dat$cntry))

# Make long format data of GDI, so that country x year has own rows

long_GDI_in_ESS_d <- GDI_in_ESS_d %>%
  pivot_longer(
    cols = starts_with("gdi_") & !matches("avg"),  # <-- exclude the "_avg" columns
    names_to = "year",
    values_to = "gdi",
    names_prefix = "gdi_"
  ) %>%
  mutate(year = as.integer(year)) %>%
  filter(year >= 2002 & year <= 2023) %>%
  dplyr::select(ISO2, iso3, country, year, gdi, gdi_2002_2023_avg)

# describe gdi average
psych::describe(GDI_in_ESS_d$gdi_2002_2023_avg)
```

```
##    vars  n mean   sd median trimmed  mad min  max range skew kurtosis se
## X1    1 34 0.98 0.03   0.98    0.98 0.02 0.9 1.03  0.13 -0.3     1.28  0
```

``` r
# get means and SDs for standardizing
gdi_mean<-mean(GDI_in_ESS_d$gdi_2002_2023_avg,na.rm=T)
gdi_sd<-sd(GDI_in_ESS_d$gdi_2002_2023_avg,na.rm=T)

# standardize 
# year specific scores
long_GDI_in_ESS_d$gdi.z<-(long_GDI_in_ESS_d$gdi-gdi_mean)/gdi_sd

# get the country means to a separate frame

GDI_country_means<-
  long_GDI_in_ESS_d %>% group_by(ISO2) %>%
  summarise(gdi.cm=mean(gdi,na.rm=T),
            gdi.z.cm=mean(gdi.z,na.rm=T))

long_GDI_in_ESS_d<-left_join(
  x=long_GDI_in_ESS_d,
  y=GDI_country_means,
  by="ISO2"
)


# center both the raw score and standardized scores
long_GDI_in_ESS_d$gdi.cmc<-(long_GDI_in_ESS_d$gdi-long_GDI_in_ESS_d$gdi.cm)
long_GDI_in_ESS_d$gdi.z.cmc<-(long_GDI_in_ESS_d$gdi.z-long_GDI_in_ESS_d$gdi.z.cm)

# link datasets by country and year

diff_dat<-
  left_join(
    x=diff_dat,
    y=long_GDI_in_ESS_d[,c("ISO2","year","gdi","gdi.cm","gdi.cmc",
                           "gdi.z","gdi.z.cm","gdi.z.cmc")],
    by=c("cntry"="ISO2","year")
  )
```

## Add GDP-variable

* log_GDP = Logarithm of GDP per capita, PPP (constant 2017 international $)



``` r
GDP<-read.csv("../data/GDP_processed.csv")

# Obtain GDP only for ESS countries
GDP_in_ESS_d<-
  GDP %>%
  filter(ISO2 %in% unique(diff_dat$cntry))

# First, rename those X2000..YR2000. style columns to simple "gdp_2000" etc.

GDP_in_ESS_d <- GDP_in_ESS_d %>%
  rename_with(
    ~ str_replace_all(., "^X(\\d{4})..YR\\1\\.$", "gdp_\\1"),
    starts_with("X")
  )

# GDP_in_ESS_d
# Make long format data of GDP, so that country x year has own rows

long_GDP_in_ESS_d <- GDP_in_ESS_d %>%
  pivot_longer(
    cols = starts_with("gdp_") & !matches("avg"),  # <-- exclude the "_avg" columns
    names_to = "year",
    values_to = "gdp",
    names_prefix = "gdp_"
  ) %>%
  mutate(year = as.integer(year)) %>%
  filter(year >= 2002 & year <= 2023) %>%
  dplyr::select(ISO2, Country.Name, year, gdp, gdp_2002_2023_avg,log_gdp_2002_2023_avg) %>%
  mutate(log_gdp=log(gdp))

#View(long_GDP_in_ESS_d)

# describe gdp average
psych::describe(GDP_in_ESS_d$gdp_2002_2023_avg)
```

```
##    vars  n     mean      sd  median  trimmed     mad      min      max    range skew kurtosis      se
## X1    1 34 43905.99 17369.7 40201.8 42845.95 15827.5 16384.61 85341.85 68957.24 0.48    -0.58 2978.88
```

``` r
psych::describe(GDP_in_ESS_d$log_gdp_2002_2023_avg)
```

```
##    vars  n  mean   sd median trimmed  mad min   max range  skew kurtosis   se
## X1    1 34 10.61 0.41   10.6   10.62 0.44 9.7 11.35  1.65 -0.27    -0.71 0.07
```

``` r
# get means and SDs for standardizing
log_gdp_mean<-mean(GDP_in_ESS_d$log_gdp_2002_2023_avg,na.rm=T)
log_gdp_sd<-sd(GDP_in_ESS_d$log_gdp_2002_2023_avg,na.rm=T)

# standardize 
# year specific scores
long_GDP_in_ESS_d$log_gdp.z<-
  (long_GDP_in_ESS_d$log_gdp-log_gdp_mean)/log_gdp_sd

# get the country means to a separate frame

GDP_country_means<-
  long_GDP_in_ESS_d %>% group_by(ISO2) %>%
  summarise(gdp.cm=mean(gdp,na.rm=T),
            #gdp.z.cm=mean(gdp.z,na.rm=T),
            log_gdp.cm=mean(log_gdp,na.rm=T),
            log_gdp.z.cm=mean(log_gdp.z,na.rm=T))

long_GDP_in_ESS_d<-left_join(
  x=long_GDP_in_ESS_d,
  y=GDP_country_means,
  by="ISO2"
)


# center both the raw score and standardized scores
long_GDP_in_ESS_d$log_gdp.cmc<-(long_GDP_in_ESS_d$log_gdp-long_GDP_in_ESS_d$log_gdp.cm)
long_GDP_in_ESS_d$log_gdp.z.cmc<-(long_GDP_in_ESS_d$log_gdp.z-long_GDP_in_ESS_d$log_gdp.z.cm)

# link datasets by country and year

diff_dat<-
  left_join(
    x=diff_dat,
    y=long_GDP_in_ESS_d[,c("ISO2","year","log_gdp","log_gdp.cm","log_gdp.cmc",
                           "log_gdp.z","log_gdp.z.cm","log_gdp.z.cmc")],
    by=c("cntry"="ISO2","year")
  )
```

## Add GGGI-variable

* GGGI = Global Gender Gap Index



``` r
GGGI<-import("../data/GGGI.xlsx")

# check if all ESS countries are in the log_GGGI data
table(unique(diff_dat$cntry) %in% GGGI$ISO2)
```

```
## 
## TRUE 
##   34
```

``` r
# Obtain GGGI only for ESS countries
GGGI_in_ESS_d<-
  GGGI %>%
  filter(ISO2 %in% unique(diff_dat$cntry))

# Make long format data of GGGI, so that country x year has own rows

long_GGGI_in_ESS_d <- GGGI_in_ESS_d %>%
  pivot_longer(
    cols = starts_with("gggi_") & !matches("avg"),  # <-- exclude the "_avg" columns
    names_to = "year",
    values_to = "gggi",
    names_prefix = "gggi_"
  ) %>%
  mutate(year = as.integer(year)) %>%
  filter(year >= 2002 & year <= 2023) %>%
  dplyr::select(ISO2, cname, year, gggi, GGGI_2002_2023_avg)

# describe gggi average
psych::describe(GGGI_in_ESS_d$GGGI_2002_2023_avg)
```

```
##    vars  n mean   sd median trimmed  mad  min  max range skew kurtosis   se
## X1    1 34 0.74 0.05   0.73    0.73 0.04 0.61 0.86  0.25  0.4     0.29 0.01
```

``` r
# get means and SDs for standardizing
gggi_mean<-mean(GGGI_in_ESS_d$GGGI_2002_2023_avg,na.rm=T)
gggi_sd<-sd(GGGI_in_ESS_d$GGGI_2002_2023_avg,na.rm=T)

# standardize 
# year specific scores
long_GGGI_in_ESS_d$gggi.z<-(long_GGGI_in_ESS_d$gggi-gggi_mean)/gggi_sd

# get the country means to a separate frame

GGGI_country_means<-
  long_GGGI_in_ESS_d %>% group_by(ISO2) %>%
  summarise(gggi.cm=mean(gggi,na.rm=T),
            gggi.z.cm=mean(gggi.z,na.rm=T))

long_GGGI_in_ESS_d<-left_join(
  x=long_GGGI_in_ESS_d,
  y=GGGI_country_means,
  by="ISO2"
)

# center both the raw score and standardized scores
long_GGGI_in_ESS_d$gggi.cmc<-
  (long_GGGI_in_ESS_d$gggi-long_GGGI_in_ESS_d$gggi.cm)

long_GGGI_in_ESS_d$gggi.z.cmc<-
  (long_GGGI_in_ESS_d$gggi.z-long_GGGI_in_ESS_d$gggi.z.cm)

# link datasets by country and year

diff_dat<-
  left_join(
    x=diff_dat,
    y=long_GGGI_in_ESS_d[,c("ISO2","year","gggi","gggi.cm","gggi.cmc",
                           "gggi.z","gggi.z.cm","gggi.z.cmc")],
    by=c("cntry"="ISO2","year")
  )
```

## Country-year-gender dataframe

Frame with one row for country-year-gender


``` r
diff_dat_cntry_year<-
  diff_dat %>% 
  group_by(cntry,year,gndr.c) %>%
  dplyr::summarise(n=n(),
                   n_wt=sum(pspwght),
                   FM.z.wt=weighted.mean(x=FM.z,w=pspwght),
                   FM.z=mean(FM.z),
                   gei.z=mean(gei.z,na.rm=T),
                   gei.z.cm=mean(gei.z.cm,na.rm=T),
                   gei.z.cmc=mean(gei.z.cmc,na.rm=T),
                   gdi.z=mean(gdi.z,na.rm=T),
                   gdi.z.cm=mean(gdi.z.cm,na.rm=T),
                   gdi.z.cmc=mean(gdi.z.cmc,na.rm=T),
                   gggi.z=mean(gggi.z,na.rm=T),
                   gggi.z.cm=mean(gggi.z.cm,na.rm=T),
                   gggi.z.cmc=mean(gggi.z.cmc,na.rm=T),
                   log_gdp.z=mean(log_gdp.z,na.rm=T),
                   log_gdp.z.cm=mean(log_gdp.z.cm,na.rm=T),
                   log_gdp.z.cmc=mean(log_gdp.z.cmc,na.rm=T),
                   gei=mean(gei,na.rm=T),
                   gei.cm=mean(gei.cm,na.rm=T),
                   gei.cmc=mean(gei.cmc,na.rm=T),
                   gdi=mean(gdi,na.rm=T),
                   gdi.cm=mean(gdi.cm,na.rm=T),
                   gdi.cmc=mean(gdi.cmc,na.rm=T),
                   gggi=mean(gggi,na.rm=T),
                   gggi.cm=mean(gggi.cm,na.rm=T),
                   gggi.cmc=mean(gggi.cmc,na.rm=T),
                   log_gdp=mean(log_gdp,na.rm=T),
                   log_gdp.cm=mean(log_gdp.cm,na.rm=T),
                   log_gdp.cmc=mean(log_gdp.cmc,na.rm=T))
```

```
## `summarise()` has regrouped the output.
## ℹ Summaries were computed grouped by cntry, year, and gndr.c.
## ℹ Output is grouped by cntry and year.
## ℹ Use `summarise(.groups = "drop_last")` to silence this message.
## ℹ Use `summarise(.by = c(cntry, year, gndr.c))` for per-operation grouping (`?dplyr::dplyr_by`) instead.
```


# Descriptive statistics

## Country-specific descriptives


``` r
# sample sizes from weights

cntry_n_frame<-
  diff_dat %>% group_by(cntry) %>%
  summarise('n ESS rounds' = mean(n_unique_essround),
            n=round(sum(pspwght),0))

# value-based male-typicality

cntry_FM_frame<-
  diff_dat %>% group_by(cntry,essround) %>%
  summarise('FM M' = weighted.mean(x=FM.z,w=pspwght),
            'FM SD' = sqrt(wtd.var(FM.z,pspwght))) %>%
  group_by(cntry) %>%
  summarise('FM M' = mean(x=`FM M`),
            'FM SD'= mean(x=`FM SD`))
```

```
## `summarise()` has regrouped the output.
## ℹ Summaries were computed grouped by cntry and essround.
## ℹ Output is grouped by cntry.
## ℹ Use `summarise(.groups = "drop_last")` to silence this message.
## ℹ Use `summarise(.by = c(cntry, essround))` for per-operation grouping (`?dplyr::dplyr_by`) instead.
```

``` r
cntry_FM_women_frame<-
  diff_dat %>%
  filter(gndr.c==-0.5) %>%
  group_by(cntry,essround) %>%
  summarise('FM M' = weighted.mean(x=FM.z,w=pspwght),
            'FM SD' = sqrt(wtd.var(FM.z,pspwght))) %>%
  group_by(cntry) %>%
  summarise('FM M Women' = mean(x=`FM M`),
            'FM SD Women'= mean(x=`FM SD`))
```

```
## `summarise()` has regrouped the output.
## ℹ Summaries were computed grouped by cntry and essround.
## ℹ Output is grouped by cntry.
## ℹ Use `summarise(.groups = "drop_last")` to silence this message.
## ℹ Use `summarise(.by = c(cntry, essround))` for per-operation grouping (`?dplyr::dplyr_by`) instead.
```

``` r
cntry_FM_men_frame<-
  diff_dat %>%
  filter(gndr.c==0.5) %>%
  group_by(cntry,essround) %>%
  summarise('FM M' = weighted.mean(x=FM.z,w=pspwght),
            'FM SD' = sqrt(wtd.var(FM.z,pspwght))) %>%
  group_by(cntry) %>%
  summarise('FM M Men' = mean(x=`FM M`),
            'FM SD Men'= mean(x=`FM SD`))
```

```
## `summarise()` has regrouped the output.
## ℹ Summaries were computed grouped by cntry and essround.
## ℹ Output is grouped by cntry.
## ℹ Use `summarise(.groups = "drop_last")` to silence this message.
## ℹ Use `summarise(.by = c(cntry, essround))` for per-operation grouping (`?dplyr::dplyr_by`) instead.
```

``` r
# link n and FM datasets

desc_frame<-
  left_join(
    x=cntry_n_frame,
    y=cntry_FM_frame,
    by="cntry"
  )

desc_frame<-
  left_join(
    x=desc_frame,
    y=cntry_FM_women_frame,
    by="cntry"
  )

desc_frame<-
  left_join(
    x=desc_frame,
    y=cntry_FM_men_frame,
    by="cntry"
  )

# Add country-specific differences
desc_frame$D<-desc_frame$`FM M Men`-desc_frame$`FM M Women`

desc_frame
```

```
## # A tibble: 34 × 10
##    cntry `n ESS rounds`     n  `FM M` `FM SD` `FM M Women` `FM SD Women` `FM M Men` `FM SD Men`     D
##    <chr>          <dbl> <dbl>   <dbl>   <dbl>        <dbl>         <dbl>      <dbl>       <dbl> <dbl>
##  1 AT                 7 14032  0.0730   1.04        -0.173         1.02      0.341        0.993 0.514
##  2 BE                11 16681 -0.0327   0.919       -0.217         0.923     0.162        0.873 0.380
##  3 BG                 7 13437  0.0337   1.09        -0.186         1.12      0.277        1.01  0.464
##  4 CH                11 15914 -0.0572   0.877       -0.272         0.862     0.169        0.835 0.441
##  5 CY                 6  4580 -0.298    1.04        -0.461         1.06     -0.118        0.986 0.342
##  6 CZ                 9 17085  0.402    1.14         0.123         1.15      0.708        1.05  0.585
##  7 DE                10 25753 -0.233    0.997       -0.476         0.976     0.0226       0.954 0.499
##  8 DK                 8 10590  0.0961   0.969       -0.142         0.927     0.341        0.948 0.484
##  9 EE                10 15958 -0.102    1.05        -0.337         1.03      0.190        0.988 0.527
## 10 ES                10 16795 -0.478    1.02        -0.646         1.03     -0.302        0.976 0.344
## # ℹ 24 more rows
```

``` r
# add country-level measures

desc_frame<-
  left_join(
    x=desc_frame,
    y=GII_country_means[,c("ISO2","gei.cm")],
    by=c("cntry"="ISO2")
  )

desc_frame<-
  left_join(
    x=desc_frame,
    y=GGGI_country_means[,c("ISO2","gggi.cm")],
    by=c("cntry"="ISO2")
  )

desc_frame<-
  left_join(
    x=desc_frame,
    y=GDI_country_means[,c("ISO2","gdi.cm")],
    by=c("cntry"="ISO2")
  )


desc_frame<-
  left_join(
    x=desc_frame,
    y=GDP_country_means[,c("ISO2","gdp.cm")],
    by=c("cntry"="ISO2")
  )

desc_frame<-left_join(
  x=desc_frame,
  y=ISO[,c("ISO2","CLDR")],
  by=c("cntry"="ISO2")
)

# finalize the formatted table

cntry_desc_tbl<-
  desc_frame %>%
  dplyr::select(
    Country = CLDR,
    `n ESS rounds`,
    n,
    `FM M`, `FM SD`,
    `FM M Women`, `FM SD Women`,
    `FM M Men`, `FM SD Men`,
    D,
    GEI = gei.cm,
    GGGI = gggi.cm,
    GDI = gdi.cm,
    GDP = gdp.cm
  ) %>%
  mutate(
    across(
      .cols = -c(Country, `n ESS rounds`, n, GDP) & where(is.numeric),
      .fns  = ~ round_tidy(.x, 2)
    )
  ) %>%
  mutate(
    across(
      .cols = GDP,
      .fns  = ~ round_tidy(.x, 0)
    )
  )
print(cntry_desc_tbl,n=35)
```

```
## # A tibble: 34 × 14
##    Country     `n ESS rounds`     n `FM M` `FM SD` `FM M Women` `FM SD Women` `FM M Men` `FM SD Men` D    
##    <chr>                <dbl> <dbl> <chr>  <chr>   <chr>        <chr>         <chr>      <chr>       <chr>
##  1 Austria                  7 14032 0.07   1.04    -0.17        1.02          0.34       0.99        0.51 
##  2 Belgium                 11 16681 -0.03  0.92    -0.22        0.92          0.16       0.87        0.38 
##  3 Bulgaria                 7 13437 0.03   1.09    -0.19        1.12          0.28       1.01        0.46 
##  4 Switzerland             11 15914 -0.06  0.88    -0.27        0.86          0.17       0.84        0.44 
##  5 Cyprus                   6  4580 -0.30  1.04    -0.46        1.06          -0.12      0.99        0.34 
##  6 Czechia                  9 17085 0.40   1.14    0.12         1.15          0.71       1.05        0.59 
##  7 Germany                 10 25753 -0.23  1.00    -0.48        0.98          0.02       0.95        0.50 
##  8 Denmark                  8 10590 0.10   0.97    -0.14        0.93          0.34       0.95        0.48 
##  9 Estonia                 10 15958 -0.10  1.05    -0.34        1.03          0.19       0.99        0.53 
## 10 Spain                   10 16795 -0.48  1.02    -0.65        1.03          -0.30      0.98        0.34 
## 11 Finland                 11 17394 -0.27  1.04    -0.53        1.01          0.03       0.98        0.56 
## 12 France                  11 18277 -0.29  1.05    -0.51        1.03          -0.05      1.02        0.46 
## 13 UK                      11 20697 -0.14  1.00    -0.36        0.98          0.10       0.96        0.47 
## 14 Greece                   6 13982 0.08   1.05    -0.07        1.10          0.24       0.97        0.31 
## 15 Croatia                  5  6895 -0.24  1.07    -0.42        1.09          -0.02      1.01        0.40 
## 16 Hungary                 11 15877 0.19   1.02    0.04         1.04          0.36       0.98        0.32 
## 17 Ireland                 11 20405 -0.01  1.06    -0.19        1.08          0.19       1.01        0.38 
## 18 Israel                   7 13446 0.30   1.01    0.16         1.03          0.45       0.96        0.29 
## 19 Iceland                  6  3446 -0.31  0.94    -0.54        0.91          -0.07      0.91        0.47 
## 20 Italy                    5 10458 0.12   0.98    -0.03        0.98          0.28       0.95        0.31 
## 21 Lithuania                7 11541 0.50   1.09    0.31         1.08          0.76       1.03        0.45 
## 22 Latvia                   3  3458 0.19   1.10    -0.05        1.07          0.52       1.05        0.57 
## 23 Montenegro               3  3422 0.35   1.19    0.17         1.22          0.54       1.11        0.37 
## 24 Netherlands             11 17478 0.22   0.90    -0.00        0.88          0.46       0.86        0.46 
## 25 Norway                  11 14317 0.08   0.95    -0.14        0.93          0.30       0.90        0.44 
## 26 Poland                  10 14740 -0.01  0.98    -0.23        0.98          0.24       0.91        0.48 
## 27 Portugal                11 16831 0.05   0.97    -0.09        0.98          0.22       0.93        0.31 
## 28 Serbia                   2  3096 -0.32  1.11    -0.55        1.11          -0.09      1.07        0.46 
## 29 Russia                   5 11105 0.30   1.14    0.17         1.17          0.48       1.08        0.31 
## 30 Sweden                  10 14058 0.04   0.99    -0.20        0.96          0.28       0.96        0.48 
## 31 Slovenia                11 12254 0.15   0.87    -0.05        0.86          0.37       0.83        0.42 
## 32 Slovakia                 8 10905 0.29   1.11    0.08         1.13          0.52       1.04        0.44 
## 33 Turkey                   2  3704 0.41   0.86    0.29         0.90          0.53       0.80        0.24 
## 34 Ukraine                  6 10778 0.18   1.20    0.03         1.22          0.39       1.14        0.37 
## # ℹ 4 more variables: GEI <chr>, GGGI <chr>, GDI <chr>, GDP <chr>
```

``` r
export(cntry_desc_tbl,"../results/cntry_desc_tbl.xlsx",overwrite=T)

# export also the entire data.frame without rounding for other analyses
export(desc_frame,"../results/cntry_lvl_data.xlsx",overwrite=T)
```

## Country-level correlations table


``` r
cor_frame<-
  desc_frame %>%
  dplyr::select(
    VBMT=`FM M`,
    VBMT_Women=`FM M Women`,
    VBMT_Men=`FM M Men`,
    D = D,
    GEI = gei.cm,
    GGGI = gggi.cm,
    GDI = gdi.cm,
    GDP = gdp.cm
  ) %>%
  mutate(
    log_GDP=log(GDP)
  ) %>%
  dplyr::select(-GDP)

apa.cor.table(
  data=cor_frame,
  filename = "../results/CorTable1.doc",
  #table.number = NA,
  show.conf.interval = FALSE,
  show.sig.stars = F,
  landscape = F
)
```

```
## The ability to suppress reporting of reporting confidence intervals has been deprecated in this version.
## The function argument show.conf.interval will be removed in a later version.
```

```
## 
## 
## Means, standard deviations, and correlations with confidence intervals
##  
## 
##   Variable      M     SD   1            2            3            4           5           6          
##   1. VBMT       0.04  0.24                                                                           
##                                                                                                      
##   2. VBMT_Women -0.16 0.26 .99                                                                       
##                            [.98, .99]                                                                
##                                                                                                      
##   3. VBMT_Men   0.26  0.24 .98          .94                                                          
##                            [.96, .99]   [.88, .97]                                                   
##                                                                                                      
##   4. D          0.42  0.09 -.19         -.34         .01                                             
##                            [-.49, .16]  [-.60, .00]  [-.33, .34]                                     
##                                                                                                      
##   5. GEI        0.87  0.07 -.34         -.40         -.28         .41                                
##                            [-.61, .01]  [-.65, -.07] [-.57, .07]  [.07, .66]                         
##                                                                                                      
##   6. GGGI       0.74  0.05 -.44         -.51         -.35         .52         .73                    
##                            [-.68, -.12] [-.72, -.21] [-.62, -.01] [.23, .73]  [.52, .86]             
##                                                                                                      
##   7. GDI        0.98  0.03 .08          .05          .18          .33         .07         .19        
##                            [-.26, .41]  [-.29, .38]  [-.17, .48]  [-.01, .60] [-.28, .41] [-.16, .50]
##                                                                                                      
##   8. log_GDP    10.61 0.41 -.27         -.31         -.24         .23         .72         .62        
##                            [-.55, .08]  [-.58, .03]  [-.54, .11]  [-.11, .53] [.50, .85]  [.36, .79] 
##                                                                                                      
##   7          
##              
##              
##              
##              
##              
##              
##              
##              
##              
##              
##              
##              
##              
##              
##              
##              
##              
##              
##              
##              
##   -.18       
##   [-.49, .17]
##              
## 
## Note. M and SD are used to represent mean and standard deviation, respectively.
## Values in square brackets indicate the 95% confidence interval.
## The confidence interval is a plausible range of population correlations 
## that could have caused the sample correlation (Cumming, 2014).
## 
```
## Intraclass coefficient for gender


``` r
fit_gndr_cntry<-
  glmer(gndr.bin~(1|cntry),
        data=diff_dat,family=binomial(link="logit"))
summary(fit_gndr_cntry)
```

```
## Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
##  Family: binomial  ( logit )
## Formula: gndr.bin ~ (1 | cntry)
##    Data: diff_dat
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##  600811.9  600833.8 -300403.9  600807.9    437741 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -1.0504 -0.9241 -0.8096  1.0774  1.4053 
## 
## Random effects:
##  Groups Name        Variance Std.Dev.
##  cntry  (Intercept) 0.03316  0.1821  
## Number of obs: 437743, groups:  cntry, 34
## 
## Fixed effects:
##             Estimate Std. Error z value Pr(>|z|)    
## (Intercept) -0.18939    0.03144  -6.024  1.7e-09 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

``` r
tau=VarCorr(fit_gndr_cntry)$cntry[1]
(ICC_gndr=tau/(tau+pi^2/3))
```

```
## [1] 0.009979834
```


# Analysis

Following preregistration: https://osf.io/7cags?view_only=f3e97a78271e46bfafb9e20ac8d35bb1 

## mod0: Random intercept model


``` r
mod0<-lmer(FM.z~(1|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod0)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: FM.z ~ (1 | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##   1326480   1326513   -663237   1326474    437740 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.3857 -0.6196 -0.0001  0.5972  8.6856 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  cntry    (Intercept) 0.06183  0.2486  
##  Residual             1.06397  1.0315  
## Number of obs: 437743, groups:  cntry, 34
## 
## Fixed effects:
##             Estimate Std. Error       df t value Pr(>|t|)
## (Intercept)  0.04819    0.04268 33.95905   1.129    0.267
```

``` r
getVC(mod0)
```

```
##        grp        var1 var2 sdcor vcov
## 1    cntry (Intercept) <NA>  0.25 0.06
## 2 Residual        <NA> <NA>  1.03 1.06
```

``` r
r2mlm(mod0,bargraph = F)
```

```
## $Decompositions
##                      total within between
## fixed, within   0.00000000      0      NA
## fixed, between  0.00000000     NA       0
## slope variation 0.00000000      0      NA
## mean variation  0.05491835     NA       1
## sigma2          0.94508165      1      NA
## 
## $R2s
##          total within between
## f1  0.00000000      0      NA
## f2  0.00000000     NA       0
## v   0.00000000      0      NA
## m   0.05491835     NA       1
## f   0.00000000     NA      NA
## fv  0.00000000      0      NA
## fvm 0.05491835     NA      NA
```

``` r
# calculate reliability of the group means

multilevel.icc(diff_dat[, c("FM.z", "cntry")],
               cluster = c("cntry"),type = "2")
```

```
## [1] 0.9987213
```

``` r
# for men
multilevel.icc(diff_dat[diff_dat$gndr.c==0.5, c("FM.z", "cntry")],
               cluster = c("cntry"),type = "2")
```

```
## [1] 0.9975226
```

``` r
# for women
multilevel.icc(diff_dat[diff_dat$gndr.c==-0.5, c("FM.z", "cntry")],
               cluster = c("cntry"),type = "2")
```

```
## [1] 0.9978874
```

``` r
## reliability across three levels

multilevel.icc(diff_dat[, c("FM.z", "cntry", "cntry_time")],
               cluster = c("cntry", "cntry_time"),type = "2")
```

```
##      cntry cntry_time 
##  0.9516789  0.9710500
```

``` r
# for men
multilevel.icc(diff_dat[diff_dat$gndr.c==0.5,
                        c("FM.z", "cntry","cntry_time")],
               cluster = c("cntry","cntry_time"),type = "2")
```

```
##      cntry cntry_time 
##  0.9526160  0.9426394
```

``` r
# for women
multilevel.icc(diff_dat[diff_dat$gndr.c==0.5, 
                        c("FM.z", "cntry","cntry_time")],
               cluster = c("cntry","cntry_time"),type = "2")
```

```
##      cntry cntry_time 
##  0.9526160  0.9426394
```

## mod1: Gender fixed effect


``` r
mod1<-lmer(FM.z~gndr.c+(1|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod1)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: FM.z ~ gndr.c + (1 | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1307163.2 1307207.2 -653577.6 1307155.2    437739 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.2198 -0.6155  0.0013  0.5996  8.5137 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  cntry    (Intercept) 0.06323  0.2514  
##  Residual             1.01803  1.0090  
## Number of obs: 437743, groups:  cntry, 34
## 
## Fixed effects:
##              Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept) 5.880e-02  4.316e-02 3.396e+01   1.362    0.182    
## gndr.c      4.286e-01  3.050e-03 4.377e+05 140.540   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##        (Intr)
## gndr.c 0.002
```

``` r
getFE(mod1,round=3)
```

```
##              Est.    SE         df       t     p     LL    UL
## (Intercept) 0.059 0.043     33.962   1.362 0.182 -0.029 0.147
## gndr.c      0.429 0.003 437710.895 140.540 0.000  0.423 0.435
```

``` r
getVC(mod1)
```

```
##        grp        var1 var2 sdcor vcov
## 1    cntry (Intercept) <NA>  0.25 0.06
## 2 Residual        <NA> <NA>  1.01 1.02
```

``` r
r2mlm(mod1,bargraph = F)
```

```
## $Decompositions
##                      total
## fixed           0.04045118
## slope variation 0.00000000
## mean variation  0.05610897
## sigma2          0.90343985
## 
## $R2s
##          total
## f   0.04045118
## v   0.00000000
## m   0.05610897
## fv  0.04045118
## fvm 0.09656015
```

## mod2: Gender fixed and random effect

* Include random effect correlation by default


``` r
mod2<-lmer(FM.z~gndr.c+(gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: FM.z ~ gndr.c + (gndr.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1306540.8 1306606.7 -653264.4 1306528.8    437737 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.3265 -0.6158  0.0016  0.6006  8.4531 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr  
##  cntry    (Intercept) 0.063185 0.25137        
##           gndr.c      0.007376 0.08589  -0.16 
##  Residual             1.016329 1.00813        
## Number of obs: 437743, groups:  cntry, 34
## 
## Fixed effects:
##             Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)  0.05878    0.04315 33.96104   1.362    0.182    
## gndr.c       0.42171    0.01515 33.33705  27.831   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##        (Intr)
## gndr.c -0.156
```

``` r
getFE(mod2,round=3)
```

```
##              Est.    SE     df      t     p     LL    UL
## (Intercept) 0.059 0.043 33.961  1.362 0.182 -0.029 0.146
## gndr.c      0.422 0.015 33.337 27.831 0.000  0.391 0.453
```

``` r
getVC(mod2)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.25 0.06
## 2    cntry      gndr.c   <NA>  0.09 0.01
## 3    cntry (Intercept) gndr.c -0.16 0.00
## 4 Residual        <NA>   <NA>  1.01 1.02
```

``` r
r2mlm(mod2,bargraph = F)
```

```
## $Decompositions
##                      total
## fixed           0.03919028
## slope variation 0.00162556
## mean variation  0.05640719
## sigma2          0.90277697
## 
## $R2s
##          total
## f   0.03919028
## v   0.00162556
## m   0.05640719
## fv  0.04081584
## fvm 0.09722303
```

``` r
anova(mod1,mod2)
```

```
## Data: diff_dat
## Models:
## mod1: FM.z ~ gndr.c + (1 | cntry)
## mod2: FM.z ~ gndr.c + (gndr.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod1    4 1307163 1307207 -653578   1307155                         
## mod2    6 1306541 1306607 -653264   1306529 626.47  2  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

``` r
lvl2_var_cond_lvl1(mod2,lvl1.var = "gndr.c",lvl1.values = c(0.5,-0.5))
```

```
##   lvl1.value lvl2.cond.var lvl2.cond.sd
## 1        0.5    0.06155836    0.2481096
## 2       -0.5    0.06849968    0.2617244
```

``` r
reliability_dms(model=mod2,diff_var="gndr.c",
              diff_var_values=c(0.5,-0.5),var = "FM.z",
              group_var = "cntry")
```

```
##              r11              r22              r12              sd1              sd2           sd_d12 
##       0.99704802       0.99748750       0.94219550       0.25664124       0.26971233       0.09040591 
##               m1               m2            m_d12 reliability_dmsa 
##       0.23929955      -0.17282443       0.41212398       0.95384905
```

* Test for random effect correlation


``` r
mod2_norecov<-lmer(FM.z~gndr.c+(gndr.c||cntry),data=diff_dat,REML=F,weights = pspwght,
                   control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2_norecov)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: FM.z ~ gndr.c + (gndr.c || cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1306539.6 1306594.5 -653264.8 1306529.6    437738 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.3259 -0.6158  0.0016  0.6007  8.4532 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  cntry    (Intercept) 0.063196 0.25139 
##  cntry.1  gndr.c      0.007371 0.08585 
##  Residual             1.016329 1.00813 
## Number of obs: 437743, groups:  cntry, 34
## 
## Fixed effects:
##             Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)  0.05879    0.04315 33.95985   1.362    0.182    
## gndr.c       0.42172    0.01515 33.34057  27.842   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##        (Intr)
## gndr.c 0.000
```

``` r
getFE(mod2_norecov,round=3)
```

```
##              Est.    SE     df      t     p     LL    UL
## (Intercept) 0.059 0.043 33.960  1.362 0.182 -0.029 0.146
## gndr.c      0.422 0.015 33.341 27.842 0.000  0.391 0.453
```

``` r
getVC(mod2_norecov)
```

```
##        grp        var1 var2 sdcor vcov
## 1    cntry (Intercept) <NA>  0.25 0.06
## 2  cntry.1      gndr.c <NA>  0.09 0.01
## 3 Residual        <NA> <NA>  1.01 1.02
```

``` r
anova(mod2_norecov,mod2)
```

```
## Data: diff_dat
## Models:
## mod2_norecov: FM.z ~ gndr.c + (gndr.c || cntry)
## mod2: FM.z ~ gndr.c + (gndr.c | cntry)
##              npar     AIC     BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)
## mod2_norecov    5 1306540 1306595 -653265   1306530                     
## mod2            6 1306541 1306607 -653264   1306529 0.8243  1     0.3639
```


## mod2 with Gender-equality index (GEI)


``` r
mod2_GEI<-lmer(FM.z~gndr.c+gei.z.cm+gndr.c:gei.z.cm+
                 (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2_GEI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: FM.z ~ gndr.c + gei.z.cm + gndr.c:gei.z.cm + (gndr.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1268297.2 1268385.0 -634140.6 1268281.2    426956 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.3751 -0.6169  0.0022  0.6023  8.4949 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr  
##  cntry    (Intercept) 0.056052 0.2368         
##           gndr.c      0.006464 0.0804   -0.02 
##  Residual             1.006181 1.0031         
## Number of obs: 426964, groups:  cntry, 33
## 
## Fixed effects:
##                 Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)      0.05443    0.04125 32.96494   1.319   0.1961    
## gndr.c           0.42221    0.01445 32.54025  29.216   <2e-16 ***
## gei.z.cm        -0.09311    0.04191 33.02177  -2.221   0.0333 *  
## gndr.c:gei.z.cm  0.03452    0.01487 34.27073   2.321   0.0264 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c g.z.cm
## gndr.c      -0.017              
## gei.z.cm    -0.001  0.000       
## gndr.c:g.z.  0.000 -0.020 -0.017
```

``` r
getFE(mod2_GEI,round=3)
```

```
##                   Est.    SE     df      t     p     LL     UL
## (Intercept)      0.054 0.041 32.965  1.319 0.196 -0.030  0.138
## gndr.c           0.422 0.014 32.540 29.216 0.000  0.393  0.452
## gei.z.cm        -0.093 0.042 33.022 -2.221 0.033 -0.178 -0.008
## gndr.c:gei.z.cm  0.035 0.015 34.271  2.321 0.026  0.004  0.065
```

``` r
getVC(mod2_GEI)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.24 0.06
## 2    cntry      gndr.c   <NA>  0.08 0.01
## 3    cntry (Intercept) gndr.c -0.02 0.00
## 4 Residual        <NA>   <NA>  1.00 1.01
```

``` r
r2mlm(mod2_GEI,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.045390777
## slope variation 0.001440176
## mean variation  0.050330231
## sigma2          0.902838817
## 
## $R2s
##           total
## f   0.045390777
## v   0.001440176
## m   0.050330231
## fv  0.046830952
## fvm 0.097161183
```

### Deconstructed associations


``` r
t1<-Sys.time()
ddsc_mod2_GEI<-
  ddsc_ml(model = mod2_GEI,
          predictor = "gei.z.cm",
          moderator = "gndr.c",moderator_values = c(-0.5,0.5),
          re_cov_test = T)
t2<-Sys.time()
t2-t1
```

```
## Time difference of 49.93352 secs
```

``` r
round(ddsc_mod2_GEI$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.068        0.262        1.016     1.085 0.063   6999.706 0.998   0.998
## 2        0.5         0.062        0.248        1.016     1.078 0.057   5875.088 0.997   0.997
```

``` r
round(ddsc_mod2_GEI$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gei.z.cm gei.z.cm_scaled
## means_y1           0.236 0.260    1.000           1.000    0.942           0.942   -0.269          -0.269
## means_y1_scaled    0.885 0.975    1.000           1.000    0.942           0.942   -0.269          -0.269
## means_y2          -0.176 0.273    0.942           0.942    1.000           1.000   -0.362          -0.362
## means_y2_scaled   -0.661 1.025    0.942           0.942    1.000           1.000   -0.362          -0.362
## gei.z.cm           0.000 1.000   -0.269          -0.269   -0.362          -0.362    1.000           1.000
## gei.z.cm_scaled    0.000 1.000   -0.269          -0.269   -0.362          -0.362    1.000           1.000
## diff_score         0.412 0.092    0.028           0.028   -0.310          -0.310    0.313           0.313
## diff_score_scaled  1.546 0.344    0.028           0.028   -0.310          -0.310    0.313           0.313
##                   diff_score diff_score_scaled
## means_y1               0.028             0.028
## means_y1_scaled        0.028             0.028
## means_y2              -0.310            -0.310
## means_y2_scaled       -0.310            -0.310
## gei.z.cm               0.313             0.313
## gei.z.cm_scaled        0.313             0.313
## diff_score             1.000             1.000
## diff_score_scaled      1.000             1.000
```

``` r
round(ddsc_mod2_GEI$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                      -0.376 0.162 34.271  -2.321   0.026   -0.705   -0.047
## w_11                         -0.110 0.043 33.054  -2.585   0.014   -0.197   -0.024
## w_21                         -0.076 0.042 33.048  -1.787   0.083   -0.162    0.011
## r_xy1                        -0.425 0.164 33.054  -2.585   0.014   -0.759   -0.091
## r_xy2                        -0.278 0.155 33.048  -1.787   0.083   -0.594    0.038
## b_11                         -0.414 0.160 33.054  -2.585   0.014   -0.740   -0.088
## b_21                         -0.285 0.159 33.048  -1.787   0.083   -0.609    0.039
## main_effect                  -0.093 0.042 33.022  -2.221   0.033   -0.178   -0.008
## moderator_effect              0.422 0.014 32.540  29.216   0.000    0.393    0.452
## interaction                   0.035 0.015 34.271   2.321   0.026    0.004    0.065
## q_b11_b21                    -0.148    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                  -0.168    NA     NA      NA      NA       NA       NA
## cross_over_point            -12.232    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.059 0.044 33.126  -1.324   0.194   -0.149    0.031
## interaction_vs_main_bscale   -0.220 0.166 33.126  -1.324   0.194   -0.558    0.118
## interaction_vs_main_rscale   -0.204 0.158 33.133  -1.288   0.207   -0.527    0.118
## dadas                        -0.152 0.085 33.048  -1.787   0.958   -0.324    0.021
## dadas_bscale                 -0.569 0.319 33.048  -1.787   0.958   -1.217    0.079
## dadas_rscale                 -0.555 0.311 33.048  -1.787   0.958   -1.188    0.077
## abs_diff                      0.035 0.015 34.271   2.321   0.013    0.004    0.065
## abs_sum                       0.186 0.084 33.022   2.221   0.017    0.016    0.357
## abs_diff_bscale               0.130 0.056 34.271   2.321   0.013    0.016    0.243
## abs_sum_bscale                0.699 0.315 33.022   2.221   0.017    0.059    1.339
## abs_diff_rscale               0.147 0.057 34.246   2.602   0.007    0.032    0.262
## abs_sum_rscale                0.703 0.315 33.022   2.232   0.016    0.062    1.343
```

``` r
round(ddsc_mod2_GEI$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.003 -0.161  0.824  1.000  0.364
```

``` r
# Structural path model
d_GEI<-ddsc_mod2_GEI$ddsc_sem_fit$data

ddsc_sem_GEI<-
  ddsc_sem(data=d_GEI,x = "gei.z.cm",y1="means_y2",
           y2="means_y1",q_sesoi = .10)
round(ddsc_sem_GEI$results,3)
```

```
##                                     est    se       z pvalue ci.lower ci.upper
## r_xy1_y2                         -0.313 0.165  -1.895  0.058   -0.637    0.011
## r_xy1                            -0.362 0.162  -2.227  0.026   -0.680   -0.043
## r_xy2                            -0.269 0.168  -1.607  0.108   -0.598    0.059
## b_11                             -0.370 0.166  -2.227  0.026   -0.696   -0.044
## b_21                             -0.263 0.163  -1.607  0.108   -0.583    0.058
## b_10                             -0.661 0.164  -4.038  0.000   -0.982   -0.340
## b_20                              0.885 0.161   5.498  0.000    0.569    1.200
## res_cov_y1_y2                     0.818 0.208   3.936  0.000    0.411    1.225
## diff_b10_b20                     -1.546 0.056 -27.570  0.000   -1.656   -1.436
## diff_b11_b21                     -0.108 0.057  -1.895  0.058   -0.219    0.004
## diff_rxy1_rxy2                   -0.092 0.057  -1.612  0.107   -0.204    0.020
## q_b11_b21                        -0.120 0.066  -1.826  0.068   -0.249    0.009
## q_rxy1_rxy2                      -0.102 0.064  -1.609  0.108   -0.227    0.022
## cross_over_point                -14.330 7.582  -1.890  0.059  -29.190    0.530
## sum_b11_b21                      -0.633 0.325  -1.949  0.051   -1.270    0.004
## main_effect                      -0.317 0.162  -1.949  0.051   -0.635    0.002
## interaction_vs_main_effect       -0.209 0.169  -1.233  0.218   -0.540    0.123
## diff_abs_b11_abs_b21              0.108 0.057   1.895  0.058   -0.004    0.219
## abs_diff_b11_b21                  0.108 0.057   1.895  0.029   -0.004    0.219
## abs_sum_b11_b21                   0.633 0.325   1.949  0.026   -0.004    1.270
## dadas                            -0.525 0.327  -1.607  0.946   -1.166    0.115
## q_r_equivalence                   0.002 0.064   0.038  0.515       NA       NA
## q_b_equivalence                   0.020 0.066   0.305  0.620       NA       NA
## cross_over_point_equivalence     14.330 7.582   1.890  0.971       NA       NA
## cross_over_point_minimal_effect  14.330 7.582   1.890  0.029       NA       NA
```

``` r
round(ddsc_sem_GEI$variance_test,3)
```

```
##             est    se      z pvalue ci.lower ci.upper
## cov_y1y2  0.912 0.232  3.939  0.000    0.458    1.366
## var_y1    1.018 0.251  4.062  0.000    0.527    1.510
## var_y2    0.921 0.227  4.062  0.000    0.477    1.366
## var_diff  0.097 0.116  0.838  0.402   -0.130    0.324
## var_ratio 1.105 0.129  8.549  0.000    0.852    1.359
## cor_y1y2  0.942 0.020 47.928  0.000    0.903    0.980
```

``` r
## random intercept mlm with double-entries for each country (men and women)

d_GEI_long <- d_GEI %>%
  # move row names into a column
  rownames_to_column("cntry") %>%
  # pivot only the means_y1 / means_y2 columns
  pivot_longer(
    cols = c(means_y1, means_y2),
    names_to = "y",
    values_to = "means"
  ) %>%
  # create a sgender column based on y1/y2
  mutate(
    gndr.c = case_when(
      y == "means_y1" ~ -0.5,
      y == "means_y2" ~ 0.5
    )
  ) %>%
  dplyr::select(cntry,gei.z.cm,means,gndr.c)


ddsc_mod2_GEI_ri<-
  ddsc_ml(data=data.frame(d_GEI_long),predictor = "gei.z.cm",
          moderator = "gndr.c",
        DV = "means",lvl2_unit = "cntry",
        moderator_values = c(-0.5,0.5))
round(ddsc_mod2_GEI_ri$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                      -0.313 0.171 31.000  -1.836   0.076   -0.661    0.035
## w_11                         -0.099 0.045 32.904  -2.178   0.037   -0.191   -0.006
## w_21                         -0.070 0.045 32.904  -1.544   0.132   -0.162    0.022
## r_xy1                        -0.362 0.166 32.904  -2.178   0.037   -0.699   -0.024
## r_xy2                        -0.269 0.175 32.904  -1.544   0.132   -0.625    0.086
## b_11                         -0.371 0.170 32.904  -2.178   0.037   -0.717   -0.024
## b_21                         -0.263 0.170 32.904  -1.544   0.132   -0.609    0.084
## main_effect                  -0.084 0.045 31.000  -1.889   0.068   -0.175    0.007
## moderator_effect              0.412 0.015 31.000  26.722   0.000    0.381    0.444
## interaction                   0.029 0.016 31.000   1.836   0.076   -0.003    0.061
## q_b11_b21                    -0.120    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                  -0.102    NA     NA      NA      NA       NA       NA
## cross_over_point            -14.330    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.056 0.047 38.510  -1.175   0.247   -0.151    0.040
## interaction_vs_main_bscale   -0.209 0.178 38.510  -1.175   0.247   -0.568    0.151
## interaction_vs_main_rscale   -0.223 0.186 37.996  -1.202   0.237   -0.600    0.153
## dadas                        -0.140 0.091 32.904  -1.544   0.934   -0.324    0.045
## dadas_bscale                 -0.525 0.340 32.904  -1.544   0.934   -1.218    0.167
## dadas_rscale                 -0.539 0.349 32.904  -1.544   0.934   -1.249    0.171
## abs_diff                      0.029 0.016 31.000   1.836   0.038   -0.003    0.061
## abs_sum                       0.169 0.089 31.000   1.889   0.034   -0.013    0.351
## abs_diff_bscale               0.108 0.059 31.000   1.836   0.038   -0.012    0.228
## abs_sum_bscale                0.633 0.335 31.000   1.889   0.034   -0.050    1.317
## abs_diff_rscale               0.092 0.059 32.264   1.551   0.065   -0.029    0.213
## abs_sum_rscale                0.631 0.335 31.001   1.881   0.035   -0.053    1.315
```

``` r
# country-time multilevel model


mod2_GEI_cntry_year<-
  lmer(FM.z.wt~gndr.c+gei.z.cm+gndr.c:gei.z.cm+
                 (gndr.c|cntry),data=diff_dat_cntry_year,REML=F,
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2_GEI_cntry_year)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: FM.z.wt ~ gndr.c + gei.z.cm + gndr.c:gei.z.cm + (gndr.c | cntry)
##    Data: diff_dat_cntry_year
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##    -364.1    -329.8     190.0    -380.1       526 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -3.6850 -0.6315 -0.0045  0.5830  4.6042 
## 
## Random effects:
##  Groups   Name        Variance  Std.Dev. Corr 
##  cntry    (Intercept) 0.0506441 0.22504       
##           gndr.c      0.0002488 0.01577  0.01 
##  Residual             0.0230735 0.15190       
## Number of obs: 534, groups:  cntry, 33
## 
## Fixed effects:
##                 Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)      0.04385    0.03986 32.81393   1.100   0.2793    
## gndr.c           0.42428    0.01381 33.52650  30.715   <2e-16 ***
## gei.z.cm        -0.08539    0.04091 34.14917  -2.087   0.0444 *  
## gndr.c:gei.z.cm  0.03272    0.01582 38.66263   2.068   0.0454 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c g.z.cm
## gndr.c       0.003              
## gei.z.cm    -0.013  0.000       
## gndr.c:g.z.  0.000 -0.223  0.003
```

``` r
getFE(mod2_GEI_cntry_year,round=3)
```

```
##                   Est.    SE     df      t     p     LL     UL
## (Intercept)      0.044 0.040 32.814  1.100 0.279 -0.037  0.125
## gndr.c           0.424 0.014 33.527 30.715 0.000  0.396  0.452
## gei.z.cm        -0.085 0.041 34.149 -2.087 0.044 -0.169 -0.002
## gndr.c:gei.z.cm  0.033 0.016 38.663  2.068 0.045  0.001  0.065
```

``` r
getVC(mod2_GEI_cntry_year)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.23 0.05
## 2    cntry      gndr.c   <NA>  0.02 0.00
## 3    cntry (Intercept) gndr.c  0.01 0.00
## 4 Residual        <NA>   <NA>  0.15 0.02
```

``` r
r2mlm(mod2_GEI,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.045390777
## slope variation 0.001440176
## mean variation  0.050330231
## sigma2          0.902838817
## 
## $R2s
##           total
## f   0.045390777
## v   0.001440176
## m   0.050330231
## fv  0.046830952
## fvm 0.097161183
```

``` r
ddsc_mod2_GEI_cntry_year<-
  ddsc_ml(model = mod2_GEI_cntry_year,
          predictor = "gei.z.cm",
          moderator = "gndr.c",moderator_values = c(-0.5,0.5),
          re_cov_test = T)

round(ddsc_mod2_GEI_cntry_year$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.059        0.243        0.023     0.082 0.720      8.029 0.996   0.954
## 2        0.5         0.055        0.234        0.023     0.077 0.704      8.029 0.996   0.950
```

``` r
round(ddsc_mod2_GEI_cntry_year$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gei.z.cm gei.z.cm_scaled
## means_y1           0.255 0.245    1.000           1.000    0.940           0.940   -0.277          -0.277
## means_y1_scaled    1.011 0.972    1.000           1.000    0.940           0.940   -0.277          -0.277
## means_y2          -0.168 0.260    0.940           0.940    1.000           1.000   -0.401          -0.401
## means_y2_scaled   -0.666 1.028    0.940           0.940    1.000           1.000   -0.401          -0.401
## gei.z.cm           0.000 1.000   -0.277          -0.277   -0.401          -0.401    1.000           1.000
## gei.z.cm_scaled    0.000 1.000   -0.277          -0.277   -0.401          -0.401    1.000           1.000
## diff_score         0.424 0.089    0.016           0.016   -0.326          -0.326    0.407           0.407
## diff_score_scaled  1.677 0.350    0.016           0.016   -0.326          -0.326    0.407           0.407
##                   diff_score diff_score_scaled
## means_y1               0.016             0.016
## means_y1_scaled        0.016             0.016
## means_y2              -0.326            -0.326
## means_y2_scaled       -0.326            -0.326
## gei.z.cm               0.407             0.407
## gei.z.cm_scaled        0.407             0.407
## diff_score             1.000             1.000
## diff_score_scaled      1.000             1.000
```

``` r
round(ddsc_mod2_GEI_cntry_year$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                      -0.370 0.179 38.663  -2.068   0.045   -0.731   -0.008
## w_11                         -0.102 0.042 34.530  -2.443   0.020   -0.186   -0.017
## w_21                         -0.069 0.042 34.549  -1.656   0.107   -0.154    0.016
## r_xy1                        -0.415 0.170 34.530  -2.443   0.020   -0.759   -0.070
## r_xy2                        -0.266 0.161 34.549  -1.656   0.107   -0.592    0.060
## b_11                         -0.403 0.165 34.530  -2.443   0.020   -0.738   -0.068
## b_21                         -0.273 0.165 34.549  -1.656   0.107   -0.609    0.062
## main_effect                  -0.085 0.041 34.149  -2.087   0.044   -0.169   -0.002
## moderator_effect              0.424 0.014 33.527  30.715   0.000    0.396    0.452
## interaction                   0.033 0.016 38.663   2.068   0.045    0.001    0.065
## q_b11_b21                    -0.147    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                  -0.169    NA     NA      NA      NA       NA       NA
## cross_over_point            -12.967    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.053 0.044 35.495  -1.200   0.238   -0.142    0.036
## interaction_vs_main_bscale   -0.209 0.174 35.495  -1.200   0.238   -0.561    0.144
## interaction_vs_main_rscale   -0.192 0.165 35.590  -1.160   0.254   -0.527    0.144
## dadas                        -0.138 0.083 34.549  -1.656   0.947   -0.307    0.031
## dadas_bscale                 -0.547 0.330 34.549  -1.656   0.947   -1.218    0.124
## dadas_rscale                 -0.532 0.321 34.549  -1.656   0.947   -1.185    0.120
## abs_diff                      0.033 0.016 38.663   2.068   0.023    0.001    0.065
## abs_sum                       0.171 0.082 34.149   2.087   0.022    0.005    0.337
## abs_diff_bscale               0.130 0.063 38.663   2.068   0.023    0.003    0.256
## abs_sum_bscale                0.676 0.324 34.149   2.087   0.022    0.018    1.335
## abs_diff_rscale               0.149 0.063 38.878   2.346   0.012    0.020    0.277
## abs_sum_rscale                0.681 0.324 34.149   2.099   0.022    0.022    1.340
```

``` r
round(ddsc_mod2_GEI_cntry_year$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.002 -0.306  0.377  1.000  0.539
```

### Comparisons across models


``` r
# full multilevel model (individuals as entries)
round(ddsc_mod2_GEI$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df t.ratio p.value ci.lower ci.upper
## r_xy1             -0.4248 0.1643 33.0535 -2.5854  0.0143  -0.7591  -0.0905
## r_xy2             -0.2777 0.1554 33.0484 -1.7870  0.0831  -0.5938   0.0384
## b_11              -0.4142 0.1602 33.0535 -2.5854  0.0143  -0.7401  -0.0883
## b_21              -0.2846 0.1593 33.0484 -1.7870  0.0831  -0.6087   0.0394
## main_effect       -0.0931 0.0419 33.0218 -2.2215  0.0333  -0.1784  -0.0078
## moderator_effect   0.4222 0.0145 32.5402 29.2163  0.0000   0.3928   0.4516
## interaction        0.0345 0.0149 34.2707  2.3207  0.0264   0.0043   0.0647
## q_b11_b21         -0.1479     NA      NA      NA      NA       NA       NA
```

``` r
# country-level structural path model (country-gender means as entries)
round(ddsc_sem_GEI$results[c("r_xy1","r_xy2","b_11","b_21","q_b11_b21","diff_b11_b21"),],4)
```

```
##                  est     se       z pvalue ci.lower ci.upper
## r_xy1        -0.3615 0.1623 -2.2274 0.0259  -0.6796  -0.0434
## r_xy2        -0.2694 0.1676 -1.6070 0.1080  -0.5980   0.0592
## b_11         -0.3705 0.1663 -2.2274 0.0259  -0.6964  -0.0445
## b_21         -0.2626 0.1634 -1.6070 0.1080  -0.5828   0.0577
## q_b11_b21    -0.1201 0.0658 -1.8260 0.0678  -0.2490   0.0088
## diff_b11_b21 -0.1079 0.0569 -1.8946 0.0581  -0.2195   0.0037
```

``` r
# multilevel model (country-gender means as entries)
round(ddsc_mod2_GEI_ri$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE     df t.ratio p.value ci.lower ci.upper
## r_xy1             -0.3615 0.1660 32.904 -2.1779  0.0367  -0.6993  -0.0238
## r_xy2             -0.2694 0.1745 32.904 -1.5437  0.1322  -0.6245   0.0857
## b_11              -0.3706 0.1701 32.904 -2.1779  0.0367  -0.7168  -0.0244
## b_21              -0.2627 0.1701 32.904 -1.5437  0.1322  -0.6089   0.0836
## main_effect       -0.0844 0.0447 31.000 -1.8892  0.0683  -0.1754   0.0067
## moderator_effect   0.4121 0.0154 31.000 26.7219  0.0000   0.3806   0.4435
## interaction        0.0288 0.0157 31.000  1.8363  0.0759  -0.0032   0.0607
## q_b11_b21         -0.1201     NA     NA      NA      NA       NA       NA
```

``` r
# multilevel model (country-year-gender means as entries)
round(ddsc_mod2_GEI_cntry_year$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df t.ratio p.value ci.lower ci.upper
## r_xy1             -0.4146 0.1697 34.5300 -2.4433  0.0198  -0.7593  -0.0699
## r_xy2             -0.2660 0.1606 34.5487 -1.6560  0.1068  -0.5923   0.0602
## b_11              -0.4030 0.1650 34.5300 -2.4433  0.0198  -0.7381  -0.0680
## b_21              -0.2734 0.1651 34.5487 -1.6560  0.1068  -0.6088   0.0619
## main_effect       -0.0854 0.0409 34.1492 -2.0874  0.0444  -0.1685  -0.0023
## moderator_effect   0.4243 0.0138 33.5265 30.7153  0.0000   0.3962   0.4524
## interaction        0.0327 0.0158 38.6626  2.0682  0.0454   0.0007   0.0647
## q_b11_b21         -0.1467     NA      NA      NA      NA       NA       NA
```


### Bootstrap and equivalence test

Takes a lot of time


``` r
t1<-Sys.time()
mod2_GEI_booted_fixef <-
  lme4::bootMer(
    x = mod2_GEI,
    FUN = lme4::fixef,
    nsim = 1000,
    use.u = FALSE,
    seed = 12345,
    type = c("parametric"),
    verbose = FALSE
  )
t2<-Sys.time()
t2-t1
```

```
## Time difference of 1.755237 hours
```



``` r
# obtain all the bootstrap estimates
mod2_GEI_boot_est <- data.frame(mod2_GEI_booted_fixef$t)

# calculate estimates
mod2_GEI_boot_est$w11<-mod2_GEI_boot_est$gei.z.cm+(-0.5)*mod2_GEI_boot_est$gndr.c.gei.z.cm
mod2_GEI_boot_est$w21<-mod2_GEI_boot_est$gei.z.cm+(0.5)*mod2_GEI_boot_est$gndr.c.gei.z.cm
mod2_GEI_boot_est$b11<-mod2_GEI_boot_est$w11/ddsc_mod2_GEI$SDs["SD_pooled"]
mod2_GEI_boot_est$b21<-mod2_GEI_boot_est$w21/ddsc_mod2_GEI$SDs["SD_pooled"]
mod2_GEI_boot_est$r_xy1<-mod2_GEI_boot_est$w11/ddsc_mod2_GEI$SDs["SD_y1"]
mod2_GEI_boot_est$r_xy2<-mod2_GEI_boot_est$w21/ddsc_mod2_GEI$SDs["SD_y2"]
mod2_GEI_boot_est$q_b<-atanh(mod2_GEI_boot_est$b11)-atanh(mod2_GEI_boot_est$b21)
mod2_GEI_boot_est$q<-atanh(mod2_GEI_boot_est$r_xy1)-atanh(mod2_GEI_boot_est$r_xy2)

# Calculate bootstrap summary statistics
mod2_GEI_boot_results <- t(as.data.frame(sapply(
  mod2_GEI_boot_est,
  function(x) {
    c(
      Estimate = mean(x, na.rm = TRUE),
      SE = stats::sd(x, na.rm = TRUE),
      stats::quantile(x, c((1 - .95) / 2,
                           1 - (1 - .95) / 2), na.rm = TRUE)
    )
  }
)))

mod2_GEI_boot_results
```

```
##                    Estimate         SE         2.5%        97.5%
## X.Intercept.     0.05347353 0.04055668 -0.026123644  0.133920406
## gndr.c           0.42296415 0.01471536  0.394396143  0.452190662
## gei.z.cm        -0.09233569 0.04180300 -0.171823287 -0.006922419
## gndr.c.gei.z.cm  0.03409613 0.01545983  0.001876137  0.062828927
## w11             -0.10938375 0.04222327 -0.189827389 -0.019443241
## w21             -0.07528762 0.04279814 -0.156601294  0.012974383
## b11             -0.41049887 0.15845685 -0.712390331 -0.072967222
## b21             -0.28254182 0.16061422 -0.587698374  0.048690682
## r_xy1           -0.42104073 0.16252612 -0.730684941 -0.074841064
## r_xy2           -0.27564044 0.15669105 -0.573343225  0.047501361
## q_b             -0.15316719 0.07420722 -0.308034726 -0.007101674
## q               -0.17618353 0.08068671 -0.344251198 -0.024550158
```

``` r
# equivalence test for q_b
tost_z(est=mod2_GEI_boot_results["q_b","Estimate"],
       se=mod2_GEI_boot_results["q_b","SE"],
       low=-.10,
       high=.10, alpha = 0.05)
```

```
## $estimate
## [1] -0.1531672
## 
## $se
## [1] 0.07420722
## 
## $low_bound
## [1] -0.1
## 
## $high_bound
## [1] 0.1
## 
## $alpha
## [1] 0.05
## 
## $z_low
## [1] -0.7164692
## 
## $p_low
## [1] 0.7631492
## 
## $z_high
## [1] -3.411625
## 
## $p_high
## [1] 0.0003228847
## 
## $ci_level
## [1] 0.9
## 
## $ci_lower
## [1] -0.2752272
## 
## $ci_upper
## [1] -0.03110717
## 
## $equivalent
## [1] FALSE
```

``` r
# equivalence test for q
tost_z(est=mod2_GEI_boot_results["q","Estimate"],
       se=mod2_GEI_boot_results["q","SE"],
       low=-.10,
       high=.10, alpha = 0.05)
```

```
## $estimate
## [1] -0.1761835
## 
## $se
## [1] 0.08068671
## 
## $low_bound
## [1] -0.1
## 
## $high_bound
## [1] 0.1
## 
## $alpha
## [1] 0.05
## 
## $z_low
## [1] -0.9441892
## 
## $p_low
## [1] 0.8274635
## 
## $z_high
## [1] -3.422912
## 
## $p_high
## [1] 0.0003097706
## 
## $ci_level
## [1] 0.9
## 
## $ci_lower
## [1] -0.3089014
## 
## $ci_upper
## [1] -0.04346569
## 
## $equivalent
## [1] FALSE
```



### Figure 


``` r
# plot the results

# start with obtaining predicted values for means and differences

# refit reduced and full models with GGGI in original scale


mod2_GEI_unstd<-lmer(FM.z~gndr.c+gei.cm+gndr.c:gei.cm+
                 (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

mod2_GEI_unstd_red<-lmer(FM.z~gndr.c+
                       (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
                     control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))


p<-
  emmip(
    mod2_GEI_unstd, 
    gndr.c ~ gei.cm,
    at=list(gndr.c = c(-0.5,0.5),
            gei.cm=
              seq(from=round(range(GII_country_means$gei.cm,na.rm=T)[1],2),
                  to=round(range(GII_country_means$gei.cm,na.rm=T)[2],2),
                  by=0.001)),
    plotit=F,CIs=T,lmerTest.limit = 1e6,disable.pbkrtest=T)

p$gndr.c<-p$tvar
levels(p$gndr.c)<-c("Women","Men")

# obtain min and max for aligned plots
min.y.comp<-min(p$LCL)
max.y.comp<-max(p$UCL)

# Men and Women mean distributions

p3<-coefficients(mod2_GEI_unstd_red)$cntry
p3<-cbind(rbind(p3,p3),weight=rep(c(-0.5,0.5),each=nrow(p3)))
p3$xvar<-p3$`(Intercept)`+p3$gndr.c*p3$weight
p3$gndr.c<-as.factor(p3$weight)
levels(p3$gndr.c)<-c("Women","Men")

# obtain min and max for aligned plots
min.y.mean.distr<-min(p3$xvar)
max.y.mean.distr<-max(p3$xvar)


# obtain the coefs for the gndr.c-effect (difference) as function of gei.cm

p2<-data.frame(
  emtrends(mod2_GEI_unstd,var="gndr.c",
           specs="gei.cm",
           at=list(#gndr.c = c(-0.5,0.5),
             gei.cm=
               seq(from=round(range(GII_country_means$gei.cm,na.rm=T)[1],2),
                   to=round(range(GII_country_means$gei.cm,na.rm=T)[2],2),
                   by=0.001)),
           lmerTest.limit = 1e6,disable.pbkrtest=T))

p2$yvar<-p2$gndr.c.trend
p2$xvar<-p2$gei.cm
p2$LCL<-p2$lower.CL
p2$UCL<-p2$upper.CL

# obtain min and max for aligned plots
min.y.diff<-min(p2$LCL)
max.y.diff<-max(p2$UCL)

# difference score distribution

p4<-coefficients(mod2_GEI_unstd_red)$cntry
p4$xvar=(+1)*p4$gndr.c

# obtain mix and max for aligned plots

min.y.diff.distr<-min(p4$xvar)
max.y.diff.distr<-max(p4$xvar)

# define mins and maxs

min.y.pred<-
  ifelse(min.y.comp<min.y.mean.distr,min.y.comp,min.y.mean.distr)

max.y.pred<-
  ifelse(max.y.comp>max.y.mean.distr,max.y.comp,max.y.mean.distr)

min.y.narrow<-
  ifelse(min.y.diff<min.y.diff.distr,min.y.diff,min.y.diff.distr)

max.y.narrow<-
  ifelse(max.y.diff>max.y.diff.distr,max.y.diff,max.y.diff.distr)

# Figures 

# p1

# scaled simple effects to the plot
pvals<-round_tidy(ddsc_mod2_GEI$results[6:7,"p.value"],3)

ests<-
  round_tidy(ddsc_mod2_GEI$results[6:7,"estimate"],2)

coef1<-paste0("std. b11 = ",ests[1],", p = ",pvals[1])
coef2<-paste0("std. b21 = ",ests[2],", p = ",pvals[2])
coefs<-data.frame(gndr.c=c("Women","Men"),
                  coef=c(coef1,coef2))

coef_q<-paste0("Cohen's q = ",round_tidy(ddsc_mod2_GEI$results["q_b11_b21","estimate"],2),", p = ",
               substr(round_tidy(ddsc_mod2_GEI$results["interaction","p.value"],3),2,5))  
# prediction plot for difference score

pvals2<-round_tidy(ddsc_mod2_GEI$results["interaction","p.value"],3)

ests2<-
  round_tidy(-1*ddsc_mod2_GEI$results["r_xy1y2","estimate"],2)

coefs2<-paste0("difference score correlation = ",ests2,", p = ",pvals2)
# attempt to make a flag-plot

flag_points_GEI<-coefficients(mod2_GEI_unstd_red)$cntry

flag_points_GEI$women<-flag_points_GEI$`(Intercept)`+(-0.5)*flag_points_GEI$gndr.c
flag_points_GEI$men<-flag_points_GEI$`(Intercept)`+(0.5)*flag_points_GEI$gndr.c
flag_points_GEI$mean_level<-flag_points_GEI$`(Intercept)`
flag_points_GEI$difference<-flag_points_GEI$men-flag_points_GEI$women
flag_points_GEI$cntry<-rownames(flag_points_GEI)

flag_points_GEI<-
  left_join(x=flag_points_GEI,
            y=GII_country_means[,c("ISO2","gei.cm")],by=c("cntry"="ISO2"))
#flag_points_GEI

flag_points_GEI_long<-
  data.frame(mean_level=c(flag_points_GEI$women,
                          flag_points_GEI$men),
             gei.cm=c(flag_points_GEI$gei.cm,
                      flag_points_GEI$gei.cm),
             cntry=c(flag_points_GEI$cntry,
                     flag_points_GEI$cntry),
             gndr.c=rep(c("Women","Men"),each=nrow(flag_points_GEI)
             ))


p1.FM.flags<-
  ggplot(p,aes(y=yvar,x=gei.cm,color=gndr.c))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("Gender Equality Index (Average 2002-2023)")+
  #ylim(c(min.y.pred,max.y.pred))+
  #xlim(c(0.60,1.00))+
  ylab("Value male-typicality (Average 2002-2023)")+
  scale_color_manual(values=met.brewer("Archambault")[c(6,2)])+
  theme(legend.position = "top",
        legend.title=element_blank(),
        text=element_text(size=16,  family="sans"),
        panel.background = element_rect(fill = "white",
                                        #colour = "black",
                                        #size = 0.5, linetype = "solid"
        ),
        panel.grid.major.x = element_line(linewidth = 0.5, linetype = 2,
                                          colour = "gray"))+
  geom_text(data = coefs,show.legend=F,
            aes(label=coef,x=0.65,
                y=c(-0.35,-0.40),size=14,hjust="left"))+
  geom_text(inherit.aes=F,aes(x=0.65,y=-0.45,
                              label=coef_q,size=14,hjust="left"),
            show.legend=F)+
  geom_point(data=flag_points_GEI_long,size=8,alpha=0.30,
             aes(x=gei.cm,y=mean_level))+
  geom_line(data=flag_points_GEI_long,aes(group = cntry,x=gei.cm,y=mean_level),
            color="black",linetype=2)+
  geom_flag(data=flag_points_GEI_long,show.legend=F,
            aes(country=tolower(cntry),size=14,x=gei.cm,y=mean_level))

p2.FM.flags<-ggplot(p2,aes(y=yvar,x=gei.cm))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("Gender Equality Index (Average 2002-2023)")+
  #xlim(c(0.60,1.00))+
  #ylim(c(min.y.narrow,max.y.narrow))+
  ylab("Difference in value male-typicality")+
  #scale_color_manual(values=met.brewer("Archambault")[c(6,2)])+
  theme(legend.position = "right",
        legend.title=element_blank(),
        text=element_text(size=16,  family="sans"),
        panel.background = element_rect(fill = "white",
                                        #colour = "black",
                                        #size = 0.5, linetype = "solid"
        ),
        panel.grid.major.x = element_line(linewidth = 0.5, linetype = 2,
                                          colour = "gray"))+
  #geom_text(coef2,aes(x=0.63,y=min(p2$LCL)))
  geom_text(data = data.frame(coefs2),show.legend=F,
            aes(label=coefs2,x=0.70,
                y=c(round(min(p2$LCL),2))+0.05,size=18,hjust="left"))+
  geom_flag(data=flag_points_GEI,show.legend=F,
            aes(country=tolower(cntry),size=18,x=gei.cm,y=difference))

pflag_comb<-
  ggarrange(p1.FM.flags,p2.FM.flags,align = "v",
            ncol=1,nrow=2,heights=c(2,1))
```

```
## Warning in geom_text(inherit.aes = F, aes(x = 0.65, y = -0.45, label = coef_q, : All aesthetics have length 1, but the data has 662 rows.
## ℹ Please consider using `annotate()` or provide this layer with data containing a single row.
```

```
## Warning: Removed 2 rows containing missing values or values outside the scale range (`geom_point()`).
```

```
## Warning: Removed 2 rows containing missing values or values outside the scale range (`geom_line()`).
```

```
## Warning: Removed 2 rows containing missing values or values outside the scale range (`geom_flag()`).
```

```
## Warning: Removed 1 row containing missing values or values outside the scale range (`geom_flag()`).
```

``` r
pflag_comb
```

![](Analysis_files/figure-html/unnamed-chunk-28-1.png)<!-- -->

``` r
png(filename = 
      "../results/GEI_flags.png",
    units = "cm",
    width = 21.0,height=29.7*(3/4),res = 300)
pflag_comb
dev.off()
```

```
## png 
##   2
```

## mod2 with Gender-equality index (GGGI)


``` r
mod2_GGGI<-lmer(FM.z~gndr.c+gggi.z.cm+gndr.c:gggi.z.cm+
                 (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2_GGGI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: FM.z ~ gndr.c + gggi.z.cm + gndr.c:gggi.z.cm + (gndr.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##  964693.5  964779.0 -482338.8  964677.5    323844 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.3637 -0.6143  0.0047  0.6004  8.2367 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr  
##  cntry    (Intercept) 0.062856 0.25071        
##           gndr.c      0.006122 0.07824  -0.04 
##  Residual             1.008538 1.00426        
## Number of obs: 323852, groups:  cntry, 34
## 
## Fixed effects:
##                  Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)       0.07381    0.04305 33.92892   1.714  0.09557 .  
## gndr.c            0.41172    0.01405 32.68018  29.294  < 2e-16 ***
## gggi.z.cm        -0.13933    0.04372 34.00190  -3.187  0.00308 ** 
## gndr.c:gggi.z.cm  0.04679    0.01454 35.04682   3.219  0.00277 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c ggg.z.
## gndr.c      -0.034              
## gggi.z.cm    0.000  0.000       
## gndr.c:gg..  0.000 -0.009 -0.034
```

``` r
getFE(mod2_GGGI,round=3)
```

```
##                    Est.    SE     df      t     p     LL     UL
## (Intercept)       0.074 0.043 33.929  1.714 0.096 -0.014  0.161
## gndr.c            0.412 0.014 32.680 29.294 0.000  0.383  0.440
## gggi.z.cm        -0.139 0.044 34.002 -3.187 0.003 -0.228 -0.050
## gndr.c:gggi.z.cm  0.047 0.015 35.047  3.219 0.003  0.017  0.076
```

``` r
getVC(mod2_GGGI)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.25 0.06
## 2    cntry      gndr.c   <NA>  0.08 0.01
## 3    cntry (Intercept) gndr.c -0.04 0.00
## 4 Residual        <NA>   <NA>  1.00 1.01
```

``` r
r2mlm(mod2_GGGI,bargraph = F)
```

```
## $Decompositions
##                      total
## fixed           0.04977292
## slope variation 0.00134440
## mean variation  0.05573267
## sigma2          0.89315002
## 
## $R2s
##          total
## f   0.04977292
## v   0.00134440
## m   0.05573267
## fv  0.05111732
## fvm 0.10684998
```


### Deconstructed associations


``` r
t1<-Sys.time()
ddsc_mod2_GGGI<-
  ddsc_ml(model = mod2_GGGI,
          predictor = "gggi.z.cm",
          moderator = "gndr.c",moderator_values = c(-0.5,0.5),
          re_cov_test = T)
t2<-Sys.time()
t2-t1
```

```
## Time difference of 47.63345 secs
```

``` r
round(ddsc_mod2_GGGI$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.068        0.262        1.016     1.085 0.063   6999.706 0.998   0.998
## 2        0.5         0.062        0.248        1.016     1.078 0.057   5875.088 0.997   0.997
```

``` r
round(ddsc_mod2_GGGI$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gggi.z.cm
## means_y1           0.251 0.284    1.000           1.000    0.951           0.951    -0.420
## means_y1_scaled    0.845 0.954    1.000           1.000    0.951           0.951    -0.420
## means_y2          -0.151 0.310    0.951           0.951    1.000           1.000    -0.525
## means_y2_scaled   -0.506 1.044    0.951           0.951    1.000           1.000    -0.525
## gggi.z.cm          0.000 1.000   -0.420          -0.420   -0.525          -0.525     1.000
## gggi.z.cm_scaled   0.000 1.000   -0.420          -0.420   -0.525          -0.525     1.000
## diff_score         0.402 0.096   -0.120          -0.120   -0.420          -0.420     0.457
## diff_score_scaled  1.351 0.324   -0.120          -0.120   -0.420          -0.420     0.457
##                   gggi.z.cm_scaled diff_score diff_score_scaled
## means_y1                    -0.420     -0.120            -0.120
## means_y1_scaled             -0.420     -0.120            -0.120
## means_y2                    -0.525     -0.420            -0.420
## means_y2_scaled             -0.525     -0.420            -0.420
## gggi.z.cm                    1.000      0.457             0.457
## gggi.z.cm_scaled             1.000      0.457             0.457
## diff_score                   0.457      1.000             1.000
## diff_score_scaled            0.457      1.000             1.000
```

``` r
round(ddsc_mod2_GGGI$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                      -0.486 0.151 35.047  -3.219   0.003   -0.792   -0.179
## w_11                         -0.163 0.045 34.023  -3.651   0.001   -0.253   -0.072
## w_21                         -0.116 0.044 33.995  -2.630   0.013   -0.206   -0.026
## r_xy1                        -0.573 0.157 34.023  -3.651   0.001   -0.892   -0.254
## r_xy2                        -0.373 0.142 33.995  -2.630   0.013   -0.662   -0.085
## b_11                         -0.548 0.150 34.023  -3.651   0.001   -0.852   -0.243
## b_21                         -0.390 0.148 33.995  -2.630   0.013   -0.692   -0.089
## main_effect                  -0.139 0.044 34.002  -3.187   0.003   -0.228   -0.050
## moderator_effect              0.412 0.014 32.680  29.294   0.000    0.383    0.440
## interaction                   0.047 0.015 35.047   3.219   0.003    0.017    0.076
## q_b11_b21                    -0.203    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                  -0.260    NA     NA      NA      NA       NA       NA
## cross_over_point             -8.800    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.093 0.046 34.013  -2.029   0.050   -0.185    0.000
## interaction_vs_main_bscale   -0.311 0.153 34.013  -2.029   0.050   -0.623    0.000
## interaction_vs_main_rscale   -0.273 0.141 34.019  -1.940   0.061   -0.560    0.013
## dadas                        -0.232 0.088 33.995  -2.630   0.994   -0.411   -0.053
## dadas_bscale                 -0.780 0.297 33.995  -2.630   0.994   -1.383   -0.177
## dadas_rscale                 -0.747 0.284 33.995  -2.630   0.994   -1.324   -0.170
## abs_diff                      0.047 0.015 35.047   3.219   0.001    0.017    0.076
## abs_sum                       0.279 0.087 34.002   3.187   0.002    0.101    0.456
## abs_diff_bscale               0.157 0.049 35.047   3.219   0.001    0.058    0.257
## abs_sum_bscale                0.938 0.294 34.002   3.187   0.002    0.340    1.536
## abs_diff_rscale               0.200 0.051 34.851   3.904   0.000    0.096    0.304
## abs_sum_rscale                0.947 0.295 34.003   3.210   0.001    0.347    1.546
```

``` r
round(ddsc_mod2_GGGI$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.003 -0.161  0.824  1.000  0.364
```

``` r
# Structural path model
d_GGGI<-ddsc_mod2_GGGI$ddsc_sem_fit$data

ddsc_sem_GGGI<-
  ddsc_sem(data=d_GGGI,x = "gggi.z.cm",y1="means_y2",
           y2="means_y1",q_sesoi = .10)
round(ddsc_sem_GGGI$results,3)
```

```
##                                    est    se       z pvalue ci.lower ci.upper
## r_xy1_y2                        -0.457 0.153  -2.997  0.003   -0.756   -0.158
## r_xy1                           -0.525 0.146  -3.601  0.000   -0.811   -0.239
## r_xy2                           -0.420 0.156  -2.695  0.007   -0.725   -0.114
## b_11                            -0.548 0.152  -3.601  0.000   -0.847   -0.250
## b_21                            -0.400 0.149  -2.695  0.007   -0.692   -0.109
## b_10                            -0.506 0.150  -3.373  0.001   -0.800   -0.212
## b_20                             0.845 0.146   5.774  0.000    0.558    1.132
## res_cov_y1_y2                    0.707 0.176   4.008  0.000    0.361    1.052
## diff_b10_b20                    -1.351 0.049 -27.762  0.000   -1.446   -1.256
## diff_b11_b21                    -0.148 0.049  -2.997  0.003   -0.245   -0.051
## diff_rxy1_rxy2                  -0.106 0.050  -2.104  0.035   -0.204   -0.007
## q_b11_b21                       -0.192 0.076  -2.520  0.012   -0.341   -0.043
## q_rxy1_rxy2                     -0.137 0.065  -2.097  0.036   -0.264   -0.009
## cross_over_point                -9.126 3.063  -2.980  0.003  -15.128   -3.123
## sum_b11_b21                     -0.949 0.297  -3.197  0.001   -1.531   -0.367
## main_effect                     -0.474 0.148  -3.197  0.001   -0.765   -0.184
## interaction_vs_main_effect      -0.326 0.153  -2.137  0.033   -0.626   -0.027
## diff_abs_b11_abs_b21             0.148 0.049   2.997  0.003    0.051    0.245
## abs_diff_b11_b21                 0.148 0.049   2.997  0.001    0.051    0.245
## abs_sum_b11_b21                  0.949 0.297   3.197  0.001    0.367    1.531
## dadas                           -0.801 0.297  -2.695  0.996   -1.383   -0.219
## q_r_equivalence                  0.037 0.065   0.562  0.713       NA       NA
## q_b_equivalence                  0.092 0.076   1.208  0.886       NA       NA
## cross_over_point_equivalence     9.126 3.063   2.980  0.999       NA       NA
## cross_over_point_minimal_effect  9.126 3.063   2.980  0.001       NA       NA
```

``` r
round(ddsc_sem_GGGI$variance_test,3)
```

```
##             est    se      z pvalue ci.lower ci.upper
## cov_y1y2  0.920 0.229  4.019  0.000    0.471    1.368
## var_y1    1.057 0.256  4.123  0.000    0.555    1.560
## var_y2    0.884 0.214  4.123  0.000    0.464    1.304
## var_diff  0.174 0.110  1.572  0.116   -0.043    0.390
## var_ratio 1.196 0.126  9.464  0.000    0.949    1.444
## cor_y1y2  0.951 0.016 58.460  0.000    0.919    0.983
```

``` r
## random intercept mlm with double-entries for each country (men and women)

d_GGGI_long <- d_GGGI %>%
  # move row names into a column
  rownames_to_column("cntry") %>%
  # pivot only the means_y1 / means_y2 columns
  pivot_longer(
    cols = c(means_y1, means_y2),
    names_to = "y",
    values_to = "means"
  ) %>%
  # create a sgender column based on y1/y2
  mutate(
    gndr.c = case_when(
      y == "means_y1" ~ -0.5,
      y == "means_y2" ~ 0.5
    )
  ) %>%
  dplyr::select(cntry,gggi.z.cm,means,gndr.c)


ddsc_mod2_GGGI_ri<-
  ddsc_ml(data=data.frame(d_GGGI_long),predictor = "gggi.z.cm",
          moderator = "gndr.c",
        DV = "means",lvl2_unit = "cntry",
        moderator_values = c(-0.5,0.5))
round(ddsc_mod2_GGGI_ri$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                      -0.457 0.157 32.000  -2.908   0.007   -0.777   -0.137
## w_11                         -0.163 0.046 33.772  -3.537   0.001   -0.257   -0.069
## w_21                         -0.119 0.046 33.772  -2.582   0.014   -0.213   -0.025
## r_xy1                        -0.525 0.149 33.772  -3.537   0.001   -0.827   -0.223
## r_xy2                        -0.420 0.163 33.772  -2.582   0.014   -0.750   -0.089
## b_11                         -0.549 0.155 33.772  -3.537   0.001   -0.865   -0.233
## b_21                         -0.401 0.155 33.772  -2.582   0.014   -0.716   -0.085
## main_effect                  -0.141 0.046 32.000  -3.102   0.004   -0.234   -0.048
## moderator_effect              0.402 0.015 32.000  26.933   0.000    0.371    0.432
## interaction                   0.044 0.015 32.000   2.908   0.007    0.013    0.075
## q_b11_b21                    -0.192    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                  -0.137    NA     NA      NA      NA       NA       NA
## cross_over_point             -9.126    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.097 0.048 39.006  -2.024   0.050   -0.194    0.000
## interaction_vs_main_bscale   -0.327 0.161 39.006  -2.024   0.050   -0.653    0.000
## interaction_vs_main_rscale   -0.367 0.175 38.185  -2.094   0.043   -0.721   -0.012
## dadas                        -0.238 0.092 33.772  -2.582   0.993   -0.426   -0.051
## dadas_bscale                 -0.802 0.310 33.772  -2.582   0.993   -1.433   -0.171
## dadas_rscale                 -0.839 0.325 33.772  -2.582   0.993   -1.500   -0.179
## abs_diff                      0.044 0.015 32.000   2.908   0.003    0.013    0.075
## abs_sum                       0.282 0.091 32.000   3.102   0.002    0.097    0.468
## abs_diff_bscale               0.148 0.051 32.000   2.908   0.003    0.044    0.252
## abs_sum_bscale                0.950 0.306 32.000   3.102   0.002    0.326    1.574
## abs_diff_rscale               0.106 0.053 36.615   2.001   0.026   -0.001    0.213
## abs_sum_rscale                0.945 0.307 32.004   3.080   0.002    0.320    1.570
```

``` r
# country-time multilevel model


mod2_GGGI_cntry_year<-
  lmer(FM.z.wt~gndr.c+gggi.z.cm+gndr.c:gggi.z.cm+
                 (gndr.c|cntry),data=diff_dat_cntry_year,REML=F,
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
```

```
## boundary (singular) fit: see help('isSingular')
```

``` r
summary(mod2_GGGI_cntry_year)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: FM.z.wt ~ gndr.c + gggi.z.cm + gndr.c:gggi.z.cm + (gndr.c | cntry)
##    Data: diff_dat_cntry_year
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##    -329.9    -298.0     173.0    -345.9       392 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -2.9614 -0.5796 -0.0102  0.5958  3.5182 
## 
## Random effects:
##  Groups   Name        Variance  Std.Dev. Corr  
##  cntry    (Intercept) 5.820e-02 0.241252       
##           gndr.c      5.755e-06 0.002399 -1.00 
##  Residual             1.823e-02 0.135034       
## Number of obs: 400, groups:  cntry, 34
## 
## Fixed effects:
##                   Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)        0.06553    0.04211  33.42891   1.556  0.12907    
## gndr.c             0.41547    0.01365 358.04154  30.432  < 2e-16 ***
## gggi.z.cm         -0.13614    0.04300  34.18432  -3.166  0.00325 ** 
## gndr.c:gggi.z.cm   0.03652    0.01467 359.46782   2.490  0.01321 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c ggg.z.
## gndr.c      -0.030              
## gggi.z.cm   -0.010  0.000       
## gndr.c:gg..  0.000 -0.144 -0.028
## optimizer (bobyqa) convergence code: 0 (OK)
## boundary (singular) fit: see help('isSingular')
```

``` r
getFE(mod2_GGGI_cntry_year,round=3)
```

```
##                    Est.    SE      df      t     p     LL     UL
## (Intercept)       0.066 0.042  33.429  1.556 0.129 -0.020  0.151
## gndr.c            0.415 0.014 358.042 30.432 0.000  0.389  0.442
## gggi.z.cm        -0.136 0.043  34.184 -3.166 0.003 -0.224 -0.049
## gndr.c:gggi.z.cm  0.037 0.015 359.468  2.490 0.013  0.008  0.065
```

``` r
getVC(mod2_GGGI_cntry_year)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.24 0.06
## 2    cntry      gndr.c   <NA>  0.00 0.00
## 3    cntry (Intercept) gndr.c -1.00 0.00
## 4 Residual        <NA>   <NA>  0.14 0.02
```

``` r
r2mlm(mod2_GGGI,bargraph = F)
```

```
## $Decompositions
##                      total
## fixed           0.04977292
## slope variation 0.00134440
## mean variation  0.05573267
## sigma2          0.89315002
## 
## $R2s
##          total
## f   0.04977292
## v   0.00134440
## m   0.05573267
## fv  0.05111732
## fvm 0.10684998
```

``` r
ddsc_mod2_GGGI_cntry_year<-
  ddsc_ml(model = mod2_GGGI_cntry_year,
          predictor = "gggi.z.cm",
          moderator = "gndr.c",moderator_values = c(-0.5,0.5),
          re_cov_test = T)

round(ddsc_mod2_GGGI_cntry_year$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.059        0.243        0.023     0.082 0.720      8.029 0.996   0.954
## 2        0.5         0.055        0.234        0.023     0.077 0.704      8.029 0.996   0.950
```

``` r
round(ddsc_mod2_GGGI_cntry_year$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gggi.z.cm
## means_y1           0.272 0.278    1.000           1.000    0.948           0.948    -0.411
## means_y1_scaled    0.934 0.955    1.000           1.000    0.948           0.948    -0.411
## means_y2          -0.140 0.304    0.948           0.948    1.000           1.000    -0.535
## means_y2_scaled   -0.479 1.043    0.948           0.948    1.000           1.000    -0.535
## gggi.z.cm          0.000 1.000   -0.411          -0.411   -0.535          -0.535     1.000
## gggi.z.cm_scaled   0.000 1.000   -0.411          -0.411   -0.535          -0.535     1.000
## diff_score         0.411 0.097   -0.104          -0.104   -0.413          -0.413     0.499
## diff_score_scaled  1.413 0.332   -0.104          -0.104   -0.413          -0.413     0.499
##                   gggi.z.cm_scaled diff_score diff_score_scaled
## means_y1                    -0.411     -0.104            -0.104
## means_y1_scaled             -0.411     -0.104            -0.104
## means_y2                    -0.535     -0.413            -0.413
## means_y2_scaled             -0.535     -0.413            -0.413
## gggi.z.cm                    1.000      0.499             0.499
## gggi.z.cm_scaled             1.000      0.499             0.499
## diff_score                   0.499      1.000             1.000
## diff_score_scaled            0.499      1.000             1.000
```

``` r
round(ddsc_mod2_GGGI_cntry_year$results,3)
```

```
##                            estimate    SE      df t.ratio p.value ci.lower ci.upper
## r_xy1y2                      -0.377 0.152 359.468  -2.490   0.013   -0.676   -0.079
## w_11                         -0.154 0.044  34.183  -3.523   0.001   -0.243   -0.065
## w_21                         -0.118 0.043  34.188  -2.715   0.010   -0.206   -0.030
## r_xy1                        -0.555 0.158  34.183  -3.523   0.001   -0.876   -0.235
## r_xy2                        -0.388 0.143  34.188  -2.715   0.010   -0.679   -0.098
## b_11                         -0.531 0.151  34.183  -3.523   0.001   -0.837   -0.225
## b_21                         -0.405 0.149  34.188  -2.715   0.010   -0.709   -0.102
## main_effect                  -0.136 0.043  34.184  -3.166   0.003   -0.224   -0.049
## moderator_effect              0.415 0.014 358.042  30.432   0.000    0.389    0.442
## interaction                   0.037 0.015 359.468   2.490   0.013    0.008    0.065
## q_b11_b21                    -0.161    NA      NA      NA      NA       NA       NA
## q_rxy1_rxy2                  -0.216    NA      NA      NA      NA       NA       NA
## cross_over_point            -11.376    NA      NA      NA      NA       NA       NA
## interaction_vs_main          -0.100 0.045  34.460  -2.212   0.034   -0.191   -0.008
## interaction_vs_main_bscale   -0.342 0.155  34.460  -2.212   0.034   -0.657   -0.028
## interaction_vs_main_rscale   -0.305 0.142  34.552  -2.137   0.040   -0.594   -0.015
## dadas                        -0.236 0.087  34.188  -2.715   0.995   -0.412   -0.059
## dadas_bscale                 -0.811 0.299  34.188  -2.715   0.995   -1.417   -0.204
## dadas_rscale                 -0.776 0.286  34.188  -2.715   0.995   -1.357   -0.195
## abs_diff                      0.037 0.015 359.468   2.490   0.007    0.008    0.065
## abs_sum                       0.272 0.086  34.184   3.166   0.002    0.098    0.447
## abs_diff_bscale               0.126 0.050 359.468   2.490   0.007    0.026    0.225
## abs_sum_bscale                0.936 0.296  34.184   3.166   0.002    0.335    1.537
## abs_diff_rscale               0.167 0.053 149.329   3.182   0.001    0.063    0.271
## abs_sum_rscale                0.943 0.296  34.184   3.184   0.002    0.341    1.546
```

``` r
round(ddsc_mod2_GGGI_cntry_year$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.002 -0.306  0.377  1.000  0.539
```

### Comparisons across models


``` r
# full multilevel model (individuals as entries)
round(ddsc_mod2_GGGI$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df t.ratio p.value ci.lower ci.upper
## r_xy1             -0.5733 0.1570 34.0227 -3.6513  0.0009  -0.8923  -0.2542
## r_xy2             -0.3734 0.1420 33.9945 -2.6302  0.0127  -0.6619  -0.0849
## b_11              -0.5476 0.1500 34.0227 -3.6513  0.0009  -0.8523  -0.2428
## b_21              -0.3901 0.1483 33.9945 -2.6302  0.0127  -0.6916  -0.0887
## main_effect       -0.1393 0.0437 34.0019 -3.1866  0.0031  -0.2282  -0.0505
## moderator_effect   0.4117 0.0141 32.6802 29.2942  0.0000   0.3831   0.4403
## interaction        0.0468 0.0145 35.0468  3.2186  0.0028   0.0173   0.0763
## q_b11_b21         -0.2029     NA      NA      NA      NA       NA       NA
```

``` r
# country-level structural path model (country-gender means as entries)
round(ddsc_sem_GGGI$results[c("r_xy1","r_xy2","b_11","b_21","q_b11_b21","diff_b11_b21"),],4)
```

```
##                  est     se       z pvalue ci.lower ci.upper
## r_xy1        -0.5254 0.1459 -3.6010 0.0003  -0.8114  -0.2395
## r_xy2        -0.4196 0.1557 -2.6954 0.0070  -0.7247  -0.1145
## b_11         -0.5484 0.1523 -3.6010 0.0003  -0.8470  -0.2499
## b_21         -0.4004 0.1485 -2.6954 0.0070  -0.6915  -0.1093
## q_b11_b21    -0.1920 0.0762 -2.5202 0.0117  -0.3414  -0.0427
## diff_b11_b21 -0.1480 0.0494 -2.9971 0.0027  -0.2449  -0.0512
```

``` r
# multilevel model (country-gender means as entries)
round(ddsc_mod2_GGGI_ri$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df t.ratio p.value ci.lower ci.upper
## r_xy1             -0.5254 0.1486 33.7716 -3.5369  0.0012  -0.8274  -0.2235
## r_xy2             -0.4196 0.1625 33.7716 -2.5821  0.0143  -0.7499  -0.0893
## b_11              -0.5490 0.1552 33.7716 -3.5369  0.0012  -0.8645  -0.2335
## b_21              -0.4008 0.1552 33.7716 -2.5821  0.0143  -0.7163  -0.0853
## main_effect       -0.1411 0.0455 32.0000 -3.1016  0.0040  -0.2338  -0.0484
## moderator_effect   0.4019 0.0149 32.0000 26.9328  0.0000   0.3715   0.4323
## interaction        0.0440 0.0151 32.0000  2.9076  0.0066   0.0132   0.0749
## q_b11_b21         -0.1923     NA      NA      NA      NA       NA       NA
```

``` r
# multilevel model (country-year-gender means as entries)
round(ddsc_mod2_GGGI_cntry_year$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE       df t.ratio p.value ci.lower ci.upper
## r_xy1             -0.5553 0.1576  34.1826 -3.5229  0.0012  -0.8756  -0.2350
## r_xy2             -0.3881 0.1430  34.1882 -2.7149  0.0103  -0.6786  -0.0977
## b_11              -0.5308 0.1507  34.1826 -3.5229  0.0012  -0.8370  -0.2247
## b_21              -0.4053 0.1493  34.1882 -2.7149  0.0103  -0.7086  -0.1020
## main_effect       -0.1361 0.0430  34.1843 -3.1658  0.0032  -0.2235  -0.0488
## moderator_effect   0.4155 0.0137 358.0415 30.4315  0.0000   0.3886   0.4423
## interaction        0.0365 0.0147 359.4678  2.4904  0.0132   0.0077   0.0654
## q_b11_b21         -0.1614     NA       NA      NA      NA       NA       NA
```

### Bootstrap and equivalence test

Takes a lot of time


``` r
t1<-Sys.time()
mod2_GGGI_booted_fixef <-
  lme4::bootMer(
    x = mod2_GGGI,
    FUN = lme4::fixef,
    nsim = 1000,
    use.u = FALSE,
    seed = 12345,
    type = c("parametric"),
    verbose = FALSE
  )
t2<-Sys.time()
t2-t1
```

```
## Time difference of 1.446369 hours
```



``` r
# obtain all the bootstrap estimates
mod2_GGGI_boot_est <- data.frame(mod2_GGGI_booted_fixef$t)

# calculate estimates
mod2_GGGI_boot_est$w11<-mod2_GGGI_boot_est$gggi.z.cm+(-0.5)*mod2_GGGI_boot_est$gndr.c.gggi.z.cm
mod2_GGGI_boot_est$w21<-mod2_GGGI_boot_est$gggi.z.cm+(0.5)*mod2_GGGI_boot_est$gndr.c.gggi.z.cm
mod2_GGGI_boot_est$b11<-mod2_GGGI_boot_est$w11/ddsc_mod2_GGGI$SDs["SD_pooled"]
mod2_GGGI_boot_est$b21<-mod2_GGGI_boot_est$w21/ddsc_mod2_GGGI$SDs["SD_pooled"]
mod2_GGGI_boot_est$r_xy1<-mod2_GGGI_boot_est$w11/ddsc_mod2_GGGI$SDs["SD_y1"]
mod2_GGGI_boot_est$r_xy2<-mod2_GGGI_boot_est$w21/ddsc_mod2_GGGI$SDs["SD_y2"]
mod2_GGGI_boot_est$q_b<-atanh(mod2_GGGI_boot_est$b11)-atanh(mod2_GGGI_boot_est$b21)
mod2_GGGI_boot_est$q<-atanh(mod2_GGGI_boot_est$r_xy1)-atanh(mod2_GGGI_boot_est$r_xy2)
```

```
## Warning in atanh(mod2_GGGI_boot_est$r_xy1): NaNs produced
```

``` r
# Calculate bootstrap summary statistics
mod2_GGGI_boot_results <- t(as.data.frame(sapply(
  mod2_GGGI_boot_est,
  function(x) {
    c(
      Estimate = mean(x, na.rm = TRUE),
      SE = stats::sd(x, na.rm = TRUE),
      stats::quantile(x, c((1 - .95) / 2,
                           1 - (1 - .95) / 2), na.rm = TRUE)
    )
  }
)))

mod2_GGGI_boot_results
```

```
##                     Estimate         SE         2.5%       97.5%
## X.Intercept.      0.07270732 0.04345944 -0.007450383  0.15892208
## gndr.c            0.41229727 0.01495826  0.381523603  0.43986940
## gggi.z.cm        -0.13790324 0.04230640 -0.227073265 -0.05762603
## gndr.c.gggi.z.cm  0.04688308 0.01485422  0.018266919  0.07796229
## w11              -0.16134478 0.04332786 -0.255024430 -0.07819334
## w21              -0.11446170 0.04257562 -0.201748804 -0.03071462
## b11              -0.54293920 0.14580200 -0.858179352 -0.26312737
## b21              -0.38517357 0.14327065 -0.678902245 -0.10335738
## r_xy1            -0.56841143 0.15264236 -0.898441200 -0.27547210
## r_xy2            -0.36865312 0.13712564 -0.649783495 -0.09892428
## q_b              -0.22059870 0.10461516 -0.465689308 -0.07736998
## q                -0.28628542 0.14054032 -0.642147209 -0.11746501
```

``` r
# equivalence test for q_b
tost_z(est=mod2_GGGI_boot_results["q_b","Estimate"],
       se=mod2_GGGI_boot_results["q_b","SE"],
       low=-.10,
       high=.10, alpha = 0.05)
```

```
## $estimate
## [1] -0.2205987
## 
## $se
## [1] 0.1046152
## 
## $low_bound
## [1] -0.1
## 
## $high_bound
## [1] 0.1
## 
## $alpha
## [1] 0.05
## 
## $z_low
## [1] -1.152784
## 
## $p_low
## [1] 0.8755005
## 
## $z_high
## [1] -3.064553
## 
## $p_high
## [1] 0.001089977
## 
## $ci_level
## [1] 0.9
## 
## $ci_lower
## [1] -0.3926753
## 
## $ci_upper
## [1] -0.04852208
## 
## $equivalent
## [1] FALSE
```

``` r
# equivalence test for q
tost_z(est=mod2_GGGI_boot_results["q","Estimate"],
       se=mod2_GGGI_boot_results["q","SE"],
       low=-.10,
       high=.10, alpha = 0.05)
```

```
## $estimate
## [1] -0.2862854
## 
## $se
## [1] 0.1405403
## 
## $low_bound
## [1] -0.1
## 
## $high_bound
## [1] 0.1
## 
## $alpha
## [1] 0.05
## 
## $z_low
## [1] -1.325495
## 
## $p_low
## [1] 0.9074964
## 
## $z_high
## [1] -2.748574
## 
## $p_high
## [1] 0.002992759
## 
## $ci_level
## [1] 0.9
## 
## $ci_lower
## [1] -0.5174537
## 
## $ci_upper
## [1] -0.05511717
## 
## $equivalent
## [1] FALSE
```



### Figure


``` r
# plot the results

# start with obtaining predicted values for means and differences

# refit reduced and full models with GGGI in original scale


mod2_GGGI_unstd<-lmer(FM.z~gndr.c+gggi.cm+gndr.c:gggi.cm+
                       (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
                     control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

mod2_GGGI_unstd_red<-lmer(FM.z~gndr.c+
                           (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
                         control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))


p<-
  emmip(
    mod2_GGGI_unstd, 
    gndr.c ~ gggi.cm,
    at=list(gndr.c = c(-0.5,0.5),
            gggi.cm=
              seq(from=round(range(GGGI_country_means$gggi.cm,na.rm=T)[1],2),
                  to=round(range(GGGI_country_means$gggi.cm,na.rm=T)[2],2),
                  by=0.001)),
    plotit=F,CIs=T,lmerTest.limit = 1e6,disable.pbkrtest=T)

p$gndr.c<-p$tvar
levels(p$gndr.c)<-c("Women","Men")

# obtain min and max for aligned plots
min.y.comp<-min(p$LCL)
max.y.comp<-max(p$UCL)

# Men and Women mean distributions

p3<-coefficients(mod2_GGGI_unstd_red)$cntry
p3<-cbind(rbind(p3,p3),weight=rep(c(-0.5,0.5),each=nrow(p3)))
p3$xvar<-p3$`(Intercept)`+p3$gndr.c*p3$weight
p3$gndr.c<-as.factor(p3$weight)
levels(p3$gndr.c)<-c("Women","Men")

# obtain min and max for aligned plots
min.y.mean.distr<-min(p3$xvar)
max.y.mean.distr<-max(p3$xvar)


# obtain the coefs for the gndr.c-effect (difference) as function of gggi.cm

p2<-data.frame(
  emtrends(mod2_GGGI_unstd,var="gndr.c",
           specs="gggi.cm",
           at=list(#gndr.c = c(-0.5,0.5),
             gggi.cm=
               seq(from=round(range(GGGI_country_means$gggi.cm,na.rm=T)[1],2),
                   to=round(range(GGGI_country_means$gggi.cm,na.rm=T)[2],2),
                   by=0.001)),
           lmerTest.limit = 1e6,disable.pbkrtest=T))

p2$yvar<-p2$gndr.c.trend
p2$xvar<-p2$gggi.cm
p2$LCL<-p2$lower.CL
p2$UCL<-p2$upper.CL

# obtain min and max for aligned plots
min.y.diff<-min(p2$LCL)
max.y.diff<-max(p2$UCL)

# difference score distribution

p4<-coefficients(mod2_GGGI_unstd_red)$cntry
p4$xvar=(+1)*p4$gndr.c

# obtain mix and max for aligned plots

min.y.diff.distr<-min(p4$xvar)
max.y.diff.distr<-max(p4$xvar)

# define mins and maxs

min.y.pred<-
  ifelse(min.y.comp<min.y.mean.distr,min.y.comp,min.y.mean.distr)

max.y.pred<-
  ifelse(max.y.comp>max.y.mean.distr,max.y.comp,max.y.mean.distr)

min.y.narrow<-
  ifelse(min.y.diff<min.y.diff.distr,min.y.diff,min.y.diff.distr)

max.y.narrow<-
  ifelse(max.y.diff>max.y.diff.distr,max.y.diff,max.y.diff.distr)

# Figures 

# p1

# scaled simple effects to the plot
pvals<-round_tidy(ddsc_mod2_GGGI$results[6:7,"p.value"],3)

ests<-
  round_tidy(ddsc_mod2_GGGI$results[6:7,"estimate"],2)

coef1<-paste0("std. b11 = ",ests[1],", p = ",pvals[1])
coef2<-paste0("std. b21 = ",ests[2],", p = ",pvals[2])
coefs<-data.frame(gndr.c=c("Women","Men"),
                  coef=c(coef1,coef2))

coef_q<-paste0("Cohen's q = ",round_tidy(ddsc_mod2_GGGI$results["q_b11_b21","estimate"],2),", p = ",
               substr(round_tidy(ddsc_mod2_GGGI$results["interaction","p.value"],3),2,5))  
# prediction plot for difference score

pvals2<-round_tidy(ddsc_mod2_GGGI$results["interaction","p.value"],3)

ests2<-
  round_tidy(-1*ddsc_mod2_GGGI$results["r_xy1y2","estimate"],2)

coefs2<-paste0("difference score correlation = ",ests2,", p = ",pvals2)
# attempt to make a flag-plot

flag_points_GGGI<-coefficients(mod2_GGGI_unstd_red)$cntry

flag_points_GGGI$women<-flag_points_GGGI$`(Intercept)`+(-0.5)*flag_points_GGGI$gndr.c
flag_points_GGGI$men<-flag_points_GGGI$`(Intercept)`+(0.5)*flag_points_GGGI$gndr.c
flag_points_GGGI$mean_level<-flag_points_GGGI$`(Intercept)`
flag_points_GGGI$difference<-flag_points_GGGI$men-flag_points_GGGI$women
flag_points_GGGI$cntry<-rownames(flag_points_GGGI)

flag_points_GGGI<-
  left_join(x=flag_points_GGGI,
            y=GGGI_country_means[,c("ISO2","gggi.cm")],by=c("cntry"="ISO2"))
#flag_points_GGGI

flag_points_GGGI_long<-
  data.frame(mean_level=c(flag_points_GGGI$women,
                          flag_points_GGGI$men),
             gggi.cm=c(flag_points_GGGI$gggi.cm,
                                 flag_points_GGGI$gggi.cm),
             cntry=c(flag_points_GGGI$cntry,
                     flag_points_GGGI$cntry),
             gndr.c=rep(c("Women","Men"),each=nrow(flag_points_GGGI)
             ))


p1.FM.flags<-
  ggplot(p,aes(y=yvar,x=gggi.cm,color=gndr.c))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("Global Gender Gap Index (Average 2002-2023)")+
  #ylim(c(min.y.pred,max.y.pred))+
  #xlim(c(0.60,1.00))+
  ylab("Value male-typicality (Average 2002-2023)")+
  scale_color_manual(values=met.brewer("Archambault")[c(6,2)])+
  theme(legend.position = "top",
        legend.title=element_blank(),
        text=element_text(size=16,  family="sans"),
        panel.background = element_rect(fill = "white",
                                        #colour = "black",
                                        #size = 0.5, linetype = "solid"
        ),
        panel.grid.major.x = element_line(linewidth = 0.5, linetype = 2,
                                          colour = "gray"))+
  geom_text(data = coefs,show.legend=F,
            aes(label=coef,x=0.61,
                y=c(-0.35,-0.40),size=14,hjust="left"))+
  geom_text(inherit.aes=F,aes(x=0.61,y=-0.45,
                              label=coef_q,size=14,hjust="left"),
            show.legend=F)+
  geom_point(data=flag_points_GGGI_long,size=8,alpha=0.30,
             aes(x=gggi.cm,y=mean_level))+
  geom_line(data=flag_points_GGGI_long,aes(group = cntry,x=gggi.cm,y=mean_level),
            color="black",linetype=2)+
  geom_flag(data=flag_points_GGGI_long,show.legend=F,
            aes(country=tolower(cntry),size=14,x=gggi.cm,y=mean_level))

p2.FM.flags<-ggplot(p2,aes(y=yvar,x=gggi.cm))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("Global Gender Gap Index (Average 2002-2023)")+
  #xlim(c(0.60,1.00))+
  #ylim(c(min.y.narrow,max.y.narrow))+
  ylab("Difference in value male-typicality")+
  #scale_color_manual(values=met.brewer("Archambault")[c(6,2)])+
  theme(legend.position = "right",
        legend.title=element_blank(),
        text=element_text(size=16,  family="sans"),
        panel.background = element_rect(fill = "white",
                                        #colour = "black",
                                        #size = 0.5, linetype = "solid"
        ),
        panel.grid.major.x = element_line(linewidth = 0.5, linetype = 2,
                                          colour = "gray"))+
  #geom_text(coef2,aes(x=0.63,y=min(p2$LCL)))
  geom_text(data = data.frame(coefs2),show.legend=F,
            aes(label=coefs2,x=0.70,
                y=c(round(min(p2$LCL),2))+0.05,size=18,hjust="left"))+
  geom_flag(data=flag_points_GGGI,show.legend=F,
            aes(country=tolower(cntry),size=18,x=gggi.cm,y=difference))

pflag_comb<-
  ggarrange(p1.FM.flags,p2.FM.flags,align = "v",
            ncol=1,nrow=2,heights=c(2,1))
```

```
## Warning in geom_text(inherit.aes = F, aes(x = 0.61, y = -0.45, label = coef_q, : All aesthetics have length 1, but the data has 502 rows.
## ℹ Please consider using `annotate()` or provide this layer with data containing a single row.
```

``` r
pflag_comb
```

![](Analysis_files/figure-html/unnamed-chunk-34-1.png)<!-- -->

``` r
png(filename = 
      "../results/GGGI_flags_new.png",
    units = "cm",
    width = 21.0,height=29.7*(3/4),res = 300)
pflag_comb
dev.off()
```

```
## png 
##   2
```


## mod2 with Gender-equality index (GDI)


``` r
mod2_GDI<-lmer(FM.z~gndr.c+gdi.z.cm+gndr.c:gdi.z.cm+
                 (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2_GDI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: FM.z ~ gndr.c + gdi.z.cm + gndr.c:gdi.z.cm + (gndr.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1306539.5 1306627.4 -653261.8 1306523.5    437735 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.3236 -0.6158  0.0017  0.6007  8.4520 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr  
##  cntry    (Intercept) 0.062450 0.24990        
##           gndr.c      0.006587 0.08116  -0.22 
##  Residual             1.016328 1.00813        
## Number of obs: 437743, groups:  cntry, 34
## 
## Fixed effects:
##                 Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)      0.05880    0.04290 33.95747   1.371   0.1794    
## gndr.c           0.42147    0.01436 34.26387  29.342   <2e-16 ***
## gdi.z.cm         0.02767    0.04356 34.02224   0.635   0.5296    
## gndr.c:gdi.z.cm  0.03119    0.01482 36.46302   2.105   0.0423 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c gd.z.c
## gndr.c      -0.208              
## gdi.z.cm     0.000  0.000       
## gndr.c:gd..  0.000 -0.008 -0.204
```

``` r
getFE(mod2_GDI,round=3)
```

```
##                  Est.    SE     df      t     p     LL    UL
## (Intercept)     0.059 0.043 33.957  1.371 0.179 -0.028 0.146
## gndr.c          0.421 0.014 34.264 29.342 0.000  0.392 0.451
## gdi.z.cm        0.028 0.044 34.022  0.635 0.530 -0.061 0.116
## gndr.c:gdi.z.cm 0.031 0.015 36.463  2.105 0.042  0.001 0.061
```

``` r
getVC(mod2_GDI)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.25 0.06
## 2    cntry      gndr.c   <NA>  0.08 0.01
## 3    cntry (Intercept) gndr.c -0.22 0.00
## 4 Residual        <NA>   <NA>  1.01 1.02
```

``` r
r2mlm(mod2_GDI,bargraph = F)
```

```
## $Decompositions
##                      total
## fixed           0.03935782
## slope variation 0.00145246
## mean variation  0.05585713
## sigma2          0.90333260
## 
## $R2s
##          total
## f   0.03935782
## v   0.00145246
## m   0.05585713
## fv  0.04081028
## fvm 0.09666740
```


### Deconstructed associations


``` r
t1<-Sys.time()
ddsc_mod2_GDI<-
  ddsc_ml(model = mod2_GDI,
          predictor = "gdi.z.cm",
          moderator = "gndr.c",moderator_values = c(-0.5,0.5),
          re_cov_test = T)
t2<-Sys.time()
t2-t1
```

```
## Time difference of 2.124499 mins
```

``` r
round(ddsc_mod2_GDI$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.068        0.262        1.016     1.085 0.063   6999.706 0.998   0.998
## 2        0.5         0.062        0.248        1.016     1.078 0.057   5875.088 0.997   0.997
```

``` r
round(ddsc_mod2_GDI$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gdi.z.cm gdi.z.cm_scaled
## means_y1           0.239 0.257    1.000           1.000    0.942           0.942    0.166           0.166
## means_y1_scaled    0.909 0.975    1.000           1.000    0.942           0.942    0.166           0.166
## means_y2          -0.173 0.270    0.942           0.942    1.000           1.000   -0.008          -0.008
## means_y2_scaled   -0.656 1.025    0.942           0.942    1.000           1.000   -0.008          -0.008
## gdi.z.cm           0.000 1.000    0.166           0.166   -0.008          -0.008    1.000           1.000
## gdi.z.cm_scaled    0.000 1.000    0.166           0.166   -0.008          -0.008    1.000           1.000
## diff_score         0.412 0.090    0.028           0.028   -0.309          -0.309    0.495           0.495
## diff_score_scaled  1.565 0.343    0.028           0.028   -0.309          -0.309    0.495           0.495
##                   diff_score diff_score_scaled
## means_y1               0.028             0.028
## means_y1_scaled        0.028             0.028
## means_y2              -0.309            -0.309
## means_y2_scaled       -0.309            -0.309
## gdi.z.cm               0.495             0.495
## gdi.z.cm_scaled        0.495             0.495
## diff_score             1.000             1.000
## diff_score_scaled      1.000             1.000
```

``` r
round(ddsc_mod2_GDI$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                      -0.345 0.164 36.463  -2.105   0.042   -0.677   -0.013
## w_11                          0.012 0.046 34.031   0.265   0.793   -0.081    0.105
## w_21                          0.043 0.043 34.060   1.014   0.318   -0.043    0.130
## r_xy1                         0.047 0.178 34.031   0.265   0.793   -0.314    0.409
## r_xy2                         0.160 0.158 34.060   1.014   0.318   -0.161    0.482
## b_11                          0.046 0.173 34.031   0.265   0.793   -0.307    0.398
## b_21                          0.164 0.162 34.060   1.014   0.318   -0.165    0.494
## main_effect                   0.028 0.044 34.022   0.635   0.530   -0.061    0.116
## moderator_effect              0.421 0.014 34.264  29.342   0.000    0.392    0.451
## interaction                   0.031 0.015 36.463   2.105   0.042    0.001    0.061
## q_b11_b21                    -0.120    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                  -0.115    NA     NA      NA      NA       NA       NA
## cross_over_point            -13.514    NA     NA      NA      NA       NA       NA
## interaction_vs_main           0.004 0.049 34.082   0.072   0.943   -0.096    0.103
## interaction_vs_main_bscale    0.013 0.185 34.082   0.072   0.943   -0.363    0.390
## interaction_vs_main_rscale    0.010 0.194 34.077   0.050   0.961   -0.384    0.403
## dadas                        -0.024 0.091 34.031  -0.265   0.604   -0.210    0.161
## dadas_bscale                 -0.092 0.347 34.031  -0.265   0.604   -0.797    0.613
## dadas_rscale                 -0.094 0.356 34.031  -0.265   0.604   -0.817    0.629
## abs_diff                      0.031 0.015 36.463   2.105   0.021    0.001    0.061
## abs_sum                       0.055 0.087 34.022   0.635   0.265   -0.122    0.232
## abs_diff_bscale               0.119 0.056 36.463   2.105   0.021    0.004    0.233
## abs_sum_bscale                0.210 0.331 34.022   0.635   0.265   -0.462    0.883
## abs_diff_rscale               0.113 0.059 36.121   1.935   0.030   -0.005    0.232
## abs_sum_rscale                0.207 0.332 34.022   0.626   0.268   -0.466    0.881
```

``` r
round(ddsc_mod2_GDI$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.003 -0.161  0.824  1.000  0.364
```

``` r
# Structural path model
d_GDI<-ddsc_mod2_GDI$ddsc_sem_fit$data

ddsc_sem_GDI<-
  ddsc_sem(data=d_GDI,x = "gdi.z.cm",y1="means_y2",
           y2="means_y1",q_sesoi = .10)
round(ddsc_sem_GDI$results,3)
```

```
##                                    est    se       z pvalue ci.lower ci.upper
## r_xy1_y2                        -0.495 0.149  -3.323  0.001   -0.787   -0.203
## r_xy1                           -0.008 0.171  -0.049  0.961   -0.344    0.328
## r_xy2                            0.166 0.169   0.980  0.327   -0.166    0.497
## b_11                            -0.009 0.176  -0.049  0.961   -0.353    0.336
## b_21                             0.161 0.165   0.980  0.327   -0.162    0.485
## b_10                            -0.656 0.173  -3.793  0.000   -0.996   -0.317
## b_20                             0.909 0.162   5.596  0.000    0.591    1.227
## res_cov_y1_y2                    0.915 0.227   4.031  0.000    0.470    1.359
## diff_b10_b20                    -1.565 0.050 -31.055  0.000   -1.664   -1.467
## diff_b11_b21                    -0.170 0.051  -3.323  0.001   -0.270   -0.070
## diff_rxy1_rxy2                  -0.174 0.050  -3.473  0.001   -0.272   -0.076
## q_b11_b21                       -0.171 0.051  -3.357  0.001   -0.272   -0.071
## q_rxy1_rxy2                     -0.176 0.051  -3.456  0.001   -0.275   -0.076
## cross_over_point                -9.206 2.786  -3.305  0.001  -14.666   -3.746
## sum_b11_b21                      0.153 0.337   0.454  0.650   -0.507    0.813
## main_effect                      0.076 0.168   0.454  0.650   -0.254    0.407
## interaction_vs_main_effect       0.094 0.186   0.503  0.615   -0.271    0.459
## diff_abs_b11_abs_b21            -0.153 0.337  -0.454  0.650   -0.813    0.507
## abs_diff_b11_b21                 0.170 0.051   3.323  0.000    0.070    0.270
## abs_sum_b11_b21                  0.153 0.337   0.454  0.325   -0.507    0.813
## dadas                            0.017 0.351   0.049  0.481   -0.672    0.706
## q_r_equivalence                  0.076 0.051   1.487  0.932       NA       NA
## q_b_equivalence                  0.071 0.051   1.399  0.919       NA       NA
## cross_over_point_equivalence     9.206 2.786   3.305  1.000       NA       NA
## cross_over_point_minimal_effect  9.206 2.786   3.305  0.000       NA       NA
```

``` r
round(ddsc_sem_GDI$variance_test,3)
```

```
##             est    se      z pvalue ci.lower ci.upper
## cov_y1y2  0.913 0.228  3.999  0.000    0.466    1.361
## var_y1    1.019 0.247  4.123  0.000    0.534    1.503
## var_y2    0.922 0.224  4.123  0.000    0.484    1.361
## var_diff  0.096 0.114  0.846  0.397   -0.127    0.319
## var_ratio 1.104 0.127  8.701  0.000    0.856    1.353
## cor_y1y2  0.942 0.019 48.936  0.000    0.904    0.980
```

``` r
## random intercept mlm with double-entries for each country (men and women)

d_GDI_long <- d_GDI %>%
  # move row names into a column
  rownames_to_column("cntry") %>%
  # pivot only the means_y1 / means_y2 columns
  pivot_longer(
    cols = c(means_y1, means_y2),
    names_to = "y",
    values_to = "means"
  ) %>%
  # create a sgender column based on y1/y2
  mutate(
    gndr.c = case_when(
      y == "means_y1" ~ -0.5,
      y == "means_y2" ~ 0.5
    )
  ) %>%
  dplyr::select(cntry,gdi.z.cm,means,gndr.c)


ddsc_mod2_GDI_ri<-
  ddsc_ml(data=data.frame(d_GDI_long),predictor = "gdi.z.cm",
          moderator = "gndr.c",
        DV = "means",lvl2_unit = "cntry",
        moderator_values = c(-0.5,0.5))
round(ddsc_mod2_GDI_ri$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                      -0.495 0.154 32.000  -3.224   0.003   -0.808   -0.182
## w_11                         -0.002 0.046 33.476  -0.049   0.961   -0.096    0.092
## w_21                          0.043 0.046 33.476   0.920   0.364   -0.051    0.137
## r_xy1                        -0.008 0.171 33.476  -0.049   0.961   -0.357    0.340
## r_xy2                         0.166 0.180 33.476   0.920   0.364   -0.201    0.532
## b_11                         -0.009 0.176 33.476  -0.049   0.961   -0.366    0.349
## b_21                          0.162 0.176 33.476   0.920   0.364   -0.196    0.519
## main_effect                   0.020 0.046 32.000   0.440   0.663   -0.073    0.113
## moderator_effect              0.412 0.014 32.000  30.128   0.000    0.384    0.440
## interaction                   0.045 0.014 32.000   3.224   0.003    0.016    0.073
## q_b11_b21                    -0.172    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                  -0.176    NA     NA      NA      NA       NA       NA
## cross_over_point             -9.206    NA     NA      NA      NA       NA       NA
## interaction_vs_main           0.025 0.048 37.856   0.516   0.609   -0.072    0.121
## interaction_vs_main_bscale    0.094 0.182 37.856   0.516   0.609   -0.274    0.461
## interaction_vs_main_rscale    0.095 0.173 38.315   0.551   0.585   -0.255    0.446
## dadas                         0.005 0.092 33.476   0.049   0.481   -0.184    0.193
## dadas_bscale                  0.017 0.351 33.476   0.049   0.481   -0.697    0.732
## dadas_rscale                  0.017 0.343 33.476   0.049   0.481   -0.680    0.714
## abs_diff                      0.045 0.014 32.000   3.224   0.001    0.016    0.073
## abs_sum                       0.040 0.091 32.000   0.440   0.331   -0.146    0.226
## abs_diff_bscale               0.170 0.053 32.000   3.224   0.001    0.063    0.278
## abs_sum_bscale                0.153 0.347 32.000   0.440   0.331   -0.555    0.861
## abs_diff_rscale               0.174 0.053 33.710   3.253   0.001    0.065    0.283
## abs_sum_rscale                0.157 0.348 32.001   0.453   0.327   -0.551    0.865
```

``` r
# country-time multilevel model


mod2_GDI_cntry_year<-
  lmer(FM.z.wt~gndr.c+gdi.z.cm+gndr.c:gdi.z.cm+
                 (gndr.c|cntry),data=diff_dat_cntry_year,REML=F,
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2_GDI_cntry_year)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: FM.z.wt ~ gndr.c + gdi.z.cm + gndr.c:gdi.z.cm + (gndr.c | cntry)
##    Data: diff_dat_cntry_year
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##    -371.2    -336.7     193.6    -387.2       538 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -3.7522 -0.6052 -0.0095  0.5865  4.5784 
## 
## Random effects:
##  Groups   Name        Variance  Std.Dev. Corr  
##  cntry    (Intercept) 0.0555416 0.23567        
##           gndr.c      0.0009091 0.03015  -0.41 
##  Residual             0.0228453 0.15115        
## Number of obs: 546, groups:  cntry, 34
## 
## Fixed effects:
##                 Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)      0.04763    0.04106 33.48256   1.160    0.254    
## gndr.c           0.42697    0.01405 30.93203  30.386   <2e-16 ***
## gdi.z.cm         0.03015    0.04215 34.97291   0.715    0.479    
## gndr.c:gdi.z.cm  0.01811    0.01697 44.75703   1.067    0.292    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c gd.z.c
## gndr.c      -0.150              
## gdi.z.cm    -0.006  0.001       
## gndr.c:gd..  0.001 -0.055 -0.126
```

``` r
getFE(mod2_GDI_cntry_year,round=3)
```

```
##                  Est.    SE     df      t     p     LL    UL
## (Intercept)     0.048 0.041 33.483  1.160 0.254 -0.036 0.131
## gndr.c          0.427 0.014 30.932 30.386 0.000  0.398 0.456
## gdi.z.cm        0.030 0.042 34.973  0.715 0.479 -0.055 0.116
## gndr.c:gdi.z.cm 0.018 0.017 44.757  1.067 0.292 -0.016 0.052
```

``` r
getVC(mod2_GDI_cntry_year)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.24 0.06
## 2    cntry      gndr.c   <NA>  0.03 0.00
## 3    cntry (Intercept) gndr.c -0.41 0.00
## 4 Residual        <NA>   <NA>  0.15 0.02
```

``` r
r2mlm(mod2_GDI,bargraph = F)
```

```
## $Decompositions
##                      total
## fixed           0.03935782
## slope variation 0.00145246
## mean variation  0.05585713
## sigma2          0.90333260
## 
## $R2s
##          total
## f   0.03935782
## v   0.00145246
## m   0.05585713
## fv  0.04081028
## fvm 0.09666740
```

``` r
ddsc_mod2_GDI_cntry_year<-
  ddsc_ml(model = mod2_GDI_cntry_year,
          predictor = "gdi.z.cm",
          moderator = "gndr.c",moderator_values = c(-0.5,0.5),
          re_cov_test = T)

round(ddsc_mod2_GDI_cntry_year$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.059        0.243        0.023     0.082 0.720      8.029 0.996   0.954
## 2        0.5         0.055        0.234        0.023     0.077 0.704      8.029 0.996   0.950
```

``` r
round(ddsc_mod2_GDI_cntry_year$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gdi.z.cm gdi.z.cm_scaled
## means_y1           0.259 0.243    1.000           1.000    0.940           0.940    0.175           0.175
## means_y1_scaled    1.036 0.970    1.000           1.000    0.940           0.940    0.175           0.175
## means_y2          -0.162 0.258    0.940           0.940    1.000           1.000    0.052           0.052
## means_y2_scaled   -0.649 1.029    0.940           0.940    1.000           1.000    0.052           0.052
## gdi.z.cm           0.000 1.000    0.175           0.175    0.052           0.052    1.000           1.000
## gdi.z.cm_scaled    0.000 1.000    0.175           0.175    0.052           0.052    1.000           1.000
## diff_score         0.422 0.088    0.005           0.005   -0.335          -0.335    0.333           0.333
## diff_score_scaled  1.685 0.350    0.005           0.005   -0.335          -0.335    0.333           0.333
##                   diff_score diff_score_scaled
## means_y1               0.005             0.005
## means_y1_scaled        0.005             0.005
## means_y2              -0.335            -0.335
## means_y2_scaled       -0.335            -0.335
## gdi.z.cm               0.333             0.333
## gdi.z.cm_scaled        0.333             0.333
## diff_score             1.000             1.000
## diff_score_scaled      1.000             1.000
```

``` r
round(ddsc_mod2_GDI_cntry_year$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                      -0.206 0.193 44.757  -1.067   0.292   -0.596    0.183
## w_11                          0.021 0.044 35.341   0.479   0.635   -0.068    0.110
## w_21                          0.039 0.042 35.579   0.935   0.356   -0.046    0.124
## r_xy1                         0.087 0.181 35.341   0.479   0.635   -0.281    0.455
## r_xy2                         0.152 0.163 35.579   0.935   0.356   -0.178    0.482
## b_11                          0.084 0.176 35.341   0.479   0.635   -0.273    0.441
## b_21                          0.157 0.168 35.579   0.935   0.356   -0.183    0.497
## main_effect                   0.030 0.042 34.973   0.715   0.479   -0.055    0.116
## moderator_effect              0.427 0.014 30.932  30.386   0.000    0.398    0.456
## interaction                   0.018 0.017 44.757   1.067   0.292   -0.016    0.052
## q_b11_b21                    -0.073    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                  -0.066    NA     NA      NA      NA       NA       NA
## cross_over_point            -23.575    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.012 0.047 36.306  -0.254   0.801   -0.108    0.084
## interaction_vs_main_bscale   -0.048 0.189 36.306  -0.254   0.801   -0.432    0.336
## interaction_vs_main_rscale   -0.054 0.199 36.212  -0.272   0.787   -0.458    0.350
## dadas                        -0.042 0.088 35.341  -0.479   0.683   -0.221    0.136
## dadas_bscale                 -0.169 0.352 35.341  -0.479   0.683   -0.883    0.545
## dadas_rscale                 -0.174 0.363 35.341  -0.479   0.683   -0.910    0.562
## abs_diff                      0.018 0.017 44.757   1.067   0.146   -0.016    0.052
## abs_sum                       0.060 0.084 34.973   0.715   0.240   -0.111    0.231
## abs_diff_bscale               0.072 0.068 44.757   1.067   0.146   -0.064    0.209
## abs_sum_bscale                0.241 0.337 34.973   0.715   0.240   -0.443    0.925
## abs_diff_rscale               0.065 0.070 43.871   0.934   0.178   -0.076    0.206
## abs_sum_rscale                0.239 0.337 34.971   0.708   0.242   -0.446    0.924
```

``` r
round(ddsc_mod2_GDI_cntry_year$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.002 -0.306  0.377  1.000  0.539
```

### Comparisons across models


``` r
# full multilevel model (individuals as entries)
round(ddsc_mod2_GDI$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df t.ratio p.value ci.lower ci.upper
## r_xy1              0.0471 0.1779 34.0312  0.2645  0.7930  -0.3144   0.4085
## r_xy2              0.1604 0.1582 34.0605  1.0139  0.3178  -0.1611   0.4819
## b_11               0.0459 0.1735 34.0312  0.2645  0.7930  -0.3066   0.3984
## b_21               0.1644 0.1621 34.0605  1.0139  0.3178  -0.1651   0.4939
## main_effect        0.0277 0.0436 34.0222  0.6352  0.5296  -0.0609   0.1162
## moderator_effect   0.4215 0.0144 34.2639 29.3417  0.0000   0.3923   0.4506
## interaction        0.0312 0.0148 36.4630  2.1047  0.0423   0.0011   0.0612
## q_b11_b21         -0.1200     NA      NA      NA      NA       NA       NA
```

``` r
# country-level structural path model (country-gender means as entries)
round(ddsc_sem_GDI$results[c("r_xy1","r_xy2","b_11","b_21","q_b11_b21","diff_b11_b21"),],4)
```

```
##                  est     se       z pvalue ci.lower ci.upper
## r_xy1        -0.0083 0.1715 -0.0487 0.9612  -0.3445   0.3278
## r_xy2         0.1657 0.1691  0.9795 0.3273  -0.1658   0.4971
## b_11         -0.0086 0.1757 -0.0487 0.9612  -0.3529   0.3358
## b_21          0.1615 0.1649  0.9795 0.3273  -0.1617   0.4847
## q_b11_b21    -0.1715 0.0511 -3.3569 0.0008  -0.2716  -0.0714
## diff_b11_b21 -0.1700 0.0512 -3.3234 0.0009  -0.2703  -0.0698
```

``` r
# multilevel model (country-gender means as entries)
round(ddsc_mod2_GDI_ri$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df t.ratio p.value ci.lower ci.upper
## r_xy1             -0.0083 0.1714 33.4756 -0.0487  0.9615  -0.3569   0.3402
## r_xy2              0.1657 0.1801 33.4756  0.9196  0.3644  -0.2006   0.5320
## b_11              -0.0086 0.1757 33.4756 -0.0487  0.9615  -0.3658   0.3487
## b_21               0.1615 0.1757 33.4756  0.9196  0.3644  -0.1957   0.5188
## main_effect        0.0201 0.0457 32.0000  0.4404  0.6626  -0.0730   0.1132
## moderator_effect   0.4121 0.0137 32.0000 30.1282  0.0000   0.3843   0.4400
## interaction        0.0448 0.0139 32.0000  3.2242  0.0029   0.0165   0.0730
## q_b11_b21         -0.1715     NA      NA      NA      NA       NA       NA
```

``` r
# multilevel model (country-year-gender means as entries)
round(ddsc_mod2_GDI_cntry_year$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df t.ratio p.value ci.lower ci.upper
## r_xy1              0.0869 0.1813 35.3406  0.4792  0.6348  -0.2811   0.4549
## r_xy2              0.1521 0.1627 35.5790  0.9350  0.3561  -0.1780   0.4822
## b_11               0.0843 0.1759 35.3406  0.4792  0.6348  -0.2727   0.4413
## b_21               0.1567 0.1676 35.5790  0.9350  0.3561  -0.1833   0.4966
## main_effect        0.0302 0.0421 34.9729  0.7154  0.4791  -0.0554   0.1157
## moderator_effect   0.4270 0.0141 30.9320 30.3864  0.0000   0.3983   0.4556
## interaction        0.0181 0.0170 44.7570  1.0673  0.2915  -0.0161   0.0523
## q_b11_b21         -0.0735     NA      NA      NA      NA       NA       NA
```

### Bootstrap and equivalence test

Takes a lot of time


``` r
t1<-Sys.time()
mod2_GDI_booted_fixef <-
  lme4::bootMer(
    x = mod2_GDI,
    FUN = lme4::fixef,
    nsim = 1000,
    use.u = FALSE,
    seed = 12345,
    type = c("parametric"),
    verbose = FALSE
  )
t2<-Sys.time()
t2-t1
```

```
## Time difference of 1.505135 hours
```



``` r
# obtain all the bootstrap estimates
mod2_GDI_boot_est <- data.frame(mod2_GDI_booted_fixef$t)

# calculate estimates
mod2_GDI_boot_est$w11<-mod2_GDI_boot_est$gdi.z.cm+(-0.5)*mod2_GDI_boot_est$gndr.c.gdi.z.cm
mod2_GDI_boot_est$w21<-mod2_GDI_boot_est$gdi.z.cm+(0.5)*mod2_GDI_boot_est$gndr.c.gdi.z.cm
mod2_GDI_boot_est$b11<-mod2_GDI_boot_est$w11/ddsc_mod2_GDI$SDs["SD_pooled"]
mod2_GDI_boot_est$b21<-mod2_GDI_boot_est$w21/ddsc_mod2_GDI$SDs["SD_pooled"]
mod2_GDI_boot_est$r_xy1<-mod2_GDI_boot_est$w11/ddsc_mod2_GDI$SDs["SD_y1"]
mod2_GDI_boot_est$r_xy2<-mod2_GDI_boot_est$w21/ddsc_mod2_GDI$SDs["SD_y2"]
mod2_GDI_boot_est$q_b<-atanh(mod2_GDI_boot_est$b11)-atanh(mod2_GDI_boot_est$b21)
mod2_GDI_boot_est$q<-atanh(mod2_GDI_boot_est$r_xy1)-atanh(mod2_GDI_boot_est$r_xy2)

# Calculate bootstrap summary statistics
mod2_GDI_boot_results <- t(as.data.frame(sapply(
  mod2_GDI_boot_est,
  function(x) {
    c(
      Estimate = mean(x, na.rm = TRUE),
      SE = stats::sd(x, na.rm = TRUE),
      stats::quantile(x, c((1 - .95) / 2,
                           1 - (1 - .95) / 2), na.rm = TRUE)
    )
  }
)))

mod2_GDI_boot_results
```

```
##                    Estimate         SE         2.5%        97.5%
## X.Intercept.     0.05769621 0.04335425 -0.022138635  0.142726873
## gndr.c           0.42212275 0.01488710  0.392756542  0.450242318
## gdi.z.cm         0.02962683 0.04500938 -0.061950976  0.113874925
## gndr.c.gdi.z.cm  0.03034924 0.01517715  0.001460714  0.061871452
## w11              0.01445221 0.04737138 -0.084209606  0.102700365
## w21              0.04480145 0.04384990 -0.043632398  0.125261193
## b11              0.05491445 0.17999833 -0.319973529  0.390233370
## b21              0.17023329 0.16661766 -0.165791208  0.475958362
## r_xy1            0.05631288 0.18458211 -0.328121876  0.400170932
## r_xy2            0.16610827 0.16258026 -0.161773831  0.464425156
## q_b             -0.12056558 0.06020384 -0.240830772 -0.006812653
## q               -0.11436730 0.06249212 -0.240562595  0.005111928
```

``` r
# equivalence test for q_b
tost_z(est=mod2_GDI_boot_results["q_b","Estimate"],
       se=mod2_GDI_boot_results["q_b","SE"],
       low=-.10,
       high=.10, alpha = 0.05)
```

```
## $estimate
## [1] -0.1205656
## 
## $se
## [1] 0.06020384
## 
## $low_bound
## [1] -0.1
## 
## $high_bound
## [1] 0.1
## 
## $alpha
## [1] 0.05
## 
## $z_low
## [1] -0.3415992
## 
## $p_low
## [1] 0.6336737
## 
## $z_high
## [1] -3.663647
## 
## $p_high
## [1] 0.0001243249
## 
## $ci_level
## [1] 0.9
## 
## $ci_lower
## [1] -0.2195921
## 
## $ci_upper
## [1] -0.02153908
## 
## $equivalent
## [1] FALSE
```

``` r
# equivalence test for q
tost_z(est=mod2_GDI_boot_results["q","Estimate"],
       se=mod2_GDI_boot_results["q","SE"],
       low=-.10,
       high=.10, alpha = 0.05)
```

```
## $estimate
## [1] -0.1143673
## 
## $se
## [1] 0.06249212
## 
## $low_bound
## [1] -0.1
## 
## $high_bound
## [1] 0.1
## 
## $alpha
## [1] 0.05
## 
## $z_low
## [1] -0.2299058
## 
## $p_low
## [1] 0.5909175
## 
## $z_high
## [1] -3.430309
## 
## $p_high
## [1] 0.0003014467
## 
## $ci_level
## [1] 0.9
## 
## $ci_lower
## [1] -0.2171577
## 
## $ci_upper
## [1] -0.01157691
## 
## $equivalent
## [1] FALSE
```



### Figure


``` r
# plot the results

# start with obtaining predicted values for means and differences


mod2_GDI_unstd<-lmer(FM.z~gndr.c+gdi.cm+gndr.c:gdi.cm+
                       (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
                     control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

mod2_GDI_unstd_red<-lmer(FM.z~gndr.c+
                           (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
                         control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))


p<-
  emmip(
    mod2_GDI_unstd, 
    gndr.c ~ gdi.cm,
    at=list(gndr.c = c(-0.5,0.5),
            gdi.cm=
              seq(from=round(range(GDI_country_means$gdi.cm,na.rm=T)[1],2),
                  to=round(range(GDI_country_means$gdi.cm,na.rm=T)[2],2),
                  by=0.001)),
    plotit=F,CIs=T,lmerTest.limit = 1e6,disable.pbkrtest=T)

p$gndr.c<-p$tvar
levels(p$gndr.c)<-c("Women","Men")

# obtain min and max for aligned plots
min.y.comp<-min(p$LCL)
max.y.comp<-max(p$UCL)

# Men and Women mean distributions

p3<-coefficients(mod2_GDI_unstd_red)$cntry
p3<-cbind(rbind(p3,p3),weight=rep(c(-0.5,0.5),each=nrow(p3)))
p3$xvar<-p3$`(Intercept)`+p3$gndr.c*p3$weight
p3$gndr.c<-as.factor(p3$weight)
levels(p3$gndr.c)<-c("Women","Men")

# obtain min and max for aligned plots
min.y.mean.distr<-min(p3$xvar)
max.y.mean.distr<-max(p3$xvar)


# obtain the coefs for the gndr.c-effect (difference) as function of gdi.cm

p2<-data.frame(
  emtrends(mod2_GDI_unstd,var="gndr.c",
           specs="gdi.cm",
           at=list(#gndr.c = c(-0.5,0.5),
             gdi.cm=
               seq(from=round(range(GDI_country_means$gdi.cm,na.rm=T)[1],2),
                   to=round(range(GDI_country_means$gdi.cm,na.rm=T)[2],2),
                   by=0.001)),
           lmerTest.limit = 1e6,disable.pbkrtest=T))

p2$yvar<-p2$gndr.c.trend
p2$xvar<-p2$gdi.cm
p2$LCL<-p2$lower.CL
p2$UCL<-p2$upper.CL

# obtain min and max for aligned plots
min.y.diff<-min(p2$LCL)
max.y.diff<-max(p2$UCL)

# difference score distribution

p4<-coefficients(mod2_GDI_unstd_red)$cntry
p4$xvar=(+1)*p4$gndr.c

# obtain mix and max for aligned plots

min.y.diff.distr<-min(p4$xvar)
max.y.diff.distr<-max(p4$xvar)

# define mins and maxs

min.y.pred<-
  ifelse(min.y.comp<min.y.mean.distr,min.y.comp,min.y.mean.distr)

max.y.pred<-
  ifelse(max.y.comp>max.y.mean.distr,max.y.comp,max.y.mean.distr)

min.y.narrow<-
  ifelse(min.y.diff<min.y.diff.distr,min.y.diff,min.y.diff.distr)

max.y.narrow<-
  ifelse(max.y.diff>max.y.diff.distr,max.y.diff,max.y.diff.distr)

# Figures 

# p1

# scaled simple effects to the plot
pvals<-round_tidy(ddsc_mod2_GDI$results[6:7,"p.value"],3)

ests<-
  round_tidy(ddsc_mod2_GDI$results[6:7,"estimate"],2)

coef1<-paste0("std. b11 = ",ests[1],", p = ",pvals[1])
coef2<-paste0("std. b21 = ",ests[2],", p = ",pvals[2])
coefs<-data.frame(gndr.c=c("Women","Men"),
                  coef=c(coef1,coef2))

coef_q<-paste0("Cohen's q = ",round_tidy(ddsc_mod2_GDI$results["q_b11_b21","estimate"],2),", p = ",
               substr(round_tidy(ddsc_mod2_GDI$results["interaction","p.value"],3),2,5))  
# prediction plot for difference score

pvals2<-round_tidy(ddsc_mod2_GDI$results["interaction","p.value"],3)

ests2<-
  round_tidy(-1*ddsc_mod2_GDI$results["r_xy1y2","estimate"],2)

coefs2<-paste0("difference score correlation = ",ests2,", p = ",pvals2)
# attempt to make a flag-plot

flag_points_GDI<-coefficients(mod2_GDI_unstd_red)$cntry

flag_points_GDI$women<-flag_points_GDI$`(Intercept)`+(-0.5)*flag_points_GDI$gndr.c
flag_points_GDI$men<-flag_points_GDI$`(Intercept)`+(0.5)*flag_points_GDI$gndr.c
flag_points_GDI$mean_level<-flag_points_GDI$`(Intercept)`
flag_points_GDI$difference<-flag_points_GDI$men-flag_points_GDI$women
flag_points_GDI$cntry<-rownames(flag_points_GDI)

flag_points_GDI<-
  left_join(x=flag_points_GDI,
            y=GDI_country_means[,c("ISO2","gdi.cm")],by=c("cntry"="ISO2"))
#flag_points_GDI

flag_points_GDI_long<-
  data.frame(mean_level=c(flag_points_GDI$women,
                          flag_points_GDI$men),
             gdi.cm=c(flag_points_GDI$gdi.cm,
                                 flag_points_GDI$gdi.cm),
             cntry=c(flag_points_GDI$cntry,
                     flag_points_GDI$cntry),
             gndr.c=rep(c("Women","Men"),each=nrow(flag_points_GDI)
             ))


p1.FM.flags<-
  ggplot(p,aes(y=yvar,x=gdi.cm,color=gndr.c))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("Gender Development Index (Average 2002-2023)")+
  #ylim(c(min.y.pred,max.y.pred))+
  #xlim(c(0.60,1.00))+
  ylab("Value male-typicality (Average 2002-2023)")+
  scale_color_manual(values=met.brewer("Archambault")[c(6,2)])+
  theme(legend.position = "top",
        legend.title=element_blank(),
        text=element_text(size=16,  family="sans"),
        panel.background = element_rect(fill = "white",
                                        #colour = "black",
                                        #size = 0.5, linetype = "solid"
        ),
        panel.grid.major.x = element_line(linewidth = 0.5, linetype = 2,
                                          colour = "gray"))+
  geom_text(data = coefs,show.legend=F,
            aes(label=coef,x=0.90,
                y=c(-0.35,-0.40),size=14,hjust="left"))+
  geom_text(inherit.aes=F,aes(x=0.90,y=-0.45,
                              label=coef_q,size=14,hjust="left"),
            show.legend=F)+
  geom_point(data=flag_points_GDI_long,size=8,alpha=0.30,
             aes(x=gdi.cm,y=mean_level))+
  geom_line(data=flag_points_GDI_long,aes(group = cntry,x=gdi.cm,y=mean_level),
            color="black",linetype=2)+
  geom_flag(data=flag_points_GDI_long,show.legend=F,
            aes(country=tolower(cntry),size=14,x=gdi.cm,y=mean_level))

#p1.FM.flags


p2.FM.flags<-ggplot(p2,aes(y=yvar,x=gdi.cm))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("Gender Development Index (Average 2002-2023)")+
  #xlim(c(0.60,1.00))+
  #ylim(c(min.y.narrow,max.y.narrow))+
  ylab("Difference in value male-typicality")+
  #scale_color_manual(values=met.brewer("Archambault")[c(6,2)])+
  theme(legend.position = "right",
        legend.title=element_blank(),
        text=element_text(size=16,  family="sans"),
        panel.background = element_rect(fill = "white",
                                        #colour = "black",
                                        #size = 0.5, linetype = "solid"
        ),
        panel.grid.major.x = element_line(size = 0.5, linetype = 2,
                                          colour = "gray"))+
  #geom_text(coef2,aes(x=0.63,y=min(p2$LCL)))
  geom_text(data = data.frame(coefs2),show.legend=F,
            aes(label=coefs2,x=0.96,
                y=c(round(min(p2$LCL),2)),size=18,hjust="left"))+
  geom_flag(data=flag_points_GDI,show.legend=F,
            aes(country=tolower(cntry),size=18,x=gdi.cm,y=difference))
```

```
## Warning: The `size` argument of `element_line()` is deprecated as of ggplot2 3.4.0.
## ℹ Please use the `linewidth` argument instead.
## This warning is displayed once per session.
## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.
```

``` r
#p2.FM.flags


pflag_comb<-
  ggarrange(p1.FM.flags,p2.FM.flags,align = "v",
            ncol=1,nrow=2,heights=c(2,1))
```

```
## Warning in geom_text(inherit.aes = F, aes(x = 0.9, y = -0.45, label = coef_q, : All aesthetics have length 1, but the data has 262 rows.
## ℹ Please consider using `annotate()` or provide this layer with data containing a single row.
```

``` r
pflag_comb
```

![](Analysis_files/figure-html/unnamed-chunk-40-1.png)<!-- -->

``` r
png(filename = 
      "../results/GDI_flags.png",
    units = "cm",
    width = 21.0,height=29.7*(3/4),res = 300)
pflag_comb
dev.off()
```

```
## png 
##   2
```


## mod2 with Gross Domestic Product (log_GDP)

* Logarithm of GDP per capita, PPP (constant 2017 international $)


``` r
mod2_log_GDP<-lmer(FM.z~gndr.c+log_gdp.z.cm+
                     gndr.c:log_gdp.z.cm+
                 (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2_log_GDP)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: FM.z ~ gndr.c + log_gdp.z.cm + gndr.c:log_gdp.z.cm + (gndr.c |      cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1306540.7 1306628.6 -653262.3 1306524.7    437735 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.3265 -0.6158  0.0016  0.6006  8.4531 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr  
##  cntry    (Intercept) 0.057500 0.23979        
##           gndr.c      0.007062 0.08403  -0.11 
##  Residual             1.016330 1.00813        
## Number of obs: 437743, groups:  cntry, 34
## 
## Fixed effects:
##                     Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)          0.05704    0.04117 33.95644   1.385    0.175    
## gndr.c               0.42184    0.01484 33.23336  28.420   <2e-16 ***
## log_gdp.z.cm        -0.07560    0.04130 33.99525  -1.830    0.076 .  
## gndr.c:log_gdp.z.cm  0.01758    0.01501 34.26466   1.171    0.250    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c lg_g..
## gndr.c      -0.103              
## lg_gdp.z.cm  0.023 -0.003       
## gndr.c:l_.. -0.003  0.004 -0.102
```

``` r
getFE(mod2_log_GDP,round=3)
```

```
##                       Est.    SE     df      t     p     LL    UL
## (Intercept)          0.057 0.041 33.956  1.385 0.175 -0.027 0.141
## gndr.c               0.422 0.015 33.233 28.420 0.000  0.392 0.452
## log_gdp.z.cm        -0.076 0.041 33.995 -1.830 0.076 -0.160 0.008
## gndr.c:log_gdp.z.cm  0.018 0.015 34.265  1.171 0.250 -0.013 0.048
```

``` r
getVC(mod2_log_GDP)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.24 0.06
## 2    cntry      gndr.c   <NA>  0.08 0.01
## 3    cntry (Intercept) gndr.c -0.11 0.00
## 4 Residual        <NA>   <NA>  1.01 1.02
```

``` r
r2mlm(mod2_log_GDP,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.042941918
## slope variation 0.001558599
## mean variation  0.051333178
## sigma2          0.904166305
## 
## $R2s
##           total
## f   0.042941918
## v   0.001558599
## m   0.051333178
## fv  0.044500517
## fvm 0.095833695
```


### Deconstructed associations


``` r
t1<-Sys.time()
ddsc_mod2_log_GDP<-
  ddsc_ml(model = mod2_log_GDP,
          predictor = "log_gdp.z.cm",
          moderator = "gndr.c",moderator_values = c(-0.5,0.5),
          re_cov_test = T)
```

```
## Warning in ddsc_ml(model = mod2_log_GDP, predictor = "log_gdp.z.cm", moderator = "gndr.c", : Predictor
## not properly standardized, SD = 1.0118945399749
```

``` r
t2<-Sys.time()
t2-t1
```

```
## Time difference of 35.13685 secs
```

``` r
round(ddsc_mod2_log_GDP$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.068        0.262        1.016     1.085 0.063   6999.706 0.998   0.998
## 2        0.5         0.062        0.248        1.016     1.078 0.057   5875.088 0.997   0.997
```

``` r
round(ddsc_mod2_log_GDP$descriptives,3)
```

```
##                          M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled log_gdp.z.cm
## means_y1             0.239 0.257    1.000           1.000    0.942           0.942       -0.248
## means_y1_scaled      0.909 0.975    1.000           1.000    0.942           0.942       -0.248
## means_y2            -0.173 0.270    0.942           0.942    1.000           1.000       -0.268
## means_y2_scaled     -0.656 1.025    0.942           0.942    1.000           1.000       -0.268
## log_gdp.z.cm        -0.024 1.012   -0.248          -0.248   -0.268          -0.268        1.000
## log_gdp.z.cm_scaled  0.000 1.000   -0.248          -0.248   -0.268          -0.268        1.000
## diff_score           0.412 0.090    0.028           0.028   -0.309          -0.309        0.093
## diff_score_scaled    1.565 0.343    0.028           0.028   -0.309          -0.309        0.093
##                     log_gdp.z.cm_scaled diff_score diff_score_scaled
## means_y1                         -0.248      0.028             0.028
## means_y1_scaled                  -0.248      0.028             0.028
## means_y2                         -0.268     -0.309            -0.309
## means_y2_scaled                  -0.268     -0.309            -0.309
## log_gdp.z.cm                      1.000      0.093             0.093
## log_gdp.z.cm_scaled               1.000      0.093             0.093
## diff_score                        0.093      1.000             1.000
## diff_score_scaled                 0.093      1.000             1.000
```

``` r
round(ddsc_mod2_log_GDP$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                      -0.194 0.166 34.265  -1.171   0.250   -0.532    0.143
## w_11                         -0.084 0.043 33.985  -1.975   0.056   -0.171    0.002
## w_21                         -0.067 0.041 33.994  -1.621   0.114   -0.151    0.017
## r_xy1                        -0.329 0.166 33.985  -1.975   0.056   -0.667    0.010
## r_xy2                        -0.248 0.153 33.994  -1.621   0.114   -0.558    0.063
## b_11                         -0.321 0.162 33.985  -1.975   0.056   -0.651    0.009
## b_21                         -0.254 0.157 33.994  -1.621   0.114   -0.572    0.064
## main_effect                  -0.076 0.041 33.995  -1.830   0.076   -0.160    0.008
## moderator_effect              0.422 0.015 33.233  28.420   0.000    0.392    0.452
## interaction                   0.018 0.015 34.265   1.171   0.250   -0.013    0.048
## q_b11_b21                    -0.073    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                  -0.089    NA     NA      NA      NA       NA       NA
## cross_over_point            -23.996    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.058 0.042 33.999  -1.366   0.181   -0.144    0.028
## interaction_vs_main_bscale   -0.220 0.161 33.999  -1.366   0.181   -0.548    0.108
## interaction_vs_main_rscale   -0.207 0.154 34.000  -1.345   0.188   -0.520    0.106
## dadas                        -0.134 0.082 33.994  -1.621   0.943   -0.301    0.034
## dadas_bscale                 -0.508 0.313 33.994  -1.621   0.943   -1.144    0.129
## dadas_rscale                 -0.495 0.306 33.994  -1.621   0.943   -1.117    0.126
## abs_diff                      0.018 0.015 34.265   1.171   0.125   -0.013    0.048
## abs_sum                       0.151 0.083 33.995   1.830   0.038   -0.017    0.319
## abs_diff_bscale               0.067 0.057 34.265   1.171   0.125   -0.049    0.183
## abs_sum_bscale                0.574 0.314 33.995   1.830   0.038   -0.063    1.212
## abs_diff_rscale               0.081 0.058 34.185   1.390   0.087   -0.037    0.200
## abs_sum_rscale                0.576 0.314 33.995   1.835   0.038   -0.062    1.215
```

``` r
round(ddsc_mod2_log_GDP$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.003 -0.161  0.824  1.000  0.364
```

``` r
# Structural path model
d_log_GDP<-ddsc_mod2_log_GDP$ddsc_sem_fit$data

ddsc_sem_log_GDP<-
  ddsc_sem(data=d_log_GDP,x = "log_gdp.z.cm",y1="means_y2",
           y2="means_y1",q_sesoi = .10)
round(ddsc_sem_log_GDP$results,3)
```

```
##                                     est     se       z pvalue ci.lower ci.upper
## r_xy1_y2                         -0.093  0.171  -0.545  0.586   -0.428    0.242
## r_xy1                            -0.268  0.165  -1.620  0.105   -0.591    0.056
## r_xy2                            -0.248  0.166  -1.496  0.135   -0.574    0.077
## b_11                             -0.274  0.169  -1.620  0.105   -0.606    0.058
## b_21                             -0.242  0.162  -1.496  0.135   -0.560    0.075
## b_10                             -0.656  0.167  -3.936  0.000   -0.983   -0.330
## b_20                              0.909  0.160   5.697  0.000    0.596    1.222
## res_cov_y1_y2                     0.849  0.213   3.990  0.000    0.432    1.266
## diff_b10_b20                     -1.565  0.058 -27.098  0.000   -1.679   -1.452
## diff_b11_b21                     -0.032  0.059  -0.545  0.586   -0.147    0.083
## diff_rxy1_rxy2                   -0.019  0.058  -0.329  0.742   -0.133    0.095
## q_b11_b21                        -0.034  0.063  -0.540  0.589   -0.158    0.090
## q_rxy1_rxy2                      -0.020  0.062  -0.329  0.743   -0.143    0.102
## cross_over_point                -49.020 90.028  -0.544  0.586 -225.471  127.431
## sum_b11_b21                      -0.516  0.326  -1.584  0.113   -1.156    0.123
## main_effect                      -0.258  0.163  -1.584  0.113   -0.578    0.061
## interaction_vs_main_effect       -0.226  0.166  -1.362  0.173   -0.552    0.099
## diff_abs_b11_abs_b21              0.032  0.059   0.545  0.586   -0.083    0.147
## abs_diff_b11_b21                  0.032  0.059   0.545  0.293   -0.083    0.147
## abs_sum_b11_b21                   0.516  0.326   1.584  0.057   -0.123    1.156
## dadas                            -0.485  0.324  -1.496  0.933   -1.119    0.150
## q_r_equivalence                  -0.080  0.062  -1.275  0.101       NA       NA
## q_b_equivalence                  -0.066  0.063  -1.039  0.149       NA       NA
## cross_over_point_equivalence     49.020 90.028   0.544  0.707       NA       NA
## cross_over_point_minimal_effect  49.020 90.028   0.544  0.293       NA       NA
```

``` r
round(ddsc_sem_log_GDP$variance_test,3)
```

```
##             est    se      z pvalue ci.lower ci.upper
## cov_y1y2  0.913 0.228  3.999  0.000    0.466    1.361
## var_y1    1.019 0.247  4.123  0.000    0.534    1.503
## var_y2    0.922 0.224  4.123  0.000    0.484    1.361
## var_diff  0.096 0.114  0.846  0.397   -0.127    0.319
## var_ratio 1.104 0.127  8.701  0.000    0.856    1.353
## cor_y1y2  0.942 0.019 48.936  0.000    0.904    0.980
```

``` r
## random intercept mlm with double-entries for each country (men and women)

d_log_GDP_long <- d_log_GDP %>%
  # move row names into a column
  rownames_to_column("cntry") %>%
  # pivot only the means_y1 / means_y2 columns
  pivot_longer(
    cols = c(means_y1, means_y2),
    names_to = "y",
    values_to = "means"
  ) %>%
  # create a sgender column based on y1/y2
  mutate(
    gndr.c = case_when(
      y == "means_y1" ~ -0.5,
      y == "means_y2" ~ 0.5
    )
  ) %>%
  dplyr::select(cntry,log_gdp.z.cm,means,gndr.c)


ddsc_mod2_log_GDP_ri<-
  ddsc_ml(data=data.frame(d_log_GDP_long),predictor = "log_gdp.z.cm",
          moderator = "gndr.c",
        DV = "means",lvl2_unit = "cntry",
        moderator_values = c(-0.5,0.5))
round(ddsc_mod2_log_GDP_ri$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                      -0.093 0.176 32.000  -0.528   0.601   -0.452    0.266
## w_11                         -0.072 0.045 34.067  -1.606   0.118   -0.164    0.019
## w_21                         -0.064 0.045 34.067  -1.419   0.165   -0.155    0.028
## r_xy1                        -0.268 0.167 34.067  -1.606   0.118   -0.606    0.071
## r_xy2                        -0.248 0.175 34.067  -1.419   0.165   -0.604    0.107
## b_11                         -0.274 0.171 34.067  -1.606   0.118   -0.621    0.073
## b_21                         -0.242 0.171 34.067  -1.419   0.165   -0.589    0.105
## main_effect                  -0.068 0.044 32.000  -1.536   0.134   -0.158    0.022
## moderator_effect              0.412 0.016 32.000  26.289   0.000    0.380    0.444
## interaction                   0.008 0.016 32.000   0.528   0.601   -0.024    0.041
## q_b11_b21                    -0.034    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                  -0.020    NA     NA      NA      NA       NA       NA
## cross_over_point            -49.020    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.060 0.047 40.142  -1.267   0.212   -0.155    0.035
## interaction_vs_main_bscale   -0.226 0.179 40.142  -1.267   0.212   -0.587    0.135
## interaction_vs_main_rscale   -0.239 0.187 39.591  -1.278   0.209   -0.617    0.139
## dadas                        -0.128 0.090 34.067  -1.419   0.917   -0.310    0.055
## dadas_bscale                 -0.485 0.342 34.067  -1.419   0.917   -1.179    0.210
## dadas_rscale                 -0.497 0.350 34.067  -1.419   0.917   -1.209    0.215
## abs_diff                      0.008 0.016 32.000   0.528   0.300   -0.024    0.041
## abs_sum                       0.136 0.088 32.000   1.536   0.067   -0.044    0.316
## abs_diff_bscale               0.032 0.060 32.000   0.528   0.300   -0.091    0.155
## abs_sum_bscale                0.517 0.336 32.000   1.536   0.067   -0.168    1.201
## abs_diff_rscale               0.019 0.061 33.220   0.313   0.378   -0.105    0.143
## abs_sum_rscale                0.516 0.336 32.001   1.534   0.067   -0.169    1.201
```

``` r
# country-time multilevel model


mod2_log_GDP_cntry_year<-
  lmer(FM.z.wt~gndr.c+log_gdp.z.cm+gndr.c:log_gdp.z.cm+
                 (gndr.c|cntry),data=diff_dat_cntry_year,REML=F,
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2_log_GDP_cntry_year)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: FM.z.wt ~ gndr.c + log_gdp.z.cm + gndr.c:log_gdp.z.cm + (gndr.c |      cntry)
##    Data: diff_dat_cntry_year
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##    -373.7    -339.3     194.9    -389.7       538 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -3.8245 -0.6131 -0.0193  0.5894  4.4991 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr  
##  cntry    (Intercept) 0.051508 0.22695        
##           gndr.c      0.000597 0.02443  -0.15 
##  Residual             0.022883 0.15127        
## Number of obs: 546, groups:  cntry, 34
## 
## Fixed effects:
##                     Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)          0.04696    0.03959 33.47498   1.186   0.2439    
## gndr.c               0.42510    0.01399 35.22414  30.397   <2e-16 ***
## log_gdp.z.cm        -0.07048    0.03991 34.14612  -1.766   0.0863 .  
## gndr.c:log_gdp.z.cm  0.01968    0.01496 38.59694   1.315   0.1961    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c lg_g..
## gndr.c      -0.044              
## lg_gdp.z.cm  0.010 -0.001       
## gndr.c:l_.. -0.001 -0.205 -0.041
```

``` r
getFE(mod2_log_GDP_cntry_year,round=3)
```

```
##                       Est.    SE     df      t     p     LL    UL
## (Intercept)          0.047 0.040 33.475  1.186 0.244 -0.034 0.127
## gndr.c               0.425 0.014 35.224 30.397 0.000  0.397 0.453
## log_gdp.z.cm        -0.070 0.040 34.146 -1.766 0.086 -0.152 0.011
## gndr.c:log_gdp.z.cm  0.020 0.015 38.597  1.315 0.196 -0.011 0.050
```

``` r
getVC(mod2_log_GDP_cntry_year)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.23 0.05
## 2    cntry      gndr.c   <NA>  0.02 0.00
## 3    cntry (Intercept) gndr.c -0.15 0.00
## 4 Residual        <NA>   <NA>  0.15 0.02
```

``` r
r2mlm(mod2_log_GDP,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.042941918
## slope variation 0.001558599
## mean variation  0.051333178
## sigma2          0.904166305
## 
## $R2s
##           total
## f   0.042941918
## v   0.001558599
## m   0.051333178
## fv  0.044500517
## fvm 0.095833695
```

``` r
ddsc_mod2_log_GDP_cntry_year<-
  ddsc_ml(model = mod2_log_GDP_cntry_year,
          predictor = "log_gdp.z.cm",
          moderator = "gndr.c",moderator_values = c(-0.5,0.5),
          re_cov_test = T)
```

```
## Warning in ddsc_ml(model = mod2_log_GDP_cntry_year, predictor = "log_gdp.z.cm", : Predictor not properly
## standardized, SD = 1.0118945399749
```

``` r
round(ddsc_mod2_log_GDP_cntry_year$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.059        0.243        0.023     0.082 0.720      8.029 0.996   0.954
## 2        0.5         0.055        0.234        0.023     0.077 0.704      8.029 0.996   0.950
```

``` r
round(ddsc_mod2_log_GDP_cntry_year$descriptives,3)
```

```
##                          M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled log_gdp.z.cm
## means_y1             0.259 0.243    1.000           1.000    0.940           0.940       -0.249
## means_y1_scaled      1.036 0.970    1.000           1.000    0.940           0.940       -0.249
## means_y2            -0.162 0.258    0.940           0.940    1.000           1.000       -0.313
## means_y2_scaled     -0.649 1.029    0.940           0.940    1.000           1.000       -0.313
## log_gdp.z.cm        -0.024 1.012   -0.249          -0.249   -0.313          -0.313        1.000
## log_gdp.z.cm_scaled  0.000 1.000   -0.249          -0.249   -0.313          -0.313        1.000
## diff_score           0.422 0.088    0.005           0.005   -0.335          -0.335        0.231
## diff_score_scaled    1.685 0.350    0.005           0.005   -0.335          -0.335        0.231
##                     log_gdp.z.cm_scaled diff_score diff_score_scaled
## means_y1                         -0.249      0.005             0.005
## means_y1_scaled                  -0.249      0.005             0.005
## means_y2                         -0.313     -0.335            -0.335
## means_y2_scaled                  -0.313     -0.335            -0.335
## log_gdp.z.cm                      1.000      0.231             0.231
## log_gdp.z.cm_scaled               1.000      0.231             0.231
## diff_score                        0.231      1.000             1.000
## diff_score_scaled                 0.231      1.000             1.000
```

``` r
round(ddsc_mod2_log_GDP_cntry_year$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                      -0.224 0.171 38.597  -1.315   0.196   -0.569    0.121
## w_11                         -0.080 0.041 34.225  -1.963   0.058   -0.163    0.003
## w_21                         -0.061 0.040 34.333  -1.505   0.141   -0.143    0.021
## r_xy1                        -0.331 0.168 34.225  -1.963   0.058   -0.673    0.012
## r_xy2                        -0.235 0.156 34.333  -1.505   0.141   -0.553    0.082
## b_11                         -0.321 0.163 34.225  -1.963   0.058   -0.653    0.011
## b_21                         -0.242 0.161 34.333  -1.505   0.141   -0.569    0.085
## main_effect                  -0.070 0.040 34.146  -1.766   0.086   -0.152    0.011
## moderator_effect              0.425 0.014 35.224  30.397   0.000    0.397    0.453
## interaction                   0.020 0.015 38.597   1.315   0.196   -0.011    0.050
## q_b11_b21                    -0.085    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                  -0.104    NA     NA      NA      NA       NA       NA
## cross_over_point            -21.606    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.051 0.042 34.756  -1.209   0.235   -0.136    0.035
## interaction_vs_main_bscale   -0.203 0.168 34.756  -1.209   0.235   -0.544    0.138
## interaction_vs_main_rscale   -0.188 0.159 34.806  -1.180   0.246   -0.510    0.135
## dadas                        -0.121 0.081 34.333  -1.505   0.929   -0.285    0.042
## dadas_bscale                 -0.485 0.322 34.333  -1.505   0.929   -1.139    0.170
## dadas_rscale                 -0.471 0.313 34.333  -1.505   0.929   -1.106    0.165
## abs_diff                      0.020 0.015 38.597   1.315   0.098   -0.011    0.050
## abs_sum                       0.141 0.080 34.146   1.766   0.043   -0.021    0.303
## abs_diff_bscale               0.079 0.060 38.597   1.315   0.098   -0.042    0.200
## abs_sum_bscale                0.563 0.319 34.146   1.766   0.043   -0.085    1.211
## abs_diff_rscale               0.095 0.061 38.291   1.567   0.063   -0.028    0.219
## abs_sum_rscale                0.566 0.319 34.145   1.773   0.043   -0.083    1.215
```

``` r
round(ddsc_mod2_log_GDP_cntry_year$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.002 -0.306  0.377  1.000  0.539
```

### Comparisons across models


``` r
# full multilevel model (individuals as entries)
round(ddsc_mod2_log_GDP$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df t.ratio p.value ci.lower ci.upper
## r_xy1             -0.3288 0.1665 33.9848 -1.9751  0.0564  -0.6671   0.0095
## r_xy2             -0.2477 0.1528 33.9945 -1.6207  0.1143  -0.5583   0.0629
## b_11              -0.3206 0.1623 33.9848 -1.9751  0.0564  -0.6506   0.0093
## b_21              -0.2538 0.1566 33.9945 -1.6207  0.1143  -0.5721   0.0645
## main_effect       -0.0756 0.0413 33.9952 -1.8303  0.0760  -0.1595   0.0083
## moderator_effect   0.4218 0.0148 33.2334 28.4198  0.0000   0.3917   0.4520
## interaction        0.0176 0.0150 34.2647  1.1714  0.2495  -0.0129   0.0481
## q_b11_b21         -0.0728     NA      NA      NA      NA       NA       NA
```

``` r
# country-level structural path model (country-gender means as entries)
round(ddsc_sem_log_GDP$results[c("r_xy1","r_xy2","b_11","b_21","q_b11_b21","diff_b11_b21"),],4)
```

```
##                  est     se       z pvalue ci.lower ci.upper
## r_xy1        -0.2676 0.1652 -1.6196 0.1053  -0.5915   0.0562
## r_xy2        -0.2485 0.1661 -1.4959 0.1347  -0.5741   0.0771
## b_11         -0.2742 0.1693 -1.6196 0.1053  -0.6060   0.0576
## b_21         -0.2423 0.1619 -1.4959 0.1347  -0.5597   0.0752
## q_b11_b21    -0.0342 0.0633 -0.5405 0.5889  -0.1583   0.0899
## diff_b11_b21 -0.0319 0.0586 -0.5446 0.5860  -0.1469   0.0830
```

``` r
# multilevel model (country-gender means as entries)
round(ddsc_mod2_log_GDP_ri$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df t.ratio p.value ci.lower ci.upper
## r_xy1             -0.2676 0.1667 34.0674 -1.6057  0.1176  -0.6063   0.0711
## r_xy2             -0.2485 0.1752 34.0674 -1.4187  0.1651  -0.6044   0.1074
## b_11              -0.2743 0.1708 34.0674 -1.6057  0.1176  -0.6214   0.0728
## b_21              -0.2423 0.1708 34.0674 -1.4187  0.1651  -0.5894   0.1048
## main_effect       -0.0680 0.0442 32.0000 -1.5364  0.1343  -0.1581   0.0221
## moderator_effect   0.4121 0.0157 32.0000 26.2890  0.0000   0.3802   0.4441
## interaction        0.0084 0.0159 32.0000  0.5283  0.6009  -0.0240   0.0408
## q_b11_b21         -0.0342     NA      NA      NA      NA       NA       NA
```

``` r
# multilevel model (country-year-gender means as entries)
round(ddsc_mod2_log_GDP_cntry_year$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df t.ratio p.value ci.lower ci.upper
## r_xy1             -0.3308 0.1685 34.2252 -1.9634  0.0578  -0.6731   0.0115
## r_xy2             -0.2353 0.1564 34.3333 -1.5049  0.1415  -0.5530   0.0823
## b_11              -0.3209 0.1635 34.2252 -1.9634  0.0578  -0.6530   0.0112
## b_21              -0.2423 0.1610 34.3333 -1.5049  0.1415  -0.5694   0.0848
## main_effect       -0.0705 0.0399 34.1461 -1.7661  0.0863  -0.1516   0.0106
## moderator_effect   0.4251 0.0140 35.2241 30.3969  0.0000   0.3967   0.4535
## interaction        0.0197 0.0150 38.5969  1.3154  0.1961  -0.0106   0.0499
## q_b11_b21         -0.0855     NA      NA      NA      NA       NA       NA
```

### Bootstrap and equivalence test

Takes a lot of time


``` r
t1<-Sys.time()
mod2_log_GDP_booted_fixef <-
  lme4::bootMer(
    x = mod2_log_GDP,
    FUN = lme4::fixef,
    nsim = 1000,
    use.u = FALSE,
    seed = 12345,
    type = c("parametric"),
    verbose = FALSE
  )
t2<-Sys.time()
t2-t1
```

```
## Time difference of 1.185969 hours
```



``` r
# obtain all the bootstrap estimates
mod2_log_GDP_boot_est <- data.frame(mod2_log_GDP_booted_fixef$t)

# calculate estimates
mod2_log_GDP_boot_est$w11<-mod2_log_GDP_boot_est$log_gdp.z.cm+(-0.5)*mod2_log_GDP_boot_est$gndr.c.log_gdp.z.cm
mod2_log_GDP_boot_est$w21<-mod2_log_GDP_boot_est$log_gdp.z.cm+(0.5)*mod2_log_GDP_boot_est$gndr.c.log_gdp.z.cm
mod2_log_GDP_boot_est$b11<-mod2_log_GDP_boot_est$w11/ddsc_mod2_log_GDP$SDs["SD_pooled"]
mod2_log_GDP_boot_est$b21<-mod2_log_GDP_boot_est$w21/ddsc_mod2_log_GDP$SDs["SD_pooled"]
mod2_log_GDP_boot_est$r_xy1<-mod2_log_GDP_boot_est$w11/ddsc_mod2_log_GDP$SDs["SD_y1"]
mod2_log_GDP_boot_est$r_xy2<-mod2_log_GDP_boot_est$w21/ddsc_mod2_log_GDP$SDs["SD_y2"]
mod2_log_GDP_boot_est$q_b<-atanh(mod2_log_GDP_boot_est$b11)-atanh(mod2_log_GDP_boot_est$b21)
mod2_log_GDP_boot_est$q<-atanh(mod2_log_GDP_boot_est$r_xy1)-atanh(mod2_log_GDP_boot_est$r_xy2)

# Calculate bootstrap summary statistics
mod2_log_GDP_boot_results <- t(as.data.frame(sapply(
  mod2_log_GDP_boot_est,
  function(x) {
    c(
      Estimate = mean(x, na.rm = TRUE),
      SE = stats::sd(x, na.rm = TRUE),
      stats::quantile(x, c((1 - .95) / 2,
                           1 - (1 - .95) / 2), na.rm = TRUE)
    )
  }
)))

mod2_log_GDP_boot_results
```

```
##                        Estimate         SE        2.5%         97.5%
## X.Intercept.         0.05597194 0.04164783 -0.02085622  0.1374367423
## gndr.c               0.42249071 0.01536107  0.39236440  0.4519651307
## log_gdp.z.cm        -0.07572176 0.04044426 -0.15493095  0.0047992536
## gndr.c.log_gdp.z.cm  0.01825901 0.01473316 -0.01008881  0.0465500903
## w11                 -0.08485126 0.04199790 -0.16759262 -0.0005271133
## w21                 -0.06659226 0.04020180 -0.14507301  0.0094079509
## b11                 -0.32241165 0.15958058 -0.63680622 -0.0020028867
## b21                 -0.25303240 0.15275587 -0.55123786  0.0357476469
## r_xy1               -0.33062208 0.16364441 -0.65302293 -0.0020538916
## r_xy2               -0.24690104 0.14905436 -0.53788052  0.0348814263
## q_b                 -0.07990374 0.06692691 -0.21663672  0.0404725835
## q                   -0.09742641 0.07172337 -0.25099206  0.0279261464
```

``` r
# equivalence test for q_b
tost_z(est=mod2_log_GDP_boot_results["q_b","Estimate"],
       se=mod2_log_GDP_boot_results["q_b","SE"],
       low=-.10,
       high=.10, alpha = 0.05)
```

```
## $estimate
## [1] -0.07990374
## 
## $se
## [1] 0.06692691
## 
## $low_bound
## [1] -0.1
## 
## $high_bound
## [1] 0.1
## 
## $alpha
## [1] 0.05
## 
## $z_low
## [1] 0.3002717
## 
## $p_low
## [1] 0.381985
## 
## $z_high
## [1] -2.688063
## 
## $p_high
## [1] 0.003593395
## 
## $ci_level
## [1] 0.9
## 
## $ci_lower
## [1] -0.1899887
## 
## $ci_upper
## [1] 0.03018123
## 
## $equivalent
## [1] FALSE
```

``` r
# equivalence test for q
tost_z(est=mod2_log_GDP_boot_results["q","Estimate"],
       se=mod2_log_GDP_boot_results["q","SE"],
       low=-.10,
       high=.10, alpha = 0.05)
```

```
## $estimate
## [1] -0.09742641
## 
## $se
## [1] 0.07172337
## 
## $low_bound
## [1] -0.1
## 
## $high_bound
## [1] 0.1
## 
## $alpha
## [1] 0.05
## 
## $z_low
## [1] 0.03588211
## 
## $p_low
## [1] 0.4856882
## 
## $z_high
## [1] -2.752609
## 
## $p_high
## [1] 0.002956121
## 
## $ci_level
## [1] 0.9
## 
## $ci_lower
## [1] -0.2154009
## 
## $ci_upper
## [1] 0.02054804
## 
## $equivalent
## [1] FALSE
```



### Figure


``` r
# plot the results

# start with obtaining predicted values for means and differences


mod2_log_GDP_unstd<-lmer(FM.z~gndr.c+log_gdp.cm+
                           gndr.c:log_gdp.cm+
                       (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
                     control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

mod2_log_GDP_unstd_red<-lmer(FM.z~gndr.c+
                           (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
                         control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))


p<-
  emmip(
    mod2_log_GDP_unstd, 
    gndr.c ~ log_gdp.cm,
    at=list(gndr.c = c(-0.5,0.5),
            log_gdp.cm=
              seq(from=round(range(GDP_country_means$log_gdp.cm,na.rm=T)[1],2),
                  to=round(range(GDP_country_means$log_gdp.cm,na.rm=T)[2],2),
                  by=0.001)),
    plotit=F,CIs=T,lmerTest.limit = 1e6,disable.pbkrtest=T)

p$gndr.c<-p$tvar
levels(p$gndr.c)<-c("Women","Men")

# obtain min and max for aligned plots
min.y.comp<-min(p$LCL)
max.y.comp<-max(p$UCL)

# Men and Women mean distributions

p3<-coefficients(mod2_log_GDP_unstd_red)$cntry
p3<-cbind(rbind(p3,p3),weight=rep(c(-0.5,0.5),each=nrow(p3)))
p3$xvar<-p3$`(Intercept)`+p3$gndr.c*p3$weight
p3$gndr.c<-as.factor(p3$weight)
levels(p3$gndr.c)<-c("Women","Men")

# obtain min and max for aligned plots
min.y.mean.distr<-min(p3$xvar)
max.y.mean.distr<-max(p3$xvar)


# obtain the coefs for the gndr.c-effect (difference) as function of log_gdp.cm

p2<-data.frame(
  emtrends(mod2_log_GDP_unstd,var="gndr.c",
           specs="log_gdp.cm",
           at=list(#gndr.c = c(-0.5,0.5),
             log_gdp.cm=
               seq(from=round(range(GDP_country_means$log_gdp.cm,na.rm=T)[1],2),
                   to=round(range(GDP_country_means$log_gdp.cm,na.rm=T)[2],2),
                   by=0.001)),
           lmerTest.limit = 1e6,disable.pbkrtest=T))

p2$yvar<-p2$gndr.c.trend
p2$xvar<-p2$log_gdp.cm
p2$LCL<-p2$lower.CL
p2$UCL<-p2$upper.CL

# obtain min and max for aligned plots
min.y.diff<-min(p2$LCL)
max.y.diff<-max(p2$UCL)

# difference score distribution

p4<-coefficients(mod2_log_GDP_unstd_red)$cntry
p4$xvar=(+1)*p4$gndr.c

# obtain mix and max for aligned plots

min.y.diff.distr<-min(p4$xvar)
max.y.diff.distr<-max(p4$xvar)

# define mins and maxs

min.y.pred<-
  ifelse(min.y.comp<min.y.mean.distr,min.y.comp,min.y.mean.distr)

max.y.pred<-
  ifelse(max.y.comp>max.y.mean.distr,max.y.comp,max.y.mean.distr)

min.y.narrow<-
  ifelse(min.y.diff<min.y.diff.distr,min.y.diff,min.y.diff.distr)

max.y.narrow<-
  ifelse(max.y.diff>max.y.diff.distr,max.y.diff,max.y.diff.distr)

# Figures 

# p1

# scaled simple effects to the plot
pvals<-round_tidy(ddsc_mod2_log_GDP$results[6:7,"p.value"],3)

ests<-
  round_tidy(ddsc_mod2_log_GDP$results[6:7,"estimate"],2)

coef1<-paste0("std. b11 = ",ests[1],", p = ",pvals[1])
coef2<-paste0("std. b21 = ",ests[2],", p = ",pvals[2])
coefs<-data.frame(gndr.c=c("Women","Men"),
                  coef=c(coef1,coef2))

coef_q<-paste0("Cohen's q = ",round_tidy(ddsc_mod2_log_GDP$results["q_b11_b21","estimate"],2),", p = ",
               substr(round_tidy(ddsc_mod2_log_GDP$results["interaction","p.value"],3),2,5))  
# prediction plot for difference score

pvals2<-round_tidy(ddsc_mod2_log_GDP$results["interaction","p.value"],3)

ests2<-
  round_tidy(-1*ddsc_mod2_log_GDP$results["r_xy1y2","estimate"],2)

coefs2<-paste0("difference score correlation = ",ests2,", p = ",pvals2)
# attempt to make a flag-plot

flag_points_log_GDP<-coefficients(mod2_log_GDP_unstd_red)$cntry

flag_points_log_GDP$women<-flag_points_log_GDP$`(Intercept)`+(-0.5)*flag_points_log_GDP$gndr.c
flag_points_log_GDP$men<-flag_points_log_GDP$`(Intercept)`+(0.5)*flag_points_log_GDP$gndr.c
flag_points_log_GDP$mean_level<-flag_points_log_GDP$`(Intercept)`
flag_points_log_GDP$difference<-flag_points_log_GDP$men-flag_points_log_GDP$women
flag_points_log_GDP$cntry<-rownames(flag_points_log_GDP)

flag_points_log_GDP<-
  left_join(x=flag_points_log_GDP,
            y=GDP_country_means[,c("ISO2","log_gdp.cm")],by=c("cntry"="ISO2"))
#flag_points_log_GDP

flag_points_log_GDP_long<-
  data.frame(mean_level=c(flag_points_log_GDP$women,
                          flag_points_log_GDP$men),
             log_gdp.cm=c(flag_points_log_GDP$log_gdp.cm,
                                 flag_points_log_GDP$log_gdp.cm),
             cntry=c(flag_points_log_GDP$cntry,
                     flag_points_log_GDP$cntry),
             gndr.c=rep(c("Women","Men"),each=nrow(flag_points_log_GDP)
             ))


p1.FM.flags<-
  ggplot(p,aes(y=yvar,x=log_gdp.cm,color=gndr.c))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("log(GDP per capita) (Average 2002-2023)")+
  #ylim(c(min.y.pred,max.y.pred))+
  #xlim(c(0.60,1.00))+
  ylab("Value male-typicality (Average 2002-2023)")+
  scale_color_manual(values=met.brewer("Archambault")[c(6,2)])+
  theme(legend.position = "top",
        legend.title=element_blank(),
        text=element_text(size=16,  family="sans"),
        panel.background = element_rect(fill = "white",
                                        #colour = "black",
                                        #size = 0.5, linetype = "solid"
        ),
        panel.grid.major.x = element_line(linewidth = 0.5, linetype = 2,
                                          colour = "gray"))+
  geom_text(data = coefs,show.legend=F,
            aes(label=coef,x=9.75,
                y=c(-0.35,-0.40),size=14,hjust="left"))+
  geom_text(inherit.aes=F,aes(x=9.75,y=-0.45,
                              label=coef_q,size=14,hjust="left"),
            show.legend=F)+
  geom_point(data=flag_points_log_GDP_long,size=8,alpha=0.30,
             aes(x=log_gdp.cm,y=mean_level))+
  geom_line(data=flag_points_log_GDP_long,aes(group = cntry,x=log_gdp.cm,y=mean_level),
            color="black",linetype=2)+
  geom_flag(data=flag_points_log_GDP_long,show.legend=F,
            aes(country=tolower(cntry),size=14,x=log_gdp.cm,y=mean_level))

p2.FM.flags<-ggplot(p2,aes(y=yvar,x=log_gdp.cm))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("log(GDP per capita) (Average 2002-2023)")+
  #xlim(c(0.60,1.00))+
  #ylim(c(min.y.narrow,max.y.narrow))+
  ylab("Difference in value male-typicality")+
  #scale_color_manual(values=met.brewer("Archambault")[c(6,2)])+
  theme(legend.position = "right",
        legend.title=element_blank(),
        text=element_text(size=16,  family="sans"),
        panel.background = element_rect(fill = "white",
                                        #colour = "black",
                                        #size = 0.5, linetype = "solid"
        ),
        panel.grid.major.x = element_line(size = 0.5, linetype = 2,
                                          colour = "gray"))+
  #geom_text(coef2,aes(x=0.63,y=min(p2$LCL)))
  geom_text(data = data.frame(coefs2),show.legend=F,
            aes(label=coefs2,x=9.75,
                y=c(round(min(p2$LCL),2)),size=18,hjust="left"))+
  geom_flag(data=flag_points_log_GDP,show.legend=F,
            aes(country=tolower(cntry),size=18,x=log_gdp.cm,y=difference))


pflag_comb<-
  ggarrange(p1.FM.flags,p2.FM.flags,align = "v",
            ncol=1,nrow=2,heights=c(2,1))
```

```
## Warning in geom_text(inherit.aes = F, aes(x = 9.75, y = -0.45, label = coef_q, : All aesthetics have length 1, but the data has 3302 rows.
## ℹ Please consider using `annotate()` or provide this layer with data containing a single row.
```

``` r
pflag_comb
```

![](Analysis_files/figure-html/unnamed-chunk-46-1.png)<!-- -->

``` r
png(filename = 
      "../results/log_GDP_flags.png",
    units = "cm",
    width = 21.0,height=29.7*(3/4),res = 300)
pflag_comb
dev.off()
```

```
## png 
##   2
```


## mod3: fixed effect of time (year)

* Year coded year-2002 (2002 is a zero-point)


``` r
mod3<-lmer(FM.z~gndr.c+essround.c+(gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod3)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: FM.z ~ gndr.c + essround.c + (gndr.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##   1305952   1306029   -652969   1305938    437736 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.2626 -0.6158  0.0017  0.6000  8.3760 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr  
##  cntry    (Intercept) 0.062785 0.25057        
##           gndr.c      0.007416 0.08611  -0.16 
##  Residual             1.014958 1.00745        
## Number of obs: 437743, groups:  cntry, 34
## 
## Fixed effects:
##               Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)  6.215e-02  4.301e-02  3.397e+01   1.445    0.158    
## gndr.c       4.218e-01  1.519e-02  3.337e+01  27.771   <2e-16 ***
## essround.c  -1.263e-02  5.194e-04  4.374e+05 -24.316   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##            (Intr) gndr.c
## gndr.c     -0.153       
## essround.c -0.003  0.000
```

``` r
getFE(mod3,round=3)
```

```
##               Est.    SE         df       t     p     LL     UL
## (Intercept)  0.062 0.043     33.966   1.445 0.158 -0.025  0.150
## gndr.c       0.422 0.015     33.370  27.771 0.000  0.391  0.453
## essround.c  -0.013 0.001 437444.841 -24.316 0.000 -0.014 -0.012
```

``` r
getVC(mod3)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.25 0.06
## 2    cntry      gndr.c   <NA>  0.09 0.01
## 3    cntry (Intercept) gndr.c -0.16 0.00
## 4 Residual        <NA>   <NA>  1.01 1.01
```

``` r
r2mlm(mod3,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.040509874
## slope variation 0.001634673
## mean variation  0.056061432
## sigma2          0.901794021
## 
## $R2s
##           total
## f   0.040509874
## v   0.001634673
## m   0.056061432
## fv  0.042144547
## fvm 0.098205979
```

``` r
anova(mod2,mod3)
```

```
## Data: diff_dat
## Models:
## mod2: FM.z ~ gndr.c + (gndr.c | cntry)
## mod3: FM.z ~ gndr.c + essround.c + (gndr.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod2    6 1306541 1306607 -653264   1306529                         
## mod3    7 1305952 1306029 -652969   1305938 590.85  1  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## mod4: random effect of time (year)


``` r
mod4<-lmer(FM.z~gndr.c+year.c+(gndr.c+year.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod4)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: FM.z ~ gndr.c + year.c + (gndr.c + year.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1301668.2 1301778.1 -650824.1 1301648.2    437733 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.5977 -0.6137  0.0028  0.5977  8.3433 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr        
##  cntry    (Intercept) 0.361543 0.60128              
##           gndr.c      0.007665 0.08755   0.23       
##           year.c      0.001759 0.04194  -0.79 -0.41 
##  Residual             1.004517 1.00226              
## Number of obs: 437743, groups:  cntry, 34
## 
## Fixed effects:
##              Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept) -0.030570   0.103357 30.782593  -0.296    0.769    
## gndr.c       0.421107   0.015423 33.290226  27.304   <2e-16 ***
## year.c       0.005096   0.007209 32.207783   0.707    0.485    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##        (Intr) gndr.c
## gndr.c  0.227       
## year.c -0.791 -0.403
```

``` r
getFE(mod4,round=3)
```

```
##               Est.    SE     df      t     p     LL    UL
## (Intercept) -0.031 0.103 30.783 -0.296 0.769 -0.241 0.180
## gndr.c       0.421 0.015 33.290 27.304 0.000  0.390 0.452
## year.c       0.005 0.007 32.208  0.707 0.485 -0.010 0.020
```

``` r
getVC(mod4)
```

```
##        grp        var1   var2 sdcor  vcov
## 1    cntry (Intercept)   <NA>  0.60  0.36
## 2    cntry      gndr.c   <NA>  0.09  0.01
## 3    cntry      year.c   <NA>  0.04  0.00
## 4    cntry (Intercept) gndr.c  0.23  0.01
## 5    cntry (Intercept) year.c -0.79 -0.02
## 6    cntry      gndr.c year.c -0.41  0.00
## 7 Residual        <NA>   <NA>  1.00  1.00
```

``` r
r2mlm(mod4,bargraph = F)
```

```
## $Decompositions
##                      total
## fixed           0.03578556
## slope variation 0.05744274
## mean variation  0.10955711
## sigma2          0.79721460
## 
## $R2s
##          total
## f   0.03578556
## v   0.05744274
## m   0.10955711
## fv  0.09322830
## fvm 0.20278540
```

``` r
anova(mod2,mod3,mod4)
```

```
## Data: diff_dat
## Models:
## mod2: FM.z ~ gndr.c + (gndr.c | cntry)
## mod3: FM.z ~ gndr.c + essround.c + (gndr.c | cntry)
## mod4: FM.z ~ gndr.c + year.c + (gndr.c + year.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L)   Chisq Df Pr(>Chisq)    
## mod2    6 1306541 1306607 -653264   1306529                          
## mod3    7 1305952 1306029 -652969   1305938  590.85  1  < 2.2e-16 ***
## mod4   10 1301668 1301778 -650824   1301648 4289.71  3  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## mod5: fixed interaction between time and gender


``` r
mod5<-lmer(FM.z~gndr.c+year.c+
             gndr.c:year.c+(gndr.c+year.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod5)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: FM.z ~ gndr.c + year.c + gndr.c:year.c + (gndr.c + year.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##   1301600   1301721   -650789   1301578    437732 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.6194 -0.6136  0.0030  0.5979  8.3154 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr        
##  cntry    (Intercept) 0.362395 0.60199              
##           gndr.c      0.007897 0.08886   0.16       
##           year.c      0.001763 0.04199  -0.79 -0.37 
##  Residual             1.004352 1.00217              
## Number of obs: 437743, groups:  cntry, 34
## 
## Fixed effects:
##                 Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)   -2.943e-02  1.035e-01  3.079e+01  -0.284    0.778    
## gndr.c         4.658e-01  1.652e-02  4.098e+01  28.194   <2e-16 ***
## year.c         4.992e-03  7.218e-03  3.221e+01   0.692    0.494    
## gndr.c:year.c -4.211e-03  5.019e-04  2.075e+05  -8.390   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c year.c
## gndr.c       0.144              
## year.c      -0.792 -0.338       
## gndr.c:yr.c -0.001 -0.322  0.002
```

``` r
getFE(mod5,round=3)
```

```
##                 Est.    SE         df      t     p     LL     UL
## (Intercept)   -0.029 0.103     30.786 -0.284 0.778 -0.241  0.182
## gndr.c         0.466 0.017     40.975 28.194 0.000  0.432  0.499
## year.c         0.005 0.007     32.207  0.692 0.494 -0.010  0.020
## gndr.c:year.c -0.004 0.001 207533.892 -8.390 0.000 -0.005 -0.003
```

``` r
getVC(mod5)
```

```
##        grp        var1   var2 sdcor  vcov
## 1    cntry (Intercept)   <NA>  0.60  0.36
## 2    cntry      gndr.c   <NA>  0.09  0.01
## 3    cntry      year.c   <NA>  0.04  0.00
## 4    cntry (Intercept) gndr.c  0.16  0.01
## 5    cntry (Intercept) year.c -0.79 -0.02
## 6    cntry      gndr.c year.c -0.37  0.00
## 7 Residual        <NA>   <NA>  1.00  1.00
```

``` r
r2mlm(mod5,bargraph = F)
```

```
## $Decompositions
##                      total
## fixed           0.03623703
## slope variation 0.05757744
## mean variation  0.10973893
## sigma2          0.79644660
## 
## $R2s
##          total
## f   0.03623703
## v   0.05757744
## m   0.10973893
## fv  0.09381447
## fvm 0.20355340
```

``` r
anova(mod4,mod5)
```

```
## Data: diff_dat
## Models:
## mod4: FM.z ~ gndr.c + year.c + (gndr.c + year.c | cntry)
## mod5: FM.z ~ gndr.c + year.c + gndr.c:year.c + (gndr.c + year.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L) Chisq Df Pr(>Chisq)    
## mod4   10 1301668 1301778 -650824   1301648                        
## mod5   11 1301600 1301721 -650789   1301578 70.13  1  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## mod6: random interaction between time and gender


``` r
mod6<-lmer(FM.z~gndr.c+year.c+
             gndr.c:year.c+(gndr.c+year.c+gndr.c:year.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: FM.z ~ gndr.c + year.c + gndr.c:year.c + (gndr.c + year.c + gndr.c:year.c |      cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1301502.6 1301667.4 -650736.3 1301472.6    437728 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.7030 -0.6137  0.0025  0.5976  8.3038 
## 
## Random effects:
##  Groups   Name          Variance Std.Dev. Corr              
##  cntry    (Intercept)   0.368778 0.60727                    
##           gndr.c        0.042779 0.20683  -0.75             
##           year.c        0.001779 0.04217  -0.79  0.56       
##           gndr.c:year.c 0.000142 0.01192   0.77 -0.87 -0.73 
##  Residual               1.003955 1.00198                    
## Number of obs: 437743, groups:  cntry, 34
## 
## Fixed effects:
##                Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)   -0.031593   0.104379 30.851697  -0.303    0.764    
## gndr.c         0.500555   0.036602 11.374954  13.676 2.06e-08 ***
## year.c         0.005133   0.007249 32.220535   0.708    0.484    
## gndr.c:year.c -0.006222   0.002146 13.533669  -2.899    0.012 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c year.c
## gndr.c      -0.728              
## year.c      -0.795  0.538       
## gndr.c:yr.c  0.731 -0.867 -0.696
```

``` r
getFE(mod6,round=3)
```

```
##                 Est.    SE     df      t     p     LL     UL
## (Intercept)   -0.032 0.104 30.852 -0.303 0.764 -0.245  0.181
## gndr.c         0.501 0.037 11.375 13.676 0.000  0.420  0.581
## year.c         0.005 0.007 32.221  0.708 0.484 -0.010  0.020
## gndr.c:year.c -0.006 0.002 13.534 -2.899 0.012 -0.011 -0.002
```

``` r
getVC(mod6)
```

```
##         grp          var1          var2 sdcor  vcov
## 1     cntry   (Intercept)          <NA>  0.61  0.37
## 2     cntry        gndr.c          <NA>  0.21  0.04
## 3     cntry        year.c          <NA>  0.04  0.00
## 4     cntry gndr.c:year.c          <NA>  0.01  0.00
## 5     cntry   (Intercept)        gndr.c -0.75 -0.09
## 6     cntry   (Intercept)        year.c -0.79 -0.02
## 7     cntry   (Intercept) gndr.c:year.c  0.77  0.01
## 8     cntry        gndr.c        year.c  0.56  0.00
## 9     cntry        gndr.c gndr.c:year.c -0.87  0.00
## 10    cntry        year.c gndr.c:year.c -0.73  0.00
## 11 Residual          <NA>          <NA>  1.00  1.00
```

``` r
r2mlm(mod6,bargraph = F)
```

```
## $Decompositions
##                      total
## fixed           0.03856817
## slope variation 0.06089767
## mean variation  0.11124088
## sigma2          0.78929328
## 
## $R2s
##          total
## f   0.03856817
## v   0.06089767
## m   0.11124088
## fv  0.09946584
## fvm 0.21070672
```

``` r
anova(mod4,mod5,mod6)
```

```
## Data: diff_dat
## Models:
## mod4: FM.z ~ gndr.c + year.c + (gndr.c + year.c | cntry)
## mod5: FM.z ~ gndr.c + year.c + gndr.c:year.c + (gndr.c + year.c | cntry)
## mod6: FM.z ~ gndr.c + year.c + gndr.c:year.c + (gndr.c + year.c + gndr.c:year.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod4   10 1301668 1301778 -650824   1301648                         
## mod5   11 1301600 1301721 -650789   1301578  70.13  1  < 2.2e-16 ***
## mod6   15 1301503 1301667 -650736   1301473 105.48  4  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

### Trends


``` r
# gender specific change over time

change_mod6<-emmeans(mod6,specs="year.c",by="gndr.c",
                     at=list(gndr.c=c(-0.5,0.5),
                             year.c=rev(range(diff_dat$year.c))),
                     disable.pbkrtest=T,
                     lmerTest.limit = 500000,infer=c(T,T))
change_mod6
```

```
## gndr.c = -0.5:
##  year.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21 -0.109 0.1000 31.5  -0.3133   0.0958  -1.084  0.2867
##       0 -0.282 0.1180 29.3  -0.5238  -0.0399  -2.381  0.0240
## 
## gndr.c =  0.5:
##  year.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21  0.261 0.0885 31.4   0.0808   0.4415   2.952  0.0059
##       0  0.219 0.0919 28.3   0.0305   0.4069   2.379  0.0244
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_mod6,adjust="none",infer=c(T,T))
```

```
## gndr.c = -0.5:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0   0.1731 0.169 31.3   -0.171    0.517   1.026  0.3127
## 
## gndr.c =  0.5:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0   0.0424 0.137 30.8   -0.238    0.323   0.309  0.7596
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# change in gender differences over time

change_in_diff_mod6<-emmeans(mod6,specs=c("gndr.c","year.c"),
                             at=list(gndr.c=c(-0.5,0.5),
                                     year.c=rev(range(diff_dat$year.c))),
                             disable.pbkrtest=T,
                             lmerTest.limit = 500000,infer=c(T,T))
change_in_diff_mod6
```

```
##  gndr.c year.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21 -0.109 0.1000 31.5  -0.3133   0.0958  -1.084  0.2867
##     0.5     21  0.261 0.0885 31.4   0.0808   0.4415   2.952  0.0059
##    -0.5      0 -0.282 0.1180 29.3  -0.5238  -0.0399  -2.381  0.0240
##     0.5      0  0.219 0.0919 28.3   0.0305   0.4069   2.379  0.0244
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod6,adjust="none",infer=c(T,T))
```

```
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21    -0.3699 0.0226 20.4   -0.417  -0.3228 -16.362 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)   0.1731 0.1690 31.3   -0.171   0.5170   1.026  0.3127
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     -0.3274 0.1500 31.7   -0.632  -0.0226  -2.189  0.0361
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)      0.5430 0.1570 31.4    0.222   0.8640   3.448  0.0016
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0         0.0424 0.1370 30.8   -0.238   0.3229   0.309  0.7596
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0      -0.5006 0.0366 11.4   -0.581  -0.4203 -13.676 <0.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# custom contrast to obtain the estimates for the entire span

changes_in_diff <- list(
  diff_ESS11 = c(-1,1,0,0),
  diff_ESS1 = c(0,0,-1,1))

diff_mod6<-contrast(change_in_diff_mod6,method = changes_in_diff,adjust="none",infer=c(T,T))

# gender differences at ESS1 and ESS11
diff_mod6
```

```
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11    0.370 0.0226 20.4    0.323    0.417  16.362 <0.0001
##  diff_ESS1     0.501 0.0366 11.4    0.420    0.581  13.676 <0.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# Test for whether and how much the differences have changed
pairs(diff_mod6,infer=c(T,T))
```

```
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1   -0.131 0.0451 13.5   -0.228  -0.0337  -2.899  0.0120
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

### Equivalence test for differences in linear change


``` r
diff_in_linear_change<-pairs(diff_mod6,infer=c(T,T))
tost_t(est = data.frame(diff_in_linear_change)[1,"estimate"],
       se = data.frame(diff_in_linear_change)[1,"SE"],
       df = data.frame(diff_in_linear_change)[1,"df"],
       low = -0.20,high = 0.20,alpha = .10)
```

```
## $estimate
## [1] -0.1306724
## 
## $se
## [1] 0.045076
## 
## $df
## [1] 13.53367
## 
## $low_bound
## [1] -0.2
## 
## $high_bound
## [1] 0.2
## 
## $alpha
## [1] 0.1
## 
## $t_low
## [1] 1.538016
## 
## $p_low
## [1] 0.07354574
## 
## $t_high
## [1] -7.335885
## 
## $p_high
## [1] 2.255241e-06
## 
## $ci_level
## [1] 0.8
## 
## $ci_lower
## [1] -0.1914046
## 
## $ci_upper
## [1] -0.06994019
## 
## $equivalent
## [1] TRUE
```


### Figure for time trends


``` r
# Figure for average patterns

# Model-based development for men and women

p_mod6<-
  emmip(
    mod6, 
    gndr.c ~ year.c,
    at=list(gndr.c = c(-0.5,0.5),
            year.c=
              unique(diff_dat$year.c)),
    plotit=F,CIs=T,lmerTest.limit = 1e6,disable.pbkrtest=T)

p_mod6$gndr.c<-p_mod6$tvar
levels(p_mod6$gndr.c)<-c("Women","Men")

p_mod6<-data.frame(p_mod6)

p_mod6$year<-p_mod6$year.c+2002

# add observed statistics as well

p_mod6$obs_mean<-NA
p_mod6$obs_sd<-NA
p_mod6$obs_n<-NA
p_mod6$obs_mean_wt<-NA
p_mod6$obs_sd_wt<-NA
p_mod6$obs_n_wt<-NA

for(i in 1:nrow(p_mod6)){
  #i=1
  #cntry_i<-p_mod6[i,"cntry"]
  year_i<-p_mod6[i,"year.c"]
  gndr_i<-as.numeric(as.character(p_mod6[i,"tvar"]))
  
  temp_diff_dat<-diff_dat %>%
    filter(#cntry==cntry_i,
      year.c==year_i,
      gndr.c==gndr_i) %>%
    dplyr::summarize(#n=n(),
      obs_n_wt=sum(pspwght),
      obs_mean_wt=weighted.mean(x=FM.z,w=pspwght),
      obs_sd_wt=sqrt(wtd.var(FM.z,pspwght)),
      obs_mean=mean(FM.z),
      obs_sd=sd(FM.z),
      obs_n=n()) 
  
  p_mod6[i,"obs_mean"]<-temp_diff_dat[,"obs_mean"]
  p_mod6[i,"obs_sd"]<-temp_diff_dat[,"obs_sd"]
  p_mod6[i,"obs_n"]<-temp_diff_dat[,"obs_n"]
  p_mod6[i,"obs_mean_wt"]<-temp_diff_dat[,"obs_mean_wt"]
  p_mod6[i,"obs_sd_wt"]<-temp_diff_dat[,"obs_sd_wt"]
  p_mod6[i,"obs_n_wt"]<-temp_diff_dat[,"obs_n_wt"]
}

# calculate confidence intervals
p_mod6$obs_mean_se<-p_mod6$obs_sd/sqrt(p_mod6$obs_n)
p_mod6$obs_mean_LL<-p_mod6$obs_mean+
  qnorm(.025)*p_mod6$obs_mean_se
p_mod6$obs_mean_UL<-p_mod6$obs_mean+
  qnorm(.975)*p_mod6$obs_mean_se

p_mod6$obs_mean_wt_se<-p_mod6$obs_sd_wt/sqrt(p_mod6$obs_n_wt)
p_mod6$obs_mean_wt_LL<-p_mod6$obs_mean_wt+
  qnorm(.025)*p_mod6$obs_mean_wt_se
p_mod6$obs_mean_wt_UL<-p_mod6$obs_mean_wt+
  qnorm(.975)*p_mod6$obs_mean_wt_se

# texts for the figure

## obtain texts for the figure

time_slopes<-data.frame(pairs(change_mod6,adjust="none",infer=c(T,T)))
time_slope_text_men<-
  paste0("Linear change from 2002 to 2023 = ",round_tidy(time_slopes[2,"estimate"],2),
         ", 95% CI [",round_tidy(time_slopes[2,"lower.CL"],2),
         ", ", round_tidy(time_slopes[2,"upper.CL"],2),"], p = ",
         round_tidy(time_slopes[2,"p.value"],3))

time_slope_text_women<-
  paste0("Linear change from 2002 to 2023 = ",round_tidy(time_slopes[1,"estimate"],2),
         ", 95% CI [",round_tidy(time_slopes[1,"lower.CL"],2),
         ", ", round_tidy(time_slopes[1,"upper.CL"],2),"], p = ",
         round_tidy(time_slopes[1,"p.value"],3))

time_slope_texts<-data.frame(
  gndr.c=c("Women","Men"),
  year=c(2002,2002),
  yvar=c(-0.45,0.55),
  label=c(time_slope_text_women,time_slope_text_men)
)

# Figure

#my_colors <- met.brewer("Cassatt2")[c(3, 8)]
my_colors <- met.brewer("Archambault")[c(6,2)]

p_time_trends<-
  ggplot(p_mod6, 
         aes(x = year, y = yvar, color = gndr.c)) +
  geom_smooth(method = "lm", se = FALSE,formula="y~x") +
  geom_point(size=8) +
  geom_point(aes(x=year,y=obs_mean_wt),size=8,shape = 1,alpha=.50)+
  geom_errorbar(aes(ymin=obs_mean_wt_LL, ymax=obs_mean_wt_UL),alpha = .50)+
  geom_text(
    data = time_slope_texts,
    aes(x = year, y = yvar, color = gndr.c,label=label),
    #inherit.aes = FALSE,
    hjust = 0,
    size = 4
  ) +
  scale_color_manual(values = my_colors) +
  scale_y_continuous(limits = c(-1.1, 1.1)) +
  scale_x_continuous(limits = c(2001, 2024),
                     breaks = c(seq(2002, 2020, 2),2023))+
  ylab("Mean-level of value male-typicality")+
  theme(legend.title=element_blank(),
        legend.text  = element_text(size = 16),
        axis.title.x = element_text(size = 18, face = "bold"),   
        axis.title.y = element_text(size = 18, face = "bold"),   
        axis.text.x  = element_text(size = 14),                  
        axis.text.y  = element_text(size = 14))

p_time_trends
```

![](Analysis_files/figure-html/unnamed-chunk-53-1.png)<!-- -->

``` r
png(filename = 
      "../results/time_trends.png",
    units = "cm",
    width = 21.0,height=29.7*(3/4),res = 300)
p_time_trends
dev.off()
```

```
## png 
##   2
```


### Figures for country-specific time trends


``` r
# COUNTRY-SPECIFIC TIME x GENDER TRENDS

# obtain country-specific coefficients

mod6_cntry_coefs<-coefficients(mod6)$cntry

# loop through each country and produce a data with all the estimated points at each wave for men and women

countries<-unique(diff_dat$cntry)
pred_list<-list()

for(unique_cntry in countries) {
  unique_cntry_rounds<-
    diff_dat %>%
    filter(cntry == unique_cntry) %>%
    pull(year.c) %>%
    unique()
  
  unique_cntry_coefs<-mod6_cntry_coefs[unique_cntry, ]
  
  unique_cntry_pred<-
    data.frame(cntry=unique_cntry,
               year.c=unique_cntry_rounds,
               gndr.c=rep(x = c(0.5,-0.5),each=length(unique_cntry_rounds)))
  pred_list[[unique_cntry]]<-unique_cntry_pred
}

pred_cntry_dat<-do.call(rbind.data.frame,pred_list)

# model based predictions for each time x country point
pred_cntry_dat$FM.z_mean<-predict(mod6,newdata=pred_cntry_dat)

pred_cntry_dat$year=pred_cntry_dat$year.c+2002

pred_cntry_dat$gender<-
  case_when(
    pred_cntry_dat$gndr.c==0.5~"men",
    pred_cntry_dat$gndr.c==-0.5~"women",
  )

range(pred_cntry_dat$FM.z_mean)
```

```
## [1] -0.9146445  0.9857775
```

``` r
# add observed statistics in addition to model-based estimates
pred_cntry_dat$obs_mean<-NA
pred_cntry_dat$obs_sd<-NA
pred_cntry_dat$obs_n<-NA
pred_cntry_dat$obs_mean_wt<-NA
pred_cntry_dat$obs_sd_wt<-NA
pred_cntry_dat$obs_n_wt<-NA

for(i in 1:nrow(pred_cntry_dat)){
  #i=1
  cntry_i<-pred_cntry_dat[i,"cntry"]
  year_i<-pred_cntry_dat[i,"year"]
  gndr_i<-pred_cntry_dat[i,"gndr.c"]
  
  temp_diff_dat<-diff_dat %>%
    filter(cntry==cntry_i,
           year==year_i,
           gndr.c==gndr_i) %>%
    dplyr::summarize(#n=n(),
      obs_n_wt=sum(pspwght),
      obs_mean_wt=weighted.mean(x=FM.z,w=pspwght),
      obs_sd_wt=sqrt(wtd.var(FM.z,pspwght)),
      obs_mean=mean(FM.z),
      obs_sd=sd(FM.z),
      obs_n=n()) 
  
  pred_cntry_dat[i,"obs_mean"]<-temp_diff_dat[,"obs_mean"]
  pred_cntry_dat[i,"obs_sd"]<-temp_diff_dat[,"obs_sd"]
  pred_cntry_dat[i,"obs_n"]<-temp_diff_dat[,"obs_n"]
  pred_cntry_dat[i,"obs_mean_wt"]<-temp_diff_dat[,"obs_mean_wt"]
  pred_cntry_dat[i,"obs_sd_wt"]<-temp_diff_dat[,"obs_sd_wt"]
  pred_cntry_dat[i,"obs_n_wt"]<-temp_diff_dat[,"obs_n_wt"]
}

# calculate confidence intervals
pred_cntry_dat$obs_mean_se<-pred_cntry_dat$obs_sd/sqrt(pred_cntry_dat$obs_n)
pred_cntry_dat$obs_mean_LL<-pred_cntry_dat$obs_mean+qnorm(.025)*pred_cntry_dat$obs_mean_se
pred_cntry_dat$obs_mean_UL<-pred_cntry_dat$obs_mean+qnorm(.975)*pred_cntry_dat$obs_mean_se

pred_cntry_dat$obs_mean_wt_se<-pred_cntry_dat$obs_sd_wt/sqrt(pred_cntry_dat$obs_n_wt)
pred_cntry_dat$obs_mean_wt_LL<-pred_cntry_dat$obs_mean_wt+qnorm(.025)*pred_cntry_dat$obs_mean_wt_se
pred_cntry_dat$obs_mean_wt_UL<-pred_cntry_dat$obs_mean_wt+qnorm(.975)*pred_cntry_dat$obs_mean_wt_se

# add GEI to figures


index_year_dat<-
  diff_dat_cntry_year %>%
  filter(gndr.c==0.5) %>%
  dplyr::select(cntry,year,gei,gggi,gdi,log_gdp)

pred_cntry_dat<-
  left_join(x=pred_cntry_dat,
            y=index_year_dat,
            by=c("cntry","year"))


# 1) Define scale for GEI
#df_ctry <- pred_cntry_dat[pred_cntry_dat$cntry == ctry, ]

gei_min <- min(pred_cntry_dat$gei, na.rm = TRUE)
gei_max <- max(pred_cntry_dat$gei, na.rm = TRUE)

# Primary axis range
y_min <- -1.1
y_max <-  1.1

# Funktion to transform GEI to primary y-axis
scale_gei_to_y <- function(x) {
  (x - gei_min) / (gei_max - gei_min) * (y_max - y_min) + y_min
}

# Reverse function: primary axis to GEI 
scale_y_to_gei <- function(x) {
  (x - y_min) / (y_max - y_min) * (gei_max - gei_min) + gei_min
}

# save the plots to same pdf-file

#my_colors <- met.brewer("Cassatt2")[c(8, 3)]
my_colors <- met.brewer("Archambault")[c(2,6)]

pdf("../results/country_specific_time_trends.pdf", width = 8, height = 6)

for (ctry in countries) {
  print(
    ggplot(pred_cntry_dat[pred_cntry_dat$cntry == ctry, ],
       aes(x = year, y = FM.z_mean, color = gender)) +
  geom_smooth(method = "lm", formula = "y ~ x", se = FALSE) +
  geom_point(size = 8) +
  geom_point(aes(x = year, y = obs_mean_wt), size = 8, shape = 1, alpha = .50) +
  geom_errorbar(aes(ymin = obs_mean_wt_LL, ymax = obs_mean_wt_UL), alpha = .50) +
  geom_line(aes(y = scale_gei_to_y(gei),
                linetype = "GEI"),   
            color = "black",         
            linewidth = 1.2) +
  scale_color_manual(values = my_colors) +
  scale_linetype_manual(          
    name   = NULL,                
    values = c("GEI" = "solid"),
    labels = c("GEI" = "GEI")
  ) +
  geom_flag(aes(country = tolower(ctry))) +
  scale_y_continuous(
    limits = c(y_min, y_max),
    name   = "Mean-level of value male-typicality",
    sec.axis = sec_axis(~ scale_y_to_gei(.),
                        name = "Gender Equality Index (GEI)")
  ) +
  scale_x_continuous(
    limits = c(2001, 2024),
    breaks = c(seq(2002, 2020, 2), 2023)
  ) +
  ggtitle(paste("Country:", ctry)) +
  theme(legend.title = element_blank())

  )
}
```

```
## Warning: Removed 12 rows containing missing values or values outside the scale range (`geom_line()`).
```

``` r
dev.off()
```

```
## png 
##   2
```

### Figure for country-specific trends in panels


``` r
#ISO<-read.csv2("../data/ISO.csv")

pred_cntry_dat<-left_join(
  x=pred_cntry_dat,
  y=ISO[,c("ISO2","Country_eng_short","CLDR")],
  by=c("cntry"="ISO2")
)

#my_colors <- met.brewer("Cassatt2")[c(8, 3)]
my_colors <- met.brewer("Archambault")[c(2,6)]

facet_plot<-
  ggplot(pred_cntry_dat, 
         aes(x = year, y = FM.z_mean, color = gender)) +
  geom_smooth(method = "lm", se = FALSE, formula="y~x") +
  geom_point(size=3.5) +
  geom_point(aes(x=year,y=obs_mean_wt),size=3.5,shape = 1,alpha=.70)+
  geom_errorbar(aes(ymin=obs_mean_wt_LL, ymax=obs_mean_wt_UL),alpha = .50)+
  geom_line(aes(y = scale_gei_to_y(gei),
                linetype = "GEI"),   
            color = "black",         
            linewidth = 1.2) +
  scale_color_manual(values = my_colors) +
  scale_linetype_manual(          
    name   = NULL,                
    values = c("GEI" = "solid"),
    labels = c("GEI" = "GEI")
  ) +
  #geom_flag(aes(country=tolower(ctry)))+
  scale_y_continuous(
    limits = c(y_min, y_max),
    name   = "Mean-level of value male-typicality",
    sec.axis = sec_axis(~ scale_y_to_gei(.),
                        name = "Gender Equality Index (GEI)")
  ) +
  scale_x_continuous(limits = c(2001, 2024),
                     breaks = c(seq(2002, 2020, 2),2023)) +
  #ggtitle(paste("Country:", ctry))+
  ylab("Mean-level of value male-typicality")+
  xlab("Year")+
  theme(legend.title=element_blank(),legend.position = "top",
        axis.text.x = element_text(angle = 45,size = 6,hjust=1))+
  facet_wrap(~CLDR,nrow=6,ncol=6)+
  #facet_wrap(~cntry,nrow=6,ncol=6)+
  geom_flag(aes(country=tolower(cntry)),size=2)

facet_plot
```

```
## Warning: Removed 12 rows containing missing values or values outside the scale range (`geom_line()`).
```

![](Analysis_files/figure-html/unnamed-chunk-55-1.png)<!-- -->

``` r
png(filename = 
      "../results/country_time_trend_facets.png",
    units = "cm",
    width = 21.0,height=29.7*(3/4),res = 600)
facet_plot
```

```
## Warning: Removed 12 rows containing missing values or values outside the scale range (`geom_line()`).
```

``` r
dev.off()
```

```
## png 
##   2
```

### Figure for country-specific estimates of change and gender difference in change


``` r
cntry_specific_changes<-
  coefficients(mod6)$cntry %>%
  mutate(change_per_year=year.c,
         gndr_change_per_year=`gndr.c:year.c`) %>%
  mutate(men_change_per_year=change_per_year+0.5*gndr_change_per_year,
         women_change_per_year=change_per_year+(-0.5)*gndr_change_per_year) %>%
  mutate(change_per_21_years=21*change_per_year,
         gndr_change_per_21_year=21*gndr_change_per_year,
         men_change_per_21_years=21*men_change_per_year,
         women_change_per_21_years=21*women_change_per_year) %>%
  dplyr::select(gndr.c,change_per_21_years,gndr_change_per_21_year,men_change_per_21_years,women_change_per_21_years) %>%
  round(.,2)

cntry_specific_changes$cntry<-rownames(cntry_specific_changes)

cntry_specific_changes<-
  left_join(x=cntry_specific_changes,
            y=n_rounds,
            by="cntry")


changes_plot_data<- data.frame(
  type=rep(c("Linear change from 2002 to 2023",
             "Gender differences in linear change"),
           each=nrow(cntry_specific_changes)),
  estimate=c(cntry_specific_changes$change_per_21_years,
             cntry_specific_changes$gndr_change_per_21_year),
  cntry=rep(cntry_specific_changes$cntry,times=2),
  nrounds=rep(cntry_specific_changes$n_unique_essround,times=2)
)

# exclude countries with less than five rounds of data
changes_plot_data_5 <- subset(changes_plot_data, nrounds > 4)
# refactor to have the correct order
changes_plot_data_5$type <- 
  factor(changes_plot_data_5$type,
         levels=c("Linear change from 2002 to 2023",
                  "Gender differences in linear change"))

# include country names
changes_plot_data_5<-left_join(
  x=changes_plot_data_5,
  y=ISO[,c("ISO2","Country_eng_short","CLDR")],
  by=c("cntry"="ISO2")
)


lambda_plot<-
  ggplot(changes_plot_data_5,aes(x=estimate,y=CLDR))+
  geom_point(size=7)+
  geom_flag(data=changes_plot_data_5,show.legend=F,
            aes(country=tolower(cntry),size=14))+
  facet_grid(~type)+
  scale_y_discrete(limits = rev)+
  ylab(NULL)+
  xlab("Estimate (SD units)")+
  geom_vline(xintercept = 0, colour = "black", linetype = "solid") +
  geom_vline(xintercept = c(-0.20, 0.20), colour = "black", linetype = "dashed")+
  scale_x_continuous(
    breaks = c(-0.5,-0.2,0,0.2,0.5, 0.8)
  )+theme(
    strip.text = element_text(size = 14),   
    axis.text  = element_text(size = 12),    
    axis.title = element_text(size = 14)  
  )

lambda_plot
```

![](Analysis_files/figure-html/unnamed-chunk-56-1.png)<!-- -->

``` r
png(filename = 
      "../results/changes_by_country.png",
    units = "cm",
    width = 21.0,height=29.7*(3/4),res = 300)
lambda_plot
dev.off()
```

```
## png 
##   2
```


### Exploratory results

Ratio of variances for time main effect and gendered time effect


``` r
getVC(mod6,round = 10)[3,"vcov"]/
  getVC(mod6,round = 10)[4,"vcov"]
```

```
## [1] 12.52648
```

Country-specific coefficients for time effect


``` r
cntry_specific_changes<-
  coefficients(mod6)$cntry %>%
  mutate(change_per_year=year.c,
         gndr_change_per_year=`gndr.c:year.c`) %>%
  mutate(men_change_per_year=change_per_year+0.5*gndr_change_per_year,
         women_change_per_year=change_per_year+(-0.5)*gndr_change_per_year) %>%
  mutate(change_per_21_years=21*change_per_year,
         gndr_change_per_21_year=21*gndr_change_per_year,
         men_change_per_21_years=21*men_change_per_year,
         women_change_per_21_years=21*women_change_per_year) %>%
  dplyr::select(gndr.c,change_per_21_years,gndr_change_per_21_year,men_change_per_21_years,women_change_per_21_years) %>%
  round(.,2)

cntry_specific_changes$cntry<-rownames(cntry_specific_changes)

cntry_specific_changes<-
  left_join(x=cntry_specific_changes,
            y=n_rounds,
            by="cntry")

cntry_specific_changes
```

```
##    gndr.c change_per_21_years gndr_change_per_21_year men_change_per_21_years women_change_per_21_years
## 1    0.57               -0.26                   -0.13                   -0.32                     -0.19
## 2    0.39               -0.07                   -0.03                   -0.08                     -0.05
## 3    0.49                0.81                   -0.04                    0.79                      0.83
## 4    0.46               -0.06                   -0.04                   -0.08                     -0.04
## 5    0.28               -0.27                    0.14                   -0.20                     -0.34
## 6    0.73                0.42                   -0.33                    0.25                      0.58
## 7    0.56               -0.50                   -0.14                   -0.57                     -0.43
## 8    0.48               -0.24                    0.00                   -0.24                     -0.24
## 9    0.53               -0.18                   -0.02                   -0.19                     -0.17
## 10   0.39               -0.50                   -0.08                   -0.54                     -0.46
## 11   0.66               -0.38                   -0.21                   -0.48                     -0.27
## 12   0.49               -0.31                   -0.04                   -0.33                     -0.30
## 13   0.53               -0.41                   -0.13                   -0.47                     -0.34
## 14   0.30                0.09                    0.03                    0.10                      0.08
## 15   0.47               -0.16                   -0.11                   -0.22                     -0.11
## 16   0.45                0.50                   -0.25                    0.37                      0.63
## 17   0.38               -0.15                    0.02                   -0.14                     -0.16
## 18   0.32                0.07                   -0.06                    0.04                      0.10
## 19   0.64               -0.23                   -0.27                   -0.37                     -0.09
## 20   0.43                0.62                   -0.18                    0.53                      0.71
## 21   0.34               -0.01                    0.17                    0.07                     -0.10
## 22   0.55               -0.45                    0.03                   -0.44                     -0.47
## 23   1.44                4.14                   -1.25                    3.52                      4.77
## 24   0.51               -0.31                   -0.11                   -0.36                     -0.25
## 25   0.52               -0.31                   -0.17                   -0.40                     -0.22
## 26   0.50               -0.16                   -0.06                   -0.18                     -0.13
## 27   0.31               -0.31                   -0.03                   -0.33                     -0.30
## 28   0.68               -0.01                   -0.25                   -0.14                      0.11
## 29   0.49                0.39                   -0.42                    0.18                      0.60
## 30   0.56               -0.44                   -0.18                   -0.53                     -0.35
## 31   0.44               -0.35                   -0.04                   -0.37                     -0.34
## 32   0.58                0.31                   -0.27                    0.18                      0.45
## 33   0.30                2.44                   -0.27                    2.30                      2.57
## 34   0.27               -0.06                    0.26                    0.07                     -0.19
##    cntry n_unique_essround
## 1     AT                 7
## 2     BE                11
## 3     BG                 7
## 4     CH                11
## 5     CY                 6
## 6     CZ                 9
## 7     DE                10
## 8     DK                 8
## 9     EE                10
## 10    ES                10
## 11    FI                11
## 12    FR                11
## 13    GB                11
## 14    GR                 6
## 15    HR                 5
## 16    HU                11
## 17    IE                11
## 18    IL                 7
## 19    IS                 6
## 20    IT                 5
## 21    LT                 7
## 22    LV                 3
## 23    ME                 3
## 24    NL                11
## 25    NO                11
## 26    PL                10
## 27    PT                11
## 28    RS                 2
## 29    RU                 5
## 30    SE                10
## 31    SI                11
## 32    SK                 8
## 33    TR                 2
## 34    UA                 6
```

``` r
# rank by overall change
cntry_specific_changes %>%
  filter(n_unique_essround>4) %>%
  dplyr::select(cntry,change_per_21_years) %>%
  arrange(change_per_21_years)
```

```
##    cntry change_per_21_years
## 1     DE               -0.50
## 2     ES               -0.50
## 3     SE               -0.44
## 4     GB               -0.41
## 5     FI               -0.38
## 6     SI               -0.35
## 7     FR               -0.31
## 8     NL               -0.31
## 9     NO               -0.31
## 10    PT               -0.31
## 11    CY               -0.27
## 12    AT               -0.26
## 13    DK               -0.24
## 14    IS               -0.23
## 15    EE               -0.18
## 16    HR               -0.16
## 17    PL               -0.16
## 18    IE               -0.15
## 19    BE               -0.07
## 20    CH               -0.06
## 21    UA               -0.06
## 22    LT               -0.01
## 23    IL                0.07
## 24    GR                0.09
## 25    SK                0.31
## 26    RU                0.39
## 27    CZ                0.42
## 28    HU                0.50
## 29    IT                0.62
## 30    BG                0.81
```

``` r
# rank by gendered change
cntry_specific_changes %>%
  filter(n_unique_essround>4) %>%
  dplyr::select(cntry,gndr_change_per_21_year) %>%
  arrange(gndr_change_per_21_year)
```

```
##    cntry gndr_change_per_21_year
## 1     RU                   -0.42
## 2     CZ                   -0.33
## 3     IS                   -0.27
## 4     SK                   -0.27
## 5     HU                   -0.25
## 6     FI                   -0.21
## 7     IT                   -0.18
## 8     SE                   -0.18
## 9     NO                   -0.17
## 10    DE                   -0.14
## 11    AT                   -0.13
## 12    GB                   -0.13
## 13    HR                   -0.11
## 14    NL                   -0.11
## 15    ES                   -0.08
## 16    IL                   -0.06
## 17    PL                   -0.06
## 18    BG                   -0.04
## 19    CH                   -0.04
## 20    FR                   -0.04
## 21    SI                   -0.04
## 22    BE                   -0.03
## 23    PT                   -0.03
## 24    EE                   -0.02
## 25    DK                    0.00
## 26    IE                    0.02
## 27    GR                    0.03
## 28    CY                    0.14
## 29    LT                    0.17
## 30    UA                    0.26
```


#### Country-level predictors: GEI


``` r
mod6_GEI<-lmer(FM.z~gndr.c+year.c+
             gndr.c:year.c+
               gei.z.cm:gndr.c+gei.z.cm:year.c+gei.z.cm:gndr.c:year.c+
               (gndr.c+year.c+gndr.c:year.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_GEI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: FM.z ~ gndr.c + year.c + gndr.c:year.c + gei.z.cm:gndr.c + gei.z.cm:year.c +  
##     gei.z.cm:gndr.c:year.c + (gndr.c + year.c + gndr.c:year.c |      cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1263210.2 1263407.6 -631587.1 1263174.2    426946 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.7505 -0.6146  0.0031  0.5989  8.3464 
## 
## Random effects:
##  Groups   Name          Variance  Std.Dev. Corr              
##  cntry    (Intercept)   0.3788676 0.61552                    
##           gndr.c        0.0420471 0.20505  -0.79             
##           year.c        0.0014005 0.03742  -0.86  0.74       
##           gndr.c:year.c 0.0001312 0.01145   0.80 -0.89 -0.83 
##  Residual               0.9935731 0.99678                    
## Number of obs: 426964, groups:  cntry, 33
## 
## Fixed effects:
##                          Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)            -0.0400356  0.1073852 30.0881088  -0.373  0.71189    
## gndr.c                  0.5081046  0.0368039 10.5010300  13.806 4.56e-08 ***
## year.c                  0.0054152  0.0065334 30.8045947   0.829  0.41357    
## gndr.c:year.c          -0.0069019  0.0021050 12.1203255  -3.279  0.00651 ** 
## gndr.c:gei.z.cm         0.0397889  0.0236767 25.2736728   1.681  0.10518    
## year.c:gei.z.cm        -0.0178981  0.0034684 33.4983633  -5.160 1.11e-05 ***
## gndr.c:year.c:gei.z.cm  0.0008915  0.0015397 29.3304715   0.579  0.56699    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c year.c gnd.:. gn.:.. yr.:..
## gndr.c      -0.767                                   
## year.c      -0.857  0.712                            
## gndr.c:yr.c  0.753 -0.888 -0.783                     
## gndr.c:g.z. -0.001 -0.025  0.002  0.037              
## yr.c:g.z.cm  0.002  0.000 -0.007  0.002  0.162       
## gndr.c:.:..  0.001  0.037  0.000 -0.079 -0.698 -0.370
```

``` r
getFE(mod6_GEI,round=2,p.round = 5)
```

```
##                         Est.   SE    df     t       p    LL    UL
## (Intercept)            -0.04 0.11 30.09 -0.37 0.71189 -0.26  0.18
## gndr.c                  0.51 0.04 10.50 13.81 0.00000  0.43  0.59
## year.c                  0.01 0.01 30.80  0.83 0.41357 -0.01  0.02
## gndr.c:year.c          -0.01 0.00 12.12 -3.28 0.00651 -0.01  0.00
## gndr.c:gei.z.cm         0.04 0.02 25.27  1.68 0.10518 -0.01  0.09
## year.c:gei.z.cm        -0.02 0.00 33.50 -5.16 0.00001 -0.02 -0.01
## gndr.c:year.c:gei.z.cm  0.00 0.00 29.33  0.58 0.56699  0.00  0.00
```

``` r
getVC(mod6_GEI)
```

```
##         grp          var1          var2 sdcor  vcov
## 1     cntry   (Intercept)          <NA>  0.62  0.38
## 2     cntry        gndr.c          <NA>  0.21  0.04
## 3     cntry        year.c          <NA>  0.04  0.00
## 4     cntry gndr.c:year.c          <NA>  0.01  0.00
## 5     cntry   (Intercept)        gndr.c -0.79 -0.10
## 6     cntry   (Intercept)        year.c -0.86 -0.02
## 7     cntry   (Intercept) gndr.c:year.c  0.80  0.01
## 8     cntry        gndr.c        year.c  0.74  0.01
## 9     cntry        gndr.c gndr.c:year.c -0.89  0.00
## 10    cntry        year.c gndr.c:year.c -0.83  0.00
## 11 Residual          <NA>          <NA>  1.00  0.99
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_GEI,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] 21.25389
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_GEI,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] 7.58362
```

``` r
# GEI-levels  
gei_mean-gei_sd
```

```
## [1] 0.7936086
```

``` r
gei_mean+gei_sd
```

```
## [1] 0.9406057
```

``` r
# Simple slopes for 21 years
change_mod6_GEI<-emmeans(mod6_GEI,specs="year.c",by="gei.z.cm",
                     at=list(gei.z.cm=c(-1,0,1),
                             gndr.c=0,
                             year.c=rev(range(diff_dat$year.c))),
                     disable.pbkrtest=T,
                     lmerTest.limit = 500000,infer=c(T,T))
```

```
## NOTE: Results may be misleading due to involvement in interactions
```

``` r
change_mod6_GEI
```

```
## gei.z.cm = -1:
##  year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21  0.4495 0.1030 33.6   0.2411   0.6580   4.385  0.0001
##       0 -0.0400 0.1070 30.1  -0.2593   0.1792  -0.373  0.7119
## 
## gei.z.cm =  0:
##  year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21  0.0737 0.0714 31.6  -0.0717   0.2191   1.032  0.3097
##       0 -0.0400 0.1070 30.1  -0.2593   0.1792  -0.373  0.7119
## 
## gei.z.cm =  1:
##  year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21 -0.3022 0.1010 32.3  -0.5087  -0.0956  -2.979  0.0055
##       0 -0.0400 0.1070 30.1  -0.2593   0.1792  -0.373  0.7119
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_mod6_GEI,adjust="none",infer=c(T,T))
```

```
## gei.z.cm = -1:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0    0.490 0.156 42.3    0.175   0.8039   3.143  0.0031
## 
## gei.z.cm =  0:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0    0.114 0.137 30.8   -0.166   0.3936   0.829  0.4136
## 
## gei.z.cm =  1:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0   -0.262 0.155 41.2   -0.575   0.0506  -1.692  0.0981
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# change in gender differences over time

change_in_diff_mod6_GEI<-emmeans(mod6_GEI,specs=c("gndr.c","year.c"),by="gei.z.cm",
                             at=list(gndr.c=c(-0.5,0.5),
                                     gei.z.cm=c(-1,0,1),
                                     year.c=rev(range(diff_dat$year.c))),
                             disable.pbkrtest=T,
                             lmerTest.limit = 500000,infer=c(T,T))
change_in_diff_mod6_GEI
```

```
## gei.z.cm = -1:
##  gndr.c year.c  emmean     SE   df  lower.CL upper.CL t.ratio p.value
##    -0.5     21  0.2972 0.1090 33.5  7.46e-02   0.5198   2.715  0.0104
##     0.5     21  0.6019 0.0978 33.5  4.03e-01   0.8007   6.155 <0.0001
##    -0.5      0 -0.2742 0.1230 29.1 -5.25e-01  -0.0233  -2.235  0.0332
##     0.5      0  0.1941 0.0948 28.6  6.24e-05   0.3882   2.047  0.0499
## 
## gei.z.cm =  0:
##  gndr.c year.c  emmean     SE   df  lower.CL upper.CL t.ratio p.value
##    -0.5     21 -0.1079 0.0762 31.3 -2.63e-01   0.0475  -1.416  0.1667
##     0.5     21  0.2553 0.0677 31.4  1.17e-01   0.3933   3.769  0.0007
##    -0.5      0 -0.2941 0.1220 28.6 -5.44e-01  -0.0443  -2.409  0.0227
##     0.5      0  0.2140 0.0940 27.6  2.13e-02   0.4067   2.276  0.0308
## 
## gei.z.cm =  1:
##  gndr.c year.c  emmean     SE   df  lower.CL upper.CL t.ratio p.value
##    -0.5     21 -0.5130 0.1080 32.0 -7.33e-01  -0.2928  -4.745 <0.0001
##     0.5     21 -0.0913 0.0965 31.9 -2.88e-01   0.1052  -0.947  0.3509
##    -0.5      0 -0.3140 0.1230 29.1 -5.65e-01  -0.0633  -2.561  0.0159
##     0.5      0  0.2339 0.0947 28.5  4.01e-02   0.4277   2.470  0.0197
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod6_GEI,adjust="none",infer=c(T,T))
```

```
## gei.z.cm = -1:
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21    -0.3047 0.0325 29.7   -0.371  -0.2382  -9.360 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)   0.5714 0.1750 40.7    0.217   0.9253   3.261  0.0022
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0      0.1031 0.1510 41.8   -0.202   0.4082   0.682  0.4990
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)      0.8761 0.1650 40.5    0.544   1.2086   5.323 <0.0001
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0         0.4077 0.1390 40.4    0.126   0.6895   2.924  0.0056
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0      -0.4683 0.0443 15.9   -0.562  -0.3744 -10.579 <0.0001
## 
## gei.z.cm =  0:
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21    -0.3632 0.0205 24.2   -0.405  -0.3210 -17.751 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)   0.1862 0.1550 29.8   -0.131   0.5031   1.200  0.2395
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     -0.3219 0.1300 29.9   -0.586  -0.0574  -2.485  0.0188
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)      0.5494 0.1470 29.9    0.249   0.8500   3.732  0.0008
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0         0.0412 0.1210 28.9   -0.206   0.2881   0.342  0.7350
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0      -0.5081 0.0368 10.5   -0.590  -0.4266 -13.806 <0.0001
## 
## gei.z.cm =  1:
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21    -0.4217 0.0292 24.5   -0.482  -0.3615 -14.449 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)  -0.1990 0.1740 39.6   -0.551   0.1529  -1.143  0.2597
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     -0.7469 0.1500 40.4   -1.050  -0.4438  -4.978 <0.0001
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)      0.2226 0.1640 39.6   -0.108   0.5536   1.360  0.1816
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0        -0.3253 0.1380 38.7   -0.605  -0.0459  -2.356  0.0237
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0      -0.5479 0.0433 14.7   -0.640  -0.4556 -12.668 <0.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# custom contrast to obtain the estimates for the entire span

changes_in_diff <- list(
  diff_ESS11 = c(-1,1,0,0),
  diff_ESS1 = c(0,0,-1,1))

diff_mod6_GEI<-contrast(change_in_diff_mod6_GEI,method = changes_in_diff,adjust="none",infer=c(T,T))

# gender differences at ESS1 and ESS11
diff_mod6_GEI
```

```
## gei.z.cm = -1:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11    0.305 0.0325 29.7    0.238    0.371   9.360 <0.0001
##  diff_ESS1     0.468 0.0443 15.9    0.374    0.562  10.579 <0.0001
## 
## gei.z.cm =  0:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11    0.363 0.0205 24.2    0.321    0.405  17.751 <0.0001
##  diff_ESS1     0.508 0.0368 10.5    0.427    0.590  13.806 <0.0001
## 
## gei.z.cm =  1:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11    0.422 0.0292 24.5    0.362    0.482  14.449 <0.0001
##  diff_ESS1     0.548 0.0433 14.7    0.456    0.640  12.668 <0.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# Test for whether and how much the differences have changed
pairs(diff_mod6_GEI,infer=c(T,T))
```

```
## gei.z.cm = -1:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1   -0.164 0.0568 22.0   -0.281  -0.0459  -2.882  0.0087
## 
## gei.z.cm =  0:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1   -0.145 0.0442 12.1   -0.241  -0.0487  -3.279  0.0065
## 
## gei.z.cm =  1:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1   -0.126 0.0527 17.2   -0.237  -0.0152  -2.397  0.0282
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```


#### Country-level predictors: GGGI


``` r
mod6_GGGI<-lmer(FM.z~gndr.c+year.c+
             gndr.c:year.c+
               gggi.z.cm:gndr.c+gggi.z.cm:year.c+gggi.z.cm:gndr.c:year.c+
               (gndr.c+year.c+gndr.c:year.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_GGGI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: FM.z ~ gndr.c + year.c + gndr.c:year.c + gggi.z.cm:gndr.c + gggi.z.cm:year.c +  
##     gggi.z.cm:gndr.c:year.c + (gndr.c + year.c + gndr.c:year.c |      cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##  961799.6  961992.0 -480881.8  961763.6    323834 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.7224 -0.6131  0.0057  0.5995  8.2896 
## 
## Random effects:
##  Groups   Name          Variance  Std.Dev. Corr              
##  cntry    (Intercept)   4.090e-01 0.639565                   
##           gndr.c        1.743e-02 0.132027 -0.53             
##           year.c        1.369e-03 0.037001 -0.92  0.56       
##           gndr.c:year.c 8.493e-05 0.009216  0.59 -0.82 -0.62 
##  Residual               9.988e-01 0.999397                   
## Number of obs: 323852, groups:  cntry, 34
## 
## Fixed effects:
##                          Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)              0.007285   0.111544 25.539420   0.065  0.94844    
## gndr.c                   0.453061   0.025077 11.481242  18.067 8.65e-10 ***
## year.c                   0.001198   0.006461 26.894875   0.185  0.85434    
## gndr.c:year.c           -0.002794   0.001774 15.615055  -1.575  0.13534    
## gndr.c:gggi.z.cm         0.071804   0.023294 36.247709   3.082  0.00391 ** 
## year.c:gggi.z.cm        -0.009470   0.002851 33.825885  -3.321  0.00216 ** 
## gndr.c:year.c:gggi.z.cm -0.002503   0.001675 31.278977  -1.495  0.14506    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c year.c gnd.:. gn.:.. yr.:..
## gndr.c      -0.482                                   
## year.c      -0.915  0.507                            
## gndr.c:yr.c  0.523 -0.836 -0.556                     
## gndr.c:gg.. -0.001 -0.007  0.001  0.019              
## yr.c:ggg.z.  0.015 -0.007 -0.043  0.015  0.164       
## gndr.c:.:..  0.000  0.018  0.004 -0.047 -0.788 -0.206
```

``` r
getFE(mod6_GGGI,round=2,p.round = 5)
```

```
##                          Est.   SE    df     t       p    LL   UL
## (Intercept)              0.01 0.11 25.54  0.07 0.94844 -0.22 0.24
## gndr.c                   0.45 0.03 11.48 18.07 0.00000  0.40 0.51
## year.c                   0.00 0.01 26.89  0.19 0.85434 -0.01 0.01
## gndr.c:year.c            0.00 0.00 15.62 -1.57 0.13534 -0.01 0.00
## gndr.c:gggi.z.cm         0.07 0.02 36.25  3.08 0.00391  0.02 0.12
## year.c:gggi.z.cm        -0.01 0.00 33.83 -3.32 0.00216 -0.02 0.00
## gndr.c:year.c:gggi.z.cm  0.00 0.00 31.28 -1.49 0.14506 -0.01 0.00
```

``` r
getVC(mod6_GGGI)
```

```
##         grp          var1          var2 sdcor  vcov
## 1     cntry   (Intercept)          <NA>  0.64  0.41
## 2     cntry        gndr.c          <NA>  0.13  0.02
## 3     cntry        year.c          <NA>  0.04  0.00
## 4     cntry gndr.c:year.c          <NA>  0.01  0.00
## 5     cntry   (Intercept)        gndr.c -0.53 -0.04
## 6     cntry   (Intercept)        year.c -0.92 -0.02
## 7     cntry   (Intercept) gndr.c:year.c  0.59  0.00
## 8     cntry        gndr.c        year.c  0.56  0.00
## 9     cntry        gndr.c gndr.c:year.c -0.82  0.00
## 10    cntry        year.c gndr.c:year.c -0.62  0.00
## 11 Residual          <NA>          <NA>  1.00  1.00
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_GGGI,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] 23.02217
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_GGGI,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] 40.18326
```

``` r
# GGGI-levels  
gggi_mean-gggi_sd
```

```
## [1] 0.6847451
```

``` r
gggi_mean+gggi_sd
```

```
## [1] 0.7878957
```

``` r
# Simple slopes for 21 years
change_mod6_GGGI<-emmeans(mod6_GGGI,specs="year.c",by="gggi.z.cm",
                     at=list(gggi.z.cm=c(-1,0,1),
                             gndr.c=0,
                             year.c=rev(range(diff_dat$year.c))),
                     disable.pbkrtest=T,
                     lmerTest.limit = 500000,infer=c(T,T))
```

```
## NOTE: Results may be misleading due to involvement in interactions
```

``` r
change_mod6_GGGI
```

```
## gggi.z.cm = -1:
##  year.c   emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21  0.23129 0.0852 34.6   0.0584  0.40423   2.716  0.0102
##       0  0.00728 0.1120 25.5  -0.2222  0.23677   0.065  0.9484
## 
## gggi.z.cm =  0:
##  year.c   emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21  0.03243 0.0562 33.0  -0.0820  0.14687   0.577  0.5681
##       0  0.00728 0.1120 25.5  -0.2222  0.23677   0.065  0.9484
## 
## gggi.z.cm =  1:
##  year.c   emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21 -0.16643 0.0790 34.1  -0.3270 -0.00583  -2.106  0.0427
##       0  0.00728 0.1120 25.5  -0.2222  0.23677   0.065  0.9484
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_mod6_GGGI,adjust="none",infer=c(T,T))
```

```
## gggi.z.cm = -1:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0   0.2240 0.151 36.9  -0.0812    0.529   1.487  0.1455
## 
## gggi.z.cm =  0:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0   0.0251 0.136 26.9  -0.2533    0.304   0.185  0.8543
## 
## gggi.z.cm =  1:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0  -0.1737 0.146 34.1  -0.4703    0.123  -1.190  0.2421
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# change in gender differences over time

change_in_diff_mod6_GGGI<-emmeans(mod6_GGGI,specs=c("gndr.c","year.c"),by="gggi.z.cm",
                             at=list(gndr.c=c(-0.5,0.5),
                                     gggi.z.cm=c(-1,0,1),
                                     year.c=rev(range(diff_dat$year.c))),
                             disable.pbkrtest=T,
                             lmerTest.limit = 500000,infer=c(T,T))
change_in_diff_mod6_GGGI
```

```
## gggi.z.cm = -1:
##  gndr.c year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21  0.0437 0.0895 34.8  -0.1380   0.2255   0.488  0.6283
##     0.5     21  0.4189 0.0837 34.4   0.2489   0.5888   5.007 <0.0001
##    -0.5      0 -0.1833 0.1190 24.9  -0.4278   0.0611  -1.545  0.1350
##     0.5      0  0.1979 0.1070 24.2  -0.0223   0.4181   1.854  0.0759
## 
## gggi.z.cm =  0:
##  gndr.c year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21 -0.1648 0.0595 32.9  -0.2858  -0.0437  -2.770  0.0091
##     0.5     21  0.2296 0.0549 33.1   0.1179   0.3414   4.180  0.0002
##    -0.5      0 -0.2192 0.1180 24.4  -0.4628   0.0243  -1.856  0.0755
##     0.5      0  0.2338 0.1060 23.6   0.0147   0.4529   2.204  0.0375
## 
## gggi.z.cm =  1:
##  gndr.c year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21 -0.3732 0.0830 34.5  -0.5418  -0.2046  -4.496 <0.0001
##     0.5     21  0.0404 0.0777 33.8  -0.1176   0.1984   0.520  0.6068
##    -0.5      0 -0.2551 0.1190 24.9  -0.4996  -0.0107  -2.150  0.0415
##     0.5      0  0.2697 0.1070 24.1   0.0496   0.4898   2.528  0.0184
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod6_GGGI,adjust="none",infer=c(T,T))
```

```
## gggi.z.cm = -1:
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21   -0.37516 0.0318 28.3 -0.44029  -0.3100 -11.795 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)  0.22706 0.1630 36.3 -0.10409   0.5582   1.390  0.1729
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0    -0.15420 0.1500 36.1 -0.45821   0.1498  -1.029  0.3105
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)     0.60222 0.1540 35.6  0.28968   0.9148   3.910  0.0004
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0        0.22096 0.1420 34.5 -0.06688   0.5088   1.559  0.1281
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0     -0.38126 0.0343 23.2 -0.45229  -0.3102 -11.099 <0.0001
## 
## gggi.z.cm =  0:
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21   -0.39439 0.0213 27.2 -0.43815  -0.3506 -18.487 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)  0.05448 0.1470 26.3 -0.24724   0.3562   0.371  0.7136
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0    -0.39858 0.1340 26.1 -0.67401  -0.1231  -2.974  0.0063
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)     0.44887 0.1390 25.7  0.16356   0.7342   3.236  0.0033
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0       -0.00419 0.1260 24.9 -0.26430   0.2559  -0.033  0.9738
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0     -0.45306 0.0251 11.5 -0.50797  -0.3981 -18.067 <0.0001
## 
## gggi.z.cm =  1:
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21   -0.41362 0.0296 28.4 -0.47420  -0.3530 -13.979 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0) -0.11809 0.1580 33.4 -0.44040   0.2042  -0.745  0.4614
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0    -0.64296 0.1450 33.4 -0.93756  -0.3484  -4.438 <0.0001
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)     0.29553 0.1500 32.6 -0.00907   0.6001   1.975  0.0568
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0       -0.22933 0.1370 31.9 -0.50809   0.0494  -1.676  0.1035
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0     -0.52486 0.0341 22.1 -0.59557  -0.4542 -15.390 <0.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# custom contrast to obtain the estimates for the entire span

changes_in_diff <- list(
  diff_ESS11 = c(-1,1,0,0),
  diff_ESS1 = c(0,0,-1,1))

diff_mod6_GGGI<-contrast(change_in_diff_mod6_GGGI,method = changes_in_diff,adjust="none",infer=c(T,T))

# gender differences at ESS1 and ESS11
diff_mod6_GGGI
```

```
## gggi.z.cm = -1:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11    0.375 0.0318 28.3    0.310    0.440  11.795 <0.0001
##  diff_ESS1     0.381 0.0343 23.2    0.310    0.452  11.099 <0.0001
## 
## gggi.z.cm =  0:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11    0.394 0.0213 27.2    0.351    0.438  18.487 <0.0001
##  diff_ESS1     0.453 0.0251 11.5    0.398    0.508  18.067 <0.0001
## 
## gggi.z.cm =  1:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11    0.414 0.0296 28.4    0.353    0.474  13.979 <0.0001
##  diff_ESS1     0.525 0.0341 22.1    0.454    0.596  15.390 <0.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# Test for whether and how much the differences have changed
pairs(diff_mod6_GGGI,infer=c(T,T))
```

```
## gggi.z.cm = -1:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1  -0.0061 0.0524 27.8   -0.114  0.10134  -0.116  0.9082
## 
## gggi.z.cm =  0:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1  -0.0587 0.0373 15.6   -0.138  0.02047  -1.575  0.1353
## 
## gggi.z.cm =  1:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1  -0.1112 0.0500 25.6   -0.214 -0.00835  -2.224  0.0352
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

#### Country-level predictors: GDI


``` r
mod6_GDI<-lmer(FM.z~gndr.c+year.c+
             gndr.c:year.c+
               gdi.z.cm:gndr.c+gdi.z.cm:year.c+gdi.z.cm:gndr.c:year.c+
               (gndr.c+year.c+gndr.c:year.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_GDI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: FM.z ~ gndr.c + year.c + gndr.c:year.c + gdi.z.cm:gndr.c + gdi.z.cm:year.c +  
##     gdi.z.cm:gndr.c:year.c + (gndr.c + year.c + gndr.c:year.c |      cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1301501.6 1301699.4 -650732.8 1301465.6    437725 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.6900 -0.6137  0.0025  0.5975  8.3042 
## 
## Random effects:
##  Groups   Name          Variance  Std.Dev. Corr              
##  cntry    (Intercept)   0.3684872 0.60703                    
##           gndr.c        0.0463758 0.21535  -0.79             
##           year.c        0.0016173 0.04022  -0.79  0.61       
##           gndr.c:year.c 0.0001409 0.01187   0.77 -0.88 -0.76 
##  Residual               1.0039556 1.00198                    
## Number of obs: 437743, groups:  cntry, 34
## 
## Fixed effects:
##                         Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)            -0.031808   0.104338 30.911697  -0.305   0.7625    
## gndr.c                  0.501977   0.038013 11.089167  13.205 3.96e-08 ***
## year.c                  0.005164   0.006914 31.177270   0.747   0.4607    
## gndr.c:year.c          -0.006433   0.002140 12.850942  -3.007   0.0102 *  
## gndr.c:gdi.z.cm         0.030884   0.025028 28.850384   1.234   0.2272    
## year.c:gdi.z.cm        -0.006911   0.004372 33.563104  -1.581   0.1233    
## gndr.c:year.c:gdi.z.cm  0.001492   0.001631 32.328053   0.915   0.3671    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c year.c gnd.:. gn.:.. yr.:..
## gndr.c      -0.766                                   
## year.c      -0.787  0.594                            
## gndr.c:yr.c  0.729 -0.881 -0.723                     
## gndr.c:gd..  0.000 -0.018  0.000  0.020              
## yr.c:gd.z.c  0.001  0.000 -0.004  0.002 -0.023       
## gndr.c:.:..  0.000  0.018  0.001 -0.036 -0.714 -0.321
```

``` r
getFE(mod6_GDI,round=2,p.round = 5)
```

```
##                         Est.   SE    df     t       p    LL   UL
## (Intercept)            -0.03 0.10 30.91 -0.30 0.76252 -0.24 0.18
## gndr.c                  0.50 0.04 11.09 13.21 0.00000  0.42 0.59
## year.c                  0.01 0.01 31.18  0.75 0.46070 -0.01 0.02
## gndr.c:year.c          -0.01 0.00 12.85 -3.01 0.01022 -0.01 0.00
## gndr.c:gdi.z.cm         0.03 0.03 28.85  1.23 0.22718 -0.02 0.08
## year.c:gdi.z.cm        -0.01 0.00 33.56 -1.58 0.12329 -0.02 0.00
## gndr.c:year.c:gdi.z.cm  0.00 0.00 32.33  0.91 0.36710  0.00 0.00
```

``` r
getVC(mod6_GDI)
```

```
##         grp          var1          var2 sdcor  vcov
## 1     cntry   (Intercept)          <NA>  0.61  0.37
## 2     cntry        gndr.c          <NA>  0.22  0.05
## 3     cntry        year.c          <NA>  0.04  0.00
## 4     cntry gndr.c:year.c          <NA>  0.01  0.00
## 5     cntry   (Intercept)        gndr.c -0.79 -0.10
## 6     cntry   (Intercept)        year.c -0.79 -0.02
## 7     cntry   (Intercept) gndr.c:year.c  0.77  0.01
## 8     cntry        gndr.c        year.c  0.61  0.01
## 9     cntry        gndr.c gndr.c:year.c -0.88  0.00
## 10    cntry        year.c gndr.c:year.c -0.76  0.00
## 11 Residual          <NA>          <NA>  1.00  1.00
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_GDI,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] 9.067384
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_GDI,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] 0.729816
```

``` r
# GDI-levels  
gdi_mean-gdi_sd
```

```
## [1] 0.9594243
```

``` r
gdi_mean+gdi_sd
```

```
## [1] 1.010451
```

``` r
# Simple slopes for 21 years
change_mod6_GDI<-emmeans(mod6_GDI,specs="year.c",by="gdi.z.cm",
                     at=list(gdi.z.cm=c(-1,0,1),
                             gndr.c=0,
                             year.c=rev(range(diff_dat$year.c))),
                     disable.pbkrtest=T,
                     lmerTest.limit = 500000,infer=c(T,T))
```

```
## NOTE: Results may be misleading due to involvement in interactions
```

``` r
change_mod6_GDI
```

```
## gdi.z.cm = -1:
##  year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21  0.2218 0.1290 33.6  -0.0407    0.484   1.718  0.0950
##       0 -0.0318 0.1040 30.9  -0.2446    0.181  -0.305  0.7625
## 
## gdi.z.cm =  0:
##  year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21  0.0766 0.0902 32.5  -0.1070    0.260   0.850  0.4018
##       0 -0.0318 0.1040 30.9  -0.2446    0.181  -0.305  0.7625
## 
## gdi.z.cm =  1:
##  year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21 -0.0685 0.1280 32.9  -0.3296    0.193  -0.534  0.5971
##       0 -0.0318 0.1040 30.9  -0.2446    0.181  -0.305  0.7625
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_mod6_GDI,adjust="none",infer=c(T,T))
```

```
## gdi.z.cm = -1:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0   0.2536 0.172 43.0  -0.0935    0.601   1.473  0.1480
## 
## gdi.z.cm =  0:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0   0.1084 0.145 31.2  -0.1876    0.405   0.747  0.4607
## 
## gdi.z.cm =  1:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0  -0.0367 0.171 42.2  -0.3826    0.309  -0.214  0.8316
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# change in gender differences over time

change_in_diff_mod6_GDI<-emmeans(mod6_GDI,specs=c("gndr.c","year.c"),by="gdi.z.cm",
                             at=list(gndr.c=c(-0.5,0.5),
                                     gdi.z.cm=c(-1,0,1),
                                     year.c=rev(range(diff_dat$year.c))),
                             disable.pbkrtest=T,
                             lmerTest.limit = 500000,infer=c(T,T))
change_in_diff_mod6_GDI
```

```
## gdi.z.cm = -1:
##  gndr.c year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21  0.0694 0.1380 33.3  -0.2113  0.35022   0.503  0.6183
##     0.5     21  0.3741 0.1220 33.2   0.1266  0.62163   3.075  0.0042
##    -0.5      0 -0.2674 0.1200 29.9  -0.5129 -0.02184  -2.224  0.0338
##     0.5      0  0.2037 0.0915 29.2   0.0166  0.39087   2.226  0.0339
## 
## gdi.z.cm =  0:
##  gndr.c year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21 -0.1068 0.0965 32.1  -0.3033  0.08966  -1.107  0.2764
##     0.5     21  0.2601 0.0848 31.9   0.0873  0.43291   3.066  0.0044
##    -0.5      0 -0.2828 0.1200 29.3  -0.5271 -0.03847  -2.366  0.0248
##     0.5      0  0.2192 0.0906 27.9   0.0335  0.40481   2.419  0.0223
## 
## gdi.z.cm =  1:
##  gndr.c year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21 -0.2831 0.1370 32.5  -0.5623 -0.00382  -2.064  0.0471
##     0.5     21  0.1461 0.1210 32.4  -0.1001  0.39219   1.208  0.2357
##    -0.5      0 -0.2982 0.1200 29.9  -0.5436 -0.05284  -2.483  0.0189
##     0.5      0  0.2346 0.0914 29.0   0.0476  0.42160   2.566  0.0157
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod6_GDI,adjust="none",infer=c(T,T))
```

```
## gdi.z.cm = -1:
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21    -0.3047 0.0329 23.9  -0.3726  -0.2367  -9.258 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)   0.3368 0.1900 42.0  -0.0471   0.7207   1.770  0.0839
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     -0.1343 0.1720 41.6  -0.4820   0.2135  -0.780  0.4401
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)      0.6415 0.1760 42.1   0.2855   0.9974   3.637  0.0007
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0         0.1704 0.1570 41.1  -0.1471   0.4878   1.084  0.2848
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0      -0.4711 0.0459 17.1  -0.5679  -0.3743 -10.266 <0.0001
## 
## gdi.z.cm =  0:
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21    -0.3669 0.0213 21.0  -0.4113  -0.3225 -17.187 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)   0.1760 0.1620 30.3  -0.1551   0.5071   1.085  0.2864
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     -0.3260 0.1400 30.6  -0.6123  -0.0397  -2.324  0.0270
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)      0.5429 0.1530 30.4   0.2308   0.8549   3.551  0.0013
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0         0.0409 0.1300 29.8  -0.2244   0.3063   0.315  0.7550
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0      -0.5020 0.0380 11.1  -0.5856  -0.4184 -13.205 <0.0001
## 
## gdi.z.cm =  1:
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21    -0.4291 0.0313 22.8  -0.4938  -0.3644 -13.722 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)   0.0152 0.1890 41.1  -0.3673   0.3976   0.080  0.9365
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     -0.5177 0.1710 40.7  -0.8639  -0.1714  -3.020  0.0044
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)      0.4443 0.1760 41.3   0.0895   0.7991   2.528  0.0154
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0        -0.0886 0.1560 40.1  -0.4046   0.2274  -0.566  0.5743
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0      -0.5329 0.0451 16.5  -0.6283  -0.4374 -11.807 <0.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# custom contrast to obtain the estimates for the entire span

changes_in_diff <- list(
  diff_ESS11 = c(-1,1,0,0),
  diff_ESS1 = c(0,0,-1,1))

diff_mod6_GDI<-contrast(change_in_diff_mod6_GDI,method = changes_in_diff,adjust="none",infer=c(T,T))

# gender differences at ESS1 and ESS11
diff_mod6_GDI
```

```
## gdi.z.cm = -1:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11    0.305 0.0329 23.9    0.237    0.373   9.258 <0.0001
##  diff_ESS1     0.471 0.0459 17.1    0.374    0.568  10.266 <0.0001
## 
## gdi.z.cm =  0:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11    0.367 0.0213 21.0    0.322    0.411  17.187 <0.0001
##  diff_ESS1     0.502 0.0380 11.1    0.418    0.586  13.205 <0.0001
## 
## gdi.z.cm =  1:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11    0.429 0.0313 22.8    0.364    0.494  13.722 <0.0001
##  diff_ESS1     0.533 0.0451 16.5    0.437    0.628  11.807 <0.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# Test for whether and how much the differences have changed
pairs(diff_mod6_GDI,infer=c(T,T))
```

```
## gdi.z.cm = -1:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1   -0.166 0.0575 22.4   -0.286  -0.0473  -2.895  0.0083
## 
## gdi.z.cm =  0:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1   -0.135 0.0449 12.8   -0.232  -0.0379  -3.007  0.0102
## 
## gdi.z.cm =  1:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1   -0.104 0.0555 20.6   -0.219   0.0118  -1.870  0.0758
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

#### Country-level predictors: log_GDP


``` r
mod6_log_GDP<-lmer(FM.z~gndr.c+year.c+
             gndr.c:year.c+
               log_gdp.z.cm:gndr.c+log_gdp.z.cm:year.c+log_gdp.z.cm:gndr.c:year.c+
               (gndr.c+year.c+gndr.c:year.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_log_GDP)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: FM.z ~ gndr.c + year.c + gndr.c:year.c + log_gdp.z.cm:gndr.c +  
##     log_gdp.z.cm:year.c + log_gdp.z.cm:gndr.c:year.c + (gndr.c +      year.c + gndr.c:year.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##   1301498   1301696   -650731   1301462    437725 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.6980 -0.6136  0.0025  0.5976  8.3025 
## 
## Random effects:
##  Groups   Name          Variance  Std.Dev. Corr              
##  cntry    (Intercept)   0.3684814 0.60703                    
##           gndr.c        0.0434460 0.20844  -0.78             
##           year.c        0.0014388 0.03793  -0.79  0.62       
##           gndr.c:year.c 0.0001453 0.01206   0.79 -0.87 -0.80 
##  Residual               1.0039600 1.00198                    
## Number of obs: 437743, groups:  cntry, 34
## 
## Fixed effects:
##                             Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)                -0.031997   0.104337 30.834748  -0.307  0.76116    
## gndr.c                      0.498614   0.036879 11.228614  13.520 2.69e-08 ***
## year.c                      0.004877   0.006524 30.566071   0.748  0.46042    
## gndr.c:year.c              -0.006121   0.002170 13.605170  -2.820  0.01394 *  
## gndr.c:log_gdp.z.cm         0.029704   0.025356 26.550536   1.171  0.25183    
## year.c:log_gdp.z.cm        -0.011343   0.003996 31.807255  -2.838  0.00783 ** 
## gndr.c:year.c:log_gdp.z.cm -0.001273   0.001505 27.392166  -0.846  0.40503    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c year.c gnd.:. g.:_.. y.:_..
## gndr.c      -0.758                                   
## year.c      -0.792  0.604                            
## gndr.c:yr.c  0.749 -0.867 -0.765                     
## gndr.c:l_.. -0.002 -0.048  0.002  0.052              
## yr.c:lg_g..  0.001  0.000  0.012 -0.006  0.008       
## gndr.:.:_..  0.002  0.051 -0.007 -0.056 -0.716 -0.407
```

``` r
getFE(mod6_log_GDP,round=2,p.round = 5)
```

```
##                             Est.   SE    df     t       p    LL   UL
## (Intercept)                -0.03 0.10 30.83 -0.31 0.76116 -0.24 0.18
## gndr.c                      0.50 0.04 11.23 13.52 0.00000  0.42 0.58
## year.c                      0.00 0.01 30.57  0.75 0.46042 -0.01 0.02
## gndr.c:year.c              -0.01 0.00 13.61 -2.82 0.01394 -0.01 0.00
## gndr.c:log_gdp.z.cm         0.03 0.03 26.55  1.17 0.25183 -0.02 0.08
## year.c:log_gdp.z.cm        -0.01 0.00 31.81 -2.84 0.00783 -0.02 0.00
## gndr.c:year.c:log_gdp.z.cm  0.00 0.00 27.39 -0.85 0.40503  0.00 0.00
```

``` r
getVC(mod6_log_GDP)
```

```
##         grp          var1          var2 sdcor  vcov
## 1     cntry   (Intercept)          <NA>  0.61  0.37
## 2     cntry        gndr.c          <NA>  0.21  0.04
## 3     cntry        year.c          <NA>  0.04  0.00
## 4     cntry gndr.c:year.c          <NA>  0.01  0.00
## 5     cntry   (Intercept)        gndr.c -0.78 -0.10
## 6     cntry   (Intercept)        year.c -0.79 -0.02
## 7     cntry   (Intercept) gndr.c:year.c  0.79  0.01
## 8     cntry        gndr.c        year.c  0.62  0.00
## 9     cntry        gndr.c gndr.c:year.c -0.87  0.00
## 10    cntry        year.c gndr.c:year.c -0.80  0.00
## 11 Residual          <NA>          <NA>  1.00  1.00
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_log_GDP,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] 19.10258
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_log_GDP,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] -2.361865
```

``` r
# log_GDP-levels  
log_gdp_mean-log_gdp_sd
```

```
## [1] 10.19598
```

``` r
log_gdp_mean+log_gdp_sd
```

```
## [1] 11.0245
```

``` r
# Simple slopes for 21 years
change_mod6_log_GDP<-emmeans(mod6_log_GDP,specs="year.c",by="log_gdp.z.cm",
                     at=list(log_gdp.z.cm=c(-1,0,1),
                             gndr.c=0,
                             year.c=rev(range(diff_dat$year.c))),
                     disable.pbkrtest=T,
                     lmerTest.limit = 500000,infer=c(T,T))
```

```
## NOTE: Results may be misleading due to involvement in interactions
```

``` r
change_mod6_log_GDP
```

```
## log_gdp.z.cm = -1:
##  year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21  0.3086 0.1170 32.0   0.0697   0.5476   2.631  0.0130
##       0 -0.0320 0.1040 30.8  -0.2448   0.1808  -0.307  0.7612
## 
## log_gdp.z.cm =  0:
##  year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21  0.0704 0.0837 31.7  -0.1002   0.2410   0.841  0.4065
##       0 -0.0320 0.1040 30.8  -0.2448   0.1808  -0.307  0.7612
## 
## log_gdp.z.cm =  1:
##  year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21 -0.1678 0.1200 31.6  -0.4118   0.0762  -1.401  0.1708
##       0 -0.0320 0.1040 30.8  -0.2448   0.1808  -0.307  0.7612
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_mod6_log_GDP,adjust="none",infer=c(T,T))
```

```
## log_gdp.z.cm = -1:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0    0.341 0.160 41.2   0.0179    0.663   2.131  0.0391
## 
## log_gdp.z.cm =  0:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0    0.102 0.137 30.6  -0.1771    0.382   0.748  0.4604
## 
## log_gdp.z.cm =  1:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0   -0.136 0.161 41.0  -0.4619    0.190  -0.841  0.4053
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# change in gender differences over time

change_in_diff_mod6_log_GDP<-emmeans(mod6_log_GDP,specs=c("gndr.c","year.c"),by="log_gdp.z.cm",
                             at=list(gndr.c=c(-0.5,0.5),
                                     log_gdp.z.cm=c(-1,0,1),
                                     year.c=rev(range(diff_dat$year.c))),
                             disable.pbkrtest=T,
                             lmerTest.limit = 500000,infer=c(T,T))
change_in_diff_mod6_log_GDP
```

```
## log_gdp.z.cm = -1:
##  gndr.c year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21  0.1251 0.1270 31.6  -0.1337   0.3839   0.985  0.3320
##     0.5     21  0.4922 0.1090 31.4   0.2696   0.7147   4.508 <0.0001
##    -0.5      0 -0.2665 0.1200 30.0  -0.5108  -0.0221  -2.227  0.0336
##     0.5      0  0.2025 0.0922 29.7   0.0141   0.3908   2.196  0.0360
## 
## log_gdp.z.cm =  0:
##  gndr.c year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21 -0.1146 0.0907 31.2  -0.2995   0.0702  -1.264  0.2155
##     0.5     21  0.2555 0.0778 31.1   0.0968   0.4141   3.283  0.0025
##    -0.5      0 -0.2813 0.1190 29.3  -0.5244  -0.0382  -2.366  0.0248
##     0.5      0  0.2173 0.0912 28.2   0.0306   0.4040   2.384  0.0241
## 
## log_gdp.z.cm =  1:
##  gndr.c year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21 -0.3543 0.1300 31.1  -0.6185  -0.0902  -2.735  0.0102
##     0.5     21  0.0187 0.1110 31.0  -0.2083   0.2458   0.168  0.8675
##    -0.5      0 -0.2962 0.1200 29.8  -0.5403  -0.0520  -2.478  0.0191
##     0.5      0  0.2322 0.0919 29.3   0.0443   0.4200   2.526  0.0172
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod6_log_GDP,adjust="none",infer=c(T,T))
```

```
## log_gdp.z.cm = -1:
##  contrast                                     estimate     SE   df  lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21    -0.3671 0.0321 21.6 -0.433686  -0.3005 -11.445 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)   0.3915 0.1790 40.3  0.029188   0.7539   2.183  0.0349
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     -0.0774 0.1610 39.9 -0.403184   0.2484  -0.480  0.6338
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)      0.7586 0.1630 40.4  0.429027   1.0882   4.650 <0.0001
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0         0.2897 0.1430 39.0 -0.000245   0.5797   2.021  0.0502
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0      -0.4689 0.0458 18.1 -0.564987  -0.3728 -10.249 <0.0001
## 
## log_gdp.z.cm =  0:
##  contrast                                     estimate     SE   df  lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21    -0.3701 0.0229 20.7 -0.417658  -0.3225 -16.188 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)   0.1667 0.1550 29.8 -0.150203   0.4836   1.075  0.2912
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     -0.3319 0.1340 30.0 -0.604942  -0.0589  -2.483  0.0188
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)      0.5368 0.1430 29.9  0.244229   0.8293   3.748  0.0008
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0         0.0381 0.1200 28.9 -0.208240   0.2845   0.317  0.7538
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0      -0.4986 0.0369 11.2 -0.579584  -0.4176 -13.520 <0.0001
## 
## log_gdp.z.cm =  1:
##  contrast                                     estimate     SE   df  lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21    -0.3731 0.0317 19.7 -0.439270  -0.3068 -11.764 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)  -0.0582 0.1810 40.1 -0.423704   0.3074  -0.322  0.7495
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     -0.5865 0.1630 39.4 -0.916104  -0.2569  -3.598  0.0009
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)      0.3149 0.1640 40.3 -0.017367   0.6471   1.915  0.0626
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0        -0.2134 0.1450 38.3 -0.505964   0.0791  -1.477  0.1480
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0      -0.5283 0.0437 15.9 -0.621096  -0.4355 -12.079 <0.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# custom contrast to obtain the estimates for the entire span

changes_in_diff <- list(
  diff_ESS11 = c(-1,1,0,0),
  diff_ESS1 = c(0,0,-1,1))

diff_mod6_log_GDP<-contrast(change_in_diff_mod6_log_GDP,method = changes_in_diff,adjust="none",infer=c(T,T))

# gender differences at ESS1 and ESS11
diff_mod6_log_GDP
```

```
## log_gdp.z.cm = -1:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11    0.367 0.0321 21.6    0.301    0.434  11.445 <0.0001
##  diff_ESS1     0.469 0.0458 18.1    0.373    0.565  10.249 <0.0001
## 
## log_gdp.z.cm =  0:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11    0.370 0.0229 20.7    0.322    0.418  16.188 <0.0001
##  diff_ESS1     0.499 0.0369 11.2    0.418    0.580  13.520 <0.0001
## 
## log_gdp.z.cm =  1:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11    0.373 0.0317 19.7    0.307    0.439  11.764 <0.0001
##  diff_ESS1     0.528 0.0437 15.9    0.436    0.621  12.079 <0.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# Test for whether and how much the differences have changed
pairs(diff_mod6_log_GDP,infer=c(T,T))
```

```
## log_gdp.z.cm = -1:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1   -0.102 0.0569 22.4   -0.220   0.0160  -1.790  0.0870
## 
## log_gdp.z.cm =  0:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1   -0.129 0.0456 13.6   -0.227  -0.0305  -2.820  0.0139
## 
## log_gdp.z.cm =  1:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1   -0.155 0.0540 18.9   -0.268  -0.0422  -2.876  0.0097
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```


## mod7_GEI and mod8_GEI: Countries' progress in GEI and convergence of gender gap

* First fit a model including country average of gender-equality and its interaction with gender and time and differential development for women and men (fixed and random)

* Then fit a model that adds time-specific within-country fluctuations from country's average gender-equality (fixed and random) and compare to a model without these parameters.


``` r
mod7_GEI<-lmer(FM.z~gndr.c+gei.z.cm+gndr.c:gei.z.cm+
                      year.c+
                      gndr.c:year.c+
                      (gndr.c+year.c+gndr.c:year.c|cntry),
               data=diff_dat,REML=F,weights = pspwght,
                    control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

summary(mod7_GEI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: FM.z ~ gndr.c + gei.z.cm + gndr.c:gei.z.cm + year.c + gndr.c:year.c +  
##     (gndr.c + year.c + gndr.c:year.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1263220.3 1263406.7 -631593.1 1263186.3    426947 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.7611 -0.6146  0.0031  0.5990  8.3472 
## 
## Random effects:
##  Groups   Name          Variance  Std.Dev. Corr              
##  cntry    (Intercept)   0.4632528 0.68063                    
##           gndr.c        0.0313432 0.17704  -0.73             
##           year.c        0.0018475 0.04298  -0.89  0.68       
##           gndr.c:year.c 0.0001026 0.01013   0.73 -0.85 -0.78 
##  Residual               0.9935763 0.99678                    
## Number of obs: 426964, groups:  cntry, 33
## 
## Fixed effects:
##                  Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)     -0.040342   0.118698 27.577285  -0.340 0.736528    
## gndr.c           0.502423   0.032045 10.799623  15.679 8.98e-09 ***
## gei.z.cm        -0.221972   0.056050 33.109521  -3.960 0.000375 ***
## year.c           0.005376   0.007499 31.361622   0.717 0.478700    
## gndr.c:gei.z.cm  0.060200   0.018197 27.817686   3.308 0.002599 ** 
## gndr.c:year.c   -0.006764   0.001880 12.570120  -3.598 0.003411 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c g.z.cm year.c gn.:..
## gndr.c      -0.702                            
## gei.z.cm     0.001 -0.005                     
## year.c      -0.886  0.652 -0.004              
## gndr.c:g.z.  0.000  0.003 -0.389  0.002       
## gndr.c:yr.c  0.684 -0.845  0.012 -0.729 -0.032
```

``` r
getFE(mod7_GEI)
```

```
##                  Est.   SE    df     t     p    LL    UL
## (Intercept)     -0.04 0.12 27.58 -0.34 0.737 -0.28  0.20
## gndr.c           0.50 0.03 10.80 15.68 0.000  0.43  0.57
## gei.z.cm        -0.22 0.06 33.11 -3.96 0.000 -0.34 -0.11
## year.c           0.01 0.01 31.36  0.72 0.479 -0.01  0.02
## gndr.c:gei.z.cm  0.06 0.02 27.82  3.31 0.003  0.02  0.10
## gndr.c:year.c   -0.01 0.00 12.57 -3.60 0.003 -0.01  0.00
```

``` r
getVC(mod7_GEI)
```

```
##         grp          var1          var2 sdcor  vcov
## 1     cntry   (Intercept)          <NA>  0.68  0.46
## 2     cntry        gndr.c          <NA>  0.18  0.03
## 3     cntry        year.c          <NA>  0.04  0.00
## 4     cntry gndr.c:year.c          <NA>  0.01  0.00
## 5     cntry   (Intercept)        gndr.c -0.73 -0.09
## 6     cntry   (Intercept)        year.c -0.89 -0.03
## 7     cntry   (Intercept) gndr.c:year.c  0.73  0.01
## 8     cntry        gndr.c        year.c  0.68  0.01
## 9     cntry        gndr.c gndr.c:year.c -0.85  0.00
## 10    cntry        year.c gndr.c:year.c -0.78  0.00
## 11 Residual          <NA>          <NA>  1.00  0.99
```

``` r
anova(mod2_GEI,mod7_GEI)
```

```
## Data: diff_dat
## Models:
## mod2_GEI: FM.z ~ gndr.c + gei.z.cm + gndr.c:gei.z.cm + (gndr.c | cntry)
## mod7_GEI: FM.z ~ gndr.c + gei.z.cm + gndr.c:gei.z.cm + year.c + gndr.c:year.c + (gndr.c + year.c + gndr.c:year.c | cntry)
##          npar     AIC     BIC  logLik -2*log(L) Chisq Df Pr(>Chisq)    
## mod2_GEI    8 1268297 1268385 -634141   1268281                        
## mod7_GEI   17 1263220 1263407 -631593   1263186  5095  9  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

### Examine trends


``` r
change_mod7_GEI<-emmeans(mod7_GEI,specs="year.c",by="gndr.c",
                              at=list(gndr.c=c(-0.5,0.5),
                                      year.c=rev(range(diff_dat$year.c)),
                                      gei.z.cmc=0,
                                      gei.z.cm=0),
                              disable.pbkrtest=T,
                              lmerTest.limit = 500000,infer=c(T,T),adjust="none")
change_mod7_GEI
```

```
## gndr.c = -0.5:
##  year.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21 -0.108 0.0811 26.5  -0.2743   0.0590  -1.326  0.1960
##       0 -0.292 0.1300 27.3  -0.5590  -0.0241  -2.235  0.0338
## 
## gndr.c =  0.5:
##  year.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21  0.253 0.0719 25.6   0.1048   0.4007   3.513  0.0017
##       0  0.211 0.1080 25.4  -0.0115   0.4332   1.951  0.0621
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_mod7_GEI,adjust="none",infer=c(T,T))
```

```
## gndr.c = -0.5:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0   0.1839 0.172 30.6   -0.168    0.536   1.067  0.2943
## 
## gndr.c =  0.5:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0   0.0419 0.144 30.2   -0.252    0.335   0.291  0.7727
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# change in gender differences over time

change_in_diff_mod7_GEI<-emmeans(mod7_GEI,specs=c("gndr.c","year.c"),
                                      at=list(gndr.c=c(-0.5,0.5),
                                              year.c=rev(range(diff_dat$year.c)),
                                              gei.z.cmc=0,
                                              gei.z.cm=0),
                                      disable.pbkrtest=T,
                                      lmerTest.limit = 500000,infer=c(T,T),adjust="none")
change_in_diff_mod7_GEI
```

```
##  gndr.c year.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21 -0.108 0.0811 26.5  -0.2743   0.0590  -1.326  0.1960
##     0.5     21  0.253 0.0719 25.6   0.1048   0.4007   3.513  0.0017
##    -0.5      0 -0.292 0.1300 27.3  -0.5590  -0.0241  -2.235  0.0338
##     0.5      0  0.211 0.1080 25.4  -0.0115   0.4332   1.951  0.0621
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod7_GEI,adjust="none",infer=c(T,T))
```

```
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21    -0.3604 0.0211 21.9   -0.404 -0.31655 -17.053 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)   0.1839 0.1720 30.6   -0.168  0.53568   1.067  0.2943
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     -0.3185 0.1520 30.6   -0.629 -0.00835  -2.096  0.0445
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)      0.5443 0.1650 30.6    0.208  0.88081   3.301  0.0025
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0         0.0419 0.1440 30.2   -0.252  0.33534   0.291  0.7727
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0      -0.5024 0.0320 10.8   -0.573 -0.43173 -15.679 <0.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# custom contrast to obtain the estimates for the entire span

changes_in_diff <- list(
  diff_ESS11 = c(-1,1,0,0),
  diff_ESS1 = c(0,0,-1,1))

diff_mod7_GEI<-contrast(change_in_diff_mod7_GEI,method = changes_in_diff,adjust="none",infer=c(T,T))

# gender differences at ESS1 and ESS10
diff_mod7_GEI
```

```
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11    0.360 0.0211 21.9    0.317    0.404  17.053 <0.0001
##  diff_ESS1     0.502 0.0320 10.8    0.432    0.573  15.679 <0.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# Test for whether and how much the differences have changed
pairs(diff_mod7_GEI,infer=c(T,T))
```

```
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1   -0.142 0.0395 12.6   -0.228  -0.0565  -3.598  0.0034
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```


``` r
mod8_GEI<-lmer(FM.z~gndr.c+gei.z.cm+gndr.c:gei.z.cm+
                      year.c+
                      gndr.c:year.c+
                      gei.z.cmc+gei.z.cmc:gndr.c+
                      (gndr.c+year.c+gndr.c:year.c+gei.z.cmc+gei.z.cmc:gndr.c|cntry),
               data=diff_dat,REML=F,weights = pspwght,
                    control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

getFE(mod8_GEI)
```

```
##                   Est.   SE    df     t     p    LL    UL
## (Intercept)       0.08 0.13 25.92  0.64 0.529 -0.19  0.36
## gndr.c            0.51 0.04 18.36 12.83 0.000  0.42  0.59
## gei.z.cm         -0.12 0.04 30.84 -2.80 0.009 -0.20 -0.03
## year.c            0.00 0.01 24.35 -0.45 0.659 -0.03  0.02
## gei.z.cmc        -0.02 0.10 24.61 -0.24 0.811 -0.23  0.19
## gndr.c:gei.z.cm   0.03 0.01 29.54  2.16 0.039  0.00  0.06
## gndr.c:year.c    -0.01 0.00 16.10 -2.51 0.023 -0.01  0.00
## gndr.c:gei.z.cmc  0.05 0.04 19.81  1.38 0.182 -0.03  0.12
```

``` r
getVC(mod8_GEI)
```

```
##         grp             var1             var2 sdcor  vcov
## 1     cntry      (Intercept)             <NA>  0.75  0.56
## 2     cntry           gndr.c             <NA>  0.20  0.04
## 3     cntry           year.c             <NA>  0.06  0.00
## 4     cntry        gei.z.cmc             <NA>  0.57  0.32
## 5     cntry    gndr.c:year.c             <NA>  0.01  0.00
## 6     cntry gndr.c:gei.z.cmc             <NA>  0.17  0.03
## 7     cntry      (Intercept)           gndr.c -0.35 -0.05
## 8     cntry      (Intercept)           year.c -0.94 -0.04
## 9     cntry      (Intercept)        gei.z.cmc  0.92  0.40
## 10    cntry      (Intercept)    gndr.c:year.c  0.48  0.01
## 11    cntry      (Intercept) gndr.c:gei.z.cmc -0.35 -0.04
## 12    cntry           gndr.c           year.c  0.33  0.00
## 13    cntry           gndr.c        gei.z.cmc -0.37 -0.04
## 14    cntry           gndr.c    gndr.c:year.c -0.94  0.00
## 15    cntry           gndr.c gndr.c:gei.z.cmc  0.91  0.03
## 16    cntry           year.c        gei.z.cmc -0.95 -0.03
## 17    cntry           year.c    gndr.c:year.c -0.50  0.00
## 18    cntry           year.c gndr.c:gei.z.cmc  0.40  0.00
## 19    cntry        gei.z.cmc    gndr.c:year.c  0.52  0.00
## 20    cntry        gei.z.cmc gndr.c:gei.z.cmc -0.46 -0.04
## 21    cntry    gndr.c:year.c gndr.c:gei.z.cmc -0.96  0.00
## 22 Residual             <NA>             <NA>  1.00  0.99
```

``` r
anova(mod2_GEI,mod7_GEI,mod8_GEI)
```

```
## Data: diff_dat
## Models:
## mod2_GEI: FM.z ~ gndr.c + gei.z.cm + gndr.c:gei.z.cm + (gndr.c | cntry)
## mod7_GEI: FM.z ~ gndr.c + gei.z.cm + gndr.c:gei.z.cm + year.c + gndr.c:year.c + (gndr.c + year.c + gndr.c:year.c | cntry)
## mod8_GEI: FM.z ~ gndr.c + gei.z.cm + gndr.c:gei.z.cm + year.c + gndr.c:year.c + gei.z.cmc + gei.z.cmc:gndr.c + (gndr.c + year.c + gndr.c:year.c + gei.z.cmc + gei.z.cmc:gndr.c | cntry)
##          npar     AIC     BIC  logLik -2*log(L)   Chisq Df Pr(>Chisq)    
## mod2_GEI    8 1268297 1268385 -634141   1268281                          
## mod7_GEI   17 1263220 1263407 -631593   1263186 5094.98  9  < 2.2e-16 ***
## mod8_GEI   30 1262891 1263220 -631416   1262831  355.23 13  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

### Examine trends


``` r
# gender specific change over time

change_mod8_GEI<-emmeans(mod8_GEI,specs="year.c",by="gndr.c",
                              at=list(gndr.c=c(-0.5,0.5),
                                      year.c=rev(range(diff_dat$year.c)),
                                      gei.z.cmc=0,
                                      gei.z.cm=0),
                              disable.pbkrtest=T,
                              lmerTest.limit = 500000,infer=c(T,T),adjust="none")
change_mod8_GEI
```

```
## gndr.c = -0.5:
##  year.c emmean    SE   df lower.CL upper.CL t.ratio p.value
##      21 -0.195 0.121 25.0  -0.4455   0.0546  -1.609  0.1201
##       0 -0.168 0.140 25.8  -0.4555   0.1188  -1.206  0.2389
## 
## gndr.c =  0.5:
##  year.c emmean    SE   df lower.CL upper.CL t.ratio p.value
##      21  0.158 0.108 24.1  -0.0642   0.3805   1.468  0.1551
##       0  0.337 0.127 25.0   0.0745   0.5992   2.644  0.0139
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_mod8_GEI,adjust="none",infer=c(T,T))
```

```
## gndr.c = -0.5:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0  -0.0271 0.245 24.8   -0.531    0.477  -0.111  0.9128
## 
## gndr.c =  0.5:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0  -0.1787 0.219 23.1   -0.632    0.274  -0.816  0.4228
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# change in gender differences over time

change_in_diff_mod8_GEI<-emmeans(mod8_GEI,specs=c("gndr.c","year.c"),
                                      at=list(gndr.c=c(-0.5,0.5),
                                              year.c=rev(range(diff_dat$year.c)),
                                              gei.z.cmc=0,
                                              gei.z.cm=0),
                                      disable.pbkrtest=T,
                                      lmerTest.limit = 500000,infer=c(T,T),adjust="none")
change_in_diff_mod8_GEI
```

```
##  gndr.c year.c emmean    SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21 -0.195 0.121 25.0  -0.4455   0.0546  -1.609  0.1201
##     0.5     21  0.158 0.108 24.1  -0.0642   0.3805   1.468  0.1551
##    -0.5      0 -0.168 0.140 25.8  -0.4555   0.1188  -1.206  0.2389
##     0.5      0  0.337 0.127 25.0   0.0745   0.5992   2.644  0.0139
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod8_GEI,adjust="none",infer=c(T,T))
```

```
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21    -0.3536 0.0272 16.2   -0.411   -0.296 -12.990 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)  -0.0271 0.2450 24.8   -0.531    0.477  -0.111  0.9128
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     -0.5323 0.2320 24.2   -1.011   -0.054  -2.296  0.0306
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)      0.3265 0.2300 24.5   -0.147    0.800   1.422  0.1676
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0        -0.1787 0.2190 23.1   -0.632    0.274  -0.816  0.4228
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0      -0.5052 0.0394 18.4   -0.588   -0.423 -12.830 <0.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# custom contrast to obtain the estimates for the entire span

changes_in_diff <- list(
  diff_ESS11 = c(-1,1,0,0),
  diff_ESS1 = c(0,0,-1,1))

diff_mod8_GEI<-contrast(change_in_diff_mod8_GEI,method = changes_in_diff,adjust="none",infer=c(T,T))

# gender differences at ESS1 and ESS10
diff_mod8_GEI
```

```
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11    0.354 0.0272 16.2    0.296    0.411  12.990 <0.0001
##  diff_ESS1     0.505 0.0394 18.4    0.423    0.588  12.830 <0.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# Test for whether and how much the differences have changed
pairs(diff_mod8_GEI,infer=c(T,T))
```

```
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1   -0.152 0.0605 16.1    -0.28  -0.0234  -2.506  0.0233
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

## mod7_GGGI and mod8_GGGI: Countries' progress in GGGI and convergence of gender gap

* First fit a model including country average of gender-equality and its interaction with gender and time and differential development for women and men (fixed and random)

* Then fit a model that adds time-specific within-country fluctuations from country's average gender-equality (fixed and random) and compare to a model without these parameters.


``` r
mod7_GGGI<-lmer(FM.z~gndr.c+gggi.z.cm+gndr.c:gggi.z.cm+
                      year.c+
                      gndr.c:year.c+
                      (gndr.c+year.c+gndr.c:year.c|cntry),
               data=diff_dat,REML=F,weights = pspwght,
                    control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

summary(mod7_GGGI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: FM.z ~ gndr.c + gggi.z.cm + gndr.c:gggi.z.cm + year.c + gndr.c:year.c +  
##     (gndr.c + year.c + gndr.c:year.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##  961804.4  961986.1 -480885.2  961770.4    323835 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.7315 -0.6131  0.0057  0.5996  8.2896 
## 
## Random effects:
##  Groups   Name          Variance  Std.Dev. Corr              
##  cntry    (Intercept)   4.438e-01 0.666160                   
##           gndr.c        1.433e-02 0.119702 -0.32             
##           year.c        1.516e-03 0.038931 -0.92  0.35       
##           gndr.c:year.c 7.386e-05 0.008594  0.40 -0.77 -0.43 
##  Residual               9.988e-01 0.999396                   
## Number of obs: 323852, groups:  cntry, 34
## 
## Fixed effects:
##                    Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)       0.0094238  0.1161894 25.6898375   0.081  0.93599    
## gndr.c            0.4487824  0.0231658 18.4016502  19.373 1.08e-13 ***
## gggi.z.cm        -0.1390368  0.0478081 35.3710265  -2.908  0.00624 ** 
## year.c            0.0008816  0.0067923 28.0966799   0.130  0.89766    
## gndr.c:gggi.z.cm  0.0471339  0.0145169 37.4966990   3.247  0.00246 ** 
## gndr.c:year.c    -0.0027156  0.0016797 19.4844697  -1.617  0.12201    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c ggg.z. year.c gn.:..
## gndr.c      -0.286                            
## gggi.z.cm    0.012 -0.004                     
## year.c      -0.922  0.310 -0.035              
## gndr.c:gg.. -0.001  0.012 -0.006  0.004       
## gndr.c:yr.c  0.347 -0.800  0.010 -0.375 -0.032
```

``` r
getFE(mod7_GGGI)
```

```
##                   Est.   SE    df     t     p    LL    UL
## (Intercept)       0.01 0.12 25.69  0.08 0.936 -0.23  0.25
## gndr.c            0.45 0.02 18.40 19.37 0.000  0.40  0.50
## gggi.z.cm        -0.14 0.05 35.37 -2.91 0.006 -0.24 -0.04
## year.c            0.00 0.01 28.10  0.13 0.898 -0.01  0.01
## gndr.c:gggi.z.cm  0.05 0.01 37.50  3.25 0.002  0.02  0.08
## gndr.c:year.c     0.00 0.00 19.48 -1.62 0.122 -0.01  0.00
```

``` r
getVC(mod7_GGGI)
```

```
##         grp          var1          var2 sdcor  vcov
## 1     cntry   (Intercept)          <NA>  0.67  0.44
## 2     cntry        gndr.c          <NA>  0.12  0.01
## 3     cntry        year.c          <NA>  0.04  0.00
## 4     cntry gndr.c:year.c          <NA>  0.01  0.00
## 5     cntry   (Intercept)        gndr.c -0.32 -0.03
## 6     cntry   (Intercept)        year.c -0.92 -0.02
## 7     cntry   (Intercept) gndr.c:year.c  0.40  0.00
## 8     cntry        gndr.c        year.c  0.35  0.00
## 9     cntry        gndr.c gndr.c:year.c -0.77  0.00
## 10    cntry        year.c gndr.c:year.c -0.43  0.00
## 11 Residual          <NA>          <NA>  1.00  1.00
```

``` r
anova(mod2_GGGI,mod7_GGGI)
```

```
## Data: diff_dat
## Models:
## mod2_GGGI: FM.z ~ gndr.c + gggi.z.cm + gndr.c:gggi.z.cm + (gndr.c | cntry)
## mod7_GGGI: FM.z ~ gndr.c + gggi.z.cm + gndr.c:gggi.z.cm + year.c + gndr.c:year.c + (gndr.c + year.c + gndr.c:year.c | cntry)
##           npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod2_GGGI    8 964694 964779 -482339    964678                         
## mod7_GGGI   17 961804 961986 -480885    961770 2907.1  9  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

### Examine trends


``` r
change_mod7_GGGI<-emmeans(mod7_GGGI,specs="year.c",by="gndr.c",
                              at=list(gndr.c=c(-0.5,0.5),
                                      year.c=rev(range(diff_dat$year.c)),
                                      gggi.z.cmc=0,
                                      gggi.z.cm=0),
                              disable.pbkrtest=T,
                              lmerTest.limit = 500000,infer=c(T,T),adjust="none")
change_mod7_GGGI
```

```
## gndr.c = -0.5:
##  year.c emmean     SE   df  lower.CL upper.CL t.ratio p.value
##      21 -0.168 0.0600 32.0 -0.290226  -0.0457  -2.797  0.0086
##       0 -0.215 0.1200 24.7 -0.462280   0.0323  -1.791  0.0855
## 
## gndr.c =  0.5:
##  year.c emmean     SE   df  lower.CL upper.CL t.ratio p.value
##      21  0.224 0.0565 30.8  0.108637   0.3390   3.964  0.0004
##       0  0.234 0.1130 23.7 -0.000415   0.4680   2.061  0.0504
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_mod7_GGGI,adjust="none",infer=c(T,T))
```

```
## gndr.c = -0.5:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0    0.047 0.150 27.1   -0.261    0.355   0.313  0.7565
## 
## gndr.c =  0.5:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0   -0.010 0.137 26.7   -0.291    0.271  -0.073  0.9424
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# change in gender differences over time

change_in_diff_mod7_GGGI<-emmeans(mod7_GGGI,specs=c("gndr.c","year.c"),
                                      at=list(gndr.c=c(-0.5,0.5),
                                              year.c=rev(range(diff_dat$year.c)),
                                              gggi.z.cmc=0,
                                              gggi.z.cm=0),
                                      disable.pbkrtest=T,
                                      lmerTest.limit = 500000,infer=c(T,T),adjust="none")
change_in_diff_mod7_GGGI
```

```
##  gndr.c year.c emmean     SE   df  lower.CL upper.CL t.ratio p.value
##    -0.5     21 -0.168 0.0600 32.0 -0.290226  -0.0457  -2.797  0.0086
##     0.5     21  0.224 0.0565 30.8  0.108637   0.3390   3.964  0.0004
##    -0.5      0 -0.215 0.1200 24.7 -0.462280   0.0323  -1.791  0.0855
##     0.5      0  0.234 0.1130 23.7 -0.000415   0.4680   2.061  0.0504
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod7_GGGI,adjust="none",infer=c(T,T))
```

```
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21     -0.392 0.0217 26.2   -0.436   -0.347 -18.020 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)    0.047 0.1500 27.1   -0.261    0.355   0.313  0.7565
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0      -0.402 0.1430 27.4   -0.694   -0.109  -2.815  0.0089
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)       0.439 0.1440 26.7    0.143    0.734   3.049  0.0051
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0         -0.010 0.1370 26.7   -0.291    0.271  -0.073  0.9424
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0       -0.449 0.0232 18.4   -0.497   -0.400 -19.373 <0.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# custom contrast to obtain the estimates for the entire span

changes_in_diff <- list(
  diff_ESS11 = c(-1,1,0,0),
  diff_ESS1 = c(0,0,-1,1))

diff_mod7_GGGI<-contrast(change_in_diff_mod7_GGGI,method = changes_in_diff,adjust="none",infer=c(T,T))

# gender differences at ESS1 and ESS10
diff_mod7_GGGI
```

```
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11    0.392 0.0217 26.2    0.347    0.436  18.020 <0.0001
##  diff_ESS1     0.449 0.0232 18.4    0.400    0.497  19.373 <0.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# Test for whether and how much the differences have changed
pairs(diff_mod7_GGGI,infer=c(T,T))
```

```
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1   -0.057 0.0353 19.5   -0.131   0.0167  -1.617  0.1220
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```


``` r
mod8_GGGI<-lmer(FM.z~gndr.c+gggi.z.cm+gndr.c:gggi.z.cm+
                      year.c+
                      gndr.c:year.c+
                      gggi.z.cmc+gggi.z.cmc:gndr.c+
                      (gndr.c+year.c+gndr.c:year.c+gggi.z.cmc+gggi.z.cmc:gndr.c|cntry),
               data=diff_dat,REML=F,weights = pspwght,
                    control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
```

```
## boundary (singular) fit: see help('isSingular')
```

``` r
getFE(mod8_GGGI)
```

```
##                    Est.   SE    df     t     p    LL    UL
## (Intercept)        0.05 0.07 17.16  0.70 0.492 -0.09  0.18
## gndr.c             0.46 0.03 24.63 15.19 0.000  0.40  0.52
## gggi.z.cm         -0.10 0.04 34.43 -2.44 0.020 -0.18 -0.02
## year.c             0.00 0.01 27.92  0.09 0.931 -0.01  0.01
## gggi.z.cmc        -0.06 0.05 21.38 -1.33 0.196 -0.16  0.04
## gndr.c:gggi.z.cm   0.03 0.01 38.85  1.92 0.062  0.00  0.06
## gndr.c:year.c      0.00 0.00 25.70 -1.90 0.069 -0.01  0.00
## gndr.c:gggi.z.cmc  0.03 0.02 21.42  1.76 0.093 -0.01  0.07
```

``` r
getVC(mod8_GGGI)
```

```
##         grp              var1              var2 sdcor  vcov
## 1     cntry       (Intercept)              <NA>  0.37  0.13
## 2     cntry            gndr.c              <NA>  0.15  0.02
## 3     cntry            year.c              <NA>  0.03  0.00
## 4     cntry        gggi.z.cmc              <NA>  0.26  0.07
## 5     cntry     gndr.c:year.c              <NA>  0.01  0.00
## 6     cntry gndr.c:gggi.z.cmc              <NA>  0.07  0.00
## 7     cntry       (Intercept)            gndr.c  0.30  0.02
## 8     cntry       (Intercept)            year.c -0.77 -0.01
## 9     cntry       (Intercept)        gggi.z.cmc  0.67  0.06
## 10    cntry       (Intercept)     gndr.c:year.c -0.08  0.00
## 11    cntry       (Intercept) gndr.c:gggi.z.cmc  0.04  0.00
## 12    cntry            gndr.c            year.c -0.24  0.00
## 13    cntry            gndr.c        gggi.z.cmc  0.37  0.01
## 14    cntry            gndr.c     gndr.c:year.c -0.83  0.00
## 15    cntry            gndr.c gndr.c:gggi.z.cmc  0.73  0.01
## 16    cntry            year.c        gggi.z.cmc -0.83 -0.01
## 17    cntry            year.c     gndr.c:year.c -0.02  0.00
## 18    cntry            year.c gndr.c:gggi.z.cmc  0.27  0.00
## 19    cntry        gggi.z.cmc     gndr.c:year.c -0.15  0.00
## 20    cntry        gggi.z.cmc gndr.c:gggi.z.cmc -0.28  0.00
## 21    cntry     gndr.c:year.c gndr.c:gggi.z.cmc -0.78  0.00
## 22 Residual              <NA>              <NA>  1.00  1.00
```

``` r
anova(mod2_GGGI,mod7_GGGI,mod8_GGGI)
```

```
## Data: diff_dat
## Models:
## mod2_GGGI: FM.z ~ gndr.c + gggi.z.cm + gndr.c:gggi.z.cm + (gndr.c | cntry)
## mod7_GGGI: FM.z ~ gndr.c + gggi.z.cm + gndr.c:gggi.z.cm + year.c + gndr.c:year.c + (gndr.c + year.c + gndr.c:year.c | cntry)
## mod8_GGGI: FM.z ~ gndr.c + gggi.z.cm + gndr.c:gggi.z.cm + year.c + gndr.c:year.c + gggi.z.cmc + gggi.z.cmc:gndr.c + (gndr.c + year.c + gndr.c:year.c + gggi.z.cmc + gggi.z.cmc:gndr.c | cntry)
##           npar    AIC    BIC  logLik -2*log(L)   Chisq Df Pr(>Chisq)    
## mod2_GGGI    8 964694 964779 -482339    964678                          
## mod7_GGGI   17 961804 961986 -480885    961770 2907.11  9  < 2.2e-16 ***
## mod8_GGGI   30 961583 961903 -480761    961523  247.85 13  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

### Examine trends


``` r
# gender specific change over time

change_mod8_GGGI<-emmeans(mod8_GGGI,specs="year.c",by="gndr.c",
                              at=list(gndr.c=c(-0.5,0.5),
                                      year.c=rev(range(diff_dat$year.c)),
                                      gggi.z.cmc=0,
                                      gggi.z.cm=0),
                              disable.pbkrtest=T,
                              lmerTest.limit = 500000,infer=c(T,T),adjust="none")
change_mod8_GGGI
```

```
## gndr.c = -0.5:
##  year.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21 -0.132 0.0724 26.4   -0.281   0.0163  -1.829  0.0787
##       0 -0.183 0.0630 15.3   -0.317  -0.0490  -2.906  0.0107
## 
## gndr.c =  0.5:
##  year.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21  0.242 0.0664 25.4    0.105   0.3788   3.645  0.0012
##       0  0.274 0.0702 19.0    0.127   0.4214   3.906  0.0009
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_mod8_GGGI,adjust="none",infer=c(T,T))
```

```
## gndr.c = -0.5:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0   0.0508 0.109 26.9   -0.173    0.275   0.466  0.6451
## 
## gndr.c =  0.5:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0  -0.0323 0.108 27.9   -0.254    0.190  -0.298  0.7680
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# change in gender differences over time

change_in_diff_mod8_GGGI<-emmeans(mod8_GGGI,specs=c("gndr.c","year.c"),
                                      at=list(gndr.c=c(-0.5,0.5),
                                              year.c=rev(range(diff_dat$year.c)),
                                              gggi.z.cmc=0,
                                              gggi.z.cm=0),
                                      disable.pbkrtest=T,
                                      lmerTest.limit = 500000,infer=c(T,T),adjust="none")
change_in_diff_mod8_GGGI
```

```
##  gndr.c year.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21 -0.132 0.0724 26.4   -0.281   0.0163  -1.829  0.0787
##     0.5     21  0.242 0.0664 25.4    0.105   0.3788   3.645  0.0012
##    -0.5      0 -0.183 0.0630 15.3   -0.317  -0.0490  -2.906  0.0107
##     0.5      0  0.274 0.0702 19.0    0.127   0.4214   3.906  0.0009
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod8_GGGI,adjust="none",infer=c(T,T))
```

```
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21    -0.3745 0.0233 26.7   -0.422   -0.327 -16.086 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)   0.0508 0.1090 26.9   -0.173    0.275   0.466  0.6451
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     -0.4068 0.1140 27.9   -0.640   -0.173  -3.571  0.0013
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)      0.4253 0.1010 25.2    0.217    0.633   4.207  0.0003
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0        -0.0323 0.1080 27.9   -0.254    0.190  -0.298  0.7680
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0      -0.4576 0.0301 24.6   -0.520   -0.395 -15.192 <0.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# custom contrast to obtain the estimates for the entire span

changes_in_diff <- list(
  diff_ESS11 = c(-1,1,0,0),
  diff_ESS1 = c(0,0,-1,1))

diff_mod8_GGGI<-contrast(change_in_diff_mod8_GGGI,method = changes_in_diff,adjust="none",infer=c(T,T))

# gender differences at ESS1 and ESS10
diff_mod8_GGGI
```

```
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11    0.374 0.0233 26.7    0.327    0.422  16.086 <0.0001
##  diff_ESS1     0.458 0.0301 24.6    0.395    0.520  15.192 <0.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# Test for whether and how much the differences have changed
pairs(diff_mod8_GGGI,infer=c(T,T))
```

```
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1  -0.0831 0.0438 25.7   -0.173  0.00688  -1.900  0.0688
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

## mod7_GDI and mod8_GDI: Countries' progress in GDI and convergence of gender gap

* First fit a model including country average of gender-equality and its interaction with gender and time and differential development for women and men (fixed and random)

* Then fit a model that adds time-specific within-country fluctuations from country's average gender-equality (fixed and random) and compare to a model without these parameters.


``` r
mod7_GDI<-lmer(FM.z~gndr.c+gdi.z.cm+gndr.c:gdi.z.cm+
                      year.c+
                      gndr.c:year.c+
                      (gndr.c+year.c+gndr.c:year.c|cntry),
               data=diff_dat,REML=F,weights = pspwght,
                    control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

summary(mod7_GDI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: FM.z ~ gndr.c + gdi.z.cm + gndr.c:gdi.z.cm + year.c + gndr.c:year.c +  
##     (gndr.c + year.c + gndr.c:year.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##   1301502   1301689   -650734   1301468    437726 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.6923 -0.6137  0.0026  0.5975  8.3043 
## 
## Random effects:
##  Groups   Name          Variance Std.Dev. Corr              
##  cntry    (Intercept)   0.375902 0.61311                    
##           gndr.c        0.048639 0.22054  -0.80             
##           year.c        0.001778 0.04217  -0.80  0.63       
##           gndr.c:year.c 0.000153 0.01237   0.78 -0.89 -0.78 
##  Residual               1.003955 1.00198                    
## Number of obs: 437743, groups:  cntry, 34
## 
## Fixed effects:
##                  Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)     -0.031761   0.105378 28.945048  -0.301   0.7653    
## gndr.c           0.502473   0.038878 11.786532  12.924 2.57e-08 ***
## gdi.z.cm        -0.018861   0.060327 30.129698  -0.313   0.7567    
## year.c           0.005139   0.007249 32.245049   0.709   0.4835    
## gndr.c:gdi.z.cm  0.041676   0.019411 28.447615   2.147   0.0405 *  
## gndr.c:year.c   -0.006437   0.002220 13.978359  -2.900   0.0117 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c gd.z.c year.c gn.:..
## gndr.c      -0.775                            
## gdi.z.cm    -0.001 -0.006                     
## year.c      -0.799  0.615 -0.002              
## gndr.c:gd..  0.000 -0.001 -0.550  0.001       
## gndr.c:yr.c  0.748 -0.888  0.013 -0.741 -0.018
```

``` r
getFE(mod7_GDI)
```

```
##                  Est.   SE    df     t     p    LL   UL
## (Intercept)     -0.03 0.11 28.95 -0.30 0.765 -0.25 0.18
## gndr.c           0.50 0.04 11.79 12.92 0.000  0.42 0.59
## gdi.z.cm        -0.02 0.06 30.13 -0.31 0.757 -0.14 0.10
## year.c           0.01 0.01 32.25  0.71 0.483 -0.01 0.02
## gndr.c:gdi.z.cm  0.04 0.02 28.45  2.15 0.040  0.00 0.08
## gndr.c:year.c   -0.01 0.00 13.98 -2.90 0.012 -0.01 0.00
```

``` r
getVC(mod7_GDI)
```

```
##         grp          var1          var2 sdcor  vcov
## 1     cntry   (Intercept)          <NA>  0.61  0.38
## 2     cntry        gndr.c          <NA>  0.22  0.05
## 3     cntry        year.c          <NA>  0.04  0.00
## 4     cntry gndr.c:year.c          <NA>  0.01  0.00
## 5     cntry   (Intercept)        gndr.c -0.80 -0.11
## 6     cntry   (Intercept)        year.c -0.80 -0.02
## 7     cntry   (Intercept) gndr.c:year.c  0.78  0.01
## 8     cntry        gndr.c        year.c  0.63  0.01
## 9     cntry        gndr.c gndr.c:year.c -0.89  0.00
## 10    cntry        year.c gndr.c:year.c -0.78  0.00
## 11 Residual          <NA>          <NA>  1.00  1.00
```

``` r
anova(mod2_GDI,mod7_GDI)
```

```
## Data: diff_dat
## Models:
## mod2_GDI: FM.z ~ gndr.c + gdi.z.cm + gndr.c:gdi.z.cm + (gndr.c | cntry)
## mod7_GDI: FM.z ~ gndr.c + gdi.z.cm + gndr.c:gdi.z.cm + year.c + gndr.c:year.c + (gndr.c + year.c + gndr.c:year.c | cntry)
##          npar     AIC     BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod2_GDI    8 1306540 1306627 -653262   1306524                         
## mod7_GDI   17 1301502 1301689 -650734   1301468 5055.6  9  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

### Examine trends


``` r
change_mod7_GDI<-emmeans(mod7_GDI,specs="year.c",by="gndr.c",
                              at=list(gndr.c=c(-0.5,0.5),
                                      year.c=rev(range(diff_dat$year.c)),
                                      gdi.z.cmc=0,
                                      gdi.z.cm=0),
                              disable.pbkrtest=T,
                              lmerTest.limit = 500000,infer=c(T,T),adjust="none")
change_mod7_GDI
```

```
## gndr.c = -0.5:
##  year.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21 -0.107 0.0994 30.1  -0.3104   0.0954  -1.082  0.2880
##       0 -0.283 0.1210 28.0  -0.5310  -0.0350  -2.337  0.0268
## 
## gndr.c =  0.5:
##  year.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21  0.260 0.0874 29.9   0.0812   0.4384   2.972  0.0058
##       0  0.219 0.0911 26.2   0.0322   0.4067   2.408  0.0234
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_mod7_GDI,adjust="none",infer=c(T,T))
```

```
## gndr.c = -0.5:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0   0.1755 0.170 31.4   -0.171    0.522   1.031  0.3104
## 
## gndr.c =  0.5:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0   0.0403 0.136 30.8   -0.237    0.317   0.297  0.7686
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# change in gender differences over time

change_in_diff_mod7_GDI<-emmeans(mod7_GDI,specs=c("gndr.c","year.c"),
                                      at=list(gndr.c=c(-0.5,0.5),
                                              year.c=rev(range(diff_dat$year.c)),
                                              gdi.z.cmc=0,
                                              gdi.z.cm=0),
                                      disable.pbkrtest=T,
                                      lmerTest.limit = 500000,infer=c(T,T),adjust="none")
change_in_diff_mod7_GDI
```

```
##  gndr.c year.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21 -0.107 0.0994 30.1  -0.3104   0.0954  -1.082  0.2880
##     0.5     21  0.260 0.0874 29.9   0.0812   0.4384   2.972  0.0058
##    -0.5      0 -0.283 0.1210 28.0  -0.5310  -0.0350  -2.337  0.0268
##     0.5      0  0.219 0.0911 26.2   0.0322   0.4067   2.408  0.0234
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod7_GDI,adjust="none",infer=c(T,T))
```

```
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21    -0.3673 0.0216 20.5   -0.412  -0.3224 -17.037 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)   0.1755 0.1700 31.4   -0.171   0.5225   1.031  0.3104
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     -0.3270 0.1470 31.6   -0.626  -0.0275  -2.225  0.0333
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)      0.5428 0.1600 31.4    0.216   0.8692   3.390  0.0019
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0         0.0403 0.1360 30.8   -0.237   0.3175   0.297  0.7686
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0      -0.5025 0.0389 11.8   -0.587  -0.4176 -12.924 <0.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# custom contrast to obtain the estimates for the entire span

changes_in_diff <- list(
  diff_ESS11 = c(-1,1,0,0),
  diff_ESS1 = c(0,0,-1,1))

diff_mod7_GDI<-contrast(change_in_diff_mod7_GDI,method = changes_in_diff,adjust="none",infer=c(T,T))

# gender differences at ESS1 and ESS10
diff_mod7_GDI
```

```
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11    0.367 0.0216 20.5    0.322    0.412  17.037 <0.0001
##  diff_ESS1     0.502 0.0389 11.8    0.418    0.587  12.924 <0.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# Test for whether and how much the differences have changed
pairs(diff_mod7_GDI,infer=c(T,T))
```

```
##  contrast               estimate     SE df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1   -0.135 0.0466 14   -0.235  -0.0352  -2.900  0.0117
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```


``` r
mod8_GDI<-lmer(FM.z~gndr.c+gdi.z.cm+gndr.c:gdi.z.cm+
                      year.c+
                      gndr.c:year.c+
                      gdi.z.cmc+gdi.z.cmc:gndr.c+
                      (gndr.c+year.c+gndr.c:year.c+gdi.z.cmc+gdi.z.cmc:gndr.c|cntry),
               data=diff_dat,REML=F,weights = pspwght,
                    control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

getFE(mod8_GDI)
```

```
##                   Est.   SE    df     t     p    LL   UL
## (Intercept)       0.09 0.05 26.13  1.97 0.060  0.00 0.18
## gndr.c            0.47 0.02 23.72 19.02 0.000  0.42 0.52
## gdi.z.cm         -0.04 0.05 32.87 -0.87 0.391 -0.13 0.05
## year.c            0.00 0.00 25.37 -1.10 0.282 -0.01 0.00
## gdi.z.cmc         0.09 0.06 31.14  1.55 0.132 -0.03 0.21
## gndr.c:gdi.z.cm   0.04 0.02 34.26  2.51 0.017  0.01 0.07
## gndr.c:year.c     0.00 0.00 22.06 -3.70 0.001 -0.01 0.00
## gndr.c:gdi.z.cmc -0.01 0.03 29.79 -0.36 0.724 -0.07 0.05
```

``` r
getVC(mod8_GDI)
```

```
##         grp             var1             var2 sdcor  vcov
## 1     cntry      (Intercept)             <NA>  0.26  0.07
## 2     cntry           gndr.c             <NA>  0.13  0.02
## 3     cntry           year.c             <NA>  0.02  0.00
## 4     cntry        gdi.z.cmc             <NA>  0.32  0.11
## 5     cntry    gndr.c:year.c             <NA>  0.00  0.00
## 6     cntry gndr.c:gdi.z.cmc             <NA>  0.15  0.02
## 7     cntry      (Intercept)           gndr.c -0.34 -0.01
## 8     cntry      (Intercept)           year.c -0.07  0.00
## 9     cntry      (Intercept)        gdi.z.cmc  0.22  0.02
## 10    cntry      (Intercept)    gndr.c:year.c  0.38  0.00
## 11    cntry      (Intercept) gndr.c:gdi.z.cmc -0.14 -0.01
## 12    cntry           gndr.c           year.c -0.33  0.00
## 13    cntry           gndr.c        gdi.z.cmc  0.31  0.01
## 14    cntry           gndr.c    gndr.c:year.c -0.74  0.00
## 15    cntry           gndr.c gndr.c:gdi.z.cmc -0.38 -0.01
## 16    cntry           year.c        gdi.z.cmc -0.03  0.00
## 17    cntry           year.c    gndr.c:year.c  0.20  0.00
## 18    cntry           year.c gndr.c:gdi.z.cmc -0.10  0.00
## 19    cntry        gdi.z.cmc    gndr.c:year.c  0.05  0.00
## 20    cntry        gdi.z.cmc gndr.c:gdi.z.cmc -0.77 -0.04
## 21    cntry    gndr.c:year.c gndr.c:gdi.z.cmc -0.02  0.00
## 22 Residual             <NA>             <NA>  1.00  1.00
```

``` r
anova(mod2_GDI,mod7_GDI,mod8_GDI)
```

```
## Data: diff_dat
## Models:
## mod2_GDI: FM.z ~ gndr.c + gdi.z.cm + gndr.c:gdi.z.cm + (gndr.c | cntry)
## mod7_GDI: FM.z ~ gndr.c + gdi.z.cm + gndr.c:gdi.z.cm + year.c + gndr.c:year.c + (gndr.c + year.c + gndr.c:year.c | cntry)
## mod8_GDI: FM.z ~ gndr.c + gdi.z.cm + gndr.c:gdi.z.cm + year.c + gndr.c:year.c + gdi.z.cmc + gdi.z.cmc:gndr.c + (gndr.c + year.c + gndr.c:year.c + gdi.z.cmc + gdi.z.cmc:gndr.c | cntry)
##          npar     AIC     BIC  logLik -2*log(L)   Chisq Df Pr(>Chisq)    
## mod2_GDI    8 1306540 1306627 -653262   1306524                          
## mod7_GDI   17 1301502 1301689 -650734   1301468 5055.58  9  < 2.2e-16 ***
## mod8_GDI   30 1301153 1301482 -650546   1301093  375.21 13  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

### Examine trends


``` r
# gender specific change over time

change_mod8_GDI<-emmeans(mod8_GDI,specs="year.c",by="gndr.c",
                              at=list(gndr.c=c(-0.5,0.5),
                                      year.c=rev(range(diff_dat$year.c)),
                                      gdi.z.cmc=0,
                                      gdi.z.cm=0),
                              disable.pbkrtest=T,
                              lmerTest.limit = 500000,infer=c(T,T),adjust="none")
change_mod8_GDI
```

```
## gndr.c = -0.5:
##  year.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21 -0.170 0.0745 28.7  -0.3225  -0.0178  -2.285  0.0299
##       0 -0.146 0.0509 25.3  -0.2502  -0.0408  -2.860  0.0084
## 
## gndr.c =  0.5:
##  year.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21  0.212 0.0705 28.9   0.0683   0.3566   3.015  0.0053
##       0  0.325 0.0435 26.6   0.2359   0.4144   7.483 <0.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_mod8_GDI,adjust="none",infer=c(T,T))
```

```
## gndr.c = -0.5:
##  contrast           estimate     SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0  -0.0247 0.0622 24.9   -0.153    0.103  -0.397  0.6950
## 
## gndr.c =  0.5:
##  contrast           estimate     SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0  -0.1127 0.0650 25.4   -0.246    0.021  -1.735  0.0949
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# change in gender differences over time

change_in_diff_mod8_GDI<-emmeans(mod8_GDI,specs=c("gndr.c","year.c"),
                                      at=list(gndr.c=c(-0.5,0.5),
                                              year.c=rev(range(diff_dat$year.c)),
                                              gdi.z.cmc=0,
                                              gdi.z.cm=0),
                                      disable.pbkrtest=T,
                                      lmerTest.limit = 500000,infer=c(T,T),adjust="none")
change_in_diff_mod8_GDI
```

```
##  gndr.c year.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21 -0.170 0.0745 28.7  -0.3225  -0.0178  -2.285  0.0299
##     0.5     21  0.212 0.0705 28.9   0.0683   0.3566   3.015  0.0053
##    -0.5      0 -0.146 0.0509 25.3  -0.2502  -0.0408  -2.860  0.0084
##     0.5      0  0.325 0.0435 26.6   0.2359   0.4144   7.483 <0.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod8_GDI,adjust="none",infer=c(T,T))
```

```
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21    -0.3826 0.0174 25.4   -0.418   -0.347 -21.978 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)  -0.0247 0.0622 24.9   -0.153    0.103  -0.397  0.6950
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     -0.4953 0.0700 25.2   -0.639   -0.351  -7.075 <0.0001
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)      0.3579 0.0595 24.8    0.235    0.480   6.020 <0.0001
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0        -0.1127 0.0650 25.4   -0.246    0.021  -1.735  0.0949
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0      -0.4707 0.0247 23.7   -0.522   -0.420 -19.021 <0.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# custom contrast to obtain the estimates for the entire span

changes_in_diff <- list(
  diff_ESS11 = c(-1,1,0,0),
  diff_ESS1 = c(0,0,-1,1))

diff_mod8_GDI<-contrast(change_in_diff_mod8_GDI,method = changes_in_diff,adjust="none",infer=c(T,T))

# gender differences at ESS1 and ESS10
diff_mod8_GDI
```

```
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11    0.383 0.0174 25.4    0.347    0.418  21.978 <0.0001
##  diff_ESS1     0.471 0.0247 23.7    0.420    0.522  19.021 <0.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# Test for whether and how much the differences have changed
pairs(diff_mod8_GDI,infer=c(T,T))
```

```
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1  -0.0881 0.0238 22.1   -0.137  -0.0387  -3.698  0.0013
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

## mod7_log_GDP and mod8_log_GDP: Countries' progress in log_GDP and convergence of gender gap

* First fit a model including country average of gender-equality and its interaction with gender and time and differential development for women and men (fixed and random)

* Then fit a model that adds time-specific within-country fluctuations from country's average gender-equality (fixed and random) and compare to a model without these parameters.


``` r
mod7_log_GDP<-lmer(FM.z~gndr.c+log_gdp.z.cm+gndr.c:log_gdp.z.cm+
                      year.c+
                      gndr.c:year.c+
                      (gndr.c+year.c+gndr.c:year.c|cntry),
               data=diff_dat,REML=F,weights = pspwght,
                    control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

summary(mod7_log_GDP)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: FM.z ~ gndr.c + log_gdp.z.cm + gndr.c:log_gdp.z.cm + year.c +  
##     gndr.c:year.c + (gndr.c + year.c + gndr.c:year.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##   1301506   1301693   -650736   1301472    437726 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.7008 -0.6137  0.0026  0.5976  8.3037 
## 
## Random effects:
##  Groups   Name          Variance  Std.Dev. Corr              
##  cntry    (Intercept)   0.3725430 0.61036                    
##           gndr.c        0.0431266 0.20767  -0.73             
##           year.c        0.0017776 0.04216  -0.80  0.53       
##           gndr.c:year.c 0.0001436 0.01198   0.77 -0.87 -0.73 
##  Residual               1.0039548 1.00198                    
## Number of obs: 437743, groups:  cntry, 34
## 
## Fixed effects:
##                      Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)         -0.031785   0.104913 28.204734  -0.303   0.7641    
## gndr.c               0.501122   0.036774 11.390832  13.627 2.11e-08 ***
## log_gdp.z.cm        -0.011265   0.058438 28.317880  -0.193   0.8485    
## year.c               0.005131   0.007247 32.240694   0.708   0.4840    
## gndr.c:log_gdp.z.cm -0.013045   0.018051 23.655594  -0.723   0.4770    
## gndr.c:year.c       -0.006263   0.002159 13.185755  -2.901   0.0122 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c lg_g.. year.c g.:_..
## gndr.c      -0.711                            
## lg_gdp.z.cm  0.009 -0.031                     
## year.c      -0.799  0.511  0.002              
## gndr.c:l_.. -0.004  0.007 -0.474 -0.002       
## gndr.c:yr.c  0.733 -0.865  0.030 -0.697 -0.011
```

``` r
getFE(mod7_log_GDP)
```

```
##                      Est.   SE    df     t     p    LL   UL
## (Intercept)         -0.03 0.10 28.20 -0.30 0.764 -0.25 0.18
## gndr.c               0.50 0.04 11.39 13.63 0.000  0.42 0.58
## log_gdp.z.cm        -0.01 0.06 28.32 -0.19 0.849 -0.13 0.11
## year.c               0.01 0.01 32.24  0.71 0.484 -0.01 0.02
## gndr.c:log_gdp.z.cm -0.01 0.02 23.66 -0.72 0.477 -0.05 0.02
## gndr.c:year.c       -0.01 0.00 13.19 -2.90 0.012 -0.01 0.00
```

``` r
getVC(mod7_log_GDP)
```

```
##         grp          var1          var2 sdcor  vcov
## 1     cntry   (Intercept)          <NA>  0.61  0.37
## 2     cntry        gndr.c          <NA>  0.21  0.04
## 3     cntry        year.c          <NA>  0.04  0.00
## 4     cntry gndr.c:year.c          <NA>  0.01  0.00
## 5     cntry   (Intercept)        gndr.c -0.73 -0.09
## 6     cntry   (Intercept)        year.c -0.80 -0.02
## 7     cntry   (Intercept) gndr.c:year.c  0.77  0.01
## 8     cntry        gndr.c        year.c  0.53  0.00
## 9     cntry        gndr.c gndr.c:year.c -0.87  0.00
## 10    cntry        year.c gndr.c:year.c -0.73  0.00
## 11 Residual          <NA>          <NA>  1.00  1.00
```

``` r
anova(mod2_log_GDP,mod7_log_GDP)
```

```
## Data: diff_dat
## Models:
## mod2_log_GDP: FM.z ~ gndr.c + log_gdp.z.cm + gndr.c:log_gdp.z.cm + (gndr.c | cntry)
## mod7_log_GDP: FM.z ~ gndr.c + log_gdp.z.cm + gndr.c:log_gdp.z.cm + year.c + gndr.c:year.c + (gndr.c + year.c + gndr.c:year.c | cntry)
##              npar     AIC     BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod2_log_GDP    8 1306541 1306629 -653262   1306525                         
## mod7_log_GDP   17 1301506 1301693 -650736   1301472 5052.7  9  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

### Examine trends


``` r
change_mod7_log_GDP<-emmeans(mod7_log_GDP,specs="year.c",by="gndr.c",
                              at=list(gndr.c=c(-0.5,0.5),
                                      year.c=rev(range(diff_dat$year.c)),
                                      log_gdp.z.cmc=0,
                                      log_gdp.z.cm=0),
                              disable.pbkrtest=T,
                              lmerTest.limit = 500000,infer=c(T,T),adjust="none")
change_mod7_log_GDP
```

```
## gndr.c = -0.5:
##  year.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21 -0.109 0.1000 26.0  -0.3144   0.0968  -1.088  0.2866
##       0 -0.282 0.1190 28.0  -0.5255  -0.0392  -2.379  0.0244
## 
## gndr.c =  0.5:
##  year.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21  0.261 0.0871 25.3   0.0815   0.4401   2.993  0.0061
##       0  0.219 0.0927 24.5   0.0275   0.4100   2.359  0.0266
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_mod7_log_GDP,adjust="none",infer=c(T,T))
```

```
## gndr.c = -0.5:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0    0.174 0.169 31.3   -0.171    0.518   1.028  0.3118
## 
## gndr.c =  0.5:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0    0.042 0.137 30.8   -0.238    0.322   0.306  0.7618
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# change in gender differences over time

change_in_diff_mod7_log_GDP<-emmeans(mod7_log_GDP,specs=c("gndr.c","year.c"),
                                      at=list(gndr.c=c(-0.5,0.5),
                                              year.c=rev(range(diff_dat$year.c)),
                                              log_gdp.z.cmc=0,
                                              log_gdp.z.cm=0),
                                      disable.pbkrtest=T,
                                      lmerTest.limit = 500000,infer=c(T,T),adjust="none")
change_in_diff_mod7_log_GDP
```

```
##  gndr.c year.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21 -0.109 0.1000 26.0  -0.3144   0.0968  -1.088  0.2866
##     0.5     21  0.261 0.0871 25.3   0.0815   0.4401   2.993  0.0061
##    -0.5      0 -0.282 0.1190 28.0  -0.5255  -0.0392  -2.379  0.0244
##     0.5      0  0.219 0.0927 24.5   0.0275   0.4100   2.359  0.0266
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod7_log_GDP,adjust="none",infer=c(T,T))
```

```
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21     -0.370 0.0229 19.8   -0.417  -0.3218 -16.149 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)    0.174 0.1690 31.3   -0.171   0.5176   1.028  0.3118
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0      -0.328 0.1510 31.6   -0.635  -0.0207  -2.175  0.0372
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)       0.543 0.1560 31.3    0.224   0.8622   3.470  0.0015
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0          0.042 0.1370 30.8   -0.238   0.3222   0.306  0.7618
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0       -0.501 0.0368 11.4   -0.582  -0.4205 -13.627 <0.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# custom contrast to obtain the estimates for the entire span

changes_in_diff <- list(
  diff_ESS11 = c(-1,1,0,0),
  diff_ESS1 = c(0,0,-1,1))

diff_mod7_log_GDP<-contrast(change_in_diff_mod7_log_GDP,method = changes_in_diff,adjust="none",infer=c(T,T))

# gender differences at ESS1 and ESS10
diff_mod7_log_GDP
```

```
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11    0.370 0.0229 19.8    0.322    0.417  16.149 <0.0001
##  diff_ESS1     0.501 0.0368 11.4    0.421    0.582  13.627 <0.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# Test for whether and how much the differences have changed
pairs(diff_mod7_log_GDP,infer=c(T,T))
```

```
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1   -0.132 0.0453 13.2   -0.229  -0.0337  -2.901  0.0122
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```


``` r
mod8_log_GDP<-lmer(FM.z~gndr.c+log_gdp.z.cm+gndr.c:log_gdp.z.cm+
                      year.c+
                      gndr.c:year.c+
                      log_gdp.z.cmc+log_gdp.z.cmc:gndr.c+
                      (gndr.c+year.c+gndr.c:year.c+log_gdp.z.cmc+log_gdp.z.cmc:gndr.c|cntry),
               data=diff_dat,REML=F,weights = pspwght,
                    control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

getFE(mod8_log_GDP)
```

```
##                       Est.   SE    df     t     p    LL   UL
## (Intercept)           0.04 0.12 28.35  0.36 0.723 -0.20 0.29
## gndr.c                0.51 0.04 10.25 11.29 0.000  0.41 0.61
## log_gdp.z.cm         -0.04 0.05 25.50 -0.65 0.521 -0.15 0.08
## year.c                0.00 0.01 30.00 -0.39 0.701 -0.02 0.01
## log_gdp.z.cmc         0.15 0.07 15.54  2.03 0.059 -0.01 0.30
## gndr.c:log_gdp.z.cm   0.00 0.02 23.69 -0.27 0.787 -0.04 0.03
## gndr.c:year.c        -0.01 0.00 11.83 -2.17 0.051 -0.01 0.00
## gndr.c:log_gdp.z.cmc  0.00 0.04 24.37 -0.04 0.968 -0.08 0.08
```

``` r
getVC(mod8_log_GDP)
```

```
##         grp                 var1                 var2 sdcor  vcov
## 1     cntry          (Intercept)                 <NA>  0.68  0.46
## 2     cntry               gndr.c                 <NA>  0.25  0.06
## 3     cntry               year.c                 <NA>  0.04  0.00
## 4     cntry        log_gdp.z.cmc                 <NA>  0.39  0.15
## 5     cntry        gndr.c:year.c                 <NA>  0.02  0.00
## 6     cntry gndr.c:log_gdp.z.cmc                 <NA>  0.16  0.02
## 7     cntry          (Intercept)               gndr.c -0.81 -0.14
## 8     cntry          (Intercept)               year.c -0.87 -0.03
## 9     cntry          (Intercept)        log_gdp.z.cmc  0.42  0.11
## 10    cntry          (Intercept)        gndr.c:year.c  0.86  0.01
## 11    cntry          (Intercept) gndr.c:log_gdp.z.cmc -0.56 -0.06
## 12    cntry               gndr.c               year.c  0.71  0.01
## 13    cntry               gndr.c        log_gdp.z.cmc -0.27 -0.03
## 14    cntry               gndr.c        gndr.c:year.c -0.92  0.00
## 15    cntry               gndr.c gndr.c:log_gdp.z.cmc  0.53  0.02
## 16    cntry               year.c        log_gdp.z.cmc -0.47 -0.01
## 17    cntry               year.c        gndr.c:year.c -0.86  0.00
## 18    cntry               year.c gndr.c:log_gdp.z.cmc  0.56  0.00
## 19    cntry        log_gdp.z.cmc        gndr.c:year.c  0.39  0.00
## 20    cntry        log_gdp.z.cmc gndr.c:log_gdp.z.cmc -0.47 -0.03
## 21    cntry        gndr.c:year.c gndr.c:log_gdp.z.cmc -0.66  0.00
## 22 Residual                 <NA>                 <NA>  1.00  1.00
```

``` r
anova(mod2_log_GDP,mod7_log_GDP,mod8_log_GDP)
```

```
## Data: diff_dat
## Models:
## mod2_log_GDP: FM.z ~ gndr.c + log_gdp.z.cm + gndr.c:log_gdp.z.cm + (gndr.c | cntry)
## mod7_log_GDP: FM.z ~ gndr.c + log_gdp.z.cm + gndr.c:log_gdp.z.cm + year.c + gndr.c:year.c + (gndr.c + year.c + gndr.c:year.c | cntry)
## mod8_log_GDP: FM.z ~ gndr.c + log_gdp.z.cm + gndr.c:log_gdp.z.cm + year.c + gndr.c:year.c + log_gdp.z.cmc + log_gdp.z.cmc:gndr.c + (gndr.c + year.c + gndr.c:year.c + log_gdp.z.cmc + log_gdp.z.cmc:gndr.c | cntry)
##              npar     AIC     BIC  logLik -2*log(L)   Chisq Df Pr(>Chisq)    
## mod2_log_GDP    8 1306541 1306629 -653262   1306525                          
## mod7_log_GDP   17 1301506 1301693 -650736   1301472 5052.70  9  < 2.2e-16 ***
## mod8_log_GDP   30 1301248 1301578 -650594   1301188  283.68 13  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

### Examine trends


``` r
# gender specific change over time

change_mod8_log_GDP<-emmeans(mod8_log_GDP,specs="year.c",by="gndr.c",
                              at=list(gndr.c=c(-0.5,0.5),
                                      year.c=rev(range(diff_dat$year.c)),
                                      log_gdp.z.cmc=0,
                                      log_gdp.z.cm=0),
                              disable.pbkrtest=T,
                              lmerTest.limit = 500000,infer=c(T,T),adjust="none")
change_mod8_log_GDP
```

```
## gndr.c = -0.5:
##  year.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21 -0.206 0.0930 21.2 -0.39882  -0.0124  -2.212  0.0381
##       0 -0.211 0.1370 27.0 -0.49145   0.0701  -1.540  0.1353
## 
## gndr.c =  0.5:
##  year.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21  0.163 0.0756 20.8  0.00557   0.3203   2.155  0.0431
##       0  0.296 0.1020 25.9  0.08559   0.5058   2.894  0.0076
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_mod8_log_GDP,adjust="none",infer=c(T,T))
```

```
## gndr.c = -0.5:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0  0.00507 0.191 29.0   -0.385    0.395   0.027  0.9790
## 
## gndr.c =  0.5:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0 -0.13279 0.141 27.3   -0.421    0.156  -0.944  0.3534
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# change in gender differences over time

change_in_diff_mod8_log_GDP<-emmeans(mod8_log_GDP,specs=c("gndr.c","year.c"),
                                      at=list(gndr.c=c(-0.5,0.5),
                                              year.c=rev(range(diff_dat$year.c)),
                                              log_gdp.z.cmc=0,
                                              log_gdp.z.cm=0),
                                      disable.pbkrtest=T,
                                      lmerTest.limit = 500000,infer=c(T,T),adjust="none")
change_in_diff_mod8_log_GDP
```

```
##  gndr.c year.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21 -0.206 0.0930 21.2 -0.39882  -0.0124  -2.212  0.0381
##     0.5     21  0.163 0.0756 20.8  0.00557   0.3203   2.155  0.0431
##    -0.5      0 -0.211 0.1370 27.0 -0.49145   0.0701  -1.540  0.1353
##     0.5      0  0.296 0.1020 25.9  0.08559   0.5058   2.894  0.0076
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod8_log_GDP,adjust="none",infer=c(T,T))
```

```
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21   -0.36853 0.0288 18.5  -0.4289   -0.308 -12.796 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)  0.00507 0.1910 29.0  -0.3850    0.395   0.027  0.9790
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0    -0.50131 0.1610 29.1  -0.8299   -0.173  -3.120  0.0041
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)     0.37360 0.1710 29.6   0.0247    0.722   2.188  0.0367
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0       -0.13279 0.1410 27.3  -0.4212    0.156  -0.944  0.3534
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0     -0.50639 0.0448 10.2  -0.6060   -0.407 -11.294 <0.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# custom contrast to obtain the estimates for the entire span

changes_in_diff <- list(
  diff_ESS11 = c(-1,1,0,0),
  diff_ESS1 = c(0,0,-1,1))

diff_mod8_log_GDP<-contrast(change_in_diff_mod8_log_GDP,method = changes_in_diff,adjust="none",infer=c(T,T))

# gender differences at ESS1 and ESS10
diff_mod8_log_GDP
```

```
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11    0.369 0.0288 18.5    0.308    0.429  12.796 <0.0001
##  diff_ESS1     0.506 0.0448 10.2    0.407    0.606  11.294 <0.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# Test for whether and how much the differences have changed
pairs(diff_mod8_log_GDP,infer=c(T,T))
```

```
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1   -0.138 0.0635 11.8   -0.276  0.00074  -2.171  0.0510
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
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
##  [1] misty_0.8.2           multilevel_2.8        MASS_7.3-65           tibble_3.3.1         
##  [5] apaTables_2.0.8       stringr_1.6.0         tidyr_1.3.2           r2mlm_0.3.8          
##  [9] nlme_3.1-169          Hmisc_5.2-5           ggpubr_0.6.3          metafor_5.0-1        
## [13] numDeriv_2016.8-1.1   metadat_1.6-0         lmerTest_3.2-1        ggflags_0.0.4        
## [17] finalfit_1.1.0        ggplot2_4.0.3         MetBrewer_0.2.0       vjihelpers_0.0.0.9000
## [21] emmeans_2.0.3         lme4_2.0-1            Matrix_1.7-5          multid_1.0.2         
## [25] dplyr_1.2.1           rio_1.3.0             knitr_1.51            rmarkdown_2.31       
## 
## loaded via a namespace (and not attached):
##   [1] mathjaxr_2.0-0     RColorBrewer_1.1-3 rstudioapi_0.18.0  jsonlite_2.0.0     shape_1.4.6.1     
##   [6] magrittr_2.0.5     estimability_1.5.1 jomo_2.7-6         farver_2.1.2       nloptr_2.2.1      
##  [11] fs_2.1.0           vctrs_0.7.3        memoise_2.0.1      minqa_1.2.8        base64enc_0.1-6   
##  [16] rstatix_0.7.3      htmltools_0.5.9    forcats_1.0.1      usethis_3.2.1      haven_2.5.5       
##  [21] broom_1.0.13       cellranger_1.1.0   Formula_1.2-5      pROC_1.19.0.1      mitml_0.4-5       
##  [26] sass_0.4.10        bslib_0.10.0       htmlwidgets_1.6.4  plyr_1.8.9         cachem_1.1.0      
##  [31] lifecycle_1.0.5    iterators_1.0.14   pkgconfig_2.0.3    R6_2.6.1           fastmap_1.2.0     
##  [36] rbibutils_2.4.1    digest_0.6.39      colorspace_2.1-2   pkgload_1.5.2      labeling_0.4.3    
##  [41] mgcv_1.9-4         abind_1.4-8        compiler_4.6.0     withr_3.0.2        htmlTable_2.5.0   
##  [46] S7_0.2.2           backports_1.5.1    carData_3.0-6      psych_2.6.3        pkgbuild_1.4.8    
##  [51] R.utils_2.13.0     ggsignif_0.6.4     pan_1.9            sessioninfo_1.2.3  pbivnorm_0.6.0    
##  [56] tools_4.6.0        foreign_0.8-91     otel_0.2.0         zip_2.3.3          nnet_7.3-20       
##  [61] quadprog_1.5-8     R.oo_1.27.1        glue_1.8.1         grid_4.6.0         checkmate_2.3.4   
##  [66] cluster_2.1.8.2    generics_0.1.4     gtable_0.3.6       tzdb_0.5.0         R.methodsS3_1.8.2 
##  [71] data.table_1.18.4  hms_1.1.4          car_3.1-5          utf8_1.2.6         foreach_1.5.2     
##  [76] pillar_1.11.1      rockchalk_1.8.164  splines_4.6.0      lattice_0.22-9     survival_3.8-6    
##  [81] kutils_1.73        tidyselect_1.2.1   reformulas_0.4.4   gridExtra_2.3      grImport2_0.3-3   
##  [86] stats4_4.6.0       xfun_0.57          devtools_2.5.2     stringi_1.8.7      yaml_2.3.12       
##  [91] boot_1.3-32        evaluate_1.0.5     codetools_0.2-20   cli_3.6.6          rpart_4.1.27      
##  [96] xtable_1.8-8       Rdpack_2.6.6       lavaan_0.6-21      jquerylib_0.1.4    readxl_1.4.5      
## [101] Rcpp_1.1.1-1.1     png_0.1-9          coda_0.19-4.1      XML_3.99-0.23      parallel_4.6.0    
## [106] ellipsis_0.3.3     readr_2.2.0        jpeg_0.1-11        glmnet_5.0         mvtnorm_1.3-7     
## [111] scales_1.4.0       writexl_1.5.4      openxlsx_4.2.8.1   purrr_1.2.2        rlang_1.2.0       
## [116] cowplot_1.2.0      mnormt_2.1.2       mice_3.19.0
```

