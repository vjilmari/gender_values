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
```

# Data

* Loads the preprocessed European Social Survey data file made within "Preparations_ESS" code


``` r
load("../data/fdat.rdata")
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
## 450838   2586    774    401    248    185    129    135    142    165  34952
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
##  AL_6  AT_1  AT_2  AT_3  AT_7  AT_8  AT_9  BE_1 BE_10  BE_2  BE_3  BE_4  BE_5  BE_6  BE_7  BE_8 
##  1117  2254  2198  2348  1795  1993  2477  1830  1334  1771  1796  1754  1699  1862  1767  1759 
##  BE_9 BG_10  BG_3  BG_4  BG_5  BG_6  BG_9  CH_1 CH_10  CH_2  CH_3  CH_4  CH_5  CH_6  CH_7  CH_8 
##  1756  2697  1295  2144  2371  2179  1926  2024  1505  2110  1780  1753  1491  1483  1521  1504 
##  CH_9  CY_3  CY_4  CY_5  CY_6  CY_9  CZ_1 CZ_10  CZ_2  CZ_4  CZ_5  CZ_6  CZ_7  CZ_8  CZ_9  DE_1 
##  1517   978  1210  1053  1110   773  1208  2369  2557  1986  2335  1973  1862  2252  2343  2819 
##  DE_2  DE_3  DE_4  DE_5  DE_6  DE_7  DE_8  DE_9  DK_1  DK_2  DK_3  DK_4  DK_5  DK_6  DK_7  DK_9 
##  2840  2884  2732  3007  2935  3006  2821  2328  1470  1458  1461  1581  1564  1621  1483  1554 
## EE_10  EE_2  EE_3  EE_4  EE_5  EE_6  EE_7  EE_8  EE_9  ES_1  ES_2  ES_3  ES_4  ES_5  ES_6  ES_7 
##  1538  1948  1466  1646  1793  2345  2036  2007  1899  1712  1623  1847  2562  1881  1871  1907 
##  ES_8  ES_9  FI_1 FI_10  FI_2  FI_3  FI_4  FI_5  FI_6  FI_7  FI_8  FI_9  FR_1 FR_10  FR_2  FR_3 
##  1929  1619  1763  1561  1701  1649  1901  1649  2158  2050  1903  1735  1355  1951  1699  1983 
##  FR_4  FR_5  FR_6  FR_7  FR_8  FR_9  GB_1 GB_10  GB_2  GB_3  GB_4  GB_5  GB_6  GB_7  GB_8  GB_9 
##  2067  1723  1960  1902  2057  1982  1798  1131  1864  2353  2311  2374  2261  2231  1942  2183 
##  GR_1 GR_10  GR_2  GR_4  GR_5 HR_10  HR_4  HR_5  HR_9  HU_1 HU_10  HU_2  HU_3  HU_4  HU_5  HU_6 
##  2551  2768  2399  2063  2669  1564  1430  1601  1781  1634  1816  1460  1462  1430  1473  1968 
##  HU_7  HU_8  HU_9  IE_1 IE_10  IE_2  IE_3  IE_4  IE_5  IE_6  IE_7  IE_8  IE_9  IL_1  IL_4  IL_5 
##  1520  1458  1643  1916  1751  1187  1589  1757  2400  2616  2380  2746  2189  2279  2382  2212 
##  IL_6  IL_7  IL_8 IS_10  IS_2  IS_6  IS_8  IS_9 IT_10  IT_6  IT_8  IT_9 LT_10  LT_5  LT_6  LT_7 
##  2378  2351  2366   886   524   739   841   844  2573   909  2531  2660  1606  1632  2108  2241 
##  LT_8  LT_9  LU_2  LV_4  LV_9 ME_10  ME_9 MK_10  NL_1 NL_10  NL_2  NL_3  NL_4  NL_5  NL_6  NL_7 
##  2079  1677  1614  1970   891  1248  1188  1400  2337  1466  1858  1860  1724  1801  1828  1823 
##  NL_8  NL_9  NO_1 NO_10  NO_2  NO_3  NO_4  NO_5  NO_6  NO_7  NO_8  NO_9  PL_1  PL_2  PL_3  PL_4 
##  1669  1657  1819  1408  1575  1550  1391  1530  1610  1423  1530  1396  2065  1683  1685  1596 
##  PL_5  PL_6  PL_7  PL_8  PL_9  PT_1 PT_10  PT_2  PT_3  PT_4  PT_5  PT_6  PT_7  PT_8  PT_9  RO_4 
##  1719  1866  1594  1675  1443  1482  1827  2024  2182  2337  2139  2138  1242  1254  1045  2104 
##  RS_9  RU_3  RU_4  RU_5  RU_6  RU_8  SE_1  SE_2  SE_3  SE_4  SE_5  SE_6  SE_7  SE_8  SE_9  SI_1 
##  1969  2339  2446  2557  2429  2374  1682  1678  1604  1556  1463  1838  1761  1526  1510  1488 
## SI_10  SI_2  SI_3  SI_4  SI_5  SI_6  SI_7  SI_8  SI_9 SK_10  SK_2  SK_3  SK_4  SK_5  SK_6  SK_9 
##  1232  1384  1465  1257  1369  1244  1189  1295  1307  1395  1425  1711  1789  1803  1827  1061 
##  TR_2  TR_4  UA_2  UA_3  UA_4  UA_5  UA_6  XK_6 
##  1790  2305  1896  1885  1766  1779  2064  1244
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
## [1] 248
```

``` r
length(unique(fdat$cntry_time))*200
```

```
## [1] 49600
```

### Coefficients for value variables in training set


``` r
round(coefficients(value_typ$cv.mod,s = "lambda.min"),2)
```

```
## 11 x 1 sparse Matrix of class "dgCMatrix"
##             lambda.min
## (Intercept)       0.64
## con               0.09
## tra              -0.03
## ben              -0.30
## uni              -0.06
## sdi               0.10
## sti               0.05
## hed               0.06
## ach               0.07
## pow               0.14
## sec              -0.16
```

``` r
plot(value_typ$cv.mod)
```

![](Analysis_youth_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

### Description of gender differences/prediction in testing dataset


``` r
# save to summary tab
sum_tab<-value_typ$D
# print the country X time -fold results
round(sum_tab,2)
```

```
##        n.1  n.0   m.1   m.0 sd.1 sd.0 pooled.sd diff    D pcc.1 pcc.0 pcc.total  auc pooled.sd.1
## AL_6   419  498  0.06 -0.01 0.34 0.36      0.35 0.07 0.20  0.56  0.50      0.53 0.55        0.37
## AT_1   941 1113  0.16 -0.04 0.42 0.44      0.43 0.19 0.45  0.66  0.53      0.59 0.63        0.37
## AT_2   910 1088  0.20 -0.02 0.41 0.43      0.42 0.21 0.50  0.69  0.52      0.60 0.64        0.37
## AT_3   994 1154  0.19 -0.02 0.43 0.40      0.41 0.21 0.51  0.69  0.52      0.60 0.65        0.37
## AT_7   753  842  0.08 -0.10 0.35 0.37      0.36 0.18 0.50  0.59  0.61      0.60 0.63        0.37
## AT_8   792 1001  0.10 -0.04 0.40 0.42      0.41 0.15 0.36  0.61  0.55      0.57 0.60        0.37
## AT_9  1043 1234  0.07 -0.12 0.36 0.39      0.38 0.18 0.49  0.58  0.62      0.61 0.64        0.37
## BE_1   845  785  0.07 -0.07 0.37 0.38      0.37 0.13 0.35  0.61  0.58      0.60 0.61        0.37
## BE_10  570  564  0.04 -0.12 0.33 0.34      0.33 0.16 0.47  0.56  0.63      0.60 0.63        0.37
## BE_2   770  801  0.09 -0.09 0.34 0.37      0.36 0.18 0.50  0.60  0.58      0.59 0.64        0.37
## BE_3   740  856  0.06 -0.08 0.35 0.35      0.35 0.14 0.40  0.61  0.56      0.58 0.61        0.37
## BE_4   760  794  0.09 -0.07 0.34 0.34      0.34 0.16 0.47  0.63  0.56      0.59 0.63        0.37
## BE_5   719  780  0.06 -0.08 0.34 0.36      0.35 0.14 0.42  0.59  0.59      0.59 0.62        0.37
## BE_6   811  851  0.07 -0.09 0.32 0.35      0.33 0.15 0.45  0.60  0.57      0.58 0.62        0.37
## BE_7   795  772  0.05 -0.08 0.33 0.33      0.33 0.13 0.40  0.56  0.56      0.56 0.61        0.37
## BE_8   784  775  0.07 -0.09 0.33 0.34      0.33 0.16 0.47  0.59  0.57      0.58 0.62        0.37
## BE_9   765  791  0.08 -0.07 0.32 0.36      0.34 0.14 0.42  0.60  0.56      0.58 0.61        0.37
## BG_10 1171 1326  0.18  0.02 0.38 0.42      0.40 0.15 0.38  0.70  0.46      0.57 0.60        0.37
## BG_3   409  686 -0.02 -0.19 0.42 0.45      0.44 0.17 0.40  0.50  0.68      0.62 0.62        0.37
## BG_4   848 1096 -0.07 -0.26 0.43 0.46      0.44 0.18 0.41  0.45  0.71      0.60 0.62        0.37
## BG_5   942 1229  0.05 -0.12 0.37 0.42      0.40 0.17 0.42  0.57  0.59      0.58 0.62        0.37
## BG_6   830 1149  0.06 -0.13 0.37 0.41      0.40 0.19 0.47  0.58  0.61      0.60 0.63        0.37
## BG_9   763  963  0.08 -0.08 0.38 0.44      0.41 0.16 0.38  0.58  0.58      0.58 0.61        0.37
## CH_1   874  950  0.06 -0.09 0.33 0.34      0.33 0.15 0.45  0.58  0.59      0.59 0.63        0.37
## CH_10  672  633  0.06 -0.11 0.32 0.33      0.32 0.17 0.53  0.60  0.61      0.60 0.65        0.37
## CH_2   834 1076  0.06 -0.10 0.34 0.35      0.35 0.16 0.46  0.58  0.61      0.60 0.63        0.37
## CH_3   702  878  0.03 -0.12 0.34 0.32      0.33 0.15 0.45  0.54  0.66      0.61 0.62        0.37
## CH_4   697  856  0.09 -0.14 0.33 0.34      0.34 0.24 0.70  0.61  0.67      0.64 0.70        0.37
## CH_5   662  629  0.04 -0.11 0.32 0.32      0.32 0.15 0.47  0.56  0.64      0.60 0.64        0.37
## CH_6   643  640  0.06 -0.10 0.30 0.33      0.32 0.15 0.49  0.59  0.60      0.59 0.63        0.37
## CH_7   659  662  0.08 -0.11 0.31 0.33      0.32 0.19 0.60  0.61  0.64      0.63 0.67        0.37
## CH_8   681  623  0.06 -0.10 0.31 0.32      0.31 0.17 0.53  0.58  0.64      0.61 0.65        0.37
## CH_9   667  650  0.05 -0.11 0.31 0.33      0.32 0.17 0.52  0.57  0.64      0.60 0.65        0.37
## CY_3   368  410 -0.01 -0.12 0.33 0.41      0.37 0.10 0.28  0.52  0.59      0.56 0.57        0.37
## CY_4   512  498 -0.01 -0.12 0.39 0.39      0.39 0.11 0.28  0.49  0.59      0.54 0.57        0.37
## CY_5   375  478 -0.07 -0.20 0.35 0.45      0.41 0.12 0.30  0.45  0.65      0.56 0.58        0.37
## CY_6   391  519 -0.12 -0.26 0.38 0.42      0.41 0.14 0.34  0.38  0.72      0.57 0.59        0.37
## CY_9   263  310 -0.18 -0.26 0.40 0.40      0.40 0.08 0.19  0.32  0.74      0.54 0.55        0.37
## CZ_1   472  536  0.11 -0.17 0.46 0.44      0.45 0.28 0.63  0.60  0.67      0.64 0.68        0.37
## CZ_10  928 1241  0.32  0.15 0.39 0.45      0.42 0.17 0.40  0.80  0.34      0.54 0.61        0.37
## CZ_2  1089 1268  0.25 -0.04 0.41 0.47      0.45 0.29 0.64  0.74  0.52      0.62 0.68        0.37
## CZ_4   863  923  0.26  0.03 0.41 0.47      0.44 0.23 0.52  0.76  0.47      0.61 0.64        0.37
## CZ_5  1069 1066  0.27  0.04 0.41 0.45      0.43 0.23 0.54  0.76  0.43      0.60 0.64        0.37
## CZ_6   875  898  0.28  0.08 0.41 0.47      0.44 0.21 0.47  0.78  0.42      0.60 0.62        0.37
## CZ_7   704  958  0.31  0.09 0.40 0.43      0.41 0.22 0.54  0.80  0.39      0.56 0.65        0.37
## CZ_8   989 1063  0.30  0.16 0.37 0.41      0.39 0.13 0.34  0.79  0.32      0.55 0.59        0.37
## CZ_9   921 1222  0.25  0.08 0.40 0.43      0.42 0.18 0.42  0.75  0.42      0.56 0.62        0.37
## DE_1  1249 1370  0.05 -0.13 0.41 0.42      0.41 0.18 0.43  0.56  0.64      0.60 0.63        0.37
## DE_2  1267 1373  0.09 -0.11 0.39 0.39      0.39 0.20 0.51  0.58  0.62      0.60 0.64        0.37
## DE_3  1319 1365  0.09 -0.10 0.38 0.39      0.39 0.19 0.49  0.60  0.59      0.60 0.64        0.37
## DE_4  1342 1190  0.04 -0.15 0.38 0.38      0.38 0.20 0.52  0.55  0.68      0.61 0.65        0.37
## DE_5  1444 1363 -0.01 -0.17 0.35 0.36      0.36 0.16 0.45  0.49  0.69      0.59 0.63        0.37
## DE_6  1385 1350 -0.02 -0.20 0.33 0.37      0.35 0.18 0.52  0.50  0.70      0.60 0.65        0.37
## DE_7  1423 1383 -0.07 -0.25 0.34 0.36      0.35 0.18 0.52  0.43  0.77      0.60 0.65        0.37
## DE_8  1392 1229 -0.08 -0.22 0.35 0.36      0.35 0.15 0.42  0.40  0.74      0.56 0.62        0.37
## DE_9  1094 1034 -0.08 -0.25 0.35 0.35      0.35 0.17 0.49  0.41  0.75      0.58 0.63        0.37
## DK_1   651  619  0.12 -0.03 0.38 0.35      0.36 0.15 0.42  0.65  0.51      0.58 0.62        0.37
## DK_2   608  650  0.15 -0.05 0.37 0.39      0.38 0.20 0.53  0.69  0.54      0.61 0.65        0.37
##       pooled.sd.0 pooled.sd.total d.sd.total
## AL_6          0.4            0.39       0.18
## AT_1          0.4            0.39       0.50
## AT_2          0.4            0.39       0.55
## AT_3          0.4            0.39       0.54
## AT_7          0.4            0.39       0.46
## AT_8          0.4            0.39       0.38
## AT_9          0.4            0.39       0.47
## BE_1          0.4            0.39       0.34
## BE_10         0.4            0.39       0.41
## BE_2          0.4            0.39       0.46
## BE_3          0.4            0.39       0.36
## BE_4          0.4            0.39       0.41
## BE_5          0.4            0.39       0.37
## BE_6          0.4            0.39       0.39
## BE_7          0.4            0.39       0.34
## BE_8          0.4            0.39       0.41
## BE_9          0.4            0.39       0.37
## BG_10         0.4            0.39       0.39
## BG_3          0.4            0.39       0.45
## BG_4          0.4            0.39       0.48
## BG_5          0.4            0.39       0.43
## BG_6          0.4            0.39       0.49
## BG_9          0.4            0.39       0.41
## CH_1          0.4            0.39       0.39
## CH_10         0.4            0.39       0.44
## CH_2          0.4            0.39       0.42
## CH_3          0.4            0.39       0.38
## CH_4          0.4            0.39       0.61
## CH_5          0.4            0.39       0.39
## CH_6          0.4            0.39       0.40
## CH_7          0.4            0.39       0.49
## CH_8          0.4            0.39       0.43
## CH_9          0.4            0.39       0.43
## CY_3          0.4            0.39       0.27
## CY_4          0.4            0.39       0.29
## CY_5          0.4            0.39       0.32
## CY_6          0.4            0.39       0.36
## CY_9          0.4            0.39       0.20
## CZ_1          0.4            0.39       0.73
## CZ_10         0.4            0.39       0.44
## CZ_2          0.4            0.39       0.74
## CZ_4          0.4            0.39       0.59
## CZ_5          0.4            0.39       0.60
## CZ_6          0.4            0.39       0.53
## CZ_7          0.4            0.39       0.58
## CZ_8          0.4            0.39       0.35
## CZ_9          0.4            0.39       0.46
## DE_1          0.4            0.39       0.46
## DE_2          0.4            0.39       0.52
## DE_3          0.4            0.39       0.48
## DE_4          0.4            0.39       0.51
## DE_5          0.4            0.39       0.41
## DE_6          0.4            0.39       0.47
## DE_7          0.4            0.39       0.47
## DE_8          0.4            0.39       0.38
## DE_9          0.4            0.39       0.44
## DK_1          0.4            0.39       0.40
## DK_2          0.4            0.39       0.52
##  [ reached 'max' / getOption("max.print") -- omitted 190 rows ]
```

``` r
# range in gender differences across folds with fold-specific SDs
range(sum_tab$D)
```

```
## [1] 0.06219156 0.70223254
```

``` r
# range in gender differences across folds with SD pooled across all folds
range(sum_tab$d.sd.total)
```

```
## [1] 0.06708605 0.73749832
```

``` r
# mean gender difference across folds with SD pooled across all folds
mean(sum_tab$d.sd.total)
```

```
## [1] 0.4158013
```

``` r
# average probability of correct classification (pcc)
mean(sum_tab$pcc.total)
```

```
## [1] 0.5805718
```

``` r
# pcc for men
mean(sum_tab$pcc.1)
```

```
## [1] 0.5964323
```

``` r
# pcc for women
mean(sum_tab$pcc.0)
```

```
## [1] 0.5720302
```

``` r
# average area under the curve
mean(sum_tab$auc)
```

```
## [1] 0.618998
```

``` r
# print smallest and largest gender differences
sum_tab[sum_tab$d.sd.total==min(sum_tab$d.sd.total),]
```

```
##       n.1 n.0        m.1         m.0      sd.1      sd.0 pooled.sd       diff          D     pcc.1
## MK_10 538 662 0.01341395 -0.01254427 0.4035074 0.4283396 0.4173914 0.02595822 0.06219156 0.5037175
##          pcc.0 pcc.total       auc pooled.sd.1 pooled.sd.0 pooled.sd.total d.sd.total
## MK_10 0.521148 0.5133333 0.5189917   0.3738373   0.3976181       0.3869392 0.06708605
```

``` r
sum_tab[sum_tab$d.sd.total==max(sum_tab$d.sd.total),]
```

```
##       n.1  n.0       m.1         m.0      sd.1      sd.0 pooled.sd     diff         D     pcc.1
## CZ_2 1089 1268 0.2469659 -0.03840109 0.4112922 0.4745104 0.4464178 0.285367 0.6392375 0.7410468
##          pcc.0 pcc.total       auc pooled.sd.1 pooled.sd.0 pooled.sd.total d.sd.total
## CZ_2 0.5157729 0.6198557 0.6766359   0.3738373   0.3976181       0.3869392  0.7374983
```

``` r
# correlation between men and women in male-typicality across the folds
cor(sum_tab$m.1,sum_tab$m.0)
```

```
## [1] 0.9199494
```

``` r
# correlation between men and women in male-typicality deviations across the folds
cor(sum_tab$sd.1,sum_tab$sd.0)
```

```
## [1] 0.8655195
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
##        AT   BE   BG   CH   CY   CZ   DE   DK   EE   ES   FI   FR   GB   GR   HR   HU   IE   IL   IS
##   1  2054 1630    0 1824    0 1008 2619 1270    0 1512 1563 1155 1598 2351    0 1434 1716 2079    0
##   2  1998 1571    0 1910    0 2357 2640 1258 1748 1423 1501 1499 1664 2199    0 1260  987    0  324
##   3  2148 1596 1095 1580  778    0 2684 1261 1266 1647 1449 1783 2153    0    0 1262 1389    0    0
##   4     0 1554 1944 1553 1010 1786 2532 1381 1446 2362 1701 1867 2111 1863 1230 1230 1557 2182    0
##   5     0 1499 2171 1291  853 2135 2807 1364 1593 1681 1449 1523 2174 2469 1401 1273 2200 2012    0
##   6     0 1662 1979 1283  910 1773 2735 1421 2145 1671 1958 1760 2061    0    0 1768 2416 2178  539
##   7  1595 1567    0 1321    0 1662 2806 1283 1836 1707 1850 1702 2031    0    0 1320 2180 2151    0
##   8  1793 1559    0 1304    0 2052 2621    0 1807 1729 1703 1857 1742    0    0 1258 2546 2166  641
##   9  2277 1556 1726 1317  573 2143 2128 1354 1699 1419 1535 1782 1983    0 1581 1443 1989    0  644
##   10    0 1134 2497 1305    0 2169    0    0 1338    0 1361 1751  931 2568 1364 1616 1551    0  686
##     
##        IT   LT   LV   ME   NL   NO   PL   PT   RU   SE   SI   SK   TR   UA
##   1     0    0    0    0 2137 1619 1865 1282    0 1482 1288    0    0    0
##   2     0    0    0    0 1658 1375 1483 1824    0 1478 1184 1225 1590 1696
##   3     0    0    0    0 1660 1350 1485 1982 2139 1404 1265 1511    0 1685
##   4     0    0 1770    0 1524 1191 1396 2137 2246 1356 1057 1589 2105 1566
##   5     0 1432    0    0 1601 1330 1519 1939 2357 1263 1169 1603    0 1579
##   6   709 1908    0    0 1628 1410 1666 1938 2229 1638 1044 1627    0 1864
##   7     0 2041    0    0 1623 1223 1394 1042    0 1561  989    0    0    0
##   8  2331 1879    0    0 1469 1330 1475 1054 2174 1326 1095    0    0    0
##   9  2460 1477  691  988 1457 1196 1243  845    0 1310 1107  861    0    0
##   10 2373 1406    0 1048 1266 1208    0 1627    0    0 1032 1195    0    0
```

``` r
# range of sample sizes
range(table(diff_dat$cntry))
```

```
## [1]  2036 23572
```

``` r
# value-based gender-typicality histogram
hist(diff_dat$pred)
```

![](Analysis_youth_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

``` r
# scale/standardize with SD pooled across all country-time-gender folds
FM_pooled_sd<-value_typ$D[1,"pooled.sd.total"]
FM_pooled_sd
```

```
## [1] 0.3869392
```

``` r
# standardized
diff_dat$FM.z<-diff_dat$pred/FM_pooled_sd
hist(diff_dat$FM.z)
```

![](Analysis_youth_files/figure-html/unnamed-chunk-8-2.png)<!-- -->

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
  pivot_longer(
    cols = starts_with("gii_") & !matches("avg"),  # <-- exclude the "_avg" columns
    names_to = "year",
    values_to = "gii",
    names_prefix = "gii_"
  ) %>%
  mutate(year = as.integer(year)) %>%
  filter(year >= 2002 & year <= 2022) %>%
  mutate(gei = 1 - gii) %>%
  select(ISO2, iso3, country, year, gii, gei, gii_2002_2022_avg, gei_2002_2022_avg)


# describe gei average
psych::describe(GII_in_ESS_d$gei_2002_2022_avg)
```

```
##    vars  n mean   sd median trimmed  mad  min  max range  skew kurtosis   se
## X1    1 32 0.87 0.07   0.87    0.87 0.06 0.62 0.96  0.34 -1.11     1.48 0.01
```

``` r
# one is missing, see which one
GII_in_ESS_d[is.na(GII_in_ESS_d$gei_2002_2022_avg),"country"]
```

```
## [1] "Ukraine"
```

``` r
# get means and SDs for standardizing
gei_mean<-mean(GII_in_ESS_d$gei_2002_2022_avg,na.rm=T)
gei_sd<-sd(GII_in_ESS_d$gei_2002_2022_avg,na.rm=T)

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

# averages across 2002-2022
#long_GII_in_ESS_d$gei_2002_2022_avg.z<-(long_GII_in_ESS_d$gei_2002_2022_avg-gei_mean)/gei_sd

#long_GII_in_ESS_d$gei.cmc<-long_GII_in_ESS_d$gei-long_GII_in_ESS_d$gei_2002_2022_avg

# add year to ESS data-frame

diff_dat$year<-
  case_when(
    diff_dat$essround.c==-4.5~2002,  
    diff_dat$essround.c==-3.5~2004,
    diff_dat$essround.c==-2.5~2006,
    diff_dat$essround.c==-1.5~2008,
    diff_dat$essround.c==-0.5~2010,
    diff_dat$essround.c==0.5~2012,
    diff_dat$essround.c==1.5~2014,
    diff_dat$essround.c==2.5~2016,
    diff_dat$essround.c==3.5~2018,
    diff_dat$essround.c==4.5~2020
  )

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
  filter(year >= 2002 & year <= 2022) %>%
  select(ISO2, iso3, country, year, gdi, gdi_2002_2022_avg)

# describe gdi average
psych::describe(GDI_in_ESS_d$gdi_2002_2022_avg)
```

```
##    vars  n mean   sd median trimmed  mad min  max range  skew kurtosis se
## X1    1 33 0.99 0.03   0.99    0.98 0.02 0.9 1.04  0.13 -0.35     1.27  0
```

``` r
# get means and SDs for standardizing
gdi_mean<-mean(GDI_in_ESS_d$gdi_2002_2022_avg,na.rm=T)
gdi_sd<-sd(GDI_in_ESS_d$gdi_2002_2022_avg,na.rm=T)

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

#GDP_in_ESS_d
# Make long format data of GDP, so that country x year has own rows

long_GDP_in_ESS_d <- GDP_in_ESS_d %>%
  pivot_longer(
    cols = starts_with("gdp_") & !matches("avg"),  # <-- exclude the "_avg" columns
    names_to = "year",
    values_to = "gdp",
    names_prefix = "gdp_"
  ) %>%
  mutate(year = as.integer(year)) %>%
  filter(year >= 2002 & year <= 2022) %>%
  select(ISO2, Country.Name, year, gdp, gdp_2002_2022_avg,log_gdp_2002_2022_avg) %>%
  mutate(log_gdp=log(gdp))

#View(long_GDP_in_ESS_d)

# describe gdp average
psych::describe(GDP_in_ESS_d$gdp_2002_2022_avg)
```

```
##    vars  n     mean      sd   median  trimmed      mad      min      max    range skew kurtosis
## X1    1 33 44309.55 16988.2 40528.07 43238.06 16282.88 16406.88 85115.96 68709.08 0.48    -0.58
##         se
## X1 2957.27
```

``` r
psych::describe(GDP_in_ESS_d$log_gdp_2002_2022_avg)
```

```
##    vars  n  mean  sd median trimmed  mad  min   max range  skew kurtosis   se
## X1    1 33 10.62 0.4  10.61   10.64 0.43 9.71 11.35  1.65 -0.25    -0.67 0.07
```

``` r
# get means and SDs for standardizing
log_gdp_mean<-mean(GDP_in_ESS_d$log_gdp_2002_2022_avg,na.rm=T)
log_gdp_sd<-sd(GDP_in_ESS_d$log_gdp_2002_2022_avg,na.rm=T)

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
##   33
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
  filter(year >= 2002 & year <= 2022) %>%
  select(ISO2, cname, year, gggi, GGGI_2002_2022_avg)

# describe gggi average
psych::describe(GGGI_in_ESS_d$GGGI_2002_2022_avg)
```

```
##    vars  n mean   sd median trimmed  mad  min  max range skew kurtosis   se
## X1    1 33 0.73 0.05   0.73    0.73 0.04 0.61 0.86  0.25 0.38     0.19 0.01
```

``` r
# get means and SDs for standardizing
gggi_mean<-mean(GGGI_in_ESS_d$GGGI_2002_2022_avg,na.rm=T)
gggi_sd<-sd(GGGI_in_ESS_d$GGGI_2002_2022_avg,na.rm=T)

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



# Analysis

## Data subset


``` r
diff_dat<-diff_dat %>%
  filter(agea>17 & agea <30)
```

## mod0: Random intercept model


``` r
mod0<-lmer(FM.z~(1|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod0)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method [lmerModLmerTest]
## Formula: FM.z ~ (1 | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##  179022.6  179049.8  -89508.3  179016.6     63044 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.1575 -0.5943  0.0206  0.5918  9.2626 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  cntry    (Intercept) 0.07234  0.269   
##  Residual             1.03704  1.018   
## Number of obs: 63047, groups:  cntry, 33
## 
## Fixed effects:
##             Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)  0.43295    0.04703 32.88379   9.206 1.27e-10 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

``` r
getVC(mod0)
```

```
##        grp        var1 var2 sdcor vcov
## 1    cntry (Intercept) <NA>  0.27 0.07
## 2 Residual        <NA> <NA>  1.02 1.04
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
## mean variation  0.06520463     NA       1
## sigma2          0.93479537      1      NA
## 
## $R2s
##          total within between
## f1  0.00000000      0      NA
## f2  0.00000000     NA       0
## v   0.00000000      0      NA
## m   0.06520463     NA       1
## f   0.00000000     NA      NA
## fv  0.00000000      0      NA
## fvm 0.06520463     NA      NA
```

## mod1: Gender fixed effect


``` r
mod1<-lmer(FM.z~gndr.c+(1|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod1)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method [lmerModLmerTest]
## Formula: FM.z ~ gndr.c + (1 | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##  176345.0  176381.2  -88168.5  176337.0     63043 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.6291 -0.5883  0.0252  0.5952  9.0799 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  cntry    (Intercept) 0.07268  0.2696  
##  Residual             0.99387  0.9969  
## Number of obs: 63047, groups:  cntry, 33
## 
## Fixed effects:
##              Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept) 4.322e-01  4.713e-02 3.289e+01   9.171  1.4e-10 ***
## gndr.c      3.828e-01  7.316e-03 6.302e+04  52.321  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##        (Intr)
## gndr.c 0.000
```

``` r
getFE(mod1,round=3)
```

```
##              Est.    SE        df      t p    LL    UL
## (Intercept) 0.432 0.047    32.886  9.171 0 0.336 0.528
## gndr.c      0.383 0.007 63015.035 52.321 0 0.368 0.397
```

``` r
getVC(mod1)
```

```
##        grp        var1 var2 sdcor vcov
## 1    cntry (Intercept) <NA>  0.27 0.07
## 2 Residual        <NA> <NA>  1.00 0.99
```

``` r
r2mlm(mod1,bargraph = F)
```

```
## $Decompositions
##                      total
## fixed           0.03316947
## slope variation 0.00000000
## mean variation  0.06588070
## sigma2          0.90094982
## 
## $R2s
##          total
## f   0.03316947
## v   0.00000000
## m   0.06588070
## fv  0.03316947
## fvm 0.09905018
```

* Alternative model with country-mean centered gender


``` r
diff_dat_gmc<-group_mean_center(data=diff_dat,
                                group.var = "cntry",vars = c("gndr.c","essround.c"))

mod1_alt<-lmer(FM.z~gndr.c.gmc+(1|cntry),data=diff_dat_gmc,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod1_alt)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method [lmerModLmerTest]
## Formula: FM.z ~ gndr.c.gmc + (1 | cntry)
##    Data: diff_dat_gmc
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##  176343.7  176379.9  -88167.9  176335.7     63043 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.6293 -0.5883  0.0252  0.5952  9.0797 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  cntry    (Intercept) 0.06995  0.2645  
##  Residual             0.99387  0.9969  
## Number of obs: 63047, groups:  cntry, 33
## 
## Fixed effects:
##              Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept) 4.237e-01  4.624e-02 3.288e+01   9.162 1.43e-10 ***
## gndr.c.gmc  3.829e-01  7.316e-03 6.302e+04  52.333  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##            (Intr)
## gndr.c.gmc -0.004
```

``` r
getFE(mod1_alt,round=3)
```

```
##              Est.    SE        df      t p    LL    UL
## (Intercept) 0.424 0.046    32.884  9.162 0 0.330 0.518
## gndr.c.gmc  0.383 0.007 63016.617 52.333 0 0.369 0.397
```

``` r
getVC(mod1_alt)
```

```
##        grp        var1 var2 sdcor vcov
## 1    cntry (Intercept) <NA>  0.26 0.07
## 2 Residual        <NA> <NA>  1.00 0.99
```

``` r
r2mlm(mod1_alt,bargraph = F)
```

```
## $Decompositions
##                      total     within between
## fixed, within   0.03312297 0.03537171      NA
## fixed, between  0.00000000         NA       0
## slope variation 0.00000000 0.00000000      NA
## mean variation  0.06357473         NA       1
## sigma2          0.90330230 0.96462829      NA
## 
## $R2s
##          total     within between
## f1  0.03312297 0.03537171      NA
## f2  0.00000000         NA       0
## v   0.00000000 0.00000000      NA
## m   0.06357473         NA       1
## f   0.03312297         NA      NA
## fv  0.03312297 0.03537171      NA
## fvm 0.09669770         NA      NA
```

## mod2: Gender fixed and random effect

* Include random effect correlation by default


``` r
mod2<-lmer(FM.z~gndr.c+(gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method [lmerModLmerTest]
## Formula: FM.z ~ gndr.c + (gndr.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##  176204.0  176258.3  -88096.0  176192.0     63041 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.6035 -0.5866  0.0247  0.5968  9.1316 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr 
##  cntry    (Intercept) 0.0725   0.2693        
##           gndr.c      0.0121   0.1100   -0.32
##  Residual             0.9906   0.9953        
## Number of obs: 63047, groups:  cntry, 33
## 
## Fixed effects:
##             Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)  0.43153    0.04707 32.88411   9.168 1.41e-10 ***
## gndr.c       0.36866    0.02085 31.15701  17.684  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##        (Intr)
## gndr.c -0.296
```

``` r
getFE(mod2,round=3)
```

```
##              Est.    SE     df      t p    LL    UL
## (Intercept) 0.432 0.047 32.884  9.168 0 0.336 0.527
## gndr.c      0.369 0.021 31.157 17.684 0 0.326 0.411
```

``` r
getVC(mod2)
```

```
##        grp        var1   var2 sdcor  vcov
## 1    cntry (Intercept)   <NA>  0.27  0.07
## 2    cntry      gndr.c   <NA>  0.11  0.01
## 3    cntry (Intercept) gndr.c -0.32 -0.01
## 4 Residual        <NA>   <NA>  1.00  0.99
```

``` r
r2mlm(mod2,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.030843151
## slope variation 0.002746729
## mean variation  0.066181156
## sigma2          0.900228964
## 
## $R2s
##           total
## f   0.030843151
## v   0.002746729
## m   0.066181156
## fv  0.033589880
## fvm 0.099771036
```

``` r
anova(mod1,mod2)
```

```
## Data: diff_dat
## Models:
## mod1: FM.z ~ gndr.c + (1 | cntry)
## mod2: FM.z ~ gndr.c + (gndr.c | cntry)
##      npar    AIC    BIC logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod1    4 176345 176381 -88168    176337                         
## mod2    6 176204 176258 -88096    176192 144.99  2  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

``` r
lvl2_var_cond_lvl1(mod2,lvl1.var = "gndr.c",lvl1.values = c(0.5,-0.5))
```

```
##   lvl1.value lvl2.cond.var lvl2.cond.sd
## 1        0.5    0.06593716    0.2567823
## 2       -0.5    0.08510945    0.2917352
```

* Test for random effect correlation


``` r
mod2_norecov<-lmer(FM.z~gndr.c+(gndr.c||cntry),data=diff_dat,REML=F,weights = pspwght,
                   control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2_norecov)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method [lmerModLmerTest]
## Formula: FM.z ~ gndr.c + (gndr.c || cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##  176204.9  176250.1  -88097.4  176194.9     63042 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.6171 -0.5866  0.0249  0.5965  9.1180 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  cntry    (Intercept) 0.07249  0.2692  
##  cntry.1  gndr.c      0.01222  0.1106  
##  Residual             0.99060  0.9953  
## Number of obs: 63047, groups:  cntry, 33
## 
## Fixed effects:
##             Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)  0.43133    0.04707 32.88524   9.164 1.42e-10 ***
## gndr.c       0.36830    0.02095 30.99607  17.582  < 2e-16 ***
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
##              Est.    SE     df      t p    LL    UL
## (Intercept) 0.431 0.047 32.885  9.164 0 0.336 0.527
## gndr.c      0.368 0.021 30.996 17.582 0 0.326 0.411
```

``` r
getVC(mod2_norecov)
```

```
##        grp        var1 var2 sdcor vcov
## 1    cntry (Intercept) <NA>  0.27 0.07
## 2  cntry.1      gndr.c <NA>  0.11 0.01
## 3 Residual        <NA> <NA>  1.00 0.99
```

``` r
anova(mod2_norecov,mod2)
```

```
## Data: diff_dat
## Models:
## mod2_norecov: FM.z ~ gndr.c + (gndr.c || cntry)
## mod2: FM.z ~ gndr.c + (gndr.c | cntry)
##              npar    AIC    BIC logLik -2*log(L)  Chisq Df Pr(>Chisq)  
## mod2_norecov    5 176205 176250 -88097    176195                       
## mod2            6 176204 176258 -88096    176192 2.9037  1    0.08838 .
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```


* Alternative model with country-mean centered gender


``` r
mod2_alt<-lmer(FM.z~gndr.c.gmc+(gndr.c.gmc|cntry),data=diff_dat_gmc,
               REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2_alt)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method [lmerModLmerTest]
## Formula: FM.z ~ gndr.c.gmc + (gndr.c.gmc | cntry)
##    Data: diff_dat_gmc
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##  176203.1  176257.4  -88095.5  176191.1     63041 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.6055 -0.5867  0.0247  0.5968  9.1296 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr 
##  cntry    (Intercept) 0.07056  0.2656        
##           gndr.c.gmc  0.01210  0.1100   -0.32
##  Residual             0.99060  0.9953        
## Number of obs: 63047, groups:  cntry, 33
## 
## Fixed effects:
##             Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)  0.42465    0.04644 32.88278   9.143  1.5e-10 ***
## gndr.c.gmc   0.36876    0.02085 31.13290  17.690  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##            (Intr)
## gndr.c.gmc -0.296
```

``` r
getFE(mod2_alt,round=3)
```

```
##              Est.    SE     df      t p    LL    UL
## (Intercept) 0.425 0.046 32.883  9.143 0 0.330 0.519
## gndr.c.gmc  0.369 0.021 31.133 17.690 0 0.326 0.411
```

``` r
getVC(mod2_alt)
```

```
##        grp        var1       var2 sdcor  vcov
## 1    cntry (Intercept)       <NA>  0.27  0.07
## 2    cntry  gndr.c.gmc       <NA>  0.11  0.01
## 3    cntry (Intercept) gndr.c.gmc -0.32 -0.01
## 4 Residual        <NA>       <NA>  1.00  0.99
```

``` r
r2mlm(mod2_alt,bargraph = F)
```

```
## $Decompositions
##                       total      within between
## fixed, within   0.030790237 0.032904818      NA
## fixed, between  0.000000000          NA       0
## slope variation 0.002739748 0.002927906      NA
## mean variation  0.064263559          NA       1
## sigma2          0.902206455 0.964167276      NA
## 
## $R2s
##           total      within between
## f1  0.030790237 0.032904818      NA
## f2  0.000000000          NA       0
## v   0.002739748 0.002927906      NA
## m   0.064263559          NA       1
## f   0.030790237          NA      NA
## fv  0.033529986 0.035832724      NA
## fvm 0.097793545          NA      NA
```

``` r
anova(mod1_alt,mod2_alt)
```

```
## Data: diff_dat_gmc
## Models:
## mod1_alt: FM.z ~ gndr.c.gmc + (1 | cntry)
## mod2_alt: FM.z ~ gndr.c.gmc + (gndr.c.gmc | cntry)
##          npar    AIC    BIC logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod1_alt    4 176344 176380 -88168    176336                         
## mod2_alt    6 176203 176257 -88096    176191 144.63  2  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## mod2 with Gender-equality index (GEI)


``` r
mod2_GEI<-lmer(FM.z~gndr.c+gei.z.cm+gndr.c:gei.z.cm+
                 (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2_GEI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method [lmerModLmerTest]
## Formula: FM.z ~ gndr.c + gei.z.cm + gndr.c:gei.z.cm + (gndr.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##  171248.7  171320.9  -85616.4  171232.7     61568 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.6538 -0.5896  0.0247  0.5985  9.1735 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr
##  cntry    (Intercept) 0.057128 0.2390       
##           gndr.c      0.006069 0.0779   0.05
##  Residual             0.979095 0.9895       
## Number of obs: 61576, groups:  cntry, 32
## 
## Fixed effects:
##                 Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)      0.42372    0.04248 31.88107   9.975 2.51e-11 ***
## gndr.c           0.37147    0.01602 31.08402  23.192  < 2e-16 ***
## gei.z.cm        -0.12727    0.04319 31.96551  -2.947  0.00595 ** 
## gndr.c:gei.z.cm  0.08098    0.01668 34.75814   4.855 2.52e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c g.z.cm
## gndr.c       0.043              
## gei.z.cm    -0.002 -0.001       
## gndr.c:g.z. -0.001 -0.048  0.043
```

``` r
getFE(mod2_GEI,round=3)
```

```
##                   Est.    SE     df      t     p     LL     UL
## (Intercept)      0.424 0.042 31.881  9.975 0.000  0.337  0.510
## gndr.c           0.371 0.016 31.084 23.192 0.000  0.339  0.404
## gei.z.cm        -0.127 0.043 31.966 -2.947 0.006 -0.215 -0.039
## gndr.c:gei.z.cm  0.081 0.017 34.758  4.855 0.000  0.047  0.115
```

``` r
getVC(mod2_GEI)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.24 0.06
## 2    cntry      gndr.c   <NA>  0.08 0.01
## 3    cntry (Intercept) gndr.c  0.05 0.00
## 4 Residual        <NA>   <NA>  0.99 0.98
```

``` r
r2mlm(mod2_GEI,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.045621882
## slope variation 0.001394007
## mean variation  0.052514715
## sigma2          0.900469396
## 
## $R2s
##           total
## f   0.045621882
## v   0.001394007
## m   0.052514715
## fv  0.047015889
## fvm 0.099530604
```

``` r
#anova(mod2,mod2_GEI)
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
## Time difference of 4.111369 secs
```

``` r
round(ddsc_mod2_GEI$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.085        0.292        0.991     1.076 0.079    987.576 0.987   0.988
## 2        0.5         0.066        0.257        0.991     1.057 0.062    922.939 0.983   0.984
```

``` r
round(ddsc_mod2_GEI$descriptives,3)
```

```
##                       M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gei.z.cm
## means_y1          0.607 0.268    1.000           1.000    0.903           0.903   -0.322
## means_y1_scaled   2.143 0.946    1.000           1.000    0.903           0.903   -0.322
## means_y2          0.240 0.298    0.903           0.903    1.000           1.000   -0.565
## means_y2_scaled   0.847 1.052    0.903           0.903    1.000           1.000   -0.565
## gei.z.cm          0.000 1.000   -0.322          -0.322   -0.565          -0.565    1.000
## gei.z.cm_scaled   0.000 1.000   -0.322          -0.322   -0.565          -0.565    1.000
## diff_score        0.367 0.128   -0.008          -0.008   -0.437          -0.437    0.641
## diff_score_scaled 1.296 0.452   -0.008          -0.008   -0.437          -0.437    0.641
##                   gei.z.cm_scaled diff_score diff_score_scaled
## means_y1                   -0.322     -0.008            -0.008
## means_y1_scaled            -0.322     -0.008            -0.008
## means_y2                   -0.565     -0.437            -0.437
## means_y2_scaled            -0.565     -0.437            -0.437
## gei.z.cm                    1.000      0.641             0.641
## gei.z.cm_scaled             1.000      0.641             0.641
## diff_score                  0.641      1.000             1.000
## diff_score_scaled           0.641      1.000             1.000
```

``` r
round(ddsc_mod2_GEI$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                      -0.633 0.130 34.758  -4.855   0.000   -0.897   -0.368
## w_11                         -0.168 0.044 31.993  -3.845   0.001   -0.257   -0.079
## w_21                         -0.087 0.044 31.977  -1.957   0.059   -0.177    0.004
## r_xy1                        -0.627 0.163 31.993  -3.845   0.001   -0.959   -0.295
## r_xy2                        -0.292 0.149 31.977  -1.957   0.059   -0.595    0.012
## b_11                         -0.594 0.154 31.993  -3.845   0.001   -0.908   -0.279
## b_21                         -0.307 0.157 31.977  -1.957   0.059   -0.627    0.012
## main_effect                  -0.127 0.043 31.966  -2.947   0.006   -0.215   -0.039
## moderator_effect              0.371 0.016 31.084  23.192   0.000    0.339    0.404
## interaction                   0.081 0.017 34.758   4.855   0.000    0.047    0.115
## q_b11_b21                    -0.366    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                  -0.436    NA     NA      NA      NA       NA       NA
## cross_over_point             -4.587    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.046 0.047 31.940  -0.986   0.332   -0.142    0.049
## interaction_vs_main_bscale   -0.164 0.166 31.940  -0.986   0.332   -0.502    0.175
## interaction_vs_main_rscale   -0.124 0.151 31.929  -0.822   0.417   -0.431    0.183
## dadas                        -0.174 0.089 31.977  -1.957   0.970   -0.354    0.007
## dadas_bscale                 -0.614 0.314 31.977  -1.957   0.970   -1.253    0.025
## dadas_rscale                 -0.583 0.298 31.977  -1.957   0.970   -1.190    0.024
## abs_diff                      0.081 0.017 34.758   4.855   0.000    0.047    0.115
## abs_sum                       0.255 0.086 31.966   2.947   0.003    0.079    0.430
## abs_diff_bscale               0.287 0.059 34.758   4.855   0.000    0.167    0.406
## abs_sum_bscale                0.901 0.306 31.966   2.947   0.003    0.278    1.523
## abs_diff_rscale               0.335 0.061 35.499   5.522   0.000    0.212    0.458
## abs_sum_rscale                0.918 0.306 31.965   2.998   0.003    0.294    1.542
```

``` r
round(ddsc_mod2_GEI$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.010 -0.324  2.904  1.000  0.088
```

``` r
d_GEI<-ddsc_mod2_GEI$ddsc_sem_fit$data

ddsc_sem_GEI<-
  ddsc_sem(data=d_GEI,x = "gei.z.cm",y1="means_y2",
           y2="means_y1",q_sesoi = .10)
round(ddsc_sem_GEI$results,3)
```

```
##                                    est    se       z pvalue ci.lower ci.upper
## r_xy1_y2                        -0.641 0.136  -4.724  0.000   -0.907   -0.375
## r_xy1                           -0.565 0.146  -3.876  0.000   -0.851   -0.279
## r_xy2                           -0.322 0.167  -1.924  0.054   -0.650    0.006
## b_11                            -0.594 0.153  -3.876  0.000   -0.895   -0.294
## b_21                            -0.304 0.158  -1.924  0.054   -0.615    0.006
## b_10                             0.847 0.151   5.612  0.000    0.551    1.143
## b_20                             2.143 0.156  13.758  0.000    1.838    2.448
## res_cov_y1_y2                    0.694 0.181   3.837  0.000    0.340    1.049
## diff_b10_b20                    -1.296 0.060 -21.457  0.000   -1.414   -1.178
## diff_b11_b21                    -0.290 0.061  -4.724  0.000   -0.410   -0.170
## diff_rxy1_rxy2                  -0.243 0.065  -3.742  0.000   -0.371   -0.116
## q_b11_b21                       -0.370 0.102  -3.644  0.000   -0.569   -0.171
## q_rxy1_rxy2                     -0.307 0.083  -3.684  0.000   -0.470   -0.143
## cross_over_point                -4.471 0.969  -4.614  0.000   -6.370   -2.571
## sum_b11_b21                     -0.899 0.306  -2.942  0.003   -1.498   -0.300
## main_effect                     -0.449 0.153  -2.942  0.003   -0.749   -0.150
## interaction_vs_main_effect      -0.160 0.169  -0.943  0.346   -0.491    0.172
## diff_abs_b11_abs_b21             0.290 0.061   4.724  0.000    0.170    0.410
## abs_diff_b11_b21                 0.290 0.061   4.724  0.000    0.170    0.410
## abs_sum_b11_b21                  0.899 0.306   2.942  0.002    0.300    1.498
## dadas                           -0.609 0.317  -1.924  0.973   -1.229    0.011
## q_r_equivalence                  0.207 0.083   2.482  0.993       NA       NA
## q_b_equivalence                  0.270 0.102   2.659  0.996       NA       NA
## cross_over_point_equivalence     4.471 0.969   4.614  1.000       NA       NA
## cross_over_point_minimal_effect  4.471 0.969   4.614  0.000       NA       NA
```

``` r
round(ddsc_sem_GEI$variance_test,3)
```

```
##             est    se      z pvalue ci.lower ci.upper
## cov_y1y2  0.870 0.229  3.791  0.000    0.420    1.319
## var_y1    1.071 0.268  4.000  0.000    0.546    1.596
## var_y2    0.866 0.217  4.000  0.000    0.442    1.291
## var_diff  0.205 0.155  1.321  0.187   -0.099    0.509
## var_ratio 1.237 0.188  6.576  0.000    0.868    1.605
## cor_y1y2  0.903 0.033 27.608  0.000    0.839    0.967
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
flag_points_GEI
```

```
##    (Intercept)    gndr.c        women        men  mean_level difference cntry    gei.cm
## 1   0.46910842 0.4485451  0.244835875 0.69338097  0.46910842  0.4485451    AT 0.9106190
## 2   0.32356214 0.3808445  0.133139904 0.51398438  0.32356214  0.3808445    BE 0.9225238
## 3   0.47387200 0.3677068  0.290018620 0.65772537  0.47387200  0.3677068    BG 0.7777143
## 4   0.22591935 0.4771312 -0.012646257 0.46448497  0.22591935  0.4771312    CH 0.9528095
## 5   0.15303188 0.2372716  0.034396071 0.27166769  0.15303188  0.2372716    CY 0.7855714
## 6   0.93057045 0.4914030  0.684868947 1.17627195  0.93057045  0.4914030    CZ 0.8596190
## 7   0.20782346 0.4492190 -0.016786057 0.43243298  0.20782346  0.4492190    DE 0.9113810
## 8   0.37107531 0.5082305  0.116960053 0.62519056  0.37107531  0.5082305    DK 0.9620476
## 9   0.38115603 0.4346623  0.163824871 0.59848719  0.38115603  0.4346623    EE 0.8486190
## 10 -0.06706344 0.3145103 -0.224318604 0.09019173 -0.06706344  0.3145103    ES 0.9090000
## 11  0.10927543 0.5032436 -0.142346355 0.36089721  0.10927543  0.5032436    FI 0.9383333
## 12  0.04684201 0.5265046 -0.216410272 0.31009429  0.04684201  0.5265046    FR 0.8897143
## 13  0.24472225 0.4807966  0.004323969 0.48512054  0.24472225  0.4807966    GB 0.8453333
## 14  0.49191071 0.2018644  0.390978498 0.59284292  0.49191071  0.2018644    GR 0.8516190
## 15  0.25062695 0.3869044  0.057174746 0.44407916  0.25062695  0.3869044    HR 0.8639048
## 16  0.55793675 0.3009417  0.407465893 0.70840761  0.55793675  0.3009417    HU 0.7514762
## 17  0.42227675 0.3629228  0.240815337 0.60373817  0.42227675  0.3629228    IE 0.8667619
## 18  0.46176273 0.2728115  0.325356964 0.59816850  0.46176273  0.2728115    IL 0.8701429
## 19 -0.11461611 0.4352629 -0.332247549 0.10301532 -0.11461611  0.4352629    IS 0.9205714
## 20  0.43510299 0.3440312  0.263087404 0.60711857  0.43510299  0.3440312    IT 0.8936190
## 21  1.03984359 0.3466629  0.866512150 1.21317503  1.03984359  0.3466629    LT 0.8493810
## 22  0.89184510 0.3648590  0.709415623 1.07427458  0.89184510  0.3648590    LV 0.7867143
## 23  0.47184627 0.2235117  0.360090429 0.58360210  0.47184627  0.2235117    ME 0.8387647
## 24  0.55600859 0.4532554  0.329380896 0.78263628  0.55600859  0.4532554    NL 0.9533333
## 25  0.34824643 0.4001341  0.148179355 0.54831350  0.34824643  0.4001341    NO 0.9506667
## 26  0.39586987 0.4233198  0.184209959 0.60752977  0.39586987  0.4233198    PL 0.8590952
## 27  0.40599286 0.2625674  0.274709176 0.53727655  0.40599286  0.2625674    PT 0.8843333
## 28  0.85454708 0.1750750  0.767009567 0.94208460  0.85454708  0.1750750    RU 0.7484286
## 29  0.25763869 0.4825205  0.016378456 0.49889893  0.25763869  0.4825205    SE 0.9611429
## 30  0.56045022 0.3714961  0.374702187 0.74619825  0.56045022  0.3714961    SI 0.9077143
## 31  0.74523002 0.3195112  0.585474428 0.90498562  0.74523002  0.3195112    SK 0.8110000
## 32  0.65783058 0.1646775  0.575491842 0.74016932  0.65783058  0.1646775    TR 0.6203810
## 33  0.68026564 0.2534444  0.553543418 0.80698786  0.68026564  0.2534444    UA       NaN
```

``` r
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
  xlab("Gender Equality Index (Average 2002-2022)")+
  #ylim(c(min.y.pred,max.y.pred))+
  #xlim(c(0.60,1.00))+
  ylab("Value male-typicality (Average 2002-2022)")+
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

p1.FM.flags
```

```
## Warning in geom_text(inherit.aes = F, aes(x = 0.65, y = -0.45, label = coef_q, : All aesthetics have length 1, but the data has 682 rows.
## ℹ Please consider using `annotate()` or provide this layer with data containing a single row.
```

```
## Warning: Removed 2 rows containing missing values or values outside the scale range
## (`geom_point()`).
```

```
## Warning: Removed 2 rows containing missing values or values outside the scale range
## (`geom_line()`).
```

```
## Warning: Removed 2 rows containing missing values or values outside the scale range
## (`geom_flag()`).
```

![](Analysis_youth_files/figure-html/unnamed-chunk-22-1.png)<!-- -->

``` r
p2.FM.flags<-ggplot(p2,aes(y=yvar,x=gei.cm))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("Gender Equality Index (Average 2002-2022)")+
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

p2.FM.flags
```

```
## Warning: Removed 1 row containing missing values or values outside the scale range (`geom_flag()`).
```

![](Analysis_youth_files/figure-html/unnamed-chunk-22-2.png)<!-- -->

``` r
pflag_comb<-
  ggarrange(p1.FM.flags,p2.FM.flags,align = "v",
            ncol=1,nrow=2,heights=c(2,1))
```

```
## Warning in geom_text(inherit.aes = F, aes(x = 0.65, y = -0.45, label = coef_q, : All aesthetics have length 1, but the data has 682 rows.
## ℹ Please consider using `annotate()` or provide this layer with data containing a single row.
```

```
## Warning: Removed 2 rows containing missing values or values outside the scale range
## (`geom_point()`).
```

```
## Warning: Removed 2 rows containing missing values or values outside the scale range
## (`geom_line()`).
```

```
## Warning: Removed 2 rows containing missing values or values outside the scale range
## (`geom_flag()`).
```

```
## Warning: Removed 1 row containing missing values or values outside the scale range (`geom_flag()`).
```

``` r
pflag_comb
```

![](Analysis_youth_files/figure-html/unnamed-chunk-22-3.png)<!-- -->

``` r
png(filename = 
      "../results/GEI_flags_youth.png",
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
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method [lmerModLmerTest]
## Formula: FM.z ~ gndr.c + gggi.z.cm + gndr.c:gggi.z.cm + (gndr.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##  124896.8  124966.5  -62440.4  124880.8     44999 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.7557 -0.5934  0.0243  0.5924  6.8789 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr 
##  cntry    (Intercept) 0.064448 0.25387       
##           gndr.c      0.007986 0.08936  -0.17
##  Residual             0.972986 0.98640       
## Number of obs: 45007, groups:  cntry, 33
## 
## Fixed effects:
##                  Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)       0.45372    0.04451 32.80692  10.193 1.07e-11 ***
## gndr.c            0.35972    0.01827 29.19456  19.693  < 2e-16 ***
## gggi.z.cm        -0.16007    0.04528 33.03273  -3.535 0.001231 ** 
## gndr.c:gggi.z.cm  0.07756    0.01934 33.93245   4.010 0.000315 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c ggg.z.
## gndr.c      -0.143              
## gggi.z.cm   -0.001 -0.002       
## gndr.c:gg.. -0.002 -0.019 -0.137
```

``` r
getFE(mod2_GGGI,round=3)
```

```
##                    Est.    SE     df      t     p     LL     UL
## (Intercept)       0.454 0.045 32.807 10.193 0.000  0.363  0.544
## gndr.c            0.360 0.018 29.195 19.693 0.000  0.322  0.397
## gggi.z.cm        -0.160 0.045 33.033 -3.535 0.001 -0.252 -0.068
## gndr.c:gggi.z.cm  0.078 0.019 33.932  4.010 0.000  0.038  0.117
```

``` r
getVC(mod2_GGGI)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.25 0.06
## 2    cntry      gndr.c   <NA>  0.09 0.01
## 3    cntry (Intercept) gndr.c -0.17 0.00
## 4 Residual        <NA>   <NA>  0.99 0.97
```

``` r
r2mlm(mod2_GGGI,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.048780368
## slope variation 0.001824793
## mean variation  0.059091651
## sigma2          0.890303188
## 
## $R2s
##           total
## f   0.048780368
## v   0.001824793
## m   0.059091651
## fv  0.050605161
## fvm 0.109696812
```

``` r
#anova(mod2,mod2_GGGI)
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
## Time difference of 3.997287 secs
```

``` r
round(ddsc_mod2_GGGI$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.085        0.292        0.991     1.076 0.079    987.576 0.987   0.988
## 2        0.5         0.066        0.257        0.991     1.057 0.062    922.939 0.983   0.984
```

``` r
round(ddsc_mod2_GGGI$descriptives,3)
```

```
##                       M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gggi.z.cm
## means_y1          0.624 0.286    1.000           1.000    0.870           0.870    -0.404
## means_y1_scaled   1.978 0.907    1.000           1.000    0.870           0.870    -0.404
## means_y2          0.280 0.342    0.870           0.870    1.000           1.000    -0.588
## means_y2_scaled   0.886 1.085    0.870           0.870    1.000           1.000    -0.588
## gggi.z.cm         0.000 1.000   -0.404          -0.404   -0.588          -0.588     1.000
## gggi.z.cm_scaled  0.000 1.000   -0.404          -0.404   -0.588          -0.588     1.000
## diff_score        0.344 0.169   -0.068          -0.068   -0.551          -0.551     0.508
## diff_score_scaled 1.091 0.536   -0.068          -0.068   -0.551          -0.551     0.508
##                   gggi.z.cm_scaled diff_score diff_score_scaled
## means_y1                    -0.404     -0.068            -0.068
## means_y1_scaled             -0.404     -0.068            -0.068
## means_y2                    -0.588     -0.551            -0.551
## means_y2_scaled             -0.588     -0.551            -0.551
## gggi.z.cm                    1.000      0.508             0.508
## gggi.z.cm_scaled             1.000      0.508             0.508
## diff_score                   0.508      1.000             1.000
## diff_score_scaled            0.508      1.000             1.000
```

``` r
round(ddsc_mod2_GGGI$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                      -0.458 0.114 33.932  -4.010   0.000   -0.691   -0.226
## w_11                         -0.199 0.048 33.051  -4.179   0.000   -0.296   -0.102
## w_21                         -0.121 0.045 33.165  -2.696   0.011   -0.213   -0.030
## r_xy1                        -0.694 0.166 33.051  -4.179   0.000   -1.033   -0.356
## r_xy2                        -0.354 0.131 33.165  -2.696   0.011   -0.621   -0.087
## b_11                         -0.633 0.151 33.051  -4.179   0.000   -0.940   -0.325
## b_21                         -0.386 0.143 33.165  -2.696   0.011   -0.677   -0.095
## main_effect                  -0.160 0.045 33.033  -3.535   0.001   -0.252   -0.068
## moderator_effect              0.360 0.018 29.195  19.693   0.000    0.322    0.397
## interaction                   0.078 0.019 33.932   4.010   0.000    0.038    0.117
## q_b11_b21                    -0.339    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                  -0.486    NA     NA      NA      NA       NA       NA
## cross_over_point             -4.638    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.083 0.047 33.374  -1.766   0.087   -0.178    0.013
## interaction_vs_main_bscale   -0.262 0.149 33.374  -1.766   0.087   -0.565    0.040
## interaction_vs_main_rscale   -0.184 0.126 33.440  -1.465   0.152   -0.440    0.071
## dadas                        -0.243 0.090 33.165  -2.696   0.995   -0.426   -0.060
## dadas_bscale                 -0.772 0.286 33.165  -2.696   0.995   -1.354   -0.189
## dadas_rscale                 -0.708 0.263 33.165  -2.696   0.995   -1.243   -0.174
## abs_diff                      0.078 0.019 33.932   4.010   0.000    0.038    0.117
## abs_sum                       0.320 0.091 33.033   3.535   0.001    0.136    0.504
## abs_diff_bscale               0.247 0.062 33.932   4.010   0.000    0.122    0.372
## abs_sum_bscale                1.018 0.288 33.033   3.535   0.001    0.432    1.604
## abs_diff_rscale               0.340 0.070 33.931   4.832   0.000    0.197    0.483
## abs_sum_rscale                1.049 0.291 33.028   3.601   0.001    0.456    1.641
```

``` r
round(ddsc_mod2_GGGI$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.010 -0.324  2.904  1.000  0.088
```

``` r
d_GGGI<-ddsc_mod2_GGGI$ddsc_sem_fit$data

ddsc_sem_GGGI<-
  ddsc_sem(data=d_GGGI,x = "gggi.z.cm",y1="means_y2",
           y2="means_y1",q_sesoi = .10)
round(ddsc_sem_GGGI$results,3)
```

```
##                                    est    se       z pvalue ci.lower ci.upper
## r_xy1_y2                        -0.508 0.150  -3.385  0.001   -0.802   -0.214
## r_xy1                           -0.588 0.141  -4.181  0.000   -0.864   -0.313
## r_xy2                           -0.404 0.159  -2.534  0.011   -0.716   -0.091
## b_11                            -0.638 0.153  -4.181  0.000   -0.938   -0.339
## b_21                            -0.366 0.144  -2.534  0.011   -0.649   -0.083
## b_10                             0.886 0.150   5.895  0.000    0.592    1.181
## b_20                             1.978 0.142  13.900  0.000    1.699    2.257
## res_cov_y1_y2                    0.604 0.162   3.733  0.000    0.287    0.921
## diff_b10_b20                    -1.091 0.079 -13.787  0.000   -1.247   -0.936
## diff_b11_b21                    -0.272 0.080  -3.385  0.001   -0.430   -0.115
## diff_rxy1_rxy2                  -0.185 0.083  -2.234  0.025   -0.347   -0.023
## q_b11_b21                       -0.371 0.144  -2.580  0.010   -0.653   -0.089
## q_rxy1_rxy2                     -0.247 0.112  -2.211  0.027   -0.466   -0.028
## cross_over_point                -4.011 1.220  -3.287  0.001   -6.402   -1.619
## sum_b11_b21                     -1.005 0.286  -3.510  0.000   -1.566   -0.444
## main_effect                     -0.502 0.143  -3.510  0.000   -0.783   -0.222
## interaction_vs_main_effect      -0.230 0.157  -1.470  0.142   -0.537    0.077
## diff_abs_b11_abs_b21             0.272 0.080   3.385  0.001    0.115    0.430
## abs_diff_b11_b21                 0.272 0.080   3.385  0.000    0.115    0.430
## abs_sum_b11_b21                  1.005 0.286   3.510  0.000    0.444    1.566
## dadas                           -0.732 0.289  -2.534  0.994   -1.299   -0.166
## q_r_equivalence                  0.147 0.112   1.317  0.906       NA       NA
## q_b_equivalence                  0.271 0.144   1.885  0.970       NA       NA
## cross_over_point_equivalence     4.011 1.220   3.287  0.999       NA       NA
## cross_over_point_minimal_effect  4.011 1.220   3.287  0.001       NA       NA
```

``` r
round(ddsc_sem_GGGI$variance_test,3)
```

```
##             est    se      z pvalue ci.lower ci.upper
## cov_y1y2  0.830 0.220  3.771  0.000    0.399    1.262
## var_y1    1.141 0.281  4.062  0.000    0.591    1.692
## var_y2    0.798 0.196  4.062  0.000    0.413    1.183
## var_diff  0.343 0.184  1.861  0.063   -0.018    0.704
## var_ratio 1.430 0.245  5.826  0.000    0.949    1.911
## cor_y1y2  0.870 0.042 20.566  0.000    0.787    0.953
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
flag_points_GGGI
```

```
##    (Intercept)    gndr.c        women        men  mean_level difference cntry   gggi.cm
## 1   0.46910842 0.4485451  0.244835875 0.69338097  0.46910842  0.4485451    AT 0.7280529
## 2   0.32356214 0.3808445  0.133139904 0.51398438  0.32356214  0.3808445    BE 0.7518824
## 3   0.47387200 0.3677068  0.290018620 0.65772537  0.47387200  0.3677068    BG 0.7206176
## 4   0.22591935 0.4771312 -0.012646257 0.46448497  0.22591935  0.4771312    CH 0.7609706
## 5   0.15303188 0.2372716  0.034396071 0.27166769  0.15303188  0.2372716    CY 0.6752824
## 6   0.93057045 0.4914030  0.684868947 1.17627195  0.93057045  0.4914030    CZ 0.6859235
## 7   0.20782346 0.4492190 -0.016786057 0.43243298  0.20782346  0.4492190    DE 0.7709588
## 8   0.37107531 0.5082305  0.116960053 0.62519056  0.37107531  0.5082305    DK 0.7700353
## 9   0.38115603 0.4346623  0.163824871 0.59848719  0.38115603  0.4346623    EE 0.7217941
## 10 -0.06706344 0.3145103 -0.224318604 0.09019173 -0.06706344  0.3145103    ES 0.7512882
## 11  0.10927543 0.5032436 -0.142346355 0.36089721  0.10927543  0.5032436    FI 0.8350824
## 12  0.04684201 0.5265046 -0.216410272 0.31009429  0.04684201  0.5265046    FR 0.7387118
## 13  0.24472225 0.4807966  0.004323969 0.48512054  0.24472225  0.4807966    GB 0.7554824
## 14  0.49191071 0.2018644  0.390978498 0.59284292  0.49191071  0.2018644    GR 0.6819118
## 15  0.25062695 0.3869044  0.057174746 0.44407916  0.25062695  0.3869044    HR 0.7096750
## 16  0.55793675 0.3009417  0.407465893 0.70840761  0.55793675  0.3009417    HU 0.6774059
## 17  0.42227675 0.3629228  0.240815337 0.60373817  0.42227675  0.3629228    IE 0.7819647
## 18  0.46176273 0.2728115  0.325356964 0.59816850  0.46176273  0.2728115    IL 0.7066000
## 19 -0.11461611 0.4352629 -0.332247549 0.10301532 -0.11461611  0.4352629    IS 0.8571471
## 20  0.43510299 0.3440312  0.263087404 0.60711857  0.43510299  0.3440312    IT 0.6919647
## 21  1.03984359 0.3466629  0.866512150 1.21317503  1.03984359  0.3466629    LT 0.7406000
## 22  0.89184510 0.3648590  0.709415623 1.07427458  0.89184510  0.3648590    LV 0.7554529
## 23  0.47184627 0.2235117  0.360090429 0.58360210  0.47184627  0.2235117    ME 0.7056111
## 24  0.55600859 0.4532554  0.329380896 0.78263628  0.55600859  0.4532554    NL 0.7530588
## 25  0.34824643 0.4001341  0.148179355 0.54831350  0.34824643  0.4001341    NO 0.8366824
## 26  0.39586987 0.4233198  0.184209959 0.60752977  0.39586987  0.4233198    PL 0.7085176
## 27  0.40599286 0.2625674  0.274709176 0.53727655  0.40599286  0.2625674    PT 0.7263294
## 28  0.85454708 0.1750750  0.767009567 0.94208460  0.85454708  0.1750750    RU 0.6969200
## 29  0.25763869 0.4825205  0.016378456 0.49889893  0.25763869  0.4825205    SE 0.8158706
## 30  0.56045022 0.3714961  0.374702187 0.74619825  0.56045022  0.3714961    SI 0.7348235
## 31  0.74523002 0.3195112  0.585474428 0.90498562  0.74523002  0.3195112    SK 0.6903294
## 32  0.65783058 0.1646775  0.575491842 0.74016932  0.65783058  0.1646775    TR 0.6111941
## 33  0.68026564 0.2534444  0.553543418 0.80698786  0.68026564  0.2534444    UA 0.6980647
```

``` r
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
  xlab("Global Gender Gap Index (Average 2002-2022)")+
  #ylim(c(min.y.pred,max.y.pred))+
  #xlim(c(0.60,1.00))+
  ylab("Value male-typicality (Average 2002-2022)")+
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

p1.FM.flags
```

```
## Warning in geom_text(inherit.aes = F, aes(x = 0.61, y = -0.45, label = coef_q, : All aesthetics have length 1, but the data has 502 rows.
## ℹ Please consider using `annotate()` or provide this layer with data containing a single row.
```

![](Analysis_youth_files/figure-html/unnamed-chunk-25-1.png)<!-- -->

``` r
p2.FM.flags<-ggplot(p2,aes(y=yvar,x=gggi.cm))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("Global Gender Gap Index (Average 2002-2022)")+
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

p2.FM.flags
```

![](Analysis_youth_files/figure-html/unnamed-chunk-25-2.png)<!-- -->

``` r
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

![](Analysis_youth_files/figure-html/unnamed-chunk-25-3.png)<!-- -->

``` r
png(filename = 
      "../results/GGGI_flags_new_youth.png",
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
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method [lmerModLmerTest]
## Formula: FM.z ~ gndr.c + gdi.z.cm + gndr.c:gdi.z.cm + (gndr.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##  176203.0  176275.4  -88093.5  176187.0     63039 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.6094 -0.5865  0.0248  0.5969  9.1258 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr 
##  cntry    (Intercept) 0.06520  0.2554        
##           gndr.c      0.01224  0.1106   -0.37
##  Residual             0.99059  0.9953        
## Number of obs: 63047, groups:  cntry, 33
## 
## Fixed effects:
##                 Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)      0.43149    0.04466 32.91465   9.661 3.91e-11 ***
## gndr.c           0.36853    0.02094 31.49110  17.597  < 2e-16 ***
## gdi.z.cm         0.08695    0.04542 33.10364   1.915   0.0642 .  
## gndr.c:gdi.z.cm  0.01273    0.02177 34.30683   0.585   0.5627    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c gd.z.c
## gndr.c      -0.341              
## gdi.z.cm     0.000 -0.001       
## gndr.c:gd.. -0.001 -0.006 -0.332
```

``` r
getFE(mod2_GDI,round=3)
```

```
##                  Est.    SE     df      t     p     LL    UL
## (Intercept)     0.431 0.045 32.915  9.661 0.000  0.341 0.522
## gndr.c          0.369 0.021 31.491 17.597 0.000  0.326 0.411
## gdi.z.cm        0.087 0.045 33.104  1.915 0.064 -0.005 0.179
## gndr.c:gdi.z.cm 0.013 0.022 34.307  0.585 0.563 -0.032 0.057
```

``` r
getVC(mod2_GDI)
```

```
##        grp        var1   var2 sdcor  vcov
## 1    cntry (Intercept)   <NA>  0.26  0.07
## 2    cntry      gndr.c   <NA>  0.11  0.01
## 3    cntry (Intercept) gndr.c -0.37 -0.01
## 4 Residual        <NA>   <NA>  1.00  0.99
```

``` r
r2mlm(mod2_GDI,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.036043088
## slope variation 0.002781442
## mean variation  0.059667223
## sigma2          0.901508246
## 
## $R2s
##           total
## f   0.036043088
## v   0.002781442
## m   0.059667223
## fv  0.038824530
## fvm 0.098491754
```

``` r
#anova(mod2,mod2_GDI)
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
## Time difference of 3.901382 secs
```

``` r
round(ddsc_mod2_GDI$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.085        0.292        0.991     1.076 0.079    987.576 0.987   0.988
## 2        0.5         0.066        0.257        0.991     1.057 0.062    922.939 0.983   0.984
```

``` r
round(ddsc_mod2_GDI$descriptives,3)
```

```
##                       M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gdi.z.cm
## means_y1          0.612 0.266    1.000           1.000    0.903           0.903    0.368
## means_y1_scaled   2.169 0.940    1.000           1.000    0.903           0.903    0.368
## means_y2          0.250 0.298    0.903           0.903    1.000           1.000    0.260
## means_y2_scaled   0.884 1.056    0.903           0.903    1.000           1.000    0.260
## gdi.z.cm          0.000 1.000    0.368           0.368    0.260           0.260    1.000
## gdi.z.cm_scaled   0.000 1.000    0.368           0.368    0.260           0.260    1.000
## diff_score        0.363 0.128   -0.030          -0.030   -0.456          -0.456    0.158
## diff_score_scaled 1.285 0.453   -0.030          -0.030   -0.456          -0.456    0.158
##                   gdi.z.cm_scaled diff_score diff_score_scaled
## means_y1                    0.368     -0.030            -0.030
## means_y1_scaled             0.368     -0.030            -0.030
## means_y2                    0.260     -0.456            -0.456
## means_y2_scaled             0.260     -0.456            -0.456
## gdi.z.cm                    1.000      0.158             0.158
## gdi.z.cm_scaled             1.000      0.158             0.158
## diff_score                  0.158      1.000             1.000
## diff_score_scaled           0.158      1.000             1.000
```

``` r
round(ddsc_mod2_GDI$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                      -0.099 0.170 34.307  -0.585   0.563   -0.445    0.246
## w_11                          0.081 0.050 33.124   1.609   0.117   -0.021    0.183
## w_21                          0.093 0.043 33.282   2.168   0.037    0.006    0.181
## r_xy1                         0.303 0.189 33.124   1.609   0.117   -0.080    0.687
## r_xy2                         0.313 0.144 33.282   2.168   0.037    0.019    0.606
## b_11                          0.286 0.178 33.124   1.609   0.117   -0.076    0.647
## b_21                          0.331 0.153 33.282   2.168   0.037    0.020    0.642
## main_effect                   0.087 0.045 33.104   1.915   0.064   -0.005    0.179
## moderator_effect              0.369 0.021 31.491  17.597   0.000    0.326    0.411
## interaction                   0.013 0.022 34.307   0.585   0.563   -0.032    0.057
## q_b11_b21                    -0.050    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                  -0.010    NA     NA      NA      NA       NA       NA
## cross_over_point            -28.958    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.074 0.057 33.238  -1.313   0.198   -0.189    0.041
## interaction_vs_main_bscale   -0.263 0.200 33.238  -1.313   0.198   -0.671    0.144
## interaction_vs_main_rscale   -0.299 0.220 33.218  -1.357   0.184   -0.747    0.149
## dadas                        -0.161 0.100 33.124  -1.609   0.941   -0.365    0.043
## dadas_bscale                 -0.572 0.355 33.124  -1.609   0.941   -1.295    0.151
## dadas_rscale                 -0.607 0.377 33.124  -1.609   0.941   -1.374    0.161
## abs_diff                      0.013 0.022 34.307   0.585   0.281   -0.032    0.057
## abs_sum                       0.174 0.091 33.104   1.915   0.032   -0.011    0.359
## abs_diff_bscale               0.045 0.077 34.307   0.585   0.281   -0.112    0.202
## abs_sum_bscale                0.617 0.322 33.104   1.915   0.032   -0.039    1.272
## abs_diff_rscale               0.009 0.086 34.052   0.110   0.457   -0.164    0.183
## abs_sum_rscale                0.616 0.325 33.100   1.897   0.033   -0.044    1.277
```

``` r
round(ddsc_mod2_GDI$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.010 -0.324  2.904  1.000  0.088
```

``` r
d_GDI<-ddsc_mod2_GDI$ddsc_sem_fit$data

ddsc_sem_GDI<-
  ddsc_sem(data=d_GDI,x = "gdi.z.cm",y1="means_y2",
           y2="means_y1",q_sesoi = .10)
round(ddsc_sem_GDI$results,3)
```

```
##                                     est     se       z pvalue ci.lower ci.upper
## r_xy1_y2                         -0.158  0.172  -0.916  0.360   -0.494    0.179
## r_xy1                             0.260  0.168   1.550  0.121   -0.069    0.590
## r_xy2                             0.368  0.162   2.277  0.023    0.051    0.686
## b_11                              0.275  0.178   1.550  0.121   -0.073    0.623
## b_21                              0.347  0.152   2.277  0.023    0.048    0.645
## b_10                              0.884  0.175   5.055  0.000    0.541    1.226
## b_20                              2.169  0.150  14.471  0.000    1.875    2.462
## res_cov_y1_y2                     0.778  0.202   3.842  0.000    0.381    1.174
## diff_b10_b20                     -1.285  0.077 -16.751  0.000   -1.435   -1.135
## diff_b11_b21                     -0.071  0.078  -0.916  0.360   -0.224    0.081
## diff_rxy1_rxy2                   -0.108  0.074  -1.455  0.146   -0.253    0.037
## q_b11_b21                        -0.079  0.084  -0.942  0.346   -0.244    0.085
## q_rxy1_rxy2                      -0.120  0.083  -1.452  0.147   -0.282    0.042
## cross_over_point                -18.003 19.677  -0.915  0.360  -56.569   20.564
## sum_b11_b21                       0.622  0.321   1.934  0.053   -0.008    1.252
## main_effect                       0.311  0.161   1.934  0.053   -0.004    0.626
## interaction_vs_main_effect       -0.239  0.201  -1.194  0.233   -0.633    0.154
## diff_abs_b11_abs_b21             -0.071  0.078  -0.916  0.360   -0.224    0.081
## abs_diff_b11_b21                  0.071  0.078   0.916  0.180   -0.081    0.224
## abs_sum_b11_b21                   0.622  0.321   1.934  0.027   -0.008    1.252
## dadas                            -0.550  0.355  -1.550  0.939   -1.246    0.146
## q_r_equivalence                   0.020  0.083   0.242  0.596       NA       NA
## q_b_equivalence                  -0.021  0.084  -0.249  0.401       NA       NA
## cross_over_point_equivalence     18.003 19.677   0.915  0.820       NA       NA
## cross_over_point_minimal_effect  18.003 19.677   0.915  0.180       NA       NA
```

``` r
round(ddsc_sem_GDI$variance_test,3)
```

```
##             est    se      z pvalue ci.lower ci.upper
## cov_y1y2  0.870 0.226  3.851  0.000    0.427    1.313
## var_y1    1.082 0.266  4.062  0.000    0.560    1.604
## var_y2    0.858 0.211  4.062  0.000    0.444    1.271
## var_diff  0.224 0.154  1.455  0.146   -0.078    0.526
## var_ratio 1.261 0.188  6.697  0.000    0.892    1.631
## cor_y1y2  0.903 0.032 28.214  0.000    0.841    0.966
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
flag_points_GDI
```

```
##    (Intercept)    gndr.c        women        men  mean_level difference cntry    gdi.cm
## 1   0.46910842 0.4485451  0.244835875 0.69338097  0.46910842  0.4485451    AT 0.9701905
## 2   0.32356214 0.3808445  0.133139904 0.51398438  0.32356214  0.3808445    BE 0.9694286
## 3   0.47387200 0.3677068  0.290018620 0.65772537  0.47387200  0.3677068    BG 0.9954286
## 4   0.22591935 0.4771312 -0.012646257 0.46448497  0.22591935  0.4771312    CH 0.9627143
## 5   0.15303188 0.2372716  0.034396071 0.27166769  0.15303188  0.2372716    CY 0.9718571
## 6   0.93057045 0.4914030  0.684868947 1.17627195  0.93057045  0.4914030    CZ 0.9764286
## 7   0.20782346 0.4492190 -0.016786057 0.43243298  0.20782346  0.4492190    DE 0.9569048
## 8   0.37107531 0.5082305  0.116960053 0.62519056  0.37107531  0.5082305    DK 0.9871905
## 9   0.38115603 0.4346623  0.163824871 0.59848719  0.38115603  0.4346623    EE 1.0277619
## 10 -0.06706344 0.3145103 -0.224318604 0.09019173 -0.06706344  0.3145103    ES 0.9820476
## 11  0.10927543 0.5032436 -0.142346355 0.36089721  0.10927543  0.5032436    FI 0.9954762
## 12  0.04684201 0.5265046 -0.216410272 0.31009429  0.04684201  0.5265046    FR 0.9865714
## 13  0.24472225 0.4807966  0.004323969 0.48512054  0.24472225  0.4807966    GB 0.9704762
## 14  0.49191071 0.2018644  0.390978498 0.59284292  0.49191071  0.2018644    GR 0.9658095
## 15  0.25062695 0.3869044  0.057174746 0.44407916  0.25062695  0.3869044    HR 0.9860952
## 16  0.55793675 0.3009417  0.407465893 0.70840761  0.55793675  0.3009417    HU 0.9923810
## 17  0.42227675 0.3629228  0.240815337 0.60373817  0.42227675  0.3629228    IE 0.9770476
## 18  0.46176273 0.2728115  0.325356964 0.59816850  0.46176273  0.2728115    IL 0.9712857
## 19 -0.11461611 0.4352629 -0.332247549 0.10301532 -0.11461611  0.4352629    IS 0.9726667
## 20  0.43510299 0.3440312  0.263087404 0.60711857  0.43510299  0.3440312    IT 0.9679524
## 21  1.03984359 0.3466629  0.866512150 1.21317503  1.03984359  0.3466629    LT 1.0292857
## 22  0.89184510 0.3648590  0.709415623 1.07427458  0.89184510  0.3648590    LV 1.0351905
## 23  0.47184627 0.2235117  0.360090429 0.58360210  0.47184627  0.2235117    ME 0.9638333
## 24  0.55600859 0.4532554  0.329380896 0.78263628  0.55600859  0.4532554    NL 0.9581905
## 25  0.34824643 0.4001341  0.148179355 0.54831350  0.34824643  0.4001341    NO 0.9973810
## 26  0.39586987 0.4233198  0.184209959 0.60752977  0.39586987  0.4233198    PL 1.0078571
## 27  0.40599286 0.2625674  0.274709176 0.53727655  0.40599286  0.2625674    PT 0.9940476
## 28  0.85454708 0.1750750  0.767009567 0.94208460  0.85454708  0.1750750    RU 1.0287500
## 29  0.25763869 0.4825205  0.016378456 0.49889893  0.25763869  0.4825205    SE 0.9892381
## 30  0.56045022 0.3714961  0.374702187 0.74619825  0.56045022  0.3714961    SI 1.0018095
## 31  0.74523002 0.3195112  0.585474428 0.90498562  0.74523002  0.3195112    SK 0.9939524
## 32  0.65783058 0.1646775  0.575491842 0.74016932  0.65783058  0.1646775    TR 0.9031905
## 33  0.68026564 0.2534444  0.553543418 0.80698786  0.68026564  0.2534444    UA 1.0176667
```

``` r
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
  xlab("Gender Development Index (Average 2002-2022)")+
  #ylim(c(min.y.pred,max.y.pred))+
  #xlim(c(0.60,1.00))+
  ylab("Value male-typicality (Average 2002-2022)")+
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

p1.FM.flags
```

```
## Warning in geom_text(inherit.aes = F, aes(x = 0.9, y = -0.45, label = coef_q, : All aesthetics have length 1, but the data has 282 rows.
## ℹ Please consider using `annotate()` or provide this layer with data containing a single row.
```

![](Analysis_youth_files/figure-html/unnamed-chunk-28-1.png)<!-- -->

``` r
p2.FM.flags<-ggplot(p2,aes(y=yvar,x=gdi.cm))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("Gender Development Index (Average 2002-2022)")+
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
## This warning is displayed once every 8 hours.
## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.
```

``` r
p2.FM.flags
```

![](Analysis_youth_files/figure-html/unnamed-chunk-28-2.png)<!-- -->

``` r
pflag_comb<-
  ggarrange(p1.FM.flags,p2.FM.flags,align = "v",
            ncol=1,nrow=2,heights=c(2,1))
```

```
## Warning in geom_text(inherit.aes = F, aes(x = 0.9, y = -0.45, label = coef_q, : All aesthetics have length 1, but the data has 282 rows.
## ℹ Please consider using `annotate()` or provide this layer with data containing a single row.
```

``` r
pflag_comb
```

![](Analysis_youth_files/figure-html/unnamed-chunk-28-3.png)<!-- -->

``` r
png(filename = 
      "../results/GDI_flags_youth.png",
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
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method [lmerModLmerTest]
## Formula: FM.z ~ gndr.c + log_gdp.z.cm + gndr.c:log_gdp.z.cm + (gndr.c |      cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##  176188.4  176260.8  -88086.2  176172.4     63039 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.6046 -0.5866  0.0251  0.5970  9.1306 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr 
##  cntry    (Intercept) 0.055971 0.23658       
##           gndr.c      0.007458 0.08636  -0.03
##  Residual             0.990587 0.99528       
## Number of obs: 63047, groups:  cntry, 33
## 
## Fixed effects:
##                     Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)          0.42882    0.04142 32.86765  10.353 7.06e-12 ***
## gndr.c               0.36792    0.01710 32.04672  21.519  < 2e-16 ***
## log_gdp.z.cm        -0.12926    0.04162 33.05324  -3.106 0.003882 ** 
## gndr.c:log_gdp.z.cm  0.07189    0.01752 33.88851   4.103 0.000242 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c lg_g..
## gndr.c      -0.028              
## lg_gdp.z.cm  0.019 -0.002       
## gndr.c:l_.. -0.002 -0.034 -0.026
```

``` r
getFE(mod2_log_GDP,round=3)
```

```
##                       Est.    SE     df      t     p     LL     UL
## (Intercept)          0.429 0.041 32.868 10.353 0.000  0.345  0.513
## gndr.c               0.368 0.017 32.047 21.519 0.000  0.333  0.403
## log_gdp.z.cm        -0.129 0.042 33.053 -3.106 0.004 -0.214 -0.045
## gndr.c:log_gdp.z.cm  0.072 0.018 33.889  4.103 0.000  0.036  0.108
```

``` r
getVC(mod2_log_GDP)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.24 0.06
## 2    cntry      gndr.c   <NA>  0.09 0.01
## 3    cntry (Intercept) gndr.c -0.03 0.00
## 4 Residual        <NA>   <NA>  1.00 0.99
```

``` r
r2mlm(mod2_log_GDP,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.045536953
## slope variation 0.001695383
## mean variation  0.050976363
## sigma2          0.901791300
## 
## $R2s
##           total
## f   0.045536953
## v   0.001695383
## m   0.050976363
## fv  0.047232337
## fvm 0.098208700
```

``` r
#anova(mod2,mod2_log_GDP)
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
## Warning in ddsc_ml(model = mod2_log_GDP, predictor = "log_gdp.z.cm", moderator = "gndr.c", :
## Predictor not properly standardized, SD = 1.01176689233303
```

``` r
t2<-Sys.time()
t2-t1
```

```
## Time difference of 4.3027 secs
```

``` r
round(ddsc_mod2_log_GDP$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.085        0.292        0.991     1.076 0.079    987.576 0.987   0.988
## 2        0.5         0.066        0.257        0.991     1.057 0.062    922.939 0.983   0.984
```

``` r
round(ddsc_mod2_log_GDP$descriptives,3)
```

```
##                          M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled log_gdp.z.cm
## means_y1             0.612 0.266    1.000           1.000    0.903           0.903       -0.346
## means_y1_scaled      2.169 0.940    1.000           1.000    0.903           0.903       -0.346
## means_y2             0.250 0.298    0.903           0.903    1.000           1.000       -0.566
## means_y2_scaled      0.884 1.056    0.903           0.903    1.000           1.000       -0.566
## log_gdp.z.cm        -0.022 1.012   -0.346          -0.346   -0.566          -0.566        1.000
## log_gdp.z.cm_scaled  0.000 1.000   -0.346          -0.346   -0.566          -0.566        1.000
## diff_score           0.363 0.128   -0.030          -0.030   -0.456          -0.456        0.601
## diff_score_scaled    1.285 0.453   -0.030          -0.030   -0.456          -0.456        0.601
##                     log_gdp.z.cm_scaled diff_score diff_score_scaled
## means_y1                         -0.346     -0.030            -0.030
## means_y1_scaled                  -0.346     -0.030            -0.030
## means_y2                         -0.566     -0.456            -0.456
## means_y2_scaled                  -0.566     -0.456            -0.456
## log_gdp.z.cm                      1.000      0.601             0.601
## log_gdp.z.cm_scaled               1.000      0.601             0.601
## diff_score                        0.601      1.000             1.000
## diff_score_scaled                 0.601      1.000             1.000
```

``` r
round(ddsc_mod2_log_GDP$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                      -0.562 0.137 33.889  -4.103   0.000   -0.840   -0.283
## w_11                         -0.165 0.043 33.115  -3.864   0.000   -0.252   -0.078
## w_21                         -0.093 0.042 33.016  -2.206   0.034   -0.179   -0.007
## r_xy1                        -0.622 0.161 33.115  -3.864   0.000   -0.950   -0.295
## r_xy2                        -0.313 0.142 33.016  -2.206   0.034   -0.601   -0.024
## b_11                         -0.586 0.152 33.115  -3.864   0.000   -0.895   -0.277
## b_21                         -0.331 0.150 33.016  -2.206   0.034   -0.636   -0.026
## main_effect                  -0.129 0.042 33.053  -3.106   0.004   -0.214   -0.045
## moderator_effect              0.368 0.017 32.047  21.519   0.000    0.333    0.403
## interaction                   0.072 0.018 33.889   4.103   0.000    0.036    0.108
## q_b11_b21                    -0.328    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                  -0.405    NA     NA      NA      NA       NA       NA
## cross_over_point             -5.118    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.057 0.045 32.816  -1.282   0.209   -0.148    0.034
## interaction_vs_main_bscale   -0.203 0.159 32.816  -1.282   0.209   -0.526    0.119
## interaction_vs_main_rscale   -0.158 0.143 32.757  -1.108   0.276   -0.449    0.132
## dadas                        -0.187 0.085 33.016  -2.206   0.983   -0.359   -0.014
## dadas_bscale                 -0.662 0.300 33.016  -2.206   0.983   -1.273   -0.051
## dadas_rscale                 -0.626 0.284 33.016  -2.206   0.983   -1.203   -0.049
## abs_diff                      0.072 0.018 33.889   4.103   0.000    0.036    0.108
## abs_sum                       0.259 0.083 33.053   3.106   0.002    0.089    0.428
## abs_diff_bscale               0.255 0.062 33.889   4.103   0.000    0.129    0.381
## abs_sum_bscale                0.917 0.295 33.053   3.106   0.002    0.316    1.518
## abs_diff_rscale               0.309 0.065 35.055   4.748   0.000    0.177    0.441
## abs_sum_rscale                0.935 0.296 33.054   3.154   0.002    0.332    1.538
```

``` r
round(ddsc_mod2_log_GDP$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.010 -0.324  2.904  1.000  0.088
```

``` r
d_log_GDP<-ddsc_mod2_log_GDP$ddsc_sem_fit$data

ddsc_sem_log_GDP<-
  ddsc_sem(data=d_log_GDP,x = "log_gdp.z.cm",y1="means_y2",
           y2="means_y1",q_sesoi = .10)
round(ddsc_sem_log_GDP$results,3)
```

```
##                                    est    se       z pvalue ci.lower ci.upper
## r_xy1_y2                        -0.601 0.139  -4.317  0.000   -0.874   -0.328
## r_xy1                           -0.566 0.144  -3.944  0.000   -0.847   -0.285
## r_xy2                           -0.346 0.163  -2.120  0.034   -0.666   -0.026
## b_11                            -0.598 0.152  -3.944  0.000   -0.895   -0.301
## b_21                            -0.326 0.154  -2.120  0.034   -0.627   -0.025
## b_10                             0.884 0.149   5.920  0.000    0.591    1.176
## b_20                             2.169 0.151  14.339  0.000    1.872    2.465
## res_cov_y1_y2                    0.681 0.176   3.877  0.000    0.337    1.026
## diff_b10_b20                    -1.285 0.062 -20.693  0.000   -1.407   -1.163
## diff_b11_b21                    -0.272 0.063  -4.317  0.000   -0.396   -0.149
## diff_rxy1_rxy2                  -0.220 0.066  -3.317  0.001   -0.350   -0.090
## q_b11_b21                       -0.352 0.105  -3.352  0.001   -0.558   -0.146
## q_rxy1_rxy2                     -0.281 0.086  -3.274  0.001   -0.448   -0.113
## cross_over_point                -4.720 1.117  -4.226  0.000   -6.909   -2.531
## sum_b11_b21                     -0.923 0.299  -3.092  0.002   -1.509   -0.338
## main_effect                     -0.462 0.149  -3.092  0.002   -0.754   -0.169
## interaction_vs_main_effect      -0.189 0.164  -1.155  0.248   -0.511    0.132
## diff_abs_b11_abs_b21             0.272 0.063   4.317  0.000    0.149    0.396
## abs_diff_b11_b21                 0.272 0.063   4.317  0.000    0.149    0.396
## abs_sum_b11_b21                  0.923 0.299   3.092  0.001    0.338    1.509
## dadas                           -0.651 0.307  -2.120  0.983   -1.253   -0.049
## q_r_equivalence                  0.181 0.086   2.107  0.982       NA       NA
## q_b_equivalence                  0.252 0.105   2.399  0.992       NA       NA
## cross_over_point_equivalence     4.720 1.117   4.226  1.000       NA       NA
## cross_over_point_minimal_effect  4.720 1.117   4.226  0.000       NA       NA
```

``` r
round(ddsc_sem_log_GDP$variance_test,3)
```

```
##             est    se      z pvalue ci.lower ci.upper
## cov_y1y2  0.870 0.226  3.851  0.000    0.427    1.313
## var_y1    1.082 0.266  4.062  0.000    0.560    1.604
## var_y2    0.858 0.211  4.062  0.000    0.444    1.271
## var_diff  0.224 0.154  1.455  0.146   -0.078    0.526
## var_ratio 1.261 0.188  6.697  0.000    0.892    1.631
## cor_y1y2  0.903 0.032 28.214  0.000    0.841    0.966
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
flag_points_log_GDP
```

```
##    (Intercept)    gndr.c        women        men  mean_level difference cntry log_gdp.cm
## 1   0.46910842 0.4485451  0.244835875 0.69338097  0.46910842  0.4485451    AT  11.008132
## 2   0.32356214 0.3808445  0.133139904 0.51398438  0.32356214  0.3808445    BE  10.946219
## 3   0.47387200 0.3677068  0.290018620 0.65772537  0.47387200  0.3677068    BG  10.024001
## 4   0.22591935 0.4771312 -0.012646257 0.46448497  0.22591935  0.4771312    CH  11.222820
## 5   0.15303188 0.2372716  0.034396071 0.27166769  0.15303188  0.2372716    CY  10.641538
## 6   0.93057045 0.4914030  0.684868947 1.17627195  0.93057045  0.4914030    CZ  10.602096
## 7   0.20782346 0.4492190 -0.016786057 0.43243298  0.20782346  0.4492190    DE  10.951411
## 8   0.37107531 0.5082305  0.116960053 0.62519056  0.37107531  0.5082305    DK  11.041181
## 9   0.38115603 0.4346623  0.163824871 0.59848719  0.38115603  0.4346623    EE  10.453078
## 10 -0.06706344 0.3145103 -0.224318604 0.09019173 -0.06706344  0.3145103    ES  10.680598
## 11  0.10927543 0.5032436 -0.142346355 0.36089721  0.10927543  0.5032436    FI  10.899727
## 12  0.04684201 0.5265046 -0.216410272 0.31009429  0.04684201  0.5265046    FR  10.820760
## 13  0.24472225 0.4807966  0.004323969 0.48512054  0.24472225  0.4807966    GB  10.795470
## 14  0.49191071 0.2018644  0.390978498 0.59284292  0.49191071  0.2018644    GR  10.458468
## 15  0.25062695 0.3869044  0.057174746 0.44407916  0.25062695  0.3869044    HR  10.312802
## 16  0.55793675 0.3009417  0.407465893 0.70840761  0.55793675  0.3009417    HU  10.328579
## 17  0.42227675 0.3629228  0.240815337 0.60373817  0.42227675  0.3629228    IE  11.199280
## 18  0.46176273 0.2728115  0.325356964 0.59816850  0.46176273  0.2728115    IL  10.568869
## 19 -0.11461611 0.4352629 -0.332247549 0.10301532 -0.11461611  0.4352629    IS  10.972691
## 20  0.43510299 0.3440312  0.263087404 0.60711857  0.43510299  0.3440312    IT  10.811351
## 21  1.03984359 0.3466629  0.866512150 1.21317503  1.03984359  0.3466629    LT  10.352204
## 22  0.89184510 0.3648590  0.709415623 1.07427458  0.89184510  0.3648590    LV  10.234123
## 23  0.47184627 0.2235117  0.360090429 0.58360210  0.47184627  0.2235117    ME   9.896653
## 24  0.55600859 0.4532554  0.329380896 0.78263628  0.55600859  0.4532554    NL  11.045016
## 25  0.34824643 0.4001341  0.148179355 0.54831350  0.34824643  0.4001341    NO  11.351260
## 26  0.39586987 0.4233198  0.184209959 0.60752977  0.39586987  0.4233198    PL  10.266271
## 27  0.40599286 0.2625674  0.274709176 0.53727655  0.40599286  0.2625674    PT  10.512324
## 28  0.85454708 0.1750750  0.767009567 0.94208460  0.85454708  0.1750750    RU  10.398951
## 29  0.25763869 0.4825205  0.016378456 0.49889893  0.25763869  0.4825205    SE  10.945268
## 30  0.56045022 0.3714961  0.374702187 0.74619825  0.56045022  0.3714961    SI  10.570158
## 31  0.74523002 0.3195112  0.585474428 0.90498562  0.74523002  0.3195112    SK  10.305981
## 32  0.65783058 0.1646775  0.575491842 0.74016932  0.65783058  0.1646775    TR  10.009841
## 33  0.68026564 0.2534444  0.553543418 0.80698786  0.68026564  0.2534444    UA   9.699690
```

``` r
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
  xlab("log(GDP per capita) (Average 2002-2022)")+
  #ylim(c(min.y.pred,max.y.pred))+
  #xlim(c(0.60,1.00))+
  ylab("Value male-typicality (Average 2002-2022)")+
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

p1.FM.flags
```

```
## Warning in geom_text(inherit.aes = F, aes(x = 9.75, y = -0.45, label = coef_q, : All aesthetics have length 1, but the data has 3302 rows.
## ℹ Please consider using `annotate()` or provide this layer with data containing a single row.
```

![](Analysis_youth_files/figure-html/unnamed-chunk-31-1.png)<!-- -->

``` r
p2.FM.flags<-ggplot(p2,aes(y=yvar,x=log_gdp.cm))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("log(GDP per capita) (Average 2002-2022)")+
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

p2.FM.flags
```

![](Analysis_youth_files/figure-html/unnamed-chunk-31-2.png)<!-- -->

``` r
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

![](Analysis_youth_files/figure-html/unnamed-chunk-31-3.png)<!-- -->

``` r
png(filename = 
      "../results/log_GDP_flags_youth.png",
    units = "cm",
    width = 21.0,height=29.7*(3/4),res = 300)
pflag_comb
dev.off()
```

```
## png 
##   2
```


## mod3: fixed effect of time (Ess round)


``` r
mod3<-lmer(FM.z~gndr.c+essround.c+(gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod3)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method [lmerModLmerTest]
## Formula: FM.z ~ gndr.c + essround.c + (gndr.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##  175790.6  175853.9  -87888.3  175776.6     63040 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.7551 -0.5862  0.0240  0.5956  9.2646 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr 
##  cntry    (Intercept) 0.07253  0.2693        
##           gndr.c      0.01239  0.1113   -0.31
##  Residual             0.98408  0.9920        
## Number of obs: 63047, groups:  cntry, 33
## 
## Fixed effects:
##               Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)  4.317e-01  4.708e-02  3.291e+01    9.17 1.39e-10 ***
## gndr.c       3.688e-01  2.105e-02  3.102e+01   17.52  < 2e-16 ***
## essround.c  -2.895e-02  1.418e-03  6.291e+04  -20.42  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##            (Intr) gndr.c
## gndr.c     -0.283       
## essround.c  0.000 -0.001
```

``` r
getFE(mod3,round=3)
```

```
##               Est.    SE        df       t p     LL     UL
## (Intercept)  0.432 0.047    32.906   9.170 0  0.336  0.528
## gndr.c       0.369 0.021    31.019  17.521 0  0.326  0.412
## essround.c  -0.029 0.001 62913.751 -20.416 0 -0.032 -0.026
```

``` r
getVC(mod3)
```

```
##        grp        var1   var2 sdcor  vcov
## 1    cntry (Intercept)   <NA>  0.27  0.07
## 2    cntry      gndr.c   <NA>  0.11  0.01
## 3    cntry (Intercept) gndr.c -0.31 -0.01
## 4 Residual        <NA>   <NA>  0.99  0.98
```

``` r
r2mlm(mod3,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.036330941
## slope variation 0.002813456
## mean variation  0.066226634
## sigma2          0.894628969
## 
## $R2s
##           total
## f   0.036330941
## v   0.002813456
## m   0.066226634
## fv  0.039144397
## fvm 0.105371031
```

``` r
anova(mod2,mod3)
```

```
## Data: diff_dat
## Models:
## mod2: FM.z ~ gndr.c + (gndr.c | cntry)
## mod3: FM.z ~ gndr.c + essround.c + (gndr.c | cntry)
##      npar    AIC    BIC logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod2    6 176204 176258 -88096    176192                         
## mod3    7 175791 175854 -87888    175777 415.41  1  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

* Alternative model with country-mean centered gender and time


``` r
mod3_alt<-lmer(FM.z~gndr.c.gmc+essround.c.gmc+(gndr.c.gmc|cntry),
               data=diff_dat_gmc,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod3_alt)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method [lmerModLmerTest]
## Formula: FM.z ~ gndr.c.gmc + essround.c.gmc + (gndr.c.gmc | cntry)
##    Data: diff_dat_gmc
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##  175789.3  175852.6  -87887.6  175775.3     63040 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.7584 -0.5863  0.0240  0.5956  9.2616 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr 
##  cntry    (Intercept) 0.07019  0.2649        
##           gndr.c.gmc  0.01232  0.1110   -0.32
##  Residual             0.98408  0.9920        
## Number of obs: 63047, groups:  cntry, 33
## 
## Fixed effects:
##                  Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)     4.264e-01  4.632e-02  3.288e+01   9.204 1.28e-10 ***
## gndr.c.gmc      3.692e-01  2.100e-02  3.110e+01  17.582  < 2e-16 ***
## essround.c.gmc -2.899e-02  1.419e-03  6.299e+04 -20.425  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr..
## gndr.c.gmc  -0.296       
## essrnd.c.gm -0.002 -0.001
```

``` r
getFE(mod3_alt,round=3)
```

```
##                  Est.    SE        df       t p     LL     UL
## (Intercept)     0.426 0.046    32.883   9.204 0  0.332  0.521
## gndr.c.gmc      0.369 0.021    31.096  17.582 0  0.326  0.412
## essround.c.gmc -0.029 0.001 62991.318 -20.425 0 -0.032 -0.026
```

``` r
getVC(mod3_alt)
```

```
##        grp        var1       var2 sdcor  vcov
## 1    cntry (Intercept)       <NA>  0.26  0.07
## 2    cntry  gndr.c.gmc       <NA>  0.11  0.01
## 3    cntry (Intercept) gndr.c.gmc -0.32 -0.01
## 4 Residual        <NA>       <NA>  0.99  0.98
```

``` r
r2mlm(mod3_alt,bargraph = F)
```

```
## $Decompositions
##                       total      within between
## fixed, within   0.035740804 0.038184945      NA
## fixed, between  0.000000000          NA       0
## slope variation 0.002794202 0.002985284      NA
## mean variation  0.064007961          NA       1
## sigma2          0.897457032 0.958829771      NA
## 
## $R2s
##           total      within between
## f1  0.035740804 0.038184945      NA
## f2  0.000000000          NA       0
## v   0.002794202 0.002985284      NA
## m   0.064007961          NA       1
## f   0.035740804          NA      NA
## fv  0.038535006 0.041170229      NA
## fvm 0.102542968          NA      NA
```

``` r
r2mlm_comp(mod2_alt,mod3_alt,bargraph = F)
```

```
## $`Model A R2s`
##           total      within between
## f1  0.030790237 0.032904818      NA
## f2  0.000000000          NA       0
## v   0.002739748 0.002927906      NA
## m   0.064263559          NA       1
## f   0.030790237          NA      NA
## fv  0.033529986 0.035832724      NA
## fvm 0.097793545          NA      NA
## 
## $`Model B R2s`
##           total      within between
## f1  0.035740804 0.038184945      NA
## f2  0.000000000          NA       0
## v   0.002794202 0.002985284      NA
## m   0.064007961          NA       1
## f   0.035740804          NA      NA
## fv  0.038535006 0.041170229      NA
## fvm 0.102542968          NA      NA
## 
## $`R2 differences, Model B - Model A`
##             total       within between
## f1   4.950567e-03 5.280127e-03      NA
## f2   0.000000e+00           NA       0
## v    5.445363e-05 5.737791e-05      NA
## m   -2.555978e-04           NA       0
## f    4.950567e-03           NA      NA
## fv   5.005021e-03 5.337505e-03      NA
## fvm  4.749423e-03           NA      NA
```

## mod4: random effect of time (Ess round)


``` r
mod4<-lmer(FM.z~gndr.c+essround.c+(gndr.c+essround.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod4)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method [lmerModLmerTest]
## Formula: FM.z ~ gndr.c + essround.c + (gndr.c + essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##  175233.5  175324.0  -87606.8  175213.5     63037 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.7953 -0.5877  0.0229  0.5914  9.2989 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr       
##  cntry    (Intercept) 0.084348 0.29043             
##           gndr.c      0.012948 0.11379  -0.36      
##           essround.c  0.002518 0.05018   0.47 -0.60
##  Residual             0.973759 0.98679             
## Number of obs: 63047, groups:  cntry, 33
## 
## Fixed effects:
##              Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)  0.433018   0.050939 32.669373   8.501 8.69e-10 ***
## gndr.c       0.369860   0.021402 31.192877  17.281  < 2e-16 ***
## essround.c  -0.021070   0.008991 26.006355  -2.344    0.027 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##            (Intr) gndr.c
## gndr.c     -0.326       
## essround.c  0.443 -0.548
```

``` r
getFE(mod4,round=3)
```

```
##               Est.    SE     df      t     p     LL     UL
## (Intercept)  0.433 0.051 32.669  8.501 0.000  0.329  0.537
## gndr.c       0.370 0.021 31.193 17.281 0.000  0.326  0.413
## essround.c  -0.021 0.009 26.006 -2.344 0.027 -0.040 -0.003
```

``` r
getVC(mod4)
```

```
##        grp        var1       var2 sdcor  vcov
## 1    cntry (Intercept)       <NA>  0.29  0.08
## 2    cntry      gndr.c       <NA>  0.11  0.01
## 3    cntry  essround.c       <NA>  0.05  0.00
## 4    cntry (Intercept)     gndr.c -0.36 -0.01
## 5    cntry (Intercept) essround.c  0.47  0.01
## 6    cntry      gndr.c essround.c -0.60  0.00
## 7 Residual        <NA>       <NA>  0.99  0.97
```

``` r
r2mlm(mod4,bargraph = F)
```

```
## $Decompositions
##                      total
## fixed           0.03348391
## slope variation 0.01937448
## mean variation  0.07327551
## sigma2          0.87386609
## 
## $R2s
##          total
## f   0.03348391
## v   0.01937448
## m   0.07327551
## fv  0.05285839
## fvm 0.12613391
```

``` r
anova(mod3,mod4)
```

```
## Data: diff_dat
## Models:
## mod3: FM.z ~ gndr.c + essround.c + (gndr.c | cntry)
## mod4: FM.z ~ gndr.c + essround.c + (gndr.c + essround.c | cntry)
##      npar    AIC    BIC logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod3    7 175791 175854 -87888    175777                         
## mod4   10 175234 175324 -87607    175214 563.05  3  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```


* Alternative model with country-mean centered gender and time


``` r
mod4_alt<-lmer(FM.z~gndr.c.gmc+essround.c.gmc+(gndr.c.gmc+essround.c.gmc|cntry),
               data=diff_dat_gmc,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod4_alt)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method [lmerModLmerTest]
## Formula: FM.z ~ gndr.c.gmc + essround.c.gmc + (gndr.c.gmc + essround.c.gmc |      cntry)
##    Data: diff_dat_gmc
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##  175230.1  175320.6  -87605.0  175210.1     63037 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.8234 -0.5877  0.0231  0.5914  9.3133 
## 
## Random effects:
##  Groups   Name           Variance Std.Dev. Corr       
##  cntry    (Intercept)    0.070628 0.26576             
##           gndr.c.gmc     0.013260 0.11515  -0.28      
##           essround.c.gmc 0.002843 0.05332   0.32 -0.62
##  Residual                0.973690 0.98676             
## Number of obs: 63047, groups:  cntry, 33
## 
## Fixed effects:
##                 Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)     0.428414   0.046464 32.868704   9.220 1.23e-10 ***
## gndr.c.gmc      0.369194   0.021629 30.696538  17.070  < 2e-16 ***
## essround.c.gmc -0.019393   0.009576 23.115103  -2.025   0.0545 .  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr..
## gndr.c.gmc  -0.260       
## essrnd.c.gm  0.305 -0.568
```

``` r
getFE(mod4_alt,round=3)
```

```
##                  Est.    SE     df      t     p     LL    UL
## (Intercept)     0.428 0.046 32.869  9.220 0.000  0.334 0.523
## gndr.c.gmc      0.369 0.022 30.697 17.070 0.000  0.325 0.413
## essround.c.gmc -0.019 0.010 23.115 -2.025 0.055 -0.039 0.000
```

``` r
getVC(mod4_alt)
```

```
##        grp           var1           var2 sdcor  vcov
## 1    cntry    (Intercept)           <NA>  0.27  0.07
## 2    cntry     gndr.c.gmc           <NA>  0.12  0.01
## 3    cntry essround.c.gmc           <NA>  0.05  0.00
## 4    cntry    (Intercept)     gndr.c.gmc -0.28 -0.01
## 5    cntry    (Intercept) essround.c.gmc  0.32  0.00
## 6    cntry     gndr.c.gmc essround.c.gmc -0.62  0.00
## 7 Residual           <NA>           <NA>  0.99  0.97
```

``` r
r2mlm(mod4_alt,bargraph = F)
```

```
## $Decompositions
##                      total     within between
## fixed, within   0.03286253 0.03511239      NA
## fixed, between  0.00000000         NA       0
## slope variation 0.01969402 0.02104233      NA
## mean variation  0.06407611         NA       1
## sigma2          0.88336735 0.94384528      NA
## 
## $R2s
##          total     within between
## f1  0.03286253 0.03511239      NA
## f2  0.00000000         NA       0
## v   0.01969402 0.02104233      NA
## m   0.06407611         NA       1
## f   0.03286253         NA      NA
## fv  0.05255655 0.05615472      NA
## fvm 0.11663265         NA      NA
```

``` r
r2mlm_comp(mod3_alt,mod4_alt,bargraph = F)
```

```
## $`Model A R2s`
##           total      within between
## f1  0.035740804 0.038184945      NA
## f2  0.000000000          NA       0
## v   0.002794202 0.002985284      NA
## m   0.064007961          NA       1
## f   0.035740804          NA      NA
## fv  0.038535006 0.041170229      NA
## fvm 0.102542968          NA      NA
## 
## $`Model B R2s`
##          total     within between
## f1  0.03286253 0.03511239      NA
## f2  0.00000000         NA       0
## v   0.01969402 0.02104233      NA
## m   0.06407611         NA       1
## f   0.03286253         NA      NA
## fv  0.05255655 0.05615472      NA
## fvm 0.11663265         NA      NA
## 
## $`R2 differences, Model B - Model A`
##             total       within between
## f1  -0.0028782783 -0.003072554      NA
## f2   0.0000000000           NA       0
## v    0.0168998184  0.018057047      NA
## m    0.0000681442           NA       0
## f   -0.0028782783           NA      NA
## fv   0.0140215400  0.014984494      NA
## fvm  0.0140896842           NA      NA
```

``` r
anova(mod3_alt,mod4_alt)
```

```
## Data: diff_dat_gmc
## Models:
## mod3_alt: FM.z ~ gndr.c.gmc + essround.c.gmc + (gndr.c.gmc | cntry)
## mod4_alt: FM.z ~ gndr.c.gmc + essround.c.gmc + (gndr.c.gmc + essround.c.gmc | cntry)
##          npar    AIC    BIC logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod3_alt    7 175789 175853 -87888    175775                         
## mod4_alt   10 175230 175321 -87605    175210 565.17  3  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## mod5: fixed interaction between time and gender


``` r
mod5<-lmer(FM.z~gndr.c+essround.c+
             gndr.c:essround.c+(gndr.c+essround.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod5)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method [lmerModLmerTest]
## Formula: FM.z ~ gndr.c + essround.c + gndr.c:essround.c + (gndr.c + essround.c |      cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##  175218.7  175318.3  -87598.3  175196.7     63036 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.8183 -0.5880  0.0223  0.5917  9.3258 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr       
##  cntry    (Intercept) 0.084390 0.29050             
##           gndr.c      0.013460 0.11602  -0.37      
##           essround.c  0.002521 0.05021   0.47 -0.60
##  Residual             0.973481 0.98665             
## Number of obs: 63047, groups:  cntry, 33
## 
## Fixed effects:
##                     Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)        4.337e-01  5.095e-02  3.268e+01   8.513 8.40e-10 ***
## gndr.c             3.689e-01  2.177e-02  3.147e+01  16.948  < 2e-16 ***
## essround.c        -2.114e-02  8.996e-03  2.606e+01  -2.350   0.0266 *  
## gndr.c:essround.c -1.143e-02  2.783e-03  3.048e+04  -4.105 4.05e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c essrn.
## gndr.c      -0.335              
## essround.c   0.449 -0.542       
## gndr.c:ssr. -0.003  0.009  0.001
```

``` r
getFE(mod5,round=3)
```

```
##                     Est.    SE        df      t     p     LL     UL
## (Intercept)        0.434 0.051    32.679  8.513 0.000  0.330  0.537
## gndr.c             0.369 0.022    31.467 16.948 0.000  0.325  0.413
## essround.c        -0.021 0.009    26.064 -2.350 0.027 -0.040 -0.003
## gndr.c:essround.c -0.011 0.003 30481.891 -4.105 0.000 -0.017 -0.006
```

``` r
getVC(mod5)
```

```
##        grp        var1       var2 sdcor  vcov
## 1    cntry (Intercept)       <NA>  0.29  0.08
## 2    cntry      gndr.c       <NA>  0.12  0.01
## 3    cntry  essround.c       <NA>  0.05  0.00
## 4    cntry (Intercept)     gndr.c -0.37 -0.01
## 5    cntry (Intercept) essround.c  0.47  0.01
## 6    cntry      gndr.c essround.c -0.60  0.00
## 7 Residual        <NA>       <NA>  0.99  0.97
```

``` r
r2mlm(mod5,bargraph = F)
```

```
## $Decompositions
##                      total
## fixed           0.03394044
## slope variation 0.01950248
## mean variation  0.07326518
## sigma2          0.87329190
## 
## $R2s
##          total
## f   0.03394044
## v   0.01950248
## m   0.07326518
## fv  0.05344292
## fvm 0.12670810
```

``` r
anova(mod4,mod5)
```

```
## Data: diff_dat
## Models:
## mod4: FM.z ~ gndr.c + essround.c + (gndr.c + essround.c | cntry)
## mod5: FM.z ~ gndr.c + essround.c + gndr.c:essround.c + (gndr.c + essround.c | cntry)
##      npar    AIC    BIC logLik -2*log(L) Chisq Df Pr(>Chisq)    
## mod4   10 175234 175324 -87607    175214                        
## mod5   11 175219 175318 -87598    175197 16.82  1  4.109e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

* Alternative model with country-mean centered gender and time


``` r
mod5_alt<-lmer(FM.z~gndr.c.gmc+essround.c.gmc+gndr.c.gmc:essround.c.gmc+
                 (gndr.c.gmc+essround.c.gmc|cntry),
               data=diff_dat_gmc,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod5_alt)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method [lmerModLmerTest]
## Formula: FM.z ~ gndr.c.gmc + essround.c.gmc + gndr.c.gmc:essround.c.gmc +  
##     (gndr.c.gmc + essround.c.gmc | cntry)
##    Data: diff_dat_gmc
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##  175214.1  175313.7  -87596.1  175192.1     63036 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.8512 -0.5881  0.0226  0.5916  9.3352 
## 
## Random effects:
##  Groups   Name           Variance Std.Dev. Corr       
##  cntry    (Intercept)    0.070624 0.26575             
##           gndr.c.gmc     0.013247 0.11510  -0.28      
##           essround.c.gmc 0.002855 0.05343   0.32 -0.63
##  Residual                0.973416 0.98662             
## Number of obs: 63047, groups:  cntry, 33
## 
## Fixed effects:
##                             Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)                4.285e-01  4.646e-02  3.287e+01   9.222 1.23e-10 ***
## gndr.c.gmc                 3.698e-01  2.162e-02  3.063e+01  17.108  < 2e-16 ***
## essround.c.gmc            -1.908e-02  9.594e-03  2.314e+01  -1.989   0.0587 .  
## gndr.c.gmc:essround.c.gmc -1.197e-02  2.823e-03  6.297e+04  -4.242 2.21e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.. essr..
## gndr.c.gmc  -0.261              
## essrnd.c.gm  0.307 -0.574       
## gndr.c.g:..  0.000 -0.007 -0.007
```

``` r
getFE(mod5_alt,round=3)
```

```
##                             Est.    SE        df      t     p     LL     UL
## (Intercept)                0.428 0.046    32.869  9.222 0.000  0.334  0.523
## gndr.c.gmc                 0.370 0.022    30.633 17.108 0.000  0.326  0.414
## essround.c.gmc            -0.019 0.010    23.139 -1.989 0.059 -0.039  0.001
## gndr.c.gmc:essround.c.gmc -0.012 0.003 62973.431 -4.242 0.000 -0.018 -0.006
```

``` r
getVC(mod5_alt)
```

```
##        grp           var1           var2 sdcor  vcov
## 1    cntry    (Intercept)           <NA>  0.27  0.07
## 2    cntry     gndr.c.gmc           <NA>  0.12  0.01
## 3    cntry essround.c.gmc           <NA>  0.05  0.00
## 4    cntry    (Intercept)     gndr.c.gmc -0.28 -0.01
## 5    cntry    (Intercept) essround.c.gmc  0.32  0.00
## 6    cntry     gndr.c.gmc essround.c.gmc -0.63  0.00
## 7 Residual           <NA>           <NA>  0.99  0.97
```

``` r
r2mlm(mod5_alt,bargraph = F)
```

```
## $Decompositions
##                      total
## fixed           0.03311114
## slope variation 0.01975860
## mean variation  0.06406875
## sigma2          0.88306151
## 
## $R2s
##          total
## f   0.03311114
## v   0.01975860
## m   0.06406875
## fv  0.05286974
## fvm 0.11693849
```

``` r
#r2mlm_comp(mod4_alt,mod5_alt,bargraph = F)
anova(mod4_alt,mod5_alt)
```

```
## Data: diff_dat_gmc
## Models:
## mod4_alt: FM.z ~ gndr.c.gmc + essround.c.gmc + (gndr.c.gmc + essround.c.gmc | cntry)
## mod5_alt: FM.z ~ gndr.c.gmc + essround.c.gmc + gndr.c.gmc:essround.c.gmc + (gndr.c.gmc + essround.c.gmc | cntry)
##          npar    AIC    BIC logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod4_alt   10 175230 175321 -87605    175210                         
## mod5_alt   11 175214 175314 -87596    175192 17.992  1  2.218e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```


## mod6: random interaction between time and gender


``` r
mod6<-lmer(FM.z~gndr.c+essround.c+
             gndr.c:essround.c+(gndr.c+essround.c+gndr.c:essround.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method [lmerModLmerTest]
## Formula: FM.z ~ gndr.c + essround.c + gndr.c:essround.c + (gndr.c + essround.c +  
##     gndr.c:essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##  175218.0  175353.8  -87594.0  175188.0     63032 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.8325 -0.5885  0.0215  0.5923  9.3452 
## 
## Random effects:
##  Groups   Name              Variance  Std.Dev. Corr             
##  cntry    (Intercept)       0.0845008 0.29069                   
##           gndr.c            0.0139007 0.11790  -0.37            
##           essround.c        0.0025323 0.05032   0.47 -0.61      
##           gndr.c:essround.c 0.0001655 0.01286  -0.45  0.06 -0.34
##  Residual                   0.9731515 0.98648                   
## Number of obs: 63047, groups:  cntry, 33
## 
## Fixed effects:
##                    Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)        0.433942   0.050985 32.673613   8.511 8.44e-10 ***
## gndr.c             0.369217   0.022171 30.097544  16.653  < 2e-16 ***
## essround.c        -0.021149   0.009015 25.999947  -2.346  0.02688 *  
## gndr.c:essround.c -0.011538   0.003696 27.132972  -3.122  0.00424 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c essrn.
## gndr.c      -0.336              
## essround.c   0.451 -0.556       
## gndr.c:ssr. -0.274  0.030 -0.205
```

``` r
getFE(mod6,round=3)
```

```
##                     Est.    SE     df      t     p     LL     UL
## (Intercept)        0.434 0.051 32.674  8.511 0.000  0.330  0.538
## gndr.c             0.369 0.022 30.098 16.653 0.000  0.324  0.414
## essround.c        -0.021 0.009 26.000 -2.346 0.027 -0.040 -0.003
## gndr.c:essround.c -0.012 0.004 27.133 -3.122 0.004 -0.019 -0.004
```

``` r
getVC(mod6)
```

```
##         grp              var1              var2 sdcor  vcov
## 1     cntry       (Intercept)              <NA>  0.29  0.08
## 2     cntry            gndr.c              <NA>  0.12  0.01
## 3     cntry        essround.c              <NA>  0.05  0.00
## 4     cntry gndr.c:essround.c              <NA>  0.01  0.00
## 5     cntry       (Intercept)            gndr.c -0.37 -0.01
## 6     cntry       (Intercept)        essround.c  0.47  0.01
## 7     cntry       (Intercept) gndr.c:essround.c -0.45  0.00
## 8     cntry            gndr.c        essround.c -0.61  0.00
## 9     cntry            gndr.c gndr.c:essround.c  0.06  0.00
## 10    cntry        essround.c gndr.c:essround.c -0.34  0.00
## 11 Residual              <NA>              <NA>  0.99  0.97
```

``` r
r2mlm(mod6,bargraph = F)
```

```
## $Decompositions
##                      total
## fixed           0.03399164
## slope variation 0.01997811
## mean variation  0.07330109
## sigma2          0.87272915
## 
## $R2s
##          total
## f   0.03399164
## v   0.01997811
## m   0.07330109
## fv  0.05396975
## fvm 0.12727085
```

``` r
anova(mod5,mod6)
```

```
## Data: diff_dat
## Models:
## mod5: FM.z ~ gndr.c + essround.c + gndr.c:essround.c + (gndr.c + essround.c | cntry)
## mod6: FM.z ~ gndr.c + essround.c + gndr.c:essround.c + (gndr.c + essround.c + gndr.c:essround.c | cntry)
##      npar    AIC    BIC logLik -2*log(L)  Chisq Df Pr(>Chisq)  
## mod5   11 175219 175318 -87598    175197                       
## mod6   15 175218 175354 -87594    175188 8.6699  4     0.0699 .
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

### Trends


``` r
# gender specific change over time

change_mod6<-emmeans(mod6,specs="essround.c",by="gndr.c",
                     at=list(gndr.c=c(-0.5,0.5),
                             essround.c=rev(range(diff_dat$essround.c))),
                     disable.pbkrtest=T,
                     lmerTest.limit = 400000,infer=c(T,T))
change_mod6
```

```
## gndr.c = -0.5:
##  essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.180 0.0869 29.6  0.00249    0.358   2.072  0.0471
##        -4.5  0.319 0.0485 24.8  0.21860    0.418   6.568  <.0001
## 
## gndr.c =  0.5:
##  essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.497 0.0711 29.5  0.35210    0.643   6.995  <.0001
##        -4.5  0.740 0.0527 26.6  0.63137    0.848  14.023  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_mod6,adjust="none",infer=c(T,T))
```

```
## gndr.c = -0.5:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   -0.138 0.0861 24.9   -0.316   0.0389  -1.608  0.1205
## 
## gndr.c =  0.5:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   -0.242 0.0794 25.6   -0.406  -0.0789  -3.050  0.0053
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# change in gender differences over time

change_in_diff_mod6<-emmeans(mod6,specs=c("gndr.c","essround.c"),
                             at=list(gndr.c=c(-0.5,0.5),
                                     essround.c=rev(range(diff_dat$essround.c))),
                             disable.pbkrtest=T,
                             lmerTest.limit = 400000,infer=c(T,T))
change_in_diff_mod6
```

```
##  gndr.c essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5        4.5  0.180 0.0869 29.6  0.00249    0.358   2.072  0.0471
##     0.5        4.5  0.497 0.0711 29.5  0.35210    0.643   6.995  <.0001
##    -0.5       -4.5  0.319 0.0485 24.8  0.21860    0.418   6.568  <.0001
##     0.5       -4.5  0.740 0.0527 26.6  0.63137    0.848  14.023  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod6,adjust="none",infer=c(T,T))
```

```
##  contrast                                                 estimate     SE   df lower.CL upper.CL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5       -0.317 0.0281 22.7  -0.3755  -0.2591
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   -0.138 0.0861 24.9  -0.3158   0.0389
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)    -0.560 0.0953 27.1  -0.7550  -0.3641
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)       0.179 0.0712 24.8   0.0321   0.3256
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)       -0.242 0.0794 25.6  -0.4056  -0.0789
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)   -0.421 0.0273 31.4  -0.4768  -0.3655
##  t.ratio p.value
##  -11.287  <.0001
##   -1.608  0.1205
##   -5.874  <.0001
##    2.511  0.0189
##   -3.050  0.0053
##  -15.417  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# custom contrast to obtain the estimates for the entire span

changes_in_diff <- list(
  diff_ESS10 = c(-1,1,0,0),
  diff_ESS1 = c(0,0,-1,1))

diff_mod6<-contrast(change_in_diff_mod6,method = changes_in_diff,adjust="none",infer=c(T,T))

# gender differences at ESS1 and ESS10
diff_mod6
```

```
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10    0.317 0.0281 22.7    0.259    0.375  11.287  <.0001
##  diff_ESS1     0.421 0.0273 31.4    0.365    0.477  15.417  <.0001
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
##  diff_ESS10 - diff_ESS1   -0.104 0.0333 27.1   -0.172  -0.0356  -3.122  0.0042
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```



* Alternative model with country-mean centered gender and time


``` r
mod6_alt<-lmer(FM.z~gndr.c.gmc+essround.c.gmc+gndr.c.gmc:essround.c.gmc+
                 (gndr.c.gmc+essround.c.gmc+gndr.c.gmc:essround.c.gmc|cntry),
               data=diff_dat_gmc,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_alt)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method [lmerModLmerTest]
## Formula: FM.z ~ gndr.c.gmc + essround.c.gmc + gndr.c.gmc:essround.c.gmc +  
##     (gndr.c.gmc + essround.c.gmc + gndr.c.gmc:essround.c.gmc |          cntry)
##    Data: diff_dat_gmc
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##  175212.1  175347.9  -87591.0  175182.1     63032 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.8711 -0.5881  0.0218  0.5918  9.3552 
## 
## Random effects:
##  Groups   Name                      Variance  Std.Dev. Corr             
##  cntry    (Intercept)               0.0706270 0.26576                   
##           gndr.c.gmc                0.0131061 0.11448  -0.28            
##           essround.c.gmc            0.0029024 0.05387   0.32 -0.64      
##           gndr.c.gmc:essround.c.gmc 0.0001774 0.01332  -0.50  0.13 -0.39
##  Residual                           0.9730847 0.98645                   
## Number of obs: 63047, groups:  cntry, 33
## 
## Fixed effects:
##                            Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)                0.428501   0.046464 32.868635   9.222 1.23e-10 ***
## gndr.c.gmc                 0.369300   0.021520 30.587321  17.161  < 2e-16 ***
## essround.c.gmc            -0.018890   0.009672 22.888249  -1.953  0.06313 .  
## gndr.c.gmc:essround.c.gmc -0.012492   0.003793 25.621670  -3.294  0.00289 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.. essr..
## gndr.c.gmc  -0.264              
## essrnd.c.gm  0.308 -0.579       
## gndr.c.g:.. -0.306  0.065 -0.244
```

``` r
getFE(mod6_alt,round=3)
```

```
##                             Est.    SE     df      t     p     LL     UL
## (Intercept)                0.429 0.046 32.869  9.222 0.000  0.334  0.523
## gndr.c.gmc                 0.369 0.022 30.587 17.161 0.000  0.325  0.413
## essround.c.gmc            -0.019 0.010 22.888 -1.953 0.063 -0.039  0.001
## gndr.c.gmc:essround.c.gmc -0.012 0.004 25.622 -3.294 0.003 -0.020 -0.005
```

``` r
getVC(mod6_alt)
```

```
##         grp                      var1                      var2 sdcor  vcov
## 1     cntry               (Intercept)                      <NA>  0.27  0.07
## 2     cntry                gndr.c.gmc                      <NA>  0.11  0.01
## 3     cntry            essround.c.gmc                      <NA>  0.05  0.00
## 4     cntry gndr.c.gmc:essround.c.gmc                      <NA>  0.01  0.00
## 5     cntry               (Intercept)                gndr.c.gmc -0.28 -0.01
## 6     cntry               (Intercept)            essround.c.gmc  0.32  0.00
## 7     cntry               (Intercept) gndr.c.gmc:essround.c.gmc -0.50  0.00
## 8     cntry                gndr.c.gmc            essround.c.gmc -0.64  0.00
## 9     cntry                gndr.c.gmc gndr.c.gmc:essround.c.gmc  0.13  0.00
## 10    cntry            essround.c.gmc gndr.c.gmc:essround.c.gmc -0.39  0.00
## 11 Residual                      <NA>                      <NA>  0.99  0.97
```

``` r
r2mlm(mod6_alt,bargraph = F)
```

```
## $Decompositions
##                      total
## fixed           0.03299417
## slope variation 0.02026957
## mean variation  0.06404235
## sigma2          0.88269391
## 
## $R2s
##          total
## f   0.03299417
## v   0.02026957
## m   0.06404235
## fv  0.05326374
## fvm 0.11730609
```

``` r
#r2mlm_comp(mod4_alt,mod6_alt,bargraph = F)
anova(mod5_alt,mod6_alt)
```

```
## Data: diff_dat_gmc
## Models:
## mod5_alt: FM.z ~ gndr.c.gmc + essround.c.gmc + gndr.c.gmc:essround.c.gmc + (gndr.c.gmc + essround.c.gmc | cntry)
## mod6_alt: FM.z ~ gndr.c.gmc + essround.c.gmc + gndr.c.gmc:essround.c.gmc + (gndr.c.gmc + essround.c.gmc + gndr.c.gmc:essround.c.gmc | cntry)
##          npar    AIC    BIC logLik -2*log(L) Chisq Df Pr(>Chisq)  
## mod5_alt   11 175214 175314 -87596    175192                      
## mod6_alt   15 175212 175348 -87591    175182 10.01  4    0.04027 *
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

### Figure for time trends


``` r
# Figure for average patterns

# Model-based development for men and women

p_mod6<-
  emmip(
    mod6, 
    gndr.c ~ essround.c,
    at=list(gndr.c = c(-0.5,0.5),
            essround.c=
              seq(from=-4.5,to=4.5,by=1)),
    plotit=F,CIs=T,lmerTest.limit = 1e6,disable.pbkrtest=T)


p_mod6$gndr.c<-p_mod6$tvar
levels(p_mod6$gndr.c)<-c("Women","Men")

p_mod6<-data.frame(p_mod6)

p_mod6$year<-
  case_when(
    p_mod6$essround.c==-4.5~2002,  
    p_mod6$essround.c==-3.5~2004,
    p_mod6$essround.c==-2.5~2006,
    p_mod6$essround.c==-1.5~2008,
    p_mod6$essround.c==-0.5~2010,
    p_mod6$essround.c==0.5~2012,
    p_mod6$essround.c==1.5~2014,
    p_mod6$essround.c==2.5~2016,
    p_mod6$essround.c==3.5~2018,
    p_mod6$essround.c==4.5~2020
  )

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
  essround_i<-p_mod6[i,"essround.c"]
  gndr_i<-as.numeric(as.character(p_mod6[i,"tvar"]))
  
  temp_diff_dat<-diff_dat %>%
    filter(#cntry==cntry_i,
      essround.c==essround_i,
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
  scale_color_manual(values = my_colors) +
  scale_y_continuous(limits = c(-1.1, 1.1)) +
  scale_x_continuous(limits = c(2001, 2021), breaks = seq(2002, 2020, 2))+
  ylab("Mean-level of value male-typicality")+
  theme(legend.title=element_blank(),
        legend.text  = element_text(size = 16),
        axis.title.x = element_text(size = 18, face = "bold"),   
        axis.title.y = element_text(size = 18, face = "bold"),   
        axis.text.x  = element_text(size = 14),                  
        axis.text.y  = element_text(size = 14))

png(filename = 
      "../results/time_trends_youth.png",
    units = "cm",
    width = 21.0,height=29.7*(3/4),res = 300)
p_time_trends
dev.off()
```

```
## png 
##   2
```
#### Add year specific predictors to the figure


``` r
# obtain predicted values from which country-variance is removed

GEI_mod<-lmer(gei~(1|ISO2)+(1|year),data=long_GII_in_ESS_d)
GEI_means<-coefficients(GEI_mod)$year
GEI_means$year<-as.numeric(rownames(GEI_means))
GEI_means$gei_year_mean<-GEI_means$`(Intercept)`

p_mod6<-left_join(
  x=p_mod6,
  y=GEI_means[,c("year","gei_year_mean")],
  by="year")

GGGI_mod<-lmer(gggi~(1|ISO2)+(1|year),data=long_GGGI_in_ESS_d)
GGGI_means<-coefficients(GGGI_mod)$year
GGGI_means$year<-as.numeric(rownames(GGGI_means))
GGGI_means$gggi_year_mean<-GGGI_means$`(Intercept)`

p_mod6<-left_join(
  x=p_mod6,
  y=GGGI_means[,c("year","gggi_year_mean")],
  by="year")

GDI_mod<-lmer(gdi~(1|ISO2)+(1|year),data=long_GDI_in_ESS_d)
GDI_means<-coefficients(GDI_mod)$year
GDI_means$year<-as.numeric(rownames(GDI_means))
GDI_means$gdi_year_mean<-GDI_means$`(Intercept)`

p_mod6<-left_join(
  x=p_mod6,
  y=GDI_means[,c("year","gdi_year_mean")],
  by="year")

GDP_mod<-lmer(log_gdp~(1|ISO2)+(1|year),data=long_GDP_in_ESS_d)
GDP_means<-coefficients(GDP_mod)$year
GDP_means$year<-as.numeric(rownames(GDP_means))
GDP_means$log_gdp_year_mean<-GDP_means$`(Intercept)`

p_mod6<-left_join(
  x=p_mod6,
  y=GDP_means[,c("year","log_gdp_year_mean")],
  by="year")


p_time_trends_w_GE<-p_time_trends+
  geom_line(data = p_mod6, aes(x = year, y = gei_year_mean),
          color = "black", linewidth = 2) +
  geom_point(data = p_mod6, aes(x = year, y = gei_year_mean),
             color = "black", size = 6) +
  geom_line(data = p_mod6, aes(x = year, y = gggi_year_mean),
            color = "darkgray", linewidth = 2) +
  geom_point(data = p_mod6, aes(x = year, y = gggi_year_mean),
             color = "darkgray", size = 6)+
  geom_text(
    data = p_mod6 %>% filter(year == min(year)),
    aes(x = year + 0.0, y = gei_year_mean+0.1, label = "GEI year average"),
    color = "black", hjust = 0, vjust = -0.5, size = 5
  ) +
  geom_text(
    data = p_mod6 %>% filter(year == 2016),
    aes(x = year - 1.0, y = gggi_year_mean-0.05, label = "GGGI year average"),
    color = "darkgray", hjust = 0, vjust = 1.5, size = 5
  ) +
  ylim(-0.6,1.1)
```

```
## Scale for y is already present.
## Adding another scale for y, which will replace the existing scale.
```

``` r
png(filename = 
      "../results/time_trends_with_GE_youth.png",
    units = "cm",
    width = 21.0,height=29.7*(3/4),res = 300)
p_time_trends_w_GE
```

```
## Warning: Removed 4 rows containing missing values or values outside the scale range
## (`geom_line()`).
```

```
## Warning: Removed 6 rows containing missing values or values outside the scale range
## (`geom_point()`).
```

``` r
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
    pull(essround.c) %>%
    unique()
  
  unique_cntry_coefs<-mod6_cntry_coefs[unique_cntry, ]
  
  unique_cntry_pred<-
    data.frame(cntry=unique_cntry,
               essround.c=unique_cntry_rounds,
               gndr.c=rep(x = c(0.5,-0.5),each=length(unique_cntry_rounds)))
  pred_list[[unique_cntry]]<-unique_cntry_pred
}

pred_cntry_dat<-do.call(rbind.data.frame,pred_list)

# model based predictions for each time x country point
pred_cntry_dat$FM.z_mean<-predict(mod6,newdata=pred_cntry_dat)

pred_cntry_dat$essround=pred_cntry_dat$essround.c+5.5
pred_cntry_dat$year<-
  case_when(
    pred_cntry_dat$essround==1~2002,  
    pred_cntry_dat$essround==2~2004,
    pred_cntry_dat$essround==3~2006,
    pred_cntry_dat$essround==4~2008,
    pred_cntry_dat$essround==5~2010,
    pred_cntry_dat$essround==6~2012,
    pred_cntry_dat$essround==7~2014,
    pred_cntry_dat$essround==8~2016,
    pred_cntry_dat$essround==9~2018,
    pred_cntry_dat$essround==10~2020
  )
pred_cntry_dat$gender<-
  case_when(
    pred_cntry_dat$gndr.c==0.5~"men",
    pred_cntry_dat$gndr.c==-0.5~"women",
  )

range(pred_cntry_dat$FM.z_mean)
```

```
## [1] -0.4766659  1.2976351
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
  essround_i<-pred_cntry_dat[i,"essround"]
  gndr_i<-pred_cntry_dat[i,"gndr.c"]
  
  temp_diff_dat<-diff_dat %>%
    filter(cntry==cntry_i,
           essround==essround_i,
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

# save the plots to same pdf-file

#my_colors <- met.brewer("Cassatt2")[c(8, 3)]
my_colors <- met.brewer("Archambault")[c(6,2)]

pdf("../results/country_specific_time_trends.pdf", width = 8, height = 6)

for (ctry in countries) {
  print(
    ggplot(pred_cntry_dat[pred_cntry_dat$cntry == ctry, ], 
           aes(x = year, y = FM.z_mean, color = gender)) +
      geom_smooth(method = "lm", se = FALSE) +
      geom_point(size=8) +
      geom_point(aes(x=year,y=obs_mean_wt),size=8,shape = 1,alpha=.50)+
      geom_errorbar(aes(ymin=obs_mean_wt_LL, ymax=obs_mean_wt_UL),alpha = .50)+
      scale_color_manual(values = my_colors) +
      geom_flag(aes(country=tolower(ctry)))+
      scale_y_continuous(limits = c(-1.1, 1.6)) +
      scale_x_continuous(limits = c(2001, 2021), breaks = seq(2002, 2020, 2)) +
      ggtitle(paste("Country:", ctry))+
    ylab("Mean-level of value male-typicality")+
      theme(legend.title=element_blank())
  )
}
```

```
## `geom_smooth()` using formula = 'y ~ x'
## `geom_smooth()` using formula = 'y ~ x'
## `geom_smooth()` using formula = 'y ~ x'
## `geom_smooth()` using formula = 'y ~ x'
## `geom_smooth()` using formula = 'y ~ x'
## `geom_smooth()` using formula = 'y ~ x'
## `geom_smooth()` using formula = 'y ~ x'
## `geom_smooth()` using formula = 'y ~ x'
## `geom_smooth()` using formula = 'y ~ x'
## `geom_smooth()` using formula = 'y ~ x'
## `geom_smooth()` using formula = 'y ~ x'
## `geom_smooth()` using formula = 'y ~ x'
## `geom_smooth()` using formula = 'y ~ x'
## `geom_smooth()` using formula = 'y ~ x'
## `geom_smooth()` using formula = 'y ~ x'
## `geom_smooth()` using formula = 'y ~ x'
## `geom_smooth()` using formula = 'y ~ x'
## `geom_smooth()` using formula = 'y ~ x'
## `geom_smooth()` using formula = 'y ~ x'
## `geom_smooth()` using formula = 'y ~ x'
## `geom_smooth()` using formula = 'y ~ x'
## `geom_smooth()` using formula = 'y ~ x'
## `geom_smooth()` using formula = 'y ~ x'
## `geom_smooth()` using formula = 'y ~ x'
## `geom_smooth()` using formula = 'y ~ x'
## `geom_smooth()` using formula = 'y ~ x'
## `geom_smooth()` using formula = 'y ~ x'
## `geom_smooth()` using formula = 'y ~ x'
## `geom_smooth()` using formula = 'y ~ x'
## `geom_smooth()` using formula = 'y ~ x'
## `geom_smooth()` using formula = 'y ~ x'
## `geom_smooth()` using formula = 'y ~ x'
## `geom_smooth()` using formula = 'y ~ x'
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
ISO<-read.csv2("../data/ISO.csv")

pred_cntry_dat<-left_join(
  x=pred_cntry_dat,
  y=ISO[,c("ISO2","Country_eng_short","CLDR")],
  by=c("cntry"="ISO2")
)

#mod6_cntry_coefs

#round((mod6_cntry_coefs$`gndr.c:essround.c`/2)*18,2)

#my_colors <- met.brewer("Cassatt2")[c(8, 3)]
my_colors <- met.brewer("Archambault")[c(2,6)]

facet_plot<-
  ggplot(pred_cntry_dat, 
         aes(x = year, y = FM.z_mean, color = gender)) +
  geom_smooth(method = "lm", se = FALSE, formula="y~x") +
  geom_point(size=3.5) +
  geom_point(aes(x=year,y=obs_mean_wt),size=3.5,shape = 1,alpha=.70)+
  geom_errorbar(aes(ymin=obs_mean_wt_LL, ymax=obs_mean_wt_UL),alpha = .50)+
  scale_color_manual(values = my_colors) +
  #geom_flag(aes(country=tolower(ctry)))+
  scale_y_continuous(limits = c(-1.1, 1.6)) +
  scale_x_continuous(limits = c(2001, 2021), breaks = seq(2002, 2020, 2)) +
  #ggtitle(paste("Country:", ctry))+
  ylab("Mean-level of value male-typicality")+
  xlab("Year")+
  theme(legend.title=element_blank(),legend.position = "top",
        axis.text.x = element_text(angle = 45,size = 6,hjust=1))+
  facet_wrap(~CLDR,nrow=6,ncol=6)+
  #facet_wrap(~cntry,nrow=6,ncol=6)+
  geom_flag(aes(country=tolower(cntry)),size=2)

#test_plot


png(filename = 
      "../results/country_time_trend_facets_youth.png",
    units = "cm",
    width = 21.0,height=29.7*(3/4),res = 600)
facet_plot
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
## [1] 15.30423
```

Country-specific coefficients for time effect


``` r
cntry_specific_changes<-
  coefficients(mod6)$cntry %>%
  mutate(change_per_year=essround.c/2,
         gndr_change_per_year=`gndr.c:essround.c`/2) %>%
  mutate(men_change_per_year=change_per_year+0.5*gndr_change_per_year,
         women_change_per_year=change_per_year+(-0.5)*gndr_change_per_year) %>%
  mutate(change_per_18_years=18*change_per_year,
         gndr_change_per_18_year=18*gndr_change_per_year,
         men_change_per_18_years=18*men_change_per_year,
         women_change_per_18_years=18*women_change_per_year) %>%
  select(gndr.c,change_per_18_years,gndr_change_per_18_year,men_change_per_18_years,women_change_per_18_years) %>%
  round(.,2)

cntry_specific_changes$cntry<-rownames(cntry_specific_changes)

cntry_specific_changes<-
  left_join(x=cntry_specific_changes,
            y=n_rounds,
            by="cntry")

cntry_specific_changes
```

```
##    gndr.c change_per_18_years gndr_change_per_18_year men_change_per_18_years
## 1    0.45               -0.47                    0.01                   -0.46
## 2    0.38               -0.25                   -0.07                   -0.28
## 3    0.36                0.24                   -0.10                    0.19
## 4    0.47               -0.21                   -0.11                   -0.26
## 5    0.26               -0.63                    0.02                   -0.62
## 6    0.50                0.34                   -0.33                    0.17
## 7    0.45               -0.70                   -0.09                   -0.74
## 8    0.52               -0.52                   -0.17                   -0.60
## 9    0.43               -0.39                    0.04                   -0.37
## 10   0.31               -0.40                    0.02                   -0.39
## 11   0.51               -0.50                   -0.07                   -0.54
## 12   0.52               -0.54                   -0.02                   -0.55
## 13   0.48               -0.61                   -0.15                   -0.68
## 14   0.19                0.26                   -0.01                    0.25
## 15   0.40               -0.34                   -0.09                   -0.38
## 16   0.29                0.27                   -0.23                    0.16
## 17   0.37               -0.27                   -0.13                   -0.34
## 18   0.27               -0.10                   -0.10                   -0.15
## 19   0.42               -0.25                   -0.01                   -0.25
## 20   0.36                0.23                   -0.12                    0.17
## 21   0.38               -0.25                   -0.16                   -0.34
## 22   0.43               -1.04                   -0.09                   -1.09
## 23   0.25                0.23                   -0.15                    0.15
## 24   0.45               -0.41                   -0.12                   -0.47
## 25   0.40               -0.51                   -0.10                   -0.56
## 26   0.41               -0.10                   -0.06                   -0.13
## 27   0.26               -0.34                   -0.08                   -0.38
## 28   0.18               -0.03                   -0.23                   -0.14
## 29   0.50               -0.56                   -0.12                   -0.62
## 30   0.37               -0.35                   -0.08                   -0.39
## 31   0.31                0.08                   -0.15                    0.01
## 32   0.07                1.36                   -0.24                    1.24
## 33   0.22                0.47                   -0.13                    0.40
##    women_change_per_18_years cntry n_unique_essround
## 1                      -0.47    AT                 6
## 2                      -0.21    BE                10
## 3                       0.29    BG                 6
## 4                      -0.15    CH                10
## 5                      -0.64    CY                 5
## 6                       0.50    CZ                 9
## 7                      -0.66    DE                 9
## 8                      -0.43    DK                 8
## 9                      -0.41    EE                 9
## 10                     -0.41    ES                 9
## 11                     -0.46    FI                10
## 12                     -0.53    FR                10
## 13                     -0.53    GB                10
## 14                      0.27    GR                 5
## 15                     -0.29    HR                 4
## 16                      0.39    HU                10
## 17                     -0.21    IE                10
## 18                     -0.04    IL                 6
## 19                     -0.25    IS                 5
## 20                      0.29    IT                 4
## 21                     -0.17    LT                 6
## 22                     -1.00    LV                 2
## 23                      0.30    ME                 2
## 24                     -0.35    NL                10
## 25                     -0.47    NO                10
## 26                     -0.07    PL                 9
## 27                     -0.30    PT                10
## 28                      0.09    RU                 5
## 29                     -0.50    SE                 9
## 30                     -0.31    SI                10
## 31                      0.15    SK                 7
## 32                      1.48    TR                 2
## 33                      0.53    UA                 5
```

``` r
# rank by overall change
cntry_specific_changes %>%
  filter(n_unique_essround>4) %>%
  select(cntry,change_per_18_years) %>%
  arrange(change_per_18_years)
```

```
##    cntry change_per_18_years
## 1     DE               -0.70
## 2     CY               -0.63
## 3     GB               -0.61
## 4     SE               -0.56
## 5     FR               -0.54
## 6     DK               -0.52
## 7     NO               -0.51
## 8     FI               -0.50
## 9     AT               -0.47
## 10    NL               -0.41
## 11    ES               -0.40
## 12    EE               -0.39
## 13    SI               -0.35
## 14    PT               -0.34
## 15    IE               -0.27
## 16    BE               -0.25
## 17    IS               -0.25
## 18    LT               -0.25
## 19    CH               -0.21
## 20    IL               -0.10
## 21    PL               -0.10
## 22    RU               -0.03
## 23    SK                0.08
## 24    BG                0.24
## 25    GR                0.26
## 26    HU                0.27
## 27    CZ                0.34
## 28    UA                0.47
```

``` r
# rank by gendered change
cntry_specific_changes %>%
  filter(n_unique_essround>4) %>%
  select(cntry,gndr_change_per_18_year) %>%
  arrange(gndr_change_per_18_year)
```

```
##    cntry gndr_change_per_18_year
## 1     CZ                   -0.33
## 2     HU                   -0.23
## 3     RU                   -0.23
## 4     DK                   -0.17
## 5     LT                   -0.16
## 6     GB                   -0.15
## 7     SK                   -0.15
## 8     IE                   -0.13
## 9     UA                   -0.13
## 10    NL                   -0.12
## 11    SE                   -0.12
## 12    CH                   -0.11
## 13    BG                   -0.10
## 14    IL                   -0.10
## 15    NO                   -0.10
## 16    DE                   -0.09
## 17    PT                   -0.08
## 18    SI                   -0.08
## 19    BE                   -0.07
## 20    FI                   -0.07
## 21    PL                   -0.06
## 22    FR                   -0.02
## 23    GR                   -0.01
## 24    IS                   -0.01
## 25    AT                    0.01
## 26    CY                    0.02
## 27    ES                    0.02
## 28    EE                    0.04
```


#### Country-level predictors: GEI


``` r
mod6_GEI<-lmer(FM.z~gndr.c+essround.c+
             gndr.c:essround.c+
               gei.z.cm:gndr.c+gei.z.cm:essround.c+gei.z.cm:gndr.c:essround.c+
               (gndr.c+essround.c+gndr.c:essround.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_GEI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method [lmerModLmerTest]
## Formula: FM.z ~ gndr.c + essround.c + gndr.c:essround.c + gei.z.cm:gndr.c +  
##     gei.z.cm:essround.c + gei.z.cm:gndr.c:essround.c + (gndr.c +  
##     essround.c + gndr.c:essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##  170265.7  170428.2  -85114.9  170229.7     61558 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.8966 -0.5899  0.0212  0.5946  9.4022 
## 
## Random effects:
##  Groups   Name              Variance  Std.Dev. Corr             
##  cntry    (Intercept)       0.0840876 0.28998                   
##           gndr.c            0.0066515 0.08156   0.14            
##           essround.c        0.0015988 0.03998   0.19 -0.23      
##           gndr.c:essround.c 0.0001819 0.01349  -0.44 -0.19 -0.35
##  Residual                   0.9613844 0.98050                   
## Number of obs: 61576, groups:  cntry, 32
## 
## Fixed effects:
##                             Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)                 0.422983   0.051669 31.679510   8.186 2.56e-09 ***
## gndr.c                      0.372023   0.016695 30.220606  22.284  < 2e-16 ***
## essround.c                 -0.023276   0.007399 22.613065  -3.146  0.00459 ** 
## gndr.c:essround.c          -0.012316   0.003860 26.373022  -3.191  0.00364 ** 
## gndr.c:gei.z.cm             0.095852   0.017583 36.511061   5.452 3.61e-06 ***
## essround.c:gei.z.cm        -0.024788   0.007684 24.723309  -3.226  0.00352 ** 
## gndr.c:essround.c:gei.z.cm  0.001640   0.004374 32.940234   0.375  0.71008    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c essrn. gnd.:. gn.:.. es.:..
## gndr.c       0.117                                   
## essround.c   0.168 -0.199                            
## gndr.c:ssr. -0.274 -0.102 -0.208                     
## gndr.c:g.z.  0.001 -0.056  0.005 -0.070              
## essrnd.c:.. -0.010  0.005 -0.054  0.007 -0.221       
## gndr.c:.:..  0.001 -0.047  0.002 -0.207  0.120 -0.150
```

``` r
getFE(mod6_GEI,round=2,p.round = 5)
```

```
##                             Est.   SE    df     t       p    LL    UL
## (Intercept)                 0.42 0.05 31.68  8.19 0.00000  0.32  0.53
## gndr.c                      0.37 0.02 30.22 22.28 0.00000  0.34  0.41
## essround.c                 -0.02 0.01 22.61 -3.15 0.00459 -0.04 -0.01
## gndr.c:essround.c          -0.01 0.00 26.37 -3.19 0.00364 -0.02  0.00
## gndr.c:gei.z.cm             0.10 0.02 36.51  5.45 0.00000  0.06  0.13
## essround.c:gei.z.cm        -0.02 0.01 24.72 -3.23 0.00352 -0.04 -0.01
## gndr.c:essround.c:gei.z.cm  0.00 0.00 32.94  0.37 0.71008 -0.01  0.01
```

``` r
getVC(mod6_GEI)
```

```
##         grp              var1              var2 sdcor vcov
## 1     cntry       (Intercept)              <NA>  0.29 0.08
## 2     cntry            gndr.c              <NA>  0.08 0.01
## 3     cntry        essround.c              <NA>  0.04 0.00
## 4     cntry gndr.c:essround.c              <NA>  0.01 0.00
## 5     cntry       (Intercept)            gndr.c  0.14 0.00
## 6     cntry       (Intercept)        essround.c  0.19 0.00
## 7     cntry       (Intercept) gndr.c:essround.c -0.44 0.00
## 8     cntry            gndr.c        essround.c -0.23 0.00
## 9     cntry            gndr.c gndr.c:essround.c -0.19 0.00
## 10    cntry        essround.c gndr.c:essround.c -0.35 0.00
## 11 Residual              <NA>              <NA>  0.98 0.96
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_GEI,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] 36.8662
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_GEI,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] -9.902397
```

``` r
# GEI-levels  
gei_mean-gei_sd
```

```
## [1] 0.7907318
```

``` r
gei_mean+gei_sd
```

```
## [1] 0.9406642
```

``` r
# Simple slopes for 18 years
change_mod6_GEI<-emmeans(mod6_GEI,specs="essround.c",by="gei.z.cm",
                     at=list(gei.z.cm=c(-1,0,1),
                             gndr.c=0,
                             essround.c=rev(range(diff_dat$essround.c))),
                     disable.pbkrtest=T,
                     lmerTest.limit = 400000,infer=c(T,T))
```

```
## NOTE: Results may be misleading due to involvement in interactions
```

``` r
change_mod6_GEI
```

```
## gei.z.cm = -1:
##  essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.430 0.0756 37.0   0.2767    0.583   5.687  <.0001
##        -4.5  0.416 0.0670 28.9   0.2792    0.553   6.214  <.0001
## 
## gei.z.cm =  0:
##  essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.318 0.0660 25.8   0.1825    0.454   4.822  0.0001
##        -4.5  0.528 0.0566 21.8   0.4103    0.645   9.327  <.0001
## 
## gei.z.cm =  1:
##  essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.207 0.0734 32.9   0.0573    0.356   2.816  0.0082
##        -4.5  0.639 0.0656 27.4   0.5047    0.774   9.740  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_mod6_GEI,adjust="none",infer=c(T,T))
```

```
## gei.z.cm = -1:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   0.0136 0.0986 26.4   -0.189   0.2161   0.138  0.8912
## 
## gei.z.cm =  0:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)  -0.2095 0.0666 22.6   -0.347  -0.0716  -3.146  0.0046
## 
## gei.z.cm =  1:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)  -0.4326 0.0934 22.1   -0.626  -0.2390  -4.634  0.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# change in gender differences over time

change_in_diff_mod6_GEI<-emmeans(mod6_GEI,specs=c("gndr.c","essround.c"),by="gei.z.cm",
                             at=list(gndr.c=c(-0.5,0.5),
                                     gei.z.cm=c(-1,0,1),
                                     essround.c=rev(range(diff_dat$essround.c))),
                             disable.pbkrtest=T,
                             lmerTest.limit = 400000,infer=c(T,T))
change_in_diff_mod6_GEI
```

```
## gei.z.cm = -1:
##  gndr.c essround.c   emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5        4.5  0.32311 0.0820 37.5   0.1570    0.489   3.940  0.0003
##     0.5        4.5  0.53647 0.0740 38.0   0.3867    0.686   7.252  <.0001
##    -0.5       -4.5  0.24669 0.0667 25.9   0.1095    0.384   3.696  0.0010
##     0.5       -4.5  0.58566 0.0721 29.5   0.4384    0.733   8.126  <.0001
## 
## gei.z.cm =  0:
##  gndr.c essround.c   emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5        4.5  0.15994 0.0697 25.3   0.0164    0.303   2.294  0.0304
##     0.5        4.5  0.47654 0.0641 26.0   0.3448    0.608   7.434  <.0001
##    -0.5       -4.5  0.31400 0.0550 19.8   0.1992    0.429   5.710  <.0001
##     0.5       -4.5  0.74145 0.0608 22.0   0.6153    0.868  12.192  <.0001
## 
## gei.z.cm =  1:
##  gndr.c essround.c   emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5        4.5 -0.00322 0.0791 32.7  -0.1641    0.158  -0.041  0.9677
##     0.5        4.5  0.41661 0.0711 32.9   0.2719    0.561   5.860  <.0001
##    -0.5       -4.5  0.38131 0.0650 24.1   0.2472    0.515   5.865  <.0001
##     0.5       -4.5  0.89723 0.0705 27.8   0.7529    1.042  12.734  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod6_GEI,adjust="none",infer=c(T,T))
```

```
## gei.z.cm = -1:
##  contrast                                                 estimate     SE   df lower.CL upper.CL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.2134 0.0393 32.3  -0.2935 -0.13325
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   0.0764 0.1070 25.2  -0.1440  0.29680
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.2626 0.1060 26.6  -0.4812 -0.04390
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.2898 0.0966 26.8   0.0915  0.48809
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      -0.0492 0.0982 27.1  -0.2506  0.15220
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.3390 0.0368 40.0  -0.4133 -0.26460
##  t.ratio p.value
##   -5.422  <.0001
##    0.714  0.4819
##   -2.465  0.0204
##    2.999  0.0058
##   -0.501  0.6203
##   -9.212  <.0001
## 
## gei.z.cm =  0:
##  contrast                                                 estimate     SE   df lower.CL upper.CL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.3166 0.0228 25.6  -0.3636 -0.26964
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  -0.1541 0.0722 21.3  -0.3041 -0.00399
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.5815 0.0718 23.0  -0.7300 -0.43299
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.1625 0.0653 22.6   0.0272  0.29785
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      -0.2649 0.0652 22.6  -0.4000 -0.12984
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.4274 0.0253 29.4  -0.4791 -0.37575
##  t.ratio p.value
##  -13.869  <.0001
##   -2.133  0.0447
##   -8.099  <.0001
##    2.487  0.0207
##   -4.062  0.0005
##  -16.900  <.0001
## 
## gei.z.cm =  1:
##  contrast                                                 estimate     SE   df lower.CL upper.CL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.4198 0.0325 24.9  -0.4867 -0.35298
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  -0.3845 0.1010 20.6  -0.5941 -0.17500
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.9005 0.1010 22.5  -1.1099 -0.69101
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.0353 0.0912 22.2  -0.1537  0.22428
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      -0.4806 0.0916 21.9  -0.6707 -0.29053
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.5159 0.0339 30.7  -0.5852 -0.44666
##  t.ratio p.value
##  -12.936  <.0001
##   -3.821  0.0010
##   -8.904  <.0001
##    0.387  0.7024
##   -5.246  <.0001
##  -15.199  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# custom contrast to obtain the estimates for the entire span

changes_in_diff <- list(
  diff_ESS10 = c(-1,1,0,0),
  diff_ESS1 = c(0,0,-1,1))

diff_mod6_GEI<-contrast(change_in_diff_mod6_GEI,method = changes_in_diff,adjust="none",infer=c(T,T))

# gender differences at ESS1 and ESS10
diff_mod6_GEI
```

```
## gei.z.cm = -1:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10    0.213 0.0393 32.3    0.133    0.293   5.422  <.0001
##  diff_ESS1     0.339 0.0368 40.0    0.265    0.413   9.212  <.0001
## 
## gei.z.cm =  0:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10    0.317 0.0228 25.6    0.270    0.364  13.869  <.0001
##  diff_ESS1     0.427 0.0253 29.4    0.376    0.479  16.900  <.0001
## 
## gei.z.cm =  1:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10    0.420 0.0325 24.9    0.353    0.487  12.936  <.0001
##  diff_ESS1     0.516 0.0339 30.7    0.447    0.585  15.199  <.0001
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
##  contrast               estimate     SE   df lower.CL  upper.CL t.ratio p.value
##  diff_ESS10 - diff_ESS1  -0.1256 0.0576 37.2   -0.242 -0.008842  -2.179  0.0357
## 
## gei.z.cm =  0:
##  contrast               estimate     SE   df lower.CL  upper.CL t.ratio p.value
##  diff_ESS10 - diff_ESS1  -0.1108 0.0347 26.4   -0.182 -0.039490  -3.191  0.0036
## 
## gei.z.cm =  1:
##  contrast               estimate     SE   df lower.CL  upper.CL t.ratio p.value
##  diff_ESS10 - diff_ESS1  -0.0961 0.0468 25.2   -0.192  0.000285  -2.053  0.0506
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```


#### Country-level predictors: GGGI


``` r
mod6_GGGI<-lmer(FM.z~gndr.c+essround.c+
             gndr.c:essround.c+
               gggi.z.cm:gndr.c+gggi.z.cm:essround.c+gggi.z.cm:gndr.c:essround.c+
               (gndr.c+essround.c+gndr.c:essround.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_GGGI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method [lmerModLmerTest]
## Formula: FM.z ~ gndr.c + essround.c + gndr.c:essround.c + gggi.z.cm:gndr.c +  
##     gggi.z.cm:essround.c + gggi.z.cm:gndr.c:essround.c + (gndr.c +  
##     essround.c + gndr.c:essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##  124468.8  124625.6  -62216.4  124432.8     44989 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.8538 -0.5954  0.0231  0.5871  6.8335 
## 
## Random effects:
##  Groups   Name              Variance  Std.Dev. Corr             
##  cntry    (Intercept)       0.0858329 0.29297                   
##           gndr.c            0.0079734 0.08929  -0.16            
##           essround.c        0.0011199 0.03347   0.41 -0.31      
##           gndr.c:essround.c 0.0003116 0.01765  -0.76  0.39 -0.57
##  Residual                   0.9618149 0.98072                   
## Number of obs: 45007, groups:  cntry, 33
## 
## Fixed effects:
##                              Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)                  0.459800   0.051477 32.556127   8.932 2.86e-10 ***
## gndr.c                       0.363310   0.018342 29.341455  19.808  < 2e-16 ***
## essround.c                  -0.025700   0.006577 24.175238  -3.907 0.000659 ***
## gndr.c:essround.c           -0.007816   0.005379 25.704440  -1.453 0.158302    
## gndr.c:gggi.z.cm             0.074092   0.019314 34.764050   3.836 0.000504 ***
## essround.c:gggi.z.cm        -0.018275   0.006907 27.756113  -2.646 0.013264 *  
## gndr.c:essround.c:gggi.z.cm -0.007892   0.005505 34.609245  -1.434 0.160674    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c essrn. gnd.:. gn.:.. es.:..
## gndr.c      -0.133                                   
## essround.c   0.357 -0.244                            
## gndr.c:ssr. -0.436  0.138 -0.313                     
## gndr.c:gg..  0.000 -0.016  0.007 -0.039              
## essrnd.c:.. -0.016  0.009 -0.082  0.022 -0.202       
## gndr.c:.:..  0.002 -0.049  0.017 -0.083  0.058 -0.187
```

``` r
getFE(mod6_GGGI,round=2,p.round = 5)
```

```
##                              Est.   SE    df     t       p    LL    UL
## (Intercept)                  0.46 0.05 32.56  8.93 0.00000  0.36  0.56
## gndr.c                       0.36 0.02 29.34 19.81 0.00000  0.33  0.40
## essround.c                  -0.03 0.01 24.18 -3.91 0.00066 -0.04 -0.01
## gndr.c:essround.c           -0.01 0.01 25.70 -1.45 0.15830 -0.02  0.00
## gndr.c:gggi.z.cm             0.07 0.02 34.76  3.84 0.00050  0.03  0.11
## essround.c:gggi.z.cm        -0.02 0.01 27.76 -2.65 0.01326 -0.03  0.00
## gndr.c:essround.c:gggi.z.cm -0.01 0.01 34.61 -1.43 0.16067 -0.02  0.00
```

``` r
getVC(mod6_GGGI)
```

```
##         grp              var1              var2 sdcor vcov
## 1     cntry       (Intercept)              <NA>  0.29 0.09
## 2     cntry            gndr.c              <NA>  0.09 0.01
## 3     cntry        essround.c              <NA>  0.03 0.00
## 4     cntry gndr.c:essround.c              <NA>  0.02 0.00
## 5     cntry       (Intercept)            gndr.c -0.16 0.00
## 6     cntry       (Intercept)        essround.c  0.41 0.00
## 7     cntry       (Intercept) gndr.c:essround.c -0.76 0.00
## 8     cntry            gndr.c        essround.c -0.31 0.00
## 9     cntry            gndr.c gndr.c:essround.c  0.39 0.00
## 10    cntry        essround.c gndr.c:essround.c -0.57 0.00
## 11 Residual              <NA>              <NA>  0.98 0.96
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_GGGI,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] 55.77553
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_GGGI,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] -88.31373
```

``` r
# GGGI-levels  
gggi_mean-gggi_sd
```

```
## [1] 0.6827146
```

``` r
gggi_mean+gggi_sd
```

```
## [1] 0.7867524
```

``` r
# Simple slopes for 18 years
change_mod6_GGGI<-emmeans(mod6_GGGI,specs="essround.c",by="gggi.z.cm",
                     at=list(gggi.z.cm=c(-1,0,1),
                             gndr.c=0,
                             essround.c=rev(range(diff_dat$essround.c))),
                     disable.pbkrtest=T,
                     lmerTest.limit = 400000,infer=c(T,T))
```

```
## NOTE: Results may be misleading due to involvement in interactions
```

``` r
change_mod6_GGGI
```

```
## gggi.z.cm = -1:
##  essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.426 0.0760 36.0    0.272    0.581   5.607  <.0001
##        -4.5  0.493 0.0592 34.0    0.373    0.614   8.332  <.0001
## 
## gggi.z.cm =  0:
##  essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.344 0.0679 26.7    0.205    0.484   5.067  <.0001
##        -4.5  0.575 0.0494 24.3    0.474    0.677  11.653  <.0001
## 
## gggi.z.cm =  1:
##  essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.262 0.0733 33.3    0.113    0.411   3.572  0.0011
##        -4.5  0.658 0.0575 33.2    0.541    0.775  11.441  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_mod6_GGGI,adjust="none",infer=c(T,T))
```

```
## gggi.z.cm = -1:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)  -0.0668 0.0893 28.3   -0.250    0.116  -0.748  0.4604
## 
## gggi.z.cm =  0:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)  -0.2313 0.0592 24.2   -0.353   -0.109  -3.907  0.0007
## 
## gggi.z.cm =  1:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)  -0.3958 0.0822 27.5   -0.564   -0.227  -4.813  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# change in gender differences over time

change_in_diff_mod6_GGGI<-emmeans(mod6_GGGI,specs=c("gndr.c","essround.c"),by="gggi.z.cm",
                             at=list(gndr.c=c(-0.5,0.5),
                                     gggi.z.cm=c(-1,0,1),
                                     essround.c=rev(range(diff_dat$essround.c))),
                             disable.pbkrtest=T,
                             lmerTest.limit = 400000,infer=c(T,T))
change_in_diff_mod6_GGGI
```

```
## gggi.z.cm = -1:
##  gndr.c essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5        4.5 0.2816 0.0877 36.5   0.1039    0.459   3.213  0.0027
##     0.5        4.5 0.5712 0.0709 35.8   0.4274    0.715   8.061  <.0001
##    -0.5       -4.5 0.3488 0.0605 32.5   0.2256    0.472   5.766  <.0001
##     0.5       -4.5 0.6377 0.0650 35.2   0.5056    0.770   9.803  <.0001
## 
## gggi.z.cm =  0:
##  gndr.c essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5        4.5 0.1801 0.0768 26.1   0.0222    0.338   2.344  0.0269
##     0.5        4.5 0.5082 0.0620 26.0   0.3807    0.636   8.192  <.0001
##    -0.5       -4.5 0.3762 0.0481 22.8   0.2767    0.476   7.824  <.0001
##     0.5       -4.5 0.7747 0.0544 24.7   0.6625    0.887  14.228  <.0001
## 
## gggi.z.cm =  1:
##  gndr.c essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5        4.5 0.0786 0.0843 33.7  -0.0927    0.250   0.932  0.3578
##     0.5        4.5 0.4453 0.0678 32.9   0.3073    0.583   6.565  <.0001
##    -0.5       -4.5 0.4036 0.0584 31.4   0.2846    0.523   6.914  <.0001
##     0.5       -4.5 0.9117 0.0637 34.4   0.7824    1.041  14.324  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod6_GGGI,adjust="none",infer=c(T,T))
```

```
## gggi.z.cm = -1:
##  contrast                                                 estimate     SE   df lower.CL upper.CL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.2896 0.0477 31.4  -0.3869  -0.1922
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  -0.0672 0.1040 29.2  -0.2803   0.1459
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.3560 0.0987 30.3  -0.5576  -0.1545
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.2224 0.0874 27.8   0.0432   0.4015
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      -0.0665 0.0877 27.5  -0.2462   0.1132
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.2889 0.0420 35.7  -0.3740  -0.2037
##  t.ratio p.value
##   -6.066  <.0001
##   -0.645  0.5243
##   -3.607  0.0011
##    2.544  0.0168
##   -0.758  0.4546
##   -6.882  <.0001
## 
## gggi.z.cm =  0:
##  contrast                                                 estimate     SE   df lower.CL upper.CL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.3281 0.0323 25.2  -0.3947  -0.2616
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  -0.1961 0.0706 23.6  -0.3420  -0.0503
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.5946 0.0661 25.3  -0.7307  -0.4585
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.1320 0.0575 24.0   0.0133   0.2507
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      -0.2665 0.0565 24.7  -0.3829  -0.1500
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.3985 0.0283 27.6  -0.4564  -0.3405
##  t.ratio p.value
##  -10.151  <.0001
##   -2.777  0.0106
##   -8.993  <.0001
##    2.295  0.0308
##   -4.715  0.0001
##  -14.091  <.0001
## 
## gggi.z.cm =  1:
##  contrast                                                 estimate     SE   df lower.CL upper.CL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.3667 0.0435 29.2  -0.4557  -0.2777
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  -0.3251 0.0960 28.2  -0.5217  -0.1285
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.8332 0.0919 29.0  -1.0211  -0.6453
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.0416 0.0805 26.8  -0.1236   0.2069
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      -0.4665 0.0807 26.5  -0.6322  -0.3008
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.5081 0.0412 31.1  -0.5922  -0.4240
##  t.ratio p.value
##   -8.424  <.0001
##   -3.386  0.0021
##   -9.069  <.0001
##    0.517  0.6094
##   -5.781  <.0001
##  -12.320  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# custom contrast to obtain the estimates for the entire span

changes_in_diff <- list(
  diff_ESS10 = c(-1,1,0,0),
  diff_ESS1 = c(0,0,-1,1))

diff_mod6_GGGI<-contrast(change_in_diff_mod6_GGGI,method = changes_in_diff,adjust="none",infer=c(T,T))

# gender differences at ESS1 and ESS10
diff_mod6_GGGI
```

```
## gggi.z.cm = -1:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10    0.290 0.0477 31.4    0.192    0.387   6.066  <.0001
##  diff_ESS1     0.289 0.0420 35.7    0.204    0.374   6.882  <.0001
## 
## gggi.z.cm =  0:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10    0.328 0.0323 25.2    0.262    0.395  10.151  <.0001
##  diff_ESS1     0.398 0.0283 27.6    0.341    0.456  14.091  <.0001
## 
## gggi.z.cm =  1:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10    0.367 0.0435 29.2    0.278    0.456   8.424  <.0001
##  diff_ESS1     0.508 0.0412 31.1    0.424    0.592  12.320  <.0001
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
##  contrast                estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10 - diff_ESS1  0.000687 0.0721 35.7   -0.146   0.1469   0.010  0.9925
## 
## gggi.z.cm =  0:
##  contrast                estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10 - diff_ESS1 -0.070344 0.0484 25.7   -0.170   0.0292  -1.453  0.1583
## 
## gggi.z.cm =  1:
##  contrast                estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10 - diff_ESS1 -0.141375 0.0663 32.2   -0.276  -0.0063  -2.131  0.0408
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

#### Country-level predictors: GDI


``` r
mod6_GDI<-lmer(FM.z~gndr.c+essround.c+
             gndr.c:essround.c+
               gdi.z.cm:gndr.c+gdi.z.cm:essround.c+gdi.z.cm:gndr.c:essround.c+
               (gndr.c+essround.c+gndr.c:essround.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_GDI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method [lmerModLmerTest]
## Formula: FM.z ~ gndr.c + essround.c + gndr.c:essround.c + gdi.z.cm:gndr.c +  
##     gdi.z.cm:essround.c + gdi.z.cm:gndr.c:essround.c + (gndr.c +  
##     essround.c + gndr.c:essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##  175209.5  175372.4  -87586.7  175173.5     63029 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.8320 -0.5886  0.0218  0.5916  9.3415 
## 
## Random effects:
##  Groups   Name              Variance  Std.Dev. Corr             
##  cntry    (Intercept)       0.0861413 0.29350                   
##           gndr.c            0.0140651 0.11860  -0.45            
##           essround.c        0.0022845 0.04780   0.64 -0.63      
##           gndr.c:essround.c 0.0001597 0.01264  -0.57  0.11 -0.38
##  Residual                   0.9731183 0.98647                   
## Number of obs: 63047, groups:  cntry, 33
## 
## Fixed effects:
##                             Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)                 0.435674   0.051424 32.596306   8.472 9.56e-10 ***
## gndr.c                      0.368193   0.022278 30.455535  16.527  < 2e-16 ***
## essround.c                 -0.020731   0.008570 28.615098  -2.419 0.022156 *  
## gndr.c:essround.c          -0.011502   0.003662 26.914066  -3.141 0.004065 ** 
## gndr.c:gdi.z.cm             0.035256   0.021299 35.276477   1.655 0.106740    
## essround.c:gdi.z.cm        -0.027909   0.007118 31.772357  -3.921 0.000441 ***
## gndr.c:essround.c:gdi.z.cm  0.006173   0.004382 38.554966   1.409 0.166914    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c essrn. gnd.:. gn.:.. es.:..
## gndr.c      -0.412                                   
## essround.c   0.616 -0.568                            
## gndr.c:ssr. -0.345  0.053 -0.222                     
## gndr.c:gd.. -0.003 -0.013  0.002 -0.006              
## essrnd.c:..  0.001  0.005 -0.024  0.001 -0.401       
## gndr.c:.:..  0.003 -0.019  0.003  0.051  0.032 -0.022
```

``` r
getFE(mod6_GDI,round=2,p.round = 5)
```

```
##                             Est.   SE    df     t       p    LL    UL
## (Intercept)                 0.44 0.05 32.60  8.47 0.00000  0.33  0.54
## gndr.c                      0.37 0.02 30.46 16.53 0.00000  0.32  0.41
## essround.c                 -0.02 0.01 28.62 -2.42 0.02216 -0.04  0.00
## gndr.c:essround.c          -0.01 0.00 26.91 -3.14 0.00407 -0.02  0.00
## gndr.c:gdi.z.cm             0.04 0.02 35.28  1.66 0.10674 -0.01  0.08
## essround.c:gdi.z.cm        -0.03 0.01 31.77 -3.92 0.00044 -0.04 -0.01
## gndr.c:essround.c:gdi.z.cm  0.01 0.00 38.55  1.41 0.16691  0.00  0.02
```

``` r
getVC(mod6_GDI)
```

```
##         grp              var1              var2 sdcor  vcov
## 1     cntry       (Intercept)              <NA>  0.29  0.09
## 2     cntry            gndr.c              <NA>  0.12  0.01
## 3     cntry        essround.c              <NA>  0.05  0.00
## 4     cntry gndr.c:essround.c              <NA>  0.01  0.00
## 5     cntry       (Intercept)            gndr.c -0.45 -0.02
## 6     cntry       (Intercept)        essround.c  0.64  0.01
## 7     cntry       (Intercept) gndr.c:essround.c -0.57  0.00
## 8     cntry            gndr.c        essround.c -0.63  0.00
## 9     cntry            gndr.c gndr.c:essround.c  0.11  0.00
## 10    cntry        essround.c gndr.c:essround.c -0.38  0.00
## 11 Residual              <NA>              <NA>  0.99  0.97
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_GDI,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] 9.785392
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_GDI,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] 3.49522
```

``` r
# GDI-levels  
gdi_mean-gdi_sd
```

```
## [1] 0.9589057
```

``` r
gdi_mean+gdi_sd
```

```
## [1] 1.011161
```

``` r
# Simple slopes for 18 years
change_mod6_GDI<-emmeans(mod6_GDI,specs="essround.c",by="gdi.z.cm",
                     at=list(gdi.z.cm=c(-1,0,1),
                             gndr.c=0,
                             essround.c=rev(range(diff_dat$essround.c))),
                     disable.pbkrtest=T,
                     lmerTest.limit = 400000,infer=c(T,T))
```

```
## NOTE: Results may be misleading due to involvement in interactions
```

``` r
change_mod6_GDI
```

```
## gdi.z.cm = -1:
##  essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.468 0.0875 38.7   0.2910    0.645   5.349  <.0001
##        -4.5  0.403 0.0527 35.5   0.2964    0.510   7.651  <.0001
## 
## gdi.z.cm =  0:
##  essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.342 0.0811 30.6   0.1769    0.508   4.223  0.0002
##        -4.5  0.529 0.0411 28.0   0.4448    0.613  12.866  <.0001
## 
## gdi.z.cm =  1:
##  essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.217 0.0869 37.6   0.0409    0.393   2.496  0.0171
##        -4.5  0.655 0.0515 35.4   0.5500    0.759  12.707  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_mod6_GDI,adjust="none",infer=c(T,T))
```

```
## gdi.z.cm = -1:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   0.0646 0.1010 37.8   -0.141   0.2700   0.637  0.5281
## 
## gdi.z.cm =  0:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)  -0.1866 0.0771 28.6   -0.344  -0.0287  -2.419  0.0222
## 
## gdi.z.cm =  1:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)  -0.4378 0.0991 35.8   -0.639  -0.2368  -4.419  0.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# change in gender differences over time

change_in_diff_mod6_GDI<-emmeans(mod6_GDI,specs=c("gndr.c","essround.c"),by="gdi.z.cm",
                             at=list(gndr.c=c(-0.5,0.5),
                                     gdi.z.cm=c(-1,0,1),
                                     essround.c=rev(range(diff_dat$essround.c))),
                             disable.pbkrtest=T,
                             lmerTest.limit = 400000,infer=c(T,T))
change_in_diff_mod6_GDI
```

```
## gdi.z.cm = -1:
##  gndr.c essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5        4.5 0.3413 0.0988 39.3  0.14147    0.541   3.454  0.0013
##     0.5        4.5 0.5947 0.0799 38.4  0.43295    0.756   7.441  <.0001
##    -0.5       -4.5 0.1971 0.0521 33.0  0.09106    0.303   3.781  0.0006
##     0.5       -4.5 0.6096 0.0599 33.9  0.48779    0.731  10.172  <.0001
## 
## gdi.z.cm =  0:
##  gndr.c essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5        4.5 0.1842 0.0904 30.3 -0.00032    0.369   2.038  0.0504
##     0.5        4.5 0.5006 0.0734 29.9  0.35075    0.650   6.824  <.0001
##    -0.5       -4.5 0.3190 0.0413 26.4  0.23417    0.404   7.725  <.0001
##     0.5       -4.5 0.7389 0.0452 28.0  0.64641    0.831  16.359  <.0001
## 
## gdi.z.cm =  1:
##  gndr.c essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5        4.5 0.0271 0.0982 38.3 -0.17160    0.226   0.276  0.7843
##     0.5        4.5 0.4065 0.0793 37.5  0.24592    0.567   5.127  <.0001
##    -0.5       -4.5 0.4408 0.0513 34.1  0.33658    0.545   8.593  <.0001
##     0.5       -4.5 0.8683 0.0589 34.6  0.74866    0.988  14.743  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod6_GDI,adjust="none",infer=c(T,T))
```

```
## gdi.z.cm = -1:
##  contrast                                                 estimate     SE   df lower.CL upper.CL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.2534 0.0409 29.6  -0.3371  -0.1697
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   0.1441 0.1070 37.1  -0.0736   0.3619
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.2683 0.1200 38.6  -0.5105  -0.0262
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.3975 0.0905 36.1   0.2141   0.5810
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      -0.0149 0.1010 35.6  -0.2207   0.1909
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.4125 0.0388 32.2  -0.4914  -0.3335
##  t.ratio p.value
##   -6.189  <.0001
##    1.341  0.1880
##   -2.242  0.0308
##    4.394  0.0001
##   -0.147  0.8837
##  -10.639  <.0001
## 
## gdi.z.cm =  0:
##  contrast                                                 estimate     SE   df lower.CL upper.CL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.3164 0.0284 22.8  -0.3752  -0.2577
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  -0.1348 0.0824 27.9  -0.3036   0.0339
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.5548 0.0916 29.4  -0.7421  -0.3675
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.1816 0.0670 27.7   0.0443   0.3190
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      -0.2383 0.0752 27.9  -0.3924  -0.0843
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.4200 0.0270 32.0  -0.4750  -0.3650
##  t.ratio p.value
##  -11.142  <.0001
##   -1.637  0.1129
##   -6.054  <.0001
##    2.710  0.0114
##   -3.169  0.0037
##  -15.553  <.0001
## 
## gdi.z.cm =  1:
##  contrast                                                 estimate     SE   df lower.CL upper.CL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.3795 0.0409 33.2  -0.4627  -0.2962
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  -0.4138 0.1050 35.9  -0.6276  -0.1999
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.8412 0.1170 36.7  -1.0791  -0.6033
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)     -0.0343 0.0879 34.2  -0.2129   0.1443
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      -0.4617 0.0995 34.9  -0.6638  -0.2597
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.4274 0.0398 39.1  -0.5080  -0.3469
##  t.ratio p.value
##   -9.273  <.0001
##   -3.924  0.0004
##   -7.166  <.0001
##   -0.390  0.6987
##   -4.641  <.0001
##  -10.729  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# custom contrast to obtain the estimates for the entire span

changes_in_diff <- list(
  diff_ESS10 = c(-1,1,0,0),
  diff_ESS1 = c(0,0,-1,1))

diff_mod6_GDI<-contrast(change_in_diff_mod6_GDI,method = changes_in_diff,adjust="none",infer=c(T,T))

# gender differences at ESS1 and ESS10
diff_mod6_GDI
```

```
## gdi.z.cm = -1:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10    0.253 0.0409 29.6    0.170    0.337   6.189  <.0001
##  diff_ESS1     0.412 0.0388 32.2    0.334    0.491  10.639  <.0001
## 
## gdi.z.cm =  0:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10    0.316 0.0284 22.8    0.258    0.375  11.142  <.0001
##  diff_ESS1     0.420 0.0270 32.0    0.365    0.475  15.553  <.0001
## 
## gdi.z.cm =  1:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10    0.379 0.0409 33.2    0.296    0.463   9.273  <.0001
##  diff_ESS1     0.427 0.0398 39.1    0.347    0.508  10.729  <.0001
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
##  diff_ESS10 - diff_ESS1   -0.159 0.0501 29.5   -0.261  -0.0567  -3.176  0.0035
## 
## gdi.z.cm =  0:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10 - diff_ESS1   -0.104 0.0330 26.9   -0.171  -0.0359  -3.141  0.0041
## 
## gdi.z.cm =  1:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10 - diff_ESS1   -0.048 0.0527 42.6   -0.154   0.0583  -0.910  0.3677
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

#### Country-level predictors: log_GDP


``` r
mod6_log_GDP<-lmer(FM.z~gndr.c+essround.c+
             gndr.c:essround.c+
               log_gdp.z.cm:gndr.c+log_gdp.z.cm:essround.c+log_gdp.z.cm:gndr.c:essround.c+
               (gndr.c+essround.c+gndr.c:essround.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_log_GDP)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method [lmerModLmerTest]
## Formula: FM.z ~ gndr.c + essround.c + gndr.c:essround.c + log_gdp.z.cm:gndr.c +  
##     log_gdp.z.cm:essround.c + log_gdp.z.cm:gndr.c:essround.c +  
##     (gndr.c + essround.c + gndr.c:essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##  175211.2  175374.1  -87587.6  175175.2     63029 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.8316 -0.5885  0.0217  0.5923  9.3395 
## 
## Random effects:
##  Groups   Name              Variance  Std.Dev. Corr             
##  cntry    (Intercept)       0.0846977 0.29103                   
##           gndr.c            0.0088696 0.09418  -0.06            
##           essround.c        0.0017459 0.04178   0.25 -0.40      
##           gndr.c:essround.c 0.0001783 0.01335  -0.54  0.10 -0.52
##  Residual                   0.9731343 0.98648                   
## Number of obs: 63047, groups:  cntry, 33
## 
## Fixed effects:
##                                 Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)                     0.430328   0.051056 32.534072   8.429 1.09e-09 ***
## gndr.c                          0.368056   0.018372 30.518370  20.034  < 2e-16 ***
## essround.c                     -0.020992   0.007588 23.425489  -2.766  0.01087 *  
## gndr.c:essround.c              -0.010503   0.003813 29.903840  -2.754  0.00991 ** 
## gndr.c:log_gdp.z.cm             0.072694   0.018820 32.243106   3.863  0.00051 ***
## essround.c:log_gdp.z.cm        -0.023175   0.007651 24.230907  -3.029  0.00576 ** 
## gndr.c:essround.c:log_gdp.z.cm -0.004642   0.003867 31.816837  -1.201  0.23877    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c essrn. gnd.:. g.:_.. e.:_..
## gndr.c      -0.052                                   
## essround.c   0.232 -0.351                            
## gndr.c:ssr. -0.332  0.061 -0.305                     
## gndr.c:l_.. -0.002 -0.034 -0.001 -0.015              
## essrnd.:_..  0.011 -0.001 -0.040 -0.005 -0.340       
## gndr.:.:_..  0.001 -0.013 -0.001 -0.212  0.068 -0.228
```

``` r
getFE(mod6_log_GDP,round=2,p.round = 5)
```

```
##                                 Est.   SE    df     t       p    LL    UL
## (Intercept)                     0.43 0.05 32.53  8.43 0.00000  0.33  0.53
## gndr.c                          0.37 0.02 30.52 20.03 0.00000  0.33  0.41
## essround.c                     -0.02 0.01 23.43 -2.77 0.01087 -0.04 -0.01
## gndr.c:essround.c              -0.01 0.00 29.90 -2.75 0.00991 -0.02  0.00
## gndr.c:log_gdp.z.cm             0.07 0.02 32.24  3.86 0.00051  0.03  0.11
## essround.c:log_gdp.z.cm        -0.02 0.01 24.23 -3.03 0.00576 -0.04 -0.01
## gndr.c:essround.c:log_gdp.z.cm  0.00 0.00 31.82 -1.20 0.23877 -0.01  0.00
```

``` r
getVC(mod6_log_GDP)
```

```
##         grp              var1              var2 sdcor vcov
## 1     cntry       (Intercept)              <NA>  0.29 0.08
## 2     cntry            gndr.c              <NA>  0.09 0.01
## 3     cntry        essround.c              <NA>  0.04 0.00
## 4     cntry gndr.c:essround.c              <NA>  0.01 0.00
## 5     cntry       (Intercept)            gndr.c -0.06 0.00
## 6     cntry       (Intercept)        essround.c  0.25 0.00
## 7     cntry       (Intercept) gndr.c:essround.c -0.54 0.00
## 8     cntry            gndr.c        essround.c -0.40 0.00
## 9     cntry            gndr.c gndr.c:essround.c  0.10 0.00
## 10    cntry        essround.c gndr.c:essround.c -0.52 0.00
## 11 Residual              <NA>              <NA>  0.99 0.97
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_log_GDP,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] 31.05432
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_log_GDP,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] -7.774407
```

``` r
# log_GDP-levels  
log_gdp_mean-log_gdp_sd
```

```
## [1] 10.22564
```

``` r
log_gdp_mean+log_gdp_sd
```

```
## [1] 11.02399
```

``` r
# Simple slopes for 18 years
change_mod6_log_GDP<-emmeans(mod6_log_GDP,specs="essround.c",by="log_gdp.z.cm",
                     at=list(log_gdp.z.cm=c(-1,0,1),
                             gndr.c=0,
                             essround.c=rev(range(diff_dat$essround.c))),
                     disable.pbkrtest=T,
                     lmerTest.limit = 400000,infer=c(T,T))
```

```
## NOTE: Results may be misleading due to involvement in interactions
```

``` r
change_mod6_log_GDP
```

```
## log_gdp.z.cm = -1:
##  essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.440 0.0763 35.1   0.2853    0.595   5.769  <.0001
##        -4.5  0.421 0.0655 30.0   0.2868    0.554   6.424  <.0001
## 
## log_gdp.z.cm =  0:
##  essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.336 0.0677 26.3   0.1968    0.475   4.962  <.0001
##        -4.5  0.525 0.0544 22.3   0.4120    0.638   9.639  <.0001
## 
## log_gdp.z.cm =  1:
##  essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.232 0.0756 33.1   0.0779    0.385   3.064  0.0043
##        -4.5  0.629 0.0634 28.1   0.4993    0.759   9.928  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_mod6_log_GDP,adjust="none",infer=c(T,T))
```

```
## log_gdp.z.cm = -1:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   0.0196 0.0989 26.5   -0.184   0.2228   0.199  0.8441
## 
## log_gdp.z.cm =  0:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)  -0.1889 0.0683 23.4   -0.330  -0.0478  -2.766  0.0109
## 
## log_gdp.z.cm =  1:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)  -0.3975 0.0950 22.8   -0.594  -0.2009  -4.184  0.0004
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# change in gender differences over time

change_in_diff_mod6_log_GDP<-emmeans(mod6_log_GDP,specs=c("gndr.c","essround.c"),by="log_gdp.z.cm",
                             at=list(gndr.c=c(-0.5,0.5),
                                     log_gdp.z.cm=c(-1,0,1),
                                     essround.c=rev(range(diff_dat$essround.c))),
                             disable.pbkrtest=T,
                             lmerTest.limit = 400000,infer=c(T,T))
change_in_diff_mod6_log_GDP
```

```
## log_gdp.z.cm = -1:
##  gndr.c essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5        4.5 0.3057 0.0854 35.1   0.1323    0.479   3.578  0.0010
##     0.5        4.5 0.5746 0.0715 35.7   0.4295    0.720   8.032  <.0001
##    -0.5       -4.5 0.2596 0.0652 27.2   0.1260    0.393   3.984  0.0005
##     0.5       -4.5 0.5814 0.0706 30.2   0.4371    0.726   8.230  <.0001
## 
## log_gdp.z.cm =  0:
##  gndr.c essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5        4.5 0.1755 0.0741 25.7   0.0231    0.328   2.369  0.0257
##     0.5        4.5 0.4963 0.0633 26.6   0.3662    0.626   7.837  <.0001
##    -0.5       -4.5 0.3171 0.0532 20.5   0.2064    0.428   5.962  <.0001
##     0.5       -4.5 0.7325 0.0583 22.3   0.6117    0.853  12.568  <.0001
## 
## log_gdp.z.cm =  1:
##  gndr.c essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5        4.5 0.0453 0.0843 32.6  -0.1264    0.217   0.537  0.5950
##     0.5        4.5 0.4179 0.0701 32.6   0.2752    0.561   5.964  <.0001
##    -0.5       -4.5 0.3746 0.0626 24.9   0.2456    0.504   5.982  <.0001
##     0.5       -4.5 0.8835 0.0682 28.1   0.7439    1.023  12.963  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod6_log_GDP,adjust="none",infer=c(T,T))
```

```
## log_gdp.z.cm = -1:
##  contrast                                                 estimate     SE   df lower.CL upper.CL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5     -0.26899 0.0392 31.9  -0.3489  -0.1890
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  0.04602 0.1080 26.4  -0.1767   0.2687
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)  -0.27572 0.1110 27.2  -0.5028  -0.0486
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)     0.31501 0.0935 26.2   0.1229   0.5071
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)     -0.00673 0.0962 25.5  -0.2048   0.1913
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5) -0.32173 0.0366 38.1  -0.3958  -0.2477
##  t.ratio p.value
##   -6.855  <.0001
##    0.424  0.6747
##   -2.490  0.0192
##    3.369  0.0023
##   -0.070  0.9448
##   -8.798  <.0001
## 
## log_gdp.z.cm =  0:
##  contrast                                                 estimate     SE   df lower.CL upper.CL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5     -0.32079 0.0259 24.3  -0.3742  -0.2674
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5) -0.14167 0.0753 22.9  -0.2975   0.0142
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)  -0.55698 0.0767 24.1  -0.7152  -0.3987
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)     0.17913 0.0642 23.2   0.0464   0.3118
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)     -0.23619 0.0651 22.5  -0.3711  -0.1012
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5) -0.41532 0.0244 30.6  -0.4650  -0.3656
##  t.ratio p.value
##  -12.388  <.0001
##   -1.881  0.0727
##   -7.262  <.0001
##    2.791  0.0103
##   -3.626  0.0015
##  -17.048  <.0001
## 
## log_gdp.z.cm =  1:
##  contrast                                                 estimate     SE   df lower.CL upper.CL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5     -0.37260 0.0347 21.6  -0.4447  -0.3005
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5) -0.32935 0.1040 22.1  -0.5445  -0.1142
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)  -0.83825 0.1070 23.6  -1.0594  -0.6171
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)     0.04325 0.0890 22.2  -0.1413   0.2278
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)     -0.46565 0.0907 21.0  -0.6542  -0.2771
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5) -0.50890 0.0328 28.4  -0.5760  -0.4418
##  t.ratio p.value
##  -10.735  <.0001
##   -3.173  0.0044
##   -7.830  <.0001
##    0.486  0.6319
##   -5.136  <.0001
##  -15.530  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# custom contrast to obtain the estimates for the entire span

changes_in_diff <- list(
  diff_ESS10 = c(-1,1,0,0),
  diff_ESS1 = c(0,0,-1,1))

diff_mod6_log_GDP<-contrast(change_in_diff_mod6_log_GDP,method = changes_in_diff,adjust="none",infer=c(T,T))

# gender differences at ESS1 and ESS10
diff_mod6_log_GDP
```

```
## log_gdp.z.cm = -1:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10    0.269 0.0392 31.9    0.189    0.349   6.855  <.0001
##  diff_ESS1     0.322 0.0366 38.1    0.248    0.396   8.798  <.0001
## 
## log_gdp.z.cm =  0:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10    0.321 0.0259 24.3    0.267    0.374  12.388  <.0001
##  diff_ESS1     0.415 0.0244 30.6    0.366    0.465  17.048  <.0001
## 
## log_gdp.z.cm =  1:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10    0.373 0.0347 21.6    0.301    0.445  10.735  <.0001
##  diff_ESS1     0.509 0.0328 28.4    0.442    0.576  15.530  <.0001
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
##  diff_ESS10 - diff_ESS1  -0.0527 0.0538 43.7   -0.161   0.0557  -0.980  0.3323
## 
## log_gdp.z.cm =  0:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10 - diff_ESS1  -0.0945 0.0343 29.9   -0.165  -0.0244  -2.754  0.0099
## 
## log_gdp.z.cm =  1:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10 - diff_ESS1  -0.1363 0.0434 23.5   -0.226  -0.0466  -3.141  0.0045
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```


## mod7_GEI and mod8_GEI: Countries' progress in GEI and convergence of gender gap

* First fit a model including country average of gender-equality and its interaction with gender and time and differential development for women and men (fixed and random)

* Then fit a model that adds time-specific within-country fluctuations from country's average gender-equality (fixed and random) and compare to a model without these parameters.


``` r
mod7_GEI<-lmer(FM.z~gndr.c+gei.z.cm+gndr.c:gei.z.cm+
                      essround.c+
                      gndr.c:essround.c+
                      (gndr.c+essround.c+gndr.c:essround.c|cntry),
               data=diff_dat,REML=F,weights = pspwght,
                    control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

summary(mod7_GEI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method [lmerModLmerTest]
## Formula: FM.z ~ gndr.c + gei.z.cm + gndr.c:gei.z.cm + essround.c + gndr.c:essround.c +  
##     (gndr.c + essround.c + gndr.c:essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##  170263.5  170417.0  -85114.8  170229.5     61559 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.8925 -0.5898  0.0216  0.5946  9.4049 
## 
## Random effects:
##  Groups   Name              Variance  Std.Dev. Corr             
##  cntry    (Intercept)       0.0579036 0.24063                   
##           gndr.c            0.0066635 0.08163   0.09            
##           essround.c        0.0025360 0.05036   0.16 -0.26      
##           gndr.c:essround.c 0.0001987 0.01409  -0.35 -0.12 -0.48
##  Residual                   0.9613422 0.98048                   
## Number of obs: 61576, groups:  cntry, 32
## 
## Fixed effects:
##                    Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)        0.421457   0.043064 31.523579   9.787 4.49e-11 ***
## gndr.c             0.372494   0.016691 29.265581  22.317  < 2e-16 ***
## gei.z.cm          -0.143711   0.042630 30.849089  -3.371  0.00203 ** 
## essround.c        -0.023702   0.009173 24.464251  -2.584  0.01615 *  
## gndr.c:gei.z.cm    0.077749   0.017233 37.097552   4.512 6.28e-05 ***
## gndr.c:essround.c -0.011981   0.003849 24.569590  -3.112  0.00466 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c g.z.cm essrn. gn.:..
## gndr.c       0.080                            
## gei.z.cm    -0.009 -0.007                     
## essround.c   0.143 -0.226 -0.007              
## gndr.c:g.z. -0.001 -0.053  0.129 -0.005       
## gndr.c:ssr. -0.229 -0.077 -0.031 -0.309 -0.057
```

``` r
getFE(mod7_GEI)
```

```
##                    Est.   SE    df     t     p    LL    UL
## (Intercept)        0.42 0.04 31.52  9.79 0.000  0.33  0.51
## gndr.c             0.37 0.02 29.27 22.32 0.000  0.34  0.41
## gei.z.cm          -0.14 0.04 30.85 -3.37 0.002 -0.23 -0.06
## essround.c        -0.02 0.01 24.46 -2.58 0.016 -0.04  0.00
## gndr.c:gei.z.cm    0.08 0.02 37.10  4.51 0.000  0.04  0.11
## gndr.c:essround.c -0.01 0.00 24.57 -3.11 0.005 -0.02  0.00
```

``` r
getVC(mod7_GEI)
```

```
##         grp              var1              var2 sdcor vcov
## 1     cntry       (Intercept)              <NA>  0.24 0.06
## 2     cntry            gndr.c              <NA>  0.08 0.01
## 3     cntry        essround.c              <NA>  0.05 0.00
## 4     cntry gndr.c:essround.c              <NA>  0.01 0.00
## 5     cntry       (Intercept)            gndr.c  0.09 0.00
## 6     cntry       (Intercept)        essround.c  0.16 0.00
## 7     cntry       (Intercept) gndr.c:essround.c -0.35 0.00
## 8     cntry            gndr.c        essround.c -0.26 0.00
## 9     cntry            gndr.c gndr.c:essround.c -0.12 0.00
## 10    cntry        essround.c gndr.c:essround.c -0.48 0.00
## 11 Residual              <NA>              <NA>  0.98 0.96
```

``` r
anova(mod2_GEI,mod7_GEI)
```

```
## Data: diff_dat
## Models:
## mod2_GEI: FM.z ~ gndr.c + gei.z.cm + gndr.c:gei.z.cm + (gndr.c | cntry)
## mod7_GEI: FM.z ~ gndr.c + gei.z.cm + gndr.c:gei.z.cm + essround.c + gndr.c:essround.c + (gndr.c + essround.c + gndr.c:essround.c | cntry)
##          npar    AIC    BIC logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod2_GEI    8 171249 171321 -85616    171233                         
## mod7_GEI   17 170264 170417 -85115    170230 1003.2  9  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

### Examine trends


``` r
change_mod7_GEI<-emmeans(mod7_GEI,specs="essround.c",by="gndr.c",
                              at=list(gndr.c=c(-0.5,0.5),
                                      essround.c=rev(range(diff_dat$essround.c)),
                                      gei.z.cmc=0,
                                      gei.z.cm=0),
                              disable.pbkrtest=T,
                              lmerTest.limit = 400000,infer=c(T,T),adjust="none")
change_mod7_GEI
```

```
## gndr.c = -0.5:
##  essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.156 0.0685 24.4   0.0143    0.297   2.271  0.0322
##        -4.5  0.315 0.0552 19.3   0.1996    0.430   5.708  <.0001
## 
## gndr.c =  0.5:
##  essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.474 0.0609 24.2   0.3484    0.600   7.783  <.0001
##        -4.5  0.741 0.0580 22.1   0.6210    0.862  12.771  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_mod7_GEI,adjust="none",infer=c(T,T))
```

```
## gndr.c = -0.5:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   -0.159 0.0894 23.7   -0.344   0.0253  -1.782  0.0875
## 
## gndr.c =  0.5:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   -0.267 0.0789 24.1   -0.430  -0.1043  -3.385  0.0024
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# change in gender differences over time

change_in_diff_mod7_GEI<-emmeans(mod7_GEI,specs=c("gndr.c","essround.c"),
                                      at=list(gndr.c=c(-0.5,0.5),
                                              essround.c=rev(range(diff_dat$essround.c)),
                                              gei.z.cmc=0,
                                              gei.z.cm=0),
                                      disable.pbkrtest=T,
                                      lmerTest.limit = 400000,infer=c(T,T),adjust="none")
change_in_diff_mod7_GEI
```

```
##  gndr.c essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5        4.5  0.156 0.0685 24.4   0.0143    0.297   2.271  0.0322
##     0.5        4.5  0.474 0.0609 24.2   0.3484    0.600   7.783  <.0001
##    -0.5       -4.5  0.315 0.0552 19.3   0.1996    0.430   5.708  <.0001
##     0.5       -4.5  0.741 0.0580 22.1   0.6210    0.862  12.771  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod7_GEI,adjust="none",infer=c(T,T))
```

```
##  contrast                                                 estimate     SE   df lower.CL upper.CL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5       -0.319 0.0231 21.4 -0.36658  -0.2706
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   -0.159 0.0894 23.7 -0.34409   0.0253
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)    -0.586 0.0878 23.7 -0.76723  -0.4044
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)       0.159 0.0805 23.9 -0.00692   0.3253
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)       -0.267 0.0789 24.1 -0.43014  -0.1043
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)   -0.426 0.0250 31.1 -0.47732  -0.3755
##  t.ratio p.value
##  -13.784  <.0001
##   -1.782  0.0875
##   -6.669  <.0001
##    1.979  0.0595
##   -3.385  0.0024
##  -17.082  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# custom contrast to obtain the estimates for the entire span

changes_in_diff <- list(
  diff_ESS10 = c(-1,1,0,0),
  diff_ESS1 = c(0,0,-1,1))

diff_mod7_GEI<-contrast(change_in_diff_mod7_GEI,method = changes_in_diff,adjust="none",infer=c(T,T))

# gender differences at ESS1 and ESS10
diff_mod7_GEI
```

```
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10    0.319 0.0231 21.4    0.271    0.367  13.784  <.0001
##  diff_ESS1     0.426 0.0250 31.1    0.376    0.477  17.082  <.0001
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
##  diff_ESS10 - diff_ESS1   -0.108 0.0346 24.6   -0.179  -0.0364  -3.112  0.0047
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```


``` r
mod8_GEI<-lmer(FM.z~gndr.c+gei.z.cm+gndr.c:gei.z.cm+
                      essround.c+
                      gndr.c:essround.c+
                      gei.z.cmc+gei.z.cmc:gndr.c+
                      (gndr.c+essround.c+gndr.c:essround.c+gei.z.cmc+gei.z.cmc:gndr.c|cntry),
               data=diff_dat,REML=F,weights = pspwght,
                    control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
```

```
## boundary (singular) fit: see help('isSingular')
```

``` r
getFE(mod8_GEI)
```

```
##                    Est.   SE    df     t     p    LL    UL
## (Intercept)        0.42 0.04 31.97  9.97 0.000  0.34  0.51
## gndr.c             0.37 0.02 29.87 20.09 0.000  0.34  0.41
## gei.z.cm          -0.16 0.04 32.67 -3.70 0.001 -0.24 -0.07
## essround.c        -0.02 0.01 25.02 -1.90 0.068 -0.05  0.00
## gei.z.cmc         -0.07 0.08 20.87 -0.85 0.407 -0.23  0.10
## gndr.c:gei.z.cm    0.07 0.01 42.80  5.00 0.000  0.04  0.10
## gndr.c:essround.c -0.01 0.01 23.39 -1.38 0.181 -0.02  0.00
## gndr.c:gei.z.cmc   0.00 0.04 32.35  0.05 0.958 -0.09  0.09
```

``` r
getVC(mod8_GEI)
```

```
##         grp              var1              var2 sdcor  vcov
## 1     cntry       (Intercept)              <NA>  0.24  0.06
## 2     cntry            gndr.c              <NA>  0.09  0.01
## 3     cntry        essround.c              <NA>  0.06  0.00
## 4     cntry         gei.z.cmc              <NA>  0.38  0.14
## 5     cntry gndr.c:essround.c              <NA>  0.02  0.00
## 6     cntry  gndr.c:gei.z.cmc              <NA>  0.12  0.01
## 7     cntry       (Intercept)            gndr.c  0.04  0.00
## 8     cntry       (Intercept)        essround.c  0.03  0.00
## 9     cntry       (Intercept)         gei.z.cmc  0.10  0.01
## 10    cntry       (Intercept) gndr.c:essround.c -0.13  0.00
## 11    cntry       (Intercept)  gndr.c:gei.z.cmc  0.02  0.00
## 12    cntry            gndr.c        essround.c -0.10  0.00
## 13    cntry            gndr.c         gei.z.cmc -0.14  0.00
## 14    cntry            gndr.c gndr.c:essround.c -0.94  0.00
## 15    cntry            gndr.c  gndr.c:gei.z.cmc  0.97  0.01
## 16    cntry        essround.c         gei.z.cmc -0.86 -0.02
## 17    cntry        essround.c gndr.c:essround.c -0.05  0.00
## 18    cntry        essround.c  gndr.c:gei.z.cmc -0.30  0.00
## 19    cntry         gei.z.cmc gndr.c:essround.c  0.26  0.00
## 20    cntry         gei.z.cmc  gndr.c:gei.z.cmc  0.03  0.00
## 21    cntry gndr.c:essround.c  gndr.c:gei.z.cmc -0.92  0.00
## 22 Residual              <NA>              <NA>  0.98  0.96
```

``` r
anova(mod2_GEI,mod7_GEI,mod8_GEI)
```

```
## Data: diff_dat
## Models:
## mod2_GEI: FM.z ~ gndr.c + gei.z.cm + gndr.c:gei.z.cm + (gndr.c | cntry)
## mod7_GEI: FM.z ~ gndr.c + gei.z.cm + gndr.c:gei.z.cm + essround.c + gndr.c:essround.c + (gndr.c + essround.c + gndr.c:essround.c | cntry)
## mod8_GEI: FM.z ~ gndr.c + gei.z.cm + gndr.c:gei.z.cm + essround.c + gndr.c:essround.c + gei.z.cmc + gei.z.cmc:gndr.c + (gndr.c + essround.c + gndr.c:essround.c + gei.z.cmc + gei.z.cmc:gndr.c | cntry)
##          npar    AIC    BIC logLik -2*log(L)    Chisq Df Pr(>Chisq)    
## mod2_GEI    8 171249 171321 -85616    171233                           
## mod7_GEI   17 170264 170417 -85115    170230 1003.193  9  < 2.2e-16 ***
## mod8_GEI   30 170229 170500 -85085    170169   60.381 13   4.49e-08 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

### Examine trends


``` r
# gender specific change over time

change_mod8_GEI<-emmeans(mod8_GEI,specs="essround.c",by="gndr.c",
                              at=list(gndr.c=c(-0.5,0.5),
                                      essround.c=rev(range(diff_dat$essround.c)),
                                      gei.z.cmc=0,
                                      gei.z.cm=0),
                              disable.pbkrtest=T,
                              lmerTest.limit = 400000,infer=c(T,T),adjust="none")
change_mod8_GEI
```

```
## gndr.c = -0.5:
##  essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.157 0.0705 28.2   0.0123    0.301   2.222  0.0345
##        -4.5  0.318 0.0710 26.1   0.1719    0.464   4.475  0.0001
## 
## gndr.c =  0.5:
##  essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.488 0.0673 28.7   0.3507    0.626   7.257  <.0001
##        -4.5  0.733 0.0726 27.1   0.5838    0.882  10.094  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_mod8_GEI,adjust="none",infer=c(T,T))
```

```
## gndr.c = -0.5:
##  contrast                         estimate    SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   -0.161 0.112 24.4   -0.393   0.0701  -1.437  0.1634
## 
## gndr.c =  0.5:
##  contrast                         estimate    SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   -0.244 0.109 25.6   -0.469  -0.0199  -2.239  0.0341
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# change in gender differences over time

change_in_diff_mod8_GEI<-emmeans(mod8_GEI,specs=c("gndr.c","essround.c"),
                                      at=list(gndr.c=c(-0.5,0.5),
                                              essround.c=rev(range(diff_dat$essround.c)),
                                              gei.z.cmc=0,
                                              gei.z.cm=0),
                                      disable.pbkrtest=T,
                                      lmerTest.limit = 400000,infer=c(T,T),adjust="none")
change_in_diff_mod8_GEI
```

```
##  gndr.c essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5        4.5  0.157 0.0705 28.2   0.0123    0.301   2.222  0.0345
##     0.5        4.5  0.488 0.0673 28.7   0.3507    0.626   7.257  <.0001
##    -0.5       -4.5  0.318 0.0710 26.1   0.1719    0.464   4.475  0.0001
##     0.5       -4.5  0.733 0.0726 27.1   0.5838    0.882  10.094  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod8_GEI,adjust="none",infer=c(T,T))
```

```
##  contrast                                                 estimate     SE   df lower.CL upper.CL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5       -0.332 0.0243 21.4  -0.3821  -0.2813
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   -0.161 0.1120 24.4  -0.3926   0.0701
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)    -0.576 0.1090 25.0  -0.8015  -0.3507
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)       0.170 0.1070 25.0  -0.0494   0.3904
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)       -0.244 0.1090 25.6  -0.4688  -0.0199
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)   -0.415 0.0438 24.9  -0.5050  -0.3247
##  t.ratio p.value
##  -13.669  <.0001
##   -1.437  0.1634
##   -5.264  <.0001
##    1.597  0.1229
##   -2.239  0.0341
##   -9.481  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# custom contrast to obtain the estimates for the entire span

changes_in_diff <- list(
  diff_ESS10 = c(-1,1,0,0),
  diff_ESS1 = c(0,0,-1,1))

diff_mod8_GEI<-contrast(change_in_diff_mod8_GEI,method = changes_in_diff,adjust="none",infer=c(T,T))

# gender differences at ESS1 and ESS10
diff_mod8_GEI
```

```
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10    0.332 0.0243 21.4    0.281    0.382  13.669  <.0001
##  diff_ESS1     0.415 0.0438 24.9    0.325    0.505   9.481  <.0001
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
##  diff_ESS10 - diff_ESS1  -0.0831 0.0602 23.4   -0.208   0.0413  -1.380  0.1806
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

## mod7_GGGI and mod8_GGGI: Countries' progress in GGGI and convergence of gender gap

* First fit a model including country average of gender-equality and its interaction with gender and time and differential development for women and men (fixed and random)

* Then fit a model that adds time-specific within-country fluctuations from country's average gender-equality (fixed and random) and compare to a model without these parameters.


``` r
mod7_GGGI<-lmer(FM.z~gndr.c+gggi.z.cm+gndr.c:gggi.z.cm+
                      essround.c+
                      gndr.c:essround.c+
                      (gndr.c+essround.c+gndr.c:essround.c|cntry),
               data=diff_dat,REML=F,weights = pspwght,
                    control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

summary(mod7_GGGI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method [lmerModLmerTest]
## Formula: FM.z ~ gndr.c + gggi.z.cm + gndr.c:gggi.z.cm + essround.c + gndr.c:essround.c +  
##     (gndr.c + essround.c + gndr.c:essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##  124468.1  124616.3  -62217.1  124434.1     44990 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.8437 -0.5941  0.0232  0.5871  6.8334 
## 
## Random effects:
##  Groups   Name              Variance  Std.Dev. Corr             
##  cntry    (Intercept)       0.0676907 0.26017                   
##           gndr.c            0.0080869 0.08993  -0.13            
##           essround.c        0.0016580 0.04072   0.33 -0.22      
##           gndr.c:essround.c 0.0002751 0.01659  -0.70  0.38 -0.38
##  Residual                   0.9617056 0.98067                   
## Number of obs: 45007, groups:  cntry, 33
## 
## Fixed effects:
##                    Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)        0.455429   0.045894 31.642786   9.923 3.09e-11 ***
## gndr.c             0.362649   0.018426 29.457905  19.681  < 2e-16 ***
## gggi.z.cm         -0.122129   0.042655 32.682658  -2.863 0.007267 ** 
## essround.c        -0.027424   0.007803 29.277979  -3.515 0.001454 ** 
## gndr.c:gggi.z.cm   0.072668   0.019240 34.803721   3.777 0.000596 ***
## gndr.c:essround.c -0.007968   0.005268 28.741444  -1.512 0.141340    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c ggg.z. essrn. gn.:..
## gndr.c      -0.106                            
## gggi.z.cm   -0.001 -0.017                     
## essround.c   0.289 -0.183  0.003              
## gndr.c:gg.. -0.004 -0.013 -0.038 -0.009       
## gndr.c:ssr. -0.386  0.124 -0.028 -0.205 -0.034
```

``` r
getFE(mod7_GGGI)
```

```
##                    Est.   SE    df     t     p    LL    UL
## (Intercept)        0.46 0.05 31.64  9.92 0.000  0.36  0.55
## gndr.c             0.36 0.02 29.46 19.68 0.000  0.32  0.40
## gggi.z.cm         -0.12 0.04 32.68 -2.86 0.007 -0.21 -0.04
## essround.c        -0.03 0.01 29.28 -3.51 0.001 -0.04 -0.01
## gndr.c:gggi.z.cm   0.07 0.02 34.80  3.78 0.001  0.03  0.11
## gndr.c:essround.c -0.01 0.01 28.74 -1.51 0.141 -0.02  0.00
```

``` r
getVC(mod7_GGGI)
```

```
##         grp              var1              var2 sdcor vcov
## 1     cntry       (Intercept)              <NA>  0.26 0.07
## 2     cntry            gndr.c              <NA>  0.09 0.01
## 3     cntry        essround.c              <NA>  0.04 0.00
## 4     cntry gndr.c:essround.c              <NA>  0.02 0.00
## 5     cntry       (Intercept)            gndr.c -0.13 0.00
## 6     cntry       (Intercept)        essround.c  0.33 0.00
## 7     cntry       (Intercept) gndr.c:essround.c -0.70 0.00
## 8     cntry            gndr.c        essround.c -0.22 0.00
## 9     cntry            gndr.c gndr.c:essround.c  0.38 0.00
## 10    cntry        essround.c gndr.c:essround.c -0.38 0.00
## 11 Residual              <NA>              <NA>  0.98 0.96
```

``` r
anova(mod2_GGGI,mod7_GGGI)
```

```
## Data: diff_dat
## Models:
## mod2_GGGI: FM.z ~ gndr.c + gggi.z.cm + gndr.c:gggi.z.cm + (gndr.c | cntry)
## mod7_GGGI: FM.z ~ gndr.c + gggi.z.cm + gndr.c:gggi.z.cm + essround.c + gndr.c:essround.c + (gndr.c + essround.c + gndr.c:essround.c | cntry)
##           npar    AIC    BIC logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod2_GGGI    8 124897 124967 -62440    124881                         
## mod7_GGGI   17 124468 124616 -62217    124434 446.65  9  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

### Examine trends


``` r
change_mod7_GGGI<-emmeans(mod7_GGGI,specs="essround.c",by="gndr.c",
                              at=list(gndr.c=c(-0.5,0.5),
                                      essround.c=rev(range(diff_dat$essround.c)),
                                      gggi.z.cmc=0,
                                      gggi.z.cm=0),
                              disable.pbkrtest=T,
                              lmerTest.limit = 400000,infer=c(T,T),adjust="none")
change_mod7_GGGI
```

```
## gndr.c = -0.5:
##  essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.169 0.0729 26.1   0.0187    0.319   2.312  0.0289
##        -4.5  0.380 0.0483 21.6   0.2794    0.480   7.865  <.0001
## 
## gndr.c =  0.5:
##  essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.495 0.0610 24.0   0.3695    0.621   8.120  <.0001
##        -4.5  0.778 0.0537 23.4   0.6671    0.889  14.483  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_mod7_GGGI,adjust="none",infer=c(T,T))
```

```
## gndr.c = -0.5:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   -0.211 0.0786 28.7   -0.372  -0.0502  -2.684  0.0119
## 
## gndr.c =  0.5:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   -0.283 0.0694 29.1   -0.425  -0.1408  -4.075  0.0003
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# change in gender differences over time

change_in_diff_mod7_GGGI<-emmeans(mod7_GGGI,specs=c("gndr.c","essround.c"),
                                      at=list(gndr.c=c(-0.5,0.5),
                                              essround.c=rev(range(diff_dat$essround.c)),
                                              gggi.z.cmc=0,
                                              gggi.z.cm=0),
                                      disable.pbkrtest=T,
                                      lmerTest.limit = 400000,infer=c(T,T),adjust="none")
change_in_diff_mod7_GGGI
```

```
##  gndr.c essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5        4.5  0.169 0.0729 26.1   0.0187    0.319   2.312  0.0289
##     0.5        4.5  0.495 0.0610 24.0   0.3695    0.621   8.120  <.0001
##    -0.5       -4.5  0.380 0.0483 21.6   0.2794    0.480   7.865  <.0001
##     0.5       -4.5  0.778 0.0537 23.4   0.6671    0.889  14.483  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod7_GGGI,adjust="none",infer=c(T,T))
```

```
##  contrast                                                 estimate     SE   df lower.CL upper.CL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5       -0.327 0.0318 27.7   -0.392  -0.2617
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   -0.211 0.0786 28.7   -0.372  -0.0502
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)    -0.609 0.0758 28.3   -0.765  -0.4543
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)       0.116 0.0693 28.1   -0.026   0.2577
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)       -0.283 0.0694 29.1   -0.425  -0.1408
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)   -0.399 0.0282 28.2   -0.456  -0.3408
##  t.ratio p.value
##  -10.285  <.0001
##   -2.684  0.0119
##   -8.041  <.0001
##    1.672  0.1056
##   -4.075  0.0003
##  -14.147  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# custom contrast to obtain the estimates for the entire span

changes_in_diff <- list(
  diff_ESS10 = c(-1,1,0,0),
  diff_ESS1 = c(0,0,-1,1))

diff_mod7_GGGI<-contrast(change_in_diff_mod7_GGGI,method = changes_in_diff,adjust="none",infer=c(T,T))

# gender differences at ESS1 and ESS10
diff_mod7_GGGI
```

```
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10    0.327 0.0318 27.7    0.262    0.392  10.285  <.0001
##  diff_ESS1     0.399 0.0282 28.2    0.341    0.456  14.147  <.0001
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
##  diff_ESS10 - diff_ESS1  -0.0717 0.0474 28.7   -0.169   0.0253  -1.512  0.1413
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```


``` r
mod8_GGGI<-lmer(FM.z~gndr.c+gggi.z.cm+gndr.c:gggi.z.cm+
                      essround.c+
                      gndr.c:essround.c+
                      gggi.z.cmc+gggi.z.cmc:gndr.c+
                      (gndr.c+essround.c+gndr.c:essround.c+gggi.z.cmc+gggi.z.cmc:gndr.c|cntry),
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
## (Intercept)        0.45 0.05 32.11  9.74 0.000  0.35  0.54
## gndr.c             0.36 0.02 27.97 19.11 0.000  0.32  0.40
## gggi.z.cm         -0.12 0.04 31.21 -2.73 0.010 -0.20 -0.03
## essround.c        -0.02 0.01 28.73 -1.93 0.064 -0.04  0.00
## gggi.z.cmc        -0.09 0.05 19.22 -1.88 0.076 -0.19  0.01
## gndr.c:gggi.z.cm   0.06 0.02 33.47  3.36 0.002  0.03  0.10
## gndr.c:essround.c -0.01 0.01 34.17 -0.99 0.330 -0.02  0.01
## gndr.c:gggi.z.cmc  0.01 0.04 19.11  0.24 0.812 -0.07  0.09
```

``` r
getVC(mod8_GGGI)
```

```
##         grp              var1              var2 sdcor  vcov
## 1     cntry       (Intercept)              <NA>  0.26  0.07
## 2     cntry            gndr.c              <NA>  0.09  0.01
## 3     cntry        essround.c              <NA>  0.05  0.00
## 4     cntry        gggi.z.cmc              <NA>  0.21  0.05
## 5     cntry gndr.c:essround.c              <NA>  0.02  0.00
## 6     cntry gndr.c:gggi.z.cmc              <NA>  0.11  0.01
## 7     cntry       (Intercept)            gndr.c -0.12  0.00
## 8     cntry       (Intercept)        essround.c  0.21  0.00
## 9     cntry       (Intercept)        gggi.z.cmc -0.01  0.00
## 10    cntry       (Intercept) gndr.c:essround.c -0.72  0.00
## 11    cntry       (Intercept) gndr.c:gggi.z.cmc  0.09  0.00
## 12    cntry            gndr.c        essround.c -0.32  0.00
## 13    cntry            gndr.c        gggi.z.cmc  0.26  0.00
## 14    cntry            gndr.c gndr.c:essround.c  0.44  0.00
## 15    cntry            gndr.c gndr.c:gggi.z.cmc -0.08  0.00
## 16    cntry        essround.c        gggi.z.cmc -0.77 -0.01
## 17    cntry        essround.c gndr.c:essround.c -0.08  0.00
## 18    cntry        essround.c gndr.c:gggi.z.cmc  0.00  0.00
## 19    cntry        gggi.z.cmc gndr.c:essround.c  0.24  0.00
## 20    cntry        gggi.z.cmc gndr.c:gggi.z.cmc -0.60 -0.01
## 21    cntry gndr.c:essround.c gndr.c:gggi.z.cmc -0.59  0.00
## 22 Residual              <NA>              <NA>  0.98  0.96
```

``` r
anova(mod2_GGGI,mod7_GGGI,mod8_GGGI)
```

```
## Data: diff_dat
## Models:
## mod2_GGGI: FM.z ~ gndr.c + gggi.z.cm + gndr.c:gggi.z.cm + (gndr.c | cntry)
## mod7_GGGI: FM.z ~ gndr.c + gggi.z.cm + gndr.c:gggi.z.cm + essround.c + gndr.c:essround.c + (gndr.c + essround.c + gndr.c:essround.c | cntry)
## mod8_GGGI: FM.z ~ gndr.c + gggi.z.cm + gndr.c:gggi.z.cm + essround.c + gndr.c:essround.c + gggi.z.cmc + gggi.z.cmc:gndr.c + (gndr.c + essround.c + gndr.c:essround.c + gggi.z.cmc + gggi.z.cmc:gndr.c | cntry)
##           npar    AIC    BIC logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod2_GGGI    8 124897 124967 -62440    124881                         
## mod7_GGGI   17 124468 124616 -62217    124434 446.65  9  < 2.2e-16 ***
## mod8_GGGI   30 124452 124714 -62196    124392  41.69 13   7.36e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

### Examine trends


``` r
# gender specific change over time

change_mod8_GGGI<-emmeans(mod8_GGGI,specs="essround.c",by="gndr.c",
                              at=list(gndr.c=c(-0.5,0.5),
                                      essround.c=rev(range(diff_dat$essround.c)),
                                      gggi.z.cmc=0,
                                      gggi.z.cm=0),
                              disable.pbkrtest=T,
                              lmerTest.limit = 400000,infer=c(T,T),adjust="none")
change_mod8_GGGI
```

```
## gndr.c = -0.5:
##  essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.190 0.0783 24.9   0.0293    0.352   2.434  0.0224
##        -4.5  0.342 0.0592 17.0   0.2175    0.467   5.785  <.0001
## 
## gndr.c =  0.5:
##  essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.522 0.0671 23.2   0.3833    0.661   7.780  <.0001
##        -4.5  0.734 0.0671 18.9   0.5936    0.875  10.940  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_mod8_GGGI,adjust="none",infer=c(T,T))
```

```
## gndr.c = -0.5:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   -0.152 0.1010 29.5   -0.357   0.0536  -1.510  0.1416
## 
## gndr.c =  0.5:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   -0.212 0.0978 26.5   -0.413  -0.0112  -2.169  0.0393
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# change in gender differences over time

change_in_diff_mod8_GGGI<-emmeans(mod8_GGGI,specs=c("gndr.c","essround.c"),
                                      at=list(gndr.c=c(-0.5,0.5),
                                              essround.c=rev(range(diff_dat$essround.c)),
                                              gggi.z.cmc=0,
                                              gggi.z.cm=0),
                                      disable.pbkrtest=T,
                                      lmerTest.limit = 400000,infer=c(T,T),adjust="none")
change_in_diff_mod8_GGGI
```

```
##  gndr.c essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5        4.5  0.190 0.0783 24.9   0.0293    0.352   2.434  0.0224
##     0.5        4.5  0.522 0.0671 23.2   0.3833    0.661   7.780  <.0001
##    -0.5       -4.5  0.342 0.0592 17.0   0.2175    0.467   5.785  <.0001
##     0.5       -4.5  0.734 0.0671 18.9   0.5936    0.875  10.940  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod8_GGGI,adjust="none",infer=c(T,T))
```

```
##  contrast                                                 estimate     SE   df lower.CL upper.CL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5       -0.332 0.0358 28.2  -0.4048  -0.2584
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   -0.152 0.1010 29.5  -0.3574   0.0536
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)    -0.544 0.1010 27.9  -0.7498  -0.3374
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)       0.180 0.0917 27.0  -0.0084   0.3678
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)       -0.212 0.0978 26.5  -0.4128  -0.0112
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)   -0.392 0.0359 26.9  -0.4655  -0.3180
##  t.ratio p.value
##   -9.273  <.0001
##   -1.510  0.1416
##   -5.401  <.0001
##    1.960  0.0604
##   -2.169  0.0393
##  -10.899  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# custom contrast to obtain the estimates for the entire span

changes_in_diff <- list(
  diff_ESS10 = c(-1,1,0,0),
  diff_ESS1 = c(0,0,-1,1))

diff_mod8_GGGI<-contrast(change_in_diff_mod8_GGGI,method = changes_in_diff,adjust="none",infer=c(T,T))

# gender differences at ESS1 and ESS10
diff_mod8_GGGI
```

```
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10    0.332 0.0358 28.2    0.258    0.405   9.273  <.0001
##  diff_ESS1     0.392 0.0359 26.9    0.318    0.465  10.899  <.0001
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
##  diff_ESS10 - diff_ESS1  -0.0601 0.0609 34.2   -0.184   0.0636  -0.987  0.3304
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

## mod7_GDI and mod8_GDI: Countries' progress in GDI and convergence of gender gap

* First fit a model including country average of gender-equality and its interaction with gender and time and differential development for women and men (fixed and random)

* Then fit a model that adds time-specific within-country fluctuations from country's average gender-equality (fixed and random) and compare to a model without these parameters.


``` r
mod7_GDI<-lmer(FM.z~gndr.c+gdi.z.cm+gndr.c:gdi.z.cm+
                      essround.c+
                      gndr.c:essround.c+
                      (gndr.c+essround.c+gndr.c:essround.c|cntry),
               data=diff_dat,REML=F,weights = pspwght,
                    control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

summary(mod7_GDI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method [lmerModLmerTest]
## Formula: FM.z ~ gndr.c + gdi.z.cm + gndr.c:gdi.z.cm + essround.c + gndr.c:essround.c +  
##     (gndr.c + essround.c + gndr.c:essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##  175212.6  175366.4  -87589.3  175178.6     63030 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.8308 -0.5885  0.0217  0.5920  9.3410 
## 
## Random effects:
##  Groups   Name              Variance  Std.Dev. Corr             
##  cntry    (Intercept)       0.0880294 0.29670                   
##           gndr.c            0.0145897 0.12079  -0.47            
##           essround.c        0.0025337 0.05034   0.67 -0.64      
##           gndr.c:essround.c 0.0001708 0.01307  -0.60  0.13 -0.39
##  Residual                   0.9731588 0.98649                   
## Number of obs: 63047, groups:  cntry, 33
## 
## Fixed effects:
##                    Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)        0.435068   0.051979 28.937993   8.370 3.22e-09 ***
## gndr.c             0.368478   0.022631 28.587167  16.282 5.37e-16 ***
## gdi.z.cm           0.139171   0.039382 30.912860   3.534  0.00131 ** 
## essround.c        -0.021074   0.008999 26.219772  -2.342  0.02705 *  
## gndr.c:gdi.z.cm   -0.005425   0.019621 36.186499  -0.276  0.78374    
## gndr.c:essround.c -0.011602   0.003706 26.737528  -3.131  0.00419 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c gd.z.c essrn. gn.:..
## gndr.c      -0.428                            
## gdi.z.cm    -0.005 -0.008                     
## essround.c   0.642 -0.586  0.012              
## gndr.c:gd.. -0.003 -0.011 -0.084 -0.009       
## gndr.c:ssr. -0.369  0.066  0.015 -0.233 -0.009
```

``` r
getFE(mod7_GDI)
```

```
##                    Est.   SE    df     t     p    LL   UL
## (Intercept)        0.44 0.05 28.94  8.37 0.000  0.33 0.54
## gndr.c             0.37 0.02 28.59 16.28 0.000  0.32 0.41
## gdi.z.cm           0.14 0.04 30.91  3.53 0.001  0.06 0.22
## essround.c        -0.02 0.01 26.22 -2.34 0.027 -0.04 0.00
## gndr.c:gdi.z.cm   -0.01 0.02 36.19 -0.28 0.784 -0.05 0.03
## gndr.c:essround.c -0.01 0.00 26.74 -3.13 0.004 -0.02 0.00
```

``` r
getVC(mod7_GDI)
```

```
##         grp              var1              var2 sdcor  vcov
## 1     cntry       (Intercept)              <NA>  0.30  0.09
## 2     cntry            gndr.c              <NA>  0.12  0.01
## 3     cntry        essround.c              <NA>  0.05  0.00
## 4     cntry gndr.c:essround.c              <NA>  0.01  0.00
## 5     cntry       (Intercept)            gndr.c -0.47 -0.02
## 6     cntry       (Intercept)        essround.c  0.67  0.01
## 7     cntry       (Intercept) gndr.c:essround.c -0.60  0.00
## 8     cntry            gndr.c        essround.c -0.64  0.00
## 9     cntry            gndr.c gndr.c:essround.c  0.13  0.00
## 10    cntry        essround.c gndr.c:essround.c -0.39  0.00
## 11 Residual              <NA>              <NA>  0.99  0.97
```

``` r
anova(mod2_GDI,mod7_GDI)
```

```
## Data: diff_dat
## Models:
## mod2_GDI: FM.z ~ gndr.c + gdi.z.cm + gndr.c:gdi.z.cm + (gndr.c | cntry)
## mod7_GDI: FM.z ~ gndr.c + gdi.z.cm + gndr.c:gdi.z.cm + essround.c + gndr.c:essround.c + (gndr.c + essround.c + gndr.c:essround.c | cntry)
##          npar    AIC    BIC logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod2_GDI    8 176203 176275 -88093    176187                         
## mod7_GDI   17 175213 175366 -87589    175179 1008.4  9  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

### Examine trends


``` r
change_mod7_GDI<-emmeans(mod7_GDI,specs="essround.c",by="gndr.c",
                              at=list(gndr.c=c(-0.5,0.5),
                                      essround.c=rev(range(diff_dat$essround.c)),
                                      gdi.z.cmc=0,
                                      gdi.z.cm=0),
                              disable.pbkrtest=T,
                              lmerTest.limit = 400000,infer=c(T,T),adjust="none")
change_mod7_GDI
```

```
## gndr.c = -0.5:
##  essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.182 0.0937 26.9  -0.0101    0.374   1.944  0.0624
##        -4.5  0.320 0.0403 27.1   0.2368    0.402   7.921  <.0001
## 
## gndr.c =  0.5:
##  essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.498 0.0758 26.5   0.3428    0.654   6.579  <.0001
##        -4.5  0.740 0.0449 28.2   0.6482    0.832  16.475  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_mod7_GDI,adjust="none",infer=c(T,T))
```

```
## gndr.c = -0.5:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   -0.137 0.0864 25.4   -0.315   0.0404  -1.591  0.1240
## 
## gndr.c =  0.5:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   -0.242 0.0788 25.7   -0.404  -0.0798  -3.070  0.0050
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# change in gender differences over time

change_in_diff_mod7_GDI<-emmeans(mod7_GDI,specs=c("gndr.c","essround.c"),
                                      at=list(gndr.c=c(-0.5,0.5),
                                              essround.c=rev(range(diff_dat$essround.c)),
                                              gdi.z.cmc=0,
                                              gdi.z.cm=0),
                                      disable.pbkrtest=T,
                                      lmerTest.limit = 400000,infer=c(T,T),adjust="none")
change_in_diff_mod7_GDI
```

```
##  gndr.c essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5        4.5  0.182 0.0937 26.9  -0.0101    0.374   1.944  0.0624
##     0.5        4.5  0.498 0.0758 26.5   0.3428    0.654   6.579  <.0001
##    -0.5       -4.5  0.320 0.0403 27.1   0.2368    0.402   7.921  <.0001
##     0.5       -4.5  0.740 0.0449 28.2   0.6482    0.832  16.475  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod7_GDI,adjust="none",infer=c(T,T))
```

```
##  contrast                                                 estimate     SE   df lower.CL upper.CL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5       -0.316 0.0290 21.1  -0.3765  -0.2560
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   -0.137 0.0864 25.4  -0.3153   0.0404
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)    -0.558 0.0960 26.9  -0.7552  -0.3611
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)       0.179 0.0702 25.1   0.0343   0.3233
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)       -0.242 0.0788 25.7  -0.4039  -0.0798
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)   -0.421 0.0272 31.2  -0.4762  -0.3652
##  t.ratio p.value
##  -10.910  <.0001
##   -1.591  0.1240
##   -5.813  <.0001
##    2.548  0.0173
##   -3.070  0.0050
##  -15.462  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# custom contrast to obtain the estimates for the entire span

changes_in_diff <- list(
  diff_ESS10 = c(-1,1,0,0),
  diff_ESS1 = c(0,0,-1,1))

diff_mod7_GDI<-contrast(change_in_diff_mod7_GDI,method = changes_in_diff,adjust="none",infer=c(T,T))

# gender differences at ESS1 and ESS10
diff_mod7_GDI
```

```
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10    0.316 0.0290 21.1    0.256    0.377  10.910  <.0001
##  diff_ESS1     0.421 0.0272 31.2    0.365    0.476  15.462  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# Test for whether and how much the differences have changed
pairs(diff_mod7_GDI,infer=c(T,T))
```

```
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10 - diff_ESS1   -0.104 0.0334 26.7   -0.173  -0.0359  -3.131  0.0042
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```


``` r
mod8_GDI<-lmer(FM.z~gndr.c+gdi.z.cm+gndr.c:gdi.z.cm+
                      essround.c+
                      gndr.c:essround.c+
                      gdi.z.cmc+gdi.z.cmc:gndr.c+
                      (gndr.c+essround.c+gndr.c:essround.c+gdi.z.cmc+gdi.z.cmc:gndr.c|cntry),
               data=diff_dat,REML=F,weights = pspwght,
                    control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

getFE(mod8_GDI)
```

```
##                    Est.   SE    df     t     p    LL    UL
## (Intercept)        0.43 0.05 31.78  8.49 0.000  0.33  0.54
## gndr.c             0.38 0.02 28.87 17.46 0.000  0.33  0.42
## gdi.z.cm           0.02 0.04 23.59  0.43 0.671 -0.07  0.10
## essround.c        -0.02 0.01 25.79 -2.85 0.008 -0.03 -0.01
## gdi.z.cmc          0.05 0.06 21.62  0.91 0.373 -0.07  0.18
## gndr.c:gdi.z.cm    0.01 0.02 31.80  0.33 0.741 -0.03  0.05
## gndr.c:essround.c -0.01 0.00 11.41 -3.68 0.003 -0.02 -0.01
## gndr.c:gdi.z.cmc  -0.04 0.05  9.26 -0.87 0.405 -0.16  0.07
```

``` r
getVC(mod8_GDI)
```

```
##         grp              var1              var2 sdcor  vcov
## 1     cntry       (Intercept)              <NA>  0.29  0.08
## 2     cntry            gndr.c              <NA>  0.11  0.01
## 3     cntry        essround.c              <NA>  0.04  0.00
## 4     cntry         gdi.z.cmc              <NA>  0.29  0.09
## 5     cntry gndr.c:essround.c              <NA>  0.01  0.00
## 6     cntry  gndr.c:gdi.z.cmc              <NA>  0.20  0.04
## 7     cntry       (Intercept)            gndr.c -0.40 -0.01
## 8     cntry       (Intercept)        essround.c  0.65  0.01
## 9     cntry       (Intercept)         gdi.z.cmc  0.38  0.03
## 10    cntry       (Intercept) gndr.c:essround.c -0.57  0.00
## 11    cntry       (Intercept)  gndr.c:gdi.z.cmc -0.10 -0.01
## 12    cntry            gndr.c        essround.c -0.76  0.00
## 13    cntry            gndr.c         gdi.z.cmc  0.12  0.00
## 14    cntry            gndr.c gndr.c:essround.c  0.19  0.00
## 15    cntry            gndr.c  gndr.c:gdi.z.cmc -0.26 -0.01
## 16    cntry        essround.c         gdi.z.cmc  0.04  0.00
## 17    cntry        essround.c gndr.c:essround.c -0.09  0.00
## 18    cntry        essround.c  gndr.c:gdi.z.cmc  0.14  0.00
## 19    cntry         gdi.z.cmc gndr.c:essround.c -0.11  0.00
## 20    cntry         gdi.z.cmc  gndr.c:gdi.z.cmc -0.86 -0.05
## 21    cntry gndr.c:essround.c  gndr.c:gdi.z.cmc -0.19  0.00
## 22 Residual              <NA>              <NA>  0.99  0.97
```

``` r
anova(mod2_GDI,mod7_GDI,mod8_GDI)
```

```
## Data: diff_dat
## Models:
## mod2_GDI: FM.z ~ gndr.c + gdi.z.cm + gndr.c:gdi.z.cm + (gndr.c | cntry)
## mod7_GDI: FM.z ~ gndr.c + gdi.z.cm + gndr.c:gdi.z.cm + essround.c + gndr.c:essround.c + (gndr.c + essround.c + gndr.c:essround.c | cntry)
## mod8_GDI: FM.z ~ gndr.c + gdi.z.cm + gndr.c:gdi.z.cm + essround.c + gndr.c:essround.c + gdi.z.cmc + gdi.z.cmc:gndr.c + (gndr.c + essround.c + gndr.c:essround.c + gdi.z.cmc + gdi.z.cmc:gndr.c | cntry)
##          npar    AIC    BIC logLik -2*log(L)    Chisq Df Pr(>Chisq)    
## mod2_GDI    8 176203 176275 -88093    176187                           
## mod7_GDI   17 175213 175366 -87589    175179 1008.412  9  < 2.2e-16 ***
## mod8_GDI   30 175176 175447 -87558    175116   62.656 13  1.751e-08 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

### Examine trends


``` r
# gender specific change over time

change_mod8_GDI<-emmeans(mod8_GDI,specs="essround.c",by="gndr.c",
                              at=list(gndr.c=c(-0.5,0.5),
                                      essround.c=rev(range(diff_dat$essround.c)),
                                      gdi.z.cmc=0,
                                      gdi.z.cm=0),
                              disable.pbkrtest=T,
                              lmerTest.limit = 400000,infer=c(T,T),adjust="none")
change_mod8_GDI
```

```
## gndr.c = -0.5:
##  essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.186 0.0805 29.7   0.0213    0.350   2.308  0.0281
##        -4.5  0.301 0.0414 21.4   0.2154    0.387   7.279  <.0001
## 
## gndr.c =  0.5:
##  essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.505 0.0665 30.0   0.3693    0.641   7.593  <.0001
##        -4.5  0.734 0.0453 24.7   0.6403    0.827  16.194  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_mod8_GDI,adjust="none",infer=c(T,T))
```

```
## gndr.c = -0.5:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   -0.116 0.0631 24.2   -0.246   0.0145  -1.833  0.0792
## 
## gndr.c =  0.5:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   -0.228 0.0613 26.1   -0.354  -0.1024  -3.726  0.0009
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# change in gender differences over time

change_in_diff_mod8_GDI<-emmeans(mod8_GDI,specs=c("gndr.c","essround.c"),
                                      at=list(gndr.c=c(-0.5,0.5),
                                              essround.c=rev(range(diff_dat$essround.c)),
                                              gdi.z.cmc=0,
                                              gdi.z.cm=0),
                                      disable.pbkrtest=T,
                                      lmerTest.limit = 400000,infer=c(T,T),adjust="none")
change_in_diff_mod8_GDI
```

```
##  gndr.c essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5        4.5  0.186 0.0805 29.7   0.0213    0.350   2.308  0.0281
##     0.5        4.5  0.505 0.0665 30.0   0.3693    0.641   7.593  <.0001
##    -0.5       -4.5  0.301 0.0414 21.4   0.2154    0.387   7.279  <.0001
##     0.5       -4.5  0.734 0.0453 24.7   0.6403    0.827  16.194  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod8_GDI,adjust="none",infer=c(T,T))
```

```
##  contrast                                                 estimate     SE   df lower.CL upper.CL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5       -0.319 0.0270 20.2   -0.376  -0.2632
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   -0.116 0.0631 24.2   -0.246   0.0145
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)    -0.548 0.0762 27.0   -0.704  -0.3914
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)       0.204 0.0488 24.6    0.103   0.3045
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)       -0.228 0.0613 26.1   -0.354  -0.1024
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)   -0.432 0.0259 23.1   -0.486  -0.3787
##  t.ratio p.value
##  -11.849  <.0001
##   -1.833  0.0792
##   -7.189  <.0001
##    4.175  0.0003
##   -3.726  0.0009
##  -16.697  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# custom contrast to obtain the estimates for the entire span

changes_in_diff <- list(
  diff_ESS10 = c(-1,1,0,0),
  diff_ESS1 = c(0,0,-1,1))

diff_mod8_GDI<-contrast(change_in_diff_mod8_GDI,method = changes_in_diff,adjust="none",infer=c(T,T))

# gender differences at ESS1 and ESS10
diff_mod8_GDI
```

```
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10    0.319 0.0270 20.2    0.263    0.376  11.849  <.0001
##  diff_ESS1     0.432 0.0259 23.1    0.379    0.486  16.697  <.0001
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
##  diff_ESS10 - diff_ESS1   -0.113 0.0307 11.4    -0.18  -0.0456  -3.678  0.0034
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

## mod7_log_GDP and mod8_log_GDP: Countries' progress in log_GDP and convergence of gender gap

* First fit a model including country average of gender-equality and its interaction with gender and time and differential development for women and men (fixed and random)

* Then fit a model that adds time-specific within-country fluctuations from country's average gender-equality (fixed and random) and compare to a model without these parameters.


``` r
mod7_log_GDP<-lmer(FM.z~gndr.c+log_gdp.z.cm+gndr.c:log_gdp.z.cm+
                      essround.c+
                      gndr.c:essround.c+
                      (gndr.c+essround.c+gndr.c:essround.c|cntry),
               data=diff_dat,REML=F,weights = pspwght,
                    control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

summary(mod7_log_GDP)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method [lmerModLmerTest]
## Formula: FM.z ~ gndr.c + log_gdp.z.cm + gndr.c:log_gdp.z.cm + essround.c +  
##     gndr.c:essround.c + (gndr.c + essround.c + gndr.c:essround.c |      cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##  175210.7  175364.6  -87588.4  175176.7     63030 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.8479 -0.5882  0.0217  0.5919  9.3620 
## 
## Random effects:
##  Groups   Name              Variance Std.Dev. Corr             
##  cntry    (Intercept)       0.061145 0.24727                   
##           gndr.c            0.009165 0.09573  -0.09            
##           essround.c        0.002524 0.05024   0.27 -0.43      
##           gndr.c:essround.c 0.000165 0.01285  -0.56  0.06 -0.28
##  Residual                   0.973136 0.98648                   
## Number of obs: 63047, groups:  cntry, 33
## 
## Fixed effects:
##                      Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)          0.432612   0.043530 31.703815   9.938 2.92e-11 ***
## gndr.c               0.367975   0.018619 27.571269  19.763  < 2e-16 ***
## log_gdp.z.cm        -0.130229   0.040983 31.490248  -3.178  0.00332 ** 
## essround.c          -0.022167   0.009016 26.507068  -2.459  0.02078 *  
## gndr.c:log_gdp.z.cm  0.052904   0.017738 34.382666   2.983  0.00523 ** 
## gndr.c:essround.c   -0.011275   0.003681 28.809935  -3.063  0.00472 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c lg_g.. essrn. g.:_..
## gndr.c      -0.077                            
## lg_gdp.z.cm  0.001 -0.009                     
## essround.c   0.254 -0.379  0.019              
## gndr.c:l_..  0.000 -0.036  0.022 -0.011       
## gndr.c:ssr. -0.345  0.037 -0.066 -0.166 -0.019
```

``` r
getFE(mod7_log_GDP)
```

```
##                      Est.   SE    df     t     p    LL    UL
## (Intercept)          0.43 0.04 31.70  9.94 0.000  0.34  0.52
## gndr.c               0.37 0.02 27.57 19.76 0.000  0.33  0.41
## log_gdp.z.cm        -0.13 0.04 31.49 -3.18 0.003 -0.21 -0.05
## essround.c          -0.02 0.01 26.51 -2.46 0.021 -0.04  0.00
## gndr.c:log_gdp.z.cm  0.05 0.02 34.38  2.98 0.005  0.02  0.09
## gndr.c:essround.c   -0.01 0.00 28.81 -3.06 0.005 -0.02  0.00
```

``` r
getVC(mod7_log_GDP)
```

```
##         grp              var1              var2 sdcor vcov
## 1     cntry       (Intercept)              <NA>  0.25 0.06
## 2     cntry            gndr.c              <NA>  0.10 0.01
## 3     cntry        essround.c              <NA>  0.05 0.00
## 4     cntry gndr.c:essround.c              <NA>  0.01 0.00
## 5     cntry       (Intercept)            gndr.c -0.09 0.00
## 6     cntry       (Intercept)        essround.c  0.27 0.00
## 7     cntry       (Intercept) gndr.c:essround.c -0.56 0.00
## 8     cntry            gndr.c        essround.c -0.43 0.00
## 9     cntry            gndr.c gndr.c:essround.c  0.06 0.00
## 10    cntry        essround.c gndr.c:essround.c -0.28 0.00
## 11 Residual              <NA>              <NA>  0.99 0.97
```

``` r
anova(mod2_log_GDP,mod7_log_GDP)
```

```
## Data: diff_dat
## Models:
## mod2_log_GDP: FM.z ~ gndr.c + log_gdp.z.cm + gndr.c:log_gdp.z.cm + (gndr.c | cntry)
## mod7_log_GDP: FM.z ~ gndr.c + log_gdp.z.cm + gndr.c:log_gdp.z.cm + essround.c + gndr.c:essround.c + (gndr.c + essround.c + gndr.c:essround.c | cntry)
##              npar    AIC    BIC logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod2_log_GDP    8 176188 176261 -88086    176172                         
## mod7_log_GDP   17 175211 175365 -87588    175177 995.73  9  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

### Examine trends


``` r
change_mod7_log_GDP<-emmeans(mod7_log_GDP,specs="essround.c",by="gndr.c",
                              at=list(gndr.c=c(-0.5,0.5),
                                      essround.c=rev(range(diff_dat$essround.c)),
                                      log_gdp.z.cmc=0,
                                      log_gdp.z.cm=0),
                              disable.pbkrtest=T,
                              lmerTest.limit = 400000,infer=c(T,T),adjust="none")
change_mod7_log_GDP
```

```
## gndr.c = -0.5:
##  essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.174 0.0729 25.7   0.0244    0.324   2.392  0.0244
##        -4.5  0.323 0.0493 22.5   0.2209    0.425   6.551  <.0001
## 
## gndr.c =  0.5:
##  essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.491 0.0624 24.8   0.3630    0.620   7.879  <.0001
##        -4.5  0.742 0.0562 24.3   0.6259    0.858  13.207  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_mod7_log_GDP,adjust="none",infer=c(T,T))
```

```
## gndr.c = -0.5:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   -0.149 0.0855 25.7   -0.325   0.0270  -1.741  0.0937
## 
## gndr.c =  0.5:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   -0.250 0.0801 25.9   -0.415  -0.0856  -3.125  0.0044
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# change in gender differences over time

change_in_diff_mod7_log_GDP<-emmeans(mod7_log_GDP,specs=c("gndr.c","essround.c"),
                                      at=list(gndr.c=c(-0.5,0.5),
                                              essround.c=rev(range(diff_dat$essround.c)),
                                              log_gdp.z.cmc=0,
                                              log_gdp.z.cm=0),
                                      disable.pbkrtest=T,
                                      lmerTest.limit = 400000,infer=c(T,T),adjust="none")
change_in_diff_mod7_log_GDP
```

```
##  gndr.c essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5        4.5  0.174 0.0729 25.7   0.0244    0.324   2.392  0.0244
##     0.5        4.5  0.491 0.0624 24.8   0.3630    0.620   7.879  <.0001
##    -0.5       -4.5  0.323 0.0493 22.5   0.2209    0.425   6.551  <.0001
##     0.5       -4.5  0.742 0.0562 24.3   0.6259    0.858  13.207  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod7_log_GDP,adjust="none",infer=c(T,T))
```

```
##  contrast                                                 estimate     SE   df lower.CL upper.CL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5       -0.317 0.0254 23.9  -0.3696  -0.2649
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   -0.149 0.0855 25.7  -0.3246   0.0270
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)    -0.567 0.0899 26.0  -0.7522  -0.3828
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)       0.168 0.0761 25.3   0.0119   0.3250
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)       -0.250 0.0801 25.9  -0.4149  -0.0856
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)   -0.419 0.0245 29.2  -0.4687  -0.3687
##  t.ratio p.value
##  -12.503  <.0001
##   -1.741  0.0937
##   -6.315  <.0001
##    2.215  0.0360
##   -3.125  0.0044
##  -17.116  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# custom contrast to obtain the estimates for the entire span

changes_in_diff <- list(
  diff_ESS10 = c(-1,1,0,0),
  diff_ESS1 = c(0,0,-1,1))

diff_mod7_log_GDP<-contrast(change_in_diff_mod7_log_GDP,method = changes_in_diff,adjust="none",infer=c(T,T))

# gender differences at ESS1 and ESS10
diff_mod7_log_GDP
```

```
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10    0.317 0.0254 23.9    0.265    0.370  12.503  <.0001
##  diff_ESS1     0.419 0.0245 29.2    0.369    0.469  17.116  <.0001
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
##  diff_ESS10 - diff_ESS1   -0.101 0.0331 28.8   -0.169  -0.0337  -3.063  0.0047
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```


``` r
mod8_log_GDP<-lmer(FM.z~gndr.c+log_gdp.z.cm+gndr.c:log_gdp.z.cm+
                      essround.c+
                      gndr.c:essround.c+
                      log_gdp.z.cmc+log_gdp.z.cmc:gndr.c+
                      (gndr.c+essround.c+gndr.c:essround.c+log_gdp.z.cmc+log_gdp.z.cmc:gndr.c|cntry),
               data=diff_dat,REML=F,weights = pspwght,
                    control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## Warning: Model failed to converge with 2 negative eigenvalues: -3.1e+03 -1.5e+04
```

``` r
getFE(mod8_log_GDP)
```

```
##                       Est.   SE    df     t     p    LL    UL
## (Intercept)           0.44 0.04 32.85 10.18 0.000  0.35  0.53
## gndr.c                0.37 0.02 26.36 20.17 0.000  0.33  0.41
## log_gdp.z.cm         -0.17 0.04 33.49 -3.88 0.000 -0.25 -0.08
## essround.c           -0.03 0.01 21.30 -3.72 0.001 -0.05 -0.01
## log_gdp.z.cmc         0.07 0.07 17.28  0.95 0.357 -0.08  0.22
## gndr.c:log_gdp.z.cm   0.04 0.02 26.19  2.69 0.012  0.01  0.08
## gndr.c:essround.c    -0.01 0.00 24.64 -1.13 0.270 -0.01  0.00
## gndr.c:log_gdp.z.cmc -0.06 0.05 24.62 -1.15 0.261 -0.16  0.05
```

``` r
getVC(mod8_log_GDP)
```

```
##         grp                 var1                 var2 sdcor  vcov
## 1     cntry          (Intercept)                 <NA>  0.25  0.06
## 2     cntry               gndr.c                 <NA>  0.09  0.01
## 3     cntry           essround.c                 <NA>  0.05  0.00
## 4     cntry        log_gdp.z.cmc                 <NA>  0.32  0.11
## 5     cntry    gndr.c:essround.c                 <NA>  0.01  0.00
## 6     cntry gndr.c:log_gdp.z.cmc                 <NA>  0.14  0.02
## 7     cntry          (Intercept)               gndr.c -0.13  0.00
## 8     cntry          (Intercept)           essround.c  0.14  0.00
## 9     cntry          (Intercept)        log_gdp.z.cmc -0.03  0.00
## 10    cntry          (Intercept)    gndr.c:essround.c -0.23  0.00
## 11    cntry          (Intercept) gndr.c:log_gdp.z.cmc -0.12  0.00
## 12    cntry               gndr.c           essround.c -0.55  0.00
## 13    cntry               gndr.c        log_gdp.z.cmc  0.25  0.01
## 14    cntry               gndr.c    gndr.c:essround.c  0.06  0.00
## 15    cntry               gndr.c gndr.c:log_gdp.z.cmc -0.21  0.00
## 16    cntry           essround.c        log_gdp.z.cmc -0.44 -0.01
## 17    cntry           essround.c    gndr.c:essround.c -0.01  0.00
## 18    cntry           essround.c gndr.c:log_gdp.z.cmc -0.57  0.00
## 19    cntry        log_gdp.z.cmc    gndr.c:essround.c  0.87  0.00
## 20    cntry        log_gdp.z.cmc gndr.c:log_gdp.z.cmc -0.24 -0.01
## 21    cntry    gndr.c:essround.c gndr.c:log_gdp.z.cmc -0.53  0.00
## 22 Residual                 <NA>                 <NA>  0.99  0.97
```

``` r
anova(mod2_log_GDP,mod7_log_GDP,mod8_log_GDP)
```

```
## Data: diff_dat
## Models:
## mod2_log_GDP: FM.z ~ gndr.c + log_gdp.z.cm + gndr.c:log_gdp.z.cm + (gndr.c | cntry)
## mod7_log_GDP: FM.z ~ gndr.c + log_gdp.z.cm + gndr.c:log_gdp.z.cm + essround.c + gndr.c:essround.c + (gndr.c + essround.c + gndr.c:essround.c | cntry)
## mod8_log_GDP: FM.z ~ gndr.c + log_gdp.z.cm + gndr.c:log_gdp.z.cm + essround.c + gndr.c:essround.c + log_gdp.z.cmc + log_gdp.z.cmc:gndr.c + (gndr.c + essround.c + gndr.c:essround.c + log_gdp.z.cmc + log_gdp.z.cmc:gndr.c | cntry)
##              npar    AIC    BIC logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod2_log_GDP    8 176188 176261 -88086    176172                         
## mod7_log_GDP   17 175211 175365 -87588    175177 995.73  9  < 2.2e-16 ***
## mod8_log_GDP   30 175183 175454 -87561    175123  53.90 13  6.297e-07 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

### Examine trends


``` r
# gender specific change over time

change_mod8_log_GDP<-emmeans(mod8_log_GDP,specs="essround.c",by="gndr.c",
                              at=list(gndr.c=c(-0.5,0.5),
                                      essround.c=rev(range(diff_dat$essround.c)),
                                      log_gdp.z.cmc=0,
                                      log_gdp.z.cm=0),
                              disable.pbkrtest=T,
                              lmerTest.limit = 400000,infer=c(T,T),adjust="none")
change_mod8_log_GDP
```

```
## gndr.c = -0.5:
##  essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.120 0.0680 27.3  -0.0191    0.260   1.770  0.0879
##        -4.5  0.394 0.0543 17.2   0.2799    0.509   7.267  <.0001
## 
## gndr.c =  0.5:
##  essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.466 0.0593 28.5   0.3444    0.587   7.854  <.0001
##        -4.5  0.787 0.0602 22.5   0.6625    0.912  13.067  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_mod8_log_GDP,adjust="none",infer=c(T,T))
```

```
## gndr.c = -0.5:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   -0.274 0.0831 20.3   -0.447   -0.101  -3.296  0.0036
## 
## gndr.c =  0.5:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   -0.321 0.0823 22.5   -0.492   -0.151  -3.908  0.0007
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# change in gender differences over time

change_in_diff_mod8_log_GDP<-emmeans(mod8_log_GDP,specs=c("gndr.c","essround.c"),
                                      at=list(gndr.c=c(-0.5,0.5),
                                              essround.c=rev(range(diff_dat$essround.c)),
                                              log_gdp.z.cmc=0,
                                              log_gdp.z.cm=0),
                                      disable.pbkrtest=T,
                                      lmerTest.limit = 400000,infer=c(T,T),adjust="none")
change_in_diff_mod8_log_GDP
```

```
##  gndr.c essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5        4.5  0.120 0.0680 27.3  -0.0191    0.260   1.770  0.0879
##     0.5        4.5  0.466 0.0593 28.5   0.3444    0.587   7.854  <.0001
##    -0.5       -4.5  0.394 0.0543 17.2   0.2799    0.509   7.267  <.0001
##     0.5       -4.5  0.787 0.0602 22.5   0.6625    0.912  13.067  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod8_log_GDP,adjust="none",infer=c(T,T))
```

```
##  contrast                                                 estimate     SE   df lower.CL upper.CL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.3455 0.0281 25.5   -0.403   -0.288
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  -0.2739 0.0831 20.3   -0.447   -0.101
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.6669 0.0897 22.4   -0.853   -0.481
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.0715 0.0736 19.8   -0.082    0.225
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      -0.3215 0.0823 22.5   -0.492   -0.151
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.3930 0.0277 26.8   -0.450   -0.336
##  t.ratio p.value
##  -12.300  <.0001
##   -3.296  0.0036
##   -7.436  <.0001
##    0.972  0.3426
##   -3.908  0.0007
##  -14.178  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# custom contrast to obtain the estimates for the entire span

changes_in_diff <- list(
  diff_ESS10 = c(-1,1,0,0),
  diff_ESS1 = c(0,0,-1,1))

diff_mod8_log_GDP<-contrast(change_in_diff_mod8_log_GDP,method = changes_in_diff,adjust="none",infer=c(T,T))

# gender differences at ESS1 and ESS10
diff_mod8_log_GDP
```

```
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10    0.345 0.0281 25.5    0.288    0.403  12.300  <.0001
##  diff_ESS1     0.393 0.0277 26.8    0.336    0.450  14.178  <.0001
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
##  diff_ESS10 - diff_ESS1  -0.0475 0.0421 24.6   -0.134   0.0393  -1.129  0.2699
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
## R version 4.5.1 (2025-06-13 ucrt)
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
##  [1] stringr_1.5.2         tidyr_1.3.1           r2mlm_0.3.8           nlme_3.1-168         
##  [5] Hmisc_5.2-4           ggpubr_0.6.2          metafor_4.8-0         numDeriv_2016.8-1.1  
##  [9] metadat_1.4-0         lmerTest_3.1-3        ggflags_0.0.4         finalfit_1.1.0       
## [13] ggplot2_4.0.0         MetBrewer_0.2.0       vjihelpers_0.0.0.9000 emmeans_1.11.2-8     
## [17] lme4_1.1-37           Matrix_1.7-3          dplyr_1.1.4           rio_1.2.4            
## [21] multid_1.0.2.9000     knitr_1.50            rmarkdown_2.30       
## 
## loaded via a namespace (and not attached):
##   [1] mnormt_2.1.1       Rdpack_2.6.4       pROC_1.19.0.1      gridExtra_2.3      readxl_1.4.5      
##   [6] rlang_1.1.6        magrittr_2.0.4     rockchalk_1.8.157  compiler_4.5.1     mgcv_1.9-3        
##  [11] png_0.1-8          vctrs_0.6.5        quadprog_1.5-8     crayon_1.5.3       pkgconfig_2.0.3   
##  [16] shape_1.4.6.1      fastmap_1.2.0      backports_1.5.0    labeling_0.4.3     pbivnorm_0.6.0    
##  [21] nloptr_2.2.1       purrr_1.1.0        xfun_0.53          glmnet_4.1-10      jomo_2.7-6        
##  [26] cachem_1.1.0       kutils_1.73        jsonlite_2.0.0     pan_1.9            jpeg_0.1-11       
##  [31] psych_2.5.6        lavaan_0.6-20      parallel_4.5.1     broom_1.0.10       cluster_2.1.8.1   
##  [36] R6_2.6.1           bslib_0.9.0        stringi_1.8.7      RColorBrewer_1.1-3 car_3.1-3         
##  [41] boot_1.3-31        rpart_4.1.24       cellranger_1.1.0   jquerylib_0.1.4    estimability_1.5.1
##  [46] Rcpp_1.1.0         iterators_1.0.14   base64enc_0.1-3    R.utils_2.13.0     splines_4.5.1     
##  [51] nnet_7.3-20        tidyselect_1.2.1   rstudioapi_0.17.1  abind_1.4-8        yaml_2.3.10       
##  [56] codetools_0.2-20   lattice_0.22-7     tibble_3.3.0       plyr_1.8.9         withr_3.0.2       
##  [61] S7_0.2.0           evaluate_1.0.5     foreign_0.8-90     survival_3.8-3     zip_2.3.3         
##  [66] pillar_1.11.1      carData_3.0-5      mice_3.18.0        stats4_4.5.1       checkmate_2.3.3   
##  [71] foreach_1.5.2      reformulas_0.4.1   generics_0.1.4     grImport2_0.3-3    mathjaxr_1.8-0    
##  [76] scales_1.4.0       minqa_1.2.8        xtable_1.8-4       glue_1.8.0         tools_4.5.1       
##  [81] data.table_1.17.8  openxlsx_4.2.8     ggsignif_0.6.4     forcats_1.0.1      XML_3.99-0.19     
##  [86] mvtnorm_1.3-3      cowplot_1.2.0      grid_4.5.1         rbibutils_2.3      colorspace_2.1-2  
##  [91] htmlTable_2.4.3    Formula_1.2-5      cli_3.6.5          gtable_0.3.6       R.methodsS3_1.8.2 
##  [96] rstatix_0.7.3      sass_0.4.10        digest_0.6.37      htmlwidgets_1.6.4  farver_2.1.2      
## [101] htmltools_0.5.8.1  R.oo_1.27.1        lifecycle_1.0.4    mitml_0.4-5        MASS_7.3-65
```

