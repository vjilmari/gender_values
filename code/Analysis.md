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
##  AL_6  AT_1  AT_2  AT_3  AT_7  AT_8  AT_9  BE_1 BE_10  BE_2  BE_3  BE_4  BE_5  BE_6  BE_7  BE_8  BE_9 
##  1117  2254  2198  2348  1795  1993  2477  1830  1334  1771  1796  1754  1699  1862  1767  1759  1756 
## BG_10  BG_3  BG_4  BG_5  BG_6  BG_9  CH_1 CH_10  CH_2  CH_3  CH_4  CH_5  CH_6  CH_7  CH_8  CH_9  CY_3 
##  2697  1295  2144  2371  2179  1926  2024  1505  2110  1780  1753  1491  1483  1521  1504  1517   978 
##  CY_4  CY_5  CY_6  CY_9  CZ_1 CZ_10  CZ_2  CZ_4  CZ_5  CZ_6  CZ_7  CZ_8  CZ_9  DE_1  DE_2  DE_3  DE_4 
##  1210  1053  1110   773  1208  2369  2557  1986  2335  1973  1862  2252  2343  2819  2840  2884  2732 
##  DE_5  DE_6  DE_7  DE_8  DE_9  DK_1  DK_2  DK_3  DK_4  DK_5  DK_6  DK_7  DK_9 EE_10  EE_2  EE_3  EE_4 
##  3007  2935  3006  2821  2328  1470  1458  1461  1581  1564  1621  1483  1554  1538  1948  1466  1646 
##  EE_5  EE_6  EE_7  EE_8  EE_9  ES_1  ES_2  ES_3  ES_4  ES_5  ES_6  ES_7  ES_8  ES_9  FI_1 FI_10  FI_2 
##  1793  2345  2036  2007  1899  1712  1623  1847  2562  1881  1871  1907  1929  1619  1763  1561  1701 
##  FI_3  FI_4  FI_5  FI_6  FI_7  FI_8  FI_9  FR_1 FR_10  FR_2  FR_3  FR_4  FR_5  FR_6  FR_7  FR_8  FR_9 
##  1649  1901  1649  2158  2050  1903  1735  1355  1951  1699  1983  2067  1723  1960  1902  2057  1982 
##  GB_1 GB_10  GB_2  GB_3  GB_4  GB_5  GB_6  GB_7  GB_8  GB_9  GR_1 GR_10  GR_2  GR_4  GR_5 HR_10  HR_4 
##  1798  1131  1864  2353  2311  2374  2261  2231  1942  2183  2551  2768  2399  2063  2669  1564  1430 
##  HR_5  HR_9  HU_1 HU_10  HU_2  HU_3  HU_4  HU_5  HU_6  HU_7  HU_8  HU_9  IE_1 IE_10  IE_2  IE_3  IE_4 
##  1601  1781  1634  1816  1460  1462  1430  1473  1968  1520  1458  1643  1916  1751  1187  1589  1757 
##  IE_5  IE_6  IE_7  IE_8  IE_9  IL_1  IL_4  IL_5  IL_6  IL_7  IL_8 IS_10  IS_2  IS_6  IS_8  IS_9 IT_10 
##  2400  2616  2380  2746  2189  2279  2382  2212  2378  2351  2366   886   524   739   841   844  2573 
##  IT_6  IT_8  IT_9 LT_10  LT_5  LT_6  LT_7  LT_8  LT_9  LU_2  LV_4  LV_9 ME_10  ME_9 MK_10  NL_1 NL_10 
##   909  2531  2660  1606  1632  2108  2241  2079  1677  1614  1970   891  1248  1188  1400  2337  1466 
##  NL_2  NL_3  NL_4  NL_5  NL_6  NL_7  NL_8  NL_9  NO_1 NO_10  NO_2  NO_3  NO_4  NO_5  NO_6  NO_7  NO_8 
##  1858  1860  1724  1801  1828  1823  1669  1657  1819  1408  1575  1550  1391  1530  1610  1423  1530 
##  NO_9  PL_1  PL_2  PL_3  PL_4  PL_5  PL_6  PL_7  PL_8  PL_9  PT_1 PT_10  PT_2  PT_3  PT_4  PT_5  PT_6 
##  1396  2065  1683  1685  1596  1719  1866  1594  1675  1443  1482  1827  2024  2182  2337  2139  2138 
##  PT_7  PT_8  PT_9  RO_4  RS_9  RU_3  RU_4  RU_5  RU_6  RU_8  SE_1  SE_2  SE_3  SE_4  SE_5  SE_6  SE_7 
##  1242  1254  1045  2104  1969  2339  2446  2557  2429  2374  1682  1678  1604  1556  1463  1838  1761 
##  SE_8  SE_9  SI_1 SI_10  SI_2  SI_3  SI_4  SI_5  SI_6  SI_7  SI_8  SI_9 SK_10  SK_2  SK_3  SK_4  SK_5 
##  1526  1510  1488  1232  1384  1465  1257  1369  1244  1189  1295  1307  1395  1425  1711  1789  1803 
##  SK_6  SK_9  TR_2  TR_4  UA_2  UA_3  UA_4  UA_5  UA_6  XK_6 
##  1827  1061  1790  2305  1896  1885  1766  1779  2064  1244
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

![](Analysis_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

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
##     
##        LT   LV   ME   NL   NO   PL   PT   RU   SE   SI   SK   TR   UA
##   1     0    0    0 2137 1619 1865 1282    0 1482 1288    0    0    0
##   2     0    0    0 1658 1375 1483 1824    0 1478 1184 1225 1590 1696
##   3     0    0    0 1660 1350 1485 1982 2139 1404 1265 1511    0 1685
##   4     0 1770    0 1524 1191 1396 2137 2246 1356 1057 1589 2105 1566
##   5  1432    0    0 1601 1330 1519 1939 2357 1263 1169 1603    0 1579
##   6  1908    0    0 1628 1410 1666 1938 2229 1638 1044 1627    0 1864
##   7  2041    0    0 1623 1223 1394 1042    0 1561  989    0    0    0
##   8  1879    0    0 1469 1330 1475 1054 2174 1326 1095    0    0    0
##   9  1477  691  988 1457 1196 1243  845    0 1310 1107  861    0    0
##   10 1406    0 1048 1266 1208    0 1627    0    0 1032 1195    0    0
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

![](Analysis_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

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

![](Analysis_files/figure-html/unnamed-chunk-9-2.png)<!-- -->

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
##    vars  n     mean      sd   median  trimmed      mad      min      max    range skew kurtosis      se
## X1    1 33 44309.55 16988.2 40528.07 43238.06 16282.88 16406.88 85115.96 68709.08 0.48    -0.58 2957.27
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
## `summarise()` has grouped output by 'cntry'. You can override using the `.groups` argument.
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
## `summarise()` has grouped output by 'cntry'. You can override using the `.groups` argument.
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
## `summarise()` has grouped output by 'cntry'. You can override using the `.groups` argument.
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

desc_frame
```

```
## # A tibble: 33 × 9
##    cntry `n ESS rounds`     n  `FM M` `FM SD` `FM M Women` `FM SD Women` `FM M Men` `FM SD Men`
##    <chr>          <dbl> <dbl>   <dbl>   <dbl>        <dbl>         <dbl>      <dbl>       <dbl>
##  1 AT                 6 11869  0.0820   1.06        -0.168         1.04      0.355        1.00 
##  2 BE                10 15312 -0.0363   0.923       -0.225         0.922     0.164        0.880
##  3 BG                 6 11429 -0.0387   1.09        -0.252         1.12      0.198        1.00 
##  4 CH                10 14723 -0.0624   0.872       -0.278         0.857     0.164        0.829
##  5 CY                 5  4116 -0.304    1.01        -0.436         1.06     -0.161        0.939
##  6 CZ                 9 17101  0.419    1.14         0.138         1.15      0.726        1.04 
##  7 DE                 9 23579 -0.225    1.000       -0.469         0.980     0.0337       0.954
##  8 DK                 8 10594  0.0664   0.960       -0.174         0.921     0.314        0.935
##  9 EE                 9 14874 -0.0956   1.05        -0.330         1.04      0.196        0.990
## 10 ES                 9 15164 -0.452    1.02        -0.625         1.03     -0.271        0.972
## # ℹ 23 more rows
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
  select(
    Country = CLDR,
    `n ESS rounds`,
    n,
    `FM M`, `FM SD`,
    `FM M Women`, `FM SD Women`,
    `FM M Men`, `FM SD Men`,
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
cntry_desc_tbl
```

```
## # A tibble: 33 × 13
##    Country    `n ESS rounds`     n `FM M` `FM SD` `FM M Women` `FM SD Women` `FM M Men` `FM SD Men` GEI  
##    <chr>               <dbl> <dbl> <chr>  <chr>   <chr>        <chr>         <chr>      <chr>       <chr>
##  1 Austria                 6 11869 0.08   1.06    -0.17        1.04          0.35       1.00        0.91 
##  2 Belgium                10 15312 -0.04  0.92    -0.23        0.92          0.16       0.88        0.92 
##  3 Bulgaria                6 11429 -0.04  1.09    -0.25        1.12          0.20       1.00        0.78 
##  4 Switzerla…             10 14723 -0.06  0.87    -0.28        0.86          0.16       0.83        0.95 
##  5 Cyprus                  5  4116 -0.30  1.01    -0.44        1.06          -0.16      0.94        0.79 
##  6 Czechia                 9 17101 0.42   1.14    0.14         1.15          0.73       1.04        0.86 
##  7 Germany                 9 23579 -0.22  1.00    -0.47        0.98          0.03       0.95        0.91 
##  8 Denmark                 8 10594 0.07   0.96    -0.17        0.92          0.31       0.94        0.96 
##  9 Estonia                 9 14874 -0.10  1.05    -0.33        1.04          0.20       0.99        0.85 
## 10 Spain                   9 15164 -0.45  1.02    -0.63        1.03          -0.27      0.97        0.91 
## # ℹ 23 more rows
## # ℹ 3 more variables: GGGI <chr>, GDI <chr>, GDP <chr>
```

``` r
export(cntry_desc_tbl,"../results/cntry_desc_tbl.xlsx",overwrite=T)
```

## Country-level correlations table


``` r
cor_frame<-
  desc_frame %>%
  select(
    VBMT=`FM M`,
    VBMT_Women=`FM M Women`,
    VBMT_Men=`FM M Men`,
    GEI = gei.cm,
    GGGI = gggi.cm,
    GDI = gdi.cm,
    GDP = gdp.cm
  ) %>%
  mutate(
    log_GDP=log(GDP)
  ) %>%
  select(-GDP)

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
##   1. VBMT       0.04  0.23                                                                           
##                                                                                                      
##   2. VBMT_Women -0.16 0.25 .99                                                                       
##                            [.97, .99]                                                                
##                                                                                                      
##   3. VBMT_Men   0.26  0.23 .98          .93                                                          
##                            [.96, .99]   [.87, .97]                                                   
##                                                                                                      
##   4. GEI        0.87  0.07 -.37         -.45         -.29                                            
##                            [-.64, -.02] [-.69, -.12] [-.58, .07]                                     
##                                                                                                      
##   5. GGGI       0.73  0.05 -.46         -.53         -.36         .73                                
##                            [-.69, -.13] [-.74, -.23] [-.62, -.01] [.52, .86]                         
##                                                                                                      
##   6. GDI        0.99  0.03 .10          .08          .17          .07         .20                    
##                            [-.25, .43]  [-.27, .41]  [-.18, .49]  [-.29, .41] [-.16, .51]            
##                                                                                                      
##   7. log_GDP    10.62 0.40 -.34         -.40         -.29         .75         .67         -.22       
##                            [-.61, .00]  [-.65, -.06] [-.58, .06]  [.55, .87]  [.42, .82]  [-.53, .13]
##                                                                                                      
## 
## Note. M and SD are used to represent mean and standard deviation, respectively.
## Values in square brackets indicate the 95% confidence interval.
## The confidence interval is a plausible range of population correlations 
## that could have caused the sample correlation (Cumming, 2014).
## 
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
## 1188657.1 1188689.7 -594325.6 1188651.1    392765 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.5011 -0.6198  0.0009  0.6002 10.1532 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  cntry    (Intercept) 0.05692  0.2386  
##  Residual             1.06236  1.0307  
## Number of obs: 392768, groups:  cntry, 33
## 
## Fixed effects:
##             Estimate Std. Error       df t value Pr(>|t|)
## (Intercept)  0.04841    0.04158 32.97671   1.164    0.253
```

``` r
getVC(mod0)
```

```
##        grp        var1 var2 sdcor vcov
## 1    cntry (Intercept) <NA>  0.24 0.06
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
## mean variation  0.05085636     NA       1
## sigma2          0.94914364      1      NA
## 
## $R2s
##          total within between
## f1  0.00000000      0      NA
## f2  0.00000000     NA       0
## v   0.00000000      0      NA
## m   0.05085636     NA       1
## f   0.00000000     NA      NA
## fv  0.00000000      0      NA
## fvm 0.05085636     NA      NA
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
## 1171109.8 1171153.3 -585550.9 1171101.8    392764 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.3384 -0.6154  0.0034  0.6020  9.8999 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  cntry    (Intercept) 0.0585   0.2419  
##  Residual             1.0159   1.0079  
## Number of obs: 392768, groups:  cntry, 33
## 
## Fixed effects:
##              Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept) 5.913e-02  4.215e-02 3.298e+01   1.403     0.17    
## gndr.c      4.308e-01  3.216e-03 3.927e+05 133.968   <2e-16 ***
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
##              Est.    SE         df       t    p     LL    UL
## (Intercept) 0.059 0.042     32.977   1.403 0.17 -0.027 0.145
## gndr.c      0.431 0.003 392736.769 133.968 0.00  0.424 0.437
```

``` r
getVC(mod1)
```

```
##        grp        var1 var2 sdcor vcov
## 1    cntry (Intercept) <NA>  0.24 0.06
## 2 Residual        <NA> <NA>  1.01 1.02
```

``` r
r2mlm(mod1,bargraph = F)
```

```
## $Decompositions
##                      total
## fixed           0.04108843
## slope variation 0.00000000
## mean variation  0.05220964
## sigma2          0.90670192
## 
## $R2s
##          total
## f   0.04108843
## v   0.00000000
## m   0.05220964
## fv  0.04108843
## fvm 0.09329808
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
## 1170511.6 1170576.9 -585249.8 1170499.6    392762 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.4447 -0.6151  0.0035  0.6022  9.8992 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr 
##  cntry    (Intercept) 0.058318 0.2415        
##           gndr.c      0.007657 0.0875   -0.17
##  Residual             1.014117 1.0070        
## Number of obs: 392768, groups:  cntry, 33
## 
## Fixed effects:
##             Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)  0.05884    0.04208 32.97483   1.398    0.171    
## gndr.c       0.41941    0.01570 32.73337  26.718   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##        (Intr)
## gndr.c -0.164
```

``` r
getFE(mod2,round=3)
```

```
##              Est.    SE     df      t     p     LL    UL
## (Intercept) 0.059 0.042 32.975  1.398 0.171 -0.027 0.144
## gndr.c      0.419 0.016 32.733 26.718 0.000  0.387 0.451
```

``` r
getVC(mod2)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.24 0.06
## 2    cntry      gndr.c   <NA>  0.09 0.01
## 3    cntry (Intercept) gndr.c -0.17 0.00
## 4 Residual        <NA>   <NA>  1.01 1.01
```

``` r
r2mlm(mod2,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.039024402
## slope variation 0.001698633
## mean variation  0.052442066
## sigma2          0.906834899
## 
## $R2s
##           total
## f   0.039024402
## v   0.001698633
## m   0.052442066
## fv  0.040723035
## fvm 0.093165101
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
## mod1    4 1171110 1171153 -585551   1171102                         
## mod2    6 1170512 1170577 -585250   1170500 602.15  2  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

``` r
lvl2_var_cond_lvl1(mod2,lvl1.var = "gndr.c",lvl1.values = c(0.5,-0.5))
```

```
##   lvl1.value lvl2.cond.var lvl2.cond.sd
## 1        0.5    0.05664492    0.2380019
## 2       -0.5    0.06381993    0.2526261
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
## 1170510.5 1170564.9 -585250.3 1170500.5    392763 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.4440 -0.6150  0.0035  0.6022  9.8974 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  cntry    (Intercept) 0.058332 0.24152 
##  cntry.1  gndr.c      0.007662 0.08754 
##  Residual             1.014117 1.00703 
## Number of obs: 392768, groups:  cntry, 33
## 
## Fixed effects:
##             Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)  0.05884    0.04209 32.97586   1.398    0.171    
## gndr.c       0.41942    0.01570 32.72455  26.708   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##        (Intr)
## gndr.c 0.001
```

``` r
getFE(mod2_norecov,round=3)
```

```
##              Est.    SE     df      t     p     LL    UL
## (Intercept) 0.059 0.042 32.976  1.398 0.171 -0.027 0.144
## gndr.c      0.419 0.016 32.725 26.708 0.000  0.387 0.451
```

``` r
getVC(mod2_norecov)
```

```
##        grp        var1 var2 sdcor vcov
## 1    cntry (Intercept) <NA>  0.24 0.06
## 2  cntry.1      gndr.c <NA>  0.09 0.01
## 3 Residual        <NA> <NA>  1.01 1.01
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
## mod2_norecov    5 1170511 1170565 -585250   1170501                     
## mod2            6 1170512 1170577 -585250   1170500 0.8959  1     0.3439
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
## 1140107.0 1140193.9 -570045.5 1140091.0    384370 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.4908 -0.6160  0.0039  0.6038  9.9445 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr
##  cntry    (Intercept) 0.049759 0.22307      
##           gndr.c      0.006046 0.07775  0.04
##  Residual             1.004645 1.00232      
## Number of obs: 384378, groups:  cntry, 32
## 
## Fixed effects:
##                 Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)      0.05305    0.03948 31.99285   1.344   0.1885    
## gndr.c           0.42160    0.01427 32.05742  29.547   <2e-16 ***
## gei.z.cm        -0.09790    0.04013 32.05339  -2.440   0.0204 *  
## gndr.c:gei.z.cm  0.04151    0.01470 33.89630   2.823   0.0079 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c g.z.cm
## gndr.c       0.037              
## gei.z.cm    -0.001  0.000       
## gndr.c:g.z.  0.000 -0.022  0.037
```

``` r
getFE(mod2_GEI,round=3)
```

```
##                   Est.    SE     df      t     p     LL     UL
## (Intercept)      0.053 0.039 31.993  1.344 0.188 -0.027  0.133
## gndr.c           0.422 0.014 32.057 29.547 0.000  0.393  0.451
## gei.z.cm        -0.098 0.040 32.053 -2.440 0.020 -0.180 -0.016
## gndr.c:gei.z.cm  0.042 0.015 33.896  2.823 0.008  0.012  0.071
```

``` r
getVC(mod2_GEI)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.22 0.05
## 2    cntry      gndr.c   <NA>  0.08 0.01
## 3    cntry (Intercept) gndr.c  0.04 0.00
## 4 Residual        <NA>   <NA>  1.00 1.00
```

``` r
r2mlm(mod2_GEI,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.046353324
## slope variation 0.001355711
## mean variation  0.044901307
## sigma2          0.907389659
## 
## $R2s
##           total
## f   0.046353324
## v   0.001355711
## m   0.044901307
## fv  0.047709035
## fvm 0.092610341
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
## Time difference of 53.36705 secs
```

``` r
round(ddsc_mod2_GEI$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.064        0.253        1.014     1.078 0.059   6470.818 0.997   0.998
## 2        0.5         0.057        0.238        1.014     1.071 0.053   5431.242 0.996   0.997
```

``` r
round(ddsc_mod2_GEI$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gei.z.cm
## means_y1           0.264 0.245    1.000           1.000    0.934           0.934   -0.315
## means_y1_scaled    1.050 0.975    1.000           1.000    0.934           0.934   -0.315
## means_y2          -0.158 0.258    0.934           0.934    1.000           1.000   -0.461
## means_y2_scaled   -0.628 1.025    0.934           0.934    1.000           1.000   -0.461
## gei.z.cm           0.000 1.000   -0.315          -0.315   -0.461          -0.461    1.000
## gei.z.cm_scaled    0.000 1.000   -0.315          -0.315   -0.461          -0.461    1.000
## diff_score         0.422 0.092    0.048           0.048   -0.311          -0.311    0.451
## diff_score_scaled  1.677 0.365    0.048           0.048   -0.311          -0.311    0.451
##                   gei.z.cm_scaled diff_score diff_score_scaled
## means_y1                   -0.315      0.048             0.048
## means_y1_scaled            -0.315      0.048             0.048
## means_y2                   -0.461     -0.311            -0.311
## means_y2_scaled            -0.461     -0.311            -0.311
## gei.z.cm                    1.000      0.451             0.451
## gei.z.cm_scaled             1.000      0.451             0.451
## diff_score                  0.451      1.000             1.000
## diff_score_scaled           0.451      1.000             1.000
```

``` r
round(ddsc_mod2_GEI$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                      -0.452 0.160 33.896  -2.823   0.008   -0.777   -0.127
## w_11                         -0.119 0.041 32.120  -2.928   0.006   -0.201   -0.036
## w_21                         -0.077 0.041 32.085  -1.879   0.069   -0.161    0.006
## r_xy1                        -0.484 0.165 32.120  -2.928   0.006   -0.821   -0.147
## r_xy2                        -0.299 0.159 32.085  -1.879   0.069   -0.624    0.025
## b_11                         -0.472 0.161 32.120  -2.928   0.006   -0.800   -0.144
## b_21                         -0.307 0.163 32.085  -1.879   0.069   -0.640    0.026
## main_effect                  -0.098 0.040 32.053  -2.440   0.020   -0.180   -0.016
## moderator_effect              0.422 0.014 32.057  29.547   0.000    0.393    0.451
## interaction                   0.042 0.015 33.896   2.823   0.008    0.012    0.071
## q_b11_b21                    -0.196    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                  -0.219    NA     NA      NA      NA       NA       NA
## cross_over_point            -10.156    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.056 0.043 32.169  -1.304   0.202   -0.144    0.032
## interaction_vs_main_bscale   -0.224 0.172 32.169  -1.304   0.202   -0.575    0.126
## interaction_vs_main_rscale   -0.207 0.164 32.177  -1.260   0.217   -0.542    0.128
## dadas                        -0.154 0.082 32.085  -1.879   0.965   -0.322    0.013
## dadas_bscale                 -0.614 0.327 32.085  -1.879   0.965   -1.279    0.052
## dadas_rscale                 -0.599 0.319 32.085  -1.879   0.965   -1.248    0.050
## abs_diff                      0.042 0.015 33.896   2.823   0.004    0.012    0.071
## abs_sum                       0.196 0.080 32.053   2.440   0.010    0.032    0.359
## abs_diff_bscale               0.165 0.058 33.896   2.823   0.004    0.046    0.284
## abs_sum_bscale                0.779 0.319 32.053   2.440   0.010    0.129    1.429
## abs_diff_rscale               0.185 0.059 34.065   3.141   0.002    0.065    0.304
## abs_sum_rscale                0.783 0.319 32.054   2.453   0.010    0.133    1.434
```

``` r
round(ddsc_mod2_GEI$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.004 -0.170  0.896  1.000  0.344
```

``` r
d_GEI<-ddsc_mod2_GEI$ddsc_sem_fit$data

ddsc_sem_GEI<-
  ddsc_sem(data=d_GEI,x = "gei.z.cm",y1="means_y2",
           y2="means_y1",q_sesoi = .10)
round(ddsc_sem_GEI$results,3)
```

```
##                                     est    se       z pvalue ci.lower ci.upper
## r_xy1_y2                         -0.451 0.158  -2.861  0.004   -0.761   -0.142
## r_xy1                            -0.461 0.157  -2.938  0.003   -0.768   -0.153
## r_xy2                            -0.315 0.168  -1.880  0.060   -0.644    0.013
## b_11                             -0.472 0.161  -2.938  0.003   -0.787   -0.157
## b_21                             -0.307 0.164  -1.880  0.060   -0.628    0.013
## b_10                             -0.628 0.158  -3.968  0.000   -0.938   -0.318
## b_20                              1.050 0.161   6.522  0.000    0.734    1.365
## res_cov_y1_y2                     0.763 0.197   3.868  0.000    0.377    1.150
## diff_b10_b20                     -1.677 0.057 -29.579  0.000   -1.789   -1.566
## diff_b11_b21                     -0.165 0.058  -2.861  0.004   -0.278   -0.052
## diff_rxy1_rxy2                   -0.146 0.059  -2.485  0.013   -0.260   -0.031
## q_b11_b21                        -0.195 0.073  -2.658  0.008   -0.339   -0.051
## q_rxy1_rxy2                      -0.172 0.070  -2.472  0.013   -0.308   -0.036
## cross_over_point                -10.175 3.572  -2.848  0.004  -17.176   -3.173
## sum_b11_b21                      -0.780 0.319  -2.443  0.015   -1.405   -0.154
## main_effect                      -0.390 0.160  -2.443  0.015   -0.703   -0.077
## interaction_vs_main_effect       -0.225 0.172  -1.306  0.192   -0.563    0.113
## diff_abs_b11_abs_b21              0.165 0.058   2.861  0.004    0.052    0.278
## abs_diff_b11_b21                  0.165 0.058   2.861  0.002    0.052    0.278
## abs_sum_b11_b21                   0.780 0.319   2.443  0.007    0.154    1.405
## dadas                            -0.615 0.327  -1.880  0.970   -1.256    0.026
## q_r_equivalence                   0.072 0.070   1.035  0.850       NA       NA
## q_b_equivalence                   0.095 0.073   1.297  0.903       NA       NA
## cross_over_point_equivalence     10.175 3.572   2.848  0.998       NA       NA
## cross_over_point_minimal_effect  10.175 3.572   2.848  0.002       NA       NA
```

``` r
round(ddsc_sem_GEI$variance_test,3)
```

```
##             est    se      z pvalue ci.lower ci.upper
## cov_y1y2  0.904 0.234  3.862  0.000    0.445    1.363
## var_y1    1.017 0.254  4.000  0.000    0.519    1.515
## var_y2    0.921 0.230  4.000  0.000    0.470    1.372
## var_diff  0.096 0.124  0.776  0.438   -0.147    0.340
## var_ratio 1.105 0.139  7.943  0.000    0.832    1.377
## cor_y1y2  0.934 0.022 41.686  0.000    0.891    0.978
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

pflag_comb<-
  ggarrange(p1.FM.flags,p2.FM.flags,align = "v",
            ncol=1,nrow=2,heights=c(2,1))
```

```
## Warning in geom_text(inherit.aes = F, aes(x = 0.65, y = -0.45, label = coef_q, : All aesthetics have length 1, but the data has 682 rows.
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

![](Analysis_files/figure-html/unnamed-chunk-22-1.png)<!-- -->

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
##  835098.3  835182.7 -417541.2  835082.3    280638 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.4749 -0.6139  0.0057  0.6032  8.0857 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr
##  cntry    (Intercept) 0.055856 0.2363       
##           gndr.c      0.006038 0.0777   0.00
##  Residual             1.007836 1.0039       
## Number of obs: 280646, groups:  cntry, 33
## 
## Fixed effects:
##                  Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)       0.07622    0.04121 32.97252   1.850  0.07335 .  
## gndr.c            0.41058    0.01426 31.97619  28.800  < 2e-16 ***
## gggi.z.cm        -0.13836    0.04188 33.07113  -3.304  0.00230 ** 
## gndr.c:gggi.z.cm  0.05014    0.01481 34.78082   3.385  0.00178 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c ggg.z.
## gndr.c      -0.002              
## gggi.z.cm    0.000  0.000       
## gndr.c:gg..  0.000 -0.011 -0.002
```

``` r
getFE(mod2_GGGI,round=3)
```

```
##                    Est.    SE     df      t     p     LL     UL
## (Intercept)       0.076 0.041 32.973  1.850 0.073 -0.008  0.160
## gndr.c            0.411 0.014 31.976 28.800 0.000  0.382  0.440
## gggi.z.cm        -0.138 0.042 33.071 -3.304 0.002 -0.224 -0.053
## gndr.c:gggi.z.cm  0.050 0.015 34.781  3.385 0.002  0.020  0.080
```

``` r
getVC(mod2_GGGI)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.24 0.06
## 2    cntry      gndr.c   <NA>  0.08 0.01
## 3    cntry (Intercept) gndr.c  0.00 0.00
## 4 Residual        <NA>   <NA>  1.00 1.01
```

``` r
r2mlm(mod2_GGGI,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.049696629
## slope variation 0.001335796
## mean variation  0.049845683
## sigma2          0.899121891
## 
## $R2s
##           total
## f   0.049696629
## v   0.001335796
## m   0.049845683
## fv  0.051032425
## fvm 0.100878109
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
## Time difference of 52.52572 secs
```

``` r
round(ddsc_mod2_GGGI$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.064        0.253        1.014     1.078 0.059   6470.818 0.997   0.998
## 2        0.5         0.057        0.238        1.014     1.071 0.053   5431.242 0.996   0.997
```

``` r
round(ddsc_mod2_GGGI$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gggi.z.cm
## means_y1           0.281 0.269    1.000           1.000    0.941           0.941    -0.418
## means_y1_scaled    0.995 0.955    1.000           1.000    0.941           0.941    -0.418
## means_y2          -0.128 0.294    0.941           0.941    1.000           1.000    -0.558
## means_y2_scaled   -0.453 1.043    0.941           0.941    1.000           1.000    -0.558
## gggi.z.cm          0.000 1.000   -0.418          -0.418   -0.558          -0.558     1.000
## gggi.z.cm_scaled   0.000 1.000   -0.418          -0.418   -0.558          -0.558     1.000
## diff_score         0.408 0.100   -0.075          -0.075   -0.408          -0.408     0.518
## diff_score_scaled  1.449 0.354   -0.075          -0.075   -0.408          -0.408     0.518
##                   gggi.z.cm_scaled diff_score diff_score_scaled
## means_y1                    -0.418     -0.075            -0.075
## means_y1_scaled             -0.418     -0.075            -0.075
## means_y2                    -0.558     -0.408            -0.408
## means_y2_scaled             -0.558     -0.408            -0.408
## gggi.z.cm                    1.000      0.518             0.518
## gggi.z.cm_scaled             1.000      0.518             0.518
## diff_score                   0.518      1.000             1.000
## diff_score_scaled            0.518      1.000             1.000
```

``` r
round(ddsc_mod2_GGGI$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                      -0.503 0.149 34.781  -3.385   0.002   -0.804   -0.201
## w_11                         -0.163 0.043 33.402  -3.842   0.001   -0.250   -0.077
## w_21                         -0.113 0.043 33.335  -2.665   0.012   -0.200   -0.027
## r_xy1                        -0.607 0.158 33.402  -3.842   0.001   -0.928   -0.286
## r_xy2                        -0.385 0.145 33.335  -2.665   0.012   -0.679   -0.091
## b_11                         -0.580 0.151 33.402  -3.842   0.001   -0.887   -0.273
## b_21                         -0.402 0.151 33.335  -2.665   0.012   -0.709   -0.095
## main_effect                  -0.138 0.042 33.071  -3.304   0.002   -0.224   -0.053
## moderator_effect              0.411 0.014 31.976  28.800   0.000    0.382    0.440
## interaction                   0.050 0.015 34.781   3.385   0.002    0.020    0.080
## q_b11_b21                    -0.236    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                  -0.298    NA     NA      NA      NA       NA       NA
## cross_over_point             -8.189    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.088 0.044 34.026  -1.987   0.055   -0.178    0.002
## interaction_vs_main_bscale   -0.313 0.158 34.026  -1.987   0.055   -0.633    0.007
## interaction_vs_main_rscale   -0.274 0.145 34.142  -1.890   0.067   -0.569    0.021
## dadas                        -0.227 0.085 33.335  -2.665   0.994   -0.400   -0.054
## dadas_bscale                 -0.804 0.302 33.335  -2.665   0.994   -1.418   -0.190
## dadas_rscale                 -0.770 0.289 33.335  -2.665   0.994   -1.358   -0.182
## abs_diff                      0.050 0.015 34.781   3.385   0.001    0.020    0.080
## abs_sum                       0.277 0.084 33.071   3.304   0.001    0.106    0.447
## abs_diff_bscale               0.178 0.053 34.781   3.385   0.001    0.071    0.285
## abs_sum_bscale                0.982 0.297 33.071   3.304   0.001    0.377    1.587
## abs_diff_rscale               0.222 0.054 35.420   4.081   0.000    0.111    0.332
## abs_sum_rscale                0.992 0.298 33.073   3.330   0.001    0.386    1.598
```

``` r
round(ddsc_mod2_GGGI$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.004 -0.170  0.896  1.000  0.344
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
## r_xy1_y2                        -0.518 0.149  -3.480  0.001   -0.810   -0.226
## r_xy1                           -0.558 0.144  -3.866  0.000   -0.841   -0.275
## r_xy2                           -0.418 0.158  -2.642  0.008   -0.728   -0.108
## b_11                            -0.582 0.151  -3.866  0.000   -0.878   -0.287
## b_21                            -0.399 0.151  -2.642  0.008   -0.695   -0.103
## b_10                            -0.453 0.148  -3.056  0.002   -0.744   -0.163
## b_20                             0.995 0.149   6.692  0.000    0.704    1.287
## res_cov_y1_y2                    0.684 0.174   3.932  0.000    0.343    1.024
## diff_b10_b20                    -1.449 0.052 -27.934  0.000   -1.550   -1.347
## diff_b11_b21                    -0.183 0.053  -3.480  0.001   -0.286   -0.080
## diff_rxy1_rxy2                  -0.140 0.055  -2.576  0.010   -0.247   -0.034
## q_b11_b21                       -0.243 0.086  -2.845  0.004   -0.411   -0.076
## q_rxy1_rxy2                     -0.185 0.072  -2.561  0.010   -0.327   -0.044
## cross_over_point                -7.904 2.289  -3.454  0.001  -12.389   -3.418
## sum_b11_b21                     -0.981 0.297  -3.304  0.001   -1.564   -0.399
## main_effect                     -0.491 0.149  -3.304  0.001   -0.782   -0.200
## interaction_vs_main_effect      -0.307 0.158  -1.946  0.052   -0.617    0.002
## diff_abs_b11_abs_b21             0.183 0.053   3.480  0.001    0.080    0.286
## abs_diff_b11_b21                 0.183 0.053   3.480  0.000    0.080    0.286
## abs_sum_b11_b21                  0.981 0.297   3.304  0.000    0.399    1.564
## dadas                           -0.798 0.302  -2.642  0.996   -1.390   -0.206
## q_r_equivalence                  0.085 0.072   1.179  0.881       NA       NA
## q_b_equivalence                  0.143 0.086   1.677  0.953       NA       NA
## cross_over_point_equivalence     7.904 2.289   3.454  1.000       NA       NA
## cross_over_point_minimal_effect  7.904 2.289   3.454  0.000       NA       NA
```

``` r
round(ddsc_sem_GGGI$variance_test,3)
```

```
##             est    se      z pvalue ci.lower ci.upper
## cov_y1y2  0.909 0.231  3.937  0.000    0.456    1.362
## var_y1    1.055 0.260  4.062  0.000    0.546    1.564
## var_y2    0.884 0.218  4.062  0.000    0.458    1.311
## var_diff  0.171 0.121  1.407  0.159   -0.067    0.408
## var_ratio 1.193 0.140  8.494  0.000    0.918    1.468
## cor_y1y2  0.941 0.020 47.281  0.000    0.902    0.980
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

![](Analysis_files/figure-html/unnamed-chunk-25-1.png)<!-- -->

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
##   1170512   1170599   -585248   1170496    392760 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.4428 -0.6151  0.0035  0.6022  9.8969 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr 
##  cntry    (Intercept) 0.057267 0.23930       
##           gndr.c      0.007183 0.08475  -0.22
##  Residual             1.014116 1.00703       
## Number of obs: 392768, groups:  cntry, 33
## 
## Fixed effects:
##                 Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)      0.05886    0.04170 32.97307   1.412    0.167    
## gndr.c           0.41932    0.01523 33.40436  27.528   <2e-16 ***
## gdi.z.cm         0.03313    0.04237 33.05236   0.782    0.440    
## gndr.c:gdi.z.cm  0.02483    0.01573 35.65795   1.578    0.123    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c gd.z.c
## gndr.c      -0.209              
## gdi.z.cm     0.000  0.000       
## gndr.c:gd..  0.001 -0.005 -0.205
```

``` r
getFE(mod2_GDI,round=3)
```

```
##                  Est.    SE     df      t     p     LL    UL
## (Intercept)     0.059 0.042 32.973  1.412 0.167 -0.026 0.144
## gndr.c          0.419 0.015 33.404 27.528 0.000  0.388 0.450
## gdi.z.cm        0.033 0.042 33.052  0.782 0.440 -0.053 0.119
## gndr.c:gdi.z.cm 0.025 0.016 35.658  1.578 0.123 -0.007 0.057
```

``` r
getVC(mod2_GDI)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.24 0.06
## 2    cntry      gndr.c   <NA>  0.08 0.01
## 3    cntry (Intercept) gndr.c -0.22 0.00
## 4 Residual        <NA>   <NA>  1.01 1.01
```

``` r
r2mlm(mod2_GDI,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.039266180
## slope variation 0.001594848
## mean variation  0.051603943
## sigma2          0.907535029
## 
## $R2s
##           total
## f   0.039266180
## v   0.001594848
## m   0.051603943
## fv  0.040861028
## fvm 0.092464971
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
## Time difference of 52.76935 secs
```

``` r
round(ddsc_mod2_GDI$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.064        0.253        1.014     1.078 0.059   6470.818 0.997   0.998
## 2        0.5         0.057        0.238        1.014     1.071 0.053   5431.242 0.996   0.997
```

``` r
round(ddsc_mod2_GDI$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gdi.z.cm
## means_y1           0.268 0.243    1.000           1.000    0.934           0.934    0.191
## means_y1_scaled    1.074 0.971    1.000           1.000    0.934           0.934    0.191
## means_y2          -0.151 0.257    0.934           0.934    1.000           1.000    0.078
## means_y2_scaled   -0.603 1.028    0.934           0.934    1.000           1.000    0.078
## gdi.z.cm           0.000 1.000    0.191           0.191    0.078           0.078    1.000
## gdi.z.cm_scaled    0.000 1.000    0.191           0.191    0.078           0.078    1.000
## diff_score         0.419 0.092    0.028           0.028   -0.332          -0.332    0.288
## diff_score_scaled  1.676 0.368    0.028           0.028   -0.332          -0.332    0.288
##                   gdi.z.cm_scaled diff_score diff_score_scaled
## means_y1                    0.191      0.028             0.028
## means_y1_scaled             0.191      0.028             0.028
## means_y2                    0.078     -0.332            -0.332
## means_y2_scaled             0.078     -0.332            -0.332
## gdi.z.cm                    1.000      0.288             0.288
## gdi.z.cm_scaled             1.000      0.288             0.288
## diff_score                  0.288      1.000             1.000
## diff_score_scaled           0.288      1.000             1.000
```

``` r
round(ddsc_mod2_GDI$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                      -0.270 0.171 35.658  -1.578   0.123   -0.617    0.077
## w_11                          0.021 0.045 33.083   0.464   0.646   -0.070    0.112
## w_21                          0.046 0.041 33.122   1.098   0.280   -0.039    0.130
## r_xy1                         0.085 0.184 33.083   0.464   0.646   -0.289    0.460
## r_xy2                         0.177 0.161 33.122   1.098   0.280   -0.151    0.505
## b_11                          0.083 0.179 33.083   0.464   0.646   -0.281    0.447
## b_21                          0.182 0.166 33.122   1.098   0.280   -0.155    0.520
## main_effect                   0.033 0.042 33.052   0.782   0.440   -0.053    0.119
## moderator_effect              0.419 0.015 33.404  27.528   0.000    0.388    0.450
## interaction                   0.025 0.016 35.658   1.578   0.123   -0.007    0.057
## q_b11_b21                    -0.101    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                  -0.093    NA     NA      NA      NA       NA       NA
## cross_over_point            -16.888    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.008 0.048 33.194  -0.172   0.864   -0.106    0.090
## interaction_vs_main_bscale   -0.033 0.193 33.194  -0.172   0.864   -0.425    0.359
## interaction_vs_main_rscale   -0.039 0.203 33.182  -0.195   0.847   -0.452    0.373
## dadas                        -0.041 0.089 33.083  -0.464   0.677   -0.223    0.140
## dadas_bscale                 -0.166 0.358 33.083  -0.464   0.677   -0.893    0.561
## dadas_rscale                 -0.171 0.368 33.083  -0.464   0.677   -0.920    0.578
## abs_diff                      0.025 0.016 35.658   1.578   0.062   -0.007    0.057
## abs_sum                       0.066 0.085 33.052   0.782   0.220   -0.106    0.239
## abs_diff_bscale               0.099 0.063 35.658   1.578   0.062   -0.028    0.227
## abs_sum_bscale                0.265 0.339 33.052   0.782   0.220   -0.425    0.955
## abs_diff_rscale               0.092 0.066 35.356   1.396   0.086   -0.042    0.225
## abs_sum_rscale                0.263 0.340 33.052   0.772   0.223   -0.429    0.954
```

``` r
round(ddsc_mod2_GDI$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.004 -0.170  0.896  1.000  0.344
```

``` r
d_GDI<-ddsc_mod2_GDI$ddsc_sem_fit$data

ddsc_sem_GDI<-
  ddsc_sem(data=d_GDI,x = "gdi.z.cm",y1="means_y2",
           y2="means_y1",q_sesoi = .10)
round(ddsc_sem_GDI$results,3)
```

```
##                                     est    se       z pvalue ci.lower ci.upper
## r_xy1_y2                         -0.288 0.167  -1.728  0.084   -0.615    0.039
## r_xy1                             0.078 0.174   0.447  0.655   -0.263    0.418
## r_xy2                             0.191 0.171   1.120  0.263   -0.144    0.526
## b_11                              0.080 0.178   0.447  0.655   -0.270    0.430
## b_21                              0.186 0.166   1.120  0.263   -0.139    0.511
## b_10                             -0.603 0.176  -3.428  0.001   -0.947   -0.258
## b_20                              1.074 0.163   6.575  0.000    0.754    1.394
## res_cov_y1_y2                     0.890 0.226   3.933  0.000    0.446    1.333
## diff_b10_b20                     -1.676 0.060 -27.750  0.000   -1.795   -1.558
## diff_b11_b21                     -0.106 0.061  -1.728  0.084   -0.226    0.014
## diff_rxy1_rxy2                   -0.114 0.060  -1.893  0.058   -0.232    0.004
## q_b11_b21                        -0.108 0.062  -1.748  0.080   -0.229    0.013
## q_rxy1_rxy2                      -0.116 0.061  -1.889  0.059   -0.237    0.004
## cross_over_point                -15.812 9.167  -1.725  0.085  -33.779    2.155
## sum_b11_b21                       0.265 0.339   0.783  0.434   -0.399    0.930
## main_effect                       0.133 0.170   0.783  0.434   -0.200    0.465
## interaction_vs_main_effect       -0.027 0.192  -0.139  0.889   -0.403    0.350
## diff_abs_b11_abs_b21             -0.106 0.061  -1.728  0.084   -0.226    0.014
## abs_diff_b11_b21                  0.106 0.061   1.728  0.042   -0.014    0.226
## abs_sum_b11_b21                   0.265 0.339   0.783  0.217   -0.399    0.930
## dadas                            -0.159 0.357  -0.447  0.672   -0.859    0.540
## q_r_equivalence                   0.016 0.061   0.261  0.603       NA       NA
## q_b_equivalence                   0.008 0.062   0.130  0.552       NA       NA
## cross_over_point_equivalence     15.812 9.167   1.725  0.958       NA       NA
## cross_over_point_minimal_effect  15.812 9.167   1.725  0.042       NA       NA
```

``` r
round(ddsc_sem_GDI$variance_test,3)
```

```
##             est    se      z pvalue ci.lower ci.upper
## cov_y1y2  0.904 0.231  3.921  0.000    0.452    1.356
## var_y1    1.026 0.253  4.062  0.000    0.531    1.521
## var_y2    0.914 0.225  4.062  0.000    0.473    1.355
## var_diff  0.112 0.124  0.906  0.365   -0.130    0.354
## var_ratio 1.123 0.140  8.030  0.000    0.849    1.397
## cor_y1y2  0.934 0.022 41.933  0.000    0.890    0.977
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

#p1.FM.flags


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
#p2.FM.flags


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

![](Analysis_files/figure-html/unnamed-chunk-28-1.png)<!-- -->

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
## 1170508.2 1170595.3 -585246.1 1170492.2    392760 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.4445 -0.6151  0.0034  0.6022  9.8993 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr 
##  cntry    (Intercept) 0.050246 0.22416       
##           gndr.c      0.006879 0.08294  -0.06
##  Residual             1.014117 1.00703       
## Number of obs: 392768, groups:  cntry, 33
## 
## Fixed effects:
##                     Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)          0.05688    0.03908 32.98150   1.456   0.1549    
## gndr.c               0.41954    0.01493 32.60491  28.107   <2e-16 ***
## log_gdp.z.cm        -0.09021    0.03923 33.02655  -2.300   0.0279 *  
## gndr.c:log_gdp.z.cm  0.02800    0.01510 33.59345   1.854   0.0725 .  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c lg_g..
## gndr.c      -0.058              
## lg_gdp.z.cm  0.021 -0.002       
## gndr.c:l_.. -0.002  0.001 -0.057
```

``` r
getFE(mod2_log_GDP,round=3)
```

```
##                       Est.    SE     df      t     p     LL     UL
## (Intercept)          0.057 0.039 32.981  1.456 0.155 -0.023  0.136
## gndr.c               0.420 0.015 32.605 28.107 0.000  0.389  0.450
## log_gdp.z.cm        -0.090 0.039 33.027 -2.300 0.028 -0.170 -0.010
## gndr.c:log_gdp.z.cm  0.028 0.015 33.593  1.854 0.073 -0.003  0.059
```

``` r
getVC(mod2_log_GDP)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.22 0.05
## 2    cntry      gndr.c   <NA>  0.08 0.01
## 3    cntry (Intercept) gndr.c -0.06 0.00
## 4 Residual        <NA>   <NA>  1.01 1.01
```

``` r
r2mlm(mod2_log_GDP,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.044858184
## slope variation 0.001528867
## mean variation  0.045112797
## sigma2          0.908500152
## 
## $R2s
##           total
## f   0.044858184
## v   0.001528867
## m   0.045112797
## fv  0.046387051
## fvm 0.091499848
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
## not properly standardized, SD = 1.01176689233303
```

``` r
t2<-Sys.time()
t2-t1
```

```
## Time difference of 52.66602 secs
```

``` r
round(ddsc_mod2_log_GDP$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.064        0.253        1.014     1.078 0.059   6470.818 0.997   0.998
## 2        0.5         0.057        0.238        1.014     1.071 0.053   5431.242 0.996   0.997
```

``` r
round(ddsc_mod2_log_GDP$descriptives,3)
```

```
##                          M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled log_gdp.z.cm
## means_y1             0.268 0.243    1.000           1.000    0.934           0.934       -0.319
## means_y1_scaled      1.074 0.971    1.000           1.000    0.934           0.934       -0.319
## means_y2            -0.151 0.257    0.934           0.934    1.000           1.000       -0.410
## means_y2_scaled     -0.603 1.028    0.934           0.934    1.000           1.000       -0.410
## log_gdp.z.cm        -0.022 1.012   -0.319          -0.319   -0.410          -0.410        1.000
## log_gdp.z.cm_scaled  0.000 1.000   -0.319          -0.319   -0.410          -0.410        1.000
## diff_score           0.419 0.092    0.028           0.028   -0.332          -0.332        0.302
## diff_score_scaled    1.676 0.368    0.028           0.028   -0.332          -0.332        0.302
##                     log_gdp.z.cm_scaled diff_score diff_score_scaled
## means_y1                         -0.319      0.028             0.028
## means_y1_scaled                  -0.319      0.028             0.028
## means_y2                         -0.410     -0.332            -0.332
## means_y2_scaled                  -0.410     -0.332            -0.332
## log_gdp.z.cm                      1.000      0.302             0.302
## log_gdp.z.cm_scaled               1.000      0.302             0.302
## diff_score                        0.302      1.000             1.000
## diff_score_scaled                 0.302      1.000             1.000
```

``` r
round(ddsc_mod2_log_GDP$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                      -0.304 0.164 33.593  -1.854   0.073   -0.638    0.029
## w_11                         -0.104 0.040 33.051  -2.581   0.014   -0.186   -0.022
## w_21                         -0.076 0.040 33.045  -1.928   0.062   -0.157    0.004
## r_xy1                        -0.430 0.166 33.051  -2.581   0.014   -0.768   -0.091
## r_xy2                        -0.296 0.154 33.045  -1.928   0.062   -0.609    0.016
## b_11                         -0.417 0.162 33.051  -2.581   0.014   -0.746   -0.088
## b_21                         -0.305 0.158 33.045  -1.928   0.062   -0.627    0.017
## main_effect                  -0.090 0.039 33.027  -2.300   0.028   -0.170   -0.010
## moderator_effect              0.420 0.015 32.605  28.107   0.000    0.389    0.450
## interaction                   0.028 0.015 33.593   1.854   0.073   -0.003    0.059
## q_b11_b21                    -0.129    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                  -0.154    NA     NA      NA      NA       NA       NA
## cross_over_point            -14.982    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.062 0.041 33.085  -1.509   0.141   -0.146    0.022
## interaction_vs_main_bscale   -0.249 0.165 33.085  -1.509   0.141   -0.585    0.087
## interaction_vs_main_rscale   -0.230 0.156 33.090  -1.470   0.151   -0.548    0.088
## dadas                        -0.152 0.079 33.045  -1.928   0.969   -0.313    0.008
## dadas_bscale                 -0.610 0.316 33.045  -1.928   0.969   -1.254    0.034
## dadas_rscale                 -0.593 0.307 33.045  -1.928   0.969   -1.219    0.033
## abs_diff                      0.028 0.015 33.593   1.854   0.036   -0.003    0.059
## abs_sum                       0.180 0.078 33.027   2.300   0.014    0.021    0.340
## abs_diff_bscale               0.112 0.060 33.593   1.854   0.036   -0.011    0.235
## abs_sum_bscale                0.722 0.314 33.027   2.300   0.014    0.083    1.361
## abs_diff_rscale               0.133 0.062 33.635   2.157   0.019    0.008    0.259
## abs_sum_rscale                0.726 0.314 33.027   2.309   0.014    0.086    1.366
```

``` r
round(ddsc_mod2_log_GDP$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.004 -0.170  0.896  1.000  0.344
```

``` r
d_log_GDP<-ddsc_mod2_log_GDP$ddsc_sem_fit$data

ddsc_sem_log_GDP<-
  ddsc_sem(data=d_log_GDP,x = "log_gdp.z.cm",y1="means_y2",
           y2="means_y1",q_sesoi = .10)
round(ddsc_sem_log_GDP$results,3)
```

```
##                                     est    se       z pvalue ci.lower ci.upper
## r_xy1_y2                         -0.302 0.166  -1.820  0.069   -0.627    0.023
## r_xy1                            -0.410 0.159  -2.579  0.010   -0.721   -0.098
## r_xy2                            -0.319 0.165  -1.936  0.053   -0.643    0.004
## b_11                             -0.421 0.163  -2.579  0.010   -0.741   -0.101
## b_21                             -0.310 0.160  -1.936  0.053   -0.624    0.004
## b_10                             -0.603 0.161  -3.747  0.000   -0.918   -0.287
## b_20                              1.074 0.158   6.810  0.000    0.765    1.383
## res_cov_y1_y2                     0.777 0.199   3.910  0.000    0.388    1.167
## diff_b10_b20                     -1.676 0.060 -27.876  0.000   -1.794   -1.558
## diff_b11_b21                     -0.111 0.061  -1.820  0.069   -0.231    0.009
## diff_rxy1_rxy2                   -0.090 0.061  -1.469  0.142   -0.210    0.030
## q_b11_b21                        -0.129 0.074  -1.740  0.082   -0.273    0.016
## q_rxy1_rxy2                      -0.104 0.071  -1.467  0.143   -0.243    0.035
## cross_over_point                -15.081 8.304  -1.816  0.069  -31.356    1.193
## sum_b11_b21                      -0.731 0.318  -2.302  0.021   -1.354   -0.109
## main_effect                      -0.366 0.159  -2.302  0.021   -0.677   -0.054
## interaction_vs_main_effect       -0.254 0.167  -1.523  0.128   -0.582    0.073
## diff_abs_b11_abs_b21              0.111 0.061   1.820  0.069   -0.009    0.231
## abs_diff_b11_b21                  0.111 0.061   1.820  0.034   -0.009    0.231
## abs_sum_b11_b21                   0.731 0.318   2.302  0.011    0.109    1.354
## dadas                            -0.620 0.320  -1.936  0.974   -1.248    0.008
## q_r_equivalence                   0.004 0.071   0.057  0.523       NA       NA
## q_b_equivalence                   0.029 0.074   0.386  0.650       NA       NA
## cross_over_point_equivalence     15.081 8.304   1.816  0.965       NA       NA
## cross_over_point_minimal_effect  15.081 8.304   1.816  0.035       NA       NA
```

``` r
round(ddsc_sem_log_GDP$variance_test,3)
```

```
##             est    se      z pvalue ci.lower ci.upper
## cov_y1y2  0.904 0.231  3.921  0.000    0.452    1.356
## var_y1    1.026 0.253  4.062  0.000    0.531    1.521
## var_y2    0.914 0.225  4.062  0.000    0.473    1.355
## var_diff  0.112 0.124  0.906  0.365   -0.130    0.354
## var_ratio 1.123 0.140  8.030  0.000    0.849    1.397
## cor_y1y2  0.934 0.022 41.933  0.000    0.890    0.977
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

![](Analysis_files/figure-html/unnamed-chunk-31-1.png)<!-- -->

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


## mod3: fixed effect of time (Ess round)


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
## 1169930.2 1170006.3 -584958.1 1169916.2    392761 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.3689 -0.6146  0.0035  0.6025  9.9546 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr 
##  cntry    (Intercept) 0.057697 0.24020       
##           gndr.c      0.007705 0.08778  -0.16
##  Residual             1.012612 1.00629       
## Number of obs: 392768, groups:  cntry, 33
## 
## Fixed effects:
##               Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)  6.060e-02  4.186e-02  3.298e+01   1.448    0.157    
## gndr.c       4.197e-01  1.574e-02  3.276e+01  26.657   <2e-16 ***
## essround.c  -1.486e-02  6.152e-04  3.924e+05 -24.163   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##            (Intr) gndr.c
## gndr.c     -0.159       
## essround.c -0.002 -0.001
```

``` r
getFE(mod3,round=3)
```

```
##               Est.    SE         df       t     p     LL     UL
## (Intercept)  0.061 0.042     32.984   1.448 0.157 -0.025  0.146
## gndr.c       0.420 0.016     32.756  26.657 0.000  0.388  0.452
## essround.c  -0.015 0.001 392435.942 -24.163 0.000 -0.016 -0.014
```

``` r
getVC(mod3)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.24 0.06
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
## fixed           0.040494433
## slope variation 0.001710019
## mean variation  0.051900657
## sigma2          0.905894891
## 
## $R2s
##           total
## f   0.040494433
## v   0.001710019
## m   0.051900657
## fv  0.042204452
## fvm 0.094105109
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
## mod2    6 1170512 1170577 -585250   1170500                         
## mod3    7 1169930 1170006 -584958   1169916 583.43  1  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## mod4: random effect of time (Ess round)


``` r
mod4<-lmer(FM.z~gndr.c+essround.c+(gndr.c+essround.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod4)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: FM.z ~ gndr.c + essround.c + (gndr.c + essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1166638.1 1166747.0 -583309.1 1166618.1    392758 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.7081 -0.6130  0.0044  0.6009  9.9630 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr       
##  cntry    (Intercept) 0.096084 0.30997             
##           gndr.c      0.007881 0.08877  -0.22      
##           essround.c  0.004408 0.06640   0.36 -0.45
##  Residual             1.003663 1.00183             
## Number of obs: 392768, groups:  cntry, 33
## 
## Fixed effects:
##             Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)  0.04508    0.05414 30.55267   0.833    0.412    
## gndr.c       0.41898    0.01590 32.73362  26.343   <2e-16 ***
## essround.c   0.00115    0.01163 27.92329   0.099    0.922    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##            (Intr) gndr.c
## gndr.c     -0.209       
## essround.c  0.354 -0.433
```

``` r
getFE(mod4,round=3)
```

```
##              Est.    SE     df      t     p     LL    UL
## (Intercept) 0.045 0.054 30.553  0.833 0.412 -0.065 0.156
## gndr.c      0.419 0.016 32.734 26.343 0.000  0.387 0.451
## essround.c  0.001 0.012 27.923  0.099 0.922 -0.023 0.025
```

``` r
getVC(mod4)
```

```
##        grp        var1       var2 sdcor  vcov
## 1    cntry (Intercept)       <NA>  0.31  0.10
## 2    cntry      gndr.c       <NA>  0.09  0.01
## 3    cntry  essround.c       <NA>  0.07  0.00
## 4    cntry (Intercept)     gndr.c -0.22 -0.01
## 5    cntry (Intercept) essround.c  0.36  0.01
## 6    cntry      gndr.c essround.c -0.45  0.00
## 7 Residual        <NA>       <NA>  1.00  1.00
```

``` r
r2mlm(mod4,bargraph = F)
```

```
## $Decompositions
##                      total
## fixed           0.03696890
## slope variation 0.02978582
## mean variation  0.08161018
## sigma2          0.85163510
## 
## $R2s
##          total
## f   0.03696890
## v   0.02978582
## m   0.08161018
## fv  0.06675472
## fvm 0.14836490
```

``` r
anova(mod2,mod3,mod4)
```

```
## Data: diff_dat
## Models:
## mod2: FM.z ~ gndr.c + (gndr.c | cntry)
## mod3: FM.z ~ gndr.c + essround.c + (gndr.c | cntry)
## mod4: FM.z ~ gndr.c + essround.c + (gndr.c + essround.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L)   Chisq Df Pr(>Chisq)    
## mod2    6 1170512 1170577 -585250   1170500                          
## mod3    7 1169930 1170006 -584958   1169916  583.43  1  < 2.2e-16 ***
## mod4   10 1166638 1166747 -583309   1166618 3298.04  3  < 2.2e-16 ***
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
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: FM.z ~ gndr.c + essround.c + gndr.c:essround.c + (gndr.c + essround.c |      cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1166568.6 1166688.3 -583273.3 1166546.6    392757 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.7348 -0.6129  0.0043  0.6006  9.9849 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr       
##  cntry    (Intercept) 0.096450 0.31056             
##           gndr.c      0.008219 0.09066  -0.27      
##           essround.c  0.004441 0.06664   0.36 -0.41
##  Residual             1.003474 1.00174             
## Number of obs: 392768, groups:  cntry, 33
## 
## Fixed effects:
##                     Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)        4.492e-02  5.424e-02  3.044e+01   0.828    0.414    
## gndr.c             4.200e-01  1.623e-02  3.257e+01  25.886   <2e-16 ***
## essround.c         9.502e-04  1.167e-02  2.778e+01   0.081    0.936    
## gndr.c:essround.c -1.034e-02  1.221e-03  1.993e+05  -8.474   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c essrn.
## gndr.c      -0.264              
## essround.c   0.347 -0.398       
## gndr.c:ssr. -0.001 -0.007  0.003
```

``` r
getFE(mod5,round=3)
```

```
##                     Est.    SE         df      t     p     LL     UL
## (Intercept)        0.045 0.054     30.444  0.828 0.414 -0.066  0.156
## gndr.c             0.420 0.016     32.568 25.886 0.000  0.387  0.453
## essround.c         0.001 0.012     27.784  0.081 0.936 -0.023  0.025
## gndr.c:essround.c -0.010 0.001 199321.161 -8.474 0.000 -0.013 -0.008
```

``` r
getVC(mod5)
```

```
##        grp        var1       var2 sdcor  vcov
## 1    cntry (Intercept)       <NA>  0.31  0.10
## 2    cntry      gndr.c       <NA>  0.09  0.01
## 3    cntry  essround.c       <NA>  0.07  0.00
## 4    cntry (Intercept)     gndr.c -0.27 -0.01
## 5    cntry (Intercept) essround.c  0.36  0.01
## 6    cntry      gndr.c essround.c -0.41  0.00
## 7 Residual        <NA>       <NA>  1.00  1.00
```

``` r
r2mlm(mod5,bargraph = F)
```

```
## $Decompositions
##                      total
## fixed           0.03734108
## slope variation 0.03003836
## mean variation  0.08197219
## sigma2          0.85064838
## 
## $R2s
##          total
## f   0.03734108
## v   0.03003836
## m   0.08197219
## fv  0.06737944
## fvm 0.14935162
```

``` r
anova(mod4,mod5)
```

```
## Data: diff_dat
## Models:
## mod4: FM.z ~ gndr.c + essround.c + (gndr.c + essround.c | cntry)
## mod5: FM.z ~ gndr.c + essround.c + gndr.c:essround.c + (gndr.c + essround.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod4   10 1166638 1166747 -583309   1166618                         
## mod5   11 1166569 1166688 -583273   1166547 71.522  1  < 2.2e-16 ***
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
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: FM.z ~ gndr.c + essround.c + gndr.c:essround.c + (gndr.c + essround.c +  
##     gndr.c:essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1166538.7 1166701.9 -583254.4 1166508.7    392753 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.8013 -0.6128  0.0043  0.6005  9.9814 
## 
## Random effects:
##  Groups   Name              Variance  Std.Dev. Corr             
##  cntry    (Intercept)       0.0970165 0.31147                   
##           gndr.c            0.0091607 0.09571  -0.37            
##           essround.c        0.0044639 0.06681   0.35 -0.39      
##           gndr.c:essround.c 0.0001306 0.01143  -0.30  0.09 -0.53
##  Residual                   1.0032910 1.00164                   
## Number of obs: 392768, groups:  cntry, 33
## 
## Fixed effects:
##                    Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)        0.044507   0.054403 30.247265   0.818     0.42    
## gndr.c             0.422406   0.017172 30.231521  24.599  < 2e-16 ***
## essround.c         0.001038   0.011704 27.498383   0.089     0.93    
## gndr.c:essround.c -0.012529   0.002449 20.797385  -5.115  4.7e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c essrn.
## gndr.c      -0.358              
## essround.c   0.338 -0.379       
## gndr.c:ssr. -0.246  0.051 -0.431
```

``` r
getFE(mod6,round=3)
```

```
##                     Est.    SE     df      t    p     LL     UL
## (Intercept)        0.045 0.054 30.247  0.818 0.42 -0.067  0.156
## gndr.c             0.422 0.017 30.232 24.599 0.00  0.387  0.457
## essround.c         0.001 0.012 27.498  0.089 0.93 -0.023  0.025
## gndr.c:essround.c -0.013 0.002 20.797 -5.115 0.00 -0.018 -0.007
```

``` r
getVC(mod6)
```

```
##         grp              var1              var2 sdcor  vcov
## 1     cntry       (Intercept)              <NA>  0.31  0.10
## 2     cntry            gndr.c              <NA>  0.10  0.01
## 3     cntry        essround.c              <NA>  0.07  0.00
## 4     cntry gndr.c:essround.c              <NA>  0.01  0.00
## 5     cntry       (Intercept)            gndr.c -0.37 -0.01
## 6     cntry       (Intercept)        essround.c  0.35  0.01
## 7     cntry       (Intercept) gndr.c:essround.c -0.30  0.00
## 8     cntry            gndr.c        essround.c -0.39  0.00
## 9     cntry            gndr.c gndr.c:essround.c  0.09  0.00
## 10    cntry        essround.c gndr.c:essround.c -0.53  0.00
## 11 Residual              <NA>              <NA>  1.00  1.00
```

``` r
r2mlm(mod6,bargraph = F)
```

```
## $Decompositions
##                      total
## fixed           0.03778764
## slope variation 0.03073979
## mean variation  0.08254090
## sigma2          0.84893167
## 
## $R2s
##          total
## f   0.03778764
## v   0.03073979
## m   0.08254090
## fv  0.06852743
## fvm 0.15106833
```

``` r
anova(mod4,mod5,mod6)
```

```
## Data: diff_dat
## Models:
## mod4: FM.z ~ gndr.c + essround.c + (gndr.c + essround.c | cntry)
## mod5: FM.z ~ gndr.c + essround.c + gndr.c:essround.c + (gndr.c + essround.c | cntry)
## mod6: FM.z ~ gndr.c + essround.c + gndr.c:essround.c + (gndr.c + essround.c + gndr.c:essround.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod4   10 1166638 1166747 -583309   1166618                         
## mod5   11 1166569 1166688 -583273   1166547 71.522  1  < 2.2e-16 ***
## mod6   15 1166539 1166702 -583254   1166509 37.917  4  1.165e-07 ***
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
##         4.5 -0.134 0.0941 31.2  -0.3257   0.0580  -1.422  0.1649
##        -4.5 -0.200 0.0632 14.5  -0.3346  -0.0645  -3.160  0.0067
## 
## gndr.c =  0.5:
##  essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.232 0.0819 31.1   0.0652   0.3992   2.836  0.0080
##        -4.5  0.279 0.0617 16.9   0.1491   0.4094   4.529  0.0003
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_mod6,adjust="none",infer=c(T,T))
```

```
## gndr.c = -0.5:
##  contrast                         estimate    SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   0.0657 0.111 27.1   -0.161    0.292   0.595  0.5570
## 
## gndr.c =  0.5:
##  contrast                         estimate    SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)  -0.0470 0.101 26.9   -0.254    0.160  -0.465  0.6454
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
##    -0.5        4.5 -0.134 0.0941 31.2  -0.3257   0.0580  -1.422  0.1649
##     0.5        4.5  0.232 0.0819 31.1   0.0652   0.3992   2.836  0.0080
##    -0.5       -4.5 -0.200 0.0632 14.5  -0.3346  -0.0645  -3.160  0.0067
##     0.5       -4.5  0.279 0.0617 16.9   0.1491   0.4094   4.529  0.0003
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod6,adjust="none",infer=c(T,T))
```

```
##  contrast                                                 estimate     SE   df lower.CL upper.CL t.ratio
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.3660 0.0209 21.4   -0.409   -0.323 -17.538
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   0.0657 0.1110 27.1   -0.161    0.292   0.595
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.4131 0.1130 28.8   -0.644   -0.182  -3.656
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.4318 0.1000 25.7    0.226    0.638   4.313
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      -0.0470 0.1010 26.9   -0.254    0.160  -0.465
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.4788 0.0199 27.8   -0.520   -0.438 -24.026
##  p.value
##   <.0001
##   0.5570
##   0.0010
##   0.0002
##   0.6454
##   <.0001
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
##  diff_ESS10    0.366 0.0209 21.4    0.323    0.409  17.538  <.0001
##  diff_ESS1     0.479 0.0199 27.8    0.438    0.520  24.026  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# Test for whether and how much the differences have changed
pairs(diff_mod6,infer=c(T,T))
```

```
##  contrast               estimate    SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10 - diff_ESS1   -0.113 0.022 20.8   -0.159  -0.0669  -5.115  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
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

p_time_trends
```

![](Analysis_files/figure-html/unnamed-chunk-37-1.png)<!-- -->

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
## [1] -0.8757598  0.8407629
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
      geom_smooth(method = "lm",formula = "y ~ x", se = FALSE) +
      geom_point(size=8) +
      geom_point(aes(x=year,y=obs_mean_wt),size=8,shape = 1,alpha=.50)+
      geom_errorbar(aes(ymin=obs_mean_wt_LL, ymax=obs_mean_wt_UL),alpha = .50)+
      scale_color_manual(values = my_colors) +
      geom_flag(aes(country=tolower(ctry)))+
      scale_y_continuous(limits = c(-1.1, 1.1)) +
      scale_x_continuous(limits = c(2001, 2021), breaks = seq(2002, 2020, 2)) +
      ggtitle(paste("Country:", ctry))+
    ylab("Mean-level of value male-typicality")+
      theme(legend.title=element_blank())
  )
}

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
  scale_y_continuous(limits = c(-1.1, 1.1)) +
  scale_x_continuous(limits = c(2001, 2021), breaks = seq(2002, 2020, 2)) +
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

![](Analysis_files/figure-html/unnamed-chunk-39-1.png)<!-- -->

``` r
png(filename = 
      "../results/country_time_trend_facets.png",
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
## [1] 34.17356
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
##    gndr.c change_per_18_years gndr_change_per_18_year men_change_per_18_years women_change_per_18_years
## 1    0.52               -0.24                   -0.09                   -0.29                     -0.19
## 2    0.39               -0.07                   -0.02                   -0.08                     -0.06
## 3    0.46                0.67                   -0.14                    0.60                      0.74
## 4    0.44               -0.03                    0.01                   -0.03                     -0.03
## 5    0.31               -0.63                   -0.03                   -0.65                     -0.62
## 6    0.59                0.36                   -0.26                    0.23                      0.49
## 7    0.50               -0.51                   -0.11                   -0.57                     -0.46
## 8    0.48               -0.21                   -0.04                   -0.23                     -0.18
## 9    0.53               -0.17                   -0.05                   -0.20                     -0.14
## 10   0.36               -0.55                   -0.03                   -0.57                     -0.54
## 11   0.56               -0.38                   -0.13                   -0.44                     -0.31
## 12   0.49               -0.33                   -0.01                   -0.33                     -0.32
## 13   0.48               -0.39                   -0.03                   -0.40                     -0.37
## 14   0.31                0.10                    0.01                    0.10                      0.09
## 15   0.44               -0.16                   -0.14                   -0.23                     -0.09
## 16   0.33                0.42                   -0.22                    0.31                      0.53
## 17   0.39               -0.15                   -0.06                   -0.18                     -0.12
## 18   0.30               -0.02                   -0.08                   -0.06                      0.02
## 19   0.49               -0.27                   -0.16                   -0.35                     -0.20
## 20   0.35                0.88                   -0.19                    0.78                      0.97
## 21   0.45               -0.05                   -0.06                   -0.08                     -0.02
## 22   0.49               -0.95                   -0.09                   -1.00                     -0.91
## 23   0.49                1.44                   -0.22                    1.33                      1.55
## 24   0.47               -0.26                   -0.09                   -0.30                     -0.21
## 25   0.44               -0.29                   -0.20                   -0.39                     -0.19
## 26   0.47               -0.10                   -0.08                   -0.14                     -0.06
## 27   0.29               -0.35                   -0.04                   -0.37                     -0.32
## 28   0.31                0.33                   -0.32                    0.17                      0.49
## 29   0.49               -0.41                   -0.12                   -0.47                     -0.35
## 30   0.41               -0.35                   -0.08                   -0.39                     -0.31
## 31   0.46                0.30                   -0.21                    0.20                      0.41
## 32   0.16                2.14                   -0.31                    1.99                      2.30
## 33   0.31                0.53                   -0.12                    0.47                      0.59
##    cntry n_unique_essround
## 1     AT                 6
## 2     BE                10
## 3     BG                 6
## 4     CH                10
## 5     CY                 5
## 6     CZ                 9
## 7     DE                 9
## 8     DK                 8
## 9     EE                 9
## 10    ES                 9
## 11    FI                10
## 12    FR                10
## 13    GB                10
## 14    GR                 5
## 15    HR                 4
## 16    HU                10
## 17    IE                10
## 18    IL                 6
## 19    IS                 5
## 20    IT                 4
## 21    LT                 6
## 22    LV                 2
## 23    ME                 2
## 24    NL                10
## 25    NO                10
## 26    PL                 9
## 27    PT                10
## 28    RU                 5
## 29    SE                 9
## 30    SI                10
## 31    SK                 7
## 32    TR                 2
## 33    UA                 5
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
## 1     CY               -0.63
## 2     ES               -0.55
## 3     DE               -0.51
## 4     SE               -0.41
## 5     GB               -0.39
## 6     FI               -0.38
## 7     PT               -0.35
## 8     SI               -0.35
## 9     FR               -0.33
## 10    NO               -0.29
## 11    IS               -0.27
## 12    NL               -0.26
## 13    AT               -0.24
## 14    DK               -0.21
## 15    EE               -0.17
## 16    IE               -0.15
## 17    PL               -0.10
## 18    BE               -0.07
## 19    LT               -0.05
## 20    CH               -0.03
## 21    IL               -0.02
## 22    GR                0.10
## 23    SK                0.30
## 24    RU                0.33
## 25    CZ                0.36
## 26    HU                0.42
## 27    UA                0.53
## 28    BG                0.67
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
## 1     RU                   -0.32
## 2     CZ                   -0.26
## 3     HU                   -0.22
## 4     SK                   -0.21
## 5     NO                   -0.20
## 6     IS                   -0.16
## 7     BG                   -0.14
## 8     FI                   -0.13
## 9     SE                   -0.12
## 10    UA                   -0.12
## 11    DE                   -0.11
## 12    AT                   -0.09
## 13    NL                   -0.09
## 14    IL                   -0.08
## 15    PL                   -0.08
## 16    SI                   -0.08
## 17    IE                   -0.06
## 18    LT                   -0.06
## 19    EE                   -0.05
## 20    DK                   -0.04
## 21    PT                   -0.04
## 22    CY                   -0.03
## 23    ES                   -0.03
## 24    GB                   -0.03
## 25    BE                   -0.02
## 26    FR                   -0.01
## 27    CH                    0.01
## 28    GR                    0.01
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
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: FM.z ~ gndr.c + essround.c + gndr.c:essround.c + gei.z.cm:gndr.c +  
##     gei.z.cm:essround.c + gei.z.cm:gndr.c:essround.c + (gndr.c +  
##     essround.c + gndr.c:essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##   1136150   1136346   -568057   1136114    384360 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.8571 -0.6139  0.0047  0.6023 10.0312 
## 
## Random effects:
##  Groups   Name              Variance  Std.Dev. Corr             
##  cntry    (Intercept)       0.1029883 0.32092                   
##           gndr.c            0.0066946 0.08182  -0.17            
##           essround.c        0.0033699 0.05805  -0.02 -0.07      
##           gndr.c:essround.c 0.0001041 0.01020  -0.10 -0.21 -0.33
##  Residual                   0.9937232 0.99686                   
## Number of obs: 384378, groups:  cntry, 32
## 
## Fixed effects:
##                              Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)                 0.0312916  0.0569342 28.3373494   0.550  0.58689    
## gndr.c                      0.4252873  0.0150716 28.9067103  28.218  < 2e-16 ***
## essround.c                  0.0002312  0.0103548 22.2727603   0.022  0.98238    
## gndr.c:essround.c          -0.0131417  0.0023423 22.1222075  -5.611 1.19e-05 ***
## gndr.c:gei.z.cm             0.0459392  0.0154397 33.7886979   2.975  0.00537 ** 
## essround.c:gei.z.cm        -0.0396317  0.0105926 23.2800845  -3.741  0.00105 ** 
## gndr.c:essround.c:gei.z.cm  0.0060235  0.0027064 25.6296013   2.226  0.03505 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c essrn. gnd.:. gn.:.. es.:..
## gndr.c      -0.159                                   
## essround.c  -0.027 -0.065                            
## gndr.c:ssr. -0.078 -0.172 -0.252                     
## gndr.c:g.z.  0.000 -0.027  0.000 -0.044              
## essrnd.c:.. -0.001  0.000 -0.012  0.001 -0.068       
## gndr.c:.:..  0.000 -0.025  0.000 -0.188 -0.036 -0.224
```

``` r
getFE(mod6_GEI,round=2,p.round = 5)
```

```
##                             Est.   SE    df     t       p    LL    UL
## (Intercept)                 0.03 0.06 28.34  0.55 0.58689 -0.09  0.15
## gndr.c                      0.43 0.02 28.91 28.22 0.00000  0.39  0.46
## essround.c                  0.00 0.01 22.27  0.02 0.98238 -0.02  0.02
## gndr.c:essround.c          -0.01 0.00 22.12 -5.61 0.00001 -0.02 -0.01
## gndr.c:gei.z.cm             0.05 0.02 33.79  2.98 0.00537  0.01  0.08
## essround.c:gei.z.cm        -0.04 0.01 23.28 -3.74 0.00105 -0.06 -0.02
## gndr.c:essround.c:gei.z.cm  0.01 0.00 25.63  2.23 0.03505  0.00  0.01
```

``` r
getVC(mod6_GEI)
```

```
##         grp              var1              var2 sdcor vcov
## 1     cntry       (Intercept)              <NA>  0.32 0.10
## 2     cntry            gndr.c              <NA>  0.08 0.01
## 3     cntry        essround.c              <NA>  0.06 0.00
## 4     cntry gndr.c:essround.c              <NA>  0.01 0.00
## 5     cntry       (Intercept)            gndr.c -0.17 0.00
## 6     cntry       (Intercept)        essround.c -0.02 0.00
## 7     cntry       (Intercept) gndr.c:essround.c -0.10 0.00
## 8     cntry            gndr.c        essround.c -0.07 0.00
## 9     cntry            gndr.c gndr.c:essround.c -0.21 0.00
## 10    cntry        essround.c gndr.c:essround.c -0.33 0.00
## 11 Residual              <NA>              <NA>  1.00 0.99
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_GEI,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] 24.5079
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_GEI,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] 20.28273
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
##  essround.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.2107 0.0872 38.5   0.0342   0.3871   2.416  0.0205
##        -4.5 -0.1481 0.0887 19.3  -0.3336   0.0375  -1.669  0.1113
## 
## gei.z.cm =  0:
##  essround.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.0323 0.0726 26.9  -0.1167   0.1813   0.445  0.6596
##        -4.5  0.0303 0.0745 14.4  -0.1292   0.1897   0.406  0.6908
## 
## gei.z.cm =  1:
##  essround.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5 -0.1460 0.0865 37.1  -0.3213   0.0292  -1.688  0.0998
##        -4.5  0.2086 0.0882 18.9   0.0239   0.3933   2.365  0.0289
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
##  essround.c4.5 - (essround.c-4.5)  0.35877 0.1340 23.3   0.0815    0.636   2.675  0.0135
## 
## gei.z.cm =  0:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)  0.00208 0.0932 22.3  -0.1911    0.195   0.022  0.9824
## 
## gei.z.cm =  1:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5) -0.35460 0.1320 22.3  -0.6292   -0.080  -2.676  0.0137
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
##  gndr.c essround.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5        4.5  0.0641 0.0912 38.3  -0.1205  0.24875   0.703  0.4864
##     0.5        4.5  0.3572 0.0852 38.5   0.1849  0.52959   4.194  0.0002
##    -0.5       -4.5 -0.3809 0.0913 17.4  -0.5731 -0.18863  -4.173  0.0006
##     0.5       -4.5  0.0847 0.0885 20.4  -0.0997  0.26911   0.957  0.3498
## 
## gei.z.cm =  0:
##  gndr.c essround.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5        4.5 -0.1507 0.0755 26.7  -0.3057  0.00417  -1.998  0.0560
##     0.5        4.5  0.2154 0.0706 26.7   0.0704  0.36041   3.050  0.0051
##    -0.5       -4.5 -0.2120 0.0763 13.2  -0.3765 -0.04738  -2.778  0.0155
##     0.5       -4.5  0.2725 0.0741 15.1   0.1147  0.43023   3.679  0.0022
## 
## gei.z.cm =  1:
##  gndr.c essround.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5        4.5 -0.3656 0.0903 36.7  -0.5487 -0.18254  -4.048  0.0003
##     0.5        4.5  0.0736 0.0842 36.7  -0.0971  0.24429   0.874  0.3879
##    -0.5       -4.5 -0.0430 0.0906 17.0  -0.2343  0.14820  -0.475  0.6410
##     0.5       -4.5  0.4602 0.0879 19.9   0.2769  0.64357   5.238  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod6_GEI,adjust="none",infer=c(T,T))
```

```
## gei.z.cm = -1:
##  contrast                                                 estimate     SE   df  lower.CL upper.CL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.2931 0.0272 35.3 -0.348342  -0.2379
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   0.4450 0.1390 22.4  0.157095   0.7329
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.0206 0.1370 24.6 -0.303601   0.2624
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.7381 0.1340 21.9  0.459082   1.0171
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)       0.2725 0.1310 23.2  0.000732   0.5443
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.4656 0.0288 29.0 -0.524527  -0.4067
##  t.ratio p.value
##  -10.769  <.0001
##    3.202  0.0040
##   -0.150  0.8821
##    5.488  <.0001
##    2.073  0.0494
##  -16.157  <.0001
## 
## gei.z.cm =  0:
##  contrast                                                 estimate     SE   df  lower.CL upper.CL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.3661 0.0168 26.6 -0.400722  -0.3316
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   0.0612 0.0964 21.5 -0.138971   0.2614
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.4232 0.0954 23.8 -0.620100  -0.2263
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.4274 0.0934 20.5  0.232795   0.6219
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      -0.0571 0.0911 21.9 -0.246032   0.1319
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.4844 0.0198 25.0 -0.525258  -0.4436
##  t.ratio p.value
##  -21.747  <.0001
##    0.635  0.5321
##   -4.438  0.0002
##    4.574  0.0002
##   -0.626  0.5376
##  -24.435  <.0001
## 
## gei.z.cm =  1:
##  contrast                                                 estimate     SE   df  lower.CL upper.CL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.4392 0.0239 27.7 -0.488235  -0.3902
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  -0.3226 0.1370 21.4 -0.607233  -0.0379
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.8258 0.1360 23.6 -1.105991  -0.5457
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.1166 0.1330 20.9 -0.159530   0.3928
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      -0.3866 0.1290 22.0 -0.655036  -0.1182
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.5033 0.0275 24.7 -0.559922  -0.4466
##  t.ratio p.value
##  -18.353  <.0001
##   -2.354  0.0282
##   -6.090  <.0001
##    0.878  0.3897
##   -2.987  0.0068
##  -18.304  <.0001
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
##  diff_ESS10    0.293 0.0272 35.3    0.238    0.348  10.769  <.0001
##  diff_ESS1     0.466 0.0288 29.0    0.407    0.525  16.157  <.0001
## 
## gei.z.cm =  0:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10    0.366 0.0168 26.6    0.332    0.401  21.747  <.0001
##  diff_ESS1     0.484 0.0198 25.0    0.444    0.525  24.435  <.0001
## 
## gei.z.cm =  1:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10    0.439 0.0239 27.7    0.390    0.488  18.353  <.0001
##  diff_ESS1     0.503 0.0275 24.7    0.447    0.560  18.304  <.0001
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
##  diff_ESS10 - diff_ESS1  -0.1725 0.0351 29.0   -0.244 -0.10074  -4.917  <.0001
## 
## gei.z.cm =  0:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10 - diff_ESS1  -0.1183 0.0211 22.1   -0.162 -0.07457  -5.611  <.0001
## 
## gei.z.cm =  1:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10 - diff_ESS1  -0.0641 0.0291 19.0   -0.125 -0.00323  -2.204  0.0401
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
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: FM.z ~ gndr.c + essround.c + gndr.c:essround.c + gggi.z.cm:gndr.c +  
##     gggi.z.cm:essround.c + gggi.z.cm:gndr.c:essround.c + (gndr.c +  
##     essround.c + gndr.c:essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##  833017.8  833207.6 -416490.9  832981.8    280628 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.8248 -0.6121  0.0063  0.6025  8.1720 
## 
## Random effects:
##  Groups   Name              Variance  Std.Dev. Corr             
##  cntry    (Intercept)       0.0766667 0.27689                   
##           gndr.c            0.0057648 0.07593  -0.04            
##           essround.c        0.0018252 0.04272   0.30  0.01      
##           gndr.c:essround.c 0.0001984 0.01409  -0.17 -0.10 -0.60
##  Residual                   0.9997799 0.99989                   
## Number of obs: 280646, groups:  cntry, 33
## 
## Fixed effects:
##                              Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)                  0.066560   0.048521 32.339522   1.372  0.17958    
## gndr.c                       0.417723   0.014114 32.854653  29.597  < 2e-16 ***
## essround.c                  -0.008746   0.007843 28.411539  -1.115  0.27413    
## gndr.c:essround.c           -0.007067   0.003246 23.488523  -2.177  0.03974 *  
## gndr.c:gggi.z.cm             0.050850   0.014756 36.565994   3.446  0.00144 ** 
## essround.c:gggi.z.cm        -0.014351   0.008157 29.749812  -1.759  0.08881 .  
## gndr.c:essround.c:gggi.z.cm -0.002337   0.003611 26.654044  -0.647  0.52301    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c essrn. gnd.:. gn.:.. es.:..
## gndr.c      -0.044                                   
## essround.c   0.285  0.006                            
## gndr.c:ssr. -0.128 -0.108 -0.458                     
## gndr.c:gg..  0.003 -0.013  0.008 -0.024              
## essrnd.c:.. -0.011  0.007 -0.078  0.028  0.001       
## gndr.c:.:..  0.005 -0.027  0.025 -0.083 -0.107 -0.421
```

``` r
getFE(mod6_GGGI,round=2,p.round = 5)
```

```
##                              Est.   SE    df     t       p    LL   UL
## (Intercept)                  0.07 0.05 32.34  1.37 0.17958 -0.03 0.17
## gndr.c                       0.42 0.01 32.85 29.60 0.00000  0.39 0.45
## essround.c                  -0.01 0.01 28.41 -1.12 0.27413 -0.02 0.01
## gndr.c:essround.c           -0.01 0.00 23.49 -2.18 0.03974 -0.01 0.00
## gndr.c:gggi.z.cm             0.05 0.01 36.57  3.45 0.00144  0.02 0.08
## essround.c:gggi.z.cm        -0.01 0.01 29.75 -1.76 0.08881 -0.03 0.00
## gndr.c:essround.c:gggi.z.cm  0.00 0.00 26.65 -0.65 0.52301 -0.01 0.01
```

``` r
getVC(mod6_GGGI)
```

```
##         grp              var1              var2 sdcor vcov
## 1     cntry       (Intercept)              <NA>  0.28 0.08
## 2     cntry            gndr.c              <NA>  0.08 0.01
## 3     cntry        essround.c              <NA>  0.04 0.00
## 4     cntry gndr.c:essround.c              <NA>  0.01 0.00
## 5     cntry       (Intercept)            gndr.c -0.04 0.00
## 6     cntry       (Intercept)        essround.c  0.30 0.00
## 7     cntry       (Intercept) gndr.c:essround.c -0.17 0.00
## 8     cntry            gndr.c        essround.c  0.01 0.00
## 9     cntry            gndr.c gndr.c:essround.c -0.10 0.00
## 10    cntry        essround.c gndr.c:essround.c -0.60 0.00
## 11 Residual              <NA>              <NA>  1.00 1.00
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_GGGI,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] 59.11212
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_GGGI,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] -51.88158
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
##  essround.c  emmean     SE   df  lower.CL upper.CL t.ratio p.value
##         4.5  0.0918 0.0785 39.0 -0.067021    0.251   1.169  0.2495
##        -4.5  0.0413 0.0643 36.9 -0.088947    0.172   0.643  0.5242
## 
## gggi.z.cm =  0:
##  essround.c  emmean     SE   df  lower.CL upper.CL t.ratio p.value
##         4.5  0.0272 0.0676 27.1 -0.111570    0.166   0.402  0.6907
##        -4.5  0.1059 0.0512 26.7  0.000759    0.211   2.068  0.0485
## 
## gggi.z.cm =  1:
##  essround.c  emmean     SE   df  lower.CL upper.CL t.ratio p.value
##         4.5 -0.0374 0.0754 36.5 -0.190178    0.115  -0.496  0.6230
##        -4.5  0.1705 0.0617 35.6  0.045274    0.296   2.762  0.0090
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
##  essround.c4.5 - (essround.c-4.5)   0.0504 0.1060 31.0   -0.165  0.26607   0.477  0.6366
## 
## gggi.z.cm =  0:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)  -0.0787 0.0706 28.4   -0.223  0.06578  -1.115  0.2741
## 
## gggi.z.cm =  1:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)  -0.2079 0.0978 30.1   -0.408 -0.00815  -2.125  0.0419
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
##    -0.5        4.5 -0.081 0.0838 38.7  -0.2506  0.08860  -0.966  0.3399
##     0.5        4.5  0.265 0.0757 38.6   0.1114  0.41778   3.494  0.0012
##    -0.5       -4.5 -0.153 0.0696 35.6  -0.2939 -0.01163  -2.196  0.0347
##     0.5       -4.5  0.235 0.0627 37.4   0.1084  0.36248   3.753  0.0006
## 
## gggi.z.cm =  0:
##  gndr.c essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5        4.5 -0.166 0.0709 26.6  -0.3113 -0.02025  -2.339  0.0271
##     0.5        4.5  0.220 0.0657 27.1   0.0854  0.35491   3.352  0.0024
##    -0.5       -4.5 -0.119 0.0540 26.2  -0.2298 -0.00786  -2.200  0.0368
##     0.5       -4.5  0.331 0.0506 26.9   0.2269  0.43450   6.536  <.0001
## 
## gggi.z.cm =  1:
##  gndr.c essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5        4.5 -0.251 0.0802 36.4  -0.4130 -0.08798  -3.125  0.0035
##     0.5        4.5  0.176 0.0728 36.0   0.0281  0.32345   2.413  0.0210
##    -0.5       -4.5 -0.085 0.0668 34.4  -0.2206  0.05074  -1.272  0.2120
##     0.5       -4.5  0.426 0.0603 36.1   0.3036  0.54831   7.059  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod6_GGGI,adjust="none",infer=c(T,T))
```

```
## gggi.z.cm = -1:
##  contrast                                                 estimate     SE   df lower.CL upper.CL t.ratio
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.3456 0.0294 30.4  -0.4057  -0.2855 -11.742
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   0.0717 0.1170 30.7  -0.1676   0.3111   0.611
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.3164 0.1080 31.2  -0.5362  -0.0967  -2.936
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.4173 0.1080 31.4   0.1979   0.6367   3.877
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)       0.0292 0.0981 30.6  -0.1710   0.2294   0.297
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.3882 0.0318 33.0  -0.4529  -0.3234 -12.202
##  p.value
##   <.0001
##   0.5454
##   0.0062
##   0.0005
##   0.7683
##   <.0001
## 
## gggi.z.cm =  0:
##  contrast                                                 estimate     SE   df lower.CL upper.CL t.ratio
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.3859 0.0192 28.3  -0.4252  -0.3466 -20.113
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  -0.0469 0.0784 28.0  -0.2074   0.1136  -0.599
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.4964 0.0719 28.6  -0.6436  -0.3493  -6.905
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.3390 0.0721 28.5   0.1915   0.4865   4.704
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      -0.1105 0.0652 27.9  -0.2441   0.0231  -1.695
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.4495 0.0214 29.3  -0.4932  -0.4058 -21.031
##  p.value
##   <.0001
##   0.5542
##   <.0001
##   0.0001
##   0.1012
##   <.0001
## 
## gggi.z.cm =  1:
##  contrast                                                 estimate     SE   df lower.CL upper.CL t.ratio
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.4263 0.0270 29.0  -0.4816  -0.3710 -15.764
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  -0.1655 0.1090 29.6  -0.3878   0.0567  -1.522
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.6764 0.0997 30.4  -0.8799  -0.4730  -6.787
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.2607 0.1000 30.1   0.0563   0.4651   2.604
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      -0.2502 0.0904 29.7  -0.4349  -0.0655  -2.767
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.5109 0.0311 29.9  -0.5744  -0.4473 -16.421
##  p.value
##   <.0001
##   0.1387
##   <.0001
##   0.0142
##   0.0096
##   <.0001
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
##  diff_ESS10    0.346 0.0294 30.4    0.286    0.406  11.742  <.0001
##  diff_ESS1     0.388 0.0318 33.0    0.323    0.453  12.202  <.0001
## 
## gggi.z.cm =  0:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10    0.386 0.0192 28.3    0.347    0.425  20.113  <.0001
##  diff_ESS1     0.450 0.0214 29.3    0.406    0.493  21.031  <.0001
## 
## gggi.z.cm =  1:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10    0.426 0.0270 29.0    0.371    0.482  15.764  <.0001
##  diff_ESS1     0.511 0.0311 29.9    0.447    0.574  16.421  <.0001
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
##  diff_ESS10 - diff_ESS1  -0.0426 0.0455 26.8   -0.136  0.05076  -0.936  0.3575
## 
## gggi.z.cm =  0:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10 - diff_ESS1  -0.0636 0.0292 23.5   -0.124 -0.00325  -2.177  0.0397
## 
## gggi.z.cm =  1:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10 - diff_ESS1  -0.0846 0.0419 23.9   -0.171  0.00177  -2.022  0.0545
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
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: FM.z ~ gndr.c + essround.c + gndr.c:essround.c + gdi.z.cm:gndr.c +  
##     gdi.z.cm:essround.c + gdi.z.cm:gndr.c:essround.c + (gndr.c +  
##     essround.c + gndr.c:essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1166534.4 1166730.3 -583249.2 1166498.4    392750 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.8036 -0.6128  0.0044  0.6005  9.9884 
## 
## Random effects:
##  Groups   Name              Variance Std.Dev. Corr             
##  cntry    (Intercept)       0.095831 0.30957                   
##           gndr.c            0.008373 0.09151  -0.38            
##           essround.c        0.003613 0.06011   0.42 -0.32      
##           gndr.c:essround.c 0.000125 0.01118  -0.31  0.07 -0.58
##  Residual                   1.003296 1.00165                   
## Number of obs: 392768, groups:  cntry, 33
## 
## Fixed effects:
##                              Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)                 0.0461379  0.0540568 30.6910423   0.854  0.39999    
## gndr.c                      0.4229073  0.0164573 30.6099945  25.697  < 2e-16 ***
## essround.c                  0.0007234  0.0105407 28.4827314   0.069  0.94577    
## gndr.c:essround.c          -0.0126719  0.0024074 20.5225962  -5.264 3.47e-05 ***
## gndr.c:gdi.z.cm             0.0267235  0.0159426 36.4517256   1.676  0.10225    
## essround.c:gdi.z.cm        -0.0287659  0.0098124 20.1032721  -2.932  0.00822 ** 
## gndr.c:essround.c:gdi.z.cm -0.0006865  0.0028421 27.8697280  -0.242  0.81089    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c essrn. gnd.:. gn.:.. es.:..
## gndr.c      -0.364                                   
## essround.c   0.412 -0.306                            
## gndr.c:ssr. -0.245  0.034 -0.467                     
## gndr.c:gd..  0.000 -0.012  0.001 -0.007              
## essrnd.c:..  0.004 -0.001 -0.012  0.007 -0.176       
## gndr.c:.:.. -0.001 -0.010  0.006 -0.006  0.035 -0.342
```

``` r
getFE(mod6_GDI,round=2,p.round = 5)
```

```
##                             Est.   SE    df     t       p    LL    UL
## (Intercept)                 0.05 0.05 30.69  0.85 0.39999 -0.06  0.16
## gndr.c                      0.42 0.02 30.61 25.70 0.00000  0.39  0.46
## essround.c                  0.00 0.01 28.48  0.07 0.94577 -0.02  0.02
## gndr.c:essround.c          -0.01 0.00 20.52 -5.26 0.00003 -0.02 -0.01
## gndr.c:gdi.z.cm             0.03 0.02 36.45  1.68 0.10225 -0.01  0.06
## essround.c:gdi.z.cm        -0.03 0.01 20.10 -2.93 0.00822 -0.05 -0.01
## gndr.c:essround.c:gdi.z.cm  0.00 0.00 27.87 -0.24 0.81089 -0.01  0.01
```

``` r
getVC(mod6_GDI)
```

```
##         grp              var1              var2 sdcor  vcov
## 1     cntry       (Intercept)              <NA>  0.31  0.10
## 2     cntry            gndr.c              <NA>  0.09  0.01
## 3     cntry        essround.c              <NA>  0.06  0.00
## 4     cntry gndr.c:essround.c              <NA>  0.01  0.00
## 5     cntry       (Intercept)            gndr.c -0.38 -0.01
## 6     cntry       (Intercept)        essround.c  0.42  0.01
## 7     cntry       (Intercept) gndr.c:essround.c -0.31  0.00
## 8     cntry            gndr.c        essround.c -0.32  0.00
## 9     cntry            gndr.c gndr.c:essround.c  0.07  0.00
## 10    cntry        essround.c gndr.c:essround.c -0.58  0.00
## 11 Residual              <NA>              <NA>  1.00  1.00
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_GDI,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] 19.06461
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_GDI,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] 4.280517
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
##  essround.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.1788 0.0963 41.4  -0.0155   0.3732   1.858  0.0703
##        -4.5 -0.0866 0.0713 19.2  -0.2356   0.0625  -1.214  0.2393
## 
## gdi.z.cm =  0:
##  essround.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.0494 0.0853 31.9  -0.1245   0.2233   0.579  0.5668
##        -4.5  0.0429 0.0553 15.5  -0.0747   0.1604   0.775  0.4499
## 
## gdi.z.cm =  1:
##  essround.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5 -0.0801 0.0959 40.7  -0.2738   0.1137  -0.835  0.4088
##        -4.5  0.1723 0.0703 18.6   0.0250   0.3196   2.452  0.0243
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
##  essround.c4.5 - (essround.c-4.5)  0.26540 0.1300 26.9 -0.00221   0.5330   2.035  0.0518
## 
## gdi.z.cm =  0:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)  0.00651 0.0949 28.5 -0.18767   0.2007   0.069  0.9458
## 
## gdi.z.cm =  1:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5) -0.25238 0.1290 26.1 -0.51711   0.0123  -1.959  0.0608
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
##  gndr.c essround.c   emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5        4.5  0.00771 0.1040 41.4 -0.20160   0.2170   0.074  0.9410
##     0.5        4.5  0.34996 0.0906 41.0  0.16699   0.5329   3.863  0.0004
##    -0.5       -4.5 -0.31162 0.0744 17.6 -0.46820  -0.1550  -4.189  0.0006
##     0.5       -4.5  0.13850 0.0708 20.8 -0.00887   0.2859   1.955  0.0641
## 
## gdi.z.cm =  0:
##  gndr.c essround.c   emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5        4.5 -0.13355 0.0913 31.6 -0.31955   0.0525  -1.463  0.1533
##     0.5        4.5  0.23234 0.0802 31.5  0.06879   0.3959   2.895  0.0068
##    -0.5       -4.5 -0.19708 0.0578 14.4 -0.32064  -0.0735  -3.412  0.0041
##     0.5       -4.5  0.28285 0.0545 16.7  0.16770   0.3980   5.189  0.0001
## 
## gdi.z.cm =  1:
##  gndr.c essround.c   emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5        4.5 -0.27481 0.1030 40.7 -0.48340  -0.0662  -2.661  0.0111
##     0.5        4.5  0.11471 0.0903 40.3 -0.06770   0.2971   1.271  0.2111
##    -0.5       -4.5 -0.08254 0.0733 17.0 -0.23722   0.0721  -1.126  0.2759
##     0.5       -4.5  0.42720 0.0699 20.3  0.28147   0.5729   6.110  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod6_GDI,adjust="none",infer=c(T,T))
```

```
## gdi.z.cm = -1:
##  contrast                                                 estimate     SE   df lower.CL upper.CL t.ratio
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.3422 0.0291 31.1  -0.4016  -0.2829 -11.758
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   0.3193 0.1380 26.5   0.0359   0.6028   2.314
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.1308 0.1380 29.4  -0.4125   0.1509  -0.949
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.6616 0.1270 24.6   0.4001   0.9230   5.216
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)       0.2115 0.1250 26.8  -0.0443   0.4672   1.697
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.4501 0.0279 28.4  -0.5073  -0.3929 -16.112
##  p.value
##   <.0001
##   0.0287
##   0.3503
##   <.0001
##   0.1013
##   <.0001
## 
## gdi.z.cm =  0:
##  contrast                                                 estimate     SE   df lower.CL upper.CL t.ratio
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.3659 0.0200 23.1  -0.4073  -0.3245 -18.283
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   0.0635 0.1000 28.0  -0.1421   0.2691   0.633
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.4164 0.1010 29.8  -0.6230  -0.2098  -4.118
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.4294 0.0912 26.6   0.2422   0.6167   4.709
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      -0.0505 0.0903 28.0  -0.2355   0.1345  -0.559
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.4799 0.0194 25.7  -0.5198  -0.4401 -24.753
##  p.value
##   <.0001
##   0.5319
##   0.0003
##   0.0001
##   0.5804
##   <.0001
## 
## gdi.z.cm =  1:
##  contrast                                                 estimate     SE   df lower.CL upper.CL t.ratio
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.3895 0.0286 31.8  -0.4478  -0.3312 -13.619
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  -0.1923 0.1360 25.7  -0.4726   0.0880  -1.411
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.7020 0.1360 28.5  -0.9809  -0.4232  -5.153
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.1972 0.1250 23.8  -0.0611   0.4556   1.576
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      -0.3125 0.1230 26.1  -0.5656  -0.0594  -2.538
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.5097 0.0279 30.8  -0.5666  -0.4528 -18.276
##  p.value
##   <.0001
##   0.1703
##   <.0001
##   0.1282
##   0.0175
##   <.0001
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
##  diff_ESS10    0.342 0.0291 31.1    0.283    0.402  11.758  <.0001
##  diff_ESS1     0.450 0.0279 28.4    0.393    0.507  16.112  <.0001
## 
## gdi.z.cm =  0:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10    0.366 0.0200 23.1    0.325    0.407  18.283  <.0001
##  diff_ESS1     0.480 0.0194 25.7    0.440    0.520  24.753  <.0001
## 
## gdi.z.cm =  1:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10    0.390 0.0286 31.8    0.331    0.448  13.619  <.0001
##  diff_ESS1     0.510 0.0279 30.8    0.453    0.567  18.276  <.0001
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
##  diff_ESS10 - diff_ESS1   -0.108 0.0336 23.3   -0.177  -0.0384  -3.209  0.0039
## 
## gdi.z.cm =  0:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10 - diff_ESS1   -0.114 0.0217 20.5   -0.159  -0.0689  -5.264  <.0001
## 
## gdi.z.cm =  1:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10 - diff_ESS1   -0.120 0.0334 27.4   -0.189  -0.0517  -3.596  0.0013
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
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: FM.z ~ gndr.c + essround.c + gndr.c:essround.c + log_gdp.z.cm:gndr.c +  
##     log_gdp.z.cm:essround.c + log_gdp.z.cm:gndr.c:essround.c +  
##     (gndr.c + essround.c + gndr.c:essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1166535.4 1166731.2 -583249.7 1166499.4    392750 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.8017 -0.6127  0.0043  0.6005  9.9812 
## 
## Random effects:
##  Groups   Name              Variance  Std.Dev. Corr             
##  cntry    (Intercept)       0.1020915 0.31952                   
##           gndr.c            0.0083521 0.09139  -0.33            
##           essround.c        0.0034054 0.05836   0.16 -0.24      
##           gndr.c:essround.c 0.0001219 0.01104  -0.21 -0.01 -0.48
##  Residual                   1.0032771 1.00164                   
## Number of obs: 392768, groups:  cntry, 33
## 
## Fixed effects:
##                                 Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)                     0.040502   0.055810 29.593934   0.726  0.47372    
## gndr.c                          0.423114   0.016444 27.700571  25.730  < 2e-16 ***
## essround.c                      0.001212   0.010245 27.376915   0.118  0.90670    
## gndr.c:essround.c              -0.012540   0.002431 22.751547  -5.159 3.25e-05 ***
## gndr.c:log_gdp.z.cm             0.021487   0.015842 33.054067   1.356  0.18420    
## essround.c:log_gdp.z.cm        -0.034685   0.010315 25.761613  -3.363  0.00242 ** 
## gndr.c:essround.c:log_gdp.z.cm  0.002736   0.002562 22.508586   1.068  0.29694    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c essrn. gnd.:. g.:_.. e.:_..
## gndr.c      -0.321                                   
## essround.c   0.148 -0.230                            
## gndr.c:ssr. -0.162 -0.027 -0.374                     
## gndr.c:l_.. -0.002 -0.007 -0.003  0.005              
## essrnd.:_..  0.011 -0.006  0.000 -0.005 -0.191       
## gndr.:.:_.. -0.003  0.007 -0.004 -0.156 -0.075 -0.339
```

``` r
getFE(mod6_log_GDP,round=2,p.round = 5)
```

```
##                                 Est.   SE    df     t       p    LL    UL
## (Intercept)                     0.04 0.06 29.59  0.73 0.47372 -0.07  0.15
## gndr.c                          0.42 0.02 27.70 25.73 0.00000  0.39  0.46
## essround.c                      0.00 0.01 27.38  0.12 0.90670 -0.02  0.02
## gndr.c:essround.c              -0.01 0.00 22.75 -5.16 0.00003 -0.02 -0.01
## gndr.c:log_gdp.z.cm             0.02 0.02 33.05  1.36 0.18420 -0.01  0.05
## essround.c:log_gdp.z.cm        -0.03 0.01 25.76 -3.36 0.00242 -0.06 -0.01
## gndr.c:essround.c:log_gdp.z.cm  0.00 0.00 22.51  1.07 0.29694  0.00  0.01
```

``` r
getVC(mod6_log_GDP)
```

```
##         grp              var1              var2 sdcor  vcov
## 1     cntry       (Intercept)              <NA>  0.32  0.10
## 2     cntry            gndr.c              <NA>  0.09  0.01
## 3     cntry        essround.c              <NA>  0.06  0.00
## 4     cntry gndr.c:essround.c              <NA>  0.01  0.00
## 5     cntry       (Intercept)            gndr.c -0.33 -0.01
## 6     cntry       (Intercept)        essround.c  0.16  0.00
## 7     cntry       (Intercept) gndr.c:essround.c -0.21  0.00
## 8     cntry            gndr.c        essround.c -0.24  0.00
## 9     cntry            gndr.c gndr.c:essround.c -0.01  0.00
## 10    cntry        essround.c gndr.c:essround.c -0.48  0.00
## 11 Residual              <NA>              <NA>  1.00  1.00
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_log_GDP,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] 23.71188
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_log_GDP,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] 6.649447
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
##         4.5  0.202 0.0900 40.4   0.0202    0.384   2.245  0.0303
##        -4.5 -0.121 0.0818 23.4  -0.2901    0.048  -1.480  0.1523
## 
## log_gdp.z.cm =  0:
##  essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.046 0.0775 30.1  -0.1123    0.204   0.593  0.5575
##        -4.5  0.035 0.0669 16.9  -0.1062    0.176   0.524  0.6072
## 
## log_gdp.z.cm =  1:
##  essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5 -0.110 0.0906 40.1  -0.2933    0.073  -1.215  0.2315
##        -4.5  0.191 0.0811 22.0   0.0230    0.359   2.357  0.0277
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
##  essround.c4.5 - (essround.c-4.5)   0.3231 0.1310 27.8   0.0549   0.5912   2.469  0.0200
## 
## log_gdp.z.cm =  0:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   0.0109 0.0922 27.4  -0.1782   0.2000   0.118  0.9067
## 
## log_gdp.z.cm =  1:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)  -0.3013 0.1310 26.0  -0.5702  -0.0323  -2.303  0.0296
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
##  gndr.c essround.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5        4.5  0.0356 0.0967 39.9  -0.1599   0.2311   0.368  0.7149
##     0.5        4.5  0.3685 0.0850 39.8   0.1966   0.5404   4.333  0.0001
##    -0.5       -4.5 -0.3562 0.0846 21.0  -0.5322  -0.1803  -4.210  0.0004
##     0.5       -4.5  0.1142 0.0816 24.6  -0.0539   0.2822   1.400  0.1741
## 
## log_gdp.z.cm =  0:
##  gndr.c essround.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5        4.5 -0.1374 0.0828 29.7  -0.3065   0.0318  -1.659  0.1076
##     0.5        4.5  0.2293 0.0731 29.7   0.0800   0.3786   3.137  0.0038
##    -0.5       -4.5 -0.2047 0.0692 15.5  -0.3518  -0.0577  -2.959  0.0095
##     0.5       -4.5  0.2748 0.0661 17.7   0.1358   0.4138   4.159  0.0006
## 
## log_gdp.z.cm =  1:
##  gndr.c essround.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5        4.5 -0.3104 0.0974 39.5  -0.5073  -0.1135  -3.187  0.0028
##     0.5        4.5  0.0901 0.0854 39.2  -0.0826   0.2629   1.055  0.2979
##    -0.5       -4.5 -0.0532 0.0837 19.8  -0.2280   0.1215  -0.636  0.5323
##     0.5       -4.5  0.4355 0.0807 23.1   0.2685   0.6024   5.395  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod6_log_GDP,adjust="none",infer=c(T,T))
```

```
## log_gdp.z.cm = -1:
##  contrast                                                 estimate     SE   df lower.CL upper.CL t.ratio
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.3329 0.0279 28.9 -0.38986  -0.2759 -11.951
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   0.3918 0.1370 27.1  0.10996   0.6737   2.852
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.0786 0.1370 29.3 -0.35943   0.2023  -0.572
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.7247 0.1280 25.9  0.46131   0.9881   5.657
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)       0.2543 0.1260 27.3 -0.00468   0.5133   2.014
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.4704 0.0293 29.4 -0.53025  -0.4105 -16.056
##  p.value
##   <.0001
##   0.0082
##   0.5719
##   <.0001
##   0.0540
##   <.0001
## 
## log_gdp.z.cm =  0:
##  contrast                                                 estimate     SE   df lower.CL upper.CL t.ratio
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.3667 0.0195 22.2 -0.40711  -0.3263 -18.800
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   0.0673 0.0968 26.8 -0.13141   0.2661   0.695
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.4122 0.0973 28.8 -0.61131  -0.2131  -4.236
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.4340 0.0899 25.3  0.24906   0.6190   4.830
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      -0.0455 0.0887 26.6 -0.22766   0.1366  -0.513
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.4795 0.0200 24.0 -0.52081  -0.4383 -23.986
##  p.value
##   <.0001
##   0.4928
##   0.0002
##   0.0001
##   0.6120
##   <.0001
## 
## log_gdp.z.cm =  1:
##  contrast                                                 estimate     SE   df lower.CL upper.CL t.ratio
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.4005 0.0264 23.5 -0.45508  -0.3459 -15.157
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  -0.2571 0.1370 25.3 -0.53958   0.0253  -1.874
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.7459 0.1380 27.6 -1.02784  -0.4639  -5.422
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.1433 0.1280 24.1 -0.12048   0.4072   1.121
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      -0.3454 0.1260 25.2 -0.60441  -0.0863  -2.745
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.4887 0.0276 24.7 -0.54567  -0.4318 -17.686
##  p.value
##   <.0001
##   0.0725
##   <.0001
##   0.2733
##   0.0110
##   <.0001
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
##  diff_ESS10    0.333 0.0279 28.9    0.276    0.390  11.951  <.0001
##  diff_ESS1     0.470 0.0293 29.4    0.410    0.530  16.056  <.0001
## 
## log_gdp.z.cm =  0:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10    0.367 0.0195 22.2    0.326    0.407  18.800  <.0001
##  diff_ESS1     0.480 0.0200 24.0    0.438    0.521  23.986  <.0001
## 
## log_gdp.z.cm =  1:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10    0.400 0.0264 23.5    0.346    0.455  15.157  <.0001
##  diff_ESS1     0.489 0.0276 24.7    0.432    0.546  17.686  <.0001
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
##  diff_ESS10 - diff_ESS1  -0.1375 0.0342 28.4   -0.207  -0.0675  -4.023  0.0004
## 
## log_gdp.z.cm =  0:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10 - diff_ESS1  -0.1129 0.0219 22.8   -0.158  -0.0676  -5.159  <.0001
## 
## log_gdp.z.cm =  1:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10 - diff_ESS1  -0.0882 0.0292 17.8   -0.150  -0.0268  -3.022  0.0074
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
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: FM.z ~ gndr.c + gei.z.cm + gndr.c:gei.z.cm + essround.c + gndr.c:essround.c +  
##     (gndr.c + essround.c + gndr.c:essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1136153.3 1136337.9 -568059.7 1136119.3    384361 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.8511 -0.6140  0.0047  0.6024 10.0275 
## 
## Random effects:
##  Groups   Name              Variance Std.Dev. Corr             
##  cntry    (Intercept)       0.083571 0.28909                   
##           gndr.c            0.006934 0.08327  -0.22            
##           essround.c        0.004981 0.07058  -0.07  0.06      
##           gndr.c:essround.c 0.000162 0.01273   0.09 -0.28 -0.61
##  Residual                   0.993717 0.99685                   
## Number of obs: 384378, groups:  cntry, 32
## 
## Fixed effects:
##                     Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)        0.0294720  0.0513430 24.6669144   0.574 0.571151    
## gndr.c             0.4279894  0.0153290 25.5049767  27.920  < 2e-16 ***
## gei.z.cm          -0.1626306  0.0520201 26.7661191  -3.126 0.004232 ** 
## essround.c         0.0005486  0.0125553 24.8397289   0.044 0.965496    
## gndr.c:gei.z.cm    0.0564025  0.0158010 30.7684227   3.570 0.001197 ** 
## gndr.c:essround.c -0.0132893  0.0026718 15.9502289  -4.974 0.000139 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c g.z.cm essrn. gn.:..
## gndr.c      -0.210                            
## gei.z.cm    -0.005  0.002                     
## essround.c  -0.077  0.062 -0.001              
## gndr.c:g.z.  0.001 -0.028 -0.200  0.000       
## gndr.c:ssr.  0.076 -0.251  0.006 -0.517 -0.053
```

``` r
getFE(mod7_GEI)
```

```
##                    Est.   SE    df     t     p    LL    UL
## (Intercept)        0.03 0.05 24.67  0.57 0.571 -0.08  0.14
## gndr.c             0.43 0.02 25.50 27.92 0.000  0.40  0.46
## gei.z.cm          -0.16 0.05 26.77 -3.13 0.004 -0.27 -0.06
## essround.c         0.00 0.01 24.84  0.04 0.965 -0.03  0.03
## gndr.c:gei.z.cm    0.06 0.02 30.77  3.57 0.001  0.02  0.09
## gndr.c:essround.c -0.01 0.00 15.95 -4.97 0.000 -0.02 -0.01
```

``` r
getVC(mod7_GEI)
```

```
##         grp              var1              var2 sdcor  vcov
## 1     cntry       (Intercept)              <NA>  0.29  0.08
## 2     cntry            gndr.c              <NA>  0.08  0.01
## 3     cntry        essround.c              <NA>  0.07  0.00
## 4     cntry gndr.c:essround.c              <NA>  0.01  0.00
## 5     cntry       (Intercept)            gndr.c -0.22 -0.01
## 6     cntry       (Intercept)        essround.c -0.07  0.00
## 7     cntry       (Intercept) gndr.c:essround.c  0.09  0.00
## 8     cntry            gndr.c        essround.c  0.06  0.00
## 9     cntry            gndr.c gndr.c:essround.c -0.28  0.00
## 10    cntry        essround.c gndr.c:essround.c -0.61  0.00
## 11 Residual              <NA>              <NA>  1.00  0.99
```

``` r
anova(mod2_GEI,mod7_GEI)
```

```
## Data: diff_dat
## Models:
## mod2_GEI: FM.z ~ gndr.c + gei.z.cm + gndr.c:gei.z.cm + (gndr.c | cntry)
## mod7_GEI: FM.z ~ gndr.c + gei.z.cm + gndr.c:gei.z.cm + essround.c + gndr.c:essround.c + (gndr.c + essround.c + gndr.c:essround.c | cntry)
##          npar     AIC     BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod2_GEI    8 1140107 1140194 -570046   1140091                         
## mod7_GEI   17 1136153 1136338 -568060   1136119 3971.7  9  < 2.2e-16 ***
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
##         4.5 -0.152 0.0766 24.1  -0.3103  0.00595  -1.986  0.0585
##        -4.5 -0.217 0.0837 11.8  -0.3996 -0.03419  -2.590  0.0239
## 
## gndr.c =  0.5:
##  essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.216 0.0710 24.7   0.0698  0.36230   3.044  0.0055
##        -4.5  0.271 0.0760 13.8   0.1077  0.43409   3.565  0.0032
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_mod7_GEI,adjust="none",infer=c(T,T))
```

```
## gndr.c = -0.5:
##  contrast                         estimate    SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   0.0647 0.120 24.1   -0.182    0.312   0.541  0.5935
## 
## gndr.c =  0.5:
##  contrast                         estimate    SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)  -0.0549 0.107 24.7   -0.276    0.166  -0.511  0.6136
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
##    -0.5        4.5 -0.152 0.0766 24.1  -0.3103  0.00595  -1.986  0.0585
##     0.5        4.5  0.216 0.0710 24.7   0.0698  0.36230   3.044  0.0055
##    -0.5       -4.5 -0.217 0.0837 11.8  -0.3996 -0.03419  -2.590  0.0239
##     0.5       -4.5  0.271 0.0760 13.8   0.1077  0.43409   3.565  0.0032
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod7_GEI,adjust="none",infer=c(T,T))
```

```
##  contrast                                                 estimate     SE   df lower.CL upper.CL t.ratio
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.3682 0.0169 23.4   -0.403   -0.333 -21.727
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   0.0647 0.1200 24.1   -0.182    0.312   0.541
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.4231 0.1130 26.1   -0.655   -0.191  -3.741
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.4329 0.1150 22.6    0.195    0.671   3.765
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      -0.0549 0.1070 24.7   -0.276    0.166  -0.511
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.4878 0.0217 15.9   -0.534   -0.442 -22.455
##  p.value
##   <.0001
##   0.5935
##   0.0009
##   0.0010
##   0.6136
##   <.0001
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
##  diff_ESS10    0.368 0.0169 23.4    0.333    0.403  21.727  <.0001
##  diff_ESS1     0.488 0.0217 15.9    0.442    0.534  22.455  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# Test for whether and how much the differences have changed
pairs(diff_mod7_GEI,infer=c(T,T))
```

```
##  contrast               estimate    SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10 - diff_ESS1    -0.12 0.024 15.9   -0.171  -0.0686  -4.974  0.0001
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
## (Intercept)        0.05 0.04 31.25  1.08 0.290 -0.04  0.14
## gndr.c             0.43 0.02 29.11 26.69 0.000  0.39  0.46
## gei.z.cm          -0.13 0.04 32.26 -3.06 0.004 -0.21 -0.04
## essround.c         0.00 0.01 26.31  0.28 0.785 -0.02  0.03
## gei.z.cmc         -0.14 0.08 22.30 -1.66 0.112 -0.31  0.03
## gndr.c:gei.z.cm    0.05 0.01 32.00  3.50 0.001  0.02  0.08
## gndr.c:essround.c -0.01 0.01 10.22 -2.54 0.029 -0.03  0.00
## gndr.c:gei.z.cmc   0.05 0.03 13.71  1.39 0.185 -0.03  0.12
```

``` r
getVC(mod8_GEI)
```

```
##         grp              var1              var2 sdcor  vcov
## 1     cntry       (Intercept)              <NA>  0.25  0.06
## 2     cntry            gndr.c              <NA>  0.09  0.01
## 3     cntry        essround.c              <NA>  0.07  0.00
## 4     cntry         gei.z.cmc              <NA>  0.45  0.20
## 5     cntry gndr.c:essround.c              <NA>  0.03  0.00
## 6     cntry  gndr.c:gei.z.cmc              <NA>  0.16  0.02
## 7     cntry       (Intercept)            gndr.c -0.04  0.00
## 8     cntry       (Intercept)        essround.c  0.01  0.00
## 9     cntry       (Intercept)         gei.z.cmc  0.18  0.02
## 10    cntry       (Intercept) gndr.c:essround.c  0.00  0.00
## 11    cntry       (Intercept)  gndr.c:gei.z.cmc -0.05  0.00
## 12    cntry            gndr.c        essround.c  0.21  0.00
## 13    cntry            gndr.c         gei.z.cmc -0.34 -0.01
## 14    cntry            gndr.c gndr.c:essround.c -0.60  0.00
## 15    cntry            gndr.c  gndr.c:gei.z.cmc  0.61  0.01
## 16    cntry        essround.c         gei.z.cmc -0.86 -0.03
## 17    cntry        essround.c gndr.c:essround.c -0.27  0.00
## 18    cntry        essround.c  gndr.c:gei.z.cmc  0.29  0.00
## 19    cntry         gei.z.cmc gndr.c:essround.c  0.44  0.01
## 20    cntry         gei.z.cmc  gndr.c:gei.z.cmc -0.50 -0.04
## 21    cntry gndr.c:essround.c  gndr.c:gei.z.cmc -0.99  0.00
## 22 Residual              <NA>              <NA>  1.00  0.99
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
##          npar     AIC     BIC  logLik -2*log(L)   Chisq Df Pr(>Chisq)    
## mod2_GEI    8 1140107 1140194 -570046   1140091                          
## mod7_GEI   17 1136153 1136338 -568060   1136119 3971.72  9  < 2.2e-16 ***
## mod8_GEI   30 1135876 1136202 -567908   1135816  303.37 13  < 2.2e-16 ***
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
##         4.5 -0.120 0.0726 29.3   -0.268   0.0285  -1.651  0.1094
##        -4.5 -0.213 0.0766 25.2   -0.371  -0.0554  -2.782  0.0101
## 
## gndr.c =  0.5:
##  essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.244 0.0703 29.9    0.101   0.3881   3.475  0.0016
##        -4.5  0.276 0.0699 24.6    0.132   0.4198   3.949  0.0006
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
##  essround.c4.5 - (essround.c-4.5)   0.0932 0.120 24.7   -0.154    0.340   0.778  0.4439
## 
## gndr.c =  0.5:
##  contrast                         estimate    SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)  -0.0314 0.109 25.7   -0.256    0.193  -0.288  0.7759
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
##    -0.5        4.5 -0.120 0.0726 29.3   -0.268   0.0285  -1.651  0.1094
##     0.5        4.5  0.244 0.0703 29.9    0.101   0.3881   3.475  0.0016
##    -0.5       -4.5 -0.213 0.0766 25.2   -0.371  -0.0554  -2.782  0.0101
##     0.5       -4.5  0.276 0.0699 24.6    0.132   0.4198   3.949  0.0006
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod8_GEI,adjust="none",infer=c(T,T))
```

```
##  contrast                                                 estimate     SE   df lower.CL upper.CL t.ratio
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.3642 0.0209 11.9   -0.410   -0.319 -17.405
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   0.0932 0.1200 24.7   -0.154    0.340   0.778
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.3956 0.1100 26.6   -0.621   -0.170  -3.598
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.4574 0.1160 25.8    0.219    0.696   3.942
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      -0.0314 0.1090 25.7   -0.256    0.193  -0.288
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.4888 0.0358 15.1   -0.565   -0.413 -13.669
##  p.value
##   <.0001
##   0.4439
##   0.0013
##   0.0006
##   0.7759
##   <.0001
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
##  diff_ESS10    0.364 0.0209 11.9    0.319    0.410  17.405  <.0001
##  diff_ESS1     0.489 0.0358 15.1    0.413    0.565  13.669  <.0001
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
##  diff_ESS10 - diff_ESS1   -0.125 0.0491 10.2   -0.234  -0.0154  -2.536  0.0291
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
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: FM.z ~ gndr.c + gggi.z.cm + gndr.c:gggi.z.cm + essround.c + gndr.c:essround.c +  
##     (gndr.c + essround.c + gndr.c:essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##  833015.9  833195.2 -416491.0  832981.9    280629 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.8273 -0.6124  0.0064  0.6025  8.1716 
## 
## Random effects:
##  Groups   Name              Variance  Std.Dev. Corr             
##  cntry    (Intercept)       0.0620275 0.24905                   
##           gndr.c            0.0058741 0.07664  -0.02            
##           essround.c        0.0021709 0.04659   0.24  0.06      
##           gndr.c:essround.c 0.0001797 0.01341  -0.14 -0.07 -0.52
##  Residual                   0.9997792 0.99989                   
## Number of obs: 280646, groups:  cntry, 33
## 
## Fixed effects:
##                    Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)        0.065206   0.043750 32.063980   1.490  0.14589    
## gndr.c             0.416798   0.014219 32.867961  29.312  < 2e-16 ***
## gggi.z.cm         -0.096960   0.043127 33.158801  -2.248  0.03133 *  
## essround.c        -0.010589   0.008507 30.422569  -1.245  0.22274    
## gndr.c:gggi.z.cm   0.052585   0.014821 36.633605   3.548  0.00108 ** 
## gndr.c:essround.c -0.007309   0.003143 25.297919  -2.326  0.02833 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c ggg.z. essrn. gn.:..
## gndr.c      -0.024                            
## gggi.z.cm   -0.006  0.001                     
## essround.c   0.233  0.055  0.003              
## gndr.c:gg..  0.004 -0.015 -0.039  0.009       
## gndr.c:ssr. -0.105 -0.093 -0.001 -0.389 -0.032
```

``` r
getFE(mod7_GGGI)
```

```
##                    Est.   SE    df     t     p    LL    UL
## (Intercept)        0.07 0.04 32.06  1.49 0.146 -0.02  0.15
## gndr.c             0.42 0.01 32.87 29.31 0.000  0.39  0.45
## gggi.z.cm         -0.10 0.04 33.16 -2.25 0.031 -0.18 -0.01
## essround.c        -0.01 0.01 30.42 -1.24 0.223 -0.03  0.01
## gndr.c:gggi.z.cm   0.05 0.01 36.63  3.55 0.001  0.02  0.08
## gndr.c:essround.c -0.01 0.00 25.30 -2.33 0.028 -0.01  0.00
```

``` r
getVC(mod7_GGGI)
```

```
##         grp              var1              var2 sdcor vcov
## 1     cntry       (Intercept)              <NA>  0.25 0.06
## 2     cntry            gndr.c              <NA>  0.08 0.01
## 3     cntry        essround.c              <NA>  0.05 0.00
## 4     cntry gndr.c:essround.c              <NA>  0.01 0.00
## 5     cntry       (Intercept)            gndr.c -0.02 0.00
## 6     cntry       (Intercept)        essround.c  0.24 0.00
## 7     cntry       (Intercept) gndr.c:essround.c -0.14 0.00
## 8     cntry            gndr.c        essround.c  0.06 0.00
## 9     cntry            gndr.c gndr.c:essround.c -0.07 0.00
## 10    cntry        essround.c gndr.c:essround.c -0.52 0.00
## 11 Residual              <NA>              <NA>  1.00 1.00
```

``` r
anova(mod2_GGGI,mod7_GGGI)
```

```
## Data: diff_dat
## Models:
## mod2_GGGI: FM.z ~ gndr.c + gggi.z.cm + gndr.c:gggi.z.cm + (gndr.c | cntry)
## mod7_GGGI: FM.z ~ gndr.c + gggi.z.cm + gndr.c:gggi.z.cm + essround.c + gndr.c:essround.c + (gndr.c + essround.c + gndr.c:essround.c | cntry)
##           npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod2_GGGI    8 835098 835183 -417541    835082                         
## mod7_GGGI   17 833016 833195 -416491    832982 2100.4  9  < 2.2e-16 ***
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
##         4.5 -0.174 0.0672 27.7  -0.3121 -0.03672  -2.596  0.0149
##        -4.5 -0.112 0.0538 27.8  -0.2223 -0.00165  -2.080  0.0469
## 
## gndr.c =  0.5:
##  essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.210 0.0632 27.0   0.0799  0.33912   3.317  0.0026
##        -4.5  0.338 0.0502 27.1   0.2348  0.44064   6.730  <.0001
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
##  essround.c4.5 - (essround.c-4.5)  -0.0624 0.0831 30.2   -0.232   0.1072  -0.751  0.4584
## 
## gndr.c =  0.5:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)  -0.1282 0.0722 30.1   -0.276   0.0193  -1.774  0.0861
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
##    -0.5        4.5 -0.174 0.0672 27.7  -0.3121 -0.03672  -2.596  0.0149
##     0.5        4.5  0.210 0.0632 27.0   0.0799  0.33912   3.317  0.0026
##    -0.5       -4.5 -0.112 0.0538 27.8  -0.2223 -0.00165  -2.080  0.0469
##     0.5       -4.5  0.338 0.0502 27.1   0.2348  0.44064   6.730  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod7_GGGI,adjust="none",infer=c(T,T))
```

```
##  contrast                                                 estimate     SE   df lower.CL upper.CL t.ratio
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.3839 0.0191 29.0   -0.423  -0.3448 -20.101
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  -0.0624 0.0831 30.2   -0.232   0.1072  -0.751
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.5121 0.0771 29.9   -0.670  -0.3546  -6.642
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.3215 0.0786 30.4    0.161   0.4820   4.088
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      -0.1282 0.0722 30.1   -0.276   0.0193  -1.774
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.4497 0.0210 30.5   -0.492  -0.4069 -21.448
##  p.value
##   <.0001
##   0.4584
##   <.0001
##   0.0003
##   0.0861
##   <.0001
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
##  diff_ESS10    0.384 0.0191 29.0    0.345    0.423  20.101  <.0001
##  diff_ESS1     0.450 0.0210 30.5    0.407    0.492  21.448  <.0001
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
##  diff_ESS10 - diff_ESS1  -0.0658 0.0283 25.3   -0.124 -0.00757  -2.326  0.0283
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
## (Intercept)        0.06 0.04 32.32  1.28 0.211 -0.03  0.15
## gndr.c             0.42 0.02 31.41 26.80 0.000  0.39  0.45
## gggi.z.cm         -0.11 0.04 33.91 -2.52 0.017 -0.20 -0.02
## essround.c        -0.01 0.01 29.32 -0.55 0.588 -0.03  0.02
## gggi.z.cmc        -0.05 0.04 28.23 -1.15 0.260 -0.13  0.04
## gndr.c:gggi.z.cm   0.05 0.01 38.55  3.59 0.001  0.02  0.08
## gndr.c:essround.c -0.01 0.00 22.07 -2.73 0.012 -0.02  0.00
## gndr.c:gggi.z.cmc  0.04 0.02 25.07  1.93 0.065  0.00  0.08
```

``` r
getVC(mod8_GGGI)
```

```
##         grp              var1              var2 sdcor  vcov
## 1     cntry       (Intercept)              <NA>  0.25  0.06
## 2     cntry            gndr.c              <NA>  0.09  0.01
## 3     cntry        essround.c              <NA>  0.06  0.00
## 4     cntry        gggi.z.cmc              <NA>  0.21  0.04
## 5     cntry gndr.c:essround.c              <NA>  0.02  0.00
## 6     cntry gndr.c:gggi.z.cmc              <NA>  0.08  0.01
## 7     cntry       (Intercept)            gndr.c -0.02  0.00
## 8     cntry       (Intercept)        essround.c  0.10  0.00
## 9     cntry       (Intercept)        gggi.z.cmc  0.17  0.01
## 10    cntry       (Intercept) gndr.c:essround.c -0.13  0.00
## 11    cntry       (Intercept) gndr.c:gggi.z.cmc  0.15  0.00
## 12    cntry            gndr.c        essround.c  0.12  0.00
## 13    cntry            gndr.c        gggi.z.cmc -0.01  0.00
## 14    cntry            gndr.c gndr.c:essround.c -0.46  0.00
## 15    cntry            gndr.c gndr.c:gggi.z.cmc  0.73  0.00
## 16    cntry        essround.c        gggi.z.cmc -0.60 -0.01
## 17    cntry        essround.c gndr.c:essround.c -0.58  0.00
## 18    cntry        essround.c gndr.c:gggi.z.cmc  0.66  0.00
## 19    cntry        gggi.z.cmc gndr.c:essround.c  0.04  0.00
## 20    cntry        gggi.z.cmc gndr.c:gggi.z.cmc -0.21  0.00
## 21    cntry gndr.c:essround.c gndr.c:gggi.z.cmc -0.90  0.00
## 22 Residual              <NA>              <NA>  1.00  1.00
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
##           npar    AIC    BIC  logLik -2*log(L)   Chisq Df Pr(>Chisq)    
## mod2_GGGI    8 835098 835183 -417541    835082                          
## mod7_GGGI   17 833016 833195 -416491    832982 2100.37  9  < 2.2e-16 ***
## mod8_GGGI   30 832789 833105 -416364    832729  252.94 13  < 2.2e-16 ***
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
##         4.5 -0.154 0.0729 27.4  -0.3031 -0.00418  -2.108  0.0443
##        -4.5 -0.156 0.0687 23.7  -0.2977 -0.01404  -2.270  0.0326
## 
## gndr.c =  0.5:
##  essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.212 0.0661 26.7   0.0764  0.34781   3.210  0.0034
##        -4.5  0.322 0.0615 23.2   0.1952  0.44969   5.240  <.0001
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
##  essround.c4.5 - (essround.c-4.5)  0.00225 0.1090 29.0   -0.221   0.2259   0.021  0.9838
## 
## gndr.c =  0.5:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5) -0.11033 0.0913 28.5   -0.297   0.0766  -1.208  0.2369
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
##    -0.5        4.5 -0.154 0.0729 27.4  -0.3031 -0.00418  -2.108  0.0443
##     0.5        4.5  0.212 0.0661 26.7   0.0764  0.34781   3.210  0.0034
##    -0.5       -4.5 -0.156 0.0687 23.7  -0.2977 -0.01404  -2.270  0.0326
##     0.5       -4.5  0.322 0.0615 23.2   0.1952  0.44969   5.240  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod8_GGGI,adjust="none",infer=c(T,T))
```

```
##  contrast                                                 estimate     SE   df lower.CL upper.CL t.ratio
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5     -0.36576 0.0200 27.6   -0.407  -0.3247 -18.266
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  0.00225 0.1090 29.0   -0.221   0.2259   0.021
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)  -0.47609 0.0981 28.8   -0.677  -0.2753  -4.852
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)     0.36800 0.1020 29.5    0.160   0.5756   3.623
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)     -0.11033 0.0913 28.5   -0.297   0.0766  -1.208
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5) -0.47834 0.0308 23.6   -0.542  -0.4147 -15.539
##  p.value
##   <.0001
##   0.9838
##   <.0001
##   0.0011
##   0.2369
##   <.0001
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
##  diff_ESS10    0.366 0.0200 27.6    0.325    0.407  18.266  <.0001
##  diff_ESS1     0.478 0.0308 23.6    0.415    0.542  15.539  <.0001
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
##  diff_ESS10 - diff_ESS1   -0.113 0.0413 22.1   -0.198   -0.027  -2.726  0.0123
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
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: FM.z ~ gndr.c + gdi.z.cm + gndr.c:gdi.z.cm + essround.c + gndr.c:essround.c +  
##     (gndr.c + essround.c + gndr.c:essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1166540.7 1166725.6 -583253.3 1166506.7    392751 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.8007 -0.6128  0.0043  0.6006  9.9773 
## 
## Random effects:
##  Groups   Name              Variance  Std.Dev. Corr             
##  cntry    (Intercept)       0.0987145 0.31419                   
##           gndr.c            0.0088768 0.09422  -0.42            
##           essround.c        0.0043750 0.06614   0.44 -0.35      
##           gndr.c:essround.c 0.0001403 0.01184  -0.30  0.10 -0.58
##  Residual                   1.0032936 1.00165                   
## Number of obs: 392768, groups:  cntry, 33
## 
## Fixed effects:
##                     Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)        0.0455594  0.0548657 29.5807169   0.830    0.413    
## gndr.c             0.4228211  0.0169207 29.0531993  24.988  < 2e-16 ***
## gdi.z.cm           0.0597107  0.0501229 21.1354502   1.191    0.247    
## essround.c         0.0007826  0.0115858 27.6530653   0.068    0.947    
## gndr.c:gdi.z.cm    0.0120219  0.0165575 29.4759242   0.726    0.474    
## gndr.c:essround.c -0.0127472  0.0025097 19.6774982  -5.079 6.01e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c gd.z.c essrn. gn.:..
## gndr.c      -0.409                            
## gdi.z.cm    -0.006  0.001                     
## essround.c   0.432 -0.337  0.006              
## gndr.c:gd..  0.002 -0.013 -0.302 -0.003       
## gndr.c:ssr. -0.247  0.062 -0.003 -0.475 -0.006
```

``` r
getFE(mod7_GDI)
```

```
##                    Est.   SE    df     t     p    LL    UL
## (Intercept)        0.05 0.05 29.58  0.83 0.413 -0.07  0.16
## gndr.c             0.42 0.02 29.05 24.99 0.000  0.39  0.46
## gdi.z.cm           0.06 0.05 21.14  1.19 0.247 -0.04  0.16
## essround.c         0.00 0.01 27.65  0.07 0.947 -0.02  0.02
## gndr.c:gdi.z.cm    0.01 0.02 29.48  0.73 0.474 -0.02  0.05
## gndr.c:essround.c -0.01 0.00 19.68 -5.08 0.000 -0.02 -0.01
```

``` r
getVC(mod7_GDI)
```

```
##         grp              var1              var2 sdcor  vcov
## 1     cntry       (Intercept)              <NA>  0.31  0.10
## 2     cntry            gndr.c              <NA>  0.09  0.01
## 3     cntry        essround.c              <NA>  0.07  0.00
## 4     cntry gndr.c:essround.c              <NA>  0.01  0.00
## 5     cntry       (Intercept)            gndr.c -0.42 -0.01
## 6     cntry       (Intercept)        essround.c  0.44  0.01
## 7     cntry       (Intercept) gndr.c:essround.c -0.30  0.00
## 8     cntry            gndr.c        essround.c -0.35  0.00
## 9     cntry            gndr.c gndr.c:essround.c  0.10  0.00
## 10    cntry        essround.c gndr.c:essround.c -0.58  0.00
## 11 Residual              <NA>              <NA>  1.00  1.00
```

``` r
anova(mod2_GDI,mod7_GDI)
```

```
## Data: diff_dat
## Models:
## mod2_GDI: FM.z ~ gndr.c + gdi.z.cm + gndr.c:gdi.z.cm + (gndr.c | cntry)
## mod7_GDI: FM.z ~ gndr.c + gdi.z.cm + gndr.c:gdi.z.cm + essround.c + gndr.c:essround.c + (gndr.c + essround.c + gndr.c:essround.c | cntry)
##          npar     AIC     BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod2_GDI    8 1170512 1170599 -585248   1170496                         
## mod7_GDI   17 1166541 1166726 -583253   1166507 3989.3  9  < 2.2e-16 ***
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
##         4.5 -0.134 0.0970 28.7   -0.332   0.0649  -1.377  0.1791
##        -4.5 -0.198 0.0597 12.7   -0.327  -0.0688  -3.318  0.0057
## 
## gndr.c =  0.5:
##  essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.232 0.0848 28.0    0.058   0.4056   2.732  0.0108
##        -4.5  0.282 0.0561 13.8    0.162   0.4026   5.030  0.0002
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
##  essround.c4.5 - (essround.c-4.5)   0.0644 0.1100 27.2   -0.161    0.290   0.585  0.5634
## 
## gndr.c =  0.5:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)  -0.0503 0.0994 27.1   -0.254    0.154  -0.506  0.6168
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
##    -0.5        4.5 -0.134 0.0970 28.7   -0.332   0.0649  -1.377  0.1791
##     0.5        4.5  0.232 0.0848 28.0    0.058   0.4056   2.732  0.0108
##    -0.5       -4.5 -0.198 0.0597 12.7   -0.327  -0.0688  -3.318  0.0057
##     0.5       -4.5  0.282 0.0561 13.8    0.162   0.4026   5.030  0.0002
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod7_GDI,adjust="none",infer=c(T,T))
```

```
##  contrast                                                 estimate     SE   df lower.CL upper.CL t.ratio
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.3655 0.0209 21.5   -0.409   -0.322 -17.473
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   0.0644 0.1100 27.2   -0.161    0.290   0.585
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.4158 0.1110 28.9   -0.643   -0.188  -3.741
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.4299 0.0998 25.5    0.224    0.635   4.305
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      -0.0503 0.0994 27.1   -0.254    0.154  -0.506
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.4802 0.0198 24.2   -0.521   -0.439 -24.308
##  p.value
##   <.0001
##   0.5634
##   0.0008
##   0.0002
##   0.6168
##   <.0001
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
##  diff_ESS10    0.365 0.0209 21.5    0.322    0.409  17.473  <.0001
##  diff_ESS1     0.480 0.0198 24.2    0.439    0.521  24.308  <.0001
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
##  diff_ESS10 - diff_ESS1   -0.115 0.0226 19.7   -0.162  -0.0676  -5.079  0.0001
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
```

```
## Warning in checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv, : Model failed to converge
## with max|grad| = 0.00529081 (tol = 0.002, component 1)
```

``` r
getFE(mod8_GDI)
```

```
##                    Est.   SE    df     t     p    LL    UL
## (Intercept)        0.06 0.05 32.41  1.16 0.254 -0.05  0.17
## gndr.c             0.42 0.02 29.26 25.34 0.000  0.39  0.46
## gdi.z.cm          -0.04 0.05 34.30 -0.86 0.393 -0.14  0.05
## essround.c        -0.01 0.01 25.36 -0.72 0.478 -0.02  0.01
## gdi.z.cmc          0.00 0.08 25.87  0.05 0.961 -0.16  0.17
## gndr.c:gdi.z.cm    0.02 0.02 34.19  1.16 0.255 -0.01  0.05
## gndr.c:essround.c -0.01 0.00 20.54 -6.14 0.000 -0.02 -0.01
## gndr.c:gdi.z.cmc   0.00 0.03 15.84  0.04 0.967 -0.07  0.07
```

``` r
getVC(mod8_GDI)
```

```
##         grp              var1              var2 sdcor  vcov
## 1     cntry       (Intercept)              <NA>  0.30  0.09
## 2     cntry            gndr.c              <NA>  0.09  0.01
## 3     cntry        essround.c              <NA>  0.04  0.00
## 4     cntry         gdi.z.cmc              <NA>  0.45  0.21
## 5     cntry gndr.c:essround.c              <NA>  0.01  0.00
## 6     cntry  gndr.c:gdi.z.cmc              <NA>  0.15  0.02
## 7     cntry       (Intercept)            gndr.c -0.32 -0.01
## 8     cntry       (Intercept)        essround.c  0.36  0.00
## 9     cntry       (Intercept)         gdi.z.cmc  0.29  0.04
## 10    cntry       (Intercept) gndr.c:essround.c -0.27  0.00
## 11    cntry       (Intercept)  gndr.c:gdi.z.cmc -0.14 -0.01
## 12    cntry            gndr.c        essround.c -0.47  0.00
## 13    cntry            gndr.c         gdi.z.cmc  0.29  0.01
## 14    cntry            gndr.c gndr.c:essround.c  0.09  0.00
## 15    cntry            gndr.c  gndr.c:gdi.z.cmc -0.26  0.00
## 16    cntry        essround.c         gdi.z.cmc -0.38 -0.01
## 17    cntry        essround.c gndr.c:essround.c -0.44  0.00
## 18    cntry        essround.c  gndr.c:gdi.z.cmc  0.32  0.00
## 19    cntry         gdi.z.cmc gndr.c:essround.c  0.14  0.00
## 20    cntry         gdi.z.cmc  gndr.c:gdi.z.cmc -0.83 -0.06
## 21    cntry gndr.c:essround.c  gndr.c:gdi.z.cmc -0.06  0.00
## 22 Residual              <NA>              <NA>  1.00  1.00
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
##          npar     AIC     BIC  logLik -2*log(L)   Chisq Df Pr(>Chisq)    
## mod2_GDI    8 1170512 1170599 -585248   1170496                          
## mod7_GDI   17 1166541 1166726 -583253   1166507 3989.34  9  < 2.2e-16 ***
## mod8_GDI   30 1166222 1166549 -583081   1166162  344.18 13  < 2.2e-16 ***
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
##         4.5 -0.150 0.0771 27.4  -0.3085  0.00782  -1.949  0.0616
##        -4.5 -0.153 0.0536 22.1  -0.2644 -0.04209  -2.858  0.0091
## 
## gndr.c =  0.5:
##  essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.221 0.0677 27.8   0.0819  0.35935   3.258  0.0030
##        -4.5  0.325 0.0530 23.9   0.2158  0.43453   6.138  <.0001
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
##  essround.c4.5 - (essround.c-4.5)  0.00291 0.0735 24.1   -0.149   0.1546   0.040  0.9688
## 
## gndr.c =  0.5:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5) -0.10457 0.0687 25.9   -0.246   0.0367  -1.521  0.1403
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
##    -0.5        4.5 -0.150 0.0771 27.4  -0.3085  0.00782  -1.949  0.0616
##     0.5        4.5  0.221 0.0677 27.8   0.0819  0.35935   3.258  0.0030
##    -0.5       -4.5 -0.153 0.0536 22.1  -0.2644 -0.04209  -2.858  0.0091
##     0.5       -4.5  0.325 0.0530 23.9   0.2158  0.43453   6.138  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod8_GDI,adjust="none",infer=c(T,T))
```

```
##  contrast                                                 estimate     SE   df lower.CL upper.CL t.ratio
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5     -0.37092 0.0190 23.9   -0.410  -0.3317 -19.518
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  0.00291 0.0735 24.1   -0.149   0.1546   0.040
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)  -0.47549 0.0793 24.6   -0.639  -0.3121  -5.999
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)     0.37383 0.0652 27.2    0.240   0.5076   5.734
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)     -0.10457 0.0687 25.9   -0.246   0.0367  -1.521
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5) -0.47840 0.0188 24.9   -0.517  -0.4397 -25.439
##  p.value
##   <.0001
##   0.9688
##   <.0001
##   <.0001
##   0.1403
##   <.0001
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
##  diff_ESS10    0.371 0.0190 23.9    0.332    0.410  19.518  <.0001
##  diff_ESS1     0.478 0.0188 24.9    0.440    0.517  25.439  <.0001
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
##  diff_ESS10 - diff_ESS1   -0.107 0.0175 20.5   -0.144   -0.071  -6.143  <.0001
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
```

```
## Warning in checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv, : Model failed to converge
## with max|grad| = 0.00248095 (tol = 0.002, component 1)
```

``` r
summary(mod7_log_GDP)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: FM.z ~ gndr.c + log_gdp.z.cm + gndr.c:log_gdp.z.cm + essround.c +  
##     gndr.c:essround.c + (gndr.c + essround.c + gndr.c:essround.c |      cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1166541.8 1166726.7 -583253.9 1166507.8    392751 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.8011 -0.6128  0.0043  0.6005  9.9819 
## 
## Random effects:
##  Groups   Name              Variance  Std.Dev. Corr             
##  cntry    (Intercept)       0.0913466 0.30224                   
##           gndr.c            0.0084489 0.09192  -0.32            
##           essround.c        0.0044919 0.06702   0.26 -0.32      
##           gndr.c:essround.c 0.0001297 0.01139  -0.27  0.05 -0.52
##  Residual                   1.0032896 1.00164                   
## Number of obs: 392768, groups:  cntry, 33
## 
## Fixed effects:
##                      Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)          0.042882   0.052810 28.165080   0.812    0.424    
## gndr.c               0.422529   0.016532 27.688784  25.559  < 2e-16 ***
## log_gdp.z.cm        -0.047245   0.051435 23.498585  -0.919    0.368    
## essround.c           0.001116   0.011743 27.090692   0.095    0.925    
## gndr.c:log_gdp.z.cm  0.014333   0.015897 28.602875   0.902    0.375    
## gndr.c:essround.c   -0.012504   0.002445 21.148615  -5.115 4.48e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c lg_g.. essrn. g.:_..
## gndr.c      -0.310                            
## lg_gdp.z.cm  0.009 -0.003                     
## essround.c   0.256 -0.311  0.013              
## gndr.c:l_.. -0.001 -0.007 -0.269 -0.007       
## gndr.c:ssr. -0.217  0.019 -0.024 -0.423 -0.014
## optimizer (bobyqa) convergence code: 0 (OK)
## Model failed to converge with max|grad| = 0.00248095 (tol = 0.002, component 1)
```

``` r
getFE(mod7_log_GDP)
```

```
##                      Est.   SE    df     t     p    LL    UL
## (Intercept)          0.04 0.05 28.17  0.81 0.424 -0.07  0.15
## gndr.c               0.42 0.02 27.69 25.56 0.000  0.39  0.46
## log_gdp.z.cm        -0.05 0.05 23.50 -0.92 0.368 -0.15  0.06
## essround.c           0.00 0.01 27.09  0.10 0.925 -0.02  0.03
## gndr.c:log_gdp.z.cm  0.01 0.02 28.60  0.90 0.375 -0.02  0.05
## gndr.c:essround.c   -0.01 0.00 21.15 -5.11 0.000 -0.02 -0.01
```

``` r
getVC(mod7_log_GDP)
```

```
##         grp              var1              var2 sdcor  vcov
## 1     cntry       (Intercept)              <NA>  0.30  0.09
## 2     cntry            gndr.c              <NA>  0.09  0.01
## 3     cntry        essround.c              <NA>  0.07  0.00
## 4     cntry gndr.c:essround.c              <NA>  0.01  0.00
## 5     cntry       (Intercept)            gndr.c -0.32 -0.01
## 6     cntry       (Intercept)        essround.c  0.26  0.01
## 7     cntry       (Intercept) gndr.c:essround.c -0.27  0.00
## 8     cntry            gndr.c        essround.c -0.32  0.00
## 9     cntry            gndr.c gndr.c:essround.c  0.05  0.00
## 10    cntry        essround.c gndr.c:essround.c -0.52  0.00
## 11 Residual              <NA>              <NA>  1.00  1.00
```

``` r
anova(mod2_log_GDP,mod7_log_GDP)
```

```
## Data: diff_dat
## Models:
## mod2_log_GDP: FM.z ~ gndr.c + log_gdp.z.cm + gndr.c:log_gdp.z.cm + (gndr.c | cntry)
## mod7_log_GDP: FM.z ~ gndr.c + log_gdp.z.cm + gndr.c:log_gdp.z.cm + essround.c + gndr.c:essround.c + (gndr.c + essround.c + gndr.c:essround.c | cntry)
##              npar     AIC     BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod2_log_GDP    8 1170508 1170595 -585246   1170492                         
## mod7_log_GDP   17 1166542 1166727 -583254   1166508 3984.5  9  < 2.2e-16 ***
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
##         4.5 -0.135 0.0896 26.3  -0.3193   0.0488  -1.510  0.1430
##        -4.5 -0.202 0.0661 12.2  -0.3454  -0.0577  -3.049  0.0100
## 
## gndr.c =  0.5:
##  essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.231 0.0787 26.0   0.0693   0.3928   2.935  0.0069
##        -4.5  0.277 0.0643 14.2   0.1396   0.4149   4.313  0.0007
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_mod7_log_GDP,adjust="none",infer=c(T,T))
```

```
## gndr.c = -0.5:
##  contrast                         estimate    SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   0.0663 0.111 26.7   -0.161    0.294   0.598  0.5546
## 
## gndr.c =  0.5:
##  contrast                         estimate    SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)  -0.0462 0.102 26.5   -0.255    0.162  -0.455  0.6526
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
##    -0.5        4.5 -0.135 0.0896 26.3  -0.3193   0.0488  -1.510  0.1430
##     0.5        4.5  0.231 0.0787 26.0   0.0693   0.3928   2.935  0.0069
##    -0.5       -4.5 -0.202 0.0661 12.2  -0.3454  -0.0577  -3.049  0.0100
##     0.5       -4.5  0.277 0.0643 14.2   0.1396   0.4149   4.313  0.0007
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod7_log_GDP,adjust="none",infer=c(T,T))
```

```
##  contrast                                                 estimate     SE   df lower.CL upper.CL t.ratio
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.3663 0.0200 19.7   -0.408   -0.324 -18.286
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   0.0663 0.1110 26.7   -0.161    0.294   0.598
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.4125 0.1120 28.2   -0.642   -0.183  -3.685
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.4326 0.1020 24.7    0.223    0.642   4.250
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      -0.0462 0.1020 26.5   -0.255    0.162  -0.455
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.4788 0.0197 26.6   -0.519   -0.438 -24.325
##  p.value
##   <.0001
##   0.5546
##   0.0010
##   0.0003
##   0.6526
##   <.0001
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
##  diff_ESS10    0.366 0.0200 19.7    0.324    0.408  18.286  <.0001
##  diff_ESS1     0.479 0.0197 26.6    0.438    0.519  24.325  <.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
# Test for whether and how much the differences have changed
pairs(diff_mod7_log_GDP,infer=c(T,T))
```

```
##  contrast               estimate    SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10 - diff_ESS1   -0.113 0.022 21.1   -0.158  -0.0668  -5.115  <.0001
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
## Warning in checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv, : Model failed to converge
## with max|grad| = 0.00737069 (tol = 0.002, component 1)
```

``` r
getFE(mod8_log_GDP)
```

```
##                       Est.   SE    df     t     p    LL    UL
## (Intercept)           0.06 0.05 31.61  1.26 0.216 -0.04  0.16
## gndr.c                0.42 0.02 25.97 25.50 0.000  0.39  0.45
## log_gdp.z.cm         -0.09 0.04 33.59 -2.24 0.032 -0.17 -0.01
## essround.c           -0.02 0.01 27.78 -1.23 0.228 -0.04  0.01
## log_gdp.z.cmc         0.08 0.08 20.99  1.01 0.323 -0.08  0.24
## gndr.c:log_gdp.z.cm   0.02 0.02 34.56  1.26 0.217 -0.01  0.05
## gndr.c:essround.c    -0.01 0.00 11.35 -2.80 0.017 -0.02  0.00
## gndr.c:log_gdp.z.cmc -0.05 0.04  8.55 -1.19 0.265 -0.13  0.04
```

``` r
getVC(mod8_log_GDP)
```

```
##         grp                 var1                 var2 sdcor  vcov
## 1     cntry          (Intercept)                 <NA>  0.27  0.07
## 2     cntry               gndr.c                 <NA>  0.09  0.01
## 3     cntry           essround.c                 <NA>  0.07  0.01
## 4     cntry        log_gdp.z.cmc                 <NA>  0.42  0.18
## 5     cntry    gndr.c:essround.c                 <NA>  0.01  0.00
## 6     cntry gndr.c:log_gdp.z.cmc                 <NA>  0.14  0.02
## 7     cntry          (Intercept)               gndr.c -0.28 -0.01
## 8     cntry          (Intercept)           essround.c -0.08  0.00
## 9     cntry          (Intercept)        log_gdp.z.cmc  0.47  0.05
## 10    cntry          (Intercept)    gndr.c:essround.c  0.17  0.00
## 11    cntry          (Intercept) gndr.c:log_gdp.z.cmc -0.27 -0.01
## 12    cntry               gndr.c           essround.c -0.13  0.00
## 13    cntry               gndr.c        log_gdp.z.cmc -0.19 -0.01
## 14    cntry               gndr.c    gndr.c:essround.c -0.07  0.00
## 15    cntry               gndr.c gndr.c:log_gdp.z.cmc  0.05  0.00
## 16    cntry           essround.c        log_gdp.z.cmc -0.75 -0.02
## 17    cntry           essround.c    gndr.c:essround.c -0.73  0.00
## 18    cntry           essround.c gndr.c:log_gdp.z.cmc  0.38  0.00
## 19    cntry        log_gdp.z.cmc    gndr.c:essround.c  0.65  0.00
## 20    cntry        log_gdp.z.cmc gndr.c:log_gdp.z.cmc -0.59 -0.04
## 21    cntry    gndr.c:essround.c gndr.c:log_gdp.z.cmc -0.51  0.00
## 22 Residual                 <NA>                 <NA>  1.00  1.00
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
##              npar     AIC     BIC  logLik -2*log(L)   Chisq Df Pr(>Chisq)    
## mod2_log_GDP    8 1170508 1170595 -585246   1170492                          
## mod7_log_GDP   17 1166542 1166727 -583254   1166508 3984.46  9  < 2.2e-16 ***
## mod8_log_GDP   30 1166319 1166646 -583130   1166259  248.34 13  < 2.2e-16 ***
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
##         4.5 -0.203 0.0790 19.9  -0.3681  -0.0382  -2.570  0.0183
##        -4.5 -0.097 0.0840 22.1  -0.2713   0.0772  -1.155  0.2606
## 
## gndr.c =  0.5:
##  essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.176 0.0693 22.0   0.0323   0.3195   2.539  0.0187
##        -4.5  0.365 0.0759 24.3   0.2083   0.5214   4.807  0.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_mod8_log_GDP,adjust="none",infer=c(T,T))
```

```
## gndr.c = -0.5:
##  contrast                         estimate    SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   -0.106 0.128 27.4   -0.369    0.157  -0.827  0.4151
## 
## gndr.c =  0.5:
##  contrast                         estimate    SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   -0.189 0.112 26.1   -0.420    0.042  -1.682  0.1046
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
##    -0.5        4.5 -0.203 0.0790 19.9  -0.3681  -0.0382  -2.570  0.0183
##     0.5        4.5  0.176 0.0693 22.0   0.0323   0.3195   2.539  0.0187
##    -0.5       -4.5 -0.097 0.0840 22.1  -0.2713   0.0772  -1.155  0.2606
##     0.5       -4.5  0.365 0.0759 24.3   0.2083   0.5214   4.807  0.0001
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod8_log_GDP,adjust="none",infer=c(T,T))
```

```
##  contrast                                                 estimate     SE   df lower.CL upper.CL t.ratio
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5       -0.379 0.0213 17.8  -0.4238   -0.334 -17.772
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   -0.106 0.1280 27.4  -0.3690    0.157  -0.827
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)    -0.568 0.1230 26.9  -0.8200   -0.316  -4.626
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)       0.273 0.1190 27.6   0.0296    0.516   2.299
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)       -0.189 0.1120 26.1  -0.4200    0.042  -1.682
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)   -0.462 0.0230 12.8  -0.5117   -0.412 -20.085
##  p.value
##   <.0001
##   0.4151
##   0.0001
##   0.0293
##   0.1046
##   <.0001
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
##  diff_ESS10    0.379 0.0213 17.8    0.334    0.424  17.772  <.0001
##  diff_ESS1     0.462 0.0230 12.8    0.412    0.512  20.085  <.0001
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
##  diff_ESS10 - diff_ESS1  -0.0829 0.0297 11.3   -0.148  -0.0179  -2.796  0.0169
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
##  [1] apaTables_2.0.8       stringr_1.5.2         tidyr_1.3.1           r2mlm_0.3.8          
##  [5] nlme_3.1-168          Hmisc_5.2-4           ggpubr_0.6.2          metafor_4.8-0        
##  [9] numDeriv_2016.8-1.1   metadat_1.4-0         lmerTest_3.1-3        ggflags_0.0.4        
## [13] finalfit_1.1.0        ggplot2_4.0.0         MetBrewer_0.2.0       vjihelpers_0.0.0.9000
## [17] emmeans_1.11.2-8      lme4_1.1-37           Matrix_1.7-3          dplyr_1.1.4          
## [21] rio_1.2.4             multid_1.0.2.9000     knitr_1.50            rmarkdown_2.30       
## 
## loaded via a namespace (and not attached):
##   [1] mnormt_2.1.1       Rdpack_2.6.4       pROC_1.19.0.1      gridExtra_2.3      writexl_1.5.4     
##   [6] readxl_1.4.5       rlang_1.1.6        magrittr_2.0.4     rockchalk_1.8.157  compiler_4.5.1    
##  [11] mgcv_1.9-3         png_0.1-8          vctrs_0.6.5        quadprog_1.5-8     crayon_1.5.3      
##  [16] pkgconfig_2.0.3    shape_1.4.6.1      fastmap_1.2.0      backports_1.5.0    labeling_0.4.3    
##  [21] pbivnorm_0.6.0     utf8_1.2.6         nloptr_2.2.1       purrr_1.1.0        xfun_0.53         
##  [26] glmnet_4.1-10      jomo_2.7-6         cachem_1.1.0       kutils_1.73        jsonlite_2.0.0    
##  [31] pan_1.9            jpeg_0.1-11        psych_2.5.6        lavaan_0.6-20      parallel_4.5.1    
##  [36] broom_1.0.10       cluster_2.1.8.1    R6_2.6.1           bslib_0.9.0        stringi_1.8.7     
##  [41] RColorBrewer_1.1-3 car_3.1-3          boot_1.3-31        rpart_4.1.24       cellranger_1.1.0  
##  [46] jquerylib_0.1.4    estimability_1.5.1 Rcpp_1.1.0         iterators_1.0.14   base64enc_0.1-3   
##  [51] R.utils_2.13.0     splines_4.5.1      nnet_7.3-20        tidyselect_1.2.1   rstudioapi_0.17.1 
##  [56] abind_1.4-8        yaml_2.3.10        codetools_0.2-20   lattice_0.22-7     tibble_3.3.0      
##  [61] plyr_1.8.9         withr_3.0.2        S7_0.2.0           coda_0.19-4.1      evaluate_1.0.5    
##  [66] foreign_0.8-90     survival_3.8-3     zip_2.3.3          pillar_1.11.1      carData_3.0-5     
##  [71] mice_3.18.0        stats4_4.5.1       checkmate_2.3.3    foreach_1.5.2      reformulas_0.4.1  
##  [76] generics_0.1.4     grImport2_0.3-3    mathjaxr_1.8-0     scales_1.4.0       minqa_1.2.8       
##  [81] xtable_1.8-4       glue_1.8.0         tools_4.5.1        data.table_1.17.8  openxlsx_4.2.8    
##  [86] ggsignif_0.6.4     forcats_1.0.1      XML_3.99-0.19      mvtnorm_1.3-3      cowplot_1.2.0     
##  [91] grid_4.5.1         rbibutils_2.3      colorspace_2.1-2   htmlTable_2.4.3    Formula_1.2-5     
##  [96] cli_3.6.5          gtable_0.3.6       R.methodsS3_1.8.2  rstatix_0.7.3      sass_0.4.10       
## [101] digest_0.6.37      htmlwidgets_1.6.4  farver_2.1.2       htmltools_0.5.8.1  R.oo_1.27.1       
## [106] lifecycle_1.0.4    mitml_0.4-5        MASS_7.3-65
```

