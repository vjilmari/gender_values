---
title: "Data preparations for UN country indices"
output: 
  html_document: 
    toc: true
    keep_md: true
---





# Packages


``` r
library(rio)
library(dplyr)
library(psych)
```

# Read Data


``` r
UN_d<-read.csv("../data/HDR25_Composite_indices_complete_time_series.csv")

# read also ISO to obtain ISO2 codes
ISO<-read.csv2("../data/ISO.csv")
```

From Fors Connolly et al. 2020

The independent variable is Gender Equality Index (Gender Inequality Index reversed). The Gender Equality Index (GEI) was collected by the United Nations Development Programme (UNDP 2018). The index is a composite of health, empowerment, and labor market participation. GEI was only available for 2000 (1995 for the Czech Republic and Slovenia), 2005, and every year from 2010 and onward. We linearly imputed values for the missing years using the adjacent values. GEI changed relatively evenly within countries over years so the imputation procedure is not likely to misrepresent the true scores. All countries increased their gender equality over time, but to varying degrees. On a 0 to 1 scale (higher scores indicating higher gender equality), countries ranged from .43 (Turkey in 2004) and .96 (Denmark, Norway, Sweden, Switzerland in 2016). The GEI grand mean was .84 in 2002 and .92 in 2016. The within-country change in GEI between 2002 and 2016 was on average .06, ranging from .01 in Sweden to .15 in Poland.

## GII/GEI (Gender inequality index, reverse coded as gender equality index)


``` r
# variable names
GII_vars<-names(UN_d)[grepl("gii",names(UN_d)) & !grepl("gii_rank",names(UN_d))]

# select GII variables
GII_d<-UN_d %>%
  select(iso3,country,all_of(GII_vars))

# add ISO2 to GII_d
GII_d <- left_join(x=GII_d,
                   y=ISO[,c("ISO2","ISO3")],
                   by=c("iso3"="ISO3"))

# calculate GII mean across the time period 2002-2023 including the middle years

GII_vars_2002_2023<-
  paste0("gii_",2002:2023)
GII_vars_2002_2023
```

```
##  [1] "gii_2002" "gii_2003" "gii_2004" "gii_2005" "gii_2006" "gii_2007" "gii_2008" "gii_2009" "gii_2010"
## [10] "gii_2011" "gii_2012" "gii_2013" "gii_2014" "gii_2015" "gii_2016" "gii_2017" "gii_2018" "gii_2019"
## [19] "gii_2020" "gii_2021" "gii_2022" "gii_2023"
```

``` r
table(GII_vars_2002_2023 %in% names(GII_d))
```

```
## 
## TRUE 
##   22
```

``` r
GII_d$gii_2002_2023_avg<-
  rowMeans(GII_d[,GII_vars_2002_2023],na.rm=T)
  

describe(GII_d$gii_2002_2023_avg)
```

```
##    vars   n mean   sd median trimmed  mad  min  max range  skew kurtosis   se
## X1    1 183 0.39 0.19   0.41    0.39 0.23 0.04 0.82  0.79 -0.17    -1.09 0.01
```

``` r
# reverse code

GII_d$gei_2002_2023_avg<-(1-GII_d$gii_2002_2023_avg)
describe(GII_d$gei_2002_2023_avg)
```

```
##    vars   n mean   sd median trimmed  mad  min  max range skew kurtosis   se
## X1    1 183 0.61 0.19   0.59    0.61 0.23 0.18 0.96  0.79 0.17    -1.09 0.01
```

``` r
GII_d
```

```
##    iso3                          country gii_1990 gii_1991 gii_1992 gii_1993 gii_1994 gii_1995 gii_1996
## 1   AFG                      Afghanistan       NA       NA       NA       NA       NA       NA       NA
## 2   ALB                          Albania       NA       NA       NA       NA       NA       NA       NA
## 3   DZA                          Algeria       NA       NA       NA       NA       NA    0.673    0.661
## 4   AND                          Andorra       NA       NA       NA       NA       NA       NA       NA
## 5   AGO                           Angola       NA       NA       NA       NA       NA       NA       NA
## 6   ATG              Antigua and Barbuda       NA       NA       NA       NA       NA       NA       NA
## 7   ARG                        Argentina    0.489    0.452    0.450    0.442    0.444    0.442    0.429
## 8   ARM                          Armenia       NA       NA       NA       NA       NA    0.468    0.465
## 9   AUS                        Australia    0.188    0.187    0.187    0.183    0.181    0.178    0.177
## 10  AUT                          Austria    0.195    0.198    0.195    0.190    0.174    0.163    0.156
## 11  AZE                       Azerbaijan       NA       NA       NA    0.416    0.409    0.393    0.389
## 12  BHS                          Bahamas    0.374    0.371    0.378    0.379    0.370    0.367    0.370
## 13  BHR                          Bahrain       NA       NA       NA       NA       NA       NA       NA
## 14  BGD                       Bangladesh    0.700    0.696    0.690    0.688    0.685    0.683    0.682
## 15  BRB                         Barbados    0.385    0.386    0.388    0.372    0.367    0.364    0.364
## 16  BLR                          Belarus       NA       NA       NA       NA       NA       NA       NA
## 17  BEL                          Belgium    0.180    0.179    0.175    0.167    0.160    0.157    0.158
## 18  BLZ                           Belize       NA       NA       NA       NA    0.540    0.532    0.526
## 19  BEN                            Benin    0.700    0.697    0.687    0.680    0.673    0.665    0.655
## 20  BTN                           Bhutan       NA       NA       NA       NA       NA       NA       NA
## 21  BOL Bolivia (Plurinational State of)       NA       NA       NA       NA       NA       NA       NA
## 22  BIH           Bosnia and Herzegovina       NA       NA       NA       NA       NA       NA       NA
## 23  BWA                         Botswana    0.595    0.590    0.581    0.568    0.559    0.546    0.541
## 24  BRA                           Brazil    0.582    0.577    0.561    0.557    0.553    0.543    0.536
## 25  BRN                Brunei Darussalam       NA       NA       NA       NA       NA       NA       NA
##    gii_1997 gii_1998 gii_1999 gii_2000 gii_2001 gii_2002 gii_2003 gii_2004 gii_2005 gii_2006 gii_2007
## 1        NA       NA       NA       NA       NA       NA       NA       NA       NA       NA       NA
## 2        NA       NA    0.307    0.291    0.237    0.262    0.278    0.266    0.260    0.252    0.253
## 3     0.643    0.627    0.624    0.610    0.597    0.564    0.559    0.554    0.552    0.546    0.528
## 4        NA       NA       NA       NA       NA       NA       NA       NA       NA       NA       NA
## 5        NA       NA       NA       NA       NA       NA       NA       NA       NA       NA       NA
## 6        NA       NA       NA       NA       NA       NA       NA       NA       NA       NA       NA
## 7     0.429    0.425    0.421    0.426    0.425    0.396    0.397    0.390    0.375    0.375    0.374
## 8     0.447    0.439    0.470    0.452    0.447    0.446    0.414    0.406    0.397    0.399    0.357
## 9     0.174    0.171    0.158    0.153    0.152    0.146    0.141    0.137    0.137    0.135    0.136
## 10    0.148    0.143    0.141    0.142    0.141    0.140    0.124    0.121    0.115    0.112    0.108
## 11    0.389    0.383    0.375    0.358    0.358    0.356    0.353    0.367    0.354    0.358    0.348
## 12    0.365    0.356    0.370    0.360    0.360    0.361    0.363    0.349    0.352    0.357    0.368
## 13       NA       NA       NA       NA       NA       NA    0.353    0.357    0.360    0.291    0.266
## 14    0.670    0.659    0.650    0.648    0.657    0.728    0.745    0.739    0.654    0.646    0.645
## 15    0.364    0.361    0.357    0.368    0.368    0.357    0.359    0.349    0.342    0.338    0.336
## 16       NA       NA       NA       NA    0.255    0.224    0.207    0.155    0.153    0.150    0.151
## 17    0.156    0.156    0.128    0.132    0.135    0.129    0.108    0.103    0.101    0.100    0.098
## 18    0.523    0.527    0.514    0.508    0.505    0.492    0.509    0.493    0.483    0.468    0.445
## 19    0.652    0.652    0.658    0.655    0.657    0.668    0.659    0.652    0.650    0.649    0.628
## 20       NA       NA       NA       NA       NA       NA       NA       NA       NA       NA       NA
## 21       NA       NA    0.582    0.581    0.569    0.540    0.539    0.536    0.541    0.533    0.535
## 22       NA       NA       NA       NA       NA       NA       NA       NA       NA       NA       NA
## 23    0.538    0.540    0.538    0.500    0.499    0.493    0.493    0.517    0.517    0.519    0.518
## 24    0.530    0.522    0.514    0.506    0.495    0.494    0.471    0.470    0.467    0.466    0.463
## 25       NA       NA       NA       NA       NA       NA       NA       NA       NA       NA       NA
##    gii_2008 gii_2009 gii_2010 gii_2011 gii_2012 gii_2013 gii_2014 gii_2015 gii_2016 gii_2017 gii_2018
## 1     0.690    0.696    0.704    0.716    0.731    0.708    0.687    0.683    0.682    0.679    0.677
## 2     0.245    0.192    0.192    0.194    0.197    0.194    0.190    0.179    0.162    0.145    0.138
## 3     0.521    0.513    0.508    0.496    0.407    0.403    0.399    0.397    0.392    0.397    0.383
## 4        NA       NA       NA       NA       NA       NA       NA       NA       NA       NA       NA
## 5        NA       NA    0.556    0.549    0.550    0.551    0.553    0.549    0.544    0.547    0.539
## 6        NA       NA       NA       NA       NA       NA       NA       NA       NA       NA       NA
## 7     0.373    0.370    0.370    0.364    0.353    0.352    0.347    0.338    0.333    0.321    0.301
## 8     0.348    0.346    0.352    0.351    0.326    0.314    0.309    0.301    0.296    0.254    0.242
## 9     0.140    0.137    0.137    0.132    0.131    0.122    0.116    0.109    0.101    0.094    0.091
## 10    0.115    0.107    0.105    0.097    0.094    0.083    0.081    0.083    0.081    0.071    0.060
## 11    0.346    0.346    0.333    0.332    0.339    0.342    0.346    0.341    0.328    0.323    0.320
## 12    0.366    0.384    0.371    0.387    0.394    0.374    0.384    0.366    0.357    0.333    0.331
## 13    0.261    0.260    0.253    0.241    0.238    0.236    0.233    0.244    0.237    0.235    0.204
## 14    0.673    0.615    0.603    0.597    0.583    0.578    0.569    0.563    0.557    0.541    0.532
## 15    0.375    0.354    0.356    0.353    0.352    0.342    0.337    0.341    0.340    0.340    0.322
## 16    0.147    0.157    0.149    0.146    0.152    0.147    0.138    0.128    0.119    0.109    0.102
## 17    0.098    0.093    0.090    0.087    0.082    0.073    0.066    0.062    0.058    0.056    0.052
## 18    0.438    0.423    0.406    0.409    0.384    0.401    0.413    0.429    0.440    0.464    0.466
## 19    0.629    0.635    0.637    0.652    0.653    0.652    0.653    0.659    0.659    0.656    0.645
## 20       NA       NA    0.442    0.413    0.392    0.430    0.404    0.385    0.363    0.347    0.277
## 21    0.531    0.489    0.480    0.476    0.479    0.473    0.435    0.441    0.433    0.431    0.430
## 22       NA       NA       NA       NA       NA       NA       NA       NA    0.169    0.167    0.170
## 23    0.518    0.522    0.509    0.530    0.538    0.498    0.502    0.505    0.480    0.459    0.492
## 24    0.460    0.456    0.447    0.446    0.441    0.445    0.441    0.436    0.432    0.422    0.401
## 25       NA       NA       NA       NA       NA       NA       NA       NA    0.312    0.286    0.286
##    gii_2019 gii_2020 gii_2021 gii_2022 gii_2023 ISO2 gii_2002_2023_avg gei_2002_2023_avg
## 1     0.676    0.674    0.642    0.665    0.661   AF        0.68568750         0.3143125
## 2     0.131    0.125    0.109    0.107    0.107   AL        0.18990909         0.8100909
## 3     0.385    0.383    0.442    0.443    0.443   DZ        0.46431818         0.5356818
## 4        NA       NA       NA       NA       NA   AD               NaN               NaN
## 5     0.536    0.533    0.532    0.524    0.515   AO        0.54128571         0.4587143
## 6        NA       NA    0.225    0.224    0.240   AG        0.22966667         0.7703333
## 7     0.283    0.285    0.272    0.262    0.264   AR        0.34068182         0.6593182
## 8     0.216    0.207    0.185    0.180    0.180   AM        0.31027273         0.6897273
## 9     0.079    0.072    0.065    0.056    0.056   AU        0.11227273         0.8877273
## 10    0.052    0.048    0.042    0.038    0.033   AT        0.08681818         0.9131818
## 11    0.314    0.332    0.315    0.317    0.315   AZ        0.33750000         0.6625000
## 12    0.335    0.326    0.329    0.327    0.325   BS        0.35631818         0.6436818
## 13    0.191    0.189    0.191    0.177    0.165   BH        0.24676190         0.7532381
## 14    0.520    0.500    0.494    0.487    0.487   BD        0.59345455         0.4065455
## 15    0.316    0.302    0.303    0.298    0.297   BB        0.33677273         0.6632273
## 16    0.096    0.094    0.092    0.088    0.080   BY        0.13563636         0.8643636
## 17    0.047    0.044    0.042    0.038    0.031   BE        0.07536364         0.9246364
## 18    0.461    0.469    0.436    0.436    0.428   BZ        0.44513636         0.5548636
## 19    0.636    0.638    0.634    0.640    0.573   BJ        0.64350000         0.3565000
## 20    0.265    0.266    0.266    0.283    0.278   BT        0.34364286         0.6563571
## 21    0.429    0.431    0.423    0.421    0.419   BO        0.47477273         0.5252273
## 22    0.153    0.144    0.138    0.158    0.157   BA        0.15700000         0.8430000
## 23    0.466    0.497    0.495    0.492    0.490   BW        0.50227273         0.4977273
## 24    0.398    0.411    0.405    0.390    0.390   BR        0.43872727         0.5612727
## 25    0.286    0.276    0.272    0.276    0.257   BN        0.28137500         0.7186250
##  [ reached 'max' / getOption("max.print") -- omitted 181 rows ]
```

``` r
rio::export(data.frame(GII_d),"../data/GII.csv",overwrite=T)
```

## GDI (Gender development index)


``` r
# variable names
GDI_vars<-names(UN_d)[grepl("gdi",names(UN_d)) & !grepl("gdi_group",names(UN_d))]

# select GDI variables
GDI_d<-UN_d %>%
  select(iso3,country,all_of(GDI_vars))

# add ISO2 to GDI_d
GDI_d <- left_join(x=GDI_d,
                   y=ISO[,c("ISO2","ISO3")],
                   by=c("iso3"="ISO3"))

# calculate GDI mean across the time period 2002-2023 including the middle years

GDI_vars_2002_2023<-
  paste0("gdi_",2002:2023)
GDI_vars_2002_2023
```

```
##  [1] "gdi_2002" "gdi_2003" "gdi_2004" "gdi_2005" "gdi_2006" "gdi_2007" "gdi_2008" "gdi_2009" "gdi_2010"
## [10] "gdi_2011" "gdi_2012" "gdi_2013" "gdi_2014" "gdi_2015" "gdi_2016" "gdi_2017" "gdi_2018" "gdi_2019"
## [19] "gdi_2020" "gdi_2021" "gdi_2022" "gdi_2023"
```

``` r
table(GDI_vars_2002_2023 %in% names(GDI_d))
```

```
## 
## TRUE 
##   22
```

``` r
GDI_d$gdi_2002_2023_avg<-
  rowMeans(GDI_d[,GDI_vars_2002_2023],na.rm=T)


describe(GDI_d$gdi_2002_2023_avg)
```

```
##    vars   n mean   sd median trimmed  mad  min  max range  skew kurtosis se
## X1    1 195 0.94 0.07   0.96    0.95 0.05 0.55 1.04  0.49 -1.64     4.57  0
```

``` r
rio::export(data.frame(GDI_d),"../data/GDI.csv",overwrite=T)
```

