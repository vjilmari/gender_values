---
title: "Analysis for Examining the Gender Equality Paradox in Values Using Universalism Value"
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

# Data for the analysis

* basic descriptives


``` r
diff_dat<-
  fdat

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
##   1  2254 1830    0 2024    0 1208 2819 1470    0 1712 1763 1355 1798 2551    0 1634 1916 2279    0    0
##   2  2198 1771    0 2110    0 2557 2840 1458 1948 1623 1701 1699 1864 2399    0 1460 1187    0  524    0
##   3  2348 1796 1295 1780  978    0 2884 1461 1466 1847 1649 1983 2353    0    0 1462 1589    0    0    0
##   4     0 1754 2144 1753 1210 1986 2732 1581 1646 2562 1901 2067 2311 2063 1430 1430 1757 2382    0    0
##   5     0 1699 2371 1491 1053 2335 3007 1564 1793 1881 1649 1723 2374 2669 1601 1473 2400 2212    0    0
##   6     0 1862 2179 1483 1110 1973 2935 1621 2345 1871 2158 1960 2261    0    0 1968 2616 2378  739  909
##   7  1795 1767    0 1521    0 1862 3006 1483 2036 1907 2050 1902 2231    0    0 1520 2380 2351    0    0
##   8  1993 1759    0 1504    0 2252 2821    0 2007 1929 1903 2057 1942    0    0 1458 2746 2366  841 2531
##   9  2477 1756 1926 1517  773 2343 2328 1554 1899 1619 1735 1982 2183    0 1781 1643 2189    0  844 2660
##   10    0 1334 2697 1505    0 2369    0    0 1538    0 1561 1951 1131 2768 1564 1816 1751    0  886 2573
##     
##        LT   LV   ME   NL   NO   PL   PT   RU   SE   SI   SK   TR   UA
##   1     0    0    0 2337 1819 2065 1482    0 1682 1488    0    0    0
##   2     0    0    0 1858 1575 1683 2024    0 1678 1384 1425 1790 1896
##   3     0    0    0 1860 1550 1685 2182 2339 1604 1465 1711    0 1885
##   4     0 1970    0 1724 1391 1596 2337 2446 1556 1257 1789 2305 1766
##   5  1632    0    0 1801 1530 1719 2139 2557 1463 1369 1803    0 1779
##   6  2108    0    0 1828 1610 1866 2138 2429 1838 1244 1827    0 2064
##   7  2241    0    0 1823 1423 1594 1242    0 1761 1189    0    0    0
##   8  2079    0    0 1669 1530 1675 1254 2374 1526 1295    0    0    0
##   9  1677  891 1188 1657 1396 1443 1045    0 1510 1307 1061    0    0
##   10 1606    0 1248 1466 1408    0 1827    0    0 1232 1395    0    0
```

``` r
# range of sample sizes
range(table(diff_dat$cntry))
```

```
## [1]  2436 25372
```

``` r
# scale/standardize with SD pooled across all country-time-gender folds
cntry.uni<-diff_dat %>% group_by(cntry,essround) %>%
  summarise(uni.ctm=mean(uni,na.rm=T),
            uni.ctsd=sd(uni,na.rm=T)) %>%
  group_by(cntry) %>%
  summarise(uni.cm=mean(uni.ctm),
            uni.csd=mean(uni.ctsd)) 
```

```
## `summarise()` has grouped output by 'cntry'. You can override using the `.groups` argument.
```

``` r
grand_mean_uni<-mean(cntry.uni$uni.cm)
grand_sd_uni<-mean(cntry.uni$uni.csd)

# standardized
diff_dat$uni.z<-(diff_dat$uni-grand_mean_uni)/grand_sd_uni
hist(diff_dat$uni.z)
```

![](Analysis_universalism_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

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

# value-based universalism

cntry_uni_frame<-
  diff_dat %>% group_by(cntry,essround) %>%
  summarise('uni M' = weighted.mean(x=uni.z,w=pspwght),
            'uni SD' = sqrt(wtd.var(uni.z,pspwght))) %>%
  group_by(cntry) %>%
  summarise('uni M' = mean(x=`uni M`),
            'uni SD'= mean(x=`uni SD`))
```

```
## `summarise()` has grouped output by 'cntry'. You can override using the `.groups` argument.
```

``` r
cntry_uni_women_frame<-
  diff_dat %>%
  filter(gndr.c==-0.5) %>%
  group_by(cntry,essround) %>%
  summarise('uni M' = weighted.mean(x=uni.z,w=pspwght),
            'uni SD' = sqrt(wtd.var(uni.z,pspwght))) %>%
  group_by(cntry) %>%
  summarise('uni M Women' = mean(x=`uni M`),
            'uni SD Women'= mean(x=`uni SD`))
```

```
## `summarise()` has grouped output by 'cntry'. You can override using the `.groups` argument.
```

``` r
cntry_uni_men_frame<-
  diff_dat %>%
  filter(gndr.c==0.5) %>%
  group_by(cntry,essround) %>%
  summarise('uni M' = weighted.mean(x=uni.z,w=pspwght),
            'uni SD' = sqrt(wtd.var(uni.z,pspwght))) %>%
  group_by(cntry) %>%
  summarise('uni M Men' = mean(x=`uni M`),
            'uni SD Men'= mean(x=`uni SD`))
```

```
## `summarise()` has grouped output by 'cntry'. You can override using the `.groups` argument.
```

``` r
# link n and uni datasets

desc_frame<-
  left_join(
    x=cntry_n_frame,
    y=cntry_uni_frame,
    by="cntry"
  )

desc_frame<-
  left_join(
    x=desc_frame,
    y=cntry_uni_women_frame,
    by="cntry"
  )

desc_frame<-
  left_join(
    x=desc_frame,
    y=cntry_uni_men_frame,
    by="cntry"
  )

# Add country-specific differences
desc_frame$D<-desc_frame$`uni M Men`-desc_frame$`uni M Women`

desc_frame
```

```
## # A tibble: 33 × 10
##    cntry `n ESS rounds`     n `uni M` `uni SD` `uni M Women` `uni SD Women` `uni M Men` `uni SD Men`
##    <chr>          <dbl> <dbl>   <dbl>    <dbl>         <dbl>          <dbl>       <dbl>        <dbl>
##  1 AT                 6 13077  0.0859    1.04         0.199           1.00     -0.0353         1.06 
##  2 BE                10 17313  0.111     0.864        0.184           0.843     0.0355         0.878
##  3 BG                 6 12641 -0.153     1.09        -0.0817          1.06     -0.230          1.12 
##  4 CH                10 16720  0.316     0.865        0.412           0.834     0.215          0.885
##  5 CY                 5  5105  0.341     0.862        0.352           0.845     0.330          0.879
##  6 CZ                 9 18934 -0.406     1.07        -0.282           1.05     -0.541          1.08 
##  7 DE                 9 25389  0.107     0.935        0.210           0.895    -0.00174        0.964
##  8 DK                 8 12198 -0.0311    1.02         0.0642          0.985    -0.129          1.05 
##  9 EE                 9 16692 -0.126     0.928       -0.0147          0.896    -0.259          0.947
## 10 ES                 9 16954  0.376     0.901        0.417           0.884     0.333          0.915
## # ℹ 23 more rows
## # ℹ 1 more variable: D <dbl>
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
    `uni M`, `uni SD`,
    `uni M Women`, `uni SD Women`,
    `uni M Men`, `uni SD Men`,
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
print(cntry_desc_tbl,n=33)
```

```
## # A tibble: 33 × 14
##    Country    `n ESS rounds`     n `uni M` `uni SD` `uni M Women` `uni SD Women` `uni M Men` `uni SD Men`
##    <chr>               <dbl> <dbl> <chr>   <chr>    <chr>         <chr>          <chr>       <chr>       
##  1 Austria                 6 13077 0.09    1.04     0.20          1.00           -0.04       1.06        
##  2 Belgium                10 17313 0.11    0.86     0.18          0.84           0.04        0.88        
##  3 Bulgaria                6 12641 -0.15   1.09     -0.08         1.06           -0.23       1.12        
##  4 Switzerla…             10 16720 0.32    0.87     0.41          0.83           0.21        0.89        
##  5 Cyprus                  5  5105 0.34    0.86     0.35          0.84           0.33        0.88        
##  6 Czechia                 9 18934 -0.41   1.07     -0.28         1.05           -0.54       1.08        
##  7 Germany                 9 25389 0.11    0.94     0.21          0.89           -0.00       0.96        
##  8 Denmark                 8 12198 -0.03   1.02     0.06          0.98           -0.13       1.05        
##  9 Estonia                 9 16692 -0.13   0.93     -0.01         0.90           -0.26       0.95        
## 10 Spain                   9 16954 0.38    0.90     0.42          0.88           0.33        0.92        
## 11 Finland                10 18050 0.16    0.95     0.35          0.88           -0.04       0.98        
## 12 France                 10 18720 0.10    1.08     0.15          1.07           0.04        1.08        
## 13 UK                     10 21456 -0.00   1.02     0.06          1.00           -0.07       1.03        
## 14 Greece                  5 12464 0.26    0.93     0.25          0.92           0.28        0.94        
## 15 Croatia                 4  6368 0.00    1.08     0.10          1.04           -0.10       1.11        
## 16 Hungary                10 16006 0.01    1.01     0.08          0.99           -0.08       1.03        
## 17 Ireland                10 20576 0.04    1.05     0.12          1.04           -0.04       1.06        
## 18 Israel                  6 13964 0.10    1.06     0.14          1.06           0.06        1.06        
## 19 Iceland                 5  3832 0.07    0.98     0.19          0.94           -0.05       1.01        
## 20 Italy                   4  8663 0.07    0.98     0.12          0.98           0.02        0.99        
## 21 Lithuania               6 11714 -0.71   1.17     -0.62         1.16           -0.82       1.18        
## 22 Latvia                  2  2866 -0.25   1.06     -0.08         1.03           -0.45       1.05        
## 23 Montenegro              2  2441 -0.06   1.19     0.03          1.16           -0.15       1.21        
## 24 Netherlan…             10 18048 -0.03   0.88     0.06          0.85           -0.12       0.90        
## 25 Norway                 10 15186 -0.22   0.99     -0.13         0.96           -0.31       1.01        
## 26 Poland                  9 15314 0.07    0.91     0.14          0.88           -0.01       0.94        
## 27 Portugal               10 17705 -0.27   1.04     -0.25         1.03           -0.28       1.04        
## 28 Russia                  5 12139 -0.19   1.10     -0.13         1.09           -0.26       1.11        
## 29 Sweden                  9 14897 -0.04   0.98     0.10          0.94           -0.20       0.99        
## 30 Slovenia               10 13238 0.28    0.84     0.36          0.82           0.20        0.86        
## 31 Slovakia                7 11132 -0.15   0.95     -0.07         0.94           -0.23       0.96        
## 32 Turkey                  2  4108 0.18    0.96     0.16          0.97           0.20        0.96        
## 33 Ukraine                 5  9454 -0.32   1.24     -0.27         1.25           -0.39       1.23        
## # ℹ 5 more variables: D <chr>, GEI <chr>, GGGI <chr>, GDI <chr>, GDP <chr>
```

``` r
export(cntry_desc_tbl,"../results/uni/cntry_desc_tbl.xlsx",overwrite=T)
```

## Country-level correlations table


``` r
cor_frame<-
  desc_frame %>%
  select(
    VBMT=`uni M`,
    VBMT_Women=`uni M Women`,
    VBMT_Men=`uni M Men`,
    D = D,
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
  filename = "../results/uni/CorTable1.doc",
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
##   Variable      M     SD   1            2            3            4            5           6          
##   1. VBMT       -0.01 0.23                                                                            
##                                                                                                       
##   2. VBMT_Women 0.07  0.22 .98                                                                        
##                            [.96, .99]                                                                 
##                                                                                                       
##   3. VBMT_Men   -0.09 0.25 .98          .92                                                           
##                            [.96, .99]   [.85, .96]                                                    
##                                                                                                       
##   4. D          -0.16 0.10 .28          .09          .46                                              
##                            [-.07, .57]  [-.26, .42]  [.14, .70]                                       
##                                                                                                       
##   5. GEI        0.87  0.07 .12          .22          .03          -.40                                
##                            [-.24, .45]  [-.14, .53]  [-.32, .38]  [-.66, -.06]                        
##                                                                                                       
##   6. GGGI       0.73  0.05 -.00         .13          -.12         -.61         .73                    
##                            [-.35, .34]  [-.22, .45]  [-.45, .23]  [-.79, -.34] [.52, .86]             
##                                                                                                       
##   7. GDI        0.99  0.03 -.54         -.48         -.59         -.42         .07         .20        
##                            [-.75, -.25] [-.71, -.17] [-.78, -.31] [-.66, -.09] [-.29, .41] [-.16, .51]
##                                                                                                       
##   8. log_GDP    10.62 0.40 .30          .37          .24          -.25         .75         .67        
##                            [-.04, .59]  [.03, .63]   [-.12, .54]  [-.54, .11]  [.55, .87]  [.42, .82] 
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
##   -.22       
##   [-.53, .13]
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
mod0<-lmer(uni.z~(1|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod0)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: uni.z ~ (1 | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1314597.5 1314630.5 -657295.8 1314591.5    441165 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -8.3900 -0.5476  0.1251  0.6611  5.0698 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  cntry    (Intercept) 0.05165  0.2273  
##  Residual             1.01552  1.0077  
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##             Estimate Std. Error       df t value Pr(>|t|)
## (Intercept) -0.01322    0.03960 33.00819  -0.334    0.741
```

``` r
getVC(mod0)
```

```
##        grp        var1 var2 sdcor vcov
## 1    cntry (Intercept) <NA>  0.23 0.05
## 2 Residual        <NA> <NA>  1.01 1.02
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
## mean variation  0.04840191     NA       1
## sigma2          0.95159809      1      NA
## 
## $R2s
##          total within between
## f1  0.00000000      0      NA
## f2  0.00000000     NA       0
## v   0.00000000      0      NA
## m   0.04840191     NA       1
## f   0.00000000     NA      NA
## fv  0.00000000      0      NA
## fvm 0.04840191     NA      NA
```

## mod1: Gender fixed effect


``` r
mod1<-lmer(uni.z~gndr.c+(1|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod1)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: uni.z ~ gndr.c + (1 | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1311682.5 1311726.4 -655837.2 1311674.5    441164 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -8.4287 -0.5515  0.0997  0.6552  4.9231 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  cntry    (Intercept) 0.0523   0.2287  
##  Residual             1.0088   1.0044  
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##               Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept) -1.634e-02  3.985e-02  3.301e+01   -0.41    0.684    
## gndr.c      -1.634e-01  3.020e-03  4.411e+05  -54.10   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##        (Intr)
## gndr.c 0.001
```

``` r
getFE(mod1,round=3)
```

```
##               Est.    SE         df       t     p     LL     UL
## (Intercept) -0.016 0.040     33.008  -0.410 0.684 -0.097  0.065
## gndr.c      -0.163 0.003 441136.096 -54.099 0.000 -0.169 -0.157
```

``` r
getVC(mod1)
```

```
##        grp        var1 var2 sdcor vcov
## 1    cntry (Intercept) <NA>  0.23 0.05
## 2 Residual        <NA> <NA>  1.00 1.01
```

``` r
r2mlm(mod1,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.006212122
## slope variation 0.000000000
## mean variation  0.048982630
## sigma2          0.944805248
## 
## $R2s
##           total
## f   0.006212122
## v   0.000000000
## m   0.048982630
## fv  0.006212122
## fvm 0.055194752
```

## mod2: Gender fixed and random effect

* Include random effect correlation by default


``` r
mod2<-lmer(uni.z~gndr.c+(gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: uni.z ~ gndr.c + (gndr.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1310992.2 1311058.1 -655490.1 1310980.2    441162 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -8.3968 -0.5510  0.1023  0.6578  4.8962 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr
##  cntry    (Intercept) 0.052550 0.22924      
##           gndr.c      0.008628 0.09289  0.28
##  Residual             1.007000 1.00349      
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##             Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept) -0.01654    0.03995 33.00841  -0.414    0.681    
## gndr.c      -0.16060    0.01655 31.81259  -9.704    5e-11 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##        (Intr)
## gndr.c 0.273
```

``` r
getFE(mod2,round=3)
```

```
##               Est.    SE     df      t     p     LL     UL
## (Intercept) -0.017 0.040 33.008 -0.414 0.681 -0.098  0.065
## gndr.c      -0.161 0.017 31.813 -9.704 0.000 -0.194 -0.127
```

``` r
getVC(mod2)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.23 0.05
## 2    cntry      gndr.c   <NA>  0.09 0.01
## 3    cntry (Intercept) gndr.c  0.28 0.01
## 4 Residual        <NA>   <NA>  1.00 1.01
```

``` r
r2mlm(mod2,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.006002895
## slope variation 0.002008170
## mean variation  0.048798874
## sigma2          0.943190061
## 
## $R2s
##           total
## f   0.006002895
## v   0.002008170
## m   0.048798874
## fv  0.008011065
## fvm 0.056809939
```

``` r
anova(mod1,mod2)
```

```
## Data: diff_dat
## Models:
## mod1: uni.z ~ gndr.c + (1 | cntry)
## mod2: uni.z ~ gndr.c + (gndr.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod1    4 1311682 1311726 -655837   1311674                         
## mod2    6 1310992 1311058 -655490   1310980 694.29  2  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

``` r
lvl2_var_cond_lvl1(mod2,lvl1.var = "gndr.c",lvl1.values = c(0.5,-0.5))
```

```
##   lvl1.value lvl2.cond.var lvl2.cond.sd
## 1        0.5    0.06065795    0.2462883
## 2       -0.5    0.04875619    0.2208080
```

* Test for random effect correlation


``` r
mod2_norecov<-lmer(uni.z~gndr.c+(gndr.c||cntry),data=diff_dat,REML=F,weights = pspwght,
                   control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2_norecov)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: uni.z ~ gndr.c + (gndr.c || cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1310992.7 1311047.7 -655491.4 1310982.7    441163 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -8.3948 -0.5510  0.1023  0.6578  4.8995 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  cntry    (Intercept) 0.052519 0.22917 
##  cntry.1  gndr.c      0.008583 0.09264 
##  Residual             1.007001 1.00349 
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##             Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept) -0.01652    0.03993 33.00669  -0.414    0.682    
## gndr.c      -0.16058    0.01651 31.81631  -9.727 4.72e-11 ***
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
##               Est.    SE     df      t     p     LL     UL
## (Intercept) -0.017 0.040 33.007 -0.414 0.682 -0.098  0.065
## gndr.c      -0.161 0.017 31.816 -9.727 0.000 -0.194 -0.127
```

``` r
getVC(mod2_norecov)
```

```
##        grp        var1 var2 sdcor vcov
## 1    cntry (Intercept) <NA>  0.23 0.05
## 2  cntry.1      gndr.c <NA>  0.09 0.01
## 3 Residual        <NA> <NA>  1.00 1.01
```

``` r
anova(mod2_norecov,mod2)
```

```
## Data: diff_dat
## Models:
## mod2_norecov: uni.z ~ gndr.c + (gndr.c || cntry)
## mod2: uni.z ~ gndr.c + (gndr.c | cntry)
##              npar     AIC     BIC  logLik -2*log(L) Chisq Df Pr(>Chisq)
## mod2_norecov    5 1310993 1311048 -655491   1310983                    
## mod2            6 1310992 1311058 -655490   1310980 2.571  1     0.1088
```


## mod2 with Gender-equality index (GEI)


``` r
mod2_GEI<-lmer(uni.z~gndr.c+gei.z.cm+gndr.c:gei.z.cm+
                 (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2_GEI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: uni.z ~ gndr.c + gei.z.cm + gndr.c:gei.z.cm + (gndr.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1275913.6 1276001.4 -637948.8 1275897.6    431770 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -8.4101 -0.5527  0.1026  0.6600  4.9238 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr
##  cntry    (Intercept) 0.050185 0.22402      
##           gndr.c      0.007471 0.08644  0.39
##  Residual             0.995011 0.99750      
## Number of obs: 431778, groups:  cntry, 32
## 
## Fixed effects:
##                  Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)     -0.006685   0.039643 32.003715  -0.169   0.8672    
## gndr.c          -0.161227   0.015688 30.704524 -10.277 1.85e-11 ***
## gei.z.cm         0.028437   0.040295 32.059068   0.706   0.4854    
## gndr.c:gei.z.cm -0.038091   0.016111 32.027798  -2.364   0.0243 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c g.z.cm
## gndr.c       0.384              
## gei.z.cm    -0.001  0.000       
## gndr.c:g.z.  0.000 -0.016  0.380
```

``` r
getFE(mod2_GEI,round=3)
```

```
##                   Est.    SE     df       t     p     LL     UL
## (Intercept)     -0.007 0.040 32.004  -0.169 0.867 -0.087  0.074
## gndr.c          -0.161 0.016 30.705 -10.277 0.000 -0.193 -0.129
## gei.z.cm         0.028 0.040 32.059   0.706 0.485 -0.054  0.111
## gndr.c:gei.z.cm -0.038 0.016 32.028  -2.364 0.024 -0.071 -0.005
```

``` r
getVC(mod2_GEI)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.22 0.05
## 2    cntry      gndr.c   <NA>  0.09 0.01
## 3    cntry (Intercept) gndr.c  0.39 0.01
## 4 Residual        <NA>   <NA>  1.00 1.00
```

``` r
r2mlm(mod2_GEI,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.007334658
## slope variation 0.001762039
## mean variation  0.047078422
## sigma2          0.943824880
## 
## $R2s
##           total
## f   0.007334658
## v   0.001762039
## m   0.047078422
## fv  0.009096697
## fvm 0.056175120
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
## Time difference of 33.91011 secs
```

``` r
round(ddsc_mod2_GEI$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.049        0.221        1.007     1.056 0.046   7204.152 0.997   0.997
## 2        0.5         0.061        0.246        1.007     1.068 0.057   6164.576 0.997   0.997
```

``` r
round(ddsc_mod2_GEI$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gei.z.cm
## means_y1          -0.088 0.250    1.000           1.000    0.918           0.918    0.039
## means_y1_scaled   -0.374 1.063    1.000           1.000    0.918           0.918    0.039
## means_y2           0.074 0.219    0.918           0.918    1.000           1.000    0.215
## means_y2_scaled    0.317 0.932    0.918           0.918    1.000           1.000    0.215
## gei.z.cm           0.000 1.000    0.039           0.039    0.215           0.215    1.000
## gei.z.cm_scaled    0.000 1.000    0.039           0.039    0.215           0.215    1.000
## diff_score        -0.162 0.099    0.489           0.489    0.103           0.103   -0.377
## diff_score_scaled -0.691 0.423    0.489           0.489    0.103           0.103   -0.377
##                   gei.z.cm_scaled diff_score diff_score_scaled
## means_y1                    0.039      0.489             0.489
## means_y1_scaled             0.039      0.489             0.489
## means_y2                    0.215      0.103             0.103
## means_y2_scaled             0.215      0.103             0.103
## gei.z.cm                    1.000     -0.377            -0.377
## gei.z.cm_scaled             1.000     -0.377            -0.377
## diff_score                 -0.377      1.000             1.000
## diff_score_scaled          -0.377      1.000             1.000
```

``` r
round(ddsc_mod2_GEI$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                       0.383 0.162 32.028   2.364   0.024    0.053    0.713
## w_11                          0.047 0.038 32.162   1.250   0.220   -0.030    0.125
## w_21                          0.009 0.044 32.067   0.214   0.832   -0.080    0.099
## r_xy1                         0.190 0.152 32.162   1.250   0.220   -0.120    0.500
## r_xy2                         0.043 0.201 32.067   0.214   0.832   -0.366    0.452
## b_11                          0.203 0.162 32.162   1.250   0.220   -0.127    0.532
## b_21                          0.040 0.188 32.067   0.214   0.832   -0.342    0.422
## main_effect                   0.028 0.040 32.059   0.706   0.485   -0.054    0.111
## moderator_effect             -0.161 0.016 30.705 -10.277   0.000   -0.193   -0.129
## interaction                  -0.038 0.016 32.028  -2.364   0.024   -0.071   -0.005
## q_b11_b21                     0.165    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.150    NA     NA      NA      NA       NA       NA
## cross_over_point             -4.233    NA     NA      NA      NA       NA       NA
## interaction_vs_main           0.010 0.049 32.094   0.198   0.844   -0.090    0.109
## interaction_vs_main_bscale    0.041 0.208 32.094   0.198   0.844   -0.382    0.465
## interaction_vs_main_rscale    0.031 0.232 32.090   0.133   0.895   -0.442    0.503
## dadas                        -0.019 0.088 32.067  -0.214   0.584   -0.198    0.160
## dadas_bscale                 -0.080 0.375 32.067  -0.214   0.584   -0.845    0.684
## dadas_rscale                 -0.086 0.402 32.067  -0.214   0.584   -0.904    0.732
## abs_diff                      0.038 0.016 32.028   2.364   0.012    0.005    0.071
## abs_sum                       0.057 0.081 32.059   0.706   0.243   -0.107    0.221
## abs_diff_bscale               0.162 0.069 32.028   2.364   0.012    0.023    0.302
## abs_sum_bscale                0.243 0.344 32.059   0.706   0.243   -0.458    0.943
## abs_diff_rscale               0.147 0.080 31.961   1.832   0.038   -0.017    0.311
## abs_sum_rscale                0.233 0.347 32.057   0.671   0.253   -0.474    0.940
```

``` r
round(ddsc_mod2_GEI$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
##  0.006  0.279  2.571  1.000  0.109
```

``` r
d_GEI<-ddsc_mod2_GEI$ddsc_sem_fit$data

ddsc_sem_GEI<-
  ddsc_sem(data=d_GEI,x = "gei.z.cm",y1="means_y2",
           y2="means_y1",q_sesoi = .10)
round(ddsc_sem_GEI$results,3)
```

```
##                                    est    se      z pvalue ci.lower ci.upper
## r_xy1_y2                         0.377 0.164  2.304  0.021    0.056    0.698
## r_xy1                            0.215 0.173  1.247  0.212   -0.123    0.554
## r_xy2                            0.039 0.177  0.219  0.827   -0.308    0.385
## b_11                             0.201 0.161  1.247  0.212   -0.115    0.516
## b_21                             0.041 0.188  0.219  0.827   -0.327    0.409
## b_10                             0.317 0.158  2.001  0.045    0.006    0.628
## b_20                            -0.374 0.185 -2.023  0.043   -0.736   -0.012
## res_cov_y1_y2                    0.874 0.227  3.858  0.000    0.430    1.318
## diff_b10_b20                     0.691 0.068 10.130  0.000    0.557    0.825
## diff_b11_b21                     0.160 0.069  2.304  0.021    0.024    0.296
## diff_rxy1_rxy2                   0.177 0.064  2.747  0.006    0.051    0.303
## q_b11_b21                        0.162 0.068  2.375  0.018    0.028    0.296
## q_rxy1_rxy2                      0.180 0.066  2.732  0.006    0.051    0.309
## cross_over_point                -4.328 1.927 -2.246  0.025   -8.104   -0.552
## sum_b11_b21                      0.242 0.343  0.706  0.480   -0.430    0.914
## main_effect                      0.121 0.171  0.706  0.480   -0.215    0.457
## interaction_vs_main_effect       0.039 0.209  0.186  0.853   -0.370    0.448
## diff_abs_b11_abs_b21             0.160 0.069  2.304  0.021    0.024    0.296
## abs_diff_b11_b21                 0.160 0.069  2.304  0.011    0.024    0.296
## abs_sum_b11_b21                  0.242 0.343  0.706  0.240   -0.430    0.914
## dadas                           -0.082 0.376 -0.219  0.587   -0.818    0.654
## q_r_equivalence                  0.080 0.066  1.215  0.888       NA       NA
## q_b_equivalence                  0.062 0.068  0.913  0.819       NA       NA
## cross_over_point_equivalence     4.328 1.927  2.246  0.988       NA       NA
## cross_over_point_minimal_effect  4.328 1.927  2.246  0.012       NA       NA
```

``` r
round(ddsc_sem_GEI$variance_test,3)
```

```
##              est    se      z pvalue ci.lower ci.upper
## cov_y1y2   0.882 0.231  3.826  0.000    0.430    1.334
## var_y1     0.842 0.211  4.000  0.000    0.430    1.255
## var_y2     1.095 0.274  4.000  0.000    0.559    1.632
## var_diff  -0.253 0.149 -1.701  0.089   -0.544    0.038
## var_ratio  0.769 0.108  7.142  0.000    0.558    0.980
## cor_y1y2   0.918 0.028 33.116  0.000    0.864    0.973
```

### Figure 


``` r
# plot the results

# start with obtaining predicted values for means and differences

# refit reduced and full models with GGGI in original scale


mod2_GEI_unstd<-lmer(uni.z~gndr.c+gei.cm+gndr.c:gei.cm+
                 (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

mod2_GEI_unstd_red<-lmer(uni.z~gndr.c+
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


p1.uni.flags<-
  ggplot(p,aes(y=yvar,x=gei.cm,color=gndr.c))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("Gender Equality Index (Average 2002-2022)")+
  #ylim(c(min.y.pred,max.y.pred))+
  #xlim(c(0.60,1.00))+
  ylab("Value universalism (Average 2002-2022)")+
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

p2.uni.flags<-ggplot(p2,aes(y=yvar,x=gei.cm))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("Gender Equality Index (Average 2002-2022)")+
  #xlim(c(0.60,1.00))+
  #ylim(c(min.y.narrow,max.y.narrow))+
  ylab("Difference in value universalism")+
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
  ggarrange(p1.uni.flags,p2.uni.flags,align = "v",
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

![](Analysis_universalism_files/figure-html/unnamed-chunk-18-1.png)<!-- -->

``` r
png(filename = 
      "../results/uni/GEI_flags.png",
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
mod2_GGGI<-lmer(uni.z~gndr.c+gggi.z.cm+gndr.c:gggi.z.cm+
                 (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2_GGGI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: uni.z ~ gndr.c + gggi.z.cm + gndr.c:gggi.z.cm + (gndr.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##  935031.9  935117.2 -467507.9  935015.9    314638 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -8.4377 -0.5514  0.1081  0.6612  4.8389 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr
##  cntry    (Intercept) 0.052598 0.2293       
##           gndr.c      0.006116 0.0782   0.37
##  Residual             1.005355 1.0027       
## Number of obs: 314646, groups:  cntry, 33
## 
## Fixed effects:
##                  Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)      -0.01113    0.03998 33.01938  -0.278 0.782399    
## gndr.c           -0.15751    0.01425 30.09696 -11.055 4.04e-12 ***
## gggi.z.cm         0.01254    0.04063 33.10244   0.309 0.759580    
## gndr.c:gggi.z.cm -0.05885    0.01474 32.25842  -3.993 0.000354 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c ggg.z.
## gndr.c       0.355              
## gggi.z.cm    0.000  0.000       
## gndr.c:gg..  0.000 -0.011  0.348
```

``` r
getFE(mod2_GGGI,round=3)
```

```
##                    Est.    SE     df       t     p     LL     UL
## (Intercept)      -0.011 0.040 33.019  -0.278 0.782 -0.092  0.070
## gndr.c           -0.158 0.014 30.097 -11.055 0.000 -0.187 -0.128
## gggi.z.cm         0.013 0.041 33.102   0.309 0.760 -0.070  0.095
## gndr.c:gggi.z.cm -0.059 0.015 32.258  -3.993 0.000 -0.089 -0.029
```

``` r
getVC(mod2_GGGI)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.23 0.05
## 2    cntry      gndr.c   <NA>  0.08 0.01
## 3    cntry (Intercept) gndr.c  0.37 0.01
## 4 Residual        <NA>   <NA>  1.00 1.01
```

``` r
r2mlm(mod2_GGGI,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.006920507
## slope variation 0.001424606
## mean variation  0.048835888
## sigma2          0.942819000
## 
## $R2s
##           total
## f   0.006920507
## v   0.001424606
## m   0.048835888
## fv  0.008345112
## fvm 0.057181000
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
## Time difference of 33.84358 secs
```

``` r
round(ddsc_mod2_GGGI$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.049        0.221        1.007     1.056 0.046   7204.152 0.997   0.997
## 2        0.5         0.061        0.246        1.007     1.068 0.057   6164.576 0.997   0.997
```

``` r
round(ddsc_mod2_GGGI$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gggi.z.cm
## means_y1          -0.090 0.253    1.000           1.000    0.908           0.908    -0.070
## means_y1_scaled   -0.376 1.057    1.000           1.000    0.908           0.908    -0.070
## means_y2           0.067 0.225    0.908           0.908    1.000           1.000     0.190
## means_y2_scaled    0.282 0.939    0.908           0.908    1.000           1.000     0.190
## gggi.z.cm          0.000 1.000   -0.070          -0.070    0.190           0.190     1.000
## gggi.z.cm_scaled   0.000 1.000   -0.070          -0.070    0.190           0.190     1.000
## diff_score        -0.157 0.106    0.461           0.461    0.046           0.046    -0.569
## diff_score_scaled -0.657 0.444    0.461           0.461    0.046           0.046    -0.569
##                   gggi.z.cm_scaled diff_score diff_score_scaled
## means_y1                    -0.070      0.461             0.461
## means_y1_scaled             -0.070      0.461             0.461
## means_y2                     0.190      0.046             0.046
## means_y2_scaled              0.190      0.046             0.046
## gggi.z.cm                    1.000     -0.569            -0.569
## gggi.z.cm_scaled             1.000     -0.569            -0.569
## diff_score                  -0.569      1.000             1.000
## diff_score_scaled           -0.569      1.000             1.000
```

``` r
round(ddsc_mod2_GGGI$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                       0.553 0.139 32.258   3.993   0.000    0.271    0.835
## w_11                          0.042 0.039 33.263   1.085   0.286   -0.037    0.121
## w_21                         -0.017 0.044 33.116  -0.386   0.702   -0.106    0.072
## r_xy1                         0.166 0.153 33.263   1.085   0.286   -0.145    0.476
## r_xy2                        -0.075 0.194 33.116  -0.386   0.702   -0.471    0.320
## b_11                          0.176 0.162 33.263   1.085   0.286   -0.154    0.505
## b_21                         -0.071 0.183 33.116  -0.386   0.702   -0.443    0.302
## main_effect                   0.013 0.041 33.102   0.309   0.760   -0.070    0.095
## moderator_effect             -0.158 0.014 30.097 -11.055   0.000   -0.187   -0.128
## interaction                  -0.059 0.015 32.258  -3.993   0.000   -0.089   -0.029
## q_b11_b21                     0.248    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.242    NA     NA      NA      NA       NA       NA
## cross_over_point             -2.676    NA     NA      NA      NA       NA       NA
## interaction_vs_main           0.046 0.048 33.162   0.969   0.340   -0.051    0.144
## interaction_vs_main_bscale    0.194 0.200 33.162   0.969   0.340   -0.213    0.600
## interaction_vs_main_rscale    0.195 0.221 33.156   0.884   0.383   -0.254    0.645
## dadas                         0.034 0.087 33.116   0.386   0.351   -0.144    0.212
## dadas_bscale                  0.141 0.366 33.116   0.386   0.351   -0.603    0.886
## dadas_rscale                  0.150 0.389 33.116   0.386   0.351   -0.641    0.941
## abs_diff                      0.059 0.015 32.258   3.993   0.000    0.029    0.089
## abs_sum                       0.025 0.081 33.102   0.309   0.380   -0.140    0.190
## abs_diff_bscale               0.246 0.062 32.258   3.993   0.000    0.121    0.372
## abs_sum_bscale                0.105 0.340 33.102   0.309   0.380   -0.587    0.796
## abs_diff_rscale               0.241 0.071 32.458   3.370   0.001    0.095    0.386
## abs_sum_rscale                0.091 0.342 33.100   0.265   0.396   -0.606    0.787
```

``` r
round(ddsc_mod2_GGGI$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
##  0.006  0.279  2.571  1.000  0.109
```

``` r
d_GGGI<-ddsc_mod2_GGGI$ddsc_sem_fit$data

ddsc_sem_GGGI<-
  ddsc_sem(data=d_GGGI,x = "gggi.z.cm",y1="means_y2",
           y2="means_y1",q_sesoi = .10)
round(ddsc_sem_GGGI$results,3)
```

```
##                                    est    se      z pvalue ci.lower ci.upper
## r_xy1_y2                         0.569 0.143  3.973  0.000    0.288    0.849
## r_xy1                            0.190 0.171  1.110  0.267   -0.145    0.525
## r_xy2                           -0.070 0.174 -0.405  0.685   -0.411    0.270
## b_11                             0.178 0.161  1.110  0.267   -0.136    0.493
## b_21                            -0.074 0.184 -0.405  0.685   -0.434    0.285
## b_10                             0.282 0.158  1.782  0.075   -0.028    0.591
## b_20                            -0.376 0.181 -2.078  0.038   -0.730   -0.021
## res_cov_y1_y2                    0.887 0.225  3.935  0.000    0.445    1.329
## diff_b10_b20                     0.657 0.063 10.498  0.000    0.535    0.780
## diff_b11_b21                     0.253 0.064  3.973  0.000    0.128    0.377
## diff_rxy1_rxy2                   0.260 0.060  4.371  0.000    0.144    0.377
## q_b11_b21                        0.255 0.063  4.027  0.000    0.131    0.379
## q_rxy1_rxy2                      0.263 0.061  4.322  0.000    0.144    0.382
## cross_over_point                -2.602 0.700 -3.716  0.000   -3.974   -1.230
## sum_b11_b21                      0.104 0.339  0.306  0.759   -0.561    0.768
## main_effect                      0.052 0.169  0.306  0.759   -0.280    0.384
## interaction_vs_main_effect       0.201 0.202  0.995  0.320   -0.195    0.596
## diff_abs_b11_abs_b21             0.104 0.339  0.306  0.759   -0.561    0.768
## abs_diff_b11_b21                 0.253 0.064  3.973  0.000    0.128    0.377
## abs_sum_b11_b21                  0.104 0.339  0.306  0.380   -0.561    0.768
## dadas                            0.149 0.367  0.405  0.343   -0.571    0.869
## q_r_equivalence                  0.163 0.061  2.676  0.996       NA       NA
## q_b_equivalence                  0.155 0.063  2.446  0.993       NA       NA
## cross_over_point_equivalence     2.602 0.700  3.716  1.000       NA       NA
## cross_over_point_minimal_effect  2.602 0.700  3.716  0.000       NA       NA
```

``` r
round(ddsc_sem_GGGI$variance_test,3)
```

```
##              est    se      z pvalue ci.lower ci.upper
## cov_y1y2   0.874 0.226  3.861  0.000    0.430    1.318
## var_y1     0.855 0.211  4.062  0.000    0.443    1.268
## var_y2     1.084 0.267  4.062  0.000    0.561    1.607
## var_diff  -0.229 0.152 -1.510  0.131   -0.526    0.068
## var_ratio  0.789 0.115  6.845  0.000    0.563    1.015
## cor_y1y2   0.908 0.031 29.611  0.000    0.848    0.968
```

### Figure


``` r
# plot the results

# start with obtaining predicted values for means and differences

# refit reduced and full models with GGGI in original scale


mod2_GGGI_unstd<-lmer(uni.z~gndr.c+gggi.cm+gndr.c:gggi.cm+
                       (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
                     control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

mod2_GGGI_unstd_red<-lmer(uni.z~gndr.c+
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


p1.uni.flags<-
  ggplot(p,aes(y=yvar,x=gggi.cm,color=gndr.c))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("Global Gender Gap Index (Average 2002-2022)")+
  #ylim(c(min.y.pred,max.y.pred))+
  #xlim(c(0.60,1.00))+
  ylab("Value universalism (Average 2002-2022)")+
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

p2.uni.flags<-ggplot(p2,aes(y=yvar,x=gggi.cm))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("Global Gender Gap Index (Average 2002-2022)")+
  #xlim(c(0.60,1.00))+
  #ylim(c(min.y.narrow,max.y.narrow))+
  ylab("Difference in value universalism")+
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
  ggarrange(p1.uni.flags,p2.uni.flags,align = "v",
            ncol=1,nrow=2,heights=c(2,1))
```

```
## Warning in geom_text(inherit.aes = F, aes(x = 0.61, y = -0.45, label = coef_q, : All aesthetics have length 1, but the data has 502 rows.
## ℹ Please consider using `annotate()` or provide this layer with data containing a single row.
```

``` r
pflag_comb
```

![](Analysis_universalism_files/figure-html/unnamed-chunk-21-1.png)<!-- -->

``` r
png(filename = 
      "../results/uni/GGGI_flags_new.png",
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
mod2_GDI<-lmer(uni.z~gndr.c+gdi.z.cm+gndr.c:gdi.z.cm+
                 (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2_GDI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: uni.z ~ gndr.c + gdi.z.cm + gndr.c:gdi.z.cm + (gndr.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1310980.8 1311068.8 -655482.4 1310964.8    441160 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -8.3977 -0.5510  0.1027  0.6578  4.8959 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr
##  cntry    (Intercept) 0.037090 0.19259      
##           gndr.c      0.007134 0.08446  0.06
##  Residual             1.006999 1.00349      
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##                 Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)     -0.01652    0.03357 33.02597  -0.492 0.625963    
## gndr.c          -0.16039    0.01512 32.73129 -10.609 3.99e-12 ***
## gdi.z.cm        -0.12640    0.03412 33.13616  -3.704 0.000769 ***
## gndr.c:gdi.z.cm -0.04078    0.01560 34.76157  -2.615 0.013104 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c gd.z.c
## gndr.c       0.062              
## gdi.z.cm     0.000  0.000       
## gndr.c:gd..  0.000 -0.005  0.062
```

``` r
getFE(mod2_GDI,round=3)
```

```
##                   Est.    SE     df       t     p     LL     UL
## (Intercept)     -0.017 0.034 33.026  -0.492 0.626 -0.085  0.052
## gndr.c          -0.160 0.015 32.731 -10.609 0.000 -0.191 -0.130
## gdi.z.cm        -0.126 0.034 33.136  -3.704 0.001 -0.196 -0.057
## gndr.c:gdi.z.cm -0.041 0.016 34.762  -2.615 0.013 -0.072 -0.009
```

``` r
getVC(mod2_GDI)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.19 0.04
## 2    cntry      gndr.c   <NA>  0.08 0.01
## 3    cntry (Intercept) gndr.c  0.06 0.00
## 4 Residual        <NA>   <NA>  1.00 1.01
```

``` r
r2mlm(mod2_GDI,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.015724723
## slope variation 0.001668331
## mean variation  0.034842815
## sigma2          0.947764131
## 
## $R2s
##           total
## f   0.015724723
## v   0.001668331
## m   0.034842815
## fv  0.017393054
## fvm 0.052235869
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
## Time difference of 33.44373 secs
```

``` r
round(ddsc_mod2_GDI$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.049        0.221        1.007     1.056 0.046   7204.152 0.997   0.997
## 2        0.5         0.061        0.246        1.007     1.068 0.057   6164.576 0.997   0.997
```

``` r
round(ddsc_mod2_GDI$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gdi.z.cm
## means_y1          -0.097 0.251    1.000           1.000    0.921           0.921   -0.587
## means_y1_scaled   -0.407 1.056    1.000           1.000    0.921           0.921   -0.587
## means_y2           0.064 0.224    0.921           0.921    1.000           1.000   -0.469
## means_y2_scaled    0.268 0.941    0.921           0.921    1.000           1.000   -0.469
## gdi.z.cm           0.000 1.000   -0.587          -0.587   -0.469          -0.469    1.000
## gdi.z.cm_scaled    0.000 1.000   -0.587          -0.587   -0.469          -0.469    1.000
## diff_score        -0.161 0.098    0.459           0.459    0.076           0.076   -0.433
## diff_score_scaled -0.676 0.413    0.459           0.459    0.076           0.076   -0.433
##                   gdi.z.cm_scaled diff_score diff_score_scaled
## means_y1                   -0.587      0.459             0.459
## means_y1_scaled            -0.587      0.459             0.459
## means_y2                   -0.469      0.076             0.076
## means_y2_scaled            -0.469      0.076             0.076
## gdi.z.cm                    1.000     -0.433            -0.433
## gdi.z.cm_scaled             1.000     -0.433            -0.433
## diff_score                 -0.433      1.000             1.000
## diff_score_scaled          -0.433      1.000             1.000
```

``` r
round(ddsc_mod2_GDI$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                       0.415 0.159 34.762   2.615   0.013    0.093    0.737
## w_11                         -0.106 0.035 33.254  -3.070   0.004   -0.176   -0.036
## w_21                         -0.147 0.035 33.251  -4.139   0.000   -0.219   -0.075
## r_xy1                        -0.422 0.137 33.254  -3.070   0.004   -0.701   -0.142
## r_xy2                        -0.655 0.158 33.251  -4.139   0.000   -0.977   -0.333
## b_11                         -0.446 0.145 33.254  -3.070   0.004   -0.741   -0.150
## b_21                         -0.617 0.149 33.251  -4.139   0.000   -0.921   -0.314
## main_effect                  -0.126 0.034 33.136  -3.704   0.001   -0.196   -0.057
## moderator_effect             -0.160 0.015 32.731 -10.609   0.000   -0.191   -0.130
## interaction                  -0.041 0.016 34.762  -2.615   0.013   -0.072   -0.009
## q_b11_b21                     0.241    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.335    NA     NA      NA      NA       NA       NA
## cross_over_point             -3.933    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.086 0.037 33.556  -2.337   0.026   -0.160   -0.011
## interaction_vs_main_bscale   -0.360 0.154 33.556  -2.337   0.026   -0.673   -0.047
## interaction_vs_main_rscale   -0.305 0.139 33.624  -2.193   0.035   -0.588   -0.022
## dadas                        -0.212 0.069 33.254  -3.070   0.998   -0.352   -0.072
## dadas_bscale                 -0.892 0.290 33.254  -3.070   0.998   -1.483   -0.301
## dadas_rscale                 -0.843 0.275 33.254  -3.070   0.998   -1.402   -0.285
## abs_diff                      0.041 0.016 34.762   2.615   0.007    0.009    0.072
## abs_sum                       0.253 0.068 33.136   3.704   0.000    0.114    0.392
## abs_diff_bscale               0.172 0.066 34.762   2.615   0.007    0.038    0.305
## abs_sum_bscale                1.063 0.287 33.136   3.704   0.000    0.479    1.647
## abs_diff_rscale               0.233 0.069 34.625   3.391   0.001    0.094    0.373
## abs_sum_rscale                1.077 0.288 33.137   3.735   0.000    0.490    1.663
```

``` r
round(ddsc_mod2_GDI$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
##  0.006  0.279  2.571  1.000  0.109
```

``` r
d_GDI<-ddsc_mod2_GDI$ddsc_sem_fit$data

ddsc_sem_GDI<-
  ddsc_sem(data=d_GDI,x = "gdi.z.cm",y1="means_y2",
           y2="means_y1",q_sesoi = .10)
round(ddsc_sem_GDI$results,3)
```

```
##                                    est    se      z pvalue ci.lower ci.upper
## r_xy1_y2                         0.433 0.157  2.760  0.006    0.125    0.741
## r_xy1                           -0.469 0.154 -3.051  0.002   -0.770   -0.168
## r_xy2                           -0.587 0.141 -4.168  0.000   -0.863   -0.311
## b_11                            -0.441 0.145 -3.051  0.002   -0.725   -0.158
## b_21                            -0.620 0.149 -4.168  0.000   -0.912   -0.328
## b_10                             0.268 0.142  1.883  0.060   -0.011    0.547
## b_20                            -0.407 0.146 -2.781  0.005   -0.695   -0.120
## res_cov_y1_y2                    0.622 0.162  3.850  0.000    0.305    0.938
## diff_b10_b20                     0.676 0.064 10.599  0.000    0.551    0.801
## diff_b11_b21                     0.179 0.065  2.760  0.006    0.052    0.306
## diff_rxy1_rxy2                   0.118 0.066  1.789  0.074   -0.011    0.248
## q_b11_b21                        0.251 0.111  2.267  0.023    0.034    0.468
## q_rxy1_rxy2                      0.165 0.092  1.781  0.075   -0.017    0.346
## cross_over_point                -3.782 1.416 -2.671  0.008   -6.557   -1.006
## sum_b11_b21                     -1.061 0.286 -3.708  0.000   -1.622   -0.500
## main_effect                     -0.531 0.143 -3.708  0.000   -0.811   -0.250
## interaction_vs_main_effect      -0.352 0.153 -2.298  0.022   -0.652   -0.052
## diff_abs_b11_abs_b21            -0.179 0.065 -2.760  0.006   -0.306   -0.052
## abs_diff_b11_b21                 0.179 0.065  2.760  0.003    0.052    0.306
## abs_sum_b11_b21                  1.061 0.286  3.708  0.000    0.500    1.622
## dadas                           -0.883 0.289 -3.051  0.999   -1.450   -0.316
## q_r_equivalence                  0.065 0.092  0.699  0.758       NA       NA
## q_b_equivalence                  0.151 0.111  1.364  0.914       NA       NA
## cross_over_point_equivalence     3.782 1.416  2.671  0.996       NA       NA
## cross_over_point_minimal_effect  3.782 1.416  2.671  0.004       NA       NA
```

``` r
round(ddsc_sem_GDI$variance_test,3)
```

```
##              est    se      z pvalue ci.lower ci.upper
## cov_y1y2   0.887 0.228  3.892  0.000    0.440    1.334
## var_y1     0.859 0.211  4.062  0.000    0.444    1.273
## var_y2     1.081 0.266  4.062  0.000    0.559    1.602
## var_diff  -0.222 0.142 -1.569  0.117   -0.500    0.055
## var_ratio  0.794 0.108  7.371  0.000    0.583    1.005
## cor_y1y2   0.921 0.026 34.845  0.000    0.869    0.973
```

### Figure


``` r
# plot the results

# start with obtaining predicted values for means and differences


mod2_GDI_unstd<-lmer(uni.z~gndr.c+gdi.cm+gndr.c:gdi.cm+
                       (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
                     control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

mod2_GDI_unstd_red<-lmer(uni.z~gndr.c+
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


p1.uni.flags<-
  ggplot(p,aes(y=yvar,x=gdi.cm,color=gndr.c))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("Gender Development Index (Average 2002-2022)")+
  #ylim(c(min.y.pred,max.y.pred))+
  #xlim(c(0.60,1.00))+
  ylab("Value universalism (Average 2002-2022)")+
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

#p1.uni.flags


p2.uni.flags<-ggplot(p2,aes(y=yvar,x=gdi.cm))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("Gender Development Index (Average 2002-2022)")+
  #xlim(c(0.60,1.00))+
  #ylim(c(min.y.narrow,max.y.narrow))+
  ylab("Difference in value universalism")+
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

#p2.uni.flags


pflag_comb<-
  ggarrange(p1.uni.flags,p2.uni.flags,align = "v",
            ncol=1,nrow=2,heights=c(2,1))
```

```
## Warning in geom_text(inherit.aes = F, aes(x = 0.9, y = -0.45, label = coef_q, : All aesthetics have length 1, but the data has 282 rows.
## ℹ Please consider using `annotate()` or provide this layer with data containing a single row.
```

``` r
pflag_comb
```

![](Analysis_universalism_files/figure-html/unnamed-chunk-24-1.png)<!-- -->

``` r
png(filename = 
      "../results/uni/GDI_flags.png",
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
mod2_log_GDP<-lmer(uni.z~gndr.c+log_gdp.z.cm+
                     gndr.c:log_gdp.z.cm+
                 (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2_log_GDP)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: uni.z ~ gndr.c + log_gdp.z.cm + gndr.c:log_gdp.z.cm + (gndr.c |      cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1310988.7 1311076.7 -655486.3 1310972.7    441160 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -8.3933 -0.5508  0.1022  0.6579  4.8954 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr
##  cntry    (Intercept) 0.047632 0.21825      
##           gndr.c      0.008088 0.08993  0.38
##  Residual             1.007001 1.00349      
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##                     Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)         -0.01505    0.03804 33.00247  -0.396    0.695    
## gndr.c              -0.16045    0.01605 31.54722  -9.999 2.65e-11 ***
## log_gdp.z.cm         0.07046    0.03819 33.04365   1.845    0.074 .  
## gndr.c:log_gdp.z.cm -0.02265    0.01620 32.29566  -1.397    0.172    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c lg_g..
## gndr.c      0.372               
## lg_gdp.z.cm 0.021  0.008        
## gndr.c:l_.. 0.008  0.006  0.370
```

``` r
getFE(mod2_log_GDP,round=3)
```

```
##                       Est.    SE     df      t     p     LL     UL
## (Intercept)         -0.015 0.038 33.002 -0.396 0.695 -0.092  0.062
## gndr.c              -0.160 0.016 31.547 -9.999 0.000 -0.193 -0.128
## log_gdp.z.cm         0.070 0.038 33.044  1.845 0.074 -0.007  0.148
## gndr.c:log_gdp.z.cm -0.023 0.016 32.296 -1.397 0.172 -0.056  0.010
```

``` r
getVC(mod2_log_GDP)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.22 0.05
## 2    cntry      gndr.c   <NA>  0.09 0.01
## 3    cntry (Intercept) gndr.c  0.38 0.01
## 4 Residual        <NA>   <NA>  1.00 1.01
```

``` r
r2mlm(mod2_log_GDP,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.009950578
## slope variation 0.001884100
## mean variation  0.044119977
## sigma2          0.944045344
## 
## $R2s
##           total
## f   0.009950578
## v   0.001884100
## m   0.044119977
## fv  0.011834679
## fvm 0.055954656
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
## Time difference of 33.97931 secs
```

``` r
round(ddsc_mod2_log_GDP$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.049        0.221        1.007     1.056 0.046   7204.152 0.997   0.997
## 2        0.5         0.061        0.246        1.007     1.068 0.057   6164.576 0.997   0.997
```

``` r
round(ddsc_mod2_log_GDP$descriptives,3)
```

```
##                          M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled log_gdp.z.cm
## means_y1            -0.097 0.251    1.000           1.000    0.921           0.921        0.240
## means_y1_scaled     -0.407 1.056    1.000           1.000    0.921           0.921        0.240
## means_y2             0.064 0.224    0.921           0.921    1.000           1.000        0.367
## means_y2_scaled      0.268 0.941    0.921           0.921    1.000           1.000        0.367
## log_gdp.z.cm        -0.022 1.012    0.240           0.240    0.367           0.367        1.000
## log_gdp.z.cm_scaled  0.000 1.000    0.240           0.240    0.367           0.367        1.000
## diff_score          -0.161 0.098    0.459           0.459    0.076           0.076       -0.224
## diff_score_scaled   -0.676 0.413    0.459           0.459    0.076           0.076       -0.224
##                     log_gdp.z.cm_scaled diff_score diff_score_scaled
## means_y1                          0.240      0.459             0.459
## means_y1_scaled                   0.240      0.459             0.459
## means_y2                          0.367      0.076             0.076
## means_y2_scaled                   0.367      0.076             0.076
## log_gdp.z.cm                      1.000     -0.224            -0.224
## log_gdp.z.cm_scaled               1.000     -0.224            -0.224
## diff_score                       -0.224      1.000             1.000
## diff_score_scaled                -0.224      1.000             1.000
```

``` r
round(ddsc_mod2_log_GDP$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                       0.230 0.165 32.296   1.397   0.172   -0.105    0.566
## w_11                          0.082 0.036 33.126   2.272   0.030    0.009    0.155
## w_21                          0.059 0.042 33.039   1.413   0.167   -0.026    0.144
## r_xy1                         0.325 0.143 33.126   2.272   0.030    0.034    0.616
## r_xy2                         0.264 0.187 33.039   1.413   0.167   -0.116    0.644
## b_11                          0.344 0.151 33.126   2.272   0.030    0.036    0.652
## b_21                          0.249 0.176 33.039   1.413   0.167   -0.109    0.607
## main_effect                   0.070 0.038 33.044   1.845   0.074   -0.007    0.148
## moderator_effect             -0.160 0.016 31.547  -9.999   0.000   -0.193   -0.128
## interaction                  -0.023 0.016 32.296  -1.397   0.172   -0.056    0.010
## q_b11_b21                     0.105    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.067    NA     NA      NA      NA       NA       NA
## cross_over_point             -7.085    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.048 0.047 33.025  -1.025   0.313   -0.143    0.047
## interaction_vs_main_bscale   -0.201 0.196 33.025  -1.025   0.313   -0.601    0.198
## interaction_vs_main_rscale   -0.233 0.216 33.029  -1.081   0.288   -0.672    0.206
## dadas                        -0.118 0.084 33.039  -1.413   0.916   -0.289    0.052
## dadas_bscale                 -0.497 0.352 33.039  -1.413   0.916   -1.214    0.219
## dadas_rscale                 -0.528 0.374 33.039  -1.413   0.916   -1.288    0.232
## abs_diff                      0.023 0.016 32.296   1.397   0.086   -0.010    0.056
## abs_sum                       0.141 0.076 33.044   1.845   0.037   -0.014    0.296
## abs_diff_bscale               0.095 0.068 32.296   1.397   0.086   -0.044    0.234
## abs_sum_bscale                0.593 0.321 33.044   1.845   0.037   -0.061    1.246
## abs_diff_rscale               0.061 0.077 32.418   0.795   0.216   -0.096    0.219
## abs_sum_rscale                0.589 0.324 33.042   1.820   0.039   -0.069    1.248
```

``` r
round(ddsc_mod2_log_GDP$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
##  0.006  0.279  2.571  1.000  0.109
```

``` r
d_log_GDP<-ddsc_mod2_log_GDP$ddsc_sem_fit$data

ddsc_sem_log_GDP<-
  ddsc_sem(data=d_log_GDP,x = "log_gdp.z.cm",y1="means_y2",
           y2="means_y1",q_sesoi = .10)
round(ddsc_sem_log_GDP$results,3)
```

```
##                                    est    se      z pvalue ci.lower ci.upper
## r_xy1_y2                         0.224 0.170  1.320  0.187   -0.109    0.556
## r_xy1                            0.367 0.162  2.268  0.023    0.050    0.685
## r_xy2                            0.240 0.169  1.419  0.156   -0.091    0.571
## b_11                             0.346 0.152  2.268  0.023    0.047    0.644
## b_21                             0.253 0.178  1.419  0.156   -0.097    0.603
## b_10                             0.268 0.150  1.788  0.074   -0.026    0.562
## b_20                            -0.407 0.176 -2.319  0.020   -0.752   -0.063
## res_cov_y1_y2                    0.802 0.206  3.895  0.000    0.399    1.206
## diff_b10_b20                     0.676 0.069  9.802  0.000    0.541    0.811
## diff_b11_b21                     0.092 0.070  1.320  0.187   -0.045    0.230
## diff_rxy1_rxy2                   0.127 0.066  1.944  0.052   -0.001    0.256
## q_b11_b21                        0.102 0.074  1.378  0.168   -0.043    0.246
## q_rxy1_rxy2                      0.141 0.073  1.937  0.053   -0.002    0.283
## cross_over_point                -7.314 5.592 -1.308  0.191  -18.275    3.647
## sum_b11_b21                      0.599 0.324  1.846  0.065   -0.037    1.234
## main_effect                      0.299 0.162  1.846  0.065   -0.019    0.617
## interaction_vs_main_effect      -0.207 0.200 -1.037  0.300   -0.598    0.184
## diff_abs_b11_abs_b21             0.092 0.070  1.320  0.187   -0.045    0.230
## abs_diff_b11_b21                 0.092 0.070  1.320  0.093   -0.045    0.230
## abs_sum_b11_b21                  0.599 0.324  1.846  0.032   -0.037    1.234
## dadas                           -0.506 0.357 -1.419  0.922   -1.206    0.193
## q_r_equivalence                  0.041 0.073  0.560  0.712       NA       NA
## q_b_equivalence                  0.002 0.074  0.021  0.509       NA       NA
## cross_over_point_equivalence     7.314 5.592  1.308  0.905       NA       NA
## cross_over_point_minimal_effect  7.314 5.592  1.308  0.095       NA       NA
```

``` r
round(ddsc_sem_log_GDP$variance_test,3)
```

```
##              est    se      z pvalue ci.lower ci.upper
## cov_y1y2   0.887 0.228  3.892  0.000    0.440    1.334
## var_y1     0.859 0.211  4.062  0.000    0.444    1.273
## var_y2     1.081 0.266  4.062  0.000    0.559    1.602
## var_diff  -0.222 0.142 -1.569  0.117   -0.500    0.055
## var_ratio  0.794 0.108  7.371  0.000    0.583    1.005
## cor_y1y2   0.921 0.026 34.845  0.000    0.869    0.973
```

### Figure


``` r
# plot the results

# start with obtaining predicted values for means and differences


mod2_log_GDP_unstd<-lmer(uni.z~gndr.c+log_gdp.cm+
                           gndr.c:log_gdp.cm+
                       (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
                     control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

mod2_log_GDP_unstd_red<-lmer(uni.z~gndr.c+
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


p1.uni.flags<-
  ggplot(p,aes(y=yvar,x=log_gdp.cm,color=gndr.c))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("log(GDP per capita) (Average 2002-2022)")+
  #ylim(c(min.y.pred,max.y.pred))+
  #xlim(c(0.60,1.00))+
  ylab("Value universalism (Average 2002-2022)")+
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

p2.uni.flags<-ggplot(p2,aes(y=yvar,x=log_gdp.cm))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("log(GDP per capita) (Average 2002-2022)")+
  #xlim(c(0.60,1.00))+
  #ylim(c(min.y.narrow,max.y.narrow))+
  ylab("Difference in value universalism")+
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
  ggarrange(p1.uni.flags,p2.uni.flags,align = "v",
            ncol=1,nrow=2,heights=c(2,1))
```

```
## Warning in geom_text(inherit.aes = F, aes(x = 9.75, y = -0.45, label = coef_q, : All aesthetics have length 1, but the data has 3302 rows.
## ℹ Please consider using `annotate()` or provide this layer with data containing a single row.
```

``` r
pflag_comb
```

![](Analysis_universalism_files/figure-html/unnamed-chunk-27-1.png)<!-- -->

``` r
png(filename = 
      "../results/uni/log_GDP_flags.png",
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
mod3<-lmer(uni.z~gndr.c+essround.c+(gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod3)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: uni.z ~ gndr.c + essround.c + (gndr.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1310791.3 1310868.3 -655388.7 1310777.3    441161 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -8.4011 -0.5506  0.1046  0.6582  4.9403 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr
##  cntry    (Intercept) 0.053792 0.2319       
##           gndr.c      0.008668 0.0931   0.29
##  Residual             1.006535 1.0033       
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##               Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept) -1.756e-02  4.041e-02  3.301e+01  -0.434    0.667    
## gndr.c      -1.608e-01  1.659e-02  3.183e+01  -9.694  5.1e-11 ***
## essround.c   8.214e-03  5.766e-04  4.408e+05  14.245  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##            (Intr) gndr.c
## gndr.c      0.280       
## essround.c -0.002 -0.001
```

``` r
getFE(mod3,round=3)
```

```
##               Est.    SE         df      t     p     LL     UL
## (Intercept) -0.018 0.040     33.005 -0.434 0.667 -0.100  0.065
## gndr.c      -0.161 0.017     31.827 -9.694 0.000 -0.195 -0.127
## essround.c   0.008 0.001 440762.784 14.245 0.000  0.007  0.009
```

``` r
getVC(mod3)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.23 0.05
## 2    cntry      gndr.c   <NA>  0.09 0.01
## 3    cntry (Intercept) gndr.c  0.29 0.01
## 4 Residual        <NA>   <NA>  1.00 1.01
```

``` r
r2mlm(mod3,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.006469549
## slope variation 0.002014932
## mean variation  0.049885275
## sigma2          0.941630244
## 
## $R2s
##           total
## f   0.006469549
## v   0.002014932
## m   0.049885275
## fv  0.008484481
## fvm 0.058369756
```

``` r
anova(mod2,mod3)
```

```
## Data: diff_dat
## Models:
## mod2: uni.z ~ gndr.c + (gndr.c | cntry)
## mod3: uni.z ~ gndr.c + essround.c + (gndr.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod2    6 1310992 1311058 -655490   1310980                         
## mod3    7 1310791 1310868 -655389   1310777 202.85  1  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## mod4: random effect of time (Ess round)


``` r
mod4<-lmer(uni.z~gndr.c+essround.c+(gndr.c+essround.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod4)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: uni.z ~ gndr.c + essround.c + (gndr.c + essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1307056.2 1307166.2 -653518.1 1307036.2    441158 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -8.4372 -0.5524  0.1071  0.6540  4.9180 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr       
##  cntry    (Intercept) 0.050262 0.22419             
##           gndr.c      0.008779 0.09369   0.25      
##           essround.c  0.001482 0.03850   0.23 -0.40
##  Residual             0.997720 0.99886             
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##              Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept) -0.011740   0.039178 33.233112  -0.300    0.766    
## gndr.c      -0.161422   0.016682 31.893759  -9.676 5.21e-11 ***
## essround.c   0.005447   0.006786 32.911780   0.803    0.428    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##            (Intr) gndr.c
## gndr.c      0.249       
## essround.c  0.220 -0.392
```

``` r
getFE(mod4,round=3)
```

```
##               Est.    SE     df      t     p     LL     UL
## (Intercept) -0.012 0.039 33.233 -0.300 0.766 -0.091  0.068
## gndr.c      -0.161 0.017 31.894 -9.676 0.000 -0.195 -0.127
## essround.c   0.005 0.007 32.912  0.803 0.428 -0.008  0.019
```

``` r
getVC(mod4)
```

```
##        grp        var1       var2 sdcor vcov
## 1    cntry (Intercept)       <NA>  0.22 0.05
## 2    cntry      gndr.c       <NA>  0.09 0.01
## 3    cntry  essround.c       <NA>  0.04 0.00
## 4    cntry (Intercept)     gndr.c  0.25 0.01
## 5    cntry (Intercept) essround.c  0.23 0.00
## 6    cntry      gndr.c essround.c -0.40 0.00
## 7 Residual        <NA>       <NA>  1.00 1.00
```

``` r
r2mlm(mod4,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.006263849
## slope variation 0.012512704
## mean variation  0.046619877
## sigma2          0.934603570
## 
## $R2s
##           total
## f   0.006263849
## v   0.012512704
## m   0.046619877
## fv  0.018776553
## fvm 0.065396430
```

``` r
anova(mod2,mod3,mod4)
```

```
## Data: diff_dat
## Models:
## mod2: uni.z ~ gndr.c + (gndr.c | cntry)
## mod3: uni.z ~ gndr.c + essround.c + (gndr.c | cntry)
## mod4: uni.z ~ gndr.c + essround.c + (gndr.c + essround.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L)   Chisq Df Pr(>Chisq)    
## mod2    6 1310992 1311058 -655490   1310980                          
## mod3    7 1310791 1310868 -655389   1310777  202.85  1  < 2.2e-16 ***
## mod4   10 1307056 1307166 -653518   1307036 3741.12  3  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## mod5: fixed interaction between time and gender


``` r
mod5<-lmer(uni.z~gndr.c+essround.c+
             gndr.c:essround.c+(gndr.c+essround.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod5)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: uni.z ~ gndr.c + essround.c + gndr.c:essround.c + (gndr.c + essround.c |      cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1307056.1 1307177.1 -653517.1 1307034.1    441157 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -8.4370 -0.5523  0.1075  0.6540  4.9210 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr       
##  cntry    (Intercept) 0.050263 0.22420             
##           gndr.c      0.008691 0.09323   0.25      
##           essround.c  0.001484 0.03852   0.23 -0.41
##  Residual             0.997716 0.99886             
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##                     Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)       -1.164e-02  3.918e-02  3.323e+01  -0.297    0.768    
## gndr.c            -1.612e-01  1.660e-02  3.189e+01  -9.712 4.78e-11 ***
## essround.c         5.394e-03  6.790e-03  3.290e+01   0.794    0.433    
## gndr.c:essround.c -1.640e-03  1.144e-03  2.427e+05  -1.433    0.152    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c essrn.
## gndr.c       0.247              
## essround.c   0.220 -0.396       
## gndr.c:ssr. -0.002 -0.007  0.005
```

``` r
getFE(mod5,round=3)
```

```
##                     Est.    SE         df      t     p     LL     UL
## (Intercept)       -0.012 0.039     33.231 -0.297 0.768 -0.091  0.068
## gndr.c            -0.161 0.017     31.889 -9.712 0.000 -0.195 -0.127
## essround.c         0.005 0.007     32.900  0.794 0.433 -0.008  0.019
## gndr.c:essround.c -0.002 0.001 242717.133 -1.433 0.152 -0.004  0.001
```

``` r
getVC(mod5)
```

```
##        grp        var1       var2 sdcor vcov
## 1    cntry (Intercept)       <NA>  0.22 0.05
## 2    cntry      gndr.c       <NA>  0.09 0.01
## 3    cntry  essround.c       <NA>  0.04 0.00
## 4    cntry (Intercept)     gndr.c  0.25 0.01
## 5    cntry (Intercept) essround.c  0.23 0.00
## 6    cntry      gndr.c essround.c -0.41 0.00
## 7 Residual        <NA>       <NA>  1.00 1.00
```

``` r
r2mlm(mod5,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.006253493
## slope variation 0.012504024
## mean variation  0.046627341
## sigma2          0.934615142
## 
## $R2s
##           total
## f   0.006253493
## v   0.012504024
## m   0.046627341
## fv  0.018757517
## fvm 0.065384858
```

``` r
anova(mod4,mod5)
```

```
## Data: diff_dat
## Models:
## mod4: uni.z ~ gndr.c + essround.c + (gndr.c + essround.c | cntry)
## mod5: uni.z ~ gndr.c + essround.c + gndr.c:essround.c + (gndr.c + essround.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)
## mod4   10 1307056 1307166 -653518   1307036                     
## mod5   11 1307056 1307177 -653517   1307034 2.0508  1     0.1521
```

## mod6: random interaction between time and gender


``` r
mod6<-lmer(uni.z~gndr.c+essround.c+
             gndr.c:essround.c+(gndr.c+essround.c+gndr.c:essround.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
```

```
## Warning in checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv, : Model failed to converge
## with max|grad| = 0.00224081 (tol = 0.002, component 1)
```

``` r
summary(mod6)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: uni.z ~ gndr.c + essround.c + gndr.c:essround.c + (gndr.c + essround.c +  
##     gndr.c:essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##   1307030   1307195   -653500   1307000    441153 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -8.4344 -0.5508  0.1053  0.6545  4.9358 
## 
## Random effects:
##  Groups   Name              Variance  Std.Dev. Corr             
##  cntry    (Intercept)       5.026e-02 0.224198                  
##           gndr.c            8.587e-03 0.092669  0.26            
##           essround.c        1.486e-03 0.038555  0.23 -0.42      
##           gndr.c:essround.c 8.402e-05 0.009166  0.05 -0.49  0.24
##  Residual                   9.976e-01 0.998790                  
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##                    Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)       -0.011694   0.039179 33.226863  -0.298    0.767    
## gndr.c            -0.161754   0.016590 31.977427  -9.750 4.23e-11 ***
## essround.c         0.005386   0.006796 32.920521   0.792    0.434    
## gndr.c:essround.c -0.002387   0.002070 27.638373  -1.153    0.259    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c essrn.
## gndr.c       0.249              
## essround.c   0.220 -0.408       
## gndr.c:ssr.  0.039 -0.392  0.189
## optimizer (bobyqa) convergence code: 0 (OK)
## Model failed to converge with max|grad| = 0.00224081 (tol = 0.002, component 1)
```

``` r
getFE(mod6,round=3)
```

```
##                     Est.    SE     df      t     p     LL     UL
## (Intercept)       -0.012 0.039 33.227 -0.298 0.767 -0.091  0.068
## gndr.c            -0.162 0.017 31.977 -9.750 0.000 -0.196 -0.128
## essround.c         0.005 0.007 32.921  0.792 0.434 -0.008  0.019
## gndr.c:essround.c -0.002 0.002 27.638 -1.153 0.259 -0.007  0.002
```

``` r
getVC(mod6)
```

```
##         grp              var1              var2 sdcor vcov
## 1     cntry       (Intercept)              <NA>  0.22 0.05
## 2     cntry            gndr.c              <NA>  0.09 0.01
## 3     cntry        essround.c              <NA>  0.04 0.00
## 4     cntry gndr.c:essround.c              <NA>  0.01 0.00
## 5     cntry       (Intercept)            gndr.c  0.26 0.01
## 6     cntry       (Intercept)        essround.c  0.23 0.00
## 7     cntry       (Intercept) gndr.c:essround.c  0.05 0.00
## 8     cntry            gndr.c        essround.c -0.42 0.00
## 9     cntry            gndr.c gndr.c:essround.c -0.49 0.00
## 10    cntry        essround.c gndr.c:essround.c  0.24 0.00
## 11 Residual              <NA>              <NA>  1.00 1.00
```

``` r
r2mlm(mod6,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.006296762
## slope variation 0.012608470
## mean variation  0.046626219
## sigma2          0.934468549
## 
## $R2s
##           total
## f   0.006296762
## v   0.012608470
## m   0.046626219
## fv  0.018905232
## fvm 0.065531451
```

``` r
anova(mod4,mod5,mod6)
```

```
## Data: diff_dat
## Models:
## mod4: uni.z ~ gndr.c + essround.c + (gndr.c + essround.c | cntry)
## mod5: uni.z ~ gndr.c + essround.c + gndr.c:essround.c + (gndr.c + essround.c | cntry)
## mod6: uni.z ~ gndr.c + essround.c + gndr.c:essround.c + (gndr.c + essround.c + gndr.c:essround.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L)   Chisq Df Pr(>Chisq)    
## mod4   10 1307056 1307166 -653518   1307036                          
## mod5   11 1307056 1307177 -653517   1307034  2.0508  1     0.1521    
## mod6   15 1307030 1307195 -653500   1307000 34.1520  4  6.936e-07 ***
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
```

```
## Note: D.f. calculations have been disabled because the number of observations exceeds 4e+05.
## To enable adjustments, add the argument 'lmerTest.limit = 441168' (or larger)
## [or, globally, 'set emm_options(lmerTest.limit = 441168)' or larger];
## but be warned that this may result in large computation time and memory use.
```

``` r
change_mod6
```

```
## gndr.c = -0.5:
##  essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##         4.5  0.0988 0.0551 Inf   -0.0092    0.2068   1.793  0.0730
##        -4.5  0.0396 0.0407 Inf   -0.0401    0.1193   0.973  0.3305
## 
## gndr.c =  0.5:
##  essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##         4.5 -0.0737 0.0555 Inf   -0.1825    0.0351  -1.328  0.1843
##        -4.5 -0.1114 0.0497 Inf   -0.2089   -0.0140  -2.242  0.0250
## 
## Degrees-of-freedom method: asymptotic 
## Confidence level used: 0.95
```

``` r
pairs(change_mod6,adjust="none",infer=c(T,T))
```

```
## gndr.c = -0.5:
##  contrast                         estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   0.0592 0.0601 Inf   -0.0586     0.177   0.985  0.3246
## 
## gndr.c =  0.5:
##  contrast                         estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   0.0377 0.0636 Inf   -0.0869     0.162   0.593  0.5529
## 
## Degrees-of-freedom method: asymptotic 
## Confidence level used: 0.95
```

``` r
# change in gender differences over time

change_in_diff_mod6<-emmeans(mod6,specs=c("gndr.c","essround.c"),
                             at=list(gndr.c=c(-0.5,0.5),
                                     essround.c=rev(range(diff_dat$essround.c))),
                             disable.pbkrtest=T,
                             lmerTest.limit = 400000,infer=c(T,T))
```

```
## Note: D.f. calculations have been disabled because the number of observations exceeds 4e+05.
## To enable adjustments, add the argument 'lmerTest.limit = 441168' (or larger)
## [or, globally, 'set emm_options(lmerTest.limit = 441168)' or larger];
## but be warned that this may result in large computation time and memory use.
```

``` r
change_in_diff_mod6
```

```
##  gndr.c essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##    -0.5        4.5  0.0988 0.0551 Inf   -0.0092    0.2068   1.793  0.0730
##     0.5        4.5 -0.0737 0.0555 Inf   -0.1825    0.0351  -1.328  0.1843
##    -0.5       -4.5  0.0396 0.0407 Inf   -0.0401    0.1193   0.973  0.3305
##     0.5       -4.5 -0.1114 0.0497 Inf   -0.2089   -0.0140  -2.242  0.0250
## 
## Degrees-of-freedom method: asymptotic 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod6,adjust="none",infer=c(T,T))
```

```
##  contrast                                                 estimate     SE  df asymp.LCL asymp.UCL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5       0.1725 0.0155 Inf    0.1421   0.20291
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   0.0592 0.0601 Inf   -0.0586   0.17701
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)    0.2102 0.0696 Inf    0.0738   0.34664
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)     -0.1133 0.0565 Inf   -0.2239  -0.00262
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)       0.0377 0.0636 Inf   -0.0869   0.16235
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)   0.1510 0.0220 Inf    0.1079   0.19410
##  z.ratio p.value
##   11.116  <.0001
##    0.985  0.3246
##    3.020  0.0025
##   -2.006  0.0448
##    0.593  0.5529
##    6.870  <.0001
## 
## Degrees-of-freedom method: asymptotic 
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
##  contrast   estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10   -0.172 0.0155 Inf    -0.203    -0.142 -11.116  <.0001
##  diff_ESS1    -0.151 0.0220 Inf    -0.194    -0.108  -6.870  <.0001
## 
## Degrees-of-freedom method: asymptotic 
## Confidence level used: 0.95
```

``` r
# Test for whether and how much the differences have changed
pairs(diff_mod6,infer=c(T,T))
```

```
##  contrast               estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10 - diff_ESS1  -0.0215 0.0186 Inf    -0.058     0.015  -1.153  0.2490
## 
## Degrees-of-freedom method: asymptotic 
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
      obs_mean_wt=weighted.mean(x=uni.z,w=pspwght),
      obs_sd_wt=sqrt(wtd.var(uni.z,pspwght)),
      obs_mean=mean(uni.z),
      obs_sd=sd(uni.z),
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
  ylab("Mean-level of value universalism")+
  theme(legend.title=element_blank(),
        legend.text  = element_text(size = 16),
        axis.title.x = element_text(size = 18, face = "bold"),   
        axis.title.y = element_text(size = 18, face = "bold"),   
        axis.text.x  = element_text(size = 14),                  
        axis.text.y  = element_text(size = 14))

p_time_trends
```

![](Analysis_universalism_files/figure-html/unnamed-chunk-33-1.png)<!-- -->

``` r
png(filename = 
      "../results/uni/time_trends.png",
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
pred_cntry_dat$uni.z_mean<-predict(mod6,newdata=pred_cntry_dat)

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

range(pred_cntry_dat$uni.z_mean)
```

```
## [1] -0.9290826  0.6262593
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
      obs_mean_wt=weighted.mean(x=uni.z,w=pspwght),
      obs_sd_wt=sqrt(wtd.var(uni.z,pspwght)),
      obs_mean=mean(uni.z),
      obs_sd=sd(uni.z),
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

pdf("../results/uni/country_specific_time_trends.pdf", width = 8, height = 6)

for (ctry in countries) {
  print(
    ggplot(pred_cntry_dat[pred_cntry_dat$cntry == ctry, ], 
           aes(x = year, y = uni.z_mean, color = gender)) +
      geom_smooth(method = "lm",formula = "y ~ x", se = FALSE) +
      geom_point(size=8) +
      geom_point(aes(x=year,y=obs_mean_wt),size=8,shape = 1,alpha=.50)+
      geom_errorbar(aes(ymin=obs_mean_wt_LL, ymax=obs_mean_wt_UL),alpha = .50)+
      scale_color_manual(values = my_colors) +
      geom_flag(aes(country=tolower(ctry)))+
      scale_y_continuous(limits = c(-1.1, 1.1)) +
      scale_x_continuous(limits = c(2001, 2021), breaks = seq(2002, 2020, 2)) +
      ggtitle(paste("Country:", ctry))+
    ylab("Mean-level of value universalism")+
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
         aes(x = year, y = uni.z_mean, color = gender)) +
  geom_smooth(method = "lm", se = FALSE, formula="y~x") +
  geom_point(size=3.5) +
  geom_point(aes(x=year,y=obs_mean_wt),size=3.5,shape = 1,alpha=.70)+
  geom_errorbar(aes(ymin=obs_mean_wt_LL, ymax=obs_mean_wt_UL),alpha = .50)+
  scale_color_manual(values = my_colors) +
  #geom_flag(aes(country=tolower(ctry)))+
  scale_y_continuous(limits = c(-1.1, 1.1)) +
  scale_x_continuous(limits = c(2001, 2021), breaks = seq(2002, 2020, 2)) +
  #ggtitle(paste("Country:", ctry))+
  ylab("Mean-level of value universalism")+
  xlab("Year")+
  theme(legend.title=element_blank(),legend.position = "top",
        axis.text.x = element_text(angle = 45,size = 6,hjust=1))+
  facet_wrap(~CLDR,nrow=6,ncol=6)+
  #facet_wrap(~cntry,nrow=6,ncol=6)+
  geom_flag(aes(country=tolower(cntry)),size=2)

facet_plot
```

![](Analysis_universalism_files/figure-html/unnamed-chunk-35-1.png)<!-- -->

``` r
png(filename = 
      "../results/uni/country_time_trend_facets.png",
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
## [1] 17.69297
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
## 1   -0.23                0.22                   -0.02                    0.21                      0.23
## 2   -0.15                0.04                   -0.03                    0.03                      0.06
## 3   -0.14               -0.37                   -0.04                   -0.39                     -0.35
## 4   -0.19                0.01                   -0.07                   -0.03                      0.04
## 5   -0.04                0.49                   -0.09                    0.45                      0.54
## 6   -0.26               -0.28                    0.02                   -0.27                     -0.29
## 7   -0.20                0.36                    0.12                    0.42                      0.31
## 8   -0.20                0.32                   -0.06                    0.29                      0.35
## 9   -0.23                0.05                   -0.07                    0.02                      0.08
## 10  -0.09                0.27                   -0.03                    0.26                      0.29
## 11  -0.38                0.28                    0.06                    0.31                      0.25
## 12  -0.10               -0.07                   -0.09                   -0.11                     -0.02
## 13  -0.12                0.35                   -0.03                    0.33                      0.37
## 14   0.03               -0.73                   -0.01                   -0.74                     -0.73
## 15  -0.21                0.56                    0.01                    0.57                      0.55
## 16  -0.16               -0.41                    0.05                   -0.38                     -0.44
## 17  -0.15               -0.01                    0.00                   -0.01                     -0.01
## 18  -0.08               -0.10                   -0.01                   -0.11                     -0.10
## 19  -0.24                0.46                    0.01                    0.47                      0.45
## 20  -0.07               -0.53                   -0.07                   -0.56                     -0.49
## 21  -0.17               -0.35                   -0.09                   -0.39                     -0.30
## 22  -0.37                0.25                    0.10                    0.30                      0.21
## 23  -0.17               -0.15                    0.00                   -0.15                     -0.15
## 24  -0.18                0.16                   -0.05                    0.14                      0.19
## 25  -0.19                0.27                    0.08                    0.31                      0.23
## 26  -0.15                0.01                   -0.06                   -0.02                      0.04
## 27  -0.03                0.18                   -0.18                    0.09                      0.27
## 28  -0.13               -0.41                   -0.01                   -0.42                     -0.41
## 29  -0.29                0.69                    0.10                    0.75                      0.64
## 30  -0.16                0.54                    0.02                    0.55                      0.53
## 31  -0.16               -0.05                   -0.08                   -0.09                     -0.01
## 32   0.01               -0.26                   -0.08                   -0.30                     -0.22
## 33  -0.13               -0.22                   -0.11                   -0.27                     -0.16
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
## 1     GR               -0.73
## 2     HU               -0.41
## 3     RU               -0.41
## 4     BG               -0.37
## 5     LT               -0.35
## 6     CZ               -0.28
## 7     UA               -0.22
## 8     IL               -0.10
## 9     FR               -0.07
## 10    SK               -0.05
## 11    IE               -0.01
## 12    CH                0.01
## 13    PL                0.01
## 14    BE                0.04
## 15    EE                0.05
## 16    NL                0.16
## 17    PT                0.18
## 18    AT                0.22
## 19    ES                0.27
## 20    NO                0.27
## 21    FI                0.28
## 22    DK                0.32
## 23    GB                0.35
## 24    DE                0.36
## 25    IS                0.46
## 26    CY                0.49
## 27    SI                0.54
## 28    SE                0.69
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
## 1     PT                   -0.18
## 2     UA                   -0.11
## 3     CY                   -0.09
## 4     FR                   -0.09
## 5     LT                   -0.09
## 6     SK                   -0.08
## 7     CH                   -0.07
## 8     EE                   -0.07
## 9     DK                   -0.06
## 10    PL                   -0.06
## 11    NL                   -0.05
## 12    BG                   -0.04
## 13    BE                   -0.03
## 14    ES                   -0.03
## 15    GB                   -0.03
## 16    AT                   -0.02
## 17    GR                   -0.01
## 18    IL                   -0.01
## 19    RU                   -0.01
## 20    IE                    0.00
## 21    IS                    0.01
## 22    CZ                    0.02
## 23    SI                    0.02
## 24    HU                    0.05
## 25    FI                    0.06
## 26    NO                    0.08
## 27    SE                    0.10
## 28    DE                    0.12
```


#### Country-level predictors: GEI


``` r
mod6_GEI<-lmer(uni.z~gndr.c+essround.c+
             gndr.c:essround.c+
               gei.z.cm:gndr.c+gei.z.cm:essround.c+gei.z.cm:gndr.c:essround.c+
               (gndr.c+essround.c+gndr.c:essround.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_GEI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: uni.z ~ gndr.c + essround.c + gndr.c:essround.c + gei.z.cm:gndr.c +  
##     gei.z.cm:essround.c + gei.z.cm:gndr.c:essround.c + (gndr.c +  
##     essround.c + gndr.c:essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1271909.1 1272106.7 -635936.6 1271873.1    431760 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -8.3918 -0.5523  0.1058  0.6568  4.9654 
## 
## Random effects:
##  Groups   Name              Variance  Std.Dev. Corr             
##  cntry    (Intercept)       4.766e-02 0.218309                  
##           gndr.c            7.767e-03 0.088129  0.39            
##           essround.c        1.176e-03 0.034286  0.14 -0.29      
##           gndr.c:essround.c 8.019e-05 0.008955  0.00 -0.50  0.26
##  Residual                   9.854e-01 0.992679                  
## Number of obs: 431778, groups:  cntry, 32
## 
## Fixed effects:
##                              Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)                -0.0008702  0.0387521 32.2454063  -0.022  0.98222    
## gndr.c                     -0.1615720  0.0160665 30.4002619 -10.056 3.47e-11 ***
## essround.c                  0.0060555  0.0061575 31.6013772   0.983  0.33286    
## gndr.c:essround.c          -0.0012301  0.0020951 29.4644962  -0.587  0.56157    
## gndr.c:gei.z.cm            -0.0447880  0.0153264 33.1164524  -2.922  0.00622 ** 
## essround.c:gei.z.cm         0.0179490  0.0063130 33.7449307   2.843  0.00753 ** 
## gndr.c:essround.c:gei.z.cm -0.0007658  0.0023891 37.9399366  -0.321  0.75032    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c essrn. gnd.:. gn.:.. es.:..
## gndr.c       0.381                                   
## essround.c   0.127 -0.280                            
## gndr.c:ssr. -0.004 -0.389  0.202                     
## gndr.c:g.z.  0.001 -0.020  0.002 -0.033              
## essrnd.c:.. -0.004  0.001 -0.025 -0.004 -0.354       
## gndr.c:.:..  0.001 -0.014 -0.005 -0.182 -0.263  0.178
```

``` r
getFE(mod6_GEI,round=2,p.round = 5)
```

```
##                             Est.   SE    df      t       p    LL    UL
## (Intercept)                 0.00 0.04 32.25  -0.02 0.98222 -0.08  0.08
## gndr.c                     -0.16 0.02 30.40 -10.06 0.00000 -0.19 -0.13
## essround.c                  0.01 0.01 31.60   0.98 0.33286 -0.01  0.02
## gndr.c:essround.c           0.00 0.00 29.46  -0.59 0.56157 -0.01  0.00
## gndr.c:gei.z.cm            -0.04 0.02 33.12  -2.92 0.00622 -0.08 -0.01
## essround.c:gei.z.cm         0.02 0.01 33.74   2.84 0.00753  0.01  0.03
## gndr.c:essround.c:gei.z.cm  0.00 0.00 37.94  -0.32 0.75032 -0.01  0.00
```

``` r
getVC(mod6_GEI)
```

```
##         grp              var1              var2 sdcor vcov
## 1     cntry       (Intercept)              <NA>  0.22 0.05
## 2     cntry            gndr.c              <NA>  0.09 0.01
## 3     cntry        essround.c              <NA>  0.03 0.00
## 4     cntry gndr.c:essround.c              <NA>  0.01 0.00
## 5     cntry       (Intercept)            gndr.c  0.39 0.01
## 6     cntry       (Intercept)        essround.c  0.14 0.00
## 7     cntry       (Intercept) gndr.c:essround.c  0.00 0.00
## 8     cntry            gndr.c        essround.c -0.29 0.00
## 9     cntry            gndr.c gndr.c:essround.c -0.50 0.00
## 10    cntry        essround.c gndr.c:essround.c  0.26 0.00
## 11 Residual              <NA>              <NA>  0.99 0.99
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_GEI,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] 20.9207
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_GEI,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] 4.557249
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
## Note: D.f. calculations have been disabled because the number of observations exceeds 4e+05.
## To enable adjustments, add the argument 'lmerTest.limit = 431778' (or larger)
## [or, globally, 'set emm_options(lmerTest.limit = 431778)' or larger];
## but be warned that this may result in large computation time and memory use.
```

```
## NOTE: Results may be misleading due to involvement in interactions
```

``` r
change_mod6_GEI
```

```
## gei.z.cm = -1:
##  essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##         4.5 -0.0544 0.0583 Inf  -0.16859   0.05981  -0.933  0.3506
##        -4.5  0.0527 0.0532 Inf  -0.05171   0.15702   0.989  0.3228
## 
## gei.z.cm =  0:
##  essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##         4.5  0.0264 0.0504 Inf  -0.07243   0.12519   0.523  0.6008
##        -4.5 -0.0281 0.0447 Inf  -0.11572   0.05948  -0.629  0.5292
## 
## gei.z.cm =  1:
##  essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##         4.5  0.1072 0.0575 Inf  -0.00547   0.21977   1.865  0.0622
##        -4.5 -0.1089 0.0527 Inf  -0.21212  -0.00566  -2.068  0.0387
## 
## Degrees-of-freedom method: asymptotic 
## Confidence level used: 0.95
```

``` r
pairs(change_mod6_GEI,adjust="none",infer=c(T,T))
```

```
## gei.z.cm = -1:
##  contrast                         estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  essround.c4.5 - (essround.c-4.5)  -0.1070 0.0803 Inf   -0.2645    0.0504  -1.332  0.1827
## 
## gei.z.cm =  0:
##  contrast                         estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   0.0545 0.0554 Inf   -0.0541    0.1631   0.983  0.3254
## 
## gei.z.cm =  1:
##  contrast                         estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   0.2160 0.0784 Inf    0.0624    0.3697   2.756  0.0058
## 
## Degrees-of-freedom method: asymptotic 
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
```

```
## Note: D.f. calculations have been disabled because the number of observations exceeds 4e+05.
## To enable adjustments, add the argument 'lmerTest.limit = 431778' (or larger)
## [or, globally, 'set emm_options(lmerTest.limit = 431778)' or larger];
## but be warned that this may result in large computation time and memory use.
```

``` r
change_in_diff_mod6_GEI
```

```
## gei.z.cm = -1:
##  gndr.c essround.c    emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##    -0.5        4.5  0.005046 0.0589 Inf   -0.1104   0.12049   0.086  0.9317
##     0.5        4.5 -0.113828 0.0600 Inf   -0.2314   0.00374  -1.898  0.0577
##    -0.5       -4.5  0.109998 0.0494 Inf    0.0132   0.20683   2.226  0.0260
##     0.5       -4.5 -0.004697 0.0608 Inf   -0.1239   0.11448  -0.077  0.9384
## 
## gei.z.cm =  0:
##  gndr.c essround.c    emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##    -0.5        4.5  0.109933 0.0493 Inf    0.0132   0.20665   2.228  0.0259
##     0.5        4.5 -0.057174 0.0526 Inf   -0.1602   0.04584  -1.088  0.2767
##    -0.5       -4.5  0.049898 0.0412 Inf   -0.0309   0.13065   1.211  0.2259
##     0.5       -4.5 -0.106138 0.0503 Inf   -0.2047  -0.00755  -2.110  0.0348
## 
## gei.z.cm =  1:
##  gndr.c essround.c    emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##    -0.5        4.5  0.214821 0.0578 Inf    0.1014   0.32819   3.714  0.0002
##     0.5        4.5 -0.000521 0.0589 Inf   -0.1160   0.11499  -0.009  0.9929
##    -0.5       -4.5 -0.010201 0.0487 Inf   -0.1056   0.08521  -0.210  0.8340
##     0.5       -4.5 -0.207580 0.0601 Inf   -0.3254  -0.08977  -3.454  0.0006
## 
## Degrees-of-freedom method: asymptotic 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod6_GEI,adjust="none",infer=c(T,T))
```

```
## gei.z.cm = -1:
##  contrast                                                 estimate     SE  df asymp.LCL asymp.UCL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      0.11887 0.0235 Inf    0.0728    0.1650
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5) -0.10495 0.0791 Inf   -0.2600    0.0501
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   0.00974 0.0899 Inf   -0.1664    0.1859
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)    -0.22383 0.0764 Inf   -0.3736   -0.0741
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)     -0.10913 0.0845 Inf   -0.2747    0.0564
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  0.11469 0.0306 Inf    0.0548    0.1746
##  z.ratio p.value
##    5.052  <.0001
##   -1.327  0.1846
##    0.108  0.9137
##   -2.930  0.0034
##   -1.292  0.1964
##    3.752  0.0002
## 
## gei.z.cm =  0:
##  contrast                                                 estimate     SE  df asymp.LCL asymp.UCL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      0.16711 0.0151 Inf    0.1374    0.1968
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  0.06004 0.0543 Inf   -0.0464    0.1665
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   0.21607 0.0619 Inf    0.0948    0.3373
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)    -0.10707 0.0532 Inf   -0.2113   -0.0028
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      0.04896 0.0581 Inf   -0.0648    0.1628
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  0.15604 0.0216 Inf    0.1138    0.1983
##  z.ratio p.value
##   11.037  <.0001
##    1.106  0.2689
##    3.492  0.0005
##   -2.013  0.0442
##    0.843  0.3990
##    7.238  <.0001
## 
## gei.z.cm =  1:
##  contrast                                                 estimate     SE  df asymp.LCL asymp.UCL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      0.21534 0.0208 Inf    0.1746    0.2561
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  0.22502 0.0768 Inf    0.0746    0.3755
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   0.42240 0.0880 Inf    0.2500    0.5948
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)     0.00968 0.0743 Inf   -0.1359    0.1552
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      0.20706 0.0820 Inf    0.0463    0.3678
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  0.19738 0.0295 Inf    0.1396    0.2552
##  z.ratio p.value
##   10.360  <.0001
##    2.931  0.0034
##    4.801  <.0001
##    0.130  0.8963
##    2.524  0.0116
##    6.693  <.0001
## 
## Degrees-of-freedom method: asymptotic 
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
##  contrast   estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10   -0.119 0.0235 Inf    -0.165   -0.0728  -5.052  <.0001
##  diff_ESS1    -0.115 0.0306 Inf    -0.175   -0.0548  -3.752  0.0002
## 
## gei.z.cm =  0:
##  contrast   estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10   -0.167 0.0151 Inf    -0.197   -0.1374 -11.037  <.0001
##  diff_ESS1    -0.156 0.0216 Inf    -0.198   -0.1138  -7.238  <.0001
## 
## gei.z.cm =  1:
##  contrast   estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10   -0.215 0.0208 Inf    -0.256   -0.1746 -10.360  <.0001
##  diff_ESS1    -0.197 0.0295 Inf    -0.255   -0.1396  -6.693  <.0001
## 
## Degrees-of-freedom method: asymptotic 
## Confidence level used: 0.95
```

``` r
# Test for whether and how much the differences have changed
pairs(diff_mod6_GEI,infer=c(T,T))
```

```
## gei.z.cm = -1:
##  contrast               estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10 - diff_ESS1 -0.00418 0.0311 Inf   -0.0651    0.0567  -0.134  0.8930
## 
## gei.z.cm =  0:
##  contrast               estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10 - diff_ESS1 -0.01107 0.0189 Inf   -0.0480    0.0259  -0.587  0.5571
## 
## gei.z.cm =  1:
##  contrast               estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10 - diff_ESS1 -0.01796 0.0259 Inf   -0.0687    0.0328  -0.694  0.4877
## 
## Degrees-of-freedom method: asymptotic 
## Confidence level used: 0.95
```


#### Country-level predictors: GGGI


``` r
mod6_GGGI<-lmer(uni.z~gndr.c+essround.c+
             gndr.c:essround.c+
               gggi.z.cm:gndr.c+gggi.z.cm:essround.c+gggi.z.cm:gndr.c:essround.c+
               (gndr.c+essround.c+gndr.c:essround.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
```

```
## Warning in checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv, : Model failed to converge
## with max|grad| = 0.00238441 (tol = 0.002, component 1)
```

``` r
summary(mod6_GGGI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: uni.z ~ gndr.c + essround.c + gndr.c:essround.c + gggi.z.cm:gndr.c +  
##     gggi.z.cm:essround.c + gggi.z.cm:gndr.c:essround.c + (gndr.c +  
##     essround.c + gndr.c:essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##  932795.8  932987.6 -466379.9  932759.8    314628 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -8.5350 -0.5526  0.1064  0.6570  4.7620 
## 
## Random effects:
##  Groups   Name              Variance Std.Dev. Corr             
##  cntry    (Intercept)       0.055756 0.23613                   
##           gndr.c            0.005482 0.07404   0.37            
##           essround.c        0.001804 0.04247   0.22 -0.20      
##           gndr.c:essround.c 0.000108 0.01039  -0.03  0.11 -0.41
##  Residual                   0.997727 0.99886                   
## Number of obs: 314646, groups:  cntry, 33
## 
## Fixed effects:
##                              Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)                 -0.001532   0.041433 33.199410  -0.037  0.97072    
## gndr.c                      -0.158722   0.013653 29.614561 -11.626 1.46e-12 ***
## essround.c                   0.001443   0.007776 29.596298   0.186  0.85409    
## gndr.c:essround.c           -0.003522   0.002658 28.234349  -1.325  0.19581    
## gndr.c:gggi.z.cm            -0.060168   0.013394 33.015755  -4.492 8.16e-05 ***
## essround.c:gggi.z.cm         0.026781   0.008226 30.991994   3.256  0.00274 ** 
## gndr.c:essround.c:gggi.z.cm  0.007689   0.003006 31.961507   2.558  0.01547 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c essrn. gnd.:. gn.:.. es.:..
## gndr.c       0.349                                   
## essround.c   0.214 -0.191                            
## gndr.c:ssr. -0.028  0.031 -0.284                     
## gndr.c:gg..  0.005 -0.015  0.014 -0.027              
## essrnd.c:.. -0.014  0.008 -0.076  0.018 -0.289       
## gndr.c:.:..  0.007 -0.025  0.018 -0.099  0.035 -0.278
## optimizer (bobyqa) convergence code: 0 (OK)
## Model failed to converge with max|grad| = 0.00238441 (tol = 0.002, component 1)
```

``` r
getFE(mod6_GGGI,round=2,p.round = 5)
```

```
##                              Est.   SE    df      t       p    LL    UL
## (Intercept)                  0.00 0.04 33.20  -0.04 0.97072 -0.09  0.08
## gndr.c                      -0.16 0.01 29.61 -11.63 0.00000 -0.19 -0.13
## essround.c                   0.00 0.01 29.60   0.19 0.85409 -0.01  0.02
## gndr.c:essround.c            0.00 0.00 28.23  -1.32 0.19581 -0.01  0.00
## gndr.c:gggi.z.cm            -0.06 0.01 33.02  -4.49 0.00008 -0.09 -0.03
## essround.c:gggi.z.cm         0.03 0.01 30.99   3.26 0.00274  0.01  0.04
## gndr.c:essround.c:gggi.z.cm  0.01 0.00 31.96   2.56 0.01547  0.00  0.01
```

``` r
getVC(mod6_GGGI)
```

```
##         grp              var1              var2 sdcor vcov
## 1     cntry       (Intercept)              <NA>  0.24 0.06
## 2     cntry            gndr.c              <NA>  0.07 0.01
## 3     cntry        essround.c              <NA>  0.04 0.00
## 4     cntry gndr.c:essround.c              <NA>  0.01 0.00
## 5     cntry       (Intercept)            gndr.c  0.37 0.01
## 6     cntry       (Intercept)        essround.c  0.22 0.00
## 7     cntry       (Intercept) gndr.c:essround.c -0.03 0.00
## 8     cntry            gndr.c        essround.c -0.20 0.00
## 9     cntry            gndr.c gndr.c:essround.c  0.11 0.00
## 10    cntry        essround.c gndr.c:essround.c -0.41 0.00
## 11 Residual              <NA>              <NA>  1.00 1.00
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_GGGI,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] -21.34537
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_GGGI,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] -28.571
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
##  essround.c   emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5 -0.11556 0.0719 42.7  -0.2606  0.02949  -1.607  0.1154
##        -4.5  0.11249 0.0620 41.4  -0.0127  0.23770   1.814  0.0769
## 
## gggi.z.cm =  0:
##  essround.c   emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.00496 0.0597 31.5  -0.1167  0.12660   0.083  0.9343
##        -4.5 -0.00802 0.0482 31.9  -0.1062  0.09011  -0.167  0.8688
## 
## gggi.z.cm =  1:
##  essround.c   emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.12548 0.0685 41.1  -0.0129  0.26381   1.832  0.0743
##        -4.5 -0.12854 0.0595 39.8  -0.2487 -0.00833  -2.162  0.0367
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
##  essround.c4.5 - (essround.c-4.5)   -0.228 0.1060 31.5  -0.4434  -0.0127  -2.158  0.0386
## 
## gggi.z.cm =  0:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)    0.013 0.0700 29.6  -0.1300   0.1560   0.186  0.8541
## 
## gggi.z.cm =  1:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)    0.254 0.0979 30.3   0.0541   0.4539   2.594  0.0145
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
##  gndr.c essround.c   emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5        4.5 -0.04106 0.0755 41.6  -0.1935   0.1114  -0.543  0.5897
##     0.5        4.5 -0.19006 0.0709 41.4  -0.3331  -0.0470  -2.682  0.0105
##    -0.5       -4.5  0.13655 0.0619 41.4   0.0115   0.2616   2.205  0.0331
##     0.5       -4.5  0.08844 0.0648 40.5  -0.0425   0.2194   1.365  0.1799
## 
## gggi.z.cm =  0:
##  gndr.c essround.c   emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5        4.5  0.09224 0.0606 31.4  -0.0313   0.2158   1.522  0.1380
##     0.5        4.5 -0.08232 0.0602 30.8  -0.2051   0.0404  -1.368  0.1811
##    -0.5       -4.5  0.06341 0.0471 32.4  -0.0324   0.1593   1.347  0.1874
##     0.5       -4.5 -0.07946 0.0508 31.3  -0.1831   0.0242  -1.563  0.1281
## 
## gggi.z.cm =  1:
##  gndr.c essround.c   emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5        4.5  0.22554 0.0716 40.4   0.0808   0.3703   3.148  0.0031
##     0.5        4.5  0.02541 0.0677 39.4  -0.1115   0.1623   0.375  0.7094
##    -0.5       -4.5 -0.00972 0.0592 40.1  -0.1294   0.1100  -0.164  0.8705
##     0.5       -4.5 -0.24736 0.0624 38.8  -0.3736  -0.1212  -3.966  0.0003
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
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      0.14900 0.0277 32.5  0.09252  0.20548   5.370
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5) -0.17760 0.1120 31.3 -0.40665  0.05145  -1.581
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)  -0.12950 0.1120 30.8 -0.35774  0.09874  -1.157
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)    -0.32660 0.1030 32.2 -0.53585 -0.11735  -3.179
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)     -0.27850 0.1020 30.8 -0.48680 -0.07020  -2.728
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  0.04811 0.0262 35.1 -0.00516  0.10137   1.833
##  p.value
##   <.0001
##   0.1240
##   0.2560
##   0.0033
##   0.0104
##   0.0753
## 
## gggi.z.cm =  0:
##  contrast                                                 estimate     SE   df lower.CL upper.CL t.ratio
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      0.17457 0.0184 29.0  0.13688  0.21226   9.473
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  0.02883 0.0743 29.3 -0.12300  0.18066   0.388
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   0.17170 0.0738 29.1  0.02075  0.32266   2.326
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)    -0.14574 0.0687 29.3 -0.28617 -0.00531  -2.122
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)     -0.00286 0.0676 28.7 -0.14110  0.13538  -0.042
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  0.14287 0.0179 28.6  0.10631  0.17944   7.995
##  p.value
##   <.0001
##   0.7007
##   0.0272
##   0.0424
##   0.9665
##   <.0001
## 
## gggi.z.cm =  1:
##  contrast                                                 estimate     SE   df lower.CL upper.CL t.ratio
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      0.20014 0.0257 29.7  0.14766  0.25261   7.793
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  0.23526 0.1040 30.2  0.02256  0.44796   2.258
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   0.47291 0.1040 29.2  0.25995  0.68586   4.541
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)     0.03513 0.0952 31.1 -0.15893  0.22919   0.369
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      0.27277 0.0944 29.4  0.07974  0.46580   2.888
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  0.23764 0.0255 30.6  0.18564  0.28965   9.325
##  p.value
##   <.0001
##   0.0313
##   0.0001
##   0.7145
##   0.0072
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
##  diff_ESS10  -0.1490 0.0277 32.5   -0.205 -0.09252  -5.370  <.0001
##  diff_ESS1   -0.0481 0.0262 35.1   -0.101  0.00516  -1.833  0.0753
## 
## gggi.z.cm =  0:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10  -0.1746 0.0184 29.0   -0.212 -0.13688  -9.473  <.0001
##  diff_ESS1   -0.1429 0.0179 28.6   -0.179 -0.10631  -7.995  <.0001
## 
## gggi.z.cm =  1:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10  -0.2001 0.0257 29.7   -0.253 -0.14766  -7.793  <.0001
##  diff_ESS1   -0.2376 0.0255 30.6   -0.290 -0.18564  -9.325  <.0001
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
##  diff_ESS10 - diff_ESS1  -0.1009 0.0378 32.0  -0.1780  -0.0238  -2.667  0.0119
## 
## gggi.z.cm =  0:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10 - diff_ESS1  -0.0317 0.0239 28.2  -0.0807   0.0173  -1.325  0.1958
## 
## gggi.z.cm =  1:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10 - diff_ESS1   0.0375 0.0343 28.2  -0.0327   0.1077   1.094  0.2834
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

#### Country-level predictors: GDI


``` r
mod6_GDI<-lmer(uni.z~gndr.c+essround.c+
             gndr.c:essround.c+
               gdi.z.cm:gndr.c+gdi.z.cm:essround.c+gdi.z.cm:gndr.c:essround.c+
               (gndr.c+essround.c+gndr.c:essround.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_GDI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: uni.z ~ gndr.c + essround.c + gndr.c:essround.c + gdi.z.cm:gndr.c +  
##     gdi.z.cm:essround.c + gdi.z.cm:gndr.c:essround.c + (gndr.c +  
##     essround.c + gndr.c:essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1307030.8 1307228.7 -653497.4 1306994.8    441150 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -8.4353 -0.5506  0.1050  0.6540  4.9389 
## 
## Random effects:
##  Groups   Name              Variance  Std.Dev. Corr             
##  cntry    (Intercept)       5.023e-02 0.224128                  
##           gndr.c            7.276e-03 0.085300  0.04            
##           essround.c        1.534e-03 0.039165  0.32 -0.45      
##           gndr.c:essround.c 7.932e-05 0.008906  0.04 -0.53  0.23
##  Residual                   9.976e-01 0.998789                  
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##                             Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)                -0.011886   0.039168 33.223393  -0.303   0.7634    
## gndr.c                     -0.161232   0.015342 32.677108 -10.509 5.17e-12 ***
## essround.c                  0.005342   0.006902 30.839440   0.774   0.4448    
## gndr.c:essround.c          -0.002344   0.002030 27.927190  -1.155   0.2580    
## gndr.c:gdi.z.cm            -0.038893   0.015837 35.033970  -2.456   0.0192 *  
## essround.c:gdi.z.cm         0.007313   0.006755 34.900646   1.083   0.2864    
## gndr.c:essround.c:gdi.z.cm -0.001138   0.002465 42.665389  -0.462   0.6467    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c essrn. gnd.:. gn.:.. es.:..
## gndr.c       0.043                                   
## essround.c   0.314 -0.434                            
## gndr.c:ssr.  0.027 -0.420  0.177                     
## gndr.c:gd..  0.000 -0.011  0.001 -0.002              
## essrnd.c:..  0.002  0.002 -0.018  0.000 -0.457       
## gndr.c:.:..  0.003 -0.006  0.000  0.001 -0.282  0.146
```

``` r
getFE(mod6_GDI,round=2,p.round = 5)
```

```
##                             Est.   SE    df      t       p    LL    UL
## (Intercept)                -0.01 0.04 33.22  -0.30 0.76342 -0.09  0.07
## gndr.c                     -0.16 0.02 32.68 -10.51 0.00000 -0.19 -0.13
## essround.c                  0.01 0.01 30.84   0.77 0.44484 -0.01  0.02
## gndr.c:essround.c           0.00 0.00 27.93  -1.15 0.25805 -0.01  0.00
## gndr.c:gdi.z.cm            -0.04 0.02 35.03  -2.46 0.01916 -0.07 -0.01
## essround.c:gdi.z.cm         0.01 0.01 34.90   1.08 0.28636 -0.01  0.02
## gndr.c:essround.c:gdi.z.cm  0.00 0.00 42.67  -0.46 0.64672 -0.01  0.00
```

``` r
getVC(mod6_GDI)
```

```
##         grp              var1              var2 sdcor vcov
## 1     cntry       (Intercept)              <NA>  0.22 0.05
## 2     cntry            gndr.c              <NA>  0.09 0.01
## 3     cntry        essround.c              <NA>  0.04 0.00
## 4     cntry gndr.c:essround.c              <NA>  0.01 0.00
## 5     cntry       (Intercept)            gndr.c  0.04 0.00
## 6     cntry       (Intercept)        essround.c  0.32 0.00
## 7     cntry       (Intercept) gndr.c:essround.c  0.04 0.00
## 8     cntry            gndr.c        essround.c -0.45 0.00
## 9     cntry            gndr.c gndr.c:essround.c -0.53 0.00
## 10    cntry        essround.c gndr.c:essround.c  0.23 0.00
## 11 Residual              <NA>              <NA>  1.00 1.00
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_GDI,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] -3.188515
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_GDI,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] 5.59408
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
## Note: D.f. calculations have been disabled because the number of observations exceeds 4e+05.
## To enable adjustments, add the argument 'lmerTest.limit = 441168' (or larger)
## [or, globally, 'set emm_options(lmerTest.limit = 441168)' or larger];
## but be warned that this may result in large computation time and memory use.
```

```
## NOTE: Results may be misleading due to involvement in interactions
```

``` r
change_mod6_GDI
```

```
## gdi.z.cm = -1:
##  essround.c   emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##         4.5 -0.02076 0.0649 Inf   -0.1480    0.1065  -0.320  0.7492
##        -4.5 -0.00302 0.0520 Inf   -0.1049    0.0988  -0.058  0.9537
## 
## gdi.z.cm =  0:
##  essround.c   emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##         4.5  0.01215 0.0571 Inf   -0.0998    0.1241   0.213  0.8315
##        -4.5 -0.03593 0.0417 Inf   -0.1176    0.0457  -0.862  0.3884
## 
## gdi.z.cm =  1:
##  essround.c   emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##         4.5  0.04506 0.0645 Inf   -0.0813    0.1714   0.699  0.4846
##        -4.5 -0.06884 0.0512 Inf   -0.1691    0.0315  -1.345  0.1785
## 
## Degrees-of-freedom method: asymptotic 
## Confidence level used: 0.95
```

``` r
pairs(change_mod6_GDI,adjust="none",infer=c(T,T))
```

```
## gdi.z.cm = -1:
##  contrast                         estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  essround.c4.5 - (essround.c-4.5)  -0.0177 0.0877 Inf   -0.1897     0.154  -0.202  0.8397
## 
## gdi.z.cm =  0:
##  contrast                         estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   0.0481 0.0621 Inf   -0.0737     0.170   0.774  0.4389
## 
## gdi.z.cm =  1:
##  contrast                         estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   0.1139 0.0861 Inf   -0.0549     0.283   1.323  0.1859
## 
## Degrees-of-freedom method: asymptotic 
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
```

```
## Note: D.f. calculations have been disabled because the number of observations exceeds 4e+05.
## To enable adjustments, add the argument 'lmerTest.limit = 441168' (or larger)
## [or, globally, 'set emm_options(lmerTest.limit = 441168)' or larger];
## but be warned that this may result in large computation time and memory use.
```

``` r
change_in_diff_mod6_GDI
```

```
## gdi.z.cm = -1:
##  gndr.c essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##    -0.5        4.5  0.0431 0.0681 Inf   -0.0903    0.1765   0.634  0.5262
##     0.5        4.5 -0.0846 0.0636 Inf   -0.2092    0.0400  -1.331  0.1831
##    -0.5       -4.5  0.0554 0.0488 Inf   -0.0402    0.1511   1.136  0.2560
##     0.5       -4.5 -0.0615 0.0589 Inf   -0.1770    0.0540  -1.043  0.2969
## 
## gdi.z.cm =  0:
##  gndr.c essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##    -0.5        4.5  0.0980 0.0586 Inf   -0.0168    0.2129   1.673  0.0943
##     0.5        4.5 -0.0737 0.0565 Inf   -0.1845    0.0370  -1.305  0.1919
##    -0.5       -4.5  0.0394 0.0396 Inf   -0.0383    0.1171   0.994  0.3200
##     0.5       -4.5 -0.1113 0.0460 Inf   -0.2015   -0.0211  -2.418  0.0156
## 
## gdi.z.cm =  1:
##  gndr.c essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##    -0.5        4.5  0.1530 0.0676 Inf    0.0206    0.2854   2.264  0.0236
##     0.5        4.5 -0.0628 0.0631 Inf   -0.1866    0.0609  -0.995  0.3196
##    -0.5       -4.5  0.0234 0.0480 Inf   -0.0707    0.1175   0.487  0.6261
##     0.5       -4.5 -0.1611 0.0582 Inf   -0.2751   -0.0470  -2.768  0.0056
## 
## Degrees-of-freedom method: asymptotic 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod6_GDI,adjust="none",infer=c(T,T))
```

```
## gdi.z.cm = -1:
##  contrast                                                 estimate     SE  df asymp.LCL asymp.UCL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5       0.1278 0.0220 Inf    0.0846   0.17089
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  -0.0123 0.0866 Inf   -0.1821   0.15745
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)    0.1046 0.0995 Inf   -0.0904   0.29961
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)     -0.1401 0.0805 Inf   -0.2978   0.01760
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      -0.0232 0.0911 Inf   -0.2017   0.15538
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)   0.1169 0.0302 Inf    0.0577   0.17610
##  z.ratio p.value
##    5.807  <.0001
##   -0.142  0.8870
##    1.051  0.2931
##   -1.741  0.0816
##   -0.254  0.7992
##    3.872  0.0001
## 
## gdi.z.cm =  0:
##  contrast                                                 estimate     SE  df asymp.LCL asymp.UCL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5       0.1718 0.0142 Inf    0.1440   0.19958
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   0.0586 0.0612 Inf   -0.0613   0.17852
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)    0.2093 0.0702 Inf    0.0718   0.34681
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)     -0.1132 0.0572 Inf   -0.2252  -0.00112
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)       0.0375 0.0644 Inf   -0.0886   0.16369
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)   0.1507 0.0209 Inf    0.1097   0.19163
##  z.ratio p.value
##   12.112  <.0001
##    0.958  0.3378
##    2.984  0.0028
##   -1.980  0.0478
##    0.583  0.5598
##    7.212  <.0001
## 
## gdi.z.cm =  1:
##  contrast                                                 estimate     SE  df asymp.LCL asymp.UCL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5       0.2158 0.0216 Inf    0.1734   0.25820
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   0.1296 0.0850 Inf   -0.0370   0.29616
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)    0.3140 0.0980 Inf    0.1220   0.50603
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)     -0.0862 0.0787 Inf   -0.2405   0.06805
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)       0.0982 0.0896 Inf   -0.0773   0.27375
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)   0.1845 0.0301 Inf    0.1254   0.24347
##  z.ratio p.value
##    9.974  <.0001
##    1.524  0.1274
##    3.206  0.0013
##   -1.095  0.2733
##    1.097  0.2727
##    6.126  <.0001
## 
## Degrees-of-freedom method: asymptotic 
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
##  contrast   estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10   -0.128 0.0220 Inf    -0.171   -0.0846  -5.807  <.0001
##  diff_ESS1    -0.117 0.0302 Inf    -0.176   -0.0577  -3.872  0.0001
## 
## gdi.z.cm =  0:
##  contrast   estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10   -0.172 0.0142 Inf    -0.200   -0.1440 -12.112  <.0001
##  diff_ESS1    -0.151 0.0209 Inf    -0.192   -0.1097  -7.212  <.0001
## 
## gdi.z.cm =  1:
##  contrast   estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10   -0.216 0.0216 Inf    -0.258   -0.1734  -9.974  <.0001
##  diff_ESS1    -0.184 0.0301 Inf    -0.243   -0.1254  -6.126  <.0001
## 
## Degrees-of-freedom method: asymptotic 
## Confidence level used: 0.95
```

``` r
# Test for whether and how much the differences have changed
pairs(diff_mod6_GDI,infer=c(T,T))
```

```
## gdi.z.cm = -1:
##  contrast               estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10 - diff_ESS1  -0.0109 0.0287 Inf   -0.0672    0.0454  -0.378  0.7055
## 
## gdi.z.cm =  0:
##  contrast               estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10 - diff_ESS1  -0.0211 0.0183 Inf   -0.0569    0.0147  -1.155  0.2483
## 
## gdi.z.cm =  1:
##  contrast               estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10 - diff_ESS1  -0.0313 0.0288 Inf   -0.0877    0.0250  -1.090  0.2758
## 
## Degrees-of-freedom method: asymptotic 
## Confidence level used: 0.95
```

#### Country-level predictors: log_GDP


``` r
mod6_log_GDP<-lmer(uni.z~gndr.c+essround.c+
             gndr.c:essround.c+
               log_gdp.z.cm:gndr.c+log_gdp.z.cm:essround.c+log_gdp.z.cm:gndr.c:essround.c+
               (gndr.c+essround.c+gndr.c:essround.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_log_GDP)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: uni.z ~ gndr.c + essround.c + gndr.c:essround.c + log_gdp.z.cm:gndr.c +  
##     log_gdp.z.cm:essround.c + log_gdp.z.cm:gndr.c:essround.c +  
##     (gndr.c + essround.c + gndr.c:essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1307029.7 1307227.6 -653496.8 1306993.7    441150 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -8.4330 -0.5508  0.1051  0.6545  4.9351 
## 
## Random effects:
##  Groups   Name              Variance  Std.Dev. Corr             
##  cntry    (Intercept)       5.032e-02 0.224313                  
##           gndr.c            8.166e-03 0.090363  0.38            
##           essround.c        1.232e-03 0.035101  0.10 -0.36      
##           gndr.c:essround.c 7.573e-05 0.008702 -0.03 -0.44  0.14
##  Residual                   9.976e-01 0.998790                  
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##                                 Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)                    -0.010782   0.039204 33.234815  -0.275   0.7850    
## gndr.c                         -0.161780   0.016195 30.747838  -9.989 3.59e-11 ***
## essround.c                      0.005351   0.006205 33.084441   0.862   0.3946    
## gndr.c:essround.c              -0.002577   0.002042 30.355114  -1.262   0.2167    
## gndr.c:log_gdp.z.cm            -0.033314   0.015288 33.071340  -2.179   0.0366 *  
## essround.c:log_gdp.z.cm         0.015366   0.006319 34.802882   2.432   0.0203 *  
## gndr.c:essround.c:log_gdp.z.cm  0.002603   0.002194 33.395674   1.186   0.2440    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c essrn. gnd.:. g.:_.. e.:_..
## gndr.c       0.371                                   
## essround.c   0.096 -0.347                            
## gndr.c:ssr. -0.027 -0.345  0.110                     
## gndr.c:l_..  0.001 -0.005 -0.009  0.006              
## essrnd.:_..  0.011 -0.005 -0.011 -0.009 -0.400       
## gndr.:.:_..  0.005  0.009 -0.007 -0.183 -0.336  0.113
```

``` r
getFE(mod6_log_GDP,round=2,p.round = 5)
```

```
##                                 Est.   SE    df     t       p    LL    UL
## (Intercept)                    -0.01 0.04 33.23 -0.28 0.78499 -0.09  0.07
## gndr.c                         -0.16 0.02 30.75 -9.99 0.00000 -0.19 -0.13
## essround.c                      0.01 0.01 33.08  0.86 0.39463 -0.01  0.02
## gndr.c:essround.c               0.00 0.00 30.36 -1.26 0.21672 -0.01  0.00
## gndr.c:log_gdp.z.cm            -0.03 0.02 33.07 -2.18 0.03656 -0.06  0.00
## essround.c:log_gdp.z.cm         0.02 0.01 34.80  2.43 0.02031  0.00  0.03
## gndr.c:essround.c:log_gdp.z.cm  0.00 0.00 33.40  1.19 0.24396  0.00  0.01
```

``` r
getVC(mod6_log_GDP)
```

```
##         grp              var1              var2 sdcor vcov
## 1     cntry       (Intercept)              <NA>  0.22 0.05
## 2     cntry            gndr.c              <NA>  0.09 0.01
## 3     cntry        essround.c              <NA>  0.04 0.00
## 4     cntry gndr.c:essround.c              <NA>  0.01 0.00
## 5     cntry       (Intercept)            gndr.c  0.38 0.01
## 6     cntry       (Intercept)        essround.c  0.10 0.00
## 7     cntry       (Intercept) gndr.c:essround.c -0.03 0.00
## 8     cntry            gndr.c        essround.c -0.36 0.00
## 9     cntry            gndr.c gndr.c:essround.c -0.44 0.00
## 10    cntry        essround.c gndr.c:essround.c  0.14 0.00
## 11 Residual              <NA>              <NA>  1.00 1.00
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_log_GDP,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] 17.11628
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_log_GDP,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] 9.864478
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
## Note: D.f. calculations have been disabled because the number of observations exceeds 4e+05.
## To enable adjustments, add the argument 'lmerTest.limit = 441168' (or larger)
## [or, globally, 'set emm_options(lmerTest.limit = 441168)' or larger];
## but be warned that this may result in large computation time and memory use.
```

```
## NOTE: Results may be misleading due to involvement in interactions
```

``` r
change_mod6_log_GDP
```

```
## log_gdp.z.cm = -1:
##  essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##         4.5 -0.0558 0.0577 Inf   -0.1689   0.05723  -0.968  0.3330
##        -4.5  0.0343 0.0544 Inf   -0.0723   0.14086   0.630  0.5284
## 
## log_gdp.z.cm =  0:
##  essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##         4.5  0.0133 0.0503 Inf   -0.0852   0.11180   0.265  0.7913
##        -4.5 -0.0349 0.0459 Inf   -0.1248   0.05510  -0.760  0.4475
## 
## log_gdp.z.cm =  1:
##  essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##         4.5  0.0824 0.0578 Inf   -0.0308   0.19573   1.427  0.1537
##        -4.5 -0.1040 0.0536 Inf   -0.2091   0.00106  -1.940  0.0524
## 
## Degrees-of-freedom method: asymptotic 
## Confidence level used: 0.95
```

``` r
pairs(change_mod6_log_GDP,adjust="none",infer=c(T,T))
```

```
## log_gdp.z.cm = -1:
##  contrast                         estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  essround.c4.5 - (essround.c-4.5)  -0.0901 0.0801 Inf   -0.2472     0.067  -1.125  0.2608
## 
## log_gdp.z.cm =  0:
##  contrast                         estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   0.0482 0.0558 Inf   -0.0613     0.158   0.862  0.3884
## 
## log_gdp.z.cm =  1:
##  contrast                         estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   0.1865 0.0793 Inf    0.0311     0.342   2.353  0.0186
## 
## Degrees-of-freedom method: asymptotic 
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
```

```
## Note: D.f. calculations have been disabled because the number of observations exceeds 4e+05.
## To enable adjustments, add the argument 'lmerTest.limit = 441168' (or larger)
## [or, globally, 'set emm_options(lmerTest.limit = 441168)' or larger];
## but be warned that this may result in large computation time and memory use.
```

``` r
change_in_diff_mod6_log_GDP
```

```
## log_gdp.z.cm = -1:
##  gndr.c essround.c   emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##    -0.5        4.5  0.02004 0.0591 Inf  -0.09579   0.13587   0.339  0.7346
##     0.5        4.5 -0.13174 0.0585 Inf  -0.24631  -0.01716  -2.253  0.0242
##    -0.5       -4.5  0.08686 0.0505 Inf  -0.01217   0.18590   1.719  0.0856
##     0.5       -4.5 -0.01829 0.0618 Inf  -0.13947   0.10288  -0.296  0.7673
## 
## log_gdp.z.cm =  0:
##  gndr.c essround.c   emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##    -0.5        4.5  0.09999 0.0499 Inf   0.00218   0.19779   2.004  0.0451
##     0.5        4.5 -0.07339 0.0518 Inf  -0.17493   0.02815  -1.417  0.1566
##    -0.5       -4.5  0.04023 0.0423 Inf  -0.04268   0.12314   0.951  0.3416
##     0.5       -4.5 -0.10996 0.0515 Inf  -0.21084  -0.00907  -2.136  0.0327
## 
## log_gdp.z.cm =  1:
##  gndr.c essround.c   emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##    -0.5        4.5  0.17994 0.0592 Inf   0.06392   0.29595   3.040  0.0024
##     0.5        4.5 -0.01504 0.0583 Inf  -0.12930   0.09921  -0.258  0.7964
##    -0.5       -4.5 -0.00641 0.0495 Inf  -0.10338   0.09057  -0.129  0.8970
##     0.5       -4.5 -0.20162 0.0610 Inf  -0.32118  -0.08206  -3.305  0.0009
## 
## Degrees-of-freedom method: asymptotic 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod6_log_GDP,adjust="none",infer=c(T,T))
```

```
## log_gdp.z.cm = -1:
##  contrast                                                 estimate     SE  df asymp.LCL asymp.UCL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      0.15177 0.0225 Inf    0.1077    0.1958
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5) -0.06683 0.0799 Inf   -0.2234    0.0898
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   0.03833 0.0907 Inf   -0.1394    0.2160
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)    -0.21860 0.0750 Inf   -0.3656   -0.0716
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)     -0.11344 0.0830 Inf   -0.2762    0.0493
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  0.10516 0.0304 Inf    0.0456    0.1647
##  z.ratio p.value
##    6.754  <.0001
##   -0.837  0.4029
##    0.423  0.6725
##   -2.915  0.0036
##   -1.366  0.1719
##    3.463  0.0005
## 
## log_gdp.z.cm =  0:
##  contrast                                                 estimate     SE  df asymp.LCL asymp.UCL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      0.17337 0.0156 Inf    0.1428    0.2040
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  0.05976 0.0556 Inf   -0.0492    0.1687
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   0.20994 0.0633 Inf    0.0859    0.3340
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)    -0.11362 0.0525 Inf   -0.2165   -0.0108
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      0.03657 0.0576 Inf   -0.0763    0.1494
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  0.15019 0.0212 Inf    0.1086    0.1917
##  z.ratio p.value
##   11.099  <.0001
##    1.075  0.2824
##    3.316  0.0009
##   -2.166  0.0303
##    0.635  0.5254
##    7.084  <.0001
## 
## log_gdp.z.cm =  1:
##  contrast                                                 estimate     SE  df asymp.LCL asymp.UCL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      0.19498 0.0210 Inf    0.1537    0.2362
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  0.18634 0.0788 Inf    0.0319    0.3408
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   0.38155 0.0901 Inf    0.2049    0.5582
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)    -0.00863 0.0737 Inf   -0.1530    0.1358
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      0.18658 0.0816 Inf    0.0267    0.3464
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  0.19521 0.0290 Inf    0.1383    0.2521
##  z.ratio p.value
##    9.267  <.0001
##    2.365  0.0180
##    4.234  <.0001
##   -0.117  0.9067
##    2.287  0.0222
##    6.728  <.0001
## 
## Degrees-of-freedom method: asymptotic 
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
##  contrast   estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10   -0.152 0.0225 Inf    -0.196   -0.1077  -6.754  <.0001
##  diff_ESS1    -0.105 0.0304 Inf    -0.165   -0.0456  -3.463  0.0005
## 
## log_gdp.z.cm =  0:
##  contrast   estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10   -0.173 0.0156 Inf    -0.204   -0.1428 -11.099  <.0001
##  diff_ESS1    -0.150 0.0212 Inf    -0.192   -0.1086  -7.084  <.0001
## 
## log_gdp.z.cm =  1:
##  contrast   estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10   -0.195 0.0210 Inf    -0.236   -0.1537  -9.267  <.0001
##  diff_ESS1    -0.195 0.0290 Inf    -0.252   -0.1383  -6.728  <.0001
## 
## Degrees-of-freedom method: asymptotic 
## Confidence level used: 0.95
```

``` r
# Test for whether and how much the differences have changed
pairs(diff_mod6_log_GDP,infer=c(T,T))
```

```
## log_gdp.z.cm = -1:
##  contrast                estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10 - diff_ESS1 -0.046614 0.0293 Inf   -0.1041    0.0109  -1.589  0.1121
## 
## log_gdp.z.cm =  0:
##  contrast                estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10 - diff_ESS1 -0.023190 0.0184 Inf   -0.0592    0.0128  -1.262  0.2071
## 
## log_gdp.z.cm =  1:
##  contrast                estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10 - diff_ESS1  0.000235 0.0244 Inf   -0.0476    0.0481   0.010  0.9923
## 
## Degrees-of-freedom method: asymptotic 
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
##   [1] mnormt_2.1.1       Rdpack_2.6.4       gridExtra_2.3      writexl_1.5.4      readxl_1.4.5      
##   [6] rlang_1.1.6        magrittr_2.0.4     rockchalk_1.8.157  compiler_4.5.1     mgcv_1.9-3        
##  [11] png_0.1-8          vctrs_0.6.5        quadprog_1.5-8     crayon_1.5.3       pkgconfig_2.0.3   
##  [16] shape_1.4.6.1      fastmap_1.2.0      backports_1.5.0    labeling_0.4.3     pbivnorm_0.6.0    
##  [21] utf8_1.2.6         nloptr_2.2.1       purrr_1.1.0        xfun_0.53          glmnet_4.1-10     
##  [26] jomo_2.7-6         cachem_1.1.0       kutils_1.73        jsonlite_2.0.0     pan_1.9           
##  [31] jpeg_0.1-11        psych_2.5.6        lavaan_0.6-20      parallel_4.5.1     broom_1.0.10      
##  [36] cluster_2.1.8.1    R6_2.6.1           bslib_0.9.0        stringi_1.8.7      RColorBrewer_1.1-3
##  [41] car_3.1-3          boot_1.3-31        rpart_4.1.24       cellranger_1.1.0   jquerylib_0.1.4   
##  [46] estimability_1.5.1 Rcpp_1.1.0         iterators_1.0.14   base64enc_0.1-3    R.utils_2.13.0    
##  [51] splines_4.5.1      nnet_7.3-20        tidyselect_1.2.1   rstudioapi_0.17.1  abind_1.4-8       
##  [56] yaml_2.3.10        codetools_0.2-20   lattice_0.22-7     tibble_3.3.0       plyr_1.8.9        
##  [61] withr_3.0.2        S7_0.2.0           coda_0.19-4.1      evaluate_1.0.5     foreign_0.8-90    
##  [66] survival_3.8-3     zip_2.3.3          pillar_1.11.1      carData_3.0-5      mice_3.18.0       
##  [71] stats4_4.5.1       checkmate_2.3.3    foreach_1.5.2      reformulas_0.4.1   generics_0.1.4    
##  [76] grImport2_0.3-3    mathjaxr_1.8-0     scales_1.4.0       minqa_1.2.8        xtable_1.8-4      
##  [81] glue_1.8.0         tools_4.5.1        data.table_1.17.8  openxlsx_4.2.8     ggsignif_0.6.4    
##  [86] forcats_1.0.1      XML_3.99-0.19      mvtnorm_1.3-3      cowplot_1.2.0      grid_4.5.1        
##  [91] rbibutils_2.3      colorspace_2.1-2   htmlTable_2.4.3    Formula_1.2-5      cli_3.6.5         
##  [96] gtable_0.3.6       R.methodsS3_1.8.2  rstatix_0.7.3      sass_0.4.10        digest_0.6.37     
## [101] htmlwidgets_1.6.4  farver_2.1.2       htmltools_0.5.8.1  R.oo_1.27.1        lifecycle_1.0.4   
## [106] mitml_0.4-5        MASS_7.3-65
```

