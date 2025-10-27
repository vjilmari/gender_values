---
title: "Analysis for Examining the Gender Equality Paradox in Values Using benevolence Value"
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
cntry.ben<-diff_dat %>% group_by(cntry,essround) %>%
  summarise(ben.ctm=mean(ben,na.rm=T),
            ben.ctsd=sd(ben,na.rm=T)) %>%
  group_by(cntry) %>%
  summarise(ben.cm=mean(ben.ctm),
            ben.csd=mean(ben.ctsd)) 
```

```
## `summarise()` has grouped output by 'cntry'. You can override using the `.groups` argument.
```

``` r
grand_mean_ben<-mean(cntry.ben$ben.cm)
grand_sd_ben<-mean(cntry.ben$ben.csd)

# standardized
diff_dat$ben.z<-(diff_dat$ben-grand_mean_ben)/grand_sd_ben
hist(diff_dat$ben.z)
```

![](Analysis_benevolence_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

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

# value-based benevolence

cntry_ben_frame<-
  diff_dat %>% group_by(cntry,essround) %>%
  summarise('ben M' = weighted.mean(x=ben.z,w=pspwght),
            'ben SD' = sqrt(wtd.var(ben.z,pspwght))) %>%
  group_by(cntry) %>%
  summarise('ben M' = mean(x=`ben M`),
            'ben SD'= mean(x=`ben SD`))
```

```
## `summarise()` has grouped output by 'cntry'. You can override using the `.groups` argument.
```

``` r
cntry_ben_women_frame<-
  diff_dat %>%
  filter(gndr.c==-0.5) %>%
  group_by(cntry,essround) %>%
  summarise('ben M' = weighted.mean(x=ben.z,w=pspwght),
            'ben SD' = sqrt(wtd.var(ben.z,pspwght))) %>%
  group_by(cntry) %>%
  summarise('ben M Women' = mean(x=`ben M`),
            'ben SD Women'= mean(x=`ben SD`))
```

```
## `summarise()` has grouped output by 'cntry'. You can override using the `.groups` argument.
```

``` r
cntry_ben_men_frame<-
  diff_dat %>%
  filter(gndr.c==0.5) %>%
  group_by(cntry,essround) %>%
  summarise('ben M' = weighted.mean(x=ben.z,w=pspwght),
            'ben SD' = sqrt(wtd.var(ben.z,pspwght))) %>%
  group_by(cntry) %>%
  summarise('ben M Men' = mean(x=`ben M`),
            'ben SD Men'= mean(x=`ben SD`))
```

```
## `summarise()` has grouped output by 'cntry'. You can override using the `.groups` argument.
```

``` r
# link n and ben datasets

desc_frame<-
  left_join(
    x=cntry_n_frame,
    y=cntry_ben_frame,
    by="cntry"
  )

desc_frame<-
  left_join(
    x=desc_frame,
    y=cntry_ben_women_frame,
    by="cntry"
  )

desc_frame<-
  left_join(
    x=desc_frame,
    y=cntry_ben_men_frame,
    by="cntry"
  )

# Add country-specific differences
desc_frame$D<-desc_frame$`ben M Men`-desc_frame$`ben M Women`

desc_frame
```

```
## # A tibble: 33 × 10
##    cntry `n ESS rounds`     n `ben M` `ben SD` `ben M Women` `ben SD Women` `ben M Men` `ben SD Men`
##    <chr>          <dbl> <dbl>   <dbl>    <dbl>         <dbl>          <dbl>       <dbl>        <dbl>
##  1 AT                 6 13077  0.166     0.983        0.283           0.955      0.0405        0.997
##  2 BE                10 17313  0.218     0.826        0.334           0.810      0.0962        0.824
##  3 BG                 6 12641  0.0194    1.05         0.0597          1.04      -0.0239        1.05 
##  4 CH                10 16720  0.350     0.810        0.458           0.786      0.237         0.819
##  5 CY                 5  5105  0.411     0.813        0.478           0.785      0.340         0.833
##  6 CZ                 9 18934 -0.483     1.12        -0.340           1.10      -0.638         1.12 
##  7 DE                 9 25389  0.214     0.870        0.319           0.852      0.104         0.875
##  8 DK                 8 12198  0.359     0.859        0.497           0.811      0.217         0.884
##  9 EE                 9 16692 -0.208     0.952       -0.0722          0.924     -0.371         0.960
## 10 ES                 9 16954  0.378     0.878        0.449           0.855      0.304         0.895
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
    `ben M`, `ben SD`,
    `ben M Women`, `ben SD Women`,
    `ben M Men`, `ben SD Men`,
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
##    Country    `n ESS rounds`     n `ben M` `ben SD` `ben M Women` `ben SD Women` `ben M Men` `ben SD Men`
##    <chr>               <dbl> <dbl> <chr>   <chr>    <chr>         <chr>          <chr>       <chr>       
##  1 Austria                 6 13077 0.17    0.98     0.28          0.95           0.04        1.00        
##  2 Belgium                10 17313 0.22    0.83     0.33          0.81           0.10        0.82        
##  3 Bulgaria                6 12641 0.02    1.05     0.06          1.04           -0.02       1.05        
##  4 Switzerla…             10 16720 0.35    0.81     0.46          0.79           0.24        0.82        
##  5 Cyprus                  5  5105 0.41    0.81     0.48          0.78           0.34        0.83        
##  6 Czechia                 9 18934 -0.48   1.12     -0.34         1.10           -0.64       1.12        
##  7 Germany                 9 25389 0.21    0.87     0.32          0.85           0.10        0.87        
##  8 Denmark                 8 12198 0.36    0.86     0.50          0.81           0.22        0.88        
##  9 Estonia                 9 16692 -0.21   0.95     -0.07         0.92           -0.37       0.96        
## 10 Spain                   9 16954 0.38    0.88     0.45          0.86           0.30        0.89        
## 11 Finland                10 18050 0.07    0.94     0.25          0.90           -0.13       0.94        
## 12 France                 10 18720 -0.02   1.14     0.10          1.13           -0.17       1.14        
## 13 UK                     10 21456 0.15    0.96     0.27          0.92           0.01        0.97        
## 14 Greece                  5 12464 0.22    0.94     0.25          0.94           0.18        0.94        
## 15 Croatia                 4  6368 0.12    1.06     0.20          1.03           0.03        1.09        
## 16 Hungary                10 16006 -0.06   1.07     0.03          1.05           -0.16       1.09        
## 17 Ireland                10 20576 0.07    1.05     0.20          1.03           -0.05       1.05        
## 18 Israel                  6 13964 0.24    1.06     0.28          1.04           0.20        1.07        
## 19 Iceland                 5  3832 0.38    0.87     0.55          0.81           0.21        0.89        
## 20 Italy                   4  8663 -0.04   1.02     0.04          1.01           -0.12       1.02        
## 21 Lithuania               6 11714 -0.91   1.27     -0.81         1.27           -1.03       1.26        
## 22 Latvia                  2  2866 -0.21   1.06     -0.06         1.01           -0.40       1.09        
## 23 Montenegro              2  2441 0.05    1.15     0.13          1.14           -0.02       1.15        
## 24 Netherlan…             10 18048 -0.02   0.88     0.11          0.85           -0.16       0.88        
## 25 Norway                 10 15186 0.02    0.94     0.15          0.93           -0.12       0.93        
## 26 Poland                  9 15314 -0.07   0.93     0.00          0.93           -0.15       0.93        
## 27 Portugal               10 17705 -0.20   1.06     -0.16         1.05           -0.25       1.06        
## 28 Russia                  5 12139 -0.31   1.16     -0.27         1.17           -0.36       1.15        
## 29 Sweden                  9 14897 0.03    0.96     0.20          0.92           -0.16       0.96        
## 30 Slovenia               10 13238 0.04    0.90     0.14          0.89           -0.07       0.89        
## 31 Slovakia                7 11132 -0.36   1.01     -0.27         1.01           -0.46       1.01        
## 32 Turkey                  2  4108 0.14    0.98     0.13          0.97           0.14        0.98        
## 33 Ukraine                 5  9454 -0.51   1.27     -0.47         1.30           -0.56       1.23        
## # ℹ 5 more variables: D <chr>, GEI <chr>, GGGI <chr>, GDI <chr>, GDP <chr>
```

``` r
export(cntry_desc_tbl,"../results/ben/cntry_desc_tbl.xlsx",overwrite=T)
```

## Country-level correlations table


``` r
cor_frame<-
  desc_frame %>%
  select(
    VBMT=`ben M`,
    VBMT_Women=`ben M Women`,
    VBMT_Men=`ben M Men`,
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
  filename = "../results/ben/CorTable1.doc",
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
##   1. VBMT       0.01  0.29                                                                            
##                                                                                                       
##   2. VBMT_Women 0.11  0.30 .99                                                                        
##                            [.98, .99]                                                                 
##                                                                                                       
##   3. VBMT_Men   -0.10 0.30 .99          .95                                                           
##                            [.97, .99]   [.90, .98]                                                    
##                                                                                                       
##   4. D          -0.21 0.09 -.00         -.15         .16                                              
##                            [-.34, .34]  [-.47, .20]  [-.19, .48]                                      
##                                                                                                       
##   5. GEI        0.87  0.07 .25          .35          .15          -.60                                
##                            [-.11, .55]  [-.00, .62]  [-.21, .47]  [-.78, -.31]                        
##                                                                                                       
##   6. GGGI       0.73  0.05 .22          .34          .10          -.75         .73                    
##                            [-.13, .53]  [-.01, .61]  [-.25, .43]  [-.87, -.54] [.52, .86]             
##                                                                                                       
##   7. GDI        0.99  0.03 -.59         -.54         -.62         -.25         .07         .20        
##                            [-.77, -.31] [-.75, -.25] [-.80, -.36] [-.55, .10]  [-.29, .41] [-.16, .51]
##                                                                                                       
##   8. log_GDP    10.62 0.40 .44          .53          .35          -.58         .75         .67        
##                            [.12, .68]   [.23, .74]   [.00, .62]   [-.77, -.30] [.55, .87]  [.42, .82] 
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
mod0<-lmer(ben.z~(1|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod0)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ben.z ~ (1 | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1313458.4 1313491.4 -656726.2 1313452.4    441165 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.6758 -0.5293  0.0725  0.6676  5.0625 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  cntry    (Intercept) 0.08253  0.2873  
##  Residual             1.01287  1.0064  
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##              Estimate Std. Error        df t value Pr(>|t|)
## (Intercept)  0.002568   0.050042 33.004611   0.051    0.959
```

``` r
getVC(mod0)
```

```
##        grp        var1 var2 sdcor vcov
## 1    cntry (Intercept) <NA>  0.29 0.08
## 2 Residual        <NA> <NA>  1.01 1.01
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
## mean variation  0.07534513     NA       1
## sigma2          0.92465487      1      NA
## 
## $R2s
##          total within between
## f1  0.00000000      0      NA
## f2  0.00000000     NA       0
## v   0.00000000      0      NA
## m   0.07534513     NA       1
## f   0.00000000     NA      NA
## fv  0.00000000      0      NA
## fvm 0.07534513     NA      NA
```

## mod1: Gender fixed effect


``` r
mod1<-lmer(ben.z~gndr.c+(1|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod1)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ben.z ~ gndr.c + (1 | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1308548.4 1308592.4 -654270.2 1308540.4    441164 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.5167 -0.5246  0.0623  0.6618  5.3502 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  cntry    (Intercept) 0.08387  0.2896  
##  Residual             1.00165  1.0008  
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##               Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept) -1.465e-03  5.044e-02  3.300e+01  -0.029    0.977    
## gndr.c      -2.115e-01  3.009e-03  4.411e+05 -70.281   <2e-16 ***
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
## (Intercept) -0.001 0.050     33.003  -0.029 0.977 -0.104  0.101
## gndr.c      -0.211 0.003 441135.679 -70.281 0.000 -0.217 -0.206
```

``` r
getVC(mod1)
```

```
##        grp        var1 var2 sdcor vcov
## 1    cntry (Intercept) <NA>  0.29 0.08
## 2 Residual        <NA> <NA>  1.00 1.00
```

``` r
r2mlm(mod1,bargraph = F)
```

```
## $Decompositions
##                      total
## fixed           0.01013549
## slope variation 0.00000000
## mean variation  0.07647681
## sigma2          0.91338769
## 
## $R2s
##          total
## f   0.01013549
## v   0.00000000
## m   0.07647681
## fv  0.01013549
## fvm 0.08661231
```

## mod2: Gender fixed and random effect

* Include random effect correlation by default


``` r
mod2<-lmer(ben.z~gndr.c+(gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ben.z ~ gndr.c + (gndr.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1307919.9 1307985.9 -653953.9 1307907.9    441162 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.4654 -0.5175  0.0655  0.6605  5.3637 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr 
##  cntry    (Intercept) 0.083738 0.2894        
##           gndr.c      0.007587 0.0871   -0.01
##  Residual             0.999979 1.0000        
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##              Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept) -0.001205   0.050405 33.002224  -0.024    0.981    
## gndr.c      -0.202281   0.015564 32.284098 -12.996 2.25e-14 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##        (Intr)
## gndr.c -0.008
```

``` r
getFE(mod2,round=3)
```

```
##               Est.    SE     df       t     p     LL     UL
## (Intercept) -0.001 0.050 33.002  -0.024 0.981 -0.104  0.101
## gndr.c      -0.202 0.016 32.284 -12.996 0.000 -0.234 -0.171
```

``` r
getVC(mod2)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.29 0.08
## 2    cntry      gndr.c   <NA>  0.09 0.01
## 3    cntry (Intercept) gndr.c -0.01 0.00
## 4 Residual        <NA>   <NA>  1.00 1.00
```

``` r
r2mlm(mod2,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.009278680
## slope variation 0.001720511
## mean variation  0.076443390
## sigma2          0.912557419
## 
## $R2s
##           total
## f   0.009278680
## v   0.001720511
## m   0.076443390
## fv  0.010999191
## fvm 0.087442581
```

``` r
anova(mod1,mod2)
```

```
## Data: diff_dat
## Models:
## mod1: ben.z ~ gndr.c + (1 | cntry)
## mod2: ben.z ~ gndr.c + (gndr.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod1    4 1308548 1308592 -654270   1308540                         
## mod2    6 1307920 1307986 -653954   1307908 632.57  2  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

``` r
lvl2_var_cond_lvl1(mod2,lvl1.var = "gndr.c",lvl1.values = c(0.5,-0.5))
```

```
##   lvl1.value lvl2.cond.var lvl2.cond.sd
## 1        0.5    0.08540870    0.2922477
## 2       -0.5    0.08586003    0.2930188
```

* Test for random effect correlation


``` r
mod2_norecov<-lmer(ben.z~gndr.c+(gndr.c||cntry),data=diff_dat,REML=F,weights = pspwght,
                   control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2_norecov)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ben.z ~ gndr.c + (gndr.c || cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1307917.9 1307972.9 -653953.9 1307907.9    441163 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.4654 -0.5175  0.0654  0.6605  5.3639 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  cntry    (Intercept) 0.083738 0.28937 
##  cntry.1  gndr.c      0.007588 0.08711 
##  Residual             0.999979 0.99999 
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##              Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept) -0.001206   0.050405 33.002501  -0.024    0.981    
## gndr.c      -0.202278   0.015565 32.282881 -12.996 2.26e-14 ***
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
##               Est.    SE     df       t     p     LL     UL
## (Intercept) -0.001 0.050 33.003  -0.024 0.981 -0.104  0.101
## gndr.c      -0.202 0.016 32.283 -12.996 0.000 -0.234 -0.171
```

``` r
getVC(mod2_norecov)
```

```
##        grp        var1 var2 sdcor vcov
## 1    cntry (Intercept) <NA>  0.29 0.08
## 2  cntry.1      gndr.c <NA>  0.09 0.01
## 3 Residual        <NA> <NA>  1.00 1.00
```

``` r
anova(mod2_norecov,mod2)
```

```
## Data: diff_dat
## Models:
## mod2_norecov: ben.z ~ gndr.c + (gndr.c || cntry)
## mod2: ben.z ~ gndr.c + (gndr.c | cntry)
##              npar     AIC     BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)
## mod2_norecov    5 1307918 1307973 -653954   1307908                     
## mod2            6 1307920 1307986 -653954   1307908 0.0025  1     0.9602
```


## mod2 with Gender-equality index (GEI)


``` r
mod2_GEI<-lmer(ben.z~gndr.c+gei.z.cm+gndr.c:gei.z.cm+
                 (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2_GEI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ben.z ~ gndr.c + gei.z.cm + gndr.c:gei.z.cm + (gndr.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1272087.7 1272175.5 -636035.8 1272071.7    431770 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.5290 -0.5160  0.0681  0.6622  5.4052 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr
##  cntry    (Intercept) 0.073106 0.27038      
##           gndr.c      0.004673 0.06836  0.29
##  Residual             0.986230 0.99309      
## Number of obs: 431778, groups:  cntry, 32
## 
## Fixed effects:
##                 Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)      0.01500    0.04783 31.99860   0.314 0.755920    
## gndr.c          -0.20433    0.01259 31.41523 -16.235  < 2e-16 ***
## gei.z.cm         0.06859    0.04861 32.03628   1.411 0.167851    
## gndr.c:gei.z.cm -0.05453    0.01300 33.46379  -4.195 0.000189 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c g.z.cm
## gndr.c       0.275              
## gei.z.cm     0.000  0.000       
## gndr.c:g.z.  0.000 -0.024  0.270
```

``` r
getFE(mod2_GEI,round=3)
```

```
##                   Est.    SE     df       t     p     LL     UL
## (Intercept)      0.015 0.048 31.999   0.314 0.756 -0.082  0.112
## gndr.c          -0.204 0.013 31.415 -16.235 0.000 -0.230 -0.179
## gei.z.cm         0.069 0.049 32.036   1.411 0.168 -0.030  0.168
## gndr.c:gei.z.cm -0.055 0.013 33.464  -4.195 0.000 -0.081 -0.028
```

``` r
getVC(mod2_GEI)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.27 0.07
## 2    cntry      gndr.c   <NA>  0.07 0.00
## 3    cntry (Intercept) gndr.c  0.29 0.01
## 4 Residual        <NA>   <NA>  0.99 0.99
```

``` r
r2mlm(mod2_GEI,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.013870400
## slope variation 0.001080745
## mean variation  0.067647200
## sigma2          0.917401655
## 
## $R2s
##           total
## f   0.013870400
## v   0.001080745
## m   0.067647200
## fv  0.014951145
## fvm 0.082598345
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
## Time difference of 40.50753 secs
```

``` r
round(ddsc_mod2_GEI$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.086        0.293            1     1.086 0.079   7204.152 0.998   0.998
## 2        0.5         0.085        0.292            1     1.085 0.079   6164.576 0.998   0.998
```

``` r
round(ddsc_mod2_GEI$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gei.z.cm
## means_y1          -0.088 0.290    1.000           1.000    0.950           0.950    0.143
## means_y1_scaled   -0.305 1.010    1.000           1.000    0.950           0.950    0.143
## means_y2           0.118 0.284    0.950           0.950    1.000           1.000    0.337
## means_y2_scaled    0.410 0.989    0.950           0.950    1.000           1.000    0.337
## gei.z.cm           0.000 1.000    0.143           0.143    0.337           0.337    1.000
## gei.z.cm_scaled    0.000 1.000    0.143           0.143    0.337           0.337    1.000
## diff_score        -0.205 0.091    0.223           0.223   -0.094          -0.094   -0.593
## diff_score_scaled -0.716 0.318    0.223           0.223   -0.094          -0.094   -0.593
##                   gei.z.cm_scaled diff_score diff_score_scaled
## means_y1                    0.143      0.223             0.223
## means_y1_scaled             0.143      0.223             0.223
## means_y2                    0.337     -0.094            -0.094
## means_y2_scaled             0.337     -0.094            -0.094
## gei.z.cm                    1.000     -0.593            -0.593
## gei.z.cm_scaled             1.000     -0.593            -0.593
## diff_score                 -0.593      1.000             1.000
## diff_score_scaled          -0.593      1.000             1.000
```

``` r
round(ddsc_mod2_GEI$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                       0.598 0.143 33.464   4.195   0.000    0.308    0.888
## w_11                          0.096 0.047 32.074   2.028   0.051    0.000    0.192
## w_21                          0.041 0.051 32.062   0.814   0.421   -0.062    0.145
## r_xy1                         0.331 0.163 32.074   2.028   0.051   -0.001    0.663
## r_xy2                         0.146 0.179 32.062   0.814   0.421   -0.218    0.510
## b_11                          0.334 0.165 32.074   2.028   0.051   -0.001    0.669
## b_21                          0.144 0.177 32.062   0.814   0.421   -0.216    0.504
## main_effect                   0.069 0.049 32.036   1.411   0.168   -0.030    0.168
## moderator_effect             -0.204 0.013 31.415 -16.235   0.000   -0.230   -0.179
## interaction                  -0.055 0.013 33.464  -4.195   0.000   -0.081   -0.028
## q_b11_b21                     0.202    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.197    NA     NA      NA      NA       NA       NA
## cross_over_point             -3.747    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.014 0.054 32.136  -0.262   0.795   -0.123    0.095
## interaction_vs_main_bscale   -0.049 0.187 32.136  -0.262   0.795   -0.429    0.331
## interaction_vs_main_rscale   -0.053 0.190 32.134  -0.279   0.782   -0.441    0.335
## dadas                        -0.083 0.102 32.062  -0.814   0.789   -0.289    0.124
## dadas_bscale                 -0.288 0.354 32.062  -0.814   0.789   -1.008    0.432
## dadas_rscale                 -0.291 0.357 32.062  -0.814   0.789   -1.019    0.437
## abs_diff                      0.055 0.013 33.464   4.195   0.000    0.028    0.081
## abs_sum                       0.137 0.097 32.036   1.411   0.084   -0.061    0.335
## abs_diff_bscale               0.190 0.045 33.464   4.195   0.000    0.098    0.282
## abs_sum_bscale                0.478 0.339 32.036   1.411   0.084   -0.212    1.168
## abs_diff_rscale               0.185 0.046 33.435   3.989   0.000    0.091    0.279
## abs_sum_rscale                0.476 0.339 32.036   1.405   0.085   -0.214    1.166
```

``` r
round(ddsc_mod2_GEI$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
##  0.000 -0.009  0.002  1.000  0.960
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
## r_xy1_y2                         0.593 0.142  4.161  0.000    0.313    0.872
## r_xy1                            0.337 0.166  2.022  0.043    0.010    0.663
## r_xy2                            0.143 0.175  0.820  0.413   -0.200    0.486
## b_11                             0.333 0.165  2.022  0.043    0.010    0.656
## b_21                             0.145 0.177  0.820  0.413   -0.202    0.491
## b_10                             0.410 0.162  2.530  0.011    0.092    0.728
## b_20                            -0.305 0.174 -1.755  0.079   -0.646    0.036
## res_cov_y1_y2                    0.873 0.222  3.933  0.000    0.438    1.308
## diff_b10_b20                     0.716 0.045 16.077  0.000    0.628    0.803
## diff_b11_b21                     0.188 0.045  4.161  0.000    0.100    0.277
## diff_rxy1_rxy2                   0.193 0.044  4.352  0.000    0.106    0.280
## q_b11_b21                        0.200 0.047  4.268  0.000    0.108    0.292
## q_rxy1_rxy2                      0.206 0.048  4.321  0.000    0.113    0.299
## cross_over_point                -3.803 0.944 -4.029  0.000   -5.653   -1.953
## sum_b11_b21                      0.478 0.339  1.411  0.158   -0.186    1.142
## main_effect                      0.239 0.169  1.411  0.158   -0.093    0.571
## interaction_vs_main_effect      -0.051 0.187 -0.272  0.786   -0.417    0.315
## diff_abs_b11_abs_b21             0.188 0.045  4.161  0.000    0.100    0.277
## abs_diff_b11_b21                 0.188 0.045  4.161  0.000    0.100    0.277
## abs_sum_b11_b21                  0.478 0.339  1.411  0.079   -0.186    1.142
## dadas                           -0.290 0.354 -0.820  0.794   -0.983    0.403
## q_r_equivalence                  0.106 0.048  2.223  0.987       NA       NA
## q_b_equivalence                  0.100 0.047  2.138  0.984       NA       NA
## cross_over_point_equivalence     3.803 0.944  4.029  1.000       NA       NA
## cross_over_point_minimal_effect  3.803 0.944  4.029  0.000       NA       NA
```

``` r
round(ddsc_sem_GEI$variance_test,3)
```

```
##              est    se      z pvalue ci.lower ci.upper
## cov_y1y2   0.920 0.236  3.896  0.000    0.457    1.383
## var_y1     0.948 0.237  4.000  0.000    0.484    1.413
## var_y2     0.989 0.247  4.000  0.000    0.504    1.474
## var_diff  -0.041 0.108 -0.378  0.705   -0.252    0.170
## var_ratio  0.959 0.106  9.039  0.000    0.751    1.167
## cor_y1y2   0.950 0.017 54.868  0.000    0.916    0.984
```

### Figure 


``` r
# plot the results

# start with obtaining predicted values for means and differences

# refit reduced and full models with GGGI in original scale


mod2_GEI_unstd<-lmer(ben.z~gndr.c+gei.cm+gndr.c:gei.cm+
                 (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

mod2_GEI_unstd_red<-lmer(ben.z~gndr.c+
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


p1.ben.flags<-
  ggplot(p,aes(y=yvar,x=gei.cm,color=gndr.c))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("Gender Equality Index (Average 2002-2022)")+
  #ylim(c(min.y.pred,max.y.pred))+
  #xlim(c(0.60,1.00))+
  ylab("Value benevolence (Average 2002-2022)")+
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

p2.ben.flags<-ggplot(p2,aes(y=yvar,x=gei.cm))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("Gender Equality Index (Average 2002-2022)")+
  #xlim(c(0.60,1.00))+
  #ylim(c(min.y.narrow,max.y.narrow))+
  ylab("Difference in value benevolence")+
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
  ggarrange(p1.ben.flags,p2.ben.flags,align = "v",
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

![](Analysis_benevolence_files/figure-html/unnamed-chunk-18-1.png)<!-- -->

``` r
png(filename = 
      "../results/ben/GEI_flags.png",
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
mod2_GGGI<-lmer(ben.z~gndr.c+gggi.z.cm+gndr.c:gggi.z.cm+
                 (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2_GGGI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ben.z ~ gndr.c + gggi.z.cm + gndr.c:gggi.z.cm + (gndr.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##  932690.1  932775.3 -466337.0  932674.1    314638 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.5346 -0.5128  0.0629  0.6661  5.3427 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr
##  cntry    (Intercept) 0.078657 0.28046      
##           gndr.c      0.003028 0.05503  0.29
##  Residual             0.997916 0.99896      
## Number of obs: 314646, groups:  cntry, 33
## 
## Fixed effects:
##                  Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)       0.01656    0.04887 33.01508   0.339   0.7369    
## gndr.c           -0.19869    0.01042 33.20533 -19.072  < 2e-16 ***
## gggi.z.cm         0.08539    0.04965 33.07033   1.720   0.0948 .  
## gndr.c:gggi.z.cm -0.06202    0.01092 37.12651  -5.678  1.7e-06 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c ggg.z.
## gndr.c       0.268              
## gggi.z.cm    0.000  0.000       
## gndr.c:gg..  0.000 -0.020  0.259
```

``` r
getFE(mod2_GGGI,round=3)
```

```
##                    Est.    SE     df       t     p     LL     UL
## (Intercept)       0.017 0.049 33.015   0.339 0.737 -0.083  0.116
## gndr.c           -0.199 0.010 33.205 -19.072 0.000 -0.220 -0.178
## gggi.z.cm         0.085 0.050 33.070   1.720 0.095 -0.016  0.186
## gndr.c:gggi.z.cm -0.062 0.011 37.127  -5.678 0.000 -0.084 -0.040
```

``` r
getVC(mod2_GGGI)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.28 0.08
## 2    cntry      gndr.c   <NA>  0.06 0.00
## 3    cntry (Intercept) gndr.c  0.29 0.00
## 4 Residual        <NA>   <NA>  1.00 1.00
```

``` r
r2mlm(mod2_GGGI,bargraph = F)
```

```
## $Decompositions
##                        total
## fixed           0.0148925892
## slope variation 0.0006880666
## mean variation  0.0716232303
## sigma2          0.9127961138
## 
## $R2s
##            total
## f   0.0148925892
## v   0.0006880666
## m   0.0716232303
## fv  0.0155806558
## fvm 0.0872038862
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
## Time difference of 53.1099 secs
```

``` r
round(ddsc_mod2_GGGI$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.086        0.293            1     1.086 0.079   7204.152 0.998   0.998
## 2        0.5         0.085        0.292            1     1.085 0.079   6164.576 0.998   0.998
```

``` r
round(ddsc_mod2_GGGI$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gggi.z.cm
## means_y1          -0.082 0.300    1.000           1.000    0.957           0.957     0.178
## means_y1_scaled   -0.274 0.996    1.000           1.000    0.957           0.957     0.178
## means_y2           0.116 0.302    0.957           0.957    1.000           1.000     0.388
## means_y2_scaled    0.384 1.003    0.957           0.957    1.000           1.000     0.388
## gggi.z.cm          0.000 1.000    0.178           0.178    0.388           0.388     1.000
## gggi.z.cm_scaled   0.000 1.000    0.178           0.178    0.388           0.388     1.000
## diff_score        -0.198 0.088    0.123           0.123   -0.170          -0.170    -0.724
## diff_score_scaled -0.658 0.293    0.123           0.123   -0.170          -0.170    -0.724
##                   gggi.z.cm_scaled diff_score diff_score_scaled
## means_y1                     0.178      0.123             0.123
## means_y1_scaled              0.178      0.123             0.123
## means_y2                     0.388     -0.170            -0.170
## means_y2_scaled              0.388     -0.170            -0.170
## gggi.z.cm                    1.000     -0.724            -0.724
## gggi.z.cm_scaled             1.000     -0.724            -0.724
## diff_score                  -0.724      1.000             1.000
## diff_score_scaled           -0.724      1.000             1.000
```

``` r
round(ddsc_mod2_GGGI$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                       0.704 0.124 37.127   5.678   0.000    0.453    0.955
## w_11                          0.116 0.049 33.144   2.399   0.022    0.018    0.215
## w_21                          0.054 0.051 33.110   1.059   0.297   -0.050    0.159
## r_xy1                         0.388 0.162 33.144   2.399   0.022    0.059    0.718
## r_xy2                         0.180 0.170 33.110   1.059   0.297   -0.166    0.526
## b_11                          0.387 0.161 33.144   2.399   0.022    0.059    0.715
## b_21                          0.181 0.171 33.110   1.059   0.297   -0.166    0.528
## main_effect                   0.085 0.050 33.070   1.720   0.095   -0.016    0.186
## moderator_effect             -0.199 0.010 33.205 -19.072   0.000   -0.220   -0.178
## interaction                  -0.062 0.011 37.127  -5.678   0.000   -0.084   -0.040
## q_b11_b21                     0.226    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.228    NA     NA      NA      NA       NA       NA
## cross_over_point             -3.204    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.023 0.054 33.239  -0.436   0.665   -0.132    0.086
## interaction_vs_main_bscale   -0.078 0.178 33.239  -0.436   0.665   -0.440    0.284
## interaction_vs_main_rscale   -0.076 0.177 33.241  -0.430   0.670   -0.436    0.284
## dadas                        -0.109 0.103 33.110  -1.059   0.851   -0.318    0.100
## dadas_bscale                 -0.362 0.341 33.110  -1.059   0.851   -1.056    0.333
## dadas_rscale                 -0.360 0.340 33.110  -1.059   0.851   -1.052    0.332
## abs_diff                      0.062 0.011 37.127   5.678   0.000    0.040    0.084
## abs_sum                       0.171 0.099 33.070   1.720   0.047   -0.031    0.373
## abs_diff_bscale               0.206 0.036 37.127   5.678   0.000    0.133    0.280
## abs_sum_bscale                0.568 0.330 33.070   1.720   0.047   -0.104    1.240
## abs_diff_rscale               0.208 0.036 37.194   5.777   0.000    0.135    0.281
## abs_sum_rscale                0.569 0.330 33.070   1.722   0.047   -0.103    1.240
```

``` r
round(ddsc_mod2_GGGI$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
##  0.000 -0.009  0.002  1.000  0.960
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
## r_xy1_y2                         0.724 0.120  6.034  0.000    0.489    0.960
## r_xy1                            0.388 0.160  2.421  0.015    0.074    0.703
## r_xy2                            0.178 0.171  1.040  0.298   -0.158    0.514
## b_11                             0.390 0.161  2.421  0.015    0.074    0.705
## b_21                             0.178 0.171  1.040  0.298   -0.157    0.512
## b_10                             0.384 0.159  2.424  0.015    0.074    0.695
## b_20                            -0.274 0.168 -1.631  0.103   -0.604    0.055
## res_cov_y1_y2                    0.861 0.214  4.019  0.000    0.441    1.281
## diff_b10_b20                     0.658 0.035 19.014  0.000    0.591    0.726
## diff_b11_b21                     0.212 0.035  6.034  0.000    0.143    0.281
## diff_rxy1_rxy2                   0.210 0.036  5.920  0.000    0.141    0.280
## q_b11_b21                        0.232 0.040  5.848  0.000    0.154    0.310
## q_rxy1_rxy2                      0.230 0.039  5.868  0.000    0.153    0.307
## cross_over_point                -3.103 0.540 -5.751  0.000   -4.161   -2.046
## sum_b11_b21                      0.567 0.330  1.720  0.086   -0.079    1.214
## main_effect                      0.284 0.165  1.720  0.086   -0.040    0.607
## interaction_vs_main_effect      -0.071 0.178 -0.402  0.688   -0.420    0.277
## diff_abs_b11_abs_b21             0.212 0.035  6.034  0.000    0.143    0.281
## abs_diff_b11_b21                 0.212 0.035  6.034  0.000    0.143    0.281
## abs_sum_b11_b21                  0.567 0.330  1.720  0.043   -0.079    1.214
## dadas                           -0.355 0.341 -1.040  0.851   -1.024    0.314
## q_r_equivalence                  0.130 0.039  3.315  1.000       NA       NA
## q_b_equivalence                  0.132 0.040  3.328  1.000       NA       NA
## cross_over_point_equivalence     3.103 0.540  5.751  1.000       NA       NA
## cross_over_point_minimal_effect  3.103 0.540  5.751  0.000       NA       NA
```

``` r
round(ddsc_sem_GGGI$variance_test,3)
```

```
##             est    se      z pvalue ci.lower ci.upper
## cov_y1y2  0.928 0.234  3.972   0.00    0.470    1.386
## var_y1    0.976 0.240  4.062   0.00    0.505    1.448
## var_y2    0.963 0.237  4.062   0.00    0.498    1.428
## var_diff  0.014 0.098  0.139   0.89   -0.178    0.205
## var_ratio 1.014 0.102  9.913   0.00    0.814    1.215
## cor_y1y2  0.957 0.015 65.496   0.00    0.928    0.986
```

### Figure


``` r
# plot the results

# start with obtaining predicted values for means and differences

# refit reduced and full models with GGGI in original scale


mod2_GGGI_unstd<-lmer(ben.z~gndr.c+gggi.cm+gndr.c:gggi.cm+
                       (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
                     control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

mod2_GGGI_unstd_red<-lmer(ben.z~gndr.c+
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


p1.ben.flags<-
  ggplot(p,aes(y=yvar,x=gggi.cm,color=gndr.c))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("Global Gender Gap Index (Average 2002-2022)")+
  #ylim(c(min.y.pred,max.y.pred))+
  #xlim(c(0.60,1.00))+
  ylab("Value benevolence (Average 2002-2022)")+
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

p2.ben.flags<-ggplot(p2,aes(y=yvar,x=gggi.cm))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("Global Gender Gap Index (Average 2002-2022)")+
  #xlim(c(0.60,1.00))+
  #ylim(c(min.y.narrow,max.y.narrow))+
  ylab("Difference in value benevolence")+
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
  ggarrange(p1.ben.flags,p2.ben.flags,align = "v",
            ncol=1,nrow=2,heights=c(2,1))
```

```
## Warning in geom_text(inherit.aes = F, aes(x = 0.61, y = -0.45, label = coef_q, : All aesthetics have length 1, but the data has 502 rows.
## ℹ Please consider using `annotate()` or provide this layer with data containing a single row.
```

``` r
pflag_comb
```

![](Analysis_benevolence_files/figure-html/unnamed-chunk-21-1.png)<!-- -->

``` r
png(filename = 
      "../results/ben/GGGI_flags_new.png",
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
mod2_GDI<-lmer(ben.z~gndr.c+gdi.z.cm+gndr.c:gdi.z.cm+
                 (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2_GDI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ben.z ~ gndr.c + gdi.z.cm + gndr.c:gdi.z.cm + (gndr.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1307907.4 1307995.4 -653945.7 1307891.4    441160 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.4672 -0.5175  0.0655  0.6605  5.3633 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr 
##  cntry    (Intercept) 0.055738 0.2361        
##           gndr.c      0.007242 0.0851   -0.19
##  Residual             0.999978 1.0000        
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##                  Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)     -0.001209   0.041137 32.995118  -0.029 0.976722    
## gndr.c          -0.202232   0.015224 32.832423 -13.284 9.31e-15 ***
## gdi.z.cm        -0.170017   0.041797 33.067960  -4.068 0.000277 ***
## gndr.c:gdi.z.cm -0.022226   0.015698 34.827370  -1.416 0.165707    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c gd.z.c
## gndr.c      -0.184              
## gdi.z.cm     0.000  0.000       
## gndr.c:gd..  0.000 -0.005 -0.181
```

``` r
getFE(mod2_GDI,round=3)
```

```
##                   Est.    SE     df       t     p     LL     UL
## (Intercept)     -0.001 0.041 32.995  -0.029 0.977 -0.085  0.082
## gndr.c          -0.202 0.015 32.832 -13.284 0.000 -0.233 -0.171
## gdi.z.cm        -0.170 0.042 33.068  -4.068 0.000 -0.255 -0.085
## gndr.c:gdi.z.cm -0.022 0.016 34.827  -1.416 0.166 -0.054  0.010
```

``` r
getVC(mod2_GDI)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.24 0.06
## 2    cntry      gndr.c   <NA>  0.09 0.01
## 3    cntry (Intercept) gndr.c -0.19 0.00
## 4 Residual        <NA>   <NA>  1.00 1.00
```

``` r
r2mlm(mod2_GDI,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.026626260
## slope variation 0.001655848
## mean variation  0.051571019
## sigma2          0.920146873
## 
## $R2s
##           total
## f   0.026626260
## v   0.001655848
## m   0.051571019
## fv  0.028282108
## fvm 0.079853127
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
## Time difference of 50.56133 secs
```

``` r
round(ddsc_mod2_GDI$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.086        0.293            1     1.086 0.079   7204.152 0.998   0.998
## 2        0.5         0.085        0.292            1     1.085 0.079   6164.576 0.998   0.998
```

``` r
round(ddsc_mod2_GDI$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gdi.z.cm
## means_y1          -0.102 0.297    1.000           1.000    0.952           0.952   -0.612
## means_y1_scaled   -0.343 0.999    1.000           1.000    0.952           0.952   -0.612
## means_y2           0.100 0.298    0.952           0.952    1.000           1.000   -0.531
## means_y2_scaled    0.335 1.001    0.952           0.952    1.000           1.000   -0.531
## gdi.z.cm           0.000 1.000   -0.612          -0.612   -0.531          -0.531    1.000
## gdi.z.cm_scaled    0.000 1.000   -0.612          -0.612   -0.531          -0.531    1.000
## diff_score        -0.202 0.092    0.149           0.149   -0.160          -0.160   -0.259
## diff_score_scaled -0.679 0.309    0.149           0.149   -0.160          -0.160   -0.259
##                   gdi.z.cm_scaled diff_score diff_score_scaled
## means_y1                   -0.612      0.149             0.149
## means_y1_scaled            -0.612      0.149             0.149
## means_y2                   -0.531     -0.160            -0.160
## means_y2_scaled            -0.531     -0.160            -0.160
## gdi.z.cm                    1.000     -0.259            -0.259
## gdi.z.cm_scaled             1.000     -0.259            -0.259
## diff_score                 -0.259      1.000             1.000
## diff_score_scaled          -0.259      1.000             1.000
```

``` r
round(ddsc_mod2_GDI$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                       0.242 0.171 34.827   1.416   0.166   -0.105    0.589
## w_11                         -0.159 0.044 33.096  -3.620   0.001   -0.248   -0.070
## w_21                         -0.181 0.041 33.163  -4.406   0.000   -0.265   -0.098
## r_xy1                        -0.534 0.148 33.096  -3.620   0.001   -0.835   -0.234
## r_xy2                        -0.608 0.138 33.163  -4.406   0.000   -0.889   -0.327
## b_11                         -0.534 0.148 33.096  -3.620   0.001   -0.834   -0.234
## b_21                         -0.609 0.138 33.163  -4.406   0.000   -0.890   -0.328
## main_effect                  -0.170 0.042 33.068  -4.068   0.000   -0.255   -0.085
## moderator_effect             -0.202 0.015 32.832 -13.284   0.000   -0.233   -0.171
## interaction                  -0.022 0.016 34.827  -1.416   0.166   -0.054    0.010
## q_b11_b21                     0.111    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.110    NA     NA      NA      NA       NA       NA
## cross_over_point             -9.099    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.148 0.047 33.178  -3.129   0.004   -0.244   -0.052
## interaction_vs_main_bscale   -0.497 0.159 33.178  -3.129   0.004   -0.819   -0.174
## interaction_vs_main_rscale   -0.498 0.159 33.178  -3.130   0.004   -0.821   -0.174
## dadas                        -0.318 0.088 33.096  -3.620   1.000   -0.496   -0.139
## dadas_bscale                 -1.068 0.295 33.096  -3.620   1.000   -1.668   -0.468
## dadas_rscale                 -1.069 0.295 33.096  -3.620   1.000   -1.670   -0.468
## abs_diff                      0.022 0.016 34.827   1.416   0.083   -0.010    0.054
## abs_sum                       0.340 0.084 33.068   4.068   0.000    0.170    0.510
## abs_diff_bscale               0.075 0.053 34.827   1.416   0.083   -0.032    0.182
## abs_sum_bscale                1.143 0.281 33.068   4.068   0.000    0.571    1.714
## abs_diff_rscale               0.074 0.053 34.816   1.395   0.086   -0.034    0.181
## abs_sum_rscale                1.142 0.281 33.068   4.067   0.000    0.571    1.714
```

``` r
round(ddsc_mod2_GDI$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
##  0.000 -0.009  0.002  1.000  0.960
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
## r_xy1_y2                         0.259 0.168  1.541  0.123   -0.070    0.589
## r_xy1                           -0.531 0.148 -3.595  0.000   -0.820   -0.241
## r_xy2                           -0.612 0.138 -4.440  0.000   -0.881   -0.342
## b_11                            -0.531 0.148 -3.595  0.000   -0.820   -0.242
## b_21                            -0.611 0.138 -4.440  0.000   -0.881   -0.341
## b_10                             0.335 0.145  2.305  0.021    0.050    0.620
## b_20                            -0.343 0.136 -2.533  0.011   -0.609   -0.078
## res_cov_y1_y2                    0.609 0.155  3.926  0.000    0.305    0.913
## diff_b10_b20                     0.679 0.051 13.279  0.000    0.578    0.779
## diff_b11_b21                     0.080 0.052  1.541  0.123   -0.022    0.182
## diff_rxy1_rxy2                   0.081 0.052  1.562  0.118   -0.021    0.183
## q_b11_b21                        0.119 0.077  1.542  0.123   -0.032    0.270
## q_rxy1_rxy2                      0.120 0.077  1.559  0.119   -0.031    0.272
## cross_over_point                -8.487 5.545 -1.531  0.126  -19.354    2.381
## sum_b11_b21                     -1.142 0.281 -4.068  0.000   -1.692   -0.592
## main_effect                     -0.571 0.140 -4.068  0.000   -0.846   -0.296
## interaction_vs_main_effect      -0.491 0.159 -3.089  0.002   -0.803   -0.179
## diff_abs_b11_abs_b21            -0.080 0.052 -1.541  0.123   -0.182    0.022
## abs_diff_b11_b21                 0.080 0.052  1.541  0.062   -0.022    0.182
## abs_sum_b11_b21                  1.142 0.281  4.068  0.000    0.592    1.692
## dadas                           -1.062 0.295 -3.595  1.000   -1.641   -0.483
## q_r_equivalence                  0.020 0.077  0.265  0.604       NA       NA
## q_b_equivalence                  0.019 0.077  0.245  0.597       NA       NA
## cross_over_point_equivalence     8.487 5.545  1.531  0.937       NA       NA
## cross_over_point_minimal_effect  8.487 5.545  1.531  0.063       NA       NA
```

``` r
round(ddsc_sem_GDI$variance_test,3)
```

```
##             est    se      z pvalue ci.lower ci.upper
## cov_y1y2  0.924 0.233  3.962  0.000    0.467    1.380
## var_y1    0.971 0.239  4.062  0.000    0.503    1.440
## var_y2    0.968 0.238  4.062  0.000    0.501    1.435
## var_diff  0.004 0.103  0.034  0.973   -0.198    0.205
## var_ratio 1.004 0.107  9.419  0.000    0.795    1.212
## cor_y1y2  0.952 0.016 58.835  0.000    0.921    0.984
```

### Figure


``` r
# plot the results

# start with obtaining predicted values for means and differences


mod2_GDI_unstd<-lmer(ben.z~gndr.c+gdi.cm+gndr.c:gdi.cm+
                       (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
                     control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

mod2_GDI_unstd_red<-lmer(ben.z~gndr.c+
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


p1.ben.flags<-
  ggplot(p,aes(y=yvar,x=gdi.cm,color=gndr.c))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("Gender Development Index (Average 2002-2022)")+
  #ylim(c(min.y.pred,max.y.pred))+
  #xlim(c(0.60,1.00))+
  ylab("Value benevolence (Average 2002-2022)")+
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

#p1.ben.flags


p2.ben.flags<-ggplot(p2,aes(y=yvar,x=gdi.cm))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("Gender Development Index (Average 2002-2022)")+
  #xlim(c(0.60,1.00))+
  #ylim(c(min.y.narrow,max.y.narrow))+
  ylab("Difference in value benevolence")+
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

#p2.ben.flags


pflag_comb<-
  ggarrange(p1.ben.flags,p2.ben.flags,align = "v",
            ncol=1,nrow=2,heights=c(2,1))
```

```
## Warning in geom_text(inherit.aes = F, aes(x = 0.9, y = -0.45, label = coef_q, : All aesthetics have length 1, but the data has 282 rows.
## ℹ Please consider using `annotate()` or provide this layer with data containing a single row.
```

``` r
pflag_comb
```

![](Analysis_benevolence_files/figure-html/unnamed-chunk-24-1.png)<!-- -->

``` r
png(filename = 
      "../results/ben/GDI_flags.png",
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
mod2_log_GDP<-lmer(ben.z~gndr.c+log_gdp.z.cm+
                     gndr.c:log_gdp.z.cm+
                 (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2_log_GDP)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ben.z ~ gndr.c + log_gdp.z.cm + gndr.c:log_gdp.z.cm + (gndr.c |      cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1307898.6 1307986.6 -653941.3 1307882.6    441160 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.4635 -0.5176  0.0654  0.6606  5.3669 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr
##  cntry    (Intercept) 0.067493 0.2598       
##           gndr.c      0.004762 0.0690   0.37
##  Residual             0.999981 1.0000       
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##                      Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)          0.001557   0.045270 32.990200   0.034 0.972775    
## gndr.c              -0.201667   0.012505 31.719025 -16.127  < 2e-16 ***
## log_gdp.z.cm         0.127859   0.045436 33.019170   2.814 0.008182 ** 
## gndr.c:log_gdp.z.cm -0.053290   0.012672 32.834415  -4.205 0.000189 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c lg_g..
## gndr.c       0.352              
## lg_gdp.z.cm  0.022  0.008       
## gndr.c:l_..  0.007 -0.004  0.349
```

``` r
getFE(mod2_log_GDP,round=3)
```

```
##                       Est.    SE     df       t     p     LL     UL
## (Intercept)          0.002 0.045 32.990   0.034 0.973 -0.091  0.094
## gndr.c              -0.202 0.013 31.719 -16.127 0.000 -0.227 -0.176
## log_gdp.z.cm         0.128 0.045 33.019   2.814 0.008  0.035  0.220
## gndr.c:log_gdp.z.cm -0.053 0.013 32.834  -4.205 0.000 -0.079 -0.028
```

``` r
getVC(mod2_log_GDP)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.26 0.07
## 2    cntry      gndr.c   <NA>  0.07 0.00
## 3    cntry (Intercept) gndr.c  0.37 0.01
## 4 Residual        <NA>   <NA>  1.00 1.00
```

``` r
r2mlm(mod2_log_GDP,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.022633365
## slope variation 0.001082653
## mean variation  0.061296200
## sigma2          0.914987782
## 
## $R2s
##           total
## f   0.022633365
## v   0.001082653
## m   0.061296200
## fv  0.023716018
## fvm 0.085012218
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
## Time difference of 31.94234 secs
```

``` r
round(ddsc_mod2_log_GDP$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.086        0.293            1     1.086 0.079   7204.152 0.998   0.998
## 2        0.5         0.085        0.292            1     1.085 0.079   6164.576 0.998   0.998
```

``` r
round(ddsc_mod2_log_GDP$descriptives,3)
```

```
##                          M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled log_gdp.z.cm
## means_y1            -0.102 0.297    1.000           1.000    0.952           0.952        0.346
## means_y1_scaled     -0.343 0.999    1.000           1.000    0.952           0.952        0.346
## means_y2             0.100 0.298    0.952           0.952    1.000           1.000        0.523
## means_y2_scaled      0.335 1.001    0.952           0.952    1.000           1.000        0.523
## log_gdp.z.cm        -0.022 1.012    0.346           0.346    0.523           0.523        1.000
## log_gdp.z.cm_scaled  0.000 1.000    0.346           0.346    0.523           0.523        1.000
## diff_score          -0.202 0.092    0.149           0.149   -0.160          -0.160       -0.577
## diff_score_scaled   -0.679 0.309    0.149           0.149   -0.160          -0.160       -0.577
##                     log_gdp.z.cm_scaled diff_score diff_score_scaled
## means_y1                          0.346      0.149             0.149
## means_y1_scaled                   0.346      0.149             0.149
## means_y2                          0.523     -0.160            -0.160
## means_y2_scaled                   0.523     -0.160            -0.160
## log_gdp.z.cm                      1.000     -0.577            -0.577
## log_gdp.z.cm_scaled               1.000     -0.577            -0.577
## diff_score                       -0.577      1.000             1.000
## diff_score_scaled                -0.577      1.000             1.000
```

``` r
round(ddsc_mod2_log_GDP$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                       0.580 0.138 32.834   4.205   0.000    0.299    0.861
## w_11                          0.155 0.044 33.032   3.541   0.001    0.066    0.243
## w_21                          0.101 0.048 33.032   2.108   0.043    0.004    0.199
## r_xy1                         0.520 0.147 33.032   3.541   0.001    0.221    0.818
## r_xy2                         0.340 0.161 33.032   2.108   0.043    0.012    0.668
## b_11                          0.519 0.147 33.032   3.541   0.001    0.221    0.817
## b_21                          0.340 0.161 33.032   2.108   0.043    0.012    0.668
## main_effect                   0.128 0.045 33.019   2.814   0.008    0.035    0.220
## moderator_effect             -0.202 0.013 31.719 -16.127   0.000   -0.227   -0.176
## interaction                  -0.053 0.013 32.834  -4.205   0.000   -0.079   -0.028
## q_b11_b21                     0.221    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.222    NA     NA      NA      NA       NA       NA
## cross_over_point             -3.784    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.075 0.051 33.064  -1.455   0.155   -0.179    0.030
## interaction_vs_main_bscale   -0.251 0.172 33.064  -1.455   0.155   -0.601    0.100
## interaction_vs_main_rscale   -0.250 0.172 33.064  -1.453   0.156   -0.600    0.100
## dadas                        -0.202 0.096 33.032  -2.108   0.979   -0.398   -0.007
## dadas_bscale                 -0.680 0.323 33.032  -2.108   0.979   -1.337   -0.024
## dadas_rscale                 -0.680 0.322 33.032  -2.108   0.979   -1.335   -0.024
## abs_diff                      0.053 0.013 32.834   4.205   0.000    0.028    0.079
## abs_sum                       0.256 0.091 33.019   2.814   0.004    0.071    0.441
## abs_diff_bscale               0.179 0.043 32.834   4.205   0.000    0.092    0.266
## abs_sum_bscale                0.859 0.305 33.019   2.814   0.004    0.238    1.480
## abs_diff_rscale               0.180 0.042 32.826   4.233   0.000    0.093    0.266
## abs_sum_rscale                0.859 0.305 33.019   2.815   0.004    0.238    1.481
```

``` r
round(ddsc_mod2_log_GDP$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
##  0.000 -0.009  0.002  1.000  0.960
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
## r_xy1_y2                         0.577 0.142  4.058  0.000    0.298    0.856
## r_xy1                            0.523 0.148  3.526  0.000    0.232    0.814
## r_xy2                            0.346 0.163  2.117  0.034    0.026    0.666
## b_11                             0.524 0.148  3.526  0.000    0.233    0.815
## b_21                             0.346 0.163  2.117  0.034    0.026    0.665
## b_10                             0.335 0.146  2.293  0.022    0.049    0.622
## b_20                            -0.343 0.161 -2.136  0.033   -0.658   -0.028
## res_cov_y1_y2                    0.748 0.188  3.988  0.000    0.380    1.116
## diff_b10_b20                     0.679 0.043 15.703  0.000    0.594    0.763
## diff_b11_b21                     0.178 0.044  4.058  0.000    0.092    0.264
## diff_rxy1_rxy2                   0.177 0.044  4.031  0.000    0.091    0.263
## q_b11_b21                        0.221 0.055  4.003  0.000    0.113    0.329
## q_rxy1_rxy2                      0.220 0.055  3.999  0.000    0.112    0.328
## cross_over_point                -3.811 0.970 -3.929  0.000   -5.712   -1.910
## sum_b11_b21                      0.869 0.309  2.813  0.005    0.264    1.475
## main_effect                      0.435 0.154  2.813  0.005    0.132    0.737
## interaction_vs_main_effect      -0.256 0.174 -1.472  0.141   -0.598    0.085
## diff_abs_b11_abs_b21             0.178 0.044  4.058  0.000    0.092    0.264
## abs_diff_b11_b21                 0.178 0.044  4.058  0.000    0.092    0.264
## abs_sum_b11_b21                  0.869 0.309  2.813  0.002    0.264    1.475
## dadas                           -0.691 0.326 -2.117  0.983   -1.331   -0.051
## q_r_equivalence                  0.120 0.055  2.180  0.985       NA       NA
## q_b_equivalence                  0.121 0.055  2.191  0.986       NA       NA
## cross_over_point_equivalence     3.811 0.970  3.929  1.000       NA       NA
## cross_over_point_minimal_effect  3.811 0.970  3.929  0.000       NA       NA
```

``` r
round(ddsc_sem_log_GDP$variance_test,3)
```

```
##             est    se      z pvalue ci.lower ci.upper
## cov_y1y2  0.924 0.233  3.962  0.000    0.467    1.380
## var_y1    0.971 0.239  4.062  0.000    0.503    1.440
## var_y2    0.968 0.238  4.062  0.000    0.501    1.435
## var_diff  0.004 0.103  0.034  0.973   -0.198    0.205
## var_ratio 1.004 0.107  9.419  0.000    0.795    1.212
## cor_y1y2  0.952 0.016 58.835  0.000    0.921    0.984
```

### Figure


``` r
# plot the results

# start with obtaining predicted values for means and differences


mod2_log_GDP_unstd<-lmer(ben.z~gndr.c+log_gdp.cm+
                           gndr.c:log_gdp.cm+
                       (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
                     control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

mod2_log_GDP_unstd_red<-lmer(ben.z~gndr.c+
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


p1.ben.flags<-
  ggplot(p,aes(y=yvar,x=log_gdp.cm,color=gndr.c))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("log(GDP per capita) (Average 2002-2022)")+
  #ylim(c(min.y.pred,max.y.pred))+
  #xlim(c(0.60,1.00))+
  ylab("Value benevolence (Average 2002-2022)")+
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

p2.ben.flags<-ggplot(p2,aes(y=yvar,x=log_gdp.cm))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("log(GDP per capita) (Average 2002-2022)")+
  #xlim(c(0.60,1.00))+
  #ylim(c(min.y.narrow,max.y.narrow))+
  ylab("Difference in value benevolence")+
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
  ggarrange(p1.ben.flags,p2.ben.flags,align = "v",
            ncol=1,nrow=2,heights=c(2,1))
```

```
## Warning in geom_text(inherit.aes = F, aes(x = 9.75, y = -0.45, label = coef_q, : All aesthetics have length 1, but the data has 3302 rows.
## ℹ Please consider using `annotate()` or provide this layer with data containing a single row.
```

``` r
pflag_comb
```

![](Analysis_benevolence_files/figure-html/unnamed-chunk-27-1.png)<!-- -->

``` r
png(filename = 
      "../results/ben/log_GDP_flags.png",
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
mod3<-lmer(ben.z~gndr.c+essround.c+(gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod3)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ben.z ~ gndr.c + essround.c + (gndr.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1306213.8 1306290.7 -653099.9 1306199.8    441161 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.5974 -0.5187  0.0684  0.6605  5.5016 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr
##  cntry    (Intercept) 0.087068 0.29507      
##           gndr.c      0.007706 0.08778  0.01
##  Residual             0.996111 0.99805      
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##               Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept) -4.141e-03  5.140e-02  3.300e+01  -0.081    0.936    
## gndr.c      -2.028e-01  1.568e-02  3.231e+01 -12.933 2.54e-14 ***
## essround.c   2.373e-02  5.737e-04  4.411e+05  41.370  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##            (Intr) gndr.c
## gndr.c      0.009       
## essround.c -0.001 -0.001
```

``` r
getFE(mod3,round=3)
```

```
##               Est.    SE         df       t     p     LL     UL
## (Intercept) -0.004 0.051     33.003  -0.081 0.936 -0.109  0.100
## gndr.c      -0.203 0.016     32.310 -12.933 0.000 -0.235 -0.171
## essround.c   0.024 0.001 441064.285  41.370 0.000  0.023  0.025
```

``` r
getVC(mod3)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.30 0.09
## 2    cntry      gndr.c   <NA>  0.09 0.01
## 3    cntry (Intercept) gndr.c  0.01 0.00
## 4 Residual        <NA>   <NA>  1.00 1.00
```

``` r
r2mlm(mod3,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.013100882
## slope variation 0.001741495
## mean variation  0.079183644
## sigma2          0.905973979
## 
## $R2s
##           total
## f   0.013100882
## v   0.001741495
## m   0.079183644
## fv  0.014842377
## fvm 0.094026021
```

``` r
anova(mod2,mod3)
```

```
## Data: diff_dat
## Models:
## mod2: ben.z ~ gndr.c + (gndr.c | cntry)
## mod3: ben.z ~ gndr.c + essround.c + (gndr.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod2    6 1307920 1307986 -653954   1307908                         
## mod3    7 1306214 1306291 -653100   1306200 1708.1  1  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## mod4: random effect of time (Ess round)


``` r
mod4<-lmer(ben.z~gndr.c+essround.c+(gndr.c+essround.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod4)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ben.z ~ gndr.c + essround.c + (gndr.c + essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##   1302876   1302986   -651428   1302856    441158 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.7780 -0.5170  0.0804  0.6594  5.2421 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr       
##  cntry    (Intercept) 0.080986 0.28458             
##           gndr.c      0.007849 0.08860  -0.03      
##           essround.c  0.001416 0.03764   0.35 -0.63
##  Residual             0.988304 0.99414             
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##              Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)  0.007380   0.049651 33.062243   0.149   0.8827    
## gndr.c      -0.202256   0.015808 32.372875 -12.795  3.3e-14 ***
## essround.c   0.013381   0.006634 31.251645   2.017   0.0523 .  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##            (Intr) gndr.c
## gndr.c     -0.032       
## essround.c  0.339 -0.607
```

``` r
getFE(mod4,round=3)
```

```
##               Est.    SE     df       t     p     LL     UL
## (Intercept)  0.007 0.050 33.062   0.149 0.883 -0.094  0.108
## gndr.c      -0.202 0.016 32.373 -12.795 0.000 -0.234 -0.170
## essround.c   0.013 0.007 31.252   2.017 0.052  0.000  0.027
```

``` r
getVC(mod4)
```

```
##        grp        var1       var2 sdcor vcov
## 1    cntry (Intercept)       <NA>  0.28 0.08
## 2    cntry      gndr.c       <NA>  0.09 0.01
## 3    cntry  essround.c       <NA>  0.04 0.00
## 4    cntry (Intercept)     gndr.c -0.03 0.00
## 5    cntry (Intercept) essround.c  0.35 0.00
## 6    cntry      gndr.c essround.c -0.63 0.00
## 7 Residual        <NA>       <NA>  0.99 0.99
```

``` r
r2mlm(mod4,bargraph = F)
```

```
## $Decompositions
##                      total
## fixed           0.01050045
## slope variation 0.01154437
## mean variation  0.07398860
## sigma2          0.90396658
## 
## $R2s
##          total
## f   0.01050045
## v   0.01154437
## m   0.07398860
## fv  0.02204482
## fvm 0.09603342
```

``` r
anova(mod2,mod3,mod4)
```

```
## Data: diff_dat
## Models:
## mod2: ben.z ~ gndr.c + (gndr.c | cntry)
## mod3: ben.z ~ gndr.c + essround.c + (gndr.c | cntry)
## mod4: ben.z ~ gndr.c + essround.c + (gndr.c + essround.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod2    6 1307920 1307986 -653954   1307908                         
## mod3    7 1306214 1306291 -653100   1306200 1708.1  1  < 2.2e-16 ***
## mod4   10 1302876 1302986 -651428   1302856 3343.7  3  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## mod5: fixed interaction between time and gender


``` r
mod5<-lmer(ben.z~gndr.c+essround.c+
             gndr.c:essround.c+(gndr.c+essround.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod5)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ben.z ~ gndr.c + essround.c + gndr.c:essround.c + (gndr.c + essround.c |      cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##   1302878   1302999   -651428   1302856    441157 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.7780 -0.5170  0.0804  0.6594  5.2421 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr       
##  cntry    (Intercept) 0.080977 0.28456             
##           gndr.c      0.007850 0.08860  -0.03      
##           essround.c  0.001416 0.03764   0.35 -0.63
##  Residual             0.988304 0.99414             
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##                     Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)        7.380e-03  4.965e-02  3.307e+01   0.149   0.8827    
## gndr.c            -2.023e-01  1.581e-02  3.236e+01 -12.794 3.32e-14 ***
## essround.c         1.338e-02  6.634e-03  3.125e+01   2.017   0.0523 .  
## gndr.c:essround.c  2.511e-06  1.138e-03  1.954e+05   0.002   0.9982    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c essrn.
## gndr.c      -0.032              
## essround.c   0.339 -0.607       
## gndr.c:ssr. -0.002 -0.007  0.006
```

``` r
getFE(mod5,round=3)
```

```
##                     Est.    SE         df       t     p     LL     UL
## (Intercept)        0.007 0.050     33.068   0.149 0.883 -0.094  0.108
## gndr.c            -0.202 0.016     32.363 -12.794 0.000 -0.234 -0.170
## essround.c         0.013 0.007     31.251   2.017 0.052  0.000  0.027
## gndr.c:essround.c  0.000 0.001 195385.476   0.002 0.998 -0.002  0.002
```

``` r
getVC(mod5)
```

```
##        grp        var1       var2 sdcor vcov
## 1    cntry (Intercept)       <NA>  0.28 0.08
## 2    cntry      gndr.c       <NA>  0.09 0.01
## 3    cntry  essround.c       <NA>  0.04 0.00
## 4    cntry (Intercept)     gndr.c -0.03 0.00
## 5    cntry (Intercept) essround.c  0.35 0.00
## 6    cntry      gndr.c essround.c -0.63 0.00
## 7 Residual        <NA>       <NA>  0.99 0.99
```

``` r
r2mlm(mod5,bargraph = F)
```

```
## $Decompositions
##                      total
## fixed           0.01050057
## slope variation 0.01154451
## mean variation  0.07398098
## sigma2          0.90397393
## 
## $R2s
##          total
## f   0.01050057
## v   0.01154451
## m   0.07398098
## fv  0.02204509
## fvm 0.09602607
```

``` r
anova(mod4,mod5)
```

```
## Data: diff_dat
## Models:
## mod4: ben.z ~ gndr.c + essround.c + (gndr.c + essround.c | cntry)
## mod5: ben.z ~ gndr.c + essround.c + gndr.c:essround.c + (gndr.c + essround.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L) Chisq Df Pr(>Chisq)
## mod4   10 1302876 1302986 -651428   1302856                    
## mod5   11 1302878 1302999 -651428   1302856     0  1     0.9982
```

## mod6: random interaction between time and gender


``` r
mod6<-lmer(ben.z~gndr.c+essround.c+
             gndr.c:essround.c+(gndr.c+essround.c+gndr.c:essround.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ben.z ~ gndr.c + essround.c + gndr.c:essround.c + (gndr.c + essround.c +  
##     gndr.c:essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1302817.2 1302982.2 -651393.6 1302787.2    441153 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.7739 -0.5184  0.0764  0.6580  5.2628 
## 
## Random effects:
##  Groups   Name              Variance  Std.Dev. Corr             
##  cntry    (Intercept)       0.0810029 0.28461                   
##           gndr.c            0.0073570 0.08577  -0.02            
##           essround.c        0.0014141 0.03760   0.35 -0.61      
##           gndr.c:essround.c 0.0001366 0.01169  -0.02 -0.47  0.06
##  Residual                   0.9880731 0.99402                   
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##                     Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)        0.0074331  0.0496558 33.0782331   0.150   0.8819    
## gndr.c            -0.2067726  0.0154541 32.0455877 -13.380 1.15e-14 ***
## essround.c         0.0133788  0.0066284 31.3350636   2.018   0.0522 .  
## gndr.c:essround.c  0.0000176  0.0024506 30.5057910   0.007   0.9943    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c essrn.
## gndr.c      -0.014              
## essround.c   0.338 -0.583       
## gndr.c:ssr. -0.019 -0.407  0.052
```

``` r
getFE(mod6,round=3)
```

```
##                     Est.    SE     df       t     p     LL     UL
## (Intercept)        0.007 0.050 33.078   0.150 0.882 -0.094  0.108
## gndr.c            -0.207 0.015 32.046 -13.380 0.000 -0.238 -0.175
## essround.c         0.013 0.007 31.335   2.018 0.052  0.000  0.027
## gndr.c:essround.c  0.000 0.002 30.506   0.007 0.994 -0.005  0.005
```

``` r
getVC(mod6)
```

```
##         grp              var1              var2 sdcor vcov
## 1     cntry       (Intercept)              <NA>  0.28 0.08
## 2     cntry            gndr.c              <NA>  0.09 0.01
## 3     cntry        essround.c              <NA>  0.04 0.00
## 4     cntry gndr.c:essround.c              <NA>  0.01 0.00
## 5     cntry       (Intercept)            gndr.c -0.02 0.00
## 6     cntry       (Intercept)        essround.c  0.35 0.00
## 7     cntry       (Intercept) gndr.c:essround.c -0.02 0.00
## 8     cntry            gndr.c        essround.c -0.61 0.00
## 9     cntry            gndr.c gndr.c:essround.c -0.47 0.00
## 10    cntry        essround.c gndr.c:essround.c  0.06 0.00
## 11 Residual              <NA>              <NA>  0.99 0.99
```

``` r
r2mlm(mod6,bargraph = F)
```

```
## $Decompositions
##                      total
## fixed           0.01091586
## slope variation 0.01164249
## mean variation  0.07394992
## sigma2          0.90349174
## 
## $R2s
##          total
## f   0.01091586
## v   0.01164249
## m   0.07394992
## fv  0.02255834
## fvm 0.09650826
```

``` r
anova(mod4,mod5,mod6)
```

```
## Data: diff_dat
## Models:
## mod4: ben.z ~ gndr.c + essround.c + (gndr.c + essround.c | cntry)
## mod5: ben.z ~ gndr.c + essround.c + gndr.c:essround.c + (gndr.c + essround.c | cntry)
## mod6: ben.z ~ gndr.c + essround.c + gndr.c:essround.c + (gndr.c + essround.c + gndr.c:essround.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod4   10 1302876 1302986 -651428   1302856                         
## mod5   11 1302878 1302999 -651428   1302856  0.000  1     0.9982    
## mod6   15 1302817 1302982 -651394   1302787 68.807  4  4.052e-14 ***
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
##         4.5  0.1710 0.0684 Inf    0.0368    0.3051   2.498  0.0125
##        -4.5  0.0507 0.0468 Inf   -0.0411    0.1425   1.081  0.2795
## 
## gndr.c =  0.5:
##  essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##         4.5 -0.0357 0.0643 Inf   -0.1618    0.0904  -0.555  0.5788
##        -4.5 -0.1562 0.0526 Inf   -0.2593   -0.0531  -2.970  0.0030
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
##  essround.c4.5 - (essround.c-4.5)     0.12 0.0601 Inf  0.002530     0.238   2.002  0.0453
## 
## gndr.c =  0.5:
##  contrast                         estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  essround.c4.5 - (essround.c-4.5)     0.12 0.0612 Inf  0.000493     0.240   1.968  0.0491
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
##    -0.5        4.5  0.1710 0.0684 Inf    0.0368    0.3051   2.498  0.0125
##     0.5        4.5 -0.0357 0.0643 Inf   -0.1618    0.0904  -0.555  0.5788
##    -0.5       -4.5  0.0507 0.0468 Inf   -0.0411    0.1425   1.081  0.2795
##     0.5       -4.5 -0.1562 0.0526 Inf   -0.2593   -0.0531  -2.970  0.0030
## 
## Degrees-of-freedom method: asymptotic 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod6,adjust="none",infer=c(T,T))
```

```
##  contrast                                                 estimate     SE  df asymp.LCL asymp.UCL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5       0.2067 0.0149 Inf  0.177502    0.2359
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   0.1203 0.0601 Inf  0.002530    0.2381
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)    0.3272 0.0698 Inf  0.190363    0.4640
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)     -0.0864 0.0522 Inf -0.188623    0.0159
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)       0.1205 0.0612 Inf  0.000493    0.2405
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)   0.2069 0.0223 Inf  0.163067    0.2506
##  z.ratio p.value
##   13.878  <.0001
##    2.002  0.0453
##    4.687  <.0001
##   -1.655  0.0979
##    1.968  0.0491
##    9.259  <.0001
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
##  diff_ESS10   -0.207 0.0149 Inf    -0.236    -0.178 -13.878  <.0001
##  diff_ESS1    -0.207 0.0223 Inf    -0.251    -0.163  -9.259  <.0001
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
##  diff_ESS10 - diff_ESS1 0.000158 0.0221 Inf   -0.0431    0.0434   0.007  0.9943
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
      obs_mean_wt=weighted.mean(x=ben.z,w=pspwght),
      obs_sd_wt=sqrt(wtd.var(ben.z,pspwght)),
      obs_mean=mean(ben.z),
      obs_sd=sd(ben.z),
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
  ylab("Mean-level of value benevolence")+
  theme(legend.title=element_blank(),
        legend.text  = element_text(size = 16),
        axis.title.x = element_text(size = 18, face = "bold"),   
        axis.title.y = element_text(size = 18, face = "bold"),   
        axis.text.x  = element_text(size = 14),                  
        axis.text.y  = element_text(size = 14))

p_time_trends
```

![](Analysis_benevolence_files/figure-html/unnamed-chunk-33-1.png)<!-- -->

``` r
png(filename = 
      "../results/ben/time_trends.png",
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
pred_cntry_dat$ben.z_mean<-predict(mod6,newdata=pred_cntry_dat)

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

range(pred_cntry_dat$ben.z_mean)
```

```
## [1] -1.0963818  0.6557916
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
      obs_mean_wt=weighted.mean(x=ben.z,w=pspwght),
      obs_sd_wt=sqrt(wtd.var(ben.z,pspwght)),
      obs_mean=mean(ben.z),
      obs_sd=sd(ben.z),
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

pdf("../results/ben/country_specific_time_trends.pdf", width = 8, height = 6)

for (ctry in countries) {
  print(
    ggplot(pred_cntry_dat[pred_cntry_dat$cntry == ctry, ], 
           aes(x = year, y = ben.z_mean, color = gender)) +
      geom_smooth(method = "lm",formula = "y ~ x", se = FALSE) +
      geom_point(size=8) +
      geom_point(aes(x=year,y=obs_mean_wt),size=8,shape = 1,alpha=.50)+
      geom_errorbar(aes(ymin=obs_mean_wt_LL, ymax=obs_mean_wt_UL),alpha = .50)+
      scale_color_manual(values = my_colors) +
      geom_flag(aes(country=tolower(ctry)))+
      scale_y_continuous(limits = c(-1.1, 1.1)) +
      scale_x_continuous(limits = c(2001, 2021), breaks = seq(2002, 2020, 2)) +
      ggtitle(paste("Country:", ctry))+
    ylab("Mean-level of value benevolence")+
      theme(legend.title=element_blank())
  )
}
```

```
## Warning: Removed 2 rows containing missing values or values outside the scale range (`geom_point()`).
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

#mod6_cntry_coefs

#round((mod6_cntry_coefs$`gndr.c:essround.c`/2)*18,2)

#my_colors <- met.brewer("Cassatt2")[c(8, 3)]
my_colors <- met.brewer("Archambault")[c(2,6)]

facet_plot<-
  ggplot(pred_cntry_dat, 
         aes(x = year, y = ben.z_mean, color = gender)) +
  geom_smooth(method = "lm", se = FALSE, formula="y~x") +
  geom_point(size=3.5) +
  geom_point(aes(x=year,y=obs_mean_wt),size=3.5,shape = 1,alpha=.70)+
  geom_errorbar(aes(ymin=obs_mean_wt_LL, ymax=obs_mean_wt_UL),alpha = .50)+
  scale_color_manual(values = my_colors) +
  #geom_flag(aes(country=tolower(ctry)))+
  scale_y_continuous(limits = c(-1.1, 1.1)) +
  scale_x_continuous(limits = c(2001, 2021), breaks = seq(2002, 2020, 2)) +
  #ggtitle(paste("Country:", ctry))+
  ylab("Mean-level of value benevolence")+
  xlab("Year")+
  theme(legend.title=element_blank(),legend.position = "top",
        axis.text.x = element_text(angle = 45,size = 6,hjust=1))+
  facet_wrap(~CLDR,nrow=6,ncol=6)+
  #facet_wrap(~cntry,nrow=6,ncol=6)+
  geom_flag(aes(country=tolower(cntry)),size=2)

facet_plot
```

```
## Warning: Removed 2 rows containing missing values or values outside the scale range (`geom_point()`).
```

![](Analysis_benevolence_files/figure-html/unnamed-chunk-35-1.png)<!-- -->

``` r
png(filename = 
      "../results/ben/country_time_trend_facets.png",
    units = "cm",
    width = 21.0,height=29.7*(3/4),res = 600)
facet_plot
```

```
## Warning: Removed 2 rows containing missing values or values outside the scale range (`geom_point()`).
```

``` r
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
## [1] 10.35406
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
## 1   -0.24                0.25                   -0.06                    0.22                      0.28
## 2   -0.23                0.20                    0.00                    0.20                      0.20
## 3   -0.09               -0.34                    0.10                   -0.29                     -0.39
## 4   -0.22                0.33                   -0.07                    0.29                      0.36
## 5   -0.16                0.21                   -0.03                    0.20                      0.23
## 6   -0.29               -0.17                    0.10                   -0.12                     -0.22
## 7   -0.22                0.56                    0.02                    0.57                      0.55
## 8   -0.28                0.32                   -0.02                    0.31                      0.33
## 9   -0.29                0.22                   -0.03                    0.20                      0.23
## 10  -0.15                0.34                   -0.08                    0.30                      0.38
## 11  -0.37                0.58                    0.12                    0.63                      0.52
## 12  -0.27                0.26                   -0.07                    0.23                      0.30
## 13  -0.25                0.38                   -0.01                    0.38                      0.39
## 14  -0.09               -0.42                   -0.19                   -0.52                     -0.33
## 15  -0.18                0.20                   -0.02                    0.19                      0.21
## 16  -0.19               -0.22                    0.16                   -0.14                     -0.30
## 17  -0.25                0.16                    0.09                    0.20                      0.11
## 18  -0.09                0.06                    0.02                    0.07                      0.05
## 19  -0.33                0.33                    0.13                    0.39                      0.27
## 20  -0.17               -0.37                    0.05                   -0.35                     -0.40
## 21  -0.23               -0.25                    0.06                   -0.22                     -0.28
## 22  -0.29               -0.01                    0.03                    0.01                     -0.03
## 23  -0.15               -0.32                    0.01                   -0.32                     -0.32
## 24  -0.27                0.38                    0.00                    0.38                      0.38
## 25  -0.27                0.56                    0.05                    0.59                      0.54
## 26  -0.16                0.07                   -0.03                    0.05                      0.09
## 27  -0.09                0.27                   -0.22                    0.16                      0.38
## 28  -0.10               -0.40                    0.01                   -0.40                     -0.41
## 29  -0.35                0.76                    0.14                    0.83                      0.69
## 30  -0.22                0.56                    0.00                    0.56                      0.56
## 31  -0.19               -0.09                    0.03                   -0.07                     -0.10
## 32  -0.02               -0.55                   -0.08                   -0.59                     -0.52
## 33  -0.13                0.12                   -0.18                    0.03                      0.21
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
## 1     GR               -0.42
## 2     RU               -0.40
## 3     BG               -0.34
## 4     LT               -0.25
## 5     HU               -0.22
## 6     CZ               -0.17
## 7     SK               -0.09
## 8     IL                0.06
## 9     PL                0.07
## 10    UA                0.12
## 11    IE                0.16
## 12    BE                0.20
## 13    CY                0.21
## 14    EE                0.22
## 15    AT                0.25
## 16    FR                0.26
## 17    PT                0.27
## 18    DK                0.32
## 19    CH                0.33
## 20    IS                0.33
## 21    ES                0.34
## 22    GB                0.38
## 23    NL                0.38
## 24    DE                0.56
## 25    NO                0.56
## 26    SI                0.56
## 27    FI                0.58
## 28    SE                0.76
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
## 1     PT                   -0.22
## 2     GR                   -0.19
## 3     UA                   -0.18
## 4     ES                   -0.08
## 5     CH                   -0.07
## 6     FR                   -0.07
## 7     AT                   -0.06
## 8     CY                   -0.03
## 9     EE                   -0.03
## 10    PL                   -0.03
## 11    DK                   -0.02
## 12    GB                   -0.01
## 13    BE                    0.00
## 14    NL                    0.00
## 15    SI                    0.00
## 16    RU                    0.01
## 17    DE                    0.02
## 18    IL                    0.02
## 19    SK                    0.03
## 20    NO                    0.05
## 21    LT                    0.06
## 22    IE                    0.09
## 23    BG                    0.10
## 24    CZ                    0.10
## 25    FI                    0.12
## 26    IS                    0.13
## 27    SE                    0.14
## 28    HU                    0.16
```


#### Country-level predictors: GEI


``` r
mod6_GEI<-lmer(ben.z~gndr.c+essround.c+
             gndr.c:essround.c+
               gei.z.cm:gndr.c+gei.z.cm:essround.c+gei.z.cm:gndr.c:essround.c+
               (gndr.c+essround.c+gndr.c:essround.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_GEI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ben.z ~ gndr.c + essround.c + gndr.c:essround.c + gei.z.cm:gndr.c +  
##     gei.z.cm:essround.c + gei.z.cm:gndr.c:essround.c + (gndr.c +  
##     essround.c + gndr.c:essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1266903.3 1267100.8 -633433.6 1266867.3    431760 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.8433 -0.5180  0.0788  0.6606  5.3095 
## 
## Random effects:
##  Groups   Name              Variance  Std.Dev. Corr             
##  cntry    (Intercept)       0.0750825 0.27401                   
##           gndr.c            0.0047654 0.06903   0.34            
##           essround.c        0.0006260 0.02502   0.23 -0.31      
##           gndr.c:essround.c 0.0001118 0.01057  -0.14 -0.61  0.20
##  Residual                   0.9740873 0.98696                   
## Number of obs: 431778, groups:  cntry, 32
## 
## Fixed effects:
##                             Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)                 0.021463   0.048546 32.082389   0.442  0.66136    
## gndr.c                     -0.207183   0.012829 30.360190 -16.149  < 2e-16 ***
## essround.c                  0.013205   0.004540 29.111941   2.909  0.00688 ** 
## gndr.c:essround.c           0.001545   0.002325 34.212501   0.665  0.51079    
## gndr.c:gei.z.cm            -0.064743   0.012564 33.297429  -5.153 1.15e-05 ***
## essround.c:gei.z.cm         0.027983   0.004648 31.500712   6.020 1.09e-06 ***
## gndr.c:essround.c:gei.z.cm -0.001889   0.002569 43.082891  -0.735  0.46606    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c essrn. gnd.:. gn.:.. es.:..
## gndr.c       0.323                                   
## essround.c   0.214 -0.292                            
## gndr.c:ssr. -0.113 -0.502  0.168                     
## gndr.c:g.z.  0.000 -0.029  0.003 -0.031              
## essrnd.c:.. -0.004  0.002 -0.041 -0.002 -0.377       
## gndr.c:.:..  0.001 -0.015 -0.005 -0.148 -0.329  0.169
```

``` r
getFE(mod6_GEI,round=2,p.round = 5)
```

```
##                             Est.   SE    df      t       p    LL    UL
## (Intercept)                 0.02 0.05 32.08   0.44 0.66136 -0.08  0.12
## gndr.c                     -0.21 0.01 30.36 -16.15 0.00000 -0.23 -0.18
## essround.c                  0.01 0.00 29.11   2.91 0.00688  0.00  0.02
## gndr.c:essround.c           0.00 0.00 34.21   0.66 0.51079  0.00  0.01
## gndr.c:gei.z.cm            -0.06 0.01 33.30  -5.15 0.00001 -0.09 -0.04
## essround.c:gei.z.cm         0.03 0.00 31.50   6.02 0.00000  0.02  0.04
## gndr.c:essround.c:gei.z.cm  0.00 0.00 43.08  -0.74 0.46606 -0.01  0.00
```

``` r
getVC(mod6_GEI)
```

```
##         grp              var1              var2 sdcor vcov
## 1     cntry       (Intercept)              <NA>  0.27 0.08
## 2     cntry            gndr.c              <NA>  0.07 0.00
## 3     cntry        essround.c              <NA>  0.03 0.00
## 4     cntry gndr.c:essround.c              <NA>  0.01 0.00
## 5     cntry       (Intercept)            gndr.c  0.34 0.01
## 6     cntry       (Intercept)        essround.c  0.23 0.00
## 7     cntry       (Intercept) gndr.c:essround.c -0.14 0.00
## 8     cntry            gndr.c        essround.c -0.31 0.00
## 9     cntry            gndr.c gndr.c:essround.c -0.61 0.00
## 10    cntry        essround.c gndr.c:essround.c  0.20 0.00
## 11 Residual              <NA>              <NA>  0.99 0.97
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_GEI,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] 55.72863
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_GEI,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] 18.17048
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
##         4.5 -0.0450 0.0607 Inf   -0.1639    0.0738  -0.743  0.4578
##        -4.5  0.0880 0.0530 Inf   -0.0160    0.1919   1.658  0.0973
## 
## gei.z.cm =  0:
##  essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##         4.5  0.0809 0.0566 Inf   -0.0300    0.1917   1.430  0.1527
##        -4.5 -0.0380 0.0485 Inf   -0.1330    0.0570  -0.783  0.4336
## 
## gei.z.cm =  1:
##  essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##         4.5  0.2068 0.0599 Inf    0.0893    0.3243   3.450  0.0006
##        -4.5 -0.1639 0.0525 Inf   -0.2668   -0.0609  -3.119  0.0018
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
##  essround.c4.5 - (essround.c-4.5)   -0.133 0.0597 Inf   -0.2500   -0.0161  -2.229  0.0258
## 
## gei.z.cm =  0:
##  contrast                         estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  essround.c4.5 - (essround.c-4.5)    0.119 0.0409 Inf    0.0388    0.1989   2.909  0.0036
## 
## gei.z.cm =  1:
##  contrast                         estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  essround.c4.5 - (essround.c-4.5)    0.371 0.0573 Inf    0.2585    0.4829   6.473  <.0001
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
##  gndr.c essround.c   emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##    -0.5        4.5  0.01845 0.0611 Inf   -0.1013    0.1382   0.302  0.7626
##     0.5        4.5 -0.10853 0.0618 Inf   -0.2297    0.0126  -1.756  0.0792
##    -0.5       -4.5  0.16691 0.0501 Inf    0.0688    0.2650   3.334  0.0009
##     0.5       -4.5  0.00902 0.0595 Inf   -0.1075    0.1256   0.152  0.8795
## 
## gei.z.cm =  0:
##  gndr.c essround.c   emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##    -0.5        4.5  0.18100 0.0559 Inf    0.0714    0.2906   3.235  0.0012
##     0.5        4.5 -0.01923 0.0578 Inf   -0.1324    0.0940  -0.333  0.7392
##    -0.5       -4.5  0.06911 0.0456 Inf   -0.0203    0.1585   1.515  0.1298
##     0.5       -4.5 -0.14503 0.0531 Inf   -0.2492   -0.0409  -2.730  0.0063
## 
## gei.z.cm =  1:
##  gndr.c essround.c   emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##    -0.5        4.5  0.34354 0.0601 Inf    0.2257    0.4614   5.713  <.0001
##     0.5        4.5  0.07007 0.0609 Inf   -0.0493    0.1894   1.151  0.2499
##    -0.5       -4.5 -0.02869 0.0494 Inf   -0.1256    0.0682  -0.580  0.5616
##     0.5       -4.5 -0.29907 0.0588 Inf   -0.4143   -0.1838  -5.086  <.0001
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
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      0.12698 0.0198 Inf   0.08813   0.16584
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5) -0.14846 0.0594 Inf  -0.26484  -0.03209
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   0.00943 0.0678 Inf  -0.12350   0.14237
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)    -0.27545 0.0564 Inf  -0.38603  -0.16486
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)     -0.11755 0.0644 Inf  -0.24385   0.00876
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  0.15790 0.0288 Inf   0.10147   0.21432
##  z.ratio p.value
##    6.406  <.0001
##   -2.500  0.0124
##    0.139  0.8894
##   -4.882  <.0001
##   -1.824  0.0681
##    5.485  <.0001
## 
## gei.z.cm =  0:
##  contrast                                                 estimate     SE  df asymp.LCL asymp.UCL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      0.20023 0.0118 Inf   0.17711   0.22335
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  0.11189 0.0404 Inf   0.03264   0.19114
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   0.32602 0.0463 Inf   0.23535   0.41669
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)    -0.08834 0.0391 Inf  -0.16496  -0.01173
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      0.12580 0.0439 Inf   0.03985   0.21174
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  0.21414 0.0202 Inf   0.17450   0.25378
##  z.ratio p.value
##   16.972  <.0001
##    2.767  0.0057
##    7.047  <.0001
##   -2.260  0.0238
##    2.869  0.0041
##   10.588  <.0001
## 
## gei.z.cm =  1:
##  contrast                                                 estimate     SE  df asymp.LCL asymp.UCL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      0.27347 0.0167 Inf   0.24084   0.30611
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  0.37224 0.0564 Inf   0.26162   0.48286
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   0.64261 0.0655 Inf   0.51426   0.77097
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)     0.09876 0.0538 Inf  -0.00671   0.20423
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      0.36914 0.0615 Inf   0.24852   0.48976
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  0.27038 0.0276 Inf   0.21622   0.32453
##  z.ratio p.value
##   16.422  <.0001
##    6.595  <.0001
##    9.813  <.0001
##    1.835  0.0665
##    5.998  <.0001
##    9.786  <.0001
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
##  diff_ESS10   -0.127 0.0198 Inf    -0.166   -0.0881  -6.406  <.0001
##  diff_ESS1    -0.158 0.0288 Inf    -0.214   -0.1015  -5.485  <.0001
## 
## gei.z.cm =  0:
##  contrast   estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10   -0.200 0.0118 Inf    -0.223   -0.1771 -16.972  <.0001
##  diff_ESS1    -0.214 0.0202 Inf    -0.254   -0.1745 -10.588  <.0001
## 
## gei.z.cm =  1:
##  contrast   estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10   -0.273 0.0167 Inf    -0.306   -0.2408 -16.422  <.0001
##  diff_ESS1    -0.270 0.0276 Inf    -0.325   -0.2162  -9.786  <.0001
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
##  diff_ESS10 - diff_ESS1   0.0309 0.0334 Inf   -0.0346    0.0964   0.925  0.3548
## 
## gei.z.cm =  0:
##  contrast               estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10 - diff_ESS1   0.0139 0.0209 Inf   -0.0271    0.0549   0.665  0.5063
## 
## gei.z.cm =  1:
##  contrast               estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10 - diff_ESS1  -0.0031 0.0288 Inf   -0.0595    0.0533  -0.108  0.9144
## 
## Degrees-of-freedom method: asymptotic 
## Confidence level used: 0.95
```


#### Country-level predictors: GGGI


``` r
mod6_GGGI<-lmer(ben.z~gndr.c+essround.c+
             gndr.c:essround.c+
               gggi.z.cm:gndr.c+gggi.z.cm:essround.c+gggi.z.cm:gndr.c:essround.c+
               (gndr.c+essround.c+gndr.c:essround.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_GGGI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ben.z ~ gndr.c + essround.c + gndr.c:essround.c + gggi.z.cm:gndr.c +  
##     gggi.z.cm:essround.c + gggi.z.cm:gndr.c:essround.c + (gndr.c +  
##     essround.c + gndr.c:essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##  930356.8  930548.7 -465160.4  930320.8    314628 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.8344 -0.5155  0.0706  0.6636  5.2893 
## 
## Random effects:
##  Groups   Name              Variance  Std.Dev. Corr             
##  cntry    (Intercept)       0.0846398 0.29093                   
##           gndr.c            0.0030318 0.05506   0.28            
##           essround.c        0.0010766 0.03281   0.24 -0.22      
##           gndr.c:essround.c 0.0001543 0.01242  -0.29 -0.10 -0.28
##  Residual                   0.9900541 0.99501                   
## Number of obs: 314646, groups:  cntry, 33
## 
## Fixed effects:
##                              Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)                  0.018558   0.050867 33.300342   0.365 0.717549    
## gndr.c                      -0.200160   0.010568 32.869415 -18.940  < 2e-16 ***
## essround.c                   0.013314   0.006069 28.907146   2.194 0.036440 *  
## gndr.c:essround.c           -0.004820   0.002918 30.556557  -1.652 0.108792    
## gndr.c:gggi.z.cm            -0.066191   0.010844 37.943686  -6.104 4.12e-07 ***
## essround.c:gggi.z.cm         0.026957   0.006430 30.477532   4.192 0.000219 ***
## gndr.c:essround.c:gggi.z.cm  0.006544   0.003185 35.523736   2.054 0.047344 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c essrn. gnd.:. gn.:.. es.:..
## gndr.c       0.254                                   
## essround.c   0.227 -0.204                            
## gndr.c:ssr. -0.216 -0.115 -0.205                     
## gndr.c:gg..  0.001 -0.023  0.013 -0.031              
## essrnd.c:.. -0.008  0.011 -0.083  0.014 -0.267       
## gndr.c:.:..  0.003 -0.031  0.013 -0.083 -0.064 -0.169
```

``` r
getFE(mod6_GGGI,round=2,p.round = 5)
```

```
##                              Est.   SE    df      t       p    LL    UL
## (Intercept)                  0.02 0.05 33.30   0.36 0.71755 -0.08  0.12
## gndr.c                      -0.20 0.01 32.87 -18.94 0.00000 -0.22 -0.18
## essround.c                   0.01 0.01 28.91   2.19 0.03644  0.00  0.03
## gndr.c:essround.c            0.00 0.00 30.56  -1.65 0.10879 -0.01  0.00
## gndr.c:gggi.z.cm            -0.07 0.01 37.94  -6.10 0.00000 -0.09 -0.04
## essround.c:gggi.z.cm         0.03 0.01 30.48   4.19 0.00022  0.01  0.04
## gndr.c:essround.c:gggi.z.cm  0.01 0.00 35.52   2.05 0.04734  0.00  0.01
```

``` r
getVC(mod6_GGGI)
```

```
##         grp              var1              var2 sdcor vcov
## 1     cntry       (Intercept)              <NA>  0.29 0.08
## 2     cntry            gndr.c              <NA>  0.06 0.00
## 3     cntry        essround.c              <NA>  0.03 0.00
## 4     cntry gndr.c:essround.c              <NA>  0.01 0.00
## 5     cntry       (Intercept)            gndr.c  0.28 0.00
## 6     cntry       (Intercept)        essround.c  0.24 0.00
## 7     cntry       (Intercept) gndr.c:essround.c -0.29 0.00
## 8     cntry            gndr.c        essround.c -0.22 0.00
## 9     cntry            gndr.c gndr.c:essround.c -0.10 0.00
## 10    cntry        essround.c gndr.c:essround.c -0.28 0.00
## 11 Residual              <NA>              <NA>  1.00 0.99
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_GGGI,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] 23.86341
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_GGGI,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] -12.9457
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
##  essround.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5 -0.0428 0.0704 42.5  -0.1848   0.0992  -0.609  0.5461
##        -4.5  0.0800 0.0604 44.0  -0.0418   0.2017   1.324  0.1924
## 
## gggi.z.cm =  0:
##  essround.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.0785 0.0630 32.2  -0.0497   0.2067   1.246  0.2216
##        -4.5 -0.0414 0.0520 32.2  -0.1472   0.0645  -0.795  0.4322
## 
## gggi.z.cm =  1:
##  essround.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.1998 0.0682 40.0   0.0620   0.3376   2.930  0.0056
##        -4.5 -0.1627 0.0586 42.3  -0.2809  -0.0445  -2.777  0.0082
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
##  essround.c4.5 - (essround.c-4.5)   -0.123 0.0828 30.8  -0.2917   0.0461  -1.483  0.1483
## 
## gggi.z.cm =  0:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)    0.120 0.0546 28.9   0.0081   0.2316   2.194  0.0364
## 
## gggi.z.cm =  1:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)    0.362 0.0762 30.0   0.2068   0.5181   4.755  <.0001
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
##  gndr.c essround.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5        4.5  0.0497 0.0736 42.6  -0.0988  0.19826   0.675  0.5032
##     0.5        4.5 -0.1354 0.0692 41.7  -0.2751  0.00434  -1.956  0.0572
##    -0.5       -4.5  0.1214 0.0595 44.3   0.0015  0.24123   2.040  0.0473
##     0.5       -4.5  0.0385 0.0640 43.9  -0.0905  0.16753   0.602  0.5502
## 
## gggi.z.cm =  0:
##  gndr.c essround.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5        4.5  0.1894 0.0646 31.8   0.0579  0.32091   2.934  0.0062
##     0.5        4.5 -0.0325 0.0623 32.0  -0.1594  0.09453  -0.521  0.6063
##    -0.5       -4.5  0.0479 0.0502 32.0  -0.0543  0.15003   0.955  0.3469
##     0.5       -4.5 -0.1306 0.0552 32.0  -0.2431 -0.01813  -2.365  0.0242
## 
## gggi.z.cm =  1:
##  gndr.c essround.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5        4.5  0.3291 0.0710 40.2   0.1855  0.47262   4.633  <.0001
##     0.5        4.5  0.0705 0.0671 39.1  -0.0652  0.20615   1.051  0.2998
##    -0.5       -4.5 -0.0256 0.0575 42.6  -0.1417  0.09045  -0.445  0.6585
##     0.5       -4.5 -0.2997 0.0622 42.2  -0.4253 -0.17412  -4.815  <.0001
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
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5       0.1851 0.0247 29.8   0.1346  0.23557   7.493
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  -0.0716 0.0888 30.4  -0.2528  0.10955  -0.807
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)    0.0112 0.0877 31.1  -0.1676  0.19001   0.128
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)     -0.2568 0.0806 31.1  -0.4211 -0.09246  -3.187
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      -0.1739 0.0815 30.5  -0.3403 -0.00752  -2.133
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)   0.0828 0.0260 39.3   0.0302  0.13545   3.183
##  p.value
##   <.0001
##   0.4259
##   0.8993
##   0.0033
##   0.0411
##   0.0028
## 
## gggi.z.cm =  0:
##  contrast                                                 estimate     SE   df lower.CL upper.CL t.ratio
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5       0.2219 0.0159 28.4   0.1893  0.25435  13.972
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   0.1415 0.0587 28.3   0.0213  0.26178   2.410
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)    0.3200 0.0577 29.2   0.2020  0.43798   5.545
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)     -0.0803 0.0535 28.8  -0.1897  0.02909  -1.502
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)       0.0981 0.0535 28.5  -0.0114  0.20765   1.834
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)   0.1785 0.0178 31.8   0.1422  0.21469  10.039
##  p.value
##   <.0001
##   0.0227
##   <.0001
##   0.1440
##   0.0771
##   <.0001
## 
## gggi.z.cm =  1:
##  contrast                                                 estimate     SE   df lower.CL upper.CL t.ratio
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5       0.2586 0.0224 28.4   0.2128  0.30439  11.560
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   0.3547 0.0818 29.6   0.1876  0.52180   4.337
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)    0.6288 0.0811 29.9   0.4632  0.79443   7.754
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.0961 0.0741 30.4  -0.0552  0.24736   1.297
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)       0.3702 0.0750 29.3   0.2169  0.52351   4.936
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)   0.2741 0.0253 33.8   0.2227  0.32556  10.830
##  p.value
##   <.0001
##   0.0002
##   <.0001
##   0.2045
##   <.0001
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
##  diff_ESS10  -0.1851 0.0247 29.8   -0.236  -0.1346  -7.493  <.0001
##  diff_ESS1   -0.0828 0.0260 39.3   -0.135  -0.0302  -3.183  0.0028
## 
## gggi.z.cm =  0:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10  -0.2219 0.0159 28.4   -0.254  -0.1893 -13.972  <.0001
##  diff_ESS1   -0.1785 0.0178 31.8   -0.215  -0.1422 -10.039  <.0001
## 
## gggi.z.cm =  1:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10  -0.2586 0.0224 28.4   -0.304  -0.2128 -11.560  <.0001
##  diff_ESS1   -0.2741 0.0253 33.8   -0.326  -0.2227 -10.830  <.0001
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
##  diff_ESS10 - diff_ESS1  -0.1023 0.0405 35.6  -0.1844  -0.0202  -2.528  0.0160
## 
## gggi.z.cm =  0:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10 - diff_ESS1  -0.0434 0.0263 30.6  -0.0970   0.0102  -1.652  0.1088
## 
## gggi.z.cm =  1:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10 - diff_ESS1   0.0155 0.0372 32.0  -0.0603   0.0914   0.417  0.6797
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

#### Country-level predictors: GDI


``` r
mod6_GDI<-lmer(ben.z~gndr.c+essround.c+
             gndr.c:essround.c+
               gdi.z.cm:gndr.c+gdi.z.cm:essround.c+gdi.z.cm:gndr.c:essround.c+
               (gndr.c+essround.c+gndr.c:essround.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_GDI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ben.z ~ gndr.c + essround.c + gndr.c:essround.c + gdi.z.cm:gndr.c +  
##     gdi.z.cm:essround.c + gdi.z.cm:gndr.c:essround.c + (gndr.c +  
##     essround.c + gndr.c:essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1302818.6 1303016.6 -651391.3 1302782.6    441150 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.7738 -0.5183  0.0765  0.6579  5.2628 
## 
## Random effects:
##  Groups   Name              Variance  Std.Dev. Corr             
##  cntry    (Intercept)       0.0809285 0.28448                   
##           gndr.c            0.0070763 0.08412  -0.22            
##           essround.c        0.0015570 0.03946   0.55 -0.61      
##           gndr.c:essround.c 0.0001351 0.01162   0.05 -0.46  0.05
##  Residual                   0.9880729 0.99402                   
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##                              Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)                 0.0067771  0.0496252 33.0780468   0.137   0.8922    
## gndr.c                     -0.2068031  0.0151755 31.3229243 -13.627 1.04e-14 ***
## essround.c                  0.0133154  0.0069452 27.4579487   1.917   0.0657 .  
## gndr.c:essround.c           0.0000825  0.0024401 30.2938781   0.034   0.9733    
## gndr.c:gdi.z.cm            -0.0306861  0.0153612 35.4814372  -1.998   0.0535 .  
## essround.c:gdi.z.cm         0.0154176  0.0060580 32.6745955   2.545   0.0158 *  
## gndr.c:essround.c:gdi.z.cm  0.0014760  0.0028447 44.9659526   0.519   0.6064    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c essrn. gnd.:. gn.:.. es.:..
## gndr.c      -0.207                                   
## essround.c   0.535 -0.588                            
## gndr.c:ssr.  0.037 -0.397  0.046                     
## gndr.c:gd..  0.000 -0.015  0.002  0.000              
## essrnd.c:..  0.001  0.003 -0.018  0.003 -0.558       
## gndr.c:.:..  0.003 -0.003  0.003 -0.010 -0.274  0.020
```

``` r
getFE(mod6_GDI,round=2,p.round = 5)
```

```
##                             Est.   SE    df      t       p    LL    UL
## (Intercept)                 0.01 0.05 33.08   0.14 0.89220 -0.09  0.11
## gndr.c                     -0.21 0.02 31.32 -13.63 0.00000 -0.24 -0.18
## essround.c                  0.01 0.01 27.46   1.92 0.06567  0.00  0.03
## gndr.c:essround.c           0.00 0.00 30.29   0.03 0.97325  0.00  0.01
## gndr.c:gdi.z.cm            -0.03 0.02 35.48  -2.00 0.05347 -0.06  0.00
## essround.c:gdi.z.cm         0.02 0.01 32.67   2.54 0.01584  0.00  0.03
## gndr.c:essround.c:gdi.z.cm  0.00 0.00 44.97   0.52 0.60639  0.00  0.01
```

``` r
getVC(mod6_GDI)
```

```
##         grp              var1              var2 sdcor  vcov
## 1     cntry       (Intercept)              <NA>  0.28  0.08
## 2     cntry            gndr.c              <NA>  0.08  0.01
## 3     cntry        essround.c              <NA>  0.04  0.00
## 4     cntry gndr.c:essround.c              <NA>  0.01  0.00
## 5     cntry       (Intercept)            gndr.c -0.22 -0.01
## 6     cntry       (Intercept)        essround.c  0.55  0.01
## 7     cntry       (Intercept) gndr.c:essround.c  0.05  0.00
## 8     cntry            gndr.c        essround.c -0.61  0.00
## 9     cntry            gndr.c gndr.c:essround.c -0.46  0.00
## 10    cntry        essround.c gndr.c:essround.c  0.05  0.00
## 11 Residual              <NA>              <NA>  0.99  0.99
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_GDI,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] -10.10851
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_GDI,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] 1.102289
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
##         4.5 -0.00268 0.0766 Inf   -0.1528    0.1475  -0.035  0.9721
##        -4.5  0.01624 0.0506 Inf   -0.0829    0.1154   0.321  0.7482
## 
## gdi.z.cm =  0:
##  essround.c   emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##         4.5  0.06670 0.0714 Inf   -0.0733    0.2066   0.934  0.3503
##        -4.5 -0.05314 0.0422 Inf   -0.1358    0.0295  -1.260  0.2078
## 
## gdi.z.cm =  1:
##  essround.c   emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##         4.5  0.13608 0.0763 Inf   -0.0134    0.2855   1.784  0.0744
##        -4.5 -0.12252 0.0499 Inf   -0.2203   -0.0247  -2.456  0.0141
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
##  essround.c4.5 - (essround.c-4.5)  -0.0189 0.0837 Inf  -0.18292     0.145  -0.226  0.8211
## 
## gdi.z.cm =  0:
##  contrast                         estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   0.1198 0.0625 Inf  -0.00267     0.242   1.917  0.0552
## 
## gdi.z.cm =  1:
##  contrast                         estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   0.2586 0.0822 Inf   0.09747     0.420   3.146  0.0017
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
##  gndr.c essround.c   emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##    -0.5        4.5  0.08851 0.0814 Inf  -0.07108    0.2481   1.087  0.2770
##     0.5        4.5 -0.09388 0.0733 Inf  -0.23746    0.0497  -1.281  0.2000
##    -0.5       -4.5  0.10116 0.0494 Inf   0.00425    0.1981   2.046  0.0408
##     0.5       -4.5 -0.06869 0.0563 Inf  -0.17904    0.0417  -1.220  0.2225
## 
## gdi.z.cm =  0:
##  gndr.c essround.c   emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##    -0.5        4.5  0.16991 0.0745 Inf   0.02387    0.3160   2.280  0.0226
##     0.5        4.5 -0.03652 0.0690 Inf  -0.17167    0.0986  -0.530  0.5964
##    -0.5       -4.5  0.05044 0.0422 Inf  -0.03229    0.1332   1.195  0.2321
##     0.5       -4.5 -0.15673 0.0449 Inf  -0.24482   -0.0686  -3.487  0.0005
## 
## gdi.z.cm =  1:
##  gndr.c essround.c   emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##    -0.5        4.5  0.25131 0.0810 Inf   0.09257    0.4101   3.103  0.0019
##     0.5        4.5  0.02084 0.0729 Inf  -0.12214    0.1638   0.286  0.7751
##    -0.5       -4.5 -0.00027 0.0488 Inf  -0.09585    0.0953  -0.006  0.9956
##     0.5       -4.5 -0.24477 0.0556 Inf  -0.35374   -0.1358  -4.403  <.0001
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
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5       0.1824 0.0228 Inf   0.13763    0.2271
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  -0.0126 0.0849 Inf  -0.17902    0.1537
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)    0.1572 0.0976 Inf  -0.03419    0.3486
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)     -0.1950 0.0736 Inf  -0.33925   -0.0508
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      -0.0252 0.0859 Inf  -0.19347    0.1431
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)   0.1698 0.0316 Inf   0.10787    0.2318
##  z.ratio p.value
##    7.986  <.0001
##   -0.149  0.8815
##    1.610  0.1074
##   -2.651  0.0080
##   -0.293  0.7692
##    5.372  <.0001
## 
## gdi.z.cm =  0:
##  contrast                                                 estimate     SE  df asymp.LCL asymp.UCL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5       0.2064 0.0148 Inf   0.17746    0.2354
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   0.1195 0.0630 Inf  -0.00394    0.2429
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)    0.3266 0.0725 Inf   0.18458    0.4687
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)     -0.0870 0.0550 Inf  -0.19469    0.0208
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)       0.1202 0.0640 Inf  -0.00515    0.2456
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)   0.2072 0.0220 Inf   0.16409    0.2503
##  z.ratio p.value
##   13.964  <.0001
##    1.897  0.0578
##    4.507  <.0001
##   -1.582  0.1136
##    1.879  0.0602
##    9.425  <.0001
## 
## gdi.z.cm =  1:
##  contrast                                                 estimate     SE  df asymp.LCL asymp.UCL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5       0.2305 0.0223 Inf   0.18668    0.2743
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   0.2516 0.0833 Inf   0.08828    0.4149
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)    0.4961 0.0962 Inf   0.30750    0.6847
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.0211 0.0719 Inf  -0.11991    0.1621
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)       0.2656 0.0845 Inf   0.10002    0.4312
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)   0.2445 0.0313 Inf   0.18308    0.3059
##  z.ratio p.value
##   10.313  <.0001
##    3.019  0.0025
##    5.156  <.0001
##    0.293  0.7692
##    3.144  0.0017
##    7.802  <.0001
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
##  diff_ESS10   -0.182 0.0228 Inf    -0.227    -0.138  -7.986  <.0001
##  diff_ESS1    -0.170 0.0316 Inf    -0.232    -0.108  -5.372  <.0001
## 
## gdi.z.cm =  0:
##  contrast   estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10   -0.206 0.0148 Inf    -0.235    -0.177 -13.964  <.0001
##  diff_ESS1    -0.207 0.0220 Inf    -0.250    -0.164  -9.425  <.0001
## 
## gdi.z.cm =  1:
##  contrast   estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10   -0.230 0.0223 Inf    -0.274    -0.187 -10.313  <.0001
##  diff_ESS1    -0.245 0.0313 Inf    -0.306    -0.183  -7.802  <.0001
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
##  contrast                estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10 - diff_ESS1 -0.012542 0.0339 Inf   -0.0790    0.0539  -0.370  0.7114
## 
## gdi.z.cm =  0:
##  contrast                estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10 - diff_ESS1  0.000742 0.0220 Inf   -0.0423    0.0438   0.034  0.9730
## 
## gdi.z.cm =  1:
##  contrast                estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10 - diff_ESS1  0.014027 0.0336 Inf   -0.0518    0.0798   0.418  0.6760
## 
## Degrees-of-freedom method: asymptotic 
## Confidence level used: 0.95
```

#### Country-level predictors: log_GDP


``` r
mod6_log_GDP<-lmer(ben.z~gndr.c+essround.c+
             gndr.c:essround.c+
               log_gdp.z.cm:gndr.c+log_gdp.z.cm:essround.c+log_gdp.z.cm:gndr.c:essround.c+
               (gndr.c+essround.c+gndr.c:essround.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
```

```
## Warning in checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv, : Model failed to converge
## with max|grad| = 0.00264536 (tol = 0.002, component 1)
```

``` r
summary(mod6_log_GDP)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ben.z ~ gndr.c + essround.c + gndr.c:essround.c + log_gdp.z.cm:gndr.c +  
##     log_gdp.z.cm:essround.c + log_gdp.z.cm:gndr.c:essround.c +  
##     (gndr.c + essround.c + gndr.c:essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1302801.8 1302999.7 -651382.9 1302765.8    441150 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.7743 -0.5184  0.0761  0.6582  5.2731 
## 
## Random effects:
##  Groups   Name              Variance  Std.Dev. Corr             
##  cntry    (Intercept)       0.0818944 0.28617                   
##           gndr.c            0.0047735 0.06909   0.44            
##           essround.c        0.0008702 0.02950   0.07 -0.36      
##           gndr.c:essround.c 0.0001330 0.01153  -0.23 -0.46 -0.09
##  Residual                   0.9880715 0.99402                   
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##                                  Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)                     0.0093235  0.0499325 33.0019876   0.187 0.853021    
## gndr.c                         -0.2059915  0.0126455 29.9599623 -16.290  < 2e-16 ***
## essround.c                      0.0134051  0.0052433 31.5998471   2.557 0.015584 *  
## gndr.c:essround.c              -0.0002627  0.0024424 32.3522730  -0.108 0.915025    
## gndr.c:log_gdp.z.cm            -0.0679745  0.0117856 32.5113447  -5.768 2.01e-06 ***
## essround.c:log_gdp.z.cm         0.0235934  0.0053827 33.1761727   4.383 0.000111 ***
## gndr.c:essround.c:log_gdp.z.cm  0.0041852  0.0025328 36.8911030   1.652 0.106943    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c essrn. gnd.:. g.:_.. e.:_..
## gndr.c       0.420                                   
## essround.c   0.059 -0.341                            
## gndr.c:ssr. -0.196 -0.388 -0.061                     
## gndr.c:l_..  0.000 -0.026 -0.005  0.009              
## essrnd.:_..  0.009 -0.001 -0.021 -0.012 -0.390       
## gndr.:.:_..  0.004  0.011 -0.009 -0.130 -0.323 -0.038
## optimizer (bobyqa) convergence code: 0 (OK)
## Model failed to converge with max|grad| = 0.00264536 (tol = 0.002, component 1)
```

``` r
getFE(mod6_log_GDP,round=2,p.round = 5)
```

```
##                                 Est.   SE    df      t       p    LL    UL
## (Intercept)                     0.01 0.05 33.00   0.19 0.85302 -0.09  0.11
## gndr.c                         -0.21 0.01 29.96 -16.29 0.00000 -0.23 -0.18
## essround.c                      0.01 0.01 31.60   2.56 0.01558  0.00  0.02
## gndr.c:essround.c               0.00 0.00 32.35  -0.11 0.91502 -0.01  0.00
## gndr.c:log_gdp.z.cm            -0.07 0.01 32.51  -5.77 0.00000 -0.09 -0.04
## essround.c:log_gdp.z.cm         0.02 0.01 33.18   4.38 0.00011  0.01  0.03
## gndr.c:essround.c:log_gdp.z.cm  0.00 0.00 36.89   1.65 0.10694  0.00  0.01
```

``` r
getVC(mod6_log_GDP)
```

```
##         grp              var1              var2 sdcor vcov
## 1     cntry       (Intercept)              <NA>  0.29 0.08
## 2     cntry            gndr.c              <NA>  0.07 0.00
## 3     cntry        essround.c              <NA>  0.03 0.00
## 4     cntry gndr.c:essround.c              <NA>  0.01 0.00
## 5     cntry       (Intercept)            gndr.c  0.44 0.01
## 6     cntry       (Intercept)        essround.c  0.07 0.00
## 7     cntry       (Intercept) gndr.c:essround.c -0.23 0.00
## 8     cntry            gndr.c        essround.c -0.36 0.00
## 9     cntry            gndr.c gndr.c:essround.c -0.46 0.00
## 10    cntry        essround.c gndr.c:essround.c -0.09 0.00
## 11 Residual              <NA>              <NA>  0.99 0.99
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_log_GDP,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] 38.45973
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_log_GDP,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] 2.631829
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
##         4.5 -0.0365 0.0615 Inf   -0.1570    0.0840  -0.594  0.5524
##        -4.5  0.0552 0.0595 Inf   -0.0615    0.1718   0.927  0.3540
## 
## log_gdp.z.cm =  0:
##  essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##         4.5  0.0696 0.0565 Inf   -0.0410    0.1803   1.233  0.2175
##        -4.5 -0.0510 0.0540 Inf   -0.1567    0.0547  -0.945  0.3445
## 
## log_gdp.z.cm =  1:
##  essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##         4.5  0.1758 0.0614 Inf    0.0554    0.2962   2.863  0.0042
##        -4.5 -0.1572 0.0588 Inf   -0.2723   -0.0420  -2.675  0.0075
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
##  essround.c4.5 - (essround.c-4.5)  -0.0917 0.0683 Inf   -0.2257    0.0423  -1.342  0.1797
## 
## log_gdp.z.cm =  0:
##  contrast                         estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   0.1206 0.0472 Inf    0.0282    0.2131   2.557  0.0106
## 
## log_gdp.z.cm =  1:
##  contrast                         estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   0.3330 0.0669 Inf    0.2019    0.4641   4.977  <.0001
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
##    -0.5        4.5  0.04249 0.0629 Inf  -0.08073   0.16571   0.676  0.4991
##     0.5        4.5 -0.11554 0.0616 Inf  -0.23636   0.00528  -1.874  0.0609
##    -0.5       -4.5  0.11417 0.0563 Inf   0.00374   0.22461   2.026  0.0427
##     0.5       -4.5 -0.00383 0.0656 Inf  -0.13245   0.12479  -0.058  0.9535
## 
## log_gdp.z.cm =  0:
##  gndr.c essround.c   emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##    -0.5        4.5  0.17323 0.0565 Inf   0.06250   0.28397   3.066  0.0022
##     0.5        4.5 -0.03394 0.0572 Inf  -0.14606   0.07818  -0.593  0.5530
##    -0.5       -4.5  0.05141 0.0505 Inf  -0.04756   0.15037   1.018  0.3087
##     0.5       -4.5 -0.15340 0.0589 Inf  -0.26881  -0.03800  -2.605  0.0092
## 
## log_gdp.z.cm =  1:
##  gndr.c essround.c   emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##    -0.5        4.5  0.30397 0.0628 Inf   0.18096   0.42699   4.843  <.0001
##     0.5        4.5  0.04766 0.0614 Inf  -0.07261   0.16793   0.777  0.4373
##    -0.5       -4.5 -0.01136 0.0554 Inf  -0.12003   0.09731  -0.205  0.8376
##     0.5       -4.5 -0.30298 0.0647 Inf  -0.42969  -0.17627  -4.686  <.0001
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
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5       0.1580 0.0197 Inf    0.1194   0.19670
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  -0.0717 0.0710 Inf   -0.2108   0.06742
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)    0.0463 0.0763 Inf   -0.1033   0.19590
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)     -0.2297 0.0643 Inf   -0.3557  -0.10372
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      -0.1117 0.0698 Inf   -0.2485   0.02510
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)   0.1180 0.0281 Inf    0.0629   0.17310
##  z.ratio p.value
##    8.010  <.0001
##   -1.010  0.3125
##    0.607  0.5439
##   -3.574  0.0004
##   -1.600  0.1095
##    4.198  <.0001
## 
## log_gdp.z.cm =  0:
##  contrast                                                 estimate     SE  df asymp.LCL asymp.UCL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5       0.2072 0.0131 Inf    0.1814   0.23293
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   0.1218 0.0491 Inf    0.0256   0.21807
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)    0.3266 0.0529 Inf    0.2230   0.43024
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)     -0.0853 0.0445 Inf   -0.1726   0.00186
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)       0.1195 0.0478 Inf    0.0258   0.21314
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)   0.2048 0.0197 Inf    0.1662   0.24345
##  z.ratio p.value
##   15.763  <.0001
##    2.481  0.0131
##    6.180  <.0001
##   -1.918  0.0551
##    2.500  0.0124
##   10.389  <.0001
## 
## log_gdp.z.cm =  1:
##  contrast                                                 estimate     SE  df asymp.LCL asymp.UCL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5       0.2563 0.0179 Inf    0.2212   0.29138
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   0.3153 0.0694 Inf    0.1792   0.45144
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)    0.6070 0.0750 Inf    0.4600   0.75395
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.0590 0.0625 Inf   -0.0635   0.18156
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)       0.3506 0.0676 Inf    0.2182   0.48307
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)   0.2916 0.0264 Inf    0.2398   0.34341
##  z.ratio p.value
##   14.325  <.0001
##    4.541  <.0001
##    8.093  <.0001
##    0.944  0.3452
##    5.189  <.0001
##   11.035  <.0001
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
##  diff_ESS10   -0.158 0.0197 Inf    -0.197   -0.1194  -8.010  <.0001
##  diff_ESS1    -0.118 0.0281 Inf    -0.173   -0.0629  -4.198  <.0001
## 
## log_gdp.z.cm =  0:
##  contrast   estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10   -0.207 0.0131 Inf    -0.233   -0.1814 -15.763  <.0001
##  diff_ESS1    -0.205 0.0197 Inf    -0.243   -0.1662 -10.389  <.0001
## 
## log_gdp.z.cm =  1:
##  contrast   estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10   -0.256 0.0179 Inf    -0.291   -0.2212 -14.325  <.0001
##  diff_ESS1    -0.292 0.0264 Inf    -0.343   -0.2398 -11.035  <.0001
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
##  contrast               estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10 - diff_ESS1 -0.04003 0.0337 Inf   -0.1060    0.0259  -1.189  0.2343
## 
## log_gdp.z.cm =  0:
##  contrast               estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10 - diff_ESS1 -0.00236 0.0220 Inf   -0.0454    0.0407  -0.108  0.9144
## 
## log_gdp.z.cm =  1:
##  contrast               estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10 - diff_ESS1  0.03530 0.0295 Inf   -0.0226    0.0932   1.195  0.2322
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

