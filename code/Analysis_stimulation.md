---
title: "Analysis for Examining the Gender Equality Paradox in Values Using stimulation Value"
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
cntry.sti<-diff_dat %>% group_by(cntry,essround) %>%
  summarise(sti.ctm=mean(sti,na.rm=T),
            sti.ctsd=sd(sti,na.rm=T)) %>%
  group_by(cntry) %>%
  summarise(sti.cm=mean(sti.ctm),
            sti.csd=mean(sti.ctsd)) 
```

```
## `summarise()` has grouped output by 'cntry'. You can override using the `.groups` argument.
```

``` r
grand_mean_ben<-mean(cntry.sti$sti.cm)
grand_sd_ben<-mean(cntry.sti$sti.csd)

# standardized
diff_dat$sti.z<-(diff_dat$sti-grand_mean_ben)/grand_sd_ben
hist(diff_dat$sti.z)
```

![](Analysis_stimulation_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

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

# value-based stimulation

cntry_ben_frame<-
  diff_dat %>% group_by(cntry,essround) %>%
  summarise('sti M' = weighted.mean(x=sti.z,w=pspwght),
            'sti SD' = sqrt(wtd.var(sti.z,pspwght))) %>%
  group_by(cntry) %>%
  summarise('sti M' = mean(x=`sti M`),
            'sti SD'= mean(x=`sti SD`))
```

```
## `summarise()` has grouped output by 'cntry'. You can override using the `.groups` argument.
```

``` r
cntry_ben_women_frame<-
  diff_dat %>%
  filter(gndr.c==-0.5) %>%
  group_by(cntry,essround) %>%
  summarise('sti M' = weighted.mean(x=sti.z,w=pspwght),
            'sti SD' = sqrt(wtd.var(sti.z,pspwght))) %>%
  group_by(cntry) %>%
  summarise('sti M Women' = mean(x=`sti M`),
            'sti SD Women'= mean(x=`sti SD`))
```

```
## `summarise()` has grouped output by 'cntry'. You can override using the `.groups` argument.
```

``` r
cntry_ben_men_frame<-
  diff_dat %>%
  filter(gndr.c==0.5) %>%
  group_by(cntry,essround) %>%
  summarise('sti M' = weighted.mean(x=sti.z,w=pspwght),
            'sti SD' = sqrt(wtd.var(sti.z,pspwght))) %>%
  group_by(cntry) %>%
  summarise('sti M Men' = mean(x=`sti M`),
            'sti SD Men'= mean(x=`sti SD`))
```

```
## `summarise()` has grouped output by 'cntry'. You can override using the `.groups` argument.
```

``` r
# link n and sti datasets

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
desc_frame$D<-desc_frame$`sti M Men`-desc_frame$`sti M Women`

desc_frame
```

```
## # A tibble: 33 × 10
##    cntry `n ESS rounds`     n `sti M` `sti SD` `sti M Women` `sti SD Women` `sti M Men` `sti SD Men`
##    <chr>          <dbl> <dbl>   <dbl>    <dbl>         <dbl>          <dbl>       <dbl>        <dbl>
##  1 AT                 6 13077  0.0293    1.01      -0.0694            1.03       0.135         0.979
##  2 BE                10 17313  0.0938    0.935      0.0223            0.938      0.169         0.926
##  3 BG                 6 12641  0.0999    1.06      -0.00645           1.08       0.215         1.03 
##  4 CH                10 16720  0.102     0.957     -0.000949          0.967      0.209         0.934
##  5 CY                 5  5105  0.154     1.09       0.0505            1.11       0.263         1.07 
##  6 CZ                 9 18934 -0.0323    1.000     -0.150             0.996      0.0955        0.986
##  7 DE                 9 25389 -0.154     0.966     -0.258             0.955     -0.0446        0.965
##  8 DK                 8 12198  0.0814    1.05      -0.00297           1.06       0.168         1.04 
##  9 EE                 9 16692 -0.0539    0.972     -0.136             0.976      0.0447        0.957
## 10 ES                 9 16954 -0.0195    1.03      -0.0684            1.03       0.0312        1.03 
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
    `sti M`, `sti SD`,
    `sti M Women`, `sti SD Women`,
    `sti M Men`, `sti SD Men`,
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
##    Country    `n ESS rounds`     n `sti M` `sti SD` `sti M Women` `sti SD Women` `sti M Men` `sti SD Men`
##    <chr>               <dbl> <dbl> <chr>   <chr>    <chr>         <chr>          <chr>       <chr>       
##  1 Austria                 6 13077 0.03    1.01     -0.07         1.03           0.13        0.98        
##  2 Belgium                10 17313 0.09    0.94     0.02          0.94           0.17        0.93        
##  3 Bulgaria                6 12641 0.10    1.06     -0.01         1.08           0.21        1.03        
##  4 Switzerla…             10 16720 0.10    0.96     -0.00         0.97           0.21        0.93        
##  5 Cyprus                  5  5105 0.15    1.09     0.05          1.11           0.26        1.07        
##  6 Czechia                 9 18934 -0.03   1.00     -0.15         1.00           0.10        0.99        
##  7 Germany                 9 25389 -0.15   0.97     -0.26         0.95           -0.04       0.96        
##  8 Denmark                 8 12198 0.08    1.05     -0.00         1.06           0.17        1.04        
##  9 Estonia                 9 16692 -0.05   0.97     -0.14         0.98           0.04        0.96        
## 10 Spain                   9 16954 -0.02   1.03     -0.07         1.03           0.03        1.03        
## 11 Finland                10 18050 0.08    0.97     0.02          1.00           0.15        0.94        
## 12 France                 10 18720 -0.09   0.99     -0.19         0.97           0.01        1.00        
## 13 UK                     10 21456 0.12    1.01     0.04          1.02           0.20        0.99        
## 14 Greece                  5 12464 0.22    1.01     0.12          1.02           0.33        0.98        
## 15 Croatia                 4  6368 -0.34   1.08     -0.48         1.07           -0.19       1.06        
## 16 Hungary                10 16006 0.10    1.00     0.02          0.99           0.19        0.99        
## 17 Ireland                10 20576 0.22    1.04     0.16          1.04           0.27        1.03        
## 18 Israel                  6 13964 0.26    1.08     0.18          1.08           0.35        1.07        
## 19 Iceland                 5  3832 0.12    1.01     0.07          1.02           0.17        0.99        
## 20 Italy                   4  8663 0.03    0.88     -0.04         0.89           0.11        0.87        
## 21 Lithuania               6 11714 -0.15   1.07     -0.29         1.07           0.02        1.04        
## 22 Latvia                  2  2866 0.22    1.01     0.13          1.00           0.33        1.02        
## 23 Montenegro              2  2441 0.14    0.99     0.03          1.00           0.25        0.96        
## 24 Netherlan…             10 18048 0.16    0.91     0.08          0.91           0.24        0.89        
## 25 Norway                 10 15186 -0.06   1.01     -0.16         1.03           0.04        0.98        
## 26 Poland                  9 15314 -0.01   0.97     -0.13         0.97           0.12        0.96        
## 27 Portugal               10 17705 -0.13   0.93     -0.25         0.94           0.00        0.89        
## 28 Russia                  5 12139 -0.07   1.09     -0.18         1.09           0.08        1.06        
## 29 Sweden                  9 14897 -0.05   0.97     -0.12         0.99           0.03        0.95        
## 30 Slovenia               10 13238 0.31    0.94     0.22          0.94           0.40        0.92        
## 31 Slovakia                7 11132 0.02    0.95     -0.08         0.96           0.14        0.93        
## 32 Turkey                  2  4108 0.30    1.02     0.23          1.02           0.36        1.01        
## 33 Ukraine                 5  9454 -0.25   1.05     -0.32         1.04           -0.17       1.04        
## # ℹ 5 more variables: D <chr>, GEI <chr>, GGGI <chr>, GDI <chr>, GDP <chr>
```

``` r
export(cntry_desc_tbl,"../results/sti/cntry_desc_tbl.xlsx",overwrite=T)
```

## Country-level correlations table


``` r
cor_frame<-
  desc_frame %>%
  select(
    VBMT=`sti M`,
    VBMT_Women=`sti M Women`,
    VBMT_Men=`sti M Men`,
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
  filename = "../results/sti/CorTable1.doc",
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
##   Variable      M     SD   1            2            3           4           5           6          
##   1. VBMT       0.04  0.15                                                                          
##                                                                                                     
##   2. VBMT_Women -0.05 0.16 .99                                                                      
##                            [.98, 1.00]                                                              
##                                                                                                     
##   3. VBMT_Men   0.14  0.14 .98          .95                                                         
##                            [.97, .99]   [.91, .98]                                                  
##                                                                                                     
##   4. D          0.19  0.05 -.44         -.55         -.27                                           
##                            [-.68, -.11] [-.75, -.26] [-.56, .08]                                    
##                                                                                                     
##   5. GEI        0.87  0.07 -.23         -.19         -.30        -.18                               
##                            [-.54, .13]  [-.50, .17]  [-.59, .06] [-.50, .18]                        
##                                                                                                     
##   6. GGGI       0.73  0.05 -.09         -.03         -.16        -.32        .73                    
##                            [-.42, .27]  [-.37, .31]  [-.47, .20] [-.60, .02] [.52, .86]             
##                                                                                                     
##   7. GDI        0.99  0.03 -.39         -.40         -.33        .36         .07         .20        
##                            [-.65, -.05] [-.66, -.07] [-.61, .01] [.02, .63]  [-.29, .41] [-.16, .51]
##                                                                                                     
##   8. log_GDP    10.62 0.40 .10          .14          .04         -.32        .75         .67        
##                            [-.25, .43]  [-.22, .46]  [-.31, .38] [-.60, .03] [.55, .87]  [.42, .82] 
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
mod0<-lmer(sti.z~(1|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod0)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: sti.z ~ (1 | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##   1308428   1308461   -654211   1308422    441165 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -4.7821 -0.7005 -0.0350  0.6509  4.8637 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  cntry    (Intercept) 0.02334  0.1528  
##  Residual             1.00148  1.0007  
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##             Estimate Std. Error       df t value Pr(>|t|)  
## (Intercept)  0.04624    0.02666 32.92490   1.735   0.0921 .
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

``` r
getVC(mod0)
```

```
##        grp        var1 var2 sdcor vcov
## 1    cntry (Intercept) <NA>  0.15 0.02
## 2 Residual        <NA> <NA>  1.00 1.00
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
## mean variation  0.02277899     NA       1
## sigma2          0.97722101      1      NA
## 
## $R2s
##          total within between
## f1  0.00000000      0      NA
## f2  0.00000000     NA       0
## v   0.00000000      0      NA
## m   0.02277899     NA       1
## f   0.00000000     NA      NA
## fv  0.00000000      0      NA
## fvm 0.02277899     NA      NA
```

## mod1: Gender fixed effect


``` r
mod1<-lmer(sti.z~gndr.c+(1|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod1)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: sti.z ~ gndr.c + (1 | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1304382.3 1304426.3 -652187.2 1304374.3    441164 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -4.7890 -0.6926 -0.0407  0.6427  4.9456 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  cntry    (Intercept) 0.0230   0.1516  
##  Residual             0.9923   0.9962  
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##              Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept) 4.988e-02  2.646e-02 3.292e+01   1.886   0.0682 .  
## gndr.c      1.910e-01  2.995e-03 4.411e+05  63.768   <2e-16 ***
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
##              Est.    SE         df      t     p     LL    UL
## (Intercept) 0.050 0.026     32.922  1.886 0.068 -0.004 0.104
## gndr.c      0.191 0.003 441137.327 63.768 0.000  0.185 0.197
```

``` r
getVC(mod1)
```

```
##        grp        var1 var2 sdcor vcov
## 1    cntry (Intercept) <NA>  0.15 0.02
## 2 Residual        <NA> <NA>  1.00 0.99
```

``` r
r2mlm(mod1,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.008849217
## slope variation 0.000000000
## mean variation  0.022447515
## sigma2          0.968703268
## 
## $R2s
##           total
## f   0.008849217
## v   0.000000000
## m   0.022447515
## fv  0.008849217
## fvm 0.031296732
```

## mod2: Gender fixed and random effect

* Include random effect correlation by default


``` r
mod2<-lmer(sti.z~gndr.c+(gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: sti.z ~ gndr.c + (gndr.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1304211.3 1304277.3 -652099.7 1304199.3    441162 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -4.8158 -0.6917 -0.0520  0.6446  5.0428 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr 
##  cntry    (Intercept) 0.022849 0.1512        
##           gndr.c      0.002285 0.0478   -0.46
##  Residual             0.991799 0.9959        
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##              Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)  0.050159   0.026373 32.921463   1.902    0.066 .  
## gndr.c       0.191504   0.008981 32.072424  21.323   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##        (Intr)
## gndr.c -0.423
```

``` r
getFE(mod2,round=3)
```

```
##              Est.    SE     df      t     p     LL    UL
## (Intercept) 0.050 0.026 32.921  1.902 0.066 -0.004 0.104
## gndr.c      0.192 0.009 32.072 21.323 0.000  0.173 0.210
```

``` r
getVC(mod2)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.15 0.02
## 2    cntry      gndr.c   <NA>  0.05 0.00
## 3    cntry (Intercept) gndr.c -0.46 0.00
## 4 Residual        <NA>   <NA>  1.00 0.99
```

``` r
r2mlm(mod2,bargraph = F)
```

```
## $Decompositions
##                        total
## fixed           0.0088943087
## slope variation 0.0005542305
## mean variation  0.0225554159
## sigma2          0.9679960449
## 
## $R2s
##            total
## f   0.0088943087
## v   0.0005542305
## m   0.0225554159
## fv  0.0094485392
## fvm 0.0320039551
```

``` r
anova(mod1,mod2)
```

```
## Data: diff_dat
## Models:
## mod1: sti.z ~ gndr.c + (1 | cntry)
## mod2: sti.z ~ gndr.c + (gndr.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod1    4 1304382 1304426 -652187   1304374                         
## mod2    6 1304211 1304277 -652100   1304199 174.99  2  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

``` r
lvl2_var_cond_lvl1(mod2,lvl1.var = "gndr.c",lvl1.values = c(0.5,-0.5))
```

```
##   lvl1.value lvl2.cond.var lvl2.cond.sd
## 1        0.5    0.02010622    0.1417964
## 2       -0.5    0.02673421    0.1635060
```

* Test for random effect correlation


``` r
mod2_norecov<-lmer(sti.z~gndr.c+(gndr.c||cntry),data=diff_dat,REML=F,weights = pspwght,
                   control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2_norecov)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: sti.z ~ gndr.c + (gndr.c || cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1304215.6 1304270.6 -652102.8 1304205.6    441163 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -4.8170 -0.6916 -0.0515  0.6450  5.0414 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  cntry    (Intercept) 0.022887 0.15128 
##  cntry.1  gndr.c      0.002276 0.04771 
##  Residual             0.991799 0.99589 
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##              Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)  0.050187   0.026395 32.919662   1.901    0.066 .  
## gndr.c       0.192222   0.008976 32.033347  21.415   <2e-16 ***
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
## (Intercept) 0.050 0.026 32.920  1.901 0.066 -0.004 0.104
## gndr.c      0.192 0.009 32.033 21.415 0.000  0.174 0.211
```

``` r
getVC(mod2_norecov)
```

```
##        grp        var1 var2 sdcor vcov
## 1    cntry (Intercept) <NA>  0.15 0.02
## 2  cntry.1      gndr.c <NA>  0.05 0.00
## 3 Residual        <NA> <NA>  1.00 0.99
```

``` r
anova(mod2_norecov,mod2)
```

```
## Data: diff_dat
## Models:
## mod2_norecov: sti.z ~ gndr.c + (gndr.c || cntry)
## mod2: sti.z ~ gndr.c + (gndr.c | cntry)
##              npar     AIC     BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)  
## mod2_norecov    5 1304216 1304271 -652103   1304206                       
## mod2            6 1304211 1304277 -652100   1304199 6.2327  1    0.01254 *
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```


## mod2 with Gender-equality index (GEI)


``` r
mod2_GEI<-lmer(sti.z~gndr.c+gei.z.cm+gndr.c:gei.z.cm+
                 (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2_GEI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: sti.z ~ gndr.c + gei.z.cm + gndr.c:gei.z.cm + (gndr.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##   1273364   1273452   -636674   1273348    431770 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -4.8209 -0.6933 -0.0506  0.6466  5.0488 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr 
##  cntry    (Intercept) 0.019284 0.1389        
##           gndr.c      0.002209 0.0470   -0.62
##  Residual             0.989319 0.9946        
## Number of obs: 431778, groups:  cntry, 32
## 
## Fixed effects:
##                  Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)      0.059322   0.024614 31.966689   2.410   0.0219 *  
## gndr.c           0.193153   0.008983 30.916566  21.501   <2e-16 ***
## gei.z.cm        -0.040317   0.025036 32.109045  -1.610   0.1171    
## gndr.c:gei.z.cm -0.009509   0.009402 34.488058  -1.011   0.3189    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c g.z.cm
## gndr.c      -0.568              
## gei.z.cm    -0.002  0.000       
## gndr.c:g.z.  0.000 -0.043 -0.551
```

``` r
getFE(mod2_GEI,round=3)
```

```
##                   Est.    SE     df      t     p     LL    UL
## (Intercept)      0.059 0.025 31.967  2.410 0.022  0.009 0.109
## gndr.c           0.193 0.009 30.917 21.501 0.000  0.175 0.211
## gei.z.cm        -0.040 0.025 32.109 -1.610 0.117 -0.091 0.011
## gndr.c:gei.z.cm -0.010 0.009 34.488 -1.011 0.319 -0.029 0.010
```

``` r
getVC(mod2_GEI)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.14 0.02
## 2    cntry      gndr.c   <NA>  0.05 0.00
## 3    cntry (Intercept) gndr.c -0.62 0.00
## 4 Residual        <NA>   <NA>  0.99 0.99
```

``` r
r2mlm(mod2_GEI,bargraph = F)
```

```
## $Decompositions
##                        total
## fixed           0.0097944122
## slope variation 0.0005387156
## mean variation  0.0192109155
## sigma2          0.9704559568
## 
## $R2s
##            total
## f   0.0097944122
## v   0.0005387156
## m   0.0192109155
## fv  0.0103331277
## fvm 0.0295440432
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
## Time difference of 30.17719 secs
```

``` r
round(ddsc_mod2_GEI$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.027        0.164        0.992     1.019 0.026   7204.152 0.995   0.995
## 2        0.5         0.020        0.142        0.992     1.012 0.020   6164.576 0.992   0.992
```

``` r
round(ddsc_mod2_GEI$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gei.z.cm
## means_y1           0.156 0.136    1.000           1.000    0.950           0.950   -0.330
## means_y1_scaled    1.044 0.912    1.000           1.000    0.950           0.950   -0.330
## means_y2          -0.037 0.162    0.950           0.950    1.000           1.000   -0.223
## means_y2_scaled   -0.248 1.081    0.950           0.950    1.000           1.000   -0.223
## gei.z.cm           0.000 1.000   -0.330          -0.330   -0.223          -0.223    1.000
## gei.z.cm_scaled    0.000 1.000   -0.330          -0.330   -0.223          -0.223    1.000
## diff_score         0.193 0.053   -0.323          -0.323   -0.602          -0.602   -0.168
## diff_score_scaled  1.292 0.357   -0.323          -0.323   -0.602          -0.602   -0.168
##                   gei.z.cm_scaled diff_score diff_score_scaled
## means_y1                   -0.330     -0.323            -0.323
## means_y1_scaled            -0.330     -0.323            -0.323
## means_y2                   -0.223     -0.602            -0.602
## means_y2_scaled            -0.223     -0.602            -0.602
## gei.z.cm                    1.000     -0.168            -0.168
## gei.z.cm_scaled             1.000     -0.168            -0.168
## diff_score                 -0.168      1.000             1.000
## diff_score_scaled          -0.168      1.000             1.000
```

``` r
round(ddsc_mod2_GEI$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                       0.178 0.176 34.488   1.011   0.319   -0.179    0.536
## w_11                         -0.036 0.028 32.171  -1.275   0.212   -0.092    0.021
## w_21                         -0.045 0.023 32.232  -1.978   0.057   -0.091    0.001
## r_xy1                        -0.261 0.204 32.171  -1.275   0.212   -0.677    0.156
## r_xy2                        -0.279 0.141 32.232  -1.978   0.057   -0.565    0.008
## b_11                         -0.238 0.187 32.171  -1.275   0.212   -0.620    0.143
## b_21                         -0.302 0.153 32.232  -1.978   0.057   -0.613    0.009
## main_effect                  -0.040 0.025 32.109  -1.610   0.117   -0.091    0.011
## moderator_effect              0.193 0.009 30.917  21.501   0.000    0.175    0.211
## interaction                  -0.010 0.009 34.488  -1.011   0.319   -0.029    0.010
## q_b11_b21                     0.069    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.019    NA     NA      NA      NA       NA       NA
## cross_over_point             20.313    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.031 0.031 32.317  -0.987   0.331   -0.094    0.033
## interaction_vs_main_bscale   -0.207 0.209 32.317  -0.987   0.331   -0.633    0.220
## interaction_vs_main_rscale   -0.252 0.241 32.282  -1.045   0.304   -0.742    0.239
## dadas                        -0.071 0.056 32.171  -1.275   0.894   -0.185    0.043
## dadas_bscale                 -0.477 0.374 32.171  -1.275   0.894   -1.239    0.285
## dadas_rscale                 -0.521 0.409 32.171  -1.275   0.894   -1.354    0.312
## abs_diff                      0.010 0.009 34.488   1.011   0.159   -0.010    0.029
## abs_sum                       0.081 0.050 32.109   1.610   0.059   -0.021    0.183
## abs_diff_bscale               0.064 0.063 34.488   1.011   0.159   -0.064    0.192
## abs_sum_bscale                0.541 0.336 32.109   1.610   0.059   -0.143    1.225
## abs_diff_rscale               0.018 0.083 33.473   0.216   0.415   -0.151    0.186
## abs_sum_rscale                0.539 0.341 32.109   1.580   0.062   -0.156    1.234
```

``` r
round(ddsc_mod2_GEI$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.003 -0.459  6.233  1.000  0.013
```

``` r
d_GEI<-ddsc_mod2_GEI$ddsc_sem_fit$data

ddsc_sem_GEI<-
  ddsc_sem(data=d_GEI,x = "gei.z.cm",y1="means_y2",
           y2="means_y1",q_sesoi = .10)
round(ddsc_sem_GEI$results,3)
```

```
##                                    est     se       z pvalue ci.lower ci.upper
## r_xy1_y2                         0.168  0.174   0.965  0.335   -0.173    0.510
## r_xy1                           -0.223  0.172  -1.293  0.196   -0.561    0.115
## r_xy2                           -0.330  0.167  -1.978  0.048   -0.657   -0.003
## b_11                            -0.241  0.186  -1.293  0.196   -0.606    0.124
## b_21                            -0.301  0.152  -1.978  0.048   -0.599   -0.003
## b_10                            -0.248  0.183  -1.353  0.176   -0.608    0.111
## b_20                             1.044  0.150   6.971  0.000    0.750    1.337
## res_cov_y1_y2                    0.837  0.215   3.901  0.000    0.416    1.257
## diff_b10_b20                    -1.292  0.061 -21.105  0.000   -1.412   -1.172
## diff_b11_b21                     0.060  0.062   0.965  0.335   -0.062    0.182
## diff_rxy1_rxy2                   0.107  0.053   2.037  0.042    0.004    0.210
## q_b11_b21                        0.065  0.064   1.014  0.311   -0.060    0.190
## q_rxy1_rxy2                      0.116  0.057   2.032  0.042    0.004    0.228
## cross_over_point                21.531 22.342   0.964  0.335  -22.258   65.321
## sum_b11_b21                     -0.542  0.334  -1.620  0.105   -1.197    0.114
## main_effect                     -0.271  0.167  -1.620  0.105   -0.599    0.057
## interaction_vs_main_effect      -0.211  0.208  -1.012  0.311   -0.619    0.197
## diff_abs_b11_abs_b21            -0.060  0.062  -0.965  0.335   -0.182    0.062
## abs_diff_b11_b21                 0.060  0.062   0.965  0.167   -0.062    0.182
## abs_sum_b11_b21                  0.542  0.334   1.620  0.053   -0.114    1.197
## dadas                           -0.482  0.373  -1.293  0.902   -1.212    0.249
## q_r_equivalence                  0.016  0.057   0.284  0.612       NA       NA
## q_b_equivalence                 -0.035  0.064  -0.551  0.291       NA       NA
## cross_over_point_equivalence    21.531 22.342   0.964  0.832       NA       NA
## cross_over_point_minimal_effect 21.531 22.342   0.964  0.168       NA       NA
```

``` r
round(ddsc_sem_GEI$variance_test,3)
```

```
##             est    se      z pvalue ci.lower ci.upper
## cov_y1y2  0.907 0.233  3.896  0.000    0.451    1.363
## var_y1    1.132 0.283  4.000  0.000    0.577    1.687
## var_y2    0.805 0.201  4.000  0.000    0.411    1.200
## var_diff  0.327 0.133  2.450  0.014    0.065    0.589
## var_ratio 1.406 0.155  9.051  0.000    1.102    1.711
## cor_y1y2  0.950 0.017 55.031  0.000    0.916    0.984
```

### Figure 


``` r
# plot the results

# start with obtaining predicted values for means and differences

# refit reduced and full models with GGGI in original scale


mod2_GEI_unstd<-lmer(sti.z~gndr.c+gei.cm+gndr.c:gei.cm+
                 (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

mod2_GEI_unstd_red<-lmer(sti.z~gndr.c+
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


p1.sti.flags<-
  ggplot(p,aes(y=yvar,x=gei.cm,color=gndr.c))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("Gender Equality Index (Average 2002-2022)")+
  #ylim(c(min.y.pred,max.y.pred))+
  #xlim(c(0.60,1.00))+
  ylab("Value stimulation (Average 2002-2022)")+
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

p2.sti.flags<-ggplot(p2,aes(y=yvar,x=gei.cm))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("Gender Equality Index (Average 2002-2022)")+
  #xlim(c(0.60,1.00))+
  #ylim(c(min.y.narrow,max.y.narrow))+
  ylab("Difference in value stimulation")+
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
  ggarrange(p1.sti.flags,p2.sti.flags,align = "v",
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

![](Analysis_stimulation_files/figure-html/unnamed-chunk-18-1.png)<!-- -->

``` r
png(filename = 
      "../results/sti/GEI_flags.png",
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
mod2_GGGI<-lmer(sti.z~gndr.c+gggi.z.cm+gndr.c:gggi.z.cm+
                 (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2_GGGI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: sti.z ~ gndr.c + gggi.z.cm + gndr.c:gggi.z.cm + (gndr.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##  931198.2  931283.5 -465591.1  931182.2    314638 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -4.9136 -0.6919 -0.0435  0.6477  5.0479 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr 
##  cntry    (Intercept) 0.027993 0.1673        
##           gndr.c      0.002107 0.0459   -0.44
##  Residual             0.993343 0.9967        
## Number of obs: 314646, groups:  cntry, 33
## 
## Fixed effects:
##                   Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)       0.075475   0.029206 32.839111   2.584   0.0144 *  
## gndr.c            0.187965   0.008936 30.700164  21.036   <2e-16 ***
## gggi.z.cm        -0.029137   0.029693 32.992322  -0.981   0.3336    
## gndr.c:gggi.z.cm -0.008035   0.009447 35.036666  -0.851   0.4008    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c ggg.z.
## gndr.c      -0.395              
## gggi.z.cm   -0.001  0.000       
## gndr.c:gg..  0.000 -0.026 -0.380
```

``` r
getFE(mod2_GGGI,round=3)
```

```
##                    Est.    SE     df      t     p     LL    UL
## (Intercept)       0.075 0.029 32.839  2.584 0.014  0.016 0.135
## gndr.c            0.188 0.009 30.700 21.036 0.000  0.170 0.206
## gggi.z.cm        -0.029 0.030 32.992 -0.981 0.334 -0.090 0.031
## gndr.c:gggi.z.cm -0.008 0.009 35.037 -0.851 0.401 -0.027 0.011
```

``` r
getVC(mod2_GGGI)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.17 0.03
## 2    cntry      gndr.c   <NA>  0.05 0.00
## 3    cntry (Intercept) gndr.c -0.44 0.00
## 4 Residual        <NA>   <NA>  1.00 0.99
```

``` r
r2mlm(mod2_GGGI,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.008880218
## slope variation 0.000507422
## mean variation  0.027411636
## sigma2          0.963200724
## 
## $R2s
##           total
## f   0.008880218
## v   0.000507422
## m   0.027411636
## fv  0.009387640
## fvm 0.036799276
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
## Time difference of 29.46783 secs
```

``` r
round(ddsc_mod2_GGGI$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.027        0.164        0.992     1.019 0.026   7204.152 0.995   0.995
## 2        0.5         0.020        0.142        0.992     1.012 0.020   6164.576 0.992   0.992
```

``` r
round(ddsc_mod2_GGGI$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gggi.z.cm
## means_y1           0.170 0.166    1.000           1.000    0.958           0.958    -0.194
## means_y1_scaled    0.966 0.944    1.000           1.000    0.958           0.958    -0.194
## means_y2          -0.018 0.185    0.958           0.958    1.000           1.000    -0.142
## means_y2_scaled   -0.101 1.053    0.958           0.958    1.000           1.000    -0.142
## gggi.z.cm          0.000 1.000   -0.194          -0.194   -0.142          -0.142     1.000
## gggi.z.cm_scaled   0.000 1.000   -0.194          -0.194   -0.142          -0.142     1.000
## diff_score         0.187 0.054   -0.208          -0.208   -0.479          -0.479    -0.112
## diff_score_scaled  1.067 0.308   -0.208          -0.208   -0.479          -0.479    -0.112
##                   gggi.z.cm_scaled diff_score diff_score_scaled
## means_y1                    -0.194     -0.208            -0.208
## means_y1_scaled             -0.194     -0.208            -0.208
## means_y2                    -0.142     -0.479            -0.479
## means_y2_scaled             -0.142     -0.479            -0.479
## gggi.z.cm                    1.000     -0.112            -0.112
## gggi.z.cm_scaled             1.000     -0.112            -0.112
## diff_score                  -0.112      1.000             1.000
## diff_score_scaled           -0.112      1.000             1.000
```

``` r
round(ddsc_mod2_GGGI$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                       0.148 0.175 35.037   0.851   0.401   -0.206    0.503
## w_11                         -0.025 0.032 33.007  -0.790   0.435   -0.090    0.040
## w_21                         -0.033 0.028 32.943  -1.174   0.249   -0.091    0.024
## r_xy1                        -0.152 0.192 33.007  -0.790   0.435   -0.542    0.239
## r_xy2                        -0.179 0.153 32.943  -1.174   0.249   -0.490    0.132
## b_11                         -0.143 0.181 33.007  -0.790   0.435   -0.512    0.226
## b_21                         -0.189 0.161 32.943  -1.174   0.249   -0.517    0.139
## main_effect                  -0.029 0.030 32.992  -0.981   0.334   -0.090    0.031
## moderator_effect              0.188 0.009 30.700  21.036   0.000    0.170    0.206
## interaction                  -0.008 0.009 35.037  -0.851   0.401   -0.027    0.011
## q_b11_b21                     0.047    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.029    NA     NA      NA      NA       NA       NA
## cross_over_point             23.393    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.021 0.034 33.007  -0.613   0.544   -0.091    0.049
## interaction_vs_main_bscale   -0.120 0.196 33.007  -0.613   0.544   -0.520    0.279
## interaction_vs_main_rscale   -0.138 0.216 33.007  -0.638   0.528   -0.576    0.301
## dadas                        -0.050 0.064 33.007  -0.790   0.782   -0.180    0.079
## dadas_bscale                 -0.287 0.363 33.007  -0.790   0.782   -1.025    0.451
## dadas_rscale                 -0.303 0.384 33.007  -0.790   0.782   -1.083    0.477
## abs_diff                      0.008 0.009 35.037   0.851   0.200   -0.011    0.027
## abs_sum                       0.058 0.059 32.992   0.981   0.167   -0.063    0.179
## abs_diff_bscale               0.046 0.054 35.037   0.851   0.200   -0.064    0.155
## abs_sum_bscale                0.332 0.339 32.992   0.981   0.167   -0.357    1.022
## abs_diff_rscale               0.028 0.063 33.955   0.440   0.331   -0.101    0.157
## abs_sum_rscale                0.331 0.341 32.994   0.971   0.169   -0.363    1.025
```

``` r
round(ddsc_mod2_GGGI$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.003 -0.459  6.233  1.000  0.013
```

``` r
d_GGGI<-ddsc_mod2_GGGI$ddsc_sem_fit$data

ddsc_sem_GGGI<-
  ddsc_sem(data=d_GGGI,x = "gggi.z.cm",y1="means_y2",
           y2="means_y1",q_sesoi = .10)
round(ddsc_sem_GGGI$results,3)
```

```
##                                    est     se       z pvalue ci.lower ci.upper
## r_xy1_y2                         0.112  0.173   0.650  0.516   -0.227    0.451
## r_xy1                           -0.142  0.172  -0.821  0.411   -0.479    0.196
## r_xy2                           -0.194  0.171  -1.139  0.255   -0.529    0.140
## b_11                            -0.149  0.181  -0.821  0.411   -0.505    0.207
## b_21                            -0.184  0.161  -1.139  0.255   -0.500    0.132
## b_10                            -0.101  0.179  -0.567  0.571   -0.451    0.249
## b_20                             0.966  0.159   6.082  0.000    0.655    1.277
## res_cov_y1_y2                    0.897  0.226   3.975  0.000    0.455    1.339
## diff_b10_b20                    -1.067  0.053 -20.323  0.000   -1.170   -0.964
## diff_b11_b21                     0.035  0.053   0.650  0.516   -0.070    0.139
## diff_rxy1_rxy2                   0.053  0.050   1.068  0.286   -0.044    0.150
## q_b11_b21                        0.036  0.054   0.659  0.510   -0.070    0.142
## q_rxy1_rxy2                      0.054  0.051   1.067  0.286   -0.046    0.154
## cross_over_point                30.803 47.435   0.649  0.516  -62.167  123.773
## sum_b11_b21                     -0.333  0.339  -0.981  0.327   -0.997    0.332
## main_effect                     -0.166  0.170  -0.981  0.327   -0.499    0.166
## interaction_vs_main_effect      -0.132  0.196  -0.671  0.502   -0.516    0.253
## diff_abs_b11_abs_b21            -0.035  0.053  -0.650  0.516   -0.139    0.070
## abs_diff_b11_b21                 0.035  0.053   0.650  0.258   -0.070    0.139
## abs_sum_b11_b21                  0.333  0.339   0.981  0.163   -0.332    0.997
## dadas                           -0.298  0.363  -0.821  0.794   -1.009    0.413
## q_r_equivalence                 -0.046  0.051  -0.893  0.186       NA       NA
## q_b_equivalence                 -0.064  0.054  -1.190  0.117       NA       NA
## cross_over_point_equivalence    30.803 47.435   0.649  0.742       NA       NA
## cross_over_point_minimal_effect 30.803 47.435   0.649  0.258       NA       NA
```

``` r
round(ddsc_sem_GGGI$variance_test,3)
```

```
##             est    se      z pvalue ci.lower ci.upper
## cov_y1y2  0.924 0.232  3.974  0.000    0.468    1.379
## var_y1    1.074 0.265  4.062  0.000    0.556    1.593
## var_y2    0.865 0.213  4.062  0.000    0.448    1.282
## var_diff  0.210 0.109  1.920  0.055   -0.004    0.423
## var_ratio 1.242 0.124 10.027  0.000    0.999    1.485
## cor_y1y2  0.958 0.014 67.073  0.000    0.930    0.986
```

### Figure


``` r
# plot the results

# start with obtaining predicted values for means and differences

# refit reduced and full models with GGGI in original scale


mod2_GGGI_unstd<-lmer(sti.z~gndr.c+gggi.cm+gndr.c:gggi.cm+
                       (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
                     control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

mod2_GGGI_unstd_red<-lmer(sti.z~gndr.c+
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


p1.sti.flags<-
  ggplot(p,aes(y=yvar,x=gggi.cm,color=gndr.c))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("Global Gender Gap Index (Average 2002-2022)")+
  #ylim(c(min.y.pred,max.y.pred))+
  #xlim(c(0.60,1.00))+
  ylab("Value stimulation (Average 2002-2022)")+
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

p2.sti.flags<-ggplot(p2,aes(y=yvar,x=gggi.cm))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("Global Gender Gap Index (Average 2002-2022)")+
  #xlim(c(0.60,1.00))+
  #ylim(c(min.y.narrow,max.y.narrow))+
  ylab("Difference in value stimulation")+
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
  ggarrange(p1.sti.flags,p2.sti.flags,align = "v",
            ncol=1,nrow=2,heights=c(2,1))
```

```
## Warning in geom_text(inherit.aes = F, aes(x = 0.61, y = -0.45, label = coef_q, : All aesthetics have length 1, but the data has 502 rows.
## ℹ Please consider using `annotate()` or provide this layer with data containing a single row.
```

``` r
pflag_comb
```

![](Analysis_stimulation_files/figure-html/unnamed-chunk-21-1.png)<!-- -->

``` r
png(filename = 
      "../results/sti/GGGI_flags_new.png",
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
mod2_GDI<-lmer(sti.z~gndr.c+gdi.z.cm+gndr.c:gdi.z.cm+
                 (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2_GDI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: sti.z ~ gndr.c + gdi.z.cm + gndr.c:gdi.z.cm + (gndr.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1304208.1 1304296.1 -652096.1 1304192.1    441160 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -4.8180 -0.6916 -0.0518  0.6447  5.0449 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr 
##  cntry    (Intercept) 0.019955 0.14126       
##           gndr.c      0.001887 0.04344  -0.36
##  Residual             0.991798 0.99589       
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##                  Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)      0.050197   0.024655 32.878009   2.036   0.0499 *  
## gndr.c           0.191435   0.008278 32.400574  23.126   <2e-16 ***
## gdi.z.cm        -0.054668   0.025075 33.077434  -2.180   0.0365 *  
## gndr.c:gdi.z.cm  0.021234   0.008794 37.919962   2.415   0.0207 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c gd.z.c
## gndr.c      -0.324              
## gdi.z.cm     0.000  0.001       
## gndr.c:gd..  0.001 -0.013 -0.309
```

``` r
getFE(mod2_GDI,round=3)
```

```
##                   Est.    SE     df      t     p     LL     UL
## (Intercept)      0.050 0.025 32.878  2.036 0.050  0.000  0.100
## gndr.c           0.191 0.008 32.401 23.126 0.000  0.175  0.208
## gdi.z.cm        -0.055 0.025 33.077 -2.180 0.036 -0.106 -0.004
## gndr.c:gdi.z.cm  0.021 0.009 37.920  2.415 0.021  0.003  0.039
```

``` r
getVC(mod2_GDI)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.14 0.02
## 2    cntry      gndr.c   <NA>  0.04 0.00
## 3    cntry (Intercept) gndr.c -0.36 0.00
## 4 Residual        <NA>   <NA>  1.00 0.99
```

``` r
r2mlm(mod2_GDI,bargraph = F)
```

```
## $Decompositions
##                        total
## fixed           0.0114071992
## slope variation 0.0004578007
## mean variation  0.0196550521
## sigma2          0.9684799479
## 
## $R2s
##            total
## f   0.0114071992
## v   0.0004578007
## m   0.0196550521
## fv  0.0118650000
## fvm 0.0315200521
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
## Time difference of 30.29069 secs
```

``` r
round(ddsc_mod2_GDI$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.027        0.164        0.992     1.019 0.026   7204.152 0.995   0.995
## 2        0.5         0.020        0.142        0.992     1.012 0.020   6164.576 0.992   0.992
```

``` r
round(ddsc_mod2_GDI$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gdi.z.cm
## means_y1           0.146 0.145    1.000           1.000    0.951           0.951   -0.301
## means_y1_scaled    0.937 0.930    1.000           1.000    0.951           0.951   -0.301
## means_y2          -0.046 0.167    0.951           0.951    1.000           1.000   -0.391
## means_y2_scaled   -0.292 1.065    0.951           0.951    1.000           1.000   -0.391
## gdi.z.cm           0.000 1.000   -0.301          -0.301   -0.391          -0.391    1.000
## gdi.z.cm_scaled    0.000 1.000   -0.301          -0.301   -0.391          -0.391    1.000
## diff_score         0.192 0.053   -0.245          -0.245   -0.532          -0.532    0.402
## diff_score_scaled  1.229 0.339   -0.245          -0.245   -0.532          -0.532    0.402
##                   gdi.z.cm_scaled diff_score diff_score_scaled
## means_y1                   -0.301     -0.245            -0.245
## means_y1_scaled            -0.301     -0.245            -0.245
## means_y2                   -0.391     -0.532            -0.532
## means_y2_scaled            -0.391     -0.532            -0.532
## gdi.z.cm                    1.000      0.402             0.402
## gdi.z.cm_scaled             1.000      0.402             0.402
## diff_score                  0.402      1.000             1.000
## diff_score_scaled           0.402      1.000             1.000
```

``` r
round(ddsc_mod2_GDI$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                      -0.400 0.166 37.920  -2.415   0.021   -0.736   -0.065
## w_11                         -0.065 0.027 33.151  -2.440   0.020   -0.120   -0.011
## w_21                         -0.044 0.024 33.135  -1.829   0.076   -0.093    0.005
## r_xy1                        -0.449 0.184 33.151  -2.440   0.020   -0.823   -0.075
## r_xy2                        -0.264 0.145 33.135  -1.829   0.076   -0.558    0.030
## b_11                         -0.418 0.172 33.151  -2.440   0.020   -0.767   -0.070
## b_21                         -0.282 0.154 33.135  -1.829   0.076   -0.596    0.032
## main_effect                  -0.055 0.025 33.077  -2.180   0.036   -0.106   -0.004
## moderator_effect              0.191 0.008 32.401  23.126   0.000    0.175    0.208
## interaction                   0.021 0.009 37.920   2.415   0.021    0.003    0.039
## q_b11_b21                    -0.156    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                  -0.212    NA     NA      NA      NA       NA       NA
## cross_over_point             -9.016    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.033 0.024 33.326  -1.400   0.171   -0.082    0.015
## interaction_vs_main_bscale   -0.214 0.153 33.326  -1.400   0.171   -0.525    0.097
## interaction_vs_main_rscale   -0.172 0.133 33.392  -1.292   0.205   -0.443    0.099
## dadas                        -0.088 0.048 33.135  -1.829   0.962   -0.186    0.010
## dadas_bscale                 -0.565 0.309 33.135  -1.829   0.962   -1.193    0.063
## dadas_rscale                 -0.529 0.289 33.135  -1.829   0.962   -1.117    0.059
## abs_diff                      0.021 0.009 37.920   2.415   0.010    0.003    0.039
## abs_sum                       0.109 0.050 33.077   2.180   0.018    0.007    0.211
## abs_diff_bscale               0.136 0.056 37.920   2.415   0.010    0.022    0.250
## abs_sum_bscale                0.701 0.321 33.077   2.180   0.018    0.047    1.355
## abs_diff_rscale               0.184 0.067 36.253   2.765   0.004    0.049    0.320
## abs_sum_rscale                0.713 0.324 33.078   2.201   0.017    0.054    1.373
```

``` r
round(ddsc_mod2_GDI$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.003 -0.459  6.233  1.000  0.013
```

``` r
d_GDI<-ddsc_mod2_GDI$ddsc_sem_fit$data

ddsc_sem_GDI<-
  ddsc_sem(data=d_GDI,x = "gdi.z.cm",y1="means_y2",
           y2="means_y1",q_sesoi = .10)
round(ddsc_sem_GDI$results,3)
```

```
##                                    est    se       z pvalue ci.lower ci.upper
## r_xy1_y2                        -0.402 0.159  -2.523  0.012   -0.715   -0.090
## r_xy1                           -0.391 0.160  -2.441  0.015   -0.705   -0.077
## r_xy2                           -0.301 0.166  -1.816  0.069   -0.627    0.024
## b_11                            -0.417 0.171  -2.441  0.015   -0.751   -0.082
## b_21                            -0.280 0.154  -1.816  0.069   -0.583    0.022
## b_10                            -0.292 0.168  -1.737  0.082   -0.621    0.037
## b_20                             0.937 0.152   6.161  0.000    0.639    1.235
## res_cov_y1_y2                    0.801 0.202   3.956  0.000    0.404    1.197
## diff_b10_b20                    -1.229 0.053 -23.085  0.000   -1.333   -1.124
## diff_b11_b21                    -0.136 0.054  -2.523  0.012   -0.242   -0.030
## diff_rxy1_rxy2                  -0.090 0.052  -1.724  0.085   -0.192    0.012
## q_b11_b21                       -0.156 0.071  -2.200  0.028   -0.294   -0.017
## q_rxy1_rxy2                     -0.102 0.059  -1.721  0.085   -0.219    0.014
## cross_over_point                -9.009 3.592  -2.508  0.012  -16.049   -1.969
## sum_b11_b21                     -0.697 0.321  -2.172  0.030   -1.326   -0.068
## main_effect                     -0.348 0.160  -2.172  0.030   -0.663   -0.034
## interaction_vs_main_effect      -0.212 0.153  -1.387  0.165   -0.512    0.088
## diff_abs_b11_abs_b21             0.136 0.054   2.523  0.012    0.030    0.242
## abs_diff_b11_b21                 0.136 0.054   2.523  0.006    0.030    0.242
## abs_sum_b11_b21                  0.697 0.321   2.172  0.015    0.068    1.326
## dadas                           -0.561 0.309  -1.816  0.965   -1.166    0.045
## q_r_equivalence                  0.002 0.059   0.036  0.514       NA       NA
## q_b_equivalence                  0.056 0.071   0.787  0.784       NA       NA
## cross_over_point_equivalence     9.009 3.592   2.508  0.994       NA       NA
## cross_over_point_minimal_effect  9.009 3.592   2.508  0.006       NA       NA
```

``` r
round(ddsc_sem_GDI$variance_test,3)
```

```
##             est    se      z pvalue ci.lower ci.upper
## cov_y1y2  0.914 0.231  3.959  0.000    0.461    1.366
## var_y1    1.100 0.271  4.062  0.000    0.569    1.631
## var_y2    0.839 0.207  4.062  0.000    0.434    1.244
## var_diff  0.261 0.122  2.149  0.032    0.023    0.500
## var_ratio 1.312 0.141  9.307  0.000    1.035    1.588
## cor_y1y2  0.951 0.017 57.365  0.000    0.919    0.984
```

### Figure


``` r
# plot the results

# start with obtaining predicted values for means and differences


mod2_GDI_unstd<-lmer(sti.z~gndr.c+gdi.cm+gndr.c:gdi.cm+
                       (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
                     control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

mod2_GDI_unstd_red<-lmer(sti.z~gndr.c+
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


p1.sti.flags<-
  ggplot(p,aes(y=yvar,x=gdi.cm,color=gndr.c))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("Gender Development Index (Average 2002-2022)")+
  #ylim(c(min.y.pred,max.y.pred))+
  #xlim(c(0.60,1.00))+
  ylab("Value stimulation (Average 2002-2022)")+
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

#p1.sti.flags


p2.sti.flags<-ggplot(p2,aes(y=yvar,x=gdi.cm))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("Gender Development Index (Average 2002-2022)")+
  #xlim(c(0.60,1.00))+
  #ylim(c(min.y.narrow,max.y.narrow))+
  ylab("Difference in value stimulation")+
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

#p2.sti.flags


pflag_comb<-
  ggarrange(p1.sti.flags,p2.sti.flags,align = "v",
            ncol=1,nrow=2,heights=c(2,1))
```

```
## Warning in geom_text(inherit.aes = F, aes(x = 0.9, y = -0.45, label = coef_q, : All aesthetics have length 1, but the data has 282 rows.
## ℹ Please consider using `annotate()` or provide this layer with data containing a single row.
```

``` r
pflag_comb
```

![](Analysis_stimulation_files/figure-html/unnamed-chunk-24-1.png)<!-- -->

``` r
png(filename = 
      "../results/sti/GDI_flags.png",
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
mod2_log_GDP<-lmer(sti.z~gndr.c+log_gdp.z.cm+
                     gndr.c:log_gdp.z.cm+
                 (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2_log_GDP)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: sti.z ~ gndr.c + log_gdp.z.cm + gndr.c:log_gdp.z.cm + (gndr.c |      cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##   1304212   1304300   -652098   1304196    441160 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -4.8155 -0.6917 -0.0517  0.6446  5.0426 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr 
##  cntry    (Intercept) 0.02277  0.15090       
##           gndr.c      0.00202  0.04494  -0.46
##  Residual             0.99180  0.99589       
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##                      Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)          0.050396   0.026334 32.910546   1.914   0.0644 .  
## gndr.c               0.192016   0.008521 31.706445  22.533   <2e-16 ***
## log_gdp.z.cm         0.008737   0.026442 32.994185   0.330   0.7432    
## gndr.c:log_gdp.z.cm -0.016081   0.008705 33.569866  -1.847   0.0735 .  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c lg_g..
## gndr.c      -0.416              
## lg_gdp.z.cm  0.021 -0.010       
## gndr.c:l_.. -0.010 -0.027 -0.409
```

``` r
getFE(mod2_log_GDP,round=3)
```

```
##                       Est.    SE     df      t     p     LL    UL
## (Intercept)          0.050 0.026 32.911  1.914 0.064 -0.003 0.104
## gndr.c               0.192 0.009 31.706 22.533 0.000  0.175 0.209
## log_gdp.z.cm         0.009 0.026 32.994  0.330 0.743 -0.045 0.063
## gndr.c:log_gdp.z.cm -0.016 0.009 33.570 -1.847 0.074 -0.034 0.002
```

``` r
getVC(mod2_log_GDP)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.15 0.02
## 2    cntry      gndr.c   <NA>  0.04 0.00
## 3    cntry (Intercept) gndr.c -0.46 0.00
## 4 Residual        <NA>   <NA>  1.00 0.99
```

``` r
r2mlm(mod2_log_GDP,bargraph = F)
```

```
## $Decompositions
##                        total
## fixed           0.0088552348
## slope variation 0.0004899747
## mean variation  0.0224663360
## sigma2          0.9681884544
## 
## $R2s
##            total
## f   0.0088552348
## v   0.0004899747
## m   0.0224663360
## fv  0.0093452095
## fvm 0.0318115456
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
## Time difference of 29.71601 secs
```

``` r
round(ddsc_mod2_log_GDP$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.027        0.164        0.992     1.019 0.026   7204.152 0.995   0.995
## 2        0.5         0.020        0.142        0.992     1.012 0.020   6164.576 0.992   0.992
```

``` r
round(ddsc_mod2_log_GDP$descriptives,3)
```

```
##                          M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled log_gdp.z.cm
## means_y1             0.146 0.145    1.000           1.000    0.951           0.951        0.002
## means_y1_scaled      0.937 0.930    1.000           1.000    0.951           0.951        0.002
## means_y2            -0.046 0.167    0.951           0.951    1.000           1.000        0.102
## means_y2_scaled     -0.292 1.065    0.951           0.951    1.000           1.000        0.102
## log_gdp.z.cm        -0.022 1.012    0.002           0.002    0.102           0.102        1.000
## log_gdp.z.cm_scaled  0.000 1.000    0.002           0.002    0.102           0.102        1.000
## diff_score           0.192 0.053   -0.245          -0.245   -0.532          -0.532       -0.314
## diff_score_scaled    1.229 0.339   -0.245          -0.245   -0.532          -0.532       -0.314
##                     log_gdp.z.cm_scaled diff_score diff_score_scaled
## means_y1                          0.002     -0.245            -0.245
## means_y1_scaled                   0.002     -0.245            -0.245
## means_y2                          0.102     -0.532            -0.532
## means_y2_scaled                   0.102     -0.532            -0.532
## log_gdp.z.cm                      1.000     -0.314            -0.314
## log_gdp.z.cm_scaled               1.000     -0.314            -0.314
## diff_score                       -0.314      1.000             1.000
## diff_score_scaled                -0.314      1.000             1.000
```

``` r
round(ddsc_mod2_log_GDP$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                       0.303 0.164 33.570   1.847   0.074   -0.031    0.637
## w_11                          0.017 0.028 32.986   0.589   0.560   -0.041    0.075
## w_21                          0.001 0.025 32.955   0.028   0.978   -0.050    0.052
## r_xy1                         0.115 0.196 32.986   0.589   0.560   -0.283    0.514
## r_xy2                         0.004 0.150 32.955   0.028   0.978   -0.301    0.309
## b_11                          0.108 0.183 32.986   0.589   0.560   -0.264    0.479
## b_21                          0.004 0.160 32.955   0.028   0.978   -0.321    0.330
## main_effect                   0.009 0.026 32.994   0.330   0.743   -0.045    0.063
## moderator_effect              0.192 0.009 31.706  22.533   0.000    0.175    0.209
## interaction                  -0.016 0.009 33.570  -1.847   0.074   -0.034    0.002
## q_b11_b21                     0.103    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.112    NA     NA      NA      NA       NA       NA
## cross_over_point             11.940    NA     NA      NA      NA       NA       NA
## interaction_vs_main           0.007 0.024 32.870   0.303   0.764   -0.042    0.057
## interaction_vs_main_bscale    0.047 0.155 32.870   0.303   0.764   -0.269    0.363
## interaction_vs_main_rscale    0.051 0.134 32.849   0.382   0.705   -0.222    0.325
## dadas                        -0.001 0.050 32.955  -0.028   0.511   -0.103    0.100
## dadas_bscale                 -0.009 0.320 32.955  -0.028   0.511   -0.660    0.643
## dadas_rscale                 -0.008 0.300 32.955  -0.028   0.511   -0.619    0.602
## abs_diff                      0.016 0.009 33.570   1.847   0.037   -0.002    0.034
## abs_sum                       0.017 0.053 32.994   0.330   0.372   -0.090    0.125
## abs_diff_bscale               0.103 0.056 33.570   1.847   0.037   -0.010    0.216
## abs_sum_bscale                0.112 0.339 32.994   0.330   0.372   -0.578    0.802
## abs_diff_rscale               0.111 0.069 33.103   1.616   0.058   -0.029    0.251
## abs_sum_rscale                0.120 0.342 32.995   0.349   0.364   -0.576    0.815
```

``` r
round(ddsc_mod2_log_GDP$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.003 -0.459  6.233  1.000  0.013
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
## r_xy1_y2                         0.314 0.165   1.901  0.057   -0.010    0.638
## r_xy1                            0.102 0.173   0.589  0.556   -0.237    0.441
## r_xy2                            0.002 0.174   0.013  0.990   -0.339    0.343
## b_11                             0.109 0.184   0.589  0.556   -0.253    0.470
## b_21                             0.002 0.162   0.013  0.990   -0.315    0.319
## b_10                            -0.292 0.182  -1.607  0.108   -0.648    0.064
## b_20                             0.937 0.159   5.875  0.000    0.624    1.249
## res_cov_y1_y2                    0.914 0.230   3.970  0.000    0.463    1.365
## diff_b10_b20                    -1.229 0.055 -22.263  0.000   -1.337   -1.121
## diff_b11_b21                     0.107 0.056   1.901  0.057   -0.003    0.216
## diff_rxy1_rxy2                   0.100 0.052   1.934  0.053   -0.001    0.201
## q_b11_b21                        0.107 0.057   1.868  0.062   -0.005    0.219
## q_rxy1_rxy2                      0.100 0.052   1.931  0.053   -0.002    0.202
## cross_over_point                11.534 6.091   1.894  0.058   -0.403   23.472
## sum_b11_b21                      0.111 0.343   0.323  0.746   -0.561    0.782
## main_effect                      0.055 0.171   0.323  0.746   -0.280    0.391
## interaction_vs_main_effect       0.051 0.157   0.326  0.745   -0.257    0.359
## diff_abs_b11_abs_b21             0.107 0.056   1.901  0.057   -0.003    0.216
## abs_diff_b11_b21                 0.107 0.056   1.901  0.029   -0.003    0.216
## abs_sum_b11_b21                  0.111 0.343   0.323  0.373   -0.561    0.782
## dadas                           -0.004 0.324  -0.013  0.505   -0.639    0.630
## q_r_equivalence                  0.000 0.052   0.001  0.500       NA       NA
## q_b_equivalence                  0.007 0.057   0.121  0.548       NA       NA
## cross_over_point_equivalence    11.534 6.091   1.894  0.971       NA       NA
## cross_over_point_minimal_effect 11.534 6.091   1.894  0.029       NA       NA
```

``` r
round(ddsc_sem_log_GDP$variance_test,3)
```

```
##             est    se      z pvalue ci.lower ci.upper
## cov_y1y2  0.914 0.231  3.959  0.000    0.461    1.366
## var_y1    1.100 0.271  4.062  0.000    0.569    1.631
## var_y2    0.839 0.207  4.062  0.000    0.434    1.244
## var_diff  0.261 0.122  2.149  0.032    0.023    0.500
## var_ratio 1.312 0.141  9.307  0.000    1.035    1.588
## cor_y1y2  0.951 0.017 57.365  0.000    0.919    0.984
```

### Figure


``` r
# plot the results

# start with obtaining predicted values for means and differences


mod2_log_GDP_unstd<-lmer(sti.z~gndr.c+log_gdp.cm+
                           gndr.c:log_gdp.cm+
                       (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
                     control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

mod2_log_GDP_unstd_red<-lmer(sti.z~gndr.c+
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


p1.sti.flags<-
  ggplot(p,aes(y=yvar,x=log_gdp.cm,color=gndr.c))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("log(GDP per capita) (Average 2002-2022)")+
  #ylim(c(min.y.pred,max.y.pred))+
  #xlim(c(0.60,1.00))+
  ylab("Value stimulation (Average 2002-2022)")+
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

p2.sti.flags<-ggplot(p2,aes(y=yvar,x=log_gdp.cm))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("log(GDP per capita) (Average 2002-2022)")+
  #xlim(c(0.60,1.00))+
  #ylim(c(min.y.narrow,max.y.narrow))+
  ylab("Difference in value stimulation")+
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
  ggarrange(p1.sti.flags,p2.sti.flags,align = "v",
            ncol=1,nrow=2,heights=c(2,1))
```

```
## Warning in geom_text(inherit.aes = F, aes(x = 9.75, y = -0.45, label = coef_q, : All aesthetics have length 1, but the data has 3302 rows.
## ℹ Please consider using `annotate()` or provide this layer with data containing a single row.
```

``` r
pflag_comb
```

![](Analysis_stimulation_files/figure-html/unnamed-chunk-27-1.png)<!-- -->

``` r
png(filename = 
      "../results/sti/log_GDP_flags.png",
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
mod3<-lmer(sti.z~gndr.c+essround.c+(gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod3)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: sti.z ~ gndr.c + essround.c + (gndr.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1304085.5 1304162.5 -652035.7 1304071.5    441161 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -4.7821 -0.6917 -0.0481  0.6453  5.0769 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr 
##  cntry    (Intercept) 0.023236 0.15243       
##           gndr.c      0.002292 0.04787  -0.46
##  Residual             0.991511 0.99575       
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##              Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept) 4.937e-02  2.659e-02 3.292e+01   1.856   0.0724 .  
## gndr.c      1.914e-01  8.992e-03 3.207e+01  21.287   <2e-16 ***
## essround.c  6.470e-03  5.721e-04 4.381e+05  11.309   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##            (Intr) gndr.c
## gndr.c     -0.427       
## essround.c -0.003 -0.001
```

``` r
getFE(mod3,round=3)
```

```
##              Est.    SE         df      t     p     LL    UL
## (Intercept) 0.049 0.027     32.922  1.856 0.072 -0.005 0.103
## gndr.c      0.191 0.009     32.067 21.287 0.000  0.173 0.210
## essround.c  0.006 0.001 438063.075 11.309 0.000  0.005 0.008
```

``` r
getVC(mod3)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.15 0.02
## 2    cntry      gndr.c   <NA>  0.05 0.00
## 3    cntry (Intercept) gndr.c -0.46 0.00
## 4 Residual        <NA>   <NA>  1.00 0.99
```

``` r
r2mlm(mod3,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.009207357
## slope variation 0.000555607
## mean variation  0.022928203
## sigma2          0.967308832
## 
## $R2s
##           total
## f   0.009207357
## v   0.000555607
## m   0.022928203
## fv  0.009762964
## fvm 0.032691168
```

``` r
anova(mod2,mod3)
```

```
## Data: diff_dat
## Models:
## mod2: sti.z ~ gndr.c + (gndr.c | cntry)
## mod3: sti.z ~ gndr.c + essround.c + (gndr.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod2    6 1304211 1304277 -652100   1304199                         
## mod3    7 1304085 1304162 -652036   1304071 127.86  1  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## mod4: random effect of time (Ess round)


``` r
mod4<-lmer(sti.z~gndr.c+essround.c+(gndr.c+essround.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod4)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: sti.z ~ gndr.c + essround.c + (gndr.c + essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1303377.1 1303487.1 -651678.6 1303357.1    441158 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -4.9142 -0.6927 -0.0431  0.6448  5.0763 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr       
##  cntry    (Intercept) 0.032297 0.17971             
##           gndr.c      0.002381 0.04880  -0.51      
##           essround.c  0.001204 0.03469   0.31 -0.34
##  Residual             0.989569 0.99477             
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##              Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)  0.057531   0.031461 31.689607   1.829   0.0769 .  
## gndr.c       0.190998   0.009138 31.749239  20.902   <2e-16 ***
## essround.c   0.012406   0.006129 25.873284   2.024   0.0534 .  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##            (Intr) gndr.c
## gndr.c     -0.472       
## essround.c  0.298 -0.308
```

``` r
getFE(mod4,round=3)
```

```
##              Est.    SE     df      t     p     LL    UL
## (Intercept) 0.058 0.031 31.690  1.829 0.077 -0.007 0.122
## gndr.c      0.191 0.009 31.749 20.902 0.000  0.172 0.210
## essround.c  0.012 0.006 25.873  2.024 0.053  0.000 0.025
```

``` r
getVC(mod4)
```

```
##        grp        var1       var2 sdcor vcov
## 1    cntry (Intercept)       <NA>  0.18 0.03
## 2    cntry      gndr.c       <NA>  0.05 0.00
## 3    cntry  essround.c       <NA>  0.03 0.00
## 4    cntry (Intercept)     gndr.c -0.51 0.00
## 5    cntry (Intercept) essround.c  0.31 0.00
## 6    cntry      gndr.c essround.c -0.34 0.00
## 7 Residual        <NA>       <NA>  0.99 0.99
```

``` r
r2mlm(mod4,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.009846247
## slope variation 0.009286491
## mean variation  0.031246650
## sigma2          0.949620612
## 
## $R2s
##           total
## f   0.009846247
## v   0.009286491
## m   0.031246650
## fv  0.019132739
## fvm 0.050379388
```

``` r
anova(mod2,mod3,mod4)
```

```
## Data: diff_dat
## Models:
## mod2: sti.z ~ gndr.c + (gndr.c | cntry)
## mod3: sti.z ~ gndr.c + essround.c + (gndr.c | cntry)
## mod4: sti.z ~ gndr.c + essround.c + (gndr.c + essround.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod2    6 1304211 1304277 -652100   1304199                         
## mod3    7 1304085 1304162 -652036   1304071 127.86  1  < 2.2e-16 ***
## mod4   10 1303377 1303487 -651679   1303357 714.36  3  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## mod5: fixed interaction between time and gender


``` r
mod5<-lmer(sti.z~gndr.c+essround.c+
             gndr.c:essround.c+(gndr.c+essround.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod5)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: sti.z ~ gndr.c + essround.c + gndr.c:essround.c + (gndr.c + essround.c |      cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1303331.4 1303452.4 -651654.7 1303309.4    441157 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -4.9358 -0.6925 -0.0435  0.6441  5.0939 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr       
##  cntry    (Intercept) 0.032320 0.17978             
##           gndr.c      0.002642 0.05140  -0.55      
##           essround.c  0.001209 0.03477   0.31 -0.36
##  Residual             0.989460 0.99472             
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##                     Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)        5.743e-02  3.147e-02  3.165e+01   1.825   0.0775 .  
## gndr.c             1.914e-01  9.561e-03  3.139e+01  20.021  < 2e-16 ***
## essround.c         1.230e-02  6.142e-03  2.583e+01   2.002   0.0559 .  
## gndr.c:essround.c -7.823e-03  1.130e-03  6.191e+04  -6.922 4.51e-12 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c essrn.
## gndr.c      -0.514              
## essround.c   0.296 -0.335       
## gndr.c:ssr.  0.000 -0.008  0.002
```

``` r
getFE(mod5,round=3)
```

```
##                     Est.    SE        df      t     p     LL     UL
## (Intercept)        0.057 0.031    31.652  1.825 0.077 -0.007  0.122
## gndr.c             0.191 0.010    31.391 20.021 0.000  0.172  0.211
## essround.c         0.012 0.006    25.835  2.002 0.056  0.000  0.025
## gndr.c:essround.c -0.008 0.001 61906.334 -6.922 0.000 -0.010 -0.006
```

``` r
getVC(mod5)
```

```
##        grp        var1       var2 sdcor  vcov
## 1    cntry (Intercept)       <NA>  0.18  0.03
## 2    cntry      gndr.c       <NA>  0.05  0.00
## 3    cntry  essround.c       <NA>  0.03  0.00
## 4    cntry (Intercept)     gndr.c -0.55 -0.01
## 5    cntry (Intercept) essround.c  0.31  0.00
## 6    cntry      gndr.c essround.c -0.36  0.00
## 7 Residual        <NA>       <NA>  0.99  0.99
```

``` r
r2mlm(mod5,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.010038856
## slope variation 0.009384796
## mean variation  0.031307875
## sigma2          0.949268472
## 
## $R2s
##           total
## f   0.010038856
## v   0.009384796
## m   0.031307875
## fv  0.019423653
## fvm 0.050731528
```

``` r
anova(mod4,mod5)
```

```
## Data: diff_dat
## Models:
## mod4: sti.z ~ gndr.c + essround.c + (gndr.c + essround.c | cntry)
## mod5: sti.z ~ gndr.c + essround.c + gndr.c:essround.c + (gndr.c + essround.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod4   10 1303377 1303487 -651679   1303357                         
## mod5   11 1303331 1303452 -651655   1303309 47.712  1  4.936e-12 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## mod6: random interaction between time and gender


``` r
mod6<-lmer(sti.z~gndr.c+essround.c+
             gndr.c:essround.c+(gndr.c+essround.c+gndr.c:essround.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: sti.z ~ gndr.c + essround.c + gndr.c:essround.c + (gndr.c + essround.c +  
##     gndr.c:essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1303309.4 1303474.4 -651639.7 1303279.4    441153 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -4.9528 -0.6930 -0.0455  0.6445  5.1192 
## 
## Random effects:
##  Groups   Name              Variance  Std.Dev. Corr             
##  cntry    (Intercept)       3.216e-02 0.179338                  
##           gndr.c            2.935e-03 0.054171 -0.56            
##           essround.c        1.198e-03 0.034606  0.33 -0.41      
##           gndr.c:essround.c 9.244e-05 0.009614  0.25 -0.17 -0.20
##  Residual                   9.893e-01 0.994645                  
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##                    Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)        0.058117   0.031392 31.730455   1.851 0.073453 .  
## gndr.c             0.191413   0.010108 27.710486  18.936  < 2e-16 ***
## essround.c         0.012141   0.006112 25.880798   1.986 0.057693 .  
## gndr.c:essround.c -0.009017   0.002114 23.150956  -4.265 0.000287 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c essrn.
## gndr.c      -0.523              
## essround.c   0.314 -0.374       
## gndr.c:ssr.  0.198 -0.149 -0.162
```

``` r
getFE(mod6,round=3)
```

```
##                     Est.    SE     df      t     p     LL     UL
## (Intercept)        0.058 0.031 31.730  1.851 0.073 -0.006  0.122
## gndr.c             0.191 0.010 27.710 18.936 0.000  0.171  0.212
## essround.c         0.012 0.006 25.881  1.986 0.058  0.000  0.025
## gndr.c:essround.c -0.009 0.002 23.151 -4.265 0.000 -0.013 -0.005
```

``` r
getVC(mod6)
```

```
##         grp              var1              var2 sdcor  vcov
## 1     cntry       (Intercept)              <NA>  0.18  0.03
## 2     cntry            gndr.c              <NA>  0.05  0.00
## 3     cntry        essround.c              <NA>  0.03  0.00
## 4     cntry gndr.c:essround.c              <NA>  0.01  0.00
## 5     cntry       (Intercept)            gndr.c -0.56 -0.01
## 6     cntry       (Intercept)        essround.c  0.33  0.00
## 7     cntry       (Intercept) gndr.c:essround.c  0.25  0.00
## 8     cntry            gndr.c        essround.c -0.41  0.00
## 9     cntry            gndr.c gndr.c:essround.c -0.17  0.00
## 10    cntry        essround.c gndr.c:essround.c -0.20  0.00
## 11 Residual              <NA>              <NA>  0.99  0.99
```

``` r
r2mlm(mod6,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.010055451
## slope variation 0.009572081
## mean variation  0.031186459
## sigma2          0.949186008
## 
## $R2s
##           total
## f   0.010055451
## v   0.009572081
## m   0.031186459
## fv  0.019627532
## fvm 0.050813992
```

``` r
anova(mod4,mod5,mod6)
```

```
## Data: diff_dat
## Models:
## mod4: sti.z ~ gndr.c + essround.c + (gndr.c + essround.c | cntry)
## mod5: sti.z ~ gndr.c + essround.c + gndr.c:essround.c + (gndr.c + essround.c | cntry)
## mod6: sti.z ~ gndr.c + essround.c + gndr.c:essround.c + (gndr.c + essround.c + gndr.c:essround.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod4   10 1303377 1303487 -651679   1303357                         
## mod5   11 1303331 1303452 -651655   1303309 47.712  1  4.936e-12 ***
## mod6   15 1303309 1303474 -651640   1303279 29.997  4  4.900e-06 ***
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
##         4.5  0.0373 0.0508 Inf   -0.0622    0.1368   0.735  0.4622
##        -4.5 -0.1125 0.0377 Inf   -0.1863   -0.0387  -2.987  0.0028
## 
## gndr.c =  0.5:
##  essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##         4.5  0.1882 0.0455 Inf    0.0990    0.2774   4.134  <.0001
##        -4.5  0.1195 0.0330 Inf    0.0547    0.1843   3.615  0.0003
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
##  essround.c4.5 - (essround.c-4.5)   0.1498 0.0573 Inf    0.0375     0.262   2.614  0.0090
## 
## gndr.c =  0.5:
##  contrast                         estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   0.0687 0.0543 Inf   -0.0377     0.175   1.265  0.2058
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
##    -0.5        4.5  0.0373 0.0508 Inf   -0.0622    0.1368   0.735  0.4622
##     0.5        4.5  0.1882 0.0455 Inf    0.0990    0.2774   4.134  <.0001
##    -0.5       -4.5 -0.1125 0.0377 Inf   -0.1863   -0.0387  -2.987  0.0028
##     0.5       -4.5  0.1195 0.0330 Inf    0.0547    0.1843   3.615  0.0003
## 
## Degrees-of-freedom method: asymptotic 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod6,adjust="none",infer=c(T,T))
```

```
##  contrast                                                 estimate     SE  df asymp.LCL asymp.UCL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.1508 0.0128 Inf   -0.1759   -0.1257
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   0.1498 0.0573 Inf    0.0375    0.2622
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.0821 0.0595 Inf   -0.1988    0.0345
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.3007 0.0521 Inf    0.1986    0.4028
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)       0.0687 0.0543 Inf   -0.0377    0.1751
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.2320 0.0149 Inf   -0.2611   -0.2028
##  z.ratio p.value
##  -11.779  <.0001
##    2.614  0.0090
##   -1.380  0.1677
##    5.773  <.0001
##    1.265  0.2058
##  -15.593  <.0001
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
##  diff_ESS10    0.151 0.0128 Inf     0.126     0.176  11.779  <.0001
##  diff_ESS1     0.232 0.0149 Inf     0.203     0.261  15.593  <.0001
## 
## Degrees-of-freedom method: asymptotic 
## Confidence level used: 0.95
```

``` r
# Test for whether and how much the differences have changed
pairs(diff_mod6,infer=c(T,T))
```

```
##  contrast               estimate    SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10 - diff_ESS1  -0.0811 0.019 Inf    -0.118   -0.0439  -4.265  <.0001
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
      obs_mean_wt=weighted.mean(x=sti.z,w=pspwght),
      obs_sd_wt=sqrt(wtd.var(sti.z,pspwght)),
      obs_mean=mean(sti.z),
      obs_sd=sd(sti.z),
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
  ylab("Mean-level of value stimulation")+
  theme(legend.title=element_blank(),
        legend.text  = element_text(size = 16),
        axis.title.x = element_text(size = 18, face = "bold"),   
        axis.title.y = element_text(size = 18, face = "bold"),   
        axis.text.x  = element_text(size = 14),                  
        axis.text.y  = element_text(size = 14))

p_time_trends
```

![](Analysis_stimulation_files/figure-html/unnamed-chunk-33-1.png)<!-- -->

``` r
png(filename = 
      "../results/sti/time_trends.png",
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
pred_cntry_dat$sti.z_mean<-predict(mod6,newdata=pred_cntry_dat)

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

range(pred_cntry_dat$sti.z_mean)
```

```
## [1] -0.5139010  0.5031949
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
      obs_mean_wt=weighted.mean(x=sti.z,w=pspwght),
      obs_sd_wt=sqrt(wtd.var(sti.z,pspwght)),
      obs_mean=mean(sti.z),
      obs_sd=sd(sti.z),
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

pdf("../results/sti/country_specific_time_trends.pdf", width = 8, height = 6)

for (ctry in countries) {
  print(
    ggplot(pred_cntry_dat[pred_cntry_dat$cntry == ctry, ], 
           aes(x = year, y = sti.z_mean, color = gender)) +
      geom_smooth(method = "lm",formula = "y ~ x", se = FALSE) +
      geom_point(size=8) +
      geom_point(aes(x=year,y=obs_mean_wt),size=8,shape = 1,alpha=.50)+
      geom_errorbar(aes(ymin=obs_mean_wt_LL, ymax=obs_mean_wt_UL),alpha = .50)+
      scale_color_manual(values = my_colors) +
      geom_flag(aes(country=tolower(ctry)))+
      scale_y_continuous(limits = c(-1.1, 1.1)) +
      scale_x_continuous(limits = c(2001, 2021), breaks = seq(2002, 2020, 2)) +
      ggtitle(paste("Country:", ctry))+
    ylab("Mean-level of value stimulation")+
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
         aes(x = year, y = sti.z_mean, color = gender)) +
  geom_smooth(method = "lm", se = FALSE, formula="y~x") +
  geom_point(size=3.5) +
  geom_point(aes(x=year,y=obs_mean_wt),size=3.5,shape = 1,alpha=.70)+
  geom_errorbar(aes(ymin=obs_mean_wt_LL, ymax=obs_mean_wt_UL),alpha = .50)+
  scale_color_manual(values = my_colors) +
  #geom_flag(aes(country=tolower(ctry)))+
  scale_y_continuous(limits = c(-1.1, 1.1)) +
  scale_x_continuous(limits = c(2001, 2021), breaks = seq(2002, 2020, 2)) +
  #ggtitle(paste("Country:", ctry))+
  ylab("Mean-level of value stimulation")+
  xlab("Year")+
  theme(legend.title=element_blank(),legend.position = "top",
        axis.text.x = element_text(angle = 45,size = 6,hjust=1))+
  facet_wrap(~CLDR,nrow=6,ncol=6)+
  #facet_wrap(~cntry,nrow=6,ncol=6)+
  geom_flag(aes(country=tolower(cntry)),size=2)

facet_plot
```

![](Analysis_stimulation_files/figure-html/unnamed-chunk-35-1.png)<!-- -->

``` r
png(filename = 
      "../results/sti/country_time_trend_facets.png",
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
## [1] 12.95533
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
## 1    0.21               -0.06                    0.04                   -0.04                     -0.08
## 2    0.15                0.13                   -0.03                    0.11                      0.14
## 3    0.22                0.08                    0.00                    0.08                      0.08
## 4    0.20                0.17                   -0.02                    0.16                      0.18
## 5    0.20               -0.11                   -0.03                   -0.13                     -0.10
## 6    0.25                0.30                   -0.14                    0.23                      0.37
## 7    0.21                0.07                   -0.09                    0.03                      0.12
## 8    0.17                0.21                   -0.06                    0.18                      0.24
## 9    0.19               -0.03                   -0.05                   -0.06                     -0.01
## 10   0.11                0.02                   -0.15                   -0.06                      0.09
## 11   0.14                0.05                   -0.08                    0.01                      0.09
## 12   0.20               -0.11                   -0.06                   -0.14                     -0.08
## 13   0.16               -0.01                   -0.06                   -0.04                      0.02
## 14   0.20               -0.12                   -0.15                   -0.20                     -0.05
## 15   0.29                0.08                   -0.12                    0.02                      0.14
## 16   0.17                0.14                   -0.05                    0.11                      0.17
## 17   0.12                0.03                    0.03                    0.04                      0.01
## 18   0.17                0.07                    0.05                    0.09                      0.04
## 19   0.15                0.05                   -0.14                   -0.02                      0.12
## 20   0.18                0.09                   -0.10                    0.03                      0.14
## 21   0.33               -0.17                   -0.12                   -0.23                     -0.11
## 22   0.20               -0.48                   -0.03                   -0.50                     -0.47
## 23   0.23                0.49                   -0.10                    0.44                      0.54
## 24   0.17                0.03                   -0.03                    0.02                      0.05
## 25   0.20                0.13                   -0.11                    0.08                      0.19
## 26   0.25               -0.14                   -0.09                   -0.19                     -0.10
## 27   0.24               -0.15                   -0.21                   -0.26                     -0.04
## 28   0.25                0.15                   -0.25                    0.03                      0.28
## 29   0.16                0.28                   -0.03                    0.27                      0.30
## 30   0.17                0.13                   -0.03                    0.12                      0.14
## 31   0.22                0.21                   -0.15                    0.13                      0.28
## 32   0.07                1.43                   -0.15                    1.36                      1.51
## 33   0.14                0.67                   -0.18                    0.58                      0.75
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
## 1     LT               -0.17
## 2     PT               -0.15
## 3     PL               -0.14
## 4     GR               -0.12
## 5     CY               -0.11
## 6     FR               -0.11
## 7     AT               -0.06
## 8     EE               -0.03
## 9     GB               -0.01
## 10    ES                0.02
## 11    IE                0.03
## 12    NL                0.03
## 13    FI                0.05
## 14    IS                0.05
## 15    DE                0.07
## 16    IL                0.07
## 17    BG                0.08
## 18    BE                0.13
## 19    NO                0.13
## 20    SI                0.13
## 21    HU                0.14
## 22    RU                0.15
## 23    CH                0.17
## 24    DK                0.21
## 25    SK                0.21
## 26    SE                0.28
## 27    CZ                0.30
## 28    UA                0.67
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
## 1     RU                   -0.25
## 2     PT                   -0.21
## 3     UA                   -0.18
## 4     ES                   -0.15
## 5     GR                   -0.15
## 6     SK                   -0.15
## 7     CZ                   -0.14
## 8     IS                   -0.14
## 9     LT                   -0.12
## 10    NO                   -0.11
## 11    DE                   -0.09
## 12    PL                   -0.09
## 13    FI                   -0.08
## 14    DK                   -0.06
## 15    FR                   -0.06
## 16    GB                   -0.06
## 17    EE                   -0.05
## 18    HU                   -0.05
## 19    BE                   -0.03
## 20    CY                   -0.03
## 21    NL                   -0.03
## 22    SE                   -0.03
## 23    SI                   -0.03
## 24    CH                   -0.02
## 25    BG                    0.00
## 26    IE                    0.03
## 27    AT                    0.04
## 28    IL                    0.05
```


#### Country-level predictors: GEI


``` r
mod6_GEI<-lmer(sti.z~gndr.c+essround.c+
             gndr.c:essround.c+
               gei.z.cm:gndr.c+gei.z.cm:essround.c+gei.z.cm:gndr.c:essround.c+
               (gndr.c+essround.c+gndr.c:essround.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_GEI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: sti.z ~ gndr.c + essround.c + gndr.c:essround.c + gei.z.cm:gndr.c +  
##     gei.z.cm:essround.c + gei.z.cm:gndr.c:essround.c + (gndr.c +  
##     essround.c + gndr.c:essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1272570.1 1272767.7 -636267.1 1272534.1    431760 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -4.9508 -0.6944 -0.0445  0.6464  5.1121 
## 
## Random effects:
##  Groups   Name              Variance  Std.Dev. Corr             
##  cntry    (Intercept)       3.211e-02 0.179187                  
##           gndr.c            3.123e-03 0.055885 -0.74            
##           essround.c        9.607e-04 0.030994  0.34 -0.40      
##           gndr.c:essround.c 8.024e-05 0.008958  0.50 -0.35  0.11
##  Residual                   9.871e-01 0.993509                  
## Number of obs: 431778, groups:  cntry, 32
## 
## Fixed effects:
##                             Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)                 0.064342   0.031845 30.715896   2.020 0.052128 .  
## gndr.c                      0.194611   0.010522 24.030195  18.496 1.02e-15 ***
## essround.c                  0.010238   0.005576 24.433933   1.836 0.078548 .  
## gndr.c:essround.c          -0.008831   0.002059 27.516094  -4.289 0.000199 ***
## gndr.c:gei.z.cm            -0.022764   0.008234 37.063839  -2.765 0.008825 ** 
## essround.c:gei.z.cm        -0.008847   0.005442 26.385526  -1.626 0.115900    
## gndr.c:essround.c:gei.z.cm  0.004439   0.002255 31.550567   1.969 0.057806 .  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c essrn. gnd.:. gn.:.. es.:..
## gndr.c      -0.694                                   
## essround.c   0.324 -0.367                            
## gndr.c:ssr.  0.388 -0.260  0.080                     
## gndr.c:g.z. -0.004 -0.055 -0.005 -0.059              
## essrnd.c:.. -0.004 -0.001 -0.028  0.005 -0.173       
## gndr.c:.:..  0.002 -0.032  0.004 -0.192  0.190 -0.058
```

``` r
getFE(mod6_GEI,round=2,p.round = 5)
```

```
##                             Est.   SE    df     t       p    LL    UL
## (Intercept)                 0.06 0.03 30.72  2.02 0.05213  0.00  0.13
## gndr.c                      0.19 0.01 24.03 18.50 0.00000  0.17  0.22
## essround.c                  0.01 0.01 24.43  1.84 0.07855  0.00  0.02
## gndr.c:essround.c          -0.01 0.00 27.52 -4.29 0.00020 -0.01  0.00
## gndr.c:gei.z.cm            -0.02 0.01 37.06 -2.76 0.00882 -0.04 -0.01
## essround.c:gei.z.cm        -0.01 0.01 26.39 -1.63 0.11590 -0.02  0.00
## gndr.c:essround.c:gei.z.cm  0.00 0.00 31.55  1.97 0.05781  0.00  0.01
```

``` r
getVC(mod6_GEI)
```

```
##         grp              var1              var2 sdcor  vcov
## 1     cntry       (Intercept)              <NA>  0.18  0.03
## 2     cntry            gndr.c              <NA>  0.06  0.00
## 3     cntry        essround.c              <NA>  0.03  0.00
## 4     cntry gndr.c:essround.c              <NA>  0.01  0.00
## 5     cntry       (Intercept)            gndr.c -0.74 -0.01
## 6     cntry       (Intercept)        essround.c  0.34  0.00
## 7     cntry       (Intercept) gndr.c:essround.c  0.50  0.00
## 8     cntry            gndr.c        essround.c -0.40  0.00
## 9     cntry            gndr.c gndr.c:essround.c -0.35  0.00
## 10    cntry        essround.c gndr.c:essround.c  0.11  0.00
## 11 Residual              <NA>              <NA>  0.99  0.99
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_GEI,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] 19.78147
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_GEI,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] 13.1949
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
##         4.5  0.1502 0.0529 Inf    0.0465    0.2540   2.838  0.0045
##        -4.5 -0.0215 0.0419 Inf   -0.1036    0.0605  -0.514  0.6070
## 
## gei.z.cm =  0:
##  essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##         4.5  0.1104 0.0465 Inf    0.0193    0.2015   2.375  0.0175
##        -4.5  0.0183 0.0336 Inf   -0.0475    0.0840   0.544  0.5861
## 
## gei.z.cm =  1:
##  essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##         4.5  0.0706 0.0522 Inf   -0.0316    0.1728   1.354  0.1759
##        -4.5  0.0581 0.0412 Inf   -0.0227    0.1389   1.409  0.1587
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
##  essround.c4.5 - (essround.c-4.5)   0.1718 0.0711 Inf   0.03244     0.311   2.416  0.0157
## 
## gei.z.cm =  0:
##  contrast                         estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   0.0921 0.0502 Inf  -0.00622     0.190   1.836  0.0663
## 
## gei.z.cm =  1:
##  contrast                         estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   0.0125 0.0692 Inf  -0.12302     0.148   0.181  0.8564
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
##    -0.5        4.5  0.07139 0.0561 Inf   -0.0386    0.1814   1.272  0.2034
##     0.5        4.5  0.22905 0.0515 Inf    0.1280    0.3301   4.444  <.0001
##    -0.5       -4.5 -0.16008 0.0455 Inf   -0.2493   -0.0709  -3.519  0.0004
##     0.5       -4.5  0.11700 0.0406 Inf    0.0375    0.1965   2.885  0.0039
## 
## gei.z.cm =  0:
##  gndr.c essround.c   emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##    -0.5        4.5  0.03297 0.0489 Inf   -0.0629    0.1289   0.674  0.5004
##     0.5        4.5  0.18785 0.0447 Inf    0.1002    0.2755   4.200  <.0001
##    -0.5       -4.5 -0.09890 0.0377 Inf   -0.1727   -0.0251  -2.625  0.0087
##     0.5       -4.5  0.13545 0.0309 Inf    0.0748    0.1961   4.379  <.0001
## 
## gei.z.cm =  1:
##  gndr.c essround.c   emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##    -0.5        4.5 -0.00544 0.0551 Inf   -0.1135    0.1026  -0.099  0.9213
##     0.5        4.5  0.14664 0.0505 Inf    0.0477    0.2455   2.906  0.0037
##    -0.5       -4.5 -0.03772 0.0447 Inf   -0.1253    0.0499  -0.844  0.3987
##     0.5       -4.5  0.15389 0.0397 Inf    0.0760    0.2318   3.874  0.0001
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
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5     -0.15766 0.0202 Inf   -0.1973  -0.11799
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  0.23148 0.0726 Inf    0.0892   0.37376
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)  -0.04561 0.0759 Inf   -0.1944   0.10319
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)     0.38914 0.0687 Inf    0.2545   0.52376
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      0.11205 0.0727 Inf   -0.0304   0.25454
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5) -0.27709 0.0204 Inf   -0.3171  -0.23711
##  z.ratio p.value
##   -7.789  <.0001
##    3.189  0.0014
##   -0.601  0.5480
##    5.665  <.0001
##    1.541  0.1233
##  -13.586  <.0001
## 
## gei.z.cm =  0:
##  contrast                                                 estimate     SE  df asymp.LCL asymp.UCL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5     -0.15487 0.0121 Inf   -0.1785  -0.13120
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  0.13188 0.0503 Inf    0.0333   0.23046
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)  -0.10247 0.0549 Inf   -0.2101   0.00518
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)     0.28675 0.0473 Inf    0.1940   0.37953
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      0.05240 0.0518 Inf   -0.0490   0.15383
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5) -0.23435 0.0157 Inf   -0.2652  -0.20353
##  z.ratio p.value
##  -12.825  <.0001
##    2.622  0.0087
##   -1.866  0.0621
##    6.057  <.0001
##    1.013  0.3113
##  -14.903  <.0001
## 
## gei.z.cm =  1:
##  contrast                                                 estimate     SE  df asymp.LCL asymp.UCL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5     -0.15208 0.0169 Inf   -0.1853  -0.11888
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  0.03227 0.0701 Inf   -0.1051   0.16961
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)  -0.15934 0.0741 Inf   -0.3045  -0.01420
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)     0.18436 0.0665 Inf    0.0541   0.31463
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)     -0.00725 0.0704 Inf   -0.1453   0.13077
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5) -0.19161 0.0189 Inf   -0.2286  -0.15460
##  z.ratio p.value
##   -8.978  <.0001
##    0.461  0.6451
##   -2.152  0.0314
##    2.774  0.0055
##   -0.103  0.9180
##  -10.148  <.0001
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
##  diff_ESS10    0.158 0.0202 Inf     0.118     0.197   7.789  <.0001
##  diff_ESS1     0.277 0.0204 Inf     0.237     0.317  13.586  <.0001
## 
## gei.z.cm =  0:
##  contrast   estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10    0.155 0.0121 Inf     0.131     0.179  12.825  <.0001
##  diff_ESS1     0.234 0.0157 Inf     0.204     0.265  14.903  <.0001
## 
## gei.z.cm =  1:
##  contrast   estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10    0.152 0.0169 Inf     0.119     0.185   8.978  <.0001
##  diff_ESS1     0.192 0.0189 Inf     0.155     0.229  10.148  <.0001
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
##  diff_ESS10 - diff_ESS1  -0.1194 0.0300 Inf    -0.178  -0.06065  -3.982  0.0001
## 
## gei.z.cm =  0:
##  contrast               estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10 - diff_ESS1  -0.0795 0.0185 Inf    -0.116  -0.04316  -4.289  <.0001
## 
## gei.z.cm =  1:
##  contrast               estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10 - diff_ESS1  -0.0395 0.0247 Inf    -0.088   0.00892  -1.599  0.1098
## 
## Degrees-of-freedom method: asymptotic 
## Confidence level used: 0.95
```


#### Country-level predictors: GGGI


``` r
mod6_GGGI<-lmer(sti.z~gndr.c+essround.c+
             gndr.c:essround.c+
               gggi.z.cm:gndr.c+gggi.z.cm:essround.c+gggi.z.cm:gndr.c:essround.c+
               (gndr.c+essround.c+gndr.c:essround.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_GGGI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: sti.z ~ gndr.c + essround.c + gndr.c:essround.c + gggi.z.cm:gndr.c +  
##     gggi.z.cm:essround.c + gggi.z.cm:gndr.c:essround.c + (gndr.c +  
##     essround.c + gndr.c:essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##  930847.9  931039.7 -465405.9  930811.9    314628 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -5.0717 -0.6922 -0.0367  0.6481  5.0954 
## 
## Random effects:
##  Groups   Name              Variance Std.Dev. Corr             
##  cntry    (Intercept)       0.028396 0.16851                   
##           gndr.c            0.002218 0.04710  -0.47            
##           essround.c        0.000348 0.01865  -0.27  0.03      
##           gndr.c:essround.c 0.000141 0.01187   0.32 -0.21  0.13
##  Residual                   0.991854 0.99592                   
## Number of obs: 314646, groups:  cntry, 33
## 
## Fixed effects:
##                               Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)                  0.0726571  0.0295645 32.5786630   2.458  0.01948 *  
## gndr.c                       0.1923752  0.0092545 29.0894389  20.787  < 2e-16 ***
## essround.c                   0.0081813  0.0035572 26.8522394   2.300  0.02946 *  
## gndr.c:essround.c           -0.0090983  0.0028223 25.8615195  -3.224  0.00341 ** 
## gndr.c:gggi.z.cm            -0.0101975  0.0090508 34.3774439  -1.127  0.26768    
## essround.c:gggi.z.cm        -0.0023124  0.0037712 31.0335224  -0.613  0.54422    
## gndr.c:essround.c:gggi.z.cm  0.0005036  0.0030518 30.2825791   0.165  0.87003    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c essrn. gnd.:. gn.:.. es.:..
## gndr.c      -0.415                                   
## essround.c  -0.258  0.026                            
## gndr.c:ssr.  0.233 -0.184  0.101                     
## gndr.c:gg.. -0.002 -0.030 -0.008 -0.032              
## essrnd.c:.. -0.008 -0.003 -0.082 -0.009 -0.071       
## gndr.c:.:.. -0.001 -0.030 -0.006 -0.079 -0.098  0.149
```

``` r
getFE(mod6_GGGI,round=2,p.round = 5)
```

```
##                              Est.   SE    df     t       p    LL   UL
## (Intercept)                  0.07 0.03 32.58  2.46 0.01948  0.01 0.13
## gndr.c                       0.19 0.01 29.09 20.79 0.00000  0.17 0.21
## essround.c                   0.01 0.00 26.85  2.30 0.02946  0.00 0.02
## gndr.c:essround.c           -0.01 0.00 25.86 -3.22 0.00341 -0.01 0.00
## gndr.c:gggi.z.cm            -0.01 0.01 34.38 -1.13 0.26768 -0.03 0.01
## essround.c:gggi.z.cm         0.00 0.00 31.03 -0.61 0.54422 -0.01 0.01
## gndr.c:essround.c:gggi.z.cm  0.00 0.00 30.28  0.17 0.87003 -0.01 0.01
```

``` r
getVC(mod6_GGGI)
```

```
##         grp              var1              var2 sdcor vcov
## 1     cntry       (Intercept)              <NA>  0.17 0.03
## 2     cntry            gndr.c              <NA>  0.05 0.00
## 3     cntry        essround.c              <NA>  0.02 0.00
## 4     cntry gndr.c:essround.c              <NA>  0.01 0.00
## 5     cntry       (Intercept)            gndr.c -0.47 0.00
## 6     cntry       (Intercept)        essround.c -0.27 0.00
## 7     cntry       (Intercept) gndr.c:essround.c  0.32 0.00
## 8     cntry            gndr.c        essround.c  0.03 0.00
## 9     cntry            gndr.c gndr.c:essround.c -0.21 0.00
## 10    cntry        essround.c gndr.c:essround.c  0.13 0.00
## 11 Residual              <NA>              <NA>  1.00 0.99
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_GGGI,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] 70.94086
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_GGGI,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] -52.5114
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
##         4.5 0.1199 0.0350 41.6   0.0491    0.191   3.421  0.0014
##        -4.5 0.0254 0.0412 40.6  -0.0578    0.109   0.617  0.5405
## 
## gggi.z.cm =  0:
##  essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5 0.1095 0.0298 28.2   0.0485    0.170   3.677  0.0010
##        -4.5 0.0358 0.0371 31.0  -0.0398    0.111   0.967  0.3411
## 
## gggi.z.cm =  1:
##  essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5 0.0991 0.0335 38.1   0.0313    0.167   2.958  0.0053
##        -4.5 0.0462 0.0403 39.2  -0.0353    0.128   1.147  0.2584
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
##  essround.c4.5 - (essround.c-4.5)   0.0944 0.0485 30.8 -0.00456    0.193   1.946  0.0608
## 
## gggi.z.cm =  0:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   0.0736 0.0320 26.9  0.00793    0.139   2.300  0.0295
## 
## gggi.z.cm =  1:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   0.0528 0.0447 29.2 -0.03858    0.144   1.182  0.2469
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
##    -0.5        4.5  0.0402 0.0364 41.3  -0.0333  0.11366   1.105  0.2757
##     0.5        4.5  0.1996 0.0371 41.5   0.1246  0.27455   5.372  <.0001
##    -0.5       -4.5 -0.0975 0.0446 40.4  -0.1875 -0.00743  -2.187  0.0346
##     0.5       -4.5  0.1483 0.0414 41.4   0.0648  0.23187   3.585  0.0009
## 
## gggi.z.cm =  0:
##  gndr.c essround.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5        4.5  0.0338 0.0306 28.2  -0.0290  0.09651   1.101  0.2800
##     0.5        4.5  0.1852 0.0306 27.8   0.1225  0.24786   6.055  <.0001
##    -0.5       -4.5 -0.0808 0.0404 31.1  -0.1632  0.00153  -2.001  0.0542
##     0.5       -4.5  0.1525 0.0355 30.1   0.0799  0.22505   4.292  0.0002
## 
## gggi.z.cm =  1:
##  gndr.c essround.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5        4.5  0.0273 0.0347 38.1  -0.0430  0.09761   0.787  0.4364
##     0.5        4.5  0.1708 0.0352 38.1   0.0996  0.24204   4.855  <.0001
##    -0.5       -4.5 -0.0642 0.0437 39.0  -0.1525  0.02418  -1.469  0.1498
##     0.5       -4.5  0.1567 0.0403 39.7   0.0751  0.23822   3.884  0.0004
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
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.1594 0.0223 26.6 -0.20516 -0.11357  -7.145
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   0.1377 0.0500 29.7  0.03553  0.23977   2.754
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.1081 0.0505 30.9 -0.21112 -0.00514  -2.142
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.2970 0.0501 32.0  0.19503  0.39900   5.932
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)       0.0512 0.0545 31.2 -0.05984  0.16231   0.941
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.2458 0.0245 34.4 -0.29565 -0.19591 -10.012
##  p.value
##   <.0001
##   0.0099
##   0.0402
##   <.0001
##   0.3542
##   <.0001
## 
## gggi.z.cm =  0:
##  contrast                                                 estimate     SE   df lower.CL upper.CL t.ratio
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.1514 0.0143 26.0 -0.18077 -0.12210 -10.612
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   0.1146 0.0332 25.1  0.04614  0.18301   3.448
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.1187 0.0331 27.6 -0.18658 -0.05090  -3.588
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.2660 0.0336 26.7  0.19712  0.33489   7.928
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)       0.0327 0.0356 27.7 -0.04030  0.10568   0.918
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.2333 0.0170 27.5 -0.26824 -0.19839 -13.695
##  p.value
##   <.0001
##   0.0020
##   0.0013
##   <.0001
##   0.3666
##   <.0001
## 
## gggi.z.cm =  1:
##  contrast                                                 estimate     SE   df lower.CL upper.CL t.ratio
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.1435 0.0201 26.1 -0.18471 -0.10229  -7.156
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   0.0915 0.0460 28.3 -0.00262  0.18561   1.990
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.1294 0.0469 28.8 -0.22523 -0.03349  -2.761
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.2350 0.0461 30.3  0.14088  0.32912   5.097
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)       0.0141 0.0503 28.6 -0.08875  0.11704   0.281
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.2209 0.0238 29.1 -0.26957 -0.17214  -9.271
##  p.value
##   <.0001
##   0.0563
##   0.0099
##   <.0001
##   0.7805
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
##  diff_ESS10    0.159 0.0223 26.6    0.114    0.205   7.145  <.0001
##  diff_ESS1     0.246 0.0245 34.4    0.196    0.296  10.012  <.0001
## 
## gggi.z.cm =  0:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10    0.151 0.0143 26.0    0.122    0.181  10.612  <.0001
##  diff_ESS1     0.233 0.0170 27.5    0.198    0.268  13.695  <.0001
## 
## gggi.z.cm =  1:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10    0.144 0.0201 26.1    0.102    0.185   7.156  <.0001
##  diff_ESS1     0.221 0.0238 29.1    0.172    0.270   9.271  <.0001
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
##  diff_ESS10 - diff_ESS1  -0.0864 0.0389 30.4   -0.166 -0.00709  -2.224  0.0338
## 
## gggi.z.cm =  0:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10 - diff_ESS1  -0.0819 0.0254 25.9   -0.134 -0.02966  -3.224  0.0034
## 
## gggi.z.cm =  1:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10 - diff_ESS1  -0.0774 0.0359 27.1   -0.151 -0.00370  -2.155  0.0402
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

#### Country-level predictors: GDI


``` r
mod6_GDI<-lmer(sti.z~gndr.c+essround.c+
             gndr.c:essround.c+
               gdi.z.cm:gndr.c+gdi.z.cm:essround.c+gdi.z.cm:gndr.c:essround.c+
               (gndr.c+essround.c+gndr.c:essround.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_GDI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: sti.z ~ gndr.c + essround.c + gndr.c:essround.c + gdi.z.cm:gndr.c +  
##     gdi.z.cm:essround.c + gdi.z.cm:gndr.c:essround.c + (gndr.c +  
##     essround.c + gndr.c:essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1303305.5 1303503.4 -651634.7 1303269.5    441150 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -4.9651 -0.6929 -0.0453  0.6444  5.1206 
## 
## Random effects:
##  Groups   Name              Variance  Std.Dev. Corr             
##  cntry    (Intercept)       3.264e-02 0.180654                  
##           gndr.c            2.502e-03 0.050024 -0.47            
##           essround.c        8.873e-04 0.029788  0.11 -0.22      
##           gndr.c:essround.c 8.797e-05 0.009379  0.21 -0.13 -0.32
##  Residual                   9.893e-01 0.994641                  
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##                             Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)                 0.056776   0.031626 31.529212   1.795 0.082205 .  
## gndr.c                      0.191460   0.009441 27.193809  20.280  < 2e-16 ***
## essround.c                  0.012361   0.005292 26.837923   2.336 0.027212 *  
## gndr.c:essround.c          -0.009114   0.002079 21.516101  -4.383 0.000248 ***
## gndr.c:gdi.z.cm             0.013299   0.009180 36.841506   1.449 0.155866    
## essround.c:gdi.z.cm        -0.017432   0.005494 29.785071  -3.173 0.003491 ** 
## gndr.c:essround.c:gdi.z.cm -0.001816   0.002498 28.739693  -0.727 0.473138    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c essrn. gnd.:. gn.:.. es.:..
## gndr.c      -0.435                                   
## essround.c   0.101 -0.198                            
## gndr.c:ssr.  0.162 -0.113 -0.250                     
## gndr.c:gd..  0.000 -0.025  0.000 -0.016              
## essrnd.c:..  0.002 -0.001 -0.027  0.011 -0.154       
## gndr.c:.:..  0.001 -0.023  0.010  0.007  0.078 -0.231
```

``` r
getFE(mod6_GDI,round=2,p.round = 5)
```

```
##                             Est.   SE    df     t       p    LL    UL
## (Intercept)                 0.06 0.03 31.53  1.80 0.08221 -0.01  0.12
## gndr.c                      0.19 0.01 27.19 20.28 0.00000  0.17  0.21
## essround.c                  0.01 0.01 26.84  2.34 0.02721  0.00  0.02
## gndr.c:essround.c          -0.01 0.00 21.52 -4.38 0.00025 -0.01  0.00
## gndr.c:gdi.z.cm             0.01 0.01 36.84  1.45 0.15587 -0.01  0.03
## essround.c:gdi.z.cm        -0.02 0.01 29.79 -3.17 0.00349 -0.03 -0.01
## gndr.c:essround.c:gdi.z.cm  0.00 0.00 28.74 -0.73 0.47314 -0.01  0.00
```

``` r
getVC(mod6_GDI)
```

```
##         grp              var1              var2 sdcor vcov
## 1     cntry       (Intercept)              <NA>  0.18 0.03
## 2     cntry            gndr.c              <NA>  0.05 0.00
## 3     cntry        essround.c              <NA>  0.03 0.00
## 4     cntry gndr.c:essround.c              <NA>  0.01 0.00
## 5     cntry       (Intercept)            gndr.c -0.47 0.00
## 6     cntry       (Intercept)        essround.c  0.11 0.00
## 7     cntry       (Intercept) gndr.c:essround.c  0.21 0.00
## 8     cntry            gndr.c        essround.c -0.22 0.00
## 9     cntry            gndr.c gndr.c:essround.c -0.13 0.00
## 10    cntry        essround.c gndr.c:essround.c -0.32 0.00
## 11 Residual              <NA>              <NA>  0.99 0.99
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_GDI,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] 25.90343
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_GDI,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] 4.833913
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
##         4.5  0.19084 0.0486 Inf   0.09564    0.2860   3.929  0.0001
##        -4.5 -0.07729 0.0454 Inf  -0.16629    0.0117  -1.702  0.0887
## 
## gdi.z.cm =  0:
##  essround.c   emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##         4.5  0.11240 0.0415 Inf   0.03115    0.1937   2.711  0.0067
##        -4.5  0.00115 0.0376 Inf  -0.07259    0.0749   0.031  0.9756
## 
## gdi.z.cm =  1:
##  essround.c   emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##         4.5  0.03396 0.0480 Inf  -0.06005    0.1280   0.708  0.4789
##        -4.5  0.07959 0.0446 Inf  -0.00788    0.1671   1.783  0.0745
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
##  essround.c4.5 - (essround.c-4.5)   0.2681 0.0696 Inf    0.1318    0.4045   3.853  0.0001
## 
## gdi.z.cm =  0:
##  contrast                         estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   0.1112 0.0476 Inf    0.0179    0.2046   2.336  0.0195
## 
## gdi.z.cm =  1:
##  contrast                         estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  essround.c4.5 - (essround.c-4.5)  -0.0456 0.0677 Inf   -0.1783    0.0871  -0.674  0.5003
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
##    -0.5        4.5  0.1182 0.0524 Inf    0.0155    0.2209   2.256  0.0241
##     0.5        4.5  0.2635 0.0466 Inf    0.1722    0.3548   5.655  <.0001
##    -0.5       -4.5 -0.1828 0.0488 Inf   -0.2784   -0.0871  -3.745  0.0002
##     0.5       -4.5  0.0282 0.0440 Inf   -0.0580    0.1144   0.641  0.5214
## 
## gdi.z.cm =  0:
##  gndr.c essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##    -0.5        4.5  0.0372 0.0440 Inf   -0.0491    0.1235   0.844  0.3986
##     0.5        4.5  0.1876 0.0397 Inf    0.1098    0.2654   4.726  <.0001
##    -0.5       -4.5 -0.1151 0.0407 Inf   -0.1948   -0.0354  -2.830  0.0047
##     0.5       -4.5  0.1174 0.0357 Inf    0.0474    0.1874   3.287  0.0010
## 
## gdi.z.cm =  1:
##  gndr.c essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##    -0.5        4.5 -0.0438 0.0517 Inf   -0.1451    0.0574  -0.848  0.3962
##     0.5        4.5  0.1117 0.0460 Inf    0.0215    0.2019   2.428  0.0152
##    -0.5       -4.5 -0.0474 0.0480 Inf   -0.1415    0.0467  -0.987  0.3237
##     0.5       -4.5  0.2066 0.0433 Inf    0.1216    0.2915   4.767  <.0001
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
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5     -0.14532 0.0199 Inf   -0.1842   -0.1064
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  0.30098 0.0745 Inf    0.1549    0.4470
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   0.08997 0.0730 Inf   -0.0532    0.2331
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)     0.44630 0.0686 Inf    0.3119    0.5807
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      0.23530 0.0675 Inf    0.1030    0.3676
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5) -0.21100 0.0197 Inf   -0.2495   -0.1725
##  z.ratio p.value
##   -7.320  <.0001
##    4.039  0.0001
##    1.232  0.2180
##    6.507  <.0001
##    3.486  0.0005
##  -10.736  <.0001
## 
## gdi.z.cm =  0:
##  contrast                                                 estimate     SE  df asymp.LCL asymp.UCL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5     -0.15045 0.0125 Inf   -0.1750   -0.1259
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  0.15226 0.0508 Inf    0.0527    0.2518
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)  -0.08021 0.0503 Inf   -0.1789    0.0185
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)     0.30271 0.0467 Inf    0.2112    0.3942
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      0.07024 0.0462 Inf   -0.0203    0.1608
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5) -0.23247 0.0140 Inf   -0.2600   -0.2050
##  z.ratio p.value
##  -12.017  <.0001
##    2.999  0.0027
##   -1.593  0.1111
##    6.484  <.0001
##    1.521  0.1283
##  -16.579  <.0001
## 
## gdi.z.cm =  1:
##  contrast                                                 estimate     SE  df asymp.LCL asymp.UCL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5     -0.15558 0.0193 Inf   -0.1934   -0.1177
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  0.00355 0.0725 Inf   -0.1386    0.1457
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)  -0.25040 0.0712 Inf   -0.3900   -0.1108
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)     0.15912 0.0666 Inf    0.0286    0.2896
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)     -0.09482 0.0659 Inf   -0.2239    0.0343
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5) -0.25394 0.0199 Inf   -0.2929   -0.2149
##  z.ratio p.value
##   -8.057  <.0001
##    0.049  0.9610
##   -3.516  0.0004
##    2.390  0.0169
##   -1.439  0.1501
##  -12.764  <.0001
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
##  diff_ESS10    0.145 0.0199 Inf     0.106     0.184   7.320  <.0001
##  diff_ESS1     0.211 0.0197 Inf     0.172     0.250  10.736  <.0001
## 
## gdi.z.cm =  0:
##  contrast   estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10    0.150 0.0125 Inf     0.126     0.175  12.017  <.0001
##  diff_ESS1     0.232 0.0140 Inf     0.205     0.260  16.579  <.0001
## 
## gdi.z.cm =  1:
##  contrast   estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10    0.156 0.0193 Inf     0.118     0.193   8.057  <.0001
##  diff_ESS1     0.254 0.0199 Inf     0.215     0.293  12.764  <.0001
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
##  diff_ESS10 - diff_ESS1  -0.0657 0.0291 Inf    -0.123  -0.00855  -2.253  0.0242
## 
## gdi.z.cm =  0:
##  contrast               estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10 - diff_ESS1  -0.0820 0.0187 Inf    -0.119  -0.04534  -4.383  <.0001
## 
## gdi.z.cm =  1:
##  contrast               estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10 - diff_ESS1  -0.0984 0.0294 Inf    -0.156  -0.04084  -3.351  0.0008
## 
## Degrees-of-freedom method: asymptotic 
## Confidence level used: 0.95
```

#### Country-level predictors: log_GDP


``` r
mod6_log_GDP<-lmer(sti.z~gndr.c+essround.c+
             gndr.c:essround.c+
               log_gdp.z.cm:gndr.c+log_gdp.z.cm:essround.c+log_gdp.z.cm:gndr.c:essround.c+
               (gndr.c+essround.c+gndr.c:essround.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_log_GDP)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: sti.z ~ gndr.c + essround.c + gndr.c:essround.c + log_gdp.z.cm:gndr.c +  
##     log_gdp.z.cm:essround.c + log_gdp.z.cm:gndr.c:essround.c +  
##     (gndr.c + essround.c + gndr.c:essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1303303.1 1303501.0 -651633.5 1303267.1    441150 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -4.9542 -0.6927 -0.0449  0.6439  5.1437 
## 
## Random effects:
##  Groups   Name              Variance  Std.Dev. Corr             
##  cntry    (Intercept)       3.242e-02 0.180060                  
##           gndr.c            2.878e-03 0.053650 -0.57            
##           essround.c        1.093e-03 0.033059  0.32 -0.50      
##           gndr.c:essround.c 6.759e-05 0.008221  0.34 -0.11 -0.07
##  Residual                   9.893e-01 0.994644                  
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##                                 Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)                     0.057328   0.031517 31.594847   1.819   0.0784 .  
## gndr.c                          0.192376   0.009998 26.141702  19.242  < 2e-16 ***
## essround.c                      0.012185   0.005847 26.671684   2.084   0.0469 *  
## gndr.c:essround.c              -0.009293   0.001935 27.282331  -4.804 5.04e-05 ***
## gndr.c:log_gdp.z.cm            -0.014892   0.008700 31.469534  -1.712   0.0968 .  
## essround.c:log_gdp.z.cm        -0.011232   0.005715 28.845041  -1.965   0.0591 .  
## gndr.c:essround.c:log_gdp.z.cm  0.004243   0.002001 25.315610   2.120   0.0440 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c essrn. gnd.:. g.:_.. e.:_..
## gndr.c      -0.535                                   
## essround.c   0.310 -0.459                            
## gndr.c:ssr.  0.250 -0.088 -0.050                     
## gndr.c:l_.. -0.005 -0.034 -0.004 -0.003              
## essrnd.:_..  0.014 -0.010 -0.015  0.001 -0.354       
## gndr.:.:_..  0.002  0.001 -0.002 -0.172  0.077 -0.128
```

``` r
getFE(mod6_log_GDP,round=2,p.round = 5)
```

```
##                                 Est.   SE    df     t       p    LL    UL
## (Intercept)                     0.06 0.03 31.59  1.82 0.07841 -0.01  0.12
## gndr.c                          0.19 0.01 26.14 19.24 0.00000  0.17  0.21
## essround.c                      0.01 0.01 26.67  2.08 0.04686  0.00  0.02
## gndr.c:essround.c              -0.01 0.00 27.28 -4.80 0.00005 -0.01 -0.01
## gndr.c:log_gdp.z.cm            -0.01 0.01 31.47 -1.71 0.09678 -0.03  0.00
## essround.c:log_gdp.z.cm        -0.01 0.01 28.85 -1.97 0.05906 -0.02  0.00
## gndr.c:essround.c:log_gdp.z.cm  0.00 0.00 25.32  2.12 0.04398  0.00  0.01
```

``` r
getVC(mod6_log_GDP)
```

```
##         grp              var1              var2 sdcor  vcov
## 1     cntry       (Intercept)              <NA>  0.18  0.03
## 2     cntry            gndr.c              <NA>  0.05  0.00
## 3     cntry        essround.c              <NA>  0.03  0.00
## 4     cntry gndr.c:essround.c              <NA>  0.01  0.00
## 5     cntry       (Intercept)            gndr.c -0.57 -0.01
## 6     cntry       (Intercept)        essround.c  0.32  0.00
## 7     cntry       (Intercept) gndr.c:essround.c  0.34  0.00
## 8     cntry            gndr.c        essround.c -0.50  0.00
## 9     cntry            gndr.c gndr.c:essround.c -0.11  0.00
## 10    cntry        essround.c gndr.c:essround.c -0.07  0.00
## 11 Residual              <NA>              <NA>  0.99  0.99
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_log_GDP,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] 8.736416
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_log_GDP,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] 26.88462
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
##  essround.c   emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##         4.5  0.16271 0.0535 Inf    0.0579    0.2675   3.043  0.0023
##        -4.5 -0.04805 0.0433 Inf   -0.1329    0.0368  -1.109  0.2673
## 
## log_gdp.z.cm =  0:
##  essround.c   emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##         4.5  0.11216 0.0469 Inf    0.0202    0.2041   2.391  0.0168
##        -4.5  0.00249 0.0342 Inf   -0.0646    0.0696   0.073  0.9419
## 
## log_gdp.z.cm =  1:
##  essround.c   emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##         4.5  0.06161 0.0535 Inf   -0.0433    0.1665   1.151  0.2496
##        -4.5  0.05304 0.0423 Inf   -0.0298    0.1359   1.254  0.2097
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
##  essround.c4.5 - (essround.c-4.5)  0.21076 0.0741 Inf   0.06544     0.356   2.843  0.0045
## 
## log_gdp.z.cm =  0:
##  contrast                         estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  essround.c4.5 - (essround.c-4.5)  0.10967 0.0526 Inf   0.00653     0.213   2.084  0.0371
## 
## log_gdp.z.cm =  1:
##  contrast                         estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  essround.c4.5 - (essround.c-4.5)  0.00857 0.0730 Inf  -0.13454     0.152   0.117  0.9065
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
##  gndr.c essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##    -0.5        4.5  0.0895 0.0573 Inf -0.022822    0.2019   1.562  0.1183
##     0.5        4.5  0.2359 0.0511 Inf  0.135692    0.3361   4.614  <.0001
##    -0.5       -4.5 -0.1821 0.0452 Inf -0.270824   -0.0935  -4.025  0.0001
##     0.5       -4.5  0.0860 0.0435 Inf  0.000849    0.1712   1.980  0.0478
## 
## log_gdp.z.cm =  0:
##  gndr.c essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##    -0.5        4.5  0.0369 0.0497 Inf -0.060579    0.1343   0.742  0.4583
##     0.5        4.5  0.1874 0.0448 Inf  0.099607    0.2753   4.183  <.0001
##    -0.5       -4.5 -0.1146 0.0367 Inf -0.186509   -0.0427  -3.124  0.0018
##     0.5       -4.5  0.1196 0.0330 Inf  0.054858    0.1843   3.621  0.0003
## 
## log_gdp.z.cm =  1:
##  gndr.c essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##    -0.5        4.5 -0.0158 0.0573 Inf -0.128088    0.0966  -0.275  0.7832
##     0.5        4.5  0.1390 0.0509 Inf  0.039196    0.2388   2.730  0.0063
##    -0.5       -4.5 -0.0471 0.0441 Inf -0.133410    0.0393  -1.068  0.2854
##     0.5       -4.5  0.1531 0.0423 Inf  0.070299    0.2360   3.623  0.0003
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
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5     -0.14636 0.0190 Inf  -0.18368   -0.1090
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  0.27167 0.0764 Inf   0.12184    0.4215
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   0.00349 0.0804 Inf  -0.15410    0.1611
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)     0.41803 0.0699 Inf   0.28094    0.5551
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      0.14985 0.0743 Inf   0.00425    0.2954
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5) -0.26818 0.0192 Inf  -0.30580   -0.2306
##  z.ratio p.value
##   -7.687  <.0001
##    3.554  0.0004
##    0.043  0.9654
##    5.976  <.0001
##    2.017  0.0437
##  -13.972  <.0001
## 
## log_gdp.z.cm =  0:
##  contrast                                                 estimate     SE  df asymp.LCL asymp.UCL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5     -0.15056 0.0127 Inf  -0.17539   -0.1257
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  0.15148 0.0538 Inf   0.04612    0.2568
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)  -0.08271 0.0579 Inf  -0.19619    0.0308
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)     0.30204 0.0488 Inf   0.20632    0.3978
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      0.06785 0.0529 Inf  -0.03585    0.1715
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5) -0.23419 0.0138 Inf  -0.26128   -0.2071
##  z.ratio p.value
##  -11.885  <.0001
##    2.818  0.0048
##   -1.429  0.1531
##    6.184  <.0001
##    1.282  0.1997
##  -16.946  <.0001
## 
## log_gdp.z.cm =  1:
##  contrast                                                 estimate     SE  df asymp.LCL asymp.UCL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5     -0.15476 0.0172 Inf  -0.18848   -0.1210
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  0.03130 0.0750 Inf  -0.11573    0.1783
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)  -0.16891 0.0795 Inf  -0.32467   -0.0132
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)     0.18606 0.0685 Inf   0.05187    0.3202
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)     -0.01415 0.0728 Inf  -0.15678    0.1285
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5) -0.20021 0.0174 Inf  -0.23433   -0.1661
##  z.ratio p.value
##   -8.995  <.0001
##    0.417  0.6765
##   -2.125  0.0335
##    2.718  0.0066
##   -0.194  0.8458
##  -11.499  <.0001
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
##  diff_ESS10    0.146 0.0190 Inf     0.109     0.184   7.687  <.0001
##  diff_ESS1     0.268 0.0192 Inf     0.231     0.306  13.972  <.0001
## 
## log_gdp.z.cm =  0:
##  contrast   estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10    0.151 0.0127 Inf     0.126     0.175  11.885  <.0001
##  diff_ESS1     0.234 0.0138 Inf     0.207     0.261  16.946  <.0001
## 
## log_gdp.z.cm =  1:
##  contrast   estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10    0.155 0.0172 Inf     0.121     0.188   8.995  <.0001
##  diff_ESS1     0.200 0.0174 Inf     0.166     0.234  11.499  <.0001
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
##  diff_ESS10 - diff_ESS1  -0.1218 0.0271 Inf   -0.1750 -0.068659  -4.491  <.0001
## 
## log_gdp.z.cm =  0:
##  contrast               estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10 - diff_ESS1  -0.0836 0.0174 Inf   -0.1178 -0.049509  -4.804  <.0001
## 
## log_gdp.z.cm =  1:
##  contrast               estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10 - diff_ESS1  -0.0454 0.0228 Inf   -0.0901 -0.000778  -1.994  0.0461
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

