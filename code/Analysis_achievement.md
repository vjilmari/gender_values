---
title: "Analysis for Examining the Gender Equality Paradox in Values Using achievement Value"
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
cntry.ach<-diff_dat %>% group_by(cntry,essround) %>%
  summarise(ach.ctm=mean(ach,na.rm=T),
            ach.ctsd=sd(ach,na.rm=T)) %>%
  group_by(cntry) %>%
  summarise(ach.cm=mean(ach.ctm),
            ach.csd=mean(ach.ctsd)) 
```

```
## `summarise()` has grouped output by 'cntry'. You can override using the `.groups` argument.
```

``` r
grand_mean_ach<-mean(cntry.ach$ach.cm)
grand_sd_ach<-mean(cntry.ach$ach.csd)

# standardized
diff_dat$ach.z<-(diff_dat$ach-grand_mean_ach)/grand_sd_ach
hist(diff_dat$ach.z)
```

![](Analysis_achievement_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

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

# value-based achievement

cntry_ach_frame<-
  diff_dat %>% group_by(cntry,essround) %>%
  summarise('ach M' = weighted.mean(x=ach.z,w=pspwght),
            'ach SD' = sqrt(wtd.var(ach.z,pspwght))) %>%
  group_by(cntry) %>%
  summarise('ach M' = mean(x=`ach M`),
            'ach SD'= mean(x=`ach SD`))
```

```
## `summarise()` has grouped output by 'cntry'. You can override using the `.groups` argument.
```

``` r
cntry_ach_women_frame<-
  diff_dat %>%
  filter(gndr.c==-0.5) %>%
  group_by(cntry,essround) %>%
  summarise('ach M' = weighted.mean(x=ach.z,w=pspwght),
            'ach SD' = sqrt(wtd.var(ach.z,pspwght))) %>%
  group_by(cntry) %>%
  summarise('ach M Women' = mean(x=`ach M`),
            'ach SD Women'= mean(x=`ach SD`))
```

```
## `summarise()` has grouped output by 'cntry'. You can override using the `.groups` argument.
```

``` r
cntry_ach_men_frame<-
  diff_dat %>%
  filter(gndr.c==0.5) %>%
  group_by(cntry,essround) %>%
  summarise('ach M' = weighted.mean(x=ach.z,w=pspwght),
            'ach SD' = sqrt(wtd.var(ach.z,pspwght))) %>%
  group_by(cntry) %>%
  summarise('ach M Men' = mean(x=`ach M`),
            'ach SD Men'= mean(x=`ach SD`))
```

```
## `summarise()` has grouped output by 'cntry'. You can override using the `.groups` argument.
```

``` r
# link n and ach datasets

desc_frame<-
  left_join(
    x=cntry_n_frame,
    y=cntry_ach_frame,
    by="cntry"
  )

desc_frame<-
  left_join(
    x=desc_frame,
    y=cntry_ach_women_frame,
    by="cntry"
  )

desc_frame<-
  left_join(
    x=desc_frame,
    y=cntry_ach_men_frame,
    by="cntry"
  )

# Add country-specific differences
desc_frame$D<-desc_frame$`ach M Men`-desc_frame$`ach M Women`

desc_frame
```

```
## # A tibble: 33 × 10
##    cntry `n ESS rounds`     n `ach M` `ach SD` `ach M Women` `ach SD Women` `ach M Men` `ach SD Men`
##    <chr>          <dbl> <dbl>   <dbl>    <dbl>         <dbl>          <dbl>       <dbl>        <dbl>
##  1 AT                 6 13077  0.209     0.999        0.0708          1.03      0.358          0.945
##  2 BE                10 17313 -0.0207    0.975       -0.0900          0.998     0.0524         0.944
##  3 BG                 6 12641  0.586     0.923        0.508           0.955     0.670          0.877
##  4 CH                10 16720 -0.0194    1.01        -0.157           1.02      0.125          0.967
##  5 CY                 5  5105  0.194     1.06         0.146           1.08      0.243          1.03 
##  6 CZ                 9 18934 -0.184     1.04        -0.329           1.07     -0.0265         0.989
##  7 DE                 9 25389 -0.149     0.999       -0.282           1.01     -0.00901        0.963
##  8 DK                 8 12198 -0.224     1.03        -0.313           1.03     -0.133          1.03 
##  9 EE                 9 16692 -0.293     1.01        -0.371           1.03     -0.199          0.988
## 10 ES                 9 16954 -0.263     1.08        -0.369           1.08     -0.152          1.06 
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
    `ach M`, `ach SD`,
    `ach M Women`, `ach SD Women`,
    `ach M Men`, `ach SD Men`,
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
##    Country    `n ESS rounds`     n `ach M` `ach SD` `ach M Women` `ach SD Women` `ach M Men` `ach SD Men`
##    <chr>               <dbl> <dbl> <chr>   <chr>    <chr>         <chr>          <chr>       <chr>       
##  1 Austria                 6 13077 0.21    1.00     0.07          1.03           0.36        0.95        
##  2 Belgium                10 17313 -0.02   0.98     -0.09         1.00           0.05        0.94        
##  3 Bulgaria                6 12641 0.59    0.92     0.51          0.95           0.67        0.88        
##  4 Switzerla…             10 16720 -0.02   1.01     -0.16         1.02           0.13        0.97        
##  5 Cyprus                  5  5105 0.19    1.06     0.15          1.08           0.24        1.03        
##  6 Czechia                 9 18934 -0.18   1.04     -0.33         1.07           -0.03       0.99        
##  7 Germany                 9 25389 -0.15   1.00     -0.28         1.01           -0.01       0.96        
##  8 Denmark                 8 12198 -0.22   1.03     -0.31         1.03           -0.13       1.03        
##  9 Estonia                 9 16692 -0.29   1.01     -0.37         1.03           -0.20       0.99        
## 10 Spain                   9 16954 -0.26   1.08     -0.37         1.08           -0.15       1.06        
## 11 Finland                10 18050 -0.49   1.03     -0.65         1.04           -0.33       0.99        
## 12 France                 10 18720 -0.58   1.09     -0.68         1.07           -0.46       1.09        
## 13 UK                     10 21456 -0.07   1.06     -0.19         1.08           0.06        1.03        
## 14 Greece                  5 12464 0.35    0.99     0.27          1.01           0.44        0.96        
## 15 Croatia                 4  6368 0.05    1.08     -0.03         1.09           0.13        1.05        
## 16 Hungary                10 16006 0.30    0.99     0.22          1.01           0.38        0.97        
## 17 Ireland                10 20576 0.08    1.05     0.01          1.08           0.15        1.01        
## 18 Israel                  6 13964 0.74    0.91     0.72          0.92           0.76        0.90        
## 19 Iceland                 5  3832 -0.57   1.05     -0.61         1.05           -0.53       1.04        
## 20 Italy                   4  8663 0.45    0.90     0.38          0.94           0.53        0.85        
## 21 Lithuania               6 11714 0.06    1.01     0.02          1.02           0.12        0.98        
## 22 Latvia                  2  2866 0.13    1.02     0.05          1.03           0.22        1.00        
## 23 Montenegro              2  2441 0.32    0.99     0.28          1.01           0.37        0.96        
## 24 Netherlan…             10 18048 -0.10   0.94     -0.22         0.98           0.02        0.89        
## 25 Norway                 10 15186 -0.35   0.98     -0.43         1.01           -0.28       0.95        
## 26 Poland                  9 15314 0.08    0.97     -0.03         1.00           0.19        0.92        
## 27 Portugal               10 17705 0.09    0.91     0.01          0.92           0.18        0.88        
## 28 Russia                  5 12139 0.23    1.05     0.16          1.07           0.31        1.01        
## 29 Sweden                  9 14897 -0.49   1.02     -0.56         1.02           -0.42       1.01        
## 30 Slovenia               10 13238 0.46    0.88     0.40          0.92           0.53        0.84        
## 31 Slovakia                7 11132 0.13    0.95     0.02          0.98           0.25        0.90        
## 32 Turkey                  2  4108 0.68    0.84     0.59          0.88           0.78        0.79        
## 33 Ukraine                 5  9454 -0.14   1.08     -0.22         1.10           -0.04       1.05        
## # ℹ 5 more variables: D <chr>, GEI <chr>, GGGI <chr>, GDI <chr>, GDP <chr>
```

``` r
export(cntry_desc_tbl,"../results/ach/cntry_desc_tbl.xlsx",overwrite=T)
```

## Country-level correlations table


``` r
cor_frame<-
  desc_frame %>%
  select(
    VBMT=`ach M`,
    VBMT_Women=`ach M Women`,
    VBMT_Men=`ach M Men`,
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
  filename = "../results/ach/CorTable1.doc",
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
##   1. VBMT       0.04  0.34                                                                           
##                                                                                                      
##   2. VBMT_Women -0.05 0.36 1.00                                                                      
##                            [.99, 1.00]                                                               
##                                                                                                      
##   3. VBMT_Men   0.13  0.33 1.00         .98                                                          
##                            [.99, 1.00]  [.97, .99]                                                   
##                                                                                                      
##   4. D          0.18  0.07 -.34         -.42         -.25                                            
##                            [-.61, .00]  [-.67, -.09] [-.54, .11]                                     
##                                                                                                      
##   5. GEI        0.87  0.07 -.60         -.60         -.60         .20                                
##                            [-.78, -.31] [-.78, -.31] [-.78, -.31] [-.16, .51]                        
##                                                                                                      
##   6. GGGI       0.73  0.05 -.70         -.69         -.72         .07         .73                    
##                            [-.84, -.48] [-.83, -.45] [-.85, -.50] [-.28, .40] [.52, .86]             
##                                                                                                      
##   7. GDI        0.99  0.03 -.23         -.20         -.25         -.15        .07         .20        
##                            [-.53, .12]  [-.51, .15]  [-.55, .10]  [-.47, .20] [-.29, .41] [-.16, .51]
##                                                                                                      
##   8. log_GDP    10.62 0.40 -.47         -.48         -.46         .24         .75         .67        
##                            [-.70, -.15] [-.71, -.17] [-.70, -.14] [-.11, .54] [.55, .87]  [.42, .82] 
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
mod0<-lmer(ach.z~(1|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod0)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ach.z ~ (1 | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1314647.1 1314680.1 -657320.6 1314641.1    441165 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -5.7669 -0.6832  0.0143  0.6571  5.2191 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  cntry    (Intercept) 0.1164   0.3411  
##  Residual             1.0156   1.0078  
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##             Estimate Std. Error       df t value Pr(>|t|)
## (Intercept)  0.04173    0.05941 32.99176   0.702    0.487
```

``` r
getVC(mod0)
```

```
##        grp        var1 var2 sdcor vcov
## 1    cntry (Intercept) <NA>  0.34 0.12
## 2 Residual        <NA> <NA>  1.01 1.02
```

``` r
r2mlm(mod0,bargraph = F)
```

```
## $Decompositions
##                     total within between
## fixed, within   0.0000000      0      NA
## fixed, between  0.0000000     NA       0
## slope variation 0.0000000      0      NA
## mean variation  0.1028004     NA       1
## sigma2          0.8971996      1      NA
## 
## $R2s
##         total within between
## f1  0.0000000      0      NA
## f2  0.0000000     NA       0
## v   0.0000000      0      NA
## m   0.1028004     NA       1
## f   0.0000000     NA      NA
## fv  0.0000000      0      NA
## fvm 0.1028004     NA      NA
```

## mod1: Gender fixed effect


``` r
mod1<-lmer(ach.z~gndr.c+(1|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod1)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ach.z ~ gndr.c + (1 | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1310442.8 1310486.8 -655217.4 1310434.8    441164 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.0345 -0.6688  0.0343  0.6622  5.0339 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  cntry    (Intercept) 0.1165   0.3413  
##  Residual             1.0059   1.0030  
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##              Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept) 4.547e-02  5.943e-02 3.299e+01   0.765     0.45    
## gndr.c      1.961e-01  3.016e-03 4.411e+05  65.011   <2e-16 ***
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
##              Est.    SE         df      t    p     LL    UL
## (Intercept) 0.045 0.059     32.991  0.765 0.45 -0.075 0.166
## gndr.c      0.196 0.003 441135.483 65.011 0.00  0.190 0.202
```

``` r
getVC(mod1)
```

```
##        grp        var1 var2 sdcor vcov
## 1    cntry (Intercept) <NA>  0.34 0.12
## 2 Residual        <NA> <NA>  1.00 1.01
```

``` r
r2mlm(mod1,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.008437946
## slope variation 0.000000000
## mean variation  0.102882338
## sigma2          0.888679716
## 
## $R2s
##           total
## f   0.008437946
## v   0.000000000
## m   0.102882338
## fv  0.008437946
## fvm 0.111320284
```

## mod2: Gender fixed and random effect

* Include random effect correlation by default


``` r
mod2<-lmer(ach.z~gndr.c+(gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ach.z ~ gndr.c + (gndr.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1310065.1 1310131.1 -655026.6 1310053.1    441162 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -5.9289 -0.6682  0.0305  0.6625  5.0663 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr 
##  cntry    (Intercept) 0.116123 0.3408        
##           gndr.c      0.004084 0.0639   -0.36
##  Residual             1.004883 1.0024        
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##             Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)  0.04525    0.05935 32.99371   0.763    0.451    
## gndr.c       0.18285    0.01165 33.05858  15.693   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##        (Intr)
## gndr.c -0.347
```

``` r
getFE(mod2,round=3)
```

```
##              Est.    SE     df      t     p     LL    UL
## (Intercept) 0.045 0.059 32.994  0.763 0.451 -0.075 0.166
## gndr.c      0.183 0.012 33.059 15.693 0.000  0.159 0.207
```

``` r
getVC(mod2)
```

```
##        grp        var1   var2 sdcor  vcov
## 1    cntry (Intercept)   <NA>  0.34  0.12
## 2    cntry      gndr.c   <NA>  0.06  0.00
## 3    cntry (Intercept) gndr.c -0.36 -0.01
## 4 Residual        <NA>   <NA>  1.00  1.00
```

``` r
r2mlm(mod2,bargraph = F)
```

```
## $Decompositions
##                        total
## fixed           0.0073461459
## slope variation 0.0008972411
## mean variation  0.1032278136
## sigma2          0.8885287994
## 
## $R2s
##            total
## f   0.0073461459
## v   0.0008972411
## m   0.1032278136
## fv  0.0082433870
## fvm 0.1114712006
```

``` r
anova(mod1,mod2)
```

```
## Data: diff_dat
## Models:
## mod1: ach.z ~ gndr.c + (1 | cntry)
## mod2: ach.z ~ gndr.c + (gndr.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod1    4 1310443 1310487 -655217   1310435                         
## mod2    6 1310065 1310131 -655027   1310053 381.62  2  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

``` r
lvl2_var_cond_lvl1(mod2,lvl1.var = "gndr.c",lvl1.values = c(0.5,-0.5))
```

```
##   lvl1.value lvl2.cond.var lvl2.cond.sd
## 1        0.5     0.1092195    0.3304837
## 2       -0.5     0.1250688    0.3536507
```

* Test for random effect correlation


``` r
mod2_norecov<-lmer(ach.z~gndr.c+(gndr.c||cntry),data=diff_dat,REML=F,weights = pspwght,
                   control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2_norecov)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ach.z ~ gndr.c + (gndr.c || cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1310067.3 1310122.3 -655028.6 1310057.3    441163 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -5.9277 -0.6685  0.0303  0.6624  5.0643 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  cntry    (Intercept) 0.11616  0.34082 
##  cntry.1  gndr.c      0.00413  0.06427 
##  Residual             1.00488  1.00244 
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##             Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)  0.04524    0.05936 32.99189   0.762    0.451    
## gndr.c       0.18304    0.01172 33.00653  15.622   <2e-16 ***
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
## (Intercept) 0.045 0.059 32.992  0.762 0.451 -0.076 0.166
## gndr.c      0.183 0.012 33.007 15.622 0.000  0.159 0.207
```

``` r
getVC(mod2_norecov)
```

```
##        grp        var1 var2 sdcor vcov
## 1    cntry (Intercept) <NA>  0.34 0.12
## 2  cntry.1      gndr.c <NA>  0.06 0.00
## 3 Residual        <NA> <NA>  1.00 1.00
```

``` r
anova(mod2_norecov,mod2)
```

```
## Data: diff_dat
## Models:
## mod2_norecov: ach.z ~ gndr.c + (gndr.c || cntry)
## mod2: ach.z ~ gndr.c + (gndr.c | cntry)
##              npar     AIC     BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)  
## mod2_norecov    5 1310067 1310122 -655029   1310057                       
## mod2            6 1310065 1310131 -655027   1310053 4.1288  1    0.04216 *
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```


## mod2 with Gender-equality index (GEI)


``` r
mod2_GEI<-lmer(ach.z~gndr.c+gei.z.cm+gndr.c:gei.z.cm+
                 (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2_GEI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ach.z ~ gndr.c + gei.z.cm + gndr.c:gei.z.cm + (gndr.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1278463.0 1278550.8 -639223.5 1278447.0    431770 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -5.9406 -0.6690  0.0315  0.6643  5.0762 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr 
##  cntry    (Intercept) 0.07407  0.27216       
##           gndr.c      0.00404  0.06356  -0.31
##  Residual             1.00091  1.00046       
## Number of obs: 431778, groups:  cntry, 32
## 
## Fixed effects:
##                 Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)      0.05066    0.04815 32.00972   1.052 0.300557    
## gndr.c           0.18304    0.01178 31.88645  15.542  < 2e-16 ***
## gei.z.cm        -0.21501    0.04893 32.04749  -4.394 0.000114 ***
## gndr.c:gei.z.cm  0.01358    0.01219 34.26986   1.114 0.272928    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c g.z.cm
## gndr.c      -0.294              
## gei.z.cm     0.000  0.000       
## gndr.c:g.z.  0.000 -0.028 -0.289
```

``` r
getFE(mod2_GEI,round=3)
```

```
##                   Est.    SE     df      t     p     LL     UL
## (Intercept)      0.051 0.048 32.010  1.052 0.301 -0.047  0.149
## gndr.c           0.183 0.012 31.886 15.542 0.000  0.159  0.207
## gei.z.cm        -0.215 0.049 32.047 -4.394 0.000 -0.315 -0.115
## gndr.c:gei.z.cm  0.014 0.012 34.270  1.114 0.273 -0.011  0.038
```

``` r
getVC(mod2_GEI)
```

```
##        grp        var1   var2 sdcor  vcov
## 1    cntry (Intercept)   <NA>  0.27  0.07
## 2    cntry      gndr.c   <NA>  0.06  0.00
## 3    cntry (Intercept) gndr.c -0.31 -0.01
## 4 Residual        <NA>   <NA>  1.00  1.00
```

``` r
r2mlm(mod2_GEI,bargraph = F)
```

```
## $Decompositions
##                        total
## fixed           0.0355601595
## slope variation 0.0008999591
## mean variation  0.0667264284
## sigma2          0.8968134531
## 
## $R2s
##            total
## f   0.0355601595
## v   0.0008999591
## m   0.0667264284
## fv  0.0364601186
## fvm 0.1031865469
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
## Time difference of 28.66401 secs
```

``` r
round(ddsc_mod2_GEI$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.125        0.354        1.005     1.130 0.111   7204.152 0.999   0.999
## 2        0.5         0.109        0.330        1.005     1.114 0.098   6164.576 0.998   0.999
```

``` r
round(ddsc_mod2_GEI$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gei.z.cm
## means_y1           0.141 0.341    1.000           1.000    0.983           0.983   -0.613
## means_y1_scaled    0.401 0.968    1.000           1.000    0.983           0.983   -0.613
## means_y2          -0.040 0.363    0.983           0.983    1.000           1.000   -0.609
## means_y2_scaled   -0.114 1.031    0.983           0.983    1.000           1.000   -0.609
## gei.z.cm           0.000 1.000   -0.613          -0.613   -0.609          -0.609    1.000
## gei.z.cm_scaled    0.000 1.000   -0.613          -0.613   -0.609          -0.609    1.000
## diff_score         0.181 0.069   -0.233          -0.233   -0.410          -0.410    0.181
## diff_score_scaled  0.515 0.197   -0.233          -0.233   -0.410          -0.410    0.181
##                   gei.z.cm_scaled diff_score diff_score_scaled
## means_y1                   -0.613     -0.233            -0.233
## means_y1_scaled            -0.613     -0.233            -0.233
## means_y2                   -0.609     -0.410            -0.410
## means_y2_scaled            -0.609     -0.410            -0.410
## gei.z.cm                    1.000      0.181             0.181
## gei.z.cm_scaled             1.000      0.181             0.181
## diff_score                  0.181      1.000             1.000
## diff_score_scaled           0.181      1.000             1.000
```

``` r
round(ddsc_mod2_GEI$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                      -0.196 0.176 34.270  -1.114   0.273   -0.552    0.161
## w_11                         -0.222 0.051 32.088  -4.347   0.000   -0.326   -0.118
## w_21                         -0.208 0.048 32.091  -4.381   0.000   -0.305   -0.111
## r_xy1                        -0.651 0.150 32.088  -4.347   0.000   -0.956   -0.346
## r_xy2                        -0.573 0.131 32.091  -4.381   0.000   -0.840   -0.307
## b_11                         -0.630 0.145 32.088  -4.347   0.000   -0.926   -0.335
## b_21                         -0.592 0.135 32.091  -4.381   0.000   -0.867   -0.317
## main_effect                  -0.215 0.049 32.047  -4.394   0.000   -0.315   -0.115
## moderator_effect              0.183 0.012 31.886  15.542   0.000    0.159    0.207
## interaction                   0.014 0.012 34.270   1.114   0.273   -0.011    0.038
## q_b11_b21                    -0.062    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                  -0.125    NA     NA      NA      NA       NA       NA
## cross_over_point            -13.474    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.201 0.047 32.204  -4.296   0.000   -0.297   -0.106
## interaction_vs_main_bscale   -0.572 0.133 32.204  -4.296   0.000   -0.844   -0.301
## interaction_vs_main_rscale   -0.534 0.125 32.219  -4.283   0.000   -0.789   -0.280
## dadas                        -0.416 0.095 32.091  -4.381   1.000   -0.610   -0.223
## dadas_bscale                 -1.183 0.270 32.091  -4.381   1.000   -1.733   -0.633
## dadas_rscale                 -1.147 0.262 32.091  -4.381   1.000   -1.680   -0.614
## abs_diff                      0.014 0.012 34.270   1.114   0.136   -0.011    0.038
## abs_sum                       0.430 0.098 32.047   4.394   0.000    0.231    0.629
## abs_diff_bscale               0.039 0.035 34.270   1.114   0.136   -0.032    0.109
## abs_sum_bscale                1.222 0.278 32.047   4.394   0.000    0.656    1.788
## abs_diff_rscale               0.078 0.038 34.288   2.034   0.025    0.000    0.155
## abs_sum_rscale                1.224 0.279 32.047   4.394   0.000    0.657    1.792
```

``` r
round(ddsc_mod2_GEI$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.008 -0.364  4.129  1.000  0.042
```

``` r
d_GEI<-ddsc_mod2_GEI$ddsc_sem_fit$data

ddsc_sem_GEI<-
  ddsc_sem(data=d_GEI,x = "gei.z.cm",y1="means_y2",
           y2="means_y1",q_sesoi = .10)
round(ddsc_sem_GEI$results,3)
```

```
##                                     est     se       z pvalue ci.lower ci.upper
## r_xy1_y2                         -0.181  0.174  -1.039  0.299   -0.521    0.160
## r_xy1                            -0.609  0.140  -4.347  0.000   -0.884   -0.335
## r_xy2                            -0.613  0.140  -4.386  0.000   -0.887   -0.339
## b_11                             -0.628  0.145  -4.347  0.000   -0.912   -0.345
## b_21                             -0.593  0.135  -4.386  0.000   -0.858   -0.328
## b_10                             -0.114  0.142  -0.798  0.425   -0.393    0.165
## b_20                              0.401  0.133   3.016  0.003    0.140    0.662
## res_cov_y1_y2                     0.589  0.149   3.943  0.000    0.296    0.882
## diff_b10_b20                     -0.515  0.034 -15.257  0.000   -0.581   -0.449
## diff_b11_b21                     -0.036  0.034  -1.039  0.299   -0.103    0.032
## diff_rxy1_rxy2                    0.003  0.033   0.104  0.917   -0.061    0.068
## q_b11_b21                        -0.057  0.061  -0.934  0.350   -0.176    0.062
## q_rxy1_rxy2                       0.005  0.053   0.104  0.917   -0.098    0.109
## cross_over_point                -14.458 13.952  -1.036  0.300  -41.802   12.887
## sum_b11_b21                      -1.221  0.278  -4.396  0.000   -1.766   -0.677
## main_effect                      -0.611  0.139  -4.396  0.000   -0.883   -0.338
## interaction_vs_main_effect       -0.575  0.134  -4.306  0.000   -0.837   -0.313
## diff_abs_b11_abs_b21              0.036  0.034   1.039  0.299   -0.032    0.103
## abs_diff_b11_b21                  0.036  0.034   1.039  0.149   -0.032    0.103
## abs_sum_b11_b21                   1.221  0.278   4.396  0.000    0.677    1.766
## dadas                            -1.186  0.270  -4.386  1.000   -1.716   -0.656
## q_r_equivalence                  -0.095  0.053  -1.795  0.036       NA       NA
## q_b_equivalence                  -0.043  0.061  -0.710  0.239       NA       NA
## cross_over_point_equivalence     14.458 13.952   1.036  0.850       NA       NA
## cross_over_point_minimal_effect  14.458 13.952   1.036  0.150       NA       NA
```

``` r
round(ddsc_sem_GEI$variance_test,3)
```

```
##             est    se       z pvalue ci.lower ci.upper
## cov_y1y2  0.950 0.240   3.965   0.00    0.480    1.420
## var_y1    1.031 0.258   4.000   0.00    0.526    1.536
## var_y2    0.907 0.227   4.000   0.00    0.463    1.351
## var_diff  0.124 0.071   1.751   0.08   -0.015    0.262
## var_ratio 1.136 0.075  15.214   0.00    0.990    1.283
## cor_y1y2  0.983 0.006 160.813   0.00    0.971    0.995
```

### Figure 


``` r
# plot the results

# start with obtaining predicted values for means and differences

# refit reduced and full models with GGGI in original scale


mod2_GEI_unstd<-lmer(ach.z~gndr.c+gei.cm+gndr.c:gei.cm+
                 (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

mod2_GEI_unstd_red<-lmer(ach.z~gndr.c+
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


p1.ach.flags<-
  ggplot(p,aes(y=yvar,x=gei.cm,color=gndr.c))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("Gender Equality Index (Average 2002-2022)")+
  #ylim(c(min.y.pred,max.y.pred))+
  #xlim(c(0.60,1.00))+
  ylab("Value achievement (Average 2002-2022)")+
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

p2.ach.flags<-ggplot(p2,aes(y=yvar,x=gei.cm))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("Gender Equality Index (Average 2002-2022)")+
  #xlim(c(0.60,1.00))+
  #ylim(c(min.y.narrow,max.y.narrow))+
  ylab("Difference in value achievement")+
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
  ggarrange(p1.ach.flags,p2.ach.flags,align = "v",
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

![](Analysis_achievement_files/figure-html/unnamed-chunk-18-1.png)<!-- -->

``` r
png(filename = 
      "../results/ach/GEI_flags.png",
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
mod2_GGGI<-lmer(ach.z~gndr.c+gggi.z.cm+gndr.c:gggi.z.cm+
                 (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2_GGGI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ach.z ~ gndr.c + gggi.z.cm + gndr.c:gggi.z.cm + (gndr.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##  932763.8  932849.1 -466373.9  932747.8    314638 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.0173 -0.6660  0.0300  0.6645  5.0749 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr 
##  cntry    (Intercept) 0.065371 0.25568       
##           gndr.c      0.003817 0.06178  -0.40
##  Residual             0.998156 0.99908       
## Number of obs: 314646, groups:  cntry, 33
## 
## Fixed effects:
##                  Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)       0.07369    0.04456 33.00420   1.654    0.108    
## gndr.c            0.17759    0.01152 31.39707  15.421 3.35e-16 ***
## gggi.z.cm        -0.24961    0.04528 33.07068  -5.513 4.03e-06 ***
## gndr.c:gggi.z.cm  0.01276    0.01201 34.55486   1.062    0.296    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c ggg.z.
## gndr.c      -0.376              
## gggi.z.cm    0.000  0.000       
## gndr.c:gg..  0.000 -0.017 -0.366
```

``` r
getFE(mod2_GGGI,round=3)
```

```
##                    Est.    SE     df      t     p     LL     UL
## (Intercept)       0.074 0.045 33.004  1.654 0.108 -0.017  0.164
## gndr.c            0.178 0.012 31.397 15.421 0.000  0.154  0.201
## gggi.z.cm        -0.250 0.045 33.071 -5.513 0.000 -0.342 -0.158
## gndr.c:gggi.z.cm  0.013 0.012 34.555  1.062 0.296 -0.012  0.037
```

``` r
getVC(mod2_GGGI)
```

```
##        grp        var1   var2 sdcor  vcov
## 1    cntry (Intercept)   <NA>  0.26  0.07
## 2    cntry      gndr.c   <NA>  0.06  0.00
## 3    cntry (Intercept) gndr.c -0.40 -0.01
## 4 Residual        <NA>   <NA>  1.00  1.00
```

``` r
r2mlm(mod2_GGGI,bargraph = F)
```

```
## $Decompositions
##                        total
## fixed           0.0479224009
## slope variation 0.0008475555
## mean variation  0.0589025782
## sigma2          0.8923274653
## 
## $R2s
##            total
## f   0.0479224009
## v   0.0008475555
## m   0.0589025782
## fv  0.0487699564
## fvm 0.1076725347
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
## Time difference of 28.55501 secs
```

``` r
round(ddsc_mod2_GGGI$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.125        0.354        1.005     1.130 0.111   7204.152 0.999   0.999
## 2        0.5         0.109        0.330        1.005     1.114 0.098   6164.576 0.998   0.999
```

``` r
round(ddsc_mod2_GGGI$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gggi.z.cm
## means_y1           0.161 0.349    1.000           1.000    0.983           0.983    -0.699
## means_y1_scaled    0.445 0.964    1.000           1.000    0.983           0.983    -0.699
## means_y2          -0.014 0.375    0.983           0.983    1.000           1.000    -0.680
## means_y2_scaled   -0.038 1.034    0.983           0.983    1.000           1.000    -0.680
## gggi.z.cm          0.000 1.000   -0.699          -0.699   -0.680          -0.680     1.000
## gggi.z.cm_scaled   0.000 1.000   -0.699          -0.699   -0.680          -0.680     1.000
## diff_score         0.175 0.072   -0.264          -0.264   -0.437          -0.437     0.149
## diff_score_scaled  0.483 0.198   -0.264          -0.264   -0.437          -0.437     0.149
##                   gggi.z.cm_scaled diff_score diff_score_scaled
## means_y1                    -0.699     -0.264            -0.264
## means_y1_scaled             -0.699     -0.264            -0.264
## means_y2                    -0.680     -0.437            -0.437
## means_y2_scaled             -0.680     -0.437            -0.437
## gggi.z.cm                    1.000      0.149             0.149
## gggi.z.cm_scaled             1.000      0.149             0.149
## diff_score                   0.149      1.000             1.000
## diff_score_scaled            0.149      1.000             1.000
```

``` r
round(ddsc_mod2_GGGI$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                      -0.178 0.168 34.555  -1.062   0.296   -0.518    0.162
## w_11                         -0.256 0.048 33.133  -5.355   0.000   -0.353   -0.159
## w_21                         -0.243 0.043 33.130  -5.600   0.000   -0.332   -0.155
## r_xy1                        -0.733 0.137 33.133  -5.355   0.000   -1.011   -0.454
## r_xy2                        -0.649 0.116 33.130  -5.600   0.000   -0.885   -0.413
## b_11                         -0.707 0.132 33.133  -5.355   0.000   -0.976   -0.439
## b_21                         -0.672 0.120 33.130  -5.600   0.000   -0.916   -0.428
## main_effect                  -0.250 0.045 33.071  -5.513   0.000   -0.342   -0.158
## moderator_effect              0.178 0.012 31.397  15.421   0.000    0.154    0.201
## interaction                   0.013 0.012 34.555   1.062   0.296   -0.012    0.037
## q_b11_b21                    -0.067    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                  -0.161    NA     NA      NA      NA       NA       NA
## cross_over_point            -13.920    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.237 0.042 33.319  -5.590   0.000   -0.323   -0.151
## interaction_vs_main_bscale   -0.654 0.117 33.319  -5.590   0.000   -0.892   -0.416
## interaction_vs_main_rscale   -0.607 0.109 33.346  -5.582   0.000   -0.829   -0.386
## dadas                        -0.486 0.087 33.130  -5.600   1.000   -0.663   -0.310
## dadas_bscale                 -1.344 0.240 33.130  -5.600   1.000   -1.832   -0.856
## dadas_rscale                 -1.298 0.232 33.130  -5.600   1.000   -1.770   -0.827
## abs_diff                      0.013 0.012 34.555   1.062   0.148   -0.012    0.037
## abs_sum                       0.499 0.091 33.071   5.513   0.000    0.315    0.683
## abs_diff_bscale               0.035 0.033 34.555   1.062   0.148   -0.032    0.103
## abs_sum_bscale                1.379 0.250 33.071   5.513   0.000    0.870    1.888
## abs_diff_rscale               0.084 0.037 34.729   2.240   0.016    0.008    0.159
## abs_sum_rscale                1.382 0.251 33.071   5.509   0.000    0.872    1.892
```

``` r
round(ddsc_mod2_GGGI$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.008 -0.364  4.129  1.000  0.042
```

``` r
d_GGGI<-ddsc_mod2_GGGI$ddsc_sem_fit$data

ddsc_sem_GGGI<-
  ddsc_sem(data=d_GGGI,x = "gggi.z.cm",y1="means_y2",
           y2="means_y1",q_sesoi = .10)
round(ddsc_sem_GGGI$results,3)
```

```
##                                     est     se       z pvalue ci.lower ci.upper
## r_xy1_y2                         -0.149  0.172  -0.865  0.387   -0.486    0.189
## r_xy1                            -0.680  0.128  -5.335  0.000   -0.930   -0.430
## r_xy2                            -0.699  0.124  -5.620  0.000   -0.943   -0.455
## b_11                             -0.704  0.132  -5.335  0.000   -0.962   -0.445
## b_21                             -0.674  0.120  -5.620  0.000   -0.910   -0.439
## b_10                             -0.038  0.130  -0.290  0.772   -0.292    0.217
## b_20                              0.445  0.118   3.767  0.000    0.214    0.677
## res_cov_y1_y2                     0.490  0.123   3.995  0.000    0.250    0.731
## diff_b10_b20                     -0.483  0.034 -14.392  0.000   -0.549   -0.417
## diff_b11_b21                     -0.029  0.034  -0.865  0.387   -0.096    0.037
## diff_rxy1_rxy2                    0.019  0.032   0.587  0.557   -0.044    0.082
## q_b11_b21                        -0.056  0.074  -0.763  0.445   -0.200    0.088
## q_rxy1_rxy2                       0.036  0.061   0.587  0.557   -0.084    0.156
## cross_over_point                -16.392 18.992  -0.863  0.388  -53.615   20.832
## sum_b11_b21                      -1.378  0.250  -5.515  0.000   -1.868   -0.888
## main_effect                      -0.689  0.125  -5.515  0.000   -0.934   -0.444
## interaction_vs_main_effect       -0.660  0.117  -5.623  0.000   -0.890   -0.430
## diff_abs_b11_abs_b21              0.029  0.034   0.865  0.387   -0.037    0.096
## abs_diff_b11_b21                  0.029  0.034   0.865  0.194   -0.037    0.096
## abs_sum_b11_b21                   1.378  0.250   5.515  0.000    0.888    1.868
## dadas                            -1.349  0.240  -5.620  1.000   -1.819   -0.878
## q_r_equivalence                  -0.064  0.061  -1.045  0.148       NA       NA
## q_b_equivalence                  -0.044  0.074  -0.596  0.276       NA       NA
## cross_over_point_equivalence     16.392 18.992   0.863  0.806       NA       NA
## cross_over_point_minimal_effect  16.392 18.992   0.863  0.194       NA       NA
```

``` r
round(ddsc_sem_GGGI$variance_test,3)
```

```
##             est    se       z pvalue ci.lower ci.upper
## cov_y1y2  0.951 0.236   4.027  0.000    0.488    1.413
## var_y1    1.038 0.255   4.062  0.000    0.537    1.538
## var_y2    0.902 0.222   4.062  0.000    0.467    1.337
## var_diff  0.136 0.071   1.923  0.055   -0.003    0.274
## var_ratio 1.150 0.074  15.564  0.000    1.006    1.295
## cor_y1y2  0.983 0.006 165.767  0.000    0.971    0.994
```

### Figure


``` r
# plot the results

# start with obtaining predicted values for means and differences

# refit reduced and full models with GGGI in original scale


mod2_GGGI_unstd<-lmer(ach.z~gndr.c+gggi.cm+gndr.c:gggi.cm+
                       (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
                     control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

mod2_GGGI_unstd_red<-lmer(ach.z~gndr.c+
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


p1.ach.flags<-
  ggplot(p,aes(y=yvar,x=gggi.cm,color=gndr.c))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("Global Gender Gap Index (Average 2002-2022)")+
  #ylim(c(min.y.pred,max.y.pred))+
  #xlim(c(0.60,1.00))+
  ylab("Value achievement (Average 2002-2022)")+
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

p2.ach.flags<-ggplot(p2,aes(y=yvar,x=gggi.cm))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("Global Gender Gap Index (Average 2002-2022)")+
  #xlim(c(0.60,1.00))+
  #ylim(c(min.y.narrow,max.y.narrow))+
  ylab("Difference in value achievement")+
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
  ggarrange(p1.ach.flags,p2.ach.flags,align = "v",
            ncol=1,nrow=2,heights=c(2,1))
```

```
## Warning in geom_text(inherit.aes = F, aes(x = 0.61, y = -0.45, label = coef_q, : All aesthetics have length 1, but the data has 502 rows.
## ℹ Please consider using `annotate()` or provide this layer with data containing a single row.
```

``` r
pflag_comb
```

![](Analysis_achievement_files/figure-html/unnamed-chunk-21-1.png)<!-- -->

``` r
png(filename = 
      "../results/ach/GGGI_flags_new.png",
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
mod2_GDI<-lmer(ach.z~gndr.c+gdi.z.cm+gndr.c:gdi.z.cm+
                 (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2_GDI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ach.z ~ gndr.c + gdi.z.cm + gndr.c:gdi.z.cm + (gndr.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##   1310066   1310154   -655025   1310050    441160 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -5.9262 -0.6682  0.0303  0.6627  5.0665 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr 
##  cntry    (Intercept) 0.111424 0.3338        
##           gndr.c      0.003944 0.0628   -0.41
##  Residual             1.004884 1.0024        
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##                  Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)      0.045240   0.058135 32.995690   0.778    0.442    
## gndr.c           0.183025   0.011468 32.730584  15.959   <2e-16 ***
## gdi.z.cm        -0.069632   0.059053 33.032367  -1.179    0.247    
## gndr.c:gdi.z.cm -0.009978   0.011950 36.031501  -0.835    0.409    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c gd.z.c
## gndr.c      -0.388              
## gdi.z.cm     0.000  0.000       
## gndr.c:gd..  0.000 -0.008 -0.378
```

``` r
getFE(mod2_GDI,round=3)
```

```
##                   Est.    SE     df      t     p     LL    UL
## (Intercept)      0.045 0.058 32.996  0.778 0.442 -0.073 0.164
## gndr.c           0.183 0.011 32.731 15.959 0.000  0.160 0.206
## gdi.z.cm        -0.070 0.059 33.032 -1.179 0.247 -0.190 0.051
## gndr.c:gdi.z.cm -0.010 0.012 36.032 -0.835 0.409 -0.034 0.014
```

``` r
getVC(mod2_GDI)
```

```
##        grp        var1   var2 sdcor  vcov
## 1    cntry (Intercept)   <NA>  0.33  0.11
## 2    cntry      gndr.c   <NA>  0.06  0.00
## 3    cntry (Intercept) gndr.c -0.41 -0.01
## 4 Residual        <NA>   <NA>  1.00  1.00
```

``` r
r2mlm(mod2_GDI,bargraph = F)
```

```
## $Decompositions
##                        total
## fixed           0.0106740568
## slope variation 0.0008673772
## mean variation  0.0991973474
## sigma2          0.8892612185
## 
## $R2s
##            total
## f   0.0106740568
## v   0.0008673772
## m   0.0991973474
## fv  0.0115414341
## fvm 0.1107387815
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
## Time difference of 28.57527 secs
```

``` r
round(ddsc_mod2_GDI$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.125        0.354        1.005     1.130 0.111   7204.152 0.999   0.999
## 2        0.5         0.109        0.330        1.005     1.114 0.098   6164.576 0.998   0.999
```

``` r
round(ddsc_mod2_GDI$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gdi.z.cm
## means_y1           0.136 0.337    1.000           1.000    0.983           0.983   -0.219
## means_y1_scaled    0.390 0.968    1.000           1.000    0.983           0.983   -0.219
## means_y2          -0.045 0.359    0.983           0.983    1.000           1.000   -0.183
## means_y2_scaled   -0.130 1.031    0.983           0.983    1.000           1.000   -0.183
## gdi.z.cm           0.000 1.000   -0.219          -0.219   -0.183          -0.183    1.000
## gdi.z.cm_scaled    0.000 1.000   -0.219          -0.219   -0.183          -0.183    1.000
## diff_score         0.181 0.068   -0.231          -0.231   -0.407          -0.407   -0.119
## diff_score_scaled  0.521 0.196   -0.231          -0.231   -0.407          -0.407   -0.119
##                   gdi.z.cm_scaled diff_score diff_score_scaled
## means_y1                   -0.219     -0.231            -0.231
## means_y1_scaled            -0.219     -0.231            -0.231
## means_y2                   -0.183     -0.407            -0.407
## means_y2_scaled            -0.183     -0.407            -0.407
## gdi.z.cm                    1.000     -0.119            -0.119
## gdi.z.cm_scaled             1.000     -0.119            -0.119
## diff_score                 -0.119      1.000             1.000
## diff_score_scaled          -0.119      1.000             1.000
```

``` r
round(ddsc_mod2_GDI$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                       0.146 0.175 36.032   0.835   0.409   -0.209    0.501
## w_11                         -0.065 0.062 33.056  -1.050   0.301   -0.190    0.061
## w_21                         -0.075 0.057 33.064  -1.308   0.200   -0.191    0.041
## r_xy1                        -0.192 0.183 33.056  -1.050   0.301   -0.564    0.180
## r_xy2                        -0.208 0.159 33.064  -1.308   0.200   -0.532    0.116
## b_11                         -0.186 0.177 33.056  -1.050   0.301   -0.546    0.174
## b_21                         -0.215 0.164 33.064  -1.308   0.200   -0.548    0.119
## main_effect                  -0.070 0.059 33.032  -1.179   0.247   -0.190    0.051
## moderator_effect              0.183 0.011 32.731  15.959   0.000    0.160    0.206
## interaction                  -0.010 0.012 36.032  -0.835   0.409   -0.034    0.014
## q_b11_b21                     0.030    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.017    NA     NA      NA      NA       NA       NA
## cross_over_point             18.342    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.060 0.065 33.133  -0.925   0.362   -0.191    0.072
## interaction_vs_main_bscale   -0.172 0.186 33.133  -0.925   0.362   -0.549    0.206
## interaction_vs_main_rscale   -0.184 0.197 33.124  -0.935   0.356   -0.584    0.216
## dadas                        -0.129 0.123 33.056  -1.050   0.849   -0.380    0.121
## dadas_bscale                 -0.372 0.354 33.056  -1.050   0.849   -1.092    0.348
## dadas_rscale                 -0.384 0.366 33.056  -1.050   0.849   -1.128    0.360
## abs_diff                      0.010 0.012 36.032   0.835   0.205   -0.014    0.034
## abs_sum                       0.139 0.118 33.032   1.179   0.123   -0.101    0.380
## abs_diff_bscale               0.029 0.034 36.032   0.835   0.205   -0.041    0.098
## abs_sum_bscale                0.400 0.340 33.032   1.179   0.123   -0.290    1.091
## abs_diff_rscale               0.016 0.040 35.730   0.404   0.344   -0.065    0.097
## abs_sum_rscale                0.400 0.340 33.032   1.175   0.124   -0.293    1.092
```

``` r
round(ddsc_mod2_GDI$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.008 -0.364  4.129  1.000  0.042
```

``` r
d_GDI<-ddsc_mod2_GDI$ddsc_sem_fit$data

ddsc_sem_GDI<-
  ddsc_sem(data=d_GDI,x = "gdi.z.cm",y1="means_y2",
           y2="means_y1",q_sesoi = .10)
round(ddsc_sem_GDI$results,3)
```

```
##                                    est     se       z pvalue ci.lower ci.upper
## r_xy1_y2                         0.119  0.173   0.688  0.492   -0.220    0.458
## r_xy1                           -0.183  0.171  -1.067  0.286   -0.518    0.153
## r_xy2                           -0.219  0.170  -1.288  0.198   -0.552    0.114
## b_11                            -0.188  0.176  -1.067  0.286   -0.534    0.158
## b_21                            -0.212  0.164  -1.288  0.198   -0.534    0.111
## b_10                            -0.130  0.174  -0.749  0.454   -0.471    0.210
## b_20                             0.390  0.162   2.412  0.016    0.073    0.708
## res_cov_y1_y2                    0.912  0.227   4.026  0.000    0.468    1.356
## diff_b10_b20                    -0.521  0.033 -15.572  0.000   -0.586   -0.455
## diff_b11_b21                     0.023  0.034   0.688  0.492   -0.043    0.090
## diff_rxy1_rxy2                   0.036  0.032   1.134  0.257   -0.026    0.098
## q_b11_b21                        0.024  0.035   0.699  0.484   -0.044    0.093
## q_rxy1_rxy2                      0.038  0.033   1.134  0.257   -0.027    0.103
## cross_over_point                22.294 32.443   0.687  0.492  -41.294   85.881
## sum_b11_b21                     -0.400  0.339  -1.179  0.239   -1.065    0.265
## main_effect                     -0.200  0.170  -1.179  0.239   -0.533    0.133
## interaction_vs_main_effect      -0.177  0.185  -0.957  0.338   -0.538    0.185
## diff_abs_b11_abs_b21            -0.023  0.034  -0.688  0.492   -0.090    0.043
## abs_diff_b11_b21                 0.023  0.034   0.688  0.246   -0.043    0.090
## abs_sum_b11_b21                  0.400  0.339   1.179  0.119   -0.265    1.065
## dadas                           -0.377  0.353  -1.067  0.857   -1.068    0.315
## q_r_equivalence                 -0.062  0.033  -1.884  0.030       NA       NA
## q_b_equivalence                 -0.076  0.035  -2.175  0.015       NA       NA
## cross_over_point_equivalence    22.294 32.443   0.687  0.754       NA       NA
## cross_over_point_minimal_effect 22.294 32.443   0.687  0.246       NA       NA
```

``` r
round(ddsc_sem_GDI$variance_test,3)
```

```
##             est    se       z pvalue ci.lower ci.upper
## cov_y1y2  0.951 0.236   4.026  0.000    0.488    1.414
## var_y1    1.031 0.254   4.062  0.000    0.533    1.528
## var_y2    0.908 0.224   4.062  0.000    0.470    1.347
## var_diff  0.122 0.069   1.766  0.077   -0.013    0.258
## var_ratio 1.135 0.073  15.496  0.000    0.991    1.278
## cor_y1y2  0.983 0.006 164.296  0.000    0.971    0.994
```

### Figure


``` r
# plot the results

# start with obtaining predicted values for means and differences


mod2_GDI_unstd<-lmer(ach.z~gndr.c+gdi.cm+gndr.c:gdi.cm+
                       (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
                     control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

mod2_GDI_unstd_red<-lmer(ach.z~gndr.c+
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


p1.ach.flags<-
  ggplot(p,aes(y=yvar,x=gdi.cm,color=gndr.c))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("Gender Development Index (Average 2002-2022)")+
  #ylim(c(min.y.pred,max.y.pred))+
  #xlim(c(0.60,1.00))+
  ylab("Value achievement (Average 2002-2022)")+
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

#p1.ach.flags


p2.ach.flags<-ggplot(p2,aes(y=yvar,x=gdi.cm))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("Gender Development Index (Average 2002-2022)")+
  #xlim(c(0.60,1.00))+
  #ylim(c(min.y.narrow,max.y.narrow))+
  ylab("Difference in value achievement")+
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

#p2.ach.flags


pflag_comb<-
  ggarrange(p1.ach.flags,p2.ach.flags,align = "v",
            ncol=1,nrow=2,heights=c(2,1))
```

```
## Warning in geom_text(inherit.aes = F, aes(x = 0.9, y = -0.45, label = coef_q, : All aesthetics have length 1, but the data has 282 rows.
## ℹ Please consider using `annotate()` or provide this layer with data containing a single row.
```

``` r
pflag_comb
```

![](Analysis_achievement_files/figure-html/unnamed-chunk-24-1.png)<!-- -->

``` r
png(filename = 
      "../results/ach/GDI_flags.png",
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
mod2_log_GDP<-lmer(ach.z~gndr.c+log_gdp.z.cm+
                     gndr.c:log_gdp.z.cm+
                 (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2_log_GDP)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ach.z ~ gndr.c + log_gdp.z.cm + gndr.c:log_gdp.z.cm + (gndr.c |      cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1310059.7 1310147.7 -655021.9 1310043.7    441160 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -5.9287 -0.6682  0.0304  0.6625  5.0661 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr 
##  cntry    (Intercept) 0.08801  0.29667       
##           gndr.c      0.00383  0.06189  -0.29
##  Residual             1.00488  1.00244       
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##                     Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)          0.04156    0.05169 33.00246   0.804  0.42708    
## gndr.c               0.18308    0.01132 33.08845  16.175  < 2e-16 ***
## log_gdp.z.cm        -0.16836    0.05187 33.02487  -3.246  0.00269 ** 
## gndr.c:log_gdp.z.cm  0.01627    0.01149 34.45007   1.415  0.16591    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c lg_g..
## gndr.c      -0.271              
## lg_gdp.z.cm  0.022 -0.006       
## gndr.c:l_.. -0.006 -0.009 -0.268
```

``` r
getFE(mod2_log_GDP,round=3)
```

```
##                       Est.    SE     df      t     p     LL     UL
## (Intercept)          0.042 0.052 33.002  0.804 0.427 -0.064  0.147
## gndr.c               0.183 0.011 33.088 16.175 0.000  0.160  0.206
## log_gdp.z.cm        -0.168 0.052 33.025 -3.246 0.003 -0.274 -0.063
## gndr.c:log_gdp.z.cm  0.016 0.011 34.450  1.415 0.166 -0.007  0.040
```

``` r
getVC(mod2_log_GDP)
```

```
##        grp        var1   var2 sdcor  vcov
## 1    cntry (Intercept)   <NA>  0.30  0.09
## 2    cntry      gndr.c   <NA>  0.06  0.00
## 3    cntry (Intercept) gndr.c -0.29 -0.01
## 4 Residual        <NA>   <NA>  1.00  1.00
```

``` r
r2mlm(mod2_log_GDP,bargraph = F)
```

```
## $Decompositions
##                        total
## fixed           0.0276799965
## slope variation 0.0008456755
## mean variation  0.0785738425
## sigma2          0.8929004855
## 
## $R2s
##            total
## f   0.0276799965
## v   0.0008456755
## m   0.0785738425
## fv  0.0285256720
## fvm 0.1070995145
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
## Time difference of 28.51246 secs
```

``` r
round(ddsc_mod2_log_GDP$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.125        0.354        1.005     1.130 0.111   7204.152 0.999   0.999
## 2        0.5         0.109        0.330        1.005     1.114 0.098   6164.576 0.998   0.999
```

``` r
round(ddsc_mod2_log_GDP$descriptives,3)
```

```
##                          M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled log_gdp.z.cm
## means_y1             0.136 0.337    1.000           1.000    0.983           0.983       -0.482
## means_y1_scaled      0.390 0.968    1.000           1.000    0.983           0.983       -0.482
## means_y2            -0.045 0.359    0.983           0.983    1.000           1.000       -0.497
## means_y2_scaled     -0.130 1.031    0.983           0.983    1.000           1.000       -0.497
## log_gdp.z.cm        -0.022 1.012   -0.482          -0.482   -0.497          -0.497        1.000
## log_gdp.z.cm_scaled  0.000 1.000   -0.482          -0.482   -0.497          -0.497        1.000
## diff_score           0.181 0.068   -0.231          -0.231   -0.407          -0.407        0.237
## diff_score_scaled    0.521 0.196   -0.231          -0.231   -0.407          -0.407        0.237
##                     log_gdp.z.cm_scaled diff_score diff_score_scaled
## means_y1                         -0.482     -0.231            -0.231
## means_y1_scaled                  -0.482     -0.231            -0.231
## means_y2                         -0.497     -0.407            -0.407
## means_y2_scaled                  -0.497     -0.407            -0.407
## log_gdp.z.cm                      1.000      0.237             0.237
## log_gdp.z.cm_scaled               1.000      0.237             0.237
## diff_score                        0.237      1.000             1.000
## diff_score_scaled                 0.237      1.000             1.000
```

``` r
round(ddsc_mod2_log_GDP$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                      -0.238 0.168 34.450  -1.415   0.166   -0.580    0.104
## w_11                         -0.176 0.054 33.047  -3.287   0.002   -0.286   -0.067
## w_21                         -0.160 0.051 33.040  -3.164   0.003   -0.263   -0.057
## r_xy1                        -0.524 0.159 33.047  -3.287   0.002   -0.848   -0.200
## r_xy2                        -0.447 0.141 33.040  -3.164   0.003   -0.734   -0.159
## b_11                         -0.508 0.154 33.047  -3.287   0.002   -0.822   -0.193
## b_21                         -0.461 0.146 33.040  -3.164   0.003   -0.757   -0.165
## main_effect                  -0.168 0.052 33.025  -3.246   0.003   -0.274   -0.063
## moderator_effect              0.183 0.011 33.088  16.175   0.000    0.160    0.206
## interaction                   0.016 0.011 34.450   1.415   0.166   -0.007    0.040
## q_b11_b21                    -0.061    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                  -0.101    NA     NA      NA      NA       NA       NA
## cross_over_point            -11.255    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.152 0.050 33.070  -3.040   0.005   -0.254   -0.050
## interaction_vs_main_bscale   -0.437 0.144 33.070  -3.040   0.005   -0.730   -0.145
## interaction_vs_main_rscale   -0.408 0.135 33.073  -3.025   0.005   -0.682   -0.134
## dadas                        -0.320 0.101 33.040  -3.164   0.998   -0.526   -0.114
## dadas_bscale                 -0.921 0.291 33.040  -3.164   0.998   -1.514   -0.329
## dadas_rscale                 -0.893 0.282 33.040  -3.164   0.998   -1.468   -0.319
## abs_diff                      0.016 0.011 34.450   1.415   0.083   -0.007    0.040
## abs_sum                       0.337 0.104 33.025   3.246   0.001    0.126    0.548
## abs_diff_bscale               0.047 0.033 34.450   1.415   0.083   -0.020    0.114
## abs_sum_bscale                0.968 0.298 33.025   3.246   0.001    0.361    1.575
## abs_diff_rscale               0.077 0.037 34.758   2.107   0.021    0.003    0.152
## abs_sum_rscale                0.971 0.299 33.025   3.247   0.001    0.363    1.579
```

``` r
round(ddsc_mod2_log_GDP$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.008 -0.364  4.129  1.000  0.042
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
## r_xy1_y2                         -0.237 0.169  -1.399  0.162   -0.568    0.095
## r_xy1                            -0.497 0.151  -3.294  0.001   -0.793   -0.201
## r_xy2                            -0.482 0.153  -3.159  0.002   -0.781   -0.183
## b_11                             -0.513 0.156  -3.294  0.001   -0.818   -0.208
## b_21                             -0.466 0.148  -3.159  0.002   -0.756   -0.177
## b_10                             -0.130 0.153  -0.849  0.396   -0.431    0.170
## b_20                              0.390 0.145   2.686  0.007    0.106    0.675
## res_cov_y1_y2                     0.719 0.179   4.015  0.000    0.368    1.070
## diff_b10_b20                     -0.521 0.033 -15.914  0.000   -0.585   -0.456
## diff_b11_b21                     -0.046 0.033  -1.399  0.162   -0.112    0.019
## diff_rxy1_rxy2                   -0.016 0.032  -0.482  0.630   -0.079    0.048
## q_b11_b21                        -0.061 0.048  -1.272  0.203   -0.155    0.033
## q_rxy1_rxy2                      -0.020 0.042  -0.482  0.630   -0.104    0.063
## cross_over_point                -11.198 8.033  -1.394  0.163  -26.942    4.546
## sum_b11_b21                      -0.979 0.302  -3.247  0.001   -1.571   -0.388
## main_effect                      -0.490 0.151  -3.247  0.001   -0.785   -0.194
## interaction_vs_main_effect       -0.443 0.146  -3.029  0.002   -0.730   -0.156
## diff_abs_b11_abs_b21              0.046 0.033   1.399  0.162   -0.019    0.112
## abs_diff_b11_b21                  0.046 0.033   1.399  0.081   -0.019    0.112
## abs_sum_b11_b21                   0.979 0.302   3.247  0.001    0.388    1.571
## dadas                            -0.933 0.295  -3.159  0.999   -1.512   -0.354
## q_r_equivalence                  -0.080 0.042  -1.872  0.031       NA       NA
## q_b_equivalence                  -0.039 0.048  -0.807  0.210       NA       NA
## cross_over_point_equivalence     11.198 8.033   1.394  0.918       NA       NA
## cross_over_point_minimal_effect  11.198 8.033   1.394  0.082       NA       NA
```

``` r
round(ddsc_sem_log_GDP$variance_test,3)
```

```
##             est    se       z pvalue ci.lower ci.upper
## cov_y1y2  0.951 0.236   4.026  0.000    0.488    1.414
## var_y1    1.031 0.254   4.062  0.000    0.533    1.528
## var_y2    0.908 0.224   4.062  0.000    0.470    1.347
## var_diff  0.122 0.069   1.766  0.077   -0.013    0.258
## var_ratio 1.135 0.073  15.496  0.000    0.991    1.278
## cor_y1y2  0.983 0.006 164.296  0.000    0.971    0.994
```

### Figure


``` r
# plot the results

# start with obtaining predicted values for means and differences


mod2_log_GDP_unstd<-lmer(ach.z~gndr.c+log_gdp.cm+
                           gndr.c:log_gdp.cm+
                       (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
                     control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

mod2_log_GDP_unstd_red<-lmer(ach.z~gndr.c+
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


p1.ach.flags<-
  ggplot(p,aes(y=yvar,x=log_gdp.cm,color=gndr.c))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("log(GDP per capita) (Average 2002-2022)")+
  #ylim(c(min.y.pred,max.y.pred))+
  #xlim(c(0.60,1.00))+
  ylab("Value achievement (Average 2002-2022)")+
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

p2.ach.flags<-ggplot(p2,aes(y=yvar,x=log_gdp.cm))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("log(GDP per capita) (Average 2002-2022)")+
  #xlim(c(0.60,1.00))+
  #ylim(c(min.y.narrow,max.y.narrow))+
  ylab("Difference in value achievement")+
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
  ggarrange(p1.ach.flags,p2.ach.flags,align = "v",
            ncol=1,nrow=2,heights=c(2,1))
```

```
## Warning in geom_text(inherit.aes = F, aes(x = 9.75, y = -0.45, label = coef_q, : All aesthetics have length 1, but the data has 3302 rows.
## ℹ Please consider using `annotate()` or provide this layer with data containing a single row.
```

``` r
pflag_comb
```

![](Analysis_achievement_files/figure-html/unnamed-chunk-27-1.png)<!-- -->

``` r
png(filename = 
      "../results/ach/log_GDP_flags.png",
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
mod3<-lmer(ach.z~gndr.c+essround.c+(gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod3)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ach.z ~ gndr.c + essround.c + (gndr.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1310029.8 1310106.8 -655007.9 1310015.8    441161 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -5.9479 -0.6685  0.0297  0.6621  5.0732 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr 
##  cntry    (Intercept) 0.116203 0.34089       
##           gndr.c      0.004088 0.06393  -0.37
##  Residual             1.004798 1.00240       
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##               Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)  4.569e-02  5.937e-02  3.299e+01   0.770    0.447    
## gndr.c       1.829e-01  1.166e-02  3.306e+01  15.691  < 2e-16 ***
## essround.c  -3.521e-03  5.762e-04  4.411e+05  -6.111  9.9e-10 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##            (Intr) gndr.c
## gndr.c     -0.352       
## essround.c -0.001 -0.001
```

``` r
getFE(mod3,round=3)
```

```
##               Est.    SE         df      t     p     LL     UL
## (Intercept)  0.046 0.059     32.992  0.770 0.447 -0.075  0.166
## gndr.c       0.183 0.012     33.059 15.691 0.000  0.159  0.207
## essround.c  -0.004 0.001 441100.530 -6.111 0.000 -0.005 -0.002
```

``` r
getVC(mod3)
```

```
##        grp        var1   var2 sdcor  vcov
## 1    cntry (Intercept)   <NA>  0.34  0.12
## 2    cntry      gndr.c   <NA>  0.06  0.00
## 3    cntry (Intercept) gndr.c -0.37 -0.01
## 4 Residual        <NA>   <NA>  1.00  1.00
```

``` r
r2mlm(mod3,bargraph = F)
```

```
## $Decompositions
##                        total
## fixed           0.0074255781
## slope variation 0.0008980208
## mean variation  0.1032986666
## sigma2          0.8883777345
## 
## $R2s
##            total
## f   0.0074255781
## v   0.0008980208
## m   0.1032986666
## fv  0.0083235989
## fvm 0.1116222655
```

``` r
anova(mod2,mod3)
```

```
## Data: diff_dat
## Models:
## mod2: ach.z ~ gndr.c + (gndr.c | cntry)
## mod3: ach.z ~ gndr.c + essround.c + (gndr.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod2    6 1310065 1310131 -655027   1310053                         
## mod3    7 1310030 1310107 -655008   1310016 37.343  1  9.908e-10 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## mod4: random effect of time (Ess round)


``` r
mod4<-lmer(ach.z~gndr.c+essround.c+(gndr.c+essround.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod4)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ach.z ~ gndr.c + essround.c + (gndr.c + essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1308828.1 1308938.1 -654404.1 1308808.1    441158 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.1211 -0.6690  0.0298  0.6608  5.1476 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr       
##  cntry    (Intercept) 0.115483 0.33983             
##           gndr.c      0.004064 0.06375  -0.36      
##           essround.c  0.001207 0.03474   0.02  0.09
##  Residual             1.001725 1.00086             
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##              Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)  0.046144   0.059279 33.060284   0.778    0.442    
## gndr.c       0.182804   0.011625 33.139316  15.725   <2e-16 ***
## essround.c  -0.003727   0.006155 28.072943  -0.606    0.550    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##            (Intr) gndr.c
## gndr.c     -0.340       
## essround.c  0.016  0.087
```

``` r
getFE(mod4,round=3)
```

```
##               Est.    SE     df      t     p     LL    UL
## (Intercept)  0.046 0.059 33.060  0.778 0.442 -0.074 0.167
## gndr.c       0.183 0.012 33.139 15.725 0.000  0.159 0.206
## essround.c  -0.004 0.006 28.073 -0.606 0.550 -0.016 0.009
```

``` r
getVC(mod4)
```

```
##        grp        var1       var2 sdcor  vcov
## 1    cntry (Intercept)       <NA>  0.34  0.12
## 2    cntry      gndr.c       <NA>  0.06  0.00
## 3    cntry  essround.c       <NA>  0.03  0.00
## 4    cntry (Intercept)     gndr.c -0.36 -0.01
## 5    cntry (Intercept) essround.c  0.02  0.00
## 6    cntry      gndr.c essround.c  0.09  0.00
## 7 Residual        <NA>       <NA>  1.00  1.00
```

``` r
r2mlm(mod4,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.007391540
## slope variation 0.008914429
## mean variation  0.102152242
## sigma2          0.881541789
## 
## $R2s
##           total
## f   0.007391540
## v   0.008914429
## m   0.102152242
## fv  0.016305969
## fvm 0.118458211
```

``` r
anova(mod2,mod3,mod4)
```

```
## Data: diff_dat
## Models:
## mod2: ach.z ~ gndr.c + (gndr.c | cntry)
## mod3: ach.z ~ gndr.c + essround.c + (gndr.c | cntry)
## mod4: ach.z ~ gndr.c + essround.c + (gndr.c + essround.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L)    Chisq Df Pr(>Chisq)    
## mod2    6 1310065 1310131 -655027   1310053                           
## mod3    7 1310030 1310107 -655008   1310016   37.343  1  9.908e-10 ***
## mod4   10 1308828 1308938 -654404   1308808 1207.690  3  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## mod5: fixed interaction between time and gender


``` r
mod5<-lmer(ach.z~gndr.c+essround.c+
             gndr.c:essround.c+(gndr.c+essround.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod5)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ach.z ~ gndr.c + essround.c + gndr.c:essround.c + (gndr.c + essround.c |      cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1308767.0 1308888.0 -654372.5 1308745.0    441157 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.1468 -0.6692  0.0306  0.6616  5.1280 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr       
##  cntry    (Intercept) 0.115434 0.33976             
##           gndr.c      0.003726 0.06104  -0.38      
##           essround.c  0.001218 0.03490   0.02  0.09
##  Residual             1.001589 1.00079             
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##                     Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)        4.593e-02  5.927e-02  3.306e+01   0.775    0.444    
## gndr.c             1.837e-01  1.117e-02  3.366e+01  16.444  < 2e-16 ***
## essround.c        -3.831e-03  6.183e-03  2.800e+01  -0.620    0.541    
## gndr.c:essround.c -9.086e-03  1.142e-03  1.170e+05  -7.955 1.81e-15 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c essrn.
## gndr.c      -0.364              
## essround.c   0.014  0.083       
## gndr.c:ssr.  0.000 -0.009  0.002
```

``` r
getFE(mod5,round=3)
```

```
##                     Est.    SE         df      t     p     LL     UL
## (Intercept)        0.046 0.059     33.061  0.775 0.444 -0.075  0.167
## gndr.c             0.184 0.011     33.660 16.444 0.000  0.161  0.206
## essround.c        -0.004 0.006     28.001 -0.620 0.541 -0.016  0.009
## gndr.c:essround.c -0.009 0.001 116984.514 -7.955 0.000 -0.011 -0.007
```

``` r
getVC(mod5)
```

```
##        grp        var1       var2 sdcor  vcov
## 1    cntry (Intercept)       <NA>  0.34  0.12
## 2    cntry      gndr.c       <NA>  0.06  0.00
## 3    cntry  essround.c       <NA>  0.03  0.00
## 4    cntry (Intercept)     gndr.c -0.38 -0.01
## 5    cntry (Intercept) essround.c  0.02  0.00
## 6    cntry      gndr.c essround.c  0.09  0.00
## 7 Residual        <NA>       <NA>  1.00  1.00
```

``` r
r2mlm(mod5,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.007607711
## slope variation 0.008915735
## mean variation  0.102117223
## sigma2          0.881359332
## 
## $R2s
##           total
## f   0.007607711
## v   0.008915735
## m   0.102117223
## fv  0.016523446
## fvm 0.118640668
```

``` r
anova(mod4,mod5)
```

```
## Data: diff_dat
## Models:
## mod4: ach.z ~ gndr.c + essround.c + (gndr.c + essround.c | cntry)
## mod5: ach.z ~ gndr.c + essround.c + gndr.c:essround.c + (gndr.c + essround.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod4   10 1308828 1308938 -654404   1308808                         
## mod5   11 1308767 1308888 -654372   1308745 63.113  1  1.952e-15 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## mod6: random interaction between time and gender


``` r
mod6<-lmer(ach.z~gndr.c+essround.c+
             gndr.c:essround.c+(gndr.c+essround.c+gndr.c:essround.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ach.z ~ gndr.c + essround.c + gndr.c:essround.c + (gndr.c + essround.c +  
##     gndr.c:essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1308766.9 1308931.9 -654368.5 1308736.9    441153 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.1234 -0.6697  0.0311  0.6615  5.1201 
## 
## Random effects:
##  Groups   Name              Variance  Std.Dev. Corr             
##  cntry    (Intercept)       1.155e-01 0.339819                  
##           gndr.c            3.746e-03 0.061206 -0.41            
##           essround.c        1.232e-03 0.035099  0.02  0.11      
##           gndr.c:essround.c 5.186e-05 0.007202  0.11 -0.09 -0.27
##  Residual                   1.002e+00 1.000759                  
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##                    Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)        0.045379   0.059280 33.058865   0.766    0.449    
## gndr.c             0.184218   0.011262 34.067039  16.358  < 2e-16 ***
## essround.c        -0.003674   0.006218 27.612440  -0.591    0.559    
## gndr.c:essround.c -0.009440   0.001782 17.503311  -5.298 5.36e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c essrn.
## gndr.c      -0.389              
## essround.c   0.012  0.104       
## gndr.c:ssr.  0.078 -0.073 -0.187
```

``` r
getFE(mod6,round=3)
```

```
##                     Est.    SE     df      t     p     LL     UL
## (Intercept)        0.045 0.059 33.059  0.766 0.449 -0.075  0.166
## gndr.c             0.184 0.011 34.067 16.358 0.000  0.161  0.207
## essround.c        -0.004 0.006 27.612 -0.591 0.559 -0.016  0.009
## gndr.c:essround.c -0.009 0.002 17.503 -5.298 0.000 -0.013 -0.006
```

``` r
getVC(mod6)
```

```
##         grp              var1              var2 sdcor  vcov
## 1     cntry       (Intercept)              <NA>  0.34  0.12
## 2     cntry            gndr.c              <NA>  0.06  0.00
## 3     cntry        essround.c              <NA>  0.04  0.00
## 4     cntry gndr.c:essround.c              <NA>  0.01  0.00
## 5     cntry       (Intercept)            gndr.c -0.41 -0.01
## 6     cntry       (Intercept)        essround.c  0.02  0.00
## 7     cntry       (Intercept) gndr.c:essround.c  0.11  0.00
## 8     cntry            gndr.c        essround.c  0.11  0.00
## 9     cntry            gndr.c gndr.c:essround.c -0.09  0.00
## 10    cntry        essround.c gndr.c:essround.c -0.27  0.00
## 11 Residual              <NA>              <NA>  1.00  1.00
```

``` r
r2mlm(mod6,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.007648028
## slope variation 0.009126964
## mean variation  0.102171812
## sigma2          0.881053196
## 
## $R2s
##           total
## f   0.007648028
## v   0.009126964
## m   0.102171812
## fv  0.016774993
## fvm 0.118946804
```

``` r
anova(mod4,mod5,mod6)
```

```
## Data: diff_dat
## Models:
## mod4: ach.z ~ gndr.c + essround.c + (gndr.c + essround.c | cntry)
## mod5: ach.z ~ gndr.c + essround.c + gndr.c:essround.c + (gndr.c + essround.c | cntry)
## mod6: ach.z ~ gndr.c + essround.c + gndr.c:essround.c + (gndr.c + essround.c + gndr.c:essround.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L)   Chisq Df Pr(>Chisq)    
## mod4   10 1308828 1308938 -654404   1308808                          
## mod5   11 1308767 1308888 -654372   1308745 63.1133  1  1.952e-15 ***
## mod6   15 1308767 1308932 -654368   1308737  8.0628  4    0.08931 .  
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
##         4.5 -0.0420 0.0679 Inf   -0.1752    0.0911  -0.619  0.5362
##        -4.5 -0.0514 0.0684 Inf   -0.1855    0.0826  -0.752  0.4520
## 
## gndr.c =  0.5:
##  essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##         4.5  0.0997 0.0644 Inf   -0.0266    0.2260   1.548  0.1217
##        -4.5  0.1753 0.0627 Inf    0.0523    0.2982   2.794  0.0052
## 
## Degrees-of-freedom method: asymptotic 
## Confidence level used: 0.95
```

``` r
pairs(change_mod6,adjust="none",infer=c(T,T))
```

```
## gndr.c = -0.5:
##  contrast                         estimate    SE  df asymp.LCL asymp.UCL z.ratio p.value
##  essround.c4.5 - (essround.c-4.5)  0.00942 0.058 Inf    -0.104    0.1231   0.162  0.8710
## 
## gndr.c =  0.5:
##  contrast                         estimate    SE  df asymp.LCL asymp.UCL z.ratio p.value
##  essround.c4.5 - (essround.c-4.5) -0.07554 0.055 Inf    -0.183    0.0323  -1.373  0.1698
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
##    -0.5        4.5 -0.0420 0.0679 Inf   -0.1752    0.0911  -0.619  0.5362
##     0.5        4.5  0.0997 0.0644 Inf   -0.0266    0.2260   1.548  0.1217
##    -0.5       -4.5 -0.0514 0.0684 Inf   -0.1855    0.0826  -0.752  0.4520
##     0.5       -4.5  0.1753 0.0627 Inf    0.0523    0.2982   2.794  0.0052
## 
## Degrees-of-freedom method: asymptotic 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod6,adjust="none",infer=c(T,T))
```

```
##  contrast                                                 estimate     SE  df asymp.LCL asymp.UCL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5     -0.14174 0.0133 Inf    -0.168   -0.1156
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  0.00942 0.0580 Inf    -0.104    0.1231
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)  -0.21728 0.0559 Inf    -0.327   -0.1077
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)     0.15116 0.0582 Inf     0.037    0.2653
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)     -0.07554 0.0550 Inf    -0.183    0.0323
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5) -0.22670 0.0143 Inf    -0.255   -0.1987
##  z.ratio p.value
##  -10.628  <.0001
##    0.162  0.8710
##   -3.886  0.0001
##    2.596  0.0094
##   -1.373  0.1698
##  -15.858  <.0001
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
##  diff_ESS10    0.142 0.0133 Inf     0.116     0.168  10.628  <.0001
##  diff_ESS1     0.227 0.0143 Inf     0.199     0.255  15.858  <.0001
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
##  diff_ESS10 - diff_ESS1   -0.085 0.016 Inf    -0.116   -0.0535  -5.298  <.0001
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
      obs_mean_wt=weighted.mean(x=ach.z,w=pspwght),
      obs_sd_wt=sqrt(wtd.var(ach.z,pspwght)),
      obs_mean=mean(ach.z),
      obs_sd=sd(ach.z),
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
  ylab("Mean-level of value achievement")+
  theme(legend.title=element_blank(),
        legend.text  = element_text(size = 16),
        axis.title.x = element_text(size = 18, face = "bold"),   
        axis.title.y = element_text(size = 18, face = "bold"),   
        axis.text.x  = element_text(size = 14),                  
        axis.text.y  = element_text(size = 14))

p_time_trends
```

![](Analysis_achievement_files/figure-html/unnamed-chunk-33-1.png)<!-- -->

``` r
png(filename = 
      "../results/ach/time_trends.png",
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
pred_cntry_dat$ach.z_mean<-predict(mod6,newdata=pred_cntry_dat)

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

range(pred_cntry_dat$ach.z_mean)
```

```
## [1] -0.7385374  0.8241162
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
      obs_mean_wt=weighted.mean(x=ach.z,w=pspwght),
      obs_sd_wt=sqrt(wtd.var(ach.z,pspwght)),
      obs_mean=mean(ach.z),
      obs_sd=sd(ach.z),
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

pdf("../results/ach/country_specific_time_trends.pdf", width = 8, height = 6)

for (ctry in countries) {
  print(
    ggplot(pred_cntry_dat[pred_cntry_dat$cntry == ctry, ], 
           aes(x = year, y = ach.z_mean, color = gender)) +
      geom_smooth(method = "lm",formula = "y ~ x", se = FALSE) +
      geom_point(size=8) +
      geom_point(aes(x=year,y=obs_mean_wt),size=8,shape = 1,alpha=.50)+
      geom_errorbar(aes(ymin=obs_mean_wt_LL, ymax=obs_mean_wt_UL),alpha = .50)+
      scale_color_manual(values = my_colors) +
      geom_flag(aes(country=tolower(ctry)))+
      scale_y_continuous(limits = c(-1.1, 1.1)) +
      scale_x_continuous(limits = c(2001, 2021), breaks = seq(2002, 2020, 2)) +
      ggtitle(paste("Country:", ctry))+
    ylab("Mean-level of value achievement")+
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
         aes(x = year, y = ach.z_mean, color = gender)) +
  geom_smooth(method = "lm", se = FALSE, formula="y~x") +
  geom_point(size=3.5) +
  geom_point(aes(x=year,y=obs_mean_wt),size=3.5,shape = 1,alpha=.70)+
  geom_errorbar(aes(ymin=obs_mean_wt_LL, ymax=obs_mean_wt_UL),alpha = .50)+
  scale_color_manual(values = my_colors) +
  #geom_flag(aes(country=tolower(ctry)))+
  scale_y_continuous(limits = c(-1.1, 1.1)) +
  scale_x_continuous(limits = c(2001, 2021), breaks = seq(2002, 2020, 2)) +
  #ggtitle(paste("Country:", ctry))+
  ylab("Mean-level of value achievement")+
  xlab("Year")+
  theme(legend.title=element_blank(),legend.position = "top",
        axis.text.x = element_text(angle = 45,size = 6,hjust=1))+
  facet_wrap(~CLDR,nrow=6,ncol=6)+
  #facet_wrap(~cntry,nrow=6,ncol=6)+
  geom_flag(aes(country=tolower(cntry)),size=2)

facet_plot
```

![](Analysis_achievement_files/figure-html/unnamed-chunk-35-1.png)<!-- -->

``` r
png(filename = 
      "../results/ach/country_time_trend_facets.png",
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
## [1] 23.7536
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
## 1    0.28                0.05                   -0.09                    0.00                      0.09
## 2    0.14                0.16                   -0.07                    0.12                      0.19
## 3    0.17               -0.22                    0.02                   -0.21                     -0.23
## 4    0.28                0.16                   -0.07                    0.13                      0.20
## 5    0.12               -0.24                   -0.10                   -0.29                     -0.19
## 6    0.30                0.29                   -0.12                    0.23                      0.35
## 7    0.27               -0.08                   -0.10                   -0.13                     -0.04
## 8    0.18                0.13                   -0.02                    0.12                      0.14
## 9    0.18               -0.13                   -0.09                   -0.18                     -0.08
## 10   0.21               -0.24                   -0.09                   -0.28                     -0.19
## 11   0.31               -0.16                   -0.04                   -0.18                     -0.14
## 12   0.23               -0.16                   -0.11                   -0.22                     -0.11
## 13   0.24               -0.11                   -0.06                   -0.14                     -0.07
## 14   0.15               -0.30                   -0.10                   -0.35                     -0.24
## 15   0.17               -0.21                   -0.11                   -0.27                     -0.16
## 16   0.16               -0.02                   -0.14                   -0.09                      0.05
## 17   0.15               -0.08                   -0.04                   -0.10                     -0.06
## 18   0.05               -0.12                   -0.02                   -0.13                     -0.12
## 19   0.13               -0.14                   -0.16                   -0.22                     -0.06
## 20   0.17                0.19                   -0.12                    0.13                      0.25
## 21   0.11               -0.31                    0.00                   -0.31                     -0.32
## 22   0.17               -1.16                   -0.03                   -1.18                     -1.15
## 23   0.17                0.57                   -0.14                    0.50                      0.65
## 24   0.24                0.03                   -0.09                   -0.01                      0.08
## 25   0.16               -0.01                   -0.09                   -0.06                      0.04
## 26   0.21               -0.01                   -0.07                   -0.05                      0.02
## 27   0.17               -0.07                   -0.15                   -0.15                      0.00
## 28   0.14               -0.29                   -0.14                   -0.36                     -0.22
## 29   0.15               -0.06                   -0.02                   -0.07                     -0.05
## 30   0.13                0.31                   -0.07                    0.28                      0.34
## 31   0.23                0.09                   -0.11                    0.03                      0.14
## 32   0.15                0.37                   -0.11                    0.32                      0.42
## 33   0.16                0.70                   -0.15                    0.63                      0.78
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
## 1     LT               -0.31
## 2     GR               -0.30
## 3     RU               -0.29
## 4     CY               -0.24
## 5     ES               -0.24
## 6     BG               -0.22
## 7     FI               -0.16
## 8     FR               -0.16
## 9     IS               -0.14
## 10    EE               -0.13
## 11    IL               -0.12
## 12    GB               -0.11
## 13    DE               -0.08
## 14    IE               -0.08
## 15    PT               -0.07
## 16    SE               -0.06
## 17    HU               -0.02
## 18    NO               -0.01
## 19    PL               -0.01
## 20    NL                0.03
## 21    AT                0.05
## 22    SK                0.09
## 23    DK                0.13
## 24    BE                0.16
## 25    CH                0.16
## 26    CZ                0.29
## 27    SI                0.31
## 28    UA                0.70
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
## 1     IS                   -0.16
## 2     PT                   -0.15
## 3     UA                   -0.15
## 4     HU                   -0.14
## 5     RU                   -0.14
## 6     CZ                   -0.12
## 7     FR                   -0.11
## 8     SK                   -0.11
## 9     CY                   -0.10
## 10    DE                   -0.10
## 11    GR                   -0.10
## 12    AT                   -0.09
## 13    EE                   -0.09
## 14    ES                   -0.09
## 15    NL                   -0.09
## 16    NO                   -0.09
## 17    BE                   -0.07
## 18    CH                   -0.07
## 19    PL                   -0.07
## 20    SI                   -0.07
## 21    GB                   -0.06
## 22    FI                   -0.04
## 23    IE                   -0.04
## 24    DK                   -0.02
## 25    IL                   -0.02
## 26    SE                   -0.02
## 27    LT                    0.00
## 28    BG                    0.02
```


#### Country-level predictors: GEI


``` r
mod6_GEI<-lmer(ach.z~gndr.c+essround.c+
             gndr.c:essround.c+
               gei.z.cm:gndr.c+gei.z.cm:essround.c+gei.z.cm:gndr.c:essround.c+
               (gndr.c+essround.c+gndr.c:essround.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_GEI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ach.z ~ gndr.c + essround.c + gndr.c:essround.c + gei.z.cm:gndr.c +  
##     gei.z.cm:essround.c + gei.z.cm:gndr.c:essround.c + (gndr.c +  
##     essround.c + gndr.c:essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1277293.7 1277491.3 -638628.9 1277257.7    431760 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.1418 -0.6704  0.0326  0.6637  5.1274 
## 
## Random effects:
##  Groups   Name              Variance  Std.Dev. Corr             
##  cntry    (Intercept)       1.183e-01 0.344012                  
##           gndr.c            3.807e-03 0.061704 -0.42            
##           essround.c        9.606e-04 0.030993  0.21  0.08      
##           gndr.c:essround.c 5.009e-05 0.007078  0.48 -0.27 -0.26
##  Residual                   9.978e-01 0.998895                  
## Number of obs: 431778, groups:  cntry, 32
## 
## Fixed effects:
##                             Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)                 0.048352   0.060923 32.044849   0.794  0.43323    
## gndr.c                      0.184195   0.011514 28.345851  15.997 1.01e-15 ***
## essround.c                 -0.006851   0.005592 24.651137  -1.225  0.23210    
## gndr.c:essround.c          -0.009780   0.001807 17.251772  -5.413 4.43e-05 ***
## gndr.c:gei.z.cm             0.003266   0.011106 36.528810   0.294  0.77038    
## essround.c:gei.z.cm         0.008919   0.005683 26.665997   1.569  0.12837    
## gndr.c:essround.c:gei.z.cm  0.005912   0.002040 22.742776   2.898  0.00817 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c essrn. gnd.:. gn.:.. es.:..
## gndr.c      -0.399                                   
## essround.c   0.201  0.083                            
## gndr.c:ssr.  0.335 -0.190 -0.185                     
## gndr.c:g.z.  0.000 -0.038 -0.002 -0.050              
## essrnd.c:.. -0.002 -0.001 -0.031  0.003  0.174       
## gndr.c:.:..  0.000 -0.024  0.004 -0.223  0.088 -0.232
```

``` r
getFE(mod6_GEI,round=2,p.round = 5)
```

```
##                             Est.   SE    df     t       p    LL    UL
## (Intercept)                 0.05 0.06 32.04  0.79 0.43323 -0.08  0.17
## gndr.c                      0.18 0.01 28.35 16.00 0.00000  0.16  0.21
## essround.c                 -0.01 0.01 24.65 -1.23 0.23210 -0.02  0.00
## gndr.c:essround.c          -0.01 0.00 17.25 -5.41 0.00004 -0.01 -0.01
## gndr.c:gei.z.cm             0.00 0.01 36.53  0.29 0.77038 -0.02  0.03
## essround.c:gei.z.cm         0.01 0.01 26.67  1.57 0.12837  0.00  0.02
## gndr.c:essround.c:gei.z.cm  0.01 0.00 22.74  2.90 0.00817  0.00  0.01
```

``` r
getVC(mod6_GEI)
```

```
##         grp              var1              var2 sdcor  vcov
## 1     cntry       (Intercept)              <NA>  0.34  0.12
## 2     cntry            gndr.c              <NA>  0.06  0.00
## 3     cntry        essround.c              <NA>  0.03  0.00
## 4     cntry gndr.c:essround.c              <NA>  0.01  0.00
## 5     cntry       (Intercept)            gndr.c -0.42 -0.01
## 6     cntry       (Intercept)        essround.c  0.21  0.00
## 7     cntry       (Intercept) gndr.c:essround.c  0.48  0.00
## 8     cntry            gndr.c        essround.c  0.08  0.00
## 9     cntry            gndr.c gndr.c:essround.c -0.27  0.00
## 10    cntry        essround.c gndr.c:essround.c -0.26  0.00
## 11 Residual              <NA>              <NA>  1.00  1.00
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_GEI,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] 22.02739
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_GEI,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] 3.415697
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
##         4.5 -0.0226 0.0752 Inf   -0.1701     0.125  -0.300  0.7638
##        -4.5  0.1193 0.0665 Inf   -0.0109     0.250   1.795  0.0726
## 
## gei.z.cm =  0:
##  essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##         4.5  0.0175 0.0704 Inf   -0.1205     0.156   0.249  0.8035
##        -4.5  0.0792 0.0611 Inf   -0.0405     0.199   1.297  0.1948
## 
## gei.z.cm =  1:
##  essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##         4.5  0.0577 0.0746 Inf   -0.0886     0.204   0.773  0.4397
##        -4.5  0.0390 0.0660 Inf   -0.0902     0.168   0.592  0.5538
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
##  essround.c4.5 - (essround.c-4.5)  -0.1419 0.0729 Inf    -0.285  0.000899  -1.948  0.0515
## 
## gei.z.cm =  0:
##  contrast                         estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  essround.c4.5 - (essround.c-4.5)  -0.0617 0.0503 Inf    -0.160  0.036980  -1.225  0.2205
## 
## gei.z.cm =  1:
##  contrast                         estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   0.0186 0.0706 Inf    -0.120  0.157030   0.264  0.7921
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
##    -0.5        4.5 -0.07777 0.0769 Inf   -0.2284    0.0729  -1.012  0.3116
##     0.5        4.5  0.03255 0.0751 Inf   -0.1146    0.1797   0.433  0.6647
##    -0.5       -4.5 -0.00646 0.0717 Inf   -0.1470    0.1341  -0.090  0.9283
##     0.5       -4.5  0.24509 0.0626 Inf    0.1224    0.3677   3.917  0.0001
## 
## gei.z.cm =  0:
##  gndr.c essround.c   emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##    -0.5        4.5 -0.05257 0.0716 Inf   -0.1929    0.0878  -0.734  0.4630
##     0.5        4.5  0.08762 0.0698 Inf   -0.0492    0.2244   1.255  0.2094
##    -0.5       -4.5 -0.03492 0.0655 Inf   -0.1634    0.0935  -0.533  0.5941
##     0.5       -4.5  0.19328 0.0573 Inf    0.0810    0.3055   3.375  0.0007
## 
## gei.z.cm =  1:
##  gndr.c essround.c   emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##    -0.5        4.5 -0.02737 0.0760 Inf   -0.1764    0.1217  -0.360  0.7189
##     0.5        4.5  0.14268 0.0743 Inf   -0.0029    0.2883   1.921  0.0547
##    -0.5       -4.5 -0.06339 0.0711 Inf   -0.2028    0.0760  -0.891  0.3728
##     0.5       -4.5  0.14148 0.0620 Inf    0.0200    0.2629   2.283  0.0224
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
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.1103 0.0211 Inf   -0.1517   -0.0689
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  -0.0713 0.0766 Inf   -0.2215    0.0788
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.3229 0.0726 Inf   -0.4652   -0.1805
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.0390 0.0766 Inf   -0.1112    0.1892
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      -0.2125 0.0715 Inf   -0.3528   -0.0723
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.2515 0.0213 Inf   -0.2933   -0.2098
##  z.ratio p.value
##   -5.226  <.0001
##   -0.931  0.3519
##   -4.444  <.0001
##    0.509  0.6108
##   -2.971  0.0030
##  -11.821  <.0001
## 
## gei.z.cm =  0:
##  contrast                                                 estimate     SE  df asymp.LCL asymp.UCL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.1402 0.0128 Inf   -0.1652   -0.1152
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  -0.0176 0.0524 Inf   -0.1204    0.0851
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.2459 0.0507 Inf   -0.3452   -0.1465
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.1225 0.0526 Inf    0.0195    0.2255
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      -0.1057 0.0495 Inf   -0.2026   -0.0087
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.2282 0.0153 Inf   -0.2582   -0.1982
##  z.ratio p.value
##  -10.979  <.0001
##   -0.336  0.7365
##   -4.851  <.0001
##    2.332  0.0197
##   -2.136  0.0327
##  -14.908  <.0001
## 
## gei.z.cm =  1:
##  contrast                                                 estimate     SE  df asymp.LCL asymp.UCL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.1701 0.0182 Inf   -0.2058   -0.1344
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   0.0360 0.0739 Inf   -0.1089    0.1809
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.1688 0.0703 Inf   -0.3067   -0.0310
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.2061 0.0743 Inf    0.0604    0.3517
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)       0.0012 0.0689 Inf   -0.1338    0.1362
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.2049 0.0199 Inf   -0.2438   -0.1659
##  z.ratio p.value
##   -9.336  <.0001
##    0.487  0.6262
##   -2.401  0.0163
##    2.773  0.0056
##    0.017  0.9861
##  -10.307  <.0001
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
##  diff_ESS10    0.110 0.0211 Inf    0.0689     0.152   5.226  <.0001
##  diff_ESS1     0.252 0.0213 Inf    0.2098     0.293  11.821  <.0001
## 
## gei.z.cm =  0:
##  contrast   estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10    0.140 0.0128 Inf    0.1152     0.165  10.979  <.0001
##  diff_ESS1     0.228 0.0153 Inf    0.1982     0.258  14.908  <.0001
## 
## gei.z.cm =  1:
##  contrast   estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10    0.170 0.0182 Inf    0.1344     0.206   9.336  <.0001
##  diff_ESS1     0.205 0.0199 Inf    0.1659     0.244  10.307  <.0001
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
##  diff_ESS10 - diff_ESS1  -0.1412 0.0271 Inf   -0.1944   -0.0881  -5.210  <.0001
## 
## gei.z.cm =  0:
##  contrast               estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10 - diff_ESS1  -0.0880 0.0163 Inf   -0.1199   -0.0562  -5.413  <.0001
## 
## gei.z.cm =  1:
##  contrast               estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10 - diff_ESS1  -0.0348 0.0216 Inf   -0.0772    0.0076  -1.609  0.1077
## 
## Degrees-of-freedom method: asymptotic 
## Confidence level used: 0.95
```


#### Country-level predictors: GGGI


``` r
mod6_GGGI<-lmer(ach.z~gndr.c+essround.c+
             gndr.c:essround.c+
               gggi.z.cm:gndr.c+gggi.z.cm:essround.c+gggi.z.cm:gndr.c:essround.c+
               (gndr.c+essround.c+gndr.c:essround.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_GGGI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ach.z ~ gndr.c + essround.c + gndr.c:essround.c + gggi.z.cm:gndr.c +  
##     gggi.z.cm:essround.c + gggi.z.cm:gndr.c:essround.c + (gndr.c +  
##     essround.c + gndr.c:essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##  932069.9  932261.8 -466017.0  932033.9    314628 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.0918 -0.6652  0.0276  0.6615  5.1294 
## 
## Random effects:
##  Groups   Name              Variance  Std.Dev. Corr             
##  cntry    (Intercept)       0.1224091 0.34987                   
##           gndr.c            0.0043457 0.06592  -0.56            
##           essround.c        0.0007192 0.02682  -0.11  0.14      
##           gndr.c:essround.c 0.0001655 0.01286   0.42 -0.34 -0.27
##  Residual                   0.9954460 0.99772                   
## Number of obs: 314646, groups:  cntry, 33
## 
## Fixed effects:
##                              Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)                  0.081476   0.061087 32.912091   1.334   0.1914    
## gndr.c                       0.182261   0.012328 23.692244  14.784 1.86e-13 ***
## essround.c                  -0.008158   0.005027 27.066890  -1.623   0.1162    
## gndr.c:essround.c           -0.008082   0.002977 21.081045  -2.715   0.0129 *  
## gndr.c:gggi.z.cm            -0.010828   0.011091 35.791182  -0.976   0.3355    
## essround.c:gggi.z.cm        -0.005531   0.005451 29.005992  -1.015   0.3186    
## gndr.c:essround.c:gggi.z.cm  0.006409   0.003159 24.479355   2.029   0.0535 .  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c essrn. gnd.:. gn.:.. es.:..
## gndr.c      -0.526                                   
## essround.c  -0.107  0.131                            
## gndr.c:ssr.  0.316 -0.280 -0.204                     
## gndr.c:gg..  0.000 -0.018 -0.004 -0.026              
## essrnd.c:.. -0.004 -0.002 -0.089  0.010  0.081       
## gndr.c:.:..  0.001 -0.026  0.011 -0.080 -0.137 -0.172
```

``` r
getFE(mod6_GGGI,round=2,p.round = 5)
```

```
##                              Est.   SE    df     t       p    LL   UL
## (Intercept)                  0.08 0.06 32.91  1.33 0.19144 -0.04 0.21
## gndr.c                       0.18 0.01 23.69 14.78 0.00000  0.16 0.21
## essround.c                  -0.01 0.01 27.07 -1.62 0.11618 -0.02 0.00
## gndr.c:essround.c           -0.01 0.00 21.08 -2.72 0.01294 -0.01 0.00
## gndr.c:gggi.z.cm            -0.01 0.01 35.79 -0.98 0.33549 -0.03 0.01
## essround.c:gggi.z.cm        -0.01 0.01 29.01 -1.01 0.31861 -0.02 0.01
## gndr.c:essround.c:gggi.z.cm  0.01 0.00 24.48  2.03 0.05347  0.00 0.01
```

``` r
getVC(mod6_GGGI)
```

```
##         grp              var1              var2 sdcor  vcov
## 1     cntry       (Intercept)              <NA>  0.35  0.12
## 2     cntry            gndr.c              <NA>  0.07  0.00
## 3     cntry        essround.c              <NA>  0.03  0.00
## 4     cntry gndr.c:essround.c              <NA>  0.01  0.00
## 5     cntry       (Intercept)            gndr.c -0.56 -0.01
## 6     cntry       (Intercept)        essround.c -0.11  0.00
## 7     cntry       (Intercept) gndr.c:essround.c  0.42  0.00
## 8     cntry            gndr.c        essround.c  0.14  0.00
## 9     cntry            gndr.c gndr.c:essround.c -0.34  0.00
## 10    cntry        essround.c gndr.c:essround.c -0.27  0.00
## 11 Residual              <NA>              <NA>  1.00  1.00
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_GGGI,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] 41.61901
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_GGGI,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] -219.115
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
##         4.5 0.0697 0.0683 34.6 -0.06897    0.208   1.020  0.3146
##        -4.5 0.0933 0.0723 34.1 -0.05362    0.240   1.290  0.2056
## 
## gggi.z.cm =  0:
##  essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5 0.0448 0.0628 26.7 -0.08421    0.174   0.713  0.4823
##        -4.5 0.1182 0.0674 27.3 -0.01999    0.256   1.754  0.0906
## 
## gggi.z.cm =  1:
##  essround.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5 0.0199 0.0666 32.4 -0.11575    0.155   0.298  0.7674
##        -4.5 0.1431 0.0711 32.8 -0.00162    0.288   2.012  0.0525
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
##  essround.c4.5 - (essround.c-4.5)  -0.0236 0.0696 28.7   -0.166  0.11882  -0.340  0.7366
## 
## gggi.z.cm =  0:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)  -0.0734 0.0452 27.1   -0.166  0.01939  -1.623  0.1162
## 
## gggi.z.cm =  1:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)  -0.1232 0.0637 28.0   -0.254  0.00729  -1.934  0.0633
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
##    -0.5        4.5  0.00571 0.0708 34.8 -0.13798   0.1494   0.081  0.9361
##     0.5        4.5  0.13359 0.0678 34.8 -0.00407   0.2713   1.971  0.0568
##    -0.5       -4.5 -0.03585 0.0792 34.4 -0.19675   0.1250  -0.453  0.6536
##     0.5       -4.5  0.22245 0.0677 34.0  0.08479   0.3601   3.284  0.0024
## 
## gggi.z.cm =  0:
##  gndr.c essround.c   emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5        4.5 -0.02818 0.0646 26.3 -0.16082   0.1045  -0.436  0.6660
##     0.5        4.5  0.11771 0.0620 26.4 -0.00965   0.2451   1.899  0.0686
##    -0.5       -4.5  0.00887 0.0735 27.1 -0.14187   0.1596   0.121  0.9048
##     0.5       -4.5  0.22750 0.0624 26.5  0.09938   0.3556   3.646  0.0011
## 
## gggi.z.cm =  1:
##  gndr.c essround.c   emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5        4.5 -0.06208 0.0689 32.5 -0.20238   0.0782  -0.901  0.3743
##     0.5        4.5  0.10182 0.0660 32.4 -0.03261   0.2363   1.542  0.1328
##    -0.5       -4.5  0.05360 0.0779 33.0 -0.10498   0.2122   0.688  0.4965
##     0.5       -4.5  0.23256 0.0665 32.5  0.09719   0.3679   3.497  0.0014
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
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.1279 0.0239 27.5 -0.17695 -0.07881  -5.343
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   0.0416 0.0760 28.5 -0.11393  0.19706   0.547
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.2167 0.0699 28.4 -0.35983 -0.07363  -3.100
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.1694 0.0733 30.2  0.01986  0.31903   2.313
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      -0.0889 0.0689 27.5 -0.23010  0.05239  -1.290
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.2583 0.0285 32.2 -0.31629 -0.20031  -9.072
##  p.value
##   <.0001
##   0.5886
##   0.0043
##   0.0277
##   0.2079
##   <.0001
## 
## gggi.z.cm =  0:
##  contrast                                                 estimate     SE   df lower.CL upper.CL t.ratio
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.1459 0.0155 26.2 -0.17766 -0.11412  -9.435
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  -0.0371 0.0497 25.9 -0.13930  0.06519  -0.745
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.2557 0.0453 27.2 -0.34861 -0.16276  -5.643
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.1088 0.0484 26.4  0.00938  0.20829   2.248
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      -0.1098 0.0445 26.1 -0.20122 -0.01837  -2.468
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.2186 0.0206 22.4 -0.26128 -0.17598 -10.621
##  p.value
##   <.0001
##   0.4629
##   <.0001
##   0.0332
##   0.0205
##   <.0001
## 
## gggi.z.cm =  1:
##  contrast                                                 estimate     SE   df lower.CL upper.CL t.ratio
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.1639 0.0217 26.9 -0.20834 -0.11947  -7.570
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  -0.1157 0.0698 27.6 -0.25871  0.02736  -1.658
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.2946 0.0640 27.6 -0.42591 -0.16337  -4.600
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.0482 0.0675 29.3 -0.08977  0.18622   0.714
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      -0.1307 0.0628 26.4 -0.25982 -0.00165  -2.080
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.1790 0.0278 28.7 -0.23587 -0.12206  -6.434
##  p.value
##   <.0001
##   0.1087
##   0.0001
##   0.4806
##   0.0473
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
##  diff_ESS10    0.128 0.0239 27.5   0.0788    0.177   5.343  <.0001
##  diff_ESS1     0.258 0.0285 32.2   0.2003    0.316   9.072  <.0001
## 
## gggi.z.cm =  0:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10    0.146 0.0155 26.2   0.1141    0.178   9.435  <.0001
##  diff_ESS1     0.219 0.0206 22.4   0.1760    0.261  10.621  <.0001
## 
## gggi.z.cm =  1:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10    0.164 0.0217 26.9   0.1195    0.208   7.570  <.0001
##  diff_ESS1     0.179 0.0278 28.7   0.1221    0.236   6.434  <.0001
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
##  diff_ESS10 - diff_ESS1  -0.1304 0.0406 25.7  -0.2139  -0.0469  -3.214  0.0035
## 
## gggi.z.cm =  0:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10 - diff_ESS1  -0.0727 0.0268 21.1  -0.1284  -0.0170  -2.715  0.0129
## 
## gggi.z.cm =  1:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10 - diff_ESS1  -0.0151 0.0375 22.7  -0.0927   0.0625  -0.402  0.6916
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

#### Country-level predictors: GDI


``` r
mod6_GDI<-lmer(ach.z~gndr.c+essround.c+
             gndr.c:essround.c+
               gdi.z.cm:gndr.c+gdi.z.cm:essround.c+gdi.z.cm:gndr.c:essround.c+
               (gndr.c+essround.c+gndr.c:essround.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_GDI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ach.z ~ gndr.c + essround.c + gndr.c:essround.c + gdi.z.cm:gndr.c +  
##     gdi.z.cm:essround.c + gdi.z.cm:gndr.c:essround.c + (gndr.c +  
##     essround.c + gndr.c:essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1308765.1 1308963.0 -654364.5 1308729.1    441150 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.1202 -0.6697  0.0307  0.6613  5.1199 
## 
## Random effects:
##  Groups   Name              Variance  Std.Dev. Corr             
##  cntry    (Intercept)       1.162e-01 0.340927                  
##           gndr.c            3.627e-03 0.060225 -0.47            
##           essround.c        1.045e-03 0.032322 -0.05  0.05      
##           gndr.c:essround.c 5.498e-05 0.007415  0.16 -0.07 -0.24
##  Residual                   1.002e+00 1.000757                  
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##                             Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)                 0.045408   0.059468 33.037575   0.764    0.451    
## gndr.c                      0.184471   0.011102 33.330099  16.616  < 2e-16 ***
## essround.c                 -0.003426   0.005743 28.099678  -0.597    0.556    
## gndr.c:essround.c          -0.009449   0.001810 17.746074  -5.219 6.06e-05 ***
## gndr.c:gdi.z.cm            -0.016243   0.010613 38.096503  -1.531    0.134    
## essround.c:gdi.z.cm        -0.015238   0.005980 31.261272  -2.548    0.016 *  
## gndr.c:essround.c:gdi.z.cm  0.001514   0.002285 26.869643   0.663    0.513    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c essrn. gnd.:. gn.:.. es.:..
## gndr.c      -0.439                                   
## essround.c  -0.057  0.048                            
## gndr.c:ssr.  0.113 -0.062 -0.173                     
## gndr.c:gd..  0.000 -0.018  0.000 -0.017              
## essrnd.c:..  0.002 -0.002 -0.029  0.008  0.027       
## gndr.c:.:.. -0.001 -0.022  0.008  0.012  0.101 -0.133
```

``` r
getFE(mod6_GDI,round=2,p.round = 5)
```

```
##                             Est.   SE    df     t       p    LL    UL
## (Intercept)                 0.05 0.06 33.04  0.76 0.45055 -0.08  0.17
## gndr.c                      0.18 0.01 33.33 16.62 0.00000  0.16  0.21
## essround.c                  0.00 0.01 28.10 -0.60 0.55553 -0.02  0.01
## gndr.c:essround.c          -0.01 0.00 17.75 -5.22 0.00006 -0.01 -0.01
## gndr.c:gdi.z.cm            -0.02 0.01 38.10 -1.53 0.13415 -0.04  0.01
## essround.c:gdi.z.cm        -0.02 0.01 31.26 -2.55 0.01596 -0.03  0.00
## gndr.c:essround.c:gdi.z.cm  0.00 0.00 26.87  0.66 0.51323  0.00  0.01
```

``` r
getVC(mod6_GDI)
```

```
##         grp              var1              var2 sdcor  vcov
## 1     cntry       (Intercept)              <NA>  0.34  0.12
## 2     cntry            gndr.c              <NA>  0.06  0.00
## 3     cntry        essround.c              <NA>  0.03  0.00
## 4     cntry gndr.c:essround.c              <NA>  0.01  0.00
## 5     cntry       (Intercept)            gndr.c -0.47 -0.01
## 6     cntry       (Intercept)        essround.c -0.05  0.00
## 7     cntry       (Intercept) gndr.c:essround.c  0.16  0.00
## 8     cntry            gndr.c        essround.c  0.05  0.00
## 9     cntry            gndr.c gndr.c:essround.c -0.07  0.00
## 10    cntry        essround.c gndr.c:essround.c -0.24  0.00
## 11 Residual              <NA>              <NA>  1.00  1.00
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_GDI,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] 15.19882
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_GDI,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] -6.010393
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
##         4.5  0.09856 0.0692 Inf  -0.03705    0.2342   1.424  0.1543
##        -4.5 -0.00774 0.0718 Inf  -0.14839    0.1329  -0.108  0.9141
## 
## gdi.z.cm =  0:
##  essround.c   emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##         4.5  0.02999 0.0635 Inf  -0.09442    0.1544   0.472  0.6366
##        -4.5  0.06083 0.0662 Inf  -0.06888    0.1905   0.919  0.3580
## 
## gdi.z.cm =  1:
##  essround.c   emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##         4.5 -0.03858 0.0687 Inf  -0.17323    0.0961  -0.562  0.5744
##        -4.5  0.12940 0.0711 Inf  -0.00999    0.2688   1.819  0.0688
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
##  essround.c4.5 - (essround.c-4.5)   0.1063 0.0757 Inf    -0.042    0.2546   1.405  0.1601
## 
## gdi.z.cm =  0:
##  contrast                         estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  essround.c4.5 - (essround.c-4.5)  -0.0308 0.0517 Inf    -0.132    0.0705  -0.597  0.5508
## 
## gdi.z.cm =  1:
##  contrast                         estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  essround.c4.5 - (essround.c-4.5)  -0.1680 0.0735 Inf    -0.312   -0.0238  -2.284  0.0224
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
##    -0.5        4.5  0.0229 0.0720 Inf   -0.1183    0.1640   0.318  0.7508
##     0.5        4.5  0.1742 0.0679 Inf    0.0413    0.3072   2.568  0.0102
##    -0.5       -4.5 -0.1328 0.0754 Inf   -0.2806    0.0151  -1.760  0.0784
##     0.5       -4.5  0.1173 0.0693 Inf   -0.0185    0.2531   1.692  0.0906
## 
## gdi.z.cm =  0:
##  gndr.c essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##    -0.5        4.5 -0.0410 0.0658 Inf   -0.1700    0.0880  -0.623  0.5335
##     0.5        4.5  0.1010 0.0618 Inf   -0.0201    0.2220   1.634  0.1022
##    -0.5       -4.5 -0.0527 0.0695 Inf   -0.1888    0.0835  -0.758  0.4483
##     0.5       -4.5  0.1743 0.0635 Inf    0.0498    0.2988   2.744  0.0061
## 
## gdi.z.cm =  1:
##  gndr.c essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##    -0.5        4.5 -0.1048 0.0715 Inf   -0.2449    0.0352  -1.467  0.1423
##     0.5        4.5  0.0277 0.0674 Inf   -0.1043    0.1597   0.411  0.6811
##    -0.5       -4.5  0.0274 0.0748 Inf   -0.1192    0.1740   0.367  0.7138
##     0.5       -4.5  0.2314 0.0687 Inf    0.0967    0.3661   3.366  0.0008
## 
## Degrees-of-freedom method: asymptotic 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod6_GDI,adjust="none",infer=c(T,T))
```

```
## gdi.z.cm = -1:
##  contrast                                                  estimate     SE  df asymp.LCL asymp.UCL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5     -0.151380 0.0207 Inf  -0.19198   -0.1108
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  0.155639 0.0788 Inf   0.00123    0.3100
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)  -0.094409 0.0767 Inf  -0.24469    0.0559
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)     0.307019 0.0778 Inf   0.15450    0.4595
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      0.056971 0.0748 Inf  -0.08955    0.2035
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5) -0.250048 0.0198 Inf  -0.28883   -0.2113
##  z.ratio p.value
##   -7.309  <.0001
##    1.976  0.0482
##   -1.231  0.2182
##    3.945  0.0001
##    0.762  0.4460
##  -12.636  <.0001
## 
## gdi.z.cm =  0:
##  contrast                                                  estimate     SE  df asymp.LCL asymp.UCL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5     -0.141950 0.0134 Inf  -0.16813   -0.1158
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  0.011685 0.0537 Inf  -0.09356    0.1169
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)  -0.215306 0.0523 Inf  -0.31789   -0.1127
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)     0.153635 0.0534 Inf   0.04902    0.2583
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)     -0.073356 0.0509 Inf  -0.17313    0.0264
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5) -0.226991 0.0142 Inf  -0.25477   -0.1992
##  z.ratio p.value
##  -10.628  <.0001
##    0.218  0.8277
##   -4.113  <.0001
##    2.878  0.0040
##   -1.441  0.1496
##  -16.016  <.0001
## 
## gdi.z.cm =  1:
##  contrast                                                  estimate     SE  df asymp.LCL asymp.UCL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5     -0.132521 0.0202 Inf  -0.17215   -0.0929
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5) -0.132269 0.0766 Inf  -0.28231    0.0178
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)  -0.336203 0.0745 Inf  -0.48230   -0.1901
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)     0.000252 0.0756 Inf  -0.14802    0.1485
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)     -0.203682 0.0728 Inf  -0.34643   -0.0609
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5) -0.203934 0.0201 Inf  -0.24326   -0.1646
##  z.ratio p.value
##   -6.554  <.0001
##   -1.728  0.0840
##   -4.510  <.0001
##    0.003  0.9973
##   -2.797  0.0052
##  -10.163  <.0001
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
##  diff_ESS10    0.151 0.0207 Inf    0.1108     0.192   7.309  <.0001
##  diff_ESS1     0.250 0.0198 Inf    0.2113     0.289  12.636  <.0001
## 
## gdi.z.cm =  0:
##  contrast   estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10    0.142 0.0134 Inf    0.1158     0.168  10.628  <.0001
##  diff_ESS1     0.227 0.0142 Inf    0.1992     0.255  16.016  <.0001
## 
## gdi.z.cm =  1:
##  contrast   estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10    0.133 0.0202 Inf    0.0929     0.172   6.554  <.0001
##  diff_ESS1     0.204 0.0201 Inf    0.1646     0.243  10.163  <.0001
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
##  diff_ESS10 - diff_ESS1  -0.0987 0.0261 Inf    -0.150   -0.0475  -3.783  0.0002
## 
## gdi.z.cm =  0:
##  contrast               estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10 - diff_ESS1  -0.0850 0.0163 Inf    -0.117   -0.0531  -5.219  <.0001
## 
## gdi.z.cm =  1:
##  contrast               estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10 - diff_ESS1  -0.0714 0.0264 Inf    -0.123   -0.0197  -2.706  0.0068
## 
## Degrees-of-freedom method: asymptotic 
## Confidence level used: 0.95
```

#### Country-level predictors: log_GDP


``` r
mod6_log_GDP<-lmer(ach.z~gndr.c+essround.c+
             gndr.c:essround.c+
               log_gdp.z.cm:gndr.c+log_gdp.z.cm:essround.c+log_gdp.z.cm:gndr.c:essround.c+
               (gndr.c+essround.c+gndr.c:essround.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
```

```
## Warning in checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv, : Model failed to converge
## with max|grad| = 0.00314846 (tol = 0.002, component 1)
```

``` r
summary(mod6_log_GDP)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ach.z ~ gndr.c + essround.c + gndr.c:essround.c + log_gdp.z.cm:gndr.c +  
##     log_gdp.z.cm:essround.c + log_gdp.z.cm:gndr.c:essround.c +  
##     (gndr.c + essround.c + gndr.c:essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1308769.2 1308967.2 -654366.6 1308733.2    441150 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.1303 -0.6696  0.0314  0.6616  5.1190 
## 
## Random effects:
##  Groups   Name              Variance  Std.Dev. Corr             
##  cntry    (Intercept)       0.1154732 0.339813                  
##           gndr.c            0.0036103 0.060086 -0.38            
##           essround.c        0.0012532 0.035400 -0.06  0.15      
##           gndr.c:essround.c 0.0000506 0.007113  0.37 -0.19 -0.34
##  Residual                   1.0015171 1.000758                  
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##                                 Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)                     0.044351   0.059286 33.094621   0.748   0.4597    
## gndr.c                          0.184351   0.011073 31.748365  16.649  < 2e-16 ***
## essround.c                     -0.003535   0.006272 27.212086  -0.564   0.5777    
## gndr.c:essround.c              -0.009865   0.001800 18.510646  -5.480 3.01e-05 ***
## gndr.c:log_gdp.z.cm             0.005430   0.010566 35.104485   0.514   0.6105    
## essround.c:log_gdp.z.cm        -0.005400   0.006426 29.417741  -0.840   0.4075    
## gndr.c:essround.c:log_gdp.z.cm  0.003955   0.001878 21.392528   2.105   0.0472 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c essrn. gnd.:. g.:_.. e.:_..
## gndr.c      -0.357                                   
## essround.c  -0.067  0.139                            
## gndr.c:ssr.  0.254 -0.143 -0.232                     
## gndr.c:l_..  0.001 -0.022  0.000  0.003              
## essrnd.:_..  0.010 -0.004 -0.018  0.004  0.122       
## gndr.:.:_.. -0.003  0.006  0.001 -0.206 -0.038 -0.211
## optimizer (bobyqa) convergence code: 0 (OK)
## Model failed to converge with max|grad| = 0.00314846 (tol = 0.002, component 1)
```

``` r
getFE(mod6_log_GDP,round=2,p.round = 5)
```

```
##                                 Est.   SE    df     t       p    LL    UL
## (Intercept)                     0.04 0.06 33.09  0.75 0.45969 -0.08  0.16
## gndr.c                          0.18 0.01 31.75 16.65 0.00000  0.16  0.21
## essround.c                      0.00 0.01 27.21 -0.56 0.57766 -0.02  0.01
## gndr.c:essround.c              -0.01 0.00 18.51 -5.48 0.00003 -0.01 -0.01
## gndr.c:log_gdp.z.cm             0.01 0.01 35.10  0.51 0.61050 -0.02  0.03
## essround.c:log_gdp.z.cm        -0.01 0.01 29.42 -0.84 0.40750 -0.02  0.01
## gndr.c:essround.c:log_gdp.z.cm  0.00 0.00 21.39  2.11 0.04725  0.00  0.01
```

``` r
getVC(mod6_log_GDP)
```

```
##         grp              var1              var2 sdcor  vcov
## 1     cntry       (Intercept)              <NA>  0.34  0.12
## 2     cntry            gndr.c              <NA>  0.06  0.00
## 3     cntry        essround.c              <NA>  0.04  0.00
## 4     cntry gndr.c:essround.c              <NA>  0.01  0.00
## 5     cntry       (Intercept)            gndr.c -0.38 -0.01
## 6     cntry       (Intercept)        essround.c -0.06  0.00
## 7     cntry       (Intercept) gndr.c:essround.c  0.37  0.00
## 8     cntry            gndr.c        essround.c  0.15  0.00
## 9     cntry            gndr.c gndr.c:essround.c -0.19  0.00
## 10    cntry        essround.c gndr.c:essround.c -0.34  0.00
## 11 Residual              <NA>              <NA>  1.00  1.00
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_log_GDP,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] -1.722099
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_log_GDP,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] 2.440637
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
##         4.5 0.05274 0.0701 Inf   -0.0847     0.190   0.752  0.4520
##        -4.5 0.03596 0.0737 Inf   -0.1086     0.180   0.488  0.6258
## 
## log_gdp.z.cm =  0:
##  essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##         4.5 0.02844 0.0639 Inf   -0.0968     0.154   0.445  0.6563
##        -4.5 0.06026 0.0674 Inf   -0.0718     0.192   0.895  0.3710
## 
## log_gdp.z.cm =  1:
##  essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##         4.5 0.00415 0.0702 Inf   -0.1334     0.142   0.059  0.9529
##        -4.5 0.08456 0.0729 Inf   -0.0583     0.227   1.160  0.2459
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
##  essround.c4.5 - (essround.c-4.5)   0.0168 0.0815 Inf    -0.143    0.1766   0.206  0.8369
## 
## log_gdp.z.cm =  0:
##  contrast                         estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  essround.c4.5 - (essround.c-4.5)  -0.0318 0.0565 Inf    -0.142    0.0788  -0.564  0.5730
## 
## log_gdp.z.cm =  1:
##  contrast                         estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  essround.c4.5 - (essround.c-4.5)  -0.0804 0.0801 Inf    -0.237    0.0765  -1.004  0.3153
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
##    -0.5        4.5 -0.00562 0.0717 Inf   -0.1462    0.1350  -0.078  0.9375
##     0.5        4.5  0.11111 0.0698 Inf   -0.0257    0.2480   1.591  0.1115
##    -0.5       -4.5 -0.08460 0.0781 Inf   -0.2376    0.0684  -1.084  0.2784
##     0.5       -4.5  0.15652 0.0707 Inf    0.0179    0.2951   2.213  0.0269
## 
## log_gdp.z.cm =  0:
##  gndr.c essround.c   emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##    -0.5        4.5 -0.04153 0.0652 Inf   -0.1693    0.0862  -0.637  0.5240
##     0.5        4.5  0.09842 0.0633 Inf   -0.0256    0.2225   1.555  0.1199
##    -0.5       -4.5 -0.05411 0.0710 Inf   -0.1933    0.0851  -0.762  0.4461
##     0.5       -4.5  0.17463 0.0643 Inf    0.0486    0.3007   2.715  0.0066
## 
## log_gdp.z.cm =  1:
##  gndr.c essround.c   emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##    -0.5        4.5 -0.07745 0.0717 Inf   -0.2180    0.0631  -1.080  0.2801
##     0.5        4.5  0.08574 0.0698 Inf   -0.0510    0.2225   1.229  0.2190
##    -0.5       -4.5 -0.02363 0.0771 Inf   -0.1747    0.1274  -0.307  0.7591
##     0.5       -4.5  0.19275 0.0697 Inf    0.0561    0.3294   2.764  0.0057
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
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.1167 0.0192 Inf   -0.1545  -0.07900
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   0.0790 0.0851 Inf   -0.0878   0.24572
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.1621 0.0810 Inf   -0.3209  -0.00337
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.1957 0.0849 Inf    0.0292   0.36220
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      -0.0454 0.0799 Inf   -0.2021   0.11129
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.2411 0.0209 Inf   -0.2822  -0.20006
##  z.ratio p.value
##   -6.064  <.0001
##    0.928  0.3533
##   -2.002  0.0453
##    2.304  0.0212
##   -0.568  0.5701
##  -11.513  <.0001
## 
## log_gdp.z.cm =  0:
##  contrast                                                 estimate     SE  df asymp.LCL asymp.UCL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.1400 0.0128 Inf   -0.1650  -0.11496
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   0.0126 0.0589 Inf   -0.1028   0.12794
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.2162 0.0560 Inf   -0.3259  -0.10641
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.1525 0.0590 Inf    0.0369   0.26820
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      -0.0762 0.0551 Inf   -0.1843   0.03186
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.2287 0.0146 Inf   -0.2574  -0.20009
##  z.ratio p.value
##  -10.974  <.0001
##    0.214  0.8308
##   -3.860  0.0001
##    2.585  0.0097
##   -1.382  0.1669
##  -15.643  <.0001
## 
## log_gdp.z.cm =  1:
##  contrast                                                 estimate     SE  df asymp.LCL asymp.UCL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.1632 0.0175 Inf   -0.1975  -0.12882
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  -0.0538 0.0833 Inf   -0.2170   0.10941
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.2702 0.0795 Inf   -0.4261  -0.11433
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.1094 0.0834 Inf   -0.0541   0.27287
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      -0.1070 0.0781 Inf   -0.2602   0.04615
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.2164 0.0192 Inf   -0.2540  -0.17874
##  z.ratio p.value
##   -9.308  <.0001
##   -0.646  0.5182
##   -3.398  0.0007
##    1.311  0.1898
##   -1.369  0.1709
##  -11.267  <.0001
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
##  diff_ESS10    0.117 0.0192 Inf     0.079     0.154   6.064  <.0001
##  diff_ESS1     0.241 0.0209 Inf     0.200     0.282  11.513  <.0001
## 
## log_gdp.z.cm =  0:
##  contrast   estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10    0.140 0.0128 Inf     0.115     0.165  10.974  <.0001
##  diff_ESS1     0.229 0.0146 Inf     0.200     0.257  15.643  <.0001
## 
## log_gdp.z.cm =  1:
##  contrast   estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10    0.163 0.0175 Inf     0.129     0.198   9.308  <.0001
##  diff_ESS1     0.216 0.0192 Inf     0.179     0.254  11.267  <.0001
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
##  diff_ESS10 - diff_ESS1  -0.1244 0.0257 Inf   -0.1748   -0.0740  -4.838  <.0001
## 
## log_gdp.z.cm =  0:
##  contrast               estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10 - diff_ESS1  -0.0888 0.0162 Inf   -0.1205   -0.0570  -5.480  <.0001
## 
## log_gdp.z.cm =  1:
##  contrast               estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10 - diff_ESS1  -0.0532 0.0209 Inf   -0.0941   -0.0123  -2.549  0.0108
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

