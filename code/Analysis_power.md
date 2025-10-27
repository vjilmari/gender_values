---
title: "Analysis for Examining the Gender Equality Paradox in Values Using power Value"
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
cntry.pow<-diff_dat %>% group_by(cntry,essround) %>%
  summarise(pow.ctm=mean(pow,na.rm=T),
            pow.ctsd=sd(pow,na.rm=T)) %>%
  group_by(cntry) %>%
  summarise(pow.cm=mean(pow.ctm),
            pow.csd=mean(pow.ctsd)) 
```

```
## `summarise()` has grouped output by 'cntry'. You can override using the `.groups` argument.
```

``` r
grand_mean_ben<-mean(cntry.pow$pow.cm)
grand_sd_ben<-mean(cntry.pow$pow.csd)

# standardized
diff_dat$pow.z<-(diff_dat$pow-grand_mean_ben)/grand_sd_ben
hist(diff_dat$pow.z)
```

![](Analysis_power_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

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

# value-based power

cntry_ben_frame<-
  diff_dat %>% group_by(cntry,essround) %>%
  summarise('pow M' = weighted.mean(x=pow.z,w=pspwght),
            'pow SD' = sqrt(wtd.var(pow.z,pspwght))) %>%
  group_by(cntry) %>%
  summarise('pow M' = mean(x=`pow M`),
            'pow SD'= mean(x=`pow SD`))
```

```
## `summarise()` has grouped output by 'cntry'. You can override using the `.groups` argument.
```

``` r
cntry_ben_women_frame<-
  diff_dat %>%
  filter(gndr.c==-0.5) %>%
  group_by(cntry,essround) %>%
  summarise('pow M' = weighted.mean(x=pow.z,w=pspwght),
            'pow SD' = sqrt(wtd.var(pow.z,pspwght))) %>%
  group_by(cntry) %>%
  summarise('pow M Women' = mean(x=`pow M`),
            'pow SD Women'= mean(x=`pow SD`))
```

```
## `summarise()` has grouped output by 'cntry'. You can override using the `.groups` argument.
```

``` r
cntry_ben_men_frame<-
  diff_dat %>%
  filter(gndr.c==0.5) %>%
  group_by(cntry,essround) %>%
  summarise('pow M' = weighted.mean(x=pow.z,w=pspwght),
            'pow SD' = sqrt(wtd.var(pow.z,pspwght))) %>%
  group_by(cntry) %>%
  summarise('pow M Men' = mean(x=`pow M`),
            'pow SD Men'= mean(x=`pow SD`))
```

```
## `summarise()` has grouped output by 'cntry'. You can override using the `.groups` argument.
```

``` r
# link n and pow datasets

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
desc_frame$D<-desc_frame$`pow M Men`-desc_frame$`pow M Women`

desc_frame
```

```
## # A tibble: 33 × 10
##    cntry `n ESS rounds`     n `pow M` `pow SD` `pow M Women` `pow SD Women` `pow M Men` `pow SD Men`
##    <chr>          <dbl> <dbl>   <dbl>    <dbl>         <dbl>          <dbl>       <dbl>        <dbl>
##  1 AT                 6 13077  0.163     0.975        0.0298          0.983      0.306         0.945
##  2 BE                10 17313 -0.156     0.942       -0.244           0.941     -0.0638        0.934
##  3 BG                 6 12641 -0.0341    1.08        -0.157           1.08       0.0979        1.07 
##  4 CH                10 16720 -0.0511    0.973       -0.167           0.980      0.0706        0.950
##  5 CY                 5  5105  0.205     1.02         0.132           1.02       0.283         1.01 
##  6 CZ                 9 18934  0.271     1.04         0.135           1.06       0.418         1.00 
##  7 DE                 9 25389 -0.287     0.990       -0.437           0.968     -0.130         0.988
##  8 DK                 8 12198 -0.212     0.982       -0.359           0.962     -0.0599        0.979
##  9 EE                 9 16692 -0.364     1.04        -0.496           1.04      -0.206         1.03 
## 10 ES                 9 16954 -0.239     1.06        -0.334           1.05      -0.140         1.05 
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
    `pow M`, `pow SD`,
    `pow M Women`, `pow SD Women`,
    `pow M Men`, `pow SD Men`,
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
##    Country    `n ESS rounds`     n `pow M` `pow SD` `pow M Women` `pow SD Women` `pow M Men` `pow SD Men`
##    <chr>               <dbl> <dbl> <chr>   <chr>    <chr>         <chr>          <chr>       <chr>       
##  1 Austria                 6 13077 0.16    0.98     0.03          0.98           0.31        0.94        
##  2 Belgium                10 17313 -0.16   0.94     -0.24         0.94           -0.06       0.93        
##  3 Bulgaria                6 12641 -0.03   1.08     -0.16         1.08           0.10        1.07        
##  4 Switzerla…             10 16720 -0.05   0.97     -0.17         0.98           0.07        0.95        
##  5 Cyprus                  5  5105 0.21    1.02     0.13          1.02           0.28        1.01        
##  6 Czechia                 9 18934 0.27    1.04     0.13          1.06           0.42        1.00        
##  7 Germany                 9 25389 -0.29   0.99     -0.44         0.97           -0.13       0.99        
##  8 Denmark                 8 12198 -0.21   0.98     -0.36         0.96           -0.06       0.98        
##  9 Estonia                 9 16692 -0.36   1.04     -0.50         1.04           -0.21       1.03        
## 10 Spain                   9 16954 -0.24   1.06     -0.33         1.05           -0.14       1.05        
## 11 Finland                10 18050 -0.60   0.98     -0.78         0.96           -0.40       0.97        
## 12 France                 10 18720 -0.63   1.02     -0.74         1.01           -0.51       1.02        
## 13 UK                     10 21456 -0.24   1.02     -0.38         1.01           -0.10       1.02        
## 14 Greece                  5 12464 0.60    1.02     0.52          1.03           0.68        1.01        
## 15 Croatia                 4  6368 -0.01   1.03     -0.10         1.04           0.08        1.01        
## 16 Hungary                10 16006 0.26    1.07     0.18          1.09           0.35        1.05        
## 17 Ireland                10 20576 0.00    1.07     -0.07         1.07           0.07        1.06        
## 18 Israel                  6 13964 0.49    1.13     0.37          1.15           0.61        1.09        
## 19 Iceland                 5  3832 -0.46   0.90     -0.60         0.87           -0.32       0.90        
## 20 Italy                   4  8663 0.26    0.93     0.23          0.93           0.30        0.94        
## 21 Lithuania               6 11714 0.26    1.03     0.17          1.03           0.37        1.02        
## 22 Latvia                  2  2866 0.24    1.04     0.15          1.02           0.34        1.05        
## 23 Montenegro              2  2441 0.20    1.05     0.12          1.08           0.29        1.01        
## 24 Netherlan…             10 18048 -0.25   0.91     -0.41         0.91           -0.09       0.89        
## 25 Norway                 10 15186 -0.36   0.91     -0.46         0.90           -0.26       0.90        
## 26 Poland                  9 15314 0.21    1.02     0.08          1.03           0.35        0.98        
## 27 Portugal               10 17705 -0.18   0.91     -0.27         0.91           -0.09       0.89        
## 28 Russia                  5 12139 0.62    1.05     0.57          1.04           0.67        1.05        
## 29 Sweden                  9 14897 -0.40   0.94     -0.52         0.92           -0.27       0.95        
## 30 Slovenia               10 13238 0.10    0.92     -0.01         0.92           0.21        0.91        
## 31 Slovakia                7 11132 0.29    1.02     0.17          1.02           0.43        0.99        
## 32 Turkey                  2  4108 0.82    0.95     0.73          0.97           0.92        0.92        
## 33 Ukraine                 5  9454 0.36    1.12     0.29          1.12           0.44        1.11        
## # ℹ 5 more variables: D <chr>, GEI <chr>, GGGI <chr>, GDI <chr>, GDP <chr>
```

``` r
export(cntry_desc_tbl,"../results/pow/cntry_desc_tbl.xlsx",overwrite=T)
```

## Country-level correlations table


``` r
cor_frame<-
  desc_frame %>%
  select(
    VBMT=`pow M`,
    VBMT_Women=`pow M Women`,
    VBMT_Men=`pow M Men`,
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
  filename = "../results/pow/CorTable1.doc",
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
##   1. VBMT       0.03  0.36                                                                           
##                                                                                                      
##   2. VBMT_Women -0.08 0.38 1.00                                                                      
##                            [.99, 1.00]                                                               
##                                                                                                      
##   3. VBMT_Men   0.14  0.34 1.00         .99                                                          
##                            [.99, 1.00]  [.97, .99]                                                   
##                                                                                                      
##   4. D          0.22  0.07 -.55         -.60         -.47                                            
##                            [-.75, -.25] [-.79, -.33] [-.70, -.15]                                    
##                                                                                                      
##   5. GEI        0.87  0.07 -.69         -.69         -.69         .40                                
##                            [-.84, -.45] [-.84, -.46] [-.84, -.45] [.06, .66]                         
##                                                                                                      
##   6. GGGI       0.73  0.05 -.79         -.79         -.78         .45         .73                    
##                            [-.89, -.61] [-.89, -.61] [-.89, -.60] [.13, .69]  [.52, .86]             
##                                                                                                      
##   7. GDI        0.99  0.03 -.11         -.09         -.11         -.06        .07         .20        
##                            [-.43, .25]  [-.42, .26]  [-.44, .24]  [-.40, .29] [-.29, .41] [-.16, .51]
##                                                                                                      
##   8. log_GDP    10.62 0.40 -.58         -.58         -.58         .32         .75         .67        
##                            [-.77, -.29] [-.77, -.30] [-.77, -.30] [-.03, .59] [.55, .87]  [.42, .82] 
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
mod0<-lmer(pow.z~(1|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod0)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: pow.z ~ (1 | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1323689.1 1323722.1 -661841.6 1323683.1    441165 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -5.9151 -0.6686 -0.0465  0.6239  6.6491 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  cntry    (Intercept) 0.1297   0.3601  
##  Residual             1.0366   1.0181  
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##             Estimate Std. Error       df t value Pr(>|t|)
## (Intercept)  0.03913    0.06271 32.99356   0.624    0.537
```

``` r
getVC(mod0)
```

```
##        grp        var1 var2 sdcor vcov
## 1    cntry (Intercept) <NA>  0.36 0.13
## 2 Residual        <NA> <NA>  1.02 1.04
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
## mean variation  0.1111749     NA       1
## sigma2          0.8888251      1      NA
## 
## $R2s
##         total within between
## f1  0.0000000      0      NA
## f2  0.0000000     NA       0
## v   0.0000000      0      NA
## m   0.1111749     NA       1
## f   0.0000000     NA      NA
## fv  0.0000000      0      NA
## fvm 0.1111749     NA      NA
```

## mod1: Gender fixed effect


``` r
mod1<-lmer(pow.z~gndr.c+(1|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod1)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: pow.z ~ gndr.c + (1 | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1317804.4 1317848.4 -658898.2 1317796.4    441164 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.2086 -0.6604 -0.0039  0.6165  6.4319 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  cntry    (Intercept) 0.1303   0.361   
##  Residual             1.0229   1.011   
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##              Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept) 4.359e-02  6.288e-02 3.299e+01   0.693    0.493    
## gndr.c      2.341e-01  3.041e-03 4.411e+05  76.982   <2e-16 ***
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
##              Est.    SE         df      t     p     LL    UL
## (Intercept) 0.044 0.063     32.992  0.693 0.493 -0.084 0.172
## gndr.c      0.234 0.003 441135.439 76.982 0.000  0.228 0.240
```

``` r
getVC(mod1)
```

```
##        grp        var1 var2 sdcor vcov
## 1    cntry (Intercept) <NA>  0.36 0.13
## 2 Residual        <NA> <NA>  1.01 1.02
```

``` r
r2mlm(mod1,bargraph = F)
```

```
## $Decompositions
##                      total
## fixed           0.01167079
## slope variation 0.00000000
## mean variation  0.11171314
## sigma2          0.87661607
## 
## $R2s
##          total
## f   0.01167079
## v   0.00000000
## m   0.11171314
## fv  0.01167079
## fvm 0.12338393
```

## mod2: Gender fixed and random effect

* Include random effect correlation by default


``` r
mod2<-lmer(pow.z~gndr.c+(gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: pow.z ~ gndr.c + (gndr.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##   1317430   1317496   -658709   1317418    441162 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.0776 -0.6584 -0.0132  0.6176  6.4333 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr 
##  cntry    (Intercept) 0.129716 0.36016       
##           gndr.c      0.004179 0.06465  -0.56
##  Residual             1.021807 1.01084       
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##             Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)  0.04318    0.06272 32.99086   0.688    0.496    
## gndr.c       0.22346    0.01178 33.32784  18.973   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##        (Intr)
## gndr.c -0.533
```

``` r
getFE(mod2,round=3)
```

```
##              Est.    SE     df      t     p     LL    UL
## (Intercept) 0.043 0.063 32.991  0.688 0.496 -0.084 0.171
## gndr.c      0.223 0.012 33.328 18.973 0.000  0.200 0.247
```

``` r
getVC(mod2)
```

```
##        grp        var1   var2 sdcor  vcov
## 1    cntry (Intercept)   <NA>  0.36  0.13
## 2    cntry      gndr.c   <NA>  0.06  0.00
## 3    cntry (Intercept) gndr.c -0.56 -0.01
## 4 Residual        <NA>   <NA>  1.01  1.02
```

``` r
r2mlm(mod2,bargraph = F)
```

```
## $Decompositions
##                        total
## fixed           0.0106414210
## slope variation 0.0008906097
## mean variation  0.1121233746
## sigma2          0.8763445947
## 
## $R2s
##            total
## f   0.0106414210
## v   0.0008906097
## m   0.1121233746
## fv  0.0115320307
## fvm 0.1236554053
```

``` r
anova(mod1,mod2)
```

```
## Data: diff_dat
## Models:
## mod1: pow.z ~ gndr.c + (1 | cntry)
## mod2: pow.z ~ gndr.c + (gndr.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L) Chisq Df Pr(>Chisq)    
## mod1    4 1317804 1317848 -658898   1317796                        
## mod2    6 1317430 1317496 -658709   1317418 378.4  2  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

``` r
lvl2_var_cond_lvl1(mod2,lvl1.var = "gndr.c",lvl1.values = c(0.5,-0.5))
```

```
##   lvl1.value lvl2.cond.var lvl2.cond.sd
## 1        0.5     0.1177492    0.3431460
## 2       -0.5     0.1437731    0.3791742
```

* Test for random effect correlation


``` r
mod2_norecov<-lmer(pow.z~gndr.c+(gndr.c||cntry),data=diff_dat,REML=F,weights = pspwght,
                   control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2_norecov)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: pow.z ~ gndr.c + (gndr.c || cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1317438.9 1317493.9 -658714.5 1317428.9    441163 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.0801 -0.6584 -0.0131  0.6174  6.4388 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  cntry    (Intercept) 0.129760 0.36022 
##  cntry.1  gndr.c      0.004137 0.06432 
##  Residual             1.021808 1.01085 
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##             Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)  0.04319    0.06273 32.99152   0.689    0.496    
## gndr.c       0.22409    0.01173 33.42924  19.098   <2e-16 ***
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
## (Intercept) 0.043 0.063 32.992  0.689 0.496 -0.084 0.171
## gndr.c      0.224 0.012 33.429 19.098 0.000  0.200 0.248
```

``` r
getVC(mod2_norecov)
```

```
##        grp        var1 var2 sdcor vcov
## 1    cntry (Intercept) <NA>  0.36 0.13
## 2  cntry.1      gndr.c <NA>  0.06 0.00
## 3 Residual        <NA> <NA>  1.01 1.02
```

``` r
anova(mod2_norecov,mod2)
```

```
## Data: diff_dat
## Models:
## mod2_norecov: pow.z ~ gndr.c + (gndr.c || cntry)
## mod2: pow.z ~ gndr.c + (gndr.c | cntry)
##              npar     AIC     BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod2_norecov    5 1317439 1317494 -658714   1317429                         
## mod2            6 1317430 1317496 -658709   1317418 10.912  1  0.0009554 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```


## mod2 with Gender-equality index (GEI)


``` r
mod2_GEI<-lmer(pow.z~gndr.c+gei.z.cm+gndr.c:gei.z.cm+
                 (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2_GEI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: pow.z ~ gndr.c + gei.z.cm + gndr.c:gei.z.cm + (gndr.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1285142.1 1285229.9 -642563.1 1285126.1    431770 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.0937 -0.6599 -0.0133  0.6189  6.4500 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr 
##  cntry    (Intercept) 0.064683 0.25433       
##           gndr.c      0.003428 0.05855  -0.39
##  Residual             1.016544 1.00824       
## Number of obs: 431778, groups:  cntry, 32
## 
## Fixed effects:
##                 Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)      0.03315    0.04500 32.01935   0.737   0.4667    
## gndr.c           0.22579    0.01094 32.45974  20.643  < 2e-16 ***
## gei.z.cm        -0.26063    0.04573 32.06333  -5.699 2.58e-06 ***
## gndr.c:gei.z.cm  0.02704    0.01136 35.26717   2.381   0.0228 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c g.z.cm
## gndr.c      -0.369              
## gei.z.cm    -0.001  0.000       
## gndr.c:g.z.  0.000 -0.032 -0.361
```

``` r
getFE(mod2_GEI,round=3)
```

```
##                   Est.    SE     df      t     p     LL     UL
## (Intercept)      0.033 0.045 32.019  0.737 0.467 -0.059  0.125
## gndr.c           0.226 0.011 32.460 20.643 0.000  0.204  0.248
## gei.z.cm        -0.261 0.046 32.063 -5.699 0.000 -0.354 -0.167
## gndr.c:gei.z.cm  0.027 0.011 35.267  2.381 0.023  0.004  0.050
```

``` r
getVC(mod2_GEI)
```

```
##        grp        var1   var2 sdcor  vcov
## 1    cntry (Intercept)   <NA>  0.25  0.06
## 2    cntry      gndr.c   <NA>  0.06  0.00
## 3    cntry (Intercept) gndr.c -0.39 -0.01
## 4 Residual        <NA>   <NA>  1.01  1.02
```

``` r
r2mlm(mod2_GEI,bargraph = F)
```

```
## $Decompositions
##                        total
## fixed           0.0518203314
## slope variation 0.0007465297
## mean variation  0.0570368204
## sigma2          0.8903963185
## 
## $R2s
##            total
## f   0.0518203314
## v   0.0007465297
## m   0.0570368204
## fv  0.0525668612
## fvm 0.1096036815
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
## Time difference of 29.39263 secs
```

``` r
round(ddsc_mod2_GEI$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.144        0.379        1.022     1.166 0.123   7204.152 0.999   0.999
## 2        0.5         0.118        0.343        1.022     1.140 0.103   6164.576 0.998   0.999
```

``` r
round(ddsc_mod2_GEI$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gei.z.cm
## means_y1           0.146 0.351    1.000           1.000    0.988           0.988   -0.706
## means_y1_scaled    0.395 0.952    1.000           1.000    0.988           0.988   -0.706
## means_y2          -0.079 0.386    0.988           0.988    1.000           1.000   -0.710
## means_y2_scaled   -0.216 1.046    0.988           0.988    1.000           1.000   -0.710
## gei.z.cm           0.000 1.000   -0.706          -0.706   -0.710          -0.710    1.000
## gei.z.cm_scaled    0.000 1.000   -0.706          -0.706   -0.710          -0.710    1.000
## diff_score         0.225 0.067   -0.447          -0.447   -0.582          -0.582    0.386
## diff_score_scaled  0.611 0.183   -0.447          -0.447   -0.582          -0.582    0.386
##                   gei.z.cm_scaled diff_score diff_score_scaled
## means_y1                   -0.706     -0.447            -0.447
## means_y1_scaled            -0.706     -0.447            -0.447
## means_y2                   -0.710     -0.582            -0.582
## means_y2_scaled            -0.710     -0.582            -0.582
## gei.z.cm                    1.000      0.386             0.386
## gei.z.cm_scaled             1.000      0.386             0.386
## diff_score                  0.386      1.000             1.000
## diff_score_scaled           0.386      1.000             1.000
```

``` r
round(ddsc_mod2_GEI$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                      -0.401 0.168 35.267  -2.381   0.023   -0.742   -0.059
## w_11                         -0.274 0.048 32.111  -5.702   0.000   -0.372   -0.176
## w_21                         -0.247 0.044 32.138  -5.616   0.000   -0.337   -0.157
## r_xy1                        -0.782 0.137 32.111  -5.702   0.000   -1.061   -0.502
## r_xy2                        -0.641 0.114 32.138  -5.616   0.000   -0.873   -0.408
## b_11                         -0.744 0.131 32.111  -5.702   0.000   -1.010   -0.479
## b_21                         -0.671 0.119 32.138  -5.616   0.000   -0.914   -0.428
## main_effect                  -0.261 0.046 32.063  -5.699   0.000   -0.354   -0.167
## moderator_effect              0.226 0.011 32.460  20.643   0.000    0.204    0.248
## interaction                   0.027 0.011 35.267   2.381   0.023    0.004    0.050
## q_b11_b21                    -0.148    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                  -0.290    NA     NA      NA      NA       NA       NA
## cross_over_point             -8.351    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.234 0.043 32.362  -5.438   0.000   -0.321   -0.146
## interaction_vs_main_bscale   -0.634 0.117 32.362  -5.438   0.000   -0.872   -0.397
## interaction_vs_main_rscale   -0.570 0.106 32.411  -5.401   0.000   -0.785   -0.355
## dadas                        -0.494 0.088 32.138  -5.616   1.000   -0.673   -0.315
## dadas_bscale                 -1.342 0.239 32.138  -5.616   1.000   -1.829   -0.855
## dadas_rscale                 -1.281 0.228 32.138  -5.616   1.000   -1.746   -0.817
## abs_diff                      0.027 0.011 35.267   2.381   0.011    0.004    0.050
## abs_sum                       0.521 0.091 32.063   5.699   0.000    0.335    0.708
## abs_diff_bscale               0.073 0.031 35.267   2.381   0.011    0.011    0.136
## abs_sum_bscale                1.416 0.248 32.063   5.699   0.000    0.910    1.921
## abs_diff_rscale               0.141 0.037 34.396   3.822   0.000    0.066    0.216
## abs_sum_rscale                1.422 0.249 32.063   5.701   0.000    0.914    1.930
```

``` r
round(ddsc_mod2_GEI$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.013 -0.559 10.912  1.000  0.001
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
## r_xy1_y2                        -0.386 0.163  -2.366  0.018   -0.705   -0.066
## r_xy1                           -0.710 0.125  -5.697  0.000   -0.954   -0.465
## r_xy2                           -0.706 0.125  -5.639  0.000   -0.951   -0.461
## b_11                            -0.742 0.130  -5.697  0.000   -0.998   -0.487
## b_21                            -0.672 0.119  -5.639  0.000   -0.905   -0.438
## b_10                            -0.216 0.128  -1.680  0.093   -0.467    0.036
## b_20                             0.395 0.117   3.372  0.001    0.166    0.625
## res_cov_y1_y2                    0.469 0.119   3.950  0.000    0.236    0.702
## diff_b10_b20                    -0.611 0.029 -20.795  0.000   -0.668   -0.553
## diff_b11_b21                    -0.071 0.030  -2.366  0.018   -0.129   -0.012
## diff_rxy1_rxy2                  -0.004 0.028  -0.129  0.897   -0.058    0.051
## q_b11_b21                       -0.142 0.092  -1.541  0.123   -0.322    0.039
## q_rxy1_rxy2                     -0.007 0.056  -0.129  0.897   -0.116    0.102
## cross_over_point                -8.651 3.680  -2.351  0.019  -15.864   -1.438
## sum_b11_b21                     -1.414 0.248  -5.705  0.000   -1.900   -0.928
## main_effect                     -0.707 0.124  -5.705  0.000   -0.950   -0.464
## interaction_vs_main_effect      -0.636 0.116  -5.486  0.000   -0.864   -0.409
## diff_abs_b11_abs_b21             0.071 0.030   2.366  0.018    0.012    0.129
## abs_diff_b11_b21                 0.071 0.030   2.366  0.009    0.012    0.129
## abs_sum_b11_b21                  1.414 0.248   5.705  0.000    0.928    1.900
## dadas                           -1.344 0.238  -5.639  1.000   -1.811   -0.877
## q_r_equivalence                 -0.093 0.056  -1.670  0.047       NA       NA
## q_b_equivalence                  0.042 0.092   0.454  0.675       NA       NA
## cross_over_point_equivalence     8.651 3.680   2.351  0.991       NA       NA
## cross_over_point_minimal_effect  8.651 3.680   2.351  0.009       NA       NA
```

``` r
round(ddsc_sem_GEI$variance_test,3)
```

```
##             est    se       z pvalue ci.lower ci.upper
## cov_y1y2  0.953 0.240   3.975  0.000    0.483    1.422
## var_y1    1.060 0.265   4.000  0.000    0.541    1.580
## var_y2    0.877 0.219   4.000  0.000    0.447    1.307
## var_diff  0.183 0.070   2.606  0.009    0.045    0.321
## var_ratio 1.209 0.067  18.076  0.000    1.078    1.340
## cor_y1y2  0.988 0.004 228.207  0.000    0.979    0.996
```

### Figure 


``` r
# plot the results

# start with obtaining predicted values for means and differences

# refit reduced and full models with GGGI in original scale


mod2_GEI_unstd<-lmer(pow.z~gndr.c+gei.cm+gndr.c:gei.cm+
                 (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

mod2_GEI_unstd_red<-lmer(pow.z~gndr.c+
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


p1.pow.flags<-
  ggplot(p,aes(y=yvar,x=gei.cm,color=gndr.c))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("Gender Equality Index (Average 2002-2022)")+
  #ylim(c(min.y.pred,max.y.pred))+
  #xlim(c(0.60,1.00))+
  ylab("Value power (Average 2002-2022)")+
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

p2.pow.flags<-ggplot(p2,aes(y=yvar,x=gei.cm))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("Gender Equality Index (Average 2002-2022)")+
  #xlim(c(0.60,1.00))+
  #ylim(c(min.y.narrow,max.y.narrow))+
  ylab("Difference in value power")+
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
  ggarrange(p1.pow.flags,p2.pow.flags,align = "v",
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

![](Analysis_power_files/figure-html/unnamed-chunk-18-1.png)<!-- -->

``` r
png(filename = 
      "../results/pow/GEI_flags.png",
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
mod2_GGGI<-lmer(pow.z~gndr.c+gggi.z.cm+gndr.c:gggi.z.cm+
                 (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2_GGGI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: pow.z ~ gndr.c + gggi.z.cm + gndr.c:gggi.z.cm + (gndr.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##  940758.9  940844.2 -470371.5  940742.9    314638 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.0722 -0.6575 -0.0041  0.6212  5.9869 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr 
##  cntry    (Intercept) 0.068952 0.26259       
##           gndr.c      0.003327 0.05768  -0.31
##  Residual             1.023849 1.01185       
## Number of obs: 314646, groups:  cntry, 33
## 
## Fixed effects:
##                  Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)       0.08476    0.04576 32.94140   1.852  0.07298 .  
## gndr.c            0.21856    0.01087 32.36252  20.110  < 2e-16 ***
## gggi.z.cm        -0.30847    0.04650 33.00595  -6.634 1.51e-07 ***
## gndr.c:gggi.z.cm  0.03119    0.01138 36.02077   2.741  0.00947 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c ggg.z.
## gndr.c      -0.282              
## gggi.z.cm    0.000  0.000       
## gndr.c:gg..  0.000 -0.019 -0.274
```

``` r
getFE(mod2_GGGI,round=3)
```

```
##                    Est.    SE     df      t     p     LL     UL
## (Intercept)       0.085 0.046 32.941  1.852 0.073 -0.008  0.178
## gndr.c            0.219 0.011 32.363 20.110 0.000  0.196  0.241
## gggi.z.cm        -0.308 0.046 33.006 -6.634 0.000 -0.403 -0.214
## gndr.c:gggi.z.cm  0.031 0.011 36.021  2.741 0.009  0.008  0.054
```

``` r
getVC(mod2_GGGI)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.26 0.07
## 2    cntry      gndr.c   <NA>  0.06 0.00
## 3    cntry (Intercept) gndr.c -0.31 0.00
## 4 Residual        <NA>   <NA>  1.01 1.02
```

``` r
r2mlm(mod2_GGGI,bargraph = F)
```

```
## $Decompositions
##                        total
## fixed           0.0699227321
## slope variation 0.0007026507
## mean variation  0.0589401721
## sigma2          0.8704344451
## 
## $R2s
##            total
## f   0.0699227321
## v   0.0007026507
## m   0.0589401721
## fv  0.0706253828
## fvm 0.1295655549
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
## Time difference of 29.55208 secs
```

``` r
round(ddsc_mod2_GGGI$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.144        0.379        1.022     1.166 0.123   7204.152 0.999   0.999
## 2        0.5         0.118        0.343        1.022     1.140 0.103   6164.576 0.998   0.999
```

``` r
round(ddsc_mod2_GGGI$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gggi.z.cm
## means_y1           0.194 0.392    1.000           1.000    0.988           0.988    -0.747
## means_y1_scaled    0.473 0.958    1.000           1.000    0.988           0.988    -0.747
## means_y2          -0.024 0.426    0.988           0.988    1.000           1.000    -0.759
## means_y2_scaled   -0.058 1.041    0.988           0.988    1.000           1.000    -0.759
## gggi.z.cm          0.000 1.000   -0.747          -0.747   -0.759          -0.759     1.000
## gggi.z.cm_scaled   0.000 1.000   -0.747          -0.747   -0.759          -0.759     1.000
## diff_score         0.217 0.072   -0.401          -0.401   -0.539          -0.539     0.428
## diff_score_scaled  0.530 0.176   -0.401          -0.401   -0.539          -0.539     0.428
##                   gggi.z.cm_scaled diff_score diff_score_scaled
## means_y1                    -0.747     -0.401            -0.401
## means_y1_scaled             -0.747     -0.401            -0.401
## means_y2                    -0.759     -0.539            -0.539
## means_y2_scaled             -0.759     -0.539            -0.539
## gggi.z.cm                    1.000      0.428             0.428
## gggi.z.cm_scaled             1.000      0.428             0.428
## diff_score                   0.428      1.000             1.000
## diff_score_scaled            0.428      1.000             1.000
```

``` r
round(ddsc_mod2_GGGI$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                      -0.432 0.158 36.021  -2.741   0.009   -0.752   -0.112
## w_11                         -0.324 0.048 33.041  -6.701   0.000   -0.422   -0.226
## w_21                         -0.293 0.045 32.986  -6.469   0.000   -0.385   -0.201
## r_xy1                        -0.826 0.123 33.041  -6.701   0.000   -1.077   -0.575
## r_xy2                        -0.687 0.106 32.986  -6.469   0.000   -0.903   -0.471
## b_11                         -0.792 0.118 33.041  -6.701   0.000   -1.032   -0.551
## b_21                         -0.715 0.111 32.986  -6.469   0.000   -0.940   -0.490
## main_effect                  -0.308 0.046 33.006  -6.634   0.000   -0.403   -0.214
## moderator_effect              0.219 0.011 32.363  20.110   0.000    0.196    0.241
## interaction                   0.031 0.011 36.021   2.741   0.009    0.008    0.054
## q_b11_b21                    -0.177    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                  -0.334    NA     NA      NA      NA       NA       NA
## cross_over_point             -7.007    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.277 0.045 32.967  -6.197   0.000   -0.368   -0.186
## interaction_vs_main_bscale   -0.677 0.109 32.967  -6.197   0.000   -0.900   -0.455
## interaction_vs_main_rscale   -0.617 0.100 32.965  -6.152   0.000   -0.821   -0.413
## dadas                        -0.586 0.091 32.986  -6.469   1.000   -0.770   -0.402
## dadas_bscale                 -1.431 0.221 32.986  -6.469   1.000   -1.881   -0.981
## dadas_rscale                 -1.374 0.212 32.986  -6.469   1.000   -1.806   -0.942
## abs_diff                      0.031 0.011 36.021   2.741   0.005    0.008    0.054
## abs_sum                       0.617 0.093 33.006   6.634   0.000    0.428    0.806
## abs_diff_bscale               0.076 0.028 36.021   2.741   0.005    0.020    0.133
## abs_sum_bscale                1.507 0.227 33.006   6.634   0.000    1.045    1.969
## abs_diff_rscale               0.139 0.032 35.338   4.381   0.000    0.075    0.204
## abs_sum_rscale                1.513 0.228 33.007   6.639   0.000    1.049    1.976
```

``` r
round(ddsc_mod2_GGGI$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.013 -0.559 10.912  1.000  0.001
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
## r_xy1_y2                        -0.428 0.157  -2.722  0.006   -0.737   -0.120
## r_xy1                           -0.759 0.113  -6.705  0.000   -0.981   -0.537
## r_xy2                           -0.747 0.116  -6.449  0.000   -0.974   -0.520
## b_11                            -0.790 0.118  -6.705  0.000   -1.021   -0.559
## b_21                            -0.715 0.111  -6.449  0.000   -0.932   -0.498
## b_10                            -0.058 0.116  -0.496  0.620   -0.285    0.170
## b_20                             0.473 0.109   4.329  0.000    0.259    0.687
## res_cov_y1_y2                    0.407 0.102   4.005  0.000    0.208    0.606
## diff_b10_b20                    -0.530 0.027 -19.443  0.000   -0.584   -0.477
## diff_b11_b21                    -0.075 0.028  -2.722  0.006   -0.130   -0.021
## diff_rxy1_rxy2                  -0.013 0.027  -0.471  0.637   -0.066    0.040
## q_b11_b21                       -0.175 0.107  -1.631  0.103   -0.385    0.035
## q_rxy1_rxy2                     -0.029 0.062  -0.471  0.637   -0.151    0.093
## cross_over_point                -7.034 2.609  -2.696  0.007  -12.147   -1.920
## sum_b11_b21                     -1.505 0.227  -6.626  0.000   -1.951   -1.060
## main_effect                     -0.753 0.114  -6.626  0.000   -0.975   -0.530
## interaction_vs_main_effect      -0.677 0.110  -6.166  0.000   -0.893   -0.462
## diff_abs_b11_abs_b21             0.075 0.028   2.722  0.006    0.021    0.130
## abs_diff_b11_b21                 0.075 0.028   2.722  0.003    0.021    0.130
## abs_sum_b11_b21                  1.505 0.227   6.626  0.000    1.060    1.951
## dadas                           -1.430 0.222  -6.449  1.000   -1.865   -0.995
## q_r_equivalence                 -0.071 0.062  -1.134  0.128       NA       NA
## q_b_equivalence                  0.075 0.107   0.700  0.758       NA       NA
## cross_over_point_equivalence     7.034 2.609   2.696  0.996       NA       NA
## cross_over_point_minimal_effect  7.034 2.609   2.696  0.004       NA       NA
```

``` r
round(ddsc_sem_GGGI$variance_test,3)
```

```
##             est    se       z pvalue ci.lower ci.upper
## cov_y1y2  0.955 0.236   4.037  0.000    0.491    1.418
## var_y1    1.050 0.259   4.062  0.000    0.544    1.557
## var_y2    0.889 0.219   4.062  0.000    0.460    1.318
## var_diff  0.161 0.066   2.462  0.014    0.033    0.290
## var_ratio 1.181 0.064  18.541  0.000    1.057    1.306
## cor_y1y2  0.988 0.004 236.468  0.000    0.980    0.996
```

### Figure


``` r
# plot the results

# start with obtaining predicted values for means and differences

# refit reduced and full models with GGGI in original scale


mod2_GGGI_unstd<-lmer(pow.z~gndr.c+gggi.cm+gndr.c:gggi.cm+
                       (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
                     control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

mod2_GGGI_unstd_red<-lmer(pow.z~gndr.c+
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


p1.pow.flags<-
  ggplot(p,aes(y=yvar,x=gggi.cm,color=gndr.c))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("Global Gender Gap Index (Average 2002-2022)")+
  #ylim(c(min.y.pred,max.y.pred))+
  #xlim(c(0.60,1.00))+
  ylab("Value power (Average 2002-2022)")+
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

p2.pow.flags<-ggplot(p2,aes(y=yvar,x=gggi.cm))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("Global Gender Gap Index (Average 2002-2022)")+
  #xlim(c(0.60,1.00))+
  #ylim(c(min.y.narrow,max.y.narrow))+
  ylab("Difference in value power")+
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
  ggarrange(p1.pow.flags,p2.pow.flags,align = "v",
            ncol=1,nrow=2,heights=c(2,1))
```

```
## Warning in geom_text(inherit.aes = F, aes(x = 0.61, y = -0.45, label = coef_q, : All aesthetics have length 1, but the data has 502 rows.
## ℹ Please consider using `annotate()` or provide this layer with data containing a single row.
```

``` r
pflag_comb
```

![](Analysis_power_files/figure-html/unnamed-chunk-21-1.png)<!-- -->

``` r
png(filename = 
      "../results/pow/GGGI_flags_new.png",
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
mod2_GDI<-lmer(pow.z~gndr.c+gdi.z.cm+gndr.c:gdi.z.cm+
                 (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2_GDI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: pow.z ~ gndr.c + gdi.z.cm + gndr.c:gdi.z.cm + (gndr.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1317433.5 1317521.5 -658708.8 1317417.5    441160 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.0764 -0.6584 -0.0133  0.6175  6.4333 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr 
##  cntry    (Intercept) 0.129031 0.3592        
##           gndr.c      0.004161 0.0645   -0.56
##  Residual             1.021808 1.0108        
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##                  Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)      0.043172   0.062556 32.992571   0.690    0.495    
## gndr.c           0.223506   0.011755 33.273376  19.014   <2e-16 ***
## gdi.z.cm        -0.026507   0.063542 33.024781  -0.417    0.679    
## gndr.c:gdi.z.cm -0.002965   0.012235 36.495667  -0.242    0.810    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c gd.z.c
## gndr.c      -0.538              
## gdi.z.cm     0.000  0.000       
## gndr.c:gd..  0.000 -0.007 -0.524
```

``` r
getFE(mod2_GDI,round=3)
```

```
##                   Est.    SE     df      t     p     LL    UL
## (Intercept)      0.043 0.063 32.993  0.690 0.495 -0.084 0.170
## gndr.c           0.224 0.012 33.273 19.014 0.000  0.200 0.247
## gdi.z.cm        -0.027 0.064 33.025 -0.417 0.679 -0.156 0.103
## gndr.c:gdi.z.cm -0.003 0.012 36.496 -0.242 0.810 -0.028 0.022
```

``` r
getVC(mod2_GDI)
```

```
##        grp        var1   var2 sdcor  vcov
## 1    cntry (Intercept)   <NA>  0.36  0.13
## 2    cntry      gndr.c   <NA>  0.06  0.00
## 3    cntry (Intercept) gndr.c -0.56 -0.01
## 4 Residual        <NA>   <NA>  1.01  1.02
```

``` r
r2mlm(mod2_GDI,bargraph = F)
```

```
## $Decompositions
##                        total
## fixed           0.0112228995
## slope variation 0.0008867402
## mean variation  0.1115392507
## sigma2          0.8763511096
## 
## $R2s
##            total
## f   0.0112228995
## v   0.0008867402
## m   0.1115392507
## fv  0.0121096397
## fvm 0.1236488904
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
## Time difference of 28.95172 secs
```

``` r
round(ddsc_mod2_GDI$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.144        0.379        1.022     1.166 0.123   7204.152 0.999   0.999
## 2        0.5         0.118        0.343        1.022     1.140 0.103   6164.576 0.998   0.999
```

``` r
round(ddsc_mod2_GDI$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gdi.z.cm
## means_y1           0.155 0.349    1.000           1.000    0.988           0.988   -0.079
## means_y1_scaled    0.421 0.950    1.000           1.000    0.988           0.988   -0.079
## means_y2          -0.068 0.385    0.988           0.988    1.000           1.000   -0.066
## means_y2_scaled   -0.186 1.048    0.988           0.988    1.000           1.000   -0.066
## gdi.z.cm           0.000 1.000   -0.079          -0.079   -0.066          -0.066    1.000
## gdi.z.cm_scaled    0.000 1.000   -0.079          -0.079   -0.066          -0.066    1.000
## diff_score         0.223 0.068   -0.463          -0.463   -0.595          -0.595   -0.032
## diff_score_scaled  0.606 0.184   -0.463          -0.463   -0.595          -0.595   -0.032
##                   gdi.z.cm_scaled diff_score diff_score_scaled
## means_y1                   -0.079     -0.463            -0.463
## means_y1_scaled            -0.079     -0.463            -0.463
## means_y2                   -0.066     -0.595            -0.595
## means_y2_scaled            -0.066     -0.595            -0.595
## gdi.z.cm                    1.000     -0.032            -0.032
## gdi.z.cm_scaled             1.000     -0.032            -0.032
## diff_score                 -0.032      1.000             1.000
## diff_score_scaled          -0.032      1.000             1.000
```

``` r
round(ddsc_mod2_GDI$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                       0.044 0.181 36.496   0.242   0.810   -0.322    0.410
## w_11                         -0.025 0.067 33.042  -0.374   0.711   -0.161    0.111
## w_21                         -0.028 0.061 33.055  -0.462   0.647   -0.151    0.095
## r_xy1                        -0.072 0.192 33.042  -0.374   0.711   -0.462    0.319
## r_xy2                        -0.073 0.157 33.055  -0.462   0.647   -0.393    0.247
## b_11                         -0.068 0.182 33.042  -0.374   0.711   -0.439    0.303
## b_21                         -0.076 0.165 33.055  -0.462   0.647   -0.412    0.259
## main_effect                  -0.027 0.064 33.025  -0.417   0.679   -0.156    0.103
## moderator_effect              0.224 0.012 33.273  19.014   0.000    0.200    0.247
## interaction                  -0.003 0.012 36.496  -0.242   0.810   -0.028    0.022
## q_b11_b21                     0.008    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.001    NA     NA      NA      NA       NA       NA
## cross_over_point             75.373    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.024 0.071 33.094  -0.333   0.741   -0.167    0.120
## interaction_vs_main_bscale   -0.064 0.193 33.094  -0.333   0.741   -0.456    0.328
## interaction_vs_main_rscale   -0.071 0.211 33.085  -0.338   0.737   -0.500    0.357
## dadas                        -0.050 0.134 33.042  -0.374   0.645   -0.322    0.222
## dadas_bscale                 -0.136 0.365 33.042  -0.374   0.645   -0.879    0.606
## dadas_rscale                 -0.143 0.384 33.042  -0.374   0.645   -0.924    0.637
## abs_diff                      0.003 0.012 36.496   0.242   0.405   -0.022    0.028
## abs_sum                       0.053 0.127 33.025   0.417   0.340   -0.206    0.312
## abs_diff_bscale               0.008 0.033 36.496   0.242   0.405   -0.059    0.076
## abs_sum_bscale                0.144 0.346 33.025   0.417   0.340   -0.560    0.849
## abs_diff_rscale               0.001 0.045 34.748   0.022   0.491   -0.090    0.092
## abs_sum_rscale                0.144 0.348 33.025   0.415   0.340   -0.563    0.852
```

``` r
round(ddsc_mod2_GDI$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.013 -0.559 10.912  1.000  0.001
```

``` r
d_GDI<-ddsc_mod2_GDI$ddsc_sem_fit$data

ddsc_sem_GDI<-
  ddsc_sem(data=d_GDI,x = "gdi.z.cm",y1="means_y2",
           y2="means_y1",q_sesoi = .10)
round(ddsc_sem_GDI$results,3)
```

```
##                                     est      se       z pvalue  ci.lower ci.upper
## r_xy1_y2                          0.032   0.174   0.182  0.855    -0.309    0.373
## r_xy1                            -0.066   0.174  -0.380  0.704    -0.406    0.274
## r_xy2                            -0.079   0.174  -0.455  0.649    -0.419    0.261
## b_11                             -0.069   0.182  -0.380  0.704    -0.426    0.288
## b_21                             -0.075   0.165  -0.455  0.649    -0.398    0.248
## b_10                             -0.186   0.179  -1.035  0.300    -0.537    0.166
## b_20                              0.421   0.162   2.593  0.010     0.103    0.739
## res_cov_y1_y2                     0.948   0.235   4.037  0.000     0.488    1.409
## diff_b10_b20                     -0.606   0.032 -19.195  0.000    -0.668   -0.544
## diff_b11_b21                      0.006   0.032   0.182  0.855    -0.057    0.069
## diff_rxy1_rxy2                    0.013   0.027   0.478  0.633    -0.040    0.066
## q_b11_b21                         0.006   0.032   0.183  0.855    -0.057    0.069
## q_rxy1_rxy2                       0.013   0.027   0.478  0.633    -0.040    0.067
## cross_over_point                103.793 569.965   0.182  0.856 -1013.319 1220.905
## sum_b11_b21                      -0.144   0.346  -0.417  0.677    -0.822    0.534
## main_effect                      -0.072   0.173  -0.417  0.677    -0.411    0.267
## interaction_vs_main_effect       -0.066   0.192  -0.345  0.730    -0.443    0.310
## diff_abs_b11_abs_b21             -0.006   0.032  -0.182  0.855    -0.069    0.057
## abs_diff_b11_b21                  0.006   0.032   0.182  0.428    -0.057    0.069
## abs_sum_b11_b21                   0.144   0.346   0.417  0.338    -0.534    0.822
## dadas                            -0.138   0.364  -0.380  0.648    -0.852    0.575
## q_r_equivalence                  -0.087   0.027  -3.186  0.001        NA       NA
## q_b_equivalence                  -0.094   0.032  -2.926  0.002        NA       NA
## cross_over_point_equivalence    103.793 569.965   0.182  0.572        NA       NA
## cross_over_point_minimal_effect 103.793 569.965   0.182  0.428        NA       NA
```

``` r
round(ddsc_sem_GDI$variance_test,3)
```

```
##             est    se       z pvalue ci.lower ci.upper
## cov_y1y2  0.953 0.236   4.037  0.000    0.490    1.416
## var_y1    1.065 0.262   4.062  0.000    0.551    1.578
## var_y2    0.875 0.215   4.062  0.000    0.453    1.297
## var_diff  0.190 0.070   2.705  0.007    0.052    0.328
## var_ratio 1.217 0.066  18.411  0.000    1.088    1.347
## cor_y1y2  0.988 0.004 233.132  0.000    0.979    0.996
```

### Figure


``` r
# plot the results

# start with obtaining predicted values for means and differences


mod2_GDI_unstd<-lmer(pow.z~gndr.c+gdi.cm+gndr.c:gdi.cm+
                       (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
                     control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

mod2_GDI_unstd_red<-lmer(pow.z~gndr.c+
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


p1.pow.flags<-
  ggplot(p,aes(y=yvar,x=gdi.cm,color=gndr.c))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("Gender Development Index (Average 2002-2022)")+
  #ylim(c(min.y.pred,max.y.pred))+
  #xlim(c(0.60,1.00))+
  ylab("Value power (Average 2002-2022)")+
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

#p1.pow.flags


p2.pow.flags<-ggplot(p2,aes(y=yvar,x=gdi.cm))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("Gender Development Index (Average 2002-2022)")+
  #xlim(c(0.60,1.00))+
  #ylim(c(min.y.narrow,max.y.narrow))+
  ylab("Difference in value power")+
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

#p2.pow.flags


pflag_comb<-
  ggarrange(p1.pow.flags,p2.pow.flags,align = "v",
            ncol=1,nrow=2,heights=c(2,1))
```

```
## Warning in geom_text(inherit.aes = F, aes(x = 0.9, y = -0.45, label = coef_q, : All aesthetics have length 1, but the data has 282 rows.
## ℹ Please consider using `annotate()` or provide this layer with data containing a single row.
```

``` r
pflag_comb
```

![](Analysis_power_files/figure-html/unnamed-chunk-24-1.png)<!-- -->

``` r
png(filename = 
      "../results/pow/GDI_flags.png",
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
mod2_log_GDP<-lmer(pow.z~gndr.c+log_gdp.z.cm+
                     gndr.c:log_gdp.z.cm+
                 (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2_log_GDP)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: pow.z ~ gndr.c + log_gdp.z.cm + gndr.c:log_gdp.z.cm + (gndr.c |      cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1317419.2 1317507.2 -658701.6 1317403.2    441160 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.0773 -0.6583 -0.0130  0.6176  6.4331 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr 
##  cntry    (Intercept) 0.082833 0.28781       
##           gndr.c      0.003738 0.06114  -0.48
##  Residual             1.021807 1.01084       
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##                     Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)          0.03842    0.05015 33.00380   0.766 0.449001    
## gndr.c               0.22395    0.01120 33.85566  20.000  < 2e-16 ***
## log_gdp.z.cm        -0.21736    0.05033 33.02790  -4.319 0.000135 ***
## gndr.c:log_gdp.z.cm  0.02103    0.01137 35.26402   1.850 0.072753 .  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c lg_g..
## gndr.c      -0.457              
## lg_gdp.z.cm  0.022 -0.010       
## gndr.c:l_.. -0.010 -0.009 -0.451
```

``` r
getFE(mod2_log_GDP,round=3)
```

```
##                       Est.    SE     df      t     p     LL     UL
## (Intercept)          0.038 0.050 33.004  0.766 0.449 -0.064  0.140
## gndr.c               0.224 0.011 33.856 20.000 0.000  0.201  0.247
## log_gdp.z.cm        -0.217 0.050 33.028 -4.319 0.000 -0.320 -0.115
## gndr.c:log_gdp.z.cm  0.021 0.011 35.264  1.850 0.073 -0.002  0.044
```

``` r
getVC(mod2_log_GDP)
```

```
##        grp        var1   var2 sdcor  vcov
## 1    cntry (Intercept)   <NA>  0.29  0.08
## 2    cntry      gndr.c   <NA>  0.06  0.00
## 3    cntry (Intercept) gndr.c -0.48 -0.01
## 4 Residual        <NA>   <NA>  1.01  1.02
```

``` r
r2mlm(mod2_log_GDP,bargraph = F)
```

```
## $Decompositions
##                        total
## fixed           0.0437375965
## slope variation 0.0008028872
## mean variation  0.0721773155
## sigma2          0.8832822008
## 
## $R2s
##            total
## f   0.0437375965
## v   0.0008028872
## m   0.0721773155
## fv  0.0445404837
## fvm 0.1167177992
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
## Time difference of 28.81606 secs
```

``` r
round(ddsc_mod2_log_GDP$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.144        0.379        1.022     1.166 0.123   7204.152 0.999   0.999
## 2        0.5         0.118        0.343        1.022     1.140 0.103   6164.576 0.998   0.999
```

``` r
round(ddsc_mod2_log_GDP$descriptives,3)
```

```
##                          M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled log_gdp.z.cm
## means_y1             0.155 0.349    1.000           1.000    0.988           0.988       -0.599
## means_y1_scaled      0.421 0.950    1.000           1.000    0.988           0.988       -0.599
## means_y2            -0.068 0.385    0.988           0.988    1.000           1.000       -0.599
## means_y2_scaled     -0.186 1.048    0.988           0.988    1.000           1.000       -0.599
## log_gdp.z.cm        -0.022 1.012   -0.599          -0.599   -0.599          -0.599        1.000
## log_gdp.z.cm_scaled  0.000 1.000   -0.599          -0.599   -0.599          -0.599        1.000
## diff_score           0.223 0.068   -0.463          -0.463   -0.595          -0.595        0.319
## diff_score_scaled    0.606 0.184   -0.463          -0.463   -0.595          -0.595        0.319
##                     log_gdp.z.cm_scaled diff_score diff_score_scaled
## means_y1                         -0.599     -0.463            -0.463
## means_y1_scaled                  -0.599     -0.463            -0.463
## means_y2                         -0.599     -0.595            -0.595
## means_y2_scaled                  -0.599     -0.595            -0.595
## log_gdp.z.cm                      1.000      0.319             0.319
## log_gdp.z.cm_scaled               1.000      0.319             0.319
## diff_score                        0.319      1.000             1.000
## diff_score_scaled                 0.319      1.000             1.000
```

``` r
round(ddsc_mod2_log_GDP$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                      -0.310 0.168 35.264  -1.850   0.073   -0.651    0.030
## w_11                         -0.228 0.053 33.045  -4.288   0.000   -0.336   -0.120
## w_21                         -0.207 0.048 33.050  -4.307   0.000   -0.305   -0.109
## r_xy1                        -0.653 0.152 33.045  -4.288   0.000   -0.963   -0.343
## r_xy2                        -0.537 0.125 33.050  -4.307   0.000   -0.791   -0.283
## b_11                         -0.621 0.145 33.045  -4.288   0.000   -0.915   -0.326
## b_21                         -0.564 0.131 33.050  -4.307   0.000   -0.830   -0.297
## main_effect                  -0.217 0.050 33.028  -4.319   0.000   -0.320   -0.115
## moderator_effect              0.224 0.011 33.856  20.000   0.000    0.201    0.247
## interaction                   0.021 0.011 35.264   1.850   0.073   -0.002    0.044
## q_b11_b21                    -0.088    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                  -0.180    NA     NA      NA      NA       NA       NA
## cross_over_point            -10.648    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.196 0.046 33.112  -4.239   0.000   -0.291   -0.102
## interaction_vs_main_bscale   -0.535 0.126 33.112  -4.239   0.000   -0.792   -0.278
## interaction_vs_main_rscale   -0.479 0.114 33.126  -4.221   0.000   -0.710   -0.248
## dadas                        -0.414 0.096 33.050  -4.307   1.000   -0.609   -0.218
## dadas_bscale                 -1.127 0.262 33.050  -4.307   1.000   -1.659   -0.595
## dadas_rscale                 -1.074 0.249 33.050  -4.307   1.000   -1.582   -0.567
## abs_diff                      0.021 0.011 35.264   1.850   0.036   -0.002    0.044
## abs_sum                       0.435 0.101 33.028   4.319   0.000    0.230    0.639
## abs_diff_bscale               0.057 0.031 35.264   1.850   0.036   -0.006    0.120
## abs_sum_bscale                1.184 0.274 33.028   4.319   0.000    0.626    1.742
## abs_diff_rscale               0.116 0.039 34.438   2.964   0.003    0.036    0.195
## abs_sum_rscale                1.190 0.276 33.028   4.318   0.000    0.629    1.751
```

``` r
round(ddsc_mod2_log_GDP$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.013 -0.559 10.912  1.000  0.001
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
## r_xy1_y2                         -0.319 0.165  -1.932  0.053   -0.642    0.005
## r_xy1                            -0.599 0.139  -4.299  0.000   -0.872   -0.326
## r_xy2                            -0.599 0.139  -4.299  0.000   -0.872   -0.326
## b_11                             -0.628 0.146  -4.299  0.000   -0.914   -0.342
## b_21                             -0.569 0.132  -4.299  0.000   -0.828   -0.310
## b_10                             -0.186 0.144  -1.290  0.197   -0.467    0.096
## b_20                              0.421 0.130   3.229  0.001    0.165    0.676
## res_cov_y1_y2                     0.607 0.151   4.023  0.000    0.311    0.902
## diff_b10_b20                     -0.606 0.030 -20.241  0.000   -0.665   -0.548
## diff_b11_b21                     -0.059 0.030  -1.932  0.053   -0.118    0.001
## diff_rxy1_rxy2                    0.000 0.027   0.001  0.999   -0.053    0.053
## q_b11_b21                        -0.092 0.062  -1.477  0.140   -0.213    0.030
## q_rxy1_rxy2                       0.000 0.042   0.001  0.999   -0.083    0.083
## cross_over_point                -10.318 5.366  -1.923  0.054  -20.835    0.198
## sum_b11_b21                      -1.197 0.277  -4.319  0.000   -1.740   -0.654
## main_effect                      -0.598 0.139  -4.319  0.000   -0.870   -0.327
## interaction_vs_main_effect       -0.540 0.128  -4.225  0.000   -0.790   -0.289
## diff_abs_b11_abs_b21              0.059 0.030   1.932  0.053   -0.001    0.118
## abs_diff_b11_b21                  0.059 0.030   1.932  0.027   -0.001    0.118
## abs_sum_b11_b21                   1.197 0.277   4.319  0.000    0.654    1.740
## dadas                            -1.138 0.265  -4.299  1.000   -1.657   -0.619
## q_r_equivalence                  -0.100 0.042  -2.352  0.009       NA       NA
## q_b_equivalence                  -0.008 0.062  -0.134  0.447       NA       NA
## cross_over_point_equivalence     10.318 5.366   1.923  0.973       NA       NA
## cross_over_point_minimal_effect  10.318 5.366   1.923  0.027       NA       NA
```

``` r
round(ddsc_sem_log_GDP$variance_test,3)
```

```
##             est    se       z pvalue ci.lower ci.upper
## cov_y1y2  0.953 0.236   4.037  0.000    0.490    1.416
## var_y1    1.065 0.262   4.062  0.000    0.551    1.578
## var_y2    0.875 0.215   4.062  0.000    0.453    1.297
## var_diff  0.190 0.070   2.705  0.007    0.052    0.328
## var_ratio 1.217 0.066  18.411  0.000    1.088    1.347
## cor_y1y2  0.988 0.004 233.132  0.000    0.979    0.996
```

### Figure


``` r
# plot the results

# start with obtaining predicted values for means and differences


mod2_log_GDP_unstd<-lmer(pow.z~gndr.c+log_gdp.cm+
                           gndr.c:log_gdp.cm+
                       (gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
                     control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

mod2_log_GDP_unstd_red<-lmer(pow.z~gndr.c+
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


p1.pow.flags<-
  ggplot(p,aes(y=yvar,x=log_gdp.cm,color=gndr.c))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("log(GDP per capita) (Average 2002-2022)")+
  #ylim(c(min.y.pred,max.y.pred))+
  #xlim(c(0.60,1.00))+
  ylab("Value power (Average 2002-2022)")+
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

p2.pow.flags<-ggplot(p2,aes(y=yvar,x=log_gdp.cm))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.2)+
  xlab("log(GDP per capita) (Average 2002-2022)")+
  #xlim(c(0.60,1.00))+
  #ylim(c(min.y.narrow,max.y.narrow))+
  ylab("Difference in value power")+
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
  ggarrange(p1.pow.flags,p2.pow.flags,align = "v",
            ncol=1,nrow=2,heights=c(2,1))
```

```
## Warning in geom_text(inherit.aes = F, aes(x = 9.75, y = -0.45, label = coef_q, : All aesthetics have length 1, but the data has 3302 rows.
## ℹ Please consider using `annotate()` or provide this layer with data containing a single row.
```

``` r
pflag_comb
```

![](Analysis_power_files/figure-html/unnamed-chunk-27-1.png)<!-- -->

``` r
png(filename = 
      "../results/pow/log_GDP_flags.png",
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
mod3<-lmer(pow.z~gndr.c+essround.c+(gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod3)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: pow.z ~ gndr.c + essround.c + (gndr.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1317337.8 1317414.7 -658661.9 1317323.8    441161 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.0469 -0.6580 -0.0133  0.6179  6.3773 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr 
##  cntry    (Intercept) 0.129160 0.35939       
##           gndr.c      0.004193 0.06475  -0.57
##  Residual             1.021590 1.01074       
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##               Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)  4.388e-02  6.259e-02  3.300e+01   0.701    0.488    
## gndr.c       2.235e-01  1.179e-02  3.332e+01  18.952   <2e-16 ***
## essround.c  -5.640e-03  5.809e-04  4.411e+05  -9.708   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##            (Intr) gndr.c
## gndr.c     -0.539       
## essround.c -0.001 -0.001
```

``` r
getFE(mod3,round=3)
```

```
##               Est.    SE         df      t     p     LL     UL
## (Intercept)  0.044 0.063     32.995  0.701 0.488 -0.083  0.171
## gndr.c       0.224 0.012     33.322 18.952 0.000  0.200  0.248
## essround.c  -0.006 0.001 441084.988 -9.708 0.000 -0.007 -0.005
```

``` r
getVC(mod3)
```

```
##        grp        var1   var2 sdcor  vcov
## 1    cntry (Intercept)   <NA>  0.36  0.13
## 2    cntry      gndr.c   <NA>  0.06  0.00
## 3    cntry (Intercept) gndr.c -0.57 -0.01
## 4 Residual        <NA>   <NA>  1.01  1.02
```

``` r
r2mlm(mod3,bargraph = F)
```

```
## $Decompositions
##                        total
## fixed           0.0108453941
## slope variation 0.0008938902
## mean variation  0.1117057647
## sigma2          0.8765549511
## 
## $R2s
##            total
## f   0.0108453941
## v   0.0008938902
## m   0.1117057647
## fv  0.0117392842
## fvm 0.1234450489
```

``` r
anova(mod2,mod3)
```

```
## Data: diff_dat
## Models:
## mod2: pow.z ~ gndr.c + (gndr.c | cntry)
## mod3: pow.z ~ gndr.c + essround.c + (gndr.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod2    6 1317430 1317496 -658709   1317418                         
## mod3    7 1317338 1317415 -658662   1317324 94.238  1  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## mod4: random effect of time (Ess round)


``` r
mod4<-lmer(pow.z~gndr.c+essround.c+(gndr.c+essround.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
```

```
## Warning in checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv, : Model failed to converge
## with max|grad| = 0.00287418 (tol = 0.002, component 1)
```

``` r
summary(mod4)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: pow.z ~ gndr.c + essround.c + (gndr.c + essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1313168.8 1313278.8 -656574.4 1313148.8    441158 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.0019 -0.6567 -0.0074  0.6164  6.2484 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr       
##  cntry    (Intercept) 0.236723 0.48654             
##           gndr.c      0.004201 0.06482  -0.31      
##           essround.c  0.009850 0.09925  -0.26 -0.17
##  Residual             1.011427 1.00570             
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##               Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept) -0.0007933  0.0848480 30.6645327  -0.009    0.993    
## gndr.c       0.2238456  0.0118072 33.3707318  18.958   <2e-16 ***
## essround.c   0.0111781  0.0173351 27.3281011   0.645    0.524    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##            (Intr) gndr.c
## gndr.c     -0.300       
## essround.c -0.258 -0.161
## optimizer (bobyqa) convergence code: 0 (OK)
## Model failed to converge with max|grad| = 0.00287418 (tol = 0.002, component 1)
```

``` r
getFE(mod4,round=3)
```

```
##               Est.    SE     df      t     p     LL    UL
## (Intercept) -0.001 0.085 30.665 -0.009 0.993 -0.174 0.172
## gndr.c       0.224 0.012 33.371 18.958 0.000  0.200 0.248
## essround.c   0.011 0.017 27.328  0.645 0.524 -0.024 0.047
```

``` r
getVC(mod4)
```

```
##        grp        var1       var2 sdcor  vcov
## 1    cntry (Intercept)       <NA>  0.49  0.24
## 2    cntry      gndr.c       <NA>  0.06  0.00
## 3    cntry  essround.c       <NA>  0.10  0.01
## 4    cntry (Intercept)     gndr.c -0.31 -0.01
## 5    cntry (Intercept) essround.c -0.26 -0.01
## 6    cntry      gndr.c essround.c -0.17  0.00
## 7 Residual        <NA>       <NA>  1.01  1.01
```

``` r
r2mlm(mod4,bargraph = F)
```

```
## $Decompositions
##                      total
## fixed           0.01003441
## slope variation 0.05636877
## mean variation  0.17787607
## sigma2          0.75572076
## 
## $R2s
##          total
## f   0.01003441
## v   0.05636877
## m   0.17787607
## fv  0.06640317
## fvm 0.24427924
```

``` r
anova(mod2,mod3,mod4)
```

```
## Data: diff_dat
## Models:
## mod2: pow.z ~ gndr.c + (gndr.c | cntry)
## mod3: pow.z ~ gndr.c + essround.c + (gndr.c | cntry)
## mod4: pow.z ~ gndr.c + essround.c + (gndr.c + essround.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L)    Chisq Df Pr(>Chisq)    
## mod2    6 1317430 1317496 -658709   1317418                           
## mod3    7 1317338 1317415 -658662   1317324   94.238  1  < 2.2e-16 ***
## mod4   10 1313169 1313279 -656574   1313149 4174.916  3  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## mod5: fixed interaction between time and gender


``` r
mod5<-lmer(pow.z~gndr.c+essround.c+
             gndr.c:essround.c+(gndr.c+essround.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod5)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: pow.z ~ gndr.c + essround.c + gndr.c:essround.c + (gndr.c + essround.c |      cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1313098.1 1313219.1 -656538.1 1313076.1    441157 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -5.9732 -0.6562 -0.0095  0.6157  6.1967 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr       
##  cntry    (Intercept) 0.235967 0.48576             
##           gndr.c      0.004019 0.06339  -0.42      
##           essround.c  0.009818 0.09908  -0.25 -0.10
##  Residual             1.011269 1.00562             
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##                     Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)       -4.915e-04  8.471e-02  3.068e+01  -0.006    0.995    
## gndr.c             2.249e-01  1.157e-02  3.344e+01  19.442   <2e-16 ***
## essround.c         1.094e-02  1.731e-02  2.732e+01   0.632    0.532    
## gndr.c:essround.c -9.820e-03  1.148e-03  1.155e+05  -8.558   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c essrn.
## gndr.c      -0.400              
## essround.c  -0.255 -0.098       
## gndr.c:ssr.  0.000 -0.008  0.001
```

``` r
getFE(mod5,round=3)
```

```
##                     Est.    SE         df      t     p     LL     UL
## (Intercept)        0.000 0.085     30.685 -0.006 0.995 -0.173  0.172
## gndr.c             0.225 0.012     33.436 19.442 0.000  0.201  0.248
## essround.c         0.011 0.017     27.318  0.632 0.532 -0.025  0.046
## gndr.c:essround.c -0.010 0.001 115522.609 -8.558 0.000 -0.012 -0.008
```

``` r
getVC(mod5)
```

```
##        grp        var1       var2 sdcor  vcov
## 1    cntry (Intercept)       <NA>  0.49  0.24
## 2    cntry      gndr.c       <NA>  0.06  0.00
## 3    cntry  essround.c       <NA>  0.10  0.01
## 4    cntry (Intercept)     gndr.c -0.42 -0.01
## 5    cntry (Intercept) essround.c -0.25 -0.01
## 6    cntry      gndr.c essround.c -0.10  0.00
## 7 Residual        <NA>       <NA>  1.01  1.01
```

``` r
r2mlm(mod5,bargraph = F)
```

```
## $Decompositions
##                      total
## fixed           0.01029285
## slope variation 0.05618049
## mean variation  0.17756661
## sigma2          0.75596004
## 
## $R2s
##          total
## f   0.01029285
## v   0.05618049
## m   0.17756661
## fv  0.06647334
## fvm 0.24403996
```

``` r
anova(mod4,mod5)
```

```
## Data: diff_dat
## Models:
## mod4: pow.z ~ gndr.c + essround.c + (gndr.c + essround.c | cntry)
## mod5: pow.z ~ gndr.c + essround.c + gndr.c:essround.c + (gndr.c + essround.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod4   10 1313169 1313279 -656574   1313149                         
## mod5   11 1313098 1313219 -656538   1313076 72.698  1  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## mod6: random interaction between time and gender


``` r
mod6<-lmer(pow.z~gndr.c+essround.c+
             gndr.c:essround.c+(gndr.c+essround.c+gndr.c:essround.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
```

```
## Warning in checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv, : Model failed to converge
## with max|grad| = 0.00325505 (tol = 0.002, component 1)
```

``` r
summary(mod6)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: pow.z ~ gndr.c + essround.c + gndr.c:essround.c + (gndr.c + essround.c +  
##     gndr.c:essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1313065.5 1313230.5 -656517.8 1313035.5    441153 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -5.9804 -0.6571 -0.0103  0.6163  6.1712 
## 
## Random effects:
##  Groups   Name              Variance Std.Dev. Corr             
##  cntry    (Intercept)       0.238759 0.48863                   
##           gndr.c            0.004816 0.06940  -0.66            
##           essround.c        0.009984 0.09992  -0.26  0.29      
##           gndr.c:essround.c 0.000171 0.01308   0.47 -0.35 -0.63
##  Residual                   1.011109 1.00554                   
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##                    Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)       -0.001278   0.085211 30.666176  -0.015  0.98813    
## gndr.c             0.229829   0.012648 23.013520  18.171 3.82e-15 ***
## essround.c         0.011194   0.017452 27.394714   0.641  0.52658    
## gndr.c:essround.c -0.010785   0.002632 13.405742  -4.097  0.00119 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c essrn.
## gndr.c      -0.634              
## essround.c  -0.264  0.274       
## gndr.c:ssr.  0.409 -0.305 -0.549
## optimizer (bobyqa) convergence code: 0 (OK)
## Model failed to converge with max|grad| = 0.00325505 (tol = 0.002, component 1)
```

``` r
getFE(mod6,round=3)
```

```
##                     Est.    SE     df      t     p     LL     UL
## (Intercept)       -0.001 0.085 30.666 -0.015 0.988 -0.175  0.173
## gndr.c             0.230 0.013 23.014 18.171 0.000  0.204  0.256
## essround.c         0.011 0.017 27.395  0.641 0.527 -0.025  0.047
## gndr.c:essround.c -0.011 0.003 13.406 -4.097 0.001 -0.016 -0.005
```

``` r
getVC(mod6)
```

```
##         grp              var1              var2 sdcor  vcov
## 1     cntry       (Intercept)              <NA>  0.49  0.24
## 2     cntry            gndr.c              <NA>  0.07  0.00
## 3     cntry        essround.c              <NA>  0.10  0.01
## 4     cntry gndr.c:essround.c              <NA>  0.01  0.00
## 5     cntry       (Intercept)            gndr.c -0.66 -0.02
## 6     cntry       (Intercept)        essround.c -0.26 -0.01
## 7     cntry       (Intercept) gndr.c:essround.c  0.47  0.00
## 8     cntry            gndr.c        essround.c  0.29  0.00
## 9     cntry            gndr.c gndr.c:essround.c -0.35  0.00
## 10    cntry        essround.c gndr.c:essround.c -0.63  0.00
## 11 Residual              <NA>              <NA>  1.01  1.01
```

``` r
r2mlm(mod6,bargraph = F)
```

```
## $Decompositions
##                      total
## fixed           0.01072724
## slope variation 0.05759819
## mean variation  0.17941387
## sigma2          0.75226070
## 
## $R2s
##          total
## f   0.01072724
## v   0.05759819
## m   0.17941387
## fv  0.06832542
## fvm 0.24773930
```

``` r
anova(mod4,mod5,mod6)
```

```
## Data: diff_dat
## Models:
## mod4: pow.z ~ gndr.c + essround.c + (gndr.c + essround.c | cntry)
## mod5: pow.z ~ gndr.c + essround.c + gndr.c:essround.c + (gndr.c + essround.c | cntry)
## mod6: pow.z ~ gndr.c + essround.c + gndr.c:essround.c + (gndr.c + essround.c + gndr.c:essround.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod4   10 1313169 1313279 -656574   1313149                         
## mod5   11 1313098 1313219 -656538   1313076 72.698  1  < 2.2e-16 ***
## mod6   15 1313066 1313230 -656518   1313036 40.602  4  3.249e-08 ***
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
##         4.5 -0.0416 0.1020 Inf   -0.2419    0.1588  -0.406  0.6844
##        -4.5 -0.1908 0.1380 Inf   -0.4606    0.0789  -1.387  0.1656
## 
## gndr.c =  0.5:
##  essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##         4.5  0.1397 0.0971 Inf   -0.0507    0.3302   1.438  0.1503
##        -4.5  0.0875 0.1230 Inf   -0.1540    0.3290   0.710  0.4775
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
##  essround.c4.5 - (essround.c-4.5)   0.1493 0.164 Inf    -0.172     0.470   0.911  0.3623
## 
## gndr.c =  0.5:
##  contrast                         estimate    SE  df asymp.LCL asymp.UCL z.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   0.0522 0.151 Inf    -0.244     0.348   0.346  0.7293
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
##    -0.5        4.5 -0.0416 0.1020 Inf   -0.2419    0.1588  -0.406  0.6844
##     0.5        4.5  0.1397 0.0971 Inf   -0.0507    0.3302   1.438  0.1503
##    -0.5       -4.5 -0.1908 0.1380 Inf   -0.4606    0.0789  -1.387  0.1656
##     0.5       -4.5  0.0875 0.1230 Inf   -0.1540    0.3290   0.710  0.4775
## 
## Degrees-of-freedom method: asymptotic 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod6,adjust="none",infer=c(T,T))
```

```
##  contrast                                                 estimate     SE  df asymp.LCL asymp.UCL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.1813 0.0145 Inf    -0.210    -0.153
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   0.1493 0.1640 Inf    -0.172     0.470
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.1291 0.1540 Inf    -0.431     0.173
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.3306 0.1610 Inf     0.015     0.646
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)       0.0522 0.1510 Inf    -0.244     0.348
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.2784 0.0198 Inf    -0.317    -0.240
##  z.ratio p.value
##  -12.541  <.0001
##    0.911  0.3623
##   -0.838  0.4022
##    2.053  0.0400
##    0.346  0.7293
##  -14.067  <.0001
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
##  diff_ESS10    0.181 0.0145 Inf     0.153     0.210  12.541  <.0001
##  diff_ESS1     0.278 0.0198 Inf     0.240     0.317  14.067  <.0001
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
##  diff_ESS10 - diff_ESS1  -0.0971 0.0237 Inf    -0.143   -0.0506  -4.097  <.0001
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
      obs_mean_wt=weighted.mean(x=pow.z,w=pspwght),
      obs_sd_wt=sqrt(wtd.var(pow.z,pspwght)),
      obs_mean=mean(pow.z),
      obs_sd=sd(pow.z),
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
  ylab("Mean-level of value power")+
  theme(legend.title=element_blank(),
        legend.text  = element_text(size = 16),
        axis.title.x = element_text(size = 18, face = "bold"),   
        axis.title.y = element_text(size = 18, face = "bold"),   
        axis.text.x  = element_text(size = 14),                  
        axis.text.y  = element_text(size = 14))

p_time_trends
```

![](Analysis_power_files/figure-html/unnamed-chunk-33-1.png)<!-- -->

``` r
png(filename = 
      "../results/pow/time_trends.png",
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
pred_cntry_dat$pow.z_mean<-predict(mod6,newdata=pred_cntry_dat)

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

range(pred_cntry_dat$pow.z_mean)
```

```
## [1] -0.8646921  1.1195709
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
      obs_mean_wt=weighted.mean(x=pow.z,w=pspwght),
      obs_sd_wt=sqrt(wtd.var(pow.z,pspwght)),
      obs_mean=mean(pow.z),
      obs_sd=sd(pow.z),
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

pdf("../results/pow/country_specific_time_trends.pdf", width = 8, height = 6)

for (ctry in countries) {
  print(
    ggplot(pred_cntry_dat[pred_cntry_dat$cntry == ctry, ], 
           aes(x = year, y = pow.z_mean, color = gender)) +
      geom_smooth(method = "lm",formula = "y ~ x", se = FALSE) +
      geom_point(size=8) +
      geom_point(aes(x=year,y=obs_mean_wt),size=8,shape = 1,alpha=.50)+
      geom_errorbar(aes(ymin=obs_mean_wt_LL, ymax=obs_mean_wt_UL),alpha = .50)+
      scale_color_manual(values = my_colors) +
      geom_flag(aes(country=tolower(ctry)))+
      scale_y_continuous(limits = c(-1.1, 1.1)) +
      scale_x_continuous(limits = c(2001, 2021), breaks = seq(2002, 2020, 2)) +
      ggtitle(paste("Country:", ctry))+
    ylab("Mean-level of value power")+
      theme(legend.title=element_blank())
  )
}
```

```
## Warning: Removed 1 row containing non-finite outside the scale range (`stat_smooth()`).
```

```
## Warning: Removed 1 row containing missing values or values outside the scale range (`geom_point()`).
## Removed 1 row containing missing values or values outside the scale range (`geom_point()`).
```

```
## Warning: Removed 1 row containing missing values or values outside the scale range (`geom_flag()`).
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
         aes(x = year, y = pow.z_mean, color = gender)) +
  geom_smooth(method = "lm", se = FALSE, formula="y~x") +
  geom_point(size=3.5) +
  geom_point(aes(x=year,y=obs_mean_wt),size=3.5,shape = 1,alpha=.70)+
  geom_errorbar(aes(ymin=obs_mean_wt_LL, ymax=obs_mean_wt_UL),alpha = .50)+
  scale_color_manual(values = my_colors) +
  #geom_flag(aes(country=tolower(ctry)))+
  scale_y_continuous(limits = c(-1.1, 1.1)) +
  scale_x_continuous(limits = c(2001, 2021), breaks = seq(2002, 2020, 2)) +
  #ggtitle(paste("Country:", ctry))+
  ylab("Mean-level of value power")+
  xlab("Year")+
  theme(legend.title=element_blank(),legend.position = "top",
        axis.text.x = element_text(angle = 45,size = 6,hjust=1))+
  facet_wrap(~CLDR,nrow=6,ncol=6)+
  #facet_wrap(~cntry,nrow=6,ncol=6)+
  geom_flag(aes(country=tolower(cntry)),size=2)

facet_plot
```

```
## Warning: Removed 1 row containing non-finite outside the scale range (`stat_smooth()`).
```

```
## Warning: Removed 1 row containing missing values or values outside the scale range (`geom_point()`).
## Removed 1 row containing missing values or values outside the scale range (`geom_point()`).
```

```
## Warning: Removed 1 row containing missing values or values outside the scale range (`geom_flag()`).
```

![](Analysis_power_files/figure-html/unnamed-chunk-35-1.png)<!-- -->

``` r
png(filename = 
      "../results/pow/country_time_trend_facets.png",
    units = "cm",
    width = 21.0,height=29.7*(3/4),res = 600)
facet_plot
```

```
## Warning: Removed 1 row containing non-finite outside the scale range (`stat_smooth()`).
```

```
## Warning: Removed 1 row containing missing values or values outside the scale range (`geom_point()`).
## Removed 1 row containing missing values or values outside the scale range (`geom_point()`).
```

```
## Warning: Removed 1 row containing missing values or values outside the scale range (`geom_flag()`).
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
## [1] 58.40003
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
## 1    0.26                0.11                   -0.17                    0.03                      0.19
## 2    0.19               -0.01                    0.00                   -0.01                     -0.01
## 3    0.26                0.77                   -0.08                    0.73                      0.81
## 4    0.24                0.20                   -0.01                    0.20                      0.20
## 5    0.17               -0.57                   -0.05                   -0.60                     -0.55
## 6    0.28                0.44                   -0.13                    0.38                      0.51
## 7    0.30               -0.30                   -0.04                   -0.32                     -0.28
## 8    0.29                0.09                   -0.09                    0.04                      0.14
## 9    0.29               -0.24                   -0.04                   -0.26                     -0.22
## 10   0.20               -0.69                   -0.09                   -0.73                     -0.65
## 11   0.38               -0.19                   -0.04                   -0.21                     -0.17
## 12   0.24               -0.20                   -0.13                   -0.26                     -0.14
## 13   0.27               -0.18                   -0.12                   -0.24                     -0.11
## 14   0.16               -0.33                    0.01                   -0.33                     -0.34
## 15   0.19               -0.43                   -0.06                   -0.46                     -0.40
## 16   0.18                0.34                   -0.16                    0.26                      0.42
## 17   0.16                0.01                    0.00                    0.01                      0.01
## 18   0.23                0.02                    0.01                    0.03                      0.01
## 19   0.29               -0.22                   -0.13                   -0.28                     -0.16
## 20   0.15                0.95                   -0.20                    0.85                      1.05
## 21   0.18               -0.52                    0.09                   -0.48                     -0.56
## 22   0.20               -1.73                    0.00                   -1.73                     -1.73
## 23   0.40                3.99                   -0.53                    3.73                      4.25
## 24   0.32               -0.08                   -0.22                   -0.19                      0.03
## 25   0.21               -0.05                   -0.24                   -0.17                      0.07
## 26   0.26                0.08                   -0.11                    0.02                      0.14
## 27   0.18               -0.58                   -0.16                   -0.67                     -0.50
## 28   0.11               -0.18                   -0.06                   -0.21                     -0.15
## 29   0.25                0.07                   -0.08                    0.03                      0.11
## 30   0.22                0.02                   -0.02                    0.01                      0.04
## 31   0.26                0.20                   -0.11                    0.15                      0.26
## 32   0.15                1.90                   -0.13                    1.84                      1.96
## 33   0.14                0.62                   -0.11                    0.56                      0.67
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
## 1     ES               -0.69
## 2     PT               -0.58
## 3     CY               -0.57
## 4     LT               -0.52
## 5     GR               -0.33
## 6     DE               -0.30
## 7     EE               -0.24
## 8     IS               -0.22
## 9     FR               -0.20
## 10    FI               -0.19
## 11    GB               -0.18
## 12    RU               -0.18
## 13    NL               -0.08
## 14    NO               -0.05
## 15    BE               -0.01
## 16    IE                0.01
## 17    IL                0.02
## 18    SI                0.02
## 19    SE                0.07
## 20    PL                0.08
## 21    DK                0.09
## 22    AT                0.11
## 23    CH                0.20
## 24    SK                0.20
## 25    HU                0.34
## 26    CZ                0.44
## 27    UA                0.62
## 28    BG                0.77
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
## 1     NO                   -0.24
## 2     NL                   -0.22
## 3     AT                   -0.17
## 4     HU                   -0.16
## 5     PT                   -0.16
## 6     CZ                   -0.13
## 7     FR                   -0.13
## 8     IS                   -0.13
## 9     GB                   -0.12
## 10    PL                   -0.11
## 11    SK                   -0.11
## 12    UA                   -0.11
## 13    DK                   -0.09
## 14    ES                   -0.09
## 15    BG                   -0.08
## 16    SE                   -0.08
## 17    RU                   -0.06
## 18    CY                   -0.05
## 19    DE                   -0.04
## 20    EE                   -0.04
## 21    FI                   -0.04
## 22    SI                   -0.02
## 23    CH                   -0.01
## 24    BE                    0.00
## 25    IE                    0.00
## 26    GR                    0.01
## 27    IL                    0.01
## 28    LT                    0.09
```


#### Country-level predictors: GEI


``` r
mod6_GEI<-lmer(pow.z~gndr.c+essround.c+
             gndr.c:essround.c+
               gei.z.cm:gndr.c+gei.z.cm:essround.c+gei.z.cm:gndr.c:essround.c+
               (gndr.c+essround.c+gndr.c:essround.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_GEI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: pow.z ~ gndr.c + essround.c + gndr.c:essround.c + gei.z.cm:gndr.c +  
##     gei.z.cm:essround.c + gei.z.cm:gndr.c:essround.c + (gndr.c +  
##     essround.c + gndr.c:essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1280854.6 1281052.2 -640409.3 1280818.6    431760 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -5.9926 -0.6590 -0.0103  0.6180  6.1859 
## 
## Random effects:
##  Groups   Name              Variance  Std.Dev. Corr             
##  cntry    (Intercept)       0.2506452 0.50064                   
##           gndr.c            0.0048325 0.06952  -0.66            
##           essround.c        0.0120698 0.10986  -0.65  0.57      
##           gndr.c:essround.c 0.0002271 0.01507   0.71 -0.53 -0.71
##  Residual                   1.0058262 1.00291                   
## Number of obs: 431778, groups:  cntry, 32
## 
## Fixed effects:
##                             Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)                -0.019478   0.088662 29.749764  -0.220 0.827614    
## gndr.c                      0.233255   0.012877 12.803666  18.115 1.67e-10 ***
## essround.c                  0.010288   0.019477 21.966584   0.528 0.602652    
## gndr.c:essround.c          -0.011063   0.003003  9.619676  -3.684 0.004501 ** 
## gndr.c:gei.z.cm             0.001244   0.010633 36.332021   0.117 0.907530    
## essround.c:gei.z.cm        -0.065748   0.015193 32.077865  -4.328 0.000138 ***
## gndr.c:essround.c:gei.z.cm  0.006886   0.002689 30.752543   2.561 0.015576 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c essrn. gnd.:. gn.:.. es.:..
## gndr.c      -0.632                                   
## essround.c  -0.649  0.547                            
## gndr.c:ssr.  0.634 -0.462 -0.635                     
## gndr.c:g.z.  0.000 -0.042 -0.001 -0.044              
## essrnd.c:.. -0.001  0.000 -0.005  0.000  0.224       
## gndr.c:.:..  0.000 -0.029  0.001 -0.132  0.072 -0.335
```

``` r
getFE(mod6_GEI,round=2,p.round = 5)
```

```
##                             Est.   SE    df     t       p    LL    UL
## (Intercept)                -0.02 0.09 29.75 -0.22 0.82761 -0.20  0.16
## gndr.c                      0.23 0.01 12.80 18.11 0.00000  0.21  0.26
## essround.c                  0.01 0.02 21.97  0.53 0.60265 -0.03  0.05
## gndr.c:essround.c          -0.01 0.00  9.62 -3.68 0.00450 -0.02  0.00
## gndr.c:gei.z.cm             0.00 0.01 36.33  0.12 0.90753 -0.02  0.02
## essround.c:gei.z.cm        -0.07 0.02 32.08 -4.33 0.00014 -0.10 -0.03
## gndr.c:essround.c:gei.z.cm  0.01 0.00 30.75  2.56 0.01558  0.00  0.01
```

``` r
getVC(mod6_GEI)
```

```
##         grp              var1              var2 sdcor  vcov
## 1     cntry       (Intercept)              <NA>  0.50  0.25
## 2     cntry            gndr.c              <NA>  0.07  0.00
## 3     cntry        essround.c              <NA>  0.11  0.01
## 4     cntry gndr.c:essround.c              <NA>  0.02  0.00
## 5     cntry       (Intercept)            gndr.c -0.66 -0.02
## 6     cntry       (Intercept)        essround.c -0.65 -0.04
## 7     cntry       (Intercept) gndr.c:essround.c  0.71  0.01
## 8     cntry            gndr.c        essround.c  0.57  0.00
## 9     cntry            gndr.c gndr.c:essround.c -0.53  0.00
## 10    cntry        essround.c gndr.c:essround.c -0.71  0.00
## 11 Residual              <NA>              <NA>  1.00  1.01
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_GEI,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] -20.88574
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_GEI,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] -32.84154
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
##         4.5  0.3227 0.1010 Inf     0.125    0.5207   3.194  0.0014
##        -4.5 -0.3616 0.1740 Inf    -0.703   -0.0203  -2.076  0.0379
## 
## gei.z.cm =  0:
##  essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##         4.5  0.0268 0.0739 Inf    -0.118    0.1716   0.363  0.7167
##        -4.5 -0.0658 0.1600 Inf    -0.380    0.2480  -0.411  0.6812
## 
## gei.z.cm =  1:
##  essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##         4.5 -0.2691 0.1000 Inf    -0.466   -0.0724  -2.682  0.0073
##        -4.5  0.2301 0.1740 Inf    -0.111    0.5710   1.323  0.1859
## 
## Degrees-of-freedom method: asymptotic 
## Confidence level used: 0.95
```

``` r
pairs(change_mod6_GEI,adjust="none",infer=c(T,T))
```

```
## gei.z.cm = -1:
##  contrast                         estimate    SE  df asymp.LCL asymp.UCL z.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   0.6843 0.223 Inf     0.248    1.1210   3.071  0.0021
## 
## gei.z.cm =  0:
##  contrast                         estimate    SE  df asymp.LCL asymp.UCL z.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   0.0926 0.175 Inf    -0.251    0.4361   0.528  0.5974
## 
## gei.z.cm =  1:
##  contrast                         estimate    SE  df asymp.LCL asymp.UCL z.ratio p.value
##  essround.c4.5 - (essround.c-4.5)  -0.4991 0.222 Inf    -0.934   -0.0644  -2.250  0.0244
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
##  gndr.c essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##    -0.5        4.5  0.2471 0.1030 Inf    0.0457    0.4484   2.405  0.0162
##     0.5        4.5  0.3983 0.1010 Inf    0.2010    0.5956   3.957  0.0001
##    -0.5       -4.5 -0.5180 0.1840 Inf   -0.8785   -0.1576  -2.817  0.0049
##     0.5       -4.5 -0.2052 0.1650 Inf   -0.5288    0.1183  -1.244  0.2137
## 
## gei.z.cm =  0:
##  gndr.c essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##    -0.5        4.5 -0.0649 0.0748 Inf   -0.2116    0.0818  -0.867  0.3858
##     0.5        4.5  0.1186 0.0736 Inf   -0.0256    0.2627   1.612  0.1070
##    -0.5       -4.5 -0.2073 0.1690 Inf   -0.5388    0.1242  -1.226  0.2203
##     0.5       -4.5  0.0757 0.1510 Inf   -0.2209    0.3724   0.501  0.6167
## 
## gei.z.cm =  1:
##  gndr.c essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##    -0.5        4.5 -0.3769 0.1020 Inf   -0.5765   -0.1773  -3.701  0.0002
##     0.5        4.5 -0.1612 0.0997 Inf   -0.3567    0.0343  -1.616  0.1061
##    -0.5       -4.5  0.1034 0.1840 Inf   -0.2565    0.4634   0.563  0.5732
##     0.5       -4.5  0.3567 0.1650 Inf    0.0338    0.6797   2.165  0.0304
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
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.1512 0.0233 Inf   -0.1969   -0.1056
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   0.7651 0.2330 Inf    0.3091    1.2211
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)    0.4523 0.2160 Inf    0.0283    0.8764
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.9163 0.2300 Inf    0.4649    1.3678
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)       0.6036 0.2140 Inf    0.1836    1.0235
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.3128 0.0280 Inf   -0.3676   -0.2579
##  z.ratio p.value
##   -6.499  <.0001
##    3.288  0.0010
##    2.091  0.0366
##    3.978  0.0001
##    2.817  0.0048
##  -11.178  <.0001
## 
## gei.z.cm =  0:
##  contrast                                                 estimate     SE  df asymp.LCL asymp.UCL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.1835 0.0137 Inf   -0.2103   -0.1566
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   0.1424 0.1840 Inf   -0.2186    0.5033
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.1407 0.1690 Inf   -0.4711    0.1898
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.3258 0.1830 Inf   -0.0321    0.6838
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)       0.0428 0.1670 Inf   -0.2846    0.3702
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.2830 0.0226 Inf   -0.3273   -0.2388
##  z.ratio p.value
##  -13.398  <.0001
##    0.773  0.4395
##   -0.834  0.4041
##    1.784  0.0744
##    0.256  0.7978
##  -12.542  <.0001
## 
## gei.z.cm =  1:
##  contrast                                                 estimate     SE  df asymp.LCL asymp.UCL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.2157 0.0197 Inf   -0.2544   -0.1770
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  -0.4803 0.2320 Inf   -0.9341   -0.0266
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.7336 0.2150 Inf   -1.1556   -0.3117
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)     -0.2646 0.2290 Inf   -0.7141    0.1848
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      -0.5179 0.2130 Inf   -0.9355   -0.1004
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.2533 0.0268 Inf   -0.3058   -0.2008
##  z.ratio p.value
##  -10.931  <.0001
##   -2.075  0.0380
##   -3.408  0.0007
##   -1.154  0.2485
##   -2.431  0.0150
##   -9.457  <.0001
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
##  diff_ESS10    0.151 0.0233 Inf     0.106     0.197   6.499  <.0001
##  diff_ESS1     0.313 0.0280 Inf     0.258     0.368  11.178  <.0001
## 
## gei.z.cm =  0:
##  contrast   estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10    0.183 0.0137 Inf     0.157     0.210  13.398  <.0001
##  diff_ESS1     0.283 0.0226 Inf     0.239     0.327  12.542  <.0001
## 
## gei.z.cm =  1:
##  contrast   estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10    0.216 0.0197 Inf     0.177     0.254  10.931  <.0001
##  diff_ESS1     0.253 0.0268 Inf     0.201     0.306   9.457  <.0001
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
##  diff_ESS10 - diff_ESS1  -0.1615 0.0386 Inf    -0.237   -0.0859  -4.187  <.0001
## 
## gei.z.cm =  0:
##  contrast               estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10 - diff_ESS1  -0.0996 0.0270 Inf    -0.153   -0.0466  -3.684  0.0002
## 
## gei.z.cm =  1:
##  contrast               estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10 - diff_ESS1  -0.0376 0.0338 Inf    -0.104    0.0287  -1.112  0.2662
## 
## Degrees-of-freedom method: asymptotic 
## Confidence level used: 0.95
```


#### Country-level predictors: GGGI


``` r
mod6_GGGI<-lmer(pow.z~gndr.c+essround.c+
             gndr.c:essround.c+
               gggi.z.cm:gndr.c+gggi.z.cm:essround.c+gggi.z.cm:gndr.c:essround.c+
               (gndr.c+essround.c+gndr.c:essround.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_GGGI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: pow.z ~ gndr.c + essround.c + gndr.c:essround.c + gggi.z.cm:gndr.c +  
##     gggi.z.cm:essround.c + gggi.z.cm:gndr.c:essround.c + (gndr.c +  
##     essround.c + gndr.c:essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##  938300.4  938492.3 -469132.2  938264.4    314628 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -5.9658 -0.6569 -0.0039  0.6189  6.0692 
## 
## Random effects:
##  Groups   Name              Variance  Std.Dev. Corr             
##  cntry    (Intercept)       0.1544820 0.39304                   
##           gndr.c            0.0032475 0.05699  -0.42            
##           essround.c        0.0021940 0.04684  -0.22  0.00      
##           gndr.c:essround.c 0.0001041 0.01020   0.29 -0.10 -0.38
##  Residual                   1.0152619 1.00760                   
## Number of obs: 314646, groups:  cntry, 33
## 
## Fixed effects:
##                              Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)                  0.077402   0.068795 32.538720   1.125   0.2688    
## gndr.c                       0.222074   0.010859 25.781972  20.450   <2e-16 ***
## essround.c                  -0.009100   0.008598 25.250063  -1.058   0.2999    
## gndr.c:essround.c           -0.007119   0.002624 22.061986  -2.713   0.0127 *  
## gndr.c:gggi.z.cm             0.012480   0.010628 36.442244   1.174   0.2479    
## essround.c:gggi.z.cm        -0.013704   0.009078 30.461256  -1.510   0.1414    
## gndr.c:essround.c:gggi.z.cm  0.005294   0.002900 26.453171   1.826   0.0792 .  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c essrn. gnd.:. gn.:.. es.:..
## gndr.c      -0.381                                   
## essround.c  -0.218  0.004                            
## gndr.c:ssr.  0.196 -0.113 -0.258                     
## gndr.c:gg..  0.000 -0.021  0.001 -0.026              
## essrnd.c:.. -0.004  0.003 -0.081  0.017 -0.081       
## gndr.c:.:..  0.002 -0.027  0.017 -0.093 -0.050 -0.222
```

``` r
getFE(mod6_GGGI,round=2,p.round = 5)
```

```
##                              Est.   SE    df     t       p    LL   UL
## (Intercept)                  0.08 0.07 32.54  1.13 0.26877 -0.06 0.22
## gndr.c                       0.22 0.01 25.78 20.45 0.00000  0.20 0.24
## essround.c                  -0.01 0.01 25.25 -1.06 0.29985 -0.03 0.01
## gndr.c:essround.c           -0.01 0.00 22.06 -2.71 0.01268 -0.01 0.00
## gndr.c:gggi.z.cm             0.01 0.01 36.44  1.17 0.24788 -0.01 0.03
## essround.c:gggi.z.cm        -0.01 0.01 30.46 -1.51 0.14145 -0.03 0.00
## gndr.c:essround.c:gggi.z.cm  0.01 0.00 26.45  1.83 0.07917  0.00 0.01
```

``` r
getVC(mod6_GGGI)
```

```
##         grp              var1              var2 sdcor  vcov
## 1     cntry       (Intercept)              <NA>  0.39  0.15
## 2     cntry            gndr.c              <NA>  0.06  0.00
## 3     cntry        essround.c              <NA>  0.05  0.00
## 4     cntry gndr.c:essround.c              <NA>  0.01  0.00
## 5     cntry       (Intercept)            gndr.c -0.42 -0.01
## 6     cntry       (Intercept)        essround.c -0.22  0.00
## 7     cntry       (Intercept) gndr.c:essround.c  0.29  0.00
## 8     cntry            gndr.c        essround.c  0.00  0.00
## 9     cntry            gndr.c gndr.c:essround.c -0.10  0.00
## 10    cntry        essround.c gndr.c:essround.c -0.38  0.00
## 11 Residual              <NA>              <NA>  1.01  1.02
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_GGGI,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] 78.02624
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_GGGI,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] 39.1326
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
##         4.5  0.0981 0.0838 30.3  -0.0729    0.269   1.171  0.2506
##        -4.5  0.0567 0.0964 30.3  -0.1401    0.253   0.588  0.5608
## 
## gggi.z.cm =  0:
##  essround.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5  0.0364 0.0712 19.1  -0.1125    0.185   0.512  0.6145
##        -4.5  0.1184 0.0860 21.9  -0.0600    0.297   1.377  0.1825
## 
## gggi.z.cm =  1:
##  essround.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##         4.5 -0.0252 0.0804 27.9  -0.1899    0.139  -0.314  0.7560
##        -4.5  0.1800 0.0939 28.6  -0.0122    0.372   1.916  0.0654
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
##  essround.c4.5 - (essround.c-4.5)   0.0414 0.1170 30.0   -0.197   0.2803   0.354  0.7257
## 
## gggi.z.cm =  0:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)  -0.0819 0.0774 25.2   -0.241   0.0774  -1.058  0.2999
## 
## gggi.z.cm =  1:
##  contrast                         estimate     SE   df lower.CL upper.CL t.ratio p.value
##  essround.c4.5 - (essround.c-4.5)  -0.2052 0.1080 28.4   -0.426   0.0156  -1.902  0.0673
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
##    -0.5        4.5  0.0212 0.0870 30.6 -0.15627   0.1988   0.244  0.8087
##     0.5        4.5  0.1750 0.0821 30.2  0.00739   0.3426   2.132  0.0413
##    -0.5       -4.5 -0.0760 0.1010 30.6 -0.28115   0.1291  -0.757  0.4551
##     0.5       -4.5  0.1894 0.0937 29.9 -0.00199   0.3808   2.021  0.0523
## 
## gggi.z.cm =  0:
##  gndr.c essround.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5        4.5 -0.0586 0.0733 18.9 -0.21193   0.0948  -0.799  0.4339
##     0.5        4.5  0.1315 0.0699 19.1 -0.01475   0.2777   1.881  0.0753
##    -0.5       -4.5 -0.0087 0.0896 21.9 -0.19453   0.1771  -0.097  0.9235
##     0.5       -4.5  0.2454 0.0831 21.4  0.07290   0.4179   2.955  0.0075
## 
## gggi.z.cm =  1:
##  gndr.c essround.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5        4.5 -0.1384 0.0833 28.2 -0.30907   0.0323  -1.660  0.1079
##     0.5        4.5  0.0880 0.0787 27.7 -0.07336   0.2493   1.117  0.2734
##    -0.5       -4.5  0.0586 0.0980 28.9 -0.14181   0.2591   0.598  0.5542
##     0.5       -4.5  0.3014 0.0913 28.3  0.11448   0.4883   3.302  0.0026
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
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.1537 0.0234 30.8  -0.2014 -0.10609  -6.583
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   0.0973 0.1230 29.9  -0.1531  0.34773   0.794
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.1682 0.1190 30.0  -0.4103  0.07400  -1.418
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.2510 0.1170 30.3   0.0114  0.49068   2.138
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      -0.0144 0.1140 30.3  -0.2473  0.21845  -0.127
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.2655 0.0245 32.2  -0.3154 -0.21547 -10.815
##  p.value
##   <.0001
##   0.4337
##   0.1664
##   0.0407
##   0.9002
##   <.0001
## 
## gggi.z.cm =  0:
##  contrast                                                 estimate     SE   df lower.CL upper.CL t.ratio
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.1900 0.0151 28.6  -0.2210 -0.15911 -12.574
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  -0.0499 0.0812 24.8  -0.2172  0.11750  -0.614
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.3040 0.0781 26.4  -0.4644 -0.14357  -3.893
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.1402 0.0782 24.3  -0.0211  0.30143   1.793
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      -0.1139 0.0752 25.9  -0.2686  0.04068  -1.515
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.2541 0.0169 22.4  -0.2892 -0.21906 -15.020
##  p.value
##   <.0001
##   0.5449
##   0.0006
##   0.0855
##   0.1419
##   <.0001
## 
## gggi.z.cm =  1:
##  contrast                                                 estimate     SE   df lower.CL upper.CL t.ratio
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.2263 0.0212 29.6  -0.2697 -0.18299 -10.667
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  -0.1970 0.1130 28.3  -0.4286  0.03459  -1.742
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.4398 0.1100 28.4  -0.6640 -0.21560  -4.016
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.0293 0.1080 28.7  -0.1923  0.25099   0.271
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      -0.2134 0.1050 28.6  -0.4285  0.00159  -2.031
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.2428 0.0237 27.6  -0.2914 -0.19409 -10.222
##  p.value
##   <.0001
##   0.0924
##   0.0004
##   0.7886
##   0.0516
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
##  diff_ESS10    0.154 0.0234 30.8    0.106    0.201   6.583  <.0001
##  diff_ESS1     0.265 0.0245 32.2    0.215    0.315  10.815  <.0001
## 
## gggi.z.cm =  0:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10    0.190 0.0151 28.6    0.159    0.221  12.574  <.0001
##  diff_ESS1     0.254 0.0169 22.4    0.219    0.289  15.020  <.0001
## 
## gggi.z.cm =  1:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10    0.226 0.0212 29.6    0.183    0.270  10.667  <.0001
##  diff_ESS1     0.243 0.0237 27.6    0.194    0.291  10.222  <.0001
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
##  diff_ESS10 - diff_ESS1  -0.1117 0.0368 26.7  -0.1872  -0.0362  -3.038  0.0053
## 
## gggi.z.cm =  0:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10 - diff_ESS1  -0.0641 0.0236 22.1  -0.1130  -0.0151  -2.713  0.0127
## 
## gggi.z.cm =  1:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS10 - diff_ESS1  -0.0164 0.0335 23.3  -0.0857   0.0529  -0.490  0.6289
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

#### Country-level predictors: GDI


``` r
mod6_GDI<-lmer(pow.z~gndr.c+essround.c+
             gndr.c:essround.c+
               gdi.z.cm:gndr.c+gdi.z.cm:essround.c+gdi.z.cm:gndr.c:essround.c+
               (gndr.c+essround.c+gndr.c:essround.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_GDI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: pow.z ~ gndr.c + essround.c + gndr.c:essround.c + gdi.z.cm:gndr.c +  
##     gdi.z.cm:essround.c + gdi.z.cm:gndr.c:essround.c + (gndr.c +  
##     essround.c + gndr.c:essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1313063.0 1313260.9 -656513.5 1313027.0    441150 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -5.9773 -0.6571 -0.0108  0.6163  6.1710 
## 
## Random effects:
##  Groups   Name              Variance Std.Dev. Corr             
##  cntry    (Intercept)       0.236919 0.48674                   
##           gndr.c            0.004579 0.06767  -0.66            
##           essround.c        0.007915 0.08897  -0.29  0.23      
##           gndr.c:essround.c 0.000153 0.01237   0.48 -0.31 -0.57
##  Residual                   1.011109 1.00554                   
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##                              Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)                -0.0003551  0.0848801 30.6762747  -0.004 0.996689    
## gndr.c                      0.2294883  0.0123634 24.3371128  18.562 7.11e-16 ***
## essround.c                  0.0110730  0.0155513 26.1562644   0.712 0.482748    
## gndr.c:essround.c          -0.0106399  0.0025266 15.2110537  -4.211 0.000734 ***
## gndr.c:gdi.z.cm            -0.0101237  0.0103682 39.9604352  -0.976 0.334730    
## essround.c:gdi.z.cm        -0.0455638  0.0151781 30.0296340  -3.002 0.005361 ** 
## gndr.c:essround.c:gdi.z.cm  0.0031484  0.0027671 32.5180740   1.138 0.263510    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c essrn. gnd.:. gn.:.. es.:..
## gndr.c      -0.626                                   
## essround.c  -0.294  0.225                            
## gndr.c:ssr.  0.415 -0.267 -0.491                     
## gndr.c:gd..  0.000 -0.021 -0.001 -0.017              
## essrnd.c:..  0.001 -0.001 -0.006  0.004  0.054       
## gndr.c:.:..  0.000 -0.019  0.003 -0.002  0.110 -0.358
```

``` r
getFE(mod6_GDI,round=2,p.round = 5)
```

```
##                             Est.   SE    df     t       p    LL    UL
## (Intercept)                 0.00 0.08 30.68  0.00 0.99669 -0.17  0.17
## gndr.c                      0.23 0.01 24.34 18.56 0.00000  0.20  0.25
## essround.c                  0.01 0.02 26.16  0.71 0.48275 -0.02  0.04
## gndr.c:essround.c          -0.01 0.00 15.21 -4.21 0.00073 -0.02 -0.01
## gndr.c:gdi.z.cm            -0.01 0.01 39.96 -0.98 0.33473 -0.03  0.01
## essround.c:gdi.z.cm        -0.05 0.02 30.03 -3.00 0.00536 -0.08 -0.01
## gndr.c:essround.c:gdi.z.cm  0.00 0.00 32.52  1.14 0.26351  0.00  0.01
```

``` r
getVC(mod6_GDI)
```

```
##         grp              var1              var2 sdcor  vcov
## 1     cntry       (Intercept)              <NA>  0.49  0.24
## 2     cntry            gndr.c              <NA>  0.07  0.00
## 3     cntry        essround.c              <NA>  0.09  0.01
## 4     cntry gndr.c:essround.c              <NA>  0.01  0.00
## 5     cntry       (Intercept)            gndr.c -0.66 -0.02
## 6     cntry       (Intercept)        essround.c -0.29 -0.01
## 7     cntry       (Intercept) gndr.c:essround.c  0.48  0.00
## 8     cntry            gndr.c        essround.c  0.23  0.00
## 9     cntry            gndr.c gndr.c:essround.c -0.31  0.00
## 10    cntry        essround.c gndr.c:essround.c -0.57  0.00
## 11 Residual              <NA>              <NA>  1.01  1.01
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_GDI,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] 20.72382
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_GDI,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] 10.51446
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
##  essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##         4.5  0.2545 0.1150 Inf    0.0283    0.4807   2.205  0.0274
##        -4.5 -0.2552 0.1430 Inf   -0.5347    0.0242  -1.790  0.0735
## 
## gdi.z.cm =  0:
##  essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##         4.5  0.0495 0.0928 Inf   -0.1324    0.2313   0.533  0.5939
##        -4.5 -0.0502 0.1250 Inf   -0.2949    0.1946  -0.402  0.6878
## 
## gdi.z.cm =  1:
##  essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##         4.5 -0.1556 0.1150 Inf   -0.3810    0.0699  -1.352  0.1762
##        -4.5  0.1549 0.1420 Inf   -0.1236    0.4333   1.090  0.2758
## 
## Degrees-of-freedom method: asymptotic 
## Confidence level used: 0.95
```

``` r
pairs(change_mod6_GDI,adjust="none",infer=c(T,T))
```

```
## gdi.z.cm = -1:
##  contrast                         estimate    SE  df asymp.LCL asymp.UCL z.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   0.5097 0.196 Inf     0.125    0.8942   2.598  0.0094
## 
## gdi.z.cm =  0:
##  contrast                         estimate    SE  df asymp.LCL asymp.UCL z.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   0.0997 0.140 Inf    -0.175    0.3740   0.712  0.4764
## 
## gdi.z.cm =  1:
##  contrast                         estimate    SE  df asymp.LCL asymp.UCL z.ratio p.value
##  essround.c4.5 - (essround.c-4.5)  -0.3104 0.195 Inf    -0.693    0.0717  -1.592  0.1114
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
##    -0.5        4.5  0.1657 0.1190 Inf   -0.0676    0.3991   1.392  0.1639
##     0.5        4.5  0.3433 0.1130 Inf    0.1222    0.5644   3.043  0.0023
##    -0.5       -4.5 -0.4061 0.1500 Inf   -0.6998   -0.1123  -2.709  0.0068
##     0.5       -4.5 -0.1044 0.1360 Inf   -0.3709    0.1621  -0.768  0.4426
## 
## gdi.z.cm =  0:
##  gndr.c essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##    -0.5        4.5 -0.0413 0.0955 Inf   -0.2285    0.1458  -0.433  0.6651
##     0.5        4.5  0.1403 0.0906 Inf   -0.0373    0.3179   1.548  0.1216
##    -0.5       -4.5 -0.1889 0.1320 Inf   -0.4468    0.0691  -1.435  0.1513
##     0.5       -4.5  0.0885 0.1180 Inf   -0.1437    0.3207   0.747  0.4551
## 
## gdi.z.cm =  1:
##  gndr.c essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##    -0.5        4.5 -0.2484 0.1190 Inf   -0.4808   -0.0159  -2.094  0.0362
##     0.5        4.5 -0.0627 0.1120 Inf   -0.2831    0.1576  -0.558  0.5768
##    -0.5       -4.5  0.0283 0.1490 Inf   -0.2644    0.3211   0.190  0.8497
##     0.5       -4.5  0.2814 0.1350 Inf    0.0158    0.5470   2.077  0.0378
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
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.1776 0.0227 Inf   -0.2220   -0.1331
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   0.5718 0.2040 Inf    0.1722    0.9714
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)    0.2701 0.1940 Inf   -0.1110    0.6512
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.7493 0.1990 Inf    0.3589    1.1398
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)       0.4477 0.1900 Inf    0.0760    0.8194
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.3017 0.0242 Inf   -0.3492   -0.2542
##  z.ratio p.value
##   -7.830  <.0001
##    2.804  0.0050
##    1.389  0.1648
##    3.762  0.0002
##    2.361  0.0182
##  -12.447  <.0001
## 
## gdi.z.cm =  0:
##  contrast                                                 estimate     SE  df asymp.LCL asymp.UCL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.1816 0.0144 Inf   -0.2098   -0.1534
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   0.1475 0.1460 Inf   -0.1384    0.4335
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.1298 0.1380 Inf   -0.3997    0.1401
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.3291 0.1430 Inf    0.0484    0.6099
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)       0.0518 0.1350 Inf   -0.2123    0.3159
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.2774 0.0189 Inf   -0.3144   -0.2403
##  z.ratio p.value
##  -12.618  <.0001
##    1.011  0.3118
##   -0.943  0.3458
##    2.298  0.0216
##    0.384  0.7008
##  -14.678  <.0001
## 
## gdi.z.cm =  1:
##  contrast                                                 estimate     SE  df asymp.LCL asymp.UCL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.1857 0.0220 Inf   -0.2287   -0.1426
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  -0.2767 0.2030 Inf   -0.6739    0.1204
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.5298 0.1930 Inf   -0.9086   -0.1510
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)     -0.0911 0.1980 Inf   -0.4791    0.2970
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      -0.3441 0.1890 Inf   -0.7136    0.0254
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.2531 0.0244 Inf   -0.3009   -0.2053
##  z.ratio p.value
##   -8.456  <.0001
##   -1.366  0.1721
##   -2.741  0.0061
##   -0.460  0.6456
##   -1.825  0.0680
##  -10.374  <.0001
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
##  diff_ESS10    0.178 0.0227 Inf     0.133     0.222   7.830  <.0001
##  diff_ESS1     0.302 0.0242 Inf     0.254     0.349  12.447  <.0001
## 
## gdi.z.cm =  0:
##  contrast   estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10    0.182 0.0144 Inf     0.153     0.210  12.618  <.0001
##  diff_ESS1     0.277 0.0189 Inf     0.240     0.314  14.678  <.0001
## 
## gdi.z.cm =  1:
##  contrast   estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10    0.186 0.0220 Inf     0.143     0.229   8.456  <.0001
##  diff_ESS1     0.253 0.0244 Inf     0.205     0.301  10.374  <.0001
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
##  diff_ESS10 - diff_ESS1  -0.1241 0.0338 Inf    -0.190   -0.0579  -3.675  0.0002
## 
## gdi.z.cm =  0:
##  contrast               estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10 - diff_ESS1  -0.0958 0.0227 Inf    -0.140   -0.0512  -4.211  <.0001
## 
## gdi.z.cm =  1:
##  contrast               estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10 - diff_ESS1  -0.0674 0.0337 Inf    -0.133   -0.0014  -2.002  0.0453
## 
## Degrees-of-freedom method: asymptotic 
## Confidence level used: 0.95
```

#### Country-level predictors: log_GDP


``` r
mod6_log_GDP<-lmer(pow.z~gndr.c+essround.c+
             gndr.c:essround.c+
               log_gdp.z.cm:gndr.c+log_gdp.z.cm:essround.c+log_gdp.z.cm:gndr.c:essround.c+
               (gndr.c+essround.c+gndr.c:essround.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
```

```
## Warning in checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv, : Model failed to converge
## with max|grad| = 0.00366423 (tol = 0.002, component 1)
```

``` r
summary(mod6_log_GDP)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: pow.z ~ gndr.c + essround.c + gndr.c:essround.c + log_gdp.z.cm:gndr.c +  
##     log_gdp.z.cm:essround.c + log_gdp.z.cm:gndr.c:essround.c +  
##     (gndr.c + essround.c + gndr.c:essround.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1313062.9 1313260.8 -656513.4 1313026.9    441150 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -5.9803 -0.6570 -0.0103  0.6162  6.1706 
## 
## Random effects:
##  Groups   Name              Variance  Std.Dev. Corr             
##  cntry    (Intercept)       0.2451565 0.49513                   
##           gndr.c            0.0050775 0.07126  -0.68            
##           essround.c        0.0093011 0.09644  -0.46  0.42      
##           gndr.c:essround.c 0.0001837 0.01355   0.59 -0.43 -0.66
##  Residual                   1.0111047 1.00554                   
## Number of obs: 441168, groups:  cntry, 33
## 
## Fixed effects:
##                                 Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)                    -0.003451   0.086343 30.710437  -0.040  0.96838    
## gndr.c                          0.230186   0.012963 17.802941  17.758 9.07e-13 ***
## essround.c                      0.010586   0.016850 27.115429   0.628  0.53508    
## gndr.c:essround.c              -0.010853   0.002727 12.433775  -3.980  0.00171 ** 
## gndr.c:log_gdp.z.cm            -0.003696   0.010081 34.869023  -0.367  0.71612    
## essround.c:log_gdp.z.cm        -0.049526   0.015060 31.777647  -3.289  0.00246 ** 
## gndr.c:essround.c:log_gdp.z.cm  0.004414   0.002475 24.618957   1.784  0.08683 .  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c essrn. gnd.:. g.:_.. e.:_..
## gndr.c      -0.653                                   
## essround.c  -0.462  0.400                            
## gndr.c:ssr.  0.512 -0.370 -0.568                     
## gndr.c:l_..  0.000 -0.027  0.002 -0.002              
## essrnd.:_..  0.004 -0.001  0.012 -0.005  0.145       
## gndr.:.:_.. -0.002  0.001 -0.006 -0.114 -0.033 -0.415
## optimizer (bobyqa) convergence code: 0 (OK)
## Model failed to converge with max|grad| = 0.00366423 (tol = 0.002, component 1)
```

``` r
getFE(mod6_log_GDP,round=2,p.round = 5)
```

```
##                                 Est.   SE    df     t       p    LL    UL
## (Intercept)                     0.00 0.09 30.71 -0.04 0.96838 -0.18  0.17
## gndr.c                          0.23 0.01 17.80 17.76 0.00000  0.20  0.26
## essround.c                      0.01 0.02 27.12  0.63 0.53508 -0.02  0.05
## gndr.c:essround.c              -0.01 0.00 12.43 -3.98 0.00171 -0.02  0.00
## gndr.c:log_gdp.z.cm             0.00 0.01 34.87 -0.37 0.71612 -0.02  0.02
## essround.c:log_gdp.z.cm        -0.05 0.02 31.78 -3.29 0.00246 -0.08 -0.02
## gndr.c:essround.c:log_gdp.z.cm  0.00 0.00 24.62  1.78 0.08683  0.00  0.01
```

``` r
getVC(mod6_log_GDP)
```

```
##         grp              var1              var2 sdcor  vcov
## 1     cntry       (Intercept)              <NA>  0.50  0.25
## 2     cntry            gndr.c              <NA>  0.07  0.01
## 3     cntry        essround.c              <NA>  0.10  0.01
## 4     cntry gndr.c:essround.c              <NA>  0.01  0.00
## 5     cntry       (Intercept)            gndr.c -0.68 -0.02
## 6     cntry       (Intercept)        essround.c -0.46 -0.02
## 7     cntry       (Intercept) gndr.c:essround.c  0.59  0.00
## 8     cntry            gndr.c        essround.c  0.42  0.00
## 9     cntry            gndr.c gndr.c:essround.c -0.43  0.00
## 10    cntry        essround.c gndr.c:essround.c -0.66  0.00
## 11 Residual              <NA>              <NA>  1.01  1.01
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_log_GDP,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] 6.844359
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_log_GDP,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] -7.419821
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
##         4.5  0.2671 0.1080 Inf     0.056    0.4781   2.481  0.0131
##        -4.5 -0.2740 0.1540 Inf    -0.576    0.0282  -1.777  0.0755
## 
## log_gdp.z.cm =  0:
##  essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##         4.5  0.0442 0.0846 Inf    -0.122    0.2100   0.522  0.6015
##        -4.5 -0.0511 0.1390 Inf    -0.323    0.2208  -0.368  0.7127
## 
## log_gdp.z.cm =  1:
##  essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##         4.5 -0.1787 0.1090 Inf    -0.393    0.0353  -1.637  0.1017
##        -4.5  0.1718 0.1550 Inf    -0.131    0.4749   1.111  0.2667
## 
## Degrees-of-freedom method: asymptotic 
## Confidence level used: 0.95
```

``` r
pairs(change_mod6_log_GDP,adjust="none",infer=c(T,T))
```

```
## log_gdp.z.cm = -1:
##  contrast                         estimate    SE  df asymp.LCL asymp.UCL z.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   0.5410 0.202 Inf     0.145    0.9373   2.676  0.0075
## 
## log_gdp.z.cm =  0:
##  contrast                         estimate    SE  df asymp.LCL asymp.UCL z.ratio p.value
##  essround.c4.5 - (essround.c-4.5)   0.0953 0.152 Inf    -0.202    0.3925   0.628  0.5298
## 
## log_gdp.z.cm =  1:
##  contrast                         estimate    SE  df asymp.LCL asymp.UCL z.ratio p.value
##  essround.c4.5 - (essround.c-4.5)  -0.3505 0.205 Inf    -0.751    0.0505  -1.713  0.0867
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
##    -0.5        4.5  0.1845 0.1110 Inf   -0.0323    0.4012   1.668  0.0954
##     0.5        4.5  0.3496 0.1060 Inf    0.1424    0.5568   3.307  0.0009
##    -0.5       -4.5 -0.4252 0.1630 Inf   -0.7445   -0.1060  -2.610  0.0090
##     0.5       -4.5 -0.1227 0.1460 Inf   -0.4090    0.1636  -0.840  0.4011
## 
## log_gdp.z.cm =  0:
##  gndr.c essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##    -0.5        4.5 -0.0465 0.0868 Inf   -0.2166    0.1237  -0.536  0.5923
##     0.5        4.5  0.1349 0.0830 Inf   -0.0278    0.2975   1.625  0.1041
##    -0.5       -4.5 -0.1906 0.1470 Inf   -0.4783    0.0971  -1.298  0.1942
##     0.5       -4.5  0.0884 0.1310 Inf   -0.1683    0.3452   0.675  0.4997
## 
## log_gdp.z.cm =  1:
##  gndr.c essround.c  emmean     SE  df asymp.LCL asymp.UCL z.ratio p.value
##    -0.5        4.5 -0.2774 0.1120 Inf   -0.4971   -0.0578  -2.476  0.0133
##     0.5        4.5 -0.0799 0.1070 Inf   -0.2898    0.1299  -0.746  0.4554
##    -0.5       -4.5  0.0440 0.1630 Inf   -0.2762    0.3643   0.270  0.7875
##     0.5       -4.5  0.2995 0.1460 Inf    0.0124    0.5866   2.044  0.0409
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
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.1652 0.0214 Inf   -0.2071  -0.12327
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   0.6097 0.2110 Inf    0.1962   1.02326
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)    0.3071 0.1980 Inf   -0.0809   0.69515
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.7749 0.2080 Inf    0.3679   1.18192
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)       0.4723 0.1950 Inf    0.0910   0.85366
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.3026 0.0266 Inf   -0.3547  -0.25046
##  z.ratio p.value
##   -7.723  <.0001
##    2.890  0.0039
##    1.551  0.1208
##    3.731  0.0002
##    2.427  0.0152
##  -11.379  <.0001
## 
## log_gdp.z.cm =  0:
##  contrast                                                 estimate     SE  df asymp.LCL asymp.UCL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.1813 0.0142 Inf   -0.2091  -0.15358
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)   0.1441 0.1590 Inf   -0.1674   0.45564
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.1349 0.1470 Inf   -0.4229   0.15308
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)      0.3255 0.1570 Inf    0.0172   0.63374
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)       0.0464 0.1450 Inf   -0.2378   0.33069
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.2790 0.0209 Inf   -0.3200  -0.23808
##  z.ratio p.value
##  -12.799  <.0001
##    0.907  0.3646
##   -0.918  0.3585
##    2.069  0.0385
##    0.320  0.7488
##  -13.356  <.0001
## 
## log_gdp.z.cm =  1:
##  contrast                                                 estimate     SE  df asymp.LCL asymp.UCL
##  (gndr.c-0.5 essround.c4.5) - gndr.c0.5 essround.c4.5      -0.1975 0.0195 Inf   -0.2358  -0.15928
##  (gndr.c-0.5 essround.c4.5) - (gndr.c-0.5 essround.c-4.5)  -0.3215 0.2130 Inf   -0.7396   0.09659
##  (gndr.c-0.5 essround.c4.5) - (gndr.c0.5 essround.c-4.5)   -0.5769 0.2000 Inf   -0.9696  -0.18426
##  gndr.c0.5 essround.c4.5 - (gndr.c-0.5 essround.c-4.5)     -0.1240 0.2100 Inf   -0.5355   0.28756
##  gndr.c0.5 essround.c4.5 - (gndr.c0.5 essround.c-4.5)      -0.3794 0.1970 Inf   -0.7650   0.00609
##  (gndr.c-0.5 essround.c-4.5) - (gndr.c0.5 essround.c-4.5)  -0.2555 0.0251 Inf   -0.3047  -0.20620
##  z.ratio p.value
##  -10.124  <.0001
##   -1.507  0.1318
##   -2.880  0.0040
##   -0.590  0.5549
##   -1.929  0.0537
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

diff_mod6_log_GDP<-contrast(change_in_diff_mod6_log_GDP,method = changes_in_diff,adjust="none",infer=c(T,T))

# gender differences at ESS1 and ESS10
diff_mod6_log_GDP
```

```
## log_gdp.z.cm = -1:
##  contrast   estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10    0.165 0.0214 Inf     0.123     0.207   7.723  <.0001
##  diff_ESS1     0.303 0.0266 Inf     0.250     0.355  11.379  <.0001
## 
## log_gdp.z.cm =  0:
##  contrast   estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10    0.181 0.0142 Inf     0.154     0.209  12.799  <.0001
##  diff_ESS1     0.279 0.0209 Inf     0.238     0.320  13.356  <.0001
## 
## log_gdp.z.cm =  1:
##  contrast   estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10    0.198 0.0195 Inf     0.159     0.236  10.124  <.0001
##  diff_ESS1     0.255 0.0251 Inf     0.206     0.305  10.163  <.0001
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
##  diff_ESS10 - diff_ESS1  -0.1374 0.0350 Inf    -0.206  -0.06886  -3.929  0.0001
## 
## log_gdp.z.cm =  0:
##  contrast               estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10 - diff_ESS1  -0.0977 0.0245 Inf    -0.146  -0.04957  -3.980  0.0001
## 
## log_gdp.z.cm =  1:
##  contrast               estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  diff_ESS10 - diff_ESS1  -0.0579 0.0312 Inf    -0.119   0.00322  -1.857  0.0633
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

