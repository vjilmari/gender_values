---
title: "Analysis for universalism values"
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


# Data for the analysis

* use the testing partition of the data set (diff_dat)
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
##   11 2314 1577 2218 1368  667    0 2381    0 1282 1833 1524 1745 1529 2745 1548 2117 1985  893  825 2783
##     
##        LT   LV   ME   NL   NO   PL   PT   RS   RU   SE   SI   SK   TR   UA
##   1     0    0    0 2337 1819 2065 1482    0    0 1682 1488    0    0    0
##   2     0    0    0 1858 1575 1683 2024    0    0 1678 1384 1425 1790 1896
##   3     0    0    0 1860 1550 1685 2182    0 2339 1604 1465 1711    0 1885
##   4     0 1970    0 1724 1391 1596 2337    0 2446 1556 1257 1789 2305 1766
##   5  1632    0    0 1801 1530 1719 2139    0 2557 1463 1369 1803    0 1779
##   6  2108    0    0 1828 1610 1866 2138    0 2429 1838 1244 1827    0 2064
##   7  2241    0    0 1823 1423 1594 1242    0    0 1761 1189    0    0    0
##   8  2079    0    0 1669 1530 1675 1254    0 2374 1526 1295    0    0    0
##   9  1677  891 1188 1657 1396 1443 1045 1969    0 1510 1307 1061    0    0
##   10 1606    0 1248 1466 1408    0 1827    0    0    0 1232 1395    0    0
##   11 1337 1225 1581 1677 1318 1423 1366 1511    0 1204 1227 1414    0 2589
```

``` r
# range of sample sizes
range(table(diff_dat$cntry))
```

```
## [1]  3480 27753
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
## `summarise()` has regrouped the output.
## ℹ Summaries were computed grouped by cntry and essround.
## ℹ Output is grouped by cntry.
## ℹ Use `summarise(.groups = "drop_last")` to silence this message.
## ℹ Use `summarise(.by = c(cntry, essround))` for per-operation grouping (`?dplyr::dplyr_by`) instead.
```

``` r
grand_mean_uni<-mean(cntry.uni$uni.cm)
grand_sd_uni<-mean(cntry.uni$uni.csd)

# standardized
diff_dat$uni.z<-(diff_dat$uni-grand_mean_uni)/grand_sd_uni
hist(diff_dat$uni.z)
```

![](Analysis_universalism_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

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
  pivot_longer(
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
                   uni.z.wt=weighted.mean(x=uni.z,w=pspwght),
                   uni.z=mean(uni.z),
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

# universalism

cntry_uni_frame<-
  diff_dat %>% group_by(cntry,essround) %>%
  summarise('uni M' = weighted.mean(x=uni.z,w=pspwght),
            'uni SD' = sqrt(wtd.var(uni.z,pspwght))) %>%
  group_by(cntry) %>%
  summarise('uni M' = mean(x=`uni M`),
            'uni SD'= mean(x=`uni SD`))
```

```
## `summarise()` has regrouped the output.
## ℹ Summaries were computed grouped by cntry and essround.
## ℹ Output is grouped by cntry.
## ℹ Use `summarise(.groups = "drop_last")` to silence this message.
## ℹ Use `summarise(.by = c(cntry, essround))` for per-operation grouping (`?dplyr::dplyr_by`) instead.
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
## `summarise()` has regrouped the output.
## ℹ Summaries were computed grouped by cntry and essround.
## ℹ Output is grouped by cntry.
## ℹ Use `summarise(.groups = "drop_last")` to silence this message.
## ℹ Use `summarise(.by = c(cntry, essround))` for per-operation grouping (`?dplyr::dplyr_by`) instead.
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
## `summarise()` has regrouped the output.
## ℹ Summaries were computed grouped by cntry and essround.
## ℹ Output is grouped by cntry.
## ℹ Use `summarise(.groups = "drop_last")` to silence this message.
## ℹ Use `summarise(.by = c(cntry, essround))` for per-operation grouping (`?dplyr::dplyr_by`) instead.
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
## # A tibble: 34 × 10
##    cntry `n ESS rounds`     n `uni M` `uni SD` `uni M Women` `uni SD Women` `uni M Men` `uni SD Men`
##    <chr>          <dbl> <dbl>   <dbl>    <dbl>         <dbl>          <dbl>       <dbl>        <dbl>
##  1 AT                 7 15400  0.133     1.02        0.247            0.982      0.0107        1.04 
##  2 BE                11 18886  0.125     0.863       0.197            0.842      0.0487        0.879
##  3 BG                 7 14857 -0.143     1.09       -0.0639           1.07      -0.227         1.12 
##  4 CH                11 18087  0.320     0.859       0.415            0.828      0.221         0.880
##  5 CY                 6  5771  0.351     0.871       0.368            0.849      0.332         0.893
##  6 CZ                 9 18934 -0.398     1.07       -0.274            1.05      -0.533         1.08 
##  7 DE                10 27753  0.128     0.925       0.232            0.886      0.0188        0.951
##  8 DK                 8 12198 -0.0247    1.02        0.0704           0.982     -0.123         1.04 
##  9 EE                10 17974 -0.117     0.932      -0.00429          0.897     -0.252         0.954
## 10 ES                10 18785  0.381     0.904       0.425            0.884      0.335         0.922
## # ℹ 24 more rows
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
  dplyr::select(
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
print(cntry_desc_tbl,n=35)
```

```
## # A tibble: 34 × 14
##    Country     `n ESS rounds`     n `uni M` `uni SD` `uni M Women` `uni SD Women` `uni M Men` `uni SD Men`
##    <chr>                <dbl> <dbl> <chr>   <chr>    <chr>         <chr>          <chr>       <chr>       
##  1 Austria                  7 15400 0.13    1.02     0.25          0.98           0.01        1.04        
##  2 Belgium                 11 18886 0.12    0.86     0.20          0.84           0.05        0.88        
##  3 Bulgaria                 7 14857 -0.14   1.09     -0.06         1.07           -0.23       1.12        
##  4 Switzerland             11 18087 0.32    0.86     0.42          0.83           0.22        0.88        
##  5 Cyprus                   6  5771 0.35    0.87     0.37          0.85           0.33        0.89        
##  6 Czechia                  9 18934 -0.40   1.07     -0.27         1.05           -0.53       1.08        
##  7 Germany                 10 27753 0.13    0.92     0.23          0.89           0.02        0.95        
##  8 Denmark                  8 12198 -0.02   1.02     0.07          0.98           -0.12       1.04        
##  9 Estonia                 10 17974 -0.12   0.93     -0.00         0.90           -0.25       0.95        
## 10 Spain                   10 18785 0.38    0.90     0.43          0.88           0.33        0.92        
## 11 Finland                 11 19568 0.18    0.95     0.37          0.88           -0.02       0.97        
## 12 France                  11 20457 0.11    1.07     0.16          1.07           0.05        1.07        
## 13 UK                      11 22979 0.02    1.02     0.09          1.00           -0.04       1.03        
## 14 Greece                   6 15212 0.20    0.93     0.19          0.92           0.21        0.95        
## 15 Croatia                  5  7914 0.02    1.06     0.12          1.02           -0.09       1.08        
## 16 Hungary                 11 18123 0.00    1.01     0.07          0.99           -0.08       1.03        
## 17 Ireland                 11 22562 0.05    1.05     0.13          1.04           -0.03       1.06        
## 18 Israel                   7 14857 0.06    1.06     0.12          1.06           0.01        1.06        
## 19 Iceland                  6  4654 0.11    0.97     0.23          0.92           -0.01       1.00        
## 20 Italy                    5 11441 0.08    0.99     0.13          0.98           0.02        0.99        
## 21 Lithuania                7 13059 -0.69   1.15     -0.60         1.14           -0.80       1.16        
## 22 Latvia                   3  4088 -0.33   1.08     -0.18         1.05           -0.53       1.09        
## 23 Montenegro               3  4028 -0.50   1.25     -0.39         1.22           -0.61       1.26        
## 24 Netherlands             11 19722 -0.01   0.87     0.08          0.85           -0.10       0.89        
## 25 Norway                  11 16505 -0.19   0.98     -0.11         0.96           -0.29       1.00        
## 26 Poland                  10 16737 0.07    0.91     0.15          0.88           -0.01       0.94        
## 27 Portugal                11 19070 -0.24   1.03     -0.22         1.03           -0.26       1.04        
## 28 Serbia                   2  3499 0.17    1.06     0.31          0.98           0.02        1.11        
## 29 Russia                   5 12139 -0.18   1.10     -0.12         1.08           -0.25       1.11        
## 30 Sweden                  10 16104 -0.01   0.97     0.13          0.93           -0.17       0.98        
## 31 Slovenia                11 14463 0.30    0.83     0.38          0.81           0.22        0.84        
## 32 Slovakia                 8 12547 -0.15   0.95     -0.08         0.94           -0.23       0.95        
## 33 Turkey                   2  4108 0.18    0.96     0.16          0.97           0.20        0.95        
## 34 Ukraine                  6 12054 -0.28   1.23     -0.22         1.22           -0.37       1.25        
## # ℹ 5 more variables: D <chr>, GEI <chr>, GGGI <chr>, GDI <chr>, GDP <chr>
```

``` r
export(cntry_desc_tbl,"../results/uni/cntry_desc_tbl.xlsx",overwrite=T)
```

## Country-level correlations table


``` r
cor_frame<-
  desc_frame %>%
  dplyr::select(
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
  dplyr::select(-GDP)

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
##   1. VBMT       -0.01 0.24                                                                            
##                                                                                                       
##   2. VBMT_Women 0.07  0.24 .98                                                                        
##                            [.97, .99]                                                                 
##                                                                                                       
##   3. VBMT_Men   -0.10 0.26 .98          .94                                                           
##                            [.97, .99]   [.88, .97]                                                    
##                                                                                                       
##   4. D          -0.17 0.09 .26          .08          .43                                              
##                            [-.08, .55]  [-.26, .41]  [.10, .67]                                       
##                                                                                                       
##   5. GEI        0.87  0.07 .15          .24          .08          -.38                                
##                            [-.20, .47]  [-.12, .54]  [-.27, .41]  [-.64, -.05]                        
##                                                                                                       
##   6. GGGI       0.74  0.05 .06          .18          -.05         -.58         .73                    
##                            [-.28, .39]  [-.17, .49]  [-.38, .30]  [-.77, -.31] [.52, .86]             
##                                                                                                       
##   7. GDI        0.98  0.03 -.48         -.43         -.53         -.39         .07         .19        
##                            [-.71, -.18] [-.67, -.11] [-.74, -.24] [-.64, -.06] [-.28, .41] [-.16, .50]
##                                                                                                       
##   8. log_GDP    10.61 0.41 .34          .38          .30          -.11         .72         .62        
##                            [.00, .61]   [.04, .63]   [-.04, .58]  [-.43, .24]  [.50, .85]  [.36, .79] 
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
## 1468575.0 1468608.3 -734284.5 1468569.0    492340 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -8.3350 -0.5462  0.1198  0.6582  5.0330 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  cntry    (Intercept) 0.06003  0.245   
##  Residual             1.01593  1.008   
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##             Estimate Std. Error       df t value Pr(>|t|)
## (Intercept) -0.01425    0.04205 33.97221  -0.339    0.737
```

``` r
getVC(mod0)
```

```
##        grp        var1 var2 sdcor vcov
## 1    cntry (Intercept) <NA>  0.25 0.06
## 2 Residual        <NA> <NA>  1.01 1.02
```

``` r
r2mlm(mod0,bargraph = F)
```

```
## $Decompositions
##                    total within between
## fixed, within   0.000000      0      NA
## fixed, between  0.000000     NA       0
## slope variation 0.000000      0      NA
## mean variation  0.055791     NA       1
## sigma2          0.944209      1      NA
## 
## $R2s
##        total within between
## f1  0.000000      0      NA
## f2  0.000000     NA       0
## v   0.000000      0      NA
## m   0.055791     NA       1
## f   0.000000     NA      NA
## fv  0.000000      0      NA
## fvm 0.055791     NA      NA
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
## 1465124.1 1465168.6 -732558.1 1465116.1    492339 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -8.4907 -0.5498  0.0972  0.6540  4.8819 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  cntry    (Intercept) 0.06068  0.2463  
##  Residual             1.00883  1.0044  
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##               Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept) -1.744e-02  4.228e-02  3.397e+01  -0.413    0.683    
## gndr.c      -1.683e-01  2.860e-03  4.923e+05 -58.864   <2e-16 ***
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
## (Intercept) -0.017 0.042     33.972  -0.413 0.683 -0.103  0.068
## gndr.c      -0.168 0.003 492310.109 -58.864 0.000 -0.174 -0.163
```

``` r
getVC(mod1)
```

```
##        grp        var1 var2 sdcor vcov
## 1    cntry (Intercept) <NA>  0.25 0.06
## 2 Residual        <NA> <NA>  1.00 1.01
```

``` r
r2mlm(mod1,bargraph = F)
```

```
## $Decompositions
##                      total
## fixed           0.00654089
## slope variation 0.00000000
## mean variation  0.05636446
## sigma2          0.93709465
## 
## $R2s
##          total
## f   0.00654089
## v   0.00000000
## m   0.05636446
## fv  0.00654089
## fvm 0.06290535
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
## 1464421.7 1464488.4 -732204.9 1464409.7    492337 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -8.5028 -0.5501  0.1016  0.6560  4.8567 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr 
##  cntry    (Intercept) 0.060960 0.24690       
##           gndr.c      0.007948 0.08915  0.27 
##  Residual             1.007158 1.00357       
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##             Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept) -0.01769    0.04238 33.96993  -0.418    0.679    
## gndr.c      -0.16933    0.01564 32.80431 -10.826 2.31e-12 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##        (Intr)
## gndr.c 0.269
```

``` r
getFE(mod2,round=3)
```

```
##               Est.    SE     df       t     p     LL     UL
## (Intercept) -0.018 0.042 33.970  -0.418 0.679 -0.104  0.068
## gndr.c      -0.169 0.016 32.804 -10.826 0.000 -0.201 -0.137
```

``` r
getVC(mod2)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.25 0.06
## 2    cntry      gndr.c   <NA>  0.09 0.01
## 3    cntry (Intercept) gndr.c  0.27 0.01
## 4 Residual        <NA>   <NA>  1.00 1.01
```

``` r
r2mlm(mod2,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.006616823
## slope variation 0.001834199
## mean variation  0.056189521
## sigma2          0.935359457
## 
## $R2s
##           total
## f   0.006616823
## v   0.001834199
## m   0.056189521
## fv  0.008451022
## fvm 0.064640543
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
## mod1    4 1465124 1465169 -732558   1465116                         
## mod2    6 1464422 1464488 -732205   1464410 706.41  2  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

``` r
lvl2_var_cond_lvl1(mod2,lvl1.var = "gndr.c",lvl1.values = c(0.5,-0.5))
```

```
##   lvl1.value lvl2.cond.var lvl2.cond.sd
## 1        0.5    0.06899616    0.2626712
## 2       -0.5    0.05689866    0.2385344
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
## 1464422.2 1464477.8 -732206.1 1464412.2    492338 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -8.5018 -0.5502  0.1017  0.6560  4.8595 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  cntry    (Intercept) 0.060933 0.24685 
##  cntry.1  gndr.c      0.007916 0.08897 
##  Residual             1.007158 1.00357 
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##             Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept) -0.01768    0.04237 33.97409  -0.417    0.679    
## gndr.c      -0.16927    0.01561 32.79954 -10.843 2.22e-12 ***
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
## (Intercept) -0.018 0.042 33.974  -0.417 0.679 -0.104  0.068
## gndr.c      -0.169 0.016 32.800 -10.843 0.000 -0.201 -0.137
```

``` r
getVC(mod2_norecov)
```

```
##        grp        var1 var2 sdcor vcov
## 1    cntry (Intercept) <NA>  0.25 0.06
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
##              npar     AIC     BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)
## mod2_norecov    5 1464422 1464478 -732206   1464412                     
## mod2            6 1464422 1464488 -732205   1464410 2.5206  1     0.1124
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
## 1420577.0 1420665.6 -710280.5 1420561.0    480356 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -8.3473 -0.5516  0.1028  0.6580  4.8873 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr 
##  cntry    (Intercept) 0.05907  0.2430        
##           gndr.c      0.00709  0.0842   0.37 
##  Residual             0.99416  0.9971        
## Number of obs: 480364, groups:  cntry, 33
## 
## Fixed effects:
##                  Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)     -0.009518   0.042341 32.967058  -0.225   0.8235    
## gndr.c          -0.168491   0.015032 31.921536 -11.209 1.34e-12 ***
## gei.z.cm         0.038560   0.043013 33.016198   0.896   0.3765    
## gndr.c:gei.z.cm -0.034805   0.015438 33.342856  -2.254   0.0309 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c g.z.cm
## gndr.c       0.360              
## gei.z.cm    -0.001  0.000       
## gndr.c:g.z.  0.000 -0.016  0.356
```

``` r
getFE(mod2_GEI,round=3)
```

```
##                   Est.    SE     df       t     p     LL     UL
## (Intercept)     -0.010 0.042 32.967  -0.225 0.824 -0.096  0.077
## gndr.c          -0.168 0.015 31.922 -11.209 0.000 -0.199 -0.138
## gei.z.cm         0.039 0.043 33.016   0.896 0.376 -0.049  0.126
## gndr.c:gei.z.cm -0.035 0.015 33.343  -2.254 0.031 -0.066 -0.003
```

``` r
getVC(mod2_GEI)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.24 0.06
## 2    cntry      gndr.c   <NA>  0.08 0.01
## 3    cntry (Intercept) gndr.c  0.37 0.01
## 4 Residual        <NA>   <NA>  1.00 0.99
```

``` r
r2mlm(mod2_GEI,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.008195253
## slope variation 0.001658234
## mean variation  0.055048171
## sigma2          0.935098342
## 
## $R2s
##           total
## f   0.008195253
## v   0.001658234
## m   0.055048171
## fv  0.009853487
## fvm 0.064901658
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
## Time difference of 35.08357 secs
```

``` r
round(ddsc_mod2_GEI$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.057        0.239        1.007     1.064 0.053   7802.647 0.998   0.998
## 2        0.5         0.069        0.263        1.007     1.076 0.064   6678.029 0.998   0.998
```

``` r
round(ddsc_mod2_GEI$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gei.z.cm gei.z.cm_scaled
## means_y1          -0.089 0.274    1.000           1.000    0.938           0.938    0.097           0.097
## means_y1_scaled   -0.343 1.053    1.000           1.000    0.938           0.938    0.097           0.097
## means_y2           0.079 0.246    0.938           0.938    1.000           1.000    0.245           0.245
## means_y2_scaled    0.303 0.944    0.938           0.938    1.000           1.000    0.245           0.245
## gei.z.cm           0.000 1.000    0.097           0.097    0.245           0.245    1.000           1.000
## gei.z.cm_scaled    0.000 1.000    0.097           0.097    0.245           0.245    1.000           1.000
## diff_score        -0.168 0.096    0.454           0.454    0.118           0.118   -0.352          -0.352
## diff_score_scaled -0.646 0.367    0.454           0.454    0.118           0.118   -0.352          -0.352
##                   diff_score diff_score_scaled
## means_y1               0.454             0.454
## means_y1_scaled        0.454             0.454
## means_y2               0.118             0.118
## means_y2_scaled        0.118             0.118
## gei.z.cm              -0.352            -0.352
## gei.z.cm_scaled       -0.352            -0.352
## diff_score             1.000             1.000
## diff_score_scaled      1.000             1.000
```

``` r
round(ddsc_mod2_GEI$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                       0.364 0.162 33.343   2.254   0.031    0.036    0.693
## w_11                          0.056 0.041 33.047   1.368   0.181   -0.027    0.139
## w_21                          0.021 0.046 33.016   0.457   0.651   -0.073    0.115
## r_xy1                         0.204 0.149 33.047   1.368   0.181   -0.099    0.508
## r_xy2                         0.086 0.188 33.016   0.457   0.651   -0.297    0.469
## b_11                          0.215 0.157 33.047   1.368   0.181   -0.105    0.535
## b_21                          0.081 0.178 33.016   0.457   0.651   -0.281    0.444
## main_effect                   0.039 0.043 33.016   0.896   0.376   -0.049    0.126
## moderator_effect             -0.168 0.015 31.922 -11.209   0.000   -0.199   -0.138
## interaction                  -0.035 0.015 33.343  -2.254   0.031   -0.066   -0.003
## q_b11_b21                     0.137    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.121    NA     NA      NA      NA       NA       NA
## cross_over_point             -4.841    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.004 0.051 33.021  -0.074   0.941   -0.107    0.099
## interaction_vs_main_bscale   -0.014 0.195 33.021  -0.074   0.941   -0.410    0.381
## interaction_vs_main_rscale   -0.027 0.213 33.020  -0.126   0.900   -0.461    0.407
## dadas                        -0.042 0.093 33.016  -0.457   0.675   -0.231    0.146
## dadas_bscale                 -0.163 0.356 33.016  -0.457   0.675   -0.887    0.562
## dadas_rscale                 -0.172 0.377 33.016  -0.457   0.675   -0.938    0.594
## abs_diff                      0.035 0.015 33.343   2.254   0.015    0.003    0.066
## abs_sum                       0.077 0.086 33.016   0.896   0.188   -0.098    0.252
## abs_diff_bscale               0.134 0.059 33.343   2.254   0.015    0.013    0.254
## abs_sum_bscale                0.296 0.331 33.016   0.896   0.188   -0.376    0.969
## abs_diff_rscale               0.118 0.068 33.135   1.736   0.046   -0.020    0.256
## abs_sum_rscale                0.290 0.333 33.016   0.871   0.195   -0.387    0.967
```

``` r
round(ddsc_mod2_GEI$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
##  0.006  0.275  2.521  1.000  0.112
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
##                                    est    se      z pvalue ci.lower ci.upper
## r_xy1_y2                         0.352 0.163  2.158  0.031    0.032    0.671
## r_xy1                            0.245 0.169  1.452  0.146   -0.086    0.576
## r_xy2                            0.097 0.173  0.563  0.574   -0.242    0.437
## b_11                             0.231 0.159  1.452  0.146   -0.081    0.544
## b_21                             0.103 0.182  0.563  0.574   -0.255    0.460
## b_10                             0.303 0.157  1.930  0.054   -0.005    0.611
## b_20                            -0.343 0.180 -1.908  0.056   -0.695    0.009
## res_cov_y1_y2                    0.882 0.223  3.952  0.000    0.444    1.319
## diff_b10_b20                     0.646 0.059 10.975  0.000    0.530    0.761
## diff_b11_b21                     0.129 0.060  2.158  0.031    0.012    0.246
## diff_rxy1_rxy2                   0.148 0.055  2.661  0.008    0.039    0.256
## q_b11_b21                        0.133 0.059  2.245  0.025    0.017    0.249
## q_rxy1_rxy2                      0.152 0.058  2.651  0.008    0.040    0.265
## cross_over_point                -5.009 2.366 -2.117  0.034   -9.646   -0.372
## sum_b11_b21                      0.334 0.337  0.990  0.322   -0.327    0.995
## main_effect                      0.167 0.169  0.990  0.322   -0.163    0.498
## interaction_vs_main_effect      -0.038 0.200 -0.191  0.849   -0.429    0.353
## diff_abs_b11_abs_b21             0.129 0.060  2.158  0.031    0.012    0.246
## abs_diff_b11_b21                 0.129 0.060  2.158  0.015    0.012    0.246
## abs_sum_b11_b21                  0.334 0.337  0.990  0.161   -0.327    0.995
## dadas                           -0.205 0.365 -0.563  0.713   -0.920    0.510
## q_r_equivalence                  0.052 0.058  0.912  0.819       NA       NA
## q_b_equivalence                  0.033 0.059  0.555  0.710       NA       NA
## cross_over_point_equivalence     5.009 2.366  2.117  0.983       NA       NA
## cross_over_point_minimal_effect  5.009 2.366  2.117  0.017       NA       NA
```

``` r
round(ddsc_sem_GEI$variance_test,3)
```

```
##              est    se      z pvalue ci.lower ci.upper
## cov_y1y2   0.905 0.230  3.931  0.000    0.454    1.356
## var_y1     0.865 0.213  4.062  0.000    0.448    1.282
## var_y2     1.075 0.265  4.062  0.000    0.556    1.593
## var_diff  -0.210 0.127 -1.651  0.099   -0.459    0.039
## var_ratio  0.805 0.097  8.306  0.000    0.615    0.995
## cor_y1y2   0.938 0.021 45.074  0.000    0.898    0.979
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
## r_xy1y2                       0.352 0.168 31.000   2.091   0.045    0.009    0.695
## w_11                          0.060 0.046 32.943   1.310   0.199   -0.033    0.154
## w_21                          0.027 0.046 32.943   0.581   0.565   -0.067    0.120
## r_xy1                         0.245 0.187 32.943   1.310   0.199   -0.136    0.626
## r_xy2                         0.097 0.168 32.943   0.581   0.565   -0.244    0.439
## b_11                          0.232 0.177 32.943   1.310   0.199   -0.128    0.592
## b_21                          0.103 0.177 32.943   0.581   0.565   -0.257    0.463
## main_effect                   0.044 0.045 31.000   0.960   0.344   -0.049    0.136
## moderator_effect             -0.168 0.016 31.000 -10.637   0.000   -0.200   -0.136
## interaction                  -0.034 0.016 31.000  -2.091   0.045   -0.066   -0.001
## q_b11_b21                     0.133    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.152    NA     NA      NA      NA       NA       NA
## cross_over_point             -5.009    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.010 0.048 38.660  -0.207   0.837   -0.107    0.087
## interaction_vs_main_bscale   -0.038 0.185 38.660  -0.207   0.837   -0.412    0.336
## interaction_vs_main_rscale   -0.024 0.167 40.063  -0.141   0.888   -0.361    0.314
## dadas                        -0.053 0.092 32.943  -0.581   0.717   -0.241    0.134
## dadas_bscale                 -0.205 0.354 32.943  -0.581   0.717   -0.926    0.515
## dadas_rscale                 -0.195 0.336 32.943  -0.581   0.717   -0.878    0.488
## abs_diff                      0.034 0.016 31.000   2.091   0.022    0.001    0.066
## abs_sum                       0.087 0.091 31.000   0.960   0.172   -0.098    0.272
## abs_diff_bscale               0.129 0.062 31.000   2.091   0.022    0.003    0.255
## abs_sum_bscale                0.335 0.349 31.000   0.960   0.172   -0.376    1.045
## abs_diff_rscale               0.148 0.065 36.760   2.281   0.014    0.016    0.279
## abs_sum_rscale                0.343 0.350 31.006   0.980   0.167   -0.370    1.055
```

``` r
# country-time multilevel model


mod2_GEI_cntry_year<-
  lmer(uni.z.wt~gndr.c+gei.z.cm+gndr.c:gei.z.cm+
                 (gndr.c|cntry),data=diff_dat_cntry_year,REML=F,
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
```

```
## boundary (singular) fit: see help('isSingular')
```

``` r
summary(mod2_GEI_cntry_year)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: uni.z.wt ~ gndr.c + gei.z.cm + gndr.c:gei.z.cm + (gndr.c | cntry)
##    Data: diff_dat_cntry_year
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##    -282.6    -248.4     149.3    -298.6       526 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -5.8697 -0.4996  0.0755  0.5738  3.0489 
## 
## Random effects:
##  Groups   Name        Variance  Std.Dev. Corr 
##  cntry    (Intercept) 0.0533356 0.2309        
##           gndr.c      0.0005615 0.0237   1.00 
##  Residual             0.0271177 0.1647        
## Number of obs: 534, groups:  cntry, 33
## 
## Fixed effects:
##                   Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)      -0.003277   0.040991  32.551253  -0.080   0.9368    
## gndr.c           -0.167367   0.015213 151.554033 -11.002   <2e-16 ***
## gei.z.cm          0.040069   0.042109  34.008721   0.952   0.3480    
## gndr.c:gei.z.cm  -0.035752   0.017318 202.024905  -2.064   0.0403 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c g.z.cm
## gndr.c       0.268              
## gei.z.cm    -0.015 -0.002       
## gndr.c:g.z. -0.002 -0.213  0.240
## optimizer (bobyqa) convergence code: 0 (OK)
## boundary (singular) fit: see help('isSingular')
```

``` r
getFE(mod2_GEI_cntry_year,round=3)
```

```
##                   Est.    SE      df       t     p     LL     UL
## (Intercept)     -0.003 0.041  32.551  -0.080 0.937 -0.087  0.080
## gndr.c          -0.167 0.015 151.554 -11.002 0.000 -0.197 -0.137
## gei.z.cm         0.040 0.042  34.009   0.952 0.348 -0.046  0.126
## gndr.c:gei.z.cm -0.036 0.017 202.025  -2.064 0.040 -0.070 -0.002
```

``` r
getVC(mod2_GEI_cntry_year)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.23 0.05
## 2    cntry      gndr.c   <NA>  0.02 0.00
## 3    cntry (Intercept) gndr.c  1.00 0.01
## 4 Residual        <NA>   <NA>  0.16 0.03
```

``` r
r2mlm(mod2_GEI,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.008195253
## slope variation 0.001658234
## mean variation  0.055048171
## sigma2          0.935098342
## 
## $R2s
##           total
## f   0.008195253
## v   0.001658234
## m   0.055048171
## fv  0.009853487
## fvm 0.064901658
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
## 1       -0.5         0.052        0.229        0.027     0.079 0.660      8.029 0.995   0.940
## 2        0.5         0.059        0.244        0.027     0.086 0.687      8.029 0.995   0.946
```

``` r
round(ddsc_mod2_GEI_cntry_year$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gei.z.cm gei.z.cm_scaled
## means_y1          -0.089 0.261    1.000           1.000    0.934           0.934    0.077           0.077
## means_y1_scaled   -0.359 1.050    1.000           1.000    0.934           0.934    0.077           0.077
## means_y2           0.083 0.235    0.934           0.934    1.000           1.000    0.237           0.237
## means_y2_scaled    0.334 0.947    0.934           0.934    1.000           1.000    0.237           0.237
## gei.z.cm           0.000 1.000    0.077           0.077    0.237           0.237    1.000           1.000
## gei.z.cm_scaled    0.000 1.000    0.077           0.077    0.237           0.237    1.000           1.000
## diff_score        -0.172 0.093    0.441           0.441    0.093           0.093   -0.384          -0.384
## diff_score_scaled -0.694 0.376    0.441           0.441    0.093           0.093   -0.384          -0.384
##                   diff_score diff_score_scaled
## means_y1               0.441             0.441
## means_y1_scaled        0.441             0.441
## means_y2               0.093             0.093
## means_y2_scaled        0.093             0.093
## gei.z.cm              -0.384            -0.384
## gei.z.cm_scaled       -0.384            -0.384
## diff_score             1.000             1.000
## diff_score_scaled      1.000             1.000
```

``` r
round(ddsc_mod2_GEI_cntry_year$results,3)
```

```
##                            estimate    SE      df t.ratio p.value ci.lower ci.upper
## r_xy1y2                       0.384 0.186 202.025   2.064   0.040    0.017    0.750
## w_11                          0.058 0.041  34.466   1.416   0.166   -0.025    0.141
## w_21                          0.022 0.045  34.183   0.493   0.625   -0.069    0.114
## r_xy1                         0.222 0.157  34.466   1.416   0.166   -0.097    0.541
## r_xy2                         0.095 0.192  34.183   0.493   0.625   -0.295    0.484
## b_11                          0.234 0.165  34.466   1.416   0.166   -0.102    0.569
## b_21                          0.090 0.182  34.183   0.493   0.625   -0.279    0.459
## main_effect                   0.040 0.042  34.009   0.952   0.348   -0.046    0.126
## moderator_effect             -0.167 0.015 151.554 -11.002   0.000   -0.197   -0.137
## interaction                  -0.036 0.017 202.025  -2.064   0.040   -0.070   -0.002
## q_b11_b21                     0.149    NA      NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.131    NA      NA      NA      NA       NA       NA
## cross_over_point             -4.681    NA      NA      NA      NA       NA       NA
## interaction_vs_main          -0.004 0.049  34.996  -0.088   0.931   -0.104    0.096
## interaction_vs_main_bscale   -0.017 0.199  34.996  -0.088   0.931   -0.421    0.386
## interaction_vs_main_rscale   -0.031 0.217  34.838  -0.141   0.889   -0.471    0.410
## dadas                        -0.044 0.090  34.183  -0.493   0.688   -0.227    0.138
## dadas_bscale                 -0.179 0.363  34.183  -0.493   0.688   -0.917    0.559
## dadas_rscale                 -0.189 0.383  34.183  -0.493   0.688   -0.967    0.589
## abs_diff                      0.036 0.017 202.025   2.064   0.020    0.002    0.070
## abs_sum                       0.080 0.084  34.009   0.952   0.174   -0.091    0.251
## abs_diff_bscale               0.144 0.070 202.025   2.064   0.020    0.006    0.282
## abs_sum_bscale                0.324 0.340  34.009   0.952   0.174   -0.367    1.015
## abs_diff_rscale               0.128 0.076  87.568   1.676   0.049   -0.024    0.280
## abs_sum_rscale                0.317 0.342  34.004   0.927   0.180   -0.378    1.012
```

``` r
round(ddsc_mod2_GEI_cntry_year$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
##  0.003  0.649  0.901  1.000  0.343
```

### Comparisons across models


``` r
# full multilevel model (individuals as entries)
round(ddsc_mod2_GEI$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df  t.ratio p.value ci.lower ci.upper
## r_xy1              0.2041 0.1492 33.0468   1.3680  0.1806  -0.0994   0.5075
## r_xy2              0.0860 0.1883 33.0157   0.4567  0.6509  -0.2971   0.4690
## b_11               0.2151 0.1573 33.0468   1.3680  0.1806  -0.1048   0.5350
## b_21               0.0813 0.1781 33.0157   0.4567  0.6509  -0.2809   0.4436
## main_effect        0.0386 0.0430 33.0162   0.8965  0.3765  -0.0490   0.1261
## moderator_effect  -0.1685 0.0150 31.9215 -11.2088  0.0000  -0.1991  -0.1379
## interaction       -0.0348 0.0154 33.3429  -2.2545  0.0309  -0.0662  -0.0034
## q_b11_b21          0.1370     NA      NA       NA      NA       NA       NA
```

``` r
# country-level structural path model (country-gender means as entries)
round(ddsc_sem_GEI$results[c("r_xy1","r_xy2","b_11","b_21","q_b11_b21","diff_b11_b21"),],4)
```

```
##                 est     se      z pvalue ci.lower ci.upper
## r_xy1        0.2451 0.1688 1.4524 0.1464  -0.0857   0.5759
## r_xy2        0.0975 0.1732 0.5625 0.5738  -0.2421   0.4370
## b_11         0.2315 0.1594 1.4524 0.1464  -0.0809   0.5439
## b_21         0.1026 0.1824 0.5625 0.5738  -0.2549   0.4600
## q_b11_b21    0.1328 0.0592 2.2449 0.0248   0.0169   0.2488
## diff_b11_b21 0.1289 0.0597 2.1575 0.0310   0.0118   0.2460
```

``` r
# multilevel model (country-gender means as entries)
round(ddsc_mod2_GEI_ri$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df  t.ratio p.value ci.lower ci.upper
## r_xy1              0.2451 0.1871 32.9433   1.3100  0.1993  -0.1356   0.6258
## r_xy2              0.0975 0.1679 32.9433   0.5806  0.5655  -0.2441   0.4390
## b_11               0.2318 0.1770 32.9433   1.3100  0.1993  -0.1282   0.5919
## b_21               0.1027 0.1770 32.9433   0.5806  0.5655  -0.2573   0.4628
## main_effect        0.0435 0.0453 31.0000   0.9600  0.3445  -0.0489   0.1360
## moderator_effect  -0.1682 0.0158 31.0000 -10.6367  0.0000  -0.2005  -0.1360
## interaction       -0.0336 0.0161 31.0000  -2.0911  0.0448  -0.0663  -0.0008
## q_b11_b21          0.1330     NA      NA       NA      NA       NA       NA
```

``` r
# multilevel model (country-year-gender means as entries)
round(ddsc_mod2_GEI_cntry_year$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE       df  t.ratio p.value ci.lower ci.upper
## r_xy1              0.2224 0.1570  34.4664   1.4165  0.1656  -0.0965   0.5413
## r_xy2              0.0945 0.1915  34.1827   0.4934  0.6249  -0.2946   0.4836
## b_11               0.2339 0.1652  34.4664   1.4165  0.1656  -0.1015   0.5694
## b_21               0.0896 0.1816  34.1827   0.4934  0.6249  -0.2793   0.4585
## main_effect        0.0401 0.0421  34.0087   0.9516  0.3480  -0.0455   0.1256
## moderator_effect  -0.1674 0.0152 151.5540 -11.0018  0.0000  -0.1974  -0.1373
## interaction       -0.0358 0.0173 202.0249  -2.0645  0.0403  -0.0699  -0.0016
## q_b11_b21          0.1485     NA       NA       NA      NA       NA       NA
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
## Time difference of 1.391205 hours
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
##                    Estimate         SE        2.5%        97.5%
## X.Intercept.    -0.01060317 0.04175051 -0.09450353  0.072995442
## gndr.c          -0.16801720 0.01484458 -0.19747661 -0.140055677
## gei.z.cm         0.03946107 0.04298525 -0.04164199  0.128397642
## gndr.c.gei.z.cm -0.03526654 0.01604701 -0.06612411 -0.004269414
## w11              0.05709434 0.04053434 -0.02131309  0.138549556
## w21              0.02182780 0.04670314 -0.06824273  0.115944793
## b11              0.21946884 0.15581270 -0.08192685  0.532580086
## b21              0.08390538 0.17952540 -0.26232288  0.445688094
## r_xy1            0.20817978 0.14779799 -0.07771269  0.505185199
## r_xy2            0.08871623 0.18981879 -0.27736359  0.471242348
## q_b              0.14234751 0.06353248  0.01697412  0.269322343
## q                0.12386307 0.07297673 -0.02453539  0.260618558
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
## [1] 0.1423475
## 
## $se
## [1] 0.06353248
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
## [1] 3.814545
## 
## $p_low
## [1] 6.821705e-05
## 
## $z_high
## [1] 0.6665489
## 
## $p_high
## [1] 0.7474698
## 
## $ci_level
## [1] 0.9
## 
## $ci_lower
## [1] 0.03784587
## 
## $ci_upper
## [1] 0.2468491
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
## [1] 0.1238631
## 
## $se
## [1] 0.07297673
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
## [1] 3.067595
## 
## $p_low
## [1] 0.001078943
## 
## $z_high
## [1] 0.3269956
## 
## $p_high
## [1] 0.6281644
## 
## $ci_level
## [1] 0.9
## 
## $ci_lower
## [1] 0.003827038
## 
## $ci_upper
## [1] 0.2438991
## 
## $equivalent
## [1] FALSE
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
  xlab("Gender Equality Index (Average 2002-2023)")+
  #ylim(c(min.y.pred,max.y.pred))+
  #xlim(c(0.60,1.00))+
  ylab("Value universalism (Average 2002-2023)")+
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
  xlab("Gender Equality Index (Average 2002-2023)")+
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

![](Analysis_universalism_files/figure-html/unnamed-chunk-23-1.png)<!-- -->

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
## 1081996.4 1082082.8 -540990.2 1081980.4    363844 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -8.5587 -0.5513  0.1057  0.6592  5.0244 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr 
##  cntry    (Intercept) 0.070835 0.26615       
##           gndr.c      0.005351 0.07315  0.37 
##  Residual             1.004785 1.00239       
## Number of obs: 363852, groups:  cntry, 34
## 
## Fixed effects:
##                  Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)      -0.01572    0.04569 33.94344  -0.344 0.732950    
## gndr.c           -0.16812    0.01313 31.82480 -12.803 4.27e-14 ***
## gggi.z.cm         0.03347    0.04639 33.99558   0.722 0.475532    
## gndr.c:gggi.z.cm -0.05130    0.01356 33.92765  -3.783 0.000601 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c ggg.z.
## gndr.c       0.352              
## gggi.z.cm    0.000  0.000       
## gndr.c:gg..  0.000 -0.011  0.346
```

``` r
getFE(mod2_GGGI,round=3)
```

```
##                    Est.    SE     df       t     p     LL     UL
## (Intercept)      -0.016 0.046 33.943  -0.344 0.733 -0.109  0.077
## gndr.c           -0.168 0.013 31.825 -12.803 0.000 -0.195 -0.141
## gggi.z.cm         0.033 0.046 33.996   0.722 0.476 -0.061  0.128
## gndr.c:gggi.z.cm -0.051 0.014 33.928  -3.783 0.001 -0.079 -0.024
```

``` r
getVC(mod2_GGGI)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.27 0.07
## 2    cntry      gndr.c   <NA>  0.07 0.01
## 3    cntry (Intercept) gndr.c  0.37 0.01
## 4 Residual        <NA>   <NA>  1.00 1.00
```

``` r
r2mlm(mod2_GGGI,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.007988779
## slope variation 0.001224962
## mean variation  0.064762144
## sigma2          0.926024115
## 
## $R2s
##           total
## f   0.007988779
## v   0.001224962
## m   0.064762144
## fv  0.009213742
## fvm 0.073975885
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
## Time difference of 38.11007 secs
```

``` r
round(ddsc_mod2_GGGI$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.057        0.239        1.007     1.064 0.053   7802.647 0.998   0.998
## 2        0.5         0.069        0.263        1.007     1.076 0.064   6678.029 0.998   0.998
```

``` r
round(ddsc_mod2_GGGI$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gggi.z.cm
## means_y1          -0.096 0.300    1.000           1.000    0.952           0.952     0.054
## means_y1_scaled   -0.335 1.044    1.000           1.000    0.952           0.952     0.054
## means_y2           0.075 0.274    0.952           0.952    1.000           1.000     0.231
## means_y2_scaled    0.262 0.954    0.952           0.952    1.000           1.000     0.231
## gggi.z.cm          0.000 1.000    0.054           0.054    0.231           0.231     1.000
## gggi.z.cm_scaled   0.000 1.000    0.054           0.054    0.231           0.231     1.000
## diff_score        -0.171 0.093    0.424           0.424    0.125           0.125    -0.507
## diff_score_scaled -0.597 0.323    0.424           0.424    0.125           0.125    -0.507
##                   gggi.z.cm_scaled diff_score diff_score_scaled
## means_y1                     0.054      0.424             0.424
## means_y1_scaled              0.054      0.424             0.424
## means_y2                     0.231      0.125             0.125
## means_y2_scaled              0.231      0.125             0.125
## gggi.z.cm                    1.000     -0.507            -0.507
## gggi.z.cm_scaled             1.000     -0.507            -0.507
## diff_score                  -0.507      1.000             1.000
## diff_score_scaled           -0.507      1.000             1.000
```

``` r
round(ddsc_mod2_GGGI$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                       0.553 0.146 33.928   3.783   0.001    0.256    0.850
## w_11                          0.059 0.045 34.006   1.328   0.193   -0.031    0.150
## w_21                          0.008 0.049 33.980   0.159   0.874   -0.092    0.108
## r_xy1                         0.197 0.149 34.006   1.328   0.193   -0.105    0.499
## r_xy2                         0.029 0.180 33.980   0.159   0.874   -0.337    0.394
## b_11                          0.206 0.155 34.006   1.328   0.193   -0.109    0.522
## b_21                          0.027 0.172 33.980   0.159   0.874   -0.321    0.376
## main_effect                   0.033 0.046 33.996   0.722   0.476   -0.061    0.128
## moderator_effect             -0.168 0.013 31.825 -12.803   0.000   -0.195   -0.141
## interaction                  -0.051 0.014 33.928  -3.783   0.001   -0.079   -0.024
## q_b11_b21                     0.182    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.171    NA     NA      NA      NA       NA       NA
## cross_over_point             -3.277    NA     NA      NA      NA       NA       NA
## interaction_vs_main           0.018 0.053 33.953   0.339   0.737   -0.089    0.125
## interaction_vs_main_bscale    0.062 0.184 33.953   0.339   0.737   -0.311    0.436
## interaction_vs_main_rscale    0.056 0.199 33.957   0.280   0.781   -0.349    0.460
## dadas                        -0.016 0.098 33.980  -0.159   0.563   -0.215    0.184
## dadas_bscale                 -0.055 0.343 33.980  -0.159   0.563   -0.752    0.643
## dadas_rscale                 -0.057 0.359 33.980  -0.159   0.563   -0.788    0.673
## abs_diff                      0.051 0.014 33.928   3.783   0.000    0.024    0.079
## abs_sum                       0.067 0.093 33.996   0.722   0.238   -0.122    0.256
## abs_diff_bscale               0.179 0.047 33.928   3.783   0.000    0.083    0.275
## abs_sum_bscale                0.234 0.324 33.996   0.722   0.238   -0.424    0.892
## abs_diff_rscale               0.169 0.054 33.750   3.107   0.002    0.058    0.279
## abs_sum_rscale                0.226 0.325 33.995   0.695   0.246   -0.435    0.887
```

``` r
round(ddsc_mod2_GGGI$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
##  0.006  0.275  2.521  1.000  0.112
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
##                                    est    se      z pvalue ci.lower ci.upper
## r_xy1_y2                         0.507 0.148  3.434  0.001    0.218    0.797
## r_xy1                            0.231 0.167  1.386  0.166   -0.096    0.558
## r_xy2                            0.054 0.171  0.315  0.753   -0.282    0.390
## b_11                             0.220 0.159  1.386  0.166   -0.091    0.532
## b_21                             0.056 0.179  0.315  0.753   -0.294    0.407
## b_10                             0.262 0.157  1.670  0.095   -0.045    0.569
## b_20                            -0.335 0.176 -1.900  0.057   -0.680    0.011
## res_cov_y1_y2                    0.908 0.224  4.053  0.000    0.469    1.347
## diff_b10_b20                     0.597 0.047 12.668  0.000    0.504    0.689
## diff_b11_b21                     0.164 0.048  3.434  0.001    0.070    0.258
## diff_rxy1_rxy2                   0.177 0.044  4.044  0.000    0.091    0.263
## q_b11_b21                        0.168 0.046  3.621  0.000    0.077    0.259
## q_rxy1_rxy2                      0.181 0.045  4.022  0.000    0.093    0.270
## cross_over_point                -3.634 1.096 -3.315  0.001   -5.783   -1.485
## sum_b11_b21                      0.277 0.335  0.826  0.409   -0.380    0.934
## main_effect                      0.138 0.168  0.826  0.409   -0.190    0.467
## interaction_vs_main_effect       0.026 0.192  0.134  0.894   -0.351    0.403
## diff_abs_b11_abs_b21             0.164 0.048  3.434  0.001    0.070    0.258
## abs_diff_b11_b21                 0.164 0.048  3.434  0.000    0.070    0.258
## abs_sum_b11_b21                  0.277 0.335  0.826  0.204   -0.380    0.934
## dadas                           -0.113 0.358 -0.315  0.624   -0.814    0.588
## q_r_equivalence                  0.081 0.045  1.806  0.965       NA       NA
## q_b_equivalence                  0.068 0.046  1.463  0.928       NA       NA
## cross_over_point_equivalence     3.634 1.096  3.315  1.000       NA       NA
## cross_over_point_minimal_effect  3.634 1.096  3.315  0.000       NA       NA
```

``` r
round(ddsc_sem_GGGI$variance_test,3)
```

```
##              est    se      z pvalue ci.lower ci.upper
## cov_y1y2   0.920 0.229  4.020  0.000    0.471    1.368
## var_y1     0.882 0.214  4.123  0.000    0.463    1.302
## var_y2     1.059 0.257  4.123  0.000    0.555    1.562
## var_diff  -0.176 0.110 -1.595  0.111   -0.393    0.040
## var_ratio  0.834 0.088  9.488  0.000    0.661    1.006
## cor_y1y2   0.952 0.016 58.766  0.000    0.920    0.983
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
## r_xy1y2                       0.507 0.152 32.000   3.332   0.002    0.197    0.818
## w_11                          0.063 0.050 33.301   1.264   0.215   -0.039    0.165
## w_21                          0.016 0.050 33.301   0.323   0.749   -0.086    0.118
## r_xy1                         0.231 0.183 33.301   1.264   0.215   -0.141    0.603
## r_xy2                         0.054 0.167 33.301   0.323   0.749   -0.286    0.394
## b_11                          0.221 0.175 33.301   1.264   0.215   -0.134    0.576
## b_21                          0.056 0.175 33.301   0.323   0.749   -0.299    0.412
## main_effect                   0.040 0.050 32.000   0.801   0.429   -0.061    0.141
## moderator_effect             -0.171 0.014 32.000 -12.290   0.000   -0.199   -0.143
## interaction                  -0.047 0.014 32.000  -3.332   0.002   -0.076   -0.018
## q_b11_b21                     0.168    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.181    NA     NA      NA      NA       NA       NA
## cross_over_point             -3.634    NA     NA      NA      NA       NA       NA
## interaction_vs_main           0.007 0.052 37.173   0.143   0.887   -0.097    0.112
## interaction_vs_main_bscale    0.026 0.180 37.173   0.143   0.887   -0.338    0.390
## interaction_vs_main_rscale    0.035 0.165 37.966   0.211   0.834   -0.299    0.368
## dadas                        -0.032 0.100 33.301  -0.323   0.626   -0.236    0.171
## dadas_bscale                 -0.113 0.349 33.301  -0.323   0.626   -0.823    0.598
## dadas_rscale                 -0.108 0.334 33.301  -0.323   0.626   -0.787    0.572
## abs_diff                      0.047 0.014 32.000   3.332   0.001    0.018    0.076
## abs_sum                       0.079 0.099 32.000   0.801   0.214   -0.122    0.281
## abs_diff_bscale               0.164 0.049 32.000   3.332   0.001    0.064    0.265
## abs_sum_bscale                0.277 0.346 32.000   0.801   0.214   -0.427    0.982
## abs_diff_rscale               0.177 0.052 38.443   3.418   0.001    0.072    0.282
## abs_sum_rscale                0.285 0.347 32.003   0.823   0.208   -0.421    0.991
```

``` r
# country-time multilevel model


mod2_GGGI_cntry_year<-
  lmer(uni.z.wt~gndr.c+gggi.z.cm+gndr.c:gggi.z.cm+
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
## Formula: uni.z.wt ~ gndr.c + gggi.z.cm + gndr.c:gggi.z.cm + (gndr.c |      cntry)
##    Data: diff_dat_cntry_year
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##    -258.1    -226.2     137.1    -274.1       392 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -5.0692 -0.4977 -0.0009  0.5918  4.5618 
## 
## Random effects:
##  Groups   Name        Variance  Std.Dev. Corr 
##  cntry    (Intercept) 0.0649176 0.25479       
##           gndr.c      0.0003761 0.01939  1.00 
##  Residual             0.0219540 0.14817       
## Number of obs: 400, groups:  cntry, 34
## 
## Fixed effects:
##                   Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)       -0.01050    0.04453  33.34416  -0.236  0.81508    
## gndr.c            -0.16750    0.01534 174.74213 -10.918  < 2e-16 ***
## gggi.z.cm          0.03238    0.04550  34.14440   0.712  0.48146    
## gndr.c:gggi.z.cm  -0.04384    0.01644 193.60984  -2.666  0.00832 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c ggg.z.
## gndr.c       0.215              
## gggi.z.cm   -0.011 -0.001       
## gndr.c:gg.. -0.001 -0.138  0.204
## optimizer (bobyqa) convergence code: 0 (OK)
## boundary (singular) fit: see help('isSingular')
```

``` r
getFE(mod2_GGGI_cntry_year,round=3)
```

```
##                    Est.    SE      df       t     p     LL     UL
## (Intercept)      -0.010 0.045  33.344  -0.236 0.815 -0.101  0.080
## gndr.c           -0.167 0.015 174.742 -10.918 0.000 -0.198 -0.137
## gggi.z.cm         0.032 0.045  34.144   0.712 0.481 -0.060  0.125
## gndr.c:gggi.z.cm -0.044 0.016 193.610  -2.666 0.008 -0.076 -0.011
```

``` r
getVC(mod2_GGGI_cntry_year)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.25 0.06
## 2    cntry      gndr.c   <NA>  0.02 0.00
## 3    cntry (Intercept) gndr.c  1.00 0.00
## 4 Residual        <NA>   <NA>  0.15 0.02
```

``` r
r2mlm(mod2_GGGI,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.007988779
## slope variation 0.001224962
## mean variation  0.064762144
## sigma2          0.926024115
## 
## $R2s
##           total
## f   0.007988779
## v   0.001224962
## m   0.064762144
## fv  0.009213742
## fvm 0.073975885
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
## 1       -0.5         0.052        0.229        0.027     0.079 0.660      8.029 0.995   0.940
## 2        0.5         0.059        0.244        0.027     0.086 0.687      8.029 0.995   0.946
```

``` r
round(ddsc_mod2_GGGI_cntry_year$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gggi.z.cm
## means_y1          -0.096 0.282    1.000           1.000    0.942           0.942     0.013
## means_y1_scaled   -0.352 1.035    1.000           1.000    0.942           0.942     0.013
## means_y2           0.074 0.262    0.942           0.942    1.000           1.000     0.221
## means_y2_scaled    0.274 0.963    0.942           0.942    1.000           1.000     0.221
## gggi.z.cm          0.000 1.000    0.013           0.013    0.221           0.221     1.000
## gggi.z.cm_scaled   0.000 1.000    0.013           0.013    0.221           0.221     1.000
## diff_score        -0.170 0.094    0.368           0.368    0.035           0.035    -0.574
## diff_score_scaled -0.625 0.347    0.368           0.368    0.035           0.035    -0.574
##                   gggi.z.cm_scaled diff_score diff_score_scaled
## means_y1                     0.013      0.368             0.368
## means_y1_scaled              0.013      0.368             0.368
## means_y2                     0.221      0.035             0.035
## means_y2_scaled              0.221      0.035             0.035
## gggi.z.cm                    1.000     -0.574            -0.574
## gggi.z.cm_scaled             1.000     -0.574            -0.574
## diff_score                  -0.574      1.000             1.000
## diff_score_scaled           -0.574      1.000             1.000
```

``` r
round(ddsc_mod2_GGGI_cntry_year$results,3)
```

```
##                            estimate    SE      df t.ratio p.value ci.lower ci.upper
## r_xy1y2                       0.464 0.174 193.610   2.666   0.008    0.121    0.808
## w_11                          0.054 0.045  34.115   1.219   0.231   -0.036    0.145
## w_21                          0.010 0.048  34.029   0.219   0.828   -0.087    0.108
## r_xy1                         0.193 0.158  34.115   1.219   0.231   -0.129    0.514
## r_xy2                         0.040 0.183  34.029   0.219   0.828   -0.331    0.411
## b_11                          0.200 0.164  34.115   1.219   0.231   -0.133    0.533
## b_21                          0.038 0.176  34.029   0.219   0.828   -0.319    0.396
## main_effect                   0.032 0.045  34.144   0.712   0.481   -0.060    0.125
## moderator_effect             -0.167 0.015 174.742 -10.918   0.000   -0.198   -0.137
## interaction                  -0.044 0.016 193.610  -2.666   0.008   -0.076   -0.011
## q_b11_b21                     0.164    NA      NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.155    NA      NA      NA      NA       NA       NA
## cross_over_point             -3.821    NA      NA      NA      NA       NA       NA
## interaction_vs_main           0.011 0.051  34.040   0.223   0.825   -0.093    0.116
## interaction_vs_main_bscale    0.042 0.189  34.040   0.223   0.825   -0.342    0.427
## interaction_vs_main_rscale    0.037 0.202  34.021   0.181   0.857   -0.373    0.446
## dadas                        -0.021 0.096  34.029  -0.219   0.586   -0.215    0.174
## dadas_bscale                 -0.077 0.352  34.029  -0.219   0.586   -0.793    0.639
## dadas_rscale                 -0.080 0.365  34.029  -0.219   0.586   -0.822    0.662
## abs_diff                      0.044 0.016 193.610   2.666   0.004    0.011    0.076
## abs_sum                       0.065 0.091  34.144   0.712   0.241   -0.120    0.250
## abs_diff_bscale               0.161 0.060 193.610   2.666   0.004    0.042    0.281
## abs_sum_bscale                0.238 0.335  34.144   0.712   0.241   -0.442    0.919
## abs_diff_rscale               0.153 0.064  95.388   2.385   0.010    0.026    0.280
## abs_sum_rscale                0.233 0.336  34.143   0.693   0.246   -0.449    0.915
```

``` r
round(ddsc_mod2_GGGI_cntry_year$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
##  0.003  0.649  0.901  1.000  0.343
```

### Comparisons across models


``` r
# full multilevel model (individuals as entries)
round(ddsc_mod2_GGGI$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df  t.ratio p.value ci.lower ci.upper
## r_xy1              0.1974 0.1486 34.0059   1.3284  0.1929  -0.1046   0.4993
## r_xy2              0.0286 0.1797 33.9799   0.1592  0.8745  -0.3366   0.3938
## b_11               0.2063 0.1553 34.0059   1.3284  0.1929  -0.1093   0.5220
## b_21               0.0273 0.1715 33.9799   0.1592  0.8745  -0.3213   0.3759
## main_effect        0.0335 0.0464 33.9956   0.7215  0.4755  -0.0608   0.1278
## moderator_effect  -0.1681 0.0131 31.8248 -12.8025  0.0000  -0.1949  -0.1414
## interaction       -0.0513 0.0136 33.9276  -3.7832  0.0006  -0.0789  -0.0237
## q_b11_b21          0.1820     NA      NA       NA      NA       NA       NA
```

``` r
# country-level structural path model (country-gender means as entries)
round(ddsc_sem_GGGI$results[c("r_xy1","r_xy2","b_11","b_21","q_b11_b21","diff_b11_b21"),],4)
```

```
##                 est     se      z pvalue ci.lower ci.upper
## r_xy1        0.2312 0.1669 1.3859 0.1658  -0.0958   0.5583
## r_xy2        0.0539 0.1712 0.3150 0.7527  -0.2817   0.3896
## b_11         0.2205 0.1591 1.3859 0.1658  -0.0913   0.5323
## b_21         0.0563 0.1789 0.3150 0.7527  -0.2942   0.4069
## q_b11_b21    0.1678 0.0463 3.6214 0.0003   0.0770   0.2586
## diff_b11_b21 0.1641 0.0478 3.4342 0.0006   0.0705   0.2578
```

``` r
# multilevel model (country-gender means as entries)
round(ddsc_mod2_GGGI_ri$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df  t.ratio p.value ci.lower ci.upper
## r_xy1              0.2312 0.1830 33.3013   1.2637  0.2151  -0.1409   0.6034
## r_xy2              0.0539 0.1671 33.3013   0.3229  0.7488  -0.2858   0.3937
## b_11               0.2207 0.1747 33.3013   1.2637  0.2151  -0.1345   0.5759
## b_21               0.0564 0.1747 33.3013   0.3229  0.7488  -0.2988   0.4116
## main_effect        0.0397 0.0495 32.0000   0.8014  0.4288  -0.0612   0.1406
## moderator_effect  -0.1711 0.0139 32.0000 -12.2899  0.0000  -0.1995  -0.1427
## interaction       -0.0471 0.0141 32.0000  -3.3316  0.0022  -0.0759  -0.0183
## q_b11_b21          0.1679     NA      NA       NA      NA       NA       NA
```

``` r
# multilevel model (country-year-gender means as entries)
round(ddsc_mod2_GGGI_cntry_year$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE       df  t.ratio p.value ci.lower ci.upper
## r_xy1              0.1929 0.1582  34.1147   1.2188  0.2313  -0.1287   0.5144
## r_xy2              0.0399 0.1826  34.0290   0.2186  0.8283  -0.3312   0.4111
## b_11               0.1998 0.1639  34.1147   1.2188  0.2313  -0.1333   0.5329
## b_21               0.0385 0.1761  34.0290   0.2186  0.8283  -0.3193   0.3963
## main_effect        0.0324 0.0455  34.1444   0.7117  0.4815  -0.0601   0.1248
## moderator_effect  -0.1675 0.0153 174.7421 -10.9176  0.0000  -0.1978  -0.1372
## interaction       -0.0438 0.0164 193.6098  -2.6663  0.0083  -0.0763  -0.0114
## q_b11_b21          0.1640     NA       NA       NA      NA       NA       NA
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
## Time difference of 1.116148 hours
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
##                      Estimate         SE        2.5%       97.5%
## X.Intercept.     -0.016746108 0.04609506 -0.10079138  0.07424214
## gndr.c           -0.167619694 0.01366399 -0.19505224 -0.14146565
## gggi.z.cm         0.034810890 0.04505140 -0.05887457  0.12153134
## gndr.c.gggi.z.cm -0.050897676 0.01392988 -0.07744833 -0.02266103
## w11               0.060259728 0.04330455 -0.02722169  0.14461008
## w21               0.009362052 0.04775975 -0.08963368  0.10080405
## b11               0.210303523 0.15113078 -0.09500237  0.50468217
## b21               0.032673106 0.16667921 -0.31281718  0.35180124
## r_xy1             0.201153906 0.14455558 -0.09086913  0.48272510
## r_xy2             0.034230084 0.17462200 -0.32772392  0.36856569
## q_b               0.185611270 0.05022991  0.08241674  0.28015715
## q                 0.173585546 0.05563629  0.06582108  0.27847286
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
## [1] 0.1856113
## 
## $se
## [1] 0.05022991
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
## [1] 5.68608
## 
## $p_low
## [1] 6.499422e-09
## 
## $z_high
## [1] 1.704388
## 
## $p_high
## [1] 0.9558457
## 
## $ci_level
## [1] 0.9
## 
## $ci_lower
## [1] 0.1029904
## 
## $ci_upper
## [1] 0.2682321
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
## [1] 0.1735855
## 
## $se
## [1] 0.05563629
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
## [1] 4.917394
## 
## $p_low
## [1] 4.385195e-07
## 
## $z_high
## [1] 1.322618
## 
## $p_high
## [1] 0.9070188
## 
## $ci_level
## [1] 0.9
## 
## $ci_lower
## [1] 0.082072
## 
## $ci_upper
## [1] 0.2650991
## 
## $equivalent
## [1] FALSE
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
  xlab("Global Gender Gap Index (Average 2002-2023)")+
  #ylim(c(min.y.pred,max.y.pred))+
  #xlim(c(0.60,1.00))+
  ylab("Value universalism (Average 2002-2023)")+
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
  xlab("Global Gender Gap Index (Average 2002-2023)")+
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

![](Analysis_universalism_files/figure-html/unnamed-chunk-29-1.png)<!-- -->

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
## 1464413.6 1464502.5 -732198.8 1464397.6    492335 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -8.5037 -0.5501  0.1018  0.6560  4.8560 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr 
##  cntry    (Intercept) 0.047671 0.21834       
##           gndr.c      0.006662 0.08162  0.10 
##  Residual             1.007158 1.00357       
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##                 Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)     -0.01766    0.03748 33.95534  -0.471  0.64046    
## gndr.c          -0.16901    0.01438 33.40369 -11.752 2.06e-13 ***
## gdi.z.cm        -0.11702    0.03807 34.03195  -3.074  0.00414 ** 
## gndr.c:gdi.z.cm -0.03743    0.01481 35.34270  -2.527  0.01612 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c gd.z.c
## gndr.c       0.098              
## gdi.z.cm     0.000  0.000       
## gndr.c:gd..  0.000 -0.008  0.097
```

``` r
getFE(mod2_GDI,round=3)
```

```
##                   Est.    SE     df       t     p     LL     UL
## (Intercept)     -0.018 0.037 33.955  -0.471 0.640 -0.094  0.059
## gndr.c          -0.169 0.014 33.404 -11.752 0.000 -0.198 -0.140
## gdi.z.cm        -0.117 0.038 34.032  -3.074 0.004 -0.194 -0.040
## gndr.c:gdi.z.cm -0.037 0.015 35.343  -2.527 0.016 -0.067 -0.007
```

``` r
getVC(mod2_GDI)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.22 0.05
## 2    cntry      gndr.c   <NA>  0.08 0.01
## 3    cntry (Intercept) gndr.c  0.10 0.00
## 4 Residual        <NA>   <NA>  1.00 1.01
```

``` r
r2mlm(mod2_GDI,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.014877014
## slope variation 0.001543776
## mean variation  0.044335949
## sigma2          0.939243261
## 
## $R2s
##           total
## f   0.014877014
## v   0.001543776
## m   0.044335949
## fv  0.016420790
## fvm 0.060756739
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
## Time difference of 37.89665 secs
```

``` r
round(ddsc_mod2_GDI$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.057        0.239        1.007     1.064 0.053   7802.647 0.998   0.998
## 2        0.5         0.069        0.263        1.007     1.076 0.064   6678.029 0.998   0.998
```

``` r
round(ddsc_mod2_GDI$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gdi.z.cm gdi.z.cm_scaled
## means_y1          -0.098 0.274    1.000           1.000    0.940           0.940   -0.486          -0.486
## means_y1_scaled   -0.374 1.052    1.000           1.000    0.940           0.940   -0.486          -0.486
## means_y2           0.071 0.246    0.940           0.940    1.000           1.000   -0.381          -0.381
## means_y2_scaled    0.273 0.945    0.940           0.940    1.000           1.000   -0.381          -0.381
## gdi.z.cm           0.000 1.000   -0.486          -0.486   -0.381          -0.381    1.000           1.000
## gdi.z.cm_scaled    0.000 1.000   -0.486          -0.486   -0.381          -0.381    1.000           1.000
## diff_score        -0.169 0.094    0.453           0.453    0.123           0.123   -0.419          -0.419
## diff_score_scaled -0.647 0.361    0.453           0.453    0.123           0.123   -0.419          -0.419
##                   diff_score diff_score_scaled
## means_y1               0.453             0.453
## means_y1_scaled        0.453             0.453
## means_y2               0.123             0.123
## means_y2_scaled        0.123             0.123
## gdi.z.cm              -0.419            -0.419
## gdi.z.cm_scaled       -0.419            -0.419
## diff_score             1.000             1.000
## diff_score_scaled      1.000             1.000
```

``` r
round(ddsc_mod2_GDI$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                       0.398 0.157 35.343   2.527   0.016    0.078    0.717
## w_11                         -0.098 0.038 34.065  -2.582   0.014   -0.176   -0.021
## w_21                         -0.136 0.039 34.062  -3.438   0.002   -0.216   -0.056
## r_xy1                        -0.358 0.139 34.065  -2.582   0.014   -0.640   -0.076
## r_xy2                        -0.551 0.160 34.062  -3.438   0.002   -0.877   -0.225
## b_11                         -0.378 0.146 34.065  -2.582   0.014   -0.675   -0.080
## b_21                         -0.521 0.152 34.062  -3.438   0.002   -0.830   -0.213
## main_effect                  -0.117 0.038 34.032  -3.074   0.004   -0.194   -0.040
## moderator_effect             -0.169 0.014 33.404 -11.752   0.000   -0.198   -0.140
## interaction                  -0.037 0.015 35.343  -2.527   0.016   -0.067   -0.007
## q_b11_b21                     0.181    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.245    NA     NA      NA      NA       NA       NA
## cross_over_point             -4.515    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.080 0.039 34.178  -2.015   0.052   -0.160    0.001
## interaction_vs_main_bscale   -0.306 0.152 34.178  -2.015   0.052   -0.614    0.002
## interaction_vs_main_rscale   -0.262 0.137 34.206  -1.911   0.064   -0.541    0.017
## dadas                        -0.197 0.076 34.065  -2.582   0.993   -0.351   -0.042
## dadas_bscale                 -0.755 0.292 34.065  -2.582   0.993   -1.349   -0.161
## dadas_rscale                 -0.717 0.278 34.065  -2.582   0.993   -1.281   -0.153
## abs_diff                      0.037 0.015 35.343   2.527   0.008    0.007    0.067
## abs_sum                       0.234 0.076 34.032   3.074   0.002    0.079    0.389
## abs_diff_bscale               0.144 0.057 35.343   2.527   0.008    0.028    0.259
## abs_sum_bscale                0.899 0.292 34.032   3.074   0.002    0.305    1.493
## abs_diff_rscale               0.193 0.061 35.031   3.176   0.002    0.069    0.316
## abs_sum_rscale                0.909 0.294 34.032   3.097   0.002    0.313    1.506
```

``` r
round(ddsc_mod2_GDI$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
##  0.006  0.275  2.521  1.000  0.112
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
##                                    est    se      z pvalue ci.lower ci.upper
## r_xy1_y2                         0.419 0.156  2.694  0.007    0.114    0.725
## r_xy1                           -0.381 0.159 -2.402  0.016   -0.692   -0.070
## r_xy2                           -0.486 0.150 -3.242  0.001   -0.780   -0.192
## b_11                            -0.360 0.150 -2.402  0.016   -0.654   -0.066
## b_21                            -0.511 0.158 -3.242  0.001   -0.820   -0.202
## b_10                             0.273 0.148  1.852  0.064   -0.016    0.563
## b_20                            -0.374 0.155 -2.408  0.016   -0.679   -0.070
## res_cov_y1_y2                    0.729 0.183  3.981  0.000    0.370    1.088
## diff_b10_b20                     0.647 0.055 11.694  0.000    0.539    0.756
## diff_b11_b21                     0.151 0.056  2.694  0.007    0.041    0.262
## diff_rxy1_rxy2                   0.105 0.056  1.861  0.063   -0.006    0.216
## q_b11_b21                        0.188 0.081  2.323  0.020    0.029    0.346
## q_rxy1_rxy2                      0.130 0.070  1.856  0.064   -0.007    0.267
## cross_over_point                -4.276 1.629 -2.625  0.009   -7.468   -1.084
## sum_b11_b21                     -0.871 0.302 -2.880  0.004   -1.464   -0.278
## main_effect                     -0.436 0.151 -2.880  0.004   -0.732   -0.139
## interaction_vs_main_effect      -0.284 0.154 -1.849  0.064   -0.585    0.017
## diff_abs_b11_abs_b21            -0.151 0.056 -2.694  0.007   -0.262   -0.041
## abs_diff_b11_b21                 0.151 0.056  2.694  0.004    0.041    0.262
## abs_sum_b11_b21                  0.871 0.302  2.880  0.002    0.278    1.464
## dadas                           -0.720 0.300 -2.402  0.992   -1.307   -0.132
## q_r_equivalence                  0.030 0.070  0.425  0.665       NA       NA
## q_b_equivalence                  0.088 0.081  1.086  0.861       NA       NA
## cross_over_point_equivalence     4.276 1.629  2.625  0.996       NA       NA
## cross_over_point_minimal_effect  4.276 1.629  2.625  0.004       NA       NA
```

``` r
round(ddsc_sem_GDI$variance_test,3)
```

```
##              est    se      z pvalue ci.lower ci.upper
## cov_y1y2   0.907 0.227  3.994  0.000    0.462    1.353
## var_y1     0.867 0.210  4.123  0.000    0.455    1.279
## var_y2     1.075 0.261  4.123  0.000    0.564    1.585
## var_diff  -0.208 0.123 -1.683  0.092   -0.450    0.034
## var_ratio  0.807 0.094  8.563  0.000    0.622    0.991
## cor_y1y2   0.940 0.020 47.290  0.000    0.901    0.979
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
## r_xy1y2                       0.419 0.160 32.000   2.614   0.014    0.093    0.746
## w_11                         -0.094 0.041 34.207  -2.270   0.030   -0.178   -0.010
## w_21                         -0.133 0.041 34.207  -3.225   0.003   -0.217   -0.049
## r_xy1                        -0.381 0.168 34.207  -2.270   0.030   -0.722   -0.040
## r_xy2                        -0.486 0.151 34.207  -3.225   0.003   -0.792   -0.180
## b_11                         -0.360 0.159 34.207  -2.270   0.030   -0.683   -0.038
## b_21                         -0.512 0.159 34.207  -3.225   0.003   -0.835   -0.189
## main_effect                  -0.114 0.041 32.000  -2.794   0.009   -0.196   -0.031
## moderator_effect             -0.169 0.015 32.000 -11.344   0.000   -0.199   -0.138
## interaction                  -0.039 0.015 32.000  -2.614   0.014   -0.070   -0.009
## q_b11_b21                     0.188    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.130    NA     NA      NA      NA       NA       NA
## cross_over_point             -4.276    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.074 0.043 40.672  -1.709   0.095   -0.162    0.013
## interaction_vs_main_bscale   -0.285 0.167 40.672  -1.709   0.095   -0.621    0.052
## interaction_vs_main_rscale   -0.328 0.183 39.493  -1.790   0.081   -0.699    0.042
## dadas                        -0.188 0.083 34.207  -2.270   0.985   -0.356   -0.020
## dadas_bscale                 -0.721 0.318 34.207  -2.270   0.985   -1.366   -0.076
## dadas_rscale                 -0.762 0.336 34.207  -2.270   0.985   -1.444   -0.080
## abs_diff                      0.039 0.015 32.000   2.614   0.007    0.009    0.070
## abs_sum                       0.227 0.081 32.000   2.794   0.004    0.062    0.393
## abs_diff_bscale               0.152 0.058 32.000   2.614   0.007    0.033    0.270
## abs_sum_bscale                0.872 0.312 32.000   2.794   0.004    0.236    1.508
## abs_diff_rscale               0.105 0.061 37.306   1.735   0.045   -0.018    0.228
## abs_sum_rscale                0.867 0.313 32.006   2.768   0.005    0.229    1.505
```

``` r
# country-time multilevel model


mod2_GDI_cntry_year<-
  lmer(uni.z.wt~gndr.c+gdi.z.cm+gndr.c:gdi.z.cm+
                 (gndr.c|cntry),data=diff_dat_cntry_year,REML=F,
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2_GDI_cntry_year)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: uni.z.wt ~ gndr.c + gdi.z.cm + gndr.c:gdi.z.cm + (gndr.c | cntry)
##    Data: diff_dat_cntry_year
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##    -298.6    -264.2     157.3    -314.6       538 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.1034 -0.4604  0.0789  0.5731  3.0369 
## 
## Random effects:
##  Groups   Name        Variance  Std.Dev. Corr 
##  cntry    (Intercept) 0.0415598 0.20386       
##           gndr.c      0.0002069 0.01438  0.36 
##  Residual             0.0269796 0.16425       
## Number of obs: 546, groups:  cntry, 34
## 
## Fixed effects:
##                 Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)     -0.01001    0.03583 33.59518  -0.279  0.78160    
## gndr.c          -0.17260    0.01432 33.45772 -12.053 1.01e-13 ***
## gdi.z.cm        -0.12256    0.03699 35.85410  -3.313  0.00212 ** 
## gndr.c:gdi.z.cm -0.02414    0.01750 48.44819  -1.379  0.17426    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c gd.z.c
## gndr.c       0.062              
## gdi.z.cm    -0.009  0.000       
## gndr.c:gd..  0.000 -0.056  0.051
```

``` r
getFE(mod2_GDI_cntry_year,round=3)
```

```
##                   Est.    SE     df       t     p     LL     UL
## (Intercept)     -0.010 0.036 33.595  -0.279 0.782 -0.083  0.063
## gndr.c          -0.173 0.014 33.458 -12.053 0.000 -0.202 -0.143
## gdi.z.cm        -0.123 0.037 35.854  -3.313 0.002 -0.198 -0.048
## gndr.c:gdi.z.cm -0.024 0.018 48.448  -1.379 0.174 -0.059  0.011
```

``` r
getVC(mod2_GDI_cntry_year)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.20 0.04
## 2    cntry      gndr.c   <NA>  0.01 0.00
## 3    cntry (Intercept) gndr.c  0.36 0.00
## 4 Residual        <NA>   <NA>  0.16 0.03
```

``` r
r2mlm(mod2_GDI,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.014877014
## slope variation 0.001543776
## mean variation  0.044335949
## sigma2          0.939243261
## 
## $R2s
##           total
## f   0.014877014
## v   0.001543776
## m   0.044335949
## fv  0.016420790
## fvm 0.060756739
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
## 1       -0.5         0.052        0.229        0.027     0.079 0.660      8.029 0.995   0.940
## 2        0.5         0.059        0.244        0.027     0.086 0.687      8.029 0.995   0.946
```

``` r
round(ddsc_mod2_GDI_cntry_year$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gdi.z.cm gdi.z.cm_scaled
## means_y1          -0.097 0.261    1.000           1.000    0.937           0.937   -0.531          -0.531
## means_y1_scaled   -0.391 1.047    1.000           1.000    0.937           0.937   -0.531          -0.531
## means_y2           0.074 0.237    0.937           0.937    1.000           1.000   -0.433          -0.433
## means_y2_scaled    0.297 0.950    0.937           0.937    1.000           1.000   -0.433          -0.433
## gdi.z.cm           0.000 1.000   -0.531          -0.531   -0.433          -0.433    1.000           1.000
## gdi.z.cm_scaled    0.000 1.000   -0.531          -0.531   -0.433          -0.433    1.000           1.000
## diff_score        -0.171 0.092    0.426           0.426    0.083           0.083   -0.392          -0.392
## diff_score_scaled -0.688 0.368    0.426           0.426    0.083           0.083   -0.392          -0.392
##                   diff_score diff_score_scaled
## means_y1               0.426             0.426
## means_y1_scaled        0.426             0.426
## means_y2               0.083             0.083
## means_y2_scaled        0.083             0.083
## gdi.z.cm              -0.392            -0.392
## gdi.z.cm_scaled       -0.392            -0.392
## diff_score             1.000             1.000
## diff_score_scaled      1.000             1.000
```

``` r
round(ddsc_mod2_GDI_cntry_year$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                       0.263 0.191 48.448   1.379   0.174   -0.120    0.646
## w_11                         -0.110 0.038 36.764  -2.941   0.006   -0.187   -0.034
## w_21                         -0.135 0.038 36.503  -3.501   0.001   -0.213   -0.057
## r_xy1                        -0.423 0.144 36.764  -2.941   0.006   -0.715   -0.132
## r_xy2                        -0.568 0.162 36.503  -3.501   0.001   -0.897   -0.239
## b_11                         -0.444 0.151 36.764  -2.941   0.006   -0.749   -0.138
## b_21                         -0.541 0.154 36.503  -3.501   0.001   -0.853   -0.228
## main_effect                  -0.123 0.037 35.854  -3.313   0.002   -0.198   -0.048
## moderator_effect             -0.173 0.014 33.458 -12.053   0.000   -0.202   -0.143
## interaction                  -0.024 0.018 48.448  -1.379   0.174   -0.059    0.011
## q_b11_b21                     0.128    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.193    NA     NA      NA      NA       NA       NA
## cross_over_point             -7.152    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.098 0.040 39.142  -2.454   0.019   -0.180   -0.017
## interaction_vs_main_bscale   -0.395 0.161 39.142  -2.454   0.019   -0.721   -0.070
## interaction_vs_main_rscale   -0.351 0.148 39.609  -2.372   0.023   -0.649   -0.052
## dadas                        -0.221 0.075 36.764  -2.941   0.997   -0.373   -0.069
## dadas_bscale                 -0.887 0.302 36.764  -2.941   0.997   -1.499   -0.276
## dadas_rscale                 -0.846 0.288 36.764  -2.941   0.997   -1.429   -0.263
## abs_diff                      0.024 0.018 48.448   1.379   0.087   -0.011    0.059
## abs_sum                       0.245 0.074 35.854   3.313   0.001    0.095    0.395
## abs_diff_bscale               0.097 0.070 48.448   1.379   0.087   -0.044    0.238
## abs_sum_bscale                0.984 0.297 35.854   3.313   0.001    0.382    1.587
## abs_diff_rscale               0.145 0.073 46.576   1.995   0.026   -0.001    0.291
## abs_sum_rscale                0.991 0.298 35.853   3.327   0.001    0.387    1.595
```

``` r
round(ddsc_mod2_GDI_cntry_year$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
##  0.003  0.649  0.901  1.000  0.343
```

### Comparisons across models


``` r
# full multilevel model (individuals as entries)
round(ddsc_mod2_GDI$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df  t.ratio p.value ci.lower ci.upper
## r_xy1             -0.3583 0.1388 34.0647  -2.5821  0.0143  -0.6404  -0.0763
## r_xy2             -0.5509 0.1602 34.0619  -3.4381  0.0016  -0.8766  -0.2253
## b_11              -0.3776 0.1462 34.0647  -2.5821  0.0143  -0.6747  -0.0804
## b_21              -0.5214 0.1516 34.0619  -3.4381  0.0016  -0.8295  -0.2132
## main_effect       -0.1170 0.0381 34.0320  -3.0740  0.0041  -0.1944  -0.0397
## moderator_effect  -0.1690 0.0144 33.4037 -11.7524  0.0000  -0.1983  -0.1398
## interaction       -0.0374 0.0148 35.3427  -2.5271  0.0161  -0.0675  -0.0074
## q_b11_b21          0.1810     NA      NA       NA      NA       NA       NA
```

``` r
# country-level structural path model (country-gender means as entries)
round(ddsc_sem_GDI$results[c("r_xy1","r_xy2","b_11","b_21","q_b11_b21","diff_b11_b21"),],4)
```

```
##                  est     se       z pvalue ci.lower ci.upper
## r_xy1        -0.3808 0.1586 -2.4017 0.0163  -0.6916  -0.0700
## r_xy2        -0.4859 0.1499 -3.2420 0.0012  -0.7797  -0.1922
## b_11         -0.3599 0.1498 -2.4017 0.0163  -0.6536  -0.0662
## b_21         -0.5113 0.1577 -3.2420 0.0012  -0.8204  -0.2022
## q_b11_b21     0.1877 0.0808  2.3233 0.0202   0.0294   0.3461
## diff_b11_b21  0.1514 0.0562  2.6942 0.0071   0.0413   0.2616
```

``` r
# multilevel model (country-gender means as entries)
round(ddsc_mod2_GDI_ri$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df  t.ratio p.value ci.lower ci.upper
## r_xy1             -0.3808 0.1678 34.2067  -2.2697  0.0296  -0.7218  -0.0399
## r_xy2             -0.4859 0.1507 34.2067  -3.2246  0.0028  -0.7921  -0.1798
## b_11              -0.3604 0.1588 34.2067  -2.2697  0.0296  -0.6830  -0.0378
## b_21              -0.5120 0.1588 34.2067  -3.2246  0.0028  -0.8347  -0.1894
## main_effect       -0.1136 0.0406 32.0000  -2.7941  0.0087  -0.1964  -0.0308
## moderator_effect  -0.1688 0.0149 32.0000 -11.3444  0.0000  -0.1991  -0.1385
## interaction       -0.0395 0.0151 32.0000  -2.6138  0.0135  -0.0702  -0.0087
## q_b11_b21          0.1881     NA      NA       NA      NA       NA       NA
```

``` r
# multilevel model (country-year-gender means as entries)
round(ddsc_mod2_GDI_cntry_year$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df  t.ratio p.value ci.lower ci.upper
## r_xy1             -0.4231 0.1439 36.7636  -2.9408  0.0056  -0.7147  -0.1315
## r_xy2             -0.5680 0.1622 36.5030  -3.5015  0.0012  -0.8969  -0.2392
## b_11              -0.4436 0.1508 36.7636  -2.9408  0.0056  -0.7493  -0.1379
## b_21              -0.5405 0.1544 36.5030  -3.5015  0.0012  -0.8534  -0.2276
## main_effect       -0.1226 0.0370 35.8541  -3.3131  0.0021  -0.1976  -0.0475
## moderator_effect  -0.1726 0.0143 33.4577 -12.0527  0.0000  -0.2017  -0.1435
## interaction       -0.0241 0.0175 48.4482  -1.3789  0.1743  -0.0593   0.0110
## q_b11_b21          0.1282     NA      NA       NA      NA       NA       NA
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
## Time difference of 1.325531 hours
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
```

```
## Warning in atanh(mod2_GDI_boot_est$r_xy2): NaNs produced
```

``` r
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
##                    Estimate         SE        2.5%        97.5%
## X.Intercept.    -0.01859990 0.03780094 -0.08829167  0.056660198
## gndr.c          -0.16861557 0.01486707 -0.19832587 -0.140064569
## gdi.z.cm        -0.11533505 0.03945175 -0.19512008 -0.041496365
## gndr.c.gdi.z.cm -0.03784031 0.01530182 -0.06787508 -0.007203178
## w11             -0.09641490 0.03947878 -0.17827119 -0.022828786
## w21             -0.13425520 0.04088250 -0.21793564 -0.060004768
## b11             -0.37034235 0.15164322 -0.68476317 -0.087688384
## b21             -0.51569197 0.15703509 -0.83711959 -0.230486240
## r_xy1           -0.35147326 0.14391694 -0.64987422 -0.083220626
## r_xy2           -0.54494788 0.16594390 -0.88461053 -0.243562039
## q_b              0.20006713 0.11578415  0.03600526  0.435452426
## q                0.27421219 0.16685944  0.08491397  0.622651863
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
## [1] 0.2000671
## 
## $se
## [1] 0.1157841
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
## [1] 2.591608
## 
## $p_low
## [1] 0.004776427
## 
## $z_high
## [1] 0.8642558
## 
## $p_high
## [1] 0.8062763
## 
## $ci_level
## [1] 0.9
## 
## $ci_lower
## [1] 0.009619151
## 
## $ci_upper
## [1] 0.3905151
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
## [1] 0.2742122
## 
## $se
## [1] 0.1668594
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
## [1] 2.242679
## 
## $p_low
## [1] 0.01245876
## 
## $z_high
## [1] 1.044066
## 
## $p_high
## [1] 0.8517725
## 
## $ci_level
## [1] 0.9
## 
## $ci_lower
## [1] -0.0002471665
## 
## $ci_upper
## [1] 0.5486716
## 
## $equivalent
## [1] FALSE
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
  xlab("Gender Development Index (Average 2002-2023)")+
  #ylim(c(min.y.pred,max.y.pred))+
  #xlim(c(0.60,1.00))+
  ylab("Value universalism (Average 2002-2023)")+
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
  xlab("Gender Development Index (Average 2002-2023)")+
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
```

```
## Warning: The `size` argument of `element_line()` is deprecated as of ggplot2 3.4.0.
## ℹ Please use the `linewidth` argument instead.
## This warning is displayed once per session.
## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.
```

``` r
#p2.uni.flags


pflag_comb<-
  ggarrange(p1.uni.flags,p2.uni.flags,align = "v",
            ncol=1,nrow=2,heights=c(2,1))
```

```
## Warning in geom_text(inherit.aes = F, aes(x = 0.9, y = -0.45, label = coef_q, : All aesthetics have length 1, but the data has 262 rows.
## ℹ Please consider using `annotate()` or provide this layer with data containing a single row.
```

``` r
pflag_comb
```

![](Analysis_universalism_files/figure-html/unnamed-chunk-35-1.png)<!-- -->

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
##   1464420   1464509   -732202   1464404    492335 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -8.5011 -0.5502  0.1017  0.6560  4.8561 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr 
##  cntry    (Intercept) 0.053449 0.23119       
##           gndr.c      0.007872 0.08872  0.33 
##  Residual             1.007158 1.00357       
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##                      Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)         -0.015708   0.039694 33.964945  -0.396   0.6948    
## gndr.c              -0.169193   0.015570 32.705398 -10.867 2.17e-12 ***
## log_gdp.z.cm         0.086897   0.039817 34.002001   2.182   0.0361 *  
## gndr.c:log_gdp.z.cm -0.007484   0.015718 33.536715  -0.476   0.6371    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c lg_g..
## gndr.c      0.319               
## lg_gdp.z.cm 0.023  0.007        
## gndr.c:l_.. 0.007  0.008  0.317
```

``` r
getFE(mod2_log_GDP,round=3)
```

```
##                       Est.    SE     df       t     p     LL     UL
## (Intercept)         -0.016 0.040 33.965  -0.396 0.695 -0.096  0.065
## gndr.c              -0.169 0.016 32.705 -10.867 0.000 -0.201 -0.138
## log_gdp.z.cm         0.087 0.040 34.002   2.182 0.036  0.006  0.168
## gndr.c:log_gdp.z.cm -0.007 0.016 33.537  -0.476 0.637 -0.039  0.024
```

``` r
getVC(mod2_log_GDP)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.23 0.05
## 2    cntry      gndr.c   <NA>  0.09 0.01
## 3    cntry (Intercept) gndr.c  0.33 0.01
## 4 Residual        <NA>   <NA>  1.00 1.01
```

``` r
r2mlm(mod2_log_GDP,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.011992088
## slope variation 0.001819673
## mean variation  0.049249695
## sigma2          0.936938544
## 
## $R2s
##           total
## f   0.011992088
## v   0.001819673
## m   0.049249695
## fv  0.013811761
## fvm 0.063061456
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
## Time difference of 38.17223 secs
```

``` r
round(ddsc_mod2_log_GDP$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.057        0.239        1.007     1.064 0.053   7802.647 0.998   0.998
## 2        0.5         0.069        0.263        1.007     1.076 0.064   6678.029 0.998   0.998
```

``` r
round(ddsc_mod2_log_GDP$descriptives,3)
```

```
##                          M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled log_gdp.z.cm
## means_y1            -0.098 0.274    1.000           1.000    0.940           0.940        0.330
## means_y1_scaled     -0.374 1.052    1.000           1.000    0.940           0.940        0.330
## means_y2             0.071 0.246    0.940           0.940    1.000           1.000        0.388
## means_y2_scaled      0.273 0.945    0.940           0.940    1.000           1.000        0.388
## log_gdp.z.cm        -0.024 1.012    0.330           0.330    0.388           0.388        1.000
## log_gdp.z.cm_scaled  0.000 1.000    0.330           0.330    0.388           0.388        1.000
## diff_score          -0.169 0.094    0.453           0.453    0.123           0.123       -0.054
## diff_score_scaled   -0.647 0.361    0.453           0.453    0.123           0.123       -0.054
##                     log_gdp.z.cm_scaled diff_score diff_score_scaled
## means_y1                          0.330      0.453             0.453
## means_y1_scaled                   0.330      0.453             0.453
## means_y2                          0.388      0.123             0.123
## means_y2_scaled                   0.388      0.123             0.123
## log_gdp.z.cm                      1.000     -0.054            -0.054
## log_gdp.z.cm_scaled               1.000     -0.054            -0.054
## diff_score                       -0.054      1.000             1.000
## diff_score_scaled                -0.054      1.000             1.000
```

``` r
round(ddsc_mod2_log_GDP$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                       0.080 0.167 33.537   0.476   0.637   -0.260    0.419
## w_11                          0.091 0.038 34.012   2.381   0.023    0.013    0.168
## w_21                          0.083 0.043 33.992   1.936   0.061   -0.004    0.170
## r_xy1                         0.330 0.139 34.012   2.381   0.023    0.048    0.612
## r_xy2                         0.338 0.174 33.992   1.936   0.061   -0.017    0.692
## b_11                          0.348 0.146 34.012   2.381   0.023    0.051    0.645
## b_21                          0.319 0.165 33.992   1.936   0.061   -0.016    0.655
## main_effect                   0.087 0.040 34.002   2.182   0.036    0.006    0.168
## moderator_effect             -0.169 0.016 32.705 -10.867   0.000   -0.201   -0.138
## interaction                  -0.007 0.016 33.537  -0.476   0.637   -0.039    0.024
## q_b11_b21                     0.032    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                  -0.008    NA     NA      NA      NA       NA       NA
## cross_over_point            -22.608    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.079 0.047 33.966  -1.682   0.102   -0.175    0.017
## interaction_vs_main_bscale   -0.305 0.181 33.966  -1.682   0.102   -0.674    0.064
## interaction_vs_main_rscale   -0.341 0.199 33.971  -1.718   0.095   -0.745    0.062
## dadas                        -0.166 0.086 33.992  -1.936   0.969   -0.341    0.008
## dadas_bscale                 -0.639 0.330 33.992  -1.936   0.969   -1.310    0.032
## dadas_rscale                 -0.675 0.349 33.992  -1.936   0.969   -1.384    0.034
## abs_diff                      0.007 0.016 33.537   0.476   0.319   -0.024    0.039
## abs_sum                       0.174 0.080 34.002   2.182   0.018    0.012    0.336
## abs_diff_bscale               0.029 0.060 33.537   0.476   0.319   -0.094    0.152
## abs_sum_bscale                0.668 0.306 34.002   2.182   0.018    0.046    1.289
## abs_diff_rscale              -0.007 0.068 33.578  -0.105   0.542   -0.145    0.130
## abs_sum_rscale                0.668 0.308 34.002   2.170   0.019    0.042    1.294
```

``` r
round(ddsc_mod2_log_GDP$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
##  0.006  0.275  2.521  1.000  0.112
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
##                                     est      se      z pvalue ci.lower ci.upper
## r_xy1_y2                          0.054   0.171  0.317  0.751   -0.281    0.390
## r_xy1                             0.388   0.158  2.456  0.014    0.078    0.698
## r_xy2                             0.330   0.162  2.038  0.042    0.013    0.647
## b_11                              0.367   0.149  2.456  0.014    0.074    0.660
## b_21                              0.347   0.170  2.038  0.042    0.013    0.681
## b_10                              0.273   0.147  1.858  0.063   -0.015    0.562
## b_20                             -0.374   0.168 -2.229  0.026   -0.703   -0.045
## res_cov_y1_y2                     0.784   0.197  3.979  0.000    0.398    1.170
## diff_b10_b20                      0.647   0.061 10.631  0.000    0.528    0.767
## diff_b11_b21                      0.020   0.062  0.317  0.751   -0.102    0.141
## diff_rxy1_rxy2                    0.058   0.058  0.996  0.319   -0.056    0.173
## q_b11_b21                         0.022   0.070  0.322  0.748   -0.115    0.160
## q_rxy1_rxy2                       0.067   0.067  0.995  0.320   -0.065    0.199
## cross_over_point                -32.991 103.969 -0.317  0.751 -236.768  170.785
## sum_b11_b21                       0.714   0.314  2.271  0.023    0.098    1.330
## main_effect                       0.357   0.157  2.271  0.023    0.049    0.665
## interaction_vs_main_effect       -0.337   0.188 -1.797  0.072   -0.705    0.031
## diff_abs_b11_abs_b21              0.020   0.062  0.317  0.751   -0.102    0.141
## abs_diff_b11_b21                  0.020   0.062  0.317  0.375   -0.102    0.141
## abs_sum_b11_b21                   0.714   0.314  2.271  0.012    0.098    1.330
## dadas                            -0.694   0.341 -2.038  0.979   -1.362   -0.027
## q_r_equivalence                  -0.033   0.067 -0.493  0.311       NA       NA
## q_b_equivalence                  -0.078   0.070 -1.109  0.134       NA       NA
## cross_over_point_equivalence     32.991 103.969  0.317  0.624       NA       NA
## cross_over_point_minimal_effect  32.991 103.969  0.317  0.376       NA       NA
```

``` r
round(ddsc_sem_log_GDP$variance_test,3)
```

```
##              est    se      z pvalue ci.lower ci.upper
## cov_y1y2   0.907 0.227  3.994  0.000    0.462    1.353
## var_y1     0.867 0.210  4.123  0.000    0.455    1.279
## var_y2     1.075 0.261  4.123  0.000    0.564    1.585
## var_diff  -0.208 0.123 -1.683  0.092   -0.450    0.034
## var_ratio  0.807 0.094  8.563  0.000    0.622    0.991
## cor_y1y2   0.940 0.020 47.290  0.000    0.901    0.979
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
## r_xy1y2                       0.054 0.177 32.000   0.308   0.760   -0.305    0.414
## w_11                          0.096 0.043 34.471   2.221   0.033    0.008    0.183
## w_21                          0.091 0.043 34.471   2.102   0.043    0.003    0.178
## r_xy1                         0.388 0.175 34.471   2.221   0.033    0.033    0.743
## r_xy2                         0.330 0.157 34.471   2.102   0.043    0.011    0.649
## b_11                          0.367 0.165 34.471   2.221   0.033    0.031    0.703
## b_21                          0.348 0.165 34.471   2.102   0.043    0.012    0.684
## main_effect                   0.093 0.042 32.000   2.203   0.035    0.007    0.179
## moderator_effect             -0.169 0.016 32.000 -10.313   0.000   -0.202   -0.135
## interaction                  -0.005 0.017 32.000  -0.308   0.760   -0.039    0.029
## q_b11_b21                     0.023    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.067    NA     NA      NA      NA       NA       NA
## cross_over_point            -32.991    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.088 0.045 41.668  -1.938   0.059   -0.180    0.004
## interaction_vs_main_bscale   -0.338 0.174 41.668  -1.938   0.059   -0.690    0.014
## interaction_vs_main_rscale   -0.301 0.158 43.381  -1.904   0.064   -0.619    0.018
## dadas                        -0.181 0.086 34.471  -2.102   0.979   -0.356   -0.006
## dadas_bscale                 -0.695 0.331 34.471  -2.102   0.979   -1.367   -0.024
## dadas_rscale                 -0.660 0.314 34.471  -2.102   0.979   -1.297   -0.022
## abs_diff                      0.005 0.017 32.000   0.308   0.380   -0.029    0.039
## abs_sum                       0.186 0.084 32.000   2.203   0.017    0.014    0.358
## abs_diff_bscale               0.020 0.064 32.000   0.308   0.380   -0.110    0.150
## abs_sum_bscale                0.715 0.325 32.000   2.203   0.017    0.054    1.376
## abs_diff_rscale               0.058 0.066 36.744   0.877   0.193   -0.076    0.193
## abs_sum_rscale                0.718 0.325 32.007   2.206   0.017    0.055    1.381
```

``` r
# country-time multilevel model


mod2_log_GDP_cntry_year<-
  lmer(uni.z.wt~gndr.c+log_gdp.z.cm+gndr.c:log_gdp.z.cm+
                 (gndr.c|cntry),data=diff_dat_cntry_year,REML=F,
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
```

```
## boundary (singular) fit: see help('isSingular')
```

``` r
summary(mod2_log_GDP_cntry_year)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: uni.z.wt ~ gndr.c + log_gdp.z.cm + gndr.c:log_gdp.z.cm + (gndr.c |      cntry)
##    Data: diff_dat_cntry_year
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##    -294.6    -260.2     155.3    -310.6       538 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -5.9143 -0.4814  0.0729  0.5733  3.1958 
## 
## Random effects:
##  Groups   Name        Variance  Std.Dev. Corr 
##  cntry    (Intercept) 0.0483078 0.21979       
##           gndr.c      0.0004487 0.02118  1.00 
##  Residual             0.0269738 0.16424       
## Number of obs: 546, groups:  cntry, 34
## 
## Fixed effects:
##                      Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)          -0.01059    0.03850  33.44535  -0.275    0.785    
## gndr.c               -0.16976    0.01487 182.27726 -11.420   <2e-16 ***
## log_gdp.z.cm          0.08634    0.03886  34.26225   2.222    0.033 *  
## gndr.c:log_gdp.z.cm  -0.01619    0.01588 211.99737  -1.020    0.309    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c lg_g..
## gndr.c       0.241              
## lg_gdp.z.cm  0.007  0.004       
## gndr.c:l_..  0.004 -0.207  0.227
## optimizer (bobyqa) convergence code: 0 (OK)
## boundary (singular) fit: see help('isSingular')
```

``` r
getFE(mod2_log_GDP_cntry_year,round=3)
```

```
##                       Est.    SE      df       t     p     LL     UL
## (Intercept)         -0.011 0.039  33.445  -0.275 0.785 -0.089  0.068
## gndr.c              -0.170 0.015 182.277 -11.420 0.000 -0.199 -0.140
## log_gdp.z.cm         0.086 0.039  34.262   2.222 0.033  0.007  0.165
## gndr.c:log_gdp.z.cm -0.016 0.016 211.997  -1.020 0.309 -0.047  0.015
```

``` r
getVC(mod2_log_GDP_cntry_year)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.22 0.05
## 2    cntry      gndr.c   <NA>  0.02 0.00
## 3    cntry (Intercept) gndr.c  1.00 0.00
## 4 Residual        <NA>   <NA>  0.16 0.03
```

``` r
r2mlm(mod2_log_GDP,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.011992088
## slope variation 0.001819673
## mean variation  0.049249695
## sigma2          0.936938544
## 
## $R2s
##           total
## f   0.011992088
## v   0.001819673
## m   0.049249695
## fv  0.013811761
## fvm 0.063061456
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
## 1       -0.5         0.052        0.229        0.027     0.079 0.660      8.029 0.995   0.940
## 2        0.5         0.059        0.244        0.027     0.086 0.687      8.029 0.995   0.946
```

``` r
round(ddsc_mod2_log_GDP_cntry_year$descriptives,3)
```

```
##                          M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled log_gdp.z.cm
## means_y1            -0.097 0.261    1.000           1.000    0.937           0.937        0.310
## means_y1_scaled     -0.391 1.047    1.000           1.000    0.937           0.937        0.310
## means_y2             0.074 0.237    0.937           0.937    1.000           1.000        0.382
## means_y2_scaled      0.297 0.950    0.937           0.937    1.000           1.000        0.382
## log_gdp.z.cm        -0.024 1.012    0.310           0.310    0.382           0.382        1.000
## log_gdp.z.cm_scaled  0.000 1.000    0.310           0.310    0.382           0.382        1.000
## diff_score          -0.171 0.092    0.426           0.426    0.083           0.083       -0.105
## diff_score_scaled   -0.688 0.368    0.426           0.426    0.083           0.083       -0.105
##                     log_gdp.z.cm_scaled diff_score diff_score_scaled
## means_y1                          0.310      0.426             0.426
## means_y1_scaled                   0.310      0.426             0.426
## means_y2                          0.382      0.083             0.083
## means_y2_scaled                   0.382      0.083             0.083
## log_gdp.z.cm                      1.000     -0.105            -0.105
## log_gdp.z.cm_scaled               1.000     -0.105            -0.105
## diff_score                       -0.105      1.000             1.000
## diff_score_scaled                -0.105      1.000             1.000
```

``` r
round(ddsc_mod2_log_GDP_cntry_year$results,3)
```

```
##                            estimate    SE      df t.ratio p.value ci.lower ci.upper
## r_xy1y2                       0.176 0.173 211.997   1.020   0.309   -0.165    0.517
## w_11                          0.094 0.038  34.595   2.494   0.018    0.018    0.171
## w_21                          0.078 0.041  34.337   1.891   0.067   -0.006    0.162
## r_xy1                         0.362 0.145  34.595   2.494   0.018    0.067    0.656
## r_xy2                         0.330 0.175  34.337   1.891   0.067   -0.025    0.685
## b_11                          0.379 0.152  34.595   2.494   0.018    0.070    0.688
## b_21                          0.314 0.166  34.337   1.891   0.067   -0.023    0.652
## main_effect                   0.086 0.039  34.262   2.222   0.033    0.007    0.165
## moderator_effect             -0.170 0.015 182.277 -11.420   0.000   -0.199   -0.140
## interaction                  -0.016 0.016 211.997  -1.020   0.309   -0.047    0.015
## q_b11_b21                     0.074    NA      NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.036    NA      NA      NA      NA       NA       NA
## cross_over_point            -10.484    NA      NA      NA      NA       NA       NA
## interaction_vs_main          -0.070 0.045  34.942  -1.552   0.130   -0.162    0.022
## interaction_vs_main_bscale   -0.282 0.181  34.942  -1.552   0.130   -0.650    0.087
## interaction_vs_main_rscale   -0.314 0.197  34.822  -1.596   0.119   -0.714    0.086
## dadas                        -0.156 0.083  34.337  -1.891   0.966   -0.325    0.012
## dadas_bscale                 -0.628 0.332  34.337  -1.891   0.966   -1.303    0.047
## dadas_rscale                 -0.660 0.349  34.337  -1.891   0.966   -1.370    0.049
## abs_diff                      0.016 0.016 211.997   1.020   0.155   -0.015    0.047
## abs_sum                       0.173 0.078  34.262   2.222   0.016    0.015    0.331
## abs_diff_bscale               0.065 0.064 211.997   1.020   0.155   -0.061    0.191
## abs_sum_bscale                0.693 0.312  34.262   2.222   0.016    0.059    1.327
## abs_diff_rscale               0.032 0.069  91.752   0.457   0.324   -0.105    0.168
## abs_sum_rscale                0.692 0.313  34.258   2.207   0.017    0.055    1.329
```

``` r
round(ddsc_mod2_log_GDP_cntry_year$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
##  0.003  0.649  0.901  1.000  0.343
```

### Comparisons across models


``` r
# full multilevel model (individuals as entries)
round(ddsc_mod2_log_GDP$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df  t.ratio p.value ci.lower ci.upper
## r_xy1              0.3304 0.1388 34.0119   2.3813  0.0230   0.0484   0.6124
## r_xy2              0.3375 0.1744 33.9923   1.9356  0.0613  -0.0168   0.6919
## b_11               0.3482 0.1462 34.0119   2.3813  0.0230   0.0510   0.6453
## b_21               0.3194 0.1650 33.9923   1.9356  0.0613  -0.0159   0.6548
## main_effect        0.0869 0.0398 34.0020   2.1824  0.0361   0.0060   0.1678
## moderator_effect  -0.1692 0.0156 32.7054 -10.8667  0.0000  -0.2009  -0.1375
## interaction       -0.0075 0.0157 33.5367  -0.4761  0.6371  -0.0394   0.0245
## q_b11_b21          0.0324     NA      NA       NA      NA       NA       NA
```

``` r
# country-level structural path model (country-gender means as entries)
round(ddsc_sem_log_GDP$results[c("r_xy1","r_xy2","b_11","b_21","q_b11_b21","diff_b11_b21"),],4)
```

```
##                 est     se      z pvalue ci.lower ci.upper
## r_xy1        0.3881 0.1581 2.4558 0.0141   0.0784   0.6979
## r_xy2        0.3299 0.1619 2.0380 0.0415   0.0126   0.6473
## b_11         0.3668 0.1494 2.4558 0.0141   0.0741   0.6595
## b_21         0.3472 0.1703 2.0380 0.0415   0.0133   0.6810
## q_b11_b21    0.0225 0.0699 0.3217 0.7477  -0.1145   0.1595
## diff_b11_b21 0.0196 0.0618 0.3175 0.7509  -0.1015   0.1408
```

``` r
# multilevel model (country-gender means as entries)
round(ddsc_mod2_log_GDP_ri$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df  t.ratio p.value ci.lower ci.upper
## r_xy1              0.3881 0.1747 34.4711   2.2213  0.0330   0.0332   0.7431
## r_xy2              0.3299 0.1569 34.4711   2.1024  0.0429   0.0112   0.6487
## b_11               0.3673 0.1654 34.4711   2.2213  0.0330   0.0314   0.7032
## b_21               0.3477 0.1654 34.4711   2.1024  0.0429   0.0118   0.6835
## main_effect        0.0931 0.0422 32.0000   2.2033  0.0349   0.0070   0.1791
## moderator_effect  -0.1688 0.0164 32.0000 -10.3135  0.0000  -0.2021  -0.1355
## interaction       -0.0051 0.0166 32.0000  -0.3080  0.7601  -0.0390   0.0287
## q_b11_b21          0.0225     NA      NA       NA      NA       NA       NA
```

``` r
# multilevel model (country-year-gender means as entries)
round(ddsc_mod2_log_GDP_cntry_year$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE       df  t.ratio p.value ci.lower ci.upper
## r_xy1              0.3616 0.1450  34.5949   2.4944  0.0175   0.0672   0.6560
## r_xy2              0.3301 0.1746  34.3370   1.8905  0.0671  -0.0246   0.6848
## b_11               0.3791 0.1520  34.5949   2.4944  0.0175   0.0704   0.6878
## b_21               0.3141 0.1662  34.3370   1.8905  0.0671  -0.0234   0.6517
## main_effect        0.0863 0.0389  34.2622   2.2219  0.0330   0.0074   0.1653
## moderator_effect  -0.1698 0.0149 182.2773 -11.4195  0.0000  -0.1991  -0.1404
## interaction       -0.0162 0.0159 211.9974  -1.0196  0.3091  -0.0475   0.0151
## q_b11_b21          0.0739     NA       NA       NA      NA       NA       NA
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
## Time difference of 1.344592 hours
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
##                         Estimate         SE         2.5%       97.5%
## X.Intercept.        -0.016703129 0.04005731 -0.090752343  0.06264347
## gndr.c              -0.168866102 0.01591149 -0.200614946 -0.13768229
## log_gdp.z.cm         0.086557127 0.03883069  0.010607348  0.16174920
## gndr.c.log_gdp.z.cm -0.006754189 0.01517979 -0.035412316  0.02183328
## w11                  0.089934222 0.03727991  0.016322506  0.16416524
## w21                  0.083180033 0.04172609  0.002088207  0.16064852
## b11                  0.345449225 0.14319705  0.062696902  0.63058037
## b21                  0.319505492 0.16027543  0.008021079  0.61707217
## r_xy1                0.327848451 0.13590111  0.059502470  0.59845205
## r_xy2                0.337631474 0.16936807  0.008476126  0.65207952
## q_b                  0.027497980 0.06947694 -0.110856415  0.15414963
## q                   -0.017695574 0.08560794 -0.186447590  0.12438602
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
## [1] 0.02749798
## 
## $se
## [1] 0.06947694
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
## [1] 1.835112
## 
## $p_low
## [1] 0.03324454
## 
## $z_high
## [1] -1.043541
## 
## $p_high
## [1] 0.148349
## 
## $ci_level
## [1] 0.9
## 
## $ci_lower
## [1] -0.08678142
## 
## $ci_upper
## [1] 0.1417774
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
## [1] -0.01769557
## 
## $se
## [1] 0.08560794
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
## [1] 0.9614112
## 
## $p_low
## [1] 0.1681727
## 
## $z_high
## [1] -1.374821
## 
## $p_high
## [1] 0.0845935
## 
## $ci_level
## [1] 0.9
## 
## $ci_lower
## [1] -0.1585081
## 
## $ci_upper
## [1] 0.123117
## 
## $equivalent
## [1] FALSE
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
  xlab("log(GDP per capita) (Average 2002-2023)")+
  #ylim(c(min.y.pred,max.y.pred))+
  #xlim(c(0.60,1.00))+
  ylab("Value universalism (Average 2002-2023)")+
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
  xlab("log(GDP per capita) (Average 2002-2023)")+
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

![](Analysis_universalism_files/figure-html/unnamed-chunk-41-1.png)<!-- -->

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


## mod3: fixed effect of time (year)

* Year coded year-2002 (2002 is a zero-point)


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
## 1464147.7 1464225.4 -732066.8 1464133.7    492336 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -8.4783 -0.5495  0.1039  0.6558  4.9069 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr 
##  cntry    (Intercept) 0.062851 0.25070       
##           gndr.c      0.007986 0.08936  0.29 
##  Residual             1.006592 1.00329       
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##               Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept) -1.993e-02  4.303e-02  3.397e+01  -0.463    0.646    
## gndr.c      -1.695e-01  1.568e-02  3.282e+01 -10.809 2.39e-12 ***
## essround.c   8.077e-03  4.861e-04  4.920e+05  16.617  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##            (Intr) gndr.c
## gndr.c      0.279       
## essround.c -0.003  0.000
```

``` r
getFE(mod3,round=3)
```

```
##               Est.    SE         df       t     p     LL     UL
## (Intercept) -0.020 0.043     33.969  -0.463 0.646 -0.107  0.068
## gndr.c      -0.169 0.016     32.818 -10.809 0.000 -0.201 -0.138
## essround.c   0.008 0.000 492025.542  16.617 0.000  0.007  0.009
```

``` r
getVC(mod3)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.25 0.06
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
## fixed           0.007176200
## slope variation 0.001839637
## mean variation  0.057816741
## sigma2          0.933167421
## 
## $R2s
##           total
## f   0.007176200
## v   0.001839637
## m   0.057816741
## fv  0.009015838
## fvm 0.066832579
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
## mod2    6 1464422 1464488 -732205   1464410                         
## mod3    7 1464148 1464225 -732067   1464134 276.03  1  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## mod4: random effect of time (year)


``` r
mod4<-lmer(uni.z~gndr.c+year.c+(gndr.c+year.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod4)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: uni.z ~ gndr.c + year.c + (gndr.c + year.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1458564.4 1458675.5 -729272.2 1458544.4    492333 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -8.4901 -0.5511  0.1054  0.6527  4.8329 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr        
##  cntry    (Intercept) 0.715164 0.84567              
##           gndr.c      0.008037 0.08965   0.05       
##           year.c      0.002656 0.05154  -0.95 -0.03 
##  Residual             0.994744 0.99737              
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##              Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)  0.097744   0.145178 32.266634   0.673    0.506    
## gndr.c      -0.169906   0.015721 32.830861 -10.807 2.39e-12 ***
## year.c      -0.004828   0.008850 32.688943  -0.545    0.589    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##        (Intr) gndr.c
## gndr.c  0.051       
## year.c -0.952 -0.030
```

``` r
getFE(mod4,round=3)
```

```
##               Est.    SE     df       t     p     LL     UL
## (Intercept)  0.098 0.145 32.267   0.673 0.506 -0.198  0.393
## gndr.c      -0.170 0.016 32.831 -10.807 0.000 -0.202 -0.138
## year.c      -0.005 0.009 32.689  -0.545 0.589 -0.023  0.013
```

``` r
getVC(mod4)
```

```
##        grp        var1   var2 sdcor  vcov
## 1    cntry (Intercept)   <NA>  0.85  0.72
## 2    cntry      gndr.c   <NA>  0.09  0.01
## 3    cntry      year.c   <NA>  0.05  0.00
## 4    cntry (Intercept) gndr.c  0.05  0.00
## 5    cntry (Intercept) year.c -0.95 -0.04
## 6    cntry      gndr.c year.c -0.03  0.00
## 7 Residual        <NA>   <NA>  1.00  0.99
```

``` r
r2mlm(mod4,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.006469517
## slope variation 0.086653347
## mean variation  0.115515910
## sigma2          0.791361225
## 
## $R2s
##           total
## f   0.006469517
## v   0.086653347
## m   0.115515910
## fv  0.093122865
## fvm 0.208638775
```

``` r
anova(mod2,mod3,mod4)
```

```
## Data: diff_dat
## Models:
## mod2: uni.z ~ gndr.c + (gndr.c | cntry)
## mod3: uni.z ~ gndr.c + essround.c + (gndr.c | cntry)
## mod4: uni.z ~ gndr.c + year.c + (gndr.c + year.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L)   Chisq Df Pr(>Chisq)    
## mod2    6 1464422 1464488 -732205   1464410                          
## mod3    7 1464148 1464225 -732067   1464134  276.03  1  < 2.2e-16 ***
## mod4   10 1458564 1458676 -729272   1458544 5589.24  3  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## mod5: fixed interaction between time and gender


``` r
mod5<-lmer(uni.z~gndr.c+year.c+
             gndr.c:year.c+(gndr.c+year.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod5)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: uni.z ~ gndr.c + year.c + gndr.c:year.c + (gndr.c + year.c |      cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1458548.2 1458670.4 -729263.1 1458526.2    492332 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -8.4995 -0.5518  0.1059  0.6527  4.8367 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr        
##  cntry    (Intercept) 0.712041 0.84383              
##           gndr.c      0.007701 0.08775   0.08       
##           year.c      0.002646 0.05144  -0.95 -0.06 
##  Residual             0.994710 0.99735              
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##                 Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)    9.816e-02  1.449e-01  3.245e+01   0.678    0.503    
## gndr.c        -1.485e-01  1.619e-02  4.007e+01  -9.172 2.16e-11 ***
## year.c        -4.868e-03  8.833e-03  3.286e+01  -0.551    0.585    
## gndr.c:year.c -2.009e-03  4.695e-04  2.821e+05  -4.278 1.89e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c year.c
## gndr.c       0.074              
## year.c      -0.952 -0.056       
## gndr.c:yr.c -0.001 -0.308  0.001
```

``` r
getFE(mod5,round=3)
```

```
##                 Est.    SE         df      t     p     LL     UL
## (Intercept)    0.098 0.145     32.447  0.678 0.503 -0.197  0.393
## gndr.c        -0.149 0.016     40.075 -9.172 0.000 -0.181 -0.116
## year.c        -0.005 0.009     32.862 -0.551 0.585 -0.023  0.013
## gndr.c:year.c -0.002 0.000 282125.887 -4.278 0.000 -0.003 -0.001
```

``` r
getVC(mod5)
```

```
##        grp        var1   var2 sdcor  vcov
## 1    cntry (Intercept)   <NA>  0.84  0.71
## 2    cntry      gndr.c   <NA>  0.09  0.01
## 3    cntry      year.c   <NA>  0.05  0.00
## 4    cntry (Intercept) gndr.c  0.08  0.01
## 5    cntry (Intercept) year.c -0.95 -0.04
## 6    cntry      gndr.c year.c -0.06  0.00
## 7 Residual        <NA>   <NA>  1.00  0.99
```

``` r
r2mlm(mod5,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.006434671
## slope variation 0.086337186
## mean variation  0.115178920
## sigma2          0.792049223
## 
## $R2s
##           total
## f   0.006434671
## v   0.086337186
## m   0.115178920
## fv  0.092771857
## fvm 0.207950777
```

``` r
anova(mod4,mod5)
```

```
## Data: diff_dat
## Models:
## mod4: uni.z ~ gndr.c + year.c + (gndr.c + year.c | cntry)
## mod5: uni.z ~ gndr.c + year.c + gndr.c:year.c + (gndr.c + year.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod4   10 1458564 1458676 -729272   1458544                         
## mod5   11 1458548 1458670 -729263   1458526 18.246  1  1.941e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## mod6: random interaction between time and gender


``` r
mod6<-lmer(uni.z~gndr.c+year.c+
             gndr.c:year.c+(gndr.c+year.c+gndr.c:year.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: uni.z ~ gndr.c + year.c + gndr.c:year.c + (gndr.c + year.c +      gndr.c:year.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1458495.2 1458661.8 -729232.6 1458465.2    492328 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -8.4516 -0.5503  0.1040  0.6530  4.8400 
## 
## Random effects:
##  Groups   Name          Variance  Std.Dev. Corr              
##  cntry    (Intercept)   7.123e-01 0.844000                   
##           gndr.c        1.473e-02 0.121356  0.12             
##           year.c        2.647e-03 0.051448 -0.95 -0.12       
##           gndr.c:year.c 2.244e-05 0.004737 -0.08 -0.80  0.11 
##  Residual               9.945e-01 0.997255                   
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##                 Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)    0.0982796  0.1448912 32.4030787   0.678   0.5024    
## gndr.c        -0.1458087  0.0221742 31.4570768  -6.576 2.25e-07 ***
## year.c        -0.0048971  0.0088347 32.8207998  -0.554   0.5831    
## gndr.c:year.c -0.0023556  0.0009872 29.5682559  -2.386   0.0236 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c year.c
## gndr.c       0.112              
## year.c      -0.952 -0.108       
## gndr.c:yr.c -0.069 -0.788  0.092
```

``` r
getFE(mod6,round=3)
```

```
##                 Est.    SE     df      t     p     LL     UL
## (Intercept)    0.098 0.145 32.403  0.678 0.502 -0.197  0.393
## gndr.c        -0.146 0.022 31.457 -6.576 0.000 -0.191 -0.101
## year.c        -0.005 0.009 32.821 -0.554 0.583 -0.023  0.013
## gndr.c:year.c -0.002 0.001 29.568 -2.386 0.024 -0.004  0.000
```

``` r
getVC(mod6)
```

```
##         grp          var1          var2 sdcor  vcov
## 1     cntry   (Intercept)          <NA>  0.84  0.71
## 2     cntry        gndr.c          <NA>  0.12  0.01
## 3     cntry        year.c          <NA>  0.05  0.00
## 4     cntry gndr.c:year.c          <NA>  0.00  0.00
## 5     cntry   (Intercept)        gndr.c  0.12  0.01
## 6     cntry   (Intercept)        year.c -0.95 -0.04
## 7     cntry   (Intercept) gndr.c:year.c -0.08  0.00
## 8     cntry        gndr.c        year.c -0.12  0.00
## 9     cntry        gndr.c gndr.c:year.c -0.80  0.00
## 10    cntry        year.c gndr.c:year.c  0.11  0.00
## 11 Residual          <NA>          <NA>  1.00  0.99
```

``` r
r2mlm(mod6,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.006507014
## slope variation 0.086470030
## mean variation  0.115205498
## sigma2          0.791817458
## 
## $R2s
##           total
## f   0.006507014
## v   0.086470030
## m   0.115205498
## fv  0.092977044
## fvm 0.208182542
```

``` r
anova(mod4,mod5,mod6)
```

```
## Data: diff_dat
## Models:
## mod4: uni.z ~ gndr.c + year.c + (gndr.c + year.c | cntry)
## mod5: uni.z ~ gndr.c + year.c + gndr.c:year.c + (gndr.c + year.c | cntry)
## mod6: uni.z ~ gndr.c + year.c + gndr.c:year.c + (gndr.c + year.c + gndr.c:year.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod4   10 1458564 1458676 -729272   1458544                         
## mod5   11 1458548 1458670 -729263   1458526 18.246  1  1.941e-05 ***
## mod6   15 1458495 1458662 -729233   1458465 60.987  4  1.799e-12 ***
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
##  year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21  0.0931 0.0652 34.0  -0.0393    0.225   1.428  0.1623
##       0  0.1712 0.1440 31.2  -0.1226    0.465   1.188  0.2437
## 
## gndr.c =  0.5:
##  year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21 -0.1022 0.0660 34.1  -0.2364    0.032  -1.548  0.1310
##       0  0.0254 0.1470 31.4  -0.2733    0.324   0.173  0.8636
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
##  year.c21 - year.c0  -0.0781 0.185 31.9   -0.455    0.298  -0.422  0.6755
## 
## gndr.c =  0.5:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0  -0.1276 0.187 32.0   -0.508    0.253  -0.683  0.4995
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
##  gndr.c year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21  0.0931 0.0652 34.0  -0.0393    0.225   1.428  0.1623
##     0.5     21 -0.1022 0.0660 34.1  -0.2364    0.032  -1.548  0.1310
##    -0.5      0  0.1712 0.1440 31.2  -0.1226    0.465   1.188  0.2437
##     0.5      0  0.0254 0.1470 31.4  -0.2733    0.324   0.173  0.8636
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod6,adjust="none",infer=c(T,T))
```

```
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21     0.1953 0.0141 31.1    0.167    0.224  13.899 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)  -0.0781 0.1850 31.9   -0.455    0.298  -0.422  0.6755
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0      0.0677 0.1880 32.4   -0.314    0.450   0.361  0.7205
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)     -0.2734 0.1850 32.2   -0.650    0.103  -1.480  0.1486
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0        -0.1276 0.1870 32.0   -0.508    0.253  -0.683  0.4995
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0       0.1458 0.0222 31.5    0.101    0.191   6.576 <0.0001
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
##  diff_ESS11   -0.195 0.0141 31.1   -0.224   -0.167 -13.899 <0.0001
##  diff_ESS1    -0.146 0.0222 31.5   -0.191   -0.101  -6.576 <0.0001
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
##  diff_ESS11 - diff_ESS1  -0.0495 0.0207 29.6  -0.0918  -0.0071  -2.386  0.0236
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
## [1] -0.04946833
## 
## $se
## [1] 0.02073204
## 
## $df
## [1] 29.56826
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
## [1] 7.260822
## 
## $p_low
## [1] 2.39552e-08
## 
## $t_high
## [1] -12.03298
## 
## $p_high
## [1] 3.209316e-13
## 
## $ci_level
## [1] 0.8
## 
## $ci_lower
## [1] -0.07664485
## 
## $ci_upper
## [1] -0.02229181
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
  scale_x_continuous(limits = c(2001, 2024),
                     breaks = c(seq(2002, 2020, 2),2023))+
  ylab("Mean-level of value universalism")+
  theme(legend.title=element_blank(),
        legend.text  = element_text(size = 16),
        axis.title.x = element_text(size = 18, face = "bold"),   
        axis.title.y = element_text(size = 18, face = "bold"),   
        axis.text.x  = element_text(size = 14),                  
        axis.text.y  = element_text(size = 14))

p_time_trends
```

![](Analysis_universalism_files/figure-html/unnamed-chunk-48-1.png)<!-- -->

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
pred_cntry_dat$uni.z_mean<-predict(mod6,newdata=pred_cntry_dat)

pred_cntry_dat$year=pred_cntry_dat$year.c+2002

pred_cntry_dat$gender<-
  case_when(
    pred_cntry_dat$gndr.c==0.5~"men",
    pred_cntry_dat$gndr.c==-0.5~"women",
  )

range(pred_cntry_dat$uni.z_mean)
```

```
## [1] -1.3839428  0.6420202
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

pdf("../results/uni/country_specific_time_trends.pdf", width = 8, height = 6)

for (ctry in countries) {
  print(
    ggplot(pred_cntry_dat[pred_cntry_dat$cntry == ctry, ],
       aes(x = year, y = uni.z_mean, color = gender)) +
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
    name   = "Mean-level of value universalism",
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

```
## Warning: Removed 2 rows containing non-finite outside the scale range (`stat_smooth()`).
```

```
## Warning: Removed 2 rows containing missing values or values outside the scale range (`geom_point()`).
## Removed 2 rows containing missing values or values outside the scale range (`geom_point()`).
```

```
## Warning: Removed 2 rows containing missing values or values outside the scale range (`geom_flag()`).
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
         aes(x = year, y = uni.z_mean, color = gender)) +
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
    name   = "Mean-level of value universalism",
    sec.axis = sec_axis(~ scale_y_to_gei(.),
                        name = "Gender Equality Index (GEI)")
  ) +
  scale_x_continuous(limits = c(2001, 2024),
                     breaks = c(seq(2002, 2020, 2),2023)) +
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

```
## Warning: Removed 2 rows containing non-finite outside the scale range (`stat_smooth()`).
```

```
## Warning: Removed 2 rows containing missing values or values outside the scale range (`geom_point()`).
## Removed 2 rows containing missing values or values outside the scale range (`geom_point()`).
```

```
## Warning: Removed 12 rows containing missing values or values outside the scale range (`geom_line()`).
```

```
## Warning: Removed 2 rows containing missing values or values outside the scale range (`geom_flag()`).
```

![](Analysis_universalism_files/figure-html/unnamed-chunk-50-1.png)<!-- -->

``` r
png(filename = 
      "../results/uni/country_time_trend_facets.png",
    units = "cm",
    width = 21.0,height=29.7*(3/4),res = 600)
facet_plot
```

```
## Warning: Removed 2 rows containing non-finite outside the scale range (`stat_smooth()`).
```

```
## Warning: Removed 2 rows containing missing values or values outside the scale range (`geom_point()`).
## Removed 2 rows containing missing values or values outside the scale range (`geom_point()`).
```

```
## Warning: Removed 12 rows containing missing values or values outside the scale range (`geom_line()`).
```

```
## Warning: Removed 2 rows containing missing values or values outside the scale range (`geom_flag()`).
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
## [1] 117.9678
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
## 1   -0.22                0.33                   -0.03                    0.32                      0.35
## 2   -0.13                0.07                   -0.04                    0.06                      0.09
## 3   -0.11               -0.26                   -0.10                   -0.30                     -0.21
## 4   -0.16                0.00                   -0.07                   -0.04                      0.03
## 5    0.02                0.29                   -0.13                    0.22                      0.36
## 6   -0.27               -0.33                    0.02                   -0.32                     -0.34
## 7   -0.25                0.36                    0.08                    0.40                      0.32
## 8   -0.16                0.38                   -0.08                    0.34                      0.41
## 9   -0.20                0.05                   -0.08                    0.01                      0.09
## 10  -0.06                0.18                   -0.07                    0.15                      0.22
## 11  -0.41                0.30                    0.06                    0.33                      0.28
## 12  -0.05               -0.04                   -0.12                   -0.10                      0.02
## 13  -0.09                0.42                   -0.09                    0.38                      0.46
## 14   0.04               -0.73                   -0.06                   -0.76                     -0.70
## 15  -0.21                0.44                    0.01                    0.44                      0.43
## 16  -0.18               -0.38                    0.06                   -0.34                     -0.41
## 17  -0.14                0.03                   -0.03                    0.02                      0.05
## 18  -0.05               -0.25                   -0.11                   -0.30                     -0.19
## 19  -0.23                0.53                    0.01                    0.53                      0.52
## 20  -0.03               -0.11                   -0.11                   -0.16                     -0.06
## 21  -0.13               -0.14                   -0.11                   -0.19                     -0.09
## 22  -0.41               -0.18                    0.10                   -0.13                     -0.24
## 23  -0.13               -6.02                   -0.09                   -6.07                     -5.98
## 24  -0.15                0.21                   -0.04                    0.19                      0.23
## 25  -0.22                0.32                    0.09                    0.36                      0.28
## 26  -0.11                0.00                   -0.11                   -0.06                      0.05
## 27   0.07                0.26                   -0.23                    0.14                      0.37
## 28  -0.26                0.33                   -0.02                    0.32                      0.34
## 29  -0.12               -0.49                   -0.03                   -0.50                     -0.47
## 30  -0.34                0.66                    0.10                    0.71                      0.61
## 31  -0.16                0.52                    0.00                    0.52                      0.51
## 32  -0.12               -0.12                   -0.08                   -0.15                     -0.08
## 33   0.06               -0.30                   -0.14                   -0.37                     -0.23
## 34  -0.05                0.16                   -0.25                    0.04                      0.29
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
## 1     GR               -0.73
## 2     RU               -0.49
## 3     HU               -0.38
## 4     CZ               -0.33
## 5     BG               -0.26
## 6     IL               -0.25
## 7     LT               -0.14
## 8     SK               -0.12
## 9     IT               -0.11
## 10    FR               -0.04
## 11    CH                0.00
## 12    PL                0.00
## 13    IE                0.03
## 14    EE                0.05
## 15    BE                0.07
## 16    UA                0.16
## 17    ES                0.18
## 18    NL                0.21
## 19    PT                0.26
## 20    CY                0.29
## 21    FI                0.30
## 22    NO                0.32
## 23    AT                0.33
## 24    DE                0.36
## 25    DK                0.38
## 26    GB                0.42
## 27    HR                0.44
## 28    SI                0.52
## 29    IS                0.53
## 30    SE                0.66
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
## 1     UA                   -0.25
## 2     PT                   -0.23
## 3     CY                   -0.13
## 4     FR                   -0.12
## 5     IL                   -0.11
## 6     IT                   -0.11
## 7     LT                   -0.11
## 8     PL                   -0.11
## 9     BG                   -0.10
## 10    GB                   -0.09
## 11    DK                   -0.08
## 12    EE                   -0.08
## 13    SK                   -0.08
## 14    CH                   -0.07
## 15    ES                   -0.07
## 16    GR                   -0.06
## 17    BE                   -0.04
## 18    NL                   -0.04
## 19    AT                   -0.03
## 20    IE                   -0.03
## 21    RU                   -0.03
## 22    SI                    0.00
## 23    HR                    0.01
## 24    IS                    0.01
## 25    CZ                    0.02
## 26    FI                    0.06
## 27    HU                    0.06
## 28    DE                    0.08
## 29    NO                    0.09
## 30    SE                    0.10
```


#### Country-level predictors: GEI


``` r
mod6_GEI<-lmer(uni.z~gndr.c+year.c+
             gndr.c:year.c+
               gei.z.cm:gndr.c+gei.z.cm:year.c+gei.z.cm:gndr.c:year.c+
               (gndr.c+year.c+gndr.c:year.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_GEI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: uni.z ~ gndr.c + year.c + gndr.c:year.c + gei.z.cm:gndr.c + gei.z.cm:year.c +  
##     gei.z.cm:gndr.c:year.c + (gndr.c + year.c + gndr.c:year.c |      cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1414626.6 1414826.1 -707295.3 1414590.6    480346 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -8.4939 -0.5519  0.1047  0.6550  4.8753 
## 
## Random effects:
##  Groups   Name          Variance  Std.Dev. Corr              
##  cntry    (Intercept)   0.7265057 0.852353                   
##           gndr.c        0.0129845 0.113949  0.11             
##           year.c        0.0025991 0.050981 -0.96 -0.06       
##           gndr.c:year.c 0.0000161 0.004013 -0.15 -0.81  0.13 
##  Residual               0.9813564 0.990634                   
## Number of obs: 480364, groups:  cntry, 33
## 
## Fixed effects:
##                          Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)             0.1127741  0.1485246 31.6709235   0.759   0.4533    
## gndr.c                 -0.1514267  0.0212441 30.0880181  -7.128 6.17e-08 ***
## year.c                 -0.0053996  0.0088866 31.9203960  -0.608   0.5477    
## gndr.c:year.c          -0.0017210  0.0009032 29.0369355  -1.905   0.0667 .  
## gndr.c:gei.z.cm        -0.0361613  0.0212417 28.5096305  -1.702   0.0996 .  
## year.c:gei.z.cm         0.0068210  0.0026080 35.5419898   2.615   0.0130 *  
## gndr.c:year.c:gei.z.cm  0.0003567  0.0010019 40.3293972   0.356   0.7237    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c year.c gnd.:. gn.:.. yr.:..
## gndr.c       0.105                                   
## year.c      -0.959 -0.055                            
## gndr.c:yr.c -0.114 -0.783  0.102                     
## gndr.c:g.z.  0.000 -0.036  0.000  0.062              
## yr.c:g.z.cm  0.001  0.001 -0.006 -0.002  0.161       
## gndr.c:.:..  0.000  0.067 -0.001 -0.162 -0.724 -0.020
```

``` r
getFE(mod6_GEI,round=2,p.round = 5)
```

```
##                         Est.   SE    df     t       p    LL    UL
## (Intercept)             0.11 0.15 31.67  0.76 0.45329 -0.19  0.42
## gndr.c                 -0.15 0.02 30.09 -7.13 0.00000 -0.19 -0.11
## year.c                 -0.01 0.01 31.92 -0.61 0.54775 -0.02  0.01
## gndr.c:year.c           0.00 0.00 29.04 -1.91 0.06668  0.00  0.00
## gndr.c:gei.z.cm        -0.04 0.02 28.51 -1.70 0.09957 -0.08  0.01
## year.c:gei.z.cm         0.01 0.00 35.54  2.62 0.01299  0.00  0.01
## gndr.c:year.c:gei.z.cm  0.00 0.00 40.33  0.36 0.72366  0.00  0.00
```

``` r
getVC(mod6_GEI)
```

```
##         grp          var1          var2 sdcor  vcov
## 1     cntry   (Intercept)          <NA>  0.85  0.73
## 2     cntry        gndr.c          <NA>  0.11  0.01
## 3     cntry        year.c          <NA>  0.05  0.00
## 4     cntry gndr.c:year.c          <NA>  0.00  0.00
## 5     cntry   (Intercept)        gndr.c  0.11  0.01
## 6     cntry   (Intercept)        year.c -0.96 -0.04
## 7     cntry   (Intercept) gndr.c:year.c -0.15  0.00
## 8     cntry        gndr.c        year.c -0.06  0.00
## 9     cntry        gndr.c gndr.c:year.c -0.81  0.00
## 10    cntry        year.c gndr.c:year.c  0.13  0.00
## 11 Residual          <NA>          <NA>  0.99  0.98
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_GEI,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] 1.805622
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_GEI,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] 28.22825
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
##  year.c    emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21 -0.143859 0.0825 39.8  -0.3106   0.0229  -1.744  0.0889
##       0  0.112774 0.1490 31.7  -0.1899   0.4154   0.759  0.4533
## 
## gei.z.cm =  0:
##  year.c    emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21 -0.000617 0.0609 32.7  -0.1246   0.1233  -0.010  0.9920
##       0  0.112774 0.1490 31.7  -0.1899   0.4154   0.759  0.4533
## 
## gei.z.cm =  1:
##  year.c    emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21  0.142625 0.0813 37.7  -0.0221   0.3073   1.753  0.0876
##       0  0.112774 0.1490 31.7  -0.1899   0.4154   0.759  0.4533
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
##  year.c21 - year.c0  -0.2566 0.195 37.2   -0.651    0.138  -1.318  0.1957
## 
## gei.z.cm =  0:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0  -0.1134 0.187 31.9   -0.494    0.267  -0.608  0.5477
## 
## gei.z.cm =  1:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0   0.0299 0.194 36.7   -0.364    0.423   0.154  0.8787
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
##  gndr.c year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21 -0.0644 0.0812 40.2  -0.2286   0.0998  -0.793  0.4326
##     0.5     21 -0.2233 0.0851 39.5  -0.3954  -0.0512  -2.624  0.0123
##    -0.5      0  0.1704 0.1480 30.9  -0.1319   0.4727   1.150  0.2590
##     0.5      0  0.0551 0.1500 31.1  -0.2516   0.3619   0.367  0.7164
## 
## gei.z.cm =  0:
##  gndr.c year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21  0.0932 0.0600 32.6  -0.0290   0.2153   1.552  0.1302
##     0.5     21 -0.0944 0.0625 32.5  -0.2217   0.0329  -1.510  0.1407
##    -0.5      0  0.1885 0.1480 30.6  -0.1131   0.4901   1.275  0.2118
##     0.5      0  0.0371 0.1500 30.7  -0.2690   0.3431   0.247  0.8065
## 
## gei.z.cm =  1:
##  gndr.c year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21  0.2507 0.0800 37.9   0.0888   0.4127   3.135  0.0033
##     0.5     21  0.0345 0.0838 37.3  -0.1353   0.2043   0.412  0.6829
##    -0.5      0  0.2066 0.1480 30.9  -0.0956   0.5088   1.394  0.1732
##     0.5      0  0.0190 0.1500 31.0  -0.2877   0.3256   0.126  0.9004
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
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21     0.1589 0.0219 40.5   0.1147  0.20306   7.269 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)  -0.2348 0.1940 36.4  -0.6291  0.15951  -1.207  0.2351
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     -0.1196 0.1950 36.4  -0.5154  0.27634  -0.612  0.5442
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)     -0.3937 0.1970 37.3  -0.7922  0.00475  -2.001  0.0527
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0        -0.2784 0.1960 36.3  -0.6764  0.11946  -1.419  0.1645
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0       0.1153 0.0306 31.7   0.0529  0.17758   3.769  0.0007
## 
## gei.z.cm =  0:
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21     0.1876 0.0134 32.1   0.1603  0.21488  13.986 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)  -0.0953 0.1860 31.1  -0.4744  0.28379  -0.513  0.6117
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0      0.0561 0.1870 31.4  -0.3259  0.43814   0.299  0.7666
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)     -0.2829 0.1870 31.4  -0.6642  0.09838  -1.513  0.1404
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0        -0.1315 0.1880 31.1  -0.5145  0.25155  -0.700  0.4892
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0       0.1514 0.0212 30.1   0.1080  0.19481   7.128 <0.0001
## 
## gei.z.cm =  1:
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21     0.2162 0.0194 32.3   0.1768  0.25567  11.167 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)   0.0442 0.1940 35.8  -0.3489  0.43724   0.228  0.8210
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0      0.2318 0.1950 35.9  -0.1631  0.62658   1.191  0.2416
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)     -0.1721 0.1960 36.8  -0.5694  0.22528  -0.878  0.3859
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0         0.0155 0.1950 35.7  -0.3810  0.41208   0.079  0.9371
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0       0.1876 0.0295 28.2   0.1272  0.24798   6.360 <0.0001
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
##  diff_ESS11   -0.159 0.0219 40.5   -0.203  -0.1147  -7.269 <0.0001
##  diff_ESS1    -0.115 0.0306 31.7   -0.178  -0.0529  -3.769  0.0007
## 
## gei.z.cm =  0:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11   -0.188 0.0134 32.1   -0.215  -0.1603 -13.986 <0.0001
##  diff_ESS1    -0.151 0.0212 30.1   -0.195  -0.1080  -7.128 <0.0001
## 
## gei.z.cm =  1:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11   -0.216 0.0194 32.3   -0.256  -0.1768 -11.167 <0.0001
##  diff_ESS1    -0.188 0.0295 28.2   -0.248  -0.1272  -6.360 <0.0001
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
##  diff_ESS11 - diff_ESS1  -0.0436 0.0305 42.3  -0.1052  0.01795  -1.430  0.1602
## 
## gei.z.cm =  0:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1  -0.0361 0.0190 29.0  -0.0749  0.00265  -1.905  0.0667
## 
## gei.z.cm =  1:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1  -0.0286 0.0259 29.0  -0.0817  0.02442  -1.104  0.2786
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```


#### Country-level predictors: GGGI


``` r
mod6_GGGI<-lmer(uni.z~gndr.c+year.c+
             gndr.c:year.c+
               gggi.z.cm:gndr.c+gggi.z.cm:year.c+gggi.z.cm:gndr.c:year.c+
               (gndr.c+year.c+gndr.c:year.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_GGGI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: uni.z ~ gndr.c + year.c + gndr.c:year.c + gggi.z.cm:gndr.c +  
##     gggi.z.cm:year.c + gggi.z.cm:gndr.c:year.c + (gndr.c + year.c +      gndr.c:year.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##   1078270   1078465   -539117   1078234    363834 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -8.4755 -0.5514  0.1011  0.6566  5.0168 
## 
## Random effects:
##  Groups   Name          Variance  Std.Dev. Corr              
##  cntry    (Intercept)   1.727e+00 1.314038                   
##           gndr.c        3.355e-02 0.183166  0.88             
##           year.c        5.374e-03 0.073305 -0.98 -0.84       
##           gndr.c:year.c 7.942e-05 0.008912 -0.86 -0.93  0.84 
##  Residual               9.938e-01 0.996890                   
## Number of obs: 363852, groups:  cntry, 34
## 
## Fixed effects:
##                          Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)              0.201899   0.227977 31.705493   0.886 0.382496    
## gndr.c                  -0.106356   0.033077  4.320796  -3.215 0.029041 *  
## year.c                  -0.010954   0.012700 32.434700  -0.863 0.394732    
## gndr.c:year.c           -0.004525   0.001694  4.636572  -2.671 0.047835 *  
## gndr.c:gggi.z.cm        -0.074371   0.018895 34.828595  -3.936 0.000378 ***
## year.c:gggi.z.cm         0.005778   0.002821 32.880750   2.048 0.048608 *  
## gndr.c:year.c:gggi.z.cm  0.002759   0.001168 35.108246   2.361 0.023881 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c year.c gnd.:. gn.:.. yr.:..
## gndr.c       0.843                                   
## year.c      -0.979 -0.811                            
## gndr.c:yr.c -0.788 -0.914  0.770                     
## gndr.c:gg..  0.001 -0.007 -0.003  0.016              
## yr.c:ggg.z.  0.007  0.001 -0.023 -0.006  0.131       
## gndr.c:.:..  0.000  0.013  0.000 -0.038 -0.761 -0.010
```

``` r
getFE(mod6_GGGI,round=2,p.round = 5)
```

```
##                          Est.   SE    df     t       p    LL    UL
## (Intercept)              0.20 0.23 31.71  0.89 0.38250 -0.26  0.67
## gndr.c                  -0.11 0.03  4.32 -3.22 0.02904 -0.20 -0.02
## year.c                  -0.01 0.01 32.43 -0.86 0.39473 -0.04  0.01
## gndr.c:year.c            0.00 0.00  4.64 -2.67 0.04783 -0.01  0.00
## gndr.c:gggi.z.cm        -0.07 0.02 34.83 -3.94 0.00038 -0.11 -0.04
## year.c:gggi.z.cm         0.01 0.00 32.88  2.05 0.04861  0.00  0.01
## gndr.c:year.c:gggi.z.cm  0.00 0.00 35.11  2.36 0.02388  0.00  0.01
```

``` r
getVC(mod6_GGGI)
```

```
##         grp          var1          var2 sdcor  vcov
## 1     cntry   (Intercept)          <NA>  1.31  1.73
## 2     cntry        gndr.c          <NA>  0.18  0.03
## 3     cntry        year.c          <NA>  0.07  0.01
## 4     cntry gndr.c:year.c          <NA>  0.01  0.00
## 5     cntry   (Intercept)        gndr.c  0.88  0.21
## 6     cntry   (Intercept)        year.c -0.98 -0.09
## 7     cntry   (Intercept) gndr.c:year.c -0.86 -0.01
## 8     cntry        gndr.c        year.c -0.84 -0.01
## 9     cntry        gndr.c gndr.c:year.c -0.93  0.00
## 10    cntry        year.c gndr.c:year.c  0.84  0.00
## 11 Residual          <NA>          <NA>  1.00  0.99
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_GGGI,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] -103.021
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_GGGI,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] -253.9522
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
##  year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21 -0.1495 0.0899 38.2  -0.3314   0.0324  -1.663  0.1045
##       0  0.2019 0.2280 31.7  -0.2626   0.6664   0.886  0.3825
## 
## gggi.z.cm =  0:
##  year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21 -0.0281 0.0635 33.2  -0.1573   0.1010  -0.443  0.6606
##       0  0.2019 0.2280 31.7  -0.2626   0.6664   0.886  0.3825
## 
## gggi.z.cm =  1:
##  year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21  0.0932 0.0837 38.3  -0.0762   0.2626   1.113  0.2725
##       0  0.2019 0.2280 31.7  -0.2626   0.6664   0.886  0.3825
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
##  year.c21 - year.c0   -0.351 0.275 36.1   -0.908    0.205  -1.280  0.2087
## 
## gggi.z.cm =  0:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0   -0.230 0.267 32.4   -0.773    0.313  -0.863  0.3947
## 
## gggi.z.cm =  1:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0   -0.109 0.272 34.9   -0.661    0.443  -0.400  0.6917
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
##    -0.5     21 -0.0570 0.0891 38.1  -0.2373  0.12327  -0.640  0.5260
##     0.5     21 -0.2419 0.0920 37.9  -0.4283 -0.05564  -2.629  0.0123
##    -0.5      0  0.2179 0.2140 30.2  -0.2199  0.65570   1.016  0.3177
##     0.5      0  0.1859 0.2420 30.6  -0.3085  0.68027   0.767  0.4487
## 
## gggi.z.cm =  0:
##  gndr.c year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21  0.0726 0.0630 33.2  -0.0556  0.20070   1.152  0.2576
##     0.5     21 -0.1288 0.0648 33.1  -0.2607  0.00303  -1.988  0.0552
##    -0.5      0  0.2551 0.2140 30.1  -0.1824  0.69251   1.191  0.2431
##     0.5      0  0.1487 0.2420 30.5  -0.3453  0.64277   0.614  0.5435
## 
## gggi.z.cm =  1:
##  gndr.c year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21  0.2021 0.0830 38.3   0.0342  0.37008   2.435  0.0196
##     0.5     21 -0.0157 0.0857 38.0  -0.1892  0.15779  -0.183  0.8556
##    -0.5      0  0.2923 0.2140 30.2  -0.1455  0.73003   1.363  0.1829
##     0.5      0  0.1115 0.2420 30.6  -0.3828  0.60590   0.460  0.6485
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod6_GGGI,adjust="none",infer=c(T,T))
```

```
## gggi.z.cm = -1:
##  contrast                                     estimate     SE    df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21     0.1850 0.0223 36.02   0.1398   0.2301   8.310 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)  -0.2749 0.2620 35.07  -0.8063   0.2566  -1.050  0.3009
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     -0.2429 0.2870 34.69  -0.8264   0.3405  -0.845  0.4037
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)     -0.4598 0.2630 35.63  -0.9934   0.0738  -1.748  0.0890
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0        -0.4279 0.2880 34.67  -1.0135   0.1578  -1.484  0.1469
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0       0.0320 0.0382  7.21  -0.0578   0.1218   0.837  0.4294
## 
## gggi.z.cm =  0:
##  contrast                                     estimate     SE    df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21     0.2014 0.0145 34.14   0.1720   0.2308  13.935 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)  -0.1825 0.2530 31.09  -0.6990   0.3339  -0.721  0.4765
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     -0.0762 0.2800 31.57  -0.6470   0.4947  -0.272  0.7875
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)     -0.3839 0.2540 31.38  -0.9015   0.1337  -1.512  0.1405
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0        -0.2775 0.2810 31.39  -0.8496   0.2945  -0.989  0.3302
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0       0.1064 0.0331  4.32   0.0171   0.1956   3.215  0.0290
## 
## gggi.z.cm =  1:
##  contrast                                     estimate     SE    df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21     0.2178 0.0207 34.87   0.1758   0.2599  10.518 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)  -0.0901 0.2590 33.75  -0.6166   0.4363  -0.348  0.7300
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0      0.0906 0.2850 33.61  -0.4886   0.6698   0.318  0.7525
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)     -0.3080 0.2600 34.23  -0.8364   0.2204  -1.184  0.2445
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0        -0.1272 0.2860 33.54  -0.7082   0.4537  -0.445  0.6590
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0       0.1807 0.0380  7.01   0.0910   0.2705   4.759  0.0021
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
##  contrast   estimate     SE    df lower.CL upper.CL t.ratio p.value
##  diff_ESS11   -0.185 0.0223 36.02   -0.230  -0.1398  -8.310 <0.0001
##  diff_ESS1    -0.032 0.0382  7.21   -0.122   0.0578  -0.837  0.4294
## 
## gggi.z.cm =  0:
##  contrast   estimate     SE    df lower.CL upper.CL t.ratio p.value
##  diff_ESS11   -0.201 0.0145 34.14   -0.231  -0.1720 -13.935 <0.0001
##  diff_ESS1    -0.106 0.0331  4.32   -0.196  -0.0171  -3.215  0.0290
## 
## gggi.z.cm =  1:
##  contrast   estimate     SE    df lower.CL upper.CL t.ratio p.value
##  diff_ESS11   -0.218 0.0207 34.87   -0.260  -0.1758 -10.518 <0.0001
##  diff_ESS1    -0.181 0.0380  7.01   -0.270  -0.0910  -4.759  0.0021
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
##  diff_ESS11 - diff_ESS1  -0.1530 0.0440 8.99   -0.252 -0.05348  -3.479  0.0070
## 
## gggi.z.cm =  0:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1  -0.0950 0.0356 4.64   -0.189 -0.00137  -2.671  0.0478
## 
## gggi.z.cm =  1:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1  -0.0371 0.0425 7.94   -0.135  0.06094  -0.874  0.4079
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

#### Country-level predictors: GDI


``` r
mod6_GDI<-lmer(uni.z~gndr.c+year.c+
             gndr.c:year.c+
               gdi.z.cm:gndr.c+gdi.z.cm:year.c+gdi.z.cm:gndr.c:year.c+
               (gndr.c+year.c+gndr.c:year.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_GDI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: uni.z ~ gndr.c + year.c + gndr.c:year.c + gdi.z.cm:gndr.c + gdi.z.cm:year.c +  
##     gdi.z.cm:gndr.c:year.c + (gndr.c + year.c + gndr.c:year.c |      cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1458485.1 1458685.1 -729224.6 1458449.1    492325 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -8.4522 -0.5502  0.1038  0.6531  4.8403 
## 
## Random effects:
##  Groups   Name          Variance  Std.Dev. Corr              
##  cntry    (Intercept)   7.141e-01 0.845021                   
##           gndr.c        1.340e-02 0.115737  0.04             
##           year.c        2.797e-03 0.052885 -0.96 -0.07       
##           gndr.c:year.c 2.103e-05 0.004586 -0.10 -0.82  0.11 
##  Residual               9.945e-01 0.997255                   
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##                          Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)             0.0982136  0.1450658 32.3612554   0.677  0.50320    
## gndr.c                 -0.1461397  0.0212540 31.6344765  -6.876 9.38e-08 ***
## year.c                 -0.0048032  0.0090810 32.5503663  -0.529  0.60044    
## gndr.c:year.c          -0.0022533  0.0009643 29.3293183  -2.337  0.02648 *  
## gndr.c:gdi.z.cm        -0.0310964  0.0218008 33.1735169  -1.426  0.16310    
## year.c:gdi.z.cm        -0.0077009  0.0024974 36.0908002  -3.084  0.00391 ** 
## gndr.c:year.c:gdi.z.cm -0.0007030  0.0010974 43.2695132  -0.641  0.52516    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c year.c gnd.:. gn.:.. yr.:..
## gndr.c       0.036                                   
## year.c      -0.964 -0.068                            
## gndr.c:yr.c -0.081 -0.808  0.091                     
## gndr.c:gd..  0.000 -0.020  0.000  0.024              
## yr.c:gd.z.c  0.001  0.000 -0.004  0.001 -0.123       
## gndr.c:.:..  0.000  0.022  0.001 -0.046 -0.772  0.042
```

``` r
getFE(mod6_GDI,round=2,p.round = 5)
```

```
##                         Est.   SE    df     t       p    LL    UL
## (Intercept)             0.10 0.15 32.36  0.68 0.50320 -0.20  0.39
## gndr.c                 -0.15 0.02 31.63 -6.88 0.00000 -0.19 -0.10
## year.c                  0.00 0.01 32.55 -0.53 0.60044 -0.02  0.01
## gndr.c:year.c           0.00 0.00 29.33 -2.34 0.02648  0.00  0.00
## gndr.c:gdi.z.cm        -0.03 0.02 33.17 -1.43 0.16310 -0.08  0.01
## year.c:gdi.z.cm        -0.01 0.00 36.09 -3.08 0.00391 -0.01  0.00
## gndr.c:year.c:gdi.z.cm  0.00 0.00 43.27 -0.64 0.52516  0.00  0.00
```

``` r
getVC(mod6_GDI)
```

```
##         grp          var1          var2 sdcor  vcov
## 1     cntry   (Intercept)          <NA>  0.85  0.71
## 2     cntry        gndr.c          <NA>  0.12  0.01
## 3     cntry        year.c          <NA>  0.05  0.00
## 4     cntry gndr.c:year.c          <NA>  0.00  0.00
## 5     cntry   (Intercept)        gndr.c  0.04  0.00
## 6     cntry   (Intercept)        year.c -0.96 -0.04
## 7     cntry   (Intercept) gndr.c:year.c -0.10  0.00
## 8     cntry        gndr.c        year.c -0.07  0.00
## 9     cntry        gndr.c gndr.c:year.c -0.82  0.00
## 10    cntry        year.c gndr.c:year.c  0.11  0.00
## 11 Residual          <NA>          <NA>  1.00  0.99
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_GDI,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] -5.667803
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_GDI,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] 6.264179
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
##  year.c   emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21  0.15907 0.0830 42.6 -0.00842  0.32655   1.916  0.0621
##       0  0.09821 0.1450 32.4 -0.19715  0.39357   0.677  0.5032
## 
## gdi.z.cm =  0:
##  year.c   emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21 -0.00265 0.0638 32.8 -0.13242  0.12711  -0.042  0.9671
##       0  0.09821 0.1450 32.4 -0.19715  0.39357   0.677  0.5032
## 
## gdi.z.cm =  1:
##  year.c   emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21 -0.16437 0.0821 40.9 -0.33019  0.00145  -2.002  0.0519
##       0  0.09821 0.1450 32.4 -0.19715  0.39357   0.677  0.5032
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
##  year.c21 - year.c0   0.0609 0.198 37.3   -0.340    0.462   0.307  0.7603
## 
## gdi.z.cm =  0:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0  -0.1009 0.191 32.5   -0.489    0.287  -0.529  0.6004
## 
## gdi.z.cm =  1:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0  -0.2626 0.198 37.0   -0.663    0.138  -1.329  0.1920
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
##    -0.5     21  0.2329 0.0845 42.3   0.0623   0.4034   2.755  0.0086
##     0.5     21  0.0853 0.0828 43.0  -0.0817   0.2522   1.030  0.3087
##    -0.5      0  0.1557 0.1460 31.5  -0.1408   0.4523   1.070  0.2926
##     0.5      0  0.0407 0.1460 31.6  -0.2573   0.3387   0.278  0.7826
## 
## gdi.z.cm =  0:
##  gndr.c year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21  0.0941 0.0645 32.8  -0.0372   0.2253   1.459  0.1542
##     0.5     21 -0.0994 0.0637 32.8  -0.2290   0.0302  -1.561  0.1282
##    -0.5      0  0.1713 0.1450 31.2  -0.1245   0.4671   1.181  0.2467
##     0.5      0  0.0251 0.1460 31.3  -0.2722   0.3225   0.172  0.8642
## 
## gdi.z.cm =  1:
##  gndr.c year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21 -0.0447 0.0835 40.5  -0.2134   0.1240  -0.535  0.5953
##     0.5     21 -0.2840 0.0818 41.2  -0.4492  -0.1188  -3.471  0.0012
##    -0.5      0  0.1868 0.1450 31.5  -0.1097   0.4833   1.284  0.2084
##     0.5      0  0.0096 0.1460 31.6  -0.2884   0.3076   0.066  0.9481
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
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21     0.1476 0.0204 34.6   0.1062   0.1890   7.245 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)   0.0771 0.1980 36.5  -0.3235   0.4777   0.390  0.6986
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0      0.1922 0.2000 37.3  -0.2133   0.5977   0.960  0.3432
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)     -0.0705 0.1980 36.6  -0.4718   0.3308  -0.356  0.7239
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0         0.0446 0.2000 36.6  -0.3600   0.4492   0.223  0.8245
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0       0.1150 0.0307 32.5   0.0525   0.1776   3.742  0.0007
## 
## gdi.z.cm =  0:
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21     0.1935 0.0129 31.2   0.1672   0.2197  15.009 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)  -0.0772 0.1900 31.7  -0.4645   0.3100  -0.406  0.6873
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0      0.0689 0.1920 32.0  -0.3217   0.4595   0.359  0.7216
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)     -0.2707 0.1910 32.0  -0.6591   0.1178  -1.419  0.1655
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0        -0.1245 0.1920 31.8  -0.5155   0.2664  -0.649  0.5210
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0       0.1461 0.0213 31.6   0.1028   0.1895   6.876 <0.0001
## 
## gdi.z.cm =  1:
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21     0.2393 0.0194 36.5   0.1999   0.2787  12.306 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)  -0.2315 0.1970 36.1  -0.6312   0.1681  -1.175  0.2478
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     -0.0543 0.2000 37.0  -0.4590   0.3503  -0.272  0.7872
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)     -0.4709 0.1980 36.3  -0.8714  -0.0704  -2.384  0.0225
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0        -0.2936 0.1990 36.3  -0.6974   0.1102  -1.474  0.1490
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0       0.1772 0.0301 32.5   0.1159   0.2386   5.879 <0.0001
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
##  diff_ESS11   -0.148 0.0204 34.6   -0.189  -0.1062  -7.245 <0.0001
##  diff_ESS1    -0.115 0.0307 32.5   -0.178  -0.0525  -3.742  0.0007
## 
## gdi.z.cm =  0:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11   -0.193 0.0129 31.2   -0.220  -0.1672 -15.009 <0.0001
##  diff_ESS1    -0.146 0.0213 31.6   -0.189  -0.1028  -6.876 <0.0001
## 
## gdi.z.cm =  1:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11   -0.239 0.0194 36.5   -0.279  -0.1999 -12.306 <0.0001
##  diff_ESS1    -0.177 0.0301 32.5   -0.239  -0.1159  -5.879 <0.0001
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
##  diff_ESS11 - diff_ESS1  -0.0326 0.0314 35.7  -0.0962  0.03108  -1.038  0.3063
## 
## gdi.z.cm =  0:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1  -0.0473 0.0202 29.3  -0.0887 -0.00592  -2.337  0.0265
## 
## gdi.z.cm =  1:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1  -0.0621 0.0300 38.2  -0.1227 -0.00142  -2.071  0.0451
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

#### Country-level predictors: log_GDP


``` r
mod6_log_GDP<-lmer(uni.z~gndr.c+year.c+
             gndr.c:year.c+
               log_gdp.z.cm:gndr.c+log_gdp.z.cm:year.c+log_gdp.z.cm:gndr.c:year.c+
               (gndr.c+year.c+gndr.c:year.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
```

```
## Warning: Model failed to converge with 1 negative eigenvalue: -2.9e+03
```

``` r
summary(mod6_log_GDP)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: uni.z ~ gndr.c + year.c + gndr.c:year.c + log_gdp.z.cm:gndr.c +  
##     log_gdp.z.cm:year.c + log_gdp.z.cm:gndr.c:year.c + (gndr.c +      year.c + gndr.c:year.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1458487.4 1458687.3 -729225.7 1458451.4    492325 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -8.4486 -0.5499  0.1037  0.6534  4.8355 
## 
## Random effects:
##  Groups   Name          Variance  Std.Dev. Corr              
##  cntry    (Intercept)   7.087e-01 0.8418                     
##           gndr.c        1.295e-02 0.1138    0.00             
##           year.c        2.401e-03 0.0490   -0.96  0.04       
##           gndr.c:year.c 1.764e-05 0.0042    0.13 -0.78 -0.18 
##  Residual               9.945e-01 0.9973                     
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##                              Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)                 0.0984568  0.1445225 32.5951483   0.681  0.50052    
## gndr.c                     -0.1435749  0.0209668 31.0915649  -6.848 1.11e-07 ***
## year.c                     -0.0047768  0.0084149 32.5592597  -0.568  0.57416    
## gndr.c:year.c              -0.0025335  0.0009149 31.2423329  -2.769  0.00937 ** 
## gndr.c:log_gdp.z.cm        -0.0329009  0.0218418 32.6171372  -1.506  0.14161    
## year.c:log_gdp.z.cm         0.0074267  0.0024366 34.3540890   3.048  0.00441 ** 
## gndr.c:year.c:log_gdp.z.cm  0.0024187  0.0009671 31.2977945   2.501  0.01783 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c year.c gnd.:. g.:_.. y.:_..
## gndr.c      -0.001                                   
## year.c      -0.958  0.041                            
## gndr.c:yr.c  0.102 -0.770 -0.139                     
## gndr.c:l_..  0.000 -0.061  0.001  0.095              
## yr.c:lg_g..  0.000  0.004  0.005 -0.004  0.134       
## gndr.:.:_..  0.000  0.092 -0.002 -0.145 -0.789 -0.136
```

``` r
getFE(mod6_log_GDP,round=2,p.round = 5)
```

```
##                             Est.   SE    df     t       p    LL    UL
## (Intercept)                 0.10 0.14 32.60  0.68 0.50052 -0.20  0.39
## gndr.c                     -0.14 0.02 31.09 -6.85 0.00000 -0.19 -0.10
## year.c                      0.00 0.01 32.56 -0.57 0.57416 -0.02  0.01
## gndr.c:year.c               0.00 0.00 31.24 -2.77 0.00937  0.00  0.00
## gndr.c:log_gdp.z.cm        -0.03 0.02 32.62 -1.51 0.14161 -0.08  0.01
## year.c:log_gdp.z.cm         0.01 0.00 34.35  3.05 0.00441  0.00  0.01
## gndr.c:year.c:log_gdp.z.cm  0.00 0.00 31.30  2.50 0.01783  0.00  0.00
```

``` r
getVC(mod6_log_GDP)
```

```
##         grp          var1          var2 sdcor  vcov
## 1     cntry   (Intercept)          <NA>  0.84  0.71
## 2     cntry        gndr.c          <NA>  0.11  0.01
## 3     cntry        year.c          <NA>  0.05  0.00
## 4     cntry gndr.c:year.c          <NA>  0.00  0.00
## 5     cntry   (Intercept)        gndr.c  0.00  0.00
## 6     cntry   (Intercept)        year.c -0.96 -0.04
## 7     cntry   (Intercept) gndr.c:year.c  0.13  0.00
## 8     cntry        gndr.c        year.c  0.04  0.00
## 9     cntry        gndr.c gndr.c:year.c -0.78  0.00
## 10    cntry        year.c gndr.c:year.c -0.18  0.00
## 11 Residual          <NA>          <NA>  1.00  0.99
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_log_GDP,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] 9.300523
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_log_GDP,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] 21.36551
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
##  year.c   emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21 -0.15782 0.0757 37.9 -0.31110 -0.00453  -2.084  0.0439
##       0  0.09846 0.1450 32.6 -0.19572  0.39263   0.681  0.5005
## 
## log_gdp.z.cm =  0:
##  year.c   emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21 -0.00186 0.0566 33.4 -0.11693  0.11322  -0.033  0.9740
##       0  0.09846 0.1450 32.6 -0.19572  0.39263   0.681  0.5005
## 
## log_gdp.z.cm =  1:
##  year.c   emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21  0.15410 0.0769 37.1 -0.00162  0.30983   2.005  0.0523
##       0  0.09846 0.1450 32.6 -0.19572  0.39263   0.681  0.5005
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
##  year.c21 - year.c0  -0.2563 0.184 37.4   -0.628    0.116  -1.395  0.1713
## 
## log_gdp.z.cm =  0:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0  -0.1003 0.177 32.6   -0.460    0.259  -0.568  0.5742
## 
## log_gdp.z.cm =  1:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0   0.0556 0.184 37.6   -0.317    0.429   0.302  0.7642
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
##    -0.5     21 -0.0505 0.0766 38.1  -0.2054    0.104  -0.659  0.5136
##     0.5     21 -0.2652 0.0762 37.5  -0.4194   -0.111  -3.482  0.0013
##    -0.5      0  0.1538 0.1450 33.3  -0.1418    0.449   1.058  0.2977
##     0.5      0  0.0431 0.1450 33.4  -0.2525    0.339   0.297  0.7686
## 
## log_gdp.z.cm =  0:
##  gndr.c year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21  0.0965 0.0573 33.4  -0.0201    0.213   1.683  0.1016
##     0.5     21 -0.1002 0.0567 33.3  -0.2155    0.015  -1.769  0.0860
##    -0.5      0  0.1702 0.1450 32.9  -0.1246    0.465   1.175  0.2485
##     0.5      0  0.0267 0.1450 33.0  -0.2681    0.321   0.184  0.8551
## 
## log_gdp.z.cm =  1:
##  gndr.c year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21  0.2435 0.0777 37.3   0.0862    0.401   3.136  0.0033
##     0.5     21  0.0647 0.0773 36.7  -0.0919    0.221   0.837  0.4080
##    -0.5      0  0.1867 0.1450 33.3  -0.1088    0.482   1.285  0.2076
##     0.5      0  0.0102 0.1450 33.3  -0.2852    0.306   0.070  0.9443
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod6_log_GDP,adjust="none",infer=c(T,T))
```

```
## log_gdp.z.cm = -1:
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21     0.2147 0.0197 35.5   0.1747   0.2547  10.887 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)  -0.2043 0.1860 38.2  -0.5807   0.1722  -1.098  0.2789
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     -0.0936 0.1850 37.6  -0.4683   0.2811  -0.506  0.6159
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)     -0.4189 0.1850 38.7  -0.7931  -0.0448  -2.265  0.0292
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0        -0.3083 0.1830 37.2  -0.6784   0.0619  -1.687  0.0999
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0       0.1107 0.0312 34.6   0.0473   0.1740   3.548  0.0011
## 
## log_gdp.z.cm =  0:
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21     0.1968 0.0137 32.3   0.1689   0.2247  14.349 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)  -0.0737 0.1780 33.0  -0.4365   0.2891  -0.413  0.6820
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0      0.0699 0.1780 32.6  -0.2921   0.4318   0.393  0.6969
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)     -0.2705 0.1770 33.1  -0.6303   0.0893  -1.529  0.1357
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0        -0.1269 0.1760 32.3  -0.4845   0.2307  -0.723  0.4751
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0       0.1436 0.0210 31.1   0.1008   0.1863   6.848 <0.0001
## 
## log_gdp.z.cm =  1:
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21     0.1789 0.0191 30.7   0.1398   0.2179   9.349 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)   0.0569 0.1860 38.3  -0.3203   0.4340   0.305  0.7619
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0      0.2333 0.1850 37.7  -0.1421   0.6088   1.258  0.2160
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)     -0.1220 0.1850 38.9  -0.4971   0.2530  -0.658  0.5143
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0         0.0544 0.1830 37.2  -0.3162   0.4251   0.298  0.7677
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0       0.1765 0.0293 29.1   0.1165   0.2365   6.016 <0.0001
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
##  diff_ESS11   -0.215 0.0197 35.5   -0.255  -0.1747 -10.887 <0.0001
##  diff_ESS1    -0.111 0.0312 34.6   -0.174  -0.0473  -3.548  0.0011
## 
## log_gdp.z.cm =  0:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11   -0.197 0.0137 32.3   -0.225  -0.1689 -14.349 <0.0001
##  diff_ESS1    -0.144 0.0210 31.1   -0.186  -0.1008  -6.848 <0.0001
## 
## log_gdp.z.cm =  1:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11   -0.179 0.0191 30.7   -0.218  -0.1398  -9.349 <0.0001
##  diff_ESS1    -0.176 0.0293 29.1   -0.236  -0.1165  -6.016 <0.0001
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
##  diff_ESS11 - diff_ESS1 -0.10400 0.0299 35.8  -0.1647  -0.0433  -3.476  0.0014
## 
## log_gdp.z.cm =  0:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1 -0.05320 0.0192 31.2  -0.0924  -0.0140  -2.769  0.0094
## 
## log_gdp.z.cm =  1:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1 -0.00241 0.0258 26.6  -0.0555   0.0507  -0.093  0.9263
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
##  [1] tibble_3.3.1          apaTables_2.0.8       stringr_1.6.0         tidyr_1.3.2          
##  [5] r2mlm_0.3.8           nlme_3.1-169          Hmisc_5.2-5           ggpubr_0.6.3         
##  [9] metafor_5.0-1         numDeriv_2016.8-1.1   metadat_1.6-0         lmerTest_3.2-1       
## [13] ggflags_0.0.4         finalfit_1.1.0        ggplot2_4.0.3         MetBrewer_0.2.0      
## [17] vjihelpers_0.0.0.9000 emmeans_2.0.3         lme4_2.0-1            Matrix_1.7-5         
## [21] dplyr_1.2.1           rio_1.3.0             multid_1.0.2          knitr_1.51           
## [25] rmarkdown_2.31       
## 
## loaded via a namespace (and not attached):
##   [1] Rdpack_2.6.6       mnormt_2.1.2       gridExtra_2.3      writexl_1.5.4      readxl_1.4.5      
##   [6] rlang_1.2.0        magrittr_2.0.5     otel_0.2.0         rockchalk_1.8.164  compiler_4.6.0    
##  [11] mgcv_1.9-4         png_0.1-9          vctrs_0.7.3        quadprog_1.5-8     pkgconfig_2.0.3   
##  [16] shape_1.4.6.1      fastmap_1.2.0      backports_1.5.1    labeling_0.4.3     pbivnorm_0.6.0    
##  [21] utf8_1.2.6         nloptr_2.2.1       purrr_1.2.2        xfun_0.57          glmnet_5.0        
##  [26] jomo_2.7-6         cachem_1.1.0       kutils_1.73        jsonlite_2.0.0     pan_1.9           
##  [31] jpeg_0.1-11        psych_2.6.3        lavaan_0.6-21      broom_1.0.13       parallel_4.6.0    
##  [36] cluster_2.1.8.2    R6_2.6.1           bslib_0.10.0       stringi_1.8.7      RColorBrewer_1.1-3
##  [41] car_3.1-5          boot_1.3-32        rpart_4.1.27       cellranger_1.1.0   jquerylib_0.1.4   
##  [46] estimability_1.5.1 Rcpp_1.1.1-1.1     iterators_1.0.14   base64enc_0.1-6    R.utils_2.13.0    
##  [51] splines_4.6.0      nnet_7.3-20        tidyselect_1.2.1   rstudioapi_0.18.0  abind_1.4-8       
##  [56] yaml_2.3.12        codetools_0.2-20   plyr_1.8.9         lattice_0.22-9     withr_3.0.2       
##  [61] S7_0.2.2           coda_0.19-4.1      evaluate_1.0.5     foreign_0.8-91     survival_3.8-6    
##  [66] zip_2.3.3          pillar_1.11.1      carData_3.0-6      mice_3.19.0        stats4_4.6.0      
##  [71] checkmate_2.3.4    foreach_1.5.2      reformulas_0.4.4   generics_0.1.4     grImport2_0.3-3   
##  [76] mathjaxr_2.0-0     scales_1.4.0       minqa_1.2.8        xtable_1.8-8       glue_1.8.1        
##  [81] tools_4.6.0        data.table_1.18.4  openxlsx_4.2.8.1   ggsignif_0.6.4     forcats_1.0.1     
##  [86] XML_3.99-0.23      mvtnorm_1.3-7      cowplot_1.2.0      grid_4.6.0         rbibutils_2.4.1   
##  [91] colorspace_2.1-2   htmlTable_2.5.0    Formula_1.2-5      cli_3.6.6          gtable_0.3.6      
##  [96] R.methodsS3_1.8.2  rstatix_0.7.3      sass_0.4.10        digest_0.6.39      htmlwidgets_1.6.4 
## [101] farver_2.1.2       htmltools_0.5.9    R.oo_1.27.1        lifecycle_1.0.5    mitml_0.4-5       
## [106] MASS_7.3-65
```

