---
title: "Analysis for benevolence values"
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
cntry.ben<-diff_dat %>% group_by(cntry,essround) %>%
  summarise(ben.ctm=mean(ben,na.rm=T),
            ben.ctsd=sd(ben,na.rm=T)) %>%
  group_by(cntry) %>%
  summarise(ben.cm=mean(ben.ctm),
            ben.csd=mean(ben.ctsd)) 
```

```
## `summarise()` has regrouped the output.
## ℹ Summaries were computed grouped by cntry and essround.
## ℹ Output is grouped by cntry.
## ℹ Use `summarise(.groups = "drop_last")` to silence this message.
## ℹ Use `summarise(.by = c(cntry, essround))` for per-operation grouping (`?dplyr::dplyr_by`) instead.
```

``` r
grand_mean_ben<-mean(cntry.ben$ben.cm)
grand_sd_ben<-mean(cntry.ben$ben.csd)

# standardized
diff_dat$ben.z<-(diff_dat$ben-grand_mean_ben)/grand_sd_ben
hist(diff_dat$ben.z)
```

![](Analysis_benevolence_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

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
                   ben.z.wt=weighted.mean(x=ben.z,w=pspwght),
                   ben.z=mean(ben.z),
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

# benevolence

cntry_ben_frame<-
  diff_dat %>% group_by(cntry,essround) %>%
  summarise('ben M' = weighted.mean(x=ben.z,w=pspwght),
            'ben SD' = sqrt(wtd.var(ben.z,pspwght))) %>%
  group_by(cntry) %>%
  summarise('ben M' = mean(x=`ben M`),
            'ben SD'= mean(x=`ben SD`))
```

```
## `summarise()` has regrouped the output.
## ℹ Summaries were computed grouped by cntry and essround.
## ℹ Output is grouped by cntry.
## ℹ Use `summarise(.groups = "drop_last")` to silence this message.
## ℹ Use `summarise(.by = c(cntry, essround))` for per-operation grouping (`?dplyr::dplyr_by`) instead.
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
## `summarise()` has regrouped the output.
## ℹ Summaries were computed grouped by cntry and essround.
## ℹ Output is grouped by cntry.
## ℹ Use `summarise(.groups = "drop_last")` to silence this message.
## ℹ Use `summarise(.by = c(cntry, essround))` for per-operation grouping (`?dplyr::dplyr_by`) instead.
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
## `summarise()` has regrouped the output.
## ℹ Summaries were computed grouped by cntry and essround.
## ℹ Output is grouped by cntry.
## ℹ Use `summarise(.groups = "drop_last")` to silence this message.
## ℹ Use `summarise(.by = c(cntry, essround))` for per-operation grouping (`?dplyr::dplyr_by`) instead.
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
## # A tibble: 34 × 10
##    cntry `n ESS rounds`     n `ben M` `ben SD` `ben M Women` `ben SD Women` `ben M Men` `ben SD Men`
##    <chr>          <dbl> <dbl>   <dbl>    <dbl>         <dbl>          <dbl>       <dbl>        <dbl>
##  1 AT                 7 15400  0.204     0.959        0.314           0.933      0.0866        0.971
##  2 BE                11 18886  0.225     0.827        0.339           0.814      0.105         0.824
##  3 BG                 7 14857  0.0249    1.06         0.0699          1.06      -0.0236        1.06 
##  4 CH                11 18087  0.352     0.807        0.458           0.785      0.242         0.815
##  5 CY                 6  5771  0.396     0.842        0.476           0.806      0.312         0.867
##  6 CZ                 9 18934 -0.481     1.12        -0.338           1.10      -0.636         1.11 
##  7 DE                10 27753  0.234     0.863        0.336           0.846      0.126         0.868
##  8 DK                 8 12198  0.360     0.858        0.498           0.809      0.218         0.883
##  9 EE                10 17974 -0.199     0.953       -0.0656          0.927     -0.358         0.958
## 10 ES                10 18785  0.388     0.875        0.459           0.853      0.314         0.891
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
print(cntry_desc_tbl,n=35)
```

```
## # A tibble: 34 × 14
##    Country     `n ESS rounds`     n `ben M` `ben SD` `ben M Women` `ben SD Women` `ben M Men` `ben SD Men`
##    <chr>                <dbl> <dbl> <chr>   <chr>    <chr>         <chr>          <chr>       <chr>       
##  1 Austria                  7 15400 0.20    0.96     0.31          0.93           0.09        0.97        
##  2 Belgium                 11 18886 0.22    0.83     0.34          0.81           0.10        0.82        
##  3 Bulgaria                 7 14857 0.02    1.06     0.07          1.06           -0.02       1.06        
##  4 Switzerland             11 18087 0.35    0.81     0.46          0.78           0.24        0.82        
##  5 Cyprus                   6  5771 0.40    0.84     0.48          0.81           0.31        0.87        
##  6 Czechia                  9 18934 -0.48   1.12     -0.34         1.10           -0.64       1.11        
##  7 Germany                 10 27753 0.23    0.86     0.34          0.85           0.13        0.87        
##  8 Denmark                  8 12198 0.36    0.86     0.50          0.81           0.22        0.88        
##  9 Estonia                 10 17974 -0.20   0.95     -0.07         0.93           -0.36       0.96        
## 10 Spain                   10 18785 0.39    0.88     0.46          0.85           0.31        0.89        
## 11 Finland                 11 19568 0.09    0.93     0.26          0.90           -0.10       0.93        
## 12 France                  11 20457 -0.02   1.14     0.11          1.12           -0.15       1.14        
## 13 UK                      11 22979 0.15    0.96     0.28          0.93           0.02        0.98        
## 14 Greece                   6 15212 0.17    0.95     0.21          0.95           0.12        0.95        
## 15 Croatia                  5  7914 0.12    1.04     0.20          1.00           0.03        1.07        
## 16 Hungary                 11 18123 -0.07   1.07     0.02          1.05           -0.16       1.09        
## 17 Ireland                 11 22562 0.08    1.05     0.20          1.02           -0.05       1.05        
## 18 Israel                   7 14857 0.18    1.06     0.23          1.05           0.13        1.07        
## 19 Iceland                  6  4654 0.41    0.85     0.57          0.79           0.25        0.88        
## 20 Italy                    5 11441 -0.05   1.02     0.02          1.01           -0.13       1.02        
## 21 Lithuania                7 13059 -0.89   1.26     -0.79         1.25           -1.02       1.25        
## 22 Latvia                   3  4088 -0.27   1.08     -0.11         1.03           -0.46       1.11        
## 23 Montenegro               3  4028 -0.40   1.24     -0.33         1.23           -0.46       1.24        
## 24 Netherlands             11 19722 -0.00   0.87     0.13          0.84           -0.14       0.87        
## 25 Norway                  11 16505 0.03    0.93     0.17          0.92           -0.11       0.93        
## 26 Poland                  10 16737 -0.07   0.93     0.01          0.93           -0.15       0.94        
## 27 Portugal                11 19070 -0.18   1.06     -0.13         1.05           -0.23       1.06        
## 28 Serbia                   2  3499 0.36    0.93     0.46          0.90           0.26        0.95        
## 29 Russia                   5 12139 -0.31   1.16     -0.27         1.17           -0.36       1.15        
## 30 Sweden                  10 16104 0.05    0.95     0.21          0.91           -0.13       0.95        
## 31 Slovenia                11 14463 0.05    0.89     0.16          0.89           -0.06       0.89        
## 32 Slovakia                 8 12547 -0.37   1.01     -0.28         1.02           -0.46       1.00        
## 33 Turkey                   2  4108 0.14    0.98     0.13          0.97           0.14        0.98        
## 34 Ukraine                  6 12054 -0.49   1.28     -0.42         1.28           -0.59       1.27        
## # ℹ 5 more variables: D <chr>, GEI <chr>, GGGI <chr>, GDI <chr>, GDP <chr>
```

``` r
export(cntry_desc_tbl,"../results/ben/cntry_desc_tbl.xlsx",overwrite=T)
```

## Country-level correlations table


``` r
cor_frame<-
  desc_frame %>%
  dplyr::select(
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
  dplyr::select(-GDP)

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
##   1. VBMT       0.01  0.30                                                                            
##                                                                                                       
##   2. VBMT_Women 0.11  0.30 .99                                                                        
##                            [.98, 1.00]                                                                
##                                                                                                       
##   3. VBMT_Men   -0.10 0.30 .99          .96                                                           
##                            [.98, .99]   [.92, .98]                                                    
##                                                                                                       
##   4. D          -0.21 0.09 -.01         -.14         .14                                              
##                            [-.34, .33]  [-.46, .21]  [-.21, .46]                                      
##                                                                                                       
##   5. GEI        0.87  0.07 .26          .35          .17          -.59                                
##                            [-.09, .55]  [.00, .62]   [-.18, .48]  [-.77, -.30]                        
##                                                                                                       
##   6. GGGI       0.74  0.05 .27          .37          .16          -.74         .73                    
##                            [-.07, .56]  [.04, .63]   [-.19, .47]  [-.86, -.53] [.52, .86]             
##                                                                                                       
##   7. GDI        0.98  0.03 -.54         -.50         -.59         -.33         .07         .19        
##                            [-.75, -.25] [-.72, -.19] [-.77, -.32] [-.60, .01]  [-.28, .41] [-.16, .50]
##                                                                                                       
##   8. log_GDP    10.61 0.41 .43          .49          .35          -.51         .72         .62        
##                            [.10, .67]   [.19, .71]   [.01, .61]   [-.72, -.21] [.50, .85]  [.36, .79] 
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
## 1468866.7 1468900.1 -734430.4 1468860.7    492340 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.6809 -0.5301  0.0655  0.6649  5.0161 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  cntry    (Intercept) 0.08982  0.2997  
##  Residual             1.01650  1.0082  
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##              Estimate Std. Error        df t value Pr(>|t|)
## (Intercept)  0.000995   0.051424 33.987671   0.019    0.985
```

``` r
getVC(mod0)
```

```
##        grp        var1 var2 sdcor vcov
## 1    cntry (Intercept) <NA>  0.30 0.09
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
## mean variation  0.08118475     NA       1
## sigma2          0.91881525      1      NA
## 
## $R2s
##          total within between
## f1  0.00000000      0      NA
## f2  0.00000000     NA       0
## v   0.00000000      0      NA
## m   0.08118475     NA       1
## f   0.00000000     NA      NA
## fv  0.00000000      0      NA
## fvm 0.08118475     NA      NA
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
## 1463409.5 1463453.9 -731700.7 1463401.5    492339 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.5220 -0.5231  0.0681  0.6609  5.3025 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  cntry    (Intercept) 0.09116  0.3019  
##  Residual             1.00529  1.0026  
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##               Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept) -3.012e-03  5.181e-02  3.399e+01  -0.058    0.954    
## gndr.c      -2.115e-01  2.855e-03  4.923e+05 -74.092   <2e-16 ***
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
## (Intercept) -0.003 0.052     33.989  -0.058 0.954 -0.108  0.102
## gndr.c      -0.212 0.003 492309.744 -74.092 0.000 -0.217 -0.206
```

``` r
getVC(mod1)
```

```
##        grp        var1 var2 sdcor vcov
## 1    cntry (Intercept) <NA>   0.3 0.09
## 2 Residual        <NA> <NA>   1.0 1.01
```

``` r
r2mlm(mod1,bargraph = F)
```

```
## $Decompositions
##                      total
## fixed           0.01003738
## slope variation 0.00000000
## mean variation  0.08230381
## sigma2          0.90765881
## 
## $R2s
##          total
## f   0.01003738
## v   0.00000000
## m   0.08230381
## fv  0.01003738
## fvm 0.09234119
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
## 1462817.8 1462884.4 -731402.9 1462805.8    492337 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.4745 -0.5144  0.0698  0.6614  5.3247 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr 
##  cntry    (Intercept) 0.09119  0.30198       
##           gndr.c      0.00656  0.08099  0.00 
##  Residual             1.00386  1.00193       
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##              Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept) -0.002897   0.051816 33.988284  -0.056    0.956    
## gndr.c      -0.204267   0.014275 33.021860 -14.310 1.03e-15 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##        (Intr)
## gndr.c 0.000
```

``` r
getFE(mod2,round=3)
```

```
##               Est.    SE     df       t     p     LL     UL
## (Intercept) -0.003 0.052 33.988  -0.056 0.956 -0.108  0.102
## gndr.c      -0.204 0.014 33.022 -14.310 0.000 -0.233 -0.175
```

``` r
getVC(mod2)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.30 0.09
## 2    cntry      gndr.c   <NA>  0.08 0.01
## 3    cntry (Intercept) gndr.c  0.00 0.00
## 4 Residual        <NA>   <NA>  1.00 1.00
```

``` r
r2mlm(mod2,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.009365625
## slope variation 0.001472454
## mean variation  0.082380556
## sigma2          0.906781365
## 
## $R2s
##           total
## f   0.009365625
## v   0.001472454
## m   0.082380556
## fv  0.010838079
## fvm 0.093218635
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
## mod1    4 1463409 1463454 -731701   1463401                         
## mod2    6 1462818 1462884 -731403   1462806 595.65  2  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

``` r
lvl2_var_cond_lvl1(mod2,lvl1.var = "gndr.c",lvl1.values = c(0.5,-0.5))
```

```
##   lvl1.value lvl2.cond.var lvl2.cond.sd
## 1        0.5    0.09283355    0.3046860
## 2       -0.5    0.09282757    0.3046762
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
## 1462815.8 1462871.3 -731402.9 1462805.8    492338 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.4745 -0.5144  0.0698  0.6614  5.3247 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  cntry    (Intercept) 0.09119  0.30198 
##  cntry.1  gndr.c      0.00656  0.08099 
##  Residual             1.00386  1.00193 
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##              Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept) -0.002897   0.051817 33.986205  -0.056    0.956    
## gndr.c      -0.204267   0.014275 33.022126 -14.310 1.03e-15 ***
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
## (Intercept) -0.003 0.052 33.986  -0.056 0.956 -0.108  0.102
## gndr.c      -0.204 0.014 33.022 -14.310 0.000 -0.233 -0.175
```

``` r
getVC(mod2_norecov)
```

```
##        grp        var1 var2 sdcor vcov
## 1    cntry (Intercept) <NA>  0.30 0.09
## 2  cntry.1      gndr.c <NA>  0.08 0.01
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
##              npar     AIC     BIC  logLik -2*log(L) Chisq Df Pr(>Chisq)
## mod2_norecov    5 1462816 1462871 -731403   1462806                    
## mod2            6 1462818 1462884 -731403   1462806     0  1     0.9995
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
## 1417460.3 1417549.0 -708722.2 1417444.3    480356 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.5496 -0.5133  0.0705  0.6633  5.3701 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr 
##  cntry    (Intercept) 0.08083  0.2843        
##           gndr.c      0.00437  0.0661   0.22 
##  Residual             0.98773  0.9938        
## Number of obs: 480364, groups:  cntry, 33
## 
## Fixed effects:
##                 Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)      0.01204    0.04952 32.98163   0.243 0.809391    
## gndr.c          -0.20364    0.01197 32.32283 -17.007  < 2e-16 ***
## gei.z.cm         0.07538    0.05030 33.01736   1.499 0.143482    
## gndr.c:gei.z.cm -0.05072    0.01237 34.50572  -4.100 0.000238 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c g.z.cm
## gndr.c       0.213              
## gei.z.cm     0.000  0.000       
## gndr.c:g.z.  0.000 -0.025  0.209
```

``` r
getFE(mod2_GEI,round=3)
```

```
##                   Est.    SE     df       t     p     LL     UL
## (Intercept)      0.012 0.050 32.982   0.243 0.809 -0.089  0.113
## gndr.c          -0.204 0.012 32.323 -17.007 0.000 -0.228 -0.179
## gei.z.cm         0.075 0.050 33.017   1.499 0.143 -0.027  0.178
## gndr.c:gei.z.cm -0.051 0.012 34.506  -4.100 0.000 -0.076 -0.026
```

``` r
getVC(mod2_GEI)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.28 0.08
## 2    cntry      gndr.c   <NA>  0.07 0.00
## 3    cntry (Intercept) gndr.c  0.22 0.00
## 4 Residual        <NA>   <NA>  0.99 0.99
```

``` r
r2mlm(mod2_GEI,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.014203172
## slope variation 0.001001734
## mean variation  0.074237698
## sigma2          0.910557396
## 
## $R2s
##           total
## f   0.014203172
## v   0.001001734
## m   0.074237698
## fv  0.015204906
## fvm 0.089442604
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
## Time difference of 36.63844 secs
```

``` r
round(ddsc_mod2_GEI$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.093        0.305        1.004     1.097 0.085   7802.647 0.998   0.999
## 2        0.5         0.093        0.305        1.004     1.097 0.085   6678.029 0.998   0.998
```

``` r
round(ddsc_mod2_GEI$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gei.z.cm gei.z.cm_scaled
## means_y1          -0.102 0.313    1.000           1.000    0.963           0.963    0.171           0.171
## means_y1_scaled   -0.329 1.009    1.000           1.000    0.963           0.963    0.171           0.171
## means_y2           0.104 0.308    0.963           0.963    1.000           1.000    0.325           0.325
## means_y2_scaled    0.336 0.991    0.963           0.963    1.000           1.000    0.325           0.325
## gei.z.cm           0.000 1.000    0.171           0.171    0.325           0.325    1.000           1.000
## gei.z.cm_scaled    0.000 1.000    0.171           0.171    0.325           0.325    1.000           1.000
## diff_score        -0.207 0.085    0.200           0.200   -0.073          -0.073   -0.545          -0.545
## diff_score_scaled -0.665 0.274    0.200           0.200   -0.073          -0.073   -0.545          -0.545
##                   diff_score diff_score_scaled
## means_y1               0.200             0.200
## means_y1_scaled        0.200             0.200
## means_y2              -0.073            -0.073
## means_y2_scaled       -0.073            -0.073
## gei.z.cm              -0.545            -0.545
## gei.z.cm_scaled       -0.545            -0.545
## diff_score             1.000             1.000
## diff_score_scaled      1.000             1.000
```

``` r
round(ddsc_mod2_GEI$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                       0.595 0.145 34.506   4.100   0.000    0.300    0.890
## w_11                          0.101 0.049 33.028   2.040   0.049    0.000    0.201
## w_21                          0.050 0.052 33.032   0.963   0.343   -0.056    0.156
## r_xy1                         0.322 0.158 33.028   2.040   0.049    0.001    0.642
## r_xy2                         0.163 0.169 33.032   0.963   0.343   -0.181    0.506
## b_11                          0.324 0.159 33.028   2.040   0.049    0.001    0.648
## b_21                          0.161 0.167 33.032   0.963   0.343   -0.179    0.501
## main_effect                   0.075 0.050 33.017   1.499   0.143   -0.027    0.178
## moderator_effect             -0.204 0.012 32.323 -17.007   0.000   -0.228   -0.179
## interaction                  -0.051 0.012 34.506  -4.100   0.000   -0.076   -0.026
## q_b11_b21                     0.174    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.169    NA     NA      NA      NA       NA       NA
## cross_over_point             -4.015    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.025 0.054 33.078  -0.455   0.652   -0.135    0.086
## interaction_vs_main_bscale   -0.079 0.175 33.078  -0.455   0.652   -0.435    0.276
## interaction_vs_main_rscale   -0.083 0.178 33.076  -0.467   0.643   -0.444    0.278
## dadas                        -0.100 0.104 33.032  -0.963   0.829   -0.311    0.111
## dadas_bscale                 -0.322 0.335 33.032  -0.963   0.829   -1.003    0.359
## dadas_rscale                 -0.325 0.338 33.032  -0.963   0.829   -1.012    0.362
## abs_diff                      0.051 0.012 34.506   4.100   0.000    0.026    0.076
## abs_sum                       0.151 0.101 33.017   1.499   0.072   -0.054    0.355
## abs_diff_bscale               0.163 0.040 34.506   4.100   0.000    0.082    0.244
## abs_sum_bscale                0.486 0.324 33.017   1.499   0.072   -0.174    1.145
## abs_diff_rscale               0.159 0.041 34.525   3.922   0.000    0.077    0.241
## abs_sum_rscale                0.484 0.324 33.017   1.494   0.072   -0.175    1.143
```

``` r
round(ddsc_mod2_GEI$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
##      0      0      0      1      1
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
## r_xy1_y2                         0.545 0.146  3.732  0.000    0.259    0.831
## r_xy1                            0.325 0.165  1.971  0.049    0.002    0.647
## r_xy2                            0.171 0.172  0.995  0.320   -0.166    0.507
## b_11                             0.322 0.163  1.971  0.049    0.002    0.641
## b_21                             0.172 0.173  0.995  0.320   -0.167    0.511
## b_10                             0.336 0.161  2.093  0.036    0.021    0.651
## b_20                            -0.329 0.170 -1.929  0.054   -0.663    0.005
## res_cov_y1_y2                    0.880 0.220  4.007  0.000    0.449    1.310
## diff_b10_b20                     0.665 0.039 16.861  0.000    0.588    0.742
## diff_b11_b21                     0.149 0.040  3.732  0.000    0.071    0.228
## diff_rxy1_rxy2                   0.154 0.039  3.903  0.000    0.077    0.231
## q_b11_b21                        0.160 0.042  3.819  0.000    0.078    0.241
## q_rxy1_rxy2                      0.164 0.042  3.885  0.000    0.081    0.247
## cross_over_point                -4.449 1.221 -3.644  0.000   -6.843   -2.056
## sum_b11_b21                      0.494 0.334  1.479  0.139   -0.161    1.148
## main_effect                      0.247 0.167  1.479  0.139   -0.080    0.574
## interaction_vs_main_effect      -0.097 0.181 -0.538  0.591   -0.452    0.258
## diff_abs_b11_abs_b21             0.149 0.040  3.732  0.000    0.071    0.228
## abs_diff_b11_b21                 0.149 0.040  3.732  0.000    0.071    0.228
## abs_sum_b11_b21                  0.494 0.334  1.479  0.070   -0.161    1.148
## dadas                           -0.344 0.346 -0.995  0.840   -1.023    0.334
## q_r_equivalence                  0.064 0.042  1.522  0.936       NA       NA
## q_b_equivalence                  0.060 0.042  1.426  0.923       NA       NA
## cross_over_point_equivalence     4.449 1.221  3.644  1.000       NA       NA
## cross_over_point_minimal_effect  4.449 1.221  3.644  0.000       NA       NA
```

``` r
round(ddsc_sem_GEI$variance_test,3)
```

```
##              est    se      z pvalue ci.lower ci.upper
## cov_y1y2   0.933 0.234  3.984  0.000    0.474    1.392
## var_y1     0.952 0.234  4.062  0.000    0.493    1.412
## var_y2     0.987 0.243  4.062  0.000    0.511    1.463
## var_diff  -0.035 0.092 -0.375  0.707   -0.215    0.146
## var_ratio  0.965 0.091 10.590  0.000    0.786    1.144
## cor_y1y2   0.963 0.013 75.156  0.000    0.937    0.988
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
## r_xy1y2                       0.545 0.151 31.000   3.617   0.001    0.238    0.852
## w_11                          0.100 0.054 31.892   1.854   0.073   -0.010    0.210
## w_21                          0.053 0.054 31.892   0.992   0.329   -0.056    0.163
## r_xy1                         0.325 0.175 31.892   1.854   0.073   -0.032    0.681
## r_xy2                         0.171 0.172 31.892   0.992   0.329   -0.180    0.521
## b_11                          0.322 0.174 31.892   1.854   0.073   -0.032    0.675
## b_21                          0.172 0.174 31.892   0.992   0.329   -0.181    0.526
## main_effect                   0.077 0.054 31.000   1.433   0.162   -0.032    0.186
## moderator_effect             -0.207 0.013 31.000 -16.342   0.000   -0.232   -0.181
## interaction                  -0.046 0.013 31.000  -3.617   0.001   -0.073   -0.020
## q_b11_b21                     0.160    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.164    NA     NA      NA      NA       NA       NA
## cross_over_point             -4.449    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.030 0.055 34.555  -0.550   0.586   -0.142    0.081
## interaction_vs_main_bscale   -0.097 0.177 34.555  -0.550   0.586   -0.457    0.262
## interaction_vs_main_rscale   -0.094 0.174 34.652  -0.538   0.594   -0.447    0.260
## dadas                        -0.107 0.108 31.892  -0.992   0.836   -0.326    0.113
## dadas_bscale                 -0.344 0.347 31.892  -0.992   0.836   -1.051    0.363
## dadas_rscale                 -0.341 0.344 31.892  -0.992   0.836   -1.042    0.360
## abs_diff                      0.046 0.013 31.000   3.617   0.001    0.020    0.073
## abs_sum                       0.153 0.107 31.000   1.433   0.081   -0.065    0.372
## abs_diff_bscale               0.149 0.041 31.000   3.617   0.001    0.065    0.234
## abs_sum_bscale                0.494 0.345 31.000   1.433   0.081   -0.209    1.197
## abs_diff_rscale               0.154 0.041 31.341   3.713   0.000    0.069    0.238
## abs_sum_rscale                0.495 0.345 31.000   1.437   0.080   -0.208    1.198
```

``` r
# country-time multilevel model


mod2_GEI_cntry_year<-
  lmer(ben.z.wt~gndr.c+gei.z.cm+gndr.c:gei.z.cm+
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
## Formula: ben.z.wt ~ gndr.c + gei.z.cm + gndr.c:gei.z.cm + (gndr.c | cntry)
##    Data: diff_dat_cntry_year
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##    -250.8    -216.6     133.4    -266.8       526 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -5.9459 -0.5916  0.1296  0.6466  3.0504 
## 
## Random effects:
##  Groups   Name        Variance  Std.Dev. Corr 
##  cntry    (Intercept) 0.0740690 0.27216       
##           gndr.c      0.0003152 0.01775  1.00 
##  Residual             0.0282961 0.16821       
## Number of obs: 534, groups:  cntry, 33
## 
## Fixed effects:
##                  Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)       0.01708    0.04808  32.64128   0.355  0.72465    
## gndr.c           -0.20949    0.01527 213.77662 -13.718  < 2e-16 ***
## gei.z.cm          0.07846    0.04926  33.77562   1.593  0.12051    
## gndr.c:gei.z.cm  -0.05024    0.01744 271.40541  -2.880  0.00429 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c g.z.cm
## gndr.c       0.200              
## gei.z.cm    -0.011 -0.001       
## gndr.c:g.z. -0.001 -0.220  0.179
## optimizer (bobyqa) convergence code: 0 (OK)
## boundary (singular) fit: see help('isSingular')
```

``` r
getFE(mod2_GEI_cntry_year,round=3)
```

```
##                   Est.    SE      df       t     p     LL     UL
## (Intercept)      0.017 0.048  32.641   0.355 0.725 -0.081  0.115
## gndr.c          -0.209 0.015 213.777 -13.718 0.000 -0.240 -0.179
## gei.z.cm         0.078 0.049  33.776   1.593 0.121 -0.022  0.179
## gndr.c:gei.z.cm -0.050 0.017 271.405  -2.880 0.004 -0.085 -0.016
```

``` r
getVC(mod2_GEI_cntry_year)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.27 0.07
## 2    cntry      gndr.c   <NA>  0.02 0.00
## 3    cntry (Intercept) gndr.c  1.00 0.00
## 4 Residual        <NA>   <NA>  0.17 0.03
```

``` r
r2mlm(mod2_GEI,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.014203172
## slope variation 0.001001734
## mean variation  0.074237698
## sigma2          0.910557396
## 
## $R2s
##           total
## f   0.014203172
## v   0.001001734
## m   0.074237698
## fv  0.015204906
## fvm 0.089442604
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
## 1       -0.5         0.086        0.294        0.028     0.115 0.752      8.029 0.995    0.96
## 2        0.5         0.085        0.292        0.028     0.113 0.749      8.029 0.995    0.96
```

``` r
round(ddsc_mod2_GEI_cntry_year$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gei.z.cm gei.z.cm_scaled
## means_y1          -0.086 0.297    1.000           1.000    0.956           0.956    0.170           0.170
## means_y1_scaled   -0.291 1.004    1.000           1.000    0.956           0.956    0.170           0.170
## means_y2           0.122 0.294    0.956           0.956    1.000           1.000    0.346           0.346
## means_y2_scaled    0.413 0.996    0.956           0.956    1.000           1.000    0.346           0.346
## gei.z.cm           0.000 1.000    0.170           0.170    0.346           0.346    1.000           1.000
## gei.z.cm_scaled    0.000 1.000    0.170           0.170    0.346           0.346    1.000           1.000
## diff_score        -0.208 0.088    0.176           0.176   -0.121          -0.121   -0.586          -0.586
## diff_score_scaled -0.704 0.297    0.176           0.176   -0.121          -0.121   -0.586          -0.586
##                   diff_score diff_score_scaled
## means_y1               0.176             0.176
## means_y1_scaled        0.176             0.176
## means_y2              -0.121            -0.121
## means_y2_scaled       -0.121            -0.121
## gei.z.cm              -0.586            -0.586
## gei.z.cm_scaled       -0.586            -0.586
## diff_score             1.000             1.000
## diff_score_scaled      1.000             1.000
```

``` r
round(ddsc_mod2_GEI_cntry_year$results,3)
```

```
##                            estimate    SE      df t.ratio p.value ci.lower ci.upper
## r_xy1y2                       0.573 0.199 271.405   2.880   0.004    0.181    0.964
## w_11                          0.104 0.048  33.976   2.137   0.040    0.005    0.202
## w_21                          0.053 0.052  33.944   1.035   0.308   -0.051    0.158
## r_xy1                         0.349 0.163  33.976   2.137   0.040    0.017    0.681
## r_xy2                         0.181 0.175  33.944   1.035   0.308   -0.175    0.537
## b_11                          0.351 0.164  33.976   2.137   0.040    0.017    0.684
## b_21                          0.181 0.174  33.944   1.035   0.308   -0.174    0.535
## main_effect                   0.078 0.049  33.776   1.593   0.121   -0.022    0.179
## moderator_effect             -0.209 0.015 213.777 -13.718   0.000   -0.240   -0.179
## interaction                  -0.050 0.017 271.405  -2.880   0.004   -0.085   -0.016
## q_b11_b21                     0.184    NA      NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.181    NA      NA      NA      NA       NA       NA
## cross_over_point             -4.170    NA      NA      NA      NA       NA       NA
## interaction_vs_main          -0.028 0.055  34.545  -0.512   0.612   -0.140    0.084
## interaction_vs_main_bscale   -0.096 0.187  34.545  -0.512   0.612   -0.474    0.283
## interaction_vs_main_rscale   -0.097 0.188  34.534  -0.518   0.608   -0.479    0.284
## dadas                        -0.107 0.103  33.944  -1.035   0.846   -0.316    0.103
## dadas_bscale                 -0.361 0.349  33.944  -1.035   0.846   -1.070    0.348
## dadas_rscale                 -0.363 0.350  33.944  -1.035   0.846   -1.075    0.349
## abs_diff                      0.050 0.017 271.405   2.880   0.002    0.016    0.085
## abs_sum                       0.157 0.099  33.776   1.593   0.060   -0.043    0.357
## abs_diff_bscale               0.170 0.059 271.405   2.880   0.002    0.054    0.286
## abs_sum_bscale                0.531 0.333  33.776   1.593   0.060   -0.147    1.209
## abs_diff_rscale               0.168 0.059 243.066   2.830   0.003    0.051    0.285
## abs_sum_rscale                0.530 0.334  33.776   1.591   0.061   -0.147    1.208
```

``` r
round(ddsc_mod2_GEI_cntry_year$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.001 -1.000  0.016  1.000  0.898
```

### Comparisons across models


``` r
# full multilevel model (individuals as entries)
round(ddsc_mod2_GEI$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df  t.ratio p.value ci.lower ci.upper
## r_xy1              0.3216 0.1576 33.0277   2.0401  0.0494   0.0009   0.6422
## r_xy2              0.1625 0.1688 33.0319   0.9629  0.3426  -0.1809   0.5059
## b_11               0.3244 0.1590 33.0277   2.0401  0.0494   0.0009   0.6480
## b_21               0.1611 0.1673 33.0319   0.9629  0.3426  -0.1793   0.5014
## main_effect        0.0754 0.0503 33.0174   1.4986  0.1435  -0.0270   0.1777
## moderator_effect  -0.2036 0.0120 32.3228 -17.0074  0.0000  -0.2280  -0.1793
## interaction       -0.0507 0.0124 34.5057  -4.0998  0.0002  -0.0758  -0.0256
## q_b11_b21          0.1741     NA      NA       NA      NA       NA       NA
```

``` r
# country-level structural path model (country-gender means as entries)
round(ddsc_sem_GEI$results[c("r_xy1","r_xy2","b_11","b_21","q_b11_b21","diff_b11_b21"),],4)
```

```
##                 est     se      z pvalue ci.lower ci.upper
## r_xy1        0.3245 0.1647 1.9710 0.0487   0.0018   0.6473
## r_xy2        0.1707 0.1715 0.9949 0.3198  -0.1655   0.5068
## b_11         0.3216 0.1632 1.9710 0.0487   0.0018   0.6415
## b_21         0.1722 0.1730 0.9949 0.3198  -0.1670   0.5113
## q_b11_b21    0.1596 0.0418 3.8190 0.0001   0.0777   0.2415
## diff_b11_b21 0.1495 0.0401 3.7318 0.0002   0.0710   0.2280
```

``` r
# multilevel model (country-gender means as entries)
round(ddsc_mod2_GEI_ri$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df  t.ratio p.value ci.lower ci.upper
## r_xy1              0.3245 0.1751 31.8915   1.8535  0.0731  -0.0322   0.6812
## r_xy2              0.1707 0.1720 31.8915   0.9922  0.3286  -0.1798   0.5211
## b_11               0.3216 0.1735 31.8915   1.8535  0.0731  -0.0319   0.6752
## b_21               0.1722 0.1735 31.8915   0.9922  0.3286  -0.1814   0.5257
## main_effect        0.0767 0.0535 31.0000   1.4330  0.1619  -0.0324   0.1858
## moderator_effect  -0.2065 0.0126 31.0000 -16.3424  0.0000  -0.2323  -0.1807
## interaction       -0.0464 0.0128 31.0000  -3.6169  0.0010  -0.0726  -0.0202
## q_b11_b21          0.1596     NA      NA       NA      NA       NA       NA
```

``` r
# multilevel model (country-year-gender means as entries)
round(ddsc_mod2_GEI_cntry_year$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE       df  t.ratio p.value ci.lower ci.upper
## r_xy1              0.3492 0.1634  33.9757   2.1372  0.0399   0.0171   0.6812
## r_xy2              0.1813 0.1752  33.9439   1.0350  0.3080  -0.1747   0.5373
## b_11               0.3506 0.1640  33.9757   2.1372  0.0399   0.0172   0.6840
## b_21               0.1805 0.1744  33.9439   1.0350  0.3080  -0.1740   0.5351
## main_effect        0.0785 0.0493  33.7756   1.5929  0.1205  -0.0217   0.1786
## moderator_effect  -0.2095 0.0153 213.7766 -13.7177  0.0000  -0.2396  -0.1794
## interaction       -0.0502 0.0174 271.4054  -2.8804  0.0043  -0.0846  -0.0159
## q_b11_b21          0.1836     NA       NA       NA      NA       NA       NA
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
## Time difference of 1.664452 hours
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
##                    Estimate         SE         2.5%      97.5%
## X.Intercept.     0.01079347 0.04881526 -0.087096215  0.1088377
## gndr.c          -0.20318771 0.01188556 -0.226870055 -0.1807550
## gei.z.cm         0.07642455 0.05024175 -0.017880465  0.1800032
## gndr.c.gei.z.cm -0.05114621 0.01287925 -0.076449245 -0.0264830
## w11              0.10199765 0.04903391  0.007556252  0.2009890
## w21              0.05085144 0.05222146 -0.049125270  0.1572154
## b11              0.32848397 0.15791395  0.024334949  0.6472863
## b21              0.16376733 0.16817949 -0.158208189  0.5063131
## r_xy1            0.32558705 0.15652130  0.024120338  0.6415778
## r_xy2            0.16523753 0.16968930 -0.159628488  0.5108585
## q_b              0.18222242 0.04715831  0.094312483  0.2787494
## q                0.17687894 0.04634879  0.090723377  0.2731036
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
## [1] 0.1822224
## 
## $se
## [1] 0.04715831
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
## [1] 5.984575
## 
## $p_low
## [1] 1.084777e-09
## 
## $z_high
## [1] 1.743541
## 
## $p_high
## [1] 0.9593804
## 
## $ci_level
## [1] 0.9
## 
## $ci_lower
## [1] 0.1046539
## 
## $ci_upper
## [1] 0.2597909
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
## [1] 0.1768789
## 
## $se
## [1] 0.04634879
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
## [1] 5.973812
## 
## $p_low
## [1] 1.158864e-09
## 
## $z_high
## [1] 1.658705
## 
## $p_high
## [1] 0.9514123
## 
## $ci_level
## [1] 0.9
## 
## $ci_lower
## [1] 0.100642
## 
## $ci_upper
## [1] 0.2531159
## 
## $equivalent
## [1] FALSE
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
  xlab("Gender Equality Index (Average 2002-2023)")+
  #ylim(c(min.y.pred,max.y.pred))+
  #xlim(c(0.60,1.00))+
  ylab("Value benevolence (Average 2002-2023)")+
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
  xlab("Gender Equality Index (Average 2002-2023)")+
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

![](Analysis_benevolence_files/figure-html/unnamed-chunk-23-1.png)<!-- -->

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
## 1081167.6 1081254.0 -540575.8 1081151.6    363844 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.5507 -0.5120  0.0748  0.6657  5.3019 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr 
##  cntry    (Intercept) 0.094489 0.3074        
##           gndr.c      0.002735 0.0523   0.32 
##  Residual             1.002523 1.0013        
## Number of obs: 363852, groups:  cntry, 34
## 
## Fixed effects:
##                   Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)       0.009966   0.052755 33.964659   0.189   0.8513    
## gndr.c           -0.201032   0.009739 32.931595 -20.643  < 2e-16 ***
## gggi.z.cm         0.105412   0.053564 34.003813   1.968   0.0573 .  
## gndr.c:gggi.z.cm -0.056294   0.010176 36.411620  -5.532 2.84e-06 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c ggg.z.
## gndr.c       0.297              
## gggi.z.cm    0.000  0.000       
## gndr.c:gg..  0.000 -0.019  0.288
```

``` r
getFE(mod2_GGGI,round=3)
```

```
##                    Est.    SE     df       t     p     LL     UL
## (Intercept)       0.010 0.053 33.965   0.189 0.851 -0.097  0.117
## gndr.c           -0.201 0.010 32.932 -20.643 0.000 -0.221 -0.181
## gggi.z.cm         0.105 0.054 34.004   1.968 0.057 -0.003  0.214
## gndr.c:gggi.z.cm -0.056 0.010 36.412  -5.532 0.000 -0.077 -0.036
```

``` r
getVC(mod2_GGGI)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.31 0.09
## 2    cntry      gndr.c   <NA>  0.05 0.00
## 3    cntry (Intercept) gndr.c  0.32 0.01
## 4 Residual        <NA>   <NA>  1.00 1.00
```

``` r
r2mlm(mod2_GGGI,bargraph = F)
```

```
## $Decompositions
##                        total
## fixed           0.0172398286
## slope variation 0.0006084813
## mean variation  0.0842610892
## sigma2          0.8978906009
## 
## $R2s
##            total
## f   0.0172398286
## v   0.0006084813
## m   0.0842610892
## fv  0.0178483099
## fvm 0.1021093991
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
## Time difference of 38.32269 secs
```

``` r
round(ddsc_mod2_GGGI$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.093        0.305        1.004     1.097 0.085   7802.647 0.998   0.999
## 2        0.5         0.093        0.305        1.004     1.097 0.085   6678.029 0.998   0.998
```

``` r
round(ddsc_mod2_GGGI$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gggi.z.cm
## means_y1          -0.103 0.344    1.000           1.000    0.973           0.973     0.223
## means_y1_scaled   -0.299 1.000    1.000           1.000    0.973           0.973     0.223
## means_y2           0.101 0.344    0.973           0.973    1.000           1.000     0.385
## means_y2_scaled    0.294 1.000    0.973           0.973    1.000           1.000     0.385
## gggi.z.cm          0.000 1.000    0.223           0.223    0.385           0.385     1.000
## gggi.z.cm_scaled   0.000 1.000    0.223           0.223    0.385           0.385     1.000
## diff_score        -0.204 0.079    0.119           0.119   -0.112          -0.112    -0.698
## diff_score_scaled -0.593 0.231    0.119           0.119   -0.112          -0.112    -0.698
##                   gggi.z.cm_scaled diff_score diff_score_scaled
## means_y1                     0.223      0.119             0.119
## means_y1_scaled              0.223      0.119             0.119
## means_y2                     0.385     -0.112            -0.112
## means_y2_scaled              0.385     -0.112            -0.112
## gggi.z.cm                    1.000     -0.698            -0.698
## gggi.z.cm_scaled             1.000     -0.698            -0.698
## diff_score                  -0.698      1.000             1.000
## diff_score_scaled           -0.698      1.000             1.000
```

``` r
round(ddsc_mod2_GGGI$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                       0.710 0.128 36.412   5.532   0.000    0.450    0.970
## w_11                          0.134 0.052 34.013   2.552   0.015    0.027    0.240
## w_21                          0.077 0.055 34.005   1.399   0.171   -0.035    0.190
## r_xy1                         0.388 0.152 34.013   2.552   0.015    0.079    0.697
## r_xy2                         0.225 0.161 34.005   1.399   0.171   -0.102    0.551
## b_11                          0.388 0.152 34.013   2.552   0.015    0.079    0.698
## b_21                          0.225 0.161 34.005   1.399   0.171   -0.102    0.551
## main_effect                   0.105 0.054 34.004   1.968   0.057   -0.003    0.214
## moderator_effect             -0.201 0.010 32.932 -20.643   0.000   -0.221   -0.181
## interaction                  -0.056 0.010 36.412  -5.532   0.000   -0.077   -0.036
## q_b11_b21                     0.181    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.181    NA     NA      NA      NA       NA       NA
## cross_over_point             -3.571    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.049 0.057 34.024  -0.857   0.398   -0.166    0.067
## interaction_vs_main_bscale   -0.143 0.167 34.024  -0.857   0.398   -0.482    0.196
## interaction_vs_main_rscale   -0.143 0.167 34.024  -0.857   0.397   -0.482    0.196
## dadas                        -0.155 0.110 34.005  -1.399   0.915   -0.379    0.070
## dadas_bscale                 -0.449 0.321 34.005  -1.399   0.915   -1.102    0.204
## dadas_rscale                 -0.450 0.321 34.005  -1.399   0.915   -1.103    0.204
## abs_diff                      0.056 0.010 36.412   5.532   0.000    0.036    0.077
## abs_sum                       0.211 0.107 34.004   1.968   0.029   -0.007    0.429
## abs_diff_bscale               0.164 0.030 36.412   5.532   0.000    0.104    0.224
## abs_sum_bscale                0.613 0.312 34.004   1.968   0.029   -0.020    1.246
## abs_diff_rscale               0.163 0.030 36.411   5.517   0.000    0.103    0.224
## abs_sum_rscale                0.613 0.312 34.004   1.968   0.029   -0.020    1.246
```

``` r
round(ddsc_mod2_GGGI$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
##      0      0      0      1      1
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
## r_xy1_y2                         0.698 0.123  5.689  0.000    0.458    0.939
## r_xy1                            0.385 0.158  2.430  0.015    0.074    0.695
## r_xy2                            0.223 0.167  1.337  0.181   -0.104    0.551
## b_11                             0.385 0.158  2.430  0.015    0.074    0.695
## b_21                             0.224 0.167  1.337  0.181   -0.104    0.551
## b_10                             0.294 0.156  1.886  0.059   -0.012    0.600
## b_20                            -0.299 0.165 -1.813  0.070   -0.622    0.024
## res_cov_y1_y2                    0.861 0.210  4.095  0.000    0.449    1.274
## diff_b10_b20                     0.593 0.028 21.259  0.000    0.538    0.647
## diff_b11_b21                     0.161 0.028  5.689  0.000    0.106    0.216
## diff_rxy1_rxy2                   0.161 0.028  5.706  0.000    0.106    0.217
## q_b11_b21                        0.178 0.031  5.677  0.000    0.117    0.239
## q_rxy1_rxy2                      0.178 0.031  5.675  0.000    0.117    0.240
## cross_over_point                -3.681 0.670 -5.496  0.000   -4.994   -2.368
## sum_b11_b21                      0.608 0.324  1.875  0.061   -0.028    1.244
## main_effect                      0.304 0.162  1.875  0.061   -0.014    0.622
## interaction_vs_main_effect      -0.143 0.173 -0.826  0.409   -0.483    0.197
## diff_abs_b11_abs_b21             0.161 0.028  5.689  0.000    0.106    0.216
## abs_diff_b11_b21                 0.161 0.028  5.689  0.000    0.106    0.216
## abs_sum_b11_b21                  0.608 0.324  1.875  0.030   -0.028    1.244
## dadas                           -0.447 0.334 -1.337  0.909   -1.103    0.208
## q_r_equivalence                  0.078 0.031  2.491  0.994       NA       NA
## q_b_equivalence                  0.078 0.031  2.487  0.994       NA       NA
## cross_over_point_equivalence     3.681 0.670  5.496  1.000       NA       NA
## cross_over_point_minimal_effect  3.681 0.670  5.496  0.000       NA       NA
```

``` r
round(ddsc_sem_GGGI$variance_test,3)
```

```
##              est    se       z pvalue ci.lower ci.upper
## cov_y1y2   0.945 0.232   4.067  0.000    0.490    1.400
## var_y1     0.970 0.235   4.123  0.000    0.509    1.431
## var_y2     0.971 0.236   4.123  0.000    0.510    1.433
## var_diff  -0.002 0.076  -0.020  0.984   -0.151    0.148
## var_ratio  0.998 0.078  12.733  0.000    0.845    1.152
## cor_y1y2   0.973 0.009 108.257  0.000    0.956    0.991
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
## r_xy1y2                       0.698 0.127 32.000   5.519   0.000    0.441    0.956
## w_11                          0.132 0.058 32.487   2.292   0.029    0.015    0.250
## w_21                          0.077 0.058 32.487   1.332   0.192   -0.041    0.194
## r_xy1                         0.385 0.168 32.487   2.292   0.029    0.043    0.726
## r_xy2                         0.223 0.168 32.487   1.332   0.192   -0.118    0.565
## b_11                          0.385 0.168 32.487   2.292   0.029    0.043    0.726
## b_21                          0.224 0.168 32.487   1.332   0.192   -0.118    0.565
## main_effect                   0.105 0.057 32.000   1.819   0.078   -0.013    0.222
## moderator_effect             -0.204 0.010 32.000 -20.624   0.000   -0.224   -0.184
## interaction                  -0.055 0.010 32.000  -5.519   0.000   -0.076   -0.035
## q_b11_b21                     0.178    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.178    NA     NA      NA      NA       NA       NA
## cross_over_point             -3.681    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.049 0.058 33.946  -0.843   0.405   -0.168    0.069
## interaction_vs_main_bscale   -0.143 0.170 33.946  -0.843   0.405   -0.488    0.202
## interaction_vs_main_rscale   -0.143 0.170 33.949  -0.843   0.405   -0.488    0.202
## dadas                        -0.154 0.115 32.487  -1.332   0.904   -0.389    0.081
## dadas_bscale                 -0.447 0.336 32.487  -1.332   0.904   -1.130    0.236
## dadas_rscale                 -0.447 0.335 32.487  -1.332   0.904   -1.130    0.236
## abs_diff                      0.055 0.010 32.000   5.519   0.000    0.035    0.076
## abs_sum                       0.209 0.115 32.000   1.819   0.039   -0.025    0.443
## abs_diff_bscale               0.161 0.029 32.000   5.519   0.000    0.102    0.220
## abs_sum_bscale                0.608 0.334 32.000   1.819   0.039   -0.073    1.289
## abs_diff_rscale               0.161 0.029 32.001   5.527   0.000    0.102    0.221
## abs_sum_rscale                0.608 0.334 32.000   1.819   0.039   -0.073    1.289
```

``` r
# country-time multilevel model


mod2_GGGI_cntry_year<-
  lmer(ben.z.wt~gndr.c+gggi.z.cm+gndr.c:gggi.z.cm+
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
## Formula: ben.z.wt ~ gndr.c + gggi.z.cm + gndr.c:gggi.z.cm + (gndr.c |      cntry)
##    Data: diff_dat_cntry_year
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##    -261.7    -229.8     138.9    -277.7       392 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -5.1109 -0.5400  0.0630  0.5796  4.5582 
## 
## Random effects:
##  Groups   Name        Variance  Std.Dev. Corr 
##  cntry    (Intercept) 0.0875934 0.29596       
##           gndr.c      0.0002899 0.01703  1.00 
##  Residual             0.0211657 0.14548       
## Number of obs: 400, groups:  cntry, 34
## 
## Fixed effects:
##                   Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)        0.01349    0.05146  33.56631   0.262  0.79482    
## gndr.c            -0.20574    0.01499 195.51562 -13.724  < 2e-16 ***
## gggi.z.cm          0.10490    0.05248  34.17835   1.999  0.05364 .  
## gndr.c:gggi.z.cm  -0.05205    0.01607 214.72803  -3.238  0.00139 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c ggg.z.
## gndr.c       0.193              
## gggi.z.cm   -0.008 -0.001       
## gndr.c:gg.. -0.001 -0.139  0.184
## optimizer (bobyqa) convergence code: 0 (OK)
## boundary (singular) fit: see help('isSingular')
```

``` r
getFE(mod2_GGGI_cntry_year,round=3)
```

```
##                    Est.    SE      df       t     p     LL     UL
## (Intercept)       0.013 0.051  33.566   0.262 0.795 -0.091  0.118
## gndr.c           -0.206 0.015 195.516 -13.724 0.000 -0.235 -0.176
## gggi.z.cm         0.105 0.052  34.178   1.999 0.054 -0.002  0.212
## gndr.c:gggi.z.cm -0.052 0.016 214.728  -3.238 0.001 -0.084 -0.020
```

``` r
getVC(mod2_GGGI_cntry_year)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.30 0.09
## 2    cntry      gndr.c   <NA>  0.02 0.00
## 3    cntry (Intercept) gndr.c  1.00 0.01
## 4 Residual        <NA>   <NA>  0.15 0.02
```

``` r
r2mlm(mod2_GGGI,bargraph = F)
```

```
## $Decompositions
##                        total
## fixed           0.0172398286
## slope variation 0.0006084813
## mean variation  0.0842610892
## sigma2          0.8978906009
## 
## $R2s
##            total
## f   0.0172398286
## v   0.0006084813
## m   0.0842610892
## fv  0.0178483099
## fvm 0.1021093991
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
## 1       -0.5         0.086        0.294        0.028     0.115 0.752      8.029 0.995    0.96
## 2        0.5         0.085        0.292        0.028     0.113 0.749      8.029 0.995    0.96
```

``` r
round(ddsc_mod2_GGGI_cntry_year$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gggi.z.cm
## means_y1          -0.087 0.324    1.000           1.000    0.968           0.968     0.227
## means_y1_scaled   -0.265 0.994    1.000           1.000    0.968           0.968     0.227
## means_y2           0.115 0.328    0.968           0.968    1.000           1.000     0.406
## means_y2_scaled    0.354 1.006    0.968           0.968    1.000           1.000     0.406
## gggi.z.cm          0.000 1.000    0.227           0.227    0.406           0.406     1.000
## gggi.z.cm_scaled   0.000 1.000    0.227           0.227    0.406           0.406     1.000
## diff_score        -0.202 0.083    0.079           0.079   -0.176          -0.176    -0.719
## diff_score_scaled -0.619 0.255    0.079           0.079   -0.176          -0.176    -0.719
##                   gggi.z.cm_scaled diff_score diff_score_scaled
## means_y1                     0.227      0.079             0.079
## means_y1_scaled              0.227      0.079             0.079
## means_y2                     0.406     -0.176            -0.176
## means_y2_scaled              0.406     -0.176            -0.176
## gggi.z.cm                    1.000     -0.719            -0.719
## gggi.z.cm_scaled             1.000     -0.719            -0.719
## diff_score                  -0.719      1.000             1.000
## diff_score_scaled           -0.719      1.000             1.000
```

``` r
round(ddsc_mod2_GGGI_cntry_year$results,3)
```

```
##                            estimate    SE      df t.ratio p.value ci.lower ci.upper
## r_xy1y2                       0.626 0.193 214.728   3.238   0.001    0.245    1.007
## w_11                          0.131 0.052  34.120   2.536   0.016    0.026    0.236
## w_21                          0.079 0.055  34.145   1.446   0.157   -0.032    0.190
## r_xy1                         0.404 0.159  34.120   2.536   0.016    0.080    0.727
## r_xy2                         0.240 0.166  34.145   1.446   0.157   -0.097    0.578
## b_11                          0.401 0.158  34.120   2.536   0.016    0.080    0.723
## b_21                          0.242 0.167  34.145   1.446   0.157   -0.098    0.581
## main_effect                   0.105 0.052  34.178   1.999   0.054   -0.002    0.212
## moderator_effect             -0.206 0.015 195.516 -13.724   0.000   -0.235   -0.176
## interaction                  -0.052 0.016 214.728  -3.238   0.001   -0.084   -0.020
## q_b11_b21                     0.179    NA      NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.183    NA      NA      NA      NA       NA       NA
## cross_over_point             -3.953    NA      NA      NA      NA       NA       NA
## interaction_vs_main          -0.053 0.058  34.167  -0.917   0.366   -0.170    0.064
## interaction_vs_main_bscale   -0.162 0.177  34.167  -0.917   0.366   -0.521    0.197
## interaction_vs_main_rscale   -0.158 0.175  34.170  -0.907   0.371   -0.513    0.196
## dadas                        -0.158 0.109  34.145  -1.446   0.921   -0.379    0.064
## dadas_bscale                 -0.484 0.334  34.145  -1.446   0.921   -1.163    0.196
## dadas_rscale                 -0.481 0.332  34.145  -1.446   0.921   -1.156    0.195
## abs_diff                      0.052 0.016 214.728   3.238   0.001    0.020    0.084
## abs_sum                       0.210 0.105  34.178   1.999   0.027   -0.003    0.423
## abs_diff_bscale               0.160 0.049 214.728   3.238   0.001    0.062    0.257
## abs_sum_bscale                0.643 0.322  34.178   1.999   0.027   -0.011    1.297
## abs_diff_rscale               0.164 0.049 254.182   3.342   0.000    0.067    0.260
## abs_sum_rscale                0.644 0.322  34.178   2.002   0.027   -0.010    1.298
```

``` r
round(ddsc_mod2_GGGI_cntry_year$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.001 -1.000  0.016  1.000  0.898
```

### Comparisons across models


``` r
# full multilevel model (individuals as entries)
round(ddsc_mod2_GGGI$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df  t.ratio p.value ci.lower ci.upper
## r_xy1              0.3882 0.1521 34.0127   2.5524  0.0154   0.0791   0.6973
## r_xy2              0.2248 0.1607 34.0047   1.3986  0.1710  -0.1018   0.5513
## b_11               0.3884 0.1522 34.0127   2.5524  0.0154   0.0791   0.6976
## b_21               0.2247 0.1606 34.0047   1.3986  0.1710  -0.1018   0.5511
## main_effect        0.1054 0.0536 34.0038   1.9679  0.0573  -0.0034   0.2143
## moderator_effect  -0.2010 0.0097 32.9316 -20.6428  0.0000  -0.2208  -0.1812
## interaction       -0.0563 0.0102 36.4116  -5.5320  0.0000  -0.0769  -0.0357
## q_b11_b21          0.1813     NA      NA       NA      NA       NA       NA
```

``` r
# country-level structural path model (country-gender means as entries)
round(ddsc_sem_GGGI$results[c("r_xy1","r_xy2","b_11","b_21","q_b11_b21","diff_b11_b21"),],4)
```

```
##                 est     se      z pvalue ci.lower ci.upper
## r_xy1        0.3847 0.1583 2.4303 0.0151   0.0745   0.6950
## r_xy2        0.2235 0.1672 1.3370 0.1812  -0.1041   0.5511
## b_11         0.3846 0.1582 2.4303 0.0151   0.0744   0.6947
## b_21         0.2236 0.1672 1.3370 0.1812  -0.1042   0.5513
## q_b11_b21    0.1780 0.0314 5.6770 0.0000   0.1165   0.2394
## diff_b11_b21 0.1610 0.0283 5.6891 0.0000   0.1055   0.2164
```

``` r
# multilevel model (country-gender means as entries)
round(ddsc_mod2_GGGI_ri$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE     df  t.ratio p.value ci.lower ci.upper
## r_xy1              0.3847 0.1679 32.487   2.2917  0.0285   0.0430   0.7265
## r_xy2              0.2235 0.1677 32.487   1.3324  0.1920  -0.1180   0.5650
## b_11               0.3846 0.1678 32.487   2.2917  0.0285   0.0430   0.7262
## b_21               0.2236 0.1678 32.487   1.3324  0.1920  -0.1180   0.5652
## main_effect        0.1046 0.0575 32.000   1.8190  0.0783  -0.0125   0.2217
## moderator_effect  -0.2038 0.0099 32.000 -20.6241  0.0000  -0.2239  -0.1837
## interaction       -0.0554 0.0100 32.000  -5.5192  0.0000  -0.0758  -0.0349
## q_b11_b21          0.1780     NA     NA       NA      NA       NA       NA
```

``` r
# multilevel model (country-year-gender means as entries)
round(ddsc_mod2_GGGI_cntry_year$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE       df  t.ratio p.value ci.lower ci.upper
## r_xy1              0.4039 0.1592  34.1197   2.5365  0.0159   0.0803   0.7274
## r_xy2              0.2403 0.1661  34.1447   1.4464  0.1572  -0.0973   0.5778
## b_11               0.4013 0.1582  34.1197   2.5365  0.0159   0.0798   0.7228
## b_21               0.2418 0.1672  34.1447   1.4464  0.1572  -0.0979   0.5815
## main_effect        0.1049 0.0525  34.1784   1.9987  0.0536  -0.0017   0.2115
## moderator_effect  -0.2057 0.0150 195.5156 -13.7239  0.0000  -0.2353  -0.1762
## interaction       -0.0520 0.0161 214.7280  -3.2381  0.0014  -0.0837  -0.0204
## q_b11_b21          0.1786     NA       NA       NA      NA       NA       NA
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
## Time difference of 1.267411 hours
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
##                      Estimate         SE        2.5%      97.5%
## X.Intercept.      0.008768654 0.05321263 -0.08814490  0.1137347
## gndr.c           -0.200601836 0.01026274 -0.22096643 -0.1817461
## gggi.z.cm         0.106986129 0.05198933 -0.00175850  0.2071856
## gndr.c.gggi.z.cm -0.056027556 0.01061425 -0.07726400 -0.0345303
## w11               0.134999908 0.05080695  0.03034573  0.2332079
## w21               0.078972351 0.05367277 -0.03278533  0.1813533
## b11               0.392545682 0.14773381  0.08823774  0.6781097
## b21               0.229631679 0.15606688 -0.09533148  0.5273296
## r_xy1             0.392392185 0.14767604  0.08820324  0.6778445
## r_xy2             0.229721543 0.15612795 -0.09536879  0.5275360
## q_b               0.187666786 0.03919845  0.11311266  0.2740097
## q                 0.187359961 0.03911585  0.11290197  0.2736032
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
## [1] 0.1876668
## 
## $se
## [1] 0.03919845
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
## [1] 7.338728
## 
## $p_low
## [1] 1.078027e-13
## 
## $z_high
## [1] 2.236486
## 
## $p_high
## [1] 0.98734
## 
## $ci_level
## [1] 0.9
## 
## $ci_lower
## [1] 0.1231911
## 
## $ci_upper
## [1] 0.2521425
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
## [1] 0.18736
## 
## $se
## [1] 0.03911585
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
## [1] 7.346382
## 
## $p_low
## [1] 1.018075e-13
## 
## $z_high
## [1] 2.233365
## 
## $p_high
## [1] 0.9872376
## 
## $ci_level
## [1] 0.9
## 
## $ci_lower
## [1] 0.1230201
## 
## $ci_upper
## [1] 0.2516998
## 
## $equivalent
## [1] FALSE
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
  xlab("Global Gender Gap Index (Average 2002-2023)")+
  #ylim(c(min.y.pred,max.y.pred))+
  #xlim(c(0.60,1.00))+
  ylab("Value benevolence (Average 2002-2023)")+
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
  xlab("Global Gender Gap Index (Average 2002-2023)")+
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

![](Analysis_benevolence_files/figure-html/unnamed-chunk-29-1.png)<!-- -->

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
## 1462805.6 1462894.4 -731394.8 1462789.6    492335 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.4766 -0.5147  0.0698  0.6614  5.3250 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr  
##  cntry    (Intercept) 0.066020 0.25694        
##           gndr.c      0.005977 0.07731  -0.22 
##  Residual             1.003860 1.00193        
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##                  Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)     -0.002889   0.044097 33.964087  -0.066  0.94815    
## gndr.c          -0.204065   0.013659 33.888597 -14.939  < 2e-16 ***
## gdi.z.cm        -0.161066   0.044778 34.019299  -3.597  0.00101 ** 
## gndr.c:gdi.z.cm -0.027875   0.014089 36.045190  -1.978  0.05556 .  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c gd.z.c
## gndr.c      -0.212              
## gdi.z.cm     0.000  0.000       
## gndr.c:gd..  0.000 -0.008 -0.209
```

``` r
getFE(mod2_GDI,round=3)
```

```
##                   Est.    SE     df       t     p     LL     UL
## (Intercept)     -0.003 0.044 33.964  -0.066 0.948 -0.093  0.087
## gndr.c          -0.204 0.014 33.889 -14.939 0.000 -0.232 -0.176
## gdi.z.cm        -0.161 0.045 34.019  -3.597 0.001 -0.252 -0.070
## gndr.c:gdi.z.cm -0.028 0.014 36.045  -1.978 0.056 -0.056  0.001
```

``` r
getVC(mod2_GDI)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.26 0.07
## 2    cntry      gndr.c   <NA>  0.08 0.01
## 3    cntry (Intercept) gndr.c -0.22 0.00
## 4 Residual        <NA>   <NA>  1.00 1.00
```

``` r
r2mlm(mod2_GDI,bargraph = F)
```

```
## $Decompositions
##                      total
## fixed           0.02479469
## slope variation 0.00135140
## mean variation  0.06039069
## sigma2          0.91346322
## 
## $R2s
##          total
## f   0.02479469
## v   0.00135140
## m   0.06039069
## fv  0.02614609
## fvm 0.08653678
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
## Time difference of 38.34109 secs
```

``` r
round(ddsc_mod2_GDI$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.093        0.305        1.004     1.097 0.085   7802.647 0.998   0.999
## 2        0.5         0.093        0.305        1.004     1.097 0.085   6678.029 0.998   0.998
```

``` r
round(ddsc_mod2_GDI$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gdi.z.cm gdi.z.cm_scaled
## means_y1          -0.117 0.321    1.000           1.000    0.965           0.965   -0.551          -0.551
## means_y1_scaled   -0.368 1.008    1.000           1.000    0.965           0.965   -0.551          -0.551
## means_y2           0.089 0.315    0.965           0.965    1.000           1.000   -0.461          -0.461
## means_y2_scaled    0.281 0.992    0.965           0.965    1.000           1.000   -0.461          -0.461
## gdi.z.cm           0.000 1.000   -0.551          -0.551   -0.461          -0.461    1.000           1.000
## gdi.z.cm_scaled    0.000 1.000   -0.551          -0.551   -0.461          -0.461    1.000           1.000
## diff_score        -0.206 0.084    0.192           0.192   -0.071          -0.071   -0.372          -0.372
## diff_score_scaled -0.649 0.264    0.192           0.192   -0.071          -0.071   -0.372          -0.372
##                   diff_score diff_score_scaled
## means_y1               0.192             0.192
## means_y1_scaled        0.192             0.192
## means_y2              -0.071            -0.071
## means_y2_scaled       -0.071            -0.071
## gdi.z.cm              -0.372            -0.372
## gdi.z.cm_scaled       -0.372            -0.372
## diff_score             1.000             1.000
## diff_score_scaled      1.000             1.000
```

``` r
round(ddsc_mod2_GDI$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                       0.332 0.168 36.045   1.978   0.056   -0.008    0.673
## w_11                         -0.147 0.047 34.021  -3.147   0.003   -0.242   -0.052
## w_21                         -0.175 0.044 34.060  -3.991   0.000   -0.264   -0.086
## r_xy1                        -0.459 0.146 34.021  -3.147   0.003   -0.755   -0.163
## r_xy2                        -0.555 0.139 34.060  -3.991   0.000   -0.838   -0.272
## b_11                         -0.463 0.147 34.021  -3.147   0.003   -0.762   -0.164
## b_21                         -0.550 0.138 34.060  -3.991   0.000   -0.831   -0.270
## main_effect                  -0.161 0.045 34.019  -3.597   0.001   -0.252   -0.070
## moderator_effect             -0.204 0.014 33.889 -14.939   0.000   -0.232   -0.176
## interaction                  -0.028 0.014 36.045  -1.978   0.056   -0.056    0.001
## q_b11_b21                     0.118    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.129    NA     NA      NA      NA       NA       NA
## cross_over_point             -7.321    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.133 0.050 34.039  -2.682   0.011   -0.234   -0.032
## interaction_vs_main_bscale   -0.419 0.156 34.039  -2.682   0.011   -0.736   -0.101
## interaction_vs_main_rscale   -0.411 0.154 34.040  -2.670   0.012   -0.724   -0.098
## dadas                        -0.294 0.094 34.021  -3.147   0.998   -0.484   -0.104
## dadas_bscale                 -0.925 0.294 34.021  -3.147   0.998   -1.523   -0.328
## dadas_rscale                 -0.918 0.292 34.021  -3.147   0.998   -1.511   -0.325
## abs_diff                      0.028 0.014 36.045   1.978   0.028   -0.001    0.056
## abs_sum                       0.322 0.090 34.019   3.597   0.001    0.140    0.504
## abs_diff_bscale               0.088 0.044 36.045   1.978   0.028   -0.002    0.178
## abs_sum_bscale                1.013 0.282 34.019   3.597   0.001    0.441    1.586
## abs_diff_rscale               0.096 0.044 36.205   2.186   0.018    0.007    0.185
## abs_sum_rscale                1.014 0.282 34.019   3.600   0.001    0.442    1.586
```

``` r
round(ddsc_mod2_GDI$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
##      0      0      0      1      1
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
## r_xy1_y2                         0.372 0.159  2.340  0.019    0.061    0.684
## r_xy1                           -0.461 0.152 -3.027  0.002   -0.759   -0.162
## r_xy2                           -0.551 0.143 -3.847  0.000   -0.831   -0.270
## b_11                            -0.457 0.151 -3.027  0.002   -0.753   -0.161
## b_21                            -0.555 0.144 -3.847  0.000   -0.838   -0.272
## b_10                             0.281 0.149  1.892  0.058   -0.010    0.573
## b_20                            -0.368 0.142 -2.589  0.010   -0.647   -0.089
## res_cov_y1_y2                    0.691 0.171  4.039  0.000    0.355    1.026
## diff_b10_b20                     0.649 0.041 15.696  0.000    0.568    0.730
## diff_b11_b21                     0.098 0.042  2.340  0.019    0.016    0.181
## diff_rxy1_rxy2                   0.090 0.042  2.121  0.034    0.007    0.173
## q_b11_b21                        0.132 0.059  2.252  0.024    0.017    0.248
## q_rxy1_rxy2                      0.121 0.057  2.116  0.034    0.009    0.234
## cross_over_point                -6.607 2.855 -2.315  0.021  -12.202   -1.013
## sum_b11_b21                     -1.012 0.292 -3.462  0.001   -1.585   -0.439
## main_effect                     -0.506 0.146 -3.462  0.001   -0.793   -0.220
## interaction_vs_main_effect      -0.408 0.158 -2.574  0.010   -0.718   -0.097
## diff_abs_b11_abs_b21            -0.098 0.042 -2.340  0.019   -0.181   -0.016
## abs_diff_b11_b21                 0.098 0.042  2.340  0.010    0.016    0.181
## abs_sum_b11_b21                  1.012 0.292  3.462  0.000    0.439    1.585
## dadas                           -0.914 0.302 -3.027  0.999   -1.506   -0.322
## q_r_equivalence                  0.021 0.057  0.370  0.644       NA       NA
## q_b_equivalence                  0.032 0.059  0.552  0.710       NA       NA
## cross_over_point_equivalence     6.607 2.855  2.315  0.990       NA       NA
## cross_over_point_minimal_effect  6.607 2.855  2.315  0.010       NA       NA
```

``` r
round(ddsc_sem_GDI$variance_test,3)
```

```
##              est    se      z pvalue ci.lower ci.upper
## cov_y1y2   0.937 0.231  4.050  0.000    0.483    1.390
## var_y1     0.955 0.232  4.123  0.000    0.501    1.409
## var_y2     0.986 0.239  4.123  0.000    0.518    1.455
## var_diff  -0.032 0.087 -0.364  0.716   -0.203    0.139
## var_ratio  0.968 0.087 11.167  0.000    0.798    1.138
## cor_y1y2   0.965 0.012 82.576  0.000    0.942    0.988
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
## r_xy1y2                       0.372 0.164 32.000   2.270   0.030    0.038    0.707
## w_11                         -0.145 0.048 33.320  -3.002   0.005   -0.244   -0.047
## w_21                         -0.177 0.048 33.320  -3.647   0.001   -0.275   -0.078
## r_xy1                        -0.461 0.153 33.320  -3.002   0.005   -0.773   -0.149
## r_xy2                        -0.551 0.151 33.320  -3.647   0.001   -0.858   -0.244
## b_11                         -0.457 0.152 33.320  -3.002   0.005   -0.767   -0.147
## b_21                         -0.555 0.152 33.320  -3.647   0.001   -0.865   -0.246
## main_effect                  -0.161 0.048 32.000  -3.359   0.002   -0.258   -0.063
## moderator_effect             -0.206 0.014 32.000 -15.227   0.000   -0.234   -0.179
## interaction                  -0.031 0.014 32.000  -2.270   0.030   -0.059   -0.003
## q_b11_b21                     0.132    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.121    NA     NA      NA      NA       NA       NA
## cross_over_point             -6.607    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.130 0.050 37.247  -2.601   0.013   -0.231   -0.029
## interaction_vs_main_bscale   -0.408 0.157 37.247  -2.601   0.013   -0.725   -0.090
## interaction_vs_main_rscale   -0.416 0.159 37.123  -2.611   0.013   -0.738   -0.093
## dadas                        -0.291 0.097 33.320  -3.002   0.997   -0.487   -0.094
## dadas_bscale                 -0.914 0.304 33.320  -3.002   0.997   -1.533   -0.295
## dadas_rscale                 -0.921 0.307 33.320  -3.002   0.997   -1.546   -0.297
## abs_diff                      0.031 0.014 32.000   2.270   0.015    0.003    0.059
## abs_sum                       0.322 0.096 32.000   3.359   0.001    0.127    0.517
## abs_diff_bscale               0.098 0.043 32.000   2.270   0.015    0.010    0.186
## abs_sum_bscale                1.012 0.301 32.000   3.359   0.001    0.398    1.626
## abs_diff_rscale               0.090 0.043 32.207   2.076   0.023    0.002    0.178
## abs_sum_rscale                1.011 0.301 32.000   3.356   0.001    0.398    1.625
```

``` r
# country-time multilevel model


mod2_GDI_cntry_year<-
  lmer(ben.z.wt~gndr.c+gdi.z.cm+gndr.c:gdi.z.cm+
                 (gndr.c|cntry),data=diff_dat_cntry_year,REML=F,
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
```

```
## boundary (singular) fit: see help('isSingular')
```

``` r
summary(mod2_GDI_cntry_year)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ben.z.wt ~ gndr.c + gdi.z.cm + gndr.c:gdi.z.cm + (gndr.c | cntry)
##    Data: diff_dat_cntry_year
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##    -261.5    -227.1     138.8    -277.5       538 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -5.9430 -0.5549  0.0986  0.6394  3.0316 
## 
## Random effects:
##  Groups   Name        Variance  Std.Dev. Corr  
##  cntry    (Intercept) 5.829e-02 0.241424       
##           gndr.c      8.695e-05 0.009324 -1.00 
##  Residual             2.843e-02 0.168622       
## Number of obs: 546, groups:  cntry, 34
## 
## Fixed effects:
##                   Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)       0.003625   0.042182  33.434784   0.086 0.932025    
## gndr.c           -0.217433   0.014545 358.437378 -14.949  < 2e-16 ***
## gdi.z.cm         -0.168584   0.043382  35.175298  -3.886 0.000431 ***
## gndr.c:gdi.z.cm  -0.012080   0.017800 425.502010  -0.679 0.497722    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c gd.z.c
## gndr.c      -0.109              
## gdi.z.cm    -0.007  0.001       
## gndr.c:gd..  0.000 -0.056 -0.090
## optimizer (bobyqa) convergence code: 0 (OK)
## boundary (singular) fit: see help('isSingular')
```

``` r
getFE(mod2_GDI_cntry_year,round=3)
```

```
##                   Est.    SE      df       t     p     LL     UL
## (Intercept)      0.004 0.042  33.435   0.086 0.932 -0.082  0.089
## gndr.c          -0.217 0.015 358.437 -14.949 0.000 -0.246 -0.189
## gdi.z.cm        -0.169 0.043  35.175  -3.886 0.000 -0.257 -0.081
## gndr.c:gdi.z.cm -0.012 0.018 425.502  -0.679 0.498 -0.047  0.023
```

``` r
getVC(mod2_GDI_cntry_year)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.24 0.06
## 2    cntry      gndr.c   <NA>  0.01 0.00
## 3    cntry (Intercept) gndr.c -1.00 0.00
## 4 Residual        <NA>   <NA>  0.17 0.03
```

``` r
r2mlm(mod2_GDI,bargraph = F)
```

```
## $Decompositions
##                      total
## fixed           0.02479469
## slope variation 0.00135140
## mean variation  0.06039069
## sigma2          0.91346322
## 
## $R2s
##          total
## f   0.02479469
## v   0.00135140
## m   0.06039069
## fv  0.02614609
## fvm 0.08653678
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
## 1       -0.5         0.086        0.294        0.028     0.115 0.752      8.029 0.995    0.96
## 2        0.5         0.085        0.292        0.028     0.113 0.749      8.029 0.995    0.96
```

``` r
round(ddsc_mod2_GDI_cntry_year$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gdi.z.cm gdi.z.cm_scaled
## means_y1          -0.101 0.304    1.000           1.000    0.959           0.959   -0.591          -0.591
## means_y1_scaled   -0.331 1.000    1.000           1.000    0.959           0.959   -0.591          -0.591
## means_y2           0.106 0.304    0.959           0.959    1.000           1.000   -0.497          -0.497
## means_y2_scaled    0.348 1.000    0.959           0.959    1.000           1.000   -0.497          -0.497
## gdi.z.cm           0.000 1.000   -0.591          -0.591   -0.497          -0.497    1.000           1.000
## gdi.z.cm_scaled    0.000 1.000   -0.591          -0.591   -0.497          -0.497    1.000           1.000
## diff_score        -0.207 0.087    0.142           0.142   -0.143          -0.143   -0.328          -0.328
## diff_score_scaled -0.678 0.285    0.142           0.142   -0.143          -0.143   -0.328          -0.328
##                   diff_score diff_score_scaled
## means_y1               0.142             0.142
## means_y1_scaled        0.142             0.142
## means_y2              -0.143            -0.143
## means_y2_scaled       -0.143            -0.143
## gdi.z.cm              -0.328            -0.328
## gdi.z.cm_scaled       -0.328            -0.328
## diff_score             1.000             1.000
## diff_score_scaled      1.000             1.000
```

``` r
round(ddsc_mod2_GDI_cntry_year$results,3)
```

```
##                            estimate    SE      df t.ratio p.value ci.lower ci.upper
## r_xy1y2                       0.139 0.205 425.502   0.679   0.498   -0.264    0.542
## w_11                         -0.163 0.045  35.499  -3.607   0.001   -0.254   -0.071
## w_21                         -0.175 0.043  35.733  -4.015   0.000   -0.263   -0.086
## r_xy1                        -0.534 0.148  35.499  -3.607   0.001   -0.834   -0.234
## r_xy2                        -0.574 0.143  35.733  -4.015   0.000   -0.863   -0.284
## b_11                         -0.534 0.148  35.499  -3.607   0.001   -0.834   -0.234
## b_21                         -0.574 0.143  35.733  -4.015   0.000   -0.863   -0.284
## main_effect                  -0.169 0.043  35.175  -3.886   0.000   -0.257   -0.081
## moderator_effect             -0.217 0.015 358.437 -14.949   0.000   -0.246   -0.189
## interaction                  -0.012 0.018 425.502  -0.679   0.498   -0.047    0.023
## q_b11_b21                     0.057    NA      NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.057    NA      NA      NA      NA       NA       NA
## cross_over_point            -18.000    NA      NA      NA      NA       NA       NA
## interaction_vs_main          -0.157 0.048  36.825  -3.237   0.003   -0.254   -0.059
## interaction_vs_main_bscale   -0.514 0.159  36.825  -3.237   0.003   -0.836   -0.192
## interaction_vs_main_rscale   -0.514 0.159  36.825  -3.237   0.003   -0.836   -0.192
## dadas                        -0.325 0.090  35.499  -3.607   1.000   -0.508   -0.142
## dadas_bscale                 -1.068 0.296  35.499  -3.607   1.000   -1.668   -0.467
## dadas_rscale                 -1.068 0.296  35.499  -3.607   1.000   -1.669   -0.467
## abs_diff                      0.012 0.018 425.502   0.679   0.249   -0.023    0.047
## abs_sum                       0.337 0.087  35.175   3.886   0.000    0.161    0.513
## abs_diff_bscale               0.040 0.058 425.502   0.679   0.249   -0.075    0.155
## abs_sum_bscale                1.107 0.285  35.175   3.886   0.000    0.529    1.686
## abs_diff_rscale               0.040 0.058 425.017   0.677   0.249   -0.075    0.155
## abs_sum_rscale                1.107 0.285  35.175   3.886   0.000    0.529    1.686
```

``` r
round(ddsc_mod2_GDI_cntry_year$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.001 -1.000  0.016  1.000  0.898
```

### Comparisons across models


``` r
# full multilevel model (individuals as entries)
round(ddsc_mod2_GDI$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df  t.ratio p.value ci.lower ci.upper
## r_xy1             -0.4590 0.1459 34.0212  -3.1466  0.0034  -0.7554  -0.1626
## r_xy2             -0.5549 0.1391 34.0598  -3.9907  0.0003  -0.8375  -0.2724
## b_11              -0.4627 0.1471 34.0212  -3.1466  0.0034  -0.7616  -0.1639
## b_21              -0.5504 0.1379 34.0598  -3.9907  0.0003  -0.8307  -0.2701
## main_effect       -0.1611 0.0448 34.0193  -3.5970  0.0010  -0.2521  -0.0701
## moderator_effect  -0.2041 0.0137 33.8886 -14.9395  0.0000  -0.2318  -0.1763
## interaction       -0.0279 0.0141 36.0452  -1.9785  0.0556  -0.0564   0.0007
## q_b11_b21          0.1182     NA      NA       NA      NA       NA       NA
```

``` r
# country-level structural path model (country-gender means as entries)
round(ddsc_sem_GDI$results[c("r_xy1","r_xy2","b_11","b_21","q_b11_b21","diff_b11_b21"),],4)
```

```
##                  est     se       z pvalue ci.lower ci.upper
## r_xy1        -0.4607 0.1522 -3.0267 0.0025  -0.7590  -0.1624
## r_xy2        -0.5507 0.1431 -3.8473 0.0001  -0.8313  -0.2702
## b_11         -0.4569 0.1510 -3.0267 0.0025  -0.7528  -0.1610
## b_21         -0.5552 0.1443 -3.8473 0.0001  -0.8381  -0.2724
## q_b11_b21     0.1325 0.0588  2.2522 0.0243   0.0172   0.2477
## diff_b11_b21  0.0983 0.0420  2.3403 0.0193   0.0160   0.1806
```

``` r
# multilevel model (country-gender means as entries)
round(ddsc_mod2_GDI_ri$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df  t.ratio p.value ci.lower ci.upper
## r_xy1             -0.4607 0.1535 33.3201  -3.0018  0.0051  -0.7729  -0.1486
## r_xy2             -0.5507 0.1510 33.3201  -3.6474  0.0009  -0.8578  -0.2437
## b_11              -0.4569 0.1522 33.3201  -3.0018  0.0051  -0.7665  -0.1474
## b_21              -0.5552 0.1522 33.3201  -3.6474  0.0009  -0.8648  -0.2456
## main_effect       -0.1609 0.0479 32.0000  -3.3587  0.0020  -0.2585  -0.0633
## moderator_effect  -0.2065 0.0136 32.0000 -15.2274  0.0000  -0.2341  -0.1789
## interaction       -0.0313 0.0138 32.0000  -2.2704  0.0301  -0.0593  -0.0032
## q_b11_b21          0.1325     NA      NA       NA      NA       NA       NA
```

``` r
# multilevel model (country-year-gender means as entries)
round(ddsc_mod2_GDI_cntry_year$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE       df  t.ratio p.value ci.lower ci.upper
## r_xy1             -0.5339 0.1480  35.4994  -3.6069  0.0009  -0.8343  -0.2336
## r_xy2             -0.5735 0.1428  35.7334  -4.0151  0.0003  -0.8633  -0.2838
## b_11              -0.5339 0.1480  35.4994  -3.6069  0.0009  -0.8342  -0.2335
## b_21              -0.5736 0.1429  35.7334  -4.0151  0.0003  -0.8634  -0.2838
## main_effect       -0.1686 0.0434  35.1753  -3.8860  0.0004  -0.2566  -0.0805
## moderator_effect  -0.2174 0.0145 358.4374 -14.9492  0.0000  -0.2460  -0.1888
## interaction       -0.0121 0.0178 425.5020  -0.6787  0.4977  -0.0471   0.0229
## q_b11_b21          0.0573     NA       NA       NA      NA       NA       NA
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
## Time difference of 1.502713 hours
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
##                     Estimate         SE         2.5%         97.5%
## X.Intercept.    -0.003982609 0.04446277 -0.085828525  0.0837675052
## gndr.c          -0.203604728 0.01429943 -0.231802503 -0.1775111097
## gdi.z.cm        -0.159100800 0.04639836 -0.252470410 -0.0725498413
## gndr.c.gdi.z.cm -0.028426512 0.01458566 -0.056381647  0.0007925688
## w11             -0.144887544 0.04845836 -0.246033246 -0.0548539822
## w21             -0.173314056 0.04542878 -0.265458058 -0.0910401446
## b11             -0.455696929 0.15241010 -0.773818036 -0.1725254672
## b21             -0.545103332 0.14288153 -0.834912501 -0.2863373424
## r_xy1           -0.452004171 0.15117504 -0.767547372 -0.1711274005
## r_xy2           -0.549593368 0.14405845 -0.841789705 -0.2886959135
## q_b              0.128254525 0.07920793 -0.004330770  0.2905013275
## q                0.141703747 0.08469238  0.009666929  0.3139173011
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
## [1] 0.1282545
## 
## $se
## [1] 0.07920793
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
## [1] 2.881713
## 
## $p_low
## [1] 0.001977599
## 
## $z_high
## [1] 0.3567133
## 
## $p_high
## [1] 0.6393468
## 
## $ci_level
## [1] 0.9
## 
## $ci_lower
## [1] -0.002030928
## 
## $ci_upper
## [1] 0.25854
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
## [1] 0.1417037
## 
## $se
## [1] 0.08469238
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
## [1] 2.853902
## 
## $p_low
## [1] 0.002159293
## 
## $z_high
## [1] 0.4924144
## 
## $p_high
## [1] 0.6887868
## 
## $ci_level
## [1] 0.9
## 
## $ci_lower
## [1] 0.002397176
## 
## $ci_upper
## [1] 0.2810103
## 
## $equivalent
## [1] FALSE
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
  xlab("Gender Development Index (Average 2002-2023)")+
  #ylim(c(min.y.pred,max.y.pred))+
  #xlim(c(0.60,1.00))+
  ylab("Value benevolence (Average 2002-2023)")+
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
  xlab("Gender Development Index (Average 2002-2023)")+
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
## Warning in geom_text(inherit.aes = F, aes(x = 0.9, y = -0.45, label = coef_q, : All aesthetics have length 1, but the data has 262 rows.
## ℹ Please consider using `annotate()` or provide this layer with data containing a single row.
```

``` r
pflag_comb
```

![](Analysis_benevolence_files/figure-html/unnamed-chunk-35-1.png)<!-- -->

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
## 1462802.7 1462891.6 -731393.3 1462786.7    492335 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.4731 -0.5144  0.0692  0.6615  5.3268 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr 
##  cntry    (Intercept) 0.074251 0.27249       
##           gndr.c      0.004817 0.06941  0.29 
##  Residual             1.003863 1.00193       
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##                       Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)          0.0001206  0.0467741 33.9697251   0.003  0.99796    
## gndr.c              -0.2040101  0.0123449 32.5635204 -16.526  < 2e-16 ***
## log_gdp.z.cm         0.1304614  0.0469156 33.9963493   2.781  0.00878 ** 
## gndr.c:log_gdp.z.cm -0.0416302  0.0125061 33.7954485  -3.329  0.00212 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c lg_g..
## gndr.c       0.283              
## lg_gdp.z.cm  0.023  0.006       
## gndr.c:l_..  0.006 -0.001  0.280
```

``` r
getFE(mod2_log_GDP,round=3)
```

```
##                       Est.    SE     df       t     p     LL     UL
## (Intercept)          0.000 0.047 33.970   0.003 0.998 -0.095  0.095
## gndr.c              -0.204 0.012 32.564 -16.526 0.000 -0.229 -0.179
## log_gdp.z.cm         0.130 0.047 33.996   2.781 0.009  0.035  0.226
## gndr.c:log_gdp.z.cm -0.042 0.013 33.795  -3.329 0.002 -0.067 -0.016
```

``` r
getVC(mod2_log_GDP)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.27 0.07
## 2    cntry      gndr.c   <NA>  0.07 0.00
## 3    cntry (Intercept) gndr.c  0.29 0.01
## 4 Residual        <NA>   <NA>  1.00 1.00
```

``` r
r2mlm(mod2_log_GDP,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.022488422
## slope variation 0.001084597
## mean variation  0.066891290
## sigma2          0.909535690
## 
## $R2s
##           total
## f   0.022488422
## v   0.001084597
## m   0.066891290
## fv  0.023573019
## fvm 0.090464310
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
## Time difference of 38.52756 secs
```

``` r
round(ddsc_mod2_log_GDP$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.093        0.305        1.004     1.097 0.085   7802.647 0.998   0.999
## 2        0.5         0.093        0.305        1.004     1.097 0.085   6678.029 0.998   0.998
```

``` r
round(ddsc_mod2_log_GDP$descriptives,3)
```

```
##                          M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled log_gdp.z.cm
## means_y1            -0.117 0.321    1.000           1.000    0.965           0.965        0.367
## means_y1_scaled     -0.368 1.008    1.000           1.000    0.965           0.965        0.367
## means_y2             0.089 0.315    0.965           0.965    1.000           1.000        0.484
## means_y2_scaled      0.281 0.992    0.965           0.965    1.000           1.000        0.484
## log_gdp.z.cm        -0.024 1.012    0.367           0.367    0.484           0.484        1.000
## log_gdp.z.cm_scaled  0.000 1.000    0.367           0.367    0.484           0.484        1.000
## diff_score          -0.206 0.084    0.192           0.192   -0.071          -0.071       -0.419
## diff_score_scaled   -0.649 0.264    0.192           0.192   -0.071          -0.071       -0.419
##                     log_gdp.z.cm_scaled diff_score diff_score_scaled
## means_y1                          0.367      0.192             0.192
## means_y1_scaled                   0.367      0.192             0.192
## means_y2                          0.484     -0.071            -0.071
## means_y2_scaled                   0.484     -0.071            -0.071
## log_gdp.z.cm                      1.000     -0.419            -0.419
## log_gdp.z.cm_scaled               1.000     -0.419            -0.419
## diff_score                       -0.419      1.000             1.000
## diff_score_scaled                -0.419      1.000             1.000
```

``` r
round(ddsc_mod2_log_GDP$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                       0.496 0.149 33.795   3.329   0.002    0.193    0.799
## w_11                          0.151 0.046 33.985   3.320   0.002    0.059    0.244
## w_21                          0.110 0.049 33.998   2.236   0.032    0.010    0.209
## r_xy1                         0.472 0.142 33.985   3.320   0.002    0.183    0.761
## r_xy2                         0.348 0.155 33.998   2.236   0.032    0.032    0.664
## b_11                          0.476 0.143 33.985   3.320   0.002    0.185    0.767
## b_21                          0.345 0.154 33.998   2.236   0.032    0.031    0.658
## main_effect                   0.130 0.047 33.996   2.781   0.009    0.035    0.226
## moderator_effect             -0.204 0.012 32.564 -16.526   0.000   -0.229   -0.179
## interaction                  -0.042 0.013 33.795  -3.329   0.002   -0.067   -0.016
## q_b11_b21                     0.158    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.150    NA     NA      NA      NA       NA       NA
## cross_over_point             -4.901    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.089 0.052 34.005  -1.714   0.096   -0.194    0.016
## interaction_vs_main_bscale   -0.279 0.163 34.005  -1.714   0.096   -0.611    0.052
## interaction_vs_main_rscale   -0.286 0.165 34.005  -1.726   0.093   -0.622    0.051
## dadas                        -0.219 0.098 33.998  -2.236   0.984   -0.419   -0.020
## dadas_bscale                 -0.690 0.308 33.998  -2.236   0.984   -1.317   -0.063
## dadas_rscale                 -0.695 0.311 33.998  -2.236   0.984   -1.327   -0.063
## abs_diff                      0.042 0.013 33.795   3.329   0.001    0.016    0.067
## abs_sum                       0.261 0.094 33.996   2.781   0.004    0.070    0.452
## abs_diff_bscale               0.131 0.039 33.795   3.329   0.001    0.051    0.211
## abs_sum_bscale                0.821 0.295 33.996   2.781   0.004    0.221    1.420
## abs_diff_rscale               0.124 0.040 33.879   3.100   0.002    0.043    0.206
## abs_sum_rscale                0.820 0.295 33.996   2.776   0.004    0.220    1.420
```

``` r
round(ddsc_mod2_log_GDP$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
##      0      0      0      1      1
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
##                                    est    se      z pvalue ci.lower ci.upper
## r_xy1_y2                         0.419 0.156  2.687  0.007    0.113    0.724
## r_xy1                            0.484 0.150  3.228  0.001    0.190    0.778
## r_xy2                            0.367 0.160  2.300  0.021    0.054    0.680
## b_11                             0.480 0.149  3.228  0.001    0.189    0.772
## b_21                             0.370 0.161  2.300  0.021    0.055    0.685
## b_10                             0.281 0.147  1.919  0.055   -0.006    0.569
## b_20                            -0.368 0.158 -2.323  0.020   -0.679   -0.057
## res_cov_y1_y2                    0.764 0.188  4.055  0.000    0.395    1.134
## diff_b10_b20                     0.649 0.040 16.039  0.000    0.570    0.729
## diff_b11_b21                     0.110 0.041  2.687  0.007    0.030    0.191
## diff_rxy1_rxy2                   0.117 0.040  2.903  0.004    0.038    0.197
## q_b11_b21                        0.135 0.049  2.771  0.006    0.040    0.231
## q_rxy1_rxy2                      0.144 0.050  2.893  0.004    0.046    0.241
## cross_over_point                -5.881 2.219 -2.650  0.008  -10.230   -1.531
## sum_b11_b21                      0.850 0.307  2.769  0.006    0.248    1.452
## main_effect                      0.425 0.154  2.769  0.006    0.124    0.726
## interaction_vs_main_effect      -0.315 0.170 -1.848  0.065   -0.648    0.019
## diff_abs_b11_abs_b21             0.110 0.041  2.687  0.007    0.030    0.191
## abs_diff_b11_b21                 0.110 0.041  2.687  0.004    0.030    0.191
## abs_sum_b11_b21                  0.850 0.307  2.769  0.003    0.248    1.452
## dadas                           -0.740 0.322 -2.300  0.989   -1.370   -0.109
## q_r_equivalence                  0.044 0.050  0.880  0.811       NA       NA
## q_b_equivalence                  0.035 0.049  0.720  0.764       NA       NA
## cross_over_point_equivalence     5.881 2.219  2.650  0.996       NA       NA
## cross_over_point_minimal_effect  5.881 2.219  2.650  0.004       NA       NA
```

``` r
round(ddsc_sem_log_GDP$variance_test,3)
```

```
##              est    se      z pvalue ci.lower ci.upper
## cov_y1y2   0.937 0.231  4.050  0.000    0.483    1.390
## var_y1     0.955 0.232  4.123  0.000    0.501    1.409
## var_y2     0.986 0.239  4.123  0.000    0.518    1.455
## var_diff  -0.032 0.087 -0.364  0.716   -0.203    0.139
## var_ratio  0.968 0.087 11.167  0.000    0.798    1.138
## cor_y1y2   0.965 0.012 82.576  0.000    0.942    0.988
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
## r_xy1y2                       0.419 0.161 32.000   2.607   0.014    0.091    0.746
## w_11                          0.153 0.051 33.146   3.008   0.005    0.049    0.256
## w_21                          0.118 0.051 33.146   2.316   0.027    0.014    0.221
## r_xy1                         0.484 0.161 33.146   3.008   0.005    0.157    0.812
## r_xy2                         0.367 0.158 33.146   2.316   0.027    0.045    0.689
## b_11                          0.480 0.160 33.146   3.008   0.005    0.156    0.805
## b_21                          0.370 0.160 33.146   2.316   0.027    0.045    0.695
## main_effect                   0.135 0.050 32.000   2.686   0.011    0.033    0.238
## moderator_effect             -0.206 0.013 32.000 -15.560   0.000   -0.234   -0.179
## interaction                  -0.035 0.013 32.000  -2.607   0.014   -0.063   -0.008
## q_b11_b21                     0.135    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.144    NA     NA      NA      NA       NA       NA
## cross_over_point             -5.881    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.100 0.052 36.561  -1.921   0.063   -0.206    0.006
## interaction_vs_main_bscale   -0.315 0.164 36.561  -1.921   0.063   -0.647    0.017
## interaction_vs_main_rscale   -0.308 0.161 36.674  -1.911   0.064   -0.635    0.019
## dadas                        -0.235 0.102 33.146  -2.316   0.987   -0.442   -0.029
## dadas_bscale                 -0.740 0.319 33.146  -2.316   0.987   -1.390   -0.090
## dadas_rscale                 -0.734 0.317 33.146  -2.316   0.987   -1.378   -0.089
## abs_diff                      0.035 0.013 32.000   2.607   0.007    0.008    0.063
## abs_sum                       0.270 0.101 32.000   2.686   0.006    0.065    0.475
## abs_diff_bscale               0.110 0.042 32.000   2.607   0.007    0.024    0.197
## abs_sum_bscale                0.850 0.317 32.000   2.686   0.006    0.205    1.495
## abs_diff_rscale               0.117 0.042 32.239   2.765   0.005    0.031    0.204
## abs_sum_rscale                0.851 0.317 32.000   2.689   0.006    0.206    1.496
```

``` r
# country-time multilevel model


mod2_log_GDP_cntry_year<-
  lmer(ben.z.wt~gndr.c+log_gdp.z.cm+gndr.c:log_gdp.z.cm+
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
## Formula: ben.z.wt ~ gndr.c + log_gdp.z.cm + gndr.c:log_gdp.z.cm + (gndr.c |      cntry)
##    Data: diff_dat_cntry_year
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##    -266.1    -231.7     141.1    -282.1       538 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -5.6974 -0.5667  0.1323  0.6531  3.3632 
## 
## Random effects:
##  Groups   Name        Variance  Std.Dev. Corr 
##  cntry    (Intercept) 0.0680884 0.26094       
##           gndr.c      0.0006171 0.02484  1.00 
##  Residual             0.0278957 0.16702       
## Number of obs: 546, groups:  cntry, 34
## 
## Fixed effects:
##                      Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)           0.00406    0.04546  33.34284   0.089  0.92937    
## gndr.c               -0.20730    0.01527 143.95319 -13.579  < 2e-16 ***
## log_gdp.z.cm          0.13124    0.04581  33.96437   2.865  0.00711 ** 
## gndr.c:log_gdp.z.cm  -0.04732    0.01629 170.48227  -2.905  0.00416 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c lg_g..
## gndr.c       0.276              
## lg_gdp.z.cm  0.011  0.005       
## gndr.c:l_..  0.005 -0.203  0.260
## optimizer (bobyqa) convergence code: 0 (OK)
## boundary (singular) fit: see help('isSingular')
```

``` r
getFE(mod2_log_GDP_cntry_year,round=3)
```

```
##                       Est.    SE      df       t     p     LL     UL
## (Intercept)          0.004 0.045  33.343   0.089 0.929 -0.088  0.097
## gndr.c              -0.207 0.015 143.953 -13.579 0.000 -0.237 -0.177
## log_gdp.z.cm         0.131 0.046  33.964   2.865 0.007  0.038  0.224
## gndr.c:log_gdp.z.cm -0.047 0.016 170.482  -2.905 0.004 -0.079 -0.015
```

``` r
getVC(mod2_log_GDP_cntry_year)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.26 0.07
## 2    cntry      gndr.c   <NA>  0.02 0.00
## 3    cntry (Intercept) gndr.c  1.00 0.01
## 4 Residual        <NA>   <NA>  0.17 0.03
```

``` r
r2mlm(mod2_log_GDP,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.022488422
## slope variation 0.001084597
## mean variation  0.066891290
## sigma2          0.909535690
## 
## $R2s
##           total
## f   0.022488422
## v   0.001084597
## m   0.066891290
## fv  0.023573019
## fvm 0.090464310
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
## 1       -0.5         0.086        0.294        0.028     0.115 0.752      8.029 0.995    0.96
## 2        0.5         0.085        0.292        0.028     0.113 0.749      8.029 0.995    0.96
```

``` r
round(ddsc_mod2_log_GDP_cntry_year$descriptives,3)
```

```
##                          M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled log_gdp.z.cm
## means_y1            -0.101 0.304    1.000           1.000    0.959           0.959        0.355
## means_y1_scaled     -0.331 1.000    1.000           1.000    0.959           0.959        0.355
## means_y2             0.106 0.304    0.959           0.959    1.000           1.000        0.500
## means_y2_scaled      0.348 1.000    0.959           0.959    1.000           1.000        0.500
## log_gdp.z.cm        -0.024 1.012    0.355           0.355    0.500           0.500        1.000
## log_gdp.z.cm_scaled  0.000 1.000    0.355           0.355    0.500           0.500        1.000
## diff_score          -0.207 0.087    0.142           0.142   -0.143          -0.143       -0.511
## diff_score_scaled   -0.678 0.285    0.142           0.142   -0.143          -0.143       -0.511
##                     log_gdp.z.cm_scaled diff_score diff_score_scaled
## means_y1                          0.355      0.142             0.142
## means_y1_scaled                   0.355      0.142             0.142
## means_y2                          0.500     -0.143            -0.143
## means_y2_scaled                   0.500     -0.143            -0.143
## log_gdp.z.cm                      1.000     -0.511            -0.511
## log_gdp.z.cm_scaled               1.000     -0.511            -0.511
## diff_score                       -0.511      1.000             1.000
## diff_score_scaled                -0.511      1.000             1.000
```

``` r
round(ddsc_mod2_log_GDP_cntry_year$results,3)
```

```
##                            estimate    SE      df t.ratio p.value ci.lower ci.upper
## r_xy1y2                       0.545 0.188 170.482   2.905   0.004    0.175    0.916
## w_11                          0.155 0.044  33.909   3.489   0.001    0.065    0.245
## w_21                          0.108 0.049  33.965   2.215   0.034    0.009    0.206
## r_xy1                         0.509 0.146  33.909   3.489   0.001    0.212    0.805
## r_xy2                         0.353 0.160  33.965   2.215   0.034    0.029    0.678
## b_11                          0.509 0.146  33.909   3.489   0.001    0.212    0.805
## b_21                          0.353 0.160  33.965   2.215   0.034    0.029    0.678
## main_effect                   0.131 0.046  33.964   2.865   0.007    0.038    0.224
## moderator_effect             -0.207 0.015 143.953 -13.579   0.000   -0.237   -0.177
## interaction                  -0.047 0.016 170.482  -2.905   0.004   -0.079   -0.015
## q_b11_b21                     0.192    NA      NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.192    NA      NA      NA      NA       NA       NA
## cross_over_point             -4.381    NA      NA      NA      NA       NA       NA
## interaction_vs_main          -0.084 0.052  34.091  -1.600   0.119   -0.191    0.023
## interaction_vs_main_bscale   -0.276 0.172  34.091  -1.600   0.119   -0.626    0.074
## interaction_vs_main_rscale   -0.276 0.172  34.091  -1.600   0.119   -0.626    0.074
## dadas                        -0.215 0.097  33.965  -2.215   0.983   -0.413   -0.018
## dadas_bscale                 -0.707 0.319  33.965  -2.215   0.983   -1.355   -0.058
## dadas_rscale                 -0.707 0.319  33.965  -2.215   0.983   -1.355   -0.058
## abs_diff                      0.047 0.016 170.482   2.905   0.002    0.015    0.079
## abs_sum                       0.262 0.092  33.964   2.865   0.004    0.076    0.449
## abs_diff_bscale               0.155 0.054 170.482   2.905   0.002    0.050    0.261
## abs_sum_bscale                0.862 0.301  33.964   2.865   0.004    0.251    1.474
## abs_diff_rscale               0.155 0.053 170.765   2.906   0.002    0.050    0.261
## abs_sum_rscale                0.862 0.301  33.964   2.865   0.004    0.251    1.474
```

``` r
round(ddsc_mod2_log_GDP_cntry_year$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.001 -1.000  0.016  1.000  0.898
```

### Comparisons across models


``` r
# full multilevel model (individuals as entries)
round(ddsc_mod2_log_GDP$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df  t.ratio p.value ci.lower ci.upper
## r_xy1              0.4719 0.1421 33.9851   3.3202  0.0022   0.1831   0.7608
## r_xy2              0.3477 0.1555 33.9982   2.2361  0.0320   0.0317   0.6637
## b_11               0.4758 0.1433 33.9851   3.3202  0.0022   0.1846   0.7670
## b_21               0.3449 0.1542 33.9982   2.2361  0.0320   0.0314   0.6583
## main_effect        0.1305 0.0469 33.9963   2.7808  0.0088   0.0351   0.2258
## moderator_effect  -0.2040 0.0123 32.5635 -16.5259  0.0000  -0.2291  -0.1789
## interaction       -0.0416 0.0125 33.7954  -3.3288  0.0021  -0.0671  -0.0162
## q_b11_b21          0.1579     NA      NA       NA      NA       NA       NA
```

``` r
# country-level structural path model (country-gender means as entries)
round(ddsc_sem_log_GDP$results[c("r_xy1","r_xy2","b_11","b_21","q_b11_b21","diff_b11_b21"),],4)
```

```
##                 est     se      z pvalue ci.lower ci.upper
## r_xy1        0.4843 0.1500 3.2281 0.0012   0.1903   0.7784
## r_xy2        0.3670 0.1595 2.3002 0.0214   0.0543   0.6796
## b_11         0.4804 0.1488 3.2281 0.0012   0.1887   0.7720
## b_21         0.3699 0.1608 2.3002 0.0214   0.0547   0.6852
## q_b11_b21    0.1351 0.0488 2.7713 0.0056   0.0396   0.2307
## diff_b11_b21 0.1104 0.0411 2.6869 0.0072   0.0299   0.1910
```

``` r
# multilevel model (country-gender means as entries)
round(ddsc_mod2_log_GDP_ri$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df  t.ratio p.value ci.lower ci.upper
## r_xy1              0.4843 0.1610 33.1456   3.0079  0.0050   0.1568   0.8119
## r_xy2              0.3670 0.1584 33.1456   2.3164  0.0269   0.0447   0.6892
## b_11               0.4804 0.1597 33.1456   3.0079  0.0050   0.1555   0.8053
## b_21               0.3700 0.1597 33.1456   2.3164  0.0269   0.0451   0.6948
## main_effect        0.1352 0.0503 32.0000   2.6859  0.0114   0.0327   0.2377
## moderator_effect  -0.2065 0.0133 32.0000 -15.5598  0.0000  -0.2335  -0.1795
## interaction       -0.0351 0.0135 32.0000  -2.6066  0.0138  -0.0625  -0.0077
## q_b11_b21          0.1351     NA      NA       NA      NA       NA       NA
```

``` r
# multilevel model (country-year-gender means as entries)
round(ddsc_mod2_log_GDP_cntry_year$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE       df  t.ratio p.value ci.lower ci.upper
## r_xy1              0.5088 0.1458  33.9087   3.4889  0.0014   0.2124   0.8052
## r_xy2              0.3533 0.1595  33.9649   2.2150  0.0336   0.0291   0.6775
## b_11               0.5088 0.1458  33.9087   3.4889  0.0014   0.2124   0.8051
## b_21               0.3534 0.1595  33.9649   2.2150  0.0336   0.0291   0.6776
## main_effect        0.1312 0.0458  33.9644   2.8648  0.0071   0.0381   0.2243
## moderator_effect  -0.2073 0.0153 143.9532 -13.5790  0.0000  -0.2375  -0.1771
## interaction       -0.0473 0.0163 170.4823  -2.9047  0.0042  -0.0795  -0.0152
## q_b11_b21          0.1918     NA       NA       NA      NA       NA       NA
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
## Time difference of 1.65398 hours
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
##                         Estimate         SE        2.5%       97.5%
## X.Intercept.        -0.001045279 0.04719321 -0.08817525  0.09181549
## gndr.c              -0.203762990 0.01267075 -0.22895525 -0.17880669
## log_gdp.z.cm         0.130078599 0.04576640  0.04167688  0.21875251
## gndr.c.log_gdp.z.cm -0.041034129 0.01212700 -0.06374782 -0.01792528
## w11                  0.150595663 0.04459521  0.06027659  0.23983428
## w21                  0.109561534 0.04768570  0.01687030  0.19998479
## b11                  0.473649972 0.14025982  0.18958054  0.75432119
## b21                  0.344590386 0.14997997  0.05306009  0.62898752
## r_xy1                0.469811731 0.13912322  0.18804427  0.74820852
## r_xy2                0.347428789 0.15121536  0.05349714  0.63416851
## q_b                  0.162610713 0.05315974  0.07100254  0.27428203
## q                    0.153287947 0.05060671  0.06240478  0.25463938
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
## [1] 0.1626107
## 
## $se
## [1] 0.05315974
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
## [1] 4.94003
## 
## $p_low
## [1] 3.905532e-07
## 
## $z_high
## [1] 1.177784
## 
## $p_high
## [1] 0.8805587
## 
## $ci_level
## [1] 0.9
## 
## $ci_lower
## [1] 0.07517072
## 
## $ci_upper
## [1] 0.2500507
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
## [1] 0.1532879
## 
## $se
## [1] 0.05060671
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
## [1] 5.005027
## 
## $p_low
## [1] 2.792712e-07
## 
## $z_high
## [1] 1.052982
## 
## $p_high
## [1] 0.8538253
## 
## $ci_level
## [1] 0.9
## 
## $ci_lower
## [1] 0.07004732
## 
## $ci_upper
## [1] 0.2365286
## 
## $equivalent
## [1] FALSE
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
  xlab("log(GDP per capita) (Average 2002-2023)")+
  #ylim(c(min.y.pred,max.y.pred))+
  #xlim(c(0.60,1.00))+
  ylab("Value benevolence (Average 2002-2023)")+
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
  xlab("log(GDP per capita) (Average 2002-2023)")+
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

![](Analysis_benevolence_files/figure-html/unnamed-chunk-41-1.png)<!-- -->

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


## mod3: fixed effect of time (year)

* Year coded year-2002 (2002 is a zero-point)


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
## 1461286.2 1461364.0 -730636.1 1461272.2    492336 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.5000 -0.5170  0.0671  0.6601  5.4507 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr 
##  cntry    (Intercept) 0.094988 0.30820       
##           gndr.c      0.006664 0.08163  0.00 
##  Residual             1.000735 1.00037       
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##               Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept) -8.149e-03  5.288e-02  3.399e+01  -0.154    0.878    
## gndr.c      -2.045e-01  1.438e-02  3.306e+01 -14.221 1.21e-15 ***
## essround.c   1.900e-02  4.847e-04  4.922e+05  39.192  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##            (Intr) gndr.c
## gndr.c      0.004       
## essround.c -0.003  0.000
```

``` r
getFE(mod3,round=3)
```

```
##               Est.    SE         df       t     p     LL     UL
## (Intercept) -0.008 0.053     33.989  -0.154 0.878 -0.116  0.099
## gndr.c      -0.205 0.014     33.060 -14.221 0.000 -0.234 -0.175
## essround.c   0.019 0.000 492249.802  39.192 0.000  0.018  0.020
```

``` r
getVC(mod3)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.31 0.09
## 2    cntry      gndr.c   <NA>  0.08 0.01
## 3    cntry (Intercept) gndr.c  0.00 0.00
## 4 Residual        <NA>   <NA>  1.00 1.00
```

``` r
r2mlm(mod3,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.012409039
## slope variation 0.001490239
## mean variation  0.085487046
## sigma2          0.900613676
## 
## $R2s
##           total
## f   0.012409039
## v   0.001490239
## m   0.085487046
## fv  0.013899278
## fvm 0.099386324
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
## mod2    6 1462818 1462884 -731403   1462806                         
## mod3    7 1461286 1461364 -730636   1461272 1533.6  1  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## mod4: random effect of time (year)


``` r
mod4<-lmer(ben.z~gndr.c+year.c+(gndr.c+year.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod4)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ben.z ~ gndr.c + year.c + (gndr.c + year.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1455970.2 1456081.3 -727975.1 1455950.2    492333 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.6627 -0.5174  0.0802  0.6604  5.3354 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr        
##  cntry    (Intercept) 0.826491 0.90912              
##           gndr.c      0.006817 0.08256   0.28       
##           year.c      0.002923 0.05407  -0.93 -0.36 
##  Residual             0.989495 0.99473              
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##              Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)  0.092582   0.156048 32.706151   0.593    0.557    
## gndr.c      -0.204779   0.014531 33.092988 -14.093 1.53e-15 ***
## year.c      -0.002696   0.009284 32.990741  -0.290    0.773    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##        (Intr) gndr.c
## gndr.c  0.273       
## year.c -0.930 -0.352
```

``` r
getFE(mod4,round=3)
```

```
##               Est.    SE     df       t     p     LL     UL
## (Intercept)  0.093 0.156 32.706   0.593 0.557 -0.225  0.410
## gndr.c      -0.205 0.015 33.093 -14.093 0.000 -0.234 -0.175
## year.c      -0.003 0.009 32.991  -0.290 0.773 -0.022  0.016
```

``` r
getVC(mod4)
```

```
##        grp        var1   var2 sdcor  vcov
## 1    cntry (Intercept)   <NA>  0.91  0.83
## 2    cntry      gndr.c   <NA>  0.08  0.01
## 3    cntry      year.c   <NA>  0.05  0.00
## 4    cntry (Intercept) gndr.c  0.28  0.02
## 5    cntry (Intercept) year.c -0.93 -0.05
## 6    cntry      gndr.c year.c -0.36  0.00
## 7 Residual        <NA>   <NA>  0.99  0.99
```

``` r
r2mlm(mod4,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.008143339
## slope variation 0.090599073
## mean variation  0.150080047
## sigma2          0.751177541
## 
## $R2s
##           total
## f   0.008143339
## v   0.090599073
## m   0.150080047
## fv  0.098742412
## fvm 0.248822459
```

``` r
anova(mod2,mod3,mod4)
```

```
## Data: diff_dat
## Models:
## mod2: ben.z ~ gndr.c + (gndr.c | cntry)
## mod3: ben.z ~ gndr.c + essround.c + (gndr.c | cntry)
## mod4: ben.z ~ gndr.c + year.c + (gndr.c + year.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod2    6 1462818 1462884 -731403   1462806                         
## mod3    7 1461286 1461364 -730636   1461272 1533.6  1  < 2.2e-16 ***
## mod4   10 1455970 1456081 -727975   1455950 5322.0  3  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## mod5: fixed interaction between time and gender


``` r
mod5<-lmer(ben.z~gndr.c+year.c+
             gndr.c:year.c+(gndr.c+year.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod5)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ben.z ~ gndr.c + year.c + gndr.c:year.c + (gndr.c + year.c |      cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1455970.6 1456092.8 -727974.3 1455948.6    492332 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.6653 -0.5177  0.0800  0.6605  5.3314 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr        
##  cntry    (Intercept) 0.826457 0.90910              
##           gndr.c      0.006802 0.08248   0.29       
##           year.c      0.002923 0.05407  -0.93 -0.37 
##  Residual             0.989492 0.99473              
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##                 Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)    9.272e-02  1.560e-01  3.270e+01   0.594    0.556    
## gndr.c        -1.985e-01  1.534e-02  4.117e+01 -12.936 4.31e-16 ***
## year.c        -2.708e-03  9.284e-03  3.299e+01  -0.292    0.772    
## gndr.c:year.c -5.939e-04  4.678e-04  2.230e+05  -1.270    0.204    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c year.c
## gndr.c       0.267              
## year.c      -0.930 -0.342       
## gndr.c:yr.c -0.001 -0.324  0.001
```

``` r
getFE(mod5,round=3)
```

```
##                 Est.    SE         df       t     p     LL     UL
## (Intercept)    0.093 0.156     32.701   0.594 0.556 -0.225  0.410
## gndr.c        -0.198 0.015     41.169 -12.936 0.000 -0.229 -0.167
## year.c        -0.003 0.009     32.986  -0.292 0.772 -0.022  0.016
## gndr.c:year.c -0.001 0.000 222995.308  -1.270 0.204 -0.002  0.000
```

``` r
getVC(mod5)
```

```
##        grp        var1   var2 sdcor  vcov
## 1    cntry (Intercept)   <NA>  0.91  0.83
## 2    cntry      gndr.c   <NA>  0.08  0.01
## 3    cntry      year.c   <NA>  0.05  0.00
## 4    cntry (Intercept) gndr.c  0.29  0.02
## 5    cntry (Intercept) year.c -0.93 -0.05
## 6    cntry      gndr.c year.c -0.37  0.00
## 7 Residual        <NA>   <NA>  0.99  0.99
```

``` r
r2mlm(mod5,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.008124025
## slope variation 0.090598564
## mean variation  0.150077288
## sigma2          0.751200123
## 
## $R2s
##           total
## f   0.008124025
## v   0.090598564
## m   0.150077288
## fv  0.098722589
## fvm 0.248799877
```

``` r
anova(mod4,mod5)
```

```
## Data: diff_dat
## Models:
## mod4: ben.z ~ gndr.c + year.c + (gndr.c + year.c | cntry)
## mod5: ben.z ~ gndr.c + year.c + gndr.c:year.c + (gndr.c + year.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)
## mod4   10 1455970 1456081 -727975   1455950                     
## mod5   11 1455971 1456093 -727974   1455949 1.6077  1     0.2048
```

## mod6: random interaction between time and gender


``` r
mod6<-lmer(ben.z~gndr.c+year.c+
             gndr.c:year.c+(gndr.c+year.c+gndr.c:year.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ben.z ~ gndr.c + year.c + gndr.c:year.c + (gndr.c + year.c +      gndr.c:year.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1455835.4 1456002.0 -727902.7 1455805.4    492328 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.6545 -0.5181  0.0762  0.6584  5.3253 
## 
## Random effects:
##  Groups   Name          Variance  Std.Dev. Corr              
##  cntry    (Intercept)   0.8272003 0.909506                   
##           gndr.c        0.0161322 0.127013 -0.12             
##           year.c        0.0029255 0.054088 -0.93  0.03       
##           gndr.c:year.c 0.0000501 0.007078  0.31 -0.83 -0.23 
##  Residual               0.9890973 0.994534                   
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##                 Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)    0.0931569  0.1561143 32.6483325   0.597    0.555    
## gndr.c        -0.2096145  0.0232804 27.0092120  -9.004 1.28e-09 ***
## year.c        -0.0027638  0.0092870 32.9394849  -0.298    0.768    
## gndr.c:year.c -0.0002652  0.0013519 23.0583548  -0.196    0.846    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c year.c
## gndr.c      -0.115              
## year.c      -0.930  0.024       
## gndr.c:yr.c  0.279 -0.837 -0.207
```

``` r
getFE(mod6,round=3)
```

```
##                 Est.    SE     df      t     p     LL     UL
## (Intercept)    0.093 0.156 32.648  0.597 0.555 -0.225  0.411
## gndr.c        -0.210 0.023 27.009 -9.004 0.000 -0.257 -0.162
## year.c        -0.003 0.009 32.939 -0.298 0.768 -0.022  0.016
## gndr.c:year.c  0.000 0.001 23.058 -0.196 0.846 -0.003  0.003
```

``` r
getVC(mod6)
```

```
##         grp          var1          var2 sdcor  vcov
## 1     cntry   (Intercept)          <NA>  0.91  0.83
## 2     cntry        gndr.c          <NA>  0.13  0.02
## 3     cntry        year.c          <NA>  0.05  0.00
## 4     cntry gndr.c:year.c          <NA>  0.01  0.00
## 5     cntry   (Intercept)        gndr.c -0.12 -0.01
## 6     cntry   (Intercept)        year.c -0.93 -0.05
## 7     cntry   (Intercept) gndr.c:year.c  0.31  0.00
## 8     cntry        gndr.c        year.c  0.03  0.00
## 9     cntry        gndr.c gndr.c:year.c -0.83  0.00
## 10    cntry        year.c gndr.c:year.c -0.23  0.00
## 11 Residual          <NA>          <NA>  0.99  0.99
```

``` r
r2mlm(mod6,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.008733771
## slope variation 0.090983813
## mean variation  0.150559621
## sigma2          0.749722794
## 
## $R2s
##           total
## f   0.008733771
## v   0.090983813
## m   0.150559621
## fv  0.099717585
## fvm 0.250277206
```

``` r
anova(mod4,mod5,mod6)
```

```
## Data: diff_dat
## Models:
## mod4: ben.z ~ gndr.c + year.c + (gndr.c + year.c | cntry)
## mod5: ben.z ~ gndr.c + year.c + gndr.c:year.c + (gndr.c + year.c | cntry)
## mod6: ben.z ~ gndr.c + year.c + gndr.c:year.c + (gndr.c + year.c + gndr.c:year.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L)    Chisq Df Pr(>Chisq)    
## mod4   10 1455970 1456081 -727975   1455950                           
## mod5   11 1455971 1456093 -727974   1455949   1.6077  1     0.2048    
## mod6   15 1455835 1456002 -727903   1455805 143.2740  4     <2e-16 ***
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
##      21  0.1427 0.0778 33.9  -0.0154   0.3008   1.835  0.0753
##       0  0.1980 0.1580 31.5  -0.1238   0.5197   1.254  0.2191
## 
## gndr.c =  0.5:
##  year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21 -0.0725 0.0750 33.9  -0.2248   0.0799  -0.967  0.3404
##       0 -0.0117 0.1550 31.6  -0.3280   0.3047  -0.075  0.9406
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
##  year.c21 - year.c0  -0.0553 0.198 32.0   -0.459    0.349  -0.278  0.7825
## 
## gndr.c =  0.5:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0  -0.0608 0.193 32.1   -0.453    0.331  -0.316  0.7542
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
##    -0.5     21  0.1427 0.0778 33.9  -0.0154   0.3008   1.835  0.0753
##     0.5     21 -0.0725 0.0750 33.9  -0.2248   0.0799  -0.967  0.3404
##    -0.5      0  0.1980 0.1580 31.5  -0.1238   0.5197   1.254  0.2191
##     0.5      0 -0.0117 0.1550 31.6  -0.3280   0.3047  -0.075  0.9406
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod6,adjust="none",infer=c(T,T))
```

```
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21     0.2152 0.0155 29.9    0.183    0.247  13.850 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)  -0.0553 0.1980 32.0   -0.459    0.349  -0.278  0.7825
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0      0.1544 0.1980 32.5   -0.248    0.557   0.780  0.4410
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)     -0.2704 0.1930 32.3   -0.664    0.123  -1.400  0.1710
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0        -0.0608 0.1930 32.1   -0.453    0.331  -0.316  0.7542
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0       0.2096 0.0233 27.0    0.162    0.257   9.004 <0.0001
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
##  diff_ESS11   -0.215 0.0155 29.9   -0.247   -0.183 -13.850 <0.0001
##  diff_ESS1    -0.210 0.0233 27.0   -0.257   -0.162  -9.004 <0.0001
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
##  diff_ESS11 - diff_ESS1 -0.00557 0.0284 23.1  -0.0643   0.0531  -0.196  0.8462
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
## [1] -0.005569526
## 
## $se
## [1] 0.02838894
## 
## $df
## [1] 23.05835
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
## [1] 6.84881
## 
## $p_low
## [1] 2.727776e-07
## 
## $t_high
## [1] -7.241183
## 
## $p_high
## [1] 1.117276e-07
## 
## $ci_level
## [1] 0.8
## 
## $ci_lower
## [1] -0.0430248
## 
## $ci_upper
## [1] 0.03188575
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
  scale_x_continuous(limits = c(2001, 2024),
                     breaks = c(seq(2002, 2020, 2),2023))+
  ylab("Mean-level of value benevolence")+
  theme(legend.title=element_blank(),
        legend.text  = element_text(size = 16),
        axis.title.x = element_text(size = 18, face = "bold"),   
        axis.title.y = element_text(size = 18, face = "bold"),   
        axis.text.x  = element_text(size = 14),                  
        axis.text.y  = element_text(size = 14))

p_time_trends
```

![](Analysis_benevolence_files/figure-html/unnamed-chunk-48-1.png)<!-- -->

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
pred_cntry_dat$ben.z_mean<-predict(mod6,newdata=pred_cntry_dat)

pred_cntry_dat$year=pred_cntry_dat$year.c+2002

pred_cntry_dat$gender<-
  case_when(
    pred_cntry_dat$gndr.c==0.5~"men",
    pred_cntry_dat$gndr.c==-0.5~"women",
  )

range(pred_cntry_dat$ben.z_mean)
```

```
## [1] -1.2428290  0.6663235
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

pdf("../results/ben/country_specific_time_trends.pdf", width = 8, height = 6)

for (ctry in countries) {
  print(
    ggplot(pred_cntry_dat[pred_cntry_dat$cntry == ctry, ],
       aes(x = year, y = ben.z_mean, color = gender)) +
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
    name   = "Mean-level of value benevolence",
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
## Warning: Removed 2 rows containing missing values or values outside the scale range (`geom_point()`).
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
         aes(x = year, y = ben.z_mean, color = gender)) +
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
    name   = "Mean-level of value benevolence",
    sec.axis = sec_axis(~ scale_y_to_gei(.),
                        name = "Gender Equality Index (GEI)")
  ) +
  scale_x_continuous(limits = c(2001, 2024),
                     breaks = c(seq(2002, 2020, 2),2023)) +
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
## Warning: Removed 2 rows containing non-finite outside the scale range (`stat_smooth()`).
```

```
## Warning: Removed 2 rows containing missing values or values outside the scale range (`geom_point()`).
```

```
## Warning: Removed 4 rows containing missing values or values outside the scale range (`geom_point()`).
```

```
## Warning: Removed 12 rows containing missing values or values outside the scale range (`geom_line()`).
```

```
## Warning: Removed 2 rows containing missing values or values outside the scale range (`geom_flag()`).
```

![](Analysis_benevolence_files/figure-html/unnamed-chunk-50-1.png)<!-- -->

``` r
png(filename = 
      "../results/ben/country_time_trend_facets.png",
    units = "cm",
    width = 21.0,height=29.7*(3/4),res = 600)
facet_plot
```

```
## Warning: Removed 2 rows containing non-finite outside the scale range (`stat_smooth()`).
```

```
## Warning: Removed 2 rows containing missing values or values outside the scale range (`geom_point()`).
```

```
## Warning: Removed 4 rows containing missing values or values outside the scale range (`geom_point()`).
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
## [1] 58.39378
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
## 1   -0.23                0.34                    0.00                    0.34                      0.34
## 2   -0.24                0.20                    0.01                    0.20                      0.19
## 3   -0.11               -0.24                    0.04                   -0.22                     -0.26
## 4   -0.19                0.29                   -0.05                    0.27                      0.32
## 5   -0.12                0.02                   -0.10                   -0.03                      0.07
## 6   -0.34               -0.20                    0.11                   -0.15                     -0.26
## 7   -0.23                0.52                    0.05                    0.55                      0.50
## 8   -0.27                0.37                   -0.03                    0.36                      0.39
## 9   -0.28                0.23                   -0.01                    0.22                      0.23
## 10  -0.12                0.29                   -0.06                    0.26                      0.32
## 11  -0.44                0.57                    0.17                    0.65                      0.49
## 12  -0.25                0.26                   -0.04                    0.24                      0.28
## 13  -0.24                0.36                   -0.02                    0.35                      0.37
## 14   0.00               -0.46                   -0.22                   -0.57                     -0.35
## 15  -0.16                0.13                   -0.03                    0.12                      0.14
## 16  -0.27               -0.25                    0.18                   -0.16                     -0.34
## 17  -0.28                0.15                    0.05                    0.17                      0.12
## 18  -0.10               -0.18                   -0.01                   -0.18                     -0.18
## 19  -0.43                0.39                    0.20                    0.48                      0.29
## 20  -0.24               -0.27                    0.12                   -0.21                     -0.33
## 21  -0.20               -0.04                   -0.04                   -0.06                     -0.03
## 22  -0.28               -0.21                   -0.07                   -0.24                     -0.17
## 23  -0.34               -6.27                    0.25                   -6.15                     -6.40
## 24  -0.27                0.42                    0.02                    0.43                      0.41
## 25  -0.28                0.56                    0.01                    0.57                      0.56
## 26  -0.13                0.06                   -0.08                    0.02                      0.10
## 27   0.02                0.34                   -0.23                    0.23                      0.46
## 28  -0.24                0.57                    0.03                    0.59                      0.55
## 29  -0.10               -0.47                   -0.02                   -0.48                     -0.46
## 30  -0.41                0.71                    0.17                    0.80                      0.63
## 31  -0.21                0.54                    0.00                    0.53                      0.54
## 32  -0.21               -0.13                    0.05                   -0.10                     -0.15
## 33   0.03               -0.70                   -0.17                   -0.78                     -0.61
## 34   0.03                0.13                   -0.47                   -0.11                      0.36
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
## 1     RU               -0.47
## 2     GR               -0.46
## 3     IT               -0.27
## 4     HU               -0.25
## 5     BG               -0.24
## 6     CZ               -0.20
## 7     IL               -0.18
## 8     SK               -0.13
## 9     LT               -0.04
## 10    CY                0.02
## 11    PL                0.06
## 12    HR                0.13
## 13    UA                0.13
## 14    IE                0.15
## 15    BE                0.20
## 16    EE                0.23
## 17    FR                0.26
## 18    CH                0.29
## 19    ES                0.29
## 20    AT                0.34
## 21    PT                0.34
## 22    GB                0.36
## 23    DK                0.37
## 24    IS                0.39
## 25    NL                0.42
## 26    DE                0.52
## 27    SI                0.54
## 28    NO                0.56
## 29    FI                0.57
## 30    SE                0.71
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
## 1     UA                   -0.47
## 2     PT                   -0.23
## 3     GR                   -0.22
## 4     CY                   -0.10
## 5     PL                   -0.08
## 6     ES                   -0.06
## 7     CH                   -0.05
## 8     FR                   -0.04
## 9     LT                   -0.04
## 10    DK                   -0.03
## 11    HR                   -0.03
## 12    GB                   -0.02
## 13    RU                   -0.02
## 14    EE                   -0.01
## 15    IL                   -0.01
## 16    AT                    0.00
## 17    SI                    0.00
## 18    BE                    0.01
## 19    NO                    0.01
## 20    NL                    0.02
## 21    BG                    0.04
## 22    DE                    0.05
## 23    IE                    0.05
## 24    SK                    0.05
## 25    CZ                    0.11
## 26    IT                    0.12
## 27    FI                    0.17
## 28    SE                    0.17
## 29    HU                    0.18
## 30    IS                    0.20
```


#### Country-level predictors: GEI


``` r
mod6_GEI<-lmer(ben.z~gndr.c+year.c+
             gndr.c:year.c+
               gei.z.cm:gndr.c+gei.z.cm:year.c+gei.z.cm:gndr.c:year.c+
               (gndr.c+year.c+gndr.c:year.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_GEI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ben.z ~ gndr.c + year.c + gndr.c:year.c + gei.z.cm:gndr.c + gei.z.cm:year.c +  
##     gei.z.cm:gndr.c:year.c + (gndr.c + year.c + gndr.c:year.c |      cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1410453.1 1410652.6 -705208.5 1410417.1    480346 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.7381 -0.5182  0.0777  0.6605  5.3876 
## 
## Random effects:
##  Groups   Name          Variance  Std.Dev. Corr              
##  cntry    (Intercept)   8.413e-01 0.917223                   
##           gndr.c        1.169e-02 0.108100 -0.25             
##           year.c        2.779e-03 0.052719 -0.95  0.30       
##           gndr.c:year.c 2.893e-05 0.005379  0.44 -0.88 -0.47 
##  Residual               9.729e-01 0.986337                   
## Number of obs: 480364, groups:  cntry, 33
## 
## Fixed effects:
##                          Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)             0.1132653  0.1598055 31.7425245   0.709 0.483641    
## gndr.c                 -0.2197321  0.0203880 15.8153158 -10.778 1.08e-08 ***
## year.c                 -0.0031331  0.0091887 31.8535698  -0.341 0.735362    
## gndr.c:year.c           0.0008770  0.0011023 12.0146801   0.796 0.441653    
## gndr.c:gei.z.cm        -0.0550004  0.0198292 29.2540925  -2.774 0.009551 ** 
## year.c:gei.z.cm         0.0116747  0.0028877 35.3952864   4.043 0.000272 ***
## gndr.c:year.c:gei.z.cm  0.0004444  0.0011138 37.9447258   0.399 0.692145    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c year.c gnd.:. gn.:.. yr.:..
## gndr.c      -0.234                                   
## year.c      -0.953  0.278                            
## gndr.c:yr.c  0.375 -0.866 -0.401                     
## gndr.c:g.z.  0.000 -0.042  0.000  0.057              
## yr.c:g.z.cm  0.001  0.000 -0.005 -0.001  0.185       
## gndr.c:.:..  0.000  0.062 -0.001 -0.118 -0.817 -0.142
```

``` r
getFE(mod6_GEI,round=2,p.round = 5)
```

```
##                         Est.   SE    df      t       p    LL    UL
## (Intercept)             0.11 0.16 31.74   0.71 0.48364 -0.21  0.44
## gndr.c                 -0.22 0.02 15.82 -10.78 0.00000 -0.26 -0.18
## year.c                  0.00 0.01 31.85  -0.34 0.73536 -0.02  0.02
## gndr.c:year.c           0.00 0.00 12.01   0.80 0.44165  0.00  0.00
## gndr.c:gei.z.cm        -0.06 0.02 29.25  -2.77 0.00955 -0.10 -0.01
## year.c:gei.z.cm         0.01 0.00 35.40   4.04 0.00027  0.01  0.02
## gndr.c:year.c:gei.z.cm  0.00 0.00 37.94   0.40 0.69215  0.00  0.00
```

``` r
getVC(mod6_GEI)
```

```
##         grp          var1          var2 sdcor  vcov
## 1     cntry   (Intercept)          <NA>  0.92  0.84
## 2     cntry        gndr.c          <NA>  0.11  0.01
## 3     cntry        year.c          <NA>  0.05  0.00
## 4     cntry gndr.c:year.c          <NA>  0.01  0.00
## 5     cntry   (Intercept)        gndr.c -0.25 -0.03
## 6     cntry   (Intercept)        year.c -0.95 -0.05
## 7     cntry   (Intercept) gndr.c:year.c  0.44  0.00
## 8     cntry        gndr.c        year.c  0.30  0.00
## 9     cntry        gndr.c gndr.c:year.c -0.88  0.00
## 10    cntry        year.c gndr.c:year.c -0.47  0.00
## 11 Residual          <NA>          <NA>  0.99  0.97
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_GEI,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] 4.998311
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_GEI,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] 42.25532
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
##      21 -0.1977 0.0882 37.6  -0.3764   -0.019  -2.241  0.0310
##       0  0.1133 0.1600 31.7  -0.2124    0.439   0.709  0.4836
## 
## gei.z.cm =  0:
##  year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21  0.0475 0.0633 32.9  -0.0814    0.176   0.750  0.4589
##       0  0.1133 0.1600 31.7  -0.2124    0.439   0.709  0.4836
## 
## gei.z.cm =  1:
##  year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21  0.2926 0.0871 35.9   0.1159    0.469   3.358  0.0019
##       0  0.1133 0.1600 31.7  -0.2124    0.439   0.709  0.4836
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
##  year.c21 - year.c0  -0.3110 0.203 37.7   -0.721   0.0992  -1.535  0.1331
## 
## gei.z.cm =  0:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0  -0.0658 0.193 31.9   -0.459   0.3273  -0.341  0.7354
## 
## gei.z.cm =  1:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0   0.1794 0.202 37.3   -0.230   0.5885   0.888  0.3802
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
##    -0.5     21 -0.1199 0.0890 38.2 -0.30008   0.0603  -1.346  0.1861
##     0.5     21 -0.2755 0.0885 37.4 -0.45470  -0.0963  -3.115  0.0035
##    -0.5      0  0.1956 0.1630 30.9 -0.13647   0.5277   1.202  0.2387
##     0.5      0  0.0309 0.1580 31.0 -0.29150   0.3533   0.195  0.8463
## 
## gei.z.cm =  0:
##  gndr.c year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21  0.1481 0.0642 32.8  0.01755   0.2787   2.308  0.0274
##     0.5     21 -0.0532 0.0630 32.9 -0.18143   0.0751  -0.844  0.4048
##    -0.5      0  0.2231 0.1620 30.7 -0.10841   0.5547   1.373  0.1796
##     0.5      0  0.0034 0.1580 30.7 -0.31842   0.3252   0.022  0.9829
## 
## gei.z.cm =  1:
##  gndr.c year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21  0.4161 0.0878 36.3  0.23807   0.5942   4.739 <0.0001
##     0.5     21  0.1691 0.0872 35.5 -0.00784   0.3461   1.939  0.0604
##    -0.5      0  0.2506 0.1630 30.9 -0.08139   0.5826   1.540  0.1338
##     0.5      0 -0.0241 0.1580 30.9 -0.34640   0.2982  -0.153  0.8798
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
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21     0.1556 0.0192 34.6   0.1167   0.1946   8.121 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)  -0.3155 0.2080 36.8  -0.7373   0.1063  -1.516  0.1381
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     -0.1508 0.2020 36.9  -0.5595   0.2580  -0.747  0.4595
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)     -0.4712 0.2050 37.3  -0.8861  -0.0562  -2.300  0.0272
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0        -0.3064 0.1980 36.7  -0.7085   0.0956  -1.545  0.1310
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0       0.1647 0.0290 23.7   0.1048   0.2247   5.675 <0.0001
## 
## gei.z.cm =  0:
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21     0.2013 0.0116 28.4   0.1776   0.2250  17.366 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)  -0.0750 0.1980 31.0  -0.4786   0.3286  -0.379  0.7073
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0      0.1447 0.1920 31.3  -0.2473   0.5368   0.753  0.4573
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)     -0.2763 0.1940 31.2  -0.6726   0.1199  -1.422  0.1650
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0        -0.0566 0.1890 30.9  -0.4413   0.3282  -0.300  0.7662
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0       0.2197 0.0204 15.8   0.1765   0.2630  10.778 <0.0001
## 
## gei.z.cm =  1:
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21     0.2470 0.0163 29.1   0.2136   0.2804  15.130 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)   0.1655 0.2070 36.3  -0.2551   0.5861   0.798  0.4301
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0      0.4402 0.2010 36.4   0.0326   0.8479   2.189  0.0351
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)     -0.0815 0.2040 36.9  -0.4954   0.3324  -0.399  0.6922
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0         0.1932 0.1980 36.1  -0.2074   0.5939   0.978  0.3346
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0       0.2747 0.0278 20.9   0.2168   0.3327   9.868 <0.0001
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
##  diff_ESS11   -0.156 0.0192 34.6   -0.195   -0.117  -8.121 <0.0001
##  diff_ESS1    -0.165 0.0290 23.7   -0.225   -0.105  -5.675 <0.0001
## 
## gei.z.cm =  0:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11   -0.201 0.0116 28.4   -0.225   -0.178 -17.366 <0.0001
##  diff_ESS1    -0.220 0.0204 15.8   -0.263   -0.176 -10.778 <0.0001
## 
## gei.z.cm =  1:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11   -0.247 0.0163 29.1   -0.280   -0.214 -15.130 <0.0001
##  diff_ESS1    -0.275 0.0278 20.9   -0.333   -0.217  -9.868 <0.0001
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
##  diff_ESS11 - diff_ESS1  0.00909 0.0348 27.4  -0.0623   0.0804   0.261  0.7960
## 
## gei.z.cm =  0:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1  0.01842 0.0231 12.0  -0.0320   0.0688   0.796  0.4417
## 
## gei.z.cm =  1:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1  0.02775 0.0309 18.9  -0.0370   0.0925   0.898  0.3805
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```


#### Country-level predictors: GGGI


``` r
mod6_GGGI<-lmer(ben.z~gndr.c+year.c+
             gndr.c:year.c+
               gggi.z.cm:gndr.c+gggi.z.cm:year.c+gggi.z.cm:gndr.c:year.c+
               (gndr.c+year.c+gndr.c:year.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_GGGI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ben.z ~ gndr.c + year.c + gndr.c:year.c + gggi.z.cm:gndr.c +  
##     gggi.z.cm:year.c + gggi.z.cm:gndr.c:year.c + (gndr.c + year.c +      gndr.c:year.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1077483.4 1077677.9 -538723.7 1077447.4    363834 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.6354 -0.5144  0.0707  0.6649  5.3324 
## 
## Random effects:
##  Groups   Name          Variance  Std.Dev. Corr              
##  cntry    (Intercept)   1.704e+00 1.30550                    
##           gndr.c        9.240e-03 0.09612   0.38             
##           year.c        5.085e-03 0.07131  -0.97 -0.35       
##           gndr.c:year.c 4.045e-05 0.00636  -0.06 -0.80  0.08 
##  Residual               9.916e-01 0.99580                    
## Number of obs: 363852, groups:  cntry, 34
## 
## Fixed effects:
##                          Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)              0.203371   0.226638 31.399555   0.897  0.37637    
## gndr.c                  -0.176469   0.019034  9.138482  -9.271 6.02e-06 ***
## year.c                  -0.009094   0.012369 31.946012  -0.735  0.46758    
## gndr.c:year.c           -0.002299   0.001305 31.834315  -1.762  0.08771 .  
## gndr.c:gggi.z.cm        -0.085820   0.019131 37.346137  -4.486 6.71e-05 ***
## year.c:gggi.z.cm         0.011021   0.003160 32.772129   3.487  0.00141 ** 
## gndr.c:year.c:gggi.z.cm  0.003002   0.001427 40.767535   2.103  0.04170 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c year.c gnd.:. gn.:.. yr.:..
## gndr.c       0.332                                   
## year.c      -0.972 -0.314                            
## gndr.c:yr.c -0.060 -0.828  0.071                     
## gndr.c:gg..  0.000 -0.011  0.000  0.021              
## yr.c:ggg.z.  0.008  0.003 -0.026 -0.007  0.036       
## gndr.c:.:..  0.001  0.019 -0.003 -0.046 -0.859  0.062
```

``` r
getFE(mod6_GGGI,round=2,p.round = 5)
```

```
##                          Est.   SE    df     t       p    LL    UL
## (Intercept)              0.20 0.23 31.40  0.90 0.37637 -0.26  0.67
## gndr.c                  -0.18 0.02  9.14 -9.27 0.00001 -0.22 -0.13
## year.c                  -0.01 0.01 31.95 -0.74 0.46758 -0.03  0.02
## gndr.c:year.c            0.00 0.00 31.83 -1.76 0.08771  0.00  0.00
## gndr.c:gggi.z.cm        -0.09 0.02 37.35 -4.49 0.00007 -0.12 -0.05
## year.c:gggi.z.cm         0.01 0.00 32.77  3.49 0.00141  0.00  0.02
## gndr.c:year.c:gggi.z.cm  0.00 0.00 40.77  2.10 0.04170  0.00  0.01
```

``` r
getVC(mod6_GGGI)
```

```
##         grp          var1          var2 sdcor  vcov
## 1     cntry   (Intercept)          <NA>  1.31  1.70
## 2     cntry        gndr.c          <NA>  0.10  0.01
## 3     cntry        year.c          <NA>  0.07  0.01
## 4     cntry gndr.c:year.c          <NA>  0.01  0.00
## 5     cntry   (Intercept)        gndr.c  0.38  0.05
## 6     cntry   (Intercept)        year.c -0.97 -0.09
## 7     cntry   (Intercept) gndr.c:year.c -0.06  0.00
## 8     cntry        gndr.c        year.c -0.35  0.00
## 9     cntry        gndr.c gndr.c:year.c -0.80  0.00
## 10    cntry        year.c gndr.c:year.c  0.08  0.00
## 11 Residual          <NA>          <NA>  1.00  0.99
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_GGGI,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] -73.7984
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_GGGI,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] 19.2624
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
##      21 -0.2190 0.0970 35.6  -0.4158  -0.0223  -2.259  0.0301
##       0  0.2034 0.2270 31.4  -0.2586   0.6654   0.897  0.3764
## 
## gggi.z.cm =  0:
##  year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21  0.0124 0.0659 33.0  -0.1216   0.1464   0.188  0.8519
##       0  0.2034 0.2270 31.4  -0.2586   0.6654   0.897  0.3764
## 
## gggi.z.cm =  1:
##  year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21  0.2438 0.0899 35.8   0.0614   0.4263   2.711  0.0102
##       0  0.2034 0.2270 31.4  -0.2586   0.6654   0.897  0.3764
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
##  year.c21 - year.c0  -0.4224 0.270 36.6   -0.969    0.124  -1.566  0.1260
## 
## gggi.z.cm =  0:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0  -0.1910 0.260 31.9   -0.720    0.338  -0.735  0.4676
## 
## gggi.z.cm =  1:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0   0.0405 0.266 35.1   -0.500    0.581   0.152  0.8801
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
##    -0.5     21 -0.1181 0.0966 36.1 -0.31403   0.0779  -1.222  0.2297
##     0.5     21 -0.3200 0.0988 35.0 -0.52055  -0.1195  -3.240  0.0026
##    -0.5      0  0.2487 0.2240 30.3 -0.20831   0.7057   1.111  0.2753
##     0.5      0  0.1580 0.2300 30.2 -0.31190   0.6280   0.687  0.4976
## 
## gggi.z.cm =  0:
##  gndr.c year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21  0.1248 0.0662 33.0 -0.00999   0.2595   1.884  0.0684
##     0.5     21 -0.1000 0.0665 32.9 -0.23525   0.0353  -1.504  0.1422
##    -0.5      0  0.2916 0.2240 30.2 -0.16505   0.7483   1.304  0.2022
##     0.5      0  0.1151 0.2300 30.1 -0.35446   0.5847   0.501  0.6203
## 
## gggi.z.cm =  1:
##  gndr.c year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21  0.3676 0.0897 36.2  0.18571   0.5495   4.098  0.0002
##     0.5     21  0.1201 0.0915 35.1 -0.06566   0.3058   1.312  0.1980
##    -0.5      0  0.3345 0.2240 30.3 -0.12249   0.7915   1.494  0.1454
##     0.5      0  0.0722 0.2300 30.2 -0.39770   0.5421   0.314  0.7558
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
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21    0.20196 0.0239 33.58   0.1534   0.2505   8.464 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0) -0.36676 0.2690 35.48  -0.9133   0.1798  -1.362  0.1819
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0    -0.27611 0.2750 35.33  -0.8332   0.2810  -1.006  0.3214
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)    -0.56872 0.2660 35.93  -1.1079  -0.0296  -2.140  0.0393
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0       -0.47807 0.2720 35.45  -1.0297   0.0735  -1.759  0.0872
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0      0.09065 0.0271 22.76   0.0345   0.1468   3.341  0.0029
## 
## gggi.z.cm =  0:
##  contrast                                     estimate     SE    df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21    0.22475 0.0158 31.71   0.1925   0.2569  14.223 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0) -0.16684 0.2590 30.90  -0.6954   0.3618  -0.644  0.5244
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     0.00963 0.2650 31.09  -0.5306   0.5499   0.036  0.9712
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)    -0.39158 0.2550 30.97  -0.9116   0.1284  -1.536  0.1347
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0       -0.21512 0.2610 30.74  -0.7478   0.3176  -0.824  0.4163
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0      0.17647 0.0190  9.14   0.1335   0.2194   9.271 <0.0001
## 
## gggi.z.cm =  1:
##  contrast                                     estimate     SE    df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21    0.24753 0.0221 33.60   0.2025   0.2925  11.185 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)  0.03308 0.2660 33.95  -0.5076   0.5737   0.124  0.9018
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     0.29537 0.2710 33.89  -0.2560   0.8467   1.089  0.2839
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)    -0.21445 0.2620 34.27  -0.7472   0.3184  -0.818  0.4192
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0        0.04784 0.2680 33.88  -0.4974   0.5931   0.178  0.8595
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0      0.26229 0.0268 21.24   0.2065   0.3181   9.773 <0.0001
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
##  diff_ESS11  -0.2020 0.0239 33.58   -0.250  -0.1534  -8.464 <0.0001
##  diff_ESS1   -0.0906 0.0271 22.76   -0.147  -0.0345  -3.341  0.0029
## 
## gggi.z.cm =  0:
##  contrast   estimate     SE    df lower.CL upper.CL t.ratio p.value
##  diff_ESS11  -0.2247 0.0158 31.71   -0.257  -0.1925 -14.223 <0.0001
##  diff_ESS1   -0.1765 0.0190  9.14   -0.219  -0.1335  -9.271 <0.0001
## 
## gggi.z.cm =  1:
##  contrast   estimate     SE    df lower.CL upper.CL t.ratio p.value
##  diff_ESS11  -0.2475 0.0221 33.60   -0.293  -0.2025 -11.185 <0.0001
##  diff_ESS1   -0.2623 0.0268 21.24   -0.318  -0.2065  -9.773 <0.0001
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
##  diff_ESS11 - diff_ESS1  -0.1113 0.0415 38.5  -0.1953 -0.02728  -2.680  0.0108
## 
## gggi.z.cm =  0:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1  -0.0483 0.0274 31.8  -0.1041  0.00755  -1.762  0.0877
## 
## gggi.z.cm =  1:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1   0.0148 0.0397 35.1  -0.0658  0.09530   0.372  0.7122
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

#### Country-level predictors: GDI


``` r
mod6_GDI<-lmer(ben.z~gndr.c+year.c+
             gndr.c:year.c+
               gdi.z.cm:gndr.c+gdi.z.cm:year.c+gdi.z.cm:gndr.c:year.c+
               (gndr.c+year.c+gndr.c:year.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_GDI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ben.z ~ gndr.c + year.c + gndr.c:year.c + gdi.z.cm:gndr.c + gdi.z.cm:year.c +  
##     gdi.z.cm:gndr.c:year.c + (gndr.c + year.c + gndr.c:year.c |      cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1455828.6 1456028.5 -727896.3 1455792.6    492325 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.6532 -0.5185  0.0761  0.6582  5.3252 
## 
## Random effects:
##  Groups   Name          Variance  Std.Dev. Corr              
##  cntry    (Intercept)   0.8288570 0.910416                   
##           gndr.c        0.0161884 0.127234 -0.21             
##           year.c        0.0031263 0.055913 -0.95  0.11       
##           gndr.c:year.c 0.0000471 0.006863  0.34 -0.84 -0.29 
##  Residual               0.9890945 0.994532                   
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##                          Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)             9.287e-02  1.563e-01  3.264e+01   0.594   0.5564    
## gndr.c                 -2.111e-01  2.330e-02  2.262e+01  -9.058 5.54e-09 ***
## year.c                 -2.664e-03  9.600e-03  3.242e+01  -0.277   0.7832    
## gndr.c:year.c          -4.535e-05  1.318e-03  2.007e+01  -0.034   0.9729    
## gndr.c:gdi.z.cm        -1.393e-02  2.334e-02  3.446e+01  -0.597   0.5545    
## year.c:gdi.z.cm        -8.520e-03  3.198e-03  3.468e+01  -2.664   0.0116 *  
## gndr.c:year.c:gdi.z.cm -1.089e-03  1.395e-03  4.373e+01  -0.781   0.4391    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c year.c gnd.:. gn.:.. yr.:..
## gndr.c      -0.200                                   
## year.c      -0.946  0.104                            
## gndr.c:yr.c  0.299 -0.841 -0.263                     
## gndr.c:gd..  0.000 -0.026  0.000  0.027              
## yr.c:gd.z.c  0.001 -0.001 -0.003  0.002 -0.262       
## gndr.c:.:..  0.000  0.026  0.001 -0.046 -0.807  0.056
```

``` r
getFE(mod6_GDI,round=2,p.round = 5)
```

```
##                         Est.   SE    df     t       p    LL    UL
## (Intercept)             0.09 0.16 32.64  0.59 0.55642 -0.23  0.41
## gndr.c                 -0.21 0.02 22.62 -9.06 0.00000 -0.26 -0.16
## year.c                  0.00 0.01 32.42 -0.28 0.78318 -0.02  0.02
## gndr.c:year.c           0.00 0.00 20.07 -0.03 0.97290  0.00  0.00
## gndr.c:gdi.z.cm        -0.01 0.02 34.46 -0.60 0.55450 -0.06  0.03
## year.c:gdi.z.cm        -0.01 0.00 34.68 -2.66 0.01163 -0.02  0.00
## gndr.c:year.c:gdi.z.cm  0.00 0.00 43.73 -0.78 0.43912  0.00  0.00
```

``` r
getVC(mod6_GDI)
```

```
##         grp          var1          var2 sdcor  vcov
## 1     cntry   (Intercept)          <NA>  0.91  0.83
## 2     cntry        gndr.c          <NA>  0.13  0.02
## 3     cntry        year.c          <NA>  0.06  0.00
## 4     cntry gndr.c:year.c          <NA>  0.01  0.00
## 5     cntry   (Intercept)        gndr.c -0.21 -0.02
## 6     cntry   (Intercept)        year.c -0.95 -0.05
## 7     cntry   (Intercept) gndr.c:year.c  0.34  0.00
## 8     cntry        gndr.c        year.c  0.11  0.00
## 9     cntry        gndr.c gndr.c:year.c -0.84  0.00
## 10    cntry        year.c gndr.c:year.c -0.29  0.00
## 11 Residual          <NA>          <NA>  0.99  0.99
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_GDI,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] -6.862678
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_GDI,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] 5.989397
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
##      21  0.2159 0.1000 38.3    0.013   0.4187   2.154  0.0376
##       0  0.0929 0.1560 32.6   -0.225   0.4109   0.594  0.5564
## 
## gdi.z.cm =  0:
##  year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21  0.0369 0.0738 32.4   -0.113   0.1872   0.500  0.6203
##       0  0.0929 0.1560 32.6   -0.225   0.4109   0.594  0.5564
## 
## gdi.z.cm =  1:
##  year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21 -0.1420 0.0994 37.1   -0.343   0.0594  -1.428  0.1615
##       0  0.0929 0.1560 32.6   -0.225   0.4109   0.594  0.5564
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
##  year.c21 - year.c0   0.1230 0.213 38.9   -0.307    0.553   0.578  0.5665
## 
## gdi.z.cm =  0:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0  -0.0559 0.202 32.4   -0.466    0.354  -0.277  0.7832
## 
## gdi.z.cm =  1:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0  -0.2349 0.212 38.6   -0.664    0.195  -1.106  0.2754
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
##    -0.5     21  0.3035 0.1040 38.3   0.0922   0.5148   2.907  0.0060
##     0.5     21  0.1283 0.0973 38.2  -0.0686   0.3251   1.318  0.1952
##    -0.5      0  0.1914 0.1590 31.8  -0.1335   0.5164   1.200  0.2388
##     0.5      0 -0.0057 0.1550 31.8  -0.3211   0.3097  -0.037  0.9709
## 
## gdi.z.cm =  0:
##  gndr.c year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21  0.1429 0.0769 32.3  -0.0137   0.2996   1.858  0.0723
##     0.5     21 -0.0691 0.0714 32.4  -0.2144   0.0762  -0.968  0.3403
##    -0.5      0  0.1984 0.1590 31.4  -0.1257   0.5225   1.248  0.2214
##     0.5      0 -0.0127 0.1540 31.4  -0.3273   0.3020  -0.082  0.9351
## 
## gdi.z.cm =  1:
##  gndr.c year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21 -0.0176 0.1040 37.1  -0.2273   0.1921  -0.170  0.8660
##     0.5     21 -0.2664 0.0964 37.0  -0.4618  -0.0710  -2.763  0.0089
##    -0.5      0  0.2054 0.1590 31.7  -0.1195   0.5302   1.288  0.2070
##     0.5      0 -0.0196 0.1550 31.8  -0.3350   0.2957  -0.127  0.8999
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
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21    0.17521 0.0236 31.6    0.127   0.2233   7.428 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)  0.11203 0.2170 37.8   -0.327   0.5511   0.517  0.6085
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     0.30916 0.2160 39.1   -0.129   0.7470   1.428  0.1612
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)    -0.06319 0.2110 37.6   -0.490   0.3638  -0.300  0.7661
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0        0.13395 0.2100 38.3   -0.292   0.5599   0.636  0.5283
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0      0.19714 0.0334 29.8    0.129   0.2654   5.901 <0.0001
## 
## gdi.z.cm =  0:
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21    0.21202 0.0150 28.4    0.181   0.2427  14.153 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0) -0.05546 0.2060 31.6   -0.475   0.3637  -0.270  0.7892
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     0.15561 0.2030 31.8   -0.259   0.5697   0.765  0.4496
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)    -0.26748 0.2010 31.8   -0.677   0.1417  -1.332  0.1924
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0       -0.05641 0.1980 31.5   -0.461   0.3480  -0.284  0.7780
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0      0.21107 0.0233 22.6    0.163   0.2593   9.058 <0.0001
## 
## gdi.z.cm =  1:
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21    0.24883 0.0222 32.6    0.204   0.2940  11.220 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0) -0.22295 0.2160 37.4   -0.661   0.2152  -1.031  0.3093
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     0.00205 0.2160 38.7   -0.435   0.4390   0.009  0.9925
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)    -0.47178 0.2100 37.3   -0.898  -0.0457  -2.243  0.0309
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0       -0.24678 0.2100 37.9   -0.672   0.1784  -1.175  0.2472
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0      0.22500 0.0326 29.4    0.158   0.2915   6.912 <0.0001
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
##  diff_ESS11   -0.175 0.0236 31.6   -0.223   -0.127  -7.428 <0.0001
##  diff_ESS1    -0.197 0.0334 29.8   -0.265   -0.129  -5.901 <0.0001
## 
## gdi.z.cm =  0:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11   -0.212 0.0150 28.4   -0.243   -0.181 -14.153 <0.0001
##  diff_ESS1    -0.211 0.0233 22.6   -0.259   -0.163  -9.058 <0.0001
## 
## gdi.z.cm =  1:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11   -0.249 0.0222 32.6   -0.294   -0.204 -11.220 <0.0001
##  diff_ESS1    -0.225 0.0326 29.4   -0.292   -0.158  -6.912 <0.0001
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
##  contrast                estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1  0.021924 0.0412 34.5  -0.0618   0.1057   0.532  0.5982
## 
## gdi.z.cm =  0:
##  contrast                estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1 -0.000952 0.0277 20.1  -0.0587   0.0568  -0.034  0.9729
## 
## gdi.z.cm =  1:
##  contrast                estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1 -0.023829 0.0394 32.6  -0.1040   0.0563  -0.605  0.5492
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

#### Country-level predictors: log_GDP


``` r
mod6_log_GDP<-lmer(ben.z~gndr.c+year.c+
             gndr.c:year.c+
               log_gdp.z.cm:gndr.c+log_gdp.z.cm:year.c+log_gdp.z.cm:gndr.c:year.c+
               (gndr.c+year.c+gndr.c:year.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_log_GDP)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ben.z ~ gndr.c + year.c + gndr.c:year.c + log_gdp.z.cm:gndr.c +  
##     log_gdp.z.cm:year.c + log_gdp.z.cm:gndr.c:year.c + (gndr.c +      year.c + gndr.c:year.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1455816.7 1456016.7 -727890.4 1455780.7    492325 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -9.6546 -0.5182  0.0760  0.6585  5.3345 
## 
## Random effects:
##  Groups   Name          Variance  Std.Dev. Corr              
##  cntry    (Intercept)   0.8278051 0.909838                   
##           gndr.c        0.0118749 0.108972 -0.29             
##           year.c        0.0025281 0.050280 -0.95  0.35       
##           gndr.c:year.c 0.0000443 0.006656  0.44 -0.81 -0.49 
##  Residual               0.9890982 0.994534                   
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##                              Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)                 0.0932864  0.1561713 32.5828384   0.597 0.554413    
## gndr.c                     -0.2044588  0.0203262 14.5634861 -10.059 6.07e-08 ***
## year.c                     -0.0024808  0.0086348 32.4094054  -0.287 0.775713    
## gndr.c:year.c              -0.0005551  0.0012846 16.3976840  -0.432 0.671295    
## gndr.c:log_gdp.z.cm        -0.0755174  0.0207169 29.3697938  -3.645 0.001025 ** 
## year.c:log_gdp.z.cm         0.0121231  0.0027582 34.3508850   4.395 0.000101 ***
## gndr.c:year.c:log_gdp.z.cm  0.0034019  0.0012417 35.5844421   2.740 0.009549 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c year.c gnd.:. g.:_.. y.:_..
## gndr.c      -0.265                                   
## year.c      -0.948  0.317                            
## gndr.c:yr.c  0.387 -0.819 -0.432                     
## gndr.c:l_..  0.000 -0.084  0.002  0.081              
## yr.c:lg_g..  0.001  0.005  0.005 -0.005  0.203       
## gndr.:.:_..  0.001  0.084 -0.002 -0.089 -0.819 -0.210
```

``` r
getFE(mod6_log_GDP,round=2,p.round = 5)
```

```
##                             Est.   SE    df      t       p    LL    UL
## (Intercept)                 0.09 0.16 32.58   0.60 0.55441 -0.22  0.41
## gndr.c                     -0.20 0.02 14.56 -10.06 0.00000 -0.25 -0.16
## year.c                      0.00 0.01 32.41  -0.29 0.77571 -0.02  0.02
## gndr.c:year.c               0.00 0.00 16.40  -0.43 0.67129  0.00  0.00
## gndr.c:log_gdp.z.cm        -0.08 0.02 29.37  -3.65 0.00102 -0.12 -0.03
## year.c:log_gdp.z.cm         0.01 0.00 34.35   4.40 0.00010  0.01  0.02
## gndr.c:year.c:log_gdp.z.cm  0.00 0.00 35.58   2.74 0.00955  0.00  0.01
```

``` r
getVC(mod6_log_GDP)
```

```
##         grp          var1          var2 sdcor  vcov
## 1     cntry   (Intercept)          <NA>  0.91  0.83
## 2     cntry        gndr.c          <NA>  0.11  0.01
## 3     cntry        year.c          <NA>  0.05  0.00
## 4     cntry gndr.c:year.c          <NA>  0.01  0.00
## 5     cntry   (Intercept)        gndr.c -0.29 -0.03
## 6     cntry   (Intercept)        year.c -0.95 -0.04
## 7     cntry   (Intercept) gndr.c:year.c  0.44  0.00
## 8     cntry        gndr.c        year.c  0.35  0.00
## 9     cntry        gndr.c gndr.c:year.c -0.81  0.00
## 10    cntry        year.c gndr.c:year.c -0.49  0.00
## 11 Residual          <NA>          <NA>  0.99  0.99
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_log_GDP,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] 13.58516
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_log_GDP,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] 11.58558
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
##      21 -0.2134 0.0825 35.8  -0.3808   -0.046  -2.586  0.0139
##       0  0.0933 0.1560 32.6  -0.2246    0.411   0.597  0.5544
## 
## log_gdp.z.cm =  0:
##  year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21  0.0412 0.0598 33.9  -0.0803    0.163   0.689  0.4954
##       0  0.0933 0.1560 32.6  -0.2246    0.411   0.597  0.5544
## 
## log_gdp.z.cm =  1:
##  year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21  0.2958 0.0840 35.1   0.1254    0.466   3.523  0.0012
##       0  0.0933 0.1560 32.6  -0.2246    0.411   0.597  0.5544
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
##  year.c21 - year.c0  -0.3067 0.190 38.1   -0.691   0.0781  -1.614  0.1149
## 
## log_gdp.z.cm =  0:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0  -0.0521 0.181 32.4   -0.421   0.3171  -0.287  0.7757
## 
## log_gdp.z.cm =  1:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0   0.2025 0.191 38.4   -0.183   0.5883   1.062  0.2948
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
##  gndr.c year.c   emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21 -0.10738 0.0845 36.1  -0.2788   0.0641  -1.270  0.2122
##     0.5     21 -0.31941 0.0819 35.5  -0.4856  -0.1532  -3.900  0.0004
##    -0.5      0  0.15776 0.1600 31.8  -0.1673   0.4828   0.989  0.3302
##     0.5      0  0.02882 0.1540 31.9  -0.2853   0.3429   0.187  0.8529
## 
## log_gdp.z.cm =  0:
##  gndr.c year.c   emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21  0.14925 0.0615 33.8   0.0242   0.2743   2.426  0.0208
##     0.5     21 -0.06687 0.0590 33.9  -0.1868   0.0531  -1.133  0.2650
##    -0.5      0  0.19552 0.1590 31.5  -0.1289   0.5199   1.228  0.2284
##     0.5      0 -0.00894 0.1540 31.6  -0.3224   0.3045  -0.058  0.9540
## 
## log_gdp.z.cm =  1:
##  gndr.c year.c   emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21  0.40587 0.0860 35.3   0.2314   0.5803   4.722 <0.0001
##     0.5     21  0.18568 0.0833 34.7   0.0166   0.3548   2.230  0.0323
##    -0.5      0  0.23327 0.1590 31.7  -0.0916   0.5582   1.463  0.1533
##     0.5      0 -0.04670 0.1540 31.8  -0.3606   0.2672  -0.303  0.7638
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
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21     0.2120 0.0220 33.2   0.1673   0.2568   9.637 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)  -0.2651 0.1970 37.4  -0.6647   0.1345  -1.344  0.1871
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     -0.1362 0.1900 37.4  -0.5208   0.2484  -0.717  0.4777
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)     -0.4772 0.1920 37.6  -0.8658  -0.0885  -2.486  0.0175
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0        -0.3482 0.1850 36.9  -0.7224   0.0260  -1.886  0.0672
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0       0.1289 0.0302 24.0   0.0666   0.1913   4.267  0.0003
## 
## log_gdp.z.cm =  0:
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21     0.2161 0.0156 30.0   0.1843   0.2480  13.858 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)  -0.0463 0.1880 31.5  -0.4285   0.3360  -0.247  0.8067
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0      0.1582 0.1810 31.8  -0.2108   0.5272   0.873  0.3890
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)     -0.2624 0.1820 31.7  -0.6339   0.1092  -1.439  0.1600
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0        -0.0579 0.1760 31.4  -0.4166   0.3007  -0.329  0.7442
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0       0.2045 0.0203 14.6   0.1610   0.2479  10.059 <0.0001
## 
## log_gdp.z.cm =  1:
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21     0.2202 0.0213 29.9   0.1768   0.2636  10.360 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)   0.1726 0.1980 37.5  -0.2279   0.5731   0.873  0.3884
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0      0.4526 0.1900 37.6   0.0671   0.8381   2.378  0.0226
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)     -0.0476 0.1920 37.9  -0.4373   0.3421  -0.247  0.8060
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0         0.2324 0.1850 37.0  -0.1424   0.6072   1.256  0.2169
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0       0.2800 0.0278 19.6   0.2219   0.3380  10.079 <0.0001
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
##  diff_ESS11   -0.212 0.0220 33.2   -0.257  -0.1673  -9.637 <0.0001
##  diff_ESS1    -0.129 0.0302 24.0   -0.191  -0.0666  -4.267  0.0003
## 
## log_gdp.z.cm =  0:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11   -0.216 0.0156 30.0   -0.248  -0.1843 -13.858 <0.0001
##  diff_ESS1    -0.204 0.0203 14.6   -0.248  -0.1610 -10.059 <0.0001
## 
## log_gdp.z.cm =  1:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11   -0.220 0.0213 29.9   -0.264  -0.1768 -10.360 <0.0001
##  diff_ESS1    -0.280 0.0278 19.6   -0.338  -0.2219 -10.079 <0.0001
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
##  diff_ESS11 - diff_ESS1  -0.0831 0.0392 30.6  -0.1630  -0.0032  -2.122  0.0420
## 
## log_gdp.z.cm =  0:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1  -0.0117 0.0270 16.4  -0.0687   0.0454  -0.432  0.6713
## 
## log_gdp.z.cm =  1:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1   0.0598 0.0358 23.5  -0.0142   0.1338   1.670  0.1083
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
##   [1] mathjaxr_2.0-0     RColorBrewer_1.1-3 rstudioapi_0.18.0  jsonlite_2.0.0     shape_1.4.6.1     
##   [6] magrittr_2.0.5     estimability_1.5.1 jomo_2.7-6         farver_2.1.2       nloptr_2.2.1      
##  [11] fs_2.1.0           vctrs_0.7.3        memoise_2.0.1      minqa_1.2.8        base64enc_0.1-6   
##  [16] rstatix_0.7.3      htmltools_0.5.9    forcats_1.0.1      usethis_3.2.1      broom_1.0.13      
##  [21] cellranger_1.1.0   Formula_1.2-5      mitml_0.4-5        sass_0.4.10        bslib_0.10.0      
##  [26] htmlwidgets_1.6.4  plyr_1.8.9         cachem_1.1.0       lifecycle_1.0.5    iterators_1.0.14  
##  [31] pkgconfig_2.0.3    R6_2.6.1           fastmap_1.2.0      rbibutils_2.4.1    digest_0.6.39     
##  [36] colorspace_2.1-2   pkgload_1.5.2      labeling_0.4.3     mgcv_1.9-4         abind_1.4-8       
##  [41] compiler_4.6.0     withr_3.0.2        htmlTable_2.5.0    S7_0.2.2           backports_1.5.1   
##  [46] carData_3.0-6      psych_2.6.3        pkgbuild_1.4.8     R.utils_2.13.0     ggsignif_0.6.4    
##  [51] pan_1.9            MASS_7.3-65        sessioninfo_1.2.3  tools_4.6.0        pbivnorm_0.6.0    
##  [56] foreign_0.8-91     otel_0.2.0         zip_2.3.3          nnet_7.3-20        R.oo_1.27.1       
##  [61] glue_1.8.1         quadprog_1.5-8     grid_4.6.0         checkmate_2.3.4    cluster_2.1.8.2   
##  [66] generics_0.1.4     gtable_0.3.6       R.methodsS3_1.8.2  data.table_1.18.4  car_3.1-5         
##  [71] utf8_1.2.6         foreach_1.5.2      pillar_1.11.1      rockchalk_1.8.164  splines_4.6.0     
##  [76] lattice_0.22-9     survival_3.8-6     kutils_1.73        tidyselect_1.2.1   reformulas_0.4.4  
##  [81] gridExtra_2.3      grImport2_0.3-3    stats4_4.6.0       xfun_0.57          devtools_2.5.2    
##  [86] stringi_1.8.7      yaml_2.3.12        boot_1.3-32        evaluate_1.0.5     codetools_0.2-20  
##  [91] cli_3.6.6          rpart_4.1.27       xtable_1.8-8       Rdpack_2.6.6       jquerylib_0.1.4   
##  [96] lavaan_0.6-21      Rcpp_1.1.1-1.1     readxl_1.4.5       png_0.1-9          coda_0.19-4.1     
## [101] XML_3.99-0.23      parallel_4.6.0     ellipsis_0.3.3     jpeg_0.1-11        glmnet_5.0        
## [106] mvtnorm_1.3-7      scales_1.4.0       openxlsx_4.2.8.1   purrr_1.2.2        writexl_1.5.4     
## [111] rlang_1.2.0        cowplot_1.2.0      mnormt_2.1.2       mice_3.19.0
```

