---
title: "Analysis for stimulation values"
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
cntry.sti<-diff_dat %>% group_by(cntry,essround) %>%
  summarise(sti.ctm=mean(sti,na.rm=T),
            sti.ctsd=sd(sti,na.rm=T)) %>%
  group_by(cntry) %>%
  summarise(sti.cm=mean(sti.ctm),
            sti.csd=mean(sti.ctsd)) 
```

```
## `summarise()` has regrouped the output.
## ℹ Summaries were computed grouped by cntry and essround.
## ℹ Output is grouped by cntry.
## ℹ Use `summarise(.groups = "drop_last")` to silence this message.
## ℹ Use `summarise(.by = c(cntry, essround))` for per-operation grouping (`?dplyr::dplyr_by`) instead.
```

``` r
grand_mean_sti<-mean(cntry.sti$sti.cm)
grand_sd_sti<-mean(cntry.sti$sti.csd)

# standardized
diff_dat$sti.z<-(diff_dat$sti-grand_mean_sti)/grand_sd_sti
hist(diff_dat$sti.z)
```

![](Analysis_stimulation_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

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
                   sti.z.wt=weighted.mean(x=sti.z,w=pspwght),
                   sti.z=mean(sti.z),
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

# stimulation

cntry_sti_frame<-
  diff_dat %>% group_by(cntry,essround) %>%
  summarise('sti M' = weighted.mean(x=sti.z,w=pspwght),
            'sti SD' = sqrt(wtd.var(sti.z,pspwght))) %>%
  group_by(cntry) %>%
  summarise('sti M' = mean(x=`sti M`),
            'sti SD'= mean(x=`sti SD`))
```

```
## `summarise()` has regrouped the output.
## ℹ Summaries were computed grouped by cntry and essround.
## ℹ Output is grouped by cntry.
## ℹ Use `summarise(.groups = "drop_last")` to silence this message.
## ℹ Use `summarise(.by = c(cntry, essround))` for per-operation grouping (`?dplyr::dplyr_by`) instead.
```

``` r
cntry_sti_women_frame<-
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
## `summarise()` has regrouped the output.
## ℹ Summaries were computed grouped by cntry and essround.
## ℹ Output is grouped by cntry.
## ℹ Use `summarise(.groups = "drop_last")` to silence this message.
## ℹ Use `summarise(.by = c(cntry, essround))` for per-operation grouping (`?dplyr::dplyr_by`) instead.
```

``` r
cntry_sti_men_frame<-
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
## `summarise()` has regrouped the output.
## ℹ Summaries were computed grouped by cntry and essround.
## ℹ Output is grouped by cntry.
## ℹ Use `summarise(.groups = "drop_last")` to silence this message.
## ℹ Use `summarise(.by = c(cntry, essround))` for per-operation grouping (`?dplyr::dplyr_by`) instead.
```

``` r
# link n and sti datasets

desc_frame<-
  left_join(
    x=cntry_n_frame,
    y=cntry_sti_frame,
    by="cntry"
  )

desc_frame<-
  left_join(
    x=desc_frame,
    y=cntry_sti_women_frame,
    by="cntry"
  )

desc_frame<-
  left_join(
    x=desc_frame,
    y=cntry_sti_men_frame,
    by="cntry"
  )

# Add country-specific differences
desc_frame$D<-desc_frame$`sti M Men`-desc_frame$`sti M Women`

desc_frame
```

```
## # A tibble: 34 × 10
##    cntry `n ESS rounds`     n `sti M` `sti SD` `sti M Women` `sti SD Women` `sti M Men` `sti SD Men`
##    <chr>          <dbl> <dbl>   <dbl>    <dbl>         <dbl>          <dbl>       <dbl>        <dbl>
##  1 AT                 7 15400  0.0328    1.00       -0.0659           1.02       0.138         0.972
##  2 BE                11 18886  0.0974    0.939       0.0253           0.942      0.173         0.931
##  3 BG                 7 14857  0.154     1.05        0.0421           1.07       0.275         1.02 
##  4 CH                11 18087  0.105     0.959       0.00465          0.970      0.209         0.936
##  5 CY                 6  5771  0.181     1.08        0.0762           1.10       0.290         1.05 
##  6 CZ                 9 18934 -0.0312    0.999      -0.149            0.995      0.0965        0.985
##  7 DE                10 27753 -0.153     0.967      -0.256            0.954     -0.0438        0.969
##  8 DK                 8 12198  0.0824    1.05       -0.00190          1.06       0.169         1.04 
##  9 EE                10 17974 -0.0454    0.971      -0.126            0.975      0.0510        0.956
## 10 ES                10 18785 -0.0298    1.03       -0.0733           1.03       0.0153        1.03 
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
print(cntry_desc_tbl,n=35)
```

```
## # A tibble: 34 × 14
##    Country     `n ESS rounds`     n `sti M` `sti SD` `sti M Women` `sti SD Women` `sti M Men` `sti SD Men`
##    <chr>                <dbl> <dbl> <chr>   <chr>    <chr>         <chr>          <chr>       <chr>       
##  1 Austria                  7 15400 0.03    1.00     -0.07         1.02           0.14        0.97        
##  2 Belgium                 11 18886 0.10    0.94     0.03          0.94           0.17        0.93        
##  3 Bulgaria                 7 14857 0.15    1.05     0.04          1.07           0.27        1.02        
##  4 Switzerland             11 18087 0.10    0.96     0.00          0.97           0.21        0.94        
##  5 Cyprus                   6  5771 0.18    1.08     0.08          1.10           0.29        1.05        
##  6 Czechia                  9 18934 -0.03   1.00     -0.15         0.99           0.10        0.99        
##  7 Germany                 10 27753 -0.15   0.97     -0.26         0.95           -0.04       0.97        
##  8 Denmark                  8 12198 0.08    1.05     -0.00         1.06           0.17        1.04        
##  9 Estonia                 10 17974 -0.05   0.97     -0.13         0.98           0.05        0.96        
## 10 Spain                   10 18785 -0.03   1.03     -0.07         1.03           0.02        1.03        
## 11 Finland                 11 19568 0.08    0.97     0.02          0.99           0.15        0.94        
## 12 France                  11 20457 -0.09   0.99     -0.19         0.97           0.02        1.00        
## 13 UK                      11 22979 0.12    1.01     0.05          1.03           0.19        0.99        
## 14 Greece                   6 15212 0.23    1.00     0.13          1.01           0.34        0.97        
## 15 Croatia                  5  7914 -0.37   1.07     -0.50         1.06           -0.23       1.05        
## 16 Hungary                 11 18123 0.11    0.99     0.04          0.99           0.20        0.99        
## 17 Ireland                 11 22562 0.21    1.04     0.16          1.04           0.28        1.03        
## 18 Israel                   7 14857 0.24    1.06     0.17          1.06           0.32        1.06        
## 19 Iceland                  6  4654 0.13    1.01     0.07          1.01           0.18        1.00        
## 20 Italy                    5 11441 0.04    0.87     -0.04         0.88           0.12        0.86        
## 21 Lithuania                7 13059 -0.11   1.06     -0.24         1.06           0.06        1.04        
## 22 Latvia                   3  4088 0.17    1.00     0.09          1.00           0.27        1.00        
## 23 Montenegro               3  4028 0.06    0.97     -0.01         0.99           0.13        0.95        
## 24 Netherlands             11 19722 0.16    0.91     0.08          0.91           0.24        0.89        
## 25 Norway                  11 16505 -0.06   1.01     -0.16         1.03           0.04        0.98        
## 26 Poland                  10 16737 -0.02   0.97     -0.14         0.97           0.11        0.96        
## 27 Portugal                11 19070 -0.12   0.93     -0.23         0.94           0.01        0.90        
## 28 Serbia                   2  3499 0.04    1.08     -0.04         1.10           0.14        1.06        
## 29 Russia                   5 12139 -0.06   1.09     -0.18         1.09           0.08        1.06        
## 30 Sweden                  10 16104 -0.04   0.97     -0.11         0.99           0.03        0.96        
## 31 Slovenia                11 14463 0.30    0.94     0.22          0.95           0.39        0.93        
## 32 Slovakia                 8 12547 0.04    0.95     -0.07         0.96           0.16        0.91        
## 33 Turkey                   2  4108 0.30    1.02     0.23          1.02           0.36        1.01        
## 34 Ukraine                  6 12054 -0.22   1.05     -0.29         1.05           -0.13       1.04        
## # ℹ 5 more variables: D <chr>, GEI <chr>, GGGI <chr>, GDI <chr>, GDP <chr>
```

``` r
export(cntry_desc_tbl,"../results/sti/cntry_desc_tbl.xlsx",overwrite=T)
```

## Country-level correlations table


``` r
cor_frame<-
  desc_frame %>%
  dplyr::select(
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
  dplyr::select(-GDP)

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
##   1. VBMT       0.05  0.15                                                                          
##                                                                                                     
##   2. VBMT_Women -0.04 0.16 .99                                                                      
##                            [.98, 1.00]                                                              
##                                                                                                     
##   3. VBMT_Men   0.14  0.14 .98          .95                                                         
##                            [.97, .99]   [.91, .98]                                                  
##                                                                                                     
##   4. D          0.19  0.05 -.47         -.58         -.31                                           
##                            [-.70, -.16] [-.77, -.30] [-.59, .03]                                    
##                                                                                                     
##   5. GEI        0.87  0.07 -.24         -.20         -.30        -.15                               
##                            [-.54, .11]  [-.51, .15]  [-.58, .05] [-.47, .20]                        
##                                                                                                     
##   6. GGGI       0.74  0.05 -.09         -.05         -.16        -.30        .73                    
##                            [-.42, .25]  [-.38, .30]  [-.47, .19] [-.58, .05] [.52, .86]             
##                                                                                                     
##   7. GDI        0.98  0.03 -.38         -.39         -.32        .37         .07         .19        
##                            [-.63, -.04] [-.64, -.06] [-.59, .02] [.04, .63]  [-.28, .41] [-.16, .50]
##                                                                                                     
##   8. log_GDP    10.61 0.41 .11          .12          .07         -.19        .72         .62        
##                            [-.24, .43]  [-.23, .44]  [-.27, .40] [-.50, .16] [.50, .85]  [.36, .79] 
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
## 1460408.5 1460441.8 -730201.2 1460402.5    492340 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -4.7663 -0.7001 -0.0296  0.6497  4.8089 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  cntry    (Intercept) 0.0218   0.1476  
##  Residual             0.9993   0.9996  
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##             Estimate Std. Error       df t value Pr(>|t|)  
## (Intercept)  0.04719    0.02538 33.96126    1.86   0.0716 .
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
## mean variation  0.02134914     NA       1
## sigma2          0.97865086      1      NA
## 
## $R2s
##          total within between
## f1  0.00000000      0      NA
## f2  0.00000000     NA       0
## v   0.00000000      0      NA
## m   0.02134914     NA       1
## f   0.00000000     NA      NA
## fv  0.00000000      0      NA
## fvm 0.02134914     NA      NA
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
## 1456015.0 1456059.4 -728003.5 1456007.0    492339 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -4.7891 -0.6932 -0.0412  0.6427  5.0062 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  cntry    (Intercept) 0.02145  0.1465  
##  Residual             0.99040  0.9952  
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##              Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept) 5.075e-02  2.517e-02 3.396e+01   2.016   0.0517 .  
## gndr.c      1.883e-01  2.834e-03 4.923e+05  66.447   <2e-16 ***
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
##              Est.    SE        df      t     p    LL    UL
## (Intercept) 0.051 0.025     33.96  2.016 0.052 0.000 0.102
## gndr.c      0.188 0.003 492312.06 66.447 0.000 0.183 0.194
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
## fixed           0.008630342
## slope variation 0.000000000
## mean variation  0.021014960
## sigma2          0.970354697
## 
## $R2s
##           total
## f   0.008630342
## v   0.000000000
## m   0.021014960
## fv  0.008630342
## fvm 0.029645303
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
## 1455821.1 1455887.8 -727904.6 1455809.1    492337 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -4.8686 -0.6919 -0.0539  0.6445  5.0879 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr  
##  cntry    (Intercept) 0.021333 0.14606        
##           gndr.c      0.002241 0.04734  -0.46 
##  Residual             0.989875 0.99492        
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##              Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)  0.051033   0.025103 33.960719   2.033   0.0499 *  
## gndr.c       0.187687   0.008721 33.382701  21.522   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##        (Intr)
## gndr.c -0.429
```

``` r
getFE(mod2,round=3)
```

```
##              Est.    SE     df      t    p   LL    UL
## (Intercept) 0.051 0.025 33.961  2.033 0.05 0.00 0.102
## gndr.c      0.188 0.009 33.383 21.522 0.00 0.17 0.205
```

``` r
getVC(mod2)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.15 0.02
## 2    cntry      gndr.c   <NA>  0.05 0.00
## 3    cntry (Intercept) gndr.c -0.46 0.00
## 4 Residual        <NA>   <NA>  0.99 0.99
```

``` r
r2mlm(mod2,bargraph = F)
```

```
## $Decompositions
##                        total
## fixed           0.0085753820
## slope variation 0.0005456516
## mean variation  0.0211458825
## sigma2          0.9697330839
## 
## $R2s
##            total
## f   0.0085753820
## v   0.0005456516
## m   0.0211458825
## fv  0.0091210336
## fvm 0.0302669161
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
## mod1    4 1456015 1456059 -728003   1456007                         
## mod2    6 1455821 1455888 -727905   1455809 197.85  2  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

``` r
lvl2_var_cond_lvl1(mod2,lvl1.var = "gndr.c",lvl1.values = c(0.5,-0.5))
```

```
##   lvl1.value lvl2.cond.var lvl2.cond.sd
## 1        0.5    0.01869327    0.1367233
## 2       -0.5    0.02509378    0.1584102
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
##   1455826   1455882   -727908   1455816    492338 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -4.8701 -0.6920 -0.0533  0.6446  5.0775 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  cntry    (Intercept) 0.021352 0.14612 
##  cntry.1  gndr.c      0.002214 0.04705 
##  Residual             0.989875 0.99492 
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##              Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)  0.051008   0.025114 33.959299   2.031   0.0501 .  
## gndr.c       0.188031   0.008682 33.456463  21.659   <2e-16 ***
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
##              Est.    SE     df      t    p   LL    UL
## (Intercept) 0.051 0.025 33.959  2.031 0.05 0.00 0.102
## gndr.c      0.188 0.009 33.456 21.659 0.00 0.17 0.206
```

``` r
getVC(mod2_norecov)
```

```
##        grp        var1 var2 sdcor vcov
## 1    cntry (Intercept) <NA>  0.15 0.02
## 2  cntry.1      gndr.c <NA>  0.05 0.00
## 3 Residual        <NA> <NA>  0.99 0.99
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
## mod2_norecov    5 1455826 1455881 -727908   1455816                        
## mod2            6 1455821 1455888 -727905   1455809 6.8403  1   0.008913 **
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
##   1416796   1416885   -708390   1416780    480356 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -4.8759 -0.6929 -0.0536  0.6461  5.1008 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr  
##  cntry    (Intercept) 0.018551 0.13620        
##           gndr.c      0.002173 0.04661  -0.60 
##  Residual             0.986529 0.99324        
## Number of obs: 480364, groups:  cntry, 33
## 
## Fixed effects:
##                  Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)      0.058614   0.023769 33.001417   2.466    0.019 *  
## gndr.c           0.189614   0.008735 32.195789  21.709   <2e-16 ***
## gei.z.cm        -0.040176   0.024166 33.155548  -1.662    0.106    
## gndr.c:gei.z.cm -0.008456   0.009134 35.793556  -0.926    0.361    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c g.z.cm
## gndr.c      -0.555              
## gei.z.cm    -0.002  0.000       
## gndr.c:g.z.  0.000 -0.042 -0.538
```

``` r
getFE(mod2_GEI,round=3)
```

```
##                   Est.    SE     df      t     p     LL    UL
## (Intercept)      0.059 0.024 33.001  2.466 0.019  0.010 0.107
## gndr.c           0.190 0.009 32.196 21.709 0.000  0.172 0.207
## gei.z.cm        -0.040 0.024 33.156 -1.662 0.106 -0.089 0.009
## gndr.c:gei.z.cm -0.008 0.009 35.794 -0.926 0.361 -0.027 0.010
```

``` r
getVC(mod2_GEI)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.14 0.02
## 2    cntry      gndr.c   <NA>  0.05 0.00
## 3    cntry (Intercept) gndr.c -0.60 0.00
## 4 Residual        <NA>   <NA>  0.99 0.99
```

``` r
r2mlm(mod2_GEI,bargraph = F)
```

```
## $Decompositions
##                        total
## fixed           0.0095226374
## slope variation 0.0005320434
## mean variation  0.0185432322
## sigma2          0.9714020870
## 
## $R2s
##            total
## f   0.0095226374
## v   0.0005320434
## m   0.0185432322
## fv  0.0100546808
## fvm 0.0285979130
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
## Time difference of 30.25013 secs
```

``` r
round(ddsc_mod2_GEI$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.025        0.158         0.99     1.015 0.025   7802.647 0.995   0.995
## 2        0.5         0.019        0.137         0.99     1.009 0.019   6678.029 0.992   0.992
```

``` r
round(ddsc_mod2_GEI$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gei.z.cm gei.z.cm_scaled
## means_y1           0.107 0.134    1.000           1.000    0.934           0.934   -0.266          -0.266
## means_y1_scaled    0.738 0.925    1.000           1.000    0.934           0.934   -0.266          -0.266
## means_y2          -0.071 0.155    0.934           0.934    1.000           1.000   -0.149          -0.149
## means_y2_scaled   -0.489 1.070    0.934           0.934    1.000           1.000   -0.149          -0.149
## gei.z.cm           0.000 1.000   -0.266          -0.266   -0.149          -0.149    1.000           1.000
## gei.z.cm_scaled    0.000 1.000   -0.266          -0.266   -0.149          -0.149    1.000           1.000
## diff_score         0.178 0.056   -0.193          -0.193   -0.530          -0.530   -0.222          -0.222
## diff_score_scaled  1.227 0.388   -0.193          -0.193   -0.530          -0.530   -0.222          -0.222
##                   diff_score diff_score_scaled
## means_y1              -0.193            -0.193
## means_y1_scaled       -0.193            -0.193
## means_y2              -0.530            -0.530
## means_y2_scaled       -0.530            -0.530
## gei.z.cm              -0.222            -0.222
## gei.z.cm_scaled       -0.222            -0.222
## diff_score             1.000             1.000
## diff_score_scaled      1.000             1.000
```

``` r
round(ddsc_mod2_GEI$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                       0.150 0.162 35.794   0.926   0.361   -0.179    0.479
## w_11                         -0.036 0.027 33.226  -1.336   0.191   -0.091    0.019
## w_21                         -0.044 0.022 33.345  -2.014   0.052   -0.089    0.000
## r_xy1                        -0.268 0.201 33.226  -1.336   0.191   -0.676    0.140
## r_xy2                        -0.286 0.142 33.345  -2.014   0.052   -0.575    0.003
## b_11                         -0.249 0.186 33.226  -1.336   0.191   -0.627    0.130
## b_21                         -0.307 0.153 33.345  -2.014   0.052   -0.617    0.003
## main_effect                  -0.040 0.024 33.156  -1.662   0.106   -0.089    0.009
## moderator_effect              0.190 0.009 32.196  21.709   0.000    0.172    0.207
## interaction                  -0.008 0.009 35.794  -0.926   0.361   -0.027    0.010
## q_b11_b21                     0.063    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.020    NA     NA      NA      NA       NA       NA
## cross_over_point             22.423    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.032 0.030 33.406  -1.054   0.299   -0.093    0.029
## interaction_vs_main_bscale   -0.219 0.208 33.406  -1.054   0.299   -0.643    0.204
## interaction_vs_main_rscale   -0.259 0.235 33.368  -1.104   0.277   -0.736    0.218
## dadas                        -0.072 0.054 33.226  -1.336   0.905   -0.181    0.038
## dadas_bscale                 -0.497 0.372 33.226  -1.336   0.905   -1.254    0.260
## dadas_rscale                 -0.536 0.401 33.226  -1.336   0.905   -1.353    0.280
## abs_diff                      0.008 0.009 35.794   0.926   0.180   -0.010    0.027
## abs_sum                       0.080 0.048 33.156   1.662   0.053   -0.018    0.179
## abs_diff_bscale               0.058 0.063 35.794   0.926   0.180   -0.070    0.187
## abs_sum_bscale                0.556 0.334 33.156   1.662   0.053   -0.124    1.236
## abs_diff_rscale               0.018 0.079 34.838   0.229   0.410   -0.143    0.179
## abs_sum_rscale                0.554 0.339 33.154   1.638   0.055   -0.134    1.243
```

``` r
round(ddsc_mod2_GEI$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.003 -0.463  6.840  1.000  0.009
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
##                                    est     se       z pvalue ci.lower ci.upper
## r_xy1_y2                         0.222  0.170   1.310  0.190   -0.110    0.555
## r_xy1                           -0.149  0.172  -0.866  0.387   -0.486    0.188
## r_xy2                           -0.266  0.168  -1.584  0.113   -0.595    0.063
## b_11                            -0.159  0.184  -0.866  0.387   -0.520    0.202
## b_21                            -0.246  0.155  -1.584  0.113   -0.550    0.058
## b_10                            -0.489  0.181  -2.698  0.007   -0.845   -0.134
## b_20                             0.738  0.153   4.828  0.000    0.438    1.037
## res_cov_y1_y2                    0.859  0.218   3.932  0.000    0.431    1.287
## diff_b10_b20                    -1.227  0.065 -18.903  0.000   -1.354   -1.100
## diff_b11_b21                     0.086  0.066   1.310  0.190   -0.043    0.216
## diff_rxy1_rxy2                   0.117  0.060   1.956  0.050    0.000    0.234
## q_b11_b21                        0.090  0.066   1.359  0.174   -0.040    0.220
## q_rxy1_rxy2                      0.122  0.063   1.951  0.051   -0.001    0.245
## cross_over_point                14.215 10.881   1.306  0.191   -7.112   35.541
## sum_b11_b21                     -0.405  0.334  -1.213  0.225   -1.060    0.250
## main_effect                     -0.203  0.167  -1.213  0.225   -0.530    0.125
## interaction_vs_main_effect      -0.116  0.205  -0.567  0.571   -0.518    0.286
## diff_abs_b11_abs_b21            -0.086  0.066  -1.310  0.190   -0.216    0.043
## abs_diff_b11_b21                 0.086  0.066   1.310  0.095   -0.043    0.216
## abs_sum_b11_b21                  0.405  0.334   1.213  0.113   -0.250    1.060
## dadas                           -0.319  0.368  -0.866  0.807   -1.041    0.403
## q_r_equivalence                  0.022  0.063   0.354  0.638       NA       NA
## q_b_equivalence                 -0.010  0.066  -0.149  0.441       NA       NA
## cross_over_point_equivalence    14.215 10.881   1.306  0.904       NA       NA
## cross_over_point_minimal_effect 14.215 10.881   1.306  0.096       NA       NA
```

``` r
round(ddsc_sem_GEI$variance_test,3)
```

```
##             est    se      z pvalue ci.lower ci.upper
## cov_y1y2  0.897 0.229  3.922  0.000    0.449    1.345
## var_y1    1.110 0.273  4.062  0.000    0.574    1.646
## var_y2    0.829 0.204  4.062  0.000    0.429    1.230
## var_diff  0.281 0.138  2.040  0.041    0.011    0.550
## var_ratio 1.338 0.166  8.064  0.000    1.013    1.664
## cor_y1y2  0.934 0.022 42.307  0.000    0.891    0.978
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
## r_xy1y2                       0.222 0.175 31.000   1.269   0.214   -0.135    0.579
## w_11                         -0.023 0.025 33.409  -0.907   0.371   -0.075    0.029
## w_21                         -0.036 0.025 33.409  -1.399   0.171   -0.087    0.016
## r_xy1                        -0.149 0.164 33.409  -0.907   0.371   -0.483    0.185
## r_xy2                        -0.266 0.190 33.409  -1.399   0.171   -0.652    0.121
## b_11                         -0.160 0.176 33.409  -0.907   0.371   -0.518    0.198
## b_21                         -0.246 0.176 33.409  -1.399   0.171   -0.605    0.112
## main_effect                  -0.029 0.025 31.000  -1.175   0.249   -0.080    0.022
## moderator_effect              0.178 0.010 31.000  18.322   0.000    0.158    0.198
## interaction                  -0.013 0.010 31.000  -1.269   0.214   -0.033    0.008
## q_b11_b21                     0.090    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.122    NA     NA      NA      NA       NA       NA
## cross_over_point             14.215    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.017 0.027 40.423  -0.627   0.534   -0.071    0.037
## interaction_vs_main_bscale   -0.117 0.186 40.423  -0.627   0.534   -0.492    0.259
## interaction_vs_main_rscale   -0.091 0.163 42.811  -0.558   0.580   -0.418    0.237
## dadas                        -0.046 0.051 33.409  -0.907   0.815   -0.150    0.057
## dadas_bscale                 -0.320 0.352 33.409  -0.907   0.815   -1.036    0.397
## dadas_rscale                 -0.298 0.328 33.409  -0.907   0.815   -0.966    0.370
## abs_diff                      0.013 0.010 31.000   1.269   0.107   -0.008    0.033
## abs_sum                       0.059 0.050 31.000   1.175   0.124   -0.043    0.161
## abs_diff_bscale               0.087 0.068 31.000   1.269   0.107   -0.053    0.226
## abs_sum_bscale                0.406 0.346 31.000   1.175   0.124   -0.299    1.111
## abs_diff_rscale               0.117 0.073 39.279   1.597   0.059   -0.031    0.264
## abs_sum_rscale                0.415 0.348 31.013   1.193   0.121   -0.294    1.124
```

``` r
# country-time multilevel model


mod2_GEI_cntry_year<-
  lmer(sti.z.wt~gndr.c+gei.z.cm+gndr.c:gei.z.cm+
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
## Formula: sti.z.wt ~ gndr.c + gei.z.cm + gndr.c:gei.z.cm + (gndr.c | cntry)
##    Data: diff_dat_cntry_year
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##    -889.8    -855.6     452.9    -905.8       526 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -3.4287 -0.6342 -0.0169  0.5649  4.1043 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr  
##  cntry    (Intercept) 0.017492 0.13226        
##           gndr.c      0.000681 0.02610  -1.00 
##  Residual             0.008682 0.09318        
## Number of obs: 534, groups:  cntry, 33
## 
## Fixed effects:
##                  Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)      0.056872   0.023463 33.007497   2.424    0.021 *  
## gndr.c           0.187918   0.009454 71.073800  19.877   <2e-16 ***
## gei.z.cm        -0.034532   0.024096 34.447684  -1.433    0.161    
## gndr.c:gei.z.cm -0.011699   0.010584 97.508748  -1.105    0.272    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c g.z.cm
## gndr.c      -0.474              
## gei.z.cm    -0.014  0.003       
## gndr.c:g.z.  0.003 -0.181 -0.432
## optimizer (bobyqa) convergence code: 0 (OK)
## boundary (singular) fit: see help('isSingular')
```

``` r
getFE(mod2_GEI_cntry_year,round=3)
```

```
##                   Est.    SE     df      t     p     LL    UL
## (Intercept)      0.057 0.023 33.007  2.424 0.021  0.009 0.105
## gndr.c           0.188 0.009 71.074 19.877 0.000  0.169 0.207
## gei.z.cm        -0.035 0.024 34.448 -1.433 0.161 -0.083 0.014
## gndr.c:gei.z.cm -0.012 0.011 97.509 -1.105 0.272 -0.033 0.009
```

``` r
getVC(mod2_GEI_cntry_year)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.13 0.02
## 2    cntry      gndr.c   <NA>  0.03 0.00
## 3    cntry (Intercept) gndr.c -1.00 0.00
## 4 Residual        <NA>   <NA>  0.09 0.01
```

``` r
r2mlm(mod2_GEI,bargraph = F)
```

```
## $Decompositions
##                        total
## fixed           0.0095226374
## slope variation 0.0005320434
## mean variation  0.0185432322
## sigma2          0.9714020870
## 
## $R2s
##            total
## f   0.0095226374
## v   0.0005320434
## m   0.0185432322
## fv  0.0100546808
## fvm 0.0285979130
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
## 1       -0.5         0.023        0.152        0.009     0.032 0.720      8.029 0.998   0.954
## 2        0.5         0.017        0.130        0.009     0.026 0.651      8.029 0.998   0.938
```

``` r
round(ddsc_mod2_GEI_cntry_year$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gei.z.cm gei.z.cm_scaled
## means_y1           0.150 0.130    1.000           1.000    0.954           0.954   -0.302          -0.302
## means_y1_scaled    1.048 0.906    1.000           1.000    0.954           0.954   -0.302          -0.302
## means_y2          -0.036 0.156    0.954           0.954    1.000           1.000   -0.203          -0.203
## means_y2_scaled   -0.254 1.086    0.954           0.954    1.000           1.000   -0.203          -0.203
## gei.z.cm           0.000 1.000   -0.302          -0.302   -0.203          -0.203    1.000           1.000
## gei.z.cm_scaled    0.000 1.000   -0.302          -0.302   -0.203          -0.203    1.000           1.000
## diff_score         0.187 0.050   -0.368          -0.368   -0.631          -0.631   -0.151          -0.151
## diff_score_scaled  1.302 0.352   -0.368          -0.368   -0.631          -0.631   -0.151          -0.151
##                   diff_score diff_score_scaled
## means_y1              -0.368            -0.368
## means_y1_scaled       -0.368            -0.368
## means_y2              -0.631            -0.631
## means_y2_scaled       -0.631            -0.631
## gei.z.cm              -0.151            -0.151
## gei.z.cm_scaled       -0.151            -0.151
## diff_score             1.000             1.000
## diff_score_scaled      1.000             1.000
```

``` r
round(ddsc_mod2_GEI_cntry_year$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                       0.232 0.210 97.509   1.105   0.272   -0.184    0.648
## w_11                         -0.029 0.027 34.716  -1.070   0.292   -0.083    0.026
## w_21                         -0.040 0.022 35.178  -1.809   0.079   -0.086    0.005
## r_xy1                        -0.221 0.206 34.716  -1.070   0.292   -0.640    0.198
## r_xy2                        -0.259 0.143 35.178  -1.809   0.079   -0.550    0.032
## b_11                         -0.201 0.188 34.716  -1.070   0.292   -0.582    0.180
## b_21                         -0.283 0.156 35.178  -1.809   0.079   -0.600    0.035
## main_effect                  -0.035 0.024 34.448  -1.433   0.161   -0.083    0.014
## moderator_effect              0.188 0.009 71.074  19.877   0.000    0.169    0.207
## interaction                  -0.012 0.011 97.509  -1.105   0.272   -0.033    0.009
## q_b11_b21                     0.087    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.041    NA     NA      NA      NA       NA       NA
## cross_over_point             16.062    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.023 0.030 35.608  -0.756   0.455   -0.084    0.038
## interaction_vs_main_bscale   -0.160 0.211 35.608  -0.756   0.455   -0.589    0.269
## interaction_vs_main_rscale   -0.201 0.245 35.357  -0.823   0.416   -0.698    0.295
## dadas                        -0.057 0.054 34.716  -1.070   0.854   -0.166    0.052
## dadas_bscale                 -0.401 0.375 34.716  -1.070   0.854   -1.163    0.361
## dadas_rscale                 -0.441 0.413 34.716  -1.070   0.854   -1.279    0.396
## abs_diff                      0.012 0.011 97.509   1.105   0.136   -0.009    0.033
## abs_sum                       0.069 0.048 34.448   1.433   0.080   -0.029    0.167
## abs_diff_bscale               0.082 0.074 97.509   1.105   0.136   -0.065    0.229
## abs_sum_bscale                0.483 0.337 34.448   1.433   0.080   -0.202    1.168
## abs_diff_rscale               0.039 0.092 51.276   0.418   0.339   -0.147    0.224
## abs_sum_rscale                0.480 0.343 34.440   1.399   0.085   -0.217    1.177
```

``` r
round(ddsc_mod2_GEI_cntry_year$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.003 -1.000  6.458  1.000  0.011
```

### Comparisons across models


``` r
# full multilevel model (individuals as entries)
round(ddsc_mod2_GEI$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df t.ratio p.value ci.lower ci.upper
## r_xy1             -0.2682 0.2007 33.2264 -1.3363  0.1905  -0.6764   0.1400
## r_xy2             -0.2863 0.1422 33.3448 -2.0141  0.0521  -0.5754   0.0028
## b_11              -0.2487 0.1861 33.2264 -1.3363  0.1905  -0.6271   0.1298
## b_21              -0.3072 0.1525 33.3448 -2.0141  0.0521  -0.6173   0.0030
## main_effect       -0.0402 0.0242 33.1555 -1.6625  0.1058  -0.0893   0.0090
## moderator_effect   0.1896 0.0087 32.1958 21.7086  0.0000   0.1718   0.2074
## interaction       -0.0085 0.0091 35.7936 -0.9258  0.3608  -0.0270   0.0101
## q_b11_b21          0.0634     NA      NA      NA      NA       NA       NA
```

``` r
# country-level structural path model (country-gender means as entries)
round(ddsc_sem_GEI$results[c("r_xy1","r_xy2","b_11","b_21","q_b11_b21","diff_b11_b21"),],4)
```

```
##                  est     se       z pvalue ci.lower ci.upper
## r_xy1        -0.1490 0.1721 -0.8658 0.3866  -0.4864   0.1883
## r_xy2        -0.2658 0.1678 -1.5836 0.1133  -0.5947   0.0632
## b_11         -0.1594 0.1842 -0.8658 0.3866  -0.5204   0.2015
## b_21         -0.2458 0.1552 -1.5836 0.1133  -0.5500   0.0584
## q_b11_b21     0.0901 0.0663  1.3589 0.1742  -0.0399   0.2200
## diff_b11_b21  0.0863 0.0659  1.3095 0.1904  -0.0429   0.2155
```

``` r
# multilevel model (country-gender means as entries)
round(ddsc_mod2_GEI_ri$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df t.ratio p.value ci.lower ci.upper
## r_xy1             -0.1490 0.1642 33.4092 -0.9075  0.3707  -0.4830   0.1849
## r_xy2             -0.2658 0.1900 33.4092 -1.3987  0.1711  -0.6521   0.1206
## b_11              -0.1599 0.1762 33.4092 -0.9075  0.3707  -0.5181   0.1984
## b_21              -0.2464 0.1762 33.4092 -1.3987  0.1711  -0.6047   0.1118
## main_effect       -0.0294 0.0250 31.0000 -1.1753  0.2488  -0.0803   0.0216
## moderator_effect   0.1779 0.0097 31.0000 18.3215  0.0000   0.1581   0.1977
## interaction       -0.0125 0.0099 31.0000 -1.2692  0.2138  -0.0326   0.0076
## q_b11_b21          0.0903     NA      NA      NA      NA       NA       NA
```

``` r
# multilevel model (country-year-gender means as entries)
round(ddsc_mod2_GEI_cntry_year$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df t.ratio p.value ci.lower ci.upper
## r_xy1             -0.2207 0.2063 34.7165 -1.0699  0.2921  -0.6395   0.1982
## r_xy2             -0.2592 0.1433 35.1777 -1.8086  0.0791  -0.5500   0.0317
## b_11              -0.2007 0.1876 34.7165 -1.0699  0.2921  -0.5817   0.1803
## b_21              -0.2826 0.1563 35.1777 -1.8086  0.0791  -0.5998   0.0346
## main_effect       -0.0345 0.0241 34.4477 -1.4331  0.1608  -0.0835   0.0144
## moderator_effect   0.1879 0.0095 71.0738 19.8769  0.0000   0.1691   0.2068
## interaction       -0.0117 0.0106 97.5087 -1.1054  0.2717  -0.0327   0.0093
## q_b11_b21          0.0870     NA      NA      NA      NA       NA       NA
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
## Time difference of 1.354424 hours
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
##                     Estimate          SE        2.5%        97.5%
## X.Intercept.     0.057952895 0.023483444  0.01001366 0.1051400238
## gndr.c           0.190064425 0.008906213  0.17299075 0.2069660694
## gei.z.cm        -0.039647715 0.024219970 -0.08616109 0.0113852963
## gndr.c.gei.z.cm -0.008857001 0.009348737 -0.02722690 0.0085019286
## w11             -0.035219215 0.026794105 -0.08738045 0.0191456627
## w21             -0.044076216 0.022338069 -0.08721230 0.0009771841
## b11             -0.243622848 0.185343600 -0.60443917 0.1324368218
## b21             -0.304889627 0.154519743 -0.60327604 0.0067595023
## r_xy1           -0.262737418 0.199885598 -0.65186327 0.1428277718
## r_xy2           -0.284212696 0.144040560 -0.56236321 0.0063010880
## q_b              0.065512688 0.071228496 -0.07602965 0.2014870963
## q                0.016923166 0.093362348 -0.18216512 0.1776011984
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
## [1] 0.06551269
## 
## $se
## [1] 0.0712285
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
## [1] 2.323686
## 
## $p_low
## [1] 0.01007115
## 
## $z_high
## [1] -0.4841786
## 
## $p_high
## [1] 0.3141296
## 
## $ci_level
## [1] 0.9
## 
## $ci_lower
## [1] -0.05164776
## 
## $ci_upper
## [1] 0.1826731
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
## [1] 0.01692317
## 
## $se
## [1] 0.09336235
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
## [1] 1.252359
## 
## $p_low
## [1] 0.1052196
## 
## $z_high
## [1] -0.8898323
## 
## $p_high
## [1] 0.186778
## 
## $ci_level
## [1] 0.9
## 
## $ci_lower
## [1] -0.1366442
## 
## $ci_upper
## [1] 0.1704906
## 
## $equivalent
## [1] FALSE
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
  xlab("Gender Equality Index (Average 2002-2023)")+
  #ylim(c(min.y.pred,max.y.pred))+
  #xlim(c(0.60,1.00))+
  ylab("Value stimulation (Average 2002-2023)")+
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
  xlab("Gender Equality Index (Average 2002-2023)")+
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

![](Analysis_stimulation_files/figure-html/unnamed-chunk-23-1.png)<!-- -->

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
## 1076597.5 1076684.0 -538290.8 1076581.5    363844 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -4.9684 -0.6919 -0.0340  0.6465  5.1233 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr  
##  cntry    (Intercept) 0.025262 0.1589         
##           gndr.c      0.002247 0.0474   -0.43 
##  Residual             0.990151 0.9951         
## Number of obs: 363852, groups:  cntry, 34
## 
## Fixed effects:
##                   Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)       0.070947   0.027331 33.914295   2.596   0.0138 *  
## gndr.c            0.184672   0.008944 31.245991  20.648   <2e-16 ***
## gggi.z.cm        -0.029227   0.027772 34.057785  -1.052   0.3000    
## gndr.c:gggi.z.cm -0.007023   0.009380 34.883359  -0.749   0.4591    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c ggg.z.
## gndr.c      -0.386              
## gggi.z.cm   -0.001  0.000       
## gndr.c:gg..  0.000 -0.021 -0.374
```

``` r
getFE(mod2_GGGI,round=3)
```

```
##                    Est.    SE     df      t     p     LL    UL
## (Intercept)       0.071 0.027 33.914  2.596 0.014  0.015 0.126
## gndr.c            0.185 0.009 31.246 20.648 0.000  0.166 0.203
## gggi.z.cm        -0.029 0.028 34.058 -1.052 0.300 -0.086 0.027
## gndr.c:gggi.z.cm -0.007 0.009 34.883 -0.749 0.459 -0.026 0.012
```

``` r
getVC(mod2_GGGI)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.16 0.03
## 2    cntry      gndr.c   <NA>  0.05 0.00
## 3    cntry (Intercept) gndr.c -0.43 0.00
## 4 Residual        <NA>   <NA>  1.00 0.99
```

``` r
r2mlm(mod2_GGGI,bargraph = F)
```

```
## $Decompositions
##                        total
## fixed           0.0086729881
## slope variation 0.0005444066
## mean variation  0.0248973790
## sigma2          0.9658852264
## 
## $R2s
##            total
## f   0.0086729881
## v   0.0005444066
## m   0.0248973790
## fv  0.0092173946
## fvm 0.0341147736
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
## Time difference of 31.44436 secs
```

``` r
round(ddsc_mod2_GGGI$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.025        0.158         0.99     1.015 0.025   7802.647 0.995   0.995
## 2        0.5         0.019        0.137         0.99     1.009 0.019   6678.029 0.992   0.992
```

``` r
round(ddsc_mod2_GGGI$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gggi.z.cm
## means_y1           0.120 0.157    1.000           1.000    0.945           0.945    -0.233
## means_y1_scaled    0.722 0.949    1.000           1.000    0.945           0.945    -0.233
## means_y2          -0.053 0.174    0.945           0.945    1.000           1.000    -0.149
## means_y2_scaled   -0.318 1.048    0.945           0.945    1.000           1.000    -0.149
## gggi.z.cm          0.000 1.000   -0.233          -0.233   -0.149          -0.149     1.000
## gggi.z.cm_scaled   0.000 1.000   -0.233          -0.233   -0.149          -0.149     1.000
## diff_score         0.172 0.057   -0.117          -0.117   -0.437          -0.437    -0.189
## diff_score_scaled  1.040 0.347   -0.117          -0.117   -0.437          -0.437    -0.189
##                   gggi.z.cm_scaled diff_score diff_score_scaled
## means_y1                    -0.233     -0.117            -0.117
## means_y1_scaled             -0.233     -0.117            -0.117
## means_y2                    -0.149     -0.437            -0.437
## means_y2_scaled             -0.149     -0.437            -0.437
## gggi.z.cm                    1.000     -0.189            -0.189
## gggi.z.cm_scaled             1.000     -0.189            -0.189
## diff_score                  -0.189      1.000             1.000
## diff_score_scaled           -0.189      1.000             1.000
```

``` r
round(ddsc_mod2_GGGI$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                       0.122 0.163 34.883   0.749   0.459   -0.209    0.454
## w_11                         -0.026 0.030 34.082  -0.862   0.395   -0.086    0.035
## w_21                         -0.033 0.026 34.104  -1.241   0.223   -0.086    0.021
## r_xy1                        -0.163 0.190 34.082  -0.862   0.395   -0.549    0.222
## r_xy2                        -0.188 0.152 34.104  -1.241   0.223   -0.497    0.120
## b_11                         -0.155 0.180 34.082  -0.862   0.395   -0.522    0.211
## b_21                         -0.198 0.159 34.104  -1.241   0.223   -0.522    0.126
## main_effect                  -0.029 0.028 34.058  -1.052   0.300   -0.086    0.027
## moderator_effect              0.185 0.009 31.246  20.648   0.000    0.166    0.203
## interaction                  -0.007 0.009 34.883  -0.749   0.459   -0.026    0.012
## q_b11_b21                     0.044    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.026    NA     NA      NA      NA       NA       NA
## cross_over_point             26.296    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.022 0.032 34.118  -0.684   0.499   -0.088    0.044
## interaction_vs_main_bscale   -0.134 0.196 34.118  -0.684   0.499   -0.533    0.264
## interaction_vs_main_rscale   -0.151 0.214 34.113  -0.707   0.484   -0.585    0.283
## dadas                        -0.051 0.060 34.082  -0.862   0.803   -0.173    0.070
## dadas_bscale                 -0.311 0.361 34.082  -0.862   0.803   -1.044    0.422
## dadas_rscale                 -0.327 0.379 34.082  -0.862   0.803   -1.098    0.444
## abs_diff                      0.007 0.009 34.883   0.749   0.230   -0.012    0.026
## abs_sum                       0.058 0.056 34.058   1.052   0.150   -0.054    0.171
## abs_diff_bscale               0.042 0.057 34.883   0.749   0.230   -0.073    0.157
## abs_sum_bscale                0.353 0.336 34.058   1.052   0.150   -0.329    1.035
## abs_diff_rscale               0.025 0.065 34.426   0.386   0.351   -0.107    0.157
## abs_sum_rscale                0.352 0.337 34.058   1.043   0.152   -0.334    1.038
```

``` r
round(ddsc_mod2_GGGI$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.003 -0.463  6.840  1.000  0.009
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
##                                    est     se       z pvalue ci.lower ci.upper
## r_xy1_y2                         0.189  0.168   1.121  0.262   -0.141    0.519
## r_xy1                           -0.149  0.170  -0.877  0.380   -0.481    0.184
## r_xy2                           -0.233  0.167  -1.398  0.162   -0.560    0.094
## b_11                            -0.156  0.178  -0.877  0.380   -0.504    0.192
## b_21                            -0.221  0.158  -1.398  0.162   -0.532    0.089
## b_10                            -0.318  0.175  -1.818  0.069   -0.662    0.025
## b_20                             0.722  0.156   4.626  0.000    0.416    1.027
## res_cov_y1_y2                    0.879  0.219   4.008  0.000    0.449    1.309
## diff_b10_b20                    -1.040  0.057 -18.086  0.000   -1.153   -0.927
## diff_b11_b21                     0.065  0.058   1.121  0.262   -0.049    0.180
## diff_rxy1_rxy2                   0.084  0.055   1.528  0.127   -0.024    0.193
## q_b11_b21                        0.068  0.059   1.146  0.252   -0.048    0.184
## q_rxy1_rxy2                      0.088  0.057   1.526  0.127   -0.025    0.200
## cross_over_point                15.889 14.196   1.119  0.263  -11.935   43.713
## sum_b11_b21                     -0.377  0.332  -1.138  0.255   -1.027    0.273
## main_effect                     -0.189  0.166  -1.138  0.255   -0.514    0.136
## interaction_vs_main_effect      -0.123  0.193  -0.637  0.524   -0.502    0.256
## diff_abs_b11_abs_b21            -0.065  0.058  -1.121  0.262   -0.180    0.049
## abs_diff_b11_b21                 0.065  0.058   1.121  0.131   -0.049    0.180
## abs_sum_b11_b21                  0.377  0.332   1.138  0.128   -0.273    1.027
## dadas                           -0.312  0.356  -0.877  0.810   -1.009    0.385
## q_r_equivalence                 -0.012  0.057  -0.215  0.415       NA       NA
## q_b_equivalence                 -0.032  0.059  -0.542  0.294       NA       NA
## cross_over_point_equivalence    15.889 14.196   1.119  0.868       NA       NA
## cross_over_point_minimal_effect 15.889 14.196   1.119  0.132       NA       NA
```

``` r
round(ddsc_sem_GGGI$variance_test,3)
```

```
##             est    se      z pvalue ci.lower ci.upper
## cov_y1y2  0.912 0.228  4.004  0.000    0.466    1.359
## var_y1    1.066 0.259  4.123  0.000    0.559    1.573
## var_y2    0.875 0.212  4.123  0.000    0.459    1.291
## var_diff  0.191 0.118  1.619  0.106   -0.040    0.423
## var_ratio 1.219 0.137  8.879  0.000    0.950    1.488
## cor_y1y2  0.945 0.018 51.082  0.000    0.908    0.981
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
## r_xy1y2                       0.189 0.174 32.000   1.088   0.285   -0.165    0.542
## w_11                         -0.026 0.029 33.981  -0.899   0.375   -0.084    0.033
## w_21                         -0.037 0.029 33.981  -1.276   0.211   -0.095    0.022
## r_xy1                        -0.149 0.166 33.981  -0.899   0.375   -0.485    0.188
## r_xy2                        -0.233 0.183 33.981  -1.276   0.211   -0.605    0.138
## b_11                         -0.156 0.174 33.981  -0.899   0.375   -0.509    0.197
## b_21                         -0.222 0.174 33.981  -1.276   0.211   -0.575    0.131
## main_effect                  -0.031 0.028 32.000  -1.104   0.278   -0.089    0.026
## moderator_effect              0.172 0.010 32.000  17.546   0.000    0.152    0.192
## interaction                  -0.011 0.010 32.000  -1.088   0.285   -0.031    0.009
## q_b11_b21                     0.068    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.088    NA     NA      NA      NA       NA       NA
## cross_over_point             15.889    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.020 0.030 39.813  -0.680   0.500   -0.081    0.040
## interaction_vs_main_bscale   -0.123 0.181 39.813  -0.680   0.500   -0.490    0.243
## interaction_vs_main_rscale   -0.107 0.165 41.100  -0.644   0.523   -0.441    0.227
## dadas                        -0.052 0.058 33.981  -0.899   0.812   -0.169    0.065
## dadas_bscale                 -0.312 0.347 33.981  -0.899   0.812   -1.018    0.394
## dadas_rscale                 -0.297 0.331 33.981  -0.899   0.812   -0.970    0.375
## abs_diff                      0.011 0.010 32.000   1.088   0.142   -0.009    0.031
## abs_sum                       0.063 0.057 32.000   1.104   0.139   -0.053    0.178
## abs_diff_bscale               0.066 0.060 32.000   1.088   0.142   -0.057    0.188
## abs_sum_bscale                0.378 0.342 32.000   1.104   0.139   -0.319    1.075
## abs_diff_rscale               0.084 0.063 37.015   1.346   0.093   -0.043    0.211
## abs_sum_rscale                0.382 0.343 32.005   1.113   0.137   -0.317    1.081
```

``` r
# country-time multilevel model


mod2_GGGI_cntry_year<-
  lmer(sti.z.wt~gndr.c+gggi.z.cm+gndr.c:gggi.z.cm+
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
## Formula: sti.z.wt ~ gndr.c + gggi.z.cm + gndr.c:gggi.z.cm + (gndr.c |      cntry)
##    Data: diff_dat_cntry_year
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##    -721.6    -689.7     368.8    -737.6       392 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -3.3535 -0.5638 -0.0065  0.5819  4.1974 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr  
##  cntry    (Intercept) 0.023765 0.15416        
##           gndr.c      0.000324 0.01800  -1.00 
##  Residual             0.006796 0.08244        
## Number of obs: 400, groups:  cntry, 34
## 
## Fixed effects:
##                    Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)        0.069656   0.026867  33.362743   2.593    0.014 *  
## gndr.c             0.184779   0.008889 100.345372  20.788   <2e-16 ***
## gggi.z.cm         -0.028287   0.027425  34.061711  -1.031    0.310    
## gndr.c:gggi.z.cm  -0.012912   0.009491 115.142980  -1.361    0.176    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c ggg.z.
## gndr.c      -0.344              
## gggi.z.cm   -0.009  0.002       
## gndr.c:gg..  0.002 -0.128 -0.328
## optimizer (bobyqa) convergence code: 0 (OK)
## boundary (singular) fit: see help('isSingular')
```

``` r
getFE(mod2_GGGI_cntry_year,round=3)
```

```
##                    Est.    SE      df      t     p     LL    UL
## (Intercept)       0.070 0.027  33.363  2.593 0.014  0.015 0.124
## gndr.c            0.185 0.009 100.345 20.788 0.000  0.167 0.202
## gggi.z.cm        -0.028 0.027  34.062 -1.031 0.310 -0.084 0.027
## gndr.c:gggi.z.cm -0.013 0.009 115.143 -1.361 0.176 -0.032 0.006
```

``` r
getVC(mod2_GGGI_cntry_year)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.15 0.02
## 2    cntry      gndr.c   <NA>  0.02 0.00
## 3    cntry (Intercept) gndr.c -1.00 0.00
## 4 Residual        <NA>   <NA>  0.08 0.01
```

``` r
r2mlm(mod2_GGGI,bargraph = F)
```

```
## $Decompositions
##                        total
## fixed           0.0086729881
## slope variation 0.0005444066
## mean variation  0.0248973790
## sigma2          0.9658852264
## 
## $R2s
##            total
## f   0.0086729881
## v   0.0005444066
## m   0.0248973790
## fv  0.0092173946
## fvm 0.0341147736
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
## 1       -0.5         0.023        0.152        0.009     0.032 0.720      8.029 0.998   0.954
## 2        0.5         0.017        0.130        0.009     0.026 0.651      8.029 0.998   0.938
```

``` r
round(ddsc_mod2_GGGI_cntry_year$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gggi.z.cm
## means_y1           0.162 0.154    1.000           1.000    0.955           0.955    -0.216
## means_y1_scaled    0.976 0.933    1.000           1.000    0.955           0.955    -0.216
## means_y2          -0.020 0.176    0.955           0.955    1.000           1.000    -0.157
## means_y2_scaled   -0.121 1.063    0.955           0.955    1.000           1.000    -0.157
## gggi.z.cm          0.000 1.000   -0.216          -0.216   -0.157          -0.157     1.000
## gggi.z.cm_scaled   0.000 1.000   -0.216          -0.216   -0.157          -0.157     1.000
## diff_score         0.182 0.054   -0.253          -0.253   -0.528          -0.528    -0.107
## diff_score_scaled  1.097 0.326   -0.253          -0.253   -0.528          -0.528    -0.107
##                   gggi.z.cm_scaled diff_score diff_score_scaled
## means_y1                    -0.216     -0.253            -0.253
## means_y1_scaled             -0.216     -0.253            -0.253
## means_y2                    -0.157     -0.528            -0.528
## means_y2_scaled             -0.157     -0.528            -0.528
## gggi.z.cm                    1.000     -0.107            -0.107
## gggi.z.cm_scaled             1.000     -0.107            -0.107
## diff_score                  -0.107      1.000             1.000
## diff_score_scaled           -0.107      1.000             1.000
```

``` r
round(ddsc_mod2_GGGI_cntry_year$results,3)
```

```
##                            estimate    SE      df t.ratio p.value ci.lower ci.upper
## r_xy1y2                       0.239 0.176 115.143   1.361   0.176   -0.109    0.588
## w_11                         -0.022 0.029  33.976  -0.744   0.462   -0.081    0.038
## w_21                         -0.035 0.026  34.061  -1.323   0.194   -0.088    0.019
## r_xy1                        -0.141 0.190  33.976  -0.744   0.462   -0.527    0.245
## r_xy2                        -0.197 0.149  34.061  -1.323   0.194   -0.501    0.106
## b_11                         -0.132 0.178  33.976  -0.744   0.462   -0.493    0.229
## b_21                         -0.210 0.159  34.061  -1.323   0.194   -0.533    0.113
## main_effect                  -0.028 0.027  34.062  -1.031   0.310   -0.084    0.027
## moderator_effect              0.185 0.009 100.345  20.788   0.000    0.167    0.202
## interaction                  -0.013 0.009 115.143  -1.361   0.176   -0.032    0.006
## q_b11_b21                     0.081    NA      NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.058    NA      NA      NA      NA       NA       NA
## cross_over_point             14.310    NA      NA      NA      NA       NA       NA
## interaction_vs_main          -0.015 0.032  33.973  -0.483   0.632   -0.080    0.049
## interaction_vs_main_bscale   -0.093 0.193  33.973  -0.483   0.632   -0.485    0.298
## interaction_vs_main_rscale   -0.113 0.215  33.956  -0.526   0.602   -0.551    0.325
## dadas                        -0.044 0.059  33.976  -0.744   0.769   -0.163    0.076
## dadas_bscale                 -0.264 0.355  33.976  -0.744   0.769   -0.986    0.457
## dadas_rscale                 -0.283 0.380  33.976  -0.744   0.769   -1.055    0.489
## abs_diff                      0.013 0.009 115.143   1.361   0.088   -0.006    0.032
## abs_sum                       0.057 0.055  34.062   1.031   0.155   -0.055    0.168
## abs_diff_bscale               0.078 0.057 115.143   1.361   0.088   -0.036    0.192
## abs_sum_bscale                0.342 0.332  34.062   1.031   0.155   -0.332    1.017
## abs_diff_rscale               0.056 0.068  51.893   0.824   0.207   -0.080    0.193
## abs_sum_rscale                0.339 0.335  34.059   1.012   0.159   -0.341    1.019
```

``` r
round(ddsc_mod2_GGGI_cntry_year$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.003 -1.000  6.458  1.000  0.011
```

### Comparisons across models


``` r
# full multilevel model (individuals as entries)
round(ddsc_mod2_GGGI$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df t.ratio p.value ci.lower ci.upper
## r_xy1             -0.1634 0.1897 34.0815 -0.8617  0.3949  -0.5489   0.2220
## r_xy2             -0.1885 0.1519 34.1042 -1.2411  0.2230  -0.4971   0.1201
## b_11              -0.1554 0.1803 34.0815 -0.8617  0.3949  -0.5218   0.2110
## b_21              -0.1978 0.1594 34.1042 -1.2411  0.2230  -0.5217   0.1261
## main_effect       -0.0292 0.0278 34.0578 -1.0524  0.3000  -0.0857   0.0272
## moderator_effect   0.1847 0.0089 31.2460 20.6479  0.0000   0.1664   0.2029
## interaction       -0.0070 0.0094 34.8834 -0.7487  0.4591  -0.0261   0.0120
## q_b11_b21          0.0438     NA      NA      NA      NA       NA       NA
```

``` r
# country-level structural path model (country-gender means as entries)
round(ddsc_sem_GGGI$results[c("r_xy1","r_xy2","b_11","b_21","q_b11_b21","diff_b11_b21"),],4)
```

```
##                  est     se       z pvalue ci.lower ci.upper
## r_xy1        -0.1487 0.1696 -0.8771 0.3804  -0.4811   0.1836
## r_xy2        -0.2332 0.1668 -1.3980 0.1621  -0.5600   0.0937
## b_11         -0.1559 0.1778 -0.8771 0.3804  -0.5043   0.1925
## b_21         -0.2214 0.1583 -1.3980 0.1621  -0.5317   0.0890
## q_b11_b21     0.0679 0.0593  1.1458 0.2519  -0.0482   0.1840
## diff_b11_b21  0.0655 0.0584  1.1214 0.2621  -0.0489   0.1798
```

``` r
# multilevel model (country-gender means as entries)
round(ddsc_mod2_GGGI_ri$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df t.ratio p.value ci.lower ci.upper
## r_xy1             -0.1487 0.1655 33.9813 -0.8986  0.3752  -0.4852   0.1877
## r_xy2             -0.2332 0.1828 33.9813 -1.2758  0.2107  -0.6046   0.1383
## b_11              -0.1561 0.1737 33.9813 -0.8986  0.3752  -0.5091   0.1970
## b_21              -0.2216 0.1737 33.9813 -1.2758  0.2107  -0.5747   0.1314
## main_effect       -0.0313 0.0283 32.0000 -1.1039  0.2779  -0.0889   0.0264
## moderator_effect   0.1723 0.0098 32.0000 17.5465  0.0000   0.1523   0.1923
## interaction       -0.0108 0.0100 32.0000 -1.0879  0.2848  -0.0312   0.0095
## q_b11_b21          0.0680     NA      NA      NA      NA       NA       NA
```

``` r
# multilevel model (country-year-gender means as entries)
round(ddsc_mod2_GGGI_cntry_year$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE       df t.ratio p.value ci.lower ci.upper
## r_xy1             -0.1414 0.1899  33.9761 -0.7444  0.4618  -0.5273   0.2446
## r_xy2             -0.1974 0.1492  34.0612 -1.3235  0.1945  -0.5005   0.1057
## b_11              -0.1321 0.1775  33.9761 -0.7444  0.4618  -0.4929   0.2286
## b_21              -0.2103 0.1589  34.0612 -1.3235  0.1945  -0.5332   0.1126
## main_effect       -0.0283 0.0274  34.0617 -1.0315  0.3096  -0.0840   0.0274
## moderator_effect   0.1848 0.0089 100.3454 20.7879  0.0000   0.1671   0.2024
## interaction       -0.0129 0.0095 115.1430 -1.3605  0.1763  -0.0317   0.0059
## q_b11_b21          0.0806     NA       NA      NA      NA       NA       NA
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
## Time difference of 1.020688 hours
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
##                      Estimate          SE        2.5%      97.5%
## X.Intercept.      0.070360784 0.027610952  0.01943797 0.12529502
## gndr.c            0.185213063 0.009596632  0.16668884 0.20274258
## gggi.z.cm        -0.028505653 0.027044580 -0.08357812 0.02383294
## gndr.c.gggi.z.cm -0.006996375 0.009753662 -0.02646536 0.01304855
## w11              -0.025007466 0.029071860 -0.08369161 0.03252338
## w21              -0.032003840 0.025791714 -0.08453287 0.01871421
## b11              -0.151089369 0.175645506 -0.50564548 0.19649881
## b21              -0.193359860 0.155827616 -0.51072821 0.11306696
## r_xy1            -0.158946011 0.184779067 -0.53193903 0.20671674
## r_xy2            -0.184252335 0.148487914 -0.48667218 0.10774135
## q_b               0.043746223 0.062685586 -0.08447818 0.16746901
## q                 0.024333260 0.071751121 -0.13016032 0.16897492
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
## [1] 0.04374622
## 
## $se
## [1] 0.06268559
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
## [1] 2.29313
## 
## $p_low
## [1] 0.01092025
## 
## $z_high
## [1] -0.8973957
## 
## $p_high
## [1] 0.1847539
## 
## $ci_level
## [1] 0.9
## 
## $ci_lower
## [1] -0.05936239
## 
## $ci_upper
## [1] 0.1468548
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
## [1] 0.02433326
## 
## $se
## [1] 0.07175112
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
## [1] 1.732841
## 
## $p_low
## [1] 0.041562
## 
## $z_high
## [1] -1.054572
## 
## $p_high
## [1] 0.1458105
## 
## $ci_level
## [1] 0.9
## 
## $ci_lower
## [1] -0.09368683
## 
## $ci_upper
## [1] 0.1423534
## 
## $equivalent
## [1] FALSE
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
  xlab("Global Gender Gap Index (Average 2002-2023)")+
  #ylim(c(min.y.pred,max.y.pred))+
  #xlim(c(0.60,1.00))+
  ylab("Value stimulation (Average 2002-2023)")+
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
  xlab("Global Gender Gap Index (Average 2002-2023)")+
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

![](Analysis_stimulation_files/figure-html/unnamed-chunk-29-1.png)<!-- -->

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
## 1455817.8 1455906.6 -727900.9 1455801.8    492335 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -4.8708 -0.6918 -0.0538  0.6445  5.0839 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr  
##  cntry    (Intercept) 0.018749 0.1369         
##           gndr.c      0.001857 0.0431   -0.37 
##  Residual             0.989874 0.9949         
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##                  Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)      0.051077   0.023541 33.956846   2.170   0.0371 *  
## gndr.c           0.187515   0.008045 34.032512  23.308   <2e-16 ***
## gdi.z.cm        -0.051677   0.023929 34.147116  -2.160   0.0379 *  
## gndr.c:gdi.z.cm  0.020966   0.008507 39.346718   2.465   0.0182 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c gd.z.c
## gndr.c      -0.338              
## gdi.z.cm    -0.001  0.001       
## gndr.c:gd..  0.001 -0.019 -0.323
```

``` r
getFE(mod2_GDI,round=3)
```

```
##                   Est.    SE     df      t     p     LL     UL
## (Intercept)      0.051 0.024 33.957  2.170 0.037  0.003  0.099
## gndr.c           0.188 0.008 34.033 23.308 0.000  0.171  0.204
## gdi.z.cm        -0.052 0.024 34.147 -2.160 0.038 -0.100 -0.003
## gndr.c:gdi.z.cm  0.021 0.009 39.347  2.465 0.018  0.004  0.038
```

``` r
getVC(mod2_GDI)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.14 0.02
## 2    cntry      gndr.c   <NA>  0.04 0.00
## 3    cntry (Intercept) gndr.c -0.37 0.00
## 4 Residual        <NA>   <NA>  0.99 0.99
```

``` r
r2mlm(mod2_GDI,bargraph = F)
```

```
## $Decompositions
##                        total
## fixed           0.0109047928
## slope variation 0.0004523474
## mean variation  0.0185432765
## sigma2          0.9700995833
## 
## $R2s
##            total
## f   0.0109047928
## v   0.0004523474
## m   0.0185432765
## fv  0.0113571402
## fvm 0.0299004167
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
## Time difference of 32.45997 secs
```

``` r
round(ddsc_mod2_GDI$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.025        0.158         0.99     1.015 0.025   7802.647 0.995   0.995
## 2        0.5         0.019        0.137         0.99     1.009 0.019   6678.029 0.992   0.992
```

``` r
round(ddsc_mod2_GDI$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gdi.z.cm gdi.z.cm_scaled
## means_y1           0.099 0.140    1.000           1.000    0.939           0.939   -0.303          -0.303
## means_y1_scaled    0.662 0.934    1.000           1.000    0.939           0.939   -0.303          -0.303
## means_y2          -0.078 0.159    0.939           0.939    1.000           1.000   -0.448          -0.448
## means_y2_scaled   -0.525 1.062    0.939           0.939    1.000           1.000   -0.448          -0.448
## gdi.z.cm           0.000 1.000   -0.303          -0.303   -0.448          -0.448    1.000           1.000
## gdi.z.cm_scaled    0.000 1.000   -0.303          -0.303   -0.448          -0.448    1.000           1.000
## diff_score         0.178 0.055   -0.173          -0.173   -0.501          -0.501    0.521           0.521
## diff_score_scaled  1.187 0.371   -0.173          -0.173   -0.501          -0.501    0.521           0.521
##                   diff_score diff_score_scaled
## means_y1              -0.173            -0.173
## means_y1_scaled       -0.173            -0.173
## means_y2              -0.501            -0.501
## means_y2_scaled       -0.501            -0.501
## gdi.z.cm               0.521             0.521
## gdi.z.cm_scaled        0.521             0.521
## diff_score             1.000             1.000
## diff_score_scaled      1.000             1.000
```

``` r
round(ddsc_mod2_GDI$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                      -0.378 0.153 39.347  -2.465   0.018   -0.688   -0.068
## w_11                         -0.062 0.026 34.236  -2.426   0.021   -0.114   -0.010
## w_21                         -0.041 0.023 34.304  -1.798   0.081   -0.088    0.005
## r_xy1                        -0.445 0.183 34.236  -2.426   0.021   -0.818   -0.072
## r_xy2                        -0.259 0.144 34.304  -1.798   0.081   -0.552    0.034
## b_11                         -0.416 0.172 34.236  -2.426   0.021   -0.765   -0.068
## b_21                         -0.276 0.153 34.304  -1.798   0.081   -0.588    0.036
## main_effect                  -0.052 0.024 34.147  -2.160   0.038   -0.100   -0.003
## moderator_effect              0.188 0.008 34.033  23.308   0.000    0.171    0.204
## interaction                   0.021 0.009 39.347   2.465   0.018    0.004    0.038
## q_b11_b21                    -0.160    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                  -0.213    NA     NA      NA      NA       NA       NA
## cross_over_point             -8.944    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.031 0.023 34.763  -1.355   0.184   -0.077    0.015
## interaction_vs_main_bscale   -0.206 0.152 34.763  -1.355   0.184   -0.514    0.102
## interaction_vs_main_rscale   -0.166 0.133 34.906  -1.250   0.220   -0.436    0.104
## dadas                        -0.082 0.046 34.304  -1.798   0.960   -0.175    0.011
## dadas_bscale                 -0.552 0.307 34.304  -1.798   0.960   -1.175    0.072
## dadas_rscale                 -0.518 0.288 34.304  -1.798   0.960   -1.104    0.067
## abs_diff                      0.021 0.009 39.347   2.465   0.009    0.004    0.038
## abs_sum                       0.103 0.048 34.147   2.160   0.019    0.006    0.201
## abs_diff_bscale               0.140 0.057 39.347   2.465   0.009    0.025    0.256
## abs_sum_bscale                0.692 0.320 34.147   2.160   0.019    0.041    1.343
## abs_diff_rscale               0.186 0.067 37.597   2.779   0.004    0.050    0.321
## abs_sum_rscale                0.704 0.323 34.146   2.180   0.018    0.048    1.360
```

``` r
round(ddsc_mod2_GDI$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.003 -0.463  6.840  1.000  0.009
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
##                                    est    se       z pvalue ci.lower ci.upper
## r_xy1_y2                        -0.521 0.146  -3.563  0.000   -0.808   -0.235
## r_xy1                           -0.448 0.153  -2.921  0.003   -0.748   -0.147
## r_xy2                           -0.303 0.163  -1.852  0.064   -0.623    0.018
## b_11                            -0.476 0.163  -2.921  0.003   -0.795   -0.157
## b_21                            -0.283 0.153  -1.852  0.064   -0.582    0.016
## b_10                            -0.525 0.160  -3.269  0.001   -0.839   -0.210
## b_20                             0.662 0.150   4.406  0.000    0.368    0.957
## res_cov_y1_y2                    0.773 0.193   4.000  0.000    0.394    1.152
## diff_b10_b20                    -1.187 0.053 -22.210  0.000   -1.292   -1.082
## diff_b11_b21                    -0.193 0.054  -3.563  0.000   -0.300   -0.087
## diff_rxy1_rxy2                  -0.145 0.054  -2.668  0.008   -0.252   -0.039
## q_b11_b21                       -0.227 0.077  -2.936  0.003   -0.379   -0.075
## q_rxy1_rxy2                     -0.170 0.064  -2.656  0.008   -0.295   -0.044
## cross_over_point                -6.141 1.746  -3.518  0.000   -9.563   -2.720
## sum_b11_b21                     -0.758 0.311  -2.439  0.015   -1.368   -0.149
## main_effect                     -0.379 0.155  -2.439  0.015   -0.684   -0.074
## interaction_vs_main_effect      -0.186 0.154  -1.204  0.229   -0.489    0.117
## diff_abs_b11_abs_b21             0.193 0.054   3.563  0.000    0.087    0.300
## abs_diff_b11_b21                 0.193 0.054   3.563  0.000    0.087    0.300
## abs_sum_b11_b21                  0.758 0.311   2.439  0.007    0.149    1.368
## dadas                           -0.565 0.305  -1.852  0.968   -1.163    0.033
## q_r_equivalence                  0.070 0.064   1.090  0.862       NA       NA
## q_b_equivalence                  0.127 0.077   1.643  0.950       NA       NA
## cross_over_point_equivalence     6.141 1.746   3.518  1.000       NA       NA
## cross_over_point_minimal_effect  6.141 1.746   3.518  0.000       NA       NA
```

``` r
round(ddsc_sem_GDI$variance_test,3)
```

```
##             est    se      z pvalue ci.lower ci.upper
## cov_y1y2  0.904 0.226  3.992  0.000    0.460    1.348
## var_y1    1.095 0.266  4.123  0.000    0.575    1.616
## var_y2    0.846 0.205  4.123  0.000    0.444    1.248
## var_diff  0.249 0.129  1.940  0.052   -0.003    0.501
## var_ratio 1.295 0.153  8.484  0.000    0.996    1.594
## cor_y1y2  0.939 0.020 46.372  0.000    0.899    0.979
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
## r_xy1y2                      -0.521 0.151 32.000  -3.457   0.002   -0.829   -0.214
## w_11                         -0.071 0.024 33.945  -2.925   0.006   -0.121   -0.022
## w_21                         -0.042 0.024 33.945  -1.737   0.091   -0.092    0.007
## r_xy1                        -0.448 0.153 33.945  -2.925   0.006   -0.759   -0.137
## r_xy2                        -0.303 0.174 33.945  -1.737   0.091   -0.657    0.051
## b_11                         -0.477 0.163 33.945  -2.925   0.006   -0.808   -0.146
## b_21                         -0.283 0.163 33.945  -1.737   0.091   -0.614    0.048
## main_effect                  -0.057 0.024 32.000  -2.366   0.024   -0.106   -0.008
## moderator_effect              0.178 0.008 32.000  21.547   0.000    0.161    0.194
## interaction                   0.029 0.008 32.000   3.457   0.002    0.012    0.046
## q_b11_b21                    -0.228    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                  -0.170    NA     NA      NA      NA       NA       NA
## cross_over_point             -6.141    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.028 0.025 39.675  -1.096   0.280   -0.079    0.024
## interaction_vs_main_bscale   -0.186 0.170 39.675  -1.096   0.280   -0.530    0.157
## interaction_vs_main_rscale   -0.230 0.191 38.444  -1.204   0.236   -0.617    0.157
## dadas                        -0.085 0.049 33.945  -1.737   0.954   -0.184    0.014
## dadas_bscale                 -0.566 0.326 33.945  -1.737   0.954   -1.229    0.096
## dadas_rscale                 -0.605 0.349 33.945  -1.737   0.954   -1.314    0.103
## abs_diff                      0.029 0.008 32.000   3.457   0.001    0.012    0.046
## abs_sum                       0.113 0.048 32.000   2.366   0.012    0.016    0.211
## abs_diff_bscale               0.194 0.056 32.000   3.457   0.001    0.080    0.308
## abs_sum_bscale                0.760 0.321 32.000   2.366   0.012    0.106    1.414
## abs_diff_rscale               0.145 0.060 40.595   2.421   0.010    0.024    0.266
## abs_sum_rscale                0.751 0.323 32.008   2.327   0.013    0.094    1.408
```

``` r
# country-time multilevel model


mod2_GDI_cntry_year<-
  lmer(sti.z.wt~gndr.c+gdi.z.cm+gndr.c:gdi.z.cm+
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
## Formula: sti.z.wt ~ gndr.c + gdi.z.cm + gndr.c:gdi.z.cm + (gndr.c | cntry)
##    Data: diff_dat_cntry_year
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##    -892.3    -857.9     454.1    -908.3       538 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -3.4540 -0.6432 -0.0165  0.5906  4.1486 
## 
## Random effects:
##  Groups   Name        Variance  Std.Dev. Corr  
##  cntry    (Intercept) 0.0172804 0.13145        
##           gndr.c      0.0003271 0.01809  -1.00 
##  Residual             0.0089878 0.09480        
## Number of obs: 546, groups:  cntry, 34
## 
## Fixed effects:
##                   Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)       0.048908   0.022995  33.884262   2.127   0.0408 *  
## gndr.c            0.183743   0.008704 108.265573  21.111   <2e-16 ***
## gdi.z.cm         -0.052885   0.023667  35.748404  -2.234   0.0318 *  
## gndr.c:gdi.z.cm   0.018958   0.010461 176.710882   1.812   0.0716 .  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c gd.z.c
## gndr.c      -0.352              
## gdi.z.cm    -0.007  0.002       
## gndr.c:gd..  0.002 -0.051 -0.297
## optimizer (bobyqa) convergence code: 0 (OK)
## boundary (singular) fit: see help('isSingular')
```

``` r
getFE(mod2_GDI_cntry_year,round=3)
```

```
##                   Est.    SE      df      t     p     LL     UL
## (Intercept)      0.049 0.023  33.884  2.127 0.041  0.002  0.096
## gndr.c           0.184 0.009 108.266 21.111 0.000  0.166  0.201
## gdi.z.cm        -0.053 0.024  35.748 -2.234 0.032 -0.101 -0.005
## gndr.c:gdi.z.cm  0.019 0.010 176.711  1.812 0.072 -0.002  0.040
```

``` r
getVC(mod2_GDI_cntry_year)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.13 0.02
## 2    cntry      gndr.c   <NA>  0.02 0.00
## 3    cntry (Intercept) gndr.c -1.00 0.00
## 4 Residual        <NA>   <NA>  0.09 0.01
```

``` r
r2mlm(mod2_GDI,bargraph = F)
```

```
## $Decompositions
##                        total
## fixed           0.0109047928
## slope variation 0.0004523474
## mean variation  0.0185432765
## sigma2          0.9700995833
## 
## $R2s
##            total
## f   0.0109047928
## v   0.0004523474
## m   0.0185432765
## fv  0.0113571402
## fvm 0.0299004167
```

``` r
ddsc_mod2_GDI_cntry_year<-
  ddsc_ml(model = mod2_GDI_cntry_year,
          predictor = "gdi.z.cm",
          moderator = "gndr.c",moderator_values = c(-0.5,0.5),
          re_cov_test = T)
```

```
## Warning in atanh(r_xy1): NaNs produced
```

```
## Warning in atanh(b_11): NaNs produced
```

```
## Warning in atanh(r_xy1): NaNs produced
```

``` r
round(ddsc_mod2_GDI_cntry_year$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.023        0.152        0.009     0.032 0.720      8.029 0.998   0.954
## 2        0.5         0.017        0.130        0.009     0.026 0.651      8.029 0.998   0.938
```

``` r
round(ddsc_mod2_GDI_cntry_year$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gdi.z.cm gdi.z.cm_scaled
## means_y1           0.142 0.137    1.000           1.000    0.955           0.955   -0.320          -0.320
## means_y1_scaled    0.958 0.920    1.000           1.000    0.955           0.955   -0.320          -0.320
## means_y2          -0.044 0.159    0.955           0.955    1.000           1.000   -0.391          -0.391
## means_y2_scaled   -0.295 1.074    0.955           0.955    1.000           1.000   -0.391          -0.391
## gdi.z.cm           0.000 1.000   -0.320          -0.320   -0.391          -0.391    1.000           1.000
## gdi.z.cm_scaled    0.000 1.000   -0.320          -0.320   -0.391          -0.391    1.000           1.000
## diff_score         0.186 0.050   -0.310          -0.310   -0.579          -0.579    0.372           0.372
## diff_score_scaled  1.253 0.337   -0.310          -0.310   -0.579          -0.579    0.372           0.372
##                   diff_score diff_score_scaled
## means_y1              -0.310            -0.310
## means_y1_scaled       -0.310            -0.310
## means_y2              -0.579            -0.579
## means_y2_scaled       -0.579            -0.579
## gdi.z.cm               0.372             0.372
## gdi.z.cm_scaled        0.372             0.372
## diff_score             1.000             1.000
## diff_score_scaled      1.000             1.000
```

``` r
round(ddsc_mod2_GDI_cntry_year$results,3)
```

```
##                            estimate    SE      df t.ratio p.value ci.lower ci.upper
## r_xy1y2                      -0.379 0.209 176.711  -1.812   0.072   -0.793    0.034
## w_11                         -0.062 0.026  36.302  -2.425   0.020   -0.114   -0.010
## w_21                         -0.043 0.023  36.682  -1.915   0.063   -0.089    0.003
## r_xy1                        -0.457 0.188  36.302  -2.425   0.020   -0.838   -0.075
## r_xy2                        -0.272 0.142  36.682  -1.915   0.063   -0.561    0.016
## b_11                         -0.422 0.174  36.302  -2.425   0.020   -0.774   -0.069
## b_21                         -0.293 0.153  36.682  -1.915   0.063   -0.604    0.017
## main_effect                  -0.053 0.024  35.748  -2.234   0.032   -0.101   -0.005
## moderator_effect              0.184 0.009 108.266  21.111   0.000    0.166    0.201
## interaction                   0.019 0.010 176.711   1.812   0.072   -0.002    0.040
## q_b11_b21                    -0.147    NA      NA      NA      NA       NA       NA
## q_rxy1_rxy2                  -0.214    NA      NA      NA      NA       NA       NA
## cross_over_point             -9.692    NA      NA      NA      NA       NA       NA
## interaction_vs_main          -0.034 0.023  40.141  -1.485   0.145   -0.080    0.012
## interaction_vs_main_bscale   -0.229 0.154  40.141  -1.485   0.145   -0.542    0.083
## interaction_vs_main_rscale   -0.180 0.133  41.738  -1.356   0.182   -0.449    0.088
## dadas                        -0.087 0.045  36.682  -1.915   0.968   -0.179    0.005
## dadas_bscale                 -0.587 0.306  36.682  -1.915   0.968   -1.208    0.034
## dadas_rscale                 -0.545 0.285  36.682  -1.915   0.968   -1.122    0.032
## abs_diff                      0.019 0.010 176.711   1.812   0.036   -0.002    0.040
## abs_sum                       0.106 0.047  35.748   2.234   0.016    0.010    0.202
## abs_diff_bscale               0.128 0.071 176.711   1.812   0.036   -0.011    0.268
## abs_sum_bscale                0.715 0.320  35.748   2.234   0.016    0.066    1.364
## abs_diff_rscale               0.184 0.082  74.429   2.247   0.014    0.021    0.347
## abs_sum_rscale                0.729 0.324  35.746   2.254   0.015    0.073    1.385
```

``` r
round(ddsc_mod2_GDI_cntry_year$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.003 -1.000  6.458  1.000  0.011
```

### Comparisons across models


``` r
# full multilevel model (individuals as entries)
round(ddsc_mod2_GDI$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df t.ratio p.value ci.lower ci.upper
## r_xy1             -0.4450 0.1834 34.2364 -2.4261  0.0207  -0.8176  -0.0723
## r_xy2             -0.2591 0.1441 34.3041 -1.7980  0.0810  -0.5519   0.0337
## b_11              -0.4163 0.1716 34.2364 -2.4261  0.0207  -0.7649  -0.0677
## b_21              -0.2759 0.1534 34.3041 -1.7980  0.0810  -0.5875   0.0358
## main_effect       -0.0517 0.0239 34.1471 -2.1596  0.0379  -0.1003  -0.0031
## moderator_effect   0.1875 0.0080 34.0325 23.3076  0.0000   0.1712   0.2039
## interaction        0.0210 0.0085 39.3467  2.4645  0.0182   0.0038   0.0382
## q_b11_b21         -0.1600     NA      NA      NA      NA       NA       NA
```

``` r
# country-level structural path model (country-gender means as entries)
round(ddsc_sem_GDI$results[c("r_xy1","r_xy2","b_11","b_21","q_b11_b21","diff_b11_b21"),],4)
```

```
##                  est     se       z pvalue ci.lower ci.upper
## r_xy1        -0.4479 0.1533 -2.9213 0.0035  -0.7485  -0.1474
## r_xy2        -0.3027 0.1635 -1.8520 0.0640  -0.6231   0.0177
## b_11         -0.4758 0.1629 -2.9213 0.0035  -0.7951  -0.1566
## b_21         -0.2826 0.1526 -1.8520 0.0640  -0.5817   0.0165
## q_b11_b21    -0.2271 0.0774 -2.9356 0.0033  -0.3787  -0.0755
## diff_b11_b21 -0.1932 0.0542 -3.5630 0.0004  -0.2996  -0.0869
```

``` r
# multilevel model (country-gender means as entries)
round(ddsc_mod2_GDI_ri$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df t.ratio p.value ci.lower ci.upper
## r_xy1             -0.4479 0.1531 33.9454 -2.9250  0.0061  -0.7592  -0.1367
## r_xy2             -0.3027 0.1743 33.9454 -1.7371  0.0914  -0.6569   0.0515
## b_11              -0.4768 0.1630 33.9454 -2.9250  0.0061  -0.8081  -0.1455
## b_21              -0.2832 0.1630 33.9454 -1.7371  0.0914  -0.6145   0.0481
## main_effect       -0.0567 0.0240 32.0000 -2.3663  0.0242  -0.1056  -0.0079
## moderator_effect   0.1776 0.0082 32.0000 21.5471  0.0000   0.1608   0.1944
## interaction        0.0289 0.0084 32.0000  3.4566  0.0016   0.0119   0.0460
## q_b11_b21         -0.2277     NA      NA      NA      NA       NA       NA
```

``` r
# multilevel model (country-year-gender means as entries)
round(ddsc_mod2_GDI_cntry_year$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE       df t.ratio p.value ci.lower ci.upper
## r_xy1             -0.4566 0.1883  36.3016 -2.4254  0.0204  -0.8384  -0.0749
## r_xy2             -0.2725 0.1423  36.6822 -1.9148  0.0633  -0.5609   0.0159
## b_11              -0.4216 0.1738  36.3016 -2.4254  0.0204  -0.7740  -0.0692
## b_21              -0.2934 0.1532  36.6822 -1.9148  0.0633  -0.6040   0.0172
## main_effect       -0.0529 0.0237  35.7484 -2.2345  0.0318  -0.1009  -0.0049
## moderator_effect   0.1837 0.0087 108.2656 21.1111  0.0000   0.1665   0.2010
## interaction        0.0190 0.0105 176.7109  1.8123  0.0716  -0.0017   0.0396
## q_b11_b21         -0.1473     NA       NA      NA      NA       NA       NA
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
## Time difference of 1.373461 hours
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
##                    Estimate          SE         2.5%        97.5%
## X.Intercept.     0.05047502 0.023766653  0.006680508  0.097926519
## gndr.c           0.18775638 0.008492713  0.171279440  0.203461548
## gdi.z.cm        -0.05060097 0.024827175 -0.100334589 -0.004015406
## gndr.c.gdi.z.cm  0.02065167 0.008916716  0.003563832  0.038343417
## w11             -0.06092680 0.026516603 -0.115562057 -0.011563545
## w21             -0.04027513 0.023862124 -0.089355793  0.002754313
## b11             -0.40800206 0.177570920 -0.773872164 -0.077436364
## b21             -0.26970622 0.159794951 -0.598379458  0.018444514
## r_xy1           -0.43614057 0.189817380 -0.827243483 -0.082776885
## r_xy2           -0.25336016 0.150110277 -0.562113543  0.017326650
## q_b             -0.17064324 0.093362407 -0.383989262 -0.026117778
## q               -0.23477693 0.151723938 -0.561302927 -0.049245804
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
## [1] -0.1706432
## 
## $se
## [1] 0.09336241
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
## [1] -0.7566562
## 
## $p_low
## [1] 0.7753721
## 
## $z_high
## [1] -2.898846
## 
## $p_high
## [1] 0.001872694
## 
## $ci_level
## [1] 0.9
## 
## $ci_lower
## [1] -0.3242107
## 
## $ci_upper
## [1] -0.01707575
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
## [1] -0.2347769
## 
## $se
## [1] 0.1517239
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
## [1] -0.8883037
## 
## $p_low
## [1] 0.8128113
## 
## $z_high
## [1] -2.206487
## 
## $p_high
## [1] 0.01367495
## 
## $ci_level
## [1] 0.9
## 
## $ci_lower
## [1] -0.4843406
## 
## $ci_upper
## [1] 0.01478674
## 
## $equivalent
## [1] FALSE
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
  xlab("Gender Development Index (Average 2002-2023)")+
  #ylim(c(min.y.pred,max.y.pred))+
  #xlim(c(0.60,1.00))+
  ylab("Value stimulation (Average 2002-2023)")+
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
  xlab("Gender Development Index (Average 2002-2023)")+
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
## Warning in geom_text(inherit.aes = F, aes(x = 0.9, y = -0.45, label = coef_q, : All aesthetics have length 1, but the data has 262 rows.
## ℹ Please consider using `annotate()` or provide this layer with data containing a single row.
```

``` r
pflag_comb
```

![](Analysis_stimulation_files/figure-html/unnamed-chunk-35-1.png)<!-- -->

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
##   1455824   1455913   -727904   1455808    492335 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -4.8687 -0.6919 -0.0539  0.6445  5.0885 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr  
##  cntry    (Intercept) 0.021224 0.14569        
##           gndr.c      0.002131 0.04616  -0.46 
##  Residual             0.989875 0.99492        
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##                      Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)          0.051291   0.025046 33.947776   2.048   0.0484 *  
## gndr.c               0.187941   0.008534 32.878035  22.022   <2e-16 ***
## log_gdp.z.cm         0.010398   0.025133 34.038510   0.414   0.6817    
## gndr.c:log_gdp.z.cm -0.009486   0.008719 35.000734  -1.088   0.2841    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c lg_g..
## gndr.c      -0.423              
## lg_gdp.z.cm  0.022 -0.010       
## gndr.c:l_.. -0.010 -0.023 -0.415
```

``` r
getFE(mod2_log_GDP,round=3)
```

```
##                       Est.    SE     df      t     p     LL    UL
## (Intercept)          0.051 0.025 33.948  2.048 0.048  0.000 0.102
## gndr.c               0.188 0.009 32.878 22.022 0.000  0.171 0.205
## log_gdp.z.cm         0.010 0.025 34.039  0.414 0.682 -0.041 0.061
## gndr.c:log_gdp.z.cm -0.009 0.009 35.001 -1.088 0.284 -0.027 0.008
```

``` r
getVC(mod2_log_GDP)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.15 0.02
## 2    cntry      gndr.c   <NA>  0.05 0.00
## 3    cntry (Intercept) gndr.c -0.46 0.00
## 4 Residual        <NA>   <NA>  0.99 0.99
```

``` r
r2mlm(mod2_log_GDP,bargraph = F)
```

```
## $Decompositions
##                        total
## fixed           0.0086158804
## slope variation 0.0005187002
## mean variation  0.0210321254
## sigma2          0.9698332940
## 
## $R2s
##            total
## f   0.0086158804
## v   0.0005187002
## m   0.0210321254
## fv  0.0091345806
## fvm 0.0301667060
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
## Time difference of 33.8762 secs
```

``` r
round(ddsc_mod2_log_GDP$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.025        0.158         0.99     1.015 0.025   7802.647 0.995   0.995
## 2        0.5         0.019        0.137         0.99     1.009 0.019   6678.029 0.992   0.992
```

``` r
round(ddsc_mod2_log_GDP$descriptives,3)
```

```
##                          M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled log_gdp.z.cm
## means_y1             0.099 0.140    1.000           1.000    0.939           0.939        0.083
## means_y1_scaled      0.662 0.934    1.000           1.000    0.939           0.939        0.083
## means_y2            -0.078 0.159    0.939           0.939    1.000           1.000        0.174
## means_y2_scaled     -0.525 1.062    0.939           0.939    1.000           1.000        0.174
## log_gdp.z.cm        -0.024 1.012    0.083           0.083    0.174           0.174        1.000
## log_gdp.z.cm_scaled  0.000 1.000    0.083           0.083    0.174           0.174        1.000
## diff_score           0.178 0.055   -0.173          -0.173   -0.501          -0.501       -0.292
## diff_score_scaled    1.187 0.371   -0.173          -0.173   -0.501          -0.501       -0.292
##                     log_gdp.z.cm_scaled diff_score diff_score_scaled
## means_y1                          0.083     -0.173            -0.173
## means_y1_scaled                   0.083     -0.173            -0.173
## means_y2                          0.174     -0.501            -0.501
## means_y2_scaled                   0.174     -0.501            -0.501
## log_gdp.z.cm                      1.000     -0.292            -0.292
## log_gdp.z.cm_scaled               1.000     -0.292            -0.292
## diff_score                       -0.292      1.000             1.000
## diff_score_scaled                -0.292      1.000             1.000
```

``` r
round(ddsc_mod2_log_GDP$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                       0.171 0.157 35.001   1.088   0.284   -0.148    0.490
## w_11                          0.015 0.027 34.048   0.556   0.582   -0.040    0.070
## w_21                          0.006 0.024 34.084   0.239   0.813   -0.042    0.054
## r_xy1                         0.108 0.195 34.048   0.556   0.582   -0.288    0.505
## r_xy2                         0.036 0.149 34.084   0.239   0.813   -0.267    0.338
## b_11                          0.101 0.182 34.048   0.556   0.582   -0.269    0.472
## b_21                          0.038 0.158 34.084   0.239   0.813   -0.284    0.360
## main_effect                   0.010 0.025 34.039   0.414   0.682   -0.041    0.061
## moderator_effect              0.188 0.009 32.878  22.022   0.000    0.171    0.205
## interaction                  -0.009 0.009 35.001  -1.088   0.284   -0.027    0.008
## q_b11_b21                     0.064    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.073    NA     NA      NA      NA       NA       NA
## cross_over_point             19.813    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.001 0.023 34.245  -0.040   0.968   -0.048    0.046
## interaction_vs_main_bscale   -0.006 0.154 34.245  -0.040   0.968   -0.318    0.306
## interaction_vs_main_rscale    0.001 0.134 34.302   0.006   0.995   -0.271    0.273
## dadas                        -0.011 0.047 34.084  -0.239   0.594   -0.107    0.085
## dadas_bscale                 -0.076 0.317 34.084  -0.239   0.594   -0.720    0.568
## dadas_rscale                 -0.071 0.298 34.084  -0.239   0.594   -0.676    0.534
## abs_diff                      0.009 0.009 35.001   1.088   0.142   -0.008    0.027
## abs_sum                       0.021 0.050 34.039   0.414   0.341   -0.081    0.123
## abs_diff_bscale               0.064 0.058 35.001   1.088   0.142   -0.055    0.182
## abs_sum_bscale                0.139 0.337 34.039   0.414   0.341   -0.545    0.823
## abs_diff_rscale               0.073 0.071 34.468   1.032   0.155   -0.070    0.216
## abs_sum_rscale                0.144 0.340 34.038   0.424   0.337   -0.546    0.834
```

``` r
round(ddsc_mod2_log_GDP$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.003 -0.463  6.840  1.000  0.009
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
##                                    est    se       z pvalue ci.lower ci.upper
## r_xy1_y2                         0.292 0.164   1.778  0.075   -0.030    0.613
## r_xy1                            0.174 0.169   1.033  0.302   -0.157    0.505
## r_xy2                            0.083 0.171   0.484  0.628   -0.252    0.418
## b_11                             0.185 0.179   1.033  0.302   -0.166    0.537
## b_21                             0.077 0.160   0.484  0.628   -0.236    0.390
## b_10                            -0.525 0.177  -2.968  0.003   -0.871   -0.178
## b_20                             0.662 0.157   4.213  0.000    0.354    0.970
## res_cov_y1_y2                    0.890 0.223   3.999  0.000    0.454    1.326
## diff_b10_b20                    -1.187 0.060 -19.814  0.000   -1.304   -1.069
## diff_b11_b21                     0.108 0.061   1.778  0.075   -0.011    0.227
## diff_rxy1_rxy2                   0.092 0.058   1.588  0.112   -0.021    0.205
## q_b11_b21                        0.110 0.064   1.724  0.085   -0.015    0.235
## q_rxy1_rxy2                      0.093 0.059   1.586  0.113   -0.022    0.209
## cross_over_point                10.979 6.200   1.771  0.077   -1.173   23.132
## sum_b11_b21                      0.262 0.334   0.786  0.432   -0.392    0.917
## main_effect                      0.131 0.167   0.786  0.432   -0.196    0.459
## interaction_vs_main_effect      -0.023 0.158  -0.147  0.883   -0.332    0.286
## diff_abs_b11_abs_b21             0.108 0.061   1.778  0.075   -0.011    0.227
## abs_diff_b11_b21                 0.108 0.061   1.778  0.038   -0.011    0.227
## abs_sum_b11_b21                  0.262 0.334   0.786  0.216   -0.392    0.917
## dadas                           -0.154 0.319  -0.484  0.686   -0.780    0.471
## q_r_equivalence                 -0.007 0.059  -0.113  0.455       NA       NA
## q_b_equivalence                  0.010 0.064   0.158  0.563       NA       NA
## cross_over_point_equivalence    10.979 6.200   1.771  0.962       NA       NA
## cross_over_point_minimal_effect 10.979 6.200   1.771  0.038       NA       NA
```

``` r
round(ddsc_sem_log_GDP$variance_test,3)
```

```
##             est    se      z pvalue ci.lower ci.upper
## cov_y1y2  0.904 0.226  3.992  0.000    0.460    1.348
## var_y1    1.095 0.266  4.123  0.000    0.575    1.616
## var_y2    0.846 0.205  4.123  0.000    0.444    1.248
## var_diff  0.249 0.129  1.940  0.052   -0.003    0.501
## var_ratio 1.295 0.153  8.484  0.000    0.996    1.594
## cor_y1y2  0.939 0.020 46.372  0.000    0.899    0.979
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
## r_xy1y2                       0.292 0.169 32.000   1.725   0.094   -0.053    0.636
## w_11                          0.028 0.026 34.118   1.059   0.297   -0.025    0.081
## w_21                          0.012 0.026 34.118   0.441   0.662   -0.042    0.065
## r_xy1                         0.174 0.165 34.118   1.059   0.297   -0.160    0.509
## r_xy2                         0.083 0.187 34.118   0.441   0.662   -0.298    0.464
## b_11                          0.186 0.175 34.118   1.059   0.297   -0.171    0.542
## b_21                          0.077 0.175 34.118   0.441   0.662   -0.279    0.434
## main_effect                   0.020 0.026 32.000   0.762   0.451   -0.033    0.072
## moderator_effect              0.178 0.009 32.000  19.222   0.000    0.159    0.196
## interaction                  -0.016 0.009 32.000  -1.725   0.094   -0.035    0.003
## q_b11_b21                     0.110    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.093    NA     NA      NA      NA       NA       NA
## cross_over_point             10.979    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.003 0.027 40.335  -0.126   0.900   -0.059    0.052
## interaction_vs_main_bscale   -0.023 0.184 40.335  -0.126   0.900   -0.394    0.348
## interaction_vs_main_rscale   -0.037 0.206 39.003  -0.179   0.859   -0.454    0.380
## dadas                        -0.023 0.052 34.118  -0.441   0.669   -0.130    0.083
## dadas_bscale                 -0.155 0.351 34.118  -0.441   0.669   -0.867    0.558
## dadas_rscale                 -0.165 0.375 34.118  -0.441   0.669   -0.927    0.596
## abs_diff                      0.016 0.009 32.000   1.725   0.047   -0.003    0.035
## abs_sum                       0.039 0.052 32.000   0.762   0.226   -0.066    0.144
## abs_diff_bscale               0.108 0.063 32.000   1.725   0.047   -0.020    0.236
## abs_sum_bscale                0.263 0.345 32.000   0.762   0.226   -0.440    0.966
## abs_diff_rscale               0.092 0.067 39.916   1.371   0.089   -0.044    0.227
## abs_sum_rscale                0.257 0.346 32.009   0.742   0.232   -0.449    0.963
```

``` r
# country-time multilevel model


mod2_log_GDP_cntry_year<-
  lmer(sti.z.wt~gndr.c+log_gdp.z.cm+gndr.c:log_gdp.z.cm+
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
## Formula: sti.z.wt ~ gndr.c + log_gdp.z.cm + gndr.c:log_gdp.z.cm + (gndr.c |      cntry)
##    Data: diff_dat_cntry_year
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##    -888.6    -854.2     452.3    -904.6       538 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -3.4102 -0.6129 -0.0070  0.5691  4.0640 
## 
## Random effects:
##  Groups   Name        Variance  Std.Dev. Corr  
##  cntry    (Intercept) 0.0196560 0.14020        
##           gndr.c      0.0004189 0.02047  -1.00 
##  Residual             0.0089784 0.09475        
## Number of obs: 546, groups:  cntry, 34
## 
## Fixed effects:
##                       Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)           0.048852   0.024468  33.595345   1.997    0.054 .  
## gndr.c                0.186848   0.009029  95.204046  20.694   <2e-16 ***
## log_gdp.z.cm          0.014374   0.024668  34.283392   0.583    0.564    
## gndr.c:log_gdp.z.cm  -0.014656   0.009592 113.973515  -1.528    0.129    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c lg_g..
## gndr.c      -0.384              
## lg_gdp.z.cm  0.010 -0.007       
## gndr.c:l_.. -0.006 -0.186 -0.363
## optimizer (bobyqa) convergence code: 0 (OK)
## boundary (singular) fit: see help('isSingular')
```

``` r
getFE(mod2_log_GDP_cntry_year,round=3)
```

```
##                       Est.    SE      df      t     p     LL    UL
## (Intercept)          0.049 0.024  33.595  1.997 0.054 -0.001 0.099
## gndr.c               0.187 0.009  95.204 20.694 0.000  0.169 0.205
## log_gdp.z.cm         0.014 0.025  34.283  0.583 0.564 -0.036 0.064
## gndr.c:log_gdp.z.cm -0.015 0.010 113.974 -1.528 0.129 -0.034 0.004
```

``` r
getVC(mod2_log_GDP_cntry_year)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.14 0.02
## 2    cntry      gndr.c   <NA>  0.02 0.00
## 3    cntry (Intercept) gndr.c -1.00 0.00
## 4 Residual        <NA>   <NA>  0.09 0.01
```

``` r
r2mlm(mod2_log_GDP,bargraph = F)
```

```
## $Decompositions
##                        total
## fixed           0.0086158804
## slope variation 0.0005187002
## mean variation  0.0210321254
## sigma2          0.9698332940
## 
## $R2s
##            total
## f   0.0086158804
## v   0.0005187002
## m   0.0210321254
## fv  0.0091345806
## fvm 0.0301667060
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
## 1       -0.5         0.023        0.152        0.009     0.032 0.720      8.029 0.998   0.954
## 2        0.5         0.017        0.130        0.009     0.026 0.651      8.029 0.998   0.938
```

``` r
round(ddsc_mod2_log_GDP_cntry_year$descriptives,3)
```

```
##                          M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled log_gdp.z.cm
## means_y1             0.142 0.137    1.000           1.000    0.955           0.955        0.065
## means_y1_scaled      0.958 0.920    1.000           1.000    0.955           0.955        0.065
## means_y2            -0.044 0.159    0.955           0.955    1.000           1.000        0.117
## means_y2_scaled     -0.295 1.074    0.955           0.955    1.000           1.000        0.117
## log_gdp.z.cm        -0.024 1.012    0.065           0.065    0.117           0.117        1.000
## log_gdp.z.cm_scaled  0.000 1.000    0.065           0.065    0.117           0.117        1.000
## diff_score           0.186 0.050   -0.310          -0.310   -0.579          -0.579       -0.195
## diff_score_scaled    1.253 0.337   -0.310          -0.310   -0.579          -0.579       -0.195
##                     log_gdp.z.cm_scaled diff_score diff_score_scaled
## means_y1                          0.065     -0.310            -0.310
## means_y1_scaled                   0.065     -0.310            -0.310
## means_y2                          0.117     -0.579            -0.579
## means_y2_scaled                   0.117     -0.579            -0.579
## log_gdp.z.cm                      1.000     -0.195            -0.195
## log_gdp.z.cm_scaled               1.000     -0.195            -0.195
## diff_score                       -0.195      1.000             1.000
## diff_score_scaled                -0.195      1.000             1.000
```

``` r
round(ddsc_mod2_log_GDP_cntry_year$results,3)
```

```
##                            estimate    SE      df t.ratio p.value ci.lower ci.upper
## r_xy1y2                       0.293 0.192 113.974   1.528   0.129   -0.087    0.674
## w_11                          0.022 0.027  34.356   0.810   0.423   -0.033    0.076
## w_21                          0.007 0.023  34.478   0.302   0.765   -0.040    0.054
## r_xy1                         0.159 0.196  34.356   0.810   0.423   -0.240    0.557
## r_xy2                         0.044 0.147  34.478   0.302   0.765   -0.254    0.342
## b_11                          0.147 0.181  34.356   0.810   0.423   -0.221    0.515
## b_21                          0.048 0.158  34.478   0.302   0.765   -0.273    0.368
## main_effect                   0.014 0.025  34.283   0.583   0.564   -0.036    0.064
## moderator_effect              0.187 0.009  95.204  20.694   0.000    0.169    0.205
## interaction                  -0.015 0.010 113.974  -1.528   0.129   -0.034    0.004
## q_b11_b21                     0.100    NA      NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.116    NA      NA      NA      NA       NA       NA
## cross_over_point             12.749    NA      NA      NA      NA       NA       NA
## interaction_vs_main           0.000 0.023  35.595   0.012   0.990   -0.046    0.047
## interaction_vs_main_bscale    0.002 0.155  35.595   0.012   0.990   -0.313    0.317
## interaction_vs_main_rscale    0.013 0.133  36.248   0.099   0.922   -0.256    0.282
## dadas                        -0.014 0.047  34.478  -0.302   0.618   -0.109    0.081
## dadas_bscale                 -0.095 0.316  34.478  -0.302   0.618   -0.737    0.546
## dadas_rscale                 -0.088 0.293  34.478  -0.302   0.618   -0.684    0.507
## abs_diff                      0.015 0.010 113.974   1.528   0.065   -0.004    0.034
## abs_sum                       0.029 0.049  34.283   0.583   0.282   -0.071    0.129
## abs_diff_bscale               0.099 0.065 113.974   1.528   0.065   -0.029    0.228
## abs_sum_bscale                0.194 0.334  34.283   0.583   0.282   -0.483    0.872
## abs_diff_rscale               0.115 0.078  52.276   1.464   0.075   -0.043    0.272
## abs_sum_rscale                0.203 0.337  34.282   0.602   0.276   -0.482    0.888
```

``` r
round(ddsc_mod2_log_GDP_cntry_year$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.003 -1.000  6.458  1.000  0.011
```

### Comparisons across models


``` r
# full multilevel model (individuals as entries)
round(ddsc_mod2_log_GDP$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df t.ratio p.value ci.lower ci.upper
## r_xy1              0.1084 0.1949 34.0481  0.5560  0.5819  -0.2878   0.5045
## r_xy2              0.0356 0.1488 34.0835  0.2390  0.8125  -0.2669   0.3380
## b_11               0.1014 0.1824 34.0481  0.5560  0.5819  -0.2692   0.4720
## b_21               0.0379 0.1584 34.0835  0.2390  0.8125  -0.2841   0.3598
## main_effect        0.0104 0.0251 34.0385  0.4137  0.6817  -0.0407   0.0615
## moderator_effect   0.1879 0.0085 32.8780 22.0217  0.0000   0.1706   0.2053
## interaction       -0.0095 0.0087 35.0007 -1.0879  0.2841  -0.0272   0.0082
## q_b11_b21          0.0639     NA      NA      NA      NA       NA       NA
```

``` r
# country-level structural path model (country-gender means as entries)
round(ddsc_sem_log_GDP$results[c("r_xy1","r_xy2","b_11","b_21","q_b11_b21","diff_b11_b21"),],4)
```

```
##                 est     se      z pvalue ci.lower ci.upper
## r_xy1        0.1744 0.1689 1.0329 0.3016  -0.1565   0.5054
## r_xy2        0.0827 0.1709 0.4839 0.6285  -0.2523   0.4177
## b_11         0.1853 0.1794 1.0329 0.3016  -0.1663   0.5369
## b_21         0.0772 0.1596 0.4839 0.6285  -0.2355   0.3899
## q_b11_b21    0.1101 0.0639 1.7240 0.0847  -0.0151   0.2353
## diff_b11_b21 0.1081 0.0608 1.7779 0.0754  -0.0111   0.2273
```

``` r
# multilevel model (country-gender means as entries)
round(ddsc_mod2_log_GDP_ri$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df t.ratio p.value ci.lower ci.upper
## r_xy1              0.1744 0.1647 34.1179  1.0589  0.2971  -0.1603   0.5091
## r_xy2              0.0827 0.1874 34.1179  0.4412  0.6619  -0.2982   0.4636
## b_11               0.1857 0.1754 34.1179  1.0589  0.2971  -0.1706   0.5420
## b_21               0.0774 0.1754 34.1179  0.4412  0.6619  -0.2789   0.4337
## main_effect        0.0196 0.0258 32.0000  0.7624  0.4514  -0.0328   0.0721
## moderator_effect   0.1776 0.0092 32.0000 19.2220  0.0000   0.1588   0.1964
## interaction       -0.0162 0.0094 32.0000 -1.7248  0.0942  -0.0353   0.0029
## q_b11_b21          0.1103     NA      NA      NA      NA       NA       NA
```

``` r
# multilevel model (country-year-gender means as entries)
round(ddsc_mod2_log_GDP_cntry_year$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE       df t.ratio p.value ci.lower ci.upper
## r_xy1              0.1589 0.1961  34.3562  0.8102  0.4234  -0.2395   0.5573
## r_xy2              0.0442 0.1466  34.4780  0.3016  0.7647  -0.2536   0.3421
## b_11               0.1467 0.1811  34.3562  0.8102  0.4234  -0.2211   0.5145
## b_21               0.0476 0.1579  34.4780  0.3016  0.7647  -0.2731   0.3683
## main_effect        0.0144 0.0247  34.2834  0.5827  0.5639  -0.0357   0.0645
## moderator_effect   0.1868 0.0090  95.2040 20.6935  0.0000   0.1689   0.2048
## interaction       -0.0147 0.0096 113.9735 -1.5279  0.1293  -0.0337   0.0043
## q_b11_b21          0.1001     NA       NA      NA      NA       NA       NA
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
## Time difference of 1.367952 hours
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
##                        Estimate          SE         2.5%       97.5%
## X.Intercept.         0.05064800 0.025297910  0.004122526 0.101600277
## gndr.c               0.18819520 0.009028393  0.170490208 0.204816200
## log_gdp.z.cm         0.01014366 0.024492871 -0.038052499 0.057338873
## gndr.c.log_gdp.z.cm -0.00906295 0.008751266 -0.024944863 0.007797803
## w11                  0.01467514 0.026648122 -0.038574801 0.067963773
## w21                  0.00561219 0.022977628 -0.038085149 0.049332232
## b11                  0.09827345 0.178451648 -0.258319777 0.455125786
## b21                  0.03758255 0.153871845 -0.255040776 0.330357924
## r_xy1                0.10505103 0.190758848 -0.276135210 0.486514256
## r_xy2                0.03530479 0.144546152 -0.239583549 0.310335959
## q_b                  0.06369926 0.062605930 -0.054561653 0.183157563
## q                    0.07397368 0.077784875 -0.072840154 0.234810467
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
## [1] 0.06369926
## 
## $se
## [1] 0.06260593
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
## [1] 2.614756
## 
## $p_low
## [1] 0.004464555
## 
## $z_high
## [1] -0.5798291
## 
## $p_high
## [1] 0.2810149
## 
## $ci_level
## [1] 0.9
## 
## $ci_lower
## [1] -0.03927833
## 
## $ci_upper
## [1] 0.1666769
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
## [1] 0.07397368
## 
## $se
## [1] 0.07778487
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
## [1] 2.2366
## 
## $p_low
## [1] 0.01265624
## 
## $z_high
## [1] -0.3345936
## 
## $p_high
## [1] 0.3689658
## 
## $ci_level
## [1] 0.9
## 
## $ci_lower
## [1] -0.05397106
## 
## $ci_upper
## [1] 0.2019184
## 
## $equivalent
## [1] FALSE
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
  xlab("log(GDP per capita) (Average 2002-2023)")+
  #ylim(c(min.y.pred,max.y.pred))+
  #xlim(c(0.60,1.00))+
  ylab("Value stimulation (Average 2002-2023)")+
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
  xlab("log(GDP per capita) (Average 2002-2023)")+
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

![](Analysis_stimulation_files/figure-html/unnamed-chunk-41-1.png)<!-- -->

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


## mod3: fixed effect of time (year)

* Year coded year-2002 (2002 is a zero-point)


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
## 1455643.7 1455721.5 -727814.9 1455629.7    492336 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -4.8715 -0.6918 -0.0481  0.6445  5.1263 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr  
##  cntry    (Intercept) 0.021893 0.14796        
##           gndr.c      0.002238 0.04731  -0.46 
##  Residual             0.989512 0.99474        
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##              Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept) 4.926e-02  2.543e-02 3.395e+01   1.937   0.0611 .  
## gndr.c      1.877e-01  8.714e-03 3.339e+01  21.536   <2e-16 ***
## essround.c  6.454e-03  4.818e-04 4.883e+05  13.396   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##            (Intr) gndr.c
## gndr.c     -0.429       
## essround.c -0.005  0.000
```

``` r
getFE(mod3,round=3)
```

```
##              Est.    SE         df      t     p     LL    UL
## (Intercept) 0.049 0.025     33.953  1.937 0.061 -0.002 0.101
## gndr.c      0.188 0.009     33.393 21.536 0.000  0.170 0.205
## essround.c  0.006 0.000 488292.082 13.396 0.000  0.006 0.007
```

``` r
getVC(mod3)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.15 0.02
## 2    cntry      gndr.c   <NA>  0.05 0.00
## 3    cntry (Intercept) gndr.c -0.46 0.00
## 4 Residual        <NA>   <NA>  0.99 0.99
```

``` r
r2mlm(mod3,bargraph = F)
```

```
## $Decompositions
##                        total
## fixed           0.0089730414
## slope variation 0.0005444409
## mean variation  0.0216840838
## sigma2          0.9687984339
## 
## $R2s
##            total
## f   0.0089730414
## v   0.0005444409
## m   0.0216840838
## fv  0.0095174822
## fvm 0.0312015661
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
## mod2    6 1455821 1455888 -727905   1455809                         
## mod3    7 1455644 1455721 -727815   1455630 179.41  1  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## mod4: random effect of time (year)


``` r
mod4<-lmer(sti.z~gndr.c+year.c+(gndr.c+year.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod4)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: sti.z ~ gndr.c + year.c + (gndr.c + year.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1454815.6 1454926.7 -727397.8 1454795.6    492333 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -5.0218 -0.6926 -0.0433  0.6435  5.0702 
## 
## Random effects:
##  Groups   Name        Variance  Std.Dev. Corr        
##  cntry    (Intercept) 0.0361634 0.19017              
##           gndr.c      0.0022694 0.04764  -0.44       
##           year.c      0.0002454 0.01567  -0.42 -0.09 
##  Residual             0.9875096 0.99374              
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##              Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)  0.021603   0.033107 29.658863   0.653   0.5191    
## gndr.c       0.187198   0.008764 33.101338  21.359   <2e-16 ***
## year.c       0.004788   0.002718 24.919164   1.762   0.0904 .  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##        (Intr) gndr.c
## gndr.c -0.405       
## year.c -0.432 -0.081
```

``` r
getFE(mod4,round=3)
```

```
##              Est.    SE     df      t     p     LL    UL
## (Intercept) 0.022 0.033 29.659  0.653 0.519 -0.046 0.089
## gndr.c      0.187 0.009 33.101 21.359 0.000  0.169 0.205
## year.c      0.005 0.003 24.919  1.762 0.090 -0.001 0.010
```

``` r
getVC(mod4)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.19 0.04
## 2    cntry      gndr.c   <NA>  0.05 0.00
## 3    cntry      year.c   <NA>  0.02 0.00
## 4    cntry (Intercept) gndr.c -0.44 0.00
## 5    cntry (Intercept) year.c -0.42 0.00
## 6    cntry      gndr.c year.c -0.09 0.00
## 7 Residual        <NA>   <NA>  0.99 0.99
```

``` r
r2mlm(mod4,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.009245721
## slope variation 0.010001176
## mean variation  0.034949457
## sigma2          0.945803647
## 
## $R2s
##           total
## f   0.009245721
## v   0.010001176
## m   0.034949457
## fv  0.019246897
## fvm 0.054196353
```

``` r
anova(mod2,mod3,mod4)
```

```
## Data: diff_dat
## Models:
## mod2: sti.z ~ gndr.c + (gndr.c | cntry)
## mod3: sti.z ~ gndr.c + essround.c + (gndr.c | cntry)
## mod4: sti.z ~ gndr.c + year.c + (gndr.c + year.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod2    6 1455821 1455888 -727905   1455809                         
## mod3    7 1455644 1455721 -727815   1455630 179.41  1  < 2.2e-16 ***
## mod4   10 1454816 1454927 -727398   1454796 834.09  3  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## mod5: fixed interaction between time and gender


``` r
mod5<-lmer(sti.z~gndr.c+year.c+
             gndr.c:year.c+(gndr.c+year.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod5)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: sti.z ~ gndr.c + year.c + gndr.c:year.c + (gndr.c + year.c |      cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##   1454772   1454894   -727375   1454750    492332 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -5.0025 -0.6926 -0.0394  0.6432  5.0861 
## 
## Random effects:
##  Groups   Name        Variance  Std.Dev. Corr        
##  cntry    (Intercept) 0.0358706 0.18940              
##           gndr.c      0.0023906 0.04889  -0.40       
##           year.c      0.0002443 0.01563  -0.42 -0.17 
##  Residual             0.9874173 0.99369              
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##                 Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)    2.171e-02  3.298e-02  2.970e+01   0.658   0.5154    
## gndr.c         2.202e-01  1.021e-02  5.442e+01  21.570  < 2e-16 ***
## year.c         4.755e-03  2.712e-03  2.489e+01   1.753   0.0919 .  
## gndr.c:year.c -3.139e-03  4.635e-04  6.489e+04  -6.773 1.27e-11 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c year.c
## gndr.c      -0.327              
## year.c      -0.429 -0.141       
## gndr.c:yr.c  0.000 -0.478  0.001
```

``` r
getFE(mod5,round=3)
```

```
##                 Est.    SE        df      t     p     LL     UL
## (Intercept)    0.022 0.033    29.696  0.658 0.515 -0.046  0.089
## gndr.c         0.220 0.010    54.425 21.570 0.000  0.200  0.241
## year.c         0.005 0.003    24.895  1.753 0.092 -0.001  0.010
## gndr.c:year.c -0.003 0.000 64893.363 -6.773 0.000 -0.004 -0.002
```

``` r
getVC(mod5)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.19 0.04
## 2    cntry      gndr.c   <NA>  0.05 0.00
## 3    cntry      year.c   <NA>  0.02 0.00
## 4    cntry (Intercept) gndr.c -0.40 0.00
## 5    cntry (Intercept) year.c -0.42 0.00
## 6    cntry      gndr.c year.c -0.17 0.00
## 7 Residual        <NA>   <NA>  0.99 0.99
```

``` r
r2mlm(mod5,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.009458900
## slope variation 0.009985817
## mean variation  0.034937430
## sigma2          0.945617854
## 
## $R2s
##           total
## f   0.009458900
## v   0.009985817
## m   0.034937430
## fv  0.019444717
## fvm 0.054382146
```

``` r
anova(mod4,mod5)
```

```
## Data: diff_dat
## Models:
## mod4: sti.z ~ gndr.c + year.c + (gndr.c + year.c | cntry)
## mod5: sti.z ~ gndr.c + year.c + gndr.c:year.c + (gndr.c + year.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod4   10 1454816 1454927 -727398   1454796                         
## mod5   11 1454772 1454894 -727375   1454750 45.627  1   1.43e-11 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## mod6: random interaction between time and gender


``` r
mod6<-lmer(sti.z~gndr.c+year.c+
             gndr.c:year.c+(gndr.c+year.c+gndr.c:year.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: sti.z ~ gndr.c + year.c + gndr.c:year.c + (gndr.c + year.c +      gndr.c:year.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1454749.3 1454915.9 -727359.7 1454719.3    492328 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -5.0313 -0.6931 -0.0412  0.6423  5.1105 
## 
## Random effects:
##  Groups   Name          Variance  Std.Dev. Corr              
##  cntry    (Intercept)   3.570e-02 0.188944                   
##           gndr.c        5.330e-03 0.073007 -0.25             
##           year.c        2.412e-04 0.015531 -0.42 -0.26       
##           gndr.c:year.c 1.573e-05 0.003966  0.02 -0.78  0.28 
##  Residual               9.873e-01 0.993623                   
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##                 Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)    0.0218306  0.0329022 29.6949240   0.663 0.512129    
## gndr.c         0.2237597  0.0142781 27.0852294  15.672 4.17e-15 ***
## year.c         0.0047227  0.0026955 24.7365399   1.752 0.092145 .  
## gndr.c:year.c -0.0033743  0.0008604 25.8874225  -3.922 0.000577 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c year.c
## gndr.c      -0.215              
## year.c      -0.428 -0.222       
## gndr.c:yr.c  0.017 -0.794  0.220
```

``` r
getFE(mod6,round=3)
```

```
##                 Est.    SE     df      t     p     LL     UL
## (Intercept)    0.022 0.033 29.695  0.663 0.512 -0.045  0.089
## gndr.c         0.224 0.014 27.085 15.672 0.000  0.194  0.253
## year.c         0.005 0.003 24.737  1.752 0.092 -0.001  0.010
## gndr.c:year.c -0.003 0.001 25.887 -3.922 0.001 -0.005 -0.002
```

``` r
getVC(mod6)
```

```
##         grp          var1          var2 sdcor vcov
## 1     cntry   (Intercept)          <NA>  0.19 0.04
## 2     cntry        gndr.c          <NA>  0.07 0.01
## 3     cntry        year.c          <NA>  0.02 0.00
## 4     cntry gndr.c:year.c          <NA>  0.00 0.00
## 5     cntry   (Intercept)        gndr.c -0.25 0.00
## 6     cntry   (Intercept)        year.c -0.42 0.00
## 7     cntry   (Intercept) gndr.c:year.c  0.02 0.00
## 8     cntry        gndr.c        year.c -0.26 0.00
## 9     cntry        gndr.c gndr.c:year.c -0.78 0.00
## 10    cntry        year.c gndr.c:year.c  0.28 0.00
## 11 Residual          <NA>          <NA>  0.99 0.99
```

``` r
r2mlm(mod6,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.009574557
## slope variation 0.009967030
## mean variation  0.034688870
## sigma2          0.945769543
## 
## $R2s
##           total
## f   0.009574557
## v   0.009967030
## m   0.034688870
## fv  0.019541587
## fvm 0.054230457
```

``` r
anova(mod4,mod5,mod6)
```

```
## Data: diff_dat
## Models:
## mod4: sti.z ~ gndr.c + year.c + (gndr.c + year.c | cntry)
## mod5: sti.z ~ gndr.c + year.c + gndr.c:year.c + (gndr.c + year.c | cntry)
## mod6: sti.z ~ gndr.c + year.c + gndr.c:year.c + (gndr.c + year.c + gndr.c:year.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod4   10 1454816 1454927 -727398   1454796                         
## mod5   11 1454772 1454894 -727375   1454750 45.627  1  1.430e-11 ***
## mod6   15 1454749 1454916 -727360   1454719 30.658  4  3.594e-06 ***
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
##      21  0.0446 0.0526 25.1  -0.0638   0.1529   0.847  0.4051
##       0 -0.0900 0.0351 30.5  -0.1618  -0.0183  -2.563  0.0156
## 
## gndr.c =  0.5:
##  year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21  0.1975 0.0517 25.3   0.0910   0.3040   3.816  0.0008
##       0  0.1337 0.0321 27.7   0.0679   0.1996   4.162  0.0003
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_mod6,adjust="none",infer=c(T,T))
```

```
## gndr.c = -0.5:
##  contrast           estimate     SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0   0.1346 0.0553 22.9   0.0201    0.249   2.433  0.0232
## 
## gndr.c =  0.5:
##  contrast           estimate     SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0   0.0637 0.0593 24.6  -0.0584    0.186   1.076  0.2925
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
##    -0.5     21  0.0446 0.0526 25.1  -0.0638   0.1529   0.847  0.4051
##     0.5     21  0.1975 0.0517 25.3   0.0910   0.3040   3.816  0.0008
##    -0.5      0 -0.0900 0.0351 30.5  -0.1618  -0.0183  -2.563  0.0156
##     0.5      0  0.1337 0.0321 27.7   0.0679   0.1996   4.162  0.0003
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod6,adjust="none",infer=c(T,T))
```

```
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21    -0.1529 0.0110 30.7  -0.1753  -0.1305 -13.914 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)   0.1346 0.0553 22.9   0.0201   0.2491   2.433  0.0232
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     -0.0892 0.0585 24.6  -0.2097   0.0314  -1.525  0.1401
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)      0.2875 0.0561 24.6   0.1718   0.4032   5.122 <0.0001
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0         0.0637 0.0593 24.6  -0.0584   0.1859   1.076  0.2925
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0      -0.2238 0.0143 27.1  -0.2531  -0.1945 -15.672 <0.0001
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
##  diff_ESS11    0.153 0.0110 30.7    0.130    0.175  13.914 <0.0001
##  diff_ESS1     0.224 0.0143 27.1    0.194    0.253  15.672 <0.0001
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
##  diff_ESS11 - diff_ESS1  -0.0709 0.0181 25.9   -0.108  -0.0337  -3.922  0.0006
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
## [1] -0.07085938
## 
## $se
## [1] 0.01806796
## 
## $df
## [1] 25.88742
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
## [1] 7.147492
## 
## $p_low
## [1] 7.013297e-08
## 
## $t_high
## [1] -14.99114
## 
## $p_high
## [1] 1.414091e-14
## 
## $ci_level
## [1] 0.8
## 
## $ci_lower
## [1] -0.09462094
## 
## $ci_upper
## [1] -0.04709783
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

# texts for the figure

## obtain texts for the figure

time_slopes<-data.frame(pairs(change_mod6,adjust="none",infer=c(T,T)))
time_slope_text_men<-
  paste0("Linear change from 2002 to 2023 = ",round_tidy(time_slopes[2,"estimate"],2),
         ", 95% CI [",round_tidy(time_slopes[2,"lower.CL"],2),
         ", ", round_tidy(time_slopes[2,"upper.CL"],2),"], p = ",
         round_tidy(time_slopes[2,"p.value"],3))

time_slope_text_women<-
  paste0("Linear change from 2002 to 2023 = ",round_tidy(time_slopes[1,"estimate"],2),
         ", 95% CI [",round_tidy(time_slopes[1,"lower.CL"],2),
         ", ", round_tidy(time_slopes[1,"upper.CL"],2),"], p = ",
         round_tidy(time_slopes[1,"p.value"],3))

time_slope_texts<-data.frame(
  gndr.c=c("Women","Men"),
  year=c(2002,2002),
  yvar=c(-0.45,0.55),
  label=c(time_slope_text_women,time_slope_text_men)
)

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
  geom_text(
    data = time_slope_texts,
    aes(x = year, y = yvar, color = gndr.c,label=label),
    #inherit.aes = FALSE,
    hjust = 0,
    size = 4
  ) +
  scale_color_manual(values = my_colors) +
  scale_y_continuous(limits = c(-1.1, 1.1)) +
  scale_x_continuous(limits = c(2001, 2024),
                     breaks = c(seq(2002, 2020, 2),2023))+
  ylab("Mean-level of value male-typicality")+
  theme(legend.title=element_blank(),
        legend.text  = element_text(size = 16),
        axis.title.x = element_text(size = 18, face = "bold"),   
        axis.title.y = element_text(size = 18, face = "bold"),   
        axis.text.x  = element_text(size = 14),                  
        axis.text.y  = element_text(size = 14))

p_time_trends
```

![](Analysis_stimulation_files/figure-html/unnamed-chunk-48-1.png)<!-- -->

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
pred_cntry_dat$sti.z_mean<-predict(mod6,newdata=pred_cntry_dat)

pred_cntry_dat$year=pred_cntry_dat$year.c+2002

pred_cntry_dat$gender<-
  case_when(
    pred_cntry_dat$gndr.c==0.5~"men",
    pred_cntry_dat$gndr.c==-0.5~"women",
  )

range(pred_cntry_dat$sti.z_mean)
```

```
## [1] -0.5043686  0.5131115
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

pdf("../results/sti/country_specific_time_trends.pdf", width = 8, height = 6)

for (ctry in countries) {
  print(
    ggplot(pred_cntry_dat[pred_cntry_dat$cntry == ctry, ],
       aes(x = year, y = sti.z_mean, color = gender)) +
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
    name   = "Mean-level of value stimulation",
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
         aes(x = year, y = sti.z_mean, color = gender)) +
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
    name   = "Mean-level of value stimulation",
    sec.axis = sec_axis(~ scale_y_to_gei(.),
                        name = "Gender Equality Index (GEI)")
  ) +
  scale_x_continuous(limits = c(2001, 2024),
                     breaks = c(seq(2002, 2020, 2),2023)) +
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

```
## Warning: Removed 12 rows containing missing values or values outside the scale range (`geom_line()`).
```

![](Analysis_stimulation_files/figure-html/unnamed-chunk-50-1.png)<!-- -->

``` r
png(filename = 
      "../results/sti/country_time_trend_facets.png",
    units = "cm",
    width = 21.0,height=29.7*(3/4),res = 600)
facet_plot
```

```
## Warning: Removed 12 rows containing missing values or values outside the scale range (`geom_line()`).
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
## [1] 15.33444
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
## 1    0.19               -0.03                    0.02                   -0.02                     -0.04
## 2    0.16                0.12                   -0.02                    0.11                      0.13
## 3    0.21                0.34                    0.03                    0.36                      0.32
## 4    0.22                0.16                   -0.04                    0.14                      0.17
## 5    0.22                0.08                   -0.04                    0.06                      0.10
## 6    0.31                0.35                   -0.14                    0.28                      0.42
## 7    0.24                0.05                   -0.07                    0.02                      0.09
## 8    0.19                0.25                   -0.05                    0.22                      0.27
## 9    0.21                0.01                   -0.05                   -0.02                      0.04
## 10   0.18               -0.07                   -0.16                   -0.14                      0.01
## 11   0.17                0.04                   -0.06                    0.01                      0.07
## 12   0.22               -0.08                   -0.04                   -0.10                     -0.06
## 13   0.21               -0.02                   -0.12                   -0.08                      0.04
## 14   0.27               -0.04                   -0.14                   -0.11                      0.03
## 15   0.33               -0.08                   -0.10                   -0.13                     -0.03
## 16   0.19                0.21                   -0.05                    0.19                      0.24
## 17   0.11                0.01                    0.04                    0.03                     -0.01
## 18   0.16               -0.01                    0.00                   -0.01                     -0.01
## 19   0.19                0.08                   -0.10                    0.03                      0.14
## 20   0.18                0.11                   -0.03                    0.10                      0.13
## 21   0.40                0.15                   -0.16                    0.07                      0.23
## 22   0.24               -0.44                   -0.12                   -0.50                     -0.38
## 23   0.25               -0.51                   -0.15                   -0.59                     -0.44
## 24   0.18                0.02                   -0.02                    0.01                      0.03
## 25   0.23                0.10                   -0.06                    0.07                      0.13
## 26   0.29               -0.17                   -0.09                   -0.21                     -0.12
## 27   0.34               -0.09                   -0.22                   -0.20                      0.02
## 28   0.22                0.13                   -0.05                    0.11                      0.16
## 29   0.35                0.17                   -0.23                    0.06                      0.29
## 30   0.17                0.23                   -0.04                    0.21                      0.26
## 31   0.19                0.09                   -0.05                    0.06                      0.11
## 32   0.27                0.21                   -0.10                    0.16                      0.27
## 33   0.11                1.62                    0.02                    1.63                      1.61
## 34   0.18                0.36                   -0.01                    0.36                      0.37
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
## 1     PL               -0.17
## 2     PT               -0.09
## 3     FR               -0.08
## 4     HR               -0.08
## 5     ES               -0.07
## 6     GR               -0.04
## 7     AT               -0.03
## 8     GB               -0.02
## 9     IL               -0.01
## 10    EE                0.01
## 11    IE                0.01
## 12    NL                0.02
## 13    FI                0.04
## 14    DE                0.05
## 15    CY                0.08
## 16    IS                0.08
## 17    SI                0.09
## 18    NO                0.10
## 19    IT                0.11
## 20    BE                0.12
## 21    LT                0.15
## 22    CH                0.16
## 23    RU                0.17
## 24    HU                0.21
## 25    SK                0.21
## 26    SE                0.23
## 27    DK                0.25
## 28    BG                0.34
## 29    CZ                0.35
## 30    UA                0.36
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
## 1     RU                   -0.23
## 2     PT                   -0.22
## 3     ES                   -0.16
## 4     LT                   -0.16
## 5     CZ                   -0.14
## 6     GR                   -0.14
## 7     GB                   -0.12
## 8     HR                   -0.10
## 9     IS                   -0.10
## 10    SK                   -0.10
## 11    PL                   -0.09
## 12    DE                   -0.07
## 13    FI                   -0.06
## 14    NO                   -0.06
## 15    DK                   -0.05
## 16    EE                   -0.05
## 17    HU                   -0.05
## 18    SI                   -0.05
## 19    CH                   -0.04
## 20    CY                   -0.04
## 21    FR                   -0.04
## 22    SE                   -0.04
## 23    IT                   -0.03
## 24    BE                   -0.02
## 25    NL                   -0.02
## 26    UA                   -0.01
## 27    IL                    0.00
## 28    AT                    0.02
## 29    BG                    0.03
## 30    IE                    0.04
```


#### Country-level predictors: GEI


``` r
mod6_GEI<-lmer(sti.z~gndr.c+year.c+
             gndr.c:year.c+
               gei.z.cm:gndr.c+gei.z.cm:year.c+gei.z.cm:gndr.c:year.c+
               (gndr.c+year.c+gndr.c:year.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_GEI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: sti.z ~ gndr.c + year.c + gndr.c:year.c + gei.z.cm:gndr.c + gei.z.cm:year.c +  
##     gei.z.cm:gndr.c:year.c + (gndr.c + year.c + gndr.c:year.c |      cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1415865.6 1416065.0 -707914.8 1415829.6    480346 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -5.0396 -0.6941 -0.0399  0.6443  5.1129 
## 
## Random effects:
##  Groups   Name          Variance  Std.Dev. Corr              
##  cntry    (Intercept)   3.158e-02 0.177695                   
##           gndr.c        5.122e-03 0.071567 -0.37             
##           year.c        1.992e-04 0.014115 -0.47 -0.36       
##           gndr.c:year.c 1.484e-05 0.003852  0.13 -0.76  0.36 
##  Residual               9.842e-01 0.992073                   
## Number of obs: 480364, groups:  cntry, 33
## 
## Fixed effects:
##                          Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)             0.0324728  0.0314694 28.2465885   1.032  0.31088    
## gndr.c                  0.2280124  0.0141724 28.1849041  16.088 9.78e-16 ***
## year.c                  0.0044560  0.0024930 25.6196073   1.787  0.08571 .  
## gndr.c:year.c          -0.0036306  0.0008618 26.4843256  -4.213  0.00026 ***
## gndr.c:gei.z.cm        -0.0149859  0.0136409 26.6886861  -1.099  0.28176    
## year.c:gei.z.cm        -0.0080235  0.0022834 29.5562690  -3.514  0.00144 ** 
## gndr.c:year.c:gei.z.cm  0.0005954  0.0009866 33.2773693   0.604  0.55026    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c year.c gnd.:. gn.:.. yr.:..
## gndr.c      -0.326                                   
## year.c      -0.486 -0.309                            
## gndr.c:yr.c  0.102 -0.775  0.274                     
## gndr.c:g.z. -0.003 -0.067  0.005  0.090              
## yr.c:g.z.cm  0.010 -0.002 -0.025  0.001 -0.557       
## gndr.c:.:..  0.000  0.083 -0.002 -0.160 -0.751  0.324
```

``` r
getFE(mod6_GEI,round=2,p.round = 5)
```

```
##                         Est.   SE    df     t       p    LL   UL
## (Intercept)             0.03 0.03 28.25  1.03 0.31088 -0.03 0.10
## gndr.c                  0.23 0.01 28.18 16.09 0.00000  0.20 0.26
## year.c                  0.00 0.00 25.62  1.79 0.08571  0.00 0.01
## gndr.c:year.c           0.00 0.00 26.48 -4.21 0.00026 -0.01 0.00
## gndr.c:gei.z.cm        -0.01 0.01 26.69 -1.10 0.28176 -0.04 0.01
## year.c:gei.z.cm        -0.01 0.00 29.56 -3.51 0.00144 -0.01 0.00
## gndr.c:year.c:gei.z.cm  0.00 0.00 33.28  0.60 0.55026  0.00 0.00
```

``` r
getVC(mod6_GEI)
```

```
##         grp          var1          var2 sdcor vcov
## 1     cntry   (Intercept)          <NA>  0.18 0.03
## 2     cntry        gndr.c          <NA>  0.07 0.01
## 3     cntry        year.c          <NA>  0.01 0.00
## 4     cntry gndr.c:year.c          <NA>  0.00 0.00
## 5     cntry   (Intercept)        gndr.c -0.37 0.00
## 6     cntry   (Intercept)        year.c -0.47 0.00
## 7     cntry   (Intercept) gndr.c:year.c  0.13 0.00
## 8     cntry        gndr.c        year.c -0.36 0.00
## 9     cntry        gndr.c gndr.c:year.c -0.76 0.00
## 10    cntry        year.c gndr.c:year.c  0.36 0.00
## 11 Residual          <NA>          <NA>  0.99 0.98
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_GEI,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] 17.4024
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_GEI,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] 5.690747
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
##      21  0.2945 0.0672 29.6   0.1571   0.4319   4.380  0.0001
##       0  0.0325 0.0315 28.2  -0.0320   0.0969   1.032  0.3109
## 
## gei.z.cm =  0:
##  year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21  0.1260 0.0461 27.3   0.0314   0.2207   2.732  0.0109
##       0  0.0325 0.0315 28.2  -0.0320   0.0969   1.032  0.3109
## 
## gei.z.cm =  1:
##  year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21 -0.0424 0.0658 27.6  -0.1774   0.0925  -0.645  0.5244
##       0  0.0325 0.0315 28.2  -0.0320   0.0969   1.032  0.3109
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_mod6_GEI,adjust="none",infer=c(T,T))
```

```
## gei.z.cm = -1:
##  contrast           estimate     SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0   0.2621 0.0719 32.0   0.1157   0.4084   3.647  0.0009
## 
## gei.z.cm =  0:
##  contrast           estimate     SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0   0.0936 0.0524 25.6  -0.0141   0.2013   1.787  0.0857
## 
## gei.z.cm =  1:
##  contrast           estimate     SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0  -0.0749 0.0701 29.0  -0.2183   0.0685  -1.068  0.2941
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
##    -0.5     21  0.2174 0.0687 28.5   0.0769  0.35796   3.166  0.0037
##     0.5     21  0.3717 0.0672 29.0   0.2342  0.50909   5.531 <0.0001
##    -0.5      0 -0.0890 0.0352 31.9  -0.1607 -0.01735  -2.530  0.0165
##     0.5      0  0.1540 0.0308 29.1   0.0910  0.21698   4.997 <0.0001
## 
## gei.z.cm =  0:
##  gndr.c year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21  0.0502 0.0471 26.3  -0.0466  0.14693   1.065  0.2966
##     0.5     21  0.2019 0.0459 26.4   0.1077  0.29613   4.403  0.0002
##    -0.5      0 -0.0815 0.0344 29.3  -0.1519 -0.01113  -2.367  0.0247
##     0.5      0  0.1465 0.0299 26.2   0.0850  0.20795   4.896 <0.0001
## 
## gei.z.cm =  1:
##  gndr.c year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21 -0.1171 0.0671 26.5  -0.2549  0.02072  -1.745  0.0926
##     0.5     21  0.0322 0.0656 26.7  -0.1024  0.16681   0.491  0.6274
##    -0.5      0 -0.0740 0.0350 31.3  -0.1455 -0.00261  -2.113  0.0427
##     0.5      0  0.1390 0.0306 28.1   0.0764  0.20157   4.549 <0.0001
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
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21    -0.1543 0.0194 34.0  -0.1936  -0.1149  -7.969 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)   0.3064 0.0693 29.9   0.1649   0.4480   4.423  0.0001
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0      0.0634 0.0771 30.6  -0.0939   0.2208   0.823  0.4171
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)      0.4607 0.0688 32.8   0.3208   0.6006   6.700 <0.0001
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0         0.2177 0.0772 31.3   0.0602   0.3752   2.819  0.0083
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0      -0.2430 0.0203 32.6  -0.2844  -0.2016 -11.957 <0.0001
## 
## gei.z.cm =  0:
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21    -0.1518 0.0114 29.6  -0.1751  -0.1284 -13.269 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)   0.1317 0.0506 23.8   0.0272   0.2362   2.601  0.0157
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     -0.0963 0.0550 25.5  -0.2095   0.0168  -1.751  0.0919
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)      0.2835 0.0512 25.4   0.1780   0.3889   5.532 <0.0001
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0         0.0555 0.0555 25.4  -0.0588   0.1697   0.999  0.3273
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0      -0.2280 0.0142 28.2  -0.2570  -0.1990 -16.088 <0.0001
## 
## gei.z.cm =  1:
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21    -0.1493 0.0164 28.4  -0.1829  -0.1157  -9.104 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)  -0.0430 0.0671 26.5  -0.1808   0.0947  -0.642  0.5264
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     -0.2561 0.0753 27.9  -0.4104  -0.1017  -3.399  0.0021
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)      0.1062 0.0669 29.3  -0.0304   0.2429   1.589  0.1227
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0        -0.1068 0.0752 28.2  -0.2608   0.0472  -1.420  0.1666
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0      -0.2130 0.0190 25.8  -0.2521  -0.1740 -11.214 <0.0001
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
##  diff_ESS11    0.154 0.0194 34.0    0.115    0.194   7.969 <0.0001
##  diff_ESS1     0.243 0.0203 32.6    0.202    0.284  11.957 <0.0001
## 
## gei.z.cm =  0:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11    0.152 0.0114 29.6    0.128    0.175  13.269 <0.0001
##  diff_ESS1     0.228 0.0142 28.2    0.199    0.257  16.088 <0.0001
## 
## gei.z.cm =  1:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11    0.149 0.0164 28.4    0.116    0.183   9.104 <0.0001
##  diff_ESS1     0.213 0.0190 25.8    0.174    0.252  11.214 <0.0001
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
##  diff_ESS11 - diff_ESS1  -0.0887 0.0296 36.6   -0.149  -0.0287  -2.997  0.0049
## 
## gei.z.cm =  0:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1  -0.0762 0.0181 26.5   -0.113  -0.0391  -4.213  0.0003
## 
## gei.z.cm =  1:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1  -0.0637 0.0252 23.9   -0.116  -0.0117  -2.527  0.0186
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```


#### Country-level predictors: GGGI


``` r
mod6_GGGI<-lmer(sti.z~gndr.c+year.c+
             gndr.c:year.c+
               gggi.z.cm:gndr.c+gggi.z.cm:year.c+gggi.z.cm:gndr.c:year.c+
               (gndr.c+year.c+gndr.c:year.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_GGGI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: sti.z ~ gndr.c + year.c + gndr.c:year.c + gggi.z.cm:gndr.c +  
##     gggi.z.cm:year.c + gggi.z.cm:gndr.c:year.c + (gndr.c + year.c +      gndr.c:year.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1075893.1 1076087.6 -537928.6 1075857.1    363834 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -5.2239 -0.6927 -0.0350  0.6466  5.0989 
## 
## Random effects:
##  Groups   Name          Variance  Std.Dev. Corr              
##  cntry    (Intercept)   0.0862456 0.293676                   
##           gndr.c        0.0050257 0.070893 -0.35             
##           year.c        0.0002155 0.014680 -0.85  0.17       
##           gndr.c:year.c 0.0000201 0.004484 -0.08 -0.74  0.23 
##  Residual               0.9877714 0.993867                   
## Number of obs: 363852, groups:  cntry, 34
## 
## Fixed effects:
##                           Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)              0.0811017  0.0515163 18.0507458   1.574  0.13278    
## gndr.c                   0.2197624  0.0150708 25.8967548  14.582 5.35e-14 ***
## year.c                   0.0003310  0.0026075 16.6163758   0.127  0.90050    
## gndr.c:year.c           -0.0030997  0.0010257 29.7855915  -3.022  0.00512 ** 
## gndr.c:gggi.z.cm        -0.0090701  0.0156086 31.9366676  -0.581  0.56525    
## year.c:gggi.z.cm        -0.0006762  0.0015148 32.6191225  -0.446  0.65826    
## gndr.c:year.c:gggi.z.cm -0.0001409  0.0011372 34.6961735  -0.124  0.90213    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c year.c gnd.:. gn.:.. yr.:..
## gndr.c      -0.285                                   
## year.c      -0.848  0.138                            
## gndr.c:yr.c -0.051 -0.796  0.165                     
## gndr.c:gg.. -0.002 -0.021  0.008  0.030              
## yr.c:ggg.z.  0.018  0.007 -0.055 -0.017 -0.201       
## gndr.c:.:..  0.001  0.027 -0.009 -0.055 -0.850  0.222
```

``` r
getFE(mod6_GGGI,round=2,p.round = 5)
```

```
##                          Est.   SE    df     t       p    LL   UL
## (Intercept)              0.08 0.05 18.05  1.57 0.13278 -0.03 0.19
## gndr.c                   0.22 0.02 25.90 14.58 0.00000  0.19 0.25
## year.c                   0.00 0.00 16.62  0.13 0.90050 -0.01 0.01
## gndr.c:year.c            0.00 0.00 29.79 -3.02 0.00512 -0.01 0.00
## gndr.c:gggi.z.cm        -0.01 0.02 31.94 -0.58 0.56525 -0.04 0.02
## year.c:gggi.z.cm         0.00 0.00 32.62 -0.45 0.65826  0.00 0.00
## gndr.c:year.c:gggi.z.cm  0.00 0.00 34.70 -0.12 0.90213  0.00 0.00
```

``` r
getVC(mod6_GGGI)
```

```
##         grp          var1          var2 sdcor  vcov
## 1     cntry   (Intercept)          <NA>  0.29  0.09
## 2     cntry        gndr.c          <NA>  0.07  0.01
## 3     cntry        year.c          <NA>  0.01  0.00
## 4     cntry gndr.c:year.c          <NA>  0.00  0.00
## 5     cntry   (Intercept)        gndr.c -0.35 -0.01
## 6     cntry   (Intercept)        year.c -0.85  0.00
## 7     cntry   (Intercept) gndr.c:year.c -0.08  0.00
## 8     cntry        gndr.c        year.c  0.17  0.00
## 9     cntry        gndr.c gndr.c:year.c -0.74  0.00
## 10    cntry        year.c gndr.c:year.c  0.23  0.00
## 11 Residual          <NA>          <NA>  0.99  0.99
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_GGGI,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] 10.66801
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_GGGI,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] -27.79134
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
##  year.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21 0.1023 0.0449 33.2   0.0110    0.194   2.279  0.0292
##       0 0.0811 0.0515 18.1  -0.0271    0.189   1.574  0.1328
## 
## gggi.z.cm =  0:
##  year.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21 0.0881 0.0295 32.1   0.0280    0.148   2.988  0.0053
##       0 0.0811 0.0515 18.1  -0.0271    0.189   1.574  0.1328
## 
## gggi.z.cm =  1:
##  year.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21 0.0739 0.0418 32.6  -0.0112    0.159   1.767  0.0866
##       0 0.0811 0.0515 18.1  -0.0271    0.189   1.574  0.1328
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_mod6_GGGI,adjust="none",infer=c(T,T))
```

```
## gggi.z.cm = -1:
##  contrast           estimate     SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0  0.02115 0.0648 26.3   -0.112    0.154   0.326  0.7468
## 
## gggi.z.cm =  0:
##  contrast           estimate     SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0  0.00695 0.0548 16.6   -0.109    0.123   0.127  0.9005
## 
## gggi.z.cm =  1:
##  contrast           estimate     SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0 -0.00725 0.0618 23.9   -0.135    0.120  -0.117  0.9076
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
##  gndr.c year.c   emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21  0.01891 0.0447 33.6  -0.0719   0.1097   0.423  0.6749
##     0.5     21  0.18560 0.0471 32.9   0.0897   0.2815   3.939  0.0004
##    -0.5      0 -0.03331 0.0547 18.7  -0.1480   0.0813  -0.609  0.5499
##     0.5      0  0.19552 0.0505 18.7   0.0896   0.3014   3.869  0.0011
## 
## gggi.z.cm =  0:
##  gndr.c year.c   emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21  0.01072 0.0297 31.9  -0.0498   0.0712   0.361  0.7205
##     0.5     21  0.16539 0.0307 32.2   0.1029   0.2279   5.388 <0.0001
##    -0.5      0 -0.02878 0.0541 17.9  -0.1426   0.0850  -0.531  0.6016
##     0.5      0  0.19098 0.0499 17.8   0.0861   0.2959   3.828  0.0013
## 
## gggi.z.cm =  1:
##  gndr.c year.c   emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21  0.00253 0.0418 33.1  -0.0824   0.0875   0.061  0.9520
##     0.5     21  0.14517 0.0438 32.3   0.0561   0.2343   3.317  0.0023
##    -0.5      0 -0.02424 0.0547 18.6  -0.1389   0.0904  -0.443  0.6627
##     0.5      0  0.18645 0.0505 18.6   0.0807   0.2922   3.695  0.0016
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod6_GGGI,adjust="none",infer=c(T,T))
```

```
## gggi.z.cm = -1:
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21   -0.16670 0.0195 34.6  -0.2063  -0.1271  -8.546 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)  0.05222 0.0639 22.9  -0.0801   0.1845   0.817  0.4224
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0    -0.17661 0.0630 26.7  -0.3059  -0.0473  -2.805  0.0093
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)     0.21892 0.0689 26.0   0.0772   0.3606   3.175  0.0038
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0       -0.00992 0.0697 30.5  -0.1522   0.1324  -0.142  0.8878
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0     -0.22883 0.0219 33.2  -0.2734  -0.1842 -10.439 <0.0001
## 
## gggi.z.cm =  0:
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21   -0.15467 0.0132 31.6  -0.1816  -0.1278 -11.714 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)  0.03950 0.0540 14.7  -0.0759   0.1549   0.731  0.4762
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0    -0.18026 0.0516 16.3  -0.2895  -0.0711  -3.494  0.0029
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)     0.19417 0.0592 16.9   0.0692   0.3191   3.280  0.0044
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0       -0.02560 0.0575 18.7  -0.1461   0.0949  -0.445  0.6615
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0     -0.21976 0.0151 25.9  -0.2507  -0.1888 -14.582 <0.0001
## 
## gggi.z.cm =  1:
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21   -0.14264 0.0181 34.9  -0.1794  -0.1059  -7.876 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)  0.02678 0.0610 20.8  -0.1002   0.1538   0.439  0.6653
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0    -0.18391 0.0598 24.3  -0.3074  -0.0605  -3.073  0.0052
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)     0.16942 0.0660 23.6   0.0330   0.3058   2.566  0.0171
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0       -0.04127 0.0664 27.8  -0.1773   0.0947  -0.622  0.5390
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0     -0.21069 0.0215 29.1  -0.2546  -0.1668  -9.813 <0.0001
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
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11    0.167 0.0195 34.6    0.127    0.206   8.546 <0.0001
##  diff_ESS1     0.229 0.0219 33.2    0.184    0.273  10.439 <0.0001
## 
## gggi.z.cm =  0:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11    0.155 0.0132 31.6    0.128    0.182  11.714 <0.0001
##  diff_ESS1     0.220 0.0151 25.9    0.189    0.251  14.582 <0.0001
## 
## gggi.z.cm =  1:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11    0.143 0.0181 34.9    0.106    0.179   7.876 <0.0001
##  diff_ESS1     0.211 0.0215 29.1    0.167    0.255   9.813 <0.0001
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
##  diff_ESS11 - diff_ESS1  -0.0621 0.0330 34.3   -0.129  0.00496  -1.882  0.0684
## 
## gggi.z.cm =  0:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1  -0.0651 0.0215 29.8   -0.109 -0.02109  -3.022  0.0051
## 
## gggi.z.cm =  1:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1  -0.0681 0.0313 30.5   -0.132 -0.00423  -2.176  0.0374
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

#### Country-level predictors: GDI


``` r
mod6_GDI<-lmer(sti.z~gndr.c+year.c+
             gndr.c:year.c+
               gdi.z.cm:gndr.c+gdi.z.cm:year.c+gdi.z.cm:gndr.c:year.c+
               (gndr.c+year.c+gndr.c:year.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_GDI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: sti.z ~ gndr.c + year.c + gndr.c:year.c + gdi.z.cm:gndr.c + gdi.z.cm:year.c +  
##     gdi.z.cm:gndr.c:year.c + (gndr.c + year.c + gndr.c:year.c |      cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1454740.6 1454940.5 -727352.3 1454704.6    492325 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -5.0303 -0.6930 -0.0409  0.6425  5.1140 
## 
## Random effects:
##  Groups   Name          Variance  Std.Dev. Corr              
##  cntry    (Intercept)   3.568e-02 0.18888                    
##           gndr.c        4.209e-03 0.06488  -0.22             
##           year.c        2.015e-04 0.01419  -0.53 -0.04       
##           gndr.c:year.c 1.398e-05 0.00374  -0.02 -0.74  0.12 
##  Residual               9.873e-01 0.99362                    
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##                          Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)             0.0213660  0.0328948 29.6720039   0.650 0.520995    
## gndr.c                  0.2228769  0.0130450 28.1475715  17.085  < 2e-16 ***
## year.c                  0.0047916  0.0024694 25.8885928   1.940 0.063295 .  
## gndr.c:year.c          -0.0033463  0.0008291 27.6663734  -4.036 0.000388 ***
## gndr.c:gdi.z.cm         0.0337978  0.0135247 32.3697608   2.499 0.017716 *  
## year.c:gdi.z.cm        -0.0081415  0.0021726 29.6496473  -3.747 0.000771 ***
## gndr.c:year.c:gdi.z.cm -0.0014569  0.0009721 39.1509726  -1.499 0.141966    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c year.c gnd.:. gn.:.. yr.:..
## gndr.c      -0.186                                   
## year.c      -0.542 -0.034                            
## gndr.c:yr.c -0.013 -0.776  0.095                     
## gndr.c:gd..  0.001 -0.015  0.000  0.011              
## yr.c:gd.z.c  0.008 -0.001 -0.022  0.001 -0.158       
## gndr.c:.:.. -0.004  0.008  0.004 -0.028 -0.775  0.092
```

``` r
getFE(mod6_GDI,round=2,p.round = 5)
```

```
##                         Est.   SE    df     t       p    LL   UL
## (Intercept)             0.02 0.03 29.67  0.65 0.52099 -0.05 0.09
## gndr.c                  0.22 0.01 28.15 17.09 0.00000  0.20 0.25
## year.c                  0.00 0.00 25.89  1.94 0.06329  0.00 0.01
## gndr.c:year.c           0.00 0.00 27.67 -4.04 0.00039 -0.01 0.00
## gndr.c:gdi.z.cm         0.03 0.01 32.37  2.50 0.01772  0.01 0.06
## year.c:gdi.z.cm        -0.01 0.00 29.65 -3.75 0.00077 -0.01 0.00
## gndr.c:year.c:gdi.z.cm  0.00 0.00 39.15 -1.50 0.14197  0.00 0.00
```

``` r
getVC(mod6_GDI)
```

```
##         grp          var1          var2 sdcor vcov
## 1     cntry   (Intercept)          <NA>  0.19 0.04
## 2     cntry        gndr.c          <NA>  0.06 0.00
## 3     cntry        year.c          <NA>  0.01 0.00
## 4     cntry gndr.c:year.c          <NA>  0.00 0.00
## 5     cntry   (Intercept)        gndr.c -0.22 0.00
## 6     cntry   (Intercept)        year.c -0.53 0.00
## 7     cntry   (Intercept) gndr.c:year.c -0.02 0.00
## 8     cntry        gndr.c        year.c -0.04 0.00
## 9     cntry        gndr.c gndr.c:year.c -0.74 0.00
## 10    cntry        year.c gndr.c:year.c  0.12 0.00
## 11 Residual          <NA>          <NA>  0.99 0.99
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_GDI,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] 16.48292
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_GDI,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] 11.09861
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
##      21  0.2930 0.0639 29.6   0.1624   0.4235   4.586 <0.0001
##       0  0.0214 0.0329 29.7  -0.0458   0.0886   0.650  0.5210
## 
## gdi.z.cm =  0:
##  year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21  0.1220 0.0438 27.4   0.0321   0.2119   2.783  0.0096
##       0  0.0214 0.0329 29.7  -0.0458   0.0886   0.650  0.5210
## 
## gdi.z.cm =  1:
##  year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21 -0.0490 0.0627 27.7  -0.1774   0.0794  -0.782  0.4409
##       0  0.0214 0.0329 29.7  -0.0458   0.0886   0.650  0.5210
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_mod6_GDI,adjust="none",infer=c(T,T))
```

```
## gdi.z.cm = -1:
##  contrast           estimate     SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0   0.2716 0.0698 32.8  0.12953   0.4137   3.890  0.0005
## 
## gdi.z.cm =  0:
##  contrast           estimate     SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0   0.1006 0.0519 25.9 -0.00599   0.2072   1.940  0.0633
## 
## gdi.z.cm =  1:
##  contrast           estimate     SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0  -0.0703 0.0683 30.0 -0.20988   0.0692  -1.030  0.3114
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
##    -0.5     21  0.2183 0.0648 28.2   0.0856  0.35090   3.370  0.0022
##     0.5     21  0.3677 0.0642 28.9   0.2363  0.49898   5.727 <0.0001
##    -0.5      0 -0.0732 0.0354 32.7  -0.1452 -0.00116  -2.068  0.0466
##     0.5      0  0.1159 0.0330 30.5   0.0485  0.18332   3.509  0.0014
## 
## gdi.z.cm =  0:
##  gndr.c year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21  0.0457 0.0445 26.2  -0.0456  0.13702   1.028  0.3134
##     0.5     21  0.1983 0.0439 26.6   0.1082  0.28843   4.517  0.0001
##    -0.5      0 -0.0901 0.0347 30.3  -0.1609 -0.01922  -2.595  0.0144
##     0.5      0  0.1328 0.0323 28.2   0.0666  0.19899   4.109  0.0003
## 
## gdi.z.cm =  1:
##  gndr.c year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21 -0.1269 0.0635 26.5  -0.2572  0.00347  -1.999  0.0560
##     0.5     21  0.0289 0.0629 27.1  -0.1002  0.15801   0.460  0.6494
##    -0.5      0 -0.1070 0.0353 32.5  -0.1789 -0.03504  -3.027  0.0048
##     0.5      0  0.1497 0.0330 30.4   0.0823  0.21707   4.535 <0.0001
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
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21    -0.1494 0.0176 31.5 -0.18518  -0.1136  -8.511 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)   0.2914 0.0699 30.1  0.14862   0.4343   4.167  0.0002
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0      0.1024 0.0714 31.6 -0.04309   0.2478   1.434  0.1614
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)      0.4408 0.0704 33.1  0.29772   0.5840   6.266 <0.0001
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0         0.2518 0.0723 32.4  0.10458   0.3989   3.483  0.0014
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0      -0.1891 0.0189 29.4 -0.22777  -0.1504  -9.988 <0.0001
## 
## gdi.z.cm =  0:
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21    -0.1526 0.0110 30.8 -0.17502  -0.1302 -13.889 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)   0.1358 0.0518 24.3  0.02900   0.2425   2.623  0.0148
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     -0.0871 0.0521 25.5 -0.19441   0.0202  -1.671  0.1070
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)      0.2884 0.0529 25.9  0.17961   0.3971   5.451 <0.0001
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0         0.0655 0.0534 25.8 -0.04431   0.1753   1.226  0.2311
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0      -0.2229 0.0130 28.1 -0.24959  -0.1962 -17.085 <0.0001
## 
## gdi.z.cm =  1:
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21    -0.1558 0.0166 36.2 -0.18954  -0.1221  -9.365 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)  -0.0199 0.0683 27.6 -0.15995   0.1201  -0.292  0.7728
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     -0.2766 0.0698 29.0 -0.41934  -0.1338  -3.963  0.0004
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)      0.1359 0.0688 30.3 -0.00464   0.2764   1.974  0.0575
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0        -0.1208 0.0708 30.0 -0.26547   0.0239  -1.705  0.0986
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0      -0.2567 0.0187 32.3 -0.29465  -0.2187 -13.763 <0.0001
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
##  diff_ESS11    0.149 0.0176 31.5    0.114    0.185   8.511 <0.0001
##  diff_ESS1     0.189 0.0189 29.4    0.150    0.228   9.988 <0.0001
## 
## gdi.z.cm =  0:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11    0.153 0.0110 30.8    0.130    0.175  13.889 <0.0001
##  diff_ESS1     0.223 0.0130 28.1    0.196    0.250  17.085 <0.0001
## 
## gdi.z.cm =  1:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11    0.156 0.0166 36.2    0.122    0.190   9.365 <0.0001
##  diff_ESS1     0.257 0.0187 32.3    0.219    0.295  13.763 <0.0001
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
##  diff_ESS11 - diff_ESS1  -0.0397 0.0272 32.1  -0.0951   0.0157  -1.459  0.1543
## 
## gdi.z.cm =  0:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1  -0.0703 0.0174 27.7  -0.1060  -0.0346  -4.036  0.0004
## 
## gdi.z.cm =  1:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1  -0.1009 0.0265 35.6  -0.1545  -0.0472  -3.813  0.0005
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

#### Country-level predictors: log_GDP


``` r
mod6_log_GDP<-lmer(sti.z~gndr.c+year.c+
             gndr.c:year.c+
               log_gdp.z.cm:gndr.c+log_gdp.z.cm:year.c+log_gdp.z.cm:gndr.c:year.c+
               (gndr.c+year.c+gndr.c:year.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_log_GDP)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: sti.z ~ gndr.c + year.c + gndr.c:year.c + log_gdp.z.cm:gndr.c +  
##     log_gdp.z.cm:year.c + log_gdp.z.cm:gndr.c:year.c + (gndr.c +      year.c + gndr.c:year.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1454748.6 1454948.5 -727356.3 1454712.6    492325 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -5.0308 -0.6928 -0.0404  0.6423  5.1077 
## 
## Random effects:
##  Groups   Name          Variance  Std.Dev. Corr              
##  cntry    (Intercept)   0.0359731 0.18967                    
##           gndr.c        0.0049151 0.07011  -0.16             
##           year.c        0.0002377 0.01542  -0.42 -0.41       
##           gndr.c:year.c 0.0000156 0.00395  -0.07 -0.78  0.42 
##  Residual               0.9872840 0.99362                    
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##                              Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)                 0.0220594  0.0330277 29.5356341   0.668 0.509379    
## gndr.c                      0.2278838  0.0138734 26.7999237  16.426 1.63e-15 ***
## year.c                      0.0046659  0.0026760 24.8220734   1.744 0.093599 .  
## gndr.c:year.c              -0.0035850  0.0008612 24.3911948  -4.163 0.000339 ***
## gndr.c:log_gdp.z.cm        -0.0285637  0.0144422 25.8607554  -1.978 0.058700 .  
## year.c:log_gdp.z.cm        -0.0019882  0.0024326 25.5478576  -0.817 0.421306    
## gndr.c:year.c:log_gdp.z.cm  0.0012126  0.0008936 26.2515153   1.357 0.186353    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c year.c gnd.:. g.:_.. y.:_..
## gndr.c      -0.132                                   
## year.c      -0.429 -0.350                            
## gndr.c:yr.c -0.059 -0.796  0.328                     
## gndr.c:l_.. -0.005 -0.101 -0.006  0.100              
## yr.c:lg_g..  0.011 -0.010  0.008  0.007 -0.433       
## gndr.:.:_.. -0.001  0.099  0.007 -0.113 -0.820  0.322
```

``` r
getFE(mod6_log_GDP,round=2,p.round = 5)
```

```
##                             Est.   SE    df     t       p    LL   UL
## (Intercept)                 0.02 0.03 29.54  0.67 0.50938 -0.05 0.09
## gndr.c                      0.23 0.01 26.80 16.43 0.00000  0.20 0.26
## year.c                      0.00 0.00 24.82  1.74 0.09360  0.00 0.01
## gndr.c:year.c               0.00 0.00 24.39 -4.16 0.00034 -0.01 0.00
## gndr.c:log_gdp.z.cm        -0.03 0.01 25.86 -1.98 0.05870 -0.06 0.00
## year.c:log_gdp.z.cm         0.00 0.00 25.55 -0.82 0.42131 -0.01 0.00
## gndr.c:year.c:log_gdp.z.cm  0.00 0.00 26.25  1.36 0.18635  0.00 0.00
```

``` r
getVC(mod6_log_GDP)
```

```
##         grp          var1          var2 sdcor vcov
## 1     cntry   (Intercept)          <NA>  0.19 0.04
## 2     cntry        gndr.c          <NA>  0.07 0.00
## 3     cntry        year.c          <NA>  0.02 0.00
## 4     cntry gndr.c:year.c          <NA>  0.00 0.00
## 5     cntry   (Intercept)        gndr.c -0.16 0.00
## 6     cntry   (Intercept)        year.c -0.42 0.00
## 7     cntry   (Intercept) gndr.c:year.c -0.07 0.00
## 8     cntry        gndr.c        year.c -0.41 0.00
## 9     cntry        gndr.c gndr.c:year.c -0.78 0.00
## 10    cntry        year.c gndr.c:year.c  0.42 0.00
## 11 Residual          <NA>          <NA>  0.99 0.99
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_log_GDP,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] 1.471294
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_log_GDP,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] 0.8257686
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
##  year.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21 0.1618 0.0720 26.6   0.0140   0.3096   2.248  0.0331
##       0 0.0221 0.0330 29.5  -0.0454   0.0896   0.668  0.5094
## 
## log_gdp.z.cm =  0:
##  year.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21 0.1200 0.0515 26.4   0.0142   0.2259   2.330  0.0277
##       0 0.0221 0.0330 29.5  -0.0454   0.0896   0.668  0.5094
## 
## log_gdp.z.cm =  1:
##  year.c emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21 0.0783 0.0731 25.9  -0.0720   0.2286   1.071  0.2942
##       0 0.0221 0.0330 29.5  -0.0454   0.0896   0.668  0.5094
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_mod6_log_GDP,adjust="none",infer=c(T,T))
```

```
## log_gdp.z.cm = -1:
##  contrast           estimate     SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0   0.1397 0.0757 28.2  -0.0152    0.295   1.847  0.0753
## 
## log_gdp.z.cm =  0:
##  contrast           estimate     SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0   0.0980 0.0562 24.8  -0.0178    0.214   1.744  0.0936
## 
## log_gdp.z.cm =  1:
##  contrast           estimate     SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0   0.0562 0.0762 26.3  -0.1004    0.213   0.738  0.4673
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
##    -0.5     21  0.08395 0.0727 25.6 -0.06565  0.23355   1.154  0.2590
##     0.5     21  0.23964 0.0721 25.9  0.09139  0.38790   3.323  0.0027
##    -0.5      0 -0.10616 0.0355 33.5 -0.17832 -0.03401  -2.992  0.0052
##     0.5      0  0.15028 0.0338 30.6  0.08129  0.21928   4.445  0.0001
## 
## log_gdp.z.cm =  0:
##  gndr.c year.c   emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21  0.04374 0.0522 25.5 -0.06357  0.15106   0.839  0.4094
##     0.5     21  0.19634 0.0515 25.6  0.09046  0.30223   3.814  0.0008
##    -0.5      0 -0.09188 0.0346 30.4 -0.16257 -0.02120  -2.653  0.0126
##     0.5      0  0.13600 0.0328 27.6  0.06869  0.20331   4.141  0.0003
## 
## log_gdp.z.cm =  1:
##  gndr.c year.c   emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21  0.00354 0.0738 24.9 -0.14853  0.15562   0.048  0.9621
##     0.5     21  0.15304 0.0732 25.1  0.00237  0.30371   2.091  0.0468
##    -0.5      0 -0.07760 0.0353 32.5 -0.14939 -0.00581  -2.200  0.0350
##     0.5      0  0.12172 0.0334 29.3  0.05336  0.19008   3.640  0.0010
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
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21    -0.1557 0.0158 35.0  -0.1877  -0.1237  -9.869 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)   0.1901 0.0727 25.8   0.0406   0.3396   2.615  0.0147
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     -0.0663 0.0801 27.8  -0.2304   0.0977  -0.828  0.4145
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)      0.3458 0.0731 28.1   0.1960   0.4956   4.728 <0.0001
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0         0.0894 0.0809 28.2  -0.0762   0.2550   1.105  0.2785
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0      -0.2564 0.0210 29.5  -0.2994  -0.2135 -12.206 <0.0001
## 
## log_gdp.z.cm =  0:
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21    -0.1526 0.0110 31.5  -0.1749  -0.1303 -13.935 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)   0.1356 0.0539 22.9   0.0241   0.2472   2.516  0.0194
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     -0.0923 0.0587 24.8  -0.2132   0.0287  -1.572  0.1286
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)      0.2882 0.0550 24.5   0.1749   0.4015   5.244 <0.0001
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0         0.0603 0.0598 24.9  -0.0628   0.1835   1.010  0.3224
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0      -0.2279 0.0139 26.8  -0.2564  -0.1994 -16.426 <0.0001
## 
## log_gdp.z.cm =  1:
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21    -0.1495 0.0149 29.8  -0.1800  -0.1190 -10.010 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)   0.0811 0.0728 23.7  -0.0693   0.2316   1.114  0.2765
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     -0.1182 0.0807 25.9  -0.2840   0.0476  -1.465  0.1549
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)      0.2306 0.0735 26.0   0.0796   0.3817   3.139  0.0042
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0         0.0313 0.0813 26.2  -0.1358   0.1985   0.385  0.7033
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0      -0.1993 0.0190 23.4  -0.2386  -0.1601 -10.495 <0.0001
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
##  diff_ESS11    0.156 0.0158 35.0    0.124    0.188   9.869 <0.0001
##  diff_ESS1     0.256 0.0210 29.5    0.214    0.299  12.206 <0.0001
## 
## log_gdp.z.cm =  0:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11    0.153 0.0110 31.5    0.130    0.175  13.935 <0.0001
##  diff_ESS1     0.228 0.0139 26.8    0.199    0.256  16.426 <0.0001
## 
## log_gdp.z.cm =  1:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11    0.149 0.0149 29.8    0.119    0.180  10.010 <0.0001
##  diff_ESS1     0.199 0.0190 23.4    0.160    0.239  10.495 <0.0001
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
##  diff_ESS11 - diff_ESS1  -0.1008 0.0275 30.0   -0.157 -0.04461  -3.665  0.0009
## 
## log_gdp.z.cm =  0:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1  -0.0753 0.0181 24.4   -0.113 -0.03799  -4.163  0.0003
## 
## log_gdp.z.cm =  1:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1  -0.0498 0.0246 20.8   -0.101  0.00126  -2.029  0.0554
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
##   [1] mnormt_2.1.2       Rdpack_2.6.6       gridExtra_2.3      writexl_1.5.4      readxl_1.4.5      
##   [6] rlang_1.2.0        magrittr_2.0.5     otel_0.2.0         rockchalk_1.8.164  compiler_4.6.0    
##  [11] mgcv_1.9-4         png_0.1-9          vctrs_0.7.3        quadprog_1.5-8     pkgconfig_2.0.3   
##  [16] shape_1.4.6.1      fastmap_1.2.0      backports_1.5.1    labeling_0.4.3     pbivnorm_0.6.0    
##  [21] utf8_1.2.6         nloptr_2.2.1       purrr_1.2.2        xfun_0.57          glmnet_5.0        
##  [26] jomo_2.7-6         cachem_1.1.0       kutils_1.73        jsonlite_2.0.0     pan_1.9           
##  [31] jpeg_0.1-11        psych_2.6.3        lavaan_0.6-21      parallel_4.6.0     broom_1.0.13      
##  [36] cluster_2.1.8.2    R6_2.6.1           bslib_0.10.0       stringi_1.8.7      RColorBrewer_1.1-3
##  [41] car_3.1-5          boot_1.3-32        rpart_4.1.27       cellranger_1.1.0   jquerylib_0.1.4   
##  [46] estimability_1.5.1 Rcpp_1.1.1-1.1     iterators_1.0.14   base64enc_0.1-6    R.utils_2.13.0    
##  [51] splines_4.6.0      nnet_7.3-20        tidyselect_1.2.1   rstudioapi_0.18.0  abind_1.4-8       
##  [56] yaml_2.3.12        codetools_0.2-20   lattice_0.22-9     plyr_1.8.9         withr_3.0.2       
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

