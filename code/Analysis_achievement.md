---
title: "Analysis for achievement values"
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
cntry.ach<-diff_dat %>% group_by(cntry,essround) %>%
  summarise(ach.ctm=mean(ach,na.rm=T),
            ach.ctsd=sd(ach,na.rm=T)) %>%
  group_by(cntry) %>%
  summarise(ach.cm=mean(ach.ctm),
            ach.csd=mean(ach.ctsd)) 
```

```
## `summarise()` has regrouped the output.
## ℹ Summaries were computed grouped by cntry and essround.
## ℹ Output is grouped by cntry.
## ℹ Use `summarise(.groups = "drop_last")` to silence this message.
## ℹ Use `summarise(.by = c(cntry, essround))` for per-operation grouping (`?dplyr::dplyr_by`) instead.
```

``` r
grand_mean_ach<-mean(cntry.ach$ach.cm)
grand_sd_ach<-mean(cntry.ach$ach.csd)

# standardized
diff_dat$ach.z<-(diff_dat$ach-grand_mean_ach)/grand_sd_ach
hist(diff_dat$ach.z)
```

![](Analysis_achievement_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

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
                   ach.z.wt=weighted.mean(x=ach.z,w=pspwght),
                   ach.z=mean(ach.z),
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

# achievement

cntry_ach_frame<-
  diff_dat %>% group_by(cntry,essround) %>%
  summarise('ach M' = weighted.mean(x=ach.z,w=pspwght),
            'ach SD' = sqrt(wtd.var(ach.z,pspwght))) %>%
  group_by(cntry) %>%
  summarise('ach M' = mean(x=`ach M`),
            'ach SD'= mean(x=`ach SD`))
```

```
## `summarise()` has regrouped the output.
## ℹ Summaries were computed grouped by cntry and essround.
## ℹ Output is grouped by cntry.
## ℹ Use `summarise(.groups = "drop_last")` to silence this message.
## ℹ Use `summarise(.by = c(cntry, essround))` for per-operation grouping (`?dplyr::dplyr_by`) instead.
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
## `summarise()` has regrouped the output.
## ℹ Summaries were computed grouped by cntry and essround.
## ℹ Output is grouped by cntry.
## ℹ Use `summarise(.groups = "drop_last")` to silence this message.
## ℹ Use `summarise(.by = c(cntry, essround))` for per-operation grouping (`?dplyr::dplyr_by`) instead.
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
## `summarise()` has regrouped the output.
## ℹ Summaries were computed grouped by cntry and essround.
## ℹ Output is grouped by cntry.
## ℹ Use `summarise(.groups = "drop_last")` to silence this message.
## ℹ Use `summarise(.by = c(cntry, essround))` for per-operation grouping (`?dplyr::dplyr_by`) instead.
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
## # A tibble: 34 × 10
##    cntry `n ESS rounds`     n `ach M` `ach SD` `ach M Women` `ach SD Women` `ach M Men` `ach SD Men`
##    <chr>          <dbl> <dbl>   <dbl>    <dbl>         <dbl>          <dbl>       <dbl>        <dbl>
##  1 AT                 7 15400  0.230     0.996        0.0989          1.03       0.372         0.941
##  2 BE                11 18886 -0.0140    0.977       -0.0785          1.00       0.0539        0.946
##  3 BG                 7 14857  0.598     0.919        0.521           0.952      0.680         0.873
##  4 CH                11 18087 -0.0278    1.00        -0.162           1.02       0.113         0.966
##  5 CY                 6  5771  0.183     1.05         0.140           1.08       0.229         1.02 
##  6 CZ                 9 18934 -0.175     1.04        -0.320           1.06      -0.0185        0.986
##  7 DE                10 27753 -0.157     0.997       -0.291           1.01      -0.0158        0.963
##  8 DK                 8 12198 -0.215     1.03        -0.304           1.02      -0.124         1.03 
##  9 EE                10 17974 -0.299     1.01        -0.380           1.02      -0.201         0.987
## 10 ES                10 18785 -0.271     1.08        -0.375           1.08      -0.162         1.06 
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
print(cntry_desc_tbl,n=35)
```

```
## # A tibble: 34 × 14
##    Country     `n ESS rounds`     n `ach M` `ach SD` `ach M Women` `ach SD Women` `ach M Men` `ach SD Men`
##    <chr>                <dbl> <dbl> <chr>   <chr>    <chr>         <chr>          <chr>       <chr>       
##  1 Austria                  7 15400 0.23    1.00     0.10          1.03           0.37        0.94        
##  2 Belgium                 11 18886 -0.01   0.98     -0.08         1.00           0.05        0.95        
##  3 Bulgaria                 7 14857 0.60    0.92     0.52          0.95           0.68        0.87        
##  4 Switzerland             11 18087 -0.03   1.00     -0.16         1.02           0.11        0.97        
##  5 Cyprus                   6  5771 0.18    1.05     0.14          1.08           0.23        1.02        
##  6 Czechia                  9 18934 -0.18   1.04     -0.32         1.06           -0.02       0.99        
##  7 Germany                 10 27753 -0.16   1.00     -0.29         1.01           -0.02       0.96        
##  8 Denmark                  8 12198 -0.22   1.03     -0.30         1.02           -0.12       1.03        
##  9 Estonia                 10 17974 -0.30   1.01     -0.38         1.02           -0.20       0.99        
## 10 Spain                   10 18785 -0.27   1.08     -0.38         1.08           -0.16       1.06        
## 11 Finland                 11 19568 -0.49   1.03     -0.64         1.03           -0.32       0.99        
## 12 France                  11 20457 -0.58   1.09     -0.69         1.07           -0.46       1.09        
## 13 UK                      11 22979 -0.07   1.07     -0.18         1.09           0.05        1.03        
## 14 Greece                   6 15212 0.31    1.01     0.23          1.03           0.39        0.97        
## 15 Croatia                  5  7914 0.03    1.06     -0.05         1.08           0.11        1.04        
## 16 Hungary                 11 18123 0.31    0.98     0.24          1.00           0.39        0.96        
## 17 Ireland                 11 22562 0.08    1.05     0.01          1.08           0.16        1.01        
## 18 Israel                   7 14857 0.71    0.91     0.69          0.91           0.72        0.89        
## 19 Iceland                  6  4654 -0.53   1.05     -0.56         1.06           -0.49       1.05        
## 20 Italy                    5 11441 0.45    0.89     0.37          0.93           0.53        0.85        
## 21 Lithuania                7 13059 0.07    1.00     0.03          1.02           0.12        0.97        
## 22 Latvia                   3  4088 0.07    1.03     -0.02         1.03           0.17        1.01        
## 23 Montenegro               3  4028 0.16    0.98     0.13          1.00           0.19        0.95        
## 24 Netherlands             11 19722 -0.10   0.94     -0.22         0.97           0.02        0.89        
## 25 Norway                  11 16505 -0.35   0.99     -0.42         1.01           -0.28       0.95        
## 26 Poland                  10 16737 0.07    0.97     -0.03         1.00           0.18        0.92        
## 27 Portugal                11 19070 0.11    0.90     0.03          0.92           0.19        0.88        
## 28 Serbia                   2  3499 0.31    1.07     0.23          1.10           0.38        1.04        
## 29 Russia                   5 12139 0.23    1.04     0.17          1.07           0.31        1.01        
## 30 Sweden                  10 16104 -0.49   1.02     -0.55         1.02           -0.43       1.01        
## 31 Slovenia                11 14463 0.46    0.88     0.40          0.92           0.53        0.84        
## 32 Slovakia                 8 12547 0.14    0.95     0.03          0.98           0.25        0.90        
## 33 Turkey                   2  4108 0.69    0.84     0.59          0.87           0.78        0.78        
## 34 Ukraine                  6 12054 -0.13   1.07     -0.20         1.09           -0.04       1.05        
## # ℹ 5 more variables: D <chr>, GEI <chr>, GGGI <chr>, GDI <chr>, GDP <chr>
```

``` r
export(cntry_desc_tbl,"../results/ach/cntry_desc_tbl.xlsx",overwrite=T)
```

## Country-level correlations table


``` r
cor_frame<-
  desc_frame %>%
  dplyr::select(
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
  dplyr::select(-GDP)

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
##   2. VBMT_Women -0.05 0.35 1.00                                                                      
##                            [.99, 1.00]                                                               
##                                                                                                      
##   3. VBMT_Men   0.13  0.33 .99          .98                                                          
##                            [.99, 1.00]  [.96, .99]                                                   
##                                                                                                      
##   4. D          0.17  0.07 -.30         -.39         -.21                                            
##                            [-.58, .04]  [-.64, -.06] [-.51, .14]                                     
##                                                                                                      
##   5. GEI        0.87  0.07 -.60         -.60         -.60         .16                                
##                            [-.78, -.33] [-.78, -.32] [-.78, -.33] [-.19, .48]                        
##                                                                                                      
##   6. GGGI       0.74  0.05 -.70         -.68         -.71         .05         .73                    
##                            [-.84, -.47] [-.83, -.44] [-.85, -.49] [-.29, .38] [.52, .86]             
##                                                                                                      
##   7. GDI        0.98  0.03 -.24         -.22         -.25         -.11        .07         .19        
##                            [-.53, .11]  [-.52, .13]  [-.54, .09]  [-.43, .24] [-.28, .41] [-.16, .50]
##                                                                                                      
##   8. log_GDP    10.61 0.41 -.46         -.47         -.45         .24         .72         .62        
##                            [-.69, -.14] [-.70, -.16] [-.68, -.13] [-.10, .54] [.50, .85]  [.36, .79] 
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
## 1467134.5 1467167.8 -733564.2 1467128.5    492340 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -5.7591 -0.6821  0.0179  0.6540  5.2297 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  cntry    (Intercept) 0.1104   0.3322  
##  Residual             1.0129   1.0064  
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##             Estimate Std. Error       df t value Pr(>|t|)
## (Intercept)  0.04332    0.05700 33.99528    0.76    0.452
```

``` r
getVC(mod0)
```

```
##        grp        var1 var2 sdcor vcov
## 1    cntry (Intercept) <NA>  0.33 0.11
## 2 Residual        <NA> <NA>  1.01 1.01
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
## mean variation  0.0982484     NA       1
## sigma2          0.9017516      1      NA
## 
## $R2s
##         total within between
## f1  0.0000000      0      NA
## f2  0.0000000     NA       0
## v   0.0000000      0      NA
## m   0.0982484     NA       1
## f   0.0000000     NA      NA
## fv  0.0000000      0      NA
## fvm 0.0982484     NA      NA
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
## 1462740.7 1462785.1 -731366.3 1462732.7    492339 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.0167 -0.6692  0.0350  0.6596  5.0393 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  cntry    (Intercept) 0.1104   0.3323  
##  Residual             1.0039   1.0020  
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##              Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept) 4.692e-02  5.701e-02 3.400e+01   0.823    0.416    
## gndr.c      1.896e-01  2.853e-03 4.923e+05  66.449   <2e-16 ***
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
## (Intercept) 0.047 0.057     33.995  0.823 0.416 -0.069 0.163
## gndr.c      0.190 0.003 492309.619 66.449 0.000  0.184 0.195
```

``` r
getVC(mod1)
```

```
##        grp        var1 var2 sdcor vcov
## 1    cntry (Intercept) <NA>  0.33 0.11
## 2 Residual        <NA> <NA>  1.00 1.00
```

``` r
r2mlm(mod1,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.007949707
## slope variation 0.000000000
## mean variation  0.098291747
## sigma2          0.893758546
## 
## $R2s
##           total
## f   0.007949707
## v   0.000000000
## m   0.098291747
## fv  0.007949707
## fvm 0.106241454
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
##   1462318   1462385   -731153   1462306    492337 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -5.9123 -0.6698  0.0267  0.6595  5.0734 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr  
##  cntry    (Intercept) 0.110115 0.33184        
##           gndr.c      0.004121 0.06419  -0.32 
##  Residual             1.002867 1.00143        
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##             Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)  0.04673    0.05693 33.99465   0.821    0.417    
## gndr.c       0.17676    0.01148 33.96039  15.396   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##        (Intr)
## gndr.c -0.307
```

``` r
getFE(mod2,round=3)
```

```
##              Est.    SE     df      t     p     LL    UL
## (Intercept) 0.047 0.057 33.995  0.821 0.417 -0.069 0.162
## gndr.c      0.177 0.011 33.960 15.396 0.000  0.153 0.200
```

``` r
getVC(mod2)
```

```
##        grp        var1   var2 sdcor  vcov
## 1    cntry (Intercept)   <NA>  0.33  0.11
## 2    cntry      gndr.c   <NA>  0.06  0.00
## 3    cntry (Intercept) gndr.c -0.32 -0.01
## 4 Residual        <NA>   <NA>  1.00  1.00
```

``` r
r2mlm(mod2,bargraph = F)
```

```
## $Decompositions
##                        total
## fixed           0.0069180539
## slope variation 0.0009123378
## mean variation  0.0985929934
## sigma2          0.8935766148
## 
## $R2s
##            total
## f   0.0069180539
## v   0.0009123378
## m   0.0985929934
## fv  0.0078303917
## fvm 0.1064233852
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
## mod1    4 1462741 1462785 -731366   1462733                         
## mod2    6 1462318 1462385 -731153   1462306 426.63  2  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

``` r
lvl2_var_cond_lvl1(mod2,lvl1.var = "gndr.c",lvl1.values = c(0.5,-0.5))
```

```
##   lvl1.value lvl2.cond.var lvl2.cond.sd
## 1        0.5     0.1043160    0.3229799
## 2       -0.5     0.1179742    0.3434738
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
## 1462319.4 1462374.9 -731154.7 1462309.4    492338 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -5.9115 -0.6698  0.0265  0.6594  5.0718 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  cntry    (Intercept) 0.11013  0.33186 
##  cntry.1  gndr.c      0.00417  0.06457 
##  Residual             1.00287  1.00143 
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##             Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)  0.04672    0.05694 33.99697   0.821    0.418    
## gndr.c       0.17692    0.01155 33.92148  15.323   <2e-16 ***
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
## (Intercept) 0.047 0.057 33.997  0.821 0.418 -0.069 0.162
## gndr.c      0.177 0.012 33.921 15.323 0.000  0.153 0.200
```

``` r
getVC(mod2_norecov)
```

```
##        grp        var1 var2 sdcor vcov
## 1    cntry (Intercept) <NA>  0.33 0.11
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
## mod2_norecov    5 1462319 1462375 -731155   1462309                       
## mod2            6 1462318 1462385 -731153   1462306 3.3035  1    0.06913 .
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
## 1422816.8 1422905.4 -711400.4 1422800.8    480356 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -5.9243 -0.6706  0.0297  0.6609  5.0836 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr  
##  cntry    (Intercept) 0.069961 0.26450        
##           gndr.c      0.004098 0.06401  -0.28 
##  Residual             0.998824 0.99941        
## Number of obs: 480364, groups:  cntry, 33
## 
## Fixed effects:
##                 Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)      0.05185    0.04607 33.01353   1.125    0.269    
## gndr.c           0.17745    0.01163 32.78606  15.261  < 2e-16 ***
## gei.z.cm        -0.20978    0.04680 33.05527  -4.482 8.39e-05 ***
## gndr.c:gei.z.cm  0.01204    0.01203 35.13897   1.001    0.324    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c g.z.cm
## gndr.c      -0.264              
## gei.z.cm     0.000  0.000       
## gndr.c:g.z.  0.000 -0.026 -0.259
```

``` r
getFE(mod2_GEI,round=3)
```

```
##                   Est.    SE     df      t     p     LL     UL
## (Intercept)      0.052 0.046 33.014  1.125 0.269 -0.042  0.146
## gndr.c           0.177 0.012 32.786 15.261 0.000  0.154  0.201
## gei.z.cm        -0.210 0.047 33.055 -4.482 0.000 -0.305 -0.115
## gndr.c:gei.z.cm  0.012 0.012 35.139  1.001 0.324 -0.012  0.036
```

``` r
getVC(mod2_GEI)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.26 0.07
## 2    cntry      gndr.c   <NA>  0.06 0.00
## 3    cntry (Intercept) gndr.c -0.28 0.00
## 4 Residual        <NA>   <NA>  1.00 1.00
```

``` r
r2mlm(mod2_GEI,bargraph = F)
```

```
## $Decompositions
##                        total
## fixed           0.0341914549
## slope variation 0.0009195811
## mean variation  0.0634522936
## sigma2          0.9014366704
## 
## $R2s
##            total
## f   0.0341914549
## v   0.0009195811
## m   0.0634522936
## fv  0.0351110360
## fvm 0.0985633296
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
## Time difference of 29.83783 secs
```

``` r
round(ddsc_mod2_GEI$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.118        0.343        1.003     1.121 0.105   7802.647 0.999   0.999
## 2        0.5         0.104        0.323        1.003     1.107 0.094   6678.029 0.998   0.999
```

``` r
round(ddsc_mod2_GEI$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gei.z.cm gei.z.cm_scaled
## means_y1           0.103 0.339    1.000           1.000    0.983           0.983   -0.595          -0.595
## means_y1_scaled    0.297 0.976    1.000           1.000    0.983           0.983   -0.595          -0.595
## means_y2          -0.067 0.356    0.983           0.983    1.000           1.000   -0.586          -0.586
## means_y2_scaled   -0.194 1.024    0.983           0.983    1.000           1.000   -0.586          -0.586
## gei.z.cm           0.000 1.000   -0.595          -0.595   -0.586          -0.586    1.000           1.000
## gei.z.cm_scaled    0.000 1.000   -0.595          -0.595   -0.586          -0.586    1.000           1.000
## diff_score         0.171 0.066   -0.164          -0.164   -0.341          -0.341    0.103           0.103
## diff_score_scaled  0.491 0.189   -0.164          -0.164   -0.341          -0.341    0.103           0.103
##                   diff_score diff_score_scaled
## means_y1              -0.164            -0.164
## means_y1_scaled       -0.164            -0.164
## means_y2              -0.341            -0.341
## means_y2_scaled       -0.341            -0.341
## gei.z.cm               0.103             0.103
## gei.z.cm_scaled        0.103             0.103
## diff_score             1.000             1.000
## diff_score_scaled      1.000             1.000
```

``` r
round(ddsc_mod2_GEI$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                      -0.183 0.183 35.139  -1.001   0.324   -0.554    0.188
## w_11                         -0.216 0.049 33.103  -4.430   0.000   -0.315   -0.117
## w_21                         -0.204 0.046 33.106  -4.467   0.000   -0.297   -0.111
## r_xy1                        -0.636 0.144 33.103  -4.430   0.000   -0.928   -0.344
## r_xy2                        -0.572 0.128 33.106  -4.467   0.000   -0.833   -0.312
## b_11                         -0.621 0.140 33.103  -4.430   0.000   -0.906   -0.336
## b_21                         -0.586 0.131 33.106  -4.467   0.000   -0.853   -0.319
## main_effect                  -0.210 0.047 33.055  -4.482   0.000   -0.305   -0.115
## moderator_effect              0.177 0.012 32.786  15.261   0.000    0.154    0.201
## interaction                   0.012 0.012 35.139   1.001   0.324   -0.012    0.036
## q_b11_b21                    -0.054    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                  -0.101    NA     NA      NA      NA       NA       NA
## cross_over_point            -14.742    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.198 0.045 33.237  -4.375   0.000   -0.290   -0.106
## interaction_vs_main_bscale   -0.569 0.130 33.237  -4.375   0.000   -0.833   -0.304
## interaction_vs_main_rscale   -0.540 0.124 33.249  -4.364   0.000   -0.792   -0.289
## dadas                        -0.408 0.091 33.106  -4.467   1.000   -0.593   -0.222
## dadas_bscale                 -1.172 0.262 33.106  -4.467   1.000   -1.706   -0.638
## dadas_rscale                 -1.144 0.256 33.106  -4.467   1.000   -1.666   -0.623
## abs_diff                      0.012 0.012 35.139   1.001   0.162   -0.012    0.036
## abs_sum                       0.420 0.094 33.055   4.482   0.000    0.229    0.610
## abs_diff_bscale               0.035 0.035 35.139   1.001   0.162   -0.036    0.105
## abs_sum_bscale                1.207 0.269 33.055   4.482   0.000    0.659    1.754
## abs_diff_rscale               0.064 0.037 35.290   1.731   0.046   -0.011    0.139
## abs_sum_rscale                1.208 0.270 33.055   4.482   0.000    0.660    1.757
```

``` r
round(ddsc_mod2_GEI$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.007 -0.321  3.304  1.000  0.069
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
##                                     est     se       z pvalue ci.lower ci.upper
## r_xy1_y2                         -0.103  0.173  -0.597  0.551   -0.443    0.236
## r_xy1                            -0.586  0.141  -4.157  0.000   -0.863   -0.310
## r_xy2                            -0.595  0.140  -4.255  0.000   -0.869   -0.321
## b_11                             -0.600  0.144  -4.157  0.000   -0.883   -0.317
## b_21                             -0.581  0.136  -4.255  0.000   -0.848   -0.313
## b_10                             -0.194  0.142  -1.364  0.173   -0.473    0.085
## b_20                              0.297  0.134   2.211  0.027    0.034    0.561
## res_cov_y1_y2                     0.614  0.153   4.009  0.000    0.314    0.915
## diff_b10_b20                     -0.491  0.032 -15.228  0.000   -0.554   -0.428
## diff_b11_b21                     -0.020  0.033  -0.597  0.551   -0.084    0.045
## diff_rxy1_rxy2                    0.009  0.032   0.281  0.778   -0.053    0.071
## q_b11_b21                        -0.030  0.053  -0.570  0.569   -0.133    0.073
## q_rxy1_rxy2                       0.014  0.049   0.281  0.778   -0.082    0.110
## cross_over_point                -25.117 42.105  -0.597  0.551 -107.641   57.407
## sum_b11_b21                      -1.181  0.279  -4.232  0.000   -1.728   -0.634
## main_effect                      -0.590  0.140  -4.232  0.000   -0.864   -0.317
## interaction_vs_main_effect       -0.571  0.135  -4.218  0.000   -0.836   -0.306
## diff_abs_b11_abs_b21              0.020  0.033   0.597  0.551   -0.045    0.084
## abs_diff_b11_b21                  0.020  0.033   0.597  0.275   -0.045    0.084
## abs_sum_b11_b21                   1.181  0.279   4.232  0.000    0.634    1.728
## dadas                            -1.161  0.273  -4.255  1.000   -1.696   -0.626
## q_r_equivalence                  -0.086  0.049  -1.765  0.039       NA       NA
## q_b_equivalence                  -0.070  0.053  -1.329  0.092       NA       NA
## cross_over_point_equivalence     25.117 42.105   0.597  0.725       NA       NA
## cross_over_point_minimal_effect  25.117 42.105   0.597  0.275       NA       NA
```

``` r
round(ddsc_sem_GEI$variance_test,3)
```

```
##             est    se       z pvalue ci.lower ci.upper
## cov_y1y2  0.952 0.236   4.028  0.000    0.489    1.416
## var_y1    1.016 0.250   4.062  0.000    0.526    1.507
## var_y2    0.923 0.227   4.062  0.000    0.478    1.368
## var_diff  0.094 0.066   1.425  0.154   -0.035    0.222
## var_ratio 1.101 0.070  15.763  0.000    0.964    1.238
## cor_y1y2  0.983 0.006 170.111  0.000    0.972    0.995
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
## r_xy1y2                      -0.103 0.179 31.000  -0.579   0.567   -0.468    0.261
## w_11                         -0.209 0.050 31.854  -4.141   0.000   -0.311   -0.106
## w_21                         -0.202 0.050 31.854  -4.006   0.000   -0.305   -0.099
## r_xy1                        -0.586 0.142 31.854  -4.141   0.000   -0.875   -0.298
## r_xy2                        -0.595 0.149 31.854  -4.006   0.000   -0.898   -0.293
## b_11                         -0.600 0.145 31.854  -4.141   0.000   -0.896   -0.305
## b_21                         -0.581 0.145 31.854  -4.006   0.000   -0.876   -0.285
## main_effect                  -0.205 0.050 31.000  -4.102   0.000   -0.307   -0.103
## moderator_effect              0.171 0.012 31.000  14.759   0.000    0.147    0.194
## interaction                   0.007 0.012 31.000   0.579   0.567   -0.017    0.031
## q_b11_b21                    -0.030    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.014    NA     NA      NA      NA       NA       NA
## cross_over_point            -25.117    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.199 0.051 34.405  -3.861   0.000   -0.303   -0.094
## interaction_vs_main_bscale   -0.571 0.148 34.405  -3.861   0.000   -0.872   -0.271
## interaction_vs_main_rscale   -0.600 0.155 34.175  -3.872   0.000   -0.914   -0.285
## dadas                        -0.404 0.101 31.854  -4.006   1.000   -0.609   -0.199
## dadas_bscale                 -1.162 0.290 31.854  -4.006   1.000   -1.752   -0.571
## dadas_rscale                 -1.190 0.297 31.854  -4.006   1.000   -1.796   -0.585
## abs_diff                      0.007 0.012 31.000   0.579   0.284   -0.017    0.031
## abs_sum                       0.411 0.100 31.000   4.102   0.000    0.206    0.615
## abs_diff_bscale               0.020 0.034 31.000   0.579   0.284   -0.049    0.088
## abs_sum_bscale                1.181 0.288 31.000   4.102   0.000    0.594    1.769
## abs_diff_rscale              -0.009 0.035 33.616  -0.259   0.601   -0.079    0.061
## abs_sum_rscale                1.181 0.288 31.000   4.100   0.000    0.594    1.769
```

``` r
# country-time multilevel model


mod2_GEI_cntry_year<-
  lmer(ach.z.wt~gndr.c+gei.z.cm+gndr.c:gei.z.cm+
                 (gndr.c|cntry),data=diff_dat_cntry_year,REML=F,
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2_GEI_cntry_year)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ach.z.wt ~ gndr.c + gei.z.cm + gndr.c:gei.z.cm + (gndr.c | cntry)
##    Data: diff_dat_cntry_year
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##    -615.6    -581.3     315.8    -631.6       526 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -3.6871 -0.6073  0.0101  0.6505  3.5345 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr  
##  cntry    (Intercept) 0.07002  0.26461        
##           gndr.c      0.00114  0.03377  -0.57 
##  Residual             0.01354  0.11638        
## Number of obs: 534, groups:  cntry, 33
## 
## Fixed effects:
##                 Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)      0.04663    0.04641 33.13399   1.005 0.322291    
## gndr.c           0.17849    0.01201 34.16921  14.866  < 2e-16 ***
## gei.z.cm        -0.20485    0.04735 33.74527  -4.326 0.000127 ***
## gndr.c:gei.z.cm  0.01562    0.01356 40.15027   1.152 0.255988    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c g.z.cm
## gndr.c      -0.277              
## gei.z.cm    -0.006  0.001       
## gndr.c:g.z.  0.001 -0.188 -0.250
```

``` r
getFE(mod2_GEI_cntry_year,round=3)
```

```
##                   Est.    SE     df      t     p     LL     UL
## (Intercept)      0.047 0.046 33.134  1.005 0.322 -0.048  0.141
## gndr.c           0.178 0.012 34.169 14.866 0.000  0.154  0.203
## gei.z.cm        -0.205 0.047 33.745 -4.326 0.000 -0.301 -0.109
## gndr.c:gei.z.cm  0.016 0.014 40.150  1.152 0.256 -0.012  0.043
```

``` r
getVC(mod2_GEI_cntry_year)
```

```
##        grp        var1   var2 sdcor  vcov
## 1    cntry (Intercept)   <NA>  0.26  0.07
## 2    cntry      gndr.c   <NA>  0.03  0.00
## 3    cntry (Intercept) gndr.c -0.57 -0.01
## 4 Residual        <NA>   <NA>  0.12  0.01
```

``` r
r2mlm(mod2_GEI,bargraph = F)
```

```
## $Decompositions
##                        total
## fixed           0.0341914549
## slope variation 0.0009195811
## mean variation  0.0634522936
## sigma2          0.9014366704
## 
## $R2s
##            total
## f   0.0341914549
## v   0.0009195811
## m   0.0634522936
## fv  0.0351110360
## fvm 0.0985633296
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
## 1       -0.5         0.116        0.340        0.014     0.129 0.895      8.029 0.998   0.986
## 2        0.5         0.100        0.316        0.014     0.114 0.880      8.029 0.998   0.983
```

``` r
round(ddsc_mod2_GEI_cntry_year$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gei.z.cm gei.z.cm_scaled
## means_y1           0.134 0.330    1.000           1.000    0.982           0.982   -0.603          -0.603
## means_y1_scaled    0.394 0.969    1.000           1.000    0.982           0.982   -0.603          -0.603
## means_y2          -0.041 0.351    0.982           0.982    1.000           1.000   -0.599          -0.599
## means_y2_scaled   -0.119 1.030    0.982           0.982    1.000           1.000   -0.599          -0.599
## gei.z.cm           0.000 1.000   -0.603          -0.603   -0.599          -0.599    1.000           1.000
## gei.z.cm_scaled    0.000 1.000   -0.603          -0.603   -0.599          -0.599    1.000           1.000
## diff_score         0.175 0.068   -0.211          -0.211   -0.393          -0.393    0.165           0.165
## diff_score_scaled  0.514 0.201   -0.211          -0.211   -0.393          -0.393    0.165           0.165
##                   diff_score diff_score_scaled
## means_y1              -0.211            -0.211
## means_y1_scaled       -0.211            -0.211
## means_y2              -0.393            -0.393
## means_y2_scaled       -0.393            -0.393
## gei.z.cm               0.165             0.165
## gei.z.cm_scaled        0.165             0.165
## diff_score             1.000             1.000
## diff_score_scaled      1.000             1.000
```

``` r
round(ddsc_mod2_GEI_cntry_year$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                      -0.228 0.198 40.150  -1.152   0.256   -0.628    0.172
## w_11                         -0.213 0.049 34.016  -4.298   0.000   -0.313   -0.112
## w_21                         -0.197 0.046 34.077  -4.272   0.000   -0.291   -0.103
## r_xy1                        -0.644 0.150 34.016  -4.298   0.000   -0.949   -0.340
## r_xy2                        -0.561 0.131 34.077  -4.272   0.000   -0.828   -0.294
## b_11                         -0.624 0.145 34.016  -4.298   0.000   -0.920   -0.329
## b_21                         -0.579 0.135 34.077  -4.272   0.000   -0.854   -0.303
## main_effect                  -0.205 0.047 33.745  -4.326   0.000   -0.301   -0.109
## moderator_effect              0.178 0.012 34.169  14.866   0.000    0.154    0.203
## interaction                   0.016 0.014 40.150   1.152   0.256   -0.012    0.043
## q_b11_b21                    -0.072    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                  -0.130    NA     NA      NA      NA       NA       NA
## cross_over_point            -11.424    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.189 0.046 34.994  -4.124   0.000   -0.282   -0.096
## interaction_vs_main_bscale   -0.556 0.135 34.994  -4.124   0.000   -0.829   -0.282
## interaction_vs_main_rscale   -0.520 0.127 35.109  -4.105   0.000   -0.777   -0.263
## dadas                        -0.394 0.092 34.077  -4.272   1.000   -0.582   -0.207
## dadas_bscale                 -1.157 0.271 34.077  -4.272   1.000   -1.708   -0.607
## dadas_rscale                 -1.123 0.263 34.077  -4.272   1.000   -1.657   -0.589
## abs_diff                      0.016 0.014 40.150   1.152   0.128   -0.012    0.043
## abs_sum                       0.410 0.095 33.745   4.326   0.000    0.217    0.602
## abs_diff_bscale               0.046 0.040 40.150   1.152   0.128   -0.035    0.126
## abs_sum_bscale                1.203 0.278 33.745   4.326   0.000    0.638    1.768
## abs_diff_rscale               0.083 0.043 41.277   1.936   0.030   -0.004    0.169
## abs_sum_rscale                1.205 0.279 33.745   4.326   0.000    0.639    1.772
```

``` r
round(ddsc_mod2_GEI_cntry_year$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.008 -0.674  4.084  1.000  0.043
```

### Comparisons across models


``` r
# full multilevel model (individuals as entries)
round(ddsc_mod2_GEI$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df t.ratio p.value ci.lower ci.upper
## r_xy1             -0.6360 0.1436 33.1034 -4.4302  0.0001  -0.9280  -0.3439
## r_xy2             -0.5722 0.1281 33.1057 -4.4669  0.0001  -0.8328  -0.3116
## b_11              -0.6206 0.1401 33.1034 -4.4302  0.0001  -0.9056  -0.3356
## b_21              -0.5860 0.1312 33.1057 -4.4669  0.0001  -0.8529  -0.3191
## main_effect       -0.2098 0.0468 33.0553 -4.4821  0.0001  -0.3050  -0.1146
## moderator_effect   0.1774 0.0116 32.7861 15.2608  0.0000   0.1538   0.2011
## interaction        0.0120 0.0120 35.1390  1.0008  0.3238  -0.0124   0.0364
## q_b11_b21         -0.0545     NA      NA      NA      NA       NA       NA
```

``` r
# country-level structural path model (country-gender means as entries)
round(ddsc_sem_GEI$results[c("r_xy1","r_xy2","b_11","b_21","q_b11_b21","diff_b11_b21"),],4)
```

```
##                  est     se       z pvalue ci.lower ci.upper
## r_xy1        -0.5863 0.1410 -4.1571 0.0000  -0.8627  -0.3099
## r_xy2        -0.5952 0.1399 -4.2550 0.0000  -0.8694  -0.3210
## b_11         -0.6002 0.1444 -4.1571 0.0000  -0.8832  -0.3172
## b_21         -0.5807 0.1365 -4.2550 0.0000  -0.8482  -0.3132
## q_b11_b21    -0.0300 0.0527 -0.5700 0.5687  -0.1332   0.0732
## diff_b11_b21 -0.0196 0.0327 -0.5970 0.5505  -0.0837   0.0446
```

``` r
# multilevel model (country-gender means as entries)
round(ddsc_mod2_GEI_ri$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df t.ratio p.value ci.lower ci.upper
## r_xy1             -0.5863 0.1416 31.8537 -4.1411  0.0002  -0.8747  -0.2978
## r_xy2             -0.5952 0.1486 31.8537 -4.0062  0.0003  -0.8979  -0.2925
## b_11              -0.6004 0.1450 31.8537 -4.1411  0.0002  -0.8958  -0.3050
## b_21              -0.5808 0.1450 31.8537 -4.0062  0.0003  -0.8762  -0.2855
## main_effect       -0.2054 0.0501 31.0000 -4.1016  0.0003  -0.3075  -0.1032
## moderator_effect   0.1708 0.0116 31.0000 14.7589  0.0000   0.1472   0.1944
## interaction        0.0068 0.0118 31.0000  0.5786  0.5670  -0.0172   0.0308
## q_b11_b21         -0.0300     NA      NA      NA      NA       NA       NA
```

``` r
# multilevel model (country-year-gender means as entries)
round(ddsc_mod2_GEI_cntry_year$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df t.ratio p.value ci.lower ci.upper
## r_xy1             -0.6442 0.1499 34.0161 -4.2978  0.0001  -0.9487  -0.3396
## r_xy2             -0.5613 0.1314 34.0774 -4.2716  0.0001  -0.8284  -0.2943
## b_11              -0.6244 0.1453 34.0161 -4.2978  0.0001  -0.9197  -0.3292
## b_21              -0.5785 0.1354 34.0774 -4.2716  0.0001  -0.8538  -0.3033
## main_effect       -0.2048 0.0474 33.7453 -4.3262  0.0001  -0.3011  -0.1086
## moderator_effect   0.1785 0.0120 34.1692 14.8664  0.0000   0.1541   0.2029
## interaction        0.0156 0.0136 40.1503  1.1524  0.2560  -0.0118   0.0430
## q_b11_b21         -0.0719     NA      NA      NA      NA       NA       NA
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
## Time difference of 1.609385 hours
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
```

```
## Warning in atanh(mod2_GEI_boot_est$b11): NaNs produced
```

``` r
mod2_GEI_boot_est$q<-atanh(mod2_GEI_boot_est$r_xy1)-atanh(mod2_GEI_boot_est$r_xy2)
```

```
## Warning in atanh(mod2_GEI_boot_est$r_xy1): NaNs produced
```

``` r
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
##                    Estimate         SE        2.5%       97.5%
## X.Intercept.     0.05067809 0.04542722 -0.04050041  0.14186737
## gndr.c           0.17800460 0.01166631  0.15552525  0.20081185
## gei.z.cm        -0.20880310 0.04676026 -0.29683578 -0.11232419
## gndr.c.gei.z.cm  0.01151192 0.01227295 -0.01219147  0.03446830
## w11             -0.21455907 0.04838349 -0.30756302 -0.11685489
## w21             -0.20304714 0.04590636 -0.29103321 -0.10993519
## b11             -0.61707005 0.13915050 -0.88454863 -0.33607367
## b21             -0.58396185 0.13202630 -0.83700904 -0.31617267
## r_xy1           -0.63232725 0.14259103 -0.90641930 -0.34438317
## r_xy2           -0.57020361 0.12891574 -0.81728896 -0.30872359
## q_b             -0.06454771 0.09353840 -0.24867970  0.05769735
## q               -0.11609586 0.09501390 -0.35583863  0.01468152
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
## [1] -0.06454771
## 
## $se
## [1] 0.0935384
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
## [1] 0.3790132
## 
## $p_low
## [1] 0.352339
## 
## $z_high
## [1] -1.759146
## 
## $p_high
## [1] 0.03927635
## 
## $ci_level
## [1] 0.9
## 
## $ci_lower
## [1] -0.2184047
## 
## $ci_upper
## [1] 0.08930926
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
## [1] -0.1160959
## 
## $se
## [1] 0.0950139
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
## [1] -0.1694053
## 
## $p_low
## [1] 0.5672611
## 
## $z_high
## [1] -2.274361
## 
## $p_high
## [1] 0.01147216
## 
## $ci_level
## [1] 0.9
## 
## $ci_lower
## [1] -0.2723798
## 
## $ci_upper
## [1] 0.04018809
## 
## $equivalent
## [1] FALSE
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
  xlab("Gender Equality Index (Average 2002-2023)")+
  #ylim(c(min.y.pred,max.y.pred))+
  #xlim(c(0.60,1.00))+
  ylab("Value achievement (Average 2002-2023)")+
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
  xlab("Gender Equality Index (Average 2002-2023)")+
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

![](Analysis_achievement_files/figure-html/unnamed-chunk-23-1.png)<!-- -->

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
## 1079668.6 1079755.0 -539826.3 1079652.6    363844 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -5.9839 -0.6669  0.0224  0.6613  5.0827 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr  
##  cntry    (Intercept) 0.059178 0.24326        
##           gndr.c      0.004012 0.06334  -0.32 
##  Residual             0.998413 0.99921        
## Number of obs: 363852, groups:  cntry, 34
## 
## Fixed effects:
##                  Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)       0.06223    0.04177 34.03451   1.490    0.145    
## gndr.c            0.17046    0.01152 32.79634  14.795 4.53e-16 ***
## gggi.z.cm        -0.23779    0.04242 34.09681  -5.606 2.77e-06 ***
## gndr.c:gggi.z.cm  0.01004    0.01195 35.47748   0.840    0.406    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c ggg.z.
## gndr.c      -0.302              
## gggi.z.cm    0.000  0.000       
## gndr.c:gg..  0.000 -0.014 -0.296
```

``` r
getFE(mod2_GGGI,round=3)
```

```
##                    Est.    SE     df      t     p     LL     UL
## (Intercept)       0.062 0.042 34.035  1.490 0.145 -0.023  0.147
## gndr.c            0.170 0.012 32.796 14.795 0.000  0.147  0.194
## gggi.z.cm        -0.238 0.042 34.097 -5.606 0.000 -0.324 -0.152
## gndr.c:gggi.z.cm  0.010 0.012 35.477  0.840 0.406 -0.014  0.034
```

``` r
getVC(mod2_GGGI)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.24 0.06
## 2    cntry      gndr.c   <NA>  0.06 0.00
## 3    cntry (Intercept) gndr.c -0.32 0.00
## 4 Residual        <NA>   <NA>  1.00 1.00
```

``` r
r2mlm(mod2_GGGI,bargraph = F)
```

```
## $Decompositions
##                        total
## fixed           0.0448037219
## slope variation 0.0008989754
## mean variation  0.0537405083
## sigma2          0.9005567945
## 
## $R2s
##            total
## f   0.0448037219
## v   0.0008989754
## m   0.0537405083
## fv  0.0457026973
## fvm 0.0994432055
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
## Time difference of 29.26521 secs
```

``` r
round(ddsc_mod2_GGGI$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.118        0.343        1.003     1.121 0.105   7802.647 0.999   0.999
## 2        0.5         0.104        0.323        1.003     1.107 0.094   6678.029 0.998   0.999
```

``` r
round(ddsc_mod2_GGGI$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gggi.z.cm
## means_y1           0.109 0.343    1.000           1.000    0.983           0.983    -0.699
## means_y1_scaled    0.310 0.974    1.000           1.000    0.983           0.983    -0.699
## means_y2          -0.052 0.361    0.983           0.983    1.000           1.000    -0.692
## means_y2_scaled   -0.147 1.025    0.983           0.983    1.000           1.000    -0.692
## gggi.z.cm          0.000 1.000   -0.699          -0.699   -0.692          -0.692     1.000
## gggi.z.cm_scaled   0.000 1.000   -0.699          -0.699   -0.692          -0.692     1.000
## diff_score         0.161 0.068   -0.174          -0.174   -0.354          -0.354     0.150
## diff_score_scaled  0.457 0.193   -0.174          -0.174   -0.354          -0.354     0.150
##                   gggi.z.cm_scaled diff_score diff_score_scaled
## means_y1                    -0.699     -0.174            -0.174
## means_y1_scaled             -0.699     -0.174            -0.174
## means_y2                    -0.692     -0.354            -0.354
## means_y2_scaled             -0.692     -0.354            -0.354
## gggi.z.cm                    1.000      0.150             0.150
## gggi.z.cm_scaled             1.000      0.150             0.150
## diff_score                   0.150      1.000             1.000
## diff_score_scaled            0.150      1.000             1.000
```

``` r
round(ddsc_mod2_GGGI$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                      -0.148 0.176 35.477  -0.840   0.406   -0.504    0.209
## w_11                         -0.243 0.045 34.177  -5.450   0.000   -0.333   -0.152
## w_21                         -0.233 0.041 34.171  -5.671   0.000   -0.316   -0.149
## r_xy1                        -0.708 0.130 34.177  -5.450   0.000   -0.971   -0.444
## r_xy2                        -0.644 0.114 34.171  -5.671   0.000   -0.875   -0.413
## b_11                         -0.689 0.126 34.177  -5.450   0.000   -0.946   -0.432
## b_21                         -0.661 0.117 34.171  -5.671   0.000   -0.898   -0.424
## main_effect                  -0.238 0.042 34.097  -5.606   0.000   -0.324   -0.152
## moderator_effect              0.170 0.012 32.796  14.795   0.000    0.147    0.194
## interaction                   0.010 0.012 35.477   0.840   0.406   -0.014    0.034
## q_b11_b21                    -0.052    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                  -0.117    NA     NA      NA      NA       NA       NA
## cross_over_point            -16.978    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.228 0.041 34.377  -5.621   0.000   -0.310   -0.145
## interaction_vs_main_bscale   -0.647 0.115 34.377  -5.621   0.000   -0.880   -0.413
## interaction_vs_main_rscale   -0.613 0.109 34.398  -5.611   0.000   -0.834   -0.391
## dadas                        -0.466 0.082 34.171  -5.671   1.000   -0.632   -0.299
## dadas_bscale                 -1.322 0.233 34.171  -5.671   1.000   -1.795   -0.848
## dadas_rscale                 -1.289 0.227 34.171  -5.671   1.000   -1.750   -0.827
## abs_diff                      0.010 0.012 35.477   0.840   0.203   -0.014    0.034
## abs_sum                       0.476 0.085 34.097   5.606   0.000    0.303    0.648
## abs_diff_bscale               0.029 0.034 35.477   0.840   0.203   -0.040    0.097
## abs_sum_bscale                1.350 0.241 34.097   5.606   0.000    0.861    1.840
## abs_diff_rscale               0.063 0.036 35.870   1.745   0.045   -0.010    0.137
## abs_sum_rscale                1.352 0.241 34.097   5.603   0.000    0.862    1.842
```

``` r
round(ddsc_mod2_GGGI$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.007 -0.321  3.304  1.000  0.069
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
##                                     est     se       z pvalue ci.lower ci.upper
## r_xy1_y2                         -0.150  0.170  -0.885  0.376   -0.482    0.182
## r_xy1                            -0.692  0.124  -5.596  0.000   -0.935   -0.450
## r_xy2                            -0.699  0.123  -5.703  0.000   -0.940   -0.459
## b_11                             -0.710  0.127  -5.596  0.000   -0.959   -0.461
## b_21                             -0.681  0.119  -5.703  0.000   -0.915   -0.447
## b_10                             -0.147  0.125  -1.174  0.240   -0.392    0.098
## b_20                              0.310  0.118   2.637  0.008    0.080    0.541
## res_cov_y1_y2                     0.483  0.119   4.052  0.000    0.249    0.717
## diff_b10_b20                     -0.457  0.032 -14.163  0.000   -0.520   -0.394
## diff_b11_b21                     -0.029  0.033  -0.885  0.376   -0.093    0.035
## diff_rxy1_rxy2                    0.007  0.032   0.213  0.831   -0.056    0.069
## q_b11_b21                        -0.056  0.070  -0.801  0.423   -0.194    0.081
## q_rxy1_rxy2                       0.013  0.062   0.213  0.831   -0.108    0.135
## cross_over_point                -15.767 17.851  -0.883  0.377  -50.754   19.221
## sum_b11_b21                      -1.391  0.244  -5.695  0.000   -1.870   -0.912
## main_effect                      -0.695  0.122  -5.695  0.000   -0.935   -0.456
## interaction_vs_main_effect       -0.666  0.119  -5.604  0.000   -0.900   -0.433
## diff_abs_b11_abs_b21              0.029  0.033   0.885  0.376   -0.035    0.093
## abs_diff_b11_b21                  0.029  0.033   0.885  0.188   -0.035    0.093
## abs_sum_b11_b21                   1.391  0.244   5.695  0.000    0.912    1.870
## dadas                            -1.362  0.239  -5.703  1.000   -1.830   -0.894
## q_r_equivalence                  -0.087  0.062  -1.402  0.080       NA       NA
## q_b_equivalence                  -0.044  0.070  -0.625  0.266       NA       NA
## cross_over_point_equivalence     15.767 17.851   0.883  0.811       NA       NA
## cross_over_point_minimal_effect  15.767 17.851   0.883  0.189       NA       NA
```

``` r
round(ddsc_sem_GGGI$variance_test,3)
```

```
##             est    se       z pvalue ci.lower ci.upper
## cov_y1y2  0.952 0.233   4.087  0.000    0.496    1.409
## var_y1    1.020 0.248   4.123  0.000    0.535    1.506
## var_y2    0.921 0.223   4.123  0.000    0.483    1.358
## var_diff  0.100 0.066   1.507  0.132   -0.030    0.230
## var_ratio 1.108 0.071  15.717  0.000    0.970    1.247
## cor_y1y2  0.983 0.006 166.517  0.000    0.971    0.994
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
## r_xy1y2                      -0.150 0.175 32.000  -0.859   0.397   -0.506    0.206
## w_11                         -0.250 0.045 33.151  -5.590   0.000   -0.341   -0.159
## w_21                         -0.240 0.045 33.151  -5.362   0.000   -0.331   -0.149
## r_xy1                        -0.692 0.124 33.151  -5.590   0.000   -0.944   -0.440
## r_xy2                        -0.699 0.130 33.151  -5.362   0.000   -0.964   -0.434
## b_11                         -0.710 0.127 33.151  -5.590   0.000   -0.969   -0.452
## b_21                         -0.681 0.127 33.151  -5.362   0.000   -0.940   -0.423
## main_effect                  -0.245 0.044 32.000  -5.525   0.000   -0.335   -0.155
## moderator_effect              0.161 0.012 32.000  13.740   0.000    0.137    0.185
## interaction                   0.010 0.012 32.000   0.859   0.397   -0.014    0.034
## q_b11_b21                    -0.056    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.013    NA     NA      NA      NA       NA       NA
## cross_over_point            -15.767    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.235 0.046 36.580  -5.114   0.000   -0.328   -0.142
## interaction_vs_main_bscale   -0.667 0.130 36.580  -5.114   0.000   -0.931   -0.402
## interaction_vs_main_rscale   -0.703 0.137 36.253  -5.135   0.000   -0.980   -0.425
## dadas                        -0.480 0.089 33.151  -5.362   1.000   -0.662   -0.298
## dadas_bscale                 -1.362 0.254 33.151  -5.362   1.000   -1.879   -0.846
## dadas_rscale                 -1.398 0.261 33.151  -5.362   1.000   -1.929   -0.868
## abs_diff                      0.010 0.012 32.000   0.859   0.198   -0.014    0.034
## abs_sum                       0.490 0.089 32.000   5.525   0.000    0.309    0.671
## abs_diff_bscale               0.029 0.034 32.000   0.859   0.198   -0.040    0.098
## abs_sum_bscale                1.391 0.252 32.000   5.525   0.000    0.878    1.904
## abs_diff_rscale              -0.007 0.034 34.352  -0.198   0.578   -0.077    0.063
## abs_sum_rscale                1.392 0.252 32.001   5.522   0.000    0.878    1.905
```

``` r
# country-time multilevel model


mod2_GGGI_cntry_year<-
  lmer(ach.z.wt~gndr.c+gggi.z.cm+gndr.c:gggi.z.cm+
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
## Formula: ach.z.wt ~ gndr.c + gggi.z.cm + gndr.c:gggi.z.cm + (gndr.c |      cntry)
##    Data: diff_dat_cntry_year
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##    -481.1    -449.1     248.5    -497.1       392 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -3.4942 -0.5982  0.0026  0.6081  3.3547 
## 
## Random effects:
##  Groups   Name        Variance  Std.Dev. Corr  
##  cntry    (Intercept) 0.0567162 0.23815        
##           gndr.c      0.0005278 0.02297  -1.00 
##  Residual             0.0121049 0.11002        
## Number of obs: 400, groups:  cntry, 34
## 
## Fixed effects:
##                   Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)        0.05795    0.04134  34.28818   1.402    0.170    
## gndr.c             0.17229    0.01180 126.85719  14.601  < 2e-16 ***
## gggi.z.cm         -0.23887    0.04214  34.84903  -5.668 2.15e-06 ***
## gndr.c:gggi.z.cm   0.01210    0.01260 144.16321   0.960    0.339    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c ggg.z.
## gndr.c      -0.332              
## gggi.z.cm   -0.007  0.001       
## gndr.c:gg..  0.001 -0.129 -0.316
## optimizer (bobyqa) convergence code: 0 (OK)
## boundary (singular) fit: see help('isSingular')
```

``` r
getFE(mod2_GGGI_cntry_year,round=3)
```

```
##                    Est.    SE      df      t     p     LL     UL
## (Intercept)       0.058 0.041  34.288  1.402 0.170 -0.026  0.142
## gndr.c            0.172 0.012 126.857 14.601 0.000  0.149  0.196
## gggi.z.cm        -0.239 0.042  34.849 -5.668 0.000 -0.324 -0.153
## gndr.c:gggi.z.cm  0.012 0.013 144.163  0.960 0.339 -0.013  0.037
```

``` r
getVC(mod2_GGGI_cntry_year)
```

```
##        grp        var1   var2 sdcor  vcov
## 1    cntry (Intercept)   <NA>  0.24  0.06
## 2    cntry      gndr.c   <NA>  0.02  0.00
## 3    cntry (Intercept) gndr.c -1.00 -0.01
## 4 Residual        <NA>   <NA>  0.11  0.01
```

``` r
r2mlm(mod2_GGGI,bargraph = F)
```

```
## $Decompositions
##                        total
## fixed           0.0448037219
## slope variation 0.0008989754
## mean variation  0.0537405083
## sigma2          0.9005567945
## 
## $R2s
##            total
## f   0.0448037219
## v   0.0008989754
## m   0.0537405083
## fv  0.0457026973
## fvm 0.0994432055
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
## 1       -0.5         0.116        0.340        0.014     0.129 0.895      8.029 0.998   0.986
## 2        0.5         0.100        0.316        0.014     0.114 0.880      8.029 0.998   0.983
```

``` r
round(ddsc_mod2_GGGI_cntry_year$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gggi.z.cm
## means_y1           0.143 0.333    1.000           1.000    0.981           0.981    -0.707
## means_y1_scaled    0.416 0.970    1.000           1.000    0.981           0.981    -0.707
## means_y2          -0.025 0.353    0.981           0.981    1.000           1.000    -0.689
## means_y2_scaled   -0.073 1.029    0.981           0.981    1.000           1.000    -0.689
## gggi.z.cm          0.000 1.000   -0.707          -0.707   -0.689          -0.689     1.000
## gggi.z.cm_scaled   0.000 1.000   -0.707          -0.707   -0.689          -0.689     1.000
## diff_score         0.168 0.070   -0.193          -0.193   -0.380          -0.380     0.114
## diff_score_scaled  0.490 0.204   -0.193          -0.193   -0.380          -0.380     0.114
##                   gggi.z.cm_scaled diff_score diff_score_scaled
## means_y1                    -0.707     -0.193            -0.193
## means_y1_scaled             -0.707     -0.193            -0.193
## means_y2                    -0.689     -0.380            -0.380
## means_y2_scaled             -0.689     -0.380            -0.380
## gggi.z.cm                    1.000      0.114             0.114
## gggi.z.cm_scaled             1.000      0.114             0.114
## diff_score                   0.114      1.000             1.000
## diff_score_scaled            0.114      1.000             1.000
```

``` r
round(ddsc_mod2_GGGI_cntry_year$results,3)
```

```
##                            estimate    SE      df t.ratio p.value ci.lower ci.upper
## r_xy1y2                      -0.173 0.180 144.163  -0.960   0.339   -0.529    0.183
## w_11                         -0.245 0.045  35.058  -5.499   0.000   -0.335   -0.155
## w_21                         -0.233 0.041  35.135  -5.735   0.000   -0.315   -0.150
## r_xy1                        -0.735 0.134  35.058  -5.499   0.000   -1.007   -0.464
## r_xy2                        -0.659 0.115  35.135  -5.735   0.000   -0.892   -0.426
## b_11                         -0.714 0.130  35.058  -5.499   0.000   -0.977   -0.450
## b_21                         -0.678 0.118  35.135  -5.735   0.000   -0.918   -0.438
## main_effect                  -0.239 0.042  34.849  -5.668   0.000   -0.324   -0.153
## moderator_effect              0.172 0.012 126.857  14.601   0.000    0.149    0.196
## interaction                   0.012 0.013 144.163   0.960   0.339   -0.013    0.037
## q_b11_b21                    -0.068    NA      NA      NA      NA       NA       NA
## q_rxy1_rxy2                  -0.149    NA      NA      NA      NA       NA       NA
## cross_over_point            -14.237    NA      NA      NA      NA       NA       NA
## interaction_vs_main          -0.227 0.040  36.162  -5.671   0.000   -0.308   -0.146
## interaction_vs_main_bscale   -0.661 0.117  36.162  -5.671   0.000   -0.897   -0.424
## interaction_vs_main_rscale   -0.621 0.110  36.311  -5.657   0.000   -0.843   -0.398
## dadas                        -0.466 0.081  35.135  -5.735   1.000   -0.630   -0.301
## dadas_bscale                 -1.357 0.237  35.135  -5.735   1.000   -1.837   -0.876
## dadas_rscale                 -1.318 0.230  35.135  -5.735   1.000   -1.784   -0.851
## abs_diff                      0.012 0.013 144.163   0.960   0.169   -0.013    0.037
## abs_sum                       0.478 0.084  34.849   5.668   0.000    0.307    0.649
## abs_diff_bscale               0.035 0.037 144.163   0.960   0.169   -0.037    0.108
## abs_sum_bscale                1.392 0.246  34.849   5.668   0.000    0.893    1.890
## abs_diff_rscale               0.076 0.040  85.911   1.928   0.029   -0.002    0.155
## abs_sum_rscale                1.394 0.246  34.849   5.664   0.000    0.894    1.894
```

``` r
round(ddsc_mod2_GGGI_cntry_year$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.008 -0.674  4.084  1.000  0.043
```

### Comparisons across models


``` r
# full multilevel model (individuals as entries)
round(ddsc_mod2_GGGI$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df t.ratio p.value ci.lower ci.upper
## r_xy1             -0.7076 0.1298 34.1773 -5.4501  0.0000  -0.9714  -0.4438
## r_xy2             -0.6443 0.1136 34.1711 -5.6711  0.0000  -0.8752  -0.4135
## b_11              -0.6894 0.1265 34.1773 -5.4501  0.0000  -0.9464  -0.4324
## b_21              -0.6609 0.1165 34.1711 -5.6711  0.0000  -0.8977  -0.4241
## main_effect       -0.2378 0.0424 34.0968 -5.6063  0.0000  -0.3240  -0.1516
## moderator_effect   0.1705 0.0115 32.7963 14.7945  0.0000   0.1470   0.1939
## interaction        0.0100 0.0120 35.4775  0.8401  0.4065  -0.0142   0.0343
## q_b11_b21         -0.0524     NA      NA      NA      NA       NA       NA
```

``` r
# country-level structural path model (country-gender means as entries)
round(ddsc_sem_GGGI$results[c("r_xy1","r_xy2","b_11","b_21","q_b11_b21","diff_b11_b21"),],4)
```

```
##                  est     se       z pvalue ci.lower ci.upper
## r_xy1        -0.6924 0.1237 -5.5955 0.0000  -0.9349  -0.4499
## r_xy2        -0.6992 0.1226 -5.7027 0.0000  -0.9395  -0.4589
## b_11         -0.7100 0.1269 -5.5955 0.0000  -0.9587  -0.4613
## b_21         -0.6810 0.1194 -5.7027 0.0000  -0.9150  -0.4469
## q_b11_b21    -0.0562 0.0702 -0.8006 0.4234  -0.1937   0.0813
## diff_b11_b21 -0.0290 0.0328 -0.8850 0.3762  -0.0932   0.0352
```

``` r
# multilevel model (country-gender means as entries)
round(ddsc_mod2_GGGI_ri$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df t.ratio p.value ci.lower ci.upper
## r_xy1             -0.6924 0.1239 33.1506 -5.5905   0.000  -0.9443  -0.4405
## r_xy2             -0.6992 0.1304 33.1506 -5.3622   0.000  -0.9644  -0.4340
## b_11              -0.7102 0.1270 33.1506 -5.5905   0.000  -0.9686  -0.4518
## b_21              -0.6812 0.1270 33.1506 -5.3622   0.000  -0.9396  -0.4228
## main_effect       -0.2450 0.0443 32.0000 -5.5254   0.000  -0.3354  -0.1547
## moderator_effect   0.1610 0.0117 32.0000 13.7400   0.000   0.1371   0.1849
## interaction        0.0102 0.0119 32.0000  0.8585   0.397  -0.0140   0.0344
## q_b11_b21         -0.0562     NA      NA      NA      NA       NA       NA
```

``` r
# multilevel model (country-year-gender means as entries)
round(ddsc_mod2_GGGI_cntry_year$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE       df t.ratio p.value ci.lower ci.upper
## r_xy1             -0.7352 0.1337  35.0581 -5.4990  0.0000  -1.0067  -0.4638
## r_xy2             -0.6588 0.1149  35.1345 -5.7352  0.0000  -0.8920  -0.4256
## b_11              -0.7135 0.1298  35.0581 -5.4990  0.0000  -0.9769  -0.4501
## b_21              -0.6783 0.1183  35.1345 -5.7352  0.0000  -0.9183  -0.4382
## main_effect       -0.2389 0.0421  34.8490 -5.6679  0.0000  -0.3244  -0.1533
## moderator_effect   0.1723 0.0118 126.8572 14.6006  0.0000   0.1489   0.1956
## interaction        0.0121 0.0126 144.1632  0.9601  0.3386  -0.0128   0.0370
## q_b11_b21         -0.0684     NA       NA      NA      NA       NA       NA
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
## Time difference of 1.125385 hours
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
```

```
## Warning in atanh(mod2_GGGI_boot_est$b11): NaNs produced
```

```
## Warning in atanh(mod2_GGGI_boot_est$b21): NaNs produced
```

``` r
mod2_GGGI_boot_est$q<-atanh(mod2_GGGI_boot_est$r_xy1)-atanh(mod2_GGGI_boot_est$r_xy2)
```

```
## Warning in atanh(mod2_GGGI_boot_est$r_xy1): NaNs produced
```

``` r
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
##                     Estimate         SE        2.5%       97.5%
## X.Intercept.      0.06129884 0.04214789 -0.01583780  0.14478525
## gndr.c            0.17109785 0.01228907  0.14692861  0.19405526
## gggi.z.cm        -0.23658568 0.04120327 -0.32207358 -0.15739765
## gndr.c.gggi.z.cm  0.01014073 0.01235281 -0.01493723  0.03628962
## w11              -0.24165604 0.04332126 -0.32976708 -0.15696176
## w21              -0.23151531 0.03993723 -0.31350724 -0.15412556
## b11              -0.68614022 0.12300318 -0.93631616 -0.44566556
## b21              -0.65734738 0.11339481 -0.89014917 -0.43761267
## r_xy1            -0.70425830 0.12625118 -0.96104033 -0.45743371
## r_xy2            -0.64086031 0.11055073 -0.86782315 -0.42663681
## q_b              -0.07089162 0.12382178 -0.32224286  0.06619446
## q                -0.14433043 0.14204683 -0.48876936  0.01662839
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
## [1] -0.07089162
## 
## $se
## [1] 0.1238218
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
## [1] 0.2350829
## 
## $p_low
## [1] 0.4070722
## 
## $z_high
## [1] -1.380142
## 
## $p_high
## [1] 0.08377148
## 
## $ci_level
## [1] 0.9
## 
## $ci_lower
## [1] -0.2745603
## 
## $ci_upper
## [1] 0.1327771
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
## [1] -0.1443304
## 
## $se
## [1] 0.1420468
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
## [1] -0.3120832
## 
## $p_low
## [1] 0.6225113
## 
## $z_high
## [1] -1.72007
## 
## $p_high
## [1] 0.0427099
## 
## $ci_level
## [1] 0.9
## 
## $ci_lower
## [1] -0.3779767
## 
## $ci_upper
## [1] 0.08931582
## 
## $equivalent
## [1] FALSE
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
  xlab("Global Gender Gap Index (Average 2002-2023)")+
  #ylim(c(min.y.pred,max.y.pred))+
  #xlim(c(0.60,1.00))+
  ylab("Value achievement (Average 2002-2023)")+
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
  xlab("Global Gender Gap Index (Average 2002-2023)")+
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

![](Analysis_achievement_files/figure-html/unnamed-chunk-29-1.png)<!-- -->

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
## 1462319.3 1462408.1 -731151.6 1462303.3    492335 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -5.9102 -0.6699  0.0261  0.6594  5.0736 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr  
##  cntry    (Intercept) 0.105111 0.32421        
##           gndr.c      0.004027 0.06346  -0.35 
##  Residual             1.002867 1.00143        
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##                  Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)      0.046731   0.055626 34.000049   0.840    0.407    
## gndr.c           0.176929   0.011360 33.698867  15.575   <2e-16 ***
## gdi.z.cm        -0.071867   0.056477 34.034735  -1.273    0.212    
## gndr.c:gdi.z.cm -0.007687   0.011793 36.680133  -0.652    0.519    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c gd.z.c
## gndr.c      -0.339              
## gdi.z.cm     0.000  0.000       
## gndr.c:gd..  0.000 -0.011 -0.332
```

``` r
getFE(mod2_GDI,round=3)
```

```
##                   Est.    SE     df      t     p     LL    UL
## (Intercept)      0.047 0.056 34.000  0.840 0.407 -0.066 0.160
## gndr.c           0.177 0.011 33.699 15.575 0.000  0.154 0.200
## gdi.z.cm        -0.072 0.056 34.035 -1.273 0.212 -0.187 0.043
## gndr.c:gdi.z.cm -0.008 0.012 36.680 -0.652 0.519 -0.032 0.016
```

``` r
getVC(mod2_GDI)
```

```
##        grp        var1   var2 sdcor  vcov
## 1    cntry (Intercept)   <NA>  0.32  0.11
## 2    cntry      gndr.c   <NA>  0.06  0.00
## 3    cntry (Intercept) gndr.c -0.35 -0.01
## 4 Residual        <NA>   <NA>  1.00  1.00
```

``` r
r2mlm(mod2_GDI,bargraph = F)
```

```
## $Decompositions
##                        total
## fixed           0.0105152932
## slope variation 0.0008923144
## mean variation  0.0942478753
## sigma2          0.8943445171
## 
## $R2s
##            total
## f   0.0105152932
## v   0.0008923144
## m   0.0942478753
## fv  0.0114076077
## fvm 0.1056554829
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
## Time difference of 29.93538 secs
```

``` r
round(ddsc_mod2_GDI$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.118        0.343        1.003     1.121 0.105   7802.647 0.999   0.999
## 2        0.5         0.104        0.323        1.003     1.107 0.094   6678.029 0.998   0.999
```

``` r
round(ddsc_mod2_GDI$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gdi.z.cm gdi.z.cm_scaled
## means_y1           0.098 0.336    1.000           1.000    0.983           0.983   -0.215          -0.215
## means_y1_scaled    0.284 0.976    1.000           1.000    0.983           0.983   -0.215          -0.215
## means_y2          -0.073 0.352    0.983           0.983    1.000           1.000   -0.213          -0.213
## means_y2_scaled   -0.212 1.024    0.983           0.983    1.000           1.000   -0.213          -0.213
## gdi.z.cm           0.000 1.000   -0.215          -0.215   -0.213          -0.213    1.000           1.000
## gdi.z.cm_scaled    0.000 1.000   -0.215          -0.215   -0.213          -0.213    1.000           1.000
## diff_score         0.171 0.065   -0.163          -0.163   -0.340          -0.340    0.044           0.044
## diff_score_scaled  0.496 0.188   -0.163          -0.163   -0.340          -0.340    0.044           0.044
##                   diff_score diff_score_scaled
## means_y1              -0.163            -0.163
## means_y1_scaled       -0.163            -0.163
## means_y2              -0.340            -0.340
## means_y2_scaled       -0.340            -0.340
## gdi.z.cm               0.044             0.044
## gdi.z.cm_scaled        0.044             0.044
## diff_score             1.000             1.000
## diff_score_scaled      1.000             1.000
```

``` r
round(ddsc_mod2_GDI$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                       0.119 0.182 36.680   0.652   0.519   -0.250    0.488
## w_11                         -0.068 0.059 34.063  -1.159   0.255   -0.187    0.051
## w_21                         -0.076 0.055 34.069  -1.381   0.176   -0.187    0.036
## r_xy1                        -0.203 0.175 34.063  -1.159   0.255   -0.558    0.153
## r_xy2                        -0.215 0.156 34.069  -1.381   0.176   -0.531    0.101
## b_11                         -0.198 0.171 34.063  -1.159   0.255   -0.544    0.149
## b_21                         -0.220 0.159 34.069  -1.381   0.176   -0.544    0.104
## main_effect                  -0.072 0.056 34.035  -1.273   0.212   -0.187    0.043
## moderator_effect              0.177 0.011 33.699  15.575   0.000    0.154    0.200
## interaction                  -0.008 0.012 36.680  -0.652   0.519   -0.032    0.016
## q_b11_b21                     0.023    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.013    NA     NA      NA      NA       NA       NA
## cross_over_point             23.017    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.064 0.061 34.150  -1.045   0.303   -0.189    0.061
## interaction_vs_main_bscale   -0.187 0.179 34.150  -1.045   0.303   -0.549    0.176
## interaction_vs_main_rscale   -0.196 0.187 34.143  -1.053   0.300   -0.575    0.183
## dadas                        -0.136 0.117 34.063  -1.159   0.873   -0.375    0.103
## dadas_bscale                 -0.395 0.341 34.063  -1.159   0.873   -1.089    0.298
## dadas_rscale                 -0.405 0.350 34.063  -1.159   0.873   -1.116    0.305
## abs_diff                      0.008 0.012 36.680   0.652   0.259   -0.016    0.032
## abs_sum                       0.144 0.113 34.035   1.273   0.106   -0.086    0.373
## abs_diff_bscale               0.022 0.034 36.680   0.652   0.259   -0.047    0.092
## abs_sum_bscale                0.418 0.328 34.035   1.273   0.106   -0.249    1.085
## abs_diff_rscale               0.012 0.038 36.726   0.329   0.372   -0.064    0.089
## abs_sum_rscale                0.418 0.329 34.035   1.270   0.106   -0.251    1.086
```

``` r
round(ddsc_mod2_GDI$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.007 -0.321  3.304  1.000  0.069
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
##                                     est      se       z pvalue ci.lower ci.upper
## r_xy1_y2                         -0.044   0.171  -0.255  0.799   -0.379    0.292
## r_xy1                            -0.213   0.168  -1.274  0.203   -0.542    0.115
## r_xy2                            -0.215   0.167  -1.286  0.198   -0.544    0.113
## b_11                             -0.218   0.171  -1.274  0.203   -0.555    0.118
## b_21                             -0.210   0.163  -1.286  0.198   -0.531    0.110
## b_10                             -0.212   0.169  -1.257  0.209   -0.544    0.119
## b_20                              0.284   0.161   1.762  0.078   -0.032    0.599
## res_cov_y1_y2                     0.909   0.222   4.087  0.000    0.473    1.345
## diff_b10_b20                     -0.496   0.032 -15.615  0.000   -0.558   -0.434
## diff_b11_b21                     -0.008   0.032  -0.255  0.799   -0.071    0.055
## diff_rxy1_rxy2                    0.002   0.031   0.065  0.948   -0.059    0.063
## q_b11_b21                        -0.009   0.034  -0.253  0.800   -0.075    0.058
## q_rxy1_rxy2                       0.002   0.033   0.065  0.948   -0.062    0.066
## cross_over_point                -60.426 237.377  -0.255  0.799 -525.677  404.825
## sum_b11_b21                      -0.429   0.333  -1.286  0.199   -1.082    0.225
## main_effect                      -0.214   0.167  -1.286  0.199   -0.541    0.112
## interaction_vs_main_effect       -0.206   0.162  -1.275  0.202   -0.523    0.111
## diff_abs_b11_abs_b21              0.008   0.032   0.255  0.799   -0.055    0.071
## abs_diff_b11_b21                  0.008   0.032   0.255  0.400   -0.055    0.071
## abs_sum_b11_b21                   0.429   0.333   1.286  0.099   -0.225    1.082
## dadas                            -0.420   0.327  -1.286  0.901   -1.061    0.220
## q_r_equivalence                  -0.098   0.033  -2.988  0.001       NA       NA
## q_b_equivalence                  -0.091   0.034  -2.690  0.004       NA       NA
## cross_over_point_equivalence     60.426 237.377   0.255  0.600       NA       NA
## cross_over_point_minimal_effect  60.426 237.377   0.255  0.400       NA       NA
```

``` r
round(ddsc_sem_GDI$variance_test,3)
```

```
##             est    se       z pvalue ci.lower ci.upper
## cov_y1y2  0.953 0.233   4.088   0.00    0.496    1.410
## var_y1    1.017 0.247   4.123   0.00    0.533    1.500
## var_y2    0.924 0.224   4.123   0.00    0.485    1.364
## var_diff  0.093 0.064   1.438   0.15   -0.034    0.219
## var_ratio 1.100 0.068  16.067   0.00    0.966    1.234
## cor_y1y2  0.983 0.006 174.148   0.00    0.972    0.994
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
## r_xy1y2                      -0.044 0.177 32.000  -0.247   0.806   -0.403    0.316
## w_11                         -0.075 0.059 32.599  -1.265   0.215   -0.196    0.046
## w_21                         -0.072 0.059 32.599  -1.218   0.232   -0.193    0.049
## r_xy1                        -0.213 0.169 32.599  -1.265   0.215   -0.557    0.130
## r_xy2                        -0.215 0.177 32.599  -1.218   0.232   -0.576    0.145
## b_11                         -0.219 0.173 32.599  -1.265   0.215   -0.570    0.133
## b_21                         -0.210 0.173 32.599  -1.218   0.232   -0.562    0.141
## main_effect                  -0.074 0.059 32.000  -1.247   0.221   -0.194    0.047
## moderator_effect              0.171 0.011 32.000  15.149   0.000    0.148    0.194
## interaction                   0.003 0.011 32.000   0.247   0.806   -0.020    0.026
## q_b11_b21                    -0.009    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.002    NA     NA      NA      NA       NA       NA
## cross_over_point            -60.426    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.071 0.060 34.392  -1.178   0.247   -0.193    0.051
## interaction_vs_main_bscale   -0.206 0.175 34.392  -1.178   0.247   -0.562    0.150
## interaction_vs_main_rscale   -0.216 0.183 34.232  -1.181   0.246   -0.589    0.156
## dadas                        -0.145 0.119 32.599  -1.218   0.884   -0.387    0.097
## dadas_bscale                 -0.421 0.345 32.599  -1.218   0.884   -1.124    0.283
## dadas_rscale                 -0.431 0.354 32.599  -1.218   0.884   -1.151    0.289
## abs_diff                      0.003 0.011 32.000   0.247   0.403   -0.020    0.026
## abs_sum                       0.148 0.118 32.000   1.247   0.111   -0.093    0.388
## abs_diff_bscale               0.008 0.033 32.000   0.247   0.403   -0.060    0.076
## abs_sum_bscale                0.429 0.344 32.000   1.247   0.111   -0.272    1.129
## abs_diff_rscale              -0.002 0.034 35.881  -0.059   0.523   -0.072    0.067
## abs_sum_rscale                0.429 0.344 32.000   1.247   0.111   -0.272    1.130
```

``` r
# country-time multilevel model


mod2_GDI_cntry_year<-
  lmer(ach.z.wt~gndr.c+gdi.z.cm+gndr.c:gdi.z.cm+
                 (gndr.c|cntry),data=diff_dat_cntry_year,REML=F,
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2_GDI_cntry_year)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ach.z.wt ~ gndr.c + gdi.z.cm + gndr.c:gdi.z.cm + (gndr.c | cntry)
##    Data: diff_dat_cntry_year
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##    -617.1    -582.7     316.5    -633.1       538 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -3.6841 -0.6040  0.0020  0.6588  3.7641 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr  
##  cntry    (Intercept) 0.102016 0.31940        
##           gndr.c      0.001007 0.03173  -0.79 
##  Residual             0.013608 0.11665        
## Number of obs: 546, groups:  cntry, 34
## 
## Fixed effects:
##                 Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)      0.04117    0.05506 34.03916   0.748    0.460    
## gndr.c           0.17921    0.01145 31.62280  15.653   <2e-16 ***
## gdi.z.cm        -0.07773    0.05610 34.55634  -1.386    0.175    
## gndr.c:gdi.z.cm -0.01280    0.01360 44.96762  -0.941    0.352    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c gd.z.c
## gndr.c      -0.374              
## gdi.z.cm    -0.002  0.001       
## gndr.c:gd..  0.000 -0.050 -0.319
```

``` r
getFE(mod2_GDI_cntry_year,round=3)
```

```
##                   Est.    SE     df      t     p     LL    UL
## (Intercept)      0.041 0.055 34.039  0.748 0.460 -0.071 0.153
## gndr.c           0.179 0.011 31.623 15.653 0.000  0.156 0.203
## gdi.z.cm        -0.078 0.056 34.556 -1.386 0.175 -0.192 0.036
## gndr.c:gdi.z.cm -0.013 0.014 44.968 -0.941 0.352 -0.040 0.015
```

``` r
getVC(mod2_GDI_cntry_year)
```

```
##        grp        var1   var2 sdcor  vcov
## 1    cntry (Intercept)   <NA>  0.32  0.10
## 2    cntry      gndr.c   <NA>  0.03  0.00
## 3    cntry (Intercept) gndr.c -0.79 -0.01
## 4 Residual        <NA>   <NA>  0.12  0.01
```

``` r
r2mlm(mod2_GDI,bargraph = F)
```

```
## $Decompositions
##                        total
## fixed           0.0105152932
## slope variation 0.0008923144
## mean variation  0.0942478753
## sigma2          0.8943445171
## 
## $R2s
##            total
## f   0.0105152932
## v   0.0008923144
## m   0.0942478753
## fv  0.0114076077
## fvm 0.1056554829
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
## 1       -0.5         0.116        0.340        0.014     0.129 0.895      8.029 0.998   0.986
## 2        0.5         0.100        0.316        0.014     0.114 0.880      8.029 0.998   0.983
```

``` r
round(ddsc_mod2_GDI_cntry_year$descriptives,3)
```

```
##                        M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled gdi.z.cm gdi.z.cm_scaled
## means_y1           0.129 0.326    1.000           1.000    0.982           0.982   -0.252          -0.252
## means_y1_scaled    0.384 0.969    1.000           1.000    0.982           0.982   -0.252          -0.252
## means_y2          -0.045 0.347    0.982           0.982    1.000           1.000   -0.216          -0.216
## means_y2_scaled   -0.135 1.030    0.982           0.982    1.000           1.000   -0.216          -0.216
## gdi.z.cm           0.000 1.000   -0.252          -0.252   -0.216          -0.216    1.000           1.000
## gdi.z.cm_scaled    0.000 1.000   -0.252          -0.252   -0.216          -0.216    1.000           1.000
## diff_score         0.175 0.067   -0.206          -0.206   -0.389          -0.389   -0.109          -0.109
## diff_score_scaled  0.518 0.200   -0.206          -0.206   -0.389          -0.389   -0.109          -0.109
##                   diff_score diff_score_scaled
## means_y1              -0.206            -0.206
## means_y1_scaled       -0.206            -0.206
## means_y2              -0.389            -0.389
## means_y2_scaled       -0.389            -0.389
## gdi.z.cm              -0.109            -0.109
## gdi.z.cm_scaled       -0.109            -0.109
## diff_score             1.000             1.000
## diff_score_scaled      1.000             1.000
```

``` r
round(ddsc_mod2_GDI_cntry_year$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                       0.190 0.202 44.968   0.941   0.352   -0.216    0.596
## w_11                         -0.071 0.059 34.780  -1.217   0.232   -0.190    0.048
## w_21                         -0.084 0.054 34.842  -1.549   0.130   -0.194    0.026
## r_xy1                        -0.218 0.180 34.780  -1.217   0.232   -0.583    0.146
## r_xy2                        -0.243 0.157 34.842  -1.549   0.130   -0.561    0.075
## b_11                         -0.212 0.174 34.780  -1.217   0.232   -0.566    0.142
## b_21                         -0.250 0.161 34.842  -1.549   0.130   -0.578    0.078
## main_effect                  -0.078 0.056 34.556  -1.386   0.175   -0.192    0.036
## moderator_effect              0.179 0.011 31.623  15.653   0.000    0.156    0.203
## interaction                  -0.013 0.014 44.968  -0.941   0.352   -0.040    0.015
## q_b11_b21                     0.040    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                   0.026    NA     NA      NA      NA       NA       NA
## cross_over_point             13.999    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.065 0.062 35.397  -1.051   0.301   -0.190    0.060
## interaction_vs_main_bscale   -0.193 0.184 35.397  -1.051   0.301   -0.565    0.180
## interaction_vs_main_rscale   -0.206 0.194 35.333  -1.064   0.294   -0.600    0.187
## dadas                        -0.143 0.117 34.780  -1.217   0.884   -0.381    0.095
## dadas_bscale                 -0.424 0.348 34.780  -1.217   0.884   -1.131    0.284
## dadas_rscale                 -0.437 0.359 34.780  -1.217   0.884   -1.166    0.292
## abs_diff                      0.013 0.014 44.968   0.941   0.176   -0.015    0.040
## abs_sum                       0.155 0.112 34.556   1.386   0.087   -0.072    0.383
## abs_diff_bscale               0.038 0.040 44.968   0.941   0.176   -0.043    0.119
## abs_sum_bscale                0.462 0.333 34.556   1.386   0.087   -0.215    1.139
## abs_diff_rscale               0.024 0.045 45.512   0.541   0.296   -0.066    0.114
## abs_sum_rscale                0.461 0.334 34.556   1.380   0.088   -0.217    1.140
```

``` r
round(ddsc_mod2_GDI_cntry_year$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.008 -0.674  4.084  1.000  0.043
```

### Comparisons across models


``` r
# full multilevel model (individuals as entries)
round(ddsc_mod2_GDI$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df t.ratio p.value ci.lower ci.upper
## r_xy1             -0.2026 0.1748 34.0626 -1.1589  0.2546  -0.5578   0.1526
## r_xy2             -0.2150 0.1556 34.0686 -1.3815  0.1761  -0.5312   0.1012
## b_11              -0.1977 0.1706 34.0626 -1.1589  0.2546  -0.5445   0.1490
## b_21              -0.2201 0.1593 34.0686 -1.3815  0.1761  -0.5438   0.1037
## main_effect       -0.0719 0.0565 34.0347 -1.2725  0.2118  -0.1866   0.0429
## moderator_effect   0.1769 0.0114 33.6989 15.5746  0.0000   0.1538   0.2000
## interaction       -0.0077 0.0118 36.6801 -0.6518  0.5186  -0.0316   0.0162
## q_b11_b21          0.0234     NA      NA      NA      NA       NA       NA
```

``` r
# country-level structural path model (country-gender means as entries)
round(ddsc_sem_GDI$results[c("r_xy1","r_xy2","b_11","b_21","q_b11_b21","diff_b11_b21"),],4)
```

```
##                  est     se       z pvalue ci.lower ci.upper
## r_xy1        -0.2134 0.1675 -1.2738 0.2027  -0.5418   0.1150
## r_xy2        -0.2154 0.1675 -1.2864 0.1983  -0.5437   0.1128
## b_11         -0.2185 0.1715 -1.2738 0.2027  -0.5546   0.1177
## b_21         -0.2102 0.1634 -1.2864 0.1983  -0.5306   0.1101
## q_b11_b21    -0.0086 0.0340 -0.2533 0.8000  -0.0752   0.0580
## diff_b11_b21 -0.0082 0.0323 -0.2546 0.7990  -0.0714   0.0550
```

``` r
# multilevel model (country-gender means as entries)
round(ddsc_mod2_GDI_ri$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df t.ratio p.value ci.lower ci.upper
## r_xy1             -0.2134 0.1687 32.5988 -1.2652  0.2148  -0.5568   0.1299
## r_xy2             -0.2154 0.1769 32.5988 -1.2176  0.2321  -0.5756   0.1447
## b_11              -0.2185 0.1727 32.5988 -1.2652  0.2148  -0.5701   0.1330
## b_21              -0.2103 0.1727 32.5988 -1.2176  0.2321  -0.5619   0.1413
## main_effect       -0.0738 0.0591 32.0000 -1.2472  0.2214  -0.1942   0.0467
## moderator_effect   0.1707 0.0113 32.0000 15.1489  0.0000   0.1478   0.1937
## interaction        0.0028 0.0114 32.0000  0.2470  0.8065  -0.0205   0.0261
## q_b11_b21         -0.0086     NA      NA      NA      NA       NA       NA
```

``` r
# multilevel model (country-year-gender means as entries)
round(ddsc_mod2_GDI_cntry_year$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df t.ratio p.value ci.lower ci.upper
## r_xy1             -0.2185 0.1796 34.7805 -1.2166  0.2319  -0.5831   0.1462
## r_xy2             -0.2426 0.1566 34.8417 -1.5490  0.1304  -0.5607   0.0754
## b_11              -0.2119 0.1742 34.7805 -1.2166  0.2319  -0.5656   0.1418
## b_21              -0.2499 0.1614 34.8417 -1.5490  0.1304  -0.5776   0.0777
## main_effect       -0.0777 0.0561 34.5563 -1.3855  0.1748  -0.1917   0.0362
## moderator_effect   0.1792 0.0114 31.6228 15.6530  0.0000   0.1559   0.2025
## interaction       -0.0128 0.0136 44.9676 -0.9412  0.3516  -0.0402   0.0146
## q_b11_b21          0.0402     NA      NA      NA      NA       NA       NA
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
## Time difference of 1.875458 hours
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
##                     Estimate         SE        2.5%      97.5%
## X.Intercept.     0.045361845 0.05607430 -0.05769494 0.15644585
## gndr.c           0.177312043 0.01194123  0.15398419 0.19937802
## gdi.z.cm        -0.069407600 0.05850602 -0.18643084 0.03934691
## gndr.c.gdi.z.cm -0.008167377 0.01225828 -0.03131668 0.01647772
## w11             -0.065323912 0.06077571 -0.19207467 0.05096969
## w21             -0.073491289 0.05680982 -0.18799816 0.03097145
## b11             -0.189895139 0.17667361 -0.55835674 0.14816774
## b21             -0.213637521 0.16514485 -0.54650641 0.09003332
## r_xy1           -0.194537668 0.18099290 -0.57200736 0.15179012
## r_xy2           -0.208658018 0.16129563 -0.53376833 0.08793480
## q_b              0.024679719 0.03872278 -0.05571630 0.09535758
## q                0.013419217 0.04306516 -0.07664176 0.09024706
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
## [1] 0.02467972
## 
## $se
## [1] 0.03872278
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
## [1] 3.219803
## 
## $p_low
## [1] 0.0006413944
## 
## $z_high
## [1] -1.945115
## 
## $p_high
## [1] 0.02588055
## 
## $ci_level
## [1] 0.9
## 
## $ci_lower
## [1] -0.03901359
## 
## $ci_upper
## [1] 0.08837303
## 
## $equivalent
## [1] TRUE
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
## [1] 0.01341922
## 
## $se
## [1] 0.04306516
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
## [1] 2.633665
## 
## $p_low
## [1] 0.004223435
## 
## $z_high
## [1] -2.01046
## 
## $p_high
## [1] 0.02219126
## 
## $ci_level
## [1] 0.9
## 
## $ci_lower
## [1] -0.05741667
## 
## $ci_upper
## [1] 0.0842551
## 
## $equivalent
## [1] TRUE
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
  xlab("Gender Development Index (Average 2002-2023)")+
  #ylim(c(min.y.pred,max.y.pred))+
  #xlim(c(0.60,1.00))+
  ylab("Value achievement (Average 2002-2023)")+
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
  xlab("Gender Development Index (Average 2002-2023)")+
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
## Warning in geom_text(inherit.aes = F, aes(x = 0.9, y = -0.45, label = coef_q, : All aesthetics have length 1, but the data has 262 rows.
## ℹ Please consider using `annotate()` or provide this layer with data containing a single row.
```

``` r
pflag_comb
```

![](Analysis_achievement_files/figure-html/unnamed-chunk-35-1.png)<!-- -->

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
## 1462312.9 1462401.8 -731148.5 1462296.9    492335 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -5.9121 -0.6698  0.0262  0.6595  5.0732 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr  
##  cntry    (Intercept) 0.085766 0.29286        
##           gndr.c      0.003833 0.06191  -0.23 
##  Residual             1.002867 1.00143        
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##                     Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)          0.04309    0.05027 34.00242   0.857  0.39734    
## gndr.c               0.17696    0.01111 33.86732  15.933  < 2e-16 ***
## log_gdp.z.cm        -0.15659    0.05042 34.02542  -3.106  0.00381 ** 
## gndr.c:log_gdp.z.cm  0.01734    0.01128 35.39170   1.537  0.13307    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c lg_g..
## gndr.c      -0.220              
## lg_gdp.z.cm  0.023 -0.005       
## gndr.c:l_.. -0.005 -0.007 -0.217
```

``` r
getFE(mod2_log_GDP,round=3)
```

```
##                       Est.    SE     df      t     p     LL     UL
## (Intercept)          0.043 0.050 34.002  0.857 0.397 -0.059  0.145
## gndr.c               0.177 0.011 33.867 15.933 0.000  0.154  0.200
## log_gdp.z.cm        -0.157 0.050 34.025 -3.106 0.004 -0.259 -0.054
## gndr.c:log_gdp.z.cm  0.017 0.011 35.392  1.537 0.133 -0.006  0.040
```

``` r
getVC(mod2_log_GDP)
```

```
##        grp        var1   var2 sdcor vcov
## 1    cntry (Intercept)   <NA>  0.29 0.09
## 2    cntry      gndr.c   <NA>  0.06 0.00
## 3    cntry (Intercept) gndr.c -0.23 0.00
## 4 Residual        <NA>   <NA>  1.00 1.00
```

``` r
r2mlm(mod2_log_GDP,bargraph = F)
```

```
## $Decompositions
##                        total
## fixed           0.0244652925
## slope variation 0.0008524127
## mean variation  0.0770608795
## sigma2          0.8976214153
## 
## $R2s
##            total
## f   0.0244652925
## v   0.0008524127
## m   0.0770608795
## fv  0.0253177052
## fvm 0.1023785847
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
## Time difference of 29.84908 secs
```

``` r
round(ddsc_mod2_log_GDP$vpc_at_moderator_values,3)
```

```
##   lvl1.value Intercept.var Intercept.sd Residual.var Total.var   VPC mean.n.obs  ICC2 ICC2_SP
## 1       -0.5         0.118        0.343        1.003     1.121 0.105   7802.647 0.999   0.999
## 2        0.5         0.104        0.323        1.003     1.107 0.094   6678.029 0.998   0.999
```

``` r
round(ddsc_mod2_log_GDP$descriptives,3)
```

```
##                          M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled log_gdp.z.cm
## means_y1             0.098 0.336    1.000           1.000    0.983           0.983       -0.442
## means_y1_scaled      0.284 0.976    1.000           1.000    0.983           0.983       -0.442
## means_y2            -0.073 0.352    0.983           0.983    1.000           1.000       -0.453
## means_y2_scaled     -0.212 1.024    0.983           0.983    1.000           1.000       -0.453
## log_gdp.z.cm        -0.024 1.012   -0.442          -0.442   -0.453          -0.453        1.000
## log_gdp.z.cm_scaled  0.000 1.000   -0.442          -0.442   -0.453          -0.453        1.000
## diff_score           0.171 0.065   -0.163          -0.163   -0.340          -0.340        0.169
## diff_score_scaled    0.496 0.188   -0.163          -0.163   -0.340          -0.340        0.169
##                     log_gdp.z.cm_scaled diff_score diff_score_scaled
## means_y1                         -0.442     -0.163            -0.163
## means_y1_scaled                  -0.442     -0.163            -0.163
## means_y2                         -0.453     -0.340            -0.340
## means_y2_scaled                  -0.453     -0.340            -0.340
## log_gdp.z.cm                      1.000      0.169             0.169
## log_gdp.z.cm_scaled               1.000      0.169             0.169
## diff_score                        0.169      1.000             1.000
## diff_score_scaled                 0.169      1.000             1.000
```

``` r
round(ddsc_mod2_log_GDP$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                      -0.268 0.174 35.392  -1.537   0.133   -0.621    0.086
## w_11                         -0.165 0.052 34.051  -3.182   0.003   -0.271   -0.060
## w_21                         -0.148 0.049 34.042  -2.988   0.005   -0.249   -0.047
## r_xy1                        -0.492 0.155 34.051  -3.182   0.003   -0.806   -0.178
## r_xy2                        -0.420 0.141 34.042  -2.988   0.005   -0.706   -0.134
## b_11                         -0.480 0.151 34.051  -3.182   0.003   -0.787   -0.174
## b_21                         -0.430 0.144 34.042  -2.988   0.005   -0.722   -0.138
## main_effect                  -0.157 0.050 34.025  -3.106   0.004   -0.259   -0.054
## moderator_effect              0.177 0.011 33.867  15.933   0.000    0.154    0.200
## interaction                   0.017 0.011 35.392   1.537   0.133   -0.006    0.040
## q_b11_b21                    -0.064    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                  -0.091    NA     NA      NA      NA       NA       NA
## cross_over_point            -10.208    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.139 0.049 34.079  -2.830   0.008   -0.239   -0.039
## interaction_vs_main_bscale   -0.405 0.143 34.079  -2.830   0.008   -0.696   -0.114
## interaction_vs_main_rscale   -0.384 0.136 34.081  -2.816   0.008   -0.661   -0.107
## dadas                        -0.296 0.099 34.042  -2.988   0.997   -0.497   -0.095
## dadas_bscale                 -0.860 0.288 34.042  -2.988   0.997   -1.445   -0.275
## dadas_rscale                 -0.840 0.281 34.042  -2.988   0.997   -1.411   -0.269
## abs_diff                      0.017 0.011 35.392   1.537   0.067   -0.006    0.040
## abs_sum                       0.313 0.101 34.025   3.106   0.002    0.108    0.518
## abs_diff_bscale               0.050 0.033 35.392   1.537   0.067   -0.016    0.117
## abs_sum_bscale                0.910 0.293 34.025   3.106   0.002    0.315    1.506
## abs_diff_rscale               0.072 0.035 35.773   2.062   0.023    0.001    0.143
## abs_sum_rscale                0.912 0.293 34.025   3.108   0.002    0.316    1.509
```

``` r
round(ddsc_mod2_log_GDP$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.007 -0.321  3.304  1.000  0.069
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
##                                     est     se       z pvalue ci.lower ci.upper
## r_xy1_y2                         -0.169  0.169  -1.000  0.318   -0.500    0.162
## r_xy1                            -0.453  0.153  -2.961  0.003   -0.752   -0.153
## r_xy2                            -0.442  0.154  -2.876  0.004   -0.744   -0.141
## b_11                             -0.463  0.157  -2.961  0.003   -0.770   -0.157
## b_21                             -0.432  0.150  -2.876  0.004   -0.726   -0.137
## b_10                             -0.212  0.154  -1.378  0.168   -0.515    0.090
## b_20                              0.284  0.148   1.919  0.055   -0.006    0.574
## res_cov_y1_y2                     0.759  0.186   4.080  0.000    0.394    1.124
## diff_b10_b20                     -0.496  0.031 -15.828  0.000   -0.558   -0.435
## diff_b11_b21                     -0.032  0.032  -1.000  0.318   -0.094    0.031
## diff_rxy1_rxy2                   -0.010  0.031  -0.335  0.737   -0.072    0.051
## q_b11_b21                        -0.040  0.042  -0.953  0.340   -0.122    0.042
## q_rxy1_rxy2                      -0.013  0.039  -0.335  0.737   -0.090    0.063
## cross_over_point                -15.600 15.637  -0.998  0.318  -46.247   15.048
## sum_b11_b21                      -0.895  0.305  -2.934  0.003   -1.493   -0.297
## main_effect                      -0.448  0.153  -2.934  0.003   -0.746   -0.149
## interaction_vs_main_effect       -0.416  0.149  -2.784  0.005   -0.708   -0.123
## diff_abs_b11_abs_b21              0.032  0.032   1.000  0.318   -0.031    0.094
## abs_diff_b11_b21                  0.032  0.032   1.000  0.159   -0.031    0.094
## abs_sum_b11_b21                   0.895  0.305   2.934  0.002    0.297    1.493
## dadas                            -0.863  0.300  -2.876  0.998   -1.452   -0.275
## q_r_equivalence                  -0.087  0.039  -2.228  0.013       NA       NA
## q_b_equivalence                  -0.060  0.042  -1.443  0.075       NA       NA
## cross_over_point_equivalence     15.600 15.637   0.998  0.841       NA       NA
## cross_over_point_minimal_effect  15.600 15.637   0.998  0.159       NA       NA
```

``` r
round(ddsc_sem_log_GDP$variance_test,3)
```

```
##             est    se       z pvalue ci.lower ci.upper
## cov_y1y2  0.953 0.233   4.088   0.00    0.496    1.410
## var_y1    1.017 0.247   4.123   0.00    0.533    1.500
## var_y2    0.924 0.224   4.123   0.00    0.485    1.364
## var_diff  0.093 0.064   1.438   0.15   -0.034    0.219
## var_ratio 1.100 0.068  16.067   0.00    0.966    1.234
## cor_y1y2  0.983 0.006 174.148   0.00    0.972    0.994
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
## r_xy1y2                      -0.169 0.174 32.000  -0.970   0.339   -0.524    0.186
## w_11                         -0.159 0.054 32.696  -2.932   0.006   -0.270   -0.049
## w_21                         -0.149 0.054 32.696  -2.731   0.010   -0.259   -0.038
## r_xy1                        -0.453 0.154 32.696  -2.932   0.006   -0.767   -0.138
## r_xy2                        -0.442 0.162 32.696  -2.731   0.010   -0.772   -0.113
## b_11                         -0.464 0.158 32.696  -2.932   0.006   -0.785   -0.142
## b_21                         -0.432 0.158 32.696  -2.731   0.010   -0.754   -0.110
## main_effect                  -0.154 0.054 32.000  -2.847   0.008   -0.264   -0.044
## moderator_effect              0.171 0.011 32.000  15.355   0.000    0.148    0.193
## interaction                   0.011 0.011 32.000   0.970   0.339   -0.012    0.034
## q_b11_b21                    -0.040    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                  -0.013    NA     NA      NA      NA       NA       NA
## cross_over_point            -15.600    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.143 0.055 34.781  -2.589   0.014   -0.255   -0.031
## interaction_vs_main_bscale   -0.416 0.161 34.781  -2.589   0.014   -0.742   -0.090
## interaction_vs_main_rscale   -0.437 0.168 34.595  -2.599   0.014   -0.779   -0.096
## dadas                        -0.297 0.109 32.696  -2.731   0.995   -0.518   -0.076
## dadas_bscale                 -0.863 0.316 32.696  -2.731   0.995   -1.507   -0.220
## dadas_rscale                 -0.885 0.324 32.696  -2.731   0.995   -1.544   -0.225
## abs_diff                      0.011 0.011 32.000   0.970   0.170   -0.012    0.034
## abs_sum                       0.308 0.108 32.000   2.847   0.004    0.088    0.528
## abs_diff_bscale               0.032 0.033 32.000   0.970   0.170   -0.035    0.099
## abs_sum_bscale                0.895 0.315 32.000   2.847   0.004    0.255    1.536
## abs_diff_rscale               0.010 0.034 35.340   0.311   0.379   -0.058    0.079
## abs_sum_rscale                0.895 0.315 32.000   2.844   0.004    0.254    1.536
```

``` r
# country-time multilevel model


mod2_log_GDP_cntry_year<-
  lmer(ach.z.wt~gndr.c+log_gdp.z.cm+gndr.c:log_gdp.z.cm+
                 (gndr.c|cntry),data=diff_dat_cntry_year,REML=F,
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2_log_GDP_cntry_year)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ach.z.wt ~ gndr.c + log_gdp.z.cm + gndr.c:log_gdp.z.cm + (gndr.c |      cntry)
##    Data: diff_dat_cntry_year
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##    -621.7    -587.3     318.9    -637.7       538 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -3.6949 -0.6187  0.0037  0.6596  3.5793 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr  
##  cntry    (Intercept) 0.084229 0.29022        
##           gndr.c      0.001049 0.03238  -0.59 
##  Residual             0.013586 0.11656        
## Number of obs: 546, groups:  cntry, 34
## 
## Fixed effects:
##                     Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)          0.03792    0.05010 33.99255   0.757  0.45427    
## gndr.c               0.17778    0.01173 36.42341  15.162  < 2e-16 ***
## log_gdp.z.cm        -0.15408    0.05034 34.25942  -3.061  0.00427 ** 
## gndr.c:log_gdp.z.cm  0.01628    0.01246 40.10777   1.307  0.19882    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c lg_g..
## gndr.c      -0.276              
## lg_gdp.z.cm  0.019 -0.006       
## gndr.c:l_.. -0.006 -0.177 -0.260
```

``` r
getFE(mod2_log_GDP_cntry_year,round=3)
```

```
##                       Est.    SE     df      t     p     LL     UL
## (Intercept)          0.038 0.050 33.993  0.757 0.454 -0.064  0.140
## gndr.c               0.178 0.012 36.423 15.162 0.000  0.154  0.202
## log_gdp.z.cm        -0.154 0.050 34.259 -3.061 0.004 -0.256 -0.052
## gndr.c:log_gdp.z.cm  0.016 0.012 40.108  1.307 0.199 -0.009  0.041
```

``` r
getVC(mod2_log_GDP_cntry_year)
```

```
##        grp        var1   var2 sdcor  vcov
## 1    cntry (Intercept)   <NA>  0.29  0.08
## 2    cntry      gndr.c   <NA>  0.03  0.00
## 3    cntry (Intercept) gndr.c -0.59 -0.01
## 4 Residual        <NA>   <NA>  0.12  0.01
```

``` r
r2mlm(mod2_log_GDP,bargraph = F)
```

```
## $Decompositions
##                        total
## fixed           0.0244652925
## slope variation 0.0008524127
## mean variation  0.0770608795
## sigma2          0.8976214153
## 
## $R2s
##            total
## f   0.0244652925
## v   0.0008524127
## m   0.0770608795
## fv  0.0253177052
## fvm 0.1023785847
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
## 1       -0.5         0.116        0.340        0.014     0.129 0.895      8.029 0.998   0.986
## 2        0.5         0.100        0.316        0.014     0.114 0.880      8.029 0.998   0.983
```

``` r
round(ddsc_mod2_log_GDP_cntry_year$descriptives,3)
```

```
##                          M    SD means_y1 means_y1_scaled means_y2 means_y2_scaled log_gdp.z.cm
## means_y1             0.129 0.326    1.000           1.000    0.982           0.982       -0.453
## means_y1_scaled      0.384 0.969    1.000           1.000    0.982           0.982       -0.453
## means_y2            -0.045 0.347    0.982           0.982    1.000           1.000       -0.474
## means_y2_scaled     -0.135 1.030    0.982           0.982    1.000           1.000       -0.474
## log_gdp.z.cm        -0.024 1.012   -0.453          -0.453   -0.474          -0.474        1.000
## log_gdp.z.cm_scaled  0.000 1.000   -0.453          -0.453   -0.474          -0.474        1.000
## diff_score           0.175 0.067   -0.206          -0.206   -0.389          -0.389        0.245
## diff_score_scaled    0.518 0.200   -0.206          -0.206   -0.389          -0.389        0.245
##                     log_gdp.z.cm_scaled diff_score diff_score_scaled
## means_y1                         -0.453     -0.206            -0.206
## means_y1_scaled                  -0.453     -0.206            -0.206
## means_y2                         -0.474     -0.389            -0.389
## means_y2_scaled                  -0.474     -0.389            -0.389
## log_gdp.z.cm                      1.000      0.245             0.245
## log_gdp.z.cm_scaled               1.000      0.245             0.245
## diff_score                        0.245      1.000             1.000
## diff_score_scaled                 0.245      1.000             1.000
```

``` r
round(ddsc_mod2_log_GDP_cntry_year$results,3)
```

```
##                            estimate    SE     df t.ratio p.value ci.lower ci.upper
## r_xy1y2                      -0.241 0.185 40.108  -1.307   0.199   -0.614    0.132
## w_11                         -0.162 0.052 34.379  -3.101   0.004   -0.268   -0.056
## w_21                         -0.146 0.049 34.364  -2.973   0.005   -0.246   -0.046
## r_xy1                        -0.497 0.160 34.379  -3.101   0.004   -0.822   -0.171
## r_xy2                        -0.421 0.142 34.364  -2.973   0.005   -0.708   -0.133
## b_11                         -0.482 0.155 34.379  -3.101   0.004   -0.798   -0.166
## b_21                         -0.434 0.146 34.364  -2.973   0.005   -0.730   -0.137
## main_effect                  -0.154 0.050 34.259  -3.061   0.004   -0.256   -0.052
## moderator_effect              0.178 0.012 36.423  15.162   0.000    0.154    0.202
## interaction                   0.016 0.012 40.108   1.307   0.199   -0.009    0.041
## q_b11_b21                    -0.061    NA     NA      NA      NA       NA       NA
## q_rxy1_rxy2                  -0.096    NA     NA      NA      NA       NA       NA
## cross_over_point            -10.920    NA     NA      NA      NA       NA       NA
## interaction_vs_main          -0.138 0.049 34.654  -2.835   0.008   -0.237   -0.039
## interaction_vs_main_bscale   -0.409 0.144 34.654  -2.835   0.008   -0.703   -0.116
## interaction_vs_main_rscale   -0.383 0.136 34.689  -2.819   0.008   -0.659   -0.107
## dadas                        -0.292 0.098 34.364  -2.973   0.997   -0.491   -0.092
## dadas_bscale                 -0.867 0.292 34.364  -2.973   0.997   -1.460   -0.275
## dadas_rscale                 -0.842 0.283 34.364  -2.973   0.997   -1.417   -0.267
## abs_diff                      0.016 0.012 40.108   1.307   0.099   -0.009    0.041
## abs_sum                       0.308 0.101 34.259   3.061   0.002    0.104    0.513
## abs_diff_bscale               0.048 0.037 40.108   1.307   0.099   -0.026    0.123
## abs_sum_bscale                0.915 0.299 34.259   3.061   0.002    0.308    1.523
## abs_diff_rscale               0.076 0.040 40.598   1.883   0.033   -0.006    0.157
## abs_sum_rscale                0.918 0.300 34.260   3.063   0.002    0.309    1.527
```

``` r
round(ddsc_mod2_log_GDP_cntry_year$re_cov_test,3)
```

```
## RE_cov RE_cor  Chisq     Df      p 
## -0.008 -0.674  4.084  1.000  0.043
```

### Comparisons across models


``` r
# full multilevel model (individuals as entries)
round(ddsc_mod2_log_GDP$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df t.ratio p.value ci.lower ci.upper
## r_xy1             -0.4922 0.1547 34.0510 -3.1821  0.0031  -0.8064  -0.1779
## r_xy2             -0.4200 0.1405 34.0422 -2.9885  0.0052  -0.7056  -0.1344
## b_11              -0.4804 0.1510 34.0510 -3.1821  0.0031  -0.7872  -0.1736
## b_21              -0.4300 0.1439 34.0422 -2.9885  0.0052  -0.7224  -0.1376
## main_effect       -0.1566 0.0504 34.0254 -3.1060  0.0038  -0.2590  -0.0541
## moderator_effect   0.1770 0.0111 33.8673 15.9334  0.0000   0.1544   0.1995
## interaction        0.0173 0.0113 35.3917  1.5375  0.1331  -0.0055   0.0402
## q_b11_b21         -0.0636     NA      NA      NA      NA       NA       NA
```

``` r
# country-level structural path model (country-gender means as entries)
round(ddsc_sem_log_GDP$results[c("r_xy1","r_xy2","b_11","b_21","q_b11_b21","diff_b11_b21"),],4)
```

```
##                  est     se       z pvalue ci.lower ci.upper
## r_xy1        -0.4528 0.1529 -2.9609 0.0031  -0.7525  -0.1531
## r_xy2        -0.4423 0.1538 -2.8756 0.0040  -0.7438  -0.1408
## b_11         -0.4634 0.1565 -2.9609 0.0031  -0.7702  -0.1567
## b_21         -0.4316 0.1501 -2.8756 0.0040  -0.7258  -0.1374
## q_b11_b21    -0.0398 0.0417 -0.9534 0.3404  -0.1216   0.0420
## diff_b11_b21 -0.0318 0.0318 -0.9996 0.3175  -0.0942   0.0306
```

``` r
# multilevel model (country-gender means as entries)
round(ddsc_mod2_log_GDP_ri$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df t.ratio p.value ci.lower ci.upper
## r_xy1             -0.4528 0.1544 32.6965 -2.9320  0.0061  -0.7670  -0.1385
## r_xy2             -0.4423 0.1620 32.6965 -2.7307  0.0101  -0.7720  -0.1127
## b_11              -0.4636 0.1581 32.6965 -2.9320  0.0061  -0.7853  -0.1418
## b_21              -0.4317 0.1581 32.6965 -2.7307  0.0101  -0.7535  -0.1100
## main_effect       -0.1540 0.0541 32.0000 -2.8467  0.0076  -0.2642  -0.0438
## moderator_effect   0.1707 0.0111 32.0000 15.3553  0.0000   0.1481   0.1934
## interaction        0.0109 0.0113 32.0000  0.9698  0.3394  -0.0120   0.0339
## q_b11_b21         -0.0398     NA      NA      NA      NA       NA       NA
```

``` r
# multilevel model (country-year-gender means as entries)
round(ddsc_mod2_log_GDP_cntry_year$results[c("r_xy1","r_xy2","b_11","b_21","main_effect","moderator_effect","interaction","q_b11_b21"),],4)
```

```
##                  estimate     SE      df t.ratio p.value ci.lower ci.upper
## r_xy1             -0.4969 0.1602 34.3790 -3.1013  0.0038  -0.8223  -0.1714
## r_xy2             -0.4209 0.1416 34.3643 -2.9733  0.0054  -0.7085  -0.1333
## b_11              -0.4819 0.1554 34.3790 -3.1013  0.0038  -0.7976  -0.1663
## b_21              -0.4335 0.1458 34.3643 -2.9733  0.0054  -0.7298  -0.1373
## main_effect       -0.1541 0.0503 34.2594 -3.0610  0.0043  -0.2564  -0.0518
## moderator_effect   0.1778 0.0117 36.4234 15.1617  0.0000   0.1540   0.2015
## interaction        0.0163 0.0125 40.1078  1.3065  0.1988  -0.0089   0.0415
## q_b11_b21         -0.0612     NA      NA      NA      NA       NA       NA
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
## Time difference of 1.765186 hours
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
##                        Estimate         SE         2.5%       97.5%
## X.Intercept.         0.04183783 0.05071365 -0.051734432  0.14175388
## gndr.c               0.17730887 0.01166166  0.154201767  0.19843873
## log_gdp.z.cm        -0.15699727 0.04918863 -0.251547478 -0.06166777
## gndr.c.log_gdp.z.cm  0.01791668 0.01112141 -0.003175667  0.03960411
## w11                 -0.16595561 0.05081812 -0.267887803 -0.06515912
## w21                 -0.14803893 0.04814981 -0.241730731 -0.05641935
## b11                 -0.48242923 0.14772713 -0.778743805 -0.18941609
## b21                 -0.43034584 0.13997041 -0.702705787 -0.16400978
## r_xy1               -0.49422359 0.15133874 -0.797782424 -0.19404690
## r_xy2               -0.42031526 0.13670795 -0.686326988 -0.16018701
## q_b                 -0.07259528 0.05313666 -0.205772112  0.01110338
## q                   -0.10452905 0.06669952 -0.265198932 -0.01225944
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
## [1] -0.07259528
## 
## $se
## [1] 0.05313666
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
## [1] 0.5157403
## 
## $p_low
## [1] 0.3030179
## 
## $z_high
## [1] -3.248139
## 
## $p_high
## [1] 0.0005808118
## 
## $ci_level
## [1] 0.9
## 
## $ci_lower
## [1] -0.1599973
## 
## $ci_upper
## [1] 0.01480675
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
## [1] -0.1045291
## 
## $se
## [1] 0.06669952
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
## [1] -0.06790234
## 
## $p_low
## [1] 0.5270683
## 
## $z_high
## [1] -3.066425
## 
## $p_high
## [1] 0.001083177
## 
## $ci_level
## [1] 0.9
## 
## $ci_lower
## [1] -0.21424
## 
## $ci_upper
## [1] 0.005181897
## 
## $equivalent
## [1] FALSE
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
  xlab("log(GDP per capita) (Average 2002-2023)")+
  #ylim(c(min.y.pred,max.y.pred))+
  #xlim(c(0.60,1.00))+
  ylab("Value achievement (Average 2002-2023)")+
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
  xlab("log(GDP per capita) (Average 2002-2023)")+
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

![](Analysis_achievement_files/figure-html/unnamed-chunk-41-1.png)<!-- -->

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


## mod3: fixed effect of time (year)

* Year coded year-2002 (2002 is a zero-point)


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
## 1462127.2 1462205.0 -731056.6 1462113.2    492336 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -5.9551 -0.6697  0.0267  0.6594  5.0899 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr  
##  cntry    (Intercept) 0.11051  0.33244        
##           gndr.c      0.00413  0.06427  -0.33 
##  Residual             1.00247  1.00124        
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##               Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)  4.860e-02  5.704e-02  3.400e+01   0.852      0.4    
## gndr.c       1.768e-01  1.149e-02  3.397e+01  15.386   <2e-16 ***
## essround.c  -6.737e-03  4.851e-04  4.923e+05 -13.887   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##            (Intr) gndr.c
## gndr.c     -0.317       
## essround.c -0.002  0.000
```

``` r
getFE(mod3,round=3)
```

```
##               Est.    SE         df       t   p     LL     UL
## (Intercept)  0.049 0.057     33.997   0.852 0.4 -0.067  0.165
## gndr.c       0.177 0.011     33.968  15.386 0.0  0.153  0.200
## essround.c  -0.007 0.000 492265.387 -13.887 0.0 -0.008 -0.006
```

``` r
getVC(mod3)
```

```
##        grp        var1   var2 sdcor  vcov
## 1    cntry (Intercept)   <NA>  0.33  0.11
## 2    cntry      gndr.c   <NA>  0.06  0.00
## 3    cntry (Intercept) gndr.c -0.33 -0.01
## 4 Residual        <NA>   <NA>  1.00  1.00
```

``` r
r2mlm(mod3,bargraph = F)
```

```
## $Decompositions
##                        total
## fixed           0.0072930545
## slope variation 0.0009140936
## mean variation  0.0989252626
## sigma2          0.8928675893
## 
## $R2s
##            total
## f   0.0072930545
## v   0.0009140936
## m   0.0989252626
## fv  0.0082071481
## fvm 0.1071324107
```

``` r
anova(mod2,mod3)
```

```
## Data: diff_dat
## Models:
## mod2: ach.z ~ gndr.c + (gndr.c | cntry)
## mod3: ach.z ~ gndr.c + essround.c + (gndr.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L) Chisq Df Pr(>Chisq)    
## mod2    6 1462318 1462385 -731153   1462306                        
## mod3    7 1462127 1462205 -731057   1462113 192.8  1  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## mod4: random effect of time (year)


``` r
mod4<-lmer(ach.z~gndr.c+year.c+(gndr.c+year.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod4)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ach.z ~ gndr.c + year.c + (gndr.c + year.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##   1460756   1460867   -730368   1460736    492333 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.0494 -0.6687  0.0282  0.6587  5.1395 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr        
##  cntry    (Intercept) 0.205322 0.45312              
##           gndr.c      0.004232 0.06505  -0.47       
##           year.c      0.000310 0.01761  -0.63  0.35 
##  Residual             0.999344 0.99967              
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##              Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)  0.144290   0.077959 31.836171   1.851   0.0735 .  
## gndr.c       0.176211   0.011619 33.784152  15.166   <2e-16 ***
## year.c      -0.006487   0.003051 27.898663  -2.126   0.0424 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##        (Intr) gndr.c
## gndr.c -0.451       
## year.c -0.631  0.338
```

``` r
getFE(mod4,round=3)
```

```
##               Est.    SE     df      t     p     LL    UL
## (Intercept)  0.144 0.078 31.836  1.851 0.073 -0.015 0.303
## gndr.c       0.176 0.012 33.784 15.166 0.000  0.153 0.200
## year.c      -0.006 0.003 27.899 -2.126 0.042 -0.013 0.000
```

``` r
getVC(mod4)
```

```
##        grp        var1   var2 sdcor  vcov
## 1    cntry (Intercept)   <NA>  0.45  0.21
## 2    cntry      gndr.c   <NA>  0.07  0.00
## 3    cntry      year.c   <NA>  0.02  0.00
## 4    cntry (Intercept) gndr.c -0.47 -0.01
## 5    cntry (Intercept) year.c -0.63 -0.01
## 6    cntry      gndr.c year.c  0.35  0.00
## 7 Residual        <NA>   <NA>  1.00  1.00
```

``` r
r2mlm(mod4,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.008099436
## slope variation 0.011691419
## mean variation  0.117431783
## sigma2          0.862777362
## 
## $R2s
##           total
## f   0.008099436
## v   0.011691419
## m   0.117431783
## fv  0.019790855
## fvm 0.137222638
```

``` r
anova(mod2,mod3,mod4)
```

```
## Data: diff_dat
## Models:
## mod2: ach.z ~ gndr.c + (gndr.c | cntry)
## mod3: ach.z ~ gndr.c + essround.c + (gndr.c | cntry)
## mod4: ach.z ~ gndr.c + year.c + (gndr.c + year.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)    
## mod2    6 1462318 1462385 -731153   1462306                         
## mod3    7 1462127 1462205 -731057   1462113  192.8  1  < 2.2e-16 ***
## mod4   10 1460756 1460867 -730368   1460736 1377.2  3  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## mod5: fixed interaction between time and gender


``` r
mod5<-lmer(ach.z~gndr.c+year.c+
             gndr.c:year.c+(gndr.c+year.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod5)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ach.z ~ gndr.c + year.c + gndr.c:year.c + (gndr.c + year.c |      cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1460670.4 1460792.5 -730324.2 1460648.4    492332 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.0787 -0.6693  0.0273  0.6590  5.1246 
## 
## Random effects:
##  Groups   Name        Variance  Std.Dev. Corr        
##  cntry    (Intercept) 0.2052466 0.45304              
##           gndr.c      0.0037951 0.06160  -0.41       
##           year.c      0.0003096 0.01759  -0.63  0.26 
##  Residual             0.9991687 0.99958              
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##                 Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)    1.450e-01  7.795e-02  3.182e+01   1.860   0.0721 .  
## gndr.c         2.230e-01  1.211e-02  4.894e+01  18.410   <2e-16 ***
## year.c        -6.565e-03  3.049e-03  2.788e+01  -2.153   0.0401 *  
## gndr.c:year.c -4.398e-03  4.686e-04  1.229e+05  -9.386   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c year.c
## gndr.c      -0.359              
## year.c      -0.631  0.226       
## gndr.c:yr.c -0.001 -0.409  0.002
```

``` r
getFE(mod5,round=3)
```

```
##                 Est.    SE         df      t     p     LL     UL
## (Intercept)    0.145 0.078     31.822  1.860 0.072 -0.014  0.304
## gndr.c         0.223 0.012     48.935 18.410 0.000  0.199  0.247
## year.c        -0.007 0.003     27.875 -2.153 0.040 -0.013  0.000
## gndr.c:year.c -0.004 0.000 122917.031 -9.386 0.000 -0.005 -0.003
```

``` r
getVC(mod5)
```

```
##        grp        var1   var2 sdcor  vcov
## 1    cntry (Intercept)   <NA>  0.45  0.21
## 2    cntry      gndr.c   <NA>  0.06  0.00
## 3    cntry      year.c   <NA>  0.02  0.00
## 4    cntry (Intercept) gndr.c -0.41 -0.01
## 5    cntry (Intercept) year.c -0.63 -0.01
## 6    cntry      gndr.c year.c  0.26  0.00
## 7 Residual        <NA>   <NA>  1.00  1.00
```

``` r
r2mlm(mod5,bargraph = F)
```

```
## $Decompositions
##                      total
## fixed           0.00837845
## slope variation 0.01158101
## mean variation  0.11731773
## sigma2          0.86272280
## 
## $R2s
##          total
## f   0.00837845
## v   0.01158101
## m   0.11731773
## fv  0.01995946
## fvm 0.13727720
```

``` r
anova(mod4,mod5)
```

```
## Data: diff_dat
## Models:
## mod4: ach.z ~ gndr.c + year.c + (gndr.c + year.c | cntry)
## mod5: ach.z ~ gndr.c + year.c + gndr.c:year.c + (gndr.c + year.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L) Chisq Df Pr(>Chisq)    
## mod4   10 1460756 1460867 -730368   1460736                        
## mod5   11 1460670 1460793 -730324   1460648 87.65  1  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## mod6: random interaction between time and gender


``` r
mod6<-lmer(ach.z~gndr.c+year.c+
             gndr.c:year.c+(gndr.c+year.c+gndr.c:year.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ach.z ~ gndr.c + year.c + gndr.c:year.c + (gndr.c + year.c +      gndr.c:year.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1460673.5 1460840.1 -730321.8 1460643.5    492328 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.0694 -0.6693  0.0266  0.6587  5.1223 
## 
## Random effects:
##  Groups   Name          Variance  Std.Dev. Corr              
##  cntry    (Intercept)   2.047e-01 0.452477                   
##           gndr.c        5.062e-03 0.071149 -0.58             
##           year.c        3.076e-04 0.017538 -0.63  0.46       
##           gndr.c:year.c 4.990e-06 0.002234  0.47 -0.48 -0.46 
##  Residual               9.991e-01 0.999567                   
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##                 Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)    0.1447148  0.0778501 31.7708013   1.859   0.0723 .  
## gndr.c         0.2185334  0.0137229 23.8163062  15.925 3.39e-14 ***
## year.c        -0.0065480  0.0030391 27.8077075  -2.155   0.0400 *  
## gndr.c:year.c -0.0041503  0.0006251 20.2744957  -6.639 1.70e-06 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c year.c
## gndr.c      -0.517              
## year.c      -0.629  0.406       
## gndr.c:yr.c  0.287 -0.569 -0.284
```

``` r
getFE(mod6,round=3)
```

```
##                 Est.    SE     df      t     p     LL     UL
## (Intercept)    0.145 0.078 31.771  1.859 0.072 -0.014  0.303
## gndr.c         0.219 0.014 23.816 15.925 0.000  0.190  0.247
## year.c        -0.007 0.003 27.808 -2.155 0.040 -0.013  0.000
## gndr.c:year.c -0.004 0.001 20.274 -6.639 0.000 -0.005 -0.003
```

``` r
getVC(mod6)
```

```
##         grp          var1          var2 sdcor  vcov
## 1     cntry   (Intercept)          <NA>  0.45  0.20
## 2     cntry        gndr.c          <NA>  0.07  0.01
## 3     cntry        year.c          <NA>  0.02  0.00
## 4     cntry gndr.c:year.c          <NA>  0.00  0.00
## 5     cntry   (Intercept)        gndr.c -0.58 -0.02
## 6     cntry   (Intercept)        year.c -0.63  0.00
## 7     cntry   (Intercept) gndr.c:year.c  0.47  0.00
## 8     cntry        gndr.c        year.c  0.46  0.00
## 9     cntry        gndr.c gndr.c:year.c -0.48  0.00
## 10    cntry        year.c gndr.c:year.c -0.46  0.00
## 11 Residual          <NA>          <NA>  1.00  1.00
```

``` r
r2mlm(mod6,bargraph = F)
```

```
## $Decompositions
##                       total
## fixed           0.008209563
## slope variation 0.011652765
## mean variation  0.117441162
## sigma2          0.862696510
## 
## $R2s
##           total
## f   0.008209563
## v   0.011652765
## m   0.117441162
## fv  0.019862328
## fvm 0.137303490
```

``` r
anova(mod4,mod5,mod6)
```

```
## Data: diff_dat
## Models:
## mod4: ach.z ~ gndr.c + year.c + (gndr.c + year.c | cntry)
## mod5: ach.z ~ gndr.c + year.c + gndr.c:year.c + (gndr.c + year.c | cntry)
## mod6: ach.z ~ gndr.c + year.c + gndr.c:year.c + (gndr.c + year.c + gndr.c:year.c | cntry)
##      npar     AIC     BIC  logLik -2*log(L)   Chisq Df Pr(>Chisq)    
## mod4   10 1460756 1460867 -730368   1460736                          
## mod5   11 1460670 1460793 -730324   1460648 87.6503  1     <2e-16 ***
## mod6   15 1460674 1460840 -730322   1460644  4.8151  4     0.3068    
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
##      21 -0.0585 0.0638 32.6  -0.1882   0.0713  -0.917  0.3657
##       0  0.0354 0.0816 31.6  -0.1309   0.2018   0.434  0.6670
## 
## gndr.c =  0.5:
##  year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21  0.0729 0.0615 32.5  -0.0523   0.1981   1.186  0.2444
##       0  0.2540 0.0745 31.3   0.1020   0.4059   3.408  0.0018
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
##  year.c21 - year.c0  -0.0939 0.0660 27.3   -0.229   0.0414  -1.424  0.1659
## 
## gndr.c =  0.5:
##  contrast           estimate     SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0  -0.1811 0.0623 26.8   -0.309  -0.0533  -2.908  0.0072
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
##    -0.5     21 -0.0585 0.0638 32.6  -0.1882   0.0713  -0.917  0.3657
##     0.5     21  0.0729 0.0615 32.5  -0.0523   0.1981   1.186  0.2444
##    -0.5      0  0.0354 0.0816 31.6  -0.1309   0.2018   0.434  0.6670
##     0.5      0  0.2540 0.0745 31.3   0.1020   0.4059   3.408  0.0018
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

``` r
pairs(change_in_diff_mod6,adjust="none",infer=c(T,T))
```

```
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21    -0.1314 0.0125 32.2   -0.157  -0.1060 -10.530 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)  -0.0939 0.0660 27.3   -0.229   0.0414  -1.424  0.1659
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     -0.3125 0.0611 27.2   -0.438  -0.1872  -5.118 <0.0001
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)      0.0374 0.0684 28.0   -0.103   0.1775   0.548  0.5883
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0        -0.1811 0.0623 26.8   -0.309  -0.0533  -2.908  0.0072
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0      -0.2185 0.0137 23.8   -0.247  -0.1902 -15.925 <0.0001
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
##  diff_ESS11    0.131 0.0125 32.2    0.106    0.157  10.530 <0.0001
##  diff_ESS1     0.219 0.0137 23.8    0.190    0.247  15.925 <0.0001
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
##  diff_ESS11 - diff_ESS1  -0.0872 0.0131 20.3   -0.115  -0.0598  -6.639 <0.0001
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
## [1] -0.08715702
## 
## $se
## [1] 0.01312713
## 
## $df
## [1] 20.2745
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
## [1] 8.596163
## 
## $p_low
## [1] 1.692104e-08
## 
## $t_high
## [1] -21.87507
## 
## $p_high
## [1] 7.073859e-16
## 
## $ci_level
## [1] 0.8
## 
## $ci_lower
## [1] -0.1045469
## 
## $ci_upper
## [1] -0.06976714
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

![](Analysis_achievement_files/figure-html/unnamed-chunk-48-1.png)<!-- -->

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
pred_cntry_dat$ach.z_mean<-predict(mod6,newdata=pred_cntry_dat)

pred_cntry_dat$year=pred_cntry_dat$year.c+2002

pred_cntry_dat$gender<-
  case_when(
    pred_cntry_dat$gndr.c==0.5~"men",
    pred_cntry_dat$gndr.c==-0.5~"women",
  )

range(pred_cntry_dat$ach.z_mean)
```

```
## [1] -0.7583059  0.8644311
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

pdf("../results/ach/country_specific_time_trends.pdf", width = 8, height = 6)

for (ctry in countries) {
  print(
    ggplot(pred_cntry_dat[pred_cntry_dat$cntry == ctry, ],
       aes(x = year, y = ach.z_mean, color = gender)) +
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
    name   = "Mean-level of value achievement",
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
         aes(x = year, y = ach.z_mean, color = gender)) +
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
    name   = "Mean-level of value achievement",
    sec.axis = sec_axis(~ scale_y_to_gei(.),
                        name = "Gender Equality Index (GEI)")
  ) +
  scale_x_continuous(limits = c(2001, 2024),
                     breaks = c(seq(2002, 2020, 2),2023)) +
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

```
## Warning: Removed 12 rows containing missing values or values outside the scale range (`geom_line()`).
```

![](Analysis_achievement_files/figure-html/unnamed-chunk-50-1.png)<!-- -->

``` r
png(filename = 
      "../results/ach/country_time_trend_facets.png",
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
## [1] 61.6387
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
## 1    0.32                0.09                   -0.10                    0.04                      0.15
## 2    0.18                0.12                   -0.10                    0.07                      0.17
## 3    0.18               -0.14                   -0.03                   -0.15                     -0.12
## 4    0.31                0.07                   -0.09                    0.02                      0.11
## 5    0.16               -0.24                   -0.10                   -0.29                     -0.19
## 6    0.35                0.34                   -0.13                    0.28                      0.40
## 7    0.30               -0.16                   -0.07                   -0.19                     -0.12
## 8    0.21                0.15                   -0.07                    0.11                      0.18
## 9    0.22               -0.19                   -0.08                   -0.23                     -0.15
## 10   0.25               -0.27                   -0.09                   -0.32                     -0.23
## 11   0.34               -0.15                   -0.07                   -0.19                     -0.12
## 12   0.28               -0.20                   -0.11                   -0.25                     -0.14
## 13   0.28               -0.12                   -0.12                   -0.18                     -0.06
## 14   0.19               -0.38                   -0.08                   -0.42                     -0.34
## 15   0.21               -0.28                   -0.09                   -0.32                     -0.24
## 16   0.21                0.03                   -0.11                   -0.03                      0.08
## 17   0.17               -0.11                   -0.05                   -0.13                     -0.08
## 18   0.07               -0.24                   -0.05                   -0.26                     -0.21
## 19   0.19                0.02                   -0.15                   -0.05                      0.10
## 20   0.21                0.03                   -0.08                   -0.01                      0.08
## 21   0.14               -0.21                   -0.06                   -0.24                     -0.19
## 22   0.20               -0.86                   -0.03                   -0.87                     -0.85
## 23   0.05               -1.79                    0.01                   -1.78                     -1.79
## 24   0.28               -0.01                   -0.10                   -0.06                      0.04
## 25   0.20               -0.03                   -0.12                   -0.09                      0.03
## 26   0.25               -0.11                   -0.10                   -0.16                     -0.06
## 27   0.24                0.00                   -0.15                   -0.08                      0.07
## 28   0.20               -0.30                   -0.07                   -0.33                     -0.26
## 29   0.19               -0.34                   -0.11                   -0.39                     -0.29
## 30   0.17               -0.13                   -0.09                   -0.18                     -0.08
## 31   0.17                0.24                   -0.08                    0.20                      0.28
## 32   0.27                0.03                   -0.10                   -0.02                      0.08
## 33   0.20                0.35                   -0.09                    0.31                      0.40
## 34   0.22                0.11                   -0.12                    0.05                      0.17
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
## 1     GR               -0.38
## 2     RU               -0.34
## 3     HR               -0.28
## 4     ES               -0.27
## 5     CY               -0.24
## 6     IL               -0.24
## 7     LT               -0.21
## 8     FR               -0.20
## 9     EE               -0.19
## 10    DE               -0.16
## 11    FI               -0.15
## 12    BG               -0.14
## 13    SE               -0.13
## 14    GB               -0.12
## 15    IE               -0.11
## 16    PL               -0.11
## 17    NO               -0.03
## 18    NL               -0.01
## 19    PT                0.00
## 20    IS                0.02
## 21    HU                0.03
## 22    IT                0.03
## 23    SK                0.03
## 24    CH                0.07
## 25    AT                0.09
## 26    UA                0.11
## 27    BE                0.12
## 28    DK                0.15
## 29    SI                0.24
## 30    CZ                0.34
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
## 1     IS                   -0.15
## 2     PT                   -0.15
## 3     CZ                   -0.13
## 4     GB                   -0.12
## 5     NO                   -0.12
## 6     UA                   -0.12
## 7     FR                   -0.11
## 8     HU                   -0.11
## 9     RU                   -0.11
## 10    AT                   -0.10
## 11    BE                   -0.10
## 12    CY                   -0.10
## 13    NL                   -0.10
## 14    PL                   -0.10
## 15    SK                   -0.10
## 16    CH                   -0.09
## 17    ES                   -0.09
## 18    HR                   -0.09
## 19    SE                   -0.09
## 20    EE                   -0.08
## 21    GR                   -0.08
## 22    IT                   -0.08
## 23    SI                   -0.08
## 24    DE                   -0.07
## 25    DK                   -0.07
## 26    FI                   -0.07
## 27    LT                   -0.06
## 28    IE                   -0.05
## 29    IL                   -0.05
## 30    BG                   -0.03
```


#### Country-level predictors: GEI


``` r
mod6_GEI<-lmer(ach.z~gndr.c+year.c+
             gndr.c:year.c+
               gei.z.cm:gndr.c+gei.z.cm:year.c+gei.z.cm:gndr.c:year.c+
               (gndr.c+year.c+gndr.c:year.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_GEI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ach.z ~ gndr.c + year.c + gndr.c:year.c + gei.z.cm:gndr.c + gei.z.cm:year.c +  
##     gei.z.cm:gndr.c:year.c + (gndr.c + year.c + gndr.c:year.c |      cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1421192.3 1421391.8 -710578.1 1421156.3    480346 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.0855 -0.6705  0.0275  0.6603  5.1314 
## 
## Random effects:
##  Groups   Name          Variance  Std.Dev. Corr              
##  cntry    (Intercept)   2.116e-01 0.459958                   
##           gndr.c        6.736e-03 0.082075 -0.72             
##           year.c        3.883e-04 0.019705 -0.74  0.60       
##           gndr.c:year.c 7.967e-06 0.002823  0.76 -0.65 -0.67 
##  Residual               9.950e-01 0.997516                   
## Number of obs: 480364, groups:  cntry, 33
## 
## Fixed effects:
##                          Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)             0.1551187  0.0803300 30.9020395   1.931   0.0627 .  
## gndr.c                  0.2182023  0.0157100 15.7715933  13.889 2.93e-10 ***
## year.c                 -0.0068685  0.0034600 23.0081424  -1.985   0.0592 .  
## gndr.c:year.c          -0.0041807  0.0007174 14.0890544  -5.828 4.28e-05 ***
## gndr.c:gei.z.cm        -0.0218479  0.0121267 35.1228732  -1.802   0.0802 .  
## year.c:gei.z.cm        -0.0062072  0.0024339 35.3076593  -2.550   0.0153 *  
## gndr.c:year.c:gei.z.cm  0.0018172  0.0007433 27.1329538   2.445   0.0213 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c year.c gnd.:. gn.:.. yr.:..
## gndr.c      -0.653                                   
## year.c      -0.741  0.540                            
## gndr.c:yr.c  0.521 -0.666 -0.460                     
## gndr.c:g.z. -0.001 -0.064  0.003  0.102              
## yr.c:g.z.cm  0.005  0.000 -0.019  0.002  0.102       
## gndr.c:.:..  0.000  0.093  0.000 -0.219 -0.478 -0.107
```

``` r
getFE(mod6_GEI,round=2,p.round = 5)
```

```
##                         Est.   SE    df     t       p    LL   UL
## (Intercept)             0.16 0.08 30.90  1.93 0.06270 -0.01 0.32
## gndr.c                  0.22 0.02 15.77 13.89 0.00000  0.18 0.25
## year.c                 -0.01 0.00 23.01 -1.99 0.05918 -0.01 0.00
## gndr.c:year.c           0.00 0.00 14.09 -5.83 0.00004 -0.01 0.00
## gndr.c:gei.z.cm        -0.02 0.01 35.12 -1.80 0.08019 -0.05 0.00
## year.c:gei.z.cm        -0.01 0.00 35.31 -2.55 0.01525 -0.01 0.00
## gndr.c:year.c:gei.z.cm  0.00 0.00 27.13  2.44 0.02128  0.00 0.00
```

``` r
getVC(mod6_GEI)
```

```
##         grp          var1          var2 sdcor  vcov
## 1     cntry   (Intercept)          <NA>  0.46  0.21
## 2     cntry        gndr.c          <NA>  0.08  0.01
## 3     cntry        year.c          <NA>  0.02  0.00
## 4     cntry gndr.c:year.c          <NA>  0.00  0.00
## 5     cntry   (Intercept)        gndr.c -0.72 -0.03
## 6     cntry   (Intercept)        year.c -0.74 -0.01
## 7     cntry   (Intercept) gndr.c:year.c  0.76  0.00
## 8     cntry        gndr.c        year.c  0.60  0.00
## 9     cntry        gndr.c gndr.c:year.c -0.65  0.00
## 10    cntry        year.c gndr.c:year.c -0.67  0.00
## 11 Residual          <NA>          <NA>  1.00  1.00
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_GEI,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] -26.23999
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_GEI,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] -59.6493
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
##      21  0.1412 0.0761 37.9 -0.01293   0.2954   1.855  0.0714
##       0  0.1551 0.0803 30.9 -0.00874   0.3190   1.931  0.0627
## 
## gei.z.cm =  0:
##  year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21  0.0109 0.0555 28.8 -0.10267   0.1244   0.196  0.8460
##       0  0.1551 0.0803 30.9 -0.00874   0.3190   1.931  0.0627
## 
## gei.z.cm =  1:
##  year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21 -0.1195 0.0748 35.6 -0.27114   0.0322  -1.598  0.1188
##       0  0.1551 0.0803 30.9 -0.00874   0.3190   1.931  0.0627
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
##  year.c21 - year.c0  -0.0139 0.0896 35.5   -0.196  0.16800  -0.155  0.8778
## 
## gei.z.cm =  0:
##  contrast           estimate     SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0  -0.1442 0.0727 23.0   -0.295  0.00607  -1.985  0.0592
## 
## gei.z.cm =  1:
##  contrast           estimate     SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0  -0.2746 0.0880 32.9   -0.454 -0.09546  -3.119  0.0038
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
##    -0.5     21  0.0842 0.0775 38.2  -0.0727   0.2410   1.086  0.2842
##     0.5     21  0.1983 0.0761 37.8   0.0441   0.3524   2.605  0.0131
##    -0.5      0  0.0351 0.0859 30.8  -0.1402   0.2103   0.408  0.6857
##     0.5      0  0.2751 0.0757 30.4   0.1206   0.4297   3.633  0.0010
## 
## gei.z.cm =  0:
##  gndr.c year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21 -0.0543 0.0566 28.7  -0.1701   0.0614  -0.960  0.3450
##     0.5     21  0.0761 0.0551 29.0  -0.0367   0.1888   1.380  0.1781
##    -0.5      0  0.0460 0.0857 30.5  -0.1288   0.2209   0.537  0.5950
##     0.5      0  0.2642 0.0754 29.9   0.1101   0.4183   3.503  0.0015
## 
## gei.z.cm =  1:
##  gndr.c year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21 -0.1928 0.0760 35.8  -0.3469  -0.0387  -2.538  0.0156
##     0.5     21 -0.0461 0.0746 35.4  -0.1976   0.1053  -0.618  0.5406
##    -0.5      0  0.0569 0.0859 30.7  -0.1182   0.2321   0.663  0.5121
##     0.5      0  0.2533 0.0756 30.2   0.0989   0.4077   3.349  0.0022
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
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21    -0.1141 0.0203 38.1  -0.1552  -0.0730  -5.618 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)   0.0491 0.0936 35.0  -0.1410   0.2392   0.524  0.6034
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     -0.1910 0.0868 34.7  -0.3673  -0.0146  -2.199  0.0346
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)      0.1632 0.0953 35.6  -0.0301   0.3564   1.713  0.0954
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0        -0.0769 0.0871 34.0  -0.2539   0.1002  -0.882  0.3838
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0      -0.2401 0.0205 30.3  -0.2818  -0.1983 -11.737 <0.0001
## 
## gei.z.cm =  0:
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21    -0.1304 0.0126 31.2  -0.1561  -0.1048 -10.364 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)  -0.1003 0.0764 22.8  -0.2585   0.0578  -1.313  0.2023
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     -0.3185 0.0685 22.4  -0.4605  -0.1766  -4.649  0.0001
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)      0.0301 0.0785 23.3  -0.1321   0.1923   0.383  0.7050
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0        -0.1881 0.0695 22.0  -0.3323  -0.0440  -2.706  0.0129
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0      -0.2182 0.0157 15.8  -0.2515  -0.1849 -13.889 <0.0001
## 
## gei.z.cm =  1:
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21    -0.1467 0.0180 32.0  -0.1834  -0.1101  -8.151 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)  -0.2498 0.0918 32.3  -0.4367  -0.0628  -2.720  0.0104
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     -0.4461 0.0850 31.8  -0.6194  -0.2729  -5.247 <0.0001
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)     -0.1031 0.0937 33.2  -0.2937   0.0876  -1.100  0.2793
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0        -0.2994 0.0852 31.1  -0.4731  -0.1257  -3.516  0.0014
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0      -0.1964 0.0192 24.3  -0.2360  -0.1567 -10.216 <0.0001
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
##  diff_ESS11    0.114 0.0203 38.1    0.073    0.155   5.618 <0.0001
##  diff_ESS1     0.240 0.0205 30.3    0.198    0.282  11.737 <0.0001
## 
## gei.z.cm =  0:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11    0.130 0.0126 31.2    0.105    0.156  10.364 <0.0001
##  diff_ESS1     0.218 0.0157 15.8    0.185    0.252  13.889 <0.0001
## 
## gei.z.cm =  1:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11    0.147 0.0180 32.0    0.110    0.183   8.151 <0.0001
##  diff_ESS1     0.196 0.0192 24.3    0.157    0.236  10.216 <0.0001
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
##  diff_ESS11 - diff_ESS1  -0.1260 0.0240 28.5  -0.1750 -0.07693  -5.259 <0.0001
## 
## gei.z.cm =  0:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1  -0.0878 0.0151 14.1  -0.1201 -0.05550  -5.828 <0.0001
## 
## gei.z.cm =  1:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1  -0.0496 0.0192 18.4  -0.0899 -0.00942  -2.589  0.0183
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```


#### Country-level predictors: GGGI


``` r
mod6_GGGI<-lmer(ach.z~gndr.c+year.c+
             gndr.c:year.c+
               gggi.z.cm:gndr.c+gggi.z.cm:year.c+gggi.z.cm:gndr.c:year.c+
               (gndr.c+year.c+gndr.c:year.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_GGGI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ach.z ~ gndr.c + year.c + gndr.c:year.c + gggi.z.cm:gndr.c +  
##     gggi.z.cm:year.c + gggi.z.cm:gndr.c:year.c + (gndr.c + year.c +      gndr.c:year.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1078175.3 1078369.8 -539069.6 1078139.3    363834 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.0865 -0.6671  0.0244  0.6578  5.1270 
## 
## Random effects:
##  Groups   Name          Variance  Std.Dev. Corr              
##  cntry    (Intercept)   4.676e-01 0.68384                    
##           gndr.c        1.752e-02 0.13235  -0.90             
##           year.c        1.081e-03 0.03287  -0.93  0.86       
##           gndr.c:year.c 2.852e-05 0.00534   0.91 -0.87 -0.90 
##  Residual               9.938e-01 0.99687                    
## Number of obs: 363852, groups:  cntry, 34
## 
## Fixed effects:
##                          Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)              0.244367   0.119273 28.019189   2.049 0.049952 *  
## gndr.c                   0.199352   0.024519  7.046063   8.131 7.93e-05 ***
## year.c                  -0.012444   0.005764 25.513475  -2.159 0.040430 *  
## gndr.c:year.c           -0.003006   0.001139  7.822828  -2.640 0.030315 *  
## gndr.c:gggi.z.cm        -0.031290   0.013998 35.787088  -2.235 0.031727 *  
## year.c:gggi.z.cm        -0.009181   0.002366 33.778779  -3.881 0.000459 ***
## gndr.c:year.c:gggi.z.cm  0.001566   0.000873 32.647757   1.793 0.082188 .  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c year.c gnd.:. gn.:.. yr.:..
## gndr.c      -0.848                                   
## year.c      -0.927  0.807                            
## gndr.c:yr.c  0.743 -0.852 -0.735                     
## gndr.c:gg.. -0.002 -0.016  0.005  0.029              
## yr.c:ggg.z.  0.012 -0.005 -0.039  0.013  0.074       
## gndr.c:.:.. -0.001  0.023  0.003 -0.060 -0.650 -0.163
```

``` r
getFE(mod6_GGGI,round=2,p.round = 5)
```

```
##                          Est.   SE    df     t       p    LL   UL
## (Intercept)              0.24 0.12 28.02  2.05 0.04995  0.00 0.49
## gndr.c                   0.20 0.02  7.05  8.13 0.00008  0.14 0.26
## year.c                  -0.01 0.01 25.51 -2.16 0.04043 -0.02 0.00
## gndr.c:year.c            0.00 0.00  7.82 -2.64 0.03032 -0.01 0.00
## gndr.c:gggi.z.cm        -0.03 0.01 35.79 -2.24 0.03173 -0.06 0.00
## year.c:gggi.z.cm        -0.01 0.00 33.78 -3.88 0.00046 -0.01 0.00
## gndr.c:year.c:gggi.z.cm  0.00 0.00 32.65  1.79 0.08219  0.00 0.00
```

``` r
getVC(mod6_GGGI)
```

```
##         grp          var1          var2 sdcor  vcov
## 1     cntry   (Intercept)          <NA>  0.68  0.47
## 2     cntry        gndr.c          <NA>  0.13  0.02
## 3     cntry        year.c          <NA>  0.03  0.00
## 4     cntry gndr.c:year.c          <NA>  0.01  0.00
## 5     cntry   (Intercept)        gndr.c -0.90 -0.08
## 6     cntry   (Intercept)        year.c -0.93 -0.02
## 7     cntry   (Intercept) gndr.c:year.c  0.91  0.00
## 8     cntry        gndr.c        year.c  0.86  0.00
## 9     cntry        gndr.c gndr.c:year.c -0.87  0.00
## 10    cntry        year.c gndr.c:year.c -0.90  0.00
## 11 Residual          <NA>          <NA>  1.00  0.99
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_GGGI,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] -251.3552
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_GGGI,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] -471.483
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
##  year.c emmean     SE   df  lower.CL upper.CL t.ratio p.value
##      21  0.176 0.0701 34.3  3.35e-02   0.3182   2.510  0.0170
##       0  0.244 0.1190 28.0  5.52e-05   0.4887   2.049  0.0500
## 
## gggi.z.cm =  0:
##  year.c emmean     SE   df  lower.CL upper.CL t.ratio p.value
##      21 -0.017 0.0460 33.3 -1.11e-01   0.0766  -0.369  0.7146
##       0  0.244 0.1190 28.0  5.52e-05   0.4887   2.049  0.0500
## 
## gggi.z.cm =  1:
##  year.c emmean     SE   df  lower.CL upper.CL t.ratio p.value
##      21 -0.210 0.0653 33.8 -3.42e-01  -0.0771  -3.214  0.0029
##       0  0.244 0.1190 28.0  5.52e-05   0.4887   2.049  0.0500
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
##  year.c21 - year.c0  -0.0685 0.133 33.8   -0.338   0.2010  -0.517  0.6087
## 
## gggi.z.cm =  0:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0  -0.2613 0.121 25.5   -0.510  -0.0123  -2.159  0.0404
## 
## gggi.z.cm =  1:
##  contrast           estimate    SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0  -0.4541 0.129 31.4   -0.717  -0.1911  -3.519  0.0013
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
##    -0.5     21  0.1085 0.0722 34.3  -0.0382  0.25526   1.503  0.1421
##     0.5     21  0.2432 0.0693 34.5   0.1025  0.38387   3.510  0.0013
##    -0.5      0  0.1290 0.1300 27.4  -0.1375  0.39563   0.992  0.3296
##     0.5      0  0.3597 0.1090 26.6   0.1352  0.58414   3.290  0.0028
## 
## gggi.z.cm =  0:
##  gndr.c year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21 -0.0851 0.0475 33.1  -0.1816  0.01146  -1.793  0.0822
##     0.5     21  0.0512 0.0455 33.5  -0.0413  0.14360   1.125  0.2686
##    -0.5      0  0.1447 0.1300 27.3  -0.1216  0.41096   1.114  0.2748
##     0.5      0  0.3440 0.1090 26.4   0.1200  0.56808   3.154  0.0040
## 
## gggi.z.cm =  1:
##  gndr.c year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21 -0.2787 0.0672 33.9  -0.4153 -0.14207  -4.146  0.0002
##     0.5     21 -0.1409 0.0646 33.8  -0.2723 -0.00947  -2.179  0.0364
##    -0.5      0  0.1603 0.1300 27.4  -0.1063  0.42692   1.233  0.2280
##     0.5      0  0.3284 0.1090 26.6   0.1040  0.55276   3.006  0.0057
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
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21    -0.1346 0.0199 33.39   -0.175 -0.09415  -6.763 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)  -0.0205 0.1420 32.56   -0.309  0.26822  -0.145  0.8859
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     -0.2512 0.1240 32.64   -0.503  0.00104  -2.027  0.0509
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)      0.1141 0.1430 32.79   -0.177  0.40574   0.796  0.4315
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0        -0.1165 0.1250 31.73   -0.370  0.13736  -0.935  0.3568
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0      -0.2306 0.0284 11.75   -0.293 -0.16857  -8.114 <0.0001
## 
## gggi.z.cm =  0:
##  contrast                                     estimate     SE    df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21    -0.1362 0.0132 30.28   -0.163 -0.10926 -10.313 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)  -0.2298 0.1300 24.83   -0.498  0.03822  -1.766  0.0896
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     -0.4291 0.1110 23.88   -0.657 -0.20081  -3.880  0.0007
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)     -0.0935 0.1330 25.39   -0.366  0.17915  -0.706  0.4867
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0        -0.2929 0.1130 23.91   -0.525 -0.06057  -2.603  0.0156
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0      -0.1994 0.0245  7.05   -0.257 -0.14145  -8.131 <0.0001
## 
## gggi.z.cm =  1:
##  contrast                                     estimate     SE    df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21    -0.1378 0.0187 32.35   -0.176 -0.09978  -7.377 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)  -0.4390 0.1380 30.18   -0.721 -0.15670  -3.175  0.0034
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     -0.6071 0.1200 30.32   -0.852 -0.36245  -5.066 <0.0001
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)     -0.3012 0.1400 30.53   -0.587 -0.01510  -2.149  0.0397
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0        -0.4693 0.1210 29.70   -0.716 -0.22226  -3.882  0.0005
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0      -0.1681 0.0280 11.09   -0.230 -0.10640  -5.994 <0.0001
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
##  diff_ESS11    0.135 0.0199 33.39   0.0942    0.175   6.763 <0.0001
##  diff_ESS1     0.231 0.0284 11.75   0.1686    0.293   8.114 <0.0001
## 
## gggi.z.cm =  0:
##  contrast   estimate     SE    df lower.CL upper.CL t.ratio p.value
##  diff_ESS11    0.136 0.0132 30.28   0.1093    0.163  10.313 <0.0001
##  diff_ESS1     0.199 0.0245  7.05   0.1415    0.257   8.131 <0.0001
## 
## gggi.z.cm =  1:
##  contrast   estimate     SE    df lower.CL upper.CL t.ratio p.value
##  diff_ESS11    0.138 0.0187 32.35   0.0998    0.176   7.377 <0.0001
##  diff_ESS1     0.168 0.0280 11.09   0.1064    0.230   5.994 <0.0001
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
##  contrast               estimate     SE    df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1  -0.0960 0.0310 18.31  -0.1610 -0.03096  -3.097  0.0061
## 
## gggi.z.cm =  0:
##  contrast               estimate     SE    df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1  -0.0631 0.0239  7.82  -0.1185 -0.00776  -2.640  0.0303
## 
## gggi.z.cm =  1:
##  contrast               estimate     SE    df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1  -0.0302 0.0292 15.19  -0.0925  0.03202  -1.034  0.3172
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

#### Country-level predictors: GDI


``` r
mod6_GDI<-lmer(ach.z~gndr.c+year.c+
             gndr.c:year.c+
               gdi.z.cm:gndr.c+gdi.z.cm:year.c+gdi.z.cm:gndr.c:year.c+
               (gndr.c+year.c+gndr.c:year.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_GDI)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ach.z ~ gndr.c + year.c + gndr.c:year.c + gdi.z.cm:gndr.c + gdi.z.cm:year.c +  
##     gdi.z.cm:gndr.c:year.c + (gndr.c + year.c + gndr.c:year.c |      cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
## 1460674.4 1460874.3 -730319.2 1460638.4    492325 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.0636 -0.6694  0.0263  0.6588  5.1222 
## 
## Random effects:
##  Groups   Name          Variance  Std.Dev. Corr              
##  cntry    (Intercept)   2.032e-01 0.450734                   
##           gndr.c        4.734e-03 0.068808 -0.59             
##           year.c        3.025e-04 0.017393 -0.66  0.41       
##           gndr.c:year.c 5.073e-06 0.002252  0.42 -0.44 -0.36 
##  Residual               9.991e-01 0.999565                   
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##                          Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)             0.1440640  0.0775534 31.8510828   1.858   0.0725 .  
## gndr.c                  0.2194762  0.0133707 24.7028643  16.415 8.53e-15 ***
## year.c                 -0.0064563  0.0030152 27.2676819  -2.141   0.0413 *  
## gndr.c:year.c          -0.0042005  0.0006284 21.1960784  -6.685 1.23e-06 ***
## gndr.c:gdi.z.cm        -0.0162782  0.0122144 38.5353623  -1.333   0.1905    
## year.c:gdi.z.cm        -0.0046808  0.0023640 34.9013987  -1.980   0.0556 .  
## gndr.c:year.c:gdi.z.cm  0.0005862  0.0007726 29.1565662   0.759   0.4541    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c year.c gnd.:. gn.:.. yr.:..
## gndr.c      -0.518                                   
## year.c      -0.662  0.356                            
## gndr.c:yr.c  0.258 -0.552 -0.222                     
## gndr.c:gd..  0.000  0.002  0.000 -0.008              
## yr.c:gd.z.c  0.005 -0.002 -0.021  0.005  0.013       
## gndr.c:.:.. -0.002 -0.011  0.004 -0.013 -0.536 -0.053
```

``` r
getFE(mod6_GDI,round=2,p.round = 5)
```

```
##                         Est.   SE    df     t       p    LL   UL
## (Intercept)             0.14 0.08 31.85  1.86 0.07249 -0.01 0.30
## gndr.c                  0.22 0.01 24.70 16.41 0.00000  0.19 0.25
## year.c                 -0.01 0.00 27.27 -2.14 0.04134 -0.01 0.00
## gndr.c:year.c           0.00 0.00 21.20 -6.68 0.00000 -0.01 0.00
## gndr.c:gdi.z.cm        -0.02 0.01 38.54 -1.33 0.19046 -0.04 0.01
## year.c:gdi.z.cm         0.00 0.00 34.90 -1.98 0.05563 -0.01 0.00
## gndr.c:year.c:gdi.z.cm  0.00 0.00 29.16  0.76 0.45411  0.00 0.00
```

``` r
getVC(mod6_GDI)
```

```
##         grp          var1          var2 sdcor  vcov
## 1     cntry   (Intercept)          <NA>  0.45  0.20
## 2     cntry        gndr.c          <NA>  0.07  0.00
## 3     cntry        year.c          <NA>  0.02  0.00
## 4     cntry gndr.c:year.c          <NA>  0.00  0.00
## 5     cntry   (Intercept)        gndr.c -0.59 -0.02
## 6     cntry   (Intercept)        year.c -0.66 -0.01
## 7     cntry   (Intercept) gndr.c:year.c  0.42  0.00
## 8     cntry        gndr.c        year.c  0.41  0.00
## 9     cntry        gndr.c gndr.c:year.c -0.44  0.00
## 10    cntry        year.c gndr.c:year.c -0.36  0.00
## 11 Residual          <NA>          <NA>  1.00  1.00
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_GDI,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] 1.640174
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_GDI,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] -1.659319
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
##      21  0.10678 0.0779 42.0  -0.0505   0.2640   1.370  0.1778
##       0  0.14406 0.0776 31.9  -0.0139   0.3021   1.858  0.0725
## 
## gdi.z.cm =  0:
##  year.c   emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21  0.00848 0.0593 33.2  -0.1122   0.1291   0.143  0.8872
##       0  0.14406 0.0776 31.9  -0.0139   0.3021   1.858  0.0725
## 
## gdi.z.cm =  1:
##  year.c   emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21 -0.08982 0.0768 40.1  -0.2450   0.0653  -1.170  0.2490
##       0  0.14406 0.0776 31.9  -0.0139   0.3021   1.858  0.0725
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
##  year.c21 - year.c0  -0.0373 0.0813 38.3   -0.202  0.12716  -0.459  0.6489
## 
## gdi.z.cm =  0:
##  contrast           estimate     SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0  -0.1356 0.0633 27.3   -0.265 -0.00572  -2.141  0.0413
## 
## gdi.z.cm =  1:
##  contrast           estimate     SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0  -0.2339 0.0797 35.4   -0.396 -0.07224  -2.936  0.0058
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
##    -0.5     21  0.0392 0.0797 42.4  -0.1217  0.20006   0.491  0.6259
##     0.5     21  0.1744 0.0772 41.6   0.0185  0.33031   2.258  0.0293
##    -0.5      0  0.0262 0.0815 32.0  -0.1397  0.19209   0.322  0.7499
##     0.5      0  0.2619 0.0746 31.7   0.1100  0.41386   3.513  0.0014
## 
## gdi.z.cm =  0:
##  gndr.c year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21 -0.0572 0.0610 33.2  -0.1811  0.06684  -0.938  0.3552
##     0.5     21  0.0741 0.0583 33.2  -0.0445  0.19271   1.271  0.2125
##    -0.5      0  0.0343 0.0812 31.7  -0.1312  0.19983   0.423  0.6754
##     0.5      0  0.2538 0.0743 31.3   0.1023  0.40530   3.416  0.0018
## 
## gdi.z.cm =  1:
##  gndr.c year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21 -0.1535 0.0786 40.4  -0.3122  0.00525  -1.954  0.0577
##     0.5     21 -0.0262 0.0761 39.7  -0.1800  0.12767  -0.344  0.7328
##    -0.5      0  0.0425 0.0814 32.0  -0.1234  0.20836   0.521  0.6057
##     0.5      0  0.2457 0.0746 31.7   0.0937  0.39759   3.295  0.0024
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
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21    -0.1352 0.0193 34.5  -0.1743  -0.0961  -7.024 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)   0.0130 0.0834 37.7  -0.1558   0.1818   0.156  0.8772
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     -0.2228 0.0803 36.6  -0.3856  -0.0599  -2.773  0.0087
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)      0.1482 0.0850 38.9  -0.0238   0.3202   1.743  0.0892
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0        -0.0875 0.0805 36.3  -0.2507   0.0756  -1.088  0.2839
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0      -0.2358 0.0181 35.9  -0.2725  -0.1991 -13.029 <0.0001
## 
## gdi.z.cm =  0:
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21    -0.1313 0.0126 31.4  -0.1569  -0.1056 -10.434 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)  -0.0915 0.0651 26.8  -0.2251   0.0421  -1.405  0.1715
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     -0.3110 0.0610 26.6  -0.4361  -0.1858  -5.100 <0.0001
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)      0.0398 0.0675 27.4  -0.0985   0.1781   0.590  0.5602
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0        -0.1797 0.0622 26.1  -0.3075  -0.0519  -2.890  0.0077
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0      -0.2195 0.0134 24.7  -0.2470  -0.1919 -16.415 <0.0001
## 
## gdi.z.cm =  1:
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21    -0.1273 0.0186 37.3  -0.1650  -0.0896  -6.845 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)  -0.1959 0.0817 34.9  -0.3618  -0.0301  -2.398  0.0219
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     -0.3991 0.0786 33.6  -0.5590  -0.2393  -5.077 <0.0001
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)     -0.0686 0.0835 36.1  -0.2379   0.1006  -0.822  0.4163
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0        -0.2718 0.0789 33.8  -0.4323  -0.1114  -3.444  0.0015
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0      -0.2032 0.0181 39.1  -0.2399  -0.1665 -11.211 <0.0001
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
##  diff_ESS11    0.135 0.0193 34.5   0.0961    0.174   7.024 <0.0001
##  diff_ESS1     0.236 0.0181 35.9   0.1991    0.272  13.029 <0.0001
## 
## gdi.z.cm =  0:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11    0.131 0.0126 31.4   0.1056    0.157  10.434 <0.0001
##  diff_ESS1     0.219 0.0134 24.7   0.1919    0.247  16.415 <0.0001
## 
## gdi.z.cm =  1:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11    0.127 0.0186 37.3   0.0896    0.165   6.845 <0.0001
##  diff_ESS1     0.203 0.0181 39.1   0.1665    0.240  11.211 <0.0001
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
##  diff_ESS11 - diff_ESS1  -0.1005 0.0211 22.7   -0.144  -0.0569  -4.775 <0.0001
## 
## gdi.z.cm =  0:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1  -0.0882 0.0132 21.2   -0.116  -0.0608  -6.685 <0.0001
## 
## gdi.z.cm =  1:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1  -0.0759 0.0208 32.3   -0.118  -0.0336  -3.653  0.0009
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

#### Country-level predictors: log_GDP


``` r
mod6_log_GDP<-lmer(ach.z~gndr.c+year.c+
             gndr.c:year.c+
               log_gdp.z.cm:gndr.c+log_gdp.z.cm:year.c+log_gdp.z.cm:gndr.c:year.c+
               (gndr.c+year.c+gndr.c:year.c|cntry),
           data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_log_GDP)
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: ach.z ~ gndr.c + year.c + gndr.c:year.c + log_gdp.z.cm:gndr.c +  
##     log_gdp.z.cm:year.c + log_gdp.z.cm:gndr.c:year.c + (gndr.c +      year.c + gndr.c:year.c | cntry)
##    Data: diff_dat
## Weights: pspwght
## Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+07))
## 
##       AIC       BIC    logLik -2*log(L)  df.resid 
##   1460678   1460878   -730321   1460642    492325 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.0740 -0.6693  0.0265  0.6587  5.1219 
## 
## Random effects:
##  Groups   Name          Variance  Std.Dev. Corr              
##  cntry    (Intercept)   2.049e-01 0.452704                   
##           gndr.c        5.548e-03 0.074482 -0.63             
##           year.c        3.214e-04 0.017929 -0.65  0.50       
##           gndr.c:year.c 5.587e-06 0.002364  0.66 -0.57 -0.60 
##  Residual               9.991e-01 0.999568                   
## Number of obs: 492343, groups:  cntry, 34
## 
## Fixed effects:
##                              Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)                 0.1447714  0.0778931 31.7581274   1.859   0.0724 .  
## gndr.c                      0.2192138  0.0143199 20.1614210  15.308 1.45e-12 ***
## year.c                     -0.0065742  0.0031056 21.9490371  -2.117   0.0458 *  
## gndr.c:year.c              -0.0042097  0.0006483 18.4033650  -6.493 3.74e-06 ***
## gndr.c:log_gdp.z.cm        -0.0097520  0.0123622 34.6780213  -0.789   0.4356    
## year.c:log_gdp.z.cm        -0.0010207  0.0023841 31.3311910  -0.428   0.6715    
## gndr.c:year.c:log_gdp.z.cm  0.0009059  0.0006231 25.1066927   1.454   0.1584    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) gndr.c year.c gnd.:. g.:_.. y.:_..
## gndr.c      -0.564                                   
## year.c      -0.649  0.446                            
## gndr.c:yr.c  0.414 -0.623 -0.375                     
## gndr.c:l_.. -0.001 -0.090  0.005  0.134              
## yr.c:lg_g..  0.008  0.000  0.001  0.000  0.119       
## gndr.:.:_.. -0.001  0.119 -0.003 -0.191 -0.550 -0.146
```

``` r
getFE(mod6_log_GDP,round=2,p.round = 5)
```

```
##                             Est.   SE    df     t       p    LL   UL
## (Intercept)                 0.14 0.08 31.76  1.86 0.07237 -0.01 0.30
## gndr.c                      0.22 0.01 20.16 15.31 0.00000  0.19 0.25
## year.c                     -0.01 0.00 21.95 -2.12 0.04584 -0.01 0.00
## gndr.c:year.c               0.00 0.00 18.40 -6.49 0.00000 -0.01 0.00
## gndr.c:log_gdp.z.cm        -0.01 0.01 34.68 -0.79 0.43556 -0.03 0.02
## year.c:log_gdp.z.cm         0.00 0.00 31.33 -0.43 0.67148 -0.01 0.00
## gndr.c:year.c:log_gdp.z.cm  0.00 0.00 25.11  1.45 0.15841  0.00 0.00
```

``` r
getVC(mod6_log_GDP)
```

```
##         grp          var1          var2 sdcor  vcov
## 1     cntry   (Intercept)          <NA>  0.45  0.20
## 2     cntry        gndr.c          <NA>  0.07  0.01
## 3     cntry        year.c          <NA>  0.02  0.00
## 4     cntry gndr.c:year.c          <NA>  0.00  0.00
## 5     cntry   (Intercept)        gndr.c -0.63 -0.02
## 6     cntry   (Intercept)        year.c -0.65 -0.01
## 7     cntry   (Intercept) gndr.c:year.c  0.66  0.00
## 8     cntry        gndr.c        year.c  0.50  0.00
## 9     cntry        gndr.c gndr.c:year.c -0.57  0.00
## 10    cntry        year.c gndr.c:year.c -0.60  0.00
## 11 Residual          <NA>          <NA>  1.00  1.00
```

``` r
# variance accounted in general slopes
100*(getVC(mod6,round=10)[3,]$vcov-getVC(mod6_log_GDP,round=10)[3,]$vcov)/
  getVC(mod6,round=10)[3,]$vcov
```

```
## [1] -4.508105
```

``` r
# variance accounted in gendered slopes
100*(getVC(mod6,round=10)[4,]$vcov-getVC(mod6_log_GDP,round=10)[4,]$vcov)/
  getVC(mod6,round=10)[4,]$vcov
```

```
## [1] -11.95391
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
##      21  0.02815 0.0785 35.6  -0.1311    0.187   0.359  0.7219
##       0  0.14477 0.0779 31.8  -0.0139    0.303   1.859  0.0724
## 
## log_gdp.z.cm =  0:
##  year.c   emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21  0.00671 0.0610 26.9  -0.1185    0.132   0.110  0.9132
##       0  0.14477 0.0779 31.8  -0.0139    0.303   1.859  0.0724
## 
## log_gdp.z.cm =  1:
##  year.c   emmean     SE   df lower.CL upper.CL t.ratio p.value
##      21 -0.01472 0.0794 35.1  -0.1759    0.146  -0.185  0.8539
##       0  0.14477 0.0779 31.8  -0.0139    0.303   1.859  0.0724
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
##  year.c21 - year.c0   -0.117 0.0822 32.2   -0.284  0.05073  -1.419  0.1655
## 
## log_gdp.z.cm =  0:
##  contrast           estimate     SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0   -0.138 0.0652 21.9   -0.273 -0.00279  -2.117  0.0458
## 
## log_gdp.z.cm =  1:
##  contrast           estimate     SE   df lower.CL upper.CL t.ratio p.value
##  year.c21 - year.c0   -0.159 0.0823 30.1   -0.327  0.00845  -1.939  0.0619
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
##    -0.5     21 -0.0326 0.0798 35.6  -0.1945   0.1292  -0.409  0.6850
##     0.5     21  0.0889 0.0781 35.5  -0.0696   0.2475   1.138  0.2627
##    -0.5      0  0.0303 0.0824 31.8  -0.1376   0.1982   0.367  0.7157
##     0.5      0  0.2593 0.0744 31.7   0.1076   0.4109   3.484  0.0015
## 
## log_gdp.z.cm =  0:
##  gndr.c year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21 -0.0587 0.0622 26.9  -0.1863   0.0689  -0.944  0.3537
##     0.5     21  0.0721 0.0604 27.0  -0.0519   0.1961   1.193  0.2431
##    -0.5      0  0.0352 0.0821 31.4  -0.1323   0.2026   0.428  0.6715
##     0.5      0  0.2544 0.0741 31.1   0.1033   0.4055   3.433  0.0017
## 
## log_gdp.z.cm =  1:
##  gndr.c year.c  emmean     SE   df lower.CL upper.CL t.ratio p.value
##    -0.5     21 -0.0848 0.0806 35.1  -0.2485   0.0789  -1.051  0.3004
##     0.5     21  0.0553 0.0790 34.8  -0.1051   0.2157   0.700  0.4884
##    -0.5      0  0.0400 0.0823 31.7  -0.1277   0.2078   0.486  0.6301
##     0.5      0  0.2495 0.0743 31.5   0.0981   0.4009   3.359  0.0021
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
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21    -0.1215 0.0174 35.5   -0.157 -0.08619  -6.976 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)  -0.0629 0.0854 31.9   -0.237  0.11099  -0.737  0.4665
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     -0.2919 0.0803 31.6   -0.456 -0.12822  -3.635  0.0010
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)      0.0586 0.0868 32.5   -0.118  0.23539   0.675  0.5043
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0        -0.1703 0.0802 31.1   -0.334 -0.00677  -2.124  0.0418
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0      -0.2290 0.0197 36.0   -0.269 -0.18893 -11.600 <0.0001
## 
## log_gdp.z.cm =  0:
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21    -0.1308 0.0122 31.1   -0.156 -0.10603 -10.765 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)  -0.0939 0.0681 21.7   -0.235  0.04741  -1.379  0.1820
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     -0.3131 0.0623 21.9   -0.442 -0.18378  -5.023 <0.0001
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)      0.0370 0.0699 22.1   -0.108  0.18184   0.529  0.6022
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0        -0.1823 0.0630 21.3   -0.313 -0.05141  -2.894  0.0086
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0      -0.2192 0.0143 20.2   -0.249 -0.18936 -15.308 <0.0001
## 
## log_gdp.z.cm =  1:
##  contrast                                     estimate     SE   df lower.CL upper.CL t.ratio p.value
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c21    -0.1401 0.0168 30.5   -0.174 -0.10570  -8.316 <0.0001
##  (gndr.c-0.5 year.c21) - (gndr.c-0.5 year.c0)  -0.1248 0.0853 29.7   -0.299  0.04939  -1.464  0.1537
##  (gndr.c-0.5 year.c21) - gndr.c0.5 year.c0     -0.3343 0.0803 29.3   -0.498 -0.17017  -4.164  0.0003
##  gndr.c0.5 year.c21 - (gndr.c-0.5 year.c0)      0.0153 0.0869 30.6   -0.162  0.19264   0.176  0.8616
##  gndr.c0.5 year.c21 - gndr.c0.5 year.c0        -0.1942 0.0800 28.8   -0.358 -0.03043  -2.426  0.0218
##  (gndr.c-0.5 year.c0) - gndr.c0.5 year.c0      -0.2095 0.0181 27.5   -0.246 -0.17244 -11.599 <0.0001
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
##  diff_ESS11    0.122 0.0174 35.5   0.0862    0.157   6.976 <0.0001
##  diff_ESS1     0.229 0.0197 36.0   0.1889    0.269  11.600 <0.0001
## 
## log_gdp.z.cm =  0:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11    0.131 0.0122 31.1   0.1060    0.156  10.765 <0.0001
##  diff_ESS1     0.219 0.0143 20.2   0.1894    0.249  15.308 <0.0001
## 
## log_gdp.z.cm =  1:
##  contrast   estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11    0.140 0.0168 30.5   0.1057    0.174   8.316 <0.0001
##  diff_ESS1     0.209 0.0181 27.5   0.1724    0.246  11.599 <0.0001
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
##  diff_ESS11 - diff_ESS1  -0.1074 0.0206 29.8   -0.150  -0.0653  -5.214 <0.0001
## 
## log_gdp.z.cm =  0:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1  -0.0884 0.0136 18.4   -0.117  -0.0598  -6.493 <0.0001
## 
## log_gdp.z.cm =  1:
##  contrast               estimate     SE   df lower.CL upper.CL t.ratio p.value
##  diff_ESS11 - diff_ESS1  -0.0694 0.0170 18.3   -0.105  -0.0337  -4.083  0.0007
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

