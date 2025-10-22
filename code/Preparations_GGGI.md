---
title: "Data preparations for GGGI country indices"
output: 
  html_document: 
    toc: true
    keep_md: true
---





# Packages


``` r
library(rio)
library(dplyr)
library(tidyr)
library(psych)
```

# Read Data


``` r
GGGI_d<-import("../data/qogdata_03_10_2025.xlsx")
```

```
## New names:
## • `` -> `...1`
```

``` r
# read also ISO to obtain ISO2 codes
ISO<-read.csv2("../data/ISO.csv")
```

## GGGI (Global Gender Gap Index)


``` r
# select relevant variables
GGGI_d<-GGGI_d %>%
  dplyr::select(cname,year,ccodealp,gggi_ggi)

# make a wide data format
GGGI_wide <- GGGI_d %>%
  arrange(year) %>%
  pivot_wider(
    id_cols = c(cname, ccodealp),
    names_from = year,
    values_from = gggi_ggi,
    names_prefix = "gggi_"
  )

# add ISO

GGGI_wide <- left_join(x=GGGI_wide,
                   y=ISO[,c("ISO2","ISO3")],
                   by=c("ccodealp"="ISO3"))


# check which years are available
table(GGGI_d$year)
```

```
## 
## 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017 2019 2020 2021 2022 2023 
##  115  128  130  133  134  135  135  136  141  144  144  144  149  153  156  146  146  146
```

``` r
names(GGGI_wide)
```

```
##  [1] "cname"     "ccodealp"  "gggi_2005" "gggi_2006" "gggi_2007" "gggi_2008" "gggi_2009" "gggi_2010" "gggi_2011" "gggi_2012" "gggi_2013"
## [12] "gggi_2014" "gggi_2015" "gggi_2016" "gggi_2017" "gggi_2019" "gggi_2020" "gggi_2021" "gggi_2022" "gggi_2023" "ISO2"
```

``` r
# calculate GGGI mean across the time period 2002-2022 including the middle years

GGGI_vars_2002_2022<-
  paste0("gggi_",2002:2022)
GGGI_vars_2002_2022
```

```
##  [1] "gggi_2002" "gggi_2003" "gggi_2004" "gggi_2005" "gggi_2006" "gggi_2007" "gggi_2008" "gggi_2009" "gggi_2010" "gggi_2011" "gggi_2012"
## [12] "gggi_2013" "gggi_2014" "gggi_2015" "gggi_2016" "gggi_2017" "gggi_2018" "gggi_2019" "gggi_2020" "gggi_2021" "gggi_2022"
```

``` r
table(GGGI_vars_2002_2022 %in% names(GGGI_wide))
```

```
## 
## FALSE  TRUE 
##     4    17
```

``` r
# which variables are not available
GGGI_vars_2002_2022[!(GGGI_vars_2002_2022 %in% names(GGGI_wide))]
```

```
## [1] "gggi_2002" "gggi_2003" "gggi_2004" "gggi_2018"
```

``` r
# code only available ones
GGGI_vars_2002_2022_available<-
  GGGI_vars_2002_2022[(GGGI_vars_2002_2022 %in% names(GGGI_wide))]

GGGI_wide$GGGI_2002_2022_avg<-
  rowMeans(GGGI_wide[,GGGI_vars_2002_2022_available],na.rm=T)
  
describe(GGGI_wide$GGGI_2002_2022_avg)
```

```
##    vars   n mean   sd median trimmed  mad  min  max range  skew kurtosis se
## X1    1 157 0.69 0.06   0.69    0.69 0.05 0.43 0.86  0.43 -0.53     1.95  0
```

``` r
rio::export(data.frame(GGGI_wide),"../data/GGGI.xlsx",overwrite=T)
```



