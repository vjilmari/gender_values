#' ---
#' title: "Examinations of VBMT reliability"
#' output: 
#'   html_document: 
#'     toc: true
#'     keep_md: true
#' ---
#' 
## ----setup, include=FALSE-------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

#' 
#' # Packages
#' 
#' Two packages are not in CRAN and need different installation
#' * devtools::install_github("vjilmari/vjihelpers")
#' * install.packages("ggflags", repos = c("https://jimjam-slam.r-universe.dev","https://cloud.r-project.org")) 
#' 
## ----message=FALSE--------------------------------------------------------------------------------------

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



#' 
#' # Custom functions
#' 
## -------------------------------------------------------------------------------------------------------
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


#' 
#' 
#' # Data
#' 
#' * Loads the preprocessed European Social Survey data file made within "Preparations_ESS" code
#' 
## -------------------------------------------------------------------------------------------------------
load("../data/fdat.rdata")


#' 
#' * Include ISO data for country-names
#' 
## -------------------------------------------------------------------------------------------------------
# read country names
ISO<-read.csv2("../data/ISO.csv")

#' 
#' 
#' # Exclusions
#' 
#' * exclude participants with missing values or gender
#' 
## -------------------------------------------------------------------------------------------------------
value.vars<-
  c("con","tra",
  "ben","uni",
  "sdi","sti",
  "hed","ach",
  "pow","sec")

fdat$miss_values<-
  rowSums(is.na(fdat[,value.vars]))
table(fdat$miss_values)

fdat<-fdat %>%
  filter(miss_values==0 & !is.na(gndr.bin))

#' 
#' # Construct value-based gender-typicality measure
#' 
#' * Use 100 men and women from each country x time data fold for training set
#' 
## -------------------------------------------------------------------------------------------------------
# set seed number for reproducibility in training/testing split
set.seed(13032023)
# check the cross-validation fold sizes
table(fdat$cntry_time)
# run the analysis
value_typ<-
  D_regularized(data=fdat,mv.vars = value.vars,
                group.var = "gndr.bin",group.values = c(1,0),
                out = T,fold = T,fold.var = "cntry_time",size = 100,
                pcc = T,auc=T,pred.prob = T,append.data=T)

#' ## Gender typicality results and description
#' 
#' 
#' ### Training sample size
#' 
## -------------------------------------------------------------------------------------------------------

# n folds in the training phase
length(unique(fdat$cntry_time))
length(unique(fdat$cntry_time))*200


#' 
#' ### Coefficients for value variables in training set
#' 
## -------------------------------------------------------------------------------------------------------
round(coefficients(value_typ$cv.mod,s = "lambda.min"),2)
plot(value_typ$cv.mod)

#' 
#' ### Description of gender differences/prediction in testing dataset
#' 
## -------------------------------------------------------------------------------------------------------
# save to summary tab
sum_tab<-value_typ$D
# print the country X time -fold results
round(sum_tab,2)
# range in gender differences across folds with fold-specific SDs
range(sum_tab$D)
# range in gender differences across folds with SD pooled across all folds
range(sum_tab$d.sd.total)
# mean gender difference across folds with SD pooled across all folds
mean(sum_tab$d.sd.total)
# average probability of correct classification (pcc)
mean(sum_tab$pcc.total)
# pcc for men
mean(sum_tab$pcc.1)
# pcc for women
mean(sum_tab$pcc.0)
# average area under the curve
mean(sum_tab$auc)
# print smallest and largest gender differences
sum_tab[sum_tab$d.sd.total==min(sum_tab$d.sd.total),]
sum_tab[sum_tab$d.sd.total==max(sum_tab$d.sd.total),]
# correlation between men and women in male-typicality across the folds
cor(sum_tab$m.1,sum_tab$m.0)
# correlation between men and women in male-typicality deviations across the folds
cor(sum_tab$sd.1,sum_tab$sd.0)

#' 
#' # Data for the analysis
#' 
#' * use the testing partition of the data set (diff_dat)
#' * basic descriptives
#' 
## -------------------------------------------------------------------------------------------------------

diff_dat<-
  value_typ$preds

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
# range of sample sizes
range(table(diff_dat$cntry))

# value-based gender-typicality histogram
hist(diff_dat$pred)
# scale/standardize with SD pooled across all country-time-gender folds
FM_pooled_sd<-value_typ$D[1,"pooled.sd.total"]
FM_pooled_sd
# standardized
diff_dat$FM.z<-diff_dat$pred/FM_pooled_sd
hist(diff_dat$FM.z)


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



#' 
#' ## Add GEI-variable
#' 
#' * GEI = Reverse-coded Gender Inequality Index (GII)
#' 
#' 
## -------------------------------------------------------------------------------------------------------
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
  select(ISO2, iso3, country, year, gii, gei, gii_2002_2023_avg, gei_2002_2023_avg)


# describe gei average
psych::describe(GII_in_ESS_d$gei_2002_2023_avg)
# one is missing, see which one
GII_in_ESS_d[is.na(GII_in_ESS_d$gei_2002_2023_avg),"country"]
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


#' 
#' 
#' ## Add GDI-variable
#' 
#' * GDI = Gender Development Index
#' 
## -------------------------------------------------------------------------------------------------------
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
  select(ISO2, iso3, country, year, gdi, gdi_2002_2023_avg)

# describe gdi average
psych::describe(GDI_in_ESS_d$gdi_2002_2023_avg)
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


#' 
#' ## Add GDP-variable
#' 
#' * log_GDP = Logarithm of GDP per capita, PPP (constant 2017 international $)
#' 
#' 
## -------------------------------------------------------------------------------------------------------
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
  select(ISO2, Country.Name, year, gdp, gdp_2002_2023_avg,log_gdp_2002_2023_avg) %>%
  mutate(log_gdp=log(gdp))

#View(long_GDP_in_ESS_d)

# describe gdp average
psych::describe(GDP_in_ESS_d$gdp_2002_2023_avg)
psych::describe(GDP_in_ESS_d$log_gdp_2002_2023_avg)

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

  

#' 
#' ## Add GGGI-variable
#' 
#' * GGGI = Global Gender Gap Index
#' 
#' 
## -------------------------------------------------------------------------------------------------------
GGGI<-import("../data/GGGI.xlsx")

# check if all ESS countries are in the log_GGGI data
table(unique(diff_dat$cntry) %in% GGGI$ISO2)

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
  select(ISO2, cname, year, gggi, GGGI_2002_2023_avg)

# describe gggi average
psych::describe(GGGI_in_ESS_d$GGGI_2002_2023_avg)
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


#' 
#' ## Country-year-gender dataframe
#' 
#' Frame with one row for country-year-gender
#' 
## -------------------------------------------------------------------------------------------------------
diff_dat_cntry_year<-
  diff_dat %>% 
  group_by(cntry,year,gndr.c) %>%
  dplyr::summarise(n=n(),
                   n_wt=sum(pspwght),
                   FM.z.wt=weighted.mean(x=FM.z,w=pspwght),
                   FM.z=mean(FM.z),
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

#' 
#' 
#' # Descriptive statistics
#' 
#' ## Country-specific descriptives
#' 
## -------------------------------------------------------------------------------------------------------

# sample sizes from weights

cntry_n_frame<-
  diff_dat %>% group_by(cntry) %>%
  summarise('n ESS rounds' = mean(n_unique_essround),
            n=round(sum(pspwght),0))

# value-based male-typicality

cntry_FM_frame<-
  diff_dat %>% group_by(cntry,essround) %>%
  summarise('FM M' = weighted.mean(x=FM.z,w=pspwght),
            'FM SD' = sqrt(wtd.var(FM.z,pspwght))) %>%
  group_by(cntry) %>%
  summarise('FM M' = mean(x=`FM M`),
            'FM SD'= mean(x=`FM SD`))

cntry_FM_women_frame<-
  diff_dat %>%
  filter(gndr.c==-0.5) %>%
  group_by(cntry,essround) %>%
  summarise('FM M' = weighted.mean(x=FM.z,w=pspwght),
            'FM SD' = sqrt(wtd.var(FM.z,pspwght))) %>%
  group_by(cntry) %>%
  summarise('FM M Women' = mean(x=`FM M`),
            'FM SD Women'= mean(x=`FM SD`))

cntry_FM_men_frame<-
  diff_dat %>%
  filter(gndr.c==0.5) %>%
  group_by(cntry,essround) %>%
  summarise('FM M' = weighted.mean(x=FM.z,w=pspwght),
            'FM SD' = sqrt(wtd.var(FM.z,pspwght))) %>%
  group_by(cntry) %>%
  summarise('FM M Men' = mean(x=`FM M`),
            'FM SD Men'= mean(x=`FM SD`))


# link n and FM datasets

desc_frame<-
  left_join(
    x=cntry_n_frame,
    y=cntry_FM_frame,
    by="cntry"
  )

desc_frame<-
  left_join(
    x=desc_frame,
    y=cntry_FM_women_frame,
    by="cntry"
  )

desc_frame<-
  left_join(
    x=desc_frame,
    y=cntry_FM_men_frame,
    by="cntry"
  )

# Add country-specific differences
desc_frame$D<-desc_frame$`FM M Men`-desc_frame$`FM M Women`

desc_frame
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
    `FM M`, `FM SD`,
    `FM M Women`, `FM SD Women`,
    `FM M Men`, `FM SD Men`,
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

export(cntry_desc_tbl,"../results/cntry_desc_tbl.xlsx",overwrite=T)

#' 
#' ## Country-level correlations table
#' 
## -------------------------------------------------------------------------------------------------------
cor_frame<-
  desc_frame %>%
  select(
    VBMT=`FM M`,
    VBMT_Women=`FM M Women`,
    VBMT_Men=`FM M Men`,
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
  filename = "../results/CorTable1.doc",
  #table.number = NA,
  show.conf.interval = FALSE,
  show.sig.stars = F,
  landscape = F
)

#' ## Intraclass coefficient for gender
#' 
## -------------------------------------------------------------------------------------------------------
fit_gndr_cntry<-
  glmer(gndr.bin~(1|cntry),
        data=diff_dat,family=binomial(link="logit"))
summary(fit_gndr_cntry)
tau=VarCorr(fit_gndr_cntry)$cntry[1]
(ICC_gndr=tau/(tau+pi^2/3))


#' 
#' 
#' # Analysis
#' 
#' Following preregistration: https://osf.io/7cags?view_only=f3e97a78271e46bfafb9e20ac8d35bb1 
#' 
#' ## mod0: Random intercept model
#' 
## -------------------------------------------------------------------------------------------------------

mod0<-lmer(FM.z~(1|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod0)
getVC(mod0)
r2mlm(mod0,bargraph = F)

#' 
#' ## mod1: Gender fixed effect
#' 
## -------------------------------------------------------------------------------------------------------

mod1<-lmer(FM.z~gndr.c+(1|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod1)
getFE(mod1,round=3)
getVC(mod1)
r2mlm(mod1,bargraph = F)

#' 
#' ## mod2: Gender fixed and random effect
#' 
#' * Include random effect correlation by default
#' 
## -------------------------------------------------------------------------------------------------------

mod2<-lmer(FM.z~gndr.c+(gndr.c|cntry),data=diff_dat,REML=F,weights = pspwght,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2)
getFE(mod2,round=3)
getVC(mod2)
r2mlm(mod2,bargraph = F)
anova(mod1,mod2)
lvl2_var_cond_lvl1(mod2,lvl1.var = "gndr.c",lvl1.values = c(0.5,-0.5))

#' 
#' * Test for random effect correlation
#' 
## -------------------------------------------------------------------------------------------------------
mod2_norecov<-lmer(FM.z~gndr.c+(gndr.c||cntry),data=diff_dat,REML=F,weights = pspwght,
                   control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2_norecov)
getFE(mod2_norecov,round=3)
getVC(mod2_norecov)

anova(mod2_norecov,mod2)

#' 
#' ## ICC2 for country mean-levels of men and women
#' 
## -------------------------------------------------------------------------------------------------------
reliab_mod2<-reliability_dms(mod2,data = diff_dat,diff_var = "gndr.c",
                diff_var_values = c(-0.5,0.5),
                group_var = "cntry",var = "FM.z")
reliab_mod2

#' 
#' 
#' # Obtain VBMT with another training sample
#' 
## -------------------------------------------------------------------------------------------------------
# rerun the analysis
value_typ2<-
  D_regularized(data=fdat,mv.vars = value.vars,
                group.var = "gndr.bin",group.values = c(1,0),
                out = T,fold = T,fold.var = "cntry_time",size = 100,
                pcc = T,auc=T,pred.prob = T,append.data=T)

#' 
#' 
#' ## Correlations between different VBMT measures
#' 
## -------------------------------------------------------------------------------------------------------
# rerun the analysis with different training set
value_typ2<-
  D_regularized(data=fdat,mv.vars = value.vars,
                group.var = "gndr.bin",group.values = c(1,0),
                out = T,fold = T,fold.var = "cntry_time",size = 100,
                pcc = T,auc=T,pred.prob = T,append.data=T)

#' 
#' ### Country-year gender mean-levels
#' 
## -------------------------------------------------------------------------------------------------------
# women
cor.test(value_typ$D$m.0,value_typ2$D$m.0)
# men
cor.test(value_typ$D$m.1,value_typ2$D$m.1)

#' 
#' ## Country-year-individual variance accounted
#' 
#' Combine the separate frames
#' 
## -------------------------------------------------------------------------------------------------------

d1<-value_typ$preds
d2<-value_typ2$preds

d3<-left_join(x=d1,y=d2[,c("cntry","essround","idno","pred")],by=c("cntry","essround","idno"))


#' 
#' 
#' Separate by gender and calculate orthogonal predictors at country, country-round, and individual level
#' 
## -------------------------------------------------------------------------------------------------------

d3_men<-d3[d3$gndr.c==0.5,]
d3_women<-d3[d3$gndr.c==-0.5,]

d3_men_cntry<-
  d3_men %>%
  group_by(cntry)%>%
  summarise(y_cntry_mean=mean(pred.y,na.rm=T),
            x_cntry_mean=mean(pred.x,na.rm=T))

d3_men_cntry_essround<-
  d3_men %>%
  group_by(cntry,essround) %>%
  summarise(y_cntry_essround_mean=mean(pred.y,na.rm=T),
            x_cntry_essround_mean=mean(pred.x,na.rm=T))

d3_men<-
  left_join(
    x=d3_men,
    y=d3_men_cntry,
    by="cntry"
  )

d3_men<-
  left_join(
    x=d3_men,
    y=d3_men_cntry_essround,
    by=c("cntry","essround")
  )

head(d3_men)

# center individual scores around cntry-essround means
d3_men$y.crc<-d3_men$pred.y-d3_men$y_cntry_essround_mean
d3_men$x.crc<-d3_men$pred.x-d3_men$x_cntry_essround_mean

# center cntry-essround means around country means

d3_men$y_cntry_essround_mean.c<-d3_men$y_cntry_essround_mean-d3_men$y_cntry_mean
d3_men$x_cntry_essround_mean.c<-d3_men$x_cntry_essround_mean-d3_men$x_cntry_mean


d3_women_cntry<-
  d3_women %>%
  group_by(cntry)%>%
  summarise(y_cntry_mean=mean(pred.y,na.rm=T),
            x_cntry_mean=mean(pred.x,na.rm=T))

d3_women_cntry_essround<-
  d3_women %>%
  group_by(cntry,essround) %>%
  summarise(y_cntry_essround_mean=mean(pred.y,na.rm=T),
            x_cntry_essround_mean=mean(pred.x,na.rm=T))

d3_women<-
  left_join(
    x=d3_women,
    y=d3_women_cntry,
    by="cntry"
  )

d3_women<-
  left_join(
    x=d3_women,
    y=d3_women_cntry_essround,
    by=c("cntry","essround")
  )

head(d3_women)

# center individual scores around cntry-essround means
d3_women$y.crc<-d3_women$pred.y-d3_women$y_cntry_essround_mean
d3_women$x.crc<-d3_women$pred.x-d3_women$x_cntry_essround_mean

# center cntry-essround means around country means

d3_women$y_cntry_essround_mean.c<-d3_women$y_cntry_essround_mean-d3_women$y_cntry_mean
d3_women$x_cntry_essround_mean.c<-d3_women$x_cntry_essround_mean-d3_women$x_cntry_mean


#' 
#' 
#' ### Variance accounted for men
#' 
## -------------------------------------------------------------------------------------------------------
## empty model

mod0_men<-lmer(pred.x~(1|cntry/essround),data=d3_men,REML=F,weights = pspwght,
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

summary(mod0_men)
getVC(mod0_men,round=8)

## predictors model

mod1_men<-lmer(pred.x~y_cntry_mean+y_cntry_essround_mean.c+y.crc+(1|cntry/essround),data=d3_men,REML=F,weights = pspwght,
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

summary(mod1_men)
getVC(mod1_men,round=8)

# variance accounted

(data.frame(getVC(mod0_men,round=8))[,"vcov"]-data.frame(getVC(mod1_men,round=8))[,"vcov"])/
  data.frame(getVC(mod0_men,round=8))[,"vcov"]

#' 
#' 
#' ### Variance accounted for women
#' 
## -------------------------------------------------------------------------------------------------------
## empty model

mod0_women<-lmer(pred.x~(1|cntry/essround),data=d3_women,REML=F,weights = pspwght,
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

summary(mod0_women)
getVC(mod0_women,round=8)

## predictors model

mod1_women<-lmer(pred.x~y_cntry_mean+y_cntry_essround_mean.c+y.crc+(1|cntry/essround),data=d3_women,REML=F,weights = pspwght,
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

summary(mod1_women)
getVC(mod1_women,round=8)

# variance accounted

(data.frame(getVC(mod0_women,round=8))[,"vcov"]-data.frame(getVC(mod1_women,round=8))[,"vcov"])/
  data.frame(getVC(mod0_women,round=8))[,"vcov"]

#' 
#' 
#' # Session information
#' 
## -------------------------------------------------------------------------------------------------------
s<-sessionInfo()
print(s,locale=FALSE)

#' 
