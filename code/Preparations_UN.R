#' ---
#' title: "Data preparations for UN country indices"
#' output: 
#'   html_document: 
#'     toc: true
#'     keep_md: true
#' ---
#' 
## ----setup, include=FALSE-----------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

#' 
#' 
#' 
#' # Packages
#' 
## ----message=FALSE------------------------------------------------------------------------------
library(rio)
library(dplyr)
library(psych)

#' 
#' # Read Data
#' 
## -----------------------------------------------------------------------------------------------
UN_d<-read.csv("../data/HDR25_Composite_indices_complete_time_series.csv")

# read also ISO to obtain ISO2 codes
ISO<-read.csv2("../data/ISO.csv")

#' 
#' From Fors Connolly et al. 2020
#' 
#' The independent variable is Gender Equality Index (Gender Inequality Index reversed). The Gender Equality Index (GEI) was collected by the United Nations Development Programme (UNDP 2018). The index is a composite of health, empowerment, and labor market participation. GEI was only available for 2000 (1995 for the Czech Republic and Slovenia), 2005, and every year from 2010 and onward. We linearly imputed values for the missing years using the adjacent values. GEI changed relatively evenly within countries over years so the imputation procedure is not likely to misrepresent the true scores. All countries increased their gender equality over time, but to varying degrees. On a 0 to 1 scale (higher scores indicating higher gender equality), countries ranged from .43 (Turkey in 2004) and .96 (Denmark, Norway, Sweden, Switzerland in 2016). The GEI grand mean was .84 in 2002 and .92 in 2016. The within-country change in GEI between 2002 and 2016 was on average .06, ranging from .01 in Sweden to .15 in Poland.
#' 
#' ## GII/GEI (Gender inequality index, reverse coded as gender equality index)
#' 
## -----------------------------------------------------------------------------------------------
# variable names
GII_vars<-names(UN_d)[grepl("gii",names(UN_d)) & !grepl("gii_rank",names(UN_d))]

# select GII variables
GII_d<-UN_d %>%
  select(iso3,country,all_of(GII_vars))

# add ISO2 to GII_d
GII_d <- left_join(x=GII_d,
                   y=ISO[,c("ISO2","ISO3")],
                   by=c("iso3"="ISO3"))

# calculate GII mean across the time period 2002-2022 including the middle years

GII_vars_2002_2022<-
  paste0("gii_",2002:2022)
GII_vars_2002_2022

table(GII_vars_2002_2022 %in% names(GII_d))

GII_d$gii_2002_2022_avg<-
  rowMeans(GII_d[,GII_vars_2002_2022],na.rm=T)
  

describe(GII_d$gii_2002_2022_avg)

# reverse code

GII_d$gei_2002_2022_avg<-(1-GII_d$gii_2002_2022_avg)
describe(GII_d$gei_2002_2022_avg)
GII_d
rio::export(data.frame(GII_d),"../data/GII.csv",overwrite=T)

#' 
#' ## GDI (Gender development index)
#' 
## -----------------------------------------------------------------------------------------------
# variable names
GDI_vars<-names(UN_d)[grepl("gdi",names(UN_d)) & !grepl("gdi_group",names(UN_d))]

# select GDI variables
GDI_d<-UN_d %>%
  select(iso3,country,all_of(GDI_vars))

# add ISO2 to GDI_d
GDI_d <- left_join(x=GDI_d,
                   y=ISO[,c("ISO2","ISO3")],
                   by=c("iso3"="ISO3"))

# calculate GDI mean across the time period 2002-2022 including the middle years

GDI_vars_2002_2022<-
  paste0("gdi_",2002:2022)
GDI_vars_2002_2022

table(GDI_vars_2002_2022 %in% names(GDI_d))

GDI_d$gdi_2002_2022_avg<-
  rowMeans(GDI_d[,GDI_vars_2002_2022],na.rm=T)


describe(GDI_d$gdi_2002_2022_avg)

rio::export(data.frame(GDI_d),"../data/GDI.csv",overwrite=T)

#' 
