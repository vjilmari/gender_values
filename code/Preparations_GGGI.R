#' ---
#' title: "Data preparations for GGGI country indices"
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
library(tidyr)
library(psych)

#' 
#' # Read Data
#' 
## -----------------------------------------------------------------------------------------------
GGGI_d<-import("../data/qogdata_03_10_2025.xlsx")

# read also ISO to obtain ISO2 codes
ISO<-read.csv2("../data/ISO.csv")

#' 
#' ## GGGI (Global Gender Gap Index)
#' 
## -----------------------------------------------------------------------------------------------
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
names(GGGI_wide)

# calculate GGGI mean across the time period 2002-2022 including the middle years

GGGI_vars_2002_2022<-
  paste0("gggi_",2002:2022)
GGGI_vars_2002_2022


table(GGGI_vars_2002_2022 %in% names(GGGI_wide))
# which variables are not available
GGGI_vars_2002_2022[!(GGGI_vars_2002_2022 %in% names(GGGI_wide))]
# code only available ones
GGGI_vars_2002_2022_available<-
  GGGI_vars_2002_2022[(GGGI_vars_2002_2022 %in% names(GGGI_wide))]

GGGI_wide$GGGI_2002_2022_avg<-
  rowMeans(GGGI_wide[,GGGI_vars_2002_2022_available],na.rm=T)
  
describe(GGGI_wide$GGGI_2002_2022_avg)
rio::export(data.frame(GGGI_wide),"../data/GGGI.xlsx",overwrite=T)


#' 
#' 
#' 
