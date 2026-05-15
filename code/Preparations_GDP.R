#' ---
#' title: "Data preparations for country GDP"
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
## -------------------------------------------------------------------------------------------------------
library(rio)
library(dplyr)
library(psych)

#' 
#' 
#' # Data
#' 
## -------------------------------------------------------------------------------------------------------
GDP<-import("../data/GDP.xlsx")

# read also ISO to obtain ISO2 codes
ISO<-read.csv2("../data/ISO.csv")

#' 
#' ## GDP
#' 
## -------------------------------------------------------------------------------------------------------
# variable names
GDP_vars<-names(GDP)[grepl("YR",names(GDP))]

# select GDP variables
GDP_d<-GDP %>%
  select('Country Name','Country Code',all_of(GDP_vars))

# transform GDP variables to numeric
GDP_d[GDP_vars] <- sapply(GDP_d[GDP_vars], as.numeric)

# add ISO2 to GII_d
GDP_d <- left_join(x=GDP_d,
                   y=ISO[,c("ISO2","ISO3")],
                   by=c("Country Code"="ISO3"))



# calculate GDP mean across the time period 2002-2023 including the middle years

GDP_vars_2002_2023<-
  paste0(2002:2023," [YR",2002:2023,"]" )

table(GDP_vars_2002_2023 %in% names(GDP_d))

GDP_d$gdp_2002_2023_avg<-
  rowMeans(GDP_d[,GDP_vars_2002_2023],na.rm=T)
  

describe(GDP_d$gdp_2002_2023_avg)

# log transform

GDP_d$log_gdp_2002_2023_avg<-log(GDP_d$gdp_2002_2023_avg)
describe(GDP_d$log_gdp_2002_2023_avg)
GDP_d
rio::export(data.frame(GDP_d),"../data/GDP_processed.csv",overwrite=T)

#' 
