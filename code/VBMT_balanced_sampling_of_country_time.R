#' ---
#' title: "Examinations of VBMT from equal number of country-time samples"
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
#' 
#' 
#' # Obtain VBMT with another training sample
#' 
#' Limit this training sample so that there is equal number of people from all participating countries. Choose 200 men and 200 women from each country. In this case, the training data is also a lot smaller.
#' 
#' 
## -------------------------------------------------------------------------------------------------------
# rerun the analysis
value_typ2<-
  D_regularized(data=fdat,mv.vars = value.vars,
                group.var = "gndr.bin",group.values = c(1,0),
                out = T,fold = T,fold.var = "cntry",size = 200,
                pcc = T,auc=T,pred.prob = T,append.data=T)

#' 
#' 
#' ## Country-year-individual variance accounted
#' 
#' Combine the separate frames
#' 
## -------------------------------------------------------------------------------------------------------

d1<-value_typ$preds
d2<-value_typ2$preds

d3<-left_join(x=d1,y=d2[,c("cntry","essround","idno","pred")],
              by=c("cntry","essround","idno"))


#' 
#' 
#' Separate by gender and calculate orthogonal predictors at country, country-round, and individual level by centering.
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

# center individual scores around cntry-essround means
d3_men$y.crc<-d3_men$pred.y-d3_men$y_cntry_essround_mean
d3_men$x.crc<-d3_men$pred.x-d3_men$x_cntry_essround_mean

# center cntry-essround means around country means

d3_men$y_cntry_essround_mean.c<-
  d3_men$y_cntry_essround_mean-d3_men$y_cntry_mean
d3_men$x_cntry_essround_mean.c<-
  d3_men$x_cntry_essround_mean-d3_men$x_cntry_mean


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

mod0_men<-lmer(pred.x~(1|cntry/essround),data=d3_men,REML=F,
               weights = pspwght,subset=!is.na(d3_men$y.crc),
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

summary(mod0_men)
getVC(mod0_men,round=8)

## predictors model

mod1_men<-lmer(pred.x~y_cntry_mean+y_cntry_essround_mean.c+y.crc+(1|cntry/essround),
               data=d3_men,REML=F,weights = pspwght,subset=!is.na(d3_men$y.crc),
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
                 subset=!is.na(d3_women$y.crc),
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

summary(mod0_women)
getVC(mod0_women,round=8)

## predictors model

mod1_women<-lmer(pred.x~y_cntry_mean+y_cntry_essround_mean.c+y.crc+(1|cntry/essround),
                 data=d3_women,REML=F,weights = pspwght,subset=!is.na(d3_women$y.crc),
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
