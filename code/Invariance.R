#' ---
#' title: "Invariance of prediction weights across countries and time"
#' output: 
#'   html_document: 
#'     toc: true
#'     keep_md: true
#' ---
#' 
## ----setup, include=FALSE-------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

#' 
#' # Pre-analysis
#' 
#' ## Packages
#' 
## -------------------------------------------------------------------------------------------------------
library(rio)
library(dplyr)
library(MetBrewer)
library(ggplot2)
library(finalfit)
library(ggflags)
library(lavaan)
library(semTools)
library(forcats)

#' 
#' ## Data
#' 
#' * Loads the preprocessed data file made within "Preparations" code
#' 
## -------------------------------------------------------------------------------------------------------
load("../data/fdat.rdata")

# get also ISO-codes for country names

ISO<-read.csv2("../data/ISO.csv")


#' 
#' ## Exclusions
#' 
#' * exclude participants with missing value variables or missing gender
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
  dplyr::filter(miss_values==0 & !is.na(gndr.bin))


#' 
#' ## Preparations
#' 
#' ### Scale gender to have SD/variance around 1
#' 
## -------------------------------------------------------------------------------------------------------
fdat$gndr.sd1<-2*fdat$gndr.c
sd(fdat$gndr.sd1)

#' 
#' # Analysis
#' 
#' ## Construct the model
#' 
#' Use a path model where a latent male-typicality (gender is the observed indicator) is predicted with the ten values
#' 
## -------------------------------------------------------------------------------------------------------

lavmod<-
  'MTYP=~gndr.sd1
  MTYP~con+tra+ben+uni+sdi+sti+hed+ach+pow+sec
  '


#' 
#' ## Fit the model as multigroup model across countries
#' 
## -------------------------------------------------------------------------------------------------------
t1<-Sys.time()

lavfit_cntry_reg<-
  sem(model=lavmod,data=fdat,
      group="cntry",group.equal="regressions")

t2<-Sys.time()
duration<-t2-t1
duration

#' 
#' ### Detect misspecifications at .10
#' 
#' * Regression coefficients from values that deviate by 0.10
#' 
## -------------------------------------------------------------------------------------------------------
t3<-Sys.time()
lavfit_cntry_reg_mi<-
  epcEquivFit(lavfit_cntry_reg,
             stdBeta = 0.1,
             cilevel = 0.95)
t4<-Sys.time()
duration2<-t4-t3
duration2

#' 
#' ### Examine misspecifications
#' 
## -------------------------------------------------------------------------------------------------------
# limit to regression parameters predicting MTYP
lavfit_cntry_reg_mi_reg<-
  lavfit_cntry_reg_mi[lavfit_cntry_reg_mi$op=="~" & 
                        lavfit_cntry_reg_mi$rhs=="MTYP",]

# decision based on modification index significance and power
table(lavfit_cntry_reg_mi_reg$decision.pow,useNA="always")
prop.table(table(lavfit_cntry_reg_mi_reg$decision.pow,useNA="always"))

# decision based on confidence intervals of expected parameter change (EPC)
table(lavfit_cntry_reg_mi_reg$decision.ci,useNA="always")
prop.table(table(lavfit_cntry_reg_mi_reg$decision.ci,useNA="always"))


#' 
#' * No clear misspecifications
#' 
#' ### Check country specific differences in the profile of coefficients
#' 
## -------------------------------------------------------------------------------------------------------
# get country labels from the model
cntry_labels<-
  data.frame(group=1:inspect(lavfit_cntry_reg,"ngroups"),
           cntry=inspect(lavfit_cntry_reg,"group.label"))

left_join(x=lavfit_cntry_reg_mi_reg,y=cntry_labels,by=c("group")) %>%
  group_by(cntry) %>%
  summarise(std.epc.M=mean(std.epc),
            std.epc.SD=sd(std.epc),
            std.epc.SQMSD=sqrt(mean(std.epc^2)),
            std.epc.MIN=min(std.epc),
            std.epc.MAX=max(std.epc))

# plot the EPCs
lavfit_cntry_reg_mi_reg<-
  left_join(x=lavfit_cntry_reg_mi_reg,y=cntry_labels,by=c("group"))


lavfit_cntry_reg_mi_reg <-
  lavfit_cntry_reg_mi_reg %>%
  mutate(Value=case_when(
    lhs=="con"~"Conformity",
    lhs=="tra"~"Tradition",
    lhs=="ben"~"Benevolence",
    lhs=="uni"~"Universalism",
    lhs=="sdi"~"Self-direction",
    lhs=="sti"~"Stimulation",
    lhs=="hed"~"Hedonism",
    lhs=="ach"~"Achievement",
    lhs=="pow"~"Power",
    lhs=="sec"~"Security"
  ))


ggplot(lavfit_cntry_reg_mi_reg,aes(x=std.epc,y=as.factor(cntry)))+
         geom_point(aes(size=abs(std.epc)))+
  facet_grid(.~Value)+
  geom_vline(xintercept=c(-0.10,0.10),linetype=2)+
  geom_vline(xintercept=c(0),linetype=1)+
  geom_errorbar(aes(xmin = lower.std.epc,xmax=upper.std.epc))+
  labs(x="Expected std. regression coefficient change if freed",
       y="Country")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
       


#' 
#' ### Detect misspecifications at .05
#' 
## -------------------------------------------------------------------------------------------------------
t5<-Sys.time()
lavfit_cntry_reg_mi_05<-
  epcEquivFit(lavfit_cntry_reg,
           stdBeta = 0.05,
           cilevel = 0.95)
t6<-Sys.time()
duration3<-t6-t5
duration3


#' 
#' ### Examine misspecifications
#' 
## -------------------------------------------------------------------------------------------------------
# limit to regression parameters predicting MTYP
lavfit_cntry_reg_mi_05_reg<-
  lavfit_cntry_reg_mi_05[lavfit_cntry_reg_mi_05$op=="~" & 
                        lavfit_cntry_reg_mi_05$rhs=="MTYP",]

# decision based on modification index significance and power
table(lavfit_cntry_reg_mi_05_reg$decision.pow,useNA="always")
proportions(table(lavfit_cntry_reg_mi_05_reg$decision.pow,useNA="always"))

# decision based on confidence intervals of expected parameter change (EPC)
table(lavfit_cntry_reg_mi_05_reg$decision.ci,useNA="always")
prop.table(table(lavfit_cntry_reg_mi_05_reg$decision.ci,useNA="always"))

# combined criteria
table(lavfit_cntry_reg_mi_05_reg$decision.pow,
      lavfit_cntry_reg_mi_05_reg$decision.ci)


lavfit_cntry_reg_mi_05_reg<-
  lavfit_cntry_reg_mi_05_reg %>%
  mutate(misspec=case_when(
    decision.pow=="EPC:M" | decision.pow=="M" | decision.ci=="M"~1,
    TRUE~0
  ))

table(lavfit_cntry_reg_mi_05_reg$misspec,useNA="always")
prop.table(table(lavfit_cntry_reg_mi_05_reg$misspec,useNA="always"))

# check if the are in specific countries
lavfit_cntry_reg_mi_05_reg<-
  left_join(x=lavfit_cntry_reg_mi_05_reg,y=cntry_labels,by=c("group"))

lavfit_cntry_reg_mi_05_reg %>% 
  group_by(cntry) %>%
  summarise(n_misspec=sum(misspec))

table(lavfit_cntry_reg_mi_05_reg$cntry,
      lavfit_cntry_reg_mi_05_reg$misspec)

# save results

export(lavfit_cntry_reg_mi_05_reg %>%
  group_by(cntry) %>%
  summarise(n_misspec=sum(misspec)) %>%
  arrange(desc(cntry)),"../results/n_misspec.xlsx",overwrite=T)


#' 
#' ### Plot country results
#' 
## -------------------------------------------------------------------------------------------------------

lavfit_cntry_reg_mi_05_reg <-
  lavfit_cntry_reg_mi_05_reg %>%
  mutate(Value=case_when(
    lhs=="con"~"Conformity",
    lhs=="tra"~"Tradition",
    lhs=="ben"~"Benevolence",
    lhs=="uni"~"Universalism",
    lhs=="sdi"~"Self-direction",
    lhs=="sti"~"Stimulation",
    lhs=="hed"~"Hedonism",
    lhs=="ach"~"Achievement",
    lhs=="pow"~"Power",
    lhs=="sec"~"Security"
  ))


# add correctly ordered country
lavfit_cntry_reg_mi_05_reg$cntry <- 
  fct_rev(fct_reorder(as.factor(lavfit_cntry_reg_mi_05_reg$cntry),
                      lavfit_cntry_reg_mi_05_reg$cntry))

# add average deviation for each country

cntry_avg_epc<-
  lavfit_cntry_reg_mi_05_reg %>%
  group_by(cntry) %>%
  summarise(std.epc.M=mean(std.epc),
            std.epc.SD=sd(std.epc),
            std.epc.MAD=mean(abs(std.epc)),
            std.epc.MIN=min(std.epc),
            std.epc.MAX=max(std.epc)) %>%
  arrange(desc(cntry))

print(cntry_avg_epc,n=100)


lavfit_cntry_reg_mi_05_reg<-
  left_join(x=lavfit_cntry_reg_mi_05_reg,y=cntry_avg_epc,by=c("cntry"))

# country label with MAD
lavfit_cntry_reg_mi_05_reg$cntry_lab<-
  paste0(lavfit_cntry_reg_mi_05_reg$cntry,
         " (",round_tidy(lavfit_cntry_reg_mi_05_reg$std.epc.MAD,2),")")


# need to reorder cntry_lab as well
lavfit_cntry_reg_mi_05_reg$cntry_lab <- 
  fct_rev(fct_reorder(as.factor(lavfit_cntry_reg_mi_05_reg$cntry_lab),
                      lavfit_cntry_reg_mi_05_reg$cntry_lab))

misspec_plot<-
  ggplot(lavfit_cntry_reg_mi_05_reg,aes(x=std.epc,y=cntry_lab,
                                      color=as.factor(misspec)))+
  geom_point(aes(size=abs(std.epc)))+
  facet_grid(.~Value)+
  geom_vline(xintercept=c(-0.05,0.05),linetype=2)+
  geom_vline(xintercept=c(0),linetype=1)+
  geom_errorbar(aes(xmin = lower.std.epc,xmax=upper.std.epc))+
  labs(x="Expected std. regression coefficient change if parameter freed (EPC)",
       y="Country (EPC Mean Absolute Deviation)",
       color=NULL,  
       size="Absolute EPC")+
  scale_color_manual(values=c("0"=met.brewer("Demuth")[8], "1"=met.brewer("Demuth")[2]),
                     labels=c("0"="Not misspecified", "1"="Misspecified"))+
  geom_flag(aes(country=tolower(cntry)),size=2.5)+
  scale_x_continuous(breaks = c(-0.05, 0.05))
  
misspec_plot

ggsave(misspec_plot,filename="misspec_plot_cntry.png",
       path="../results",
       device = "png",width = 29.7,height=18,dpi = 600,units = "cm")                     


#' 
#' ## Fit the model as multigroup model across ess rounds (time)
#' 
## -------------------------------------------------------------------------------------------------------
t7<-Sys.time()

lavfit_essround_reg<-
  sem(model=lavmod,data=fdat,
      group="essround",group.equal="regressions")

t8<-Sys.time()
duration<-t8-t7
duration

#' 
#' ### Detect misspecifications at .10
#' 
#' * Regression coefficients from values that deviate by 0.10
#' 
## -------------------------------------------------------------------------------------------------------
t9<-Sys.time()
lavfit_essround_reg_mi<-
  epcEquivFit(lavfit_essround_reg,
           stdBeta = 0.1,
           cilevel = 0.95)
t10<-Sys.time()
duration2<-t10-t9
duration2


#' 
#' ### Examine misspecifications
#' 
## -------------------------------------------------------------------------------------------------------
# limit to regression parameters predicting MTYP
lavfit_essround_reg_mi_reg<-
  lavfit_essround_reg_mi[lavfit_essround_reg_mi$op=="~" & 
                        lavfit_essround_reg_mi$rhs=="MTYP",]

# decision based on modification index significance and power
table(lavfit_essround_reg_mi_reg$decision.pow,useNA="always")
prop.table(table(lavfit_essround_reg_mi_reg$decision.pow,useNA="always"))

# decision based on confidence intervals of expected parameter change (EPC)
table(lavfit_essround_reg_mi_reg$decision.ci,useNA="always")
prop.table(table(lavfit_essround_reg_mi_reg$decision.ci,useNA="always"))


#' 
#' * No misspecifications
#' 
#' ### Check country specific differences in the profile of coefficients
#' 
## -------------------------------------------------------------------------------------------------------
# get country labels from the model
essround_labels<-
  data.frame(group=1:inspect(lavfit_essround_reg,"ngroups"),
           essround=inspect(lavfit_essround_reg,"group.label"))


left_join(x=lavfit_essround_reg_mi_reg,y=essround_labels,by=c("group")) %>%
  group_by(essround) %>%
  summarise(std.epc.M=mean(std.epc),
            std.epc.SD=sd(std.epc),
            std.epc.SQMSD=sqrt(mean(std.epc^2)),
            std.epc.MIN=min(std.epc),
            std.epc.MAX=max(std.epc))

# plot the EPCs
lavfit_essround_reg_mi_reg<-
  left_join(x=lavfit_essround_reg_mi_reg,y=essround_labels,by=c("group"))


lavfit_essround_reg_mi_reg <-
  lavfit_essround_reg_mi_reg %>%
  mutate(Value=case_when(
    lhs=="con"~"Conformity",
    lhs=="tra"~"Tradition",
    lhs=="ben"~"Benevolence",
    lhs=="uni"~"Universalism",
    lhs=="sdi"~"Self-direction",
    lhs=="sti"~"Stimulation",
    lhs=="hed"~"Hedonism",
    lhs=="ach"~"Achievement",
    lhs=="pow"~"Power",
    lhs=="sec"~"Security"
  ))

lavfit_essround_reg_mi_reg$essround <- factor(
  lavfit_essround_reg_mi_reg$essround,
  levels = sort(as.numeric(unique(lavfit_essround_reg_mi_reg$essround)))
)

ggplot(lavfit_essround_reg_mi_reg, aes(x = std.epc, y = essround)) +
  geom_point(aes(size = abs(std.epc))) +
  facet_grid(. ~ Value) +
  geom_vline(xintercept = c(-0.10, 0.10), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 1) +
  geom_errorbar(aes(xmin = lower.std.epc, xmax = upper.std.epc)) +
  labs(
    x = "Expected std. regression coefficient change if freed",
    y = "Essround"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



#' 
#' 
#' ### Detect misspecifications at .05
#' 
## -------------------------------------------------------------------------------------------------------
t5<-Sys.time()
lavfit_essround_reg_mi_05<-
  epcEquivFit(lavfit_essround_reg,
           stdBeta = 0.05,
           cilevel = 0.95)
t6<-Sys.time()
duration3<-t6-t5
duration3


#' 
#' ### Examine misspecifications
#' 
## -------------------------------------------------------------------------------------------------------
# limit to regression parameters predicting MTYP
lavfit_essround_reg_mi_05_reg<-
  lavfit_essround_reg_mi_05[lavfit_essround_reg_mi_05$op=="~" & 
                        lavfit_essround_reg_mi_05$rhs=="MTYP",]

# decision based on modification index significance and power
table(lavfit_essround_reg_mi_05_reg$decision.pow,useNA="always")
proportions(table(lavfit_essround_reg_mi_05_reg$decision.pow,useNA="always"))


# decision based on confidence intervals of expected parameter change (EPC)
table(lavfit_essround_reg_mi_05_reg$decision.ci,useNA="always")
prop.table(table(lavfit_essround_reg_mi_05_reg$decision.ci,useNA="always"))

# combined criteria
table(lavfit_essround_reg_mi_05_reg$decision.pow,
      lavfit_essround_reg_mi_05_reg$decision.ci)


table(lavfit_essround_reg_mi_05_reg$decision.ci)


lavfit_essround_reg_mi_05_reg<-
  lavfit_essround_reg_mi_05_reg %>%
  mutate(misspec=case_when(
    decision.pow=="EPC:M" | decision.pow=="M" | decision.ci=="M"~1,
    TRUE~0
  ))

table(lavfit_essround_reg_mi_05_reg$misspec,useNA="always")


# check if these are in specific essrounds
lavfit_essround_reg_mi_05_reg<-
  left_join(x=lavfit_essround_reg_mi_05_reg,y=essround_labels,by=c("group"))

lavfit_essround_reg_mi_05_reg %>% 
  group_by(essround) %>%
  summarise(n_misspec=sum(misspec))

table(lavfit_essround_reg_mi_05_reg$essround,
      lavfit_essround_reg_mi_05_reg$misspec)




#' 
#' ### Plot essround results
#' 
## -------------------------------------------------------------------------------------------------------

lavfit_essround_reg_mi_05_reg <-
  lavfit_essround_reg_mi_05_reg %>%
  mutate(Value=case_when(
    lhs=="con"~"Conformity",
    lhs=="tra"~"Tradition",
    lhs=="ben"~"Benevolence",
    lhs=="uni"~"Universalism",
    lhs=="sdi"~"Self-direction",
    lhs=="sti"~"Stimulation",
    lhs=="hed"~"Hedonism",
    lhs=="ach"~"Achievement",
    lhs=="pow"~"Power",
    lhs=="sec"~"Security"
  ))


# add average deviation for each essround

essround_avg_epc<-
  lavfit_essround_reg_mi_05_reg %>%
  group_by(essround) %>%
  summarise(std.epc.M=mean(std.epc),
            std.epc.SD=sd(std.epc),
            std.epc.MAD=mean(abs(std.epc)),
            std.epc.MIN=min(std.epc),
            std.epc.MAX=max(std.epc)) %>%
  arrange(desc(essround))

print(essround_avg_epc,n=100)

lavfit_essround_reg_mi_05_reg<-
  left_join(x=lavfit_essround_reg_mi_05_reg,
            y=essround_avg_epc,by=c("essround"))


# make a factor of essround variable
lavfit_essround_reg_mi_05_reg$essround<-
  factor(lavfit_essround_reg_mi_05_reg$essround)

lavfit_essround_reg_mi_05_reg$essround<-
  factor(lavfit_essround_reg_mi_05_reg$essround,
         levels = sort(as.numeric(levels(lavfit_essround_reg_mi_05_reg$essround))))



# essround label with MAD
lavfit_essround_reg_mi_05_reg$essround_lab<-
  paste0(lavfit_essround_reg_mi_05_reg$essround,
         " (",round_tidy(lavfit_essround_reg_mi_05_reg$std.epc.MAD,2),")")

# need to reorder essround_lab as well
lavfit_essround_reg_mi_05_reg$essround_lab <- 
  fct_rev(fct_reorder(as.factor(lavfit_essround_reg_mi_05_reg$essround_lab),
                      as.numeric(lavfit_essround_reg_mi_05_reg$essround)))


misspec_plot_essround<-
  ggplot(lavfit_essround_reg_mi_05_reg,aes(x=std.epc,y=essround_lab,
                                      color=as.factor(misspec)))+
  geom_point(aes(size=abs(std.epc)))+
  facet_grid(.~Value)+
  geom_vline(xintercept=c(-0.05,0.05),linetype=2)+
  geom_vline(xintercept=c(0),linetype=1)+
  geom_errorbar(aes(xmin = lower.std.epc,xmax=upper.std.epc))+
  labs(x="Expected std. regression coefficient change if parameter freed (EPC)",
       y="ESS round (EPC Mean Absolute Deviation)",
       color=NULL,  
       size="Absolute EPC")+
  scale_color_manual(values=c("0"=met.brewer("Demuth")[8], "1"=met.brewer("Demuth")[2]),
                     labels=c("0"="Not misspecified", "1"="Misspecified"))+
  #geom_flag(aes(country=tolower(essround)),size=2.5)+
  scale_x_continuous(breaks = c(-0.05, 0.05),limits = c(-0.10,0.10))
  

misspec_plot_essround

ggsave(misspec_plot_essround,filename="misspec_plot_essround.png",
       path="../results",
       device = "png",width = 29.7,height=18,dpi = 600,units = "cm")                     


#' 
#' # Session information
#' 
## -------------------------------------------------------------------------------------------------------
s<-sessionInfo()
print(s,locale=F)

