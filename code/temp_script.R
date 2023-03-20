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

fdat<-import("data/fdat.xlsx")

# exclude participants with missing values or gender
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

# set seed number 

set.seed(13032023)
table(fdat$cntry_time)
value_typ<-
  D_regularized(data=fdat,mv.vars = value.vars,
                group.var = "gndr.bin",group.values = c(1,0),
                out = T,fold = T,fold.var = "cntry_time",size = 100,
                pcc = T,auc=T,pred.prob = T,append.data=T)
round(value_typ$D,2)
round(coefficients(value_typ$cv.mod,s = "lambda.min"),2)
plot(value_typ$cv.mod)

# export summary table
sum_tab<-value_typ$D
range(sum_tab$D)
range(sum_tab$d.sd.total)
mean(sum_tab$d.sd.total)
mean(sum_tab$pcc.total)
mean(sum_tab$pcc.1)
mean(sum_tab$pcc.0)
sum_tab[sum_tab$d.sd.total==min(sum_tab$d.sd.total),]
sum_tab[sum_tab$d.sd.total==max(sum_tab$d.sd.total),]

cor(sum_tab$m.1,sum_tab$m.0)
cor(sum_tab$sd.1,sum_tab$sd.0)

# export testing data
test_dat<-
  value_typ$preds
str(test_dat)

test_dat$age.c<-as.numeric(test_dat$age.c)
test_dat$rlgdgr.c<-as.numeric(test_dat$rlgdgr.c)
test_dat$eduyrs.c<-as.numeric(test_dat$eduyrs.c)
#export(value_typ$preds,"data/testdat.xlsx",overwrite=T)

# data exclusions

test_dat<-test_dat %>%
  dplyr::filter(age_included==1 & waves>1 &
                  same_gndr_partner==0 &
                  #essround < 9 &
                  !is.na(childless3) &
                  !is.na(married.c) &
                  !is.na(eduyrs.c) &
                  !is.na(rlgdgr.c) &
                  !is.na(age.c))

table(test_dat$waves)

# scale FM with pooled SD 

hist(test_dat$pred)
FM_pooled_sd<-value_typ$D[1,"pooled.sd.total"]
test_dat$FM.z<-test_dat$pred/FM_pooled_sd
hist(test_dat$FM.z)

# center around country mean-levels
cntry_FM<-test_dat %>%
  group_by(cntry) %>%
  summarise(FM.z.cm=mean(FM.z,na.rm=T))

test_dat<-left_join(
  x=test_dat,
  y=cntry_FM,
  by="cntry"
)  


test_dat$FM.z.cmc<-test_dat$FM.z-test_dat$FM.z.cm

# Analysis

## null model

fit0<-glmer(childless3~(1|cntry),data=test_dat,
            weights = pspwght,family=binomial(link="logit"),
            control=glmerControl(optimizer="bobyqa"))
summary(fit0)

## add gender fixed effect
fit1<-glmer(childless3~gndr.c+(1|cntry),data=test_dat,
            weights = pspwght,family=binomial(link="logit"),
            control=glmerControl(optimizer="bobyqa"))
summary(fit1)

## add FM fixed
fit2<-glmer(childless3~gndr.c+FM.z.cmc+(1|cntry),data=test_dat,
            weights = pspwght,family=binomial(link="logit"),
            control=glmerControl(optimizer="bobyqa"))
summary(fit2)
getFE_glmer(fit2)
getVC(fit2)

## add FM random
fit3<-glmer(childless3~gndr.c+FM.z.cmc+(FM.z.cmc|cntry),data=test_dat,
            weights = pspwght,family=binomial(link="logit"),
            control=glmerControl(optimizer="bobyqa"))
summary(fit3)
anova(fit2,fit3)

# add gender random
fit4<-glmer(childless3~gndr.c+FM.z.cmc+(gndr.c+FM.z.cmc|cntry),data=test_dat,
            weights = pspwght,family=binomial(link="logit"),
            control=glmerControl(optimizer="bobyqa"))
summary(fit4)
anova(fit3,fit4)

# add FM by gender interaction fixed

fit5<-glmer(childless3~gndr.c+FM.z.cmc+gndr.c:FM.z.cmc+
              (gndr.c+FM.z.cmc|cntry),data=test_dat,
            weights = pspwght,family=binomial(link="logit"),
            control=glmerControl(optimizer="bobyqa"))
summary(fit5)
anova(fit4,fit5)
getFE_glmer(fit5)
getVC(fit5)

FM.slopes.fit5<-
  emtrends(fit5,var="FM.z.cmc",specs="gndr.c",
           at=list(gndr.c=c(0.5,-0.5)),infer=c(T,T))
FM.slopes.fit5

# add FM by gender interaction random

fit6<-glmer(childless3~gndr.c+FM.z.cmc+gndr.c:FM.z.cmc+
              (gndr.c+FM.z.cmc+gndr.c:FM.z.cmc|cntry),data=test_dat,
            weights = pspwght,family=binomial(link="logit"),
            control=glmerControl(optimizer="bobyqa",
                                 optCtrl=list(maxfun=2e7)))
summary(fit6)
anova(fit5,fit6)

FM.slopes.fit6<-
  emtrends(fit6,var="FM.z.cmc",specs="gndr.c",
           at=list(gndr.c=c(0.5,-0.5)),infer=c(T,T))
FM.slopes.fit6



# test with covariates

fit6_covariates<-glmer(childless3~gndr.c+FM.z.cmc+gndr.c:FM.z.cmc+
                         married.c+eduyrs.c+rlgdgr.c+age.c+essround.c+
                         (gndr.c+FM.z.cmc+gndr.c:FM.z.cmc|cntry),
                       data=test_dat,
                       weights = pspwght,family=binomial(link="logit"),
                       control=glmerControl(optimizer="bobyqa",
                                            optCtrl=list(maxfun=2e7)))
summary(fit6_covariates)
getFE_glmer(fit6_covariates)
getVC(fit6_covariates)

FM.slopes.fit6_covariates<-
  emtrends(fit6_covariates,var="FM.z.cmc",specs="gndr.c",
           at=list(gndr.c=c(0.5,-0.5)),infer=c(T,T))
FM.slopes.fit6_covariates
exp(c(0.0891,0.0605,0.118))
exp(c(0.2580,0.2244,0.292))
pairs(FM.slopes.fit6_covariates,adjust="none",infer=c(T,T))
exp(c(-0.169,-0.207,-0.131))
# aggregate so that the question is about gender-typicality and childlessness

(0.0891+(-0.2580))/2

mlist <- list(
  gender_typicality = c(0.5, -0.5))

contrast(FM.slopes.fit6_covariates,method=mlist,adjust="none",infer=c(T,T))
exp(c(-0.0845,-0.103,-0.0656))

# add time component to the model

fit7<-glmer(childless3~gndr.c+FM.z.cmc+gndr.c:FM.z.cmc+
                         married.c+eduyrs.c+rlgdgr.c+age.c+
              essround.c+essround.c:gndr.c+
              essround.c:FM.z.cmc+essround.c:FM.z.cmc:gndr.c+
                         (gndr.c+FM.z.cmc+gndr.c:FM.z.cmc|cntry),
                       data=test_dat,
                       weights = pspwght,family=binomial(link="logit"),
                       control=glmerControl(optimizer="bobyqa",
                                            optCtrl=list(maxfun=2e7)))
summary(fit7)
getFE_glmer(fit7)
anova(fit6_covariates,fit7)

# general time trend
ref_grid(fit7)

emtrends(fit7,var="essround.c",specs=c("essround.c"),
         at=list(essround.c=0),infer=c(T,T))
exp(c(0.0288,0.0231,0.0345))
exp(c(10*0.0288,10*0.0231,10*0.0345))

# change of typicality-childlessness link over time

time_main_fit7<-
  emtrends(fit7,var="FM.z.cmc",specs=c("essround.c"),
         at=list(essround.c=seq(from=-4.5,to=4.5,by=1)),infer=c(T,T))
time_main_fit7
round(exp(c(0.182,0.147,0.218)),2)
round(exp(c(0.165,0.129 ,0.201)),2)

contrast(time_main_fit7,
         method=list(eff=c(-1,1,0,0,0,0,0,0,0,0)),
         adjust="none",infer=c(T,T))
exp(c(-0.00195,-0.00773,0.00382))

contrast(time_main_fit7,
         method=list(eff=c(-1,0,0,0,0,0,0,0,0,1)),
         adjust="none",infer=c(T,T))
exp(c(-0.0176,-0.0695,0.0344))


FM.slopes.fit7<-
  emtrends(fit7,var="FM.z.cmc",specs=c("essround.c","gndr.c"),
           at=list(gndr.c=c(0.5,-0.5),
                   essround.c=seq(from=-4.5,to=4.5,by=1)),infer=c(T,T))
FM.slopes.fit7
exp(c(-0.0176,-0.0695,0.0344))


FM.contrast.women.fit7<-
  emtrends(fit7,var="FM.z.cmc",specs=c("essround.c","gndr.c"),
                                 at=list(gndr.c=c(-0.5),
                                         essround.c=seq(from=-4.5,to=4.5,by=1)),
           infer=c(T,T))

FM.contrast.women.fit7

exp(c(0.0932,0.0476,0.139))
exp(c(0.0840,0.0367,0.131))

contrast(FM.contrast.women.fit7,
         method=list(eff=c(-1,1,0,0,0,0,0,0,0,0)),
         adjust="none",infer=c(T,T))
exp(c(-0.00288,-0.011,0.0052))

contrast(FM.contrast.women.fit7,
         method=list(eff=c(-1,0,0,0,0,0,0,0,0,1)),
         adjust="none",infer=c(T,T))
exp(c(-0.026,-0.0987,0.0468))



FM.contrast.men.fit7<-
  emtrends(fit7,var="FM.z.cmc",specs=c("essround.c","gndr.c"),
           at=list(gndr.c=c(0.5),
                   essround.c=seq(from=-4.5,to=4.5,by=1)),infer=c(T,T))
FM.contrast.men.fit7

exp(c(0.0932,0.0476,0.139))
exp(c(0.0840,0.0367,0.131))

contrast(FM.contrast.men.fit7,
         method=list(eff=c(-1,1,0,0,0,0,0,0,0,0)),
         adjust="none",infer=c(T,T))
exp(c(-0.00102,-0.00914,0.00711))

contrast(FM.contrast.men.fit7,
         method=list(eff=c(-1,0,0,0,0,0,0,0,0,1)),
         adjust="none",infer=c(T,T))
exp(c(-0.00916,-0.0823,0.064))


FM.contrast.time.fit7<-
  emtrends(fit7,var="FM.z.cmc",specs=c("essround.c"),
           at=list(
                   essround.c=c(0.5,1.5),infer=c(T,T)))
FM.contrast.time.fit7
pairs(FM.contrast.time.fit7,adjust="none",reverse=T)

# plotting
range(test_dat$FM.z.cmc)

p<-
  emmip(fit6_covariates, gndr.c ~ FM.z.cmc,
        at=list(gndr.c = c(-0.5,0.5),
                FM.z.cmc=seq(from=-3,to=3,by=0.01)),
        plotit=F,CIs=TRUE,type="response")
head(p)

p$gender<-p$tvar

levels(p$gender)<-c("Women","Men")

min.women<-
  min(test_dat[test_dat$gndr.c==(-0.5),
               "FM.z.cmc"])
max.women<-
  max(test_dat[test_dat$gndr.c==(-0.5),
               "FM.z.cmc"])

min.men<-
  min(test_dat[test_dat$gndr.c==(0.5),
               "FM.z.cmc"])
max.men<-
  max(test_dat[test_dat$gndr.c==(0.5),
               "FM.z.cmc"])

p$filter.women<-
  ifelse(p$gender=="Women" &
           (p$FM.z.cmc<min.women | 
              p$FM.z.cmc>max.women),0,1)
table(p$filter.women)

p$filter.men<-
  ifelse(p$gender=="Men" &
           (p$FM.z.cmc<min.men | 
              p$FM.z.cmc>max.men),0,1)
table(p$filter.men)

p.ex<-p[p$filter.men!=0 & p$filter.women!=0,]

met.brewer("Cassatt2")
met.brewer("Cassatt2")[c(2,9)]

p1<-ggplot(p.ex,aes(y=yvar,x=xvar,color=gender))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.5)+
  xlab("Male-typicality of values")+
  ylab("P(Childlessness)")+
  scale_color_manual(values=met.brewer("Cassatt2")[c(2,9)])+
  theme(legend.position = "top",
        legend.title=element_blank(),
        text=element_text(size=16,  family="sans"),
        panel.background = element_rect(fill = "white",
                                        #colour = "black",
                                        #size = 0.5, linetype = "solid"
        ),
        panel.grid.major.y = element_line(linewidth = 0.5, linetype = 2,
                                          colour = "black"))

png(filename = 
      "results/adjusted_main_effects.png",
    units = "cm",
    width = 21.0,height=29.7*(3/4),res = 300)
p1
dev.off()

# flip the x-axis to gender-typicality

head(p)
p$gender_typicality<-
  ifelse(p$gender=="Women",-1*p$xvar,p$xvar)
c(min.women,max.women)

p$filter2.women<-
  ifelse(p$gender=="Women" &
           (p$gender_typicality>(-1)*min.women | 
              p$gender_typicality<(-1)*max.women),0,1)
table(p$filter2.women)

p.ex2<-p[p$filter.men!=0 & p$filter2.women!=0,]


p2<-ggplot(p.ex2,aes(y=yvar,x=gender_typicality,color=gender))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LCL, ymax=UCL),alpha=0.5)+
  xlab("Gender-typicality of values")+
  ylab("P(Childlessness)")+
  scale_color_manual(values=met.brewer("Cassatt2")[c(2,9)])+
  theme(legend.position = "top",
        legend.title=element_blank(),
        text=element_text(size=16,  family="sans"),
        panel.background = element_rect(fill = "white",
                                        #colour = "black",
                                        #size = 0.5, linetype = "solid"
        ),
        panel.grid.major.y = element_line(linewidth = 0.5, linetype = 2,
                                          colour = "black"))

png(filename = 
      "results/adjusted_main_effects_2.png",
    units = "cm",
    width = 21.0,height=29.7*(3/4),res = 300)
p2
dev.off()


table(test_dat$childless3,test_dat$gndr.bin)
table(fdat$childless3,fdat$gndr.bin)

# obtain random effects for countries
coefficients(fit6_covariates)$cntry

# check if the same pattern persists without covariates
fit_gndr<-glmer(childless3~gndr.c+(gndr.c|cntry),data=test_dat,
                weights = pspwght,family=binomial(link="logit"),
                control=glmerControl(optimizer="bobyqa"))
summary(fit_gndr)
getFE_glmer(fit_gndr)
getVC(fit_gndr)

round(exp(coefficients(fit_gndr)$cntry["gndr.c"]),2)

# check childlessness associations for each country

fit6_coefs<-coefficients(fit6_covariates)$cntry
fit6_coefs$male_typicality_OR=exp(fit6_coefs$FM.z.cmc)
fit6_coefs$male_typicality_OR_women=exp(fit6_coefs$FM.z.cmc+(-0.5)*fit6_coefs$'gndr.c:FM.z.cmc')
fit6_coefs$male_typicality_OR_men=exp(fit6_coefs$FM.z.cmc+(0.5)*fit6_coefs$'gndr.c:FM.z.cmc')
fit6_coefs

cntry_OR<-data.frame(
  cntry=rownames(fit6_coefs),
  male_typicality_OR=round_tidy(fit6_coefs$male_typicality_OR,2),
  male_typicality_OR_women=round_tidy(fit6_coefs$male_typicality_OR_women,2),
  male_typicality_OR_men=round_tidy(fit6_coefs$male_typicality_OR_men,2))
cntry_OR
fit6_coefs

cntry.plot<-
  ggplot(fit6_coefs,aes(y=reorder(rownames(fit6_coefs),male_typicality_OR),x=male_typicality_OR))+
  geom_point(color="black",size=8)+
  geom_flag(aes(country=tolower(rownames(fit6_coefs))))+
  xlab("Odds ratio of Male-typicality on Childlessness")+
  ylab("Country")+
  geom_point(inherit.aes=F,color=met.brewer("Cassatt2")[3],size=8,
             aes(y=rownames(fit6_coefs),x=male_typicality_OR_women))+
  geom_flag(inherit.aes=F,aes(country=tolower(rownames(fit6_coefs)),
                y=rownames(fit6_coefs),x=male_typicality_OR_women))+
  geom_point(inherit.aes=F,color=met.brewer("Cassatt2")[8],size=8,
             aes(y=rownames(fit6_coefs),x=male_typicality_OR_men))+
  geom_flag(inherit.aes=F,aes(country=tolower(rownames(fit6_coefs)),
                              y=rownames(fit6_coefs),x=male_typicality_OR_men))
  
png(filename = 
      "results/cntry_male_typicality_OR_plot.png",
    units = "cm",
    width = 21.0,height=29.7*(3/4),res = 300)
cntry.plot
dev.off()

# country-moderation associations

# religiosity

cntry.relig<-
  test_dat %>%
  group_by(cntry) %>%
  summarise(rlgdgr.c.cm=mean(rlgdgr.c))
cntry.relig


test_dat<-left_join(
  x=test_dat,
  y=cntry.relig,
  by=c("cntry"))

# this model does not converge, remove the RE_covs
fit8<-glmer(childless3~gndr.c+FM.z.cmc+gndr.c:FM.z.cmc+
                         married.c+eduyrs.c+rlgdgr.c+age.c+essround.c+
              rlgdgr.c.cm+rlgdgr.c.cm:FM.z.cmc+
              rlgdgr.c.cm:gndr.c+rlgdgr.c.cm:FM.z.cmc:gndr.c+
                         (gndr.c+FM.z.cmc+gndr.c:FM.z.cmc|cntry),
                       data=test_dat,
                       weights = pspwght,family=binomial(link="logit"),
                       control=glmerControl(optimizer="bobyqa",
                                            optCtrl=list(maxfun=2e7)))
summary(fit8)
getFE_glmer(fit8)
getVC(fit8)

# recovs removed

fit8_norecov<-glmer(childless3~gndr.c+FM.z.cmc+gndr.c:FM.z.cmc+
              married.c+eduyrs.c+rlgdgr.c+age.c+essround.c+
              rlgdgr.c.cm+rlgdgr.c.cm:FM.z.cmc+
              rlgdgr.c.cm:gndr.c+rlgdgr.c.cm:FM.z.cmc:gndr.c+
              (gndr.c+FM.z.cmc+gndr.c:FM.z.cmc||cntry),
            data=test_dat,
            weights = pspwght,family=binomial(link="logit"),
            control=glmerControl(optimizer="bobyqa",
                                 optCtrl=list(maxfun=2e7)))
summary(fit8_norecov)
getFE_glmer(fit8_norecov)
getVC(fit8_norecov)

points<-c(mean(cntry.relig$rlgdgr.c.cm)-sd(cntry.relig$rlgdgr.c.cm),
          mean(cntry.relig$rlgdgr.c.cm),
          mean(cntry.relig$rlgdgr.c.cm)+sd(cntry.relig$rlgdgr.c.cm))
          

FM.slopes.by.rlgdgr.c.cm.fit8<-
  emtrends(fit8_norecov,var="FM.z.cmc",specs="rlgdgr.c.cm",
           at=list(rlgdgr.c.cm=points),infer=c(T,T))
FM.slopes.by.rlgdgr.c.cm.fit8
pairs(FM.slopes.by.rlgdgr.c.cm.fit8,
      adjust="none")

FM.gndr.slopes.by.rlgdgr.c.cm.fit8<-
  emtrends(fit8_norecov,var="FM.z.cmc",specs=c("rlgdgr.c.cm","gndr.c"),
           at=list(rlgdgr.c.cm=points,gndr.c=c(-0.5,0.5)),infer=c(T,T))
FM.gndr.slopes.by.rlgdgr.c.cm.fit8

FM.women.slopes.by.rlgdgr.c.cm.fit8<-
  emtrends(fit8_norecov,var="FM.z.cmc",specs=c("rlgdgr.c.cm","gndr.c"),
           at=list(rlgdgr.c.cm=points,gndr.c=c(-0.5)),infer=c(T,T))
FM.women.slopes.by.rlgdgr.c.cm.fit8

pairs(FM.women.slopes.by.rlgdgr.c.cm.fit8,
      adjust="none")

FM.men.slopes.by.rlgdgr.c.cm.fit8<-
  emtrends(fit8_norecov,var="FM.z.cmc",specs=c("rlgdgr.c.cm","gndr.c"),
           at=list(rlgdgr.c.cm=points,gndr.c=c(0.5)),infer=c(T,T))
FM.men.slopes.by.rlgdgr.c.cm.fit8

pairs(FM.men.slopes.by.rlgdgr.c.cm.fit8,
      adjust="none")

#library(remotes)
#remotes::install_github("sebastiansauer/pradadata")
#library(pradadata)
#data(cult_values)
#export(cult_values,"data/cult_values.xlsx",overwrite=T)

cult_values<-import("data/cult_values.xlsx")
cult_values

test_dat<-
  left_join(
    x=test_dat,
    y=cult_values,
    by=c("cntry"="ISO2")
  )

names(test_dat)

# center the country-level variables

cult_values_ESS<-cult_values[cult_values$ISO2 %in% unique(test_dat$cntry),]

test_dat$embedded.c<-test_dat$embedded-mean(cult_values_ESS$embedded,na.rm=T)
embedded_points<-
  c(-1*sd(cult_values_ESS$embedded,na.rm=T),
    0,
    sd(cult_values_ESS$embedded,na.rm=T))


fit9_norecov<-glmer(childless3~gndr.c+FM.z.cmc+gndr.c:FM.z.cmc+
                      married.c+eduyrs.c+rlgdgr.c+age.c+essround.c+
                      embedded.c+embedded.c:FM.z.cmc+
                      embedded.c:gndr.c+embedded.c:FM.z.cmc:gndr.c+
                      (gndr.c+FM.z.cmc+gndr.c:FM.z.cmc||cntry),
                    data=test_dat,
                    weights = pspwght,family=binomial(link="logit"),
                    control=glmerControl(optimizer="bobyqa",
                                         optCtrl=list(maxfun=2e7)))
summary(fit9_norecov)
getFE_glmer(fit9_norecov)
getVC(fit9_norecov)

FM.slopes.by.embedded.c.fit9<-
  emtrends(fit9_norecov,var="FM.z.cmc",specs="embedded.c",
           at=list(embedded.c=points),infer=c(T,T))
FM.slopes.by.embedded.c.fit9
pairs(FM.slopes.by.embedded.c.fit9,
      adjust="none")

FM.gndr.slopes.by.embedded.c.fit9<-
  emtrends(fit9_norecov,var="FM.z.cmc",specs=c("embedded.c","gndr.c"),
           at=list(embedded.c=points,gndr.c=c(-0.5,0.5)),infer=c(T,T))
FM.gndr.slopes.by.embedded.c.fit9

FM.women.slopes.by.embedded.c.fit9<-
  emtrends(fit9_norecov,var="FM.z.cmc",specs=c("embedded.c","gndr.c"),
           at=list(embedded.c=points,gndr.c=c(-0.5)),infer=c(T,T))
FM.women.slopes.by.embedded.c.fit9

pairs(FM.women.slopes.by.embedded.c.fit9,
      adjust="none")

FM.men.slopes.by.embedded.c.fit9<-
  emtrends(fit9_norecov,var="FM.z.cmc",specs=c("embedded.c","gndr.c"),
           at=list(embedded.c=points,gndr.c=c(0.5)),infer=c(T,T))
FM.men.slopes.by.embedded.c.fit9

pairs(FM.men.slopes.by.embedded.c.fit9,
      adjust="none")

# test region domicile rurality
table(test_dat$rural.c,useNA="always")
t1<-Sys.time()
fit10<-glmer(childless3~gndr.c+FM.z.cmc+gndr.c:FM.z.cmc+
                         married.c+eduyrs.c+rlgdgr.c+age.c+essround.c+
               rural.c+rural.c:FM.z.cmc+rural.c:gndr.c+rural.c:FM.z.cmc:gndr.c+
                         (gndr.c+FM.z.cmc+gndr.c:FM.z.cmc|cntry),
                       data=test_dat,
                       weights = pspwght,family=binomial(link="logit"),
                       control=glmerControl(optimizer="bobyqa",
                                            optCtrl=list(maxfun=2e7)))
t2<-Sys.time()
t2-t1
summary(fit10)
getFE_glmer(fit10)
getVC(fit10)

rural.slopes.fit10<-
  emtrends(fit10,var="rural.c",specs="rural.c",
           at=list(rural.c=0),infer=c(T,T))
rural.slopes.fit10
exp(c(-0.15,-0.181,-0.118))

FM.rural.slopes.fit10<-
  emtrends(fit10,var="FM.z.cmc",specs="rural.c",
           at=list(rural.c=c(-0.5,0.5)),infer=c(T,T))
FM.rural.slopes.fit10
contrast(FM.rural.slopes.fit10,method=list(eff=c(-1,1)),infer=c(T,T))
exp(c(-0.00945,-0.0429,0.024))

gndr.rural.slopes.fit10<-
  emtrends(fit10,var="rural.c",specs="gndr.c",
           at=list(gndr.c=c(-0.5,0.5)),infer=c(T,T))
gndr.rural.slopes.fit10
contrast(gndr.rural.slopes.fit10,method=list(eff=c(-1,1)),infer=c(T,T))
exp(c(0.132,0.0691,0.195))
exp(c(-0.2157,-0.261,-0.1705))
exp(c(-0.0834,-0.128,-0.0388))

FM.gndr.rural.slopes.fit10<-
  emtrends(fit10,var="FM.z.cmc",specs=c("rural.c","gndr.c"),
           at=list(rural.c=c(-0.5,0.5),
                   gndr.c=c(-0.5,0.5)),infer=c(T,T))

contr.FM.gndr.rural.slopes.fit10<-
  contrast(FM.gndr.rural.slopes.fit10,
         method=list(women=c(-1,1,0,0),
                     men=c(0,0,-1,1)),
                     adjust="none")
contrast(contr.FM.gndr.rural.slopes.fit10,
        method=list(eff=c(-1,1)),adjust="none",infer=c(T,T))
exp(c(0.00891,-0.0575,0.0753))

# variance explained

performance::r2(fit0)

#library(remotes)
#install_github("timnewbold/StatisticalModels")
#install.packages("MuMIn")

library(StatisticalModels)
library(MuMIn)
library(sjstats)

fit11_only_time<-glmer(childless3~essround.c+
                         (1|cntry),
                       data=test_dat,
                       weights = pspwght,family=binomial(link="logit"),
                       control=glmerControl(optimizer="bobyqa",
                                            optCtrl=list(maxfun=2e7)))

R2GLMER(fit11_only_time)
r.squaredGLMM(fit11_only_time)
performance::r2(fit11_only_time)
str(performance::r2(fit11_only_time))

fit11_only_covariates<-glmer(childless3~essround.c+
                               married.c+eduyrs.c+rlgdgr.c+age.c+gndr.c+
                         (1|cntry),
                       data=test_dat,
                       weights = pspwght,family=binomial(link="logit"),
                       control=glmerControl(optimizer="bobyqa",
                                            optCtrl=list(maxfun=2e7)))

R2GLMER(fit11_only_covariates)$marginal-R2GLMER(fit11_only_time)$marginal
performance::r2(fit11_only_covariates)
performance::r2(fit11_only_covariates)$R2_marginal-
  performance::r2(fit11_only_time)$R2_marginal


fit11_FM_gndr<-glmer(childless3~essround.c+FM.z.cmc+FM.z.cmc:gndr.c+
                       married.c+eduyrs.c+rlgdgr.c+age.c+gndr.c+
                       (1|cntry),
                     data=test_dat,
                     weights = pspwght,family=binomial(link="logit"),
                     control=glmerControl(optimizer="bobyqa",
                                          optCtrl=list(maxfun=2e7)))

fit11_FM_main<-glmer(childless3~essround.c+FM.z.cmc+
                       married.c+eduyrs.c+rlgdgr.c+age.c+gndr.c+
                       (1|cntry),
                     data=test_dat,
                     weights = pspwght,family=binomial(link="logit"),
                     control=glmerControl(optimizer="bobyqa",
                                          optCtrl=list(maxfun=2e7)))

fit11_FM_main_noFM<-glmer(childless3~essround.c+
                       married.c+eduyrs.c+rlgdgr.c+age.c+gndr.c+
                       (1|cntry),
                     data=test_dat,
                     weights = pspwght,family=binomial(link="logit"),
                     control=glmerControl(optimizer="bobyqa",
                                          optCtrl=list(maxfun=2e7)))

fit11_FM_main_nomarriage<-glmer(childless3~essround.c+FM.z.cmc+
                       eduyrs.c+rlgdgr.c+age.c+gndr.c+
                       (1|cntry),
                     data=test_dat,
                     weights = pspwght,family=binomial(link="logit"),
                     control=glmerControl(optimizer="bobyqa",
                                          optCtrl=list(maxfun=2e7)))

fit11_FM_main_noedu<-glmer(childless3~essround.c+FM.z.cmc+
                       married.c+rlgdgr.c+age.c+gndr.c+
                       (1|cntry),
                     data=test_dat,
                     weights = pspwght,family=binomial(link="logit"),
                     control=glmerControl(optimizer="bobyqa",
                                          optCtrl=list(maxfun=2e7)))

fit11_FM_main_norelig<-glmer(childless3~essround.c+FM.z.cmc+
                       married.c+eduyrs.c+age.c+gndr.c+
                       (1|cntry),
                     data=test_dat,
                     weights = pspwght,family=binomial(link="logit"),
                     control=glmerControl(optimizer="bobyqa",
                                          optCtrl=list(maxfun=2e7)))

fit11_FM_main_noage<-glmer(childless3~essround.c+FM.z.cmc+
                       married.c+eduyrs.c+rlgdgr.c+gndr.c+
                       (1|cntry),
                     data=test_dat,
                     weights = pspwght,family=binomial(link="logit"),
                     control=glmerControl(optimizer="bobyqa",
                                          optCtrl=list(maxfun=2e7)))



library(sjPlot)
tab_model(fit0,fit11_only_time,
          fit11_FM_gndr,fit11_FM_main,
          fit11_FM_main_noage,fit11_FM_main_noedu,
          fit11_FM_main_noFM,fit11_FM_main_nomarriage,
          fit11_FM_main_norelig)

# try single added variables instead

fit11_gndr<-glmer(childless3~essround.c+gndr.c+
                       (1|cntry),
                     data=test_dat,
                     weights = pspwght,family=binomial(link="logit"),
                     control=glmerControl(optimizer="bobyqa",
                                          optCtrl=list(maxfun=2e7)))

fit11_married<-glmer(childless3~essround.c+married.c+
                    (1|cntry),
                  data=test_dat,
                  weights = pspwght,family=binomial(link="logit"),
                  control=glmerControl(optimizer="bobyqa",
                                       optCtrl=list(maxfun=2e7)))

fit11_edu<-glmer(childless3~essround.c+eduyrs.c+
                    (1|cntry),
                  data=test_dat,
                  weights = pspwght,family=binomial(link="logit"),
                  control=glmerControl(optimizer="bobyqa",
                                       optCtrl=list(maxfun=2e7)))

fit11_FM<-glmer(childless3~essround.c+FM.z.cmc+
                    (1|cntry),
                  data=test_dat,
                  weights = pspwght,family=binomial(link="logit"),
                  control=glmerControl(optimizer="bobyqa",
                                       optCtrl=list(maxfun=2e7)))

fit11_relig<-glmer(childless3~essround.c+rlgdgr.c+
                  (1|cntry),
                data=test_dat,
                weights = pspwght,family=binomial(link="logit"),
                control=glmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e7)))

fit11_age<-glmer(childless3~essround.c+age.c+
                     (1|cntry),
                   data=test_dat,
                   weights = pspwght,family=binomial(link="logit"),
                   control=glmerControl(optimizer="bobyqa",
                                        optCtrl=list(maxfun=2e7)))

fit11_FM_gndr<-glmer(childless3~essround.c+gndr.c+FM.z.cmc+
                   (1|cntry),
                 data=test_dat,
                 weights = pspwght,family=binomial(link="logit"),
                 control=glmerControl(optimizer="bobyqa",
                                      optCtrl=list(maxfun=2e7)))

fit11_FM_gndr_int<-glmer(childless3~essround.c+gndr.c+FM.z.cmc+gndr.c:FM.z.cmc-
                       (1|cntry),
                     data=test_dat,
                     weights = pspwght,family=binomial(link="logit"),
                     control=glmerControl(optimizer="bobyqa",
                                          optCtrl=list(maxfun=2e7)))

fit11_full_fixed<-glmer(childless3~essround.c+FM.z.cmc+FM.z.cmc:gndr.c+
                       married.c+eduyrs.c+rlgdgr.c+age.c+gndr.c+
                       (1|cntry),
                     data=test_dat,
                     weights = pspwght,family=binomial(link="logit"),
                     control=glmerControl(optimizer="bobyqa",
                                          optCtrl=list(maxfun=2e7)))

tab_model(fit0,fit11_only_time,
          fit11_gndr,fit11_edu,fit11_age,
          fit11_relig,fit11_married,
          fit11_FM,fit11_FM_gndr,fit11_FM_gndr_int,fit11_full_fixed,
          file="results/R2_tabs.html")
table(test_dat$agea)

# Education interactions

fit12<-glmer(childless3~essround.c+gndr.c+FM.z.cmc+FM.z.cmc:gndr.c+
                          married.c+eduyrs.c+rlgdgr.c+age.c+gndr.c+
               eduyrs.c:gndr.c+eduyrs.c:FM.z.cmc+eduyrs.c:gndr.c:FM.z.cmc+
                          (gndr.c+FM.z.cmc+FM.z.cmc:gndr.c||cntry),
                        data=test_dat,
                        weights = pspwght,family=binomial(link="logit"),
                        control=glmerControl(optimizer="bobyqa",
                                             optCtrl=list(maxfun=2e7)))


getFE_glmer(fit12)

edutrends_fit12<-
  emtrends(fit12,var="eduyrs.c",specs="gndr.c",
           at=list(gndr.c=c(-0.5,0.5)),infer=c(T,T))
edutrends_fit12
contrast(edutrends_fit12,method=list(diff=c(-1,1)),infer=c(T,T))
round(exp(c(-0.274,-0.309,-0.24)),2)
round(exp(c(0.529,0.504,0.554)),2)
round(exp(c(0.254,0.230,0.279)),2)

# analyze within nordic countries

unique(test_dat$cntry)

fit12_nordics<-glmer(childless3~essround.c+gndr.c+FM.z.cmc+FM.z.cmc:gndr.c+
               married.c+eduyrs.c+rlgdgr.c+age.c+gndr.c+
               eduyrs.c:gndr.c+eduyrs.c:FM.z.cmc+eduyrs.c:gndr.c:FM.z.cmc+
               (gndr.c+FM.z.cmc+FM.z.cmc:gndr.c||cntry),
             data=test_dat[
               test_dat$cntry=="FI" | test_dat$cntry=="SE" | test_dat$cntry=="NO" |
                 test_dat$cntry=="DK" | test_dat$cntry=="IS",
             ],
             weights = pspwght,family=binomial(link="logit"),
             control=glmerControl(optimizer="bobyqa",
                                  optCtrl=list(maxfun=2e7)))

summary(fit12_nordics)
getFE_glmer(fit12_nordics)

edutrends_fit12_nordics<-
  emtrends(fit12_nordics,var="eduyrs.c",specs="gndr.c",
           at=list(gndr.c=c(-0.5,0.5)),infer=c(T,T))
edutrends_fit12_nordics
contrast(edutrends_fit12_nordics,method=list(diff=c(-1,1)),infer=c(T,T))
round(exp(c(-0.156,-0.246,-0.0657)),2)
round(exp(c(0.331,0.264,0.399)),2)
round(exp(c(0.176,0.113,0.238)),2)

# analyze three-way between time, education, and gender

unique(test_dat$cntry)

fit13_time_edu_gndr<-glmer(childless3~essround.c+gndr.c+FM.z.cmc+FM.z.cmc:gndr.c+
                       married.c+eduyrs.c+rlgdgr.c+age.c+gndr.c+
                         essround.c*eduyrs.c*gndr.c+
                       (1|cntry),
                     data=test_dat,
                     weights = pspwght,family=binomial(link="logit"),
                     control=glmerControl(optimizer="bobyqa",
                                          optCtrl=list(maxfun=2e7)))

summary(fit13_time_edu_gndr)
getFE_glmer(fit13_time_edu_gndr)

edutrends_fit13_time_edu<-
  emtrends(fit13_time_edu_gndr,var="eduyrs.c",specs="essround.c",
           at=list(essround.c=seq(from=-4.5,to=4.5,by=1)),infer=c(T,T))
edutrends_fit13_time_edu
contrast(edutrends_fit13_time_edu,
         method=list(diff=c(-1,0,0,0,0,0,0,0,0,1)),infer=c(T,T))
round(exp(c(-0.00143,-0.00747,0.0046)),2)
round(exp(c(-0.0129,-0.0672,0.0414)),2)

edutrends_fit13_time_edu_gndr<-
  emtrends(fit13_time_edu_gndr,var="eduyrs.c",specs=c("essround.c","gndr.c"),
           at=list(essround.c=seq(from=-4.5,to=4.5,by=1),
                   gndr.c=c(-0.5,0.5)),infer=c(T,T))
edutrends_fit13_time_edu_gndr

edutrends_fit13_time_edu_men<-
  emtrends(fit13_time_edu_gndr,var="eduyrs.c",specs=c("essround.c","gndr.c"),
           at=list(essround.c=seq(from=-4.5,to=4.5,by=1),
                   gndr.c=c(0.5)),infer=c(T,T))
  
edutrends_fit13_time_edu_women<-
  emtrends(fit13_time_edu_gndr,var="eduyrs.c",specs=c("essround.c","gndr.c"),
           at=list(essround.c=seq(from=-4.5,to=4.5,by=1),
                   gndr.c=c(-0.5)),infer=c(T,T))

edutrends_fit13_time_edu_men
edutrends_fit13_time_edu_women

contrast_fit13_time_edu_men<-
  contrast(edutrends_fit13_time_edu_men,
           method=list(trend_men=c(-1,1,0,0,0,0,0,0,0,0)),infer=c(T,T))

contrast_fit13_time_edu_women<-
  contrast(edutrends_fit13_time_edu_women,
           method=list(trend_women=c(-1,1,0,0,0,0,0,0,0,0)),infer=c(T,T))
contrast_fit13_time_edu_women

rbind(contrast_fit13_time_edu_women,
    contrast_fit13_time_edu_men,adjust="none")

contrast(rbind(contrast_fit13_time_edu_women,
               contrast_fit13_time_edu_men,adjust="none"),
         method=list(c(-1,1)),infer=c(T,T))

round(exp(c(0.0108,-0.00107,0.0227)),2)
round(exp(c(0.331,0.264,0.399)),2)
round(exp(c(0.176,0.113,0.238)),2)


# Study 1

diff_dat<-
  value_typ$preds

diff_dat$age.c<-as.numeric(diff_dat$age.c)
diff_dat$rlgdgr.c<-as.numeric(diff_dat$rlgdgr.c)
diff_dat$eduyrs.c<-as.numeric(diff_dat$eduyrs.c)

# data exclusions

diff_dat<-diff_dat %>%
  dplyr::filter(waves>1 &
                  !is.na(gndr))

table(diff_dat$waves)

hist(diff_dat$pred)
FM_pooled_sd<-value_typ$D[1,"pooled.sd.total"]
diff_dat$FM.z<-diff_dat$pred/FM_pooled_sd

# modeling

mod0<-lmer(FM.z~(1|cntry),data=diff_dat,REML=F,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod0)
getVC(mod0)

mod1<-lmer(FM.z~gndr.c+(1|cntry),data=diff_dat,REML=F,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod1)
getFE(mod1)
getVC(mod1)

mod2<-lmer(FM.z~gndr.c+(gndr.c|cntry),data=diff_dat,REML=F,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2)
getFE(mod2)
getVC(mod2)

mean(sum_tab$d.sd.total)

lvl2_var_cond_lvl1(mod2,lvl1.var = "gndr.c",lvl1.values = c(0.5,-0.5))

mod2_norecov<-lmer(FM.z~gndr.c+(gndr.c||cntry),data=diff_dat,REML=F,
                   control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod2_norecov)
getFE(mod2_norecov)
getVC(mod2_norecov)

anova(mod2_norecov,mod2)

# add time

mod3<-lmer(FM.z~gndr.c+essround.c+(gndr.c|cntry),data=diff_dat,REML=F,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod3)
getFE(mod3)
getVC(mod3)

# time as random

mod4<-lmer(FM.z~gndr.c+essround.c+(gndr.c+essround.c|cntry),data=diff_dat,REML=F,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod4)
getFE(mod4)
getVC(mod4)

anova(mod3,mod4)
coefficients(mod4)$cntry

# interaction between gender and time

mod5<-lmer(FM.z~gndr.c+essround.c+
             gndr.c:essround.c+(gndr.c+essround.c|cntry),data=diff_dat,REML=F,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod5)
getFE(mod5)
getVC(mod5)

# interaction between gender and time as random

mod6<-lmer(FM.z~gndr.c+essround.c+
             gndr.c:essround.c+(gndr.c+essround.c+gndr.c:essround.c|cntry),
           data=diff_dat,REML=F,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6)
getFE(mod6)
getVC(mod6)
anova(mod5,mod6)

# obtain predicted values

time_trends_mod6<-
  emtrends(mod6,var="essround.c",specs="gndr.c",
           at=list(gndr.c=c(-0.5,0.5)),
           disable.pbkrtest=T,
           lmerTest.limit = 400000,infer=c(T,T))

time_trends_mod6
pairs(time_trends_mod6,adjust="none")



# fit a fixed model for dadas

# scale the essround.c to obtain a standardized coef
diff_dat$essround.z<-(diff_dat$essround.c-mean(diff_dat$essround.c))/
  sd(rep(unique(diff_dat$essround.c),1000))

mod7<-lmer(FM.z~gndr.c+essround.z+
             gndr.c:essround.z+(gndr.c+essround.z+
                                  gndr.c:essround.z|cntry),
           data=diff_dat,REML=F,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod7)
round(getFE(mod7),3)
round(getFE(mod5),3)
getVC(mod7)


# scaled estimates

lvl2_var_cond_lvl1(mod2,lvl1.var = "gndr.c",lvl1.values = c(0.5,-0.5))

time_trends_mod7<-
  emtrends(mod7,var="essround.z",specs="gndr.c",
           at=list(gndr.c=c(0.5,-0.5)),
           disable.pbkrtest=T,
           lmerTest.limit = 400000,infer=c(T,T))
time_trends_mod7
-0.0029/0.2377514
0.0246/0.2497135

change_mod6<-emmeans(mod6,specs="essround.c",by="gndr.c",
                     at=list(gndr.c=c(-0.5,0.5),
                             essround.c=rev(range(diff_dat$essround.c))),
                     disable.pbkrtest=T,
                     lmerTest.limit = 400000,infer=c(T,T))
change_mod6
pairs(change_mod6,adjust="none")

# change in differences

change_in_diff_mod6<-emmeans(mod6,specs=c("gndr.c","essround.c"),
                     at=list(gndr.c=c(-0.5,0.5),
                             essround.c=rev(range(diff_dat$essround.c))),
                     disable.pbkrtest=T,
                     lmerTest.limit = 400000,infer=c(T,T))
change_in_diff_mod6
pairs(change_mod6,adjust="none")

changes_in_diff <- list(
  diff_ESS10 = c(-1,1,0,0),
  diff_ESS1 = c(0,0,-1,1))

diff_mod6<-contrast(change_in_diff_mod6,method = changes_in_diff,adjust="none")
diff_mod6
pairs(diff_mod6)


mod7_dadas<-ml_dadas(model=mod7,predictor = "essround.z",
                     diff_var = "gndr.c",diff_var_values = c(0.5,-0.5),
                     scaled_estimates = T,re_cov_test = T)
round(mod7_dadas$dadas,3)
round(mod7_dadas$scaled_estimates,3)
round(mod7_dadas$vpc_at_reduced,3)

# fit a mod6 without random effects


mod6_nore<-lmer(FM.z~gndr.c+essround.c+
             gndr.c:essround.c+(gndr.c|cntry),
           data=diff_dat,REML=F,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_nore)
getFE(mod6_nore)
getVC(mod6_nore)

time_trends_mod6_nore<-
  emtrends(mod6_nore,var="essround.c",specs="gndr.c",
           at=list(gndr.c=c(-0.5,0.5)),
           disable.pbkrtest=T,
           lmerTest.limit = 400000,infer=c(T,T))

getVC(mod6)
time_trends_mod6_nore
time_trends_mod6
pairs(time_trends_mod6_nore,adjust="none")

change_in_diff_mod6_nore<-
  emmeans(mod6_nore,specs=c("gndr.c","essround.c"),
          at=list(gndr.c=c(-0.5,0.5),
                  essround.c=rev(range(diff_dat$essround.c))),
          disable.pbkrtest=T,
          lmerTest.limit = 400000,infer=c(T,T))
change_in_diff_mod6_nore

changes_in_diff <- list(
  diff_ESS10 = c(-1,1,0,0),
  diff_ESS1 = c(0,0,-1,1))

diff_mod6_nore<-
  contrast(change_in_diff_mod6_nore,
           method = changes_in_diff,adjust="none")
diff_mod6_nore
pairs(diff_mod6_nore)


# try model categorical time
table(diff_dat$essround)
diff_dat$essround.f<-as.factor(diff_dat$essround)


mod8<-lmer(FM.z~gndr.c+essround.f+
             gndr.c:essround.f+(gndr.c+essround.f+
                                  gndr.c:essround.f|cntry),
           data=diff_dat,REML=F,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod8)
round(getFE(mod8),3)
getVC(mod8)

# remove the recovs

mod8_norecov<-lmer(FM.z~gndr.c+essround.f+
             gndr.c:essround.f+(gndr.c+essround.f+
                                  gndr.c:essround.f||cntry),
           data=diff_dat,REML=F,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))

# 



# test divergence for values

value.vars

diff_dat$resp<-rowMeans(diff_dat[,value.vars])
diff_dat$pow.c<-diff_dat$pow-diff_dat$resp
diff_dat$pow.z<-diff_dat$pow.c/sd(diff_dat$pow.c)

mod6_pow<-lmer(pow.z~gndr.c+essround.c+
             gndr.c:essround.c+(gndr.c+essround.c+gndr.c:essround.c||cntry),
           data=diff_dat,REML=F,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_pow)
getFE(mod6_pow)

time_trends_mod6_pow<-
  emtrends(mod6_pow,var="essround.c",specs="gndr.c",
           at=list(gndr.c=c(-0.5,0.5)),
           disable.pbkrtest=T,
           lmerTest.limit = 400000,infer=c(T,T))

time_trends_mod6_pow
pairs(time_trends_mod6_pow,adjust="none")

change_in_diff_mod6_pow<-
  emmeans(mod6_pow,specs=c("gndr.c","essround.c"),
          at=list(gndr.c=c(-0.5,0.5),
                  essround.c=rev(range(diff_dat$essround.c))),
          disable.pbkrtest=T,
          lmerTest.limit = 400000,infer=c(T,T))
change_in_diff_mod6_pow

changes_in_diff <- list(
  diff_ESS10 = c(-1,1,0,0),
  diff_ESS1 = c(0,0,-1,1))

diff_mod6_pow<-
  contrast(change_in_diff_mod6_pow,
           method = changes_in_diff,adjust="none")
diff_mod6_pow
pairs(diff_mod6_pow)

mod6_pow_nore<-lmer(pow.z~gndr.c+essround.c+
                 gndr.c:essround.c+(gndr.c|cntry),
               data=diff_dat,REML=F,
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_pow_nore)
getFE(mod6_pow_nore)

time_trends_mod6_pow_nore<-
  emtrends(mod6_pow_nore,var="essround.c",specs="gndr.c",
           at=list(gndr.c=c(-0.5,0.5)),
           disable.pbkrtest=T,
           lmerTest.limit = 400000,infer=c(T,T))

time_trends_mod6_pow_nore
pairs(time_trends_mod6_pow_nore,adjust="none")

change_in_diff_mod6_pow_nore<-
  emmeans(mod6_pow_nore,specs=c("gndr.c","essround.c"),
                             at=list(gndr.c=c(-0.5,0.5),
                                     essround.c=rev(range(diff_dat$essround.c))),
                             disable.pbkrtest=T,
                             lmerTest.limit = 400000,infer=c(T,T))
change_in_diff_mod6_pow_nore

changes_in_diff <- list(
  diff_ESS10 = c(-1,1,0,0),
  diff_ESS1 = c(0,0,-1,1))

diff_mod6_pow_nore<-
  contrast(change_in_diff_mod6_pow_nore,
           method = changes_in_diff,adjust="none")
diff_mod6_pow_nore
pairs(diff_mod6_pow_nore)

# achievement

diff_dat$ach.c<-diff_dat$ach-diff_dat$resp
diff_dat$ach.z<-diff_dat$ach.c/sd(diff_dat$ach.c)

mod6_ach<-lmer(ach.z~gndr.c+essround.c+
                 gndr.c:essround.c+(gndr.c+essround.c+gndr.c:essround.c||cntry),
               data=diff_dat,REML=F,
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_ach)
getFE(mod6_ach)

time_trends_mod6_ach<-
  emtrends(mod6_ach,var="essround.c",specs="gndr.c",
           at=list(gndr.c=c(-0.5,0.5)),
           disable.pbkrtest=T,
           lmerTest.limit = 400000,infer=c(T,T))

time_trends_mod6_ach
pairs(time_trends_mod6_ach,adjust="none")

change_in_diff_mod6_ach<-
  emmeans(mod6_ach,specs=c("gndr.c","essround.c"),
          at=list(gndr.c=c(-0.5,0.5),
                  essround.c=rev(range(diff_dat$essround.c))),
          disable.pbkrtest=T,
          lmerTest.limit = 400000,infer=c(T,T))
change_in_diff_mod6_ach

changes_in_diff <- list(
  diff_ESS10 = c(-1,1,0,0),
  diff_ESS1 = c(0,0,-1,1))

diff_mod6_ach<-
  contrast(change_in_diff_mod6_ach,
           method = changes_in_diff,adjust="none")
diff_mod6_ach
pairs(diff_mod6_ach)

mod6_ach_nore<-lmer(ach.z~gndr.c+essround.c+
                      gndr.c:essround.c+(gndr.c|cntry),
                    data=diff_dat,REML=F,
                    control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_ach_nore)
getFE(mod6_ach_nore)

time_trends_mod6_ach_nore<-
  emtrends(mod6_ach_nore,var="essround.c",specs="gndr.c",
           at=list(gndr.c=c(-0.5,0.5)),
           disable.pbkrtest=T,
           lmerTest.limit = 400000,infer=c(T,T))

time_trends_mod6_ach_nore
pairs(time_trends_mod6_ach_nore,adjust="none")

change_in_diff_mod6_ach_nore<-
  emmeans(mod6_ach_nore,specs=c("gndr.c","essround.c"),
          at=list(gndr.c=c(-0.5,0.5),
                  essround.c=rev(range(diff_dat$essround.c))),
          disable.pbkrtest=T,
          lmerTest.limit = 400000,infer=c(T,T))
change_in_diff_mod6_ach_nore

changes_in_diff <- list(
  diff_ESS10 = c(-1,1,0,0),
  diff_ESS1 = c(0,0,-1,1))

diff_mod6_ach_nore<-
  contrast(change_in_diff_mod6_ach_nore,
           method = changes_in_diff,adjust="none")
diff_mod6_ach_nore
pairs(diff_mod6_ach_nore)

# benevolence

diff_dat$ben.c<-diff_dat$ben-diff_dat$resp
diff_dat$ben.z<-diff_dat$ben.c/sd(diff_dat$ben.c)

mod6_ben<-lmer(ben.z~gndr.c+essround.c+
                 gndr.c:essround.c+(gndr.c+essround.c+gndr.c:essround.c||cntry),
               data=diff_dat,REML=F,
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_ben)
getFE(mod6_ben)

time_trends_mod6_ben<-
  emtrends(mod6_ben,var="essround.c",specs="gndr.c",
           at=list(gndr.c=c(-0.5,0.5)),
           disable.pbkrtest=T,
           lmerTest.limit = 400000,infer=c(T,T))

time_trends_mod6_ben
pairs(time_trends_mod6_ben,adjust="none")

change_in_diff_mod6_ben<-
  emmeans(mod6_ben,specs=c("gndr.c","essround.c"),
          at=list(gndr.c=c(-0.5,0.5),
                  essround.c=rev(range(diff_dat$essround.c))),
          disable.pbkrtest=T,
          lmerTest.limit = 400000,infer=c(T,T))
change_in_diff_mod6_ben

changes_in_diff <- list(
  diff_ESS10 = c(-1,1,0,0),
  diff_ESS1 = c(0,0,-1,1))

diff_mod6_ben<-
  contrast(change_in_diff_mod6_ben,
           method = changes_in_diff,adjust="none")
diff_mod6_ben
pairs(diff_mod6_ben)

mod6_ben_nore<-lmer(ben.z~gndr.c+essround.c+
                      gndr.c:essround.c+(gndr.c|cntry),
                    data=diff_dat,REML=F,
                    control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_ben_nore)
getFE(mod6_ben_nore)

time_trends_mod6_ben_nore<-
  emtrends(mod6_ben_nore,var="essround.c",specs="gndr.c",
           at=list(gndr.c=c(-0.5,0.5)),
           disable.pbkrtest=T,
           lmerTest.limit = 400000,infer=c(T,T))

time_trends_mod6_ben_nore
pairs(time_trends_mod6_ben_nore,adjust="none")

change_in_diff_mod6_ben_nore<-
  emmeans(mod6_ben_nore,specs=c("gndr.c","essround.c"),
          at=list(gndr.c=c(-0.5,0.5),
                  essround.c=rev(range(diff_dat$essround.c))),
          disable.pbkrtest=T,
          lmerTest.limit = 400000,infer=c(T,T))
change_in_diff_mod6_ben_nore

changes_in_diff <- list(
  diff_ESS10 = c(-1,1,0,0),
  diff_ESS1 = c(0,0,-1,1))

diff_mod6_ben_nore<-
  contrast(change_in_diff_mod6_ben_nore,
           method = changes_in_diff,adjust="none")
diff_mod6_ben_nore
pairs(diff_mod6_ben_nore)

# stimulation

diff_dat$sti.c<-diff_dat$sti-diff_dat$resp
diff_dat$sti.z<-diff_dat$sti.c/sd(diff_dat$sti.c)

mod6_sti<-lmer(sti.z~gndr.c+essround.c+
                 gndr.c:essround.c+(gndr.c+essround.c+gndr.c:essround.c||cntry),
               data=diff_dat,REML=F,
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_sti)
getFE(mod6_sti)

time_trends_mod6_sti<-
  emtrends(mod6_sti,var="essround.c",specs="gndr.c",
           at=list(gndr.c=c(-0.5,0.5)),
           disable.pbkrtest=T,
           lmerTest.limit = 400000,infer=c(T,T))

time_trends_mod6_sti
pairs(time_trends_mod6_sti,adjust="none")

change_in_diff_mod6_sti<-
  emmeans(mod6_sti,specs=c("gndr.c","essround.c"),
          at=list(gndr.c=c(-0.5,0.5),
                  essround.c=rev(range(diff_dat$essround.c))),
          disable.pbkrtest=T,
          lmerTest.limit = 400000,infer=c(T,T))
change_in_diff_mod6_sti

changes_in_diff <- list(
  diff_ESS10 = c(-1,1,0,0),
  diff_ESS1 = c(0,0,-1,1))

diff_mod6_sti<-
  contrast(change_in_diff_mod6_sti,
           method = changes_in_diff,adjust="none")
diff_mod6_sti
pairs(diff_mod6_sti)

mod6_sti_nore<-lmer(sti.z~gndr.c+essround.c+
                      gndr.c:essround.c+(gndr.c|cntry),
                    data=diff_dat,REML=F,
                    control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_sti_nore)
getFE(mod6_sti_nore)

time_trends_mod6_sti_nore<-
  emtrends(mod6_sti_nore,var="essround.c",specs="gndr.c",
           at=list(gndr.c=c(-0.5,0.5)),
           disable.pbkrtest=T,
           lmerTest.limit = 400000,infer=c(T,T))

time_trends_mod6_sti_nore
pairs(time_trends_mod6_sti_nore,adjust="none")

change_in_diff_mod6_sti_nore<-
  emmeans(mod6_sti_nore,specs=c("gndr.c","essround.c"),
          at=list(gndr.c=c(-0.5,0.5),
                  essround.c=rev(range(diff_dat$essround.c))),
          disable.pbkrtest=T,
          lmerTest.limit = 400000,infer=c(T,T))
change_in_diff_mod6_sti_nore

changes_in_diff <- list(
  diff_ESS10 = c(-1,1,0,0),
  diff_ESS1 = c(0,0,-1,1))

diff_mod6_sti_nore<-
  contrast(change_in_diff_mod6_sti_nore,
           method = changes_in_diff,adjust="none")
diff_mod6_sti_nore
pairs(diff_mod6_sti_nore)


# universalism

diff_dat$uni.c<-diff_dat$uni-diff_dat$resp
diff_dat$uni.z<-diff_dat$uni.c/sd(diff_dat$uni.c)

mod6_uni<-lmer(uni.z~gndr.c+essround.c+
                 gndr.c:essround.c+(gndr.c+essround.c+gndr.c:essround.c||cntry),
               data=diff_dat,REML=F,
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_uni)
getFE(mod6_uni)

time_trends_mod6_uni<-
  emtrends(mod6_uni,var="essround.c",specs="gndr.c",
           at=list(gndr.c=c(-0.5,0.5)),
           disable.pbkrtest=T,
           lmerTest.limit = 400000,infer=c(T,T))

time_trends_mod6_uni
pairs(time_trends_mod6_uni,adjust="none")

change_in_diff_mod6_uni<-
  emmeans(mod6_uni,specs=c("gndr.c","essround.c"),
          at=list(gndr.c=c(-0.5,0.5),
                  essround.c=rev(range(diff_dat$essround.c))),
          disable.pbkrtest=T,
          lmerTest.limit = 400000,infer=c(T,T))
change_in_diff_mod6_uni

changes_in_diff <- list(
  diff_ESS10 = c(-1,1,0,0),
  diff_ESS1 = c(0,0,-1,1))

diff_mod6_uni<-
  contrast(change_in_diff_mod6_uni,
           method = changes_in_diff,adjust="none")
diff_mod6_uni
pairs(diff_mod6_uni)

mod6_uni_nore<-lmer(uni.z~gndr.c+essround.c+
                      gndr.c:essround.c+(gndr.c|cntry),
                    data=diff_dat,REML=F,
                    control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod6_uni_nore)
getFE(mod6_uni_nore)

time_trends_mod6_uni_nore<-
  emtrends(mod6_uni_nore,var="essround.c",specs="gndr.c",
           at=list(gndr.c=c(-0.5,0.5)),
           disable.pbkrtest=T,
           lmerTest.limit = 400000,infer=c(T,T))

time_trends_mod6_uni_nore
pairs(time_trends_mod6_uni_nore,adjust="none")

change_in_diff_mod6_uni_nore<-
  emmeans(mod6_uni_nore,specs=c("gndr.c","essround.c"),
          at=list(gndr.c=c(-0.5,0.5),
                  essround.c=rev(range(diff_dat$essround.c))),
          disable.pbkrtest=T,
          lmerTest.limit = 400000,infer=c(T,T))
change_in_diff_mod6_uni_nore

changes_in_diff <- list(
  diff_ESS10 = c(-1,1,0,0),
  diff_ESS1 = c(0,0,-1,1))

diff_mod6_uni_nore<-
  contrast(change_in_diff_mod6_uni_nore,
           method = changes_in_diff,adjust="none")
diff_mod6_uni_nore
pairs(diff_mod6_uni_nore)



# GGGI

GGGI<-import("data/GGGI_2006_2022.xlsx")
str(GGGI)

table(GGGI$ISO2 %in% unique(diff_dat$cntry))

GGGI.vars<-
  c(paste0("X",2006:2018),paste0("X",2020:2022))
GGGI.vars



GGGI<-GGGI %>%
  dplyr::filter(Indicator=="Overall Global Gender Gap Index" &
                  Subindicator.Type == "Index")

GGGI$GGGI.m<-
  rowMeans(GGGI[,GGGI.vars],na.rm=T)
names(GGGI)
GGGI<-GGGI[GGGI$ISO2 %in% unique(diff_dat$cntry),]

GGGI$GGGI.m.z<-scale(GGGI$GGGI.m,center=T,scale=T)

GEP_dat<-
  left_join(x=diff_dat,
            y=GGGI,
            by=c("cntry"="ISO2"))


mod_GEP<-lmer(FM.z~gndr.c+essround.c+GGGI.m.z+GGGI.m.z:gndr.c+(gndr.c|cntry),
           data=GEP_dat,REML=F,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
summary(mod_GEP)
getFE(mod_GEP)
getVC(mod_GEP)

GEP_dadas<-
  ml_dadas(model=mod_GEP,predictor = "GGGI.m.z",
           diff_var = "gndr.c",diff_var_values = c(-0.5,0.5),
           scaled_estimates = T,re_cov_test = T)
round(GEP_dadas$dadas,3)
round(GEP_dadas$vpc_at_reduced,3)
round(GEP_dadas$scaled_estimates,3)
GEP_dadas$re_cov_test
