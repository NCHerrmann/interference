###R Script for Culbertson and Herrmann###
#TABLE OF CONTENTS#
#Line 19 = data summary and prep
#Line 122 = main analyses of male-male trials
#Line 296 = main analyses of green anole focal lizards (comparing sexes)
#Line 418 = supplemental analyses

#Set working directory, clean environment, and load pacakges
setwd("")
rm(list=ls())
options(stringsAsFactors = FALSE)
library(rstanarm)
library(ggplot2)
library(grid)
#Import data
dat<-read.csv(file="DATA.csv",header=TRUE)

###########################
###DATA SUMMARY AND PREP###
###########################

#TABLE 1 - SUMMARY OF TRIALS
sagrei=subset(dat,dat$focspec=="SA") #34 trials with focal A sagrei males
carmales=subset(dat,dat$focspec=="CA" & dat$focsex=="M") #28 trials with focal A carolinensis males
carfemales=subset(dat,dat$focspec=="CA" & dat$focsex=="F") #30 trials with focal A carolinensis females
#trial duration (minutes)
mean(sagrei$duration) #3.3
sd(sagrei$duration) #1.3
mean(carmales$duration) #3.3
sd(carmales$duration) #1.4
mean(carfemales$duration) #3.6
sd(carfemales$duration) #1.6
#perch height (cm)
mean(sagrei$height) #84
sd(sagrei$height) #58
mean(carmales$height) #274
sd(carmales$height) #133
mean(carfemales$height) #205
sd(carfemales$height) #125
#was intruder larger? (binary)
mean(sagrei$intlarger,na.rm = TRUE) #18%
mean(carmales$intlarger,na.rm = TRUE) #68%
#did intruder display or move during trial? (binary)
mean(sagrei$intreact) #68%
mean(carmales$intreact) #54%
mean(carfemales$intreact) #50%
#was at least one conspecific present?
mean(sagrei$con) #44%
mean(carmales$con) #11%
mean(carfemales$con) #10%
#did focal lizard attack? (binary)
mean(sagrei$attack) #6%
mean(carmales$attack) #0%
mean(carfemales$attack) #0%
#did focal lizard display? (binary)
mean(sagrei$display) #91%
mean(carmales$display) #57%
mean(carfemales$display) #50%
#focal lizard displays/minute
mean(sagrei$disprate) #2.1
mean(carmales$disprate) #1.4
mean(carfemales$disprate) #0.5
#did focal lizard retreat? (binary)
mean(sagrei$scatter) #9%
mean(carmales$scatter) #79%
mean(carfemales$scatter) #77%
#if focal lizardretreated, did it retreat upward? (binary)
sagrei$up[sagrei$scatter==1] #1 of 3
mean(sagrei$up[sagrei$scatter==1]) #33%
carmales$up[carmales$scatter==1] #15 of 22
mean(carmales$up[carmales$scatter==1]) #68%
carfemales$up[carfemales$scatter==1] #12 of 23
mean(carfemales$up[carfemales$scatter==1]) #52%

#Subset data for analyses
males=subset(dat,dat$focsex=="M" & !is.na(dat$intlarger)) #56 male-male trials in which larger lizard identified with certainty
malessvldiff=subset(males,!is.na(males$svldiff)) #50 male-male trials in which exact svl difference is known
malesfirst=males[males$firsttrial==1,] #30 male-male trials including only the first trial for each tethered individual
car=subset(dat,dat$focspec=="CA") #58 trials with A carolinensis focal lizards
carfirst=subset(car, car$firsttrial==1) #40 carolinesis trials including only the first trial for each tethered individual
carmoved=subset(car,car$scatter==1) #45 trials in which focal A carolinensis retreated
carfirstmoved=subset(carfirst,carfirst$scatter==1) #29 trials in which focal A carolinensis retreated and including only the first trial for each tethered individual

#Prep data for analysis
#This section standardizes all predcitor variables by centering and dividing by 2SD
#Standarized data subsets are labeled with a "Z" even though variables are technically not Z-scored (Z-scoring requires dividng by only 1 SD)
malesZ=males
malesZsvldiff=malessvldiff
malesZfirst=malesfirst
carZ=car
carZfirst=carfirst
carZmoved=carmoved
carZfirstmoved=carfirstmoved

#functions to center predictor variables and dividing by 2 SD
prep=function(x){as.numeric(as.factor(x))-1} #required for binary predictors only
doubleZ=function(x){(x-mean(x, na.rm=TRUE))/(2*sd(x, na.rm=TRUE))} #required for all predictors

#predictors that are binary except focspec and focsex
predsBIN=c("site","placement","intlarger","intreact","con","firsttrial")

#all predictors
predsALL=c(predsBIN, "svldiff", "height", "duration")

#applies functions to appropriate predictors across all data subsets
malesZ[predsBIN] <- lapply(malesZ[predsBIN],prep)
malesZ[predsALL] <- lapply(malesZ[predsALL],doubleZ)
malesZsvldiff[predsBIN] <- lapply(malesZsvldiff[predsBIN],prep)
malesZsvldiff[predsALL] <- lapply(malesZsvldiff[predsALL],doubleZ)
malesZfirst[predsBIN] <- lapply(malesZfirst[predsBIN],prep)
malesZfirst[predsALL] <- lapply(malesZfirst[predsALL],doubleZ)
carZ[predsBIN] <- lapply(carZ[predsBIN],prep)
carZ[predsALL] <- lapply(carZ[predsALL],doubleZ)
carZfirst[predsBIN] <- lapply(carZfirst[predsBIN],prep)
carZfirst[predsALL] <- lapply(carZfirst[predsALL],doubleZ)
carZmoved[predsBIN] <- lapply(carZmoved[predsBIN],prep)
carZmoved[predsALL] <- lapply(carZmoved[predsALL],doubleZ)
carZfirstmoved[predsBIN] <- lapply(carZfirstmoved[predsBIN],prep)
carZfirstmoved[predsALL] <- lapply(carZfirstmoved[predsALL],doubleZ)

###############################
###MAIN ANALYSES - male-male###
###############################

#Do green and brown males differ in THE PROBABILITY OF DISPLAY?

#Fit model
displaymodel <- stan_glm(display ~ focspec + site + intlarger +
                           intreact + height + placement + con + duration
                         ,family=binomial(link='logit'),data=malesZ, iter=4000)
summary(displaymodel)
plot(displaymodel)
prior_summary(displaymodel)

#Create data frame to hold posterior predictions
preds1display <- data.frame(focspec=rep(unique(malesZ$focspec)),
                            site=rep(0,2),
                            intlarger=rep(0,2),
                            height=rep(0,2),
                            con=rep(0,2),
                            intreact=rep(0,2),
                            placement=rep(0,2),
                            duration=rep(0,2))
preds1display

#Generate posterior predictions
lpdisplay <- posterior_linpred(displaymodel,newdata=preds1display,transform=TRUE)
quantfun <- function(x){quantile(x,probs = c(0.025,0.25,0.5,0.75,0.975))}
lpdisplay_quants <- t(apply(lpdisplay,FUN=quantfun,MARGIN=2))
preds1display$lwr_025 <- lpdisplay_quants[,1]
preds1display$lwr_25 <- lpdisplay_quants[,2]
preds1display$median <- lpdisplay_quants[,3]
preds1display$upr_75 <- lpdisplay_quants[,4]
preds1display$upr_975 <- lpdisplay_quants[,5]
preds1display$focspec <- c("CA (n = 22)", "SA (n = 34)")

#Plot posterior predictions
fig1display <- ggplot(preds1display,aes(colour=focspec))+
  geom_linerange(aes(x=focspec,ymin=lwr_025,ymax=upr_975))+
  geom_linerange(aes(x=focspec,ymin=lwr_25,ymax=upr_75),lwd=1.5)+
  geom_point(aes(x=focspec,y=median),shape=21,size=2,fill="white")+
  scale_color_manual(values=c("forestgreen","saddlebrown"))+
  coord_flip()+
  scale_y_continuous("Probability of display",limits=c(0,1))+
  theme_bw()
fig1display +
  theme(legend.position="none") +
  theme(axis.title.y=element_blank(),
        axis.title.x = element_text(size=8, margin = margin(t = 8, r = 0, b = 0, l = 0)))



#Do green and brown males differ in THE RATE OF DISPLAY?

#Fit model
displayratemodel <- stan_glm(disprate ~ focspec + site + intlarger +
                               intreact + height + placement + con + duration
                             ,family=gaussian(link = "identity"),data=malesZ, iter=4000)
summary(displayratemodel)
plot(displayratemodel)
prior_summary(displayratemodel)

#Create data frame to hold posterior predictions
preds1numdisplay <- data.frame(focspec=rep(unique(malesZ$focspec)),
                               site=rep(0,2),
                               intlarger=rep(0,2),
                               height=rep(0,2),
                               con=rep(0,2),
                               intreact=rep(0,2),
                               placement=rep(0,2),
                               duration=rep(0,2))

preds1numdisplay

#Generate posterior predictions
lpnumdisplay <- posterior_linpred(displayratemodel,newdata=preds1numdisplay,transform=TRUE)
quantfun <- function(x){quantile(x,probs = c(0.025,0.25,0.5,0.75,0.975))}
lpnumdisplay_quants <- t(apply(lpnumdisplay,FUN=quantfun,MARGIN=2))
preds1numdisplay$lwr_025 <- lpnumdisplay_quants[,1]
preds1numdisplay$lwr_25 <- lpnumdisplay_quants[,2]
preds1numdisplay$median <- lpnumdisplay_quants[,3]
preds1numdisplay$upr_75 <- lpnumdisplay_quants[,4]
preds1numdisplay$upr_975 <- lpnumdisplay_quants[,5]
preds1numdisplay$focspec <- c("", " ")
preds1numdisplay

#Plot posterior predictions
fig1displayrate <- ggplot(preds1numdisplay,aes(colour=focspec))+
  geom_linerange(aes(x=focspec,ymin=lwr_025,ymax=upr_975))+
  geom_linerange(aes(x=focspec,ymin=lwr_25,ymax=upr_75),lwd=1.5)+
  geom_point(aes(x=focspec,y=median),shape=21,size=2,fill="white")+
  scale_color_manual(values=c("forestgreen","saddlebrown"))+
  coord_flip()+
  scale_y_continuous("Displays per min",limits=c(0,2.7))+
  theme_bw()
fig1displayrate +
  theme(legend.position="none") +
  theme(axis.title.y=element_blank(),
        axis.title.x = element_text(size=8, margin = margin(t = 8, r = 0, b = 0, l = 0)))

#Do green and brown males differ in PROBABILITY OF RETREATING?
scatmodel <- stan_glm(scatter ~ focspec + site + intlarger +
                        intreact + height + placement
                      ,family=binomial(link='logit'),data=malesZ, iter=4000)
summary(scatmodel)
plot(scatmodel)
prior_summary(scatmodel)

#Create data frame to hold posterior predictions
preds1scat <- data.frame(focspec=rep(unique(malesZ$focspec)),
                         site=rep(0,2),
                         intlarger=rep(0,2),
                         intreact=rep(0,2),
                         height=rep(0,2),
                         placement=rep(0,2))
preds1scat

#Generate posterior predictions
lpscat <- posterior_linpred(scatmodel,newdata=preds1scat,transform=TRUE)
quantfun <- function(x){quantile(x,probs = c(0.025,0.25,0.5,0.75,0.975))}
lpscat_quants <- t(apply(lpscat,FUN=quantfun,MARGIN=2))
preds1scat$lwr_025 <- lpscat_quants[,1]
preds1scat$lwr_25 <- lpscat_quants[,2]
preds1scat$median <- lpscat_quants[,3]
preds1scat$upr_75 <- lpscat_quants[,4]
preds1scat$upr_975 <- lpscat_quants[,5]
preds1scat$focspec <- c("", " ")


#Plot posterior predictions
fig1scat <- ggplot(preds1scat,aes(colour=focspec))+
  geom_linerange(aes(x=focspec,ymin=lwr_025,ymax=upr_975))+
  geom_linerange(aes(x=focspec,ymin=lwr_25,ymax=upr_75),lwd=1.5)+
  geom_point(aes(x=focspec,y=median),shape=21,size=2,fill="white")+
  scale_color_manual(values=c("forestgreen","saddlebrown"))+
  coord_flip()+
  scale_y_continuous("Probability of retreat",limits=c(0,1))+
  theme_bw()
fig1scat +
  theme(legend.position="none") +
  theme(axis.title.y=element_blank(),
        axis.title.x = element_text(size=8, margin = margin(t = 8, r = 0, b = 0, l = 0)))

#FIGURE 2 - SPECIES-SPECIFIC DIFFERENCES (males)
jpeg("Figure2.jpeg",width=9.5,height=2, units="in",res=1000)
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,3, widths=c(3.5,3,3))))

vplayout<-function(x,y) viewport(layout.pos.row=x,layout.pos.col=y)

print(fig1display +
        theme(legend.position="none") +
        theme(axis.title.y=element_blank(),
              axis.title.x = element_text(size=9, margin = margin(t = 9, r = 0, b = 0, l = 0))),
      vp=vplayout(1,1))

print(fig1displayrate +
        theme(legend.position="none") +
        theme(axis.title.y=element_blank(),
              axis.title.x = element_text(size=9, margin = margin(t = 9, r = 0, b = 0, l = 0))),
      vp=vplayout(1,2))

print(fig1scat +
        theme(legend.position="none") +
        theme(axis.title.y=element_blank(),
              axis.title.x = element_text(size=9, margin = margin(t = 9, r = 0, b = 0, l = 0))),
      vp=vplayout(1,3))

dev.off()

#Are male green anoles that retreat more likely than chance to retreat upward?
carmalesretreat=subset(car,car$focsex=="M" & car$scatter==1) #22 trials
binom.test(sum(carmalesretreat$up),length(carmalesretreat$up),p=0.5,alternative="greater") #p-value = 0.0669

##################################
###MAIN ANALYSES - green anoles###
##################################

#Do A. carolinensis males and females differ in the PROBABILITY OF RETREAT?

#Fit model
carscatmodel <- stan_glm(scatter ~ focsex + site + 
                           intreact + height + placement
                         ,family=binomial(link='logit'),data=carZ, iter=4000)
summary(carscatmodel)
plot(carscatmodel)
prior_summary(carscatmodel)

#Create data frame to hold posterior predictions
preds1scatCAR <- data.frame(focsex=rep(unique(carZ$focsex)),
                            site=rep(0,2),
                            intreact=rep(0,2),
                            height=rep(0,2),
                            placement=rep(0,2))
preds1scatCAR

#Generate posterior predictions
lp <- posterior_linpred(carscatmodel,newdata=preds1scatCAR,transform=TRUE)
quantfun <- function(x){quantile(x,probs = c(0.025,0.25,0.5,0.75,0.975))}
lp_quants <- t(apply(lp,FUN=quantfun,MARGIN=2))
preds1scatCAR$lwr_025 <- lp_quants[,1]
preds1scatCAR$lwr_25 <- lp_quants[,2]
preds1scatCAR$median <- lp_quants[,3]
preds1scatCAR$upr_75 <- lp_quants[,4]
preds1scatCAR$upr_975 <- lp_quants[,5]
preds1scatCAR$focsex <- c("Females (n = 30)","Males (n = 28)")

#Plot posterior predictions
carfig1scat <- ggplot(preds1scatCAR,aes(colour=focsex))+
  geom_linerange(aes(x=focsex,ymin=lwr_025,ymax=upr_975))+
  geom_linerange(aes(x=focsex,ymin=lwr_25,ymax=upr_75),lwd=1.5)+
  geom_point(aes(x=focsex,y=median),shape=21,size=2,fill="white")+
  coord_flip()+
  scale_color_manual(values=c("forestgreen","green3"))+
  scale_y_continuous("Probability of retreat",limits=c(0.5,1))+
  theme_bw()
carfig1scat +
  theme(legend.position="none") +
  theme(axis.title.y=element_blank(),
        axis.title.x = element_text(size=10, margin = margin(t = 14, r = 0, b = 0, l = 0)))

#Do carolinensis males and females that retreat differ in the PROBABILITY OF MOVING UP?
length(carmoved$up[carmoved$focsex=="M"]) #22
length(carmoved$up[carmoved$focsex=="F"]) #23
mean(carmoved$up[carmoved$focsex=="M"]) #68%
mean(carmoved$up[carmoved$focsex=="F"]) #52%

#Fit model
carupmodel <- stan_glm(up ~ focsex + site +
                         intreact + height + placement
                       ,family=binomial(link='logit'),data=carZmoved, iter=4000)

summary(carupmodel)
plot(carupmodel)
prior_summary(carupmodel)

#Create data frame to hold posterior predictions
preds1upCAR <- data.frame(focsex=rep(unique(carZmoved$focsex)),
                          site=rep(0,2),
                          intreact=rep(0,2),
                          height=rep(0,2),
                          placement=rep(0,2))
preds1upCAR

#Generate posterior predictions
lpup <- posterior_linpred(carupmodel,newdata=preds1upCAR,transform=TRUE)
quantfun <- function(x){quantile(x,probs = c(0.025,0.25,0.5,0.75,0.975))}
lpup_quants <- t(apply(lpup,FUN=quantfun,MARGIN=2))
preds1upCAR$lwr_025 <- lpup_quants[,1]
preds1upCAR$lwr_25 <- lpup_quants[,2]
preds1upCAR$median <- lpup_quants[,3]
preds1upCAR$upr_75 <- lpup_quants[,4]
preds1upCAR$upr_975 <- lpup_quants[,5]
preds1upCAR$focsex <- c("", " ")

#Plot posterior predictions
carfig1up <- ggplot(preds1upCAR,aes(colour=focsex))+
  geom_linerange(aes(x=focsex,ymin=lwr_025,ymax=upr_975))+
  geom_linerange(aes(x=focsex,ymin=lwr_25,ymax=upr_75),lwd=1.5)+
  geom_point(aes(x=focsex,y=median),shape=21,size=2,fill="white")+
  coord_flip()+
  scale_color_manual(values=c("forestgreen","green3"))+
  scale_y_continuous("Probability of moving up (if retreating)",limits=c(0.3,0.9))+
  theme_bw()
carfig1up  +
  theme(legend.position="none") +
  theme(axis.title.y=element_blank(),
        axis.title.x = element_text(size=10, margin = margin(t = 14, r = 0, b = 0, l = 0)))

#Are female green anoles that retreat more likely than chance to retreat upward?
carfemalesretreat=subset(car,car$focsex=="F" & car$scatter==1) #23 trials
binom.test(sum(carfemalesretreat$up),length(carfemalesretreat$up),p=0.5,alternative="greater") #p-value = 0.5

#FIGURE 3 - sex-specific probabilities (greens)

jpeg("Figure3.jpeg",width=6.85,height=2, units="in",res=1000)
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,2, widths=c(3.85,3))))

vplayout<-function(x,y) viewport(layout.pos.row=x,layout.pos.col=y)

print(carfig1scat +
        theme(legend.position="none") +
        theme(axis.title.y=element_blank(),
              axis.title.x = element_text(size=9, margin = margin(t = 9, r = 0, b = 0, l = 0))),
      vp=vplayout(1,1))

print(carfig1up +
        theme(legend.position="none") +
        theme(axis.title.y=element_blank(),
              axis.title.x = element_text(size=9, margin = margin(t = 9, r = 0, b = 0, l = 0))),
      vp=vplayout(1,2))

dev.off()


############################
###SUPPLEMENTAL ANALYSES ###
############################

###
#Rechecking male-male models using SVL difference instead of "intlarger" to account for body size difference
###

#Probability of display
displaymodelSVLdiff <- stan_glm(display ~ focspec + site + svldiff +
                                  intreact + height + placement + con + duration
                                ,family=binomial(link='logit'),data=malesZsvldiff, iter=4000)
summary(displaymodelSVLdiff)
plot(displaymodelSVLdiff)
preds1displaySVLdiff <- data.frame(focspec=rep(unique(malesZsvldiff$focspec)),
                                   site=rep(0,2),
                                   svldiff=rep(0,2),
                                   height=rep(0,2), 
                                   con=rep(0,2),
                                   intreact=rep(0,2),
                                   placement=rep(0,2),
                                   duration=rep(0,2))
lp <- posterior_linpred(displaymodelSVLdiff,newdata=preds1displaySVLdiff,transform=TRUE)
quantfun <- function(x){quantile(x,probs = c(0.025,0.25,0.5,0.75,0.975))}
lp_quants <- t(apply(lp,FUN=quantfun,MARGIN=2))
preds1displaySVLdiff$lwr_025 <- lp_quants[,1]
preds1displaySVLdiff$lwr_25 <- lp_quants[,2]
preds1displaySVLdiff$median <- lp_quants[,3]
preds1displaySVLdiff$upr_75 <- lp_quants[,4]
preds1displaySVLdiff$upr_975 <- lp_quants[,5]
preds1displaySVLdiff$focspec <- c("CA (n = 18)", "SA (n = 32)")
fig1displaySVLdiff <- ggplot(preds1displaySVLdiff,aes(colour=focspec))+
  geom_linerange(aes(x=focspec,ymin=lwr_025,ymax=upr_975))+
  geom_linerange(aes(x=focspec,ymin=lwr_25,ymax=upr_75),lwd=1.5)+
  geom_point(aes(x=focspec,y=median),shape=21,size=2,fill="white")+
  scale_color_manual(values=c("forestgreen","saddlebrown"))+
  coord_flip()+
  scale_y_continuous("Probability of Display",limits=c(0,1))+
  theme_bw()

#Display rate
displayratemodelSVLdiff <- stan_glm(disprate ~ focspec + site + svldiff +
                                      intreact + height + placement + con + duration
                                    ,family=gaussian(link = "identity"),data=malesZsvldiff, iter=4000)
summary(displayratemodelSVLdiff)
plot(displayratemodelSVLdiff)
preds1displayrateSVLdiff <- data.frame(focspec=rep(unique(malesZsvldiff$focspec)),
                                       site=rep(0,2),
                                       svldiff=rep(0,2),
                                       height=rep(0,2),
                                       con=rep(0,2),
                                       intreact=rep(0,2),
                                       placement=rep(0,2),
                                       duration=rep(0,2))

preds1displayrateSVLdiff
lp <- posterior_linpred(displayratemodelSVLdiff,newdata=preds1displayrateSVLdiff,transform=TRUE)
quantfun <- function(x){quantile(x,probs = c(0.025,0.25,0.5,0.75,0.975))}
lp_quants <- t(apply(lp,FUN=quantfun,MARGIN=2))
preds1displayrateSVLdiff$lwr_025 <- lp_quants[,1]
preds1displayrateSVLdiff$lwr_25 <- lp_quants[,2]
preds1displayrateSVLdiff$median <- lp_quants[,3]
preds1displayrateSVLdiff$upr_75 <- lp_quants[,4]
preds1displayrateSVLdiff$upr_975 <- lp_quants[,5]
preds1displayrateSVLdiff$focspec <- c("CA (n = 18)", "SA (n = 32)")
preds1displayrateSVLdiff
fig1displayrateSVLdiff <- ggplot(preds1displayrateSVLdiff,aes(colour=focspec))+
  geom_linerange(aes(x=focspec,ymin=lwr_025,ymax=upr_975))+
  geom_linerange(aes(x=focspec,ymin=lwr_25,ymax=upr_75),lwd=1.5)+
  geom_point(aes(x=focspec,y=median),shape=21,size=2,fill="white")+
  coord_flip()+
  scale_color_manual(values=c("forestgreen","saddlebrown"))+
  scale_y_continuous("Displays/min",limits=c(0,3))+
  theme_bw()

#Probably of retreat
scatmodelsvl <- stan_glm(scatter ~ focspec + site + svldiff +
                           intreact + height + placement
                         ,family=binomial(link='logit'),data=malesZsvldiff, iter=4000)
summary(scatmodelsvl)
plot(scatmodelsvl)
preds1scatsvl <- data.frame(focspec=rep(unique(malesZsvldiff$focspec)),
                            site=rep(0,2),
                            intreact=rep(0,2),
                            height=rep(0,2),
                            svldiff=rep(0,2),
                            placement=rep(0,2))
lpsvl <- posterior_linpred(scatmodelsvl,newdata=preds1scatsvl,transform=TRUE)
quantfun <- function(x){quantile(x,probs = c(0.025,0.25,0.5,0.75,0.975))}
lpsvl_quants <- t(apply(lpsvl,FUN=quantfun,MARGIN=2))
preds1scatsvl$lwr_025 <- lpsvl_quants[,1]
preds1scatsvl$lwr_25 <- lpsvl_quants[,2]
preds1scatsvl$median <- lpsvl_quants[,3]
preds1scatsvl$upr_75 <- lpsvl_quants[,4]
preds1scatsvl$upr_975 <- lpsvl_quants[,5]
preds1scatsvl$focspec <- c("CA (n = 18)", "SA (n = 32)")
fig1scatSVLdiff <- ggplot(preds1scatsvl,aes(colour=focspec))+
  geom_linerange(aes(x=focspec,ymin=lwr_025,ymax=upr_975))+
  geom_linerange(aes(x=focspec,ymin=lwr_25,ymax=upr_75),lwd=1.5)+
  geom_point(aes(x=focspec,y=median),shape=21,size=2,fill="white")+
  scale_color_manual(values=c("forestgreen","saddlebrown"))+
  coord_flip()+
  scale_y_continuous("Probability of Retreat",limits=c(0,1))+
  theme_bw()

###
#Rechecking male-male models using only trials in which an individual intruder was used for the first time
###

#Probability of display
displaymodelft <- stan_glm(display ~ focspec + site + intlarger +
                             intreact + height + placement + con + duration
                           ,family=binomial(link='logit'),data=malesZfirst, iter=4000)
summary(displaymodelft)
plot(displaymodelft)
preds1displayft <- data.frame(focspec=rep(unique(malesZfirst$focspec)),
                              site=rep(0,2),
                              intlarger=rep(0,2),
                              height=rep(0,2), 
                              con=rep(0,2),
                              intreact=rep(0,2),
                              placement=rep(0,2),
                              duration=rep(0,2))
lp <- posterior_linpred(displaymodelft,newdata=preds1displayft,transform=TRUE)
quantfun <- function(x){quantile(x,probs = c(0.025,0.25,0.5,0.75,0.975))}
lp_quants <- t(apply(lp,FUN=quantfun,MARGIN=2))
preds1displayft$lwr_025 <- lp_quants[,1]
preds1displayft$lwr_25 <- lp_quants[,2]
preds1displayft$median <- lp_quants[,3]
preds1displayft$upr_75 <- lp_quants[,4]
preds1displayft$upr_975 <- lp_quants[,5]
preds1displayft$focspec <- c("CA (n = 16)", "SA (n = 14)")
fig1displayft <- ggplot(preds1displayft,aes(colour=focspec))+
  geom_linerange(aes(x=focspec,ymin=lwr_025,ymax=upr_975))+
  geom_linerange(aes(x=focspec,ymin=lwr_25,ymax=upr_75),lwd=1.5)+
  geom_point(aes(x=focspec,y=median),shape=21,size=2,fill="white")+
  scale_color_manual(values=c("forestgreen","saddlebrown"))+
  coord_flip()+
  scale_y_continuous("Probability of Display",limits=c(0,1))+
  theme_bw()

#Display rate
displayratemodelft <- stan_glm(disprate ~ focspec + site + intlarger +
                                 intreact + height + placement + con + duration
                               ,family=gaussian(link = "identity"),data=malesZfirst, iter=4000)
summary(displayratemodelft)
plot(displayratemodelft)
preds1displayrateft <- data.frame(focspec=rep(unique(malesZfirst$focspec)),
                                  site=rep(0,2),
                                  intlarger=rep(0,2),
                                  height=rep(0,2),
                                  con=rep(0,2),
                                  intreact=rep(0,2),
                                  placement=rep(0,2),
                                  duration=rep(0,2))

preds1displayrateft
lp <- posterior_linpred(displayratemodelft,newdata=preds1displayrateft,transform=TRUE)
quantfun <- function(x){quantile(x,probs = c(0.025,0.25,0.5,0.75,0.975))}
lp_quants <- t(apply(lp,FUN=quantfun,MARGIN=2))
preds1displayrateft$lwr_025 <- lp_quants[,1]
preds1displayrateft$lwr_25 <- lp_quants[,2]
preds1displayrateft$median <- lp_quants[,3]
preds1displayrateft$upr_75 <- lp_quants[,4]
preds1displayrateft$upr_975 <- lp_quants[,5]
preds1displayrateft$focspec <- c("CA (n = 16)", "SA (n = 14)")
preds1displayrateft$lwr_025[1]=0
fig1displayrateft <- ggplot(preds1displayrateft,aes(colour=focspec))+
  geom_linerange(aes(x=focspec,ymin=lwr_025,ymax=upr_975))+
  geom_linerange(aes(x=focspec,ymin=lwr_25,ymax=upr_75),lwd=1.5)+
  geom_point(aes(x=focspec,y=median),shape=21,size=2,fill="white")+
  coord_flip()+
  scale_color_manual(values=c("forestgreen","saddlebrown"))+
  scale_y_continuous("Displays/min",limits=c(0,4))+
  theme_bw()

#Probability of retreat
scatmodelft <- stan_glm(scatter ~ focspec + site 
                        + intlarger +
                          intreact + height + placement
                        ,family=binomial(link='logit'),data=malesZfirst, iter=4000)
summary(scatmodelft)
plot(scatmodelft)
preds1scatft <- data.frame(focspec=rep(unique(malesZfirst$focspec)),
                           site=rep(0,2),
                           intreact=rep(0,2),
                           height=rep(0,2),
                           intlarger=rep(0,2),
                           placement=rep(0,2))
lpft <- posterior_linpred(scatmodelft,newdata=preds1scatft,transform=TRUE)
quantfun <- function(x){quantile(x,probs = c(0.025,0.25,0.5,0.75,0.975))}
lpft_quants <- t(apply(lpft,FUN=quantfun,MARGIN=2))
preds1scatft$lwr_025 <- lpft_quants[,1]
preds1scatft$lwr_25 <- lpft_quants[,2]
preds1scatft$median <- lpft_quants[,3]
preds1scatft$upr_75 <- lpft_quants[,4]
preds1scatft$upr_975 <- lpft_quants[,5]
preds1scatft$focspec <- c("CA (n = 16)", "SA (n = 14)")
fig1scatft <- ggplot(preds1scatft,aes(colour=focspec))+
  geom_linerange(aes(x=focspec,ymin=lwr_025,ymax=upr_975))+
  geom_linerange(aes(x=focspec,ymin=lwr_25,ymax=upr_75),lwd=1.5)+
  geom_point(aes(x=focspec,y=median),shape=21,size=2,fill="white")+
  scale_color_manual(values=c("forestgreen","saddlebrown"))+
  coord_flip()+
  scale_y_continuous("Probability of Retreat",limits=c(0,1))+
  theme_bw()

###
#Rechecking green anole (sex difference) models using only trials in which an individual intruder was used for the first time
###

#Probability of retreat
carscatmodelft <- stan_glm(scatter ~ focsex + site + 
                             intreact + height + placement
                           ,family=binomial(link='logit'),data=carZfirst, iter=4000)
summary(carscatmodelft)
plot(carscatmodelft)
preds1scatCARft <- data.frame(focsex=rep(unique(carZfirst$focsex)),
                              site=rep(0,2),
                              intreact=rep(0,2),
                              height=rep(0,2),
                              placement=rep(0,2))
preds1scatCARft
lp <- posterior_linpred(carscatmodelft,newdata=preds1scatCARft,transform=TRUE)
quantfun <- function(x){quantile(x,probs = c(0.025,0.25,0.5,0.75,0.975))}
lp_quants <- t(apply(lp,FUN=quantfun,MARGIN=2))
preds1scatCARft$lwr_025 <- lp_quants[,1]
preds1scatCARft$lwr_25 <- lp_quants[,2]
preds1scatCARft$median <- lp_quants[,3]
preds1scatCARft$upr_75 <- lp_quants[,4]
preds1scatCARft$upr_975 <- lp_quants[,5]
preds1scatCARft$focsex <- c("Females (n = 22)", "Males (n = 18)")
carfig1scatft <- ggplot(preds1scatCARft,aes(colour=focsex))+
  geom_linerange(aes(x=focsex,ymin=lwr_025,ymax=upr_975))+
  geom_linerange(aes(x=focsex,ymin=lwr_25,ymax=upr_75),lwd=1.5)+
  geom_point(aes(x=focsex,y=median),shape=21,size=2,fill="white")+
  coord_flip()+
  scale_color_manual(values=c("forestgreen","green3"))+
  scale_y_continuous("Probability of Retreat",limits=c(0.5,1))+
  theme_bw()
carfig1scatft +
  theme(legend.position="none") +
  theme(axis.title.y=element_blank(),
        axis.title.x = element_text(size=10, margin = margin(t = 14, r = 0, b = 0, l = 0)))

#Probaility of moving upward (if retreating)
length(carfirstmoved$up[carfirstmoved$focsex=="M"]) #13
length(carfirstmoved$up[carfirstmoved$focsex=="F"]) #16
carupmodelft <- stan_glm(up ~ focsex + site +
                           intreact + height + placement
                         ,family=binomial(link='logit'),data=carZfirstmoved, iter=4000)
summary(carupmodelft)
plot(carupmodelft)
preds1upCARft <- data.frame(focsex=rep(unique(carZfirstmoved$focsex)),
                            site=rep(0,2),
                            intreact=rep(0,2),
                            height=rep(0,2),
                            placement=rep(0,2))
lpup <- posterior_linpred(carupmodelft,newdata=preds1upCARft,transform=TRUE)
quantfun <- function(x){quantile(x,probs = c(0.025,0.25,0.5,0.75,0.975))}
lpup_quants <- t(apply(lpup,FUN=quantfun,MARGIN=2))
preds1upCARft$lwr_025 <- lpup_quants[,1]
preds1upCARft$lwr_25 <- lpup_quants[,2]
preds1upCARft$median <- lpup_quants[,3]
preds1upCARft$upr_75 <- lpup_quants[,4]
preds1upCARft$upr_975 <- lpup_quants[,5]
preds1upCARft$focsex <- c("Females (n = 16)", "Males (n= 13)")
carfig1upft <- ggplot(preds1upCARft,aes(colour=focsex))+
  geom_linerange(aes(x=focsex,ymin=lwr_025,ymax=upr_975))+
  geom_linerange(aes(x=focsex,ymin=lwr_25,ymax=upr_75),lwd=1.5)+
  geom_point(aes(x=focsex,y=median),shape=21,size=2,fill="white")+
  coord_flip()+
  scale_color_manual(values=c("forestgreen","green3"))+
  scale_y_continuous("Probability of Moving Up (if retreating)",limits=c(0,1))+
  theme_bw()

#SUPPLEMENATRY FIGURES

#Supplmentary Figure - male-male retreat model posterior predictions based on species and size

preds2scat <- data.frame(focspec=rep(unique(malesZ$focspec),2),
                         site=rep(0,4),
                         intreact=rep(0,4),
                         height=rep(0,4),
                         placement=rep(0,4),
                         intlarger=c(rep(min(malesZ$intlarger),2),
                                     rep(max(malesZ$intlarger),2)))
preds2scat

lp <- posterior_linpred(scatmodel,newdata=preds2scat,transform=TRUE)
quantfun <- function(x){quantile(x,probs = c(0.025,0.25,0.5,0.75,0.975))}
lp_quants <- t(apply(lp,FUN=quantfun,MARGIN=2))
preds2scat$lwr_025 <- lp_quants[,1]
preds2scat$lwr_25 <- lp_quants[,2]
preds2scat$median <- lp_quants[,3]
preds2scat$upr_75 <- lp_quants[,4]
preds2scat$upr_975 <- lp_quants[,5]
preds2scat$FocalLizard <- c("larger CA (n = 7)", "larger SA (n = 28)", "smaller CA (n = 15)", "smaller SA (n = 6)")

fig2scat <- ggplot(preds2scat)+
  geom_linerange(aes(x=FocalLizard,ymin=lwr_025,ymax=upr_975))+
  geom_linerange(aes(x=FocalLizard,intlarger,ymin=lwr_25,ymax=upr_75),lwd=1.5)+
  geom_point(aes(x=FocalLizard,y=median),shape=21,size=2,fill="white")+
  coord_flip()+
  scale_y_continuous("Probability of Retreat",limits=c(0,1))+
  theme_bw()

pdf("Supplement Fig 1 - retreat by spec and size.pdf",width=4,height=2)
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,1)))

vplayout<-function(x,y) viewport(layout.pos.row=x,layout.pos.col=y)

print(fig2scat +
        theme(legend.position="none") +
        theme(axis.title.y=element_blank(),
              axis.title.x = element_text(size=8, margin = margin(t = 5, r = 0, b = 0, l = 0))),
      vp=vplayout(1,1))

dev.off()

#Supplmentary Figure - green anole up model posterior predictions based on sex and intruder placement

sum(carmoved$focsex=="F" & carmoved$placement=="B") #11
sum(carmoved$focsex=="M" & carmoved$placement=="B") #12
sum(carmoved$focsex=="F" & carmoved$placement=="A") #12
sum(carmoved$focsex=="M" & carmoved$placement=="A") #10

preds2upCAR <- data.frame(focsex=rep(unique(carZmoved$focsex),2),
                          site=rep(0,4),
                          intreact=rep(0,4),
                          height=rep(0,4),
                          placement=c(rep(min(carZmoved$placement),2),
                                      rep(max(carZmoved$placement),2)))
preds2upCAR

lp2up <- posterior_linpred(carupmodel,newdata=preds2upCAR,transform=TRUE)
quantfun <- function(x){quantile(x,probs = c(0.025,0.25,0.5,0.75,0.975))}
lp2up_quants <- t(apply(lp2up,FUN=quantfun,MARGIN=2))
preds2upCAR$lwr_025 <- lp2up_quants[,1]
preds2upCAR$lwr_25 <- lp2up_quants[,2]
preds2upCAR$median <- lp2up_quants[,3]
preds2upCAR$upr_75 <- lp2up_quants[,4]
preds2upCAR$upr_975 <- lp2up_quants[,5]
preds2upCAR$focsex <- c("F w/ intruder below (n = 11)", "M w/ intruder below (n = 12)",
                        "F w/ intruder above (n = 12)","M w/ intruder above (n = 10)")

carfig2up <- ggplot(preds2upCAR)+
  geom_linerange(aes(x=focsex,ymin=lwr_025,ymax=upr_975))+
  geom_linerange(aes(x=focsex,ymin=lwr_25,ymax=upr_75),lwd=1.5)+
  geom_point(aes(x=focsex,y=median),shape=21,size=2,fill="white")+
  coord_flip()+
  scale_y_continuous("Probability of Moving Up (if retreating)",limits=c(0,1))+
  theme_bw()


pdf("Supplement Fig 2 - direction by sex and intruder placement.pdf",width=4,height=2)
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,1)))

vplayout<-function(x,y) viewport(layout.pos.row=x,layout.pos.col=y)

print(carfig2up +
        theme(legend.position="none") +
        theme(axis.title.y=element_blank(),
              axis.title.x = element_text(size=8, margin = margin(t = 5, r = 0, b = 0, l = 0))),
      vp=vplayout(1,1))

dev.off()

#Supplmentary Figure - species-specific probabilities (males) when using SVLdiff instead of intlarger

pdf("Supplement Fig 3 - SVLdiff.pdf",width=4,height=6)
grid.newpage()
pushViewport(viewport(layout=grid.layout(3,1)))

vplayout<-function(x,y) viewport(layout.pos.row=x,layout.pos.col=y)

print(fig1displaySVLdiff +
        theme(legend.position="none") +
        theme(axis.title.y=element_blank(),
              axis.title.x = element_text(size=8, margin = margin(t = 5, r = 0, b = 0, l = 0))),
      vp=vplayout(1,1))

print(fig1displayrateSVLdiff +
        theme(legend.position="none") +
        theme(axis.title.y=element_blank(),
              axis.title.x = element_text(size=8, margin = margin(t = 5, r = 0, b = 0, l = 0))),
      vp=vplayout(2,1))

print(fig1scatSVLdiff +
        theme(legend.position="none") +
        theme(axis.title.y=element_blank(),
              axis.title.x = element_text(size=8, margin = margin(t = 5, r = 0, b = 0, l = 0))),
      vp=vplayout(3,1))

dev.off()


#Supplmentary Figure - species-specific probabilities (males) only first use of each intruder

pdf("Supplement Fig 4 - firsttime.pdf",width=4,height=6)
grid.newpage()
pushViewport(viewport(layout=grid.layout(3,1)))

vplayout<-function(x,y) viewport(layout.pos.row=x,layout.pos.col=y)

print(fig1displayft +
        theme(legend.position="none") +
        theme(axis.title.y=element_blank(),
              axis.title.x = element_text(size=8, margin = margin(t = 5, r = 0, b = 0, l = 0))),
      vp=vplayout(1,1))

print(fig1displayrateft +
        theme(legend.position="none") +
        theme(axis.title.y=element_blank(),
              axis.title.x = element_text(size=8, margin = margin(t = 5, r = 0, b = 0, l = 0))),
      vp=vplayout(2,1))

print(fig1scatft +
        theme(legend.position="none") +
        theme(axis.title.y=element_blank(),
              axis.title.x = element_text(size=8, margin = margin(t = 5, r = 0, b = 0, l = 0))),
      vp=vplayout(3,1))

dev.off()

#Supplementary Figure - sex-specific probabilities (greens) with only first trial for each intruder

pdf("Supplement Fig 5 - firsttimeCAR.pdf",width=4,height=4)
grid.newpage()
pushViewport(viewport(layout=grid.layout(2,1)))

vplayout<-function(x,y) viewport(layout.pos.row=x,layout.pos.col=y)

print(carfig1scatft +
        theme(legend.position="none") +
        theme(axis.title.y=element_blank(),
              axis.title.x = element_text(size=8, margin = margin(t = 5, r = 0, b = 0, l = 0))),
      vp=vplayout(1,1))

print(carfig1upft +
        theme(legend.position="none") +
        theme(axis.title.y=element_blank(),
              axis.title.x = element_text(size=8, margin = margin(t = 5, r = 0, b = 0, l = 0))),
      vp=vplayout(2,1))

dev.off()
