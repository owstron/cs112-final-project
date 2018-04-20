### The following portion of code has been replicated and modified from code provided by Lueders et al
### They used it for the following paper.
### Providing driver’s licenses to unauthorized immigrants in California improves traffic safety
### PNAS
### Hans Lueders, Jens Hainmueller, and Duncan Lawrence
### Replication files
### February 9, 2017


# Preliminaries -----------------------------------------------------------

### Load packages
require(foreign)
require(plm)
require(lmtest)
require(ggplot2)
require(ggthemes)
require(stargazer)
require(Rmisc)
require(dplyr)
require(gridExtra)
require(grid)
require(ggrepel)
require(lmtest)
require(stargazer)

### empty environment
rm(list=ls())

### Load vcovCluster function to compute clustered SEs
### Provided in replication set
source("vcovCluster.R")

### load dataset, which includes variables on:
# Population
# Outstanding Driver's Licenses
# Median Income
# Unemployment
# Vehicle Registrations
# PPIC estimated number of unauthorized immigrants

## Accidents  (compiled from FARS) 
### acc_total: total number of accident by county and year
### hr_total: total number of hit and run accidents by county and year
### acc_killed: total number of people killed by county and year

# See Stata 13 package for data sources.
install.packages("readstata13")

library(readstata13)
d <- read.dta13("Lueders_et al_2017_driverslicenses_replication.dta")
acc.rep <- read.csv("ACC_CALIFORNIA.csv", header=TRUE) #This dataset was created using FARS_data_preparer.R


# Prepare dataset for analysis --------------------------------------------

### generate treatment variable
d$treat <- ifelse(d$year==2015, 1, 0)

d <- subset(d, month == 12)
acc.rep <- subset(ACC_cali, YEAR < 2016)

d$acc_total <- acc.rep$CRASH
d$hr_total <- acc.rep$HR


### generate exposure variables
# predict the number of AB60 licenses issued in 2015
d.dl <- subset(d[,c("county", "year", "dl")], d$month==12 & d$year>=2011 & d$year<=2015)
reg.dl <- lm(dl ~ year*as.factor(county), data=d.dl, 
             d.dl$year>=2011 & d.dl$year<=2014)
d.dl$predict <- predict(reg.dl, newdata=d.dl)
d.dl$resid <- d.dl$dl - d.dl$predict
d.dl$resid <- ifelse(d.dl$resid<0,0,d.dl$resid)

# create thresholds
d.ab60 <- subset(d.dl, d.dl$year==2015)
d.ab60$exposure <- d.ab60$resid / d.ab60$dl

ab60.mean <- mean(d.ab60$exposure)
ab60.median <- median(d.ab60$exposure)
ab60.terciles <- quantile(d.ab60$exposure, probs=c((1/3), 2/3))

# merge continuous exposure measure to dataset
d.ab60 <- d.ab60[,c("county", "exposure", "resid")]
d <- merge(d, d.ab60, by=c("county"), all=T)

# create other exposure measures
d$exposure.mean <- ifelse(d$exposure < ab60.mean, 0, 1)
d$exposure.median <- ifelse(d$exposure < ab60.median, 0, 1)
d$exposure.tercile <- ifelse(d$exposure < ab60.terciles[1], 1,
                             ifelse(d$exposure >=ab60.terciles[1] & d$exposure < ab60.terciles[2], 2,
                                    ifelse(d$exposure >= ab60.terciles[2], 3, NA)))


### re-code year variable. To ease interpretation, 2015 is set equal to the reference year. 
d$year2 <- ifelse(d$year==2015, 2006, d$year)

### Create main dependent variables
### 
# accidents per 1,000 capita
d$acc.pc <- (d$acc_total / d$population) * 1000

# hit and run accidents as share of all accidents (in %)
d$hr.share <- (d$hr_total / d$acc_total) * 100


### supplemntary DVs:

# total accidents that were not hit-and-run
d$acc.nonhr    <- d$acc_total - d$hr_total
d$acc.nonhr.pc <- (d$acc.nonhr / d$population) * 1000

# hit and run accidents per 1,000 capita
d$hr.pc <- (d$hr_total / d$population) * 1000

### Rescale income variable: median income in 1,000 USD
d$income <- d$income / 1000

### Exposure variable should be in percent
d$exposure <- d$exposure * 100

### discard data frames that are not needed
rm(d.ab60, d.dl)

### create the samples of interest
# restrict sample to 2006-2015
d1 <- subset(d, d$year >= 2006 & d$year<=2015)

# There is one observation with no accidents. Consequently, you cannot define % fatal or HR accidents:
# d2 <- d1[!is.na(d1$hr.share),]
d2 <- na.omit(d1)

# Replication of Figure 2 ----------------------------------------------------------------

d6 <- d1
d6$after <- 0 
d6$after[d6$year==2015] <- 1

d6 <- na.omit(d6[,c("acc.pc", "hr.share","population","exposure",
                    "county","after")])

d6 <- aggregate(d6[,c("acc.pc", "hr.share","population","exposure") ],
                by=as.list(d6[,c("after", "county")]), median, na.rm=T)

d6 <- reshape(d6, timevar = "after", 
              idvar = "county", 
              direction = "wide",
              v.names = c("acc.pc","hr.share","population","exposure"))


d6$acc.pc.delta <- (d6$acc.pc.1 - d6$acc.pc.0) 
d6$hr.share.delta <- (d6$hr.share.1 - d6$hr.share.0) 

# convert the data into long format to facilitate faceting
d6 <- d6[,c("county", "exposure.1", "acc.pc.delta", "hr.share.delta")]
colnames(d6) <- c("county", "exposure", "delta1", "delta2")

d6 <- reshape(d6, 
              varying = c("delta1", "delta2"), 
              v.names = "delta",
              timevar = "outcome", 
              times = c("1", "2"), 
              direction = "long")

d6$outcome <- factor(d6$outcome, 
                     levels=c(1,2),
                     labels=c("(a) Accidents per 1,000 capita",
                              "(b) % Hit-and-Run Accidents"))

p <- ggplot(d6, aes(x= exposure, y= delta, label=county)) +
  geom_point(size=2, col="gray50") + 
  geom_text_repel(aes(label=county),size=6) +
  stat_smooth(method="lm", se=T,col="firebrick3") +
  facet_wrap(~ outcome, ncol=1, scales="free") +  
  xlab("Share of AB60 licenses") + 
  ylab("Change between pre- and post-treatment period") + 
  theme(axis.title = element_text(size=22, face="bold"),
        axis.text=element_text(size=20),
        strip.text.x = element_text(size = 26, face ="bold"))

ggsave(p, file="changes_pre_post.pdf", width=13, height=29) 


# Main Models: Table 1 ------------------------------------

### Regressions
# a. accidents per capita x terciles
m1 <- lm(acc.pc ~ treat * as.factor(exposure.tercile) + as.factor(year2) +  
           as.factor(county), data=d1)
c1 <- coeftest(m1, vcov = vcovCluster(m1, factor(d1$county)))

# b. accidents per capita x continuous
m2 <- lm(acc.pc ~ treat * exposure + as.factor(year2) + as.factor(county), 
         data=d1)
c2 <- coeftest(m2, vcov = vcovCluster(m2, factor(d1$county)))

# e. % hit-and-run accidents x terciles
m3 <- lm(hr.share ~ treat * as.factor(exposure.tercile) + as.factor(year2) +  
           as.factor(county), data=d2)
c3 <- coeftest(m3, vcov = vcovCluster(m3, factor(d2$county)))

# f. % hit-and-run accidents x continuous
m4 <- lm(hr.share ~ treat * exposure + as.factor(year2) +  as.factor(county), 
         data=d2)
c4 <- coeftest(m4, vcov = vcovCluster(m4, factor(d2$county)))


### output
stargazer(m1, m2, m3, m4,
          se=list(c1[,2], c2[,2], 
                  c3[,2], c4[,2]),
          omit.stat = c("rsq", "f", "ser"), df=F, notes.align="c", no.space=T,
          font.size="scriptsize",
          model.names = F,
          omit=c("year", "county", "Constant"),
          dep.var.caption = "",
          dep.var.labels = c("Accidents per 1,000 capita", "Share Hit-and-Run Accidents"),
          star.cutoffs = c(0.05),
          add.lines=list(c("Implied percent change in outcome compared to 2014", 
                           paste(round(m1$coefficients[71] / mean(d$acc.pc[d$year==2014], na.rm=T) * 100, digits=1),"%",sep=""),
                           paste(round(m2$coefficients[69] / mean(d$acc.pc[d$year==2014], na.rm=T) * 100 * mean(d$exposure, na.rm=T), 
                                       digits=1),"%",sep=""),
                           paste(round(m5$coefficients[71] / mean(d$hr.share[d$year==2014], na.rm=T) * 100, digits=1), "%", sep=""),
                           paste(round(m6$coefficients[69] / mean(d$hr.share[d$year==2014], na.rm=T) * 100 * mean(d$exposure, na.rm=T),
                                       digits=1), "%", sep="")),
                         c("",
                           paste("(",round(c1[69,2] / mean(d$acc.pc[d$year==2014], na.rm=T) * 100, digits=1),"%)",sep=""),
                           paste("(",round(c2[68,2] / mean(d$acc.pc[d$year==2014], na.rm=T) * 100 * mean(d$exposure, na.rm=T), 
                                           digits=1),"%)",sepdi=""),
                           paste("(",round(c3[69,2] / mean(d$hr.share[d$year==2014], na.rm=T) * 100, digits=1),"%)",sep=""),
                           paste("(",round(c4[68,2] / mean(d$hr.share[d$year==2014], na.rm=T) * 100 * mean(d$exposure, na.rm=T), 
                                           digits=1),"%)",sep=""))))

### mean and SDs for DVs
# accidents per capita
paste("Mean = ", round(mean(d1$acc.pc, na.rm=T), digits=3), 
      ", SD = ",round(sd(d1$acc.pc, na.rm=T), digits=3), sep="")

# % hit-and-run accidents
paste("Mean = ", round(mean(d2$hr.share, na.rm=T), digits=3), 
      ", SD = ",round(sd(d2$hr.share, na.rm=T), digits=3), sep="")



