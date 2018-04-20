### Providing driver’s licenses to unauthorized immigrants in California improves traffic safety
### PNAS
### Hans Lueders, Jens Hainmueller, and Duncan Lawrence
### Replication files
### February 9, 2017

### Note:
### This script contains the code for our analysis in the paper and SM
### The data on accidents is controlled by the California Highway Patrol and we legally not allowed to post it. 
### Therefore we can only share the limited data that excludes these variables

###########################################################################
######## FOR PAPER ########################################################
###########################################################################

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

### empty environment
rm(list=ls())

### Load vcovCluster function to compute clustered SEs
source("vcovCluster.R")

### load dataset, which includes variables on:
# Population
# Outstanding Driver's Licenses
# Median Income
# Unemployment
# Vehicle Registrations
# PPIC estimated number of unauthorized immigrants

## Accidents  (posting restricted by CHP)
### acc_total: total number of accident by county and month
### acc_fatal: total number of fatal accidents by county and month
### hr_total: total number of hit and run accidents by county and month
### acc_killed: total number of people killed by county and month

# See SM for data sources.
d <- read.dta("Lueders_et al_2017_driverslicenses_replication.dta")

# Prepare dataset for analysis --------------------------------------------

### generate treatment variable
d$treat <- ifelse(d$year==2015, 1, 0)

### generate exposure variables
# predict the number of AB60 licenses issued in 2015
d.dl <- subset(d[,c("county", "year", "dl")], d$month==1 & d$year>=2011 & d$year<=2015)
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

# share of all accidents that were fatal
d$fatal.share <- (d$acc_fatal / d$acc_total) * 100

# hit and run accidents as share of all accidents (in %)
d$hr.share <- (d$hr_total / d$acc_total) * 100


### supplemntary DVs:

# total accidents that were not hit-and-run
d$acc.nonhr    <- d$acc_total - d$hr_total
d$acc.nonhr.pc <- (d$acc.nonhr / d$population) * 1000

# fatal accidents per 1,000 capita
d$acc.fatal.pc <- (d$acc_fatal / d$population) * 1000

# people killed in accidents per 1,000 capita
d$acc.killed.pc <- (d$acc_killed / d$population) * 1000

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
d2 <- d1[!is.na(d1$hr.share),]

# Main Models: Table 1 ------------------------------------

### Regressions
# a. accidents per capita x terciles
m1 <- lm(acc.pc ~ treat * as.factor(exposure.tercile) + as.factor(year2) + as.factor(month) + 
           as.factor(county), data=d1)
c1 <- coeftest(m1, vcov = vcovCluster(m1, factor(d1$county)))

# b. accidents per capita x continuous
m2 <- lm(acc.pc ~ treat * exposure + as.factor(year2) + as.factor(month) + as.factor(county), 
         data=d1)
c2 <- coeftest(m2, vcov = vcovCluster(m2, factor(d1$county)))

# c. % fatal accidents x terciles
m3 <- lm(fatal.share ~ treat * as.factor(exposure.tercile) + as.factor(year2) + 
           as.factor(month) + as.factor(county), data=d2)
c3 <- coeftest(m3, vcov = vcovCluster(m3, factor(d2$county)))

# d. % fatal accidents x continuous
m4 <- lm(fatal.share ~ treat * exposure + as.factor(year2) + as.factor(month) + as.factor(county), 
         data=d2)
c4 <- coeftest(m4, vcov = vcovCluster(m4, factor(d2$county)))

# e. % hit-and-run accidents x terciles
m5 <- lm(hr.share ~ treat * as.factor(exposure.tercile) + as.factor(year2) + as.factor(month) + 
           as.factor(county), data=d2)
c5 <- coeftest(m5, vcov = vcovCluster(m5, factor(d2$county)))

# f. % hit-and-run accidents x continuous
m6 <- lm(hr.share ~ treat * exposure + as.factor(year2) + as.factor(month) + as.factor(county), 
         data=d2)
c6 <- coeftest(m6, vcov = vcovCluster(m6, factor(d2$county)))


### output
stargazer(m1, m2, m3, m4, m5, m6,
          se=list(c1[,2], c2[,2], c3[,2], c4[,2], 
                  c5[,2], c6[,2]),
          omit.stat = c("rsq", "f", "ser"), df=F, notes.align="c", no.space=T,
          font.size="scriptsize",
          model.names = F,
          omit=c("month", "year", "county", "Constant"),
          dep.var.caption = "",
          dep.var.labels = c("Accidents per 1,000 capita", "Share Fatal Accidents", "Share Hit-and-Run Accidents"),
          star.cutoffs = c(0.05),
          add.lines=list(c("Implied percent change in outcome compared to 2014", 
                           paste(round(m1$coefficients[82] / mean(d$acc.pc[d$year==2014], na.rm=T) * 100, digits=1),"%",sep=""),
                           paste(round(m2$coefficients[80] / mean(d$acc.pc[d$year==2014], na.rm=T) * 100 * mean(d$exposure, na.rm=T), 
                                       digits=1),"%",sep=""),
                           paste(round(m3$coefficients[82] / mean(d$fatal.share[d$year==2014], na.rm=T) * 100, digits=1),"%",sep=""),
                           paste(round(m4$coefficients[80] / mean(d$fatal.share[d$year==2014], na.rm=T) * 100 * mean(d$exposure, na.rm=T), 
                                       digits=1),"%",sep=""),
                           paste(round(m5$coefficients[82] / mean(d$hr.share[d$year==2014], na.rm=T) * 100, digits=1), "%", sep=""),
                           paste(round(m6$coefficients[80] / mean(d$hr.share[d$year==2014], na.rm=T) * 100 * mean(d$exposure, na.rm=T),
                                       digits=1), "%", sep="")),
                         c("",
                           paste("(",round(c1[80,2] / mean(d$acc.pc[d$year==2014], na.rm=T) * 100, digits=1),"%)",sep=""),
                           paste("(",round(c2[79,2] / mean(d$acc.pc[d$year==2014], na.rm=T) * 100 * mean(d$exposure, na.rm=T), 
                                           digits=1),"%)",sep=""),
                           paste("(",round(c3[80,2] / mean(d$fatal.share[d$year==2014], na.rm=T) * 100, digits=1),"%)",sep=""),
                           paste("(",round(c4[79,2] / mean(d$fatal.share[d$year==2014], na.rm=T) * 100 * mean(d$exposure, na.rm=T), 
                                           digits=1),"%)",sep=""),
                           paste("(",round(c5[80,2] / mean(d$hr.share[d$year==2014], na.rm=T) * 100, digits=1),"%)",sep=""),
                           paste("(",round(c6[79,2] / mean(d$hr.share[d$year==2014], na.rm=T) * 100 * mean(d$exposure, na.rm=T), 
                                           digits=1),"%)",sep=""))))

### mean and SDs for DVs
# accidents per capita
paste("Mean = ", round(mean(d1$acc.pc, na.rm=T), digits=3), 
      ", SD = ",round(sd(d1$acc.pc, na.rm=T), digits=3), sep="")

# % fatal accidents
paste("Mean = ", round(mean(d2$fatal.share, na.rm=T), digits=3), 
      ", SD = ",round(sd(d2$fatal.share, na.rm=T), digits=3), sep="")

# % hit-and-run accidents
paste("Mean = ", round(mean(d2$hr.share, na.rm=T), digits=3), 
      ", SD = ",round(sd(d2$hr.share, na.rm=T), digits=3), sep="")




# Figure 1 ----------------------------------------------------------------

# load required packages
library(maptools)
library(mapproj)
library(rgeos)
library(rgdal)
library(ggmap)
library(ggalt)

# select data of interest
plottingdata <- subset(d, d$year==2010 & d$month==1)[,c("county", "year")]
plottingdata$id <- plottingdata$county

# if required: install albersusa package, which requires the "sf" and "devtools" packages
# library(devtools)
# devtools::install_github("hrbrmstr/albersusa")
library(albersusa)

# plot(usa_composite(proj="laea"))

us <- usa_composite()
us_map <- fortify(us, region="name")

# coding the policy
us_map$policy <- ifelse(us_map$id=="California" | us_map$id=="Nevada" | us_map$id=="Utah" | us_map$id=="Colorado" |
                          us_map$id=="Illinois" | us_map$id=="Vermont" | us_map$id=="Connecticut" | us_map$id=="Maryland" | 
                          us_map$id=="District of Columbia" | us_map$id=="Delaware" | us_map$id=="Hawaii", 1, 
                        ifelse(us_map$id=="Washington" | us_map$id == "New Mexico", 2, 0))


# state abbreviations
us_map$id2 <- NA
states <- unique(us_map$id)
stateid <- c("AL", 	"AK", 	"AZ", 	"AR", 	"CA", 	
             "CO", 	"CT", 	"DE", 	"DC", 	"FL", 	"GA", 	
             "HI", 	"ID", 	"IL", 	"IN", 	"IA", 	"KS", 
             "KY", 	"LA", 	"ME", 	"MD", 	"MA", 	"MI", 	
             "MN", 	"MS", 	"MO", 	"MT", 	"NE", 	"NV", 	
             "NH", 	"NJ", 	"NM", 	"NY", 	"NC", 	"ND", 	
             "OH", 	"OK", 	"OR", 	"PA", 	"RI", 	"SC", 	
             "SD", 	"TN", 	"TX", 	"UT", 	"VT", 	"VA", 	
             "WA", 	"WV", 	"WI", 	"WY")

for(i in 1:length(states)){
  us_map$id2[us_map$id==states[i]] <- stateid[i]
}

# for state labels
snames <- aggregate(cbind(long, lat) ~ id2, data=us_map, 
                    FUN=function(x)mean(range(x)))

# recode the position of state labels to center them better
snames$long[snames$id2=="AK"] <- -112.147686
snames$lat[snames$id2=="AK"] <- 27.503678
snames$long[snames$id2=="CA"] <- -120.27054
snames$lat[snames$id2=="CA"] <- 37.27184
snames$long[snames$id2=="DE"] <- -75.548163
snames$lat[snames$id2=="DE"] <- 38.620091
snames$long[snames$id2=="FL"] <- -81.417685
snames$lat[snames$id2=="FL"] <- 27.889527
snames$long[snames$id2=="ID"] <- -114.500689
snames$lat[snames$id2=="ID"] <- 43.572862
snames$long[snames$id2=="IL"] <- -89.20409
snames$lat[snames$id2=="IL"] <- 39.73930
snames$long[snames$id2=="KY"] <- -85.36990
snames$lat[snames$id2=="KY"] <- 37.42229
snames$long[snames$id2=="LA"] <- -91.63008
snames$lat[snames$id2=="LA"] <- 30.47219
snames$long[snames$id2=="MA"] <- -71.732516
snames$lat[snames$id2=="MA"] <- 42.318379
snames$long[snames$id2=="MD"] <- -76.966901
snames$lat[snames$id2=="MD"] <- 39.457513
snames$long[snames$id2=="MI"] <- -84.606075
snames$lat[snames$id2=="MI"] <- 43.383979
snames$long[snames$id2=="MN"] <- -94.695559
snames$lat[snames$id2=="MN"] <- 45.984195
snames$long[snames$id2=="NC"] <- -79.29026
snames$lat[snames$id2=="NC"] <- 35.51632
snames$long[snames$id2=="NH"] <- -71.543811
snames$lat[snames$id2=="NH"] <- 43.201985
snames$long[snames$id2=="NV"] <- -116.544210
snames$lat[snames$id2=="NV"] <- 39.550948
snames$long[snames$id2=="OK"] <- -97.305605
snames$lat[snames$id2=="OK"] <- 35.437769
snames$long[snames$id2=="TX"] <- -98.989946
snames$lat[snames$id2=="TX"] <- 31.466855
snames$long[snames$id2=="VA"] <- -78.472453
snames$lat[snames$id2=="VA"] <- 37.303534
snames$long[snames$id2=="VT"] <- -72.497788
snames$lat[snames$id2=="VT"] <- 44.520066
snames$long[snames$id2=="WV"] <- -81.227773
snames$lat[snames$id2=="WV"] <- 38.297669


# Plot
map.license.overview <- ggplot() + 
  geom_map(data=us_map, map=us_map,
           aes(x=long, y=lat, map_id=id), fill=ifelse(us_map$policy==1 | us_map$policy==2, "firebrick3", "gray80"),
           color="gray40", size=0.25) + 
  geom_text(data=snames, aes(long, lat, label = id2), size=4, col="black")  +
  theme_nothing(legend=T) +
  theme(legend.text = element_text(size = 17, face="bold"), legend.position="bottom") +
  theme(plot.title=element_text(size = 25, face="bold"))
ggsave(map.license.overview, file="map_license_overview.pdf", width=13, height=10)





# Figure 2 ----------------------------------------------------------------

d6 <- d1
d6$after <- 0 
d6$after[d6$year==2015] <- 1

d6 <- na.omit(d6[,c("acc.pc", "fatal.share","hr.share","population","exposure",
                    "county","after")])

d6 <- aggregate(d6[,c("acc.pc", "fatal.share","hr.share","population","exposure") ],
                by=as.list(d6[,c("after", "county")]), mean, na.rm=T)

d6 <- reshape(d6, timevar = "after", 
              idvar = "county", 
              direction = "wide",
              v.names = c("acc.pc","fatal.share","hr.share","population","exposure"))


d6$acc.pc.delta <- (d6$acc.pc.1 - d6$acc.pc.0) 
d6$fatal.share.delta <- (d6$fatal.share.1 - d6$fatal.share.0) 
d6$hr.share.delta <- (d6$hr.share.1 - d6$hr.share.0) 

# convert the data into long format to facilitate faceting
d6 <- d6[,c("county", "exposure.1", "acc.pc.delta", "fatal.share.delta", "hr.share.delta")]
colnames(d6) <- c("county", "exposure", "delta1", "delta2", "delta3")

d6 <- reshape(d6, 
              varying = c("delta1", "delta2", "delta3"), 
              v.names = "delta",
              timevar = "outcome", 
              times = c("1", "2", "3"), 
              direction = "long")

d6$outcome <- factor(d6$outcome, 
                   levels=c(1,2,3),
                   labels=c("(a) Accidents per 1,000 capita",
                            "(b) % Fatal Accidents",
                            "(c) % Hit-and-Run Accidents"))

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




# Figure 3 ----------------------------------------------------------------

# Computing increases
dta <- subset(d, d$month==1 & d$year<=2015)[,c("county", "year", "exposure.median", "autos", "dl")]

dta <- ddply(dta, .(county), transform, autos.lag = lag(autos))
dta <- ddply(dta, .(county), transform, dl.lag = lag(dl))

dta$autos.delta <- ((dta$autos - dta$autos.lag) / dta$autos.lag)*100
dta$dl.delta <- ((dta$dl - dta$dl.lag) / dta$dl.lag)*100

aux <- seq(2010, 2015, 1)

dta$exposure.median.factor <- factor(dta$exposure.median, 
                                     levels=c(0,1),
                                     labels=c("AB60 Low", "AB60 High"))

change <- ggplot() + 
  geom_point(data=subset(dta, dta$year<2015 & dta$county != "Alpine" & 
                           dta$county != "Mono" & dta$county != "Sierra"), aes(x=dl.delta, y=autos.delta), 
             col="gray30", size=1) + 
  geom_point(data=subset(dta, dta$year==2015 & dta$county != "Alpine" & 
                           dta$county != "Mono" & dta$county != "Sierra"), aes(x=dl.delta, y=autos.delta), 
             col="firebrick3", size=1) +   
  geom_abline(intercept=0, slope=1, col="grey70", lty=2) +
  stat_smooth(data=subset(dta, dta$year==2010 & dta$county != "Alpine" & 
                            dta$county != "Mono" & dta$county != "Sierra"), aes(x=dl.delta, y=autos.delta), 
              method="lm", se=F,col="gray10") +
  stat_smooth(data=subset(dta, dta$year==2011 & dta$county != "Alpine" & 
                            dta$county != "Mono" & dta$county != "Sierra"), aes(x=dl.delta, y=autos.delta), 
              method="lm", se=F,col="gray10") +
  stat_smooth(data=subset(dta, dta$year==2012 & dta$county != "Alpine" & 
                            dta$county != "Mono" & dta$county != "Sierra"), aes(x=dl.delta, y=autos.delta), 
              method="lm", se=F,col="gray10") +
  stat_smooth(data=subset(dta, dta$year==2013 & dta$county != "Alpine" & 
                            dta$county != "Mono" & dta$county != "Sierra"), aes(x=dl.delta, y=autos.delta), 
              method="lm", se=F,col="gray10") +  
  stat_smooth(data=subset(dta, dta$year==2014 & dta$county != "Alpine" & 
                            dta$county != "Mono" & dta$county != "Sierra"), aes(x=dl.delta, y=autos.delta), 
              method="lm", se=F,col="gray10") +
  stat_smooth(data=subset(dta, dta$year==2015 & dta$county != "Alpine" & 
                            dta$county != "Mono" & dta$county != "Sierra"), aes(x=dl.delta, y=autos.delta), 
              method="lm", se=F,col="firebrick3") +  
  facet_grid(. ~ exposure.median.factor) +
  theme(strip.text.x = element_text(size = 22, face ="bold")) +
  xlab(expression(bold(Delta~Licenses~(percent)))) + 
  ylab(expression(bold(Delta~Registered~Autos~(percent)))) + 
  theme(axis.title = element_text(size=22, face="bold"),
        axis.text=element_text(size=20))
ggsave(change, file="change_autos_dl1.pdf", width=15, height=10)  






###########################################################################
######## FOR SM ###########################################################
###########################################################################



# Summary Statistics ------------------------------------------------------

summ.stats <- matrix(ncol=7, nrow=10)
colnames(summ.stats) <- c("Mean", "SD", "Min", "1st Quartile", "Median", "3rd Quartile", "Max")
rownames(summ.stats) <- c("Accidents p.c.", "Share Fatal Accidents", "Share Hit-and-Run Accidents", 
                          "Pct. AB60 Licenses", "Median Income", "Unemployment Rate","Hit-and-Run Accidents p.c.",
                          "Fatal Accidents p.c.", "People Killed in Accidents p.c.", "Non-Hit-and-Run Accidents p.c.")
summ.stats[1,] <- c(mean(d1$acc.pc, na.rm=T), sd(d1$acc.pc, na.rm=T), fivenum(d1$acc.pc))
summ.stats[2,] <- c(mean(d1$fatal.share, na.rm=T), sd(d1$fatal.share, na.rm=T), fivenum(d1$fatal.share))
summ.stats[3,] <- c(mean(d1$hr.share, na.rm=T), sd(d1$hr.share, na.rm=T), fivenum(d1$hr.share))
summ.stats[4,] <- c(mean(d1$exposure, na.rm=T), sd(d1$exposure, na.rm=T), fivenum(d1$exposure))
summ.stats[5,] <- c(mean(d1$income, na.rm=T), sd(d1$income, na.rm=T), fivenum(d1$income))
summ.stats[6,] <- c(mean(d1$unemploymentrate, na.rm=T), sd(d1$unemploymentrate, na.rm=T), fivenum(d1$unemploymentrate))
summ.stats[7,] <- c(mean(d1$hr.pc, na.rm=T), sd(d1$hr.pc, na.rm=T), fivenum(d1$hr.pc))
summ.stats[8,] <- c(mean(d1$acc.fatal.pc, na.rm=T), sd(d1$acc.fatal.pc, na.rm=T), fivenum(d1$acc.fatal.pc))
summ.stats[9,] <- c(mean(d1$acc.killed.pc, na.rm=T), sd(d1$acc.killed.pc, na.rm=T), fivenum(d1$acc.killed.pc))
summ.stats[10,] <- c(mean(d1$acc.nonhr.pc, na.rm=T), sd(d1$acc.nonhr.pc, na.rm=T), fivenum(d1$acc.nonhr.pc))
xtable(summ.stats, digits=2)



# Boxplot: yearly changes in driver's licenses ----------------------------

box <- ggplot() + 
  geom_boxplot(data=dta[dta$year>=2010 & dta$year<=2014,], aes(factor(year),dl.delta), color="black") + 
  geom_boxplot(data=dta[dta$year==2015,], aes(factor(year),dl.delta), color="firebrick2") +  
  facet_grid(. ~ exposure.median.factor) +
  xlab("Year") + 
  ylab(expression(bold(Delta~Licenses~(percent)))) +
  theme(axis.title = element_text(size=22, face="bold"),
        axis.text=element_text(size=20)) +
  theme(strip.text.x = element_text(size = 18, face ="bold"))
ggsave(box, file="changes_licenses_boxplot.pdf", width=15, height=10)  




# Illustration of AB60 estimation strategy --------------------------------

d.license <- subset(d, d$month==1 & d$year>=2008 & d$year<=2015)
d.license <- d.license[,c("county", "year","dl","exposure", "population")]
d.license$ab60 <- ifelse(d.license$year<2015, 0, d.license$exposure/100*d.license$dl)

license.reg <- lm(dl ~ year*as.factor(county), data=d.license, year>=2011 & year<=2014)

d.license$prediction <- predict(license.reg, newdata=d.license)

d.license.state <- aggregate(d.license$dl, by=list(d.license$year), FUN=sum)
d.license.state2 <- aggregate(d.license$prediction, by=list(d.license$year), FUN=sum)

d.license.state$predict <- d.license.state2$x
colnames(d.license.state) <- c("year", "actual", "predict")

temp <- paste("Delta ==", 559200)

data.segment <- data.frame(x1 = 2015, x2=2015, y1=d.license.state[8,3], y2=d.license.state[8,2])


modelab60 <- ggplot() +
  geom_line(data=d.license.state, 
            aes(x=year, y=actual, col="Actual number"), size=1.1)  +
  geom_line(data=d.license.state[d.license.state$year>=2011,], 
            aes(x=year, y=predict, col="Prediction"), lty=2, size=1.1)  +
  ylab("Number of Licenses") + 
  xlab("Year") + 
  scale_x_continuous(lim=c(2008,2016), 
                     breaks=c(2008,2009,2010,2011,2012,2013,2014,2015),
                     labels=c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015")) +  
  scale_y_continuous(lim=c(2.25e+07, 2.6e+07),
                     breaks=c(2.3e+07, 2.4e+07, 2.5e+07, 2.6e+07),
                     labels=c("23 million", "24 million", "25 million", "26 million")) +
  theme(axis.title = element_text(size=22, face="bold"),
        axis.text=element_text(size=20)) +
  theme(plot.title = element_text(size=20, face="bold", colour="black",
                                  hjust=.5)) + 
  scale_color_manual(values = c("Actual number" = "steelblue",
                                "Prediction" = "firebrick"),
                     name="") + 
  theme(legend.title = element_text(face="bold", size=22),
        legend.text = element_text(size = 20), legend.position="bottom") +
  guides(col = guide_legend(nrow=1, override.aes = list(linetype=c(1,2))))  +
  geom_segment(aes(x=x1, xend=x2, y=y1, yend=y2), data=data.segment, lty=3, col="black", size=1.1) +
  annotate("text", x=2015.6, y=(data.segment[,4] - (data.segment[,4] - data.segment[,3])/2), parse=T,
           label=temp, size=7)
ggsave(modelab60, file="modelab60.pdf", width=15, height=12) 



# Estimating the counterfactual reduction in Number of Hit and Run Accidents --------

# everything is based on this model:
m6 <- lm(hr.share ~ treat * exposure + as.factor(year2) + as.factor(month) + 
           as.factor(county), data=d2)
c6 <- coeftest(m6, vcov = vcovCluster(m6, factor(d2$county)))

# extract coefficient on interaction term
beta <- c6[79,1]

# create empty vector to store each county's prediction
sums <- vector()

# loop through all 58 counties to compute each county's predicted 
# number of accidents that were not hit-and-run
for(i in 1:length(unique(d$county))) { 
  sums[i] <- (beta/100) * d$exposure[d$year==2015 & d$month==1 & d$county==unique(d$county)[i]] * 
    sum(d$acc_total[d$year==2015 & d$county == unique(d$county)[i]], na.rm=T)
}

# sum up all 58 values to compute the total number of hit-and-runs that did not happen
sum(sums) 

# To compare: hit and run accidents in 2014:
sum(d$hr_total[d$year==2014], na.rm=T)

# To compare: hit and run accidents in 2015:
sum(d$hr_total[d$year==2015], na.rm=T)



# Year-to-year increases auto registration --------------------------------

### computing increases
d.vehicles <- subset(d, d$month==1 & d$year<=2015)[,c("county", "year", "exposure.median", "autos", "dl")]

d.vehicles <- ddply(d.vehicles, .(county), transform, autos.lag = lag(autos))
d.vehicles <- ddply(d.vehicles, .(county), transform, dl.lag = lag(dl))

d.vehicles$autos.delta <- ((d.vehicles$autos - d.vehicles$autos.lag) / d.vehicles$autos.lag)*100
d.vehicles$dl.delta <- ((d.vehicles$dl - d.vehicles$dl.lag) / d.vehicles$dl.lag)*100


### Regressions
# 2010
increase2010 <- lm(autos.delta ~ dl.delta, data=d.vehicles, year==2010 & 
                     county != "Alpine" & county != "Mono" & county != "Sierra")
increase2010.coef <- coeftest(increase2010, vcov=vcovHC(increase2010, type="HC2"))

# 2011
increase2011 <- lm(autos.delta ~ dl.delta, data=d.vehicles, year==2011 & 
                     county != "Alpine" & county != "Mono" & county != "Sierra")
increase2011.coef <- coeftest(increase2011, vcov=vcovHC(increase2011, type="HC2"))

# 2012
increase2012 <- lm(autos.delta ~ dl.delta, data=d.vehicles, year==2012 & 
                     county != "Alpine" & county != "Mono" & county != "Sierra")
increase2012.coef <- coeftest(increase2012, vcov=vcovHC(increase2012, type="HC2"))

# 2013
increase2013 <- lm(autos.delta ~ dl.delta, data=d.vehicles, year==2013 & 
                     county != "Alpine" & county != "Mono" & county != "Sierra")
increase2013.coef <- coeftest(increase2013, vcov=vcovHC(increase2013, type="HC2"))

# 2014
increase2014 <- lm(autos.delta ~ dl.delta, data=d.vehicles, year==2014 & 
                     county != "Alpine" & county != "Mono" & county != "Sierra")
increase2014.coef <- coeftest(increase2014, vcov=vcovHC(increase2014, type="HC2"))

# 2015
increase2015 <- lm(autos.delta ~ dl.delta, data=d.vehicles, year==2015 & 
                     county != "Alpine" & county != "Mono" & county != "Sierra")
increase2015.coef <- coeftest(increase2015, vcov=vcovHC(increase2015, type="HC2"))

# 2015: low
increase2015l <- lm(autos.delta ~ dl.delta, data=d.vehicles, year==2015 & 
                     county != "Alpine" & county != "Mono" & county != "Sierra" & 
                      exposure.median==0)
increase2015l.coef <- coeftest(increase2015l, vcov=vcovHC(increase2015l, type="HC2"))

# 2015: high
increase2015h <- lm(autos.delta ~ dl.delta, data=d.vehicles, year==2015 & 
                      county != "Alpine" & county != "Mono" & county != "Sierra" & 
                      exposure.median==1)
increase2015h.coef <- coeftest(increase2015h, vcov=vcovHC(increase2015h, type="HC2"))


### output
stargazer(increase2010, increase2011, increase2012, increase2013, increase2014, increase2015, increase2015l, increase2015h,
          se=list(increase2010.coef[,2], increase2011.coef[,2], increase2012.coef[,2], increase2013.coef[,2], 
                  increase2014.coef[,2], increase2015.coef[,2], increase2015l.coef[,2], increase2015h.coef[,2]),
          omit.stat = c("rsq", "f", "ser"), df=F, notes.align="c", no.space=T,
          font.size="scriptsize",
          star.cutoffs = c(0.05),        
          model.names = F,
          dep.var.caption = "",
          dep.var.labels = "")



# Weighted regressions ------------------------------------------------

### need to code country's average population--will be weight
d.pop <- as.data.frame(matrix(nrow=length(unique(d$county)), ncol=2))
d.pop[,1] <- unique(d$county)
colnames(d.pop) <- c("county", "pop.mean")

for(i in 1:length(unique(d$county))){
  d.pop[i,2] <- mean(d$population[d$month==1 & d$county==unique(d$county)[i]], na.rm=T)
}

d <- merge(d, d.pop, by="county")
d1 <- merge(d1, d.pop, by="county")
d2 <- merge(d2, d.pop, by="county")


### Regressions
# a. accidents per capita x terciles
m1 <- lm(acc.pc ~ treat * as.factor(exposure.tercile) + as.factor(year2) + as.factor(month) + 
           as.factor(county), data=d1, weights = log(pop.mean))
c1 <- coeftest(m1, vcov = vcovCluster(m1, factor(d1$county)))

# b. accidents per capita x continuous
m2 <- lm(acc.pc ~ treat * exposure + as.factor(year2) + as.factor(month) + as.factor(county), 
         data=d1, weights = log(pop.mean))
c2 <- coeftest(m2, vcov = vcovCluster(m2, factor(d1$county)))

# c. % fatal.accidents x terciles
m3 <- lm(fatal.share ~ treat * as.factor(exposure.tercile) + as.factor(year2) + 
           as.factor(month) + as.factor(county), data=d2, weights = log(pop.mean))
c3 <- coeftest(m3, vcov = vcovCluster(m3, factor(d2$county)))

# d. % fatal.accidents x continuous
m4 <- lm(fatal.share ~ treat * exposure + as.factor(year2) + as.factor(month) + as.factor(county), 
         data=d2, weights = log(pop.mean))
c4 <- coeftest(m4, vcov = vcovCluster(m4, factor(d2$county)))

# e. % hit-and-run accidents x terciles
m5 <- lm(hr.share ~ treat * as.factor(exposure.tercile) + as.factor(year2) + as.factor(month) + 
           as.factor(county), data=d2, weights = log(pop.mean))
c5 <- coeftest(m5, vcov = vcovCluster(m5, factor(d2$county)))

# f. % hit-and-run accidents x continuous
m6 <- lm(hr.share ~ treat * exposure + as.factor(year2) + as.factor(month) + as.factor(county), 
         data=d2, weights = log(pop.mean))
c6 <- coeftest(m6, vcov = vcovCluster(m6, factor(d2$county)))

### output
stargazer(m1, m2, m3, m4, m5, m6,
          se=list(c1[,2], c2[,2], c3[,2], c4[,2], 
                  c5[,2], c6[,2]),
          omit.stat = c("rsq", "f", "ser"), df=F, notes.align="c", no.space=T,
          font.size="scriptsize",
          model.names = F,
          omit=c("month", "year", "county", "Constant"),
          dep.var.caption = "",
          dep.var.labels = c("Accidents per 1,000 capita", "Share Fatal Accidents", "Share Hit-and-Run Accidents"),
          star.cutoffs = c(0.05),
          add.lines=list(c("Implied percent change in outcome compared to 2014", 
                           paste(round(m1$coefficients[82] / mean(d$acc.pc[d$year==2014], na.rm=T) * 100, digits=1),"%",sep=""),
                           paste(round(m2$coefficients[80] / mean(d$acc.pc[d$year==2014], na.rm=T) * 100 * mean(d$exposure, na.rm=T), 
                                       digits=1),"%",sep=""),
                           paste(round(m3$coefficients[82] / mean(d$fatal.share[d$year==2014], na.rm=T) * 100, digits=1),"%",sep=""),
                           paste(round(m4$coefficients[80] / mean(d$fatal.share[d$year==2014], na.rm=T) * 100 * mean(d$exposure, na.rm=T), 
                                       digits=1),"%",sep=""),
                           paste(round(m5$coefficients[82] / mean(d$hr.share[d$year==2014], na.rm=T) * 100, digits=1), "%", sep=""),
                           paste(round(m6$coefficients[80] / mean(d$hr.share[d$year==2014], na.rm=T) * 100 * mean(d$exposure, na.rm=T),
                                       digits=1), "%", sep="")),
                         c("",
                           paste("(",round(c1[80,2] / mean(d$acc.pc[d$year==2014], na.rm=T) * 100, digits=1),"%)",sep=""),
                           paste("(",round(c2[79,2] / mean(d$acc.pc[d$year==2014], na.rm=T) * 100 * mean(d$exposure, na.rm=T), 
                                           digits=1),"%)",sep=""),
                           paste("(",round(c3[80,2] / mean(d$fatal.share[d$year==2014], na.rm=T) * 100, digits=1),"%)",sep=""),
                           paste("(",round(c4[79,2] / mean(d$fatal.share[d$year==2014], na.rm=T) * 100 * mean(d$exposure, na.rm=T), 
                                           digits=1),"%)",sep=""),
                           paste("(",round(c5[80,2] / mean(d$hr.share[d$year==2014], na.rm=T) * 100, digits=1),"%)",sep=""),
                           paste("(",round(c6[79,2] / mean(d$hr.share[d$year==2014], na.rm=T) * 100 * mean(d$exposure, na.rm=T), 
                                           digits=1),"%)",sep=""))))




# Using Alternative Thresholds of the AB60 exposure measure -----------
### Regressions
# a. accidents per capita x mean
m1 <- lm(acc.pc ~ treat * exposure.mean + as.factor(year2) + as.factor(month) + 
           as.factor(county), data=d1)
c1 <- coeftest(m1, vcov = vcovCluster(m1, factor(d1$county)))

# b. accidents per capita x median
m2 <- lm(acc.pc ~ treat * exposure.median + as.factor(year2) + as.factor(month) + as.factor(county), 
         data=d1)
c2 <- coeftest(m2, vcov = vcovCluster(m2, factor(d1$county)))

# c. % fatal.accidents x mean
m3 <- lm(fatal.share ~ treat * exposure.mean + as.factor(year2) + 
           as.factor(month) + as.factor(county), data=d2)
c3 <- coeftest(m3, vcov = vcovCluster(m3, factor(d2$county)))

# d. % fatal.accidents x median
m4 <- lm(fatal.share ~ treat * exposure.median + as.factor(year2) + as.factor(month) + as.factor(county), 
         data=d2)
c4 <- coeftest(m4, vcov = vcovCluster(m4, factor(d2$county)))

# e. % hit-and-run accidents x mean
m5 <- lm(hr.share ~ treat * exposure.mean + as.factor(year2) + as.factor(month) + 
           as.factor(county), data=d2)
c5 <- coeftest(m5, vcov = vcovCluster(m5, factor(d2$county)))

# f. % hit-and-run accidents x median
m6 <- lm(hr.share ~ treat * exposure.median + as.factor(year2) + as.factor(month) + as.factor(county), 
         data=d2)
c6 <- coeftest(m6, vcov = vcovCluster(m6, factor(d2$county)))


### output
stargazer(m1, m2, m3, m4, m5, m6,
          se=list(c1[,2], c2[,2], c3[,2], c4[,2], 
                  c5[,2], c6[,2]),
          omit.stat = c("rsq", "f", "ser"), df=F, notes.align="c", no.space=T,
          font.size="scriptsize",
          model.names = F,
          omit=c("month", "year", "county", "Constant"),
          dep.var.caption = "",
          dep.var.labels = c("Accidents per 1,000 capita", "Share Fatal Accidents", "Share Hit-and-Run Accidents"),
          star.cutoffs = c(0.05),
          add.lines=list(c("Implied percent change in outcome compared to 2014", 
                           paste(round(m1$coefficients[80] / mean(d$acc.pc[d$year==2014], na.rm=T) * 100, digits=1),"%",sep=""),
                           paste(round(m2$coefficients[80] / mean(d$acc.pc[d$year==2014], na.rm=T) * 100, digits=1),"%",sep=""),
                           paste(round(m3$coefficients[80] / mean(d$fatal.share[d$year==2014], na.rm=T) * 100, digits=1),"%",sep=""),
                           paste(round(m4$coefficients[80] / mean(d$fatal.share[d$year==2014], na.rm=T) * 100, digits=1),"%",sep=""),
                           paste(round(m5$coefficients[80] / mean(d$hr.share[d$year==2014], na.rm=T) * 100, digits=1), "%", sep=""),
                           paste(round(m6$coefficients[80] / mean(d$hr.share[d$year==2014], na.rm=T) * 100, digits=1), "%", sep="")),
                         c("",
                           paste("(",round(c1[79,2] / mean(d$acc.pc[d$year==2014], na.rm=T) * 100, digits=1),"%)",sep=""),
                           paste("(",round(c2[79,2] / mean(d$acc.pc[d$year==2014], na.rm=T) * 100, digits=1),"%)",sep=""),
                           paste("(",round(c3[79,2] / mean(d$fatal.share[d$year==2014], na.rm=T) * 100, digits=1),"%)",sep=""),
                           paste("(",round(c4[79,2] / mean(d$fatal.share[d$year==2014], na.rm=T) * 100, digits=1),"%)",sep=""),
                           paste("(",round(c5[79,2] / mean(d$hr.share[d$year==2014], na.rm=T) * 100, digits=1),"%)",sep=""),
                           paste("(",round(c6[79,2] / mean(d$hr.share[d$year==2014], na.rm=T) * 100, digits=1),"%)",sep=""))))




# Adding control variables --------------------------------------------

### Regressions
# a. accidents per capita x terciles
m1 <- lm(acc.pc ~ treat * as.factor(exposure.tercile) + unemploymentrate + log(income) + 
           as.factor(year2) + as.factor(month) + as.factor(county), data=d1)
c1 <- coeftest(m1, vcov = vcovCluster(m1, factor(d1$county)))

# b. accidents per capita x continuous
m2 <- lm(acc.pc ~ treat * exposure + unemploymentrate + log(income) + 
           as.factor(year2) + as.factor(month) + as.factor(county), data=d1)
c2 <- coeftest(m2, vcov = vcovCluster(m2, factor(d1$county)))

# c. % fatal accidents x terciles
m3 <- lm(fatal.share ~ treat * as.factor(exposure.tercile) + unemploymentrate + log(income) +
           as.factor(year2) + as.factor(month) + as.factor(county), data=d2)
c3 <- coeftest(m3, vcov = vcovCluster(m3, factor(d2$county)))

# d. % fatal accidents x continuous
m4 <- lm(fatal.share ~ treat * exposure + unemploymentrate + log(income) + 
           as.factor(year2) + as.factor(month) + as.factor(county), data=d2)
c4 <- coeftest(m4, vcov = vcovCluster(m4, factor(d2$county)))

# e. % hit-and-run accidents x terciles
m5 <- lm(hr.share ~ treat * as.factor(exposure.tercile) + unemploymentrate + log(income) +
           as.factor(year2) + as.factor(month) + as.factor(county), data=d2)
c5 <- coeftest(m5, vcov = vcovCluster(m5, factor(d2$county)))

# f. % hit-and-run accidents x continuous
m6 <- lm(hr.share ~ treat * exposure + unemploymentrate + log(income) + 
           as.factor(year2) + as.factor(month) + as.factor(county), data=d2)
c6 <- coeftest(m6, vcov = vcovCluster(m6, factor(d2$county)))


### output
stargazer(m1, m2, m3, m4, m5, m6,
          se=list(c1[,2], c2[,2], c3[,2], c4[,2], 
                  c5[,2], c6[,2]),
          omit.stat = c("rsq", "f", "ser"), df=F, notes.align="c", no.space=T,
          font.size="scriptsize",
          model.names = F,
          omit=c("month", "year", "county", "Constant"),
          dep.var.caption = "",
          dep.var.labels = c("Accidents per 1,000 capita", "Share Fatal Accidents", "Share Hit-and-Run Accidents"),
          star.cutoffs = c(0.05),
          add.lines=list(c("Implied percent change in outcome compared to 2014", 
                           paste(round(m1$coefficients[84] / mean(d$acc.pc[d$year==2014], na.rm=T) * 100, digits=1),"%",sep=""),
                           paste(round(m2$coefficients[82] / mean(d$acc.pc[d$year==2014], na.rm=T) * 100 * mean(d$exposure, na.rm=T), 
                                       digits=1),"%",sep=""),
                           paste(round(m3$coefficients[84] / mean(d$fatal.share[d$year==2014], na.rm=T) * 100, digits=1),"%",sep=""),
                           paste(round(m4$coefficients[82] / mean(d$fatal.share[d$year==2014], na.rm=T) * 100 * mean(d$exposure, na.rm=T), 
                                       digits=1),"%",sep=""),
                           paste(round(m5$coefficients[84] / mean(d$hr.share[d$year==2014], na.rm=T) * 100, digits=1), "%", sep=""),
                           paste(round(m6$coefficients[82] / mean(d$hr.share[d$year==2014], na.rm=T) * 100 * mean(d$exposure, na.rm=T),
                                       digits=1), "%", sep="")),
                         c("",
                           paste("(",round(c1[82,2] / mean(d$acc.pc[d$year==2014], na.rm=T) * 100, digits=1),"%)",sep=""),
                           paste("(",round(c2[81,2] / mean(d$acc.pc[d$year==2014], na.rm=T) * 100 * mean(d$exposure, na.rm=T), 
                                           digits=1),"%)",sep=""),
                           paste("(",round(c3[82,2] / mean(d$fatal.share[d$year==2014], na.rm=T) * 100, digits=1),"%)",sep=""),
                           paste("(",round(c4[81,2] / mean(d$fatal.share[d$year==2014], na.rm=T) * 100 * mean(d$exposure, na.rm=T), 
                                           digits=1),"%)",sep=""),
                           paste("(",round(c5[82,2] / mean(d$hr.share[d$year==2014], na.rm=T) * 100, digits=1),"%)",sep=""),
                           paste("(",round(c6[81,2] / mean(d$hr.share[d$year==2014], na.rm=T) * 100 * mean(d$exposure, na.rm=T), 
                                           digits=1),"%)",sep=""))))



# Restricting the sample period to 14/15 ----------------------------------

### Regressions
# a. accidents per capita x terciles
m1 <- lm(acc.pc ~ treat * as.factor(exposure.tercile) + as.factor(year2) + as.factor(month) + 
           as.factor(county), data=d1, d1$year>=2014 & d1$year <= 2015)
c1 <- coeftest(m1, vcov = vcovCluster(m1, factor(d1$county[d1$year>=2014 & d1$year <= 2015])))

# b. accidents per capita x continuous
m2 <- lm(acc.pc ~ treat * exposure + as.factor(year2) + as.factor(month) + as.factor(county), 
         data=d1, d1$year>=2014 & d1$year <= 2015)
c2 <- coeftest(m2, vcov = vcovCluster(m2, factor(d1$county[d1$year>=2014 & d1$year <= 2015])))

# c. % fatal accidents x terciles
m3 <- lm(fatal.share ~ treat * as.factor(exposure.tercile) + as.factor(year2) + 
           as.factor(month) + as.factor(county), data=d2, d2$year>=2014 & d2$year <= 2015)
c3 <- coeftest(m3, vcov = vcovCluster(m3, factor(d2$county[d2$year>=2014 & d2$year <= 2015])))

# d. % fatal accidents x continuous
m4 <- lm(fatal.share ~ treat * exposure + as.factor(year2) + as.factor(month) + as.factor(county), 
         data=d2, d2$year>=2014 & d2$year <= 2015)
c4 <- coeftest(m4, vcov = vcovCluster(m4, factor(d2$county[d2$year>=2014 & d2$year <= 2015])))

# e. % hit-and-run accidents x terciles
m5 <- lm(hr.share ~ treat * as.factor(exposure.tercile) + as.factor(year2) + as.factor(month) + 
           as.factor(county), data=d2, d2$year>=2014 & d2$year <= 2015)
c5 <- coeftest(m5, vcov = vcovCluster(m5, factor(d2$county[d2$year>=2014 & d2$year <= 2015])))

# f. % hit-and-run accidents x continuous
m6 <- lm(hr.share ~ treat * exposure + as.factor(year2) + as.factor(month) + as.factor(county), 
         data=d2, d2$year>=2014 & d2$year <= 2015)
c6 <- coeftest(m6, vcov = vcovCluster(m6, factor(d2$county[d2$year>=2014 & d2$year <= 2015])))


### output
stargazer(m1, m2, m3, m4, m5, m6,
          se=list(c1[,2], c2[,2], c3[,2], c4[,2], 
                  c5[,2], c6[,2]),
          omit.stat = c("rsq", "f", "ser"), df=F, notes.align="c", no.space=T,
          font.size="scriptsize",
          model.names = F,
          omit=c("month", "year", "county", "Constant"),
          dep.var.caption = "",
          dep.var.labels = c("Accidents per 1,000 capita", "Share Fatal Accidents", "Share Hit-and-Run Accidents"),
          star.cutoffs = c(0.05),
          add.lines=list(c("Implied percent change in outcome compared to 2014", 
                           paste(round(m1$coefficients[75] / mean(d$acc.pc[d$year==2014], na.rm=T) * 100, digits=1),"%",sep=""),
                           paste(round(m2$coefficients[73] / mean(d$acc.pc[d$year==2014], na.rm=T) * 100 * mean(d$exposure, na.rm=T), 
                                       digits=1),"%",sep=""),
                           paste(round(m3$coefficients[75] / mean(d$fatal.share[d$year==2014], na.rm=T) * 100, digits=1),"%",sep=""),
                           paste(round(m4$coefficients[73] / mean(d$fatal.share[d$year==2014], na.rm=T) * 100 * mean(d$exposure, na.rm=T), 
                                       digits=1),"%",sep=""),
                           paste(round(m5$coefficients[75] / mean(d$hr.share[d$year==2014], na.rm=T) * 100, digits=1), "%", sep=""),
                           paste(round(m6$coefficients[73] / mean(d$hr.share[d$year==2014], na.rm=T) * 100 * mean(d$exposure, na.rm=T),
                                       digits=1), "%", sep="")),
                         c("",
                           paste("(",round(c1[72,2] / mean(d$acc.pc[d$year==2014], na.rm=T) * 100, digits=1),"%)",sep=""),
                           paste("(",round(c2[71,2] / mean(d$acc.pc[d$year==2014], na.rm=T) * 100 * mean(d$exposure, na.rm=T), 
                                           digits=1),"%)",sep=""),
                           paste("(",round(c3[72,2] / mean(d$fatal.share[d$year==2014], na.rm=T) * 100, digits=1),"%)",sep=""),
                           paste("(",round(c4[71,2] / mean(d$fatal.share[d$year==2014], na.rm=T) * 100 * mean(d$exposure, na.rm=T), 
                                           digits=1),"%)",sep=""),
                           paste("(",round(c5[72,2] / mean(d$hr.share[d$year==2014], na.rm=T) * 100, digits=1),"%)",sep=""),
                           paste("(",round(c6[71,2] / mean(d$hr.share[d$year==2014], na.rm=T) * 100 * mean(d$exposure, na.rm=T), 
                                           digits=1),"%)",sep=""))))



### mean and SDs for DVs
# accidents per capita
paste("Mean = ", round(mean(d1$acc.pc[d1$year>=2014 & d1$year<=2015], na.rm=T), digits=3), 
      ", SD = ",round(sd(d1$acc.pc[d1$year>=2014 & d1$year<=2015], na.rm=T), digits=3), sep="")

# % fatal accidents
paste("Mean = ", round(mean(d2$fatal.share[d2$year>=2014 & d2$year<=2015], na.rm=T), digits=3), 
      ", SD = ",round(sd(d2$fatal.share[d2$year>=2014 & d2$year<=2015], na.rm=T), digits=3), sep="")

# % hit-and-run accidents
paste("Mean = ", round(mean(d2$hr.share[d2$year>=2014 & d2$year<=2015], na.rm=T), digits=3), 
      ", SD = ",round(sd(d2$hr.share[d2$year>=2014 & d2$year<=2015], na.rm=T), digits=3), sep="")


# Other Traffic Safety Outcomes -------------------------------------------

### Regressions
# a. fatal accidents per capita x terciles
m1 <- lm(acc.fatal.pc ~ treat * as.factor(exposure.tercile) + as.factor(year2) + as.factor(month) + 
           as.factor(county), data=d2)
c1 <- coeftest(m1, vcov = vcovCluster(m1, factor(d2$county)))

# b. fatal accidents per capita x continuous
m2 <- lm(acc.fatal.pc ~ treat * exposure + as.factor(year2) + as.factor(month) + as.factor(county), 
         data=d2)
c2 <- coeftest(m2, vcov = vcovCluster(m2, factor(d2$county)))

# c. people killed in accidents per capita x terciles
m3 <- lm(acc.killed.pc ~ treat * as.factor(exposure.tercile) + as.factor(year2) + 
           as.factor(month) + as.factor(county), data=d2)
c3 <- coeftest(m3, vcov = vcovCluster(m3, factor(d2$county)))

# d. people killed in accidents per capita x continuous
m4 <- lm(acc.killed.pc ~ treat * exposure + as.factor(year2) + as.factor(month) + as.factor(county), 
         data=d2)
c4 <- coeftest(m4, vcov = vcovCluster(m4, factor(d2$county)))

# e. people killed in accidents per capita x terciles
m5 <- lm(hr.pc ~ treat * as.factor(exposure.tercile) + as.factor(year2) + 
           as.factor(month) + as.factor(county), data=d2)
c5 <- coeftest(m5, vcov = vcovCluster(m5, factor(d2$county)))

# f. people killed in accidents per capita x continuous
m6 <- lm(hr.pc ~ treat * exposure + as.factor(year2) + as.factor(month) + as.factor(county), 
         data=d2)
c6 <- coeftest(m6, vcov = vcovCluster(m6, factor(d2$county)))

# g. people killed in accidents per capita x terciles
m7 <- lm(acc.nonhr.pc ~ treat * as.factor(exposure.tercile) + as.factor(year2) + 
           as.factor(month) + as.factor(county), data=d1)
c7 <- coeftest(m7, vcov = vcovCluster(m7, factor(d1$county)))

# h. people killed in accidents per capita x continuous
m8 <- lm(acc.nonhr.pc ~ treat * exposure + as.factor(year2) + as.factor(month) + as.factor(county), 
         data=d1)
c8 <- coeftest(m8, vcov = vcovCluster(m8, factor(d1$county)))



### output
stargazer(m1, m2, m3, m4,m5,m6,m7,m8,
          se=list(c1[,2], c2[,2], c3[,2], c4[,2], c5[,2], c6[,2], c7[,2], c8[,2]),
          omit.stat = c("rsq", "f", "ser"), df=F, notes.align="c", no.space=T,
          font.size="scriptsize",
          model.names = F,
          omit=c("month", "year", "county", "Constant"),
          dep.var.caption = "",
          dep.var.labels = c("Fatal accidents per 1,000 capita", "People killed in accidents per 1,000 capita",
                             "Hit-and-Run accidents per 1,000 capita", "Non-Hit-and-Run accidents per 1,000 capita"),
          star.cutoffs = c(0.05),
          add.lines=list(c("Implied percent change in outcome compared to 2014", 
                           paste(round(m1$coefficients[82] / mean(d$acc.fatal.pc[d$year==2014], na.rm=T) * 100, digits=1),"%",sep=""),
                           paste(round(m2$coefficients[80] / mean(d$acc.fatal.pc[d$year==2014], na.rm=T) * 100 * mean(d$exposure, na.rm=T), 
                                       digits=1),"%",sep=""),
                           paste(round(m3$coefficients[82] / mean(d$acc.killed.pc[d$year==2014], na.rm=T) * 100, digits=1),"%",sep=""),
                           paste(round(m4$coefficients[80] / mean(d$acc.killed.pc[d$year==2014], na.rm=T) * 100 * mean(d$exposure, na.rm=T), 
                                       digits=1),"%",sep=""),
                           paste(round(m5$coefficients[82] / mean(d$hr.pc[d$year==2014], na.rm=T) * 100, digits=1),"%",sep=""),
                           paste(round(m6$coefficients[80] / mean(d$hr.pc[d$year==2014], na.rm=T) * 100 * mean(d$exposure, na.rm=T), 
                                       digits=1),"%",sep=""),
                           paste(round(m7$coefficients[82] / mean(d$acc.nonhr.pc[d$year==2014], na.rm=T) * 100, digits=1),"%",sep=""),
                           paste(round(m8$coefficients[80] / mean(d$acc.nonhr.pc[d$year==2014], na.rm=T) * 100 * mean(d$exposure, na.rm=T), 
                                       digits=1),"%",sep="")),
                         c("",
                           paste("(",round(c1[80,2] / mean(d$acc.fatal.pc[d$year==2014], na.rm=T) * 100, digits=1),"%)",sep=""),
                           paste("(",round(c2[79,2] / mean(d$acc.fatal.pc[d$year==2014], na.rm=T) * 100 * mean(d$exposure, na.rm=T), 
                                           digits=1),"%)",sep=""),
                           paste("(",round(c3[80,2] / mean(d$acc.killed.pc[d$year==2014], na.rm=T) * 100, digits=1),"%)",sep=""),
                           paste("(",round(c4[79,2] / mean(d$acc.killed.pc[d$year==2014], na.rm=T) * 100 * mean(d$exposure, na.rm=T), 
                                           digits=1),"%)",sep=""),
                           paste("(",round(c5[80,2] / mean(d$hr.pc[d$year==2014], na.rm=T) * 100, digits=1),"%)",sep=""),
                           paste("(",round(c6[79,2] / mean(d$hr.pc[d$year==2014], na.rm=T) * 100 * mean(d$exposure, na.rm=T), 
                                           digits=1),"%)",sep=""),
                           paste("(",round(c7[80,2] / mean(d$acc.nonhr.pc[d$year==2014], na.rm=T) * 100, digits=1),"%)",sep=""),
                           paste("(",round(c8[79,2] / mean(d$acc.nonhr.pc[d$year==2014], na.rm=T) * 100 * mean(d$exposure, na.rm=T), 
                                           digits=1),"%)",sep=""))))


### mean and SDs for DVs
# fatal accidents per capita
paste("Mean = ", round(mean(d1$acc.fatal.pc, na.rm=T), digits=3), 
      ", SD = ",round(sd(d1$acc.fatal.pc, na.rm=T), digits=3), sep="")

# killed in accidents per capita
paste("Mean = ", round(mean(d1$acc.killed.pc, na.rm=T), digits=3), 
      ", SD = ",round(sd(d1$acc.killed.pc, na.rm=T), digits=3), sep="")

# hit-and-run accidents per capita
paste("Mean = ", round(mean(d1$hr.pc, na.rm=T), digits=3), 
      ", SD = ",round(sd(d1$hr.pc, na.rm=T), digits=3), sep="")

# non-hit-and-run accidents per capita
paste("Mean = ", round(mean(d1$acc.nonhr.pc, na.rm=T), digits=3), 
      ", SD = ",round(sd(d1$acc.nonhr.pc, na.rm=T), digits=3), sep="")



# Placebo tests: AB353 and AB4 --------------------------------------------

# code new treatment periods
d$impoundlaw <- ifelse(d$year<2012, 0, 1)
d$trustact <- ifelse(d$year<2014, 0, 1)

d3 <- subset(d, d$year>=2011 & d$year<=2012)
d4 <- subset(d, d$year>=2013 & d$year<=2014)
d5 <- d4[!is.na(d4$hr.share),]


### 1: Anti-Impound Law, 2011/12
## Regressions
# a. accidents per capita x terciles
m1 <- lm(acc.pc ~ impoundlaw * as.factor(exposure.tercile) + as.factor(year2) + as.factor(month) + 
           as.factor(county), data=d3)
c1 <- coeftest(m1, vcov = vcovCluster(m1, factor(d3$county)))

# b. accidents per capita x continuous
m2 <- lm(acc.pc ~ impoundlaw * exposure + as.factor(year2) + as.factor(month) + as.factor(county), 
         data=d3)
c2 <- coeftest(m2, vcov = vcovCluster(m2, factor(d3$county)))

# c. % fatal accidents x terciles
m3 <- lm(fatal.share ~ impoundlaw * as.factor(exposure.tercile) + as.factor(year2) + 
           as.factor(month) + as.factor(county), data=d3)
c3 <- coeftest(m3, vcov = vcovCluster(m3, factor(d3$county)))

# d. % fatal accidents x continuous
m4 <- lm(fatal.share ~ impoundlaw * exposure + as.factor(year2) + as.factor(month) + as.factor(county), 
         data=d3)
c4 <- coeftest(m4, vcov = vcovCluster(m4, factor(d3$county)))

# e. % hit-and-run accidents x terciles
m5 <- lm(hr.share ~ impoundlaw * as.factor(exposure.tercile) + as.factor(year2) + as.factor(month) + 
           as.factor(county), data=d3)
c5 <- coeftest(m5, vcov = vcovCluster(m5, factor(d3$county)))

# f. % hit-and-run accidents x continuous
m6 <- lm(hr.share ~ impoundlaw * exposure + as.factor(year2) + as.factor(month) + as.factor(county), 
         data=d3)
c6 <- coeftest(m6, vcov = vcovCluster(m6, factor(d3$county)))


## output
stargazer(m1, m2, m3, m4, m5, m6,
          se=list(c1[,2], c2[,2], c3[,2], c4[,2], 
                  c5[,2], c6[,2]),
          omit.stat = c("rsq", "f", "ser"), df=F, notes.align="c", no.space=T,
          font.size="scriptsize",
          model.names = F,
          omit=c("month", "year", "county", "Constant"),
          dep.var.caption = "",
          dep.var.labels = c("Accidents per 1,000 capita", "Share Fatal Accidents", "Share Hit-and-Run Accidents"),
          star.cutoffs = c(0.05),
          add.lines=list(c("Implied percent change in outcome compared to 2011", 
                           paste(round(m1$coefficients[75] / mean(d$acc.pc[d$year==2011], na.rm=T) * 100, digits=1),"%",sep=""),
                           paste(round(m2$coefficients[73] / mean(d$acc.pc[d$year==2011], na.rm=T) * 100 * mean(d$exposure, na.rm=T), 
                                       digits=1),"%",sep=""),
                           paste(round(m3$coefficients[75] / mean(d$fatal.share[d$year==2011], na.rm=T) * 100, digits=1),"%",sep=""),
                           paste(round(m4$coefficients[73] / mean(d$fatal.share[d$year==2011], na.rm=T) * 100 * mean(d$exposure, na.rm=T), 
                                       digits=1),"%",sep=""),
                           paste(round(m5$coefficients[75] / mean(d$hr.share[d$year==2011], na.rm=T) * 100, digits=1), "%", sep=""),
                           paste(round(m6$coefficients[73] / mean(d$hr.share[d$year==2011], na.rm=T) * 100 * mean(d$exposure, na.rm=T),
                                       digits=1), "%", sep="")),
                         c("",
                           paste("(",round(c1[72,2] / mean(d$acc.pc[d$year==2011], na.rm=T) * 100, digits=1),"%)",sep=""),
                           paste("(",round(c2[71,2] / mean(d$acc.pc[d$year==2011], na.rm=T) * 100 * mean(d$exposure, na.rm=T), 
                                           digits=1),"%)",sep=""),
                           paste("(",round(c3[72,2] / mean(d$fatal.share[d$year==2011], na.rm=T) * 100, digits=1),"%)",sep=""),
                           paste("(",round(c4[71,2] / mean(d$fatal.share[d$year==2011], na.rm=T) * 100 * mean(d$exposure, na.rm=T), 
                                           digits=1),"%)",sep=""),
                           paste("(",round(c5[72,2] / mean(d$hr.share[d$year==2011], na.rm=T) * 100, digits=1),"%)",sep=""),
                           paste("(",round(c6[71,2] / mean(d$hr.share[d$year==2011], na.rm=T) * 100 * mean(d$exposure, na.rm=T), 
                                           digits=1),"%)",sep=""))))


### mean and SDs for DVs
# accidents per capita
paste("Mean = ", round(mean(d3$acc.pc, na.rm=T), digits=3), 
      ", SD = ",round(sd(d3$acc.pc, na.rm=T), digits=3), sep="")

# % fatal accidents
paste("Mean = ", round(mean(d3$fatal.share, na.rm=T), digits=3), 
      ", SD = ",round(sd(d3$fatal.share, na.rm=T), digits=3), sep="")

# % hit-and-run accidents
paste("Mean = ", round(mean(d3$hr.share, na.rm=T), digits=3), 
      ", SD = ",round(sd(d3$hr.share, na.rm=T), digits=3), sep="")




### 2: TRUST Act, 2013/14
## Regressions
# a. accidents per capita x terciles
m1 <- lm(acc.pc ~ trustact * as.factor(exposure.tercile) + as.factor(year2) + as.factor(month) + 
           as.factor(county), data=d4)
c1 <- coeftest(m1, vcov = vcovCluster(m1, factor(d4$county)))

# b. accidents per capita x continuous
m2 <- lm(acc.pc ~ trustact * exposure + as.factor(year2) + as.factor(month) + as.factor(county), 
         data=d4)
c2 <- coeftest(m2, vcov = vcovCluster(m2, factor(d4$county)))

# c. % fatal accidents x terciles
m3 <- lm(fatal.share ~ trustact * as.factor(exposure.tercile) + as.factor(year2) + 
           as.factor(month) + as.factor(county), data=d5)
c3 <- coeftest(m3, vcov = vcovCluster(m3, factor(d5$county)))

# d. % fatal accidents x continuous
m4 <- lm(fatal.share ~ trustact * exposure + as.factor(year2) + as.factor(month) + as.factor(county), 
         data=d5)
c4 <- coeftest(m4, vcov = vcovCluster(m4, factor(d5$county)))

# e. % hit-and-run accidents x terciles
m5 <- lm(hr.share ~ trustact * as.factor(exposure.tercile) + as.factor(year2) + as.factor(month) + 
           as.factor(county), data=d5)
c5 <- coeftest(m5, vcov = vcovCluster(m5, factor(d5$county)))

# f. % hit-and-run accidents x continuous
m6 <- lm(hr.share ~ trustact * exposure + as.factor(year2) + as.factor(month) + as.factor(county), 
         data=d5)
c6 <- coeftest(m6, vcov = vcovCluster(m6, factor(d5$county)))


## output
stargazer(m1, m2, m3, m4, m5, m6,
          se=list(c1[,2], c2[,2], c3[,2], c4[,2], 
                  c5[,2], c6[,2]),
          omit.stat = c("rsq", "f", "ser"), df=F, notes.align="c", no.space=T,
          font.size="scriptsize",
          model.names = F,
          omit=c("month", "year", "county", "Constant"),
          dep.var.caption = "",
          dep.var.labels = c("Accidents per 1,000 capita", "Share Fatal Accidents", "Share Hit-and-Run Accidents"),
          star.cutoffs = c(0.05),
          add.lines=list(c("Implied percent change in outcome compared to 2013", 
                           paste(round(m1$coefficients[75] / mean(d$acc.pc[d$year==2013], na.rm=T) * 100, digits=1),"%",sep=""),
                           paste(round(m2$coefficients[73] / mean(d$acc.pc[d$year==2013], na.rm=T) * 100 * mean(d$exposure, na.rm=T), 
                                       digits=1),"%",sep=""),
                           paste(round(m3$coefficients[75] / mean(d$fatal.share[d$year==2013], na.rm=T) * 100, digits=1),"%",sep=""),
                           paste(round(m4$coefficients[73] / mean(d$fatal.share[d$year==2013], na.rm=T) * 100 * mean(d$exposure, na.rm=T), 
                                       digits=1),"%",sep=""),
                           paste(round(m5$coefficients[75] / mean(d$hr.share[d$year==2013], na.rm=T) * 100, digits=1), "%", sep=""),
                           paste(round(m6$coefficients[73] / mean(d$hr.share[d$year==2013], na.rm=T) * 100 * mean(d$exposure, na.rm=T),
                                       digits=1), "%", sep="")),
                         c("",
                           paste("(",round(c1[72,2] / mean(d$acc.pc[d$year==2013], na.rm=T) * 100, digits=1),"%)",sep=""),
                           paste("(",round(c2[71,2] / mean(d$acc.pc[d$year==2013], na.rm=T) * 100 * mean(d$exposure, na.rm=T), 
                                           digits=1),"%)",sep=""),
                           paste("(",round(c3[72,2] / mean(d$fatal.share[d$year==2013], na.rm=T) * 100, digits=1),"%)",sep=""),
                           paste("(",round(c4[71,2] / mean(d$fatal.share[d$year==2013], na.rm=T) * 100 * mean(d$exposure, na.rm=T), 
                                           digits=1),"%)",sep=""),
                           paste("(",round(c5[72,2] / mean(d$hr.share[d$year==2013], na.rm=T) * 100, digits=1),"%)",sep=""),
                           paste("(",round(c6[71,2] / mean(d$hr.share[d$year==2013], na.rm=T) * 100 * mean(d$exposure, na.rm=T), 
                                           digits=1),"%)",sep=""))))


### mean and SDs for DVs
# accidents per capita
paste("Mean = ", round(mean(d4$acc.pc, na.rm=T), digits=3), 
      ", SD = ",round(sd(d4$acc.pc, na.rm=T), digits=3), sep="")

# % fatal accidents
paste("Mean = ", round(mean(d5$fatal.share, na.rm=T), digits=3), 
      ", SD = ",round(sd(d5$fatal.share, na.rm=T), digits=3), sep="")

# % hit-and-run accidents
paste("Mean = ", round(mean(d5$hr.share, na.rm=T), digits=3), 
      ", SD = ",round(sd(d5$hr.share, na.rm=T), digits=3), sep="")



