## Performing synthetic control on CALIFORNIA 
ACC_usa <- read.csv("ACC_usa_data_compiled.csv", header=TRUE)

#### The data set contains following Variables
# State.Code
# State.Name
# YEAR
# acc_total: total number of accident by state and year (From FARS)
# hr_total: total number of hit and run accidents by state and year (From FARS)
# pop: Total Population by year (from US Census Bureau)
# income: Total Median Income by year (from US Census Bureau)
# unemp: Rate of unemployment (From Labor Market Information)


# Removing data for 12 states other than California that gives driving license to unauthorized immigrants
# Colorado, Connecticut, Delaware, DC, Hawaii, Illinois, Maryland, Nevada, New Mexico, Utah, Vermont, Washington
ls.skip_state_code = c(8,9,10,11,15, 17,24,32,35,49,50,53)
ls.skip_state_name = c("COLORADO",
                       "CONNECTICUT",
                       "DELAWARE",
                       "DISTRICT OF COLUMBIA",
                       "HAWAII",
                       "ILLINOIS",
                       "MARYLAND",
                       "NEVADA",
                       "NEW MEXICO",
                       "UTAH",
                       "VERMONT",
                       "WASSHINGTON")
for (i in ls.skip_state_code) {
  ACC_usa = subset(ACC_usa, State.Code != i)
}

# Calculating Number of Accidents per 1000 population
ACC_usa$hr_pc <- (ACC_usa$hr_total / ACC_usa$pop)*1000
ACC_usa$hr_share <- (ACC_usa$hr_total / ACC_usa$acc_total)*100
ACC_usa$acc_pc <- (ACC_usa$acc_total / ACC_usa$pop)*1000

# Changing State.Name variable to character
ACC_usa$State.Name <- as.character(ACC_usa$State.Name)

## Performing Synthetic Control Number of accidents per 1000 population
dataprep.out <- dataprep(foo = ACC_usa,
                         predictors = c("pop", "income", "unemp"),
                         predictors.op = c("median"),
                         dependent = c("acc_pc"),
                         unit.variable = c("State.Code"),
                         time.variable = "YEAR",
                         treatment.identifier = 6,
                         controls.identifier = unique(ACC_usa$State.Code)[-5],
                         time.predictors.prior = c(2006:2014),
                         time.optimize.ssr = c(2010:2014),
                         time.plot = c(2006:2016),
                         unit.names.variable = c("State.Name")
)

synth.out <- synth(dataprep.out)
round(synth.out$solution.w,2)
synth.out$solution.v

gaps<- dataprep.out$Y1plot-(
  dataprep.out$Y0plot%*%synth.out$solution.w
) ; gaps

synth.tables <- synth.tab(
  dataprep.res = dataprep.out,
  synth.res = synth.out)
print(synth.tables)

path.plot(dataprep.res = dataprep.out,
          synth.res = synth.out,
          Ylab = c("Number of Accidents per 1000 population"),
          Xlab = c("Year"),
          Main = "Synthetic Control for California to find the impact of AB60")

# Line for Year of Treatment
abline(v=2015, lwd=1, lty=2)

## plot the gaps (treated - synthetic)
gaps.plot(dataprep.res = dataprep.out,synth.res = synth.out)


## Predicting Sythetic Control for share of hit and run accidents
dataprep2.out <- dataprep(foo = ACC_usa,
                         predictors = c("pop", "income", "unemp"),
                         predictors.op = c("mean"),
                         dependent = c("hr_share"),
                         unit.variable = c("State.Code"),
                         time.variable = "YEAR",
                         treatment.identifier = 6,
                         controls.identifier = unique(ACC_usa$State.Code)[-5],
                         time.predictors.prior = c(2006:2014),
                         time.optimize.ssr = c(2006:2014),
                         time.plot = c(2006:2016),
                         unit.names.variable = c("State.Name")
)

synth2.out <- synth(dataprep2.out)
round(synth2.out$solution.w,2)
synth2.out$solution.v


gaps2<- dataprep2.out$Y1plot-(
  dataprep2.out$Y0plot%*%synth2.out$solution.w
) ; gaps2

synth2.tables <- synth.tab(
  dataprep.res = dataprep2.out,
  synth.res = synth2.out)
print(synth2.tables)

path.plot(dataprep.res = dataprep2.out,
          synth.res = synth2.out,
          Ylab = c("Percentage of Hit and run cases"),
          Xlab = c("Year"),
          Main = "Synthetic Control for California to find the impact of AB60")

# Line for Year of Treatment
abline(v=2015, lwd=1, lty=2)

## plot the gaps (treated - synthetic)
gaps.plot(dataprep.res = dataprep2.out,synth.res = synth2.out)




### PLacebo In time
## Performing Synthetic Control Number of accidents per 1000 population
dataprep_pt.out <- dataprep(foo = ACC_usa,
                         predictors = c("pop", "income", "unemp"),
                         predictors.op = c("median"),
                         dependent = c("acc_pc"),
                         unit.variable = c("State.Code"),
                         time.variable = "YEAR",
                         treatment.identifier = 6,
                         controls.identifier = unique(ACC_usa$State.Code)[-5],
                         time.predictors.prior = c(2006:2013),
                         time.optimize.ssr = c(2006:2013),
                         time.plot = c(2006:2016),
                         unit.names.variable = c("State.Name")
)

synth_pt.out <- synth(dataprep_pt.out)
round(synth.out$solution.w,2)
synth_pt.out$solution.v


gaps_pt<- dataprep_pt.out$Y1plot-(
  dataprep_pt.out$Y0plot%*%synth_pt.out$solution.w
) ; gaps_pt

synth_pt.tables <- synth.tab(
  dataprep.res = dataprep_pt.out,
  synth.res = synth_pt.out)
print(synth.tables)

path.plot(dataprep.res = dataprep_pt.out,
          synth.res = synth_pt.out,
          Ylab = c("Number of Accidents per 1000 population"),
          Xlab = c("Year"),
          Main = "Sythetic Control with Placebo in Time")

# Line for Year of Placebo in Time
abline(v=2013)
# Line for Year of Treatment
abline(v=2015, lwd=1, lty=2)

## plot the gaps (treated - synthetic)
gaps.plot(dataprep.res = dataprep_pt.out,synth.res = synth_pt.out)


### PLacebo In Space ---------
## Performing Synthetic Control Number of accidents per 1000 population
controls = unique(ACC_usa$State.Code)[-24]
controls = controls[-5]
# 24 and 5 are indexes of NEW YORK and CALIFORNIA
dataprep_ps.out <- dataprep(foo = ACC_usa,
                         predictors = c("pop", "income", "unemp"),
                         predictors.op = c("median"),
                         dependent = c("acc_pc"),
                         unit.variable = c("State.Code"),
                         time.variable = "YEAR",
                         treatment.identifier = 36,
                         controls.identifier = unique(ACC_usa$State.Code)[-24][-5],
                         time.predictors.prior = c(2006:2014),
                         time.optimize.ssr = c(2006:2014),
                         time.plot = c(2006:2016),
                         unit.names.variable = c("State.Name")
)

synth_ps.out <- synth(dataprep_ps.out)
round(synth_ps.out$solution.w,2)
synth_ps.out$solution.v

gaps_ps<- dataprep_ps.out$Y1plot-(
  dataprep_ps.out$Y0plot%*%synth_ps.out$solution.w
) ; gaps_ps

synth_ps.tables <- synth.tab(
  dataprep.res = dataprep_ps.out,
  synth.res = synth_ps.out)
print(synth_ps.tables)

path.plot(dataprep.res = dataprep_ps.out,
          synth.res = synth_ps.out,
          Ylab = c("Number of Accidents per 1000 population"),
          Xlab = c("Year"),
          Main = "Synthetic Control: Placebo in Place on New York")

# Line for Year of Treatment
abline(v=2015, lwd=1, lty=2)

## plot the gaps (treated - synthetic)
gaps.plot(dataprep.res = dataprep_ps.out,synth.res = synth_ps.out)

