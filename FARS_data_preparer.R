setwd("~/Desktop/Minerva/CS112/Final_Project/")

#### Preparing Data for Replication ------------------------
#ACC_AUX.csv is provided by FARS and obtained from the following link
ACC_AUX <- read.csv("https://query.data.world/s/gsycqrxu62oelfgy744o3ernkm3hei", header=TRUE);
View(ACC_AUX)
#FRPP_GLC_US.csv give the code for states and county that are used in FARS and other governmental dataset
frpp_glc <- read.csv("FRPP_GLC_US.csv", header = TRUE, stringsAsFactors = FALSE, skip = 1)
View(frpp_glc)

dim(frpp_glc)
dim(ACC_AUX)

# Creating data set of State Name and respective Code
df.state = data.frame(unique(frpp_glc$State.Code),unique(frpp_glc$State.Name))

# The dataset with demographics provided with the main paper starts from 2006
ACC_subset <- subset(ACC_AUX, YEAR > 2005)

# The state code of CALIFORNIA is 6
ACC_cali <- subset(ACC_AUX, STATE == 6)
# Creating a data frame for California with only required information like type of crash, hit and run
ACC_cali_cleared <- data.frame(ACC_cali$A_CRAINJ, ACC_cali$A_HR, ACC_cali$YEAR, ACC_cali$STATE, ACC_cali$COUNTY)

# Getting and creating datagrame with all the Counties in CALIFORNIA
df.state.cali = subset(frpp_glc, frpp_glc$State.Name == "CALIFORNIA")
df <- data.frame(code = df.state.cali$County.Code, name = df.state.cali$County.Name)
df.state.cali.county_data = unique(df[c("code", "name")])
df.state.cali.county_data <- df.state.cali.county_data[order(df.state.cali.county_data$code),]

# Picking only data above 2005 because the initial dataset only has data from 2005
ACC_cali_cleared <- subset(ACC_cali_cleared, ACC_cali.YEAR > 2005)
num_years = length(unique(ACC_cali_cleared$ACC_cali.YEAR)) #number of years we are replicating

# Repeating the years and counties to prepare data in same format as the initial dataset
# Repeating years based on number of counties
year_repeat <- data.frame(
  rep(unique(ACC_cali_cleared$ACC_cali.YEAR), length(df.state.cali.county_data$code)))
df = df.state.cali.county_data
# Repeating Counties based on number of years
df2 = df[rep(seq_len(nrow(df)), each=num_years),]
row.names(df2) <- c(1:nrow(df2))
df2 = cbind(df2, year_repeat)

# Creating empty columns to store data about crash and hit and runs
ACC_cali <- data.frame(df2,
                rep(0, nrow(df2)),
                rep(0, nrow(df2))
                )
names(ACC_cali) <- c("County.Code", "County.Name","YEAR", "CRASH","HR")
for (row_index in c(1:nrow(ACC_cali))) {
  row = ACC_cali[row_index, ]
  ACC_cali[row_index,4] <- nrow(subset(ACC_cali_cleared, ACC_cali.YEAR == row$YEAR & ACC_cali.COUNTY == row$County.Code))
  ACC_cali[row_index,5] <- nrow(subset(ACC_cali_cleared, ACC_cali.YEAR == row$YEAR & ACC_cali.COUNTY == row$County.Code & ACC_cali.A_HR == 1))
}

View(ACC_cali)
write.csv(ACC_cali, "ACC_CALIFORNIA.csv")


#############################################
#######Synthetic Control#####################
#############################################
#### Cleaning for Synthetic Control ------------------

# Creating data set with State.Name State.Code
df.state = data.frame(unique(frpp_glc$State.Code),unique(frpp_glc$State.Name))
names(df.state) = c("State.Code", "State.Name")

# Performing synthetic control on data set with data > 2005
ACC_rec <- subset(ACC_AUX, YEAR > 2005)
num_years = length(unique(ACC_rec$YEAR)) #number of years

# Repeating rows for both YEAR and State.Code for having multiple years for same State
year_repeat <- data.frame(
  rep(unique(ACC_rec$YEAR), length(df.state$State.Code))
)
df2.state = df.state[rep(seq_len(nrow(df.state)),each=num_years),]
row.names(df2.state) <- c(1:nrow(df2.state))
df2.state = cbind(df2.state, year_repeat)

# Creating empty columns to store more informaiton about demographics and accidents
ACC_usa <- data.frame(df2.state, 
                      rep(0, nrow(df2.state)), 
                      rep(0, nrow(df2.state)),
                      rep(0, nrow(df2.state)),
                      rep(0, nrow(df2.state)),
                      rep(0, nrow(df2.state)))
names(ACC_usa) <- c("State.Code", "State.Name", "YEAR", "acc_total", "hr_total","pop", "income","unemp")

# Adding information about accidents and hit and run cases
for (row_index in c(1:nrow(ACC_usa))) {
  row = ACC_usa[row_index, ]
  ACC_usa[row_index, 4] <- nrow(subset(ACC_rec, YEAR == row$YEAR & STATE == row$State.Code))
  ACC_usa[row_index, 5] <- nrow(subset(ACC_rec, YEAR == row$YEAR & STATE == row$State.Code & A_HR == 1))
}
 

## Adding data about demographics (Population, Median Income and Unemployment rate) -----
# The POP_DATA.csv file contains data about median income, population and unemployment rate of each state from 2006 to 2016
# The data is in a different format than ACC_usa, so they are transposed to fit in the vector
POP_DATA <- read.csv("POP_DATA.csv", header = TRUE, stringsAsFactors = FALSE, skip = 0)
View(POP_DATA)

# Compiling columns with incomes to one single vector to add with ACC_usa dataframe
# The structure or the two data frame are different
median_income = POP_DATA[,c(3:13)]
median_vector = c()
for (row_index in c(1:nrow(median_income))) {
  for (col_index in c(1:length(median_income[row_index, ]))) {
    median_vector = c(median_vector, median_income[row_index, col_index])
  }
}
ACC_usa$income = median_vector

# Same process of compiling data from POP_DATA to ACC_usa
popn = POP_DATA[,c(25:35)]
popn_vector = c()
for (row_index in c(1:nrow(popn))) {
  for (col_index in c(1:length(popn[row_index, ]))) {
    popn_vector = c(popn_vector, popn[row_index, col_index])
  }
}
ACC_usa$pop = popn_vector

# Same process of compiling data from POP_DATA to ACC_usa
unemp = POP_DATA[,c(14:24)]
unemp_vector = c()
for (row_index in c(1:nrow(unemp))) {
  for (col_index in c(1:length(unemp[row_index, ]))) {
    unemp_vector = c(unemp_vector, unemp[row_index, col_index])
  }
}
ACC_usa$unemp = unemp_vector

write.csv(ACC_usa, "ACC_usa_data_compiled.csv")
