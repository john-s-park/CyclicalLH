##############

## Code to: calculate 'mu' from life history measurements in all populations

##############

# read raw life history data
clutch <- read.csv("2017_clutch_data.csv")

# packages 
library(plyr)
library(reshape2)
library(ggplot2)
library(Rmisc)


# birth date of clutch that was used to sib-mate 
clutch$clutch_1_date <- as.Date(clutch$clutch_1_date, "%m/%d/%y")
clutch$julian_1 <- julian(clutch$clutch_1_date) #add column corresponding to julian date of first clutch (origin: 1970 Jan 01)
clutch$daynight_1 <- as.numeric(as.character(mapvalues(clutch$AMPM_1, c("AM", "PM"), c("0.0","0.5")))) #convert AM/PM into numerical
clutch$timing_1 <- as.numeric(clutch$julian_1 + clutch$daynight_1) #numerical value of birth timing

# numerical values for dates when gravid females were seen from sib-mating (in order, up to five)
clutch$maturity_1_date <- as.Date(clutch$inbred_gravid_1_date, "%m/%d/%y")
clutch$maturity_julian_1 <- julian(clutch$maturity_1_date) #add column corresponding to julian date of first clutch (origin: 1970 Jan 01)
clutch$maturity_daynight_1 <- as.numeric(as.character(mapvalues(clutch$inb_AMPM_1, c("AM", "PM"), c("0.0","0.5")))) #convert AM/PM into numerical
clutch$maturity_timing_1 <- as.numeric(clutch$maturity_julian_1 + clutch$maturity_daynight_1) 

clutch$maturity_2_date <- as.Date(clutch$inbred_gravid_2_date, "%m/%d/%y")
clutch$maturity_julian_2 <- julian(clutch$maturity_2_date) #add column corresponding to julian date of first clutch (origin: 1970 Jan 01)
clutch$maturity_daynight_2 <- as.numeric(as.character(mapvalues(clutch$inb_AMPM_2, c("AM", "PM"), c("0.0","0.5")))) #convert AM/PM into numerical
clutch$maturity_timing_2 <- as.numeric(clutch$maturity_julian_2 + clutch$maturity_daynight_2) 

clutch$maturity_3_date <- as.Date(clutch$inbred_gravid_3_date, "%m/%d/%y")
clutch$maturity_julian_3 <- julian(clutch$maturity_3_date) #add column corresponding to julian date of first clutch (origin: 1970 Jan 01)
clutch$maturity_daynight_3 <- as.numeric(as.character(mapvalues(clutch$inb_AMPM_3, c("AM", "PM"), c("0.0","0.5")))) #convert AM/PM into numerical
clutch$maturity_timing_3 <- as.numeric(clutch$maturity_julian_3 + clutch$maturity_daynight_3)

clutch$maturity_4_date <- as.Date(clutch$inbred_gravid_4_date, "%m/%d/%y")
clutch$maturity_julian_4 <- julian(clutch$maturity_4_date) #add column corresponding to julian date of first clutch (origin: 1970 Jan 01)
clutch$maturity_daynight_4 <- as.numeric(as.character(mapvalues(clutch$inb_AMPM_4, c("AM", "PM"), c("0.0","0.5")))) #convert AM/PM into numerical
clutch$maturity_timing_4 <- as.numeric(clutch$maturity_julian_4 + clutch$maturity_daynight_4)

clutch$maturity_5_date <- as.Date(clutch$inbred_gravid_5_date, "%m/%d/%y")
clutch$maturity_julian_5 <- julian(clutch$maturity_5_date) #add column corresponding to julian date of first clutch (origin: 1970 Jan 01)
clutch$maturity_daynight_5 <- as.numeric(as.character(mapvalues(clutch$inb_AMPM_5, c("AM", "PM"), c("0.0","0.5")))) #convert AM/PM into numerical
clutch$maturity_timing_5 <- as.numeric(clutch$maturity_julian_5 + clutch$maturity_daynight_5)

# numerical values for time it took each female (up to five) to reach maturity
clutch$age_maturity1 <- as.numeric(clutch$maturity_timing_1 - clutch$timing_1)
clutch$age_maturity2 <- as.numeric(clutch$maturity_timing_2 - clutch$timing_1)
clutch$age_maturity3 <- as.numeric(clutch$maturity_timing_3 - clutch$timing_1)
clutch$age_maturity4 <- as.numeric(clutch$maturity_timing_4 - clutch$timing_1)
clutch$age_maturity5 <- as.numeric(clutch$maturity_timing_5 - clutch$timing_1)
clutch$age_maturity_mean <- rowMeans(clutch[,c("age_maturity1","age_maturity2","age_maturity3","age_maturity4","age_maturity5")], na.rm = TRUE)


# create dataframe of individual mean age at maturity (first five to become mature) per each mother's (N=~12) clutch 
mu_individualmeans <- as.data.frame(cbind(clutch$`Population..`, clutch$age_maturity_mean))

# turn 'age at maturity' into rate (mu) by taking inverse
mu_individualmeans$V2 <- 1/mu_individualmeans$V2
colnames(mu_individualmeans) <- c("Population", "mu")

# create csv of individual mean mu
write.csv(mu_individualmeans, "mu_individualmeans.csv")

# create csv of POPULATION MEANS of individual mean mu
mu_pops_summ <- summarySE(mu_individualmeans, measurevar = "mu", groupvars = "Population", na.rm = TRUE)
write.csv(mu_pops_summ, "mu_pops_summ.csv")
