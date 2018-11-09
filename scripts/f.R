##############

## Code to: calculate 'f' from life history measurements in all populations

##############

# read raw life history data
clutch <- read.csv("2017_clutch_data.csv")

# packages 
library(plyr)
library(reshape2)
library(ggplot2)
library(Rmisc)


# taking out outliers (as identified in later steps) -- also crosschecked with data notes for individuals that seemed to be unhealthy or otherwise unusual
clutch <- clutch[-c(24,49,103),]

# converting all check dates into julian values and adding numerical AM / PM increments
clutch$clutch_1_date <- as.Date(clutch$clutch_1_date, "%m/%d/%y")
clutch$julian_1 <- julian(clutch$clutch_1_date) #add column corresponding to julian date of first clutch (origin: 1970 Jan 01)
clutch$daynight_1 <- as.numeric(as.character(mapvalues(clutch$AMPM_1, c("AM", "PM"), c("0.0","0.5")))) #convert AM/PM into numerical
clutch$timing_1 <- as.numeric(clutch$julian_1 + clutch$daynight_1) #numerical value of first clutch timing

clutch$clutch_2_date <- as.Date(clutch$clutch_2_date, "%m/%d/%y")
clutch$julian_2 <- julian(clutch$clutch_2_date) 
clutch$daynight_2 <- as.numeric(as.character(mapvalues(clutch$AMPM_2, c("AM", "PM"), c("0.0","0.5")))) 
clutch$timing_2 <- as.numeric(clutch$julian_2 + clutch$daynight_2) #numerical value of second clutch timing

clutch$clutch_3_date <- as.Date(clutch$clutch_3_date, "%m/%d/%y")
clutch$julian_3 <- julian(clutch$clutch_3_date) 
clutch$daynight_3 <- as.numeric(as.character(mapvalues(clutch$AMPM_3, c("AM", "PM"), c("0.0","0.5")))) 
clutch$timing_3 <- as.numeric(clutch$julian_3 + clutch$daynight_3) #numerical value of third clutch timing

clutch$clutch_4_date <- as.Date(clutch$clutch_4_date, "%m/%d/%y")
clutch$julian_4 <- julian(clutch$clutch_4_date) 
clutch$daynight_4 <- as.numeric(as.character(mapvalues(clutch$AMPM_4, c("AM", "PM"), c("0.0","0.5")))) 
clutch$timing_4 <- as.numeric(clutch$julian_4 + clutch$daynight_4) #numerical value of fourth clutch timing

clutch$clutch_5_date <- as.Date(clutch$clutch_5_date, "%m/%d/%y")
clutch$julian_5 <- julian(clutch$clutch_5_date) 
clutch$daynight_5 <- as.numeric(as.character(mapvalues(clutch$AMPM_5, c("AM", "PM"), c("0.0","0.5")))) 
clutch$timing_5 <- as.numeric(clutch$julian_5 + clutch$daynight_5) #numerical value of fifth clutch timing

clutch$clutch_6_date <- as.Date(clutch$clutch_6_date, "%m/%d/%y")
clutch$julian_6 <- julian(clutch$clutch_6_date) 
clutch$daynight_6 <- as.numeric(as.character(mapvalues(clutch$AMPM_6, c("AM", "PM"), c("0.0","0.5")))) 
clutch$timing_6 <- as.numeric(clutch$julian_6 + clutch$daynight_6) #numerical value of sixth clutch timing

clutch$clutch_7_date <- as.Date(clutch$clutch_7_date, "%m/%d/%y")
clutch$julian_7 <- julian(clutch$clutch_7_date) 
clutch$daynight_7 <- as.numeric(as.character(mapvalues(clutch$AMPM_7, c("AM", "PM"), c("0.0","0.5")))) 
clutch$timing_7 <- as.numeric(clutch$julian_7 + clutch$daynight_7) #numerical value of seventh clutch timing

# adding columns for interval calculations
clutch$first_int <- as.numeric(clutch$timing_2 - clutch$timing_1)
clutch$second_int <- as.numeric(clutch$timing_3 - clutch$timing_2)
clutch$third_int <- as.numeric(clutch$timing_4 - clutch$timing_3)
clutch$fourth_int <- as.numeric(clutch$timing_5 - clutch$timing_4)
clutch$fifth_int <- as.numeric(clutch$timing_6 - clutch$timing_5)
clutch$sixth_int <- as.numeric(clutch$timing_7 - clutch$timing_6)
clutch$mean_int <- rowMeans(clutch[,c("first_int", "second_int", "third_int","fourth_int","fifth_int","sixth_int")], na.rm = TRUE)


# create dataframe of individual mean clutch interval
f_individualmeans <- as.data.frame(cbind(clutch$`Population..`, clutch$mean_int))

# turn 'clutch interval' into rate (f) by taking inverse AND multiplying by mean clutch SIZE
f_individualmeans$V2 <- 1/f_individualmeans$V2 * 47.32
colnames(f_individualmeans) <- c("Population", "f")

# create csv of individual mean f
write.csv(f_individualmeans, "f_individualmeans.csv")

# create csv of POPULATION MEANS of individual mean f
f_pops_summ <- summarySE(f_individualmeans, measurevar = "f", groupvars = "Population", na.rm = TRUE)
write.csv(f_pops_summ, "f_pops_summ.csv")
