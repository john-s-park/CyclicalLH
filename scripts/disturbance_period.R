##############

## Code to calculate mean disturbance period in each pool

##############

# function to:
# 1) take temperature timeseries and calculate derivatives, 
# 2) cluster vector of derivatives into days, 
# 3) flag days that contain derivatives that are below an arbitrary threshold (disturbance signal), and 
# 4) calculate mean intervals between disturbance days

dist_days <- function(pool_ts){
  ts_diff <- diff(pool_ts$Temp)/diff(pool_ts$X.)
  days <- split(ts_diff, ceiling(seq_along(ts_diff)/288)) # 12 * 24 = 288 time points in each day
  
  days_disturbed <- NA
  for (i in seq_along(days)){
    days_disturbed[i] <- ifelse(any(days[[i]] < mean(days[[i]]) - 7.15 * sd(days[[i]])), 1, 0)
    days_disturbed <- as.numeric(days_disturbed)
  }
  
  gaps <- c()
  c <- 1
  for (i in 1:length(days_disturbed)){
    if (days_disturbed[i] == 0) {
      c <- c + 1
    } else {
      gaps = c(gaps, c)
      c = 1
    }
  }
  
  print(mean(gaps))
}


# list of raw temperature timeseries datasets from all populations

# trim all time series so they start at midnight and end at midnight 
df_1181396 <- read.csv("1181396_editted.csv"); df_1181396 <- df_1181396[-c(1:144, 24050:nrow(df_1181396)),] 
df_20148877 <- read.csv("20148877_editted.csv"); df_20148877 <- df_20148877[-c(1:144, 24050:nrow(df_20148877)),]
df_20148879 <- read.csv("20148879_editted.csv"); df_20148879 <- df_20148879[-c(1:144, 24050:nrow(df_20148879)),]
df_2276020 <- read.csv("2276020_editted.csv"); df_2276020 <- df_2276020[-c(1:144, 24050:nrow(df_2276020)),]
df_2278765 <- read.csv("2278765_editted.csv"); df_2278765 <- df_2278765[-c(1:252, 35966:nrow(df_2278765)),]
df_2292396 <- read.csv("2292396_editted.csv"); df_2292396 <- df_2292396[-c(1:252, 35966:nrow(df_2292396)),]
df_2401944 <- read.csv("2401944_editted.csv"); df_2401944 <- df_2401944[-c(1:252, 35966:nrow(df_2401944)),]
df_906044 <- read.csv("906044_editted.csv"); df_906044 <- df_906044[-c(1:144, 24050:nrow(df_906044)),]
df_9742335 <- read.csv("9742335_editted.csv"); df_9742335 <- df_9742335[-c(1:252, 35966:nrow(df_9742335)),]
df_20148874 <- read.csv("20148874_editted.csv"); df_20148874 <- df_20148874[-c(1:144, 24050:nrow(df_20148874)),]
df_20148876 <- read.csv("20148876_editted.csv"); df_20148876 <- df_20148876[-c(1:144, 24050:nrow(df_20148876)),]
df_20148878 <- read.csv("20148878_editted.csv"); df_20148878 <- df_20148878[-c(1:144, 24050:nrow(df_20148878)),]
df_2276018 <- read.csv("2276018_editted.csv"); df_2276018 <- df_2276018[-c(1:144, 24050:nrow(df_2276018)),]
df_2276206 <- read.csv("2276206_editted.csv"); df_2276206 <- df_2276206[-c(1:252, 35966:nrow(df_2276206)),]
df_2278766 <- read.csv("2278766_editted.csv"); df_2278766 <- df_2278766[-c(1:252, 35966:nrow(df_2278766)),]
df_2382989 <- read.csv("2382989_editted.csv"); df_2382989 <- df_2382989[-c(1:252),] # logger presumably got cut short; did not trim at the tail because it almost goes up to midnight of 08/18
df_906035 <- read.csv("906035_editted.csv"); df_906035 <- df_906035[-c(1:144, 24050:nrow(df_906035)),]
df_9742332 <- read.csv("9742332_editted.csv"); df_9742332 <- df_9742332[-c(1:252, 35966:nrow(df_9742332)),]
df_9742339 <- read.csv("9742339_editted.csv"); df_9742339 <- df_9742339[-c(1:144, 24050:nrow(df_9742339)),]

# concatenate temp timeseries datasets into list
rawdf.ls <- list(df_1181396,
                 df_20148877, 
                 df_20148879, 
                 df_2276020,
                 df_2278765, 
                 df_2292396, 
                 df_2401944, 
                 df_906044, 
                 df_9742335, 
                 df_20148874, 
                 df_20148876, 
                 df_20148878, 
                 df_2276018, 
                 df_2276206, 
                 df_2278766, 
                 df_2382989, 
                 df_906035, 
                 df_9742332, 
                 df_9742339)

# create results dataframe including population number (temp logger id#'s that were used to hide regional identity during life history measurements) and region label
period.df <- data.frame(rbind(c(1181396, "FH"),
                              c(20148877, "FH"),
                              c(20148879, "FH"),
                              c(2276020, "FH"),
                              c(2278765, "NB"),
                              c(2292396, "NB"),
                              c(2401944, "NB"),
                              c(906044, "FH"),
                              c(9742335, "NB"),
                              c(20148874, "FH"),
                              c(20148876, "FH"),
                              c(20148878, "FH"),
                              c(2276018, "FH"),
                              c(2276206, "NB"),
                              c(2278766, "NB"),
                              c(2382989, "NB"),
                              c(906035, "FH"),
                              c(9742332, "NB"),
                              c(9742339, "FH") 
)
)
colnames(period.df) <- c("Population", "Region")

# create column in period.df that contains mean period per population
for (i in 1:length(rawdf.ls)){ # fill in result dataframe column with mean disturbance intervals associated with each population
  period.df$Period[i] <- dist_days(rawdf.ls[[i]])
}