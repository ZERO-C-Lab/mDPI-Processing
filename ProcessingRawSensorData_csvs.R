##Used for processing raw sensor data from the mDPI using .csv files (Dol_13 and later for DolMICROBE)

#   STEPS      #
# 1) fill out metadata: cruise, rvsav, station, diurnal, wd, and intended save path 
# 2) Read in required packages and functions
# 3) Read in different sensor data and combine to the highest resolution sensor (CTD)
# 4) Optional plotting throughout to check results and can clean sensor data at the end (i.e., remove erroneous values)

library(tidyverse) #pipelines
library(magrittr)  #Reverse pipeline operators %<>% - will overwrite the object
library(readr)     #file loading/reading
library(plyr)      #dataframe manipulation
library(lubridate) #date-time data handling functions 

cruise <- 'Dol_20'
rvsav <- 'SAV22_16'
station <- '45m'
diurnal <- "day"
wd <- "E:/GitHub-stuff/mDPI-ProcessingCode/ExampleData/45m/" #Directory of sensor data .csv files
save_path <- "E:/test-process/" #Directory for the processed data output
setwd(wd) 

cruise_name <- paste0(cruise,"_",rvsav)
#just saving some strings for use in the graphs and file-saving later, using the input strings above

tow_save_name <- paste0(cruise_name,"_",station,"_",diurnal)

#Merge the different sensors using the unix time - then generate a julian time at the end that 
#you can use to associate with the images from the mDPI
#-------------------------------------------------------------------------------------
#functions needed

#Different than julianmaker for the .txt output files from old mDPI due to data formatting
#modified julianmaker
#slight change to the function to make hhmmss longer to include the extra digits
julianmaker_unix <- function (t) {
  datetime <- as.POSIXlt(t, origin="1970-01-01", tz = "EST")
  hhmmss <- substr(datetime,12,24)
  hr <- substr(hhmmss, 1,2)
  hr <- as.numeric(hr)
  minu <- as.numeric(substr(hhmmss, 4,5))
  sec <- as.numeric(substr(hhmmss, 7,13)) 
  julian <- (hr/24)+(minu/(24*60))+(sec/(24*60*60))
  return(julian)
}

options(digits.secs=4) #need this to avoid repeated julians

#Make an index of the nearest time point for merging
nearphys <- function(x,y) {
  n <- length(x) 
  phystime <- vector(length = n)
  for(i in 1:n) {
    phystime[i] <- which.min(abs(x[i]-y))
  }
  return(phystime)
}

#function to convert micromole liter from optode to mg liter at the end
O2mol_to_mg <- function(O2umoles){
  O2moles = O2umoles/1000000
  O2mg= as.numeric((O2moles*32000))
  return(O2mg)
}

#Use at the end for the final data frame
hhmmss_fromUnixTime <- function(t){
  datetime <- as.POSIXlt(t, origin="1970-01-01", tz = "EST")
  hhmmss <- substr(datetime,12,19)
  return(hhmmss)
}


##-----------------------------------------------------------------------------------
#reading in CTD data
CTDfiles <- list.files(wd,pattern=glob2rx("CTD*.csv"))

#We choose to use read.delim as out input function due to issues in read.csv skipping numbers at the end of decimals
CTDdf <- ldply(CTDfiles,read.delim ,skip =1,  sep = ",", header=FALSE) #includes 1000ths of a second - 

#Rename the columns from the CTD files
colnames(CTDdf) <- c("Time","Temp","Cond","Pressure","Depth","Salinity","Sound Velocity")
#All columns should be numeric

#Exclude the NA values
CTDdf_cleaned <- na.exclude(CTDdf)
#head(CTDdf_cleaned)

#-------------------------------------------------------------------------------------
#reading in chlorophyll data and cleaning it
Fluorfiles <- list.files(wd,pattern=glob2rx("Fluo*.csv"))

Chla_df <- ldply(Fluorfiles, read.delim, skip = 1, sep = ",", header = FALSE)
#Renaming the columns
colnames(Chla_df) <- c("Time", "Wavelength", "ChlaCounts", "Signal")

Chla_df <- Chla_df[,c(-2,-4)] #Removing the Wavelength and signal columns
#head(Chla_df) #check and see if there are 2 columns - time and chl-a counts

Chla_df$ChlaCounts <- as.numeric(Chla_df$ChlaCounts)
Chla_df_cleaned <- Chla_df %>% na.exclude(ChlaCounts)
#check
#head(Chla_df_cleaned)

#-------------------------------------------------------------------------------------------------
#oxygen data loading
Optodefiles <- list.files(wd,pattern=glob2rx("Oxy*.csv"))
Opt_df <- ldply(Optodefiles, read.delim, skip = 0,  sep = ",", header = TRUE)

#readjusting column names 
colnames(Opt_df) <- c("Time","Model", "Serial Number", "DOmolar","AirSat","Temp","raw data1",
                      "raw data2", "raw data3", "raw data4","raw data5","raw data6","raw data7")
#removing some junk columns
Opt_df <- Opt_df[,c(-2, -3, -7:-13)]

#Making sure all the proper columns are formatted as numeric data
Opt_df$DOmolar <- as.numeric(Opt_df$DOmolar)
Opt_df$AirSat <- as.numeric(Opt_df$AirSat)
Opt_df$Temp <- as.numeric(Opt_df$Temp)

Opt_df_cleaned <- Opt_df %>% 
  na.exclude(Temp) %>% na.exclude(AirSat) %>% na.exclude(DOmolar)

#MERGING
#----------------------------------------------------------------
#using nearphys function to merge data from the 3 sensors, must do once to merge CTD and chla 
#and again to merge with oxygen optode

#merging CTD and Chla - find where they matched unix time
mergephys1 <- nearphys(CTDdf_cleaned$Time, Chla_df_cleaned$Time)
#Creating a new array of 'nearest points for each row'
nearunix <- Chla_df_cleaned[mergephys1,1] #pull the unix Time column (1st one)
#adding the nearjul array as a column onto the CTD dataframe before the merge
#The nearestjul is the nearest julian point on the CTD dataframe for each Chla-df row
CTD2 <- cbind(nearunix, CTDdf_cleaned)
colnames(Chla_df_cleaned)[1] <- "nearunix" #change the name of julian to nearjul for merger

#Making the merge b/w CTD and Chla dataframes via the nearest julian-timestamped sample
CTDChla <- merge(CTD2, Chla_df_cleaned, by = "nearunix", all.x = TRUE)
head(CTDChla)

###check
#ggplot(CTDChla, aes(Temp, -Depth, colour = ChlaCounts))+geom_point() #looks good!

#merging CTD and chla with oxygen Optode
mergephys2 <- nearphys(CTDChla$nearunix, Opt_df_cleaned$Time)
nearunix2 <- Opt_df_cleaned[mergephys2,1] #pull the unix Time
CTD3 <- cbind(nearunix2, CTDChla)
colnames(Opt_df_cleaned)[1] <- "nearunix2" #change the name of Time to nearunix2 for merger
#Making the merge b/w CTD+Chla and Optode
CTDOpt <- merge(CTD3, Opt_df_cleaned, by = "nearunix2", all.x = TRUE)
head(CTDOpt)

#can plot to check
#ggplot(CTDOpt, aes(Temp.x, -Depth, colour = DOmolar))+geom_point() #looks good!

#drop the first 2 columns then make a julian time
CTDOpt <- CTDOpt[,-1:-2]
CTDOpt$julian <- julianmaker_unix(CTDOpt$Time)
CTDOpt <- CTDOpt[,-11] #get rid of the Temp.y (the 11th column) which is from the o2 sensor, and we do not use
#head(CTDOpt)

#Clean up the allsensors dataframe
#Remove junk columns (nearjuls, duplicate Times/date/temps)

allsensors <- CTDOpt

#rename some of the columns
allsensors %<>% dplyr::rename(Temp.CTD = 'Temp.x') %<>%
  mutate(O2mgL = O2mol_to_mg(DOmolar)) %<>% dplyr::rename(Conductivity ='Cond') %<>%
  mutate(Station = as.factor(station)) %<>% mutate(TOD = as.factor(diurnal)) %<>%
  mutate(Cruise =as.factor(cruise)) %<>%
  mutate(Time = hhmmss_fromUnixTime(Time)) %<>% dplyr::rename(S.Velocity ='Sound Velocity')

#write your new data frame to the disk - will use the information you gave it
#at the beginning to make the file name
write.csv(allsensors, paste0(save_path,tow_save_name,'.csv'), row.names = FALSE)

#can filter to remove obviously erroneous values and drop variables so it fits with 
#sensor data collected in the past (only pressure was recorded, for example)
allsensors %<>%
  filter(Salinity>15) %<>% filter(Pressure>1) %<>% filter(Temp.CTD > 3) %<>% filter(Conductivity > 0) %<>%
  select(-c(DOmolar,Depth))


