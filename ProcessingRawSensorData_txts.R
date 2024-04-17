#TROUBLESHOOTING:
#If files fail to load properly using ldply(), check raw data and delete unusual characters.
#The CTD.txt files routinely write characters in rows that prevent as.numeric() from working

#Dol01 --> Dol11 use the old .txt raw data for sensors

#   STEPS      #
# 1) fill out metadata: cruise, rvsav, station, diurnal, wd, and intended save path 
# 2) Read in required packages and functions
# 3) Read in different sensor data and combine to the highest resolution sensor (CTD)
# 4) Optional cleaning of sensor data at the end (i.e., remove erroneous values)

library(tidyverse) #Code organization and pipe operators
library(readr) #File loading and importing
library(plyr)  #Array manipulation
library(lubridate) #Makes dealing with date-time data easier
library(magrittr) #reverse pipe operator in tidyverse system %<>% - overwrites the initial object

cruise <- 'Dol_10'
rvsav <- 'SAV21_16'
tow_day <- 'Aug_10'
station <- '25m'
diurnal <- "day"
wd <-"E:/GitHub-stuff/mDPI-ProcessingCode/ExampleData/Station1-25m/"
save_path <- "E:/test-process/"
setwd(wd)

cruise_name <- paste0(cruise,'_',rvsav)

#just saving some strings for use in the graphs and file-saving later, using the input strings above
tow_save_name <-paste0(cruise_name,"_",station,"_",diurnal)
tow_name <- paste0(tow_day,': ',tow_save_name)

#modified julianmaker, no slices needed for physical data, use seconds to a decimal place
julianmaker <- function (t) {
  jul <- as.character(t)
  hr <- substr(jul, 1,2)
  hr <- as.numeric(hr)
  minu <- as.numeric(substr(jul, 4,5))
  sec <- as.numeric(substr(jul, 7,12)) 
  julian <- (hr/24)+(minu/(24*60))+(sec/(24*60*60))
  return(julian)
}

#Make an index of the nearest time point for merging
nearphys <- function(x,y) {
  n <- length(x) 
  phystime <- vector(length = n)
  for(i in 1:n) {
    phystime[i] <- which.min(abs(x[i]-y))
  }
  return(phystime)
}

#function to convert micromole liter from optode to mg liter
O2mol_to_mg <- function(O2umoles){
  O2moles = O2umoles/1000000
  O2mg= as.numeric((O2moles*32000))
  return(O2mg)
}

#-------------------------------------------------------------------------------------------------
#CTD data loading
CTDfiles <- list.files(wd,pattern=glob2rx("CTD*.txt"))

CTDdf <- ldply(CTDfiles,read_csv2, skip = 1) #keeps the fractions of a second
head(CTDdf)

CTDdf <- as.data.frame(str_split_fixed(CTDdf[,1], ",", 7))
head(CTDdf)
#tail(CTDdf) #bunch of junk at the bottom from blank data reads

colnames(CTDdf) <- c("Date", "Time", "Temp", "Conductivity", "Pressure",
                     "Salinity", "Velocity")

#need to make sure columns 3-7 are numeric, not characters
for (i in 3:7){
  CTDdf[,i] <- as.numeric(as.character(CTDdf[,i]))
}

CTDdf_naexclude <- na.exclude(CTDdf) #excludes NA values from data-frame

#get rid of NAs and convert Time to julian
CTDdf_cleaned <- CTDdf_naexclude %>% 
  mutate(julian = julianmaker(Time)) 

#-------------------------------------------------------------------------------------------------
#Chl-a data loading
Fluorfiles <- list.files(wd,pattern=glob2rx("ECO*.txt"))

Chla_df <- ldply(Fluorfiles,read.delim,skip = 1)
#Have to remove 2 columns before na.exclude as the ECO-FL produces
#2 columns with basic sensor information
Chla_df <- Chla_df[,c(-3,-4)]

#rename the chl-a counts column then remove the wavelength and signal columns
colnames(Chla_df)[4] <- "ChlaCounts" 
Chla_df %<>% select(-c('Wavelength','Signal'))

Chla_df$ChlaCounts <- as.numeric(as.character(Chla_df$ChlaCounts))

Chla_df_cleaned <- Chla_df %>% 
  na.exclude(ChlaCounts) %>%
  mutate(julian = julianmaker(Time))

#-------------------------------------------------------------------------------------------------
#Dissolved oxygen data loading
Optodefiles <- list.files(wd,pattern=glob2rx("Opt*.txt"))
Opt_df <- ldply(Optodefiles, read.delim, skip = 2,  sep = "\t", header = FALSE)

colnames(Opt_df) <- c("Date","Time","Model", "Serial Number", "DOmolar","AirSat","Temp","raw data1",
                      "raw data2", "raw data3", "raw data4","raw data5","raw data6","raw data7")
#remove junk columns
Opt_df <- Opt_df[,c(-3, -4, -8:-14)]

#make sure all relevant data are numeric
Opt_df$DOmolar <- as.numeric(Opt_df$DOmolar)
Opt_df$AirSat <- as.numeric(Opt_df$AirSat)
Opt_df$Temp <- as.numeric(Opt_df$Temp)

#remove NA values from Optode data
Opt_df_cleaned <- Opt_df %>% 
  na.exclude(Temp) %>% na.exclude(AirSat) %>% na.exclude(DOmolar)

Opt_df_cleaned$julian <- julianmaker(Opt_df_cleaned$Time)
head(Opt_df_cleaned)

#----------------------------------------------------------------
#Using nearphys function to merge data from the 3 sensors.
#Must do once to merge CTD abd chla and again to merge with oxygen optode.

#merging CTD and Chla
mergephys1 <- nearphys(CTDdf_cleaned$julian, Chla_df_cleaned$julian)


nearjul <- Chla_df_cleaned[mergephys1,4] #pull the julians
CTD2 <- cbind(nearjul, CTDdf_cleaned)
colnames(Chla_df_cleaned)[4] <- "nearjul" #change the name of julian to nearjul for merger

CTDChla <- merge(CTD2, Chla_df_cleaned, by = "nearjul", all.x = TRUE)
head(CTDChla)

#merging CTD and chla with oxygen Optode
mergephys2 <- nearphys(CTDChla$julian, Opt_df_cleaned$julian)

nearjul2 <- Opt_df_cleaned[mergephys2, 6] #pull the julians
CTD3 <- cbind(nearjul2, CTDChla)
colnames(Opt_df_cleaned)[6] <- "nearjul2" #change the name of julian to nearjul for merger

CTDOpt <- merge(CTD3, Opt_df_cleaned, by = "nearjul2", all.x = TRUE)
head(CTDOpt)
allsensors <- CTDOpt

#Rename some columns and make calculations
allsensors  %<>% dplyr::rename(Temp.Opt = 'Temp.y') %<>% dplyr::rename(Temp.CTD = 'Temp.x')%<>% 
  mutate(O2mgL = O2mol_to_mg(DOmolar)) %<>%
  mutate(Station = as.factor(station)) %<>% mutate(TOD = as.factor(diurnal)) %<>%
  select(-DOmolar)

#Ged rid of the useless columns that were created for the merging
drop.cols <- c('nearjul2','nearjul','Date.y','Time.y','Date.x','Time')
allsensors_cleaned <- allsensors %>% select(-drop.cols) %>% dplyr::rename(Time = 'Time.x') %>% mutate(Cruise = cruise)

write.csv(allsensors_cleaned,paste0(save_path,tow_save_name,'.csv'), row.names = FALSE)

#can filter to remove obviously erroneous values and drop variables 
allsensors_cleaned %<>%
  filter(Salinity>15) %<>% filter(Pressure>1) %<>% filter(Temp.CTD > 3) %<>% filter(Conductivity > 0)
