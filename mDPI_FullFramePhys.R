#Match sensor data to each full frame image from the mDPI

#This script requires: 
#1) A directory of full frame mDPI images
#2) A .csv file of sensor data with a julian time for each row

#The output is a .csv file with one set of sensor data for each full frame image
#The time key is a unique identifier that can be used to link any segment or full frame image to
#sensor data collected on the mDPI compact vehicle. Seconds are to the thousandths
#Format of timekey: YYYYMMDD_HHMMSSSSS_FRAME
#Example: 20210810_143626523_1 = first image taken on August 10, 2021 at 14:36:26.523 local time (Eastern US) 

#required functions to make julian time and find nearest julian time between 2 data frames
julianmakerPD <- function (t,s) {
  if(missing(s)){
    #Handles time without any image number (slices) in the format hh:mm:ss.sss
    #This input without slices must be a string
    jul <- as.character(t)
    hr <- as.numeric(substr(jul, 1,2))
    minu <- as.numeric(substr(jul, 4,5))
    sec <- as.numeric(substr(jul, 7,12))
    julian <- (hr/24)+(minu/(24*60))+(sec/(24*60*60))
  } else {
    #Handles the timestamps from the ROI which are formatted without decimals
    #Uses the image number (slice) to add seconds onto julian calculation
    jul <- as.character(t)
    hr <- as.numeric(substr(jul, 1,2))
    minu <- as.numeric(substr(jul, 3,4))
    sec <- as.numeric(paste0(substr(jul, 5,6), ".", substr(jul,7,9)))
    #4.924 images per second, so this must be added to total seconds
    julian <- (hr/24)+(minu/(24*60))+((sec+(s/4.924))/(24*60*60))
  }
  return(julian)
}

nearphys <- function(x,y) {
  n <- length(x) #x will be the longer data frame
  phystime <- vector(length = n)
  for(i in 1:n) {
    phystime[i] <- which.min(abs(x[i]-y))
  }
  return(phystime)
}

#Read in the .csv file with all of the sensor data combined
#Better to not filter the sensor data so that every image has corresponding sensor data
phys <- read.csv("e:/test-process/Dol_10_SAV21_16_25m_day.csv", header = TRUE)

#Read in the list of full frame file names
#frame_list <-  list.files("Z:/2020_DolMicrobe/Cruise_Data/SAV 21-16/Station1-25m/")
#Read from a .csv file in the example data
frame_list <- read.csv("E:/GitHub-stuff/mDPI-ProcessingCode/ExampleData/Station1-ListofFullFrames.csv", header = TRUE)

#Get rid any weird files or folders in in there with the image files
#And get rid of .tiff in the strings
frame_list <- frame_list[!frame_list %in% c("data","Thumbs.db")]
#frame_list <- gsub(".tiff", "", frame_list)
frame_list <- gsub(".tiff", "", frame_list[,1])

#split by the underscore and create a matrix
#with the values read in rows for 8 columns
imgf <- unlist(strsplit(frame_list, "[_]"))
imgf2 <- matrix(imgf, ncol = 8, byrow = TRUE)
imgf2 <- as.data.frame(imgf2, stringsAsFactors = FALSE)

#isolate the time component and the frame number
ftime <- data.frame(imgf2$V6, imgf2$V7, imgf2$V8)
colnames(ftime) <- c("date","timestamp", "frame")

#keep timestamp as character to avoid dropping a leading 0
ftime$timestamp <- as.character(ftime$timestamp)
ftime$frame <- as.numeric(as.character(ftime$frame))

#make the julian time using timestamp and frame #
jul <- julianmakerPD(ftime$timestamp, ftime$frame)

#make the time key
#and attach the corresponding julian time
tkey <- paste0(ftime$date,"_", ftime$timestamp, "_", ftime$frame)
nf <- data.frame(jul, tkey)

#merge each tkey to the nearest data point from phys
ind1 <- nearphys(nf$jul, phys$julian)
julph <- phys$julian[ind1]
nf2 <- cbind(julph, nf)
colnames(nf2)[1] <- "julian"
head(nf2)

tkeyphys <- merge(nf2, phys, by = "julian", all.x = TRUE)

#get rid of the 2 julian columns and save
tkeyphys <- tkeyphys[,c(-1, -2)]
write.csv(tkeyphys, file = "e:/test-process/Dol_10_SAV21_16_25m_day-frames.csv", row.names = FALSE)

#Once you have this .csv file with the sensor data for each image, it is very easy to use the first column (tkey) to
#merge with any segment or ROI or particle detected in that frame to calculate abundances or vertical distributions.
