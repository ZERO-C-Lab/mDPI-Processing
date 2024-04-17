#Merge identified image segments with the full frame sensor data

#This script requires: 
#1) A .csv file with the sensor data associated with each full frame image
#2) Folder of image segments with the ID as the lead portion of the file name

#The script shows how you can make vertical distribution plots and then run some further data exporation for 
#determining if some groups are correlated with relative chl-a fluorescence. This is just an example.
#The main benefit of this approach is that no julian time is needed to merge sensor data with identified image segments.
#And calculating sample volumes is easy because you know the volume of each image and the time spent in each part of the
#water column.

library(ggplot2)
library(plyr)
library(reshape2) #for making vertical distributions

#read in one required function that makes a vertical distribution
#knowing the characteristics of the mDPI (volume per image)
#also requires a plankton environment "pe" that is the nearest sensor data
#for each identified segment.
qvertplot <- function(ph, pe, bwidth) {
  bins <- seq(0, round_any(max(ph$Pressure)+(2*bwidth), bwidth, ceiling), bwidth)
  numobs <- hist(ph$Pressure, breaks = bins, plot = FALSE)$counts
  volsamp <- numobs * .102 * .077 *.25 #field of view and depth of field in m
  planknum <- ddply(pe, .(ID),
                    function(x){a <- hist(x$Pressure, breaks = bins, plot = F)$counts
                    return(a)})
  planknum <- melt(planknum)[,-2]
  numtaxa <- length(unique(planknum$ID))
  plak <- cbind(planknum, rep(bins[-1], each = numtaxa))
  plak2 <- cbind(plak, rep(volsamp[seq(1,(max(bins)/bwidth),1)], each = numtaxa))
  colnames(plak2) <- c("ID", "count", "depth", "volume.sampled")
  plak2$conc <- as.numeric(plak2$count)/as.numeric(plak2$volume.sampled)
  return(plak2)
}

#read in the frame-key sensor data
tkeyphys <- read.csv("E:/GitHub-stuff/mDPI-ProcessingCode/ExampleData/Dol_10_SAV21_16_25m_day-frames.csv", header = TRUE)

#read in a list of classified segments, make a frame key for these segments (same format as phys), merge, and plot
imglab <- list.files("E:/GitHub-stuff/mDPI-ProcessingCode/ExampleData/Station1_classified/")
imglab <- gsub(".tiff.tif", "", imglab)

#Get rid of Thumbs.db if it exists
imglab <- imglab[! imglab %in% c("data", "Thumbs.db")]
imgf2 <- unlist(strsplit(imglab, "[_]"))
imgf2 <- matrix(imgf2, ncol = 14, byrow = TRUE)
tkey <- paste0(imgf2[,8], "_", imgf2[,9], "_", imgf2[,10])

#combine the segment ID, timekey, and the full file name
biodat <- data.frame(imgf2[,1], tkey, imglab)
colnames(biodat) <- c("ID", "timekey", "filename")
colnames(tkeyphys)[1] <- "timekey"

#merge the segments to the time-key phys using the time-key
#YYYYMMDD_HHMMSSSSS_FRAME
plankenv <- merge(biodat, tkeyphys, by = "timekey", all.x = TRUE)
levels(as.factor(plankenv$ID))

alltaxa <- qvertplot(tkeyphys, plankenv, 1) #1 m vertical depth bins

#Plot the vertical distribution of all identified segments
ggplot(alltaxa, aes(-depth, conc, fill = ID))+geom_bar(stat = "identity")+coord_flip()

#Look at the thaliaceans only
thals <- subset(alltaxa, ID %in% c("dolio", "dolnurse", "dolphor", "app", "salp", "salpchain"))
ggplot(thals, aes(-depth, conc, fill = ID))+geom_bar(stat = "identity")+coord_flip()+
  labs(x = "Depth (m)", y = expression(paste("Abundance"~"(ind."~m^{-3},")")))+theme_bw()

#which zoops are the most common in the transect?
totvol <- nrow(tkeyphys)*.102 * .077 *.25
ddply(plankenv, .(ID), summarise, conc = length(ID)/totvol)

#look at the hydros, chaetos, and ctenos only (common zoops)
czoops <- subset(alltaxa, ID %in% c("cteno", "chaeto", "hydro"))
ggplot(czoops, aes(-depth, conc, fill = ID))+geom_bar(stat = "identity")+coord_flip()+
  labs(x = "Depth (m)", y = expression(paste("Abundance"~"(ind."~m^{-3},")")))+theme_bw()

#where were the marine snow and the diatoms and tricho?
snow <- subset(alltaxa, ID %in% c("snow", "diatom", "tricho"))
ggplot(snow, aes(-depth, conc, fill = ID))+geom_bar(stat = "identity")+coord_flip()+
  labs(x = "Depth (m)", y = expression(paste("Abundance"~"(ind."~m^{-3},")")))+theme_bw()

#shrimp and fish were somewhat common
fishy <- subset(alltaxa, ID %in% c("shrimp", "stoma", "fish"))
ggplot(fishy, aes(-depth, conc, fill = ID))+geom_bar(stat = "identity")+coord_flip()+
  labs(x = "Depth (m)", y = expression(paste("Abundance"~"(ind."~m^{-3},")")))+theme_bw()

#does the marine snow or diatom abundance correlate with chl-a averaged in 0.5 m bins?
tkeyphys$rdepth <- round_any(tkeyphys$Pressure, 0.5, ceiling)
chlaround <- ddply(tkeyphys, .(rdepth), summarise, mchla = mean(ChlaCounts))

diat <- subset(plankenv, ID %in% c("snow", "diatom"))
diavert <- qvertplot(tkeyphys, diat, 0.5)
diavert <- subset(diavert, volume.sampled > 0)


#combine
diaconc <- subset(diavert, ID == "diatom")[,5]
snowconc <- subset(diavert, ID == "snow")[,5]

nozoop <- data.frame(chlaround, diaconc, snowconc)

mnz <- melt(nozoop, id = c("rdepth", "mchla"))

ggplot(mnz, aes(mchla, value))+geom_point()+facet_wrap(~variable)+geom_smooth(method = "lm")
 
#what about both diatoms and snow abundances together?
nozoop$both <- nozoop$diaconc + nozoop$snowconc
#remake the melted data frame and re plot
mnz <- melt(nozoop, id = c("rdepth", "mchla"))
ggplot(mnz, aes(mchla, value, colour = -rdepth))+geom_point()+facet_wrap(~variable)+geom_smooth(method = "lm")+
  labs(x = "Mean Chl-a (counts)", y = expression(paste('Abundance'~'(ind.'~m^{-3},')')), colour = "Depth (m)")

#slope is steepest with both added together, but the deeper waters have high chl-a (relationship appears nonlinear). 
#Deep waters have low abundances of these large phyto or detritus  - so 
#perhaps there are smaller phytoplankton down there that are fluorescing?