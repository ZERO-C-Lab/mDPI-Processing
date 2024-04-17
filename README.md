# mDPI-Processing
 Code and example data for processing mDPI sensor data and plotting distributions of zooplankton

This repository includes code for:
1) Reading in raw sensor data on the mDPI compact vehicle, in .txt or .csv format, and combining into one large .csv file (ProcessingRawSensorData_txts.R or ProcessingRawSensorData_csvs.R).
2) Matching sensor or physical oceanographic data from the large .csv file to a list of full frame images from the mDPI. This process generates oceanographic variables associated with each image (mDPI_FullFramePhys.R).
3) Merging identified image segments (i.e., zooplankton) to their nearest sensor data and plotting vertical distributions (mDPI_VerticalDistributionPlots.R).
4) A .pdf that explains how the mDPI works and goes through the data processing: from the large sensor data .csv file to plotting vertical distributions of zooplankton.

The folder ExampleData.zip has the files needed to run a test and look at the data from one station sampled on August 10, 2021. 
These example data can be downloaded at https://outlookuga-my.sharepoint.com/:u:/g/personal/atgreer_uga_edu/EbSmr0JAAXxJrpWok1-Ca6cBDg5OQPq_U_vFWHDX9IzBUA?e=7GQzKX 

