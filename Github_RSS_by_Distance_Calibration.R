###########################################################################################################################################################
##
##  Kristina L Paxton
##
##  April 1 2021
##
##    Code to examine the relationship between RSS values and distance for a given study area using an exponential decay model
##      -- Based on data from a test dataset of a radio transmitter at 135 known locations distributed throughout the Guam Network 
##         -- At each test location a transmitter was held stationary for a 5-minute time period 
##            -- Removed first and last minute of each test to ensure that times matched between tests and the node network
##            -- For each test, calculated an average RSS value for the 3-min middle time period individually for each node that detected the transmitter
##      -- Data for this analysis published at: https://doi.org/10.5066/P94LQWIE 
##          -- Datafile called: RSS.Localization.PaperDataset.csv 
##             -- Data associated with this test is coded as 'A' in the column 'DataSet'

##
##    1st step: Data preparation - isolates raw RSS data from node network that is associated with test data - creates the published data file - RSS.Localization.PaperDataset.csv
##    2nd step: Exponential Decay Model - uses the dataset created in step 1 to examine the relationship between RSS values and distance
##
##
##    Files Needed
##        1. TestInfo.csv: file with test information that is saved in the working directory defined below
##            - For each unique test location there should a row of information associated with each of the 3 inner minutes of that test
##                - meaning that each test will have 3 rows of information with each row identifying one of the 3 inner minutes of the test (see TestInfo_Example.csv)
##            - Columns
##                -- TagId: unique identifier of the transmitter used during the specified test
##                -- TestId: unique identifier given to each of the unique locations where a test was conducted
##                -- Date: date when the specified test was conducted
##                -- Start.Time: time when the specified test was started
##                -- End.Time: time when the specified test was ended
##                -- Min: 1-minute time period of the specified test
##                -- Hour: hour of the the specified minute of the test
##                -- TestUTMx: Easting location of the specified test 
##                -- TestUTMy: Northing location of the specified test 
##
##       2. BeepData.rds: file with the raw RSS values collected by the node network during the time period of the test. The file should be saved in the working directory defined below
##            -- output file from Import_beep.data.Github.R script (see BeepData_Example.rds)
##            -- If using a file created from another source the following columns are needed in the specified format: 
##                -- TagId: Factor identifying the unique code of a tag
##                -- Time.local: POSIXct value in the format: "2021-12-29 09:02:51" that identifies the unique datetime stamp (IN THE LOCAL TIME ZONE OF THE STUDY) when the tag was detected by the specified node
##                    ***** If you don't have a column with the local time, but only UTC time - then change lines 42-44 in Functions_RSS.Based.Localizations.R from 'Time.local' to 'Time' ************
##                    ***** BUT be sure that Dates and Times in TestInfo.csv are also in UTC time and not local time *********
##                -- NodeId: Factor identifying the node in the network that detected the specified tag
##                -- TagRSSI: Integer indicating the RSS value (in dB) of the signal of the specified tag received by the specified node
##
##       3. Nodes.csv: file with a list of all nodes in your network and their UTM locations. The file should be saved in the working directory defined below
##            - Columns
##                -- NodeId: unique identifier of given to each node 
#                 -- NodeUTMx: Easting location of the specified node
##                -- NodeUTMy: Northing location of the specified node 
##            - Other columns can also be in the file for your reference 
##            - If Node names are being converted to scientific notation in Excel open the file with a text editor (e.g. BBedit) to change names to the correct format and save the file  
## 
##       4. Functions_RSS.Based.Localizations.R - R script with functions needed to run the code below - file should be saved in the working directory defined below    
##      
##  Important Output
##      1. Calibration Dataset that contains the average RSS values for a given node that detected the signal of a test transmitter at a known location
##      2. Model coefficients from an exponential decay model that show the relationship between RSS and Distance for a node network
##            exponential model formula: avgRSS ~ a * exp(-S * distance) + K
##                  a = intercept
##                  S = decay factor
##                  K = horizontal asymptote
##
##  
##########################################################################################################################################################



# Packages needed
library(dplyr)
library(lubridate)
library(ggplot2)

# Reset R's brain - removes all previous objects
rm(list=ls())

## Set by User
# Working Directory - Provide/path/on/your/computer/where/master/csv/file/of/nodes/is/found/and/where/Functions_CTT.Network.R/is/located
working.directory <- " add here "

# Directory for Output data - Provide/path/where/you/want/output/data/to/be/stored/
outpath <- " add here - make sure it ends with a /"



## Bring in functions 
setwd(working.directory)
source("Functions_RSS.Based.Localizations.R")



## Bring in 3 Needed files - Test Information, RSS values, and Node Information - change file names in " " as needed
test.info <- read.csv("Test.Info.csv", header = T)
str(test.info) # check that data imported properly

beep.dat <- readRDS("BeepData.rds") 
str(beep.dat) # check that data imported properly

nodes <- read.csv("Nodes.csv", header = T)
str(nodes)





#################################################################
#   Step 1
# Run function to isolate raw RSS data from a node network that
# is associated with test data  
#################################################################

    ## Variables to define for function
        ## DATE.FORMAT = format of the date column in Test.Info.csv file using R standard date expressions (e.g., "%m-%d-%y", "%Y-%m-%d")
        ## TIMEZONE = Time zone where data was collected, use grep("<insert location name here>", OlsonNames(), value=TRUE) to find valid time zone name


    ## Output in R environment
        ## combined.data - dataframe that contains the average RSS values for a given node associated with each unique test


    ## Output saved
        ## RSS.Localization.PaperDataset.csv - same dataframe saved in the R environment is also saved in the folder specified by the outpath


##******* Define Variables - replace values below with user specified values *******## 
DATE.FORMAT <- "%m/%d/%y"
TIME.ZONE <- "Pacific/Guam"


# Combine RSS data from a node network with test information into a dataframe
combined.data <- data.setup(DATE.FORMAT, TIME.ZONE)





##########################################################
#    Step 2
# Exponential Decay Function to Examine Relationship 
# between distance and Tag RSS values 
##########################################################

# Use the dataset created above or bring in a file with the data
#combined.data <- read.csv("RSS.Localization.PaperDataset.csv")


## Visualize data

  # Plot of the relationship between RSS and distance
ggplot(data = combined.data, aes(x = distance, y = avgRSS, color = NodeId)) +
  geom_point(size = 2)



## Preliminary Exponential Decay Model to determine starting values for the final model 

    # SSasymp - self start for exponential model to finding the starting values of the data
    # Asym - horizontal asymptote (when large values ) - y values decay to this value 
    # R0 - numeric value when avgRSS (i.e., response variable) = 0 
    # lrc - natural logarithm of the rate constant (rate of decay)


  # Preliminary Model
exp.mod <- nls(avgRSS ~ SSasymp(distance, Asym, R0, lrc), data = combined.data)
  # Summary of Model
summary(exp.mod)
  # rate of decay
exp(coef(exp.mod)[["lrc"]])




## Exponential Decay Final Model with user provided self-starting values 
## based on visualization of the data and values in the Preliminary Model Output

    # exponential model formula: avgRSS ~ a * exp(-S * distance) + K
      # a = intercept
      # S = decay factor
      # K = horizontal asymptote


##  ***** Variables to define for final model below - replace values below with values from exp.mod ****  ## 
a <- -59
S <- 0.004561831
K <- -104.63033

  
  # Final Model
nls.mod <- nls(avgRSS ~ a * exp(-S * distance) + K, start = list(a = a, S = S, K= K), 
               data = combined.data)
  # Model Summary
summary(nls.mod)
  # Model Coefficients - *you will use these coefficient values to estimate distance based on RSS value in Github_Simulations.R*
coef(nls.mod)


## Check the fit of the model and get predicted values

  # Get residuals and fit of model and add variables to main table
combined.data$E <- residuals(nls.mod)
combined.data$fit <- fitted(nls.mod)

  # Plot residuals by fit or distance
ggplot(combined.data, aes(x = distance, y = E, color = NodeId)) +
         geom_point(size = 2)
ggplot(combined.data, aes(x = fit, y = E, color = NodeId)) +
  geom_point(size = 2)

  # Get model predictions
combined.data$pred <- predict(nls.mod)


## Save Final Dataset with Residuals and Predictions
write.csv(combined.data, paste0(outpath, "RSS.Localization.PaperDataset_withResiduals_Predictions.csv"),
          row.names = F)


## Plot with predicted line

pdf(paste0(outpath, "Relationship_Distance~RSSI.pdf"), width = 8, height = 6)

ggplot(combined.data, aes(x = distance, y = avgRSS)) + 
  geom_point() +
  geom_line(aes(y = pred), color = "#689FBB", lwd = 1.25) +
  scale_y_continuous(name = "RSS (dB)") +
  scale_x_continuous(name = "Distance (m)") +
  theme_classic()

dev.off()




