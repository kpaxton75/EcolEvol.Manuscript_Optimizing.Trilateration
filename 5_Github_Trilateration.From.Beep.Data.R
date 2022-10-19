###############################################################################################################################################
# Kristina L Paxton
#  March 29 2022
#
#   Code to estimate locations of tagged individuals based on Trilateration
#
#
#  Files Required to be in the working directory
#
#  1. Beep data -- exported from Import_beep.data.Github.R and save in the working directory indicated below
      # Name and type of file - RDS file named beep_data_START.DATE_END.DATE.rds (e.g., beep_data_2022-03-26_2022-03-27.rds)
          # START.DATE and END.DATE are replaced with the start and end date of the beep data
          # Columns 
              # SensorId - name of the Sensor Station where the data was gathered
              # Time - GMT time of the detection
              # RadioId - unique identifier of radio port
              # TagId - unique identifier of the tag detected
              # TagRSSI - signal strenght of the tag detected in dB
              # NodeId - unique identifier of the node that detected the signal
              # Validated - CTT field for Hybrid and LifeTags
              # Time.local - time for the local time zone defined by user
              # v - version of the CTT node

# 2. Node data - File that contains a list of all nodes in your network that is saved in the working directory defined below
      # Type of file - CSV file (ex. Nodes_trilateration.csv)
          # Columns REQUIRED - MUST BE IN THE FOLLOWING ORDER: 
              # NodeId - unique identifier of each node in your network
              # UTMx - Easting location of the specified node
              # UTMy - Northing location of the specified node
              # Other columns may be added for your reference after the 3 specified columns

# 3. Tag data - file tat contains a list of all tags that you want to calculate location estimates based on trilateration
      # Type of file - CSV file (ex. Tags_trilateration.csv)
          # Columns required:
              # TagId - unique identifier of the tag detected
              # StartDate - date when tag was deployed

# 4. Functions_RSS.Based.Localizations.R - R script containing functions to run the trilateration script
#
#
######################################################################################################################################
  

### packages needed
library(tidyverse) # organization of data
library(lubridate) # dates
library(slider) # sliding window
library(padr) # thicken function
library(nlstools) # trilateration
library(raster) # calculate distance between nodes



### Reset R's brain - removes all previous objects
rm(list=ls())



## Set by User 
  # Working Directory  - Provide/path/on/your/computer/where/master/csv/file/of/nodes/is/found/and/where/Functions_CTT.Network.R/is/located
working.directory <- " add here "

# Directory for Output data - Provide/path/where/you/want/output/data/to/be/stored/
outpath <- " add here - make sure it ends with a /"



## Bring in Functions
setwd(working.directory)
source("Functions_RSS.Based.Localizations.R")





#################################################################
#   Step 1
#
# Bring in 3 files needed 
#   - Beep Data, Node Information, Tag Information
#################################################################

    ## Variables to define
        ## DATE.FORMAT = format of the date column in Beep Data file using R standard date expressions (e.g., "%m-%d-%y", "%Y/%m/%d")
        ## TIMEZONE = Time zone where data was collected, use grep("<insert location name here>", OlsonNames(), value=TRUE) to find valid time zone name
        ## START = Start date of the beep data in format YYYY-MM-DD
        ## END = End date of the beep data in format YYYY-MM-DD

    ## Output in R environment
        ## beep.data - data frame that contains beep data for the specified time period
        ## nodes - data frame that contains information about the nodes in your network 
        ## tags - data frame that contains information about the tags that you want to calculate location data

    

    ##******* Define Variables - replace values below with user specified values *******## 

TIMEZONE <- "Pacific/Guam"
DATE.FORMAT <- "%m/%d/%y"
START <- "2022-03-26"
END <- "2022-03-27"


    
    ## Bring in 3 files needed for functions below and format date columns where needed
beep.data <- readRDS(paste0("beep_data_", START, "_", END, ".rds")) %>%
  dplyr::mutate(Date = as.Date.POSIXct(Time.local, tz = TIMEZONE))
str(beep.data) # check that data imported properly


nodes <- read.csv("Nodes_trilateration.csv", header = T)
str(nodes) 


tags <- read.csv("Tags_trilateration.csv", header = T) %>%
  dplyr::mutate(StartDate = as.Date(StartDate, format = DATE.FORMAT))
str(tags)



######################################################################################
# Step 2
#   Function to filter beep data that occurs prior to the Start Date of each TagId
#     -- If you do not want to filter beep data then skip 
#        this step and replace input file 'beep.filtered' with 'beep.data' in Step 3
######################################################################################


beep.filtered <- filter.dates(beep.data,tags)





###################################################################################################
#   Step 3
#
#  Prepare data for trilateration
#     -- A sliding window will be applied to RSSI values for a specified number of minutes  
#        for each TagId and NodeId pair to help remove the effects of aberrant RSSI values
#     -- RSSI values for a given TagId and NodeId will be averaged across a specified time period
#     -- Estimate distance of a tag signal based on RSSI value and the exponential decay between
#         RSS values and distance
#
#   ** This step can take a long time depending on the number of observations in beep.data file and
#      the speed of your computer -- recommend that you test script with a small data set **
#
######################################################################################################

    ## Variables to define for function
        ## Coefficients in the exponential model: use output from Github_RSS_by_Distance_Calibration (e.g., coef(nls.mod))
              ## a = intercept
              ## S = decay factor
              ## K = horizontal asymptote
        ## SLIDE.TIME = number of minutes before or after an RSS value to include in the sliding window 
              ## If you do not want to run a sliding window across RSS values indicate 0
        ## GROUP.TIME = number of minutes to average RSS values across
        
    ## Output in R environment
        ## beep.grouped - dataframe that contains average RSSI value for a given time period (e.g., 1 min) along with the
                          ## estimated distance of the signal from the receiving node
            ## Columns:
                ## TagId - unique identifier of the transmitter
                ## NodeId - unique identifier of the node receiving a signal from the transmitter
                ## Time.group - Date and start time of the 1-min time period (or other specified time period) when transmitter signal was received
                ## Date - Date when signal was received
                ## mean_rssi - average RSSI value across a 1-min time period (or other specified time period) for a given transmitter and node
                ## beep_count - number of signals received by the node during the 1-min time period (or other specified time period)
                ## e.dist - estimated distance (in meters) between the specified node and the transmitter signal
                ## UTMx - Easting location of the specified node
                ## UTMy - Northing location of the specified node



  ##******* Define Variables - replace values below with user specified values *******##

a <- 47.23 
S <- 0.005
K <- -105.16
SLIDE.TIME <- 2
GROUP.TIME <- "1 min"


# Function to prepare beep data for trilateration 
# by estimating distance of a signal based on RSS values
beep.grouped <- prep.data(beep.filtered,nodes) 






########################################################################
#   Step 4
#
#   Estimate locations of tagged individuals based on trilateration 
#
########################################################################


    ## Variables to define for function
        ## DIST.filter = distances (in meters) that will be used to filter data prior to trilateration analysis
              # For each time period the node with the strongest signal is identified and then only nodes
              # within the specified distance from the strongest node are kept for trilateration analysis
                  ## If you do not want to filter by distance then select a distance that is greater than the
                  ## largest distance between nodes in your network (e.g., 2000)
        ## RSS.FILTER = RSS value that will be used to filter data prior to trilateration analysis
              # All RSS values less than the specified value will be removed from the beep dataset
                  ## If you do not want to filter by RSS values then select a value greater than your horizontal asymptope, K (e.g., -106)

  ## Output
      ## While the function is running it will indicate which tag it is processing and the number of unique
      ## time periods were a location will be estimated from trilateration

      ## Location estimates for the specified tags within the time period of the beep data will be:
          ## Saved as an RDS file in the outpath specified with the name 'Estimated.Locations_StartDate_EndDate.rds' (e.g., Estimated.Locations_2022-03-26_2022-03_27.rds)
          ## Dataframe in the R Environment named 'location.estimates'
              ## Columns
                  # TagId - unique identifier of the transmitter
                  # Time.group - Date and start time of the 1-min time period (or other specified time period) when transmitter signal was received
                  # Hour - Hour of the day (in military time) when the transmitter signal was received
                  # No.Nodes - number of nodes used to estimate the location
                  # UTMx_est = Estimated Easting location
                  # UTMy_est = Estimated Northing location
                  # x.LCI = Lower 95% confidence interval around the estimated Easting location
                  # x.UCI = Upper 95% confidence interval around the estimated Easting location
                  # y.LCI = Lower 95% confidence interval around the estimated Northing location
                  # y.UCI = Upper 95% confidence interval around the estimated Northing location
      



##******* Define Variables - replace values below with user specified values *******##

DIST.filter <- 350
RSS.filter <- -90


location.estimates <- trilateration(beep.grouped)





