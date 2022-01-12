###########################################################################################################################################################
##
##  Kristina L Paxton
##
##  April 1 2021
##
##    -- Code to assess the localization error associated with RSS-based localization estimates for a node network based on a test dataset
##        -- Test dataset should be locations that differ from the dataset used to calibrate the relationship between RSS and distance for the node network
##        -- localization error = difference in distance (m) between true and estimated location
##
##    -- Analysis in Ecology and Evolution Paper is based on data from a test with a radio transmitter at 54 known locations randomly distributed throughout the Guam Network 
##        -- At each test location a transmitter was held stationary for a 5-minute time period 
##            -- Removed first and last minute of each test to ensure that times matched between tests and the node network
##            -- For each test, calculated an average RSS value for the 3-min middle time period individually for each node that detected the transmitter
##
##
##    1st step: Data preparation - isolates raw RSS data from node network that is associated with test data 
##       -- Creates the data file that was published with this paper at: https://doi.org/10.5066/P94LQWIE 
##            -- Data associated with this test is coded as 'B' in the column 'DataSet'
##    2nd step: Estimate the distance of the test signal from each node that detected the signal based on the exponential decay relationship between RSS and Distance 
##    3rd step: Calculate the error associated with RSS-based localization estimates when filters are applied prior to trilateration analysis
##
##
##    Files Needed
##        Files Needed
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
##            -- Columns
##                -- NodeId: unique identifier of given to each node 
#                 -- NodeUTMx: Easting location of the specified node
##                -- NodeUTMy: Northing location of the specified node 
##            -- Other columns can also be in the file for your reference 
##            -- If Node names are being converted to scientific notation in Excel open the file with a text editor (e.g. BBedit) to change names to the correct format and save the file  
##
##       4. Model coefficients from exponential decay model that describes the relationship between RSS and Distance for the node network where test data was collected
##            -- values are from final model calculated in: Github_RSS-by_Distance_Calibration.R script
##            -- exponential model formula: avgRSS ~ a * exp(-S * distance) + K
##                     -- a = intercept
##                     -- S = decay factor
##                     -- K = horizontal asymptote 
##
##       5. Functions_RSS.Based.Localizations.R - R script with functions needed to run the code below - file should be saved in the working directory defined below    
##
##    Important Output
##        1. Dataframes with summary statistics of localization error when different filters (No Filter, RSS filters, Distance filters) 
##           are applied prior to trilateration analysis
##  
##########################################################################################################################################################


# Packages needed
library(tidyverse)
library(raster)
library(nlstools)

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
test.info <- read.csv("Test.Info_Example.csv", header = T)
str(test.info) # check that data imported properly

beep.dat <- readRDS("BeepData_Example.rds") 
str(beep.dat) # check that data imported properly

nodes <- read.csv("Nodes_Example.csv", header = T)
str(nodes)



#################################################################
#   Step 1
# Run function to isolate raw RSS data from a node network that
# is associated with test data  
#################################################################

    ## Variables to define for function
        ## TEST.TYPE = category indicating the type of dataset (Calibration or LocError)
                ## "Calibration" - test dataset used for calibrating the relationship between RSS and Distance (purpose of this R script)
                ## "LocError" - test data used to determine localization error associated with RSS-based localization estimates for a node network
        ## DATE.FORMAT = format of the date column in Test.Info.csv file using R standard date expressions (e.g., "%m-%d-%y", "%Y-%m-%d")
        ## TIMEZONE = Time zone where data was collected, use grep("<insert location name here>", OlsonNames(), value=TRUE) to find valid time zone name


    ## Output in R environment
        ## combined.data - dataframe that contains the average RSS values for a given node associated with each unique test
            ## Columns:
                ## NodeId - unique identifier of the specified node
                ## TestId - unique identifier of the specified test
                ## avgRSS - average RSS value across a 3-min time period for the given node and test
                ## sdRSS - standard deviation of RSS values across the 3-min time period for the given node and test
                ## distance - true distance (in meters) between the specified node and test location
                ## NodeUTMx - Easting location of the specified node
                ## NodeUTMy - Northing location of the specified node
                ## TestUTMx - Easting location of the specified test
                ## TestUTMy - Northing location of the specified test

    
##******* Define Variables - replace values below with user specified values *******## 
TEST.TYPE <- "LocError"
DATE.FORMAT <- "%m/%d/%y"
TIME.ZONE <- "Pacific/Guam"


# Function to combine RSS data from a node network with test information into a dataframe
combined.data <- data.setup(TEST.TYPE, DATE.FORMAT, TIME.ZONE)




##########################################################################################################
# Step 2 - Estimate Distance
#     Function to estimate the distance of each test location based on the RSS value detected by a node
#     Use the model coefficients from exponential decay model created with calibration data
###########################################################################################################


    ## Variables to define for function  
        ## Coefficients in the exponential model: use output from Github_RSS_by_Distance_Calibration (e.g., coef(nls.mod))
            ## a = intercept
            ## S = decay factor
            ## K = horizontal asymptote

    ## Output in R environment
        ## combined.data - with the added column:
            ## e.dist = estimated distance of the test signal from the specified node based on RSS values

    ## Output saved
        ## LocError_Dataset.csv - .csv file of the dataframe 'combine.data' saved in the folder specified by the outpath


##******** Define Variables - replace values below with user specified values ********##  
a <- 47.23 
S <- 0.005
K <- -105.16


# Function to estimate the distance of each test signal from the node based on the RSS values detected
combined.data <- estimate.distance(combined.data)




####################################
# Step 3 - Trilateration Estimates
####################################

########################################################################################
# Trilateration of Test Data - No filter
#    Function that will estimate the location of each test location in the test
#    dataset with no filters applied prior to trilateration analysis and then 
#    summarize localization error
#######################################################################################

    ## Output in R environment
        ## no.filters - dataframe that contains the average localization error when no filters are applied prior to trilateration analysis
            ## Columns
                ## n.pts - number of random points a location estimate could be calculated based on trilateration 
                ## avg.no.nodes - average number of nodes used in each trilateration analysis
                ## avg.diff - average localization error of trilateration estimates when the specified filter was applied
                             ## localization error = difference in distance (m) between true and estimated location
                ## sd.diff - standard deviation of localization error of trilateration estimates when the specified filter was applied
                ## lower.ci - lower 95% confidence interval of localization error of trilateration estimates when the specified filter was applied
                ## upper.ci - upper 95% confidence interval of localization error of trilateration estimates when the specified filter was applied
                ## med.diff - median localization error of trilateration estimates when the specified filter was applied
                ## min.diff - minimum localization error of trilateration estimates when the specified filter was applied
                ## max.dist - maximum localization error of trilateration estimates when the specified filter was applied
                ## filter - filter applied prior to trilateration analysis


    ## Output saved
        ## Trilateration.TestData_NoFilter_Summary.Stats.csv - .csv file of the dataframe 'no.filters' saved in the folder specified by the outpath
        ## Trilateration.TestData_TestData_NoFilter_Summary.Results.csv - one .csv file that has the summary statistics of localization error when the specified
        ##                                                                filter was applied prior to trilateration. Each row represents the summary statistics 
        ##                                                                for a given test location. Files are for your reference and not used in downstream functions. 



# Calculate error of location estimates of each test location when no filters are applied prior to trilateration 
no.filters <- trilateration.TestData.NoFilter(combined.data)




#####################################################################################
# Trilateration of Test Data - RSS filters
#    Function that will estimate the location of each test location in the test
#    dataset when RSS filters are applied prior to trilateration and then 
#    summarize localization error
######################################################################################

    ## Variable to define for function 
        ## RSS.FILTER = RSS values that will be used to filter data prior to trilateration analysis
            ## All RSS values less than the specified value will be removed from the test data


    ## Output in R environment
        ## RSS.filters - dataframe that contains the overall localization error when RSS filters are applied prior to trilateration analysis
            ## columns are the same as described in the no.filter table
                  
    ## Output saved
        ## Trilateration.TestData_Filters.RSS_Summary.Stats.csv - .csv file of the dataframe 'RSS.filters' saved in the folder specified by the outpath
        ## _Summary.Results.csv - 4 .csv files each with the RSS filter name in the title that has the summary statistics of localization error when the specified
        ##                        filter was applied prior to trilateration. Each row represents the summary statistics for a given test location. 
        ##                        Files are for your reference and not used in down stream functions. 




##****** Define variables - indicate the RSS values to filter data prior to trilateration ************#
RSS.FILTER <- c(-80, -85, -90, -95)


# Calculate error of location estimates of each test location when RSS filters are applied prior to trilateration 
RSS.filters <- trilateration.TestData.RSS.Filter(combined.data)






###############################################################################################
# Trilateration - Distance from strongest node filter 
#   Function that will estimate the location of test location in the test dataset when
#   distance filters are applied prior to trilateration and then summarize localization error 
################################################################################################

    ## Variable to define for function
        ## DIST.FILTER = distances that are a multiple (x1.25, x2, x3, x4) of the average grid spacing
            # For each random test location the node with the strongest signal is identified and then only nodes
            # within the specified distance from the strongest node are kept for trilateration analysis


    ## Output in R environment
        ## Dist.filters - dataframe that contains the overall localization error when Distance filters are applied prior to trilateration analysis
            ## columns are the same as described in the no.filter table

    ## Output saved
        ## Trilateration.TestData_Filters.Distance_Summary.Stats.csv - .csv file of the dataframe 'Dist.filters' saved in the folder specified by the outpath
        ## _Summary.Results.csv - 4 .csv files each with the Distance filter name in the title that has the summary statistics of localization error when the specified
        ##                        filter was applied prior to trilateration. Each row represents the summary statistics for a given test location. 
        ##                        Files are for your reference and not used in down stream functions. 



##****** Define variables - indicate the Distance values to filter data prior to trilateration ************#
DIST.FILTER <- c(315,500,750,1000)


# Calculate error of location estimates of each test location when Distance filters are applied prior to trilateration 
Dist.filters <- trilateration.TestData.Distance.Filter(combined.data)



