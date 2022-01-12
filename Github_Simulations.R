###########################################################################################################################################################
##
##  Kristina L Paxton
##
##  April 1 2021
##
##    -- Code to simulate 3 node configurations within a 12.5 km2 grid - node spacing: 100m (n = 100 nodes), 175m (n = 64 nodes), or 250m (n = 36 nodes)
##    -- Create simulated dataset of RSS values for a given node configuration incorporating variability in RSS values based on noise in your node network
##    -- Assess the error associated with RSS-based localization estimates when filters are applied prior to trilateration analysis 
##        -- localization error = difference in distance (m) between true and estimated location
##      
##
##    Files Needed
##        1. Calibration_Dataset.csv that is saved in the working directory defined below
##            - Test dataset used to calibrate RSS by Distance relationship
##            - File was created in Github_RSS_by_Distance_Calibration.R script 
##            - File MUST have the two following olumns:
##                -- avgRSS: average RSS value for a 1-min time period for the specified node during the specified test
##                -- distance: true distance between the specified node and specified test
##
##        2. Functions_RSS.Based.Localizations.R - R script with functions needed to run the code below - file should be saved in the working directory defined below    
##      
##    Important Output
##        1. Dataset of simulated RSS values for a given node configuration
##        2. Dataframes with summary statistics of localization error when different filters (No Filter, RSS filters, Distance filters) 
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


## Bring Needed File - RSS by Distance Calibration Data - - change file name in " " as needed
combined.data <- read.csv("Calibration_Dataset.csv", header = T)





###########################################################################
# Create Grid of Nodes with Specified Spacing & Generate 100 Random Points
##########################################################################

    ## Variables to define for function
        ## SIZE = spacing of the nodes - options: "n100", "n175", "n250"

    ## Output in R environment
        ## nodes = dataframe containing x and y locations of nodes with the designated spacing within a 12.5 km2 grid
        ## random.points = dataframe with x and y locations of 100 random points within the node configuration
        ## DIST.FILTER = vector with the distances that are a multiple (x1.25, x2, x3, x4) of the average grid spacing indicated
        ## plot of nodes (black circles) and random points (grey triangles) for the designated spacing


##******** Define Variable - replace value below with "n100", "n175", or "n250" ********##
SIZE <- "n250"


# Create node configuration and random points
create.node.config(SIZE)





##################################################################################################################################
# Function to simulate RSS values for a given node configuration incorporating variability based on noise in your node network
# Step 1: Use the true distance between a random point and a node to randomly select an RSS value
#         associated with that distance in the RSS by Distance Calibration Dataset
# Step 2: Estimate the distance between the random point and node based on the randomly selected RSS value 
#         and the relationship between RSS values and distance determined based on an exponential decay function 
################################################################################################################################

    ## Variables to define for function  
        ## NUM.SIM = number of simulations to run (recommend between 100 to 1,000)
            # Simulations in the paper were based on 1,000 - the more simulations the longer time to generate data, but the more robust results
        ## Coefficients in the exponential model: use output from Github_RSS_by_Distance_Calibration (e.g., coef(nls.mod))
            ## a = intercept
            ## S = decay factor
            ## K = horizontal asymptote
  

    ## Output in Console
        ## Iteration of the simulation will be printed at the end of each cycle to help you keep track of how fast the simulation is running

    ## Output in R environment
        ## random.rssi_results.all.sim - dataframe that contains the simulated values that will be needed for trilateration analysis
            ## Columns
                  ## NodeId - unique identifier of the specified node
                  ## PointId - random point id of the specified test location
                  ## t.dist - true distance between the specified node and test location
                  ## r.rss - random RSS value selected from the dataset for the given true distance
                  ## n.rss - number of potential RSS values for the given true distance
                  ## sd.rss - standard deviation of the RSS values for the given true distance
                  ## e.dist - estimated distance between the specified node and test location based on the random RSS value selected (e.g., relationship between RSS values and distance) 
                  ## SimId - unique identifier given to each simulation
                  ## x - location on the x-axis of the random point in the 12.5 km2 grid
                  ## y - location on the y-axis of the random point in the 12.5 km2 grid
                  ## TestId - unique identifier indicating the random point location and the simulation

    ## Output saved
        ## Simulated.rss.values_ - .rds file of 'random.rssi_results.all.sim' file in specified outpath with the spacing of the nodes in the name
        ##                          indicated at the end of the file name (e.g., Simulated.rss.values_250m.rds)


##******** Define Variables - replace values below with user specified values ********##  
NUM.SIM <- 100
a <- 47.23 
S <- 0.005
K <- -105.16


# Create simulated dataset to be used for trilateration analysis
random.rss_results.all.sim <- get.RSS.values(NUM.SIM, a, S, K)






########################################################################################
# Trilateration of Simulated Data - No filter
    # Function that will estimate the location of each random point in the
    # simulated dataset with no filters applied prior to trilateration analysis
#######################################################################################

    ## Output in Console
        ## Iteration of the simulation will be printed at the end of each cycle to help you keep track of how fast the simulation is running
            ## Total number of potential location estimates = 100 random points x NUM.SIM 

    ## Output in R environment
        ## no.filters - dataframe that contains the average localization error when no filters are applied prior to trilateration 
            ## Columns
                ## n.pts - number of random points a location estimate could be calculated based on trilateration 
                ## avg.no.sim - average number of simulations for each random point when a location estimate could be calculated based on trilateration 
                ## avg.no.nodes - average number of nodes used in each trilateration analysis
                ## avg.diff - average localization error of trilateration estimates when the specified filter was applied
                ##            localization error = difference in distance (m) between true and estimated location
                ## sd.diff - standard deviation of localization error of trilateration estimates when the specified filter was applied
                ## lower.ci - lower 95% confidence interval of localization error of trilateration estimates when the specified filter was applied
                ## upper.ci - upper 95% confidence interval of localization error of trilateration estimates when the specified filter was applied
                ## med.diff - median localization error of trilateration estimates when the specified filter was applied
                ## min.diff - minimum localization error of trilateration estimates when the specified filter was applied
                ## max.dist - maximum localization error of trilateration estimates when the specified filter was applied
                ## filter - filter applied prior to trilateration analysis

    ## Output saved
        ## Trilateration.Simulation_NoFilters_Summary.Stats.csv - .csv file of the dataframe 'No.filters' saved in the folder specified by the outpath
        ## Trilateration.Simulation_NoFilter_Summary.Results.csv - .csv file that has the summary statistics of localization error of when no filter  
        ##                                                          was applied prior to trilateration. Each row represents the summary statistics 
        ##                                                          for all simulations for a given random point. Files are for your reference 




# Estimate locations of 100 random points when no filters are applied prior to trilateration 
no.filters <- trilateration.Sim.NoFilter(random.rss_results.all.sim)


 



#####################################################################################
# Trilateration - RSS filters
    # Function that will estimate the location of each random point
    # in the simulated dataset when RSS filters are applied prior to trilateration 
######################################################################################

    ## Variable to define for function 
        ## RSS.FILTER = RSS values that will be used to filter data prior to trilateration analysis
                        ## All RSS values below the specified value will be removed from the simulated data


    ## Output in Console
        ## Iteration will be printed indicating which RSS filter is currently being processed
      

    ## Output in R environment
        ## RSS.filters - dataframe that contains the overall localization error when RSS filters are applied prior to trilateration analysis
            ## columns are the same as described in the no.filter table


   ## Output saved
        ## Trilateration.Simulation_Filters.RSS_Summary.Stats.csv - .csv file of the dataframe 'RSS.filters' saved in the folder specified by the outpath 
        ## _Filters.RSS_Summary.Results.csv - 4 .csv file each with the RSS filter name in the title that has the summary statistics of localization error 
        ##                                    when the specified filter was applied prior to trilateration. Each row represents the summary statistics
        ##                                    for all simulations for a given random point. Files are for your reference. 




##****** Define variables - indicate the RSS values to filter data prior to trilateration ************#
RSS.FILTER <- c(-80, -85, -90, -95)


# Estimate locations of 100 random points when RSS filters are applied prior to trilateration 
RSS.filters <- trilateration.Sim.RSS.Filter(random.rss_results.all.sim)






###########################################################################################
# Trilateration - Distance from strongest node filter 
    # Function that will estimate the location of each random point
    # in the simulated dataset when distance filters are applied prior to trilateration 
###########################################################################################


    ## Output in Console
        ## Iteration will be printed indicating which distance filter is currently being processed


    ## Output in R environment
        ## Dist.filters - dataframe that contains the overall localization error when distance filters are applied prior to trilateration analysis
            ## columns are the same as described in the no.filter table


   ## Output saved
        ## Trilateration.Simulation_Filters.Distance_Summary.Stats.csv - .csv file of the dataframe 'Dist.filters' saved in the folder specified by the outpath
        ## Filters.Distance_Summary.Results.csv - 4 .csv file each with the Distance filter name in the title that has the summary statistics of  
        ##                                        localization error when the specified filter was applied prior to trilateration. Each row represents  
        ##                                        the summary statistics for all simulations for a given random point. Files are for your reference. 




# Estimate locations of 100 random points when distance filters are applied prior to trilateration 
    # Distance filter was created in previous step when you defined the node configuration
        # DIST.FILTER = distances that are a multiple (x1.25, x2, x3, x4) of the average grid spacing
            # For each random test point the node with the strongest signal is identified and then only nodes
            # within the specified distance from the strongest node are kept for trilateration analysis
Dist.filters <- trilateration.Sim.Distance.Filter(random.rss_results.all.sim)





