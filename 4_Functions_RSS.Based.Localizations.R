###################################################################################
## Kristina L Paxton
##
## April 26 2021
##
## Functions used to calculate data generated in RSS-based Localization Paper
##
##
###################################################################################



################################################################################ 
# Run function to bring together test information and raw RSS values 
# and create dataset used in analysis - either from calibration or test data 
################################################################################


data.setup <- function(TEST.TYPE, DATE.FORMAT, TIME.ZONE) {

## Prepare test info dataset

# Format date
test.info$Date <- as.Date(test.info$Date, format = DATE.FORMAT)

# Bring together tag, date, hour, min
test.info <- test.info %>%
  dplyr::mutate(Tag.Date.HM = paste(TagId, Date, Hour, Min, sep = "_"))

# Make TestId a factor
test.info$TestId <- as.factor(test.info$TestId)

# Make a dataframe with only 1 row per test
test.UTM <- test.info %>%
  dplyr::group_by(TestId) %>%
  dplyr::slice_head(n=1)


## Prepare the beep dataset

# Add column for date, hour and min and then bring together with TagId
beep.dat$Hour <- lubridate::hour(beep.dat$Time.local)
beep.dat$Min <- lubridate::minute(beep.dat$Time.local)
beep.dat$Date <- as.Date(beep.dat$Time.local, tz = TIME.ZONE)
beep.dat <- beep.dat %>%
  dplyr::mutate(Tag.Date.HM = paste(TagId, Date, Hour, Min, sep = "_"))


## Isolate raw RSS values associated with the test data
test.dat <- beep.dat %>%
  dplyr::filter(Tag.Date.HM %in% test.info$Tag.Date.HM) %>%
  dplyr::left_join(test.info[,c("Tag.Date.HM", "TestId")])
test.dat$NodeId <- droplevels(test.dat$NodeId)
test.dat$TagId <- droplevels(test.dat$TagId)


## Calculate average RSS value for each unique test and node 
summary.test.tags <- test.dat %>%
  dplyr::group_by(NodeId, TestId) %>%
  dplyr::summarise(avgRSS = mean(TagRSSI),
                   sdRSS = sd(TagRSSI),
                   n.det = n())


## Calculate euclidean distance between each test location and node

# Calculate distance between nodes and test locations
dst <- raster::pointDistance(test.UTM[,c("TestUTMx", "TestUTMy")], nodes[,c("NodeUTMx", "NodeUTMy")], lonlat = F, allpairs = T)

# Make matrix into a dataframe
dist_df <- data.frame(dst, row.names = test.UTM$TestId)
colnames(dist_df) <- nodes$NodeId
dist_df$TestId <- rownames(dist_df)

# rearrange data
dist.gather <- dist_df %>%
  tidyr::gather(key = "NodeId", value = "distance", -TestId)


## Combine distances with summary data 
summary.dist <- summary.test.tags %>%
  dplyr::left_join(dist.gather) 

# Add UTMs of nodes and test locations to dataset
summary.dist <- summary.dist %>%
  dplyr::left_join(nodes[, c("NodeId", "NodeUTMx", "NodeUTMy")]) %>%
  dplyr::left_join(test.UTM[, c("TestId", "TestUTMx", "TestUTMy")]) 


## save file
write.csv(summary.dist, paste0(outpath, TEST.TYPE, "_Dataset.csv"),
          row.names = F)


## Add to R environment
return(summary.dist)


}





#########################################################
# Function to estimate distance of a signal based
# on RSS value and the exponential decay between
# RSS values and distance
########################################################

estimate.distance <- function(x) {

# supress warnings
options(warn = -1)

# Calculate estimated distance based on RSSI~Distance relationship and indicate simulation round
combined.data <- x %>%
  dplyr::mutate(e.dist = (log(avgRSS - K) - log(a)) / -S)
          
# Remove rows with NAs
combined.data <- combined.data[complete.cases(combined.data),]

# Change negative distances to 10 (rss values > intercept of exponential curve and thus negative) - indicates very close to the node
combined.data <- combined.data %>%
  dplyr::mutate(e.dist = dplyr::case_when(e.dist < 0 ~ 10,
                                          e.dist >=0 ~ e.dist))

# Save File to outpath
write.csv(combined.data, paste0(outpath, "LocError_Dataset.csv"),
          row.names = F)

return(combined.data)

}




###############################################
# Trilateration of Test Data - No filter
################################################


trilateration.TestData.NoFilter <- function(x) {
  
  # supress warnings
  options(warn = -1)
  
  # make a vector of unique trilaterations to run
  tests = unique(x$TestId)
  
  # Make a dataframe with only 1 row per test
  test.UTM <- x %>%
    dplyr::group_by(TestId) %>%
    dplyr::slice_head(n=1) %>%
    dplyr::select(TestId, TestUTMx, TestUTMy)
  
  # Create a dataframe for output estimates
  estimated.location_results <- data.frame(TestId=character(), No.Nodes = numeric(), x.est=numeric(), y.est=numeric(), 
                                           x.ci.lower =numeric(), x.ci.upper =numeric(), y.ci.lower = numeric(), y.ci.upper = numeric())
  
  
  for(j in 1:length(tests)) {
    
    # Isolate the test 
    sub.test <- x %>% dplyr::filter(TestId == tests[j]) 
    
    # Determine the node with the strongest RSS value
    max.RSS <- sub.test[which.max(sub.test$avgRSS),]
    
    # Calculate no nodes for the test
    no.nodes <- dplyr::n_distinct(sub.test$NodeId)
    
    
    # To deal with potential errors where the model fails due to bad starting values using tryCatch everything you want evaluated by tryCatch goes inside {},
    # then the error will be printed but the loop will continue
    
    # Non-linear test to optimize the location of unknown signal by looking at the radius around each Node based on estimated distance and the pairwise distance between all nodes
    tryCatch( {nls.test <- nls(e.dist ~ raster::pointDistance(data.frame(NodeUTMx, NodeUTMy), c(NodeUTMx_solution, NodeUTMy_solution), lonlat = F, allpairs = T),
                               data = sub.test, start=list(NodeUTMx_solution=max.RSS$NodeUTMx, NodeUTMy_solution=max.RSS$NodeUTMy),
                               control=nls.control(warnOnly = T, minFactor=1/30000, maxiter = 100)) # gives a warning, but doesn't stop the test from providing an estimate based on the last itteration before the warning
    
    
    
    # Determine an error around the point location estimate
    par.est = cbind(coef(nls.test), confint2(nls.test))
    lng.ci.upper =  par.est[1,3] 
    lng.ci.lower =  par.est[1,2]
    lat.ci.upper =  par.est[2,3] 
    lat.ci.lower =  par.est[2,2]}
    
    ,error = function(e)  {cat("ERROR :",conditionMessage(e), "\n")})
    
    # estimated location of the point and error
    estimated.loc <- data.frame(TestId = tests[j], No.Nodes = no.nodes, x.est = par.est[1,1], y.est = par.est[2,1], 
                                x.ci.lower = lng.ci.lower, x.ci.upper = lng.ci.upper,  y.ci.lower = lat.ci.lower, y.ci.upper = lat.ci.upper)
    
    # Populate dataframe with results
    estimated.location_results <- rbind(estimated.location_results, estimated.loc)
    
  }
  
  # combine estimated locations with true locations
  combined_results <- estimated.location_results %>%
    dplyr::left_join(test.UTM)
  
  # Calculate difference distance between estimated and true location    
  dst <- raster::pointDistance(combined_results[,9:10], combined_results[,c(3:4)], lonlat = F, allpairs = F)
  
  # bring all together
  combined_results_final <- dplyr::bind_cols(combined_results, data.frame(diff.dist = dst)) %>%
    dplyr::mutate(filter = "No Filter")
  
  # save file
  write.csv(combined_results_final, paste0(outpath, "Trilateration.TestData_NoFilter.RSS_Results.csv"),  row.names = F)
  
  # summarize statitics for a given filter
  summary.stats <- combined_results_final %>%
    dplyr::summarise(n.est.tests = dplyr::n_distinct(TestId),
                     avg.no.nodes = mean(No.Nodes),
                     avg.diff = mean(diff.dist),
                     sd.diff = sd(diff.dist),
                     lower.ci = avg.diff - qnorm(0.975)*sd.diff/sqrt(n.est.tests),
                     upper.ci = avg.diff + qnorm(0.975)*sd.diff/sqrt(n.est.tests),
                     med.diff = median(diff.dist),
                     min.diff = min(diff.dist),
                     max.dist = max(diff.dist)) %>%
    dplyr::mutate(filter = "No Filter")
  
  # save file
  write.csv(summary.stats, paste0(outpath, "Trilateration.TestData_NoFilter_Summary.Stats.csv"),  row.names = F)
  
  return(summary.stats)
}






###############################################
# Trilateration of Test Data - RSS filters
################################################


trilateration.TestData.RSS.Filter <- function(x) {
  
  # supress warnings
  options(warn = -1)
  
  # Empty data frame to populate with summary results
  summary.stats_results <- data.frame(n.est.tests = numeric(), avg.no.nodes = numeric(), avg.diff = numeric(), sd.diff = numeric(),
                                      lower.ci = numeric(), upper.ci = numeric(), med.diff = numeric(), 
                                      min.diff = numeric(), max.diff = numeric(), filter = character())
  
  # Make a dataframe with only 1 row per test
  test.UTM <- x %>%
    dplyr::group_by(TestId) %>%
    dplyr::slice_head(n=1) %>%
    dplyr::select(TestId, TestUTMx, TestUTMy)
  
  for (i in 1:length(RSS.FILTER)) {
    
    # Remove RSSI <= filter
    combined.data <- x %>%
      dplyr::filter(avgRSS >= RSS.FILTER[i])
    
    # Filter list so only contains tests which had <3 nodes detecting the transmitter
    sample.size <- combined.data %>%
      dplyr::group_by(TestId) %>%
      dplyr::summarise(n.nodes = n()) %>%
      dplyr::filter(n.nodes < 3)
    
    # Remove Tests with <3 nodes picking up a signal > value indicated
    combined.data.red <- combined.data %>%
      dplyr::filter(!(TestId %in% sample.size$TestId)) 
  
    # make a vector of unique trilaterations to run
    tests = unique(combined.data.red$TestId)
    
   # Create a dataframe for output estimates
    estimated.location_results <- data.frame(TestId=character(), No.Nodes = numeric(), x.est=numeric(), y.est=numeric(), 
                                             x.ci.lower =numeric(), x.ci.upper =numeric(), y.ci.lower = numeric(), y.ci.upper = numeric())
    
    
    for(j in 1:length(tests)) {
      
      # Isolate the test 
      sub.test <- combined.data.red %>% dplyr::filter(TestId == tests[j]) 
      
      # Determine the node with the strongest RSS value
      max.RSS <- sub.test[which.max(sub.test$avgRSS),]
      
      # Calculate no nodes for the test
      no.nodes <- dplyr::n_distinct(sub.test$NodeId)
      
      
      # To deal with potential errors where the model fails due to bad starting values using tryCatch everything you want evaluated by tryCatch goes inside {},
      # then the error will be printed but the loop will continue
      
      # Non-linear test to optimize the location of unknown signal by looking at the radius around each Node based on estimated distance and the pairwise distance between all nodes
      tryCatch( {nls.test <- nls(e.dist ~ raster::pointDistance(data.frame(NodeUTMx, NodeUTMy), c(NodeUTMx_solution, NodeUTMy_solution), lonlat = F, allpairs = T),
                                 data = sub.test, start=list(NodeUTMx_solution=max.RSS$NodeUTMx, NodeUTMy_solution=max.RSS$NodeUTMy),
                                 control=nls.control(warnOnly = T, minFactor=1/30000, maxiter = 100)) # gives a warning, but doesn't stop the test from providing an estimate based on the last itteration before the warning
      
      
      
      # Determine an error around the point location estimate
      par.est = cbind(coef(nls.test), confint2(nls.test))
      lng.ci.upper =  par.est[1,3] 
      lng.ci.lower =  par.est[1,2]
      lat.ci.upper =  par.est[2,3] 
      lat.ci.lower =  par.est[2,2]}
      
      ,error = function(e)  {cat("ERROR :",conditionMessage(e), "\n")})
      
      # estimated location of the point and error
      estimated.loc <- data.frame(TestId = tests[j], No.Nodes = no.nodes, x.est = par.est[1,1], y.est = par.est[2,1], 
                                  x.ci.lower = lng.ci.lower, x.ci.upper = lng.ci.upper,  y.ci.lower = lat.ci.lower, y.ci.upper = lat.ci.upper)
      
      # Populate dataframe with results
      estimated.location_results <- rbind(estimated.location_results, estimated.loc)
      
    }
    
    # combine estimated locations with true locations
    combined_results <- estimated.location_results %>%
      dplyr::left_join(test.UTM)
    
    # Calculate difference distance between estimated and true location    
    dst <- raster::pointDistance(combined_results[,9:10], combined_results[,c(3:4)], lonlat = F, allpairs = F)
    
    # bring all together
    combined_results_final <- dplyr::bind_cols(combined_results, data.frame(diff.dist = dst)) %>%
      dplyr::mutate(filter = paste0("RSS", " ", RSS.FILTER[i]))
    
    # save file
    write.csv(combined_results_final, paste0(outpath, "Trilateration.Test.Data_Filter.RSS.", RSS.FILTER[i], "_Results.csv"),  row.names = F)
    
    # summarize statitics for a given filter
    summary.stats <- combined_results_final %>%
      dplyr::summarise(n.est.tests = dplyr::n_distinct(TestId),
                       avg.no.nodes = mean(No.Nodes),
                       avg.diff = mean(diff.dist),
                       sd.diff = sd(diff.dist),
                       lower.ci = avg.diff - qnorm(0.975)*sd.diff/sqrt(n.est.tests),
                       upper.ci = avg.diff + qnorm(0.975)*sd.diff/sqrt(n.est.tests),
                       med.diff = median(diff.dist),
                       min.diff = min(diff.dist),
                       max.dist = max(diff.dist)) %>%
      dplyr::mutate(filter = paste0("RSS", " ", RSS.FILTER[i]))
    
    # write summary to empty dataframe
    summary.stats_results <- rbind(summary.stats_results, summary.stats)
    
    
  }
  
  # save file
  write.csv(summary.stats_results, paste0(outpath, "Trilateration.Test.Data_Filters.RSSI_Summary.Stats.csv"),  row.names = F)
  
  return(summary.stats_results)

}




###############################################
# Trilateration of Test Data - Distance Filters
################################################

trilateration.TestData.Distance.Filter <- function(x){
  
  # supress warnings
  options(warn = -1)
  
  # Make a dataframe with only 1 row per test
  test.UTM <- x %>%
    dplyr::group_by(TestId) %>%
    dplyr::slice_head(n=1) %>%
    dplyr::select(TestId, TestUTMx, TestUTMy)

  # make a vector of unique trilaterations to run
  tests = unique(x$TestId)
  
  # Calculate distance between nodes
  dist.nodes <- raster::pointDistance(nodes[,c("NodeUTMx", "NodeUTMy")], nodes[,c("NodeUTMx", "NodeUTMy")], lonlat = F, allpairs = T)
  
  # Make matrix into a dataframe with a row for NodeId
  dist.nodes_df <- data.frame(dist.nodes, row.names = nodes$NodeId)
  colnames(dist.nodes_df) <- nodes$NodeId
  dist.nodes_df$NodeId <- rownames(dist.nodes_df)

  
  # Empty data frame to populate with summary results
  summary.stats_results <- data.frame(n.est.tests = numeric(), avg.no.nodes = numeric(), avg.diff = numeric(), sd.diff = numeric(),
                                      lower.ci = numeric(), upper.ci = numeric(), med.diff = numeric(), 
                                      min.diff = numeric(), max.diff = numeric(), filter = character())
  

# loop through distances

for(i in 1:length(DIST.FILTER)) {
  
  # Create a dataframe for output estimates
  estimated.location_results <- data.frame(TestId=character(), No.Nodes = numeric(), x.est=numeric(), y.est=numeric(), 
                                           x.ci.lower =numeric(), x.ci.upper =numeric(), y.ci.lower = numeric(), y.ci.upper = numeric())
  
  # identify the distance
  dist.no <- DIST.FILTER[i]
  

  # Loop through random points
  for(j in 1:length(tests)) {
    

    # Isolate the test 
    sub.test <- x %>% dplyr::filter(TestId == tests[j]) 
    
    # Determine the node with the strongest RSSI value
    max.RSS <- sub.test[which.max(sub.test$avgRSS),]
    
    # Filter matrix of node distances by node with strongest avg.RSSI
    # and get all nodes within a particular distance of the node
    nodes.test <- dist.nodes_df %>%
      dplyr::filter(NodeId %in% max.RSS$NodeId) %>%
      tidyr::gather(key = "NodeId", value = "distance", -NodeId) %>%
      dplyr::filter(distance <= dist.no)
    
    # Only keep nodes that are within the specified distance of node with strongest avg.RSSI
    sub.test.red <- sub.test %>%
      dplyr::filter(NodeId %in% nodes.test$NodeId)
    
    # Calculate no nodes for the test
    no.nodes <- dplyr::n_distinct(sub.test.red$NodeId)
    
    # If the number of nodes is not greater than 3 the rest of the loop will not be continued and the next
    # iteration of the loop is started
    if(no.nodes < 3) {
      next
    }
    
    # To deal with potential errors where the model fails due to bad starting values using tryCatch everything you want evaluated by tryCatch goes inside {},
    # then the error will be printed but the loop will continue
      
    # Non-linear test to optimize the location of unknown signal by looking at the radius around each Node based on estimated distance and the pairwise distance between all nodes
    tryCatch( {nls.test <- nls(e.dist ~ raster::pointDistance(data.frame(NodeUTMx, NodeUTMy), c(NodeUTMx_solution, NodeUTMy_solution), lonlat = F, allpairs = T),
                               data = sub.test.red, start=list(NodeUTMx_solution=max.RSS$NodeUTMx, NodeUTMy_solution=max.RSS$NodeUTMy),
                               control=nls.control(warnOnly = T, minFactor=1/30000, maxiter = 100)) # gives a warning, but doesn't stop the test from providing an estimate based on the last itteration before the warning
    
    
    # Determine an error around the point location estimate
    par.est = cbind(coef(nls.test), confint2(nls.test))
    lng.ci.upper =  par.est[1,3] 
    lng.ci.lower =  par.est[1,2]
    lat.ci.upper =  par.est[2,3] 
    lat.ci.lower =  par.est[2,2] }
    
    ,error = function(e)  {cat("ERROR :",conditionMessage(e), "\n")})
    
    # estimated location of the point and error
    estimated.loc <- data.frame(TestId = tests[j], No.Nodes = no.nodes, x.est = par.est[1,1], y.est = par.est[2,1], 
                                x.ci.lower = lng.ci.lower, x.ci.upper = lng.ci.upper, y.ci.lower = lat.ci.lower, y.ci.upper = lat.ci.upper)
    
    # Populate dataframe with results
    estimated.location_results <- rbind(estimated.location_results, estimated.loc)
    
  }
  
  
  # combine estimated locations with true locations
  combined_results <- estimated.location_results %>%
    dplyr::left_join(test.UTM)
  
  # Calculate difference distance between estimated and true location    
  dst <- raster::pointDistance(combined_results[,9:10], combined_results[,c(3,4)], lonlat = F, allpairs = F)
  
  # bring all together
  combined_results_final <- dplyr::bind_cols(combined_results, data.frame(diff.dist = dst)) %>%
    dplyr::mutate(filter = paste0("Distance", " ", DIST.FILTER[i]))
  
  # save file
  write.csv(combined_results_final, paste0(outpath, "Trilateration.Test.Data_Filter.Distance.", DIST.FILTER[i], "_Results.csv"),  row.names = F)
  
  
  # Summary Stats of data for a given filter
  summary.stats <- combined_results_final %>%
    dplyr::summarise(n.est.tests = dplyr::n_distinct(TestId),
                     avg.no.nodes = mean(No.Nodes),
                     avg.diff = mean(diff.dist),
                     sd.diff = sd(diff.dist),
                     lower.ci = avg.diff - qnorm(0.975)*sd.diff/sqrt(n.est.tests),
                     upper.ci = avg.diff + qnorm(0.975)*sd.diff/sqrt(n.est.tests),
                     med.diff = median(diff.dist),
                     min.diff = min(diff.dist),
                     max.dist = max(diff.dist)) %>%
    dplyr::mutate(filter = paste0("Distance", " ", DIST.FILTER[i]))
  
  # populate empty dataframe
  summary.stats_results <- rbind(summary.stats_results, summary.stats)
  
  
}

  # save file
  write.csv(summary.stats_results, paste0(outpath, "Trilateration.Simulation_Filters.Distance_Summary.Stats.csv"),  row.names = F)
  
  return(summary.stats_results)
  
}





######################################## 
# Create Simulation Node Configuration 
########################################

## Create a grid of nodes, random points, and plot both
create.node.config <- function(SIZE) { 
  
  if(SIZE == "n100") {
    
    # Create nodes for a 100m spacing for a Grid 12 x 12
    node.x <- rep(c(0,100,200,300,400,500,600,700,800,900,1000,1100,1200), each = 13)
    node.y <- rep(c(0,100,200,300,400,500,600,700,800,900,1000,1100,1200), 13)
    NodeId <- paste0(rep("N",169), 1:169)
    nodes <- data.frame(x= node.x, y = node.y, NodeId = NodeId)
    
    # Create a set of 100 random points for both x and y within the grid
    ## set seed for reproducible results of the random sample of points
    set.seed(40549)
    random.x <- sample(0:1250, 100, replace = TRUE)
    random.y <- sample(0:1250, 100, replace = TRUE)
    PointId <- paste0(rep("P",100), 1:100)
    random.points <- data.frame(x = random.x, y = random.y, PointId = PointId)
    
    # Distance filter for 100m nodes
    DIST.FILTER <- c(125, 200, 300, 400) 
        
   # Create a list of the node configuration and random points and then add it to R environment
      new.list <- list("nodes" = nodes, "random.points" = random.points, "DIST.FILTER" = DIST.FILTER)
      list2env(new.list, .GlobalEnv)
        
    # Look at configuration of nodes and random points
      plot(y~x, data =nodes, bg = "black", pch = 21,
           xlim = c(-50, 1250), ylim = c(-50, 1250),
           xlab = "", ylab = "")
      points(y~x, data = random.points, bg = "grey", pch = 24)
      legend(-70, -5,
             legend = c("nodes", "random points"),
             pch = c(19, 17),
             col = c("black", "grey"),
             bty = "n")
  }
 
   
  if(SIZE == "n175") {
    
    # Create nodes for a 175m spacing for a
    node.x <- rep(c(0,175,350,525,700,875,1050,1225), each = 8)
    node.y <- rep(c(0,175,350,525,700,875,1050,1225), 8)
    NodeId <- paste0(rep("N",64), 1:64)
    nodes <- data.frame(x= node.x, y = node.y, NodeId = NodeId)
    
    # Create a set of 100 random points for both x and y within the grid
    ## set seed for reproducible results of the random sample of points
    set.seed(40549)
    random.x <- sample(0:1250, 100, replace = TRUE)
    random.y <- sample(0:1250, 100, replace = TRUE)
    PointId <- paste0(rep("P",100), 1:100)
    random.points <- data.frame(x = random.x, y = random.y, PointId = PointId)
    
    # Distance filter for 175 m nodes
    DIST.FILTER <- c(220,350,525,700) 
    
    # Create a list of the node configuration and random points and then add it to R environment
    new.list <- list("nodes" = nodes, "random.points" = random.points, "DIST.FILTER" = DIST.FILTER)
    list2env(new.list, .GlobalEnv)
    
    # Look at configuration of nodes and random points
    plot(y~x, data =nodes, bg = "black", pch = 21,
         xlim = c(-50, 1250), ylim = c(-50, 1250),
         xlab = "", ylab = "")
    points(y~x, data = random.points, bg = "grey", pch = 24)
    legend(-70, -5,
           legend = c("nodes", "random points"),
           pch = c(19, 17),
           col = c("black", "grey"),
           bty = "n")
    
   }
  



  if(SIZE == "n250") {
  
  # Create nodes for 250 m spacing for a Grid 12.5 km x 12.5 km
  node.x <- rep(c(0,250,500,750,1000,1250), each = 6)
  node.y <- rep(c(0,250,500,750,1000,1250), 6)
  NodeId <- paste0(rep("N",36), 1:36)
  nodes <- data.frame(x = node.x, y = node.y, NodeId = NodeId)
  
  # Create a set of 100 random points for both x and y within the grid
  ## set seed for reproducible results of the random sample of points
  set.seed(40549)
  random.x <- sample(0:1250, 100, replace = TRUE)
  random.y <- sample(0:1250, 100, replace = TRUE)
  PointId <- paste0(rep("P",100), 1:100)
  random.points <- data.frame(x = random.x, y = random.y, PointId = PointId)
  
  # Distance filter for 250 m nodes
  DIST.FILTER <- c(315,500,750,1000)
  
  
  # Create a list of the node configuration and random points and then add it to R environment
  new.list <- list("nodes" = nodes, "random.points" = random.points, "DIST.FILTER" = DIST.FILTER)
  list2env(new.list, .GlobalEnv)
  
  # Look at configuration of nodes and random points
  plot(y~x, data =nodes, bg = "black", pch = 21,
       xlim = c(-50, 1250), ylim = c(-50, 1250),
       xlab = "", ylab = "")
  points(y~x, data = random.points, bg = "grey", pch = 24)
  legend(-70, -5,
         legend = c("nodes", "random points"),
         pch = c(19, 17),
         col = c("black", "grey"),
         bty = "n")
  
  }

}





#############################################################################################
# Randomly Select RSSI value for a given distance between pairs of nodes and random points  
#    based on test data of RSSI~Distance
# Estimate distance between pairs of nodes and random points for a randomly selected RSSI
#    based on relationship between RSSI~Distance calculated from test data
#############################################################################################

get.RSS.values <- function(NUM.SIM, a, S, K){

## Isolate RSS and distance values from dataset and round to tenth place
combined.data <- combined.data  %>% 
    dplyr::select(avgRSS, distance)
combined.data$avgRSS <- round(combined.data$avgRSS, digits = 0)
combined.data$distance <- round(combined.data$distance, digits = 0)


## Calculate true distance between each random point and each node

  # Calculate euclidean distance between nodes and random points
dst <- raster::pointDistance(nodes[,1:2], random.points[,1:2], lonlat = F, allpairs = F)

  # Make matrix into a dataframe
dist_df <- data.frame(dst, row.names = nodes$NodeId)
colnames(dist_df) <- random.points$PointId 
dist_df$NodeId <- rownames(dist_df)

  # rearrange data so one column of distances for each NodeId and PointId pair
dist.gather <- dist_df %>%
  tidyr::gather(key = "PointId", value = "distance", -NodeId)

  # Round distance to nearest whole number
dist.gather$distance <- round(dist.gather$distance, digits = 0)


## Create empty dataframe for all simulations data
random.rss_results.all.sim <- data.frame(NodeId = character(), PointId = character(), SimId = character(), 
                                          t.dist = numeric(), e.dist = numeric(), r.RSS = numeric(), n.RSS = numeric(), sd.RSS = numeric())

# Repeat the following steps 1,000 times to simulate 1,000 draws of RSS for a given distance

for (j in c(1:NUM.SIM)) {
  
  
  # Create an empty dataframe to populate
  random.rss_results <- data.frame(NodeId = character(), PointId = character(), t.dist = numeric(), r.RSS = numeric(), 
                                    n.RSS = numeric(), sd.RSS = numeric())
  
  # Based on distance between nodes and random points get a Random RSS value based on RSS~ Distance Relationship
  for (i in 1:nrow(dist.gather)) {
    
    # isolate row of data in the simulated distance between nodes and random points
    sub <- dist.gather[i,]
    bandwidth = 0
    
    # In test data select all RSS values for a given distance in simulated dataset 
    pool <- combined.data$avgRSS[dplyr::between(combined.data$distance, sub$distance,
                                   sub$distance)]
    bw <- bandwidth
    # If there is less than 4 RSS value for a given distance try +-1 of that distance
    while (length(pool) < 4) {
      bw <- bw + 1
      pool <- combined.data$avgRSS[between(combined.data$distance, sub$distance - bw,
                              sub$distance + bw)]
    }
    # If there is less than 4 RSS value for a given distance try +-3 of that distance
    while (length(pool) < 4) {
      bw <- bw + 3
      pool <- combined.data$avgRSS[between(combined.data$distance, sub$distance - bw,
                              sub$distance + bw)]
    }
    # If there is less than 4 RSS value for a given distance try +-5 of that distance
    while (length(pool) < 4) {
      bw <- bw + 5
      pool <- combined.data$avgRSS[between(combined.data$distance, sub$distance - bw,
                              sub$distance + bw)]
    }
    # If there is less than 4 RSS value for a given distance try +-5 of that distance
    while (length(pool) < 4) {
      bw <- bw + 5
      pool <- combined.data$avgRSS[between(combined.data$distance, sub$distance - bw,
                              sub$distance + bw)]
    }
    # If there is less than 4 RSS value for a given distance try +-10 of that distance
    while (length(pool) < 4) {
      bw <- bw + 10
      pool <- combined.data$avgRSS[between(combined.data$distance, sub$distance - bw,
                              sub$distance + bw)]
    }
    
    
    # Take a random sample of the RSS values identified for a given distance 
    # and indicate the no. of RSS values and sd of RSS values for the given distance
    rRSS <- sample(pool, 1)
    n <- length(pool)
    sd <- sd(pool)
    
    # Make a dataframe of results 
    random.rss <- data.frame(NodeId = sub$NodeId, PointId = sub$PointId, t.dist = sub$distance, 
                              r.RSS = rRSS, n.RSS = as.numeric(n), sd.RSS = sd)
    
    # Populate dataframe with results
    random.rss_results <- rbind(random.rss_results, random.rss)
    
    
  }
  
  # supress warnings
  options(warn = -1)
  
  # Calculate estimated distance based on RSSI~Distance relationship and indicate simulation round
  random.rss_results <- random.rss_results %>%
    dplyr::mutate(e.dist = (log(r.RSS - K) - log(a)) / -S,
                  SimId = paste0("S", j))
  
  # Remove rows with NAs
  random.rss_results <- random.rss_results[complete.cases(random.rss_results),]
  
  # Populate dataframe with results
  random.rss_results.all.sim <- rbind(random.rss_results.all.sim, random.rss_results)
  
  
  print(j)
  
    }


  # Change negative distances to 10 (rss values > intercept of exponential curve and thus negative) - indicates very close to the node
  random.rss_results.all.sim <- random.rss_results.all.sim %>%
    dplyr::mutate(e.dist = dplyr::case_when(e.dist < 0 ~ 10,
                                          e.dist >=0 ~ e.dist))

  # Add Node locations to dataframe and create a new column to indicate random point and simulation
  random.rss_results.all.sim <- random.rss_results.all.sim %>%
    dplyr::left_join(nodes) %>% 
    dplyr::mutate(TestId = paste0(PointId, "_", SimId))


    # save results
    saveRDS(random.rss_results.all.sim, paste0(outpath, "Simulated.rss.values_", SIZE, ".rds"))

    return(random.rss_results.all.sim)
    
}






  
  


###############################################
# Trilateration of Simulated Data - No filter
################################################


trilateration.Sim.NoFilter <- function(x) {

# make a vector unique trilaterations to run
tests = unique(x$TestId)

# Create a dataframe for output estimates
estimated.location_results <- data.frame(TestId=character(), PointId = character(), No.Nodes = numeric(), x.est=numeric(), y.est=numeric(), 
                                         x.ci.lower =numeric(), x.ci.upper =numeric(), y.ci.lower = numeric(), y.ci.upper = numeric())


for(j in 1:length(tests)) {
  
  
  # Isolate the test 
  sub.test <- x %>% dplyr::filter(TestId == tests[j]) 
  
  # Determine the node with the strongest RSSI value
  max.RSS <- sub.test[which.max(sub.test$r.RSS),]
  
  # Identify the random point
  rd.pt <- unique(sub.test$PointId)
  
  # Calculate no nodes for the test
  no.nodes <- dplyr::n_distinct(sub.test$NodeId)
  
  
  # To deal with potential errors where the model fails due to bad starting values using tryCatch everything you want evaluated by tryCatch goes inside {},
  # then the error will be printed but the loop will continue
  
  # Non-linear test to optimize the location of unknown signal by looking at the radius around each Node based on estimated distance and the pairwise distance between all nodes
  tryCatch( {nls.test <- nls(e.dist ~ raster::pointDistance(data.frame(x, y), c(x_solution, y_solution), lonlat = F, allpairs = T),
                             data = sub.test, start=list(x_solution=max.RSS$x, y_solution=max.RSS$y),
                             control=nls.control(warnOnly = T, minFactor=1/30000, maxiter = 100)) # gives a warning, but doesn't stop the test from providing an estimate based on the last itteration before the warning
  
  
  # Determine an error around the point location estimate
  par.est = cbind(coef(nls.test), confint2(nls.test))
  lng.ci.upper =  par.est[1,3] 
  lng.ci.lower =  par.est[1,2]
  lat.ci.upper =  par.est[2,3] 
  lat.ci.lower =  par.est[2,2]}
  
  ,error = function(e)  {cat("ERROR :",conditionMessage(e), "\n")})
  
  
  # estimated location of the point and error
  estimated.loc <- data.frame(TestId = tests[j], PointId = rd.pt, No.Nodes = no.nodes, x.est = par.est[1,1], y.est = par.est[2,1], 
                              x.ci.lower = lng.ci.lower, x.ci.upper = lng.ci.upper, y.ci.lower = lat.ci.lower, y.ci.upper = lat.ci.upper)
  
  # Populate dataframe with results
  estimated.location_results <- rbind(estimated.location_results, estimated.loc)
  
  print(j)
  
}

# Summarize simulations for each random point
sim.summary <- estimated.location_results %>%
  dplyr::group_by(PointId) %>%
  dplyr::summarise(avg.x.est = mean(x.est),
                   sd.x.est = sd(x.est),
                   med.x.est = median(x.est),
                   min.x.est = min(x.est),
                   max.x.est = max(x.est),
                   avg.y.est = mean(y.est),
                   sd.y.est = sd(y.est),
                   med.y.est = median(y.est),
                   min.y.est = min(y.est),
                   max.y.est = max(y.est),
                   mean.no.nodes = mean(No.Nodes),
                   n.sim = n(),
                   filter = paste0("No Filter"))


# combine estimated locations with true locations
combined_results <- sim.summary %>%
  dplyr::left_join(random.points)

# Calculate difference distance between estimated and true location    
dst <- raster::pointDistance(combined_results[,15:16], combined_results[,c(2,7)], lonlat = F, allpairs = F)

# bring all together
combined_results_final <- dplyr::bind_cols(combined_results, data.frame(diff.dist = dst))

# save file
write.csv(combined_results_final, paste0(outpath, "Trilateration.Simulation_NoFilter_Summary.Results.csv"),  row.names = F)

# summarize statistics for a given filter
summary.stats <- combined_results_final %>%
  dplyr::summarise(n.pts = dplyr::n_distinct(PointId),
                   avg.no.sim = mean(n.sim),
                   avg.no.nodes = mean(mean.no.nodes),
                   avg.diff = mean(diff.dist),
                   sd.diff = sd(diff.dist),
                   lower.ci = avg.diff - qnorm(0.975)*sd.diff/sqrt(n.pts),
                   upper.ci = avg.diff + qnorm(0.975)*sd.diff/sqrt(n.pts),
                   med.diff = median(diff.dist),
                   min.diff = min(diff.dist),
                   max.dist = max(diff.dist)) %>%
  dplyr::mutate(filter = paste0("No Filter"))

# save file
write.csv(summary.stats, paste0(outpath, "Trilateration.Simulation_NoFilters_Summary.Stats.csv"),  row.names = F)

return(summary.stats)

}








###############################################
# Trilateration of Simulated Data - RSS filter
###############################################


trilateration.Sim.RSS.Filter <- function(x) {

# Create an empty dataframe to populate
summary.stats_results <- data.frame(n.pts = numeric(), avg.no.sim = numeric(), avg.no.nodes = numeric(), avg.diff = numeric(), sd.diff = numeric(),
                                    lower.ci = numeric(), upper.ci = numeric(), 
                                    med.diff = numeric(), min.diff = numeric(), max.dist = numeric(),
                                    filter = character())


for (i in 1:length(RSS.FILTER)) {
  
  # Filter data based on RSS value
  sub.rss <- x %>%
    dplyr::filter(r.RSS >= RSS.FILTER[i])
  
  # Filter list so only contains tests which had <3 nodes detecting the transmitter
  sample.size <- sub.rss %>%
    dplyr::group_by(TestId) %>%
    dplyr::summarise(n.nodes = n()) %>%
    dplyr::filter(n.nodes < 3)
  
  # Remove Tests with <3 nodes picking up a signal 
  sub.rss.red <- sub.rss %>%
    dplyr::filter(!(TestId %in% sample.size$TestId)) 
  
  # make a vector of unique trilaterations to run
  tests = unique(sub.rss.red$TestId)
  
  # Create a dataframe for output estimates
  estimated.location_results <- data.frame(TestId=character(), PointId = character(), No.Nodes = numeric(), x.est=numeric(), y.est=numeric(), 
                                           x.ci.lower =numeric(), x.ci.upper =numeric(), y.ci.lower = numeric(), y.ci.upper = numeric())
  
    for(j in 1:length(tests)) {
    
    
    # Isolate the test 
    sub.test <- sub.rss.red %>% dplyr::filter(TestId == tests[j]) 
    
    # Determine the node with the strongest RSSI value
    max.RSS <- sub.test[which.max(sub.test$r.RSS),]
    
    # Identify the random point
    rd.pt <- unique(sub.test$PointId)
    
    # Calculate no nodes for the test
    no.nodes <- dplyr::n_distinct(sub.test$NodeId)
    
    
    # To deal with potential errors where the model fails due to bad starting values using tryCatch everything you want evaluated by tryCatch goes inside {},
    # then the error will be printed but the loop will continue
    
    # Non-linear test to optimize the location of unknown signal by looking at the radius around each Node based on estimated distance and the pairwise distance between all nodes
    tryCatch( {nls.test <- nls(e.dist ~ raster::pointDistance(data.frame(x, y), c(x_solution, y_solution), lonlat = F, allpairs = T),
                               data = sub.test, start=list(x_solution=max.RSS$x, y_solution=max.RSS$y),
                               control=nls.control(warnOnly = T, minFactor=1/30000, maxiter = 100)) # gives a warning, but doesn't stop the test from providing an estimate based on the last itteration before the warning
    
    
    # Determine an error around the point location estimate
    par.est = cbind(coef(nls.test), confint2(nls.test))
    lng.ci.upper =  par.est[1,3] 
    lng.ci.lower =  par.est[1,2]
    lat.ci.upper =  par.est[2,3] 
    lat.ci.lower =  par.est[2,2]}

    ,error = function(e)  {cat("ERROR :",conditionMessage(e), "\n")})
    
    # estimated location of the point and error
    estimated.loc <- data.frame(TestId = tests[j], PointId = rd.pt, No.Nodes = no.nodes, x.est = par.est[1,1], y.est = par.est[2,1], 
                                x.ci.lower = lng.ci.lower, x.ci.upper = lng.ci.upper, y.ci.lower = lat.ci.lower, y.ci.upper = lat.ci.upper)
    
    # Populate dataframe with results
    estimated.location_results <- rbind(estimated.location_results, estimated.loc)
    
   
    
  }
  
  # Summarize simulations for each random point
  sim.summary <- estimated.location_results %>%
    dplyr::group_by(PointId) %>%
    dplyr::summarise(avg.x.est = mean(x.est),
                     sd.x.est = sd(x.est),
                     med.x.est = median(x.est),
                     min.x.est = min(x.est),
                     max.x.est = max(x.est),
                     avg.y.est = mean(y.est),
                     sd.y.est = sd(y.est),
                     med.y.est = median(y.est),
                     min.y.est = min(y.est),
                     max.y.est = max(y.est),
                     mean.no.nodes = mean(No.Nodes),
                     n.sim = n(),
                     filter = paste0("RSS"," ",RSS.FILTER[i]))
  
  
  # combine estimated locations with true locations
  combined_results <- sim.summary %>%
    dplyr::left_join(random.points)
  
  # Calculate difference distance between estimated and true location    
  dst <- raster::pointDistance(combined_results[,15:16], combined_results[,c(2,7)], lonlat = F, allpairs = F)
  
  # bring all together
  combined_results_final <- dplyr::bind_cols(combined_results, data.frame(diff.dist = dst))
  
  # save file
  write.csv(combined_results_final, paste0(outpath, "Trilateration.Simulation_Filter.RSS.", RSS.FILTER[i], "_Summary.Results.csv"),  row.names = F)
  
  # summarize statitics for a given filter
  summary.stats <- combined_results_final %>%
    dplyr::summarise(n.pts = dplyr::n_distinct(PointId),
                     avg.no.sim = mean(n.sim),
                     avg.no.nodes = mean(mean.no.nodes),
                     avg.diff = mean(diff.dist),
                     sd.diff = sd(diff.dist),
                     lower.ci = avg.diff - qnorm(0.975)*sd.diff/sqrt(n.diff),
                     upper.ci = avg.diff + qnorm(0.975)*sd.diff/sqrt(n.diff),
                     med.diff = median(diff.dist),
                     min.diff = min(diff.dist),
                     max.dist = max(diff.dist)) %>%
    dplyr::mutate(filter = paste0("RSS", " ", RSS.FILTER[i]))
  
  # write summary to empty dataframe
  summary.stats_results <- rbind(summary.stats_results, summary.stats)
  
  # save file
  write.csv(summary.stats_results, paste0(outpath, "Trilateration.Simulation_Filters.RSS_Summary.Stats.csv"),  row.names = F)
  
  print(i)
  
  
  
}

return(summary.stats_results)

}







#########################################################
# Trilateration of Simulated Data- Distance from strongest node filter 
########################################################


trilateration.Sim.Distance.Filter <- function(x) {
  
 
# Calculate distance between nodes
dist.nodes <- raster::pointDistance(nodes[,1:2], nodes[,1:2], lonlat = F, allpairs = T)

# Make matrix into a dataframe with a row for NodeId
dist.nodes_df <- data.frame(dist.nodes, row.names = nodes$NodeId)
colnames(dist.nodes_df) <- nodes$NodeId
dist.nodes_df$NodeId <- rownames(dist.nodes_df)

# make a vector unique trilaterations to run
tests = unique(x$TestId)

# Create an empty dataframe to populate
summary.stats_results <- data.frame(n.pts = numeric(), avg.no.sim = numeric(), avg.no.nodes = numeric(), avg.diff = numeric(), sd.diff = numeric(),
                                    lower.ci = numeric(), upper.ci = numeric(), 
                                    med.diff = numeric(), min.diff = numeric(), max.dist = numeric(),
                                    filter = character())

# loop through distances

for(i in 1:length(DIST.FILTER)) {
  
  # Create a dataframe for output estimates
  estimated.location_results <- data.frame(TestId=character(), PointId = character(), No.Nodes = numeric(), x.est=numeric(), y.est=numeric(), 
                                           x.ci.lower =numeric(), x.ci.upper =numeric(), y.ci.lower = numeric(), y.ci.upper = numeric())
  
  # identify the distance
  dist.no <- DIST.FILTER[i]
  
  print(i)
  
  # Loop through random points
  
  for(j in 1:length(tests)) {
    
    # Isolate the test 
    sub.test <- x %>% dplyr::filter(TestId == tests[j]) 
    
    # Determine the node with the strongest RSSI value
    max.RSS <- sub.test[which.max(sub.test$r.RSS),]
    
    # Identify the random point
    rd.pt <- unique(sub.test$PointId)
    
    # Filter matrix of node distances by node with strongest avg.RSSI
    # and get all nodes within a particular distance of the node
    nodes.test <- dist.nodes_df %>%
      dplyr::filter(NodeId %in% max.RSS$NodeId) %>%
      tidyr::gather(key = "NodeId", value = "distance", -NodeId) %>%
      dplyr::filter(distance <= dist.no)
    
    # Only keep nodes that are within the specified distance of node with strongest avg.RSSI
    sub.test.red <- sub.test %>%
      dplyr::filter(NodeId %in% nodes.test$NodeId)
    
    # Calculate no nodes for the test
    no.nodes <- dplyr::n_distinct(sub.test.red$NodeId)
    
    # If the number of nodes is not greater than 3 the rest of the loop will not be continued and the next
    # iteration of the loop is started
    if(no.nodes < 3) {
      next
    }
    
    # To deal with potential errors where the model fails due to bad starting values using tryCatch everything you want evaluated by tryCatch goes inside {},
    # then the error will be printed but the loop will continue
    
    # Non-linear test to optimize the location of unknown signal by looking at the radius around each Node based on estimated distance and the pairwise distance between all nodes
    tryCatch( {nls.test <- nls(e.dist ~ raster::pointDistance(data.frame(x, y), c(x_solution, y_solution), lonlat = F, allpairs = T),
                               data = sub.test.red, start=list(x_solution=max.RSS$x, y_solution=max.RSS$y),
                               control=nls.control(warnOnly = T, minFactor=1/30000, maxiter = 100)) # gives a warning, but doesn't stop the test from providing an estimate based on the last itteration before the warning
    
    
    # Determine an error around the point location estimate
    par.est = cbind(coef(nls.test), confint2(nls.test))
    lng.ci.upper =  par.est[1,3] 
    lng.ci.lower =  par.est[1,2]
    lat.ci.upper =  par.est[2,3] 
    lat.ci.lower =  par.est[2,2] }
    
    
    ,error = function(e)  {cat("ERROR :",conditionMessage(e), "\n")})
    
    # estimated location of the point and error
    estimated.loc <- data.frame(TestId = tests[j], PointId = rd.pt, No.Nodes = no.nodes, x.est = par.est[1,1], y.est = par.est[2,1], 
                                x.ci.lower = lng.ci.lower, x.ci.upper = lng.ci.upper, y.ci.lower = lat.ci.lower, y.ci.upper = lat.ci.upper)
    
    # Populate dataframe with results
    estimated.location_results <- rbind(estimated.location_results, estimated.loc)
    
  
  }
  
   # Summarize simulations for each random point
  sim.summary <- estimated.location_results %>%
    dplyr::group_by(PointId) %>%
    dplyr::summarise(avg.x.est = mean(x.est),
                     sd.x.est = sd(x.est),
                     med.x.est = median(x.est),
                     min.x.est = min(x.est),
                     max.x.est = max(x.est),
                     avg.y.est = mean(y.est),
                     sd.y.est = sd(y.est),
                     med.y.est = median(y.est),
                     min.y.est = min(y.est),
                     max.y.est = max(y.est),
                     mean.no.nodes = mean(No.Nodes),
                     n.sim = n(),
                     filter = paste0("Distance"," ",  DIST.FILTER[i]))
  
  # combine estimated locations with true locations
  combined_results <- sim.summary %>%
    dplyr::left_join(random.points)
  
  # Calculate difference distance between estimated and true location    
  dst <- raster::pointDistance(combined_results[,15:16], combined_results[,c(2,7)], lonlat = F, allpairs = F)
  
  # bring all together
  combined_results_final <- dplyr::bind_cols(combined_results, data.frame(diff.dist = dst))
  
  # save file
  write.csv(combined_results_final, paste0(outpath, "Trilateration.Simulation_Filter.Distance.", DIST.FILTER[i], "_Summary.Results.csv"),  row.names = F)
  
  # summarize statitics for a given filter
  summary.stats <- combined_results_final %>%
    dplyr::summarise(n.pts = dplyr::n_distinct(PointId),
                     avg.no.sim = mean(n.sim),
                     avg.no.nodes = mean(mean.no.nodes),
                     avg.diff = mean(diff.dist),
                     sd.diff = sd(diff.dist),
                     lower.ci = avg.diff - qnorm(0.975)*sd.diff/sqrt(n.pts),
                     upper.ci = avg.diff + qnorm(0.975)*sd.diff/sqrt(n.pts),
                     med.diff = median(diff.dist),
                     min.diff = min(diff.dist),
                     max.dist = max(diff.dist)) %>%
    dplyr::mutate(filter = paste0("Distance_", DIST.FILTER[i]))
  
  # populate empty dataframe
  summary.stats_results <- rbind(summary.stats_results, summary.stats)
  
  # save file
  write.csv(summary.stats_results, paste0(outpath, "Trilateration.Simulation_Filters.Distance_Summary.Stats.csv"),  row.names = F)
  

}

return(summary.stats_results)

}




