EcolEvol.Manuscript_Optimizing.Trilateration

R code used in Ecology and Evoultion Manuscript - Optimizing trilateration estimates for tracking fine-scale movement of wildlife using automated radio telemetry networks. Scripts include:

1. Github_RSS_by_Distance_Calibration.R - contains code to examine the relationship between RSS values and distance for a given study area using an exponential decay model - see the beginning of the script for detailed information about files needed to run the code and output generated from the code

2. Github_TestDataset_LocalizationError.R - contains code to determine the localization error associated with RSS-based localization estimates for a node network based on a test dataset. Information generated in Github_RSS_by_Distance_Calibration.R script is necessary for this script. See the beginning of the script for detailed information about files needed to run the code and output generated from the code.

3. Github_Simulations.R - contains code to 1) simulate 3 node configurations within a 12.5 km2 grid - node spacing: 100m (n = 100 nodes), 175m (n = 64 nodes), or 250m (n = 36 nodes), 2) create simulated dataset of RSS values for a given node configuration incorporating variability in RSS values based on noise in your node network, 3) assess the error associated with RSS-based localization estimates when filters are applied prior to trilateration analysis for different node configurations. Information generated in Github_RSS_by_Distance_Calibration.R script is necessary for this script. See the beginning of the script for detailed information about files needed to run the code and output generated from the code.

4. Functions_RSS.Based.Localizations.R - Functions used to calculate code in the scripts above associated with the RSS-based Localization Paper

