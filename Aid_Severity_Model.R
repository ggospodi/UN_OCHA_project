# This is a RF Forest Exploring the Predictors to the Variable Degree
# author: Georgi D. Gospodinov
# date: "September 20, 2015"
# 
# Data Sources:
#
# https://data.hdx.rwlabs.org/dataset/scnepal-agency-data
# master_hlcit.csv
# 
# Relevant materials and data can be found at:
# 
# https://www.dropbox.com/sh/tb9854hzcof7x23/AACEDTGk8EmYQ6r4ukSFLBspa?dl = 0
#
# in the folder /Displacement Network Model/
#
# 
#
# LOAD PACKAGES
library(plyr)
library(dplyr)
library(igraph)
library(RColorBrewer)
#
#
#
# SET FILE SOURCE PATH
DIR <- "/Users/ggospodinov/Desktop/UN_OCHA_project/data/"
#
#
# DEFINE FUNCTIONS
#
#
#
# WRITE OBJECT FUNCTION
writeObj <- function(obj, file_name) {
  save(obj, file=file_name)
  return(file_name)
}

# READ OBJECT FUNCTION
readObj <- function(file_name) {
  obj_name <- load(file_name)
  obj <- get(obj_name)
  return(obj)
}



#
#
#
#
#
#
# READ IN THE DATA
#
#
#
#
#
#


# READ IN AID AND SEVERITY COMPLETE DATA TABLE
aid_sev <- readObj(file_name = paste0(DIR,"aid_and_severity.df"))

# READ MODELING TABLE
aid_sev_modeling <- readObj(file_name=paste0(DIR,"aid_sev_modeling.df"))
