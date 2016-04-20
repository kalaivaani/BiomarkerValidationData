################################### Project: add survey data to biomarker validation data #####################################
# date: 2016-04-14
# written by: Kala Sundararajan
# biobank sample data from Rajiv Gandhi
# biobank survey data from Daniel Antflek
# LEAP BL and 3m survey data from DADOS


########################### do.R: load, clean/merge data ###############################
rm(list = ls())
library(plyr)
library(lubridate)

# load functions
source("Y:/LEAP/23. LEAP OA Data Quality/scripts/functions.R")
source("scripts/functions.R")

# load data
source("scripts/load.R")

# clean data
source("scripts/clean.R")