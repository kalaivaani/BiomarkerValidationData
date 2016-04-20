################################### Project: add survey data to biomarker validation data #####################################
# date: 2016-04-14
# written by: Kala Sundararajan
# biobank sample data from Rajiv Gandhi
# biobank survey data from Daniel Antflek
# LEAP BL and 3m survey data from DADOS

########################### load.R: load raw data ###############################

# serum data
Serum <- read.csv("raw_data/SerumData.csv", stringsAsFactors=FALSE)

# SF data
SF <- read.csv("raw_data/SFData.csv", stringsAsFactors=FALSE)

# biobank master list
BioMaster <- read.csv("raw_data/BiobankKneeMaster.csv", stringsAsFactors=FALSE)

# biobank data
Bio_BL <- read.csv("raw_data/BiobankKneeBL.csv", stringsAsFactors=FALSE)
Bio_3m <- read.csv("raw_data/BiobankKnee3m.csv", stringsAsFactors=FALSE)

# DADOS data
LEAP_BL <- read.csv("raw_data/LEAP_BL_20160415.csv", stringsAsFactors=FALSE)
LEAP_3m <- read.csv("raw_data/LEAP_3m_20160415.csv", stringsAsFactors=FALSE)

