################################### Project: add survey data to biomarker validation data #####################################
# date: 2016-04-14
# written by: Kala Sundararajan
# biobank sample data from Rajiv Gandhi
# biobank survey data from Daniel Antflek
# LEAP BL and 3m survey data from DADOS


########################### functions.R: clean/merge data ###############################

## clean biobank sample IDs
clean_bioID <- function(ID) {
  # remove 'B' and hyphens
  ID <- gsub("B|-", "", ID)
  
  # remove leading zeros
  ID <- sub("^0+", "", ID)
  
  return(ID)
}

# remove prefixes from variable names
remove_prefixes <- function(dat) {
  names(dat) <- sub("^X.+?_","",names(dat))
  return(dat)
}

