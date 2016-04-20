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

# mode function
modestat <- function(x) {
  return(unique(x)[which.max(tabulate(match(x, unique(x))))])
}

#### BMI calculation function ####
# argument: data.frame with input height (ft, inches, m, cm) and weight (kg, lb)
# validate all inputs
# return calculated height (m), weight (kg), BMI (kg/m2)

calc_BMI <- function(HWdata) {
  HWdata <- as.data.frame(HWdata)
  names(HWdata) <- c("ft", "inc", "m", "cm", "lb", "kg")
  HWdata <- lapply(HWdata, function(x) as.numeric(as.character(x)))
  
  # validate inches
  HWdata$inc <- ifelse(HWdata$inc <= 0 & HWdata$inc < 12, HWdata$inc, NA)
  
  # calculate height in inches
  HWdata$inches <- ifelse(!is.na(HWdata$ft), HWdata$ft*12, 0) + ifelse(!is.na(HWdata$inc), HWdata$inc, 0)
  
  # validate inches: 42-96 inches
  HWdata$inches <- ifelse(HWdata$inches > 42 & HWdata$inches < 96, HWdata$inches, NA)
  
  # calculate height in m
  HWdata$metres <- ifelse(!is.na(HWdata$m), HWdata$m, 0) + ifelse(!is.na(HWdata$cm), HWdata$cm/100, 0)
  
  # validate metres: 1.0668-2.4384 m
  HWdata$metres <- ifelse(HWdata$metres > 1.0668 & HWdata$metres < 2.4384, HWdata$metres, 0)
  
  # calculate height (metres)
  HWdata$height <- ifelse(!is.na(HWdata$inches), HWdata$inches*0.0254, HWdata$metres)
  
  # valid kg: >22, <227
  HWdata$kg <- ifelse(HWdata$kg > 22 & HWdata$kg < 227, HWdata$kg, NA)
  
  # valid lb: >50, <500
  HWdata$lb <- ifelse(HWdata$lb > 50 & HWdata$lb < 500, HWdata$lb, NA)
  
  # calculate weight (kg)
  HWdata$weight <- ifelse(!is.na(HWdata$lb), HWdata$lb/2.2, HWdata$kg)
  
  ### calculate BMI ###
  BMI <- ifelse(!is.na(HWdata$height) & !is.na(HWdata$weight), HWdata$weight/(HWdata$height^2), NA)
  BMI <- ifelse(BMI < 15 | BMI > 75, NA, BMI)
    
  return(cbind(HEIGHT=ifelse(HWdata$height>0,HWdata$height,NA), WEIGHT=ifelse(HWdata$weight>0, HWdata$weight, NA), BMI=BMI))
}

# WOMAC scoring function: give KOOS responses with number suffixes
WOMACscoring <- function(KOOS) {
  # pain: sum p5-p9, 2+ missing responses -> invalid
  # stiffness: sum s6-s7, 2 missing responses -> invalid
  # function: sum a1-a17, 5+ missing responses -> invalid
  
  PainItems <- KOOS[,c(grep("_P5", names(KOOS)):grep("_P9", names(KOOS)))]
  StiffItems <- KOOS[,c(grep("_S6", names(KOOS)):grep("_S7", names(KOOS)))]
  FuncItems <- KOOS[,c(grep("_A1$", names(KOOS)):grep("_A17", names(KOOS)))]
  
  WOMAC_PAIN <- ifelse(apply(PainItems, 1, function(x) table(is.na(x))["FALSE"])==5, apply(PainItems,1,sum),
                       ifelse(apply(PainItems, 1, function(x) table(is.na(x))["FALSE"])==4, apply(PainItems,1,sum, na.rm=TRUE) + apply(PainItems,1,mean, na.rm=TRUE),
                              NA))
  
  WOMAC_STIFFNESS <- ifelse(apply(StiffItems, 1, function(x) table(is.na(x))["FALSE"])==2, apply(StiffItems,1,sum),
                            ifelse(apply(StiffItems, 1, function(x) table(is.na(x))["FALSE"])==1, apply(StiffItems,1,sum, na.rm=TRUE)*2,
                                   NA))
  
  WOMAC_FUNCTION <- ifelse(apply(FuncItems, 1, function(x) table(is.na(x))["FALSE"])==17, apply(FuncItems,1,sum),
                           ifelse(apply(FuncItems, 1, function(x) table(is.na(x))["FALSE"])>=13, apply(FuncItems,1,sum, na.rm=TRUE) + (17-apply(FuncItems, 1, function(x) table(is.na(x))["FALSE"])) * apply(FuncItems,1,mean, na.rm=TRUE),
                                  NA))
  
  return(cbind(WOMAC_PAIN,WOMAC_STIFFNESS,WOMAC_FUNCTION))
}