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

score_sf36 <- function(responses,version) {
  ### copy variables
  responses$GH1 <- responses[,2]
  responses$HT <- responses[,3]
  responses$PF01 <- responses[,4]
  responses$PF02 <- responses[,5]
  responses$PF03 <- responses[,6]
  responses$PF04 <- responses[,7]
  responses$PF05 <- responses[,8]
  responses$PF06 <- responses[,9]
  responses$PF07 <- responses[,10]
  responses$PF08 <- responses[,11]
  responses$PF09 <- responses[,12]
  responses$PF10 <- responses[,13]
  responses$RP01 <- responses[,14]
  responses$RP02 <- responses[,15]
  responses$RP03 <- responses[,16]
  responses$RP04 <- responses[,17]
  responses$RE01 <- responses[,18]
  responses$RE02 <- responses[,19]
  responses$RE03 <- responses[,20]
  responses$SF1 <- responses[,21]
  responses$BP1 <- responses[,22]
  responses$BP2 <- responses[,23]
  responses$VT1 <- responses[,24]
  responses$MH1 <- responses[,25]
  responses$MH2 <- responses[,26]
  responses$MH3 <- responses[,27]
  responses$VT2 <- responses[,28]
  responses$MH4 <- responses[,29]
  responses$VT3 <- responses[,30]
  responses$MH5 <- responses[,31]
  responses$VT4 <- responses[,32]
  responses$SF2 <- responses[,33]
  responses$GH2 <- responses[,34]
  responses$GH3 <- responses[,35]
  responses$GH4 <- responses[,36]
  responses$GH5 <- responses[,37]
  
  #### version 2 coding
  if (version==2) {
    ### recode out of range values to missing
    responses$GH1[responses$GH1 < 1 | responses$GH1 > 5] <- NA
    responses$HT[responses$HT < 1 | responses$HT > 5] <- NA
    responses$PF01[responses$PF01 < 1 | responses$PF01 > 3] <- NA
    responses$PF02[responses$PF02 < 1 | responses$PF02 > 3] <- NA
    responses$PF03[responses$PF03 < 1 | responses$PF03 > 3] <- NA
    responses$PF04[responses$PF04 < 1 | responses$PF04 > 3] <- NA
    responses$PF05[responses$PF05 < 1 | responses$PF05 > 3] <- NA
    responses$PF06[responses$PF06 < 1 | responses$PF06 > 3] <- NA
    responses$PF07[responses$PF07 < 1 | responses$PF07 > 3] <- NA
    responses$PF08[responses$PF08 < 1 | responses$PF08 > 3] <- NA
    responses$PF09[responses$PF09 < 1 | responses$PF09 > 3] <- NA
    responses$PF10[responses$PF10 < 1 | responses$PF10 > 3] <- NA
    responses$RP01[responses$RP01 < 1 | responses$RP01 > 5] <- NA
    responses$RP02[responses$RP02 < 1 | responses$RP02 > 5] <- NA
    responses$RP03[responses$RP03 < 1 | responses$RP03 > 5] <- NA
    responses$RP04[responses$RP04 < 1 | responses$RP04 > 5] <- NA
    responses$RE01[responses$RE01 < 1 | responses$RE01 > 5] <- NA
    responses$RE02[responses$RE02 < 1 | responses$RE02 > 5] <- NA
    responses$RE03[responses$RE03 < 1 | responses$RE03 > 5] <- NA
    responses$SF1[responses$SF1 < 1 | responses$SF1 > 5] <- NA
    responses$BP1[responses$BP1 < 1 | responses$BP1 > 6] <- NA
    responses$BP2[responses$BP2 < 1 | responses$BP2 > 5] <- NA
    responses$VT1[responses$VT1 < 1 | responses$VT1 > 5] <- NA
    responses$MH1[responses$MH1 < 1 | responses$MH1 > 5] <- NA
    responses$MH2[responses$MH2 < 1 | responses$MH2 > 5] <- NA
    responses$MH3[responses$MH3 < 1 | responses$MH3 > 5] <- NA
    responses$VT2[responses$VT2 < 1 | responses$VT2 > 5] <- NA
    responses$MH4[responses$MH4 < 1 | responses$MH4 > 5] <- NA
    responses$VT3[responses$VT3 < 1 | responses$VT3 > 5] <- NA
    responses$MH5[responses$MH5 < 1 | responses$MH5 > 5] <- NA
    responses$VT4[responses$VT4 < 1 | responses$VT4 > 5] <- NA
    responses$SF2[responses$SF2 < 1 | responses$SF2 > 5] <- NA
    responses$GH2[responses$GH2 < 1 | responses$GH2 > 5] <- NA
    responses$GH3[responses$GH3 < 1 | responses$GH3 > 5] <- NA
    responses$GH4[responses$GH4 < 1 | responses$GH4 > 5] <- NA
    responses$GH5[responses$GH5 < 1 | responses$GH5 > 5] <- NA
    
    ## invert reversed scores
    responses$GH1 <- 6 - responses$GH1
    responses$HT <- 6 - responses$HT
    responses$BP1 <- 7 - responses$BP1
    responses$BP2 <- 6 - responses$BP2
    responses$VT1 <- 6 - responses$VT1
    responses$MH3 <- 6 - responses$MH3
    responses$VT2 <- 6 - responses$VT2
    responses$MH5 <- 6 - responses$MH5
    responses$GH3 <- 6 - responses$GH3
    responses$GH5 <- 6 - responses$GH5
    
    ### convert all to percentages
    responses$GH1 <- (responses$GH1-1)/4*100
    responses$HT <- (responses$HT-1)/4*100
    responses$PF01 <- (responses$PF01-1)/2*100
    responses$PF02 <- (responses$PF02-1)/2*100
    responses$PF03 <- (responses$PF03-1)/2*100
    responses$PF04 <- (responses$PF04-1)/2*100
    responses$PF05 <- (responses$PF05-1)/2*100
    responses$PF06 <- (responses$PF06-1)/2*100
    responses$PF07 <- (responses$PF07-1)/2*100
    responses$PF08 <- (responses$PF08-1)/2*100
    responses$PF09 <- (responses$PF09-1)/2*100
    responses$PF10 <- (responses$PF10-1)/2*100
    responses$RP01 <- (responses$RP01-1)/4*100
    responses$RP02 <- (responses$RP02-1)/4*100
    responses$RP03 <- (responses$RP03-1)/4*100
    responses$RP04 <- (responses$RP04-1)/4*100
    responses$RE01 <- (responses$RE01-1)/4*100
    responses$RE02 <- (responses$RE02-1)/4*100
    responses$RE03 <- (responses$RE03-1)/4*100
    responses$SF1 <- (responses$SF1-1)/4*100
    responses$BP1 <- (responses$BP1-1)/5*100
    responses$BP2 <- (responses$BP2-1)/4*100
    responses$VT1 <- (responses$VT1-1)/4*100
    responses$MH1 <- (responses$MH1-1)/4*100
    responses$MH2 <- (responses$MH2-1)/4*100
    responses$MH3 <- (responses$MH3-1)/4*100
    responses$VT2 <- (responses$VT2-1)/4*100
    responses$MH4 <- (responses$MH4-1)/4*100
    responses$VT3 <- (responses$VT3-1)/4*100
    responses$MH5 <- (responses$MH5-1)/4*100
    responses$VT4 <- (responses$VT4-1)/4*100
    responses$SF2 <- (responses$SF2-1)/4*100
    responses$GH2 <- (responses$GH2-1)/4*100
    responses$GH3 <- (responses$GH3-1)/4*100
    responses$GH4 <- (responses$GH4-1)/4*100
    responses$GH5 <- (responses$GH5-1)/4*100
  }
  
  ### v1 coding
  if (version == 1) {
    ### recode out of range values to missing
    responses$GH1[responses$GH1 < 1 | responses$GH1 > 5] <- NA
    responses$HT[responses$HT < 1 | responses$HT > 5] <- NA
    responses$PF01[responses$PF01 < 1 | responses$PF01 > 3] <- NA
    responses$PF02[responses$PF02 < 1 | responses$PF02 > 3] <- NA
    responses$PF03[responses$PF03 < 1 | responses$PF03 > 3] <- NA
    responses$PF04[responses$PF04 < 1 | responses$PF04 > 3] <- NA
    responses$PF05[responses$PF05 < 1 | responses$PF05 > 3] <- NA
    responses$PF06[responses$PF06 < 1 | responses$PF06 > 3] <- NA
    responses$PF07[responses$PF07 < 1 | responses$PF07 > 3] <- NA
    responses$PF08[responses$PF08 < 1 | responses$PF08 > 3] <- NA
    responses$PF09[responses$PF09 < 1 | responses$PF09 > 3] <- NA
    responses$PF10[responses$PF10 < 1 | responses$PF10 > 3] <- NA
    responses$RP01[responses$RP01 < 1 | responses$RP01 > 2] <- NA
    responses$RP02[responses$RP02 < 1 | responses$RP02 > 2] <- NA
    responses$RP03[responses$RP03 < 1 | responses$RP03 > 2] <- NA
    responses$RP04[responses$RP04 < 1 | responses$RP04 > 2] <- NA
    responses$RE01[responses$RE01 < 1 | responses$RE01 > 2] <- NA
    responses$RE02[responses$RE02 < 1 | responses$RE02 > 2] <- NA
    responses$RE03[responses$RE03 < 1 | responses$RE03 > 2] <- NA
    responses$SF1[responses$SF1 < 1 | responses$SF1 > 5] <- NA
    responses$BP1[responses$BP1 < 1 | responses$BP1 > 6] <- NA
    responses$BP2[responses$BP2 < 1 | responses$BP2 > 5] <- NA
    responses$VT1[responses$VT1 < 1 | responses$VT1 > 6] <- NA
    responses$MH1[responses$MH1 < 1 | responses$MH1 > 6] <- NA
    responses$MH2[responses$MH2 < 1 | responses$MH2 > 6] <- NA
    responses$MH3[responses$MH3 < 1 | responses$MH3 > 6] <- NA
    responses$VT2[responses$VT2 < 1 | responses$VT2 > 6] <- NA
    responses$MH4[responses$MH4 < 1 | responses$MH4 > 6] <- NA
    responses$VT3[responses$VT3 < 1 | responses$VT3 > 6] <- NA
    responses$MH5[responses$MH5 < 1 | responses$MH5 > 6] <- NA
    responses$VT4[responses$VT4 < 1 | responses$VT4 > 5] <- NA
    responses$SF2[responses$SF2 < 1 | responses$SF2 > 5] <- NA
    responses$GH2[responses$GH2 < 1 | responses$GH2 > 5] <- NA
    responses$GH3[responses$GH3 < 1 | responses$GH3 > 5] <- NA
    responses$GH4[responses$GH4 < 1 | responses$GH4 > 5] <- NA
    responses$GH5[responses$GH5 < 1 | responses$GH5 > 5] <- NA
    
    ## invert reversed scores
    responses$GH1 <- 6 - responses$GH1
    responses$HT <- 6 - responses$HT
    responses$SF1 <- 6 - responses$SF1
    responses$BP1 <- 7 - responses$BP1
    responses$BP2 <- 6 - responses$BP2
    responses$VT1 <- 7 - responses$VT1
    responses$MH3 <- 7 - responses$MH3
    responses$VT2 <- 7 - responses$VT2
    responses$MH5 <- 7 - responses$MH5
    responses$GH3 <- 6 - responses$GH3
    responses$GH5 <- 6 - responses$GH5
    
    ### convert all to percentages
    responses$GH1 <- (responses$GH1-1)/4*100
    responses$HT <- (responses$HT-1)/4*100
    responses$PF01 <- (responses$PF01-1)/2*100
    responses$PF02 <- (responses$PF02-1)/2*100
    responses$PF03 <- (responses$PF03-1)/2*100
    responses$PF04 <- (responses$PF04-1)/2*100
    responses$PF05 <- (responses$PF05-1)/2*100
    responses$PF06 <- (responses$PF06-1)/2*100
    responses$PF07 <- (responses$PF07-1)/2*100
    responses$PF08 <- (responses$PF08-1)/2*100
    responses$PF09 <- (responses$PF09-1)/2*100
    responses$PF10 <- (responses$PF10-1)/2*100
    responses$RP01 <- (responses$RP01-1)*100
    responses$RP02 <- (responses$RP02-1)*100
    responses$RP03 <- (responses$RP03-1)*100
    responses$RP04 <- (responses$RP04-1)*100
    responses$RE01 <- (responses$RE01-1)*100
    responses$RE02 <- (responses$RE02-1)*100
    responses$RE03 <- (responses$RE03-1)*100
    responses$SF1 <- (responses$SF1-1)/4*100
    responses$BP1 <- (responses$BP1-1)/5*100
    responses$BP2 <- (responses$BP2-1)/4*100
    responses$VT1 <- (responses$VT1-1)/5*100
    responses$MH1 <- (responses$MH1-1)/5*100
    responses$MH2 <- (responses$MH2-1)/5*100
    responses$MH3 <- (responses$MH3-1)/5*100
    responses$VT2 <- (responses$VT2-1)/5*100
    responses$MH4 <- (responses$MH4-1)/5*100
    responses$VT3 <- (responses$VT3-1)/5*100
    responses$MH5 <- (responses$MH5-1)/5*100
    responses$VT4 <- (responses$VT4-1)/4*100
    responses$SF2 <- (responses$SF2-1)/4*100
    responses$GH2 <- (responses$GH2-1)/4*100
    responses$GH3 <- (responses$GH3-1)/4*100
    responses$GH4 <- (responses$GH4-1)/4*100
    responses$GH5 <- (responses$GH5-1)/4*100
  } 
  
  ### calculate raw subscores
  responses$PF <- rowMeans(responses[,c("PF01","PF02","PF03","PF04","PF05","PF06","PF07","PF08","PF09","PF10")], na.rm=TRUE)
  
  responses$RP <- rowMeans(responses[,c("RP01","RP02","RP03","RP04")], na.rm=TRUE)
  
  responses$BP <- rowMeans(responses[,c("BP1","BP2")], na.rm=TRUE)
  
  responses$GH <- rowMeans(responses[,c("GH1","GH2","GH3","GH4","GH5")], na.rm=TRUE)
  
  responses$VT <- rowMeans(responses[,c("VT1","VT2","VT3","VT4")], na.rm=TRUE)
  
  responses$SF <- rowMeans(responses[,c("SF1","SF2")], na.rm=TRUE)
  
  responses$RE <- rowMeans(responses[,c("RE01","RE02","RE03")], na.rm=TRUE)
  
  responses$MH <- rowMeans(responses[,c("MH1","MH2","MH3","MH4","MH5")], na.rm=TRUE)
  
  ## recode NaN to NA
  responses$PF[is.nan(responses$PF)] <- NA
  responses$RP[is.nan(responses$RP)] <- NA
  responses$BP[is.nan(responses$BP)] <- NA
  responses$GH[is.nan(responses$GH)] <- NA
  responses$VT[is.nan(responses$VT)] <- NA
  responses$SF[is.nan(responses$SF)] <- NA
  responses$RE[is.nan(responses$RE)] <- NA
  responses$MH[is.nan(responses$MH)] <- NA
  
  ## calculate norm-based scores
  
  responses$PF_Z <- (responses$PF - 81.16997) / 29.10361
  responses$RP_Z <- (responses$RP - 80.50434) / 27.11497
  responses$BP_Z <- (responses$BP - 81.74681) / 24.53386
  responses$GH_Z <- (responses$GH - 72.19753) / 23.19289
  responses$VT_Z <- (responses$VT - 55.59788) / 24.84266
  responses$SF_Z <- (responses$SF - 83.73762) / 24.75248
  responses$RE_Z <- (responses$RE - 86.4124) / 22.34803
  responses$MH_Z <- (responses$MH - 70.18443) / 20.4918
  
  responses$PCS_Z = (responses$PF_Z * 0.42402) + 
    (responses$RP_Z * 0.35119) + 
    (responses$BP_Z * 0.31754) +
    (responses$GH_Z * 0.24954) + 
    (responses$VT_Z * 0.02877) + 
    (responses$SF_Z * -.00753) +
    (responses$RE_Z * -.19206) + 
    (responses$MH_Z * -.22069) 
  
  responses$MCS_Z = (responses$PF_Z * -.22999) + 
    (responses$RP_Z * -.12329) + 
    (responses$BP_Z * -.09731) +
    (responses$GH_Z * -.01571) + 
    (responses$VT_Z * 0.23534) + 
    (responses$SF_Z * 0.26876) +
    (responses$RE_Z * 0.43407) + 
    (responses$MH_Z * 0.48581) 
  
  responses$PF_NBS <- 50 + (responses$PF_Z * 10)
  responses$RP_NBS <- 50 + (responses$RP_Z * 10)
  responses$BP_NBS <- 50 + (responses$BP_Z * 10)
  responses$GH_NBS <- 50 + (responses$GH_Z * 10)
  responses$VT_NBS <- 50 + (responses$VT_Z * 10)
  responses$SF_NBS <- 50 + (responses$SF_Z * 10)
  responses$RE_NBS <- 50 + (responses$RE_Z * 10)
  responses$MH_NBS <- 50 + (responses$MH_Z * 10)
  responses$PCS_NBS <- 50 + (responses$PCS_Z * 10)
  responses$MCS_NBS <- 50 + (responses$MCS_Z * 10)
  
  x1 <- responses[c(1:37)]
  x2 <- responses[c("PF","RP","BP","GH","VT","SF","RE","MH","PF_NBS","RP_NBS","BP_NBS","GH_NBS","VT_NBS","SF_NBS","RE_NBS","MH_NBS","PCS_NBS","MCS_NBS")]
  
  return(cbind(ID=x1[,1],x2))
}