################################### Project: add survey data to biomarker validation data #####################################
# date: 2016-04-14
# written by: Kala Sundararajan
# biobank sample data from Rajiv Gandhi
# biobank survey data from Daniel Antflek
# LEAP BL and 3m survey data from DADOS


########################### clean.R: clean/merge data ###############################

# create clean IDs for biobank datasets
BioMaster$Study.ID <- clean_bioID(BioMaster$Study.ID)
Serum$ID <- clean_bioID(Serum$ID)
SF$ID <- clean_bioID(SF$ID)
Bio_BL$Study.ID <- clean_bioID(Bio_BL$Study.ID)
Bio_3m$Study.ID <- clean_bioID(Bio_3m$Study.ID)

# prefix serum and SF variables
names(Serum)[-1] <- paste0("Serum_",names(Serum)[-1])
names(SF)[-1] <- paste0("SF_",names(SF)[-1])

# merge Serum and SF data
biomarkers <- merge(Serum, SF, by="ID", all=T)

# subset Knee data
biomarkers <- subset(biomarkers, biomarkers$ID %in% BioMaster$Study.ID)

# merge LEAP IDs from biobank master list
BioMaster[BioMaster == "" | BioMaster=="NA" | BioMaster=="?"| BioMaster=="Pending"] <- NA
biomarkers <- merge(biomarkers, BioMaster[,c("Study.ID", "LEAP.ID")], by.x="ID", by.y="Study.ID", all.x=TRUE, all.y=FALSE)

# convert date variables in LEAP
LEAP_BL[,grep("DOB$|DOE$|DOC$", names(LEAP_BL))] <- lapply(LEAP_BL[,grep("DOB$|DOE$|DOC$", names(LEAP_BL))], ymd)
LEAP_3m[,grep("DOB$|DOE$|DOC$", names(LEAP_3m))] <- lapply(LEAP_3m[,grep("DOB$|DOE$|DOC$", names(LEAP_3m))], ymd)

# remove prefixes from LEAP data
LEAP_BL <- remove_prefixes(LEAP_BL)
LEAP_3m <- remove_prefixes(LEAP_3m)

########################### LEAP BL calculations ##############################
# calculate baseline age
LEAP_BL$BASELINE_AGE <- interval(LEAP_BL$DOB, LEAP_BL$DOE)/years(1)
LEAP_BL$BASELINE_AGE[LEAP_BL$BASELINE_AGE>100 | LEAP_BL$BASELINE_AGE < 16] <- NA

# calculate height, weight, BMI
HWBMI <- calc_BMI(LEAP_BL[grep("INPUT", names(LEAP_BL))])
colnames(HWBMI) <- paste0("BASELINE_", colnames(HWBMI))

# calculate BL WOMAC
BL_WOMAC <- WOMACscoring(LEAP_BL[,grep("KOOS", names(LEAP_BL))])
colnames(BL_WOMAC) <- paste0("BASELINE_", colnames(BL_WOMAC))

LEAP_BL <- cbind(LEAP_BL, HWBMI, BL_WOMAC)

######################## LEAP 3m calculations #################################
# calculate 3m WOMAC
x3m_WOMAC <- WOMACscoring(LEAP_3m[,grep("KOOS", names(LEAP_3m))])
colnames(x3m_WOMAC) <- paste0("3MONTH_", colnames(x3m_WOMAC))

LEAP_3m <- cbind(LEAP_3m, x3m_WOMAC)


######################## LEAP BL dataset #######################################
# subject ID, age, sex, baseline DOC (expectation), BASELINE_ETHNIC --> BASELINE_EMPLOYMENT_OTHER
BLdata <- LEAP_BL[c("SUBJECT_ID", "DOE", "GENDER", "BASELINE_AGE","BASELINE_EXPECTATIONDOC")]
names(BLdata)[names(BLdata)=="DOE"] <- "DOS"
names(BLdata)[names(BLdata)=="BASELINE_EXPECTATIONDOC"] <- "BASELINE_DOC"

# add demographics
BLdata <- cbind(BLdata, LEAP_BL[,grep("BASELINE_ETHNIC", names(LEAP_BL)):grep("BASELINE_EMPLOYMENT_OTHER", names(LEAP_BL))])

# add height, weight, BMI
BLdata <- cbind(BLdata, HWBMI)

# add CHA + smoking
BLdata <- cbind(BLdata, LEAP_BL[,grep("CHA1A", names(LEAP_BL)):grep("SMOKING", names(LEAP_BL))])

# add PCS scores
BLdata <- cbind(BLdata, LEAP_BL[,grep("BASELINE_RUMINATION", names(LEAP_BL)):grep("BASELINE_PCSTOTAL", names(LEAP_BL))])

# add expectations
BLdata <- cbind(BLdata, LEAP_BL[,grep("BASELINE_EXPECTATION1", names(LEAP_BL)):grep("BASELINE_EXPECTATION8", names(LEAP_BL))])

# add homunculus count, waist, hip, algometer variables (all blank)
BLdata[,c("BASELINE_HOMUNC_COUNT","BASELINE_CLINICAL_WAIST","BASELINE_CLINICAL_HIP",
          "BASELINE_CLINICAL_PAIN_FOREARM","BASELINE_CLINICAL_PAIN_LST","BASELINE_CLINICAL_PAIN_LJL","BASELINE_CLINICAL_PAIN_MJL",
          "BASELINE_CLINICAL_PAIN_MST")] <- NA

# add SF-12 scores
BLdata <- cbind(BLdata, LEAP_BL[,grep("BASELINE_SF12BP", names(LEAP_BL)):grep("BASELINE_SF12PCS", names(LEAP_BL))])

# add KOOS scores
BLdata <- cbind(BLdata, LEAP_BL[,grep("BASELINE_SYMPTOMSCORE", names(LEAP_BL)):grep("BASELINE_QOLSCORE", names(LEAP_BL))])

# add WOMAC scores
BLdata <- cbind(BLdata, BL_WOMAC)


################################ LEAP 3m dataset ########################################
x3mdata <- LEAP_3m[,c("SUBJECT_ID","3MONTH_SATISDOC")]
names(x3mdata)[names(x3mdata)=="3MONTH_SATISDOC"] <- "3MONTH_DOC"

# SF-12 scores
x3mdata <- cbind(x3mdata, LEAP_3m[,grep("3MONTH_SF12BP", names(LEAP_3m)):grep("3MONTH_SF12PCS", names(LEAP_3m))])

# KOOS scores
x3mdata <- cbind(x3mdata, LEAP_3m[,grep("3MONTH_SYMPTOMSCORE", names(LEAP_3m)):grep("3MONTH_QOLSCORE", names(LEAP_3m))])

# WOMAC scores
x3mdata <- cbind(x3mdata, x3m_WOMAC)

# satisfaction
x3mdata <- cbind(x3mdata, LEAP_3m[,grep("3MONTH_SATISFACTION1", names(LEAP_3m)):grep("3MONTH_SATISFACTION3C", names(LEAP_3m))])


########################### merge LEAP data ############################################
LEAPdata <- merge(BLdata, x3mdata, by="SUBJECT_ID", all=T)

########################### merge biomarker and LEAP data ##################################
biomarkers <- merge(biomarkers, LEAPdata, by.x="LEAP.ID", by.y="SUBJECT_ID", all.x=TRUE, all.y=FALSE)



########################## get BL and 3m data from biobank dataset ################################
Bio_BL[Bio_BL == "." |Bio_BL == "" | Bio_BL == " " | Bio_BL == "NA"] <- NA
Bio_BL[,c(which(names(Bio_BL)=="C18"):length(Bio_BL))] <- lapply(Bio_BL[,c(which(names(Bio_BL)=="C18"):length(Bio_BL))], function(x) as.numeric(as.character(x)))

Bio_3m[Bio_3m == "." |Bio_3m == "" | Bio_3m == " " | Bio_3m == "NA"] <- NA
Bio_3m[,c(7:55,57:118,119:194)] <- lapply(Bio_3m[,c(7:55,57:118,119:194)], function(x) as.numeric(as.character(x)))

# baseline: DOS (from master), gender, age, ethnicity, country, height, weight, BMI, education, marital status, living arrangement, 
  # employment?, CHA, sf scores (calculate), WOMAC scores (calculate), expectations, PCS scores (calculate)

# score SF-36
Bio_BL_SF36 <- score_sf36(Bio_BL[,c(1, grep("SF", names(Bio_BL)))], 1)
names(Bio_BL_SF36)[2:9] <- paste0("BASELINE_SF36",names(Bio_BL_SF36)[2:9]) 
  
Bio_3m_SF36 <- score_sf36(Bio_3m[,c(1, grep("SF", names(Bio_3m)))], 1)
names(Bio_3m_SF36)[2:9] <- paste0("BASELINE_SF36",names(Bio_3m_SF36)[2:9]) 


# score WOMAC
names(Bio_BL) <- sub("^X0WOMAC([1-5])$", "WOMAC_P\\1")

# score PCS


# BB_BL <- BioMaster[,c("Study.ID", "LEAP.ID", "MRN", "Surgery.Date")]
# BB_BL[,c("Surgery.Date", "Baseline.Date")] <- lapply(BB_BL[,c("Surgery.Date", "Baseline.Date")], dmy)
# names(BB_BL) <- c("ID", "LEAP.ID", "MRN", "DOS", "BASELINE_DOC")
# 
# # demographics and CHA
# BB_BL <- merge(BB_BL, Bio_BL[,c(1, 7:70)], by.x="ID", by.y="Study.ID", all.x=T, all.y=F)
# BB_BL$Language <- BB_BL$C17 <- NULL
# names(BB_BL)[6:19] <- c("GENDER", "BASELINE_AGE", "BASELINE_ETHNIC", "BASELINE_COUNTRY",
#                         "BASELINE_HEIGHT", "BASELINE_WEIGHT","BASELINE_BMI", "BASELINE_WAIST", "BASELINE_HIP",
#                         "BASELINE_SCHOOL", "BASELINE_LIVING", "BASELINE_MARITAL", "BASELINE_EMPLOYMENT")
# 
# 
# names(BB_BL) <- sub("^C(.+)(A|B|C)$", "BASELINE_CHA\\1\\2", names(BB_BL))
# names(BB_BL)[names(BB_BL)=="C18"] <- "BASELINE_SMOKING"
# 
# #baseline sf-36
# Bio_BL <- score_sf36(Bio_BL)
# 
# # 3 month: SF scores (calculate), WOMAC scores (calculate), satisfaction




















write.csv(biomarkers, file="processed_data/biomarkers_LEAP.csv", row.names=F, na="")
