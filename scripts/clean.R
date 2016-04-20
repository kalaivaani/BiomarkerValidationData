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

#### LEAP BL calculations ####
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

#### LEAP 3m calculations #####
# calculate 3m WOMAC
x3m_WOMAC <- WOMACscoring(LEAP_3m[,grep("KOOS", names(LEAP_3m))])
colnames(x3m_WOMAC) <- paste0("3MONTH_", colnames(x3m_WOMAC))

LEAP_3m <- cbind(LEAP_3m, x3m_WOMAC)


##### LEAP BL dataset ######
# subject ID, age, sex, baseline DOC (expectation), BASELINE_ETHNIC --> BASELINE_EMPLOYMENT_OTHER
BLdata <- LEAP_BL[,c("SUBJECT_ID")]