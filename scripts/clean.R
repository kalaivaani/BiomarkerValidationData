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

# select BL and 3m data from LEAP
LEAP_BL <- LEAPKnee[,c(grep("SUBJECT_ID", names(LEAPKnee)),grep("BASELINE_", names(LEAPKnee)))]
LEAP_3m <- LEAPKnee[,c(grep("SUBJECT_ID", names(LEAPKnee)),grep("X3MONTH_", names(LEAPKnee)))]