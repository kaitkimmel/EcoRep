library(here)
#### FULL DATASET WITH ESTIMATES
RepDat <- read.csv(here("Data/RepDataFeb21.csv"), na.strings=c("","NA"))
RepDat$coefficient <- as.numeric(RepDat$coefficient) #NAs introduced because of coefficients reported as "<0.01" for example
RepDat$sample_size <- as.numeric(RepDat$sample_size) #NAs where sample size was not assigned because could not determine from paper 
RepDat$std_error <- as.numeric(RepDat$std_error) #NAs introduced because of error reported as "<0.01" for example

#### DATA OF ALL PAPERS CONSIDERED
papers <- read.csv(here("Data/Paper_trackerFeb21.csv"))

# making entries congruent in categorical columns

RepDat$table_loc[RepDat$table_loc == "main "|RepDat$table_loc == "Main"|RepDat$table_loc == "Main "] <- "main"
RepDat$table_loc[RepDat$table_loc == "supplemental"|RepDat$table_loc == "Supplemental"|
                   RepDat$table_loc == "Supplement "|RepDat$table_loc == "Supplement"|
                   RepDat$table_loc == "Extended Data"|RepDat$table_loc == "Supplement S2"|
                   RepDat$table_loc == "Supplement S3"|RepDat$table_loc == "Supplement S4"|
                   RepDat$table_loc == "Supplement S5"|RepDat$table_loc == "Supplement S6" |
                   RepDat$table_loc == "extended data" | RepDat$table_loc == "Supplemental "] <- "supplement"

RepDat$error_type[RepDat$error_type == "CI-upper limit"] <- "CI"

RepDat$study_type[RepDat$study_type == "Experimental" | RepDat$study_type == "Experimantal"|
                    RepDat$study_type == "Exp"|RepDat$study_type == "exp"] <- "experimental"
RepDat$study_type[RepDat$study_type == "Observational"|RepDat$study_type == "Obs"|
                    RepDat$study_type == "obs"] <- "observational"
RepDat$study_type[RepDat$study_type == "observational and experimental" | 
                    RepDat$study_type == "Observational and Experimental" |
                    RepDat$study_type == "Observational/Experimental" |
                    RepDat$study_type == "Experimental and observational" |
                    RepDat$study_type == "Experimental and Observational" |
                    RepDat$study_type == "Experimental & observational" |
                    RepDat$study_type == "Experimantal & observational " |
                    RepDat$study_type == "experimental & observational" | 
                    RepDat$study_type == "experimental and observational"] <- "combined"
RepDat$initials[which(is.na(RepDat$initials))] <- "UNK"


Kicked_out <- RepDat[RepDat$Keep != "yes",]
Kicked_out <- Kicked_out[-which(is.na(Kicked_out$checked)),] #several rows of all NA values
write.csv(Kicked_out, here("Data/Kicked.csv"))

RepDat <- RepDat[RepDat$Keep == "yes",]
RepDat <- RepDat[-which(is.na(RepDat$journal_id)),]

## Converting CI to standard errors
RepDat$error <- RepDat$std_error

for (i in 1:nrow(RepDat)){
  if(RepDat$error_type[i] == "CI"){
    RepDat$std_error[i] = RepDat$error[i]/1.96
  }
}

RepDat <- RepDat[,c(2:18,23)]

write.csv(RepDat, here("Data/CleanedDat.csv"))


# Standardize reasons to get total amount of papers used. 
papers$include[papers$include == "0" | papers$include == "No " |
                 papers$include == "N" |papers$include == "N " |
                 papers$include == "no" | papers$include == "no "] <- "No"
papers$include[papers$include == "Yes " | papers$include == "Y" |
                 papers$include == "yes" | papers$include == "yes " |
                 papers$include == "1"] <- "Yes"
papers$journal_id[papers$journal_id == "E "] <- "E"
papers$journal_id[papers$journal_id == "N "] <- "N"


write.csv(papers, here("Data/CleanedPapers.csv"))
