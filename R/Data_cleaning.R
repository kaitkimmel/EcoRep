############################################################
#### Data pre-processing for Replicability in Ecology ####
##########################################################
## Last updated: Feb 6, 2023

## This code is used to clean 2 datasets used for analyses in Replicability in Ecology
## 1. RepDat is the file of all the estimates from papers - it needs to be cleaned to remove data
## that was later excluded after double checking and to standardize how categorical variables were
## coded
## 2. papers is the file that recorded EVERY paper that was in the journals for our sample period
## it needs to be cleaned to standardize how categorical variables were coded

## I used the here::here() function to call my data and save files. I have a folder "Data" in which all the data is stored. 
## The here function makes it so that a new working directory for other users does not have to be set. 
## The code should run easily if a Data folder is made within a new project and here() is used. 

##### Load Libraries
library(here)

#### Load data
# FULL DATASET WITH ESTIMATES
RepDat <- read.csv(here("Data/RepDataJan23.csv"), na.strings=c("","NA"))
RepDat$coefficient <- as.numeric(RepDat$coefficient) #NAs introduced because of coefficients reported as "<0.01" for example
RepDat$sample_size <- as.numeric(RepDat$sample_size) #NAs where sample size was not assigned because could not determine from paper 
RepDat$std_error <- as.numeric(RepDat$std_error) #NAs introduced because of error reported as "<0.01" for example

# DATA OF ALL PAPERS CONSIDERED
papers <- read.csv(here("Data/Paper_trackerFeb21.csv"))


##### Cleaning RepDat file
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

# get a dataframe of all the data that was eventually removed after checking data
# keeping this so we can calculate how many papers & estimates we had to remove because
# sample size was unclear
Kicked_out <- RepDat[RepDat$Keep != "yes",]
Kicked_out$reason <- Kicked_out$Checking.Notes
Kicked_out$reason[Kicked_out$Checking.Notes == "These values are leave-out-one cross validation, not model coefficients"] <- "Not model coef"
Kicked_out$reason[Kicked_out$Checking.Notes == "random effects"] <-  "Not model coef"
Kicked_out$reason[Kicked_out$Checking.Notes == "model selection"] <-  "Not model coef"
Kicked_out$reason[Kicked_out$Checking.Notes == "Random effects - remove"] <-  "Not model coef"
Kicked_out$reason[Kicked_out$Checking.Notes == "Random effects- remove"] <-  "Not model coef"
Kicked_out$reason[Kicked_out$Checking.Notes == "Random effects - remove"] <-  "Not model coef"
Kicked_out$reason[Kicked_out$Checking.Notes == "summary stats"] <-  "Not model coef"
Kicked_out$reason[Kicked_out$Checking.Notes == "unstandardized path estimate"] <-  "Not model coef"
Kicked_out$reason[Kicked_out$Checking.Notes == "unsure if stats were done on these... or just the Fis values"] <-  "Not model coef"
Kicked_out$reason[Kicked_out$Checking.Notes == "unstandardized path estimate"] <-  "Not model coef"
Kicked_out$reason[Kicked_out$Checking.Notes == "model selection - sort of, reversed order of parameters with type I sum of squares"] <- "Estimates not part of study"
Kicked_out$reason[Kicked_out$Checking.Notes == "zero-inflation parameter?"] <-  "Not model coef"
Kicked_out$reason[Kicked_out$Checking.Notes == "Top five models presented"] <-  "Not model coef"
Kicked_out$reason[Kicked_out$Checking.Notes == "random"] <-  "Not model coef"
Kicked_out$reason[Kicked_out$Checking.Notes == "not estimates and SE"] <-  "Not model coef"
Kicked_out$reason[Kicked_out$Checking.Notes == "LOO estimates for model comparison"] <-  "Not model coef"
Kicked_out$reason[Kicked_out$Checking.Notes == "full model"] <-  "Reduced model used in manuscript - full model would repeat coefficients"
Kicked_out$reason[Kicked_out$Checking.Notes == "intercepts"] <-  "Not model coef"
Kicked_out$reason[Kicked_out$Checking.Notes == "summary statistics"] <-  "Not model coef"
Kicked_out$reason[Kicked_out$Checking.Notes == "not stastics"] <-  "Not model coef"
Kicked_out$reason[Kicked_out$Checking.Notes == "not a statistical test"] <-  "Not model coef"
Kicked_out$reason[Kicked_out$Checking.Notes == "covariance part of glmm"] <-  "Not model coef"
Kicked_out$reason[Kicked_out$Checking.Notes == "rarefied dietary breadth - not statistics"] <-  "Not model coef"
Kicked_out$reason[Kicked_out$Checking.Notes == "random effect"] <-  "Not model coef"
Kicked_out$reason[Kicked_out$Checking.Notes == "model averaging"] <-  "Not model coef"
Kicked_out$reason[Kicked_out$Checking.Notes == "not sure that we can use this one either. Each coefficient is from a replicate over time."] <-  "Not model coef"


Kicked_out$reason[Kicked_out$Checking.Notes == "metaanalysis"] <- "meta-analysis"
Kicked_out$reason[Kicked_out$Checking.Notes == "meta-analysis "] <- "meta-analysis"
Kicked_out$reason[Kicked_out$Checking.Notes == "meta-analysis"] <- "meta-analysis"

Kicked_out$reason[Kicked_out$Checking.Notes == "unsure of n"]  <- "unclear n"
Kicked_out$reason[Kicked_out$Checking.Notes== "unsure of n (do not know number of points per time period)" ]  <- "unclear n"
Kicked_out$reason[Kicked_out$Checking.Notes == "unsure of n "]  <- "unclear n"
Kicked_out$reason[Kicked_out$Checking.Notes == "unsure of n for shade and open plots" ]  <- "unclear n"
Kicked_out$reason[Kicked_out$Checking.Notes == "can't find sample size" ]  <- "unclear n"
Kicked_out$reason[Kicked_out$Checking.Notes == "unsure of n for Jaccard" ]  <- "unclear n"
Kicked_out$reason[Kicked_out$Checking.Notes == "does not give number of cells used" ]  <- "unclear n"
Kicked_out$reason[Kicked_out$Checking.Notes == "unsure of n - gives total n, not split by taxa" ]  <- "unclear n"
Kicked_out$reason[Kicked_out$Checking.Notes == "unsure of n - does not split down observations by sex" ]  <- "unclear n"
Kicked_out$reason[Kicked_out$Checking.Notes == "would need to go into data file and subset out lat and long with values associated with them" ]  <- "unclear n"
Kicked_out$reason[Kicked_out$Checking.Notes == "unsure of observations per density" ]  <- "unclear n"
Kicked_out$reason[Kicked_out$Checking.Notes == "only gives average sample size per treatment, not actual sample size" ]  <- "unclear n"
Kicked_out$reason[Kicked_out$Checking.Notes == "unsure of N by date" ]  <- "unclear n"
Kicked_out$reason[Kicked_out$Checking.Notes== "unsure ofn by species and date" ]  <- "unclear n"
Kicked_out$reason[Kicked_out$Checking.Notes == "n unclear, currently assuming analyses at site level instead of plot level" ]  <- "unclear n"
Kicked_out$reason[Kicked_out$Checking.Notes == "n unclear"]  <- "unclear n"
Kicked_out$reason[Kicked_out$Checking.Notes == "only gives average sample size per treatment"]  <- "unclear n"
Kicked_out$reason[Kicked_out$Checking.Notes == "unsure of n for mortality " ]  <- "unclear n"
Kicked_out$reason[Kicked_out$Checking.Notes == "unsure of n - gives total n, not split by taxa" ]  <- "unclear n"
Kicked_out$reason[Kicked_out$Checking.Notes == "unsure of n - do not know individual trees" ]  <- "unclear n"
Kicked_out$reason[Kicked_out$Checking.Notes == "rounded effecive N values from bayesian analyses" ]  <- "unclear n"
Kicked_out$reason[Kicked_out$Checking.Notes == "unsure of n - does not split down observations by sex" ]  <- "unclear n"
Kicked_out$reason[Kicked_out$Checking.Notes == "number of radars"]  <- "unclear n"
Kicked_out$reason[Kicked_out$Checking.Notes == "unsure of n - do not give sample sizes per time period and treatment"]  <- "unclear n"
Kicked_out$reason[Kicked_out$Checking.Notes == "does not give n for growth rates"]  <- "unclear n"
Kicked_out$reason[Kicked_out$Checking.Notes == "unsure of N - does not give number of saplings or total plots sampled"]  <- "unclear n"
Kicked_out$reason[Kicked_out$Checking.Notes == "unsure of N"]  <- "unclear n"

write.csv(Kicked_out, here("Data/Kicked.csv"))

# get the dataframe with all the estimates we are keep
RepDat <- RepDat[RepDat$Keep == "yes",]

# Converting CI to standard errors
RepDat$error <- RepDat$std_error
for (i in 1:nrow(RepDat)){
  if(RepDat$error_type[i] == "CI"){
    RepDat$std_error[i] = RepDat$error[i]/1.96
  }
}

# pull out only necessary columns for data analysis
RepDat <- RepDat[,c(2:18,23:25)]
# save data
write.csv(RepDat, here("Data/CleanedDat.csv"), row.names = FALSE)

#### Cleaning papers file
# Standardize reasons to get total amount of papers used. 
papers$include[papers$include == "0" ] <- "No"
papers$include[papers$include == "No " ] <- "No"
papers$include[papers$include == "N" ] <- "No"
papers$include[papers$include == "N " ] <- "No"
papers$include[papers$include == "no" ] <- "No" 
papers$include[papers$include == "no " ] <- "No" 
papers$include[papers$include == 'maybe'] <- "No"
papers$include[papers$include == "Yes " ] <- "Yes"
papers$include[papers$include == "Y" ] <- "Yes"
papers$include[papers$include == "yes"] <- "Yes"
papers$include[papers$include == "yes " ] <- "Yes"
papers$include[papers$include == "1"] <- "Yes"
# Fixing typos in journal_id
papers$journal_id[papers$journal_id == "E "] <- "E"
papers$journal_id[papers$journal_id == "N "] <- "N"
papers <- papers[,c(1:8)]

papers$reason[papers$reason == "no_tables"] <- "No Tables" 
papers$reason[papers$reason == "no tables"]<- "No Tables"
papers$reason[papers$reason == "no tables with error"] <- "No Tables"
papers$reason[papers$reason == "no tables with coefficients"]<- "No Tables"
papers$reason[papers$reason == "no tables with coefficients "]<- "No Tables"
papers$reason[papers$reason == "no tables "]<- "No Tables"
papers$reason[papers$reason == "no_tables "] <- "No Tables"
papers$reason[papers$reason == "No table"] <- "No Tables"
papers$reason[papers$reason == "no table with coefficients"]<- "No Tables"
papers$reason[papers$reason == "no table"] <- "No Tables"


papers$reason[papers$reason == "paper_type"] <- "Paper Type"
papers$reason[papers$reason =="paper type"]<- "Paper Type"
papers$reason[papers$reason =="Paper type"]<- "Paper Type"
papers$reason[papers$reason =="Simulated Data"]<- "Paper Type"
papers$reason[papers$reason =="meta-analysis"]<- "Paper Type"
papers$reason[papers$reason =="modeling paper"]<- "Paper Type"
papers$reason[papers$reason =="Modeling"]<- "Paper Type"
papers$reason[papers$reason =="Meta-analysis"]<- "Paper Type"

papers$reason[papers$reason == "only reports standard deviations"] <- "No Error Reported"
papers$reason[papers$reason == "credible interval not included"] <- "No Error Reported"
papers$reason[papers$reason == "credible interval not included"] <- "No Error Reported"
papers$reason[papers$reason == "No table with error"] <- "No Error Reported"
papers$reason[papers$reason == "no table with error"]<- "No Error Reported"
papers$reason[papers$reason == "No table with error"]<- "No Error Reported"
papers$reason[papers$resaon == "credible intervals"]<- "No Error Reported"


# save
write.csv(papers, here("Data/CleanedPapers.csv"))

