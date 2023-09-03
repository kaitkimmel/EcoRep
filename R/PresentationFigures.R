##############################
### PRESENTATION FIGURES ####
#############################


#####################################
## Ecology Replicability Analyses ##
####################################
## Authors: Kaitlin Kimmel, Meghan Avolio, and Paul Ferraro
## Code written by: Kaitlin Kimmel (kaitlinakimmel@gmail.com)
## Code adapted from Pallavi Shukla's Stata code

#### Load libraries
library(here)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)

#### set seed
set.seed(2617)

### load in  cleaned data 
# PRIMARY DATASET WITH ESTIMATES
df <- read.csv(here("Data", "CleanedDat.csv")) 
# ALL PAPERS EXPLORED
papers <- read.csv(here("Data", "CleanedPapers.csv"), row.names = 1)
# DATA THAT WAS REMOVED FROM THE ESTIMATES FILE AFTER CHECKING
kicked <- read.csv(here("Data", "Kicked.csv"), row.names = 1)
# DATA FROM SURVEY RESULTS
survey <- read.csv(here("Data","RepEco_Survey.csv"))


##### GETTING DATA READY FOR ANALYSES
# create a unique for each paper in the dataset (uid)
uid_papers <- unique(df[,c(1:5)])
uid_papers$uid <- seq(1:nrow(uid_papers))
df <- merge(df,uid_papers) # merge uids with full dataset


# Get rid of estimates with SE of 0
nrow(df[df$std_error == 0,]) #881
df <- df[which(df$std_error != 0),]
# Get rid of coefficients with NA values - these were reported in the manuscripts as <0.001
df <- df[-which(is.na(df$coefficient)),]

################################################
#### De-round estimates and create weights ####
###############################################
# Derounding as in Brodeur et al 2016
## function to find number of decimal places for derounding
## from stackoverflow user: fvfaleiro

countDecimalPlaces <- function(x) {
  if ((x %% 1) != 0) {
    strs <- strsplit(as.character(format(x, scientific = F)), "\\.")
    n <- nchar(strs[[1]][2])
  } else {
    n <- 0
  }
  return(n) 
}

# derounding recorded coefficients and standard errors
for(i in 1:nrow(df)){
  y <- countDecimalPlaces(df$coefficient[i])
  df$coefficient_sm[i] <- runif(1,df$coefficient[i] - (0.5* 10^(-y)), df$coefficient[i] + (0.5* 10^(-y)))
  z <- countDecimalPlaces(df$std_error[i])
  df$std_error_sm[i] <- runif(1,df$std_error[i] - (0.5* 10^(-z)), df$std_error[i] + (0.5* 10^(-z)))
}


# calculate t-stats
df$tstat_sm <- df$coefficient_sm/df$std_error_sm
df$abs_tstat_sm <- abs(df$tstat_sm)

quantile(df$abs_tstat_sm, c(.01,.95,.99))

## Get rid of entries with t-stats above the 99th percentile (255 ESTIMATES TOTAL)
df <-df[df$abs_tstat_sm < 97,]

# generate weights
# count the number of observations per article
obs_by_article <- plyr::count(df, vars = "uid")
names(obs_by_article)[2] <- "obs_by_article"
# count the number of observations per table per article
obs_by_table <- plyr::count(df, vars = c("uid", "table_no")) 
names(obs_by_table)[3] <- "obs_by_table"

# count the number of tables in each article
temp_df <- df[,c("uid", "table_no")]
temp_df <- unique(temp_df)
temp_df <- temp_df %>% group_by(uid) %>% summarize(count=n())
names(temp_df)[2] <- "tab_count"

# merge observations by article, observations by table, and number of tables into main dataframe
df<- merge(df, obs_by_article)
df <- merge(df, obs_by_table)
df <- merge(df, temp_df)

# calculate weights - we want to weight articles by the inverse of the number of observations
df$weight_article <- 1/df$obs_by_article
df$weight_table <- (1/df$obs_by_table) * (1/df$tab_count)

## CALCULATE PARTIAL CORRELATION COEFFICIENTS
## PCC = TSTAT/SQRT(TSTAT^2 + RES. DF)
df$pcc <- df$tstat_sm/sqrt((df$tstat_sm^2) + df$sample_size)
df$abs_pcc <- abs(df$pcc)
df$SE_pcc <- df$pcc/df$tstat_sm
df$SE_sq <- df$SE_pcc^2
df$precision_sq <- 1/df$SE_sq


#### WLS_FE estimation and under-powered calculations ####
# get rid of one article with more than 6000 estimates
df1 <- df[df$obs_by_article <6000,]
# calculate threshold and values
# threshold for 80% power
df1$WLS_threshold <- weighted.mean(df1$abs_pcc, df1$precision_sq)/2.8 #2.8 for .8, 2.63 for 0.75 and 2.21 for 0.6
df1$WLS_value <- weighted.mean(df1$abs_pcc, df1$precision_sq)
# threshold for 75% power
df1$WLS_threshold.75 <- weighted.mean(df1$abs_pcc, df1$precision_sq)/2.63
# threshold for 60% power
df1$WLS_threshold.60 <- weighted.mean(df1$abs_pcc, df1$precision_sq)/2.21

# create columns to say if estimate is adequately powered at that cutoff
df1$powered <- NA
df1$powered.75 <- NA
df1$powered.6 <- NA
for(i in 1:nrow(df1)){
  if(df1$SE_pcc[i] <= df1$WLS_threshold[i]){
    df1$powered[i] = 1
  } else{
    df1$powered[i] = 0
  }
}
sum(df1$powered)/nrow(df1) # proportion of estimates that make the 80% cutoff

for(i in 1:nrow(df1)){
  if(df1$SE_pcc[i] <= df1$WLS_threshold.75[i]){
    df1$powered.75[i] = 1
  } else{
    df1$powered.75[i] = 0
  }
}
for(i in 1:nrow(df1)){
  if(df1$SE_pcc[i] <= df1$WLS_threshold.6[i]){
    df1$powered.6[i] = 1
  } else{
    df1$powered.6[i] = 0
  }
}
sum(df1$powered.75)/nrow(df1)# proportion of estimates that make the 75% cutoff
sum(df1$powered.6)/nrow(df1)# proportion of estimates that make the 60% cutoff


####################
##### Figure 1 ####
###################
ggsave(
  here("Figures", "Power_full_Sortee.png"),
  ggplot(df1, aes(x = SE_pcc)) +
    geom_histogram(aes(y = (..count..)/sum(..count..)), bins = 50,color = "black", fill = "gray") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), n.breaks = 6) + 
    geom_vline(xintercept = weighted.mean(df1$abs_pcc, df1$precision_sq)/2.8, color = "red", lty = "dashed") +
    annotate(geom="text", x=(weighted.mean(df1$abs_pcc, df1$precision_sq)/2.8)+.2, y=.115, label="Under-powered estimates", size = 3) +
    geom_segment(x =  weighted.mean(df1$abs_pcc, df1$precision_sq)/2.8, y = .115,
                 xend =  weighted.mean(df1$abs_pcc, df1$precision_sq)/2.8 +.05, yend = .115,
                 lineend = "round", linejoin = "round", size = 1, arrow = arrow(length = unit(0.07, "inches"))) +
    labs(x = "Standard Error of Partial Correlation Coefficient", y = "Percentage of Estimates") + 
    theme_pubclean()+ 
    theme(axis.title = element_text(face = "bold"), text = element_text(size = 14)), 
  width = 5,
  height = 3,
  dpi = 300
)


####################################
##### Range of PCC Thresholds ######
###### Supplemental Figure ########
###################################

# Create a dataframe 
newdf <- data.frame("threshold"= seq(0.01, 0.2, by = 0.01), "no_estimates" = NA)
# count the number of estimates that would fall into the hypothetical threshold
for(i in 1:nrow(newdf)){
  newdf$no_estimates[i] <- sum(df1$SE_pcc >= newdf$threshold[i])
}

# calculate the frequency 
newdf$freq <- newdf$no_estimates/19006
newdf$threshold <- as.factor(newdf$threshold)

ggsave(
  here("Figures", "Supp_thresholdrange.png"),
  ggplot(aes(x = threshold, y = freq), data = newdf) + 
    geom_bar(stat = "identity", color = "black") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(y = "Underpowered (% of estimates)", x = "Value of WLS-FE PCC") +
    theme_pubclean() + 
    theme(axis.title = element_text(face = "bold"), text = element_text(size = 12)), 
  width = 8,
  height = 4,
  dpi = 300
)

######################
#### Median power ####
#####################
#"Median power for a given area of research can then be calculated from the cumulative normal 
#probability of the difference between 1.96 and the absolute value of an estimate of the true 
#effect divided by the median standard error

# full 
WLS_threshold_median = unique(df1$WLS_value)/median(df1$SE_pcc)
1-pnorm(1.96-WLS_threshold_median)

#############################
##### Exaggeration Bias ####
############################

#step 1: calculate WLS-FE on that sub-sample of the research record that is adequately powered 
#This yields: weighted average of the adequately powered estimator (WAAP) 

# select the estimates that are adequately powered
df1.p <- df1[df1$powered == 1,]
# calculate WAAP
df1$WAAP <- weighted.mean(df1.p$abs_pcc, df1.p$precision_sq)
# select estimates that are not adequately powered
df1.u <- df1[df1$powered == 0,]

# step 2: calculate the exaggeration bias of underpowered estimates
df1.u$main_pcc_waap_diff = (df1.u$abs_pcc/df1.u$WAAP) 

# categorize exaggeration bias
df1.u$category <- NA
for (i in 1:nrow(df1.u)){
  if(df1.u$main_pcc_waap_diff[i] < 1){
    df1.u$category[i] <- "Deflation"
  }
  if(df1.u$main_pcc_waap_diff[i] >= 1 & df1.u$main_pcc_waap_diff[i]< 2){
    df1.u$category[i] <- "0-100%"
  }
  # if(df1.u$main_pcc_waap_diff[i] >= 2 & df1.u$main_pcc_waap_diff[i]<3){
  #   df1.u$category[i] <- "100-300%"
  # }
  if(df1.u$main_pcc_waap_diff[i] >= 2 & df1.u$main_pcc_waap_diff[i]<4){
    df1.u$category[i] <- "100-300%"
  }
  if(df1.u$main_pcc_waap_diff[i] >=4){
    df1.u$category[i] <- "300%+"
  }
  
}

#### manipulating data for creating bar graphs
category_counts1 <- plyr::count(df1.u$category)
names(category_counts1) <- c("category", "frequency")
category_counts1$percentage_of_studies <- category_counts1$frequency/nrow(df1.u)
category_counts1$category <- factor(category_counts1$category, 
                                    levels = c("Deflation", "none", "0-100%", "100-300%", "300%+"))


#####################
##### Figure 2 #####
####################
ggsave(here("Figures", "Exaggeration_full.png"),
       ggplot(aes(x = category, y = percentage_of_studies), data = category_counts1) + 
         geom_bar(stat = "identity", color = "black", fill = "gray") + 
         scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
         labs(y = "Percentage of Estimates", x = "Exaggeration Bias") +
         theme_pubclean() +  theme(text = element_text(size = 14), axis.title = element_text(face = "bold")),
       width = 5,
       height = 3,
       dpi = 300
)


###########################
#### p-hacking curves ####
##########################

##########################
#### Figure 3 graphs ####
#########################

# full sample - weighted by estimates per table per article
ggsave(here("Figures", "full_sample_phacking_SORTEE.png"),
ggplot(data = df1[df1$abs_tstat_sm <10,]) +
  geom_histogram(aes(x = abs_tstat_sm,  y = ..density.., weight = weight_table), 
                 bins = 50, fill = "gray") + 
  geom_density(aes(x = abs_tstat_sm, weight = weight_table), color = "black", kernel = "epanechnikov", bw = 0.2) + 
  scale_x_continuous(breaks = c(0,1.96,10)) +
  ylim(c(0,.35)) + 
  geom_segment(x =  1.9, y = .3, xend =  1.9, yend = .25,
               lineend = "round", linejoin = "round", size = .5, arrow = arrow(length = unit(0.07, "inches")), color = "red") +
  theme_pubclean() + theme(text = element_text(size = 12), axis.title = element_text(face = "bold")) + 
  labs(x = "t-statistic", y = "Density"))


nrow(df1[df1$abs_tstat_sm <10,])
# estimates presented in main text - weighted by estimates per table per article
ggsave(here("Figures", "maintext_phacking_SORTEE.png"),
ggplot(data = df3[df3$abs_tstat_sm <10,]) +
  geom_histogram(aes(x = abs_tstat_sm,  y = ..density.., weight = weight_table), 
                 bins = 50, fill = "gray") + 
  geom_density(aes(x = abs_tstat_sm, weight = weight_table), color = "black", kernel = "epanechnikov", bw = 0.2) + 
  scale_x_continuous(breaks = c(0, 1.96,10)) +
  ylim(c(0,.35)) + 
  geom_segment(x =  1.9, y = .3, xend =  1.9, yend = .25,
               lineend = "round", linejoin = "round", size = .5, arrow = arrow(length = unit(0.07, "inches")), color = "red") +
  theme_pubclean() + theme(text = element_text(size = 12), axis.title = element_text(face = "bold")) + 
  labs(x = "t-statistic", y = "Density"))

nrow(df3[df3$abs_tstat_sm <10,])
# estimates in supplemental text - weighted by estimates per table per article
ggsave(here("Figures", "supp_phacking_SORTEE.png"),
ggplot(data = df1[df1$abs_tstat_sm <10 & df1$table_loc == "supplement",]) +
  geom_histogram(aes(x = abs_tstat_sm,  y = ..density.., weight = weight_table), 
                 bins = 50, fill = "gray") + 
  geom_density(aes(x = abs_tstat_sm, weight = weight_table), color = "black", kernel = "epanechnikov") + 
  scale_x_continuous(breaks = c(0,1.96,10)) +
  ylim(c(0,.35)) + 
  geom_segment(x =  1.9, y = .3, xend =  1.9, yend = .25,
               lineend = "round", linejoin = "round", size = .5, arrow = arrow(length = unit(0.07, "inches")), color = "red") +
  theme_pubclean() + theme(text = element_text(size = 12), axis.title = element_text(face = "bold")) + 
  labs(x = "t-statistic", y = "Density"))
nrow(df1[df1$abs_tstat_sm <10 & df1$table_loc == "supplement",])

####################################
## Multiple hypothesis testing #####
####################################
# Create a smaller dataframe with only multiple hypothesis testing, code availability, and data ability 
multhyp <- unique(df[,c(1,4,7,3,14,15,17,18,19)]) 
names(multhyp)<- c("uid", "journal_id", "first_author", "year_pub", "mult_hyp", "corr", "data_avail", "code", "bayesian")

sum(multhyp$mult_hyp)/nrow(multhyp) #85% of studies use multiple hypothesis testing
sum(multhyp$corr)/sum(multhyp$mult_hyp) # of those studies that use multiple hypothesis testing, 14% use a correction


# setting up data for graphing
multhyp$mult_hyp_YN <- "Yes"
multhyp$mult_hyp_YN[multhyp$mult_hyp == 0] <- "No"
multhyp$correction <- "Corrected"
multhyp$correction[multhyp$corr == 0] <- "Not corrected"
multhyp$correction <- factor(multhyp$correction, levels = c("Not corrected", "Corrected"))
multhyp$data_avail[multhyp$data_avail == 0] <- "No"
multhyp$data_avail[multhyp$data_avail == 1] <- "Yes"
multhyp$code[multhyp$code == 0] <- "No"
multhyp$code[multhyp$code == 1] <- "Yes"
multhyp$bayesian[multhyp$bayesian == "Yes"] <- 1
multhyp$bayesian[multhyp$bayesian == "Yes"] <- 0

nrow(multhyp[which(multhyp$bayesian == 1 & multhyp$mult_hyp == 0),]) # bayesian no multiple hypothesis
nrow(multhyp[which(multhyp$bayesian == 1 & multhyp$mult_hyp == 1),]) # bayesian multiple hypothesis
nrow(multhyp[which(multhyp$bayesian == 1 & multhyp$mult_hyp == 1 & multhyp$corr == 1),]) # bayesian multiple hypothesis & correction

###################
#### Figure 4 ####
##################
ggsave(here("Figures", "Multiple_hypothesis_testing.png"),
       ggplot(aes(x = mult_hyp_YN), data = multhyp) + 
         geom_bar(aes(fill = correction,y = (..count..)/sum(..count..)), color = "black")+
         scale_fill_manual(values = c("white", "gray")) +
         scale_y_continuous(labels = scales::percent_format(accuracy = 1), n.breaks = 5) +
         labs(x = "Multiple Hypothesis Testing", fill= "Correction Used", y = "Percent of Studies") + 
         theme_pubclean() +theme(axis.title = element_text(face = "bold"), legend.title = element_blank()),
       width = 3,
       height = 3,
       dpi = 300
)

nrow(multhyp[multhyp$data_avail =="Yes",])/nrow(multhyp) # 78% of studies have data available
gr1 <- ggplot(aes(x = as.factor(data_avail)), data = multhyp) +
  geom_bar(aes(fill = journal_id, y = (..count..)/sum(..count..)), color = "black") + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), n.breaks = 6, limits = c(0,.9)) +
  scale_fill_brewer(palette = "Dark2") + 
  labs(x = "Data Available", fill= "Journal", y = "Percent of Studies") + 
  theme_pubclean() +  
  guides(fill = FALSE) +
  theme(axis.title = element_text(face = "bold"))

nrow(multhyp[multhyp$code =="Yes",])/nrow(multhyp) # 17.7% of studies have code available

gr2 <- ggplot(aes(x = as.factor(code)), data = multhyp) +
  geom_bar(aes(fill = journal_id, y = (..count..)/sum(..count..)), color = "black") + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), n.breaks = 6, limits = c(0,.9)) +
  scale_fill_brewer(palette = "Dark2") + 
  labs(x = "Code Available", fill= "Journal", y = "Percent of Studies") + 
  theme_pubclean() + 
  theme(axis.title = element_text(face = "bold"), legend.text = element_text(size = 10), 
        legend.title = element_text(size = 10), legend.key.size = unit(.5, units = "line"), legend.position = "right") 

##################
#### Figure 5 ####
###################
plots.mult <- (gr1 + gr2)
ggsave(here("Figures","data_code_avail.png"), plots.mult + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(face = "bold")), height = 6, width = 8, units = "in")


nrow(multhyp[multhyp$code == "Yes" & multhyp$data_avail == "Yes",])/nrow(multhyp) # 17.7% of studies with code and data

# graph of number of studies included by journal
ggsave(here("Figures", "Papers_SORTEE.png"),
ggplot(aes(x = include), data = papers[papers$include %in% c("No", "Yes"),]) + 
  geom_bar(aes(fill = journal_id,), color = "black",position = position_dodge()) + 
  #scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_brewer(palette = "Dark2") + 
  labs(x = "Included", fill= "Journal", y = "Number of Studies") + 
  theme_pubclean()+theme(axis.title = element_text(face = "bold")),
width = 6,
height = 6,
dpi = 300
)




################################
#### Why data were kicked ####
##############################

# pulling those removed because of unclear sample size 
notused <- kicked[kicked$Checking.Notes == "unsure of n" | 
                    kicked$Checking.Notes == "unsure of n (do not know number of points per time period)" |
                    kicked$Checking.Notes == "unsure of n " |
                    kicked$Checking.Notes == "unsure of n for shade and open plots" |
                    kicked$Checking.Notes == "can't find sample size" |
                    kicked$Checking.Notes == "unsure of n for Jaccard" |
                    kicked$Checking.Notes == "does not give number of cells used" |
                    kicked$Checking.Notes == "unsure of n - gives total n, not split by taxa" |
                    kicked$Checking.Notes == "unsure of n - does not split down observations by sex" |
                    kicked$Checking.Notes == "would need to go into data file and subset out lat and long with values associated with them" |
                    kicked$Checking.Notes == "unsure of observations per density" |
                    kicked$Checking.Notes == "only gives average sample size per treatment, not actual sample size" |
                    kicked$Checking.Notes == "unsure of N by date" |
                    kicked$Checking.Notes == "unsure ofn by species and date" |
                    kicked$Checking.Notes == "n unclear, currently assuming analyses at site level instead of plot level" |
                    kicked$Checking.Notes == "n unclear" |
                    kicked$Checking.Notes == "only gives average sample size per treatment" |
                    kicked$Checking.Notes == "unsure of n for mortality " |
                    kicked$Checking.Notes == "unsure of N - does not give number of saplings or total plots sampled",
]
notused <- notused[-which(is.na(notused$initials)),] # Lots of all NA rows - deleting
articlesnotused <- unique(notused[,c(2:6)])
nrow(notused) # 5411 estimates not used
nrow(articlesnotused) # 29 articles 


#########################
### Survey Results #####
#######################

# delete first 2 rows, not actual data
survey <- survey[-c(1,2),]
# take out responses that were not complete
survey <- survey[which(survey$Finished == "True"),]
# take out "do not agree"
survey <- survey[which(survey$Q1 == "I agree to participate"),]

# List of survey questions:
# Q2 = On average, what percentage of tests do you think passed the conventional target of 80% power?
# Q3 = Do you conduct ecological experiments?
# Q4 = If so, do you perform power analyses before starting a new experiment?
# Q5 = What is your career stage?

# creating dataframes for each question summarizing the number of answers in each category
Q2 <- as.data.frame(table(survey[,12]))
Q2$percentage <- Q2$Freq/238
Q3 <- as.data.frame(table(survey[,13]))
Q3$percentage <- Q3$Freq/238
Q4 <- as.data.frame(table(survey[,14]))
Q4 <- Q4[-which(Q4$Var1 == ""),]
Q4$percentage <- Q4$Freq/192
Q5 <- as.data.frame(table(survey[,15]))
Q5$percentage <- Q5$Freq/238
Q2.clicked <- as.data.frame(table(survey[,c(12,17)])) # did respondants click on the link before answering question 2

###############################
#### Supplemental figures ####
#############################

ggsave(here("Figures", "Survey_Q2_SORTEE.png"),
       ggplot(data = survey) + 
  geom_histogram(aes(y = (..count..)/sum(..count..), x = Q2), stat = "count", color = "black") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), n.breaks = 8) +  
  xlab("What percentage of tests do you think passed \nthe conventional target of 80% power?") +
  ylab("Percentage of respondents") +
  theme_pubclean() + 
  theme(axis.title = element_text(face = "bold"), legend.text = element_text(size = 7), 
        legend.title = element_text(size = 7), axis.text.x = element_text(angle = 90)),
  width = 6,
  height = 4,
  dpi = 300
)

q3 <- ggplot(data = survey) + 
  geom_histogram(aes(y = Q3), stat = "count", color = "black") +
  ylab("Do you conduct ecological experiments?") +
  theme_pubclean() + 
  theme(axis.title = element_text(face = "bold"), legend.text = element_text(size = 7), 
        legend.title = element_text(size = 7), legend.key.size = unit(.5, units = "line"))

survey$Q4 <- factor(survey$Q4, levels = c("Never", "Less than 25% of the time", "25 - 50% of the time",
                                          "50 - 75% of the time", "75% or more of the time", "Always", "" ))
ggsave(here("Figures", "Survey_Q4_SORTEE.png"),
ggplot(data = survey[-which(survey$Q4 == ""),]) + 
  geom_histogram(aes(x = (..count..)/sum(..count..), y = Q4), stat = "count", color = "black") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), n.breaks = 8, limits = c(0,.35)) +  
  ylab("Do you perform power anaylses \nbefore starting a new experiment?") +
  xlab("Percentage of respondents") +
  theme_pubclean() + 
  theme(axis.title = element_text(face = "bold"), legend.text = element_text(size = 7), 
        legend.title = element_text(size = 7), legend.key.size = unit(.5, units = "line")),
width = 6,
height = 6,
dpi = 300
)

survey$Q5 <- factor(survey$Q5, levels =c("Graduate Student", "Post-doc", "Faculty", "Researcher outside academia", "Other"))
q5 <- ggplot(data = survey) + 
  geom_histogram(aes(x = (..count..)/sum(..count..), y = Q5), stat = "count", color = "black") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), n.breaks = 8) +  
  ylab("Current position") +
  xlab("Percentage of respondents") +
  theme_pubclean() + 
  theme(axis.title = element_text(face = "bold"), legend.text = element_text(size = 7), 
        legend.title = element_text(size = 7), legend.key.size = unit(.5, units = "line"))


#### END OF CODE ####