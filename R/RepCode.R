#####################################
## Ecology Replicability Analyses ##
####################################


#### Load libraries
library(here)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(scales)
#library(devtools) only load if patchwork not installed
#install_github("thomasp85/patchwork")
library(patchwork)

#### set seed
set.seed(2617)

### load in  cleaned data 
# PRIMARY DATASET WITH ESTIMATES
df <- read.csv(here("Data", "CleanedDat.csv")) 
# ALL PAPERS EXPLORED
papers <- read.csv(here("Data", "CleanedPapers.csv"), row.names = 1)
# DATA THAT WAS REMOVED FROM THE ESTIMATES FILE AFTER CHECKING
kicked <- read.csv(here("Data", "Kicked.csv"), row.names = 1)
# DATA FROM SURVEY RESULTS NOTE: This data was requested to be removed from the 
# final manuscript by the Nature Ecology & Evolution Editor and a reviewer
# The analyses are still in this file, just commented out. Please uncomment them
# if you wish to see the output.
# survey <- read.csv(here("Data","RepEco_Survey.csv"))


##### GETTING DATA READY FOR ANALYSES
# create a unique for each paper in the dataset (uid)
uid_papers <- unique(df[,c(1:5)])
uid_papers$uid <- seq(1:nrow(uid_papers))
df <- merge(df,uid_papers) # merge uids with full dataset

# summary stats
median(df$sample_size) #79  
min(df$sample_size) #2
max(df$sample_size) #580978
mean(df$sample_size) #2999


# Get rid of estimates with SE of 0
nrow(df[which(df$std_error == 0),]) #810 estimates with SE of 0
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

## Get rid of entries with t-stats above the 99th percentile = estimates > 96.4,
# using 97 as the cutoff for simplicity (253 ESTIMATES TOTAL)
nrow(df[df$abs_tstat_sm > 97,])
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
median(df1$abs_pcc)
# threshold for 75% power
df1$WLS_threshold.75 <- weighted.mean(df1$abs_pcc, df1$precision_sq)/2.63
# threshold for 60% power
df1$WLS_threshold.60 <- weighted.mean(df1$abs_pcc, df1$precision_sq)/2.21

# create columns to say if estimate is adequately powered at that cutoff
df1$powered <- NA
df1$powered.75 <- NA
df1$powered.6 <- NA
# Conventional 80% power
for(i in 1:nrow(df1)){
  if(df1$SE_pcc[i] <= df1$WLS_threshold[i]){
    df1$powered[i] = 1
  } else{
    df1$powered[i] = 0
  }
}
sum(df1$powered)/nrow(df1) # proportion of estimates that make the 80% cutoff - 13.3%

# 75% power
for(i in 1:nrow(df1)){
  if(df1$SE_pcc[i] <= df1$WLS_threshold.75[i]){
    df1$powered.75[i] = 1
  } else{
    df1$powered.75[i] = 0
  }
}
# 60% power
for(i in 1:nrow(df1)){
  if(df1$SE_pcc[i] <= df1$WLS_threshold.6[i]){
    df1$powered.6[i] = 1
  } else{
    df1$powered.6[i] = 0
  }
}
sum(df1$powered.75)/nrow(df1)# proportion of estimates that make the 75% cutoff - 14%
sum(df1$powered.6)/nrow(df1)# proportion of estimates that make the 60% cutoff - 18%


####################
##### Figure 1a ####
###################
# Histogram of SE of PCC to show underpowered estimates
f1a <- ggplot(df1, aes(x = SE_pcc)) +
  geom_histogram(aes(y = (..count..)/sum(..count..)), bins = 50,color = "black", fill = "gray") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), n.breaks = 6) + 
  scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)) +
  geom_vline(xintercept = weighted.mean(df1$abs_pcc, df1$precision_sq)/2.8, color = "red", lty = "dashed") +
  annotate(geom="text", x=(weighted.mean(df1$abs_pcc, df1$precision_sq)/2.8)+.2, y=.115, label="Under-powered estimates", size = 4) +
  geom_segment(x =  weighted.mean(df1$abs_pcc, df1$precision_sq)/2.8, y = .115,
    xend =  weighted.mean(df1$abs_pcc, df1$precision_sq)/2.8 +.05, yend = .115,
    lineend = "round", linejoin = "round", size = 1, arrow = arrow(length = unit(0.07, "inches"))) +
  labs(x = "Standard Error of Partial Correlation Coefficients (PCC)", y = "Percentage of Estimates") + 
    theme_pubclean()+ 
  theme(axis.title = element_text(face = "bold"), text = element_text(size = 12))


#####################
#####Figure 1b######
####################
# Range of PCC Thresholds - moved to main text for R1 & 4 
# Create a dataframe 
newdf <- data.frame("pcc_values"= seq(0.01, 0.3, by = 0.01), "no_estimates" = NA)
newdf$threshold <- newdf$pcc_values/2.8
# count the number of estimates that would fall into the hypothetical threshold
for(i in 1:nrow(newdf)){
  newdf$no_estimates[i] <- sum(df1$SE_pcc >= newdf$threshold[i])
}

# calculate the frequency 
newdf$freq <- newdf$no_estimates/nrow(df1)
newdf$pcc_values <- as.factor(newdf$pcc_values)

f1b <- ggplot(aes(x = pcc_values, y = freq), data = newdf) + 
  geom_bar(stat = "identity", color = "black") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(y = "Underpowered \n(% of estimates)", x = "Partial Correlation Coefficient (PCC)") +
  geom_curve(x =10, y = .9,
               xend =  6, yend = .9,
               lineend = "round", size = 1, arrow = arrow()) +
  annotate(geom="text", x=14, y=.90, label="Weighted PCC value \ncalculated in this study", size = 4) + 
  theme_pubclean() + 
  theme(axis.title = element_text(face = "bold"), text = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust = 1))

fig1 <- f1a/f1b
ggsave(here("Figures","figure1.pdf"), fig1 + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(face = "bold")), height = 150, width = 180, units = "mm", dpi = 300)


################################################
#### Power analyses for main estimates only ####
###############################################
# data frame with main results (N.B. results classified as main not only those presented in the main text)
df2 <- df1[df1$main_results..Y.1..N.0. == 1,]
# calculate threshold and value
df2$WLS_threshold <- weighted.mean(df2$abs_pcc, df2$precision_sq)/2.8
df2$WLS_value <- weighted.mean(df2$abs_pcc, df2$precision_sq)
# count estimates that pass threshould
df2$powered <- NA
for(i in 1:nrow(df2)){
  if(df2$SE_pcc[i] <= df2$WLS_threshold[i]){
    df2$powered[i] = 1
  } else{
    df2$powered[i] = 0
  }
}

sum(df2$powered)/nrow(df2) # 11% reach 80% threshold

#############################################################
#### Power analyses for estimates presented in main text ####
############################################################
## tables in main text
df3 <- df1[df1$table_loc == "main",]

df3$WLS_threshold <- weighted.mean(df3$abs_pcc, df3$precision_sq)/2.8
df3$WLS_value <- weighted.mean(df3$abs_pcc, df3$precision_sq)
df3$powered <- NA

for(i in 1:nrow(df3)){
  if(df3$SE_pcc[i] <= df3$WLS_threshold[i]){
    df3$powered[i] = 1
  } else{
    df3$powered[i] = 0
  }
}
sum(df3$powered)/nrow(df3) # 18% meet 80% threshold

#############################################################
#### Power analyses for estimates presented in supp text ####
############################################################
# supplemental tables
df4 <- df1[df1$table_loc == "supplement",]

df4$WLS_threshold <- weighted.mean(df4$abs_pcc, df4$precision_sq)/2.8
df4$WLS_value <- weighted.mean(df4$abs_pcc, df4$precision_sq)
df4$powered <- NA

for(i in 1:nrow(df4)){
  if(df4$SE_pcc[i] <= df4$WLS_threshold[i]){
    df4$powered[i] = 1
  } else{
    df4$powered[i] = 0
  }
}

sum(df4$powered)/nrow(df4) #13% meet threshold

#############################################################
#### Power analyses for observational study estimates   ####
############################################################
df5 <- df1[df1$study_type == "observational",]

df5$WLS_threshold <- weighted.mean(df5$abs_pcc, df5$precision_sq)/2.8
df5$WLS_value <- weighted.mean(df5$abs_pcc, df5$precision_sq)
df5$powered <- NA

for(i in 1:nrow(df5)){
  if(df5$SE_pcc[i] <= df5$WLS_threshold[i]){
    df5$powered[i] = 1
  } else{
    df5$powered[i] = 0
  }
}
sum(df5$powered)/nrow(df5) # 12% meet threshold
#############################################################
#### Power analyses for experimental study estimates    ####
############################################################
df6 <- df1[df1$study_type == "experimental",]

df6$WLS_threshold <- weighted.mean(df6$abs_pcc, df6$precision_sq)/2.8
df6$WLS_value <- weighted.mean(df6$abs_pcc, df6$precision_sq)
df6$powered <- NA

for(i in 1:nrow(df6)){
  if(df6$SE_pcc[i] <= df6$WLS_threshold[i]){
    df6$powered[i] = 1
  } else{
    df6$powered[i] = 0
  }
}
sum(df6$powered)/nrow(df6) # 43% meet threshold

######################
#### Median power ####
#####################
#"Median power for a given area of research can then be calculated from the cumulative normal 
#probability of the difference between 1.96 and the absolute value of an estimate of the true 
#effect divided by the median standard error"

# full 
WLS_threshold_median = unique(df1$WLS_value)/median(df1$SE_pcc)
1-pnorm(1.96-WLS_threshold_median) # 13.4

# main
WLS_threshold_median = unique(df2$WLS_value)/median(df2$SE_pcc)
1-pnorm(1.96-WLS_threshold_median) #10.9

# in main text
WLS_threshold_median = unique(df3$WLS_value)/median(df3$SE_pcc)
1-pnorm(1.96-WLS_threshold_median) # 16.3

# in supplemental text
WLS_threshold_median = unique(df4$WLS_value)/median(df4$SE_pcc)
1-pnorm(1.96-WLS_threshold_median) #13.4

# in observational studies
WLS_threshold_median = unique(df5$WLS_value)/median(df5$SE_pcc)
1-pnorm(1.96-WLS_threshold_median) # 10.9

# in experimental studies
WLS_threshold_median = unique(df6$WLS_value)/median(df6$SE_pcc)
1-pnorm(1.96-WLS_threshold_median) #68

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

########################################
##### Figure 2a - Exaggeration Bias #####
#######################################
f2a <- ggplot(aes(x = category, y = percentage_of_studies), data = category_counts1) + 
  geom_bar(stat = "identity", color = "black", fill = "gray") + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(y = "Percentage of Estimates", x = "Exaggeration Bias") +
  theme_pubclean() +  theme(text = element_text(size = 12), axis.title = element_text(face = "bold"))


###########################################
#### Figure 2b - Range of WAAP values ####
##########################################
waapdf<- data.frame("waap_values"= seq(0.01, 0.3, by = 0.01), "no_estimates" = NA)
# count the number of estimates that would fall into the hypothetical threshold
for(i in 1:nrow(waapdf)){
  waapdf$no_estimates[i] <- sum(df1.u$abs_pcc/waapdf$waap_values[i] >= 2)
}

# calculate the frequency 
waapdf$freq <- waapdf$no_estimates/nrow(df1.u)
waapdf$waap_values <- as.factor(waapdf$waap_values)


f2b <- ggplot(aes(x = waap_values, y = freq), data = waapdf) + 
         geom_bar(stat = "identity", color = "black") +
         scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
         labs(y = "Exaggerated by 100% or more \n(% of estimates)", x = "Weighted average of the adequately powered estimator (WAAP)") +
         geom_curve(x =10, y = .8,
                    xend =  5, yend = .67,
                    lineend = "round", size = 1, arrow = arrow()) +
         annotate(geom="text", x=11, y=.72, label="approx. WAAP value \ncalculated in this study", size = 4) + 
         theme_pubclean() + 
         theme(axis.title = element_text(face = "bold"), text = element_text(size = 12), axis.text = element_text(angle = 45, hjust = 1))


fig2 <- f2a/f2b
ggsave(here("Figures","figure2.pdf"), fig2 + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(face = "bold")), height = 150, width = 180, units = "mm", dpi = 300)

###########################
#### p-hacking curves ####
##########################
# curves are plotted several different ways
# using 1000 bins and an epanechnikov kernal for plotting
# Focusing on t-stats <10 to exclude long tail of distribution

# 1. unweighted
ggplot(data = df1[df1$abs_tstat_sm<10,]) +
  geom_histogram(aes(x = abs_tstat_sm,  y = ..density..), bins = 1000, fill = "gray") + 
  geom_density(aes(x = abs_tstat_sm), color = "black", kernel = "epanechnikov", bw = 0.2) + 
  scale_x_continuous(breaks = c(0,1.96,10)) +
  theme_pubclean() +
  labs(title = "Rounded T-stats", x = "|t-stat|")

# 2. weighted by number of estimates per articles
ggplot(data = df1[df1$abs_tstat_sm <10,]) +
  geom_histogram(aes(x = abs_tstat_sm,  y = ..density.., weight = weight_article), 
                 bins = 1000, fill = "gray") + 
  geom_density(aes(x = abs_tstat_sm, weight = weight_article), color = "black", kernel = "epanechnikov", bw = 0.2) + 
  scale_x_continuous(breaks = c(0,1.96,10)) +
  geom_segment(x =  1.94, y = .4,
               xend =  1.94, yend = .2,
               lineend = "round", linejoin = "round", size = .5, arrow = arrow(length = unit(0.07, "inches")), color = "red") +
  theme_pubclean() + 
  labs(title ="Weighted by Article", x= "|t-stat|")

# 3. Weighted by number of estimates per table per article
  ggplot(data = df1[df1$abs_tstat_sm <10,]) +
  geom_histogram(aes(x = abs_tstat_sm,  y = ..density.., weight = weight_table), 
                 bins = 1000, fill = "gray") + 
  geom_density(aes(x = abs_tstat_sm, weight = weight_table), color = "black", kernel = "epanechnikov", bw = 0.2) + 
  scale_x_continuous(breaks = c(0,1.96,10)) +
  geom_segment(x =  1.9, y = .4, xend =  1.9, yend = .25,
               lineend = "round", linejoin = "round", size = .5, arrow = arrow(length = unit(0.07, "inches")), color = "red") +
  theme_pubclean() + theme(text = element_text(size = 14), axis.title = element_text(face = "bold")) + 
  labs(x = "t-statistic", y = "Density")

# 4. Weighted by number of estimates per table per article for "main" estimates
  ggplot(data = df2[df2$abs_tstat_sm <10,]) +
    geom_histogram(aes(x = abs_tstat_sm,  y = ..density.., weight = weight_table), 
                   bins = 1000, fill = "gray") + 
    geom_density(aes(x = abs_tstat_sm, weight = weight_table), color = "black", kernel = "epanechnikov", bw = 0.2) + 
    scale_x_continuous(breaks = c(0,1.96,10)) +
    geom_segment(x =  1.9, y = .4, xend =  1.9, yend = .25,
                 lineend = "round", linejoin = "round", size = .5, arrow = arrow(length = unit(0.07, "inches")), color = "red") +
    theme_pubclean() + theme(text = element_text(size = 14), axis.title = element_text(face = "bold")) + 
    labs(x = "t-statistic", y = "Density")
  
  
  
##########################
#### Figure 3 graphs ####
#########################

# Figure 3c
# full sample - weighted by estimates per table per article
gra1 <- ggplot(data = df1[df1$abs_tstat_sm <10,]) +
  geom_histogram(aes(x = abs_tstat_sm,  y = ..density.., weight = weight_table), 
                 bins = 50, fill = "gray") + 
  geom_density(aes(x = abs_tstat_sm, weight = weight_table), color = "black", kernel = "epanechnikov", bw = 0.2) + 
  scale_x_continuous(breaks = c(0,1.96,10)) +
  ylim(c(0,.35)) + 
  geom_segment(x =  1.9, y = .3, xend =  1.9, yend = .25,
               lineend = "round", linejoin = "round", size = .5, arrow = arrow(length = unit(0.07, "inches")), color = "red") +
  theme_pubclean() + theme(text = element_text(size = 12), axis.title = element_text(face = "bold"), plot.title = element_text(hjust = 0.5)) + 
  labs(x = "t-statistic", y = "Density", title= "All Tables")

  nrow(df1[df1$abs_tstat_sm <10,]) #16,950
# Figure 3a
# estimates presented in main text - weighted by estimates per table per article
gra2 <- ggplot(data = df3[df3$abs_tstat_sm <10,]) +
  geom_histogram(aes(x = abs_tstat_sm,  y = ..density.., weight = weight_table), 
                 bins = 50, fill = "gray") + 
  geom_density(aes(x = abs_tstat_sm, weight = weight_table), color = "black", kernel = "epanechnikov", bw = 0.2) + 
  scale_x_continuous(breaks = c(0, 1.96,10)) +
  ylim(c(0,.35)) + 
  geom_segment(x =  1.9, y = .3, xend =  1.9, yend = .25,
               lineend = "round", linejoin = "round", size = .5, arrow = arrow(length = unit(0.07, "inches")), color = "red") +
  theme_pubclean() + theme(text = element_text(size = 12), axis.title = element_text(face = "bold"), plot.title = element_text(hjust = 0.5)) + 
  labs(x = "t-statistic", y = "Density", title = "Main Text Tables")
nrow(df3[df3$abs_tstat_sm <10,]) # 2,278
# Figure 3b
# estimates in supplemental text - weighted by estimates per table per article
gra3 <- ggplot(data = df1[df1$abs_tstat_sm <10 & df1$table_loc == "supplement",]) +
  geom_histogram(aes(x = abs_tstat_sm,  y = ..density.., weight = weight_table), 
                 bins = 50, fill = "gray") + 
  geom_density(aes(x = abs_tstat_sm, weight = weight_table), color = "black", kernel = "epanechnikov") + 
  scale_x_continuous(breaks = c(0,1.96,10)) +
  ylim(c(0,.35)) + 
  geom_segment(x =  1.9, y = .3, xend =  1.9, yend = .25,
               lineend = "round", linejoin = "round", size = .5, arrow = arrow(length = unit(0.07, "inches")), color = "red") +
  theme_pubclean() + theme(text = element_text(size = 12), axis.title = element_text(face = "bold"), plot.title = element_text(hjust = 0.5)) + 
  labs(x = "t-statistic", y = "Density", title = "Supplement Tables")
nrow(df1[df1$abs_tstat_sm <10 & df1$table_loc == "supplement",]) #14,672

plots <- gra2/gra3/gra1
ggsave(here("Figures","figure3.pdf"), plots + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(face = "bold")), height = 210, width = 150, units = "mm", dpi = 300)


#####################################################################
### p-hacking curves for experimental and observational studies ####
### not presented in manuscript - looking at for Reviewer ##########
####################################################################

ggplot(data = df5[df5$abs_tstat_sm <10 & df5$table_loc == "main",]) +
  geom_histogram(aes(x = abs_tstat_sm,  y = ..density.., weight = weight_table), 
                 bins = 50, fill = "gray") + 
  geom_density(aes(x = abs_tstat_sm, weight = weight_table), color = "black", kernel = "epanechnikov") + 
  scale_x_continuous(breaks = c(0,1.96,10)) +
  ylim(c(0,.35)) + 
  geom_segment(x =  1.9, y = .3, xend =  1.9, yend = .25,
               lineend = "round", linejoin = "round", size = .5, arrow = arrow(length = unit(0.07, "inches")), color = "red") +
  theme_pubclean() + theme(text = element_text(size = 12), axis.title = element_text(face = "bold"), plot.title = element_text(hjust = 0.5)) + 
  labs(x = "t-statistic", y = "Density", title = "Observational Tables")

ggplot(data = df6[df6$abs_tstat_sm <10 & df6$table_loc == 'main',]) +
  geom_histogram(aes(x = abs_tstat_sm,  y = ..density.., weight = weight_table), 
                 bins = 50, fill = "gray") + 
  geom_density(aes(x = abs_tstat_sm, weight = weight_table), color = "black", kernel = "epanechnikov") + 
  scale_x_continuous(breaks = c(0,1.96,10)) +
  ylim(c(0,.35)) + 
  geom_segment(x =  1.9, y = .3, xend =  1.9, yend = .25,
               lineend = "round", linejoin = "round", size = .5, arrow = arrow(length = unit(0.07, "inches")), color = "red") +
  theme_pubclean() + theme(text = element_text(size = 12), axis.title = element_text(face = "bold"), plot.title = element_text(hjust = 0.5)) + 
  labs(x = "t-statistic", y = "Density", title = "Experimental")
nrow(df1[df1$abs_tstat_sm <10 & df1$study_type == "experimental" & df1$table_loc == "main",]) #1069
nrow(df1[df1$abs_tstat_sm <10 & df1$study_type == "observational" & df1$table_loc == "supplement",]) #9393

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

nrow(multhyp[which(multhyp$bayesian == 1 & multhyp$mult_hyp == 0),]) # bayesian no multiple hypothesis = 9
nrow(multhyp[which(multhyp$bayesian == 1 & multhyp$mult_hyp == 1),]) # bayesian multiple hypothesis = 21
nrow(multhyp[which(multhyp$bayesian == 1 & multhyp$mult_hyp == 1 & multhyp$corr == 1),]) # bayesian multiple hypothesis & correction = 1

###################
#### Figure 4 ####
##################

# bar graph with multiple hypothesis testing, yes/no 

ggsave(here("Figures", "figure4.pdf"),
       ggplot(aes(x = mult_hyp_YN), data = multhyp) + 
         geom_bar(aes(fill = correction,y = (..count..)/sum(..count..)), color = "black")+
         scale_fill_manual(values = c("white", "gray")) +
         scale_y_continuous(labels = scales::percent_format(accuracy = 1), n.breaks = 5) +
         labs(x = "Multiple Hypothesis Testing", fill= "Correction Used", y = "Percent of Studies") + 
         theme_pubclean() +theme(axis.title = element_text(face = "bold"), legend.title = element_blank()),
       width = 88,
       height =100,
       units = 'mm',
       dpi = 300
)


##################
#### Figure 5 ####
##################
# percentage of studies with data available Figure 5a
nrow(multhyp[multhyp$data_avail =="Yes",])/nrow(multhyp) # 78% of studies have data available
gr1 <- ggplot(aes(x = as.factor(data_avail)), data = multhyp) +
  geom_bar(aes(fill = journal_id, y = (..count..)/sum(..count..)), color = "black") + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), n.breaks = 6, limits = c(0,.9)) +
  scale_fill_brewer(palette = "Dark2") + 
  labs(x = "Data Available", fill= "Journal", y = "Percent of Studies") + 
  theme_pubclean() +  
  guides(fill = FALSE) +
  theme(axis.title = element_text(face = "bold"))


# percentage of studies with code available Figure 5b
nrow(multhyp[multhyp$code =="Yes",])/nrow(multhyp) # 18% of studies have code available
gr2 <- ggplot(aes(x = as.factor(code)), data = multhyp) +
  geom_bar(aes(fill = journal_id, y = (..count..)/sum(..count..)), color = "black") + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), n.breaks = 6, limits = c(0,.9)) +
  scale_fill_brewer(palette = "Dark2") + 
  labs(x = "Code Available", fill= "Journal", y = "Percent of Studies") + 
  theme_pubclean() + 
  theme(axis.title = element_text(face = "bold"), legend.text = element_text(size = 10), 
        legend.title = element_text(size = 10), legend.key.size = unit(.5, units = "line"), legend.position = "right") 


plots.mult <- (gr1 + gr2)
ggsave(here("Figures","figure5.pdf"), plots.mult + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(face = "bold")), height = 88, width = 130, units = "mm", dpi = 300)


nrow(multhyp[multhyp$code == "Yes" & multhyp$data_avail == "Yes",])/nrow(multhyp) # 17.7% of studies with code and data

# graph of number of studies included by journal
# not included in manuscript
ggplot(aes(x = include), data = papers[papers$include %in% c("No", "Yes"),]) + 
  geom_bar(aes(fill = journal_id,), color = "black",position = position_dodge()) + 
  #scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_brewer(palette = "Dark2") + 
  labs(x = "Included", fill= "Journal", y = "Number of Studies") + 
  theme_pubclean()


################################
#### Why data were kicked ####
##############################
# First look at papers that never had values recorded.
removed <- papers[which(papers$include == "No"),]
nrow(removed[removed$reason== "No Tables",]) #1024 
nrow(removed[removed$reason == "Paper Type",]) #132
nrow(removed[removed$reason == "No Error Reported",]) #14 + 1


# get papers in kicked which match uid of papers in df1 - these are included in the final sample
# but may have had some estimates removed for some reson or the other

df1_papers <- unique(df1[,c(1,3:7)])
kicked_sub <- kicked[-which(is.na(kicked$reason)),]
kicked_sub <- unique(kicked[,c(2:6,26)])
kicked_sub <- merge(kicked_sub, df1_papers, all.x = TRUE)

# pulling those removed because of unclear sample size 
notused <- kicked[which(kicked$reason == "unclear n"),]
articlesnotused <- unique(notused[,c(2:6)])
nrow(notused) # 5412 estimates not used
nrow(articlesnotused) # 29 articles - includes some that were in the final dataset.
articles <- merge(articlesnotused, df1_papers, all.x = TRUE)
articles <- articles[which(is.na(articles$uid)),]
articles$reason <- "unsure N"
# Now looking at other reasons
kicked_sub <- kicked_sub[which(is.na(kicked_sub$uid)),] #any papers without an UID were not in the final dataset
unique(kicked_sub$reason)
nrow(kicked_sub[which(kicked_sub$reason == "Not model coef"),]) # 13
nrow(kicked_sub[which(kicked_sub$reason == "meta-analysis"),]) # 4
nrow(kicked_sub[which(kicked_sub$reason == "unclear n"),])
x <- merge(articles, kicked_sub, all = TRUE)

#########################
### Survey Results #####
#######################
# 
# # delete first 2 rows, not actual data
# survey <- survey[-c(1,2),]
# # take out responses that were not complete
# survey <- survey[which(survey$Finished == "True"),]
# # take out "do not agree"
# survey <- survey[which(survey$Q1 == "I agree to participate"),]
# 
# # List of survey questions:
# # Q2 = On average, what percentage of tests do you think passed the conventional target of 80% power?
# # Q3 = Do you conduct ecological experiments?
# # Q4 = If so, do you perform power analyses before starting a new experiment?
# # Q5 = What is your career stage?
# 
# # creating dataframes for each question summarizing the number of answers in each category
# Q2 <- as.data.frame(table(survey[,12]))
# Q2$percentage <- Q2$Freq/238
# Q3 <- as.data.frame(table(survey[,13]))
# Q3$percentage <- Q3$Freq/238
# Q4 <- as.data.frame(table(survey[,14]))
# Q4 <- Q4[-which(Q4$Var1 == ""),]
# Q4$percentage <- Q4$Freq/192
# Q5 <- as.data.frame(table(survey[,15]))
# Q5$percentage <- Q5$Freq/238
# Q2.clicked <- as.data.frame(table(survey[,c(12,17)])) # did respondents click on the link before answering question 2
# 
# #############################################
# #### Supplemental figures Survey Results ####
# #############################################
# q2 <- ggplot(data = survey) + 
#   geom_histogram(aes(x = (..count..)/sum(..count..), y = Q2), stat = "count", color = "black") +
#   scale_x_continuous(labels = scales::percent_format(accuracy = 1), n.breaks = 8) +  
#   ylab("What percentage of tests do you think passed \nthe conventional target of 80% power?") +
#   xlab("Percentage of repsondents") +
#   theme_pubclean() + 
#   theme(axis.title = element_text(face = "bold"), legend.text = element_text(size = 7), 
#         legend.title = element_text(size = 7), legend.key.size = unit(.5, units = "line"))
# 
# q3 <- ggplot(data = survey) + 
#   geom_histogram(aes(y = Q3), stat = "count", color = "black") +
#   ylab("Do you conduct ecological experiments?") +
#   theme_pubclean() + 
#   theme(axis.title = element_text(face = "bold"), legend.text = element_text(size = 7), 
#         legend.title = element_text(size = 7), legend.key.size = unit(.5, units = "line"))
# 
# survey$Q4 <- factor(survey$Q4, levels = c("Never", "Less than 25% of the time", "25 - 50% of the time",
#                                           "50 - 75% of the time", "75% or more of the time", "Always", "" ))
# q4 <- ggplot(data = survey[-which(survey$Q4 == ""),]) + 
#   geom_histogram(aes(x = (..count..)/sum(..count..), y = Q4), stat = "count", color = "black") +
#   scale_x_continuous(labels = scales::percent_format(accuracy = 1), n.breaks = 8, limits = c(0,.35)) +  
#   ylab("Do you perform power anaylses \nbefore starting a new experiment?") +
#   xlab("Percentage of respondents") +
#   theme_pubclean() + 
#   theme(axis.title = element_text(face = "bold"), legend.text = element_text(size = 7), 
#         legend.title = element_text(size = 7), legend.key.size = unit(.5, units = "line"),
#         axis.text.y = element_text(angle = 45))
# survey$Q5 <- factor(survey$Q5, levels =c("Graduate Student", "Post-doc", "Faculty", "Researcher outside academia", "Other"))
# q5 <- ggplot(data = survey) + 
#   geom_histogram(aes(x = (..count..)/sum(..count..), y = Q5), stat = "count", color = "black") +
#   scale_x_continuous(labels = scales::percent_format(accuracy = 1), n.breaks = 8) +  
#   ylab("Current position") +
#   xlab("Percentage of respondents") +
#   theme_pubclean() + 
#   theme(axis.title = element_text(face = "bold"), legend.text = element_text(size = 7), 
#         legend.title = element_text(size = 7), legend.key.size = unit(.5, units = "line"))
# 
# # saving figures for supplement
# png(here("Figures/surv_results.png"), width = 9, height = 6, units = "in", res = 300)
# ggarrange(plotlist = list(q2,q4), ncol = 2, legend = "right", labels = c("A", "B"), widths = c(1,1.25))
# dev.off()


###############################################
## Supplemental Figures - PCC range for R4 ####
###############################################

sup1 <- ggplot(aes(x = abs_pcc), data = df1) + 
  geom_histogram(aes(y = (..count..)/sum(..count..), fill = study_type), bins = 20, color = "black") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(breaks = c(0, 0.1, .2, .3, .4, .5, .6, .7, .8, .9, 1))+
  scale_fill_manual(values = c("black", "gray", "forestgreen")) + 
  geom_vline(xintercept = 0.14, lty = "dashed") + 
  annotate(geom="text", x=median(df1$abs_pcc) + .2, y=.17, label="un-weighted median PCC", size = 4) +
  geom_segment(x =  0.2, y = .17,
               xend =  .14, yend = .17,
               lineend = "round", linejoin = "round", size = 1, arrow = arrow(length = unit(0.07, "inches"))) + 
  xlab("Absolute Value of Partial Correlation Coefficient (PCC)")+
  ylab("Percentage of Estimates") + 
  theme_pubclean()

sup2<- ggplot(aes(x = SE_pcc), data = df1) + 
  geom_histogram(aes(y = (..count..)/sum(..count..), fill = study_type), bins = 20, color = "black") + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7))+
  scale_fill_manual(values = c("black","gray", "forestgreen")) +
  xlab("Standard Error of Partial Correlation Coefficient (PCC)")+
  ylab("Percantage of Estimates") + 
  guides(fill = FALSE) + 
  theme_pubclean()

sup <- sup1/sup2
ggsave(here("Figures/supplemental_figure.png"), sup + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(face = "bold")), dpi = 300)


#####################################################
## For R3 - Number of mixed effects estimates ####
###################################################

nrow(df[df$mixed_effects ==1,])/nrow(df)  ## 42% of estimates are from mixed effect models

#####################################################
### For R4 - Number of studies that report power ####
#####################################################

power <- unique(df[,c(1,20)])
power <- power[-which(is.na(power$report_power)),]
nrow(power[power$report_power ==1,]) # one study reports power & it is only post hoc






#### END OF CODE ####
