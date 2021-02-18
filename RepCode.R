## Ecology Replicability Analyses ##
## Authors: Kaitlin Kimmel, Meghan Avolio, and Paul Ferraro
## Code adapted from Pallavi Shukla's Stata code

#### load libraries ####
library(here)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)

#### set seed ####
set.seed(2617)

### load in  cleaned data ###
df <- read.csv(here("Data", "CleanedDat.csv"), row.names = 1)
papers <- read.csv(here("Data", "CleanedPapers.csv"), row.names = 1)
kicked <- read.csv(here("Data", "Kicked.csv"), row.names = 1)
survey <- read.csv(here("Data","RepEco_Survey.csv"))

# create a unique for each paper in the dataset (uid)
uid_papers <- unique(df[,c(1:5)])
uid_papers$uid <- seq(1:nrow(uid_papers))
df <- merge(df,uid_papers) # merge uids with full dataset

# summary stats
median(df$sample_size[-which(is.na(df$sample_size))]) #79
min(df$sample_size[-which(is.na(df$sample_size))]) #2
max(df$sample_size[-which(is.na(df$sample_size))]) #580978
mean(df$sample_size[-which(is.na(df$sample_size))]) #2963.5

## Get rid of estimates with SE of 0

nrow(df[df$std_error == 0,]) #881
df <- df[df$std_error != 0,]
df <- df[-which(is.na(df$sample_size)),]
df <- df[-which(is.na(df$coefficient)),]

################################################
#### De-round estimates and create weights ####
###############################################
# Derounding as in Brodeur et al 2016
## function to find number of decimal places for derounding
## from stackoverflow user: daroczig
decimalplaces <- function(x) {
  if (abs(x - round(x)) > .Machine$double.eps^0.5) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed = TRUE)[[1]])
  } else {
    return(0)
  }
} 

# derounding recorded coefficients and standard errors
for(i in 1:nrow(df)){
  y <- decimalplaces(df$coefficient[i])
  df$coefficient_sm[i] <- runif(1,df$coefficient[i] - (0.5* 10^(-y)), df$coefficient[i] + (0.5* 10^(-y)))
  z <- decimalplaces(df$std_error[i])
  df$std_error_sm[i] <- runif(1,df$std_error[i] - (0.5* 10^(-y)), df$std_error[i] + (0.5* 10^(-y)))
}


# calculate t-stats
df$tstat_sm <- df$coefficient_sm/df$std_error_sm
df$abs_tstat_sm <- abs(df$tstat_sm)

quantile(df$abs_tstat_sm, c(.01,.95,.99))

## Get rid of entries with t-stats above the 99th percentile (256 ESTIMATES TOTAL)
df <- df[df$abs_tstat_sm < 97,]

# generate weights
obs_by_article <- plyr::count(df, vars = "uid")
names(obs_by_article)[2] <- "obs_by_article"
obs_by_table <- plyr::count(df, vars = c("uid", "table_no"))
names(obs_by_table)[3] <- "obs_by_table"

# tables in each article
temp_df <- df[,c("uid", "table_no")]
temp_df <- unique(temp_df)
temp_df <- temp_df %>% group_by(uid) %>% summarize(count=n())
names(temp_df)[2] <- "tab_count"

df<- merge(df, obs_by_article)
df <- merge(df, obs_by_table)
df <- merge(df, temp_df)

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
# subset data for pcc values that exist with full data set
df1 <- df[df$obs_by_article <6000,]
# calculate threshold and values
df1$WLS_threshold <- weighted.mean(df1$abs_pcc, df1$precision_sq)/2.8 #2.8 for .8, 2.63 for 0.75 and 2.21 for 0.6
df1$WLS_value <- weighted.mean(df1$abs_pcc, df1$precision_sq)
df1$WLS_threshold.75 <- weighted.mean(df1$abs_pcc, df1$precision_sq)/2.63
df1$WLS_threshold.60 <- weighted.mean(df1$abs_pcc, df1$precision_sq)/2.21
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

sum(df1$powered)/nrow(df1)
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
sum(df1$powered.75)/nrow(df1)
sum(df1$powered.6)/nrow(df1)

ggsave(
  here("Figures", "Power_full.png"),
ggplot(df1, aes(x = SE_pcc)) +
  geom_histogram(aes(y = (..count..)/sum(..count..)), bins = 50,color = "black", fill = "gray") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
  geom_vline(xintercept = weighted.mean(df1$abs_pcc, df1$precision_sq)/2.8, color = "red", lty = "dashed") +
  annotate(geom="text", x=(weighted.mean(df1$abs_pcc, df1$precision_sq)/2.8)+.2, y=.21, label="Under-powered estimates", size = 3) +
  geom_segment(x =  weighted.mean(df1$abs_pcc, df1$precision_sq)/2.8, y = .21,
    xend =  weighted.mean(df1$abs_pcc, df1$precision_sq)/2.8 +.05, yend = .21,
    lineend = "round", linejoin = "round", size = 1, arrow = arrow(length = unit(0.07, "inches"))) +
  labs(x = "Standard Error of PCC", y = "Percentage of Estimates") + 
    theme_pubclean()+ 
  theme(axis.title = element_text(face = "bold"), text = element_text(size = 14)), 
width = 5,
height = 3,
dpi = 300
)

## main estimates

df2 <- df1[df1$main_results..Y.1..N.0. == 1,]

df2$WLS_threshold <- weighted.mean(df2$abs_pcc, df2$precision_sq)/2.8
df2$WLS_value <- weighted.mean(df2$abs_pcc, df2$precision_sq)
df2$powered <- NA

for(i in 1:nrow(df2)){
  if(df2$SE_pcc[i] <= df2$WLS_threshold[i]){
    df2$powered[i] = 1
  } else{
    df2$powered[i] = 0
  }
}

sum(df2$powered)/nrow(df2)


ggsave(
  here("Figures", "Power_main.png"),
ggplot(df2, aes(x = SE_pcc)) +
  geom_histogram(aes(y = (..count..)/sum(..count..)), bins = 50, color = "black", fill = "gray") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
  geom_vline(xintercept = weighted.mean(df2$abs_pcc, df2$precision_sq)/2.8, color = "red", lty = "dashed") +
  annotate(geom="text", x=(weighted.mean(df2$abs_pcc, df2$precision_sq)/2.8)+.2, y=.2, label="Under-powered estimates", size = 3) + 
  geom_segment(x =  weighted.mean(df2$abs_pcc, df2$precision_sq)/2.8, y = .2,
               xend =  weighted.mean(df2$abs_pcc, df2$precision_sq)/2.8 +.04, yend = .2,
               lineend = "round", linejoin = "round", size = 1, arrow = arrow(length = unit(0.07, "inches"))) +
  labs(x = "Standard Error of PCC", y = "Percentage of Studies") + 
  theme_pubclean() + 
  theme(axis.title = element_text(face = "bold"), text = element_text(size = 14)),
width = 5,
height = 3,
dpi = 300
)
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

sum(df3$powered)/nrow(df3)

ggplot(df3, aes(x = SE_pcc)) +
  geom_histogram(aes(y = (..count..)/sum(..count..)), bins = 30, color = "black", fill = "gray") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
  geom_vline(xintercept = weighted.mean(df3$abs_pcc, df3$precision_sq)/2.8, color = "red", lty = "dashed") +
  geom_vline(xintercept = weighted.mean(df1$abs_pcc, df1$precision_sq)/2.8, color = "pink", lty = "dashed") +
  annotate(geom="text", x=(weighted.mean(df3$abs_pcc, df3$precision_sq)/2.8)+.17, y=.21, label="Under-powered estimates") + 
  geom_segment(x =  weighted.mean(df3$abs_pcc, df3$precision_sq)/2.8, y = .21,
               xend =  weighted.mean(df3$abs_pcc, df3$precision_sq)/2.8 +.04, yend = .21,
               lineend = "round", linejoin = "round", size = .5, arrow = arrow(length = unit(0.07, "inches"))) +
  labs(x = "Standard Error of PCC", y = "Percentage of Studies", title = "Experimental") + 
  theme_pubclean()

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

sum(df4$powered)/nrow(df4)

ggplot(df4, aes(x = SE_pcc)) +
  geom_histogram(aes(y = (..count..)/sum(..count..)), bins = 50, color = "black", fill = "gray") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
  geom_vline(xintercept = weighted.mean(df4$abs_pcc, df4$precision_sq)/2.8, color = "red", lty = "dashed") +
  geom_vline(xintercept = weighted.mean(df1$abs_pcc, df1$precision_sq)/2.8, color = "pink", lty = "dashed") +
  annotate(geom="text", x=(weighted.mean(df4$abs_pcc, df4$precision_sq)/2.8)+.2, y=.35, label="Under-powered estimates") + 
  geom_segment(x =  weighted.mean(df4$abs_pcc, df4$precision_sq)/2.8, y = .35,
               xend =  weighted.mean(df4$abs_pcc, df4$precision_sq)/2.8 +.04, yend = .35,
               lineend = "round", linejoin = "round", size = .5, arrow = arrow(length = unit(0.07, "inches"))) +
  labs(x = "Standard Error of PCC", y = "Percentage of Studies", title = "Observational") + 
  theme_pubclean()

####################
## Median power ###
###################
#"Median power for a given area of research can then be calculated from the cumulative normal 
#probability of the difference between 1.96 and the absolute value of an estimate of the true 
#effect divided by the median standard error

# full 
WLS_threshold_median = unique(df1$WLS_value)/median(df1$SE_pcc)
1-pnorm(1.96-WLS_threshold_median)

# main
WLS_threshold_median = unique(df2$WLS_value)/median(df2$SE_pcc)
1-pnorm(1.96-WLS_threshold_median)

# in main text
WLS_threshold_median = unique(df3$WLS_value)/median(df3$SE_pcc)
1-pnorm(1.96-WLS_threshold_median)

# in supplemental text
WLS_threshold_median = unique(df4$WLS_value)/median(df4$SE_pcc)
1-pnorm(1.96-WLS_threshold_median)



#############################
##### Exaggeration Bias ####
############################

#step 1: calculate WLS-FE on that sub-sample of the research record that is adequately powered 
#This yields: weighted average of the adequately powered estimator (WAAP) 

df1.p <- df1[df1$powered == 1,]
df1$WAAP <- weighted.mean(df1.p$abs_pcc, df1.p$precision_sq)
df1.u <- df1[df1$powered == 0,]
df1.u$main_pcc_waap_diff = 100*(df1.u$abs_pcc/df1.u$WAAP) 
df1.u$category <- NA
# categorize exaggeration bias
for (i in 1:nrow(df1.u)){
  if(df1.u$main_pcc_waap_diff[i] < 100){
    df1.u$category[i] <- "Deflation"
  }
  if(df1.u$main_pcc_waap_diff[i] >= 100 & df1.u$main_pcc_waap_diff[i]< 200){
    df1.u$category[i] <- "0-100%"
  }
  if(df1.u$main_pcc_waap_diff[i] >= 200 & df1.u$main_pcc_waap_diff[i]<300){
    df1.u$category[i] <- "100-300%"
  }
  if(df1.u$main_pcc_waap_diff[i] >= 200 & df1.u$main_pcc_waap_diff[i]<400){
    df1.u$category[i] <- "100-300%"
  }
  if(df1.u$main_pcc_waap_diff[i] >=400){
    df1.u$category[i] <- "300%+"
  }
  
}

nrow(df1.u)
mean(df1.u$main_pcc_waap_diff)
category_counts1 <- plyr::count(df1.u$category)
names(category_counts1) <- c("category", "frequency")
category_counts1$percentage_of_studies <- category_counts1$frequency/nrow(df1.u)
category_counts1$category <- factor(category_counts1$category, 
                                      levels = c("Deflation", "none", "0-100%", "100-300%", "300%+"))


ggsave(here("Figures", "Exaggeration_full.png"),
ggplot(aes(x = category, y = percentage_of_studies), data = category_counts1) + 
  geom_bar(stat = "identity", color = "black", fill = "gray") + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(y = "Percentage of Estiamtes", x = "Exaggeration Bias") +
  theme_pubclean() +  theme(text = element_text(size = 14), axis.title = element_text(face = "bold")),
width = 5,
height = 3,
dpi = 300
)
#######################
#### Main results ####
#######################
df2.p <- df2[df2$powered == 1,]
df2$WAAP <- weighted.mean(df2.p$abs_pcc, df2.p$precision_sq)
df2.u <- df2[df2$powered == 0,]
df2.u$main_pcc_waap_diff = 100*(df2.u$abs_pcc/df2.u$WAAP) 
df2.u$category <- NA
# categorize exaggeration bias
for (i in 1:nrow(df2.u)){
  if(df2.u$main_pcc_waap_diff[i] < 100){
    df2.u$category[i] <- "Deflation"
  }
  if(df2.u$main_pcc_waap_diff[i] >= 100 & df2.u$main_pcc_waap_diff[i]< 200){
    df2.u$category[i] <- "0-100%"
  }
  if(df2.u$main_pcc_waap_diff[i] >= 200 & df2.u$main_pcc_waap_diff[i]<300){
    df2.u$category[i] <- "100-300%"
  }
  if(df2.u$main_pcc_waap_diff[i] >= 200 & df2.u$main_pcc_waap_diff[i]<400){
    df2.u$category[i] <- "100-300%"
  }
  if(df2.u$main_pcc_waap_diff[i] >=400){
    df2.u$category[i] <- "300%+"
  }
  
}

nrow(df2.u)
mean(df2.u$main_pcc_waap_diff)
category_counts2 <- plyr::count(df2.u$category)
names(category_counts2) <- c("category", "frequency")
category_counts2$percentage_of_studies <- category_counts2$frequency/nrow(df2.u)
category_counts2$category <- factor(category_counts2$category, 
                                    levels = c("Deflation", "none", "0-100%", "100-300%", "300%+"))
ggsave(here("Figures", "Exaggeration_main.png"),
ggplot(aes(x = category, y = percentage_of_studies), data = category_counts2) + 
  geom_bar(stat = "identity", color = "black", fill = "gray") + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), n.breaks = 7) +
  labs(y = "Percentage of Estiamtes", x = "Exaggeration") +
  theme_pubclean() +  theme(text = element_text(size = 14), axis.title = element_text(face = "bold")),
width = 5,
height = 3,
dpi = 300
)

df3.p <- df3[df3$powered == 1,]
df3$WAAP <- weighted.mean(df3.p$abs_pcc, df3.p$precision_sq)
df3.u <- df3[df3$powered == 0,]
df3.u$main_pcc_waap_diff = 100*(df3.u$abs_pcc/df3.u$WAAP) 
df3.u$category <- NA
# categorize exaggeration bias
for (i in 1:nrow(df3.u)){
  if(df3.u$main_pcc_waap_diff[i] < 100){
    df3.u$category[i] <- "Deflation"
  }
  if(df3.u$main_pcc_waap_diff[i] >= 100 & df3.u$main_pcc_waap_diff[i]< 200){
    df3.u$category[i] <- "0-100%"
  }
  if(df3.u$main_pcc_waap_diff[i] > 200){
    df3.u$category[i] <- "100-300%"
  }
}


category_counts3 <- plyr::count(df3.u$category)
names(category_counts3) <- c("category", "frequency")
category_counts3$percentage_of_studies <- category_counts3$frequency/nrow(df3.u)
category_counts3$category <- factor(category_counts3$category, 
                                   levels = c("Deflation", "0-100%", "100-300%"))

ggplot(aes(x = category, y = percentage_of_studies), data = category_counts3) + 
  geom_bar(stat = "identity", color = "black", fill = "gray") + 
  scale_y_continuous(labels = scales::percent_format()) +
  labs(y = "Percentage of Estimates") + 
  theme_pubclean()

df4.p <- df4[df4$powered == 1,]
df4$WAAP <- weighted.mean(df4.p$abs_pcc, df4.p$precision_sq)
df4.u <- df4[df4$powered == 0,]
df4.u$main_pcc_waap_diff = 100*(df4.u$abs_pcc/df4.u$WAAP) 
df4.u$category <- NA
# categorize exaggeration bias
for (i in 1:nrow(df4.u)){
  if(df4.u$main_pcc_waap_diff[i] < 100){
    df4.u$category[i] <- "Deflation"
  }
  if(df4.u$main_pcc_waap_diff[i] >= 100 & df4.u$main_pcc_waap_diff[i]< 200){
    df4.u$category[i] <- "0-100%"
  }
  if(df4.u$main_pcc_waap_diff[i] > 200){
    df4.u$category[i] <- "100-300%"
  }
}


category_counts4 <- plyr::count(df4.u$category)
names(category_counts4) <- c("category", "frequency")
category_counts4$percentage_of_studies <- category_counts4$frequency/nrow(df4.u)
category_counts4$category <- factor(category_counts4$category, 
                                   levels = c("Deflation", "0-100%", "100-300%"))

ggplot(aes(x = category, y = percentage_of_studies), data = category_counts4) + 
  geom_bar(stat = "identity", color = "black", fill = "gray") + 
  scale_y_continuous(labels = scales::percent_format()) +
  labs(y = "Percentage of Estimates") + 
  theme_pubclean()

######################
## p-hacking curves ##
#####################
ggsave(here("Figures", "Unweighted_phacking.png"),
ggplot(data = df1[df1$abs_tstat_sm<10,]) +
  geom_histogram(aes(x = abs_tstat_sm,  y = ..density..), bins = 1000, fill = "gray") + 
  geom_density(aes(x = abs_tstat_sm), color = "black", kernel = "epanechnikov", bw = 0.2) + 
  scale_x_continuous(breaks = c(0,1.96,10)) +
  theme_pubclean() +
  labs(title = "Rounded T-stats", x = "|t-stat|"),
width = 5,
height = 3,
dpi = 300
)

ggsave(here("Figures", "Article_phacking.png"),
ggplot(data = df1[df1$abs_tstat_sm <10,]) +
  geom_histogram(aes(x = abs_tstat_sm,  y = ..density.., weight = weight_article), 
                 bins = 1000, fill = "gray") + 
  geom_density(aes(x = abs_tstat_sm, weight = weight_article), color = "black", kernel = "epanechnikov", bw = 0.2) + 
  scale_x_continuous(breaks = c(0,1.96,10)) +
  theme_pubclean() + 
  labs(title ="Weighted by Article", x= "|t-stat|"),
width = 5,
height = 3,
dpi = 300
)

ggsave(here("Figures", "TableWeight_phacking.png"),
  ggplot(data = df1[df1$abs_tstat_sm <10,]) +
  geom_histogram(aes(x = abs_tstat_sm,  y = ..density.., weight = weight_table), 
                 bins = 1000, fill = "gray") + 
  geom_density(aes(x = abs_tstat_sm, weight = weight_table), color = "black", kernel = "epanechnikov", bw = 0.2) + 
  scale_x_continuous(breaks = c(0,1.96,10)) +
  theme_pubclean() + theme(text = element_text(size = 14), axis.title = element_text(face = "bold")) + 
  labs(x = "Selective Reporting of t-statistic", y = "Density"), 
  width = 5,
  height = 3,
  dpi = 300
)

###################################
## Multiple hypothesis testing ##
################################

multhyp <- unique(df[,c(1,4,7,3,14,15,17,18)])
names(multhyp)<- c("uid", "journal_id", "first_author", "year_pub", "mult_hyp", "corr", "data_avail", "code")

sum(multhyp$mult_hyp)/nrow(multhyp) #84.84% of studies use multiple hypothesis testing
sum(multhyp$corr)/sum(multhyp$mult_hyp) # of those studies that use multiple hypothesis testing, 13.96% use a correction

ggsave(here("Figures", "MultHyp.png"),
ggplot(aes(x = as.factor(mult_hyp)), data = multhyp) + 
  geom_bar(aes(fill = as.factor(corr),y = (..count..)/sum(..count..)), color = "black")+
  scale_fill_manual(values = c("white", "gray")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), n.breaks = 5) +
  labs(x = "Multiple Hypothesis Testing", fill= "Correction Used", y = "Percent of Studies") + 
  theme_pubclean() +theme(axis.title = element_text(face = "bold")),
width = 3,
height = 3,
dpi = 300
)

nrow(multhyp[multhyp$data_avail ==1,])/nrow(multhyp) # 78.51% of studies have data available
ggsave(here("Figures", "DataAvail.png"),
ggplot(aes(x = as.factor(data_avail)), data = multhyp) +
  geom_bar(aes(fill = journal_id, y = (..count..)/sum(..count..)), color = "black") + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), n.breaks = 10) +
  scale_fill_brewer(palette = "Dark2") + 
  labs(x = "Data Available", fill= "Journal", y = "Percent of Studies") + 
  theme_pubclean() +   
  theme(axis.title = element_text(face = "bold"), legend.text = element_text(size = 7), 
        legend.title = element_text(size = 7), legend.key.size = unit(.5, units = "line")),
width = 3,
height = 3,
dpi = 300
)

nrow(multhyp[multhyp$code ==1,])/nrow(multhyp) # 17.5% of studies have code available

ggsave(here("Figures", "CodeAvail.png"),
ggplot(aes(x = as.factor(code)), data = multhyp) +
  geom_bar(aes(fill = journal_id, y = (..count..)/sum(..count..)), color = "black") + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), n.breaks = 8) +
  scale_fill_brewer(palette = "Dark2") + 
  labs(x = "Code Available", fill= "Journal", y = "Percent of Studies") + 
  theme_pubclean() + 
  theme(axis.title = element_text(face = "bold"), legend.text = element_text(size = 7), 
        legend.title = element_text(size = 7), legend.key.size = unit(.5, units = "line")), 
width = 3,
height = 3,
dpi = 300
)

#### multi-panel figure ####
gr1 <-ggplot(aes(x = as.factor(mult_hyp)), data = multhyp) + 
  geom_bar(aes(fill = as.factor(corr),y = (..count..)/sum(..count..)), color = "black")+
  scale_fill_manual(values = c("white", "gray")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), n.breaks = 5) +
  labs(x = "Multiple Hypothesis Testing", fill= "Correction \nUsed", y = "Percent of Studies") + 
  theme_pubclean() +theme(axis.title = element_text(face = "bold"))
gr2 <- ggplot(aes(x = as.factor(data_avail)), data = multhyp) +
  geom_bar(aes(fill = journal_id, y = (..count..)/sum(..count..)), color = "black") + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), n.breaks = 10) +
  scale_fill_brewer(palette = "Dark2") + 
  labs(x = "Data Available", fill= "Journal", y = "Percent of Studies") + 
  theme_pubclean() +   
  theme(axis.title = element_text(face = "bold"), legend.text = element_text(size = 7), 
        legend.title = element_text(size = 7), legend.key.size = unit(.5, units = "line"))
gr3 <- ggplot(aes(x = as.factor(code)), data = multhyp) +
  geom_bar(aes(fill = journal_id, y = (..count..)/sum(..count..)), color = "black") + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), n.breaks = 8) +
  scale_fill_brewer(palette = "Dark2") + 
  labs(x = "Code Available", fill= "Journal", y = "Percent of Studies") + 
  theme_pubclean() + 
  theme(axis.title = element_text(face = "bold"), legend.text = element_text(size = 7), 
        legend.title = element_text(size = 7), legend.key.size = unit(.5, units = "line"))
png(here("Figures/mult_pan.png"), width = 4, height = 9, units = "in", res = 300)
ggarrange(plotlist = list(gr1,gr2,gr3),ncol = 1,nrow =3, legend = "right", labels = c("A", "B", "C"))
dev.off()

nrow(multhyp[multhyp$code ==1 & multhyp$data_avail == 1,])/nrow(multhyp)

ggplot(aes(x = include), data = papers[papers$include %in% c("No", "Yes"),]) + 
  geom_bar(aes(fill = journal_id,), color = "black",position = position_dodge()) + 
  #scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_brewer(palette = "Dark2") + 
  labs(x = "Included", fill= "Journal", y = "Number of Studies") + 
  theme_pubclean()


#### graphs for prospectus ####
jpeg(here("Figures", "prospectus_figs.jpg"), width = 8, height = 5, units = "in", res = 300)
ggarrange(plotlist = list(gr1, gr2), labels = "AUTO")
dev.off()


unique(kicked$Checking.Notes)

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
notused <- notused[-which(is.na(notused$initials)),]
articlesnotused <- unique(notused[,c(2:6)])



#########################
### Survey Results #####
#######################

# delete first 2 rows
survey <- survey[-c(1,2),]
# take out responses that were not complete
survey <- survey[which(survey$Finished == "True"),]
# take out "do not agree"
survey <- survey[which(survey$Q1 == "I agree to participate"),]

Q2 <- as.data.frame(table(survey[,12]))
Q2$percentage <- Q2$Freq/238
Q3 <- as.data.frame(table(survey[,13]))
Q3$percentage <- Q3$Freq/238
Q4 <- as.data.frame(table(survey[,14]))
Q4 <- Q4[-which(Q4$Var1 == ""),]
Q4$percentage <- Q4$Freq/192
Q5 <- as.data.frame(table(survey[,15]))
Q5$percentage <- Q5$Freq/238
Q2.clicked <- as.data.frame(table(survey[,c(12,17)]))

ggplot(data = survey) + 
  geom_histogram(aes(y = Q2), stat = "count", color = "black") +
  ylab("What percentage of tests do you think passed \nthe conventional target of 80% power?") +
  theme_pubclean() + 
  theme(axis.title = element_text(face = "bold"), legend.text = element_text(size = 7), 
        legend.title = element_text(size = 7), legend.key.size = unit(.5, units = "line"))

ggplot(data = survey) + 
  geom_histogram(aes(y = Q3), stat = "count", color = "black") +
  ylab("Do you conduct ecological experiments?") +
  theme_pubclean() + 
  theme(axis.title = element_text(face = "bold"), legend.text = element_text(size = 7), 
        legend.title = element_text(size = 7), legend.key.size = unit(.5, units = "line"))

survey$Q4 <- factor(survey$Q4, levels = c("Never", "Less than 25% of the time", "25 - 50% of the time",
                                          "50 - 75% of the time", "75% or more of the time", "Always", "" ))
ggplot(data = survey[-which(survey$Q4 == ""),]) + 
  geom_histogram(aes(y = Q4), stat = "count", color = "black") +
  ylab("Do you perform power anaylses \nbefore starting a new experiment?") +
  theme_pubclean() + 
  theme(axis.title = element_text(face = "bold"), legend.text = element_text(size = 7), 
        legend.title = element_text(size = 7), legend.key.size = unit(.5, units = "line"))
survey$Q5 <- factor(survey$Q5, levels =c("Graduate Student", "Post-doc", "Faculty", "Researcher outside academia", "Other"))
ggplot(data = survey) + 
  geom_histogram(aes(y = Q5), stat = "count", color = "black") +
  ylab("Current position") +
  theme_pubclean() + 
  theme(axis.title = element_text(face = "bold"), legend.text = element_text(size = 7), 
        legend.title = element_text(size = 7), legend.key.size = unit(.5, units = "line"))
