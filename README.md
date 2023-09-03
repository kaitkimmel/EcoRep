---
output:
  html_document: default

---
# EcoRep
This R code and data produces the analyses in Kimmel, Avolio, & Ferraro (2023). Please contact kaitlinakimmel@gmail.com for questions regarding the code and data files. See https://osf.io/9yd2b/ for the corresponding Open Science Framework page which also includes a preprint version of the final manuscript.

## <u>Code</u>
Note: I used the here::here() function to call my data and save files. The here function makes it so that a new working directory for other users does not have to be set. When setting up the .Rproj create the following folders: "Data", "Figures", and "R". The code should run easily if the folder structure is made within a new project and here() is used. 

## <u>Code</u>
Note: I used the here::here() function to call my data and save files. The here function makes it so that a new working directory for other users does not have to be set. When setting up the .Rproj create the following folders: "Data", "Figures", and "R". I have a folder "Data" in which all the data is stored. The code should run easily if the folder structure is made within a new project and here() is used. 

+ Data will house the data files used in the analysis
+ R will house the code for analysis
+ Figures will house the output from the analysis

1. <b> Data_Cleaning.R </b>
 + inputs : RepDataJan23.csv, Paper_trackerFeb21.csv
 + outputs :  
      + Kicked.csv | .csv file of the papers that were removed from the RepDataJan23.csv after checking ;
      + CleanedDat.csv| .csv of the final cleaned RepDataJan23.csv including coefficients; 
      + CleanedPapers.csv | .csv of the final cleaned Paper_trackerFeb21.csv
 + purpose : standarizes columns for data analysis
2. <b> RepCode.R </b>
 + inputs: CleanedDat.csv, CleanedPapers.csv, Kicked.csv, RepEco_Survey.csv
 + outputs: 
      + Figure 1 - Figure 5 
      + Figure S1 & S2 [note Figure S2 is removed from final analysis as required by editor and R4]
 + purpose: the code for all the analyses presented in the manuscript. Code is commented throughout to explain the steps and purpose of each section. 

## <u>Raw Data Inputs</u>
1. <b>Paper_trackerFeb21.csv</b> is the file that recorded EVERY paper that was in the journals for our sample period it needs to be cleaned to standardize how categorical variables were coded

+ initials : initials of the research assistant who first looked at the paper
+ year_pub : publication year
+ journal_id : N = Nature, S = Science, E = Ecology, EL = Ecology Letters, JOE = Journal of Ecology
+ volume : volume number
+ issue	: issue number
+ first_author : first author last name
+ include	: Yes = Coefficients pulled in RepDataJan23.csv, No = not included
+ reason : why the paper was not included - some RAs included why they included them. This column in standardized in the Data_Cleaning.R code.

2. <b> RepDataJan23.csv</b> is the file of all the estimates from papers - it needs to be cleaned to remove data that was later excluded after double checking and to standardize how categorical variables were coded
+ initials : initials of the research assistant who first looked at the paper
+ year_pub : publication year
+ journal_id : N = Nature, S = Science, E = Ecology, EL = Ecology Letters, JOE = Journal of Ecology
+ volume : volume number
+ issue	: issue number
+ first_author : first author last name
+ table_no : table number in the publication
+ table_loc	: main text or in supplemental information
+ coefficient	: coefficient from the table
+ std_error	: error recorded from the table associated with the coefficient
+ error_type : confidence interval or standard error	
+ sample_size	: estimated sample size
+ study_type : experimental, observational, experimental and observational, or meta-analysis - note that these are standardized in the Data_cleaning.R code
+ mult_hyp (Y = 1; N = 0)	: 1 if there was multiple hypothesis testing, 0 if there was not
+ mult_hyp_corr(Y =1; N = 0) : If there was multiple hypothesis testing, 1 if there was corrections, 0 if there was not
+ main_results (Y=1; N=0)	: Was the coefficient mentioned in the main text
+ data_available : 1 if there was a data availability statement, 0 if there is not
+ code_available : 1 is there is a code availability statement, 0 if there is not
+ notes	: any notes when looking at paper
+ Checking Notes : notes from the person who checked the data, includes reason why the data was not kept if applicable
+ Keep : should the estimate be kept after checking? 1 = yes, 0 = no
+ checked	: 1 = checked, 0 = not
+ Bayesian (Y=1, N =0): 1 = bayesian stats, 0 = frequentist
+ report_power : 1 = statistical power reported, 0 = no mention of statistical power	
+ mixed_effects : was a mixed effects model used? 1 = yes, 0 = no

3. RepEco_Survey.csv is the file downloaded from Qualtrics with the survey response data

## <u>Cleaned Data </u>##
1. Kicked.csv
 + same headings as RepDataJan23.csv - see above
2. CleanedDat.csv
 + same headings as RepDataJan23.csv - see above; no initials column
3. CleanedPapers.csv
 + same headings as PaperTracker_Feb21.csv - see above
