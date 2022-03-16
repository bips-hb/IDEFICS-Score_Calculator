##########################################################################
# Program: Example_Calculator.R
# Project: IDEFICS
# R-Version: 4.1.2
#
# Purpose: Example use of IDEFICS_Score_Calculator.R to calculate the metabolic 
#         syndrome score for external study data 
#
# Input data files:  Example_Data.csv (Example study data)
#                    all_para_tables.RData 
#                   (The provided file includes the reference parameters)
#
# Remarks: For a step-by-step description see the readme 
#          "Example_Calculator-readme.pdf"
# 
# Author: D.Thies
# Date: 17.12.2021
##########################################################################

# install the required packages
install.packages(c("gamlss", "plyr", "dplyr"))

# set your working directory; make sure all files are in this folder
setwd("C:/Users/Person/Folder") 
#(See "Example_Calculator-readme.pdf", STEP 2)

# read the file of the study data set
# make sure that your variables have the correct names and units 
# (you can find the exact requirements in "ScoreCalc-help.pdf")
data_input <- read.csv(file = "Example_Data.csv")#(STEP 3)

# define the path of the file 'all_para_tables.RData' which is used in 
# the "ScoreCalc" function 
file_path <- "./all_para_tables.RData"#(STEP 3)

# load required functions
source("./IDEFICS_Score_Calculator.R")#(STEP 3)

# use the function ScoreCalc to get a new data set with calculated z-scores, 
# percentiles and MetS scores
z_scored_data <- ScoreCalc(data_input, file_path)

# export the table as csv-file
write.csv(z_scored_data, "./z_Scored_Data.csv")#(STEP 3)
