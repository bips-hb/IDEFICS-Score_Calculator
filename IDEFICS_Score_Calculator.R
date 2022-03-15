##########################################################################
# Program: IDEFICS_Score_Calculator.R
# Project: IDEFICS
# R-Version: 4.1.2
#
# Purpose: Functions to calculate percentiles, z-scores, IDEFICS metabolic 
#           syndrome score classification levels based publication 
#           Ahrens et al. Metabolic syndrome in young children: Definitions 
#           and results of the IDEFICS study. International Journal of Obesity. 
#           2014;38(S2):S4-S14. https://doi.org/10.1038/ijo.2014.130
#
# Input: none
# Output: none
#
# Author: D. Thies
# Date: 17.12.2021
##########################################################################

###############################################################################
# Calculate percentiles and z-scores for a given data set
###############################################################################
# The function calculates the percentiles and z-scores for each value of a 
# given clinical parameter (waist', 'sbp', 'dbp', 'trg', 'hdl' and 'homa') 
# using the respective sex-specifc reference table 
#
# Arguments: 
## data_input: data set with study data including the clinical parameters 
## sex: sex of the subjects in the data set ('m' or 'f')
## p: clinical parameter that shall be investigated (use specific abbreviation)
## tablepath: the path of the file 'all_para_tables.RData'
###############################################################################
paravalues <- function(data_input, sex, p, tablepath) { 
    
    # load all necessary parameter tables
    load(file = tablepath)

    # define names dependent on sex
    if (sex == "f") {
        s <- "girls"
    }
    if (sex == "m") {
        s <- "boys"
    }

    # save the current data frame and parameter table in a temporary
    # variable
    temp <- eval(parse(text = paste("par", p, s, sep = "_")))

    # merge the j-th parameter table to the sex-specific input table
    data_input <- merge(data_input, temp, 
                        by = intersect(colnames(data_input), colnames(temp)),
                        all.x = TRUE, suffixes = NULL)
    
    # define new columns to fill them with values if possible 
    data_input$percentile <- NA   
    data_input$z_score <- NA 

    # row-by-row calculation of the percentiles and z-scores 
    for (i in 1:nrow(data_input)) {
        if (!is.na(eval(parse(text = paste0("data_input$", p, "[i]"))))) {
            if (data_input$dist[i] == "BCCG" & !is.na(data_input$dist[i])) {

                data_input$percentile[i] <- 
                  pBCCG(q = eval(parse(text = paste0("data_input$", p, "[i]"))), 
                                                  mu = data_input$mu[i], 
                                                  sigma = data_input$sigma[i], 
                                                  nu = data_input$nu[i])
                data_input$z_score[i] <- qNO(data_input$percentile[i])
            } else {
                if (data_input$dist[i] == "BCT" & !is.na(data_input$dist[i])) {

                  data_input$percentile[i] <- 
                    pBCT(q = eval(parse(text = paste0("data_input$",p, "[i]"))), 
                                                   mu = data_input$mu[i], 
                                                   sigma = data_input$sigma[i], 
                                                   nu = data_input$nu[i],
                    tau = data_input$tau[i])
                  data_input$z_score[i] <- qNO(data_input$percentile[i])
                } else {
                  if (data_input$dist[i] == "BCPE" & 
                      !is.na(data_input$dist[i])) {

                    data_input$percentile[i] <- 
                      pBCPE(q = eval(parse(text = 
                                             paste0("data_input$",p, "[i]"))), 
                                                      mu = data_input$mu[i], 
                                                    sigma = data_input$sigma[i],
                                                      nu = data_input$nu[i], 
                                                      tau = data_input$tau[i])
                    data_input$z_score[i] <- qNO(data_input$percentile[i])
                  } else {
                    if (data_input$dist[i] == "LO" & 
                        !is.na(data_input$dist[i])) {

                      data_input$percentile[i] <- 
                        pLO(q = eval(parse(text = 
                                             paste0("data_input$",p, "[i]"))), 
                                                      mu = data_input$mu[i], 
                                                    sigma = data_input$sigma[i])
                      data_input$z_score[i] <- qNO(data_input$percentile[i])
                    }
                  }
                }
            }
        }
    }

    # remove the columns that came from the parameter table
    data_input$dist <- NULL
    data_input$mu <- NULL
    data_input$sigma <- NULL
    data_input$nu <- NULL
    data_input$tau <- NULL

    # rename the columns appropriately
    colnames(data_input)[colnames(data_input) == "percentile"] <- 
      paste("perc", p, sep = ".")
    colnames(data_input)[colnames(data_input) == "z_score"] <- 
      paste("z", p, sep = ".")

    return(data_input)
}


###############################################################################
# Calculate the MetS-score for a given data set
###############################################################################
# The function takes the z-scores of the clinical parameters 'waist', 'sbp', 
# 'dbp', 'trg', 'hdl' and 'homa' to calculate the individual MetS-score
#
# Arguments: 
## data_input: data set with z-scores of of the clinical parameters
###############################################################################
MetSScore <- function(data_input) {
    data_input$MetS <- data_input$z.waist + data_input$z.homa + 
      (data_input$z.sbp + data_input$z.dbp)/2 + 
      (data_input$z.trg - data_input$z.hdl)/2

    # for distribution purposes, we also need a shifted value to
    # calculate percentile and z-score
    data_input$MetS_shifted <- data_input$MetS + 100 
    return(data_input)
}


###############################################################################
# Derive the monitoring/action level status for component 'adiposity'
###############################################################################
# The function derive for the 'adiposity' component whether a certain limit of 
# the percentile of 'waist' has been exceeded to classify a child (none, 
# monitoring, action)
#
# Arguments: 
## data_input: data set with percentile values for adiposity component
## lvl_name: name of the classification level ('monit' or 'action')
## perc_level: the appropriate cut-off percentile for the the corresponding 
##           classification (0.9 or 0.95)
###############################################################################
waistlvl <- function(data_input, lvl_name, perc_level) { 
    data_input$adiposity <- 
      (select(data_input, matches("perc.waist")) >= perc_level)[,1]
    colnames(data_input)[colnames(data_input) == "adiposity"] <- 
      paste("adiposity", lvl_name, sep = ".")
    return(data_input)
}


###############################################################################
# Derive the monitoring/action level status for component 'blood_pressure'
###############################################################################
# The function derives for category 'blood_pressure' whether a certain limit 
# of the percentiles of 'dbp' or 'sbp' has been exceeded to classify a child 
# (none, monitoring, action)
#
# Arguments: 
## data_input: data set with percentile values for blood pressure component
## lvl_name: name of the classification level ('monit' or 'action')
## perc_level: the appropriate cut-off percentile for the the corresponding 
##           classification (0.9 or 0.95)
###############################################################################
bloodlvl <- function(data_input, lvl_name, perc_level) {
    data_input$blood_pressure <- 
      (select(data_input, matches("perc.dbp|perc.sbp")) >= perc_level)[,1]
    colnames(data_input)[colnames(data_input) == "blood_pressure"] <- 
      paste("blood_pressure", lvl_name, sep = ".")
    return(data_input)
}


###############################################################################
# Derive the monitoring/action level status for component 'blood_lipids'
###############################################################################
# The function derives for the 'blood_lipids' component whether a certain limit  
# of the percentiles of 'trg' or 'hdl' has been exceeded classifiy a child  
# (none, monitoring, action)
#
# Arguments: 
## data_input: data set with percentile values for blood lipid component
## lvl_name: name of the classification level ('monit' or 'action')
## perc_level: the appropriate cut-off percentile for the corresponding 
##           classification (0.9 or 0.95)
###############################################################################
lipidlvl <- function(data_input, lvl_name, perc_level) {
    data_input$blood_lipids <- 
      (select(data_input, matches("perc.trg")) >= perc_level | 
         select(data_input, matches("perc.hdl")) <= 1 - perc_level)[,1]
    colnames(data_input)[colnames(data_input) == "blood_lipids"] <- 
      paste("blood_lipids", lvl_name, sep = ".")
    return(data_input)
}


###############################################################################
# Derive the monitoring/action level status for component 'blood_glu_insu'
###############################################################################
# The function derives for the 'blood_glu_insu' component whether a certain  
# limit of the percentiles of 'glu' or 'insu' has been exceeded to classifiy a  
# child (none, moitoring, action)
#
# Arguments: 
## data_input: data set with percentile values for blood_glu_insu'  component
## lvl_name: name of the classification level ('monit' or 'action')
## perc_level: the appropriate cut-off percentile for the corresponding 
##           classification (0.9 or 0.95)
###############################################################################
glulvl <- function(data_input, lvl_name, perc_level) {
    data_input$blood_glu_insu <- 
      (select(data_input, matches("perc.homa|perc.glu")) >= perc_level)[,1]
    colnames(data_input)[colnames(data_input) == "blood_glu_insu"] <- 
      paste("blood_glu_insu", lvl_name, sep = ".")
    return(data_input)
}


###############################################################################
# Derive the overall monitoring/action level status
###############################################################################
# The function derives the overall monitoring/action level status based on
# on the monitoring/action level limit of the single components (adiposity etc.)
#
# Arguments: 
## data_input: data set with monitoring/action level variables of the components 
## lvl_name: name of the classification level ('monit' or 'action')
###############################################################################
showlvl <- function(data_input, lvl_name) {
    for (i in 1:nrow(data_input)) {

        # return True/False for whether level is reached
        data_input$lvl_indicator[i] <- 
          (sum(eval(parse(text = 
                            paste0("data_input$adiposity.", lvl_name, "[i]"))), 
               eval(parse(text = 
                            paste0("data_input$blood_pressure.", lvl_name,
                                   "[i]"))), 
               eval(parse(text = 
                            paste0("data_input$blood_lipids.", lvl_name, 
                                   "[i]"))),
               eval(parse(text = 
                            paste0("data_input$blood_glu_insu.", lvl_name, 
                                   "[i]"))), 
               na.rm = TRUE) >= 3)

        # avoid 'False' if there are too many missing values
        if (sum(is.na(eval(parse(text = paste0("data_input$adiposity.", 
                                               lvl_name, "[i]")))),
            is.na(eval(parse(text = paste0("data_input$blood_pressure.", 
                                           lvl_name, "[i]")))),
            is.na(eval(parse(text = paste0("data_input$blood_lipids.", 
                                           lvl_name, "[i]")))),
            is.na(eval(parse(text = paste0("data_input$blood_glu_insu.", 
                                           lvl_name, "[i]"))))) >= 2) {
            data_input$lvl_indicator[i] <- NA
        }

    }
    colnames(data_input)[colnames(data_input) == "lvl_indicator"] <- 
      paste(lvl_name, "level", sep = ".")
    return(data_input)
}



###############################################################################
# Calculate the metabolic syndrome score for a whole data set 
###############################################################################
# The function calculates the metabolic syndrome score for a whole data set 
# using all previously defined functions.
# Requires the packages 'gamlss', 'plyr' and 'dplyr'.
#
# Arguments: 
## data_set: data set including the study data with clinical parameters
## tablepath: the path of the file 'all_para_tables.RData'
###############################################################################
ScoreCalc <- function(data_set, tablepath) { 
    
    #install.packages(c("gamlss", "plyr", "dplyr"))
    library(gamlss)
    library(plyr)
    library(dplyr)
    
    #add a temporary enumeration to sort by at the end
    data_set$Temp_Number <- seq(1:nrow(data_set))
    
    # split the data input by sex
    datad <- split(data_set, data_set$sex)$d
    dataf <- split(data_set, data_set$sex)$f
    datam <- split(data_set, data_set$sex)$m

    # define character vector to use in the loop
    allparas <- c("bmi", "glu", "hdl", "height", "homa", "insu", "trg",
        "waist", "sbp", "dbp")
    paranames <- intersect(colnames(data_set), allparas)

    # calculation of percentiles and z-scores
    for (j in 1:length(paranames)) {
        if(length(datam) > 0){
            datam <- paravalues(datam, "m", paranames[j], tablepath)
        }
        if(length(dataf) > 0){
            dataf <- paravalues(dataf, "f", paranames[j], tablepath)
        }
    }
    
    # define vector with parameters necessary to calculate the
    # MetS-score
    necparas <- c("waist", "homa", "sbp", "dbp", "trg", "hdl")

    # MetS-score calculation if all necessary parameters exist in the
    # data set
    if (sum(necparas %in% paranames) == 6) {
        
        if(length(datam) > 0){
            
            # MetS-score calculation
            datam <- MetSScore(datam)
            
            # percentiles and z-scores for the MetS-score
            datam <- paravalues(datam, "m", "MetS_shifted", tablepath)
            
            # deletion of the shifted value
            datam$MetS_shifted <- NULL 
            
            #rename the percentile and z-score columns of MetS_shifted to MetS
            colnames(datam)[colnames(datam) == "perc.MetS_shifted"] <- 
              "perc.MetS"
            colnames(datam)[colnames(datam) == "z.MetS_shifted"] <- "z.MetS"
        }
        
        if(length(dataf) > 0){
            
            # MetS-score calculation
            dataf <- MetSScore(dataf)
            
            # percentiles and z-scores for the MetS-score
            dataf <- paravalues(dataf, "f", "MetS_shifted", tablepath)
            
            # deletion of the shifted value
            dataf$MetS_shifted <- NULL
            
            #rename the percentile and z-score columns of MetS_shifted to MetS
            colnames(dataf)[colnames(dataf) == "perc.MetS_shifted"] <- 
              "perc.MetS"
            colnames(dataf)[colnames(dataf) == "z.MetS_shifted"] <- "z.MetS"
        }
        
    }

    # reunite the (edited) data into one data frame
    output <- rbind.fill(datam, dataf, datad) 

    # calculate monitoring and action level for the parameters
    lvl <- c("monit", "action")
    p <- c(0.9, 0.95)
    for (i in 1:2) {
        output <- waistlvl(output, lvl[i], p[i])
        output <- bloodlvl(output, lvl[i], p[i])
        output <- lipidlvl(output, lvl[i], p[i])  
        #(for 'hdl' the percentile level is then 1 - perc_level)
        output <- glulvl(output, lvl[i], p[i])
        output <- showlvl(output, lvl[i])
    }
    
    #sort by temporary enumeration and delete it afterwards
    output <- output[order(output$Temp_Number),]
    output$Temp_Number <- NULL
    
    return(output)
}


