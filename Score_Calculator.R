##########################################################################
#
# Project: Z-Score and MetS-Score Calculator for whole data sets
# File: Score_Calculator.R
# Author: D.Thies
# Date: 28.09.2021
#
# Purpose: Function to calculate percentiles, z-scores, MetS-scores and 
#          classification levels based on the percentiles calculated from 
#          IDEFICS study
#
# Input data files:  all_para_tables.RData
# Output data files: --
#
##########################################################################

setwd("H:/Documents/Programs/Score-Calculator/Final")

library(gamlss)
library(plyr)
library(dplyr)


# REQUIRED FUNCTIONS

# function to calculate the percentiles and z-scores. enter 'm' or 'f'
# for sex and one of the parameter abbreviations for p
paravalues <- function(df, sex, p) {

    # load all necessary parameter tables
    load(file = "all_para_tables.RData")

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
    df <- merge(df, temp, by = intersect(colnames(df), colnames(temp)),
        all.x = TRUE, suffixes = NULL)
    
    # define new columns to fill them with values if possible 
    df$new <- NA
    df$new1 <- NA

    # row-by-row calculation of the percentiles and z-scores of the
    # parameter
    for (i in 1:nrow(df)) {
        if (!is.na(eval(parse(text = paste0("df$", p, "[i]"))))) {
            if (df$dist[i] == "BCCG" & !is.na(df$dist[i])) {

                df$new[i] <- pBCCG(q = eval(parse(text = paste0("df$",
                  p, "[i]"))), mu = df$mu[i], sigma = df$sigma[i], nu = df$nu[i])
                df$new1[i] <- qNO(df$new[i])
            } else {
                if (df$dist[i] == "BCT" & !is.na(df$dist[i])) {

                  df$new[i] <- pBCT(q = eval(parse(text = paste0("df$",
                    p, "[i]"))), mu = df$mu[i], sigma = df$sigma[i], nu = df$nu[i],
                    tau = df$tau[i])
                  df$new1[i] <- qNO(df$new[i])
                } else {
                  if (df$dist[i] == "BCPE" & !is.na(df$dist[i])) {

                    df$new[i] <- pBCPE(q = eval(parse(text = paste0("df$",
                      p, "[i]"))), mu = df$mu[i], sigma = df$sigma[i],
                      nu = df$nu[i], tau = df$tau[i])
                    df$new1[i] <- qNO(df$new[i])
                  } else {
                    if (df$dist[i] == "LO" & !is.na(df$dist[i])) {

                      df$new[i] <- pLO(q = eval(parse(text = paste0("df$",
                        p, "[i]"))), mu = df$mu[i], sigma = df$sigma[i])
                      df$new1[i] <- qNO(df$new[i])
                    }
                  }
                }
            }
        }
    }

    # remove the columns that came from the parameter table
    df$dist <- NULL
    df$mu <- NULL
    df$sigma <- NULL
    df$nu <- NULL
    df$tau <- NULL

    # rename the columns appropriately
    colnames(df)[colnames(df) == "new"] <- paste("perc", p, sep = ".")
    colnames(df)[colnames(df) == "new1"] <- paste("z", p, sep = ".")

    return(df)
}

# function to calculate the MetS-score for each data point in the
# data frame. enter a data frame containing the z-scores for: waist,
# homa, sbp, dbp, trg, hdl
MetSScore <- function(df) {
    df$MetS <- df$z.waist + df$z.homa + (df$z.sbp + df$z.dbp)/2 + (df$z.trg -
        df$z.hdl)/2

    # for distribution purposes, we also need a shifted value to
    # calculate percentile and z-score
    df$mets <- df$MetS + 100
    return(df)
}

# the following functions each calculate whether a certain level is
# reached for a certain parameter. enter the data frame, the level
# name and the confidence level for each function
waistlvl <- function(df, name, c) {
    df$new <- (select(df, matches("perc.waist")) >= c)
    colnames(df)[colnames(df) == "new"] <- paste("adiposity", name, sep = ".")
    return(df)
}

bloodlvl <- function(df, name, c) {
    df$new <- (select(df, matches("perc.dbp|perc.sbp")) >= c)
    colnames(df)[colnames(df) == "new"] <- paste("blood_pressure", name,
        sep = ".")
    return(df)
}

lipidlvl <- function(df, name, c) {
    df$new <- (select(df, matches("perc.trg")) >= c | select(df, matches("perc.hdl")) <=
        1 - c)
    colnames(df)[colnames(df) == "new"] <- paste("blood_lipids", name,
        sep = ".")
    return(df)
}

glulvl <- function(df, name, c) {
    df$new <- (select(df, matches("perc.homa|perc.glu")) >= c)
    colnames(df)[colnames(df) == "new"] <- paste("blood_glu/insu", name,
        sep = ".")
    return(df)
}

# function to define for each row of a data frame whether
# monit/action level applies. enter data frame and the level to be
# checked
showlvl <- function(df, l) {
    for (i in 1:nrow(df)) {

        # return True/False for whether level is reached
        df$new[i] <- (sum(eval(parse(text = paste0("df$adiposity.", l,
            "[i]"))), eval(parse(text = paste0("df$blood_pressure.", l,
            "[i]"))), eval(parse(text = paste0("df$blood_lipids.", l, "[i]"))),
            eval(parse(text = paste0("df$`blood_glu/insu.", l, "`[i]"))),
            na.rm = TRUE) >= 3)

        # avoid 'False' if there are too many missing values
        if (sum(is.na(eval(parse(text = paste0("df$adiposity.", l, "[i]")))),
            is.na(eval(parse(text = paste0("df$blood_pressure.", l, "[i]")))),
            is.na(eval(parse(text = paste0("df$blood_lipids.", l, "[i]")))),
            is.na(eval(parse(text = paste0("df$`blood_glu/insu.", l, "`[i]"))))) >=
            2) {
            df$new[i] <- NA
        }

    }
    colnames(df)[colnames(df) == "new"] <- paste(l, "level", sep = ".")
    return(df)
}




# THE MAIN COMPOSITE FUNCTION
# function that unites the upper functions and calculates all the
# required values. enter a data frame
ScoreCalc <- function(input) {
    
    # Split the data input by sex
    datad <- split(input, input$sex)$d
    dataf <- split(input, input$sex)$f
    datam <- split(input, input$sex)$m

    # define character vector to use in the loop
    allparas <- c("bmi", "glu", "hdl", "height", "homa", "insu", "trg",
        "waist", "sbp", "dbp")
    paranames <- intersect(colnames(input), allparas)

    # calculation of percentiles and z-scores
    for (j in 1:length(paranames)) {
        if(length(datam) > 0){
            datam <- paravalues(datam, "m", paranames[j])
        }
        if(length(dataf) > 0){
            dataf <- paravalues(dataf, "f", paranames[j])
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
            datam <- paravalues(datam, "m", "mets")
            
            # deletion of the shifted value
            datam$mets <- NULL
            colnames(datam)[colnames(datam) == "MetS"] <- "mets"
        }
        
        if(length(dataf) > 0){
            
            # MetS-score calculation
            dataf <- MetSScore(dataf)
            
            # percentiles and z-scores for the MetS-score
            dataf <- paravalues(dataf, "f", "mets")
            
            # deletion of the shifted value
            dataf$mets <- NULL
            colnames(dataf)[colnames(dataf) == "MetS"] <- "mets"
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
        output <- glulvl(output, lvl[i], p[i])
        output <- showlvl(output, lvl[i])
    }

    return(output)
}

df <- read.table(file="../Example_Data3.txt", header=T)
df$height[1] <- 20

lala <- ScoreCalc(df)
