################################################################################
# Author: Benjamin Stucky
# Year:   2024
# Aim:    Calculate the MatchingFrontier Algorithm in the HypnoLaus data-set
# Input:  HypnoLaus/CoLaus dataset and additional information
# Output: matching p-value, confidence interval, cohen's d
################################################################################

# Libraries
library(tidyverse)        # easy calculations
library(MatchIt)          # matching toolkit
library(MatchingFrontier) # MatchingFrontier Algorithm, install via remotes::install_github('IQSS/MatchingFrontier')
library(openxlsx)         # read xlsx data
library(haven)            # read dta data
library(effsize)          # cohen d function
library(ggplot2)          # plotting

# set working directory
setwd("...")

# load the data set
load(file = "data_UKB_CL_v2020_coffcups.RData")

# add caffeine and other information
coffad <- read.xlsx(xlsxFile = "additional_data.xlsx")
which_add <- c("F1pt", colnames(coffad)[!colnames(coffad) %in% colnames(data_colaus)])
data_colaus <- merge(data_colaus, coffad[,which_add], by = "F1pt")

# add socio-demographic-economic and health related data
socioeco_add <- read_dta("data/Data_20231106_additional_matchingfrontier/Caffeine_genetics_data2.dta")
which_add <- c("F1pt", colnames(socioeco_add)[!colnames(socioeco_add) %in% colnames(data_colaus)])
data_colaus <- merge(data_colaus, socioeco_add[,which_add], by = "F1pt")

# which variables to use for matching individuals
matching_vars <- c("F0sex", # gender
                   "F1age", # age
                   
                   "F1conso_hebdo", # alcohol units per week
                   "F1equiv",       # smoking cigarettes equivalent
                   
                   # ethnicity
                   "ethori_self",
                   
                   # education
                   "edtyp",
                   
                   # martial status
                   "F1mrtsts",
                   "F1dmst",
                   "F1nochd",
                   "F1agechd1",
                   "F1agechd2",
                   "F1sclhlp",
                   
                   # JOB STATUS
                   "F1job_curr1",
                   "F1job_curr4b",
                   "F1job_curr8",
                   "F1job_curr8a",
                   
                   # personal history (cardiovascular illnesses, strokes,...)
                   "F1cmp",
                   "F1hdv",
                   "F1chf",
                   "F1artm",
                   "F1cad",
                   "F1angn",
                   "F1miac",
                   "F1strk",
                   "F1vslg",
                   "F1ccth",
                   "F1cabg",
                   
                   # health status
                   "F1Quest1",
                   
                   # physical activity weekly excluding sleep
                   "F1etsemns",
                   
                   # mental health
                   "F1SC",
                   "F1DA",
                   "F1IR",
                   "F1PA",
                   "F1CESD",
                   
                   # BMI, waist, hip
                   "F1waist",
                   "F1hip",
                   "f1bmi", # BMI
                   
                   # resting blood pressure, HR
                   "F1SBP",
                   "F1DBP",
                   "F1HRTRTE",
                   
                   # cholesterol
                   "F1chol",
                   "F1hdlch",
                   "F1ldlch",
                   "F1trig",
                   
                   # inflammation
                   "F1crpu"
                   )

# create caffeine intake variable
data_colaus$F1cafuse <- factor(data_colaus$F1cafuse, 
                               levels = 0:3, labels = c("0", "1-3", "4-6", ">6"))
data_colaus$caf <- NA
data_colaus$caf[data_colaus$F1cafuse == "0" |data_colaus$F1cafuse == "1-3"] <- "moderate"
data_colaus$caf[data_colaus$F1cafuse == "4-6" |data_colaus$F1cafuse == ">6"] <- "high"

# total sleep time from minutes in hours
data_colaus$tsth <- data_colaus$tst/60

# outcome variables and their corresponding names
naming <- c("Total sleep time [min]", 
            "Log sleep latency [min]", 
            "Logit sleep efficiency [%]", 
            "Log # awakenings",
            "NREM Delta [%]", 
            "Logit NREM Sigma [%]",
            "REM sleep [%]",
            "Log PSQI global score", 
            "Log ESS score",
            "MEQ score")
outcoms <- c("tst",
             "latencytos2lightoff", 
             "sleepefficiency", 
             "nrofawakenings",
             "EEG_deltaR_nrem", 
             "EEG_sigmaR_nrem",
             "rem",
             "F1psqi", 
             "F1epworth", 
             "F1horne")

# transform outcome variables if necessary
data_colaus[, c("latencytos2lightoff","nrofawakenings","F1psqi","F1epworth")] <- log(data_colaus[, c("latencytos2lightoff","nrofawakenings","F1psqi","F1epworth")])
data_colaus[, c("sleepefficiency", "EEG_sigmaR_nrem")] <- car::logit(data_colaus[, c("sleepefficiency", "EEG_sigmaR_nrem")]/100)

# convert variables to factor if necessary
data_colaus[,setdiff(matching_vars, "F1nochd")] <- lapply(data_colaus[,setdiff(matching_vars, "F1nochd")],\(x){
  if("character" %in% class(x)){
    as.factor(x)
  } else if("haven_labelled" %in% class(x)){
    factor(x, levels = attr(x, "labels"), labels = names(attr(x, "labels")))
  } else {
    x
  }
})



# MatchingFrontier

# initialize
matchingsave <- list()

# loop over all outcome variables
for(outcome in 1:length(outcoms)){
  # omit na's and non-finite values
  data_match <- na.omit(data_colaus[, c("caf", outcoms[outcome],matching_vars)])
  data_match <- data_match[is.finite(data_match[,outcoms[outcome]]), ]
  
  # caffeine as the main factor, make sure moderate corresponds to 0
  data_match$caf <- factor(data_match$caf, levels = c("moderate", "high"))
  
  # make the frontier
  formula_frontier <- as.formula(paste("caf", paste(matching_vars, collapse = " + "), sep = " ~ "))
  mahal.frontier <- makeFrontier(formula_frontier,
                                 data = data_match, 
                                 QOI = "FSATE", 
                                 metric = "dist",
                                 distance.mat = "mahalanobis", 
                                 verbose = FALSE)

  
  # create the estimates (only for plotting purpose)
  basef <- as.formula(paste(outcoms[outcome], "~", "caf"))
  mahal.estimates <- estimateEffects(mahal.frontier, 
                                     base.form = basef,
                                     verbose = FALSE)
  
  
  # save "check" plots
  p1 <- plot(mahal.frontier) +
    ggtitle(paste(naming[outcome], ": frontier plot"))
  p2 <- plot(mahal.estimates, band = "confidence") +
    geom_hline(yintercept = mahal.estimates$un$coef, linetype = "dashed")+
    geom_hline(yintercept = 0)+
    ggtitle(paste(naming[outcome], ": effects plot"))
  ggsave(plot = p1,
         filename = paste0("plots/matchingfrontier_frontier_new_", outcoms[outcome], ".jpg"))
  ggsave(plot = p2,
         filename = paste0("plots/matchingfrontier_cibands_new_", outcoms[outcome], ".jpg"))
  
  # generate matching data
  datamatch <- generateDataset(mahal.frontier, 
                               Ndrop = floor(dim(mahal.frontier$data)[1]/2), # drop the worst 50% matches
                               dup = T)                                      # to retain the matching id's (subclass)

  
  # merge the correct matches for easy t-test calculations
  caf1 <- datamatch[datamatch$caf == "1", c("subclass", outcoms[outcome])]
  caf0 <- datamatch[datamatch$caf == "0", c("subclass", outcoms[outcome])]
  caf10 <- merge(caf1,caf0, by = "subclass", suffixes = c("_1", "_0"))

  
  # calculate the paired two-sided t-tests
  tout <- t.test(caf10[, paste(outcoms[outcome], "1", sep = "_")], 
                 caf10[, paste(outcoms[outcome], "0", sep = "_")], 
         paired = T, alternative = "two.sided")
  
  # calculate the according cohen's d value
  dout <- cohen.d(caf10[, paste(outcoms[outcome], "1", sep = "_")], 
                  caf10[, paste(outcoms[outcome], "0", sep = "_")],
                  paired = T)
  
  # store all information
  matchingsave[[outcome]] <- data.frame(variable = outcoms[outcome],
                                        var_name = naming[outcome],
                                        estimate = tout$estimate,
                                        pvalue = tout$p.value,
                                        cohend = dout$estimate,
                                        ci_low = tout$conf.int[1],
                                        ci_high = tout$conf.int[2],
                                        panel = if(outcoms[outcome] %in% c("F1psqi", "F1epworth", "F1horne")){
                                          "subjective"
                                        } else {"objective"},
                                        method = "Matching")
  
}
# bind the rows of the for loops together
matchings <- bind_rows(matchingsave)
rownames(matchings) <- NULL

