################################################################################
# Author: Leonard Henckel
# Year:   2024
# Aim:    Calculate two-sample MR method estimates for the causal effect
# of habitual caffeine consumption on measures of sleep quality
# Input:  HypnoLaus and UK Biobank dataset and additional information
# Output: estimates, p-values, confidence intervals, cohen's d
################################################################################


library(AER)
library(MendelianRandomization)
library(gtools)
library(stats)
library(outliers)
library(lmtest)
library(sandwich)

source("Functions/betaComp.R")
source("Functions/MainFunction.R")

## running data preparation
source("DataPreparation/data_prep.R")

## computing SNP-treatment association to reduce computational burden
## as they do not depend on the choice of outcome variable
k <- length(snip_names)
BetaXG <- seBetaXG <- pBetaXG <- R2XG <- c()
for(i in 1:k){
  res_ukb <- beta.comp(t(snps_ukb_reorder_matrix),k+1,i,c(k+3,k+4)) 
  ## adjusting for age and sex
  BetaXG[i] <- res_ukb[1]
  seBetaXG[i] <- res_ukb[2]
  pBetaXG[i] <- res_ukb[3]
  R2XG[i] <- res_ukb[4]
}
names(BetaXG) <- names_colaus_snps
names(seBetaXG) <- names_colaus_snps
names(pBetaXG) <- names_colaus_snps


## available SNPs
index <- c("rs2472297",
  "rs2470893",
  "rs35107470",
  "rs4410790",
  "rs6968865",
  "rs2472304",
  "rs12909047",
  "rs6968554",
  "rs1992145",
  "rs62005807",
  "rs10275488",
  "rs2892838",
  "rs762551",
  "rs56113850",
  "rs7800944",
  "rs17685",
  "rs10516471",
  "rs7605062",
  "rs1800498",
  "rs2668822",
  "rs767778",
  "rs10007278",
  "rs6575353",
  "rs347306",
  "rs6279",
  "rs1571536",
  "rs66500423"
  )

# SNPs for inverse variance weighted estimator
setIVW <- index[c(2:12,16,18,22)]
# SNPs for MR-Egger estimator
setMREgger <- index[c(2,4,14:24,26,27)]
# all SNPs used for median estimator
setmedian <- index
# individual SNPs
ins1 <- index[2]
ins2 <- index[4]
ins3 <- index[16]

# risk-score SNPs
setRS <- c(ins1,ins2,ins3)

ins <- c(ins1,ins2,ins3)

results <- list()

## transforming some outcome variables to ensure constant treatment effect
## assumption is plausible
adj_snips_colaus_clean <- snips_colaus_clean
## logit for percentage data with thresholds influencing distribution
adj_snips_colaus_clean[83,] <- logit(snips_colaus_clean[83,]/100)
adj_snips_colaus_clean[84,] <- logit(snips_colaus_clean[84,]/100)
## log for strictly positive data
adj_snips_colaus_clean[86,] <- log(snips_colaus_clean[86,]+1)
adj_snips_colaus_clean[87,] <- log(snips_colaus_clean[87,]+1)
adj_snips_colaus_clean[89,] <- log(snips_colaus_clean[89,]+1)
adj_snips_colaus_clean[90,] <- log(snips_colaus_clean[90,]+1)


## running main analysis
for(i in 1:10){
  results[[i]] <- main(data=adj_snips_colaus_clean,data2=t(snps_ukb_reorder_matrix),BetaXG,seBetaXG,
                       Y=81+i,X=78,
                       W=c(80,81),
                       # W=c(),
                       ins1=ins1,ins2=ins2,ins3=ins3, # single instruments
                       setRS=setRS, # have to be valid individually 
                       setIVW=setIVW, # have to be valid individually
                       setmedian = setmedian, # at least half valid
                       setMREgger=setMREgger, # INSIDE has to hold
                       names_colaus_snps)
}

## running analysis where we consider all SNPs individually
## separated for computational reasons
source("Functions/CompareSNPsExecution.R")

## running function that prepares a results object to save
source("Functions/ResultsStorage.R")
source("Functions/CohensD.R")

## save results
save(results,summary_res,MR_IVW_Q,MR_egger_Q,tranform_flag,results_snps,CohensD, file="MR_results.RData")

