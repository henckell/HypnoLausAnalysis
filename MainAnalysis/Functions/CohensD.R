## code to compute CohensD statistics for each of the five main methods
## and the ten outcome variables 

names_outcomes <- names(Y)

Y_cohens <- snips_colaus_clean[82:91,]

Y_cohens[which(tranform_flag=="log"),] <- log(Y_cohens[which(tranform_flag=="log"),]+1)
Y_cohens[which(tranform_flag=="logit"),] <- logit(Y_cohens[which(tranform_flag=="logit"),]/100)

var_rm <- function(X){var(X,na.rm=TRUE)}

SE <- sqrt(apply(Y_cohens,1,var_rm))

est <- c()

for(i in 1:10){
  est <- cbind(est,summary_res[[i]][,1])
}

CohensD <- est

for(i in 1:10){
  CohensD[,i] <- est[,i]/SE[i]
}

colnames(CohensD) <- names_outcomes

