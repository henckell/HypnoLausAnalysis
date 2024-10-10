## function that given the data, SNP-treatment associations and SNP choices
## runs analysis for a given outcome for a variety of possible MR-methods
library(AER)
library(MendelianRandomization)
library(gtools)
library(stats)


main <- function(data,data2,BetaXG,seBetaXG,
                 Y,X,W=numeric(),ins1,ins2,ins3,
                 setRS, 
                 setIVW,
                 setmedian,
                 setMREgger,
                 names_colaus_snps){
  
  
  BetaYG <- seBetaYG <- pBetaYG <- c()
  
  for(i in 1:77){
    res <- beta.comp(data,Y,i,W)
    BetaYG[i] <- res[1]
    seBetaYG[i] <- res[2]
    pBetaYG[i] <- res[3]
  }
  
  sample_y <- sum(!is.na(data[Y,]))
  
  names(BetaYG) <- names_colaus_snps
  names(seBetaYG) <- names_colaus_snps
  names(pBetaYG) <- names_colaus_snps
  
  RSins1 <- RSins2 <- c()
  # compute genetic riskscore (RS) and add corresponding beta coefficient
  EffectDirection <- sign(BetaXG[setRS])
  RSins1 <- EffectDirection[1] * data[setRS[1],]
  RSins2 <- EffectDirection[1] * data2[setRS[1],]
  for(i in 2:length(setRS)){
    RSins1 <- RSins1 + EffectDirection[i] * data[setRS[i],]
    RSins2 <- RSins2 + EffectDirection[i] * data2[setRS[i],]
  }
  
  # now add to compute beta coefficients for genetic riskscore
  RS <- 78
  if(length(W)==0){
    res_1 <- beta.comp(rbind(RSins1,data[Y,],data[W,]),2,1,W) 
  } else{
    res_1 <- beta.comp(rbind(RSins1,data[Y,],data[W,]),2,1,3:(2+length(W))) 
  }
  
  ## make sure X and Y match the right coefficients in data!!!
  BetaYG[RS] <- res_1[1]
  seBetaYG[RS] <- res_1[2]
  if(length(W)==0){
    res_2 <- beta.comp(rbind(RSins2,data2[X,],data2[W-7,]),2,1,W) 
  } else{
    res_2 <- beta.comp(rbind(RSins2,data2[X,],data2[W-7,]),2,1,3:(2+length(W))) 
  }
  rbind(RSins2,data2[X,])
  BetaXG[RS] <- res_2[1]
  seBetaXG[RS] <- res_2[2]
  
  
  # maybe remove correlation for now if it does weird stuff
  Z <- data2[1:77,1:77]
  corr <- cor(t(Z),t(Z))
  
  MRInputObjectIns1 <- mr_input(bx = BetaXG[ins1],
                               bxse = seBetaXG[ins1],
                               by = BetaYG[ins1],
                               byse = seBetaYG[ins1])
  
  MRInputObjectIns2 <- mr_input(bx = BetaXG[ins2],
                                bxse = seBetaXG[ins2],
                                by = BetaYG[ins2],
                                byse = seBetaYG[ins2])
  
  MRInputObjectIns3 <- mr_input(bx = BetaXG[ins3],
                                bxse = seBetaXG[ins3],
                                by = BetaYG[ins3],
                                byse = seBetaYG[ins3])
  
  MRInputObjectRS <- mr_input(bx = BetaXG[RS],
                                bxse = seBetaXG[RS],
                                by = BetaYG[RS],
                                byse = seBetaYG[RS])
  
  MRInputObjectIVW <- mr_input(bx = BetaXG[setIVW],
                               bxse = seBetaXG[setIVW],
                               by = BetaYG[setIVW],
                               byse = seBetaYG[setIVW],
                               correlation=corr[setIVW,setIVW])
  
  MRInputObjectmedian <- mr_input(bx = BetaXG[setmedian],
                                  bxse = seBetaXG[setmedian],
                                  by = BetaYG[setmedian],
                                  byse = seBetaYG[setmedian],
                                  correlation=corr[setmedian,setmedian])
  
  MRInputObjectMREgger <- mr_input(bx = BetaXG[setMREgger],
                                   bxse = seBetaXG[setMREgger],
                                   by = BetaYG[setMREgger],
                                   byse = seBetaYG[setMREgger],
                                   correlation=corr[setMREgger,setMREgger])
  
  # add corrleation and use a subset 
  fit1 <- mr_ivw(MRInputObjectIns1,
                      model = "default",
                      robust = FALSE,
                      penalized = FALSE,
                      correl = FALSE,
                      weights = "simple", # "delta"
                      psi = 0, # cor(data[X,],data[Y,],"na.or.complete")
                      distribution = "normal",
                      alpha = 0.05)
  
  fit2 <- mr_ivw(MRInputObjectIns2,
                      model = "default",
                      robust = FALSE,
                      penalized = FALSE,
                      correl = FALSE,
                      weights = "simple", # "delta"
                      psi = 0, # cor(data[X,],data[Y,],"na.or.complete")
                      distribution = "normal",
                      alpha = 0.05)
  
  fit3 <- mr_ivw(MRInputObjectIns3,
                 model = "default",
                 robust = FALSE,
                 penalized = FALSE,
                 correl = FALSE,
                 weights = "simple", # "delta"
                 psi = 0, # cor(data[X,],data[Y,],"na.or.complete")
                 distribution = "normal",
                 alpha = 0.05)
  
  fitRS <- mr_ivw(MRInputObjectRS,
                      model = "default",
                      robust = FALSE,
                      penalized = FALSE,
                      correl = FALSE,
                      weights = "simple", # "delta"
                      psi = 0, # cor(data[X,],data[Y,],"na.or.complete")
                      distribution = "normal",
                      alpha = 0.05)
  
  IVWObject <- mr_ivw(MRInputObjectIVW,
                      model = "default",
                      robust = FALSE,
                      penalized = FALSE,
                      correl = TRUE,
                      weights = "simple", # "delta"
                      psi = 0, # cor(data[X,],data[Y,],"na.or.complete")
                      distribution = "normal",
                      alpha = 0.05)
  
  # # median
  WeightedMedianObject <- mr_median(MRInputObjectmedian,
                                    weighting = "weighted",
                                    distribution = "normal",
                                    alpha = 0.05,
                                    iterations = 10000,
                                    seed = 314159265)
  
  # MR-egger
  EggerObject <- mr_egger(MRInputObjectMREgger,
                          robust = FALSE,
                          penalized = FALSE,
                          correl = TRUE, # add correlation matrix
                          distribution = "normal", # might change
                          alpha = 0.05)
  
  # naive OLS 
  if(length(W==0)==0){
    fit_yx <- lm(data[Y,]~data[X,])
  } else{
    fit_yx <- lm(data[Y,]~.,data=as.data.frame(t(as.matrix(data[c(X,W),]))))
  }
  
  beta <- c()
  beta[1] <- fit1$Estimate
  beta[2] <- fit2$Estimate
  beta[3] <- fit3$Estimate
  beta[4] <- fitRS$Estimate
  beta[5] <- IVWObject$Estimate
  beta[6] <- WeightedMedianObject$Estimate
  beta[7] <- EggerObject$Estimate[1] ## make sure this is not intercept
  beta[8] <- fit_yx$coefficients[2]
  
  pvalue <- c()
  pvalue[1] <- fit1$Pvalue
  pvalue[2] <- fit2$Pvalue
  pvalue[3] <- fit3$Pvalue
  pvalue[4] <- fitRS$Pvalue
  pvalue[5] <- IVWObject$Pvalue
  pvalue[6] <- WeightedMedianObject$Pvalue
  pvalue[7] <- EggerObject$Pvalue.Est ## make sure this is not intercept
  pvalue[8] <- summary(fit_yx)$coefficients[2,4]
  
  CI <- matrix(numeric(16),ncol=2,nrow=8)
  CI[1,1] <- fit1$CILower
  CI[1,2] <- fit1$CIUpper
  CI[2,1] <- fit2$CILower
  CI[2,2] <- fit2$CIUpper
  CI[3,1] <- fit3$CILower
  CI[3,2] <- fit3$CIUpper
  CI[4,1] <- fitRS$CILower
  CI[4,2] <- fitRS$CIUpper
  CI[5,1] <- IVWObject$CILower
  CI[5,2] <- IVWObject$CIUpper 
  CI[6,1] <- WeightedMedianObject$CILower
  CI[6,2] <- WeightedMedianObject$CIUpper
  CI[7,1] <- EggerObject$CILower.Est
  CI[7,2] <- EggerObject$CIUpper.Est
  CI[8,1] <- confint(fit_yx)[2,1]
  CI[8,2] <- confint(fit_yx)[2,2]
  
  return(list(beta=beta, pvalue=pvalue, CI=CI,fit1=fit1,fit2=fit2,fit3=fit3,fitRS=fitRS,
              IVW=IVWObject,
              Median=WeightedMedianObject,
              MRegger=EggerObject,sample_y=sample_y))
  # TSHT=TSHTObject))
}


