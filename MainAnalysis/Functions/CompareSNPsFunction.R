main_2 <- function(data,data2,
                 BetaXG,seBetaXG,
                 Y,W,
                 names_colaus_snps,index){
  
  
BetaYG <- seBetaYG <- beta <- c()
CI <- matrix(numeric(2*length(index)),ncol=2)

  for(i in 1:length(BetaXG)){
    res <- beta.comp(data,Y,i,W)
    BetaYG[i] <- res[1]
    seBetaYG[i] <- res[2]
  }

names(BetaYG) <- names_colaus_snps
names(seBetaYG) <- names_colaus_snps

for(j in 1:length(index)){
    MRInputObject <- mr_input(bx = BetaXG[index[j]],
                                  bxse = seBetaXG[index[j]],
                                  by = BetaYG[index[j]],
                                  byse = seBetaYG[index[j]])
    fit <- mr_ivw(MRInputObject,
                   model = "default",
                   robust = FALSE,
                   penalized = FALSE,
                   correl = FALSE,
                   weights = "simple", # "delta"
                   psi = 0, # cor(data[X,],data[Y,],"na.or.complete")
                   distribution = "normal",
                   alpha = 0.05)
    
    beta[j] <- fit$Estimate
    CI[j,1] <- fit$CILower
    CI[j,2] <- fit$CIUpper
    
}
names(beta) <- index
rownames(CI) <- index
return(list(beta=beta,CI=CI))
}

CIplot <- function(summary,title,clustering){
  beta <- summary$beta
  CI <- summary$CI
  method <- as.factor(1:27)
  
  data.frame <- data.frame(beta=beta,CIlow = CI[,1], CIhigh = CI[,2],method=method)
  
  g <- ggplot(data.frame, aes(y=beta,x=method,color=clustering))  +
    geom_pointrange(aes(ymin=CIlow, ymax=CIhigh)) + 
    ggtitle(title) + xlab("Method") + ylab("Est. effect") 
  # + scale_color_manual(values = c("#FFDB6D", "#C4961A", "#F4EDCA","#D16103", "#C3D7A4"))
  g
}
