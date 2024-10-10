## function to compute adjusted SNP-outcome and SNP-treatment associations
beta.comp <- function(data1,y,x,Z=numeric(0),errors="robust"){
  BetaYG    <- c()
  seBetaYG    <- c()
  pBetaYG <- c()
  R2YG <- c()
  
  missing <- which(is.na(data1[y,]))
  
  if(length(missing)==0){
    if(length(Z)==0){
      fitYZ <- lm(data1[y,]~data1[x,])
    } else{
    fitYZ <- lm(data1[y,]~.,data=as.data.frame(t(as.matrix(data1[c(x,Z),]))))
    }
  } else{
    if(length(Z)==0){
      fitYZ <- lm(data1[y,-missing]~data1[x,-missing])
    } else{
    fitYZ <- lm(data1[y,-missing]~.,data=as.data.frame(t(as.matrix(data1[c(x,Z),-missing]))))
    }
  }
  BetaYG <- fitYZ$coefficients[2]
  R2YG <- summary(fitYZ)$r.squared
  
  if(errors=="standard"){
  seBetaYG <- sqrt(diag(vcov(fitYZ)))[2]
  pBetaYG <- summary(fitYZ)[[4]][2,4]
  } else{ 
    res <- coeftest(fitYZ, vcov = vcovHC(fitYZ))
    seBetaYG<-res[2,2]
    pBetaYG<-res[2,4]
  }
  return(c(BetaYG, seBetaYG,pBetaYG,R2YG))
}