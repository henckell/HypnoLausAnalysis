library(ggplot2)
source("Functions/CompareSNPsFunction.R")
results_snps <- list()

k <- length(snip_names)
BetaXG <- seBetaXG <- c()
for(i in 1:k){
  res_ukb <- beta.comp(t(snps_ukb_reorder_matrix),k+1,i,c(k+3,k+4))
  BetaXG[i] <- res_ukb[1]
  seBetaXG[i] <- res_ukb[2]
}

names(BetaXG) <- names_colaus_snps
names(seBetaXG) <- names_colaus_snps

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

results_snps <- list()
for(i in 1:10){
  results_snps[[rownames(snips_colaus_clean)[81+i]]] <- main_2(data=adj_snips_colaus_clean,data2=t(snps_ukb_reorder_matrix),
                            BetaXG,seBetaXG,
                            Y=81+i,W=c(80,81),
                            names_colaus_snps=names_colaus_snps,index=index)
}



