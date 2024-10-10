### small file to compute Holm-Bonferroni adjusted p-values 

## only accounting for estimates included in paper
p_1 <- c()
for(i in c(1,2,4:10)){
  p_1 <- c(p_1,summary_res[[i]][2:4,4])
}
sort(p.adjust(p_1,method="holm"))

## also accounting for sleep efficiency as additional outcome and RS as estimator
p_2 <- c()
for(i in 1:10){
  p_2 <- c(p_2,summary_res[[i]][1:4,4])
}
sort(p.adjust(p_2,method="holm"))


