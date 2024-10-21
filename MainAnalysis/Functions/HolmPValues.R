### small file to compute Holm-Bonferroni adjusted p-values 

p_1 <- c()
for(i in c(1,2,4:10)){
  p_1 <- c(p_1,summary_res[[i]][2:4,4])
}
sort(p.adjust(p_1,method="holm"))


