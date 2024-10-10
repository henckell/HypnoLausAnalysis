## Code to extract core results from general results list 
## and store as a matrix
summary_res <- list()
MR_IVW_Q <- MR_egger_Q <- tranform_flag <- list()

for(i in 1:10){
  summary_res[[rownames(snips_colaus_clean)[81+i]]] <- results[[i]]$beta
  summary_res[[rownames(snips_colaus_clean)[81+i]]] <- cbind(summary_res[[i]],results[[i]]$CI)
  summary_res[[rownames(snips_colaus_clean)[81+i]]] <- cbind(summary_res[[i]],results[[i]]$pvalue)
  
  summary_res[[rownames(snips_colaus_clean)[81+i]]] <-  summary_res[[i]][4:8,]
  
  colnames(summary_res[[rownames(snips_colaus_clean)[81+i]]])<- c("beta","CI_lower","CI_higher","pvalue")
  rownames(summary_res[[rownames(snips_colaus_clean)[81+i]]])<- c("MR_RS","MR_IVW","MR_median","MR_egger","obervational")
  
  MR_IVW_Q[[rownames(snips_colaus_clean)[81+i]]] <- results[[i]][[8]]$Heter.Stat[2]
  MR_egger_Q[[rownames(snips_colaus_clean)[81+i]]] <- results[[i]][[10]]$Heter.Stat[2]
  
  tranform_flag[[rownames(snips_colaus_clean)[81+i]]] <- "none" 
}

## logit
tranform_flag[["EEG_sigmaR_nrem"]] <- "logit"
tranform_flag[["sleepefficiency"]] <- "logit"

##
tranform_flag[["latencytos2lightoff"]] <- "log"
tranform_flag[["nrofawakenings"]] <- "log"
tranform_flag[["F1psqi"]] <- "log"
tranform_flag[["F1epworth"]] <- "log"



