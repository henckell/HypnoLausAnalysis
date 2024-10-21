################################################################################
# Author: Benjamin Stucky
# Year:   2024
# Input:  UKBiobank and HypnoLaus/CoLaus datasets
# Output: 
#         1) objective and subjective sleep variable planel plot
#         2) distribution of cups per day plot
#         3) demographic tables
################################################################################

# libraries
library("ggplot2")      # plotting libraries
library("ggbeeswarm")
library("gridExtra")
library("grid")
library("gridtext")
library("scales")
library("tidyr")        # tidy code
library("dplyr")
library("openxlsx")     # read xlsx data
library(patchwork)


# load CoLaus data: data_colaus, UKBiobank data: data_ukb and SNP information: snp_combined
load(file = "data_UKB_Hypno.RData")

# Take out participants that wish to be withdrawn from the UK Biobank 
# They have to be omitted for the analysis
withdraw <- read.csv(file = "withdrawdata.csv",header = F)
data_ukb <- data_ukb[!(data_ukb$eid %in% as.character(as.vector(t(withdraw)))),]


# set working directory for storing plots and tables
setwd("")

# create caffeine variable, in accordance with UK Biobank
data_colaus$F1cafuse <- factor(data_colaus$F1cafuse, 
                               levels = 0:3, labels = c("0", "1-3", "4-6", ">6"))

# add additional data: Cigarettes, Alcohol
coffad <- read.xlsx(xlsxFile = "additional_data.xlsx")
data_colaus <- merge(data_colaus, coffad[,c("F1pt", "F1cgtno", "F1alcfrq")], by = "F1pt")

# function to extract number of observations
give.n <- function(x, rg){
  return(data.frame(y = Inf,
                    label = paste("",length(x), sep = ""))) 
}

# function to extract mean labels
mean.n <- function(x, rg){
  return(data.frame(y = -Inf,
                    label = paste("", round(mean(x),1), sep = "")))
}



############
# 1) objective and subjective violin panel

# total sleep time in hours
data_colaus$tsth <- data_colaus$tst/60

# transform colaus data to long data format for easy plotting
colausl <- gather(data_colaus, key = "variable", value = "yy", 
                  c("latencytos2lightoff", "sleepefficiency", 
                    "rem","nrofawakenings",
                    "EEG_deltaR_nrem", "EEG_sigmaR_nrem", 
                    "tst", "F1psqi", 
                    "F1epworth", "F1horne")
                  )
colausl$variable <- factor(colausl$variable, 
       levels = c("tst",
                  "latencytos2lightoff", 
                  "sleepefficiency", 
                  "nrofawakenings",
                  "EEG_deltaR_nrem", 
                  "EEG_sigmaR_nrem",
                  "rem",
                  "F1psqi", 
                  "F1epworth", 
                  "F1horne"),
       labels = 
         c("Total sleep time [min]",
           "Sleep latency [min]",
           "Sleep efficiency [%]",
           "# Awakenings",
           "NREM Delta [%]",
           "NREM Sigma [%]",
           "REM sleep [%]",
           "PSQI global score",
           "ESS score",
           "MEQ score")
       
       )

# create two panels: objective and subjective, through a factor
colausl$panel <- "objective"
colausl$panel[colausl$variable %in% c("PSQI global score", 
                                      "ESS score",
                                      "MEQ score")] <- "subjective"
colausl$variable <- factor(colausl$variable, levels = c("Total sleep time [min]", "Sleep latency [min]", 
                                                        "# Awakenings", "REM sleep [%]",
                                                        "NREM Delta [%]", "NREM Sigma [%]", 
                                                        "Sleep efficiency [%]", 
                                                        "PSQI global score", "ESS score", "MEQ score"),
                           labels = c("Total sleep time [min]", "Sleep latency [min]", 
                                      "# Awakenings", "REM sleep [%]",
                                      "NREM delta power [%]", "NREM sigma power [%]", 
                                      "Sleep efficiency [%]", 
                                      "PSQI global score", "ESS score", "MEQ score"))
# objective panel
g1 <- ggplot(data = colausl[!is.na(colausl$F1cafuse) &
                             colausl$panel == "objective" &
                              !(colausl$variable  %in% "Sleep efficiency [%]"),], aes(x = F1cafuse, y = yy)) +
  geom_violin(trim = T) +
  #geom_quasirandom(color = "black", alpha = 0.8) + 
  geom_boxplot(outlier.alpha = 0.3, width = 0.3) +
  ggtitle("Objective") +
  xlab("Caffeinated beverages per day") +
  ylab(NULL) +
  stat_summary(fun.data = give.n, vjust = 1.3,
               geom = "text", colour = "blue", size = 6) +
  stat_summary(fun.data = mean.n, vjust=-0.5,
              geom = "text", colour = "darkorange3", size = 6) +
  facet_wrap(~variable, ncol = 2, scales = "free_y") +
  scale_y_continuous(expand = c(0.2,0.2), breaks = pretty_breaks(n = 4))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        strip.text = element_text(size=14),
        panel.background = element_blank(),
        axis.line = element_line(),
        plot.title = element_text(size=18))+
  theme_nice(basesize = 18)

# subjective panel
g2 <- ggplot(data = colausl[!is.na(colausl$F1cafuse) &
                             colausl$panel == "subjective",], 
             aes(x = F1cafuse, y = yy)) +
  geom_violin(trim = T) +
  #geom_quasirandom(color = "black", alpha = 0.8) + 
  geom_boxplot(outlier.alpha = 0.3, width = 0.3) +
  ggtitle("Subjective") +
  xlab("Caffeinated beverages per day") +
  ylab(NULL) +
  stat_summary(fun.data = give.n, vjust = 1.3,
               geom = "text", colour = "blue", size = 6) +
  stat_summary(fun.data = mean.n, vjust=-0.5,
               geom = "text", colour = "darkorange3", size = 6) +
  facet_wrap(~variable, ncol = 1, scales = "free_y") +
  scale_y_continuous(expand = c(0.2,0.2), breaks = pretty_breaks(n = 4))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        strip.text = element_text(size=14),
        panel.background = element_blank(),
        axis.line = element_line(),
        plot.title = element_text(size=18)) +
  theme_nice(basesize = 18)

g <- g1 + g2 + plot_layout(ncol = 2, guides = "collect", widths = c(2.2,1), axis_titles = "collect")&
  theme(legend.position='bottom')&
  guides(color = guide_legend(nrow = 1, byrow = TRUE))

# save plot
ggsave(plot = g, filename = "violin_paper.png", 
       height = 9, width = 12, dpi = 900)



############
# 2) Distribution plot
# create a combined caffeine data-set
colc <- data_colaus$F1cafuse
ukb <- na.omit(as.character(data_ukb$coffeecups))
hyp <- na.omit(colc)
cuf <- factor(c(ukb,as.character(hyp)), levels = c("0","1-3", "4-6", ">6"))
dtp <- data.frame(coffeecups = cuf, 
                  cohort = c(rep("UK Biobank", times = length(ukb)),
                             rep("HypnoLaus", times = length(hyp))))
dtp$cohort <- factor(dtp$cohort, levels = c("UK Biobank", "HypnoLaus"))
dtp$Intake[dtp$coffeecups %in% c("0", "1-3")] <- "moderate"
dtp$Intake[dtp$coffeecups %in% c("4-6", ">6")] <- "high"

dis <- ggplot(dtp, aes(x= coffeecups, fill = Intake)) + geom_bar()+ 
  facet_wrap(~cohort, scales = "free_y") + 
  scale_y_continuous(labels = comma)+
  scale_x_discrete(na.translate = FALSE) + 
  ylab("Count") + 
  xlab("Caffeinated beverages per day") +
  theme(panel.background = element_blank(), 
        text = element_text(size=20),
        legend.position = "bottom")

ggsave(width =8, height =5, plot = dis,
       file = "distribution_paper.png", dpi = 900)



############
# 3) Demographics table

# calculation of table is split off into HypnoLaus and UKBiobank datasets

# part 1: HypnoLaus
# -----------------

# new variable caf
data_colaus$caf <- NA
data_colaus$caf[data_colaus$F1cafuse == "0" |data_colaus$F1cafuse == "1-3"] <- "moderate"
data_colaus$caf[data_colaus$F1cafuse == "4-6" |data_colaus$F1cafuse == ">6"] <- "high"

# total sleep time in hours
data_colaus$tsth <- data_colaus$tst/60

# levels smoking
data_colaus$F1sbsmk <- factor(data_colaus$F1sbsmk, 
                              levels = c(0,1,2), 
                              labels = c("non-smoker", "current", "former"))

# cigarettes per day, only smokers
data_colaus$cigd <- data_colaus$F1cgtno
data_colaus$cigd[data_colaus$F1sbsmk == "non-smoker"] <- NA

# alcohol use
data_colaus$F1alcfrq <- factor(data_colaus$F1alcfrq, 
                               levels = 1:7,
                               labels = c(">2 / day",
                                          "2 / day", 
                                          "1 / day",
                                          "3 - 6 / week",
                                          "1 - 2 / week",
                                          "less frequent",
                                          "never"))

# sleep, genetic and caffeine subpopulation
data_colaus <- data_colaus[!is.na(data_colaus$tst)&
                            !is.na(data_colaus$caf),]

vars <- c("F1age", 
          "F0sex", 
          "f1bmi", 
          "tsth",
          "F1horne", 
          "F1psqi", 
          "F1epworth",
          "F1sbsmk",
          "cigd",
          "F1alcfrq"
          )
varn <- c("Age [y]", 
          "Gender", 
          "BMI",
          "Total sleep Time [h]",
          "Horne-Oestberg Score", 
          "PSQI Score", 
          "Epwort Score",
          "Smoking status",
          "Cigarettes per day (only smokers)",
          "Alcohol intake frequency"
          )
# initialize
democ <- list()
for(i in 1:length(vars)){
  # extract data for moderate and high and overall
  moderate <- data_colaus[data_colaus$caf == "moderate",vars[i]]
  high <- data_colaus[data_colaus$caf == "high",vars[i]]
  ovr <- data_colaus[,vars[i]]

  # numerics and characters have a different summary procedures
  if(vars[i] %in% c("F1age", 
                    "f1bmi", 
                    "F1horne", 
                    "F1psqi", 
                    "F1epworth", 
                    "tsth",
                    "cigd")){
    # calculate t-test and wilcoxon test, chisquared test is NA
    pvl <- round(t.test(moderate, high, alternative = "two.sided")$p.value,4)
    wvl <- round(wilcox.test(moderate, high, alternative = "two.sided")$p.value,4)
    chivl <- NA
    
    # store information
    democ[[i]] <-  data.frame(variable = varn[i], 
                              frequency = "",
                              overall = paste(round(mean(ovr, na.rm= T),2), 
                                              " (", round(sd(ovr, na.rm= T),2),")",sep = ""),
                              moderate = paste(round(mean(moderate, na.rm= T),2), 
                                               " (", round(sd(moderate, na.rm= T),2),")",sep = ""),
                              high = paste(round(mean(high, na.rm= T),2), 
                                           " (", round(sd(high, na.rm= T),2),")",sep = ""),
                              difference = as.character(round(mean(high, na.rm= T)-
                                                   mean(moderate, na.rm = T),2)),
                              t_pvalue = pvl,
                              wilcox_pvalue = wvl,
                              chi_pvalue = chivl
                              )
  } else {
    # t-test and wilcoxon tests are NA
    pvl <- NA
    wvl <- NA
    
    
    # create comüarison table for chisq.test()
    gd <- expand.grid(caf = c("moderate", "high"), 
                      var = levels(data_colaus[,vars[i]]))
    tdel <- data_colaus[,c("caf", vars[i])] %>%
      group_by(eval(as.name(vars[i])) ,caf) %>%
      tally
    
    sdel <- data_colaus[,c("caf", vars[i])]  %>%
      group_by( eval(as.name(vars[i]))) %>%
      tally
    tdel <- aggregate(tdel$n, 
                      by=list(caf= tdel$caf, 
                              var = tdel$"eval(as.name(vars[i]))"), FUN=sum)
    
    tdel <- merge(tdel, gd, by = c("var", "caf"),all = T)
    tdel$x[is.na(tdel$x)] <- 0
    
    ldel <- levels(tdel$var)
    
    # calculate chisquared test
    chivl <- chisq.test(cbind(tdel$x[tdel$var == ldel[1]], 
                        tdel$x[tdel$var == ldel[2]]))$p.value
    
    mdel <- tdel[tdel$caf == "moderate", c("var", "x")]
    hdel <- tdel[tdel$caf == "high", c("var", "x")]
    
    
    # store information
    democ[[i]] <-  data.frame(variable = c(varn[i], 
                                           rep("", times = length(hdel$x)-1)), 
                              frequency = as.character(mdel$var),
                              overall = as.character(mdel$x + hdel$x),
                              moderate = as.character(mdel$x),
                              high = as.character(hdel$x),
                              difference = paste(as.character(round(100*hdel$x/sum(hdel$x)-
                                                                    100*mdel$x/sum(mdel$x),1)), "%"),
                              t_pvalue = rep(pvl, times = length(hdel$x)),
                              wilcox_pvalue = rep(wvl, times = length(hdel$x)),
                              chi_pvalue = c(round(chivl,4), 
                                             rep("", times = length(hdel$x)-1))
    )
  }
  print(varn[i])
}
# bind for loops together
democ <- do.call(rbind, democ)


# part 2: ukbiobank
# -----------------

# same caf variable as before
data_ukb$caf <- NA
data_ukb$caf[data_ukb$coffeecups == "0" |data_ukb$coffeecups == "1-3"] <- "moderate"
data_ukb$caf[data_ukb$coffeecups == "4-6"|data_ukb$coffeecups == ">6"] <- "high"

vars <- c("Age at recruitment_0.0", 
          "Sex_0.0", 
          "Body mass index (BMI)_0.0", 
          "Sleep duration_0.0", 
          "Sleeplessness / insomnia_0.0",
          "Smoking status_0.0",
          "Number of cigarettes currently smoked daily (current cigarette smokers)_0.0", 
          "Alcohol intake frequency._0.0")
varn <- c("Age [y]", 
          "Gender", 
          "BMI", 
          "Subjective sleep duration", 
          "Subjective sleeplessness",
          "Smoking status", 
          "Cigarettes per day (only smokers)",
          "Alcohol intake frequency")
# initialize
demou <- list()
for(i in 1:length(vars)){
  # extract moderate, high and overall data-sets
  moderate <- data_ukb[data_ukb$caf == "moderate",vars[i]]
  high <- data_ukb[data_ukb$caf == "high",vars[i]]
  ovr <- data_ukb[,vars[i]]
  
  # numeric and character variables have a different summary procedure
  if(vars[i] %in% c("Age at recruitment_0.0", 
                    "Body mass index (BMI)_0.0", 
                    "Number of cigarettes currently smoked daily (current cigarette smokers)_0.0", 
                    "Sleep duration_0.0")){
    # calculate t-test and wilcoxon test, chisquared test is NA
    pvl <- round(t.test(moderate, high, alternative = "two.sided")$p.value,8)
    wvl <- round(wilcox.test(moderate, high, alternative = "two.sided")$p.value,8)
    chivl <- NA
    
    # store information
    demou[[i]] <-  data.frame(variable = varn[i],
                              frequency = "",
                              overall = paste(round(mean(ovr, na.rm= T),2), 
                                              " (", round(sd(ovr, na.rm= T),2),")",sep = ""),
                              moderate = paste(round(mean(moderate, na.rm= T),2), 
                                               " (", round(sd(moderate, na.rm= T),2),")",sep = ""),
                              high = paste(round(mean(high, na.rm= T),2), 
                                           " (", round(sd(high, na.rm= T),2),")",sep = ""),
                              difference = as.character(round(mean(high, na.rm= T)-
                                                              mean(moderate, na.rm = T),2)),
                              t_pvalue = pvl,
                              wilcox_pvalue = wvl,
                              chi_pvalue = chivl
    )
  } else {
    # t-test and wilcoxon tests are NA
    pvl <- NA
    wvl <- NA
    
    # create table for chisquared test
    gd <- expand.grid(caf = c("moderate", "high"), 
                      var = levels(data_ukb[,vars[i]]))
    tdel <- data_ukb[,c("caf", vars[i])] %>%
      group_by(eval(as.name(vars[i])) ,caf) %>%
      tally
    
    sdel <- data_ukb[,c("caf", vars[i])]  %>%
      group_by( eval(as.name(vars[i]))) %>%
      tally
    tdel <- aggregate(tdel$n, 
                      by=list(caf= tdel$caf, 
                              var = tdel$"eval(as.name(vars[i]))"), FUN=sum)
    
    tdel <- merge(tdel, gd, by = c("var", "caf"),all = T)
    tdel$x[is.na(tdel$x)] <- 0
    
    ldel <- levels(tdel$var)
    
    # calculate chisquared test
    chivl <- chisq.test(cbind(tdel$x[tdel$var == ldel[1]], 
                        tdel$x[tdel$var == ldel[2]]))$p.value
    
    mdel <- tdel[tdel$caf == "moderate", c("var", "x")]
    hdel <- tdel[tdel$caf == "high", c("var", "x")]
    
    # store information
    demou[[i]] <-  data.frame(variable = c(varn[i], 
                                           rep("", times = length(hdel$x)-1)), 
                              frequency = as.character(mdel$var),
                              overall = as.character(mdel$x + hdel$x),
                              moderate = as.character(mdel$x),
                              high = as.character(hdel$x),
                              difference = paste(as.character(round(100*mdel$x/sum(mdel$x)-
                                                                    100*hdel$x/sum(hdel$x),1)), "%"),
                              t_pvalue = rep(pvl, times = length(hdel$x)),
                              wilcox_pvalue = rep(wvl, times = length(hdel$x)),
                              chi_pvalue = c(round(chivl,8), 
                                             rep("", times = length(hdel$x)-1))
    )
  }
  print(varn[i])
}
# rbind for loops
demou <- do.call(rbind, demou)


# denote the difference as high - moderate
colnames(democ)[colnames(democ) == "difference"] <- "high - moderate"
colnames(demou)[colnames(demou) == "difference"] <- "high - moderate"

# save tables
write.xlsx(democ, "demographic_hypnolaus.xlsx")
write.xlsx(demou, "demographic_ukb.xlsx")