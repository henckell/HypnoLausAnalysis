################################################################################
# Author: Benjamin Stucky
# Year:   2024
# Input:  UKBiobank and HypnoLaus/CoLaus datasets
# Output: 
#         1) SNP estimated effects plot 
#         2) plot the causal estimates
#         3) create the cohen's d table
################################################################################

# libraries
library(ggplot2)                 # plotting data
library(ggpubr)
library(gridExtra)
library(grid)
library(scales)
library(dplyr)                   # tidy code
library(stringr)
library(openxlsx)                # read xlsx data
library(patchwork)


# working directory
setwd("...")

# load the MR results
load(file = "MR_results.RData")

# load the matching results
load(file ="matching_results.RData")

# load cohen's d
load(file = "CohensD.RData")


# naming
naming <- c("Total sleep time [min]", 
            "Log sleep latency [min]", 
            "Logit sleep efficiency [%]", 
            "Log # awakenings",
            "NREM Delta [%]", 
            "Logit NREM Sigma [%]",
            "REM sleep [%]",
            "Log PSQI global score", 
            "Log ESS score",
            "MEQ score")
ordr <- c(7,5,3,6,1,2,4,8,9,10)

############
# 1) SNP estimated effects plot 
# cluster coloring and naming
clustering <- rep("none", times = length(results_snps[[1]]$beta))

clustering[which(names(results_snps[[1]]$beta) %in% c("rs2472297",
                                                      "rs2470893",
                                                      "rs2472304",
                                                      "rs12909047",
                                                      
                                                      "rs35107470",
                                                      "rs1992145",
                                                      "rs62005807",
                                                      "rs762551"))] <- "CYP1A1/A2"

clustering[which(names(results_snps[[1]]$beta) %in% c("rs4410790",
                                                      "rs6968554",
                                                      "rs6968865",
                                                      "rs10275488",
                                                      "rs2892838"))] <- "AHR"

clustering[which(names(results_snps[[1]]$beta) %in% c("rs4822492",
                                                      "rs2298383",
                                                      "rs5751876"))] <- "ADORA2A"

clustering[which(names(results_snps[[1]]$beta) %in% c("rs1800498",
                                                      "rs6279"))] <- "DRD"
clustering <- factor(clustering, levels = c("CYP1A1/A2",
                                            "AHR",
                                            "DRD",
                                            "none"))


# extract the MR data from list()
out <- lapply(1:10, FUN = function(x){
  data.frame(estimate = results_snps[[x]]$beta, 
            snp_order = 1:length(results_snps[[x]]$beta),
            snp_name = names(results_snps[[x]]$beta),
            ci_low = results_snps[[x]]$CI[,1],
            ci_high = results_snps[[x]]$CI[,2],
            variable = names(results_snps)[x],
            var_name = naming[which(ordr == x)],
            clustering = clustering)
            
            
})
dat <- do.call(rbind, out)
rownames(dat) <- NULL


# set the factor order for the panel
dat$var_name <- factor(dat$var_name, levels = naming)

# create two panels: objective and subjective, through a factor
dat$panel <- "objective"
dat$panel[dat$var_name %in% c("Log PSQI global score", 
                              "Log ESS score",
                              "MEQ score")] <- "subjective"

dat$var_name <- factor(dat$var_name, levels = c("Total sleep time [min]", "Log sleep latency [min]", 
                          "Log # awakenings", "REM sleep [%]",
                          "NREM Delta [%]", "Logit NREM Sigma [%]", 
                          "Logit sleep efficiency [%]", 
                           "Log PSQI global score", "Log ESS score", "MEQ score"),
                       labels = c("Total sleep time [min]", "Sleep latency [min]", 
                                  "# Awakenings", "REM sleep [%]",
                                  "NREM delta power [%]", "NREM sigma power [%]", 
                                  "Sleep efficiency [%]", 
                                  "PSQI global score", "ESS score", "MEQ score"))



manual_colors <- c("cadetblue",
                   "mediumpurple3",
                   "salmon2",
                   "black")

# Objective panel
g1 <- ggplot(dat[dat$panel == "objective" & !(dat$var_name  %in% "Logit sleep efficiency [%]"),], 
             aes(y = estimate,
                     x = snp_order,
                     color = clustering))  +
  geom_hline(yintercept = 0, col = "gray60")+
  geom_pointrange(aes(ymin = ci_low, 
                      ymax = ci_high), show.legend = F, size = 0.7) + 
  ggtitle("Objective") +
  xlab(NULL) + 
  ylab("Estimated effect") + 
  facet_wrap(~var_name, ncol = 2, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_continuous(breaks= unique(dat[,c("snp_order", "snp_name")])[,1],
                     labels= unique(dat[,c("snp_order", "snp_name")])[,2])+
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=16,face="bold"),
        strip.text = element_text(size=14),
        panel.grid.major.x =  element_line(),
        panel.grid.minor.x =  element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.title = element_text(size=18)) +
  # scale_color_manual(values = c(gg_color_hue(4)[1:3],"black")) + 
  scale_color_manual(values = manual_colors) + 
  theme_nice(basesize = 18) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
# Subjective panel
g2 <- ggplot(dat[dat$panel == "subjective",], aes(y = estimate,
                                                 x = snp_order,
                                                 color = clustering))  +
  
  geom_hline(yintercept = 0, col = "gray60")+
  geom_pointrange(aes(ymin = ci_low, 
                      ymax = ci_high), size = 0.7)+ #, show.legend = F) + 
  ggtitle("Subjective") +
  xlab(NULL) + 
  ylab(NULL) + 
  facet_wrap(~var_name, ncol = 1, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_continuous(breaks= unique(dat[,c("snp_order", "snp_name")])[,1],
                     labels= unique(dat[,c("snp_order", "snp_name")])[,2])+
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=16,face="bold"),
        strip.text = element_text(size=14),
        panel.grid.major.x =  element_line(),
        panel.grid.minor.x =  element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.title = element_text(size=18))+
  # scale_color_manual(values = c(gg_color_hue(4)[1:3],"black")) + 
  scale_color_manual(values = manual_colors, name = "Linkage") + 
  theme_nice(basesize = 18) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

g <- g1 + g2 + plot_layout(ncol = 2, guides = "collect", widths = c(2.2,1))&
  theme(legend.position='bottom')&
  guides(color = guide_legend(nrow = 1, byrow = TRUE))

# save plot
ggsave(plot = g, filename = "tree_paper.png", 
       height = 9, width = 12, dpi = 900)




####
# create the general plotting data.frame
# extract MR information
out2 <- lapply(1:10, FUN = function(x){
  cbind(summary_res[[x]], 
        data.frame(method = rownames(summary_res[[x]]), 
                   variable = names(summary_res)[x],
                   var_name = naming[which(ordr == x)]))
  
})
dat2 <- do.call(rbind, out2)
rownames(dat2) <- NULL

# factor order for panel
dat2$var_name <- factor(dat2$var_name, levels = naming)

# create two panels: objective and subjective, through a factor
dat2$panel <- "objective"
dat2$panel[dat2$var_name %in% c("Log PSQI global score", 
                              "Log ESS score",
                              "MEQ score")] <- "subjective"

# set the method's name
dat2$method <- factor(dat2$method, levels = c("MR_egger",
                                              "MR_IVW", 
                                              "MR_median",
                                              "MR_RS", 
                                              "obervational"),
                      labels = c("MR Egger", 
                                 "IVW",
                                 "Median", 
                                 "Risk Score", 
                                 "Observational"))
colnames(dat2)[1:3] <- c("estimate", 
                    "ci_low", 
                    "ci_high")

# add the cohen's d
colnames(CohensD) <- naming[c(5,6,3,7,2,4,1,8,9,10)]
CohensD <- as.data.frame(CohensD)
CohensD <- CohensD[,naming]
CohensD$method <- c("Risk Score", "IVW", "Median", "MR Egger", "Observational")
cod <- pivot_longer(CohensD, cols = 1:10, names_to = "var_name", values_to = "cohend")
dat2 <- merge(dat2, cod, by = c("method", "var_name"))

dat2 <- bind_rows(dat2, matchings)

# omit risk score
dat2 <- dat2[dat2$method != "Risk Score",]

# set levels
dat2$method <- factor(dat2$method, levels = c("MR Egger", "IVW", "Median", "Matching", "Observational"))

dat2$cohend[dat2$method %in% c("IVW", "Median", "MR Egger", "Observational")] <- NA

# merge SE
dat2$se[dat2$method %in% "Matching"] <- dat2$SE[dat2$method %in% "Matching"]
dat2$SE <- NULL

# save the full results
write.xlsx(x = dat2, file = "results_paper.xlsx")



# split off into observational and full
dat2$var_name <- factor(dat2$var_name, levels = naming)
dat4 <- dat2
dat3 <- dat2[which(dat2$method == "Observational"), ]
# dat2 <- dat2[which(dat2$method != "Observational"), ]


######
# 2) plot the causal estimates
# Objective panel
dat2$var_name <- factor(dat2$var_name, levels = c("Total sleep time [min]", "Log sleep latency [min]", 
                                                "Log # awakenings", "REM sleep [%]",
                                                "NREM Delta [%]", "Logit NREM Sigma [%]", 
                                                "Logit sleep efficiency [%]", 
                                                "Log PSQI global score", "Log ESS score", "MEQ score"),
                        labels = c("Total sleep time [min]", "Sleep latency [min]", 
                                   "# Awakenings", "REM sleep [%]",
                                   "NREM delta power [%]", "NREM sigma power [%]", 
                                   "Sleep efficiency [%]", 
                                   "PSQI global score", "ESS score", "MEQ score"))
levels(dat2$method)[1] <- "MR-Egger"
g1 <- ggplot(dat2[dat2$panel == "objective"& !(dat2$var_name  %in% "Sleep efficiency [%]"),], 
             aes(y = estimate,
                 x = method))  +
  geom_hline(yintercept = 0, col = "red", linetype = "dashed")+
  geom_pointrange(aes(ymin = ci_low, 
                      ymax = ci_high), show.legend = F, size = 0.7) +
  ggtitle("Objective") +
  xlab(NULL) + 
  ylab("Estimated causal effect") + 
  facet_wrap(~var_name, ncol = 2, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=16,face="bold"),
        strip.text = element_text(size=14),
        panel.grid.major.x =  element_line(),
        panel.grid.minor.x =  element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.title = element_text(size=18)) + 
  theme_nice(basesize = 18) + 
  scale_y_continuous(breaks = breaks_pretty(4))
# Subjective panel
g2 <- ggplot(dat2[dat2$panel == "subjective",], 
             aes(y = estimate,
                 x = method))  +
  geom_hline(yintercept = 0, col = "red", linetype = "dashed")+
  geom_pointrange(aes(ymin = ci_low, 
                      ymax = ci_high), show.legend = F, size = 0.7) + 
  ggtitle("Subjective") +
  xlab(NULL) + 
  ylab(NULL) + 
  facet_wrap(~var_name, ncol = 1, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=16,face="bold"),
        strip.text = element_text(size=14),
        panel.grid.major.x =  element_line(),
        panel.grid.minor.x =  element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.title = element_text(size=18)) + 
  theme_nice(basesize = 18) + 
  scale_y_continuous(breaks = breaks_pretty(4))

g <- g1 + g2 + plot_layout(ncol = 2, guides = "collect", widths = c(2.15,1))&
  theme(legend.position='bottom')&
  guides(color = guide_legend(nrow = 1, byrow = TRUE))

# save plot
ggsave(plot = g, filename = "results_paper.png", 
       height = 9, width = 12, dpi = 900)
