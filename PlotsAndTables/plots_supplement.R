################################################################################
# Author: Benjamin Stucky
# Year:   2024
# Input:  R2 and D' data files
# Output: D' heatmap and D' vs R2 comparison plots
################################################################################

# libraries 
library(ggplot2)          # for plotting
library(gridExtra)
library(viridis)
library(corrplot)
library(tidyverse)        # tidy code

# directory
setwd("...")

# list all the R2 and D' data files
files <- list.files(pattern = "D\\.txt$|R2\\.txt$")

# read all files at once
dat <- lapply(files, function(x)read.table(x, header=T)) 
dat <- lapply(dat, FUN = function(x){data.frame(x[,-1],row.names = x[,1])})
names(dat)<- gsub("\\..*","",files)


# long format for heatmap
dat2 <- dat %>%
  lapply(. %>% rownames_to_column) %>%
  lapply(. %>% gather(colname, value, -rowname))

# plot individual heatmaps, assign automatically to p1, p2, ....
for(i in 1:length(files)){
titles <- paste("Chromosome",gsub("[^0-9.]", "", strsplit(gsub("\\..*","",files)[i],"_")[[1]][1]),
                if(strsplit(gsub("\\..*","",files)[i],"_")[[1]][3]== "D"){"D'"}else{"R2"}, collapse = " ")

assign(paste("p",i,sep =""), ggplot(dat2[[i]], aes(x = reorder(rowname, desc(rowname)), 
                                                   y = colname, 
                                                   fill = value)) +
  geom_tile(colour="gray30",linewidth = 0.25) +
    scale_fill_gradientn(limits = c(0,1), colours = rev(viridis(10))) +
    #scale_fill_gradient(low = "skyblue3", high = "brown3")+
  xlab("") + ylab("") + labs(title = titles) + 
  scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) + # ,position = "top"
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
)
}

p <- grid.arrange(p1,p2,p11,p12,p17,p18,p19,p20,p21,p22,p23,
                 p24,p3,p4,p5,p6,p7,p8,p9,p10,p13,p14,p15,p16, 
            ncol = 2)
ggsave("LD_plot.png", plot = p ,dpi = 300, width = 10, height = 50, units = "in", limitsize = F)



# distribution d' and r2
dat22 <- list()
for(i in names(dat2)){
  tmp <- dat2[[i]]
  nams <- strsplit(i, "_")[[1]]
  tmp$measure <- nams[length(nams)]
  dat22[[i]] <- tmp
}
dat22 <- do.call(rbind, dat22)
dat22$link <- paste(dat22$rowname, dat22$colname, sep = "_")
dat22$link_opposite <- paste(dat22$colname, dat22$rowname, sep = "_")

# remove duplicates
k <- 1
while(k <= dim(dat22)[1]){
  matchs <- which(dat22$link_opposite %in% dat22$link[k])
  matchs <- matchs[matchs != k]
  if(length(matchs)>0){
    dat22 <- dat22[-matchs,]
  }
  k <- k + 1
}
dat22$measure <- factor(dat22$measure, levels = c("D", "R2"))

# plot
g <- ggplot(dat22[!dat22$link == dat22$link_opposite,], aes(x = measure, y = value)) +
  geom_violin() +
  geom_point(alpha = 0.2, size = 1) +
  geom_line(aes(group = link),
            alpha = 0.05) +
  theme_nice +
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=16,face="bold"),
        strip.text = element_text(size=14),
        panel.grid.major.x =  element_line(),
        panel.grid.minor.x =  element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.title = element_text(size=18)) +
  scale_y_continuous(breaks = c(seq(0,1, by = 0.2), 0.5)) +
  scale_x_discrete(labels = c(D = "D'", R2 = expression(r^2))) +
  stat_summary(fun.y=mean, geom="point", shape=20, size=4, color="red", fill="red", alpha = 0.9)+
  stat_summary(fun.y=median, geom="point", shape=20, size=4, color="darkorange", fill="darkorange", alpha = 0.9) +
  annotate('segment',x = 0.75, xend = 1.25, y = 0.5, yend = 0.5, linetype = "longdash", color = "blue", alpha = 0.5) +
  xlab("Linkage Disequlibrium measure") +
  ylab("Value")
ggsave(plot = g, filename = "dprime_r2_distribution_paper.png", 
       height = 7, width = 6, dpi = 900)
