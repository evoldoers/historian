#!/usr/bin/env R

# install.packages("ggplot2")
# install.packages("gridExtra")
library(ggplot2)
library(gridExtra)

rms.all <- read.csv(file="gp120sim/all.1-100.summary",sep=" ")
rms.sym <- read.csv(file="gp120sim/allsym.1-100.summary",sep=" ")
rms.tree <- read.csv(file="gp120sim/alltree.1-100.summary",sep=" ")

data.all <- read.csv(file="gp120sim/results.1-100.dat",sep=" ")
data.sym <- data.all[data.all$tree=="sym",]
data.tree <- data.all[data.all$tree=="tree",]

quartiles.all <- read.csv(file="gp120sim/quartiles.1-100.dat",sep=" ")
quartiles.sym <- quartiles.all[quartiles.all$tree=="sym",]
quartiles.tree <- quartiles.all[quartiles.all$tree=="tree",]

simplot <- function(method,title) {
	rms <- rms.tree
	data <- data.tree
	quartiles <- quartiles.tree
	err.mean <- rms$insdelmean[rms$method==method]
	err.rms <- rms$insdelrmse[rms$method==method]
	err.lq <- quartiles$lower[quartiles$method==method]
	err.uq <- quartiles$upper[quartiles$method==method]
	err.median <- quartiles$median[quartiles$method==method]
	return (ggplot(data[data$method==method,],aes(x=bin,weight=count))
		+geom_bar()
		+annotate("text",label=paste0(title,"\nMean error: ",sprintf("%.2f",round(err.mean-1,digits=2)),"\nRMS error: ",sprintf("%.2f",round(err.rms,digits=2))),x=1.6,y=1.75,size=3)
		+geom_segment(aes(y=3,yend=3,x=err.lq,xend=err.uq))
		+annotate("point",y=3,x=err.median)
		+xlab("Estimated indel rate / true indel rate")
		+scale_x_continuous(limits = c(0.25,1.75))
		+scale_y_continuous(limits = c(0,3))
		+theme(axis.title.y=element_blank(),
		       axis.text.y=element_blank(),
		       axis.ticks.y=element_blank()))
}

ma <- simplot("ma.tree","(a) True alignment")
histslow <- simplot("histslow","(b) Historian (default)")
hist <- simplot("hist","(c) Historian (fast)")
prank <- simplot("prank","(d) Prank")
prob <- simplot("probcons","(e) ProbCons")
muscle <- simplot("muscle","(f) Muscle")

g <- arrangeGrob(ma,histslow,hist,prank,prob,muscle,ncol=1)
ggsave("results.pdf",g)

warnings()
