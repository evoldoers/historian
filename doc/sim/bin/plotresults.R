install.packages("ggplot2")
install.packages("gridExtra")
library(ggplot2)
library(gridExtra)

data.raw <- read.csv(file="gp120sim/results.1-100.dat",sep=" ")

data <- data.raw[data$tree="tree"]

maIns <- ggplot(data[data$method=="ma.tree"&data$event=="ins",],aes(x=bin,weight=count))+geom_bar()
prankIns <- ggplot(data[data$method=="prank"&data$event=="ins",],aes(x=bin,weight=count))+geom_bar()
histIns <- ggplot(data[data$method=="hist"&data$event=="ins",],aes(x=bin,weight=count))+geom_bar()
muscleIns <- ggplot(data[data$method=="muscle"&data$event=="ins",],aes(x=bin,weight=count))+geom_bar()
probIns <- ggplot(data[data$method=="probcons"&data$event=="ins",],aes(x=bin,weight=count))+geom_bar()

maDel <- ggplot(data[data$method=="ma.tree"&data$event=="del",],aes(x=bin,weight=count))+geom_bar()
prankDel <- ggplot(data[data$method=="prank"&data$event=="del",],aes(x=bin,weight=count))+geom_bar()
histDel <- ggplot(data[data$method=="hist"&data$event=="del",],aes(x=bin,weight=count))+geom_bar()
muscleDel <- ggplot(data[data$method=="muscle"&data$event=="del",],aes(x=bin,weight=count))+geom_bar()
probDel <- ggplot(data[data$method=="probcons"&data$event=="del",],aes(x=bin,weight=count))+geom_bar()

g <- arrangeGrob(maIns,maDel,histIns,histDel,prankIns,prankDel,muscleIns,muscleDel,probIns,probDel,ncol=2)
ggsave("together.pdf",g)
