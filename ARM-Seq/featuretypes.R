

library(ggplot2)
library(reshape2)
library(scales)

args <- commandArgs(trailingOnly = TRUE)

#args <- c("hg19-genomictrailertable.txt","hg19-genomicbarplot.png")
#args =c("repfragtypenormcounts.txt","barplot.png")


#Rscript trailerbarplot.R hg19-trailertable.txt hg19-barplot.png
counts <- read.table(args[1])

selectcounts = counts





temp = cbind(selectcounts, seq = factor(rownames(selectcounts),rev(rownames(selectcounts)), ordered = TRUE))

#levels(temp$seq) <- rev(rownames(selectcounts))

countsmelt = melt(temp, id.vars = c('seq'))


countsmelt = within(countsmelt, seq <- factor(seq, 
    rev(rownames(selectcounts)), 
    ordered = TRUE))

ggplot(countsmelt,aes(x = variable, y = value,fill = seq, stat="identity")) +
	geom_bar(position = "fill",stat="identity") + 
    geom_bar(position = "fill",stat="identity",color="black",show_guide=FALSE) + 
    scale_y_continuous(labels = percent_format()) +
    xlab("Sample") +
    ylab("Percentage of Reads") + 
    theme(axis.title.x = element_text(face="bold", size=15), axis.text.x = element_text(face="bold", size=9,angle = 90, vjust = .5))

ggsave(filename=args[2])
    