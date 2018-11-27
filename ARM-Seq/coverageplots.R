library(ggplot2)
library(reshape2)
library(scales)
library(plyr)
library(gridExtra)

args <- commandArgs(trailingOnly = TRUE)

#args <- c("aging-coverage.txt","sacCer3-trnatable.txt", "YeastAging.txt", "aging-coverage.pdf")
#args <- c("ExosomeData-coverage.txt","/soe/holmes/pythonsource/trnatest/hgtrnadb/hg19-trnatable.txt","exosomesamples.txt","ExosomeData-SizeFactors.txt","ExosomeData-coverage.pdf")
#args

#hcvm.i <- melt(hcvm.i, id.vars=c(grep("^readC", names(hcvm.i), value=TRUE, invert=TRUE)), variable.name="tRNA.basePosition", value.name="read.density")

coverages <- read.table(args[1], header = TRUE)
trnatable <- read.table(args[2])
sampletable <- read.table(args[3])
normalizationtable <- read.table(args[4], header = TRUE)
outputfile <- args[5]
combinedfile <- "allcoverage.pdf"
multipage <- "multipage.pdf"
combinedfile = "coverage.pdf"
#colnames(coverages)

myBreaks <- function(x){
    breaks <- c(min(0),floor(max(x)))
    names(breaks) <- attr(breaks,"labels")
    breaks
}


#colnames(coverages) <- c("Feature", "Sample",1:(length(colnames(coverages)) - 2))
#trnatable[coveragemelt[,1],c(3,4)]


#remove columns with too many NAs
coverages <- coverages[ , colSums(is.na(coverages)) < nrow(coverages)/8]

#aggregate(coverages, by=list(trnatable[,3]), FUN=sum)[2]

#1:(length(colnames(coverages)) - 2)
#coveragemelt = melt(coverages, id.vars = c("Feature", "Sample"), measure.vars = as.character(colnames(coverages)[3:length(colnames(coverages)) - 2]))
coveragemelt = melt(coverages, id.vars = c("Feature", "Sample"))

#trnatable[coveragemelt[,1],c(3,4)]

#colnames(coveragemelt)

#aggregate(coveragemelt, by=list(trnatable[,3]), FUN=mean)[2]


#ddply(coveragemelt, ,summarise,value = sum(!is.na(value)))

#coveragemelt$value
coveragemelt[is.na(coveragemelt)] <- 0

#ddply(coveragemelt, c("Feature", "sample", "variable"),summarise,value = sum(!is.na(value)))

#aggregate(coveragemelt, data=dat, FUN = function(x) c(M=mean(x), SD=sd(x)))

#tapply(coveragemelt$value,sampletable[coveragemelt$Sample,2])
#( normalizationtable[2,coveragemelt$Sample])

#Normalization
#This normalization takes way too long
#coveragemelt$value <- as.vector(coveragemelt$value / as.vector( normalizationtable[1,coveragemelt$Sample]), mode = "numeric")
#coveragetest <- coveragemelt

#coveragemelt
#out <- as.vector(coveragetest$value / as.vector( normalizationtable[1,coveragetest$Sample]), mode = "numeric")

print("***")
coveragemeltagg <- aggregate(coveragemelt$value, by=list(Feature = coveragemelt$Feature, Sample = sampletable[match(coveragemelt$Sample,sampletable[,1]),2], variable = coveragemelt$variable), FUN=mean)
colnames(coveragemeltagg)[colnames(coveragemeltagg) == "x"]  <- "value"
coveragemeltagg$Sample <- factor(coveragemeltagg$Sample,levels = unique(sampletable[,2]), ordered = TRUE)
coveragemelt <- coveragemeltagg
#coveragemelt[coveragemelt$tRNA == "tRNA-Phe-GAA-2" & coveragemelt$Sample == "dmSCDd12",]
#coveragemelt[coveragemelt$tRNA == "tRNA-Phe-GAA-2" & coveragemelt$Sample == "dmSCDd1",]
#coveragemelt[coveragemelt$tRNA == "tRNA-Phe-GAA-2",]
#"41390"	"tRNA-Leu-TAA-1"	"dmMet_Amino"	"105"	546.017749321

#write.table(coveragemelt,"aggtables.txt" ,sep = "\t")


acceptorType = trnatable[match(coveragemelt$Feature, trnatable[,1]),3]
#acceptorType

#coveragemeltagg <- aggregate(coveragemelt$value, by=list(Sample = coveragemelt$Sample, variable = coveragemelt$variable), FUN=mean)
#colnames(coveragemeltagg)[colnames(coveragemeltagg) == "x"]  <- "value"
#coveragemeltagg$Sample <- factor(coveragemeltagg$Sample,levels = unique(sampletable[,2]), ordered = TRUE)

#colnames(coveragemeltagg)
#unique(coveragemelt$tRNA) 



covsummary <- ggplot(coveragemelt,aes(x=variable,y=value)) + facet_grid( ~ Sample, scales="free") + geom_bar(aes(fill = factor(acceptorType)),stat="identity")+theme(axis.title.x=element_blank(), axis.text.y=element_text(colour="black",size=8),axis.text.x=element_blank(), strip.text.y = element_text(angle=0,size=4),strip.text.x = element_text(angle=0,size=5))+ ylab("Normalized read count") +   scale_y_continuous(breaks=myBreaks) +scale_fill_discrete(drop=FALSE, name="Acceptor\ntype") #+scale_x_discrete("Position") 
ggsave(filename=combinedfile,covsummary)

#axis.text.y=element_text is y-axis coverage labels
#strip.text.y is trna names
#strip.text.x is experiment names
#ggplot(coveragemelt,aes(x=variable,y=value)) + facet_grid(Feature ~ Sample, scales="free") + geom_bar(stat="identity")+theme(axis.title.x=element_blank(), axis.text.y=element_text(colour="black",size=8),axis.text.x=element_blank(), strip.text.y = element_text(angle=0,size=2),strip.text.x = element_text(angle=0,size=8))+ ylab("read count") +   scale_y_continuous(breaks=myBreaks) #+scale_x_discrete("Position")

allcoverages <- ggplot(coveragemelt,aes(x=variable,y=value), size = 2) + facet_grid(Feature ~ Sample, scales="free") + geom_bar(stat="identity")+theme(axis.title.x=element_blank(), axis.text.y=element_text(colour="black",size=8),axis.text.x=element_blank(), strip.text.y = element_text(angle=0),strip.text.x = element_text(angle=0, size = 4))+ ylab("read count") +   scale_y_continuous(breaks=myBreaks) #+scale_x_discrete("Position") 
#ggsave(filename=outputfile, width = 30, height = 30)
#3*length(unique(coveragemelt$Feature))
#5*length(unique(coveragemelt$Sample))

plots = dlply(coveragemelt , "Feature", `%+%`, e1 = allcoverages)
ml = do.call(marrangeGrob, c(plots, list(nrow=4, ncol=1)))
ggsave(multipage, ml)

#ggsave(filename="testwrap.pdf", allcoverages,scale=2)


allcoverages <- ggplot(coveragemelt,aes(x=variable,y=value), size = 2) + facet_grid(Feature ~ Sample, scales="free") + geom_bar(stat="identity")+theme(axis.title.x=element_blank(), axis.text.y=element_text(colour="black",size=8),axis.text.x=element_blank(), strip.text.y = element_text(angle=0),strip.text.x = element_text(angle=0))+ ylab("read count") +   scale_y_continuous(breaks=myBreaks) #+scale_x_discrete("Position") 
#ggsave(filename=outputfile, width = 30, height = 30)
3*length(unique(coveragemelt$Feature))
5*length(unique(coveragemelt$Sample))
ggsave(filename=outputfile, allcoverages,scale=2,height=.5*length(unique(coveragemelt$Feature)),width=1.5*length(unique(coveragemelt$Sample)), limitsize=FALSE)
#set dpi




#ggplot(coveragemelt,aes(x=variable,y=value)) + facet_grid(Feature ~ Sample, scales="free") + geom_bar(stat="identity")+theme(axis.title.x=element_blank(), axis.text.x=element_text(colour="black",size=2), strip.text.y = element_text(angle=0))+ ylab("Normalized read count")+scale_x_discrete("Position")

 #scale_x_discrete()


#theme(axis.title.x=element_blank(), axis.text.x=element_text(colour="black"), strip.text.y = element_text(angle=0))
#+ scale_x_discrete()
#ggplot(scv.i, aes(x=paste(V1,V2,sep = "_" ), y=value)) + geom_bar(aes(fill = factor(baseMod, levels=c("m1A", "m1G", "m3C", "other bases", "not documented"))), stat="identity") + facet_grid(Feature ~ treatment, scales="free", labeller=sig_labeller_WT) + ylab("read count") + scale_x_discrete(breaks=grep("other_9|AntiC.loop_34|AntiC.loop_35|AntiC.loop_36|T.loop_58", unique(scv.i$tRNA.basePosition), value=TRUE), labels=c("res.9", "", "anti\ncodon", "","A58")) + theme(axis.title.x=element_blank(), axis.text.x=element_text(colour="black"), strip.text.y = element_text(angle=0)) + scale_fill_manual(limits=c("m1A", "m1G", "m3C", "other bases", "not documented"), values=profilePlotColors, name="Modomics\nbase\nmodifications")
