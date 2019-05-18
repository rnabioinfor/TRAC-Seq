

suppressPackageStartupMessages(library("DESeq2"))



colgetlogname =  function(currtable, newname){

newtable = currtable[,c(2),drop = FALSE]
colnames(newtable)[1] <- newname
return(newtable)
}

colgetavgname =  function(currtable, newname){

newtable = currtable[,c(2),drop = FALSE]
colnames(newtable)[1] <- newname
return(newtable)
}

colrename =  function(currtable, newname){

newtable = currtable[,c(5),drop = FALSE]
colnames(newtable)[1] <- newname
return(newtable)
}


args = commandArgs(trailingOnly = TRUE)

#args =c("YeastAging","featurecounts.txt","agingshort.txt", "dmStat_Amino:dmAll_Amino","dmStat_Amino:dmMet_Amino","dmStat_Amino:dmLeu_Amino")
experimentname = args[1]
inputtable = args[2]
samplefile = args[3]


readcounts = read.table(inputtable)
sampledata = read.table(samplefile)

#sampleinfo = as.character(sampledata[colnames(readcounts) == sampledata[,1],2])
sampleinfo = as.character(sampledata[colnames(readcounts) == gsub("-", ".", sampledata[,1]) ,2])


#gsub("-", ".", sampledata[,1]) 
#sampledata[,2]
#colnames(readcounts)
sampledata[,1]
if (length(args) > 3){
comparisons = strsplit(args[4:length(args)], ":", fixed = TRUE)

}else{
comparisons = combn(unique(sampleinfo),2,simplify = FALSE)
}


#coldata
#condition


coldata = data.frame(condition=factor(sampleinfo))

coldata
cds = DESeqDataSetFromMatrix(countData = readcounts,coldata  ,design = ~ condition)
cds = DESeq(cds,betaPrior=TRUE)



names = lapply(comparisons, function(currcompare){ })

compareresults = lapply(comparisons, function(currcompare){ list(paste(currcompare[[1]],currcompare[[2]] ,sep= ":"),results( cds, contrast=c("condition", currcompare[[1]] ,currcompare[[2]]),cooksCutoff  =TRUE))})

reslist = lapply(compareresults, function(currresult){colrename(currresult[[2]],currresult[[1]])})

resloglist = lapply(compareresults, function(currresult){colgetlogname(currresult[[2]],currresult[[1]])})

resavglist = lapply(compareresults, function(currresult){colgetlogname(currresult[[2]],currresult[[1]])})                                


dds = cds



#print adjusted p-values
allprobs = Reduce(function(x,y) cbind(x,y), reslist)
write.table(allprobs,paste(experimentname,"-padjs.txt", sep = ""),sep="	")
                                                                   
#Print log values
alllogvals = Reduce(function(x,y) cbind(x,y), resloglist)
write.table(alllogvals,paste(experimentname,"-logvals.txt", sep = ""),sep="	")

#Print out the size factors
write.table(rbind(colnames(readcounts),dds$sizeFactor),file=paste(experimentname,"-SizeFactors.txt", sep = ""), row.names=FALSE,col.names=FALSE)

#get deseq normalized  raw counts
normalizedrnas = sweep(readcounts,2,dds$sizeFactor, "/" )
write.table(normalizedrnas,paste(experimentname,"-normalized.txt", sep = ""), sep = "\t")


