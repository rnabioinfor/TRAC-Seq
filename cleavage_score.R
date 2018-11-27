library("data.table")
library("GenomicAlignments")
library("GenomicFeatures")
library("signal")
library("seqinr")
library("cowplot")
library(gridExtra)


args <- commandArgs(TRUE)
trnaSeqFile <- args[1]
trnaStrucFile <- args[2]
inputFile <- args[3]
tracFile <- args[4]

transcript_seqs <- read.fasta(trnaSeqFile, seqtype = 'DNA', as.string = T)
transcript_seqs <- data.table(tRNAid=names(transcript_seqs), seq=as.character(transcript_seqs))
transcript_seqs[,tRNAid:=sub("\\|.*","",tRNAid)]

data_summary <- function(tRNAid, seq) {
        len = length(unlist(strsplit(seq,"")))
        cc = data.table(tRNA=tRNAid,start=1:len)
        return(cc)
}
cct = mapply(data_summary,transcript_seqs$tRNAid,transcript_seqs$seq,SIMPLIFY = FALSE)
allsites = rbindlist(cct)

structures = setDT(read.table(trnaStrucFile, header = F))

get_coverage <- function(bamfile) {
    filename <- gsub('\\..*', '', bamfile)
    param <- ScanBamParam(flag = scanBamFlag(), simpleCigar = FALSE, reverseComplement = FALSE, tag = 'nM', tagFilter = list(), what = character(0), mapqFilter=1)
    reads <- as.data.table(readGAlignments(bamfile, index=paste(bamfile, '.bai', sep = ''), use.names=FALSE, param=param, with.which_label=FALSE))
    reads = reads[strand=="+",]
    coverage_table <- reads[, .N, by=list(seqnames,strand,start)]
    setnames(coverage_table, c('N'), c('count'))
    coverage_table[,end:=start]
    # Release memory
    rm(reads)
    gc()
    countfile = paste0(bamfile,".bg")
    bg = setDT(read.table(countfile, head=F,sep="\t",as.is=T))
    setnames(bg,c("seqnames","start","end","coverage"))
    bg[,start := start+1]
    bg[,tRNA:=gsub("\\|.*$", "", seqnames)]
    setkey(coverage_table,seqnames,start,end)
    setkey(bg,seqnames,start,end)
    tc = foverlaps(coverage_table, bg, nomatch=0L)
    tc[,c("seqnames","start","end","strand"):=NULL]
    setnames(tc,c("i.start","i.end"),c("start","end"))
    tc[,faivalue:=count/coverage]
    tc = merge(allsites,tc,by=c("tRNA","start"),all.x=T)
    tc[is.na(faivalue),faivalue := 0]
    if(filename == "WTDN2") {
      tc[count<3,faivalue := 0]
      tc[faivalue<0.1,faivalue := 0]
    }
    tc[,sample:=filename]
    return(tc)
}

tc1 = get_coverage(inputFile)
tc2 = get_coverage(tracFile)
dy = merge(tc1,tc2,by=c("tRNA","start"))
dy[,faifc:=log2(faivalue.y/faivalue.x)]
dy[is.infinite(faifc),faifc := 0]
dy[is.na(faifc),faifc := 0]
dy[,type:=gsub("chr\\S+?-(\\w+)$", "\\1", tRNA)]
dy = dy[!grep("Undet",type),]

color = c("#FF0000FF","#9932CC","#00FF66FF","#0066FFFF")
names(color) = c("codon","D_loop","A_loop","T_loop")

plot_usage2 <- function(objanticodon) {
    #objanticodon = "chr1.trna680-GluTTC"
    dat = dy[tRNA==objanticodon,]
    nucs = toupper(unlist(strsplit(transcript_seqs[tRNAid==objanticodon,seq],"")))
    my_cols = rep("#000000",length(nucs))
    objstu = structures[V1==objanticodon,]
    for(st in c("D_loop","A_loop","T_loop","codon")) {
       if(st %in% objstu$V3){
           sl = objstu[V3==st,V4]
           el = objstu[V3==st,V5]
           my_cols[sl:el] = color[st]
       }
    }
    cols <- c("codon"="#FF0000FF","D_loop"="#9932CC","A_loop"="#00FF66FF","T_loop"="#0066FFFF")
    rp = ggplot(dat, aes(x=start, y=faifc)) + geom_line() +scale_x_continuous(breaks=1:length(nucs), labels=nucs,expand = c(0, 0)) + labs(x="",y="Cleavage score",title=objanticodon) +
    theme(axis.text.x = element_text(color = my_cols, size = 6))+
    annotate("text", x = 59, y = 0, label = "D_loop",colour="#9932CC") +
    annotate("text", x = 63, y = 0, label = "A_loop",colour="#00FF66FF") +
    annotate("text", x = 67, y = 0, label = "T_loop",colour="#0066FFFF") +
    annotate("text", x = 71, y = 0, label = "Codon",colour="#FF0000FF") +
    theme(panel.background = element_blank(), axis.line = element_line(colour ="black", size = 0.8), panel.grid.major.y = element_line(colour = "gray90", size = 0.5, linetype='solid'),
    panel.grid.major.x = element_line(colour = "gray90", size = 0.5, linetype='33'),
    axis.ticks = element_line(colour = "black", size = 1.2, linetype='solid'),
    axis.text=element_text(size=12),axis.text.x=element_text(size=12),legend.position = "top",legend.title = element_blank())
    return(rp)
}

plot_tRNA <- function(tRNAtype) {
    objtRNAs = unique(dy[type==tRNAtype,tRNA])
    aa = mapply(plot_usage2,objtRNAs,SIMPLIFY=F)
    glist <- lapply(aa, ggplotGrob)
    cat(length(objtRNAs))
    cat("\n")
    if(length(glist)>15){
       glist = glist[1:15]
    }
    ggsave(paste0("tRNA_coverage_",tRNAtype,".pdf"), width = 12, height = 3*length(glist), marrangeGrob(grobs = glist, layout_matrix =matrix(1:length(glist),  ncol = 1, nrow=length(glist), byrow=TRUE)))
}
sapply(unique(dy$type),plot_tRNA)

#list output
dy = dy[faifc>6,]
setnames(dy,"tRNA","tRNAid")
dy[transcript_seqs, trnaseq:=i.seq, on=.(tRNAid)]
dy[, peakseq := toupper(substr(trnaseq, start-4, start+2))]
write.table(dy, file="trimvalue.fc.lt6.list.txt",quote=F,row.names=F,col.names=T,sep="\t")

#seq output for meme
dy[, motifseq := toupper(substr(trnaseq, start-10, start+10))]
write.fasta(sequences = as.list(dy$motifseq), names = dy$tRNAid, file.out = "for_motif.fasta") 
