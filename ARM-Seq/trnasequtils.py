#!/usr/bin/env python

import pysam
import sys
import tempfile
from collections import defaultdict

def tempmultifasta(allseqs):
    fafile = tempfile.NamedTemporaryFile(suffix=".fa")


    for seqname, seq in allseqs:
        fafile.write(">"+seqname+"\n")
        fafile.write(seq+"\n")
    fafile.flush()
    return fafile
def invertstrand(strand):
    if strand == "+":
        return "-"
    elif strand == "-":
        return "+"

class alignment:
    def __init__(self, sequences):
        if len(sequences.keys()) < 1:
            raise EmptyAlignException() 
        self.aligns = sequences

        if max(len(curr) for curr in self.aligns.itervalues()) != min(len(curr) for curr in self.aligns.itervalues()):
            print >>sys.stderr, "Non-matching sequence lengths in multiple alignment"
        #self.alignlength = len(self.aligns[self.aligns.keys()[0]])
        self.alignlength = max(len(curr) for curr in self.aligns.itervalues())
    def getseqlength(self, seqname):
        return len(self.aligns[seqname].replace("-",""))
    def numseqs(self):
        return len(self.aligns.keys())
    def getseqnames(self):
        return self.aligns.keys()
    def getalignseq(self ,seqname):
        return self.aligns[seqname]
    def toseqdict(self):
        seqs = dict()
        for currname, currseq in self.aligns.iteritems():
            seqs[currname] = string.translate(currseq, None, "-.~")
        return seqs
    def getsubset(self, subset):
        newaligns = dict()
        for curr in subset:
            newaligns[curr] = ""
        for i in range(self.alignlength):
            if any(self.aligns[curr][i] not in gapchars for curr in subset):
                for curr in subset:
                    newaligns[curr] += self.aligns[curr][i]
        return alignment(newaligns)
    def prottocodonalignment(self, nucseqs):
        #I don't do any checking to make sure that the nucleotide sequence is equivalent here
        newalign = defaultdict(str)
        currpoint = defaultdict(int)
        for i in range(self.alignlength):
            for currgene in nucseqs.iterkeys():
                if self.aligns[currgene][i] != "-":
                    newalign[currgene] += nucseqs[currgene][currpoint[currgene]:currpoint[currgene] + 3]
                    currpoint[currgene] += 3
                else:
                    newalign[currgene] += "---"
         
        #print "|".join(str(curr) for curr in currpoint.values())
        #print "|".join(str(len(curr)) for curr in nucseqs.values())
        #print "|".join(curr for curr in newalign.values())
        
        #print str(self.alignlength) pyroIsla1  
        return alignment(newalign, aligned = True)
    def tempmultifasta(self):
        return tempmultifasta(self.aligns)
    def fastaformat(self, prefix = ""):
        output = ""
        for seqname, seq in self.aligns.iteritems():
            output += ">"+prefix+seqname+"\n"+seq+"\n"
        return output
    def getsegment(self,start, end):
        newalign = dict()
        for name, seq in self.aligns.iteritems():
            newalign[name] = seq[start:end]
        return alignment(newalign)
    def getsubsets(self,windowsize, stepsize):
        for curr in range(0, self.alignlength - windowsize, stepsize):
            start = curr
            end = curr + windowsize
            yield self.getsegment(start,end)
    def removeemptyseqs(self):
        newalign = dict()
        for name, seq in self.aligns.iteritems():
            if len(string.translate(seq, None, "-.")) != 0:
                newalign[name] = seq
        return alignment(newalign)
    #http://bioinformatics.oxfordjournals.org/content/25/5/668.full

    def phylipformat(self):
        #print ",".join(str(curr)for curr in self.aligns.keys())
        output = str(len(self.aligns.keys()))+ " "+str(self.alignlength)+"\n"
        for currgene in self.aligns.iterkeys():
            #sys.stderr.write("**"+str(currgene).ljust( 11)+"**\n")
            output += str(currgene).ljust( 14)
            output +=  self.aligns[currgene]+"\n"
        #print output
        return output


    def nexusformat(self):
        output = "#NEXUS\nBegin data;\nDimensions ntax="+str(len(self.aligns.keys()))+" nchar="+str(self.alignlength)+";\n"
        output += "Format datatype=dna symbols=\""+"ATCG"+"\" missing=? gap=-;\n"
        output += "Matrix\n"
        for currgene in self.aligns.iterkeys():
            output += currgene.ljust( 11)
            output +=  self.aligns[currgene]+"\n"
        output+=";\nEnd;"    
    def getrealcoordinate(self, alignname, coord):
        currreal = 0
        #print self.aligns[alignname]
        #print coord
        for i in range(0, coord):
            if self.aligns[alignname][i] not in  set("-.~"):
                currreal += 1
        return min([currreal, len(list(curr for curr in self.aligns[alignname] if curr not in set("-.~"))) - 1])
    def getseqrange(self, alignname, start, end):
        return self.getsegment(self.getaligncoordinate(alignname, start),self.getaligncoordinate(alignname, end))

            
    def printstk(self, name = None):
        print "# STOCKHOLM 1.0"
        if name is not None:
            print "#=GF ID "+name
        padlength = max(len(curr) for curr in self.aligns.iterkeys()) + 4
        for currname, currseq in self.aligns.iteritems():
            print string.ljust(currname,padlength ) +currseq
        print "//"
    def printhtml(self, name = None):
        print "<CODE>"
        for currname, currseq in self.aligns.iteritems():
            print currname + "\t"+currseq +"<BR/>"
        print "</CODE>"
    def clustalwformat(self):
        output = "CLUSTAL W 2.1 multiple sequence alignment\n\n"
        conservestring = ""
        for i in range(0, self.alignlength):
            conservestring += ":"                     
        for currpos in range(0, self.alignlength, 60):
            for seqname, seq in self.aligns.iteritems():
                output += seqname+"\t"+seq[currpos:min([currpos+60,self.alignlength])]+"\n"
            output += "\t"+conservestring[currpos:min([currpos+60,self.alignlength])] + "\n\n"
            
        return output
        

class RnaAlignment(alignment):
    def __init__(self,alignseqs, structure, consensus = None,energies = None):
        self.aligns = alignseqs
        self.currstruct = structure
        self.energies = energies 
        self.consensus = consensus
        self.alignlength = max(len(curr) for curr in alignseqs.values())
    def addupstream(self, seqs, struct = None):
        newseqs = dict()
        #print >>sys.stderr, seqs.keys()
        for curr in self.aligns.iterkeys():
            newseqs[curr] = seqs[curr] + self.aligns[curr]
        if struct is None:
            newstruct = (max(len(curr) for curr in seqs.itervalues()) * ":") + self.currstruct
        else:
            newstruct = struct + self.currstruct
        return RnaAlignment(newseqs, newstruct)
    def adddownstream(self, seqs, struct = None):
        newseqs = dict()
        for curr in self.aligns.iterkeys():
            newseqs[curr] =  self.aligns[curr]+  seqs[curr]
        if struct is None:
            newstruct = self.currstruct + (max(len(curr) for curr in seqs.itervalues()) * ":")
        else:
            newstruct = self.currstruct + struct
        return RnaAlignment(newseqs, newstruct)
    def addmargin(self, length):
        newseqs = dict()
        for curr in self.aligns.iterkeys():
            newseqs[curr] =  length*"N" + self.aligns[curr]+ length*"N"
        newstruct = length * ":" + self.currstruct + length * ":"

        return RnaAlignment(newseqs, newstruct)

    def viennaformat(self):
        output = ""
        for currseq in self.aligns.iterkeys():
            output += ">"+currseq+"\n"
            output += self.aligns[currseq]+"\n"
            output += self.currstruct+"\n"
        return output
    def viennatempfile(self):
        viennafile = tempfile.NamedTemporaryFile()
        viennafile.write(self.viennaformat())
        viennafile.flush()
        return viennafile
    def printstk(self, name = None):
        print "# STOCKHOLM 1.0"
        if name is not None:
            print "#=GF ID "+name
        for currname, currseq in self.aligns.iteritems():
            print currname + "\t"+currseq
        structline = ""
        structpos = 0
        secpos = 0
        print "#=GC SS_cons\t"+self.currstruct
        print "//"
    def printhtml(self, name = None):
        for currname, currseq in self.aligns.iteritems():
            print currname + "\t"+currseq +"</BR>"
        structline = ""
        structpos = 0
        secpos = 0
        print "#=GC SS_cons\t"+self.currstruct+"</BR>"
        print "//"+"</BR>"
def convertmaturealign(rnaalign):
    newseqs = dict()
    newalign = rnaalign.struct
    for name, seq in rnaalign.aligns.iteritems():
        newseqs[name]
    
def readrnastk(stk):
    seqs = defaultdict(str)
    struct = ""
    consensus = ""
    energyscore = None
    for line in stk:
        line = line.rstrip()
        if line.startswith("//"):
            if consensus == "":
                consensus = None
            yield RnaAlignment(seqs, struct, consensus = consensus)
            seqs = defaultdict(str)
            struct = ""
            consensus = ""
            energyscore = None
        elif not line.startswith("#") and len(line.split()) > 1:
            currname = line.split()[0]
            currseq = line.split()[1]
            seqs[currname] += currseq
        elif line.startswith("#=GC SS_cons"):
            struct += line.split()[2]
        elif line.startswith("#=GC RF"):
            consensus += line.split()[2]
            

class transcriptfile:
    def __init__(self, trnafilename):
        trnafile = open(trnafilename)
        locustranscript = dict()
        trnatranscripts = list()
        for i, line in enumerate(trnafile):
            fields = line.split()
            if len(fields) < 2:
                continue
            trnatranscripts.append(fields[0])
            for currlocus in fields[1].split(','):
                locustranscript[currlocus] = fields[0]

        
        self.locustranscript = locustranscript
        self.transcripts = trnatranscripts
    def gettranscripts(self):
        return set(self.transcripts)
        


class samplefile:
    def __init__(self, samplefilename):
        samplefile = open(samplefilename)
        samplelist = list()
        samplefiles = dict()
        replicatename = dict()
        
        replicatelist = list()
        for i, line in enumerate(samplefile):
            fields = line.split()
            if len(fields) < 2:
                continue
            samplefiles[fields[0]] = fields[2]
            replicatename[fields[0]] = fields[1]
            
            samplelist.append(fields[0])
            if fields[1] not in set(replicatelist):
                replicatelist.append(fields[1])
        
        #bamlist = list(curr + "_sort.bam" for curr in samplefiles.iterkeys())
        samplenames = list(curr  for curr in samplefiles.iterkeys())
        self.samplelist = samplelist
        self.samplefiles = samplefiles
        self.replicatename = replicatename
        self.replicatelist = replicatelist
        #self.bamlist = list(curr+ "_sort.bam" for curr in samplelist)
    def getsamples(self):
        return self.samplelist
    def getbamlist(self):
        return list(curr+ ".bam" for curr in samplelist)
    def getbam(self, sample):
        #return sample+"/"+sample+ ".bam"
        return sample+ ".bam"
    def getfastq(self, sample):
        return self.samplefiles[sample]
    def getreplicatename(self, sample):
        return self.replicatename[sample]
    def allreplicates(self):
        return self.replicatelist
    def getrepsamples(self, replicate):
        return set(currsample for currsample in self.replicatename.iterkeys() if self.replicatename[currsample] == replicate)
        
def getsizefactors( sizefactorfile):
    sizefactorfile = open(sizefactorfile)
    sizefactors = dict()
    bamheaders =list()
    sizes = list()
    for i, line in enumerate(sizefactorfile):
        if i == 0:
            bamheaders = list(curr.strip("\"\n") for curr in line.split())
        elif i == 1:
            sizes = list(float(curr.strip("\"\n")) for curr in line.split())
    for i in range(0, len(bamheaders)):
        sizefactors[bamheaders[i]] = sizes[i]
        #print >>sys.stderr, bamheaders[i]+":"+ str(sizes[i])
    return sizefactors
    
class GenomeRange:
    __slots__ = "dbname", "chrom", "strand","name", "fastafile"
    def __eq__(self, other):
        return self.strand == other.strand and self.chrom == other.chrom and self.start == other.start and self.end == other.end
    def __hash__(self):
        return  self.start + self.end + hash(self.chrom) + hash(self.strand)
    def __init__(self, dbname, chrom, start, end, strand = None,name = None, orderstrand = False, data = None, fastafile = None):
        self.dbname =dbname
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.fastafile = fastafile
        if orderstrand and self.start > self.end:
            temp = self.start
            self.start = self.end
            self.end = temp
            self.strand = invertstrand(strand)
        self.data = data
        self.name = name        
    def coverage(self, other):
        if self.strand == other.strand and self.chrom == other.chrom:
            start = max([self.start,other.start])
            end = min([self.end,other.end])
            
            if end - start < 0:
                return 0
            else:
                return end - start
        else:
            return 0
            
    def bedstring(self, name = None,score = 1000):
        if self.strand == None:
            self.strand = "+"
        if name is None and self.name is None:
            name = "FEAT"
        elif name is None:
            name = self.name
        return "\t".join([self.chrom,str(self.start),str(self.end),name,str(score), self.strand])
    def length(self):
        return self.end - self.start
    def addmargin(self, dist = 50, name = None):
        newname = name
        if name is None:
            newname = self.name
        return GenomeRange(self.dbname, self.chrom, self.start - dist,self.end + dist,self.strand, name = newname, fastafile = self.fastafile)


'''
Still need to add trailer fragment code, both here and elsewhere
'''
def getfragtype(currfeat, currread, maxoffset = 10):
    if currread.start < currfeat.start + maxoffset and currread.end > currfeat.end - maxoffset:
        return "Whole"
    elif currread.start < currfeat.start + maxoffset:
        if currfeat.strand == "+":
            return "Fiveprime"
        else:
            return "Threeprime"
    elif currread.end > currfeat.end - maxoffset:
        if currfeat.strand == "+":
            return "Threeprime"
        else:
            return "Fiveprime"
def readfeatures(filename, orgdb="genome", seqfile= None):
    if filename.endswith(".bed") or filename.endswith(".bed.gz"):
        return readbed(filename, orgdb, seqfile)
    elif filename.endswith(".gtf") or filename.endswith(".gtf.gz") or filename.endswith(".gff") or filename.endswith(".gff.gz"):
        
        return readgtf(filename, orgdb, seqfile)
    else:
        print >>sys.stderr, filename+" not valid feature file"
        sys.exit()


def readgtf(filename, orgdb="genome", seqfile= None):
    bedfile = None
    #print >>sys.stderr, "****"
    if filename == "stdin":
        bedfile = sys.stdin
    elif filename.endswith(".gz"):
        bedfile = gzip.open(filename, 'rb')
    else:
        bedfile = open(filename, "r")
    
    for currline in bedfile:
        #print currline
        if currline.startswith('track') or currline.startswith('#'):
            continue
        fields = currline.rstrip().split("\t")
        if len(fields) > 2:
            featname = None
            #print >>sys.stderr, len(fields)
            if fields[2] != "transcript":
                continue
                
                
            for currattr in fields[8].rstrip(";").split(";"):
                #print >>sys.stderr,  currattr
                currname = currattr.strip().split()[0]
                currvalue = currattr.strip().split()[1]
                if currname == "name" or currname == "gene_id":
                    featname = currvalue.strip('"')
            #print >>sys.stderr,  GenomeRange( orgdb, fields[0],fields[3],fields[4],fields[6], name = featname, fastafile = seqfile).bedstring()                    
            yield GenomeRange( orgdb, fields[0],fields[3],fields[4],fields[6], name = featname, fastafile = seqfile)
            
def readbed(filename, orgdb="genome", seqfile= None):
    bedfile = None
    if filename == "stdin":
        bedfile = sys.stdin
    elif filename.endswith(".gz"):
        bedfile = gzip.open(filename, 'rb')
    else:
        bedfile = open(filename, "r")
    for currline in bedfile:
        #print currline
        if currline.startswith('track') or currline.startswith('#'):
            continue
        fields = currline.rstrip().split()
        if len(fields) > 2:
            #print >>sys.stderr, len(fields)
            if len(fields) < 5:
                strand = "+"
            else:
                #print >>sys.stderr, len(fields)
                strand = fields[5]
                
            yield GenomeRange( orgdb, fields[0],fields[1],fields[2],strand, name = fields[3], fastafile = seqfile)

def ifelse(arg, trueres,falseres):
    if arg:
        return trueres
    else:
        return falseres
        
        
def isprimarymapping(mapping):
    return not (mapping.flag & 0x0100 > 0)        
def getbamrange(bamfile, chromrange = None, primaryonly = False):
    bamiter = None
    if chromrange is not None:
        bamiter = bamfile.fetch(chromrange.chrom, chromrange.start, chromrange.end)
    else:
        bamiter = bamfile.fetch()
    
    for currline in bamiter:
        rname = bamfile.getrname(currline.rname)
        #print rname
        #need to fix this with cigar stuff
        #len(currline.pos)
        strand = "+"
        strand = ifelse(currline.is_reverse, '-','+')
        seq = None
        #print currline.cigar
        if primaryonly and not isprimarymapping(currline):
            continue
        if strand == "-":
            pass
        yield GenomeRange( "genome",rname,currline.pos,currline.aend,strand, name = currline.qname , data = {"score":currline.mapq, "CIGAR":currline.cigar, "seq":seq, "flags": currline.flag, "qual":currline.qual})


def getpileuprange(bamfile, chromrange = None):
    bamiter = None
    if chromrange is not None:
        bamiter = bamfile.pileup(chromrange.chrom, chromrange.start, chromrange.end)
    else:
        bamiter = bamfile.pileup()
    
    for currpos in bamiter:
        readcounts = defaultdict(int)
        reference = 'N'
        for read in currpos.pileups:
            aln = read.alignment
            
            
            if not read.indel and not read.is_del:
                readcounts[aln.seq[read.qpos]] += 1
        yield read.qpos, readcounts
                
        
def getseqdict(genelist, faifiles = None):
    namedict = getnamedict(genelist)
    allorgs = set(currgene.dbname for currgene in genelist)
    dbdict = dict()
    fastafiles = dict()
    for currorg in allorgs:
        dbdict[currorg] = dict()
        
    
    for currgene in genelist:
        dbdict[currgene.dbname][currgene.name] = currgene
        if currgene.fastafile is not None:
            fastafiles[currgene.dbname] =  currgene.fastafile
            
        else:
            #print >>sys.stderr, currgene.dbname+":"+genomefile(currorg)
            fastafiles[currgene.dbname] =  genomefile(currgene.dbname)        
    seqdict = dict()
    
    for currorg in allorgs:
        
        if faifiles is not None:
            
            currseqs = getseqs(fastafiles[currorg], dbdict[currorg], faindex = faifiles[currorg])
        else:
            currseqs = getseqs(fastafiles[currorg], dbdict[currorg])
        seqdict.update(currseqs)
    return seqdict
    
def getnamedict(genelist):
    namedict = dict()
    for currgene in genelist:
        namedict[currgene.name] = currgene
    return namedict
        
        
def getseqs(fafile,rangedict, faindex = None):
    #print >>sys.stderr, rangedict
    if faindex is not None:
        try:
            faifile = fastaindex(fafile, faindex)
        except IOError as e:
            print >>sys.stderr, "Cannot read fasta file "+fafile
            print >>sys.stderr, "Ensure that file "+fafile +" exits and generate fastaindex "+faindex+" with samtools faidx"
            sys.exit(1)
        return faifile.getseqs(rangedict)
    genomefile = open(fafile, "r")
    reheader = re.compile( r"\>([^\s\,]+)")
    allseqs = defaultdict(str)
    currloc = 0
    for line in genomefile:
        line = line.rstrip("\n")
        currheader = reheader.match(line)
        if currheader: #sequence += line[chromstart - currloc:chromend - currloc]
            currseq = currheader.groups(1)[0]
            #print >>sys.stderr, currseq
            currloc = 0
        else:
            for currname, location in rangedict.iteritems():
                if currseq == location.chrom:
                    
                    chromstart = location.start
                    chromend = location.end
                    if location.dbname == 'eschColi_K12':
                        pass
                        #print >>sys.stderr, genomefile
                        #print >>sys.stderr, currseq+":"+str(chromstart)+"-"+str(chromend)+":"+str(currloc)
                    if currloc <= chromstart <= currloc + len(line) and currloc <= chromend <= currloc + len(line):
                        allseqs[currname] += line[chromstart - currloc:chromend - currloc]
                    elif currloc <= chromstart <= currloc + len(line):
                        allseqs[currname] += line[chromstart - currloc:]
                    elif currloc <= chromend <= currloc + len(line):
                        allseqs[currname] += line[:chromend - currloc]
                    elif  currloc < chromstart <chromend < currloc + len(line):
                        pass
                    elif chromstart <= currloc < currloc + len(line) < chromend:
                        allseqs[currname] += line
            currloc += len(line)
    genomefile.close()
    finalseqs = dict()
    for currname in allseqs.iterkeys():
        #allseqs[currname] = allseqs[currname].upper()
        if (rangedict[currname].strand == "-"):
            seq = list(allseqs[currname].upper())
            seq.reverse()
            comp = {"A":"T","T":"A", "C":"G","G":"C","N":"N","R":"Y","Y":"R","S":"W","W":"S", "K":"M", "M":"K"}
            finalseqs[currname]  = ''.join(comp[base] for base in seq)
        else:
            finalseqs[currname] = allseqs[currname].upper()
    for currseq in rangedict.iterkeys():
        if currseq not in finalseqs:
            print >>sys.stderr, "No sequence extracted for "+rangedict[currseq].dbname+"."+rangedict[currseq].chrom+":"+str(rangedict[currseq].start)+"-"+str(rangedict[currseq].end)
    return finalseqs        
    
class fastaindex:
    def __init__(self, fafile, faifile):
        self.fafile = fafile
        fai = open(faifile)
        self.chromsize = dict()
        self.chromoffset = dict()
        self.seqlinesize = dict()
        self.seqlinebytes = dict()
        
        for line in fai:
            fields = line.split("\t")
            self.chromsize[fields[0]] = int(fields[1])
            self.chromoffset[fields[0]] = int(fields[2])
            self.seqlinesize[fields[0]] = int(fields[3])
            self.seqlinebytes[fields[0]] = int(fields[4])
    def getchrombed(self, dbname = 'genome'):
        for curr in self.chromsize.iterkeys():
            yield GenomeRange(dbname,curr,0,self.chromsize[curr],name=curr, strand = "+")
    def getseek(self, currchrom,loc):
        #print >>sys.stderr, (self.seqlinebytes[currchrom] - self.seqlinesize[currchrom])
        return self.chromoffset[currchrom] + loc + int(loc/(self.seqlinesize[currchrom]))*(self.seqlinebytes[currchrom] - self.seqlinesize[currchrom])
    def getfullseqs(self, names):
        genomefile = open(self.fafile, "r")
        for currchrom in names:
            #print >>sys.stderr, currchrom+":"+str(self.chromsize[currchrom])
            genomefile.seek(self.getseek(currchrom,0))
            #seq = genomefile.read(self.getseek(currchrom,self.chromsize[currchrom]))
            seq = genomefile.read(self.getseek(currchrom,self.chromsize[currchrom]) - self.getseek(currchrom,0))
            #seq = seq.replace("\n","")
            yield currchrom, seq 
    def getseqs(self,  rangedict):
        genomefile = open(self.fafile, "r")
        allseqs = dict()
        for currname, currregion in rangedict.iteritems():
            currchrom = currregion.chrom
            #faskip = 
            #print >>sys.stderr, int(currregion.start/(self.seqlinebytes[currchrom] - self.seqlinesize[currchrom]))
            genomefile.seek(self.getseek(currchrom,currregion.start))
            seq = genomefile.read(self.getseek(currchrom,currregion.end) - self.getseek(currchrom,currregion.start))
            seq = seq.replace("\n","")
            allseqs[currname] = seq
            #print >>sys.stderr, len(seq)
            #print >>sys.stderr, str(currregion.end - currregion.start)
            
        genomefile.close()
        finalseqs = dict()
        for currname in allseqs.iterkeys():
            #allseqs[currname] = allseqs[currname].upper()
            if (rangedict[currname].strand == "-"):
                seq = list(allseqs[currname].upper())
                seq.reverse()
                comp = {"A":"T","T":"A", "C":"G","G":"C","N":"N","R":"Y","Y":"R","S":"W","W":"S", "K":"M", "M":"K"}
                finalseqs[currname]  = ''.join(comp[base] for base in seq)
            else:
                finalseqs[currname] = allseqs[currname].upper()
        return finalseqs    
        
        
        
#the object I use for storing groups of genomeranges.  Useful for finding overlaps
class RangeBin:     
    def __init__(self,rangelist, binfactor = 10000):
        self.binfactor = binfactor
        self.bins = []
        self.length = 0
        for curr in rangelist:
            self.additem(curr)
    def __len__(self):
        return self.length
    def __iter__(self):
        for currbin in self.bins:
            #print currbin
            for currgene in currbin:
                yield currgene
    def additem(self, item):
        binstart = int(item.start / self.binfactor)
        #print "**"+str(binstart)
        binend = int(item.end / self.binfactor) + 1
        while (binstart + 2 >= len(self.bins)):
            self.bins.append(set())
        self.bins[binstart].add(item)
        self.length += 1
        #print self.bins[binstart]
    def getrange(self, item):
        for i in range(int(item.start / self.binfactor)-1,int(item.end / self.binfactor)+1):
            if i < len(self.bins):
                for currrange in self.bins[i]:
                    if currrange.start >= item.start and currrange.end <= item.end:
                            yield currrange
                            
    def getbin(self, item):
        for i in range(int(item.start / self.binfactor)-1,int(item.end / self.binfactor)+1):
            if i < len(self.bins) and i >= 0:
                
                for currrange in self.bins[i]:
                    yield currrange
    def getbinpos(self, item):
        for i in range(int(item / self.binfactor)-1,int(item / self.binfactor)+1):
            if i < len(self.bins) and i >= 0:
                
                for currrange in self.bins[i]:
                    yield currrange
      
            
