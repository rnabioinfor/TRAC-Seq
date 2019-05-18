#!/usr/bin/env python

import sys
import subprocess
import argparse
from tempfile import NamedTemporaryFile
import os
import threading
from trnasequtils import *

MAXMAPS = 100

'''
~/pythonsource/trnaseq/mapreads.py --samplefile=/projects/lowelab/users/holmes/pythonsource/trnatest/testcomp.txt --bowtiedb=hg19-trnagenome >counts.txt

~/pythonsource/trnaseq/mapreads.py --samplefile=agingshort.txt --trnafile=sacCer3-trnatable.txt --bowtiedb=sacCer3-trnagenome >aging-counts.txt

'''
import Queue
'''
def getbesttrnamappings(inputbam, bamout, trnafile, logfile = sys.stderr):
    
    trnatranscripts = transcriptfile(trnafile)
    alltranscripts = trnatranscripts.gettranscripts()
    readtargets = set()

    hitscores = dict()
    curreadname = None
    readlength = None
    totalmatch = re.compile(r"^(?:(?P<startclip>\d+)S)?(?P<matchlength>\d+)M(?:(?P<endclip>\d+)S)?$")
    #19S17M
    totalclips = 0
    totalreads = 0
    multimaps = 0
    shortened = 0
    mapsremoved = 0
    totalmaps = 0
    #prints the sam header that samtools need to convert to bam
    #grab and iterate through all mappings of a single read
    #most mappers will ensure that read mappings are in the same order as they were in the fastq file, with all mappings of the same read together
    #once the file has been sorted using "samtools sort", this doesn't work anymore.  I can't detect that here, so no error message will be output
    trnareads = 0
    maxreads = 0
    diffreads = 0
    #gzsam = gzip.open(samfile, "rb")
    bamfile = pysam.Samfile(inputbam, "r" )
    sort = True
    sortjob = None
    if bamout:
        outfile = pysam.Samfile( "-", "wb", template = bamfile )
    else:
        outfile = pysam.Samfile( "-", "w", template = bamfile )
    for pairedname, allmaps in itertools.groupby(bamfile,lambda x: x.qname):
        totalreads += 1
        readlength = None
        hitscores = dict()
        readtargets = set()
        clipsize = 50
        mappings = 0
        currscore = None
        newset = set()
        readlength = None
        #iterate through all mappings of the current read
        for currmap in allmaps:
            tagdict = dict()
            for curr in currmap.tags:
                tagdict[curr[0]] = curr[1]
            totalmaps += 1
            if currmap.tid is -1:
                continue
            #print >>sys.stderr, bamfile.getrname(currmap.tid)
            chromname = bamfile.getrname(currmap.tid)
            #sys.exit()
    
            readlength = len(currmap.seq)
            mappings += 1
            #if this is the best mapping of this read, discard any previous mappings
            if currscore is None or currscore < tagdict["AS"]:
                newset = set()
                newset.add(currmap)
                currscore = tagdict["AS"]
            #if this mappings is the same as the previous best mapping, add it to the list
            elif currscore == tagdict["AS"]:
                newset.add(currmap)
            else:
                pass
        #here is where I count a bunch of thing so I can report them at the end
        if mappings > 1:
            multimaps += 1
        if len(newset) < mappings:
            #print  >>sys.stderr, pairedname
            #print >>sys.stderr, str(len(newset))+"/"+str(mappings)
            mapsremoved += mappings - len(newset)
            shortened += 1
        #print str(len(newset))+"\t"+str(readlength)
        #best mappings are printed out here
        if len(newset) >= 50:
            maxreads += 1
            #print >>sys.stderr, len(newset)
        finalset = list()
        if sum(bamfile.getrname(curr.tid) in alltranscripts for curr in newset) > 0:
            trnareads += 1
            diff = len(newset) - sum(bamfile.getrname(curr.tid) in alltranscripts for curr in newset)
            if diff > 0:
                diffreads += 1
            for curr in newset:
                if bamfile.getrname(curr.tid) in alltranscripts:
                    pass
                    #outfile.write(curr)
                    finalset.append(curr)
                else:
                    #print curr.data["bamline"].rstrip()
                    pass
        else:
            for curr in newset:
                pass
                finalset.append(curr)
        #print >>sys.stderr,  sum(isprimarymapping(curr) for curr in finalset)
        #This bit is for ensuring that, if I remove the old primary mapping, a new one is chosen
        #Nesecarry for calculating read proportions
        if sum(isprimarymapping(curr) for curr in finalset) < 1:
            
            
            for i, curr in enumerate(finalset):
                #This
                if i == 0:
                   
                    #print >>sys.stderr, "fixed "+str(len(finalset))+"/"+ str(len(newset))
                    curr.flag &= ~ 0x0100
                    outfile.write(curr)
                else:
                    outfile.write(curr)
        else:
            #print >>sys.stderr, "**"
            for curr in finalset:
                outfile.write(curr)
    outfile.close()        
'''        
        
        
class SamReader(threading.Thread):
    def __init__(self, fifo,queue ):
        #assert callable(fd.readline)
        threading.Thread.__init__(self)
        self._fifo = fifo
        self._queue = queue
    def run(self):
        outputsam = pysam.Samfile(fifoname, "r" )
        for samentry in iter(self._fifo.readline, ''):
            self._queue.put(line)
        
        
 
 
bowtiecommand = ['/projects/lowelab/share/bin/x86_64/bowtie2','-x','sacCer3','-k',"100",'-U',' /projects/lowelab/users/holmes/pythonsource/seqqa/YeastAminoRestrict_sample4_Merge_seqprep.fq.gz']
print >>sys.stderr, " ".join(bowtiecommand)

#tmpdir = tempfile.mkdtemp()
#fifoname = os.path.join(tmpdir, 'samfifo.sam')
fifoname = 'samfifo.sam'
print fifoname
try:
    os.mkfifo(fifoname)
except OSError, e:
    print "Failed to create FIFO: %s" % e
    sys.exit()
    
#getbesttrnamappings(outfile)
#fifoname



sam_queue = Queue.Queue()
bowtiereader = SamReader(fifoname, sam_queue)
bowtiereader.start()

fifo = open(fifoname, "w")

bowtierun = subprocess.Popen(bowtiecommand, stdout = fifo)

print >>sys.stderr, "Started bowtie2"

while not sam_queue.empty():
    bamline = sam_queue.get()
    print bamline.qname


        