#!/usr/bin/env python

import re
import sys
import os.path
import itertools
import pysam
import subprocess

from trnasequtils import *


from collections import defaultdict

'''
Here is where I need to use the tRNA ontology between mature tRNAs and chromosomes


'''

def isprimarymapping(mapping):
    return not (mapping.flag & 0x0100 > 0)

def getbesttrnamappings(trnafile, bamout = True, logfile = sys.stderr):
    
    trnadata = transcriptfile(trnafile)
    trnatranscripts = set(trnadata.gettranscripts())
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
    bamfile = pysam.Samfile("-", "r" )
    sort = True
    sortjob = None
    if bamout:
        outfile = pysam.Samfile( "-", "wb", template = bamfile )
    else:
        outfile = pysam.Samfile( "-", "w", template = bamfile )
    for pairedname, allmaps in itertools.groupby(bamfile,lambda x: x.qname):
        allmaps = list(allmaps)
        if sum(curr.flag & 0x004 > 0 for curr in allmaps):
            continue
        totalreads += 1
        #print >>sys.stderr, "**"+pairedname
        readlength = None
        hitscores = dict()
        readtargets = set()
        clipsize = 50
        mappings = 0
        currscore = None
        newset = set()
        readlength = None
        
        #iterate through all mappings of the current read
        #print >>sys.stderr, "**"+str(len(list(allmaps)))
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
            
            shortened += 1
        #print str(len(newset))+"\t"+str(readlength)
        #best mappings are printed out here
        if len(newset) >= 50:
            maxreads += 1
            #print >>sys.stderr, len(newset)
        finalset = list()
        #
        #print >>sys.stderr, ",".join(bamfile.getrname(curr.tid)  for curr in newset)
        #print >>sys.stderr, trnatranscripts
        if sum(bamfile.getrname(curr.tid) in trnatranscripts for curr in newset) > 0:
            
            trnareads += 1
            diff = len(newset) - sum(bamfile.getrname(curr.tid) in trnatranscripts for curr in newset)
            if diff > 0:
                diffreads += 1
            for curr in newset:
                if bamfile.getrname(curr.tid) in trnatranscripts:
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
        mapsremoved += mappings - len(finalset)
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
        

    #print >>logfile, str(diffreads)+"/"+str(trnareads)
    #print >>logfile, str(trnareads)+"/"+str(totalreads)
    #print >>logfile, str(maxreads)+"/"+str(totalreads)
    #print >>logfile, str(multimaps)+"/"+str(totalreads)
    #print >>logfile, str(shortened)+"/"+str(multimaps)
    print >>logfile, "Mappings Removed:"+str(mapsremoved)+"/"+str(totalmaps)
    
    
    

getbesttrnamappings(sys.argv[1])