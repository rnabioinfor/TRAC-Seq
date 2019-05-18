#!/usr/bin/env python

import sys
import subprocess
import argparse
from tempfile import NamedTemporaryFile
import os
from trnasequtils import *

MAXMAPS = 100

THREADS = 1

'''
~/pythonsource/trnaseq/mapreads.py --samplefile=/projects/lowelab/users/holmes/pythonsource/trnatest/testcomp.txt --bowtiedb=hg19-trnagenome >counts.txt

~/pythonsource/trnaseq/mapreads.py --samplefile=agingshort.txt --trnafile=sacCer3-trnatable.txt --bowtiedb=sacCer3-trnagenome >aging-counts.txt

cat <(samtools view -H HBV10_1.bam) <(samtools view HBV10_1.bam | grep UNC11-SN627:301:C21RKACXX:3:1101:20489:57142) | ~/pythonsource/trnaseq/choosemappings.py hg19-trnatable.txt | samtools view -

'''

def wrapbowtie2(bowtiedb, unpaired, outfile, threads, maxmaps = MAXMAPS,program = 'bowtie2', logfile = None):
    bowtiecommand = program+' -x '+bowtiedb+' -k '+str(maxmaps)+' -U '+unpaired+' -p '+str(threads)
    print >>sys.stderr, bowtiecommand
    
    bowtiecommand = bowtiecommand + ' | '+scriptdir+'choosemappings.py '+trnafile+' | samtools sort - '+outfile
    print bowtiecommand 
    bowtierun = None
    if logfile is not None:
        #print >>sys.stderr, "***LOG"
        bowtierun = subprocess.Popen(bowtiecommand, shell = True, stderr = logfile)
    else:                                                       
        bowtierun = subprocess.Popen(bowtiecommand, shell = True)
        
    #bowtierun = subprocess.popen(bowtiecommand, shell = True, stderr = subprocess.PIPE)
    #errorcode = bowtierun.wait()
    #bowtierun.stdout.close()
    errorcode = bowtierun.wait()
    
    
    if errorcode:
        print >>sys.stderr, "Failure to Bowtie2 map"

            
        sys.exit(1)
    
'''
./trnaseq/mapreads.py --samplefile=HumanTrnas.txt --bowtiedb=/scratch/encodeshortrna/hg19
nohup ../trnaseq/mapreads.py --samplefile=SaraShortSamples.txt --bowtiedb=/scratch/encodeshortrna/hg19



nohup ../trnaseq/mapreads.py --samplefile=SaraSamples.txt --bowtiedb=/projects/lowelab/users/holmes/pythonsource/seqqa/combinedb/hg19-pad
'''

parser = argparse.ArgumentParser(description='Map reads with bowtie2 and process mappings')

parser.add_argument('--samplefile',
                   help='Sample file in format')
parser.add_argument('--trnafile',
                   help='tRNA file in format')
parser.add_argument('--logfile',
                   help='optional log file for error messages and mapping stats')
parser.add_argument('--bowtiedb',
                   help='Location of Bowtie2 database')
parser.add_argument('--threads',
                   help='Number of threads for running Bowtie2') 

args = parser.parse_args()

#print >>sys.stderr, os.path.dirname(os.path.realpath(sys.argv[0]))
#print >>sys.stderr, os.path.abspath(__file__)
scriptdir = os.path.dirname(os.path.realpath(sys.argv[0]))+"/"
#sys.exit()

workingdir = './'
#samplefile = open(args.samplefile)
bowtiedb = args.bowtiedb
sampledata = samplefile(args.samplefile)
samples = sampledata.getsamples()

#scriptdir = '/projects/lowelab/users/holmes/pythonsource/trnaseq/'
trnafile = args.trnafile

threads = THREADS
if args.threads:
    threads = args.threads
 
if args.logfile:
    logfile = open(args.logfile,'w')
else:
    logfile = sys.stderr
for samplename in samples:
    print >>logfile, "Mapping "+samplename
    #put stuff in this so that it returns the err if bowtie2 fails instead of logging it
    
    #print >>sys.stderr, sampledata.getfastq(samplename)
    bamfile = workingdir+samplename
    #print >>sys.stderr, bamfile
    #sys.exit()
    
    
    wrapbowtie2(bowtiedb, sampledata.getfastq(samplename),bamfile,threads,logfile=logfile)
    

    
    
    #print >>logfile, "Processing "+samplename +" mappings"

    
    
    #result = subprocess.call(scriptdir+'choosemappings.py '+trnafile+' <'+bamfile +' | samtools view -F 4 -b - | samtools sort - '+workingdir+samplename+'_sort', shell = True)


