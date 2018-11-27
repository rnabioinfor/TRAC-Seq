#!/usr/bin/env bash

#Remove adapters from small RNA sequencing studies
echo "Removing sequencing adapters from reads"
cutadapt -m 15 --adapter='TCGTATGCCGTCTTCT' SRR029131.fastq |  gzip -c  >SRR029131_trimmed.fastq.gz
cutadapt -m 15 --adapter='TCGTATGCCGTCTTCT' SRR029124.fastq |  gzip -c  >SRR029124_trimmed.fastq.gz
cutadapt -m 15 --adapter='CGTATGCCGTCT' SRR207111.fastq |  gzip -c  >SRR207111_trimmed.fastq.gz 
cutadapt -m 15 --adapter='CGTATGCCGTCT' SRR207116.fastq |  gzip -c  >SRR207116_trimmed.fastq.gz


REALNAME=$(readlink -f $0)
SCRIPTDIR=$( cd "$( dirname "$REALNAME" )" && pwd )

#Create the tRNA database
# Params
# 1. tRNA database name
# 2. tRNAscan-SE output file
# 3. Fasta file of reference genome
echo "Creating tRNA database"
"$SCRIPTDIR/maketrnadb.bash" hg19 hg19-tRNAs.out hg19.fa

#Map the tRNAreads
# Params
# 1. Name of experiments
# 2. tRNA database name
# 3. Tab-delimited file specifying the samples and fastq files
# 4. Number of threads for running bowtie2 
echo "Mapping reads to tRNA database"
"$SCRIPTDIR/mapreads.bash" TestTrnas hg19 TrnaSamples.txt hg19-nontRNA-ncRNA.gtf 4
