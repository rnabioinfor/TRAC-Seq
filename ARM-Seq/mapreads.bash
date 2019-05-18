#!/usr/bin/env bash

#$1 is experiment name
#$2 is database name
#$3 is sample file
#$4 is bed feature for other sRNAs
#$5 is number of threads for running bowtie2

function print_usage() {
  echo "USAGE: $0 experimentname databasename samplefile.txt otherfeatures.bed threads" >&2
  echo "    experimentname: Name of experiment that will be given to output files " >&2
  echo "    databasename: Name of database created by maketrnadb.bash " >&2
  echo "    samplefile.txt: tRNAscan-SE file containing tRNAs to be used " >&2
  echo "    otherfeatures.bed:  Bed file containing non-tRNA features" >&2
  echo "    threads: Number of threads for running bowtie2" >&2
}



REALNAME=$(readlink -f $0)
SCRIPTDIR=$( cd "$( dirname "$REALNAME" )" && pwd )


#"$SCRIPTDIR/mapreads.py" --samplefile=$3 --trnafile=$2-trnatable.txt --bowtiedb=${2}-tRNAgenome --threads=$5 

"$SCRIPTDIR/countreads.py" --samplefile=$3 --bedfile=$4 --maturetrnas=$2-maturetRNAs.bed --trnaloci=${2}-trnaloci.bed >$1-counts.txt
