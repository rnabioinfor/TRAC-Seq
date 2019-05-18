#!/usr/bin/env bash

#$1 is database name
#$2 is trnascan file
#$3 is fasta file of genome

function print_usage() {
  echo "USAGE: $0 databasename tRNAscan.txt genome.fa" >&2
  echo "    databasename: Name of database that will be given to output files " >&2
  echo "    tRNAscan.txt: tRNAscan-SE file containing tRNAs to be used " >&2
  echo "    genome.fa:  Fasta file that contains genome of organism" >&2
 
}

#echo "The script you are running has basename `basename $0`, dirname `dirname $0`"

#exit 1

#`dirname $0`

REALNAME=$(readlink -f $0)
SCRIPTDIR=$( cd "$( dirname "$REALNAME" )" && pwd )

#echo $SCRIPTDIR
#echo ${3} 
	
samtools faidx ${3}
#echo "samtools faidx ${3}"
#exit
"$SCRIPTDIR/getmaturetrnas.py" --trnascan $2  --genome $3  --bedfile=${1}-maturetRNAs.bed --maturetrnatable=${1}-trnatable.txt --trnaalignment=${1}-trnaalign.stk >${1}-maturetRNAs.fa
"$SCRIPTDIR/gettrnabed.py" --trnascan $2  --genome $3  >${1}-trnaloci.bed

"$SCRIPTDIR/aligntrnalocus.py" --genomefile $3 --stkfile=${1}-trnaloci.stk  --trnaloci=${1}-trnaloci.bed

cat ${1}-maturetRNAs.fa $3 >${1}-tRNAgenome.fa
samtools faidx ${1}-tRNAgenome.fa
bowtie2-build --threads 6 ${1}-tRNAgenome.fa ${1}-tRNAgenome
