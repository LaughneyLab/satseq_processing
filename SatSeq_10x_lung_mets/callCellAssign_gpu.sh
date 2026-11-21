#!/bin/bash 
infile=$1
outroot=$2
markerFile=$3
seed=$4
modval=$5
n=$6
seed2=$7

echo "USER=$USER"
echo "infile=$infile"
echo "outroot=$outroot"
echo "markerFile=$markerFile"

if [[ ! -e $infile ]]; then echo "$infile does not exist"; exit 1; fi
if [[ ! -e $markerFile ]]; then echo "$markerFile does not exist"; exit 1; fi

R --vanilla --args --infile $infile --outroot $outroot --markers $markerFile --seed $seed --seed2 $seed2 --sample-cells $modval $n  < /workdir/analysis_pipeline/scripts/cellAssign.R > ${outroot}_stdout.txt

