#!/bin/bash

# For debugging purposes
#set -x

#################
#  Environment  #
#################

heredir=$(pwd)
outdir=$heredir/out
wrkdir=$heredir/wrk
exedir=$heredir/../../exe
auxdir=$heredir/data

code=taurus_mix.exe

if [ ! -d $outdir ]; then mkdir $outdir; fi 
if [ ! -d $wrkdir ]; then mkdir $wrkdir; fi 

#################
#  Calculation  #
#################

cd $wrkdir 

cp $exedir/$code . 
cp $auxdir/template_input.txt input.txt
cp $auxdir/projmatelem.tar.gz .

tar -xf projmatelem.tar.gz

# Running the script      
./$code < input.txt > results
    
echo "The results can be found in the files: results, out/*txt"
echo "They can be compared to the benchmark calculation: diff results data/results_benchmark"

# Clean up
rm -f $code input.txt projmatelem_states.bin projmatelem_E2.bin projmatelem.tar.gz

mv results $heredir/
mv *txt $outdir/
