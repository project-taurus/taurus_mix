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

# Gatherinr the files 
for file in projmatelem_states.bin projmatelem_M1.bin projmatelem_E2.bin
do 
  if [ -f $file ]; then rm $file; fi
  touch $file
done

for file in matelem*states.bin
do 
  cat $file >> projmatelem_states.bin
done

for file in matelem*M1.bin
do 
  cat $file >> projmatelem_M1.bin
done

for file in matelem*E2.bin
do 
  cat $file >> projmatelem_E2.bin
done

# Running the script      
./$code < input.txt > results
    
echo "The results can be found in the files: results, out/*txt"

# Clean up
rm -f $code input.txt projmatelem.tar.gz *.bin

mv results $heredir/
mv *txt $outdir/
