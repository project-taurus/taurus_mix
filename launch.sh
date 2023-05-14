#!/bin/bash

# For debugging purposes
#set -x

#################
#  Environment  #
#################

heredir=$(pwd)
exedir=$heredir/exe
wrkdir=$heredir/wrk

code=taurus_mix.exe

input=input.txt

#################
#  Calculation  #
#################

if [ ! -d $wrkdir ]; then mkdir $wrkdir; fi 

cd $wrkdir

cp $exedir/$code .
cp $heredir/$input .        

./$code < $input

# Some cleaning
rm $code $input
