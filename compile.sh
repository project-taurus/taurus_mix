#!/bin/bash

#################
#  Description  #
#################
 
# This is an example of script to compile TAURUS_mix. The code requires the
# BLAS/LAPACK libraries. When using the intel compiler "ifort", we recommend
# to use of their specific Math Kernel Library (MKL).

# The script takes one argument: 
#  FC = $1 (Fortran compiler)
#     = "gfortran", "ifort"

# This script is only given as an example and we do not guarantee that it will
# work on your system. In particular, check the version of your compiler and
# the directories of the libraries.

#################
#  Directories  #
#################

heredir=$(pwd)
srcdir=$heredir/src
wrkdir=$heredir/wrk
exedir=$heredir/exe

#############################################
#  Fortran compiler, options and libraries  #
#############################################

FC=$1

# By default, use gfortran
if [ -z $FC ]; then
 FC="gfortran"
fi   

if [ $FC = "ifort" ]; then
  LIB=""
  OPT="-O3 -mkl" 
elif [ $FC = "gfortran" ]; then
  LIB=" "
  LIB="-L/usr/lib -llapack -lblas"
  OPT="-O3" 
  OPT="-O3 -Wall -Wno-maybe-uninitialized"
else
  echo "Wrong compiler ('gfortran' or 'ifort'). Exiting."
  exit
fi

#################
#  Compilation  #
#################

echo "Starting the compilation process with $FC $OPT"
 
code=taurus_mix

# Creates a list of file to be compiled in the correct order
filelist="module_constants.xx module_mathmethods.xx \
          module_parameters.xx module_cutoffs.xx module_projmatelem.xx \
          module_spectroscopy.xx module_initialization.xx"

# The final list of .f90, .f and .o files
filef90=$(echo "$filelist" | sed "s/.xx/.f90/g") 
fileo=$(echo "$filelist" | sed "s/.xx/.o/g" | sed "s/.yy/.o/g") 
filemod=$(echo "$filelist" | sed "s/.xx/.mod/g" | sed "s/module//g" \
                           | sed "s/\_//g")
         
# Creates wrk directory
wrkc=0
if [ ! -d $wrkdir ]; then 
  mkdir $wrkdir; echo "directory '$wrkdir' created"
  wrkc=1
fi

# Copy the files and removes mpi flag if necessary
for file in $filef90
do 
  cp $srcdir/$file $wrkdir/
done

cp $srcdir/${code}.f90 $wrkdir/

echo "source files copied"

# Changes directory and performs the compilation
cd $wrkdir
         
for file in $filef90
do 
  echo "compiling ${file}"
  $FC $OPT -c $file
done

echo "compiling ${code}.f90"
$FC $OPT -o ${code}.exe ${code}.f90 $fileo $LIB

# Creates exe directory and move the exe file
if [ ! -d $exedir ]; then 
  mkdir $exedir; echo "directory '$exedir' created"
fi

if [ -f ${code}.exe ]; then mv ${code}.exe $exedir/; fi

##############
#  Clean up  #
##############

echo "cleaning up" 

cd $heredir

# Removes the wrkdir if not existing prior to the compilation
if [ $wrkc = 1 ]; then
 rm -rf $wrkdir
 echo "directory '$wrkdir' deleted"
else 
  for file in $filef90 $fileo $filemod
  do 
    rm -f $wrkdir/$file
  done
fi

# Final check to see if the exe was produced and move to exedir
if [ -f $exedir/${code}.exe ]; then 
  echo "compilation successful."
else
  echo "compilation failed."
fi
