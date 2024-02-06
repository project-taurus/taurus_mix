#################
#  Description  #
#################
 
# This is an example of makefile to compile TAURUS_hwg. The code requires the
# BLAS/LAPACK libraries. When using the intel compiler "ifort", we recommend
# to use of their specific Math Kernel Library (MKL).

# and the directories of the libraries.
# This script is only given as an example and we do not guarantee that it will
# work on your system. In particular, check the version of your compiler and
# the directories of the libraries.

#################
#  Directories  #
#################

SRCDIR=src
OBJNAM=obj
OBJDIR=$(SRCDIR)/${OBJNAM}
MODDIR=$(SRCDIR)/mod
EXEDIR=exe

#############################################
#  Fortran compiler, options and libraries  #
#############################################

# default 
FC=gfortran

ifeq ($(FC),$(filter $(FC), gfortran))
 OPT=-O3 -J$(MODDIR)
 LIB=-L/usr/lib -llapack -lblas 
else ifeq ($(FC),$(filter $(FC), ifort))
 OPT=-O3 -module $(MODDIR) -qmkl -lmkl_sequential
 LIB=
endif

#######################
#  Files definitions  #
#######################

code=taurus_mix
exec=$(code).exe

wSRC77=$(wildcard $(SRCDIR)/*.f)
wSRC90=$(wildcard $(SRCDIR)/*.f90)

SRC77=$(shell echo $(wSRC77) | sed "s/$(SRCDIR)/$(SRCDIR)\/${OBJNAM}/g")
SRC90=$(shell echo $(wSRC90) | sed "s/$(SRCDIR)/$(SRCDIR)\/${OBJNAM}/g")

OBJ90=$(SRC90:.f90=.o)
OBJ77=$(SRC77:.f=.o)

################
#  Make rules  #
################

.PHONY: all clean deepclean

# Main file
#==========

all: $(code)

$(code): $(OBJ90) $(OBJ77) | $(OBJDIR)/ $(MODDIR)/ $(EXEDIR)/
	$(FC) $(OPT) -o $(EXEDIR)/$(exec) $^ $(LIB)
	@echo "compilation finished."

# General rules 
#==============

$(OBJDIR)/%.o: $(SRCDIR)/%.f90 | $(OBJDIR)/ $(MODDIR)/
	@cp $< $(OBJDIR)/tmp.f90
	$(FC) $(OPT) -o $@ -c $(OBJDIR)/tmp.f90
	@rm -f $(OBJDIR)/tmp.f90

$(OBJDIR)/%.o: $(SRCDIR)/%.f | $(OBJDIR)/
	$(FC) $(OPT) -o $@ -c $<

$(OBJDIR)/:
	mkdir -p $(OBJDIR)

$(MODDIR)/:
	mkdir -p $(MODDIR)

$(EXEDIR)/:
	mkdir -p $(EXEDIR)

# Dependencies
#=============

$(OBJDIR)/module_mathmethods.o: $(OBJDIR)/module_constants.o

$(OBJDIR)/module_parameters.o: $(OBJDIR)/module_constants.o 

$(OBJDIR)/module_cutoffs.o: $(OBJDIR)/module_parameters.o

$(OBJDIR)/module_projmatelem.o: $(OBJDIR)/module_cutoffs.o

$(OBJDIR)/module_spectroscopy.o: $(OBJDIR)/module_projmatelem.o $(OBJDIR)/module_mathmethods.o

$(OBJDIR)/module_initialization.o: $(OBJDIR)/module_cutoffs.o

# Debug
#======

debug:
	@echo "SRC90 = $(SRC90)"
	@echo "SRC77 = $(SRC77)"
	@echo "OBJ90 = $(OBJ90)"
	@echo "OBJ77 = $(OBJ77)"
	@echo "code  = $(code)"

# Clean up
#=========

clean:
	rm -f $(MODDIR)/*.mod
	rm -f $(OBJDIR)/*.o  

deepclean:
	rm -rf $(MODDIR)/
	rm -rf $(OBJDIR)/
