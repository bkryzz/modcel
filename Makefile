#****************************#
#                            #
# Makefile for Modlog v.1.0  #
#                            #
# F Cleri 04/2014            #
#                            #
#****************************#

SHELL=/bin/sh

.SUFFIXES:
.SUFFIXES: .f90 .o

FC=gfortran
# change the fortran compiler according to your machine
FFLAGS= -O2
LFLAGS= $(FFLAGS)

EXE=modlog
MOD=modif
MF=Makefile
VPATH=source
RME=README

SRC=    modlog.f setup.f build_mat.f      \
	duplicate.f space.f ssbdsb.f mapping.f

SRX=    modif.f

_OBJ=   modlog.o setup.o build_mat.o      \
	duplicate.o space.o ssbdsb.o mapping.o
_OBI=   modif.o
OBJ=   $(patsubst %,$(VPATH)/%,$(_OBJ))
OBI=   $(patsubst %,$(VPATH)/%,$(_OBI))

$(VPATH)/%.o: %.f 
	$(FC) $(FFLAGS) -c -o $@ $<

all:    $(EXE) $(MOD)

$(EXE): $(OBJ)
	$(FC) $(LFLAGS) -o $@ $(OBJ)

$(MOD): $(OBI)
	$(FC) $(LFLAGS) -o $@ $(OBI)

tar:
	tar -cvf $(EXE).tar $(MF) $(patsubst %,$(VPATH)/%,$(SRC)) $(RME)

clean:
	rm -f $(OBJ) $(OBI) $(EXE) $(MOD) fort.* core

modlog.o: setup.o build_mat.o       \
	  duplicate.o space.o ssbdsb.o mapping.o 



