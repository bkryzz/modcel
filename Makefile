#********************************#
#                                #
#    Makefile for Modlog v.2.0   #
#                                #
#    F Cleri 02/2017             #
#                                #
#********************************#

SHELL=/bin/sh

.SUFFIXES:
.SUFFIXES: .f90 .o

FC=gfortran
FFLAGS= -O3 -ffree-form
LFLAGS= $(FFLAGS)

EXE=modlog
MOD=modif
MF=Makefile
VPATH=source
RME=README

SRC=    modlog.f setup.f build_mat.f diffusion.f  \
	duplicate.f space.f ssbdsb.f mapping.f  \
	brid3d.f f25b.f f35b.f phase.f transform.f

SRX=    modif.f

_OBJ=   modlog.o setup.o build_mat.o diffusion.o  \
	duplicate.o space.o ssbdsb.o mapping.o  \
	brid3d.o f25.o f35b.o phase.o transform.o

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

void:
	rm fort.* maps*

modlog_s.o: setup.o build_mat.o diffusion.o      \
	  duplicate.o space.o ssbdsb.o mapping.o    \
	  brid3d.o f25.o f35b.o phase.o transform.o



