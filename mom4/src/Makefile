#
# Makefile for mom4 processing
#
# --- usage: make mom4nc2archv ARCH=intelGF >& Make_mom4nc2archv.log
#
# --- Tunable parameters in ../../config/$(ARCH)_setup
#
# --- assumes netcdf module and library are in this directory
#

.SUFFIXES: 
.SUFFIXES: .c .F .f .o

.F:
	@echo "Must have an explicit rule for" $*
	@exit 1
.f:
	@echo "Must have an explicit rule for" $*
	@exit 1
.c:
	@echo "Must have an explicit rule for" $*
	@exit 1

include ../../config/$(ARCH)_setup

#
# ---------------------------------------------------------------------
# Standard part of the Makefile starts here
# ---------------------------------------------------------------------
#

all:
	/bin/csh Make_all.csh

clean:
	/bin/rm -f *.o *.a *.mod M*log *2d *3z *restart *hv *mom4 *arctic

MODS     =	mod_mom4.o mod_xc.o mod_za.o wtime.o

ARCHV_OBJ = 	bigrid.o blkin.o  \
		forday.o  \
                getdat.o putdat.o \
		zh.o    dum_gks.o 

#
# --- executables.
#

mom4nc2archv:       $(MODS) mom4nc2archv.o $(ARCHV_OBJ)
	$(LD) $(LDFLAGS) -o mom4nc2archv mom4nc2archv.o $(ARCHV_OBJ) $(MODS) $(EXTRANCDF)

#
# --- explicit dependencies for each .[fF] file.
#

mom4nc2archv.o:      mom4nc2archv.f      mod_mom4.o mod_za.o

bigrid.o:   bigrid.f      mod_mom4.o
blkin.o:    blkin.f
dum_gks.o:  dum_gks.f
forday.o:   forday.f
getdat.o:   getdat.f      mod_mom4.o mod_za.o
putdat.o:   putdat.f      mod_mom4.o mod_za.o
wtime.o:    wtime.F

mod_mom4.o: mod_mom4.F
mod_xc.o:   mod_xc.F   
mod_za.o:   mod_za.F   mod_xc.o

zh.o:       zh.F zh_sun.f zh_t3e.f zh_ymp.f
