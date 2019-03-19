#
# Makefile for ncom processing
#
# --- usage: make ncom2archv ARCH=intelGF >& Make_ncom2archv.log
#
# --- Tunable parameters in ../../config/$(ARCH)_setup
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
	/bin/rm -f *.o *.a *.mod M*log *2d *3z *restart *hv *ncom *arctic

MODS     =	mod_ncom.o mod_xc.o mod_za.o wtime.o

ARCHV_OBJ = 	bigrid.o blkin.o  \
		forday.o  \
                getdat.o putdat.o \
		zh.o    dum_gks.o 

#
# --- executables.
#

ncom2archv:         $(MODS) ncom2archv.o   $(ARCHV_OBJ)
	$(LD) $(LDFLAGS) -o ncom2archv  ncom2archv.o  $(ARCHV_OBJ) $(MODS)

ncomc2archv:        $(MODS) ncomc2archv.o  $(ARCHV_OBJ)
	$(LD) $(LDFLAGS) -o ncomc2archv ncomc2archv.o $(ARCHV_OBJ) $(MODS)

grid2ncom:          $(MODS) grid2ncom.o    $(ARCHV_OBJ)
	$(LD) $(LDFLAGS) -o grid2ncom   grid2ncom.o  $(ARCHV_OBJ) $(MODS)

grid2ncom_arctic:   $(MODS) grid2ncom_arctic.o    $(ARCHV_OBJ)
	$(LD) $(LDFLAGS) -o grid2ncom_arctic   grid2ncom_arctic.o  $(ARCHV_OBJ) $(MODS)

#
# --- explicit dependencies for each .[fF] file.
#

ncom2archv.o:       ncom2archv.f       mod_ncom.o mod_za.o
ncomc2archv.o:      ncomc2archv.f      mod_ncom.o mod_za.o
grid2ncom.o:        grid2ncom.f        mod_ncom.o mod_za.o
grid2ncom_arctic.o: grid2ncom_arctic.f mod_ncom.o mod_za.o

bigrid.o:   bigrid.f      mod_ncom.o
blkin.o:    blkin.f
dum_gks.o:  dum_gks.f
forday.o:   forday.f
getdat.o:   getdat.f      mod_ncom.o mod_za.o
putdat.o:   putdat.f      mod_ncom.o mod_za.o
wtime.o:    wtime.F

mod_ncom.o: mod_ncom.F
mod_xc.o:   mod_xc.F   
mod_za.o:   mod_za.F   mod_xc.o

zh.o:       zh.F zh_sun.f zh_t3e.f zh_ymp.f
