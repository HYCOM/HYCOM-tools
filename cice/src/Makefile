#
# Makefile for cice processing
#
# --- usage: make grid2cice ARCH=intelGF >& Make_grid2cice.log
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
	/bin/rm -f *.o *.a *.mod M*log *cice *hycom *tiles

MODS     =	mod_cice.o mod_xc.o mod_za.o wtime.o

GRID_OBJ = 	blkin.o  \
		zh.o    dum_gks.o 


ARCHV_OBJ = 	bigrid.o blkin.o  \
		forday.o  \
                putdat.o \
		zh.o    dum_gks.o 

#
# --- executables.
#

cice2archv:         $(MODS) cice2archv.o   $(ARCHV_OBJ)
	$(LD) $(LDFLAGS) -o cice2archv  cice2archv.o  $(ARCHV_OBJ) $(MODS)

cice2hycom:         $(MODS) cice2hycom.o   $(GRID_OBJ)
	$(LD) $(LDFLAGS) -o cice2hycom  cice2hycom.o  $(GRID_OBJ) $(MODS) $(EXTRANCDF)

cice_tiles:         $(MODS) cice_tiles.o   blkin.o
	$(LD) $(LDFLAGS) -o cice_tiles  cice_tiles.o  blkin.o


grid2cice:          $(MODS) grid2cice.o    $(GRID_OBJ)
	$(LD) $(LDFLAGS) -o grid2cice   grid2cice.o   $(GRID_OBJ)  $(MODS)

#
# --- explicit dependencies for each .[fF] file.
#

cice_tiles.o:       cice_tiles.F
cice2hycom.o:       cice2hycom.f                  mod_za.o
cice2archv.o:       cice2archv.f       mod_cice.o mod_za.o
grid2cice.o:        grid2cice.f        mod_cice.o mod_za.o

bigrid.o:   bigrid.f      mod_cice.o
blkin.o:    blkin.f
dum_gks.o:  dum_gks.f
forday.o:   forday.f
putdat.o:   putdat.f      mod_cice.o mod_za.o
wtime.o:    wtime.F

mod_cice.o: mod_cice.F
mod_xc.o:   mod_xc.F   
mod_za.o:   mod_za.F   mod_xc.o

zh.o:       zh.F zh_sun.f zh_t3e.f zh_ymp.f
