#
# Makefile for HYCOM T&S climatology interpolation program
#
# --- usage: make nemo2grid ARCH=intelGF >& Make_nemo2grid.log
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
	/bin/rm -f *.o *.a *.mod M*log *grid

MODS = mod_xc.o mod_za.o wtime.o
LIBS = zh.o blkin.o
LIBN = zh.o blkin.o $(EXTRANCDF)

# select one of the programs to make,
#  e.g. make nemo2grid ARCH=sun >& Make_nemo2grid

nemo2grid:       $(MODS) nemo2grid.o             $(LIBS)
	$(LD) $(LDFLAGS) nemo2grid.o     $(MODS) $(LIBN) -o nemo2grid

mod_xc.o:   mod_xc.F

mod_za.o:   mod_za.F mod_xc.o

zh.o:       zh.F zh_sun.f zh_t3e.f zh_ymp.f

blkin.o:    blkin.f
