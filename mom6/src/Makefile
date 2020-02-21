#
# Makefile for mom6 processing
#
# --- usage: make mom6nc2archv ARCH=intelGF >& Make_mom6nc2archv.log
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
	/bin/rm -f *.o *.a *.mod M*log *2vgrid *2archv *2field *2field3d *2tide *2uvfield

MODS      =	mod_mom6.o mod_xc.o mod_za.o wtime.o
MODM      =	mod_mom6.o mod_xc.o mod_za.o wtime.o mom_eos_wright.o

ARCHV_OBJ = 	bigrid.o blkin.o  \
		forday.o  \
                getdat.o putdat.o \
		zh.o    dum_gks.o 

FIELD_OBJ = 	blkin.o  \
		forday.o  \
                getdat.o putdat.o \
		zh.o    dum_gks.o 

VGRID_OBJ = 	blkin.o  \
		zh.o    dum_gks.o 

#
# --- executables.
#

hycom2vgrid:        $(MODS) hycom2vgrid.o  $(VGRID_OBJ)
	$(LD) $(LDFLAGS) -o hycom2vgrid   hycom2vgrid.o   $(VGRID_OBJ) $(MODS) $(EXTRANCDF)

mom6nc2archv:       $(MODM) mom6nc2archv.o $(ARCHV_OBJ)
	$(LD) $(LDFLAGS) -o mom6nc2archv  mom6nc2archv.o  $(ARCHV_OBJ) $(MODM) $(EXTRANCDF)

mom6nc2tide:        $(MODS) mom6nc2tide.o  $(ARCHV_OBJ)
	$(LD) $(LDFLAGS) -o mom6nc2tide   mom6nc2tide.o   $(ARCHV_OBJ) $(MODS) $(EXTRANCDF)

mom6nc2field:       $(MODS) mom6nc2field.o    $(FIELD_OBJ)
	$(LD) $(LDFLAGS) -o mom6nc2field   mom6nc2field.o   $(FIELD_OBJ) $(MODS) $(EXTRANCDF)

mom6nc2field3d:     $(MODS) mom6nc2field3d.o  $(FIELD_OBJ)
	$(LD) $(LDFLAGS) -o mom6nc2field3d mom6nc2field3d.o $(FIELD_OBJ) $(MODS) $(EXTRANCDF)

mom6nc2uvfield:     $(MODS) mom6nc2uvfield.o  $(FIELD_OBJ)
	$(LD) $(LDFLAGS) -o mom6nc2uvfield mom6nc2uvfield.o $(FIELD_OBJ) $(MODS) $(EXTRANCDF)

mom6nc82field:      $(MODS) mom6nc82field.o  $(FIELD_OBJ)
	$(LD) $(LDFLAGS) -o mom6nc82field  mom6nc82field.o  $(FIELD_OBJ) $(MODS) $(EXTRANCDF)

#
# --- explicit dependencies for each .[fF] file.
#

hycom2vgrid.o:       hycom2vgrid.f       mod_xc.o

mom6nc2archv.o:      mom6nc2archv.f      mod_mom6.o mod_za.o mom_eos_wright.o

mom6nc2tide.o:       mom6nc2tide.f       mod_mom6.o mod_za.o

mom6nc2field.o:      mom6nc2field.f      mod_mom6.o mod_za.o

mom6nc2field3d.o:    mom6nc2field3d.f    mod_mom6.o mod_za.o

mom6nc2uvfield.o:    mom6nc2uvfield.f    mod_mom6.o mod_za.o

mom6nc82field.o:     mom6nc82field.f     mod_mom6.o mod_za.o

bigrid.o:   bigrid.f      mod_mom6.o
blkin.o:    blkin.f
dum_gks.o:  dum_gks.f
forday.o:   forday.f
getdat.o:   getdat.f      mod_mom6.o mod_za.o
putdat.o:   putdat.f      mod_mom6.o mod_za.o
wtime.o:    wtime.F

mod_mom6.o: mod_mom6.F
mod_xc.o:   mod_xc.F   
mod_za.o:   mod_za.F   mod_xc.o

mom_eos_wright.o: mom_eos_wright.f

zh.o:       zh.F zh_sun.f zh_t3e.f zh_ymp.f
