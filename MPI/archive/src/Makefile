#
# Makefile for archive processing
#
# --- usage: make hycomproc ARCH=sun >& Make_hycomproc
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
#include ../../config/$(ARCH)-mpi2io_mpi

#
# ---------------------------------------------------------------------
# Standard part of the Makefile starts here
# ---------------------------------------------------------------------
#

all:	archv2data2d archv2data2t archv2data3z archv2strmf archv2ncdf2d archv2ncdf2t archv2ncdf3z archv2restart ncoda_archv ncoda_archv_inc ncoda_archv_inc2 ncoda_archv_new restart2archv

clean:
	/bin/rm -f *.o *.a *.mod M*log *strmf *2archv *2d *2data *2t *3z *sf *sfl *sfz *sig2a *restart *hv *ncombc *_inc *_inc2 *_new *9t *7t


MODS     =	mod_plot.o mod_xc.o mod_dimensions.o mod_za.o wtime.o
MODZ     =	mod_plot.o mod_xc.o mod_dimensions.o mod_za.o wtime.o mod_ppsw.o
MODR     =	mod_plot.o mod_xc.o mod_dimensions.o mod_za.o wtime.o mod_restart.o

DATAX_OBJ = 	bigrid.o           blkmpi2.o     buoflx.o  fordate.o   getdat.o  getdepth.o \
                              psmoo.o        zh_sun.o 

DATAH_OBJ = 	$(DATAX_OBJ) horout.o
DATAR_OBJ = 	$(DATAX_OBJ) putdat.o
DATAZ_OBJ = 	$(DATAH_OBJ) layer2z.o mixlay.o

NCDFH_OBJ = 	$(DATAX_OBJ) horout_nc.o
NCDFZ_OBJ = 	$(NCDFH_OBJ) layer2z.o mixlay.o

#
# --- executables.
#
  

archv2data2d:       $(MODZ) archv2data2d.o  $(DATAZ_OBJ)
	$(LD) $(LDFLAGS) -o archv2data2d  archv2data2d.o  $(MODZ) $(DATAZ_OBJ)

archv2data2t:       $(MODZ) archv2data2t.o  $(DATAZ_OBJ)
	$(LD) $(LDFLAGS) -o archv2data2t  archv2data2t.o  $(MODZ) $(DATAZ_OBJ)

archv2data3z:       $(MODZ) archv2data3z.o  $(DATAZ_OBJ)
	$(LD) $(LDFLAGS) -o archv2data3z  archv2data3z.o  $(MODZ) $(DATAZ_OBJ)

archv2ncdf2d:       $(MODZ) archv2data2d.o  $(NCDFZ_OBJ)
	$(LD) $(LDFLAGS) -o archv2ncdf2d  archv2data2d.o  $(MODZ) $(NCDFZ_OBJ) $(EXTRANCDF)

archv2ncdf2t:       $(MODZ) archv2data2t.o  $(NCDFZ_OBJ)
	$(LD) $(LDFLAGS) -o archv2ncdf2t  archv2data2t.o  $(MODZ) $(NCDFZ_OBJ) $(EXTRANCDF)

archv2ncdf3z:       $(MODZ) archv2data3z.o  $(NCDFZ_OBJ)
	$(LD) $(LDFLAGS) -o archv2ncdf3z  archv2data3z.o  $(MODZ) $(NCDFZ_OBJ) $(EXTRANCDF)

archv2strmf:        $(MODZ) archv2strmf.o   $(NCDFZ_OBJ)
	$(LD) $(LDFLAGS) -o archv2strmf   archv2strmf.o   $(MODZ) $(NCDFZ_OBJ) $(EXTRANCDF)

archv2restart:      $(MODR) archv2restart.o $(DATAX_OBJ)
	$(LD) $(LDFLAGS) -o archv2restart archv2restart.o $(MODR) $(DATAX_OBJ)

ncoda_archv:        $(MODS) ncoda_archv.o $(DATAR_OBJ)
	$(LD) $(LDFLAGS) -o ncoda_archv ncoda_archv.o $(MODS) $(DATAR_OBJ)

ncoda_archv_inc:    $(MODS) ncoda_archv_inc.o  $(DATAR_OBJ)
	$(LD) $(LDFLAGS) -o ncoda_archv_inc ncoda_archv_inc.o $(MODS) $(DATAR_OBJ)

ncoda_archv_inc2:   $(MODS) ncoda_archv_inc2.o $(DATAR_OBJ)
	$(LD) $(LDFLAGS) -o ncoda_archv_inc2 ncoda_archv_inc2.o $(MODS) $(DATAR_OBJ)

ncoda_archv_new:    $(MODS) ncoda_archv_new.o  $(DATAR_OBJ)
	$(LD) $(LDFLAGS) -o ncoda_archv_new ncoda_archv_new.o $(MODS) $(DATAR_OBJ)

restart2archv:      $(MODR) restart2archv.o $(DATAR_OBJ)
	$(LD) $(LDFLAGS) -o restart2archv restart2archv.o $(MODR) $(DATAR_OBJ)

#
# --- explicit dependencies for each .[fF] file.
#

archv2data2d.o:  archv2data2d.f  mod_plot.o mod_za.o
archv2data2t.o:  archv2data2t.f  mod_plot.o mod_za.o
archv2data3z.o:  archv2data3z.f  mod_plot.o mod_za.o

archv2restart.o:    archv2restart.f    mod_plot.o mod_za.o mod_restart.o
ncoda_archv.o:      ncoda_archv.f      mod_plot.o mod_za.o
ncoda_archv_inc.o:  ncoda_archv_inc.f  mod_plot.o mod_za.o
ncoda_archv_inc2.o: ncoda_archv_inc2.f mod_plot.o mod_za.o
ncoda_archv_new.o:  ncoda_archv_new.f  mod_plot.o mod_za.o

bigrid.o:   bigrid.f      mod_plot.o
blkmpi2.o:  blkmpi2.f
buoflx.o:   buoflx.f
fordate.o:  fordate.f
getdat.o:   getdat.f      mod_plot.o mod_za.o
getdepth.o: getdepth.f    mod_plot.o mod_za.o
horout.o:   horout.f      mod_plot.o mod_xc.o
horout_nc.o:   horout_nc.f      mod_plot.o mod_xc.o
indxj.o:     indxj.f
layer2z.o:  layer2z.f
mixlay.o:   mixlay.f
psmoo.o:    psmoo.f       mod_plot.o
wtime.o:    wtime.F
putdat.o:   putdat.f      mod_plot.o mod_za.o

mod_xc.o:   mod_xc.F mod_dimensions.o  
mod_za.o:   mod_za.F   mod_xc.o
mod_dimensions.o:  mod_dimensions.f

mod_plot.o: mod_plot.f mod_xc.o
mod_ppsw.o: mod_ppsw.f mod_xc.o
mod_restart.o:  mod_restart.f mod_plot.o mod_za.o mod_xc.o
zh_sun.o:       zh_sun.f
