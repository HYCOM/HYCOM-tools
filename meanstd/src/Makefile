#
# Makefile for meanstd processing
#
# --- usage: make mnsq ARCH=intelGF >& Make_mnsq.log
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
	/bin/rm -f *.o *.a *.mod M*log *mnsq *std *diff *wsum

MODS     =	mod_mean.o mod_mean_esmf.o mod_xc.o mod_za.o wtime.o zh.o

MEAN_OBJ = 	bigrid.o blkin.o \
		extrct.o \
                getdat.o putdat.o

MEAN_SRC = 	putdat.f

#
# --- executables.
#

hycom_diff: $(MODS) hycom_diff.o $(MEAN_OBJ)
	$(LD) $(LDFLAGS)  -o hycom_diff  hycom_diff.o $(MODS) $(MEAN_OBJ)

hycom_mnsq: $(MODS) hycom_mnsq.o $(MEAN_OBJ)
	$(LD) $(LDFLAGS)  -o hycom_mnsq  hycom_mnsq.o $(MODS) $(MEAN_OBJ)

hycom_std:  $(MODS) hycom_std.o $(MEAN_OBJ)
	$(LD) $(LDFLAGS)  -o hycom_std   hycom_std.o  $(MODS) $(MEAN_OBJ)

hycom_wsum: $(MODS) hycom_wsum.o $(MEAN_OBJ)
	$(LD) $(LDFLAGS)  -o hycom_wsum  hycom_wsum.o $(MODS) $(MEAN_OBJ)

hesmf_mean: $(MODS) hesmf_mean.o $(MEAN_OBJ)
	$(LD) $(LDFLAGS)  -o hesmf_mean  hesmf_mean.o $(MODS) $(MEAN_OBJ)

hesmf_std:  $(MODS) hesmf_std.o $(MEAN_OBJ)
	$(LD) $(LDFLAGS)  -o hesmf_std   hesmf_std.o  $(MODS) $(MEAN_OBJ)

#
# --- explicit dependencies for each .[fF] file.
#

hesmf_mean.o: hesmf_mean.f            mod_mean_esmf.o mod_za.o
hesmf_std.o:  hesmf_std.f             mod_mean_esmf.o mod_za.o

hycom_diff.o: hycom_diff.f mod_mean.o                 mod_za.o
hycom_mnsq.o: hycom_mnsq.f mod_mean.o                 mod_za.o
hycom_std.o:  hycom_std.f  mod_mean.o                 mod_za.o
hycom_wsum.o: hycom_wsum.f mod_mean.o                 mod_za.o

bigrid.o:  bigrid.f        mod_mean.o mod_mean_esmf.o
blkin.o:   blkin.f
extrct.o:  extrct.f
getdat.o:  getdat.f        mod_mean.o mod_mean_esmf.o mod_za.o
putdat.o:  putdat.f        mod_mean.o mod_mean_esmf.o mod_za.o
wtime.o:   wtime.F

mod_mean.o: mod_mean.F
mod_mean_esmf.o: mod_mean_esmf.F
mod_xc.o:   mod_xc.F
mod_za.o:   mod_za.F	mod_xc.o

zh.o:       zh.F zh_sun.f zh_t3e.f zh_ymp.f
