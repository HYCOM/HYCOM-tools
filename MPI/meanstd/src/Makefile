#
# Makefile for meanstd processing
#
# --- usage: make mean ARCH=xt5 >& Make_mean.log
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
	/bin/rm -f *.o *.a *.mod M*.log *mean *std *diff

MODS     =	mod_mean.o mod_xc.o mod_za.o mod_dimensions.o wtime.o

MEAN_OBJ = 	bigrid.o blkmpi2.o \
                getdat.o putdat.o

#
# --- executables.
#

diff:   $(MODS) hycom_diff.o $(MEAN_OBJ)
	$(LD) $(LDFLAGS)  -o hycom_diff  hycom_diff.o $(MODS) $(MEAN_OBJ)

mean:   $(MODS) hycom_mean.o $(MEAN_OBJ)
	$(LD) $(LDFLAGS)  -o hycom_mean  hycom_mean.o $(MODS) $(MEAN_OBJ)

std:    $(MODS) hycom_std.o $(MEAN_OBJ)
	$(LD) $(LDFLAGS)  -o hycom_std   hycom_std.o  $(MODS) $(MEAN_OBJ)

hycom_diff: diff
hycom_mean: mean
hycom_std:  std

#
# --- explicit dependencies for each .[fF] file.
#

hycom_diff.o: hycom_diff.f $(MODS) $(MEAN_OBJ)
hycom_mean.o: hycom_mean.f $(MODS) $(MEAN_OBJ)
hycom_std.o:  hycom_std.f  $(MODS) $(MEAN_OBJ)

bigrid.o:  bigrid.f      mod_mean.o mod_za.o
blkmpi2.o: blkmpi2.f     mod_xc.o
getdat.o:  getdat.f      mod_mean.o mod_za.o mod_xc.o
putdpt.o:  putdat.f      mod_mean.o mod_za.o mod_xc.o
wtime.o:   wtime.F

mod_dimensions.o: mod_dimensions.F
mod_mean.o: mod_mean.F mod_xc.o
mod_xc.o:   mod_xc.F mod_dimensions.o unit_offset.h
mod_za.o:   mod_za.F mod_xc.o
