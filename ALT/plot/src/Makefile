#
# Makefile for hycomproc
#
# --- usage: make hycomproc ARCH=gfortran >& Make_hycomproc.log
#
# --- Tunable parameters in ../../config/$(ARCH)_setup
#
# --- ncargf90 is used in place of $(LD).
# --- make sure this is in your path and consistent with $(LDFLAGS).
# --- if ncargf90 does not exist, create it by editing ncargf77.
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
# output type (one of: metafile, X11, PostScript portrait or landscape)
#
#OPNGKS_OBJ = opngks_meta.o
#OPNGKS_OBJ = opngks_x11.o
#OPNGKS_OBJ = opngks_psl.o
OPNGKS_OBJ = opngks_psp.o

#
# ---------------------------------------------------------------------
# Standard part of the Makefile starts here
# ---------------------------------------------------------------------
#

all:
	/bin/csh Make_all.csh

clean:
	/bin/rm -f *ps? *_x11 *meta *proc *test *cell *ctst *.a *.o *.mod *.inc *.log M*.log


MODS     =	mod_plot.o mod_xc.o mod_za.o wtime.o

PROC_OBJ = 	arrow1.o bigrid.o bigrd1.o blkin.o bord.o \
		colbar.o colors.o conrec.o contur.o \
		extrct.o filtr1.o fram.o fordate.o \
                getdat.o          getdepth.o grdtrns.o \
		horplt.o indxi.o  indxj.o legend1.o \
		         poisnd.o prtmsk.o psmoo.o \
		sbmerg.o xsecij.o ssort.o  zebra.o \
		tracks.o \
		zh.o

PROC_SRC = 	arrow1.f bigrid.f bigrd1.f bklin.f bord.f \
		colbar.f colors.f conrec.f contur.f \
		extrct.f filtr1.f fram.f fordate.f \
		getdat.f          getdepth.f grdtrns.f \
		horplt.f indxi.f  indxj.f legend1.f \
		         poisnd.f prtmsk.f psmoo.f \
		sbmerg.f xsecij.f ssort.f  zebra.f \
		tracks.f \
		mod_xc.F mod_za.F mod_za_sm.F wtime.F

NAVO_OBJ = 	$(PROC_OBJ) getdat_nc.o libnetcdf.a $(EXTRANCDF)

NAVO_LIB = 	libnetcdf.a $(EXTRANCDF)

NAVO_SRC = 	$(PROC_SRC) getdat_nc.f

FIELD_OBJ = 	blkin.o bord.o \
		colbar.o colors.o conrec.o contur.o \
		getdepth.o extrct.o fram.o \
		horplt.o \
		prtmsk.o \
		arrow1.o carrow.o legend1.o \
		tracks.o \
		         zebra.o zh.o
CELL_OBJ = 	blkin.o bord.o \
		colbar.o colors.o conrec.o contur.o \
		getdepth.o extrct.o fram.o \
		prtmsk.o \
		tracks.o \
		cpmpxy_cell.o \
		         zebra.o zh.o

FIELD_SRC = 	fieldproc.f carrow.f fieldcell.f cpmpxy_cell.f

#
# --- default executables.
#
hycomproc: $(MODS) hycomproc.o $(PROC_OBJ) $(OPNGKS_OBJ)
	ncargf90 $(LDFLAGS)  -o hycomproc  hycomproc.o $(MODS) $(PROC_OBJ) $(OPNGKS_OBJ)

hycomnavo: $(MODS) hycomnavo.o $(NAVO_OBJ) $(OPNGKS_OBJ)
	ncargf90 $(LDFLAGS)  -o hycomnavo  hycomnavo.o $(MODS) $(NAVO_OBJ) $(NAVO_LIB) $(OPNGKS_OBJ)

fieldproc: $(MODS) fieldproc.o $(FIELD_OBJ) $(OPNGKS_OBJ)
	ncargf90 $(LDFLAGS)  -o fieldproc  fieldproc.o $(MODS) $(FIELD_OBJ) $(OPNGKS_OBJ)

fieldcell: $(MODS) fieldcell.o $(CELL_OBJ) $(OPNGKS_OBJ)
	ncargf90 $(LDFLAGS)  -o fieldcell  fieldcell.o $(MODS) $(CELL_OBJ) $(OPNGKS_OBJ)

#
# --- executables for specific output media.
#
hp_meta: $(MODS) hycomproc.o $(PROC_OBJ) opngks_meta.o
	ncargf90 $(LDFLAGS)  -o hp_meta  hycomproc.o $(MODS) $(PROC_OBJ) opngks_meta.o

hp_psp: $(MODS) hycomproc.o $(PROC_OBJ) opngks_psp.o
	ncargf90 $(LDFLAGS)  -o hp_psp   hycomproc.o $(MODS) $(PROC_OBJ) opngks_psp.o

hp_psl: $(MODS) hycomproc.o $(PROC_OBJ) opngks_psl.o
	ncargf90 $(LDFLAGS)  -o hp_psl   hycomproc.o $(MODS) $(PROC_OBJ) opngks_psl.o

hp_x11: $(MODS) hycomproc.o $(PROC_OBJ) opngks_x11.o
	ncargf90 $(LDFLAGS)  -o hp_x11   hycomproc.o $(MODS) $(PROC_OBJ) opngks_x11.o

fp_meta: $(MODS) fieldproc.o $(FIELD_OBJ) opngks_meta.o
	ncargf90 $(LDFLAGS)  -o fp_meta    fieldproc.o $(MODS) $(FIELD_OBJ) opngks_meta.o

fp_psp: $(MODS) fieldproc.o $(FIELD_OBJ) opngks_psp.o
	ncargf90 $(LDFLAGS)  -o fp_psp     fieldproc.o $(MODS) $(FIELD_OBJ) opngks_psp.o

fp_psl: $(MODS) fieldproc.o $(FIELD_OBJ) opngks_psl.o
	ncargf90 $(LDFLAGS)  -o fp_psl     fieldproc.o $(MODS) $(FIELD_OBJ) opngks_psl.o

fp_x11: $(MODS) fieldproc.o $(FIELD_OBJ) opngks_x11.o
	ncargf90 $(LDFLAGS)  -o fp_x11     fieldproc.o $(MODS) $(FIELD_OBJ) opngks_x11.o

fc_meta: $(MODS) fieldcell.o $(CELL_OBJ) opngks_meta.o
	ncargf90 $(LDFLAGS)  -o fc_meta    fieldcell.o $(MODS) $(CELL_OBJ) opngks_meta.o

fc_psp:  $(MODS) fieldcell.o $(CELL_OBJ) opngks_psp.o
	ncargf90 $(LDFLAGS)  -o fc_psp     fieldcell.o $(MODS) $(CELL_OBJ) opngks_psp.o

fc_psl:  $(MODS) fieldcell.o $(CELL_OBJ) opngks_psl.o
	ncargf90 $(LDFLAGS)  -o fc_psl     fieldcell.o $(MODS) $(CELL_OBJ) opngks_psl.o

fc_x11:  $(MODS) fieldcell.o $(CELL_OBJ) opngks_x11.o
	ncargf90 $(LDFLAGS)  -o fc_x11     fieldcell.o $(MODS) $(CELL_OBJ) opngks_x11.o

#
# --- test executables.
#
hycomtest: $(MODS) hycomtest.o $(PROC_OBJ) $(OPNGKS_OBJ)
	ncargf90 $(LDFLAGS)  -o hycomtest  hycomtest.o $(MODS) $(PROC_OBJ) $(OPNGKS_OBJ)

fieldtest: $(MODS) fieldtest.o $(FIELD_OBJ) $(OPNGKS_OBJ)
	ncargf90 $(LDFLAGS)  -o fieldtest  fieldtest.o $(MODS) $(FIELD_OBJ) $(OPNGKS_OBJ)

fieldctst: $(MODS) fieldctst.o $(CELL_OBJ) $(OPNGKS_OBJ)
	ncargf90 $(LDFLAGS)  -o fieldctst  fieldctst.o $(MODS) $(CELL_OBJ) $(OPNGKS_OBJ)

#
# --- explicit dependencies for each .[fF] file.
#
hycomnavo.o: hycomnavo.f mod_plot.o mod_za.o
hycomproc.o: hycomproc.f mod_plot.o mod_za.o
fieldproc.o: fieldproc.f mod_plot.o mod_za.o
fieldcell.o: fieldcell.f mod_plot.o mod_za.o

hycomtest.o: hycomtest.f mod_plot.o mod_za.o
fieldtest.o: fieldtest.f mod_plot.o mod_za.o
fieldctst.o: fieldctst.f mod_plot.o mod_za.o

cpmpxy_cell.o: cpmpxy_cell.f mod_plot.o

arrow1.o:   arrow1.f
bigrd1.o:   bigrd1.f      mod_plot.o
bigrid.o:   bigrid.f      mod_plot.o
blkin.o:    blkin.f
bord.o:     bord.f        mod_plot.o
carrow.o:   carrow.f
colbar.o:   colbar.f
colors.o:   colors.f
conrec.o:   conrec.f
contur.o:   contur.f
extrct.o:   extrct.f
filtr1.o:   filtr1.f
fram.o:     fram.f
fordate.o:  fordate.f
getdat.o:   getdat.f      mod_plot.o mod_za.o
getdepth.o: getdepth.f    mod_plot.o mod_za.o
grdtrns.o:  grdtrns.f
horplt.o:   horplt.f
indxi.o:    indxi.f       mod_plot.o
indxj.o:    indxj.f       mod_plot.o
legend1.o:  legend1.f
poisnd.o:   poisnd.f      mod_plot.o
prtmsk.o:   prtmsk.f
psmoo.o:    psmoo.f       mod_plot.o
sbmerg.o:   sbmerg.f      mod_plot.o
ssort.o:    ssort.f
tracks.o:   tracks.f
wtime.o:    wtime.F
xsecij.o:   xsecij.f
zebra.o:    zebra.f

opngks_meta.o: opngks_meta.f
opngks_psp.o:  opngks_psp.f
opngks_psl.o:  opngks_psl.f
opngks_x11.o:  opngks_x11.f

mod_plot.o: mod_plot.F
mod_xc.o:   mod_xc.F   
mod_za.o:   mod_za.F   mod_xc.o

zh.o:       zh.F zh_sun.f zh_t3e.f zh_ymp.f

#
# --- create a .tar file.
#
tar:  *comproc.f $(PROC_SRC) $(ARCHV_SRC) $(FIELD_SRC) opngks*.f Makefile README.*
	tar -cvf hycomproc.tar *comproc.f $(PROC_SRC) $(FIELD_SRC) opngks*.f Makefile README.*
