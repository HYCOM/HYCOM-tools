#
# Makefile for archive processing
#
# --- usage: make hycomproc ARCH=intelGF >& Make_hycomproc.log
#
# --- Tunable parameters in ../../config/$(ARCH)_setup
#

.SUFFIXES: 
.SUFFIXES: .c .F .f .f90 .o

.f90:
	@echo "Must have an explicit rule for" $*
	@exit 1

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
	/bin/rm -f *.o *.a *.mod *.inc M*log *2archv *2d *2data *2b *2t *3z *sf *sfl *sfz *sig2a *restart *hv *ncombc *_inc *_inc2 *gen *ncdf *oam *oam_m *_nc *_temp *_bottom *_mixlay

MODS     =	mod_plot.o mod_xc.o mod_za.o               wtime.o
MODM     =      mod_plot.o mod_xc.o mod_za.o               wtime.o mom_eos_wright.o
MODB     =	mod_plot.o mod_xc.o mod_za.o mod_zb.o      wtime.o
MODZ     =	mod_plot.o mod_xc.o mod_za.o mod_zb.o      wtime.o mod_ppsw.o
MODG     =	mod_plot.o mod_xc.o mod_za.o mod_zb.o      wtime.o mod_ppsw.o
MODR     =	mod_plot.o mod_xc.o mod_za.o mod_restart.o wtime.o
MODN     =	mod_plot.o mod_xc.o mod_za.o mod_ncio.o    wtime.o

DATAX_OBJ = 	bigrid.o           blkin.o    \
		buoflx.o                      \
		extrct.o           fordate.o  \
                getdat.o           getdepth.o \
		          indxi.o  indxj.o    \
                                   psmoo.o    \
		          zebra.o  zh.o       \
                dum_gks.o

DATAH_OBJ = 	$(DATAX_OBJ) horout.o
DATAR_OBJ = 	$(DATAX_OBJ) putdat.o
DATAZ_OBJ = 	$(DATAH_OBJ) layer2z.o mixlay.o

NCDFH_OBJ = 	$(DATAX_OBJ) horout_nc.o
NCDFZ_OBJ = 	$(NCDFH_OBJ) layer2z.o mixlay.o

ARCHV_OBJ = 	bigrid.o blkin.o  \
		extrct.o          fordate.o \
                getdat.o          getdepth.o putdat.o \
		indxi.o  indxj.o             psmoo.o  \
		zebra.o   zh.o    dum_gks.o 

ARCHT_OBJ = 	          zh.o    dum_gks.o 

ARCHG_OBJ = 	bigrid.o blkin.o  \
		extrct.o          fordate.o \
                getdat_nc.o       getdepth.o putdat.o \
		indxi.o  indxj.o             psmoo.o  \
		zebra.o   zh.o    dum_gks.o 

NCDFG_OBJ = 	$(ARCHG_OBJ) horout_nc.o layer2z.o mixlay.o

#
# --- executables.
#

archv2data2d:       $(MODB) archv2data2d.o  $(DATAH_OBJ)
	$(LD) $(LDFLAGS) -o archv2data2d  archv2data2d.o  $(MODB) $(DATAH_OBJ)

archv2data2b:       $(MODB) archv2data2b.o  $(DATAH_OBJ)
	$(LD) $(LDFLAGS) -o archv2data2b  archv2data2b.o  $(MODB) $(DATAH_OBJ)

archv2data2t:       $(MODB) archv2data2t.o  $(DATAH_OBJ)
	$(LD) $(LDFLAGS) -o archv2data2t  archv2data2t.o  $(MODB) $(DATAH_OBJ)

archv2data3z:       $(MODZ) archv2data3z.o  $(DATAZ_OBJ)
	$(LD) $(LDFLAGS) -o archv2data3z  archv2data3z.o  $(MODZ) $(DATAZ_OBJ)

archv2datasf:       $(MODB) archv2datasf.o  $(DATAZ_OBJ)
	$(LD) $(LDFLAGS) -o archv2datasf  archv2datasf.o  $(MODB) $(DATAZ_OBJ)

archv2datasfl:      $(MODB) archv2datasfl.o $(DATAZ_OBJ)
	$(LD) $(LDFLAGS) -o archv2datasfl archv2datasfl.o $(MODB) $(DATAZ_OBJ)

archv2datasfz:      $(MODB) archv2datasfz.o $(DATAZ_OBJ)
	$(LD) $(LDFLAGS) -o archv2datasfz archv2datasfz.o $(MODB) $(DATAZ_OBJ)

field2data:         $(MODB) field2data.o    $(DATAH_OBJ)
	$(LD) $(LDFLAGS) -o field2data    field2data.o    $(MODB) $(DATAH_OBJ)

field2data3z:       $(MODB) field2data3z.o  $(DATAH_OBJ)
	$(LD) $(LDFLAGS) -o field2data3z  field2data3z.o  $(MODB) $(DATAH_OBJ)

field2ncdf:         $(MODB) field2data.o    $(NCDFH_OBJ)
	$(LD) $(LDFLAGS) -o field2ncdf    field2data.o    $(MODB) $(NCDFH_OBJ) $(EXTRANCDF)

field2ncdf3z:       $(MODB) field2data3z.o  $(NCDFH_OBJ)
	$(LD) $(LDFLAGS) -o field2ncdf3z  field2data3z.o  $(MODB) $(NCDFH_OBJ) $(EXTRANCDF)

archv2ncdf2d:       $(MODB) archv2data2d.o  $(NCDFH_OBJ)
	$(LD) $(LDFLAGS) -o archv2ncdf2d  archv2data2d.o  $(MODB) $(NCDFH_OBJ) $(EXTRANCDF)

archv2ncdf2b:       $(MODB) archv2data2b.o  $(NCDFH_OBJ)
	$(LD) $(LDFLAGS) -o archv2ncdf2b  archv2data2b.o  $(MODB) $(NCDFH_OBJ) $(EXTRANCDF)

archv2ncdf2t:       $(MODB) archv2data2t.o  $(NCDFH_OBJ)
	$(LD) $(LDFLAGS) -o archv2ncdf2t  archv2data2t.o  $(MODB) $(NCDFH_OBJ) $(EXTRANCDF)

archv2ncdf3z:       $(MODZ) archv2data3z.o  $(NCDFZ_OBJ)
	$(LD) $(LDFLAGS) -o archv2ncdf3z  archv2data3z.o  $(MODZ) $(NCDFZ_OBJ) $(EXTRANCDF)

archv2ncdfsf:       $(MODB) archv2datasf.o  $(NCDFZ_OBJ)
	$(LD) $(LDFLAGS) -o archv2ncdfsf  archv2datasf.o  $(MODB) $(NCDFZ_OBJ) $(EXTRANCDF)

archv2ncdfsfl:      $(MODB) archv2datasfl.o $(NCDFZ_OBJ)
	$(LD) $(LDFLAGS) -o archv2ncdfsfl archv2datasfl.o $(MODB) $(NCDFZ_OBJ) $(EXTRANCDF)

archv2ncdfsfz:      $(MODB) archv2datasfz.o $(NCDFZ_OBJ)
	$(LD) $(LDFLAGS) -o archv2ncdfsfz archv2datasfz.o $(MODB) $(NCDFZ_OBJ) $(EXTRANCDF)

archv2strmf_nc:     $(MODZ) archv2strmf_nc.o $(NCDFZ_OBJ)
	$(LD) $(LDFLAGS) -o archv2strmf_nc archv2strmf_nc.o $(MODB) $(NCDFZ_OBJ) $(EXTRANCDF)

archv2mom6res:      $(MODN) archv2mom6res.o $(ARCHV_OBJ)
	$(LD) $(LDFLAGS) -o archv2mom6res archv2mom6res.o $(MODN) $(ARCHV_OBJ) $(EXTRANCDF)

archv2ncombc:       $(MODB) archv2ncombc.o  $(DATAH_OBJ)
	$(LD) $(LDFLAGS) -o archv2ncombc  archv2ncombc.o  $(MODB) $(DATAH_OBJ)

archv2oam:          $(MODB) archv2oam.o     $(DATAH_OBJ)
	$(LD) $(LDFLAGS) -o archv2oam     archv2oam.o     $(MODB) $(DATAH_OBJ)

archv2oam_m:        $(MODB) archv2oam_m.o   $(DATAH_OBJ)
	$(LD) $(LDFLAGS) -o archv2oam_m   archv2oam_m.o   $(MODB) $(DATAH_OBJ)

archv2potden:       $(MODS) archv2potden.o $(ARCHV_OBJ)
	$(LD) $(LDFLAGS) -o archv2potden archv2potden.o $(MODS) $(ARCHV_OBJ)

archv2restart:      $(MODR) archv2restart.o $(DATAX_OBJ)
	$(LD) $(LDFLAGS) -o archv2restart archv2restart.o $(MODR) $(DATAX_OBJ)

addbot_archv:       $(MODS) addbot_archv.o  $(ARCHV_OBJ)
	$(LD) $(LDFLAGS) -o addbot_archv addbot_archv.o $(MODS) $(ARCHV_OBJ)

cooperhaines_archv: $(MODS) cooperhaines_archv.o  $(ARCHV_OBJ)
	$(LD) $(LDFLAGS) -o cooperhaines_archv cooperhaines_archv.o $(MODS) $(ARCHV_OBJ)

archt2archv:        $(MODS) archt2archv.o  $(ARCHT_OBJ)
	$(LD) $(LDFLAGS) -o archt2archv archt2archv.o $(MODS) $(ARCHT_OBJ)

hycomarchv:         $(MODS) hycomarchv.o   $(ARCHV_OBJ)
	$(LD) $(LDFLAGS) -o hycomarchv  hycomarchv.o $(MODS) $(ARCHV_OBJ)

botlay_archv:       $(MODS) botlay_archv.o $(ARCHV_OBJ)
	$(LD) $(LDFLAGS) -o botlay_archv  botlay_archv.o $(MODS) $(ARCHV_OBJ)

conv_archv:         $(MODS) conv_archv.o   $(ARCHV_OBJ)
	$(LD) $(LDFLAGS) -o conv_archv  conv_archv.o $(MODS) $(ARCHV_OBJ)

hybgen_archv:       $(MODS) hybgen_archv.o     $(ARCHV_OBJ)
	$(LD) $(LDFLAGS) -o hybgen_archv      hybgen_archv.o     $(MODS) $(ARCHV_OBJ)

hybgen_archv_old:   $(MODS) hybgen_archv_old.o $(ARCHV_OBJ)
	$(LD) $(LDFLAGS) -o hybgen_archv_old  hybgen_archv_old.o $(MODS) $(ARCHV_OBJ)

mrgl_archv:         $(MODS) mrgl_archv.o   $(ARCHV_OBJ)
	$(LD) $(LDFLAGS) -o mrgl_archv  mrgl_archv.o $(MODS) $(ARCHV_OBJ)

ncoda_archv:        $(MODS) ncoda_archv.o      $(ARCHV_OBJ)
	$(LD) $(LDFLAGS) -o ncoda_archv      ncoda_archv.o      $(MODS) $(ARCHV_OBJ)

ncoda_archv_inc:    $(MODS) ncoda_archv_inc.o  $(ARCHV_OBJ)
	$(LD) $(LDFLAGS) -o ncoda_archv_inc  ncoda_archv_inc.o  $(MODS) $(ARCHV_OBJ)

ncoda_archv_inc2:   $(MODS) ncoda_archv_inc2.o $(ARCHV_OBJ)
	$(LD) $(LDFLAGS) -o ncoda_archv_inc2 ncoda_archv_inc2.o $(MODS) $(ARCHV_OBJ)

trim_archv:         $(MODS) trim_archv.o   $(ARCHV_OBJ)
	$(LD) $(LDFLAGS) -o trim_archv  trim_archv.o $(MODS) $(ARCHV_OBJ)

tsis_archv:         $(MODS) tsis_archv.o   $(ARCHV_OBJ)
	$(LD) $(LDFLAGS) -o tsis_archv   tsis_archv.o   $(MODS) $(ARCHV_OBJ)

remap_archv:        $(MODS) remap_archv.o  $(ARCHV_OBJ)
	$(LD) $(LDFLAGS) -o remap_archv  remap_archv.o  $(MODS) $(ARCHV_OBJ)

remaph_archv:       $(MODS) remaph_archv.o $(ARCHV_OBJ)
	$(LD) $(LDFLAGS) -o remaph_archv remaph_archv.o $(MODS) $(ARCHV_OBJ)

remapi_archv:       $(MODS) remapi_archv.o $(ARCHV_OBJ)
	$(LD) $(LDFLAGS) -o remapi_archv remapi_archv.o $(MODS) $(ARCHV_OBJ)

remap_ncombc:       $(MODS) remap_ncombc.o $(ARCHV_OBJ)
	$(LD) $(LDFLAGS) -o remap_ncombc remap_ncombc.o $(MODS) $(ARCHV_OBJ)

restart2archv:      $(MODR) restart2archv.o $(DATAR_OBJ)
	$(LD) $(LDFLAGS) -o restart2archv restart2archv.o $(MODR) $(DATAR_OBJ)

restart2archv_ssh:  $(MODR) restart2archv_ssh.o   $(DATAR_OBJ)
	$(LD) $(LDFLAGS)  -o restart2archv_ssh   restart2archv_ssh.o   $(MODR) $(DATAR_OBJ)

restart2archv_sig2a: $(MODR) restart2archv_sig2a.o $(DATAR_OBJ)
	$(LD) $(LDFLAGS)  -o restart2archv_sig2a restart2archv_sig2a.o $(MODR) $(DATAR_OBJ)

f2d_archv:          $(MODS) f2d_archv.o    $(ARCHV_OBJ)
	$(LD) $(LDFLAGS) -o f2d_archv    f2d_archv.o    $(MODS) $(ARCHV_OBJ)

insitu_archv:       $(MODM) insitu_archv.o $(ARCHV_OBJ)
	$(LD) $(LDFLAGS) -o insitu_archv insitu_archv.o $(MODM) $(ARCHV_OBJ)

potden_archv:       $(MODM) potden_archv.o $(ARCHV_OBJ)
	$(LD) $(LDFLAGS) -o potden_archv potden_archv.o $(MODS) $(ARCHV_OBJ)

steric_archv:       $(MODS) steric_archv.o $(ARCHV_OBJ)
	$(LD) $(LDFLAGS) -o steric_archv steric_archv.o $(MODS) $(ARCHV_OBJ)

stats_archv:        $(MODS) stats_archv.o  $(ARCHV_OBJ)
	$(LD) $(LDFLAGS) -o stats_archv  stats_archv.o  $(MODS) $(ARCHV_OBJ)

gofs_ptemp:         $(MODS) gofs_ptemp.o   $(NCDFG_OBJ)
	$(LD) $(LDFLAGS) -o gofs_ptemp  gofs_ptemp.o  $(MODG) $(NCDFG_OBJ) $(EXTRANCDF)

espc2archv:         $(MODS) espc2archv.o   $(NCDFG_OBJ)
	$(LD) $(LDFLAGS) -o espc2archv espc2archv.o   $(MODG) $(NCDFG_OBJ) $(EXTRANCDF)

espc_bottom:        $(MODS) espc_bottom.o  $(NCDFG_OBJ)
	$(LD) $(LDFLAGS) -o espc_bottom espc_bottom.o $(MODG) $(NCDFG_OBJ) $(EXTRANCDF)

espc_mixlay:        $(MODS) espc_mixlay.o  $(NCDFG_OBJ)
	$(LD) $(LDFLAGS) -o espc_mixlay espc_mixlay.o $(MODG) $(NCDFG_OBJ) $(EXTRANCDF)

#
# --- explicit dependencies for each .[fF] file.
#

archv2data2d.o:  archv2data2d.f  mod_plot.o mod_za.o
archv2data3z.o:  archv2data3z.f  mod_plot.o mod_za.o
archv2datasf.o:  archv2datasf.f  mod_plot.o mod_za.o
archv2datasfl.o: archv2datasfl.f mod_plot.o mod_za.o
archv2datasfz.o: archv2datasfz.f mod_plot.o mod_za.o
archv2ncombc.o:  archv2ncombc.f  mod_plot.o mod_za.o
archv2restart.o: archv2restart.f mod_plot.o mod_za.o mod_restart.o
archv2mom6res.o: archv2mom6res.f mod_plot.o mod_za.o mod_ncio.o

hycomarchv.o:    hycomarchv.f    mod_plot.o mod_za.o
botlay_archv.o:  botlay_archv.f  mod_plot.o mod_za.o
gofs_ptemp.o:    gofs_ptemp.f    mod_plot.o mod_za.o
espc_bottom.o:   espc_bottom.f   mod_plot.o mod_za.o
espc_mixlay.o:   espc_mixlay.f   mod_plot.o mod_za.o
mrgl_archv.o:    mrgl_archv.f    mod_plot.o mod_za.o
trim_archv.o:    trim_archv.f    mod_plot.o mod_za.o
remap_archv.o:   remap_archv.f   mod_plot.o mod_za.o
stats_archv.o:   stats_archv.f   mod_plot.o mod_za.o

hybgen_archv.o:     hybgen_archv.f     mod_plot.o mod_za.o
hybgen_archv_old.o: hybgen_archv_old.f mod_plot.o mod_za.o

ncoda_archv.o:      ncoda_archv.f      mod_plot.o mod_za.o
ncoda_archv_inc.o:  ncoda_archv_inc.f  mod_plot.o mod_za.o
ncoda_archv_inc2.o: ncoda_archv_inc2.f mod_plot.o mod_za.o

insitu_archv.o:     insitu_archv.f     mod_plot.o mod_za.o mom_eos_wright.o
potden_archv.o:     potden_archv.f     mod_plot.o mod_za.o
steric_archv.o:     steric_archv.f     mod_plot.o mod_za.o

restart2archv.o: restart2archv.f mod_plot.o mod_za.o mod_restart.o

bigrid.o:   bigrid.f      mod_plot.o
blkin.o:    blkin.f
buoflx.o:   buoflx.f
dum_gks.o:  dum_gks.f
extrct.o:   extrct.f
fordate.o:  fordate.f
getdat.o:   getdat.f      mod_plot.o mod_za.o
getdepth.o: getdepth.f    mod_plot.o mod_za.o
horout.o:   horout.f      mod_plot.o mod_xc.o mod_zb.o
horout_nc.o: horout_nc.f  mod_plot.o mod_xc.o mod_zb.o
indxi.o:    indxi.f       mod_plot.o
indxj.o:    indxj.f       mod_plot.o
layer2z.o:  layer2z.f
mixlay.o:   mixlay.f
prtmsk.o:   prtmsk.f
psmoo.o:    psmoo.f       mod_plot.o
putdat.o:   putdat.f      mod_plot.o mod_za.o
wtime.o:    wtime.F
zebra.o:    zebra.f

mod_ncio.o: mod_ncio.f90
mod_plot.o: mod_plot.F
mod_ppsw.o: mod_ppsw.F
mod_xc.o:   mod_xc.F   
mod_za.o:   mod_za.F   mod_xc.o
mod_zb.o:   mod_zb.F

mom_eos_wright.o: mom_eos_wright.f

zh.o:       zh.F zh_sun.f zh_t3e.f zh_ymp.f
