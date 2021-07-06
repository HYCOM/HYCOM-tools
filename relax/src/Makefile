#
# Makefile for HYCOM T&S climatology interpolation program
#
# --- usage: make z_levitus ARCH=intelGF >& Make_z_levitus.log
#
# --- Tunable parameters in ../../config/$(ARCH)_setup
#

.SUFFIXES: 
.SUFFIXES: .c .F .f .o .f90

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
	/bin/rm -f *.o *.a *.mod *.inc M*log *at *pf *noaa *s *u *u2 *gdem3 *gdem4 *gdem42 *woa *woa13 *density *cer *st *ve *fy *km *field *xi *xv bottom *_linear

MODS = mod_xc.o mod_za.o wtime.o
MODA = $(MODS) mod_ppsw.o
LIBS = blkdat.o interp.o zh.o
LIBN = blkdat.o interp.o zh.o $(EXTRANCDF)

# select one of the programs to make,
#  e.g. make z_levitus ARCH=sun >& Make_z_levitus

tracer_const:    $(MODS) tracer_const.o          $(LIBS)
	$(LD) $(LDFLAGS) tracer_const.o  $(MODS) $(LIBS) -o tracer_const

iso_density:     $(MODS) iso_density.o           $(LIBS)
	$(LD) $(LDFLAGS) iso_density.o   $(MODS) $(LIBS) -o iso_density

z_levitus:       $(MODS) z_levitus.o             $(LIBS)
	$(LD) $(LDFLAGS) z_levitus.o     $(MODS) $(LIBS) -o z_levitus

z_lev_tracer:    $(MODS) z_lev_tracer.o          $(LIBS)
	$(LD) $(LDFLAGS) z_lev_tracer.o  $(MODS) $(LIBS) -o z_lev_tracer

z_caspian:       $(MODS) z_caspian.o             $(LIBS)
	$(LD) $(LDFLAGS) z_caspian.o     $(MODS) $(LIBN) -o z_caspian

z_medatlas:      $(MODS) z_medatlas.o            $(LIBS)
	$(LD) $(LDFLAGS) z_medatlas.o    $(MODS) $(LIBS) -o z_medatlas

z_woa:           $(MODS) z_woa.o                 $(LIBS)
	$(LD) $(LDFLAGS) z_woa.o         $(MODS) $(LIBN) -o z_woa

z_woa_med:       $(MODS) z_woa_med.o             $(LIBS)
	$(LD) $(LDFLAGS) z_woa_med.o     $(MODS) $(LIBN) -o z_woa_med

z_woa_med78:     $(MODS) z_woa_med78.o           $(LIBS)
	$(LD) $(LDFLAGS) z_woa_med78.o   $(MODS) $(LIBN) -o z_woa_med78

z_woa_tracer:    $(MODS) z_woa_tracer.o          $(LIBS)
	$(LD) $(LDFLAGS) z_woa_tracer.o  $(MODS) $(LIBN) -o z_woa_tracer

woa18_to_woa13_format:   woa18_to_woa13_format.o
	$(LD) $(LDFLAGS) woa18_to_woa13_format.o $(EXTRANCDF) -o woa18_to_woa13_format

z_woa13:         $(MODA) z_woa13.o               $(LIBS)
	$(LD) $(LDFLAGS) z_woa13.o       $(MODA) $(LIBN) -o z_woa13

sst_woa:         $(MODS) sst_woa.o               $(LIBS)
	$(LD) $(LDFLAGS) sst_woa.o       $(MODS) $(LIBN) -o sst_woa

z_gdem3:         $(MODS) z_gdem3.o               $(LIBS)
	$(LD) $(LDFLAGS) z_gdem3.o       $(MODS) $(LIBN) -o z_gdem3

sst_gdem3:       $(MODS) sst_gdem3.o             $(LIBS)
	$(LD) $(LDFLAGS) sst_gdem3.o $(MODS)     $(LIBN) -o sst_gdem3

z_gdem4:         $(MODS) z_gdem4.o               $(LIBS)
	$(LD) $(LDFLAGS) z_gdem4.o       $(MODS) $(LIBN) -o z_gdem4

sst_gdem4:       $(MODS) sst_gdem4.o             $(LIBS)
	$(LD) $(LDFLAGS) sst_gdem4.o $(MODS)     $(LIBN) -o sst_gdem4

z_gdem42:        $(MODS) z_gdem42.o              $(LIBS)
	$(LD) $(LDFLAGS) z_gdem42.o      $(MODS) $(LIBN) -o z_gdem42

mom6_gdem42:     $(MODS) mom6_gdem42.o           $(LIBS)
	$(LD) $(LDFLAGS) mom6_gdem42.o   $(MODS) $(LIBN) -o mom6_gdem42

ssh_modas:       $(MODS) ssh_modas.o             $(LIBS)
	$(LD) $(LDFLAGS) ssh_modas.o     $(MODS) $(LIBN) -o ssh_modas

z_modas:         $(MODS) z_modas.o               $(LIBS)
	$(LD) $(LDFLAGS) z_modas.o       $(MODS) $(LIBS) \
	                 -L${OTKHOME}/lib -lsds -ldf -o z_modas

z_const:         $(MODS) z_const.o               $(LIBS)
	$(LD) $(LDFLAGS) z_const.o       $(MODS) $(LIBS) -o z_const

z_modify:        $(MODS) z_modify.o              $(LIBS)
	$(LD) $(LDFLAGS) z_modify.o      $(MODS) $(LIBS) -o z_modify

sst_modas:       $(MODS) sst_modas.o              $(LIBS)
	$(LD) $(LDFLAGS) sst_modas.o      $(MODS) $(LIBN) -o sst_modas

sst_pf:          $(MODS) sst_pf.o                 $(LIBS)
	$(LD) $(LDFLAGS) sst_pf.o         $(MODS) $(LIBS) -o sst_pf

sst_pf_4km:      $(MODS) sst_pf_4km.o             $(LIBS)
	$(LD) $(LDFLAGS) sst_pf_4km.o     $(MODS) $(LIBS) -o sst_pf_4km

sst_noaa:        $(MODS) sst_noaa.o               $(LIBS)
	$(LD) $(LDFLAGS) sst_noaa.o       $(MODS) $(LIBS) -o sst_noaa

sst_rs:          $(MODS) sst_rs.o                 $(LIBS)
	$(LD) $(LDFLAGS) sst_rs.o         $(MODS) $(LIBS) \
	                  -L${OTKHOME}/lib -lsds -ldf -o sst_rs

rmu:             $(MODS) rmu.o          zh.o
	$(LD) $(LDFLAGS) rmu.o  $(MODS) zh.o -o rmu

rmu2:            $(MODS) rmu2.o         zh.o
	$(LD) $(LDFLAGS) rmu2.o $(MODS) zh.o -o rmu2

rmus:            $(MODS) rmus.o         zh.o
	$(LD) $(LDFLAGS) rmus.o $(MODS) zh.o -o rmus

rmu_linear:      $(MODS)   rmu_linear.o                $(LIBS)
	$(LD) $(LDFLAGS)   rmu_linear.o        $(MODS) $(LIBS) -o rmu_linear

zncoda_density:   $(MODS)  zncoda_density.o            $(LIBS)
	$(LD) $(LDFLAGS)   zncoda_density.o    $(MODS) $(LIBS) -o zncoda_density

zncoda_archive:   $(MODS)  zncoda_archive.o            $(LIBS)
	$(LD) $(LDFLAGS)   zncoda_archive.o    $(MODS) $(LIBS) -o zncoda_archive

z_archive:         $(MODS) z_archive.o                 $(LIBS)
	$(LD) $(LDFLAGS)   z_archive.o         $(MODS) $(LIBS) -o z_archive

relax_flat_rivers: $(MODS) relax_flat_rivers.o         $(LIBS)
	$(LD) $(LDFLAGS)   relax_flat_rivers.o $(MODS) $(LIBS) -o relax_flat_rivers

relax_tracer:      $(MODS) relax_tracer.o              $(LIBS)
	$(LD) $(LDFLAGS)   relax_tracer.o      $(MODS) $(LIBS) -o relax_tracer

bottom:          $(MODS) bottom.o                 $(LIBS)
	$(LD) $(LDFLAGS) bottom.o         $(MODS) $(LIBS) -o bottom

relaxv:          $(MODS) relaxv.o                 $(LIBS)
	$(LD) $(LDFLAGS) relaxv.o         $(MODS) $(LIBS) -o relaxv

relaxi:          $(MODS) relaxi.o                 $(LIBS)
	$(LD) $(LDFLAGS) relaxi.o         $(MODS) $(LIBS) -o relaxi

relaxi_dens:     $(MODS) relaxi_dens.o            $(LIBS)
	$(LD) $(LDFLAGS) relaxi_dens.o    $(MODS) $(LIBS) -o relaxi_dens

zonal_field:     $(MODS) zonal_field.o            $(LIBS)
	$(LD) $(LDFLAGS) zonal_field.o    $(MODS) $(LIBS) -o zonal_field

mod_xc.o:   mod_xc.F

mod_za.o:   mod_za.F mod_xc.o

zh.o:       zh.F zh_sun.f zh_t3e.f zh_ymp.f
