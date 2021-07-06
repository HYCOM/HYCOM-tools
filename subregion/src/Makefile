#
# Makefile for HYCOM full region to subregion archive program
#
# --- usage: make subregion ARCH=intelGF >& Make_subregion.log
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
	/bin/rm -f *.o *.a *.mod *.inc M*log *pog *ion *nam *grid *arche *field *gmap *gmapi *count

MODS = mod_xc.o mod_za.o mod_zb.o wtime.o zh.o

MODX = mod_xc.o mod_za.o mod_zb.o wtime.o zh.o mod_scrip.o

# select one of the programs to make,
#  e.g. make subregion ARCH=sun >& Make_subregion

isubaregion:     $(MODS) isubaregion.o        
	$(LD) $(LDFLAGS) isubaregion.o $(MODS) -o isubaregion

isubregion:      $(MODS) isubregion.o        
	$(LD) $(LDFLAGS) isubregion.o  $(MODS) -o isubregion

isub3region:     $(MODS) isub3region.o extrct.o
	$(LD) $(LDFLAGS) isub3region.o extrct.o $(MODS) -o isub3region

isubp3region:    $(MODS) isubp3region.o extrct.o
	$(LD) $(LDFLAGS) isubp3region.o extrct.o $(MODS) -o isubp3region

isuba_arche:     $(MODS) isuba_arche.o        
	$(LD) $(LDFLAGS) isuba_arche.o $(MODS) -o isuba_arche

isuba_field:     $(MODS) isuba_field.o        
	$(LD) $(LDFLAGS) isuba_field.o $(MODS) -o isuba_field

isubm_field:     $(MODS) isubm_field.o        
	$(LD) $(LDFLAGS) isubm_field.o $(MODS) -o isubm_field

isubs_count:     $(MODX) isubs_count.o        
	$(LD) $(LDFLAGS) isubs_count.o $(MODX) $(EXTRANCDF) -o isubs_count

isubs_field:     $(MODX) isubs_field.o        
	$(LD) $(LDFLAGS) isubs_field.o $(MODX) $(EXTRANCDF) -o isubs_field

isuba_gmap:      $(MODS) isuba_gmap.o        
	$(LD) $(LDFLAGS) isuba_gmap.o  $(MODS) -o isuba_gmap

isuba_gmapi:     $(MODS) isuba_gmapi.o        
	$(LD) $(LDFLAGS) isuba_gmapi.o $(MODS) -o isuba_gmapi

isuba_topog:     $(MODS) isuba_topog.o        
	$(LD) $(LDFLAGS) isuba_topog.o $(MODS) -o isuba_topog

isub_topog:      $(MODS) isub_topog.o        
	$(LD) $(LDFLAGS) isub_topog.o  $(MODS) -o isub_topog

isubx_gmap:      $(MODS) isubx_gmap.o        
	$(LD) $(LDFLAGS) isubx_gmap.o  $(MODS) -o isubx_gmap

sub_grid:        $(MODS) sub_grid.o  extrct.o
	$(LD) $(LDFLAGS) sub_grid.o  extrct.o $(MODS)  -o sub_grid

half_topog:      $(MODS) half_topog.o        
	$(LD) $(LDFLAGS) half_topog.o  $(MODS) -o half_topog

mod_scrip.o: mod_scrip.f

mod_xc.o:   mod_xc.F

mod_za.o:   mod_za.F mod_xc.o

mod_zb.o:   mod_za.F mod_xc.o

zh.o:       zh.F zh_sun.f zh_t3e.f zh_ymp.f
