#
# Makefile for bathymetery processing (e.g. MICOM to HYCOM conversion).
#
# --- usage: make batrop ARCH=intelGF >& Make_batrop.log
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
	/bin/rm -f *.o *.a *.mod netcdf.inc M*log *hgram *min *1d *2d *2h *ppmX *sea *oth *clip *ff *flat *fill *mask *set *shrink *sub *map *nam *sby *advect *trop *tor *rge *fy *it *noshrink *son *ude *zcells *zthin *arctic *_test *_ij *_resize *_ns *1deg  *lat *ts *ed *pe *gh *ip *15sec *30sec *bedmap *RTopo *sindhu *pi *sea-b *360 *island *vrtcmprs *lon *find *mom6 *dx *tiles *ic *1st

MODS = mod_xc.o mod_za.o wtime.o
MODN = mod_xc.o mod_za.o wtime.o 
MODG = $(MODS) spherdist.o rotang.o
MOD6 = $(MODN) spherdist.o rotang.o geodesic.o

# update ARCH in Make_all.csh and 
# make everything with:  csh Make_all.csh
#
# or select one of the programs to make,
#  e.g. make m2h ARCH=sun >& Make_m2h

bathy_15sec:     $(MODS) bathy_15sec.o   interp.o zh.o
	$(LD) $(LDFLAGS) bathy_15sec.o   interp.o zh.o $(MODN) $(EXTRANCDF) -o bathy_15sec

bathy_RTopo:     $(MODS) bathy_RTopo.o   interp.o zh.o
	$(LD) $(LDFLAGS) bathy_RTopo.o   interp.o zh.o $(MODN) $(EXTRANCDF) -o bathy_RTopo

bathy_30sec:     $(MODS) bathy_30sec.o   interp.o zh.o
	$(LD) $(LDFLAGS) bathy_30sec.o   interp.o zh.o $(MODN) $(EXTRANCDF) -o bathy_30sec

bathy_01min:     $(MODS) bathy_01min.o   interp.o zh.o
	$(LD) $(LDFLAGS) bathy_01min.o   interp.o zh.o $(MODN) $(EXTRANCDF) -o bathy_01min

bathy_02min:     $(MODS) bathy_02min.o   interp.o zh.o
	$(LD) $(LDFLAGS) bathy_02min.o   interp.o zh.o $(MODS) -o bathy_02min

bathy_05min:     $(MODS) bathy_05min.o   interp.o zh.o
	$(LD) $(LDFLAGS) bathy_05min.o   interp.o zh.o $(MODS) -o bathy_05min

bathy_bedmap:    $(MODS) bathy_bedmap.o  interp.o zh.o
	$(LD) $(LDFLAGS) bathy_bedmap.o  interp.o zh.o $(MODN) $(EXTRANCDF) -o bathy_bedmap

bathy_caspi:     $(MODS) bathy_caspi.o   interp.o zh.o
	$(LD) $(LDFLAGS) bathy_caspi.o   interp.o zh.o $(MODS) -o bathy_caspi

bathy_sindhu:    $(MODS) bathy_sindhu.o  interp.o zh.o
	$(LD) $(LDFLAGS) bathy_sindhu.o  interp.o zh.o $(MODN) $(EXTRANCDF) -o bathy_sindhu


landsea_15sec:   $(MODS) landsea_15sec.o interp.o zh.o
	$(LD) $(LDFLAGS) landsea_15sec.o interp.o zh.o $(MODN) $(EXTRANCDF) -o landsea_15sec

landsea_30sec:   $(MODS) landsea_30sec.o interp.o zh.o
	$(LD) $(LDFLAGS) landsea_30sec.o interp.o zh.o $(MODN) $(EXTRANCDF) -o landsea_30sec

landsea_01min:   $(MODS) landsea_01min.o interp.o zh.o
	$(LD) $(LDFLAGS) landsea_01min.o interp.o zh.o $(MODN) $(EXTRANCDF) -o landsea_01min

landsea_02min:   $(MODS) landsea_02min.o interp.o zh.o
	$(LD) $(LDFLAGS) landsea_02min.o interp.o zh.o $(MODS) -o landsea_02min

landsea_05min:   $(MODS) landsea_05min.o interp.o zh.o
	$(LD) $(LDFLAGS) landsea_05min.o interp.o zh.o $(MODS) -o landsea_05min

landsea_caspi:   $(MODS) landsea_caspi.o interp.o zh.o
	$(LD) $(LDFLAGS) landsea_caspi.o interp.o zh.o $(MODS) -o landsea_caspi

1d:              $(MODS) topo_1d.o         zh.o
	$(LD) $(LDFLAGS) topo_1d.o         zh.o $(MODS) -o topo_1d

2d:              $(MODS) topo_2d.o         zh.o
	$(LD) $(LDFLAGS) topo_2d.o         zh.o $(MODS) -o topo_2d

advect:          $(MODG) topo_advect.o     zh.o
	$(LD) $(LDFLAGS) topo_advect.o     zh.o $(MODS) -o topo_advect

biharmonic:      $(MODG) topo_biharmonic.o zh.o
	$(LD) $(LDFLAGS) topo_biharmonic.o zh.o $(MODS) -o topo_biharmonic

batrop:          $(MODG) topo_batrop.o     zh.o
	$(LD) $(LDFLAGS) topo_batrop.o     zh.o $(MODS) -o topo_batrop

clip:            $(MODS) topo_clip.o       zh.o
	$(LD) $(LDFLAGS) topo_clip.o       zh.o $(MODS) -o topo_clip

diff:            $(MODS) topo_diff.o       zh.o
	$(LD) $(LDFLAGS) topo_diff.o       zh.o $(MODS) -o topo_diff

edit:            $(MODS) topo_edit.o       zh.o
	$(LD) $(LDFLAGS) topo_edit.o       zh.o $(MODS) -o topo_edit

flat:            $(MODS) topo_flat.o       zh.o
	$(LD) $(LDFLAGS) topo_flat.o       zh.o $(MODS) -o topo_flat

grid_360:        $(MODS) grid_360.o        zh.o
	$(LD) $(LDFLAGS) grid_360.o        zh.o $(MODS) -o grid_360

shrink:          $(MODG) topo_shrink.o     zh.o
	$(LD) $(LDFLAGS) topo_shrink.o     zh.o $(MODS) -o topo_shrink

slope:           $(MODG) topo_slope.o      zh.o
	$(LD) $(LDFLAGS) topo_slope.o      zh.o $(MODS) -o topo_slope

ellipsoid:       $(MODG) grid_ellipsoid.o  zh.o
	$(LD) $(LDFLAGS) grid_ellipsoid.o  zh.o $(MODG) -o grid_ellipsoid

hudson:          $(MODG) grid_hudson.o     zh.o
	$(LD) $(LDFLAGS) grid_hudson.o     zh.o $(MODG) -o grid_hudson

island:          $(MODS) topo_island.o     zh.o
	$(LD) $(LDFLAGS) topo_island.o     zh.o $(MODS) -o topo_island

landfill:        $(MODS) topo_landfill.o   zh.o
	$(LD) $(LDFLAGS) topo_landfill.o   zh.o $(MODS) -o topo_landfill

landmask:        $(MODS) topo_landmask.o   zh.o
	$(LD) $(LDFLAGS) topo_landmask.o   zh.o $(MODS) -o topo_landmask

latitude:        $(MODG) grid_latitude.o   zh.o
	$(LD) $(LDFLAGS) grid_latitude.o   zh.o $(MODG) -o grid_latitude

lonlat:          $(MODG) grid_lonlat.o     zh.o
	$(LD) $(LDFLAGS) grid_lonlat.o     zh.o $(MODG) -o grid_lonlat

lonlat_2d:       $(MODG) grid_lonlat_2d.o  zh.o
	$(LD) $(LDFLAGS) grid_lonlat_2d.o  zh.o $(MODG) -o grid_lonlat_2d

lpanam:          $(MODG) grid_lpanam.o     zh.o
	$(LD) $(LDFLAGS) grid_lpanam.o     zh.o $(MODG) -o grid_lpanam

map:             $(MODS) topo_map.o        zh.o
	$(LD) $(LDFLAGS) topo_map.o        zh.o $(MODS) -o topo_map

mapsub:          $(MODS) topo_mapsub.o     zh.o
	$(LD) $(LDFLAGS) topo_mapsub.o     zh.o $(MODS) -o topo_mapsub

mercator:        $(MODG) grid_mercator.o   zh.o
	$(LD) $(LDFLAGS) grid_mercator.o   zh.o $(MODG) -o grid_mercator

mom6:        $(MOD6)     grid_mom6.o       zh.o
	$(LD) $(LDFLAGS) grid_mom6.o       zh.o $(MOD6) $(EXTRANCDF) -o grid_mom6

2mom6:       $(MOD6)     topo_2mom6.o      zh.o
	$(LD) $(LDFLAGS) topo_2mom6.o      zh.o $(MOD6) $(EXTRANCDF) -o topo_2mom6

hycom2mom6:  $(MOD6)     grid_hycom2mom6.o   zh.o
	$(LD) $(LDFLAGS) grid_hycom2mom6.o   zh.o $(MOD6) $(EXTRANCDF) -o grid_hycom2mom6

hycom2mom6dx: $(MOD6)    grid_hycom2mom6dx.o zh.o
	$(LD) $(LDFLAGS) grid_hycom2mom6dx.o zh.o $(MOD6) $(EXTRANCDF) -o grid_hycom2mom6dx

modify:          $(MODS) mask_modify.o     zh.o
	$(LD) $(LDFLAGS) mask_modify.o     zh.o $(MODS) -o mask_modify

panam:           $(MODG) grid_panam.o      zh.o
	$(LD) $(LDFLAGS) grid_panam.o      zh.o $(MODG) -o grid_panam

ports:           $(MODS) topo_ports.o      zh.o
	$(LD) $(LDFLAGS) topo_ports.o      zh.o $(MODS) -o topo_ports

ports_find:      $(MODS) topo_ports_find.o zh.o
	$(LD) $(LDFLAGS) topo_ports_find.o zh.o $(MODS) -o topo_ports_find

ports_latlon:    $(MODS) topo_ports_latlon.o zh.o
	$(LD) $(LDFLAGS) topo_ports_latlon.o zh.o $(MODS) -o topo_ports_latlon

olddist:         $(MODS) grid_olddist.o    zh.o
	$(LD) $(LDFLAGS) grid_olddist.o    zh.o $(MODS) -o grid_olddist

rosby:           $(MODS) grid_rosby.o      zh.o
	$(LD) $(LDFLAGS) grid_rosby.o      zh.o $(MODS) -o grid_rosby

rough:           $(MODG) topo_rough.o      zh.o
	$(LD) $(LDFLAGS) topo_rough.o      zh.o $(MODS) -o topo_rough

rotated:         $(MODG) grid_rotated.o    zh.o
	$(LD) $(LDFLAGS) grid_rotated.o    zh.o $(MODG) -o grid_rotated

merge:           $(MODS) topo_merge.o      zh.o
	$(LD) $(LDFLAGS) topo_merge.o      zh.o $(MODS) -o topo_merge

onesea:          $(MODS) topo_onesea.o     zh.o
	$(LD) $(LDFLAGS) topo_onesea.o     zh.o $(MODS) -o topo_onesea

onesea-b:        $(MODS) topo_onesea-b.o   zh.o
	$(LD) $(LDFLAGS) topo_onesea-b.o   zh.o $(MODS) -o topo_onesea-b

onesea_fill:     $(MODS) topo_onesea_fill.o zh.o
	$(LD) $(LDFLAGS) topo_onesea_fill.o zh.o $(MODS) -o topo_onesea_fill

onesea-b_fill:   $(MODS) topo_onesea-b_fill.o zh.o
	$(LD) $(LDFLAGS) topo_onesea-b_fill.o zh.o $(MODS) -o topo_onesea-b_fill

onesea-b_mask:   $(MODS) topo_onesea-b_mask.o zh.o
	$(LD) $(LDFLAGS) topo_onesea-b_mask.o zh.o $(MODS) -o topo_onesea-b_mask

partit:          $(MODS) partit.o          zh.o
	$(LD) $(LDFLAGS) partit.o          zh.o $(MODS) -o partit

partit_arctic:   $(MODS) partit_arctic.o   zh.o
	$(LD) $(LDFLAGS) partit_arctic.o   zh.o $(MODS) -o partit_arctic

partit_arctic_ns: $(MODS) partit_arctic_ns.o zh.o
	$(LD) $(LDFLAGS)  partit_arctic_ns.o zh.o $(MODS) -o partit_arctic_ns

partit_noshrink: $(MODS) partit_noshrink.o zh.o
	$(LD) $(LDFLAGS) partit_noshrink.o zh.o $(MODS) -o partit_noshrink

partit_landsea:  $(MODS) partit_landsea.o  zh.o
	$(LD) $(LDFLAGS) partit_landsea.o  zh.o $(MODS) -o partit_landsea

partit_mom6:     $(MODS) partit_mom6.o zh.o
	$(LD) $(LDFLAGS) partit_mom6.o     zh.o $(MODS) -o partit_mom6

partit_resize:   $(MODS) partit_resize.o   zh.o
	$(LD) $(LDFLAGS) partit_resize.o   zh.o $(MODS) -o partit_resize

partit_test:     $(MODS) partit_test.o     zh.o
	$(LD) $(LDFLAGS) partit_test.o     zh.o $(MODS) -o partit_test

partit_ij:       $(MODS) partit_ij.o       zh.o
	$(LD) $(LDFLAGS) partit_ij.o       zh.o $(MODS) -o partit_ij

partit_row1st:   $(MODS) partit_row1st.o   zh.o
	$(LD) $(LDFLAGS) partit_row1st.o   zh.o $(MODS) -o partit_row1st

hgram:           $(MODS) topo_hgram.o      zh.o
	$(LD) $(LDFLAGS) topo_hgram.o      zh.o $(MODS) -o topo_hgram

ppmX:            $(MODS) topo_ppmX.o       zh.o
	$(LD) $(LDFLAGS) topo_ppmX.o       zh.o $(MODS) -o topo_ppmX

diff_ppmX:       $(MODS) topo_diff_ppmX.o  zh.o
	$(LD) $(LDFLAGS) topo_diff_ppmX.o  zh.o $(MODS) -o topo_diff_ppmX

smallsea:        $(MODS) topo_smallsea.o   zh.o
	$(LD) $(LDFLAGS) topo_smallsea.o   zh.o $(MODS) -o topo_smallsea

smooth:          $(MODS) topo_smooth.o     zh.o
	$(LD) $(LDFLAGS) topo_smooth.o     zh.o $(MODS) -o topo_smooth

smooth_skip:     $(MODS) topo_smooth_skip.o  zh.o
	$(LD) $(LDFLAGS) topo_smooth_skip.o  zh.o $(MODS) -o topo_smooth_skip

smooth_slope:    $(MODS) topo_smooth_slope.o zh.o
	$(LD) $(LDFLAGS) topo_smooth_slope.o zh.o $(MODS) -o topo_smooth_slope

subset:          $(MODS) topo_subset.o     zh.o
	$(LD) $(LDFLAGS) topo_subset.o     zh.o $(MODS) -o topo_subset

tiles:           $(MODS) topo_tiles.o      zh.o
	$(LD) $(LDFLAGS) topo_tiles.o      zh.o $(MODS) -o topo_tiles

vrtcmprs:        $(MODS) topo_vrtcmprs.o   zh.o
	$(LD) $(LDFLAGS) topo_vrtcmprs.o   zh.o $(MODS) -o topo_vrtcmprs

zcells:          $(MODS) topo_zcells.o     zh.o
	$(LD) $(LDFLAGS) topo_zcells.o     zh.o $(MODS) -o topo_zcells

zthin:           $(MODS) topo_zthin.o      zh.o
	$(LD) $(LDFLAGS) topo_zthin.o      zh.o $(MODS) -o topo_zthin

latlonmask:      $(MODN) espc_latlonmask.o zh.o
	$(LD) $(LDFLAGS) espc_latlonmask.o zh.o $(MODN) $(EXTRANCDF) -o espc_latlonmask

mod_xc.o:   mod_xc.F

mod_za.o:   mod_za.F mod_xc.o

zh.o:       zh.F zh_sun.f zh_t3e.f zh_ymp.f
