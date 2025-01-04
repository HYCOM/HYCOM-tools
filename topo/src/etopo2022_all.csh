#
set echo
#
# --- Run all scripts needed to download ETOPO 2022 15sec SSH and
# --- convert it to a single netcdf file that can be input to
# --- bathy_15sec, see the example script at the end of this one.
#
# --- Requires the NCO and gdal packages.
# --- A simple way to install gdal is via conda, but note that we
# --- also require the netCDF plugin libgdal-netcdf
# --- https://gdal.org/en/stable/download.html#cross-platform-package-managers
#
# --- Alan J. Wallcraft, COAPS/FSU, January 2025.
#
csh etopo2022_wget.csh
#
csh etopo2022_merge.csh
#
csh etopo2022_rename.csh
#
# --- example script
#
exit
#
cat <<'E-o-D' >! depth_ATLc0.02_22E15s.csh
#! /bin/csh
#
set echo
set time=1
#
# --- interpolate 15-sec ETOPO_2022 to hycom bathymetry.
# --- ETOPO_2022_v1_15s.nc created with gdal_merge.
#
cd ~/hycom/ATLc0.02/topo
#
setenv FOR061  fort.61
setenv FOR061A fort.61A
#
/bin/rm -f $FOR061 $FOR061A
#
setenv CDF_GEBCO /p/work1/wallcraf/topo_ieee/ETOPO_2022_v1_15s/ETOPO_2022_v1_15s.nc
#
~/HYCOM-tools/topo/src/bathy_15sec <<'E-o-D'
 &TOPOG
  CTITLE = 'bathymetery from 15-second ETOPO_2022 global dataset',
  COAST  =     0.02,! DEPTH OF MODEL COASTLINE (-ve keeps orography)
  FLAND  =     0.0, ! FAVOR LAND VALUES
  INTERP = -3,      ! =-N; AVERAGE OVER (2*N+1)x(2*N+1) GRID PATCH
                    ! = 0; PIECEWISE LINEAR.  = 1; CUBIC SPLINE.
  MTYPE  =  0,      ! = 0; CLOSED DOMAIN. = 1; NEAR GLOBAL. = 2; FULLY GLOBAL.
 /
'E-o-D'
#
/bin/mv fort.61  depth_ATLc0.02_22E15s.b
/bin/mv fort.61A depth_ATLc0.02_22E15s.a
'E-o-D'
