#
set echo
#
# --- merge ETOPO 2022 ssh 15s tiles into a global array
# --- this takes 5 hours on a fast system and a unknown amount of memory
#
date
/bin/rm -f    ETOPO_2022_v1_15s.nc
gdal_merge -of netCDF -co "FORMAT=NC4" -co "BAND_NAMES=z" -o ETOPO_2022_v1_15s.nc ETOPO_2022_v1_15s_*_surface.nc
date
