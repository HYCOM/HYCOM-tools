#!/bin/csh
#
#set echo
#
# --- Usage:  ./Make_ncdf.csh >& Make_ncdf.log
#
# --- make all netCDF relax executables
#
source ../../Make_ncdf.src
#
# --- set ARCH to the correct value for this machine.
#
source ../../Make_all.src
#
echo "NCDF = " $NCDF
echo "ARCH = " $ARCH
#
if (! -e ../../config/${ARCH}_setup) then
  echo "ARCH = " $ARCH "  is not supported"
  exit 1
endif
#
# --- softlink to netCDF module and library
#
/bin/rm -f netcdf.mod libnetcdf.a libnetcdff.a
/bin/rm -f netcdf.inc
/bin/rm -f typesizes.mod
#
ln -s ${NCDFC}/lib/libnetcdf.a .
ln -s ${NCDF}/lib/libnetcdff.a .
ln -s ${NCDF}/include/netcdf.mod .
ln -s ${NCDF}/include/netcdf.inc .
ln -s ${NCDF}/include/typesizes.mod .
#
# --- netCDF programs
#
foreach m ( bathy_01min landsea_01min bathy_30sec landsea_30sec bathy_15sec landsea_15sec bathy_bedmap bathy_sindhu  latlonmask hycom2mom6 hycom2mom6dx 2mom6 mom6 )
  make ${m} ARCH=${ARCH} >&! Make_${m}.log
  if ($status) then
    echo "Make failed:" ${m} " - see Make_${m}.log"
  else
    echo "Make worked:" ${m}
  endif
end
