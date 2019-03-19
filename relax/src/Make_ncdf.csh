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
foreach m ( ssh_modas sst_modas sst_gdem3 sst_gdem4 sst_woa z_gdem3 z_gdem4 z_gdem42 z_woa13 z_woa_tracer mom6_gdem42 )
  make ${m} ARCH=${ARCH} >&! Make_${m}.log
  if ($status) then
    echo "Make failed:" ${m} " - see Make_${m}.log"
  else
    echo "Make worked:" ${m}
  endif
end
