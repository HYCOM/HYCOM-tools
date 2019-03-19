#!/bin/csh
#
#set echo
#
# --- Usage:  ./Make_ncdf.com >& Make_ncdf.log
#
# --- make all netCDF archive executables
#
# --- set NCDF to the root directory for netCDF version 4.3 or higher
# --- available from: http://www.unidata.ucar.edu/packages/netcdf/
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
# --- netCDF archive programs
#
foreach m ( archv2strmf archv2ncdf2d archv2ncdf2t archv2ncdf3z )
  make ${m} ARCH=${ARCH} >&! Make_${m}.log
  if ($status) then
    echo "Make failed:" ${m} " - see Make_${m}.log"
  else
    echo "Make worked:" ${m}
  endif
  if (-e /usr/bin/ldedit) then
#   try to set medium pages on POWER5+
    /usr/bin/ldedit -bdatapsize=64K -bstackpsize=64K ${m}
  endif
end
