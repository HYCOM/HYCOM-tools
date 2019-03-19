#!/bin/csh
#
#set echo
#
# --- Usage:  ./Make_scrip.csh >& Make_scrip.log
#
# --- make all SCRIP-based executables
#
# --- set SCRIP to the root directory for SCRIP
# --- available from: http://climate.lanl.gov/Software/SCRIP
#
setenv SCRIP /u/data/wallcraf/GNU/SCRIP
#
# --- set NCDF to the root directory for netCDF
# --- available from: http://www.unidata.ucar.edu/packages/netcdf/
#
source ../../Make_ncdf.src
#
# --- set ARCH to the correct value for this machine.
#
source Make_all.src
#
echo "NCDF = " $NCDF
echo "ARCH = " $ARCH
#
if (! -e ../../config/${ARCH}_setup) then
  echo "ARCH = " $ARCH "  is not supported"
  exit 1
endif
#
# --- softlink to SCRIP modules and object files
#
ln -sf ${SCRIP}/source/*.mod .
ln -sf ${SCRIP}/source/*.o   .
#
# --- softlink to netCDF module and library (and typesizes.mod for OSF1 only)
#
/bin/rm -f netcdf.mod libnetcdf.a
/bin/rm -f typesizes.mod
#
ln -s ${NCDF}/include/*.mod   .
ln -s ${NCDF}/lib/libnetcdf.a .
#
# --- SCRIP programs
#
foreach m ( isubs_field )
  make ${m} ARCH=${ARCH} >&! Make_${m}.log
  if ($status) then
    echo "Make failed:" ${m} " - see Make_${m}.log"
  else
    echo "Make worked:" ${m}
  endif
end
