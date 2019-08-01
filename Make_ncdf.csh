#!/bin/csh
#set echo
#
# --- Usage:  ./Make_ncdf.csh >& Make_ncdf.log
# ---         but first edit ./Make_ncdf.src for this machine.
#
# --- make all HYCOM pre/post processing executables that require NetCDF.
#
# --- set ARCH to the correct value for this machine.
#
source Make_all.src
#
printenv ARCH
#
if (! -e ./config/${ARCH}_setup) then
  echo "ARCH = " $ARCH "  is not supported"
  exit 1
endif
#
# --- set NCDF to the correct value for this machine.
#
source Make_ncdf.src
#
printenv NCDF
#
#if (! -e ${NCDF}/lib/libnetcdf.a) then
#  echo "NCDF = " $NCDF "  is not correct"
#  exit 2
#endif
#
setenv A $cwd
#
foreach d ( archive force relax topo )
  echo "PROCESSING ${d}/src:"
  cd ${A}/${d}/src
  csh Make_ncdf.csh >&! Make_ncdf.log
  cat Make_ncdf.log
end
#
# ./bin/Make_ncdf.csh does not use Make_ncdf.src, but
#  should not need editing for any supported machine type.
#
echo "PROCESSING bin:"
cd ${A}/bin
csh Make_ncdf.csh >&! Make_ncdf.log
cat Make_ncdf.log
