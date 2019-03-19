#!/bin/csh
#set echo
#
# --- Usage:  ./Make_all.csh >& Make_all.log
# ---         but first edit ./Make_all.src for this machine.
#
# --- make all HYCOM pre/post processing executables,
# --- except those in archive/src/Make_ncdf.csh that require NetCDF
# --- plot is not included: make if from the plot or ALT/plot directory
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
setenv A $cwd
#
foreach d ( archive cice force meanstd mom4 mom6 ncom relax roms sample subregion topo )
  echo "PROCESSING ${d}/src:"
  cd ${A}/${d}/src
  csh Make_all.csh >&! Make_all.log
  cat Make_all.log
end
#
# ./bin/Make_all.csh does not use Make_all.src, but
#  should not need editing for any supported machine type.
#
echo "PROCESSING bin:"
cd ${A}/bin
csh Make_all.csh >&! Make_all.log
cat Make_all.log
