#!/bin/csh
#set echo
#
# --- Usage:  ./Make_ncdf.csh >& Make_ncdf.log
# ---         but first edit ./Make_ncdf.src for this machine.
#
# --- make all MPI-based HYCOM pre/post processing executables that require NetCDF.
#
# --- set ARCH to the correct value for this machine.
#
#module swap mpi mpi/intel/impi/4.1.0
#module list
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
setenv A $cwd
#
foreach d ( archive/src )
  echo "PROCESSING ${d}:"
  cd ${A}/${d}
  csh Make_ncdf.csh >&! Make_ncdf.log
  cat Make_ncdf.log
end
