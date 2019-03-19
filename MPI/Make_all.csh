#!/bin/csh
#set echo
#
# --- Usage:  ./Make_all.csh >& Make_all.log
# ---         but first edit ./Make_all.src for this machine.
#
# --- make all HYCOM MPI-based pre/post processing executables,
# --- except those in archive/src/Make_ncdf.csh that require NetCDF
#
# --- set ARCH to the correct value for this machine.
#
#module swap mpi mpi/intel/impi/4.1.0
module list
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
foreach d ( archive/src  meanstd/src )
  echo "PROCESSING ${d}:"
  cd ${A}/${d}
  csh Make_all.csh >&! Make_all.log
  cat Make_all.log
end
