#!/bin/csh
#
#set echo
#
# --- Usage:  ./Make_all.csh >& Make_all.log
#
# --- make all archive executables (except netCDF)
#
# --- set ARCH to the correct value for this machine.
#
module swap mpi mpi/intel/impi/4.1.0
module list
source Make_all.src
#
printenv ARCH
#
if (! -e ../../config/${ARCH}_setup) then
  echo "ARCH = " $ARCH "  is not supported"
  exit 1
endif
#
# --- archive modifying programs
#
foreach m ( ncoda_archv ncoda_archv_inc ncoda_archv_inc2 ncoda_archv_new archv2data2d archv2data2t archv2data3z archv2restart restart2archv )
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
