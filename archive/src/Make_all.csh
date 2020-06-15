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
source ../../Make_all.src
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
foreach m ( hycomarchv botlay_archv conv_archv hybgen_archv mrgl_archv \
            ncoda_archv ncoda_archv_inc ncoda_archv_inc2 \
            remap_archv remaph_archv remapi_archv \
            insitu_archv potden_archv steric_archv stats_archv trim_archv \
            archt2archv archv2data2d archv2data2b archv2data2t archv2data3z \
            archv2datasf archv2datasfl archv2datasfz \
            archv2ncombc archv2restart \
            field2data field2data3z restart2archv archv2oam archv2oam_m )
  make ${m} ARCH=${ARCH} >&! Make_${m}.log
  if ($status) then
    echo "Make failed:" ${m} " - see Make_${m}.log"
#   cat Make_${m}.log
  else
    echo "Make worked:" ${m}
  endif
  if (-e /usr/bin/ldedit) then
#   try to set medium pages on POWER5+
    /usr/bin/ldedit -bdatapsize=64K -bstackpsize=64K ${m}
  endif
end
