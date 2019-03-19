#!/bin/csh
#
#set echo
#
# --- Usage:  ./Make_all.csh >& Make_all.log
#
# --- make all transport sampling executables
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
# --- programs
#
foreach m ( baro_vel transport transport2 transport3 transport3_lm transp_mn transp_mn2 transp_mn3 transp_mn_2p0 meantspt mergetspt )
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
