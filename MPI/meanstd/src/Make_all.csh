#! /bin/csh
#set echo
#
# --- Usage:  ./Make_all.com >& Make_all.log
#
# --- make all meanstd executables
#
# --- set ARCH to the correct value for this machine.
#
module list
source Make_all.src
#
printenv ARCH
#
if (! -e ../../config/${ARCH}_setup) then
  echo "ARCH = " $ARCH_setup" is not supported"
  exit 1
endif
#
# --- meanstd programs
#
foreach m ( mean std diff )
  make ${m} ARCH=${ARCH} >&! Make_${m}.log
  if ($status) then
    echo "Make failed:" ${m} " - see Make_${m}.log"
#   cat Make_${m}.log
  else
    echo "Make worked:" ${m}
  endif
  if (-e /usr/bin/ldedit) then
#   try to set medium pages on POWER5+
    /usr/bin/ldedit -bdatapsize=64K -bstackpsize=64K hycom_${m}
  endif
end
