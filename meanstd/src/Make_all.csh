#
#set echo
#
# --- Usage:  ./Make_all.csh >& Make_all.log
#
# --- make all meanstd executables
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
# --- meanstd programs
#
foreach m ( hycom_mean hycom_std hycom_diff hycom_wsum hesmf_mean hesmf_std )
  make ${m} ARCH=${ARCH} >&! Make_${m}.log
  if ($status) then
    echo "Make failed:" ${m} " - see Make_${m}.log"
  else
    echo "Make worked:" ${m}
  endif
  if (-e /usr/bin/ldedit) then
#   try to set medium pages on POWER5+ and POWER6
    /usr/bin/ldedit -bdatapsize=64K -bstackpsize=64K ${m}
  endif
end
