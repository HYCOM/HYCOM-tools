#!/bin/csh
#
#set echo
#
# --- Usage:  ./Make_all.csh >& Make_all.log
#
# --- make all nemo executables
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
# --- nemo programs
#foreach m ( NONE )
#  make ${m} ARCH=${ARCH} >&! Make_${m}.log
#  if ($status) then
#    echo "Make failed:" ${m} " - see Make_${m}.log"
#  else
#    echo "Make worked:" ${m}
#  endif
#end
echo "No progrsams to make"
