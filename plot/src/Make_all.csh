#!/bin/csh 
#set echo
#
# --- Usage:  ./Make_all.com >& Make_all.log
#
# --- make all plot executables
#
# --- set ARCH to the correct value for this machine.
#
#setenv NCARG_ROOT ~/pkgs/ncarg/gnu/5.2.1
#set path=($NCARG_ROOT/bin $path)
which ncargf90
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
# --- standard plot programs
#
foreach m ( hycomproc fieldproc fieldcell )
  make ${m} ARCH=${ARCH} >&! Make_${m}.log
  if ($status) then
    echo "Make failed:" ${m} "- see Make_${m}.log"
  else
    echo "Make worked:" ${m}
  endif
end
#
# --- executables for specific output media.
#
foreach p ( hp fp fc )
  foreach m ( meta psp psl x11 )
    make ${p}_${m} ARCH=${ARCH} >&! Make_${p}_${m}.log
    if ($status) then
      echo "Make failed:" ${p}_${m} " - see Make_${p}_${m}.log"
    else
      echo "Make worked:" ${p}_${m}
    endif
  end
end
