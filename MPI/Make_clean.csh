#!/bin/csh
#
#set echo
#
# --- Usage:  ./Make_clean.com
#
# --- clean all HYCOM pre/post processing source directories
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
foreach d ( archive meanstd )
  echo "CLEANING ${d}/src:"
  cd ${A}/${d}/src
  make clean ARCH=$ARCH
end
