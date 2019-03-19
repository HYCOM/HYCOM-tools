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
foreach d ( archive/src_* meanstd/src_* )
  echo "CLEANING ${d}"
  cd ${A}/${d}
  make clean
end
