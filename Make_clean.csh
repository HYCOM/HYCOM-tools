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
foreach d ( archive cice force meanstd mom4 mom6 ncom relax roms sample subregion topo )
  echo "CLEANING ${d}/src:"
  cd ${A}/${d}/src
  make clean ARCH=$ARCH
end
#
# ./bin/Make_clean.com does not use Make_clean.src, but
#  should not need editing for any supported machine type.
#
echo "PROCESSING bin:"
cd ${A}/bin
csh Make_clean.csh >&! Make_clean.log
cat Make_clean.log
