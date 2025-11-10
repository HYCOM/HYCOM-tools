#
#set echo
#
# --- Usage:  ./Make_all.csh >& Make_all.log
#
# --- make all topo setup executables
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
foreach m ( cfl2spd spd2cfl bathy_02min bathy_05min hudson landsea_02min landsea_05min \
            1d 2d batrop biharmonic clip diff distland dry edit flat grid_360 hgram \
            landfill landmask latitude lonlat lonlat_2d lpanam map mapsub \
            mercator merge modify onesea onesea-b onesea_fill onesea-b_fill \
            panam partit partit_noshrink partit_arctic partit_arctic_ns \
            partit_mom6 partit_landsea partit_resize partit_row1st partit_test \
            ppmX diff_ppmX ports ports_latlon rotated rough shrink slope \
            smallsea smooth smooth_skip smooth_slope subset tid_mask tid_merge \
            tiles vrtcmprs zcells zthin )
  make ${m} ARCH=${ARCH} >&! Make_${m}.log
  if ($status) then
    echo "Make failed:" ${m} " - see Make_${m}.log"
  else
    echo "Make worked:" ${m}
  endif
end
