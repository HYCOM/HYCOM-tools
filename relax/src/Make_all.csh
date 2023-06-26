#
#set echo
#
# --- Usage:  ./Make_all.csh >& Make_all.log
#
# --- make all relax setup executables, except those needing netCDF.
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
foreach m ( bottom iso_density relax_flat_rivers relax_tracer relaxi relaxi_dens relaxv relax_archive rmu rmu2 rmu_linear tracer_const z_archive z_const z_levitus z_medatlas z_modify sst_pf sst_pf_4km )
  make ${m} ARCH=${ARCH} >&! Make_${m}.log
  if ($status) then
    echo "Make failed:" ${m} " - see Make_${m}.log"
#   cat Make_${m}.log
  else
    echo "Make worked:" ${m}
  endif
end
