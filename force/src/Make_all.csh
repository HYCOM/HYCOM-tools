#
#set echo
#
# --- Usage:  ./Make_all.csh >& Make_all.log
#
# --- make all force setup executables
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
foreach m ( ap diff kp kpc kphfc pzero riv_mon riv_hf runoff stoch time_interp time_shift tp wi wi_curl wi_ewdnwd wi_xwdywd wi_magstress wi_meanfit wc zero )
  make ${m} ARCH=${ARCH} >&! Make_${m}.log
  if ($status) then
    echo "Make failed:" ${m} " - see Make_${m}.log"
  else
    echo "Make worked:" ${m}
  endif
end
#
foreach m ( aphf_add aphf_climo aphf_diurnal aphf_extend aphf_flcorr aphf_margin aphf_meanfit aphf_monthly aphf_offset aphf_scale aphf_tacorr kphf_table tp_sal )
  make ${m} ARCH=${ARCH} >&! Make_${m}.log
  if ($status) then
    echo "Make failed:" ${m} " - see Make_${m}.log"
  else
    echo "Make worked:" ${m}
  endif
end
