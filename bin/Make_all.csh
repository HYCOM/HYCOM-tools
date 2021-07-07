#
# --- create HYCOM related executables.
#
set echo
#
#setenv OS `/bin/uname`
setenv OS `uname`
if ($OS == "Linux") then
# setenv OS XC30
# setenv OS LinuxPGF
# setenv OS LinuxGF
  setenv OS LinuxIF
endif
#
# --- the following are extracted from HYCOM-tools/config/*_setup
#
switch ($OS)
case 'LinuxPGF':
#       compile for Portland Group compiler
	setenv FC	"pgf90"
	setenv FFLAGS	"-g -fast -byteswapio -mcmodel=medium -Mnolarge_arrays"
	setenv FLIBS	""
	setenv CC	"gcc"
	setenv CFLAGS	"-O -m64 -mcmodel=medium"
	breaksw
case 'LinuxIF':
#       compile for Intel compiler
	setenv FC	"ifort"
	setenv FFLAGS	"-traceback -xHost -O3 -fp-model precise -ftz -align array64byte -convert big_endian -assume byterecl -warn nogeneral -diag-disable 10212"
	setenv FLIBS	""
	setenv CC	"icc"
	setenv CFLAGS	"-traceback -xHost -O"
	breaksw
case 'LinuxGF':
#       compile for gfortran
	setenv FC	"gfortran"
	setenv FFLAGS	"-fPIC -m64 -fno-second-underscore -fconvert=big-endian -O"
#	setenv FFLAGS	"-fPIC -m64 -fno-second-underscore -fconvert=big-endian -O -I/usr/include -I/usr/lib64/gfortran/modules"
	setenv FLIBS	""
	setenv CC	"gcc"
	setenv CFLAGS	"-fPIC -m64 -O"
	breaksw
case 'XC30':
#       compile for XC30/XC40 with ifort via aprun
        module switch PrgEnv-cray PrgEnv-intel
        setenv FC       "ftn"
        setenv FFLAGS   "-traceback -xHost -O3 -fp-model precise -ftz -align array64byte -convert big_endian -warn nogeneral -diag-disable 10212"
        setenv FLIBS    ""
        setenv CC       "cc"
        setenv CFLAGS   "-traceback -xHost -O"
        breaksw
case 'AIX':
#       compile for IBM power, note that pwr4 is likely out of date
	setenv FC	"xlf95"
	setenv FFLAGS	"-g -qfixed -O3 -qstrict -qarch=pwr4 -qtune=pwr4 -qcache=auto -qspillsize=32000 -q64 -qwarn64  -qflttrap=overflow:zerodivide:invalid:enable:imprecise -qsigtrap"
	setenv FLIBS	""
	setenv CC	"cc"
	setenv CFLAGS	"-O -q64 -DAIX"
	breaksw
default:
	echo 'Unknown Operating System: ' $OS
	exit (1)
endsw
#
# --- *.c programs
#
foreach f ( echo2 endian )
  if ( ! -e ${f}_${OS} ) then
    $CC $CFLAGS ${f}.c -o ${f}_${OS}
  else if ( -f `find ${f}.c -prune -newer ${f}_${OS}` ) then
    $CC $CFLAGS ${f}.c -o ${f}_${OS}
  else
    echo "${f}_${OS} is already up to date"
  endif
  touch       ${f}
  /bin/rm -f  ${f}
  chmod a+rx  ${f}_${OS}
  /bin/ln -s  ${f}_${OS} ${f}
end
#
# --- *.f programs
#
foreach f ( atmos_gaussian wind_to_cd hycom_sigma )
  if ( ! -e ${f}_${OS} ) then
    $FC $FFLAGS ${f}.f $FLIBS -o ${f}_${OS}
  else if ( -f `find ${f}.f -prune -newer ${f}_${OS}` ) then
    $FC $FFLAGS ${f}.f $FLIBS -o ${f}_${OS}
  else
    echo "${f}_${OS} is already up to date"
  endif
  touch       ${f}
  /bin/rm -f  ${f}
  chmod a+rx  ${f}_${OS}
  /bin/ln -s  ${f}_${OS} ${f}
end
foreach f ( clim_stat wind_stat wind_stat_check \
            wind_stat_range wind_stat_range2 wind_stat_range4 wind_stat_range5 \
            wind_stat_raw )
  if ( ! -e ${f}_${OS} ) then
    $FC $FFLAGS ${f}.f $FLIBS -o ${f}_${OS}
  else if ( -f `find ${f}.f -prune -newer ${f}_${OS}` ) then
    $FC $FFLAGS ${f}.f $FLIBS -o ${f}_${OS}
  else
    echo "${f}_${OS} is already up to date"
  endif
  touch       ${f}.exe
  /bin/rm -f  ${f}.exe
  chmod a+rx  ${f}_${OS}
  /bin/ln -s  ${f}_${OS} ${f}.exe
end
foreach f ( cice_restart cice_restart_mask cice_restart_range cice_restart_superset \
            cice_stat cice_wind_ymdh cice_x1 hycom_palette lonlat_dist \
            hycom_alat hycom_archm_dates hycom_archv_dates hycom_depth hycom_depth_40 \
            hycom_nest_dates hycom_profile+sig hycom_profile+thstar hycom_profile2pcm \
            hycom_profile2z hycom_profile2zi \
            hycom_profile_insitu_25t hycom_profile_locsig_25t \
            hycom_profile_layprs sigma2s_to_locsig \
            hycom_profile_ncoda_inc hycom_profile_mld hycom_profile_remap \
            hycom_sigma hycom_river_anom hycom_ts hycom_wind_date \
            hycom_wind_ymdh hycom_ymdh_wind hycom_yoflat \
            rhos_to_t sigma0_to_sigma2 sigma2_to_sigma0 ts_to_sigma z2zi zi2z \
            hycom_date_wind hycom_month_day hycom_profile2plm hycom_record_size \
            hycom_subset_xy hycom_dp0k hycom_dp0k_cm hycom_dp0k_sigma \
            hycom_tideport_diff hycom_tideport_scale )
  if ( ! -e ${f}_${OS} ) then
    $FC $FFLAGS ${f}.f $FLIBS -o ${f}_${OS}
  else if ( -f `find ${f}.f -prune -newer ${f}_${OS}` ) then
    $FC $FFLAGS ${f}.f $FLIBS -o ${f}_${OS}
  else
    echo "${f}_${OS} is already up to date"
  endif
  touch       ${f}
  /bin/rm -f  ${f}
  chmod a+rx  ${f}_${OS}
  /bin/ln -s  ${f}_${OS} ${f}
end
#
$FC $FFLAGS -c mom_eos_wright.f
#
foreach f ( tsp_to_insitu hycom_profile_insitu_wright )
  if ( ! -e ${f}_${OS} ) then
    $FC $FFLAGS ${f}.f  mom_eos_wright.o $FLIBS -o ${f}_${OS}
  else if ( -f `find ${f}.f -prune -newer ${f}_${OS}` ) then
    $FC $FFLAGS ${f}.f  mom_eos_wright.o $FLIBS -o ${f}_${OS}
  else
    echo "${f}_${OS} is already up to date"
  endif
  touch       ${f}
  /bin/rm -f  ${f}
  chmod a+rx  ${f}_${OS}
  /bin/ln -s  ${f}_${OS} ${f}
end
#
# --- *.F programs, may need hycom_endian_io.o and/or parse.o.
#
$FC $FFLAGS -c hycom_endian_io.F
$CC $CFLAGS -c parse.c
#
foreach f ( cice_range hycom_accumulate hycom_cfl hycom_compare \
            hycom_crosscorr hycom_crosscorr_lag hycom_join \
            track_histogram unf42hycom unf82hycom hycom2raw hycom2raw8 hycom2unf4 hycom2unf8 \
            hycom_1st_isopyc hycom_arctic hycom_arctic_ok hycom_bandmask \
            hycom_binning hycom_binning_fld hycom_bouflx hycom_clip hycom_count \
            hycom_diurnal hycom_eddy_center hycom_expr hycom_extent2ice \
            hycom_extract hycom_fill hycom_fill_sm \
            hycom_halfmask hycom_halfsm hycom_histogram hycom_ice_blend hycom_icefreeday \
            hycom_ij2lonlat hycom_islands hycom_iselect hycom_larger \
            hycom_lonlat2ij hycom_lonlat2xy hycom_mask hycom_mass hycom_mean \
            hycom_meanfit hycom_median hycom_meridional hycom_meridional_lon \
            hycom_mixlay hycom_mxthrd hycom_NaN \
            hycom_pad_fix hycom_pad_ok hycom_potdens \
            hycom_perturbation hycom_perturbation_scale hycom_print hycom_diff_print \
            hycom_random hycom_range hycom_range_ij hycom_rivers hycom_rotate \
            hycom_runmean hycom_sample hycom_sample_list hycom_sea_ok hycom_shift \
            hycom_skill hycom_slopefit hycom_smooth hycom_speed hycom_stericssh \
            hycom_subset hycom_superset hycom_swfrac hycom_sym2d hycom_thirdsm \
            hycom_tidelat hycom_transpose hycom_triple hycom_void hycom_xy2lonlat \
            hycom_zonal hycom_zonal_lat \
            ascii2hycom byte2hycom char2hycom raw2hycom raw82hycom \
            hycom_2d_ok hycom_autocorr hycom_autocorr_lag hycom_boxmean hycom_boxtime \
            hycom_index_sort hycom_mask_ok hycom_mask2_ok hycom_mass_corr hycom_newzi hycom_quadlsq \
            hycom_regression hycom_sstice hycom_botfric hycom_botslope hycom_boxsmooth \
            hycom_diflat hycom_merge hycom_sample_xy hycom_scatter hycom_tidal_ap2dist \
            hycom_tidal_ap2ri hycom_tidal_best2 hycom_tidal_foreman \
            hycom_tidal_ri2ap hycom_tidal_ri2port \
            hycom_tidal_rms hycom_tideReIm8 hycom_tidebody hycom_tideport \
            hycom_vmean hycom_xward )
  if ( ! -e ${f}_${OS} ) then
    $FC $FFLAGS ${f}.F $FLIBS hycom_endian_io.o parse.o -o ${f}_${OS}
  else if ( -f `find ${f}.F -prune -newer ${f}_${OS}` ) then
    $FC $FFLAGS ${f}.F $FLIBS hycom_endian_io.o parse.o -o ${f}_${OS}
  else
    echo "${f}_${OS} is already up to date"
  endif
  touch       ${f}
  /bin/rm -f  ${f}
  chmod a+rx  ${f}_${OS}
  /bin/ln -s  ${f}_${OS} ${f}
end
#
# --- archive reading programs, may additionally need hycom_profile_lib.o
#
if (! -e hycom_profile_lib.o  ) then
  $FC $FFLAGS -c hycom_profile_lib.F
  foreach f ( hycom_profile_list hycom_profile_argo hycom_profile_isop hycom_profile_newsig \
              hycom_profile_obs hycom_profile_offset hycom_profile_stericsshanom \
              hycom_profile_stokes hycom_bad_velocity )
    $FC $FFLAGS ${f}.F $FLIBS hycom_profile_lib.o hycom_endian_io.o parse.o -o ${f}_${OS}
  end
else if ( -f `find hycom_profile_lib.F -prune -newer hycom_profile_lib.o` ) then
  $FC $FFLAGS -c hycom_profile_lib.F
  foreach f ( hycom_profile_list hycom_profile_argo hycom_profile_isop hycom_profile_newsig \
              hycom_profile_obs hycom_profile_offset hycom_profile_stericsshanom \
              hycom_profile_stokes hycom_bad_velocity )
    $FC $FFLAGS ${f}.F $FLIBS hycom_profile_lib.o hycom_endian_io.o parse.o -o ${f}_${OS}
  end
else
  foreach f ( hycom_profile_list hycom_profile_argo hycom_profile_isop hycom_profile_newsig \
              hycom_profile_obs hycom_profile_offset hycom_profile_stericsshanom \
              hycom_profile_stokes hycom_bad_velocity )
    if ( ! -e ${f}_${OS} ) then
      $FC $FFLAGS ${f}.F $FLIBS hycom_profile_lib.o hycom_endian_io.o parse.o -o ${f}_${OS}
    else if ( -f `find ${f}.F -prune -newer ${f}_${OS}` ) then
      $FC $FFLAGS ${f}.F $FLIBS hycom_profile_lib.o hycom_endian_io.o parse.o -o ${f}_${OS}
    else
      echo "${f}_${OS} is already up to date"
    endif
  end
  touch       ${f}
  /bin/rm -f  ${f}
  chmod a+rx  ${f}_${OS}
  /bin/ln -s  ${f}_${OS} ${f}
endif
