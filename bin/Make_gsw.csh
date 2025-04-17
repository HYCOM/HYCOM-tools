#
# --- create HYCOM related executables that use the GSW TEOS-10 toolkit.
# --- To compile with ifort, update gsw_saar.f90 and gsw_mod_saar_data.f90
# --- https://github.com/jamesorr/mocsy/commit/65ac49e307cc596dd2e965eb3c1606d7025a3c65
#
set echo
#
setenv GSW_PATH /p/home/wallcraf/pkgs/GSW-Fortran-3.05-6-ifort
#
setenv GSW_OBJM  "$GSW_PATH/modules/gsw_mod_kinds.o $GSW_PATH/modules/gsw_mod_teos10_constants.o $GSW_PATH/modules/gsw_mod_toolbox.o $GSW_PATH/modules/gsw_mod_error_functions.o $GSW_PATH/modules/gsw_mod_baltic_data.o $GSW_PATH/modules/gsw_mod_saar_data.o $GSW_PATH/modules/gsw_mod_specvol_coefficients.o"
setenv GSW_OBJT "$GSW_PATH/toolbox/gsw_sa_from_sp.o $GSW_PATH/toolbox/gsw_ct_from_pt.o $GSW_PATH/toolbox/gsw_rho.o $GSW_PATH/toolbox/gsw_pt_from_ct.o $GSW_PATH/toolbox/gsw_ct_from_rho.o $GSW_PATH/toolbox/gsw_sa_from_sp_baltic.o $GSW_PATH/toolbox/gsw_saar.o $GSW_PATH/toolbox/gsw_specvol.o $GSW_PATH/toolbox/gsw_add_barrier.o $GSW_PATH/toolbox/gsw_add_mean.o $GSW_PATH/toolbox/gsw_util_indx.o $GSW_PATH/toolbox/gsw_util_xinterp1.o $GSW_PATH/toolbox/gsw_alpha.o"
setenv GSW_OBJX "$GSW_PATH/toolbox/gsw_ct_freezing_poly.o $GSW_PATH/toolbox/gsw_ct_maxdensity.o $GSW_PATH/toolbox/gsw_gibbs_pt0_pt0.o $GSW_PATH/toolbox/gsw_rho_alpha_beta.o"
#
foreach f ( $GSW_PATH/modules/gsw_mod_kinds.mod $GSW_PATH/modules/gsw_mod_teos10_constants.mod $GSW_PATH/modules/gsw_mod_toolbox.mod $GSW_PATH/modules/gsw_mod_error_functions.mod $GSW_PATH/modules/gsw_mod_baltic_data.mod $GSW_PATH/modules/gsw_mod_saar_data.mod $GSW_PATH/modules/gsw_mod_specvol_coefficients.mod )
  ln -sf $f .
end
#
#setenv OS `/bin/uname`
setenv OS `uname`
if ($OS == "Linux") then
# setenv OS XC30
# setenv OS LinuxPGF
# setenv OS LinuxGF
# setenv OS LinuxIF
  setenv OS LinuxAIF
endif
#
# --- the following are extracted from HYCOM-tools/config/*_setup
#
switch ($OS)
case 'LinuxPGF':
#       compile for Portland Group compiler
	setenv FC	"pgf90"
	setenv FFLAGS	"-g -fast -byteswapio -mcmodel=medium -Mnolarge_arrays"
	setenv FLIBS	"$GSW_OBJM $GSW_OBJT $GSW_OBJX mom_eos_wright.o"
	setenv CC	"gcc"
	setenv CFLAGS	"-O -m64 -mcmodel=medium"
	breaksw
case 'LinuxIF':
#       compile for Intel compiler
	setenv FC	"ifort"
	setenv FFLAGS	"-traceback -xHost -O3 -fp-model precise -ftz -align array64byte -convert big_endian -assume byterecl -warn nogeneral -diag-disable 10212"
	setenv FLIBS	"$GSW_OBJM $GSW_OBJT $GSW_OBJX mom_eos_wright.o"
	setenv CC	"icc"
	setenv CFLAGS	"-traceback -xHost -O"
	breaksw
case 'LinuxAIF':
#       compile for Intel compiler and AMD processors
	setenv FC	"ifort"
	setenv FFLAGS	"-traceback -march=core-avx2 -O3 -fp-model precise -ftz -align array64byte -convert big_endian -assume byterecl -warn nogeneral -diag-disable 10212"
	setenv FLIBS	"$GSW_OBJM $GSW_OBJT $GSW_OBJX mom_eos_wright.o"
	setenv CC	"icc"
	setenv CFLAGS	"-traceback -O"
	breaksw
case 'LinuxGF':
#       compile for gfortran
	setenv FC	"gfortran"
	setenv FFLAGS	"-fPIC -m64 -fno-second-underscore -fconvert=big-endian -O"
#	setenv FFLAGS	"-fPIC -m64 -fno-second-underscore -fconvert=big-endian -O -I/usr/include -I/usr/lib64/gfortran/modules"
	setenv FLIBS	"$GSW_OBJM $GSW_OBJT $GSW_OBJX mom_eos_wright.o"
	setenv CC	"gcc"
	setenv CFLAGS	"-fPIC -m64 -O"
	breaksw
case 'XC30':
#       compile for XC30/XC40 with ifort via aprun
        module switch PrgEnv-cray PrgEnv-intel
        setenv FC       "ftn"
        setenv FFLAGS   "-traceback -xHost -O3 -fp-model precise -ftz -align array64byte -convert big_endian -warn nogeneral -diag-disable 10212"
	setenv FLIBS	"$GSW_OBJM $GSW_OBJT $GSW_OBJX mom_eos_wright.o"
        setenv CC       "cc"
        setenv CFLAGS   "-traceback -xHost -O"
        breaksw
case 'AIX':
#       compile for IBM power, note that pwr4 is likely out of date
	setenv FC	"xlf95"
	setenv FFLAGS	"-g -qfixed -O3 -qstrict -qarch=pwr4 -qtune=pwr4 -qcache=auto -qspillsize=32000 -q64 -qwarn64  -qflttrap=overflow:zerodivide:invalid:enable:imprecise -qsigtrap"
	setenv FLIBS	"$GSW_OBJM $GSW_OBJT $GSW_OBJX mom_eos_wright.o"
	setenv CC	"cc"
	setenv CFLAGS	"-O -q64 -DAIX"
	breaksw
default:
	echo 'Unknown Operating System: ' $OS
	exit (1)
endsw
#
# --- *.f programs
#
$FC $FFLAGS -c mom_eos_wright.f
#
foreach f ( tsp_to_insitu_gsw rsp_to_insitu_gsw rsp_to_fcompr_gsw rsp_to_sigma2_gsw )
  if ( ! -e ${f}_${OS} ) then
    $FC $FFLAGS -c ${f}.f -o ${f}.o
    $FC $FFLAGS $FLIBS ${f}.o -o ${f}_${OS}
  else if ( -f `find ${f}.f -prune -newer ${f}_${OS}` ) then
    $FC $FFLAGS -c ${f}.f -o ${f}.o
    $FC $FFLAGS $FLIBS ${f}.o -o ${f}_${OS}
  else
    echo "${f}_${OS} is already up to date"
  endif
  touch       ${f}
  /bin/rm -f  ${f}
  chmod a+rx  ${f}_${OS}
  /bin/ln -s  ${f}_${OS} ${f}
end
#
/bin/rm gsw_mod*.mod
