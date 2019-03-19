#
# --- create HYCOM related netCDF executables.
# --- run after Make_all.com
#
set echo
#
# --- set NCDFC to the root directory for netCDF version 4.3.
# --- set NCDF  to the root directory for netCDF version 4.3 Fortran.
# --- Use EXTRANCDF for the libraries needed for NCDF v4.3 or later
#
source ../Make_ncdf.src
#
# --- softlink to netCDF module and library (and typesizes.mod if needed)
#
/bin/rm -f netcdf.mod libnetcdf.a libnetcdff.a
/bin/rm -f netcdf.inc
/bin/rm -f typesizes.mod
#
ln -s ${NCDFC}/lib/libnetcdf.a .
ln -s ${NCDF}/lib/libnetcdff.a .
ln -s ${NCDF}/include/netcdf.mod .
ln -s ${NCDF}/include/netcdf.inc .
ln -s ${NCDF}/include/typesizes.mod .
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
# --- the following are extracted from hycom/ALL/config/*_setup
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
        setenv FFLAGS   "-traceback -xHost -O3 -fp-model precise -ftz -align array64 byte -convert big_endian -warn nogeneral -diag-disable 10212"
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
foreach f ( hycom2nc hycom_binning_nc hycom_scrip_nc )
  if ( ! -e ${f}_${OS} ) then
    $FC $FFLAGS ${f}.F hycom_endian_io.o parse.o ${EXTRANCDF} -o ${f}_${OS}
  else if ( -f `find ${f}.F -prune -newer ${f}_${OS}` ) then
    $FC $FFLAGS ${f}.F hycom_endian_io.o parse.o ${EXTRANCDF} -o ${f}_${OS}
  else
    echo "${f}_${OS} is already up to date"
  endif
  touch       ${f}
  /bin/rm -f  ${f}
  chmod a+rx  ${f}_${OS}
  /bin/ln -s  ${f}_${OS} ${f}
end
#
foreach f ( hycom_profile2z_nc hycom_profile2s_nc hycom_seaice_nc )
  if ( ! -e ${f}_${OS} ) then
    $FC $FFLAGS ${f}.F hycom_profile_lib.o hycom_endian_io.o parse.o ${EXTRANCDF} -o ${f}_${OS}
  else if ( -f `find ${f}.F -prune -newer ${f}_${OS}` ) then
    $FC $FFLAGS ${f}.F hycom_profile_lib.o hycom_endian_io.o parse.o ${EXTRANCDF} -o ${f}_${OS}
  else
    echo "${f}_${OS} is already up to date"
  endif
  touch       ${f}
  /bin/rm -f  ${f}
  chmod a+rx  ${f}_${OS}
  /bin/ln -s  ${f}_${OS} ${f}
end
#
foreach f ( wind_stat_nc wind_stat_range_nc )
  if ( ! -e ${f}_${OS} ) then
    $FC $FFLAGS ${f}.f ${EXTRANCDF} -o ${f}_${OS}
  else if ( -f `find ${f}.f -prune -newer ${f}_${OS}` ) then
    $FC $FFLAGS ${f}.f ${EXTRANCDF} -o ${f}_${OS}
  else
    echo "${f}_${OS} is already up to date"
  endif
  touch       ${f}
  /bin/rm -f  ${f}
  chmod a+rx  ${f}_${OS}
end
