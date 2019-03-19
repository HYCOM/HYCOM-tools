#
# --- create HYCOM related executables.
#
set echo
#
#setenv OS `/bin/uname`
setenv OS `uname`
if ($OS == "Linux") then
  if (`/bin/uname -m` == "x86_64") then
	setenv OS Linux64
  endif
# setenv OS LinuxIFC
  setenv OS LinuxICE
# setenv OS LinuxGF
# setenv OS XT5
endif
#if ($OS == "SunOS") then
#  setenv OS SunOS64
#endif
if ($OS == "UNICOS/mp") then
  setenv OS X1
endif
#
# --- the following are extracted from hycom/ALL/config/*_setup
#
switch ($OS)
case 'Linux64':
#       compile for 64-bit AMD64
	setenv FC	"pgf90"
	setenv FFLAGS	"-g -fast -byteswapio -tp k8-64 -mcmodel=medium -Mnolarge_arrays"
	setenv FLIBS	""
	setenv CC	"gcc"
	setenv CFLAGS	"-O -march=k8 -m64 -mcmodel=medium"
	breaksw
case 'LinuxICE':
#       compile for SGI Altix ICE, Intel compiler
#       compile for XC30 login nodes with ifort
        module switch PrgEnv-cray PrgEnv-intel
	setenv FC	"ifort"
#	setenv FFLAGS	"-g -O3 -fp-model source -convert big_endian"
	setenv FFLAGS   "-traceback -O3 -fp-model source -convert big_endian  -assume byterecl -mcmodel=medium -fPIC"
	setenv FLIBS	"-shared-intel"
	setenv CC	"icc"
	setenv CFLAGS	"-O"
	breaksw
case 'LinuxGF':
#       compile for gfortran
	setenv FC	"gfortran"
	setenv FFLAGS	"-fPIC -fno-second-underscore -fconvert=big-endian -O"
	setenv FLIBS	""
	setenv CC	"gcc"
	setenv CFLAGS	"-fPIC -fno-second-underscore -O"
	breaksw
case 'XT5':
#       compile for XT5 via aprun
        setenv FC       "env XTPE_INFO_MESSAGE_OFF=1 ftn"
        setenv FFLAGS   "-g -fast -byteswapio -tp barcelona-64 -mcmodel=medium -Mnolarge_arrays"
        setenv FLIBS    ""
        setenv CC       "env XTPE_INFO_MESSAGE_OFF=1 cc"
        setenv CFLAGS   "-g -O -tp barcelona-64 -mcmodel=medium -Mnolarge_arrays"
        breaksw
case 'LinuxIFC':
#       compile for Pentium 4, Intel compiler
	setenv FC	"ifort"
	setenv FFLAGS	"-g -tpp7 -O3 -convert big_endian"
	setenv FLIBS	"-Vaxlib"
	setenv CC	"gcc"
	setenv CFLAGS	"-O -march=pentium4 -m32"
	breaksw
case 'Linux':
#       compile for Pentium 4 (also 32-bit AMD64)
	setenv FC	"pgf90"
	setenv FFLAGS	"-g -fast -byteswapio -tp p7"
	setenv FLIBS	"-Mlfs"
	setenv CC	"gcc"
	setenv CFLAGS	"-O -march=pentium4 -m32"
	breaksw
case 'AIX':
	setenv FC	"xlf95"
	setenv FFLAGS	"-qfixed -O3 -qstrict -qarch=pwr3 -qtune=pwr3 -qcache=auto -qspillsize=32000 -q64"
	setenv FLIBS	""
	setenv CC	"cc"
	setenv CFLAGS	"-O -q64 -DAIX"
	breaksw
case 'IRIX64':
	setenv FC	"f90"
	setenv FFLAGS	"-g3 -64 -O3 -macro_expand"
	setenv FLIBS	""
	setenv CC	"cc"
	setenv CFLAGS	"g3 -64 -O3"
	breaksw
case 'SunOS64':
	setenv FC	"f95"
	setenv FFLAGS	"-g -fast -xarch=native64 -nodepend -xvector=no -O2 -xpp=cpp"
	setenv FLIBS	""
	setenv CC	"cc"
	setenv CFLAGS	"-g -fast -xarch=native64"
	breaksw
case 'SunOS':
	setenv FC	"f95"
	setenv FFLAGS	"-g -fast -nodepend -xvector=no -O2 -xpp=cpp"
	setenv FLIBS	""
	setenv CC	"cc"
	setenv CFLAGS	"-g -fast"
	breaksw
case 'OSF1':
	setenv FC	"f90"
	setenv FFLAGS	"-g3 -fpe1 -fast -O5 -convert big_endian -assume byterecl"
	setenv FLIBS	""
	setenv CC	"cc"
	setenv CFLAGS	"-g3 -fast"
	breaksw
case 'unicosmk':
	setenv FC	"f90"
	setenv FFLAGS	"-X 1 -V -f fixed -O scalar2,unroll2,pipeline1,vector3 -d p -M 801"
	setenv FLIBS	""
	setenv CC	"cc"
	setenv CFLAGS	""
	breaksw
case 'X1':
	setenv FC	"ftn"
	setenv FFLAGS	"-Ossp"
	setenv FLIBS	"x1_sys.o"
	setenv CC	"cc"
	setenv CFLAGS	"-hssp -UCRAY"
	$FC $FFLAGS -c x1_sys.f
	breaksw
default:
	echo 'Unknown Operating System: ' $OS
	exit (1)
endsw
#
# --- *.F programs, may need hycom_endian_io.o and/or parse.o.
#
$FC $FFLAGS -c hycom_endian_io.F
$CC $CFLAGS -c parse.c
#
# --- mixed-layer archive reading programs, may additionally need hycom_profile_lib.o
#
$FC $FFLAGS -c hycom_profile_lib.F
#
foreach f ( hycom_profile_list )
  if ( ! -e ${f}_${OS} ) then
    $FC $FFLAGS ${f}.F $FLIBS hycom_profile_lib.o hycom_endian_io.o parse.o -o ${f}_${OS}
  else if ( -f `find ${f}.F -prune -newer ${f}_${OS}` ) then
    $FC $FFLAGS ${f}.F $FLIBS hycom_profile_lib.o hycom_endian_io.o parse.o -o ${f}_${OS}
  else
    echo "${f}_${OS} is already up to date"
  endif
  touch       ${f}
  /bin/rm -f  ${f}
  chmod a+rx  ${f}_${OS}
  /bin/ln -s  ${f}_${OS} ${f}
end
