      PROGRAM HYCOM_DIFLAT
      IMPLICIT NONE
C
C  hycom_diflat - Usage:  hycom_diflat diflat.a difbac [grid.a]
C
C                 outputs latitudinaly varying diffusion
C
C                 difbac is the background/internal wave diffusivity (m^2/s)
C
C                 grid.a is a hycom grid file, default regional.grid.a.
C                 Note that the corresponding grid.b must also exist.
C
C                 idm,jdm are taken from grid.a
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  February 2009.
C
      REAL*4,  ALLOCATABLE :: D(:,:),PLAT(:,:)
      REAL*4               :: PAD(4096)
      INTEGER IOS
      INTEGER      IARGC
      INTEGER      NARG
      CHARACTER*240 CARG
C
      INTEGER       IDM,JDM,NPAD
      REAL          DIFBAC
      CHARACTER*6   CVARIN
      CHARACTER*240 CFILE,CFILEB,CFILEG
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.2) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) DIFBAC
        CFILEG = 'regional.grid.a'
      ELSEIF (NARG.EQ.3) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) DIFBAC
        CALL GETARG(3,CFILEG)
      ELSE
        WRITE(6,*) 
     +   'Usage:  hycom_diflat difout.a difbac [grid.a]'
        CALL EXIT(1)
      ENDIF
C
C     GET IDM,JDM FROM grid.b.
C
      CFILEB = CFILEG(1:LEN_TRIM(CFILEG)-1) // 'b'
C
      OPEN(UNIT=11,FILE=CFILEB,FORM='FORMATTED',
     &     STATUS='OLD',ACTION='READ')
C
      READ( 11,*) IDM,CVARIN
      IF (CVARIN.NE.'idm   ') THEN
        WRITE(6,*) 'hycom_diflat: bad header file ',
     &             CFILEB(1:LEN_TRIM(CFILEB))
        CALL EXIT(2)
      ENDIF
      READ( 11,*) JDM,CVARIN
      IF (CVARIN.NE.'jdm   ') THEN
        WRITE(6,*) 'hycom_diflat: bad header file ',
     &             CFILEB(1:LEN_TRIM(CFILEB))
        CALL EXIT(2)
      ENDIF
C
      CLOSE(UNIT=11)
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( D(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_diflat: could not allocate ',
     +             IDM*JDM,' words for D'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( PLAT(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_diflat: could not allocate ',
     +             IDM*JDM,' words for PLAT'
        CALL EXIT(2)
      ENDIF
C
      CALL DIFLAT(D,PLAT,IDM,JDM,PAD,NPAD,
     +            DIFBAC, CFILE,CFILEG)
      CALL EXIT(0)
 5000 FORMAT(I4)
      END
      SUBROUTINE DIFLAT(D,PLAT,IDM,JDM,PAD,NPAD,
     +                  DIFBAC, CFILE,CFILEG)
      IMPLICIT NONE
C
      CHARACTER*240 CFILE,CFILEG
      INTEGER      IDM,JDM,NPAD
      REAL         DIFBAC
      REAL*4       D(IDM,JDM),PLAT(IDM,JDM),PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      CHARACTER*18 CASN
      INTEGER      I,J,IOS,NRECL
      REAL         APLAT
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
C     INPUT PLAT ARRAY.
C
      INQUIRE( IOLENGTH=NRECL) PLAT,PAD
C
#ifdef CRAY
#ifdef t3e
      IF     (MOD(NRECL,4096).EQ.0) THEN
        WRITE(CASN,8000) NRECL/4096
 8000   FORMAT('-F cachea:',I4.4,':1:0')
        IU8 = 11
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          write(6,*) 'Error: can''t asnunit 11'
          write(6,*) 'ios  = ',ios8
          write(6,*) 'casn = ',casn
          CALL EXIT(5)
        ENDIF
        IU8 = 21
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          write(6,*) 'Error: can''t asnunit 21'
          write(6,*) 'ios  = ',ios8
          write(6,*) 'casn = ',casn
          CALL EXIT(5)
        ENDIF
      ENDIF
#else
      CALL ASNUNIT(11,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 11'
        write(6,*) 'ios = ',ios
        CALL EXIT(5)
      ENDIF
      CALL ASNUNIT(21,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 21'
        write(6,*) 'ios = ',ios
        CALL EXIT(5)
      ENDIF
#endif
#endif
C
      OPEN(UNIT=11, FILE=CFILEG, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEG)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
      READ(11,REC=2,IOSTAT=IOS) PLAT
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(PLAT,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read ',TRIM(CFILEG)
        CALL EXIT(4)
      ENDIF
C
      CLOSE(UNIT=11)
C
C     OUTPUT FILE.
C
      OPEN(UNIT=21, FILE=CFILE, FORM='UNFORMATTED', STATUS='NEW',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILE)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
C     CREATE LATITUDINALLY DEPENDENT DIFFUSION.
C
      DO J= 1,JDM
        DO I= 1,IDM
          APLAT = ABS(PLAT(I,J))
          IF     (APLAT.LE.28.9) THEN
            D(I,J) = 1.e-6 + 2.0*DIFBAC*APLAT/28.9
          ELSE
            D(I,J) = max( DIFBAC,
     &                    1.e-6 + 2.0*DIFBAC*(2.0*28.9-APLAT)/28.9 )
          ENDIF
        ENDDO
      ENDDO
C
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(D,IDM*JDM)
#endif
      IF     (NPAD.EQ.0) THEN
        WRITE(21,REC=1,IOSTAT=IOS) D
      ELSE
        WRITE(21,REC=1,IOSTAT=IOS) D,PAD
      ENDIF
      CLOSE(21)
C
      RETURN
      END
