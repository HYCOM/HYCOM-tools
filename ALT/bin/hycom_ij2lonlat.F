      PROGRAM IJ2LONLAT
      IMPLICIT NONE
C
C  hycom_ij2lonlat_ne or
C  hycom_ij2lonlat_en or
C  hycom_ij2lonlat - Usage:  hycom_ij2lonlat ip jp [grid.a]
C                            hycom_ij2lonlat [grid.a] < ij.txt
C
C     Prints the lon,lat of an ip,jp point on the HYCOM pressure grid.
C
C     A single i,j can be specified on the command line,
C     or a sequence of i,j pairs can be input from stdin.
C
C  grid.a is a hycom grid file, default regional.grid.a.
C  Note that the corresponding grid.b must also exist.
C
C  if hycom_ij2lonlat_en is used the output is degE,degN as numbers,
C  with degN between 0 and 360, which can be input to hycom_lonlat2ij
C
C  if hycom_ij2lonlat_ne is used the output is degN,degE as numbers,
C  which can be input to TPXO tidal model
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  November 2001.
C
      REAL*4, ALLOCATABLE :: PLAT(:,:),PLON(:,:)
      REAL*4              :: PAD(4096)
      INTEGER       IOS
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      LOGICAL       NEONLY,ENONLY
      INTEGER       IDM,JDM,NPAD,IP,JP,L
      CHARACTER*6   CVARIN
      CHARACTER*240 CFILEA,CFILEB
C
C     READ ARGUMENTS.
C
      CALL GETARG(0,CARG)
      L = LEN_TRIM(CARG)
*     WRITE(6,"(4a)") TRIM(CARG),'"',CARG(L-2:L),'"'
      NEONLY = CARG(L-2:L).EQ.'_ne'
      ENONLY = CARG(L-2:L).EQ.'_en'
C
      NARG = IARGC()
C
      IF     (NARG.EQ.3) THEN
        CALL GETARG(1,CARG)
        READ(CARG,*) IP
        CALL GETARG(2,CARG)
        READ(CARG,*) JP
        CALL GETARG(3,CFILEA)
      ELSEIF (NARG.EQ.2) THEN
        CALL GETARG(1,CARG)
        READ(CARG,*) IP
        CALL GETARG(2,CARG)
        READ(CARG,*) JP
        CFILEA = 'regional.grid.a'
      ELSEIF (NARG.EQ.1) THEN
        CALL GETARG(1,CFILEA)
        IP = 0
        JP = 0  ! stdin flag
      ELSEIF (NARG.EQ.0) THEN
        CFILEA = 'regional.grid.a'
        IP = 0
        JP = 0  ! stdin flag
      ELSE
        WRITE(6,*) 'Usage: hycom_ij2lonlat ip jp [grid.a]'
        CALL EXIT(1)
      ENDIF
C
C     GET IDM,JDM FROM grid.b.
C
      CFILEB = CFILEA(1:LEN_TRIM(CFILEA)-1) // 'b'
C
      OPEN(UNIT=11,FILE=CFILEB,FORM='FORMATTED',
     &     STATUS='OLD',ACTION='READ')
C
      READ( 11,*) IDM,CVARIN
      IF (CVARIN.NE.'idm   ') THEN
        WRITE(6,*) 'hycom_ij2lonlat: bad header file ',
     &             CFILEB(1:LEN_TRIM(CFILEB))
        CALL EXIT(2)
      ENDIF
      READ( 11,*) JDM,CVARIN
      IF (CVARIN.NE.'jdm   ') THEN
        WRITE(6,*) 'hycom_ij2lonlat: bad header file ',
     &             CFILEB(1:LEN_TRIM(CFILEB))
        CALL EXIT(2)
      ENDIF
C
      CLOSE(UNIT=11)
C
      IF     (IP.LT.0 .OR. IP.GT.IDM .OR.
     &        JP.LT.0 .OR. JP.GT.JDM     ) THEN
        WRITE(6,*) 'Error in hycom_ij2lonlat: ip or jp out of range'
        WRITE(6,*) 'idm,jdm = ',idm,jdm
        WRITE(6,*) 'ip,jp,  = ',ip,jp
        CALL EXIT(5)
      ENDIF
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( PLON(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_ij2lonlat: could not allocate ',
     +             IDM*JDM,' words for PLON'
        CALL EXIT(3)
      ENDIF
      ALLOCATE( PLAT(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_ij2lonlat: could not allocate ',
     +             IDM*JDM,' words for PLAT'
        CALL EXIT(3)
      ENDIF
C
      CALL LONLAT(PLON,PLAT,IDM,JDM,PAD,NPAD,
     &            IP,JP,NEONLY,ENONLY,  CFILEA)
      CALL EXIT(0)
      END
      SUBROUTINE LONLAT(PLON,PLAT,IDM,JDM,PAD,NPAD,
     +                  IP,JP, NEONLY,ENONLY, CFILEA)
      IMPLICIT NONE
C
      CHARACTER*240 CFILEA
      LOGICAL       NEONLY,ENONLY
      INTEGER       IDM,JDM,NPAD,IP,JP
      REAL*4        PLON(IDM,JDM),PLAT(IDM,JDM),PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C
      CHARACTER*1  CEW,CNS
      CHARACTER*18 CASN
      REAL*4       X,Y
      INTEGER      IOS,NRECL
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
C     READ IN THE P-GRID LON/LAT ARRAYS.
C
      INQUIRE( IOLENGTH=NRECL) PLON,PAD
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
      ENDIF
#else
      CALL ASNUNIT(11,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 11'
        write(6,*) 'ios = ',ios
        CALL EXIT(5)
      ENDIF
#endif
#endif
      OPEN(UNIT=11, FILE=CFILEA, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',CFILEA(1:LEN_TRIM(CFILEA))
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
      READ(11,REC=1,IOSTAT=IOS) PLON
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(PLON,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read ',CFILEA(1:LEN_TRIM(CFILEA))
        CALL EXIT(4)
      ENDIF
      READ(11,REC=2,IOSTAT=IOS) PLAT
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(PLAT,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read ',CFILEA(1:LEN_TRIM(CFILEA))
        CALL EXIT(4)
      ENDIF
C
C     SELECT THE DESIRED POINT.
C
      IF     (JP.NE.0) THEN  ! single i,j
        X = PLON(IP,JP)
        Y = PLAT(IP,JP)
C
C       OUTPUT X,Y.
C
        X = MOD(X+180.0,360.0)-180.0
        IF     (NEONLY) THEN
          WRITE(6,'(F9.3,F9.3)') Y,X
        ELSEIF (ENONLY) THEN
          IF     (X.LT.0.0) THEN
            WRITE(6,'(F9.3,F9.3)') X+360.0,Y
          ELSE
            WRITE(6,'(F9.3,F9.3)') X,Y
          ENDIF
        ELSE
          IF     (X.GE.0.0) THEN
            CEW = 'E'
          ELSE
            CEW = 'W'
            X   = -X
          ENDIF
          IF     (Y.GE.0.0) THEN
            CNS = 'N'
          ELSE
            CNS = 'S'
            Y   = -Y
          ENDIF
          WRITE(6,'(F8.3,A1,F8.3,A1)') X,CEW,Y,CNS
        ENDIF
      ELSE
        DO  !input loop
          READ(5,*,IOSTAT=IOS) IP,JP
          IF     (IOS.NE.0) THEN
            EXIT
          ENDIF
C
          X = PLON(IP,JP)
          Y = PLAT(IP,JP)
C
C         OUTPUT X,Y.
C
          X = MOD(X+180.0,360.0)-180.0
          IF     (NEONLY) THEN
            WRITE(6,'(F9.3,F9.3)') Y,X
          ELSEIF (ENONLY) THEN
            IF     (X.LT.0.0) THEN
              WRITE(6,'(F9.3,F9.3)') X+360.0,Y
            ELSE
              WRITE(6,'(F9.3,F9.3)') X,Y
            ENDIF
          ELSE
            IF     (X.GE.0.0) THEN
              CEW = 'E'
            ELSE
              CEW = 'W'
              X   = -X
            ENDIF
            IF     (Y.GE.0.0) THEN
              CNS = 'N'
            ELSE
              CNS = 'S'
              Y   = -Y
            ENDIF
            WRITE(6,'(F8.3,A1,F8.3,A1)') X,CEW,Y,CNS
          ENDIF
        ENDDO  !input loop
      ENDIF
      RETURN
      END
