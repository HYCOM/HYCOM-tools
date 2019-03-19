      PROGRAM XY2LONLAT
      IMPLICIT NONE
C
C  hycom_xy2lonlat - Usage:  hycom_xy2lonlat xp yp [grid.a]
C                            hycom_xy2lonlat [grid.a] < xy.txt
C
C     Prints the lon,lat of an xp,yp point on the HYCOM pressure grid.
C
C     A single x,y can be specified on the command line,
C     or a sequence of x,y pairs can be input from stdin.
C
C  grid.a is a hycom grid file, default regional.grid.a.
C  Note that the corresponding grid.b must also exist.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  November 2001.
C
      REAL*4, ALLOCATABLE :: PLAT(:,:),PLON(:,:)
      REAL*4              :: PAD(4096)
      INTEGER      IOS
      INTEGER      IARGC
      INTEGER      NARG
      CHARACTER*240 CARG
C
      INTEGER       IDM,JDM,NPAD
      REAL*4        XP,YP
      CHARACTER*6   CVARIN
      CHARACTER*240 CFILEA,CFILEB
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.3) THEN
        CALL GETARG(1,CARG)
        READ(CARG,*) XP
        CALL GETARG(2,CARG)
        READ(CARG,*) YP
        CALL GETARG(3,CFILEA)
      ELSEIF (NARG.EQ.2) THEN
        CALL GETARG(1,CARG)
        READ(CARG,*) XP
        CALL GETARG(2,CARG)
        READ(CARG,*) YP
        CFILEA = 'regional.grid.a'
      ELSEIF (NARG.EQ.1) THEN
        CALL GETARG(1,CFILEA)
        XP = 0.0
        YP = 0.0  ! stdin flag
      ELSEIF (NARG.EQ.0) THEN
        CFILEA = 'regional.grid.a'
        XP = 0.0
        YP = 0.0  ! stdin flag
      ELSE
        WRITE(6,*) 'Usage: hycom_xy2lonlat xp yp [grid.a]'
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
        WRITE(6,*) 'hycom_xy2lonlat: bad header file ',
     &             CFILEB(1:LEN_TRIM(CFILEB))
        CALL EXIT(2)
      ENDIF
      READ( 11,*) JDM,CVARIN
      IF (CVARIN.NE.'jdm   ') THEN
        WRITE(6,*) 'hycom_xy2lonlat: bad header file ',
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
      ALLOCATE( PLON(IDM+1,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_xy2lonlat: could not allocate ',
     +             (IDM+1)*JDM,' words for PLON'
        CALL EXIT(3)
      ENDIF
      ALLOCATE( PLAT(IDM+1,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_xy2lonlat: could not allocate ',
     +             (IDM+1)*JDM,' words for PLAT'
        CALL EXIT(3)
      ENDIF
C
      CALL LONLAT(PLON,PLAT,IDM,JDM,PAD,NPAD, XP,YP, CFILEA)
      CALL EXIT(0)
      END
      SUBROUTINE LONLAT(PLON,PLAT,IDM,JDM,PAD,NPAD, XP,YP, CFILEA)
      IMPLICIT NONE
C
      CHARACTER*240 CFILEA
      INTEGER       IDM,JDM,NPAD
      REAL*4        XP,YP
      REAL*4        PLON(IDM+1,JDM),PLAT(IDM+1,JDM),PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C
      CHARACTER*1  CEW,CNS
      CHARACTER*18 CASN
      REAL*4       X,Y,DX,DY
      INTEGER      IP,JP,J
      INTEGER      IOS,NRECL
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
C     READ IN THE P-GRID LON/LAT ARRAYS.
C
      INQUIRE( IOLENGTH=NRECL) PLON(1:IDM,1:JDM),PAD
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
      PLON(IDM+1,:) = 0.0
      READ(11,REC=1,IOSTAT=IOS) PLON(1:IDM,1:JDM)
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(PLON,(IDM+1)*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read ',CFILEA(1:LEN_TRIM(CFILEA))
        CALL EXIT(4)
      ENDIF
      PLON(IDM+1,:) = 0.0
      READ(11,REC=2,IOSTAT=IOS) PLAT(1:IDM,1:JDM)
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(PLAT,(IDM+1)*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read ',CFILEA(1:LEN_TRIM(CFILEA))
        CALL EXIT(4)
      ENDIF
C
C     PERIODIC OVERLAP
C
      IF     (YP.EQ.0.0 .OR. INT(XP).EQ.IDM) THEN
        DO J= 1,JDM
          X = MOD(PLON(1,J) - PLON(IDM,J),360.0)
          IF     (X.LT.-180.0) THEN
            X = X + 360.0
          ELSEIF (X.GT. 180.0) THEN
            X = X - 360.0
          ENDIF
          PLON(IDM+1,J) = PLON(IDM,J) + X
          PLAT(IDM+1,J) = PLAT(1,  J)
        ENDDO
      ENDIF
C
C     SELECT THE DESIRED POINT.
C
      IF     (YP.NE.0.0) THEN  ! single x,y
        IF     (XP.LT.1 .OR. INT(XP).GT.IDM .OR.
     &          YP.LT.1 .OR.     YP .GT.JDM     ) THEN
          WRITE(6,*) 'Error in hycom_xy2lonlat: xp or yp out of range'
          WRITE(6,*) 'idm,jdm = ',idm,jdm
          WRITE(6,*) 'xp,yp,  = ',xp,yp
          CALL EXIT(5)
        ENDIF
C
        IP = XP
        JP = YP
        DX = XP - IP
        DY = YP - JP
        X = (1.0-DX)*(1.0-DY)*PLON(IP,  JP)   + 
     +      (1.0-DX)*     DY *PLON(IP,  JP+1) + 
     +           DX *(1.0-DY)*PLON(IP+1,JP)   + 
     +           DX *     DY *PLON(IP+1,JP+1)
        Y = (1.0-DX)*(1.0-DY)*PLAT(IP,  JP)   + 
     +      (1.0-DX)*     DY *PLAT(IP,  JP+1) + 
     +           DX *(1.0-DY)*PLAT(IP+1,JP)   + 
     +           DX *     DY *PLAT(IP+1,JP+1)
C
C       OUTPUT X,Y.
C
        X = MOD(X+180.0,360.0)-180.0
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
      ELSE
        DO  !input loop
          READ(5,*,IOSTAT=IOS) XP,YP
          IF     (IOS.NE.0) THEN
            EXIT
          ENDIF
C
          IF     (XP.LT.1 .OR. INT(XP).GT.IDM .OR.
     &            YP.LT.1 .OR.     YP .GT.JDM     ) THEN
            WRITE(6,*) 'Error in hycom_xy2lonlat: xp or yp out of range'
            WRITE(6,*) 'idm,jdm = ',idm,jdm
            WRITE(6,*) 'xp,yp,  = ',xp,yp
            CALL EXIT(5)
          ENDIF
C
          IP = XP
          JP = YP
          DX = XP - IP
          DY = YP - JP
          X = (1.0-DX)*(1.0-DY)*PLON(IP,  JP)   + 
     +        (1.0-DX)*     DY *PLON(IP,  JP+1) + 
     +             DX *(1.0-DY)*PLON(IP+1,JP)   + 
     +             DX *     DY *PLON(IP+1,JP+1)
          Y = (1.0-DX)*(1.0-DY)*PLAT(IP,  JP)   + 
     +        (1.0-DX)*     DY *PLAT(IP,  JP+1) + 
     +             DX *(1.0-DY)*PLAT(IP+1,JP)   + 
     +             DX *     DY *PLAT(IP+1,JP+1)
C
C         OUTPUT X,Y.
C
          X = MOD(X+180.0,360.0)-180.0
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
        ENDDO  !input loop
      ENDIF
      RETURN
      END
