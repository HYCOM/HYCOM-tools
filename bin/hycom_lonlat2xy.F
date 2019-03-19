      PROGRAM LONLAT2XY
      IMPLICIT NONE
C
C  hycom_lonlat2xy - Usage:  hycom_lonlat2xy lon lat [grid.a]
C                            hycom_lonlat2xy [grid.a]  < lonlat.txt
C                            hycom_lonlat2xy_mask lon lat mask.a [grid.a]
C                            hycom_lonlat2xy_mask mask.a [grid.a] < lonlat.txt
C                            hycom_lonlat2ij lon lat [grid.a]
C                            hycom_lonlat2ij [grid.a]  < lonlat.txt
C                            hycom_lonlat2ij_mask lon lat mask.a [grid.a]
C                            hycom_lonlat2ij_mask mask.a [grid.a] < lonlat.txt
C
C     Prints the HYCOM p-grid location of lon,lat or 0,0 if out of range.
C     The "2ij" versions print the nearest array point location.
C
C     A single lon,lat can be specified on the command line,
C     or a sequence of lon,lat pairs can be input from stdin.
C
C     A "track" file can also be input, since lines without a lon-lat
C     pair are passed through to the output.
C
C  grid.a is a hycom grid file, default regional.grid.a.
C  Note that the corresponding grid.b must also exist.
C  for efficiency, distance is dx+dy in m, rather that sqrt(dx**2+dy**2)
C
C  mask.a is a hycom mask file.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  November 2001.
C
      REAL*4, ALLOCATABLE :: PLAT(:,:),PLON(:,:),AMSK(:,:)
      REAL*4                 PAD(4096)
      INTEGER      IOS,L
      INTEGER      IARGC
      INTEGER      NARG
      CHARACTER*240 CARG
C
      LOGICAL       L2IJ
      INTEGER       IDM,JDM,NPAD
      REAL*4        XP,YP
      CHARACTER*6   CVARIN
      CHARACTER*240 CFILEA,CFILEB,CFILEM
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      CALL GETARG(0,CARG)
      L = LEN_TRIM(CARG)
*     WRITE(6,"(4a)") TRIM(CARG),'"',CARG(L-4:L),'"'
      IF     (CARG(L-4:L).EQ.'_mask') THEN
        L2IJ = CARG(L-7:L-5).EQ.'2ij'
        IF     (NARG.EQ.4) THEN
          CALL GETARG(1,CARG)
          READ(CARG,*) XP
          CALL GETARG(2,CARG)
          READ(CARG,*) YP
          CALL GETARG(3,CFILEM)
          CALL GETARG(4,CFILEA)
        ELSEIF (NARG.EQ.3) THEN
          CALL GETARG(1,CARG)
          READ(CARG,*) XP
          CALL GETARG(2,CARG)
          CALL GETARG(3,CFILEM)
          READ(CARG,*) YP
          CFILEA = 'regional.grid.a'
        ELSEIF (NARG.EQ.2) THEN
          CALL GETARG(1,CFILEM)
          CALL GETARG(2,CFILEA)
          XP = 0.0
          YP = -999.0  ! stdin flag
        ELSEIF (NARG.EQ.1) THEN
          CALL GETARG(1,CFILEM)
          CFILEA = 'regional.grid.a'
          XP = 0.0
          YP = -999.0  ! stdin flag
        ELSEIF (L2IJ) THEN
          WRITE(6,*) 
     &      'Usage: hycom_lonlat2ij_mask lon lat mask.a [grid.a]'
          CALL EXIT(1)
        ELSE
          WRITE(6,*) 
     &      'Usage: hycom_lonlat2xy_mask lon lat mask.a [grid.a]'
          CALL EXIT(1)
        ENDIF
      ELSE !hycom_lonlat2xy
        L2IJ   = CARG(L-2:L).EQ.'2ij'
        CFILEM = ""
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
          YP = -999.0  ! stdin flag
        ELSEIF (NARG.EQ.0) THEN
          CFILEA = 'regional.grid.a'
          XP = 0.0
          YP = -999.0  ! stdin flag
        ELSEIF (L2IJ) THEN
          WRITE(6,*) 'Usage: hycom_lonlat2ij lon lat [grid.a]'
          CALL EXIT(1)
        ELSE
          WRITE(6,*) 'Usage: hycom_lonlat2xy lon lat [grid.a]'
          CALL EXIT(1)
        ENDIF
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
        WRITE(6,*) 'hycom_lonlat2xy: bad header file ',
     &             CFILEB(1:LEN_TRIM(CFILEB))
        CALL EXIT(2)
      ENDIF
      READ( 11,*) JDM,CVARIN
      IF (CVARIN.NE.'jdm   ') THEN
        WRITE(6,*) 'hycom_lonlat2xy: bad header file ',
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
      ALLOCATE( PLON(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_lonlat2xy: could not allocate ',
     +             IDM*JDM,' words for PLON'
        CALL EXIT(3)
      ENDIF
      ALLOCATE( PLAT(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_lonlat2xy: could not allocate ',
     +             IDM*JDM,' words for PLAT'
        CALL EXIT(3)
      ENDIF
      ALLOCATE( AMSK(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_lonlat2xy: could not allocate ',
     +             IDM*JDM,' words for AMSK'
        CALL EXIT(3)
      ENDIF
C
      CALL LONLAT(PLON,PLAT,AMSK,IDM,JDM,PAD,NPAD, XP,YP,L2IJ,
     &            CFILEA,CFILEM)
      CALL EXIT(0)
      END
      SUBROUTINE LONLAT(PLON,PLAT,AMSK,IDM,JDM,PAD,NPAD, XP,YP,L2IJ,
     &                  CFILEA,CFILEM)
      IMPLICIT NONE
C
      CHARACTER*240 CFILEA,CFILEM
      LOGICAL       L2IJ
      INTEGER       IDM,JDM,NPAD
      REAL*4        XP,YP
      REAL*4        PLON(IDM,JDM),PLAT(IDM,JDM),AMSK(IDM,JDM),PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C
      CHARACTER*240 CLINE
      CHARACTER*18  CASN
      REAL*8        X,Y
      REAL*4        DX,DY,DEG2RAD,DIST,DISTJ,DIST_MAX,D180,D360,QDX
      REAL*4        PLAT_MIN(JDM),PLAT_MAX(JDM)
      LOGICAL       L_IPJP,LPERIOD
      INTEGER       I,II,IP,J,JP,IOS,L,NRECL
#ifdef CRAY
      INTEGER*8     IU8,IOS8
#endif
C
      REAL*8         ZTECNF
      EXTERNAL       ZTECNF,ZTECNG,ZTECNP,ZTECNB
C
      REAL*8         A,B
      COMMON/ZAECNB/ A(0:2,0:2),B(0:2,0:2)
      SAVE  /ZAECNB/
C
      REAL*8         ACC,ERR,STEP
      REAL*8         X2(2),W(6),D(0:2,0:2)
      INTEGER        ITS
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
      CLOSE(UNIT=11)
C
      IF     (CFILEM.EQ." ") THEN
        AMSK(:,:) = 0.0
      ELSE
        OPEN(UNIT=11, FILE=CFILEM, FORM='UNFORMATTED', STATUS='OLD',
     +           ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
        IF     (IOS.NE.0) THEN
          write(6,*) 'Error: can''t open ',CFILEM(1:LEN_TRIM(CFILEM))
          write(6,*) 'ios   = ',ios
          write(6,*) 'nrecl = ',nrecl
          CALL EXIT(3)
        ENDIF
C
        READ(11,REC=1,IOSTAT=IOS) AMSK
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(AMSK,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'can''t read ',CFILEM(1:LEN_TRIM(CFILEM))
          CALL EXIT(4)
        ENDIF
        CLOSE(UNIT=11)
      ENDIF !cfilem:else
C
      LPERIOD = MAXVAL(PLON(:,:))-MINVAL(PLON(:,:)) .GT. 350.0
      DO J= 1,JDM
        IF     (.NOT.LPERIOD) THEN
          EXIT
        ENDIF
        DY = MOD( ABS(PLON(1,J) - PLON(  3,J)), 360.0 )
        IF     (DY.GT.180.0) THEN
          DY = 360.0 - DY  !abs distance
        ENDIF
        DX = MOD( ABS(PLON(1,J) - PLON(IDM,J)), 360.0 )
        IF     (DX.GT.180.0) THEN
          DX = 360.0 - DX  !abs distance
        ENDIF
        LPERIOD = LPERIOD .AND. DX.LT.DY  !1 and idm closer than 1 and 3
      ENDDO
C
C     FIND THE NEAREST POINT BY EXHAUSTIVE SEARCH (NOT EFFICIENT)
C     OPTIMIZED FOR MULTIPLE POINT CASE.
C
      D180    = 180.0
      D360    = 360.0
      DEG2RAD = 4.D0*ATAN(1.D0)/180.D0  !PI/180
C
      DO J= 1,JDM
        PLAT_MIN(J) = MINVAL(PLAT(:,J))
        PLAT_MAX(J) = MAXVAL(PLAT(:,J))
      ENDDO !j
      DIST_MAX = 0.0
      DO J= 1,JDM-1
        DO I= 1,IDM-1
          DIST_MAX = MAX( ABS(PLAT(I,J)-PLAT(I+1,J)),
     &                    ABS(PLAT(I,J)-PLAT(I,J+1)),
     &                    DIST_MAX )
        ENDDO !i
      ENDDO !j
      DIST_MAX = 2*DIST_MAX  !target must be at least this close in latitude
C
      IF     (YP.GT.-900.0) THEN  ! single lon,lat
        DIST   = 1.E30
        L_IPJP = .FALSE.
        DO J= 1,JDM
          DISTJ = MIN(DIST,DIST_MAX)
          IF     (YP.LT.PLAT_MIN(J)-DISTJ .OR.
     &            YP.GT.PLAT_MAX(J)+DISTJ     ) THEN
            CYCLE  ! far away row
          ENDIF
          L_IPJP = .TRUE.  !have a valid ip,jp
          IF     (DIST.EQ.0.0) THEN
            EXIT   ! found exact location
          ENDIF
          DO I= 1,IDM
            IF (AMSK(I,J).LT.2.0**99) THEN
              QDX = MAX(0.001,ABS(COS(YP*DEG2RAD)))
              DY  =      ABS(PLAT(I,J) - YP)
              DX  = MOD( ABS(PLON(I,J) - XP), D360 )
              IF     (DX.GT.D180) THEN
                DX = D360 - DX
              ENDIF
              IF     (QDX*DX+DY.LE.DIST) THEN
                IP   = I
                JP   = J
                DIST = QDX*DX+DY
              ENDIF
            ENDIF !amsk
          ENDDO !i
        ENDDO !j
        IF     (.NOT.L_IPJP) THEN
          IP = 0
          JP = 0
        ENDIF
C
        IF     (L2IJ) THEN
          WRITE(6,'(I6,I10)') IP,JP
        ELSEIF (JP.LE.1 .OR. JP.EQ.JDM) THEN
          X = IP
          Y = JP
          WRITE(6,'(2F10.3)') X,Y
        ELSE !2xy
C
C         FIND EXACT LOCATION WITH NAPACK ROUTINE(S).
C
          DO J= 0,2
            DO I= 0,2
              II = IP+I-1
              IF     (LPERIOD) THEN
                IF     (II.EQ.0) THEN
                  II = IDM
                ELSEIF (II.EQ.IDM+1) THEN
                  II = 1
                ENDIF   
              ENDIF
              B(I,J) =      PLAT(II,JP+J-1) - YP
              A(I,J) = MOD( PLON(II,JP+J-1) - XP, D360 )
              IF     (A(I,J).LT.-180.0) THEN
                A(I,J) = 360.0 + A(I,J)
              ELSEIF (A(I,J).GT. 180.0) THEN
                A(I,J) = A(I,J) - 360.0
              ENDIF
            ENDDO !i
          ENDDO !j
          STEP   = 0.0                        
          X2(1)  = 1.0
          X2(2)  = 1.0
          ACC    = 1.E-3
          CALL CG(X2,ERR,ITS,STEP,ACC,10,2,2,
     &            ZTECNF,ZTECNG,ZTECNB,ZTECNP,W)
          IF     (ITS.LT.0) THEN  !very flat extrema
            X2(1)  = 1.0
            X2(2)  = 1.0
          ELSEIF (MIN(X2(1),X2(2)).LT.-1.0 .OR.
     &            MAX(X2(1),X2(2)).GT. 3.0     ) THEN  !very bad CG result
            X2(1)  = 1.0
            X2(2)  = 1.0
          ENDIF
          X = IP + X2(1)-1.0
          Y = JP + X2(2)-1.0
          IF     (LPERIOD) THEN
            IF     (X.LT.1.0 .AND. X.GE.0.0) THEN
              X = X + IDM
            ELSEIF (X.GE.REAL(IDM+1) .AND. X.LT.REAL(IDM+2)) THEN
              X = X - IDM
            ENDIF
          ENDIF
          IF     (X.LT.1.0 .OR. X.GE.REAL(IDM+1) .OR.
     &            Y.LT.1.0 .OR. Y.GE.REAL(JDM+1)     ) THEN
            WRITE(6,*) "*** ",X,Y
            WRITE(6,'(2F9.3)') 0.0,0.0  !out of range
          ELSEIF (AMSK(MOD(NINT(X)-1,IDM)+1,
     &                 MOD(NINT(Y)-1,JDM)+1).GT.2.0**99) THEN
            WRITE(6,*) "*** ",X,Y
            WRITE(6,*) "*** ",MOD(NINT(X)-1,IDM)+1,
     &                        MOD(NINT(Y)-1,JDM)+1
            WRITE(6,'(2F9.3)') 0.0,0.0  !masked
          ELSE
            WRITE(6,'(2F9.3)') X,Y
          ENDIF
        ENDIF !2ij;2xy
C
      ELSE  !multiple points
        IP = IDM/2
        JP = JDM/2
        DO  !input loop
          READ(5,'(A)',IOSTAT=IOS) CLINE
          IF     (IOS.NE.0) THEN
            EXIT
          ENDIF
          IF     (CLINE     .EQ." "   .OR.
     &            CLINE(1:1).EQ."#"   .OR.
     &            CLINE(1:3).EQ.">>>" .OR.
     &            CLINE(1:3).EQ."+++" .OR.
     &            CLINE(1:3).EQ."***"     ) THEN
            WRITE(6,'(A)') TRIM(CLINE)
            CYCLE
          ENDIF
          READ(CLINE,*) XP,YP
C
C         USE THE LAST LOCATION AS A FIRST GUESS.
C
          QDX = MAX(0.001,ABS(COS(YP*DEG2RAD)))
          DY  =      ABS(PLAT(IP,JP) - YP)
          DX  = MOD( ABS(PLON(IP,JP) - XP), D360 )
          IF     (DX.GT.D180) THEN
            DX = D360 - DX
          ENDIF
          DIST = QDX*DX+DY
C
          L_IPJP = .FALSE.
          DO J= 1,JDM
            DISTJ = MIN(DIST,DIST_MAX)
            IF     (YP.LT.PLAT_MIN(J)-DISTJ .OR.
     &              YP.GT.PLAT_MAX(J)+DISTJ     ) THEN
              CYCLE  ! far away row
            ENDIF
            L_IPJP = .TRUE.  !have a valid ip,jp
            IF     (DIST.EQ.0.0) THEN
              EXIT   ! found exact location
            ENDIF
            DO I= 1,IDM
              IF (AMSK(I,J).LT.2.0**99) THEN
                DY =      ABS(PLAT(I,J) - YP)
                DX = MOD( ABS(PLON(I,J) - XP), D360 )
                IF     (DX.GT.D180) THEN
                  DX = D360 - DX
                ENDIF
                IF     (QDX*DX+DY.LE.DIST) THEN
                  IP   = I
                  JP   = J
                  DIST = QDX*DX+DY
                ENDIF
              ENDIF !amsk
            ENDDO !i
          ENDDO !j
          IF     (.NOT.L_IPJP) THEN  !not in the domain
            IP = 0
            JP = 0
          ENDIF
C
          IF     (L2IJ) THEN
            WRITE(6,'(I6,I10)') IP,JP
          ELSEIF (JP.LE.1 .OR. JP.EQ.JDM) THEN
            X = IP
            Y = JP
            WRITE(6,'(2F10.3)') X,Y
          ELSE !2xy
C
C           FIND EXACT LOCATION WITH NAPACK ROUTINE(S).
C
            DO J= 0,2
              DO I= 0,2
                II = IP+I-1
                IF     (LPERIOD) THEN
                  IF     (II.EQ.0) THEN
                    II = IDM
                  ELSEIF (II.EQ.IDM+1) THEN
                    II = 1
                  ENDIF   
                ENDIF
                B(I,J) =      PLAT(II,JP+J-1) - YP
                A(I,J) = MOD( PLON(II,JP+J-1) - XP, D360 )
                IF     (A(I,J).LT.-180.0) THEN
                  A(I,J) = 360.0 + A(I,J)
                ELSEIF (A(I,J).GT. 180.0) THEN
                  A(I,J) = A(I,J) - 360.0
                ENDIF
*               D(I,J) = SQRT( A(I,J)**2 +B(I,J)**2 )
              ENDDO !i
            ENDDO !j
*
*             WRITE(6,*) '***** A.2 = ',A(0:2,2)
*             WRITE(6,*) '***** A.1 = ',A(0:2,1)
*             WRITE(6,*) '***** A.0 = ',A(0:2,0)
*             WRITE(6,*) '***** B.2 = ',B(0:2,2)
*             WRITE(6,*) '***** B.1 = ',B(0:2,1)
*             WRITE(6,*) '***** B.0 = ',B(0:2,0)
*             WRITE(6,*) '***** D.2 = ',D(0:2,2)
*             WRITE(6,*) '***** D.1 = ',D(0:2,1)
*             WRITE(6,*) '***** D.0 = ',D(0:2,0)
*           
            STEP   = 0.0
            X2(1)  = 1.0
            X2(2)  = 1.0
            ACC    = 1.E-3
            CALL CG(X2,ERR,ITS,STEP,ACC,10,2,2,
     &              ZTECNF,ZTECNG,ZTECNB,ZTECNP,W)
C         
*           WRITE(6,*) '***** ITS = ',ITS
*           WRITE(6,*) '***** X2  = ',X2(:)
            IF     (ITS.LT.0) THEN  !very flat extrema
              X2(1)  = 1.0
              X2(2)  = 1.0
*             WRITE(6,*) '***** X2  = ',X2(:)
            ELSEIF (MIN(X2(1),X2(2)).LT.-1.0 .OR.
     &              MAX(X2(1),X2(2)).GT. 3.0     ) THEN  !very bad CG result
              X2(1)  = 1.0
              X2(2)  = 1.0
*             WRITE(6,*) '***** X2  = ',X2(:)
            ENDIF
            X = IP + X2(1)-1.0
            Y = JP + X2(2)-1.0
C
            IF     (LPERIOD) THEN
              IF     (X.LT.1.0 .AND. X.GE.0.0) THEN
                X = X + IDM
              ELSEIF (X.GE.REAL(IDM+1) .AND. X.LT.REAL(IDM+2)) THEN
                X = X - IDM
              ENDIF
            ENDIF
C
            IF     (X.LT.1.0 .OR. X.GE.REAL(IDM+1) .OR.
     &              Y.LT.1.0 .OR. Y.GE.REAL(JDM+1)     ) THEN
*             WRITE(6,*) "*** ",X,Y
              WRITE(6,'(2F10.3)') 0.0,0.0  !out of range
            ELSEIF (AMSK(MOD(NINT(X)-1,IDM)+1,
     &                   MOD(NINT(Y)-1,JDM)+1).GT.2.0**99) THEN
*             WRITE(6,*) "*** ",X,Y
*             WRITE(6,*) "*** ",MOD(NINT(X)-1,IDM)+1,
*    &                          MOD(NINT(Y)-1,JDM)+1
              WRITE(6,'(2F10.3)') 0.0,0.0  !masked
            ELSE
              WRITE(6,'(2F10.3)') X,Y
            ENDIF
            CALL FLUSH(6)
          ENDIF !2ij;2xy
        ENDDO  !input loop
      ENDIF
      RETURN
      END
C
C --- USER-LEVEL ROUTINES FOR NAPACK'S CG.
C
      REAL*8           FUNCTION ZTECNF(X)
      IMPLICIT NONE
C
      REAL*8           X(2)
C
C     WRAPPER FOR ZTECMB.
C
      REAL*8           F,G(2)
C
      CALL ZTECNB(F,G,X)
      ZTECNF = F
      RETURN
C     END OF ZTECNF.
      END
      SUBROUTINE ZTECNG(G,X)
      IMPLICIT NONE
C
      REAL*8           G(2),X(2)
C
C     WRAPPER FOR ZTECMB.
C
      REAL*8           F
C
      CALL ZTECNB(F,G,X)
      RETURN
C     END OF ZTECNG.
      END
      SUBROUTINE ZTECNP(Y,Z)
      IMPLICIT NONE
C
      REAL*8           Y(2),Z(2)
C
C     NULL PRECONDITIONER
C
      Y(1) = Z(1)
      Y(2) = Z(2)
      RETURN
C     END OF ZTECNP.
      END
      SUBROUTINE ZTECNB(F,G,X)
CFPP$ NOCONCUR R
      IMPLICIT NONE
C
      REAL*8         X(2),F,G(2)
C
      REAL*8         A,B
      COMMON/ZAECNB/ A(0:2,0:2),B(0:2,0:2)
      SAVE  /ZAECNB/
C
C**********
C*
C  1) CALCULATES FUNCTION (F) AND ITS GRADIENT (G) AT A POINT (X).
C
C  2) FUNCTION DEFINED IN [0.,2.]*[0.,2.] VIA BI-LINEAR FITS TO
C      A AND B (PASSED VIA /ZAECNB/) WITH THE RESULT ABS(A)+ABS(B).
C
C     THIS FUNCTION IS USED FOR COMPATIBILITY WITH BI-LINEAR
C      INTERPOLATION FROM ARRAY INDEX TO LON,LAT SPACE.
C
C  3) PASSED TO THE MINIMIZATION ROUTINE 'CG'.
C*
C**********
C
      INTEGER IP,JP
      REAL*8  D1,D2,DX,DY,FX(2),FY(2)
C
C     CHOOSE THE QUADRENT.
C
      IF     (X(1).GE.1.0) THEN
        IP = 1
      ELSE
        IP = 0
      ENDIF
      IF     (X(2).GE.1.0) THEN
        JP = 1
      ELSE
        JP = 0
      ENDIF
C
C     F  AT  X(1),X(2)
C
      DX = X(1)-IP
      DY = X(2)-JP
      D1 = (1.d0-DX)*(1.d0-DY)*A(IP,  JP  ) +
     &     (1.d0-DX)*      DY *A(IP,  JP+1) +
     &           DX *(1.d0-DY)*A(IP+1,JP  ) +
     &           DX *      DY *A(IP+1,JP+1)
      D2 = (1.d0-DX)*(1.d0-DY)*B(IP,  JP  ) +
     &     (1.d0-DX)*      DY *B(IP,  JP+1) +
     &           DX *(1.d0-DY)*B(IP+1,JP  ) +
     &           DX *      DY *B(IP+1,JP+1)
      F  = SQRT( D1**2 + D2**2 )
C
C     1ST DERIVATIVES.
C
      DX = X(1)-IP + 0.01
      DY = X(2)-JP
      D1 = (1.d0-DX)*(1.d0-DY)*A(IP,  JP  ) +
     &     (1.d0-DX)*      DY *A(IP,  JP+1) +
     &           DX *(1.d0-DY)*A(IP+1,JP  ) +
     &           DX *      DY *A(IP+1,JP+1)
      D2 = (1.d0-DX)*(1.d0-DY)*B(IP,  JP  ) +
     &     (1.d0-DX)*      DY *B(IP,  JP+1) +
     &           DX *(1.d0-DY)*B(IP+1,JP  ) +
     &           DX *      DY *B(IP+1,JP+1)
      FX(1) = SQRT( D1**2 + D2**2 )
C
      DX = X(1)-IP - 0.01
      DY = X(2)-JP
      D1 = (1.d0-DX)*(1.d0-DY)*A(IP,  JP  ) +
     &     (1.d0-DX)*      DY *A(IP,  JP+1) +
     &           DX *(1.d0-DY)*A(IP+1,JP  ) +
     &           DX *      DY *A(IP+1,JP+1)
      D2 = (1.d0-DX)*(1.d0-DY)*B(IP,  JP  ) +
     &     (1.d0-DX)*      DY *B(IP,  JP+1) +
     &           DX *(1.d0-DY)*B(IP+1,JP  ) +
     &           DX *      DY *B(IP+1,JP+1)
      FX(2) = SQRT( D1**2 + D2**2 )
C
      DX = X(1)-IP
      DY = X(2)-JP + 0.01
      D1 = (1.d0-DX)*(1.d0-DY)*A(IP,  JP  ) +
     &     (1.d0-DX)*      DY *A(IP,  JP+1) +
     &           DX *(1.d0-DY)*A(IP+1,JP  ) +
     &           DX *      DY *A(IP+1,JP+1)
      D2 = (1.d0-DX)*(1.d0-DY)*B(IP,  JP  ) +
     &     (1.d0-DX)*      DY *B(IP,  JP+1) +
     &           DX *(1.d0-DY)*B(IP+1,JP  ) +
     &           DX *      DY *B(IP+1,JP+1)
      FY(1) = SQRT( D1**2 + D2**2 )
C
      DX = X(1)-IP
      DY = X(2)-JP - 0.01
      D1 = (1.d0-DX)*(1.d0-DY)*A(IP,  JP  ) +
     &     (1.d0-DX)*      DY *A(IP,  JP+1) +
     &           DX *(1.d0-DY)*A(IP+1,JP  ) +
     &           DX *      DY *A(IP+1,JP+1)
      D2 = (1.d0-DX)*(1.d0-DY)*B(IP,  JP  ) +
     &     (1.d0-DX)*      DY *B(IP,  JP+1) +
     &           DX *(1.d0-DY)*B(IP+1,JP  ) +
     &           DX *      DY *B(IP+1,JP+1)
      FY(2) = SQRT( D1**2 + D2**2 )
C
      G(1) = (FX(1)-FX(2))/0.02
      G(2) = (FY(1)-FY(2))/0.02
*
*     WRITE(6,*) '***** X,  = ',X(1),X(2)
*     WRITE(6,*) '***** FX  = ',F,FX
*     WRITE(6,*) '***** FY  = ',F,FY
*     WRITE(6,*) '***** F,G = ',F,G(1),G(2)
      RETURN
C     END OF ZTECNB.
      END
C
C      ________________________________________________________
C     |                                                        |
C     |   MINIMIZE A FUNCTION USING THE FLETCHER-REEVES FORM   |
C     |            OF THE CONJUGATE GRADIENT METHOD            |
C     |            WITH (OR WITHOUT) PRECONDITIONING           |
C     |                                                        |
C     |    INPUT:                                              |
C     |                                                        |
C     |         X     --ARRAY CONTAINING STARTING GUESS        |
C     |                                                        |
C     |         STEP  --STARTING GUESS FOR MINIMIZER IN DIREC- |
C     |                 TION OF NEGATIVE GRADIENT DURING FIRST |
C     |                 ITERATION (E. G. STEP=1) WHEN STEP=0,  |
C     |                 THE PROGRAM SELECTS A STARTING GUESS   |
C     |                                                        |
*     |         T     --COMPUTING TOLERANCE (ITERATIONS STOP   |
*     |                 WHEN MAX-NORM OF GRADIENT .LE. T)      |
C     |         TT    --COMPUTING TOLERANCE (ITERATIONS STOP   |
C     |                 WHEN FUNCTION RESULT .LE. T)           |
C     |                                                        |
C     |         LIMIT --MAXIMUM NUMBER OF ITERATIONS           |
C     |                                                        |
C     |         N     --NUMBER OF UNKNOWNS                     |
C     |                                                        |
C     |         M     --NUMBER OF ITERATIONS UNTIL THE SEARCH  |
C     |                 DIRECTIONS ARE RENORMALIZED ALONG THE  |
C     |                 NEGATIVE GRADIENT (TYPICALLY, M = N)   |
C     |                                                        |
C     |         VALUE --NAME OF COST EVALUATION FUNC. ROUTINE  |
C     |                 (EXTERNAL IN MAIN PROGRAM)             |
C     |                 VALUE(X) IS VALUE OF COST AT X         |
C     |                                                        |
C     |         GRAD  --NAME OF GRADIENT EVALUATION SUBROUTINE |
C     |                 (EXTERNAL IN MAIN PROGRAM)             |
C     |                 GRAD(G,X) PUTS IN G THE GRADIENT AT X  |
C     |                                                        |
C     |         BOTH  --NAME SUBROUTINE TO EVALUATE BOTH COST  |
C     |                 AND ITS GRADIENT (EXTERNAL IN MAIN     |
C     |                 PROGRAM) BOTH(V,G,X) PUTS THE VALUE IN |
C     |                 V AND THE GRADIENT IN G FOR THE POINT X|
C     |                                                        |
C     |         PRE   --NAME OF PRECONDITIONING SUBROUTINE     |
C     |                 (EXTERNAL IN MAIN PROGRAM)             |
C     |                 PRE(Y,Z) APPLIES THE PRECONDITIONER TO |
C     |                 Z, STORING THE RESULT IN Y.            |
C     |                 IF PRECONDITIONING NOT USED SET Y = Z  |
C     |                                                        |
C     |         H     --WORK ARRAY (LENGTH AT LEAST 3N)        |
C     |                                                        |
C     |    OUTPUT:                                             |
C     |                                                        |
C     |         X     --MINIMIZER                              |
C     |                                                        |
*     |         E     --MAX-NORM OF GRADIENT                   |
C     |         EE    --FUNCTION RESULT
C     |                                                        |
C     |         IT    --NUMBER OF ITERATIONS PERFORMED         |
C     |                                                        |
C     |         STEP  --STEP SIZE ALONG SEARCH DIRECTION FOR   |
C     |                 FINAL ITERATION                        |
C     |                                                        |
C     |    BUILTIN FUNCTIONS: DABS,DEXP,IDINT,DLOG,DSQRT,DMAX1,|
C     |                         DMIN1,DSIGN                    |
C     |    PACKAGE ROUTINES: CUB,FD,FV,FVD,INS                 |
C     |________________________________________________________|
C
*     SUBROUTINE CG(X,E,IT,STEP,T,LIMIT,N,M,VALUE,GRAD,BOTH,PRE,H)
      SUBROUTINE CG(X,EE,IT,STEP,TT,LIMIT,N,M,VALUE,GRAD,BOTH,PRE,H)
      IMPLICIT NONE
      INTEGER I,IQ,IT,J,K,L,LIMIT,M,N,NA,NB,NC,ND
      REAL*8 H(N,*),X(*),Y(50),Z(50),A1,A2,A3,A4,A5,A6,A7,A8,A,B,C,C0,C1
      REAL*8 D,D0,DA,DB,E,F,F0,F1,FA,FB,FC,G,L3,P,Q,R,S,STEP,T,V,W
      REAL*8 TT,EE
      REAL*8 FV,FD,VALUE
      EXTERNAL BOTH,GRAD,PRE,VALUE
      DATA A1/.1D0/,A2/.9D0/,A3/5.D0/,A4/.2D0/,A5/10.D0/,A6/.9D0/
      DATA A7/.3D0/
      A8 = A3 + .01D0
      IT = 0
      CALL BOTH(F,H(1,3),X)
      E = 0.
      DO 10 I = 1,N
10         IF ( DABS(H(I,3)) .GT. E ) E = DABS(H(I,3))
*     IF ( E .LE. T ) RETURN
      EE = F
      IF (EE .LE. TT) RETURN
      L3 = 1./DLOG(A3)
      CALL PRE(H(1,2),H(1,3))
      A = STEP
      IF ( A .GT. 0. ) GOTO 30
      DO 20 I = 1,N
20         IF ( DABS(X(I)) .GT. A ) A = DABS(X(I))
      A = .01*A/E
      IF ( A .EQ. 0. ) A = 1.
30    G = 0.
      DO 40 I = 1,N
40         G = G + H(I,2)*H(I,3)
      IF ( G .LT. 0. ) GOTO 620
50    L = 0
      DO 60 I = 1,N
60         H(I,1) = -H(I,2)
      D = -G
70    FA = FV(A,X,H,N,VALUE)
      C0 = A
      F0 = FA
      J = 2
      Y(1) = 0.
      Z(1) = F
      Y(2) = A
      Z(2) = FA
      V = A1*D
      W = A2*D
      IQ = 0
      IF ( FA .LE. F ) GOTO 80
      C = A
      B = 0.
      A = 0.
      FC = FA
      FB = F
      FA = F
      GOTO 90
80    C = 0.
      B = 0.
      FC = F
      FB = F
      IQ = 1
90    NA = 0
      NB = 0
      NC = 0
      ND = 0
      Q = (D+(F-F0)/C0)/C0
      IF ( Q .LT. 0. ) GOTO 110
      Q = A
100   ND = ND + 1
      IF ( ND .GT. 25 ) GOTO 610
      Q = A3*Q
      P = FV(Q,X,H,N,VALUE)
      CALL INS(Q,P,A,B,C,FA,FB,FC,J,Y,Z)
      IF ( P-F .LT. W*Q ) GOTO 100
      GOTO 260
110   Q = .5*D/Q
      IF ( Q .LT. .01*C0 ) Q = .01*C0
      P = FV(Q,X,H,N,VALUE)
      IF ( P .LE. F0 ) GOTO 120
      F1 = F0
      C1 = C0
      F0 = P
      C0 = Q
      GOTO 130
120   F1 = P
      C1 = Q
130   CALL INS(Q,P,A,B,C,FA,FB,FC,J,Y,Z)
135   IF ( A .EQ. 0. ) GOTO 140
      IF ( FA-F .GE. V*A ) GOTO 160
      IF ( FA-F .LT. W*A ) GOTO 210
      GOTO 280
140   Q = C0
      IF ( C1 .LT. Q ) Q = C1
150   NA = NA + 1
      IF ( NA .GT. 25 ) GOTO 630
      Q = A4*Q
      P = FV(Q,X,H,N,VALUE)
      CALL INS(Q,P,A,B,C,FA,FB,FC,J,Y,Z)
      IF ( P-F .GE. V*Q ) GOTO 150
      GOTO 250
160   IF ( C0 .GT. C1 ) GOTO 200
      IF ( F0-F .GT. V*C0 ) GOTO 180
      IF ( F0-F .GE. W*C0 ) GOTO 320
      IF ( C1 .LE. A5*C0 ) GOTO 320
      R = DLOG(C1/C0)
      S = -IDINT(R*L3+.999)
      R = .999*DEXP(R/S)
      Q = C1
170   Q = Q*R
      IF ( Q .LT. C0 ) GOTO 320
      P = FV(Q,X,H,N,VALUE)
      CALL INS(Q,P,A,B,C,FA,FB,FC,J,Y,Z)
      NA = NA + 1
      IF ( P-F .GT. V*Q ) GOTO 170
      GOTO 320
180   Q = C0
190   NA = NA + 1
      IF ( NA .GT. 25 ) GOTO 630
      Q = A4*Q
      P = FV(Q,X,H,N,VALUE)
      CALL INS(Q,P,A,B,C,FA,FB,FC,J,Y,Z)
      IF ( P-F .GE. V*Q ) GOTO 190
      GOTO 250
200   Q = A
      GOTO 190
210   IF ( C0 .LT. C1 ) GOTO 290
      IF ( F0-F .GE. V*C0 ) GOTO 230
      IF ( F0-F .GE. W*C0 ) GOTO 250
      Q = C0
220   ND = ND  + 1
      IF ( ND .GT. 25 ) GOTO 610
      Q = A3*Q
      P = FV(Q,X,H,N,VALUE)
      CALL INS(Q,P,A,B,C,FA,FB,FC,J,Y,Z)
      IF ( P-F .LT. W*Q ) GOTO 220
      GOTO 250
230   IF ( C0 .LE. A5*C1 ) GOTO 250
      R = DLOG(C0/C1)
      S = IDINT(R*L3+.999)
      R = 1.001*DEXP(R/S)
      Q = A
240   Q = Q*R
      IF ( Q .GT. C0 ) GOTO 250
      ND = ND + 1
      P = FV(Q,X,H,N,VALUE)
      CALL INS(Q,P,A,B,C,FA,FB,FC,J,Y,Z)
      IF ( P-F .LT. W*Q ) GOTO 240
250   IF ( IQ .EQ. 1 ) GOTO 320
260   IF ( B .EQ. 0. ) GOTO 280
      IF ( C .EQ. 0. ) GOTO 270
      V = C - A
      W = A - B
      R = 1./V
      S = 1./W
      P = FC - FA
      Q = FB - FA
      E = P*R + Q*S
      IF ( DSIGN(E,C-B) .NE. E ) GOTO 320
      IF ( E .EQ. 0. ) GOTO 320
      Q = (P*R)*W - (Q*S)*V
      Q = A - .5*Q/E
      P = FV(Q,X,H,N,VALUE)
      CALL INS(Q,P,A,B,C,FA,FB,FC,J,Y,Z)
      GOTO 320
270   R = 1./A
      S = 1./B
      P = R*(FA-F) - D
      Q = S*(FB-F) - D
      E = A - B
      V = (R*P-S*Q)/E
      W = (A*Q*S-B*P*R)/E
      V = W*W-3.*V*D
      IF ( V .LT. 0. ) V = 0.
      V = DSQRT(V)
      IF ( W+V .EQ. 0. ) GOTO 320
      Q = -D/(W+V)
      IF ( Q .LE. 0. ) GOTO 320
      P = FV(Q,X,H,N,VALUE)
      CALL INS(Q,P,A,B,C,FA,FB,FC,J,Y,Z)
      GOTO 320
280   IF ( IQ .EQ. 1 ) GOTO  320
      Q = (D+(F-FA)/A)/A
      IF ( Q .GE. 0. ) GOTO 320
      Q = .5*D/Q
      P = FV(Q,X,H,N,VALUE)
      CALL INS(Q,P,A,B,C,FA,FB,FC,J,Y,Z)
      GOTO 320
290   IF ( F0-F .GT. V*C0 ) GOTO 300
      IF ( F0-F .GT. W*C0 ) GOTO 320
300   Q = A
310   ND = ND + 1
      IF ( ND .GT. 25 ) GOTO 610
      Q = A3*Q
      P = FV(Q,X,H,N,VALUE)
      CALL INS(Q,P,A,B,C,FA,FB,FC,J,Y,Z)
      IF ( P-F .LT. W*Q ) GOTO 310
      GOTO 250
320   DA = FD(A,X,H,N,GRAD)
      IF ( DA .GT. A6*G ) GOTO 410
      IF ( DA .GE. 0. ) GOTO 560
      R = A
      Q = 0.
      DO 330 I = 1,J
           IF ( Y(I) .GT. A ) GOTO 370
           IF ( Y(I) .LE. Q ) GOTO 330
           IF ( Y(I) .EQ. A ) GOTO 330
           Q = Y(I)
330   CONTINUE
      IF ( A .LE. A8*Q ) GOTO 560
      Q = A
340   ND = ND + 1
      IF ( ND .GT. 25 ) GOTO 610
      Q = A3*Q
      P = FV(Q,X,H,N,VALUE)
      F1 = FA
      CALL INS(Q,P,A,B,C,FA,FB,FC,J,Y,Z)
      IF ( P .LT. F1 ) GOTO 340
      IF ( A .GT. R ) GOTO 360
      DO 350 I = 1,N
350        H(I,2) = X(I) + A*H(I,1)
      GOTO 560
360   DA = FD(A,X,H,N,GRAD)
      IF ( DA .GT. A6*G ) GOTO 410
      GOTO 560
370   Q = Y(I)
      DO 380 K = I,J
           IF ( Y(K) .LE. A ) GOTO 380
           IF ( Y(K) .LT. Q ) Q = Y(K)
380   CONTINUE
      IF ( Q .LE. A5*A ) GOTO 560
      F0 = DLOG(Q/A)
      S = IDINT(F0*L3+.999)
      F0 = 1.001*DEXP(F0/S)
      S = A
390   S = S*F0
      IF ( S .GE. Q ) GOTO 320
      P = FV(S,X,H,N,VALUE)
      F1 = FA
      CALL INS(S,P,A,B,C,FA,FB,FC,J,Y,Z)
      IF ( P .LT. F1 ) GOTO 390
      IF ( A .GT. R ) GOTO 320
      DO 400 I = 1,N
400        H(I,2) = X(I) + A*H(I,1)
      GOTO 560
410   B = 0.
      K = 1
      I = K
420   I = I + 1
      IF ( I .GT. J ) GOTO 430
      IF ( Y(I) .GE. A ) GOTO 420
      IF ( Y(I) .LT. B ) GOTO 420
      B = Y(I)
      K = I
      GOTO 420
430   FB = Z(K)
      DB = D
      IF ( B .NE. 0. ) DB = FD(B,X,H,N,GRAD)
440   W = 2.*DABS(B-A)
      CALL CUB(C,A,B,FA,FB,DA,DB)
      NC = 1
      GOTO 480
450   W = .5*W
      IF ( W .LT. DABS(C0-C) ) GOTO 550
      IF ( C0 .LT. C ) GOTO 460
      IF ( D0 .GE. D ) GOTO 470
      GOTO 550
460   IF ( D0 .GT. D ) GOTO 550
470   CALL CUB(C,C,C0,F,F0,D,D0)
      NC = NC + 1
      IF ( NC .GT. 30 ) GOTO 600
480   R = DMAX1(A,B)
      S = DMIN1(A,B)
      IF ( C .GT. R ) GOTO 490
      IF ( C .GT. S ) GOTO 500
      C = S + (S-C)
      S = .5*(A+B)
      IF ( C .GT. S ) C = S
      GOTO 500
490   C = R - (C-R)
      S = .5*(A+B)
      IF ( C .LT. S ) C = S
500   C0 = A
      F0 = FA
      D0 = DA
      CALL FVD(F,D,C,X,H,N,BOTH)
      IF ( F .LT. FA ) GOTO 510
      B = C
      FB = F
      DB = D
      GOTO 450
510   IF ( C .LT. A ) GOTO 540
      IF ( D .LT. 0. ) GOTO 530
520   B = A
      FB = FA
      DB = DA
530   A = C
      FA = F
      DA = D
      IF ( D .GT. A6*G ) GOTO 450
      GOTO 560
540   IF ( D .LT. 0. ) GOTO 520
      GOTO 530
550   C = .5*(A+B)
      NB = NB + 1
      W = DABS(B-A)
      GOTO 500
560   E = 0.
      DO 570 I = 1,N
           IF ( DABS(H(I,3)) .GT. E ) E = DABS(H(I,3))
570        X(I) = H(I,2)
      IT = IT + 1
*     IF ( E .LE. T ) GOTO 660
      EE = F
      IF ( EE .LE. TT ) GOTO 660
      IF ( IT .GE. LIMIT ) GOTO 660
      F = FA
      D = DA
      A = A7*A
      CALL PRE(H(1,2),H(1,3))
      R = 0.
      DO 580 I = 1,N
580        R = R + H(I,2)*H(I,3)
      IF ( R .LT. 0. ) GOTO 620
      S = R/G
      G = R
      L = L + 1
      IF ( L .GE. M ) GOTO 50
      D = 0.
      DO 590 I = 1,N
           H(I,1) = -H(I,2) + S*H(I,1)
590        D = D + H(I,1)*H(I,3)
      GOTO 70
600   IF ( D .LT. G ) GOTO 560
*       WRITE(6,*) 'UNABLE TO OBTAIN DESCENT DIRECTION'
*       STOP
        IT = -1
        RETURN
610   CONTINUE
*       WRITE(6,*) 'THE FUNCTION DECREASES WITH NO MINIMUM'
*       STOP
        IT = -1
        RETURN
620   CONTINUE
*       WRITE(6,*) 'PRECONDITIONER NOT POSITIVE DEFINITE'
*       STOP
        IT = -1
        RETURN
630   CONTINUE
      Q = Q*A3**25
      ND = 0
640   ND = ND + 1
      IF ( ND .GT. 25 ) GOTO 650
      Q = A3*Q
      P = FV(Q,X,H,N,VALUE)
      CALL INS(Q,P,A,B,C,FA,FB,FC,J,Y,Z)
      IF ( P-F .GT. V*Q ) GOTO 640
      GOTO 135
650   CONTINUE
*       WRITE(6,*) 'UNABLE TO SATISFY ARMIJO CONDITION'
        IT = -1
        RETURN
660   CONTINUE
      STEP = A
      RETURN
      END
      REAL*8 FUNCTION FV(A,X,H,N,VALUE)
      REAL*8 H(N,*),X(*),A,VALUE
      EXTERNAL VALUE
      DO 10 I = 1 , N
10         H(I,2) = X(I) + A*H(I,1)
      FV = VALUE(H(1,2))
      RETURN
      END
      REAL*8 FUNCTION FD(A,X,H,N,GRAD)
      REAL*8 H(N,*),X(*),A,D
      EXTERNAL GRAD
      DO 10 I = 1 , N
10         H(I,2) = X(I) + A*H(I,1)
      CALL GRAD(H(1,3),H(1,2))
      D = 0.
      DO 20 I = 1,N
20         D = D + H(I,1)*H(I,3)
      FD = D
      RETURN
      END
      SUBROUTINE FVD(V,D,A,X,H,N,BOTH)
      IMPLICIT NONE
      INTEGER  N
      REAL*8   H(N,*),X(*),A,D,V
      EXTERNAL BOTH
      INTEGER  I
      DO 10 I = 1 , N
10         H(I,2) = X(I) + A*H(I,1)
      CALL BOTH(V,H(1,3),H(1,2))
      D = 0.
      DO 20 I = 1,N
20         D = D + H(I,1)*H(I,3)
      RETURN
      END
      SUBROUTINE CUB(X,A,B,C,D,E,F)
      IMPLICIT NONE
      REAL*8 A,B,C,D,E,F,G,V,W,X,Y,Z
      G = B - A
      IF ( G .EQ. 0. ) GOTO 50
      V = E + F - 3*(D-C)/G
      W = V*V-E*F
      IF ( W .LT. 0. ) W = 0.
      W = DSIGN(DSQRT(W),G)
      Y = E + V
      Z = F + V
      IF ( DSIGN(Y,G) .NE. Y ) GOTO 30
      IF ( DSIGN(Z,G) .NE. Z ) GOTO 20
      IF ( Z .EQ. 0. ) GOTO 20
10    X = B - G*F/(Z+W)
      RETURN
20    IF ( C .LT. D ) X = A
      IF ( C .GE. D ) X = B
      RETURN
30    IF ( DSIGN(Z,G) .NE. Z ) GOTO 40
      IF ( DABS(E) .GT. DABS(F) ) GOTO 10
40    X = A + G*E/(Y-W)
      RETURN
50    X = A
      RETURN
      END
      SUBROUTINE INS(S,F,A,B,C,FA,FB,FC,J,Y,Z)
      IMPLICIT NONE
      REAL*8 A,B,C,F,FA,FB,FC,S,Y(*),Z(*)
      INTEGER J
      J = J + 1
      Y(J) = S
      Z(J) = F
      IF ( F .LE. FA ) GOTO 20
      IF ( F .LE. FB ) GOTO 10
      IF ( F .GT. FC ) RETURN
      C = S
      FC = F
      RETURN
10    C = B
      B = S
      FC = FB
      FB = F
      RETURN
20    C = B
      B = A
      A = S
      FC = FB
      FB = FA
      FA = F
      RETURN
      END
