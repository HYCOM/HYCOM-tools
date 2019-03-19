      PROGRAM HYCOM_EDDY_CENTER
      IMPLICIT NONE
C
C  hycom_eddy_center - Usage:  hycom_eddy_center     ssh.a lon lat [grid.a]
C                              hycom_eddy_center_low ssh.a lon lat [grid.a]
C
C     Prints the location of the eddy center closest to lon lat in the
C     first SSH record and tracks this eddy center through subsequent SSH's.
C
C     Use hycom_eddy_center to track a "high" (maximum SSH) eddy and 
C     hycom_eddy_center_low to track a "low"  (minimum SSH) eddy.
C
C  grid.a is a hycom grid file, default regional.grid.a.
C  Note that the corresponding grid.b must also exist.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  November 2004.
C
      REAL*4, ALLOCATABLE :: PLAT(:,:),PLON(:,:),SSH(:,:)
      REAL*4              :: PAD(4096)
      INTEGER      IOS
      INTEGER      IARGC
      INTEGER      NARG
      CHARACTER*240 CARG
C
      LOGICAL       LDEBUG
      INTEGER       IDM,JDM,NPAD
      REAL*4        XP,YP,ETYPE
      CHARACTER*6   CVARIN
      CHARACTER*240 CFILES,CFILEG
C
C     READ ARGUMENTS.
C
      CALL GETARG(0,CARG)
      IF     (CARG.EQ.'hycom_eddy_center') THEN
        ETYPE  =  1.0
        LDEBUG = .FALSE.
      ELSEIF (CARG.EQ.'hycom_eddy_center_low') THEN
        ETYPE  = -1.0
        LDEBUG = .FALSE.
      ELSEIF (CARG.EQ.'hycom_eddy_center_debug') THEN
        ETYPE  =  1.0
        LDEBUG = .TRUE.
      ELSEIF (CARG.EQ.'hycom_eddy_center_low_debug') THEN
        ETYPE  = -1.0
        LDEBUG = .TRUE.
      ELSE
        WRITE(6,'(2a)')
     &    'Usage:  ',
     &    'hycom_eddy_center or hycom_eddy_center_low ...'
        CALL EXIT(1)
      ENDIF
C
      NARG = IARGC()
C
      IF     (NARG.EQ.3) THEN
        CALL GETARG(1,CFILES)
        CALL GETARG(2,CARG)
        READ(CARG,*) XP
        CALL GETARG(3,CARG)
        READ(CARG,*) YP
        CFILEG = 'regional.grid.a'
      ELSEIF (NARG.EQ.4) THEN
        CALL GETARG(1,CFILES)
        CALL GETARG(2,CARG)
        READ(CARG,*) XP
        CALL GETARG(3,CARG)
        READ(CARG,*) YP
        CALL GETARG(4,CFILEG)
      ELSE
        WRITE(6,*) 'Usage: hycom_eddy_center ssh.a lon lat [grid.a]'
        CALL EXIT(1)
      ENDIF
C
C     GET IDM,JDM FROM grid.b.
C
      OPEN(UNIT=11,FILE=CFILEG(1:LEN_TRIM(CFILEG)-1) // 'b',
     &     FORM='FORMATTED',STATUS='OLD',ACTION='READ')
C
      READ( 11,*) IDM,CVARIN
      IF (CVARIN.NE.'idm   ') THEN
        WRITE(6,*) 'hycom_eddy_center: bad header file ',
     &             CFILEG(1:LEN_TRIM(CFILEG)-1) // 'b'
        CALL EXIT(2)
      ENDIF
      READ( 11,*) JDM,CVARIN
      IF (CVARIN.NE.'jdm   ') THEN
        WRITE(6,*) 'hycom_eddy_center: bad header file ',
     &             CFILEG(1:LEN_TRIM(CFILEG)-1) // 'b'
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
      ALLOCATE(  SSH(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_eddy_center: could not allocate ',
     +             IDM*JDM,' words for  SSH'
        CALL EXIT(3)
      ENDIF
      ALLOCATE( PLON(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_eddy_center: could not allocate ',
     +             IDM*JDM,' words for PLON'
        CALL EXIT(3)
      ENDIF
      ALLOCATE( PLAT(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_eddy_center: could not allocate ',
     +             IDM*JDM,' words for PLAT'
        CALL EXIT(3)
      ENDIF
C
      CALL CENTER(SSH,PLON,PLAT,IDM,JDM,PAD,NPAD,
     &            XP,YP,ETYPE, LDEBUG, CFILES,CFILEG)
      CALL EXIT(0)
      END
      SUBROUTINE CENTER(SSH,PLON,PLAT,IDM,JDM,PAD,NPAD,
     &                  XP,YP,ETYPE, LDEBUG, CFILES,CFILEG)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILES,CFILEG
      LOGICAL       LDEBUG
      INTEGER       IDM,JDM,NPAD
      REAL*4        XP,YP,ETYPE
      REAL*4        SSH(IDM,JDM),PLON(IDM,JDM),PLAT(IDM,JDM),PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C
      REAL*8         ZTECNF
      EXTERNAL       ZTECNF,ZTECNG,ZTECNP,ZTECNB
C
      REAL*8         A
      COMMON/ZAECNB/ A(0:2,0:2)
C
      CHARACTER*18 CASN
      REAL*8       XC,YC
      REAL*8       ACC,ERR,STEP           
      REAL*8       X2(2),W(6)      
      REAL*4       SSHMAX,DX,DY
      INTEGER      I,IC,ICP1,II,ITS,J,JC,JJ,K,IOS,NRECL
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
      OPEN(UNIT=11, FILE=CFILEG, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEG)
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
        WRITE(6,*) 'can''t read ',TRIM(CFILEG)
        CALL EXIT(4)
      ENDIF
      READ(11,REC=2,IOSTAT=IOS) PLAT
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(PLAT,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read ',TRIM(CFILEG)
        CALL EXIT(4)
      ENDIF
C
C     FIND THE ARRAY COORDINATES OF THE STARTING POINT.
C
      CALL LONLAT2XY(PLON,PLAT,IDM,JDM, XP,YP, XC,YC)
C
      IF     (LDEBUG) THEN
        WRITE(6,*) '***** LONLAT2XY - XC,YC = ',XC,YC
      ENDIF
C
C     LOOP OVER ALL INPUT SSH ARRAYS.
C
      OPEN(UNIT=12, FILE=CFILES, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILES)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
      DO K= 1,999999
        READ(12,REC=K,IOSTAT=IOS) SSH
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(SSH,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          EXIT
        ENDIF
C
C       FIND ARRAY POINT NEAREST TO THE EDDY CENTER.
C
        II   = MAX( 2, MIN( IDM-1, INT(XC + 0.5D0) ) ) 
        JJ   = MAX( 2, MIN( JDM-1, INT(YC + 0.5D0) ) )
        IF     (SSH(II,JJ).EQ.SPVAL) THEN              
          SSHMAX = SPVAL
        ELSE                         
          SSHMAX = ETYPE*SSH(II,JJ)
        ENDIF                      
        IF     (LDEBUG) THEN
          WRITE(6,*) '***** NP - II,JJ,SSHMAX = ',
     +                           II,JJ,SSHMAX,ETYPE
        ENDIF
C
        DO !eddy-center
          IC = II         
          JC = JJ         
C              
          DO J= MAX(2,JC-2),MIN(JDM-1,JC+2)
            DO I= MAX(2,IC-2),MIN(IDM-1,IC+2)
              IF     (ETYPE*SSH(I,J).GT.SSHMAX .AND.
     +                      SSH(I,J).NE.SPVAL       ) THEN
                II     = I                                
                JJ     = J                                
                SSHMAX = ETYPE*SSH(I,J)
                IF     (LDEBUG) THEN
                  WRITE(6,*) '***** EC - II,JJ,SSHMAX = ',
     +                                   II,JJ,SSHMAX
                ENDIF
              ENDIF                    
            ENDDO !i
          ENDDO !j
          IF     (LDEBUG) THEN
            WRITE(6,*) '***** EC *****'
          ENDIF
C                 
          IF     (II.EQ.IC .AND. JJ.EQ.JC) THEN
            EXIT
          ENDIF                                
        ENDDO !eddy-center
C             
C       INITIALISE /ZATRAF/ AND CHECK FOR LAND.
C                                            
        DO J= 0,2                          
          DO I= 0,2
            IF     (SSH(IC-1+I,JC-1+J).EQ.SPVAL) THEN
              ETYPE = 0.0
            ELSEIF (SSHMAX.GT.0.0) THEN
              A(I,J) = -SSH(IC-1+I,JC-1+J)/SSH(IC,JC)
            ELSEIF (SSHMAX.LT.0.0) THEN
              A(I,J) =  SSH(IC-1+I,JC-1+J)/SSH(IC,JC)
            ELSEIF (SSHMAX.EQ.0.0) THEN
              ETYPE = 0.0  !safest option
            ENDIF                                    
          ENDDO !i
        ENDDO !j
C
        IF     (ETYPE.NE.0.0) THEN
C                                
C         FIND EDDY CENTER WITH NAPACK ROUTINE(S).
C                                               
          ACC    = 1.E-4*MAX( ABS(A(1,1) - A(1,0)),
     +                        ABS(A(1,1) - A(1,2)),
     +                        ABS(A(1,1) - A(0,1)),
     +                        ABS(A(1,1) - A(0,2)) )
C
          IF     (LDEBUG) THEN
            WRITE(6,*) '***** S.2 = ',SSH(IC-1:IC+1,JC+1)
            WRITE(6,*) '***** S.1 = ',SSH(IC-1:IC+1,JC  )
            WRITE(6,*) '***** S.0 = ',SSH(IC-1:IC+1,JC-1)
            WRITE(6,*) '***** A.2 = ',A(0:2,2)
            WRITE(6,*) '***** A.1 = ',A(0:2,1)
            WRITE(6,*) '***** A.0 = ',A(0:2,0)
            WRITE(6,*) '***** ACC = ',ACC
          ENDIF
C
          STEP   = 0.0
          X2(1)  = 1.0
          X2(2)  = 1.0
          CALL CG(X2,ERR,ITS,STEP,ACC,10,2,2,
     &            ZTECNF,ZTECNG,ZTECNB,ZTECNP,W)
C       
          IF     (ITS.LT.0) THEN  !very flat extrema
            X2(1)  = 1.0
            X2(2)  = 1.0
          ELSEIF (MIN(X2(1),X2(2)).LT.0.0 .OR.
     &            MAX(X2(1),X2(2)).GT.2.0     ) THEN  !bad CG result
            X2(1)  = 1.0
            X2(2)  = 1.0
          ENDIF
          XC = IC-1 + X2(1)
          YC = JC-1 + X2(2)
C
          IF     (LDEBUG) THEN
            IF     (X2(1).EQ.1.0 .AND. X2(2).EQ.1.0) THEN
              WRITE(6,*) '***** CG -   XC,YC = ',XC,YC
            ELSE
              WRITE(6,*) '***** CG - X,XC,YC = ',
     +                               X2(1),X2(2),XC,YC
            ENDIF
          ENDIF
        ELSEIF (MAX(ABS(XC-IC),ABS(YC-JC)).GT.1.0) THEN
          XC = IC  !best we can do
          YC = JC
          IF     (LDEBUG) THEN
            WRITE(6,*) '***** LAND - XC,YC = ',XC,YC
          ENDIF
        ELSE  !otherwise don't change xc,yc
          IF     (LDEBUG) THEN
            WRITE(6,*) '***** OLD  - XC,YC = ',XC,YC
          ENDIF
        ENDIF
C       
C       WRITE OUT EDDY CENTER (LAT-LON) POSITION.
C
        IC = XC
        JC = YC
        DX = XC - IC
        DY = YC - JC
        IF     (IC.EQ.IDM) THEN
          ICP1 = 1  !periodic  
        ELSE                   
          ICP1 = IC+1        
        ENDIF        
        XP = (1.0-DX)*(1.0-DY)*PLON(IC,  JC)   +
     +       (1.0-DX)*     DY *PLON(IC,  JC+1) +
     +            DX *(1.0-DY)*PLON(ICP1,JC)   +
     +            DX *     DY *PLON(ICP1,JC+1)  
        YP = (1.0-DX)*(1.0-DY)*PLAT(IC,  JC)   +
     +       (1.0-DX)*     DY *PLAT(IC,  JC+1) +
     +            DX *(1.0-DY)*PLAT(ICP1,JC)   +
     +            DX *     DY *PLAT(ICP1,JC+1)  
        WRITE(6,'(2F9.3)') XP,YP
      ENDDO !k
      CLOSE(UNIT=12)
C     END CENTER.
      END
      SUBROUTINE LONLAT2XY(PLON,PLAT,IDM,JDM, XP,YP, X,Y)
      IMPLICIT NONE
C
      INTEGER       IDM,JDM
      REAL*4        PLON(IDM,JDM),PLAT(IDM,JDM)
      REAL*4        XP,YP
      REAL*8        X,Y
C
C     CONVERT XP,YP (LON,LAT SPACE) TO X,Y (ARRAY INDEX SPACE).
C
      REAL*4       DX,DY,DIST,D180,D360
      REAL*4       PLAT_MIN(JDM),PLAT_MAX(JDM)
      LOGICAL      LPERIOD
      INTEGER      I,IP,J,JP
C
C     FIND THE NEAREST POINT BY EXHAUSTIVE SEARCH (NOT EFFICIENT)
C
      D180 = 180.0
      D360 = 360.0
C
      LPERIOD = MAXVAL(PLON(:,:))-MINVAL(PLON(:,:)) .GT. 350.0
C
        DIST = 1.E30
        DO J= 1,JDM
          DO I= 1,IDM
            DY =      ABS(PLAT(I,J) - YP)
            DX = MOD( ABS(PLON(I,J) - XP), D360 )
            IF     (DX.GT.D180) THEN
              DX = D360 - DX
            ENDIF
            IF     (DX+DY.LE.DIST) THEN
              IP   = I
              JP   = J
              DIST = DX+DY
            ENDIF
          ENDDO
        ENDDO
C
C       CONVERT NEAREST POINT TO EXACT LOCATION IN GRID INDEX SPACE.
C
        IF     (JP.EQ.JDM) THEN
          JP = JP-1
        ELSEIF (J.NE.1)    THEN
          DY =            (YP-PLAT(IP,JP))/
     &         (PLAT(IP,JP+1)-PLAT(IP,JP))
          IF     (DY.LT.0.0) THEN
            JP = JP-1
          ENDIF
        ENDIF
        IF     (I.EQ.IDM) THEN
          IF     (LPERIOD) THEN
            DX = MOD( PLON(IP,JP)-XP, 360.0 )
            IF     (DX.LT.-180.0) THEN
              DX = 360.0 + DX
            ELSEIF (DX.GT. 180.0) THEN
              DX = 360.0 - DX
            ENDIF
            DY = MOD( PLON(IP,JP)-PLON(   1,JP), 360.0 )
            IF     (DY.LT.-180.0) THEN
              DY = 360.0 + DY
            ELSEIF (DY.GT. 180.0) THEN
              DY = 360.0 - DY
            ENDIF
            IF     (ABS(DY).GT.0.001) THEN
              IF (DX/DY.LT.0.0) THEN
                IP = IP-1
              ENDIF
            ENDIF  !ABS(DY).GT.0.001
          ELSE
            IP = IP-1
          ENDIF  !LPERIOD:ELSE
        ELSE
          DX = MOD( PLON(IP,JP)-XP, 360.0 )
          IF     (DX.LT.-180.0) THEN
            DX = 360.0 + DX
          ELSEIF (DX.GT. 180.0) THEN
            DX = 360.0 - DX
          ENDIF
          DY = MOD( PLON(IP,JP)-PLON(IP+1,JP), 360.0 )
          IF     (DY.LT.-180.0) THEN
            DY = 360.0 + DY
          ELSEIF (DY.GT. 180.0) THEN
            DY = 360.0 - DY
          ENDIF
          IF     (ABS(DY).GT.0.001) THEN
            IF     (DX/DY.LT.0.0) THEN
              IF     (IP.NE.1) THEN
                IP = IP-1
              ELSEIF (LPERIOD) THEN
                IP = IDM
              ENDIF
            ENDIF  !DX/DY.LT.0.0
          ENDIF  !ABS(DY).GT.0.001
        ENDIF  !I.EQ.IDM:ELSE
C
        DX = MOD( PLON(IP,JP)-XP, 360.0 )
        IF     (DX.LT.-180.0) THEN
          DX = 360.0 + DX
        ELSEIF (DX.GT. 180.0) THEN
          DX = 360.0 - DX
        ENDIF
        DY = MOD( PLON(IP,JP)-PLON(IP+1,JP), 360.0 )
        IF     (DY.LT.-180.0) THEN
          DY = 360.0 + DY
        ELSEIF (DY.GT. 180.0) THEN
          DY = 360.0 - DY
        ENDIF
        IF     (ABS(DY).GT.0.001) THEN
          DX = DX/DY
        ELSE
          DX = 0.0
        ENDIF
        IF     (ABS(DX).LT.0.001) THEN
          DX = 0.0
        ELSEIF (ABS(DX-1.0).LT.0.001) THEN
          DX = 1.0
        ELSEIF (ABS(DX+1.0).LT.0.001) THEN
          DX = -1.0
        ENDIF
        X = IP
        X = X +DX
C
        DY =            (YP-PLAT(IP,JP))/
     &       (PLAT(IP,JP+1)-PLAT(IP,JP))
        IF     (ABS(DY).LT.0.001) THEN
          DY = 0.0
        ELSEIF (ABS(DY-1.0).LT.0.001) THEN
          DY = 1.0
        ELSEIF (ABS(DY+1.0).LT.0.001) THEN
          DY = -1.0
        ENDIF
        Y = JP
        Y = Y + DY
*
*       IF     (X.LT.1.0 .OR. X.GT.REAL(IDM) .OR.
*    &          Y.LT.1.0 .OR. Y.GT.REAL(JDM)     ) THEN
*         WRITE(6,'(2F9.3)') 0.0,0.0  !out of range
*       ELSE
*         WRITE(6,'(2F9.3)') X,Y
*       ENDIF
      RETURN
C     END LONLAT2XY.
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
      REAL*8           X(2),F,G(2)
C
      REAL*8           A
      COMMON/ZAECNB/   A(0:2,0:2)
C
C**********
C*
C  1) CALCULATES FUNCTION (F) AND ITS GRADIENT (G) AT A POINT (X).
C
C  2) FUNCTION DEFINED IN [0.,2.]*[0.,2.] VIA BI-QUADRATIC THAT FITS
C      A (PASSED VIA /ZAECNB/).
C
C  3) PASSED TO THE MINIMIZATION ROUTINE 'CG'.
C*
C**********
C
      INTEGER          I,J
      REAL*8           AX(0:2),AY(0:2)
C
C     F  AT  XX = X(1), YY = 0, 1, AND 2.
C
      DO 110 J= 0,2
        AX(J) = A(0,J) + X(1)*( (A(1,J) - A(0,J)) +
     +                          0.5*(X(1) - 1.0)*
     +                              (A(0,J) - 2.*A(1,J) + A(2,J)) )
  110 CONTINUE
C
C     F  AT  XX = 0, 1, AND 2, YY = X(2).
C
      DO 210 I= 0,2
        AY(I) = A(I,0) + X(2)*( (A(I,1) - A(I,0)) +
     +                          0.5*(X(2) - 1.0)*
     +                              (A(I,0) - 2.*A(I,1) + A(I,2)) )
  210 CONTINUE
C
C     F  AT  X(1),X(2)  (USING AX).
C
      F = AX(0) + X(2)*( (AX(1) - AX(0)) +
     +                   0.5*(X(2) - 1.0)*(AX(0) - 2.*AX(1) + AX(2)) )
C
C     1ST DERIVATIVES.
C
      G(1) = X(1)*(AY(0) - 2.*AY(1) + AY(2)) - 1.5*AY(0) + 
     +                                         2.0*AY(1) - 
     +                                         0.5*AY(2)
      G(2) = X(2)*(AX(0) - 2.*AX(1) + AX(2)) - 1.5*AX(0) + 
     +                                         2.0*AX(1) - 
     +                                         0.5*AX(2)
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
C     |         T     --COMPUTING TOLERANCE (ITERATIONS STOP   |
C     |                 WHEN MAX-NORM OF GRADIENT .LE. T)      |
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
C     |         E     --MAX-NORM OF GRADIENT                   |
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
      SUBROUTINE CG(X,E,IT,STEP,T,LIMIT,N,M,VALUE,GRAD,BOTH,PRE,H)
      IMPLICIT NONE
      INTEGER I,IQ,IT,J,K,L,LIMIT,M,N,NA,NB,NC,ND
      REAL*8 H(N,*),X(*),Y(50),Z(50),A1,A2,A3,A4,A5,A6,A7,A8,A,B,C,C0,C1
      REAL*8 D,D0,DA,DB,E,F,F0,F1,FA,FB,FC,G,L3,P,Q,R,S,STEP,T,V,W
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
      IF ( E .LE. T ) RETURN
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
      IF ( E .LE. T ) GOTO 660
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
