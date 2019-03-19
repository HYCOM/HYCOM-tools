      PROGRAM BOTTOM
      USE MOD_ZA  ! HYCOM array I/O interface
      IMPLICIT NONE
C
C     DIAGNOSTIC/DEBUGGING VARIABLES.
C
      INTEGER        ITEST,JTEST
      COMMON/DEBUGI/ ITEST,JTEST
      SAVE  /DEBUGI/
C
C     BLKDAT VARIABLES.
C
      CHARACTER*79 CTITLE
      CHARACTER*40 SIGFMT
      INTEGER      IVERSN,IEXPT,YRFLAG,MONTH,THFLAG,SIGVER
      INTEGER      JDW
      REAL*4       THBASE
C
C     I/O VARIABLES.
C
      CHARACTER PREAMBL(5)*79,CLINE*80
      REAL*4    HMINA,HMINB,HMAXA,HMAXB
C
C     INPUT CLIM ARRAYS.
C
      REAL*4               :: ZLEV(999)
      REAL*4,  ALLOCATABLE :: TZ(:,:,:),RZ(:,:,:)
C
C     OUTPUT ARRAYS
C
      REAL*4,  ALLOCATABLE :: PM(:,:),TM(:,:),SM(:,:),RM(:,:)
C
C     OTHER ARRAYS.
C
      INTEGER, ALLOCATABLE :: MSK(:,:)
      REAL*4,  ALLOCATABLE :: DEPTH(:,:),PLAT(:,:),
     +                        WORK(:,:)
C
C**********
C*
C 1)  FROM A Z-LEVEL CLIMATOLOGY ON THE HYCOM REGION GRID,
C      CREATE A DUMMY 1-LAYER 'ISOPYCNAL' CLIMATOLOGY
C      REPRESENTING THE BOTTOM VALUES.
C
C      ONLY FOR USE WITH HYCOM 2.0.00 OR LATER.
C
C 2)  PARAMETERS:
C
C     MODEL GRID SPECIFICATION (W.R.T. PRESSURE GRID):
C
C        IDM    = 1ST DIMENSION OF MAJOR (HYCOM) MODEL ARRAYS
C        JDM    = 2ND DIMENSION OF MAJOR (HYCOM) MODEL ARRAYS
C        KZ     = NUMBER OF Z-LEVELS IN INPUT CLIMATOLOGY
C
C 3)  INPUT:
C        ON UNIT 51:  BATHYMETRY FILE
C        ON UNIT 51A: BATHYMETRY FILE
C        ON UNIT 71:  DENS. Z-LEVEL CLIM FILE
C        ON UNIT 72:  TEMP. Z-LEVEL CLIM FILE
C        ON UNIT 99:  A SUBSET OF blkdat.input FROM TARGET SIMULATION
C     OUTPUT:
C        ON UNIT 21:  DUMMY HYCOM ARCHIVE FILE
C        ON UNIT 21A: DUMMY HYCOM ARCHIVE FILE
C
C 4)  ALAN J. WALLCRAFT,  NAVAL RESEARCH LABORATORY, MARCH 2001
C                                                 AND JUNE 2009.
C*
C**********
C
      REAL*4     ZERO,ONE,RADIAN
      PARAMETER (ZERO=0.0, ONE=1.0, RADIAN=57.2957795)
C
      INTEGER I,J,K,KK,KZ,KZTOP,L
      REAL*4  DP0KF,DPMS,DS0KF,DSMS,SZK,RZK,TZK,TIME,THK,THIKMN,
     +        RZLOC,DMIN,
     +        PAVE,XAVE,XMAX,XMIN,TMIN,SMIN,PMIN,ZJ,ZZ,Q,
     +        PINTEG,SIGMAA,SIGMAB,ZBOT,ZTOP
      REAL*4  SIG_V,SOFSIG_V
C
      CALL XCSPMD
C
      CALL ZHOPEN(71, 'FORMATTED', 'OLD', 0)
      READ(71,*) !5-line header
      READ(71,*)
      READ(71,*)
      READ(71,*)
      READ(71,*)
      DO K= 1,9999
        READ(71,*,END=100) !one line per level
      ENDDO
  100 CONTINUE
      KZ = K-1
      write(6,*) 'kz = ',kz
      CLOSE(71)
C
      ALLOCATE( TZ(KZ+1,IDM,JDM) )
      ALLOCATE( RZ(KZ+1,IDM,JDM) )
      ALLOCATE(      PM(IDM,JDM) )
      ALLOCATE(      TM(IDM,JDM) )
      ALLOCATE(      SM(IDM,JDM) )
      ALLOCATE(      RM(IDM,JDM) )
      ALLOCATE(   DEPTH(IDM,JDM) )
      ALLOCATE(    PLAT(IDM,JDM) )
      ALLOCATE(    WORK(IDM,JDM) )
      ALLOCATE(     MSK(IDM,JDM) )
C
C     BLKDAT INPUT.
C
      CALL ZHOPEN(99, 'FORMATTED', 'OLD', 0)
C
C --- 'ctitle' = climatology title
C
      WRITE(6,*)
      READ(99,'(A79)') CTITLE
      WRITE(6,'(A79)') CTITLE
C
C --- 'month'  = month (1 to 12)
C --- 'sigver' = version of the equation of state
C --- 'iversn' = hycom version number x10
C --- 'iexpt'  = experiment number x10
C --- 'yrflag' = days in year flag (0=360,1=366,2=366Jan1)
C --- 'idm   ' = longitudinal array size
C --- 'jdm   ' = latitudinal  array size
C --- 'jdw   ' = width of zonal average (optional, default 0)
C --- 'itest ' = grid point where detailed diagnostics are desired
C --- 'jtest ' = grid point where detailed diagnostics are desired
C
      WRITE(6,*)
      CALL BLKINI(MONTH, 'month ')
      CALL BLKINI(SIGVER, 'sigver')
      CALL BLKINI(IVERSN,'iversn')
      CALL BLKINI(IEXPT, 'iexpt ')
      CALL BLKINI(YRFLAG,'yrflag')
      WRITE(6,*)
      CALL BLKINI(I,     'idm   ')
      CALL BLKINI(J,     'jdm   ')
C
      IF     (I.NE.IDM) THEN
        WRITE(6,*)
        WRITE(6,*) 'ERROR - wrong IDM'
        WRITE(6,*)
        CALL ZHFLSH(6)
        STOP
      ELSEIF (J.NE.JDM) THEN
        WRITE(6,*)
        WRITE(6,*) 'ERROR - wrong JDM'
        WRITE(6,*)
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
      CALL BLKINI2(I,J, 'jdw   ','itest ')
      IF     (J.EQ.1) THEN !jdw
        JDW = I
        CALL BLKINI(ITEST, 'itest ')
      ELSE !itest
        ITEST = I
        JDW   = 0
      ENDIF
      CALL BLKINI(JTEST, 'jtest ')
C
C --- 'thflag' = reference pressure flag (0=Sigma-0,2=Sigma-2,4=Sigma-4)
      WRITE(6,*)
      CALL BLKINI(THFLAG,'thflag')
      IF     (THFLAG.EQ.0) THEN
        WRITE(SIGFMT,'(A,I2,A)')
     &    '(a6," =",f10.4," sigma-0 (sigver=',SIGVER,')")'
        IF     (MOD(SIGVER,2).NE.1 .OR. SIGVER.GE.40) THEN
          WRITE(6,*)
          WRITE(6,*) 'ERROR - thflag=0 not compatible with sigver =',
     &               SIGVER
          WRITE(6,*)
          CALL ZHFLSH(6)
          STOP
        ENDIF
      ELSEIF (THFLAG.EQ.2) THEN
        WRITE(SIGFMT,'(A,I2,A)')
     &    '(a6," =",f10.4," sigma-2 (sigver=',SIGVER,')")'
        IF     (MOD(SIGVER,2).NE.0 .OR. SIGVER.GE.40) THEN
          WRITE(6,*)
          WRITE(6,*) 'ERROR - thflag=2 not compatible with sigver =',
     &               SIGVER
          WRITE(6,*)
          CALL ZHFLSH(6)
          STOP
        ENDIF
      ELSEIF (THFLAG.EQ.4) THEN
        WRITE(SIGFMT,'(A,I2,A)')
     &    '(a6," =",f10.4," sigma-4 (sigver=',SIGVER,')")'
        IF     (SIGVER.LT.40) THEN
          WRITE(6,*)
          WRITE(6,*) 'ERROR - thflag=4 not compatible with sigver =',
     &               SIGVER
          WRITE(6,*)
          CALL ZHFLSH(6)
          STOP
        ENDIF
      ELSE
        WRITE(6,*)
        WRITE(6,*) 'ERROR - thflag must be 0 or 2 or 4'
        WRITE(6,*)
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
C --- 'thbase' = reference density (sigma units)
      WRITE(6,*)
      CALL BLKINR(THBASE,'thbase',SIGFMT)
      WRITE(6,*)
      CLOSE(UNIT=99)
C
C     TOPOGRAPHY INPUT.
C
      CALL ZAIOST
C
      CALL ZHOPEN(51, 'FORMATTED', 'OLD', 0)
      READ (51,'(A79)') PREAMBL
      READ (51,'(A)')   CLINE
      CLOSE(UNIT=51)
      WRITE(6,'(/(1X,A79))') PREAMBL,CLINE
C
      I = INDEX(CLINE,'=')
      READ (CLINE(I+1:),*)   HMINB,HMAXB
C
      CALL ZAIOPN('OLD', 51)
      CALL ZAIORD(DEPTH,MSK,.FALSE., HMINA,HMAXA, 51)
      CALL ZAIOCL(51)
C
      IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &        ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
        WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b topography files not consistent:',
     &    '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &    '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
C     INITIALIZE LAND MASK.
C
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (DEPTH(I,J).LT.2.0**99) THEN
            MSK(  I,J) = 1
          ELSE
            MSK(  I,J) = 0
            DEPTH(I,J) = ZERO
          ENDIF
        ENDDO
      ENDDO
C
C     CHECK ITEST,JTEST.
C
      IF     (MIN(ITEST,JTEST).GT.0) THEN
        IF     (ITEST.GT.IDM) THEN
          WRITE(6,'(/ a /)') 'error - itest > idm'
          CALL ZHFLSH(6)
          STOP
        ELSEIF (JTEST.GT.JDM) THEN
          WRITE(6,'(/ a /)') 'error - jtest > jdm'
          CALL ZHFLSH(6)
          STOP
        ELSEIF(MSK(ITEST,JTEST).EQ.0) THEN
          WRITE(6,'(/ a /)') 'error - itest,jtest is a land point'
          CALL ZHFLSH(6)
          STOP
        ENDIF
      ENDIF
C
C     LATITUDE GRID INPUT.
C
      CALL ZHOPNC(31, 'regional.grid.b', 'FORMATTED', 'OLD', 0)
      CALL ZAIOPF('regional.grid.a', 'OLD', 31)
C
      READ(31,*) ! skip idm
      READ(31,*) ! skip jdm
      READ(31,*) ! skip mapflg
      READ(31,*) ! skip plon
      CALL ZAIOSK(31)
      READ(31,'(A)') CLINE
      I = INDEX(CLINE,'=')
      READ (CLINE(I+1:),*)   HMINB,HMAXB
      CALL ZAIORD(PLAT,MSK,.FALSE., HMINA,HMAXA, 31)
      IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &        ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
        WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b grid files not consistent (plat):',
     &    '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &    '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
C     Z-LEVEL CLIMATOLOGY INPUT.
C
      CALL ZAIOPN('OLD', 71)
      CALL ZHOPEN(71, 'FORMATTED', 'OLD', 0)
      READ (71,'(A79)') PREAMBL
      WRITE(6, '(/(1X,A79))') PREAMBL
      DO K= 1,KZ
        READ (71,'(A)') CLINE
        I = INDEX(CLINE,'=')
        READ (CLINE(I+1:),*) ZLEV(K),HMINB,HMAXB
        IF     (HMINB.EQ.HMAXB) THEN  ! constant field
          DO J= 1,JDM
            DO I= 1,IDM
              RZ(K,I,J) = HMINB
            ENDDO
          ENDDO
          CALL ZAIOSK(71)
        ELSE
          CALL ZAIORD(WORK,MSK,.FALSE., HMINA,HMAXA, 71)
          IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &            ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
            WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &       'error - .a and .b density files not consistent:',
     &       '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &       '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
            CALL ZHFLSH(6)
            STOP
          ENDIF
          DO J= 1,JDM
            DO I= 1,IDM
              RZ(K,I,J) = WORK(I,J)
            ENDDO
          ENDDO
        ENDIF
*       if     (min(itest,jtest).gt.0) then
*         write(6,'(a,2i5,i3,a,f7.3)')
*    +     'i,j,k =',itest,jtest,k,
*    +     '   r =',rz(k,itest,jtest)
*       endif
C
      ENDDO
      CLOSE(UNIT=71)
      CALL ZAIOCL(71)
C
      CALL ZAIOPN('OLD', 72)
      CALL ZHOPEN(72, 'FORMATTED', 'OLD', 0)
      READ (72,'(A79)') PREAMBL
      WRITE(6, '(/(1X,A79))') PREAMBL
      DO K= 1,KZ
        READ (72,'(A)') CLINE
        I = INDEX(CLINE,'=')
        READ (CLINE(I+1:),*) ZLEV(K),HMINB,HMAXB
        IF     (HMINB.EQ.HMAXB) THEN  ! constant field
          DO J= 1,JDM
            DO I= 1,IDM
              WORK(I,J) = HMINB
            ENDDO
          ENDDO
          CALL ZAIOSK(72)
        ELSE
          CALL ZAIORD(WORK,MSK,.FALSE., HMINA,HMAXA, 72)
          IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &            ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
            WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &        'error - .a and .b temperature files not consistent:',
     &        '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &        '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
            CALL ZHFLSH(6)
            STOP
          ENDIF
        ENDIF
C
        DO J= 1,JDM
          DO I= 1,IDM
            TZ(K,I,J) = WORK(I,J)
          ENDDO
        ENDDO
*       if     (min(itest,jtest).gt.0) then
*         write(6,'(a,2i5,i3,a,f7.3)')
*    +     'i,j,k =',itest,jtest,k,
*    +     '   r =',rz(k,itest,jtest)
*         write(6,'(a,2i5,i3,a,f7.3)')
*    +     'i,j,k =',itest,jtest,k,
*    +     '   t =',tz(k,itest,jtest)
*         write(6,'(a,2i5,i3,a,f7.3)')
*    +     'i,j,k =',itest,jtest,k,
*    +     '   w =',work(itest,jtest)
*       endif
      ENDDO
      CLOSE (UNIT=72)
      CALL ZAIOCL(72)
C
C     ADD LEVEL KZ+1
C
      ZLEV(KZ+1) = ZLEV(KZ) + 10000.0
      DO J= 1,JDM
        DO I= 1,IDM
          RZ(KZ+1,I,J) = RZ(KZ,I,J) + 0.001
          TZ(KZ+1,I,J) = TZ(KZ,I,J)
        ENDDO
      ENDDO
C
C     DIAGNOSTIC PRINTOUT.
C
      IF     (MIN(ITEST,JTEST).GT.0) THEN
        WRITE(6,*)
        DO K= 1,KZ
          SZK = SOFSIG_V( RZ(K,ITEST,JTEST), TZ(K,ITEST,JTEST), SIGVER )
          WRITE(6,'(A,2I5,I3,A,F8.2,A,3F7.3)')
     +     'I,J,K =',ITEST,JTEST,K,
     +     '  ZLEV =',ZLEV(K),
     +     '   R,T,S =',RZ(K,ITEST,JTEST),
     +                  TZ(K,ITEST,JTEST),SZK
        ENDDO
      ENDIF
C
C     INITIALIZE CLIMATOLOGY OUTPUT.
C
      PREAMBL(1) = CTITLE
C
      WRITE(PREAMBL(2),4000) IEXPT/10,MOD(IEXPT,10)
C
      WRITE(PREAMBL(3),'(A,I1,A,I2,A)')
     +        'Layered averages w.r.t. Sigma-',THFLAG,
     +        ' (sigver=',SIGVER,')'
C
      WRITE(PREAMBL(5),'(A,2I5,I3,F9.3,F9.2,2F6.3)')
     +        'i/jdm =',
     +       IDM,JDM
C
      WRITE(6,*)
      PREAMBL(4) = 'Potential Temperature, or ' //
     +             'Salinity, or ' //
     +             'Interface Depths'
      WRITE(6, 4101) PREAMBL
      WRITE(6,*)
C
C     INITIALIZE DUMMY ARCHIVE OUTPUT.
C
      CALL ZAIOPN('NEW', 21)
      CALL ZHOPEN(21, 'FORMATTED', 'NEW', 0)
      YRFLAG = MAX(0,MIN(2,YRFLAG))
      WRITE(21,4200) MONTH,PREAMBL(1),PREAMBL(2),
     &               IVERSN,IEXPT,YRFLAG,IDM,JDM

      IF     (YRFLAG.EQ.0) THEN  ! 360 days, starting Jan 16
        TIME = (MONTH-1)*30.0
      ELSEIF (YRFLAG.EQ.1) THEN  ! 366 days, starting Jan 16
        TIME = (MONTH-1)*30.5
      ELSEIF (YRFLAG.EQ.2) THEN  ! 366 days, starting Jan 1
        TIME = (MONTH-1)*30.5 + 15.0
      ENDIF
C
C     FIRST INTERFACE IS AT SURFACE.
C     MIXED LAYER AVERAGE TEMPERATURE IS NOMINALLY SST-0.25 DEG.
C     THIS IS EQUIVALENT TO A 0.5 DEG CHANGE ACROSS THE MIXED
C     LAYER, WHICH IS LARGER THAN TYPICALLY USED FOR ACTUAL
C     PROFILES BECAUSE THIS IS A SMOOTH CLIMATOLOGICAL PROFILE.
C
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (DEPTH(I,J).GT.ZERO) THEN
            TM(  I,J) = TZ(1,I,J)-0.25
            SM(  I,J) = SOFSIG_V( RZ(1,I,J),
     +                            TZ(1,I,J), SIGVER )
            SIGMAA    = SIG_V(    TZ(1,I,J)-0.5,
     +                            SM(  I,J)    , SIGVER )
            CALL FIND_DENSITY(RZLOC,SIGMAA,RZ(1,I,J),KZ+1,1)
            IF     (RZLOC.EQ.0.0) THEN
              PM(I,J) = ZLEV(2)
            ELSEIF (RZLOC.EQ.KZ+1) THEN
              PM(I,J) = DEPTH(I,J)
            ELSE
              L = RZLOC
              Q = RZLOC - L
              PM(I,J) = (1.0-Q)*ZLEV(L) + Q*ZLEV(L+1) 
            ENDIF
            PM(I,J) = MIN( PM(I,J), DEPTH(I,J) )
              if (i.eq.itest .and. j.eq.jtest) then
                WRITE(6,'(A,2F10.3)')
     +            'PMIX,RZLOC =',PM(I,J),RZLOC
              endif
          ELSE
            PM(I,J) = ONE
            TM(I,J) = ZERO
            SM(I,J) = ZERO
            RM(I,J) = ZERO
C
            PM(I,J) = ZERO
C
            DEPTH(I,J) = ZERO  ! should already be zero
          ENDIF
        ENDDO
      ENDDO
C
C     PUT "SURFACE" FIELDS IN THE THE DUMMY ARCHIVE.
C
      DO J= 1,JDM
        DO I= 1,IDM
          WORK(I,J) = ZERO
        ENDDO
      ENDDO
      CALL ZAIOWR(WORK,MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'montg1  ',MONTH,TIME,0,ZERO,XMIN,XMAX
      CALL ZAIOWR(WORK,MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'srfhgt  ',MONTH,TIME,0,ZERO,XMIN,XMAX
      CALL ZAIOWR(WORK,MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'surflx  ',MONTH,TIME,0,ZERO,XMIN,XMAX
      CALL ZAIOWR(WORK,MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'salflx  ',MONTH,TIME,0,ZERO,XMIN,XMAX
C
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (DEPTH(I,J).GT.ZERO) THEN
            WORK(I,J) = PM(I,J)*9806.0
          ENDIF
        ENDDO
      ENDDO
      CALL ZAIOWR(WORK,MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'bl_dpth ',MONTH,TIME,0,ZERO,XMIN,XMAX
      CALL ZAIOWR(WORK,MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'mix_dpth',MONTH,TIME,0,ZERO,XMIN,XMAX
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (DEPTH(I,J).GT.ZERO) THEN
            TM(  I,J) = TZ(1,I,J)-0.25
            SM(  I,J) = SOFSIG_V( RZ(1,I,J),
     +                            TZ(1,I,J), SIGVER )
            WORK(I,J) = SIG_V(    TM(  I,J),
     +                            SM(  I,J), SIGVER ) - THBASE
          ELSE
            TM(  I,J) = ZERO
            SM(  I,J) = ZERO
            WORK(I,J) = ZERO
          ENDIF
        ENDDO
      ENDDO
      CALL ZAIOWR(TM,  MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'tmix    ',MONTH,TIME,0,  ZERO,XMIN,XMAX
      CALL ZAIOWR(SM,  MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'smix    ',MONTH,TIME,0,  ZERO,XMIN,XMAX
      CALL ZAIOWR(WORK,MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'thmix   ',MONTH,TIME,0,THBASE,XMIN,XMAX
C
      DO J= 1,JDM
        DO I= 1,IDM
          WORK(I,J) = ZERO
        ENDDO
      ENDDO
      CALL ZAIOWR(WORK,MSK,.FALSE., XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'umix    ',MONTH,TIME,0,ZERO,XMIN,XMAX
      CALL ZAIOWR(WORK,MSK,.FALSE., XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'vmix    ',MONTH,TIME,0,ZERO,XMIN,XMAX
      CALL ZAIOWR(WORK,MSK,.FALSE., XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'u_btrop ',MONTH,TIME,0,ZERO,XMIN,XMAX
      CALL ZAIOWR(WORK,MSK,.FALSE., XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'v_btrop ',MONTH,TIME,0,ZERO,XMIN,XMAX
C
C     ONLY LAYER IS AT DEEPEST TEMP AND SALINITY.
C
      K = 1
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (DEPTH(I,J).GT.ZERO) THEN
              PM(I,J) = DEPTH(I,J)
              CALL FIND_DEPTH(RZLOC,
     +                        PM(I,J), ZLEV,KZ+1,1)
              L = RZLOC
              Q = RZLOC - L
              RM(I,J) = (1.0-Q)*RZ(L,I,J) + Q*RZ(L+1,I,J)
              TM(I,J) = (1.0-Q)*TZ(L,I,J) + Q*TZ(L+1,I,J)
              SM(I,J) = SOFSIG_V(RM(I,J),TM(I,J),SIGVER)
            ENDIF  !DEPTH>0
          ENDDO  !J=1,JDM
        ENDDO  !I=1,IDM
C
C       STATISTICS
C
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (DEPTH(I,J).LE.ZERO) THEN
                RM(I,J) = RM(ITEST,JTEST)
              WORK(I,J) = 0.0
            ELSE
              WORK(I,J) = PM(I,J)
            ENDIF
          ENDDO
        ENDDO
        CALL LAYSTAT(TM,  WORK,IDM,JDM, TMIN,XAVE,XMAX)
        WRITE(6,8100) ' tem', TMIN,XAVE,XMAX, K,0.1
        CALL LAYSTAT(SM,  WORK,IDM,JDM, SMIN,XAVE,XMAX)
        WRITE(6,8100) ' sal', SMIN,XAVE,XMAX, K,0.1
        CALL LAYSTAT(RM,  WORK,IDM,JDM, XMIN,XAVE,XMAX)
        WRITE(6,8100) ' den', XMIN,XAVE,XMAX, K,0.1
        CALL LAYSTAT(PM,  WORK,IDM,JDM, PMIN,XAVE,XMAX)
        WRITE(6,8100) ' inf', PMIN,XAVE,XMAX, K,0.1
        CALL LAYSTAT(WORK, WORK,IDM,JDM, XMIN,XAVE,XMAX)
        WRITE(6,8100) ' thk', XMIN,XAVE,XMAX, K,0.1
C
C       DIAGNOSTIC PRINTOUT.
C
        IF     (MIN(ITEST,JTEST).GT.0) THEN
          WRITE(6,*)
          WRITE(6,'(A,2I5,I3,A,F9.2,F9.2,A,3F7.3)')
     +     'I,J,K =',ITEST,JTEST,K,
     +     '   INF,THK =',  PM(ITEST,JTEST),
     +                    WORK(ITEST,JTEST),
     +     '   R,T,S =',    RM(ITEST,JTEST),
     +                      TM(ITEST,JTEST),
     +                      SM(ITEST,JTEST)
          WRITE(6,*)
          CALL ZHFLSH(6)
        ENDIF
C
C       ZONAL INTERFACE DEPTHS.
C
        IF     (K.GT.1 .AND. JDW.NE.0) THEN
          XAVE = ZERO
          DO J= 1,JDM-1
            CALL LAYSTAT(PLAT(1,J), WORK(1,J),IDM,MIN(JDW,JDM-J),
     +                   XMIN,PAVE,XMAX)
            CALL LAYSTAT(PM(  1,J), WORK(1,J),IDM,MIN(JDW,JDM-J),
     +                   XMIN,XAVE,XMAX)
            WRITE(6,8200) K-1,J,MIN(J+JDW-1,JDM-1),
     +                    PAVE,XAVE
          ENDDO
        ENDIF
C
C       WRITE OUT DUMMY ARCHIVE.
C
        DO J= 1,JDM
          DO I= 1,IDM
            WORK(I,J) = ZERO
          ENDDO
        ENDDO
        CALL ZAIOWR(WORK,MSK,.FALSE., XMIN,XMAX, 21, .FALSE.)
        WRITE(21,4201) 'u-vel.  ',MONTH,TIME,K,0.1,XMIN,XMAX
        CALL ZAIOWR(WORK,MSK,.FALSE., XMIN,XMAX, 21, .FALSE.)
        WRITE(21,4201) 'v-vel.  ',MONTH,TIME,K,0.1,XMIN,XMAX
        DO J= 1,JDM
          DO I= 1,IDM
            WORK(I,J) = PM(I,J)*9806.0
          ENDDO
        ENDDO
        CALL ZAIOWR(WORK,MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
        WRITE(21,4201) 'thknss  ',MONTH,TIME,K,0.1,XMIN,XMAX
        CALL ZAIOWR(TM,  MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
        WRITE(21,4201) 'temp    ',MONTH,TIME,K,0.1,XMIN,XMAX
        CALL ZAIOWR(SM,  MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
        WRITE(21,4201) 'salin   ',MONTH,TIME,K,0.1,XMIN,XMAX
        DO J= 1,JDM
          DO I= 1,IDM
            WORK(I,J) = RM(I,J) - THBASE
          ENDDO
        ENDDO
        CALL ZAIOWR(WORK,MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
        WRITE(21,4201) 'density ',MONTH,TIME,K,0.1,XMIN,XMAX
C
        WRITE(6,6300) K,0.1
        CALL ZHFLSH(6)
      CLOSE (UNIT=21)
      CALL ZAIOCL(21)
      STOP
C
 4000 FORMAT('Expt ',I2.2,'.',I1.1)
 4101 FORMAT(A79)
 4102 FORMAT(A,': month,layer,dens,range = ',I2.2,I4.2,F7.3,1P2E16.7)
 4200 FORMAT(
     + 'Dummy HYCOM archive from climatology for month ',I2.2,'.' /
     + A80/A80/
     + '1234567890123456789012345678901234567890',
     + '1234567890123456789012345678901234567890'/
     & I5,4X,'''iversn'' = hycom version number x10'/
     & I5,4X,'''iexpt '' = experiment number x10'/
     & I5,4x,'''yrflag'' = days in year flag'/
     & I5,4x,'''idm   '' = longitudinal array size'/
     & I5,4x,'''jdm   '' = latitudinal  array size'/
     & 'field       time step  model day',
     & '  k  dens        min              max')
 4201 FORMAT(a8,' =',i11,f11.2,i3,f7.3,1p2e16.7)
 5000 FORMAT(A40)
 5500 FORMAT(6E13.6)
 6200 FORMAT(/ 1X,'MIN,MAX I COORDS = ',F8.2,',',F8.2 
     +       / 1X,'MIN,MAX J COORDS = ',F8.2,',',F8.2 /)
 6300 FORMAT(10X,'WRITING CLIM LAYER',I5,'    SIGMA =',F7.3 /)
 8100 FORMAT(1X,A,': min=',F9.2,' ave=',F9.2,' max=',F9.2,
     +   '   (k,sigma =',i3,F7.2,')')
 8200 FORMAT(' k =',I3.2,' j = ',I4.4,' to ',I4.4,
     +       ' lat =',F6.1,' inf =',F8.2)
C     END OF PROGRAM WNDINT.
      END
      SUBROUTINE LAYSTAT(PM, THICK, IDM,JDM, PMIN,PAVE,PMAX)
      IMPLICIT NONE
C
      INTEGER IDM,JDM
      REAL*4  PM(IDM,JDM),THICK(IDM,JDM), PMIN,PAVE,PMAX
C
C --- CALCULATE STATISTICS FOR PM.
C --- ONLY WHERE LAYER THICKNESS IS AT LEAST 10 CM.
C --- AVERAGE DOES NOT ALLOW FOR VARIATION IN GRID CELL SIZE.
C
      REAL*4     TENCM
      PARAMETER (TENCM=0.1)
C
      INTEGER I,J,IJSUM
      REAL*4  PPMIN,PPMAX
      REAL*8  PPSUM
C
      PPMIN =  1.E10
      PPMAX = -1.E10
      IJSUM =  0
      PPSUM =  0.0D0
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (THICK(I,J).GT.TENCM) THEN
            PPMIN = MIN( PPMIN, PM(I,J) )
            PPMAX = MAX( PPMAX, PM(I,J) )
            IJSUM = IJSUM + 1
            PPSUM = PPSUM + PM(I,J)
          ENDIF
        ENDDO
      ENDDO
      IF     (IJSUM.NE.0) THEN
        PMIN = PPMIN
        PMAX = PPMAX
        PAVE = PPSUM / IJSUM
      ELSE
        PMIN = 99.9
        PMAX = 99.9
        PAVE = 99.9
      ENDIF
      RETURN
C     END OF LAYSTAT.
      END
      SUBROUTINE FIND_DENSITY(DENLOC,DENTARG,RZ,KZ,MINZ)
      IMPLICIT NONE
C
      INTEGER KZ,MINZ
      REAL*4  DENLOC,DENTARG,RZ(KZ)
C
C     FIND EXACT LOCATION IN LAYER SPACE OF DENTARG
C     SEARCHING WITHIN RZ(MINZ:KZ)
C
C     ASSUME RZ IS MONOTONICALLY NON-DECREASING.
C
C     RETURN 0.0 IF DENTARG < RZ(MINZ)
C     RETURN KZ  IF DENTARG > RZ(KZ)
C
      INTEGER K
C
      IF     (DENTARG.LT.RZ(MINZ)) THEN
        DENLOC = 0.0
      ELSEIF (DENTARG.EQ.RZ(MINZ)) THEN
        DENLOC = MINZ
      ELSEIF (DENTARG.GE.RZ(KZ)) THEN
        DENLOC = KZ
      ELSE
        DO K= MINZ+1,KZ
          IF     (DENTARG.LE.RZ(K)) THEN !dentarg>rz(k-1)
            DENLOC = K-1 + (DENTARG-RZ(K-1))/(RZ(K)-RZ(K-1))
            EXIT
          ENDIF
        ENDDO !k
      ENDIF
      RETURN
      END
      SUBROUTINE FIND_DEPTH(ZLOC,ZTARG,Z,KZ,MINZ)
      IMPLICIT NONE
C
      INTEGER KZ,MINZ
      REAL*4  ZLOC,ZTARG,Z(KZ)
C
C     FIND EXACT LOCATION IN LAYER SPACE OF ZTARG
C     SEARCHING WITHIN Z(MINZ:KZ)
C
C     RETURN 0.0 IF ZTARG < Z(MINZ)
C     RETURN KZ  IF ZTARG > Z(KZ)
C
      INTEGER K
C
      IF     (ZTARG.LT.Z(MINZ)) THEN
        ZLOC = 0.0
      ELSEIF (ZTARG.EQ.Z(MINZ)) THEN
        ZLOC = MINZ
      ELSEIF (ZTARG.GE.Z(KZ)) THEN
        ZLOC = KZ
      ELSE
        DO K= MINZ+1,KZ
          IF     (ZTARG.LE.Z(K)) THEN !ztarg>z(k-1)
            ZLOC = K-1 + (ZTARG-Z(K-1))/(Z(K)-Z(K-1))
            EXIT
          ENDIF
        ENDDO !k
      ENDIF
      RETURN
      END
      REAL*4 FUNCTION SIG_V(TT,SS,SIGVER)
      IMPLICIT NONE
      INTEGER SIGVER
      REAL*4  TT,SS
C
C     SIGVER WRAPPER FOR SIG
C
      REAL*8 SS8,TT8
      REAL*8 SIG_1,SIG_2,SIG_3,SIG_4,SIG_5,SIG_6,SIG_7,SIG_8,
     &       SIG_46,SIG_48
C
      TT8 = TT
      SS8 = SS
      IF     (SIGVER.GT.40) THEN
        IF     (SIGVER.EQ.46) THEN
          SIG_V = SIG_46(TT8,SS8)
        ELSEIF (SIGVER.EQ.48) THEN
          SIG_V = SIG_48(TT8,SS8)
        ENDIF
      ELSEIF (MOD(SIGVER,2).EQ.1) THEN
        IF     (SIGVER.EQ.1) THEN
          SIG_V = SIG_1(TT8,SS8)
        ELSEIF (SIGVER.EQ.3) THEN
          SIG_V = SIG_3(TT8,SS8)
        ELSEIF (SIGVER.EQ.5) THEN
          SIG_V = SIG_5(TT8,SS8)
        ELSEIF (SIGVER.EQ.7) THEN
          SIG_V = SIG_7(TT8,SS8)
        ENDIF
      ELSE
        IF     (SIGVER.EQ.2) THEN
          SIG_V = SIG_2(TT8,SS8)
        ELSEIF (SIGVER.EQ.4) THEN
          SIG_V = SIG_4(TT8,SS8)
        ELSEIF (SIGVER.EQ.6) THEN
          SIG_V = SIG_6(TT8,SS8)
        ELSEIF (SIGVER.EQ.8) THEN
          SIG_V = SIG_8(TT8,SS8)
        ENDIF
      ENDIF
      RETURN
      END
      REAL*4 FUNCTION SOFSIG_V(RR,TT,SIGVER)
      IMPLICIT NONE
      INTEGER SIGVER
      REAL*4  RR,TT
C
C     SIGVER WRAPPER FOR SOFSIG
C
      REAL*8 RR8,TT8
      REAL*8 SOFSIG_1,SOFSIG_2,SOFSIG_3,SOFSIG_4,
     &       SOFSIG_5,SOFSIG_6,SOFSIG_7,SOFSIG_8,
     &       SOFSIG_46,SOFSIG_48
C
      RR8 = RR
      TT8 = TT
      IF     (SIGVER.GT.40) THEN
        IF     (SIGVER.EQ.46) THEN
          SOFSIG_V = SOFSIG_46(RR8,TT8)
        ELSEIF (SIGVER.EQ.48) THEN
          SOFSIG_V = SOFSIG_48(RR8,TT8)
        ENDIF
      ELSEIF (MOD(SIGVER,2).EQ.1) THEN
        IF     (SIGVER.EQ.1) THEN
          SOFSIG_V = SOFSIG_1(RR8,TT8)
        ELSEIF (SIGVER.EQ.3) THEN
          SOFSIG_V = SOFSIG_3(RR8,TT8)
        ELSEIF (SIGVER.EQ.5) THEN
          SOFSIG_V = SOFSIG_5(RR8,TT8)
        ELSEIF (SIGVER.EQ.7) THEN
          SOFSIG_V = SOFSIG_7(RR8,TT8)
        ENDIF
      ELSE
        IF     (SIGVER.EQ.2) THEN
          SOFSIG_V = SOFSIG_2(RR8,TT8)
        ELSEIF (SIGVER.EQ.4) THEN
          SOFSIG_V = SOFSIG_4(RR8,TT8)
        ELSEIF (SIGVER.EQ.6) THEN
          SOFSIG_V = SOFSIG_6(RR8,TT8)
        ELSEIF (SIGVER.EQ.8) THEN
          SOFSIG_V = SOFSIG_8(RR8,TT8)
        ENDIF
      ENDIF
      RETURN
      END
      REAL*8 FUNCTION SIG_1(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA0_7term.h'
      SIG_1 = SIG(TT8,SS8)
      END
      REAL*8 FUNCTION SOFSIG_1(RR8,TT8)
      IMPLICIT NONE
      REAL*8  RR8,TT8
      INCLUDE '../../include/stmt_fns_SIGMA0_7term.h'
      SOFSIG_1 = SOFSIG(RR8,TT8)
      END
      REAL*8 FUNCTION SIG_3(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA0_9term.h'
      SIG_3 = SIG(TT8,SS8)
      END
      REAL*8 FUNCTION SOFSIG_3(RR8,TT8)
      IMPLICIT NONE
      REAL*8  RR8,TT8
      INCLUDE '../../include/stmt_fns_SIGMA0_9term.h'
      SOFSIG_3 = SOFSIG(RR8,TT8)
      END
      REAL*8 FUNCTION SIG_5(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA0_17term.h'
      SIG_5 = SIG(TT8,SS8)
      END
      REAL*8 FUNCTION SOFSIG_5(RR8,TT8)
      IMPLICIT NONE
      REAL*8  RR8,TT8
      REAL*8, PARAMETER :: TOL=1.D-6
      INTEGER NN
      REAL*8  SN,SO
      REAL*8  SOFSIG_7
      INCLUDE '../../include/stmt_fns_SIGMA0_17term.h'
C     sofsig via Newton iteration from a 12-term 1st guess
      SN = SOFSIG_7(RR8,TT8)  !non-negative
      DO NN= 1,10
        SO = SN
        SN = SO - (SIG(TT8,SO)-RR8)/DSIGDS(TT8,SO)
        IF     (NN.EQ.10 .OR. ABS(SN-SO).LT.TOL) THEN
          EXIT
        ENDIF
      ENDDO !nn
      SOFSIG_5 = SN
      END
      REAL*8 FUNCTION SIG_7(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA0_12term.h'
      SIG_7 = SIG(TT8,SS8)
      END
      REAL*8 FUNCTION SOFSIG_7(RR8,TT8)
      IMPLICIT NONE
      REAL*8  RR8,TT8
      INCLUDE '../../include/stmt_fns_SIGMA0_12term.h'
      SOFSIG_7 = SOFSIG(RR8,TT8)
      END
      REAL*8 FUNCTION SIG_2(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA2_7term.h'
      SIG_2 = SIG(TT8,SS8)
      END
      REAL*8 FUNCTION SOFSIG_2(RR8,TT8)
      IMPLICIT NONE
      REAL*8  RR8,TT8
      INCLUDE '../../include/stmt_fns_SIGMA2_7term.h'
      SOFSIG_2 = SOFSIG(RR8,TT8)
      END
      REAL*8 FUNCTION SIG_4(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA2_9term.h'
      SIG_4 = SIG(TT8,SS8)
      END
      REAL*8 FUNCTION SOFSIG_4(RR8,TT8)
      IMPLICIT NONE
      REAL*8  RR8,TT8
      INCLUDE '../../include/stmt_fns_SIGMA2_9term.h'
      SOFSIG_4 = SOFSIG(RR8,TT8)
      END
      REAL*8 FUNCTION SIG_6(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA2_17term.h'
      SIG_6 = SIG(TT8,SS8)
      END
      REAL*8 FUNCTION SOFSIG_6(RR8,TT8)
      IMPLICIT NONE
      REAL*8  RR8,TT8
      REAL*8, PARAMETER :: TOL=1.D-6
      INTEGER NN
      REAL*8  SN,SO
      REAL*8  SOFSIG_8
      INCLUDE '../../include/stmt_fns_SIGMA2_17term.h'
C     sofsig via Newton iteration from a 12-term 1st guess
      SN = SOFSIG_8(RR8,TT8)  !non-negative
      DO NN= 1,10
        SO = SN
        SN = SO - (SIG(TT8,SO)-RR8)/DSIGDS(TT8,SO)
        IF     (NN.EQ.10 .OR. ABS(SN-SO).LT.TOL) THEN
          EXIT
        ENDIF
      ENDDO !nn
      SOFSIG_6 = SN
      END
      REAL*8 FUNCTION SIG_8(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA2_12term.h'
      SIG_8 = SIG(TT8,SS8)
      END
      REAL*8 FUNCTION SOFSIG_8(RR8,TT8)
      IMPLICIT NONE
      REAL*8  RR8,TT8
      INCLUDE '../../include/stmt_fns_SIGMA2_12term.h'
      SOFSIG_8 = SOFSIG(RR8,TT8)
      END
      REAL*8 FUNCTION SIG_46(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA4_17term.h'
      SIG_46 = SIG(TT8,SS8)
      END
      REAL*8 FUNCTION SOFSIG_46(RR8,TT8)
      IMPLICIT NONE
      REAL*8  RR8,TT8
      REAL*8, PARAMETER :: TOL=1.D-6
      INTEGER NN
      REAL*8  SN,SO
      REAL*8  SOFSIG_48
      INCLUDE '../../include/stmt_fns_SIGMA4_17term.h'
C     sofsig via Newton iteration from a 12-term 1st guess
      SN = SOFSIG_48(RR8,TT8)  !non-negative
      DO NN= 1,10
        SO = SN
        SN = SO - (SIG(TT8,SO)-RR8)/DSIGDS(TT8,SO)
        IF     (NN.EQ.10 .OR. ABS(SN-SO).LT.TOL) THEN
          EXIT
        ENDIF
      ENDDO !nn
      SOFSIG_46 = SN
      END
      REAL*8 FUNCTION SIG_48(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA4_12term.h'
      SIG_48 = SIG(TT8,SS8)
      END
      REAL*8 FUNCTION SOFSIG_48(RR8,TT8)
      IMPLICIT NONE
      REAL*8  RR8,TT8
      INCLUDE '../../include/stmt_fns_SIGMA4_12term.h'
      SOFSIG_48 = SOFSIG(RR8,TT8)
      END
