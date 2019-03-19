      PROGRAM PRIVER
      USE MOD_ZA  ! HYCOM array I/O interface
      IMPLICIT NONE
C
      INTEGER MXRIV
      PARAMETER (MXRIV=9000)
C
C     PRECIP ARRAY.
C
      INTEGER, ALLOCATABLE :: MSK(:,:)
      REAL*4,  ALLOCATABLE :: PLON(:,:),PLAT(:,:)
      REAL*4,  ALLOCATABLE :: DXDY(:,:)
      REAL*4,  ALLOCATABLE :: PCM(:,:),DEPTH(:,:)
C
      CHARACTER PREAMBL(5)*79,CLINE*80,CMOVE*11
C
      CHARACTER*79     CTITLE(3)
      NAMELIST/RTITLE/ CTITLE
      LOGICAL          IGNORE
      INTEGER          NRIVERS,NSMOOTH,MAXDIST,IJRIVER(2,MXRIV)
      REAL*8           TRIVER(12,MXRIV),TSCALE
      NAMELIST/RIVERS/ NRIVERS,NSMOOTH,MAXDIST,IJRIVER,TRIVER,TSCALE,
     &                 IGNORE
C
C**********
C*
C 1)  GENERATE A MONTHLY PRECIPITATION BOGAS FIELD
C      PARAMETERIZING RIVER INFLOW.
C
C 3)  NAMELIST INPUT:
C
C     /RTITLE/
C        CTITLE - THREE (79-CHARACTER) LINE TITLE.
C
C     /RIVERS/
C        NRIVERS - NUMBER OF RIVERS
C        NSMOOTH - NUMBER OF SMOOTHING PASSES.
C        MAXDIST - MAXIMUM NUMBER OF GRID POINTS TO MOVE RIVERS
C        IGNORE  - IGNORE "BAD" RIVERS (DEFAULT .FALSE.)
C        IJRIVER - I,J LOCATION OF EACH RIVER MOUTH
C        TRIVER  - MONTHLY RIVER INFLOW (UNITS OF SV/TSCALE)
C        TSCALE  - SCALE FACTOR TO CONVERT TRIVER TO SV (1.E6 M^3/S)
C
C 4)  INPUT:
C        ON UNIT  5: NAMELIST /RTITLE/, /RIVERN/, /RIVERS/
C        ON UNIT 51: HYCOM BATHYMETRY FILE (.[ab])
C     OUTPUT:
C        ON UNIT 14: HYCOM RIVER PRECIP. FILE (.[ab])
C*
C**********
C
      LOGICAL LARCTIC,LPERIOD,LFATAL
      INTEGER I,II,IJD,J,JJ,K,KR,MONTH
      REAL*4  X,Y,HMINA,HMINB,HMAXA,HMAXB,XMIN,XMAX,PLONWE
C
C     ARRAYS.
C
      CALL XCSPMD
      ALLOCATE(   MSK(IDM,JDM) )
      ALLOCATE(  PLON(IDM,JDM) )
      ALLOCATE(  PLAT(IDM,JDM) )
      ALLOCATE(  DXDY(IDM,JDM) )
      ALLOCATE(   PCM(IDM,JDM) )
      ALLOCATE( DEPTH(IDM,JDM) )
C
C     NAMELIST INPUT.
C
      CALL ZHOPEN(6, 'FORMATTED', 'UNKNOWN', 0)
C
      IGNORE    = .FALSE.
      CTITLE(1) = ' '
      CTITLE(2) = ' '
      CTITLE(3) = ' '
      WRITE(6,*) 'READING /RTITLE/'
      CALL ZHFLSH(6)
      READ( 5,RTITLE)
      WRITE(6,RTITLE)
C
      NRIVERS      = 0
      NSMOOTH      = 0
      MAXDIST      = 0
      IJRIVER(:,:) = 0
      TRIVER(:,:)  = 0.0
      TSCALE       = 1.0  ! default is TRIVER in SV
      WRITE(6,*) 'READING /RIVERS/'
      CALL ZHFLSH(6)
      READ( 5,RIVERS)
      WRITE(6,RIVERS)
      WRITE(6,*)
      CALL ZHFLSH(6)
C
      IF     (NRIVERS.LE.0 .OR. NRIVERS.GT.MXRIV) THEN
        WRITE(6,*) 'ERROR - NRIVERS OUT OF RANGE'
        WRITE(6,*)
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
C     GRID INPUT.
C
      CALL ZAIOST
C
      CALL ZHOPNC(21, 'regional.grid.b', 'FORMATTED', 'OLD', 0)
      CALL ZAIOPF('regional.grid.a', 'OLD', 21)
C
      READ(21,*) ! skip idm
      READ(21,*) ! skip jdm
      READ(21,*) ! skip mapflg
      READ(21,'(A)') CLINE
      I = INDEX(CLINE,'=')
      READ (CLINE(I+1:),*)   HMINB,HMAXB
      CALL ZAIORD(PLON,MSK,.FALSE., HMINA,HMAXA, 21)
      IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &        ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
        WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b grid files not consistent (plon):',
     &    '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &    '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
        CALL ZHFLSH(6)
        STOP
      ENDIF
      LPERIOD = HMAXB-HMINB .GT. 350.0
C
      READ(21,'(A)') CLINE
      I = INDEX(CLINE,'=')
      READ (CLINE(I+1:),*)   HMINB,HMAXB
      CALL ZAIORD(PLAT,MSK,.FALSE., HMINA,HMAXA, 21)
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
      READ(21,*) ! skip qlon
      CALL ZAIOSK(21)
      READ(21,*) ! skip qlat
      CALL ZAIOSK(21)
      READ(21,*) ! skip ulon
      CALL ZAIOSK(21)
      READ(21,*) ! skip ulat
      CALL ZAIOSK(21)
      READ(21,*) ! skip vlon
      CALL ZAIOSK(21)
      READ(21,*) ! skip vlat
      CALL ZAIOSK(21)
      READ(21,*) ! skip pang
      CALL ZAIOSK(21)
      READ(21,'(A)') CLINE
      I = INDEX(CLINE,'=')
      READ (CLINE(I+1:),*)   HMINB,HMAXB
      CALL ZAIORD(PCM,  MSK,.FALSE., HMINA,HMAXA, 21)
      IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &        ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
        WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b grid files not consistent (pscx):',
     &    '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &    '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
      READ(21,'(A)') CLINE
      I = INDEX(CLINE,'=')
      READ (CLINE(I+1:),*)   HMINB,HMAXB
      CALL ZAIORD(DXDY,  MSK,.FALSE., HMINA,HMAXA, 21)
      IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &        ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
        WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b grid files not consistent (pscy):',
     &    '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &    '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
        CALL ZHFLSH(6)
        STOP
      ENDIF
      DXDY(:,:) = PCM(:,:)*DXDY(:,:)  !dx*dy
C
      CLOSE(UNIT=21)
      CALL ZAIOCL(21)
C
C     CHECK FOR AN ARCTIC BIPOLAR PATCH
C
      LARCTIC = PLON(3*IDM/4+1,JDM-1).EQ.PLON(IDM/4,JDM) .AND.
     &          PLAT(3*IDM/4+1,JDM-1).EQ.PLAT(IDM/4,JDM)
C
      WRITE(6,*)
      WRITE(6,*) 'LPERIOD = ',LPERIOD
      WRITE(6,*) 'LARCTIC = ',LARCTIC
      WRITE(6,*)
      CALL ZHFLSH(6)
C
C     TOPOGRAPHY INPUT.
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
C     INITIALIZE LAND/SEA MASK.
C
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (DEPTH(I,J).GT.  0.0 .AND.
     &            DEPTH(I,J).LE.HMAXA      ) THEN
            MSK(I,J) = 1
          ELSE
            MSK(I,J) = 0
          ENDIF
        ENDDO
      ENDDO
C
C     INITIALIZE OUTPUT.
C
      CALL ZAIOPN('NEW', 14)
C
      CALL ZHOPEN(14, 'FORMATTED', 'NEW', 0)
C
      PREAMBL(1) = 
     +   'Monthly river inflow bogused to monthly precipitation (m/s)'
      PREAMBL(2) = CTITLE(1)
      PREAMBL(3) = CTITLE(2)
      PREAMBL(4) = CTITLE(3)
      WRITE(PREAMBL(5),'(I4,A,2I6)')
     +        NSMOOTH,' smoothing passes.  i/jdm =',IDM,JDM
C
      WRITE(14,4101) PREAMBL
      WRITE(6,*)
      WRITE(6, 4101) PREAMBL
      WRITE(6,*)
      CALL ZHFLSH(6)
C
C     LOOP THROUGH ALL MONTHS.
C
      DO MONTH= 1,12
C
C       ADD RIVER POINT SOURCES
C
        DO J= 1,JDM
          DO I= 1,IDM
            PCM(I,J) = 0.0
          ENDDO
        ENDDO
C
        LFATAL = .FALSE.
        DO KR= 1,NRIVERS
          I = IJRIVER(1,KR)
          J = IJRIVER(2,KR)
C
          IF     (MONTH.EQ.1) THEN
            CMOVE = '           '
            IF     (MSK(I,J).NE.1) THEN
C
C             MOVE RIVER TO CLOSE SEA POINT, IF THERE IS ONE.
C
              IJD = MAXDIST**2 + 1
              DO K= 1,MAXDIST
                IF     (IJD.LT.K**2) THEN
                  EXIT
                ENDIF  !early exit
                DO II= MAX(I-K,1),MIN(I+K,IDM)
                  JJ=MAX(J-K,1)
                    IF     (MSK(II,JJ).EQ.1) THEN
                      IF     ((I-II)**2+(J-JJ)**2.LT.IJD) THEN
                        IJD = (I-II)**2+(J-JJ)**2  !distance squared to I,J
                        IJRIVER(1,KR) = II
                        IJRIVER(2,KR) = JJ
                      ENDIF
                    ENDIF
                  JJ=MIN(J+K,JDM)
                    IF     (MSK(II,JJ).EQ.1) THEN
                      IF     ((I-II)**2+(J-JJ)**2.LT.IJD) THEN
                        IJD = (I-II)**2+(J-JJ)**2  !distance squared to I,J
                        IJRIVER(1,KR) = II
                        IJRIVER(2,KR) = JJ
                      ENDIF
                    ENDIF
                ENDDO !ii
                DO JJ= MAX(J-K+1,1),MIN(J+K-1,JDM)
                  II=MAX(I-K,1)
                    IF     (MSK(II,JJ).EQ.1) THEN
                      IF     ((I-II)**2+(J-JJ)**2.LT.IJD) THEN
                        IJD = (I-II)**2+(J-JJ)**2  !distance squared to I,J
                        IJRIVER(1,KR) = II
                        IJRIVER(2,KR) = JJ
                      ENDIF
                    ENDIF
                  II=MIN(I+K,IDM)
                    IF     (MSK(II,JJ).EQ.1) THEN
                      IF     ((I-II)**2+(J-JJ)**2.LT.IJD) THEN
                        IJD = (I-II)**2+(J-JJ)**2  !distance squared to I,J
                        IJRIVER(1,KR) = II
                        IJRIVER(2,KR) = JJ
                      ENDIF
                    ENDIF
                ENDDO !jj
              ENDDO !k
              IF     (IJD.NE.MAXDIST**2+1) THEN
*               write(6,'(a,4i5)')
*    &            'lnd,new ij = ',I,J,IJRIVER(1,KR),IJRIVER(2,KR)
                II  = IJRIVER(1,KR)
                JJ  = IJRIVER(2,KR)
                IJD = (I-II)**2+(J-JJ)**2  !distance squared to I,J
                I  = II
                J  = JJ
                WRITE(CMOVE,'(a,f6.1)') '(new)',SQRT(REAL(IJD))
              ELSE
                CMOVE = '(bad)'
              ENDIF
            ENDIF
            IF     (MSK(I,J).EQ.1 .AND.
     &              MINVAL(MSK(MAX(I-1,1):MIN(I+1,IDM),
     &                         MAX(J-1,1):MIN(J+1,JDM) )).NE.0) THEN
C
C             MOVE RIVER CLOSER TO LAND, IF POSSIBLE
C
              IJD = MAXDIST**2 + 1
              DO K= 1,MAXDIST
                IF     (IJD.LT.K**2) THEN
                  EXIT
                ENDIF  !early exit
                DO II= MAX(I-K,1),MIN(I+K,IDM)
                  JJ=MAX(J-K,1)
                    IF     (MSK(II,JJ).EQ.1 .AND.
     +                      MINVAL(MSK(MAX(II-1,1):MIN(II+1,IDM),
     &                                 MAX(JJ-1,1):MIN(JJ+1,JDM) ))
     &                      .EQ.0) THEN
                      IF     ((I-II)**2+(J-JJ)**2.LT.IJD) THEN
                        IJD = (I-II)**2+(J-JJ)**2  !distance squared to I,J
                        IJRIVER(1,KR) = II
                        IJRIVER(2,KR) = JJ
                      ENDIF
                    ENDIF
                  JJ=MIN(J+K,JDM)
                    IF     (MSK(II,JJ).EQ.1 .AND.
     +                      MINVAL(MSK(MAX(II-1,1):MIN(II+1,IDM),
     &                                 MAX(JJ-1,1):MIN(JJ+1,JDM) ))
     &                      .EQ.0) THEN
                      IF     ((I-II)**2+(J-JJ)**2.LT.IJD) THEN
                        IJD = (I-II)**2+(J-JJ)**2  !distance squared to I,J
                        IJRIVER(1,KR) = II
                        IJRIVER(2,KR) = JJ
                      ENDIF
                    ENDIF
                ENDDO !ii
                DO JJ= MAX(J-K+1,1),MIN(J+K-1,JDM)
                  II=MAX(I-K,1)
                    IF     (MSK(II,JJ).EQ.1 .AND.
     +                      MINVAL(MSK(MAX(II-1,1):MIN(II+1,IDM),
     &                                 MAX(JJ-1,1):MIN(JJ+1,JDM) ))
     &                      .EQ.0) THEN
                      IF     ((I-II)**2+(J-JJ)**2.LT.IJD) THEN
                        IJD = (I-II)**2+(J-JJ)**2  !distance squared to I,J
                        IJRIVER(1,KR) = II
                        IJRIVER(2,KR) = JJ
                      ENDIF
                    ENDIF
                  II=MIN(I+K,IDM)
                    IF     (MSK(II,JJ).EQ.1 .AND.
     +                      MINVAL(MSK(MAX(II-1,1):MIN(II+1,IDM),
     &                                 MAX(JJ-1,1):MIN(JJ+1,JDM) ))
     &                      .EQ.0) THEN
                      IF     ((I-II)**2+(J-JJ)**2.LT.IJD) THEN
                        IJD = (I-II)**2+(J-JJ)**2  !distance squared to I,J
                        IJRIVER(1,KR) = II
                        IJRIVER(2,KR) = JJ
                      ENDIF
                    ENDIF
                ENDDO !jj
              ENDDO !k
              IF     (IJD.NE.MAXDIST**2+1) THEN
*               write(6,'(a,4i5)')
*    &            'sea,new ij = ',I,J,IJRIVER(1,KR),IJRIVER(2,KR)
                II  = IJRIVER(1,KR)
                JJ  = IJRIVER(2,KR)
                IJD = (I-II)**2+(J-JJ)**2  !distance squared to I,J
                I   = II
                J   = JJ
                WRITE(CMOVE,'(a,f6.1)') '(new)',SQRT(REAL(IJD))
              ELSE
                CMOVE = '(sea)'
              ENDIF
            ENDIF
C
            PLONWE = 360.0
            PLONWE = MOD(PLON(I,J),PLONWE)
            IF     (PLONWE.LT.-180.0) THEN  !-360:-180 to  0:180
              PLONWE = PLONWE + 360.0
            ELSEIF (PLONWE.GT. 180.0) THEN  ! 180: 360 to -180:0
              PLONWE = PLONWE - 360.0
            ENDIF
            WRITE(6,'(A,I5,A,2I6,A,F9.3,F8.3,3X,A)') 
     &        'river',KR,
     &        ' at i,j =',I,J,
     &        ' or lon,lat =',PLONWE,PLAT(I,J),TRIM(CMOVE)
C
            IF     (MSK(I,J).NE.1) THEN
              IF     (IGNORE) THEN
                WRITE(6,'(A)') 
     &            'ERROR - RIVER OVER LAND (SKIPPED)'
              ELSE
                WRITE(6,'(A)') 
     &            'ERROR - RIVER OVER LAND (FATAL)'
              ENDIF
              CALL ZHFLSH(6)
              LFATAL = .TRUE.
              CYCLE
            ELSEIF (MIN(I,J).NE.1 .AND.
     &              MINVAL(MSK(MAX(I-1,1):MIN(I+1,IDM),
     &                         MAX(J-1,1):MIN(J+1,JDM) )).NE.0) THEN
              IF     (IGNORE) THEN
                WRITE(6,'(A)') 
     &            'ERROR - RIVER NOT CLOSE TO LAND (SKIPPED)'
              ELSE
                WRITE(6,'(A)') 
     &            'ERROR - RIVER NOT CLOSE TO LAND (FATAL)'
              ENDIF
              CALL ZHFLSH(6)
              LFATAL = .TRUE.
              CYCLE
            ENDIF
            CALL ZHFLSH(6)
          ENDIF
C
          PCM(I,J) = PCM(I,J) + TRIVER(MONTH,KR)*TSCALE*1.E6/
     &                                              DXDY(I,J)
        ENDDO
C
        IF     (LFATAL .AND. .NOT.IGNORE) THEN
          WRITE(6,*) 
          CALL ZHFLSH(6)
          STOP
        ENDIF
C
C       SMOOTH POINT SOURCES OVER NEARBY OCEAN.
C
        CALL SMOOTH1(PCM,MSK,IDM,JDM, 1,NSMOOTH,
     &               LPERIOD,LARCTIC,.FALSE.)
C
C       WRITE OUT HYCOM PRECIP.
C
        CALL ZAIOWR(PCM,MSK,.FALSE., XMIN,XMAX, 14, .FALSE.)
        WRITE(14,4102) '  river precip',MONTH,XMIN,XMAX
        WRITE( 6,4102) '  river precip',MONTH,XMIN,XMAX
        CALL ZHFLSH(6)
      ENDDO !months
C
      CALL ZAIOCL(14)
      CLOSE( UNIT=14)
      STOP
C
 4101 FORMAT(A79)
 4102 FORMAT(A,': month,range =',I3,1X,1P2E16.7)
C     END OF PROGRAM PRIVER.
      END
