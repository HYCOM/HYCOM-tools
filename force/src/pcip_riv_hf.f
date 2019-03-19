      PROGRAM PRIVER
      USE MOD_ZA  ! HYCOM array I/O interface
      IMPLICIT NONE
C
      INTEGER, PARAMETER   ::  MXRIV =99    !only the high frequency rivers
      INTEGER, PARAMETER   ::  MXTIME=366*5 !about 5 years
C
C     PRECIP ARRAY.
C
      INTEGER, ALLOCATABLE ::  MSK(:,:)
      REAL*4,  ALLOCATABLE :: PLON(:,:), PLAT(:,:)
      REAL*4,  ALLOCATABLE :: DXDY(:,:)
      REAL*4,  ALLOCATABLE ::  PCM(:,:),DEPTH(:,:)
      REAL*4,  ALLOCATABLE :: TRIV(:,:)
C
      CHARACTER PREAMBL(5)*79,CLINE*80,CMOVE*11
C
      CHARACTER*79     CTITLE(3)
      NAMELIST/RTITLE/ CTITLE
      INTEGER          NRIVERS,NSMOOTH,MAXDIST,IJRIVER(2,MXRIV)
      REAL*8           WSTART,FSTART,TSTART,TMAX,RSTART,RMAX,RINC
      INTEGER          KRIVER
      REAL*8           TRIVER(MXTIME),TSCALE
      NAMELIST/RIVLOC/ NRIVERS,NSMOOTH,MAXDIST,IJRIVER
      NAMELIST/RIVTIM/ WSTART,FSTART,TSTART,TMAX,RSTART,RMAX,RINC
      NAMELIST/RIVER/  KRIVER,TSCALE,TRIVER
C
C**********
C*
C 1)  GENERATE A HIGH (TIME) FREQUENCY PRECIPITATION BOGAS FIELD
C      PARAMETERIZING RIVER INFLOW.
C     THIS FILE IS USUALLY ADDED TO THE STANDARD PRECIP FORCING FILE.
C     USE time_interp TO BOTH CONVERT RIVERS TO THE PRECIP SAMPLING
C      INTERVAL (IF NECESSARY) AND ADD PRECIP AND RIVERS TOGETHER.
C     NOTE THAT THIS CAN EITHER BE USED IN CONJUCTION WITH pcip_riv_mon,
C     E.G. (A) SOME RIVERS MONTHLY CLIMO AND OTHERS HIGH FREQUENCY OR
C     (B) ALL RIVERS MONTHLY CLIMO WITH (SOME) HIGH FREQUENCY ANOMALIES,
C     OR IT CAN BE USED IN PLACE OF pcip_riv_mon AND forcing.precip.[ab].
C
C 3)  NAMELIST INPUT:
C
C     /RTITLE/
C        CTITLE - THREE (79-CHARACTER) LINE TITLE.
C
C     /RIVLOC/
C        NRIVERS - NUMBER OF RIVERS
C        NSMOOTH - NUMBER OF SMOOTHING PASSES.
C        MAXDIST - MAXIMUM NUMBER OF GRID POINTS TO MOVE RIVERS
C        IJRIVER - I,J LOCATION OF EACH RIVER MOUTH
C
C     /RIVTIM/
C        WSTART - TIME OF WIND START, IGNORED HERE         (DAYS)
C        FSTART - TIME OF HEAT FLUX START                  (DAYS)
C        TSTART - TIME OF START OF CURRENT INTEGRATION     (DAYS)
C        TMAX   - TIME OF END   OF CURRENT INTEGRATION     (DAYS)
C        RSTART - TIME OF RIVER START                      (DAYS)
C        RMAX   - TIME OF RIVER END                        (DAYS)
C        RINC   - TIME INCREMENT BETWEEN INFLOW VALUES     (DAYS)
C
C     /RIVER/  (INVOKED NRIVERS TIMES, WITH KRIVER=1...NRIVER)
C        KRIVER  - RIVER NUMBER
C        TSCALE  - SCALE FACTOR TO CONVERT TRIVER TO SV (1.E6 M^3/S)
C        TRIVER  - RIVER INFLOW (UNITS OF SV/TSCALE)
C                   TRIVER(L) IS FOR TIME RSTART+(L-1)*RINC
C
C 4)  INPUT:
C        ON UNIT  5: NAMELIST /RTITLE/, /RIVLOC/, /RIVTIM/, (multiple) /RIVER/
C        ON UNIT 51: HYCOM BATHYMETRY FILE (.[ab])
C     OUTPUT:
C        ON UNIT 14: HYCOM RIVER BOGAS PRECIP. FILE (.[ab])
C
C 5)  TIME IN DAYS IS W.R.T. JAN 1ST 1901.  BY CONVENTION, FLUX DATES 
C     BEFORE JAN 1ST 1905 INDICATE A CLIMATOLOGY.
C
C 6)  OUTPUT TIMES WILL BE OF THE FORM RSTART+(L-1)*RINC, COVERING
C      FSTART TO FSTART+(TMAX-TSTART).
*
C**********
C
      LOGICAL LARCTIC,LPERIOD,LFATAL
      INTEGER I,II,IJD,J,JJ,K,KR,KREC,NREC
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
      WRITE(6,*) 'READING /RIVLOC/'
      CALL ZHFLSH(6)
      READ( 5,RIVLOC)
      WRITE(6,RIVLOC)
      WRITE(6,*)
      CALL ZHFLSH(6)
C
      IF     (NRIVERS.LE.0) THEN
        WRITE(6,*) 'ERROR - NON-POSITIVE NRIVERS'
        WRITE(6,*)
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
      IF     (NRIVERS.GT.MXRIV) THEN
        WRITE(6,*) 'ERROR - NRIVERS > MXRIV = ',MXRIV
        WRITE(6,*)
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
      FSTART = 0.0
      WSTART = 0.0
      TSTART = 0.0
      TMAX   = 0.0
      RSTART = 0.0
      RMAX   = 0.0
      RINC   = 0.0
      WRITE(6,*) 'READING /RIVTIM/'
      CALL ZHFLSH(6)
      READ( 5,RIVTIM)
      WRITE(6,RIVTIM)
      WRITE(6,*)
      CALL ZHFLSH(6)
C
      NREC = NINT((RMAX-RSTART)/RINC) + 1
C
      IF     (NREC.LE.0 .OR. NREC.GT.MXTIME) THEN
        WRITE(6,*) 'ERROR - NREC OUT OF RANGE'
        WRITE(6,*)
        WRITE(6,*) 'NREC = NINT((RMAX-RSTART)/RINC) + 1'
        WRITE(6,*) 'NREC = ',NREC
        WRITE(6,*)
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
      ALLOCATE( TRIV(NRIVERS,NREC) )
C
      TSCALE = 1.0  ! default is TRIVER in SV, and then last value
      DO KR= 1,NRIVERS
        KRIVER     =    0
        TRIVER(:)  = -1.0E20  !data void
        WRITE(6,*) 'READING /RIVER/'
        CALL ZHFLSH(6)
        READ( 5,RIVER)
        WRITE(6,RIVER)
        CALL ZHFLSH(6)
        IF     (KR.NE.KRIVER) THEN
          WRITE(6,*)
          WRITE(6,*) 'ERROR - KRIVER OUT OF ORDER, SHOULD BE ',KR
          WRITE(6,*)
          CALL ZHFLSH(6)
          STOP
        ENDIF
        IF     (TRIVER(NREC).EQ.-1.0E20) THEN
          WRITE(6,*)
          WRITE(6,*) 'ERROR - TRIVER DOES NOT HAVE ENOUGH VALUES'
          WRITE(6,*) 'NREC,TRIVER(NREC)=',NREC,TRIVER(NREC)
          WRITE(6,*)
          CALL ZHFLSH(6)
          STOP
        ENDIF
        TRIV(KR,1:NREC) = TRIVER(1:NREC)
      ENDDO
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
      CALL ZAIORD(PCM, MSK,.FALSE., HMINA,HMAXA, 21)
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
      CALL ZAIORD(DXDY,MSK,.FALSE., HMINA,HMAXA, 21)
      IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &        ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
        WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b grid files not consistent (pscy):',
     &    '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &    '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
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
      WRITE(PREAMBL(1),'(A,F7.3,A)')
     +   'River inflow as precipitation (m/s), every',RINC,' days'
      PREAMBL(2) = ' '
      PREAMBL(3) = ' '
      PREAMBL(4) = ' '
      PREAMBL(5) = ' '
      I = 1
      DO J= 1,3
        IF     (LEN_TRIM(CTITLE(J)).EQ.0) THEN
          EXIT
        ENDIF
        I = I + 1
        PREAMBL(I) = CTITLE(J)
      ENDDO
      WRITE(PREAMBL(I+1),'(I4,A,2I6)')
     +        NSMOOTH,' smoothing passes.  i/jdm =',IDM,JDM
C
      WRITE(14,4101) PREAMBL
      CALL ZHFLSH(14)
      WRITE(6,*)
      WRITE(6, 4101) PREAMBL
      WRITE(6,*)
      CALL ZHFLSH(6)
C
C     LOOP THROUGH ALL RECORDS.
C
      DO KREC= 1,NREC
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
          IF     (KREC.EQ.1) THEN
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
            WRITE(6,'(A,I3,A,2I6,A,2F9.3,3X,A)') 
     &        'river no.',KR,
     &        ' at i,j =',I,J,
     &        ' or lon,lat =',PLONWE,PLAT(I,J),TRIM(CMOVE)
C
            IF     (MSK(I,J).NE.1) THEN
              WRITE(6,'(A)') 
     &          'ERROR - RIVER OVER LAND'
              CALL ZHFLSH(6)
              LFATAL = .TRUE.
            ELSEIF (MIN(I,J).NE.1 .AND.
     &              MINVAL(MSK(MAX(I-1,1):MIN(I+1,IDM),
     &                         MAX(J-1,1):MIN(J+1,JDM) )).NE.0) THEN
              WRITE(6,'(A)') 
     &          'ERROR - RIVER NOT CLOSE TO LAND'
              CALL ZHFLSH(6)
              LFATAL = .TRUE.
            ENDIF
            CALL ZHFLSH(6)
          ENDIF
C
          PCM(I,J) = PCM(I,J) + TRIV(KR,KREC)*TSCALE*1.E6/
     &                                     (DXDY(I,J))
        ENDDO !kr
C
        IF     (LFATAL) THEN
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
        WRITE(14,4122) 'riverp',RSTART+(KREC-1)*RINC,RINC,XMIN,XMAX
        CALL ZHFLSH(14)
        WRITE( 6,4122) 'riverp',RSTART+(KREC-1)*RINC,RINC,XMIN,XMAX
        CALL ZHFLSH(6)
      ENDDO !krec
C
      CALL ZAIOCL(14)
      CLOSE( UNIT=14)
      STOP
C
 4101 FORMAT(A79)
 4122 FORMAT(2X,A,': day,span,range =',F12.5,F10.6,1P2E16.7)
C     END OF PROGRAM PRIVER.
      END
