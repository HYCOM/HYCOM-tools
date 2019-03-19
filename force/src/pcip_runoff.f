      PROGRAM PRUNOFF
      USE MOD_ZA  ! HYCOM array I/O interface
      IMPLICIT NONE
C
      INTEGER MXRIV
      PARAMETER (MXRIV=9000)
C
C     PRECIP ARRAY.
C
      INTEGER, ALLOCATABLE :: MSK(:,:),ISEA(:,:),JSEA(:,:)
      REAL*4,  ALLOCATABLE :: DXDY(:,:)
      REAL*4,  ALLOCATABLE :: ROI(:,:),PCM(:,:),DEPTH(:,:)
C
      CHARACTER PREAMBL(5)*79,CLINE*240
C
      LOGICAL          IGNORE
      INTEGER          NSMOOTH,MAXDIST
      NAMELIST/RUNOFF/ NSMOOTH,MAXDIST,
     &                 IGNORE
C
C**********
C*
C 1)  GENERATE A MONTHLY PRECIPITATION BOGAS FIELD
C      PARAMETERIZING RUNOFF INFLOW.
C
C 3)  NAMELIST INPUT:
C
C     /RUNOFF/
C        NSMOOTH - NUMBER OF SMOOTHING PASSES.
C        MAXDIST - MAXIMUM NUMBER OF GRID POINTS TO MOVE RUNOFF
C        IGNORE  - IGNORE "BAD" RUNOFF LOCATIONS (DEFAULT .FALSE.)
C
C 4)  INPUT:
C        ON UNIT  5: NAMELIST /RTITLE/, /RUNOFF/
C        ON UNIT 11: HYCOM UNMASKED RUNOFF FILE (.[ab])
C        ON UNIT 51: HYCOM BATHYMETRY      FILE (.[ab])
C     OUTPUT:
C        ON UNIT 14: HYCOM RUNOFF PRECIP.   FILE (.[ab])
C*
C**********
C
      LOGICAL LARCTIC,LPERIOD,LFATAL
      INTEGER I,II,IJD,IJDMAX,IOS,J,JJ,K,MONTH
      REAL*4  X,Y,HMINA,HMINB,HMAXA,HMAXB,XMIN,XMAX
C
C     ARRAYS.
C
      CALL XCSPMD
      ALLOCATE(   MSK(IDM,JDM) )
      ALLOCATE(  ISEA(IDM,JDM) )
      ALLOCATE(  JSEA(IDM,JDM) )
      ALLOCATE(  DXDY(IDM,JDM) )
      ALLOCATE(   ROI(IDM,JDM) )
      ALLOCATE(   PCM(IDM,JDM) )
      ALLOCATE( DEPTH(IDM,JDM) )
C
C     NAMELIST INPUT.
C
      CALL ZHOPEN(6, 'FORMATTED', 'UNKNOWN', 0)
C
      IGNORE  = .FALSE.
      NSMOOTH = 0
      MAXDIST = 0
      WRITE(6,*) 'READING /RUNOFF/'
      CALL ZHFLSH(6)
      READ( 5,RUNOFF)
      WRITE(6,RUNOFF)
      WRITE(6,*)
      CALL ZHFLSH(6)
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
C
      READ(21,*) ! skip plon
      CALL ZAIOSK(21)
      READ(21,*) ! skip plat
      CALL ZAIOSK(21)
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
C     INITIALIZE LAND/SEA MASK AND LOCATION OF NEAREST COASTLINE.
C
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (DEPTH(I,J).GT.  0.0 .AND.
     &            DEPTH(I,J).LE.HMAXA      ) THEN
            MSK(I,J) = 1
          ELSE
            MSK(I,J) = 0
          ENDIF
          ISEA(I,J) = 0
          JSEA(I,J) = 0
        ENDDO
      ENDDO
C
C     CHECK FOR AN ARCTIC BIPOLAR PATCH
C
      LARCTIC = MAXVAL(MSK(:,JDM)).EQ.1  !open top   boundary
      LPERIOD = MAXVAL(MSK(IDM,:)).EQ.1  !open right boundary
C
      WRITE(6,*)
      WRITE(6,*) 'LPERIOD = ',LPERIOD
      WRITE(6,*) 'LARCTIC = ',LARCTIC
      WRITE(6,*)
      CALL ZHFLSH(6)
C
C     INITIALIZE OUTPUT.
C
      CALL ZAIOPN('OLD', 11)
      CALL ZAIOPN('NEW', 14)
C
      CALL ZHOPEN(11, 'FORMATTED', 'OLD', 0)
      CALL ZHOPEN(14, 'FORMATTED', 'NEW', 0)
C
      READ( 11,'(A79)') PREAMBL
      WRITE(PREAMBL(5),'(I4,A)') NSMOOTH,' smoothing passes.'
      WRITE(14,4101) PREAMBL
      WRITE(6,*)
      WRITE(6, 4101) PREAMBL
      WRITE(6,*)
      CALL ZHFLSH(6)
C
C     LOOP THROUGH ALL MONTHS.
C
      IJDMAX = 0
      DO MONTH= 1,HUGE(MONTH)
        READ( 11,          '(A)',IOSTAT=IOS) CLINE
        IF     (IOS.NE.0) THEN
          EXIT
        ENDIF
        CALL ZAIORD(ROI,MSK,.FALSE., XMIN,XMAX, 11)
C
        DO J= 1,JDM
          DO I= 1,IDM
            PCM(I,J) = 0.0
          ENDDO
        ENDDO
C
        LFATAL = .FALSE.
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (ROI(I,J).EQ.0.0) THEN
              CYCLE  !nothing to do
            ENDIF
            IF     (MSK(I,J).NE.1) THEN
              IF     (ISEA(I,J).EQ.0) THEN
C
C               FIND A CLOSE SEA POINT, IF THERE IS ONE.
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
                          ISEA(I,J) = II
                          JSEA(I,J) = JJ
                        ENDIF
                      ENDIF
                    JJ=MIN(J+K,JDM)
                      IF     (MSK(II,JJ).EQ.1) THEN
                        IF     ((I-II)**2+(J-JJ)**2.LT.IJD) THEN
                          IJD = (I-II)**2+(J-JJ)**2  !distance squared to I,J
                          ISEA(I,J) = II
                          JSEA(I,J) = JJ
                        ENDIF
                      ENDIF
                  ENDDO !ii
                  DO JJ= MAX(J-K+1,1),MIN(J+K-1,JDM)
                    II=MAX(I-K,1)
                      IF     (MSK(II,JJ).EQ.1) THEN
                        IF     ((I-II)**2+(J-JJ)**2.LT.IJD) THEN
                          IJD = (I-II)**2+(J-JJ)**2  !distance squared to I,J
                          ISEA(I,J) = II
                          JSEA(I,J) = JJ
                        ENDIF
                      ENDIF
                    II=MIN(I+K,IDM)
                      IF     (MSK(II,JJ).EQ.1) THEN
                        IF     ((I-II)**2+(J-JJ)**2.LT.IJD) THEN
                          IJD = (I-II)**2+(J-JJ)**2  !distance squared to I,J
                          ISEA(I,J) = II
                          JSEA(I,J) = JJ
                        ENDIF
                      ENDIF
                  ENDDO !jj
                ENDDO !k
              ENDIF !ISEA==0
              IF     (ISEA(I,J).NE.0) THEN
                II  = ISEA(I,J)
                JJ  = JSEA(I,J)
                IJD = (I-II)**2+(J-JJ)**2  !distance squared to I,J
                IF     (IJD.GT.IJDMAX) THEN
                  IJDMAX = IJD
                  WRITE(6,'(a,2i6,f12.2)')
     &              'I,J,DIST =',I,J,SQRT(REAL(IJD))
                ENDIF
                PCM(II,JJ) = PCM(II,JJ) + ROI(I,J)*DXDY(I,J)
                ROI(I, J ) = 0.0
              ELSE
                WRITE(6,*) "can't move land point",i,j," to the coast"
                LFATAL = .TRUE.
                II = ISEA(I,J)
                JJ = JSEA(I,J)
              ENDIF !ISEA/=0
            ELSE !MSK(I,J).EQ.1
              IF     (MINVAL(MSK(MAX(I-1,1):MIN(I+1,IDM),
     &                           MAX(J-1,1):MIN(J+1,JDM) )).NE.0) THEN
                IF     (ISEA(I,J).EQ.0) THEN
C
C                 FIND A SEA POINT  CLOSER TO LAND, IF THERE IS ONE.
C
                  IJD = MAXDIST**2 + 1
                  DO K= 1,MAXDIST
                    IF     (IJD.LT.K**2) THEN
                      EXIT
                    ENDIF  !early exit
                    DO II= MAX(I-K,1),MIN(I+K,IDM)
                      JJ=MAX(J-K,1)
                        IF     (MSK(II,JJ).EQ.1 .AND.
     &                          MINVAL(MSK(MAX(II-1,1):MIN(II+1,IDM),
     &                                     MAX(JJ-1,1):MIN(JJ+1,JDM) ))
     &                          .EQ.0) THEN
                          IF     ((I-II)**2+(J-JJ)**2.LT.IJD) THEN
                            IJD = (I-II)**2+(J-JJ)**2  !distance squared to I,J
                            ISEA(I,J) = II
                            JSEA(I,J) = JJ
                          ENDIF
                        ENDIF
                      JJ=MIN(J+K,JDM)
                        IF     (MSK(II,JJ).EQ.1 .AND.
     &                          MINVAL(MSK(MAX(II-1,1):MIN(II+1,IDM),
     &                                     MAX(JJ-1,1):MIN(JJ+1,JDM) ))
     &                          .EQ.0) THEN
                          IF     ((I-II)**2+(J-JJ)**2.LT.IJD) THEN
                            IJD = (I-II)**2+(J-JJ)**2  !distance squared to I,J
                            ISEA(I,J) = II
                            JSEA(I,J) = JJ
                          ENDIF
                        ENDIF
                    ENDDO !ii
                    DO JJ= MAX(J-K+1,1),MIN(J+K-1,JDM)
                      II=MAX(I-K,1)
                        IF     (MSK(II,JJ).EQ.1 .AND.
     &                          MINVAL(MSK(MAX(II-1,1):MIN(II+1,IDM),
     &                                     MAX(JJ-1,1):MIN(JJ+1,JDM) ))
     &                          .EQ.0) THEN
                          IF     ((I-II)**2+(J-JJ)**2.LT.IJD) THEN
                            IJD = (I-II)**2+(J-JJ)**2  !distance squared to I,J
                            ISEA(I,J) = II
                            JSEA(I,J) = JJ
                          ENDIF
                        ENDIF
                      II=MIN(I+K,IDM)
                        IF     (MSK(II,JJ).EQ.1 .AND.
     &                          MINVAL(MSK(MAX(II-1,1):MIN(II+1,IDM),
     &                                     MAX(JJ-1,1):MIN(JJ+1,JDM) ))
     &                          .EQ.0) THEN
                          IF     ((I-II)**2+(J-JJ)**2.LT.IJD) THEN
                            IJD = (I-II)**2+(J-JJ)**2  !distance squared to I,J
                            ISEA(I,J) = II
                            JSEA(I,J) = JJ
                          ENDIF
                        ENDIF
                    ENDDO !jj
                  ENDDO !k
                ENDIF !isea==0
                IF     (ISEA(I,J).EQ.0) THEN  !no better sea point
                  ISEA(I,J) = I
                  JSEA(I,J) = J
                ENDIF !isea==0
              ELSE  !sea next to land
                ISEA(I,J) = I
                JSEA(I,J) = J
              ENDIF !sea away from land:else
              II  = ISEA(I,J)
              JJ  = JSEA(I,J)
              IJD = (I-II)**2+(J-JJ)**2  !distance squared to I,J
              IF     (IJD.GT.IJDMAX) THEN
                IJDMAX = IJD
                WRITE(6,'(a,2i6,f12.2)')
     &            'I,J,DIST =',I,J,SQRT(REAL(IJD))
              ENDIF
              PCM(II,JJ) = PCM(II,JJ) + ROI(I,J)*DXDY(I,J)
              ROI(I, J ) = 0.0
            ENDIF  !msk.ne.1:else
          ENDDO !j
        ENDDO !i
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
        DO J= 1,JDM
          DO I= 1,IDM
            PCM(I,J) = PCM(I,J)/DXDY(I,J)
          ENDDO
        ENDDO
C
C       WRITE OUT HYCOM PRECIP.
C
        CALL ZAIOWR(PCM,MSK,.FALSE., XMIN,XMAX, 14, .FALSE.)
        I = INDEX(CLINE,'=') + 4
        WRITE(14,'(A,1P2E16.7)') CLINE(1:I),XMIN,XMAX
        WRITE( 6,'(A,1P2E16.7)') CLINE(1:I),XMIN,XMAX
        CALL ZHFLSH(6)
      ENDDO !months
C
      CALL ZAIOCL(14)
      CLOSE( UNIT=14)
      STOP
C
 4101 FORMAT(A79)
C     END OF PROGRAM PRUNOFF.
      END
