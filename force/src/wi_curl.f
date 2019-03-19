      PROGRAM W_CURL
      USE MOD_ZA  ! HYCOM array I/O interface
      IMPLICIT NONE
C
C     WIND ARRAYS.
C
      INTEGER, ALLOCATABLE :: MSK(:,:)
      REAL*4,  ALLOCATABLE :: RSCX(:,:),RSCY(:,:),PANG(:,:)
      REAL*4,  ALLOCATABLE :: TXM(:,:),TYM(:,:),CURL(:,:),QCURL(:,:),
     &                        TXP(:),TYP(:)
C
      REAL*4,  PARAMETER   :: ONEM=1.0
C
      CHARACTER PREAMBL(5)*79
C
C**********
C*
C 1)  CALCULATE THE CURL OF EXISTING HYCOM WIND STRESS FILES.
C     ALSO WORKS ON 10M WIND VELOCITY FILES.
C
C 2)  INPUT
C        ON UNIT 10:    .a/.b FORMAT MODEL TAUEWD  FILE, SEE (3).
C        ON UNIT 11:    .a/.b FORMAT MODEL TAUNWD  FILE, SEE (3).
C        ON UNIT 21:    .a/.b FORMAT regional.grid FILE.
C     OUTPUT:
C        ON UNIT 12:    .a/.b FORMAT MODEL CURL    FILE, SEE (4).
C
C 3)  THE INPUT WIND STRESSES HAVE THEIR COMPONENTS ON EITHER
C      EVERY POINT OF THE MODEL'S 'U' AND 'V' GRIDS RESPECTIVELY,
C      OR BOTH ON EVERY POINT OF THE MODEL'S 'P' GRID.  IF ON THE
C      P-GRID, IT IS FIRST INTERPOLATED TO THE U AND V GRIDS AND
C      THEN THE CURL CALCULATED FROM THE STRESSES ON THESE GRIDS.
C
C 4)  THE CURL IS CALCULATED ON THE STREAMFUNCTION-GRID, AS IS
C      NATURAL FOR STRESSES ON THE U AND V GRIDS, AND INTERPOLATED
C      TO THE P-GRID FOR OUTPUT.  THIS MAY REDUCE ACCURACY.
C      ARRAY SIZE IS 'IDM' BY 'JDM', AND THE DATA IS OUTPUT .a/.b
C      FORMAT TOGETHER WITH EITHER (A) THE MONTH, OR (B) THE DAY THE
C      WIND REPRESENTS AND THE INCREMENT IN DAYS TO THE NEXT WIND
C      RECORD.  
C
C 5)  ALAN J. WALLCRAFT,  NRL,  NOVEMBER 2002.
C*
C**********
C
      LOGICAL      PGRID,GLOBAL,ARCTIC,ROTATE
      CHARACTER*80 CLINE
      INTEGER      I,II,IOS,IM1,IP1,J,JJ,JM1,JP1,KREC
      REAL*4       TXMIJ,TYMIJ,COSPANG,SINPANG
      REAL*4       HMINA,HMINB,HMAXA,HMAXB
      REAL*4       XMIN,XMAX
C
C --- MODEL ARRAYS.
C
      CALL XCSPMD  !define idm,jdm
      ALLOCATE(   MSK(IDM,JDM) )
      ALLOCATE(  RSCX(IDM,JDM) )
      ALLOCATE(  RSCY(IDM,JDM) )
      ALLOCATE(  PANG(IDM,JDM) )
      ALLOCATE(   TXM(IDM,JDM) )
      ALLOCATE(   TYM(IDM,JDM) )
      ALLOCATE(  CURL(IDM,JDM) )
      ALLOCATE( QCURL(IDM,JDM) )
      ALLOCATE(   TXP(0:IDM) )
      ALLOCATE(   TYP(0:JDM) )
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
      READ(21,'(A)') CLINE  !read plon range
      WRITE(6,'(A)') CLINE
      I = INDEX(CLINE,'=')
      READ(CLINE(I+1:),*) HMINB,HMAXB
      GLOBAL = HMAXB-HMINB .GT. 350.0
      CALL ZAIORD(RSCX,MSK,.FALSE., HMINA,HMAXA, 21)  !temporarily store plon
      READ(21,*) ! skip plat
      CALL ZAIORD(RSCY,MSK,.FALSE., HMINA,HMAXA, 21)  !temporarily store plat
      ARCTIC = RSCX(3*IDM/4+1,JDM-1).EQ.RSCX(IDM/4,JDM) .AND.
     &         RSCY(3*IDM/4+1,JDM-1).EQ.RSCY(IDM/4,JDM)
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
C
      READ(21,'(A)') CLINE
      WRITE(6,'(A)') CLINE
      I = INDEX(CLINE,'=')
      READ(CLINE(I+1:),*) HMINB,HMAXB
      CALL ZAIORD(PANG,MSK,.FALSE., HMINA,HMAXA, 21)
      IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &        ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
        WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b grid files not consistent (pang):',
     &    '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &    '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
      READ(21,*) ! skip pscx
      CALL ZAIOSK(21)
      READ(21,*) ! skip pscy
      CALL ZAIOSK(21)
C
      READ(21,'(A)') CLINE
      WRITE(6,'(A)') CLINE
      I = INDEX(CLINE,'=')
      READ(CLINE(I+1:),*) HMINB,HMAXB
      CALL ZAIORD(RSCX,MSK,.FALSE., HMINA,HMAXA, 21)
      IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &        ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
        WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b grid files not consistent (qscx):',
     &    '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &    '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
      READ(21,'(A)') CLINE
      WRITE(6,'(A)') CLINE
      I = INDEX(CLINE,'=')
      READ(CLINE(I+1:),*) HMINB,HMAXB
      CALL ZAIORD(RSCY,MSK,.FALSE., HMINA,HMAXA, 21)
      IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &        ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
        WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b grid files not consistent (qscy):',
     &    '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &    '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
      CLOSE(UNIT=21)
      CALL ZAIOCL(21)
C
      DO I= 1,IDM
        DO J= 1,JDM
          RSCX(I,J) = 1.0/MAX(ONEM,RSCX(I,J)) !1/qscx
          RSCY(I,J) = 1.0/MAX(ONEM,RSCY(I,J)) !1/qscy
        ENDDO
      ENDDO
C
C     INITIALIZE HYCOM INPUT AND OUTPUT.
C
      CALL ZAIOPN('OLD', 10)
      CALL ZAIOPN('OLD', 11)
      CALL ZAIOPN('NEW', 12)
      CALL ZHOPEN(10, 'FORMATTED', 'OLD', 0)
      CALL ZHOPEN(11, 'FORMATTED', 'OLD', 0)
      CALL ZHOPEN(12, 'FORMATTED', 'NEW', 0)
C
      READ( 10,'(A79)') PREAMBL
      READ( 11,'(A79)') PREAMBL
      PGRID  = INDEX(PREAMBL(2),'p-grid') .NE.0
      ROTATE = INDEX(PREAMBL(2),'E-wards').NE.0
      IF     (INDEX(PREAMBL(2),'Velocities').EQ.0) THEN
        PREAMBL(2) = 'Wind Stress Curl on p-grid'
      ELSE
        PREAMBL(2) = 'Wind Velocity Curl on p-grid'
      ENDIF
      WRITE(12,'(A79)') PREAMBL
      WRITE( 6,'(A79)') PREAMBL
C
      IF     (PGRID) THEN
        WRITE(6,*) 
        IF     (ROTATE) THEN
          WRITE(6,*) 
     &      'Winds first rotated and interpolated to u and v grids.'
        ELSE
          WRITE(6,*) 
     &      'Winds first interpolated to u and v grids.'
        ENDIF
        WRITE(6,*) 
      ELSE
        IF     (ARCTIC) THEN
          WRITE(6,'(/ a /)')
     &    'error - wind stresses must be on p-grid for tri-pole regions'
          CALL ZHFLSH(6)
          STOP
        ELSEIF (ROTATE) THEN
          WRITE(6,'(/ a /)')
     &    'error - E-ward,N-ward wind stresses must be on p-grid'
          CALL ZHFLSH(6)
          STOP
        ENDIF
        WRITE(6,*) 
        WRITE(6,*) 'Wind stresses input on u and v grids.'
        WRITE(6,*) 
      ENDIF
C
C     PROCESS ALL THE WIND RECORDS.
C
      DO 810 KREC= 1,999999
C
C       READ THE INPUT WIND STRESSES.
C
        READ(10,'(A)',IOSTAT=IOS) CLINE
        IF     (IOS.NE.0) THEN
          EXIT !probably end of file
        ENDIF
        WRITE(6,'(A)') CLINE
        READ(11,'(A)') CLINE
        WRITE(6,'(A)') CLINE
        CALL ZAIORD(TXM,MSK,.FALSE., XMIN,XMAX, 10)
        CALL ZAIORD(TYM,MSK,.FALSE., XMIN,XMAX, 11)
        IF     (ROTATE) THEN
C         rotate from e-ward,n-ward to x-ward,y-ward
          DO J= 1,JDM
            DO I= 1,IDM
              TXMIJ    = TXM(I,J)
              TYMIJ    = TYM(I,J)
              COSPANG  = COS(PANG(I,J))
              SINPANG  = SIN(PANG(I,J))
              TXM(I,J) = COSPANG*TXMIJ + SINPANG*TYMIJ
              TYM(I,J) = COSPANG*TYMIJ - SINPANG*TXMIJ
            ENDDO
          ENDDO
        ENDIF !rotate
        IF     (PGRID) THEN
C         interpolate from p ro u and v grids
          DO J= 1,JDM
            DO I= 1,IDM
              TXP(I) = TXM(I,J)
            ENDDO
            IF     (GLOBAL) THEN
              TXP(0) = TXP(IDM)  ! assume     periodic boundary
            ELSE
              TXP(0) = TXP(  1)  ! assume non-periodic boundary
            ENDIF
            DO I= 1,IDM
              TXM(I,J) = 0.5*(TXP(I) + TXP(I-1))
            ENDDO
          ENDDO
          DO I= 1,IDM
            DO J= 1,JDM
              TYP(J) = TYM(I,J)
            ENDDO
            TYP(0) = TYP(1)  ! assume no change across southern boundary
            DO J= 1,JDM
              TYM(I,J) = 0.5*(TYP(J) + TYP(J-1))
            ENDDO
          ENDDO
        ENDIF !pgrid
C
C       CALCULATE CURL ON Q-GRID.
C
        DO J= 1,JDM
          IF     (J.NE.1) THEN
            JM1 = J-1
          ELSE
            JM1 = 1    ! assume no change across southern boundary
          ENDIF
          DO I= 1,IDM
            IF     (I.NE.1) THEN
              IM1 = I-1
            ELSEIF (GLOBAL) THEN
              IM1 = IDM  ! assume     periodic boundary
            ELSE
              IM1 =   1  ! assume non-periodic boundary
            ENDIF
            QCURL(I,J) = RSCX(I,J)*(TYM(I,J) - TYM(IM1,J)) -
     &                   RSCY(I,J)*(TXM(I,J) - TXM(I,JM1))
          ENDDO
        ENDDO
        IF     (ARCTIC) THEN !correct northern boundary, q-grid
          J  = JDM
          JJ = JDM-1
          DO I= 1,IDM
            II = MOD(IDM-(I-1),IDM)+1
            QCURL(I,J) = QCURL(II,JJ)
          ENDDO !i
        ENDIF
C
C       INTERPOLATE CURL TO P-GRID.
C
        DO J= 1,JDM
          IF     (J.NE.JDM) THEN
            JP1 = J+1
          ELSE
            JP1 = JDM  ! assume no change across northern boundary
          ENDIF
          DO I= 1,IDM
            IF     (I.NE.IDM) THEN
              IP1 = I+1
            ELSEIF (GLOBAL) THEN
              IP1 =   1  ! assume     periodic boundary
            ELSE
              IP1 = IDM  ! assume non-periodic boundary
            ENDIF
            CURL(I,J) = 0.25*(QCURL(I,  J)   +
     &                        QCURL(IP1,J)   +
     &                        QCURL(I,  JP1) +
     &                        QCURL(IP1,JP1)  )
          ENDDO
        ENDDO
        IF     (ARCTIC) THEN !correct northern boundary, p-grid
          J  = JDM
          JJ = JDM-1
          DO I= 1,IDM
            II = IDM-MOD(I-1,IDM)
            CURL(I,J) = CURL(II,JJ)
          ENDDO !i
        ENDIF
C
C       WRITE OUT HYCOM CURL.
C
        CALL ZAIOWR(CURL,MSK,.FALSE., XMIN,XMAX, 12, .FALSE.)
        IF     (INDEX(CLINE,'month').NE.0) THEN
          WRITE(12,'(2A,1P2E16.7)') ' wndcurl',CLINE(9:26),XMIN,XMAX
          WRITE( 6,'(2A,1P2E16.7)') ' wndcurl',CLINE(9:26),XMIN,XMAX
        ELSE
          WRITE(12,'(2A,1P2E16.7)') ' wndcurl',CLINE(9:45),XMIN,XMAX
          WRITE( 6,'(2A,1P2E16.7)') ' wndcurl',CLINE(9:45),XMIN,XMAX
        ENDIF
        CALL ZHFLSH(12)
        CALL ZHFLSH( 6)
  810 CONTINUE
C
      CALL ZAIOCL(10)
      CLOSE( UNIT=10)
      CALL ZAIOCL(11)
      CLOSE( UNIT=11)
      CALL ZAIOCL(12)
      CLOSE( UNIT=12)
      STOP
C
C     END OF PROGRAM WNDCUR.
      END
