      PROGRAM WI_MAGSTRESS
      USE MOD_ZA  ! HYCOM array I/O interface
      IMPLICIT NONE
C
C     WIND ARRAYS.
C
      INTEGER, ALLOCATABLE :: MSK(:,:)
      REAL*4,  ALLOCATABLE :: TXM(:,:),TYM(:,:),TMAG(:,:),CORI(:,:),
     &                        TXU(:),TYV(:)
C
      CHARACTER PREAMBL(5)*79
C
C     NAMELIST.
C
      INTEGER          ITYPE
      REAL*4           CEKMAN,CORMIN
      NAMELIST/MAGSTR/ ITYPE,CEKMAN,CORMIN
C
C**********
C*
C 1)  CALCULATE THE MAGNITUDE OF EXISTING HYCOM WIND STRESS FILES.
C     OF THE EKMAN DEPTH.
C
C 2)  NAMELIST INPUT:
C
C     /MAGSTR/
C        ITYPE   - TYPE OF OUTPUT
C                   =0; WIND STRESS MAGNITUDE (DEFAULT)
C                   =1; EKMAN DEPTH
C        CEKMAN  - SCALE FACTOR FOR EKMAN DEPTH (DEFAULT 0.7)
C        CORMIN  - MINIMUM CORIOLIS MAGNITUDE   (DEFAULT 1.0E-5, OR 4N)
C
C 3)  INPUT
C        ON UNIT 10:    .a/.b FORMAT MODEL TAUEWD  FILE, SEE (3).
C        ON UNIT 11:    .a/.b FORMAT MODEL TAUNWD  FILE, SEE (3).
C        ON UNIT 21:    .a/.b FORMAT regional.grid FILE.
C     OUTPUT:
C        ON UNIT 12:    .a/.b FORMAT MODEL STRMAG  FILE, SEE (4).
C
C 4)  THE INPUT WIND STRESSES HAVE THEIR COMPONENTS ON EITHER
C      EVERY POINT OF THE MODEL'S 'U' AND 'V' GRIDS RESPECTIVELY,
C      OR BOTH ON EVERY POINT OF THE MODEL'S 'P' GRID.  IF ON THE
C      U AND V GRIDS, IT IS FIRST INTERPOLATED TO THE P-GRID AND
C      THEN THE MAGNITUDE CALCULATED FROM THE STRESSES ON THIS GRID.
C
C 5)  THE STRESS MAGNITUDE IS CALCULATED AND OUTPUT ON THE P-GRID.
C      ARRAY SIZE IS 'IDM' BY 'JDM', AND THE DATA IS OUTPUT .a/.b
C      FORMAT TOGETHER WITH EITHER (A) THE MONTH, OR (B) THE DAY THE
C      WIND REPRESENTS AND THE INCREMENT IN DAYS TO THE NEXT WIND
C      RECORD.  
C
C 6)  ALAN J. WALLCRAFT,  NRL,  NOVEMBER 2004 AND JUNE 2006.
C*
C**********
C
      LOGICAL      PGRID,GLOBAL
      CHARACTER*80 CLINE
      INTEGER      I,IOS,IM1,IP1,J,JM1,JP1,KREC
      REAL*4       HMINA,HMINB,HMAXA,HMAXB
      REAL*4       XMIN,XMAX
C
C --- MODEL ARRAYS.
C
      CALL XCSPMD  !define idm,jdm
      ALLOCATE(   MSK(IDM,JDM) )
      ALLOCATE(   TXM(IDM,JDM) )
      ALLOCATE(   TYM(IDM,JDM) )
      ALLOCATE(  TMAG(IDM,JDM) )
      ALLOCATE(  CORI(IDM,JDM) )
      ALLOCATE(   TXU(1:IDM+1) )
      ALLOCATE(   TYV(1:JDM+1) )
C
C     GRID INPUT.
C
      CALL ZAIOST
C
C     NAMELIST INPUT.
C
      CALL ZHOPEN(6, 'FORMATTED', 'UNKNOWN', 0)
C
      ITYPE   = 0
      CEKMAN  = 0.7
      CORMIN  = 1.0E-5 !4N
      WRITE(6,*) 'READING /MAGSTR/'
      CALL ZHFLSH(6)
      READ( 5,MAGSTR)
      WRITE(6,MAGSTR)
      WRITE(6,*) 
      CALL ZHFLSH(6)
C
      CALL ZHOPNC(21, 'regional.grid.b', 'FORMATTED', 'OLD', 0)
      CALL ZAIOPF('regional.grid.a', 'OLD', 21)
C
      READ(21,*) ! skip idm
      READ(21,*) ! skip jdm
      READ(21,*) ! skip mapflg
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
      READ(21,*) ! skip pscx
      CALL ZAIOSK(21)
      READ(21,*) ! skip pscy
      CALL ZAIOSK(21)
      READ(21,*) ! skip qscx
      CALL ZAIOSK(21)
      READ(21,*) ! skip qscy
      CALL ZAIOSK(21)
      READ(21,*) ! skip uscx
      CALL ZAIOSK(21)
      READ(21,*) ! skip uscy
      CALL ZAIOSK(21)
      READ(21,*) ! skip vscx
      CALL ZAIOSK(21)
      READ(21,*) ! skip vscy
      CALL ZAIOSK(21)
C
      READ(21,'(A)') CLINE
      WRITE(6,'(A)') CLINE
      I = INDEX(CLINE,'=')
      READ(CLINE(I+1:),*) HMINB,HMAXB
      CALL ZAIORD(CORI,MSK,.FALSE., HMINA,HMAXA, 21)
      IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &        ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
        WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b grid files not consistent (cori):',
     &    '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &    '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
      CLOSE(UNIT=21)
      CALL ZAIOCL(21)
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
      PGRID = INDEX(PREAMBL(2),'p-grid').NE.0
      IF     (ITYPE.EQ.0) THEN
        PREAMBL(2) = 'Magnitude of Wind Stress on p-grid'
      ELSE
        PREAMBL(2) = 'Ekman Depth on p-grid'
      ENDIF
      WRITE(12,'(A79)') PREAMBL
      WRITE( 6,'(A79)') PREAMBL
C
      IF     (PGRID) THEN
        WRITE(6,*) 
        WRITE(6,*) 'Wind stresses input on p-grid.'
        WRITE(6,*) 
      ELSE
        WRITE(6,*) 
        WRITE(6,*) 'Wind stresses first interpolated to p-grid.'
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
        ELSEIF (CLINE.EQ.' ') THEN
          EXIT !probably the effective end of file
        ENDIF
        READ(11,'(A)') CLINE
        CALL ZHFLSH( 6)
        CALL ZAIORD(TXM,MSK,.FALSE., XMIN,XMAX, 10)
        CALL ZAIORD(TYM,MSK,.FALSE., XMIN,XMAX, 11)
        IF     (.NOT.PGRID) THEN
          DO J= 1,JDM
            DO I= 1,IDM
              TXU(I) = TXM(I,J)
            ENDDO
            IF     (GLOBAL) THEN
              TXU(IDM+1) = TXU(  1)  ! assume     periodic boundary
            ELSE
              TXU(IDM+1) = TXU(IDM)  ! assume non-periodic boundary
            ENDIF
            DO I= 1,IDM
              TXM(I,J) = 0.5*(TXU(I) + TXU(I+1))
            ENDDO
          ENDDO
          DO I= 1,IDM
            DO J= 1,JDM
              TYV(J) = TYM(I,J)
            ENDDO
            TYV(JDM+1) = TYV(JDM)  ! assume no change across southern boundary
            DO J= 1,JDM
              TYM(I,J) = 0.5*(TYV(J) + TYV(J+1))
            ENDDO
          ENDDO
        ENDIF
C
C       CALCULATE MAGNITUDE ON P-GRID.
C
        DO J= 1,JDM
          DO I= 1,IDM
            TMAG(I,J) = SQRT(TXM(I,J)**2 + TYM(I,J)**2)
            IF     (ITYPE.EQ.1) THEN
              TMAG(I,J) = SQRT(1.e-3*TMAG(I,J))  !u-star
              TMAG(I,J) = TMAG(I,J)*(CEKMAN*4.0)/
     &                    MAX( CORMIN*4.0,
     &                         ABS(CORI(I,J  ))+ABS(CORI(I+1,J  ))+
     &                         ABS(CORI(I,J+1))+ABS(CORI(I+1,J+1)) )
            ENDIF
          ENDDO
        ENDDO
C
        IF     (ITYPE.NE.1) THEN
C
C         WRITE OUT HYCOM STRESS MAGNITUDE.
C
          CALL ZAIOWR(TMAG,MSK,.FALSE., XMIN,XMAX, 12, .FALSE.)
          IF     (INDEX(CLINE,'month').NE.0) THEN
            WRITE(12,'(2A,1P2E16.7)') ' str_mag',CLINE(9:26),XMIN,XMAX
            WRITE( 6,'(2A,1P2E16.7)') ' str_mag',CLINE(9:26),XMIN,XMAX
          ELSEIF (INDEX(CLINE,'=').EQ.10) THEN
            WRITE(12,'(A,1P2E16.7)') CLINE(1:11),XMIN,XMAX
            WRITE( 6,'(A,1P2E16.7)') CLINE(1:11),XMIN,XMAX
          ELSE
            WRITE(12,'(2A,1P2E16.7)') ' str_mag',CLINE(9:45),XMIN,XMAX
            WRITE( 6,'(2A,1P2E16.7)') ' str_mag',CLINE(9:45),XMIN,XMAX
          ENDIF
          CALL ZHFLSH(12)
        ELSE
C
C         WRITE OUT HYCOM STRESS MAGNITUDE.
C
          CALL ZAIOWR(TMAG,MSK,.FALSE., XMIN,XMAX, 12, .FALSE.)
          IF     (INDEX(CLINE,'month').NE.0) THEN
            WRITE(12,'(2A,1P2E16.7)') '  hekman',CLINE(9:26),XMIN,XMAX
            WRITE( 6,'(2A,1P2E16.7)') '  hekman',CLINE(9:26),XMIN,XMAX
          ELSEIF (INDEX(CLINE,'=').EQ.10) THEN
            WRITE(12,'(A,1P2E16.7)') CLINE(1:11),XMIN,XMAX
            WRITE( 6,'(A,1P2E16.7)') CLINE(1:11),XMIN,XMAX
          ELSE
            WRITE(12,'(2A,1P2E16.7)') '  hekman',CLINE(9:45),XMIN,XMAX
            WRITE( 6,'(2A,1P2E16.7)') '  hekman',CLINE(9:45),XMIN,XMAX
          ENDIF
        ENDIF !itype
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
C     END OF PROGRAM WI_MAGSTRESS.
      END
