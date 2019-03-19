      PROGRAM WI_XWDYWD
      USE MOD_ZA  ! HYCOM array I/O interface
      IMPLICIT NONE
C
C     WIND ARRAYS.
C
      INTEGER, ALLOCATABLE :: MSK(:,:)
      REAL*4,  ALLOCATABLE :: PANG(:,:)
      REAL*4,  ALLOCATABLE :: TXM(:,:),TYM(:,:)
C
      CHARACTER PREAMBL(5)*79
C
C**********
C*
C 1)  ROTATE HYCOM WIND STRESS FROM E-WARD,N-WARD TO X-WARD,Y-WARD.
C     NOTE THAT USUALLY HYCOM WIND STRESS IS X-WARD,Y-WARD, SO THIS
C     PROGRAM IS PRIMARILY FOR TESTING.
C
C 2)  INPUT
C        ON UNIT 10:    .a/.b FORMAT MODEL TAUXWD  FILE, SEE (3).
C        ON UNIT 11:    .a/.b FORMAT MODEL TAUYWD  FILE, SEE (3).
C        ON UNIT 21:    .a/.b FORMAT regional.grid FILE.
C     OUTPUT:
C        ON UNIT 12:    .a/.b FORMAT MODEL TAUEWD   FILE, SEE (7).
C        ON UNIT 13:    .a/.b FORMAT MODEL TAUNWD   FILE, SEE (7).
C
C 3)  THE INPUT (AND OUTPUT) WIND STRESSES HAVE THEIR COMPONENTS 
C      BOTH ON EVERY POINT OF THE MODEL'S 'P' GRID.  THE FIELDS
C      MAY INCLUDE 2.0**100 TO INDICATE A DATA VOID.
C
C 4)  ALAN J. WALLCRAFT,  NRL,  FEBRUARY, 2015.
C*
C**********
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*80 CLINE,CLINE_X,CLINE_Y
      LOGICAL      LARCTIC
      INTEGER      I,IOS,J,KREC
      REAL*4       TXIJ,TYIJ
      REAL*4       HMINA,HMINB,HMAXA,HMAXB
      REAL*4       XMIN,XMAX
C
C --- MODEL ARRAYS.
C
      CALL XCSPMD  !define idm,jdm
      ALLOCATE(   MSK(IDM,JDM) )
      ALLOCATE(  PANG(IDM,JDM) )
      ALLOCATE(   TXM(IDM,JDM) )
      ALLOCATE(   TYM(IDM,JDM) )
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
C
      I = INDEX(CLINE,'=')
      READ(CLINE(I+1:),*) HMINB,HMAXB
      CALL ZAIORD(TXM,MSK,.FALSE., HMINA,HMAXA, 21)  !plon
      IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &        ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
        WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b grid files not consistent (plon):',
     &    '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &    '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
      READ(21,'(A)') CLINE
      I = INDEX(CLINE,'=')
      READ(CLINE(I+1:),*) HMINB,HMAXB
      CALL ZAIORD(TYM,MSK,.FALSE., HMINA,HMAXA, 21) !plat
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
C
      READ(21,'(A)') CLINE
      WRITE(6,'(A)') CLINE
      I = INDEX(CLINE,'=')
      READ(CLINE(I+1:),*) HMINB,HMAXB
      CALL ZAIORD(PANG,MSK,.FALSE., HMINA,HMAXA, 21)  !pang
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
      CLOSE(UNIT=21)
      CALL ZAIOCL(21)
C
C     CHECK FOR AN ARCTIC BIPOLAR PATCH
C
      LARCTIC =  TXM(3*IDM/4+1,JDM-1).EQ. TXM(IDM/4,JDM) .AND.
     &           TYM(3*IDM/4+1,JDM-1).EQ. TYM(IDM/4,JDM)
C
      WRITE(6,*)
      WRITE(6,*) 'LARCTIC = ',LARCTIC
      WRITE(6,*)
      CALL ZHFLSH(6)
C
C     INITIALIZE HYCOM INPUT AND OUTPUT.
C
      CALL ZAIOPN('OLD', 10)
      CALL ZAIOPN('OLD', 11)
      CALL ZAIOPN('NEW', 12)
      CALL ZAIOPN('NEW', 13)
      CALL ZHOPEN(10, 'FORMATTED', 'OLD', 0)
      CALL ZHOPEN(11, 'FORMATTED', 'OLD', 0)
      CALL ZHOPEN(12, 'FORMATTED', 'NEW', 0)
      CALL ZHOPEN(13, 'FORMATTED', 'NEW', 0)
C
      READ( 10,'(A79)') PREAMBL
      PREAMBL(3) = 'rotated from e-ward to x-ward by wi_xwdywd'
      WRITE(12,'(A79)') PREAMBL
      WRITE( 6,'(A79)') PREAMBL
C
      READ( 11,'(A79)') PREAMBL
      PREAMBL(3) = 'rotated from n-ward to y-ward by wi_xwdywd'
      WRITE(13,'(A79)') PREAMBL
      WRITE( 6,'(A79)') PREAMBL
C
C     PROCESS ALL THE WIND RECORDS.
C
      DO 810 KREC= 1,999999
C
C       READ THE INPUT WIND STRESSES.
C
        READ(10,'(A)',IOSTAT=IOS) CLINE_X
        IF     (IOS.NE.0) THEN
          EXIT !probably end of file
        ENDIF
        WRITE(6,'(A)') CLINE_X
        READ(11,'(A)') CLINE_Y
        WRITE(6,'(A)') CLINE_Y
        CALL ZAIORD(TXM,MSK,.FALSE., XMIN,XMAX, 10)
        CALL ZAIORD(TYM,MSK,.FALSE., XMIN,XMAX, 11)
C
C ---   Rotate from Xward and Yward to Eastward and Northward
C
        DO J= 1,JDM
          DO I= 1,IDM
            IF (TXM(I,J).NE.SPVAL) THEN
              TXIJ = TXM(I,J)
              TYIJ = TYM(I,J)
              TXM(I,J) = COS(PANG(I,J))*TXIJ +
     &                   SIN(PANG(I,J))*TYIJ
              TYM(I,J) = COS(PANG(I,J))*TYIJ -
     &                   SIN(PANG(I,J))*TXIJ
            ENDIF
          ENDDO
        ENDDO
C
        CALL ARCUPD(TXM,  IDM,JDM, LARCTIC,.FALSE.)
        CALL ARCUPD(TYM,  IDM,JDM, LARCTIC,.FALSE.)
C
C       WRITE OUT ROTATED FIELDS.
C
        CALL ZAIOWR(TXM,MSK,.FALSE., XMIN,XMAX, 12, .FALSE.)
        WRITE(12,4112) CLINE_X(1:48),XMIN,XMAX
        CALL ZHFLSH(12)
        WRITE( 6,4112) CLINE_X(1:48),XMIN,XMAX
C
        CALL ZAIOWR(TYM,MSK,.FALSE., XMIN,XMAX, 13, .FALSE.)
        WRITE(13,4112) CLINE_Y(1:48),XMIN,XMAX
        CALL ZHFLSH(13)
        WRITE( 6,4112) CLINE_Y(1:48),XMIN,XMAX
        CALL ZHFLSH( 6)
  810 CONTINUE
C
      CALL ZAIOCL(10)
      CLOSE( UNIT=10)
      CALL ZAIOCL(11)
      CLOSE( UNIT=11)
      CALL ZAIOCL(12)
      CLOSE( UNIT=12)
      CALL ZAIOCL(13)
      CLOSE( UNIT=13)
      STOP
C
 4112 FORMAT(A,1P2E16.7)
C
C     END OF PROGRAM WI_EWDNWD
      END
