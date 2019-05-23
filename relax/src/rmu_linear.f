      PROGRAM RMU_LINEAR
      USE MOD_ZA  ! HYCOM array I/O interface
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
C     OUTPUT ARRAYS
C
      REAL*4,  ALLOCATABLE :: RMU(:,:)
C
C     OTHER ARRAYS.
C
      INTEGER, ALLOCATABLE :: IRMU(:,:)
      REAL*4,  ALLOCATABLE :: DEPTH(:,:)
C
C     OTHER VARIABLES.
C
      REAL*4    RMUMIN,RMUMAX,HMINA,HMINB,HMAXA,HMAXB
      CHARACTER PREAMBL(5)*79,CLINE*80
C
C     BLKDAT VARIABLES.
C
      CHARACTER*79 CTITLE
      INTEGER      IF,IL,JF,JL
      REAL*4       EFOLD,EFOLDQ(4)
C
C**********
C*
C 1)  CREATE A HYCOM RELAXATION MASK
C
C 2)  PARAMETERS:
C
C     MODEL GRID SPECIFICATION (W.R.T. PRESSURE GRID):
C
C        IDM    = 1ST DIMENSION OF MAJOR (HYCOM) MODEL ARRAYS
C        JDM    = 2ND DIMENSION OF MAJOR (HYCOM) MODEL ARRAYS
C
C 3)  NAMELIST PARAMETERS
C
C     /MASK/
C        CTITLE      - 2ND LINE OF OUTPUT PREAMBL
C                       (1ST LINE IS ALWAYS 'Relaxation Mask')
C        IF,IL,JF,JL - ARRAY BOX WHERE EFOLD RELAXATION IS APPLIED
C                       = 4*0.0; END OF BOX LIST
C        EFOLD       - RELAXATION E-FOLDING TIME IN DAYS
C                       = 0.0; NO RELAXATION
C
C     THE ALLOWED RANGE FOR IF,IL IS 1 TO IDM
C     THE ALLOWED RANGE FOR JF,JL IS 1 TO JDM
C
C 4)  INPUT:
C        ON UNIT  5: NAMELIST INPUT
C        ON UNIT 51: BATHYMETRY FILE
C     OUTPUT:
C        ON UNIT 21: RELAXATION MASK
C
C 5)  ALAN J. WALLCRAFT,  NAVAL RESEARCH LABORATORY, JUNE 2000.
C*
C**********
C
      INTEGER I,II,ISF,ISL,ISEC,J,K
      REAL*4  BLK,SI,SJ,RVAL,RVMIN
C
      CHARACTER*1 C(-1:9)
      DATA C / '*', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' /
C
      RVMIN = HUGE(RVMIN)
C
      CALL XCSPMD
      ALLOCATE(   RMU(IDM,JDM) )
      ALLOCATE(  IRMU(IDM,JDM) )
      ALLOCATE( DEPTH(IDM,JDM) )
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
      CALL ZAIORD(DEPTH,IRMU,.FALSE., HMINA,HMAXA, 51)
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
C     INITIALIZE MASK TO ZERO
C
      DO J= 1,JDM
        DO I= 1,IDM
          RMU(I,J) = 0.0
        ENDDO
      ENDDO
C
C     BLKDAT INPUT.
C
      CALL ZHOPEN(99, 'FORMATTED', 'OLD', 0)
      READ(99,'(a)') CTITLE
C
      DO !sub-region loop
C
C ---   'if    ' = first i point of sub-region (<=0 to end)
C ---   'il    ' = last  i point of sub-region
C ---   'jf    ' = last  j point of sub-region
C ---   'jl    ' = last  j point of sub-region
C
        CALL BLKINI(IF,    'if    ')
        IF     (IF.LE.0) THEN
          EXIT !sub-region loop
        ENDIF
        CALL BLKINI(IL,    'il    ')
        CALL BLKINI(JF,    'jf    ')
        CALL BLKINI(JL,    'jl    ')
C
        IF     (IF.GT.IL .OR.
     +          IF.LT.1  .OR.
     +          IL.GT.IDM    ) THEN
          WRITE(6,9000) K,IF,IL,IDM
          CALL ZHFLSH(6)
          STOP
        ENDIF
        IF     (JF.GT.JL .OR.
     +          JF.LT.1  .OR.
     +          JL.GT.JDM    ) THEN
          WRITE(6,9100) K,JF,JL,JDM
          CALL ZHFLSH(6)
          STOP
        ENDIF
C
C ---   'efold ' = relaxation e-folding time in days (0.0=no relax)
C ---   'efoldA' = bottom left  e-folding time in days, optional
C ---   'efoldB' = bottom right e-folding time in days, optional
C ---   'efoldC' = top    right e-folding time in days, optional
C ---   'efoldD' = top    left  e-folding time in days, optional
C
        CALL BLKINR2(BLK,II,'efold ','(a6," =",f10.4," days")',
     &                      'efoldA','(a6," =",f10.4," days")')
        IF     (II.EQ.1) THEN !efold
          EFOLD = BLK
          RVAL  = 1.0/(EFOLD*86400.0)
          DO J= JF,JL
            DO I= IF,IL
              IF     (DEPTH(I,J).NE.SPVAL) THEN
                RMU(I,J) = MAX(RMU(I,J),RVAL)
              ENDIF
            ENDDO !i
          ENDDO !j
          RVMIN = MIN( RVMIN, RVAL )
        ELSE !efoldA
          EFOLDQ(1) = BLK
          CALL BLKINR(EFOLDQ(2),'efoldB','(a6," =",f10.4," days")')
          CALL BLKINR(EFOLDQ(3),'efoldC','(a6," =",f10.4," days")')
          CALL BLKINR(EFOLDQ(4),'efoldD','(a6," =",f10.4," days")')
          DO J= JF,JL
            SJ = REAL(J-JF)/REAL(JL-JF)  !0.0 to 1.0
            DO I= IF,IL
              IF     (DEPTH(I,J).NE.SPVAL) THEN
                SI = REAL(I-IF)/REAL(IL-IF)  !0.0 to 1.0
                EFOLD    = EFOLDQ(1) + 
     &                     (1.0-SI)*(1.0-SJ)*(EFOLDQ(1)-EFOLDQ(1)) +
     &                          SI *(1.0-SJ)*(EFOLDQ(2)-EFOLDQ(1)) +
     &                          SI *     SJ *(EFOLDQ(3)-EFOLDQ(1)) +
     &                     (1.0-SI)*     SJ *(EFOLDQ(4)-EFOLDQ(1)) 
                RVAL     = 1.0/(EFOLD*86400.0)
                RMU(I,J) = MAX(RMU(I,J),RVAL)
                RVMIN    = MIN(RVMIN,   RVAL)
              ENDIF
            ENDDO !i
          ENDDO !j
        ENDIF !ii
      ENDDO
C
C     WRITE OUT THE MASK
C
      PREAMBL(1) = 'Relaxation Mask'
      PREAMBL(2) = CTITLE
      PREAMBL(3) = ' '
      PREAMBL(4) = ' '
      WRITE(PREAMBL(5),'(A,2I5)')
     +        'i/jdm =',
     +       IDM,JDM
C
      CALL ZAIOPN('NEW', 21)
      CALL ZAIOWR(RMU,IRMU,.FALSE., RMUMIN,RMUMAX, 21, .FALSE.)
      CALL ZAIOCL(21)
C
      CALL ZHOPEN(21, 'FORMATTED', 'NEW', 0)
      WRITE(21,4101) PREAMBL
      WRITE(21,4102) '     rmu',RMUMIN,RMUMAX
      CLOSE(UNIT=21)
C
      WRITE(6, *)
      WRITE(6, 4101) PREAMBL
      WRITE(6, 4102) ' rmu',RMUMIN,RMUMAX
      WRITE(6, *)
      CALL ZHFLSH(6)
C
C     PRINTOUT THE MASK
C
      WRITE(6,6000) IDM,JDM,RVMIN
      ISEC = (IDM-1)/100 + 1
      DO K= 1,ISEC
        ISF = (K-1)*100 + 1
        ISL = MIN(IDM, ISF+100-1)
        WRITE(6,6050) ISF,ISL
        DO J= JDM,1,-1
          DO I= ISF,ISL
            IF     (DEPTH(I,J).EQ.SPVAL) THEN
              IRMU(I,J) = -1
            ELSEIF (RMU(I,J).EQ.0.0) THEN
              IRMU(I,J) =  0
            ELSE
              IRMU(I,J) =  MIN( 9, NINT( RMU(I,J)/RVMIN ) )
            ENDIF
          ENDDO
          WRITE(6,6100) J,(C(IRMU(I,J)),I=ISF,ISL)
        ENDDO
      ENDDO
      STOP
C
 4101 FORMAT(A79)
 4102 FORMAT(A,': range = ',1P2E16.7)
 6000 FORMAT(1H1 / 30X,'RELAXATION MASK FOR AN',I5,'  BY',I5,' MESH.'
     +           / 31X,'LAND = *,  OCEAN = NINT(RMU /',1PE9.2,')' )
 6050 FORMAT(/ / / 21X,'I =',I5,'  TO',I5,'  :' / /)
 6100 FORMAT(4X,'J =',I5,5X,10(10A1,1X))
 9000 FORMAT('ERROR - ILLEGAL IF OR IL FOR K =',I3,
     +       '   MUST HAVE 1<=IF<=IL<=IDM-1' /
     +       'IF,IL,IDM-1 = ',3I5 /)
 9100 FORMAT('ERROR - ILLEGAL JF OR JL FOR K =',I3,
     +       '   MUST HAVE 1<=JF<=JL<=JDM-1' /
     +       'JF,JL,JDM-1 = ',3I5 /)
      END
