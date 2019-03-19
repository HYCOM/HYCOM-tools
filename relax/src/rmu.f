      PROGRAM RLXMSK
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
C     NAMELIST INPUT.
C
      CHARACTER*79   CTITLE
      INTEGER        IF(999),IL(999),JF(999),JL(999)
      REAL*4         EFOLD(999)
      NAMELIST/MASK/ CTITLE,IF,IL,JF,JL,EFOLD
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
C     THE ALLOWED RANGE FOR IF,IL IS 1 TO IDM-1
C     THE ALLOWED RANGE FOR JF,JL IS 1 TO JDM-1
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
      INTEGER I,ISF,ISL,ISEC,J,K
      REAL*4  RVAL,RVMIN
C
      CHARACTER*1 C(-1:9)
      DATA C / '*', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' /
C
      CALL XCSPMD
      ALLOCATE(   RMU(IDM,JDM) )
      ALLOCATE(  IRMU(IDM,JDM) )
      ALLOCATE( DEPTH(IDM,JDM) )
C
C     NAMELIST INPUT.
C
      DO K= 1,999
        IF(   K) = 0
        IL(   K) = 0
        JF(   K) = 0
        JL(   K) = 0
        EFOLD(K) = 0.0
      ENDDO
      READ( 5,MASK)
      WRITE(6,MASK)
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
C     INITIALIZE SPECIFIED PATCHES.
C
      RVMIN = 1.E20
      DO K= 1,999
        IF     (MIN(IF(K),IL(K),JF(K),JL(K)).GT.0   .AND.
     +                              EFOLD(K).GT.0.0      ) THEN
          IF     (IF(K).GT.IL(K) .OR.
     +            IF(K).LT.1     .OR.
     +            IL(K).GT.IDM-1     ) THEN
            WRITE(6,9000) K,IF(K),IL(K),IDM-1
            CALL ZHFLSH(6)
            STOP
          ENDIF
          IF     (JF(K).GT.JL(K) .OR.
     +            JF(K).LT.1     .OR.
     +            JL(K).GT.JDM-1     ) THEN
            WRITE(6,9100) K,JF(K),JL(K),JDM-1
            CALL ZHFLSH(6)
            STOP
          ENDIF
C
          RVAL = 1.0/(EFOLD(K)*86400.0)
          DO J= JF(K),JL(K)
            DO I= IF(K),IL(K)
              IF     (DEPTH(I,J).NE.SPVAL) THEN
                RMU(I,J) = MAX(RMU(I,J),RVAL)
              ENDIF
            ENDDO
          ENDDO
          RVMIN = MIN( RVMIN, RVAL )
        ENDIF
      ENDDO
C
      WRITE(6,*)
      DO K= 1,999
        IF     (MIN(IF(K),IL(K),JF(K),JL(K)).GT.0   .AND.
     +                              EFOLD(K).GT.0.0      ) THEN
          RVAL = 1.0/(EFOLD(K)*86400.0)
          WRITE(6,6200) IF(K),IL(K),JF(K),JL(K),
     +                  EFOLD(K),RVAL,RVAL/RVMIN
          CALL ZHFLSH(6)
        ENDIF
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
 6200 FORMAT('IF,IL,JF,JL =',4I5,'   EFOLD,RVAL =',F6.1,1PE12.3,
     +   '  (',0PF4.1,')')
 9000 FORMAT('ERROR - ILLEGAL IF OR IL FOR K =',I3,
     +       '   MUST HAVE 1<=IF(K)<=IL(K)<=IDM-1' /
     +       'IF(K),IL(K),IDM-1 = ',3I5 /)
 9100 FORMAT('ERROR - ILLEGAL JF OR JL FOR K =',I3,
     +       '   MUST HAVE 1<=JF(K)<=JL(K)<=JDM-1' /
     +       'JF(K),JL(K),JDM-1 = ',3I5 /)
      END
