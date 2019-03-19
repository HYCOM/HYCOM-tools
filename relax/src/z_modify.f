      PROGRAM MODIFY
      USE MOD_ZA  ! HYCOM array I/O interface
      IMPLICIT NONE
C
C     I/O VARIABLES.
C
      CHARACTER PREAMBL_T(5)*79,PREAMBL_R(5)*79,CTITLE*79,CLINE*80
      INTEGER   IF,IL,JF,JL
      REAL*4    RZK,TZK,HMINA,HMINB,HMAXA,HMAXB
C
C     CLIM ARRAYS.
C
      REAL*4               :: ZLEV(999)
      REAL*4,  ALLOCATABLE :: TZ(:,:,:),RZ(:,:,:)
      INTEGER, ALLOCATABLE :: MSK(:,:)
C
C**********
C*
C 1)  MODIFY A Z-LEVEL CLIMATOLOGY ON THE HYCOM REGION GRID.
C
C 2)  PARAMETERS:
C
C     MODEL GRID SPECIFICATION (W.R.T. PRESSURE GRID):
C
C        IDM    = 1ST DIMENSION OF MAJOR (HYCOM) MODEL ARRAYS
C        JDM    = 2ND DIMENSION OF MAJOR (HYCOM) MODEL ARRAYS
C        KZ     = NUMBER OF Z-LEVELS IN CLIMATOLOGY
C
C 3)  INPUT:
C        ON UNIT 61:  DENS. Z-LEVEL CLIM FILE
C        ON UNIT 62:  TEMP. Z-LEVEL CLIM FILE
C        ON UNIT 99:  THE MODIFIED CLIMATOLOGICAL VALUES.
C     OUTPUT:
C        ON UNIT 71:  DENS. Z-LEVEL CLIM FILE
C        ON UNIT 72:  TEMP. Z-LEVEL CLIM FILE
C
C 4)  ALAN J. WALLCRAFT,  NAVAL RESEARCH LABORATORY, MAY 2004.
C*
C**********
C
      INTEGER I,J,K,KSIGMA,KZ
      REAL*4  XAVE,XMAX,XMIN
C
      CALL XCSPMD  !input idm,jdm
C
      CALL ZHOPEN(61, 'FORMATTED', 'OLD', 0)
      READ(61,*) !5-line header
      READ(61,*) 
      READ(61,*) 
      READ(61,*) 
      READ(61,*) 
      DO K= 1,999
        READ(61,*,END=100) !one line per level
      ENDDO
  100 CONTINUE
      KZ = K-1
      write(6,*) 'kz = ',kz
      CLOSE(61)
C
      ALLOCATE(  TZ(IDM,JDM,KZ+1) )
      ALLOCATE(  RZ(IDM,JDM,KZ+1) )
      ALLOCATE( MSK(IDM,JDM) )
C
C     Z-LEVEL CLIMATOLOGY INPUT.
C
      CALL ZAIOST
C
      CALL ZAIOPN('OLD', 61)
      CALL ZHOPEN(61, 'FORMATTED', 'OLD', 0)
      READ (61,'(A79)') PREAMBL_R
      WRITE(6, '(/(1X,A79))') PREAMBL_R
      DO K= 1,KZ
        READ (61,'(A)') CLINE
        I = INDEX(CLINE,'=')
        READ (CLINE(I+1:),*) ZLEV(K),HMINB,HMAXB
        IF     (HMINB.EQ.HMAXB) THEN  ! constant field
          DO J= 1,JDM
            DO I= 1,IDM
              RZ(I,J,K) = HMINB
            ENDDO
          ENDDO
          CALL ZAIOSK(61)
        ELSE
          CALL ZAIORD(RZ(1,1,K),MSK,.FALSE., HMINA,HMAXA, 61)
          IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &            ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
            WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &       'error - .a and .b density files not consistent:',
     &       '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &       '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
            CALL ZHFLSH(6)
            STOP
          ENDIF
        ENDIF
      ENDDO
      CLOSE(UNIT=61)
      CALL ZAIOCL(61)
C
      I = INDEX(CLINE,'sigma-0')
      IF     (I.GT.0) THEN
        KSIGMA = 0
      ELSE
        KSIGMA = 2
      ENDIF
C
      CALL ZAIOPN('OLD', 62)
      CALL ZHOPEN(62, 'FORMATTED', 'OLD', 0)
      READ (62,'(A79)') PREAMBL_T
      WRITE(6, '(/(1X,A79))') PREAMBL_T
      DO K= 1,KZ
        READ (62,'(A)') CLINE
        I = INDEX(CLINE,'=')
        READ (CLINE(I+1:),*) ZLEV(K),HMINB,HMAXB
        IF     (HMINB.EQ.HMAXB) THEN  ! constant field
          DO J= 1,JDM
            DO I= 1,IDM
              TZ(I,J,K) = HMINB
            ENDDO
          ENDDO
          CALL ZAIOSK(62)
        ELSE
          CALL ZAIORD(TZ(1,1,K),MSK,.FALSE., HMINA,HMAXA, 62)
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
      ENDDO
      CLOSE (UNIT=62)
      CALL ZAIOCL(62)
C
C     BLKDAT INPUT.
C
      CALL ZHOPEN(99, 'FORMATTED', 'OLD', 0)
C
C --- 'ctitle' = climatology modification title
C
      WRITE(6,*)
      READ(99,'(A79)') CTITLE
      WRITE(6,'(A79)') CTITLE
C            
C --- LOOP THROUGH SUB-REGIONS WITH DIFFERENT TARGET DENSITIES
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
        DO K= 1,KZ 
          CALL BLKINR(RZK,'rzk   ','(a6," =",f10.4)')
          CALL BLKINR(TZK,'tzk   ','(a6," =",f10.4)')
          RZ(IF:IL,JF:JL,K) = RZK
          TZ(IF:IL,JF:JL,K) = TZK
        ENDDO !k
      ENDDO !sub-region loop
      CLOSE(UNIT=99)        
C
C     CLIMATOLOGY OUTPUT.
C
      CALL ZHOPEN(71, 'FORMATTED', 'NEW', 0)
      CALL ZHOPEN(72, 'FORMATTED', 'NEW', 0)
      PREAMBL_R(4) = CTITLE
      PREAMBL_T(4) = CTITLE
      WRITE(71,4101) PREAMBL_R
      WRITE(72,4101) PREAMBL_T
C
      CALL ZAIOPN('NEW', 71)
      CALL ZAIOPN('NEW', 72)
C
      DO K= 1,KZ
        CALL ZAIOWR(RZ(1,1,K),MSK,.FALSE., XMIN,XMAX, 71, .FALSE.)
        IF     (KSIGMA.EQ.0) THEN
          WRITE(71,4102) '              sigma-0',ZLEV(K),XMIN,XMAX
        ELSE
          WRITE(71,4102) '              sigma-2',ZLEV(K),XMIN,XMAX
        ENDIF
C
        CALL ZAIOWR(TZ(1,1,K),MSK,.FALSE., XMIN,XMAX, 72, .FALSE.)
        WRITE(72,4102) 'potential temperature',ZLEV(K),XMIN,XMAX
C
        WRITE(6,6300) K,ZLEV(K)
        CALL ZHFLSH(6)
      ENDDO !k
C
      CALL ZAIOCL(71)
      CLOSE( UNIT=71)
      CALL ZAIOCL(72)
      CLOSE( UNIT=72)
C
C     SUMMARY.
C
      WRITE(6,6400) KZ
      CALL ZHFLSH(6)
      STOP
C
 4101 FORMAT(A79)
 4102 FORMAT(A,': depth,range = ',F7.1,1P2E16.7)
 5000 FORMAT(A40)
 5500 FORMAT(6E13.6)
 6200 FORMAT(/ 1X,'MIN,MAX I COORDS = ',F8.2,',',F8.2 
     +       / 1X,'MIN,MAX J COORDS = ',F8.2,',',F8.2 /)
 6300 FORMAT(10X,'WRITING CLIM RECORD',I3,'     ZLEV =',F7.1 /)
 6400 FORMAT(I5,' LEVEL CLIMATOLOGY COMPLETED.')
 8100 FORMAT(1X,A,': min=',F9.2,' ave=',F9.2,' max=',F9.2,
     +   '   (k,sigma =',i3,F7.2,')')
C     END OF PROGRAM MODIFY.
      END
