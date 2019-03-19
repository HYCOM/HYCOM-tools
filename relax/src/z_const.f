      PROGRAM Z_CONST
      USE MOD_ZA  ! HYCOM array I/O interface
      IMPLICIT NONE
C
C     NAMELIST
C
      CHARACTER CTITLE*79,CTYPE*21
      INTEGER   KZ
      REAL*4    ZLEV(1000),TRCZ(1000)
      NAMELIST/ZCONST/ CTITLE,CTYPE,KZ,ZLEV,TRCZ
C
C**********
C*
C 1)  CREATE A CONSTANT Z-LEVEL TRACER FIELD.
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
C     /ZCONST/             
C        CTITLE - 1ST LINE OF OUTPUT PREAMBL
C        KZ     - NUMBER OF Z-LEVELS IN CLIMATOLOGY
C        ZLEV   - Z-LEVEL DEPTHS (1ST 0, ZLEV(K+1)>ZLEV(K))
C        TRCZ   - TRACER VALUES (ARBITRARY)
C
C        NOTE THAT THE PROFILE WILL BE TREATED AS A FINITE VOLUME,
C        WITH ZLEV THE "REPRESENTATIVE" DEPTH OF EACH CELL.
C
C 3)  INPUT:
C        ON UNIT  5:  NAMELIST INPUT
C     OUTPUT:
C        ON UNIT 71:  TRACER Z-LEVEL CLIM FILE
C
C 4)  ALAN J. WALLCRAFT,  NAVAL RESEARCH LABORATORY, NOVEMBER 2004.
C*
C**********
C
      CHARACTER PREAMBL(5)*79,CSAVE*21
      INTEGER   K,L
      REAL*4    ZOLD
C
      CALL XCSPMD  !input idm,jdm
C                               
C     NAMELIST INPUT.           
C                    
      KZ = 0
      CTITLE = ' '
*     CTYPE  = '123456789012345678901'
      CTYPE  = '               tracer'
      DO K= 1,1000    
        ZLEV(K) = 0
        TRCZ(K) = 0.0
      ENDDO         
      READ( 5,ZCONST)
      WRITE(6,ZCONST)
      CALL ZHFLSH(6)
C
      IF     (KZ.LE.0 .OR. KZ.GT.1000) THEN
        WRITE(6,'(/A/)') 'ERROR - KZ MUST BE BETWEEN 1 AND 1000'
        CALL ZHFLSH(6)
        STOP
      ENDIF
      ZOLD = -1.0
      DO K= 1,KZ
        IF     (ZLEV(K).LE.ZOLD) THEN
          WRITE(6,'(/A,I4/)') 'ERROR  -  ZLEV(K) < ZLEV(K-1);  K =',K
          CALL ZHFLSH(6)
          STOP
        ENDIF
        ZOLD = ZLEV(K)
      ENDDO         
C
      L     = LEN_TRIM(CTYPE)
      CSAVE = CTYPE
      CTYPE = ' '
      CTYPE(22-L:21) = CSAVE(1:L)
C
C     Z-LEVEL CLIMATOLOGY OUTPUT.
C
      CALL ZHOPEN(71, 'FORMATTED', 'NEW', 0)
      PREAMBL(1) = CTITLE
      PREAMBL(2) = ' '
      PREAMBL(3) = ' '
      PREAMBL(4) = ' '
      WRITE(PREAMBL(5),'(A,2I5)')
     +        'i/jdm =',
     +       IDM,JDM
      WRITE(71,4101) PREAMBL
      WRITE( 6,*)
      WRITE( 6,4101) PREAMBL
      CALL ZHFLSH(6)
      DO K= 1,KZ
        WRITE(71,4102) TRIM(CTYPE),ZLEV(K),TRCZ(K),TRCZ(K)
        WRITE( 6,4102) TRIM(CTYPE),ZLEV(K),TRCZ(K),TRCZ(K)
        CALL ZHFLSH(6)
      ENDDO !k
      CLOSE(UNIT=71)
C
C     SUMMARY.
C
      WRITE(6,6400) KZ
      CALL ZHFLSH(6)
      STOP
C
 4101 FORMAT(A79)
 4102 FORMAT(A,': depth,range = ',F7.1,1P2G16.6)
 6400 FORMAT(I5,' LEVEL CLIMATOLOGY COMPLETED.')
C     END OF PROGRAM Z_CONST.
      END
