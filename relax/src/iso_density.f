      PROGRAM ISO_DENISITY
      USE MOD_ZA  ! HYCOM array I/O interface
      IMPLICIT NONE
C
C     BLKDAT VARIABLES.
C
      INTEGER      KDM,IF,IL,JF,JL
      REAL*4       SIGMA(99)
C
C     OUTPUT ARRAYS
C
      REAL*4,  ALLOCATABLE :: AMSK(:,:),R(:,:,:)
C
C     OTHER ARRAYS.
C
      INTEGER, ALLOCATABLE :: MSK(:,:)
C
C**********
C*
C 1)  CREATE A SPACIALLY VARYING SET OF TARGET DENSITIES FOR HYCOM.
C
C      ONLY FOR USE WITH HYCOM 2.1.09 OR LATER.
C
C 2)  PARAMETERS:
C
C     MODEL GRID SPECIFICATION (W.R.T. PRESSURE GRID):
C
C        IDM    = 1ST DIMENSION OF MAJOR (HYCOM) MODEL ARRAYS
C        JDM    = 2ND DIMENSION OF MAJOR (HYCOM) MODEL ARRAYS
C
C 3)  INPUT:
C        ON UNIT 61A: MASK FILE
C        ON UNIT 99:  THE REQUIRED PROFILES.
C                     FOLLOWED BY THE REQUIRED PROFILE
C     OUTPUT:
C        ON UNIT 10:  TARGET DENSITIES FOR EACH LAYER
C        ON UNIT 10A: TARGET DENSITIES FOR EACH LAYER
C
C 4)  THE NON-DEFAULT PROFILES WILL ONLY BE USED WHERE MASK IS NOT
C     DATA-VOID.  IF NO SPECIAL MASK IS REQUIRED, JUST INPUT A COPY
C     OF THE BATHYMETRY ON UNIT 61A.
C
C 5)  ALAN J. WALLCRAFT,  NAVAL RESEARCH LABORATORY, MARCH 2004.
C*
C**********
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      LOGICAL LCONST
      INTEGER I,II,J,K
      REAL*4  BLK,SI,SJ,SIGMAQ(4),XMAX,XMIN
C
      CALL XCSPMD
      CALL ZAIOST
C
C     BLKDAT INPUT.
C
      CALL ZHOPEN(99, 'FORMATTED', 'OLD', 0)
C
C --- 'idm   ' = longitudinal array size
C --- 'jdm   ' = latitudinal  array size
C --- 'kdm   ' = layer        array size
C
      WRITE(6,*)
      CALL BLKINI(I,     'idm   ')
      CALL BLKINI(J,     'jdm   ')
      CALL BLKINI(KDM,   'kdm   ')
C
      IF     (I.NE.IDM) THEN
        WRITE(6,*)
        WRITE(6,*) 'ERROR - wrong IDM'
        WRITE(6,*)
        CALL ZHFLSH(6)
        STOP
      ELSEIF (J.NE.JDM) THEN
        WRITE(6,*)
        WRITE(6,*) 'ERROR - wrong JDM'
        WRITE(6,*)
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
C ---   'sigma ' = default layer densities (sigma units)
C
      DO K=1,KDM
        CALL BLKINR(SIGMA(K),'sigma ','(a6," =",f10.4," sigma")')
      ENDDO
C
C --- LOOP THROUGH SUB-REGIONS WITH DIFFERENT TARGET DENSITIES
C
      LCONST = .TRUE.
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
        IF     (LCONST) THEN  !first loop
          ALLOCATE(   R(IDM,JDM,KDM) )
          ALLOCATE( AMSK(IDM,JDM)    )
          ALLOCATE(  MSK(IDM,JDM)    )
          DO K=1,KDM
            R(:,:,K) = SIGMA(K)
          ENDDO
          CALL ZAIOPN('OLD', 61)
          CALL ZAIORD(AMSK,MSK,.FALSE., XMIN,XMAX, 61)
          CALL ZAIOCL(61)
          LCONST = .FALSE.
        ENDIF
C
C ---   'sigma ' = sub-region layer densities (sigma units), or
C ---   'sigmaA' = bottom left  layer density (sigma units), optional
C ---   'sigmaB' = bottom right layer density (sigma units), optional
C ---   'sigmaC' = top    right layer density (sigma units), optional
C ---   'sigmaD' = top    left  layer density (sigma units), optional
C
        DO K=1,KDM
          CALL BLKINR2(BLK,II,'sigma ','(a6," =",f10.4," sigma")',
     &                        'sigmaA','(a6," =",f10.4," sigma")')
          IF     (II.EQ.1) THEN !sigma
            SIGMA(K) = BLK
            DO J= JF,JL
              DO I= IF,IL
                IF     (AMSK(I,J).NE.SPVAL) THEN
                  R(I,J,K) = SIGMA(K)
                ENDIF
              ENDDO !i
            ENDDO !j
          ELSE !sigma1
            SIGMAQ(1) = BLK
            CALL BLKINR(SIGMAQ(2),'sigmaB','(a6," =",f10.4," sigma")')
            CALL BLKINR(SIGMAQ(3),'sigmaC','(a6," =",f10.4," sigma")')
            CALL BLKINR(SIGMAQ(4),'sigmaD','(a6," =",f10.4," sigma")')
            SIGMA(K) = 0.25*SUM(SIGMAQ(:))
            DO J= JF,JL
              SJ = REAL(J-JF)/REAL(JL-JF)  !0.0 to 1.0
              DO I= IF,IL
                IF     (AMSK(I,J).NE.SPVAL) THEN
                  SI = REAL(I-IF)/REAL(IL-IF)  !0.0 to 1.0
                  R(I,J,K) = SIGMAQ(1) + 
     &                       (1.0-SI)*(1.0-SJ)*(SIGMAQ(1)-SIGMAQ(1)) + 
     &                            SI *(1.0-SJ)*(SIGMAQ(2)-SIGMAQ(1)) + 
     &                            SI *     SJ *(SIGMAQ(3)-SIGMAQ(1)) + 
     &                       (1.0-SI)*     SJ *(SIGMAQ(4)-SIGMAQ(1))
                ENDIF
              ENDDO !i
            ENDDO !j
          ENDIF !ii
        ENDDO !k
      ENDDO !sub-region loop
      CLOSE(UNIT=99)
C
      IF     (LCONST) THEN
C
C ---   CONSTANT TARGET DENSITIES.
C
        CALL ZHOPEN(10, 'FORMATTED', 'NEW', 0)
        DO K=1,KDM
          WRITE(10,4102) K,SIGMA(K),SIGMA(K)
          WRITE( 6,4102) K,SIGMA(K),SIGMA(K)
        ENDDO
        CLOSE(UNIT=10)
      ELSE
C
C ---   VARYING TARGET DENSITIES.
C
        CALL ZHOPEN(10, 'FORMATTED', 'NEW', 0)
        CALL ZAIOPN('NEW', 10)
        DO K=1,KDM
          CALL ZAIOWR(R(1,1,K),MSK,.FALSE., XMIN,XMAX, 10, .FALSE.)
          WRITE(10,4102) K,XMIN,XMAX
          WRITE( 6,4102) K,XMIN,XMAX
        ENDDO
        CALL ZAIOCL(10)
        CLOSE(UNIT=10)
      ENDIF
      STOP
 4102 FORMAT('target_density: layer,range = ',I2.2,2F10.4)
      END
