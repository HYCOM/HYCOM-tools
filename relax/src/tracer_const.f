      PROGRAM TRACER_CONST
      USE MOD_ZA  ! HYCOM array I/O interface
      IMPLICIT NONE
C
C     BLKDAT VARIABLES.
C
      INTEGER      MONTH,KDM,NTRACR,IF,IL,JF,JL
      REAL*4       TRACER(99,19)
C
C     OUTPUT ARRAYS
C
      REAL*4,  ALLOCATABLE :: R(:,:,:,:)
C
C     OTHER ARRAYS.
C
      INTEGER, ALLOCATABLE :: MSK(:,:)
C
C**********
C*
C 1)  CREATE A SPACIALLY VARYING SET OF CLIMATOLOGICAL TRACER FIELDS.
C
C 2)  INPUT:
C        ON UNIT 99:  THE REQUIRED TRACER VALUES.
C     OUTPUT:
C        ON UNIT 10:  TRACER FOR EACH LAYER AND TRACER NUMBER
C        ON UNIT 10A: TRACER FOR EACH LAYER AND TRACER NUMBER
C
C 3)  ALAN J. WALLCRAFT,  NAVAL RESEARCH LABORATORY, DEEMBER 2003.
C*
C**********
C
      LOGICAL   LCONST
      INTEGER   I,J,K,KTR
      REAL*4    XMAX,XMIN
      CHARACTER PREAMBL(5)*79
C
      CALL XCSPMD
C
C     BLKDAT INPUT.
C
      CALL ZHOPEN(99, 'FORMATTED', 'OLD', 0)
C
C --- 'month'  = month (1 to 12)
C --- 'idm   ' = longitudinal array size
C --- 'jdm   ' = latitudinal  array size
C --- 'kdm   ' = layer        array size
C --- 'ntracr' = number of tracers
C
      WRITE(6,*)
      CALL BLKINI(MONTH, 'month ')
      WRITE(6,*)
      CALL BLKINI(I,     'idm   ')
      CALL BLKINI(J,     'jdm   ')
      CALL BLKINI(KDM,   'kdm   ')
      CALL BLKINI(NTRACR,'ntracr')
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
      ELSEIF (NTRACR.LT.1 .OR. NTRACR.GT.19) THEN
        WRITE(6,*)
        WRITE(6,*) 'ERROR - NTRACR out of range 1-19'
        WRITE(6,*)
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
C ---   'tracer' = default layer tracer value
C
      DO KTR= 1,NTRACR
        DO K=1,KDM
          CALL BLKINR(TRACER(K,KTR),'tracer','(a6," =",f10.4)')
        ENDDO
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
          ALLOCATE(   R(IDM,JDM,KDM,NTRACR) )
          ALLOCATE( MSK(IDM,JDM)            )
          DO KTR= 1,NTRACR
            DO K=1,KDM
              R(:,:,K,KTR) = TRACER(K,KTR)
            ENDDO
          ENDDO
          LCONST = .FALSE.
        ENDIF
C
C ---   'sigma ' = sub-region layer densities (sigma units)
C
        DO KTR= 1,NTRACR
          DO K=1,KDM
            CALL BLKINR(TRACER(K,KTR),'tracer','(a6," =",f10.4)')
            R(IF:IL,JF:JL,K,KTR) = TRACER(K,KTR)
          ENDDO
        ENDDO
      ENDDO !sub-region loop
      CLOSE(UNIT=99)
C
      IF     (LCONST) THEN
C
C ---   CONSTANT TRACERS
C
        CALL ZHOPEN(10, 'FORMATTED', 'NEW', 0)
C
        PREAMBL(1) = 'Tracer Climatology via tracer_const'
        PREAMBL(2) = ' '
        PREAMBL(3) = ' '
        PREAMBL(4) = ' '
        WRITE(PREAMBL(5),'(A,2I5,2I3)')
     +          'idm,jdm,kdm,ntracr =',
     +         IDM,JDM,KDM,NTRACR
        WRITE(10,4101) PREAMBL
        WRITE(6,*)
        WRITE(6, 4101) PREAMBL
        WRITE(6,*)
C
        DO KTR= 1,NTRACR
          DO K=1,KDM
            WRITE(10,4102) MONTH,K,REAL(KTR),TRACER(K,KTR),
     &                                       TRACER(K,KTR)
            WRITE( 6,4102) MONTH,K,REAL(KTR),TRACER(K,KTR),
     &                                       TRACER(K,KTR)
          ENDDO
        ENDDO
        CLOSE(UNIT=10)
      ELSE
C
C ---   VARYING TARGET DENSITIES.
C
        CALL ZHOPEN(10, 'FORMATTED', 'NEW', 0)
        CALL ZAIOST
        CALL ZAIOPN('NEW', 10)
C
        PREAMBL(1) = 'Tracer Climatology via tracer_const'
        PREAMBL(2) = ' '
        PREAMBL(3) = ' '
        PREAMBL(4) = ' '
        WRITE(PREAMBL(5),'(A,2I5,2I3)')
     +          'idm,jdm,kdm,ntracr =',
     +         IDM,JDM,KDM,NTRACR
        WRITE(10,4101) PREAMBL
        WRITE(6,*)
        WRITE(6, 4101) PREAMBL
        WRITE(6,*)
C
        DO KTR= 1,NTRACR
          DO K=1,KDM
            CALL ZAIOWR(R(1,1,K,KTR),MSK,.FALSE.,XMIN,XMAX,10,.FALSE.)
            WRITE(10,4102) MONTH,K,REAL(KTR),XMIN,XMAX
            WRITE( 6,4102) MONTH,K,REAL(KTR),XMIN,XMAX
          ENDDO
        ENDDO
        CALL ZAIOCL(10)
        CLOSE(UNIT=10)
      ENDIF
      STOP
 4101 FORMAT(A79)
 4102 FORMAT('trcr: month,layer,ktrc,range = ',I2.2,I4.2,F7.3,1P2E16.7)
      END
