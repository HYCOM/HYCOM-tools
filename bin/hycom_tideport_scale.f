      PROGRAM HYCOM_TIDEPORT_SCALE
      IMPLICIT NONE 
C
C  hycom_tideport_scale - Usage:  
C            hycom_tideport_scale ports_u.input ports.depth ports_unew.input
C                 on output ports_unew.input is a copy of ports_u.input
C                 with velocities scaled by inverse of the 6th column of 
C                 ports.depth which is typically model_depth/tpxo_depth
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraf, NRL, July 2013.
C
      INTEGER       IARGC
      INTEGER       NARG
C
      CHARACTER*240 CFILEI,CFILED,CFILEO
      CHARACTER*240 CLINE
      REAL          DUM,SCL,U(16)
      INTEGER       I,IDUM,IOS
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.3) THEN
        CALL GETARG(1,CFILEI)
        CALL GETARG(2,CFILED)
        CALL GETARG(3,CFILEO)
      ELSE
        WRITE(6,'(2a)')
     +  'Usage: hycom_tideport_scale ',
     +     'ports_u.input ports.depth ports_unew.input'
        CALL EXIT(1)
      ENDIF
C
C     OPEN ALL FILES.
C
      OPEN(UNIT=11, FILE=CFILEI, FORM='FORMATTED', STATUS='OLD',
     +     IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',TRIM(CFILEI)
        WRITE(6,*) 'ios   = ',ios
        CALL EXIT(3)
      ENDIF
      OPEN(UNIT=12, FILE=CFILED, FORM='FORMATTED', STATUS='OLD',
     +     IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',TRIM(CFILED)
        WRITE(6,*) 'ios   = ',ios
        CALL EXIT(4)
      ENDIF
      OPEN(UNIT=21, FILE=CFILEO, FORM='FORMATTED', STATUS='NEW',
     +     IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',TRIM(CFILEO)
        WRITE(6,*) 'ios   = ',ios
        CALL EXIT(5)
      ENDIF
C
C     COPY FIRST THREE LINES
C
      READ( 11,'(a)') CLINE
      WRITE(21,'(a)') TRIM(CLINE)
      READ( 11,'(a)') CLINE
      WRITE(21,'(a)') TRIM(CLINE)//'    scaled by:  '//TRIM(CFILED)
      READ( 11,'(a)') CLINE
      WRITE(21,'(a)') TRIM(CLINE)
C
C     SKIP HEADER OF CFILED
C
      DO I= 1,5
        READ(12,'(a)') CLINE
        IF     (CLINE(1:1).NE.'#') THEN
          WRITE(6,'(2a)') 'ERROR READING HEADER OF ',TRIM(CFILED)
          CALL EXIT(6)
        ENDIF
      ENDDO !i
C
C     READ THE PORT LINES
C
      DO I= 1,999999
        READ(11,'(a)',IOSTAT=IOS) CLINE
        IF     (IOS.NE.0) THEN
          IF     (I.EQ.1) THEN
            WRITE(6,'(2a)') 'ERROR READING ',TRIM(CFILEI)
            CALL EXIT(6)
          ENDIF
          EXIT
        ENDIF
        READ(CLINE,*) DUM,DUM,U(:)
C
        READ(12,*) IDUM,DUM,DUM,DUM,DUM,SCL
        SCL = MAX( 0.25, MIN( 2.0, SCL ) )
C
        WRITE(21,'(a,16f8.3)') CLINE(1:19),U(:)/SCL
      ENDDO !i
      CLOSE(21)
      CLOSE(11)
      CLOSE(12)
      END
