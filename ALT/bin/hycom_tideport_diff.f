      PROGRAM HYCOM_TIDEPORT_DIFF
      IMPLICIT NONE 
C
C  hycom_tideport_diff - Usage:  
C            hycom_tideport_diff ports_1.input ports_2.input ports_new.input
C                 on output ports_new.input is the difference between
C                 ports_1.input and ports_2.input
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraf, NRL, august 2016.
C
      INTEGER       IARGC
      INTEGER       NARG
C
      CHARACTER*240 CFILE1,CFILE2,CFILEO
      CHARACTER*240 CLINE,CLINE2
      REAL          DUM,U1(16),U2(16)
      INTEGER       I,IDUM,IOS
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.3) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CFILE2)
        CALL GETARG(3,CFILEO)
      ELSE
        WRITE(6,'(2a)')
     +  'Usage: hycom_tideport_diff ',
     +     'ports_1.input ports_2.input ports_new.input'
        CALL EXIT(1)
      ENDIF
C
C     OPEN ALL FILES.
C
      OPEN(UNIT=11, FILE=CFILE1, FORM='FORMATTED', STATUS='OLD',
     +     IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',TRIM(CFILE1)
        WRITE(6,*) 'ios   = ',ios
        CALL EXIT(3)
      ENDIF
      OPEN(UNIT=12, FILE=CFILE2, FORM='FORMATTED', STATUS='OLD',
     +     IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',TRIM(CFILE2)
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
      READ( 12,'(a)') CLINE2
      WRITE(21,'(a)') TRIM(CLINE)//' - '//TRIM(CLINE2)
      READ( 11,'(a)') CLINE
      READ( 12,'(a)') CLINE2
      WRITE(21,'(a)') TRIM(CLINE)
      READ( 11,'(a)') CLINE
      READ( 12,'(a)') CLINE2
      WRITE(21,'(a)') TRIM(CLINE)
C
C     READ THE PORT LINES
C
      DO I= 1,999999
        READ(11,'(a)',IOSTAT=IOS) CLINE
        IF     (IOS.NE.0) THEN
          IF     (I.EQ.1) THEN
            WRITE(6,'(2a)') 'ERROR READING ',TRIM(CFILE1)
            CALL EXIT(6)
          ENDIF
          EXIT
        ENDIF
        READ(12,'(a)',IOSTAT=IOS) CLINE2
        IF     (IOS.NE.0) THEN
          IF     (I.EQ.1) THEN
            WRITE(6,'(2a)') 'ERROR READING ',TRIM(CFILE1)
            CALL EXIT(6)
          ENDIF
          EXIT
        ENDIF
        IF     (INDEX(CLINE, 'Site').NE.0 .OR.
     &          INDEX(CLINE2,'Site').NE.0     ) THEN
          WRITE(6,'(a,i5)') 'ERROR UNDEFINED PORT LOCATION',I
          CALL EXIT(9)
        ENDIF
        IF     (LEN_TRIM(CLINE) .EQ.147) THEN
          READ(CLINE,'(F9.4,F10.4,16F8.3)') DUM,DUM,U1(:)
          IF     (LEN_TRIM(CLINE2).EQ.147) THEN
            READ(CLINE2,'(F9.4,F10.4,16F8.3)') DUM,DUM,U2(:)
          ELSEIF (LEN_TRIM(CLINE2).EQ.146) THEN
            READ(CLINE2,     '(2F9.3,16F8.3)') DUM,DUM,U2(:)
          ELSE
            WRITE(6,'(2a)') 'ERROR UNKNOWN FORMAT ',TRIM(CFILE2)
            CALL EXIT(7)
          ENDIF
          WRITE(21,'(a,16f8.3)') CLINE(1:19),U1(:)-U2(:)
        ELSEIF (LEN_TRIM(CLINE).EQ.146) THEN
          READ(CLINE,      '(2F9.3,16F8.3)') DUM,DUM,U1(:)
          IF     (LEN_TRIM(CLINE2).EQ.147) THEN
            READ(CLINE2,'(F9.4,F10.4,16F8.3)') DUM,DUM,U2(:)
          ELSEIF (LEN_TRIM(CLINE2).EQ.146) THEN
            READ(CLINE2,     '(2F9.3,16F8.3)') DUM,DUM,U2(:)
          ELSE
            WRITE(6,'(2a)') 'ERROR UNKNOWN FORMAT ',TRIM(CFILE2)
            CALL EXIT(6)
          ENDIF
          WRITE(21,'(a,16f8.3)') CLINE(1:18),U1(:)-U2(:)
        ELSE
          WRITE(6,'(2a)') 'ERROR UNKNOWN FORMAT ',TRIM(CFILE1)
          CALL EXIT(8)
        ENDIF
      ENDDO !i
      CLOSE(21)
      CLOSE(11)
      CLOSE(12)
      END
