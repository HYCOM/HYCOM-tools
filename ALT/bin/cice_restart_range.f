      PROGRAM CICE_RESTART_RANGE
      IMPLICIT NONE
C
C  cice_restart_range - Usage:  cice_restart_range restart idm jdm [itest jtest ]
C
C  prints min and max of each field in a CICE restart
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  April 2015.
C
      REAL*8,  ALLOCATABLE :: A8(:,:)
C
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      INTEGER       IDM,JDM,ITEST,JTEST
      CHARACTER*240 CFILE1
      INTEGER       I,J,K,N,ISTEP,IOS
      REAL*8        RDAY,FDAY
      REAL*8        A8MAX,A8MIN
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.3) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        ITEST = 0
        JTEST = 0
      ELSEIF (NARG.EQ.5) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CARG)
        READ(CARG,*) ITEST
        CALL GETARG(5,CARG)
        READ(CARG,*) JTEST
      ELSE
        WRITE(6,*)
     &    'Usage:  cice_restart_range restart idm jdm [itest jtest]'
        CALL EXIT(1)
      ENDIF
C
      ALLOCATE( A8(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in cice_restart_range: could not allocate ',
     +             IDM*JDM,' 8-byte words'
        CALL EXIT(2)
      ENDIF
C
C --- RESTART
C
      OPEN(UNIT=21, FILE=CFILE1, FORM='UNFORMATTED', STATUS='OLD',
     +         IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILE1)
        write(6,*) 'ios   = ',ios
        CALL EXIT(3)
      ENDIF
C
C     HEADER RECORD
C
      READ( 21) ISTEP,RDAY,FDAY
      WRITE(6,'(A,I9,5X,2F20.6)')
     &  'step,time,time_forc = ',ISTEP,RDAY/86400.d0,FDAY/86400.d0
C
C     ARRAY RECORDS
C
      DO N= 1,HUGE(N)
        READ(21,IOSTAT=IOS) A8
        IF     (IOS.NE.0) THEN
          EXIT !n
        ENDIF
        A8MIN = A8(1,1)
        A8MAX = A8(1,1)
        DO J= 1,JDM
          DO I= 1,IDM
            A8MIN = MIN( A8MIN, A8(I,J) )
            A8MAX = MAX( A8MAX, A8(I,J) )
          ENDDO !j
        ENDDO !i
        WRITE(6,*) A8MIN,A8MAX
        IF     (MAX(ITEST,JTEST).GT.0) THEN
          WRITE(6,*) ITEST,JTEST,A8(ITEST,JTEST)
        ENDIF
      ENDDO !n
      CLOSE(21)
C
      END
