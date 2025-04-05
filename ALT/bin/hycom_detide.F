      PROGRAM HYCOM_DETIDE
      IMPLICIT NONE
C
C  hycom_detide - Usage:  hycom_detide fhr.a idm jdm fout.a
C
C                 Outputs a 49-hour filtered hourly time series
C
C  The filter is the convolution of a 21 hr (10 point) 2nd-order
C   Savitzky-Golay smoothing filter and a 24.842 hr boxcar filter.
C   it passes 0.02% of semi-diurnal and 3.2% of diurnal tides
C   (1.2% of the total tides).
C
C  The first 48 output fields are zero as the filter initializes.
C
C  fhr.a is assumed to be a HYCOM data file containing a sequence of
C   the same field one hour apart, with companion text file fhr.b.
C   It contains idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  COAPS/FSU,  November 2024.
C
      REAL*4, ALLOCATABLE :: A49(:,:,:),B(:,:)
      REAL*4              :: PAD(4096)
      INTEGER IOS
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      INTEGER       IDM,JDM,ITEST,JTEST,NPAD
      CHARACTER*240 CFILE,CFILEO
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.4) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        ITEST = 0
        JTEST = 0
        CALL GETARG(4,CFILEO)
      ELSEIF (NARG.EQ.6) THEN  !undocumented, for debugging
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CARG)
        READ(CARG,*) ITEST
        CALL GETARG(5,CARG)
        READ(CARG,*) JTEST
        CALL GETARG(6,CFILEO)
      ELSE
        WRITE(6,*)
     &  'Usage: hycom_detide fhr.a idm jdm fout.a'
        CALL EXIT(1)
      ENDIF
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( A49(IDM,JDM,49), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_detide: could not allocate 1st ',
     +             49*IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( B(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_detide: could not allocate 2nd ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
C
      CALL DETIDE(A49,B,IDM,JDM,PAD,NPAD,
     &            ITEST,JTEST, CFILE,CFILEO)
      CALL EXIT(0)
      END
      SUBROUTINE DETIDE(A49,B,IDM,JDM,PAD,NPAD,
     &                  ITEST,JTEST, CFILE,CFILEO)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE,CFILEO
      INTEGER      IDM,JDM,NPAD,ITEST,JTEST
      REAL*4       A49(IDM,JDM,49),
     +             B(IDM,JDM),PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
      CHARACTER*240 CFILEB,CLINE
      CHARACTER*18  CASN
      INTEGER       I,IHR,J,K,IOS,NRECL,N,NH
      REAL*4        AMN,AMX
      REAL*8        DAY,DAYI
C
      real*4  fhrly(49)
      real*4  f2hrly(0:24)
      real*4  fc
      real*8  fsum 
C
      data    f2hrly / 9.297873 ,
     &                 9.297873 , 9.338932 , 9.835878     ,
     &                 10.04647 , 10.00111 , 9.730179     ,
     &                 9.264084 , 8.633219 , 7.867976     ,
     &                 6.998751 , 6.055939 , 5.069937     ,
     &                 4.071137 , 3.089936 , 2.156729     ,
     &                 1.301912 , 0.5558784,-5.0975680E-02,
     &                -0.4882553,-0.7255653,-0.7325106    ,
     &                -0.4786960, 0.0      , 0.0           /
!
! --- normalize the filter weights
      fsum = f2hrly(0)
      do ihr=1,24
        fsum = fsum + 2.0*f2hrly(ihr)
      enddo
      fc = 1.0/fsum
      fhrly(25) = fc*f2hrly(0)
      do ihr=1,24
        fhrly(25+ihr) = fc*f2hrly(ihr)
        fhrly(25-ihr) = fc*f2hrly(ihr)
      enddo
C
      INQUIRE( IOLENGTH=NRECL) B,PAD
#ifdef CRAY
#ifdef t3e
      IF     (MOD(NRECL,4096).EQ.0) THEN
        WRITE(CASN,8000) NRECL/4096
 8000   FORMAT('-F cachea:',I4.4,':1:0')
        IU8 = 11
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          write(6,*) 'Error: can''t asnunit 11'
          write(6,*) 'ios  = ',ios8
          write(6,*) 'casn = ',casn
          CALL EXIT(5)
        ENDIF
        IU8 = 21
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          write(6,*) 'Error: can''t asnunit 21'
          write(6,*) 'ios  = ',ios8
          write(6,*) 'casn = ',casn
          CALL EXIT(5)
        ENDIF
      ENDIF
#else
      CALL ASNUNIT(11,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 11'
        write(6,*) 'ios = ',ios
        CALL EXIT(5)
      ENDIF
      CALL ASNUNIT(21,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 21'
        write(6,*) 'ios = ',ios
        CALL EXIT(5)
      ENDIF
#endif
#endif
      OPEN(UNIT=11, FILE=CFILE, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',CFILE(1:LEN_TRIM(CFILE))
        write(6,*) 'unit  = ',11
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
      OPEN(UNIT=21, FILE=CFILEO, FORM='UNFORMATTED', STATUS='NEW',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',CFILEO(1:LEN_TRIM(CFILEO))
        write(6,*) 'unit  = ',21
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
C --- .b header
C
      CFILEB = CFILE(1:LEN_TRIM(CFILE)-1) // 'b'
      OPEN(UNIT=12, FILE=CFILEB, FORM='FORMATTED', STATUS='OLD',
     +         IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',CFILEB(1:LEN_TRIM(CFILEB))
        write(6,*) 'unit  = ',12
        WRITE(6,*) 'ios   = ',ios
        CALL EXIT(3)
      ENDIF
C
      CFILEB = CFILEO(1:LEN_TRIM(CFILEO)-1) // 'b'
      OPEN(UNIT=22, FILE=CFILEB, FORM='FORMATTED', STATUS='NEW',
     +         IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',CFILEB(1:LEN_TRIM(CFILEB))
        write(6,*) 'unit  = ',22
        WRITE(6,*) 'ios   = ',ios
        CALL EXIT(3)
      ENDIF
      READ( 12,'(a)')      CLINE  !1st
      WRITE(22,'(a)') TRIM(CLINE)
      READ( 12,'(a)')      CLINE  !2nd
      WRITE(22,'(a)') TRIM(CLINE)
      READ( 12,'(a)')      CLINE  !3rd
      WRITE(22,'(a)') TRIM(CLINE)
      READ( 12,'(a)')      CLINE  !4th
      CLINE = 'hycom_detide ' // TRIM(CLINE)
      WRITE(22,'(a)') TRIM(CLINE)
      READ( 12,'(a)')      CLINE  !5th
      WRITE(22,'(a)') TRIM(CLINE)
C
C --- TIME SERIES PROCESSING
C
      AMN = 0.0
      AMX = 0.0
        B(:,:)   = 0.0
      A49(:,:,:) = 0.0
      DO K= 1,9999999
C ---   PHYSICALLY SHIFT THE ARRAYS (NOT MOST EFFICIENT APPROACH)
        DO J= 1,JDM
          DO I= 1,IDM
            DO N= 1,48
              A49(I,J,N) = A49(I,J,N+1)
            ENDDO !n
          ENDDO !i
        ENDDO !j
        READ(11,REC=K,IOSTAT=IOS) A49(:,:,49)
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A49(1,1,49),IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          EXIT
        ENDIF
        READ( 12,'(a)') CLINE
C
        IF     (K.GE.49) THEN
          AMN =  SPVAL
          AMX = -SPVAL
          DO J= 1,JDM
            DO I= 1,IDM
              IF     (A49(I,J,1) .NE. SPVAL) THEN
                FSUM = 0.0
                DO N=1,49
                  FSUM = FSUM + FHRLY(N)*A49(I,J,N)
                ENDDO !n
                B(I,J) = FSUM
                AMX = MAX( AMX, B(I,J) )
                AMN = MIN( AMN, B(I,J) )
              ELSE
                B(I,J) = SPVAL
              ENDIF
            ENDDO !i
          ENDDO !j
#ifdef ENDIAN_IO
          CALL ENDIAN_SWAP(B,IDM*JDM)
#endif
        ENDIF !>=49
C
        WRITE(21,REC=K,IOSTAT=IOS) B
        N = INDEX(CLINE,'=')
        READ(CLINE(N+1:),*) DAY,DAYI
        WRITE(22,'(a,1x,2F11.6,1p2g16.8)')
     &     CLINE(1:N),DAY,DAYI,AMN,AMX
        IF     (MIN(ITEST,JTEST).GT.0) THEN
          WRITE(6,'(a,F12.6,F10.6,1pg16.8)')
     &       CLINE(1:N),DAY,DAYI,B(ITEST,JTEST)
        ELSE
          WRITE(6,'(a,F12.6,F10.6,1p2g16.8)')
     &       CLINE(1:N),DAY,DAYI,AMN,AMX
        ENDIF
      ENDDO !k
C
      CLOSE(11)
      CLOSE(21)
C
      RETURN
      END
