      PROGRAM HYCOM_TIDAL_FILTER
      IMPLICIT NONE
C
C  hycom_tidal_filter - Usage:  hycom_tidal_filter fhr.a idm jdm bwM2 bwS2 bwK1 bwO1  fout.a
C
C                 Outputs a streaming tidal filtered hourly time series
C
C  The filter detects instantaneous tidal signals in the hourly input using
C   up to four coupled ODEs (the filter equations), given the target 
C   frequencies and the bandwidth parameters (bwXX) of the filters.
C   For each hourly timestep, the filter returns sinusoidal motions near
C   its target frequency.
C
C  Xu, C., & Zaron, E. D. (2024). Detecting instantaneous tidal signals in
C   ocean models utilizing streaming band-pass filters. JAMES.
C
C  fhr.a is assumed to be a HYCOM data file containing a sequence of
C   the same field one hour apart, with companion text file fhr.b.
C   It contains idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  0 < bwXX < 1 is a parameter that determines the bandwidth, with smaller
C   values blocking more of the off-frequency signal.  bwM2 (bwK1) of
C   0.1 passes 0.25% of the diurnal (semi‐diurnal) signal.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  COAPS/FSU,  November 2024.
C
      REAL*4, ALLOCATABLE :: A(:,:),AF(:,:),S(:,:,:),B(:,:,:)
      REAL*4              :: PAD(4096)
      INTEGER IOS
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      REAL*4        BW(4),FT(4)
      REAL*4        tidal_frequency
      INTEGER       IDM,JDM,ITEST,JTEST,NPAD
      CHARACTER*240 CFILE,CFILEO
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.8) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        ITEST = 0
        JTEST = 0
        CALL GETARG(4,CARG)
        READ(CARG,*) BW(1)
        CALL GETARG(5,CARG)
        READ(CARG,*) BW(2)
        CALL GETARG(6,CARG)
        READ(CARG,*) BW(3)
        CALL GETARG(7,CARG)
        READ(CARG,*) BW(4)
        CALL GETARG(8,CFILEO)
      ELSEIF (NARG.EQ.10) THEN  !undocumented, for debugging
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CARG)
        READ(CARG,*) ITEST
        CALL GETARG(5,CARG)
        READ(CARG,*) JTEST
        CALL GETARG(6,CARG)
        READ(CARG,*) BW(1)
        CALL GETARG(7,CARG)
        READ(CARG,*) BW(2)
        CALL GETARG(8,CARG)
        READ(CARG,*) BW(3)
        CALL GETARG(9,CARG)
        READ(CARG,*) BW(4)
        CALL GETARG(10,CFILEO)
      ELSE
        WRITE(6,'(2a)')
     &  'Usage: hycom_tidal_filter ',
     &          'fhr.a idm jdm bwM2 bwS2 bwK1 bwO1 fout.a'
        CALL EXIT(1)
      ENDIF
C
      FT(1) = tidal_frequency('M2')
      FT(2) = tidal_frequency('S2')
      FT(3) = tidal_frequency('K1')
      FT(4) = tidal_frequency('Q1')
      WRITE(6,*) '# FT = ',FT(:)
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE(  A(IDM,JDM), 
     &          AF(IDM,JDM), 
     &           S(IDM,JDM,4), 
     &           B(IDM,JDM,4), 
     &          STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_tidal_filter: could not allocate ',
     &             10*IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
C
      CALL FILTER(A,AF,S,B,IDM,JDM,PAD,NPAD,
     &            ITEST,JTEST, FT,BW, CFILE,CFILEO)
      CALL EXIT(0)
      END
      SUBROUTINE FILTER(A,AF,S,B,IDM,JDM,PAD,NPAD,
     &                  ITEST,JTEST, FT,BW, CFILE,CFILEO)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE,CFILEO
      INTEGER      IDM,JDM,NPAD,ITEST,JTEST
      REAL*4       A(IDM,JDM),AF(IDM,JDM)
      REAL*4       S(IDM,JDM,4),B(IDM,JDM,4),PAD(NPAD)
      REAL*4       FT(4),BW(4)
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
      CHARACTER*240 CFILEB,CLINE
      CHARACTER*18  CASN
      INTEGER       I,J,L,KR,IOS,N,NRECL
      REAL*4        AMN,AMX,C1,C2
      REAL*8        DAY,DAYI
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
      CLINE = 'hycom_tidal_filter ' // TRIM(CLINE)
      WRITE(22,'(a)') TRIM(CLINE)
      READ( 12,'(a)')      CLINE  !5th
      WRITE(22,'(a)') TRIM(CLINE)
C
C --- TIME SERIES PROCESSING
C
      DO KR= 1,9999999
        READ(11,REC=KR,IOSTAT=IOS) A
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          EXIT
        ENDIF
        READ( 12,'(a)') CLINE
C
        IF     (KR.EQ.1) THEN
           S(:,:,:) = 0.0
           DO L= 1,4
             IF     (BW(L).EQ.0.0) THEN
               B(:,:,L) = 0.0
             ELSE
               B(:,:,L) = A(:,:)
             ENDIF
           ENDDO !l
        ELSE
          DO L= 1,4
            IF     (BW(L).NE.0.0) THEN
              C1 = FT(L) * 300.0  !5 minute time step, 12 steps per hour
              C2 = 1.0 - BW(L) * C1
!!!!!         IF (KR.LE.3) THEN
!!!!!           WRITE(6,*) '# C1,C2,FT,BW = ',C1,C2,FT(L),BW(L)
!!!!!         ENDIF
              DO J= 1,JDM
                DO I= 1,IDM
                  IF     (A(I,J) .NE. SPVAL) THEN
                    DO N= 1,12
                      S(I,J,L) =      S(I,J,L) +  C1 *  B(I,J,L)
                      B(I,J,L) = C2 * B(I,J,L) + -C1 * (S(I,J,L) -
     &                                         BW(L) *  A(I,J))
!!!!!                 IF (KR.LE.3 .AND. I.EQ.ITEST .AND. J.EQ.JTEST) THEN
!!!!!                   WRITE(6,*) '# S,B,A =',L,S(I,J,L),B(I,J,L),A(I,J)
!!!!!                 ENDIF
                    ENDDO
                  ENDIF
                ENDDO !i
              ENDDO !j
            ENDIF !bw
          ENDDO !l
        ENDIF !kr
C
        AMN =  SPVAL
        AMX = -SPVAL
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (A(I,J) .NE. SPVAL) THEN
              AF(I,J) = B(I,J,1) + B(I,J,2) + B(I,J,3) + B(I,J,4)
!!!!!           IF (KR.LE.3 .AND. I.EQ.ITEST .AND. J.EQ.JTEST) THEN
!!!!!             WRITE(6,*) '# AF,A =',AF(I,J),A(I,J)
!!!!!           ENDIF
              AMX = MAX( AMX, AF(I,J) )
              AMN = MIN( AMN, AF(I,J) )
            ELSE
              AF(I,J) = SPVAL
            ENDIF
          ENDDO !i
        ENDDO !j
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(AF,IDM*JDM)
#endif
C
        WRITE(21,REC=KR,IOSTAT=IOS) AF
        N = INDEX(CLINE,'=')
        READ(CLINE(N+1:),*) DAY,DAYI
        WRITE(22,'(a,F13.6,F10.6,1p2g16.8)')
     &     CLINE(1:N),DAY,DAYI,AMN,AMX
        IF     (MAX(ITEST,JTEST).EQ.0) THEN
          WRITE(6,'(a,F13.6,F10.6,1p2g16.8)')
     &       CLINE(1:N),DAY,DAYI,AMN,AMX
        ELSE
          WRITE(6,'(a,F13.6,F10.6,1pg16.8)')
     &       CLINE(1:N),DAY,DAYI,AF(ITEST,JTEST)
        ENDIF
      ENDDO !k
C
      CLOSE(11)
      CLOSE(12)
      CLOSE(21)
      CLOSE(22)
C
      RETURN
      END
      function tidal_frequency(constit)
      real*4            :: tidal_frequency !> Angular frequency [rad s-1]
      character (len=2) :: constit         !> Constituent to look up
c
      select case (constit)
        case ("M2","m2")
          tidal_frequency = 1.4051890e-4
        case ("S2","s2")
          tidal_frequency = 1.4544410e-4
        case ("N2","n2")
          tidal_frequency = 1.3787970e-4
        case ("K2","k2")
          tidal_frequency = 1.4584234e-4 
        case ("K1","k1") 
          tidal_frequency = 0.7292117e-4
        case ("O1","o1") 
          tidal_frequency = 0.6759774e-4
        case ("P1","p1") 
          tidal_frequency = 0.7252295e-4
        case ("Q1","q1")
          tidal_frequency = 0.6495854e-4
        case default
          write(6,*) 'tidal_frequency: unrecognized constituent ',
     &               constit
      end select
      end function tidal_frequency
