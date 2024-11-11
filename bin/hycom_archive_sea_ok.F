      PROGRAM ASEAOK
      IMPLICIT NONE
C
C  hycom_archive_sea_ok - Usage:  hycom_archive_sea_ok archv.a depth.a anom.a
C
C                 checks that the archive is consistent with the depth file
C
C   archv.a is assumed to be an HYCOM archive data file, with companion
C   header file archv.b.  Both standard and mean archive files are allowed.
C
C   anom.a will contain difference between the sum of the layers and depth
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  COAPS/FSU, January 2024.
C
      REAL*4     QONEM,SPVAL
      PARAMETER (QONEM=1.0/9806.0, SPVAL=2.0**100)
C
      REAL*4, ALLOCATABLE :: DEPTH(:,:),DP(:,:),PK(:,:)
      REAL*4              :: PAD(4096)
C
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      INTEGER       IDM,JDM,KDM,NSURF,NLAY,IEXPT,YRFLAG
      INTEGER       IBADS,NBAD1,NBAD2
      INTEGER       NPAD,ITYPE,ITEST,JTEST
      REAL          THBASE,SIGMA(99),TIME
      CHARACTER*240 CFILEA,CFILEB,CFILED,CFILEM
C
      INTEGER       I,J,K,KREC,KREC0,IOS,NRECL
      REAL          ZMIN,ZMAX
#ifdef CRAY
      INTEGER*8     IU8,IOS8
#endif
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.3) THEN
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CFILED)
        CALL GETARG(3,CFILEM)
        ITEST = 0
        JTEST = 0
      ELSEIF (NARG.EQ.5) THEN  !undocumented, for debugging
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CFILED)
        CALL GETARG(3,CFILEM)
        CALL GETARG(4,CARG)
        READ(CARG,*) ITEST
        CALL GETARG(5,CARG)
        READ(CARG,*) JTEST
      ELSE
        WRITE(6,*) 
     +    'Usage: hycom_archive_sea_ok archv.a depth.a anom.a'
        CALL EXIT(1)
      ENDIF
C
C     EXTRACT MODEL PARAMETERS FROM ".b" FILE.
C
      CFILEB = CFILEA(1:LEN_TRIM(CFILEA)-1) // 'b'
      CALL READ_B(CFILEB,
     +            IEXPT,YRFLAG,IDM,JDM,KDM,NSURF,NLAY,
     +            THBASE,SIGMA,TIME)
C
C     OPEN ".a" FILE.
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( DEPTH(IDM,JDM),
     +             DP(IDM,JDM),
     +             PK(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_archive_sea_ok: could not allocate ',
     +             3*IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
C
      IF     (NPAD.EQ.0) THEN
        INQUIRE( IOLENGTH=NRECL) DP
      ELSE
        INQUIRE( IOLENGTH=NRECL) DP,PAD(1:NPAD)
      ENDIF
*     write(6,*) 'nrecl = ',nrecl
C
      OPEN(UNIT=11, FILE=CFILEA, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',CFILEA(1:LEN_TRIM(CFILEA))
        WRITE(6,*) 'ios   = ',ios
        WRITE(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
      OPEN(UNIT=51, FILE=CFILED, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',CFILED(1:LEN_TRIM(CFILEA))
        WRITE(6,*) 'ios   = ',ios
        WRITE(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
C     OPEN OUTPUT UNITS
C
      OPEN(UNIT=21, FILE=CFILEM, FORM='UNFORMATTED', STATUS='NEW',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',CFILEM(1:LEN_TRIM(CFILEM))
        WRITE(6,*) 'ios   = ',ios
        WRITE(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
C     DEPTH
C
      CALL DAREAD(DEPTH,IDM,JDM, 1, 51,CFILED)
      CLOSE(51)
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (DEPTH(I,J).LT.0.0) THEN
            DEPTH(I,J) = SPVAL
          ENDIF
        ENDDO
      ENDDO
C
C     ALL LAYERS
C
      PK(:,:) = 0.0
      DO K= 1,KDM
        KREC0 = NSURF+NLAY*(K-1)
        CALL DAREAD(DP,IDM,JDM, KREC0+1, 11,CFILEA)
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (DP(I,J).NE.SPVAL) THEN
              PK(I,J) = PK(I,J) + DP(I,J)
              if     (i.eq.itest .and. j.eq.jtest) then
                write(6,'(a,i3,f14.7)') 'k,p  = ',k,pk(i,j)*qonem
              endif
            ELSE
              PK(I,J) = SPVAL
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      CLOSE(11)
C
      IBADS = 0
      NBAD1 = 0
      NBAD2 = 0
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (DEPTH(I,J).EQ.SPVAL .AND. PK(I,J).NE.SPVAL) THEN
            NBAD1 = NBAD1 + 1
            IF     (NBAD1.LE.100 .OR. MOD(NBAD1,100).EQ.1) THEN
              WRITE(6,"(a,2i6)") 'sea over land:',I,J
            ENDIF
          ENDIF
          IF     (DEPTH(I,J).NE.SPVAL .AND. PK(I,J).EQ.SPVAL) THEN
            NBAD2 = NBAD2 + 1
            IF     (NBAD2.LE.100 .OR. MOD(NBAD2,100).EQ.1) THEN
              WRITE(6,"(a,2i6)") 'land over sea:',I,J
            ENDIF
          ENDIF
C
          IF     (PK(I,J).NE.SPVAL) THEN
            DP(I,J) = DEPTH(I,J) - PK(I,J)*QONEM
            ZMIN = MIN( ZMIN, DP(I,J) )
            ZMAX = MAX( ZMAX, DP(I,J) )
            if     (abs(dp(i,j)).gt.0.01) then  !> 1cm
              ibads = ibads + 1
              IF     (IBADS.LE.100 .OR. MOD(IBADS,100).EQ.1) THEN
                write(6,'(a,2i6,2f14.8)')
     &            'large archive depth mismatch: ',
     &            i,j,dp(i,j),depth(i,j)
              ENDIF
            endif !large mismatch
          ELSE
            DP(I,J) = SPVAL
          ENDIF
        ENDDO
      ENDDO
      write(6,*)
      IF     (NBAD1+NBAD2.EQ.0) THEN
        WRITE(6,'(A)') 'ARCHIVE LAND/SEA is OK'
      ENDIF
      IF     (NBAD1.NE.0) THEN
        WRITE(6,'(A,A,I9,A)') 'ARCHIVE',
     &                           ' has',NBAD1,' sea values over land'
      ENDIF
      IF     (NBAD2.NE.0) THEN
        WRITE(6,'(A,A,I9,A)') 'ARCHIVE',
     &                           ' has',NBAD2,' land values over sea'
      ENDIF
      if     (ibads.ne.0) then
        write(6,*)
        write(6,*) 'error - wrong bathymetry for this archive file'
        write(6,*) 'number of depth mismatches = ',ibads
      endif !ibads.ne.0
      write(6,*)
C
C     OUTPUT THE DEPTH ANOMALY
C
      WRITE(6,'(A,2F14.8)') 
     +  'Depth anomaly: min,max =',
     +  ZMIN,ZMAX
      IF     (NPAD.EQ.0) THEN
        WRITE(21,REC=1) DP
      ELSE
        PAD(1:NPAD) = SPVAL
        WRITE(21,REC=1) DP,PAD(1:NPAD)
      ENDIF
      CLOSE(21)
      END
      SUBROUTINE READ_B(CFILEB,
     &                  IEXPT,YRFLAG,IDM,JDM,KDM,NSURF,NLAY,
     &                  THBASE,SIGMA,TIME)
      IMPLICIT NONE
C
      INTEGER       IEXPT,YRFLAG,IDM,JDM,KDM,NSURF,NLAY
      REAL          THBASE,SIGMA(99),TIME
      CHARACTER*240 CFILEB
C
C     EXTRACT NEEDED MODEL PARAMETERS FROM ARCHIVE .b FILE.
C
      INTEGER      IDUM,IOS,K,L,NSTEP
      REAL         THBASE_IN
      CHARACTER*6  CVARIN*6
      CHARACTER*240 CLINE
C
      OPEN(UNIT=12, FILE=CFILEB, FORM='FORMATTED', STATUS='OLD',
     +         IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',CFILEB(1:LEN_TRIM(CFILEB))
        WRITE(6,*) 'ios   = ',ios
        CALL EXIT(3)
      ENDIF
      READ(12,*)  ! skip title(1)
      READ(12,*)  ! skip title(2)
      READ(12,*)  ! skip title(3)
      READ(12,*)  ! skip title(4)
      READ(12,*)  ! skip iversn
      READ(12,*) IEXPT,CVARIN
      IF     (CVARIN.NE.'iexpt ') THEN
        WRITE(6,*) 'Error in hycom_profile: bad .b file'
        WRITE(6,*) 'filename: ',CFILEB(1:LEN_TRIM(CFILEB))
        CALL EXIT(4)
      endif
      READ(12,*) YRFLAG
      READ(12,*) IDM
      READ(12,*) JDM
C
C     FIND NSURF
C
      READ(12,'(a)') CLINE
      DO L= 1,99
        READ(12,'(a)',IOSTAT=IOS) CLINE
        IF     (IOS.NE.0) THEN
          EXIT
        ELSEIF (CLINE(1:8).EQ.'thknss  ') THEN
*         write(6,*) trim(cline)
          EXIT
        ENDIF
      ENDDO
      NSURF = L-1
C
C     FIND NLAY (ALLOWING FOR TRACERS)
C
      DO L= 2,99
        READ(12,'(a)',IOSTAT=IOS) CLINE
        IF     (IOS.NE.0) THEN
          EXIT
        ELSEIF (CLINE(1:8).EQ.'thknss  ') THEN
*         write(6,*) trim(cline)
          EXIT
        ENDIF
      ENDDO
      NLAY = L-1
C
C     FIND KDM.
C
      REWIND(UNIT=12)
      DO K= 1,NSURF+10
        READ(12,'(a)') CLINE
*       write(6,*) trim(cline)
      ENDDO
*     write(6,*) '-----------------------------------------------'
C
      DO K= 1,999
        READ(12,'(a)',IOSTAT=IOS) CLINE
        IF     (IOS.NE.0) THEN
          EXIT
        ELSEIF (CLINE(1:8).NE.'thknss  ') THEN
*         write(6,*) trim(cline)
          EXIT
        ENDIF
C
        DO L= 2,NLAY
          READ(12,'(a)',IOSTAT=IOS) CLINE
          IF     (IOS.NE.0) THEN
            EXIT
          ENDIF
        ENDDO
      ENDDO
      KDM = K-1
*     write(6,*) 'kdm = ',kdm
      CLOSE(UNIT=12)
      RETURN
      END
      SUBROUTINE DAREAD(A,IDM,JDM, KREC, IUNIT, CFILEA)
      IMPLICIT NONE
C
      CHARACTER*240 CFILEA
      INTEGER       IDM,JDM,KREC,IUNIT
      REAL*4        A(IDM,JDM)
C
C --- READ ONE RECORD ON UNIT IUNIT
C
      INTEGER IOS
C
      READ(IUNIT,REC=KREC,IOSTAT=IOS) A
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read record ',KREC,
     +             ' from file ',TRIM(CFILEA)
        CALL EXIT(4)
        STOP
      ENDIF
      END
