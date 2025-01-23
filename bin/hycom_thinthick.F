      PROGRAM THINTHICK
      IMPLICIT NONE
C
C  hycom_thinthick - Usage:  hycom_thinthick archv.a thinQthick.a
C
C              generates the thickness ratio layer.k/max(layer.k-1,layer.k+1)
C
C   archv.a is assumed to be an HYCOM archive data file, with companion
C   header file archv.b.  Both standard and mean archive files are allowed.
C
C   thinQthick.[ab] will contain the resulting ratio fields.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  COAPS/FSU,  August 2024.
C
      REAL*4     ONEM,QONEM,SPVAL
      PARAMETER (ONEM=9806.0, QONEM=1.0/9806.0, SPVAL=2.0**100)
C
      REAL*4, ALLOCATABLE :: DP(:,:),PK(:,:,:)
      REAL*4              :: PAD(4096)
C
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      INTEGER       IDM,JDM,KDM,NSURF,NLAY,NTHK,IEXPT,YRFLAG
      INTEGER       NPAD,ITYPE,ITEST,JTEST
      REAL          THBASE,SIGMA(99),TIME
      REAL          DP_TOP,DP_CEN,DP_BOT,DP_INC
      CHARACTER*30  CMTYPE
      CHARACTER*240 CFILEA,CFILEB,CFILEM
C
      CHARACTER*18  CASN
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
      IF     (NARG.EQ.2) THEN
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CFILEM)
        ITEST = 0
        JTEST = 0
      ELSEIF (NARG.EQ.4) THEN  !undocumented, for debugging
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CFILEM)
        CALL GETARG(3,CARG)
        READ(CARG,*) ITEST
        CALL GETARG(4,CARG)
        READ(CARG,*) JTEST
      ELSE
        WRITE(6,*) 
     +    'Usage: hycom_thinthick archv.a thinQthick.a'
        CALL EXIT(1)
      ENDIF
C
C     EXTRACT MODEL PARAMETERS FROM ".b" FILE.
C
      CFILEB = CFILEA(1:LEN_TRIM(CFILEA)-1) // 'b'
      CALL READ_B(CFILEB,
     +            IEXPT,YRFLAG,IDM,JDM,KDM,NSURF,NLAY,NTHK,
     +            THBASE,SIGMA,TIME)
C
C     OPEN ".a" FILE.
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( DP(IDM,JDM),
     +          PK(IDM,JDM,KDM+1), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_thinthick: could not allocate ',
     +             8*IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
C
      IF     (NPAD.EQ.0) THEN
        INQUIRE( IOLENGTH=NRECL) DP
      ELSE
        INQUIRE( IOLENGTH=NRECL) DP,PAD(1:NPAD)
      ENDIF
*     write(6,*) 'nrecl = ',nrecl
#ifdef CRAY
#ifdef t3e
      IF     (MOD(NRECL,4096).EQ.0) THEN
        WRITE(CASN,8000) NRECL/4096
 8000   FORMAT('-F cachea:',I4.4,':1:0')
        IU8 = 11
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          WRITE(6,*) 'Error: can''t asnunit 11'
          WRITE(6,*) 'ios  = ',ios8
          WRITE(6,*) 'casn = ',casn
          CALL EXIT(5)
        ENDIF
        IU8 = 21
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          WRITE(6,*) 'Error: can''t asnunit 21'
          WRITE(6,*) 'ios  = ',ios8
          WRITE(6,*) 'casn = ',casn
          CALL EXIT(5)
        ENDIF
      ENDIF
#else
      CALL ASNUNIT(11,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t asnunit 11'
        WRITE(6,*) 'ios = ',ios
        CALL EXIT(5)
      ENDIF
      CALL ASNUNIT(21,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t asnunit 21'
        WRITE(6,*) 'ios = ',ios
        CALL EXIT(5)
      ENDIF
#endif
#endif
      OPEN(UNIT=11, FILE=CFILEA, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',CFILEA(1:LEN_TRIM(CFILEA))
        WRITE(6,*) 'ios   = ',ios
        WRITE(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
C     OPEN OUTPUT UNITS (20 AND 21).
C
      OPEN(UNIT=21, FILE=CFILEM, FORM='UNFORMATTED', STATUS='NEW',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',CFILEM(1:LEN_TRIM(CFILEM))
        WRITE(6,*) 'ios   = ',ios
        WRITE(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
      CFILEB = CFILEM(1:LEN_TRIM(CFILEM)-1) // 'b'
      OPEN(UNIT=20, FILE=CFILEB, FORM='FORMATTED', STATUS='NEW',
     +         IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',CFILEB(1:LEN_TRIM(CFILEB))
        WRITE(6,*) 'ios   = ',ios
        CALL EXIT(3)
      ENDIF
C
C     ALL LAYERS
C
      PK(:,:,1) = 0.0
      DO K= 1,KDM
        KREC0 = NSURF+NLAY*(K-1)
        CALL DAREAD(DP,IDM,JDM, KREC0+NTHK, CFILEA)
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (DP(I,J).NE.SPVAL) THEN
              PK(I,J,K+1) = PK(I,J,K) + DP(I,J)
              if     (i.eq.itest .and. j.eq.jtest) then
                write(6,'(a,i3,f14.7)') 'k,dp = ',k,dp(i,j)*qonem
                write(6,'(a,i3,f14.7)') 'k,p  = ',k,pk(i,j,k+1)*qonem
              endif
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      CLOSE(11)
C
C     OUTPUT THE RATIO
C
      DO K= 1,KDM
        ZMIN =  1.E10
        ZMAX = -1.E10
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (DP(I,J).NE.SPVAL) THEN
              if     (k.eq.1) then
                dp_cen = pk(i,j,k+1) - pk(i,j,k)   + 0.001*onem
                dp_bot = pk(i,j,k+2) - pk(i,j,k+1) + 0.001*onem
                dp_top = dp_bot
              elseif (k.eq.kdm) then
                dp_top = pk(i,j,k)   - pk(i,j,k-1) + 0.001*onem
                dp_cen = pk(i,j,k+1) - pk(i,j,k)   + 0.001*onem
                dp_bot = dp_top
              else
                dp_top = pk(i,j,k)   - pk(i,j,k-1) + 0.001*onem
                dp_cen = pk(i,j,k+1) - pk(i,j,k)   + 0.001*onem
                dp_bot = pk(i,j,k+2) - pk(i,j,k+1) + 0.001*onem
              endif
              if     (dp_cen.gt.onem) then
                DP(I,J) = MIN( 9.0, DP_CEN / MIN( DP_TOP, DP_BOT ) )
              else
                DP(I,J) = 1.0
              endif
              ZMIN = MIN( ZMIN, DP(I,J) )
              ZMAX = MAX( ZMAX, DP(I,J) )
              if     (dp(i,j).lt.0.001 .or.
     &                i.eq.itest .and. j.eq.jtest) then
                write(6,'(a,3x,2i5,i4.2,3f10.3,f10.5)')
     &                'ratio:',i,j,k,
     &                 dp_top*qonem,dp_cen*qonem,dp_bot*qonem,
     &                 DP(I,J)
              endif !debug
            ENDIF
          ENDDO
        ENDDO
        WRITE(20,'(A,F12.2,2F10.5)') 
     +    'layer ratio day,min,max =',
     +    TIME,ZMIN,ZMAX
        IF     (NPAD.EQ.0) THEN
          WRITE(21,REC=K) DP
        ELSE
          PAD(1:NPAD) = SPVAL
          WRITE(21,REC=K) DP,PAD(1:NPAD)
        ENDIF
      ENDDO
      CLOSE(20)
      CLOSE(21)
      END
      SUBROUTINE READ_B(CFILEB,
     &                  IEXPT,YRFLAG,IDM,JDM,KDM,NSURF,NLAY,NTHK,
     &                  THBASE,SIGMA,TIME)
      IMPLICIT NONE
C
      INTEGER       IEXPT,YRFLAG,IDM,JDM,KDM,NSURF,NLAY,NTHK
      REAL          THBASE,SIGMA(99),TIME
      CHARACTER*240 CFILEB
C
C     EXTRACT NEEDED MODEL PARAMETERS FROM ARCHIVE .b FILE.
C
      INTEGER      IDUM,IOS,K,KR,L,NSTEP
      LOGICAL      LSTERIC,LONETA,LWTRFLX,LSEAICE
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
      READ(12,*)  ! skip field line
C
C     FIND KDM.
C
      NSURF = 13  ! number of surface arrays
      NLAY  =  6  ! number of arrays per layer
      write(6,*) 'initial nsurf =',NSURF
C
      READ(12,'(a)') CLINE  !montg1
      READ(12,'(a)') CLINE  !srfhgt
      K = INDEX(CLINE,'=')
      READ(CLINE(K+1:),*) NSTEP,TIME
      READ(12,'(a)') CLINE  !steric or oneta or surflx
      LSTERIC = CLINE(1:8).EQ.'steric  '
      write(6,*) trim(cline)
      write(6,*) 'lsteric= ',lsteric
      IF     (LSTERIC) THEN
        NSURF = NSURF + 1
        write(6,*) 'steric  nsurf =',NSURF
        READ(12,'(a)') CLINE  !oneta or surflx 
      ENDIF
      LONETA  = CLINE(1:8).EQ.'oneta   '
      write(6,*) trim(cline)
      write(6,*) 'loneta = ',loneta
      IF     (LONETA) THEN
        NSURF = NSURF + 1
        write(6,*) 'oneta   nsurf =',NSURF
        READ(12,'(a)') CLINE  !surflx 
      ENDIF
      READ(12,'(a)') CLINE  !wtrflx or salflx
      LWTRFLX = CLINE(1:8).EQ.'wtrflx  '
      write(6,*) trim(cline)
      write(6,*) 'lwtrflx= ',lwtrflx
      IF     (LWTRFLX) THEN
        NSURF = NSURF + 1
        write(6,*) 'wtrflx  nsurf =',NSURF
        READ(12,'(a)') CLINE  !salflx
      ENDIF
         DO KR= 6,10
          READ(12,'(a)') CLINE
        ENDDO
*       write(6,*) trim(cline)
       IF     (CLINE(1:8).EQ.'surtx   ') THEN
          NSURF = NSURF + 6
        write(6,*) 'surtx   nsurf =',NSURF
          NLAY  = -1  !set below
          DO KR= 1,6
            READ(12,'(a)') CLINE
          ENDDO
*         write(6,*) trim(cline)
        ENDIF
*       write(6,*) trim(cline)
        IF     (CLINE(1:8).EQ.'thmix   ') THEN
          READ(CLINE(36:42),*) THBASE_IN
          IF     (THBASE_IN.NE.0.0) THEN
            THBASE = THBASE_IN
          ENDIF
        ELSE
          WRITE(6,*)
          WRITE(6,*) 'Expected thmix but got:'
          WRITE(6,*) CLINE(1:LEN_TRIM(CLINE))
          WRITE(6,*)
          CALL EXIT(2)
          STOP
        ENDIF
      DO K= 11,13
        READ(12,'(a)') CLINE
      ENDDO
*     write(6,*) trim(cline)
      IF     (CLINE(1:8).EQ.'kemix   ') THEN
        NLAY  = 7  ! mean archive
        NSURF = NSURF + 1
        write(6,*) 'kemix   nsurf =',NSURF
        READ(12,'(a)') CLINE
      ENDIF
*     write(6,*) trim(cline)
      IF     (CLINE(1:8).EQ.'covice  ') THEN
        NSURF = NSURF + 3
        write(6,*) 'covice  nsurf =',NSURF
        READ(12,'(a)') CLINE
        READ(12,'(a)') CLINE
        READ(12,'(a)') CLINE
      ENDIF
*     write(6,*) trim(cline)
      READ(12,'(a)') CLINE
      IF     (NLAY.EQ.7) THEN
        NSURF = NSURF + 1
        write(6,*) 'NLAY=7  nsurf =',NSURF
        READ(12,'(a)') CLINE  ! kebtrop
      ENDIF
*     write(6,*) trim(cline)
C
      L = INDEX(CLINE,'=')
      READ(CLINE(L+1:),*) NSTEP,TIME
C
C     FIND NLAY (ALLOWING FOR TRACERS) AND NTHK
C
      NTHK = -1
      READ(12,'(a)') CLINE
      DO L= 2,99
        READ(12,'(a)',IOSTAT=IOS) CLINE
        IF     (IOS.NE.0) THEN
          EXIT
        ELSEIF (CLINE(1:8).EQ.'thknss  ') THEN
          NTHK = L
          write(6,*) trim(cline)
        ELSEIF (CLINE(1:8).EQ.'u-vel.  ') THEN
*         write(6,*) trim(cline)
          EXIT
        ENDIF
      ENDDO
      NLAY = L-1
      write(6,*) 'FIND     nlay =',NLAY
      write(6,*) 'FIND     nthk =',NTHK
C
      REWIND(UNIT=12)
      DO K= 1,NSURF+10
        READ(12,'(a)') CLINE
        write(6,*) trim(cline)
      ENDDO
*     write(6,*) '-----------------------------------------------'
C
      DO K= 1,999
        READ(12,'(a)',IOSTAT=IOS) CLINE
        IF     (IOS.NE.0) THEN
          EXIT
        ELSEIF (CLINE(1:8).NE.'u-vel.  ') THEN
          write(6,*) trim(cline)
          EXIT
        ENDIF
        L = INDEX(CLINE,'=')
        READ(CLINE(L+1:),*) NSTEP,TIME,IDUM,SIGMA(K)
*       write(6,*) trim(cline),"  (1)",sigma(k)
*       write(6,*) "k,sigma",k,sigma(k)
C
        DO L= 2,NLAY
          READ(12,'(a)') CLINE
        ENDDO
      ENDDO
      KDM = K-1
*     write(6,*) 'kdm = ',kdm
      CLOSE(UNIT=12)
      RETURN
      END
      SUBROUTINE DAREAD(A,IDM,JDM, KREC, CFILEA)
      IMPLICIT NONE
C
      CHARACTER*240 CFILEA
      INTEGER       IDM,JDM,KREC
      REAL*4        A(IDM,JDM)
C
C --- READ ONE RECORD ON UNIT 11
C
      INTEGER IOS
C
      READ(11,REC=KREC,IOSTAT=IOS) A
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
