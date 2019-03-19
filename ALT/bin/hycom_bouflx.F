      PROGRAM BOUFLX
      IMPLICIT NONE
C
C  hycom_bouflx - Usage:  hycom_bouflx archv.a bouflx.a [grid.a]
C
C                 generates surface bouancy flux from a HYCOM archive file.
C
C   archv.a is assumed to be an HYCOM archive data file, with companion
C   header file archv.b.  Both standard and mean archive files are allowed.
C
C   bouflx.a and bouflx.b will contain the resulting bouancy field and
C   the corresponding Golystin depth.
C
C   grid.a is a hycom grid file, default regional.grid.a.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  April 2003.
C
      REAL*4     QONEM,SPVAL
      PARAMETER (QONEM=1.0/9806.0, SPVAL=2.0**100)
C
      REAL*4, ALLOCATABLE :: ST(:,:),SS(:,:),GD(:,:),F(:,:),
     +                       HF(:,:),SF(:,:),BF(:,:)
      REAL*4              :: PAD(4096)
C
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      INTEGER       IDM,JDM,KDM,NSURF,NLAY,IEXPT,YRFLAG
      INTEGER       NPAD,ITYPE,ITEST,JTEST
      REAL          THBASE,TIME
      CHARACTER*240 CFILEA,CFILEB,CFILEM,CFILEG
C
      CHARACTER*18  CASN
      INTEGER       I,J,IOS,NRECL
      REAL          BMIN,BMAX,GMIN,GMAX
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
        CFILEG = 'regional.grid.a'
        ITEST = 0
        JTEST = 0
      ELSEIF (NARG.EQ.3) THEN
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CFILEM)
        CALL GETARG(3,CFILEG)
        ITEST = 0
        JTEST = 0
      ELSEIF (NARG.EQ.4) THEN  !undocumented, for debugging
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CFILEM)
        CALL GETARG(3,CARG)
        READ(CARG,*) ITEST
        CALL GETARG(4,CARG)
        READ(CARG,*) JTEST
        CFILEG = 'regional.grid.a'
      ELSEIF (NARG.EQ.5) THEN  !undocumented, for debugging
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CFILEM)
        CALL GETARG(3,CFILEG)
        CALL GETARG(4,CARG)
        READ(CARG,*) ITEST
        CALL GETARG(5,CARG)
        READ(CARG,*) JTEST
      ELSE
        WRITE(6,*) 
     +    'Usage: hycom_bouflx archv.a bouflx.a'
        CALL EXIT(1)
      ENDIF
C
C     EXTRACT MODEL PARAMETERS FROM ".b" FILE.
C
      CFILEB = CFILEA(1:LEN_TRIM(CFILEA)-1) // 'b'
      CALL READ_B(CFILEB,
     +            IEXPT,YRFLAG,IDM,JDM,KDM,NSURF,NLAY, THBASE,TIME)
C
C     OPEN ".a" FILE.
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( ST(IDM,JDM), SS(IDM,JDM), GD(IDM,JDM),
     +           F(IDM,JDM),
     +          HF(IDM,JDM), SF(IDM,JDM), BF(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_bouflx: could not allocate ',
     +             7*IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
C
      IF     (NPAD.EQ.0) THEN
        INQUIRE( IOLENGTH=NRECL) ST
      ELSE
        INQUIRE( IOLENGTH=NRECL) ST,PAD(1:NPAD)
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
C     ONLY NEED FIRST LAYER.
C
      CALL DAREAD(HF, IDM,JDM,            3, CFILEA)  !heat flux
      CALL DAREAD(SF, IDM,JDM,            4, CFILEA)  !salinity flux
      CALL DAREAD(ST, IDM,JDM, NSURF+NLAY-2, CFILEA)  !SST
      CALL DAREAD(SS, IDM,JDM, NSURF+NLAY-1, CFILEA)  !SSS
      CLOSE(11)
      if     (itest.ne.0) then
        write(6,'(a,g20.7)')  'sst  = ',st(itest,jtest)
        write(6,'(a,g20.7)')  'sss  = ',ss(itest,jtest)
        write(6,'(a,g20.7)')  'hflx = ',hf(itest,jtest)
        write(6,'(a,g20.7)')  'sflx = ',sf(itest,jtest)
        call flush(6)
      endif
C
C     BOUANCY FLUX.
C
      IF     (THBASE.LT.30.0) THEN
        CALL BOUFLX0(BF, HF,SF, ST,SS, IDM,JDM) !sigma0
      ELSE
        CALL BOUFLX2(BF, HF,SF, ST,SS, IDM,JDM) !sigma2
      ENDIF
      if     (itest.ne.0) then
        write(6,'(a,g20.7)')  'bflx = ',bf(itest,jtest)
        call flush(6)
      endif
C
C     Golystin scale, sqrt(buoflx/f**3)
C
      OPEN(UNIT=11, FILE=CFILEG, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',CFILEG(1:LEN_TRIM(CFILEG))
        WRITE(6,*) 'ios   = ',ios
        WRITE(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
      CALL DAREAD(F, IDM,JDM, 18, CFILEG)  !Coriolis
      CLOSE(11)
      if     (itest.ne.0) then
        write(6,'(a,g20.7)')  'f    = ',f(itest,jtest)
        call flush(6)
      endif
C
      BMIN =  HUGE(BMIN)
      BMAX = -HUGE(BMAX)
      GMIN =  HUGE(GMIN)
      GMAX = -HUGE(GMAX)
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (ST(I,J).NE.SPVAL) THEN
            GD(I,J) = SQRT(ABS(BF(I,J))/
     +                     MAX(ABS(F(I,J)),1.E-5)**3)  !f(4N) is about 1.e-5
            BMIN = MIN( BMIN, BF(I,J) )
            BMAX = MAX( BMAX, BF(I,J) )
            GMIN = MIN( GMIN, GD(I,J) )
            GMAX = MAX( GMAX, GD(I,J) )
          ELSE
            GD(I,J) = SPVAL
          ENDIF
        ENDDO
      ENDDO
      if     (itest.ne.0) then
        write(6,'(a,g20.7)')  'goly = ',gd(itest,jtest)
        call flush(6)
      endif
C
C     OUTPUT THE BUOYANCY FLUX AND THE GOLYSTIN DEPTH.
C
      WRITE(20,'(A,F12.2,2G20.7)') 
     +    'bouflx: day,min,max =',
     +    TIME,BMIN,BMAX
      WRITE(20,'(A,F12.2,2G20.7)') 
     +    'l(rot): day,min,max =',
     +    TIME,GMIN,GMAX
      IF     (NPAD.EQ.0) THEN
        WRITE(21,REC=1) BF
        WRITE(21,REC=2) GD
      ELSE
        PAD(1:NPAD) = SPVAL
        WRITE(21,REC=1) BF,PAD(1:NPAD)
        WRITE(21,REC=2) GD,PAD(1:NPAD)
      ENDIF
      CLOSE(20)
      CLOSE(21)
      END
      SUBROUTINE READ_B(CFILEB,
     &                  IEXPT,YRFLAG,IDM,JDM,KDM,NSURF,NLAY,
     &                  THBASE,TIME)
      IMPLICIT NONE
C
      INTEGER      IEXPT,YRFLAG,IDM,JDM,KDM,NSURF,NLAY
      REAL         THBASE,TIME
      CHARACTER*240 CFILEB
C
C     EXTRACT NEEDED MODEL PARAMETERS FROM ARCHIVE .b FILE.
C
      INTEGER      IOS,K,NSTEP
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
C     FIND KDM.
C
      NSURF = 13  ! number of surface arrays
      NLAY  =  6  ! number of arrays per layer
C
      DO K= 1,10
        READ(12,'(a)') CLINE
      ENDDO
*     write(6,*) cline(1:len_trim(cline))
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
      ENDIF
      DO K= 11,13
        READ(12,'(a)') CLINE
      ENDDO
*     write(6,*) cline(1:len_trim(cline))
      IF     (CLINE(1:8).EQ.'kemix   ') THEN
        NLAY  = 7  ! mean archive
        NSURF = NSURF + 1
        READ(12,'(a)') CLINE
      ENDIF
*     write(6,*) cline(1:len_trim(cline))
      IF     (CLINE(1:8).EQ.'covice  ') THEN
        NSURF = NSURF + 3
        READ(12,'(a)') CLINE
        READ(12,'(a)') CLINE
        READ(12,'(a)') CLINE
      ENDIF
*     write(6,*) cline(1:len_trim(cline))
      READ(12,'(a)') CLINE
      IF     (NLAY.EQ.7) THEN
        NSURF = NSURF + 1
        READ(12,'(a)') CLINE  ! kebtrop
      ENDIF
*     write(6,*) cline(1:len_trim(cline))
C
      K = INDEX(CLINE,'=')
      READ(CLINE(K+1:),*) NSTEP,TIME
C
      DO K= 1,999
        READ(12,'(a)',IOSTAT=IOS) CLINE
        IF     (IOS.NE.0) THEN
          EXIT
        ELSEIF (CLINE(1:8).NE.'u-vel.  ') THEN
*         write(6,*) cline(1:len_trim(cline))
          EXIT
        ENDIF
*       write(6,*) cline(1:len_trim(cline)),"  (1)"
C
        READ(12,'(a)') CLINE
        READ(12,'(a)') CLINE
        READ(12,'(a)') CLINE
        READ(12,'(a)') CLINE
        READ(12,'(a)') CLINE
        IF     (NLAY.EQ.7) THEN
          READ(12,'(a)') CLINE
        ENDIF
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
      SUBROUTINE BOUFLX0(BF, HF,SF, ST,SS, IDM,JDM)
      IMPLICIT NONE
C
      INTEGER       IDM,JDM
      REAL*4        BF(IDM,JDM),
     +              HF(IDM,JDM),SF(IDM,JDM),
     +              ST(IDM,JDM),SS(IDM,JDM)
C
C --- CONVERT HEAT AND SALINTY FLUX INTO A BOUANCY FLUX.
C
      REAL*4     SPVAL
      REAL*4     G,THREF,SPCIFH
      PARAMETER (SPVAL=2.0**100)
      PARAMETER (G=9.806, THREF=1.0e-3, SPCIFH=3990.0)
C
      INTEGER I,J
C
c-----------------------------------------------------------------------------
      real*4 dsigdt,dsigds
      real*4 s,t
c
c --- coefficients for sigma-0 (based on Brydon & Sun fit)
      real*4     c1,c2,c3,c4,c5,c6,c7
      parameter (c1=-1.36471E-01, c2= 4.68181E-02, c3= 8.07004E-01,
     &           c4=-7.45353E-03, c5=-2.94418E-03,
     &           c6= 3.43570E-05, c7= 3.48658E-05)
c
c --- sigma-theta as a function of temp (deg c) and salinity (mil)
c --- (friedrich-levitus 3rd degree polynomial fit)
c
c     sig(t,s)=(c1+c3*s+t*(c2+c5*s+t*(c4+c7*s+c6*t)))
c
c --- d(sig)/dt
      dsigdt(t,s)=(c2+c5*s+2.*t*(c4+c7*s+1.5*c6*t))
c
c --- d(sig)/ds
      dsigds(t,s)=(c3+t*(c5+t*c7))
c-----------------------------------------------------------------------------
C
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (ST(I,J).NE.SPVAL) THEN
            BF(I,J) = G*THREF*(
     &                 DSIGDT(ST(I,J),SS(I,J))*HF(I,J)*THREF/SPCIFH +
     &                 DSIGDS(ST(I,J),SS(I,J))*SF(I,J)*THREF         )
          ELSE
            BF(I,J) = 0.0
          ENDIF
        ENDDO
      ENDDO
      END
      SUBROUTINE BOUFLX2(BF, HF,SF, ST,SS, IDM,JDM)
      IMPLICIT NONE
C
      INTEGER       IDM,JDM
      REAL*4        BF(IDM,JDM),
     +              HF(IDM,JDM),SF(IDM,JDM),
     +              ST(IDM,JDM),SS(IDM,JDM)
C
C --- CONVERT HEAT AND SALINTY FLUX INTO A BOUANCY FLUX.
C
      REAL*4     SPVAL
      REAL*4     G,THREF,SPCIFH
      PARAMETER (SPVAL=2.0**100)
      PARAMETER (G=9.806, THREF=1.0e-3, SPCIFH=3990.0)
C
      INTEGER I,J
C
c-----------------------------------------------------------------------------
      real*4 dsigdt,dsigds
      real*4 s,t
c
c --- coefficients for sigma-2 (based on Brydon & Sun fit)
      real*4     c1,c2,c3,c4,c5,c6,c7
      parameter (c1= 9.77093E+00, c2=-2.26493E-02, c3= 7.89879E-01,
     &           c4=-6.43205E-03, c5=-2.62983E-03,
     &           c6= 2.75835E-05, c7= 3.15235E-05)
c
c --- sigma-theta as a function of temp (deg c) and salinity (mil)
c --- (friedrich-levitus 3rd degree polynomial fit)
c
c     sig(t,s)=(c1+c3*s+t*(c2+c5*s+t*(c4+c7*s+c6*t)))
c
c --- d(sig)/dt
      dsigdt(t,s)=(c2+c5*s+2.*t*(c4+c7*s+1.5*c6*t))
c
c --- d(sig)/ds
      dsigds(t,s)=(c3+t*(c5+t*c7))
c-----------------------------------------------------------------------------
C
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (ST(I,J).NE.SPVAL) THEN
            BF(I,J) = G*THREF*(
     &                 DSIGDT(ST(I,J),SS(I,J))*HF(I,J)*THREF/SPCIFH +
     &                 DSIGDS(ST(I,J),SS(I,J))*SF(I,J)*THREF         )
          ELSE
            BF(I,J) = 0.0
          ENDIF
        ENDDO
      ENDDO
      END
