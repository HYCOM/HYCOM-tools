      PROGRAM PROFILE
      IMPLICIT NONE
C
C  hycom_profile_argo - Usage:  hycom_profile_argo argoT.txt argoS.txt archz.txt
C
C                 converts an Argo T&S profile to a hycom profile file.
C
C   argoT.txt is assumed to be a Argo temperature text profile file
C   argoS.txt is assumed to be a Argo salinity    text profile file
C   archz.txt will be the output hycom text profile file
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  Oct. 2008
C
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      CHARACTER*240 CFILEA,CFILEB,CFILEC
      CHARACTER*240 CHEAD(20),CLINE
      REAL          ALAT,ALON,THK,BOT,PZK,PZKP,RZK,ZZK
      INTEGER       IOS,ITYPE,K,KZ
      INTEGER       I,KMAX
C
      REAL, ALLOCATABLE :: ZZ(:),TZ(:),SZ(:)
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.3) THEN
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CFILEB)
        CALL GETARG(3,CFILEC)
      ELSE
        WRITE(6,*)
     +  'Usage: hycom_profile_argo argoT.txt argoS.txt archz.txt'
        CALL EXIT(1)
      ENDIF
C
C     OPEN ALL FILES.
C
      OPEN(UNIT=11, FILE=CFILEA, FORM='FORMATTED', STATUS='OLD',
     +     IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',TRIM(CFILEA)
        WRITE(6,*) 'ios   = ',ios
        CALL EXIT(3)
      ENDIF
      OPEN(UNIT=12, FILE=CFILEB, FORM='FORMATTED', STATUS='OLD',
     +     IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',TRIM(CFILEB)
        WRITE(6,*) 'ios   = ',ios
        CALL EXIT(4)
      ENDIF
      OPEN(UNIT=21, FILE=CFILEC, FORM='FORMATTED', STATUS='NEW',
     +     IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',TRIM(CFILEC)
        WRITE(6,*) 'ios   = ',ios
        CALL EXIT(5)
      ENDIF
C
C     READ THE ARGO PROFILE HEADER, TO GET KZ.
C
      DO K= 1,20
        READ(11,'(a)',IOSTAT=IOS) CHEAD(K)
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'Error: bad Argo T file'
          CALL EXIT(6)
        ENDIF
        READ(12,'(a)',IOSTAT=IOS) CLINE
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'Error: bad Argo S file'
          CALL EXIT(6)
        ENDIF
        IF     (K.NE.18 .AND. CHEAD(K).NE.CLINE) THEN
          WRITE(6,*) 'Error: bad Argo files (line',k,')'
          WRITE(6,*) TRIM(CHEAD(K))
          WRITE(6,*) TRIM(CLINE)
          CALL EXIT(6)
        ENDIF
      ENDDO !k
      IF     (CHEAD( 3)(1:32).NE.
     &        "%  profile latitude            :") THEN
        WRITE(6,*) 'Error: bad Argo file (line  3)'
        WRITE(6,*) TRIM(CHEAD( 3))
        CALL EXIT(6)
      ENDIF
      READ(CHEAD( 3)(33:),*) ALAT
      IF     (CHEAD( 4)(1:32).NE.
     &        "%  profile longitude           :") THEN
        WRITE(6,*) 'Error: bad Argo file (line  4)'
        WRITE(6,*) TRIM(CHEAD( 4))
        CALL EXIT(6)
      ENDIF
      READ(CHEAD( 4)(33:),*) ALON
      IF     (CHEAD( 5)(1:32).NE.
     &        "%  profile observed DTG        :") THEN
        WRITE(6,*) 'Error: bad Argo file (line  5)'
        WRITE(6,*) TRIM(CHEAD( 5))
        CALL EXIT(6)
      ENDIF
      IF     (CHEAD( 7)(1:32).NE.
     &        "%  DBDBV bottom depth          :") THEN
        WRITE(6,*) 'Error: bad Argo file (line  7)'
        WRITE(6,*) TRIM(CHEAD( 7))
        CALL EXIT(6)
      ENDIF
      READ(CHEAD( 7)(33:),*) BOT
      IF     (CHEAD(11)(1:32).NE.
     &        "%  observed temperature levels :") THEN
        WRITE(6,*) 'Error: bad Argo file (line 11)'
        WRITE(6,*) TRIM(CHEAD(11))
        CALL EXIT(6)
      ENDIF
      READ(CHEAD(11)(33:),*) KZ
      IF     (CHEAD(12)(1:32).NE.
     &        "%  observed salinity levels    :") THEN
        WRITE(6,*) 'Error: bad Argo file (line 12)'
        WRITE(6,*) TRIM(CHEAD(12))
        CALL EXIT(6)
      ENDIF
      READ(CHEAD(12)(33:),*) K
      IF     (K.NE.KZ) THEN
        WRITE(6,*) 'Error: Argo does not have both T&S'
        WRITE(6,*) TRIM(CHEAD(11))
        WRITE(6,*) TRIM(CHEAD(12))
        CALL EXIT(6)
      ENDIF
C
C     INPUT ALL LEVELS
C
      ALLOCATE( ZZ(KZ), TZ(KZ), SZ(KZ) )
C
      DO K= 1,KZ
        READ(11,*,IOSTAT=IOS) ZZ(K),TZ(K)
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'Error: bad Argo T file (level',K,')'
          CALL EXIT(6)
        ENDIF
        READ(12,*,IOSTAT=IOS) ZZK,SZ(K)
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'Error: bad Argo S file (level',K,')'
          CALL EXIT(6)
        ENDIF
        IF     (ZZ(K).NE.ZZK) THEN
          WRITE(6,*) 'Error: Argo does not have both T&S'
          WRITE(6,*) 'K,Z,T = ',K,ZZ(K),TZ(K)
          WRITE(6,*) 'K,Z,S = ',K,ZZK,  SZ(K)
          CALL EXIT(6)
        ENDIF
      ENDDO !k
      CLOSE(11)
      CLOSE(12)
C
C     OUTPUT.
C
      CALL SIG_I(2)  !7-term sigma-2
      CALL SIG_P(TZ(1),SZ(1),RZK)
      WRITE(21,'(2a/a,4i7,2f7.1,i7,a)')
     &  '##   expt    idm    jdm    kdm    lon    lat',
     &  ' yrflag   DTG',
     &  '##',    0,     1,     1,    KZ,  ALON,  ALAT,
     &  3,CHEAD(5)(35:46)
      WRITE(21,'(3a/a,f11.2,f8.2,f8.1,2f9.3,3f8.3,4f8.2)')
     &  '## model-day  srfhgt  surflx',
     &  '     dpbl   dpmixl    tmix    smix   thmix',
     &  '    umix    vmix   ubavg   vbavg',
     &  '#',
     &  0.0,   ! model day
     &  0.0,   ! cm
     &  0.0,   ! W/m**2
     &  0.0,   ! m
     &  0.0,   ! m
     &  TZ(1), ! degC
     &  SZ(1), ! psu
     &  RZK,   ! SigmaT
     &  0.0,   ! cm/s
     &  0.0,   ! cm/s
     &  0.0,   ! cm/s
     &  0.0    ! cm/s
      WRITE(21,'(a,a)')
     &  '#   k',
     &  '    utot    vtot    temp    saln    dens    thkns      dpth'
      PZKP = 0.0
      DO K= 1,KZ
        PZK  = PZKP
        IF     (K.LT.KZ) THEN
          PZKP = 0.5*(ZZ(K)+ZZ(K+1))
        ELSE
          PZKP = MAX(ZZ(KZ),BOT)
        ENDIF
        CALL POT_T(TZ(K),SZ(K),ZZ(K),ALAT)
        CALL SIG_P(TZ(K),SZ(K),RZK)
        WRITE(21,'(i4,1x,2f8.2,3f8.3,f9.3,f10.3)')
     &    K,
     &    0.0,               !cm/s
     &    0.0,               !cm/s
     &    TZ(K),             !degC
     &    SZ(K),             !psu
     &    RZK,               !SigmaT
     &    PZKP-PZK,          !m
     &    ZZ(K)              !m
      ENDDO
      CLOSE(21)
      END
      SUBROUTINE POT_T(T,S,Z,ALAT)
      IMPLICIT NONE
C
      REAL T,S,Z,ALAT
C
C     CONVERT TEMPERATURE TO POTENTIAL TEMPERATURE
C
      REAL*8 Z8,T8,S8,P8,PT8,L8
C
      Z8 = Z
      T8 = T
      S8 = S
      P8 = 0.0
      L8 = ALAT
      CALL PTEM85(Z8,T8,S8,P8,PT8,L8)
      T = PT8
*     WRITE(6,'(a,5f10.3)') 'pot_t - ',PT8,T8,S8,Z8,L8
      END
        subroutine ptem85(zobs,tobs,sobs,pref,ptem,alat)
        implicit real*8 (a-h,o-z)
c
c       *********************************************
c
c       this subroutine calculates potential temperature on an arbi-
c       trary reference pressure level, based on a fourth-order
c       runge-kutta integration of the adiabatic lapse rate.
c       numerical errors are 0.1 mdeg c or less for a single step
c       integration.  overall accuracy is determined by errors of
c       adiabatic lapse rate estimates.  see fofonoff (1977) dsr.,
c       vol. 24, pp. 489-491 for runge-kutta scheme and overal
c       error estimates.  note: see dsr vol. 25 for correction to
c       vol. 24 scheme.  in addition see bryden (1973) dsr. vol. 20,
c       pp. 401-408 for the empirical formulation of adiabatic
c       lapse rate.
c
c       note:  the classical definition of potential temperature is
c       that temperature a water parcel would have if raised
c       adiabatically and isentropically to the ocean surface.
c       potential temperature may be computed here by using a
c        reference pressure (pref) = 0. dbars.
c
c         variable definitions:
c            zobs = observed depth(m)
c            tobs = observed temperature (degrees c)
c            sobs = observed salinity (ppt)
c            pref = reference pressure level (db)
c            ptem = output variable, potential temperaure (deg. c)
c              on the reference pressure level.
c            alat = latitude of the observation (deg)
c
c            written october 1979, steven worley, tamu
c
c      december 1985 - this subroutine was renamed to ptem85 from
c                      ptem77 - when it was up dated with a new
c                      depth to pressure algorithm based on fofonoff's
c                      latest equation of state.
c
c      ************************************************************
cc      common/conre1/ioffp,spval1
c
c       bryden's empirical function for adiabatic lapse rate
c       gamma (deg c /1000 db)
c
      gamma(p,t,s)=0.35803e-1 + p*(0.18741e-4 - 0.67795e-6*t +
     &    0.87330e-8*t*t - 0.54481e-10*t*t*t + p*(-0.46206e-9 +
     &    0.18676e-10*t - 0.21687e-12*t*t)) + (s-35.)*(0.18932e-2 -
     &    0.42393e-4*t - 0.11351e-6*p + 0.27759e-8*t*p) +
     &    0.85258e-2*t - 0.68360e-4*t*t + 0.66228e-6*t*t*t
c
      data spval/-.03125/
      data spval1/-.03125/
c
      if (sobs.eq.spval1.or.tobs.eq.spval1) then
      ptem=spval1
      return
      endif
      s=sobs
      p=pres85(zobs,alat)
      t=tobs
      pr=pref
      delp=pr-p
      delt1=delp*gamma(p,t,s)/1000.
      t1=t+.5*delt1
      p2=p+.5*delp
      delt2=delp*gamma(p2,t1,s)/1000.
      q1=delt1
      sqrt2=sqrt(2.)
      t2=t1+(1.-1./sqrt2)*(delt2-q1)
      delt3=delp*gamma(p2,t2,s)/1000.
      q2=(2.-sqrt2)*delt2+(-2.+3./sqrt2)*q1
      t3=t2+(1.+1./sqrt2)*(delt3-q2)
      p3=delp+p
      delt4=delp*gamma(p3,t3,s)/1000.
      q3=(2.+sqrt2)*delt3+(-2-3/sqrt2)*q2
      t4=t3+(1./6.)*(delt4-2.*q3)
      ptem=t4
      return
      end
      real*8 function pres85(z,lat)
      implicit real*8 (a-h,o-z)
c
c     a function obtained by inverting the empirical depth
c     function given by fofonoff in 1985
c
c     the accuracy of this approximate algorthim is .1 db or less
c     at all depths less than 7500 meters
c     see testing information in file for specifices
c     -- worley december 1985 ------
c
c     input depth in meters and latitude in degrees
c
      real*8 lat
      data spval/-.03125/
      data spval1/-.03125/
c
      if(z.eq.spval)go to 99
c
      x=sin(lat/57.29578)
      x=x*x
c
c     the gravity anomaly using z instead of p for an approximation
c
      gr=9.780318*(1.0+(5.2788e-3+2.36e-5*x)*x) + 1.092e-6*z
      dp=gr*z
c
c     scale dp to account for the scaling done in 4th order linear
c     model used to obtain the coefficients
c
      dp=dp/1.0e3
c
      pres85=(((-2.641561e-7*dp+2.524876e-5)*dp+2.267149e-2)*dp +
     *       1.028361e2)*dp
c
      return
99    pres85=spval
      return
      end
