      PROGRAM HYCOM_PROFILE_THSTAR
      IMPLICIT NONE
C
C  hycom_profile+thstar - Usage:  hycom_profile+thstar archv.txt archs.txt kapref
C
C                 add virtual potential density (thstar) and
C                 locally referenced potential density (sigloc)
C                 to a HYCOM isopycnal text profile file
C
C   archv.txt is assumed to be an HYCOM archive text profile file
C   archs.txt will be the output text profile file, with thstar added
C   kapref is thermobaric reference state (1=Arctic; 2=Atlantic; 3=Med.)
C
C   the output profile file will not include viscosity or tracers.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  march 2002.
C
      INTEGER      IARGC
      INTEGER      NARG
      CHARACTER*240 CARG
C
      CHARACTER*240 CFILEA,CFILEC
      CHARACTER*240 CLINE
      REAL          UK,VK,TK,SK,RK,PK,THK,CENTER
      REAL          TKM,SKM,RKL,alfadt,betads
      INTEGER       IOS,K,KDM,KAPREF
C
c------------------------------------------------------------------------
      real    kappaf,sigloc,dsiglocdt,dsiglocds
      real    c1p,c2p,c3p,c4p,c5p,c6p,c7p,kappaf1
      real    r,s,t,prs
      integer kkf
c
      real       pref
      parameter (pref=2000.e4)
c
c --- coefficients for kappa^(theta)
c --- new values (w.r.t. t-toff,s-soff,prs) from Shan Sun, Sep.2004
c --- 1=Arctic/Antarctic; 2=Atlantic; 3=Mediterranean
      real, parameter, dimension(3) ::
     &  toff = (/  0.0,             3.0,            13.0 /)
     & ,soff = (/ 34.5,            35.0,            38.5 /)
     & ,qttt = (/ -3.03869354E-05, -3.03869352E-05, -3.03869353E-05 /)
     & ,qtt  = (/  4.56625601E-03,  4.29277358E-03,  3.38116552E-03 /)
     & ,qt   = (/ -2.88801209E-01, -2.61828868E-01, -1.81335007E-01 /)
     & ,qs   = (/ -1.08670290E-01, -1.05131061E-01, -9.33336309E-02 /)
     & ,qst  = (/  7.90503772E-04,  7.71096940E-04,  1.07270585E-03 /)
     & ,qpt  = (/  1.07813750E-09,  1.00638435E-09,  7.57239852E-10 /)
     & ,qpst = (/  1.41541548E-11,  1.48598578E-11,  3.89226107E-12 /)
     & ,qptt = (/ -1.31383708E-11, -1.31383707E-11, -1.31383708E-11 /)
c
c --- sub-coefficients for locally referenced sigma
c --- a fit towards Jackett & McDougall (1995)
      real, parameter, dimension(7) ::
     &  alphap = (/ -0.1364705627213484   , 0.04681812123458564,
     &               0.80700383913187     ,-0.007453530323180844,
     &              -0.002944183249153631 , 0.00003435702568990446,
     &               0.0000348657661057688 /)
     & ,betap  = (/  0.05064226654169138  ,-0.0003571087848996894,
     &              -0.0000876148051892879, 5.252431910751829e-6,
     &               1.579762259448864e-6 ,-3.466867400295792e-8,
     &              -1.687643078774232e-8 /)
     & ,gammap = (/ -5.526396144304812e-6 , 4.885838128243163e-8,
     &               9.96026931578033e-9  ,-7.251389796582352e-10,
     &              -3.987360250058777e-11, 4.006307891935698e-12,
     &               8.26367520608008e-13 /)
c --- thermobaric compressibility coefficient (integral from prs to pref)
c ---     Sun et.al. (1999) JPO 29 pp 2719-2729.
c --- kappaf1 used internally to simplify offsetting T and S,
c --- always invoke via kappaf.
c --- offset limits based on stability estimates from:
c ---     Hallberg (2005) Ocean Modelling 8 pp 279-300.
c --- t: potential temperature; s: psu; prs: pressure; kkf: ref.state
c ---     example: kappaf(4.5,34.5,1.e7,1) =  0.11411243 
c ---     example: kappaf(4.5,34.5,1.e7,2) =  0.03091669 
c ---     example: kappaf(4.5,34.5,1.e7,3) = -0.06423524
*     kappaf1(t,s,prs,kkf)=(1.e-11*qthref)*(prs-pref)*
      kappaf1(t,s,prs,kkf)=(1.e-11*1.e3)*(prs-pref)*
     &  ( s*( qs(kkf)+t* qst(kkf) ) +
     &    t*( qt(kkf)+t*(qtt(kkf)+t*qttt(kkf))+
     &        0.5*(prs+pref)*(qpt(kkf)+s*qpst(kkf)+t*qptt(kkf)) ) )
      kappaf(t,s,prs,kkf)=
     &     kappaf1(max(-1.5,         t-toff(kkf) ),
     &             max(-4.0,min(2.0, s-soff(kkf))),
     &             prs,kkf)
c
c --- locally referenced sigma, a fit towards Jackett & McDougall (1995)
c --- t: potential temperature; s: psu; prs: pressure
      c1p(prs)=alphap(1)+1.e-5*prs*(betap(1)+1.e-5*prs*gammap(1))
      c2p(prs)=alphap(2)+1.e-5*prs*(betap(2)+1.e-5*prs*gammap(2))
      c3p(prs)=alphap(3)+1.e-5*prs*(betap(3)+1.e-5*prs*gammap(3))
      c4p(prs)=alphap(4)+1.e-5*prs*(betap(4)+1.e-5*prs*gammap(4))
      c5p(prs)=alphap(5)+1.e-5*prs*(betap(5)+1.e-5*prs*gammap(5))
      c6p(prs)=alphap(6)+1.e-5*prs*(betap(6)+1.e-5*prs*gammap(6))
      c7p(prs)=alphap(7)+1.e-5*prs*(betap(7)+1.e-5*prs*gammap(7))
      sigloc(t,s,prs)=c1p(prs)+c3p(prs)*s+
     &       t*(c2p(prs)+c5p(prs)*s+t*(c4p(prs)+c7p(prs)*s+c6p(prs)*t))
      dsiglocdt(t,s,prs)=(c2p(prs)+c5p(prs)*s+
     &       2.*t*(c4p(prs)+c7p(prs)*s+1.5*c6p(prs)*t))
      dsiglocds(t,s,prs)=(c3p(prs)+t*(c5p(prs)+t*c7p(prs)))
c
c------------------------------------------------------------------------
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.3) THEN
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CFILEC)
        CALL GETARG(3,CARG)
        READ(CARG,*) KAPREF
      ELSE
        WRITE(6,*)
     +    'Usage:  hycom_profile+thstar archv.txt archs.txt kapref'
        CALL EXIT(1)
      ENDIF
C
C     OPEN ALL FILES.
C
      OPEN(UNIT=11, FILE=CFILEA, FORM='FORMATTED', STATUS='OLD',
     +     IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',CFILEA(1:LEN_TRIM(CFILEA))
        WRITE(6,*) 'ios   = ',ios
        CALL EXIT(3)
      ENDIF
      OPEN(UNIT=21, FILE=CFILEC, FORM='FORMATTED', STATUS='NEW',
     +     IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',CFILEC(1:LEN_TRIM(CFILEC))
        WRITE(6,*) 'ios   = ',ios
        CALL EXIT(5)
      ENDIF
C
C     COPY PROFILE HEADER TO OUTPUT.
C
      DO K= 1,4
        READ( 11,'(a)') CLINE
        WRITE(21,'(a)') CLINE(1:LEN_TRIM(CLINE))
      ENDDO
      READ( 11,'(a)') CLINE
      WRITE(21,'(a,a,a)')
     &  '#   k',
     &  '    utot    vtot    temp    saln    dens    thkns      dpth',
     &  '  thstar  sigloc'
C
C     READ AND WRITE THE ISOPYCNAL PROFILE.
C
      KDM  = -1
      DO K= 1,99
        READ(11,'(a)',IOSTAT=IOS) CLINE
        IF     (IOS.NE.0) THEN
          IF     (K.NE.KDM+1) THEN
            WRITE(6,*) 'Error: inconsistent input profile'
            CALL EXIT(6)
          ENDIF
          EXIT
        ENDIF
        READ( CLINE,*) KDM,UK,VK,TK,SK,RK,THK,CENTER
        PK = 9806.0*(CENTER + 0.5*THK)
        IF     (K.EQ.1) THEN
          RKL = RK
        ELSE
          alfadt=0.5*
     &             (dsiglocdt(TKM,SKM,PK)+
     &              dsiglocdt(TK, SK, PK) )*
     &             (TKM-TK)
          betads=0.5*
     &             (dsiglocds(TKM,SKM,PK)+
     &              dsiglocds(TK, SK, PK) )*
     &             (SKM-SK)
          RKL=RKL-alfadt-betads
        ENDIF
        WRITE(21,'(A,2F8.3)') CLINE(1:64),
     &                        RK+KAPPAF(TK,SK,PK,KAPREF),
     &                        RKL
        TKM = TK
        SKM = SK
      ENDDO !k
      CLOSE(11)
      CLOSE(21)
      END
