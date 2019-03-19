      PROGRAM PROFILE
      IMPLICIT NONE
C
C  hycom_profile+sig - Usage:  hycom_profile+sig archv.txt archs.txt [type]
C
C                 add potential density to a HYCOM text profile file
C
C   archv.txt is assumed to be an HYCOM archive text profile file
C   archs.txt will be the output text profile file, with sig[02] added
C   type      is 0 for sig0 (the default) and 2 for sig2
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
      REAL*8        T8,S8,R8
      INTEGER       IOS,ISIG,K,KDM
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.2) THEN
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CFILEC)
        ISIG = 0  ! default sig0
      ELSEIF (NARG.EQ.3) THEN
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CFILEC)
        CALL GETARG(3,CLINE)
        READ(CLINE,*) ISIG
      ELSE
        WRITE(6,*)
     +    'Usage:  hycom_profile+sig archv.txt archs.txt [type]'
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
      IF     (ISIG.EQ.0) THEN
        WRITE(21,'(a,a,a)')
     &    '#  k',
     &    '    utot    vtot    temp    saln    dens    thkns      dpth',
     &    '     sig0'
      ELSE
        WRITE(21,'(a,a,a)')
     &    '#  k',
     &    '    utot    vtot    temp    saln    dens    thkns      dpth',
     &    '     sig2'
      ENDIF
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
        T8 = TK
        S8 = SK
        IF     (ISIG.EQ.0) THEN
          CALL SIGMA0(R8,T8,S8)
        ELSE
          CALL SIGMA2(R8,T8,S8)
        ENDIF
        WRITE(21,'(A,F8.3)') CLINE(1:64),R8
      ENDDO
      CLOSE(11)
      CLOSE(21)
      END
      SUBROUTINE SIGMA0(RSEA,TSEA,SSEA)
      IMPLICIT NONE
C
      REAL*8  TSEA,SSEA,RSEA
C
C     CONVERT FROM SIGMA-0 TO SIGMA-2.
C
C     STATEMENT FUNCTIONS.
C
      REAL*8     DZERO,DTHIRD
      PARAMETER (DZERO=0.D0,DTHIRD=1.D0/3.D0)
C
      REAL*8     C1,C2,C3,C4,C5,C6,C7
C --- coefficients for sigma-0 (based on Brydon & Sun fit)
      PARAMETER (C1=-1.36471e-01, C2= 4.68181e-02, C3= 8.07004E-01,
     &           C4=-7.45353e-03, C5=-2.94418E-03,
     &           C6= 3.43570e-05, C7= 3.48658e-05)
C
      REAL*8  S,T
      REAL*8  SIG
C
C --- sigma-theta as a function of temp (deg c) and salinity (mil)
C --- (friedrich-levitus 3rd degree polynomial fit)
      SIG(T,S)=(C1+C3*S+T*(C2+C5*S+T*(C4+C7*S+C6*T)))
C
C     T AND S TO SIG-2,
C
      RSEA = SIG( TSEA, SSEA )
      RETURN
      END
      SUBROUTINE SIGMA2(RSEA,TSEA,SSEA)
      IMPLICIT NONE
C
      REAL*8  TSEA,SSEA,RSEA
C
C     CONVERT FROM SIGMA-0 TO SIGMA-2.
C
C     STATEMENT FUNCTIONS.
C
      REAL*8     DZERO,DTHIRD
      PARAMETER (DZERO=0.D0,DTHIRD=1.D0/3.D0)
C
      REAL*8     C1,C2,C3,C4,C5,C6,C7
C --- coefficients for sigma-2 (based on Brydon & Sun fit)
      PARAMETER (C1= 9.77093e+00, C2=-2.26493e-02, C3= 7.89879E-01,
     &           C4=-6.43205E-03, C5=-2.62983E-03,
     &           C6= 2.75835E-05, C7= 3.15235E-05)
C
      REAL*8  S,T
      REAL*8  SIG
C
C --- sigma-theta as a function of temp (deg c) and salinity (mil)
C --- (friedrich-levitus 3rd degree polynomial fit)
      SIG(T,S)=(C1+C3*S+T*(C2+C5*S+T*(C4+C7*S+C6*T)))
C
C     T AND S TO SIG-2,
C
      RSEA = SIG( TSEA, SSEA )
      RETURN
      END
