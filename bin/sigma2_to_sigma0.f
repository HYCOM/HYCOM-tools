      program sigma2_to_sigma0
      implicit none
c
c     usage:  echo dens2 saln | sigma2_to_sigma0
c
c     input:  dens2 saln
c     output: saln pot.temp dens2 dens0
c
c             dens2  = Sigma-2 potential density
c             saln   = Salinity (psu)
c
c     alan j. wallcraft, naval research laboratory, may 2005.
c
      real*8  dens0,dens2,saln,temp
      integer ios
c
      write(6,'(a)') '#     saln  pot.temp     dens2     dens0'
      do
        read(5,*,iostat=ios) dens2,saln
        if     (ios.ne.0) then
          exit
        endif
        call tofsig2(temp,  dens2,saln)
        call  sigma0(dens0, temp, saln)
        write(6,'(4f10.3)') saln,temp,dens2,dens0
      enddo
      end
      SUBROUTINE TOFSIG2(TSEA, RSEA,SSEA)
      IMPLICIT NONE
C
      REAL*8  TSEA,SSEA,RSEA
C
C     CALCULATE POTENTIAL TEMPERATURE FROM SIGMA-2 DENSITY AND SALINITY.
C
C     STATEMENT FUNCTIONS.
C
      REAL*8     DZERO,DTHIRD
      PARAMETER (DZERO=0.D0,DTHIRD=1.D0/3.D0)
C
      REAL*8 C1,C2,C3,C4,C5,C6,C7
C --- coefficients for sigma-2 (based on Brydon & Sun fit)
      PARAMETER (C1= 9.77093e+00, C2=-2.26493e-02, C3= 7.89879E-01,
     .           C4=-6.43205E-03, C5=-2.62983E-03,
     .           C6= 2.75835E-05, C7= 3.15235E-05)
C
      REAL*8  R,S,T
      REAL*8  A0,A1,A2,CUBQ,CUBR,CUBAN,CUBRL,CUBIM,TOFSIG
C
C --- auxiliary statements for finding root of 3rd degree polynomial
      A0(S)=(C1+C3*S)/C6
      A1(S)=(C2+C5*S)/C6
      A2(S)=(C4+C7*S)/C6
      CUBQ(S)=DTHIRD*A1(S)-(DTHIRD*A2(S))**2
      CUBR(R,S)=DTHIRD*(0.5D0*A1(S)*A2(S)-1.5D0*(A0(S)-R/C6))
     .   -(DTHIRD*A2(S))**3
C --- if q**3+r**2>0, water is too dense to yield real root at given
C --- salinitiy. setting q**3+r**2=0 in that case is equivalent to
C --- lowering sigma until a double real root is obtained.
      CUBAN(R,S)=DTHIRD*ATAN2(SQRT(MAX(DZERO,
     .   -(CUBQ(S)**3+CUBR(R,S)**2))),CUBR(R,S))
      CUBRL(R,S)=SQRT(-CUBQ(S))*COS(CUBAN(R,S))
      CUBIM(R,S)=SQRT(-CUBQ(S))*SIN(CUBAN(R,S))
C
C --- temp (deg c) as a function of sigma and salinity (mil)
      TOFSIG(R,S)=-CUBRL(R,S)+SQRT(3.)*CUBIM(R,S)-DTHIRD*A2(S)
C
      TSEA = TOFSIG( RSEA,SSEA )
      RETURN
      END
      SUBROUTINE SIGMA0(RSEA,TSEA,SSEA)
      IMPLICIT NONE
C
      REAL*8  TSEA,SSEA,RSEA
C
C     CONVERT FROM SIGMA-2 TO SIGMA-0.
C
C     STATEMENT FUNCTIONS.
C
      REAL*8     DZERO,DTHIRD
      PARAMETER (DZERO=0.D0,DTHIRD=1.D0/3.D0)
C
      REAL*8     C1,C2,C3,C4,C5,C6,C7
C --- coefficients for sigma-0 (based on Brydon & Sun fit)
      DATA C1,C2,C3,C4,C5,C6,C7/
     . -1.36471E-01, 4.68181E-02, 8.07004E-01,-7.45353E-03,-2.94418E-03,
     .  3.43570E-05, 3.48658E-05/
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
