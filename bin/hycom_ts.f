      PROGRAM TSCURVE
      IMPLICIT NONE
C
C  hycom_ts - Usage:  hycom_ts density
C
C                 generates a text file table of T and S at density
C                 for HYCOM's equation of state (sigma0 version).
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  August 2001.
C
      INTEGER      IARGC
      INTEGER      NARG
      CHARACTER*240 CARG
C
      INTEGER      ITT
      REAL*8       RR,TT,SS
c
c --- coefficients for sigma-0 (based on Brydon & Sun fit)
      real*8     c1,c2,c3,c4,c5,c6,c7
      parameter (c1=-1.36471E-01, c2= 4.68181E-02, c3= 8.07004E-01,
     .           c4=-7.45353E-03, c5=-2.94418E-03,
     .           c6= 3.43570E-05, c7= 3.48658E-05)
c
c --- coefficients for sigma-2 (based on Brydon & Sun fit)
csig2 real       c1,c2,c3,c4,c5,c6,c7
csig2 parameter (c1= 9.77093E+00, c2=-2.26493E-02, c3= 7.89879E-01,
csig2.           c4=-6.43205E-03, c5=-2.62983E-03,
csig2.           c6= 2.75835E-05, c7= 3.15235E-05)
c
c --- salinity (mil) as a function of sigma and temperature (deg c)
      real*8 sofsig,r,t
      sofsig(r,t)=(r-c1-t*(c2+t*(c4+c6*t)))/(c3+t*(c5+c7*t))
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.1) THEN
        CALL GETARG(1,CARG)
        READ(CARG,*) RR
      ELSE
        WRITE(6,*) 'Usage: hycom_ts density'
        CALL EXIT(1)
      ENDIF
C
      WRITE(6,'(a,f7.3 / a)')
     &  '#   HYCOM equation of state, potential density =',RR,
     &  '#      T          S        RHO'
C
      DO ITT= -50,350
        TT = ITT*0.1
        SS = SOFSIG(RR,TT)
        WRITE(6,'(3f11.6)') TT,SS,RR
      ENDDO
      END
