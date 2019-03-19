      program atmos_gaussian
      implicit none
c
c --- Usage: echo JH | atmos_gaussian
c ---   print the 2*JH latitudes of the specified gaussian grid
c
c --- A "Gaussian" grid, is regular in longitude and
c ---   almost regular in latitude (Hortal and Simmons, 1991). 
c --- https://en.wikipedia.org/wiki/Gaussian_grid
c
      integer              :: jh,j
      real*8               :: pscale
      real*4,  allocatable :: plat(:)
      real*8,  allocatable :: slat(:)
c
      read(5,*) jh
c
      allocate( plat(2*jh) )
      allocate( slat(  jh) )
c
      call glat8(jh,slat)
      pscale = 180.d0/acos(-1.d0)
      do j= 1,jh
        plat(2*jh+1-j) =  pscale*asin(slat(j))
        plat(       j) = -plat(2*jh+1-j)
      enddo !j
      do j= 1,2*jh
        write(6,'(f10.5)') plat(j)
      enddo
      end
      SUBROUTINE GLAT8(JH,SLAT)
      IMPLICIT NONE
C
      INTEGER JH
      REAL*8  SLAT(JH)
C
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    GLAT        COMPUTE GAUSSIAN LATITUDE FUNCTIONS
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 92-10-31
C
C ABSTRACT: COMPUTES SINES OF GAUSSIAN LATITUDE BY ITERATION.
C
C PROGRAM HISTORY LOG:
C   91-10-31  MARK IREDELL
C
C USAGE:    CALL GLAT(JH,SLAT)
C
C   INPUT ARGUMENT LIST:
C     JH       - INTEGER NUMBER OF GAUSSIAN LATITUDES IN A HEMISPHERE
C
C   OUTPUT ARGUMENT LIST:
C     SLAT     - REAL (JH) SINES OF (POSITIVE) GAUSSIAN LATITUDE
C
C ATTRIBUTES:
C   LANGUAGE: CRAY FORTRAN
C
C$$$
C
      INTEGER   JBZ
      PARAMETER(JBZ=50)
      REAL*8    PI,C,EPS
      PARAMETER(PI=3.14159265358979,C=(1.-(2./PI)**2)*0.25,EPS=1.E-14)
C
      INTEGER ITS,J,N
      REAL*8  PKM2,R,SP,SPMAX
      REAL*8  PK(JH),PKM1(JH)
C
      REAL*8  BZ(JBZ)
      DATA BZ        / 2.4048255577,  5.5200781103,
     $  8.6537279129, 11.7915344391, 14.9309177086, 18.0710639679,
     $ 21.2116366299, 24.3524715308, 27.4934791320, 30.6346064684,
     $ 33.7758202136, 36.9170983537, 40.0584257646, 43.1997917132,
     $ 46.3411883717, 49.4826098974, 52.6240518411, 55.7655107550,
     $ 58.9069839261, 62.0484691902, 65.1899648002, 68.3314693299,
     $ 71.4729816036, 74.6145006437, 77.7560256304, 80.8975558711,
     $ 84.0390907769, 87.1806298436, 90.3221726372, 93.4637187819,
     $ 96.6052679510, 99.7468198587, 102.888374254, 106.029930916,
     $ 109.171489649, 112.313050280, 115.454612653, 118.596176630,
     $ 121.737742088, 124.879308913, 128.020877005, 131.162446275,
     $ 134.304016638, 137.445588020, 140.587160352, 143.728733573,
     $ 146.870307625, 150.011882457, 153.153458019, 156.295034268 /
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  ESTIMATE LATITUDES USING BESSEL FUNCTION
      R=1./SQRT((2*JH+0.5)**2+C)
      DO J=1,MIN(JH,JBZ)
        SLAT(J)=COS(BZ(J)*R)
      ENDDO
      DO J=JBZ+1,JH
        SLAT(J)=COS((BZ(JBZ)+(J-JBZ)*PI)*R)
      ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  CONVERGE UNTIL ALL SINES OF GAUSSIAN LATITUDE ARE WITHIN EPS
      DO 110 ITS= 1,999
        SPMAX=0.
        DO J=1,JH
          PKM1(J)=1.
          PK(J)=SLAT(J)
        ENDDO
        DO N=2,2*JH
          DO J=1,JH
            PKM2=PKM1(J)
            PKM1(J)=PK(J)
            PK(J)=((2*N-1)*SLAT(J)*PKM1(J)-(N-1)*PKM2)/N
          ENDDO
        ENDDO
        DO J=1,JH
          SP=PK(J)*(1.-SLAT(J)**2)/(2*JH*(PKM1(J)-SLAT(J)*PK(J)))
          SLAT(J)=SLAT(J)-SP
          SPMAX=MAX(SPMAX,ABS(SP))
        ENDDO
        IF     (SPMAX.LE.EPS) THEN
          GOTO 1110
        ENDIF
  110 CONTINUE
 1110 CONTINUE
      RETURN
C     END OF GLAT8.
      END
