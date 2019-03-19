      SUBROUTINE AVERMS(X,IW,JW, XAVE,XRMS)
      IMPLICIT NONE
C
      INTEGER IW,JW
      REAL*4  X(IW,JW),XAVE,XRMS
C
C     ROUTINE TO CALCULATE THE MEAN AND RMS OF X.
C     ***NOTE*** THIS CALCULATES STATISTICS OVER THE ENTIRE ARRAY, 
C     THUS IF THE INPUT WIND DATA SET HAS NON-REALISTIC VALUES OVER 
C     LAND THE NUMBERS MAY BE MEANINGLESS, BUT THEY CAN STILL BE 
C     USED AS A CHECK BETWEEN DIFFERENT MACHINES.
C
      INTEGER I,J
      REAL*8  XSUM,XSUMSQ
C
      XSUM   = 0.D0
      XSUMSQ = 0.D0
      DO 110 J=1,JW
        DO 111 I=1,IW
           XSUM   = XSUM   + X(I,J)
           XSUMSQ = XSUMSQ + X(I,J)**2
  111   CONTINUE
  110 CONTINUE
C
      XAVE =       XSUM   / (IW*JW)
      XRMS = SQRT( XSUMSQ / (IW*JW) ) 
C
      RETURN
C     END OF AVERMS.
      END
      SUBROUTINE MINMAX(X,IW,JW,XMIN,XMAX)
      IMPLICIT NONE
C
      INTEGER IW,JW
      REAL*4  X(IW,JW), XMIN,XMAX
C
C     FIND THE MINUMUM AND MAXIMUM OF X.
C     ***NOTE*** THIS CALCULATES STATISTICS OVER THE ENTIRE ARRAY, 
C     THUS IF THE INPUT WIND DATA SET HAS NON-REALISTIC VALUES OVER 
C     LAND THE NUMBERS MAY BE MEANINGLESS, BUT THEY CAN STILL BE 
C     USED AS A CHECK BETWEEN DIFFERENT MACHINES.
C
      INTEGER I,J
C
      XMIN =  1.0E10
      XMAX = -1.0E10
      DO 110 J=1,JW
        DO 111 I=1,IW
          XMIN = MIN(XMIN,X(I,J))
          XMAX = MAX(XMAX,X(I,J))
  111   CONTINUE
  110 CONTINUE
      RETURN
C     END OF MINMAX.
      END
      SUBROUTINE GAUSS(YAF,PLAT,IW,JW, YAG,JWI,YFIN,WK)
      IMPLICIT NONE
C
      INTEGER IW,JW,JWI
      REAL*4  YAF(IW,JW),PLAT(IW,JW), YAG(JWI),YFIN,WK(JWI)
C
C     MAP FROM GAUSSIAN TO MODEL GRIDS.
C     ON ENTRY: PLAT CONTAIN MODEL LATITUDES
C     ON EXIT:  YAF  CONTAIN CORRESPONDING GAUSSIAN GRID COORDS
C
      INTEGER I,J
      REAL*4  PSCALE
C
      CALL GLAT(JWI/2,WK, YAG(1),YAG(JWI/2+1))
      PSCALE = 180.0/ACOS(-1.0)
      DO 110 J= 1,JWI/2
        YAG(JWI+1-J) =  PSCALE*ASIN(WK(J))
        YAG(      J) = -YAG(JWI+1-J)
*       WRITE(6,*) 'J,YAG = ',J,YAG(J),YAG(JWI+1-J)
  110 CONTINUE
C
      IF     (MOD(JWI,2).NE.0) THEN
        WRITE(6,*)
        WRITE(6,*) 'ERROR IN GAUSS - JWI MUST BE EVEN'
        WRITE(6,*)
        WRITE(6,*) 'JWI = ',JWI
        WRITE(6,*)
        STOP
      ELSEIF (ABS(YAG(1)-YFIN).GT.0.1) THEN
        WRITE(6,*)
        WRITE(6,*) 'ERROR IN GAUSS - YFIN NOT AT START OF GRID'
        WRITE(6,*)
        WRITE(6,*) 'YFIN,YAG = ',YFIN,YAG(1),YFIN-YAG(1)
        WRITE(6,*)
        STOP
      ENDIF
C
      DO I= 1,IW
        DO J= 1,JW
          CALL GAUSSI(YAF(I,J), PLAT(I,J),YAG,JWI)
        ENDDO
      ENDDO
      RETURN
C     END OF GAUSS.
      END
      SUBROUTINE GAUSSI(YA, YO,YAG,JWI)
      IMPLICIT NONE
C
      INTEGER JWI
      REAL*4  YA,YO, YAG(JWI)
C
C     MAP ONE POINT FROM MODEL TO GAUSSIAN GRID.
C     NOTE THE ADDITIONAL +2.0 OFFSET FOR THE EXTENDED GRID.
C
      INTEGER J
C
      IF     (YAG(1)  .GE.YO) THEN
        YA = 3.0          ! south pole maps to first grid point
      ELSEIF (YAG(JWI).LE.YO) THEN
        YA = 1.999 + JWI  ! north pole maps to last  grid point
      ELSE
        DO 110 J= 1,JWI-1
          IF     (YAG(J).LE.YO .AND. YO.LT.YAG(J+1)) THEN
            YA = 2.0 + J + (YO-YAG(J))/(YAG(J+1)-YAG(J))
            GOTO 1110
          ENDIF
  110   CONTINUE
          WRITE(6,*)
          WRITE(6,*) 'ERROR - MODEL GRID POINT NOT IN GAUSSIAN GRID.'
          WRITE(6,*) 'MODEL LATITUDE = ',YO
          WRITE(6,*) 'GAUSSIAN GRID  = ',YAG(1),' TO ',YAG(JWI)
          WRITE(6,*)
          STOP
 1110   CONTINUE
      ENDIF
      RETURN
C     END OF GAUSSI.
      END
      SUBROUTINE GLAT(JH,SLAT,PK,PKM1)
      IMPLICIT NONE
C
      INTEGER JH
      REAL*4  SLAT(JH),PK(JH),PKM1(JH)
C
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    GLAT        COMPUTE GAUSSIAN LATITUDE FUNCTIONS
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 92-10-31
C
C ABSTRACT: COMPUTES SINES OF GAUSSIAN LATITUDE BY ITERATION.
C           THE GAUSSIAN WEIGHTS ARE ALSO COMPUTED.
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
      REAL*4    PI,C,EPS
*     PARAMETER(PI=3.14159265358979,C=(1.-(2./PI)**2)*0.25,EPS=1.E-14)
      PARAMETER(PI=3.14159265358979,C=(1.-(2./PI)**2)*0.25,EPS=1.E-7)
C
      INTEGER ITS,J,N
      REAL*4  PKM2,R,SP,SPMAX
C
      REAL*4  BZ(JBZ)
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
      DO 110 ITS= 1,99
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
C     END OF GLAT.
      END
      SUBROUTINE PATCH(FLD,FX,FY,NDX,NX,NY,
     +                 FLDI,NDXI,NXI,NYI, NPATCH, FLAND)
      IMPLICIT NONE
C
      INTEGER NDX,NX,NY, NDXI,NXI,NYI,NPATCH
      REAL*4  FLD(NDX,NY),FX(NDX,NY),FY(NDX,NY)
      REAL*4  FLDI(NDXI,NYI),FLAND
C
C**********
C*
C  1) INTERPOLATE FROM THE ARRAY FLDI TO THE ARRAY FLD, WHERE 
C      FLD(I,J) IS AT (FX(I,J),FY(I,J)) W.R.T. THE FLDI GRID
C      (1:NXI,1:NYI).
C
C     AVERAGE OVER THE NEAREST (2*NPATCH+1)x(2*NPATCH+1) GRID
C      POINT PATCH.
C
C     FOR SIMPLICITY, IT IS ASSUMED THAT FX LIES BETWEEN 3 AND NXI-2 
C      AND THAT FY LIES BETWEEN 3 AND NYI-2.  THIS AVOIDS SPECIAL 
C      CASES DURING THE INTERPOLATION, BUT MAY REQUIRE THE CONSTRUCTION 
C      OF A LARGER FLDI ARRAY FROM THE ORIGINAL BEFORE CALLING BESSEL.
C      IF SUCH A LARGER RECTANGLE IS REQUIRED IT SHOULD BE FILLED WITH
C      THE APPROPRIATE VALUES FROM THE ARRAY IN PERIODIC CASES, AND
C      OTHERWISE WITH VALUES LINEARLY EXTRAPOLATED FROM THE ARRAY.
C
C  2) ARGUMENT LIST:
C       FLD     - INTERPOLATED ARRAY ON EXIT
C       FX      - MAPPING OF 1ST DIMENSION OF FLD INTO THAT OF FLDI
C                  FX(I,J).GE.3 .AND. FX(I,J).LE.NXI-2
C       FY      - MAPPING OF 2ND DIMENSION OF FLD INTO THAT OF FLDI
C                  FY(I,J).GE.3 .AND. FY(I,J).LE.NYI-2
C       NDX     - ACTUAL 1ST DIMENSION OF FLD (.GE.NX)
C       NX,NY   - SIZE OF FLD ARRAY
C       FLDI    - ARRAY OF VALUES FROM WHICH TO INTERPOLATE
C       NDXI    - ACTUAL 1ST DIMENSION OF FLDI (.GE.NXI)
C       NXI,NYI - SIZE OF FLDI ARRAY
C       FLAND   - IF NEAREST POINT IS < FLAND, USE NEAREST POINT
C                  I.E. FAVOR LAND
C
C  4) ALAN J. WALLCRAFT, NRL, JUNE 2000.
C*
C**********
C
C     LOCAL VARIABLES.
C
      REAL, ALLOCATABLE :: FP(:,:),SP(:,:)
C
      INTEGER I,II,IP,IQ,J,JJ,JP,JQ
      REAL*4  FPMIN,FPMAX,FPSUM,FS,RNP
C
      ALLOCATE( FP(-NPATCH:NPATCH,-NPATCH:NPATCH),
     &          SP(-NPATCH:NPATCH,-NPATCH:NPATCH) )
C
      SP( -NPATCH:NPATCH,   -NPATCH:NPATCH)   = 0.5
      SP(1-NPATCH:NPATCH-1,1-NPATCH:NPATCH-1) = 1.0
      RNP = 1.0/SUM(SP(:,:))
C
      DO J= 1,NY
        DO I= 1,NX
          II = NINT(FX(I,J))
          JJ = NINT(FY(I,J))
          IF     (FLDI(II,JJ).LT.FLAND) THEN
            FLD(I,J) = FLDI(II,JJ)  !favor land
          ELSE
            FPMIN =  HUGE(FPMIN)
            FPMAX = -HUGE(FPMAX)
            DO JP = -NPATCH,NPATCH
              JQ = JJ+JP
              IF     (JQ.LT.1 .OR. JQ.GT.NYI) THEN
                WRITE(6,*) 'ERROR - BAD PATCH'
                WRITE(6,*) 'I,J,II,JJ,JQ = ',I,J,II,JJ,JQ
                WRITE(6,*) 
                STOP
              ENDIF
              DO IP = -NPATCH,NPATCH
                IQ = II+IP
                IF     (IQ.LT.1) THEN
                  IQ = NXI-4+IQ  !0 goes to nxi-4
                ELSEIF (IQ.GT.NXI) THEN
                  IQ = IQ-NXI+4  !NXI+1 goes to 5
                ENDIF
                FP(IP,JP) = FLDI(IQ,JQ)
C
                FPMIN = MIN( FPMIN, FP(IP,JP) )
                FPMAX = MAX( FPMAX, FP(IP,JP) )
              ENDDO !ip
            ENDDO !jp
C
            IF     (FPMAX.GT.0.0 .AND. FPMIN.LT.-FPMAX) THEN
C ---         REDUCE THE EFFECT OF HIGH LAND VALUES
              FPSUM = 0.0
              DO JP = -NPATCH,NPATCH
                DO IP = -NPATCH,NPATCH
                  FPSUM = FPSUM + SP(IP,JP)*MAX( FP(IP,JP), -FPMAX )
                ENDDO !ip
              ENDDO !jp
            ELSE
              FPSUM = 0.0
              DO JP = -NPATCH,NPATCH
                DO IP = -NPATCH,NPATCH
                  FPSUM = FPSUM + SP(IP,JP)*FP(IP,JP)
                ENDDO !ip
              ENDDO !jp
            ENDIF
            FLD(I,J) = FPSUM*RNP
          ENDIF
        ENDDO
      ENDDO
      RETURN
C     END OF PATCH.
      END
      SUBROUTINE LINEAR(FLD,FX,FY,NDX,NX,NY, FLDI,NDXI,NXI,NYI)
      IMPLICIT NONE
C
      INTEGER NDX,NX,NY, NDXI,NXI,NYI
      REAL*4  FLD(NDX,NY),FX(NDX,NY),FY(NDX,NY)
      REAL*4  FLDI(NDXI,NYI)
C
C**********
C*
C  1) INTERPOLATE FROM THE ARRAY FLDI TO THE ARRAY FLD, WHERE 
C      FLD(I,J) IS AT (FX(I,J),FY(I,J)) W.R.T. THE FLDI GRID
C      (1:NXI,1:NYI).
C
C     BI-LINEAR INTERPOLATION IS USED.
C     THE GRIDS NEED NOT HAVE IDENTICAL ORIENTATIONS IN SPACE.
C
C     FOR SIMPLICITY, IT IS ASSUMED THAT FX LIES BETWEEN 3 AND NXI-2 
C      AND THAT FY LIES BETWEEN 3 AND NYI-2.  THIS AVOIDS SPECIAL 
C      CASES DURING THE INTERPOLATION, BUT MAY REQUIRE THE CONSTRUCTION 
C      OF A LARGER FLDI ARRAY FROM THE ORIGINAL BEFORE CALLING BESSEL.
C      IF SUCH A LARGER RECTANGLE IS REQUIRED IT SHOULD BE FILLED WITH
C      THE APPROPRIATE VALUES FROM THE ARRAY IN PERIODIC CASES, AND
C      OTHERWISE WITH VALUES LINEARLY EXTRAPOLATED FROM THE ARRAY.
C
C  2) ARGUMENT LIST:
C       FLD     - INTERPOLATED ARRAY ON EXIT
C       FX      - MAPPING OF 1ST DIMENSION OF FLD INTO THAT OF FLDI
C                  FX(I,J).GE.3 .AND. FX(I,J).LE.NXI-2
C       FY      - MAPPING OF 2ND DIMENSION OF FLD INTO THAT OF FLDI
C                  FY(I,J).GE.3 .AND. FY(I,J).LE.NYI-2
C       NDX     - ACTUAL 1ST DIMENSION OF FLD (.GE.NX)
C       NX,NY   - SIZE OF FLD ARRAY
C       FLDI    - ARRAY OF VALUES FROM WHICH TO INTERPOLATE
C       NDXI    - ACTUAL 1ST DIMENSION OF FLDI (.GE.NXI)
C       NXI,NYI - SIZE OF FLDI ARRAY
C
C  4) ALAN J. WALLCRAFT, NRL, JUNE 2000.
C*
C**********
C
C     LOCAL VARIABLES.
C
      INTEGER I,II,J,JJ
      REAL*4  SI,SJ
C
      DO J= 1,NY
        DO I= 1,NX
          II = FX(I,J)
          SI = FX(I,J) - II
          JJ = FY(I,J)
          SJ = FY(I,J) - JJ
          FLD(I,J) = (1.0-SI)*(1.0-SJ)*FLDI(II,  JJ  ) +
     +                    SI *(1.0-SJ)*FLDI(II+1,JJ  ) +
     +               (1.0-SI)*     SJ *FLDI(II,  JJ+1) +
     +                    SI *     SJ *FLDI(II+1,JJ+1)
        ENDDO
      ENDDO
      RETURN
C     END OF LINEAR.
      END
      SUBROUTINE CUBSPL(FLD,FX,FY,NDX,NX,NY,
     +                  FLDI,NDXI,NXI,NYI, IBD, FXI,FYI,WKI,WK)
      IMPLICIT NONE
C
      INTEGER NDX,NX,NY, NXI,NDXI,NYI, IBD(4)
      REAL*4  FLD( NDX, NY), FX(NDX,NY), FY(NDX,NY)
      REAL*4  FLDI(NDXI,NYI),FXI(NXI),FYI(NYI),
     +        WKI(NDXI,NYI,3),WK(3*(NXI+NYI)+1)
C
C**********
C*
C  1) INTERPOLATE FROM THE ARRAY FLDI TO THE ARRAY FLD, WHERE 
C      FLD(I,J) IS AT (FX(I,J),FY(I,J)) W.R.T. THE FLDI GRID (1:NXI,1:NYI).
C
C     CUBIC SPLINE INTERPOLATION IS USED.
C     THE GRIDS NEED NOT HAVE IDENTICAL ORIENTATIONS IN SPACE.
C
C     FOR COMPATABILITY WITH SUBROUTINE 'BESSEL', IT IS ASSUMED THAT
C      FX(I,J) LIES BETWEEN 3 AND NXI-2 AND THAT FY(I,J) LIES BETWEEN
C      3 AND NYI-2.
C
C  2) ARGUMENT LIST:
C       FLD     - INTERPOLATED ARRAY ON EXIT
C       FX      - MAPPING OF 1ST DIMENSION OF FLD INTO THAT OF FLDI
C                  FX(I,J).GE.3 .AND. FX(I,J).LE.NXI-2
C       FY      - MAPPING OF 2ND DIMENSION OF FLD INTO THAT OF FLDI
C                  FY(I,J).GE.3 .AND. FY(I,J).LE.NYI-2
C       NDX     - ACTUAL 1ST DIMENSION OF FLD (.GE.NX)
C       NX,NY   - SIZE OF FLD ARRAY
C       INDX    - ARRAY FOR INTEGER CONSTANTS
C       FLDI    - ARRAY OF VALUES FROM WHICH TO INTERPOLATE
C       NDXI    - ACTUAL 1ST DIMENSION OF FLDI (.GE.NXI)
C       NXI,NYI - SIZE OF FLDI ARRAY
C       IBD     - SPLINE BOUNDARY CONDITION FLAGS
C                  IBD(1:2)=3 FOR PERIODIC REGION, OTHERWISE IBD(:)=2.
C       FXI     - WORKSPACE ARRAY, OVERWRITTEN ON EXIT
C       FYI     - WORKSPACE ARRAY, OVERWRITTEN ON EXIT
C       WKI     - WORKSPACE ARRAY, OVERWRITTEN ON EXIT
C       WK      - WORKSPACE ARRAY, OVERWRITTEN ON EXIT
C
C  3) ALAN J. WALLCRAFT, NRL, JULY 1998.
C*
C**********
C
      INTEGER I,J,NXII,NYII
      REAL*4  TERP2
C
C     BOUNDARY CONDITIONS.
C
      IF     (IBD(1).EQ.3) THEN
        NXII = NXI - 3
        NYII = NYI - 4
      ELSE
        NXII = NXI - 4
        NYII = NYI - 4
        DO J= 1,NYI
          WKI(    3,J,1) = 0.0  !1st derivative assumed zero
          WKI(NXI-2,J,1) = 0.0  !1st derivative assumed zero
        ENDDO
      ENDIF
      DO I= 1,NXI
        WKI(I,    3,2) = 0.0  !1st derivative assumed zero
        WKI(I,NYI-2,2) = 0.0  !1st derivative assumed zero
      ENDDO
C
C     GRID COORDINATES.
C
      DO I= 1,NXI
        FXI(I) = I
      ENDDO
      DO J= 1,NYI
        FYI(J) = J
      ENDDO
C
C     SPLINE COEFFICIENTS.
C
      CALL COEFF2(NXII,FXI(3),NYII,FYI(3),FLDI(3,3),
     +            WKI(3,3,1),WKI(3,3,2),WKI(3,3,3),NDXI, IBD,WK)
C
C     INTERPOLATE.
C
      DO J= 1,NY
        DO I= 1,NX
          FLD(I,J) = TERP2(FX(I,J),FY(I,J),
     +                     NXII,FXI(3),NYII,FYI(3),FLDI(3,3),
     +                     WKI(3,3,1),WKI(3,3,2),WKI(3,3,3),NDXI, 0,0)
        ENDDO
      ENDDO
      RETURN
C     END OF CUBSPL.
      END
C
C CUBSPL     FROM NSSL                                     08/15/79    
C
      SUBROUTINE COEFF1 (N,X,F,W,IOP,INT,WK)    
      IMPLICIT REAL*4 (A-H,O-Z)
C    
C PACKAGE CUBSPL         NOTE--DOCUMENTATION FOR INDIVIDUAL ROUTINES    
C                              FOLLOWS THE GENERAL PACKAGE INFORMATION  
C    
C    
C LATEST REVISION        JANUARY 1978    
C    
C PURPOSE                TO PERFORM ONE AND TWO-DIMENSIONAL CUBIC SPLINE
C                        INTERPOLATION WITH CHOICE OF BOUNDARY    
C                        CONDITIONS.  THE FUNCTION AND SELECTED    
C                        DERIVATIVES MAY BE EVALUATED AT ANY POINT WHERE
C                        INTERPOLATION IS REQUIRED.  IN THE    
C                        TWO-DIMENSIONAL CASE, THE GIVEN DATA POINTS    
C                        MUST BE ON A RECTANGULAR GRID, WHICH NEED NOT  
C                        BE EQUALLY SPACED.  THE PACKAGE CUBSPL CONTAINS
C                        SEVEN ROUTINES.    
C    
C                        SUBROUTINE COEFF1    
C                          COMPUTES THE COEFFICIENTS FOR ONE-DIMENSIONAL
C                          CUBIC SPLINE INTERPOLATION USING ONE OF    
C                          THE FOLLOWING BOUNDARY CONDITIONS AT    
C                          EACH END OF THE RANGE.    
C                            . SECOND DERIVATIVE GIVEN AT BOUNDARY.    
C                            . FIRST DERIVATIVE GIVEN AT BOUNDARY.    
C                            . PERIODIC BOUNDARY CONDITION.    
C                            . FIRST DERIVATIVE DETERMINED BY FITTING A 
C                              CUBIC TO THE FOUR POINTS NEAREST TO THE  
C                              BOUNDARY.    
C    
C                        SUBROUTINE TERP1    
C                          USING THE COEFFICIENTS COMPUTED BY COEFF1,   
C                          THIS ROUTINE EVALUATES THE FUNCTION AND/OR   
C                          FIRST AND SECOND DERIVATIVES AT ANY POINT    
C                          WHERE INTERPOLATION IS REQUIRED.    
C    
C                        SUBROUTINE COEFF2    
C                          COMPUTES THE COEFFICIENTS FOR TWO-DIMENSIONAL
C                          BICUBIC SPLINE INTERPOLATION WITH THE SAME   
C                          CHOICE OF BOUNDARY CONDITIONS AS FOR COEFF1. 
C    
C                        FUNCTION TERP2    
C                          USING THE COEFFICIENTS PRODUCED BY COEFF2,   
C                          THIS SUBROUTINE EVALUATES THE FUNCTION OR A  
C                          SELECTED DERIVATIVE AT ANY POINT WHERE    
C                          TWO-DIMENSIONAL INTERPOLATION IS REQUIRED.   
C    
C                        SUBROUTINE TRIP    
C                          A SIMPLE PERIODIC, TRIDIAGONAL LINEAR    
C                          EQUATION SOLVER USED BY COEFF1.    
C    
C                        SUBROUTINE SEARCH    
C                          PERFORMS A BINARY SEARCH IN A ONE-DIMENSIONAL
C                          FLOATING POINT TABLE ARRANGED IN ASCENDING   
C                          ORDER.  THIS ROUTINE IS CALLED BY TERP1 AND  
C                          TERP2.    
C    
C                        SUBROUTINE INTERP    
C                          GIVEN COEFFICIENTS PROVIDED BY COEFF1 AND THE
C                          POSITION OF THE INTERPOLATION POINT IN THE   
C                          INDEPENDENT VARIABLE TABLE, THIS ROUTINE    
C                          PERFORMS ONE-DIMENSIONAL INTERPOLATION FOR   
C                          THE FUNCTION VALUE, FIRST AND SECOND    
C                          DERIVATIVE, AS DESIRED.  THIS ROUTINE IS    
C                          CALLED BY TERP1 AND TERP2.    
C    
C ACCESS CARDS           *FORTRAN,S=ULIB,N=CUBSPL    
C    
C                          THESE CARDS ACCESS ALL SEVEN ROUTINES OF THE 
C                          PACKAGE, CUBSPL.    
C    
C USAGE                  FOR ONE-DIMENSIONAL CUBIC SPLINE INTERPOLATION,
C                        THE USER FIRST CALLS COEFF1 BY    
C    
C                          CALL COEFF1 (N,X,F,W,IOP,INT,WK)    
C    
C                        THIS SUBROUTINE RETURNS THE COEFFICIENTS    
C                        NEEDED FOR THE SUBSEQUENT INTERPOLATION    
C                        IN THE ARRAY W.  THE USER THEN CALLS    
C                        SUBROUTINE TERP1 BY    
C    
C                          CALL TERP1 (N,X,F,W,Y,INT,TAB,ITAB)    
C    
C                        IN ORDER TO COMPUTE THE VALUE OF THE    
C                        FUNCTION AND/OR ITS DERIVATIVES.  THE USER    
C                        SPECIFIES Y, THE VALUE OF THE INDEPENDENT    
C                        VARIABLE WHERE THE INTERPOLATION IS TO BE    
C                        PERFORMED.  THE INTERPOLATED VALUES ARE    
C                        RETURNED IN TAB.  THE PARAMETERS    
C                        N,X,F,W, AND INT MUST NOT BE CHANGED    
C                        BETWEEN THE SUBROUTINE CALLS.    
C    
C                        FOR TWO-DIMENSIONAL CUBIC SPLINE INTERPOLATION 
C                        THE USER FIRST CALLS COEFF2 BY    
C    
C                          CALL COEFF2 (NX,X,NY,Y,F,FXX,FYY,FXXYY,    
C                                       IDM,IBD,WK)    
C    
C                        THIS SUBROUTINE RETURNS THE COEFFICIENTS    
C                        NEEDED FOR THE SUBSEQUENT INTERPOLATION IN    
C                        THE ARRAYS FXX, FYY, FXXYY.  THE USER THEN    
C                        CALLS THE ROUTINE TERP2 BY    
C    
C                          R = TERP2 (XB,YB,NX,X,NY,Y,F,FXX,FYY,FXXYY,  
C                                     IDM,IXD,IYD)    
C    
C                        IN ORDER TO PERFORM THE INTERPOLATION AT THE   
C                        POINT SPECIFIED BY THE VALUES OF XB AND YB.    
C                        DEPENDING ON THE VALUES INPUT IN IXD AND IYD,  
C                        THE ROUTINE RETURNS AN INTERPOLATED VALUE    
C                        FOR THE FUNCTION OR ONE OF ITS PARTIAL    
C                        DERIVATIVES.  THE PARAMETERS NX,X,NY,Y,F,FXX,  
C                        FYY,FXXYY, AND IDM MUST NOT BE CHANGED    
C                        BETWEEN THE CALLS TO COEFF2 AND TERP2.    
C    
C SPACE REQUIRED         TOTAL SPACE FOR ALL SEVEN ROUTINES IS    
C                        2243 (OCTAL) = 1187 (DECIMAL).    
C    
C ENTRY POINTS           COEFF1, TERP1, COEFF2, TERP2, TRIP, SEARCH,    
C                        INTERP    
C    
C SPECIAL CONDITIONS     TABLES OF INDEPENDENT VARIABLES MUST BE    
C                        ARRANGED IN ASCENDING ORDER.  FOR    
C                        TWO-DIMENSIONAL INTERPOLATION, THE DATA POINTS 
C                        MUST LIE ON A RECTANGULAR MESH.    
C    
C COMMON BLOCKS          NONE    
C    
C I/O                    NONE    
C    
C PRECISION              SINGLE    
C    
C REQUIRED ULIB          NONE    
C ROUTINES    
C    
C SPECIALIST             CICELY RIDLEY, NCAR, BOULDER, COLORADO 80307   
C    
C LANGUAGE               FORTRAN    
C    
C HISTORY                THIS PACKAGE IS BASED ON THE ROUTINES    
C                            LA E 102A, SPL1D1    
C                            LA E 103A, SPL1D2    
C                            LA E 104A, SPL2D1    
C                            LA E 105A, SPL2D2    
C                        OF THE LOS ALAMOS CUBIC SPLINE PACKAGE WRITTEN 
C                        BY THOMAS J. JORDAN AND BERTHA FAGAN, 1968.    
C                        THE ROUTINES HAVE BEEN STREAMLINED AND    
C                        STANDARDIZED.  THE ALGORITHM FOR    
C                        TWO-DIMENSIONAL INTERPOLATION IS CONSIDERABLY  
C                        MODIFIED.    
C    
C ALGORITHM              FOR ONE-DIMENSIONAL INTERPOLATION, THE CUBIC   
C                        SPLINE IS EXPRESSED IN TERMS OF THE FUNCTION   
C                        VALUES AND SECOND DERIVATIVES AT THE DATA    
C                        POINTS.  THE SECOND DERIVATIVES ARE    
C                        DETERMINED FROM A TRIDIAGONAL LINEAR SYSTEM    
C                        WHICH DESCRIBES THE CONTINUITY OF THE FIRST    
C                        DERIVATIVE AND INCORPORATES THE GIVEN    
C                        BOUNDARY CONDITIONS.  COEFF1  SETS UP THIS    
C                        SYSTEM AND CALLS SUBROUTINE TRIP TO SOLVE IT.  
C    
C                        THE CUBIC SEGMENT BETWEEN TWO ADJACENT    
C                        TABULAR POINTS IS CONSTRUCTED FROM THE    
C                        FUNCTION VALUES AND SECOND DERIVATIVES AT    
C                        THESE POINTS.  THESE PROVIDE THE FOUR    
C                        CONSTANTS NEEDED TO DEFINE THE CUBIC    
C                        UNIQUELY.  FROM THIS CUBIC, VALUES OF THE    
C                        FUNCTION AND ITS FIRST AND SECOND    
C                        DERIVATIVES ARE READILY DETERMINED AT ANY    
C                        INTERMEDIATE POINT.  ONE-DIMENSIONAL    
C                        INTERPOLATION IS PERFORMED BY THE ROUTINE    
C                        TERP1.  FOR TWO-DIMENSIONAL INTERPOLATION,    
C                        THE BICUBIC SPLINE IS DESCRIBED IN TERMS OF    
C                        VALUES OF F,FXX,FYY, AND FXXYY  AT EACH    
C                        POINT ON THE GIVEN TWO-DIMENSIONAL    
C                        RECTANGULAR GRID OF DATA POINTS.  HERE F    
C                        IS THE FUNCTION VALUE,    
C    
C                          FXX = (D/DX)**2*F    
C    
C                        AND SO ON.  THE COEFFICIENTS ARE DETERMINED    
C                        BY COEFF2, WHICH USES SUCCESSIVE    
C                        APPLICATIONS OF COEFF1.    
C    
C                        1.  THE ARRAY FXX IS DETERMINED BY APPLYING    
C                            COEFF1 TO F ALONG EACH LINE IN THE    
C                            X-DIRECTION.    
C    
C                        2.  THE ARRAY FYY IS DETERMINED BY APPLYING    
C                            COEFF1 TO F ALONG EACH LINE IN THE    
C                            Y-DIRECTION.    
C    
C                        3.  FXXYY IS DETERMINED ON THE CONSTANT Y    
C                            BOUNDARIES BY APPLYING COEFF1 TO FYY    
C                            ALONG THESE BOUNDARIES.    
C    
C                        4.  THE REMAINDER OF THE ARRAY FXXYY IS    
C                            DETERMINED BY APPLYING COEFF1 TO FXX    
C                            ALONG EACH LINE IN THE Y-DIRECTION.    
C    
C                        THE BICUBIC WITHIN ANY RECTANGULAR ELEMENT    
C                        OF THE GRID IS CONSTRUCTED FROM THE VALUES    
C                        OF F,FXX,FYY,FXXYY AT THE FOUR CORNERS.    
C                        THESE PROVIDE THE 16 CONSTANTS NECESSARY    
C                        TO DEFINE THE BICUBIC UNIQUELY.  TO FIND    
C                        THE VALUE OF F CORRESPONDING TO A POINT    
C                        (XB,YB) WITHIN THE ELEMENTARY RECTANGLE,    
C                        (X(I),Y(J)),(X(I+1),Y(J)),(X(I),Y(J+1)),    
C                        (X(I+1),Y(J+1)), FIVE ONE DIMENSIONAL    
C                        INTERPOLATIONS ARE PERFORMED.    
C    
C                        1.  F AT (XB,Y(J)) AND (XB,Y(J+1)) ARE    
C                            FOUND BY INTERPOLATING F IN THE    
C                            X-DIRECTION USING FXX. (TWO INTERPOLATIONS)
C    
C                        2.  FYY AT (XB,Y(J)) AND (XB,Y(J+1)) ARE    
C                            FOUND BY INTERPOLATING FYY IN THE    
C                            X-DIRECTION USING FXXYY. (TWO    
C                            INTERPOLATIONS.)    
C    
C                        3.  FINALLY F AT (XB,YB) IS FOUND BY    
C                        INTERPOLATING BETWEEN F(XB,Y(J)) AND    
C                        F(XB,Y(J+1)) IN THE Y-DIRECTION USING    
C                        VALUES OF FYY(XB,Y(J)) AND FYY(XB,Y(J+1))    
C                        OBTAINED ABOVE. (ONE INTERPOLATION).    
C    
C                        TWO-DIMENSIONAL INTERPOLATION IS PERFORMED    
C                        IN TERP2.    
C    
C                        FOR GREATER DETAIL, SEE J.L.WALSH,    
C                        J.H.AHLBERG, E.N.NILSEN, BEST APPROXIMATION    
C                        PROPERTIES OF THE SPLINE FIT, JOURNAL OF    
C                        MATHEMATICS AND MECHANICS, VOL.II(1962),    
C                        225-234.    
C    
C                        T.L. JORDAN, SMOOTHING AND MULTIVARIABLE    
C                        INTERPOLATION WITH SPLINES, LOS ALAMOS    
C                        REPORT, LA-3137, 1965.    
C    
C ACCURACY               NEAR MACHINE ACCURACY WAS OBTAINED WHEN    
C                        INTERPOLATING A CUBIC IN ONE DIMENSION    
C                        OR A BICUBIC IN TWO DIMENSIONS.    
C    
C PORTABILITY            FULLY PORTABLE    
C    
C REQUIRED RESIDENT      NONE    
C ROUTINES    
C    
C    
C    
C    
C    
C    
C    
C    
C    
C    
C    
C    
C SUBROUTINE COEFF1 (N,X,F,W,IOP,INT,WK)    
C    
C DIMENSION OF           X(N),F(INT*(N-1)+1),W(INT*(N-1)+1),IOP(2),    
C ARGUMENTS              WK(3*N+1)    
C    
C LATEST REVISION        JANUARY 1978    
C    
C PURPOSE                SUBROUTINE COEFF1 COMPUTES THE COEFFICIENTS    
C                        FOR ONE-DIMENSIONAL CUBIC SPLINE    
C                        INTERPOLATION USING ONE OF THE FOLLOWING    
C                        BOUNDARY CONDITIONS AT EACH END OF THE    
C                        RANGE    
C    
C                        .  SECOND DERIVATIVES GIVEN AT BOUNDARY    
C                        .  FIRST DERIVATIVE GIVEN AT BOUNDARY    
C                        .  PERIODIC BOUNDARY CONDITIONS    
C                        .  FIRST DERIVATIVE CALCULATED BY FITTING A    
C                           CUBIC TO THE FOUR POINTS NEAREST TO THE    
C                           BOUNDARY    
C    
C                        NOTE THAT TERP1 MUST BE CALLED TO PERFORM    
C                        THE INTERPOLATION.    
C    
C ACCESS CARD            *FORTRAN,S=ULIB,N=CUBSPL    
C    
C USAGE                  CALL COEFF1 (N,X,F,W,IOP,INT,WK)    
C    
C ARGUMENTS    
C    
C ON INPUT               N    
C                          THE NUMBER OF DATA POINTS.  N MUST BE AT    
C                          LEAST 4.    
C    
C                        X    
C                          TABLE OF N INDEPENDENT VARIABLE VALUES IN    
C                          ASCENDING ORDER.  DIMENSION OF X IN CALLING  
C                          PROGRAM MUST BE AT LEAST N.    
C    
C                        F    
C                          TABLE OF N CORRESPONDING DEPENDENT VARIABLE  
C                          VALUES.  THE VALUES ARE SEPARATED BY INTERVAL
C                          INT.  THIS IS USUALLY UNITY FOR    
C                          ONE-DIMENSIONAL INTERPOLATION.  OTHER VALUES 
C                          ARE USEFUL WHEN COEFF1 IS INCORPORATED IN A  
C                          TWO-DIMENSIONAL INTERPOLATION SCHEME (SEE    
C                          COEFF2).  DIMENSION OF F IN THE CALLING    
C                          PROGRAM MUST BE AT LEAST (INT*(N-1)+1).    
C    
C                        IOP    
C                          TWO ELEMENT INTEGER ARRAY DEFINING BOUNDARY  
C                          CONDITIONS AT X(1) AND X(N) ACCORDING TO THE 
C                          FOLLOWING CODE.    
C    
C                          FOR IOP(1)    
C                          = 1  SECOND DERIVATIVE GIVEN AT X(1).  PLACE 
C                               VALUE OF SECOND DERIVATIVE IN W(1)    
C                               BEFORE CALL TO COEFF1.    
C                          = 2  FIRST DERIVATIVE GIVEN AT X(1).  PLACE  
C                               VALUE OF FIRST DERIVATIVE IN W(1) BEFORE
C                               CALL.    
C                          = 3  PERIODIC BOUNDARY CONDITION.  X(1) AND  
C                               X(N) ARE EQUIVALENT POINTS.  F(1) AND   
C                               F(INT*(N-1)+1) ARE EQUAL.    
C                          = 4  THE FIRST DERIVATIVE AT X(1) IS    
C                               CALCULATED BY FITTING A CUBIC TO POINTS 
C                               X(1) THROUGH X(4).    
C                          SIMILARLY, IOP(2) DEFINES THE BOUNDARY    
C                          CONDITION AT X(N).  WHEN IOP(2) = 1 (OR 2),  
C                          THE VALUE OF THE SECOND (OR FIRST) DERIVATIVE
C                          MUST BE PLACED IN W(INT*(N-1)+1).  NOTE THAT 
C                          IF IOP(1) = 3, CONSISTENCY DEMANDS THAT    
C                          IOP(2) = 3 ALSO.    
C    
C                        INT    
C                          SPACING IN F AND W TABLES.  FOR    
C                          ONE-DIMENSIONAL INTERPOLATION THIS WILL    
C                          USUALLY BE UNITY.    
C    
C                        WK    
C                          WORK AREA OF DIMENSION AT LEAST (3*N+1).    
C    
C ON OUTPUT              W    
C                          TABLE OF SECOND DERIVATIVES CORRESPONDING TO 
C                          GIVEN X AND F VALUES.  THE SEPARATION OF    
C                          TABULAR ENTRIES IS INT (SEE ABOVE).    
C                          DIMENSION OF W IN THE CALLING PROGRAM MUST BE
C                          AT LEAST (INT*(N-1)+1).    
C    
C                          THE ARRAYS X, F, W ARE USED AS INPUT FOR THE 
C                          ROUTINE TERP1, WHICH PERFORMS INTERPOLATION  
C                          AT A GIVEN VALUE OF THE INDEPENDENT VARIABLE.
C    
C SPECIAL CONDITIONS     NONE    
C    
C COMMON BLOCKS          NONE    
C    
C I/O                    NONE    
C    
C PRECISION              SINGLE    
C    
C SPECIALIST             CICELY RIDELY, NCAR, BOULDER, COLORADO 80307   
C    
C LANGUAGE               FORTRAN    
C    
C HISTORY                THIS ROUTINE IS A STREAMLINED AND    
C                        STANDARDIZED VERSION OF A SIMILAR ROUTINE    
C                        FROM THE LOS ALAMOS CUBIC SPLINE PACKAGE.    
C    
C SPACE REQUIRED         654 (OCTAL) = 428 (DECIMAL)    
C    
C TIMING                 THE TIMING IS LINEARLY PROPORTIONAL TO N, THE  
C                        NUMBER OF DATA POINTS.  FOR 21 DATA POINTS, THE
C                        TIME ON THE NCAR CDC 7600 WAS APPROXIMATELY    
C                        .24 MILLISECONDS.    
C    
C PORTABILITY            FULLY PORTABLE    
C    
C REQUIRED RESIDENT      NONE    
C ROUTINES    
C    
C    
C    
C    
C    
C    
C    
C    
      DIMENSION       X(*)       ,F(*)       ,W(*)       ,IOP(2)     ,  
     1                WK(N,*)    
C    
C ARITHMETIC STATENENT FUNCTION USED TO LOCATE ENTRIES IN F AND W ARRAYS
C    
      II(INDEX)=(INDEX-1)*INT+1    
C    
C    
C    
C    
C    
C    
C THE FOLLOWING CALL IS FOR GATHERING STATISTICS ON LIBRARY USE AT NCAR 
C.....CALL Q8QST4( 4HNSSL      , 6HCUBSPL    ,6HCOEFF1 ,10HVERSION 11)  
C    
C START TO SET UP TRIDIAGONAL SYSTEM    
C    
      J0 = 1    
      DO 101 I=2,N    
         JM = J0    
         J0 = J0+INT    
         WK(I,1) = X(I)-X(I-1)    
         WK(I,2) = (F(J0)-F(JM))/WK(I,1)    
         WK(I,3) = WK(I,1)/6.    
         WK(I,1) = WK(I,1)/3.    
  101 CONTINUE    
      NN = N    
      MK = IOP(1)    
      ML = IOP(2)    
C    
C APPLY BOUNDARY CONDITIONS AT BOUNDARY 1    
C    
      GO TO (102,103,104,105),MK    
C    
C SECOND DERIVATIVE GIVEN AT BOUNDARY 1    
C    
  102 CONTINUE    
      WK(2,2) = WK(3,2)-WK(2,2)-WK(2,3)*W(1)    
      WK(2,3) = 0.    
      WK(2,1) = WK(2,1)+WK(3,1)    
      I1 = 2    
      NN = NN-1    
      GO TO 106    
C    
C FIRST DERIVATIVE GIVEN AT BOUNDARY 1    
C    
  103 CONTINUE    
      WK(1,2) = WK(2,2)-W(1)    
      WK(2,2) = WK(3,2)-WK(2,2)    
      WK(1,3) = 0.    
      WK(1,1) = WK(2,1)    
      WK(2,1) = WK(2,1)+WK(3,1)    
      I1 = 1    
      GO TO 106    
C    
C PERIODIC BOUNDARY CONDITION    
C    
  104 CONTINUE    
      Y2 = WK(2,2)    
      B2 = WK(2,1)    
      WK(2,2) = WK(3,2)-WK(2,2)    
      WK(2,1) = WK(3,1)+WK(2,1)    
      I1 = 2    
      NN = NN-1    
      GO TO 106    
C    
C FIRST DERIVATIVE AT BOUNDARY 1 FROM 4 POINT INTERPOLATION.    
C    
  105 CONTINUE    
      A12 = X(1)-X(2)    
      A13 = X(1)-X(3)    
      A14 = X(1)-X(4)    
      A23 = X(2)-X(3)    
      A24 = X(2)-X(4)    
      A34 = X(3)-X(4)    
      J1 = 1    
      J2 = J1+INT    
      J3 = J2+INT    
      J4 = J3+INT    
      W(1)    = (1./A12+1./A13+1./A14)*F(J1)-    
     1          A13*A14/(A12*A23*A24)*F(J2)+A12*A14/(A13*A23*A34)*F(J3)-
     2          A12*A13/(A14*A24*A34)*F(J4)    
      GO TO 103    
C COMPUTE TRIDIAGONAL ARRAYS    
  106 CONTINUE    
      I2 = N-2    
      DO 107 I=3,I2    
         WK(I,2) = WK(I+1,2)-WK(I,2)    
         WK(I,1) = WK(I+1,1)+WK(I,1)    
  107 CONTINUE    
C    
C APPLY BOUNDARY CONDITIONS AT BOUNDARY 2.    
C    
      IN = II(N)    
      GO TO (108,109,110,111),ML    
C    
C SECOND DERIVATIVE GIVEN AT BOUNDARY 2.    
C    
  108 CONTINUE    
      WK(N-1,2) = WK(N,2)-WK(N-1,2)-WK(N,3)*W(IN)    
      WK(N,3) = 0.    
      WK(N-1,1) = WK(N-1,1)+WK(N,1)    
      NN = NN-1    
      GO TO 112    
C    
C FIRST DERIVATIVE GIVEN AT BOUNDARY 2.    
C    
  109 CONTINUE    
      WK(N-1,2) = WK(N,2)-WK(N-1,2)    
      WK(N,2) = -WK(N,2)+W(IN)    
      WK(N-1,1) = WK(N-1,1)+WK(N,1)    
      WK(1,4) = 0.    
      GO TO 112    
C    
C PERIODIC BOUNDARY CONDITION    
C    
  110 CONTINUE    
      WK(N-1,2) = WK(N,2)-WK(N-1,2)    
      WK(N,2) = Y2-WK(N,2)    
      WK(N-1,1) = WK(N-1,1)+WK(N,1)    
      WK(N,1) = WK(N,1)+B2    
      WK(1,4) = WK(2,3)    
      GO TO 112    
C    
C FIRST DERIVATIVE AT BOUNDARY 2 FROM 4 POINT INTERPOLATION.    
C    
  111 CONTINUE    
      A12 = X(N)-X(N-1)    
      A13 = X(N)-X(N-2)    
      A14 = X(N)-X(N-3)    
      A23 = X(N-1)-X(N-2)    
      A24 = X(N-1)-X(N-3)    
      A34 = X(N-2)-X(N-3)    
      J1 = IN    
      J2 = J1-INT    
      J3 = J2-INT    
      J4 = J3-INT    
      W(IN)   = (1./A12+1./A13+1./A14)*F(J1)-    
     1          A13*A14/(A12*A23*A24)*F(J2)+A12*A14/(A13*A23*A34)*F(J3)-
     2          A12*A13/(A14*A24*A34)*F(J4)    
      GO TO 109    
  112 CONTINUE    
      IW1 = II(I1)    
      CALL TRIP (NN,WK(I1,3),WK(I1,1),WK(I1+1,3),WK(I1,2),W(IW1),INT)   
      GO TO (114,114,113,114),MK    
  113 CONTINUE    
      W(1) = W(IN)    
  114 CONTINUE    
      RETURN    
      END    
      SUBROUTINE TRIP (N,A,B,C,Y,Z,INT)    
      IMPLICIT REAL*4 (A-H,O-Z)
      DIMENSION       A(*)       ,B(*)       ,C(*)       ,Y(*)       ,  
     1                Z(*)    
C    
C ARITHMETIC STATEMENT FUNCTION USED TO LOCATE ENTRIES IN ARRAY Z.    
C    
      II(INDEX)=(INDEX-1)*INT+1    
C    
C GAUSSIAN ELIMINATION    
C    
      BN = B(N)    
      YN = Y(N)    
      V = C(N)    
      Y(1) = Y(1)/B(1)    
      A(1) = A(1)/B(1)    
      B(1) = C(1)/B(1)    
      NM2 = N-2    
      DO 101 J=2,NM2    
         DEN = B(J)-A(J)*B(J-1)    
         B(J) = C(J)/DEN    
         Y(J) = (Y(J)-A(J)*Y(J-1))/DEN    
         A(J) = -A(J)*A(J-1)/DEN    
         BN = BN-V*A(J-1)    
         YN = YN-V*Y(J-1)    
         V = -V*B(J-1)    
  101 CONTINUE    
      DEN = B(N-1)-A(N-1)*B(N-2)    
      B(N-1) = (C(N-1)-A(N-1)*A(N-2))/DEN    
      Y(N-1) = (Y(N-1)-A(N-1)*Y(N-2))/DEN    
      BN = BN-V*A(N-2)    
      YN = YN-V*Y(N-2)    
      V = A(N)-V*B(N-2)    
C BACK SUBSTITUTION    
      IIN = II(N)    
      Z(IIN) = (YN-V*Y(N-1))/(BN-V*B(N-1))    
      IIN2 = II(N-1)    
      Z(IIN2) = Y(N-1)-B(N-1)*Z(IIN)    
      NM1 = N-1    
      IN = II(N)    
      DO 102 J=2,NM1    
         K = N-J    
         IK = II(K)    
      IKT = IK+INT    
  102 Z(IK) = Y(K)-B(K)*Z(IKT)-A(K)*Z(IN)    
      RETURN    
      END    
      SUBROUTINE SEARCH (XBAR,X,N,I)    
      IMPLICIT REAL*4 (A-H,O-Z)
      DIMENSION       X(*)    
      DATA B/.69314718/    
C    
C IF XBAR IS OUTSIDE RANGE OF X TABLE EXTRAPOLATE    
C    
      IF (XBAR .GT. X(2)) GO TO 101    
      I = 1    
      RETURN    
  101 CONTINUE    
      IF (XBAR .LT. X(N-1)) GO TO 102    
      I = N-1    
      RETURN    
  102 CONTINUE    
C    
C FIND MAXIMUM POWER OF TWO LESS THAN N    
C    
      M = INT((ALOG(FLOAT(N)))/B)    
      I = 2**M    
      IF (I .GE. N) I = I/2    
      K = I    
      NM1 = N-1    
C    
C CONDUCT BINARY SEARCH.    
C    
  103 CONTINUE    
      K = K/2    
      IF (XBAR .GE. X(I)) GO TO 104    
      I = I-K    
      GO TO 103    
  104 CONTINUE    
      IF (XBAR .LE. X(I+1)) RETURN    
      I = MIN0(I+K,NM1)    
      GO TO 103    
      END    
      SUBROUTINE INTERP (N,X,F,W,Y,I,INT,TAB,ITAB)    
      IMPLICIT REAL*4 (A-H,O-Z)
      DIMENSION       X(*)       ,F(*)       ,W(*)       ,TAB(3)    
     -                ,ITAB(3)    
C    
C ARITHMETIC STATEMENT FUNCTION USED TO LOCATE ENTRIES IN F AND W ARRAYS
C    
      II(INDEX)=(INDEX-1)*INT+1    
C    
C PERFORM INTERPOLATION OR EXTRAPOLATION    
C    
      FLK = X(I+1)-X(I)    
      FLP = X(I+1)-Y    
      FL0 = Y-X(I)    
      I0 = II(I)    
      IP = I0+INT    
      IF (ITAB(1) .NE. 0) GO TO 101    
      GO TO 102    
  101 CONTINUE    
C    
C CALCULATE F(Y)    
C    
      A = (W(I0)*FLP**3+W(IP)*FL0**3)/(6.*FLK)    
      B = (F(IP)/FLK-W(IP)*FLK/6.)*FL0    
      C = (F(I0)/FLK-W(I0)*FLK/6.)*FLP    
      TAB(1) = A+B+C    
  102 CONTINUE    
      IF (ITAB(2) .NE. 0) GO TO 103    
      GO TO 104    
  103 CONTINUE    
C    
C CALCULATE FIRST DERIVATIVE AT Y    
C    
      A = (W(IP)*FL0**2-W(I0)*FLP**2)/(2.*FLK)    
      B = (F(IP)-F(I0))/FLK    
      C = (W(I0)-W(IP))*FLK/6.    
      TAB(2) = A+B+C    
  104 CONTINUE    
      IF (ITAB(3) .NE. 0) GO TO 105    
      GO TO 106    
  105 CONTINUE    
C    
C CALCULATE SECOND DERIVATIVE AT Y    
C    
      TAB(3) = (W(I0)*FLP+W(IP)*FL0)/FLK    
  106 CONTINUE    
      RETURN    
      END    
      SUBROUTINE TERP1 (N,X,F,W,Y,INT,TAB,ITAB)    
      IMPLICIT REAL*4 (A-H,O-Z)
C    
C    
C DIMENSION OF           X(N),F(INT*(N-1)+1),W(INT*(N-1)+1),TAB(3),    
C ARGUMENTS              ITAB(3)    
C    
C LATEST REVISION        JANUARY 1978    
C    
C PURPOSE                USING THE COEFFICIENTS COMPUTED BY COEFF1,    
C                        THIS ROUTINE EVALUATES THE FUNCTION AND/OR    
C                        FIRST AND SECOND DERIVATIVES AT ANY POINT.    
C    
C ACCESS CARD            *FORTRAN,S=ULIB,N=CUBSPL    
C    
C USAGE                  CALL TERP1 (N,X,F,W,Y,INT,TAB,ITAB)    
C    
C ARGUMENTS    
C    
C ON INPUT               N    
C                          THE NUMBER OF DATA POINTS.  N MUST BE AT    
C                          LEAST 4.    
C    
C                        X    
C                          TABLE OF N INDEPENDENT VARIABLE VALUES IN    
C                          ASCENDING ORDER.  DIMENSION OF X IN THE    
C                          CALLING PROGRAM MUST BE AT LEAST N.    
C    
C                        F    
C                          TABLE OF N CORRESPONDING DEPENDENT VARIABLE  
C                          VALUES SEPARATED BY INTERVAL INT, USUALLY    
C                          UNITY FOR ONE-DIMENSIONAL INTERPOLATION.    
C                          DIMENSION OF F IN THE CALLING PROGRAM MUST BE
C                          AT LEAST (INT*(N-1)+1).    
C    
C                        W    
C                          TABLE OF SECOND DERIVATIVES COMPUTED BY    
C                          COEFF1.  THE SEPARATION OF TABULAR ENTRIES IS
C                          INT.  DIMENSION OF W IN THE CALLING PROGRAM  
C                          MUST BE AT LEAST (INT*(N-1)+1).    
C    
C                        Y    
C                          VALUE OF THE INDEPENDENT VARIABLE AT WHICH   
C                          INTERPOLATION IS REQUIRED.  IF Y LIES OUTSIDE
C                          THE RANGE OF THE TABLE, EXTRAPOLATION TAKES  
C                          PLACE.    
C    
C                        INT    
C                          SPACING OF TABULAR ENTRIES IN F AND W ARRAYS.
C                          THIS IS USUALLY UNITY FOR ONE-DIMENSIONAL    
C                          INTERPOLATION.    
C    
C                        ITAB    
C                          THREE ELEMENT INTEGER ARRAY DEFINING    
C                          INTERPOLATION TO BE PERFORMED AT Y.    
C                            IF ITAB(1) = 1, THE FUNCTION VALUE IS    
C                                            RETURNED IN TAB(1).    
C                            IF ITAB(2) = 1, THE FIRST DERIVATIVE IS    
C                                            RETURNED IN TAB(2).    
C                            IF ITAB(3) = 1, THE SECOND DERIVATIVE IS   
C                                            RETURNED IN TAB(3).    
C                          IF ITAB(I) = 0 FOR I = 1, 2 OR 3, THE    
C                          CORRESPONDING FUNCTION VALUE OR DERIVATIVE IS
C                          NOT COMPUTED AND TAB(I) IS NOT REFERENCED.   
C    
C ON OUTPUT              TAB    
C                          THREE ELEMENT ARRAY IN WHICH INTERPOLATED    
C                          FUNCTION VALUE, FIRST AND SECOND DERIVATIVES 
C                          ARE RETURNED AS DICTATED BY ITAB (SEE ABOVE).
C    
C SPECIAL CONDITIONS     NONE    
C    
C COMMON BLOCKS          NONE    
C    
C I/O                    NONE    
C    
C PRECISION              SINGLE    
C    
C SPECIALIST             CICELY RIDELY, NCAR, BOULDER, COLORADO 80307   
C    
C LANGUAGE               FORTRAN    
C    
C HISTORY                THIS ROUTINE IS A STREAMLINED AND    
C                        STANDARDIZED VERSION OF A SIMILAR ROUTINE    
C                        FROM THE LOS ALAMOS CUBIC SPLINE PACKAGE.    
C    
C SPACE REQUIRED         47 (OCTAL) = 39 (DECIMAL)    
C    
C TIMING                 THIS PROCEDURE IS FAST.  THE MAXIMUM TIME FOR  
C                        THE BINARY SEARCH IS PROPORTIONAL TO ALOG(N).  
C                        THE TIME FOR FUNCTION AND DERIVATIVE EVALUATION
C                        IS INDEPENDENT OF N.  FOR 21 DATA POINTS THE   
C                        TIME TAKEN ON THE NCAR 7600 IS ABOUT    
C                        80 MICROSECONDS.    
C    
C PORTABILITY            FULLY PORTABLE    
C    
C REQUIRED RESIDENT      ALOG    
C ROUTINES    
C    
C    
C    
C    
C    
C    
C    
      DIMENSION       X(*)       ,F(*)       ,W(*)       ,TAB(3)     ,  
     1                ITAB(3)    
C THE FOLLOWING CALL IS FOR GATHERING STATISTICS ON LIBRARY USE AT NCAR 
C.....CALL Q8QST4( 4HNSSL      , 6HCUBSPL    ,5HTERP1  ,10HVERSION 11)  
C    
C PERFORM SEARCH    
C    
      CALL SEARCH (Y,X,N,I)    
C    
C CARRY OUT INTERPOLATION (OR EXTRAPOLATION)    
C    
      CALL INTERP (N,X,F,W,Y,I,INT,TAB,ITAB)    
      RETURN    
      END    
      SUBROUTINE COEFF2 (NX,X,NY,Y,F,FXX,FYY,FXXYY,IDM,IBD,WK)    
      IMPLICIT REAL*4 (A-H,O-Z)
C    
C    
C DIMENSION OF           X(NX),Y(NY),F(IDM,NY),FXX(IDM,NY),FYY(IDM,NY), 
C ARGUMENTS              FXXYY(IDM,NY),IBD(4),WK(3*MAX0(NX,NY)+1)    
C                        (IDM MUST BE .GE. NX)    
C    
C LATEST REVISION        JANUARY 1978    
C    
C PURPOSE                SUBROUTINE COEFF2 COMPUTES THE COEFFICIENTS    
C                        FOR TWO-DIMENSIONAL BICUBIC SPLINE    
C                        INTERPOLATION WITH THE SAME CHOICE OF    
C                        BOUNDARY CONDITIONS AS FOR COEFF1.  TERP2    
C                        IS CALLED TO PERFORM INTERPOLATION.    
C    
C ACCESS CARD            *FORTRAN,S=ULIB,N=CUBSPL    
C    
C USAGE                  CALL COEFF2 (NX,X,NY,Y,F,FXX,FYY,FXXYY,IDM,IBD,
C                                     WK)    
C    
C ARGUMENTS    
C    
C ON INPUT               NX    
C                          NUMBER OF GRID POINTS IN THE X-DIRECTION.  NX
C                          MUST BE AT LEAST 4.    
C    
C                        X    
C                          TABLE OF NX VALUES OF THE FIRST INDEPENDENT  
C                          VARIABLE ARRANGED IN ASCENDING ORDER.    
C                          DIMENSION OF X IN THE CALLING PROGRAM MUST BE
C                          AT LEAST NX.    
C    
C                        NY    
C                          NUMBER OF GRID POINTS IN THE Y-DIRECTION.  NY
C                          MUST BE AT LEAST 4.    
C    
C                        Y    
C                          TABLE OF NY VALUES OF THE SECOND INDEPENDENT 
C                          VARIABLE ARRANGED IN ASCENDING ORDER.    
C                          DIMENSION OF Y IN THE CALLING PROGRAM MUST BE
C                          AT LEAST NY.    
C    
C                        F    
C                          TWO DIMENSIONAL ARRAY OF FUNCTION VALUES AT  
C                          THE GRID POINTS DEFINED BY THE ARRAYS X AND  
C                          Y.  DIMENSION OF F IN THE CALLING PROGRAM IS 
C                          (IDM, NYY) WHERE    
C                              IDM .GE. NX    
C                              NYY .GE. NY    
C    
C                        IDM    
C                          FIRST DIMENSION IN THE CALLING PROGRAM OF    
C                          ARRAYS F, FXX, FYY, FXXYY.  IDM MUST BE AT   
C                          LEAST NX.    
C    
C                        IBD    
C                          FOUR ELEMENT INTEGER ARRAY DEFINING BOUNDARY 
C                          CONDITIONS ACCORDING TO THE FOLLOWING CODE.  
C                          FOR IBD(1)    
C                          = 1  THE SECOND DERIVATIVE OF F WITH RESPECT 
C                               TO X IS GIVEN AT (X(1),Y(J)) FOR    
C                               J = 1,NY,1.  VALUES OF THIS SECOND    
C                               DERIVATIVE MUST BE PLACED IN FXX(1,J)   
C                               FOR J = 1,NY,1, BEFORE CALLING COEFF2.  
C                          = 2  THE FIRST DERIVATIVE OF F WITH RESPECT  
C                               TO X IS GIVEN AT (X(1),Y(J)) FOR    
C                               J = 1,NY,1.  VALUES OF THE DERIVATIVE   
C                               MUST BE PLACED IN FXX(1,J) FOR    
C                               J = 1,NY,1 BEFORE CALLING COEFF2.    
C                          = 3  PERIODIC BOUNDARY CONDITION IN THE    
C                               X-DIRECTION.  (X(1),Y(J)) AND    
C                               AND (X(NX),Y(J)) ARE EQUIVALENT POINTS  
C                               FOR J = 1,NY,1.  F(1,J) AND F(NX,J) ARE 
C                               EQUAL.    
C                          = 4  THE FIRST DERIVATIVE OF F WITH RESPECT  
C                               TO X AT (X(1),Y(J)) IS COMPUTED BY    
C                               FITTING A CUBIC TO F(1,J) THROUGH F(4,J)
C                               FOR J = 1,NY,1.    
C    
C                          SIMILARLY, IBD(2) DEFINES THE BOUNDARY    
C                          CONDITION AT (X(NX),Y(J)) FOR J = 1,NY,1.    
C                          WHEN IBD(2) = 1 (OR 2) THE VALUES OF THE    
C                          SECOND (OR FIRST) DERIVATIVE OF F WITH    
C                          RESPECT TO X ARE PLACED IN FXX(NX,J) FOR    
C                          J = 1,NY,1.    
C                            NOTE THAT IF IBD(1) = 3, CONSISTENCY    
C                            REQUIRES THAT IBD(2) = 3 ALSO.    
C                          FOR IBD(3)    
C                          = 1  THE SECOND DERIVATIVE OF F WITH RESPECT 
C                               TO Y IS GIVEN AT (X(I),Y(1)).  PLACE    
C                               VALUES OF THE DERIVATIVE IN FYY(I,1) FOR
C                               I = 1,NX,1 BEFORE CALLING COEFF2.    
C                          = 2  THE FIRST DERIVATIVE OF F WITH RESPECT  
C                               TO Y IS GIVEN AT (X(I),Y(1)).  VALUES OF
C                               THIS DERIVATIVE MUST BE PLACED IN    
C                               FYY(I,1) FOR I = 1,NX,1 BEFORE CALLING  
C                               COEFF2.    
C                          = 3  PERIODIC BOUNDARY CONDITION IN THE    
C                               Y-DIRECTION.  (X(I),Y(1)) AND    
C                               (X(I),Y(NY)) ARE EQUIVALENT POINTS.    
C                               F(I,1) AND F(I,NY) ARE EQUAL.    
C                          = 4  THE FIRST DERIVATIVE OF F WITH RESPECT  
C                               TO Y AT (X(I),Y(1)) IS COMPUTED BY    
C                               FITTING A CUBIC TO F(I,1) THROUGH F(I,4)
C                               FOR I = 1,NX,1.    
C    
C                          SIMILARY, IBD(4) DEFINES THE BOUNDARY    
C                          CONDITION AT (X(I),Y(NY)) FOR I = 1,NX,1 AND 
C                          GIVEN DERIVATIVE VALUES ARE PLACED IN    
C                          FYY(I,NY).    
C                            NOTE THAT CONSISTENCY DEMANDS THAT IF    
C                            IBD(3) = 3, THEN IBD(4) = 3 ALSO.    
C    
C                        WK    
C                          WORK AREA OF DIMENSION AT LEAST    
C                          (3*MAX0(NX,NY)+1)    
C    
C ON OUTPUT              FXX    
C                          ARRAY OF SECOND DERIVATIVES OF F WITH RESPECT
C                          TO X COMPUTED BY COEFF2.  FXX(I,J) IS    
C                          DERIVATIVE AT (X(I),Y(J)).  AS FOR F,    
C                          DIMENSION OF FXX IN THE CALLING PROGRAM IS   
C                          (IDM,NYY).    
C    
C                        FYY    
C                          ARRAY OF SECOND DERIVATIVES OF F WITH RESPECT
C                          TO Y COMPUTED BY COEFF2.  DIMENSION OF FYY IN
C                          THE CALLING PROGRAM IS (IDM,NYY).    
C    
C                        FXXYY    
C                          ARRAY OF FOURTH DERIVATIVES    
C                          (D/DX)**2*(D/DY)**2*F, COMPUTED BY COEFF2.   
C                          DIMENSION OF FXXYY IN THE CALLING PROGRAM IS 
C                          (IDM,NYY).    
C    
C                        THE ARRAYS X, Y, F, FXX, FYY, FXXYY ARE USED AS
C                        INPUT FOR THE ROUTINE TERP2 WHICH PERFORMS    
C                        INTERPOLATION AT REQUIRED VALUES OF THE TWO    
C                        INDEPENDENT VARIABLES.    
C    
C SPECIAL CONDITIONS     NONE    
C    
C COMMON BLOCKS          NONE    
C    
C I/O                    NONE    
C    
C PRECISION              SINGLE    
C    
C SPECIALIST             CICELY RIDLEY, NCAR, BOULDER, COLORADO 80307   
C    
C LANGUAGE               FORTRAN    
C    
C HISTORY                THIS ROUTINE IS A STREAMLINED AND    
C                        STANDARDIZED VERSION OF A SIMILAR ROUTINE    
C                        FROM THE LOS ALAMOS CUBIC SPLINE PACKAGE.    
C    
C SPACE REQUIRED         370 (OCTAL) = 248 (DECIMAL)    
C    
C TIMING                 THE TIMING IS PROPORTIONAL TO NX*NY.  FOR A    
C                        21 X 21 GRID OF DATA POINTS, THE TIME TAKEN ON 
C                        THE NCAR 7600 WAS 19.5 MILLISECONDS.    
C    
C PORTABILITY            FULLY PORTABLE    
C    
C REQUIRED RESIDENT      NONE    
C ROUTINES    
C    
C    
C    
C    
C    
C    
C    
C    
C    
C    
C    
      DIMENSION       X(*)       ,Y(*)       ,F(IDM,*)   ,FXX(IDM,*) ,  
     1                FYY(IDM,*) ,FXXYY(IDM,*)           ,IBD(4)     ,  
     2                ILOC(2)    ,JLOC(2)    ,WK(*)
      DATA ILOC(1),ILOC(2),JLOC(1),JLOC(2)/1,1,4,4/    
C THE FOLLOWING CALL IS FOR GATHERING STATISTICS ON LIBRARY USE AT NCAR 
C.....CALL Q8QST4( 4HNSSL      , 6HCUBSPL    ,6HCOEFF2 ,10HVERSION 11)  
C    
C COMPUTE FXX    
C    
      DO 101 J=1,NY    
         CALL COEFF1 (NX,X,F(1,J),FXX(1,J),IBD(1),1,WK)    
  101 CONTINUE    
C    
C COMPUTE FYY    
C    
      DO 102 I=1,NX    
         CALL COEFF1 (NY,Y,F(I,1),FYY(I,1),IBD(3),IDM,WK)    
  102 CONTINUE    
C    
C CHECK FOR PERIODIC BOUNDARY CONDITION IN BOTH DIRECTIONS    
C    
      IF (IBD(1) .EQ. 3) GO TO 103    
      IF (IBD(3) .EQ. 3) GO TO 105    
C    
C CALCULATE FXXYY ALONG LEFT AND RIGHT BOUNDARIES    
C    
      CALL COEFF1 (NY,Y,FXX(1,1),FXXYY(1,1),JLOC,IDM,WK)    
      CALL COEFF1 (NY,Y,FXX(NX,1),FXXYY(NX,1),JLOC,IDM,WK)    
      GO TO 106    
  103 CONTINUE    
C    
C PERIODIC IN X DIRECTION . CALCULATE FXXYY ALONG LOWER AND UPPER    
C BOUNDARIES.    
C    
      CALL COEFF1 (NX,X,FYY(1,1),FXXYY(1,1),IBD(1),1,WK)    
      CALL COEFF1 (NX,X,FYY(1,NY),FXXYY(1,NY),IBD(1),1,WK)    
C    
C CALCULATE REMAINING FXXYY    
C    
      DO 104 I=1,NX    
         CALL COEFF1 (NY,Y,FXX(I,1),FXXYY(I,1),ILOC,IDM,WK)    
  104 CONTINUE    
      GO TO 108    
  105 CONTINUE    
C    
C PERIODIC IN Y DIRECTION. CALCULATE FXXYY ALONG LEFT AND RIGHT    
C BOUNDARIES.    
C    
      CALL COEFF1 (NY,Y,FXX(1,1),FXXYY(1,1),IBD(3),IDM,WK)    
      CALL COEFF1 (NY,Y,FXX(NX,1),FXXYY(NX,1),IBD(3),IDM,WK)    
  106 CONTINUE    
C    
C CALCULATE REMAINING FXXYY    
C    
      DO 107 J=1,NY    
         CALL COEFF1 (NX,X,FYY(1,J),FXXYY(1,J),ILOC,1,WK)    
  107 CONTINUE    
  108 CONTINUE    
      RETURN    
      END    
      FUNCTION TERP2 (XB,YB,NX,X,NY,Y,F,FXX,FYY,FXXYY,IDM,IXD,IYD)    
      IMPLICIT REAL*4 (A-H,O-Z)
C    
C    
C DIMENSION OF           X(NX),Y(NY),F(IDM,NY),FXX(IDM,NY),FYY(IDM,NY), 
C ARGUMENTS              FXXYY(IDM,NY))    
C                        (IDM MUST BE .GE. NX)    
C    
C LATEST REVISION        JANUARY 1978    
C    
C PURPOSE                USING THE COEFFICIENTS PRODUCED BY COEFF2,    
C                        THIS ROUTINE EVALUATES THE FUNCTION ON A    
C                        SELECTED DERIVATIVE OF ANY POINT WHERE    
C                        TWO-DIMENSIONAL INTERPOLATION IS REQUIRED.    
C    
C ACCESS CARD            *FORTRAN,S=ULIB,N=CUBSPL    
C    
C USAGE                  R = TERP2 (XB,YB,NX,X,NY,Y,F,FXX,FYY,FXXYY,IDM,
C                                   IXD,IYD)    
C    
C ARGUMENTS    
C    
C ON INPUT               XB, YB    
C                          VALUES OF THE INDEPENDENT VARIABLES, X AND Y,
C                          AT WHICH INTERPOLATION IS REQUIRED.    
C    
C                        NX    
C                          NUMBER OF GRID POINTS IN THE X-DIRECTION.  NX
C                          MUST BE AT LEAST 4.    
C    
C                        X    
C                          TABLE OF NX VALUES OF THE INDEPENDENT    
C                          VARIABLE, X, ARRANGED IN ASCENDING ORDER.    
C                          DIMENSION OF X IN THE CALLING PROGRAM MUST BE
C                          AT LEAST NX.    
C    
C                        NY    
C                          NUMBER OF GRID POINTS IN THE Y-DIRECTION.  NY
C                          MUST BE AT LEAST 4.    
C    
C                        Y    
C                          TABLE OF NY VALUES OF THE INDEPENDENT    
C                          VARIABLE, Y, ARRANGED IN ASCENDING ORDER.    
C                          DIMENSION OF Y IN THE CALLING PROGRAM MUST BE
C                          AT LEAST NY.    
C    
C                        F    
C                          TWO-DIMENSIONAL ARRAY OF FUNCTION VALUES AT  
C                          GRID POINTS DEFINED BY THE ARRAYS X AND Y.   
C                          DIMENSION OF F IN THE CALLING PROGRAM IS    
C                          (IDM,NYY), WHERE    
C                              IDM .GE. NX    
C                              NYY .GE. NY    
C    
C                        FXX    
C                          ARRAY OF SECOND DERIVATIVES OF F WITH RESPECT
C                          TO X COMPUTED BY COEFF2.  DIMENSION OF FXX IN
C                          THE CALLING PROGRAM IS (IDM,NYY).  SEE UNDER 
C                          F ABOVE.    
C    
C                        FYY    
C                          ARRAY OF SECOND DERIVATIVES OF F WITH RESPECT
C                          TO Y COMPUTED BY COEFF2.  DIMENSION OF FYY IN
C                          THE CALLING PROGRAM IS (IDM,NYY).    
C    
C                        FXXYY    
C                          ARRAY OF FOURTH DERIVATIVES,    
C                          (D/DX)**2*(D/DY)**2*F, COMPUTED BY COEFF2.   
C                          DIMENSION OF FXXYY IN THE CALLING PROGRAM IS 
C                          (IDM,NYY).    
C    
C                        IDM    
C                          FIRST DIMENSION IN THE CALLING PROGRAM OF    
C                          ARRAYS F, FXX, FYY AND FXXYY,    
C                              IDM .GE. NX    
C    
C                        IXD, IYD    
C                          DEFINE DERIVATIVE TO BE RETURNED BY THE    
C                          FUNCTION TERP2.  IXD, IYD MAY EACH TAKE THE  
C                          THE VALUES 0, 1, 2.  THE DERIVATIVE RETURNED 
C                          IS (D/DX)**IXD*(D/DY)**IYD*F.    
C                            NOTE THAT IF IXD = IYD = 0, THE FUNCTION   
C                            VALUE ITSELF IS RETURNED.    
C    
C SPECIAL CONDITIONS     NONE    
C    
C COMMON BLOCKS          NONE    
C    
C I/O                    NONE    
C    
C PRECISION              SINGLE    
C    
C SPECIALIST             CICELY RIDLEY, NCAR, BOULDER, COLORADO 80307   
C    
C LANGUAGE               FORTRAN    
C    
C HISTORY                THIS ROUTINE IS A STREAMLINED AND    
C                        STANDARDIZED VERSION OF A SIMILAR ROUTINE    
C                        FROM THE LOS ALAMOS CUBIC SPLINE PACKAGE.    
C    
C SPACE REQUIRED         220 (OCTAL) = 144 (DECIMAL)    
C    
C TIMING                 THIS PROCEDURE IS FAST.  THE MAXIMUM TIME FOR  
C                        THE BINARY SEARCH IS PROPORTIONAL TO    
C                        ALOG(NX*NY).  THE TIME FOR FUNCTION EVALUATION 
C                        IS INDEPENDENT OF N.  FOR A 21 X 21 GRID OF    
C                        DATA POINTS, AN AVERAGE TIME FOR AN    
C                        INTERPOLATION ON THE NCAR CDC 7600 IS ABOUT    
C                        .29 MILLISECONDS.    
C    
C PORTABILITY            FULLY PORTABLE    
C    
C REQUIRED RESIDENT      ALOG    
C ROUTINES    
C    
C    
C    
C    
C    
C    
C    
C    
C    
C    
C    
      DIMENSION       X(*)       ,Y(*)       ,F(IDM,*)   ,FXX(IDM,*) ,  
     1                FYY(IDM,*) ,FXXYY(IDM,*)           ,FF(2)      ,  
     2                WW(2)      ,TAB(3)     ,ITAB(3)    
C THE FOLLOWING CALL IS FOR GATHERING STATISTICS ON LIBRARY USE AT NCAR 
C.....CALL Q8QST4( 4HNSSL      , 6HCUBSPL    ,5HTERP2  ,10HVERSION 11)  
C    
C SEARCH IN X AND Y ARRAYS.    
C    
      CALL SEARCH (XB,X,NX,I)    
      CALL SEARCH (YB,Y,NY,J)    
C    
C INTERPOLATE IN X DIRECTION    
C    
      DO 101 I1=1,3    
         ITAB(I1) = 0    
  101 CONTINUE    
      I1 = IXD+1    
      ITAB(I1) = 1    
      DO 102 J1=1,2    
         JJ = J+J1-1    
         CALL INTERP (N,X,F(1,JJ),FXX(1,JJ),XB,I,1,TAB,ITAB)    
         FF(J1) = TAB(I1)    
         CALL INTERP (N,X,FYY(1,JJ),FXXYY(1,JJ),XB,I,1,TAB,ITAB)    
         WW(J1) = TAB(I1)    
  102 CONTINUE    
C    
C INTERPOLATE IN Y DIRECTION    
C    
      DO 103 J1=1,3    
         ITAB(J1) = 0    
  103 CONTINUE    
      J1 = IYD+1    
      ITAB(J1) = 1    
      CALL INTERP (2,Y(J),FF,WW,YB,1,1,TAB,ITAB)    
      TERP2 = TAB(J1)    
      RETURN    
C    
C REVISION HISTORY---    
C    
C JUNE 1977        REPLACED NON-STANDARD STATEMENT FUNCTIONS AND    
C                  SUBSCRIPTS TO ENHANCE PORTABILITY.    
C    
C JANUARY 1978     DELETED REFERENCES TO THE  *COSY  CARDS, MOVED    
C                  THE REVISION HISTORIES TO APPEAR BEFORE THE    
C                  FINAL END CARD, AND MOVED THE INITIAL COMMENT    
C                  CARDS TO APPEAR AFTER THE FIRST SUBROUTINE CARD    
C                  AND CHANGED  ITAB  FROM LOGICAL TO INTEGER IN    
C                  SUBROUTINE INTERP AND CORRECTED PROBLEM WITH    
C                  VERSION NUMBERS IN ONE STATISTICS CALL    
C-----------------------------------------------------------------------
      END    
C    
      subroutine landfill1(a1,mask,m,n, npass, ldebug)
      implicit none
c
      logical ldebug
      integer m,n,npass
      integer mask(m,n)
      real*4  a1(m,n)
c
c --- creeping extrapolation into the land mask,
c --- using npass's of a 9-point smoother based extrapolation scheme.
c --- mask == 1 for ocean.
c
      integer, allocatable :: mm(:,:,:),mmsum(:)
c
      integer i,ii,ip0,ip1,ipass,j,jj,ki,kj,nup
      real*4  sa,ss
c
      real*4  s(-1:1,-1:1)
      data s / 1.0,2.0,1.0, 2.0,4.0,2.0, 1.0,2.0,1.0 /
c
      allocate( mm(0:m+1,0:n+1,0:1) )
      allocate( mmsum(   1:n      ) )
c
      mm( : ,0,  0) = 0
      mm( : ,n+1,0) = 0
      mm(1:m,1:n,0) = mask(:,:)
      do j= 1,n
        mmsum(j) = m - sum( mm(1:m,j,0) )
      enddo
c
      do ipass= 1,npass
        ip0 = mod(ipass+1,2)
        ip1 = mod(ipass,  2)
        mm(0,  :,ip0) = mm(m,:,ip0)  !assume periodic
        mm(m+1,:,ip0) = mm(1,:,ip0)  !assume periodic
        mm(:,  :,ip1) = mm(:,:,ip0)
        nup = 0
        do j= 1,n
          if     (mmsum(j).ne.0) then  !if land left in this row
            do i= 1,m
              if     (mm(i,j,ip0).eq.0) then
                sa = 0.0
                ss = 0.0
                do kj= -1,1
                  jj = j+kj
                  do ki= -1,1
                    ii = i+ki
                    if     (ii.lt.1) then !assume periodic
                      ii = m
                    elseif (ii.gt.m) then !assume periodic
                      ii = 1
                    endif
                    if     (mm(ii,jj,ip0).eq.1) then
                      sa = sa + s(ki,kj)*a1(ii,jj)
                      ss = ss + s(ki,kj)
                    endif
                  enddo !ki
                enddo !kj
                if     (ss.gt.1.5) then  !not just a single diagonal 
                  a1(i,j)     = sa/ss
                  mm(i,j,ip1) = 1
                  nup         = nup + 1
                  mmsum(j)    = mmsum(j) - 1
                  if     (mmsum(j).eq.0) then
                    exit  !i-loop
                  endif
*                 if     (ldebug .and. mod(nup,1000).eq.1) then
*                   write(6,'(a,2i5,f5.1,f10.3)') 
*    &                '   i,j,ss,a = ',i,j,ss,a1(i,j)
*                 endif
                endif !ss>1.5
              endif !mm.eq.0
            enddo !i
          endif !mmsum.ne.0
        enddo !j
        if     (ldebug) then
          write(6,'(a,i4,a,i6,a)') 'landfill1: pass',ipass,
     &                             ' filled in',nup,' points'
        endif
        if     (nup.eq.0) then
          exit
        endif
      enddo  ! ipass=1,npass
      if     (ldebug) then
        write(6,*)
      endif
c
      deallocate( mm, mmsum )
c
      return
      end
C    
      subroutine landfill2(a1,a2,mask,m,n, npass, ldebug)
      implicit none
c
      logical ldebug
      integer m,n,npass
      integer mask(m,n)
      real*4  a1(m,n),a2(m,n)
c
c --- creeping extrapolation into the land mask,
c --- using npass's of a 9-point smoother based extrapolation scheme.
c --- mask == 1 for ocean.
c
      integer, allocatable :: mm(:,:,:),mmsum(:)
c
      integer i,ii,ip0,ip1,ipass,j,jj,ki,kj,nup
      real*4  sa(2),ss
c
      real*4  s(-1:1,-1:1)
      data s / 1.0,2.0,1.0, 2.0,4.0,2.0, 1.0,2.0,1.0 /
c
      allocate( mm(0:m+1,0:n+1,0:1) )
      allocate( mmsum(   1:n      ) )
c
      mm( : ,0,  0) = 0
      mm( : ,n+1,0) = 0
      mm(1:m,1:n,0) = mask(:,:)
      do j= 1,n
        mmsum(j) = m - sum( mm(1:m,j,0) )
      enddo
c
      do ipass= 1,npass
        ip0 = mod(ipass+1,2)
        ip1 = mod(ipass,  2)
        mm(0,  :,ip0) = mm(m,:,ip0)  !assume periodic
        mm(m+1,:,ip0) = mm(1,:,ip0)  !assume periodic
        mm(:,  :,ip1) = mm(:,:,ip0)
        nup = 0
        do j= 1,n
          if     (mmsum(j).ne.0) then  !if land left in this row
            do i= 1,m
              if     (mm(i,j,ip0).eq.0) then
                sa(1) = 0.0
                sa(2) = 0.0
                ss    = 0.0
                do kj= -1,1
                  jj = j+kj
                  do ki= -1,1
                    ii = i+ki
                    if     (ii.lt.1) then !assume periodic
                      ii = m
                    elseif (ii.gt.m) then !assume periodic
                      ii = 1
                    endif
                    if     (mm(ii,jj,ip0).eq.1) then
                      sa(1) = sa(1) + s(ki,kj)*a1(ii,jj)
                      sa(2) = sa(2) + s(ki,kj)*a2(ii,jj)
                      ss    = ss    + s(ki,kj)
                    endif
                  enddo !ki
                enddo !kj
                if     (ss.gt.1.5) then  !not just a single diagonal 
                  a1(i,j)     = sa(1)/ss
                  a2(i,j)     = sa(2)/ss
                  mm(i,j,ip1) = 1
                  nup         = nup + 1
                  mmsum(j)    = mmsum(j) - 1
                  if     (mmsum(j).eq.0) then
                    exit  !i-loop
                  endif
*                 if     (ldebug .and. mod(nup,1000).eq.1) then
*                   write(6,'(a,2i5,f5.1,f10.3)') 
*    &                '   i,j,ss,a = ',i,j,ss,a1(i,j)
*                 endif
                endif !ss>1.5
              endif !mm.eq.0
            enddo !i
          endif !mmsum.ne.0
        enddo !j
        if     (ldebug) then
          write(6,'(a,i4,a,i6,a)') 'landfill2: pass',ipass,
     &                             ' filled in',nup,' points'
        endif
        if     (nup.eq.0) then
          exit
        endif
      enddo  ! ipass=1,npass
      if     (ldebug) then
        write(6,*)
      endif
c
      deallocate( mm, mmsum )
c
      return
      end
C    
      subroutine landfill5(a1,a2,a3,a4,a5,mask,m,n, npass, ldebug)
      implicit none
c
      logical ldebug
      integer m,n,npass
      integer mask(m,n)
      real*4  a1(m,n),a2(m,n),a3(m,n),a4(m,n),a5(m,n)
c
c --- creeping extrapolation into the land mask,
c --- using npass's of a 9-point smoother based extrapolation scheme.
c --- mask == 1 for ocean.
c
      integer, allocatable :: mm(:,:,:),mmsum(:)
c
      integer i,ii,ip0,ip1,ipass,j,jj,ki,kj,nup
      real*4  sa(5),ss
c
      real*4  s(-1:1,-1:1)
      data s / 1.0,2.0,1.0, 2.0,4.0,2.0, 1.0,2.0,1.0 /
c
      allocate( mm(0:m+1,0:n+1,0:1) )
      allocate( mmsum(   1:n      ) )
c
      mm( : ,0,  0) = 0
      mm( : ,n+1,0) = 0
      mm(1:m,1:n,0) = mask(:,:)
      do j= 1,n
        mmsum(j) = m - sum( mm(1:m,j,0) )
      enddo
c
      do ipass= 1,npass
        ip0 = mod(ipass+1,2)
        ip1 = mod(ipass,  2)
        mm(0,  :,ip0) = mm(m,:,ip0)  !assume periodic
        mm(m+1,:,ip0) = mm(1,:,ip0)  !assume periodic
        mm(:,  :,ip1) = mm(:,:,ip0)
        nup = 0
        do j= 1,n
          if     (mmsum(j).ne.0) then  !if land left in this row
            do i= 1,m
              if     (mm(i,j,ip0).eq.0) then
                sa(1) = 0.0
                sa(2) = 0.0
                sa(3) = 0.0
                sa(4) = 0.0
                sa(5) = 0.0
                ss    = 0.0
                do kj= -1,1
                  jj = j+kj
                  do ki= -1,1
                    ii = i+ki
                    if     (ii.lt.1) then !assume periodic
                      ii = m
                    elseif (ii.gt.m) then !assume periodic
                      ii = 1
                    endif
                    if     (mm(ii,jj,ip0).eq.1) then
                      sa(1) = sa(1) + s(ki,kj)*a1(ii,jj)
                      sa(2) = sa(2) + s(ki,kj)*a2(ii,jj)
                      sa(3) = sa(3) + s(ki,kj)*a3(ii,jj)
                      sa(4) = sa(4) + s(ki,kj)*a4(ii,jj)
                      sa(5) = sa(5) + s(ki,kj)*a5(ii,jj)
                      ss    = ss    + s(ki,kj)
                    endif
                  enddo !ki
                enddo !kj
                if     (ss.gt.1.5) then  !not just a single diagonal 
                  a1(i,j)     = sa(1)/ss
                  a2(i,j)     = sa(2)/ss
                  a3(i,j)     = sa(3)/ss
                  a4(i,j)     = sa(4)/ss
                  a5(i,j)     = sa(5)/ss
                  mm(i,j,ip1) = 1
                  nup         = nup + 1
                  mmsum(j)    = mmsum(j) - 1
                  if     (mmsum(j).eq.0) then
                    exit  !i-loop
                  endif
*                 if     (ldebug .and. mod(nup,1000).eq.1) then
*                   write(6,'(a,2i5,f5.1,f10.3)') 
*    &                '   i,j,ss,a = ',i,j,ss,a1(i,j)
*                 endif
                endif !ss>1.5
              endif !mm.eq.0
            enddo !i
          endif !mmsum.ne.0
        enddo !j
        if     (ldebug) then
          write(6,'(a,i4,a,i6,a)') 'landfill5: pass',ipass,
     &                             ' filled in',nup,' points'
        endif
        if     (nup.eq.0) then
          exit
        endif
      enddo  ! ipass=1,npass
      if     (ldebug) then
        write(6,*)
      endif
c
      deallocate( mm, mmsum )
c
      return
      end
C
      subroutine smooth1(a1,mask,m,n, mask_smooth,npass, ldebug)
      implicit none
c
      logical ldebug
      integer m,n,mask_smooth,npass
      integer mask(m,n)
      real*4  a1(m,n)
c
c --- smooth npass times where mask == mask_smooth.
c
      integer, allocatable :: mm(:,:)
      real*4,  allocatable :: sm(:,:,:)
c
      integer i,ii,ip0,ip1,ipass,j,jj,ki,kj,npts
      real*4  sa,qss
c
      real*4  s(-1:1,-1:1)
      data s / 1.0,2.0,1.0, 2.0,4.0,2.0, 1.0,2.0,1.0 /
c
      qss = 1.0/16.0 !1.0/sum(s(:,:))
c
      allocate( mm(1:m,0:n+1)     )
      allocate( sm(1:m,1:n,  0:1) )
c
      mm( : ,0  ) = mask_smooth + 1  !not smoothed
      mm( : ,n+1) = mask_smooth + 1  !not smoothed
      do j= 1,n
        do i= 1,m
          mm(i,j)   = mask(i,j)
          sm(i,j,0) =   a1(i,j)
          sm(i,j,1) =   a1(i,j)
        enddo
      enddo
c
      do ipass= 1,npass
        ip0 = mod(ipass+1,2) !last pass
        ip1 = mod(ipass,  2) !this pass
        npts = 0
        do j= 1,n
          do i= 1,m
            if     (mm(i,j).eq.mask_smooth) then !smooth
              sa = 0.0
              do kj= -1,1
                jj = j+kj
                do ki= -1,1
                  ii = i+ki
                  if     (ii.lt.1) then !assume periodic
                    ii = m
                  elseif (ii.gt.m) then !assume periodic
                    ii = 1
                  endif
                  if     (mm(ii,jj).eq.mask_smooth) then
                    sa = sa + s(ki,kj)*sm(ii,jj,ip0)
                  else
                    sa = sa + s(ki,kj)*sm(i, j ,ip0)
                  endif
                enddo !ki
              enddo !kj
              sm(i,j,ip1) = sa*qss
              npts = npts + 1
            endif !mm.eq.mask_smooth
          enddo !i
        enddo !j
        if     (ldebug) then
          write(6,'(a,2i12)') 'smooth1 - ipass,npts =',ipass,npts
        endif
      enddo  ! ipass=1,npass
*     write(6,'(a,2i12)') 'smooth1 - ipass,npts =',npass,npts
c
      do j= 1,n
        do i= 1,m
          a1(i,j) = sm(i,j,ip1)
        enddo
      enddo
c
      deallocate( mm, sm )
c
      return
      end
