      program sigma2s_to_locsig
      implicit none
c
c     usage:  echo pot.density saln depth | sigma2s_to_locsig
c
c     input:  p.density saln depth (sigma2 psu m)
c     output: p.density saln p.temp depth density
c
c             p.density  = potential density (sigma2)
c             saln       = salinity (psu)
c             p.temp     = potential temperature, 17-term (degC)
c             depth      = depth (m)
c             density    = in-situ density (for p.T,S,depth)
c
c     potential temperature is always w.r.t the surface
c
c     alan j. wallcraft, COAPS, January 2018.
c
      integer ios
      real*8  temp,saln,depth,pden,dens
      real*8  tofsig_6,sigloc_6
c
      do
        read(5,*,iostat=ios) pden,saln,depth
        if     (ios.ne.0) then
          exit
        endif
        temp = tofsig_6(pden,saln)
        dens = sigloc_6(temp,saln,depth*1.D4)
        write(6,'(5f10.3)') pden,saln,temp,depth,dens
      enddo
      end

      REAL*8 FUNCTION SIGLOC_6(TT8,SS8,PRS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8,PRS8
      INCLUDE '../include/stmt_fns_SIGMA2_17term.h'
      SIGLOC_6 = SIGLOC(TT8,SS8,PRS8)
      END
      REAL*8 FUNCTION TOFSIG_6(RR8,SS8)
      IMPLICIT NONE
      REAL*8  RR8,SS8
      REAL*8, PARAMETER :: TOL=1.D-6
      INTEGER NN
      REAL*8  TN,TO
      REAL*8  TOFSIG_8
      INCLUDE '../include/stmt_fns_SIGMA2_17term.h'
C     sofsig via Newton iteration from a 12-term 1st guess
      TN = TOFSIG_8(RR8,SS8)  !non-negative
      DO NN= 1,10
        TO = TN
        TN = TO - (SIG(TO,SS8)-RR8)/DSIGDT(TO,SS8)
        IF     (NN.EQ.10 .OR. ABS(TN-TO).LT.TOL) THEN
          EXIT  
        ENDIF 
      ENDDO !nn
      TOFSIG_6 = TN
      END
      REAL*8 FUNCTION TOFSIG_8(RR8,SS8)
      IMPLICIT NONE
      REAL*8  RR8,SS8
      INCLUDE '../include/stmt_fns_SIGMA2_12term.h'
      TOFSIG_8 = TOFSIG(RR8,SS8)
      END
