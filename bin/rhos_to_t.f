      program rhos_to_t
      implicit none
c
c     usage:  echo pot.density saln sigtyp | rhos_to_t
c
c     input:  p.density saln {0,2}
c     output: p.density saln p.temp_7t p.temp_9t p.temp_17t p.temp_12t
c
c             p.density  = potential density (sigma0 or sigma2)
c             saln       = Salinity (psu)
c             p.temp_7t  = potential temperature,  7-term
c             p.temp_9t  = potential temperature,  9-term
c             p.temp_17t = potential temperature, 17-term
c             p.temp_12t = potential temperature, 12-term
c
c     potential temperature is always w.r.t the surface
c
c     alan j. wallcraft, naval research laboratory, september 2013.
c
      real*8  temp(4),saln,dens
      integer ios,n,sigtyp,st_old
      real*8  tofsig_v
c
      st_old = -1
      do
        read(5,*,iostat=ios) dens,saln,sigtyp
        if     (ios.ne.0) then
          exit
        endif
        if     (sigtyp.ne.st_old) then
          write(6,'(2a)') '#    p.dens   saln   p.temp_7t  p.temp_9t',
     &                                        ' p.temp_17t p.temp_12t'
          st_old = sigtyp
        endif
        if     (sigtyp.eq.0) then
          do n= 1,4
            temp(n) = tofsig_v(dens,saln,2*n-1)
          enddo
        else
          do n= 1,4
            temp(n) = tofsig_v(dens,saln,2*n)
          enddo
        endif
        write(6,'(6f10.3)') dens,saln,temp(1:4)
      enddo
      end

      REAL*8 FUNCTION TOFSIG_V(RR,SS,SIGVER)
      IMPLICIT NONE
      INTEGER SIGVER
      REAL*8  RR,SS
C
C     SIGVER WRAPPER FOR TOFSIG
C
      REAL*8 RR8,SS8
      REAL*8 TOFSIG_1,TOFSIG_2,TOFSIG_3,TOFSIG_4,
     &       TOFSIG_5,TOFSIG_6,TOFSIG_7,TOFSIG_8,
     &       TOFSIG_46,TOFSIG_48
C
      RR8 = RR
      SS8 = SS
      IF     (SIGVER.GT.40) THEN
        IF     (SIGVER.EQ.46) THEN
          TOFSIG_V = TOFSIG_46(RR8,SS8)
        ELSEIF (SIGVER.EQ.48) THEN
          TOFSIG_V = TOFSIG_48(RR8,SS8)
        ENDIF
      ELSEIF (MOD(SIGVER,2).EQ.1) THEN
        IF     (SIGVER.EQ.1) THEN
          TOFSIG_V = TOFSIG_1(RR8,SS8)
        ELSEIF (SIGVER.EQ.3) THEN
          TOFSIG_V = TOFSIG_3(RR8,SS8)
        ELSEIF (SIGVER.EQ.5) THEN
          TOFSIG_V = TOFSIG_5(RR8,SS8)
        ELSEIF (SIGVER.EQ.7) THEN
          TOFSIG_V = TOFSIG_7(RR8,SS8)
        ENDIF
      ELSE
        IF     (SIGVER.EQ.2) THEN
          TOFSIG_V = TOFSIG_2(RR8,SS8)
        ELSEIF (SIGVER.EQ.4) THEN
          TOFSIG_V = TOFSIG_4(RR8,SS8)
        ELSEIF (SIGVER.EQ.6) THEN
          TOFSIG_V = TOFSIG_6(RR8,SS8)
        ELSEIF (SIGVER.EQ.8) THEN
          TOFSIG_V = TOFSIG_8(RR8,SS8)
        ENDIF
      ENDIF
      RETURN
      END
      REAL*8 FUNCTION TOFSIG_1(RR8,SS8)
      IMPLICIT NONE
      REAL*8  RR8,SS8
      INCLUDE '../include/stmt_fns_SIGMA0_7term.h'
      TOFSIG_1 = TOFSIG(RR8,SS8)
      END
      REAL*8 FUNCTION TOFSIG_3(RR8,SS8)
      IMPLICIT NONE
      REAL*8  RR8,SS8
      INCLUDE '../include/stmt_fns_SIGMA0_9term.h'
      TOFSIG_3 = TOFSIG(RR8,SS8)
      END
      REAL*8 FUNCTION TOFSIG_5(RR8,SS8)
      IMPLICIT NONE
      REAL*8  RR8,SS8
      REAL*8, PARAMETER :: TOL=1.D-6
      INTEGER NN
      REAL*8  TN,TO
      REAL*8  TOFSIG_7
      INCLUDE '../include/stmt_fns_SIGMA0_17term.h'
C     sofsig via Newton iteration from a 12-term 1st guess
      TN = TOFSIG_7(RR8,SS8)  !non-negative
      DO NN= 1,10
        TO = TN
        TN = TO - (SIG(TO,SS8)-RR8)/DSIGDT(TO,SS8)
        IF     (NN.EQ.10 .OR. ABS(TN-TO).LT.TOL) THEN
          EXIT  
        ENDIF 
      ENDDO !nn
      TOFSIG_5 = TN
      END
      REAL*8 FUNCTION TOFSIG_7(RR8,SS8)
      IMPLICIT NONE
      REAL*8  RR8,SS8
      INCLUDE '../include/stmt_fns_SIGMA0_12term.h'
      TOFSIG_7 = TOFSIG(RR8,SS8)
      END
      REAL*8 FUNCTION TOFSIG_2(RR8,SS8)
      IMPLICIT NONE
      REAL*8  RR8,SS8
      INCLUDE '../include/stmt_fns_SIGMA2_7term.h'
      TOFSIG_2 = TOFSIG(RR8,SS8)
      END
      REAL*8 FUNCTION TOFSIG_4(RR8,SS8)
      IMPLICIT NONE
      REAL*8  RR8,SS8
      INCLUDE '../include/stmt_fns_SIGMA2_9term.h'
      TOFSIG_4 = TOFSIG(RR8,SS8)
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
      REAL*8 FUNCTION TOFSIG_46(RR8,SS8)
      IMPLICIT NONE
      REAL*8  RR8,SS8
      REAL*8, PARAMETER :: TOL=1.D-6
      INTEGER NN
      REAL*8  TN,TO
      REAL*8  TOFSIG_48
      INCLUDE '../include/stmt_fns_SIGMA4_17term.h'
C     sofsig via Newton iteration from a 12-term 1st guess
      TN = TOFSIG_48(RR8,SS8)  !non-negative
      DO NN= 1,10
        TO = TN
        TN = TO - (SIG(TO,SS8)-RR8)/DSIGDT(TO,SS8)
        IF     (NN.EQ.10 .OR. ABS(TN-TO).LT.TOL) THEN
          EXIT  
        ENDIF 
      ENDDO !nn
      TOFSIG_46 = TN
      END
      REAL*8 FUNCTION TOFSIG_48(RR8,SS8)
      IMPLICIT NONE
      REAL*8  RR8,SS8
      INCLUDE '../include/stmt_fns_SIGMA4_12term.h'
      TOFSIG_48 = TOFSIG(RR8,SS8)
      END
