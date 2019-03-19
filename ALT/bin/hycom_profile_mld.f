      PROGRAM HYCOM_PROFILE_MLD
      IMPLICIT NONE
C
C  hycom_profile_mld - Usage: hycom_profile_mld archv.txt tmljmp [itype]
C
C                 calculates MLD from a HYCOM isopycnal text profile file
C
C   archv.txt is assumed to be an HYCOM archive text profile file
C   tmljmp    is the equivalent temperature jump across mixed-layer (degC)
C   itype     is the density profile method (default 1)
C                =-1; interpolation between cell centers
C                = 0; picewise  constant  method (PCM)
C                = 1; piecewise linear    method (PLM)
C                = 2; piecewise linear    method (PLM), locsig
C                = 3; piecewise parabolic method (PQM), locsig
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  September 2002.
C
      INTEGER      IARGC
      INTEGER      NARG
      CHARACTER*240 CARG
C
      CHARACTER*240 CFILEA
      CHARACTER*240 CLINE
      REAL          DPMIXL,RMLJMP,TMLJMP,THK(99),FLAG
      REAL          U,V,T(99),S(99),R(99),P(99+1), R0(99)
      INTEGER       I,IOS,ITYPE,K,KDM
      REAL          DUM,DVM,DTM,DSM,DRM
C
      INTEGER       YRFLAG,IYR,IDY,IHR
      REAL          ALON,ALAT,DAYM,YEARM
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.3) THEN
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CLINE)
        READ(CLINE,*) TMLJMP
        CALL GETARG(3,CLINE)
        READ(CLINE,*) ITYPE
      ELSEIF (NARG.EQ.2) THEN
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CLINE)
        READ(CLINE,*) TMLJMP
        ITYPE = 1
      ELSE
        WRITE(6,*)
     +    'Usage: hycom_profile_mld archv.txt tmljmp [itype]'
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
C
C     READ THE ISOPYCNAL PROFILE.
C
      READ( 11,'(a)') CLINE
      READ( 11,'(2x,6i7,2f8.1,4i7)') 
     +  I,I,I,I,I,I,ALON,ALAT,YRFLAG,IYR,IDY,IHR
      READ( 11,'(a)') CLINE
      READ( 11,'(2x,f11.2)')
     +  DAYM
      READ( 11,'(a)') CLINE
      IF     (YRFLAG.EQ.0) THEN
        YEARM = IYR + (IDY-1.0+IHR/24.0)/360.0
      ELSEIF (YRFLAG.LE.2) THEN
        YEARM = IYR + (IDY-1.0+IHR/24.0)/366.0
      ELSE
        IF     (MOD(IYR,4).NE.0) THEN
          YEARM = IYR + (IDY-1.0+IHR/24.0)/365.0
        ELSE
          YEARM = IYR + (IDY-1.0+IHR/24.0)/366.0
        ENDIF
      ENDIF
C
      P(1) =  0.0
      KDM  = -1
      DO K= 1,99
        DO !extra header record loop
          READ(11,'(a)',IOSTAT=IOS) CLINE
          IF     (IOS.NE.0 .OR. CLINE(1:1).NE."#") THEN
            EXIT
          ENDIF
        ENDDO !extra header record loop
        IF     (IOS.NE.0) THEN
          IF     (K.NE.KDM+1) THEN
            WRITE(6,*) 'Error: inconsistent input profile'
            CALL EXIT(6)
          ENDIF
          EXIT
        ENDIF
        READ(CLINE,*) KDM,U,V,T(K),S(K),R(K),THK(K)
        P(K+1) = P(K) + THK(K)
        IF     (THK(K).EQ.0.0) THEN
          T(K) = T(K-1)
          S(K) = S(K-1)
          R(K) = R(K-1)
        ENDIF
      ENDDO
      CLOSE(11)
C
C     FIND THE MIXED LAYER DEPTH.
C
      IF     (ITYPE.GE.2) THEN
        CALL RMLJMPL(RMLJMP, T(1),S(1), 1,1, TMLJMP)
      ELSEIF (R(KDM).LT.32.0) THEN
        CALL RMLJMP0(RMLJMP, T(1),S(1), 1,1, TMLJMP)
      ELSE
        CALL RMLJMP2(RMLJMP, T(1),S(1), 1,1, TMLJMP)
      ENDIF
*
      write(6,"('jmp =',2f9.5)") RMLJMP,TMLJMP*0.03
*
      RMLJMP = MAX(RMLJMP,TMLJMP*0.03)  !cold-water fix
C
      IF     (ITYPE.EQ.-1) THEN
        CALL MLD_CEN(R, P,   KDM, RMLJMP, DPMIXL)
      ELSEIF (ITYPE.EQ.0) THEN
        CALL MLD_PCM(R, P,   KDM, RMLJMP, DPMIXL)
      ELSEIF (ITYPE.EQ.1) THEN
        CALL MLD_PLM(R, THK, KDM, RMLJMP, DPMIXL)
      ELSEIF (ITYPE.EQ.2) THEN
        CALL MLD_LOC(R,T,S,P, THK, KDM, RMLJMP, DPMIXL)
      ELSEIF (ITYPE.EQ.3) THEN
        CALL MLD_PQM(R,T,S,P, THK, KDM, RMLJMP, DPMIXL)
      ELSE
        WRITE(6,*)
     +    'Usage: hycom_profile_mld archv.txt tmljmp [itype]'
        WRITE(6,*) 'unsupported itype value'
        CALL EXIT(1)
      ENDIF
C
C     OUTPUT.
C
      WRITE(6,'(a)')
     +  '#       day         MLD     SST     SSS        year  RMLJMP'
      WRITE(6,'(f11.2,f12.2,f8.2,f8.2,f12.5,f8.4)')
     +  DAYM,DPMIXL,T(1),S(1),YEARM,RMLJMP
      END

      SUBROUTINE RMLJMP0(RMLJMP, SST,SSS, IDM,JDM, TMLJMP)
      IMPLICIT NONE
C
      INTEGER       IDM,JDM
      REAL          RMLJMP(IDM,JDM),SST(IDM,JDM),SSS(IDM,JDM),TMLJMP
C
C --- CONVERT A TEMPERATURE JUMP INTO A DENSITY JUMP
C
      REAL       SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      INTEGER I,J
C
c-----------------------------------------------------------------------------
      real   dsigdt
      real   s,t
c
c --- coefficients for sigma-0 (based on Brydon & Sun fit)
      real       c1,c2,c3,c4,c5,c6,c7
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
c-----------------------------------------------------------------------------
C
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (SST(I,J).NE.SPVAL) THEN
            RMLJMP(I,J) = -TMLJMP*DSIGDT(SST(I,J),SSS(I,J))
          ELSE
            RMLJMP(I,J) = 0.0
          ENDIF
        ENDDO
      ENDDO
      END
      SUBROUTINE RMLJMP2(RMLJMP, SST,SSS, IDM,JDM, TMLJMP)
      IMPLICIT NONE
C
      INTEGER       IDM,JDM
      REAL          RMLJMP(IDM,JDM),SST(IDM,JDM),SSS(IDM,JDM),TMLJMP
C
C --- CONVERT A TEMPERATURE JUMP INTO A DENSITY JUMP
C
      REAL       SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      INTEGER I,J
C
c-----------------------------------------------------------------------------
      real   dsigdt
      real   s,t
c
c --- coefficients for sigma-2 (based on Brydon & Sun fit)
      real       c1,c2,c3,c4,c5,c6,c7
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
c-----------------------------------------------------------------------------
C
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (SST(I,J).NE.SPVAL) THEN
            RMLJMP(I,J) = -TMLJMP*DSIGDT(SST(I,J),SSS(I,J))
          ELSE
            RMLJMP(I,J) = 0.0
          ENDIF
        ENDDO
      ENDDO
      END
      SUBROUTINE RMLJMPL(RMLJMP, SST,SSS, IDM,JDM, TMLJMP)
      IMPLICIT NONE
C
      INTEGER       IDM,JDM
      REAL          RMLJMP(IDM,JDM),SST(IDM,JDM),SSS(IDM,JDM),TMLJMP
C
C --- CONVERT A TEMPERATURE JUMP INTO A DENSITY JUMP
C
      REAL       SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      INTEGER I,J
C
c-----------------------------------------------------------------------------
      real    sigloc,dsiglocdt,dsiglocds
      real    s,t,prs
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
      real    c1p,c2p,c3p,c4p,c5p,c6p,c7p
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
c-----------------------------------------------------------------------------
C
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (SST(I,J).NE.SPVAL) THEN
            RMLJMP(I,J) = -TMLJMP*dsiglocdt(SST(I,J),SSS(I,J),0.0)
          ELSE
            RMLJMP(I,J) = 0.0
          ENDIF
        ENDDO
      ENDDO
      END

      SUBROUTINE MLD_CEN(R, P, KDM, RMLJMP, DPMIXL)
      IMPLICIT NONE
C
      INTEGER KDM
      REAL    R(KDM),P(KDM+1), RMLJMP,DPMIXL
C
C     ORIGINAL MIXED LAYER DEPTH ALGORITHM, BASED ON CELL CENTERS.
C
      INTEGER K
      REAL    Q,RM,ZO,ZN
C
      RM = R(1) + RMLJMP
      ZN = 0.5*(P(1) + P(2))
      DO K= 2,KDM
        ZO = ZN
        ZN = 0.5*(P(K) + P(K+1))
        IF     (R(K).GE.RM) THEN
C
C          MLD IS BETWEEN ZO AND ZN
C
           Q      = (R(K)-RM)/(R(K)-R(K-1))   !R(K-1)<R(K)
           DPMIXL = Q*ZO + (1.0-Q)*ZN
           EXIT
         ELSEIF (P(K+1)+1.0.GT.P(KDM+1)) THEN  !.true. for k=kdm
C
C          MLD IS AT THE BOTTOM
C
           DPMIXL = P(KDM+1)
           EXIT
        ENDIF
      ENDDO
      END

      SUBROUTINE MLD_PCM(R, P, KDM, RMLJMP, DPMIXL)
      IMPLICIT NONE
C
      INTEGER KDM
      REAL    R(KDM),P(KDM+1), RMLJMP,DPMIXL
C
C     NEW MIXED LAYER DEPTH ALGORITHM, BASED ON PCM
C
      INTEGER K
      REAL    RM
C
      RM = R(1) + RMLJMP
      DO K= 2,KDM
        IF     (R(K).GE.RM) THEN
C
C          MLD IS AT THE INTERFACE ABOVE LAYER K
C
           DPMIXL = P(K)
           EXIT
         ELSEIF (P(K+1)+1.0.GT.P(KDM+1)) THEN  !.true. for k=kdm
C
C          MLD IS AT THE BOTTOM
C
           DPMIXL = P(KDM+1)
           EXIT
        ENDIF
      ENDDO
      END

      SUBROUTINE MLD_PLM(R, THK, KDM, SIGMLJ, DPMIXL)
      IMPLICIT NONE
C
      INTEGER KDM
      REAL    R(KDM),THK(KDM), SIGMLJ,DPMIXL
C
C     NEW MIXED LAYER DEPTH ALGORITHM, BASED ON SIMPLIFIED PLM.
C
      INTEGER K
      REAL    ZC(KDM+1)
      REAL    THSUR,THTOP,THBOT,THJMP(KDM)
C
      REAL       EPSIL
      PARAMETER (EPSIL=1.0e-11)
C
      ZC(1) = -0.5*THK(1)
      DO K= 2,KDM
        ZC(K) = ZC(K-1) - 0.5*(THK(K) + THK(K-1))
      ENDDO !K
      ZC(KDM+1) = ZC(KDM) - 0.5*THK(K)
C
      DPMIXL   = -ZC(KDM+1)  !bottom
      THJMP(1) = 0.0
      THSUR    = R(1)
      DO K= 2,KDM
        thsur    = min(R(k),thsur)  !ignore surface inversion
        thjmp(k) = max(R(k)-thsur,
     &                 thjmp(k-1)) !stable profile simplifies the code
*
        write(6,"('th,thsur,jmp =',2f9.3,f9.5)") 
     &    R(k),thsur,thjmp(k)
*
c
        if (thjmp(k).ge.sigmlj) then
c
c ---     find the two densities on the interface between layers
c ---     k-1 and k, using PLM assuming layers k-2 and k+1 are PCM.
c
          if     (k.eq.2) then
            thtop = thjmp(1)  
          else              
            thtop = thjmp(k-1) +
     &                min(thjmp(k-1)-thjmp(k-2),
     &                    thjmp(k)  -thjmp(k-1) )
          endif !k.eq.2:else
          if     (k.eq.kdm) then
            thbot = thjmp(k)
          else
            thsur      = min(R(k+1),thsur)
            thjmp(k+1) = max(R(k+1)-thsur,
     &                       thjmp(k))
            thbot = thjmp(k) -
     &                min(thjmp(k+1)-thjmp(k),
     &                    thjmp(k)  -thjmp(k-1) )
          endif !k.eq.klist:else
          if     (thtop.gt.thbot) then
            thtop = 0.5*(thtop+thbot)  !make stable at interface
            thbot = thtop
          endif
c                   
          if      (thtop.gt.sigmlj) then
c         
c ---       in bottom half of layer k-1
c
            dpmixl =
     &        -ZC(k-1) +
     &                 0.5*THK(k-1)*
     &                 (sigmlj+epsil-thjmp(k-1))/
     &                 (thtop +epsil-thjmp(k-1))
          elseif (thbot.ge.sigmlj) then
c         
c ---       at layer interface
c
            dpmixl =
     &        -ZC(k-1) + 0.5*THK(k-1)
          else
c         
c ---       in top half of layer k
c
            dpmixl =
     &        -ZC(k) -
     &                 0.5*THK(k)*
     &                 (1.0-(sigmlj  +epsil-thbot)/
     &                      (thjmp(k)+epsil-thbot) )
          endif !part of layer
*                 
            write (6,'(i3,a,3f7.3,f9.2)')
     &        k,
     &        '   thsur,top,bot,dpmixl =',
     &        thsur,thtop,thbot,dpmixl
*                 
          exit  !calculated dpmixl
        endif  !found dpmixl layer (thjmp(k).ge.sigmlj)
      enddo  !k
      END

      SUBROUTINE MLD_LOC(R,TT,SS,P, THK, KDM, SIGMLJ, DPMIXL)
      IMPLICIT NONE
C
      INTEGER KDM
      REAL    R(KDM),TT(KDM),SS(KDM),P(KDM),THK(KDM), SIGMLJ,DPMIXL
C
C     NEW MIXED LAYER DEPTH ALGORITHM, BASED ON SIMPLIFIED PLM.
C
      INTEGER K
      REAL    ZC(KDM+1)
      REAL    THSUR,THTOP,THBOT,THJMP(KDM)
      REAL    alfadt,betads,prsk,q,RKL
C
      REAL       EPSIL
      PARAMETER (EPSIL=1.0e-11)
C
c-----------------------------------------------------------------------------
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
c
      real    sigloc,dsiglocdt,dsiglocds
      real    s,t,prs
      real    c1p,c2p,c3p,c4p,c5p,c6p,c7p
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
c-----------------------------------------------------------------------------
C
      ZC(1) = -0.5*THK(1)
      DO K= 2,KDM
        ZC(K) = ZC(K-1) - 0.5*(THK(K) + THK(K-1))
      ENDDO !K
      ZC(KDM+1) = ZC(KDM) - 0.5*THK(K)
C
      DPMIXL   = -ZC(KDM+1)  !bottom
      THJMP(1) = 0.0
      THSUR    = R(1)
      RKL      = R(1)
      DO K= 2,KDM
        q     = (p(k+1)-p(k))/(p(k+1)-p(k-1)+0.1)
        prsk  = 9806.0*P(K)
        alfadt=((1.0-q)*dsiglocdt(TT(K-1),SS(K-1),prsk) +
     &               q *dsiglocdt(TT(K),  SS(K),  prsk)  )*
     &             (TT(K-1)-TT(K))
        betads=((1.0-q)*dsiglocds(TT(K-1),SS(K-1),prsk) +
     &               q *dsiglocds(TT(K),  SS(K),  prsk)  )*
     &             (SS(K-1)-SS(K))
        RKL=RKL-alfadt-betads
        thsur    = min(RKL,thsur)  !ignore surface inversion
        thjmp(k) = max(RKL-thsur,
     &                 thjmp(k-1)) !stable profile simplifies the code
*
        write(6,"('th,thsur,jmp,q,k =',2f9.3,f9.5,f6.3,i3)") 
     &    RKL,thsur,thjmp(k),q,k
*
c
        if (thjmp(k).ge.sigmlj) then
c
c ---     find the two densities on the interface between layers
c ---     k-1 and k, using PLM assuming layers k-2 and k+1 are PCM.
c
          if     (k.eq.2) then
            thtop = thjmp(1)  
          else              
            thtop = thjmp(k-1) +
     &                min(thjmp(k-1)-thjmp(k-2),
     &                    thjmp(k)  -thjmp(k-1) )
          endif !k.eq.2:else
          if     (k.eq.kdm) then
            thbot = thjmp(k)
          else
            prsk   = 9806.0*P(K+1)
            alfadt=0.5*(dsiglocdt(TT(K+1),SS(K+1),prsk  )+
     &                  dsiglocdt(TT(K),  SS(K),  prsk  ) )*
     &                 (TT(K)-TT(K+1))
            betads=0.5*(dsiglocds(TT(K+1),SS(K+1),prsk  )+
     &                  dsiglocds(TT(K),  SS(K),  prsk  ) )*
     &                 (SS(K)-SS(K+1))
            RKL=RKL-alfadt-betads
            thsur      = min(RKL,thsur)
            thjmp(k+1) = max(RKL-thsur,
     &                       thjmp(k))
            thbot = thjmp(k) -
     &                min(thjmp(k+1)-thjmp(k),
     &                    thjmp(k)  -thjmp(k-1) )
          endif !k.eq.klist:else
          if     (thtop.gt.thbot) then
            thtop = 0.5*(thtop+thbot)  !make stable at interface
            thbot = thtop
          endif
c                   
          if      (thtop.gt.sigmlj) then
c         
c ---       in bottom half of layer k-1
c
            dpmixl =
     &        -ZC(k-1) +
     &                 0.5*THK(k-1)*
     &                 (sigmlj+epsil-thjmp(k-1))/
     &                 (thtop +epsil-thjmp(k-1))
          elseif (thbot.ge.sigmlj) then
c         
c ---       at layer interface
c
            dpmixl =
     &        -ZC(k-1) + 0.5*THK(k-1)
          else
c         
c ---       in top half of layer k
c
            dpmixl =
     &        -ZC(k) -
     &                 0.5*THK(k)*
     &                 (1.0-(sigmlj  +epsil-thbot)/
     &                      (thjmp(k)+epsil-thbot) )
          endif !part of layer
*                 
          write (6,'(i3,a,2f9.5,2f9.2)')
     &      k,
     &      '   top,bot,intf,dpmixl =',
     &      thtop,thbot,-ZC(k-1)+0.5*THK(k-1),dpmixl
*                 

          exit  !calculated dpmixl
        endif  !found dpmixl layer (thjmp(k).ge.sigmlj)
      enddo  !k
      END

      SUBROUTINE MLD_PQM(R,TT,SS,P, THK, KDM, SIGMLJ, DPMIXL)
      IMPLICIT NONE
C
      INTEGER KDM
      REAL    R(KDM),TT(KDM),SS(KDM),P(KDM),THK(KDM), SIGMLJ,DPMIXL
C
C     NEW MIXED LAYER DEPTH ALGORITHM, BASED ON SIMPLIFIED PQM.
C
      INTEGER K
      REAL    ZC(KDM+1)
      REAL    THSUR,THTOP,THJMP(KDM)
      REAL    alfadt,betads,prsk,q,RKL,z
C
      REAL       EPSIL
      PARAMETER (EPSIL=1.0e-11)
C
c-----------------------------------------------------------------------------
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
c
      real    sigloc,dsiglocdt,dsiglocds
      real    s,t,prs
      real    c1p,c2p,c3p,c4p,c5p,c6p,c7p
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
c-----------------------------------------------------------------------------
C
      ZC(1) = -0.5*THK(1)
      DO K= 2,KDM
        ZC(K) = ZC(K-1) - 0.5*(THK(K) + THK(K-1))
      ENDDO !K
      ZC(KDM+1) = ZC(KDM) - 0.5*THK(K)
C
      DPMIXL   = -ZC(KDM+1)  !bottom
      THJMP(1) = 0.0
      THSUR    = R(1)
      RKL      = R(1)
      DO K= 2,KDM
        q     = (p(k+1)-p(k))/(p(k+1)-p(k-1)+0.1)
        prsk  = 9806.0*P(K)
        alfadt=((1.0-q)*dsiglocdt(TT(K-1),SS(K-1),prsk) +
     &               q *dsiglocdt(TT(K),  SS(K),  prsk)  )*
     &             (TT(K-1)-TT(K))
        betads=((1.0-q)*dsiglocds(TT(K-1),SS(K-1),prsk) +
     &               q *dsiglocds(TT(K),  SS(K),  prsk)  )*
     &             (SS(K-1)-SS(K))
        RKL=RKL-alfadt-betads
        thsur    = min(RKL,thsur)  !ignore surface inversion
        thjmp(k) = max(RKL-thsur,
     &                 thjmp(k-1)) !stable profile simplifies the code
*
        write(6,"('th,thsur,jmp,q,k =',2f9.3,f9.5,f6.3,i3)") 
     &    RKL,thsur,thjmp(k),q,k
*
c
        if (thjmp(k).ge.sigmlj) then
c
c ---     find the density on the interface between layers
c ---     k-1 and k, using the same cubic polynominal as PQM
c
          if     (k.eq.2) then
c ---       linear between cell centers
            thtop = thjmp(1) + (thjmp(2)-thjmp(1))*
     &                         thk(1)/(thk(1)+thk(2))
          elseif (k.eq.kdm) then
c ---       linear between cell centers
            thtop = thjmp(kdm) + (thjmp(kdm-1)-thjmp(kdm))*
     &                           thk(kdm)/(thk(kdm)+thk(kdm-1))
          else
            prsk   = 9806.0*P(K+1)
            alfadt=0.5*(dsiglocdt(TT(K+1),SS(K+1),prsk)+
     &                  dsiglocdt(TT(K),  SS(K),  prsk) )*
     &                 (TT(K)-TT(K+1))
            betads=0.5*(dsiglocds(TT(K+1),SS(K+1),prsk)+
     &                  dsiglocds(TT(K),  SS(K),  prsk) )*
     &                 (SS(K)-SS(K+1))
            RKL=RKL-alfadt-betads
            thsur      = min(RKL,thsur)
            thjmp(k+1) = max(RKL-thsur,
     &                       thjmp(k))
c
            z     = zc(k-1) - 0.5*thk(k-1)
            thtop = thjmp(k-2)*(                  (z      -zc(k-1))*
     &                          (z      -zc(k  ))*(z      -zc(k+1)) )/
     &                         (                  (zc(k-2)-zc(k-1))*
     &                          (zc(k-2)-zc(k  ))*(zc(k-2)-zc(k+1)) ) +
     &              thjmp(k-1)*((z      -zc(k-2))*                  
     &                          (z      -zc(k  ))*(z      -zc(k+1)) )/
     &                         ((zc(k-1)-zc(k-2))*                  
     &                          (zc(k-1)-zc(k  ))*(zc(k-1)-zc(k+1)) ) +
     &              thjmp(k  )*((z      -zc(k-2))*(z      -zc(k-1))*
     &                                            (z      -zc(k+1)) )/
     &                         ((zc(k  )-zc(k-2))*(zc(k  )-zc(k-1))*
     &                                            (zc(k  )-zc(k+1)) ) +
     &              thjmp(k+1)*((z      -zc(k-2))*(z      -zc(k-1))*
     &                          (z      -zc(k  ))                   )/
     &                         ((zc(k+1)-zc(k-2))*(zc(k+1)-zc(k-1))*
     &                          (zc(k+1)-zc(k  ))                   )
            thtop = max( thjmp(k-1), min( thjmp(k), thtop ) )
          endif !k.eq.2:k.eq.kdm:else
c                   
          if      (thtop.ge.sigmlj) then
c         
c ---       in bottom half of layer k-1, use linear interpolation
c
            dpmixl =
     &        -ZC(k-1) +
     &                 0.5*THK(k-1)*
     &                 (sigmlj+epsil-thjmp(k-1))/
     &                 (thtop +epsil-thjmp(k-1))
          else
c         
c ---       in top half of layer k, use linear interpolation
c
            dpmixl =
     &        -ZC(k) -
     &                 0.5*THK(k)*
     &                 (1.0-(sigmlj  +epsil-thtop)/
     &                      (thjmp(k)+epsil-thtop) )
          endif !part of layer
*                 
          write (6,'(i3,a,f9.5,2f9.2)')
     &      k,
     &      '   top,intf,dpmixl =',
     &      thtop,-ZC(k-1)+0.5*THK(k-1),dpmixl
*                 
          exit  !calculated dpmixl
        endif  !found dpmixl layer (thjmp(k).ge.sigmlj)
      enddo  !k
      END
