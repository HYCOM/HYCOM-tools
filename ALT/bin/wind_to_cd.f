      program wind_to_cd
      implicit none
c
c     usage:  echo wind vpmx airt pair sst | wind_to_cd
c
c     input:  wind vpmx airt pair sst
c     output: cd_coarep cd_coare cd_core2
c
c             wind = wind speed   (m/s)
c             vpmx = mixing ratio (kg/kg)
c             airt = air   temp   (degC)
c             pair = air pressure (Pa)
c             sst  = ocean temp   (degC)
c
c     alan j. wallcraft, naval research laboratory, april 2015
c
      real    cd_coarep,cd_coare,cd_core2  !functions
      real    wind,vpmx,airt,pair,sst
      integer ios
c
      write(6,'(2a)')
     &  '#cd_coarep  cd_coare  cd_core2',
     &  '      wind      vpmx      airt  pair-1e5       sst'
      do
        read(5,*,iostat=ios) wind,vpmx,airt,pair,sst
        if     (ios.ne.0) then
          exit
        endif
        write(6,'(8f10.4)') 
     &    cd_coarep(wind,vpmx,airt,pair,sst),  !wind in place of samo
     &    cd_coare(wind,vpmx,airt,sst),
     &    cd_core2(wind,vpmx,airt,sst),
     &    wind,vpmx,airt,pair-1.e5,sst
      enddo
      end

      real function cd_coare(wind,vpmx,airt,sst)
      implicit none
c
      real    wind,vpmx,airt,sst
c
c --- Wind stress drag coefficient * 10^3 from an approximation
c --- to the COARE 3.0 bulk algorithm (Fairall et al. 2003).
c
c --- wind = wind speed (m/s)
c --- vpmx = water vapor mixing ratio (kg/kg)
c --- airt = air temperature (C)
c --- sst  = sea temperature (C)
c
c ---              Ta-Ts
c ---           ==============
c ---   STABLE:  7    to  0.75 degC
c ---  NEUTRAL:  0.75 to -0.75 degC
c --- UNSTABLE: -0.75 to -8    degC
c
c ---              Va
c ---           ==============
c ---   Low:     1    to   5   m/s
c ---   High:    5    to  34   m/s
c
c --- vamax of 34 m/s from Sraj et al, 2013 (MWR-D-12-00228.1).
c
      real    tamts,q,qva,va
c
      real, parameter :: vamin=  1.0,  vamax=34.0
      real, parameter :: tdmin= -8.0,  tdmax= 7.0
      real, parameter :: tzero=273.16
c
      real, parameter ::
     &  as0_00=-0.06695,   as0_10= 0.09966,  as0_20=-0.02477,
     &  as0_01= 0.3133,    as0_11=-2.116,    as0_21= 0.2726,
     &  as0_02=-0.001473,  as0_12= 4.626,    as0_22=-0.5558,
     &  as0_03=-0.004056,  as0_13=-2.680,    as0_23= 0.3139

      real, parameter ::
     &  as5_00= 0.55815,   as5_10=-0.005593, as5_20= 0.0006024,
     &  as5_01= 0.08174,   as5_11= 0.2096,   as5_21=-0.02629,
     &  as5_02=-0.0004472, as5_12=-8.634,    as5_22= 0.2121,
     &  as5_03= 2.666e-6,  as5_13= 18.63,    as5_23= 0.7755

      real, parameter ::
     &  au0_00= 1.891,     au0_10=-0.006304, au0_20= 0.0004406,
     &  au0_01=-0.7182,    au0_11=-0.3028,   au0_21=-0.01769,
     &  au0_02= 0.1975,    au0_12= 0.3120,   au0_22= 0.01303,
     &  au0_03=-0.01790,   au0_13=-0.1210,   au0_23=-0.003394

      real, parameter ::
     &  au5_00= 0.6497,    au5_10= 0.003827, au5_20=-4.83e-5,
     &  au5_01= 0.06993,   au5_11=-0.2756,   au5_21= 0.007710,
     &  au5_02= 3.541e-5,  au5_12=-1.091,    au5_22=-0.2555,
     &  au5_03=-3.428e-6,  au5_13= 4.946,    au5_23= 0.7654

      real, parameter ::
     &  an0_00= 1.057,     an5_00= 0.6825,
     &  an0_01=-0.06949,   an5_01= 0.06945,
     &  an0_02= 0.01271,   an5_02=-0.0001029

      real, parameter ::
     &  ap0_10= as0_00 + as0_10*0.75 + as0_20*0.75**2,
     &  ap0_11=          as0_11*0.75 + as0_21*0.75**2,
     &  ap0_12=          as0_12*0.75 + as0_22*0.75**2,
     &  ap0_13=          as0_13*0.75 + as0_23*0.75**2

      real, parameter ::
     &  ap5_10= as5_00 + as5_10*0.75 + as5_20*0.75**2,
     &  ap5_11=          as5_11*0.75 + as5_21*0.75**2,
     &  ap5_12=          as5_12*0.75 + as5_22*0.75**2,
     &  ap5_13=          as5_13*0.75 + as5_23*0.75**2

      real, parameter ::
     &  am0_10= au0_00 - au0_10*0.75 + au0_20*0.75**2,
     &  am0_11=        - au0_11*0.75 + au0_21*0.75**2,
     &  am0_12=        - au0_12*0.75 + au0_22*0.75**2,
     &  am0_13=        - au0_13*0.75 + au0_23*0.75**2

      real, parameter ::
     &  am5_10= au5_00 - au5_10*0.75 + au5_20*0.75**2,
     &  am5_11=        - au5_11*0.75 + au5_21*0.75**2,
     &  am5_12=        - au5_12*0.75 + au5_22*0.75**2,
     &  am5_13=        - au5_13*0.75 + au5_23*0.75**2
c
c --- saturation specific humidity (lowe, j.appl.met., 16, 100-103, 1976)
      real qsatur,t
      qsatur(t)=.622e-3*(6.107799961e+00+t*(4.436518521e-01
     &               +t*(1.428945805e-02+t*(2.650648471e-04
     &               +t*(3.031240396e-06+t*(2.034080948e-08
     &               +t* 6.136820929e-11))))))
c
          tamts = airt-sst - 0.61*(airt+tzero)*(qsatur(airt)-vpmx)
          tamts = min( tdmax, max( tdmin, tamts ) )
          va    = max(vamin,min(vamax,wind))
          qva   = 1.0/va
          if     (va.le.5.0) then
            if     (tamts.ge. 0.75) then
              cd_coare = 
     &           (as0_00 + as0_01* va + as0_02* va**2 + as0_03* va**3)
     &         + (as0_10 + as0_11*qva + as0_12*qva**2 + as0_13*qva**3)
     &           *tamts
     &         + (as0_20 + as0_21*qva + as0_22*qva**2 + as0_23*qva**3)
     &           *tamts**2 
            elseif (tamts.le.-0.75) then
              cd_coare = 
     &           (au0_00 + au0_01* va + au0_02* va**2 + au0_03* va**3)
     &         + (au0_10 + au0_11*qva + au0_12*qva**2 + au0_13*qva**3)
     &           *tamts
     &         + (au0_20 + au0_21*qva + au0_22*qva**2 + au0_23*qva**3)
     &           *tamts**2 
            elseif (tamts.ge. -0.098)  then
              q =  (tamts+0.098)/0.848  !linear between  0.75 and -0.098
              cd_coare = q*
     &        (  (         as0_01* va + as0_02* va**2 + as0_03* va**3)
     &         + (ap0_10 + ap0_11*qva + ap0_12*qva**2 + ap0_13*qva**3)
     &        ) + (1.0-q)*
     &           (an0_00 + an0_01* va + an0_02* va**2)
            else
              q = (-tamts-0.098)/0.652  !linear between -0.75 and -0.098
              cd_coare = q*
     &        (  (         au0_01* va + au0_02* va**2 + au0_03* va**3)
     &         + (am0_10 + am0_11*qva + am0_12*qva**2 + am0_13*qva**3)
     &        ) + (1.0-q)*
     &           (an0_00 + an0_01* va + an0_02* va**2)
            endif !tamts
          else !va>5
            if     (tamts.ge. 0.75) then
              cd_coare = 
     &           (as5_00 + as5_01* va + as5_02* va**2 + as5_03* va**3)
     &         + (as5_10 + as5_11*qva + as5_12*qva**2 + as5_13*qva**3)
     &           *tamts
     &         + (as5_20 + as5_21*qva + as5_22*qva**2 + as5_23*qva**3)
     &           *tamts**2 
            elseif (tamts.le.-0.75) then
              cd_coare = 
     &           (au5_00 + au5_01* va + au5_02* va**2 + au5_03* va**3)
     &         + (au5_10 + au5_11*qva + au5_12*qva**2 + au5_13*qva**3)
     &           *tamts
     &         + (au5_20 + au5_21*qva + au5_22*qva**2 + au5_23*qva**3)
     &           *tamts**2 
            elseif (tamts.ge. -0.098)  then
              q =  (tamts+0.098)/0.848  !linear between  0.75 and -0.098
              cd_coare = q*
     &        (  (         as5_01* va + as5_02* va**2 + as5_03* va**3)
     &         + (ap5_10 + ap5_11*qva + ap5_12*qva**2 + ap5_13*qva**3)
     &        ) + (1.0-q)*
     &           (an5_00 + an5_01* va + an5_02* va**2)
            else
              q = (-tamts-0.098)/0.652  !linear between -0.75 and -0.098
              cd_coare = q*
     &        (  (         au5_01* va + au5_02* va**2 + au5_03* va**3)
     &         + (am5_10 + am5_11*qva + am5_12*qva**2 + am5_13*qva**3)
     &        ) + (1.0-q)*
     &           (an5_00 + an5_01* va + an5_02* va**2)
            endif !tamts
          endif !va
c
      end function cd_coare
c
      real function cd_core2(wind,vpmx,airt,sst)
      implicit none
c
      real    wind,vpmx,airt,sst
c
c --- Wind stress drag coefficient * 10^3 from
c --- the CORE v2 bulk algorithm (Large and Yeager, 2009).
c
c --- wind = wind speed (m/s)
c --- vpmx = water vapor mixing ratio (kg/kg)
c --- airt = air temperature (C)
c --- sst  = sea temperature (C)
c
c --- Added to HYCOM by Alexandra Bozec, FSU.
c
      integer it_a
      real    u10,v10,uw10,uw
      real    cd_n10,cd_n10_rt,ce_n10,ch_n10,cd_rt,stab,
     &        tv,tstar,qstar,bstar,zeta,x2,x,xx,
     &        psi_m,psi_h,z0,rair,qrair,zi
      real    cd10,ce10,ch10,ustar
c
      real, parameter :: vonkar=  0.4        !Von Karmann constant
      real, parameter ::  tzero=273.16       !celsius to kelvin offset
      real, parameter ::      g=  9.806      !same as HYCOM's g
c
c --- saturation specific humidity
      real qsatur5,t,qra
      qsatur5(t,qra)= 0.98*qra*6.40380e5*exp(-5107.4/(t+tzero))
c
c --- CORE v2 Large and Yeager 2009 Clim. Dyn.: The global climatology
c ---  of an interannually varying air-sea flux dataset.
c --- The bulk formulae effectively transform the problem of specifying
c --- the turbulent surface fluxes (at zi=10m) into one of describing 
c --- the near surface atmospheric state (wind, temperature and humidity).
c
*     write(6,'(a,1p4g14.5)')
*    &   '   wind,vpmx,airt,sst =',wind,vpmx,airt,sst
*
      rair = 1.22
      rair = 1.0/rair

      zi = 10.0
      tv = (airt+tzero)*(1.0+0.608*vpmx)  !in Kelvin
      uw = max(wind, 0.5)  !0.5 m/s floor on wind (undocumented NCAR)
      uw10 = uw            !first guess 10m wind

      cd_n10 = (2.7/uw10+0.142+0.0764*uw10)*1.0e-3         !L-Y eqn. 6a
*     write(6,'(a,1p4g14.5)')
*    &   'cd_n10 =',cd_n10,
*    &     (2.7/uw10)*1.0e-3,0.142*1.0e-3,(0.0764*uw10)*1.0e-3
      cd_n10_rt = sqrt(cd_n10)
      ce_n10 =  34.6 *cd_n10_rt*1.0e-3                     !L-Y eqn. 6b
      stab   = 0.5 + sign(0.5,airt-sst)
      ch_n10 = (18.0*stab+32.7*(1-stab))*cd_n10_rt*1.0e-3  !L-Y eqn. 6c

      cd10 = cd_n10  !first guess for exchange coeff's at z
      ch10 = ch_n10
      ce10 = ce_n10
*     write(6,'(a,1p4g14.5)')
*    &   'uw10,psi_m,cd_n10,cd10=',uw10,0.0,cd_n10,cd10
*
      do it_a= 1,2  !Monin-Obukhov iteration
        cd_rt = sqrt(cd10)
        ustar = cd_rt*uw                              !L-Y eqn. 7a
        tstar = (ch10/cd_rt)*(airt-sst)               !L-Y eqn. 7b
        qstar = (ce10/cd_rt)*
     &          (vpmx-qsatur5(sst,qrair))             !L-Y eqn. 7c
        bstar = g*(tstar/tv+qstar/(vpmx+1.0/0.608))
        zeta  = vonkar*bstar*zi/(ustar*ustar)         !L-Y eqn. 8a
        zeta  = sign( min(abs(zeta),10.0), zeta )     !undocumented NCAR
        x2 = sqrt(abs(1.0-16.0*zeta))                 !L-Y eqn. 8b
        x2 = max(x2, 1.0)                             !undocumented NCAR
        x  = sqrt(x2)

        if (zeta > 0.0) then
            psi_m = -5.0*zeta                         !L-Y eqn. 8c
            psi_h = -5.0*zeta                         !L-Y eqn. 8c
        else
            psi_m = log((1.0+2.0*x+x2)*(1.0+x2)/8.0)
     &            - 2.0*(atan(x)-atan(1.0))           !L-Y eqn. 8d
            psi_h = 2.0*log((1.0+x2)/2.0)             !L-Y eqn. 8e
        end if

        uw10 = uw/(1.0+cd_n10_rt*(log(zi/10.0)-psi_m) !L-Y eqn. 9
     &           /vonkar)
        cd_n10 = (2.7/uw10+0.142+0.0764*uw10)*1.0e-3  !L-Y eqn. 6a again
        cd_n10_rt = sqrt(cd_n10)
        ce_n10 = 34.6*cd_n10_rt*1.0e-3                !L-Y eqn. 6b again
        stab   = 0.5 + sign(0.5,zeta)
        ch_n10 =(18.0*stab+32.7*(1.0-stab))*cd_n10_rt*1.0e-3  !L-Y eqn. 6c again
        z0     = 10.0*exp(-vonkar/cd_n10_rt)          !diagnostic
c

        xx   = (log(zi/10.0)-psi_m)/vonkar
        cd10 = cd_n10/(1.0+cd_n10_rt*xx)**2           !L-Y 10a
        xx   = (log(zi/10.0)-psi_h)/vonkar
        ch10 = ch_n10/(1.0+ch_n10*xx/cd_n10_rt) *
     &                 sqrt(cd10/cd_n10)              !L-Y 10b
        ce10 = ce_n10/(1.0+ce_n10*xx/cd_n10_rt) *
     &                 sqrt(cd10/cd_n10)              !L-Y 10c
*       write(6,'(a,1p4g14.5)')
*    &     'uw10,psi_m,cd_n10,cd10=',uw10,psi_m,cd_n10,cd10
      enddo
c
      cd_core2 = cd10*1.0e3
      return
      end function cd_core2

      real function cd_coarep(samo,vpmx,airt,pair,sst)
      implicit none
c
      real    samo,vpmx,airt,pair,sst
c
c --- Wind stress drag coefficient * 10^3 from an approximation
c --- to the COARE 3.0 bulk algorithm (Fairall et al. 2003).
c
c --- samo = wind-ocean speed (m/s)
c --- vpmx = water vapor mixing ratio (kg/kg)
c --- airt = air temperature (C)
c --- pair = air pressure (Pa)
c --- sst  = sea temperature (C)
c
c ---              Ta-Ts
c ---           ==============
c ---   STABLE:  7    to  0.75 degC
c ---  NEUTRAL:  0.75 to -0.75 degC
c --- UNSTABLE: -0.75 to -8    degC
c
c ---              Va
c ---           ==============
c ---   Low:     1    to   5   m/s
c ---   High:    5    to  34   m/s
c
c --- vamax of 34 m/s from Sraj et al, 2013 (MWR-D-12-00228.1).
c
      real    tamts,q,qva,va
c
      real, parameter :: vamin=  1.0,  vamax=34.0
      real, parameter :: tdmin= -8.0,  tdmax= 7.0
      real, parameter :: tzero=273.16
c
      real, parameter ::
     &  as0_00=-0.06695,   as0_10= 0.09966,  as0_20=-0.02477,
     &  as0_01= 0.3133,    as0_11=-2.116,    as0_21= 0.2726,
     &  as0_02=-0.001473,  as0_12= 4.626,    as0_22=-0.5558,
     &  as0_03=-0.004056,  as0_13=-2.680,    as0_23= 0.3139

      real, parameter ::
     &  as5_00= 0.55815,   as5_10=-0.005593, as5_20= 0.0006024,
     &  as5_01= 0.08174,   as5_11= 0.2096,   as5_21=-0.02629,
     &  as5_02=-0.0004472, as5_12=-8.634,    as5_22= 0.2121,
     &  as5_03= 2.666e-6,  as5_13= 18.63,    as5_23= 0.7755

      real, parameter ::
     &  au0_00= 1.891,     au0_10=-0.006304, au0_20= 0.0004406,
     &  au0_01=-0.7182,    au0_11=-0.3028,   au0_21=-0.01769,
     &  au0_02= 0.1975,    au0_12= 0.3120,   au0_22= 0.01303,
     &  au0_03=-0.01790,   au0_13=-0.1210,   au0_23=-0.003394

      real, parameter ::
     &  au5_00= 0.6497,    au5_10= 0.003827, au5_20=-4.83e-5,
     &  au5_01= 0.06993,   au5_11=-0.2756,   au5_21= 0.007710,
     &  au5_02= 3.541e-5,  au5_12=-1.091,    au5_22=-0.2555,
     &  au5_03=-3.428e-6,  au5_13= 4.946,    au5_23= 0.7654

      real, parameter ::
     &  an0_00= 1.057,     an5_00= 0.6825,
     &  an0_01=-0.06949,   an5_01= 0.06945,
     &  an0_02= 0.01271,   an5_02=-0.0001029

      real, parameter ::
     &  ap0_10= as0_00 + as0_10*0.75 + as0_20*0.75**2,
     &  ap0_11=          as0_11*0.75 + as0_21*0.75**2,
     &  ap0_12=          as0_12*0.75 + as0_22*0.75**2,
     &  ap0_13=          as0_13*0.75 + as0_23*0.75**2

      real, parameter ::
     &  ap5_10= as5_00 + as5_10*0.75 + as5_20*0.75**2,
     &  ap5_11=          as5_11*0.75 + as5_21*0.75**2,
     &  ap5_12=          as5_12*0.75 + as5_22*0.75**2,
     &  ap5_13=          as5_13*0.75 + as5_23*0.75**2

      real, parameter ::
     &  am0_10= au0_00 - au0_10*0.75 + au0_20*0.75**2,
     &  am0_11=        - au0_11*0.75 + au0_21*0.75**2,
     &  am0_12=        - au0_12*0.75 + au0_22*0.75**2,
     &  am0_13=        - au0_13*0.75 + au0_23*0.75**2

      real, parameter ::
     &  am5_10= au5_00 - au5_10*0.75 + au5_20*0.75**2,
     &  am5_11=        - au5_11*0.75 + au5_21*0.75**2,
     &  am5_12=        - au5_12*0.75 + au5_22*0.75**2,
     &  am5_13=        - au5_13*0.75 + au5_23*0.75**2
c
      real satvpr,qsaturp,t,t6,p6
c
c --- saturation vapor pressure (Pa),
c --- from a polynominal approximation (lowe, j.appl.met., 16, 100-103, 1976)
      satvpr(t)=  100.0*(6.107799961e+00+t*(4.436518521e-01
     &               +t*(1.428945805e-02+t*(2.650648471e-04
     &               +t*(3.031240396e-06+t*(2.034080948e-08
     &               +t* 6.136820929e-11))))))
c
c --- pressure dependent saturation mixing ratio (kg/kg)
c --- p6 is pressure in Pa
      qsaturp(t6,p6)=0.622*(satvpr(t6)/(p6-satvpr(t6)))
c
c ---     correct tamts to 100% humidity
          tamts = airt-sst -
     &            0.608*(airt+tzero)*(qsaturp(airt,pair)-vpmx)
          tamts = min( tdmax, max( tdmin, tamts ) )
          va    = max( vamin, min( vamax, samo  ) )
          qva   = 1.0/va
          if     (va.le.5.0) then
            if     (tamts.ge. 0.75) then
              cd_coarep =
     &           (as0_00 + as0_01* va + as0_02* va**2 + as0_03* va**3)
     &         + (as0_10 + as0_11*qva + as0_12*qva**2 + as0_13*qva**3)
     &           *tamts
     &         + (as0_20 + as0_21*qva + as0_22*qva**2 + as0_23*qva**3)
     &           *tamts**2 
            elseif (tamts.le.-0.75) then
              cd_coarep =
     &           (au0_00 + au0_01* va + au0_02* va**2 + au0_03* va**3)
     &         + (au0_10 + au0_11*qva + au0_12*qva**2 + au0_13*qva**3)
     &           *tamts
     &         + (au0_20 + au0_21*qva + au0_22*qva**2 + au0_23*qva**3)
     &           *tamts**2 
            elseif (tamts.ge. -0.098)  then
              q =  (tamts+0.098)/0.848  !linear between  0.75 and -0.098
              cd_coarep = q*
     &        (  (         as0_01* va + as0_02* va**2 + as0_03* va**3)
     &         + (ap0_10 + ap0_11*qva + ap0_12*qva**2 + ap0_13*qva**3)
     &        ) + (1.0-q)*
     &           (an0_00 + an0_01* va + an0_02* va**2)
            else
              q = (-tamts-0.098)/0.652  !linear between -0.75 and -0.098
              cd_coarep = q*
     &        (  (         au0_01* va + au0_02* va**2 + au0_03* va**3)
     &         + (am0_10 + am0_11*qva + am0_12*qva**2 + am0_13*qva**3)
     &        ) + (1.0-q)*
     &           (an0_00 + an0_01* va + an0_02* va**2)
            endif !tamts
          else !va>5
            if     (tamts.ge. 0.75) then
              cd_coarep = 
     &           (as5_00 + as5_01* va + as5_02* va**2 + as5_03* va**3)
     &         + (as5_10 + as5_11*qva + as5_12*qva**2 + as5_13*qva**3)
     &           *tamts
     &         + (as5_20 + as5_21*qva + as5_22*qva**2 + as5_23*qva**3)
     &           *tamts**2 
            elseif (tamts.le.-0.75) then
              cd_coarep = 
     &           (au5_00 + au5_01* va + au5_02* va**2 + au5_03* va**3)
     &         + (au5_10 + au5_11*qva + au5_12*qva**2 + au5_13*qva**3)
     &           *tamts
     &         + (au5_20 + au5_21*qva + au5_22*qva**2 + au5_23*qva**3)
     &           *tamts**2 
            elseif (tamts.ge. -0.098)  then
              q =  (tamts+0.098)/0.848  !linear between  0.75 and -0.098
              cd_coarep = q*
     &        (  (         as5_01* va + as5_02* va**2 + as5_03* va**3)
     &         + (ap5_10 + ap5_11*qva + ap5_12*qva**2 + ap5_13*qva**3)
     &        ) + (1.0-q)*
     &           (an5_00 + an5_01* va + an5_02* va**2)
            else
              q = (-tamts-0.098)/0.652  !linear between -0.75 and -0.098
              cd_coarep = q*
     &        (  (         au5_01* va + au5_02* va**2 + au5_03* va**3)
     &         + (am5_10 + am5_11*qva + am5_12*qva**2 + am5_13*qva**3)
     &        ) + (1.0-q)*
     &           (an5_00 + an5_01* va + an5_02* va**2)
            endif !tamts
          endif !va
c
      end function cd_coarep
