      PROGRAM HYCOM_PROFILE_REMAP
      IMPLICIT NONE
C
C  hycom_profile_remap - Usage: hycom_profile_remap archv.txt map.txt archr.txt [itype]
C
C                 converts an HYCOM isopycnal text profile file to
C                 a text profile file across a new set of layers
C
C   archv.txt is assumed to be an HYCOM archive text profile file
C   map.txt   is assumed to contain a HYCOM hybrid grid definition.
C   archr.txt will be the output, remapped, text profile file
C   itype     is the input interpolation type (default 0)
C                =0; piecewise constant  method (PCM) or donor cell
C                =1; piecewise linear    method (PLM) or VanLeer
C
C  all input interpolation schemes conserve the original cell average,
C  and the output is the average of the interpolation profile across
C  each new cell so the total depth average is conserved.  This is
C  "remapping" because we choose the location of the output layers.
C
C  If the output layer locations are known, use hycom_profile2zi.
C
C  this version for "serial" Unix systems.
C
C  Tim Campbell,       Mississippi State University
C  Alan J. Wallcraft,  Naval Research Laboratory,  March 2003.
C
      INTEGER      IARGC
      INTEGER      NARG
      CHARACTER*240 CARG
C
      CHARACTER*240 CFILEA,CFILEB,CFILEC
      CHARACTER*240 CLINE
      CHARACTER*26  SIGFMT
      REAL          THK,FLAG
      REAL          U(9999),V(9999),T(9999),S(9999),R(9999),P(9999+1)
      REAL          U0(9999),V0(9999),T0(9999),S0(9999),R0(9999)
      REAL          UZ(9999),VZ(9999),TZ(9999),SZ(9999),RZ(9999),
     +              DP0(9999),ZZ(9999),PZ(9999+1)
      INTEGER       IOS,ITYPE,K,KDM,KZ
      INTEGER       I
      REAL          UT,VT,TT,ST,RT
      REAL          UZT,VZT,TZT,SZT,RZT
      REAL          DUT,DVT,DTT,DST,DRT
      REAL          DUM,DVM,DTM,DSM,DRM
C
      INTEGER       KDMNEW,NHYBRD,NSIGMA,THFLAG,HYBFLG
      REAL          DP00,DP00X,DP00F,DS00,DS00X,DS00F,SIGMA(9999)
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.4) THEN
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CFILEB)
        CALL GETARG(3,CFILEC)
        CALL GETARG(4,CLINE)
        READ(CLINE,*) ITYPE
      ELSEIF (NARG.EQ.3) THEN
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CFILEB)
        CALL GETARG(3,CFILEC)
        ITYPE = 0
      ELSE
        WRITE(6,*)
     +  'Usage: hycom_profile_remap archv.txt map.txt archr.txt [itype]'
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
      OPEN(UNIT=12, FILE=CFILEB, FORM='FORMATTED', STATUS='OLD',
     +     IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',CFILEB(1:LEN_TRIM(CFILEB))
        WRITE(6,*) 'ios   = ',ios
        CALL EXIT(4)
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
C
C     READ THE ISOPYCNAL PROFILE.
C
      P(1) =  0.0
      KDM  = -1
      DO K= 1,9999
        READ(11,'(a)',IOSTAT=IOS) CLINE
        IF     (IOS.NE.0) THEN
          IF     (K.NE.KDM+1) THEN
            WRITE(6,*) 'Error: inconsistent input profile'
            CALL EXIT(6)
          ENDIF
          EXIT
        ENDIF
        READ(CLINE,*) KDM,U(K),V(K),T(K),S(K),R(K),THK
        P(K+1) = P(K) + THK
        IF     (THK.EQ.0.0) THEN
          U(K) = U(K-1)
          V(K) = V(K-1)
          T(K) = T(K-1)
          S(K) = S(K-1)
          R(K) = R(K-1)
        ENDIF
      ENDDO
      CLOSE(11)
C
C     STORE COPY OF INPUT AVERAGES AND
C     COMPUTE INPUT TOTAL COLUMN AVERAGES
C
      UT=0.0
      VT=0.0
      TT=0.0
      ST=0.0
      RT=0.0
      DO K=1,KDM
        THK=P(K+1)-P(K)
        UT=UT+THK*U(K)
        VT=VT+THK*V(K)
        TT=TT+THK*T(K)
        ST=ST+THK*S(K)
        RT=RT+THK*R(K)
        U0(K)=U(K)
        V0(K)=V(K)
        T0(K)=T(K)
        S0(K)=S(K)
        R0(K)=R(K)
      ENDDO
      UT=UT/P(KDM+1)
      VT=VT/P(KDM+1)
      TT=TT/P(KDM+1)
      ST=ST/P(KDM+1)
      RT=RT/P(KDM+1)
C
C     READ THE MAP FILE, ON UNIT 12.
C
C --- 'kdm   ' = number of output layers
C --- 'nhybrd' = number of output hybrid levels (0=all isopycnal)
C --- 'nsigma' = number of output sigma  levels (nhybrd-nsigma z-levels)
C --- 'dp00'   = deep    z-level spacing minimum thickness (m)
C --- 'dp00x'  = deep    z-level spacing maximum thickness (m)
C --- 'dp00f'  = deep    z-level spacing stretching factor (1.0=const.z)
C --- 'ds00'   = shallow z-level spacing minimum thickness (m)
C --- 'ds00x'  = shallow z-level spacing maximum thickness (m)
C --- 'ds00f'  = shallow z-level spacing stretching factor (1.0=const.z)
C
C --- the above specifies a vertical coord. that is isopycnal or:
C ---     z in    deep water, based on dp00,dp00x,dp00f
C ---     z in shallow water, based on ds00,ds00x,ds00f and nsigma
C ---     sigma between them, based on ds00,ds00x,ds00f and nsigma
C --- for z-only set nsigma=0 (and ds00,ds00x,ds00f=dp00,dp00x,dp00f)
C --- for sigma-z (shallow-deep) use a very small ds00
C ---  (pure sigma-z also has ds00f=dp00f and ds00x=dp00x*ds00/dp00)
C --- for z-sigma (shallow-deep) use a very large dp00 (not recommended)
C --- for sigma-only set nsigma=kdm, dp00 large, and ds00 small
C
      call blkini(kdmnew,'kdm   ')
      call blkini(nhybrd,'nhybrd')
      call blkini(nsigma,'nsigma')
      call blkinr(dp00,  'dp00  ','(a6," =",f10.4," m")')
      call blkinr(dp00x, 'dp00x ','(a6," =",f10.4," m")')
      call blkinr(dp00f, 'dp00f ','(a6," =",f10.4," ")')
      call blkinr(ds00,  'ds00  ','(a6," =",f10.4," m")')
      call blkinr(ds00x, 'ds00x ','(a6," =",f10.4," m")')
      call blkinr(ds00f, 'ds00f ','(a6," =",f10.4," ")')
C
C --- 'thflag' = reference pressure flag (0=Sigma-0, 2=Sigma-2, 4=Sigma-4)
      call blkini(thflag,'thflag')
c
      if     (thflag.eq.0) then
        sigfmt = '(a6," =",f10.4," sigma-0")'
      elseif (thflag.eq.2) then
        sigfmt = '(a6," =",f10.4," sigma-2")'
      elseif (thflag.eq.4) then
        sigfmt = '(a6," =",f10.4," sigma-4")'
      endif
C
C --- 'sigma ' = layer densities (sigma units)
      do k=1,kdmnew
        call blkinr(sigma(k),'sigma ',sigfmt)
C
        if     (k.gt.1) then
          if      (sigma(k).le.sigma(k-1)) then
            WRITE(6,*)
     &  'Usage: hycom_profile_remap archv.txt map.txt archr.txt [itype]'
            WRITE(6,*)
     &  'sigma (in map.txt) is not stabaly stratified'
            CALL EXIT(2)
          endif
        endif
      enddo
C
C --- 'hybflg' = hybrid generator  flag (0=T&S, 1=th&S, 2=th&T)
      call blkini(hybflg,'hybflg')
C
      CLOSE(12)
C
C     FIND THE MINIMUM OUTPUT LAYER THICKNESSES.
C
      KZ = KDMNEW
      CALL MIN_THICK(DP0,KZ,
     &               P(KDM+1), !depth
     &               NHYBRD,NSIGMA,DP00,DP00X,DP00F,DS00,DS00X,DS00F)
C
C     FIND THE NEW INTERFACE DEPTHS.
C
      IF     (ITYPE.EQ.0 .OR. ITYPE.EQ.1) THEN
        CALL LAYERS_PXM(PZ, DP0,SIGMA,KZ, R,P,KDM, ITYPE)
      ELSE
        WRITE(6,*)
     +  'Usage: hycom_profile_remap archv.txt map.txt archr.txt [itype]'
        WRITE(6,*) 'unsupported itype value'
        CALL EXIT(1)
      ENDIF
C
C     FORM NEW CELL AVERAGES
C
      FLAG = 0.0
      IF     (ITYPE.EQ.0) THEN
        CALL LAYER2ZI_PCM(U, V, T, S, R, P, KDM, 
     &                    UZ,VZ,TZ,SZ,RZ,PZ,KZ, FLAG)
      ELSEIF (ITYPE.EQ.1) THEN
        CALL LAYER2ZI_PLM(U, V, T, S, R, P, KDM, 
     &                    UZ,VZ,TZ,SZ,RZ,PZ,KZ, FLAG)
      ELSE
        WRITE(6,*)
     +  'Usage: hycom_profile_remap archv.txt map.txt archr.txt [itype]'
        WRITE(6,*) 'unsupported itype value'
        CALL EXIT(1)
      ENDIF
C
C     COMPUTE AND PRINT NEW TOTAL COLUMN AVERAGE AND DEVIATION METRIC
C
      UZT=0.0
      VZT=0.0
      TZT=0.0
      SZT=0.0
      RZT=0.0
      DUT=0.0
      DVT=0.0
      DTT=0.0
      DST=0.0
      DRT=0.0
      DUM=0.0
      DVM=0.0
      DTM=0.0
      DSM=0.0
      DRM=0.0
      DO K=1,KZ
        THK=PZ(K+1)-PZ(K)
        UZT=UZT+THK*UZ(K)
        VZT=VZT+THK*VZ(K)
        TZT=TZT+THK*TZ(K)
        SZT=SZT+THK*SZ(K)
        RZT=RZT+THK*RZ(K)
      ENDDO
      UZT=UZT/PZ(KZ+1)
      VZT=VZT/PZ(KZ+1)
      TZT=TZT/PZ(KZ+1)
      SZT=SZT/PZ(KZ+1)
      RZT=RZT/PZ(KZ+1)
      WRITE(*,'(/A)')       '        LABEL:  U  V  T  S  R'
      WRITE(*,'(A,5e14.6)') ' INPUT TOTALS: ',UT ,VT ,TT ,ST ,RT
      WRITE(*,'(A,5e14.6)') 'OUTPUT TOTALS: ',UZT,VZT,TZT,SZT,RZT
C
C     OUTPUT.
C
      WRITE(21,'(a,a)')
     &  '#   k',
     &  '    utot    vtot    temp    saln    dens    thkns      dpth'
      DO K= 1,KZ
        WRITE(21,'(i4,1x,2f8.2,3f8.3,f9.3,f10.3)')
     &    K,
     &    UZ(K),               !cm/s
     &    VZ(K),               !cm/s
     &    TZ(K),               !degC
     &    SZ(K),               !psu
     &    RZ(K),               !SigmaT
     &    PZ(K+1)-PZ(K),       !m
     &    ZZ(K)                !m
      ENDDO
      END

      SUBROUTINE MIN_THICK(DP0,KZ,
     &                     DEPTH,
     &                     NHYBRD,NSIGMA,
     &                     DP00,DP00X,DP00F,DS00,DS00X,DS00F)
      IMPLICIT NONE
c
      INTEGER KZ,NHYBRD,NSIGMA
      REAL    DP0(KZ), DEPTH, DP00,DP00X,DP00F,DS00,DS00X,DS00F
c
c**********
c*
c  1) return minimum layer thicknesses for hybrid coordinates.
c
c  2) input arguments:
c      kz     - number of layers
c      depth  - total depth
c      nhybrd - number of hybrid levels (0=all isopycnal)
c      nsigma - number of sigma  levels (nhybrd-nsigma z-levels)
c      dp00   - deep    z-level spacing minimum thickness (m)
c      dp00x  - deep    z-level spacing maximum thickness (m)
c      dp00f  - deep    z-level spacing stretching factor (1.0=const.z)
c      ds00   - shallow z-level spacing minimum thickness (m)
c      ds00x  - shallow z-level spacing maximum thickness (m)
c      ds00f  - shallow z-level spacing stretching factor (1.0=const.z)
c
c  3) output arguments:
c       dp0   - required minimum layer thicknesses
c
c  5) the minimum thicknesses are from the following fixed coordinates:
c       z in    deep water, based on dp00,dp00x,dp00f
c       z in shallow water, based on ds00,ds00x,ds00f and nsigma
c       sigma between them, based on ds00,ds00x,ds00f and nsigma
c
c  6) Alan J. Wallcraft, Naval Research Laboratory, March 2003.
c*
c**********
c
      INTEGER K
      REAL    DP0K(KZ),DP0KF,DS0K(KZ),DS0KF,DSMS,DSSK(KZ)
C
C     DEEP MINIMUM THICKNESSES.
C
      IF     (NHYBRD.GT.0) THEN
        DP0K(1)=DP00
      ELSE
        DP0K(1)=0.0
      ENDIF
C
      DP0KF=1.0
      DO K= 2,KZ
        DP0KF=DP0KF*DP00F
        IF     (K.LE.NHYBRD) THEN
          DP0K(K)=MIN(DP00*DP0KF,DP00X)
        ELSE
          DP0K(K)=0.0
        ENDIF
      ENDDO
C
C     SHALLOW MINIMUM THICKNESSES.
C
      IF     (NHYBRD.GT.0) THEN
        DS0K(1)=DS00
      ELSE
        DS0K(1)=0.0
      ENDIF
      DSMS = DS0K(1)
C
      DS0KF=1.0
      DO K= 2,KZ
        DS0KF=DS0KF*DS00F
        IF     (K.LE.NHYBRD) THEN
          DS0K(K)=MIN(DS00*DS0KF,DS00X)
        ELSE
          DS0K(K)=0.0
        ENDIF
        DSMS = DSMS + DS0K(K)
      ENDDO
C
C     SIGMA-DEPTH SCALE FACTORS
C
      DO K= 1,NSIGMA
        DSSK(K)=DS0K(K)/DSMS  ! fraction of depths in sigma layer k
      ENDDO
      DO K= NSIGMA+1,KZ
        DS0K(K)=DP0K(K)
        DSSK(K)=0.0           ! these layers are zero in sigma mode
      ENDDO
C
C     ESTABLISH ACTUAL MINIMUM LAYER THICKNESSES
C
      DO K= 1,KZ
        DP0(K) = MIN( DP0K(K),
     &                MAX( DS0K(K),
     &                     DSSK(K)*DEPTH ) )
      ENDDO
      RETURN
      END

      SUBROUTINE LAYERS_PXM(PZ, DP0,SIGMA,KZ, R,P,KDM, ITYPE)
      IMPLICIT NONE
c
      INTEGER KDM,KZ,ITYPE
      REAL    PZ(KZ+1), DP0(KZ),SIGMA(KZ), R(KDM),P(KDM+1)
c
c**********
c*
c  1) Choose the location in z-space of a hybrid (isopycnal/z)
c     coordinate system applied against an input density profile.
c
c  2) input arguments:
c      dp0    - output layer minimum thicknesses
c      sigma  - output layer target densities
c      kz     - number of output layers
c      r      - input density in p-layer space
c      p      - input layer interface depths (non-negative m)
c                 p(    1) is the surface
c                 p(  k+1) >= p(k)
c                 p(kdm+1) is the bathymetry
c      kdm    - number of  input layers
c      itype  - input interpolation type
c                =0; piecewise constant  method (PCM) or donor cell
c                =1; piecewise linear    method (PLM) or VanLeer
c
c  3) output arguments:
c       pz    - required interface depths (non-negative m)
c                 pz(   1)  = p(1)
c                 pz(k +1) >= pz(k)
c                 pz(kz+1)  = p(kdm+1)
c
c  5) Alan J. Wallcraft, Naval Research Laboratory, March 2003.
c*
c**********
c
      INTEGER K,KK
      REAL    DEPTH, DP,DPN,PSMX,PSUM,RSMX,RSUM,RS(KDM)
C
C     MAKE THE DENSITY PROFILE STABLE.
C     NOTE THAT THIS IS ONLY USED IN CALCULATING PZ.
C
      RS(1) = R(1)
      DO K= 2,KDM
        RS(K) = MAX( R(K), RS(K-1) )
      ENDDO
C
C     FIND EACH LAYER'S LOCATION.
C     GREEDY FROM TOP FOR ISOPYCNALS (Rainer Bleck's REGRID)
C
      DEPTH = P(KDM+1)
C
      PZ(   1) = 0.0
      PZ(KZ+1) = DEPTH
C
      DO KK= 1,KZ-1
C
C       FIND THE THICKNESS OF AN ISOPYCNAL LAYER STARTING AT PZ(KK)
C
        RSUM = 0.0
        PSUM = 0.0
        DO K= 1,KDM
          IF     (P(K+1).GT.PZ(KK)) THEN  !input layer is deep enough
            IF     (PSUM.EQ.0.0) THEN  !1st layer in sum
              IF     (RS(K).GT.SIGMA(KK)) THEN  !no isopycnal layer
                PSUM = 0.0
                EXIT
              ELSEIF (RS(K).EQ.SIGMA(KK)) THEN  !exactly the needed layer
                PSUM = P(K+1) - PZ(KK)
                RSUM = RS(K)*PSUM
                EXIT
              ELSE  !inside the needed layer, start the sum
                PSUM = P(K+1) - PZ(KK)
                IF     (ITYPE.EQ.0) THEN  !pcm
                  RSUM = RS(K)*PSUM
                ELSE  !plm
                  write(6,*) 'ERROR in LAYERS_PXM',
     &                       ' - PLM not implemented'
                  STOP
                ENDIF
              ENDIF
            ELSE  !2nd or deeper layer in sum
              DP   = P(K+1) - P(K)
              IF     (DP.GT.1.E-6) THEN
                PSMX = PSUM + DP
                RSMX = RSUM + RS(K)*DP
                IF     (RSMX/PSMX.GE.SIGMA(KK)) THEN  !last layer in sum
                  !find layer thickness to complete an isopycnal layer
                  IF     (ITYPE.EQ.0) THEN  !pcm
                    DPN  = (SIGMA(KK)*PSUM - RSUM) /
     &                         (RS(K) - SIGMA(KK))
                  ELSE  !plm
                    write(6,*) 'ERROR in LAYERS_PXM',
     &                         ' - PLM not implemented'
                    STOP
                  ENDIF
                  PSUM = PSUM + MIN( DPN, DP )  !should have dpn<=dp
                  RSUM = SIGMA(KK)*PSUM  ! found an isopycnal layer
                  EXIT
                ELSE  !inside the needed layer, add to sum
                  PSUM = PSMX
                  RSUM = RSMX
                ENDIF
              ENDIF  !non-zero layer
            ENDIF
          ENDIF  !deep enough layer
        ENDDO  !k
C
        PZ(KK+1) = PZ(KK) + MAX( PSUM, DP0(KK) )
C
        IF     (PZ(KK+1).GE.DEPTH) THEN  !at the bottom
          DO K= KK+1,KZ
            PZ(K) = DEPTH
          ENDDO
          EXIT
        ENDIF
      ENDDO  !kk
C
      RETURN
      END

      subroutine layer2zi_pcm(u, v, t, s, r, p, kk,
     &                        uz,vz,tz,sz,rz,pz,kz, flag)
      implicit none
c
      integer kk,kz
      real    u( kk),v( kk),t( kk),s( kk),r( kk),p( kk+1),
     &        uz(kz),vz(kz),tz(kz),sz(kz),rz(kz),pz(kz+1),flag
c
c**********
c*
c  1) remap from one set of vertical cells to another.
c     method: piecewise constant across each input cell
c             the output is the average of the interpolation
c             profile across each output cell.
c
c  2) input arguments:
c       u,v   - scalar fields in p-layer space
c       t,s,r - scalar fields in p-layer space
c       p     - layer interface depths (non-negative m)
c                 p(   1) is the surface
c                 p(kk+1) is the bathymetry
c       kk    - dimension of a  (number of  input layers)
c       pz    - target interface depths (non-negative m)
c                 pz(k+1) >= pz(k)
c       flag  - data void (land) marker
c       kz    - dimension of az (number of output layers)
c
c  3) output arguments:
c       uz    - scalar field in pz-layer space
c       vz    - scalar field in pz-layer space
c       tz    - scalar field in pz-layer space
c       sz    - scalar field in pz-layer space
c       rz    - scalar field in pz-layer space
c
c  4) except at data voids, must have:
c           p(   1) == zero (surface)
c           p( l+1) >= p(:,:,l)
c           p(kk+1) == bathymetry
c           0 <= pz(k) <= pz(k+1)
c      output layers completely below the bathymetry inherit values
c      from the layer above.
c
c  5) Alan J. Wallcraft, Naval Research Laboratory, September 2002.
c*
c**********
c
      real       thin 
      parameter (thin=1.e-6)  ! minimum layer thickness (no division by 0.0)
c
      integer k,l,lf
      real    q,zb,zt,uzk,vzk,tzk,szk,rzk
c
      if     (r(1).eq.flag) then
        do k= 1,kz
          uz(k) = flag  ! land
          vz(k) = flag  ! land
          tz(k) = flag  ! land
          sz(k) = flag  ! land
          rz(k) = flag  ! land
        enddo
      else
        lf=1
        zb=pz(1)
        do k= 1,kz
          zt = zb
          zb = pz(k+1)
*         WRITE(6,*) 'k,zt,zb = ',k,zt,zb
          if     (zb-zt.lt.thin .or. zt.ge.p(kk+1)) then
c
c ---       thin or bottomed layer, values taken from layer above
c
            uz(k) = uz(k-1)
            vz(k) = vz(k-1)
            tz(k) = tz(k-1)
            sz(k) = sz(k-1)
            rz(k) = rz(k-1)
          else
c
c           form layer averages.
c
            if     (p(lf).gt.zt) then
              WRITE(6,*) 'bad lf = ',lf
              stop
            endif
            uzk = 0.0
            vzk = 0.0
            tzk = 0.0
            szk = 0.0
            rzk = 0.0
            do l= lf,kk
              if     (p(l).gt.zb) then
*               WRITE(6,*) 'l,lf= ',l,lf,l-1
                lf = l-1
                exit
              elseif (p(l).ge.zt .and. p(l+1).le.zb) then
c
c               the input layer is completely inside the output layer
c
                q   = max(p(l+1)-p(l),thin)/(zb-zt)
                uzk = uzk + q*u(l)
                vzk = vzk + q*v(l)
                tzk = tzk + q*t(l)
                szk = szk + q*s(l)
                rzk = rzk + q*r(l)
*               WRITE(6,*) 'L,q = ',l,q
              else
c
c               the input layer is partially inside the output layer
c
                q   = max(min(p(l+1),zb)-max(p(l),zt),thin)/(zb-zt)
                uzk = uzk + q*u(l)
                vzk = vzk + q*v(l)
                tzk = tzk + q*t(l)
                szk = szk + q*s(l)
                rzk = rzk + q*r(l)
*               WRITE(6,*) 'l,q = ',l,q
              endif
            enddo !l
            uz(k) = uzk
            vz(k) = vzk
            tz(k) = tzk
            sz(k) = szk
            rz(k) = rzk
          endif
        enddo !k
      endif
      return
      end subroutine layer2zi_pcm

      subroutine layer2zi_plm(u, v, t, s, r, p, kk,
     &                        uz,vz,tz,sz,rz,pz,kz, flag)
      implicit none
c
      integer kk,kz
      real    u( kk),v( kk),t( kk),s( kk),r( kk),p( kk+1),
     &        uz(kz),vz(kz),tz(kz),sz(kz),rz(kz),pz(kz+1),flag
c
c**********
c*
c  1) remap from one set of vertical cells to another.
c     method: piecewise linear across each input cell
c             the output is the average of the interpolation
c             profile across each output cell.
c
c  2) input arguments:
c       u,v   - scalar fields in p-layer space
c       t,s,r - scalar fields in p-layer space
c       p     - layer interface depths (non-negative m)
c                 p(   1) is the surface
c                 p(kk+1) is the bathymetry
c       kk    - dimension of a  (number of  input layers)
c       pz    - target interface depths (non-negative m)
c                 pz(k+1) >= pz(k)
c       flag  - data void (land) marker
c       kz    - dimension of az (number of output layers)
c
c  3) output arguments:
c       uz    - scalar field in pz-layer space
c       vz    - scalar field in pz-layer space
c       tz    - scalar field in pz-layer space
c       sz    - scalar field in pz-layer space
c       rz    - scalar field in pz-layer space
c
c  4) except at data voids, must have:
c           p(   1) == zero (surface)
c           p( l+1) >= p(:,:,l)
c           p(kk+1) == bathymetry
c           0 <= pz(k) <= pz(k+1)
c      output layers completely below the bathymetry inherit values
c      from the layer above.
c
c  5) Tim Campbell, Mississippi State University, October 2002.
c*
c**********
c
      real,parameter :: thin=1.e-6  !minimum layer thickness
c
c
      integer k,l,lf
      real    q,qc,zb,zc,zt,uzk,vzk,tzk,szk,rzk
      real    us(kk),vs(kk),ts(kk),ss(kk),rs(kk),pt(kk+1)
c
      if     (r(1).eq.flag) then
        do k= 1,kz
          uz(k) = flag  ! land
          vz(k) = flag  ! land
          tz(k) = flag  ! land
          sz(k) = flag  ! land
          rz(k) = flag  ! land
        enddo
      else
c ---   compute PLM slopes for input layers
        do k=1,kk
          pt(k)=max(p(k+1)-p(k),thin)
        enddo
        call plm(pt, u,v,t,s,r, us,vs,ts,ss,rs, kk)
c ---   compute output layer averages
        lf=1
        zb=pz(1)
        do k= 1,kz
          zt = zb
          zb = pz(k+1)
*         WRITE(6,*) 'k,zt,zb = ',k,zt,zb
          if     (zb-zt.lt.thin .or. zt.ge.p(kk+1)) then
c
c ---       thin or bottomed layer, values taken from layer above
c
            uz(k) = uz(k-1)
            vz(k) = vz(k-1)
            tz(k) = tz(k-1)
            sz(k) = sz(k-1)
            rz(k) = rz(k-1)
          else
c
c           form layer averages.
c
            if     (p(lf).gt.zt) then
              WRITE(6,*) 'bad lf = ',lf
              stop
            endif
            uzk = 0.0
            vzk = 0.0
            tzk = 0.0
            szk = 0.0
            rzk = 0.0
            do l= lf,kk
              if     (p(l).gt.zb) then
*               WRITE(6,*) 'l,lf= ',l,lf,l-1
                lf = l-1
                exit
              elseif (p(l).ge.zt .and. p(l+1).le.zb) then
c
c               the input layer is completely inside the output layer
c
                q   = max(p(l+1)-p(l),thin)/(zb-zt)
                uzk = uzk + q*u(l)
                vzk = vzk + q*v(l)
                tzk = tzk + q*t(l)
                szk = szk + q*s(l)
                rzk = rzk + q*r(l)
*               WRITE(6,*) 'L,q = ',l,q
              else
c
c               the input layer is partially inside the output layer
c               average of linear profile is its center value
c
                q   = max( min(p(l+1),zb)-max(p(l),zt), thin )/(zb-zt)
                zc  = 0.5*(min(p(l+1),zb)+max(p(l),zt))
                qc  = (zc-p(l))/pt(l) - 0.5
                uzk = uzk + q*(u(l) + qc*us(l))
                vzk = vzk + q*(v(l) + qc*vs(l))
                tzk = tzk + q*(t(l) + qc*ts(l))
                szk = szk + q*(s(l) + qc*ss(l))
                rzk = rzk + q*(r(l) + qc*rs(l))
*               WRITE(6,*) 'l,q,qc = ',l,q,qc
              endif
            enddo !l
            uz(k) = uzk
            vz(k) = vzk
            tz(k) = tzk
            sz(k) = szk
            rz(k) = rzk
          endif
        enddo !k
      endif
      return
      end subroutine layer2zi_plm

      subroutine plm(pt, u, v, t, s, r,
     &                   us,vs,ts,ss,rs,kk)
      implicit none
c
      integer kk
      real    u(kk),v(kk),t(kk),s(kk),r(kk),pt(kk),
     &        us(kk),vs(kk),ts(kk),ss(kk),rs(kk)
c
c**********
c*
c  1) generate a monotonic PLM interpolation of a layered field
c
c  2) input arguments:
c       pt    - layer interface thicknesses (non-zero)
c       u,v   - scalar fields in layer space
c       t,s,r - scalar fields in layer space
c       kk    - dimension of a  (number of layers)
c
c  3) output arguments:
c       us    - scalar field slopes for PLM interpolation
c       vs    - scalar field slopes for PLM interpolation
c       ts    - scalar field slopes for PLM interpolation
c       ss    - scalar field slopes for PLM interpolation
c       rs    - scalar field slopes for PLM interpolation
c
c  4) except at data voids, must have:
c           p(   1) == zero (surface)
c           p( l+1) >= p(:,:,l)
c           p(kk+1) == bathymetry
c
c  5) Tim Campbell, Mississippi State University, September 2002.
c*
c**********
c
      integer l
      real    ql(kk),qc(kk),qr(kk)
c
      !compute grid spacing ratios for slope computations
      ql(1)=0.0
      qc(1)=0.0
      qr(1)=0.0
      do l=2,kk-1
        ql(l)=2.0*pt(l)/(pt(l-1)+pt(l))
        qc(l)=2.0*pt(l)/(pt(l-1)+2.0*pt(l)+pt(l+1))
        qr(l)=2.0*pt(l)/(pt(l)+pt(l+1))
      enddo
      ql(kk)=0.0
      qc(kk)=0.0
      qr(kk)=0.0
      !compute normalized layer slopes
      call slope(ql,qc,qr,u,us,kk)
      call slope(ql,qc,qr,v,vs,kk)
      call slope(ql,qc,qr,t,ts,kk)
      call slope(ql,qc,qr,s,ss,kk)
      call slope(ql,qc,qr,r,rs,kk)
      return
      end subroutine plm

      subroutine slope(rl,rc,rr,a,s,n)
      implicit none
c
      integer,intent(in)  :: n
      real,   intent(in)  :: rl(n),rc(n),rr(n),a(n)
      real,   intent(out) :: s(n)
c
c**********
c*
c  1) generate slopes for monotonic piecewise linear distribution
c
c  2) input arguments:
c       rl   - left grid spacing ratio
c       rc   - center grid spacing ratio
c       rr   - right grid spacing ratio
c       a    - scalar field zone averages
c       n    - number of zones
c
c  3) output arguments:
c       s    - zone slopes
c
c  4) Tim Campbell, Mississippi State University, September 2002.
c*
c**********
c
      integer,parameter :: ic=2, im=1, imax=100
      real,parameter :: fracmin=1e-6, dfac=0.5
c
      integer i,j
      real    sl,sc,sr
      real    dnp,dnn,dl,dr,ds,frac
c
c Compute zone slopes
c Campbell Eq(15) -- nonuniform grid
c
      s(1)=0.0
      do j=2,n-1
        sl=rl(j)*(a(j)-a(j-1))
        sr=rr(j)*(a(j+1)-a(j))
        if (sl*sr.gt.0.) then
          s(j)=sign(min(abs(sl),abs(sr)),sl)
        else
          s(j)=0.0
        endif
      enddo
      s(n)=0.0
c
c Minimize discontinuities between zones
c Apply single pass discontinuity minimization: Campbell Eq(19)
c
      do j=2,n-1
        if(s(j).ne.0.0) then
          dl=-0.5*(s(j)+s(j-1))+a(j)-a(j-1)
          dr=-0.5*(s(j+1)+s(j))+a(j+1)-a(j)
          ds=sign(min(abs(dl),abs(dr)),dl)
          s(j)=s(j)+2.0*ds
        endif
      enddo
      return
      end subroutine slope

      subroutine blkinr(rvar,cvar,cfmt)
      implicit none
c
      real      rvar
      character cvar*6,cfmt*(*)
c
c     read in one real value, on unit 12
c
      character*6 cvarin
c
      read(12,*) rvar,cvarin
*     write(6,cfmt) cvarin,rvar
c
      if     (cvar.ne.cvarin) then
        write(6,*) 
        write(6,*) 'error in blkinr - input ',cvarin,
     &                      ' but should be ',cvar
        write(6,*) 
        stop
      endif
      return
      end

      subroutine blkini(ivar,cvar)
      implicit none
c
      integer     ivar
      character*6 cvar
c
c     read in one integer value, on unit 12
c
      character*6 cvarin
c
      read(12,*) ivar,cvarin
*     write(6,6000) cvarin,ivar
c
      if     (cvar.ne.cvarin) then
        write(6,*) 
        write(6,*) 'error in blkini - input ',cvarin,
     &                      ' but should be ',cvar
        write(6,*) 
        stop
      endif
      return
 6000 format('# ',a6,' =',i6)
      end
