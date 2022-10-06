      PROGRAM PROFILE
      IMPLICIT NONE
C
C  hycom_profile2z - Usage:  hycom_profile2z archv.txt z.txt archz.txt [itype [bot [plat]]]
C
C                 converts an HYCOM isopycnal text profile file to
C                 a z-level text profile file.
C
C   archv.txt is assumed to be an HYCOM archive text profile file
C   z.txt     is assumed to contain a list on z-depths, one per line
C   archz.txt will be the output z-level text profile file
C   itype     is the interpolation type (default 1)
C                =-4; piecewise quadratic across cell, WENO-like, new extrema
C                =-3; piecewise quadratic across cell, WENO-like
C                =-2; piecewise quadratic across cell
C                =-1; piecewise linear    across cell
C                = 0; piecewise constant  across cell
C                = 1; linear interpolation between cell centers
C                = 3; piecewise cubic hermite interpolation between cell centers
C   bot       ignore layers within bot of the bottom (default 0.0)
C                < 0; split lowest layer into two thk+bot,-bot
C   plat      latitude of the location (write out temperature)
C
C  this version for "serial" Unix systems.
C
C  Tim Campbell,       Mississippi State University
C  Alan J. Wallcraft,  Naval Research Laboratory,  Mar. 2002 (Aug. 2007).
C
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      CHARACTER*240 CFILEA,CFILEB,CFILEC
      CHARACTER*240 CLINE
      REAL          THK,DUM(5),FLAG,BOT,PLAT,PZK,PZKP,TEMP,ZZ(9999)
      INTEGER       IOS,ITYPE,K,KA,KK,KDM,KZ
      INTEGER       I,KMAX
      LOGICAL       LDEBUG
C
      REAL, ALLOCATABLE :: SI(:,:),P(:),SZ(:,:)
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      CALL GETARG(0,CLINE)
      K = LEN_TRIM(CLINE)
      LDEBUG = CLINE(K-5:K) .EQ. "_debug"
C
      IF     (NARG.EQ.6) THEN
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CFILEB)
        CALL GETARG(3,CFILEC)
        CALL GETARG(4,CLINE)
        READ(CLINE,*) ITYPE
        CALL GETARG(5,CLINE)
        READ(CLINE,*) BOT
        CALL GETARG(6,CLINE)
        READ(CLINE,*) PLAT
      ELSEIF (NARG.EQ.5) THEN
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CFILEB)
        CALL GETARG(3,CFILEC)
        CALL GETARG(4,CLINE)
        READ(CLINE,*) ITYPE
        CALL GETARG(5,CLINE)
        READ(CLINE,*) BOT
        PLAT  = -99.99
      ELSEIF (NARG.EQ.4) THEN
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CFILEB)
        CALL GETARG(3,CFILEC)
        CALL GETARG(4,CLINE)
        READ(CLINE,*) ITYPE
        BOT   =  0.0
        PLAT  = -99.99
      ELSEIF (NARG.EQ.3) THEN
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CFILEB)
        CALL GETARG(3,CFILEC)
        ITYPE =  1
        BOT   =  0.0
        PLAT  = -99.99
      ELSE
        WRITE(6,'(2a)')
     +  'Usage: hycom_profile2z archv.txt z.txt archz.txt',
     +  ' [itype [bot [plat]]]'
        CALL EXIT(1)
      ENDIF
C
C     OPEN ALL FILES.
C
      OPEN(UNIT=11, FILE=CFILEA, FORM='FORMATTED', STATUS='OLD',
     +     IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',TRIM(CFILEA)
        WRITE(6,*) 'ios   = ',ios
        CALL EXIT(3)
      ENDIF
      OPEN(UNIT=12, FILE=CFILEB, FORM='FORMATTED', STATUS='OLD',
     +     IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',TRIM(CFILEB)
        WRITE(6,*) 'ios   = ',ios
        CALL EXIT(4)
      ENDIF
      OPEN(UNIT=21, FILE=CFILEC, FORM='FORMATTED', STATUS='NEW',
     +     IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',TRIM(CFILEC)
        WRITE(6,*) 'ios   = ',ios
        CALL EXIT(5)
      ENDIF
C
C     COPY PROFILE HEADER TO OUTPUT.
C     HEADER CONTAINS MULTIPLE PAIRS OF LINES STARTING WITH "##"
C     THEN ONE LINE STARTING WITH "#  k" (WHICH IS DISCARDED)
C
      DO 
        READ( 11,'(a)') CLINE
        IF     (CLINE(1:2).NE.'##') THEN
          EXIT
        ENDIF
        WRITE(21,'(a)') TRIM(CLINE)
        READ( 11,'(a)')      CLINE
        WRITE(21,'(a)') TRIM(CLINE)
      ENDDO
C
C     READ THE ISOPYCNAL PROFILE, TO GET KDM.
C
      KDM   = -1
      DO K= 1,99
        READ(11,'(a)',IOSTAT=IOS) CLINE
        IF     (IOS.NE.0) THEN
          IF     (K.NE.KDM+1) THEN
            WRITE(6,*) 'Error: inconsistent input profile'
            CALL EXIT(6)
          ENDIF
          EXIT
        ENDIF
        READ(CLINE,*) KDM
      ENDDO
C
      IF     (LDEBUG) THEN
        WRITE(6,*) 'KDM = ',KDM
      ENDIF
C
C     RE-READ THE ISOPYCNAL PROFILE.
C
      ALLOCATE( P(KDM+1), SI(KDM,5) )
C
      REWIND(11)
      DO 
        READ( 11,'(a)') CLINE
        IF     (CLINE(1:2).NE.'##') THEN
          EXIT
        ENDIF
        READ( 11,'(a)') CLINE
      ENDDO
      P(1) =  0.0
      DO K= 1,KDM
        READ(11,'(a)',IOSTAT=IOS) CLINE
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'Error: inconsistent input profile'
          CALL EXIT(6)
        ENDIF
        READ(CLINE,*) KK,SI(K,1),SI(K,2),SI(K,3),SI(K,4),SI(K,5),THK
        IF     (KK.NE.K) THEN
          WRITE(6,*) 'Error: inconsistent input profile - k,kin =',
     &                K,KK
          CALL EXIT(6)
        ENDIF
        READ(CLINE,*) KK,SI(K,1),SI(K,2),SI(K,3),SI(K,4),SI(K,5),THK
        P(K+1) = P(K) + THK
        IF     (THK.EQ.0.0) THEN
          DO KK= 1,5
            SI(K,KK)=SI(K-1,KK)
          ENDDO !kk
        ENDIF
C
        IF     (LDEBUG) THEN
          WRITE(6,'(a,i3,7f10.2)') 'K,SI,THK = ',K,SI(K,1:5),THK,P(K+1)
        ENDIF
      ENDDO
      CLOSE(11)
C
      IF     (BOT.GT.0.0) THEN
        DO K= 2,KDM
          IF      (P(K).GE.P(KDM+1)-BOT) THEN
            DO KK= 1,5
              SI(K,KK)=SI(K-1,KK)
            ENDDO !kk
            IF     (LDEBUG) THEN
              WRITE(6,'(a,i3,6f10.2)') 'k,si,thk = ',K,SI(K,1:5),
     &                                 P(K+1)-P(K)
            ENDIF
          ENDIF
        ENDDO
      ELSEIF (BOT.LT.0.0) THEN
        BOT = ABS(BOT)
        DO K= 1,KDM-1
          IF     (P(K+1).EQ.P(KDM+1)) THEN
            IF     (P(K+1)-P(K).GT.3.0*BOT) THEN
              P(K+1) = P(K+1)-BOT
              SI(K,  1) = (P(K+1)-P(K))/(P(K+1)-P(K)-BOT)*SI(K,1)
              SI(K,  2) = (P(K+1)-P(K))/(P(K+1)-P(K)-BOT)*SI(K,2)
              SI(K+1,1:2) = 0.0        !zero velocity at bottom
              SI(K+1,3:5) = SI(K,3:5)  !no change in scalars
              IF     (LDEBUG) THEN
                WRITE(6,'(a,i3,6f10.2)') 'K,si,thk = ',K,  SI(K,1:5),
     &                                   P(K+1)-P(K)
                WRITE(6,'(a,i3,6f10.2)') 'K,si,thk = ',K+1,SI(K+1,1:5),
     &                                   P(K+2)-P(K+1)
              ENDIF
              DO KA= K+2,KDM
                SI(KA,1:5)=SI(KA-1,1:5)
                IF     (LDEBUG) THEN
                  WRITE(6,'(a,i3,6f10.2)') 'K,si,thk = ',KA,SI(KA,1:5),
     &                                     P(KA+1)-P(KA)
                ENDIF
              ENDDO !ka
              EXIT
            ELSE
              EXIT
            ENDIF !>3bot
          ENDIF !==p(kdm+1)
        ENDDO !k
      ENDIF !bot
C
C     READ THE Z FILE.
C
      DO K= 1,9999
        READ(12,'(a)',IOSTAT=IOS) CLINE
        IF     (IOS.NE.0) THEN
          KZ = K-1
          EXIT
        ENDIF
        READ(CLINE,*) ZZ(K)
        ZZ(K) =MIN( ZZ(K), P(KDM+1) )
        IF     (LDEBUG) THEN
          WRITE(6,'(a,i3,f10.2)') 'k,zz = ',K,ZZ(K)
        ENDIF
      ENDDO
      CLOSE(12)
C
      ALLOCATE( SZ(KZ,5) )
C
C     INTERPOLATE TO Z
C
      FLAG = 0.0
      IF     (ITYPE.EQ.-4) THEN
        CALL LAYER2Z_WNX(SI,P, KDM,5,
     &                   SZ,ZZ,KZ,   FLAG)
      ELSEIF (ITYPE.EQ.-3) THEN
        CALL LAYER2Z_WNO(SI,P, KDM,5,
     &                   SZ,ZZ,KZ,   FLAG)
      ELSEIF (ITYPE.EQ.-2) THEN
        CALL LAYER2Z_PPM(SI,P, KDM,5,
     &                   SZ,ZZ,KZ,   FLAG)
      ELSEIF (ITYPE.EQ.-1) THEN
        CALL LAYER2Z_PLM(SI,P, KDM,5,
     &                   SZ,ZZ,KZ,   FLAG)
      ELSEIF (ITYPE.EQ. 0) THEN
        CALL LAYER2Z_PCM(SI,P, KDM,5,
     &                   SZ,ZZ,KZ,   FLAG)
      ELSEIF (ITYPE.EQ. 1) THEN
        CALL LAYER2Z_LIN(SI,P, KDM,5,
     &                   SZ,ZZ,KZ,   FLAG)
      ELSEIF (ITYPE.EQ. 3) THEN
        CALL LAYER2Z_PCH(SI,P, KDM,5,
     &                   SZ,ZZ,KZ,   FLAG)
      ELSE
        WRITE(6,*)
     +  'Usage: hycom_profile2z archv.txt z.txt archz.txt [itype [bot]]'
        WRITE(6,*) 'unsupported itype value'
        CALL EXIT(1)
      ENDIF
C
C     OUTPUT.
C
      IF     (PLAT.LT.-90.0) THEN
        WRITE(21,'(a,a)')
     &  '#   k',
     &  '    utot    vtot  p.temp    saln  p.dens    thkns      dpth'
        PZKP = 0.0
        DO K= 1,KZ
          PZK  = PZKP
          IF     (K.LT.KZ) THEN
            PZKP = 0.5*(ZZ(K)+ZZ(K+1))
          ELSE
            PZKP = MAX(ZZ(KZ),P(KDM+1))
          ENDIF
          WRITE(21,'(i5,2f8.2,3f8.4,f9.3,f10.3)')
     &      K,
     &      SZ(K,1),             !cm/s
     &      SZ(K,2),             !cm/s
     &      SZ(K,3),             !degC
     &      SZ(K,4),             !psu
     &      SZ(K,5),             !SigmaT
     &      PZKP-PZK,            !m
     &      ZZ(K)                !m
          IF     (ZZ(K).EQ.P(KDM+1)) THEN
            EXIT !hit the bottom
          ENDIF
        ENDDO !k
      ELSE
        WRITE(21,'(a,a,a)')
     &  '#   k',
     &  '    utot    vtot  p.temp    saln  p.dens    thkns      dpth',
     &  '    temp'
        PZKP = 0.0
        DO K= 1,KZ
          PZK  = PZKP
          IF     (K.LT.KZ) THEN
            PZKP = 0.5*(ZZ(K)+ZZ(K+1))
          ELSE
            PZKP = MAX(ZZ(KZ),P(KDM+1))
          ENDIF
          CALL PT2T(TEMP,SZ(K,3),SZ(K,4),ZZ(K),PLAT)
          WRITE(21,'(i4,1x,2f8.2,3f8.4,f9.3,f10.3,f8.4)')
     &      K,
     &      SZ(K,1),             !cm/s
     &      SZ(K,2),             !cm/s
     &      SZ(K,3),             !degC
     &      SZ(K,4),             !psu
     &      SZ(K,5),             !SigmaT
     &      PZKP-PZK,            !m
     &      ZZ(K),               !m
     &      TEMP                 !degC
          IF     (ZZ(K).EQ.P(KDM+1)) THEN
            EXIT !hit the bottom
          ENDIF
        ENDDO !k
      ENDIF !plat
      END

      subroutine layer2z_pch(si,p,kk,ks,
     &                       sz,z,kz,   flag)
      implicit none
c
      integer kk,ks,kz
      real    si(kk,ks),p(kk+1),
     &        sz(kz,ks),z(kz),flag
c
c**********
c*
c  1) interpolate a set of layered fields to fixed z depths.
c     method: piecewise cubic hermite interpolation (PCHIP) between cell centers
c
c  2) input arguments:
c       si    - scalar fields in layer space
c       p     - layer interface depths (non-negative m)
c                 p(   1) is the surface
c                 p(kk+1) is the bathymetry
c       kk    - dimension of si (number of layers)
c       ks    - dimension of si (number of fields)
c       z     - target z-level  depths (non-negative m)
c       flag  - data void (land) marker
c       kz    - dimension of sz (number of levels)
c
c  3) output arguments:
c       sz    - scalar fields in z-space
c
c  4) except at data voids, must have:
c           p(   1) == zero (surface)
c           p( l+1) >= p(:,:,l)
c           p(kk+1) == bathymetry
c           0 <= z(k) <= z(k+1)
c     note that z(k) > p(kk+1) implies that az(k)=flag,
c      since the z-level is then below the bathymetry.
c
c  5) Alan J. Wallcraft, COAPS, January 2018.
c*
c**********
c
      integer i,k,l,lf,ngood
      real    oldz(kk+1),olds(kk+1)
c
      if     (si(1,ks).eq.flag) then
        do k= 1,kz
          do i= 1,ks
            sz(k,i) = flag  !land
          enddo !i
        enddo !k
      else
        do k= 1,kk
          oldz(k) = 0.5*(p(k+1)+p(k))
          if     (p(k+1).eq.p(kk+1)) then  !.true. for k==kk
            ngood = k
            exit
          endif
        enddo !k
        oldz(ngood+1) = p(kk+1)
c
c ---   loop through scalar fields
        do i= 1,ks
          do k= 1,ngood
            olds(k) = si(k,i)
          enddo !k
          olds(ngood+1) = olds(ngood)
          call pchip(ngood+1,oldz,olds, flag, kz,z,sz(1,i))
        enddo !i
      endif
      return
      end

      subroutine pchip(nin,xin,yin, flag, nout,xout,yout)
! 20040924 Rowley, C.
!   modified somewhat from pchiptest.f from Mike Carnes Sep 2004
!   AJW: removed gapsiz added flag
      implicit none
c
      integer, intent(in)  :: nin
      real,    intent(in)  :: xin(nin),yin(nin)
      real,    intent(in)  :: flag
      integer, intent(in)  :: nout
      real,    intent(in)  :: xout(nout)
      real,    intent(out) :: yout(nout)
c
      integer :: k,kk
      real    :: del(nin),b(nin),c(nin),d(nin)
      logical :: keepgoing
      real    :: hh,s
!
! PCHIP  Piecewise Cubic Hermite Interpolating Polynomial.
! X is a row or column vector.  Y is a row or column vector of the same
! length as X, or a matrix with length(X) columns.
! YI = PCHIP(X,Y,XI) evaluates an interpolant at the elements of XI.
! PP = PCHIP(X,Y) returns a piecewise polynomial structure for use by PPVAL.
!
! The PCHIP interpolating function, P(x), satisfies:
!   On each subinterval,  x(k) <= x <= x(k+1),  P(x) is the cubic Hermite
!     interpolant to the given values and certain slopes at the two endpoints.
!   Therefore, P(x) interpolates y, i.e., P(x(j)) = y(j), and the first
!     derivative, P(x), is continuous, but P''(x) is probably not continuous;
!     there may be jumps at the x(j).
!   The slopes at the x(j) are chosen in such a way that P(x) is
!     "shape preserving" and "respects monotonicity". This means that,
!     on intervals where the data are monotonic, so is P(x);
!     at points where the data have a local extremum, so does P(x).
!
! References:
!   F. N. Fritsch and R. E. Carlson, "Monotone Piecewise Cubic
!   Interpolation", SIAM J. Numerical Analysis 17, 1980, 238-246.
!   David Kahaner, Cleve Moler and Stephen Nash, Numerical Methods
!   and Software, Prentice Hall, 1988.
!
! Find indices of subintervals, x(k) <= u < x(k+1).
! there must be at least two x values
!
c      kk=1
c      keepgoing=.true.
c      do k=1,nout
c        if(u(k).lt.x(1)) then
c          kk=1
c        else
c          dowhile((x(kk+1).lt.u(k)).and.keepgoing) then
c            if(kk.lt.nin-1) then
c              kk=kk+1
c            else
c              keepgoing=.false.
c            endif
c          enddo
c        endif
c        ku(k)=kk 
c        s(k) = u(k) - x(ku(k));
c      enddo

! Compute slopes and other coefficients.
!   
c      write(10,'(a)') 'k,del(k)'
      do k=1,nin-1
        del(k)=(yin(k+1)-yin(k))/(xin(k+1)-xin(k));  ! write(10,*) k,del(k)
      enddo
c
      call pchipslopes(nin,xin,yin,del,d);           ! write(10,'(a)') 'k,c,b,d'
c
      do k=1,nin-1 
        hh=xin(k+1)-xin(k)
        c(k)=(3.*del(k)-2.*d(k)-d(k+1))/hh
        b(k)=(d(k)-2.*del(k)+d(k+1))/(hh*hh);        ! write(10,*) k,c(k),b(k),d(k)
      enddo

! Evaluate interpolant.
c     write(10,'(a)') 'evaluate interpolant'
c
      kk=1
      do k=1,nout
        if(xout(k).lt.xin(1)) then
          yout(k)=yin(1)
        elseif(xout(k).ge.xin(1).and.xout(k).le.xin(nin)) then
          if(xout(k).eq.xin(nin)) then
            yout(k)=yin(nin)
          else
            do while(xin(kk+1).lt.xout(k))
              kk=kk+1
            enddo 
            if(abs(xout(k)-xin(kk)).lt.0.001) then
              yout(k)=yin(kk)
            else
              s=xout(k)-xin(kk)
              yout(k)=yin(kk)+s*(d(kk)+s*(c(kk)+s*b(kk)));
              ! write(10,*) k,kk,s,y(kk),d(kk),c(kk),b(kk)
            endif
          endif
        else
          yout(k)=flag
        endif
      enddo
c
      return
      end subroutine pchip

      subroutine pchipslopes(n,x,y,del,d)
      implicit none
c
      integer, intent(in)    :: n
      real,    intent(in)    :: x(n),y(n),del(n)
      real,    intent(inout) :: d(n)
c      
      real    :: h(n)
      real    :: hs,w1,w2,dmax,dmin
      integer :: k
      integer :: isign1,isign2
!
! PCHIPSLOPES  Derivative values for Piecewise Hermite Cubic Interpolation.
! d = pchipslopes(x,y,del) computes the first derivatives, d(k) = P'(x(k))'.
!
      do k=1,n
        d(k)=0.
      enddo
!
!     Special case n=2, use linear interpolation.
      if(n.eq.2) then
        d(1) = del(1);
        d(2) = del(1);
        return
      endif
!
!  Slopes at interior points.
!  d(k) = weighted average of del(k-1) and del(k) when they have the same sign.
!  d(k) = 0 when del(k-1) and del(k) have opposites signs or either is zero.
!
      do k=1,n
        d(k)=0.
      enddo
c
      do k=1,n-1
        h(k)=x(k+1)-x(k)
      enddo
c
      do k=1,n-2
        if ((del(k).gt.0..and.del(k+1).gt.0.).or.
     &      (del(k).lt.0..and.del(k+1).lt.0.)    ) then
          hs=h(k)+h(k+1)
          w1=(h(k)+hs)/(3.*hs)
          w2=(hs+h(k+1))/(3.*hs)
          dmax=max(abs(del(k)),abs(del(k+1)))
          dmin=min(abs(del(k)),abs(del(k+1)))
          d(k+1)=dmin/(w1*(del(k)/dmax)+w2*(del(k+1)/dmax))
        endif
      enddo
!
!  Slopes at end points.
!  Set d(1) and d(n) via non-centered, shape-preserving three-point formulae.
      d(1)=((2.*h(1)+h(2))*del(1)-h(1)*del(2))/(h(1)+h(2))
c
      isign1=0
      if(d(1).gt.0..and.del(1).gt.0.) then
        isign1=1
      elseif(d(1).lt.0..and.del(1).lt.0.) then
        isign1=1
      elseif(d(1).eq.0..and.del(1).eq.0.) then
        isign1=1
      endif
c
      isign2=0
      if(del(1).gt.0..and.del(2).gt.0.) then
        isign2=1
      elseif(del(1).lt.0..and.del(2).lt.0.) then
        isign2=1
      elseif(del(1).eq.0..and.del(2).eq.0.) then
        isign2=1
      endif
c
      if(isign1.eq.0) then
        d(1)=0.
      elseif((isign2.eq.0).and.(abs(d(1)).gt.abs(3.*del(1)))) then
        d(1)=3.*del(1)
      endif
c
      d(n)=((2.*h(n-1)+h(n-2))*del(n-1)-h(n-1)*del(n-2))/(h(n-1)+h(n-2))
c      
      isign1=0
      if(d(n).gt.0..and.del(n-1).gt.0.) then
        isign1=1
      elseif(d(n).lt.0..and.del(n-1).lt.0.) then
        isign1=1
      elseif(d(n).eq.0..and.del(n-1).eq.0.) then
        isign1=1
      endif
c
      isign2=0
      if(del(n-1).gt.0..and.del(n-2).gt.0.) then
        isign2=1
      elseif(del(n-1).lt.0..and.del(n-2).lt.0.) then
        isign2=1
      elseif(del(n-1).eq.0..and.del(n-2).eq.0.) then
        isign2=1
      endif
c
      if(isign1.eq.0) then
        d(n)=0.
      elseif((isign2.eq.0).and.(abs(d(n)).gt.abs(3.*del(n-1)))) then
        d(n)=3.*del(n-1)
      endif 
c
      return
      end subroutine pchipslopes

      subroutine layer2z_lin(si,p,kk,ks,
     &                       sz,z,kz,   flag)
      implicit none
c
      integer kk,ks,kz
      real    si(kk,ks),p(kk+1),
     &        sz(kz,ks),z(kz),flag
c
c**********
c*
c  1) interpolate a set of layered fields to fixed z depths.
c     method: linear interpolation to cell centers
c
c  2) input arguments:
c       si    - scalar fields in layer space
c       p     - layer interface depths (non-negative m)
c                 p(   1) is the surface
c                 p(kk+1) is the bathymetry
c       kk    - dimension of si (number of layers)
c       ks    - dimension of si (number of fields)
c       z     - target z-level  depths (non-negative m)
c       flag  - data void (land) marker
c       kz    - dimension of sz (number of levels)
c
c  3) output arguments:
c       sz    - scalar fields in z-space
c
c  4) except at data voids, must have:
c           p(   1) == zero (surface)
c           p( l+1) >= p(:,:,l)
c           p(kk+1) == bathymetry
c           0 <= z(k) <= z(k+1)
c     note that z(k) > p(kk+1) implies that az(k)=flag,
c      since the z-level is then below the bathymetry.
c
c  5) Alan J. Wallcraft, Naval Research Laboratory, February 2002.
c*
c**********
c
      integer i,k,l,lf
      real    q,zk,z0,zm,zp
c
      if     (si(1,ks).eq.flag) then
        do k= 1,kz
          do i= 1,ks
            sz(k,i) = flag  !land
          enddo !i
        enddo !k
      else
        lf=1
        do k= 1,kz
          zk=z(k)
          do l= lf,kk
            if     (p(l).le.zk .and. p(l+1).ge.zk) then
c
c             z(k) is in layer l.
c             linear interpolation between layer centers
c
              z0 = 0.5*(p(l)+p(l+1))
              if     (zk.le.z0) then
c
c               z(k) is in the upper half of the layer
c
                if     (l.eq.1) then
                  do i= 1,ks
                    sz(k,i) = si(1,i)
                  enddo !i
                else
                  zm = 0.5*(p(l-1)+p(l))
                  q  = (z0 - zk)/(z0 - zm)
                  do i= 1,ks
                    sz(k,i) = q*si(l-1,i) + (1.0-q)*si(l,i)
                  enddo !i
                endif
              else
c
c               z(k) is in the lower half of the layer
c
                if     (p(l+1).eq.p(kk+1)) then
                  do i= 1,ks
                    sz(k,i) = si(kk,i)
                  enddo !i
                else
                  zp = 0.5*(p(l+1)+p(l+2))
                  q  = (zk - z0)/(zp - z0)
                  do i= 1,ks
                    sz(k,i) = q*si(l+1,i) + (1.0-q)*si(l,i)
                  enddo !i
                endif
              endif
              lf = l  ! z monotonic increasing, so z(k+1) in layers l:kk
              exit
            elseif (l.eq.kk) then
              do i= 1,ks
                sz(k,i) = flag  ! below the bottom
              enddo !i
              lf = l
              exit
            endif
          enddo !l
        enddo !k
      endif
      return
      end

      subroutine layer2z_pcm(si,p,kk,ks,
     &                       sz,z,kz,   flag)
      implicit none
c
      integer kk,ks,kz
      real    si(kk,ks),p(kk+1),
     &        sz(kz,ks),z(kz),flag
c
c**********
c*
c  1) interpolate a set of layered field to fixed z depths.
c     method: piecewise constant across each cell
c             i.e. sample the layer spaning each depth
c
c  2) input arguments:
c       si    - scalar fields in layer space
c       p     - layer interface depths (non-negative m)
c                 p(   1) is the surface
c                 p(kk+1) is the bathymetry
c       kk    - dimension of si (number of layers)
c       ks    - dimension of si (number of fields)
c       z     - target z-level  depths (non-negative m)
c       flag  - data void (land) marker
c       kz    - dimension of sz (number of levels)
c
c  3) output arguments:
c       sz    - scalar fields in z-space
c
c  4) except at data voids, must have:
c           p(   1) == zero (surface)
c           p( l+1) >= p(:,:,l)
c           p(kk+1) == bathymetry
c           0 <= z(k) <= z(k+1)
c     note that z(k) > p(kk+1) implies that az(k)=flag,
c      since the z-level is then below the bathymetry.
c
c  5) Alan J. Wallcraft, Naval Research Laboratory, February 2002.
c*
c**********
c
      integer i,k,l,lf
      real    zk
c
      if     (si(1,ks).eq.flag) then
        do k= 1,kz
          do i= 1,ks
           sz(k,i) = flag  !land
          enddo !i
        enddo !k
      else
        lf=1
        do k= 1,kz
          zk=z(k)
          do l= lf,kk
            if     (p(l).le.zk .and. p(l+1).ge.zk) then
c
c             z(k) is in layer l.
c             return cell average.
c
              do i= 1,ks
                sz(k,i) = si(l,i)
              enddo !i
              lf = l  ! z monotonic increasing, so z(k+1) in layers l:kk
              exit
            elseif (l.eq.kk) then
              do i= 1,ks
                sz(k,i) = flag  ! below the bottom
              enddo !i
              lf = l
              exit
            endif
          enddo !l
        enddo !k
      endif
      return
      end

      subroutine layer2z_plm(si,p,kk,ks,
     &                       sz,z,kz,   flag)
      implicit none
c
      integer kk,ks,kz
      real    si(kk,ks),p(kk+1),
     &        sz(kz,ks),z(kz),flag
c
c**********
c*
c  1) interpolate a set of layered field to fixed z depths.
c     method: piecewise linear across each cell
c
c  2) input arguments:
c       si    - scalar fields in layer space
c       p     - layer interface depths (non-negative m)
c                 p(   1) is the surface
c                 p(kk+1) is the bathymetry
c       kk    - dimension of si (number of layers)
c       ks    - dimension of si (number of fields)
c       z     - target z-level  depths (non-negative m)
c       flag  - data void (land) marker
c       kz    - dimension of sz (number of levels)
c
c  3) output arguments:
c       sz    - scalar fields in z-space
c
c  4) except at data voids, must have:
c           p(   1) == zero (surface)
c           p( l+1) >= p(:,:,l)
c           p(kk+1) == bathymetry
c           0 <= z(k) <= z(k+1)
c     note that z(k) > p(kk+1) implies that az(k)=flag,
c      since the z-level is then below the bathymetry.
c
c  5) Alan J. Wallcraft, Naval Research Laboratory, February 2002.
c*
c**********
c
      real, parameter :: thin=1.e-6  !minimum layer thickness
c
      integer i,k,l,lf
      real    q,zk
      real    sis(kk,ks),pt(kk+1)
c
      if     (si(1,ks).eq.flag) then
        do k= 1,kz
          do i= 1,ks
            sz(k,i) = flag  !land
          enddo !i
        enddo
      else
c ---   compute PLM slopes for input layers
        do k=1,kk
          pt(k)=max(p(k+1)-p(k),thin)
        enddo
        call plm(pt,si,sis,kk,ks)
c
        lf=1
        do k= 1,kz
          zk=z(k)
          do l= lf,kk
            if     (p(l).le.zk .and. p(l+1).ge.zk) then
c
c             z(k) is in layer l, sample the linear profile at zk.
c
              q = (zk-p(l))/pt(l) - 0.5
              do i= 1,ks
                sz(k,i) = si(l,i) + q*sis(l,i)
              enddo !i
              lf = l  ! z monotonic increasing, so z(k+1) in layers l:kk
              exit
            elseif (l.eq.kk) then
              do i= 1,ks
                sz(k,i) = flag  ! below the bottom
              enddo !i
              lf = l
              exit
            endif
          enddo !l
        enddo !k
      endif
      return
      end

      subroutine plm(pt, s,ss,ki,ks)
      implicit none
c
      integer ki,ks
      real    pt(ki+1),s(ki,ks),ss(ki,ks)
c
c**********
c*
c  1) generate a monotonic PLM interpolation of a layered field
c
c  2) input arguments:
c       pt    - layer interface thicknesses (non-zero)
c       s     - scalar fields in layer space
c       ki    - 1st dimension of s (number of layers)
c       ks    - 2nd dimension of s (number of fields)
c
c  3) output arguments:
c       ss    - scalar field slopes for PLM interpolation
c
c  4) except at data voids, must have:
c           pi(   1) == zero (surface)
c           pi( l+1) >= pi(:,:,l)
c           pi(ki+1) == bathymetry
c
c  5) Tim Campbell, Mississippi State University, September 2002.
c*
c**********
c
      integer l
      real    ql(ki),qc(ki),qr(ki)
c
      !compute grid spacing ratios for slope computations
      ql(1)=0.0
      qc(1)=0.0
      qr(1)=0.0
      do l=2,ki-1
        ql(l)=2.0*pt(l)/(pt(l-1)+pt(l))
        qc(l)=2.0*pt(l)/(pt(l-1)+2.0*pt(l)+pt(l+1))
        qr(l)=2.0*pt(l)/(pt(l)+pt(l+1))
      enddo
      ql(ki)=0.0
      qc(ki)=0.0
      qr(ki)=0.0
      !compute normalized layer slopes
      do l=1,ks
        call slope(ql,qc,qr,s(1,l),ss(1,l),ki)
      enddo
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

      subroutine layer2z_ppm(si,p,kk,ks,
     &                       sz,z,kz,   flag)
      implicit none
c
      integer kk,ks,kz
      real    si(kk,ks),p(kk+1),
     &        sz(kz,ks),z(kz),flag
c
c**********
c*
c  1) interpolate a set of layered field to fixed z depths.
c     method: piecewise parabolic method across each cell
c
c  2) input arguments:
c       si    - scalar fields in layer space
c       p     - layer interface depths (non-negative m)
c                 p(   1) is the surface
c                 p(kk+1) is the bathymetry
c       kk    - dimension of si (number of layers)
c       ks    - dimension of si (number of fields)
c       z     - target z-level  depths (non-negative m)
c       flag  - data void (land) marker
c       kz    - dimension of az (number of levels)
c
c  3) output arguments:
c       sz    - scalar fields in z-space
c
c  4) except at data voids, must have:
c           p(   1) == zero (surface)
c           p( l+1) >= p(:,:,l)
c           p(kk+1) == bathymetry
c           0 <= z(k) <= z(k+1)
c     note that z(k) > p(kk+1) implies that az(k)=flag,
c      since the z-level is then below the bathymetry.
c
c  5) Alan J. Wallcraft, Naval Research Laboratory, February 2002.
c*
c**********
c
      real, parameter :: thin=1.e-6  !minimum layer thickness
c
      integer i,k,l,lf
      real    q,zk
      real    sic(kk,ks,3),pt(kk+1)
c
      if     (si(1,ks).eq.flag) then
        do k= 1,kz
          do i= 1,ks
            sz(k,i) = flag  !land
          enddo !i
        enddo
      else
c ---   compute PPM coefficients for input layers
        do k=1,kk
          pt(k)=max(p(k+1)-p(k),thin)
        enddo
        call ppm(pt,si,sic,kk,ks)
c
        lf=1
        do k= 1,kz
          zk=z(k)
          do l= lf,kk
            if     (p(l).le.zk .and. p(l+1).ge.zk) then
c
c             z(k) is in layer l, sample the quadratic profile at zk.
c
              q = (zk-p(l))/pt(l)
              do i= 1,ks
                sz(k,i) =    sic(l,i,1) +
     &                    q*(sic(l,i,2) +
     &                       sic(l,i,3)*(1.0-q))
              enddo !i
              lf = l  ! z monotonic increasing, so z(k+1) in layers l:kk
              exit
            elseif (l.eq.kk) then
              do i= 1,ks
                sz(k,i) = flag  ! below the bottom
              enddo !i
              lf = l
              exit
            endif
          enddo !l
        enddo !k
      endif
      return
      end

      subroutine ppm(pt, s,sc,ki,ks)
      implicit none
c
      integer ki,ks
      real    pt(ki+1),s(ki,ks),sc(ki,ks,3)
c
c**********
c*
c  1) generate a monotonic PPM interpolation of a layered field:
c     Colella, P. & P.R. Woodward, 1984, J. Comp. Phys., 54, 174-201.
c
c  2) input arguments:
c       pt    - layer interface thicknesses (non-zero)
c       s     - scalar fields in layer space
c       ki    - 1st dimension of s (number of layers)
c       ks    - 2nd dimension of s (number of fields)
c
c  3) output arguments:
c       sc    - scalar field coefficients for PPM interpolation
c
c  4) except at data voids, must have:
c           pi(   1) == zero (surface)
c           pi( l+1) >= pi(:,:,l)
c           pi(ki+1) == bathymetry
c
c  5) Tim Campbell, Mississippi State University, September 2002;
C     Alan J. Wallcraft,  Naval Research Laboratory,  August 2007.
c*
c**********
c
      integer j,k,l, jm1,jp1
      real    da,a6,slj,scj,srj
      real    as(ki),al(ki),ar(ki)
      real     ptjp(ki), pt2jp(ki), ptj2p(ki),
     &        qptjp(ki),qpt2jp(ki),qptj2p(ki),ptq3(ki),qpt4(ki)
c
      !compute grid metrics
      do j=1,ki
         jm1=max(1,j-1)
         jp1=min(ki,j+1)
         ptq3( j) = pt(j)/(pt(jm1)+pt(j)+pt(jp1))
         ptjp( j) = pt(j)   + pt(jp1)
         pt2jp(j) = pt(j)   + ptjp(j)
         ptj2p(j) = ptjp(j) + pt(jp1)
        qptjp( j) = 1.0/ptjp( j)
        qpt2jp(j) = 1.0/pt2jp(j)
        qptj2p(j) = 1.0/ptj2p(j)
      enddo !j
         ptq3(2) = pt(2)/(pt(1)+ptjp(2))
      do j=3,ki-1
         ptq3(j) = pt(j)/(pt(j-1)+ptjp(j))  !pt(j)/      (pt(j-1)+pt(j)+pt(j+1))
         qpt4(j) = 1.0/(ptjp(j-2)+ptjp(j))  !1.0/(pt(j-2)+pt(j-1)+pt(j)+pt(j+1))
      enddo !j
c
      do l= 1,ks
        !Compute average slopes: Colella, Eq. (1.8)
        as(1)=0.
        do j=2,ki-1
          slj=s(j,  l)-s(j-1,l)
          srj=s(j+1,l)-s(j,  l)
          if (slj*srj.gt.0.) then
            scj=ptq3(j)*( pt2jp(j-1)*srj*qptjp(j)
     &                   +ptj2p(j)  *slj*qptjp(j-1) )
            as(j)=sign(min(abs(2.0*slj),abs(scj),abs(2.0*srj)),scj)
          else
            as(j)=0.
          endif
        enddo !j
        as(ki)=0.
        !Compute "first guess" edge values: Colella, Eq. (1.6)
        al(1)=s(1,l)  !1st layer PCM
        ar(1)=s(1,l)  !1st layer PCM
        al(2)=s(1,l)  !1st layer PCM
        do j=3,ki-1
          al(j)=s(j-1,l)+pt(j-1)*(s(j,l)-s(j-1,l))*qptjp(j-1)
     &         +qpt4(j)*(
     &            2.*pt(j)*pt(j-1)*qptjp(j-1)*(s(j,l)-s(j-1,l))*
     &            ( ptjp(j-2)*qpt2jp(j-1)
     &             -ptjp(j)  *qptj2p(j-1) )
     &            -pt(j-1)*as(j)  *ptjp(j-2)*qpt2jp(j-1)
     &            +pt(j)  *as(j-1)*ptjp(j)  *qptj2p(j-1)
     &              )
          ar(j-1)=al(j)
        enddo !j
        ar(ki-1)=s(ki,l)  !last layer PCM
        al(ki)  =s(ki,l)  !last layer PCM
        ar(ki)  =s(ki,l)  !last layer PCM
        !Impose monotonicity: Colella, Eq. (1.10)
        do j=2,ki-1
          if ((s(j+1,l)-s(j,l))*(s(j,l)-s(j-1,l)).le.0.) then !local extremum
            al(j)=s(j,l)
            ar(j)=s(j,l)
          else
            da=ar(j)-al(j)
            a6=6.0*s(j,l)-3.0*(al(j)+ar(j))
            if     (da*a6 .gt.  da*da) then !peak in right half of zone
              al(j)=3.0*s(j,l)-2.0*ar(j)
            elseif (da*a6 .lt. -da*da) then !peak in left half of zone
              ar(j)=3.0*s(j,l)-2.0*al(j)
            endif
          endif
        enddo !j
        !Set coefficients
        do j=1,ki
          sc(j,l,1)=al(j)
          sc(j,l,2)=ar(j)-al(j)
          sc(j,l,3)=6.0*s(j,l)-3.0*(al(j)+ar(j))
        enddo !j
      enddo !l
      return
      end subroutine ppm

      subroutine layer2z_wno(si,p,kk,ks,
     &                       sz,z,kz,   flag)
      implicit none
c
      integer kk,ks,kz
      real    si(kk,ks),p(kk+1),
     &        sz(kz,ks),z(kz),flag
c
c**********
c*
c  1) interpolate a set of layered field to fixed z depths.
c     method: monotonic WENO-like alternative to PPM across each input cell
c             a second order polynomial approximation of the profiles
c             using a WENO reconciliation of the slopes to compute the
c             interfacial values
c
c  2) input arguments:
c       si    - scalar fields in layer space
c       p     - layer interface depths (non-negative m)
c                 p(   1) is the surface
c                 p(kk+1) is the bathymetry
c       kk    - dimension of si (number of layers)
c       ks    - dimension of si (number of fields)
c       z     - target z-level  depths (non-negative m)
c       flag  - data void (land) marker
c       kz    - dimension of az (number of levels)
c
c  3) output arguments:
c       sz    - scalar fields in z-space
c
c  4) except at data voids, must have:
c           p(   1) == zero (surface)
c           p( l+1) >= p(:,:,l)
c           p(kk+1) == bathymetry
c           0 <= z(k) <= z(k+1)
c     note that z(k) > p(kk+1) implies that az(k)=flag,
c      since the z-level is then below the bathymetry.
c
c  5) Alan J. Wallcraft, Naval Research Laboratory, August 2010.
c*
c**********
c
      real, parameter :: thin=1.e-6  !minimum layer thickness
c
      integer i,k,l,lf
      real    q,q0,q1,q2,zk
      real    c1d(kk,ks,2),dpi(kk)
c
      if     (si(1,ks).eq.flag) then
        do k= 1,kz
          do i= 1,ks
            sz(k,i) = flag  !land
          enddo !i
        enddo
      else
c ---   compute WENO coefficients for input layers
        do k=1,kk
          dpi(k) = max( p(k+1) - p(k), thin )
        enddo
        call weno_coefs(si,dpi,c1d,kk,ks)
c
        lf=1
        do k= 1,kz
          zk=z(k)
          do l= lf,kk
            if     (p(l).le.zk .and. p(l+1).ge.zk) then
c
c             z(k) is in layer l, sample the quadratic profile at zk.
c
              q = (zk-p(l))/dpi(l)
              do i= 1,ks
                q1 = 3.0*q**2 - 4.0*q + 1.0   !1 at q=0, 0 at q=1
                q2 =     q1   + 2.0*q - 1.0   !0 at q=0, 1 at q=1
                q0 = 1.0 - q1 - q2            !0 at q=0, 0 at q=1
                sz(k,i) =q0*si( l,i)   +
     &                   q1*c1d(l,i,1) +
     &                   q2*c1d(l,i,2)
              enddo !i
              lf = l  ! z monotonic increasing, so z(k+1) in layers l:kk
              exit
            elseif (l.eq.kk) then
              do i= 1,ks
                sz(k,i) = flag  ! below the bottom
              enddo !i
              lf = l
              exit
            endif
          enddo !l
        enddo !k
      endif
      return
      end

      subroutine weno_coefs(s,dp,ci,kk,ks)
      implicit none
c
      integer kk,ks
      real    s(kk,ks),dp(kk),ci(kk,ks,2)
c
c-----------------------------------------------------------------------
c  1) coefficents for remaping from one set of vertical cells to another.
c     method: monotonic WENO-like alternative to PPM across each input cell
c             a second order polynomial approximation of the profiles
c             using a WENO reconciliation of the slopes to compute the
c             interfacial values
c
c     REFERENCE: Alexander F. Shchepetkin - Personal Communication
c
c  2) input arguments:
c       s     - initial scalar fields in pi-layer space
c       dp    - initial layer thicknesses (>=thin)
c       kk    - number of layers
c       ks    - number of fields
c
c  3) output arguments:
c       ci    - coefficents for weno_remap
c                ci.1 is value at interface above
c                ci.2 is value at interface below
c
c  4) Laurent Debreu, Grenoble.
c     Alan J. Wallcraft,  Naval Research Laboratory,  July 2008.
c     Method by Alexander F. Shchepetkin.
c-----------------------------------------------------------------------
c
      real, parameter :: dsmll=1.0e-8
c
      integer j,i
      real    q,q01,q02,q001,q002
      real    qdpjm(kk),qdpjmjp(kk),dpjm2jp(kk)
      real    zw(kk+1,3)

      !compute grid metrics
      do j=2,kk-1
        qdpjm(  j) = 1.0/(dp(j-1) +     dp(j))
        qdpjmjp(j) = 1.0/(dp(j-1) +     dp(j) + dp(j+1))
        dpjm2jp(j) =      dp(j-1) + 2.0*dp(j) + dp(j+1)
      enddo !j
      j=kk
        qdpjm(  j) = 1.0/(dp(j-1) +     dp(j))
c
      do i= 1,ks
        do j=2,kk
          zw(j,3) = qdpjm(j)*(s(j,i)-s(j-1,i))
        enddo !j
          j = 1  !PCM first layer
            ci(j,i,1) = s(j,i)
            ci(j,i,2) = s(j,i)
            zw(j,  1) = 0.0
            zw(j,  2) = 0.0
        do j=2,kk-1
          q001 = dp(j)*zw(j+1,3)
          q002 = dp(j)*zw(j,  3)
          if (q001*q002 < 0.0) then
            q001 = 0.0
            q002 = 0.0
          endif
          q01 = dpjm2jp(j)*zw(j+1,3)
          q02 = dpjm2jp(j)*zw(j,  3)
          if     (abs(q001) > abs(q02)) then
            q001 = q02
          endif
          if     (abs(q002) > abs(q01)) then
            q002 = q01
          endif
          q    = (q001-q002)*qdpjmjp(j)
          q001 = q001-q*dp(j+1)
          q002 = q002+q*dp(j-1)

          ci(j,i,2) = s(j,i)+q001
          ci(j,i,1) = s(j,i)-q002
          zw(  j,1) = (2.0*q001-q002)**2
          zw(  j,2) = (2.0*q002-q001)**2
        enddo !j
          j = kk  !PCM last layer
            ci(j,i,1) = s(j,i)
            ci(j,i,2) = s(j,i)
            zw(j,  1) = 0.0
            zw(j,  2) = 0.0

        do j=2,kk
          q002 = max(zw(j-1,2),dsmll)
          q001 = max(zw(j,  1),dsmll)
          zw(j,3) = (q001*ci(j-1,i,2)+q002*ci(j,i,1))/(q001+q002)
        enddo !j
          zw(   1,3) = 2.0*s( 1,i)-zw( 2,3)  !not used?
          zw(kk+1,3) = 2.0*s(kk,i)-zw(kk,3)  !not used?

        do j=2,kk-1
          q01  = zw(j+1,3)-s(j,i)
          q02  = s(j,i)-zw(j,3)
          q001 = 2.0*q01
          q002 = 2.0*q02
          if     (q01*q02 < 0.0) then
            q01 = 0.0
            q02 = 0.0
          elseif (abs(q01) > abs(q002)) then
            q01 = q002
          elseif (abs(q02) > abs(q001)) then
            q02 = q001
          endif
          ci(j,i,1) = s(j,i)-q02
          ci(j,i,2) = s(j,i)+q01
        enddo !j
      enddo !i
      return
      end subroutine weno_coefs


      subroutine layer2z_wnx(si,p,kk,ks,
     &                       sz,z,kz,   flag)
      implicit none
c
      integer kk,ks,kz
      real    si(kk,ks),p(kk+1),
     &        sz(kz,ks),z(kz),flag
c
c**********
c*
c  1) interpolate a set of layered field to fixed z depths.
c     method: WENO-like alternative to PPM across each input cell
c             a second order polynomial approximation of the profiles
c             using a WENO reconciliation of the slopes to compute the
c             interfacial values and allowing new extrema
c
c  2) input arguments:
c       si    - scalar fields in layer space
c       p     - layer interface depths (non-negative m)
c                 p(   1) is the surface
c                 p(kk+1) is the bathymetry
c       kk    - dimension of si (number of layers)
c       ks    - dimension of si (number of fields)
c       z     - target z-level  depths (non-negative m)
c       flag  - data void (land) marker
c       kz    - dimension of az (number of levels)
c
c  3) output arguments:
c       sz    - scalar fields in z-space
c
c  4) except at data voids, must have:
c           p(   1) == zero (surface)
c           p( l+1) >= p(:,:,l)
c           p(kk+1) == bathymetry
c           0 <= z(k) <= z(k+1)
c     note that z(k) > p(kk+1) implies that az(k)=flag,
c      since the z-level is then below the bathymetry.
c
c  5) Alan J. Wallcraft, Naval Research Laboratory, August 2010.
c*
c**********
c
      real, parameter :: thin=1.e-6  !minimum layer thickness
c
      integer i,k,l,lf
      real    q,q0,q1,q2,zk
      real    c1d(kk,ks,2),dpi(kk)
c
      if     (si(1,ks).eq.flag) then
        do k= 1,kz
          do i= 1,ks
            sz(k,i) = flag  !land
          enddo !i
        enddo
      else
c ---   compute WENO coefficients for input layers
        do k=1,kk
          dpi(k) = max( p(k+1) - p(k), thin )
        enddo
        call weno_ext_coefs(si,dpi,c1d,kk,ks)
c
        lf=1
        do k= 1,kz
          zk=z(k)
          do l= lf,kk
            if     (p(l).le.zk .and. p(l+1).ge.zk) then
c
c             z(k) is in layer l, sample the quadratic profile at zk.
c
              q = (zk-p(l))/dpi(l)
              do i= 1,ks
                q1 = 3.0*q**2 - 4.0*q + 1.0   !1 at q=0, 0 at q=1
                q2 =     q1   + 2.0*q - 1.0   !0 at q=0, 1 at q=1
                q0 = 1.0 - q1 - q2            !0 at q=0, 0 at q=1
                sz(k,i) =q0*si( l,i)   +
     &                   q1*c1d(l,i,1) +
     &                   q2*c1d(l,i,2)
              enddo !i
              lf = l  ! z monotonic increasing, so z(k+1) in layers l:kk
              exit
            elseif (l.eq.kk) then
              do i= 1,ks
                sz(k,i) = flag  ! below the bottom
              enddo !i
              lf = l
              exit
            endif
          enddo !l
        enddo !k
      endif
      return
      end

      subroutine weno_ext_coefs(s,dp,ci,kk,ks)
      implicit none
c
      integer kk,ks
      real    s(kk,ks),dp(kk),ci(kk,ks,2)
c
c-----------------------------------------------------------------------
c  1) coefficents for remaping from one set of vertical cells to another.
c     method: WENO-like alternative to PPM across each input cell
c             a second order polynomial approximation of the profiles
c             using a WENO reconciliation of the slopes to compute the
c             interfacial values and allowing new extrema
c
c     REFERENCE: Alexander F. Shchepetkin - Personal Communication
c
c  2) input arguments:
c       s     - initial scalar fields in pi-layer space
c       dp    - initial layer thicknesses (>=thin)
c       kk    - number of layers
c       ks    - number of fields
c
c  3) output arguments:
c       ci    - coefficents for weno_remap
c                ci.1 is value at interface above
c                ci.2 is value at interface below
c
c  4) Laurent Debreu, Grenoble.
c     Alan J. Wallcraft,  Naval Research Laboratory,  July 2008.
c     Method by Alexander F. Shchepetkin.
c-----------------------------------------------------------------------
c
      real, parameter :: dsmll=1.0e-8
c
      integer j,i
      real    q,q01,q02,q001,q002
      real    qdpjm(kk),qdpjmjp(kk),dpjm2jp(kk)
      real    zw(kk+1,3)

      !compute grid metrics
      do j=2,kk-1
        qdpjm(  j) = 1.0/(dp(j-1) +     dp(j))
        qdpjmjp(j) = 1.0/(dp(j-1) +     dp(j) + dp(j+1))
        dpjm2jp(j) =      dp(j-1) + 2.0*dp(j) + dp(j+1)
      enddo !j
      j=kk
        qdpjm(  j) = 1.0/(dp(j-1) +     dp(j))
c
      do i= 1,ks
        do j=2,kk
          zw(j,3) = qdpjm(j)*(s(j,i)-s(j-1,i))
        enddo !j
          j = 1  !PCM first layer
            ci(j,i,1) = s(j,i)
            ci(j,i,2) = s(j,i)
            zw(j,  1) = 0.0
            zw(j,  2) = 0.0
        do j=2,kk-1
          q001 = dp(j)*zw(j+1,3)
          q002 = dp(j)*zw(j,  3)
          if (q001*q002 < 0.0) then
            q001 = 0.0
            q002 = 0.0
          endif
          q01 = dpjm2jp(j)*zw(j+1,3)
          q02 = dpjm2jp(j)*zw(j,  3)
          if     (abs(q001) > abs(q02)) then
            q001 = q02
          endif
          if     (abs(q002) > abs(q01)) then
            q002 = q01
          endif
          q    = (q001-q002)*qdpjmjp(j)
          q001 = q001-q*dp(j+1)
          q002 = q002+q*dp(j-1)

          ci(j,i,2) = s(j,i)+q001
          ci(j,i,1) = s(j,i)-q002
          zw(  j,1) = (2.0*q001-q002)**2
          zw(  j,2) = (2.0*q002-q001)**2
        enddo !j
          j = kk  !PCM last layer
            ci(j,i,1) = s(j,i)
            ci(j,i,2) = s(j,i)
            zw(j,  1) = 0.0
            zw(j,  2) = 0.0

*       do j=2,kk
*         q002 = max(zw(j-1,2),dsmll)
*         q001 = max(zw(j,  1),dsmll)
*         zw(j,3) = (q001*ci(j-1,i,2)+q002*ci(j,i,1))/(q001+q002)
*       enddo !j
*         zw(   1,3) = 2.0*s( 1,i)-zw( 2,3)  !not used?
*         zw(kk+1,3) = 2.0*s(kk,i)-zw(kk,3)  !not used?

*       do j=2,kk-1
*         q01  = zw(j+1,3)-s(j,i)
*         q02  = s(j,i)-zw(j,3)
*         q001 = 2.0*q01
*         q002 = 2.0*q02
*         if     (q01*q02 < 0.0) then
*           q01 = 0.0
*           q02 = 0.0
*         elseif (abs(q01) > abs(q002)) then
*           q01 = q002
*         elseif (abs(q02) > abs(q001)) then
*           q02 = q001
*         endif
*         ci(j,i,1) = s(j,i)-q02
*         ci(j,i,2) = s(j,i)+q01
*       enddo !j
        ci(2,i,1) = ci(1,i,2)
        do j=3,kk-1
          q01 = 0.5*(ci(j-1,i,2)+ci(j,i,1))
          ci(j-1,i,2) = q01
          ci(j,  i,1) = q01
        enddo !j
        ci(kk-1,i,2) = ci(kk,i,1)
      enddo !i
      return
      end subroutine weno_ext_coefs

      subroutine pt2t(temp, ptemp,saln,z,plat)
      implicit none
c
      real temp, ptemp,saln,z,plat
c
c --- calculate temperature from potential temperature
c
      external p80,theta
      real     p80,theta
      real     dbar
c
      dbar  = p80(z,plat)
      temp  = theta(saln,ptemp, 0.0,dbar)
      end

c separate set of subroutines, adapted from WHOI CTD group
c
      real function atg(s,t,p)
c ****************************
c adiabatic temperature gradient deg c per decibar
c ref: bryden,h.,1973,deep-sea res.,20,401-408
c units:
c       pressure        p        decibars
c       temperature     t        deg celsius (ipts-68)
c       salinity        s        (ipss-78)
c       adiabatic       atg      deg. c/decibar
c checkvalue: atg=3.255976e-4 c/dbar for s=40 (ipss-78),
c t=40 deg c,p0=10000 decibars
      ds = s - 35.0
      atg = (((-2.1687e-16*t+1.8676e-14)*t-4.6206e-13)*p
     x+((2.7759e-12*t-1.1351e-10)*ds+((-5.4481e-14*t
     x+8.733e-12)*t-6.7795e-10)*t+1.8741e-8))*p
     x+(-4.2393e-8*t+1.8932e-6)*ds
     x+((6.6228e-10*t-6.836e-8)*t+8.5258e-6)*t+3.5803e-5
      return
      end
      real function theta(s,t0,p0,pr)
c ***********************************
c to compute local potential temperature at pr
c using bryden 1973 polynomial for adiabatic lapse rate
c and runge-kutta 4-th order integration algorithm.
c ref: bryden,h.,1973,deep-sea res.,20,401-408
c fofonoff,n.,1977,deep-sea res.,24,489-491
c units:
c       pressure        p0       decibars
c       temperature     t0       deg celsius (ipts-68)
c       salinity        s        (ipss-78)
c       reference prs   pr       decibars
c       potential tmp.  theta    deg celsius
c checkvalue: theta= 36.89073 c,s=40 (ipss-78),t0=40 deg c,
c p0=10000 decibars,pr=0 decibars
c
c      set-up intermediate temperature and pressure variables
      p=p0
      t=t0
c**************
      h = pr - p
      xk = h*atg(s,t,p)
      t = t + 0.5*xk
      q = xk
      p = p + 0.5*h
      xk = h*atg(s,t,p)
      t = t + 0.29289322*(xk-q)
      q = 0.58578644*xk + 0.121320344*q
      xk = h*atg(s,t,p)
      t = t + 1.707106781*(xk-q)
      q = 3.414213562*xk - 4.121320344*q
      p = p + 0.5*h
      xk = h*atg(s,t,p)
      theta = t + (xk-2.0*q)/6.0
      return
      end
c
c pressure from depth from saunder's formula with eos80.
c reference: saunders,peter m., practical conversion of pressure
c            to depth., j.p.o. , april 1981.
c r millard
c march 9, 1983
c check value: p80=7500.004 dbars;for lat=30 deg., depth=7321.45 meters
      real function p80(dpth,xlat)
      real, parameter :: pi=3.141592654
      plat=abs(xlat*pi/180.)
      d=sin(plat)
      c1=5.92e-3+d**2*5.25e-3
      p80=((1-c1)-sqrt(((1-c1)**2)-(8.84e-6*dpth)))/4.42e-6
      return
      end
