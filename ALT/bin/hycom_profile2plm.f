      PROGRAM HYCOM_PROFILE2PLM
      IMPLICIT NONE
C
C  hycom_profile2plm - Usage: hycom_profile2plm archv.txt archp.txt [hybgen]
C
C                 converts a HYCOM text profile file to one
C                 which will plot as the PLM profile
C
C   archv.txt is assumed to be an HYCOM archive text profile file
C   archp.txt will be the output, PLM plot, text profile file
C   if hybgen is present, the hybgen PLM will be used
C
C  the input and the output are equivalent HYCOM profiles, but
C  the ouput adds two zero thickness layers at each interface
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  March 2002
C
      REAL          THIN 
      PARAMETER    (THIN=1.E-6)  !minimum layer thickness
C
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      CHARACTER*240 CFILEA,CFILEC
      CHARACTER*240 CLINE
      REAL          THK,FLAG
      LOGICAL       HYBGEN
      INTEGER       IOS,K,KDM,KK
C
      REAL, ALLOCATABLE :: SI(:,:),SS(:,:),PI(:),DP(:)
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.3) THEN
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CFILEC)
        HYBGEN = .true.
      ELSEIF (NARG.EQ.2) THEN
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CFILEC)
        HYBGEN = .false.
      ELSE
        WRITE(6,*)
     +    'Usage: hycom_profile2plm archv.txt archz.txt [hybgen]'
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
      OPEN(UNIT=21, FILE=CFILEC, FORM='FORMATTED', STATUS='NEW',
     +     IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',TRIM(CFILEC)
        WRITE(6,*) 'ios   = ',ios
        CALL EXIT(5)
      ENDIF
C
C     COPY PROFILE HEADER TO OUTPUT.
C
      DO K= 1,99
        READ( 11,'(a)')      CLINE
        IF     (CLINE(1:5).EQ.'#  k ') then
          EXIT
        ENDIF
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
C     RE-READ THE ISOPYCNAL PROFILE.
C
      ALLOCATE( PI(KDM+1), DP(KDM), SI(KDM,5), SS(KDM,5) )
C
      REWIND(11)
      DO K= 1,99
        READ( 11,'(a)') CLINE
        IF     (CLINE(1:5).EQ.'#  k ') then
          EXIT
        ENDIF
      ENDDO
      PI(1) =  0.0
      DO K= 1,KDM
        READ(11,'(a)',IOSTAT=IOS) CLINE
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'Error: inconsistent input profile'
          CALL EXIT(6)
        ENDIF
        READ(CLINE,*) KDM,(SI(K,KK),KK=1,5),THK
        DP(K)   = MAX( THK, THIN )
        PI(K+1) = PI(K) + THK
        IF     (THK.EQ.THIN) THEN
          DO KK= 1,5
            SI(K,KK)=SI(K-1,KK)
          ENDDO !kk
        ENDIF
      ENDDO
      CLOSE(11)
C
C     CALCULATE PLM SLOPES.
C
      IF     (HYBGEN) THEN
        CALL PLM_HYBGEN(DP, SI,SS,KDM,5, THIN)
      ELSE
        CALL PLM_INTERP(DP, SI,SS,KDM,5)
      ENDIF
C
C     OUTPUT, 3 COPIES OF EACH LAYER.
C
      WRITE(21,'(3a)')
     &      '#  k',
     &      '    utot    vtot  p.temp    saln  p.dens',
     &      '    thkns      dpth'
      DO K= 1,KDM
        THK = PI(K+1) - PI(K)
        WRITE(21,'(i4,2f8.2,3f8.3,f9.3,f10.3)')
     &    3*K-2,(SI(K,KK)-0.5*SS(K,KK),KK=1,5),0.0,     PI(K)
        WRITE(21,'(i4,2f8.2,3f8.3,f9.3,f10.3)')
     &    3*K-1,(SI(K,KK),             KK=1,5),THK,0.5*(PI(K)+PI(K+1))
        WRITE(21,'(i4,2f8.2,3f8.3,f9.3,f10.3)')
     &    3*K  ,(SI(K,KK)+0.5*SS(K,KK),KK=1,5),0.0,           PI(K+1)
      ENDDO
      END

      subroutine plm_interp(pt, s,ss,ki,ks)
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
      end subroutine plm_interp

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

      subroutine plm_hybgen(dpi, si,ci,kk,ks, thin)
      implicit none
c
      integer kk,ks
      real    si(kk,ks),dpi(kk),ci(kk,ks),thin
c
c-----------------------------------------------------------------------
c  1) coefficents for remaping from one set of vertical cells to another.
c     method: piecewise linear across each input cell
c             the output is the average of the interpolation
c             profile across each output cell.
c
c  2) input arguments:
c       si    - initial scalar fields in pi-layer space
c       dpi   - initial layer thicknesses (>= thin)
c       kk    - number of layers
c       ks    - number of fields
c       thin  - layer thickness (>0) that can be ignored
c
c  3) output arguments:
c       ci    - coefficents (slopes) for hybgen_plm_remap
c
c  4) Tim Campbell, Mississippi State University, October 2002.
c     Alan J. Wallcraft,  Naval Research Laboratory,  Aug. 2007.
c-----------------------------------------------------------------------
c
      integer k,i
      real    qcen,zbot,zcen,ztop
c
      do i= 1,ks
        ci(1, i) = 0.0
        ci(kk,i) = 0.0
      enddo !i
      do k= 2,kk-1
        if     (dpi(k).le.thin) then  !use PCM
          do i= 1,ks
            ci(k,i) = 0.0
          enddo !i
        else
c ---     use qcen in place of 0.5 to allow for non-uniform grid
          qcen = dpi(k)/(dpi(k)+0.5*(dpi(k-1)+dpi(k+1)))
          do i= 1,ks
c ---       PLM (non-zero slope, but no new extrema)
c ---       layer value is si-0.5*ci at top    interface,
c ---                  and si+0.5*ci at bottom interface.
c
c ---       monotonized central-difference limiter (van Leer, 1977,
c ---       JCP 23 pp 276-299).  For a discussion of PLM limiters, see
c ---       Finite Volume Methods for Hyperbolic Problems by R.J. Leveque.
            ztop = 2.0*(si(k,  i)-si(k-1,i))
            zbot = 2.0*(si(k+1,i)-si(k,  i))
            zcen =qcen*(si(k+1,i)-si(k-1,i))
            if     (ztop*zbot.gt.0.0) then !ztop,zbot are the same sign
              ci(k,i)=sign(min(abs(zcen),abs(zbot),abs(ztop)),zbot)
            else
              ci(k,i)=0.0  !local extrema, so no slope
            endif
          enddo !i
        endif  !PCM:PLM
      enddo !k
      return
      end subroutine plm_hybgen
