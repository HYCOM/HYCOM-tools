      PROGRAM PROFILEI
      IMPLICIT NONE
C
C  hycom_profile2zi - Usage: hycom_profile2zi archv.txt zi.txt archz.txt [itype]
C
C                 converts an HYCOM isopycnal text profile file to
C                 a finite volume z-level text profile file.
C
C   archv.txt is assumed to be an HYCOM archive text profile file
C   zi.txt    is assumed to contain a list of interface depths, one per line
C               the first layer is assumed to be from the surface.
C               the last  layer need not reach the bottom, so set the
C               last interface depth very deep to guarentee full coverage.
C               z-levels are nominally half way between interface depths.
C   archz.txt will be the output, new cell, text profile file
C   itype     is the input interpolation type (default 1)
C                =0; piecewise constant  method (PCM) or donor cell
C                =1; piecewise linear    method (PLM) or VanLeer
C                =2; piecewise parabolic method (PPM)
C
C  all input interpolation schemes conserve the original cell average,
C  and the output is the average of the interpolation profile across
C  each new cell.  So the total depth average is conserved (provided
C  the last interface is at or below the total depth).
C
C  If some or all of the target coordinate system is in density space,
C  use hycom_profile_remap to both choose the interface depths and
C  perform the averaging over the new cells.
C
C  this version for "serial" Unix systems.
C
C  Tim Campbell,       Mississippi State University
C  Alan J. Wallcraft,  Naval Research Laboratory,  Sep. 2002 (Aug. 2005).
C
      INTEGER      IARGC
      INTEGER      NARG
      CHARACTER*240 CARG
C
      CHARACTER*240 CFILEA,CFILEB,CFILEC
      CHARACTER*240 CLINE
      REAL          THK,DUM5(5),FLAG
      REAL          PO(9999)
      REAL          SIT(5),SOT(5)
      INTEGER       IOS,ITYPE,K,KK,KDM,KO
      INTEGER       I,KMAX
C
      REAL, ALLOCATABLE :: SI(:,:),PI(:),SO(:,:),ZO(:)
C
      INTEGER       NSAMP
      REAL          PMAX
      PARAMETER(NSAMP=1,PMAX=500.0)
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
        ITYPE = 1
      ELSE
        WRITE(6,*)
     +    'Usage: hycom_profile2zi archv.txt zi.txt archz.txt [itype]'
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
C
      DO K= 1,4
        READ( 11,'(a)')      CLINE
        WRITE(21,'(a)') TRIM(CLINE)
      ENDDO
      READ( 11,'(a)') CLINE
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
      ALLOCATE( PI(KDM+1), SI(KDM,5) )
C
      REWIND(11)
      DO K= 1,5
        READ(11,*)
      ENDDO
      PI(1) =  0.0
      DO K= 1,KDM
        READ(11,'(a)',IOSTAT=IOS) CLINE
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'Error: inconsistent input profile'
          CALL EXIT(6)
        ENDIF
        READ(CLINE,*) KDM,SI(K,1),SI(K,2),SI(K,3),SI(K,4),SI(K,5),THK
        PI(K+1) = PI(K) + THK
        IF     (THK.EQ.0.0) THEN
          DO KK= 1,5
            SI(K,KK)=SI(K-1,KK)
          ENDDO !kk
        ENDIF
      ENDDO
      CLOSE(11)
C
C     STORE COPY OF INPUT AVERAGES AND
C     COMPUTE INPUT TOTAL COLUMN AVERAGES
C
      KMAX=0
      DO KK= 1,5
        SIT(KK)=0.0
      ENDDO !kk
      DO K=1,KDM
        IF (PI(K).LE.PMAX) KMAX=KMAX+1
        THK=PI(K+1)-PI(K)
        DO KK= 1,5
          SIT(KK)=SIT(KK)+THK*SI(K,KK)
        ENDDO !kk
      ENDDO
      DO KK= 1,5
        SIT(KK)=SIT(KK)/PI(KDM+1)
      ENDDO !kk
C
C     READ THE ZI FILE.
C
      PO(1) = 0.0
      DO K= 2,99999
        READ(12,'(a)',IOSTAT=IOS) CLINE
        IF     (IOS.NE.0) THEN
          KO = K-2
          EXIT
        ENDIF
        READ(CLINE,*) PO(K)
        IF     (K.EQ.2 .AND. PO(K).LE.0.0) THEN
          WRITE(6,*) 'Error: inconsistent output profile'
          WRITE(6,*) '1st value must be > 0.0'
          WRITE(6,*) PO(2)
          CALL EXIT(7)
        ELSEIF (PO(K).LT.PO(K-1)) THEN
          WRITE(6,*) 'Error: inconsistent output profile'
          WRITE(6,*) k-1,'-th value < previous value'
          WRITE(6,*) PO(2:K)
          CALL EXIT(8)
        ENDIF
        PO(K) =MIN( PO(K), PI(KDM+1) )
      ENDDO
      CLOSE(12)
C
      ALLOCATE( ZO(KO), SO(KO,5) )
C
      DO K= 1,KO
        ZO(K) = 0.5*(PO(K)+PO(K+1))
      ENDDO
C
C     FORM NEW CELL AVERAGES
C
      FLAG = 9999.99
      IF     (ITYPE.EQ.0) THEN
        CALL LAYER2ZI_PCM(SI,PI,KDM,5,
     &                    SO,PO,KO,   FLAG)
      ELSEIF (ITYPE.EQ.1) THEN
        CALL LAYER2ZI_PLM(SI,PI,KDM,5,
     &                    SO,PO,KO,   FLAG)
      ELSEIF (ITYPE.EQ.2) THEN
        CALL LAYER2ZI_PPM(SI,PI,KDM,5,
     &                    SO,PO,KO,   FLAG)
      ELSE
        WRITE(6,*)
     +    'Usage: hycom_profile2zi archv.txt zi.txt archz.txt [itype]'
        WRITE(6,*) 'unsupported itype value'
        CALL EXIT(1)
      ENDIF
C
C     COMPUTE AND PRINT NEW TOTAL COLUMN AVERAGE AND DEVIATION METRIC
C
      DO KK= 1,5
        SOT(KK)=0.0
        DO K=1,KO
          THK=PO(K+1)-PO(K)
          SOT(KK)=SOT(KK)+THK*SO(K,KK)
        ENDDO !k
        SOT(KK)=SOT(KK)/PO(KO+1)
      ENDDO !kk
      WRITE(*,'(/A)')       '        LABEL:  U  V  T  S  R'
      WRITE(*,'(A,5e14.6)') ' INPUT TOTALS: ',SIT(:)
      WRITE(*,'(A,5e14.6)') 'OUTPUT TOTALS: ',SOT(:)
C
C     OUTPUT.
C
      WRITE(21,'(a,a)')
     &  '#   k',
     &  '    utot    vtot    temp    saln    dens    thkns      dpth'
      DO K= 1,KO
        WRITE(21,'(i4,1x,2f8.2,3f8.3,f9.3,f10.3)')
     &    K,
     &    SO(K,1),             !cm/s
     &    SO(K,2),             !cm/s
     &    SO(K,3),             !degC
     &    SO(K,4),             !psu
     &    SO(K,5),             !SigmaT
     &    PO(K+1)-PO(K),       !m
     &    ZO(K)                !m
      ENDDO
      END

      subroutine layer2zi_pcm(si,pi,ki,ks,
     &                        so,po,ko,   flag)
      implicit none
c
      integer ki,ks,ko
      real    si(ki,ks),pi(ki+1),
     &        so(ko,ks),po(ko+1),flag
c
c**********
c*
c  1) remap from one set of vertical cells to another.
c     method: piecewise constant across each input cell
c             the output is the average of the interpolation
c             profile across each output cell.
c
c  2) input arguments:
c       si    - scalar fields in pi-layer space
c       pi    - layer interface depths (non-negative m)
c                 pi(   1) is the surface
c                 pi(ki+1) is the bathymetry
c       ki    - 1st dimension of si     (number of  input layers)
c       ks    - 2nd dimension of si,so  (number of fields)
c       po    - target interface depths (non-negative m)
c                 po(k+1) >= po(k)
c       ko    - 1st dimension of so     (number of output layers)
c       flag  - data void (land) marker
c
c  3) output arguments:
c       so    - scalar fields in po-layer space
c
c  4) except at data voids, must have:
c           pi(   1) == zero (surface)
c           pi( l+1) >= pi(l)
c           pi(ki+1) == bathymetry
c           0 <= po(k) <= po(k+1)
c      output layers completely below the bathymetry inherit values
c      from the layer above.
c
c  5) Alan J. Wallcraft,  Naval Research Laboratory,  Sep. 2002 (Aug. 2005).
c*
c**********
c
      real       thin 
      parameter (thin=1.e-6)  ! minimum layer thickness (no division by 0.0)
c
      integer i,k,l,lf
      real    q,zb,zt,sok(ks)
c
      if     (si(1,1).eq.flag) then
        do k= 1,ko
          do i= 1,ks
            so(k,i) = flag  !land
          enddo !i
        enddo !k
      else
        lf=1
        zb=po(1)
        do k= 1,ko
          zt = zb
          zb = po(k+1)
*         WRITE(6,*) 'k,zt,zb = ',k,zt,zb
          if     (zb-zt.lt.thin .or. zt.ge.pi(ki+1)) then
c
c ---       thin or bottomed layer, values taken from layer above
c
            do i= 1,ks
              so(k,i) = so(k-1,i)
            enddo !i
          else
c
c           form layer averages.
c
            if     (pi(lf).gt.zt) then
              WRITE(6,*) 'bad lf = ',lf
              stop
            endif
            do i= 1,ks
              sok(i) = 0.0
            enddo !i
            do l= lf,ki
              if     (pi(l).gt.zb) then
*               WRITE(6,*) 'l,lf= ',l,lf,l-1
                lf = l-1
                exit
              elseif (pi(l).ge.zt .and. pi(l+1).le.zb) then
c
c               the input layer is completely inside the output layer
c
                q   = max(pi(l+1)-pi(l),thin)/(zb-zt)
                do i= 1,ks
                  sok(i) = sok(i) + q*si(l,i)
                enddo !i
*               WRITE(6,*) 'L,q = ',l,q
              else
c
c               the input layer is partially inside the output layer
c
                q   = max(min(pi(l+1),zb)-max(pi(l),zt),thin)/(zb-zt)
                do i= 1,ks
                  sok(i) = sok(i) + q*si(l,i)
                enddo !i
*               WRITE(6,*) 'l,q = ',l,q
              endif
            enddo !l
            do i= 1,ks
              so(k,i) = sok(i)
            enddo !i
          endif
        enddo !k
      endif
      return
      end subroutine layer2zi_pcm

      subroutine layer2zi_plm(si,pi,ki,ks,
     &                        so,po,ko,   flag)
      implicit none
c
      integer ki,ks,ko
      real    si(ki,ks),pi(ki+1),
     &        so(ko,ks),po(ko+1),flag
c
c**********
c*
c  1) remap from one set of vertical cells to another.
c     method: piecewise linear across each input cell
c             the output is the average of the interpolation
c             profile across each output cell.
c
c  2) input arguments:
c       si    - scalar fields in pi-layer space
c       pi    - layer interface depths (non-negative m)
c                 pi(   1) is the surface
c                 pi(ki+1) is the bathymetry
c       ki    - 1st dimension of si     (number of  input layers)
c       ks    - 2nd dimension of si,so  (number of fields)
c       po    - target interface depths (non-negative m)
c                 po(k+1) >= po(k)
c       ko    - 1st dimension of so     (number of output layers)
c       flag  - data void (land) marker
c
c  3) output arguments:
c       so    - scalar fields in po-layer space
c
c  4) except at data voids, must have:
c           pi(   1) == zero (surface)
c           pi( l+1) >= pi(l)
c           pi(ki+1) == bathymetry
c           0 <= po(k) <= po(k+1)
c      output layers completely below the bathymetry inherit values
c      from the layer above.
c
c  5) Tim Campbell, Mississippi State University, October 2002.
C     Alan J. Wallcraft,  Naval Research Laboratory,  Aug. 2005.
c*
c**********
c
      real,parameter :: thin=1.e-6  !minimum layer thickness
c
      integer i,k,l,lf
      real    q,qc,zb,zc,zt,sok(ks)
      real    sis(ki,ks),pit(ki+1)
c
      if     (si(1,1).eq.flag) then
        do k= 1,ko
          do i= 1,ks
            so(k,i) = flag  !land
          enddo !i
        enddo
      else
c ---   compute PLM slopes for input layers
        do k=1,ki
          pit(k)=max(pi(k+1)-pi(k),thin)
        enddo
        call plm(pit,si,sis,ki,ks)
c ---   compute output layer averages
        lf=1
        zb=po(1)
        do k= 1,ko
          zt = zb
          zb = po(k+1)
*         WRITE(6,*) 'k,zt,zb = ',k,zt,zb
          if     (zb-zt.lt.thin .or. zt.ge.pi(ki+1)) then
c
c ---       thin or bottomed layer, values taken from layer above
c
            do i= 1,ks
              so(k,i) = so(k-1,i)
            enddo !i
          else
c
c           form layer averages.
c
            if     (pi(lf).gt.zt) then
              WRITE(6,*) 'bad lf = ',lf
              stop
            endif
            do i= 1,ks
              sok(i) = 0.0
            enddo !i
            do l= lf,ki
              if     (pi(l).gt.zb) then
*               WRITE(6,*) 'l,lf= ',l,lf,l-1
                lf = l-1
                exit
              elseif (pi(l).ge.zt .and. pi(l+1).le.zb) then
c
c               the input layer is completely inside the output layer
c
                q   = max(pi(l+1)-pi(l),thin)/(zb-zt)
                do i= 1,ks
                  sok(i) = sok(i) + q*si(l,i)
                enddo !i
*               WRITE(6,*) 'L,q = ',l,q
              else
c
c               the input layer is partially inside the output layer
c               average of linear profile is its center value
c
                q   = max( min(pi(l+1),zb)-max(pi(l),zt), thin )/(zb-zt)
                zc  = 0.5*(min(pi(l+1),zb)+max(pi(l),zt))
                qc  = (zc-pi(l))/pit(l) - 0.5
                do i= 1,ks
                  sok(i) = sok(i) + q*(si(l,i) + qc*sis(l,i))
                enddo !i
*               WRITE(6,*) 'l,q,qc = ',l,q,qc
              endif
            enddo !l
            do i= 1,ks
              so(k,i) = sok(i)
            enddo !i
          endif
        enddo !k
      endif
      return
      end subroutine layer2zi_plm

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

      subroutine layer2zi_ppm(si,pi,ki,ks,
     &                        so,po,ko,   flag)
      implicit none
c
      integer ki,ks,ko
      real    si(ki,ks),pi(ki+1),
     &        so(ko,ks),po(ko+1),flag
c
c**********
c*
c  1) remap from one set of vertical cells to another.
c     method: piecewise parabolic method across each input cell
c             the output is the average of the interpolation
c             profile across each output cell.
c
c  2) input arguments:
c       si    - scalar fields in pi-layer space
c       pi    - layer interface depths (non-negative m)
c                 pi(   1) is the surface
c                 pi(ki+1) is the bathymetry
c       ki    - 1st dimension of si     (number of  input layers)
c       ks    - 2nd dimension of si,so  (number of fields)
c       po    - target interface depths (non-negative m)
c                 po(k+1) >= po(k)
c       ko    - 1st dimension of so     (number of output layers)
c       flag  - data void (land) marker
c
c  3) output arguments:
c       so    - scalar fields in po-layer space
c
c  4) except at data voids, must have:
c           pi(   1) == zero (surface)
c           pi( l+1) >= pi(l)
c           pi(ki+1) == bathymetry
c           0 <= po(k) <= po(k+1)
c      output layers completely below the bathymetry inherit values
c      from the layer above.
c
c  5) Tim Campbell, Mississippi State University, October 2002.
C     Alan J. Wallcraft,  Naval Research Laboratory,  Aug. 2005.
c*
c**********
c
      real,parameter :: thin=1.e-6  !minimum layer thickness
c
      integer i,k,l,lb,lt
      real    sz,xb,xt,zb,zt
      real    sic(ki,ks,3),pit(ki+1)
c
      if     (si(1,1).eq.flag) then
        do k= 1,ko
          do i= 1,ks
            so(k,i) = flag  !land
          enddo !i
        enddo
      else
c ---   compute PPM coefficients for input layers
        do k=1,ki
          pit(k)=max(pi(k+1)-pi(k),thin)
        enddo
        call ppm(pit,si,sic,ki,ks)
c ---   compute output layer averages
        lb=1
        zb=po(1)
        do k= 1,ko
          zt = zb
          zb = po(k+1)
*         WRITE(6,*) 'k,zt,zb = ',k,zt,zb
          if     (zb-zt.lt.thin .or. zt.ge.pi(ki+1)) then
c
c ---       thin or bottomed layer, values taken from layer above
c
            do i= 1,ks
              so(k,i) = so(k-1,i)
            enddo !i
          else
c
c           form layer averages.
c
            if     (pi(lb).gt.zt) then
              WRITE(6,*) 'bad lb = ',lb
              stop
            endif
            lt=lb !top will always correspond to bottom of previous
            lb=lt !find input layer containing bottom output interface
            do while (pi(lb+1).lt.zb.and.lb.lt.ki)
              lb=lb+1
            enddo
            xt=(zt-pi(lt))/pit(lt)
            xb=(zb-pi(lb))/pit(lb)
            do i= 1,ks
              if     (lt.ne.lb) then
                sz=   pit(lt)*(      sic(lt,i,1)  *(1.-xt)
     &                         +0.5*(sic(lt,i,2)+
     &                               sic(lt,i,3) )*(1.-xt**2)
     &                              -sic(lt,i,3)  *(1.-xt**3)/3.0)
                do l=lt+1,lb-1
                  sz=sz+pit(l)*si(l,i)
                enddo !l
                sz=sz+pit(lb)*(      sic(lb,i,1)  *    xb
     &                         +0.5*(sic(lb,i,2)+
     &                               sic(lb,i,3) )*    xb**2
     &                              -sic(lb,i,3)  *    xb**3 /3.0)
              else
                sz=pit(lt)*(      sic(lt,i,1)  *(xb-xt)
     &                      +0.5*(sic(lt,i,2)+
     &                            sic(lt,i,3) )*(xb**2-xt**2)
     &                           -sic(lt,i,3)  *(xb**3-xt**3)/3.0)
              endif
              so(k,i) = sz/(zb-zt)
            enddo !i
          endif !thin:std layer
        enddo !k
      endif
      return
      end subroutine layer2zi_ppm

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
      integer j,k,l
      real    da,a6,slj,scj,srj
      real    as(ki),al(ki),ar(ki)
      real     ptjp(ki), pt2jp(ki), ptj2p(ki),
     &        qptjp(ki),qpt2jp(ki),qptj2p(ki),ptq3(ki),qpt4(ki)
c
      !compute grid metrics
      do j=1,ki
         ptq3( j) = pt(j)/(pt(j-1)+pt(j)+pt(j+1))
         ptjp( j) = pt(j)   + pt(j+1)
         pt2jp(j) = pt(j)   + ptjp(j)
         ptj2p(j) = ptjp(j) + pt(j+1)
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
