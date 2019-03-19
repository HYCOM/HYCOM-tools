      PROGRAM FNEWZI
      IMPLICIT NONE
C
C  hycom_newzi - Usage:  hycom_newzi fin.a idm jdm kz zi.txt pnew.a kdm fout.a [itype]
C
C                 Outputs kdm (1:idm,1:jdm) fields, representing the
C                 vertical remapping from zi(1:kz) to pnew(:,:,1:kdm).
C
C  fin.a and pnew.a are assumed to contain idm*jdm 32-bit IEEE real values
C   for each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C   fin.a     is the layer average across layers defined by zi(1:kz) 
C             if itype=4xx the fields are sampled at "cell centers"
C   zi.txt    is a list of kz interface depths, one per line, 
C               the first layer is assumed to be from the surface and 
C               the last  layer need not reach the bottom, so set the
C               last interface depth very deep to guarentee full coverage.
C             if itype=4xx the input values are "cell centers"
C   pnew.a    is the target kdm interface depth arrays
C   fout.a    will be the average across layers defined by pnew(:,:,1:kdm)
C             if itype=xx3 it will be sampled at pnew's cell centers
C   itype     is the 3-digit input interpolation type (default 001)
C                =xx0; piecewise constant  method (PCM) or donor cell
C                =xx1; piecewise linear    method (PLM) or VanLeer
C                =xx2; linear between cell centers
C                =xx3; linear between cell centers, at cell center
C                =x0x; pnew is interface depth in m
C                =x1x; pnew is interface depth in pressure units
C                =x2x; pnew is layer thickness in m
C                =x3x; pnew is layer thickness in pressure units
C                =0xx; zi   is interface depth in m
C                =4xx; zi   is cell center     in m
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  March 2005.
C
      REAL*4, ALLOCATABLE :: AZ(:,:,:),P(:,:,:),AP(:,:,:)
      REAL*4              :: PAD(4096)
      INTEGER       IOS,L
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      INTEGER      IDM,JDM,KZ,KDM,NPAD,ITYPE,ITEST,JTEST
      CHARACTER*240 CFILE1,CFILEZ,CFILEP,CFILEO
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.8) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CARG)
        READ(CARG,*) KZ
        CALL GETARG(5,CFILEZ)
        CALL GETARG(6,CFILEP)
        CALL GETARG(7,CARG)
        READ(CARG,*) KDM
        CALL GETARG(8,CFILEO)
        ITYPE = 001
        ITEST = 0
        JTEST = 0
      ELSEIF (NARG.EQ.9) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CARG)
        READ(CARG,*) KZ
        CALL GETARG(5,CFILEZ)
        CALL GETARG(6,CFILEP)
        CALL GETARG(7,CARG)
        READ(CARG,*) KDM
        CALL GETARG(8,CFILEO)
        CALL GETARG(9,CARG)
        READ(CARG,*) ITYPE
        ITEST = 0
        JTEST = 0
      ELSEIF (NARG.EQ.11) THEN  !undocumented debug option
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CARG)
        READ(CARG,*) KZ
        CALL GETARG(5,CFILEZ)
        CALL GETARG(6,CFILEP)
        CALL GETARG(7,CARG)
        READ(CARG,*) KDM
        CALL GETARG(8,CFILEO)
        CALL GETARG(9,CARG)
        READ(CARG,*) ITYPE
        CALL GETARG(10,CARG)
        READ(CARG,*) ITEST
        CALL GETARG(11,CARG)
        READ(CARG,*) JTEST
      ELSE
        WRITE(6,'(3a)')
     &    'Usage:  ',
     &    'hycom_newzi ',
     &    'fin.a idm jdm kz zi.txt pnew.a kdm fout.a [itype]'
        CALL EXIT(1)
      ENDIF
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( AZ(IDM,JDM,KZ), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_newzi: could not allocate ',
     +             IDM*JDM*KZ,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( P(IDM,JDM,KDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_newzi: could not allocate ',
     +             IDM*JDM*KDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( AP(IDM,JDM,KDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_newzi: could not allocate ',
     +             IDM*JDM*KDM,' words'
        CALL EXIT(2)
      ENDIF
C
      CALL NEWZI(AZ,P,AP,IDM,JDM,KZ,KDM,PAD,NPAD,
     &          ITYPE, ITEST,JTEST, CFILE1,CFILEZ,CFILEP,CFILEO)
      CALL EXIT(0)
      END
      SUBROUTINE NEWZI(AZ,P,AP,IDM,JDM,KZ,KDM,PAD,NPAD,
     &                ITYPE, ITEST,JTEST, CFILE1,CFILEZ,CFILEP,CFILEO)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE1,CFILEZ,CFILEP,CFILEO
      INTEGER      IDM,JDM,KZ,KDM,NPAD,ITYPE,ITEST,JTEST
      REAL*4       AZ(IDM,JDM,KZ),
     +             P(IDM,JDM,KDM),AP(IDM,JDM,KDM),PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      CHARACTER*18 CASN
      INTEGER      ITYPE_I,ITYPE_P,ITYPE_Z
      INTEGER      I,J,K,IOS,IR,NR,NRECL,NUMR
      REAL*4       AMN,AMX,RNUMR
      REAL         SCALE,FLAG
      REAL         RP(KDM),PP(KDM+1)
      REAL         RZ(KZ),ZZ(KZ),PZ(KZ+1)
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
      ITYPE_I = MOD(ITYPE,  10)
      ITYPE_P = MOD(ITYPE, 100)/10
      ITYPE_Z = MOD(ITYPE,1000)/100
C
      INQUIRE( IOLENGTH=NRECL) AP(:,:,1),PAD
#ifdef CRAY
#ifdef t3e
      IF     (MOD(NRECL,4096).EQ.0) THEN
        WRITE(CASN,8000) NRECL/4096
 8000   FORMAT('-F cachea:',I4.4,':1:0')
        IU8 = 11
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          write(6,*) 'Error: can''t asnunit 11'
          write(6,*) 'ios  = ',ios8
          write(6,*) 'casn = ',casn
          CALL EXIT(5)
        ENDIF
        IU8 = 12
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          write(6,*) 'Error: can''t asnunit 12'
          write(6,*) 'ios  = ',ios8
          write(6,*) 'casn = ',casn
          CALL EXIT(5)
        ENDIF
        IU8 = 21
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          write(6,*) 'Error: can''t asnunit 21'
          write(6,*) 'ios  = ',ios8
          write(6,*) 'casn = ',casn
          CALL EXIT(5)
        ENDIF
      ENDIF
#else
      CALL ASNUNIT(11,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 11'
        write(6,*) 'ios = ',ios
        CALL EXIT(5)
      ENDIF
      CALL ASNUNIT(12,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 12'
        write(6,*) 'ios = ',ios
        CALL EXIT(5)
      ENDIF
      CALL ASNUNIT(21,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 21'
        write(6,*) 'ios = ',ios
        CALL EXIT(5)
      ENDIF
#endif
#endif
      OPEN(UNIT=11, FILE=CFILE1, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILE1)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
      OPEN(UNIT=12, FILE=CFILEP, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEP)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
      OPEN(UNIT=13, FILE=CFILEZ, FORM='FORMATTED', STATUS='OLD',
     +     IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',TRIM(CFILEZ)
        WRITE(6,*) 'ios   = ',ios
        CALL EXIT(3)
      ENDIF
      OPEN(UNIT=21, FILE=CFILEO, FORM='UNFORMATTED', STATUS='NEW',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEO)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
C     READ THE ZI FILE.
C
      IF     (ITYPE_Z.EQ.0) THEN  !interface depths in m
        PZ(1) = 0.0
        DO K= 2,KZ+1
          READ(13,*) PZ(K)
          IF     (K.EQ.2 .AND. PZ(K).LE.0.0) THEN
            WRITE(6,*) 'Error: inconsistent input profile'
            WRITE(6,*) '1st value must be > 0.0'
            WRITE(6,*) PZ(2)
            CALL EXIT(7)
          ELSEIF (PZ(K).LT.PZ(K-1)) THEN
            WRITE(6,*) 'Error: inconsistent input profile'
            WRITE(6,*) k-1,'-th value < previous value'
            WRITE(6,*) PZ(2:K)
            CALL EXIT(8)
          ENDIF
        ENDDO
        CLOSE(13)
C
C       CONVERT TO PRESSURE UNITS
C
        DO K= 2,KZ+1
          PZ(K) = 9806.0*PZ(K)
        ENDDO
      ELSEIF (ITYPE_Z.EQ.4) THEN  !cell centers in m
        DO K= 1,KZ
          READ(13,*) ZZ(K)
          IF     (K.EQ.2 .AND. ZZ(K).LT.0.0) THEN
            WRITE(6,*) 'Error: inconsistent input profile'
            WRITE(6,*) '1st value must be >= 0.0'
            WRITE(6,*) ZZ(2)
            CALL EXIT(7)
          ELSEIF (ZZ(K).LT.ZZ(K-1)) THEN
            WRITE(6,*) 'Error: inconsistent output profile'
            WRITE(6,*) k-1,'-th value < previous value'
            WRITE(6,*) ZZ(2:K)
            CALL EXIT(8)
          ENDIF
        ENDDO
        CLOSE(13)
        DO K= 1,KZ
          ZZ(K) = 9806.0*ZZ(K)
        ENDDO
        PZ(1) = 0.0
        DO K= 1,KZ-1
          PZ(K+1) = 0.5*(ZZ(K)+ZZ(K+1))
          if     (itest.ne.0) then
            write(6,'(a,i3,2f10.3)') 
     &          'k,zz,pz = ',k+1,zz(k)/9806.0,pz(k+1)/9806.0
          endif
        ENDDO
        PZ(KZ+1) = 9806.0*9999.0
        if     (itest.ne.0) then
          write(6,'(a,i3,2f10.3)') 
     &        'k,zz,pz = ',kz+1,zz(kz)/9806.0,pz(kz+1)/9806.0
        endif
      ELSE  !illegal
        WRITE(6,*) 'Error: inconsistent itype (not 0xx or 4xx)'
        CALL EXIT(8)
      ENDIF !itype_z
C
C     READ IN AZ
C
      DO K= 1,KZ
        READ(11,REC=K,IOSTAT=IOS) AZ(:,:,K)
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(AZ(1,1,K),IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          IF     (K.EQ.1) THEN
              WRITE(6,*) 'can''t read ',TRIM(CFILE1)
              CALL EXIT(4)
          ELSE
            WRITE(6,*) TRIM(CFILE1),' is too short'
            CALL EXIT(4)
          ENDIF !k==1:else
        ENDIF !ios
      ENDDO !k
      CLOSE(11)
C
C     READ IN PNEW
C
      DO K= 1,KDM
        READ(12,REC=K,IOSTAT=IOS) P(:,:,K)
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(P(1,1,K),IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          IF     (K.EQ.1) THEN
              WRITE(6,*) 'can''t read ',TRIM(CFILEP)
              CALL EXIT(4)
          ELSE
            WRITE(6,*) TRIM(CFILEP),' is too short'
            CALL EXIT(4)
          ENDIF !k==1:else
        ENDIF !ios
      ENDDO !k
      CLOSE(12)
C
C     REMAP (DEPTHS IN PRESSURE UNITS)
C
      IF     (ITYPE_P.EQ.1 .OR. ITYPE_P.EQ.3) THEN
        SCALE = 1.0    !already in pressure units
      ELSE
        SCALE = 9806.0 !g/thref
      ENDIF
      PP(1) = 0.0
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (P(I,J,1).NE.SPVAL) THEN
            DO K= 1,KZ
              RZ(K) = AZ(I,J,K)
            ENDDO
            IF     (ITYPE_P.EQ.0 .OR. ITYPE_P.EQ.2) THEN  !interface depths
              DO K= 1,KDM
                PP(K+1) =         SCALE*P(I,J,K)
              ENDDO
            ELSE  !layer thicknesses
              DO K= 1,KDM
                PP(K+1) = PP(K) + SCALE*P(I,J,K)
              ENDDO
            ENDIF !itype_p
            FLAG = 0.0
C           note the inversion of "p" and "z" in the layer2zi routines
            IF     (ITYPE_I.EQ.0) THEN
              CALL LAYER2ZI_PCM(RZ,PZ,KZ, RP,PP,KDM, FLAG)
            ELSEIF (ITYPE_I.EQ.1) THEN
              CALL LAYER2ZI_PLM(RZ,PZ,KZ, RP,PP,KDM, FLAG)
*           ELSEIF (ITYPE_I.EQ.2 .AND.
*    &              i.eq.itest .and. j.eq.jtest) then
*             CALL LAYER2ZI_LCC_DEBUG(RZ,ZZ,KZ, RP,PP,KDM, FLAG)
            ELSEIF (ITYPE_I.EQ.2) THEN
              CALL LAYER2ZI_LCC(RZ,ZZ,KZ, RP,PP,KDM, FLAG)
*           ELSEIF (ITYPE_I.EQ.3 .AND.
*    &              i.eq.itest .and. j.eq.jtest) then
*             CALL LAYER2ZI_ACC_DEBUG(RZ,ZZ,KZ, RP,PP,KDM, FLAG)
            ELSEIF (ITYPE_I.EQ.3) THEN
              CALL LAYER2ZI_ACC(RZ,ZZ,KZ, RP,PP,KDM, FLAG)
            ENDIF
            DO K= 1,KDM
              AP(I,J,K) = RP(K)
            ENDDO
c
c ---       debugging printout
c
            if     (i.eq.itest .and. j.eq.jtest) then
              if     (itype_i.lt.2) then
                write(6,'(a,i3,f10.3,10x,f10.3)')
     &            'k,pz,rz,pp,rp =',
     &             0,pz(1)/9806.0,pp(1)/9806.0
                do k= 1,kz !assume kdm<kz
                  if     (k.le.kdm) then
                    write(6,'(a,i3,4f10.3)')
     &                'k,pz,rz,pp,rp =',
     &                 k,pz(k+1)/9806.0,rz(k),pp(k+1)/9806.0,rp(k)
                  else
                    write(6,'(a,i3,2f10.3)')
     &                'k,pz,rz       =',
     &                 k,pz(k+1)/9806.0,rz(k)
                  endif
                enddo !k
              elseif (itype_i.eq.2) then
                write(6,'(a,i3,20x,f10.3)')
     &            'k,      pp    =',
     &             0,pp(1)/9806.0
                do k= 1,kz !assume kdm<kz
                  if     (k.le.kdm) then
                    write(6,'(a,i3,4f10.3)')
     &                'k,zz,rz,pp,rp =',
     &                 k,zz(k)/9806.0,rz(k),pp(k+1)/9806.0,rp(k)
                  elseif (k.le.kz) then
                    write(6,'(a,i3,2f10.3)')
     &                'k,zz,rz       =',
     &                 k,zz(k)/9806.0,rz(k)
                  endif
                enddo !k
              else
                do k= 1,kz !assume kdm<kz
                  if     (k.le.kdm) then
                    write(6,'(a,i3,4f10.3)')
     &                'k,zz,rz,pc,rp =',
     &                 k,zz(k)/9806.0,rz(k),
     &                 0.5*(pp(k)+pp(k+1))/9806.0,rp(k)
                  elseif (k.le.kz) then
                    write(6,'(a,i3,2f10.3)')
     &                'k,zz,rz       =',
     &                 k,zz(k)/9806.0,rz(k)
                  endif
                enddo !k
              endif !itype_i
            endif !test point
c
c ---       end debugging printout
c
          ELSE
            DO K= 1,KDM
              AP(I,J,K) = SPVAL
            ENDDO
          ENDIF
        ENDDO !j
      ENDDO !i
C
C     WRITE OUT AP.
C
      DO K= 1,KDM
        AMN =  SPVAL
        AMX = -SPVAL
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (AP(I,J,K).NE.SPVAL) THEN
              AMN = MIN( AMN, AP(I,J,K) )
              AMX = MAX( AMX, AP(I,J,K) )
            ENDIF
          ENDDO !i
        ENDDO !j
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(AP(1,1,K),IDM*JDM)
#endif
        WRITE(21,REC=K,IOSTAT=IOS) AP(:,:,K)
        WRITE(6,'(a,1p2g16.8)') 'min, max = ',AMN,AMX
        IF     (IOS.NE.0) THEN
           WRITE(6,*) 'can''t write ',TRIM(CFILEO)
           CALL EXIT(4)
        ENDIF !ios
      ENDDO !k
      CLOSE(21)
C
      RETURN
      END

      subroutine layer2zi_pcm(r,p,kk, rz,pz,kz, flag)
      implicit none
c
      integer kk,kz
      real    r( kk),p( kk+1),
     &        rz(kz),pz(kz+1),flag
c
c**********
c*
c  1) remap from one set of vertical cells to another.
c     method: piecewise constant across each input cell.
c             the output is the average of the interpolation
c             profile across each output cell.
c
c  2) input arguments:
c       r     - scalar field in p-layer space
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
      real    q,zb,zt,rzk
c
      if     (r(1).eq.flag) then
        do k= 1,kz
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
            rz(k) = rz(k-1)
          else
c
c           form layer averages.
c
            if     (p(lf).gt.zt) then
              WRITE(6,*) 'bad lf = ',lf
              stop
            endif
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
                rzk = rzk + q*r(l)
*               WRITE(6,*) 'l,q = ',l,q
              else
c
c               the input layer is partially inside the output layer
c
                q   = max(min(p(l+1),zb)-max(p(l),zt),thin)/(zb-zt)
                rzk = rzk + q*r(l)
*               WRITE(6,*) 'l,q = ',l,q
              endif
            enddo !l
            rz(k) = rzk
          endif
        enddo !k
      endif
      return
      end subroutine layer2zi_pcm

      subroutine layer2zi_plm(r,p,kk, rz,pz,kz, flag)
      implicit none
c
      integer kk,kz
      real    r( kk),p( kk+1),
     &        rz(kz),pz(kz+1),flag
c
c**********
c*
c  1) remap from one set of vertical cells to another.
c     method: piecewise linear across each input cell.
c             the output is the average of the interpolation
c             profile across each output cell.
c
c  2) input arguments:
c       r     - scalar field in p-layer space
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
      real    q,qc,zb,zc,zt,rzk
      real    rs(kk),pt(kk+1)
c
      if     (r(1).eq.flag) then
        do k= 1,kz
          rz(k) = flag  ! land
        enddo
      else
c ---   compute PLM slopes for input layers
        do k=1,kk
          pt(k)=max(p(k+1)-p(k),thin)
        enddo
        call plm(pt, r, rs, kk)
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
            rz(k) = rz(k-1)
          else
c
c           form layer averages.
c
            if     (p(lf).gt.zt) then
              WRITE(6,*) 'bad lf = ',lf
              stop
            endif
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
                rzk = rzk + q*r(l)
*               WRITE(6,*) 'l,q = ',l,q
              else
c
c               the input layer is partially inside the output layer
c               average of linear profile is its center value
c
                q   = max( min(p(l+1),zb)-max(p(l),zt), thin )/(zb-zt)
                zc  = 0.5*(min(p(l+1),zb)+max(p(l),zt))
                qc  = (zc-p(l))/pt(l) - 0.5
                rzk = rzk + q*(r(l) + qc*rs(l))
*               WRITE(6,*) 'l,q,qc = ',l,q,qc
              endif
            enddo !l
            rz(k) = rzk
          endif
        enddo !k
      endif
      return
      end subroutine layer2zi_plm

      subroutine plm(pt, r, rs,kk)
      implicit none
c
      integer kk
      real    r(kk),pt(kk), rs(kk)
c
c**********
c*
c  1) generate a monotonic PLM interpolation of a layered field
c
c  2) input arguments:
c       pt    - layer interface thicknesses (non-zero)
c       r     - scalar field in layer space
c       kk    - dimension of a  (number of layers)
c
c  3) output arguments:
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

      subroutine layer2zi_lcc(r,pc,kk, rz,pz,kz, flag)
      implicit none
c
      integer kk,kz
      real    r( kk),pc( kk),
     &        rz(kz),pz(kz+1),flag
c
c**********
c*
c  1) remap from one set of vertical cells to another.
c     method: linear between input cell centers.
c             the output is the average of the interpolation
c             profile across each output cell.
c
c  2) input arguments:
c       r     - scalar field in p-layer space
c       pc    - layer center depths (non-negative m)
c                 pc(   1) can be the surface
c                 pc(kk+1) can be the bathymetry
c       kk    - dimension of a  (number of  input layers)
c       pz    - target interface depths (non-negative m)
c                 pz(k+1) >= pz(k)
c       flag  - data void (land) marker
c       kz    - dimension of az (number of output layers)
c
c  3) output arguments:
c       rz    - scalar field in pz-layer space
c
c  4) except at data voids, must have:
c           pc(   1) >= zero
c           pc( l+1) >= pc(:,:,l)
c           pc(kk+1) == bathymetry
c           0 <= pz(k) <= pz(k+1)
c      output layers completely below the bathymetry inherit values
c      from the layer above.
c
c  5) Alan J. Wallcraft, Naval Research Laboratory, March 2005.
c*
c**********
c
      real       thin 
      parameter (thin=1.e-6)  ! minimum layer thickness (no division by 0.0)
c
      integer k,l,lf
      real    q,zb,zt,rzk,rb,rt
c
      if     (r(1).eq.flag) then
        do k= 1,kz
          rz(k) = flag  ! land
        enddo
      else
        lf=1
        zb=pz(1)
        do k= 1,kz
          zt = zb
          zb = pz(k+1)
*         WRITE(6,*) 'k,zt,zb = ',k,zt,zb
          if     (zb-zt.lt.thin .or. zt.ge.pc(kk)) then
c
c ---       thin or bottomed layer, values taken from layer above
c
            rz(k) = rz(k-1)
          else
c
c           form layer averages, treat pc(l) to pc(l+1) as a "layer".
c
            if     (pc(lf).gt.zt) then
              WRITE(6,*) 'bad lf = ',lf
              stop
            endif
            rzk = 0.0
            do l= lf,kk-1
              if     (pc(l).gt.zb) then
*               WRITE(6,*) 'l,lf= ',l,lf,l-1
                lf = l-1
                exit
              elseif (pc(l).ge.zt .and. pc(l+1).le.zb) then
c
c               the input layer is completely inside the output layer
c
                q   = max(pc(l+1)-pc(l),thin)/(zb-zt)
                rzk = rzk + q*0.5*(r(l)+r(l+1))
*               WRITE(6,*) 'l,q = ',l,q
              else
c
c               the input layer is partially inside the output layer
c
                q   = max(min(pc(l+1),zb)-max(pc(l),zt),thin)/(zb-zt)
                rt  = r(l)   + (r(l+1)-r(l))*
     &                         (max(pc(l)  ,zt)-pc(l)  )/(pc(l+1)-pc(l))
                rb  = r(l+1) + (r(l)-r(l+1))*
     &                         (min(pc(l+1),zb)-pc(l+1))/(pc(l)-pc(l+1))
                rzk = rzk + q*0.5*(rt+rb)
*               WRITE(6,*) 'l,qrt,rb = ',l,q,rt,rb
              endif
            enddo !l
            rz(k) = rzk
          endif
        enddo !k
      endif
      return
      end subroutine layer2zi_lcc

      subroutine layer2zi_lcc_debug(r,pc,kk, rz,pz,kz, flag)
      implicit none
c
      integer kk,kz
      real    r( kk),pc( kk),
     &        rz(kz),pz(kz+1),flag
c
c**********
c*
c  1) remap from one set of vertical cells to another.
c     method: linear between input cell centers.
c             the output is the average of the interpolation
c             profile across each output cell.
c
c  2) input arguments:
c       r     - scalar field in p-layer space
c       pc    - layer center depths (non-negative m)
c                 pc(   1) can be the surface
c                 pc(kk+1) can be the bathymetry
c       kk    - dimension of a  (number of  input layers)
c       pz    - target interface depths (non-negative m)
c                 pz(k+1) >= pz(k)
c       flag  - data void (land) marker
c       kz    - dimension of az (number of output layers)
c
c  3) output arguments:
c       rz    - scalar field in pz-layer space
c
c  4) except at data voids, must have:
c           pc(   1) >= zero
c           pc( l+1) >= pc(:,:,l)
c           pc(kk+1) == bathymetry
c           0 <= pz(k) <= pz(k+1)
c      output layers completely below the bathymetry inherit values
c      from the layer above.
c
c  5) Alan J. Wallcraft, Naval Research Laboratory, March 2005.
c*
c**********
c
      real       thin 
      parameter (thin=1.e-6)  ! minimum layer thickness (no division by 0.0)
c
      integer k,l,lf
      real    q,zb,zt,rzk,rb,rt
c
      if     (r(1).eq.flag) then
        do k= 1,kz
          rz(k) = flag  ! land
        enddo
      else
        lf=1
        zb=pz(1)
        do k= 1,kz
          zt = zb
          zb = pz(k+1)
          WRITE(6,*) 'k,zt,zb = ',k,zt,zb
          if     (zb-zt.lt.thin .or. zt.ge.pc(kk)) then
c
c ---       thin or bottomed layer, values taken from layer above
c
            rz(k) = rz(k-1)
          else
c
c           form layer averages, treat pc(l) to pc(l+1) as a "layer".
c
            if     (pc(lf).gt.zt) then
              WRITE(6,*) 'bad lf = ',lf
              stop
            endif
            rzk = 0.0
            do l= lf,kk-1
              if     (pc(l).gt.zb) then
                WRITE(6,*) 'l,lf= ',l,lf,l-1
                lf = l-1
                exit
              elseif (pc(l).ge.zt .and. pc(l+1).le.zb) then
c
c               the input layer is completely inside the output layer
c
                q   = max(pc(l+1)-pc(l),thin)/(zb-zt)
                rzk = rzk + q*0.5*(r(l)+r(l+1))
                WRITE(6,*) 'l,q = ',l,q
              else
c
c               the input layer is partially inside the output layer
c
                q   = max(min(pc(l+1),zb)-max(pc(l),zt),thin)/(zb-zt)
                rt  = r(l)   + (r(l+1)-r(l))*
     &                         (max(pc(l)  ,zt)-pc(l)  )/(pc(l+1)-pc(l))
                rb  = r(l+1) + (r(l)-r(l+1))*
     &                         (min(pc(l+1),zb)-pc(l+1))/(pc(l)-pc(l+1))
                rzk = rzk + q*0.5*(rt+rb)
                WRITE(6,*) 'l,qrt,rb = ',l,q,rt,rb
              endif
            enddo !l
            rz(k) = rzk
          endif
        enddo !k
      endif
      return
      end subroutine layer2zi_lcc_debug

      subroutine layer2zi_acc(r,pc,kk, rz,pz,kz, flag)
      implicit none
c
      integer kk,kz
      real    r( kk),pc( kk),
     &        rz(kz),pz(kz+1),flag
c
c**********
c*
c  1) remap from one set of vertical cells to another.
c     method: linear between input cell centers.
c             sampled at output cell center (not conservative).
c
c  2) input arguments:
c       r     - scalar field in p-layer space
c       pc    - layer center depths (non-negative m)
c                 pc( 1) can be the surface
c                 pc(kk) can be the bathymetry
c       kk    - dimension of a  (number of  input layers)
c       pz    - target interface depths (non-negative m)
c                 pz(k+1) >= pz(k)
c       flag  - data void (land) marker
c       kz    - dimension of az (number of output layers)
c
c  3) output arguments:
c       rz    - scalar field in pz-layer space
c
c  4) except at data voids, must have:
c           0 <= pc(l) <= pc(l+1)
c           0 <= pz(k) <= pz(k+1)
c      output layers completely below the bathymetry inherit values
c      from the layer above.
c
c  5) Alan J. Wallcraft, Naval Research Laboratory, March 2005.
c*
c**********
c
      integer k,l,lf
      real    q,zb,zc,zt,rzk
c
      if     (r(1).eq.flag) then
        do k= 1,kz
          rz(k) = flag  ! land
        enddo
      else
        lf=1
        zb=pz(1)
        do k= 1,kz
          zt = zb
          zb = pz(k+1)
          zc = 0.5*(zb+zt)  !target cell center
*         WRITE(6,*) 'k,zt,zb,zc = ',k,zt,zb,zc
          if     (zc.gt.pc(kk)) then
c
c ---       bottomed layer, value taken from lowest input level
c
            rz(k) = r(kk)
          else
            if     (pc(lf).gt.zc) then
              WRITE(6,*) 'bad lf = ',lf
              stop
            endif
            do l= lf,kk-1
              if     (pc(l).le.zc .and. pc(l+1).ge.zc) then
c
c ---           the input layer spans the output cell center
c
                q   = (zc-pc(l))/(pc(l+1)-pc(l))
                rzk = r(l) + q*(r(l+1)-r(l))
*               WRITE(6,*) 'l,q = ',l,q
                lf = l
                exit
              endif
            enddo !l
            rz(k) = rzk
          endif
        enddo !k
      endif
      return
      end subroutine layer2zi_acc

      subroutine layer2zi_acc_debug(r,pc,kk, rz,pz,kz, flag)
      implicit none
c
      integer kk,kz
      real    r( kk),pc( kk),
     &        rz(kz),pz(kz+1),flag
c
c**********
c*
c  1) remap from one set of vertical cells to another.
c     method: linear between input cell centers.
c             sampled at output cell center (not conservative).
c
c  2) input arguments:
c       r     - scalar field in p-layer space
c       pc    - layer center depths (non-negative m)
c                 pc( 1) can be the surface
c                 pc(kk) can be the bathymetry
c       kk    - dimension of a  (number of  input layers)
c       pz    - target interface depths (non-negative m)
c                 pz(k+1) >= pz(k)
c       flag  - data void (land) marker
c       kz    - dimension of az (number of output layers)
c
c  3) output arguments:
c       rz    - scalar field in pz-layer space
c
c  4) except at data voids, must have:
c           0 <= pc(l) <= pc(l+1)
c           0 <= pz(k) <= pz(k+1)
c      output layers completely below the bathymetry inherit values
c      from the layer above.
c
c  5) Alan J. Wallcraft, Naval Research Laboratory, March 2005.
c*
c**********
c
      integer k,l,lf
      real    q,zb,zc,zt,rzk
c
      if     (r(1).eq.flag) then
        do k= 1,kz
          rz(k) = flag  ! land
        enddo
      else
        lf=1
        zb=pz(1)
        do k= 1,kz
          zt = zb
          zb = pz(k+1)
          zc = 0.5*(zb+zt)  !target cell center
          WRITE(6,*) 'k,zt,zb,zc = ',k,zt,zb,zc
          if     (zc.gt.pc(kk)) then
c
c ---       bottomed layer, value taken from lowest input level
c
            rz(k) = r(kk)
          else
            if     (pc(lf).gt.zc) then
              WRITE(6,*) 'bad lf = ',lf
              stop
            endif
            do l= lf,kk-1
              if     (pc(l).le.zc .and. pc(l+1).ge.zc) then
c
c ---           the input layer spans the output cell center
c
                q   = (zc-pc(l))/(pc(l+1)-pc(l))
                rzk = r(l) + q*(r(l+1)-r(l))
                WRITE(6,*) 'l,q = ',l,q
                lf = l
                exit
              endif
            enddo !l
            rz(k) = rzk
          endif
        enddo !k
      endif
      return
      end subroutine layer2zi_acc_debug

