      PROGRAM HYCOM_PERTURBATION
      IMPLICIT NONE
C
C  hycom_perturbation - Usage:  hycom_perturbation hscl nsample nfld [grid.a [mask.a]] fout.a
C  hycom_perturbation_1st - Usage:  hycom_perturbation_1st hscl nsample nfld [grid.a [mask.a]] fout.a
C
C                 creates a HYCOM .[ab] file containing scalar
C                 perturbation fields with given gaussian length scale
C
C                 a positive hscl is the gaussian length scale in m
C                 a negative hscl is minus the gaussian length scale in deg
C                 when converting from degrees to m, the minimum length
C                 is clipped half the maximum length scale
C
C                 nsample is the number of random impulse samples
C
C                 nfld is the number of scalar perturbation fields
C
C                 grid.a is a hycom grid file, default regional.grid.a.
C                 Note that the corresponding grid.b must also exist.
C
C                 mask.a contains an optional mask array, that turns
C                 off the field anywhere it is 2.0**100 (data void).
C                 Note that grid.a must be specified if mask.a is used.
C
C                 idm,jdm are taken from grid.a
C
C  For hycom_perturbation the range of each field is scaled to [-1.0,1.0].
C
C  For hycom_perturbation_1st the range of the 1st field is scaled to
C  [-1.0,1.0] and the same scale factor is then applied to all subsequent
C  fields.  Fields with an absmax of more than 3.0 are discarded.
C
C  the grid is assumed to be p-grid global with an arctic bi-polar patch.
C  always use mask.a for closed domains.
C
C  fout.a will contain idm*jdm 32-bit IEEE real values for the array,
C   in standard f77 element order, followed by padding to a multiple
C   of 4096 32-bit words, but otherwise with no control bytes/words,
C   and values of 2.0**100 indicating a data void.
C
C  the random number generator is initialized from the real time clock
C   unless the command name is hycom_perturbation_debug, which is always
C   initialized with the same seed.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  March 2013.
C  Based on a program by Matthew J. Carrier.
C
      REAL*4,  ALLOCATABLE :: R2D(:,:),AMSK(:,:),
     +                        PSCX(:,:),PSCY(:,:),HSCL(:,:)
      REAL*4               :: PAD(4096)
      INTEGER       IOS
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      LOGICAL       LDEBUG,LFIRST
      INTEGER       IDM,JDM,NPAD,NSAMPLE,NFLD
      REAL*4        HHSCL
      CHARACTER*6   CVARIN
      CHARACTER*240 CFILEG,CFILEM,CFILEB,CFILEO
C
C     READ ARGUMENTS.
C
      CALL GETARG(0,CARG)
      LFIRST = CARG.EQ.'hycom_perturbation_1st' .or.
     &         CARG.EQ.'hycom_perturbation_1st_debug'
      LDEBUG = CARG.EQ.'hycom_perturbation_debug' .or.
     &         CARG.EQ.'hycom_perturbation_1st_debug'
C
      NARG = IARGC()
C
      IF     (NARG.EQ.4) THEN
        CALL GETARG(1,CARG)
        READ(CARG,*)    HHSCL
        CALL GETARG(2,CARG)
        READ(CARG,*)    NSAMPLE
        CALL GETARG(3,CARG)
        READ(CARG,*)    NFLD
        CFILEG = 'regional.grid.a'
        CFILEM = 'NONE'
        CALL GETARG(3,CFILEO)
      ELSEIF (NARG.EQ.5) THEN
        CALL GETARG(1,CARG)
        READ(CARG,*)    HHSCL
        CALL GETARG(2,CARG)
        READ(CARG,*)    NSAMPLE
        CALL GETARG(3,CARG)
        READ(CARG,*)    NFLD
        CALL GETARG(4,CFILEG)
        CFILEM = 'NONE'
        CALL GETARG(5,CFILEO)
      ELSEIF (NARG.EQ.6) THEN
        CALL GETARG(1,CARG)
        READ(CARG,*)    HHSCL
        CALL GETARG(2,CARG)
        READ(CARG,*)    NSAMPLE
        CALL GETARG(3,CARG)
        READ(CARG,*)    NFLD
        CALL GETARG(4,CFILEG)
        CALL GETARG(5,CFILEM)
        CALL GETARG(6,CFILEO)
      ELSE
        WRITE(6,*) 'Usage:  hycom_perturbation '//
     +    'hscl nsample nfld [grid.a [mask.a]] fout.a'
        CALL EXIT(1)
      ENDIF
C
C     GET IDM,JDM FROM grid.b.
C
      CFILEB = CFILEG(1:LEN_TRIM(CFILEG)-1) // 'b'
C
      OPEN(UNIT=11,FILE=CFILEB,FORM='FORMATTED',
     &     STATUS='OLD',ACTION='READ')
C
      READ( 11,*) IDM,CVARIN
      IF (CVARIN.NE.'idm   ') THEN
        WRITE(6,*) 'hycom_perturbation: bad header file ',
     &             TRIM(CFILEB)
        CALL EXIT(2)
      ENDIF
      READ( 11,*) JDM,CVARIN
      IF (CVARIN.NE.'jdm   ') THEN
        WRITE(6,*) 'hycom_perturbation: bad header file ',
     &             TRIM(CFILEB)
        CALL EXIT(2)
      ENDIF
C
      CLOSE(UNIT=11)
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( R2D(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_perturbation: could not allocate ',
     +             IDM*JDM,' words for R2D'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( AMSK(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_perturbation: could not allocate ',
     +             IDM*JDM,' words for AMSK'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( PSCX(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_perturbation: could not allocate ',
     +             IDM*JDM,' words for PSCX'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( PSCY(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_perturbation: could not allocate ',
     +             IDM*JDM,' words for PSCY'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( HSCL(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_perturbation: could not allocate ',
     +             IDM*JDM,' words for HSCL'
        CALL EXIT(2)
      ENDIF
C
      CALL PERTURBATION(R2D,AMSK,PSCX,PSCY,HSCL,IDM,JDM,PAD,NPAD,
     +                  HHSCL,NSAMPLE,NFLD, CFILEG,CFILEM,CFILEO,
     &                  LDEBUG,LFIRST)
      CALL EXIT(0)
 5000 FORMAT(I4)
      END
      SUBROUTINE PERTURBATION(R2D,AMSK,PSCX,PSCY,HSCL,IDM,JDM, PAD,NPAD,
     +                        HHSCL,NSAMPLE,NFLD, CFILEG,CFILEM,CFILEO,
     +                        LDEBUG,LFIRST)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE,CFILEG,CFILEM,CFILEO
      LOGICAL       LDEBUG,LFIRST
      INTEGER       IDM,JDM,NPAD,NSAMPLE,NFLD
      REAL*4        R2D(IDM,JDM),AMSK(IDM,JDM),
     +              PSCX(IDM,JDM),PSCY(IDM,JDM),
     +              HSCL(IDM,JDM),PAD(NPAD),
     +              HHSCL,DEG2RAD
C
C     MOST OF WORK IS DONE HERE.
C
C
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      INTEGER, ALLOCATABLE :: NEW_SEED(:)
      CHARACTER*10 DATE, TIME
      CHARACTER*18 CASN
      INTEGER      DATES(8),SEED_SIZE
      INTEGER      I,II,IOS,J,K,KS,NRECL
      REAL         SCL,V,X,Y
      REAL*4       AMX,AMN
      REAL*8       FSUM,ASUM
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
C     INPUT MASK ARRAY A.
C
      INQUIRE( IOLENGTH=NRECL) AMSK,PAD
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
      ENDIF
#else
      CALL ASNUNIT(11,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 11'
        write(6,*) 'ios = ',ios
        CALL EXIT(5)
      ENDIF
#endif
#endif
C
      IF     (CFILEM.EQ.'NONE') THEN
        AMSK(:,:) = 1.0  !all land
      ELSE
        OPEN(UNIT=11, FILE=CFILEM, FORM='UNFORMATTED', STATUS='OLD',
     +           ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
        IF     (IOS.NE.0) THEN
          write(6,*) 'Error: can''t open ',TRIM(CFILEM)
          write(6,*) 'ios   = ',ios
          write(6,*) 'nrecl = ',nrecl
          CALL EXIT(3)
        ENDIF
C
        READ(11,REC=1,IOSTAT=IOS) AMSK
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(AMSK,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'can''t read record 1 of ',TRIM(CFILEM)
          CALL EXIT(4)
        ENDIF
C
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (AMSK(I,J).NE.SPVAL) THEN
              AMSK(I,J) = 1.0
            ELSE
              AMSK(I,J) = 0.0
            ENDIF
          ENDDO !i
        ENDDO !j
        CLOSE(UNIT=11)
      ENDIF
C
C     INPUT GRID ARRAYS.
C
      OPEN(UNIT=11, FILE=CFILEG, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEG)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
      READ(11,REC=10,IOSTAT=IOS) PSCX  ! pscx
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(PSCX,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read record 10 of ',TRIM(CFILEG)
        CALL EXIT(4)
      ENDIF
C
      READ(11,REC=11,IOSTAT=IOS) PSCY  ! pscy
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(PSCY,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read record 11 of ',TRIM(CFILEG)
        CALL EXIT(4)
      ENDIF
C
      IF     (HHSCL.GE.0.0) THEN
        HSCL(:,:) = HHSCL*AMSK(:,:)
        AMN = HHSCL
        AMX = HHSCL
      ELSE
        READ(11,REC=2,IOSTAT=IOS) HSCL  ! plat
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(HSCL,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'can''t read record 2 of ',TRIM(CFILEG)
          CALL EXIT(4)
        ENDIF
        DEG2RAD = 4.D0*ATAN(1.D0)/180.D0  !PI/180
        AMN =  SPVAL
        AMX = -SPVAL
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (AMSK(I,J).EQ.1.0) THEN
              HSCL(I,J) = -HHSCL * 111.2E3 * 
     &                      SQRT(ABS(COS(HSCL(I,J)*DEG2RAD)))
              AMX = MAX( AMX, HSCL(I,J) )
              AMN = MIN( AMN, HSCL(I,J) )
            ENDIF
          ENDDO !i
        ENDDO !j
        AMN = MAX( AMN, 0.5*AMX )  !clip any small values
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (AMSK(I,J).EQ.1.0) THEN
              HSCL(I,J) = MAX( HSCL(I,J), AMN )
              if     (i.eq.idm/4) then
                write(6,*) 'j,hscl = ',j,hscl(i,j),
     &                     hscl(i,j)/sqrt(pscx(i,j)*pscy(i,j))
              endif
            ELSE
              HSCL(I,J) = 0.0
            ENDIF
          ENDDO !i
        ENDDO !j
      ENDIF !hhscl>0:else
        write(6,*) 'hscl = ',AMN,AMX
        call flush(6)
C
      CLOSE(UNIT=11)
C
C     OUTPUT FILE.
C
      OPEN(UNIT=21, FILE=CFILEO, FORM='UNFORMATTED', STATUS='NEW',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEO)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
C     LOOP THROUGH ALL PERTURBATION FIELDS
C
      CALL RANDOM_SEED(SIZE=SEED_SIZE) ! Determine seed array size.
      ALLOCATE(NEW_SEED(SEED_SIZE))     ! Allocate the seed array.
      NEW_SEED(:) = [(k,k=1,SEED_SIZE)]
      IF     (.NOT.LDEBUG) THEN
C       SEED FROM REAL TIME CLOCK
        CALL DATE_AND_TIME(DATE,TIME,VALUES=DATES)
        NEW_SEED(:)=NEW_SEED(:)+
     &              DATES(1)+DATES(2)+DATES(3)+
     &              DATES(5)+DATES(6)+DATES(7)+DATES(8)
      ENDIF
        WRITE(6,*) 'new_seed=',NEW_SEED(:),' is of size=',SEED_SIZE
      CALL RANDOM_SEED(PUT=NEW_SEED)    ! set Random_number with the new seed.
C
      K = 1
      DO 
C
C       CREATE NSAMPLE NON-OVERLAPPING IMPULSES OF RANDOM STRENGTH
C
        R2D(:,:) = 0.0
        KS = 0
        DO
          CALL RANDOM_NUMBER(X)
          CALL RANDOM_NUMBER(Y)
          CALL RANDOM_NUMBER(V)
          I = NINT(0.5001 + (IDM-0.5001)*X)  !1:IDM
          J = NINT(0.5001 + (JDM-1.5001)*Y)  !1:JDM-1
*              WRITE(6,*) 'i,j,y = ',i,j,y
*              CALL FLUSH(6)
          IF     (AMSK(I,J).NE.0.0 .AND. R2D(I,J).EQ.0.0) THEN
*              WRITE(6,*) 'i,j,v = ',i,j,v-0.5,KS+1
*              CALL FLUSH(6)
            R2D(I,J) = V - 0.5
            KS = KS + 1
            IF     (KS.GT.NSAMPLE) THEN
              EXIT
            ENDIF
          ENDIF
        ENDDO !ks
        DO I= 1,IDM
          II = IDM-MOD(I-1,IDM)
          R2D(I,JDM) = R2D(II,JDM-1)
        ENDDO !i
C
C       ERROR CORRELATION MODEL BASED ON THE IMPLICIT SOLUTION
C       OF A DIFFUSION EQUATION
C
        CALL SP_CORR2D(IDM,JDM,AMSK,PSCX,PSCY,HSCL,R2D)
C
C       SCALE TO APPROXIMATE RANGE -1.0:1.0
C
        IF     (LFIRST) THEN
          IF     (K.EQ.1) THEN
            SCL = 1.0/MAX( ABS(MINVAL(R2D(:,:))),
     &                     ABS(MAXVAL(R2D(:,:))) )
          ENDIF !k==1
        ELSE
          SCL = 1.0/MAX( ABS(MINVAL(R2D(:,:))),
     &                   ABS(MAXVAL(R2D(:,:))) )
        ENDIF !lfirst:else
C
        AMN =  SPVAL
        AMX = -SPVAL
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (AMSK(I,J).NE.0.0) THEN
              R2D(I,J) = SCL*R2D(I,J)
              AMX = MAX( AMX, R2D(I,J) )
              AMN = MIN( AMN, R2D(I,J) )
            ELSE
              R2D(I,J) = SPVAL
            ENDIF
          ENDDO !j
        ENDDO !i
C
        IF     (MAX(ABS(AMN),ABS(AMX)).GT.3.0) THEN
          WRITE(6,'(a,1p3g16.8,a)')
     &       'scl, min, max = ',SCL,AMN,AMX," *DISCARDED*"
          CYCLE !redo k
        ENDIF
C
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(R2D,IDM*JDM)
#endif
        WRITE(21,REC=K,IOSTAT=IOS) R2D
        WRITE(6,'(a,1p3g16.8)')
     &     'scl, min, max = ',SCL,AMN,AMX
        CALL FLUSH(6)
        K = K + 1
        IF     (K.GT.NFLD) THEN
          EXIT
        ENDIF
      ENDDO !do forever
      WRITE(6,*)
      WRITE(6,*) K-1,' FIELDS PROCESSED'
      WRITE(6,*)
C
      CLOSE(21)
C
      RETURN
      END
      subroutine sp_corr2d(n,m,msk,dx,dy,hscl,r2d)
      implicit none
c
      integer     n,m
      real        msk( n,m)
      real        dx(  n,m)
      real        dy(  n,m)
      real        hscl(n,m)
      real        r2d( n,m)
c
c --- the error correlation model based on the implicit solution
c --- of a diffusion equation (properly normalized by norm2d)
c --- 30 June, 2008 (Matthew J. Carrier)
c --- Updated for HYCOM arctic bi-polar patch grid, March 2013.
c
c --- solve x.0 = r2d; x.n+1 = x.n + DT K DEL.SQ x.n+1 
c --- using preconditioned congugate gradients:
c --- http://en.wikipedia.org/wiki/Conjugate_gradient_method
c --- http://en.wikipedia.org/wiki/Conjugate_gradient_method#The_preconditioned_conjugate_gradient_method
c
c local parameter list
c
      real,        parameter :: dtol   = 2.0
      real,        parameter :: p5     = 0.5
      real,        parameter :: zero   = 0.0
      integer,     parameter :: itermn =  10
      integer,     parameter :: itermx = 900
      integer,     parameter :: Ntcov  = 2
c
c CG variables
c
      integer itr,itrst,n_it
      real dtcov(n,m)
      real ax(   n,m)
      real ay(   n,m)
      real zz(   n,m)
      real xx(   n,m)
      real rr(   n,m)
      real pp(   n,m)
      real qq(   n,m)
      real bb
      real bbr(  n,m)
      real ppc,ppe,ppn,pps,ppw
      real rtz1,ptq1,cgs,rtol
      real rtz,ptq,alpha,beta,rtz_old,a,cg_sum
      real rrnorm,rrtol,grad,grad_test,descent
      integer i,ii,j,k,im,ip,jm,jp,p
c
c compute the necessary matrix elements and time step
c
      do j=1,m
      do i=1,n
         ax(i,j)=1./(dx(i,j)*dx(i,j))
         ay(i,j)=1./(dy(i,j)*dy(i,j))
         dtcov(i,j) = (0.25*(hscl(i,j)**2))
      enddo
      enddo 
c
c calculate diagonal for preconditioner
c
      bbr(:,:) = 0.0
      do j=1,m-1
        jp=j+1
        do i=1,n
          if(msk(i,j).eq.1.0) then
            if     (i.ne.n) then
              ip=i+1
            else
              ip=1
            endif
            bb = 1.0+dtcov(i,j)*(ax(ip,j)+ax(i,j)+ay(i,jp)+ay(i,j))
            bbr(i,j) = 1.0/bb
          endif
        enddo !i
      enddo !j
c ---   arctic patch
        do i= 1,n
          ii = n-mod(i-1,n)
          bbr(i,m) = bbr(ii,m-1)
        enddo !i
*       write(6,*) 'bbr = ',minval(bbr),maxval(bbr)
        call flush(6)
c
c===============================================c
c         PRECONDITIONED CG GRADIENT LOOP       c
c===============================================c
c
      do p=1,Ntcov
         xx(:,:) = 0.0  !zero 1st guess simplifies the initial pass
         rr(:,:) = r2d(:,:)
         rrtol   = 0.1*maxval(abs(rr(:,:)))
         zz(:,:) = rr(:,:)*bbr(:,:)
         pp(:,:) = zz(:,:)
         rtz = sum(rr(:,:)*zz(:,:))
*           write(6,*) 'rrtol  = ',rrtol
*           write(6,*) 'xx  = ',minval(xx(:,:)),maxval(xx(:,:))
*           write(6,*) 'rr  = ',minval(rr(:,:)),maxval(rr(:,:))
*           write(6,*) 'zz  = ',minval(qq(:,:)),maxval(qq(:,:))
*           write(6,*) 'pp  = ',minval(pp(:,:)),maxval(pp(:,:))
*           call flush(6)
         itr  = 0
         do !CG loop
            itr  = itr + 1
c ---       qq = A * pp
            do j=1,m-1
              jm=j-1; jm=max(1,jm)
              jp=j+1
              do i=1,n
                if(msk(i,j).eq.1.0) then
                  if     (i.ne.1) then
                    im=i-1
                  else
                    im=n
                  endif
                  if     (i.ne.n) then
                    ip=i+1
                  else
                    ip=1
                  endif
                  ppc = pp(i, j)
                  if (msk(im,j).eq.1.0) then
                    ppw = pp(im,j)
                  else
                    ppw = pp(i ,j)
                  endif
                  if (msk(ip,j).eq.1.0) then
                    ppe = pp(ip,j)
                  else
                    ppe = pp(i ,j)
                  endif
                  if (msk(i,jm).eq.1.0) then
                    pps = pp(i,jm)
                  else
                    pps = pp(i,j)
                  endif
                  if (msk(i,jp).eq.1.0) then
                    ppn = pp(i,jp)
                  else
                    ppn = pp(i,j)
                  endif
                  qq(i,j)=ppc+dtcov(i,j)*(  ax(ip,j) *(ppc-ppe)
     &                                     +ax(i, j) *(ppc-ppw)
     &                                     +ay(i, jp)*(ppc-ppn)
     &                                     +ay(i, j) *(ppc-pps) )
                else
                  qq(i,j) = 0.0
                endif
              enddo !i
            enddo !j
            do i= 1,n !arctic patch
              ii = n-mod(i-1,n)
              qq(i,m) = qq(ii,m-1)
            enddo !i
c ---       qq = A * pp --- END
*              write(6,*) 'qq  = ',minval(qq(:,:)),maxval(qq(:,:))
*              call flush(6)
c
            ptq = sum(pp(:,:)*qq(:,:))
            if(ptq.ne.zero) then
               alpha = rtz/ptq
            else
               alpha = zero
            endif
            xx(:,:) = xx(:,:) + alpha*pp(:,:)
            rr(:,:) = rr(:,:) - alpha*qq(:,:)
            do i= 1,n !arctic patch
              ii = n-mod(i-1,n)
              xx(i,m) = xx(ii,m-1)
              rr(i,m) = rr(ii,m-1)
            enddo !i
*              write(6,*) 'xx  = ',minval(xx(:,:)),maxval(xx(:,:))
*              write(6,*) 'rr  = ',minval(rr(:,:)),maxval(rr(:,:))
*              call flush(6)
c
c ---       convergence test
            if     (itr.ge.itermn) then
              rrnorm   = maxval(abs(rr(:,:)))
*                write(6,*) 'rrnrom = ',rrnorm
*                call flush(6)
              if(itr.gt.itermx .or. rrnorm.le.rrtol) then
                 r2d(:,:) = xx(:,:)*msk(:,:)
                 exit
              endif !rrnorm
            endif !itermn
c
            zz(:,:) = rr(:,:)*bbr(:,:)
            rtz_old = rtz
            rtz     = sum(rr(:,:)*zz(:,:))
            beta    = rtz/rtz_old
            pp(:,:) = zz(:,:)+beta*pp(:,:)
c
*              write(6,*) 'zz  = ',minval(zz(:,:)),maxval(zz(:,:))
*              write(6,*) 'pp  = ',minval(pp(:,:)),maxval(pp(:,:))
*              call flush(6)
c
c end CG loop
c
         enddo !CG loop
         if(itr.gt.itermx) then
            write(*,*)'WARNING:  Spatial correlation operator reached' 
            write(*,*)'maximum iterations.  Check values to ensure'
            write(*,*)'code is funcitioning correctly'
         endif
c
c end loop on Ntcov
c
      enddo !Ntcov
c
c program termination
c
      return
      end
