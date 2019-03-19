      PROGRAM HYCOM_PROFILE_LAYPRS
      IMPLICIT NONE
C
C  hycom_profile_layprs - Usage: hycom_profile_layprs archva.txt archvb.txt layprs.txt [flag]
C
C                 calculates the NCODA-style displacement of isopycnals
C                 between two with the same number of layers.
C
C   archva.txt is assumed to be an HYCOM archive text profile file
C   archvb.txt is assumed to be an HYCOM archive text profile file
C   layprs.txt will be the output displacement in meters from a to b
C   flag selects options (default 0)
C     =0; no smooth, and no limiter
C     =1; no smooth, and    limiter
C     =2;    smooth, and    limiter
C   under the limiter option the displacements are limited as done in
C   ncoda_archv_inc
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  May 2012.
C
      INTEGER      IARGC
      INTEGER      NARG
      CHARACTER*240 CARG
C
      CHARACTER*240 CFILEA,CFILEB,CFILEC,CFORMAT
      CHARACTER*240 CLINE
      REAL          HYBISO,NOHISO,THK,DEPTH,DUM5(5),FLAG,
     &              DZMAX,DZMIN
      INTEGER       IOS,K,KI,KK,KDM,KO,KP,KT
      INTEGER       I,IFLAG,KMAX
      LOGICAL       LIMITER,LSMOOTH
C
      REAL, ALLOCATABLE :: SI(:,:),PI(:),ZI(:),
     &                     SO(:,:),PO(:),ZO(:),DZ(:),DS(:)
C
      INTEGER       NSAMP
      REAL          PMAX
      PARAMETER(NSAMP=1,PMAX=500.0)
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.3) THEN
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CFILEB)
        CALL GETARG(3,CFILEC)
        IFLAG = 0
      ELSEIF (NARG.EQ.4) THEN
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CFILEB)
        CALL GETARG(3,CFILEC)
        CALL GETARG(4,CARG)
        READ(CARG,*) IFLAG
      ELSE
        WRITE(6,"(a,a)")
     +    'Usage: hycom_profile_layprs archva.txt archvb.txt ',
     +                                'layprs.txt [flag]'
        CALL EXIT(1)
      ENDIF
      LIMITER = IFLAG.GT.0
      LSMOOTH = IFLAG.EQ.2
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
      DO K= 1,9
        READ( 11,'(a)') CLINE
        IF     (CLINE(1:4).EQ.'#  k') THEN
          EXIT
        ENDIF
        WRITE(21,'(a)') TRIM(CLINE)
      ENDDO
C
C     READ 1ST THE ISOPYCNAL PROFILE, TO GET KDM.
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
C     RE-READ THE 1ST ISOPYCNAL PROFILE.
C
      ALLOCATE( PI(KDM+1), ZI(KDM), SI(KDM,5) )
C
      REWIND(11)
      DO K= 1,9
        READ(11,'(a)') CLINE
        IF     (CLINE(1:4).EQ.'#  k') THEN
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
        READ(CLINE,*) KI,(SI(K,KK),KK=1,5),THK,DEPTH
        PI(K+1) = PI(K) + THK
        ZI(K)   = PI(K) + THK*0.5
        IF     (THK.EQ.0.0) THEN
          DO KK= 1,5
            SI(K,KK)=SI(K-1,KK)
          ENDDO !kk
        ENDIF
      ENDDO
      CLOSE(11)
C
C     READ THE 2ND ISOPYCNAL PROFILE.
C
      ALLOCATE( PO(KDM+1), ZO(KDM), SO(KDM,5), DZ(KDM), DS(KDM) )
C
      REWIND(12)
      DO K= 1,9
        READ(12,'(a)') CLINE
        IF     (CLINE(1:4).EQ.'#  k') THEN
          EXIT
        ENDIF
      ENDDO
      PO(1) =  0.0
      DO K= 1,KDM
        READ(12,'(a)',IOSTAT=IOS) CLINE
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'Error: inconsistent input profile'
          CALL EXIT(6)
        ENDIF
        READ(CLINE,*) KI,(SO(K,KK),KK=1,5),THK,DEPTH
        PO(K+1) = PO(K) + THK
        ZO(K)   = PO(K) + THK*0.5
        IF     (THK.EQ.0.0) THEN
          DO KK= 1,5
            SO(K,KK)=SO(K-1,KK)
          ENDDO !kk
        ENDIF
        SO(K,5) = MAX( SO(K,5), SO(MAX(K-1,1),5) )  !stable profile
      ENDDO
      CLOSE(12)
C
C     CALCULATE THE DISPLACEMENTS
C
      CALL LAYPRS(SI(1,5),ZI,
     &            SO(1,5),ZO, DZ, KDM, PI(KDM+1))
C
C     APPLY THE SMOOTHER (TWICE)?
C
      IF     (LSMOOTH) THEN
        DS(:) = DZ(:)
        DO K= 2,KDM-1
          DZ(K) = 0.5*DS(K) + 0.25*(DS(K-1)+DS(K+1))
*                       write(6,*) 'k,ds,dz  = ',k,ds(k),dz(k)
        ENDDO !k
        DS(:) = DZ(:)
        DO K= 2,KDM-1
          DZ(K) = 0.5*DS(K) + 0.25*(DS(K-1)+DS(K+1))
*                       write(6,*) 'k,ds,dz  = ',k,ds(k),dz(k)
        ENDDO !k
      ENDIF
C
C     APPLY LIMITER?
C
      IF     (LIMITER) THEN
        DO K= 2,KDM
          DZMIN = (ZI(K-1)-ZI(K)) !    negative
          DZMAX = (ZI(KDM)-ZI(K)) !non-negative
*                       write(6,*) 'k,dz     = ',k,dz(k)
*                       write(6,*) 'k,dzmin  = ',k,dz(k-1)+dzmin
*                       write(6,*) 'k,dzmax  = ',k,dzmax
C         surface moves from zi(k) to zi(k)+dz(k), which
C         must be below zi(k-1)+dz(k-1) but above zi(kdm)
          DZ(K) = MAX( DZ(K),
     &                 DZ(K-1)+DZMIN )
          DZ(K) = MIN( DZ(K),  DZMAX )
*                       write(6,*) 'k,dz_lim = ',k,dz(k)
        ENDDO !k
      ENDIF
C
C     OUTPUT THE DISPLACEMENTS
C
      WRITE(21,'(4a)')
     &      '#  k',
     &      ' p.densA p.densB',
     &      '    layprs      zprs',
     &      '        zA        zB'
      DO K= 1,KDM
        WRITE(21,'(i4,2f8.3,2f10.3,2f10.3)')
     &    K,
     &    SI(K,5),SO(K,5),
     &    DZ(K),ZI(K)+DZ(K),
     &    ZI(K),ZO(K)
      ENDDO !k
      CLOSE(21)
      END

      subroutine layprs(ri,zi,ro,zo,dz,kk, depth)
      implicit none
c
      integer kk
      real    ri(kk),zi(kk),
     &        ro(kk),zo(kk),dz(kk),depth
c
c**********
c*
c  1) calculate the isopycnal displacement between two density profiles.
c
c  2) input arguments:
c       ri    - 1st density profile values
c       zi    - 1st density profile depths
c       ro    - 2nd density profile values
c       zo    - 2nd density profile depths
c       kk    - number of levels
c       depth - maximum depth, can be less than zi(kk) and zo(kk)
c
c  3) output arguments:
c       dz    - displacement between the two density profiles
c
c  4) assumes that ro is a non-decreasing profile.
c     zi(:)+dz(:) will be non-decreasing, between zo(1) and depth.
c     ro.z at z=zi(k)+dz(k) should be ri(k), but this will only be
c     the case where ri is locally non-decreasing.
c
c  5) Alan J. Wallcraft,  Naval Research Laboratory,  May 2012.
c*
c**********
c
      real       thin 
      parameter (thin=1.e-6)  ! minimum layer thickness
c
      integer k,ko,kom1
      real    q,z,zold
c
      kom1 = 1
      zold = zo(1)
      do k= 1,kk
c ---   find the location of ri(k) in ro(:).
        if     (ro(kom1).ge.ri(k)) then
          dz(k) = zold - zi(k)   !zold and kom1 unchanged
        else
          do ko= kom1+1,kk
            if     (ro(ko).ge.ri(k)) then  !also ro(ko-1).lt.ri(k)
              q     = (ro(ko) - ri(k))/max(ro(ko) - ro(ko-1), thin)
              z     = q*zo(ko-1) + (1.0-q)*zo(ko)
              z     = min( z, depth )
              zold  = max( z, zold  )
              dz(k) = zold - zi(k)
              kom1  = ko-1
              exit
            endif
          enddo !ko
        endif !ro.kom1:else
      enddo !k
      return
      end subroutine layprs
