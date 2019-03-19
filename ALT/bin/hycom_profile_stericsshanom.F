      PROGRAM HYCOM_PROFILE_STERICSSHANOM
      IMPLICIT NONE
C
C  hycom_profile_stericsshanom - Usage: hycom_profile_stericsshanom archva.txt archvb.txt z.txt [zref] sssha.txt [sssha_t.txt sssha_s.txt]
C
C                 calculates the density anomaly and 
C                 the cumalative steric SSH anomaly at fixed depths
C                 between two profiles with the same number of layers.
C                 if zref is present also write out sssha(z)/sssha(zref).
C
C   archva.txt is assumed to be an HYCOM archive text profile file
C   archvb.txt is assumed to be an HYCOM archive text profile file
C        z.txt is the required depths in m
C         zref is a  reference depth  in m
C    sssha.txt will be the output steric SSH anomaly vs depth
C  sssha_t.txt will be the output steric SSH anomaly vs depth from T only
C  sssha_s.txt will be the output steric SSH anomaly vs depth from S only
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  Dec. 2015.
C
      INTEGER      IARGC
      INTEGER      NARG
      CHARACTER*240 CARG
C
      CHARACTER*240 CFILEA,CFILEB,CFILEZ,CFILEC,CFILET,CFILES,CFORMAT
      CHARACTER*240 CLINE
      REAL          HYBISO,NOHISO,THK,DEPTH,DUM5(5),FLAG,AREF,ZREF
      REAL          SSHI,SSHO
      INTEGER       IOS,K,KI,KK,KDM,KO,KP,KT,NZ,SIGVER
      INTEGER       I,KMAX
C
      REAL, ALLOCATABLE ::  SI(:,:),PI(:),
     &                      SO(:,:),PO(:),
     &                     SSSHA(:),ZZ(:),RA(:)
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
        CALL GETARG(3,CFILEZ)
        CALL GETARG(4,CFILEC)
        ZREF   = 0
        CFILET = ' '
        CFILES = ' '
      ELSEIF (NARG.EQ.5) THEN
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CFILEB)
        CALL GETARG(3,CFILEZ)
        CALL GETARG(4,CLINE)
        READ(CLINE,*) ZREF
        CALL GETARG(5,CFILEC)
        CFILET = ' '
        CFILES = ' '
      ELSEIF (NARG.EQ.6) THEN
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CFILEB)
        CALL GETARG(3,CFILEZ)
        CALL GETARG(4,CFILEC)
        CALL GETARG(5,CFILET)
        CALL GETARG(6,CFILES)
        ZREF   = 0
      ELSEIF (NARG.EQ.7) THEN
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CFILEB)
        CALL GETARG(3,CFILEZ)
        CALL GETARG(4,CLINE)
        READ(CLINE,*) ZREF
        CALL GETARG(5,CFILEC)
        CALL GETARG(6,CFILET)
        CALL GETARG(7,CFILES)
      ELSE
        WRITE(6,"(2a)")
     &    'Usage: hycom_profile_stericsshanom archva.txt archvb.txt ',
     &             'z.txt [zref] sssha.txt [sssha_t.txt sssha_s.txt]'
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
      OPEN(UNIT=13, FILE=CFILEZ, FORM='FORMATTED', STATUS='OLD',
     +     IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',TRIM(CFILEZ)
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
      IF     (CFILET.NE.' ') THEN
        OPEN(UNIT=22, FILE=CFILET, FORM='FORMATTED', STATUS='NEW',
     +       IOSTAT=IOS)
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'Error: can''t open ',TRIM(CFILET)
          WRITE(6,*) 'ios   = ',ios
          CALL EXIT(5)
        ENDIF
        OPEN(UNIT=23, FILE=CFILES, FORM='FORMATTED', STATUS='NEW',
     +       IOSTAT=IOS)
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'Error: can''t open ',TRIM(CFILES)
          WRITE(6,*) 'ios   = ',ios
          CALL EXIT(5)
        ENDIF
      ENDIF
C
C     READ 1ST THE ISOPYCNAL PROFILE, TO GET KDM.
C
      DO K= 1,9
        READ(11,'(a)') CLINE
        IF     (CLINE(1:4).EQ.'#  k') THEN
          EXIT
        ENDIF
      ENDDO
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
        IF     (CLINE(1:1).NE.'#') THEN
          READ(CLINE,*) KDM
        ENDIF
      ENDDO
C
C     RE-READ THE 1ST ISOPYCNAL PROFILE.
C
      ALLOCATE( PI(KDM+1), SI(KDM,5) )
C
      REWIND(11)
      DO K= 1,9
        READ(11,'(a)') CLINE
        IF     (K.EQ.4) THEN
          CLINE = CLINE(13:)
          READ(CLINE,*) SSHI
          SSHI = SSHI*0.01 !cm -> m
        ELSEIF (CLINE(1:4).EQ.'#  k') THEN
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
        IF     (THK.EQ.0.0) THEN
          DO KK= 1,5
            SI(K,KK)=SI(K-1,KK)
          ENDDO !kk
        ENDIF
      ENDDO
      CLOSE(11)
      IF     (SI(KDM,5).LT.30.0) THEN
        SIGVER = 1  ! 7-term sigma-0
      ELSE
        SIGVER = 6  !17-term sigma-2
      ENDIF
      CALL SIG_I(SIGVER)
C
C     READ THE 2ND ISOPYCNAL PROFILE.
C
      ALLOCATE( PO(KDM+1), SO(KDM,7) )
C
      REWIND(12)
      DO K= 1,9
        READ(12,'(a)') CLINE
        IF     (K.EQ.4) THEN
          CLINE = CLINE(13:)
          READ(CLINE,*) SSHO
          SSHO = SSHO*0.01 !cm -> m
        ELSEIF (CLINE(1:4).EQ.'#  k') THEN
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
        IF     (THK.EQ.0.0) THEN
          DO KK= 1,5
            SO(K,KK)=SO(K-1,KK)
          ENDDO !kk
        ENDIF
        SO(K,5) = MAX( SO(K,5), SO(MAX(K-1,1),5) )  !stable profile
        IF     (CFILET.NE.' ') THEN
          CALL SIG_P(SO(K,3),SI(K,4), SO(K,6))  !new  T orig S
          CALL SIG_P(SI(K,3),SO(K,4), SO(K,7))  !orig T new  S
          SO(K,6) = MAX( SO(K,6), SO(MAX(K-1,1),6) )  !stable profile
          SO(K,7) = MAX( SO(K,7), SO(MAX(K-1,1),7) )  !stable profile
*         write(6,'(a,i3,3f12.4)') 'k,t=',k,SO(K,3),SO(K,3),SI(K,3)
*         write(6,'(a,i3,3f12.4)') 'k,s=',k,SO(K,4),SI(K,4),SO(K,4)
*         write(6,'(a,i3,3f12.4)') 'k,r=',k,SO(K,5),SO(K,6),SO(K,7)
        ENDIF
      ENDDO
      CLOSE(12)
C
C     READ DEPTHS, TO GET NZ
C
      DO K= 1,99999
        READ(13,'(a)',IOSTAT=IOS) CLINE
        IF     (IOS.NE.0) THEN
          EXIT
        ENDIF
      ENDDO
      NZ = K-1
C
C     RE-READ THE DEPTHS
C
      ALLOCATE( ZZ(NZ), SSSHA(NZ), RA(NZ) )
C
      REWIND(13)
      DO K= 1,NZ
        READ(13,*) ZZ(K)
      ENDDO
C
C     CALCULATE THE FULL SSSH ANOMALY
C
      CALL STERIC(SI(1,5),PI,
     &            SO(1,5),PO, KDM,
     &              ZREF, ZZ,  NZ, AREF,SSSHA,RA)
C
C     OUTPUT THE FULL SSSH ANOMALY
C
      IF     (ZREF.GT.0.0) THEN
        WRITE(21,'(a,10x,f8.4,a)')
     &        '#  input ssha = ',SSHO-SSHI,' (m)'
        WRITE(21,'(a,f10.4,f8.4)')
     &        '#  zref,sssha = ',zref,aref
        WRITE(21,'(a)')
     &        '#        z   sssha  a.z/a.zref    rhoa'
        DO K= 1,NZ
          WRITE(21,'(f10.4,f8.4,f12.5,f8.4)')
     &      ZZ(K),SSSHA(K),SSSHA(K)/AREF,RA(K)
        ENDDO !k
      ELSE
        WRITE(21,'(a,10x,f8.4,a)')
     &        '#  ssha = ',SSHO-SSHI,' (m)'
        WRITE(21,'(a)')
     &        '#        z   sssha    rhoa'
        DO K= 1,NZ
          WRITE(21,'(f10.4,2f8.4)')
     &      ZZ(K),SSSHA(K),RA(K)
        ENDDO !k
      ENDIF
      CLOSE(21)
C
      IF     (CFILET.NE.' ') THEN
C
C       CALCULATE THE T-ONLY SSSH ANOMALY
C
        CALL STERIC(SI(1,5),PI,
     &              SO(1,6),PO, KDM,
     &                ZREF, ZZ,  NZ, AREF,SSSHA,RA)
C
C       OUTPUT THE T-ONLY SSSH ANOMALY
C
        IF     (ZREF.GT.0.0) THEN
          WRITE(22,'(a,10x,f8.4,a)')
     &          '#  input ssha = ',SSHO-SSHI,' (m)'
          WRITE(22,'(a,f10.4,f8.4)')
     &          '#  zref,sssha = ',zref,aref
          WRITE(22,'(a)')
     &          '#        z   sssha  a.z/a.zref    rhoa  (T-only)'
          DO K= 1,NZ
            WRITE(22,'(f10.4,f8.4,f12.5,f8.4)')
     &        ZZ(K),SSSHA(K),SSSHA(K)/AREF,RA(K)
          ENDDO !k
        ELSE
          WRITE(22,'(a,10x,f8.4,a)')
     &          '#  ssha = ',SSHO-SSHI,' (m)'
          WRITE(22,'(a)')
     &          '#        z   sssha    rhoa'
          DO K= 1,NZ
            WRITE(22,'(f10.4,2f8.4)')
     &        ZZ(K),SSSHA(K),RA(K)
          ENDDO !k
        ENDIF
        CLOSE(22)
C
C       CALCULATE THE S-ONLY SSSH ANOMALY
C
        CALL STERIC(SI(1,5),PI,
     &              SO(1,7),PO, KDM,
     &                ZREF, ZZ,  NZ, AREF,SSSHA,RA)
C
C       OUTPUT THE S-ONLY SSSH ANOMALY
C
        IF     (ZREF.GT.0.0) THEN
          WRITE(23,'(a,10x,f8.4,a)')
     &          '#  input ssha = ',SSHO-SSHI,' (m)'
          WRITE(23,'(a,f10.4,f8.4)')
     &          '#  zref,sssha = ',zref,aref
          WRITE(23,'(a)')
     &          '#        z   sssha  a.z/a.zref    rhoa  (S-only)'
          DO K= 1,NZ
            WRITE(23,'(f10.4,f8.4,f12.5,f8.4)')
     &        ZZ(K),SSSHA(K),SSSHA(K)/AREF,RA(K)
          ENDDO !k
        ELSE
          WRITE(23,'(a,10x,f8.4,a)')
     &          '#  ssha = ',SSHO-SSHI,' (m)'
          WRITE(23,'(a)')
     &          '#        z   sssha    rhoa'
          DO K= 1,NZ
            WRITE(23,'(f10.4,2f8.4)')
     &        ZZ(K),SSSHA(K),RA(K)
          ENDDO !k
        ENDIF
        CLOSE(23)
      ENDIF
      END
      subroutine steric(ri,pi, ro,po,  kk,
     &                  zref,zz, nz, aref,sssha,ra)
      implicit none
c
      integer kk,nz
      real    zref,aref
      real    ri(kk),pi(kk+1),
     &        ro(kk),po(kk+1),
     &        zz(nz),sssha(nz),ra(nz)
c
c**********
c*
c  1) calculate the steric SSH anomaly between two density profiles.
c
c  2) input arguments:
c       ri    - 1st density profile values
c       pi    - 1st density profile layer depths
c       ro    - 2nd density profile values
c       po    - 2nd density profile layer depths
c       kk    - number of layers
c       zref  - reference depth
c       zz    - sample depths
c       nz    - number of sample depths
c
c  3) output arguments:
c       aref  - steric SSH anomaly at zref
c       sssha - steric SSH anomaly at zz(:)
c
c  4) Alan J. Wallcraft,  Naval Research Laboratory,  Dec. 2015.
c*
c**********
c
      integer k,kz
      real    r0,riz,roz,q,z0,zm,zp
      real*8  risum(0:kk),rosum(0:kk),risumz,rosumz
c
      r0 = 0.5*(ri(1)+ri(kk))
      
      risum(0) = 0.d0
      rosum(0) = 0.d0
      do k= 1,kk
c ---   density * thickness at layer depths
        risum(k) = risum(k-1) + (ri(k)-r0)*(pi(k+1)-pi(k))
        rosum(k) = rosum(k-1) + (ro(k)-r0)*(po(k+1)-po(k))
      enddo !k
c
      do kz= 1,nz
c
c ---   density * thickness at zz depths
c
        if     (zz(kz).eq.0.0) then
          risumz = 0.d0
             riz = ri(1)
        elseif (zz(kz).gt.pi(kk+1)) then
          risumz = risum(kk)
             riz = ri(kk)
        else
          do k= 1,kk
            if     (pi(k+1).ge.zz(kz)) then  !pi(k).le.zz(kz)
              exit
            endif
          enddo
          risumz = risum(k-1) + (ri(k)-r0)*(zz(kz)-pi(k))
          z0 = 0.5*(pi(k+1)+pi(k))
          if     (z0.ge.zz(kz)) then
c           zz(k) is in the upper half of the layer
            if     (k.eq.1) then
              riz = ri(k)
            else
              zm  = 0.5*(pi(k)+pi(k-1))
              q   = (z0 - zz(kz))/(z0 - zm)
              riz = q*ri(k-1) + (1.0-q)*ri(k)
            endif
          else
c           zz(k) is in the upper half of the layer
            if     (k.eq.kk) then
              riz = ri(k)
            else
              zp  = 0.5*(pi(k+2)+pi(k+1))
              q   = (zz(kz) - z0)/(zp - z0)
              riz = q*ri(k+1) + (1.0-q)*ri(k)
            endif
          endif !top:bottom half
        endif
c
        if     (zz(kz).eq.0.0) then
          rosumz = 0.d0
             roz = ro(1)
        elseif (zz(kz).gt.po(kk+1)) then
          rosumz = rosum(kk)
             roz =    ro(kk)
        else
          do k= 1,kk
            if     (po(k+1).ge.zz(kz)) then  !po(k).le.zz(kz)
              exit
            endif
          enddo
          rosumz = rosum(k-1) + (ro(k)-r0)*(zz(kz)-po(k))
          z0 = 0.5*(po(k+1)+po(k))
          if     (z0.ge.zz(kz)) then
c           zz(k) is in the upper half of the layer
            if     (k.eq.1) then
              roz = ro(k)
            else
              zm  = 0.5*(po(k)+po(k-1))
              q   = (z0 - zz(kz))/(z0 - zm)
              roz = q*ro(k-1) + (1.0-q)*ro(k)
            endif
          else
c           zz(k) is in the upper half of the layer
            if     (k.eq.kk) then
              roz = ro(k)
            else
              zp  = 0.5*(po(k+2)+po(k+1))
              q   = (zz(kz) - z0)/(zp - z0)
              roz = q*ro(k+1) + (1.0-q)*ro(k)
            endif
          endif !top:bottom half
        endif
c
c ---   sssha at zz depths
c
           ra(kz) =     riz -    roz
        sssha(kz) = (risumz - rosumz)/
     &     (1000.d0 + r0 + rosumz/min(zz(kz),po(kk+1)))
      enddo !kz
c
      if     (zref.ne.0.0) then
c
c ---   density * thickness at zref depth
c
        if     (zref.gt.pi(kk+1)) then
          risumz = risum(kk)
        else
          do k= 1,kk
            if     (pi(k+1).ge.zref) then  !pi(k).le.zref
              exit
            endif
          enddo
          risumz = risum(k-1) + (ri(k)-r0)*(zref-pi(k))
        endif
c
        if     (zref.gt.po(kk+1)) then
          rosumz = rosum(kk)
        else
          do k= 1,kk
            if     (po(k+1).ge.zref) then  !po(k).le.zref
              exit
            endif
          enddo
          rosumz = rosum(k-1) + (ro(k)-r0)*(zref-po(k))
        endif
c
c ---   sssha at zz depths
c
        aref = (risumz - rosumz)/
     &     (1000.d0 + r0 + rosumz/min(zref,po(kk+1)))
      endif !zref
      return
      end subroutine steric
