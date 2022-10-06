      PROGRAM HYCOM_PROFILE_LOCSIG
      IMPLICIT NONE
C
C  hycom_profile_locsig_25t - Usage:  hycom_profile_locsig_25t archv.txt archs.txt [mean]
C
C                 add locally referenced potential density (sigloc)
C                 as two tracers to a HYCOM text profile file
C                 if mean is present, add a vertical mean at the end
C                 version for 25-term equation of state
C
C   archv.txt is assumed to be an HYCOM archive text profile file
C   archs.txt will be the output text profile file, with sigloc and
C   sigloc from dsiglocdt and dsiglocds as tracers
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  December 2006.
C
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      CHARACTER*240 CFILEA,CFILEC,CFORMAT
      CHARACTER*240 CLINE
      LOGICAL       LMEAN
      REAL          THK,DEPTH,FLAG,alfadt,betads,ROFF
      REAL*8        THKMN,SM(7),KEMN
      INTEGER       IOS,K,KDM,KI,KK,KP
C
c------------------------------------------------------------------------
      include '../include/stmt_fns_SIGMA0_17term.h'
c------------------------------------------------------------------------
      REAL, ALLOCATABLE :: SI(:,:),P(:)
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      LMEAN = NARG.EQ.3
      IF     (NARG.EQ.2 .OR. LMEAN) THEN
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CFILEC)
      ELSE
        WRITE(6,*)
     +    'Usage:  hycom_profile_locsig_25t archv.txt archs.txt [mean]'
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
      DO K= 1,99999
        READ(11,'(a)',IOSTAT=IOS) CLINE
        IF     (IOS.NE.0) THEN
          EXIT
        ENDIF
      ENDDO
      KDM = K-1
C
C     RE-READ THE ISOPYCNAL PROFILE.
C
      ALLOCATE( P(KDM+1), SI(KDM,7) )
C
      REWIND(11)
      DO K= 1,99
        READ( 11,'(a)') CLINE
        IF     (CLINE(1:5).EQ.'#  k ') then
          EXIT
        ENDIF
      ENDDO
      P(1) =  0.0
      DO K= 1,KDM
        READ(11,'(a)',IOSTAT=IOS) CLINE
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'Error: inconsistent input profile'
          CALL EXIT(6)
        ENDIF
        READ(CLINE,*) KI,(SI(K,KK),KK=1,5),THK,DEPTH
        P(K+1) = P(K) + THK
        IF     (THK.EQ.0.0) THEN
          DO KK= 1,5
            SI(K,KK)=SI(K-1,KK)
          ENDDO !kk
        ENDIF
      ENDDO
      CLOSE(11)
C
C     SIGLOC, AS TRACERS
C
      K=1
        SI(K,6) = sigloc(r8(SI(1,3)),
     &                   r8(SI(1,4)),
     &                   r8(9806.0*0.5*(P(1)+P(2))))
        SI(K,7) = SI(1,5)
      DO K= 2,KDM
        SI(K,6) = sigloc(r8(SI(K,3)),
     &                   r8(SI(K,4)),
     &                   r8(9806.0*0.5*(P(K)+P(K+1))))
        alfadt = dsiglocdt(r8(0.5*(SI(K-1,3)+SI(K,3))),
     &                     r8(0.5*(SI(K-1,4)+SI(K,4))),
     &                     r8(9806.0*P(K)))*
     &                         (SI(K-1,3)-SI(K,3))
        betads = dsiglocds(r8(0.5*(SI(K-1,3)+SI(K,3))),
     &                     r8(0.5*(SI(K-1,4)+SI(K,4))),
     &                     r8(9806.0*P(K)))*
     &                         (SI(K-1,4)-SI(K,4))
        SI(K,7) = SI(K-1,7)-alfadt-betads
      ENDDO
C
C     OUTPUT
C
        IF     (.NOT. LMEAN) THEN
          WRITE(CFORMAT,'(a)')
     &      '(3a)'
          WRITE(21,CFORMAT)
     &      '#  k',
     &      '    utot    vtot  p.temp    saln  p.dens',
     &      '    thkns      dpth  tracer  tracer'
C
          WRITE(CFORMAT,'(a)')
     &      '(i4,2f8.2,3f8.4,f9.3,f10.3,2f8.4)'
        ELSE
          WRITE(CFORMAT,'(a)')
     &      '(3a)'
          WRITE(21,CFORMAT)
     &      '#  k',
     &      '    utot    vtot  p.temp    saln  p.dens',
     &      '    thkns      dpth  tracer  tracer            k.e.'
C
          WRITE(CFORMAT,'(a)')
     &      '(i4,2f8.2,3f8.4,f9.3,f10.3,2f8.4,,f16.4)'
        ENDIF
C
        IF     (LMEAN) THEN
          SM(1:7) = 0.d0
          THKMN   = 0.d0
           KEMN   = 0.d0
        ENDIF
        DO K= 1,KDM
          THK = P(K+1) - P(K)
          IF     (.NOT. LMEAN) THEN
            WRITE(21,CFORMAT)
     &        K,(SI(K,KK),KK=1,5),THK,0.5*(P(K)+P(K+1)),
     &          (SI(K,KK),KK=6,7)
          ELSE !lmean
            SM(1:7) = SM(1:7) + THK*SI(K,1:7)
            THKMN   = THKMN   + THK
            KEMN    =  KEMN   + THK*(1000.0+SI(K,5))*         !density
     &                          0.5*(SI(K,1)**2 + SI(K,2)**2) !u^2+v^2
            WRITE(21,CFORMAT)
     &        K,(SI(K,KK),KK=1,5),THK,0.5*(P(K)+P(K+1)),
     &          (SI(K,KK),KK=6,7),KEMN
          ENDIF
        ENDDO !k
        IF     (LMEAN) THEN
          WRITE(CFORMAT,'(a)')
     &      '(a4,2f8.2,3f8.4,f9.3,f10.3,2f8.4,f16.4)'
          WRITE(21,CFORMAT)
     &      '#SUM',(SM(KK)/THKMN,KK=1,5),THKMN,0.5*THKMN,
     &             (SM(KK)/THKMN,KK=6,7),KEMN
        ENDIF
      CLOSE(21)
      END
