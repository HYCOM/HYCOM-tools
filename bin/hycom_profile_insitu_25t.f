      PROGRAM HYCOM_PROFILE_INSITU_25T
C
      IMPLICIT NONE
C
C  hycom_profile_insitu_25t - Usage:  hycom_profile_insitu_25t archv.txt archs.txt [mean]
C
C                 add in-situ density (locsig)
C                 as a tracer to a HYCOM text profile file
C                 if mean is present, add a vertical mean at the end
C                 version for HYCOM 25-term rational function equation of state
C
C   archv.txt is assumed to be an HYCOM archive text profile file
C   archs.txt will be the output text profile file, with in-situ density
C   as a tracer
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  August 2018.
C
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      CHARACTER*240 CFILEA,CFILEC,CFORMAT
      CHARACTER*240 CLINE
      LOGICAL       LMEAN
      REAL          THK,DEPTH,DUMMY,SSH
      REAL*8        THKMN,SM(6),KEMN
      REAL*8        T5(5),S5(5),P5(5),R5(5),SCLP
      INTEGER       IOS,K,KDM,KI,KK,KP
C
c------------------------------------------------------------------------
      include '../include/stmt_fns_SIGMA2_17term.h'
c------------------------------------------------------------------------
c
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
     +    'Usage:  '//
     +    'hycom_profile_insitu_25t archv.txt archs.txt [mean]'
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
C
        IF     (K.EQ.4) THEN
          READ(CLINE(2:),*) DUMMY,SSH
          SSH = 0.01*SSH  !cm to m
*               write(6,'(a,f8.5)') 'ssh:',SSH
        ENDIF
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
      ALLOCATE( P(KDM+1), SI(KDM,6) )
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
C     IN-SITU DENISTY, AS A TRACER
C
      DO K= 1,KDM
        t5(:) = r8(SI(K,3))
        s5(:) = r8(SI(K,4))
        p5(1) = r8(9806.0*(1.00*P(K)+0.00*P(K+1)))  !9806.0 for m to pressure
        p5(2) = r8(9806.0*(0.75*P(K)+0.25*P(K+1)))
        p5(3) = r8(9806.0*(0.50*P(K)+0.50*P(K+1)))
        p5(4) = r8(9806.0*(0.25*P(K)+0.75*P(K+1)))
        p5(5) = r8(9806.0*(0.00*P(K)+1.00*P(K+1)))
*             write(6,'(a,i3,5f8.2)') 'p5:',k,p5(:)*1.d-4
        r5(1) = sigloc(t5(1),s5(1),p5(1))
        r5(2) = sigloc(t5(2),s5(2),p5(2))
        r5(3) = sigloc(t5(3),s5(3),p5(3))
        r5(4) = sigloc(t5(4),s5(4),p5(4))
        r5(5) = sigloc(t5(5),s5(5),p5(5))
c ---   Bode's (Boole's) Rule for integration
*             write(6,'(a,i3,5f8.2)') 'r5:',k,r5(:)
        SI(K,6) = 1.0d0/90.0d0*( 7.0d0*(r5(1)+r5(5))+
     &                          32.0d0*(r5(2)+r5(4))+
     &                          12.0d0* r5(3)        )
      ENDDO !k
C
C     OUTPUT
C
        IF     (.NOT. LMEAN) THEN
          WRITE(CFORMAT,'(a)')
     &      '(3a)'
          WRITE(21,CFORMAT)
     &      '#  k',
     &      '    utot    vtot  p.temp    saln  p.dens',
     &      '    thkns      dpth  tracer'
C
          WRITE(CFORMAT,'(a)')
     &      '(i4,2f8.2,3f8.4,f9.3,f10.3,f8.4)'
        ELSE
          WRITE(CFORMAT,'(a)')
     &      '(3a)'
          WRITE(21,CFORMAT)
     &      '#  k',
     &      '    utot    vtot  p.temp    saln  p.dens',
     &      '    thkns      dpth  tracer            k.e.'
C
          WRITE(CFORMAT,'(a)')
     &      '(i4,2f8.2,3f8.4,f9.3,f10.3,f8.4,f16.4)'
        ENDIF
C
        IF     (LMEAN) THEN
          SM(1:6) = 0.d0
          THKMN   = 0.d0
           KEMN   = 0.d0
        ENDIF
        DO K= 1,KDM
          THK = P(K+1) - P(K)
          IF     (.NOT. LMEAN) THEN
            WRITE(21,CFORMAT)
     &        K,(SI(K,KK),KK=1,5),THK,0.5*(P(K)+P(K+1)),
     &          (SI(K,KK),KK=6,6)
          ELSE !lmean
            SM(1:6) = SM(1:6) + THK*SI(K,1:6)
            THKMN   = THKMN   + THK
            KEMN    =  KEMN   + THK*(1000.0+SI(K,5))*         !density
     &                          0.5*(SI(K,1)**2 + SI(K,2)**2) !u^2+v^2
            WRITE(21,CFORMAT)
     &        K,(SI(K,KK),KK=1,5),THK,0.5*(P(K)+P(K+1)),
     &          (SI(K,KK),KK=6,6),KEMN
          ENDIF
        ENDDO !k
        IF     (LMEAN) THEN
          WRITE(CFORMAT,'(a)')
     &      '(a4,2f8.2,3f8.4,f9.3,f10.3,f8.4,f16.4)'
          WRITE(21,CFORMAT)
     &      '#SUM',(SM(KK)/THKMN,KK=1,5),THKMN,0.5*THKMN,
     &             (SM(KK)/THKMN,KK=6,6),KEMN
        ENDIF
      CLOSE(21)
      END
