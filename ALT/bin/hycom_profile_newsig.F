      PROGRAM HYCOM_PROFILE_NEWSIG
      IMPLICIT NONE
C
C  hycom_profile_newsig - Usage:  hycom_profile_newsig archv.txt sigver archo.txt [stable]
C
C                 convert potenital density to sigver
C
C   archv.txt  is an HYCOM archive text profile file
C   sigver     equation of state type
C   archo.txt  will be the output text profile file, with sigver applied
C
C   if stable is present, the delta density is included as a tracer
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  December 2015.
C
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      CHARACTER*240 CFILEA,CFILEO,CFILEC,CFORMAT
      CHARACTER*240 CLINE
      LOGICAL       LSTABLE
      REAL          THK,DEPTH,FLAG,ROFF
      INTEGER       IOS,K,KDM,KI,KK,KP,SIGVER
C
      REAL, ALLOCATABLE :: SI(:,:),P(:)
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      LSTABLE = NARG.EQ.4
C
      IF     (NARG.EQ.3 .OR. NARG.EQ.4) THEN
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CLINE)
        READ(CLINE,*) SIGVER
        CALL GETARG(3,CFILEC)
      ELSE
        WRITE(6,*)
     +    'Usage:  hycom_profile_newsig archv.txt sigver archo.txt'
     +    // ' [stable]'
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
      ALLOCATE( P(KDM+1), SI(KDM,5) )
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
      ENDDO !k
C
      CALL SIG_I(SIGVER)
      DO K= 1,KDM
        CALL SIG_P(SI(K,3),SI(K,4), SI(K,5))
      ENDDO !k
C
C     OUTPUT
C
      IF     (.NOT.LSTABLE) THEN  !default
        WRITE(CFORMAT,'(a)')
     &    '(3a)'
        WRITE(21,CFORMAT)
     &      '#  k',
     &      '    utot    vtot  p.temp    saln  p.dens',
     &      '    thkns      dpth'
C
          WRITE(CFORMAT,'(a)')
     &      '(i4,2f8.2,3f8.4,f9.3,f10.3)'
C
        DO K= 1,KDM
          THK = P(K+1) - P(K)
          WRITE(21,CFORMAT)
     &      K,(SI(K,KK),KK=1,5),THK,0.5*(P(K)+P(K+1))
        ENDDO !k
      ELSE !lstable
        WRITE(CFORMAT,'(a)')
     &    '(3a)'
        WRITE(21,CFORMAT)
     &      '#  k',
     &      '    utot    vtot  p.temp    saln  p.dens',
     &      '    thkns      dpth  tracer'
C
          WRITE(CFORMAT,'(a)')
     &      '(i4,2f8.2,3f8.4,f9.3,f10.3,f8.4)'
C
        DO K= 1,KDM
          THK = P(K+1) - P(K)
          WRITE(21,CFORMAT)
     &      K,(SI(K,KK),KK=1,5),THK,0.5*(P(K)+P(K+1)),
     &      SI(MIN(K+1,KDM),5)-SI(K,5)
        ENDDO !k
      ENDIF !.not.lstable:else
      CLOSE(21)
      END
