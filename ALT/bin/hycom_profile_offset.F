      PROGRAM HYCOM_PROFILE_LOCSIG
      IMPLICIT NONE
C
C  hycom_profile_offset - Usage:  hycom_profile_offset archv.txt offset.txt archo.txt
C
C                 add u,v,t,s offsets to a HYCOM text profile file
C
C   archv.txt  is an HYCOM archive text profile file
C   offset.txt is a set of u,v,t,s offsets for all input layers
C   archo.txt  will be the output text profile file, with offsets applied
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
      REAL          THK,DEPTH,FLAG,ROFF
      INTEGER       IOS,K,KDM,KI,KK,KP,SIGVER
C
      REAL, ALLOCATABLE :: SI(:,:),OI(:,:),P(:)
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.3) THEN
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CFILEO)
        CALL GETARG(3,CFILEC)
      ELSE
        WRITE(6,*)
     +    'Usage:  hycom_profile_offset archv.txt offset.txt archo.txt'
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
      OPEN(UNIT=12, FILE=CFILEO, FORM='FORMATTED', STATUS='OLD',
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
      ALLOCATE( P(KDM+1), SI(KDM,5), OI(KDM,4))
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
C
        READ(12,'(a)',IOSTAT=IOS) CLINE
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'Error: inconsistent input offset'
          CALL EXIT(6)
        ENDIF
        READ(CLINE,*) (OI(K,KK),KK=1,4)
      ENDDO
      CLOSE(11)
      IF     (SI(KDM,5).LT.30.0) THEN
        SIGVER = 1  ! 7-term sigma-0
      ELSE
        SIGVER = 6  !17-term sigma-2
      ENDIF
      CALL SIG_I(SIGVER)
C
C     ADD OFFSET
C
      K=1
      DO K= 1,KDM
        SI(K,1) = SI(K,1) + OI(K,1) !u
        SI(K,2) = SI(K,2) + OI(K,2) !v
        SI(K,3) = SI(K,3) + OI(K,3) !t
        SI(K,4) = SI(K,4) + OI(K,4) !s
        CALL SIG_P(SI(K,3),SI(K,4), SI(K,5))
      ENDDO
C
C     OUTPUT
C
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
      CLOSE(21)
      END
