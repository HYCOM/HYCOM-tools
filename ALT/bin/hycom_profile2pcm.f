      PROGRAM HYCOM_PROFILE2PCM
      IMPLICIT NONE
C
C  hycom_profile2pcm - Usage: hycom_profile2pcm archv.txt archp.txt [[nK] nT]
C
C                 converts a HYCOM text profile file to one
C                 which will plot as the PCM profile
C
C   archv.txt is assumed to be an HYCOM archive text profile file
C   archp.txt will be the output, PCM plot, text profile file
C   nK is the number of K-profiles (visc or visc,tdif or visc,tdif,sdif)
C   nK is the number of tracers
C
C  the input and the output are equivalent HYCOM profiles, but
C  the ouput adds two zero thickness layers at each interface
C  and (if K-profiles are output) at each cell center.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  March 2002
C
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      CHARACTER*240 CFILEA,CFILEC,CFORMAT
      CHARACTER*240 CLINE
      REAL          THK,DEPTH,FLAG,DIFF(3)
      INTEGER       IOS,K,KDM,KI,KK,KP,NDIF,KT,NTRC
C
      REAL, ALLOCATABLE :: SI(:,:),SK(:,:),P(:)
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.2) THEN
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CFILEC)
        NDIF = 0
        NTRC = 0
      ELSEIF (NARG.EQ.3) THEN
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CFILEC)
        CALL GETARG(3,CARG)
        READ(CARG,*) NDIF
        NTRC = 0
      ELSEIF (NARG.EQ.4) THEN
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CFILEC)
        CALL GETARG(3,CARG)
        READ(CARG,*) NDIF
        CALL GETARG(4,CARG)
        READ(CARG,*) NTRC
      ELSE
        WRITE(6,*)
     +    'Usage: hycom_profile2pcm archv.txt archz.txt [[nK] nT]'
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
      ALLOCATE( P(KDM+1), SI(KDM,5+NTRC), SK(KDM,3) )
      SK(:,:) = 0.0
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
        READ(CLINE,*) KI,(SI(K,KK),KK=1,5     ),THK,DEPTH,
     &                   (SK(K,KP),KP=1,  NDIF),
     &                   (SI(K,KT),KT=6,5+NTRC)
        P(K+1) = P(K) + THK
        IF     (THK.EQ.0.0) THEN
          DO KK= 1,5
            SI(K,KK)=SI(K-1,KK)
          ENDDO !kk
          DO KT= 6,5+NTRC
            SI(K,KT)=SI(K-1,KT)
          ENDDO !kk
        ENDIF
      ENDDO
      CLOSE(11)
C
C     OUTPUT, 3 OR 6 COPIES OF EACH LAYER.
C
      IF     (NTRC.EQ.0) THEN
        WRITE(CFORMAT,'(a)')
     &    '(3a)'
      ELSE
        WRITE(CFORMAT,'(a,i2,a)')
     &    '(3a,',NTRC,'a)'
      ENDIF
      IF     (NDIF.EQ.0) THEN
        WRITE(21,CFORMAT)
     &      '#  k',
     &      '    utot    vtot  p.temp    saln  p.dens',
     &      '    thkns      dpth',
     &      ('  tracer',KT=1,NTRC)
      ELSEIF (NDIF.EQ.1) THEN
        WRITE(21,CFORMAT)
     &      '#  k',
     &      '    utot    vtot  p.temp    saln  p.dens',
     &      '    thkns      dpth  viscty',
     &      ('  tracer',KT=1,NTRC)
      ELSEIF (NDIF.EQ.2) THEN
        WRITE(21,CFORMAT)
     &      '#  k',
     &      '    utot    vtot  p.temp    saln  p.dens',
     &      '    thkns      dpth  viscty  t-diff',
     &      ('  tracer',KT=1,NTRC)
      ELSEIF (NDIF.EQ.3) THEN
        WRITE(21,CFORMAT)
     &      '#  k',
     &      '    utot    vtot  p.temp    saln  p.dens',
     &      '    thkns      dpth  viscty  t-diff  s-diff',
     &      ('  tracer',KT=1,NTRC)
      ENDIF
C
      IF     (NTRC.EQ.0) THEN
        IF     (NDIF.EQ.0) THEN
          WRITE(CFORMAT,'(a)')
     &    '(i4,2f8.2,3f8.3,f9.3,f10.3)'
        ELSE
          WRITE(CFORMAT,'(a,i1,a)')
     &    '(i4,2f8.2,3f8.3,f9.3,f10.3,',NDIF,'f8.2)'
        ENDIF !ndifo
      ELSE
        IF     (NDIF.EQ.0) THEN
          WRITE(CFORMAT,'(a,i2,a)')
     &    '(i4,2f8.2,3f8.3,f9.3,f10.3,',             NTRC,'f8.4)'
        ELSE
          WRITE(CFORMAT,'(a,i1,a,i2,a)')
     &    '(i4,2f8.2,3f8.3,f9.3,f10.3,',NDIF,'f8.2,',NTRC,'f8.4)'
        ENDIF !ndifo
      ENDIF !ntrc
C
      IF     (NDIF.EQ.0) THEN
C
C       ZERO THICKNESS LAYERS AT INTERFACES
C
        DO K= 1,KDM
          THK = P(K+1) - P(K)
          WRITE(21,CFORMAT)
     &      3*K-2,(SI(K,KK),KK=1,5     ),0.0,     P(K),
     &            (SI(K,KT),KT=6,5+NTRC)
          WRITE(21,CFORMAT)
     &      3*K-1,(SI(K,KK),KK=1,5     ),THK,0.5*(P(K)+P(K+1)),
     &            (SI(K,KT),KT=6,5+NTRC)
          WRITE(21,CFORMAT)
     &      3*K  ,(SI(K,KK),KK=1,5     ),0.0,          P(K+1),
     &            (SI(K,KT),KT=6,5+NTRC)
        ENDDO !k
      ELSE
C
C       ZERO THICKNESS LAYERS AT INTERFACES AND THE LAYER CENTER
C
        DIFF(:) = SK(1,:)
        DO K= 1,KDM
          THK = 0.5*(P(K+1) - P(K))  !two half layers
          WRITE(21,CFORMAT)
     &      6*K-5,(SI(K,KK),KK=1,5     ),0.0,    P(K),
     &            (DIFF(KP),KP=1,  NDIF),  !SK(K-1,:)
     &            (SI(K,KT),KT=6,5+NTRC)
          WRITE(21,CFORMAT)
     &      6*K-4,(SI(K,KK),KK=1,5     ),THK,0.75*P(K)+0.25*P(K+1),
     &            (DIFF(KP),KP=1,  NDIF),  !SK(K-1,:)
     &            (SI(K,KT),KT=6,5+NTRC)
          WRITE(21,CFORMAT)
     &      6*K-3,(SI(K,KK),KK=1,5     ),0.0,0.50*P(K)+0.50*P(K+1),
     &            (DIFF(KP),KP=1,  NDIF),  !SK(K-1,:)
     &            (SI(K,KT),KT=6,5+NTRC)
          DIFF(:) = SK(K,:)
          WRITE(21,CFORMAT)
     &      6*K-2,(SI(K,KK),KK=1,5     ),0.0,0.50*P(K)+0.50*P(K+1),
     &            (DIFF(KP),KP=1,  NDIF),  !SK(K,:)
     &            (SI(K,KT),KT=6,5+NTRC)
          WRITE(21,CFORMAT)
     &      6*K-1,(SI(K,KK),KK=1,5     ),THK,0.25*P(K)+0.75*P(K+1),
     &            (DIFF(KP),KP=1,  NDIF),  !SK(K,:)
     &            (SI(K,KT),KT=6,5+NTRC)
          WRITE(21,CFORMAT)
     &      6*K  ,(SI(K,KK),KK=1,5     ),0.0,               P(K+1),
     &            (DIFF(KP),KP=1,  NDIF),  !SK(K,:)
     &            (SI(K,KT),KT=6,5+NTRC)
        ENDDO !k
      ENDIF !ndif==0:else
      CLOSE(21)
      END
