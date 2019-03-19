      PROGRAM PROFILE
      IMPLICIT NONE
C
C  hycom_profile_stokes - Usage:  hycom_profile_stokes stokes.txt archz.txt
C
C                 converts an Stokes drift profile to a hycom profile file.
C
C   stokes.txt is assumed to be a Stokes drift text profile file
C   archz.txt will be the output hycom text profile file
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  Oct. 2008
C
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      CHARACTER*240 CFILEA,CFILEC
      REAL          THK,BOT,PZK,PZKP,ZZK,UB,US,VB,VS
      INTEGER       IOS,ITYPE,K,KS,KZ
      INTEGER       I
C
      REAL, ALLOCATABLE :: DZ(:),PZ(:),UZ(:),VZ(:)
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.2) THEN
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CFILEC)
      ELSE
        WRITE(6,*)
     +  'Usage: hycom_profile_stokes stokes.txt archz.txt'
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
C     READ THE STOKES PROFILE HEADER, TO GET KZ.
C
      DO 
        READ(11,*,IOSTAT=IOS) KZ
        IF     (IOS.NE.0) THEN
          write(6,*) 'KZ = ',KZ
          EXIT
        ENDIF
      ENDDO
C
C     INPUT ALL LEVELS
C
      ALLOCATE( DZ(KZ), PZ(KZ), UZ(KZ), VZ(KZ) )
C
      REWIND(11)
      DO K= 1,KZ
        READ(11,*,IOSTAT=IOS) KS,UZ(K),VZ(K),DZ(K),PZ(K)
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'Error: bad Stokes file (level',K,')'
          CALL EXIT(6)
        ENDIF
        IF     (K.EQ.2) THEN
          UB = UZ(2)*DZ(2)
          VB = VZ(2)*DZ(2)
        ELSEIF (K.GT.2) THEN
          UB = UB+ UZ(K)*DZ(K)
          VB = VB+ VZ(K)*DZ(K)
        ENDIF
      ENDDO !k
      CLOSE(11)
      UB = UB/PZ(KZ)
      VB = VB/PZ(KZ)
C
C     LEVEL 1 IS THE SURFACE
C
      US = UZ(1)
      VS = VZ(1)
      DZ(1) = DZ(2)*0.01
      PZ(1) = DZ(1)
      DZ(2) = DZ(2) - DZ(1)
      UZ(2) = ((DZ(1) + DZ(2))*UZ(2) - DZ(1)*UZ(1))/DZ(2)
C
C     OUTPUT.
C
      WRITE(21,'(2a/a,4i7,2f7.1,i7,a)')
     &  '##   expt    idm    jdm    kdm    lon    lat',
     &  ' yrflag   DTG',
     &  '##',    0,     1,     1,    KZ,   0.0,   0.0,
     &  3," 1900/12/31"
      WRITE(21,'(3a/a,f11.2,f8.2,f8.1,2f9.3,3f8.3,4f8.2)')
     &  '## model-day  srfhgt  surflx',
     &  '     dpbl   dpmixl    tmix    smix   thmix',
     &  '    usur    vsur   ubavg   vbavg',
     &  '#',
     &  0.0,   ! model day
     &  0.0,   ! cm
     &  0.0,   ! W/m**2
     &  0.0,   ! m
     &  0.0,   ! m
     &  0.0,   ! degC
     &  0.0,   ! psu
     &  0.0,   ! SigmaT
     &  US*100.0, !cm/s
     &  VS*100.0, ! cm/s
     &  UB*100.0, ! cm/s
     &  VB*100.0  ! cm/s
      WRITE(21,'(a,a)')
     &  '#   k',
     &  '    utot    vtot    temp    saln    dens    thkns      dpth'
      PZKP = 0.0
      DO K= 1,KZ
        PZK  = PZKP
        PZKP = PZ(K)
        WRITE(21,'(i4,1x,2f8.2,3f8.3,f9.3,f10.3)')
     &    K,
     &    UZ(K)*100.0,       !cm/s
     &    VZ(K)*100.0,       !cm/s
     &    0.0,               !degC
     &    0.0,               !psu
     &    0.0,               !SigmaT
     &    PZKP-PZK,          !m
     &    0.5*(PZKP+PZK)     !m
      ENDDO
      CLOSE(21)
      END
