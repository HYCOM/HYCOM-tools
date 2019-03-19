      PROGRAM PROFILE
      IMPLICIT NONE
C
C  hycom_profile - Usage:  hycom_profile archv.a ip jp
C
C                 generates a text profile file at point (ip,jp)
C                 from a HYCOM archive file.
C
C   archv.a is assumed to be an HYCOM archive data file, with companion
C   header file archv.b.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  August 2001.
C
      REAL*4     ONEM
      PARAMETER (ONEM=9806.0)
C
      REAL*4, ALLOCATABLE :: A(:,:)
      REAL*4              :: PAD(4096)
C
      INTEGER      IARGC
      INTEGER      NARG
      CHARACTER*240 CARG
C
      INTEGER      IDM,JDM,KDM,KTR,NSURF,NLAY,NDIF,NTRC,IEXPT,YRFLAG
      INTEGER      NPAD,IP,JP
      REAL         P(999),PROF(106),SURF(18),
     +             UB,VB,UK,VK,VISC,TDIF,THBASE,TIME
      CHARACTER*240 CFILEA,CFILEB,CFORMAT
C
      CHARACTER*18 CASN
      INTEGER      I,J,K,KREC,KREC0,IOS,NRECL
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.3) THEN
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CARG)
        READ(CARG,'(i4)') IP
        CALL GETARG(3,CARG)
        READ(CARG,'(i4)') JP
        THBASE = 25.0  !usual value for sigma-theta cases
      ELSEIF (NARG.EQ.4) THEN  ! undocumented, specify THBASE
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CARG)
        READ(CARG,'(i4)') IP
        CALL GETARG(3,CARG)
        READ(CARG,'(i4)') JP
        CALL GETARG(4,CARG)
        READ(CARG,*) THBASE
      ELSE
        WRITE(6,*) 'Usage: hycom_profile archv.a ip jp'
        CALL EXIT(1)
      ENDIF
C
C     EXTRACT MODEL PARAMETERS FROM ".b" FILE.
C
      CFILEB = CFILEA(1:LEN_TRIM(CFILEA)-1) // 'b'
      CALL READ_B(CFILEB,
     +            IEXPT,YRFLAG,IDM,JDM,KDM,NSURF,NLAY,NDIF,NTRC,
     +            THBASE,TIME)
C
C     OPEN ".a" FILE.
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( A(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_profile: could not allocate ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
C
      IF     (NPAD.EQ.0) THEN
        INQUIRE( IOLENGTH=NRECL) A
      ELSE
        INQUIRE( IOLENGTH=NRECL) A,PAD(1:NPAD)
      ENDIF
*     write(6,*) 'nrecl = ',nrecl
#ifdef CRAY
#ifdef t3e
      IF     (MOD(NRECL,4096).EQ.0) THEN
        WRITE(CASN,8000) NRECL/4096
 8000   FORMAT('-F cachea:',I4.4,':1:0')
        IU8 = 11
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          WRITE(6,*) 'Error: can''t asnunit 11'
          WRITE(6,*) 'ios  = ',ios8
          WRITE(6,*) 'casn = ',casn
          CALL EXIT(5)
        ENDIF
      ENDIF
#else
      CALL ASNUNIT(11,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t asnunit 11'
        WRITE(6,*) 'ios = ',ios
        CALL EXIT(5)
      ENDIF
#endif
#endif
      OPEN(UNIT=11, FILE=CFILEA, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',CFILEA(1:LEN_TRIM(CFILEA))
        WRITE(6,*) 'ios   = ',ios
        WRITE(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
C     SURFACE FIELDS.
C
      DO KREC= 1,NSURF
        READ(11,REC=KREC,IOSTAT=IOS) A
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'can''t read ',CFILEA(1:LEN_TRIM(CFILEA))
          CALL EXIT(4)
        ENDIF
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A(IP,JP),1)
#endif
*       write(6,*) 'REC,A.MIN = ',KREC,MINVAL(A),A(IP,JP)
        SURF(KREC) = A(IP,JP)
      ENDDO
      write(6,*) 'SURF = ',SURF(1:NSURF)
C
      IF     (NLAY.EQ.6) THEN
        UB = SURF(NSURF-1)
        VB = SURF(NSURF)
        UK = SURF(10)+UB
        VK = SURF(11)+VB
      ELSE !mean archive
        UB = SURF(NSURF-2)
        VB = SURF(NSURF-1)
        UK = SURF(10)  !already full mixed-layer velocity
        VK = SURF(11)  !already full mixed-layer velocity
      ENDIF
      UB = MAX( -9.99, MIN( 9.99, UB ) )
      VB = MAX( -9.99, MIN( 9.99, VB ) )
      UK = MAX( -9.99, MIN( 9.99, UK ) )
      VK = MAX( -9.99, MIN( 9.99, VK ) )
C
      K = LEN_TRIM(CFILEA)
      WRITE(6,'(2a/a,7i7,3a7)')
     &  '##   expt    idm    jdm    kdm  itest  jtest',
     &  ' yrflag   year    day     hr',
     &  '##',IEXPT,   IDM,   JDM,   KDM,    IP,    JP,
     &  YRFLAG,CFILEA(K-12:K-9),CFILEA(K-7:K-5),CFILEA(K-3:K-2)
      WRITE(6,'(3a/a,f11.2,f8.2,f8.1,2f9.3,3f8.3,4f8.2)')
     &  '## model-day  srfhgt  surflx',
     &  '     dpbl   dpmixl    tmix    smix   thmix    umix    vmix',
     &  '   ubavg   vbavg',
     &  '#',TIME,
     &  SURF( 2)/(9.806*0.01),  ! cm
     &  SURF( 3),               ! W/m**2
     &  SURF( 5)/ONEM,          ! m
     &  SURF( 6)/ONEM,          ! m
     &  SURF( 7),               ! degC
     &  SURF( 8),               ! psu
     &  SURF( 9)+THBASE,        ! SigmaT
     &  UK*100.0,               ! cm/s
     &  VK*100.0,               ! cm/s
     &  UB*100.0,               ! cm/s
     &  VB*100.0                ! cm/s
C
C     PER LAYER FIELDS.
C
      IF     (NTRC.EQ.0) THEN
        WRITE(CFORMAT,'(a)')
     &    '(3a)'
      ELSE
        WRITE(CFORMAT,'(a,i2,a)')
     &    '(3a,',NTRC,'a)'
      ENDIF
      WRITE(6,CFORMAT)
     &      '#  k',
     &      '    utot    vtot    temp    saln    dens',
     &      '    thkns      dpth  viscty  t-diff',
     &      ('  tracer',KTR=1,NTRC)
C
      IF     (NTRC.EQ.0) THEN
        WRITE(CFORMAT,'(a)')
     &    '(i4,2f8.2,3f8.3,f9.3,f10.3,2f8.2)'
      ELSE
        WRITE(CFORMAT,'(a,i2,a)')
     &    '(i4,2f8.2,3f8.3,f9.3,f10.3,2f8.2,',NTRC,'f8.4)'
      ENDIF
C
      P(1) = 0.0
      IF     (NLAY.EQ.6) THEN
        UB = SURF(NSURF-1)
        VB = SURF(NSURF)
      ELSE  !mean archive
        UB = 0.0  !already full layer velocity
        VB = 0.0  !already full layer velocity
      ENDIF
      VISC = 0.0
      TDIF = 0.0
      DO K= 1,KDM
        KREC0 = NSURF+(NLAY+NDIF+NTRC)*(K-1)
        DO KREC= KREC0+1,KREC0+NLAY+NDIF+NTRC
          READ(11,REC=KREC,IOSTAT=IOS) A
          IF     (IOS.NE.0) THEN
            WRITE(6,*) 'can''t read ',CFILEA(1:LEN_TRIM(CFILEA))
            CALL EXIT(4)
          ENDIF
#ifdef ENDIAN_IO
          CALL ENDIAN_SWAP(A(IP,JP),1)
#endif
          PROF(KREC-KREC0) = A(IP,JP)
        ENDDO
        write(6,*) 'PROF = ',PROF
        P(K+1) = P(K) + PROF(NLAY-3)
        UK = PROF(1)+UB
        VK = PROF(2)+VB
        UK = MAX( -9.99, MIN( 9.99, UK ) )
        VK = MAX( -9.99, MIN( 9.99, VK ) )
        IF     (NDIF.EQ.3) THEN
          VISC = MIN(9999.99, PROF(NLAY+NTRC+1)*1.E4)
          TDIF = MIN(9999.99, PROF(NLAY+NTRC+2)*1.E4)
        ENDIF
        WRITE(6,CFORMAT)
     &    K,
     &    UK*100.0,                   !cm/s
     &    VK*100.0,                   !cm/s
     &    PROF(NLAY-2),               !degC
     &    PROF(NLAY-1),               !psu
     &    PROF(NLAY)  +THBASE,        !SigmaT
     &    PROF(NLAY-3)/ONEM,          !m
     &    (P(K+1)+P(K))*0.5/ONEM,     !m
     &    VISC,TDIF,                  !cm^2/s
     &    (PROF(NLAY+KTR),KTR=1,NTRC) !0-99
      ENDDO
      END
      SUBROUTINE READ_B(CFILEB,
     &                  IEXPT,YRFLAG,IDM,JDM,KDM,NSURF,NLAY,NDIF,NTRC,
     &                  THBASE,TIME)
      IMPLICIT NONE
C
      INTEGER      IEXPT,YRFLAG,IDM,JDM,KDM,NSURF,NLAY,NDIF,NTRC
      REAL         THBASE,TIME
      CHARACTER*240 CFILEB
C
C     EXTRACT NEEDED MODEL PARAMETERS FROM ARCHIVE .b FILE.
C
      INTEGER      IOS,K,KTR,NSTEP
      REAL         THBASE_IN
      CHARACTER*6  CVARIN*6
      CHARACTER*240 CLINE
C
      OPEN(UNIT=12, FILE=CFILEB, FORM='FORMATTED', STATUS='OLD',
     +         IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',CFILEB(1:LEN_TRIM(CFILEB))
        WRITE(6,*) 'ios   = ',ios
        CALL EXIT(3)
      ENDIF
      READ(12,*)  ! skip title(1)
      READ(12,*)  ! skip title(2)
      READ(12,*)  ! skip title(3)
      READ(12,*)  ! skip title(4)
      READ(12,*)  ! skip iversn
      READ(12,*) IEXPT,CVARIN
      IF     (CVARIN.NE.'iexpt ') THEN
        WRITE(6,*) 'Error in hycom_profile: bad .b file'
        WRITE(6,*) 'filename: ',CFILEB(1:LEN_TRIM(CFILEB))
        CALL EXIT(4)
      endif
      READ(12,*) YRFLAG
      READ(12,*) IDM
      READ(12,*) JDM
C
C     FIND KDM.
C
      NSURF = 13  ! number of surface arrays
      NLAY  =  6  ! number of arrays per layer
C
      DO K= 1,10
        READ(12,'(a)') CLINE
      ENDDO
*     write(6,*) cline(1:len_trim(cline))
      IF     (CLINE(1:8).EQ.'thmix   ') THEN
        READ(CLINE(36:42),*) THBASE_IN
        IF     (THBASE_IN.NE.0.0) THEN
          THBASE = THBASE_IN
        ENDIF
      ELSE
        WRITE(6,*) 
        WRITE(6,*) 'Expected thmix but got:'
        WRITE(6,*) CLINE(1:LEN_TRIM(CLINE))
        WRITE(6,*) 
        CALL EXIT(2)
      ENDIF
      DO K= 11,13
        READ(12,'(a)') CLINE
      ENDDO
*     write(6,*) cline(1:len_trim(cline))
      IF     (CLINE(1:8).EQ.'kemix   ') THEN
        NLAY  = 7  ! mean archive
        NSURF = NSURF + 1
        READ(12,'(a)') CLINE
      ENDIF
*     write(6,*) cline(1:len_trim(cline))
      IF     (CLINE(1:8).EQ.'covice  ') THEN
        NSURF = NSURF + 3
        READ(12,'(a)') CLINE
        READ(12,'(a)') CLINE
        READ(12,'(a)') CLINE
      ENDIF
*     write(6,*) cline(1:len_trim(cline))
      READ(12,'(a)') CLINE
      IF     (NLAY.EQ.7) THEN
        NSURF = NSURF + 1
        READ(12,'(a)') CLINE  ! kebtrop
      ENDIF
*     write(6,*) cline(1:len_trim(cline))
C
      K = INDEX(CLINE,'=')
      READ(CLINE(K+1:),*) NSTEP,TIME
C
C     FIRST LAYER (FIND NUMBER OF LAYERS)
C
      READ(12,'(a)') CLINE
      READ(12,'(a)') CLINE
      READ(12,'(a)') CLINE
      READ(12,'(a)') CLINE
      READ(12,'(a)') CLINE
      READ(12,'(a)') CLINE
      IF     (NLAY.EQ.7) THEN
        READ(12,'(a)') CLINE
      ENDIF
      NDIF = 0
      DO KTR= 1,9999
        READ(12,'(a)',IOSTAT=IOS) CLINE
        IF     (IOS.NE.0) THEN
C
C         EXACTLY ONE LAYER
C
          NTRC = KTR-NDIF-1
          KDM  = 1
*         write(6,*) 'kdm,ndif,ntrc = ',kdm,ndif,ntrc
          CLOSE(UNIT=12)
          RETURN
        ELSEIF (CLINE(1:8).EQ.'viscty  ') THEN
          NDIF = 3
        ELSEIF (CLINE(1:8).EQ.'t-diff  ') THEN
          CYCLE
        ELSEIF (CLINE(1:8).EQ.'s-diff  ') THEN
          CYCLE
        ELSEIF (CLINE(1:8).NE.'tracer  ') THEN
          EXIT
        ENDIF
      ENDDO
      NTRC = KTR-NDIF-1
C
C     MORE THAN 1 LAYER
C
*     READ(12,'(a)') CLINE  ! done above
      READ(12,'(a)') CLINE
      READ(12,'(a)') CLINE
      READ(12,'(a)') CLINE
      READ(12,'(a)') CLINE
      READ(12,'(a)') CLINE
      IF     (NLAY.EQ.7) THEN
        READ(12,'(a)') CLINE
      ENDIF
      DO KTR= 1,NDIF+NTRC
        READ(12,'(a)') CLINE
      ENDDO
      DO K= 3,999
        READ(12,'(a)',IOSTAT=IOS) CLINE
        IF     (IOS.NE.0) THEN
          EXIT
        ELSEIF (CLINE(1:8).NE.'u-vel.  ') THEN
*         write(6,*) cline(1:len_trim(cline))
          EXIT
        ENDIF
*       write(6,*) cline(1:len_trim(cline)),"  (1)"
C
        READ(12,'(a)') CLINE
        READ(12,'(a)') CLINE
        READ(12,'(a)') CLINE
        READ(12,'(a)') CLINE
        READ(12,'(a)') CLINE
        IF     (NLAY.EQ.7) THEN
          READ(12,'(a)') CLINE
        ENDIF
        DO KTR= 1,NDIF+NTRC
          READ(12,'(a)') CLINE
        ENDDO
      ENDDO
      KDM = K-1
      write(6,*) 'kdm,ndif,ntrc = ',kdm,ndif,ntrc
      CLOSE(UNIT=12)
      RETURN
      END
