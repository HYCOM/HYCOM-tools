      PROGRAM HYCOM_PROFILE_OBS
      IMPLICIT NONE
C
C  hycom_profile_obs - Usage:  hycom_profile_obs archv.a list.txt obs.txt [itype [bot]]
C
C                 produce archive field values at a list of points (x,y,z)
C
C   archv.a   is assumed to be an HYCOM archive data file, with companion
C             header file archv.b.
C   list.txt  contains a list of points,   one per line
C   obs.txt   will be the output samples,  one per line
C   itype     is the interpolation type (default 1)
C                =-2; piecewise quadratic cross cell
C                =-1; piecewise linear   across cell
C                = 0; piecewise constant across cell
C                = 1; linear interpolation between cell centers
C   bot       ignore layers within bot of the bottom (default 0.0)
C
C  Each line of list.txt should contain three values: x,y,z
C  The x,y values are the p-grid location of the point,
C  note that hycom/bin/hycom_lonlat2xy will convert lon,lat to x,y.
C  The z value is the depth of the point in m.
C
C  It is (slightly) more efficient to group points with the same
C  x,y together with z increasing, i.e. with z(i)<=z(i+1).
C
C  A value of 99.9999 in the output obs.txt indicates a data void.
C  If SSH (and U,V,T,S) is a data void, the x,y location is over land,
C  otherwise U,V,T,S are data voids when the z depth is below the bottom.
C
C  The files regional.grid.[ab] for this domain must be in the current
C  directory.  They are used to calculate lon,lat values for each point.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  December 2015.
C
      REAL*4     ONEM,SPVAL,OSPVAL
      PARAMETER (ONEM=9806.0, SPVAL=2.0**100, OSPVAL=99.9999)
C
      REAL*4, ALLOCATABLE :: A(:,:),AA(:,:)
      REAL*4, ALLOCATABLE :: XP(:),YP(:)
      REAL*4, ALLOCATABLE :: PLON(:),PLAT(:),PANG(:)
      REAL*4, ALLOCATABLE :: P(:,:),TK(:,:),SK(:,:),
     +                       UK(:,:),VK(:,:),
     +                       UB(:),VB(:),FSD(:)
      REAL*4, ALLOCATABLE :: ZZ(:),TZ(:),SZ(:),
     +                       UZ(:),VZ(:)
      REAL*4              :: PAD(4096)
C
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      LOGICAL       LSTERIC,LSEAICE,LFATAL
      INTEGER       IDM,JDM,KDM,KTR,NSURF,NLAY,NDIF,NTRC,
     +              ARTYPE,IEXPT,YRFLAG,SIGVER,I_ARCH(17)
      INTEGER       ITYPE,IDEBUG,KZ,NPAD,M,N,NX
      REAL*4        DUMZ,FLAG,THBASE,BOT
      INTEGER       NSTEP
      REAL*8        TIME3(3)
      CHARACTER*240 CFILEA,CFILEB,CFILEL,CFILEO,CFILEVS,CFORMAT
C
      CHARACTER*18 CASN
      LOGICAL      LDONE
      INTEGER      I,II,IP,IPP1,IPP2,J,JJ,JP,JPP1,JPP2,K,LANDF,LANDL,
     +             KREC,KREC0,IOS,NRECL
      REAL*4       A00,A01,A10,A11,AXY,DX,DY
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
        CALL GETARG(2,CFILEL)
        CALL GETARG(3,CFILEO)
        ITYPE  = 1
        IDEBUG = 0
        BOT    = 0.0
      ELSEIF (NARG.EQ.4) THEN
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CFILEL)
        CALL GETARG(3,CFILEO)
        CALL GETARG(4,CARG)
        READ(CARG,*) ITYPE
        IDEBUG = 0
        BOT    = 0.0
      ELSEIF (NARG.EQ.5) THEN 
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CFILEL)
        CALL GETARG(3,CFILEO)
        CALL GETARG(4,CARG)
        READ(CARG,*) ITYPE
        CALL GETARG(5,CARG)
        READ(CARG,*) BOT
        IDEBUG = 0
      ELSEIF (NARG.EQ.6) THEN !undocumented debug option
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CFILEL)
        CALL GETARG(3,CFILEO)
        CALL GETARG(4,CARG)
        READ(CARG,*) ITYPE
        CALL GETARG(5,CARG)
        READ(CARG,*) BOT
        CALL GETARG(6,CARG)
        READ(CARG,*) IDEBUG  !>0 print this station
      ELSE
        WRITE(6,"(2a)") 
     +    'Usage: hycom_profile_obs',
     +    ' archv.a list.txt obs.txt [itype [bot]]'
        CALL EXIT(1)
        STOP
      ENDIF
C
C     EXTRACT MODEL PARAMETERS FROM ".b" FILE.
C
      cfilevs = ' '
      call getenv('ARCHVS',cfilevs)
      if     (cfilevs.eq.' ') then
        CFILEB = CFILEA(1:LEN_TRIM(CFILEA)-1) // 'b'
        CALL READ_BS(CFILEB,
     +               IEXPT,YRFLAG,IDM,JDM,KDM,NSURF,NLAY,NDIF,NTRC,
     +               LSTERIC,ARTYPE,SIGVER,THBASE,TIME3,NSTEP)
        I_ARCH( :) = -1  !turn off i_arch
      else
        CFILEB = CFILEA(1:LEN_TRIM(CFILEA)-1) // 'b'
        CALL READ_BSCS(CFILEB,CFILEVS,I_ARCH,
     +                 IEXPT,YRFLAG,IDM,JDM,KDM,NSURF,NLAY,NDIF,NTRC,
     +                 LSEAICE,LSTERIC,ARTYPE,SIGVER,THBASE,TIME3,NSTEP)
      endif
*     write(6,*) 'IDM = ',idm
*     write(6,*) 'JDM = ',jdm
*     write(6,*) 'KDM = ',kdm
C
      CALL SIG_I(SIGVER)
C
C     FIND THE NUMBER OF PROFILES (SAMPLE LOCATIONS)
C
      OPEN(UNIT=31, FILE=CFILEL, FORM='FORMATTED', STATUS='OLD',
     +         IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',TRIM(CFILEL)
        WRITE(6,*) 'ios   = ',ios
        CALL EXIT(3)
        STOP
      ENDIF
C
      DO N= 1,99999
        READ(31,*,IOSTAT=IOS) DX,DY,DUMZ
        IF     (IOS.NE.0) THEN
          EXIT
        ENDIF
      ENDDO !n
      N = N-1
      REWIND(31)
*     write(6,*) 'N = ',n
C
C     ALLOCATE ARRAYS.
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( A(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_profile_obs: could not allocate ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
        STOP
      ENDIF
      ALLOCATE( AA(IDM+1,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_profile_obs: could not allocate ',
     +             (IDM+1)*JDM,' words'
        CALL EXIT(2)
        STOP
      ENDIF
      ALLOCATE(    P(0:KDM,N),
     +            TK(1:KDM,N),
     +            SK(1:KDM,N),
     +            UK(1:KDM,N),
     +            VK(1:KDM,N), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_profile_obs: could not allocate ',
     +             '(1:kdm,1:n) K arrays'
        CALL EXIT(2)
        STOP
      ENDIF
      ALLOCATE( ZZ(      N),
     +          TZ(      N),
     +          SZ(      N),
     +          UZ(      N),
     +          VZ(      N), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_profile_obs: could not allocate ',
     +             '(1:n) Z arrays'
        CALL EXIT(2)
        STOP
      ENDIF
      ALLOCATE(   UB(N),
     +            VB(N),
     +           FSD(N),
     +          PLON(N),
     +          PLAT(N),
     +          PANG(N),
     +            XP(N),
     +            YP(N), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_profile_obs: could not allocate ',
     +             '(1:n) P arrays'
        CALL EXIT(2)
        STOP
      ENDIF
C
C     READ PROFILE LOCATIONS
C
      DO I= 1,N
        READ(31,*,IOSTAT=IOS) XP(I),YP(I),ZZ(I)
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'Error: can''t read ',TRIM(CFILEL)
          WRITE(6,*) 'ios   = ',ios
          WRITE(6,*) 'i,n   = ',i,n
          CALL EXIT(3)
          STOP
        ENDIF
      ENDDO !i
      CLOSE(31)
C
C     OPEN "regional.grid.a" FILE.
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
          STOP
        ENDIF
      ENDIF
#else
      CALL ASNUNIT(11,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t asnunit 11'
        WRITE(6,*) 'ios = ',ios
        CALL EXIT(5)
        STOP
      ENDIF
#endif
#endif
      OPEN(UNIT=11, FILE="regional.grid.a",
     +         FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open regional.grid.a'
        WRITE(6,*) 'ios   = ',ios
        WRITE(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
        STOP
      ENDIF
C
      READ(11,REC=1,IOSTAT=IOS) A  !plon
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read regional.grid.a'
        WRITE(6,*) 'ios  = ',ios
        CALL EXIT(5)
        STOP
      ENDIF
      AA(1:IDM,:) = A(:,:)
      DO J= 1,JDM
        A00 = MOD(AA(1,J) - AA(IDM,J),360.0)
        IF     (A00.LT.-180.0) THEN
          A00 = A00 + 360.0
        ELSEIF (A00.GT. 180.0) THEN
          A00 = A00 - 360.0
        ENDIF
        AA(IDM+1,J) = AA(IDM,J) + A00
      ENDDO !j
      DO I= 1,N
        IP = XP(I)
        JP = YP(I)
        DX = XP(I) - IP
        DY = YP(I) - JP
        IF     (DX.EQ.0.0) THEN
          IPP1 = IP
        ELSEIF (DX.EQ.1.0) THEN
          IPP1 = IP+1
          IP   = IPP1
        ELSE
          IPP1 = IP+1
        ENDIF
        IF     (DY.EQ.0.0) THEN
          JPP1 = JP
        ELSEIF (DY.EQ.1.0) THEN
          JPP1 = JP+1
          JP   = JPP1
        ELSE
          JPP1 = JP+1
        ENDIF
        AXY = (1.0-DX)*(1.0-DY)*AA(IP,  JP)   +
     +        (1.0-DX)*     DY *AA(IP,  JPP1) +
     +             DX *(1.0-DY)*AA(IPP1,JP)   +
     +             DX *     DY *AA(IPP1,JPP1)
        PLON(I) = AXY
        PLON(I) = MOD( PLON(I) + 1080.0, 360.0 )
        IF     (PLON(I).GT.180.0) THEN
          PLON(I) = PLON(I) - 360.0
        ENDIF
        IF     (IDEBUG.EQ.I) THEN
          WRITE(6,*) 'debug: i,plon = ',I,PLON(I)
        ENDIF
      ENDDO !i
C
      READ(11,REC=2,IOSTAT=IOS) A  !plat
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read regional.grid.a'
        WRITE(6,*) 'ios  = ',ios
        CALL EXIT(5)
        STOP
      ENDIF
      AA(1:IDM,:) = A(:,:)
      DO J= 1,JDM
        AA(IDM+1,J) = AA(IDM,J)
      ENDDO !j
      DO I= 1,N
        IP = XP(I)
        JP = YP(I)
        DX = XP(I) - IP
        DY = YP(I) - JP
        IF     (DX.EQ.0.0) THEN
          IPP1 = IP
        ELSEIF (DX.EQ.1.0) THEN
          IPP1 = IP+1
          IP   = IPP1
        ELSEIF (IP.EQ.IDM) THEN
          IPP1 = 1
        ELSE
          IPP1 = IP+1
        ENDIF
        IF     (DY.EQ.0.0) THEN
          JPP1 = JP
        ELSEIF (DY.EQ.1.0) THEN
          JPP1 = JP+1
          JP   = JPP1
        ELSE
          JPP1 = JP+1
        ENDIF
        AXY = (1.0-DX)*(1.0-DY)*AA(IP,  JP)   +
     +        (1.0-DX)*     DY *AA(IP,  JPP1) +
     +             DX *(1.0-DY)*AA(IPP1,JP)   +
     +             DX *     DY *AA(IPP1,JPP1)
        PLAT(I) = AXY
        IF     (IDEBUG.EQ.I) THEN
          WRITE(6,*) 'debug: i,plat = ',I,PLAT(I)
        ENDIF
      ENDDO !i
C
      READ(11,REC=9,IOSTAT=IOS) A  !pang
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read regional.grid.a'
        WRITE(6,*) 'ios  = ',ios
        CALL EXIT(5)
        STOP
      ENDIF
      DO I= 1,N
        IP = NINT(XP(I))
        JP = NINT(YP(I))
        IF     (IP.GT.IDM) THEN
          IP = 1 !periodic wrap
        ENDIF
        IF     (JP.GT.JDM) THEN
          JP = JDM
        ENDIF
        PANG(I) = A(IP,JP)  !no interpolation of angle, use nearest point
        IF     (IDEBUG.EQ.I) THEN
          WRITE(6,*) 'debug: i,pang = ',I,PANG(I)
        ENDIF
      ENDDO !i
      CLOSE(UNIT=11)
      DEALLOCATE( AA )
C
C     OPEN ARCHIVE ".a" FILE.
C
      OPEN(UNIT=11, FILE=CFILEA, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',TRIM(CFILEA)
        WRITE(6,*) 'ios   = ',ios
        WRITE(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
        STOP
      ENDIF
C
C     SURFACE FIELDS (JUST FSD,UB,VB).
C
      IOS = 0
C
      IF     (NLAY.LE.6) THEN
        LFATAL = .FALSE.
        IF     (I_ARCH(2).LT.0) THEN
          READ(11,REC=2,IOSTAT=IOS) A
#ifdef ENDIAN_IO
          CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        ELSEIF (I_ARCH(2).GT.0) THEN
          READ(11,REC=I_ARCH(2),IOSTAT=IOS) A
#ifdef ENDIAN_IO
          CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        ELSE
          A(:,:) = 0.0
        ENDIF
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'can''t read ',TRIM(CFILEA)
          CALL EXIT(4)
          STOP
        ENDIF
        LANDF = 0
        LANDL = 0
        DO I= 1,N
          IP = XP(I)
          JP = YP(I)
          DX = XP(I) - IP
          DY = YP(I) - JP
          IF     (DX.EQ.0.0) THEN !1-d interpolation
            IPP1 = IP
          ELSEIF (DX.EQ.1.0) THEN !1-d interpolation
            IPP1 = IP+1
            IP   = IPP1
          ELSEIF (IP.EQ.IDM) THEN
            IPP1 = 1
          ELSE
            IPP1 = IP+1
          ENDIF
          IF     (DY.EQ.0.0) THEN !1-d interpolation
            JPP1 = JP
          ELSEIF (DY.EQ.1.0) THEN !1-d interpolation
            JPP1 = JP+1
            JP   = JPP1
          ELSE
            JPP1 = JP+1
          ENDIF
          IF     (MAX( A(IP,  JP),
     +                 A(IP,  JPP1),
     +                 A(IPP1,JP),
     +                 A(IPP1,JPP1) ).GE.SPVAL) THEN
C
C           NEAR LAND, IS NEAREST POINT ACTUALLY LAND?
C
            IF     (DX.LE.0.5 .AND. DY.LE.0.5) THEN
              IP = IP
              JP = JP
            ELSEIF (DX.LE.0.5 .AND. DY.GT.0.5) THEN
              IP = IP
              JP = JPP1
            ELSEIF (DX.GT.0.5 .AND. DY.LE.0.5) THEN
              IP = IPP1
              JP = JP
            ELSE ! (DX.GT.0.5 .AND. DY.GT.0.5) THEN
              IP = IPP1
              JP = JPP1
            ENDIF
*           IF     (.TRUE.) THEN !require 4 good points
            IF     (A(IP,JP).GE.SPVAL) THEN
C
C             NEAREST POINT IS LAND.
C
              FSD(I) = OSPVAL
              IF     (LANDF.EQ.0) THEN
                LANDF = I
              ENDIF
              LANDL = I
              IF     (IDEBUG.LT.0) THEN
                WRITE(6,'(a,2i6)') 'nearest point is land: ',IP,JP
                LDONE = .FALSE.
                DO JJ= MAX(1,JP-1),MIN(JDM,JP+1)
                  DO II= MAX(1,IP-1),MIN(IDM,IP+1)
                    IF     (A(II,JJ).LT.SPVAL) THEN
                      LDONE = .TRUE.
                      WRITE(6,'(a,2i6)') '  but sea at: ',II,JJ
                    ENDIF
                  ENDDO
                ENDDO
                IF     (.NOT.LDONE) THEN !25-pt search
                  DO JJ= MAX(1,JP-2),MIN(JDM,JP+2)
                    DO II= MAX(1,IP-2),MIN(IDM,IP+2)
                      IF     (A(II,JJ).LT.SPVAL) THEN
                        LDONE = .TRUE.
                        WRITE(6,'(a,2i6)') '  but sea at: ',II,JJ
                      ENDIF
                    ENDDO
                  ENDDO
                ENDIF
              ENDIF !IDEBUG
            ELSE
C
C             USE NEAREST POINT AS THE EXACT VALUE.
C
              IF     (LANDF.NE.0) THEN
                WRITE(6,'(a,i6,a,i6,a,i6,a)')
     +            'WARNING - LOCATIONS',LANDF,
     +                            ' TO',LANDL,
     +                            ' OF',N,    ' ARE TREATED AS LAND'
                LANDF = 0
              ENDIF
              FSD(I) = A(IP,JP)/9.806
               XP(I) = IP  !for all subsequent interpolations
               YP(I) = JP  !for all subsequent interpolations
            ENDIF
          ELSE  !all 4 points ok
            IF     (LANDF.NE.0) THEN
              WRITE(6,'(a,i6,a,i6,a,i6,a)')
     +          'WARNING - LOCATIONS',LANDF,
     +                          ' TO',LANDL,
     +                          ' OF',N,    ' ARE TREATED AS LAND'
              LANDF = 0
            ENDIF
            AXY = (1.0-DX)*(1.0-DY)*A(IP,  JP)   +
     +            (1.0-DX)*     DY *A(IP,  JPP1) +
     +                 DX *(1.0-DY)*A(IPP1,JP)   +
     +                 DX *     DY *A(IPP1,JPP1)
            FSD(I) = AXY/9.806
          ENDIF
          IF     (IDEBUG.EQ.I) THEN
            WRITE(6,*) 'debug: i,k,fsd= ',I,K,FSD(I)
          ENDIF
        ENDDO !i
        IF     (LANDF.NE.0) THEN
          WRITE(6,'(a,i6,a,i6,a,i6,a)')
     +      'WARNING - LOCATIONS',LANDF,
     +                      ' TO',LANDL,
     +                      ' OF',N,    ' ARE TREATED AS LAND'
        ENDIF
        IF     (I_ARCH(11).LT.0) THEN
          READ(11,REC=NSURF-1,IOSTAT=IOS) A
#ifdef ENDIAN_IO
          CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        ELSEIF (I_ARCH(11).GT.0) THEN
          READ(11,REC=I_ARCH(11),IOSTAT=IOS) A
#ifdef ENDIAN_IO
          CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        ELSE
          A(:,:) = 0.0
        ENDIF
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'can''t read ',TRIM(CFILEA)
          CALL EXIT(4)
          STOP
        ENDIF
        DO I= 1,N
          IF     (FSD(I).EQ.OSPVAL) THEN
            CYCLE
          ENDIF
          IP = XP(I)
          JP = YP(I)
          DX = XP(I) - IP
          DY = YP(I) - JP
          IF     (DX.EQ.0.0) THEN
            IPP1 = IP
          ELSEIF (DX.EQ.1.0) THEN
            IPP1 = IP+1
            IP   = IPP1
          ELSEIF (IP.EQ.IDM) THEN
            IPP1 = 1
          ELSE
            IPP1 = IP+1
          ENDIF
          IF     (IPP1.EQ.IDM) THEN
            IPP2 = 1
          ELSE
            IPP2 = IPP1+1
          ENDIF
          IF     (DY.EQ.0.0) THEN
            JPP1 = JP
          ELSEIF (DY.EQ.1.0) THEN
            JPP1 = JP+1
            JP   = JPP1
          ELSE
            JPP1 = JP+1
          ENDIF
          IF (A(IP,  JP)  .EQ.SPVAL) A(IP,  JP)   = 0.0
          IF (A(IPP1,JP)  .EQ.SPVAL) A(IPP1,JP)   = 0.0
          IF (A(IPP2,JP)  .EQ.SPVAL) A(IPP2,JP)   = 0.0
          IF (A(IP,  JPP1).EQ.SPVAL) A(IP,  JPP1) = 0.0
          IF (A(IPP1,JPP1).EQ.SPVAL) A(IPP1,JPP1) = 0.0
          IF (A(IPP2,JPP1).EQ.SPVAL) A(IPP2,JPP1) = 0.0
          IF     (ARTYPE.GT.0) THEN
            A00 = 0.5*(A(IP,  JP)   + A(IP+1,JP)   )
            A01 = 0.5*(A(IP,  JPP1) + A(IP+1,JPP1) )
            A10 = 0.5*(A(IPP1,JP)   + A(IPP2,JP)   )
            A11 = 0.5*(A(IPP1,JPP1) + A(IPP2,JPP1) )
          ELSE !p-grid
            A00 =      A(IP,  JP)
            A01 =      A(IP,  JPP1)
            A10 =      A(IPP1,JP)
            A11 =      A(IPP1,JPP1)
          ENDIF !artype
          AXY = (1.0-DX)*(1.0-DY)*A00 +
     +          (1.0-DX)*     DY *A01 +
     +               DX *(1.0-DY)*A10 +
     +               DX *     DY *A11
          UB( I) = AXY
          IF     (A00.EQ.0.0 .OR. A01.EQ.0.0 .OR.
     +            A10.EQ.0.0 .OR. A11.EQ.0.0     ) THEN
*           LFATAL = .TRUE.
*           WRITE(6,*) 'ERROR - STATION ',I,' HAS BAD U',
*    +                 IP,IPP1,IPP2,JP,JPP1
            WRITE(6,*) 'WARNING - STATION ',I,' HAS BAD U',
     +                 IP,IPP1,IPP2,JP,JPP1
            WRITE(6,*) 
     +        'DX,DY = ',DX,DY
            WRITE(6,*) 
     +        'U.0 = ',A(IP,  JP),  A(IPP1,JP),  A(IPP2,JP)
            WRITE(6,*) 
     +        'U.1 = ',A(IP,  JPP1),A(IPP1,JPP1),A(IPP2,JPP1)
          ENDIF
        ENDDO !i
C
        IF     (I_ARCH(12).LT.0) THEN
          READ(11,REC=NSURF,  IOSTAT=IOS) A
#ifdef ENDIAN_IO
          CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        ELSEIF (I_ARCH(12).GT.0) THEN
          READ(11,REC=I_ARCH(12),IOSTAT=IOS) A
#ifdef ENDIAN_IO
          CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        ELSE
          A(:,:) = 0.0
        ENDIF
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'can''t read ',TRIM(CFILEA)
          CALL EXIT(4)
          STOP
        ENDIF
        DO I= 1,N
          IF     (FSD(I).EQ.OSPVAL) THEN
            CYCLE
          ENDIF
          IP = XP(I)
          JP = YP(I)
          DX = XP(I) - IP
          DY = YP(I) - JP
          IF     (DX.EQ.0.0) THEN
            IPP1 = IP
          ELSEIF (DX.EQ.1.0) THEN
            IPP1 = IP+1
            IP   = IPP1
          ELSEIF (IP.EQ.IDM) THEN
            IPP1 = 1
          ELSE
            IPP1 = IP+1
          ENDIF
          IF     (DY.EQ.0.0) THEN
            JPP1 = JP
          ELSEIF (DY.EQ.1.0) THEN
            JPP1 = JP+1
            JP   = JPP1
          ELSE
            JPP1 = JP+1
          ENDIF
          IF     (JPP1.NE.JDM) THEN
            JPP2 = JPP1 + 1
          ELSE
            JPP2 = JPP1
          ENDIF
          IF (A(IP,  JP)  .EQ.SPVAL) A(IP,  JP)   = 0.0
          IF (A(IP,  JP+1).EQ.SPVAL) A(IP,  JP+1) = 0.0
          IF (A(IP,  JPP2).EQ.SPVAL) A(IP,  JPP2) = 0.0
          IF (A(IPP1,JP)  .EQ.SPVAL) A(IPP1,JP)   = 0.0
          IF (A(IPP1,JP+1).EQ.SPVAL) A(IPP1,JP+1) = 0.0
          IF (A(IPP1,JPP2).EQ.SPVAL) A(IPP1,JPP2) = 0.0
          IF     (ARTYPE.GT.0) THEN
            A00 = 0.5*(A(IP,  JP)   + A(IP,  JP+1) )
            A01 = 0.5*(A(IP,  JPP1) + A(IP,  JPP2) )
            A10 = 0.5*(A(IPP1,JP)   + A(IPP1,JP+1) )
            A11 = 0.5*(A(IPP1,JPP1) + A(IPP1,JPP2) )
          ELSE !p-grid
            A00 =      A(IP,  JP)
            A01 =      A(IP,  JPP1)
            A10 =      A(IPP1,JP)
            A11 =      A(IPP1,JPP1)
          ENDIF !artype
          AXY = (1.0-DX)*(1.0-DY)*A00 +
     +          (1.0-DX)*     DY *A01 +
     +               DX *(1.0-DY)*A10 +
     +               DX *     DY *A11
          VB( I) = AXY
          P(0,I) = 0.0
          IF     (A00.EQ.0.0 .OR. A01.EQ.0.0 .OR.
     +            A10.EQ.0.0 .OR. A11.EQ.0.0     ) THEN
*           LFATAL = .TRUE.
*           WRITE(6,*) 'ERROR - STATION ',I,' HAS BAD V',
*    +                 IP,IPP1,JP,JPP1,JPP2
            WRITE(6,*) 'WARNING - STATION ',I,' HAS BAD V',
     +                 IP,IPP1,JP,JPP1,JPP2
            WRITE(6,*) 
     +        'DX,DY = ',DX,DY
            WRITE(6,*) 
     +        'V.0 = ',A(IP,  JP),A(IP,  JP+1),A(IP,  JPP2)
            WRITE(6,*) 
     +        'V.1 = ',A(IPP1,JP),A(IPP1,JP+1),A(IPP1,JPP2)
          ENDIF
        ENDDO !i
        IF     (LFATAL) THEN
          CALL EXIT(9)
          STOP
        ENDIF
      ELSE  !mean archive
        DO I= 1,N
          UB( I) = 0.0
          VB( I) = 0.0
          P(0,I) = 0.0
        ENDDO !i
      ENDIF
C
C     LAYER FIELDS.
C
      DO K= 1,KDM
        KREC0 = NSURF+(NLAY+NDIF+NTRC)*(K-1)
c
c ---   u
        IF     (I_ARCH(13).LT.0) THEN
          READ(11,REC=KREC0+1,IOSTAT=IOS) A
#ifdef ENDIAN_IO
          CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        ELSEIF (I_ARCH(13).GT.0) THEN
          READ(11,REC=I_ARCH(13),IOSTAT=IOS) A
#ifdef ENDIAN_IO
          CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        ELSE
          A(:,:) = 0.0
        ENDIF
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'can''t read ',TRIM(CFILEA)
          CALL EXIT(4)
          STOP
        ENDIF
        DO I= 1,N
          IF     (FSD(I).EQ.OSPVAL) THEN
            CYCLE
          ENDIF
          IP = XP(I)
          JP = YP(I)
          DX = XP(I) - IP
          DY = YP(I) - JP
          IF     (DX.EQ.0.0) THEN
            IPP1 = IP
          ELSEIF (DX.EQ.1.0) THEN
            IPP1 = IP+1
            IP   = IPP1
          ELSEIF (IP.EQ.IDM) THEN
            IPP1 = 1
          ELSE
            IPP1 = IP+1
          ENDIF
          IF     (IPP1.EQ.IDM) THEN
            IPP2 = 1
          ELSE
            IPP2 = IPP1+1
          ENDIF
          IF     (DY.EQ.0.0) THEN
            JPP1 = JP
          ELSEIF (DY.EQ.1.0) THEN
            JPP1 = JP+1
            JP   = JPP1
          ELSE
            JPP1 = JP+1
          ENDIF
          IF (A(IP,  JP)  .EQ.SPVAL) A(IP,  JP)   = 0.0
          IF (A(IPP1,JP)  .EQ.SPVAL) A(IPP1,JP)   = 0.0
          IF (A(IPP2,JP)  .EQ.SPVAL) A(IPP2,JP)   = 0.0
          IF (A(IP,  JPP1).EQ.SPVAL) A(IP,  JPP1) = 0.0
          IF (A(IPP1,JPP1).EQ.SPVAL) A(IPP1,JPP1) = 0.0
          IF (A(IPP2,JPP1).EQ.SPVAL) A(IPP2,JPP1) = 0.0
          IF     (ARTYPE.GT.0) THEN
            A00 = 0.5*(A(IP,  JP)   + A(IP+1,  JP)   )
            A01 = 0.5*(A(IP,  JPP1) + A(IP+1,  JPP1) )
            A10 = 0.5*(A(IPP1,JP)   + A(IPP2,  JP)   )
            A11 = 0.5*(A(IPP1,JPP1) + A(IPP2,  JPP1) )
          ELSE !p-grid
            A00 =      A(IP,  JP)
            A01 =      A(IP,  JPP1)
            A10 =      A(IPP1,JP)
            A11 =      A(IPP1,JPP1)
          ENDIF !artype
          AXY = (1.0-DX)*(1.0-DY)*A00 +
     +          (1.0-DX)*     DY *A01 +
     +               DX *(1.0-DY)*A10 +
     +               DX *     DY *A11
          UK(K,I) = UB(I) + AXY
        ENDDO !i
c
c ---   v
        IF     (I_ARCH(14).LT.0) THEN
          READ(11,REC=KREC0+2,IOSTAT=IOS) A
#ifdef ENDIAN_IO
          CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        ELSEIF (I_ARCH(14).GT.0) THEN
          READ(11,REC=I_ARCH(14),IOSTAT=IOS) A
#ifdef ENDIAN_IO
          CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        ELSE
          A(:,:) = 0.0
        ENDIF
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'can''t read ',TRIM(CFILEA)
          CALL EXIT(4)
          STOP
        ENDIF
        DO I= 1,N
          IF     (FSD(I).EQ.OSPVAL) THEN
            CYCLE
          ENDIF
          IP = XP(I)
          JP = YP(I)
          DX = XP(I) - IP
          DY = YP(I) - JP
          IF     (DX.EQ.0.0) THEN
            IPP1 = IP
          ELSEIF (DX.EQ.1.0) THEN
            IPP1 = IP+1
            IP   = IPP1
          ELSEIF (IP.EQ.IDM) THEN
            IPP1 = 1
          ELSE
            IPP1 = IP+1
          ENDIF
          IF     (DY.EQ.0.0) THEN
            JPP1 = JP
          ELSEIF (DY.EQ.1.0) THEN
            JPP1 = JP+1
            JP   = JPP1
          ELSE
            JPP1 = JP+1
          ENDIF
          IF     (JPP1.NE.JDM) THEN
            JPP2 = JPP1 + 1
          ELSE
            JPP2 = JPP1
          ENDIF
          IF (A(IP,  JP)  .EQ.SPVAL) A(IP,  JP)   = 0.0
          IF (A(IP,  JP+1).EQ.SPVAL) A(IP,  JP+1) = 0.0
          IF (A(IP,  JPP2).EQ.SPVAL) A(IP,  JPP2) = 0.0
          IF (A(IPP1,JP)  .EQ.SPVAL) A(IPP1,JP)   = 0.0
          IF (A(IPP1,JP+1).EQ.SPVAL) A(IPP1,JP+1) = 0.0
          IF (A(IPP1,JPP2).EQ.SPVAL) A(IPP1,JPP2) = 0.0
          IF     (ARTYPE.GT.0) THEN
            A00 = 0.5*(A(IP,  JP)   + A(IP,  JP+1) )
            A01 = 0.5*(A(IP,  JPP1) + A(IP,  JPP2) )
            A10 = 0.5*(A(IPP1,JP)   + A(IPP1,JP+1) )
            A11 = 0.5*(A(IPP1,JPP1) + A(IPP1,JPP2) )
          ELSE !p-grid
            A00 =      A(IP,  JP)
            A01 =      A(IP,  JPP1)
            A10 =      A(IPP1,JP)
            A11 =      A(IPP1,JPP1)
          ENDIF !artype
          AXY = (1.0-DX)*(1.0-DY)*A00 +
     +          (1.0-DX)*     DY *A01 +
     +               DX *(1.0-DY)*A10 +
     +               DX *     DY *A11
          VK(K,I) = VB(I) + AXY
        ENDDO !i
c
c ---   p
        IF     (ABS(ARTYPE).EQ.1) THEN
          IF     (I_ARCH(15).LT.0) THEN
            READ(11,REC=KREC0+3,IOSTAT=IOS) A
#ifdef ENDIAN_IO
            CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
          ELSEIF (I_ARCH(15).GT.0) THEN
            READ(11,REC=I_ARCH(15),IOSTAT=IOS) A
#ifdef ENDIAN_IO
            CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
          ELSE
            A(:,:) = 0.0
          ENDIF
        ELSE  !mean archive
          READ(11,REC=KREC0+4,IOSTAT=IOS) A
#ifdef ENDIAN_IO
          CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        ENDIF
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'can''t read ',TRIM(CFILEA)
          CALL EXIT(4)
          STOP
        ENDIF
        DO I= 1,N
          IF     (FSD(I).EQ.OSPVAL) THEN
            CYCLE
          ENDIF
          IP = XP(I)
          JP = YP(I)
          DX = XP(I) - IP
          DY = YP(I) - JP
          IF     (DX.EQ.0.0) THEN
            IPP1 = IP
          ELSEIF (DX.EQ.1.0) THEN
            IPP1 = IP+1
            IP   = IPP1
          ELSEIF (IP.EQ.IDM) THEN
            IPP1 = 1
          ELSE
            IPP1 = IP+1
          ENDIF
          IF     (DY.EQ.0.0) THEN
            JPP1 = JP
          ELSEIF (DY.EQ.1.0) THEN
            JPP1 = JP+1
            JP   = JPP1
          ELSE
            JPP1 = JP+1
          ENDIF
          AXY = (1.0-DX)*(1.0-DY)*A(IP,  JP)   +
     +          (1.0-DX)*     DY *A(IP,  JPP1) +
     +               DX *(1.0-DY)*A(IPP1,JP)   +
     +               DX *     DY *A(IPP1,JPP1)
          P(K,I) = P(K-1,I) + AXY/ONEM  !interface depth in m
          IF     (IDEBUG.EQ.I) THEN
            WRITE(6,*) 'debug: i,k,p  = ',I,K,P(K,I)
          ENDIF
        ENDDO !i
c
c ---   t
        IF     (ABS(ARTYPE).EQ.1) THEN
          IF     (I_ARCH(16).LT.0) THEN
            READ(11,REC=KREC0+4,IOSTAT=IOS) A
#ifdef ENDIAN_IO
            CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
          ELSEIF (I_ARCH(16).GT.0) THEN
            READ(11,REC=I_ARCH(16),IOSTAT=IOS) A
#ifdef ENDIAN_IO
            CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
          ELSE
            A(:,:) = 0.0
          ENDIF
        ELSE  !mean archive
          READ(11,REC=KREC0+5,IOSTAT=IOS) A
#ifdef ENDIAN_IO
          CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        ENDIF
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'can''t read ',TRIM(CFILEA)
          CALL EXIT(4)
          STOP
        ENDIF
        DO I= 1,N
          IF     (FSD(I).EQ.OSPVAL) THEN
            CYCLE
          ENDIF
          IP = XP(I)
          JP = YP(I)
          DX = XP(I) - IP
          DY = YP(I) - JP
          IF     (DX.EQ.0.0) THEN
            IPP1 = IP
          ELSEIF (DX.EQ.1.0) THEN
            IPP1 = IP+1
            IP   = IPP1
          ELSEIF (IP.EQ.IDM) THEN
            IPP1 = 1
          ELSE
            IPP1 = IP+1
          ENDIF
          IF     (DY.EQ.0.0) THEN
            JPP1 = JP
          ELSEIF (DY.EQ.1.0) THEN
            JPP1 = JP+1
            JP   = JPP1
          ELSE
            JPP1 = JP+1
          ENDIF
          AXY = (1.0-DX)*(1.0-DY)*A(IP,  JP)   +
     +          (1.0-DX)*     DY *A(IP,  JPP1) +
     +               DX *(1.0-DY)*A(IPP1,JP)   +
     +               DX *     DY *A(IPP1,JPP1)
          TK(K,I) = AXY
          IF     (IDEBUG.EQ.I) THEN
            WRITE(6,*) 'debug: i,k,tk = ',I,K,TK(K,I)
          ENDIF
        ENDDO !i
c
c ---   s
        IF     (ABS(ARTYPE).EQ.1) THEN
          IF     (I_ARCH(17).LT.0) THEN
            READ(11,REC=KREC0+5,IOSTAT=IOS) A
#ifdef ENDIAN_IO
            CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
          ELSEIF (I_ARCH(17).GT.0) THEN
            READ(11,REC=I_ARCH(17),IOSTAT=IOS) A
#ifdef ENDIAN_IO
            CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
          ELSE
            A(:,:) = 0.0
          ENDIF
        ELSE  !mean archive
          READ(11,REC=KREC0+6,IOSTAT=IOS) A
#ifdef ENDIAN_IO
          CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        ENDIF
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'can''t read ',TRIM(CFILEA)
          CALL EXIT(4)
          STOP
        ENDIF
        DO I= 1,N
          IF     (FSD(I).EQ.OSPVAL) THEN
            CYCLE
          ENDIF
          IP = XP(I)
          JP = YP(I)
          DX = XP(I) - IP
          DY = YP(I) - JP
          IF     (DX.EQ.0.0) THEN
            IPP1 = IP
          ELSEIF (DX.EQ.1.0) THEN
            IPP1 = IP+1
            IP   = IPP1
          ELSEIF (IP.EQ.IDM) THEN
            IPP1 = 1
          ELSE
            IPP1 = IP+1
          ENDIF
          IF     (DY.EQ.0.0) THEN
            JPP1 = JP
          ELSEIF (DY.EQ.1.0) THEN
            JPP1 = JP+1
            JP   = JPP1
          ELSE
            JPP1 = JP+1
          ENDIF
          AXY = (1.0-DX)*(1.0-DY)*A(IP,  JP)   +
     +          (1.0-DX)*     DY *A(IP,  JPP1) +
     +               DX *(1.0-DY)*A(IPP1,JP)   +
     +               DX *     DY *A(IPP1,JPP1)
          SK(K,I) = AXY
        ENDDO !i
      ENDDO !k
C
C     INTERPOLATE ALL PROFILES TO Z
C
      M = 1 !start of the 1st profile
      DO I= 1,N
        IF     (FSD(I).EQ.OSPVAL) THEN
          UZ( I) = OSPVAL
          VZ( I) = OSPVAL
          TZ( I) = OSPVAL
          SZ( I) = OSPVAL
          M = M+1  !start of the next z-profile
          CYCLE
        ENDIF
        IF     (I.LT.N) THEN
          IF     (XP(I).EQ.XP(I+1) .AND.
     &            YP(I).EQ.YP(I+1) .AND.
     &            ZZ(I).LE.ZZ(I+1)       ) THEN
C           I AND I+1 FORM PART OF A SINGLE PROFILE
            CYCLE
          ENDIF
        ENDIF
C
C       POINTS M TO I ARE A SINGLE PROFILE
C
        KZ = I-M+1
C
        IF     (BOT.NE.0.0) THEN
C
C         IGNORE LAYERS WITHIN BOT OF THE BOTTOM
C
          DO K= 2,KDM
            IF      (P(K-1,I).GE.P(KDM,I)-BOT) THEN
              UK(K,I) = UK(K-1,I)
              VK(K,I) = VK(K-1,I)
              TK(K,I) = TK(K-1,I)
              SK(K,I) = SK(K-1,I)
            ENDIF
          ENDDO !k
        ENDIF
C
        FLAG = OSPVAL
        IF     (ITYPE.EQ.-2) THEN
          CALL LAYER2Z_PPM(UK(1,I),VK(1,I),TK(1,I),SK(1,I),P(0,I), KDM,
     &                     UZ(  M),VZ(  M),TZ(  M),SZ(  M),ZZ( M),  KZ,
     &                     FLAG)
        ELSEIF (ITYPE.EQ.-1) THEN
          CALL LAYER2Z_PLM(UK(1,I),VK(1,I),TK(1,I),SK(1,I),P(0,I), KDM,
     &                     UZ(  M),VZ(  M),TZ(  M),SZ(  M),ZZ( M),  KZ,
     &                     FLAG)
        ELSEIF (ITYPE.EQ. 0) THEN
          CALL LAYER2Z_PCM(UK(1,I),VK(1,I),TK(1,I),SK(1,I),P(0,I), KDM,
     &                     UZ(  M),VZ(  M),TZ(  M),SZ(  M),ZZ( M),  KZ,
     &                     FLAG)
        ELSEIF (ITYPE.EQ. 1 .AND. IDEBUG.GE.M .AND. IDEBUG.LE.I) THEN
          WRITE(6,*) 'debug: m,i  = ',m,i,kz
          CALL LAYER2Z_LIN_DEBUG(
     &                     UK(1,I),VK(1,I),TK(1,I),SK(1,I),P(0,I), KDM,
     &                     UZ(  M),VZ(  M),TZ(  M),SZ(  M),ZZ( M),  KZ,
     &                     FLAG)
        ELSEIF (ITYPE.EQ. 1) THEN
          CALL LAYER2Z_LIN(UK(1,I),VK(1,I),TK(1,I),SK(1,I),P(0,I), KDM,
     &                     UZ(  M),VZ(  M),TZ(  M),SZ(  M),ZZ( M),  KZ,
     &                     FLAG)
        ELSE
          WRITE(6,"(2a)") 
     +      'Usage: hycom_profile_obs',
     +      ' archv.a list.txt obs.txt [itype [bot]]'
          WRITE(6,*) 'unsupported itype value'
          CALL EXIT(1)
          STOP
        ENDIF
        IF     (PANG(I).NE.0.0) THEN
C
C         ROTATE TO EASTWARDS,NORTHWARDS
C
          DO K= M,I
            IF     (UZ(K).NE.SPVAL) THEN
              DX = UZ(K)
              DY = VZ(K)
              UZ(K) = COS(PANG(I))*DX + SIN(-PANG(I))*DY
              VZ(K) = COS(PANG(I))*DY - SIN(-PANG(I))*DX
            ENDIF
          ENDDO !k
        ENDIF
        IF     (IDEBUG.GE.M .AND. IDEBUG.LE.I) THEN
          WRITE(6,*) 'debug: m,i  = ',m,i
          WRITE(6,*) 'debug: i,tz = ',IDEBUG,TZ(IDEBUG)
        ENDIF
        M = I+1  !start of the next z-profile
      ENDDO !i
C
C     OUTPUT.
C
      OPEN(UNIT=21, FILE=CFILEO, FORM='FORMATTED', STATUS='NEW',
     +         IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',TRIM(CFILEO)
        WRITE(6,*) 'ios   = ',ios
        CALL EXIT(3)
        STOP
      ENDIF
C
      WRITE(21,'(2a)')
     &  '#     LAT     LON       Z',
     &    '       U       V       T       S     SSH'
      DO I= 1,N
        WRITE(21,'(f9.3,f8.3,f8.2,5f8.4)')
     &    PLON(I),PLAT(I),ZZ(I),
     &    UZ(I),VZ(I),TZ(I),SZ(I),FSD(I)
      ENDDO !i
C
      CLOSE(21)
      END
