      PROGRAM PROFILE
      IMPLICIT NONE
C
C  hycom_profile      - Usage:  hycom_profile archv.a ip jp [grid.a]
C  hycom_profile_uv   - Usage:  hycom_profile_uv   archv.a ip jp [grid.a]
C  hycom_profile_all  - Usage:  hycom_profile_all archv.a if il ii jf jl ji prof [grid.a]
C  hycom_profile_list - Usage:  hycom_profile_list archv.a n list.txt prof [grid.a]
C
C                 generates a text profile file at points (ip,jp)
C                 from a HYCOM archive file.
C
C   hycom_profile:      single point, output on stdout
C   hycom_profile_uv:   single point, output on stdout, actual u and v
C   hycom_profile_all:  each point where ip=if:il:ii and jp=jf:jl:ji
C   hycom_profile_list: the n points, one per line, in list.txt
C
C   Append "1" or "2" or "3" to the command name (e.g. hycom_profile3)
C    to include visc or visc,tdif or visc,tdif,sdif in the output.
C   Note that previous versions always included visc,tdif - so use
C    the "2" version to recover the old default behaviour.
C
C   archv.a is assumed to be an HYCOM archive data file, with companion
C   header file archv.b.
C
C   grid.a is a hycom grid file, default regional.grid.a.
C   Note that the corresponding grid.b must also exist.
C
C   For hycom_profile the profile is printed to standard out.  Otherwise,
C   prof is the output filname suffix, the actual filenames are of the
C   form prof_IIII_JJJJ.txt, where IIII is ip and JJJJ is jp as 4 digits.
C   If prof ends with X or XX or XXX or XXXX the X's are replaced with
C   the location within the list, instead of adding _IIII_JJJJ, e.g.
C   if prof is archv.2003_133_BXX the 19th point in the list will be
C   written to the file archv.2003_133_B19.txt and it would be an error
C   for the list to contain 100 or more points.  For hycom_profile_all
C   the list order has i varing fastest from if,jf to il,jl.
C   For hycom_profile_list, if prof ends with the string NAME this is
C   replaced by the string that follows ip,jp on each line of list.txt.
C
C  Except when using hycom_profile_uv, velocity is interpolated to the p-grid.
C  This will return half the correct value at an open nested boundary.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  May 2003 and May 2007.
C
      REAL*4     ONEM,SPVAL
      PARAMETER (ONEM=9806.0, SPVAL=2.0**100)
C
      REAL*4, ALLOCATABLE :: A(:,:)
      INTEGER,ALLOCATABLE :: IP(:),JP(:),IP1(:),JP1(:)
      CHARACTER*240, ALLOCATABLE ::SP(:)
      REAL*4, ALLOCATABLE :: P(:,:),PLOC(:,:),
     +                       PROF(:,:,:),PROU(:,:,:),PROV(:,:,:),
     +                       SURF(:,:),  SURU(:,:),  SURV(:,:),
     +                       UK(:,:),VK(:,:),RK(:,:),
     +                       UB(:),VB(:),SSHS(:),ONETA(:),WFLX(:),
     +                       CICE(:),HICE(:),TICE(:)
      REAL*4              :: PAD(4096)
      REAL*4              :: UBIJ,VBIJ,UKIJ,VKIJ,DIF(3)
C
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      LOGICAL       LSTERIC,LONETA,LWTRFLX,LSEAICE,LUV
      INTEGER       IDM,JDM,KDM,KTR,NSURF,NLAY,NDIF,NTRC,IEXPT,YRFLAG
      INTEGER       NTRC1
      INTEGER       ITYPE,NPAD,NX
      INTEGER       IF,IL,II,JF,JL,JI,N,NM
      INTEGER       ARTYPE,SIGVER,NSTEP,NDIFO
      REAL          THBASE,TIME
      REAL*8        TIME3(3)
      CHARACTER*240 CFILEA,CFILEB,CFILEL,CFILEP,CFILER,CFORMAT
C
      CHARACTER*18 CASN
      INTEGER      I,IU,J,K,KREC,KREC0,L,IOS,NRECL
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
C     READ ARGUMENTS.
C
      CALL GETARG(0,CARG)
      L = LEN_TRIM(CARG)
*     WRITE(6,"(4a)") TRIM(CARG),'"',CARG(L-4:L),'"'
      IF     (CARG(L:L).EQ."3") THEN
        NDIFO = 3
        L = L-1
      ELSEIF (CARG(L:L).EQ."2") THEN
        NDIFO = 2
        L = L-1
      ELSEIF (CARG(L:L).EQ."1") THEN
        NDIFO = 1
        L = L-1
      ELSE
        NDIFO = 0
      ENDIF
      LUV = .FALSE.  !usual case
      IF     (CARG(L-4:L).EQ.'_list') THEN
        ITYPE=3
      ELSEIF (CARG(L-3:L).EQ.'_all') THEN
        ITYPE=2
      ELSEIF (CARG(L-6:L).EQ.'profile') THEN
        ITYPE=1
      ELSEIF (CARG(L-9:L).EQ.'profile_uv') THEN
        ITYPE=1
        LUV  =.TRUE.
      ELSE
        WRITE(6,'(2a)')
     &    'Usage:  ',
     &    'hycom_profile[_uv]; hycom_profile_all; hycom_profile_list'
        CALL EXIT(1)
      ENDIF
C
      NARG = IARGC()
C
      IF     (ITYPE.EQ.1) THEN
        IF     (NARG.EQ.3) THEN
          CALL GETARG(1,CFILEA)
          CALL GETARG(2,CARG)
          READ(CARG,*) IF
          CALL GETARG(3,CARG)
          READ(CARG,*) JF
          CFILER = 'regional.grid.a'
          THBASE = 25.0  !usual value for sigma-theta cases
        ELSEIF (NARG.EQ.4) THEN  ! undocumented, specify THBASE
          CALL GETARG(1,CFILEA)
          CALL GETARG(2,CARG)
          READ(CARG,*) IF
          CALL GETARG(3,CARG)
          READ(CARG,*) JF
          CALL GETARG(4,CFILER)
          THBASE = 25.0  !usual value for sigma-theta cases
        ELSEIF (NARG.EQ.5) THEN  ! undocumented, specify THBASE
          CALL GETARG(1,CFILEA)
          CALL GETARG(2,CARG)
          READ(CARG,*) IF
          CALL GETARG(3,CARG)
          READ(CARG,*) JF
          CALL GETARG(4,CFILER)
          CALL GETARG(5,CARG)
          READ(CARG,*) THBASE
        ELSE
          WRITE(6,*) 'Usage: hycom_profile archv.a ip jp [grid.a]'
          CALL EXIT(1)
        ENDIF
        N = 1
        CFILEP = "stdout"
      ELSEIF (ITYPE.EQ.2) THEN
        IF     (NARG.EQ.8) THEN
          CALL GETARG(1,CFILEA)
          CALL GETARG(2,CARG)
          READ(CARG,*) IF
          CALL GETARG(3,CARG)
          READ(CARG,*) IL
          CALL GETARG(4,CARG)
          READ(CARG,*) II
          CALL GETARG(5,CARG)
          READ(CARG,*) JF
          CALL GETARG(6,CARG)
          READ(CARG,*) JL
          CALL GETARG(7,CARG)
          READ(CARG,*) JI
          CALL GETARG(8,CFILEP)
          CFILER = 'regional.grid.a'
          THBASE = 25.0  !usual value for sigma-theta cases
        ELSEIF (NARG.EQ.9) THEN
          CALL GETARG(1,CFILEA)
          CALL GETARG(2,CARG)
          READ(CARG,*) IF
          CALL GETARG(3,CARG)
          READ(CARG,*) IL
          CALL GETARG(4,CARG)
          READ(CARG,*) II
          CALL GETARG(5,CARG)
          READ(CARG,*) JF
          CALL GETARG(6,CARG)
          READ(CARG,*) JL
          CALL GETARG(7,CARG)
          READ(CARG,*) JI
          CALL GETARG(8,CFILEP)
          CALL GETARG(9,CFILER)
          READ(CARG,*) THBASE
        ELSEIF (NARG.EQ.10) THEN  ! undocumented, specify THBASE
          CALL GETARG(1,CFILEA)
          CALL GETARG(2,CARG)
          READ(CARG,*) IF
          CALL GETARG(3,CARG)
          READ(CARG,*) IL
          CALL GETARG(4,CARG)
          READ(CARG,*) II
          CALL GETARG(5,CARG)
          READ(CARG,*) JF
          CALL GETARG(6,CARG)
          READ(CARG,*) JL
          CALL GETARG(7,CARG)
          READ(CARG,*) JI
          CALL GETARG(8,CFILEP)
          CALL GETARG(9,CFILER)
          CALL GETARG(10,CARG)
          READ(CARG,*) THBASE
        ELSE
          WRITE(6,*) 
     +'Usage: hycom_profile_all archv.a if il ii jf jl ji prof [grid.a]'
          CALL EXIT(1)
        ENDIF
        N = ((IL-IF)/II + 1)*((JL-JF)/JI + 1)
      ELSEIF (ITYPE.EQ.3) THEN
        IF     (NARG.EQ.4) THEN
          CALL GETARG(1,CFILEA)
          CALL GETARG(2,CARG)
          READ(CARG,*) N
          CALL GETARG(3,CFILEL)
          CALL GETARG(4,CFILEP)
          CFILER = 'regional.grid.a'
          THBASE = 25.0  !usual value for sigma-theta cases
        ELSEIF (NARG.EQ.5) THEN  ! undocumented, specify THBASE
          CALL GETARG(1,CFILEA)
          CALL GETARG(2,CARG)
          READ(CARG,*) N
          CALL GETARG(3,CFILEL)
          CALL GETARG(4,CFILEP)
          CALL GETARG(5,CFILER)
          THBASE = 25.0  !usual value for sigma-theta cases
        ELSEIF (NARG.EQ.6) THEN  ! undocumented, specify THBASE
          CALL GETARG(1,CFILEA)
          CALL GETARG(2,CARG)
          READ(CARG,*) N
          CALL GETARG(3,CFILEL)
          CALL GETARG(4,CFILEP)
          CALL GETARG(5,CFILER)
          CALL GETARG(6,CARG)
          READ(CARG,*) THBASE
        ELSE
          WRITE(6,*) 
     +      'Usage: hycom_profile_list archv.a n list.txt prof [grid.a]'
          CALL EXIT(1)
        ENDIF
      ENDIF
C
C     EXTRACT MODEL PARAMETERS FROM ".b" FILE.
C
      CFILEB = CFILEA(1:LEN_TRIM(CFILEA)-1) // 'b'
      CALL READ_BSCW(CFILEB,
     +               IEXPT,YRFLAG,IDM,JDM,KDM,NSURF,NLAY,NDIF,NTRC,
     +               LSTERIC,LONETA,LWTRFLX,LSEAICE,
     +               ARTYPE,SIGVER,THBASE,TIME3,NSTEP)
      TIME = TIME3(1)
C
      IF     (ARTYPE.EQ.3 .OR. ARTYPE.EQ.4) THEN
        NTRC1 = NTRC + 1  !allow for thknss
      ELSE
        NTRC1 = NTRC
      ENDIF
C
      CALL SIG_I(SIGVER)
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
        WRITE(6,*) 'Error in hycom_profile_list: could not allocate ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( PLOC(              2,      N),
     +          SURF(             18,      N),
     +          SURU(             18,      N),
     +          SURV(             18,      N),
     +          PROF(NLAY+NDIF+NTRC1,1:KDM,N),
     +          PROU(NLAY+NDIF+NTRC1,1:KDM,N),
     +          PROV(NLAY+NDIF+NTRC1,1:KDM,N),
     +            RK(                1:KDM,N),
     +             P(                0:KDM,N),
     +            UK(                0:KDM,N),
     +            VK(                0:KDM,N),
     +            UB(                      N),
     +            VB(                      N),
     +          SSHS(                      N),
     +         ONETA(                      N),
     +          WFLX(                      N),
     +          CICE(                      N),
     +          HICE(                      N),
     +          TICE(                      N),
     +            SP(                      N),
     +            IP(                      N),
     +           IP1(                      N),
     +            JP(                      N),
     +           JP1(                      N), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_profile_list: could not allocate ',
     +             '(...,1:n) arrays'
        CALL EXIT(2)
      ENDIF
C
C     GENERATE IP,JP
C
      NX =  0
      IF     (ITYPE.EQ.1) THEN
        IP(1) = IF
        JP(1) = JF
      ELSEIF (ITYPE.EQ.2) THEN
        NM = 0
        DO J= JF,JL,JI
          DO I= IF,IL,II
            NM = NM + 1
            IP(NM) = I
            JP(NM) = J
          ENDDO !i
        ENDDO !j
      ELSEIF (ITYPE.EQ.3) THEN
C
C       DOES OUTPUT FILENAME END WITH "NAME"?
C
        I = LEN_TRIM(CFILEP)
        IF     (CFILEP(I-3:I).EQ.'NAME') THEN
          NX = -1
        ENDIF
C
C       READ PROFILE LOCATIONS
C
        OPEN(UNIT=99, FILE=CFILEL, FORM='FORMATTED', STATUS='OLD',
     +           IOSTAT=IOS)
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'Error: can''t open ',TRIM(CFILEL)
          WRITE(6,*) 'ios   = ',ios
          CALL EXIT(3)
        ENDIF
C
        DO I= 1,N
          IF     (NX.NE.-1) THEN
            READ(99,*,IOSTAT=IOS) IP(I),JP(I)
          ELSE
            READ(99,*,IOSTAT=IOS) IP(I),JP(I),SP(I)
          ENDIF
          IF     (IOS.NE.0) THEN
            WRITE(6,*) 'Error: can''t read ',TRIM(CFILEL)
            WRITE(6,*) 'ios   = ',ios
            WRITE(6,*) 'i,n   = ',i,n
            CALL EXIT(3)
          ENDIF
        ENDDO
      ENDIF !itype
C
      IF     (LUV) THEN  !velocity at ip,jp
        DO I= 1,N
          IP1(I) = IP(I)
          JP1(I) = JP(I)
        ENDDO
      ELSE !usual case, interpolated velocity
        DO I= 1,N
          IP1(I) = MOD(IP(I),  IDM) + 1
          JP1(I) = MIN(JP(I)+1,JDM)
        ENDDO
      ENDIF
C
      DO I= 1,N
        IP1(I) = MOD(IP(I),  IDM) + 1
        JP1(I) = MIN(JP(I)+1,JDM)
      ENDDO
C
C     CHECK FOR THE XXX SYNTAX.
C
      IF     (NX.NE.-1) THEN
        I = LEN_TRIM(CFILEP)
        NX = I
        DO K= I,1,-1
          IF     (CFILEP(K:K).NE.'X') THEN
            NX = I-K
            EXIT
          ENDIF
        ENDDO
      ENDIF
      IF     (NX.EQ.1) THEN
        IF     (N.GT.   9) THEN
          WRITE(6,*) 'Error: 1 digit number but more than 9 points'
          CALL EXIT(4)
        ENDIF
      ELSEIF (NX.EQ.2) THEN
        IF     (N.GT.  99) THEN
          WRITE(6,*) 'Error: 2 digit number but more than 99 points'
          CALL EXIT(4)
        ENDIF
      ELSEIF (NX.EQ.3) THEN
        IF     (N.GT. 999) THEN
          WRITE(6,*) 'Error: 3 digit number but more than 999 points'
          CALL EXIT(4)
        ENDIF
      ELSEIF (NX.EQ.4) THEN
        IF     (N.GT.9999) THEN
          WRITE(6,*) 'Error: 1 digit number but more than 9999 points'
          CALL EXIT(4)
        ENDIF
      ELSEIF (NX.GT.4) THEN
        IF     (N.GT.9999) THEN
          WRITE(6,*) 'Error: output filename ends in XXXXX'
          CALL EXIT(4)
        ENDIF
      ENDIF
C
C     OPEN ".a" FILES.
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
      OPEN(UNIT=11, FILE=CFILER, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',TRIM(CFILER)
        WRITE(6,*) 'ios   = ',ios
        WRITE(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
C     LON,LAT FIELDS.
C
      DO KREC= 1,2
        READ(11,REC=KREC,IOSTAT=IOS) A
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'can''t read ',TRIM(CFILEA)
          CALL EXIT(4)
        ENDIF
*       write(6,*) 'REC,A.MIN = ',KREC,MINVAL(A),A(IF,JF)
        DO I= 1,N
          PLOC(KREC,I) = A(IP(I),JP(I))
        ENDDO
      ENDDO !krec
*     write(6,*) 'PLOC.1 = ',PLOC(1:2,1)
      CLOSE(11)
C
      OPEN(UNIT=11, FILE=CFILEA, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',TRIM(CFILEA)
        WRITE(6,*) 'ios   = ',ios
        WRITE(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
C     SURFACE FIELDS.
C
      DO KREC= 1,NSURF
        READ(11,REC=KREC,IOSTAT=IOS) A
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'can''t read ',TRIM(CFILEA)
          CALL EXIT(4)
        ENDIF
*       write(6,*) 'REC,A.MIN = ',KREC,MINVAL(A),A(IF,JF)
        DO I= 1,N
          SURF(KREC,I) = A(IP( I),JP( I))
          SURU(KREC,I) = A(IP1(I),JP( I))
          SURV(KREC,I) = A(IP( I),JP1(I))
          IF     (SURF(KREC,I).EQ.SPVAL) THEN
*           write(6,*) 'i,kr,s.p =',I,KREC,SURF(KREC,I)
            SURF(KREC,I) = 0.0  !assume this is a velocity field
          ENDIF
          IF     (SURU(KREC,I).EQ.SPVAL) THEN
*           write(6,*) 'i,kr,s.u =',I,KREC,SURU(KREC,I)
            SURU(KREC,I) = 0.0  !assume this is a velocity field
          ENDIF
          IF     (SURV(KREC,I).EQ.SPVAL) THEN
*           write(6,*) 'i,kr,s.v =',I,KREC,SURV(KREC,I)
            SURV(KREC,I) = 0.0  !assume this is a velocity field
          ENDIF
        ENDDO
      ENDDO !krec
*     write(6,*) 'NSURF  = ',NSURF
*     write(6,*) 'SURF.1 = ',SURF(1:NSURF,1)
*     write(6,*) 'SURU.1 = ',SURU(1:NSURF,1)
*     write(6,*) 'SURV.1 = ',SURV(1:NSURF,1)
C
      IF     (ARTYPE.GT.0) THEN
        IF     (ARTYPE.EQ.1) THEN
          UB(:) = 0.5*(SURF(NSURF-1,:)+SURU(NSURF-1,:))
          VB(:) = 0.5*(SURF(NSURF,  :)+SURV(NSURF,  :))
        ELSE  !mean archive
          UB(:) = 0.5*(SURF(NSURF-2,:)+SURU(NSURF-2,:))
          VB(:) = 0.5*(SURF(NSURF-1,:)+SURV(NSURF-1,:))
        ENDIF
      ELSE !p-grid archive
        IF     (ARTYPE.EQ.-1) THEN
          UB(:) = SURF(NSURF-1,:)
          VB(:) = SURF(NSURF,  :)
        ELSE  !mean archive
          UB(:) = SURF(NSURF-2,:)
          VB(:) = SURF(NSURF-1,:)
        ENDIF
      ENDIF
C
      P(1,:) = 0.0
      DO K= 1,KDM
        KREC0 = NSURF+(NLAY+NDIF+NTRC)*(K-1)
        DO KREC= KREC0+1,KREC0+NLAY+NDIF+NTRC
          READ(11,REC=KREC,IOSTAT=IOS) A
#ifdef ENDIAN_IO
          CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
          IF     (IOS.NE.0) THEN
            WRITE(6,*) 'can''t read ',TRIM(CFILEA)
            CALL EXIT(4)
          ENDIF
          DO I= 1,N
            PROF(KREC-KREC0,K,I) = A(IP( I),JP( I))
            PROU(KREC-KREC0,K,I) = A(IP1(I),JP( I))
            PROV(KREC-KREC0,K,I) = A(IP( I),JP1(I))
            IF     (PROF(KREC-KREC0,K,I).EQ.SPVAL) THEN
              PROF(KREC-KREC0,K,I) = 0.0  !assume this is a velocity field
            ENDIF
            IF     (PROU(KREC-KREC0,K,I).EQ.SPVAL) THEN
              PROU(KREC-KREC0,K,I) = 0.0  !assume this is a velocity field
            ENDIF
            IF     (PROV(KREC-KREC0,K,I).EQ.SPVAL) THEN
              PROV(KREC-KREC0,K,I) = 0.0  !assume this is a velocity field
            ENDIF
          ENDDO
        ENDDO
        IF     (ARTYPE.GT.0) THEN
          IF     (ARTYPE.EQ.1) THEN
            UK(K, :) =  UB(:) + 0.5*(PROF(1,K,:)+PROU(1,K,:))
            VK(K, :) =  VB(:) + 0.5*(PROF(2,K,:)+PROV(2,K,:))
            P(K+1,:) = P(K,:) + PROF(3,K,:)
          ELSEIF (ARTYPE.EQ.2) THEN  !mean archive
            UK(K, :) = 0.5*(PROF(1,K,:)+PROU(1,K,:))  !already full layer velocity
            VK(K, :) = 0.5*(PROF(2,K,:)+PROV(2,K,:))  !already full layer velocity
            P(K+1,:) = P(K,:) + PROF(NLAY-3,K,:)
          ELSE  !std. or diff archive
            UK(K, :) = 0.5*(PROF(1,K,:)+PROU(1,K,:))  !already full layer velocity
            VK(K, :) = 0.5*(PROF(2,K,:)+PROV(2,K,:))  !already full layer velocity
            P(K+1,:)         = P(K,:) + PROF(NLAY-4,K,:)     !mnthknss
            PROF(NLAY+NDIF+NTRC1,K,:) = PROF(NLAY-3,K,:)/ONEM  !thknss
          END IF
        ELSE !p-grid archive
          IF     (ARTYPE.EQ.-1) THEN
            UK(K, :) =  UB(:) + PROF(1,K,:)
            VK(K, :) =  VB(:) + PROF(2,K,:)
            P(K+1,:) = P(K,:) + PROF(3,K,:)
          ELSE  !mean archive
            UK(K, :) =          PROF(1,K,:)  !already full layer velocity
            VK(K, :) =          PROF(2,K,:)  !already full layer velocity
            P(K+1,:) = P(K,:) + PROF(NLAY-3,K,:)
          ENDIF
        ENDIF
C
        IF     (SIGVER.GT.0) THEN
          DO I= 1,N
            CALL SIG_P(PROF(NLAY-1,K,I),  !T
     &                 PROF(NLAY,  K,I),  !S
     &                          RK(K,I))
          ENDDO !I
        ENDIF
      ENDDO !K
C
      IF     (LSEAICE) THEN
        IF     (ABS(ARTYPE).EQ.1) THEN
          CICE(:) = SURF(NSURF-4,:)
          HICE(:) = SURF(NSURF-3,:)
          TICE(:) = SURF(NSURF-2,:)
        ELSE  !mean archive
          CICE(:) = SURF(NSURF-5,:)
          HICE(:) = SURF(NSURF-4,:)
          TICE(:) = SURF(NSURF-3,:)
        ENDIF
      ENDIF !lseaice
C
      IF     (SIGVER.NE.0) THEN
        IF     (LSTERIC .AND. LONETA .AND. LWTRFLX) THEN
C         remove steric SSH,one+eta,water flux
          SSHS(  :) = SURF( 3,:) !steric
          ONETA( :) = SURF( 4,:) !oneta 
          WFLX(  :) = SURF( 6,:) !wtrflx
          SURF(3,:) = SURF( 5,:) !surflx
          SURF(4,:) = SURF( 7,:) !salflx
          SURF(5,:) = SURF( 8,:) !bl_dpth
          SURF(6,:) = SURF( 9,:) !mix_dpth
          SURF(7,:) = SURF(10,:) !tmix
          SURF(8,:) = SURF(11,:) !smix
          SURF(9,:) = SURF(12,:) !thmix
        ELSEIF (LONETA .AND. LWTRFLX) THEN
C         remove one+eta,water flux
          SSHS(  :) = SURF( 1,:) !steric = montg1
          ONETA( :) = SURF( 3,:) !oneta 
          WFLX(  :) = SURF( 5,:) !wtrflx
          SURF(3,:) = SURF( 4,:) !surflx
          SURF(4,:) = SURF( 6,:) !salflx
          SURF(5,:) = SURF( 7,:) !bl_dpth
          SURF(6,:) = SURF( 8,:) !mix_dpth
          SURF(7,:) = SURF( 9,:) !tmix
          SURF(8,:) = SURF(10,:) !smix
          SURF(9,:) = SURF(11,:) !thmix
        ELSEIF (LSTERIC .AND. LONETA) THEN
C         remove steric SSH,one+eta
          SSHS(  :) = SURF( 3,:) !steric
          ONETA( :) = SURF( 4,:) !oneta 
          SURF(3,:) = SURF( 5,:) !surflx
          SURF(4,:) = SURF( 6,:) !salflx
          SURF(5,:) = SURF( 7,:) !bl_dpth
          SURF(6,:) = SURF( 8,:) !mix_dpth
          SURF(7,:) = SURF( 9,:) !tmix
          SURF(8,:) = SURF(10,:) !smix
          SURF(9,:) = SURF(11,:) !thmix
          WFLX(  :) = -SURF( 4,I)*1.E-3/PROF(NLAY,1,:)
        ELSEIF (LSTERIC) THEN
C         remove steric SSH
          SSHS(  :) = SURF( 3,:) !steric
          SURF(3,:) = SURF( 4,:) !surflx
          SURF(4,:) = SURF( 5,:) !salflx
          SURF(5,:) = SURF( 6,:) !bl_dpth
          SURF(6,:) = SURF( 7,:) !mix_dpth
          SURF(7,:) = SURF( 8,:) !tmix
          SURF(8,:) = SURF( 9,:) !smix
          SURF(9,:) = SURF(10,:) !thmix
          ONETA( :) = 1.0
          WFLX(  :) = -SURF( 4,I)*1.E-3/PROF(NLAY,1,:)
        ELSE
          SSHS(  :) = SURF(1,:) !steric = montg1
          ONETA( :) = 1.0
          WFLX(  :) = -SURF( 4,I)*1.E-3/PROF(NLAY,1,:)
        ENDIF  !steric/oneta/wtrflx
        SURF(7,:) = PROF(NLAY-1,1,:)         !T
        SURF(8,:) = PROF(NLAY,  1,:)         !S
        SURF(9,:) =          RK(1,:)-THBASE  !TH
          UK(0,:) =          UK(1,:)
          VK(0,:) =          VK(1,:)
      ELSEIF (ARTYPE.EQ. 1) THEN
        IF     (LSTERIC .AND. LONETA .AND. LWTRFLX) THEN
C         remove steric SSH,one+eta,water flux
          SSHS(  :) = SURF( 3,:) !steric
          ONETA( :) = SURF( 4,:) !oneta 
          WFLX(  :) = SURF( 6,:) !wtrflx
          SURF(3,:) = SURF( 5,:) !surflx
          SURF(4,:) = SURF( 7,:) !salflx
          SURF(5,:) = SURF( 8,:) !bl_dpth
          SURF(6,:) = SURF( 9,:) !mix_dpth
          SURF(7,:) = SURF(10,:) !tmix
          SURF(8,:) = SURF(11,:) !smix
          SURF(9,:) = SURF(12,:) !thmix
          UK(0,  :) = 0.5*(SURF(13,:)+SURU(13,:))+UB(:)
          VK(0,  :) = 0.5*(SURF(14,:)+SURV(14,:))+VB(:)
        ELSEIF (LONETA .AND. LWTRFLX) THEN
C         remove one+eta,water flux
          SSHS(  :) = SURF( 1,:) !steric = montg1
          ONETA( :) = SURF( 3,:) !oneta 
          WFLX(  :) = SURF( 5,:) !wtrflx
          SURF(3,:) = SURF( 4,:) !surflx
          SURF(4,:) = SURF( 6,:) !salflx
          SURF(5,:) = SURF( 7,:) !bl_dpth
          SURF(6,:) = SURF( 8,:) !mix_dpth
          SURF(7,:) = SURF( 9,:) !tmix
          SURF(8,:) = SURF(10,:) !smix
          SURF(9,:) = SURF(11,:) !thmix
          UK(0,  :) = 0.5*(SURF(12,:)+SURU(12,:))+UB(:)
          VK(0,  :) = 0.5*(SURF(13,:)+SURV(13,:))+VB(:)
        ELSEIF (LSTERIC .AND. LONETA) THEN
C         remove steric SSH,one+eta
          SSHS(  :) = SURF( 3,:) !steric
          ONETA( :) = SURF( 4,:) !oneta 
          SURF(3,:) = SURF( 5,:) !surflx
          SURF(4,:) = SURF( 6,:) !salflx
          SURF(5,:) = SURF( 7,:) !bl_dpth
          SURF(6,:) = SURF( 8,:) !mix_dpth
          SURF(7,:) = SURF( 9,:) !tmix
          SURF(8,:) = SURF(10,:) !smix
          SURF(9,:) = SURF(11,:) !thmix
          WFLX(  :) = -SURF( 4,I)*1.E-3/PROF(NLAY,1,:)
          UK(0,  :) = 0.5*(SURF(12,:)+SURU(12,:))+UB(:)
          VK(0,  :) = 0.5*(SURF(13,:)+SURV(13,:))+VB(:)
        ELSEIF (LSTERIC) THEN
C         remove steric SSH
          SSHS(  :) = SURF( 3,:) !steric
          SURF(3,:) = SURF( 4,:) !surflx
          SURF(4,:) = SURF( 5,:) !salflx
          SURF(5,:) = SURF( 6,:) !bl_dpth
          SURF(6,:) = SURF( 7,:) !mix_dpth
          SURF(7,:) = SURF( 8,:) !tmix
          SURF(8,:) = SURF( 9,:) !smix
          SURF(9,:) = SURF(10,:) !thmix
          ONETA( :) = 1.0
          WFLX(  :) = -SURF( 4,I)*1.E-3/PROF(NLAY,1,:)
          UK(0,  :) = 0.5*(SURF(11,:)+SURU(11,:))+UB(:)
          VK(0,  :) = 0.5*(SURF(12,:)+SURV(12,:))+VB(:)
        ELSE
          SSHS(  :) = SURF(1,:) !steric = montg1
          ONETA( :) = 1.0
          WFLX(  :) = -SURF( 4,I)*1.E-3/PROF(NLAY,1,:)
          UK(0,  :) = 0.5*(SURF(10,:)+SURU(10,:))+UB(:)
          VK(0,  :) = 0.5*(SURF(11,:)+SURV(11,:))+VB(:)
        ENDIF  !steric/oneta/wtrflx
      ELSEIF (ARTYPE.EQ.-1) THEN  !p-grid archive
          UK(0,:) = SURF(10,:)+UB(:)
          VK(0,:) = SURF(11,:)+VB(:)
      ELSE  !mean archive
        IF     (LSTERIC .AND. LONETA .AND. LWTRFLX) THEN
C         remove steric SSH,one+eta,water flux
          SSHS(  :) = SURF( 3,:) !steric
          ONETA( :) = SURF( 4,:) !oneta 
          WFLX(  :) = SURF( 6,:) !wtrflx
          SURF(3,:) = SURF( 5,:) !surflx
          SURF(4,:) = SURF( 7,:) !salflx
          SURF(5,:) = SURF( 8,:) !bl_dpth
          SURF(6,:) = SURF( 9,:) !mix_dpth
          SURF(7,:) = SURF(10,:) !tmix
          SURF(8,:) = SURF(11,:) !smix
          SURF(9,:) = SURF(12,:) !thmix
          IF (ARTYPE.GT.0) THEN
            UK(0,:) = 0.5*(SURF(13,:)+SURU(13,:))  !already full m.l. velocity
            VK(0,:) = 0.5*(SURF(14,:)+SURV(14,:))  !already full m.l. velocity
          ELSE !p-grid archive
            UK(0,:) = SURF(13,:)  !already full m.l. velocity
            VK(0,:) = SURF(14,:)  !already full m.l. velocity
          ENDIF
        ELSEIF (LONETA .AND. LWTRFLX) THEN
C         remove one+eta,water flux
          SSHS(  :) = SURF( 1,:) !steric = montg1
          ONETA( :) = SURF( 3,:) !oneta 
          WFLX(  :) = SURF( 5,:) !wtrflx
          SURF(3,:) = SURF( 4,:) !surflx
          SURF(4,:) = SURF( 6,:) !salflx
          SURF(5,:) = SURF( 7,:) !bl_dpth
          SURF(6,:) = SURF( 8,:) !mix_dpth
          SURF(7,:) = SURF( 9,:) !tmix
          SURF(8,:) = SURF(10,:) !smix
          SURF(9,:) = SURF(11,:) !thmix
          IF (ARTYPE.GT.0) THEN
            UK(0,:) = 0.5*(SURF(12,:)+SURU(12,:))  !already full m.l. velocity
            VK(0,:) = 0.5*(SURF(13,:)+SURV(13,:))  !already full m.l. velocity
          ELSE !p-grid archive
            UK(0,:) = SURF(12,:)  !already full m.l. velocity
            VK(0,:) = SURF(13,:)  !already full m.l. velocity
          ENDIF
        ELSEIF (LSTERIC .AND. LONETA) THEN
C         remove steric SSH,one+eta
          SSHS(  :) = SURF( 3,:) !steric
          ONETA( :) = SURF( 4,:) !oneta 
          SURF(3,:) = SURF( 5,:) !surflx
          SURF(4,:) = SURF( 6,:) !salflx
          SURF(5,:) = SURF( 7,:) !bl_dpth
          SURF(6,:) = SURF( 8,:) !mix_dpth
          SURF(7,:) = SURF( 9,:) !tmix
          SURF(8,:) = SURF(10,:) !smix
          SURF(9,:) = SURF(11,:) !thmix
          WFLX(  :) = -SURF( 4,I)*1.E-3/PROF(NLAY,1,:)
          IF (ARTYPE.GT.0) THEN
            UK(0,:) = 0.5*(SURF(12,:)+SURU(12,:))  !already full m.l. velocity
            VK(0,:) = 0.5*(SURF(13,:)+SURV(13,:))  !already full m.l. velocity
          ELSE !p-grid archive
            UK(0,:) = SURF(12,:)  !already full m.l. velocity
            VK(0,:) = SURF(13,:)  !already full m.l. velocity
          ENDIF
        ELSEIF (LSTERIC) THEN
C         remove steric SSH
          SSHS(  :) = SURF( 3,:) !steric
          SURF(3,:) = SURF( 4,:) !surflx
          SURF(4,:) = SURF( 5,:) !salflx
          SURF(5,:) = SURF( 6,:) !bl_dpth
          SURF(6,:) = SURF( 7,:) !mix_dpth
          SURF(7,:) = SURF( 8,:) !tmix
          SURF(8,:) = SURF( 9,:) !smix
          SURF(9,:) = SURF(10,:) !thmix
          ONETA( :) = 1.0
          WFLX(  :) = -SURF( 4,I)*1.E-3/PROF(NLAY,1,:)
          IF (ARTYPE.GT.0) THEN
            UK(0,:) = 0.5*(SURF(11,:)+SURU(11,:))  !already full m.l. velocity
            VK(0,:) = 0.5*(SURF(12,:)+SURV(12,:))  !already full m.l. velocity
          ELSE !p-grid archive
            UK(0,:) = SURF(11,:)  !already full m.l. velocity
            VK(0,:) = SURF(12,:)  !already full m.l. velocity
          ENDIF
        ELSE
          SSHS(  :) = SURF(1,:) !steric = montg1
          ONETA( :) = 1.0
          WFLX(  :) = -SURF( 4,I)*1.E-3/PROF(NLAY,1,:)
          IF (ARTYPE.GT.0) THEN
            UK(0,:) = 0.5*(SURF(10,:)+SURU(10,:))  !already full m.l. velocity
            VK(0,:) = 0.5*(SURF(11,:)+SURV(11,:))  !already full m.l. velocity
          ELSE !p-grid archive
            UK(0,:) = SURF(10,:)  !already full m.l. velocity
            VK(0,:) = SURF(11,:)  !already full m.l. velocity
          ENDIF
        ENDIF  !steric/oneta/wtrflx
      ENDIF  !ARTYPE
C
C     LOOP THROUGH OUTPUT PROFILE FILES.
C
      IF     (ITYPE.EQ.1) THEN
        IU =  6
      ELSE
        IU = 21
      ENDIF
      DO I= 1,N
          IF     (IU.NE.6) THEN
            IF     (NX.EQ.-1) THEN
C   from profNAME, form profA.txt, where A is sp.
              WRITE(CFILEB,'(A,A,".txt")')
     +          CFILEP(1:LEN_TRIM(CFILEP)-4),TRIM(SP(I))
            ELSEIF (NX.EQ.0) THEN
C   form prof_IIII_JJJJ.txt, where IIII is ip and JJJJ is jp as 4 digits.
              WRITE(CFILEB,'(A,"_",I4.4,"_",I4.4,".txt")')
     +          TRIM(CFILEP),IP(I),JP(I)
            ELSEIF (NX.EQ.1) THEN
C   from profX, form profI.txt, where I is i as 1 digit.
              WRITE(CFILEB,'(A,I1.1,".txt")')
     +          CFILEP(1:LEN_TRIM(CFILEP)-1),I
            ELSEIF (NX.EQ.2) THEN
C   from profX, form profII.txt, where II is i as 2 digits.
              WRITE(CFILEB,'(A,I2.2,".txt")')
     +          CFILEP(1:LEN_TRIM(CFILEP)-2),I
            ELSEIF (NX.EQ.3) THEN
C   from profX, form profIII.txt, where III is i as 3 digits.
              WRITE(CFILEB,'(A,I3.3,".txt")')
     +          CFILEP(1:LEN_TRIM(CFILEP)-3),I
            ELSEIF (NX.EQ.4) THEN
C   from profX, form profIIII.txt, where IIII is i as 4 digits.
              WRITE(CFILEB,'(A,I4.4,".txt")')
     +          CFILEP(1:LEN_TRIM(CFILEP)-4),I
            ENDIF
            IF     (SURF(2,I).GT.1.E20) THEN
              WRITE(6,'(2a)') 'Over Land : ',TRIM(CFILEB)
              CYCLE  !land point
            ENDIF
            OPEN(UNIT=IU, FILE=CFILEB, FORM='FORMATTED', STATUS='NEW',
     +           IOSTAT=IOS)
            IF     (IOS.NE.0) THEN
              WRITE(6,'(2a)') 'Can''t open: ',TRIM(CFILEB)
              CYCLE
*           ELSE
*             WRITE(6,'(2a)') 'Writing to: ',TRIM(CFILEB)
            ENDIF
          ELSE
            IF     (SURF(2,I).GT.1.E20) THEN
              WRITE(6,'(a)') 'Over Land'
              CYCLE  !land point
            ENDIF
          ENDIF !IU.ne.6
C
          UBIJ = MAX( -9.99, MIN( 9.99, UB(I) ) )
          VBIJ = MAX( -9.99, MIN( 9.99, VB(I) ) )
          UKIJ = MAX( -9.99, MIN( 9.99, UK(0,I) ) )
          VKIJ = MAX( -9.99, MIN( 9.99, VK(0,I) ) )
C
          K = LEN_TRIM(CFILEA)
          IF     (CFILEA(K-8:K-8).EQ.'_' .AND.
     &            CFILEA(K-4:K-4).EQ.'_'      ) THEN !YYYY_JJJ_HH
            IF     (LUV) THEN  !ILOC,JLOC indicate u&v at location
              WRITE(IU,'(3a/a,6i7,f9.3,f8.3,i7,a5,a4,a3)')
     &          '##   expt    idm    jdm    kdm',
     &          '   ILOC   JLOC   lonloc  latloc',
     &          ' yrflag year day hr',
     &          '##',IEXPT,   IDM,   JDM,   KDM,
     &          IP(I),JP(I),PLOC(1,I),PLOC(2,I),
     &          YRFLAG,CFILEA(K-12:K-9),CFILEA(K-7:K-5),CFILEA(K-3:K-2)
            ELSE !usual case
              WRITE(IU,'(3a/a,6i7,f9.3,f8.3,i7,a5,a4,a3)')
     &            '##   expt    idm    jdm    kdm',
     &          '   iloc   jloc   lonloc  latloc',
     &          ' yrflag year day hr',
     &          '##',IEXPT,   IDM,   JDM,   KDM,
     &          IP(I),JP(I),PLOC(1,I),PLOC(2,I),
     &          YRFLAG,CFILEA(K-12:K-9),CFILEA(K-7:K-5),CFILEA(K-3:K-2)
            ENDIF !luv:else
          ELSE
           IF     (LUV) THEN  !ILOC,JLOC indicate u&v at location
              WRITE(IU,'(3a/a,6i7,2f8.1,i7)')
     &          '##   expt    idm    jdm    kdm',
     &          '   ILOC   JLOC  lonloc  latloc',
     &          ' yrflag',
     &          '##',IEXPT,   IDM,   JDM,   KDM,
     &          IP(I),JP(I),PLOC(1,I),PLOC(2,I),
     &          YRFLAG
            ELSE !usual case
              WRITE(IU,'(3a/a,6i7,2f8.1,i7)')
     &          '##   expt    idm    jdm    kdm',
     &          '   iloc   jloc  lonloc  latloc',
     &          ' yrflag',
     &          '##',IEXPT,   IDM,   JDM,   KDM,
     &          IP(I),JP(I),PLOC(1,I),PLOC(2,I),
     &          YRFLAG
            ENDIF !luv:else
          ENDIF
          WRITE(IU,
     & '(4a/a,f11.4,f8.2,f8.1,2f9.3,3f8.4,6f8.2,f8.4,f12.3)')
     &      '## model-day  srfhgt  surflx',
     &      '     dpbl   dpmixl    tmix    smix   thmix',
     &      '    umix    vmix   ubavg   vbavg  steric  nsterc',
     &      '   oneta         E-P',
     &      '#',TIME,
     &      SURF( 2,I)/(9.806*0.01),  ! cm
     &      SURF( 3,I),               ! W/m**2
     &      SURF( 5,I)/ONEM,          ! m
     &      SURF( 6,I)/ONEM,          ! m
     &      SURF( 7,I),               ! degC
     &      SURF( 8,I),               ! psu
     &      SURF( 9,I)+THBASE,        ! SigmaT
     &      UKIJ*100.0,               ! cm/s
     &      VKIJ*100.0,               ! cm/s
     &      UBIJ*100.0,               ! cm/s
     &      VBIJ*100.0,               ! cm/s
     &      SSHS(   I)/(9.806*0.01),  ! cm
     &      SURF( 2,I)/(9.806*0.01)-
     &      SSHS(   I)/(9.806*0.01),  ! cm
     &     ONETA(   I),               ! unitless
     &     -WFLX(   I)*86400.0*1.E3   ! mm/day
C
          IF     (LSEAICE) THEN
            WRITE(IU,'(2a / a,f11.4, 3f8.2)')
     &        '## model-day',
     &        '  covice  thkice  temice',
     &        '#',TIME,
     &          CICE(I)*100.0,        ! %
     &          HICE(I),              ! m
     &          TICE(I)               ! degC
          ENDIF !lseaice
C
C         PER LAYER FIELDS.
C
          IF     (NTRC1.EQ.NTRC) THEN  !the usual case
            IF     (NTRC.EQ.0) THEN
              WRITE(CFORMAT,'(a)')
     &          '(3a)'
            ELSE
              WRITE(CFORMAT,'(a,i2,a)')
     &          '(3a,',NTRC,'(a,i2.2))'
            ENDIF
            IF     (NDIFO.EQ.0) THEN
              WRITE(IU,CFORMAT)
     &            '#  k',
     &            '    utot    vtot  p.temp    saln  p.dens',
     &            '    thkns      dpth',
     &            ('  trcr',KTR, KTR=1,NTRC)
            ELSEIF (NDIFO.EQ.1) THEN
              WRITE(IU,CFORMAT)
     &            '#  k',
     &            '    utot    vtot  p.temp    saln  p.dens',
     &            '    thkns      dpth  viscty',
     &            ('  trcr',KTR, KTR=1,NTRC)
            ELSEIF (NDIFO.EQ.2) THEN
              WRITE(IU,CFORMAT)
     &            '#  k',
     &            '    utot    vtot  p.temp    saln  p.dens',
     &            '    thkns      dpth  viscty  t-diff',
     &            ('  trcr',KTR, KTR=1,NTRC)
            ELSEIF (NDIFO.EQ.3) THEN
              WRITE(IU,CFORMAT)
     &            '#  k',
     &            '    utot    vtot  p.temp    saln  p.dens',
     &            '    thkns      dpth  viscty  t-diff  s-diff',
     &            ('  trcr',KTR, KTR=1,NTRC)
            ENDIF
C
            IF     (NTRC.EQ.0) THEN
              IF     (NDIFO.EQ.0) THEN
                WRITE(CFORMAT,'(a)')
     &          '(i4,2f8.2,3f8.4,f9.3,f10.3,a)'  !last a not used
              ELSE
                WRITE(CFORMAT,'(a,i1,a)')
     &          '(i4,2f8.2,3f8.4,f9.3,f10.3,',NDIFO,'f8.2)'
              ENDIF !ndifo
            ELSE
              IF     (NDIFO.EQ.0) THEN
                WRITE(CFORMAT,'(a,i2,a)')
     &        '(i4,2f8.2,3f8.4,f9.3,f10.3,',NTRC,'f8.4,a)'  !last a not used
              ELSE
                WRITE(CFORMAT,'(a,i1,a,i2,a)')
     &        '(i4,2f8.2,3f8.4,f9.3,f10.3,',NDIFO,'f8.2,',NTRC,'f8.4)'
              ENDIF !ndifo
            ENDIF !ntrc
C
            DIF(1:3) = 0.0
            DO K= 1,KDM
*             write(IU,*) 'PROF = ',PROF
              UKIJ = MAX( -9.99, MIN( 9.99, UK(K,I) ) )
              VKIJ = MAX( -9.99, MIN( 9.99, VK(K,I) ) )
              IF     (NDIF.EQ.3) THEN
                DIF(1) = MIN(9999.99, PROF(NLAY+NTRC+1,K,I)*1.E4)
                DIF(2) = MIN(9999.99, PROF(NLAY+NTRC+2,K,I)*1.E4)
                DIF(3) = MIN(9999.99, PROF(NLAY+NTRC+3,K,I)*1.E4)
              ENDIF
              IF     (SIGVER.EQ.0) THEN
                WRITE(IU,CFORMAT)
     &            K,
     &            UKIJ*100.0,                          !cm/s
     &            VKIJ*100.0,                          !cm/s
     &            PROF(NLAY-2,K,I),                    !degC
     &            PROF(NLAY-1,K,I),                    !psu
     &            PROF(NLAY,  K,I)  +THBASE,           !SigmaT
     &            (P(K+1,I)-P(K,I))    /ONEM,          !m
     &            (P(K+1,I)+P(K,I))*0.5/ONEM,          !m
     &            (DIF(KTR),KTR=1,NDIFO),              !cm^2/s
     &            (PROF(NLAY+KTR,K,I),KTR=1,NTRC)      !0-99
              ELSE
                WRITE(IU,CFORMAT)
     &            K,
     &            UKIJ*100.0,                          !cm/s
     &            VKIJ*100.0,                          !cm/s
     &            PROF(NLAY-1,K,I),                    !degC
     &            PROF(NLAY,  K,I),                    !psu
     &                     RK(K,I),                    !SigmaT
     &            PROF(NLAY-2,K,I)/ONEM,               !m
     &            (P(K+1,I)+P(K,I))*0.5/ONEM,          !m
     &            (DIF(KTR),KTR=1,NDIFO),              !cm^2/s
     &            (PROF(NLAY+KTR,K,I),KTR=1,NTRC)      !0-99
              ENDIF
            ENDDO !k
          ELSE  !NTRC1.NE.NTRC
            IF     (NTRC.EQ.0) THEN
              WRITE(CFORMAT,'(a)')
     &          '(4a)'
            ELSE
              WRITE(CFORMAT,'(a,i2,a)')
     &          '(4a,',NTRC,'a)'
            ENDIF
            IF     (NDIFO.EQ.0) THEN
              WRITE(IU,CFORMAT)
     &            '#  k',
     &            '    utot    vtot  p.temp    saln  p.dens',
     &            '    thkns      dpth',
     &            ('  tracer',KTR=1,NTRC),
     &             '    thknss'
            ELSEIF (NDIFO.EQ.1) THEN
              WRITE(IU,CFORMAT)
     &            '#  k',
     &            '    utot    vtot  p.temp    saln  p.dens',
     &            '    thkns      dpth  viscty',
     &            ('  tracer',KTR=1,NTRC),
     &             '    thknss'
            ELSEIF (NDIFO.EQ.2) THEN
              WRITE(IU,CFORMAT)
     &            '#  k',
     &            '    utot    vtot  p.temp    saln  p.dens',
     &            '    thkns      dpth  viscty  t-diff',
     &            ('  tracer',KTR=1,NTRC),
     &             '    thknss'
            ELSEIF (NDIFO.EQ.3) THEN
              WRITE(IU,CFORMAT)
     &            '#  k',
     &            '    utot    vtot  p.temp    saln  p.dens',
     &            '    thkns      dpth  viscty  t-diff  s-diff',
     &            ('  tracer',KTR=1,NTRC),
     &             '    thknss'
            ENDIF
C
            IF     (NTRC.EQ.0) THEN
              IF     (NDIFO.EQ.0) THEN
                WRITE(CFORMAT,'(a)')
     &          '(i4,2f8.2,3f8.4,f9.3,f10.3,f10.3)'
              ELSE
                WRITE(CFORMAT,'(a,i1,a)')
     &          '(i4,2f8.2,3f8.4,f9.3,f10.3,',
     &                          NDIFO,'f8.2,f10.3)'
              ENDIF !ndifo
            ELSE !NTRC.ne.0
              IF     (NDIFO.EQ.0) THEN
                WRITE(CFORMAT,'(a,i2,a)')
     &          '(i4,2f8.2,3f8.4,f9.3,f10.3,',NTRC,'f8.4,f10.3)'
              ELSE
                WRITE(CFORMAT,'(a,i1,a,i2,a)')
     &          '(i4,2f8.2,3f8.4,f9.3,f10.3,',
     &                          NDIFO,'f8.2,',NTRC,'f8.4,f10.3)'
              ENDIF !ndifo
            ENDIF !ntrc
C
            DIF(1:3) = 0.0
            DO K= 1,KDM
*             write(IU,*) 'PROF = ',PROF
              UKIJ = MAX( -9.99, MIN( 9.99, UK(K,I) ) )
              VKIJ = MAX( -9.99, MIN( 9.99, VK(K,I) ) )
              IF     (NDIF.EQ.3) THEN
                DIF(1) = MIN(9999.99, PROF(NLAY+NTRC+1,K,I)*1.E4)
                DIF(2) = MIN(9999.99, PROF(NLAY+NTRC+2,K,I)*1.E4)
                DIF(3) = MIN(9999.99, PROF(NLAY+NTRC+3,K,I)*1.E4)
              ENDIF
              IF     (SIGVER.EQ.0) THEN
                WRITE(IU,CFORMAT)
     &            K,
     &            UKIJ*100.0,                          !cm/s
     &            VKIJ*100.0,                          !cm/s
     &            PROF(NLAY-2,K,I),                    !degC
     &            PROF(NLAY-1,K,I),                    !psu
     &            PROF(NLAY,  K,I)  +THBASE,           !SigmaT
     &            (P(K+1,I)-P(K,I))    /ONEM,          !m
     &            (P(K+1,I)+P(K,I))*0.5/ONEM,          !m
     &            (DIF(KTR),KTR=1,NDIFO),              !cm^2/s
     &            (PROF(NLAY+KTR,K,I),KTR=1,NTRC),     !0-99
     &             PROF(NLAY+NDIF+NTRC1,K,I)           !thknss
              ELSE
                WRITE(IU,CFORMAT)
     &            K,
     &              UKIJ*100.0,                          !cm/s
     &            VKIJ*100.0,                          !cm/s
     &            PROF(NLAY-1,K,I),                    !degC
     &            PROF(NLAY,  K,I),                    !psu
     &                     RK(K,I),                    !SigmaT
     &            PROF(NLAY-2,K,I)/ONEM,               !m
     &            (P(K+1,I)+P(K,I))*0.5/ONEM,          !m
     &            (DIF(KTR),KTR=1,NDIFO),              !cm^2/s
     &            (PROF(NLAY+KTR,K,I),KTR=1,NTRC),     !0-99
     &             PROF(NLAY+NDIF+NTRC1,K,I)           !thknss
              ENDIF
            ENDDO !k
          ENDIF !ntrc1
          IF     (IU.NE.6) THEN
            CLOSE(IU)
          ENDIF
      ENDDO !I
      CALL EXIT(0)
      END
