      PROGRAM HYCOM_BINNING_NC
      IMPLICIT NONE
C
C  hycom_binning_nc - Usage:  hycom_binning_nc file.a cell [lonmin lonmax latmin latmax] [grid.a] [mask.a]
C                             hycom_binning_nc SAMPLE cell [lonmin lonmax latmin latmax] [grid.a]
C
C                 creates a netDCF file containing the array of means of
C                 lon/lat bins from all the (1:idm,1:jdm) arrays in file.a
C                 if file.a ls "SAMPLE" the input is a sequence of plain
C                 text lon,lat,value triplets in a separate file per field
C                 with the filenames in environment variable CDF_SAMPLE
C
C                 The input array sequence is assumed to contain
C                 CDF_NLAY layers from each time, or if CDF_NLAY is
C                 0 (the default) the arrays represent the same
C                 field over time.
C
C                 use hycom_binning for plain text output (single array).
C
C                 NetCDF environment variables:
C                    CDF_FILE   ndcf filename
C                    CDF_TITLE  title
C                    CDF_INST   institution
C                    CDF_CNAME  field name (short, no spaces)
C                    CDF_PNAME  plot  name (longer,   spaces)
C                    CDF_SNAME  field standard name
C                    CDF_UNAME  field units
C                    CDF_NLAY   number of layers (default 0)
C                    CDF_DAY1   first day
C                    CDF_DAYI   increment between fields in days
C                    CDF_YEAR   type of year (0=360,1=366,2=366Jan1,3=actual)
C                    CDF_ZFILE  filename, containing list of z-levels
C                    CDF_DFILE  filename, containing list of days
C                    CDF_SAMPLE filename, containing list of sample filenames
C
C                 longitude bands are centered on 0E and latitude bands
C                 are centered on the equator, both are nominally cell
C                 degrees wide, e.g. (N-0.5)*cell to (N+0.5)*cell.
C                 However, the cell sample width is (approximately) 
C                 scaled by 1/cos(lat) in longitude, i.e. the cell
C                 width in meters is aproximately constant.
C
C                 lonmin lonmax latmin latmax define the extent of
C                 the sampled subregion, default the entire region.
C
C                 grid.a is a hycom grid file, default regional.grid.a.
C                 Note that the corresponding grid.b must also exist.
C
C                 mask.a contains an optional mask array, that turns
C                 off sampling anywhere it is 2.0**100 (data void).
C                 Note that grid.a must be specified if mask.a is used.
C
C                 idm,jdm are taken from grid.a and the input array
C                 is assumed to be on the p-grid (i.e. plat from grid.a
C                 represents the field's latitude).
C
C  file.a is assumed to contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  this version for "serial" Unix systems.
C  Based on hycom_binning.F.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  August 2004.
C
      REAL*4, ALLOCATABLE :: A(:,:),AMSK(:,:),
     +                       PLAT(:,:),PLON(:,:),PSC2(:,:)
      REAL*4              :: PAD(4096)
C
      INTEGER       IOS
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      INTEGER       IDM,JDM,NPAD
      REAL*4        CELL,LONMIN,LONMAX,LATMIN,LATMAX
      CHARACTER*6   CVARIN
      CHARACTER*240 CFILE,CFILEG,CFILEM,CFILEB
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.2) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*)    CELL
        CFILEG = 'regional.grid.a'
        CFILEM = 'NONE'
        LATMAX = 999.0  ! indicates entire region
      ELSEIF (NARG.EQ.3) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*)    CELL
        CALL GETARG(3,CFILEG)
        CFILEM = 'NONE'
        LATMAX = 999.0  ! indicates entire region
      ELSEIF (NARG.EQ.4) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*)    CELL
        CALL GETARG(3,CFILEG)
        CALL GETARG(4,CFILEM)
        LATMAX = 999.0  ! indicates entire region
      ELSEIF (NARG.EQ.6) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*)    CELL
        CALL GETARG(3,CARG)
        READ(CARG,*)    LONMIN
        CALL GETARG(4,CARG)
        READ(CARG,*)    LONMAX
        CALL GETARG(5,CARG)
        READ(CARG,*)    LATMIN
        CALL GETARG(6,CARG)
        READ(CARG,*)    LATMAX
        CFILEG = 'regional.grid.a'
        CFILEM = 'NONE'
      ELSEIF (NARG.EQ.7) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*)    CELL
        CALL GETARG(3,CARG)
        READ(CARG,*)    LONMIN
        CALL GETARG(4,CARG)
        READ(CARG,*)    LONMAX
        CALL GETARG(5,CARG)
        READ(CARG,*)    LATMIN
        CALL GETARG(6,CARG)
        READ(CARG,*)    LATMAX
        CALL GETARG(7,CFILEG)
        CFILEM = 'NONE'
      ELSEIF (NARG.EQ.8) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*)    CELL
        CALL GETARG(3,CARG)
        READ(CARG,*)    LONMIN
        CALL GETARG(4,CARG)
        READ(CARG,*)    LONMAX
        CALL GETARG(5,CARG)
        READ(CARG,*)    LATMIN
        CALL GETARG(6,CARG)
        READ(CARG,*)    LATMAX
        CALL GETARG(7,CFILEG)
        CALL GETARG(8,CFILEM)
      ELSE
        WRITE(6,*) 
     +   'Usage:  hycom_binning_nc file.a cell ' //
     +   '[lonmin lonmax latmin latmax] [grid.a] [mask.a]'
        CALL EXIT(1)
      ENDIF
C
C     GET IDM,JDM FROM grid.b.
C
      CFILEB = CFILEG(1:LEN_TRIM(CFILEG)-1) // 'b'
C
      OPEN(UNIT=11,FILE=CFILEB,FORM='FORMATTED',
     &     STATUS='OLD',ACTION='READ')
C
      READ( 11,*) IDM,CVARIN
      IF (CVARIN.NE.'idm   ') THEN
        WRITE(6,*) 'hycom_binning_nc: bad header file ',
     &             TRIM(CFILEB)
        CALL EXIT(2)
      ENDIF
      READ( 11,*) JDM,CVARIN
      IF (CVARIN.NE.'jdm   ') THEN
        WRITE(6,*) 'hycom_binning_nc: bad header file ',
     &             TRIM(CFILEB)
        CALL EXIT(2)
      ENDIF
C
      CLOSE(UNIT=11)
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( A(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_zonal: could not allocate ',
     +             IDM*JDM,' words for A'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( AMSK(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_zonal: could not allocate ',
     +             IDM*JDM,' words for AMSK'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( PLAT(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_zonal: could not allocate ',
     +             IDM*JDM,' words for PLAT'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( PLON(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_zonal: could not allocate ',
     +             IDM*JDM,' words for PLON'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( PSC2(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_zonal: could not allocate ',
     +             IDM*JDM,' words for PSC2'
        CALL EXIT(2)
      ENDIF
C
      CALL BINNING(A,AMSK,PLAT,PLON,PSC2,IDM,JDM,PAD,NPAD,
     +             CELL,LONMIN,LONMAX,LATMIN,LATMAX,
     +             CFILE,CFILEG,CFILEM)
      CALL EXIT(0)
 5000 FORMAT(I4)
      END
      SUBROUTINE BINNING(A,AMSK,PLAT,PLON,PSC2,IDM,JDM, PAD,NPAD,
     +                   CELL,LONMIN,LONMAX,LATMIN,LATMAX,
     +                   CFILE,CFILEG,CFILEM)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE,CFILEG,CFILEM
      INTEGER      IDM,JDM,NPAD
      REAL*4       A(IDM,JDM),AMSK(IDM,JDM),
     +             PLAT(IDM,JDM),PLON(IDM,JDM),PSC2(IDM,JDM),PAD(NPAD),
     +             CELL,LONMIN,LONMAX,LATMIN,LATMAX
C
C     MOST OF WORK IS DONE HERE.
C
      REAL*8,  ALLOCATABLE :: SUMA(:,:),SUMB(:,:)
      REAL,    ALLOCATABLE :: AMEAN(:,:),ALON(:),ALAT(:)
      INTEGER, ALLOCATABLE :: IW(:)
C
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      CHARACTER*240 CNAMES,CFILES
      CHARACTER*18  CASN
      LOGICAL       LSUBREG,LGLOBAL
      INTEGER       I,IC,IF,IL,IZ,IP,IQ,J,JC,JF,JL,JZ,IOS,K,NRECL
      REAL*8        SUMAT,SUMBT
      REAL*4        QA,QLAT,QLON,XLON,DEG2RAD,MINCOS
#ifdef CRAY
      INTEGER*8     IU8,IOS8
#endif
C
C     INPUT ARRAY A.
C
      INQUIRE( IOLENGTH=NRECL) A,PAD
#ifdef CRAY
#ifdef t3e
      IF     (MOD(NRECL,4096).EQ.0) THEN
        WRITE(CASN,8000) NRECL/4096
 8000   FORMAT('-F cachea:',I4.4,':1:0')
        IU8 = 11
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          write(6,*) 'Error: can''t asnunit 11'
          write(6,*) 'ios  = ',ios8
          write(6,*) 'casn = ',casn
          CALL EXIT(5)
        ENDIF
      ENDIF
#else
      CALL ASNUNIT(11,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 11'
        write(6,*) 'ios = ',ios
        CALL EXIT(5)
      ENDIF
#endif
#endif
C
C     INPUT MASK ARRAY.
C
      IF     (CFILEM.EQ.'NONE') THEN
        AMSK(:,:) = 0.0
      ELSE
        OPEN(UNIT=11, FILE=CFILEM, FORM='UNFORMATTED', STATUS='OLD',
     +           ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
        IF     (IOS.NE.0) THEN
          write(6,*) 'Error: can''t open ',TRIM(CFILEM)
          write(6,*) 'ios   = ',ios
          write(6,*) 'nrecl = ',nrecl
          CALL EXIT(3)
        ENDIF
C
        READ(11,REC=1,IOSTAT=IOS) AMSK
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(AMSK,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'can''t read ',TRIM(CFILEM)
          CALL EXIT(4)
        ENDIF
      ENDIF
C
      CLOSE(UNIT=11)
C
C     INPUT GRID ARRAYS.
C
      OPEN(UNIT=11, FILE=CFILEG, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEG)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
      READ(11,REC=10,IOSTAT=IOS) PLAT  ! pscx
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(PLAT,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read ',TRIM(CFILEG)
        CALL EXIT(4)
      ENDIF
C
      READ(11,REC=11,IOSTAT=IOS) PSC2  ! pscy
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(PSC2,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read ',TRIM(CFILEG)
        CALL EXIT(4)
      ENDIF
C
      PSC2(:,:) = PSC2(:,:)*PLAT(:,:)  ! psc2=pscx*pscy
C
      READ(11,REC= 1,IOSTAT=IOS) PLON  ! plon
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(PLON,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read ',TRIM(CFILEG)
        CALL EXIT(4)
      ENDIF
C
      READ(11,REC= 2,IOSTAT=IOS) PLAT  ! plat
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(PLAT,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read ',TRIM(CFILEG)
        CALL EXIT(4)
      ENDIF
C
      CLOSE(UNIT=11)
C
C     SUBREGION?
C
      LSUBREG = LATMAX.LT.91.0  ! sub-region specified
      LGLOBAL = .FALSE.         ! over-ride below if necessary
      IF     (.NOT.LSUBREG) THEN
C
C      ASSUME THAT LONGITUDES ARE EITHER GLOBAL OR COMPACT
C
        IF     (PLON(IDM,1)-PLON(1,1).GT.350.0) THEN  !global
          LONMIN  = -0.5*CELL
          LONMAX  = (NINT(360.0/CELL)-1.5)*CELL 
          LGLOBAL = .TRUE.
        ELSE
          LONMIN = MINVAL(PLON(:,:))
          LONMAX = MAXVAL(PLON(:,:))
        ENDIF
      ELSE
C
C       MOVE LONMAX WITHIN +360 OF LONMIN
C
        IF     (LONMAX.LT.LONMIN) THEN
          DO
            LONMAX = LONMAX + 360.0
            IF     (LONMAX.GT.LONMIN) THEN
              EXIT
            ENDIF
          ENDDO
        ELSEIF (LONMAX.GE.LONMIN+360.0) THEN
          DO
            LONMAX = LONMAX - 360.0
            IF     (LONMAX.LT.LONMIN+360.0) THEN
              EXIT
            ENDIF
          ENDDO
        ENDIF
      ENDIF
C
C     MOVE PLON WITHIN +360 OF LONMIN, AND
C     CONVERT LONGITUDE TO CELL SPACE.
C
      DO J= 1,JDM
        DO I = 1,IDM
          XLON = PLON(I,J)
          IF     (XLON.LT.LONMIN) THEN
            DO
              XLON = XLON + 360.0
              IF     (XLON.GT.LONMIN) THEN
                EXIT
              ENDIF
            ENDDO
          ELSEIF (XLON.GE.LONMIN+360.0) THEN
            DO
              XLON = XLON - 360.0
              IF     (XLON.LT.LONMIN+360.0) THEN
                EXIT
              ENDIF
            ENDDO
          ENDIF
          PLON(I,J) = NINT(XLON/CELL)
        ENDDO !i
      ENDDO !j
      IF = MINVAL(PLON(:,:))
      IL = MAXVAL(PLON(:,:))
C
      IF     (LSUBREG) THEN
        IF = MAX( IF, NINT(LONMIN/CELL) )
        IL = MIN( IL, NINT(LONMAX/CELL) )
      ENDIF
C
C     CONVERT LATITUDE TO CELL SPACE.
C
      PLAT(:,:) = NINT(PLAT(:,:)/CELL)
      JF = MINVAL(PLAT(:,:))
      JL = MAXVAL(PLAT(:,:))
C
      IF     (LSUBREG) THEN
        JF = MAX( JF, NINT(LATMIN/CELL) )
        JL = MIN( JL, NINT(LATMAX/CELL) )
      ENDIF
C
      IC = (IF+IL)/2
      JC = (JF+JL)/2
C
C     ALLOCATE SUM ARRAYS IN CELL SPACE.
C
      ALLOCATE(  SUMA(IF:IL,JF:JL),
     +           SUMB(IF:IL,JF:JL) )
      ALLOCATE( AMEAN(IF:IL,JF:JL) )
      ALLOCATE( ALON(IF:IL) )
      ALLOCATE( ALAT(JF:JL) )
      ALLOCATE(   IW(JF:JL) )
C
C     PRINT A HEADER.
C
      IF     (.NOT. LSUBREG) THEN ! full region
        WRITE(6,6000)
     +    TRIM(CFILE),
     +    CELL,
     +    TRIM(CFILEG),
     +    TRIM(CFILEM)
      ELSE
        WRITE(6,6100)
     +    TRIM(CFILE),
     +    CELL,
     +    LONMIN,LONMAX,LATMIN,LATMAX,
     +    TRIM(CFILEG),
     +    TRIM(CFILEM)
      ENDIF
C
C     INPUT FILE
C
      IF     (CFILE.NE."SAMPLE") THEN
        OPEN(UNIT=11, FILE=CFILE, FORM='UNFORMATTED', STATUS='OLD',
     +           ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
        IF     (IOS.NE.0) THEN
          write(6,*) 'Error: can''t open ',TRIM(CFILE)
          write(6,*) 'ios   = ',ios
          write(6,*) 'nrecl = ',nrecl
          CALL EXIT(3)
        ENDIF
      ELSE
        CNAMES = ' '
        CALL GETENV('CDF_SAMPLE',CNAMES)
        IF     (CNAMES.EQ.' ') THEN
          write(6,*) 'Error: CDF_SAMPLE is undefined'
          CALL EXIT(3)
        ENDIF
        OPEN(UNIT=12, FILE=CNAMES, FORM='FORMATTED', STATUS='OLD',
     +           IOSTAT=IOS)
        IF     (IOS.NE.0) THEN
          write(6,*) 'Error: can''t open ',TRIM(CNAMES)
          CALL EXIT(3)
        ENDIF
      ENDIF
C
C     LONGITUDE CELL WIDTH (1+2*IW)
C
      DEG2RAD = 4.D0*ATAN(1.D0)/180.D0  !PI/180
      MINCOS  = 1.0/REAL(NINT(360.0/CELL)-1)
      DO JZ= JF,JL
        IW(JZ) = INT( 0.5/MAX(  MINCOS,
     &                         ABS(COS(JZ*CELL*DEG2RAD)) ) )
      ENDDO !jz
C
      DO IZ= IF,IL
        ALON(IZ) = IZ*CELL
      ENDDO
      DO JZ= JF,JL
        ALAT(JZ) = JZ*CELL
      ENDDO
      write(6,*) 
      write(6,*) 'progrm - alon = ',alon
      write(6,*) 'progrm - alat = ',alat
      write(6,*) 
C
C     ALL RECORDS IN THE FILE
C
      DO K= 1,HUGE(K)/2  !safe limit
        IF     (CFILE.NE."SAMPLE") THEN
          READ(11,REC=K,IOSTAT=IOS) A
#ifdef ENDIAN_IO
          CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
          IF     (K.EQ.1 .AND. IOS.NE.0) THEN
            WRITE(6,*) 'can''t read ',TRIM(CFILE)
            CALL EXIT(4)
          ELSEIF (IOS.NE.0) THEN
            CLOSE(UNIT=11)
            EXIT
          ENDIF
C
C         DO ALL BINS.
C
          SUMAT   = 0.0
          SUMBT   = 0.0
          SUMA(:,:) = 0.0
          SUMB(:,:) = 0.0
          DO J= 1,JDM
            DO I= 1,IDM
              IF     (   A(I,J).NE.SPVAL .AND.
     +                AMSK(I,J).NE.SPVAL      ) THEN
                IZ = PLON(I,J)
                JZ = PLAT(I,J)
                IF     (IZ.LT.IF .OR. IZ.GT.IL .OR.
     +                  JZ.LT.JF .OR. JZ.GT.JL     ) THEN
                  CYCLE !i
                ENDIF
C
                SUMAT = SUMAT + A(I,J)*PSC2(I,J)
                SUMBT = SUMBT +        PSC2(I,J)
C
C               ALLOW FOR CELL WIDTH
C
                DO IP= IZ-IW(JZ),IZ+IW(JZ)
                  IF     (LGLOBAL) THEN
                    IQ = MOD(IP,IL+1)
                  ELSEIF (IP.LT.IF .OR. IP.GT.IL) THEN
                    CYCLE !ip
                  ELSE
                    IQ = IP
                  ENDIF
                  SUMA(IQ,JZ) = SUMA(IQ,JZ) + A(I,J)*PSC2(I,J)
                  SUMB(IQ,JZ) = SUMB(IQ,JZ) +        PSC2(I,J)
*
*                 IF     (IQ.EQ.IC .AND. JZ.EQ.JC) THEN
*                   WRITE(6,'(a,2i5,1pg16.8)') 'debug: ',I,J,A(I,J)
*                 ENDIF
                ENDDO !ip
              ENDIF
            ENDDO !jz
          ENDDO !iz
        ELSE !SAMPLE
          READ(12,'(a)',IOSTAT=IOS) CFILES
          IF     (K.EQ.1 .AND. IOS.NE.0) THEN
            WRITE(6,*) 'can''t read ',TRIM(CNAMES)
            CALL EXIT(4)
          ELSEIF (IOS.NE.0) THEN
            CLOSE(UNIT=12)
            EXIT
          ENDIF
          OPEN(UNIT=11, FILE=CFILES, FORM='FORMATTED', STATUS='OLD',
     +             IOSTAT=IOS)
          IF     (IOS.NE.0) THEN
            write(6,*) 'Error: can''t open ',TRIM(CFILES)
            CALL EXIT(3)
          ENDIF
C
C         DO ALL BINS.
C
          SUMAT     = 0.0
          SUMBT     = 0.0
          SUMA(:,:) = 0.0
          SUMB(:,:) = 0.0
          DO I= 1,HUGE(K)/2  !safe limit
            READ(11,*,IOSTAT=IOS) QLON,QLAT,QA
            IF     (IOS.NE.0) THEN
              CLOSE(UNIT=11)
              EXIT
            ENDIF
            IF     (QLON.LT.LONMIN) THEN
              DO
                QLON = QLON + 360.0
                IF     (QLON.GT.LONMIN) THEN
                  EXIT
                ENDIF
              ENDDO
            ELSEIF (QLON.GE.LONMIN+360.0) THEN
              DO
                QLON = QLON - 360.0
                IF     (QLON.LT.LONMIN+360.0) THEN
                  EXIT
                ENDIF
              ENDDO
            ENDIF
            IZ = NINT(QLON/CELL)
            JZ = NINT(QLAT/CELL)
            IF     (IZ.LT.IF .OR. IZ.GT.IL .OR.
     +              JZ.LT.JF .OR. JZ.GT.JL     ) THEN
              CYCLE !i
            ENDIF
C
            SUMAT = SUMAT + QA
            SUMBT = SUMBT + 1.0
C
C           ALLOW FOR CELL WIDTH
C
            DO IP= IZ-IW(JZ),IZ+IW(JZ)
              IF     (LGLOBAL) THEN
                IQ = MOD(IP,IL+1)
              ELSEIF (IP.LT.IF .OR. IP.GT.IL) THEN
                CYCLE !ip
              ELSE
                IQ = IP
              ENDIF
              SUMA(IQ,JZ) = SUMA(IQ,JZ) + QA
              SUMB(IQ,JZ) = SUMB(IQ,JZ) + 1.0
*
*             IF     (IQ.EQ.IC .AND. JZ.EQ.JC) THEN
*               WRITE(6,'(a,i6,1pg16.8)') 'debug: ',I,QA
*             ENDIF
            ENDDO !ip
          ENDDO !iz
        ENDIF !.not.SAMPLE;else
C
        DO JZ= JF,JL
          DO IZ= IF,IL
            IF     (SUMB(IZ,JZ).NE.0.0) THEN
              AMEAN(IZ,JZ) = SUMA(IZ,JZ)/SUMB(IZ,JZ)
            ELSE
              AMEAN(IZ,JZ) = 2.0**100
            ENDIF
          ENDDO !jz
        ENDDO !iz
c
        CALL HOROUT(AMEAN,ALON,ALAT,IL-IF+1,JL-JF+1,CFILE)
C
        SUMAT = SUMAT/SUMBT
        WRITE(6,'(a/a,1pg16.8)')
     +     '#    LON    LAT        MEAN',
     +     '#    ALL    ALL',SUMAT
      ENDDO !k
      RETURN
 6000 FORMAT(
     +   '# hycom_binning_nc file.a cell grid.a mask.a' /
     +   '# hycom_binning_nc ',A       /
     +   '#               ',F7.2 /
     +   '#               ',A       /
     +   '#               ',A       /
     +   '#    LON    LAT        MEAN')
 6100 FORMAT(
     +   '# hycom_binning_nc file.a cell ',
     +   'lonmin lonmax latmin latmax grid.a mask.a' /
     +   '# hycom_binning_nc ',A       /
     +   '#               ',F7.2 /
     +   '#               ',2F9.2,2F8.2 /
     +   '#               ',A       /
     +   '#               ',A       /
     +   '#    LON    LAT        MEAN')
      END

      subroutine horout(array,plon,plat,ii,jj,cfile)
      use netcdf   ! NetCDF fortran 90 interface
      implicit none
c
      integer       ii,jj
      real          array(ii,jj),plon(ii),plat(jj)
      character*(*) cfile
c
c     NetCDF environment variables:
c       CDF_FILE  ndcf filename
c       CDF_TITLE title
c       CDF_INST  institution
c       CDF_CNAME field name (short, no spaces)
c       CDF_PNAME plot  name (longer,   spaces)
c       CDF_SNAME field standard name
c       CDF_UNAME field units
c       CDF_NLAY  number of layers (default 0)
c       CDF_DAY1  first day
c       CDF_DAYI  increment between fields in days
c       CDF_YEAR  type of year
c                 (HYCOM yrflag: 0=360,1=366,2=366Jan1,3=actual,4=365Jan1)
c       CDF_ZFILE filename, containing list of z-levels
c       CDF_DFILE filename, containing list of days
c
c     This routine needs version 3.5 of the NetCDF library, from: 
c     http://www.unidata.ucar.edu/packages/netcdf/
c
      integer          :: ncfileID, status, varID
      integer          :: pLatDimID,pLonDimID,LayDimID,
     &                    pLatVarID,pLonVarID,LayVarID
      integer          :: MTDimID,MTVarID,datVarID
      character*240    :: ncfile,ncenv
      character*240    :: name,namec,names,units
c
      logical          :: lopen,lexist
      integer          :: i,j,k,kk,l,iyear,month,iday,ihour,
     &                               iyrms,monms,idms,ihrms,yrflag
      real             :: hmin,hmax,hrange(2)
      real             :: zz(999)
      double precision :: time,year,wday1,wdayi,date
c
      integer          :: kk_rec     = 0  ! saved below
      integer          :: mt_rec     = 0  ! saved below
      real,    parameter :: fill_value = 2.0**100
c
      save  !save all local variables
c
      if     (mt_rec.eq.0) then
c
c       initialization.
c
        write(6,*) 
        write(6,*) 'horout - plon = ',plon
        write(6,*) 'horout - plat = ',plat
        write(6,*) 
c
        kk_rec = 1
        mt_rec = 1
c
        namec = ' '
        call getenv('CDF_CNAME',namec)
        if     (namec.eq.' ') then
          namec = "fld"
        endif
        name  = ' '
        call getenv('CDF_PNAME',name)
        if     (name .eq.' ') then
          name  = "field"
        endif
        names = ' '
        call getenv('CDF_SNAME',names)
        if     (names.eq.' ') then
          names = ' '
        endif
        units = ' '
        call getenv('CDF_UNAME',units)
        if     (units.eq.' ') then
          units = ' '
        endif
        write(6,'(2a)') 'horout - namec=',trim(namec)
        write(6,'(2a)') 'horout - name =',trim( name)
        write(6,'(2a)') 'horout - names=',trim(names)
        write(6,'(2a)') 'horout - units=',trim(units)
c
        ncenv = ' '
        call getenv('CDF_NLAY',ncenv)
        if     (ncenv.eq.' ') then
          kk   = 0
        else
          read(ncenv,*) kk
        endif
        write(6,'(a,i4)') 'horout - kk   =',kk
c
        if     (kk.gt.0) then
          ncenv = ' '
          call getenv('CDF_ZFILE',ncenv)
          if     (ncenv.eq.' ') then
            zz(1) = -1.0
          else
            open(unit=21,file=ncenv,form='formatted', status='old')
            do k= 1,kk
              read(21,*) zz(k)
              write(6,'(a,f12.3)') 'horout - zz   =',zz(k)
            enddo
            close(21)
          endif
        endif
c
        ncenv = ' '
        call getenv('CDF_YEAR',ncenv)
        if     (ncenv.eq.' ') then
          yrflag = -1.0
        else
          read(ncenv,*) yrflag
        endif
c
        ncenv = ' '
        call getenv('CDF_DFILE',ncenv)
        if     (ncenv.eq.' ') then
          ncenv = ' '
          call getenv('CDF_DAY1',ncenv)
          if     (ncenv.eq.' ') then
            wday1 = 1.0
          else
            read(ncenv,*) wday1
          endif
          ncenv = ' '
          call getenv('CDF_DAYI',ncenv)
          if     (ncenv.eq.' ') then
            wdayi = 1.0
          else
            read(ncenv,*) wdayi
          endif
        else
          open(unit=22,file=ncenv,form='formatted', status='old')
          read(22,*) wday1
          wdayi = 0.0
        endif
        write(6,'(a,f12.3)') 'horout - wday1=',wday1
        write(6,'(a,f12.3)') 'horout - wdayi=',wday1
        write(6,'(a,i9)')    'horout - yrflg=',yrflag
c
        write( 6,'(/2a/)') 'horout - NetCDF I/O (lat/lon axes)'
c
c       NetCDF I/O
c
        time = wday1
        if     (yrflag.ge.0) then
          call fordate(time,yrflag, iyear,month,iday,ihour)
          date = (iday + 100 * month + 10000 * iyear) + ihour/24.d0
        else
          date = mt_rec
        endif
c
        write(6,6300) kk_rec,mt_rec,time,date
c
        ncfile = ' '
        call getenv('CDF_FILE',ncfile)
        if     (ncfile.eq.' ') then
          write( 6,'(/a/)')  'error in horout - CDF_FILE not defined'
          stop
        endif
c
        call ncrange(array,ii,jj,1, fill_value, hmin,hmax)
c
        inquire(file= ncfile, exist=lexist)
        if (lexist) then
          write( 6,'(/2a/a/)') 'error in horout - ',
     &                        'CDF_FILE is an existing file',
     &                        trim(ncfile)
          stop
        else
c
c         create a new NetCDF and write data to it
c         netcdf-4 classic model, netcdf version 4.3 and later
c
          call nchek('nf90_create',
     &                nf90_create(trim(ncfile),
     &                            or(nf90_clobber,
     &                               or(nf90_hdf5,
     &                                  nf90_classic_model)),
     &                            ncfileID))
          ! define the dimensions
          call nchek("nf90_def_dim-MT",
     &                nf90_def_dim(ncfileID,
     &                             "MT", nf90_unlimited,MTDimID))
          if     (kk.ne.0) then
            call nchek("nf90_def_dim-Layer",
     &                  nf90_def_dim(ncfileID,"Layer",kk,LayDimID))
          endif
            call nchek("nf90_def_dim-Latitude",
     &                  nf90_def_dim(ncfileID,
     &                               "Latitude",  jj,pLatDimID))
            call nchek("nf90_def_dim-Longitude",
     &                  nf90_def_dim(ncfileID,
     &                               "Longitude", ii,pLonDimID))
          ! create the global attributes
          call nchek("nf90_put_att-Conventions",
     &                nf90_put_att(ncfileID,nf90_global,
     &                             "Conventions",
     &                             "CF-1.0"))
            ncenv = ' '
            call getenv('CDF_TITLE',ncenv)
            if     (ncenv.eq.' ') then
              ncenv = "HYCOM"
            endif
            call nchek("nf90_put_att-title",
     &                  nf90_put_att(ncfileID,nf90_global,
     &                               "title",
     &                               trim(ncenv)))
            ncenv = ' '
            call getenv('CDF_INST',ncenv)
            if     (ncenv.ne.' ') then
              call nchek("nf90_put_att-institution",
     &                    nf90_put_att(ncfileID,nf90_global,
     &                                 "institution",
     &                                 trim(ncenv)))
            endif
            call nchek("nf90_put_att-source",
     &                  nf90_put_att(ncfileID,nf90_global,
     &                               "source",
     &                               trim(cfile)))
            call nchek("nf90_put_att-history",
     &                  nf90_put_att(ncfileID,nf90_global,
     &                               "history",
     &                               "hycom_binning_nc"))
          ! create the variables and attributes
            call nchek("nf90_def_var-MT",
     &                  nf90_def_var(ncfileID,"MT",  nf90_double,
     &                               MTDimID,MTVarID))
            if     (yrflag.eq.0) then
              call nchek("nf90_put_att-long_name",
     &                    nf90_put_att(ncfileID,MTVarID,
     &                                 "long_name",
     &                                 "model time"))
              call nchek("nf90_put_att-units",
     &                    nf90_put_att(ncfileID,MTVarID,
     &                                 "units",
     &                            "days since 0001-01-16 00:00:00"))
              call nchek("nf90_put_att-calendar",
     &                    nf90_put_att(ncfileID,MTVarID,
     &                                 "calendar",
     &                                 "360_day"))
            elseif (yrflag.eq.1) then
              call nchek("nf90_put_att-long_name",
     &                    nf90_put_att(ncfileID,MTVarID,
     &                                 "long_name",
     &                                 "model time"))
              call nchek("nf90_put_att-units",
     &                    nf90_put_att(ncfileID,MTVarID,
     &                                 "units",
     &                            "days since 0001-01-16 00:00:00"))
              call nchek("nf90_put_att-calendar",
     &                    nf90_put_att(ncfileID,MTVarID,
     &                                 "calendar",
     &                                 "366_day"))
            elseif (yrflag.eq.2) then
              call nchek("nf90_put_att-long_name",
     &                    nf90_put_att(ncfileID,MTVarID,
     &                                 "long_name",
     &                                 "model time"))
              call nchek("nf90_put_att-units",
     &                    nf90_put_att(ncfileID,MTVarID,
     &                                 "units",
     &                            "days since 0001-01-01 00:00:00"))
              call nchek("nf90_put_att-calendar",
     &                    nf90_put_att(ncfileID,MTVarID,
     &                                 "calendar",
     &                                 "366_day"))
            elseif (yrflag.eq.4) then
              call nchek("nf90_put_att-long_name",
     &                    nf90_put_att(ncfileID,MTVarID,
     &                                 "long_name",
     &                                 "model time"))
              call nchek("nf90_put_att-units",
     &                    nf90_put_att(ncfileID,MTVarID,
     &                                 "units",
     &                            "days since 0001-01-01 00:00:00"))
              call nchek("nf90_put_att-calendar",
     &                    nf90_put_att(ncfileID,MTVarID,
     &                                 "calendar",
     &                                 "365_day"))
            elseif (yrflag.eq.3) then
              call nchek("nf90_put_att-long_name",
     &                    nf90_put_att(ncfileID,MTVarID,
     &                                 "long_name",
     &                                 "time"))
              call nchek("nf90_put_att-units",
     &                    nf90_put_att(ncfileID,MTVarID,
     &                                 "units",
     &                            "days since 1900-12-31 00:00:00"))
              call nchek("nf90_put_att-calendar",
     &                    nf90_put_att(ncfileID,MTVarID,
     &                                 "calendar",
     &                                 "standard"))
            else
              call nchek("nf90_put_att-long_name",
     &                    nf90_put_att(ncfileID,MTVarID,
     &                                 "long_name",
     &                                 "record number"))
            endif !yrflag
            call nchek("nf90_put_att-axis",
     &                  nf90_put_att(ncfileID,MTVarID,
     &                               "axis","T"))
            if     (kk.gt.0 .and. zz(1).eq.-1.0) then
              call nchek("nf90_def_var-Layer",
     &                    nf90_def_var(ncfileID,"Layer", nf90_int,
     &                                 LayDimID,LayVarID))
              call nchek("nf90_put_att-units",
     &                    nf90_put_att(ncfileID,LayVarID,
     &                                 "units","layer"))
              call nchek("nf90_put_att-positive",
     &                    nf90_put_att(ncfileID,LayVarID,
     &                                 "positive","down"))
              call nchek("nf90_put_att-axis",
     &                    nf90_put_att(ncfileID,LayVarID,
     &                                 "axis","Z"))
            elseif (kk.gt.0) then !z-depths
              call nchek("nf90_def_var-Depth",
     &                    nf90_def_var(ncfileID,"Depth", nf90_float,
     &                                 LayDimID,LayVarID))
              call nchek("nf90_put_att-units",
     &                    nf90_put_att(ncfileID,LayVarID,
     &                                 "units","m"))
              call nchek("nf90_put_att-positive",
     &                    nf90_put_att(ncfileID,LayVarID,
     &                                 "positive","down"))
              call nchek("nf90_put_att-axis",
     &                    nf90_put_att(ncfileID,LayVarID,
     &                                 "axis","Z"))
            endif !kk>0
            if     (yrflag.ge.0) then
              call nchek("nf90_def_var-Date",
     &                    nf90_def_var(ncfileID,"Date", nf90_double,
     &                                 MTDimID,datVarID))
              call nchek("nf90_put_att-long_name",
     &                    nf90_put_att(ncfileID,datVarID,
     &                                 "long_name",
     &                                 "date"))
              call nchek("nf90_put_att-units",
     &                    nf90_put_att(ncfileID,datVarID,
     &                                 "units",
     &                                 "day as %Y%m%d.%f"))
              call nchek("nf90_put_att-C_format",
     &                    nf90_put_att(ncfileID,datVarID,
     &                                 "C_format",
     &                                 "%13.4f"))
              call nchek("nf90_put_att-FORTRAN_format",
     &                    nf90_put_att(ncfileID,datVarID,
     &                                 "FORTRAN_format",
     &                                 "(f13.4)"))
            endif !yrflag>0
              call nchek("nf90_def_var-Latitude",
     &                    nf90_def_var(ncfileID,"Latitude",  nf90_float,
     &                                 pLatDimID,pLatVarID))
            call nchek("nf90_put_att-standard_name",
     &                  nf90_put_att(ncfileID,pLatVarID,
     &                               "standard_name","latitude"))
            call nchek("nf90_put_att-units",
     &                  nf90_put_att(ncfileID,pLatVarID,
     &                               "units","degrees_north"))
            if     (abs((plat(jj)-plat(1))-
     &                  (plat( 2)-plat(1))*(jj-1)).lt.1.e-2) then
              call nchek("nf90_put_att-point_spacing",
     &                    nf90_put_att(ncfileID,pLatVarID,
     &                                 "point_spacing","even"))  !ferret
            endif
            call nchek("nf90_put_att-axis",
     &                  nf90_put_att(ncfileID,pLatVarID,
     &                               "axis","Y"))
              call nchek("nf90_def_var-Longitude",
     &                    nf90_def_var(ncfileID,"Longitude", nf90_float,
     &                                 pLonDimID,pLonVarID))
            call nchek("nf90_put_att-standard_name",
     &                  nf90_put_att(ncfileID,pLonVarID,
     &                               "standard_name","longitude"))
            call nchek("nf90_put_att-units",
     &                  nf90_put_att(ncfileID,pLonVarID,
     &                               "units","degrees_east"))
            if     (abs((plon(ii)-plon(1))-
     &                  (plon( 2)-plon(1))*(ii-1)).lt.1.e-2) then
              call nchek("nf90_put_att-point_spacing",
     &                    nf90_put_att(ncfileID,pLonVarID,
     &                                 "point_spacing","even"))  !ferret
            endif
            if     (abs((plon(ii)+(plon(2)-plon(1)))-
     &                  (plon( 1)+ 360.0) ).lt.1.e-2) then
              call nchek("nf90_put_att-modulo",
     &                    nf90_put_att(ncfileID,pLonVarID,
     &                                 "modulo","360 degrees"))  !ferret
            endif
            call nchek("nf90_put_att-axis",
     &                  nf90_put_att(ncfileID,pLonVarID,
     &                               "axis","X"))
          if     (kk.le.0) then
            ! model 2d variable
              call nchek("nf90_def_var-namec",
     &                    nf90_def_var(ncfileID,trim(namec),nf90_float,
     &                               (/pLonDimID, pLatDimID, MTDimID/),
     &                                 varID))
          else
            ! model 3d variable
              call nchek("nf90_def_var-namec",
     &                    nf90_def_var(ncfileID,trim(namec),nf90_float,
     &                     (/pLonDimID, pLatDimID, LayDimID, MTDimID/),
     &                                 varID))
          endif
            if     (yrflag.ge.0) then
              call nchek("nf90_put_att-coordinates",
     &                    nf90_put_att(ncfileID,varID,
     &                                 "coordinates",
     &                                 "Date"))
            endif
          if     (name .ne." ") then
            call nchek("nf90_put_att-long_name",
     &                  nf90_put_att(ncfileID,varID,
     &                               "long_name",trim(name)))
          endif
          if     (names.ne." ") then
            call nchek("nf90_put_att-standard_name",
     &                  nf90_put_att(ncfileID,varID,
     &                               "standard_name",trim(names)))
          endif
          call nchek("nf90_put_att-units",
     &                nf90_put_att(ncfileID,varID,"units",trim(units)))
          call nchek("nf90_put_att-_FillValue",
     &                nf90_put_att(ncfileID,varID,
     &                             "_FillValue",fill_value))
          call nchek("nf90_put_att-valid_range",
     &                nf90_put_att(ncfileID,varID,
     &                             "valid_range",
     &                             (/hmin, hmax/)))
          ! leave def mode
          call nchek("nf90_enddef",
     &                nf90_enddef(ncfileID))
          ! write data into coordinate variables
            call nchek("nf90_put_var-time",
     &                  nf90_put_var(ncfileID,MTVarID, time))
            if     (yrflag.ge.0) then
              call nchek("nf90_put_var-date",
     &                    nf90_put_var(ncfileID,datVarID,date))
            endif
            if     (kk.gt.0 .and. zz(1).eq.-1.0) then
              call nchek("nf90_put_var-LayVarID",
     &                    nf90_put_var(ncfileID,LayVarID,
     &                                 (/(k, k=1,kk)/)))
            elseif (kk.gt.0) then !z-depths
              call nchek("nf90_put_var-LayVarID",
     &                    nf90_put_var(ncfileID,LayVarID,zz(1:kk)))
            endif
            call nchek("nf90_put_var-pLatVarID",
     &                  nf90_put_var(ncfileID,pLatVarID,
     &                               plat(:)))     !1-d Latitudes
            call nchek("nf90_put_var-pLonVarID",
     &                  nf90_put_var(ncfileID,pLonVarID,
     &                               plon(:)))     !1-d Longtudes
          ! write to model variable
          if     (kk.le.0) then
            call nchek("nf90_put_var-array",
     &                  nf90_put_var(ncfileID,varID,array(:,:)))
          else
            call nchek("nf90_put_var-array",
     &                  nf90_put_var(ncfileID,varID,array(:,:)))
          endif
          ! close NetCDF file
          call nchek("nf90_close",
     &                nf90_close(ncfileID))
        endif !lexist
        return  !from first call
      endif  !initialization
c
c     Append data to the NetCDF file
c
      kk_rec = kk_rec + 1
      if     (kk_rec.gt.kk) then  !always .true. for k==0.
        kk_rec = 1
        mt_rec = mt_rec + 1
c
        if     (wdayi.ne.0.0) then 
          time = wday1 + (mt_rec-1)*wdayi
        else
          read(22,*) time
        endif
        if     (yrflag.ge.0) then
          call fordate(time,yrflag, iyear,month,iday,ihour)
          date = (iday + 100 * month + 10000 * iyear) + ihour/24.d0
        else
          date = mt_rec
        endif
      endif  !new time
c
      write(6,6300) kk_rec,mt_rec,time,date
c
      ! open NetCDF file
      call nchek("nf90_open",
     &            nf90_open(trim(ncfile),nf90_write, ncfileID))
      !append values
      if     (kk_rec.eq.1) then
        call nchek("nf90_put_var-time",
     &              nf90_put_var(ncfileID,MTVarID, time,
     &                           start=(/mt_rec/)))
        call nchek("nf90_put_var-date",
     &              nf90_put_var(ncfileID,datVarID,date,
     &                           start=(/mt_rec/)))
      endif
      if     (kk.le.0) then
        call nchek("nf90_put_var-array",
     &              nf90_put_var(ncfileID,varID,array(:,:),
     &                           start=(/1,1,mt_rec/)))
      else
        call nchek("nf90_put_var-array",
     &              nf90_put_var(ncfileID,varID,array(:,:),
     &                           start=(/1,1,kk_rec,mt_rec/)))
      endif
      !update valid_range
      call ncrange(array,ii,jj,1, fill_value, hmin,hmax)
      call nchek("nf90_get_att-valid_range",
     &            nf90_get_att(ncfileID,varID,
     &                         "valid_range",
     &                         hrange(1:2)))
      hrange(1) = min( hrange(1), hmin )
      hrange(2) = max( hrange(2), hmax )
      call nchek("nf90_put_att-valid_range",
     &            nf90_put_att(ncfileID,varID,
     &                         "valid_range",
     &                         hrange(1:2)))
      ! close file 
      call nchek("nf90_close",
     &            nf90_close(ncfileID))
      return
 6300 FORMAT(10X,'WRITING RECORD:',
     +                 '  K =',I5,
     +                '  MT =',I5,
     +           '     FDAY =',F9.2,
     +            '   FDATE =',F13.4 )
      end

      subroutine nchek(cnf90,status)
      use netcdf   ! NetCDF fortran 90 interface
      implicit none
c
      character*(*), intent(in) :: cnf90
      integer,       intent(in) :: status
c
c     subroutine to handle NetCDF errors
c
      if     (.FALSE.) then !nodebug
*     if     (.TRUE. ) then !debug
        write(6,'(a)') trim(cnf90)
      endif

      if (status /= nf90_noerr) then
        write(6,'(/a)')   'error in profout - from NetCDF library'
        write(6,'(a/)')   trim(cnf90)
        write(6,'(a/)')   trim(nf90_strerror(status))
        stop
      end if
      end subroutine nchek

      subroutine ncrange(h,ii,jj,kk, fill_value, hmin,hmax)
      implicit none
c
      integer, intent(in ) :: ii,jj,kk
      real,    intent(in ) :: h(ii,jj,kk),fill_value
      real,    intent(out) :: hmin,hmax
c
c     return range of array, ignoring fill_value
c
      integer i,j,k
      real    hhmin,hhmax
c
      hhmin =  abs(fill_value)
      hhmax = -abs(fill_value)
      do k= 1,kk
        do j= 1,jj
          do i= 1,ii
            if     (h(i,j,k).ne.fill_value) then
              hhmin = min(hhmin,h(i,j,k))
              hhmax = max(hhmax,h(i,j,k))
            endif
          enddo
        enddo
      enddo
      hmin = hhmin
      hmax = hhmax
      end subroutine ncrange

      subroutine fordate(dtime,yrflag, iyear,month,iday,ihour)
      implicit none
c
      double precision dtime
      integer          yrflag, iyear,month,iday,ihour
c
c --- converts model day to "calendar" date (year,month,day,hour).
c
      integer          jday,k,m
c
      integer month0(13,3)
      data month0 / 1,  31,  61,  91, 121, 151, 181,
     +                 211, 241, 271, 301, 331, 361,
     +              1,  32,  60,  91, 121, 152, 182,
     +                 213, 244, 274, 305, 335, 366,
     +              1,  32,  61,  92, 122, 153, 183,
     +                 214, 245, 275, 306, 336, 367 /
c
      call forday(dtime,yrflag, iyear,jday,ihour)
c
      if (yrflag.eq.3) then
        if     (mod(iyear,4).eq.0) then
          k = 3
        else
          k = 2
        endif
      elseif (yrflag.eq.0) then
        k = 1
      elseif (yrflag.eq.4) then
        k = 2
      else
        k = 3
      endif
      do m= 1,12
        if     (jday.ge.month0(m,  k) .and.
     +          jday.lt.month0(m+1,k)      ) then
          month = m
          iday  = jday - month0(m,k) + 1
        endif
      enddo
      return
      end
      subroutine forday(dtime,yrflag, iyear,iday,ihour)
      implicit none
c
      double precision dtime
      integer          yrflag, iyear,iday,ihour
c
c --- converts model day to "calendar" date (year,julian-day,hour).
c
      integer        lp
      common/linepr/ lp
      save  /linepr/
c
      double precision dtim1,day
      integer          iyr,nleap
c
      if     (yrflag.eq.0) then
c ---   360 days per model year, starting Jan 16
        iyear =  int((dtime+15.001d0)/360.d0) + 1
        iday  =  mod( dtime+15.001d0 ,360.d0) + 1
        ihour = (mod( dtime+15.001d0 ,360.d0) + 1.d0 - iday)*24.d0
c
      elseif (yrflag.eq.1) then
c ---   366 days per model year, starting Jan 16
        iyear =  int((dtime+15.001d0)/366.d0) + 1
        iday  =  mod( dtime+15.001d0 ,366.d0) + 1
        ihour = (mod( dtime+15.001d0 ,366.d0) + 1.d0 - iday)*24.d0
c
      elseif (yrflag.eq.2) then
c ---   366 days per model year, starting Jan 01
        iyear =  int((dtime+ 0.001d0)/366.d0) + 1
        iday  =  mod( dtime+ 0.001d0 ,366.d0) + 1
        ihour = (mod( dtime+ 0.001d0 ,366.d0) + 1.d0 - iday)*24.d0
c
      elseif (yrflag.eq.4) then
c ---   365 days per model year, starting Jan 01
        iyear =  int((dtime+ 0.001d0)/365.d0) + 1
        iday  =  mod( dtime+ 0.001d0 ,365.d0) + 1
        ihour = (mod( dtime+ 0.001d0 ,365.d0) + 1.d0 - iday)*24.d0
c
      elseif (yrflag.eq.3) then
c ---   model day is calendar days since 01/01/1901
        iyr   = (dtime-1.d0)/365.25d0
        nleap = iyr/4
        dtim1 = 365.d0*iyr + nleap + 1.d0
        day   = dtime - dtim1 + 1.d0
        if     (dtim1.gt.dtime) then
          iyr = iyr - 1
        elseif (day.ge.367.d0) then
          iyr = iyr + 1
        elseif (day.ge.366.d0 .and. mod(iyr,4).ne.3) then
          iyr = iyr + 1
        endif
        nleap = iyr/4
        dtim1 = 365.d0*iyr + nleap + 1.d0
c
        iyear =  1901 + iyr
        iday  =  dtime - dtim1 + 1
        ihour = (dtime - dtim1 + 1.d0 - iday)*24.d0
c
      else
        write(lp,*)
        write(lp,*) 'error in forday - unsupported yrflag value'
        write(lp,*)
*       call flush(lp)
        stop '(forday)'
      endif
      return
      end
