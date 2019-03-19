      PROGRAM HYCOM_BINNING
      IMPLICIT NONE
C
C  hycom_binning - Usage:  hycom_binning file.a k cell [lonmin lonmax latmin latmax] [grid.a] [mask.a]
C
C                 prints the mean, min, max of lon/lat bins
C                 from the k-th (1:idm,1:jdm) array in file.a
C                 data voids have mean,min,max=1.e10.
C
C                 longitude bands are centered on 0E and latitude bands
C                 are centered on the equator, both are nominally cell
C                 degrees wide, e.g. (N-0.5)*cell to (N+0.5)*cell.
C                 However, the cell sample width is (approximately) 
C                 scaled by 1/cos(lat) in longitude, i.e. the cell
C                 width in meters is aproximately constant.  The width 
C                 in longitude cells is the last column of the printout.
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
C  For latitude bands, use the similar hycom_zonal_lat command.
C
C  file.a is assumed to contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  January 2001.
C
      REAL*4, ALLOCATABLE :: A(:,:),AMSK(:,:),
     +                       PLAT(:,:),PLON(:,:),PSC2(:,:)
      REAL*4              :: PAD(4096)
      INTEGER IOS
      INTEGER      IARGC
      INTEGER      NARG
      CHARACTER*240 CARG
C
      INTEGER       IDM,JDM,K,NPAD
      REAL*4        CELL,LONMIN,LONMAX,LATMIN,LATMAX
      CHARACTER*6   CVARIN
      CHARACTER*240 CFILE,CFILEG,CFILEM,CFILEB
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.3) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) K
        CALL GETARG(3,CARG)
        READ(CARG,*)    CELL
        CFILEG = 'regional.grid.a'
        CFILEM = 'NONE'
        LATMAX = 999.0  ! indicates entire region
      ELSEIF (NARG.EQ.4) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) K
        CALL GETARG(3,CARG)
        READ(CARG,*)    CELL
        CALL GETARG(4,CFILEG)
        CFILEM = 'NONE'
        LATMAX = 999.0  ! indicates entire region
      ELSEIF (NARG.EQ.5) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) K
        CALL GETARG(3,CARG)
        READ(CARG,*)    CELL
        CALL GETARG(4,CFILEG)
        CALL GETARG(5,CFILEM)
        LATMAX = 999.0  ! indicates entire region
      ELSEIF (NARG.EQ.7) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) K
        CALL GETARG(3,CARG)
        READ(CARG,*)    CELL
        CALL GETARG(4,CARG)
        READ(CARG,*)    LONMIN
        CALL GETARG(5,CARG)
        READ(CARG,*)    LONMAX
        CALL GETARG(6,CARG)
        READ(CARG,*)    LATMIN
        CALL GETARG(7,CARG)
        READ(CARG,*)    LATMAX
        CFILEG = 'regional.grid.a'
        CFILEM = 'NONE'
      ELSEIF (NARG.EQ.8) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) K
        CALL GETARG(3,CARG)
        READ(CARG,*)    CELL
        CALL GETARG(4,CARG)
        READ(CARG,*)    LONMIN
        CALL GETARG(5,CARG)
        READ(CARG,*)    LONMAX
        CALL GETARG(6,CARG)
        READ(CARG,*)    LATMIN
        CALL GETARG(7,CARG)
        READ(CARG,*)    LATMAX
        CALL GETARG(8,CFILEG)
        CFILEM = 'NONE'
      ELSEIF (NARG.EQ.9) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) K
        CALL GETARG(3,CARG)
        READ(CARG,*)    CELL
        CALL GETARG(4,CARG)
        READ(CARG,*)    LONMIN
        CALL GETARG(5,CARG)
        READ(CARG,*)    LONMAX
        CALL GETARG(6,CARG)
        READ(CARG,*)    LATMIN
        CALL GETARG(7,CARG)
        READ(CARG,*)    LATMAX
        CALL GETARG(8,CFILEG)
        CALL GETARG(9,CFILEM)
      ELSE
        WRITE(6,*) 
     +   'Usage:  hycom_binning file.a k cell ' //
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
        WRITE(6,*) 'hycom_binning: bad header file ',
     &             CFILEB(1:LEN_TRIM(CFILEB))
        CALL EXIT(2)
      ENDIF
      READ( 11,*) JDM,CVARIN
      IF (CVARIN.NE.'jdm   ') THEN
        WRITE(6,*) 'hycom_binning: bad header file ',
     &             CFILEB(1:LEN_TRIM(CFILEB))
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
     +             K,CELL,LONMIN,LONMAX,LATMIN,LATMAX,
     +             CFILE,CFILEG,CFILEM)
      CALL EXIT(0)
      END
      SUBROUTINE BINNING(A,AMSK,PLAT,PLON,PSC2,IDM,JDM, PAD,NPAD,
     +                   K,CELL,LONMIN,LONMAX,LATMIN,LATMAX,
     +                   CFILE,CFILEG,CFILEM)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE,CFILEG,CFILEM
      INTEGER      IDM,JDM,NPAD,K
      REAL*4       A(IDM,JDM),AMSK(IDM,JDM),
     +             PLAT(IDM,JDM),PLON(IDM,JDM),PSC2(IDM,JDM),PAD(NPAD),
     +             CELL,LONMIN,LONMAX,LATMIN,LATMAX
C
C     MOST OF WORK IS DONE HERE.
C
      REAL*8,  ALLOCATABLE :: SUMA(:,:),SUMB(:,:)
      REAL*4,  ALLOCATABLE :: MINA(:,:),MAXA(:,:)
      INTEGER, ALLOCATABLE ::   IW(:)
C
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      CHARACTER*18 CASN
      INTEGER      LEN_TRIM
      LOGICAL      LSUBREG,LGLOBAL
      INTEGER      I,IC,IF,IL,IZ,IP,IQ,J,JC,JF,JL,JZ,IOS,NRECL
      REAL*8       SUMAT,SUMBT
      REAL*4       MINAT,MAXAT,XLON,DEG2RAD,MINCOS
#ifdef CRAY
      INTEGER*8    IU8,IOS8
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
      OPEN(UNIT=11, FILE=CFILE, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',CFILE(1:LEN_TRIM(CFILE))
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
      READ(11,REC=K,IOSTAT=IOS) A
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read ',CFILE(1:LEN_TRIM(CFILE))
        CALL EXIT(4)
      ENDIF
C
      CLOSE(UNIT=11)
C
C     INPUT MASK ARRAY.
C
      IF     (CFILEM.EQ.'NONE') THEN
        AMSK(:,:) = 0.0
      ELSE
        OPEN(UNIT=11, FILE=CFILEM, FORM='UNFORMATTED', STATUS='OLD',
     +           ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
        IF     (IOS.NE.0) THEN
          write(6,*) 'Error: can''t open ',CFILEM(1:LEN_TRIM(CFILEM))
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
          WRITE(6,*) 'can''t read ',CFILEM(1:LEN_TRIM(CFILEM))
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
        write(6,*) 'Error: can''t open ',CFILEG(1:LEN_TRIM(CFILEG))
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
        WRITE(6,*) 'can''t read ',CFILEG(1:LEN_TRIM(CFILEG))
        CALL EXIT(4)
      ENDIF
C
      READ(11,REC=11,IOSTAT=IOS) PSC2  ! pscy
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(PSC2,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read ',CFILEG(1:LEN_TRIM(CFILEG))
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
        WRITE(6,*) 'can''t read ',CFILEG(1:LEN_TRIM(CFILEG))
        CALL EXIT(4)
      ENDIF
C
      READ(11,REC= 2,IOSTAT=IOS) PLAT  ! plat
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(PLAT,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read ',CFILEG(1:LEN_TRIM(CFILEG))
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
          LONMIN  = MINVAL(PLON(:,:))
          LONMAX  = MAXVAL(PLON(:,:))
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
      ALLOCATE( SUMA(IF:IL,JF:JL),
     +          SUMB(IF:IL,JF:JL) )
      ALLOCATE( MINA(IF:IL,JF:JL),
     +          MAXA(IF:IL,JF:JL) )
      ALLOCATE(         IW(JF:JL) )
C
C     PRINT A HEADER.
C
      IF     (.NOT. LSUBREG) THEN ! full region
        WRITE(6,6000)
     +    TRIM(CFILE),
     +    K,CELL,
     +    TRIM(CFILEG),
     +    TRIM(CFILEM)
      ELSE
        WRITE(6,6100)
     +    TRIM(CFILE),
     +    K,CELL,
     +    LONMIN,LONMAX,LATMIN,LATMAX,
     +    TRIM(CFILEG),
     +    TRIM(CFILEM)
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
C     DO ALL BINS.
C
      SUMAT     = 0.0
      SUMBT     = 0.0
      SUMA(:,:) = 0.0
      SUMB(:,:) = 0.0
      MINA(:,:) =  HUGE(MINA(IF,JF))
      MAXA(:,:) = -HUGE(MAXA(IF,JF))
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (   A(I,J).NE.SPVAL .AND.
     +            AMSK(I,J).NE.SPVAL      ) THEN
            IZ = PLON(I,J)
            JZ = PLAT(I,J)
            IF     (IZ.LT.IF .OR. IZ.GT.IL .OR.
     +              JZ.LT.JF .OR. JZ.GT.JL     ) THEN
              CYCLE !i
            ENDIF
C
            SUMAT = SUMAT + A(I,J)*PSC2(I,J)
            SUMBT = SUMBT +        PSC2(I,J)
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
              MINA(IQ,JZ) = MIN( MINA(IQ,JZ), A(I,J) )
              MAXA(IQ,JZ) = MAX( MAXA(IQ,JZ), A(I,J) )
C
              SUMA(IQ,JZ) = SUMA(IQ,JZ) + A(I,J)*PSC2(I,J)
              SUMB(IQ,JZ) = SUMB(IQ,JZ) +        PSC2(I,J)
*
*             IF     (IQ.EQ.IC .AND. JZ.EQ.JC) THEN
*               WRITE(6,'(a,2i5,1pg16.8)') 'debug: ',I,J,A(I,J)
*             ENDIF
            ENDDO !ip
          ENDIF
        ENDDO !i
      ENDDO !j
C
C     PRINTOUT
C
      DO JZ= JF,JL
        DO IZ= IF,IL
          IF     (SUMB(IZ,JZ).NE.0.0) THEN
            WRITE(6,'(F8.2,F7.2,1p3g16.8,I8)') IZ*CELL,JZ*CELL,
     +                                      SUMA(IZ,JZ)/SUMB(IZ,JZ),
     +                                      MINA(IZ,JZ),MAXA(IZ,JZ),
     +                                      1+2*IW(JZ)
          ELSE
            WRITE(6,'(F8.2,F7.2,1p3g16.8,I8)') IZ*CELL,JZ*CELL,
     +                                      1.0e10,1.0e10,1.0e10,
     +                                      1+2*IW(JZ)
          ENDIF
        ENDDO !iz
      ENDDO !jz
      SUMAT = SUMAT/SUMBT
      MINAT = MINVAL(MINA(:,:))
      MAXAT = MAXVAL(MAXA(:,:))
      WRITE(6,'(a,a/a,1p3g16.8,i8)')
     +   '#    LON    LAT        MEAN             MIN             MAX',
     +   '       WIDTH',
     +   '#    ALL    ALL',SUMAT,MINAT,MAXAT,0
      RETURN
 6000 FORMAT(
     +   '# hycom_binning file.a k cell grid.a mask.a' /
     +   '# hycom_binning ',A       /
     +   '#               ',I5,F7.2 /
     +   '#               ',A       /
     +   '#               ',A       /
     +   '#    LON    LAT        MEAN             MIN             MAX',
     +   '       WIDTH')
 6100 FORMAT(
     +   '# hycom_binning file.a k cell ',
     +   'lonmin lonmax latmin latmax grid.a mask.a' /
     +   '# hycom_binning ',A       /
     +   '#               ',I5,F7.2 /
     +   '#               ',2F9.2,2F8.2 /
     +   '#               ',A       /
     +   '#               ',A       /
     +   '#    LON    LAT        MEAN             MIN             MAX',
     +   '       WIDTH')
      END
