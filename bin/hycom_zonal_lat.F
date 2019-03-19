      PROGRAM FZONAL
      IMPLICIT NONE
C
C  hycom_zonal_lat - Usage:  hycom_zonal_lat file.a k cell [lonmin lonmax latmin latmax] [grid.a] [mask.a]
C
C                 prints the min, max, mean and rms of zonal latitudinal
C                 extents from the k-th (1:idm,1:jdm) array in file.a
C
C                 latitude bands are centered on the equator and are
C                 cell degrees wide, e.g. (N-0.5)*cell to (N+0.5)*cell.
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
C
C  For cell by cell bins, use the similar hycom_binning command.
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
        READ(CARG,5000) K
        CALL GETARG(3,CARG)
        READ(CARG,*)    CELL
        CFILEG = 'regional.grid.a'
        CFILEM = 'NONE'
        LATMAX = 999.0  ! indicates entire region
      ELSEIF (NARG.EQ.4) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,5000) K
        CALL GETARG(3,CARG)
        READ(CARG,*)    CELL
        CALL GETARG(4,CFILEG)
        CFILEM = 'NONE'
        LATMAX = 999.0  ! indicates entire region
      ELSEIF (NARG.EQ.5) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,5000) K
        CALL GETARG(3,CARG)
        READ(CARG,*)    CELL
        CALL GETARG(4,CFILEG)
        CALL GETARG(5,CFILEM)
        LATMAX = 999.0  ! indicates entire region
      ELSEIF (NARG.EQ.7) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,5000) K
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
        READ(CARG,5000) K
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
        READ(CARG,5000) K
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
     +   'Usage:  hycom_zonal_lat file.a k cell ' //
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
        WRITE(6,*) 'hycom_zonal_lat: bad header file ',
     &             CFILEB(1:LEN_TRIM(CFILEB))
        CALL EXIT(2)
      ENDIF
      READ( 11,*) JDM,CVARIN
      IF (CVARIN.NE.'jdm   ') THEN
        WRITE(6,*) 'hycom_zonal_lat: bad header file ',
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
      CALL ZONAL(A,AMSK,PLAT,PLON,PSC2,IDM,JDM,PAD,NPAD,
     +           K,CELL,LONMIN,LONMAX,LATMIN,LATMAX,
     +           CFILE,CFILEG,CFILEM)
      CALL EXIT(0)
 5000 FORMAT(I4)
      END
      SUBROUTINE ZONAL(A,AMSK,PLAT,PLON,PSC2,IDM,JDM, PAD,NPAD,
     +                 K,CELL,LONMIN,LONMAX,LATMIN,LATMAX,
     +                 CFILE,CFILEG,CFILEM)
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
      REAL*8, ALLOCATABLE :: SUMA(:),SUMB(:),SUMR(:)
      REAL*4, ALLOCATABLE :: MINA(:),MAXA(:)
C
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      CHARACTER*18 CASN
      INTEGER      LEN_TRIM
      LOGICAL      OUT_OF_LON_BAND
      LOGICAL      LSUBREG
      INTEGER      I,J,JF,JL,JZ,IOS,NRECL
      REAL*8       SUMAT,SUMBT,SUMRT
      REAL*4       MINAT,MAXAT
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
C     SUBREGION?
C
      LSUBREG = LATMAX.LT.91.0  ! sub-region specified
C
      IF     (LSUBREG) THEN
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
        ELSEIF (LONMAX.GT.LONMIN+360.0) THEN
          DO
            LONMAX = LONMAX - 360.0
            IF     (LONMAX.LE.LONMIN+360.0) THEN
              EXIT
            ENDIF
          ENDDO
        ENDIF !lonmax
        IF     (LONMAX.LT.LONMIN) THEN  !full region
          LONMAX = LONMIN + 360.0
        ENDIF !lonmax
      ENDIF !lsubreg
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
C     ALLOCATE SUM ARRAYS IN CELL SPACE.
C
      ALLOCATE( SUMA(JF:JL),
     +          SUMB(JF:JL),
     +          SUMR(JF:JL) )
      ALLOCATE( MINA(JF:JL),
     +          MAXA(JF:JL) )
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
C     DO ALL LATITUDE BAND SUMS.
C
      SUMAT   = 0.0
      SUMBT   = 0.0
      SUMA(:) = 0.0
      SUMB(:) = 0.0
      MINA(:) =  HUGE(MINA(JF))
      MAXA(:) = -HUGE(MAXA(JF))
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (   A(I,J).NE.SPVAL .AND.
     +            AMSK(I,J).NE.SPVAL      ) THEN
            JZ = PLAT(I,J)
            IF     (JZ.LT.JF .OR. JZ.GT.JL) THEN
              CYCLE
            ENDIF
            IF     (LSUBREG .AND.
     +              OUT_OF_LON_BAND(PLON(I,J), LONMIN,LONMAX)) THEN
              CYCLE
            ENDIF
C
            MINA(JZ) = MIN( MINA(JZ), A(I,J) )
            MAXA(JZ) = MAX( MAXA(JZ), A(I,J) )
C
            SUMA(JZ) = SUMA(JZ) + A(I,J)*PSC2(I,J)
            SUMAT    = SUMAT    + A(I,J)*PSC2(I,J)
            SUMB(JZ) = SUMB(JZ) +        PSC2(I,J)
            SUMBT    = SUMBT    +        PSC2(I,J)
          ENDIF
        ENDDO
      ENDDO
      DO JZ= JF,JL
        IF     (SUMB(JZ).NE.0.0) THEN
          SUMA(JZ) = SUMA(JZ)/SUMB(JZ)
        ENDIF
      ENDDO
      SUMAT = SUMAT/SUMBT
C
C     DO ALL LATITUDE BAND SUM OF SQUARES.
C
      SUMRT   = 0.0
      SUMR(:) = 0.0
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (   A(I,J).NE.SPVAL .AND.
     +            AMSK(I,J).NE.SPVAL      ) THEN
            JZ = PLAT(I,J)
            IF     (JZ.LT.JF .OR. JZ.GT.JL) THEN
              CYCLE
            ENDIF
            IF     (LSUBREG .AND.
     +              OUT_OF_LON_BAND(PLON(I,J), LONMIN,LONMAX)) THEN
              CYCLE
            ENDIF
            SUMR(JZ) = SUMR(JZ) + (A(I,J)-SUMA(JZ))**2 * PSC2(I,J)
            SUMRT    = SUMRT    + (A(I,J)-SUMAT   )**2 * PSC2(I,J)
          ENDIF
        ENDDO
      ENDDO
      DO JZ= JF,JL
        IF     (SUMB(JZ).NE.0.0) THEN
          SUMR(JZ) = SQRT(SUMR(JZ)/SUMB(JZ))
          WRITE(6,'(f7.2,1p4g16.8,2pf9.3)') JZ*CELL,
     +                                      SUMA(JZ),SUMR(JZ),
     +                                      MINA(JZ),MAXA(JZ),
     +                                      SUMB(JZ)/SUMBT
        ELSE
*         WRITE(6,'(f7.2,64x,2pf9.3)')      JZ*CELL,0.0
          WRITE(6,*)
        ENDIF
      ENDDO
      SUMRT = SQRT(SUMRT/SUMBT)
      MINAT = MINVAL(MINA(:))
      MAXAT = MAXVAL(MAXA(:))
      WRITE(6,'(a,a/a,1p4g16.8,2pf9.3)')
     +   '#   LAT        MEAN             RMS',
     +      '             MIN             MAX      PERCENT',
     +   '#   ALL',SUMAT,SUMRT,MINAT,MAXAT,1.0
      RETURN
 6000 FORMAT(
     +   '# hycom_zonal_lat file.a k cell grid.a mask.a' /
     +   '# hycom_zonal_lat ',A       /
     +   '#                 ',I5,F7.2 /
     +   '#                 ',A       /
     +   '#                 ',A       /
     +   '#   LAT        MEAN             RMS',
     +      '             MIN             MAX      PERCENT')
 6100 FORMAT(
     +   '# hycom_zonal_lat file.a k cell ',
     +   'lonmin lonmax latmin latmax grid.a mask.a' /
     +   '# hycom_zonal_lat ',A       /
     +   '#                 ',I5,F7.2 /
     +   '#                 ',2F9.2,2F8.2 /
     +   '#                 ',A       /
     +   '#                 ',A       /
     +   '#   LAT        MEAN             RMS',
     +      '             MIN             MAX      PERCENT')
      END
      LOGICAL FUNCTION OUT_OF_LON_BAND(LON, LONMIN,LONMAX)
      IMPLICIT NONE
      REAL*4 LON, LONMIN,LONMAX
C
C     IS LON OUTSIDE THE RANGE LONMIN TO LONMAX?
C
      REAL*4 XLON
C
      XLON = LON
      IF     (LON.LT.LONMIN) THEN
        DO
          XLON = XLON + 360.0
          IF     (XLON.GT.LONMIN) THEN
            EXIT
          ENDIF
        ENDDO
      ELSEIF (LON.GE.LONMIN+360.0) THEN
        DO
          XLON = XLON - 360.0
          IF     (XLON.LT.LONMIN+360.0) THEN
            EXIT
          ENDIF
        ENDDO
      ENDIF
C
      OUT_OF_LON_BAND = XLON .GT. LONMAX
*
*     IF     (OUT_OF_LON_BAND) THEN
*       write(6,'(a,4f10.3)') 'lon,xlon,min,max = ',
*    &                         LON,XLON,LONMIN,LONMAX
*     ENDIF
      END
