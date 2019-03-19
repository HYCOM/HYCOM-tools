      PROGRAM HYCOM_SCRIP_NC
      IMPLICIT NONE
C
C  hycom_scrip_nc - Usage:  hycom_scrip_nc [[grid.a] mask.a]
C
C                 creates a netDCF file describing the input region
C                 suitable for use with SCRIP and ESMF REGRID.
C                 http://climate.lanl.gov/Software/SCRIP/
C
C                 NetCDF environment variables:
C                    CDF_FILE   ndcf filename
C                    CDF_TITLE  title
C
C                 grid.a is a hycom grid file, default regional.grid.a.
C                 Note that the corresponding grid.b must also exist.
C
C                 mask.a contains an optional mask array, that turns
C                 off sampling anywhere it is 2.0**100 (data void).
C                 Note that grid.a must be specified if mask.a is used.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  June 2010.
C
      REAL*4,  ALLOCATABLE          :: A(:,:)
      REAL*4                        :: PAD(4096)
      INTEGER, ALLOCATABLE          :: PMSK(:)
      DOUBLE PRECISION, ALLOCATABLE :: PLAT(:),PLON(:),
     +                                 Q4LAT(:,:),Q4LON(:,:)
C
      INTEGER       IOS
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      INTEGER       IDM,JDM,NPAD
      CHARACTER*240 CFILEB,CFILEG,CFILEM,CVARIN
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.0) THEN
        CFILEG = 'regional.grid.a'
        CFILEM = 'NONE'
      ELSEIF (NARG.EQ.1) THEN
        CALL GETARG(1,CFILEG)
        CFILEM = 'NONE'
      ELSEIF (NARG.EQ.2) THEN
        CALL GETARG(1,CFILEG)
        CALL GETARG(2,CFILEM)
      ELSE
        WRITE(6,*) 
     +   'Usage:  hycom_scrip_nc [[grid.a [mask.a]'
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
        WRITE(6,*) 'hycom_scrip_nc: bad header file ',
     &             TRIM(CFILEB)
        CALL EXIT(2)
      ENDIF
      READ( 11,*) JDM,CVARIN
      IF (CVARIN.NE.'jdm   ') THEN
        WRITE(6,*) 'hycom_scrip_nc: bad header file ',
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
        WRITE(6,*) 'Error in hycom_scrip_nc: could not allocate ',
     +             IDM*JDM,' words for A'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( PMSK(IDM*JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_scrip_nc: could not allocate ',
     +             IDM*JDM,' words for PMSK'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( PLAT(IDM*JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_scrip_nc: could not allocate ',
     +             IDM*JDM,' r*8 words for PLAT'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( PLON(IDM*JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_scrip_nc: could not allocate ',
     +             IDM*JDM,' r*8 words for PLON'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( Q4LAT(4,IDM*JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_scrip_nc: could not allocate ',
     +             4*IDM*JDM,' r*8 words for Q4LAT'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( Q4LON(4,IDM*JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_scrip_nc: could not allocate ',
     +             4*IDM*JDM,' r*8 words for Q4LON'
        CALL EXIT(2)
      ENDIF
C
      CALL SCRIP_NC(A,PMSK,PLAT,PLON,Q4LAT,Q4LON,IDM,JDM,PAD,NPAD,
     +              CFILEG,CFILEM)
      CALL EXIT(0)
 5000 FORMAT(I4)
      END
      SUBROUTINE SCRIP_NC(A,PMSK,PLAT,PLON,Q4LAT,Q4LON,IDM,JDM,
     +                    PAD,NPAD, CFILEG,CFILEM)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240    CFILEG,CFILEM
      INTEGER          IDM,JDM,NPAD
      REAL*4           A(IDM,JDM),PAD(NPAD)
      INTEGER           PMSK(  IDM*JDM)
      DOUBLE PRECISION  PLAT(  IDM*JDM), PLON(  IDM*JDM),
     +                 Q4LAT(4,IDM*JDM),Q4LON(4,IDM*JDM)
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      CHARACTER*240 CNAMES,CFILES
      CHARACTER*18  CASN
      LOGICAL       LPERIOD,LARCTIC
      INTEGER       I,II,IP1,J,JJ,IOS,NRECL
      REAL*8        DEG2RAD
#ifdef CRAY
      INTEGER*8     IU8,IOS8
#endif
C
      DEG2RAD = 4.D0*ATAN(1.D0)/180.D0  !PI/180
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
        PMSK(:) = 1  !all sea points
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
        READ(11,REC=1,IOSTAT=IOS) A
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'can''t read ',TRIM(CFILEM)
          CALL EXIT(4)
        ENDIF
C
        CLOSE(UNIT=11)
C
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (A(I,J).NE.SPVAL) THEN
              PMSK(I+(J-1)*IDM) = 1  ! sea point
            ELSE
              PMSK(I+(J-1)*IDM) = 0  !land point
            ENDIF
          ENDDO !i
        ENDDO !j
      ENDIF
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
      READ(11,REC=1,IOSTAT=IOS) A  ! plon
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read ',TRIM(CFILEG)
        CALL EXIT(4)
      ENDIF
      DO J= 1,JDM
        DO I= 1,IDM
          PLON(I+(J-1)*IDM) = A(I,J)*DEG2RAD
        ENDDO !i
      ENDDO !j
C --- LPERIOD and LARCTIC are guesses, LARCTIC is really .not. recti-linear
      LPERIOD = ABS(A(IDM,1)-A(1,1)).GT.350.0
      LARCTIC = LPERIOD .AND. (A(1,    1).NE.A(1,    JDM) .OR.
     &                         A(IDM/2,1).NE.A(IDM/2,JDM)     )
C
      WRITE(6,*) 
      WRITE(6,*) 'LPERIOD = ',LPERIOD
      WRITE(6,*) 'LARCTIC = ',LARCTIC
      WRITE(6,*) 
C
      READ(11,REC= 2,IOSTAT=IOS) A  ! plat
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read ',TRIM(CFILEG)
        CALL EXIT(4)
      ENDIF
      DO J= 1,JDM
        DO I= 1,IDM
          PLAT(I+(J-1)*IDM) = A(I,J)*DEG2RAD
        ENDDO !i
      ENDDO !j
C
      READ(11,REC=3,IOSTAT=IOS) A  ! qlon
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read ',TRIM(CFILEG)
        CALL EXIT(4)
      ENDIF
      DO J= 1,JDM
        DO I= 1,IDM
          Q4LON(1,I+(J-1)*IDM) = A(I,J)*DEG2RAD
        ENDDO !i
      ENDDO !j
C
      READ(11,REC= 4,IOSTAT=IOS) A  ! qlat
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read ',TRIM(CFILEG)
        CALL EXIT(4)
      ENDIF
      DO J= 1,JDM
        DO I= 1,IDM
          Q4LAT(1,I+(J-1)*IDM) = A(I,J)*DEG2RAD
        ENDDO !i
      ENDDO !j
C
      CLOSE(UNIT=11)
C
C     CORNERS IN COUNTER-CLOCKWISE ORDER.
C
      IF     (LPERIOD) THEN
        DO J= 1,JDM-1
          DO I= 1,IDM
            IP1 = MOD(I,IDM) + 1
            Q4LAT(2,I+(J-1)*IDM) = Q4LAT(1,IP1+(J-1)*IDM) !i+1,j
            Q4LAT(3,I+(J-1)*IDM) = Q4LAT(1,IP1+ J   *IDM) !i+1,j+1
            Q4LAT(4,I+(J-1)*IDM) = Q4LAT(1,I  + J   *IDM) !i,  j+1
            Q4LON(2,I+(J-1)*IDM) = Q4LON(1,IP1+(J-1)*IDM) !i+1,j
            Q4LON(3,I+(J-1)*IDM) = Q4LON(1,IP1+ J   *IDM) !i+1,j+1
            Q4LON(4,I+(J-1)*IDM) = Q4LON(1,I  + J   *IDM) !i,  j+1
          ENDDO !i
        ENDDO !j
        IF     (LARCTIC) THEN
C ---     NO NEED FOR J=JDM.
        ELSE
          J=JDM
          DO I= 1,IDM
            IF     (ABS(PLAT(I)+
     &                  PLAT(I+(J-1)*IDM)).LT.0.01*DEG2RAD) THEN
C ---         ANTI-SYMETRIC ABOUT THE EQUATOR.
              Q4LAT(2,I+(J-1)*IDM) = -Q4LAT(3,I)
              Q4LAT(3,I+(J-1)*IDM) = -Q4LAT(2,I)
              Q4LAT(4,I+(J-1)*IDM) = -Q4LAT(1,I)
            ELSE
              Q4LAT(2,I+(J-1)*IDM) = 2.0*Q4LAT(2,I+(J-2)*IDM) -
     &                                   Q4LAT(2,I+(J-3)*IDM)
              Q4LAT(3,I+(J-1)*IDM) = 2.0*Q4LAT(3,I+(J-2)*IDM) -
     &                                   Q4LAT(3,I+(J-3)*IDM)
              Q4LAT(4,I+(J-1)*IDM) = 2.0*Q4LAT(4,I+(J-2)*IDM) -
     &                                   Q4LAT(4,I+(J-3)*IDM)
            ENDIF
              Q4LON(2,I+(J-1)*IDM) = 2.0*Q4LON(2,I+(J-2)*IDM) -
     &                                   Q4LON(2,I+(J-3)*IDM)
              Q4LON(3,I+(J-1)*IDM) = 2.0*Q4LON(3,I+(J-2)*IDM) -
     &                                   Q4LON(3,I+(J-3)*IDM)
              Q4LON(4,I+(J-1)*IDM) = 2.0*Q4LON(4,I+(J-2)*IDM) -
     &                                   Q4LON(4,I+(J-3)*IDM)
          ENDDO !i
        ENDIF
      ELSE
        DO J= 1,JDM-1
          DO I= 1,IDM-1
            Q4LAT(2,I+(J-1)*IDM) = Q4LAT(1,I+1+(J-1)*IDM) !i+1,j
            Q4LAT(3,I+(J-1)*IDM) = Q4LAT(1,I+1+ J   *IDM) !i+1,j+1
            Q4LAT(4,I+(J-1)*IDM) = Q4LAT(1,I  + J   *IDM) !i,  j+1
            Q4LON(2,I+(J-1)*IDM) = Q4LON(1,I+1+(J-1)*IDM) !i+1,j
            Q4LON(3,I+(J-1)*IDM) = Q4LON(1,I+1+ J   *IDM) !i+1,j+1
            Q4LON(4,I+(J-1)*IDM) = Q4LON(1,I  + J   *IDM) !i,  j+1
          ENDDO !i
            I=IDM
            Q4LAT(2,I+(J-1)*IDM) = 2.0*Q4LAT(2,I-1+(J-1)*IDM) -
     &                                 Q4LAT(2,I-2+(J-1)*IDM)
            Q4LAT(3,I+(J-1)*IDM) = 2.0*Q4LAT(3,I-1+(J-1)*IDM) -
     &                                 Q4LAT(3,I-2+(J-1)*IDM)
            Q4LAT(4,I+(J-1)*IDM) = 2.0*Q4LAT(4,I-1+(J-1)*IDM) -
     &                                 Q4LAT(4,I-2+(J-1)*IDM)
            Q4LON(2,I+(J-1)*IDM) = 2.0*Q4LON(2,I-1+(J-1)*IDM) -
     &                                 Q4LON(2,I-2+(J-1)*IDM)
            Q4LON(3,I+(J-1)*IDM) = 2.0*Q4LON(3,I-1+(J-1)*IDM) -
     &                                 Q4LON(3,I-2+(J-1)*IDM)
            Q4LON(4,I+(J-1)*IDM) = 2.0*Q4LON(4,I-1+(J-1)*IDM) -
     &                                 Q4LON(4,I-2+(J-1)*IDM)
        ENDDO !j
          J=JDM
          DO I= 1,IDM
            Q4LAT(2,I+(J-1)*IDM) = 2.0*Q4LAT(2,I+(J-2)*IDM) -
     &                                 Q4LAT(2,I+(J-3)*IDM)
            Q4LAT(3,I+(J-1)*IDM) = 2.0*Q4LAT(3,I+(J-2)*IDM) -
     &                                 Q4LAT(3,I+(J-3)*IDM)
            Q4LAT(4,I+(J-1)*IDM) = 2.0*Q4LAT(4,I+(J-2)*IDM) -
     &                                 Q4LAT(4,I+(J-3)*IDM)
            Q4LON(2,I+(J-1)*IDM) = 2.0*Q4LON(2,I+(J-2)*IDM) -
     &                                 Q4LON(2,I+(J-3)*IDM)
            Q4LON(3,I+(J-1)*IDM) = 2.0*Q4LON(3,I+(J-2)*IDM) -
     &                                 Q4LON(3,I+(J-3)*IDM)
            Q4LON(4,I+(J-1)*IDM) = 2.0*Q4LON(4,I+(J-2)*IDM) -
     &                                 Q4LON(4,I+(J-3)*IDM)
          ENDDO !i
      ENDIF
C
C
C
      WRITE(6,*)
      WRITE(6,*) 'plon(  1,  1)   = ', plon(  1)      /deg2rad
      WRITE(6,*) 'plat(  1,  1)   = ', plat(  1)      /deg2rad
      WRITE(6,*) 'qlon(  1,  1).1 = ',q4lon(1,1)      /deg2rad
      WRITE(6,*) 'qlat(  1,  1).1 = ',q4lat(1,1)      /deg2rad
      WRITE(6,*) 'qlon(  1,  1).2 = ',q4lon(2,1)      /deg2rad
      WRITE(6,*) 'qlat(  1,  1).2 = ',q4lat(2,1)      /deg2rad
      WRITE(6,*) 'qlon(  1,  1).3 = ',q4lon(3,1)      /deg2rad
      WRITE(6,*) 'qlat(  1,  1).3 = ',q4lat(3,1)      /deg2rad
      WRITE(6,*) 'qlon(  1,  1).4 = ',q4lon(4,1)      /deg2rad
      WRITE(6,*) 'qlat(  1,  1).4 = ',q4lat(4,1)      /deg2rad
      WRITE(6,*)
      WRITE(6,*) 'plon(idm,  1)   = ', plon(  idm)    /deg2rad
      WRITE(6,*) 'plat(idm,  1)   = ', plat(  idm)    /deg2rad
      WRITE(6,*) 'qlon(idm,  1).1 = ',q4lon(1,idm)    /deg2rad
      WRITE(6,*) 'qlat(idm,  1).1 = ',q4lat(1,idm)    /deg2rad
      WRITE(6,*) 'qlon(idm,  1).2 = ',q4lon(2,idm)    /deg2rad
      WRITE(6,*) 'qlat(idm,  1).2 = ',q4lat(2,idm)    /deg2rad
      WRITE(6,*) 'qlon(idm,  1).3 = ',q4lon(3,idm)    /deg2rad
      WRITE(6,*) 'qlat(idm,  1).3 = ',q4lat(3,idm)    /deg2rad
      WRITE(6,*) 'qlon(idm,  1).4 = ',q4lon(4,idm)    /deg2rad
      WRITE(6,*) 'qlat(idm,  1).4 = ',q4lat(4,idm)    /deg2rad
      WRITE(6,*)
      IF     (LARCTIC) THEN
        JJ= JDM-1
      ELSE
        JJ= JDM
      ENDIF
      WRITE(6,*) 'plon(idm, jj)   = ', plon(  idm* jj)/deg2rad
      WRITE(6,*) 'plat(idm, jj)   = ', plat(  idm* jj)/deg2rad
      WRITE(6,*) 'qlon(idm, jj).1 = ',q4lon(1,idm* jj)/deg2rad
      WRITE(6,*) 'qlat(idm, jj).1 = ',q4lat(1,idm* jj)/deg2rad
      WRITE(6,*) 'qlon(idm, jj).2 = ',q4lon(2,idm* jj)/deg2rad
      WRITE(6,*) 'qlat(idm, jj).2 = ',q4lat(2,idm* jj)/deg2rad
      WRITE(6,*) 'qlon(idm, jj).3 = ',q4lon(3,idm* jj)/deg2rad
      WRITE(6,*) 'qlat(idm, jj).3 = ',q4lat(3,idm* jj)/deg2rad
      WRITE(6,*) 'qlon(idm, jj).4 = ',q4lon(4,idm* jj)/deg2rad
      WRITE(6,*) 'qlat(idm, jj).4 = ',q4lat(4,idm* jj)/deg2rad
      WRITE(6,*)
C
C     CREATE NCDF FILE.
C
      CALL SCRIP_OUT(PMSK,PLON,PLAT,Q4LON,Q4LAT,IDM,JJ)
      END

      subroutine scrip_out(pmsk,plon,plat,q4lon,q4lat,ii,jj)
      use netcdf   ! NetCDF fortran 90 interface
      implicit none
c
      integer          ii,jj
      integer           pmsk(  ii*jj)
      double precision  plat(  ii*jj), plon(  ii*jj),
     &                 q4lat(4,ii*jj),q4lon(4,ii*jj)
c
c     NetCDF environment variables:
c       CDF_FILE  ndcf filename
c       CDF_TITLE title
c
c     This routine needs version 3.5 of the NetCDF library, from: 
c     http://www.unidata.ucar.edu/packages/netcdf/
c
      integer          :: ncfileID, status, varID
      integer          :: GSDimID,GRDimID,GCDimID
      character*240    :: ncfile,ncenv
      character*240    :: name,namec,names,units
c
      logical          :: lopen,lexist
      integer          :: grid_size, grid_rank, grid_corners,
     &                    grid_dims(2)
c
c       initialization.
c
        write(6,*) 
        write(6,*) 'scrip_out - plon = ',plon(1)
        write(6,*) 'scrip_out - plat = ',plat(1)
        write(6,*) 
c
        ncfile = ' '
        call getenv('CDF_FILE',ncfile)
        if     (ncfile.eq.' ') then
          write( 6,'(/2a/)')  'error in scrip_out - ',
     &                        'CDF_FILE not defined'
          stop
        endif
c
        inquire(file= ncfile, exist=lexist)
        if (lexist) then
          write( 6,'(/2a/a/)') 'error in scrip_out - ',
     &                        'CDF_FILE is an existing file',
     &                        trim(ncfile)
          stop
        endif
c
c       create a new NetCDF and write data to it
c       netcdf-4 classic model, netcdf version 4.3 and later
c
        call nchek('nf90_create',
     &              nf90_create(trim(ncfile),
     &                          or(nf90_clobber,
     &                             or(nf90_hdf5,
     &                                nf90_classic_model)),
     &                          ncfileID))
c       define the dimensions
c
        grid_size    = ii*jj
        grid_rank    = 2
        grid_corners = 4
c
        call nchek("nf90_def_dim-grid_size",
     &              nf90_def_dim(ncfileID,
     &                           "grid_size",
     &                            grid_size,
     &                            GSDimID))
        call nchek("nf90_def_dim-grid_size",
     &              nf90_def_dim(ncfileID,
     &                           "grid_rank",
     &                            grid_rank,
     &                            GRDimID))
        call nchek("nf90_def_dim-grid_size",
     &              nf90_def_dim(ncfileID,
     &                           "grid_corners",
     &                            grid_corners,
     &                            GCDimID))
          ! create the global attributes
            ncenv = ' '
            call getenv('CDF_TITLE',ncenv)
            if     (ncenv.eq.' ') then
              write(ncenv,'(i5,a,i5,a)') ii,' by',jj,' HYCOM region'
            endif
            call nchek("nf90_put_att-title",
     &                  nf90_put_att(ncfileID,nf90_global,
     &                               "title",
     &                               trim(ncenv)))
            ncenv = ' '
            call nchek("nf90_put_att-history",
     &                  nf90_put_att(ncfileID,nf90_global,
     &                               "history",
     &                               "hycom_scrip_nc"))
          ! leave def mode
          call nchek("nf90_enddef",
     &                nf90_enddef(ncfileID))
c
c ---     grid_dims
c
          call nchek("nf90_redef",
     &                nf90_redef(ncfileID))
          call nchek("nf90_def_var-grid_dims",
     &                nf90_def_var(ncfileID,"grid_dims",nf90_int,
     &                             GRDimID,
     &                             varID))
          call nchek("nf90_enddef",
     &                nf90_enddef(ncfileID))
          grid_dims(1) = ii
          grid_dims(2) = jj
          write(6,*) 'grid_dims = ',grid_dims(:)
          call nchek("nf90_put_var-grid_dims",
     &                nf90_put_var(ncfileID,varID,grid_dims(:)))
c
c ---     grid_center_lat
c
          call nchek("nf90_redef",
     &                nf90_redef(ncfileID))
          call nchek("nf90_def_var-grid_center_lat",
     &                nf90_def_var(ncfileID,"grid_center_lat",
     &                             nf90_double,
     &                             GSDimID,
     &                             varID))
          call nchek("nf90_put_att-units",
     &                nf90_put_att(ncfileID,varID,
     &                             "units","radians"))
          call nchek("nf90_enddef",
     &                nf90_enddef(ncfileID))
          call nchek("nf90_put_var-grid_center_lat",
     &                nf90_put_var(ncfileID,varID,plat(:)))
c
c ---     grid_center_lon
c
          call nchek("nf90_redef",
     &                nf90_redef(ncfileID))
          call nchek("nf90_def_var-grid_center_lon",
     &                nf90_def_var(ncfileID,"grid_center_lon",
     &                             nf90_double,
     &                             GSDimID,
     &                             varID))
          call nchek("nf90_put_att-units",
     &                nf90_put_att(ncfileID,varID,
     &                             "units","radians"))
          call nchek("nf90_enddef",
     &                nf90_enddef(ncfileID))
          call nchek("nf90_put_var-grid_center_lon",
     &                nf90_put_var(ncfileID,varID,plon(:)))
c
c ---     grid_imask
c
          call nchek("nf90_redef",
     &                nf90_redef(ncfileID))
          call nchek("nf90_def_var-grid_imask",
     &                nf90_def_var(ncfileID,"grid_imask",
     &                             nf90_int,
     &                             GSDimID,
     &                             varID))
          call nchek("nf90_put_att-units",
     &                nf90_put_att(ncfileID,varID,
     &                             "units","unitless"))
          call nchek("nf90_enddef",
     &                nf90_enddef(ncfileID))
          call nchek("nf90_put_var-grid_imask",
     &                nf90_put_var(ncfileID,varID,pmsk(:)))
c
c ---     grid_corner_lat
c
          call nchek("nf90_redef",
     &                nf90_redef(ncfileID))
          call nchek("nf90_def_var-grid_corner_lat",
     &                nf90_def_var(ncfileID,"grid_corner_lat",
     &                             nf90_double,
     &                             (/GCDimID, GSDimID/),
     &                             varID))
          call nchek("nf90_put_att-units",
     &                nf90_put_att(ncfileID,varID,
     &                             "units","radians"))
          call nchek("nf90_enddef",
     &                nf90_enddef(ncfileID))
          call nchek("nf90_put_var-grid_corner_lat",
     &                nf90_put_var(ncfileID,varID,q4lat(:,:)))
c
c ---     grid_corner_lon
c
          call nchek("nf90_redef",
     &                nf90_redef(ncfileID))
          call nchek("nf90_def_var-grid_corner_lon",
     &                nf90_def_var(ncfileID,"grid_corner_lon",
     &                             nf90_double,
     &                             (/GCDimID, GSDimID/),
     &                             varID))
          call nchek("nf90_put_att-units",
     &                nf90_put_att(ncfileID,varID,
     &                             "units","radians"))
          call nchek("nf90_enddef",
     &                nf90_enddef(ncfileID))
          call nchek("nf90_put_var-grid_corner_lon",
     &                nf90_put_var(ncfileID,varID,q4lon(:,:)))
c
          ! close NetCDF file
          call nchek("nf90_close",
     &                nf90_close(ncfileID))
      return
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
*     if     (.TRUE. ) then !debug
      if     (.FALSE.) then !nodebug
        write(6,'(a)') trim(cnf90)
      endif

      if (status /= nf90_noerr) then
        write(6,'(/a)')   'error in profout - from NetCDF library'
        write(6,'(a/)')   trim(cnf90)
        write(6,'(a/)')   trim(nf90_strerror(status))
        stop
      end if
      end subroutine nchek
