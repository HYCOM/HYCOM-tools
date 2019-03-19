      PROGRAM HYCOM_SCATTER
      IMPLICIT NONE
C
C  hycom_scatter - Usage:  hycom_scatter points.txt field.a dist [grid.a] [mask.a]
C                 
C                 maps the values in points.txt onto a field
C                 the nearest input point is used at each grid point,
C                 but only values within dist m are included
C
C                 points.txt contains value lon lat triplets
C
C                 field.a will contain the generated field
C
C                 grid.a is a hycom grid file, default regional.grid.a.
C                 Note that the corresponding grid.b must also exist.
C
C                 mask.a contains an optional mask array, that turns
C                 off mapping anywhere it is 2.0**100 (data void).
C                 Note that grid.a must be specified if mask.a is used.
C
C  file.a will contain idm*jdm 32-bit IEEE real values for
C   the array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and values of 2.0**100 indicating a data void.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  March 2009.
C
      REAL*4, ALLOCATABLE :: A(:,:),W(:,:),
     +                       PLAT(:,:),PLON(:,:)
      REAL*4              :: PAD(4096)
      INTEGER IOS
      INTEGER      IARGC
      INTEGER      NARG
      CHARACTER*240 CARG
C
      INTEGER       IDM,JDM,NPAD
      REAL*4        DIST
      CHARACTER*6   CVARIN
      CHARACTER*240 CFILEP,CFILE,CFILEG,CFILEM,CFILEB
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.3) THEN
        CALL GETARG(1,CFILEP)
        CALL GETARG(2,CFILE)
        CALL GETARG(3,CARG)
        READ(CARG,*)    DIST
        CFILEG = 'regional.grid.a'
        CFILEM = 'NONE'
      ELSEIF (NARG.EQ.4) THEN
        CALL GETARG(1,CFILEP)
        CALL GETARG(2,CFILE)
        CALL GETARG(3,CARG)
        READ(CARG,*)    DIST
        CALL GETARG(4,CFILEG)
        CFILEM = 'NONE'
      ELSEIF (NARG.EQ.5) THEN
        CALL GETARG(1,CFILEP)
        CALL GETARG(2,CFILE)
        CALL GETARG(3,CARG)
        READ(CARG,*)    DIST
        CALL GETARG(4,CFILEG)
        CALL GETARG(5,CFILEM)
      ELSE
        WRITE(6,*) 
     +   'Usage:  hycom_scatter points.txt field.a dist ' //
     +   '[grid.a] [mask.a]'
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
        WRITE(6,*) 'hycom_scatter: bad header file ',
     &             CFILEB(1:LEN_TRIM(CFILEB))
        CALL EXIT(2)
      ENDIF
      READ( 11,*) JDM,CVARIN
      IF (CVARIN.NE.'jdm   ') THEN
        WRITE(6,*) 'hycom_scatter: bad header file ',
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
        WRITE(6,*) 'Error in hycom_scatter: could not allocate ',
     +             IDM*JDM,' words for A'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( W(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_scatter: could not allocate ',
     +             IDM*JDM,' words for W'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( PLAT(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_scatter: could not allocate ',
     +             IDM*JDM,' words for PLAT'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( PLON(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_scatter: could not allocate ',
     +             IDM*JDM,' words for PLON'
        CALL EXIT(2)
      ENDIF
C
      CALL SCATTER(A,W,PLAT,PLON,IDM,JDM,PAD,NPAD,
     +             DIST, CFILEP,CFILE,CFILEG,CFILEM)
      CALL EXIT(0)
 5000 FORMAT(I4)
      END
      SUBROUTINE SCATTER(A,W,PLAT,PLON,IDM,JDM, PAD,NPAD,
     +                   DIST, CFILEP,CFILE,CFILEG,CFILEM)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILEP,CFILE,CFILEG,CFILEM
      INTEGER      IDM,JDM,NPAD
      REAL*4       A(IDM,JDM),W(IDM,JDM),
     +             PLAT(IDM,JDM),PLON(IDM,JDM),PAD(NPAD),
     +             DIST
C
C     MOST OF WORK IS DONE HERE.
C
      CHARACTER*80 CLINE
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      INTEGER      I,J,K,IOS,NRECL
      REAL*4       AMX,AMN,D,DI,VALUE,VLAT,VLON
      REAL*4       SPHERDIST  !function
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
C     INPUT MASK ARRAY.
C
      INQUIRE( IOLENGTH=NRECL) W,PAD
C
      IF     (CFILEM.EQ.'NONE') THEN
        W(:,:) = 0.0
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
        READ(11,REC=1,IOSTAT=IOS) W
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(W,IDM*JDM)
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
C     INITIALIZE WEIGHTS (DISTANCE TO SAMPLED POINT)
C
      DO J= 1,JDM
        DO I= 1,IDM
          A(I,J) = SPVAL
          IF     (W(I,J).NE.SPVAL) THEN
            W(I,J) = DIST
          ENDIF
        ENDDO
      ENDDO
C
C     LOOP THROUGH ALL POINTS
C
      OPEN(UNIT=12,FILE=CFILEP,FORM='FORMATTED',
     &     STATUS='OLD',ACTION='READ')
      DO K= 1,9999999
        READ(12,'(a)',IOSTAT=IOS) CLINE
        IF     (IOS.NE.0) THEN
          EXIT
        ENDIF
        READ(CLINE,*) VALUE,VLON,VLAT
*       IF     (MOD(K,1000).EQ.1) THEN
*         WRITE(6,*) 'V = ',VALUE,VLON,VLAT
*       ENDIF
C
C       LOOP THROUGH ALL GRID POINTS
C
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (W(I,J).NE.SPVAL) THEN
              D = SPHERDIST(VLON,VLAT, PLON(I,J),PLAT(I,J))
              IF     (D.LT.W(I,J)) THEN  !nearest input point so far
                W(I,J) = D
                A(I,J) = VALUE
*               IF     (MOD(K,1000).EQ.1) THEN
*                 WRITE(6,*) 'D = ',D,I,J,A(I,J)
*               ENDIF
              ENDIF
            ENDIF
          ENDDO !i
        ENDDO !j
      ENDDO !points
      CLOSE(12)
C
C     FINAL RESULT
C
      AMN =  SPVAL
      AMX = -SPVAL
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (A(I,J).NE.SPVAL) THEN
            AMX = MAX( AMX, A(I,J) )
            AMN = MIN( AMN, A(I,J) )
          ENDIF
        ENDDO !i
      ENDDO !j
C
C     OUTPUT FILE.
C
      OPEN(UNIT=21, FILE=CFILE, FORM='UNFORMATTED', STATUS='NEW',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILE)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
      WRITE(21,REC=1,IOSTAT=IOS) A
      WRITE(6,'(a,1p2g16.8)')
     &   'min, max = ',AMN,AMX
      CLOSE(21)

      RETURN
      END
      real*4 function spherdist(lon1,lat1,lon2,lat2)
      implicit none
      real*4, intent(in) :: lon1,lat1,lon2,lat2 ! Pos. in degrees
c
c --- ------------------------------------------------
c --- Computes the distance between geo. pos.
c --- lon1,lat1 and lon2,lat2. 
c --- input is in degrees.
c
c --- output is real*4 for better global consistancy,
c --- by truncating double precision roundoff errors.
c --- real*4 is not in f90, but is widely supported.
c
c --- Based on m_spherdist.F90 from Geir Evanson.
c --- ------------------------------------------------
c
      double precision, parameter :: invradian=0.017453292d0
      double precision, parameter ::    rearth=6371001.0d0  ! Radius of earth
c
      double precision  dlon1,dlon2
      double precision  rlon1,rlat1,rlon2,rlat2           ! Pos. in radians
      double precision  x1,y1,z1,x2,y2,z2                 ! Cartesian position
      double precision  dr                                ! Arc length
c
c     ensure that spherdist(ax,ay,bx,by) == spherdist(bx,by,ax,ay)
c
      dlon1 = lon1
      dlon1 = mod(dlon1,360.d0)
      if     (dlon1.lt.0.d0) then
        dlon1 = dlon1 + 360.d0
      endif
      dlon2 = lon2
      dlon2 = mod(dlon2,360.d0)
      if     (dlon2.lt.0.d0) then
        dlon2 = dlon2 + 360.d0
      endif
      if     (lat1.lt.lat2) then
        rlon1=dlon1*invradian            !lon1 in rad
        rlat1=(90.d0-lat1)*invradian     !90-lat1 in rad 
        rlon2=dlon2*invradian            !lon2 in rad
        rlat2=(90.d0-lat2)*invradian     !90-lat2 in rad 
      elseif (lat1.eq.lat2 .and. dlon1.le.dlon2) then
        rlon1=dlon1*invradian            !lon1 in rad
        rlat1=(90.d0-lat1)*invradian     !90-lat1 in rad 
        rlon2=dlon2*invradian            !lon2 in rad
        rlat2=(90.d0-lat2)*invradian     !90-lat2 in rad 
      else
        rlon2=dlon1*invradian            !lon1 in rad
        rlat2=(90.d0-lat1)*invradian     !90-lat1 in rad 
        rlon1=dlon2*invradian            !lon2 in rad
        rlat1=(90.d0-lat2)*invradian     !90-lat2 in rad 
      endif
c
      x1= sin(rlat1)*cos(rlon1)        !x,y,z of pos 1.
      y1= sin(rlat1)*sin(rlon1)
      z1= cos(rlat1) 
c
      x2= sin(rlat2)*cos(rlon2)        !x,y,z of pos 2.
      y2= sin(rlat2)*sin(rlon2)
      z2= cos(rlat2) 
c
      dr=acos(min(1.d0,x1*x2+y1*y2+z1*z2))  ! Arc length
c
      spherdist=dr*rearth
c
      end function spherdist
