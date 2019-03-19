c
c --- main program below module
c
      MODULE EULER_ROTATION
      IMPLICIT NONE

!-----------------------------------------------------------------------
! euler rotation of model geometry
!-----------------------------------------------------------------------
! based on code fragments provided for AOMIP by David Holland.
! http://fish.cims.nyu.edu/project_aomip/model_grid/overview.html
!-----------------------------------------------------------------------
!
! C_ALPHA : rotation about z-axis
! C_BETA  : rotation about new y-axis 
! C_GAMMA : rotation about new z-axis 
!
! C_FORWARD : rotation matricies to geographic coordinates
! C_REVERSE : rotation matricies to model coordinates
!
!-----------------------------------------------------------------------

      REAL*8 :: 
     &      C_ALPHA  = 000.0   

      REAL*8 ::
     &      C_BETA   = 000.0   

      REAL*8 ::
     &      C_GAMMA  = 000.0   

      REAL*8,
     &      DIMENSION (1:3, 1:3) ::
     &      C_FORWARD = 0.0

      REAL*8,
     &      DIMENSION (1:3, 1:3) ::
     &      C_REVERSE = 0.0

      REAL*8, PRIVATE, PARAMETER :: C_EPSILON=1.0d-6
      REAL*8, PRIVATE, PARAMETER ::    RADIAN=57.29577951D0
      REAL*8, PRIVATE, PARAMETER :: INVRADIAN=1.d0/RADIAN

!-----------------------------------------------------------------------
      CONTAINS
!-----------------------------------------------------------------------

      REAL*8 FUNCTION C_SIND(X)
      REAL*8 X
      C_SIND = SIN(X*INVRADIAN)
      END FUNCTION C_SIND

      REAL*8 FUNCTION C_COSD(X)
      REAL*8 X
      C_COSD = COS(X*INVRADIAN)
      END FUNCTION C_COSD

      REAL*8 FUNCTION C_INV_SIND(X)
      REAL*8 X
      C_INV_SIND = RADIAN*ASIN(X)
      END FUNCTION C_INV_SIND

      REAL*8 FUNCTION C_INV_COSD(X)
      REAL*8 X
      C_INV_COSD = RADIAN*ACOS(X)
      END FUNCTION C_INV_COSD

      SUBROUTINE C_EULER_ROTATION_MATRIX 
 
!-----------------------------------------------------------------------
! purpose: compute euler rotation matrix 
!-----------------------------------------------------------------------
!
! C_ALPHA   : rotation about z-axis  
! C_BETA    : rotation about new y-axis  
! C_GAMMA   : rotation about new z-axis      
! C_REVERSE : rotation to model coordinates 
! C_FORWARD : rotation to geographic coordinates
!
!-----------------------------------------------------------------------

      !local vars
      REAL*8 :: 
     &
     &      COSA, 
     &      SINA, 
     &      COSB, 
     &      SINB, 
     &      COSG, 
     &      SING

!-----------------------------------------------------------------------

      !compute sin and cos of all angles in degrees
      COSA  = C_COSD (C_ALPHA )
      SINA  = C_SIND (C_ALPHA )
      COSB  = C_COSD (C_BETA  )
      SINB  = C_SIND (C_BETA  )
      COSG  = C_COSD (C_GAMMA )
      SING  = C_SIND (C_GAMMA )
 
      !compute rotation matrix into model coordinates
      C_REVERSE (1,1) =   COSA * COSB * COSG  -  SINA * SING
      C_REVERSE (1,2) =   SINA * COSB * COSG  +  COSA * SING
      C_REVERSE (1,3) =        - SINB * COSG
      C_REVERSE (2,1) = - COSA * COSB * SING  -  SINA * COSG
      C_REVERSE (2,2) = - SINA * COSB * SING  +  COSA * COSG
      C_REVERSE (2,3) =          SINB * SING
      C_REVERSE (3,1) =   COSA * SINB
      C_REVERSE (3,2) =   SINA * SINB
      C_REVERSE (3,3) =          COSB
    
      !compute rotation matrix into geographical coordinates
      C_FORWARD (1,1) =   COSG * COSB * COSA  -  SING * SINA
      C_FORWARD (1,2) = - SING * COSB * COSA  -  COSG * SINA
      C_FORWARD (1,3) =          SINB * COSA
      C_FORWARD (2,1) =   COSG * COSB * SINA  +  SING * COSA
      C_FORWARD (2,2) = - SING * COSB * SINA  +  COSG * COSA
      C_FORWARD (2,3) =          SINB * SINA
      C_FORWARD (3,1) = - COSG * SINB
      C_FORWARD (3,2) =   SING * SINB
      C_FORWARD (3,3) =          COSB
   
!-----------------------------------------------------------------------
      END SUBROUTINE C_EULER_ROTATION_MATRIX


      SUBROUTINE C_MODEL_TO_GEO (X_IN_ANGLE, 
     &                           Y_IN_ANGLE, 
     &                           X_OUT_ANGLE, 
     &                           Y_OUT_ANGLE)

!---------------------------------------------------------------------
! purpose: a coordinate transformation from model coordinates
!          to geographic coordinates
!--------------------------------------------------------------------
!
! variables
!
! X_IN_ANGLE  : X-COORDINATE IN MODEL SYSTEM.
! Y_IN_ANGLE  : Y-COORDINATE IN MODEL SYSTEM.
! X_OUT_ANGLE : X-COORDINATE IN GEOGRAPHICAL SYSTEM.
! Y_OUT_ANGLE : Y-COORDINATE IN GEOGRAPHICAL SYSTEM.
!
! THE ROTATION FROM THE COORDINATE SYSTEM,
! X''',Y''',Z''' BACK TO THE ORIGINAL COORDINATE SYTEM
! CAN BE PERFORMED BY DOING THE ROTATION IN
! REVERSED ORDER AND WITH OPPOSITE SIGN ON THE ROTATED ANGLES.
!
!---------------------------------------------------------------------

      !passed arguments
      REAL*8,
     &      INTENT (IN) :: 
     &
     &      X_IN_ANGLE,
     &      Y_IN_ANGLE

      !passed arguments
      REAL*8,
     &      INTENT (OUT) :: 
     &
     &      X_OUT_ANGLE,
     &      Y_OUT_ANGLE

      !local vars
      REAL*8 :: 
     &
     &      X_TMP_ANGLE,
     &      Y_TMP_ANGLE,
     &      XX,
     &      YY,
     &      ZZ, 
     &      XNEW,
     &      YNEW,
     &      ZNEW,
     &      COSTN,
     &      PHI_NEW,
     &      THETA_NEW

!---------------------------------------------------------------------

      !initial angles of rotation
      X_TMP_ANGLE = X_IN_ANGLE
      Y_TMP_ANGLE = Y_IN_ANGLE

      !avoid trouble of an exactly zero angle by adding offset
      X_TMP_ANGLE = X_TMP_ANGLE + C_EPSILON
      Y_TMP_ANGLE = Y_TMP_ANGLE + C_EPSILON

      !spherical coordiantes to cartesian coordinates
      XX = C_COSD (Y_TMP_ANGLE) * C_COSD (X_TMP_ANGLE)
      YY = C_COSD (Y_TMP_ANGLE) * C_SIND (X_TMP_ANGLE)
      ZZ = C_SIND (Y_TMP_ANGLE)
 
      !new cartesian coordinates are given by
      XNEW = C_FORWARD (1,1) * XX 
     &     + C_FORWARD (1,2) * YY 
     &     + C_FORWARD (1,3) * ZZ
      YNEW = C_FORWARD (2,1) * XX 
     &     + C_FORWARD (2,2) * YY 
     &     + C_FORWARD (2,3) * ZZ
      ZNEW = C_FORWARD (3,1) * XX 
     &     + C_FORWARD (3,2) * YY 
     &     + C_FORWARD (3,3) * ZZ
 
      !obtain new angles THETA_NEW,COSTN,PHI_NEW
      THETA_NEW  = C_INV_SIND (ZNEW)
      COSTN = SQRT (1.0 - ZNEW**2)

      IF    ((XNEW > 0.0) .AND. (YNEW > 0.0)) THEN

         IF (XNEW < YNEW) THEN
            PHI_NEW = C_INV_COSD (XNEW/COSTN)
         ELSE
            PHI_NEW = C_INV_SIND (YNEW/COSTN)
         ENDIF

      ELSEIF ((XNEW < 0.0) .AND. (YNEW > 0.0)) THEN

         IF (ABS (XNEW) < YNEW) THEN
            PHI_NEW = 180.0 - C_INV_COSD (ABS (XNEW)/COSTN)
         ELSE
            PHI_NEW = 180.0 - C_INV_SIND (YNEW/COSTN)
         ENDIF

      ELSEIF ((XNEW < 0.0) .AND. (YNEW < 0.0)) THEN

         IF (ABS (XNEW) < ABS (YNEW)) THEN
            PHI_NEW =-180.0 + C_INV_COSD (ABS (XNEW)/COSTN)
         ELSE
            PHI_NEW =-180.0 + C_INV_SIND (ABS (YNEW)/COSTN)
         ENDIF

      ELSEIF ((XNEW > 0.0) .AND. (YNEW < 0.0)) THEN

        IF(    XNEW  < ABS (YNEW)) THEN
           PHI_NEW =    - C_INV_COSD (ABS (XNEW)/COSTN)
        ELSE
           PHI_NEW =    - C_INV_SIND (ABS (YNEW)/COSTN)
        ENDIF

      ENDIF

      !new spherical coordinates 
      X_OUT_ANGLE = PHI_NEW
      Y_OUT_ANGLE = THETA_NEW

*     IF (X_OUT_ANGLE < 0.0) X_OUT_ANGLE = X_OUT_ANGLE + 360.0

      !avoid trouble of an exactly zero angle by subtracting offset
      X_OUT_ANGLE = X_OUT_ANGLE - C_EPSILON
      Y_OUT_ANGLE = Y_OUT_ANGLE - C_EPSILON

!-----------------------------------------------------------------------
      END SUBROUTINE C_MODEL_TO_GEO


!-----------------------------------------------------------------------
      END MODULE EULER_ROTATION


      program grid_rotated
      use mod_za  ! HYCOM array I/O interface
      use euler_rotation
      implicit none
c
      logical          lperiod
      integer          i,j,mapflg
      double precision ralpha,rbeta,rgamma
      real*8           hmaxa,hmina
      real*8           xi,xo,yi,yo
c
c --- from an existing grid definition file, produce a new grid
c --- definition file after rotating the coordinates.  The existing
c --- file is taken to identify the grids in the rotated frame, and
c --- the new file will identify these same grids in true lon-lat.
c
c --- the new grid will have the same grid sizes in meters, but 
c --- a) the true lon and lat will replace the pre-rotation lon and lat,
c --- b) the coriolis   will be w.r.t. true lon and lat, and
c --- c) the grid angle will be w.r.t. true W-E,S-N
c
c --- the model, and pre-processing programs that involve interpolation,
c --- must use the new grid definition file, but post-processing 
c --- diagnostic programs that don't use coriollis or the grid angle
c --- can either use the true lon-lats (new grid definition file) or
c --- stay in "rotated space" by using the original grid definition
c --- file which has lon-lats w.r.t. the rotated coordinates.
c
c --- based on the AOMIP Euler rotation code fragments at
c --- http://fish.cims.nyu.edu/project_aomip/model_grid/overview.html
c
      double precision, parameter :: half   = 0.5d0
      double precision, parameter :: halfpi = 1.5707963268d0
      double precision, parameter :: radian = 57.29578d0
c
      real*8 spherdist  ! fn for distance between geo. pos
c
      integer, allocatable :: ip(:,:)
      real*8,  allocatable :: pang(:,:),cori(:,:)
      real*8,  allocatable :: plon(:,:),qlon(:,:),ulon(:,:),vlon(:,:)
      real*8,  allocatable :: plat(:,:),qlat(:,:),ulat(:,:),vlat(:,:)
      real*8,  allocatable :: pscx(:,:),qscx(:,:),uscx(:,:),vscx(:,:)
      real*8,  allocatable :: pscy(:,:),qscy(:,:),uscy(:,:),vscy(:,:)
c
      call xcspmd  !input idm,jdm
      allocate(   ip(idm,jdm) )
      allocate( pang(idm,jdm), cori(idm,jdm) )
      allocate( plat(idm,jdm), plon(idm,jdm) )
      allocate( qlat(idm,jdm), qlon(idm,jdm) )
      allocate( ulat(idm,jdm), ulon(idm,jdm) )
      allocate( vlat(idm,jdm), vlon(idm,jdm) )
      allocate( pscy(idm,jdm), pscx(idm,jdm) )
      allocate( qscy(idm,jdm), qscx(idm,jdm) )
      allocate( uscy(idm,jdm), uscx(idm,jdm) )
      allocate( vscy(idm,jdm), vscx(idm,jdm) )
c
c --- read in the map projection.
c ---   'idm   ' = longitudinal array size
c ---   'jdm   ' = latitudinal  array size
c ---   'mapflg' = map flag (0=mercator,2=uniform,4=f-plane)
c ---   'ralpha' = 1st euler rotation angle
c ---   'rbeta ' = 2nd euler rotation angle
c ---   'rgamma' = 3rd euler rotation angle
c 
      call blkini(i,      'idm   ')
      call blkini(j,      'jdm   ')
      call blkini(mapflg, 'mapflg')
      call blkinr(ralpha,
     &           'ralpha','("blkinr: ",a6," =",f11.4," deg")')
      call blkinr(rbeta, 
     &           'rbeta ','("blkinr: ",a6," =",f11.4," deg")')
      call blkinr(rgamma,
     &           'rgamma','("blkinr: ",a6," =",f11.4," deg")')
c
      if     (i.ne.idm .or. j.ne.jdm) then
        write(lp,'(/a,a/)') 'stdin and regional.grid.b have',
     &                      ' different idm,jdm values'
        call flush(lp)
        stop
      endif
c
c --- initialize rotated grid
c
      c_alpha = ralpha
      c_beta  = rbeta 
      c_gamma = rgamma
c
      call c_euler_rotation_matrix
c
      write(6,*) 'c_angles  = ',c_alpha,c_beta,c_gamma
      write(6,*)
      write(6,*) 'c_forward = ',C_FORWARD
      write(6,*)
      write(6,*) 'c_reverse = ',C_REVERSE
      write(6,*)
c
c --- read in the original, un-rotated, grid.
c
      call zaiost
c
      call zaiopn( 'old', 51)
      call zaiord8(plon, ip,.false., hmina,hmaxa, 51)
      write(6, 6100) 'r-plon',hmina,hmaxa
      call zaiord8(plat, ip,.false., hmina,hmaxa, 51)
      write(6, 6100) 'r-plat',hmina,hmaxa
      call zaiord8(qlon, ip,.false., hmina,hmaxa, 51)
      write(6, 6100) 'r-qlon',hmina,hmaxa
      call zaiord8(qlat, ip,.false., hmina,hmaxa, 51)
      write(6, 6100) 'r-qlat',hmina,hmaxa
      call zaiord8(ulon, ip,.false., hmina,hmaxa, 51)
      write(6, 6100) 'r-ulon',hmina,hmaxa
      call zaiord8(ulat, ip,.false., hmina,hmaxa, 51)
      write(6, 6100) 'r-ulat',hmina,hmaxa
      call zaiord8(vlon, ip,.false., hmina,hmaxa, 51)
      write(6, 6100) 'r-vlon',hmina,hmaxa
      call zaiord8(vlat, ip,.false., hmina,hmaxa, 51)
      write(6, 6100) 'r-vlat',hmina,hmaxa
      call zaiocl( 51)
c
      lperiod = maxval(plon(:,:))-minval(plon(:,:)) .gt. 350.0
c
c --- rotate each grid.
c
      do j= 1,jdm
        do i= 1,idm
          xi = plon(i,j)
          yi = plat(i,j)
          call c_model_to_geo(xi,yi, xo,yo)
          plon(i,j) = xo
          plat(i,j) = yo
          if     (i.eq.idm/2 .and. j.eq.jdm/2) then
            write(6,*) 'i,j,p.rot = ',i,j,xi,yi
            write(6,*) 'i,j,p.geo = ',i,j,xo,yo
          endif
c
          xi = qlon(i,j)
          yi = qlat(i,j)
          call c_model_to_geo(xi,yi, xo,yo)
          qlon(i,j) = xo
          qlat(i,j) = yo
          if     (i.eq.idm/2 .and. j.eq.jdm/2) then
            write(6,*) 'i,j,q.rot = ',i,j,xi,yi
            write(6,*) 'i,j,q.geo = ',i,j,xo,yo
          endif
c
          xi = ulon(i,j)
          yi = ulat(i,j)
          call c_model_to_geo(xi,yi, xo,yo)
          ulon(i,j) = xo
          ulat(i,j) = yo
          if     (i.eq.idm/2 .and. j.eq.jdm/2) then
            write(6,*) 'i,j,u.rot = ',i,j,xi,yi
            write(6,*) 'i,j,u.geo = ',i,j,xo,yo
          endif
c
          xi = vlon(i,j)
          yi = vlat(i,j)
          call c_model_to_geo(xi,yi, xo,yo)
          vlon(i,j) = xo
          vlat(i,j) = yo
          if     (i.eq.idm/2 .and. j.eq.jdm/2) then
            write(6,*) 'i,j,v.rot = ',i,j,xi,yi
            write(6,*) 'i,j,v.geo = ',i,j,xo,yo
          endif
        enddo
      enddo
c
      write(6, *)
      do j= 1,jdm
        write(6,'(a,i5,2f10.3)')
     &    'j,qlat =',j,minval(qlat(:,j)),maxval(qlat(:,j))
        write(6,'(a,i5,2f10.3)')
     &    'j,plat =',j,minval(plat(:,j)),maxval(plat(:,j))
      enddo
      write(6, *)
      do i= 1,idm
        write(6,'(a,i5,2f10.3)')
     &    'i,qlon =',i,minval(qlon(i,:)),maxval(qlon(i,:))
        write(6,'(a,i5,2f10.3)')
     &    'i,plon =',i,minval(plon(i,:)),maxval(plon(i,:))
      enddo
      call zhflsh(6)
c
c --- finite difference approximation for pang.
c
      call rotang(plat,plon,idm,jdm,pang)
c
c --- write header.
c
      call zhopen(61, 'formatted', 'new', 0)
      write(61,'(i5,a)')
     &  idm,   "    'idm   ' = longitudinal array size"
      write(61,'(i5,a)')
     &  jdm,   "    'jdm   ' = latitudinal  array size"
      write(61,'(i5,a,a)')
     &  mapflg,"    'mapflg' = map flag",
     &         " (-1=unknown,0=mercator,2=uniform,4=f-plane)"
      call zhflsh(61)
c
      write(6, *)
      write(6,'(i5,a)')
     &  idm,   "    'idm   ' = longitudinal array size"
      write(6,'(i5,a)')
     &  jdm,   "    'jdm   ' = latitudinal  array size"
      write(6,'(i5,a)')
     &  mapflg,"    'mapflg' = map flag (-1=unknown,0=mercator,...)"
      call zhflsh(6)
c
c --- write grid arrays.
c
      call zaiopn('new', 61)
      call zaiowr8(plon, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6100) 'plon',hmina,hmaxa
      write(6, 6100) 'plon',hmina,hmaxa
      call zaiowr8(plat, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6100) 'plat',hmina,hmaxa
      write(6, 6100) 'plat',hmina,hmaxa
      call zaiowr8(qlon, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6100) 'qlon',hmina,hmaxa
      write(6, 6100) 'qlon',hmina,hmaxa
      call zaiowr8(qlat, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6100) 'qlat',hmina,hmaxa
      write(6, 6100) 'qlat',hmina,hmaxa
      call zaiowr8(ulon, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6100) 'ulon',hmina,hmaxa
      write(6, 6100) 'ulon',hmina,hmaxa
      call zaiowr8(ulat, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6100) 'ulat',hmina,hmaxa
      write(6, 6100) 'ulat',hmina,hmaxa
      call zaiowr8(vlon, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6100) 'vlon',hmina,hmaxa
      write(6, 6100) 'vlon',hmina,hmaxa
      call zaiowr8(vlat, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6100) 'vlat',hmina,hmaxa
      write(6, 6100) 'vlat',hmina,hmaxa
      call zaiowr8(pang, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6100) 'pang',hmina,hmaxa
      write(6, 6100) 'pang',hmina,hmaxa
      write(6, *)
 6100 format(a,':  min,max = ',2f15.5)
      call zhflsh(61)
      call zhflsh(6)
c
c --- grid spacing in meters.
c
      if     (lperiod) then
        write(lp,'(a)') 'domain assumed to be periodic'
      else
        write(lp,'(a)') 'domain assumed to be non-periodic'
      endif
c
      do j= 1,jdm
        do i= 2,idm
          uscx(i,j) = spherdist(plon(i,  j),plat(i,  j),
     &                          plon(i-1,j),plat(i-1,j) )
          qscx(i,j) = spherdist(vlon(i,  j),vlat(i,  j),
     &                          vlon(i-1,j),vlat(i-1,j) )
        enddo
        i=1
        if     (lperiod) then
          uscx(i,j) = spherdist(plon(i,  j),plat(i,  j),
     &                          plon(idm,j),plat(idm,j) )
          qscx(i,j) = spherdist(vlon(i,  j),vlat(i,  j),
     &                          vlon(idm,j),vlat(idm,j) )
        else
          uscx(i,j) = uscx(i+1,j)  ! updated below except in corner
          qscx(i,j) = qscx(i+1,j)  ! updated below except in corner
        endif
        do i= 1,idm-1
          vscx(i,j) = spherdist(qlon(i,  j),qlat(i,  j),
     &                          qlon(i+1,j),qlat(i+1,j) )
          pscx(i,j) = spherdist(ulon(i,  j),ulat(i,  j),
     &                          ulon(i+1,j),ulat(i+1,j) )
          if     (pscx(i,j).gt.100000.0) then
            write(6,*) 'spherdist ',ulon(i,  j),ulat(i,  j),
     &                              ulon(i+1,j),ulat(i+1,j),
     &                              pscx(i,j)
          endif
        enddo
        i=idm
        if     (lperiod) then
          vscx(i,j) = spherdist(qlon(i,  j),qlat(i,  j),
     &                          qlon(1,  j),qlat(1,  j) )
          pscx(i,j) = spherdist(ulon(i,  j),ulat(i,  j),
     &                          ulon(1,  j),ulat(1,  j) )
        else
          vscx(i,j) = vscx(i-1,j)  ! updated below except in corner
          pscx(i,j) = pscx(i-1,j)  ! updated below except in corner
        endif
      enddo
c
      do j= 1,jdm
        if     (j.ne.1) then
          do i= 1,idm
            vscy(i,j) = spherdist(plon(i,  j),plat(i,  j),
     &                            plon(i,j-1),plat(i,j-1) )
            qscy(i,j) = spherdist(ulon(i,  j),ulat(i,  j),
     &                            ulon(i,j-1),ulat(i,j-1) )
          enddo
        endif
        if     (j.ne.jdm) then
          do i= 1,idm
            uscy(i,j) = spherdist(qlon(i,  j),qlat(i,  j),
     &                            qlon(i,j+1),qlat(i,j+1) )
            pscy(i,j) = spherdist(vlon(i,  j),vlat(i,  j),
     &                            vlon(i,j+1),vlat(i,j+1) )
          enddo
        endif
      enddo
c
c     fill in the edges, assuming constant grid aspect ratio
c
      j=1
      do i= 1,idm
        vscy(i,j) = vscx(i,j)*vscy(i,j+1)/vscx(i,j+1)
        qscy(i,j) = qscx(i,j)*qscy(i,j+1)/qscx(i,j+1)
      enddo
      j=jdm
      do i= 1,idm
        uscy(i,j) = uscx(i,j)*uscy(i,j-1)/uscx(i,j-1)
        pscy(i,j) = pscx(i,j)*pscy(i,j-1)/pscx(i,j-1)
      enddo
c
      if     (.not.lperiod) then
        i=1
        do j= 1,jdm
          uscx(i,j) = uscy(i,j)*uscx(i+1,j)/uscy(i+1,j)
          qscx(i,j) = qscy(i,j)*qscx(i+1,j)/qscy(i+1,j)
        enddo
        i=idm
        do j= 1,jdm
          vscx(i,j) = vscy(i,j)*vscx(i-1,j)/vscy(i-1,j)
          pscx(i,j) = pscy(i,j)*pscx(i-1,j)/pscy(i-1,j)
        enddo
      endif
c
      write(lp,*)
      do j= 1,jdm
        write(lp,'(a,i5,3f10.2)')
     &    'j,vy =',j,minval(vscy(:,j)),maxval(vscy(:,j)),
     &               maxval(vscy(:,j))-minval(vscy(:,j))
        write(lp,'(a,i5,3f10.2)')
     &    'j,vx =',j,minval(vscx(:,j)),maxval(vscx(:,j)),
     &               maxval(vscx(:,j))-minval(vscx(:,j))
        write(lp,'(a,i5,3f10.2)')
     &    'j,uy =',j,minval(uscy(:,j)),maxval(uscy(:,j)),
     &               maxval(uscy(:,j))-minval(uscy(:,j))
        write(lp,'(a,i5,3f10.2)')
     &    'j,ux =',j,minval(uscx(:,j)),maxval(uscx(:,j)),
     &               maxval(uscx(:,j))-minval(uscx(:,j))
      enddo
      write(lp,*)
c
c --- write grid arrays.
c
      call zaiowr8(pscx, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6110) 'pscx',hmina,hmaxa
      write(lp,6110) 'pscx',hmina,hmaxa
      call zaiowr8(pscy, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6110) 'pscy',hmina,hmaxa
      write(lp,6110) 'pscy',hmina,hmaxa
      call zaiowr8(qscx, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6110) 'qscx',hmina,hmaxa
      write(lp,6110) 'qscx',hmina,hmaxa
      call zaiowr8(qscy, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6110) 'qscy',hmina,hmaxa
      write(lp,6110) 'qscy',hmina,hmaxa
      call zaiowr8(uscx, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6110) 'uscx',hmina,hmaxa
      write(lp,6110) 'uscx',hmina,hmaxa
      call zaiowr8(uscy, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6110) 'uscy',hmina,hmaxa
      write(lp,6110) 'uscy',hmina,hmaxa
      call zaiowr8(vscx, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6110) 'vscx',hmina,hmaxa
      write(lp,6110) 'vscx',hmina,hmaxa
      call zaiowr8(vscy, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6110) 'vscy',hmina,hmaxa
      write(lp,6110) 'vscy',hmina,hmaxa
      write(lp,*)
 6110 format(a,':  min,max = ',2f15.5)
c
c --- coriolis
c
      do j= 1,jdm
        do i= 1,idm
          cori(i,j)=sin(qlat(i,j)/radian)*
     &              8.d0*halfpi/86164.0d0  ! sidereal day
        enddo
      enddo
c
      call zaiowr8(cori, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6120) 'cori',hmina,hmaxa
      write(lp,6120) 'cori',hmina,hmaxa
      write(lp,*)
 6120 format(a,':  min,max = ',2f15.10)
c
c --- grid aspect ratios
c
      do j= 1,jdm
        do i= 1,idm
          if     (pscy(i,j).eq.0.0) then
            pscx(i,j) = 99.0
          elseif (pscx(i,j).ge.99.0*pscy(i,j)) then
            pscx(i,j) = 99.0
          else
            pscx(i,j) = pscx(i,j)/pscy(i,j)
          endif
        enddo
      enddo
c
      call zaiowr8(pscx, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6130) 'pasp',hmina,hmaxa
      write(lp,6130) 'pasp',hmina,hmaxa
      write(lp,*)
 6130 format(a,':  min,max = ',2f15.5)
c
      close(unit=61)
      call zaiocl(61)
c
      end
      subroutine blkinr(rvar,cvar,cfmt)
      implicit none
c
      double precision rvar
      character        cvar*6,cfmt*(*)
c
      integer       lp
      common/linepr/lp
c
c     read in one real value from stdin
c
      character*6 cvarin
c
      read(*,*) rvar,cvarin
      write(lp,cfmt) cvarin,rvar
      call flush(lp)
c
      if     (cvar.ne.cvarin) then
        write(lp,*) 
        write(lp,*) 'error in blkinr - input ',cvarin,
     +                      ' but should be ',cvar
        write(lp,*) 
        call flush(lp)
        stop
      endif
      return
      end
      subroutine blkini(ivar,cvar)
      implicit none
c
      integer     ivar
      character*6 cvar
c
      integer       lp
      common/linepr/lp
c
c     read in one integer value from stdin
c
      character*6 cvarin
c
      read(*,*) ivar,cvarin
      write(lp,6000) cvarin,ivar
      call flush(lp)
c
      if     (cvar.ne.cvarin) then
        write(lp,*) 
        write(lp,*) 'error in blkini - input ',cvarin,
     +                      ' but should be ',cvar
        write(lp,*) 
        call flush(lp)
        stop
      endif
      return
 6000 format('blkini: ',a6,' =',i6)
      end
