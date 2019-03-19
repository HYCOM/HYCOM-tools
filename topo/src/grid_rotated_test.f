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

      REAL*8, PRIVATE, PARAMETER :: C_EPSILON=1.0d-5
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
      write(6,*) 'X_TMP_ANGLE,Y_TMP_ANGLE = ',X_TMP_ANGLE,Y_TMP_ANGLE
      call flush(6)

      !spherical coordiantes to cartesian coordinates
      XX = C_COSD (Y_TMP_ANGLE) * C_COSD (X_TMP_ANGLE)
      YY = C_COSD (Y_TMP_ANGLE) * C_SIND (X_TMP_ANGLE)
      ZZ = C_SIND (Y_TMP_ANGLE)
      write(6,*) 'XX,YY,ZZ = ',XX,YY,ZZ
      call flush(6)
 
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
      write(6,*) 'XNEW,YNEW,ZNEW = ',XNEW,YNEW,ZNEW
      call flush(6)
 
      !obtain new angles THETA_NEW,COSTN,PHI_NEW
      THETA_NEW  = C_INV_SIND (ZNEW)
      write(6,*) 'THETA_NEW = ',THETA_NEW
      call flush(6)
      COSTN = SQRT (1.0 - ZNEW**2)
      write(6,*) 'COSTN     = ',COSTN
      call flush(6)

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
      write(6,*) 'PHI_NEW   = ',PHI_NEW
      call flush(6)

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


      program test
      use euler_rotation
      implicit none
c
      integer i,j
      real*8  xi,xo,yi,yo
c
c --- based on the AOMIP Euler rotation code fragments at
c --- http://fish.cims.nyu.edu/project_aomip/model_grid/overview.html
c
c --- initialize rotated grid
c
      c_alpha = -30.0
      c_beta  = -90.0
      c_gamma =   0.0
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
      do j= -30,30,10
        do i= -30,30,10
          xi = i
          yi = j
          call c_model_to_geo(xi,yi, xo,yo)
          write(6,*) 'xi,yi = ',xi,yi
          write(6,*) 'xo,yo = ',xo,yo
        enddo
      enddo
      end
