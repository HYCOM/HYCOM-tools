      SUBROUTINE SCRIP_READ
      USE mod_scrip  ! SCRIP remapping arrays and routines
      IMPLICIT NONE
C
C**********
C*
C  1) READ IN THE SCRIP WEIGHTS
C
C  2) THE PRECALCULATED SCRIP REGRIDDING WEIGHTS FILE IS FROM
C      ENVIRONMENT VARIABLE CDFSCRIP.
C     GENERATE THESE USING SCRIP, OR ESMF_RegridWeightGen, FROM
C      SCRIP SOURCE AND TARGET GRID FILES THAT CAN IN TURN BE
C      GENERATED FROM regional.grid USING hycom/ALL/bin/hycom_scrip_nc.
C     IN ORDER TO USE hycom_scrip_nc ON THE SCOURCE 'NATIVE' WIND
C      LAT-LON GRID THIS MUST BE MAPPED TO THE CORRESPONDING HYCOM
C      REGION.
C     THE SCRIP SOURCE AND TARGET GRIDS MUST BOTH BE UNMASKED AND ALL
C      TARGET GRID POINTS MUST GET AN INTERPOLATED VALUE.
C     FOR SCRIP, SEE http://climate.lanl.gov/Software/SCRIP/
C
C  3) ALAN J. WALLCRAFT, NRL, FEBRUARY 2012.
C*
C**********
C
      CHARACTER(CHAR_LEN) :: CSCRIP,MAP_NAME
C
      CSCRIP = ' '
      CALL GETENV('CDFSCRIP',CSCRIP)
      IF     (CSCRIP.EQ.' ') THEN
        WRITE(0,*) 'scrip_read: no CDFSCRIP environment variable'
        CALL EXIT(1)
        STOP
      ENDIF
      CALL READ_REMAP(MAP_NAME, CSCRIP)
C     END OF SCRIP_READ.
      END
      SUBROUTINE SCRIP(FLD,NX,NY, FLDI,NXI,NYI)
      USE mod_scrip  ! SCRIP remapping arrays and routines
      IMPLICIT NONE
C
      INTEGER NX,NY, NXI,NYI
      REAL*4  FLD(NX,NY)
      REAL*4  FLDI(NXI,NYI)
C
C**********
C*
C  1) INTERPOLATE FROM THE ARRAY FLDI TO THE ARRAY FLD.
C
C     PRECALCULATED WEIGHTS, IN mod_scrip, ARE USED FOR THE
C      INTERPOLATION.  INPUT THEM USING SCRIP_READ.
C
C  2) ARGUMENT LIST:
C       FLD     - INTERPOLATED ARRAY ON EXIT
C       NX,NY   - SIZE OF FLD ARRAY
C       FLDI    - ARRAY OF VALUES FROM WHICH TO INTERPOLATE
C       NXI,NYI - SIZE OF FLDI ARRAY
C
C  3) IF THE OUTPUT ARRAY IS ON A TRIPOLE GRID ITS TOP ROW
C      WILL BE ZERO ON EXIT.  CALL ARCUPD TO UPDATE IT.
C
C  4) ALAN J. WALLCRAFT, NRL, FEBRUARY 2012.
C*
C**********
C
C     LOCAL VARIABLES.
C
      REAL (KIND=DBL_KIND), DIMENSION(:), ALLOCATABLE ::
     &    G1_A,
     &    G2_A
      INTEGER I,II,J,NYA
c
c --- error check
c
      if     (grid1_size.ne.nxi*nyi) then
        write(6,*)
        write(6,*) 'scrip error: wrong grid1_size'
        write(6,*) 'nxi,nyi    = ',nxi,nyi
        write(6,*) 'nxi*nyi    = ',nxi*nyi
        write(6,*) 'grid1_size = ',grid1_size
        write(6,*)
        stop
      endif
c
      if     (grid2_size.eq.nx*ny) then
        nya = ny
      elseif (grid2_size.eq.nx*(ny-1)) then
        nya = ny-1
        fld(:,ny) = 0.0
      else
        write(6,*)
        write(6,*) 'scrip error: wrong grid2_size'
        write(6,*) 'nx,ny      = ',nx,ny
        write(6,*) 'nx*ny      = ',nx*ny
        write(6,*) 'grid2_size = ',grid2_size
        write(6,*)
        stop
      endif
C
      allocate (g1_a(grid1_size),
     &          g2_a(grid2_size) )
C
      do j= 1,nyi
        do i= 1,nxi
          ii = i + (j-1)*nxi
          g1_a(ii) = fldi(i,j)
        enddo !i
      enddo !j
c
      call remap(g2_a, wts_map1, grid2_add_map1, grid1_add_map1, g1_a)
c
      do j= 1,nya
        do i= 1,nx
          ii = i + (j-1)*nx
          if     (grid2_frac(ii).gt.0.d0) then
            fld(i,j) = g2_a(ii)/grid2_frac(ii)
          else
            write(6,*)
            write(6,*) 'scrip error: interpolation does not cover array'
            write(6,*) 'fld undefined at i,j = ',i,j
            write(6,*)
            stop
          endif
        enddo !i
      enddo !j
      if     (nya.ne.ny) then
        do i= 1,nx
          fld(i,ny) = 0.0  !defined later by arcupd
        enddo !i
      endif
c
      deallocate (g1_a, g2_a)
      return
C     END OF SCRIP.
      end
