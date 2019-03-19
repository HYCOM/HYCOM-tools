      program vgrd2ncom
      implicit none
c
c --- convert a HYCOM sigma-Z grid to NCOM
c
      integer          k,kdm,nsigma
      real             al,als,dp00,dp00f,z(0:999)
      double precision zz,dz
c
      integer       lp
      common/linepr/lp
c
      lp = 6
c
c --- 'kdm   ' = number of layers
c --- 'nsigma' = number of sigma  levels (kdm-nsigma z-levels)
c
      call blkini(kdm,   'kdm   ')
      call blkini(nsigma,'nsigma')
c
c --- 'dp00  ' = z-level spacing minimum thickness (m)
c --- 'dp00f ' = z-level spacing stretching factor 
c
      call blkinr(dp00,  'dp00  ','(a6," =",f10.4," m")')
      call blkinr(dp00f, 'dp00f ','(a6," =",f10.4,"  ")')
c
      z(0) = 0.0
      zz   = 0.0d0
      dz   = dp00
      do k= 1,kdm
        zz = zz + dz
        z(k) = -zz
        dz = dz * dp00f
        write(6,'(a,i3,f10.3,f10.5)') 
     &    'k,zz,dz = ',k,zz,dz
      enddo
      open(unit=21,file='ovgrd_1.D',form='unformatted')
      al  = kdm+1
      als = nsigma+1
      write(21) al,als,z(0:kdm)
      end
