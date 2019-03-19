      program mom6nc82field
      use mod_mom6  ! HYCOM mom6 array interface
      use mod_za    ! HYCOM array I/O interface
c
      implicit none
c
c --- convert MOM6 p-grid 2-D output to a HYCOM .[ab] file.
c --- requires the HYCOM regional.grid for the MOM6 domain.
c --- this version reads MOM6 real*8 and writes HYCOM real*4
c
      character*256    flnm_o,flnm_t,name_t
      logical          larctic
      integer          i,ia,irec,j,jja,mro
      integer          inirec,increc,maxrec
      integer          yrflag
      real             xmin,xmax
      double precision time3(3)
c
c --- spval  = hycom data void marker, 2^100 or about 1.2676506e30
      real, parameter :: spval=2.0**100
c
      call xcspmd  !input the HYCOM array dimensions
      call zaiost  !initialize HYCOM I/O
      lp     = 6
      ii     = idm
      jj     = jdm
      yrflag = 3
c
c --- 'name_t' = name of mom6  scalar  field
c --- 'flnm_t' = name of mom6  scalars file  (input)
c ---             ends in ".nc" for a single netCDF file, or
c ---                     ".nc.DDDD-EEEE" for subregion netCDF files
c --- 'flnm_o' = name of hycom scalars file (output)
c --- 'inirec' = initial time record to read in
c --- 'maxrec' = final   time record to read in, 0 for to the end
c --- 'increc' =         time record increment
c
      read (*,'(a)') name_t
      write (lp,'(2a)') ' input MOM6 field: ',trim(name_t)
      call flush(lp)
      read (*,'(a)') flnm_t
      i = len_trim(flnm_t)
      if     (flnm_t(i-2:i).eq.'.nc') then
        write (lp,'(2a)') ' input MOM6  file: ',trim(flnm_t)
      else
        write (lp,'(2a)') ' input MOM6 files: ',flnm_t(1:i-5)
        write (lp,'(3a)') '               to: ',flnm_t(1:i-9),
     &                                          flnm_t(i-3:i)
      endif
      call flush(lp)
      read (*,'(a)') flnm_o
      write (lp,'(2a)') 'output HYCOM file: ',trim(flnm_o)
      call flush(lp)
c
      call blkini(inirec,  'inirec')
      call blkini(maxrec,  'maxrec')
      call blkini(increc,  'increc')
c
c --- mom6 dimensions
c
      call rd_dimen2(nto,mto,mro, flnm_t,name_t)
c
      larctic = mto .eq. jdm-1
      if     (larctic) then
        jja = jj-1
      else
        jja = jj
      endif
c
      write(lp,*) 
      write(lp,*) 'nto,mto = ',nto,mto
      write(lp,*) 'ii,jj   = ',ii, jj
      write(lp,*) 'larctic = ',larctic
      write(lp,*) 
      call zhflsh(lp)
c
      if     (maxrec.eq.0) then
        maxrec = mro
      elseif (maxrec.gt.mro) then
        write(lp,*) 
        write(lp,*) 'maxrec .gt. time dimension = ',mro
        write(lp,*) 
        call zhflsh(lp)
        stop
      endif
c
c --- array allocation
c
      call mom6_alloc_field8
c
c     open the output .[ab] files.
c
      i = len_trim(flnm_o)
      call zaiopf(      flnm_o(1:i-2)//'.a',       'new',21)
      open(unit=21,file=flnm_o(1:i-2)//'.b',status='new',
     &     form='formatted',action='write')
c
c --- read the mom6 file.
c
      do irec= inirec,maxrec,inirec
        f8_nc(:,:) = 0.0
        call rd_out2nc8(nto,mto,irec,
     &                  f8_nc,
     &                  time3,  !HYCOM time
     &                  name_t,flnm_t)
        call zhflsh(lp)
        nstep = int(time3(3))*24 + nint((time3(3)-int(time3(3)))*24.0)  !number of hours
        write(lp,*) 
        write(lp,*) 'rd_out2nc8, time = ',time3(3),nstep
        write(lp,*) 'rd_out2nc8, fld  = ',minval(f8_nc(:,:)),
     &                                    maxval(f8_nc(:,:))
        call zhflsh(lp)
c
c ---   convert to hycom arrays.
c
        do j= 1,jja
          do i= 1,ii
            field(i,j) = f8_nc(i,j)
          enddo !i
        enddo !j
        if     (larctic) then  !assume p-grid scalar field
          do i= 1,ii
            ia = ii-mod(i-1,ii)
            field(i,jj) = field(ia,jj-1)
          enddo !i
        endif
        write(lp,*) 'convert    fld  = ',minval(field(:,:)),
     &                                   maxval(field(:,:))
c
c ---   write the .[ab] files
c
        call zaiowr(field, ip,.false., xmin,xmax, 21, .false.)
        write(21,'(a,a,f14.4,2g15.6)') 
     &    trim(name_t),': day,min,max =',time3(3),xmin,xmax
        call flush(21)
        write(lp,'(a,a,f14.4,2g15.6)') 
     &    trim(name_t),': day,min,max =',time3(3),xmin,xmax
        call flush(lp)
      enddo !irec
      end program mom6nc82field
