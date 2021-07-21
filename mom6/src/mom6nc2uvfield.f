      program mom6nc2uvfield
      use mod_mom6  ! HYCOM mom6 array interface
      use mod_za    ! HYCOM array I/O interface
c
      implicit none
c
c --- convert MOM6 u- and v-grid 2-D output to HYCOM p-grid .[ab] files.
c --- requires the HYCOM regional.grid for the MOM6 domain.
c
      character*256    flnm_u,flnm_i,name_u,
     &                 flnm_v,       name_v
      logical          larctic,lsymetr
      integer          i,ia,irec,j,nuo,mro
      integer          inirec,increc,maxrec
      integer          yrflag
      real             xmin,xmax,misval
      double precision time3(3)
c
      call xcspmd  !input the HYCOM array dimensions
      call zaiost  !initialize HYCOM I/O
      lp     = 6
      ii     = idm
      jj     = jdm
      yrflag = 3
c
c --- 'name_u' = name of mom6  u-vel field
c --- 'name_v' = name of mom6  v-vel field
c --- 'flnm_i' = name of mom6    vel file  (input)
c ---              ends in ".nc" for a single netCDF file, or
c ---                      ".nc.DDDD-EEEE" for subregion netCDF files
c --- 'flnm_u' = name of hycom u-vel file (output)
c --- 'flnm_v' = name of hycom v-vel file (output)
c --- 'inirec' = initial time record to read in
c --- 'maxrec' = final   time record to read in, 0 for to the end
c --- 'increc' =         time record increment
c
      read (*,'(a)') name_u
      write (lp,'(2a)') ' input MOM6 u-veloc: ',trim(name_u)
      call flush(lp)
      read (*,'(a)') name_v
      write (lp,'(2a)') ' input MOM6 v-veloc: ',trim(name_v)
      call flush(lp)
      read (*,'(a)') flnm_i
      i = len_trim(flnm_i)
      if     (flnm_i(i-2:i).eq.'.nc') then
        write (lp,'(2a)') ' input MOM6    file: ',trim(flnm_i)
      else
        write (lp,'(2a)') ' input MOM6   files: ',flnm_i(1:i-5)
        write (lp,'(3a)') '                 to: ',flnm_i(1:i-9),
     &                                          flnm_i(i-3:i)
      endif
      call flush(lp)
      read (*,'(a)') flnm_u
      write (lp,'(2a)') 'output HYCOM u-file: ',trim(flnm_u)
      call flush(lp)
      read (*,'(a)') flnm_v
      write (lp,'(2a)') 'output HYCOM v-file: ',trim(flnm_v)
      call flush(lp)
c
      call blkini(inirec,  'inirec')
      call blkini(maxrec,  'maxrec')
      call blkini(increc,  'increc')
c
c --- mom6 dimensions
c
      call rd_dimen2(nuo,mto,mro, flnm_i,name_u)
      call rd_dimen2(nto,mvo,mro, flnm_i,name_v)
      call rd_missing(misval,     flnm_i,name_u)
c
      lsymetr = mto .lt. mvo
      larctic = mto .eq. jdm-1 .and. nto .eq. idm
c
      write(lp,*) 
      write(lp,*) 'nto,mto = ',nto,mto
      write(lp,*) 'nuo,mvo = ',nuo,mvo
      write(lp,*) 'ii,jj   = ',ii, jj
      write(lp,*) 'lsymetr = ',lsymetr
      write(lp,*) 'larctic = ',larctic
      write(lp,*) 
      write(lp,*) 'misval  = ',misval
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
      call mom6_alloc_field
c
c     open the output .[ab] files.
c
      i = len_trim(flnm_u)
      call zaiopf(      flnm_u(1:i-2)//'.a',       'new',21)
      open(unit=21,file=flnm_u(1:i-2)//'.b',status='new',
     &     form='formatted',action='write')
      i = len_trim(flnm_v)
      call zaiopf(      flnm_v(1:i-2)//'.a',       'new',22)
      open(unit=22,file=flnm_v(1:i-2)//'.b',status='new',
     &     form='formatted',action='write')
c
c --- read the mom6 file.
c
      do irec= inirec,maxrec,inirec
c
c ---   u-vel
c
        f_nc(:,:) = misval
        call rd_out2nc(nto,mto,irec,
     &                 f_nc,
     &                 time3,  !HYCOM time
     &                 name_u,flnm_i)
        call zhflsh(lp)
        nstep = int(time3(3))*24 + nint((time3(3)-int(time3(3)))*24.0)  !number of hours
        write(lp,*) 
        write(lp,*) 'rd_out2nc, time = ',time3(3),nstep
        write(lp,*) 'rd_out2nc, fld  = ',minval(f_nc(:,:)),
     &                                   maxval(f_nc(:,:))
        call zhflsh(lp)
C
        call m2h_u(f_nc,nto,mto,1,misval, field,idm,jdm,lsymetr,larctic)
        write(lp,*) 'convert    fld  = ',minval(field(:,:)),
     &                                   maxval(field(:,:))
c
c ---   write the .[ab] files
c
        call zaiowr(field, ip,.false., xmin,xmax, 21, .false.)
        write(21,'(a,a,f14.4,2g15.6)') 
     &    trim(name_u),': day,min,max =',time3(3),xmin,xmax
        call flush(21)
        write(lp,'(a,a,f14.4,2g15.6)') 
     &    trim(name_u),': day,min,max =',time3(3),xmin,xmax
        call flush(lp)
c
c ---   v-vel
c
        v_nc(:,:,1) = misval
        call rd_out2nc(nto,mvo,irec,
     &                 v_nc,
     &                 time3,  !HYCOM time
     &                 name_v,flnm_i)
        call zhflsh(lp)
        nstep = int(time3(3))*24 + nint((time3(3)-int(time3(3)))*24.0)  !number of hours
        write(lp,*) 
        write(lp,*) 'rd_out2nc, time = ',time3(3),nstep
        write(lp,*) 'rd_out2nc, fld  = ',minval(v_nc(:,:,1)),
     &                                   maxval(v_nc(:,:,1))
        call zhflsh(lp)
C
        call m2h_v(f_nc,nto,mvo,1,misval, field,idm,jdm,lsymetr,larctic)
        write(lp,*) 'convert    fld  = ',minval(field(:,:)),
     &                                   maxval(field(:,:))
c
c ---   write the .[ab] files
c
        call zaiowr(field, ip,.false., xmin,xmax, 22, .false.)
        write(22,'(a,a,f14.4,2g15.6)') 
     &    trim(name_v),': day,min,max =',time3(3),xmin,xmax
        call flush(22)
        write(lp,'(a,a,f14.4,2g15.6)') 
     &    trim(name_v),': day,min,max =',time3(3),xmin,xmax
        call flush(lp)
      enddo !irec
      close (unit=21)
      call zaiocl(21)
      close (unit=22)
      call zaiocl(22)
      end program mom6nc2uvfield

      subroutine m2h_u(f_nc,nto,mto,kk,misval,
     &                 field,ii,jj,lsymetr,larctic)
      implicit none
c
      logical lsymetr,larctic
      integer nto,mto,kk,ii,jj
      real    f_nc(nto,mto,kk),misval,field(ii,jj,kk)
c
c --- convert u-grid mom6 array to hycom p-grid.
c --- nonsymetric mom6  has "q" at i+0.5,j+0.5 w.r.t. p.ij
c ---    symetric mom6  has "q" at i-0.5,j-0.5 w.r.t. p.ij
c ---             hycom has "q" at i-0.5,j-0.5 w.r.t. p.ij
c
c --- spval  = hycom data void marker, 2^100 or about 1.2676506e30
      real, parameter :: spval=2.0**100
c
      integer i,ia,ip1,j,k
c
      if     (lsymetr) then
        do k= 1,kk
          do j= 1,mto
            do i= 1,nto
              ip1 = mod(i,nto)+1
              if     (f_nc(i,  j,k).ne.misval .and.
     &                f_nc(ip1,j,k).ne.misval      ) then
                 field(i,j,k) = 0.5*(f_nc(i,j,k)+f_nc(ip1,j,k))
              else
                 field(i,j,k) = spval
              endif
            enddo !i
          enddo !j
          if     (larctic) then  !p-grid vector field, mto=jj-1
            do i= 1,nto
              ia = nto-mod(i-1,nto)
              if     (field(ia,jj-1,k).ne.spval) then
                field(i,jj,k) = -field(ia,jj-1,k)
              else
                field(i,jj,k) = spval
              endif
            enddo !i
          else !lsymetr
            do j= mto+1,jj
              do i= 1,nto  !must be land
                field(i,j,k) = spval
              enddo !i
            enddo !j
            do i= nto+1,ii
              do j= 1,jj   !must be land
                field(i,j,k) = spval
              enddo !j
            enddo !i
          endif !larctic:else
        enddo !k
      else
        do k= 1,kk
          do j= 1,mto
            do i= 1,nto
              ia  = mod(i,nto)+1
              ip1 = mod(i,nto)+1
              if     (f_nc(i,  j,k).ne.misval .and.
     &                f_nc(ip1,j,k).ne.misval      ) then
                 field(ia,j,k) = 0.5*(f_nc(i,j,k)+f_nc(ip1,j,k))
              else
                 field(ia,j,k) = spval
              endif
            enddo !i
          enddo !j
          if     (larctic) then  !p-grid vector field, mto=jj-1
            do i= 1,nto
              ia = nto-mod(i-1,nto)
              if     (field(ia,jj-1,k).ne.spval) then
                field(i,jj,k) = -field(ia,jj-1,k)
              else
                field(i,jj,k) = spval
              endif
            enddo !i
          endif
        enddo !k
      endif !lsymetr:else
      return
      end

      subroutine m2h_v(f_nc,nto,mto,kk,misval,
     &                 field,ii,jj,lsymetr,larctic)
      implicit none
c
      logical lsymetr,larctic
      integer nto,mto,kk,ii,jj
      real    f_nc(nto,mto,kk),misval,field(ii,jj,kk)
c
c --- convert v-grid mom6 array to hycom.
c --- mom6  standard has "q" at i+0.5,j+0.5 w.r.t. p.ij
c --- mom6  symetric has "q" at i-0.5,j-0.5 w.r.t. p.ij
c --- hycom          has "q" at i-0.5,j-0.5 w.r.t. p.ij
c
c --- spval  = hycom data void marker, 2^100 or about 1.2676506e30
      real, parameter :: spval=2.0**100
c
      integer i,ia,j,k
c
      if     (lsymetr) then
        do k= 1,kk
          do j= 1,min(mto-1,jj)
            do i= 1,nto
              if     (f_nc(i,j,  k).ne.misval .and.
     &                f_nc(i,j+1,k).ne.misval      ) then
                 field(i,j,k) = 0.5*(f_nc(i,j,k)+f_nc(i,j+1,k))
              else
                 field(i,j,k) = spval
              endif
            enddo !i
          enddo !j
          if     (larctic) then  !p-grid vector field, mto=jj
            do i= 1,nto
              ia = nto-mod(i-1,nto)
              if     (field(ia,jj-1,k).ne.spval) then
                field(i,jj,k) = -field(ia,jj-1,k)
              else
                field(i,jj,k) = spval
              endif
            enddo !i
          else !lsymetr
            do j= mto,jj
              do i= 1,nto  !must be land
                field(i,j,k) = spval
              enddo !i
            enddo !j
            do i= nto+1,ii
              do j= 1,jj   !must be land
                field(i,j,k) = spval
              enddo !j
            enddo !i
          endif !larctic:else
        enddo !k
      else
        do k= 1,kk
          do j= 1,min(mto,jj-1)  !mto if larctic
            do i= 1,nto
              if     (f_nc(i,j,  k).ne.misval .and.
     &                f_nc(i,j+1,k).ne.misval     ) then
                 field(i,j+1,k) = 0.5*(f_nc(i,j,k)+f_nc(i,j+1,k))
              else
                 field(i,j+1,k) = spval
              endif
            enddo !i
          enddo !j
          do i= 1,nto  !must be land
            field(i,1,k) = spval
          enddo !i
          if     (larctic) then  !p-grid vector field, mto=jj-1
            do i= 1,nto
              ia = nto-mod(i-1,nto)
              if     (field(ia,jj-1,k).ne.spval) then
                field(i,jj,k) = -field(ia,jj-1,k)
              else
                field(i,jj,k) = spval
              endif
            enddo !i
          endif
        enddo !k
      endif !lsymetr:else
      return
      end
