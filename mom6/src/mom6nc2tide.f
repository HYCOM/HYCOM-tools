      program mom6nc2tide
      use mod_mom6  ! HYCOM mom6 array interface
      use mod_za    ! HYCOM array I/O interface
c
      implicit none
c
c --- convert MOM6 p-grid tide-related 2-D output to a HYCOM .[ab] file.
c --- requires the HYCOM regional.grid for the MOM6 domain.
c
c --- if environment variable MOM6_STERIC is an .a file, MOM6 input must
c --- contain col_height and bottom pressure (pbo) and the HYCOM output
c --- archive file will contain bottom pressure anomaly as a surface height
c --- in montg and steric ssh, based on column density anomaly, in steric.
c
      character*256    flnm_o,
     &                 flnm_e,flnm_r,
     &                 name_e,name_b,name_u,name_v
      character*14     c_ydh
c
      logical          larctic,lsymetr
      logical          lsteric,icegln,trcout
      integer          i,im1,ip1,irec,j,jja,jm1,jp1,k,l,mro,nvo
      integer          artype,iexpt,iversn,yrflag,itest,jtest
      integer          jday,ihour,iyear
      real             onem,thref,thbase,denij
      double precision time3(3)
      double precision area,onemm,sum1,sum2,sum3,sum4
      real             smn1,smx1,smn3,smx3,smn2,smx2
c
      real, allocatable :: scp2(:,:),plat(:,:),den_mn(:,:),ssh_mn(:,:)
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
      flnm_r = ' '
      call getenv('MOM6_STERIC',flnm_r)
      lsteric   = flnm_r.ne.' '
c
c --- 'flnm_X' ends in ".nc" for a single netCDF file, or
c ---                  ".nc.DDDD-EEEE" for subregion netCDF files
c
c --- 'flnm_e' = name of mom6  ssh  file
c --- 'name_e' = name of mom6  ssh  field (lsteric: hbt field)
c --- 'name_b' = name of mom6  pbt  field (lsteric only)
c --- 'name_u' = name of mom6  ubt  field (enter ZERO for zero field)
c --- 'name_v' = name of mom6  vbt  field (enter ZERO for zero field)
c --- 'flnm_o' = name of hycom surface archive file (output)
c --- 'in_rec' = time record to read in
c --- 'iexpt ' = experiment number x10
c --- 'itest ' = i-index for debugging printout
c --- 'jtest ' = j-index for debugging printout
c
      read (*,'(a)') flnm_e
      i = len_trim(flnm_e)
      if     (flnm_e(i-2:i).eq.'.nc') then
        write (lp,'(2a)') ' input MOM6  ssh  file: ',trim(flnm_e)
      else
        write (lp,'(2a)') ' input MOM6  ssh  files: ',flnm_e(1:i-5)
        write (lp,'(3a)') '                     to: ',flnm_e(1:i-9),
     &                                                flnm_e(i-3:i)
      endif
      call flush(lp)
      if     (lsteric) then
        read (*,'(a)') name_e
        write (lp,'(2a)') ' input MOM6  hbt  field: ',trim(name_e)
        read (*,'(a)') name_b
        write (lp,'(2a)') ' input MOM6  pbt  field: ',trim(name_b)
      else
        read (*,'(a)') name_e
        write (lp,'(2a)') ' input MOM6  ssh  field: ',trim(name_e)
      endif
      read (*,'(a)') name_u
      write (lp,'(2a)') ' input MOM6  ubt  field: ',trim(name_u)
      read (*,'(a)') name_v
      write (lp,'(2a)') ' input MOM6  vbt  field: ',trim(name_v)
      call flush(lp)
c
      read (*,'(a)') flnm_o
      write (lp,'(2a)') 'output HYCOM archive file: ',trim(flnm_o)
      call flush(lp)
c
      call blkini(irec,  'in_rec')
      call blkini(iexpt, 'iexpt ')
      call blkini(itest, 'itest ')
      call blkini(jtest, 'jtest ')
c
      thbase = 34.0
      yrflag = 3
      sigver = 6    !equation of state is (nominally) 17-term sigma-2
      artype = 1    !snapshot archive
c
c --- mom6 dimensions
c
      kk = 0
      call rd_dimen2(nto,mto,mro, flnm_e,name_e)
      if     (trim(name_v) .eq. 'ZERO') then
        mvo = nto !mvo is not used
      else
        call rd_dimen2(nvo,mvo,mro, flnm_e,name_v)
      endif
c
      lsymetr = mto .lt. mvo
      larctic = mto .eq. jdm-1 .and. nto .eq. idm
c
      write(lp,*) 
      write(lp,*) 'nto,mto,mvo,kk = ',nto,mto,mvo,kk
      write(lp,*) 'ii, jj,     kk = ',ii, jj,   0,kk
      write(lp,*) 'lsymetr = ',lsymetr
      write(lp,*) 'larctic = ',larctic
      write(lp,*) 
      call zhflsh(lp)
c
      if     (irec.gt.mro) then
        write(lp,*) 
        write(lp,*) 'irec .gt. time dimension = ',mro
        write(lp,*) 
        call zhflsh(lp)
        stop
      endif
c
c --- array allocation
c
      call mom6_alloc_tide
c
c --- cell area
c
      allocate( scp2(ii,jj), plat(ii,jj) )
      call rd_scp2(ii,jj,scp2,plat, steric)  !steric is workspace
c
c --- bathymetry.
c
      call rd_bathy(ii,jj,steric)  !steric is workspace
c
      depths(1:ii,1:jj) = steric(1:ii,1:jj)
      depths(0,   1:jj) = steric(  ii,1:jj)  !periodic
      depths(0:ii,0)    = 0.0                !land to south
c
      call bigrid
c
      if     (lsteric) then
c
c ---    ssh_mn and den_mn
c
        allocate( den_mn(ii,jj), ssh_mn(ii,jj) )
        call rd_steric(ii,jj,ssh_mn,den_mn, flnm_r)
      endif
c
c --- read the mom6 file, convert to HYCOM
c
      f_nc(:,:) = 1.0e20
      call rd_out2nc(nto,mto,irec,
     &               f_nc,
     &               time3,  !HYCOM time
     &               name_e,flnm_e)
      call zhflsh(lp)
      call m2h_p(f_nc,nto,mto,1, srfht,idm,jdm,lsymetr,larctic)  !ssh or hbt
c
      nstep =                int(time3(3)) *24 +
     &        nint((time3(3)-int(time3(3)))*24.0)  !number of hours
      write(lp,*)
      write(lp,*) 'rd_out2nc, time,nstep = ',time3(3),nstep
      call zhflsh(lp)
c
      if     (lsteric) then
        f_nc(:,:) = 1.0e20
        call rd_out2nc(nto,mto,irec,
     &                 f_nc,
     &                 time3,  !HYCOM time
     &                 name_b,flnm_e)
        call zhflsh(lp)
        call m2h_p(f_nc,nto,mto,1, steric,idm,jdm,lsymetr,larctic)  !pbt
      endif
c
      if     (trim(name_u) .eq. 'ZERO') then
        f_nc(:,:) =  0.0
      else
        f_nc(:,:) = 1.0e20
        call rd_out2nc(nto,mto,irec,
     &                 f_nc,
     &                 time3,  !HYCOM time
     &                 name_u,flnm_e)
        call zhflsh(lp)
      endif
      call m2h_u(f_nc,nto,mto,1, ubaro,idm,jdm,lsymetr,larctic)
c
      if     (trim(name_v) .eq. 'ZERO') then
        f_nc(:,:) =  0.0
      else
        f_nc(:,:) = 1.0e20
        call rd_out2nc(nto,mvo,irec,
     &                 f_nc,
     &                 time3,  !HYCOM time
     &                 name_v,flnm_e)
        call zhflsh(lp)
      endif
      call m2h_v(f_nc,nto,mvo,1, vbaro,idm,jdm,lsymetr,larctic)
c
      onem   = 9806.0  ! HYCOM mks pressure units
      thref  = 1.0e-3  ! HYCOM 1/rho_0
      do j= 1,jj
        do i= 1,ii
          if     (ip(i,j).eq.1) then
            if     (lsteric) then
c ---         convert hbt and pbt to ssh and steric ssh (m)
c ---         in addition replace montg with bottom pressure anomaly (m)
c ---         both assume that ssh_mn and den_mn are entirely steric
c
*
               montg(i,j) = steric(i,j)/(9.8*den_mn(i,j)) -  !bpa=pbt/(g*rho)-h.mn
     &                     (depths(i,j) + ssh_mn(i,j))
                    denij = steric(i,j)/(9.8* srfht(i,j))    !den=pbt/(g*hbt)
*
                  if     (i.eq.itest .and. j.eq.jtest) then
                    write(lp,'(a,2i5,3g19.8)')
     &                'i,j,_mn  = ',i,j,ssh_mn(i,j),den_mn(i,j),
     &                                              depths(i,j)
                    write(lp,'(a,2i5,3g19.8)')
     &                'i,j,hbt  = ',i,j,steric(i,j),denij,
     &                                   srfht(i,j)
                  endif !debug
*
              steric(i,j) = (denij      - den_mn(i,j))*         !sssh=(den.a*
     &                     (depths(i,j) + ssh_mn(i,j))/denij +  !      h.mn)/den+
     &                      ssh_mn(i,j)                         !     ssh.mn
               srfht(i,j) =  srfht(i,j) - depths(i,j)           ! ssh=hbt-D
*
                  if     (i.eq.itest .and. j.eq.jtest) then
                    write(lp,'(a,2i5,2g19.8)')
     &                'i,j,sssh = ',i,j,steric(i,j),denij-den_mn(i,j)
                    write(lp,'(a,2i5,2g19.8)')
     &                'i,j, ssh = ',i,j, srfht(i,j),montg(i,j)
                  endif !debug
*
            endif
             montg(i,j)   =  montg(i,j) *onem*thref
            steric(i,j)   = steric(i,j) *onem*thref
             srfht(i,j)   =  srfht(i,j) *onem*thref
          else
             montg(i,j)   = spval
            steric(i,j)   = spval
             srfht(i,j)   = spval
          endif !ip
          if     (iu(i,j).eq.0) then
            ubaro(i,j) = 0.0
          endif !iu
          if     (iv(i,j).eq.0) then
            vbaro(i,j) = 0.0
          endif !iu
        enddo !i
      enddo !j
*
      if     (min(itest,jtest).gt.0) then
        i = itest
        j = jtest
        write(lp,'(a,2i5,3g19.8)')
     &    'i,j,ssh  = ',i,j, srfht(i,j)/(onem*thref),
     &                       montg(i,j)/(onem*thref),
     &                      steric(i,j)/(onem*thref)
      endif !debug
*
c
c --- write the archive file
c
      iversn    = 21
      ctitle(1) = 'from MOM6 netCDF files'
      ctitle(2) = 'converted by mom6nc2tide'
      if     (lsteric) then
        ctitle(3) = 'montg1 is bottom pressure anomaly'
        ctitle(4) = ' '
      else
        ctitle(3) = ' '
        ctitle(4) = ' '
      endif
      call putdat_tide(flnm_o,artype,time3,lsteric,
     &                 iexpt,iversn,yrflag, thbase)
c
c --- write out area-mean statistics
c
      l = len_trim(flnm_o)
      open (unit=11,file=flnm_o(1:l-2)//'.stat',
     &      form='formatted',status='new',action='write')
      write(lp,'(a,a)') 'open: ',flnm_o(1:l-2)//'.stat'
c
      call forday(time3(3),yrflag, iyear,jday,ihour)
      write(c_ydh,'('' ('',i4.4,''/'',i3.3,1x,i2.2,'')'')')
     &  iyear,jday,ihour
c
      if     (larctic) then
        jja = jj-1
      else
        jja = jj
      endif
      area =  0.d0
      sum1 =  0.d0
      smn1 =  huge(smn1)
      smx1 = -huge(smx1)
      sum2 =  0.d0
      smn2 =  huge(smn2)
      smx2 = -huge(smx2)
      sum3 =  0.d0
      smn3 =  huge(smn3)
      smx3 = -huge(smx3)
      do j= 1,jja
        do i= 1,ii
          if     (ip(i,j).eq.1) then
            area = area + scp2(i,j)
            sum1 = sum1 + scp2(i,j)*srfht(i,j)
            smn1 = min( smn1,       srfht(i,j) )
            smx1 = max( smx1,       srfht(i,j) )
            if     (lsteric) then
              sum2 = sum2 + scp2(i,j)* montg(i,j)
              smn2 = min( smn2,        montg(i,j) )
              smx2 = max( smx2,        montg(i,j) )
              sum3 = sum3 + scp2(i,j)*steric(i,j)
              smn3 = min( smn3,       steric(i,j) )
              smx3 = max( smx3,       steric(i,j) )
            endif !lsteric
          endif !ip
        enddo !i
      enddo !j
      onemm = 9806.0d0*0.001d0  ! HYCOM mks pressure units to mm
      write (11,'(i9,a,
     &              '' mean      SSH (mm):'',f8.2,
     &              ''  ('',1pe8.1,'' to '',e8.1,'')'')')
     &  nstep,c_ydh,
     &  sum1/(area*thref*onemm),smn1/(thref*onemm),smx1/(thref*onemm)
      write (lp,'(i9,a,
     &              '' mean      SSH (mm):'',f8.2,
     &              ''  ('',1pe8.1,'' to '',e8.1,'')'')')
     &  nstep,c_ydh,
     &  sum1/(area*thref*onemm),smn1/(thref*onemm),smx1/(thref*onemm)
      call flush(lp)
      if     (lsteric) then
        write (11,'(i9,a,
     &                '' mean BOTPRS.A (mm):'',f8.2,
     &                ''  ('',1pe8.1,'' to '',e8.1,'')'')')
     &    nstep,c_ydh,
     &    sum2/(area*thref*onemm),smn2/(thref*onemm),smx2/(thref*onemm)
        write (lp,'(i9,a,
     &                '' mean BOTPRS.A (mm):'',f8.2,
     &                ''  ('',1pe8.1,'' to '',e8.1,'')'')')
     &    nstep,c_ydh,
     &    sum2/(area*thref*onemm),smn2/(thref*onemm),smx2/(thref*onemm)
        call flush(lp)
        write (11,'(i9,a,
     &                '' mean SSH-BP.A (mm):'',f8.2)')
     &    nstep,c_ydh,
     &    (sum1-sum2)/(area*1.d-3*onemm)
        write (lp,'(i9,a,
     &                '' mean SSH-BP.A (mm):'',f8.2)')
     &    nstep,c_ydh,
     &    (sum1-sum2)/(area*1.d-3*onemm)
        call flush(lp)
        write (11,'(i9,a,
     &                '' mean    S.SSH (mm):'',f8.2,
     &                ''  ('',1pe8.1,'' to '',e8.1,'')'')')
     &    nstep,c_ydh,
     &    sum3/(area*1.d-3*onemm),smn3/(1.d-3*onemm),smx3/(1.d-3*onemm)
        write (lp,'(i9,a,
     &                '' mean    S.SSH (mm):'',f8.2,
     &                ''  ('',1pe8.1,'' to '',e8.1,'')'')')
     &    nstep,c_ydh,
     &    sum3/(area*1.d-3*onemm),smn3/(1.d-3*onemm),smx3/(1.d-3*onemm)
        call flush(lp)
        write (11,'(i9,a,
     &                '' mean SSH-SSSH (mm):'',f8.2)')
     &    nstep,c_ydh,
     &    (sum1-sum3)/(area*1.d-3*onemm)
        write (lp,'(i9,a,
     &                '' mean SSH-SSSH (mm):'',f8.2)')
     &    nstep,c_ydh,
     &    (sum1-sum3)/(area*1.d-3*onemm)
        call flush(lp)
      endif !lsteric
c
      close(unit=11)
c
      end program mom6nc2tide

      subroutine m2h_p(f_nc,nto,mto,kk, field,ii,jj,lsymetr,larctic)
      implicit none
c
      logical lsymetr,larctic
      integer nto,mto,kk,ii,jj
      real    f_nc(nto,mto,kk),field(ii,jj,kk)
c
c --- spval  = hycom data void marker, 2^100 or about 1.2676506e30
      real, parameter :: spval=2.0**100
c
c --- convert p-grid mom6 array to hycom.
c
      integer i,ia,j,k
c
      do k= 1,kk
        do j= 1,mto
          do i= 1,nto
            if     (f_nc(i,j,k).ne.1.0e20) then
               field(i,j,k) = f_nc(i,j,k)
            else
               field(i,j,k) = spval
            endif
          enddo !i
        enddo !j
        if     (larctic) then  !p-grid scalar field, mto=jj-1
          do i= 1,nto
            ia = nto-mod(i-1,nto)
            field(i,jj,k) = field(ia,jj-1,k)
          enddo !i
        elseif (lsymetr) then
          do i= 1,nto
            field(i,jj,k) = spval
          enddo !i
          do j= 1,jj
            field(ii,j,k) = spval
          enddo !j
        endif !larctic:lsymetr
      enddo !k
      return
      end

      subroutine m2h_u(f_nc,nto,mto,kk, field,ii,jj,lsymetr,larctic)
      implicit none
c
      logical lsymetr,larctic
      integer nto,mto,kk,ii,jj
      real    f_nc(nto,mto,kk),field(ii,jj,kk)
c
c --- convert u-grid mom6 array to hycom.
c --- nonsymetric mom6  has "q" at i+0.5,j+0.5 w.r.t. p.ij
c ---    symetric mom6  has "q" at i-0.5,j-0.5 w.r.t. p.ij
c ---             hycom has "q" at i-0.5,j-0.5 w.r.t. p.ij
c
c --- spval  = hycom data void marker, 2^100 or about 1.2676506e30
      real, parameter :: spval=2.0**100
c
      integer i,ia,j,k
c
      if     (lsymetr) then
        do k= 1,kk
          do j= 1,mto
            do i= 1,nto
              if     (f_nc(i,j,k).ne.1.0e20) then
                 field(i,j,k) = f_nc(i,j,k)
              else
                 field(i,j,k) = spval
              endif
            enddo !i
          enddo !j
          if     (larctic) then  !u-grid vector field, mto=jj-1
            do i= 1,nto
              ia = mod(nto-(i-1),nto)+1 
              if     (field(ia,jj-1,k).ne.spval) then
                field(i,jj,k) = -field(ia,jj-1,k)
              else
                 field(i,j,k) = spval
              endif
            enddo !i
          else
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
              ia = mod(i,nto)+1
              if     (f_nc(i,j,k).ne.1.0e20) then
                 field(ia,j,k) = f_nc(i,j,k)
              else
                 field(ia,j,k) = spval
              endif
            enddo !i
          enddo !j
          if     (larctic) then  !u-grid vector field, mto=jj-1
            do i= 1,nto
              ia = mod(nto-(i-1),nto)+1 
              if     (field(ia,jj-1,k).ne.spval) then
                field(i,jj,k) = -field(ia,jj-1,k)
              else
                 field(i,j,k) = spval
              endif
            enddo !i
          endif
        enddo !k
      endif !lsymetr:else
      return
      end

      subroutine m2h_v(f_nc,nto,mvo,kk, field,ii,jj,lsymetr,larctic)
      implicit none
c
      logical lsymetr,larctic
      integer nto,mvo,kk,ii,jj
      real    f_nc(nto,mvo,kk),field(ii,jj,kk)
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
          do j= 1,min(mvo,jj)
            do i= 1,nto
              if     (f_nc(i,j,k).ne.1.0e20) then
                 field(i,j,k) = f_nc(i,j,k)
              else
                 field(i,j,k) = spval
              endif
            enddo !i
          enddo !j
          do j= mvo+1,jj
            do i= 1,nto  !must be land
              field(i,j,k) = spval
            enddo !i
          enddo !j
          do i= nto+1,ii
            do j= 1,jj   !must be land
              field(i,j,k) = spval
            enddo !j
          enddo !i
        enddo !k
      else
        do k= 1,kk
          do j= 1,min(mvo,jj-1)
            do i= 1,nto
              if     (f_nc(i,j,k).ne.1.0e20) then
                 field(i,j+1,k) = f_nc(i,j,k)
              else
                 field(i,j+1,k) = spval
              endif
            enddo !i
          enddo !j
          do i= 1,nto  !must be land
            field(i,1,k) = spval
          enddo !i
        enddo !k
      endif !lsymetr:else
      return
      end
