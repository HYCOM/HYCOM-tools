      subroutine forday(dtime,yrflag, iyear,iday,ihour)
      implicit none
c
      double precision dtime
      integer          yrflag, iyear,iday,ihour
c
c --- converts model day to "calendar" date (year,julian-day,hour).
c
      double precision dtim1,day
      integer          iyr,nleap
c
      if     (yrflag.eq.0) then
c ---   360 days per model year, starting Jan 16
        iyear =  int((dtime+15.001d0)/360.d0) + 1
        iday  =  mod( dtime+15.001d0 ,360.d0) + 1
        ihour = (mod( dtime+15.001d0 ,360.d0) + 1.d0 - iday)*24.d0
c
      elseif (yrflag.eq.1) then
c ---   366 days per model year, starting Jan 16
        iyear =  int((dtime+15.001d0)/366.d0) + 1
        iday  =  mod( dtime+15.001d0 ,366.d0) + 1
        ihour = (mod( dtime+15.001d0 ,366.d0) + 1.d0 - iday)*24.d0
c
      elseif (yrflag.eq.2) then
c ---   366 days per model year, starting Jan 01
        iyear =  int((dtime+ 0.001d0)/366.d0) + 1
        iday  =  mod( dtime+ 0.001d0 ,366.d0) + 1
        ihour = (mod( dtime+ 0.001d0 ,366.d0) + 1.d0 - iday)*24.d0
c
      elseif (yrflag.eq.3) then
c ---   model day is calendar days since 01/01/1901
        iyr   = (dtime-1.d0)/365.25d0
        nleap = iyr/4
        dtim1 = 365.d0*iyr + nleap + 1.d0
        day   = dtime - dtim1 + 1.d0
        if     (dtim1.gt.dtime) then
          iyr = iyr - 1
        elseif (day.ge.367.d0) then
          iyr = iyr + 1
        elseif (day.ge.366.d0 .and. mod(iyr,4).ne.3) then
          iyr = iyr + 1
        endif
        nleap = iyr/4
        dtim1 = 365.d0*iyr + nleap + 1.d0
c
        iyear =  1901 + iyr
        iday  =  dtime - dtim1 + 1
        ihour = (dtime - dtim1 + 1.d0 - iday)*24.d0
c
      else
        stop '(forday)'
      endif
      return
      end

      subroutine getdat_name(flnm, dtime,yrflag)
      use mod_za     ! HYCOM array I/O interface
c
      implicit none
c
      character        flnm*(*)
      double precision dtime
      integer          yrflag
c
c --- modify the archive filename for the specified model day.
c --- flnm is assumed to be of the form */arch*_????_???_??.b.
c
      integer   iyear,iday,ihour,l
c
      call forday(dtime,yrflag, iyear,iday,ihour)
c
      l = len_trim(flnm)
      write(flnm(l-12:l-9),'(i4.4)') iyear
      write(flnm(l- 7:l-5),'(i3.3)') iday
      write(flnm(l- 3:l-2),'(i2.2)') ihour
      return
      end

      subroutine getdat_time(flnm, time,iexpt)
      use mod_za     ! HYCOM array I/O interface
c
      implicit none
c
      character        flnm*(*)
      double precision time
      integer          iexpt
c
c --- extract the model time (and iexpt) from an archive.
c
      character cline*80,cvarin*6
      integer   i,irec,l,ni,nstep
c
      data ni/14/
c
      l = len_trim(flnm)
      open (unit=ni,file=flnm(1:l-2)//'.b',form='formatted',
     &      status='old',action='read')
      do irec= 1,5
        read (ni,'(a)',end=6) cline
      enddo
      read( ni,*) iexpt,cvarin
      if (cvarin.ne.'iexpt ') then
        write(lp,*)
        write(lp,*) 'error in getdat - input ',cvarin,
     &                        ' but should be iexpt '
        write(lp,*)
        call flush(lp)
        stop
      endif
      do irec= 1,7  !ok for standard or mean file
        read (ni,'(a)',end=6) cline
        write(lp,*) trim(cline)
      enddo
      i = index(cline,'=')
      read (cline(i+1:),*)  nstep,time
      write(lp,*) 'getdat_time = ',time
      call flush(lp)
      close(unit=ni)
      return
c
c --- unexpected end of file
 6    continue
        write (lp,*) '***** unexpected end of archive file *****'
        write (lp,*) flnm(1:l-2)//'.b'
        call flush(lp)
        stop '(e-o-f)'
      end

      subroutine getdat(flnm,time, iexpt,yrflag, trcout)
      use mod_trans  ! HYCOM transport section archive array interface
      use mod_za     ! HYCOM array I/O interface
c
      implicit none
c
      character        flnm*(*)
      double precision time
      integer          iexpt,yrflag
      logical          trcout  !ignored
c
c --- read model fields for transport sections.
c --- HYCOM 2.0 array I/O archive file.
c
      character cline*80
      character preambl(5)*79
      character cvarin*6
      real      hminb,hmaxb,thet
      integer   i,idmtst,irec,iversn,ios,jdmtst,j,k,l,layer,ni,nstep
      integer   artype,ktr,ntr
      logical   icegln,nodens
c
      data ni/14/
c
      l = len_trim(flnm)
      open (unit=ni,file=flnm(1:l-2)//'.b',form='formatted',
     &      status='old',action='read')
      call zaiopf(flnm(1:l-2)//'.a','old', ni)
c
      read( ni,'(a)') cline
      write(lp,'(a)') trim(cline)
      read( ni,'(a)') cline
      write(lp,'(a)') trim(cline)
      read( ni,'(a)') cline
      write(lp,'(a)') trim(cline)
      read( ni,'(a)') cline
      write(lp,'(a)') trim(cline)
      call flush(lp)
      read( ni,*) iversn,cvarin
      write(lp,*) cvarin,' = ',iversn
      if (cvarin.ne.'iversn') then
        write(lp,*)
        write(lp,*) 'error in getdat - input ',cvarin,
     &                          ' but should be iversn'
        write(lp,*)
        call flush(lp)
        stop
      endif
      read( ni,*) iexpt,cvarin
      write(lp,*) cvarin,' = ',iexpt
      if (cvarin.ne.'iexpt ') then
        write(lp,*)
        write(lp,*) 'error in getdat - input ',cvarin,
     &                        ' but should be iexpt '
        write(lp,*)
        call flush(lp)
        stop
      endif
      read( ni,*) yrflag,cvarin
      write(lp,*) cvarin,' = ',yrflag
      if (cvarin.ne.'yrflag') then
        write(lp,*)
        write(lp,*) 'error in getdat - input ',cvarin,
     &                        ' but should be yrflag'
        write(lp,*)
        call flush(lp)
        stop
      endif
      read( ni,*) idmtst,cvarin
      write(lp,*) cvarin,' = ',idmtst
      if (cvarin.ne.'idm   ') then
        write(lp,*)
        write(lp,*) 'error in getdat - input ',cvarin,
     &                        ' but should be idm   '
        write(lp,*)
        call flush(lp)
        stop
      endif
      read( ni,*) jdmtst,cvarin
      write(lp,*) cvarin,' = ',jdmtst
      if (cvarin.ne.'jdm   ') then
        write(lp,*)
        write(lp,*) 'error in getdat - input ',cvarin,
     &                        ' but should be jdm   '
        write(lp,*)
        call flush(lp)
        stop
      endif
c
      if (idmtst.ne.idm .or. jdmtst.ne.jdm) then
        write(lp,*)
        write(lp,*) 'error in getdat - input idm,jdm',
     &                        ' not consistent with parameters'
        write(lp,*) 'idm,jdm = ',idm,   jdm,   '  (REGION.h)'
        write(lp,*) 'idm,jdm = ',idmtst,jdmtst,'  (input)'
        write(lp,*)
        call flush(lp)
        stop
      endif
c
c --- artype
c
      read( ni,'(a)') cline
      write(lp,'(a)') trim(cline)
      if     (cline(24:28).eq.'model') then
        artype = 1  ! standard archive
      elseif (cline(24:28).eq.' mean') then
        artype = 2  ! mean archive, see below
      else
        write(lp,*)
        write(lp,*) 'error in getdat - wrong archive type.'
        write(lp,*) 'only "model", or " mean" allowed.'
        write(lp,*)
        call flush(lp)
        stop
      endif 
      call flush(lp)
c
c --- detect version 2.2 normal archive files
c
      read (ni,'(a)',end=6) cline
      write(lp,'(a)')  trim(cline)
      call flush(lp)
      i = index(cline,'=')
      read (cline(i+1:),*)  nstep,time,layer,thet,hminb,hmaxb
      nodens = layer.ne.0
      call zaiosk(ni)
c
c --- skip surface fields, except ubaro and vbaro.
c
      do irec= 2,9999
        read (ni,'(a)',end=6) cline
        write(lp,'(a)')  trim(cline)
        if     (cline(1:8).eq.'u_btrop') then
          call flush(lp)
          i = index(cline,'=')
          read (cline(i+1:),*)  nstep,time,layer,thet,hminb,hmaxb
          call getfld(ubaro,    ni, hminb,hmaxb)
          write(lp,'("input  ",a," into ",a)') cline(1:8),'ubaro'
c
          read (ni,'(a)',end=6) cline
          write(lp,'(a)')  trim(cline)
          call flush(lp)
          i = index(cline,'=')
          read (cline(i+1:),*)  nstep,time,layer,thet,hminb,hmaxb
          call getfld(vbaro,    ni, hminb,hmaxb)
          write(lp,'("input  ",a," into ",a)') cline(1:8),'vbaro'
c
          if     (artype.eq.2) then  !kebaro
            read (ni,'(a)',end=6) cline
            write(lp,'(a)')  trim(cline)
            call zaiosk(ni)
          endif
c
          exit
        else
          call zaiosk(ni)
        endif
      enddo
c
      do k=1,kk
        if     (k.eq.2) then
c ---     already input at end of k=1 loop.
        else
          read (ni,'(a)',end=6)   cline
        endif
        write(lp,'(a)')    trim(cline)
        call flush(lp)
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time,layer,thet,hminb,hmaxb
        if     (cline(1:8).ne.'u-vel.  ') then
          write(lp,*)
          write(lp,*) 'error in getdat - layer ',k,
     &               ' does not exist (kk= ',kk,')'
          write(lp,*)
          call flush(lp)
          stop
        endif
        call getfld(u(1,1,k), ni, hminb,hmaxb)
        if     (k.eq.1 .or. k.eq.kk) then
          write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'u ',k
        endif
        call flush(lp)
c
        read (ni,'(a)',end=6) cline
        write(lp,'(a)')  trim(cline)
        call flush(lp)
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time,layer,thet,hminb,hmaxb
        call getfld(v(1,1,k), ni, hminb,hmaxb)
        if     (k.eq.1 .or. k.eq.kk) then
          write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'v ',k
        endif
        call flush(lp)
c
        if     (artype.eq.2) then
c
c ---     mean archive, skip ke and convert from total to baroclinic velocity
c
          read (ni,'(a)',end=6) cline
          call zaiosk(ni)
c
          do j= 1,jj
            do i= 1,ii
               u(i,j,k) = u(i,j,k) - ubaro(i,j)
               v(i,j,k) = v(i,j,k) - vbaro(i,j)
            enddo
          enddo
        endif
c
        read (ni,'(a)',end=6) cline
        write(lp,'(a)')  trim(cline)
        call flush(lp)
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time,layer,thet,hminb,hmaxb
        call getfld(dp(1,1,k), ni, hminb,hmaxb)
        if     (k.eq.1 .or. k.eq.kk) then
          write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'dp',k
        endif
        call flush(lp)
c
        if     (nodens) then
          do irec= 1,2 !temp,saln
            read (ni,'(a)',end=6) cline
            call zaiosk(ni)
          enddo
        else
          do irec= 1,3 !temp,saln,th3d
            read (ni,'(a)',end=6) cline
            call zaiosk(ni)
          enddo
        endif !nodens:else
        if     (kk.ne.1) then
          if     (k.eq.1) then
            do ktr= 1,999
              read (ni,'(a)',end=6) cline
              if     (cline(1:8).ne.'tracer  ' .and.
     &                cline(1:8).ne.'viscty  ' .and.
     &                cline(1:8).ne.'t-diff  ' .and.
     &                cline(1:8).ne.'s-diff  '      ) then
                exit !end of tracers and visc/diff
              else
                call zaiosk(ni)
              endif
            enddo !ktr
            ntr=ktr-1
          else
            do ktr= 1,ntr
              read (ni,'(a)',end=6) cline
              call zaiosk(ni)
            enddo
          endif  !k.eq.1:else
        endif  !kk.ne.1
        theta(k)=thet
      enddo  !k=1,kk
c
      close( unit=14)
      call zaiocl(14)
      return
c
c --- unexpected end of file
 6    continue
        write (lp,*) '***** unexpected end of archive file *****'
        write (lp,*) flnm(1:l-2)//'.b'
        call flush(lp)
        stop '(e-o-f)'
      end

      subroutine getdat_layer(flnm,time, kn, iexpt,yrflag, trcout)
      use mod_trans  ! HYCOM transport section archive array interface
      use mod_za     ! HYCOM array I/O interface
c
      implicit none
c
      character        flnm*(*)
      double precision time
      integer          kn,iexpt,yrflag
      logical          trcout  !ignored
c
c --- read single layer model fields for transport sections.
c --- HYCOM 2.0 array I/O archive file.
c
      character cline*80
      character preambl(5)*79
      character cvarin*6
      real      hminb,hmaxb,thet
      integer   i,idmtst,irec,iversn,ios,jdmtst,j,k,l,layer,ni,nstep
      integer   artype,ktr,ntr
      logical   icegln,nodens
c
      data ni/14/
c
      l = len_trim(flnm)
      open (unit=ni,file=flnm(1:l-2)//'.b',form='formatted',
     &      status='old',action='read')
      call zaiopf(flnm(1:l-2)//'.a','old', ni)
c
      read( ni,'(a)') cline
      if(kn.eq.1) write(lp,'(a)') trim(cline)
      read( ni,'(a)') cline
      if(kn.eq.1) write(lp,'(a)') trim(cline)
      read( ni,'(a)') cline
      if(kn.eq.1) write(lp,'(a)') trim(cline)
      read( ni,'(a)') cline
      if(kn.eq.1) write(lp,'(a)') trim(cline)
      call flush(lp)
      read( ni,*) iversn,cvarin
      if(kn.eq.1) write(lp,*) cvarin,' = ',iversn
      if (cvarin.ne.'iversn') then
        write(lp,*)
        write(lp,*) 'error in getdat - input ',cvarin,
     &                          ' but should be iversn'
        write(lp,*)
        call flush(lp)
        stop
      endif
      read( ni,*) iexpt,cvarin
      if(kn.eq.1) write(lp,*) cvarin,' = ',iexpt
      if (cvarin.ne.'iexpt ') then
        write(lp,*)
        write(lp,*) 'error in getdat - input ',cvarin,
     &                        ' but should be iexpt '
        write(lp,*)
        call flush(lp)
        stop
      endif
      read( ni,*) yrflag,cvarin
      if(kn.eq.1) write(lp,*) cvarin,' = ',yrflag
      if (cvarin.ne.'yrflag') then
        write(lp,*)
        write(lp,*) 'error in getdat - input ',cvarin,
     &                        ' but should be yrflag'
        write(lp,*)
        call flush(lp)
        stop
      endif
      read( ni,*) idmtst,cvarin
      if(kn.eq.1) write(lp,*) cvarin,' = ',idmtst
      if (cvarin.ne.'idm   ') then
        write(lp,*)
        write(lp,*) 'error in getdat - input ',cvarin,
     &                        ' but should be idm   '
        write(lp,*)
        call flush(lp)
        stop
      endif
      read( ni,*) jdmtst,cvarin
      if(kn.eq.1) write(lp,*) cvarin,' = ',jdmtst
      if (cvarin.ne.'jdm   ') then
        write(lp,*)
        write(lp,*) 'error in getdat - input ',cvarin,
     &                        ' but should be jdm   '
        write(lp,*)
        call flush(lp)
        stop
      endif
c
      if (idmtst.ne.idm .or. jdmtst.ne.jdm) then
        write(lp,*)
        write(lp,*) 'error in getdat - input idm,jdm',
     &                        ' not consistent with parameters'
        write(lp,*) 'idm,jdm = ',idm,   jdm,   '  (REGION.h)'
        write(lp,*) 'idm,jdm = ',idmtst,jdmtst,'  (input)'
        write(lp,*)
        call flush(lp)
        stop
      endif
c
c --- artype
c
      read( ni,'(a)') cline
      if(kn.eq.1) write(lp,'(a)') trim(cline)
      if     (cline(24:28).eq.'model') then
        artype = 1  ! standard archive
      elseif (cline(24:28).eq.' mean') then
        artype = 2  ! mean archive, see below
      else
        write(lp,*)
        write(lp,*) 'error in getdat - wrong archive type.'
        write(lp,*) 'only "model", or " mean" allowed.'
        write(lp,*)
        call flush(lp)
        stop
      endif 
      call flush(lp)
c
c --- detect version 2.2 normal archive files
c
      read (ni,'(a)',end=6) cline
      if(kn.eq.1) write(lp,'(a)')  trim(cline)
      call flush(lp)
      i = index(cline,'=')
      read (cline(i+1:),*)  nstep,time,layer,thet,hminb,hmaxb
      nodens = layer.ne.0
      call zaiosk(ni)
c
c --- skip surface fields, except ubaro and vbaro for k=1.
c
      do irec= 2,9999
        read (ni,'(a)',end=6) cline
        if     (kn.eq.1) then
          write(lp,'(a)')  trim(cline)
        endif
        if     (cline(1:8).eq.'u_btrop') then
          if     (kn.eq.1) then
            call flush(lp)
            i = index(cline,'=')
            read (cline(i+1:),*)  nstep,time,layer,thet,hminb,hmaxb
            call getfld(ubaro,    ni, hminb,hmaxb)
            write(lp,'("input  ",a," into ",a)') cline(1:8),'ubaro'
c
            read (ni,'(a)',end=6) cline
            write(lp,'(a)')  trim(cline)
            call flush(lp)
            i = index(cline,'=')
            read (cline(i+1:),*)  nstep,time,layer,thet,hminb,hmaxb
            call getfld(vbaro,    ni, hminb,hmaxb)
            write(lp,'("input  ",a," into ",a)') cline(1:8),'vbaro'
            if     (artype.eq.2) then  !kebaro
              read (ni,'(a)',end=6) cline
              write(lp,'(a)')  trim(cline)
              call zaiosk(ni)
            endif
          else !k>1
            call zaiosk(ni) !ubaro
            read (ni,'(a)',end=6) cline
            call zaiosk(ni) !vbaro
            if     (artype.eq.2) then  !kebaro
              read (ni,'(a)',end=6) cline
              call zaiosk(ni)
            endif
          endif
c
          exit
        else
          call zaiosk(ni)
        endif
      enddo
c
      do k=1,kn
        if     (k.eq.2) then
c ---     already input at end of k=1 loop.
        else
          read (ni,'(a)',end=6)   cline
        endif
        if     (k.eq.kn) then
          write(lp,'(a)')    trim(cline)
          call flush(lp)
        endif
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time,layer,thet,hminb,hmaxb
        if     (cline(1:8).ne.'u-vel.  ') then
          write(lp,*)
          write(lp,*) 'error in getdat - layer ',k,
     &               ' does not exist (kk= ',kk,')'
          write(lp,*)
          call flush(lp)
          stop
        endif
        if     (k.eq.kn) then
          call getfld(u(1,1,1), ni, hminb,hmaxb)  !single layer input
          if     (k.eq.1 .or. k.eq.kk) then
            write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'u ',k
            call flush(lp)
          endif
        else
          call zaiosk(ni)
        endif
c
        read (ni,'(a)',end=6) cline
        if     (k.eq.kn) then
          write(lp,'(a)')  trim(cline)
          call flush(lp)
          i = index(cline,'=')
          read (cline(i+1:),*)  nstep,time,layer,thet,hminb,hmaxb
          call getfld(v(1,1,1), ni, hminb,hmaxb)  !single layer input
          if     (k.eq.1 .or. k.eq.kk) then
            write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'v ',k
            call flush(lp)
          endif
        else
          call zaiosk(ni)
        endif
c
        if     (artype.eq.2) then
c
c ---     mean archive, skip ke and convert from total to baroclinic velocity
c
          read (ni,'(a)',end=6) cline
          call zaiosk(ni)
c
          if     (k.eq.kn) then
            do j= 1,jj
              do i= 1,ii
                 u(i,j,1) = u(i,j,1) - ubaro(i,j)
                 v(i,j,1) = v(i,j,1) - vbaro(i,j)
              enddo
            enddo
          endif
        endif !artype
c
        read (ni,'(a)',end=6) cline
        if     (k.eq.kn) then
          write(lp,'(a)')  trim(cline)
          call flush(lp)
          i = index(cline,'=')
          read (cline(i+1:),*)  nstep,time,layer,thet,hminb,hmaxb
          call getfld(dp(1,1,1), ni, hminb,hmaxb)  !single layer input
          if     (k.eq.1 .or. k.eq.kk) then
            write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'dp',k
            call flush(lp)
          endif
        else
          call zaiosk(ni)
        endif
c
        if     (nodens) then
          do irec= 1,2 !temp,saln
            read (ni,'(a)',end=6) cline
            call zaiosk(ni)
          enddo
        else
          do irec= 1,3 !temp,saln,th3d
            read (ni,'(a)',end=6) cline
            call zaiosk(ni)
          enddo
        endif !nodens:else
        if     (kk.ne.1) then
          if     (k.eq.1) then
            do ktr= 1,999
              read (ni,'(a)',end=6) cline
              if     (cline(1:8).ne.'tracer  ' .and.
     &                cline(1:8).ne.'viscty  ' .and.
     &                cline(1:8).ne.'t-diff  ' .and.
     &                cline(1:8).ne.'s-diff  '      ) then
                exit !end of tracers and visc/diff
              else
                call zaiosk(ni)
              endif
            enddo !ktr
            ntr=ktr-1
          else
            do ktr= 1,ntr
              read (ni,'(a)',end=6) cline
              call zaiosk(ni)
            enddo
          endif  !k.eq.1:else
        endif  !kk.ne.1
        theta(k)=thet
      enddo  !k=1,kk
c
      close( unit=14)
      call zaiocl(14)
      return
c
c --- unexpected end of file
 6    continue
        write (lp,*) '***** unexpected end of archive file *****'
        write (lp,*) flnm(1:l-2)//'.b'
        call flush(lp)
        stop '(e-o-f)'
      end

      subroutine getdat_mean(flnm,timef,timel, iexpt,yrflag)
      use mod_trans  ! HYCOM transport section archive array interface
      use mod_za     ! HYCOM array I/O interface
c
      implicit none
c
      character        flnm*(*)
      double precision timef,timel
      integer          iexpt,yrflag
c
c --- read mean model fields for transport sections.
c --- HYCOM 2.0/1 array I/O mean archive file.
c
      character        cline*80
      character        preambl(5)*79
      character        cvarin*6
      real             hminb,hmaxb,thet
      double precision time
      integer          i,idmtst,irec,iversn,ios,
     &                 j,jdmtst,k,l,layer,ni,nstep,
     &                 ktr,ntr
      logical          icegln
c
      data ni/14/
c
      l = len_trim(flnm)
      open (unit=ni,file=flnm(1:l-2)//'.b',form='formatted',
     &      status='old',action='read')
      call zaiopf(flnm(1:l-2)//'.a','old', ni)
c
      read( ni,'(a)') cline
      write(lp,'(a)') trim(cline)
      read( ni,'(a)') cline
      write(lp,'(a)') trim(cline)
      read( ni,'(a)') cline
      write(lp,'(a)') trim(cline)
      read( ni,'(a)') cline
      write(lp,'(a)') trim(cline)
      call flush(lp)
      read( ni,*) iversn,cvarin
      write(lp,*) cvarin,' = ',iversn
      if (cvarin.ne.'iversn') then
        write(lp,*)
        write(lp,*) 'error in getdat - input ',cvarin,
     &                          ' but should be iversn'
        write(lp,*)
        call flush(lp)
        stop
      endif
      read( ni,*) iexpt,cvarin
      write(lp,*) cvarin,' = ',iexpt
      if (cvarin.ne.'iexpt ') then
        write(lp,*)
        write(lp,*) 'error in getdat - input ',cvarin,
     &                        ' but should be iexpt '
        write(lp,*)
        call flush(lp)
        stop
      endif
      read( ni,*) yrflag,cvarin
      write(lp,*) cvarin,' = ',yrflag
      if (cvarin.ne.'yrflag') then
        write(lp,*)
        write(lp,*) 'error in getdat - input ',cvarin,
     &                        ' but should be yrflag'
        write(lp,*)
        call flush(lp)
        stop
      endif
      read( ni,*) idmtst,cvarin
      write(lp,*) cvarin,' = ',idmtst
      if (cvarin.ne.'idm   ') then
        write(lp,*)
        write(lp,*) 'error in getdat - input ',cvarin,
     &                        ' but should be idm   '
        write(lp,*)
        call flush(lp)
        stop
      endif
      read( ni,*) jdmtst,cvarin
      write(lp,*) cvarin,' = ',jdmtst
      if (cvarin.ne.'jdm   ') then
        write(lp,*)
        write(lp,*) 'error in getdat - input ',cvarin,
     &                        ' but should be jdm   '
        write(lp,*)
        call flush(lp)
        stop
      endif
c
      if (idmtst.ne.idm .or. jdmtst.ne.jdm) then
        write(lp,*)
        write(lp,*) 'error in getdat - input idm,jdm',
     &                        ' not consistent with parameters'
        write(lp,*) 'idm,jdm = ',idm,   jdm,   '  (REGION.h)'
        write(lp,*) 'idm,jdm = ',idmtst,jdmtst,'  (input)'
        write(lp,*)
        call flush(lp)
        stop
      endif
c
c --- artype==2 for mean archive files
c
      read( ni,'(a)') cline
      write(lp,'(a)') trim(cline)
      if     (cline(25:28).ne.'mean') then
        write(lp,*)
        write(lp,*) 'error in getdat - ',
     &              'archive is not a mean file'
        write(lp,*)
        call flush(lp)
        stop
      endif
c
c --- skip surface fields, but extract timef and timel.
c
      read (ni,'(a)',end=6) cline
      call zaiosk(ni)
      write(lp,'(a)')  trim(cline)
      i = index(cline,'=')
      read (cline(i+1:),*)  nstep,timef
c
      read (ni,'(a)',end=6) cline
      call zaiosk(ni)
      write(lp,'(a)')  trim(cline)
      i = index(cline,'=')
      read (cline(i+1:),*)  nstep,timel
c
c --- only need ubaro and vbaro.
c
      do irec= 1,9999
        read (ni,'(a)',end=6) cline
        write(lp,'(a)')  trim(cline)
        if     (cline(1:8).eq.'u_btrop') then
          i = index(cline,'=')
          read (cline(i+1:),*)  nstep,time,layer,thet,hminb,hmaxb
          call getfld(ubaro,    ni, hminb,hmaxb)
          write(lp,'("input  ",a," into ",a)') cline(1:8),'ubaro'
c
          read (ni,'(a)',end=6) cline
          write(lp,'(a)')  trim(cline)
          i = index(cline,'=')
          read (cline(i+1:),*)  nstep,time,layer,thet,hminb,hmaxb
          call getfld(vbaro,    ni, hminb,hmaxb)
          write(lp,'("input  ",a," into ",a)') cline(1:8),'vbaro'
          exit
        else
          call zaiosk(ni)
        endif
      enddo
      read (ni,'(a)',end=6) cline
      write(lp,'(a)')  trim(cline)
      call zaiosk(ni)
      call flush(lp)
c
      do k=1,kk
        if     (k.eq.2) then
c ---     already input at end of k=1 loop.
        else
          read (ni,'(a)',iostat=ios)   cline
        endif
        if     (ios.ne.0) then
          if     (k.ne.2) then
            goto 6
          endif
c
          write (lp,*) 'surface archive (barotropic transport only)'
          call flush(lp)
          do j= 1,jj
            do i= 1,ii
               u(i,j,1:kk-1) = 0.0
               v(i,j,1:kk-1) = 0.0
              dp(i,j,1:kk-1) = 0.0
               u(i,j,1:kk)   = ubaro(i,j)
               v(i,j,1:kk)   = vbaro(i,j)
              dp(i,j,  kk)   = depths(i,j)*9806.0
            enddo
          enddo
          exit
        endif
        if     (cline(1:8).ne.'u-vel.  ') then
          write(lp,'(a)') trim(cline)
          write(lp,*)
          write(lp,*) 'error in getdat - layer ',k,
     &               ' does not exist (kk= ',kk,')'
          write(lp,*)
          call flush(lp)
          stop
        endif
        write(lp,'(a)')         trim(cline)
        call flush(lp)
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time,layer,thet,hminb,hmaxb
        call getfld(u(1,1,k), ni, hminb,hmaxb)
        if     (k.eq.1 .or. k.eq.kk) then
          write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'u ',k
        endif
        call flush(lp)
c
        read (ni,'(a)',end=6) cline
        write(lp,'(a)')  trim(cline)
        call flush(lp)
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time,layer,thet,hminb,hmaxb
        call getfld(v(1,1,k), ni, hminb,hmaxb)
        if     (k.eq.1 .or. k.eq.kk) then
          write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'v ',k
        endif
        call flush(lp)
c
        read (ni,'(a)',end=6) cline
        call zaiosk(ni)
c
        read (ni,'(a)',end=6) cline
        write(lp,'(a)')  trim(cline)
        call flush(lp)
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time,layer,thet,hminb,hmaxb
        call getfld(dp(1,1,k), ni, hminb,hmaxb)
        if     (k.eq.1 .or. k.eq.kk) then
          write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'dp',k
        endif
        call flush(lp)
c
        do irec= 1,3
          read (ni,'(a)',end=6) cline
          call zaiosk(ni)
        enddo
        if     (kk.ne.1) then
          if     (k.eq.1) then
            do ktr= 1,999
              read (ni,'(a)',iostat=ios)   cline
              if     (ios.ne.0) then
                exit !end of input, but loop once more
              elseif (cline(1:8).ne.'tracer  ' .and.
     &                cline(1:8).ne.'viscty  ' .and.
     &                cline(1:8).ne.'t-diff  ' .and.
     &                cline(1:8).ne.'s-diff  '      ) then
                exit !end of tracers and visc/diff
              else
                call zaiosk(ni)
              endif
            enddo !ktr
            ntr=ktr-1
          else
            do ktr= 1,ntr
              read (ni,'(a)',end=6) cline
              call zaiosk(ni)
            enddo
          endif  !k.eq.1:else
        endif  !kk.ne.1
        theta(k)=thet
      enddo  !k=1,kk
      call flush(lp)
c
      close( unit=14)
      call zaiocl(14)
      return
c
c --- unexpected end of file
 6    continue
        write (lp,*) '***** unexpected end of archive file *****'
        write (lp,*) flnm(1:l-2)//'.b'
        call flush(lp)
        stop '(e-o-f)'
      end

      subroutine getfld(field, iunit, hminb,hmaxb)
      use mod_za ! HYCOM array I/O interface
c
      implicit none
c
      integer iunit
      real    field(idm,jdm), hminb,hmaxb
c
c --- read a single array
c
      integer i,j,mask(idm,jdm)
      real    hmina,hmaxa
c
      call zaiord(field,mask,.false., hmina,hmaxa, iunit)
c
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b files not consistent:',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        call flush(lp)
        stop
      endif
c
      do j= 1,jdm
        do i= 1,idm
          if     (field(i,j).gt.2.0**99) then
            field(i,j) = 0.0
          endif
        enddo
      enddo
      return
      end
