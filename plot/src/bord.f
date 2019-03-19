      subroutine bordplt(lalogr)
      use mod_plot  ! HYCOM plot array interface
      implicit none
c
c     line contour land/sea boundary
c
c     lalogr is the plotted lat/lon grid spacing
c       >0;    plot grid over land only
c       =0; no plot grid
c       <0;    plot grid over land and sea
c
c     if abs(lalogr) > 1000, it is +/- (1000 + grid*100).
c
c     lalogr is only significant on the first call,
c      i.e. land contouring is identical on all horrizonal plots.
c
      integer      lalogr
c
      integer         ipalet,nbase,ibase
      common /colopt/ ipalet,nbase,ibase(2,99)
c
      integer                    :: i,ikeep,j,nspval
      real                       :: ploni,plonim1,ll360
      real,    save, allocatable :: landlat(:,:),landlon(:,:),
     &                               mcoast(:,:)
      real,    save              :: qlatlo,qlathi,qlonlo,qlonhi,
     &                              qlatlon
      logical, save              :: lmcoast = .false.
      logical, save              :: lfirst  = .true.
      real,    parameter         :: spval   = -.03125
c
      if     (lfirst) then
c ---   prepare for lat/lon grid over land
        allocate( landlat(ii,jj), landlon(ii,jj), mcoast(ii,jj) )
        if     (abs(lalogr).lt.1000) then
          qlatlon = abs(lalogr)
        else
          qlatlon = (abs(lalogr) - 1000)/100.0
        endif
        qlatlo =  999.0
        qlathi = -999.0
        qlonlo =  999.0
        qlonhi = -999.0
        do j= 1,jj
          plonim1 = plon(1,j)
          nspval  = 0
*         write(6,*) 'j = ',j
          do i= 1,ii
            ploni = plon(i,j)
            if     (ploni-plonim1.ge. 180.0+qlatlon) then
              ploni = ploni - 360.0
            elseif (ploni-plonim1.le.-180.0+qlatlon) then
              ploni = ploni + 360.0
            endif
            plonim1 = ploni
            if     (lalogr.lt.0 .or. coast(i,j).eq.0.0) then
              landlat(i,j) = plat(i,j)
              qlatlo  = min( landlat(i,j), qlatlo )
              qlathi  = max( landlat(i,j), qlathi )
              if     (j.eq.1) then
                landlon(i,j) = ploni
                qlonlo  = min( landlon(i,j), qlonlo )
                qlonhi  = max( landlon(i,j), qlonhi )
              elseif (landlon(i,j-1).eq.spval .or.
     &                abs(ploni-landlon(i,j-1)).lt.350.0) then
                landlon(i,j) = ploni
                qlonlo  = min( landlon(i,j), qlonlo )
                qlonhi  = max( landlon(i,j), qlonhi )
              elseif (nspval.gt.0) then
                if     (ploni-landlon(i,j-1).ge.350.0) then
                  landlon(i,j) = ploni - 360.0
*                 write(6,'(a,2f10.2)') 
*    &              'bordplt - SEAM llon,plon= ',landlon(i,j),ploni
                else
                  landlon(i,j) = ploni + 360.0
*                 write(6,'(a,2f10.2)') 
*    &              'bordplt - SEAM llon,plon= ',landlon(i,j),ploni
                endif
                qlonlo  = min( landlon(i,j), qlonlo )
                qlonhi  = max( landlon(i,j), qlonhi )
              else
                landlon(i,j) = spval  ! seam in longitude
                nspval = nspval + 1
*               write(6,'(a,2f10.2)') 
*    &            'bordplt - seam llon,plon= ',landlon(i,j-1),ploni
              endif
              if     (landlon(i,j).ne.spval .and.
     &                abs(plat(i,j)).gt.90.0-qlatlon) then
                ll360 = mod(abs(landlon(i,j)),360.0)
                if     (    ll360       .ge.qlatlon .and.
     &                  abs(ll360-180.0).ge.qlatlon .and.
     &                  abs(ll360-360.0).ge.qlatlon      ) then
*                 write(6,'(a,2f10.2)') 
*    &              'bordplt - near pole llon =',landlon(i,j),plat(i,j)
                  landlon(i,j) = spval  ! near pole and not near date line
                endif  !not near date line
              endif  !near pole
            else  !no lon,lat contour
              landlat(i,j) = spval
              landlon(i,j) = spval
            endif !lalogr or coast:else
            if     (depths(i,j).gt.0.0) then
              mcoast(i,j) = 1.0  ! model-sea
            else
              mcoast(i,j) = 0.0  ! model-land
            endif
            lmcoast = lmcoast .or. mcoast(i,j).ne.coast(i,j)
          enddo !i
        enddo !j
        if     (lalogr.eq.0) then
          write(6,'(a)')
     &      'bordplt - no grid over land'
        else
          qlatlo  = nint(qlatlo/qlatlon)*qlatlon - qlatlon
          qlathi  = nint(qlathi/qlatlon)*qlatlon + qlatlon
          qlonlo  = nint(qlonlo/qlatlon)*qlatlon - qlatlon
          qlonhi  = nint(qlonhi/qlatlon)*qlatlon + qlatlon
          write(6,'(a,3f10.2)') 
     &      'bordplt - qlatlo,qlathi,qlatlon = ',qlatlo,qlathi,qlatlon
          write(6,'(a,3f10.2)') 
     &      'bordplt - qlonlo,qlonhi,qlatlon = ',qlonlo,qlonhi,qlatlon
        endif
        lfirst  = .false.
      endif
c
c     color-fill land and contour land/sea boundary
c
      ikeep  = ipalet
      if     (lmcoast) then
        ipalet = 101  !outside conventional palette range
        call conrec(mcoast,ii,ii1,jj1,-0.5,0.5,1.0,1,-1,0)
      endif
      ipalet = 100  !outside conventional palette range
      call conrec(coast,ii,ii1,jj1,-0.5,0.5,1.0,1,-1,0)
      ipalet = -9
      call conrec(coast,ii,ii1,jj1,-0.5,1.5,1.0,1,-1,0)
      ipalet = ikeep
c
c     mark lat/lon lines over land
c
      if     (qlatlon.ne.0.0) then
        ikeep  = ipalet
        ipalet = -9
        call conrec(landlat,ii,ii1,jj1,qlatlo,qlathi,qlatlon,1,-1,0)
        call conrec(landlon,ii,ii1,jj1,qlonlo,qlonhi,qlatlon,1,-1,0)
        ipalet = ikeep
      endif
      return
      end 
