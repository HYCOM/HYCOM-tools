      program topo_merge
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      character*79    ctitle(5)
      integer         if(999),il(999),jf(999),jl(999)
      real*4          scale(999),boxscl(4,999)
      namelist/merge/ ctitle,if,il,jf,jl,scale,boxscl
c
c --- merge two hycom topography files.
c 
      integer   i,im2,ip2,j,ja,jm2,jp2,k,nfill,nzero
      real      hmaxa,hmaxb,hmina,hminb,dhsave,dhmin,dhmax,si,sj
      character preambl(5)*79,cline*80
c
      integer, allocatable :: ip(:,:)
      real,    allocatable :: dh(:,:),dh2(:,:),s2(:,:)
c
      call xcspmd  !input idm,jdm
      allocate( ip( idm,jdm) )
      allocate( dh( idm,jdm) )
      allocate( dh2(idm,jdm) )
      allocate( s2( idm,jdm) )
c
c ---      ctitle      - HYCOM 5-line title for merged bathymetry
c ---      if,il,jf,jl - array box where scale is applied
c ---                     the allowed range for if,il is 1 to idm
c ---                     the allowed range for jf,jl is 1 to jdm-1
c ---                     = 4*0; end of box list
c ---      scale       - multiplier for 2nd bathymetry in box
c ---                     =-9.0; use boxscl linearly varying multiplier
c ---                     < 0.0; use 1st bathymetry, but merge land
c ---                     = 0.0; use 1st bathymetry only (default)
c ---                     = 1.0; use 2nd bathymetry only
c ---                     = 0.0-1.0; use fraction of each bathymetry
c ---                     = 2.0; use 2nd bathymetry where 1st is land
c ---                     = 3.0; use 2nd bathymetry where 1st is near land
c ---      boxscl      - alternative multiplier for 2nd bathymetry in box
c ---      boxscl(1,:) - scale factor for (if,jf) between 0.0 and 1.0
c ---      boxscl(2,:) - scale factor for (il,jf) between 0.0 and 1.0
c ---      boxscl(3,:) - scale factor for (il,jl) between 0.0 and 1.0
c ---      boxscl(4,:) - scale factor for (if,jl) between 0.0 and 1.0
c
      ctitle(:) = ' '
      if(      :) = 0
      il(      :) = 0
      jf(      :) = 0
      jl(      :) = 0
      scale(   :) = 0.0
      boxscl(:,:) = 0.0
      read( 5,merge)
      write(6,merge)
c
      call zhopen(51, 'formatted', 'old', 0)
      read (51,'(a79)') preambl
      read (51,'(a)')   cline
      close(unit=51)
      write(6,'(/a/(a))') '1st header:',
     &                    preambl,cline(1:len_trim(cline))
c
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
c
      call zaiost
      call zaiopn('old', 51)
      call zaiord(dh,ip,.false., hmina,hmaxa, 51)
      call zaiocl(51)
c
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        write(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - 1st .a and .b topography files not consistent:',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        call zhflsh(6)
        stop
      endif
c
      call zhopen(61, 'formatted', 'old', 0)
      read (61,'(a79)') preambl
      read (61,'(a)')   cline
      close(unit=61)
      write(6,'(/a/(a))') '2nd header:',
     &                     preambl,cline(1:len_trim(cline))
c
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
c
      call zaiopn('old', 61)
      call zaiord(dh2,ip,.false., hmina,hmaxa, 61)
      call zaiocl(61)
c
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        write(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - 2nd .a and .b topography files not consistent:',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        call zhflsh(6)
        stop
      endif
c
c     initialize s2 from /merge/
c
      s2(:,:) = 0.0
      do k= 1,999
        if     (min(if(k),il(k),jf(k),jl(k)).gt.0   .and.
     +                              scale(k).ne.0.0      ) then
          if     (if(k).gt.il(k) .or.
     +            if(k).lt.1     .or.
     +            il(k).gt.idm       ) then
            write(6,9000) k,if(k),il(k),idm
            stop
          endif
          if     (jf(k).gt.jl(k) .or.
     +            jf(k).lt.1     .or.
     +            jl(k).gt.jdm-1     ) then
            write(6,9100) k,jf(k),jl(k),jdm-1
            stop
          endif
C
          if     (scale(k).gt.0.0) then
            do j= jf(k),jl(k)
              do i= if(k),il(k)
                s2(i,j) = max( s2(i,j), scale(k) )
              enddo
            enddo
            write(6,*) 'scale> 0: ',
     &                 (if(k)+il(k))/2,(jf(k)+jl(k))/2,
     &              s2((if(k)+il(k))/2,(jf(k)+jl(k))/2)
          elseif (scale(k).eq.-9.0) then  !use boxscl
            do j= jf(k),jl(k)
              sj = real(j-jf(k))/real(jl(k)-jf(k))  !0.0 to 1.0
              do i= if(k),il(k)
                si = real(i-if(k))/real(il(k)-if(k))  !0.0 to 1.0
                s2(i,j) = max( s2(i,j),
     &                         boxscl(1,k) +
     &                    (1.0-si)*(1.0-sj)*(boxscl(1,k)-boxscl(1,k)) +
     &                         si *(1.0-sj)*(boxscl(2,k)-boxscl(1,k)) +
     &                         si *     sj *(boxscl(3,k)-boxscl(1,k)) +
     &                    (1.0-si)*     sj *(boxscl(4,k)-boxscl(1,k))  )
              enddo
            enddo
            write(6,*) 'scale=-9: ',if(k),jf(k),s2(if(k),jf(k))
            write(6,*) 'scale=-9: ',if(k),jl(k),s2(if(k),jl(k))
            write(6,*) 'scale=-9: ',il(k),jf(k),s2(il(k),jf(k))
            write(6,*) 'scale=-9: ',il(k),jl(k),s2(il(k),jl(k))
            write(6,*) 'scale=-9: ',
     &                 (if(k)+il(k))/2,(jf(k)+jl(k))/2,
     &              s2((if(k)+il(k))/2,(jf(k)+jl(k))/2)
          else
            do j= jf(k),jl(k)
              do i= if(k),il(k)
                s2(i,j) = -1.0  !merge land masks
              enddo
            enddo
            write(6,*) 'scale= X: ',
     &                 (if(k)+il(k))/2,(jf(k)+jl(k))/2,
     &              s2((if(k)+il(k))/2,(jf(k)+jl(k))/2)
          endif
        endif
      enddo
c
c     merge the bathymetries, result in dh.
c
      do j= 1,jdm
        jm2 = max(j-2,   1)
        jp2 = min(j+2, jdm)
        do i= 1,idm
          if (dh(i,j).gt.2.0**99) then
            dh(i,j) = 0.0
          endif
          if (dh2(i,j).gt.2.0**99) then
            dh2(i,j) = 0.0
          endif
          dhsave = dh(i,j)
          if     (s2(i,j).eq.0.0) then  ! 1st bathymetry only
            if (dh(i,j).le.0.0) then
              ip(i,j) = 0
              dh(i,j) = 0.0
              dhsave  = 0.0
            else
              ip(i,j) = 1
            endif
          elseif (s2(i,j).eq.1.0) then  ! 2nd bathymetry only
            if (dh2(i,j).le.0.0) then
              ip(i,j) = 0
              dh(i,j) = 0.0
            else
              ip(i,j) = 1
              dh(i,j) = dh2(i,j)
            endif
          elseif (s2(i,j).eq.2.0) then  ! 2nd bathymetry where 1st is land
            if (dh(i,j).le.0.0) then
              if (dh2(i,j).le.0.0) then
                ip(i,j) = 0
                dh(i,j) = 0.0
                dhsave  = 0.0
              else
                ip(i,j) = 1
                dh(i,j) = dh2(i,j)
              endif
            else
              ip(i,j) = 1
            endif
          elseif (s2(i,j).eq.3.0) then  ! 2nd bathymetry where 1st is near land
            im2   = max(i-2,   1)
            ip2   = min(i+2, idm)
            dhmin = minval(dh(im2:ip2,jm2:jp2))
            dhmax = maxval(dh(im2:ip2,jm2:jp2))
            if (dhmin.le.0.0     .or.
     &          dhmax.gt.2.0**99     ) then
              if (dh2(i,j).le.0.0) then
                ip(i,j) = 0
                dh(i,j) = 0.0
                dhsave  = 0.0
              else
                ip(i,j) = 1
                dh(i,j) = dh2(i,j)
              endif
            else
              ip(i,j) = 1
            endif
          elseif (s2(i,j).lt.0.0) then  ! merge land masks only
            if (min(dh(i,j),dh2(i,j)).le.0.0) then
              ip(i,j) = 0
              dh(i,j) = 0.0
            else
              ip(i,j) = 1
            endif
          else  ! merge bathymetries
            if (min(dh(i,j),dh2(i,j)).le.0.0) then
              ip(i,j) = 0
              dh(i,j) = 0.0
            else
              ip(i,j) = 1
              dh(i,j) = (1.0-s2(i,j))*dh(i,j) + s2(i,j)*dh2(i,j)
            endif
          endif
          if     (i.eq.idm/2 .or. j.eq.jdm/2) then
            if     (s2(i,j).gt.0.0 .and. s2(i,j).lt.1.0) then
              write(6,'(a,2i6,4f10.3)') 'I,J,dh,s,dh1,dh2 =',
     &          i,j,dh(i,j),s2(i,j),dhsave,dh2(i,j)
            elseif (mod(i,100).eq.1 .or. mod(j,100).eq.1) then
              write(6,'(a,2i6,4f10.3)') 'i,j,dh,s,dh1,dh2 =',
     &          i,j,dh(i,j),s2(i,j),dhsave,dh2(i,j)
            endif
          endif
        enddo
      enddo
c
c --- warn about single-width inlets and 1-point seas.
c
      write (6,'(/a)') 'single-width inlet check:'
      do j=1,jdm-1  !safe for closed and arctic
        do i=1,idm
          nzero=0
          if (dh(i,j).gt.0.0) then
            if     (i.eq.  1) then
              if (dh(idm,j).le.0.0) nzero=nzero+1  !assuming periodic is safe
              if (dh(i+1,j).le.0.0) nzero=nzero+1
            elseif (i.eq.idm) then
              if (dh(i-1,j).le.0.0) nzero=nzero+1
              if (dh(  1,j).le.0.0) nzero=nzero+1  !assuming periodic is safe
            else
              if (dh(i-1,j).le.0.0) nzero=nzero+1
              if (dh(i+1,j).le.0.0) nzero=nzero+1
            endif
            if (j.eq.  1.or.dh(i,j-1).le.0.0) nzero=nzero+1
            if (j.eq.jdm.or.dh(i,j+1).le.0.0) nzero=nzero+1
            if (nzero.ge.3) then
              write (6,'(a,i6,a,i6,a,i2,a)')
     +          ' dh(',i,',',j,') has',
     +          nzero,' land nieghbours'
              nfill=nfill+1
            end if
          end if
        enddo
      enddo
      if (nfill.gt.0) then
        write (6,'(/a/)')
     &   'WARNING - single-width inlets and/or 1-point seas exist'
      endif
c
c --- write out the merged hycom topography file.
c
      call zaiopn('new', 71)
      call zaiowr(dh, ip,.true., hmina,hmaxa, 71, .false.)
c
      call zhopen(71, 'formatted', 'new', 0)
      write(71,'(a79)') ctitle
      write(71,7100) hmina,hmaxa
c
      write(6, *)
      write(6, '(a79)') ctitle
      write(6, 7100) hmina,hmaxa
      write(6, *)
c
 7100 format('min,max depth =',2f12.5)
 9000 format('error - illegal if or il for k =',i3,
     +       '   must have 1<=if(k)<=il(k)<=idm' /
     +       'if(k),il(k),idm = ',3i6 /)
 9100 format('error - illegal jf or jl for k =',i3,
     +       '   must have 1<=jf(k)<=jl(k)<=jdm-1' /
     +       'jf(k),jl(k),jdm-1 = ',3i6 /)
      end
