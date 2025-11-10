      program grid_partition_mom6
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
c --- maximum number of tiles along one dimension
      integer    maxpe1
      parameter (maxpe1=512)
c
c --- halo size
      integer    mbdy
      parameter (mbdy=4)
c
      integer   iipe(maxpe1,maxpe1),ispt(maxpe1,maxpe1)
      integer   iipx(maxpe1,maxpe1),ispx(maxpe1,maxpe1)
      integer   jjpe(maxpe1),jspt(maxpe1)  ! always separable
      integer   jjnp(maxpe1)
      integer   ibig,jbig,minsea,maxsea,mavsea,nsea,nseao,nsea2,nreg
c
      real      hmaxa,hmaxb,hmina,hminb
      character preambl(5)*79,cline*80
c
      integer, allocatable :: map(:,:),ip(:,:),mapsum(:)
      real,    allocatable :: depths(:,:)
c
c --- This program reads in a standard HYCOM depth file and writes out 
c --- the patch distribution for a given npe x mpe 2-d node grid.
c --- This version for MOM6 layouts, i.e. constant-sized tiles.
c
      logical ldiscard,luniform
      integer i,ic,iipxnj,j,jpe_lim,n,nl,npe,mpe,mpe_old,nmpe
      integer i_w,i_e,j_s,j_n,ifrst,ilast,isec,nchar
      real    sfudge
c
      character char3*3
      character fmt*13
      data fmt/'(i4,1x,120i1)'/
c
      call xcspmd  !input idm,jdm
      allocate( map(1-mbdy:idm+mbdy,1-mbdy:jdm+mbdy) )
      allocate( ip(      0:idm+1,        0:jdm+1)    )
      allocate( mapsum(1-mbdy:idm+jdm+mbdy) )
      allocate( depths(idm,jdm) )
c
c --- read in 2-d node grid dimensions.
c ---   npe,mpe: 2-d node grid (-ve for constant-sized tiles)
c ---   sfudge:  size fudge factor (0.5 to 9.9, larger for more variation)
c ---              < 1.0 to keep all  constant-sized tiles
c
      read(5,*) npe,mpe,sfudge
c
      luniform = .true.
      ldiscard = sfudge.ge.1.0
c
      if     (min(npe,mpe).lt.1 .or. max(npe,mpe).gt.maxpe1) then
        write(6,'(/ a,i3 /)') 'error - npe,mpe must be between 1 and ',
     &                        maxpe1
        call zhflsh(6)
        stop
      endif
c
      if     (sfudge.lt.0.5 .or. sfudge.gt.9.9) then
        write(6,'(/ a /)')
     &     'error - sfudge must be between 0.5 and 9.9'
        call zhflsh(6)
        stop
      endif
c
c --- acquire basin depths from unit 51.
c
      call zhopen(51, 'formatted', 'old', 0)
      read (51,'(a79)') preambl
      read (51,'(a)')   cline
      close(unit=51)
      write(6,'(/(a))') preambl,cline(1:len_trim(cline))
c
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
c
      call zaiost
      call zaiopn('old', 51)
      call zaiord(depths,map,.false., hmina,hmaxa, 51)
      call zaiocl(51)
c
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        write(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b topography files not consistent:',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        call zhflsh(6)
        stop
      endif
c
c --- land mask
c
      do j= 1,jdm
        do i= 1,idm
          if     (depths(i,j).lt.2.0**99) then
            ip(i,j) = 1
          else
            ip(i,j) = 0
          endif
        enddo
      enddo
c
c --- handle all region types (nreg = closed:0; periodic:2; arctic:3)
c
      if     (sum(ip(1:idm-1,jdm)).ne.0) then
c
c ---   arctic dipole patch.
c
        if     (luniform) then
          if     (npe.ne.1 .and. mod(npe,2).ne.0) then
            write(6,'(/ a /)')
     &        'error - pan-am grid must have npe even'
            call zhflsh(6)
            stop
          endif
        else
          write(6,'(/ a /)')
     &      'error - pan-am grid must use equal-size tiles'
          call zhflsh(6)
          stop
        endif
        do i= 1,idm
          ic = idm+1-i
          ip(i,    0) = 0
          ip(i,jdm+1) = ip(ic,jdm-2)
        enddo
        nreg = 3  ! arctic
      else
        ip(:,    0) = 0
        ip(:,jdm+1) = 0
        nreg = 0  ! closed or periodic
      endif  !luniform
c
      if     (sum(ip(idm,:)).ne.0) then
c
c ---   periodic domain.
c
        ip(    0,:) = ip(idm,:)
        ip(idm+1,:) = ip(  1,:)
        nreg = max( nreg, 2 )  ! periodic or arctic
      else
c
c ---   treatable as a periodic domain.
c
        ip(    0,:) = ip(idm,:)
        ip(idm+1,:) = ip(  1,:)
        nreg = max( nreg, 2 )  ! periodic or arctic
      endif
c
c --- allow for all grids (p,q,u,v).
c
      do j= 1,jdm
        do i= 1,idm
          map(i,j) = ip(i,j)
          if (ip(i-1,j).gt.0.and.ip(i,j).gt.0) then
            map(i,j)=1  ! not needed, subset of ip
          endif
          if (ip(i,j-1).gt.0.and.ip(i,j).gt.0) then
            map(i,j)=1  ! not needed, subset of ip
          endif
          if (min(ip(i,j),ip(i-1,j),ip(i,j-1),ip(i-1,j-1)).gt.0) then
            map(i,j)=1  ! not needed, subset of ip
          elseif ((ip(i  ,j).gt.0.and.ip(i-1,j-1).gt.0).or.
     &            (ip(i-1,j).gt.0.and.ip(i  ,j-1).gt.0)    ) then
            map(i,j)=1  !     needed
          endif
          if (ip(i-1,j-1).gt.0.and.ip(i,j-1).gt.0) then
            map(i,j)=1  ! pvtrop for u
          endif
          if (ip(i-1,j-1).gt.0.and.ip(i-1,j).gt.0) then
            map(i,j)=1  ! pvtrop for v
          endif
        enddo
      enddo
c
c --- map halo.
c
*     if     (nreg.lt.3) then  ! not arctic
      if     (nreg.eq.0) then  ! assume arctic unless closed
        map(1:idm,1-mbdy:0)        = 0
        map(1:idm, jdm+1:jdm+mbdy) = 0
      else                     ! arctic
        do j= 1,mbdy
          do i= 1,idm
            ic = idm+1-i
            map(i,  1-j) = 0
            map(i,jdm+j) = map(ic,jdm-1-j)
          enddo
        enddo
      endif
      if     (nreg.eq.0) then  ! closed
        map(1-mbdy:0,       :) = 0
        map( idm+1:idm+mbdy,:) = 0
      else                     ! periodic or arctic
        map(1-mbdy:0,       :) = map( idm-mbdy+1:idm, :)
        map( idm+1:idm+mbdy,:) = map(          1:mbdy,:)
      endif
c
      if     (.FALSE.) then
c
c ---   write out -map-  array, with 0 or * for land and 1 or 2 for sea.
c ---   data are written in strips Nchar points wide
c
        do j= 1,jdm
          do i= 1,idm
            if     (ip(i,j).eq.0 .and. map(i,j).eq.1) then
              ip(i,j) = -1
            endif
          enddo
        enddo
        nchar=100
        isec =(idm-1)/nchar
        do ifrst=0,nchar*isec,nchar
          ilast=min(idm,ifrst+nchar)
          write (char3,'(i3)') ilast-ifrst
          fmt(8:10)=char3
          write(6,'(//''ip array, cols'',i6,'' --'',i6)') ifrst+1,ilast
          write(6,fmt) -9999,(mod(i,10),i=ifrst+1,ilast)
          write(6,*) 
          do j= jdm,1,-1
            write(6,fmt) j,(ip(i,j),i=ifrst+1,ilast)
          enddo
          write(6,*)
          write(6,fmt) -9999,(mod(i,10),i=ifrst+1,ilast)
        enddo
        write(6,*)
      endif  !write -map- 
c
      if     (luniform) then
c
c ---   uniform 2-d distribution, mom6 sizes.
c ---   discard tiles over land.
c
        iipe(1:npe,1:mpe) = idm/npe
        n = idm - iipe(1,1)*npe
        do i= 1,n
          iipe(i,1:mpe) = idm/npe + 1  ! 1st edge tiles larger
        enddo
        ispt(1,1:mpe) = 1
        do i= 2,npe
          ispt(i,1:mpe) = ispt(i-1,1:mpe) + iipe(i-1,1:mpe)
        enddo
c
        if     (jdm/mpe.lt.mbdy) then
          write(6,'(/ a,i3 /)') 'error - mpe is too large'
          call zhflsh(6)
          stop
        endif
c
        jjpe(1:mpe) = jdm/mpe
        n = jdm - jjpe(1)*mpe
        do j= 1,n
          jjpe(j) = jdm/mpe + 1  ! 1st edge tiles larger
        enddo
        if     (jjpe(mpe).le.mbdy) then
          if     (n.lt.1) then
            write(6,'(/ a,i3 /)') 'error - mpe is too large'
            call zhflsh(6)
            stop
          endif
          jjpe(mpe)   = jjpe(mpe)  +1
          jjpe(mpe/2) = jjpe(mpe/2)-1
        endif
        jspt(1) = 1
        do j= 2,mpe
          jspt(j) = jspt(j-1) + jjpe(j-1)
        enddo
c
c ---   optionally discard all-land tiles.
c
        nmpe = 0
        do j= 1,mpe
          do i= 1,npe
            i_w = ispt(i,j)
            i_e = ispt(i,j) + iipe(i,j)-1
            j_s = jspt(  j)
            j_n = jspt(  j) + jjpe(  j)-1
c ---       assume N-S exchange is before E-W exchange.
            nsea = sum(map(i_w:i_e,j_s:j_n)) +
     &             sum(map(i_w:i_w+mbdy-1,j_s-mbdy:j_s-1)) +
     &             sum(map(i_e-mbdy+1:i_e,j_s-mbdy:j_s-1)) +
     &             sum(map(i_w:i_w+mbdy-1,j_n+1:j_n+mbdy)) +
     &             sum(map(i_e-mbdy+1:i_e,j_n+1:j_n+mbdy))
            if     (ldiscard .and. nsea.eq.0) then
              ispx(i,j) = 0
              iipx(i,j) = 0
            else
              nmpe = nmpe + 1
              ispx(i,j) = ispt(i,j)
              iipx(i,j) = iipe(i,j)
*             write(6,'(a,5i6,i8)') 'pe,tile,sea = ',
*    &          nmpe,ispt(i,j),ispt(i,j)+iipe(i,j)-1,
*    &               jspt(  j),jspt(  j)+jjpe(  j)-1,nsea
            endif
          enddo !i
          if     (nreg.ne.0 .and. j.eq.mpe) then
c ---       all arctic patch tiles must be the same size or empty,
c ---       and empty tiles must be "twinned" across the top boundary.
            do i= 1,npe
              if      (iipx(      i,j).eq.0 .and. 
     &                 iipx(npe+1-i,j).ne.0      ) then
                nmpe = nmpe + 1
                ispx(i,j) = ispt(i,j)
                iipx(i,j) = iipe(i,j)
              endif
            enddo !i
          endif
        enddo !j
      endif !luniform
c
      ibig = maxval(iipx(1:npe,1:mpe))
      jbig = maxval(jjpe(      1:mpe))
c
      minsea = idm*jdm
      maxsea = 0
      mavsea = 0
      do j= 1,mpe
        do i= 1,npe
          if     (iipx(i,j).ne.0) then
            nsea   = sum( map(ispx(i,j):ispx(i,j)+iipx(i,j)-1,
     &                        jspt(  j):jspt(  j)+jjpe(  j)-1) )
            minsea = min( minsea, nsea )
            maxsea = max( maxsea, nsea )
            mavsea = mavsea + nsea
          endif
        enddo
      enddo
c
      write( 6,'(/8a6,3a8)')   '  npes','   npe','   mpe',
     &                         '   idm','   jdm','  ibig','  jbig',
     &                         '  nreg',
     &                         '  minsea','  maxsea','  avesea'
      write( 6,'(8i6,3i8/)') nmpe,npe,mpe,idm,jdm,
     &                       ibig,jbig,nreg,
     &                       minsea,maxsea,mavsea/nmpe
c
      call zhopen(21, 'formatted', 'new', 0)
      write(21,'(8a6,3a8)')    '  npes','   npe','   mpe',
     &                         '   idm','   jdm','  ibig','  jbig',
     &                         '  nreg',
     &                         '  minsea','  maxsea','  avesea'
      write(21,'(8i6,3i8/)') nmpe,npe,mpe,idm,jdm,
     &                       ibig,jbig,nreg,
     &                       minsea,maxsea,mavsea/nmpe
      do j= 1,mpe
        if     (npe.le.8) then
          write(21,'(a5,i3,a4,8i6)') 
     &         'ispt(',j,') = ',(ispx(i,j),i=1,npe)
          write(21,'(a5,i3,a4,8i6)') 
     &         'iipe(',j,') = ',(iipx(i,j),i=1,npe)
        else
          write(21,'(a5,i3,a4,8i6)') 
     &         'ispt(',j,') = ',(ispx(i,j),i=1,min(npe,8))
          write(21,'(12x,8i6)') (ispx(i,j),i=9,npe)
          write(21,'(a5,i3,a4,8i6)') 
     &         'iipe(',j,') = ',(iipx(i,j),i=1,min(npe,8))
          write(21,'(12x,8i6)') (iipx(i,j),i=9,npe)
        endif
      enddo
      write(21,*)
      if     (mpe.le.8) then
        write(21,'(a12,8i6)') 'jspt(  1) = ',(jspt(j),j=1,mpe)
        write(21,'(a12,8i6)') 'jjpe(  1) = ',(jjpe(j),j=1,mpe)
      else
        write(21,'(a12,8i6)') 'jspt(  1) = ',(jspt(j),j=1,min(mpe,8))
        write(21,'(12x,8i6)')                (jspt(j),j=9,mpe)
        write(21,'(a12,8i6)') 'jjpe(  1) = ',(jjpe(j),j=1,min(mpe,8))
        write(21,'(12x,8i6)')                (jjpe(j),j=9,mpe)
      endif
      close(21)
      end
      subroutine split1d(mapsum,jdm, jspt,jjpe,mpe)
      implicit none
c
      integer jdm,mpe
      integer mapsum(jdm),jspt(mpe),jjpe(mpe)
c
c --- load balance across mpe processors.
c
      logical update
      integer jsum(mpe), isum,jmsum,jpsum,m,maxjj,maxsum,minsum
c
c --- first guess
c
      call r_split1d(mapsum,jdm, 0,jspt,jjpe,mpe)
c
c     tile sums.
c
      do m= 1,mpe
        jsum(m) = sum( mapsum(jspt(m):jspt(m)+jjpe(m)-1) )
      enddo
c
c     try to reduce the maximum tile sum
c
      do
        maxsum = maxval(jsum)
        update = .false.
        do m= 1,mpe
          if     (jsum(m).eq.maxsum) then
            jmsum = jsum(max(m-1,  1))
            jpsum = jsum(min(m+1,mpe))
            if     (jmsum.lt.min(jsum(m),jpsum)) then
              isum = mapsum(jspt(m))
              if     (jmsum+isum.lt.maxsum) then
                jsum(m-1) = jsum(m-1) + isum
                jsum(m  ) = jsum(m  ) - isum
                jjpe(m-1) = jjpe(m-1) + 1
                jjpe(m  ) = jjpe(m  ) - 1
                jspt(m  ) = jspt(m  ) + 1
                update = .true.
*               write(6,'(a,4i8)') 'split1d: m-',
*    &                             m,jsum(m),jsum(m-1),isum
              else
*               write(6,'(a,4i8)') 'split1d: M-',
*    &                             m,jsum(m),jsum(m-1),isum
              endif
            elseif (jpsum.lt.min(jsum(m),jmsum)) then
              isum = mapsum(jspt(m)+jjpe(m)-1)
              if     (jpsum+isum.lt.maxsum) then
                jsum(m+1) = jsum(m+1) + isum
                jsum(m  ) = jsum(m  ) - isum
                jspt(m+1) = jspt(m+1) - 1
                jjpe(m+1) = jjpe(m+1) + 1
                jjpe(m  ) = jjpe(m  ) - 1
                update = .true.
*               write(6,'(a,4i8)') 'split1d: m+',
*    &                             m,jsum(m),jsum(m+1),isum
              else
*               write(6,'(a,4i8)') 'split1d: M+',
*    &                             m,jsum(m),jsum(m+1),isum
              endif
            endif
          endif
        enddo
        if (.not. update) then
          exit
        endif
      enddo
c
c     try to reduce the maximum tile size
c
      maxsum = maxval(jsum)
c
      do
        maxjj  = maxval(jjpe)
        update = .false.
        do m= 1,mpe
          if     (jjpe(m).eq.maxjj) then
            if     (m.gt.1 .and. jjpe(m-1).le.maxjj-2) then
              jmsum = jsum(m-1)
              isum  = mapsum(jspt(m))
              if     (jmsum+isum.lt.maxsum) then
                jsum(m-1) = jsum(m-1) + isum
                jsum(m  ) = jsum(m  ) - isum
                jjpe(m-1) = jjpe(m-1) + 1
                jjpe(m  ) = jjpe(m  ) - 1
                jspt(m  ) = jspt(m  ) + 1
                update = .true.
*               write(6,'(a,4i8)') 'split1D: m-',
*    &                             m,maxsum,jsum(m-1),isum
              else
*               write(6,'(a,4i8)') 'split1D: M-',
*    &                             m,maxsum,jsum(m-1),isum
              endif
            endif
            if     (m.lt.mpe .and. jjpe(m+1).le.maxjj-2) then
              jpsum = jsum(m+1)
              isum  = mapsum(jspt(m)+jjpe(m)-1)
              if     (jpsum+isum.lt.maxsum) then
                jsum(m+1) = jsum(m+1) + isum
                jsum(m  ) = jsum(m  ) - isum
                jspt(m+1) = jspt(m+1) - 1
                jjpe(m+1) = jjpe(m+1) + 1
                jjpe(m  ) = jjpe(m  ) - 1
                update = .true.
*               write(6,'(a,4i8)') 'split1D: m+',
*    &                             m,maxsum,jsum(m+1),isum
              else
*               write(6,'(a,4i8)') 'split1D: M+',
*    &                             m,maxsum,jsum(m+1),isum
              endif
            endif
          endif
        enddo
        if (.not. update) then
          exit
        endif
      enddo
      return
      end
      recursive subroutine r_split1d(mapsum,jdm, j0,jspt,jjpe,mpe)
      implicit none
c
      integer jdm,mpe,j0
      integer mapsum(jdm),jspt(mpe),jjpe(mpe)
c
c --- load balance across mpe processors.
c
      integer ipcum,iptot,iptgt,j,j01,j02,j1,j2,jd1,jd2,jp1,jp2
c
*     write(6,'(a,3i6)') 'r_split1d jdm,mpe,j0 = ',jdm,mpe,j0 
      if     (mpe.eq.1) then
        jspt(1) = j0+1
        jjpe(1) = jdm
*       write(6,'(a,2i6)') 'r_split1d jspt,jjpe  = ',jspt(1),jjpe(1)
      elseif (mpe.gt.jdm) then
c
c ---   too few rows.
c
        write(6,'(/ a /)') 'error in split2x - not enough rows'
        call zhflsh(6)
        stop
      elseif (mod(mpe,2).eq.0) then
c
c ---   even: split into two equal pieces and recurse.
c
        iptot = sum(mapsum(1:jdm))
        iptgt = iptot/2
        ipcum = 0
        do j= 1,jdm
          ipcum = ipcum + mapsum(j)
          if     (ipcum.ge.iptgt) then
c ---       must get here, since iptgt .le. ipcum(j=jdm)
*           write(6,'(a,i6,2i8,2i6)') 'r_split1d j,iptgt,ipcum  =',
*    &                                          j,iptgt,ipcum,
*    &                                          mapsum(j)-(ipcum-iptgt),
*    &                                          ipcum-iptgt
            if     (j.eq.jdm .or.
     &              ipcum-iptgt.gt.(mapsum(j)+1)/2) then
              j2 = j
              ipcum = ipcum-mapsum(j)
*             write(6,'(a,4i8)') 'r_split1d: j2=j  ',
*    &                           ipcum,iptot-ipcum,mapsum(j),iptgt-ipcum
            else
              j2 = j+1
*             write(6,'(a,4i8)') 'r_split1d: j2=j+1',
*    &                           ipcum,iptot-ipcum,mapsum(j),ipcum-iptgt
            endif
            exit
          endif
        enddo
        j1  = 1
        jd1 = j2-1
        j01 = j0
        jp1 = 1
        jd2 = jdm-j2+1
        j02 = j0+ j2-1
        jp2 = mpe/2+1
        call r_split1d(mapsum(j1),jd1, j01,jspt(jp1),jjpe(jp1),mpe/2)
        call r_split1d(mapsum(j2),jd2, j02,jspt(jp2),jjpe(jp2),mpe/2)
      else
c
c ---   odd: peel off the first processor and recurse.
c
        iptot = sum(mapsum(1:jdm))
        iptgt = iptot/mpe
        ipcum = 0
        do j= 1,jdm
          ipcum = ipcum + mapsum(j)
          if     (ipcum.ge.iptgt) then
c ---       must get here, since iptgt .le. ipcum(j=jdm)
*           write(6,'(a,i6,2i8,2i6)') 'r_split1d j,iptgt,ipcum  =',
*    &                                          j,iptgt,ipcum,
*    &                                          mapsum(j)-(ipcum-iptgt),
*    &                                          ipcum-iptgt
            if     (j.eq.jdm .or.
     &              ipcum-iptgt.gt.(3*mapsum(j))/4) then  ! greedy bias
              j2 = j
              ipcum = ipcum-mapsum(j)
*             write(6,'(a,4i8)') 'r_split1d: j2=j  ',
*    &                           ipcum,iptot-ipcum,mapsum(j),iptgt-ipcum
            else
              j2 = j+1
*             write(6,'(a,4i8)') 'r_split1d: j2=j+1',
*    &                           ipcum,iptot-ipcum,mapsum(j),ipcum-iptgt
            endif
            exit
          endif
        enddo
        j1  = 1
        jd1 = j2-1
        j01 = j0
        jp1 = 1
        jd2 = jdm-j2+1
        j02 = j0+ j2-1
        jp2 = 2
        call r_split1d(mapsum(j1),jd1, j01,jspt(jp1),jjpe(jp1),    1)
        call r_split1d(mapsum(j2),jd2, j02,jspt(jp2),jjpe(jp2),mpe-1)
      endif
      return
      end
