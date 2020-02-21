      program archv2data2d
      use mod_plot         ! HYCOM plot array interface
      use mod_za           ! HYCOM array I/O interface
c
c --- hycom to 2-d diagnostic field extractor
c
      real,    allocatable, dimension (:,:)   ::
     &   uflux,vflux, strmf,dpdx,dpdy, util1,work
      real,    allocatable, dimension (:,:,:) ::
     &   utilk,w
      real   depthu,depthv,ubi,vbi,s1,s2
      real*8 strmft,strmfu,strmfv
c
      common/conrng/ amn,amx
c
      character flnm*240,frmt*80,cline*240
      character ctrc_title(99)*80,ctrc_units(99)*80,
     &          ctrc_lname(99)*80,ctrc_sname(99)*80
      logical   smooth,mthin,lsteric,icegln,lperiod,baclin,xyward,
     &          barouv,gstruv
c
      logical          ltheta
      integer          artype,iexpt,iversn,kkin,yrflag,mxlflg
      real             dudxdn,dudxup,dvdydn,dvdyup,dp00f


      double precision time3(3)
      double precision dsumth,dsumdp
c
      real, parameter :: flag = 2.0**100
c
c --- 'trcout' -- tracer input
      logical   trcout,dbg
      data      trcout/.false./
c
      real      tenm,onem,temcm,onecm,onemm
      data      tenm/10./,onem/1./,tencm/.1/,onecm/.01/,onemm/.001/
c
      logical   initl
      data      initl /.true. /
      real      thref,spcifh
      data      thref/1.e-3/,spcifh/3990./
      character blank*40
      data      blank/'                                        '/
c
      call xcspmd
      call zaiost
      call blkinit
      lp=6
c
c --- read model data
c ---   'flnm  ' = name of file containing the actual data
c ---   'frmt  ' = output format or type (HYCOM, BINARY, netCDF)
c ---                see horout for more details on frmt
c ---   'iexpt ' = experiment number x10  (000=from archive file)
c ---   'yrflag' = days in year flag (0=360J16,1=366J16,2=366J01,3=actual)
c ---   'ntracr' = number of tracers (to output, optional with default 0)
c ---    one name line per tracer: 8-letter plot and units, field, standard_
c ---      or separated by "|" (i.e. plot|units|field|standard).
c ---      the field name must only contain alphanumerics and "_", and 
c ---      the standard_name is either blank or from the CF 1.0 conventions
c ---   'idm   ' = longitudinal array size
c ---   'jdm   ' = latitudinal  array size
c ---   'kdm   ' = number of layers
        read (uoff+99,'(a)') flnm
        if(mnproc.eq.1)then
          write (lp,'(2a)') ' input file: ',trim(flnm)
          call flush(lp)
        endif
        read (uoff+99,'(a)') frmt
        if(mnproc.eq.1)then
        write (lp,'(2a)') 'output type: ',trim(frmt)
        call flush(lp)
        endif
        call blkini(iexpt, 'iexpt ')
        call blkini(yrflag,'yrflag')
        call blkini2(i,j,  'ntracr','idm   ')  !read ntracr or idm
        if (j.eq.1) then
          ntracr = i
          do ktr= 1,ntracr
            read(uoff+99,'(a)') cline
            i = index(cline,'|')
            if     (i.eq.0) then  !8-letter plot and units, field has no spaces
              ctrc_title(ktr) = cline(1:8)
              ctrc_units(ktr) = cline(9:16)
              cline = cline(17:)
              do
                i = index(cline,' ')
                if     (i.ne.1) then
                  exit
                endif
                cline = cline(2:) !remove a leading space
              enddo
              ctrc_lname(ktr) = cline(1:i-1)
              ctrc_sname(ktr) = cline(i+1:)
            else  !separated by "|" 
              ctrc_title(ktr) = cline(1:i-1)
              cline = cline(i+1:)
              i = index(cline,'|')
              ctrc_units(ktr) = cline(1:i-1)
              cline = cline(i+1:)
              i = index(cline,'|')
              ctrc_lname(ktr) = cline(1:i-1)
              ctrc_sname(ktr) = cline(i+1:)
            endif
        if(mnproc.eq.1)then
              write (lp,'(2x,i2,3a)')
     &        ktr,' title  = "',trim(ctrc_title(ktr)),'"',
     &        ktr,' units  = "',trim(ctrc_units(ktr)),'"',
     &        ktr,' l.name = "',trim(ctrc_lname(ktr)),'"',
     &        ktr,' s.name = "',trim(ctrc_sname(ktr)),'"'
            call flush(lp)
      endif
            if     (   index(ctrc_lname(ktr),' ').ne.0 .and.
     &                 index(ctrc_lname(ktr),' ').le.
     &              len_trim(ctrc_lname(ktr))                ) then
              ! does not catch all illegal l.names.
        if(mnproc.eq.1)then
              write(lp,*)
              write(lp,*) 'error - l.name contains spaces'
              write(lp,*)
              call flush(lp)
        endif
              call xcstop('(archv - l.name wrong)')
            elseif (   index(ctrc_lname(ktr),'-').ne.0) then
              ! still does not catch all illegal l.names.
        if(mnproc.eq.1)then
              write(lp,*)
              write(lp,*) 'error - l.name contains "-"'
              write(lp,*)
              call flush(lp)
      endif
              call xcstop('(archv - l.name wrong)')
            endif !l.name check
          enddo
          
          call blkini(iit,  'idm   ')
        else
          ntracr = 0
          iit     = i
        endif
        call blkini(jjt,    'jdm   ')
        call blkini(kkt,    'kdm   ')
        kk=kkt
        kkmax=kk
        if     (iit.ne.itdm .or. jjt.ne.jtdm) then
          if(mnproc.eq.1)then
            write(lp,*)
            write(lp,*) 'error - wrong idm or jdm (should be:',
     &                                            itdm,jtdm,')'                  
            write(lp,*)
            call flush(lp)
          endif
          call xcstop('(archv - wrong idm or jdm)')
        endif

c
c ---   'thbase' = reference density (sigma units)
        call blkinr(thbase,
     &             'thbase','("blkinr: ",a6," =",f11.4," sig")')
c
c ---   'smooth' = smooth fields before output
c ---   'mthin ' = mask thin layers from output
        call blkinl(smooth,'smooth')
        call blkinl(mthin, 'mthin ')
c
c ---   'baclin' = extract baroclinic velocity (0=total:DEFAULT,1=baroclinic)
c ---   'xyward' = output original unrotated velocities (0=no:DEFAULT,1=yes)
        call blkini3(i,j,  'baclin','xyward','iorign')  !read one of three
        if (j.eq.1) then
          baclin = i.eq.1  !baroclinic, vs total, velocity
          call blkini2(i,j,  'xyward','iorign')  !read one of two
          if (j.eq.1) then
            xyward = i.eq.1  !x-y, vs east-north, velocity
            call blkini(iorign,'iorign')
          else
            xyward = .false. !eastward,northward velocity
            iorign = i
          endif
        elseif (j.eq.2) then
          baclin = .false. !total velocity
          xyward = i.eq.1  !x-y, vs east-north, velocity
          call blkini(iorign,'iorign')
        else
          baclin = .false. !total velocity
          xyward = .false. !eastward,northward velocity
          iorign = i
        endif
c
c ---   'iorign' = i-origin of sampled subregion
c ---   'jorign' = j-origin of sampled subregion
c ---   'idmp  ' = i-extent of sampled subregion (<=idm; 0 implies idm)
c ---   'jdmp  ' = j-extent of sampled subregion (<=jdm; 0 implies jdm)
c----    exit if  idmp .ne. jtdm .ne. 0
*       call blkini(iorign,'iorign')  !see above
        call blkini(jorign,'jorign')
        call blkini(iit,    'idmp  ')
        call blkini(jjt,    'jdmp  ')
        if     (iit.eq.0) then
c          ii=idm
        else
           if(mnproc.eq.1)then
        write(lp,*)'subregion version not supported: idmp .ne. 0,=',iit
        call flush(lp) 
           endif
         call xcstop('(archv2data2d_mpi - Subregion error)')
        endif
        if     (jjt.eq.0) then
c          jj=jdm
        else
           if(mnproc.eq.1)then
        write(lp,*)'subregion version not supported: jdmp .ne. 0,=',jjt
        call flush(lp) 
           endif
         call xcstop('(archv2data2d_mpi - Subregion error )')
        endif

c ---   'iorign,jorign' denote the origin of the subgrid to be extracted 
c ---   from the full history grid (dimensioned idm x jdm). 
c ---   The size of the subgrid is determined by ii,jj.

        if(mnproc.eq.1)then
          write (lp,'(/ 2(a,i5),9x,2(a,i5) /)') 'extracting i =',iorign,
     &    ' ...',iorign+itdm-1,'j =',jorign,' ...',jorign+jtdm-1
          call flush(lp)
        endif
c
        call getartype(flnm,artype)
c
c --- array allocation
c
      call plot_alloc
c
      allocate(  uflux(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
      allocate(  vflux(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
      allocate(  strmf(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
      allocate(   work(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
      allocate(  util1(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
      allocate(  dpdx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
      allocate(  dpdy(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
c
      allocate(  utilk(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,1:kk) )
c
      dpthfil = 'regional.depth'
c

      do j=1-nbdy,jj+nbdy
        do i=1-nbdy,ii+nbdy
          p(i,j,1)=0.
        enddo
      enddo
c
c --- read the archive file.
c
        if     (artype.ne.3) then
          call getdat( flnm,time3,iweight,mntype,lsteric,icegln,trcout,
     &                 iexpt,iversn,yrflag,kkin)     ! hycom input
          time = time3(3)
        else
          if(mnproc.eq.1)then
            write(lp,*)'Attempted to read stddev Archive - Line 265'
            call flush(lp)
          endif
          call xcstop('(archv3data2d_mpi stddev Archive error)')
        endif
        if (kkin.ne.kk) then
          if(mnproc.eq.1)then
            write(lp,*)
            write(lp,*) 'error - kkin must be kdm'
            write(lp,*)
            call flush(lp)
          endif  
          call xcstop('(archv2data2d - kk)')
        endif
c
      if     (yrflag.eq.0) then
        year  = 360.0d0
      elseif (yrflag.lt.3) then
        year  = 366.0d0
      else
        year  = 365.25d0
      endif
      
      
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C  Serial getdat  called getdepth
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call getdepth(dpthfil)
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C   Ensure Geometrical factor arrays are tiled!
c
      call xctilr(scvx,1,1,nbdy,nbdy,halo_vv) 
      call xctilr(scpx,1,1,nbdy,nbdy,halo_ps)         
c
c --- define grid scale
c
c*****************************************************************
C  Patch to get global max and min for plat  and plon
c  Dan Moore QinetiQ  --  July 2010
c
      platmax=maxval(plat(1:ii,1:jj))
      platmin=minval(plat(1:ii,1:jj))
      plonmax=maxval(plon(1:ii,1:jj))
      plonmin=minval(plon(1:ii,1:jj))

      call xcmaxr(platmax)
      call xcminr(platmin)
      call xcmaxr(plonmax)
      call xcminr(plonmin)

        if(mnproc.eq.1)then
          write(lp,'(/a,2f8.2/a,2f8.2)') 
     &     'sub-domain longitude range = ',
     &    plonmin,plonmax,
     &     'sub-domain latitude  range = ',
     &    platmin,platmax
        endif
c
      lperiod = ii.eq.idm .and.
     &       plonmax-plonmin .gt. 350.0
      if     (lperiod) then
        if(mnproc.eq.1)then
          write(lp,'(/a/)') 'sub-domain assumed to be periodic'
        endif
      else
        if(mnproc.eq.1)then
          write(lp,'(/a/)') 'sub-domain assumed to be non-periodic'
        endif
      endif
c*******************************************************************
c
      call bigrid(depths)
c
c --- check that bathymetry is consistent with this archive.
c --- only possible with hycom .[ab] file input.
c
      if     (iversn.ge.20) then
        ibadl = 0
        ibads = 0
        do j= 1,jj
          do i= 1,ii
            if     (ip(i,j).eq.1) then
              if     (srfht(i,j).gt.2.0**99) then
                ibads = ibads + 1   ! topo sea, srfht land
              endif
            else
              if     (srfht(i,j).lt.2.0**99) then
                ibadl = ibadl + 1   ! topo land, srfht sea
               endif
            endif
          enddo !i
        enddo !j
        if     (ibads.ne.0) then
        if(mnproc.eq.1)then
          write(lp,*)
          write(lp,*) 'error - wrong bathymetry for this archive file'
          write(lp,*) 'number of topo sea  mismatches = ',ibads
          write(lp,*) 'number of topo land mismatches = ',ibadl
          write(lp,*)
          call flush(lp)
        endif
          call xcstop('(archv2data2d - topo mismatch)')
        endif !ibads.ne.0
        if     (ibadl.ne.0) then
        if(mnproc.eq.1)then
          write(lp,*)
*         write(lp,*) 'error - wrong bathymetry for this archive file'
          write(lp,*) 'warning - wrong bathymetry for this archive file'
          write(lp,*) 'number of topo sea  mismatches = ',ibads
          write(lp,*) 'number of topo land mismatches = ',ibadl
          write(lp,*)
          call flush(lp)
        endif
*         call xcstop
        endif !ibadl.ne.0
      endif !iversn.ge.20  
c
      do 3 k=1,kkin
      do 3 j=1,jj
      do 3 i=1,ii
c
c --- convert baroclinic to total velocities by adding barotropic component
c --- note that mean archives already contain total velocity
      if (k.eq.1) then
        if     (iu(i,j).eq.1 .and. artype.eq.1) then
          umix(i,j)=umix(i,j)+ubaro(i,j)  !ignore baclin, always total
        elseif (iu(i,j).ne.1) then
          umix(i,j)=0.
        end if
        if     (iv(i,j).eq.1 .and. artype.eq.1) then
          vmix(i,j)=vmix(i,j)+vbaro(i,j)  !ignore baclin, always total
        elseif (iv(i,j).ne.1) then
          vmix(i,j)=0.
        end if
      endif
      if     (iu(i,j).eq.1) then
        if     (artype.eq.1 .and. .not.baclin) then
          u(i,j,k)=u(i,j,k)+ubaro(i,j)  !total velocity
        elseif (artype.eq.2 .and.      baclin) then
          u(i,j,k)=u(i,j,k)-ubaro(i,j)  !baroclinic velocity
        end if
      else !iu(i,j).ne.1
        u(i,j,k)=0.
      end if
      if     (iv(i,j).eq.1) then
        if     (artype.eq.1 .and. .not.baclin) then
          v(i,j,k)=v(i,j,k)+vbaro(i,j)  !total velocity
        elseif (artype.eq.2 .and.      baclin) then
          v(i,j,k)=v(i,j,k)-vbaro(i,j)  !baroclinic velocity
        end if
      else !iv(i,j).ne.1
        v(i,j,k)=0.
      end if
c
c --- convert layer thickness to meters
      if (depths(i,j).gt.0.) then
        dp(i,j,k)=dp(i,j,k)/9806.
        p(i,j,k+1)=p(i,j,k)+dp(i,j,k)
        if     (artype.eq.3) then
          dpsd(i,j,k)=dpsd(i,j,k)/9806.
        else
          th3d(i,j,k)=th3d(i,j,k)+thbase
        endif
      else
        saln(i,j,k)=flag
        temp(i,j,k)=flag
        th3d(i,j,k)=flag
        dp(i,j,k)=flag
        p(i,j,k+1)=flag
        if     (artype.eq.3) then
          dpsd(i,j,k)=flag
          ke(  i,j,k)=flag
        elseif (artype.eq.2) then
          ke(  i,j,k)=flag
        endif

        do ktr= 1,ntracr
          trcr(i,j,k,ktr)=flag
        enddo
      endif
 3    continue

        call xctilr( p(1-nbdy,1-nbdy,1),1,kk+1,nbdy,nbdy, halo_ps)
        call xctilr(dp(1-nbdy,1-nbdy,1),1,kk  ,nbdy,nbdy, halo_ps)
        call xctilr( u(1-nbdy,1-nbdy,1),1,kk  ,nbdy,nbdy, halo_uv)
        call xctilr( v(1-nbdy,1-nbdy,1),1,kk  ,nbdy,nbdy, halo_vv)
        call xctilr(th3d(1-nbdy,1-nbdy,1),1,kk,   nbdy,nbdy, halo_ps)
        call xctilr(temp(1-nbdy,1-nbdy,1),1,kk,   nbdy,nbdy, halo_ps)
        call xctilr(saln(1-nbdy,1-nbdy,1),1,kk,   nbdy,nbdy, halo_ps)
                     
c
c --- eddy kinetic energy
      if     (artype.eq.3) then
        do k=1,kkin
          do j= 1,jj
c            jp1 = min(j+1,jj)
            jp1 = min(j+1+j0,jtdm)-j0
            do i= 1,ii
              if (depths(i,j).gt.0.) then
                if     (i.ne.ii) then
                  ip1 = i+1
                elseif (lperiod) then !i=ii
                  ip1 = 1
                else !i=ii (non-periodic)
                  ip1 = ii
                endif
c ---           ke = 0.5*( std(u)**2 + std(v)**2 )
                ke(i,j,k)= 0.125*((u(i,j,k)+u(ip1,j,k))**2 +
     &                            (v(i,j,k)+v(i,jp1,k))**2  )
              endif
            enddo !i
          enddo !j
        enddo !k
      endif !std archive
c
c    Set limits for j for somothings and scalings
c
      jstop=jj
      if(j0+jj.ge.jtdm)jstop=jj-1
      jstart=1
      if(mnproc.eq.1)jstart=2

c --- fluid vertical velocity
c

      if     (artype.ne.3) then  !not available for std. case
      allocate(   w(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,1:kk) )
      do j= 1,jstop
        do i= 1,ii-1
          if     (ip(i,j).eq.1) then
            dudxdn=
     &            (u(i+1,j  ,1)*scuy(i+1,j  )-u(i,j,1)*scuy(i,j))
     &            /(scpx(i,j)*scpy(i,j))
            dvdydn=
     &            (v(i  ,j+1,1)*scvx(i  ,j+1)-v(i,j,1)*scvx(i,j))
     &        /(scpx(i,j)*scpy(i,j))
            w(i,j,1)=0.5*dp(i,j,1)*(dudxdn+dvdydn)
          else
            w(i,j,1) = flag
          endif
        enddo
      enddo
        
      do k= 2,kkin
        do j= 1,jstop
          do i= 1,ii-1
            if     (iu(i,j).eq.1) then
              dpdx(i,j)=
     &                (p(i,j,k)*scpy(i,j)-p(i-1,j  ,k)*scpy(i-1,j  ))
     &                /(scux(i,j)*scuy(i,j))
            endif
          enddo
        enddo
        do j= 1,jstop
          do i= 1,ii-1   
            if     (iv(i,j).eq.1) then
              dpdy(i,j)=
     &                (p(i,j,k)*scpx(i,j)-p(i  ,j-1,k)*scpx(i  ,j-1))
     &                /(scvx(i,j)*scvy(i,j))
            endif
          enddo
        enddo
        do j=1,jj
          if     (iu(2   ,j).eq.1) then
            dpdx(1 ,j)=dpdx(2   ,j)
          endif
          if     (iu(ii-1,j).eq.1) then
            dpdx(ii,j)=dpdx(ii-1,j)
          endif
        enddo
        do i=1,ii
          if(iv(i,2).eq.1.and.mnproc.eq.1) then
            dpdy(i,1 )=dpdy(i,2   )
          endif
          if(iv(i,jj-1).eq.1.and.jj+j0.ge.jtdm) then
            dpdy(i,jj)=dpdy(i,jj-1)
          endif
        enddo
        call xctilr(dpdy,1,1,nbdy,nbdy,halo_ps)
        do j= 1,jstop
          do i= 1,ii-1
            if     (ip(i,j).eq.1) then
              dudxup=
     &              (u(i+1,j  ,k-1)*scuy(i+1,j  )-
     &               u(i  ,j  ,k-1)*scuy(i  ,j  ))
     &              /(scpx(i,j)*scpy(i,j))
              dvdyup=
     &              (v(i  ,j+1,k-1)*scvx(i  ,j+1)-
     &               v(i  ,j  ,k-1)*scvx(i  ,j  ))
     &              /(scpx(i,j)*scpy(i,j))
              dudxdn=
     &              (u(i+1,j  ,k  )*scuy(i+1,j  )-
     &               u(i  ,j  ,k  )*scuy(i  ,j  ))
     &              /(scpx(i,j)*scpy(i,j))
              dvdydn=
     &              (v(i  ,j+1,k  )*scvx(i  ,j+1)-
     &               v(i  ,j  ,k  )*scvx(i  ,j  ))
     &              /(scpx(i,j)*scpy(i,j))
              w(i,j,k)=w(i,j,k-1)+0.5*(dp(i,j,k-1)*(dudxup+dvdyup)+
     &                                 dp(i,j,k  )*(dudxdn+dvdydn)-
     &                  (u(i  ,j  ,k)-u(i  ,j  ,k-1))*dpdx(i  ,j  )-
     &                  (u(i+1,j  ,k)-u(i+1,j  ,k-1))*dpdx(i+1,j  )-
     &                  (v(i  ,j  ,k)-v(i  ,j  ,k-1))*dpdy(i  ,j  )-
     &                  (v(i  ,j+1,k)-v(i  ,j+1,k-1))*dpdy(i  ,j+1))
            else
              w(i,j,k) = flag
            endif
          enddo
        enddo
        if(jj+j0.eq.jtdm)then
          do i= 1,ii
            w(i ,jj,k) = flag
          enddo
        endif  
        do j= 1,jj
          w(ii,j ,k) = flag
        enddo
      enddo  ! k
c
c --- w is noisy - smooth at least once.
c
      if     (smooth) then
        do k= 1,kkin
          call xctilr(w(1-nbdy,1-nbdy,k),1,1,nbdy,nbdy, halo_ps)
          call psmoo(w(1-nbdy,1-nbdy,k),work)
          call xctilr(w(1-nbdy,1-nbdy,k),1,1,nbdy,nbdy, halo_ps)
          call psmoo(w(1-nbdy,1-nbdy,k),work)
          call xctilr(w(1-nbdy,1-nbdy,k),1,1,nbdy,nbdy, halo_ps)
        enddo
      else
        do k= 1,kkin 
          call xctilr(w(1-nbdy,1-nbdy,k),1,1,nbdy,nbdy, halo_ps)
          call psmoo(w(1-nbdy,1-nbdy,k),work)
          call xctilr(w(1-nbdy,1-nbdy,k),1,1,nbdy,nbdy, halo_ps)
        enddo
      endif
      endif !artype.ne.3

c
      do 7 j=1,jj
c      jp1 = min(j+1,jj)
      jp1 = min(j+1+j0,jtdm)-j0
      do 7 i=1,ii
      if     (ip(i,j).eq.1) then
        if     (artype.eq.3) then
          if     (i.ne.ii) then
            ip1 = i+1
          elseif (lperiod) then !i=ii
            ip1 = 1
          else !i=ii (non-periodic)
            ip1 = ii
          endif     
c ---     ke = 0.5*( std(u)**2 + std(v)**2 )
          kemix( i,j)=0.125*((umix( i,j)+umix( ip1,j))**2 +
     &                       (vmix( i,j)+vmix( i,jp1))**2  )
          kebaro(i,j)=0.125*((ubaro(i,j)+ubaro(ip1,j))**2 +
     &                       (vbaro(i,j)+vbaro(i,jp1))**2  )
        endif !std archive
        dpbl(  i,j)=dpbl(  i,j)/9806.          ! m
        dpmixl(i,j)=dpmixl(i,j)/9806.          ! m
        thmix( i,j)=thmix( i,j)+thbase         ! SigmaT
        if     (artype.ne.3) then
          ttrend(i,j)=surflx(i,j)*thref*8.64E4
     &                    /spcifh/max(0.1,dpbl(i,j))    ! deg/day
          strend(i,j)=salflx(i,j)*thref*8.64E4
     &                           /max(0.1,dpbl(i,j))    ! psu/day
          emnp(  i,j)=-wtrflx(i,j)                      ! kg/m^2/s into the atmos
        else  ! std.dev, archive
          ttrend(i,j)=flag
          strend(i,j)=flag
          emnp(  i,j)=-wtrflx(i,j)                      ! kg/m^2/s into the atmos
        endif
        if     (covice(i,j).eq.0.0) then
          thkice(i,j)= 0.0
          temice(i,j)=-1.8
        endif
      else
        steric(i,j)=flag
        srfht( i,j)=flag
        montg( i,j)=flag
        surflx(i,j)=flag
        salflx(i,j)=flag
        ttrend(i,j)=flag
        strend(i,j)=flag
          emnp(i,j)=flag
        covice(i,j)=flag
        thkice(i,j)=flag
        temice(i,j)=flag
        dpbl(  i,j)=flag
        dpmixl(i,j)=flag
        tmix( i,j)=flag
        smix( i,j)=flag
        thmix(i,j)=flag
        if     (artype.gt.1) then
          kemix( i,j)=flag
          kebaro(i,j)=flag
        endif
      end if
      if (kkin.eq.1 .or. kkin.lt.kk) then
        if     (ip(i,j).eq.1) then
          p(i,j,kk+1)=depths(i,j)
        else
          p(i,j,kk+1)=flag
        endif
      endif
 7    continue

c
      dpth=0.5*onecm

c        
      if (smooth) then
c
c --- smooth mass field variables
c
      call xctilr(temp,1,1,nbdy,nbdy, halo_ps)
      call xctilr(saln,1,1,nbdy,nbdy, halo_ps)
      call xctilr(th3d,1,1,nbdy,nbdy, halo_ps)
      call psmoo(temp,work)
      call psmoo(saln,work)
      call psmoo(th3d,work)

      do ktr= 1,ntracr      
      call xctilr(trcr(1-nbdy,1-nbdy,1,ktr),1,1,nbdy,nbdy, halo_ps)
        call psmoo(trcr(1-nbdy,1-nbdy,1,ktr),work) 
      enddo

      do 3800 k=2,kkin
c
      do 76 j=1,jj
      do 76 i=1,ii
      if (depths(i,j).gt.0.) then
        util1(i,j)=max(onemm,dp(i,j,k))
        temp(i,j,k)=temp(i,j,k)*util1(i,j)
        saln(i,j,k)=saln(i,j,k)*util1(i,j)
        th3d(i,j,k)=th3d(i,j,k)*util1(i,j)
        do ktr= 1,ntracr
          trcr(i,j,k,ktr)=trcr(i,j,k,ktr)*util1(i,j)
        enddo
      else
        temp(i,j,k)=flag
        saln(i,j,k)=flag
        th3d(i,j,k)=flag
        do ktr= 1,ntracr
          trcr(i,j,k,ktr)=flag
        enddo
      end if
 76   continue
c  
      call xctilr(util1,1,1,nbdy,nbdy, halo_ps)
      call xctilr(temp(1-nbdy,1-nbdy,k),1,1,nbdy,nbdy, halo_ps)
      call xctilr(saln(1-nbdy,1-nbdy,k),1,1,nbdy,nbdy, halo_ps)
      call xctilr(th3d(1-nbdy,1-nbdy,k),1,1,nbdy,nbdy, halo_ps)

      call psmoo(util1,work)
      call psmoo(temp(1-nbdy,1-nbdy,k),work)
      call psmoo(saln(1-nbdy,1-nbdy,k),work)
      call psmoo(th3d(1-nbdy,1-nbdy,k),work)

      do ktr= 1,ntracr    
        call xctilr(trcr(1-nbdy,1-nbdy,k,ktr),1,1,nbdy,nbdy, halo_ps)
        call psmoo(trcr(1-nbdy,1-nbdy,k,ktr),work) 
      enddo
c
      do 38 j=1,jj
      do 38 i=1,ii
      if (depths(i,j).gt.0.) then
        temp(i,j,k)=temp(i,j,k)/util1(i,j)
        saln(i,j,k)=saln(i,j,k)/util1(i,j)
        th3d(i,j,k)=th3d(i,j,k)/util1(i,j)
        do ktr= 1,ntracr
          trcr(i,j,k,ktr)=trcr(i,j,k,ktr)/util1(i,j)
        enddo
      end if

   38  continue
 3800   continue  ! K loop
   
   

c
c --- smooth velocity and layer thickness fields
c
      do 30 k=1,kkin
c
      do 31 j=1,jstop
      do 31 i=2,ii1
 31   uflux(i,j)=u(i,j,k)*max(onecm,dp(i,j,k)+dp(i-1,j,k))
c
      do 32 j=jstart,jstop
      do 32 i=1,ii1
 32   vflux(i,j)=v(i,j,k)*max(onecm,dp(i,j,k)+dp(i,j-1,k))
      
      call xctilr(uflux,1,1,nbdy,nbdy, halo_uv)
      call usmoo(uflux,work)

      call xctilr(vflux,1,1,nbdy,nbdy, halo_vv)
      call vsmoo(vflux,work)

      call xctilr(dp(1-nbdy,1-nbdy,k),1,1,nbdy,nbdy, halo_ps)
      call psmoo(dp(1-nbdy,1-nbdy,k),work)
      call xctilr(dp(1-nbdy,1-nbdy,k),1,1,nbdy,nbdy, halo_ps)
c      
c --- (warning: smoothed -dp- field unsuitable for deriving interface depths)
c   
      do 33 j=1,jstop
      do 33 i=2,ii1
 33   u(i,j,k)=uflux(i,j)/max(onecm,dp(i,j,k)+dp(i-1,j,k))
c     
      do 34 j=jstart,jstop
      do 34 i=1,ii1
 34   v(i,j,k)=vflux(i,j)/max(onecm,dp(i,j,k)+dp(i,j-1,k)) 
c
c --- now smooth layer interfaces and find corresponding -dp- field
c
      if (k.lt.kkin) then
        call xctilr(p(1-nbdy,1-nbdy,k+1),1,1,nbdy,nbdy, halo_ps)
        call psmo1(p(1-nbdy,1-nbdy,k+1),work,p(1-nbdy,1-nbdy,kk+1))
        call xctilr(p(1-nbdy,1-nbdy,k+1),1,1,nbdy,nbdy, halo_ps)
      endif
      
c --- now smooth mixed layer base
      if (k.eq.1) then
        call xctilr(dpbl,1,1,nbdy,nbdy, halo_ps)
        call psmo1(dpbl,work,p(1-nbdy,1-nbdy,kk+1))
        call xctilr(dpbl,1,1,nbdy,nbdy, halo_ps)
        call xctilr(dpmixl,1,1,nbdy,nbdy, halo_ps)
        call psmo1(dpmixl,work,p(1-nbdy,1-nbdy,kk+1))
        call xctilr(dpmixl,1,1,nbdy,nbdy, halo_ps)
      end if
      
c  
      do 35 j=1,jstop
      do 35 i=1,ii1
      if (depths(i,j).gt.0.) dp(i,j,k)=p(i,j,k+1)-p(i,j,k)
 35   continue
c
 30   continue
c
      end if  !smooth = .true.
c
c --- -------------------------
c --- output non-layered fields
c --- -------------------------
c
      k=0
      ltheta=.false.
c
c --- 'areaio' = grid area  I/O unit (0 no I/O) - optional
c --- 'zlayio' = z-layer    I/O unit (0 no I/O) - optional
c --- 'dp00f ' = z-level spacing stretching factor (only for zlayo/=0)
c --- 'botio ' = bathymetry I/O unit (0 no I/O)
      call blkini3(ioin,j,  'areaio','zlayio','botio ')  !read one of three
      if (j.eq.1) then
        if (ioin.gt.0) then
          do j=1,jj
            do i=1,ii
              if (ip(i,j).ne.0) then
                util1(i,j)=scpx(i,j)*scpy(i,j)
              else
                util1(i,j)=flag
              endif
            enddo
          enddo
          call horout(util1, artype,yrflag,time3,iexpt,.true.,
     &                ' grid cell area   ',       ! plot name
     &                'area',                     ! ncdf name
     &                'sea_area',                 ! ncdf standard_name
     &                'm2',                       ! units
     &                k,ltheta, frmt,ioin)
        endif
        call blkini(ioin,'botio ')  !can't have areaio and zlayio
      elseif (j.eq.2) then
        if (ioin.gt.0) then
          call blkinr(dp00f,  
     &               'dp00f ','("blkinr: ",a6," =",f11.4," ")')
          do j=1,jj
            do i=1,ii
              if (ip(i,j).ne.0) then
                util1(i,j)=kk
                 work(i,j)=p(i,j,kk+1)
                do k=2,kk
                  if     (dp(i,j,k).gt.1.1*dp00f*dp(i,j,k-1)) then
                    util1(i,j)=k-1
                     work(i,j)=p(i,j,k)
                    exit
                  endif
                enddo !k
              else
                util1(i,j)=flag
                 work(i,j)=flag
              endif
            enddo
          enddo
          k = 0
          call horout(util1, artype,yrflag,time3,iexpt,.true.,
     &                'no. fixed layers  ',       ! plot name
     &                'zlay',                     ! ncdf name
     &                ' ',                        ! ncdf standard_name
     &                ' ',                        ! units
     &                k,ltheta, frmt,ioin)
          call horout(work,  artype,yrflag,time3,iexpt,.true.,
     &                'depth fixed layers',       ! plot name
     &                'zldepth',                  ! ncdf name
     &                ' ',                        ! ncdf standard_name
     &                'm',                        ! units
     &                k,ltheta, frmt,ioin)
        endif
        call blkini(ioin,'botio ')  !can't have zlayio and areaio
      endif !areaio:zlayio
      if (ioin.gt.0) then
        do j=1,jj
          do i=1,ii
            util1(i,j)=p(i,j,kk+1)
          enddo
        enddo
        call horout(util1, artype,yrflag,time3,iexpt,.true.,
     &              ' bathymetry       ',       ! plot name
     &              'bathymetry',               ! ncdf name
     &              'sea_floor_depth',          ! ncdf standard_name
     &              'm',                        ! units
     &              k,ltheta, frmt,ioin)
      endif
c
c --- -------------------
c ---  surface fluxes
c --- -------------------
c
c --- 'flxio ' = surf. heat  flux I/O unit (0 no I/O)
      call blkini(ioin,'flxio ')
      if (ioin.gt.0) then
        call horout(surflx,artype,yrflag,time3,iexpt,.true.,
     &              ' surf. heat flux  ',                ! plot name
     &              'qtot',                              ! ncdf name (mersea)
     &              'surface_downward_heat_flux_in_air', ! ncdf standard_name
     &              'w/m2',                              ! units
     &              k,ltheta, frmt,ioin)
      endif
c --- 'empio ' = surf. evap-pcip I/O unit (0 no I/O)
      call blkini(ioin,'empio ')
      if (ioin.gt.0) then
        do j=1,jj
          do i=1,ii
            if (ip(i,j).ne.0) then
              util1(i,j)=-emnp(i,j) !water_flux_into_ocean
            else
              util1(i,j)=flag
            endif
          enddo
        enddo
        call horout( util1,artype,yrflag,time3,iexpt,.true.,
     &              ' surf. water flux ',                 ! plot name
     &              'emp',                                ! ncdf name (mersea)
     &              'water_flux_into_ocean',              ! ncdf standard_name
     &              'kg/m2/s',                            ! units
     &              k,ltheta, frmt,ioin)
      endif
      call blkini3(ioin,i, 'ttrio ','tbfio ','txio  ')  !read one of three
      if     (i.eq.1) then
c ---   'ttrio ' = surf. temp trend I/O unit (0 no I/O)
        if (ioin.gt.0) then
          call horout(ttrend,artype,yrflag,time3,iexpt,.true.,
     &                ' surf. temp. trend',         ! plot name
     &                'surface_temperature_trend',  ! ncdf name
     &                ' ',                          ! ncdf standard_name
     &                'degC/day',                   ! units
     &                k,ltheta, frmt,ioin)
        endif
c ---   'strio ' = surf. saln trend I/O unit (0 no I/O)
        call blkini(ioin,'strio ')
        if (ioin.gt.0) then
          call horout(strend,artype,yrflag,time3,iexpt,.true.,
     &                ' surf. saln. trend',       ! plot name
     &                'surface_salinity_trend',   ! ncdf name
     &                ' ',                        ! ncdf standard_name
     &                'psu/day',                  ! units
     &                k,ltheta, frmt,ioin)
        endif
      elseif (i.eq.2) then
c ---   'tbfio ' = temp buoyancy flux I/O unit (0 no I/O)
        if (ioin.gt.0) then
          do j=1,jj
            do i=1,ii
              if (ip(i,j).ne.0) then
c               using the serial buoflx for each i,j is not efficient
                call buoflx(util1(i,j),ip(i,j),surflx(i,j),wtrflx(i,j),
     &                      temp(i,j,1),saln(i,j,1),1,1, 1)
              else
                util1(i,j)=flag
              endif
            enddo
          enddo
          call horout(util1, artype,yrflag,time3,iexpt,.true.,
     &                ' temp. bouy. flux ',         ! plot name
     &                'surface_t_bouyancy_flux',    ! ncdf name
     &                ' ',                          ! ncdf standard_name
     &                'm2/s3',                      ! units
     &                k,ltheta, frmt,ioin)
        endif
c ---   'sbfio ' = saln buoyancy flux I/O unit (0 no I/O)
        call blkini(ioin,'sbfio ')
        if (ioin.gt.0) then
          do j=1,jj
            do i=1,ii
              if (ip(i,j).ne.0) then
c               using the serial buoflx for each i,j is not efficient
                call buoflx(util1(i,j),ip(i,j),surflx(i,j),wtrflx(i,j),
     &                      temp(i,j,1),saln(i,j,1),1,1, 2)
              else
                util1(i,j)=flag
              endif
            enddo
          enddo
          call horout(util1, artype,yrflag,time3,iexpt,.true.,
     &                ' saln. bouy. flux ',         ! plot name
     &                'surface_s_bouyancy_flux',    ! ncdf name
     &                ' ',                          ! ncdf standard_name
     &                'm2/s3',                      ! units
     &                k,ltheta, frmt,ioin)
        endif
c ---   'abfio ' = tot. buoyancy flux I/O unit (0 no I/O)
        call blkini(ioin,'abfio ')
        if (ioin.gt.0) then
          do j=1,jj
            do i=1,ii
              if (ip(i,j).ne.0) then
c               using the serial buoflx for each i,j is not efficient
                call buoflx(util1(i,j),ip(i,j),surflx(i,j),wtrflx(i,j),
     &                      temp(i,j,1),saln(i,j,1),1,1, 3)
              else
                util1(i,j)=flag
              endif
            enddo
          enddo
          call horout(util1, artype,yrflag,time3,iexpt,.true.,
     &                ' total bouy. flux ',         ! plot name
     &                'surface_bouyancy_flux',      ! ncdf name
     &                ' ',                          ! ncdf standard_name
     &                'm2/s3',                      ! units
     &                k,ltheta, frmt,ioin)
        endif
      else
c ---   'txio  ' = surf. x-stress I/O unit (0 no I/O)
        if (ioin.gt.0) then
          do j=1,jj
            do i=1,ii
              if     (xyward .or. pang(i,j).eq.0.0) then
                util1(i,j)=surtx(i,j)
              else
c ---           Rotate from Xward and Yward to Eastward
                util1(i,j)=cos( pang(i,j))*surtx(i,j) +
     &                     sin(-pang(i,j))*surty(i,j)
              endif !pang
            enddo
          enddo
          if     (xyward) then
          call horout(util1,artype,yrflag,time3,iexpt,.true.,
     &                '  surf. x-stress  ',            ! plot name
     &                'surface_x_stress',              ! ncdf name
     &                'surface_downward_xward_stress', ! ncdf standard_name
     &                'Pa',                            ! units
     &                k,ltheta, frmt,ioin)
          else
          call horout(util1,artype,yrflag,time3,iexpt,.true.,
     &                '  surf. e-stress  ',               ! plot name
     &                'surface_e_stress',                 ! ncdf name
     &                'surface_downward_eastward_stress', ! ncdf standard_name
     &                'Pa',                               ! units
     &                k,ltheta, frmt,ioin)
          endif !xyward:else
        endif !txio
c ---   'tyio  ' = surf. y-stress I/O unit (0 no I/O)
        call blkini(ioin,'tyio  ')
        if (ioin.gt.0) then
          do j=1,jj
            do i=1,ii
              if     (xyward .or. pang(i,j).eq.0.0) then
                util1(i,j)=surty(i,j)
              else
c ---           Rotate from Xward and Yward to Northward
                util1(i,j)=cos( pang(i,j))*surty(i,j) -
     &                     sin(-pang(i,j))*surtx(i,j)
              endif !pang
            enddo
          enddo
          if     (xyward) then
          call horout(util1,artype,yrflag,time3,iexpt,.true.,
     &                '  surf. y-stress  ',            ! plot name
     &                'surface_y_stress',              ! ncdf name
     &                'surface_downward_yward_stress', ! ncdf standard_name
     &                'Pa',                            ! units
     &                k,ltheta, frmt,ioin)
          else
          call horout(util1,artype,yrflag,time3,iexpt,.true.,
     &                '  surf. n-stress  ',                ! plot name
     &                'surface_n_stress',                  ! ncdf name
     &                'surface_downward_northward_stress', ! ncdf standard_name
     &                'Pa',                                ! units
     &                k,ltheta, frmt,ioin)
          endif !xyward:else
        endif !tyio
c ---   'curlio' = curl of surf. stress I/O unit (0 no I/O)
c ---   using 2dx and 2dy stencil
        call blkini(ioin,'curlio')
        if (ioin.gt.0) then
          do j=1,jj
c            jm1 = max(j-1, 1)
c            jp1 = min(j+1,jj)
            jm1=j-1
            if(j0.eq.0)jm1=max(j-1,1)
            jp1=j+1
            if(jj+j0.eq.jtdm)jp1=min(j+1,jj)
            do i=1,ii
              if     (i.ne.1) then
                im1 = i-1
              elseif (lperiod) then !i=1
                im1 = ii
              else !i=1 (non-periodic)
                im1 =  1
              endif
              if     (i.ne.ii) then
                ip1 = i+1
              elseif (lperiod) then !i=ii
                ip1 = 1
              else !i=ii (non-periodic)
                ip1 = ii
              endif
              if     (ip(i,j).ne.0 .and.
     &                jm1.ne.j .and. jp1.ne.j .and.
     &                im1.ne.i .and. ip1.ne.i      ) then
                util1(i,j)=(surtx(ip1,j)-surtx(im1,j))/(2.0*scpx(i,j)) -
     &                     (surty(i,jp1)-surty(i,jm1))/(2.0*scpy(i,j))
              else
                util1(i,j)=flag
              endif
            enddo
          enddo
          call horout(util1,artype,yrflag,time3,iexpt,.true.,
     &                'curl surf. stress ',            ! plot name
     &                'surface_stress_curl',           ! ncdf name
     &                ' ',                             ! ncdf standard_name
     &                'Pa/m',                          ! units
     &                k,ltheta, frmt,ioin)
        endif !curlio
      endif !blkini3
c
c --- -----------------
c --- output ice fields
c --- -----------------
c
c --- 'icvio ' = ice coverage I/O unit (0 no I/O)
      call blkini(ioin,'icvio ')
      if (ioin.gt.0) then
        if     (.not.icegln) then
          if(mnproc.eq.1)then
          write(lp,'(a)') 'error - no Sea Ice available'
          call flush(lp)
          endif  
          call xcstop('(archv2data2d - kk)')
          stop
        endif
        call horout(covice,artype,yrflag,time3,iexpt,.true.,
     &              '     ice coverage ',       ! plot name
     &              'ice_coverage',             ! ncdf name
     &              'sea_ice_area_fraction',    ! ncdf standard_name
     &              ' ',                        ! units
     &              k,ltheta, frmt,ioin)
      endif
c
c --- 'ithio ' = ice thickness I/O unit (0 no I/O)
      call blkini(ioin,'ithio ')
      if (ioin.gt.0) then
        if     (.not.icegln) then
          if(mnproc.eq.1)then
          write(lp,'(a)') 'error - no Sea Ice available'
          call flush(lp)
          endif  
          call xcstop('(archv2data2d - kk)')
          stop
        endif
        call horout(thkice,artype,yrflag,time3,iexpt,.true.,
     &              '    ice thickness ',       ! plot name
     &              'ice_thickness',            ! ncdf name
     &              'sea_ice_thickness',        ! ncdf standard_name
     &              'm',                        ! units
     &              k,ltheta, frmt,ioin)
      endif
c
c --- 'ictio ' = ice temperature I/O unit (0 no I/O)
      call blkini(ioin,'ictio ')
      if (ioin.gt.0) then
        if     (.not.icegln) then
          if(mnproc.eq.1)then
          write(lp,'(a)') 'error - no Sea Ice available'
          call flush(lp)
          endif  
          call xcstop('(archv2data2d - kk)')
          stop
        endif
        call horout(temice,artype,yrflag,time3,iexpt,.true.,
     &              '  ice temperature ',       ! plot name
     &              'ice_temperature',          ! ncdf name
     &              ' ',                        ! ncdf standard_name
     &              'degC',                     ! units
     &              k,ltheta, frmt,ioin)
      endif
c
c --- ---------------------
c --- output surface fields
c --- ---------------------
c
c --- 'oetaio' = one+eta          I/O unit (0 no I/O), OPTIONAL, first
c --- 'atthio' = average density  I/O unit (0 no I/O), OPTIONAL
c --- 'bsshio' = bot. pres. SSH   I/O unit (0 no I/O), OPTIONAL, not with atthio
c --- 'dsshio' = non-b.pres SSH   I/O unit (0 no I/O), always after bsshio
c --- 'ssshio' = steric     SSH   I/O unit (0 no I/O), OPTIONAL, not with atthio
c ---                                                  can be after dsshio
c --- 'nsshio' = non-steric SSH   I/O unit (0 no I/O), OPTIONAL, after ssshio
c --- 'montio' = Montgomery Pot.  I/O unit (0 no I/O), OPTIONAL
c --- 'sshio ' = total      SSH   I/O unit (0 no I/O, -ve MKS for NCOM)
      call blkini9(ioin,j, 'atthio','ssshio','montio',
     &                     'bsshio','sshio ','oetaio',  !read 1 of 6
     &                     'XXXXXX','XXXXXX','XXXXXX')
      if (j.eq.6) then !oetaio
        if (ioin.gt.0) then
          if     (.not.loneta) then
            write(lp,'(a)') 'error - no one+eta available'
            call flush(lp)
            stop
          endif
          do j=1,jj
            do i=1,ii
              if (ip(i,j).ne.0) then
                if     (artype.eq.3) then
                  util1(i,j)=onetas(i,j)  !unitless
                else
                  util1(i,j)= oneta(i,j)  !unitless
                endif
              else
                util1(i,j)=flag
              endif
            enddo
          enddo
          k = 0
          call horout(util1, artype,yrflag,time3,iexpt,.true.,
     &                ' one + eta        ',       ! plot name
     &                'one_eta',                  ! ncdf
     &                ' ',                        ! ncdf standard_name
     &                ' ',                        ! units
     &                k,ltheta, frmt,ioin)
        endif !ioin
        call blkini9(ioin,j, 'atthio','ssshio','montio',
     &                       'bsshio','sshio ','XXXXXX',  !read 1 of 5
     &                       'XXXXXX','XXXXXX','XXXXXX')
      endif  !j==5 (oetaio)
      if (j.eq.1) then !atthio
        if (ioin.gt.0) then
          do j=1,jj
            do i=1,ii
              if (ip(i,j).ne.0) then
                dsumth = 0.d0
                dsumdp = 0.d0
                do k=1,kk
                  dsumth = dsumth + dp(i,j,k)*th3d(i,j,k)
                  dsumdp = dsumdp + dp(i,j,k)
                enddo !k
                util1(i,j)=dsumth/dsumdp - thbase
              else
                util1(i,j)=flag
              endif
            enddo
          enddo
          k = 0
          call horout(util1, artype,yrflag,time3,iexpt,.true.,
     &                'average HYCOM th3d',       ! plot name
     &                'avth',                     ! ncdf name
     &                ' ',                        ! ncdf standard_name
     &                'kg/m^3',                   ! units
     &                k,ltheta, frmt,ioin)
        endif
        call blkini(ioin,'sshio ')
        j = 1
      elseif (j.eq.3) then !montio
        if (ioin.gt.0) then
          do j=1,jj
            do i=1,ii
              if (ip(i,j).ne.0) then
                util1(i,j)=montg(i,j)/(thref*9806.0)  !MKS
              else
                util1(i,j)=flag
              endif
            enddo
          enddo
          call horout(util1, artype,yrflag,time3,iexpt,.true.,
     &                ' montgomery pot.  ',       ! plot name
     &                'mont_pot',                 ! ncdf name
     &                ' ',                        ! ncdf standard_name
     &                'm',                        ! units
     &                k,ltheta, frmt,ioin)
        endif
        call blkini(ioin,'sshio ')
        j = 1
      elseif (j.eq.4) then !bsshio
        if (ioin.gt.0) then
          do j=1,jj
            do i=1,ii
              if (ip(i,j).ne.0) then
                util1(i,j)=(srfht(i,j)-montg(i,j))/(thref*9806.0)  !MKS
              else
                util1(i,j)=flag
              endif
            enddo
          enddo
          k = 0
          call horout(util1, artype,yrflag,time3,iexpt,.true.,
     &                ' bot.prs.anom SSH ',       ! plot name
     &                'bp_ssh',                   ! ncdf name
     &                ' ',                        ! ncdf standard_name
     &                'm',                        ! units
     &                k,ltheta, frmt,ioin)
        endif !ioin
        call blkini(ioin,'dsshio')
        if (ioin.gt.0) then
          do j=1,jj
            do i=1,ii
              if (ip(i,j).ne.0) then
                util1(i,j)=montg(i,j)/(thref*9806.0)  !MKS
              else
                util1(i,j)=flag
              endif
            enddo
          enddo
          k = 0
          call horout(util1, artype,yrflag,time3,iexpt,.true.,
     &                ' non-b.p.a  SSH   ',       ! plot name
     &                'nonbp_ssh',                ! ncdf name
     &                ' ',                        ! ncdf standard_name
     &                'm',                        ! units
     &                k,ltheta, frmt,ioin)
        endif !ioin
        call blkini2(ioin,j, 'sshio ','ssshio')  !read one of two
      endif !atthio:montio:bsshio
c --- ssshio?
      if     (j.eq.2) then !ssshio
        if (ioin.gt.0) then
          if     (.not.lsteric) then
            if(mnproc.eq.1)then
            write(lp,'(a)') 'error - no Steric SSH available'
            call flush(lp)
            endif  
            call xcstop('(archv2data2d - kk)')
            stop
          endif
          do j=1,jj
            do i=1,ii
              if (ip(i,j).ne.0) then
                util1(i,j)=steric(i,j)/(thref*9806.0)  !MKS
              else
                util1(i,j)=flag
              endif
            enddo
          enddo
          k = 0
          call horout(util1, artype,yrflag,time3,iexpt,.true.,
     &                ' steric SSH       ',       ! plot name
     &                'steric_ssh',               ! ncdf name
     &                ' ',                        ! ncdf standard_name
     &                'm',                        ! units
     &                k,ltheta, frmt,ioin)
        endif !ioin
        call blkini2(ioin,j, 'nsshio','sshio ')  !read one of two
        if (j.eq.1) then !nsshio
          if (ioin.gt.0) then
            if     (.not.lsteric) then
              if(mnproc.eq.1)then
              write(lp,'(a)') 'error - no Steric SSH available'
              call flush(lp)
              endif  
              call xcstop('(archv2data2d - kk)')
              stop
            endif
            do j=1,jj
              do i=1,ii
                if (ip(i,j).ne.0) then
                  util1(i,j)=(srfht(i,j)-steric(i,j))/(thref*9806.0)  !MKS
                else
                  util1(i,j)=flag
                endif
              enddo
            enddo
            k = 0
            call horout(util1, artype,yrflag,time3,iexpt,.true.,
     &                  ' non-steric SSH   ',       ! plot name
     &                  'nonsteric_ssh',            ! ncdf name
     &                  ' ',                        ! ncdf standard_name
     &                  'm',                        ! units
     &                  k,ltheta, frmt,ioin)
          endif !ioin
          call blkini(ioin,'sshio ')
        endif !nsshio
      endif !ssshio
c --- sshio
      if (ioin.gt.0) then
        do j=1,jj
          do i=1,ii
            if (ip(i,j).ne.0) then
              util1(i,j)=srfht(i,j)/(thref*9806.0)  !MKS
            else
              util1(i,j)=flag
            endif
          enddo
        enddo
        k = 0
        call horout(util1, artype,yrflag,time3,iexpt,.true.,
     &              ' sea surf. height ',       ! plot name
     &              'ssh',                      ! ncdf name (mersea)
     &              'sea_surface_elevation',    ! ncdf standard_name
     &              'm',                        ! units
     &              k,ltheta, frmt,ioin)
      elseif (ioin.lt.0) then
        do j=1,jj
          do i=1,ii
            if (ip(i,j).ne.0) then
              util1(i,j)=srfht(i,j)/(thref*9806.0)  !MKS
            else
              util1(i,j)=0.0
            endif
          enddo
        enddo
        k = 0
        call horout(util1, artype,yrflag,time3,iexpt,.true.,
     &              ' sea surf. height ',       ! plot name
     &              'ssh',                      ! ncdf name (mersea)
     &              'sea_surface_elevation',    ! ncdf standard_name
     &              'm',                        ! units
     &              k,ltheta, frmt,-ioin)
      endif
c
      if     (artype.gt.1) then  ! mean or std. archive
c ---   'bkeio ' = baro. kinetic energy I/O unit (0 no I/O)
        call blkini(ioin,'bkeio ')
        if (ioin.gt.0) then
          k = 0
          call horout(kebaro,artype,yrflag,time3,iexpt,.true.,
     &                '    baro.k.e./mass',                  ! plot name
     &                'barotropic_kinetic_energy_per_mass',  ! ncdf name
     &                ' ',                          ! ncdf standard_name
     &                'm2/s2',                               ! units
     &                k,ltheta, frmt,ioin)
        endif
      endif  ! mean or std. archive
c
c --- 'guvio ' = geosr. u-velocity I/O unit (0 no I/O)
      call blkini3(ioin,i,'guvio ','buvio ','bsfio ')  !read one of three
      if     (i.eq.1) then
        gstruv = .true.
      elseif (i.eq.2) then
        gstruv = .false.
        barouv = .true.
      else
        gstruv = .false.
        barouv = .false.
        iobsf = ioin
      endif
      if     (gstruv) then
        if (ioin.gt.0) then
c ---     use steric SSH for geostrophic currents, if available
          if     (.not.lsteric) then
            steric(:,:) = srfht(:,:)
          endif
          do j=1,jj
c            jm1 = max(j-1, 1)
c            jp1 = min(j+1,jj)
            jm1=j-1
            if(j0.eq.0)jm1=max(j-1,1)
            jp1=j+1
            if(jj+j0.eq.jtdm)jp1=min(j+1,jj)
            do i=1,ii
              if (ip(i,j).ne.0) then  !may be approximate at i=ii and j=jj
                if     (i.ne.1) then
                  im1 = i-1
                elseif (lperiod) then !i=1
                  im1 = ii
                else !i=1 (non-periodic)
                  im1 =  1
                endif
                if     (i.ne.ii) then
                  ip1 = i+1
                elseif (lperiod) then !i=ii
                  ip1 = 1
                else !i=ii (non-periodic)
                  ip1 = ii
                endif
                if     (ip(im1,j).ne.0) then
                  sshw =     steric(im1,j)
                elseif (ip(ip1,j).ne.0) then
                  sshw = 2.0*steric(i,j) - steric(ip1,j)
                else
                  sshw =     steric(i,j)
                endif
                if     (ip(ip1,j).ne.0) then
                  sshe =     steric(ip1,j)
                elseif (ip(im1,j).ne.0) then
                  sshe = 2.0*steric(i,j) - steric(im1,j)
                else
                  sshe =     steric(i,j)
                endif
                if     (ip(i,jm1).ne.0) then
                  sshs =     steric(i,jm1)
                elseif (ip(i,jp1).ne.0) then
                  sshs = 2.0*steric(i,j) - steric(i,jp1)
                else
                  sshs =     steric(i,j)
                endif
                if     (ip(i,jp1).ne.0) then
                  sshn =     steric(i,jp1)
                elseif (ip(i,jm1).ne.0) then
                  sshn = 2.0*steric(i,j) - steric(i,jm1)
                else
                  sshn =     steric(i,j)
                endif
                if     (plat(i,j).ge.0.0) then
                  cori = sin(max( 5.0,plat(i,j))/57.29578)*
     &                   8.d0*1.5707963268/86164.0d0  ! sidereal day
                else
                  cori = sin(min(-5.0,plat(i,j))/57.29578)*
     &                   8.d0*1.5707963268/86164.0d0  ! sidereal day
                endif
                ugi  = -(sshn-sshs)/(cori*2.0*scpy(i,j))
                vgi  =  (sshe-sshw)/(cori*2.0*scpx(i,j))
                if     (xyward .or. pang(i,j).eq.0.0) then
                  util1(i,j)=ugi
                else
c ---             Rotate from Xward and Yward to Eastward
                  util1(i,j)=
     &               cos( pang(i,j))*ugi
     &              +sin(-pang(i,j))*vgi
                endif !pang
              else
                util1(i,j)=flag
              endif
            enddo !i
          enddo !j
          if     (xyward) then
          k = 0
          call horout(util1,artype,yrflag,time3,iexpt,.true.,
     &                '  geostr. x-vel.  ',         ! plot name
     &               'u_geostrophic_velocity',      ! ncdf name
     & 'surface_xward_geostrophic_sea_water_velocity',! ncdf standard_name
     &                'm/s',                        ! units
     &                k,ltheta, frmt,ioin)
          else
          k = 0
          call horout(util1,artype,yrflag,time3,iexpt,.true.,
     &                '   geostr u-vel.  ',         ! plot name
     &               'u_geostrophic_velocity',      ! ncdf name
     & 'surface_eastward_geostrophic_sea_water_velocity',! ncdf standard_name
     &                'm/s',                        ! units
     &                k,ltheta, frmt,ioin)
          endif !xyward:else
        endif !ioin
c
c ---   'gvvio ' = geostr. v-velocity I/O unit (0 no I/O)
        call blkini(ioin,'gvvio ')
        if (ioin.gt.0) then
          do j=1,jj
c            jm1 = max(j-1, 1)
c            jp1 = min(j+1,jj)
            jm1=j-1
            if(j0.eq.0)jm1=max(j-1,1)
            jp1=j+1
            if(jj+j0.eq.jtdm)jp1=min(j+1,jj)
            do i=1,ii
              if (ip(i,j).ne.0) then  !may be approximate at i=ii and j=jj
                if     (i.ne.1) then
                  im1 = i-1
                elseif (lperiod) then !i=1
                  im1 = ii
                else !i=1 (non-periodic)
                  im1 =  1
                endif
                if     (i.ne.ii) then
                  ip1 = i+1
                elseif (lperiod) then !i=ii
                  ip1 = 1
                else !i=ii (non-periodic)
                  ip1 = ii
                endif
                if     (ip(im1,j).ne.0) then
                  sshw =     steric(im1,j)
                elseif (ip(ip1,j).ne.0) then
                  sshw = 2.0*steric(i,j) - steric(ip1,j)
                else
                  sshw =     steric(i,j)
                endif
                if     (ip(ip1,j).ne.0) then
                  sshe =     steric(ip1,j)
                elseif (ip(im1,j).ne.0) then
                  sshe = 2.0*steric(i,j) - steric(im1,j)
                else
                  sshe =     steric(i,j)
                endif
                if     (ip(i,jm1).ne.0) then
                  sshs =     steric(i,jm1)
                elseif (ip(i,jp1).ne.0) then
                  sshs = 2.0*steric(i,j) - steric(i,jp1)
                else
                  sshs =     steric(i,j)
                endif
                if     (ip(i,jp1).ne.0) then
                  sshn =     steric(i,jp1)
                elseif (ip(i,jm1).ne.0) then
                  sshn = 2.0*steric(i,j) - steric(i,jm1)
                else
                  sshn =     steric(i,j)
                endif
                if     (plat(i,j).ge.0.0) then
                  cori = sin(max( 5.0,plat(i,j))/57.29578)*
     &                   8.d0*1.5707963268/86164.0d0  ! sidereal day
                else
                  cori = sin(min(-5.0,plat(i,j))/57.29578)*
     &                   8.d0*1.5707963268/86164.0d0  ! sidereal day
                endif
                ugi  = -(sshn-sshs)/(cori*2.0*scpy(i,j))
                vgi  =  (sshe-sshw)/(cori*2.0*scpx(i,j))
                if     (xyward .or. pang(i,j).eq.0.0) then
                  util1(i,j)=vgi
                else
c ---             Rotate from Xward and Yward to Northward
                  util1(i,j)=
     &               cos( pang(i,j))*vgi
     &              -sin(-pang(i,j))*ugi
                endif !pang
              else
                util1(i,j)=flag
              endif
            enddo !i
          enddo !j
          if     (xyward) then
          k = 0
          call horout(util1,artype,yrflag,time3,iexpt,.true.,
     &                '  geostr. y-vel.  ',         ! plot name
     &               'v_geostrophic_velocity',      ! ncdf name
     & 'surface_yward_geostrophic_sea_water_velocity',! ncdf standard_name
     &                'm/s',                        ! units
     &                k,ltheta, frmt,ioin)
          else
          k = 0
          call horout(util1,artype,yrflag,time3,iexpt,.true.,
     &                '  geostr. v-vel.  ',         ! plot name
     &               'v_geostrophic_velocity',      ! ncdf name
     & 'surface_northward_geostrophic_sea_water_velocity',! ncdf standard_name
     &                'm/s',                        ! units
     &                k,ltheta, frmt,ioin)
          endif !xyward:else
        endif !ioin
c
c ---   'gspio ' = geostr. speed I/O unit (0 no I/O)
        call blkini(ioin, 'gspio ')
        if (ioin.gt.0) then
          do j=1,jj
c            jm1 = max(j-1, 1)
c            jp1 = min(j+1,jj)
            jm1=j-1
            if(j0.eq.0)jm1=max(j-1,1)
            jp1=j+1
            if(jj+j0.eq.jtdm)jp1=min(j+1,jj)
            do i=1,ii
              if (ip(i,j).ne.0) then  !may be approximate at i=ii and j=jj
                if     (i.ne.1) then
                  im1 = i-1
                elseif (lperiod) then !i=1
                  im1 = ii
                else !i=1 (non-periodic)
                  im1 =  1
                endif
                if     (i.ne.ii) then
                  ip1 = i+1
                elseif (lperiod) then !i=ii
                  ip1 = 1
                else !i=ii (non-periodic)
                  ip1 = ii
                endif
                if     (ip(im1,j).ne.0) then
                  sshw =     steric(im1,j)
                elseif (ip(ip1,j).ne.0) then
                  sshw = 2.0*steric(i,j) - steric(ip1,j)
                else
                  sshw =     steric(i,j)
                endif
                if     (ip(ip1,j).ne.0) then
                  sshe =     steric(ip1,j)
                elseif (ip(im1,j).ne.0) then
                  sshe = 2.0*steric(i,j) - steric(im1,j)
                else
                  sshe =     steric(i,j)
                endif
                if     (ip(i,jm1).ne.0) then
                  sshs =     steric(i,jm1)
                elseif (ip(i,jp1).ne.0) then
                  sshs = 2.0*steric(i,j) - steric(i,jp1)
                else
                  sshs =     steric(i,j)
                endif
                if     (ip(i,jp1).ne.0) then
                  sshn =     steric(i,jp1)
                elseif (ip(i,jm1).ne.0) then
                  sshn = 2.0*steric(i,j) - steric(i,jm1)
                else
                  sshn =     steric(i,j)
                endif
                if     (plat(i,j).ge.0.0) then
                  cori = sin(max( 5.0,plat(i,j))/57.29578)*
     &                   8.d0*1.5707963268/86164.0d0  ! sidereal day
                else
                  cori = sin(min(-5.0,plat(i,j))/57.29578)*
     &                   8.d0*1.5707963268/86164.0d0  ! sidereal day
                endif
                ugi  = -(sshn-sshs)/(cori*2.0*scpy(i,j))
                vgi  =  (sshe-sshw)/(cori*2.0*scpx(i,j))
                util1(i,j)=sqrt( ugi**2 + vgi**2 )
              else
                util1(i,j)=flag
              endif
            enddo !i
          enddo !j
          k = 0
          call horout(util1, artype,yrflag,time3,iexpt,.true.,
     &               'geostrophic speed ',        ! plot name
     &               'geostrophic_speed',         ! ncdf name
     &    'surface_geostrophic_sea_water_speed',  ! ncdf standard_name
     &                'm/s',                      ! units
     &                k,ltheta, frmt,ioin)
        endif !ioin
c
        call blkini2(ioin,i,'buvio ','bsfio ')
        if     (i.eq.1) then
          barouv = .true.
        else
          barouv = .false.
          iobsf = ioin
        endif
      endif !gstruv
c
c --- 'buvio ' = baro.  u-velocity I/O unit (0 no I/O, -ve on u-grid for NCOM)
      if     (barouv) then
        if (ioin.gt.0) then
          do j=1,jj
c            jm1 = max(j-1, 1)
c            jp1 = min(j+1,jj)
            jm1=j-1
            if(j0.eq.0)jm1=max(j-1,1)
            jp1=j+1
            if(jj+j0.eq.jtdm)jp1=min(j+1,jj)
            do i=1,ii
              if (ip(i,j).ne.0) then  !may be approximate at i=ii and j=jj
                if     (i.ne.1) then
                  im1 = i-1
                elseif (lperiod) then !i=1
                  im1 = ii
                else !i=1 (non-periodic)
                  im1 =  1
                endif
                if     (i.ne.ii) then
                  ip1 = i+1
                elseif (lperiod) then !i=ii
                  ip1 = 1
                else !i=ii (non-periodic)
                  ip1 = ii
                endif
c ---           Make transport accurate at the expense of velocity
                if     (xyward .or. pang(i,j).eq.0.0) then
                  util1(i,j)=(ubaro(i,  j)*min(p(i,  j,kk+1),
     &                                         p(im1,j,kk+1) )+
     &                        ubaro(ip1,j)*min(p(ip1,j,kk+1),
     &                                         p(i,  j,kk+1) ) )
     &                       /max(2.0*p(i,j,kk+1),onemm)
                else
c ---             Rotate from Xward and Yward to Eastward
                  util1(i,j)=
     &              (cos( pang(i,j))*(ubaro(i,  j)*min(p(i,  j,kk+1),
     &                                                 p(im1,j,kk+1) )+
     &                                ubaro(ip1,j)*min(p(ip1,j,kk+1),
     &                                                 p(i,  j,kk+1) ) )
     &              +sin(-pang(i,j))*(vbaro(i,j  )*min(p(i,j  ,kk+1),
     &                                                 p(i,jm1,kk+1) )+
     &                                vbaro(i,jp1)*min(p(i,jp1,kk+1),
     &                                                 p(i,j  ,kk+1) ) )
     &              )/max(2.0*p(i,j,kk+1),onemm)
                endif !pang
              else
                util1(i,j)=flag
              endif
            enddo !i
          enddo !j
          if     (xyward) then
          k = 0
          call horout(util1,artype,yrflag,time3,iexpt,.true.,
     &                '    baro. x-vel.  ',         ! plot name
     &                'u_barotropic_velocity',      ! ncdf name
     &            'barotropic_sea_water_x_velocity',! ncdf standard_name
     &                'm/s',                        ! units
     &                k,ltheta, frmt,ioin)
          else
          k = 0
          call horout(util1,artype,yrflag,time3,iexpt,.true.,
     &                '    baro. u-vel.  ',         ! plot name
     &                'u_barotropic_velocity',      ! ncdf name
     &     'barotropic_eastward_sea_water_velocity',! ncdf standard_name
     &                'm/s',                        ! units
     &                k,ltheta, frmt,ioin)
          endif !xyward:else
        elseif (ioin.lt.0) then
          do j=1,jj
            do i=1,ii
              if (iu(i,j).ne.0) then
                util1(i,j)=ubaro(i,j)
              else
                util1(i,j)=0.0
              endif
            enddo
          enddo
          k = 0
          call horout(util1,artype,yrflag,time3,iexpt,.true.,
     &                '    baro. u-vel.  ',         ! plot name
     &                'u_barotropic_velocity',      ! ncdf name
     &                ' ',                          ! ncdf standard_name
     &                'm/s',                        ! units
     &                k,ltheta, frmt,-ioin)
        endif !ioin
c
c ---   'bvvio ' = baro. v-velocity I/O unit (0 no I/O, -ve on v-grid for NCOM)
        call blkini(ioin,'bvvio ')
        if (ioin.gt.0) then
          do j=1,jj
c            jm1 = max(j-1, 1)
c            jp1 = min(j+1,jj)
            jm1=j-1
            if(j0.eq.0)jm1=max(j-1,1)
            jp1=j+1
            if(jj+j0.eq.jtdm)jp1=min(j+1,jj)
            do i=1,ii
              if (ip(i,j).ne.0) then  !may be approximate at i=ii and j=jj
                if     (xyward .or. pang(i,j).eq.0.0) then
                  util1(i,j)=(vbaro(i,j  )*min(p(i,j  ,kk+1),
     &                                         p(i,jm1,kk+1) )+
     &                        vbaro(i,jp1)*min(p(i,jp1,kk+1),
     &                                         p(i,j  ,kk+1) ) )
     &                       /max(2.0*p(i,j,kk+1),onemm)
                else
c ---             Rotate from Xward and Yward to Eastward
                  if     (i.ne.ii) then
                    ip1 = i+1
                  elseif (lperiod) then !i=ii
                    ip1 = 1
                  else !i=ii (non-periodic)
                    ip1 = ii
                  endif
                  if     (i.ne.1) then
                    im1 = i-1
                  elseif (lperiod) then !i=1
                    im1 = ii
                  else !i=1 (non-periodic)
                    im1 =  1
                  endif
                  util1(i,j)=
     &              (cos( pang(i,j))*(vbaro(i,j  )*min(p(i,j  ,kk+1),
     &                                                 p(i,jm1,kk+1) )+
     &                                vbaro(i,jp1)*min(p(i,jp1,kk+1),
     &                                                 p(i,j  ,kk+1) ) )
     &              -sin(-pang(i,j))*(ubaro(i,  j)*min(p(i,  j,kk+1),
     &                                                 p(im1,j,kk+1) )+
     &                                ubaro(ip1,j)*min(p(ip1,j,kk+1),
     &                                                 p(i,  j,kk+1) ) )
     &              )/max(2.0*p(i,j,kk+1),onemm)
                endif !pang
              else
                util1(i,j)=flag
              endif
            enddo
          enddo
          if     (xyward) then
          k = 0
          call horout(util1,artype,yrflag,time3,iexpt,.true.,
     &                '    baro. y-vel.  ',         ! plot name
     &                'v_barotropic_velocity',      ! ncdf name
     &            'barotropic_sea_water_y_velocity',! ncdf standard_name
     &                'm/s',                        ! units
     &                k,ltheta, frmt,ioin)
          else
          k = 0
          call horout(util1,artype,yrflag,time3,iexpt,.true.,
     &                '    baro. v-vel.  ',         ! plot name
     &                'v_barotropic_velocity',      ! ncdf name
     &    'barotropic_northward_sea_water_velocity',! ncdf standard_name
     &                'm/s',                        ! units
     &                k,ltheta, frmt,ioin)
          endif !xyward:else
        elseif (ioin.lt.0) then
          do j=1,jj
            do i=1,ii
              if (iv(i,j).ne.0) then
                util1(i,j)=vbaro(i,j)
              else
                util1(i,j)=0.0
              endif
            enddo
          enddo
          k = 0
          call horout(util1,artype,yrflag,time3,iexpt,.true.,
     &                '    baro. v-vel.  ',         ! plot name
     &                'v_barotropic_velocity',      ! ncdf name
     &                ' ',                          ! ncdf standard_name
     &                'm/s',                        ! units
     &                k,ltheta, frmt,-ioin)
        endif !ioin
c
c ---   'bspio ' = baro. speed I/O unit (0 no I/O)
        call blkini(ioin, 'bspio ')
        if (ioin.gt.0) then
          do j=1,jj
c            jm1 = max(j-1, 1)
c            jp1 = min(j+1,jj)
            jm1=j-1
            if(j0.eq.0)jm1=max(j-1,1)
            jp1=j+1
            if(jj+j0.eq.jtdm)jp1=min(j+1,jj)
            do i=1,ii
              if (ip(i,j).ne.0) then  !may be approximate at i=ii and j=jj
                if     (i.ne.1) then
                  im1 = i-1
                elseif (lperiod) then !i=1
                  im1 = ii
                else !i=1 (non-periodic)
                  im1 =  1
                endif
                if     (i.ne.ii) then
                  ip1 = i+1
                elseif (lperiod) then !i=ii
                  ip1 = 1
                else !i=ii (non-periodic)
                  ip1 = ii
                endif
c ---           Make transport accurate at the expense of velocity
                ubi=(ubaro(i,  j)*min(p(i,  j,kk+1),
     &                                p(im1,j,kk+1) )+
     &               ubaro(ip1,j)*min(p(ip1,j,kk+1),
     &                                p(i,  j,kk+1) ) )
     &              /max(2.0*p(i,j,kk+1),onemm)
                vbi=(vbaro(i,j  )*min(p(i,j  ,kk+1),
     &                                p(i,jm1,kk+1) )+
     &               vbaro(i,jp1)*min(p(i,jp1,kk+1),
     &                                p(i,j  ,kk+1) ) )
     &              /max(2.0*p(i,j,kk+1),onemm)
                util1(i,j)=sqrt( ubi**2 + vbi**2 )
              else
                util1(i,j)=flag
              endif
            enddo
          enddo
          k = 0
          call horout(util1, artype,yrflag,time3,iexpt,.true.,
     &                'barotropic speed ',        ! plot name
     &                'barotropic_speed',         ! ncdf name
     &     'barotropic_sea_water_speed',          ! ncdf standard_name
     &                'm/s',                      ! units
     &                k,ltheta, frmt,ioin)
        endif !ioin
c
      endif !barouv
c
c --- 'bsfio ' = barotropic strmfn. I/O unit (0 no I/O)
      if     (barouv) then
        call blkini(ioin,'bsfio ')
      else
        ioin = iobsf
      endif !barouv:else
      if (ioin.gt.0) then
c
c ---   calculate strmf on the Q grid, mostly using vbaro.
c ---   note that strmf(ii,jj) is always 0.0, but streamfunctions
c ---    are arbitraty w.r.t. a constant offset anyway
        do j=jj,1,-1
c          jm1 = max(1,j-1)  !not accurate for j=1, when using a subregion.
          jm1=j-1
          if(j0.eq.0)jm1=max(j-1,1)
          i   = ii
          im1 = ii-1
          if (iu(i,j).eq.1) then
            ubi = ubaro(i,j)*min(p(i,j,kk+1),p(im1,j,kk+1))
          else
            ubi = 0.0
          endif
          if     (j.eq.jj) then
            strmfu = 0.0
          else
            strmfu = strmfu + ubi*scux(i,j)
          endif
          strmf(i,j) = strmfu
          strmft     = strmfu
          do i=ii-1,1,-1
            if (iv(i,j).eq.1) then
              vbi = vbaro(i,j)*min(p(i,j,kk+1),p(i,jm1,kk+1))
            else
              vbi = 0.0
            endif
            strmft = strmft - vbi*scvy(i,j)
            strmf(i,j) = strmft  !m^3/s
          enddo !i
        enddo !j
c
c ---   interpolate from q-grid to p-grid.
        do j=1,jj
          do i=1,ii
            if (ip(i,j).ne.0 .and. i.lt.ii .and. j.lt.jj) then
              s1=0.0
              s2=0.0
              if     (iq(i,  j)  .eq.1) then
                s1 = s1 + 1.0
                s2 = s2 + strmf(i,  j)
              endif
              if     (iq(i+1,j)  .eq.1) then
                s1 = s1 + 1.0
                s2 = s2 + strmf(i+1,j)
              endif
              if     (iq(i,  j+1).eq.1) then
                s1 = s1 + 1.0
                s2 = s2 + strmf(i,  j+1)
              endif
              if     (iq(i+1,j+1).eq.1) then
                s1 = s1 + 1.0
                s2 = s2 + strmf(i+1,j+1)
              endif
              if     (s1.ne.0.0) then
                util1(i,j)=s2/s1
              else
                util1(i,j)=flag
              endif
            else
              util1(i,j)=flag
            endif
          enddo
        enddo
c
ccc     call zebra(util1,ii,ii,jj)
ccc     write (*,'('' shown above: barotropic stream function'')')
c
        k = 0
        call horout(util1, artype,yrflag,time3,iexpt,.true.,
     &              'barotr. strmf. (V)',              ! plot name
     &              'bfsd_v',                          ! ncdf name
     &              'ocean_barotropic_streamfunction', ! ncdf standard_name
     &              'm3/s',                            ! units
     &              k,ltheta, frmt,ioin)
c
c ---   calculate strmf on the Q grid, mostly using ubaro.
c ---   note that strmf(ii,jj) is always 0.0, but streamfunctions
c ---    are arbitraty w.r.t. a constant offset anyway
        do i=ii,1,-1
          im1 = max(1,i-1) !not accurate for i=1, when using a subregion.
          j   = jj
          jm1 = jj-1
          if (iv(i,j).eq.1) then
            vbi = vbaro(i,j)*min(p(i,j,kk+1),p(i,jm1,kk+1))
          else
            vbi = 0.0
          endif
          if     (i.eq.ii) then
            strmfv = 0.0
          else
            strmfv = strmfv - vbi*scvy(i,j)
          endif
          strmf(i,j) = strmfv
          strmft     = strmfv
          do j=jj-1,1,-1
            if (iu(i,j).eq.1) then
              ubi = ubaro(i,j)*min(p(i,j,kk+1),p(im1,j,kk+1))
            else
              ubi = 0.0
            endif
            strmft     = strmft + ubi*scux(i,j)
            strmf(i,j) = strmft  !m^3/s
          enddo !i
        enddo !j
c
c ---   interpolate from q-grid to p-grid.
        do j=1,jj
          do i=1,ii
            if (ip(i,j).ne.0 .and. i.lt.ii .and. j.lt.jj) then
              s1=0.0
              s2=0.0
              if     (iq(i,  j)  .eq.1) then
                s1 = s1 + 1.0
                s2 = s2 + strmf(i,  j)
              endif
              if     (iq(i+1,j)  .eq.1) then
                s1 = s1 + 1.0
                s2 = s2 + strmf(i+1,j)
              endif
              if     (iq(i,  j+1).eq.1) then
                s1 = s1 + 1.0
                s2 = s2 + strmf(i,  j+1)
              endif
              if     (iq(i+1,j+1).eq.1) then
                s1 = s1 + 1.0
                s2 = s2 + strmf(i+1,j+1)
              endif
              if     (s1.ne.0.0) then
                util1(i,j)=s2/s1
              else
                util1(i,j)=flag
              endif
            else
              util1(i,j)=flag
            endif
          enddo
        enddo
c
ccc     call zebra(util1,ii,ii,jj)
ccc     write (*,'('' shown above: barotropic stream function'')')
c
        k = 0
        call horout(util1, artype,yrflag,time3,iexpt,.true.,
     &              'barotr. strmf. (U)',              ! plot name
     &              'bfsd_u',                          ! ncdf name
     &              'ocean_barotropic_streamfunction', ! ncdf standard_name
     &              'm3/s',                            ! units
     &              k,ltheta, frmt,ioin)
      endif
c
c --- -------------------------
c --- output mixed layer fields
c --- -------------------------
c
c --- 'uvmio ' = mixed layer u-velocity I/O unit (0 no I/O)
      call blkini(ioin, 'uvmio ')
      if (ioin.gt.0) then
        do j=1,jj
          do i=1,ii
            if (ip(i,j).ne.0 .and. i.lt.ii) then
              util1(i,j)=0.5*(umix(i,j)+umix(i+1,j))
            else
              util1(i,j)=flag
            endif
          enddo
        enddo
        k = 0
        call horout(util1, artype,yrflag,time3,iexpt,.true.,
     &              'mix.l. u-velocity ',       ! plot name
     &              'mixed_layer_u_velocity',   ! ncdf name
     &              ' ',                        ! ncdf standard_name
     &              'm/s',                      ! units
     &              k,ltheta, frmt,ioin)
      endif
c
c --- 'vvmio ' = mixed layer v-velocity I/O unit (0 no I/O)
      call blkini(ioin, 'vvmio ')
      if (ioin.gt.0) then
        do j=1,jj
          do i=1,ii
            if (ip(i,j).ne.0 .and. j.lt.jj) then
              util1(i,j)=0.5*(vmix(i,j)+vmix(i,j+1))
            else
              util1(i,j)=flag
            endif
          enddo
        enddo
        k = 0
        call horout(util1, artype,yrflag,time3,iexpt,.true.,
     &              ' mixl. v-velocity ',       ! plot name
     &              'mixed_layer_v_velocity',   ! ncdf name
     &              ' ',                        ! ncdf standard_name
     &              'm/s',                      ! units
     &              k,ltheta, frmt,ioin)
      endif
c
c --- 'spmio ' = mixed layer speed I/O unit (0 no I/O)
      call blkini(ioin, 'spmio ')
      if (ioin.gt.0) then
        do j=1,jj
          do i=1,ii
            if (ip(i,j).ne.0 .and. i.lt.ii .and. j.lt.jj) then
              util1(i,j)=0.5*sqrt( (umix(i,j)+umix(i+1,j))**2 +
     &                             (vmix(i,j)+vmix(i,j+1))**2  )
            else
              util1(i,j)=flag
            endif
          enddo
        enddo
        k = 0
        call horout(util1, artype,yrflag,time3,iexpt,.true.,
     &              'mixed-layer speed ',       ! plot name
     &              'mixed_layer_speed',        ! ncdf name
     &              ' ',                        ! ncdf standard_name
     &              'm/s',                      ! units
     &              k,ltheta, frmt,ioin)
      endif
c
c --- 'bltio ' = bnd. lay. thick. I/O unit (0 no I/O)
      call blkini(ioin,'bltio ')
      if (ioin.gt.0) then
        k = 0
        call horout(dpbl,  artype,yrflag,time3,iexpt,.true.,
     &              'bnd.layr.thickness',                ! plot name
     &              'surface_boundary_layer_thickness',  ! ncdf name
     &              ' ',                        ! ncdf standard_name
     &              'm',                                 ! units
     &              k,ltheta, frmt,ioin)
      endif
c
c --- 'mltio ' = mix. lay. thick. I/O unit (0 no I/O)
      call blkini(ioin,'mltio ')
      if (ioin.gt.0) then
        k = 0
        call horout(dpmixl,artype,yrflag,time3,iexpt,.true.,
     &              'mix.layr.thickness',          ! plot name
     &              'mixed_layer_thickness',       ! ncdf name
     &              'ocean_mixed_layer_thickness', ! ncdf standard_name
     &              'm',                           ! units
     &              k,ltheta, frmt,ioin)
      endif
c
c --- 'sstio ' = mix. lay. temp.  I/O unit (0 no I/O)
      call blkini(ioin,'sstio ')
      if (ioin.gt.0) then
        k = 0
        call horout(tmix,  artype,yrflag,time3,iexpt,.true.,
     &              'mix.layr.temp     ',       ! plot name
     &              'mixed_layer_temperature',  ! ncdf name
     &              ' ',                        ! ncdf standard_name
     &              'degC',                     ! units
     &              k,ltheta, frmt,ioin)
      endif
c
c --- 'sssio ' = mix. lay. saln.  I/O unit (0 no I/O)
      call blkini(ioin,'sssio ')
      if (ioin.gt.0) then
        k = 0
        call horout(smix,  artype,yrflag,time3,iexpt,.true.,
     &              'mix.layr.saln     ',       ! plot name
     &              'mixed_layer_salinity',     ! ncdf name
     &              ' ',                        ! ncdf standard_name
     &              'psu',                      ! units
     &              k,ltheta, frmt,ioin)
      endif
c
c --- 'ssdio ' = mix. lay. dens.  I/O unit (0 no I/O)
      call blkini(ioin,'ssdio ')
      if (ioin.gt.0) then
        k = 0
        call horout(thmix, artype,yrflag,time3,iexpt,.true.,
     &              'mix.layr.dens     ',       ! plot name
     &              'mixed_layer_density',      ! ncdf name
     &              ' ',                        ! ncdf standard_name
     &              'sigma',                    ! units
     &              k,ltheta, frmt,ioin)
      endif
c
      if     (artype.gt.1) then  ! mean or std. archive
c ---   'mkeio ' = m.l. kinetic energy I/O unit (0 no I/O)
        call blkini(ioin,'mkeio ')
        if (ioin.gt.0) then
          k = 0
          call horout(kemix, artype,yrflag,time3,iexpt,.true.,
     &                '    mixl.k.e./mass',                   ! plot name
     &                'mixed_layer_kinetic_energy_per_mass',  ! ncdf name
     &                ' ',                           ! ncdf standard_name
     &                'm2/s2',                                ! units
     &                k,ltheta, frmt,ioin)
        endif
      endif  ! mean or std. archive
c
c --- ----------------------
c --- output selected layers
c --- ----------------------
c
      do  !layer loop
c
c ---   'kf    ' = first output layer (=0 end output; <0 label with layer #)
c ---   'kl    ' = last  output layer
        call blkini(kin,'kf    ')
        ltheta = kin.gt.0
        kf     = abs(kin)
        if     (kf.eq.0) then
          exit
        endif
        call blkini(kl, 'kl    ')
        kl = abs(kl)
        if     (kl.gt.kkin) then
          write(lp,'(a)') 'error - kl larger than kdm'
          exit
        elseif (kl.lt.kf) then
          write(lp,'(a)') 'error - kl smaller than kf'
          exit
        endif
c
c ---   -------------------
c ---   output layer velocity
c ---   -------------------
c
c ---   'uvlio ' = u-velocity I/O unit (0 no I/O, -ve on u-grid for NCOM)
        call blkini(ioin, 'uvlio ')
        if (ioin.gt.0) then
          do k= kf,kl
          do j=1,jj
c            jp1 = min(j+1,jj)
            jp1=j+1
            if(jj+j0.eq.jtdm)jp1=min(j+1,jj)
            do i=1,ii
              if (ip(i,j).ne.0) then  !may be approximate at i=ii and j=jj
                if     (i.ne.ii) then
                  ip1 = i+1
                elseif (lperiod) then !i=ii
                  ip1 = 1
                else !i=ii (non-periodic)
                  ip1 = ii
                endif
                if     (xyward .or. pang(i,j).eq.0.0) then
                  utilk(i,j,k)=     0.5*(u(i,j,k)+u(ip1,j,k))
                else
c ---             Rotate from Xward and Yward to Eastward
                  utilk(i,j,k)=
     &              cos( pang(i,j))*0.5*(u(i,j,k)+u(ip1,j,k)) +
     &              sin(-pang(i,j))*0.5*(v(i,j,k)+v(i,jp1,k))
                endif !pang
                if (mthin) then
                  if (p(i,j,k)+onecm.gt.p(i,j,k+1)) then
                    utilk(i,j,k)=flag
                  endif
                else
                  if (p(i,j,k)+onecm.gt.p(i,j,kk+1)) then
                    utilk(i,j,k)=flag
                  endif
                endif
              else
                utilk(i,j,k)=flag
              endif
            enddo
          enddo
          enddo
          if     (baclin) then
          if     (xyward) then
          k = 0
          call horout_3d(utilk, artype,yrflag,time3,iexpt,.true.,
     &                ' x-b.vel.',                   ! plot name
     &                'u_baroclinic_velocity',       ! ncdf name
     &            'baroclinic_sea_water_x_velocity', ! ncdf standard_name
     &                'm/s',                         ! units
     &                kf,kl,ltheta, frmt,ioin)
          else
          k = 0
          call horout_3d(utilk, artype,yrflag,time3,iexpt,.true.,
     &                ' u-b.vel.',                   ! plot name
     &                'u_baroclinic_velocity',       ! ncdf name
     &     'baroclinic_eastward_sea_water_velocity', ! ncdf standard_name
     &                'm/s',                         ! units
     &                kf,kl,ltheta, frmt,ioin)
          endif !xyward:else
          else !total velocity
          if     (xyward) then
          k = 0
          call horout_3d(utilk, artype,yrflag,time3,iexpt,.true.,
     &                ' x-veloc.',                   ! plot name
     &                'u_velocity',                  ! ncdf name
     &                'sea_water_x_velocity',        ! ncdf standard_name
     &                'm/s',                         ! units
     &                kf,kl,ltheta, frmt,ioin)
          else
          k = 0
          call horout_3d(utilk, artype,yrflag,time3,iexpt,.true.,
     &                ' u-veloc.',                   ! plot name
     &                'u_velocity',                  ! ncdf name
     &                'eastward_sea_water_velocity', ! ncdf standard_name
     &                'm/s',                         ! units
     &                kf,kl,ltheta, frmt,ioin)
          endif !xyward:else
          endif !baclin:else
        elseif (ioin.lt.0) then
          do k= kf,kl
          do j=1,jj
            do i=1,ii
              if (iu(i,j).ne.0) then
                utilk(i,j,k)=u(i,j,k)
              else
                utilk(i,j,k)=0.0
              endif
            enddo
          enddo
          enddo
          if     (baclin) then
          call horout_3d(utilk, artype,yrflag,time3,iexpt,.true.,
     &                ' u-b.vel.',                   ! plot name
     &                'u_baroclinic_velocity',       ! ncdf name
     &                ' ',                           ! ncdf standard_name
     &                'm/s',                         ! units
     &                kf,kl,ltheta, frmt,-ioin)
          else !total velocity
          call horout_3d(utilk, artype,yrflag,time3,iexpt,.true.,
     &                ' u-veloc.',                   ! plot name
     &                'u_velocity',                  ! ncdf name
     &                ' ',                           ! ncdf standard_name
     &                'm/s',                         ! units
     &                kf,kl,ltheta, frmt,-ioin)
          endif
        endif
c
c ---   'vvlio ' = v-velocity I/O unit (0 no I/O, -ve on v-grid for NCOM)
        call blkini(ioin, 'vvlio ')
        if (ioin.gt.0) then
          do k= kf,kl
          do j=1,jj
c           jp1 = min(j+1,jj)
            jp1=j+1
            if(jj+j0.eq.jtdm)jp1=min(j+1,jj)
            do i=1,ii
              if (ip(i,j).ne.0) then  !may be approximate at i=ii and j=jj
                if     (xyward .or. pang(i,j).eq.0.0) then
                  utilk(i,j,k)=     0.5*(v(i,j,k)+v(i,jp1,k))
                else
c ---             Rotate from Xward and Yward to Northward
                  if     (i.ne.ii) then
                    ip1 = i+1
                  elseif (lperiod) then !i=ii
                    ip1 = 1
                  else !i=ii (non-periodic)
                    ip1 = ii
                  endif
                  utilk(i,j,k)=
     &              cos( pang(i,j))*0.5*(v(i,j,k)+v(i,jp1,k)) -
     &              sin(-pang(i,j))*0.5*(u(i,j,k)+u(ip1,j,k))
                endif !pang
                if (mthin) then
                  if (p(i,j,k)+onecm.gt.p(i,j,k+1)) then
                    utilk(i,j,k)=flag
                  endif
                else
                  if (p(i,j,k)+onecm.gt.p(i,j,kk+1)) then
                    utilk(i,j,k)=flag
                  endif
                endif
              else
                utilk(i,j,k)=flag
              endif
            enddo
          enddo
          enddo
          if     (baclin) then
          if     (xyward) then
          call horout_3d(utilk, artype,yrflag,time3,iexpt,.true.,
     &                ' y-b.vel.',                    ! plot name
     &                'v_baroclinic_velocity',        ! ncdf name
     &             'baroclinic_sea_water_y_velocity', ! ncdf standard_name
     &                'm/s',                          ! units
     &                kf,kl,ltheta, frmt,ioin)
          else
          call horout_3d(utilk, artype,yrflag,time3,iexpt,.true.,
     &                ' v-b.vel.',                    ! plot name
     &                'v_baroclinic_velocity',        ! ncdf name
     &     'baroclinic_northward_sea_water_velocity', ! ncdf standard_name
     &                'm/s',                          ! units
     &                kf,kl,ltheta, frmt,ioin)
          endif !xyward:else
          else !total velocity
          if     (xyward) then
          call horout_3d(utilk, artype,yrflag,time3,iexpt,.true.,
     &                ' y-veloc.',                    ! plot name
     &                'v_velocity',                   ! ncdf name
     &                'sea_water_y_velocity',         ! ncdf standard_name
     &                'm/s',                          ! units
     &                kf,kl,ltheta, frmt,ioin)
          else
          call horout_3d(utilk, artype,yrflag,time3,iexpt,.true.,
     &                ' v-veloc.',                    ! plot name
     &                'v_velocity',                   ! ncdf name
     &                'northward_sea_water_velocity', ! ncdf standard_name
     &                'm/s',                          ! units
     &                kf,kl,ltheta, frmt,ioin)
          endif !xyward:else
          endif !baclin:else
        elseif (ioin.lt.0) then
          do k= kf,kl
          do j=1,jj
            do i=1,ii
              if (iv(i,j).ne.0) then
                utilk(i,j,k)=v(i,j,k)
              else
                utilk(i,j,k)=0.0
              endif
            enddo
          enddo
          enddo
          if     (baclin) then
          call horout_3d(utilk, artype,yrflag,time3,iexpt,.true.,
     &                ' v-b.vel.',                ! plot name
     &                'v_baroclinic_velocity',    ! ncdf name
     &                ' ',                        ! ncdf standard_name
     &                'm/s',                      ! units
     &                kf,kl,ltheta, frmt,-ioin)
          else !total velocity
          call horout_3d(utilk, artype,yrflag,time3,iexpt,.true.,
     &                ' v-veloc.',                ! plot name
     &                'v_velocity',               ! ncdf name
     &                ' ',                        ! ncdf standard_name
     &                'm/s',                      ! units
     &                kf,kl,ltheta, frmt,-ioin)
          endif
        endif
c
c ---   'splio ' = speed I/O unit (0 no I/O)
        call blkini(ioin, 'splio ')
        if (ioin.gt.0) then
          do k= kf,kl
          do j=1,jj
c            jp1 = min(j+1,jj)
            jp1=j+1
            if(jj+j0.eq.jtdm)jp1=min(j+1,jj)
            do i=1,ii
              if (ip(i,j).ne.0) then  !may be approximate at i=ii and j=jj
                if     (i.ne.ii) then
                  ip1 = i+1
                elseif (lperiod) then !i=ii
                  ip1 = 1
                else !i=ii (non-periodic)
                  ip1 = ii
                endif
                utilk(i,j,k)=0.5*sqrt( (u(i,j,k)+u(ip1,j,k))**2 +
     &                                 (v(i,j,k)+v(i,jp1,k))**2  )
                if (mthin) then
                  if (p(i,j,k)+onecm.gt.p(i,j,k+1)) then
                    utilk(i,j,k)=flag
                  endif
                else
                  if (p(i,j,k)+onecm.gt.p(i,j,kk+1)) then
                    utilk(i,j,k)=flag
                  endif
                endif
              else
                utilk(i,j,k)=flag
              endif
            enddo
          enddo
          enddo
          if     (baclin) then
          call horout_3d(utilk, artype,yrflag,time3,iexpt,.true.,
     &                ' bcl.spd.',                ! plot name
     &      'baroclinic_speed',                   ! ncdf name
     &     'baroclinic_sea_water_speed',          ! ncdf standard_name
     &                'm/s',                      ! units
     &                kf,kl,ltheta, frmt,ioin)
          else !total velocity
          call horout_3d(utilk, artype,yrflag,time3,iexpt,.true.,
     &                ' speed   ',                ! plot name
     &                'speed',                    ! ncdf name
     &                'sea_water_speed',          ! ncdf standard_name
     &                'm/s',                      ! units
     &                kf,kl,ltheta, frmt,ioin)
          endif
        endif
c
c --- -----------------------
c --- w-velocity
c --- -----------------------
c
c --- 'wvlio ' = w-velocity I/O unit (0 no I/O)
        call blkini2(ioin,i,'wvlio ','infio ')
        if     (i.eq.1) then
          infio = -99
        else
          infio = max(ioin,-1)
          ioin  = 0
        endif
        if (ioin.gt.0 .and. artype.ne.3) then  !not available for std. case
          do k= kf,kl
          do j=1,jj
            do i=1,ii
              utilk(i,j,k)=w(i,j,k)
              if (mthin) then
                if (p(i,j,k)+onecm.gt.p(i,j,k+1)) then
                  utilk(i,j,k)=flag
                endif
              else
                if (p(i,j,k)+onecm.gt.p(i,j,kk+1)) then
                  utilk(i,j,k)=flag
                endif
              endif
            enddo
          enddo
          enddo
          call horout_3d(utilk, artype,yrflag,time3,iexpt,.true.,
     &                ' w-veloc.',                   ! plot name
     &                'w_velocity',                  ! ncdf name
     &                'upward_sea_water_velocity',   ! ncdf standard_name
     &                'm/s',                         ! units
     &                kf,kl,ltheta, frmt,ioin)
        endif
c
c ---   --------------------
c ---   interface depth
c ---   --------------------
c
c ---   'infio ' = intf. k depth  I/O unit (0 no I/O)
        if     (infio.eq.-99) then
          call blkini(ioin,'infio ')
        else
          ioin = infio
        endif
        if (ioin.gt.0 .and. artype.ne.3) then    !not available for std. case
          do k= kf,kl
          do j=1,jj
            do i=1,ii
              utilk(i,j,k)=p(i,j,k+1)
            enddo
          enddo
          enddo
          call horout_3d(utilk, artype,yrflag,time3,iexpt,.true.,
     &                '  i.depth',                ! plot name
     &                'interface_depth',          ! ncdf name
     &                ' ',                        ! ncdf standard_name
     &                'm',                        ! units
     &                kf,kl,ltheta, frmt,ioin)
        endif
c
c ---   --------------------
c ---   layer thickness
c ---   --------------------
c
c ---   'thkio ' = lay.  k thick. I/O unit (0 no I/O)
        call blkini(ioin,'thkio ')
        if (ioin.gt.0) then
          do k= kf,kl
          do j=1,jj
            do i=1,ii
              if     (artype.ne.3) then
                utilk(i,j,k)=dp(  i,j,k)
              else
                utilk(i,j,k)=dpsd(i,j,k)  !artype==3
              endif
              if (p(i,j,k)+onecm.gt.p(i,j,kk+1)) then
                utilk(i,j,k)=flag
              endif
            enddo
          enddo
          enddo
          call horout_3d(utilk, artype,yrflag,time3,iexpt,.true.,
     &                '  thknss ',                ! plot name
     &                'layer_thickness',          ! ncdf name
     &                ' ',                        ! ncdf standard_name
     &                'm',                        ! units
     &                kf,kl,ltheta, frmt,ioin)
        endif
c
c ---   ----------------
c ---   temperature
c ---   ----------------
c
c ---   'temio ' = layer k temp I/O unit (0 no I/O, -ve MKS for NCOM)
        call blkini(ioin,'temio ')
        if (ioin.gt.0) then
          do k= kf,kl
          do j=1,jj
            do i=1,ii
              utilk(i,j,k)=temp(i,j,k)
              if (mthin) then
                if (p(i,j,k)+onecm.gt.p(i,j,k+1)) then
                  utilk(i,j,k)=flag
                endif
              else
                if (p(i,j,k)+onecm.gt.p(i,j,kk+1)) then
                  utilk(i,j,k)=flag
                endif
              endif
            enddo
          enddo
          enddo
          if     (kf.ne.1 .or. kf.ne.kl) then
            call horout_3d(utilk, artype,yrflag,time3,iexpt,.true.,
     &                  '  temp   ',                       ! plot name
     &                  'layer_temperature',               ! ncdf name
     &                  'sea_water_potential_temperature', ! ncdf standard_name
     &                  'degC',                            ! units
     &                  kf,kl,ltheta, frmt,ioin)
          else !kf==kl==1
            call horout(utilk, artype,yrflag,time3,iexpt,.true.,
     &                  ' sea surf. temp.  ',       ! plot name
     &                  'sst',                      ! ncdf name (mersea)
     &                  'sea_surface_temperature',  ! ncdf standard_name
     &                  'degC',                     ! units
     &                  0,ltheta, frmt,ioin)
          endif
        elseif (ioin.lt.0) then
          do k= kf,kl
          do j=1,jj
            do i=1,ii
              if (ip(i,j).ne.0) then
                utilk(i,j,k)=temp(i,j,k)
*               if (p(i,j,k).ge.p(i,j,kk+1)) then
*                 utilk(i,j,k)=0.0
*               endif
              else
                utilk(i,j,k)=0.0
              endif
            enddo
          enddo
          enddo
          call horout_3d(utilk, artype,yrflag,time3,iexpt,.true.,
     &                '  temp   ',                       ! plot name
     &                'layer_temperature',               ! ncdf name
     &                'sea_water_potential_temperature', ! ncdf standard_name
     &                'degC',                            ! units
     &                kf,kl,ltheta, frmt,-ioin)
        endif
c
c ---   -------------
c ---   salinity
c ---   -------------
c
c ---   'salio ' = lay.  k saln. I/O unit (0 no I/O)
        call blkini(ioin,'salio ')
        if (ioin.gt.0) then
          do k= kf,kl
          do j=1,jj
            do i=1,ii
              utilk(i,j,k)=saln(i,j,k)
              if (mthin) then
                if (p(i,j,k)+onecm.gt.p(i,j,k+1)) then
                  utilk(i,j,k)=flag
                endif
              else
                if (p(i,j,k)+onecm.gt.p(i,j,kk+1)) then
                  utilk(i,j,k)=flag
                endif
              endif
            enddo
          enddo
          enddo
          if     (kf.ne.1 .or. kf.ne.kl) then
            call horout_3d(utilk, artype,yrflag,time3,iexpt,.true.,
     &                  ' salinity',                ! plot name
     &                  'layer_salinity',           ! ncdf name
     &                  'sea_water_salinity',       ! ncdf standard_name
     &                  'psu',                      ! units
     &                  kf,kl,ltheta, frmt,ioin)
          else !kf==kl==1
            call horout(utilk, artype,yrflag,time3,iexpt,.true.,
     &                  'sea surf. salnity ',       ! plot name
     &                  'sss',                      ! ncdf name (mersea)
     &                  'sea_surface_salinity',     ! ncdf standard_name
     &                  'psu',                      ! units
     &                  0,ltheta, frmt,ioin)
          endif
        elseif (ioin.lt.0) then
          do k= kf,kl
          do j=1,jj
            do i=1,ii
              if (ip(i,j).ne.0) then
                utilk(i,j,k)=saln(i,j,k)
*               if (p(i,j,k).ge.p(i,j,kk+1)) then
*                 utilk(i,j,k)=0.0
*               endif
              else
                utilk(i,j,k)=0.0
              endif
            enddo
          enddo
          enddo
          call horout_3d(utilk, artype,yrflag,time3,iexpt,.true.,
     &                ' salinity',                ! plot name
     &                'layer_salinity',           ! ncdf name
     &                'sea_water_salinity',       ! ncdf standard_name
     &                'psu',                      ! units
     &                kf,kl,ltheta, frmt,-ioin)
        endif
c
c ---   -------------
c ---   density
c ---   -------------
c
c ---   'tthio ' = layer k density I/O unit (0 no I/O)
        call blkini(ioin,'tthio ')
        if (ioin.gt.0) then
          do k= kf,kl
          do j=1,jj
            do i=1,ii
              utilk(i,j,k)=th3d(i,j,k)
              if (mthin) then
                if (p(i,j,k)+onecm.gt.p(i,j,k+1)) then
                  utilk(i,j,k)=flag
                endif
              else
                if (p(i,j,k)+onecm.gt.p(i,j,kk+1)) then
                  utilk(i,j,k)=flag
                endif
              endif
            enddo
          enddo
          enddo
          call horout_3d(utilk, artype,yrflag,time3,iexpt,.true.,
     &                ' density ',                   ! plot name
     &                'layer_density',               ! ncdf name
     &                'sea_water_potential_density', ! ncdf standard_name
     &                'sigma',                       ! units
     &                kf,kl,ltheta, frmt,ioin)
        endif
c
c ---   -------------
c ---   tracers
c ---   -------------
c
        do ktr= 1,ntracr
c ---   'trcio ' = layer k tracer I/O unit (0 no I/O)
        call blkini(ioin,'trcio ')
        if (ioin.gt.0) then
          do k= kf,kl
          do j=1,jj
            do i=1,ii
              utilk(i,j,k)=trcr(i,j,k,ktr)
              if (mthin) then
                if (p(i,j,k)+onecm.gt.p(i,j,k+1)) then
                  utilk(i,j,k)=flag
                endif
              else
                if (p(i,j,k)+onecm.gt.p(i,j,kk+1)) then
                  utilk(i,j,k)=flag
                endif
              endif
            enddo
          enddo
          enddo
          call horout_3d(utilk, artype,yrflag,time3,iexpt,.true.,
     &                   trim(ctrc_title(ktr)),         ! plot name
     &                   trim(ctrc_lname(ktr)),         ! ncdf name
     &                   trim(ctrc_sname(ktr)),         ! ncdf standard_name
     &                   trim(ctrc_units(ktr)),         ! units
     &                   kf,kl,ltheta, frmt,ioin)
        endif
        enddo  !ktr= 1,ntracr
c
c ---   --------------------------
c ---   layer kinetic energy
c ---   --------------------------
c
        if     (artype.gt.1) then  ! mean or std. archive
c ---     'keio  ' = kinetic energy I/O unit (0 no I/O)
          call blkini(ioin,'keio  ')
          if (ioin.gt.0) then
            do k= kf,kl
            do j=1,jj
              do i=1,ii
                utilk(i,j,k)=ke(i,j,k)
                if (mthin) then
                  if (p(i,j,k)+onecm.gt.p(i,j,k+1)) then
                    utilk(i,j,k)=flag
                  endif
                else
                  if (p(i,j,k)+onecm.gt.p(i,j,kk+1)) then
                    utilk(i,j,k)=flag
                  endif
                endif
              enddo
            enddo
            enddo
            call horout_3d(utilk, artype,yrflag,time3,iexpt,.true.,
     &                  ' ke/mass ',                      ! plot name
     &                  'layer_kinetic_energy_per_mass',  ! ncdf name
     &                  ' ',                              ! ncdf standard_name
     &                  'm2/s2',                          ! units
     &                  kf,kl,ltheta, frmt,ioin)
          endif
        endif  ! mean or std. archive
c
c ---   --------------------------
c ---   layer stream function
c ---   --------------------------
c
c ---   'sfnio ' = layer k strmfn. I/O unit (0 no I/O)
        call blkini(ioin,'sfnio ')
        if (ioin.gt.0) then
          do k= kf,kl
c
c ---      calculate strmf on the Q grid, mostly using v.k.
c ---      note that strmf(ii,jj) is always 0.0, but streamfunctions
c ---       are arbitraty w.r.t. a constant offset anyway
           do j=jj,1,-1
c            jm1 = max(1,j-1)  !not accurate for j=1, when using a subregion.
            jm1=j-1
            if(j0.eq.0)jm1=max(j-1,1)
             i   = ii
             im1 = ii-1
             if (iu(i,j).eq.1) then
c ---          accurate transport calculation, allowing for depthu
               depthu = min(p(i,j,kk+1),p(im1,j,kk+1))
               ubi = u(i,j,k)*
     &                 max(0.0,
     &                     min(depthu,0.5*(p(i,j,k+1)+p(im1,j,k+1)))-
     &                     min(depthu,0.5*(p(i,j,k  )+p(im1,j,k  ))) )
             else
               ubi = 0.0
             endif
             if     (j.eq.jj) then
               strmfu = 0.0
             else
               strmfu = strmfu + ubi*scux(i,j)
             endif
             strmf(i,j) = strmfu
             strmft     = strmfu
             do i=ii-1,1,-1
               if (iv(i,j).eq.1) then
c ---            accurate transport calculation, allowing for depthv
                 depthv = min(p(i,j,kk+1),p(i,jm1,kk+1))
                 vbi = v(i,j,k)*
     &                   max(0.0,
     &                       min(depthv,0.5*(p(i,j,k+1)+p(i,jm1,k+1)))-
     &                       min(depthv,0.5*(p(i,j,k  )+p(i,jm1,k  ))) )
               else
                 vbi = 0.0
               endif
               strmft = strmft - vbi*scvy(i,j)
               strmf(i,j) = strmft  !m^3/s
             enddo !i
           enddo !j
c
c ---     interpolate from q-grid to p-grid.
          do j=1,jj
            do i=1,ii
              if (ip(i,j).ne.0 .and. i.lt.ii .and. j.lt.jj) then
                s1=0.0
                s2=0.0
                if     (iq(i,  j)  .eq.1) then
                  s1 = s1 + 1.0
                  s2 = s2 + strmf(i,  j)
                endif
                if     (iq(i+1,j)  .eq.1) then
                  s1 = s1 + 1.0
                  s2 = s2 + strmf(i+1,j)
                endif
                if     (iq(i,  j+1).eq.1) then
                  s1 = s1 + 1.0
                  s2 = s2 + strmf(i,  j+1)
                endif
                if     (iq(i+1,j+1).eq.1) then
                  s1 = s1 + 1.0
                  s2 = s2 + strmf(i+1,j+1)
                endif
                if     (s1.ne.0.0) then
                  utilk(i,j,k)=s2/s1
                else
                  utilk(i,j,k)=flag
                endif
              else
                utilk(i,j,k)=flag
              endif
            enddo
          enddo
c
ccc       call zebra(utilk(1,1,k),ii,ii,jj)
ccc       write (*,'('' shown above: layer''i3'' stream function'')') k
          enddo !k=kf,kl
c
          call horout_3d(utilk, artype,yrflag,time3,iexpt,.true.,
     &                'strmf (V)',                ! plot name
     &                'layer_streamfunction_v',   ! ncdf name
     &                'ocean_streamfunction',     ! ncdf standard_name
     &                'm3/s',                     ! units
     &                kf,kl,ltheta, frmt,ioin)
c
          do k= kf,kl
c
c ---     calculate strmf on the Q grid, mostly using u.k.
c ---     note that strmf(ii,jj) is always 0.0, but streamfunctions
c ---      are arbitraty w.r.t. a constant offset anyway
          do i=ii,1,-1
            im1 = max(1,i-1) !not accurate for i=1, when using a subregion.
            j   = jj
c            jm1 = jj-1
            jm1=j-1
            if(j0.eq.0)jm1=max(j-1,1)
            if (iv(i,j).eq.1) then
c ---         accurate transport calculation, allowing for depthv
              depthv = min(p(i,j,kk+1),p(i,jm1,kk+1))
              vbi = v(i,j,k)*
     &                max(0.0,
     &                    min(depthv,0.5*(p(i,j,k+1)+p(i,jm1,k+1)))-
     &                    min(depthv,0.5*(p(i,j,k  )+p(i,jm1,k  ))) )
            else
              vbi = 0.0
            endif
            if     (i.eq.ii) then
              strmfv = 0.0
            else
              strmfv = strmfv - vbi*scvy(i,j)
            endif
            strmf(i,j) = strmfv
            strmft     = strmfv
            do j=jj-1,1,-1
              if (iu(i,j).eq.1) then
c ---           accurate transport calculation, allowing for depthu
                depthu = min(p(i,j,kk+1),p(im1,j,kk+1))
                ubi = u(i,j,k)*
     &                  max(0.0,
     &                      min(depthu,0.5*(p(i,j,k+1)+p(im1,j,k+1)))-
     &                      min(depthu,0.5*(p(i,j,k  )+p(im1,j,k  ))) )
              else
                ubi = 0.0
              endif
              strmft     = strmft + ubi*scux(i,j)
              strmf(i,j) = strmft  !m^3/s
            enddo !i
          enddo !j
c
c ---     interpolate from q-grid to p-grid.
          do j=1,jj
            do i=1,ii
              if (ip(i,j).ne.0 .and. i.lt.ii .and. j.lt.jj) then
                s1=0.0
                s2=0.0
                if     (iq(i,  j)  .eq.1) then
                  s1 = s1 + 1.0
                  s2 = s2 + strmf(i,  j)
                endif
                if     (iq(i+1,j)  .eq.1) then
                  s1 = s1 + 1.0
                  s2 = s2 + strmf(i+1,j)
                endif
                if     (iq(i,  j+1).eq.1) then
                  s1 = s1 + 1.0
                  s2 = s2 + strmf(i,  j+1)
                endif
                if     (iq(i+1,j+1).eq.1) then
                  s1 = s1 + 1.0
                  s2 = s2 + strmf(i+1,j+1)
                endif
                if     (s1.ne.0.0) then
                  utilk(i,j,k)=s2/s1
                else
                  utilk(i,j,k)=flag
                endif
              else
                utilk(i,j,k)=flag
              endif
            enddo
          enddo
c
ccc       call zebra(utilk(1,1,k),ii,ii,jj)
ccc       write (*,'('' shown above: layer''i3'' stream function'')') k
          enddo !k=kf,kl
c
          call horout_3d(utilk, artype,yrflag,time3,iexpt,.true.,
     &                'strmf (U)',                ! plot name
     &                'layer_streamfunction_u',   ! ncdf name
     &                'ocean_streamfunction',     ! ncdf standard_name
     &                'm3/s',                     ! units
     &                kf,kl,ltheta, frmt,ioin)
        endif
c
      enddo  !layer loop
c
      call xcstop('(normal)')
      end
