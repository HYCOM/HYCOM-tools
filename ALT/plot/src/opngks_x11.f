      subroutine opngks
c
c     X11 graphics output
c
c     Add something like the followng to your .Xdefaults file:
c     Xgks*geometry:  612x612-13-292
c
      external utilbd
c
      write(6,'(a)')   ' calling X11 version of opngks'
c
      call gopks(6,idummy)
      call ngseti('PC',-1)
      call gopwk(1,2,8)
      call gacwk(1)
*     call gstxfp(-7,2)  ! Hershey Complex Roman
      call gstxfp(-4,2)  ! Hershey Simplex Roman
      END
