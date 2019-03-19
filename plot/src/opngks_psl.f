      subroutine opngks
      implicit none
c
c     landscape color postscript graphics output
c     output file from environment variable NCARG_GKS_PS
c     or NCARG_GKS_PS_CYMK with default of 'gmeta1.ps'.
c
c     if NCARG_GKS_PS_CYMK is defined then the output is in the 
c     CMYK color model, otherwise in the default RGB color model.
c
      character cfile*256
      logical   lcymk
      integer   idummy
c
      external utilbd
c
      cfile = ' '
      call getenv('NCARG_GKS_PS_CYMK',cfile)
      lcymk = cfile.ne.' '
      if     (lcymk) then
        write(6,'(a)')
     &    ' calling landscape CYMK PostScript version of opngks'
        write(6,'(a,a)') ' filename: ',cfile(1:len_trim(cfile))
      else
        call getenv('NCARG_GKS_PS',cfile)
        if     (cfile.eq.' ') then
          cfile = 'gmeta1.ps'
        endif
        write(6,'(a)')
     &    ' calling landscape RGB PostScript version of opngks'
        write(6,'(a,a)') ' filename: ',cfile(1:len_trim(cfile))
      endif
c
      call gopks(6,idummy)
      call ngsetc('ME',cfile)
      if     (lcymk) then
        call ngseti('CM',0)
      endif
      call gopwk(1,2,26)
      call gacwk(1)
      call guwk( 1,1)
      call gstxfp(-22,2)  ! Helvetica-Bold
      END
