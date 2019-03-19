      subroutine fram(ncount)
      ncount=ncount+1
      write (*,'('' frame advance no.'',i4)') ncount
      call frame
      return
      end
