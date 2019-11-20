      program hycom_month_day
      implicit none
c
c     input:  ddd (ordinal day)
c     output: ddd, months that span the day, and their weights (linear interpolation)
c
c     example: echo 360 | hycom_month_day
c     result:  360 12 01 0.6967211 0.3032789
c
      integer jday,mon0,mon1
      real*4  w0,w1,x
c
      read(5,*) jday
      x    = 1.0 + mod( jday+366.0-15.25,366.0)/30.5
      mon0 = x
      mon1 = mod(mon0,12)+1
      w1   = mod(x,1.0) !fraction
      w0   = 1.0-w1
      write(6,'(i3.3,2i3.2,2f10.7)') jday,mon0,mon1,w0,w1
      call exit(0)
      end
