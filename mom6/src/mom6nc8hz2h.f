      program mom6nc8hz2h
      use mod_mom6r8 ! mom6 real*8 array interface
c
      implicit none
c
c --- Convert MOM6 hz (m) to h (kg m^-2)
c --- Rename MOM6 Boussinesq h in restart file to hz before running and
c --- the corresponding non-Boussinesq h will be output to the same file.
c --- Also reads temp and salt.
c
      character*256    flnm_r
      integer          i,j,k,lp,mro
      integer          itest,jtest
      double precision xmin,xmax
      double precision time3(3)
c
      lp = 6
c
c --- 'flnm_r' = name of mom6 file (input and output)
c --- 'itest ' = i-index for debugging printout (0 no debug)
c --- 'jtest ' = j-index for debugging printout (0 no debug)
c
      read (*,'(a)') flnm_r
      i = len_trim(flnm_r)
      write (lp,'(2a)') ' input MOM6 file: ',trim(flnm_r)
      call flush(lp)
c
      call blkini(itest, 'itest ')
      call blkini(jtest, 'jtest ')
c
c --- mom6 dimensions
c
      call rd_dimen(nto,mto,kk,mro, flnm_r,'Temp')
c
      write(lp,*) 
      write(lp,*) 'nto,mto,kk = ',nto,mto,kk
      write(lp,*) 
      call flush(lp)
c
c --- array allocation
c
      call mom6r8_alloc
c
c --- read the mom6 file.
c
        call rd_out3nc8(nto,mto,kk,1,
     &                  hz,
     &                  time3,  !HYCOM time
     &                  "hz",flnm_r)
        call flush(lp)
        write(lp,*) 
        write(lp,*) 'rd_out3nc8, hz,   time = ',time3(3)
        call flush(lp)
        call rd_out3nc8(nto,mto,kk,1,
     &                  temp,
     &                  time3,  !HYCOM time
     &                  "Temp",flnm_r)
        call flush(lp)
        write(lp,*) 'rd_out3nc8, temp, time = ',time3(3)
        call flush(lp)
        call rd_out3nc8(nto,mto,kk,1,
     &                  saln,
     &                  time3,  !HYCOM time
     &                  "Salt",flnm_r)
        call flush(lp)
        write(lp,*) 'rd_out3nc8, salt, time = ',time3(3)
        write(lp,*) 
        call flush(lp)
c
c ---   convert from h to hz
c
       call dens_wright(temp,saln,hz,rho,nto,mto,kk, itest,jtest)
        xmin =  1.e20
        xmax = -1.e20
        do j= 1,mto
          do i= 1,nto
            if     (hz(i,j,1).eq.1.0e20) then
              h(i,j,:) = 1.e20
            else
              do k= 1,kk
                h(i,j,k) = hz(i,j,k) * rho(i,j,k)
                xmin = min( xmin, h(i,j,k) )
                xmax = max( xmax, h(i,j,k) )
              enddo
            endif
          enddo !i
        enddo !j
        if     (itest.ne.0) then
          do k= 1,kk
            write(lp,'(a,i3,a,f9.3,f9.3,f14.3)') 
     &        'k =',k,
     &        ' R,hz,h =',rho(itest,jtest,k),
     &                     hz(itest,jtest,k),
     &                      h(itest,jtest,k)
          enddo
        endif
c
c ---   write out h
c
        call wr_out3nc8(nto,mto,kk,1,
     &                  h,time3, 
     &                  "h","Layer Thickness","kg m-2","Temp", flnm_r)
        write(lp,'(a,a,f14.4,2g15.6)') 
     &    'h',': day,min,max =',time3(3),xmin,xmax
        call flush(21)
c
      end program mom6nc8hz2h

      subroutine dens_wright(temp,saln,hz,rho,no,mo,kk, itest,jtest)
      use MOM_EOS_Wright ! MOM6 Wright equation of state
      implicit none
c
      integer no,mo,kk,itest,jtest
      real*8  temp(no,mo,kk),saln(no,mo,kk),hz(no,mo,kk),
     &         rho(no,mo,kk)
c
c --- calculate in-situ density using MOM6's Wright EOS
c
c     spval  = data void marker, 1.e20
      real, parameter :: spval=1.e20
c
      integer i,j,k,it
      real    qonem
      real*8  t5(5),s5(5),p5(5),r5(5)
c
      real,   allocatable ::  pij(:)
      real*8, allocatable :: prij(:)
c
      allocate(  pij(kk+1) )
      allocate( prij(kk+1) )
c
      do j= 1,mo
        do i= 1,no
          if     (temp(i,j,1).ne.spval) then
            pij(1) = 0.0
            do k=1,kk
              pij(k+1) = pij(k) + hz(i,j,k)
            enddo !k
                if     (i.eq.itest .and. j.eq.jtest) then
                  write(6,'(a,f12.5)') 'bot: ',pij(kk+1)
                endif !test
            prij(1) = 0.0d0
            do k= 1,kk
c ---         MOM6 uses int_density_dz_wright, which is faster for Wright,
c ---         but calculate_density works for all equations of state.
c
              t5(:) = temp(i,j,k)
              s5(:) = saln(i,j,k)
c
c ---         den at top of layer and initial estimate of prij(k+1)
c
              p5(1) = prij(k)
              call calculate_density_wright(t5,s5,p5, r5, 1,1)  !r5(1) only
              prij(k+1) = prij(k) + 9.8d0*r5(1)*(pij(k+1)-pij(k))
c
c ---         iterate to get layer in-situ density and pressure
c
              do it= 1,3
                p5(4) = 0.25*prij(k) + 0.75*prij(k+1)
                p5(5) =                     prij(k+1)
                call calculate_density_wright(t5,s5,p5, r5, 2,4)  !r5(2:5)
                    if     (i.eq.itest .and. j.eq.jtest) then
                      write(6,'(a,i3,i2,5f8.2)') 'p5:',k,it,p5(:)*1.d-4
                      write(6,'(a,i3,i2,5f8.2)') 'r5:',k,it,r5(:)
                    endif !test
c ---           Bode's (Boole's) Rule for integration
                r5(3) = 1.0d0/90.0d0*( 7.0d0*(r5(1)+r5(5))+
     &                                32.0d0*(r5(2)+r5(4))+
     &                                12.0d0* r5(3)        )
                prij(k+1) = prij(k) + 9.8d0*r5(3)*(pij(k+1)-pij(k))
              enddo !it
              rho(i,j,k) = r5(3)
            enddo !k
          else
            rho(i,j,:) = spval
          endif
        enddo !i
      enddo !j
      return
      end
