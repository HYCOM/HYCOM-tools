      module mod_za
      use mod_xc  ! HYCOM communication API
c
      implicit none
c
c --- HYCOM I/O interface.
c
c --- See README.src.mod_za for more details.
c
#if defined(SERIAL_IO)
      private zaiordd,zaiowrd
#else
      integer, save, private ::
     &  file_info_zaiord,file_count_zaiord,
     &  file_info_zaiowr,file_count_zaiowr
#endif
c
      integer, private, parameter :: uaoff = 1000 + uoff  !array I/O unit offset
      integer, private, save      :: n2drec
c
      real*4,private,save,allocatable::w(:),wminy(:),wmaxy(:),htmp(:)

      contains

#if defined(SERIAL_IO)
# include "mod_za_mp1.h"
#else
# include "mod_za_mp.h"
#endif

      end module mod_za

#if defined(ENDIAN_IO)
      subroutine zaio_endian(a,n)
      implicit none
c
      integer,         intent(in)    :: n
      real*4, intent(inout) :: a(n)  ! 4-bytes
c
c**********
c*
c 1)  swap the endian-ness of the array.
c
c 2)  assumes integer(kind=1) and integer(kind=4) ocupy one and four
c     bytes respectively.
c*
c**********
c
      integer         k
      real*4          b
      integer(kind=4) ib
      equivalence     (b,ib)
c
      integer(kind=4) ii4,   io4     ! 4-bytes
      integer(kind=1) ii1(4),io1(4)  ! 1-byte
      equivalence    (ii4,ii1(1)), (io4,io1(1))  ! non-standard f90
c
      do k= 1,n
        b   = a(k)
        ii4 = ib
        io1(1) = ii1(4)
        io1(2) = ii1(3)
        io1(3) = ii1(2)
        io1(4) = ii1(1)
        ib = io4
        a(k) = b
      enddo
      return
      end subroutine zaio_endian
#endif
