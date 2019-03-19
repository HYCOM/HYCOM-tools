      module mod_xc
      use mod_dimensions  !include 'dimensions.h'
      implicit none
c
c --- HYCOM communication interface.
c --- see README.src.mod_xc for more details.
c
      include 'unit_offset.h'
c
c --- tile dimensions and tile numbers (counting from 1), see xcspmd
      integer, public, save      :: ipr,  jpr,  ijpr,
     &                              mproc,nproc,mnproc,
     &                              mp_1st
#if defined(MPI)
c
c --- needed for some versions of mod_za
      integer, public, save      :: group_1st_in_row
#endif
c
c --- region type (-1==unknown,
c ---               0==  closed/closed,
c ---               1==periodic/closed,
c ---               2==periodic/arctic,
c ---               3==periodic/fplane
c ---               4==  closed/fplane)
      integer, public, save      :: nreg
c
c --- timers on, usually and default .true.
      logical, public, save      :: timer_on=.true.
c
c --- fill value for land, usually 0.0
      real,    public, save      :: vland
      real*4,  public, save      :: vland4  !xcget4 only
c
c --- xctilr halo options
      integer, public, parameter :: halo_ps=1, halo_pv=11,
     &                              halo_qs=2, halo_qv=12,
     &                              halo_us=3, halo_uv=13,
     &                              halo_vs=4, halo_vv=14
c
c --- xcsync stdout flushing options
      logical, public, parameter :: flush_lp=.true.,
     &                              no_flush=.false.
c
c --- generic subroutine names
      interface xcmaxr
         module procedure xcmaxr_0  ! rank 0 array (i.e. scalar)
         module procedure xcmaxr_1  ! rank 1 array
         module procedure xcmaxr_0o ! rank 0 array, old interface
         module procedure xcmaxr_1o ! rank 1 array, old interface
      end interface

      interface xcminr
         module procedure xcminr_0  ! rank 0 array (i.e. scalar)
         module procedure xcminr_1  ! rank 1 array
         module procedure xcminr_0o ! rank 0 array, old interface
         module procedure xcminr_1o ! rank 1 array, old interface
      end interface
#if defined(USE_ESMF)
c
c --- public data structures for ESMF, see xcspmd
      integer, public,  save :: countde1(iqr),
     &                          countde2(jqr)
#endif
c
c --- private timer variables, see xctmri
      character*6, private, dimension(97), save :: cc
      integer,     private,                save :: nxc
      integer,     private, dimension(97), save :: nc
      real*8,      private, dimension(97), save :: tc,t0
      real*8,      private, dimension(2),  save :: tcxc,tcxl
#if defined(MPI) || defined(SHMEM)
c
c --- private message passing data structures, see xcspmd
      integer, private, save :: idproc( 0: iqr+1,0:jqr+1),
     &                          idproc1(0:ijqr+1),idhalo(2),
     &                          i0_pe(iqr,jqr),ii_pe(iqr,jqr),
     &                          j0_pe(iqr,jqr),jj_pe(iqr,jqr),
     &                          mpe_1(     jqr),
     &                          mpe_e(     jqr)

      integer, private, allocatable :: mpe_i(:,:),npe_j(:),
     &                                 sum8t(:),sum8j(:)
  
      real, private, allocatable    :: ai(:,:),aj(:,:)

      integer, private, save      :: ilen,jlen
      integer, private, parameter :: maxld=64

      integer, private, save :: i1sum(iqr,jqr),iisum(iqr,jqr)
      integer, private, save :: m0_top,i0_st(iqr),ii_st(iqr),
     &                          mm_top,i0_gt(iqr),ii_gt(iqr),
     &                          m0_bot,i0_sb(iqr),ii_sb(iqr),
     &                          mm_bot,i0_gb(iqr),ii_gb(iqr)
      integer, private, save :: null_tile
#endif
#if defined(MPI)
      integer, private, save :: mpi_comm_hycom
#endif
c
c --- actual module subroutines
      contains


/* BARRIER       set a barrier; for SPMD versions      */
#if defined(MPI)
# define BARRIER call mpi_barrier(mpi_comm_hycom,mpierr)
#elif defined(SHMEM)
# define BARRIER call shmem_barrier_all()
#endif

/* BARRIER_MP    halo synchronization; SHMEM only      */
/* BARRIER_NP    halo synchronization; SHMEM only      */
#if defined(RINGB)
#define BARRIER_MP call xctbar(idproc(mproc-1,nproc),idproc(mproc+1,nproc))
#define BARRIER_NP call xctbar(idproc(mproc,nproc-1),idproc(mproc,nproc+1))
#else
#define BARRIER_MP BARRIER
#define BARRIER_NP BARRIER
#endif

#if defined(MPI)
/* #define MPISR */
/* MTYPE4        mpi type for real*4                   */
/* MTYPER        mpi type for real                     */
/* MTYPED        mpi type for real*8                   */
/* MTYPEI        mpi type for integer                  */
/* MPI_SEND      either mpi_send  or mpi_ssend         */
/* MPI_ISEND     either mpi_isend or mpi_issend        */
#if defined(NOMPIR8) /* LAM does not support mpi_real[48] */
#if defined(REAL4)
# define MTYPE4 mpi_real
# define MTYPER mpi_real
# define MTYPED mpi_double_precision
# define MTYPEI mpi_integer
#else /* REAL8 */
# define MTYPE4 mpi_real
# define MTYPER mpi_double_precision
# define MTYPED mpi_double_precision
# define MTYPEI mpi_integer
#endif
#else /* most MPI's allow mpi_real[48] */
#if defined(REAL4)
# define MTYPE4 mpi_real4
# define MTYPER mpi_real4
# define MTYPED mpi_real8
# define MTYPEI mpi_integer
#else /* REAL8 */
# define MTYPE4 mpi_real4
# define MTYPER mpi_real8
# define MTYPED mpi_real8
# define MTYPEI mpi_integer
#endif
#endif
#if defined(SSEND)
# define MPI_SEND mpi_ssend
# define MPI_ISEND mpi_issend
#else
# define MPI_SEND mpi_send
# define MPI_ISEND mpi_isend
#endif /* SSEND:else */
#endif /* MPI */

#if defined(SHMEM)
/* SHMEM_GET4    get real*4  variables                 */
/* SHMEM_GETR    get real    variables                 */
/* SHMEM_GETD    get real*8  variables                 */
/* SHMEM_GETI    get integer variables                 */
/* SHMEM_MYPE    return number of this PE (0...npes-1) */
/* SHMEM_NPES    return number of PEs                  */
#if defined(REAL4)
# define SHMEM_GET4 shmem_get32
# define SHMEM_GETR shmem_get32
# define SHMEM_GETD shmem_get64
# define SHMEM_GETI shmem_integer_get
#else /* REAL8 */
# define SHMEM_GET4 shmem_get32
# define SHMEM_GETR shmem_get64
# define SHMEM_GETD shmem_get64
# define SHMEM_GETI shmem_integer_get
#endif
# define SHMEM_MYPE shmem_my_pe
# define SHMEM_NPES shmem_n_pes
#endif /* SHMEM */

c
c-----------------------------------------------------------------------
c
c     auxillary routines that involve off-processor communication.
c     message passing version, contained in module mod_xc.
c
c     author:  Alan J. Wallcraft,  NRL.
c
c-----------------------------------------------------------------------
c
      subroutine xcastr(a, mnflg)
      implicit none
c
      real,    intent(inout) :: a(:)
      integer, intent(in)    :: mnflg
c
c**********
c*
c  1) broadcast array a to all tiles.
c
c  2) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    a               real           in/out    target array
c    mnflg           integer        input     node originator flag
c*
c**********
c
#if defined(MPI)
      include 'mpif.h'
      integer        mpierr,mpireq,mpistat
      common/xcmpii/ mpierr,mpireq(4),
     &               mpistat(mpi_status_size,4*iqr)
      save  /xcmpii/
#endif
c
      integer    nmax
      parameter (nmax=1024)
c
      real            b,c
      common/xcmaxr4/ b(nmax),c(nmax)
      save  /xcmaxr4/
c
      integer i,is0,isl,mn,n,nn
#if defined(TIMER)
c
      if     (nxc.eq.0) then
        call xctmr0( 9)
        nxc =  9
      endif
#endif
c
c     stripmine a.
c
      n = size(a)
c
      do is0= 0,n-1,nmax
        isl = min(is0+nmax,n)
        nn = isl - is0
        if     (mnproc.eq.mnflg) then
          do i= 1,nn
            b(i) = a(is0+i)
          enddo
        endif
#if defined(MPI)
        call mpi_bcast(b,nn,MTYPER,
     &                 idproc1(mnflg),
     &                 mpi_comm_hycom,mpierr)
#elif defined(SHMEM)
        BARRIER
        if     (mnproc.ne.mnflg) then
c         get from source processor
          call SHMEM_GETR(b,b,nn, idproc1(mnflg))
        endif
        BARRIER
#endif
        if     (mnproc.ne.mnflg) then
          do i= 1,nn
            a(is0+i) = b(i)
          enddo
        endif
      enddo  ! stripmine loop
#if defined(TIMER)
c
      if     (nxc.eq. 9) then
        call xctmr1( 9)
        nxc = 0
      endif
#endif
      return
      end subroutine xcastr

      subroutine xceget(aelem, a, ia,ja)
      implicit none
c
      real,    intent(out)   :: aelem
      real,    intent(in)    :: a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      integer, intent(in)    :: ia,ja
c
c**********
c*
c  1) find the value of a(ia,ja) on the non-tiled 2-D grid.
c
c  2) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    aelem           real           output    required element
c    a               real           input     source array
c    ia              integer        input     1st index into a
c    ja              integer        input     2nd index into a
c
c  3) the global variable vland is returned when outside active tiles.
c*
c**********
c
#if defined(MPI)
      include 'mpif.h'
      integer        mpierr,mpireq,mpistat
      common/xcmpii/ mpierr,mpireq(4),
     &               mpistat(mpi_status_size,4*iqr)
      save  /xcmpii/
#endif
c
c     double buffer to reduce the number of required barriers.
c
      real            elem(0:1)
      common/xcegetr/ elem
      save  /xcegetr/
c
      integer, save :: kdb = 0
c
      integer i,j,mp,np
#if defined(TIMER)
c
*     if     (nxc.eq.0) then
*       call xctmr0( 2)
*       nxc = 2
*     endif
#endif
c
      kdb = mod(kdb+1,2)  ! the least recently used of the two buffers
c
c     find the host tile.
c
      np = npe_j(ja)
      mp = mpe_i(ia,np)
c
      if      (mp.le.0) then
c
c       no tile.
c
        elem(kdb) = vland
      elseif  (mp.eq.mproc .and. np.eq.nproc) then
c
c       this tile.
c
        i = ia - i0
        j = ja - j0
c
        elem(kdb) = a(i,j)
#if defined(MPI)
        call mpi_bcast(elem(kdb),1,MTYPER,
     &                 idproc(mp,np),mpi_comm_hycom,mpierr)
#elif defined(SHMEM)
        BARRIER
#endif
      else
c
c       another tile.
c
#if defined(MPI)
        call mpi_bcast(elem(kdb),1,MTYPER,
     &                 idproc(mp,np),mpi_comm_hycom,mpierr)
#elif defined(SHMEM)
        BARRIER
        call SHMEM_GETR(elem(kdb),
     &                  elem(kdb),1,idproc(mp,np))
        ! no barrier needed here because of double buffering
#endif
      endif
      aelem = elem(kdb)
#if defined(TIMER)
c
*     if     (nxc.eq. 2) then
*       call xctmr1( 2)
*       nxc = 0
*     endif
#endif
      return
      end subroutine xceget

      subroutine xceput(aelem, a, ia,ja)
      implicit none
c
      integer, intent(in)    :: ia,ja
      real,    intent(in)    :: aelem
      real,    intent(inout) :: a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
c
c**********
c*
c  1) fill a single element in the non-tiled 2-D grid.
c
c  2) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    aelem           real           input     element value
c    a               real           in/out    target array
c    ia              integer        input     1st index into a
c    ja              integer        input     2nd index into a
c*
c**********
c
#if defined(MPI)
      include 'mpif.h'
      integer        mpierr,mpireq,mpistat
      common/xcmpii/ mpierr,mpireq(4),
     &               mpistat(mpi_status_size,4*iqr)
      save  /xcmpii/
#endif
c
      integer mp,np
#if defined(TIMER)
c
*     if     (nxc.eq.0) then
*       call xctmr0( 5)
*       nxc = 5
*     endif
#endif
      if     (i0.lt.ia .and. ia.le.i0+ii .and.
     &        j0.lt.ja .and. ja.le.j0+jj      ) then
c
c       this tile.
c
        a(ia-i0,ja-j0) = aelem
      endif
#if defined(TIMER)
c
*     if     (nxc.eq. 5) then
*       call xctmr1( 5)
*       nxc = 0
*     endif
#endif
      return
      end subroutine xceput

      subroutine xcgetc(iline)
      implicit none
c
      integer, intent(inout) :: iline(81)
c
c**********
c*
c  1) machine specific routine for broadcasting iline.
c
c  2) only use in zagetc (hence the name).
c*
c**********
c
#if defined(MPI)
      include 'mpif.h'
      integer        mpierr,mpireq,mpistat
      common/xcmpii/ mpierr,mpireq(4),
     &               mpistat(mpi_status_size,4*iqr)
      save  /xcmpii/
#endif
c
c     broadcast to all processors
c
#if defined(MPI)
      call mpi_bcast(iline,81,MTYPEI,
     &               idproc1(1),mpi_comm_hycom,mpierr)
#elif defined(SHMEM)
      BARRIER
      if     (mnproc.ne.1) then
        call SHMEM_GETI(iline,
     &                  iline,81, idproc1(1))
      endif
      ! no barrier needed here because zagetc is using two buffers
#endif
      return
      end subroutine xcgetc

      subroutine xchalt(cerror)
      implicit none
c
      character*(*), intent(in) :: cerror
c
c**********
c*
c  1) stop all processes.
c
c  2) only one processes need call this routine, i.e. it is for
c      emergency stops.  use 'xcstop' for ordinary stops called
c      by all processes.
c
c  3) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    cerror          char*(*)       input     error message
c*
c**********
c
#if defined(MPI)
      include 'mpif.h'
      integer        mpierr,mpireq,mpistat
      common/xcmpii/ mpierr,mpireq(4),
     &               mpistat(mpi_status_size,4*iqr)
      save  /xcmpii/
#endif
c
c     message passing version.
c
      if     (cerror.ne.' ') then
        write(lp,*) '**************************************************'
        write(lp,*) cerror
      endif
      write(lp,*) '**************************************************'
      write(lp,*) 'XCHALT CALLED ON PROC = ',mnproc,mproc,nproc
      write(lp,*) '**************************************************'
      call flush(lp)
c
#if defined(MPI)
      call mpi_abort(mpi_comm_hycom,9)
#else
      call abort()
#endif
      stop '(xchalt)'
      end subroutine xchalt

 
      subroutine xclput(aline,nl, a, i1,j1,iinc,jinc)
      implicit none
c
      integer, intent(in)    ::  nl,i1,j1,iinc,jinc
      real,    intent(in)    ::  aline(nl)
      real,    intent(inout) ::  a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
c
c**********
c*
c  1) fill a line of elements in the non-tiled 2-D grid.
c
c  2) aline(i) == aa(i1+i1*(i-1),j1+j1*(i-1)), for i=1...nl.
c     where aa is the non-tiled representation of a, and
c     one of iinc and jinc must be 0, and the other must be 1.
c     also updates the halo.
c
c  3) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    aline           real           input     line of element values
c    nl              integer        input     dimension of aline
c    a               real           in/out    target array
c    i1              integer        input     1st index into a
c    j1              integer        input     2nd index into a
c    iinc            integer        input     1st index increment
c    jinc            integer        input     2nd index increment
c*
c**********
c
      integer i,j
#if defined(TIMER)
c
      if     (nxc.eq.0) then
        call xctmr0( 5)
        nxc = 5
      endif
#endif
c
      if     (jinc.eq.0) then
        if     (j1-j0.ge.1-nbdy .and. j1-j0.le.jj+nbdy) then
          do i= max(1-nbdy,i1-i0),min(i1-i0+nl-1,ii+nbdy)
            a(i,j1-j0) = aline(i+i0-i1+1)
          enddo
          if     (nreg.ne.0 .and.
     &            i0.eq.0 .and. i1+nl-1.ge.itdm-nbdy+1) then  ! periodic
            do i= max(itdm-nbdy+1,i1),i1+nl-1
              a(i-itdm,j1-j0) = aline(i)
            enddo
          endif
          if     (nreg.ne.0 .and.
     &            i0+ii.eq.itdm .and. i1.le.nbdy) then        ! periodic
            do i= i1,min(nbdy,i1+nl-1)
              a(ii+i,j1-j0) = aline(i)
            enddo
          endif
        endif
      elseif (iinc.eq.0) then
        if     (i1-i0.ge.1-nbdy .and. i1-i0.le.ii+nbdy) then
          do j= max(1-nbdy,j1-j0),min(j1-j0+nl-1,jj+nbdy)
            a(i1-i0,j) = aline(j+j0-j1+1)
          enddo
        endif
        if     (nreg.ne.0 .and.
     &          i0.eq.0 .and. i1.ge.itdm-nbdy+1) then       ! periodic
          do j= max(1-nbdy,j1-j0),min(j1-j0+nl-1,jj+nbdy)
            a(i1-itdm,j) = aline(j+j0-j1+1)
          enddo
        endif
        if     (nreg.ne.0 .and.
     &          i0+ii.eq.itdm .and. i1.le.nbdy) then        ! periodic
          do j= max(1-nbdy,j1-j0),min(j1-j0+nl-1,jj+nbdy)
            a(ii+i1,j) = aline(j+j0-j1+1)
          enddo
        endif
      endif
#if defined(TIMER)
c
      if     (nxc.eq. 5) then
        call xctmr1( 5)
        nxc = 0
      endif
#endif
      return
      end subroutine xclput

      subroutine xcmaxr_0(a, mnflg)
      implicit none
c
      real,    intent(inout) :: a
      integer, intent(in)    :: mnflg
c
c**********
c*
c  1) replace scalar a with its element-wise maximum over all tiles.
c
c  2) mnflg selects which nodes must return the minimum
c        = 0; all nodes
c        = n; node number n (mnproc=n)
c
c  3) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    a               real           in/out    target variable
c    mnflg           integer        input     node return flag
c*
c**********
c
      real a1(1)
c
      a1(1) = a
      call xcmaxr_1(a1, mnflg)
      a = a1(1)
      return
      end subroutine xcmaxr_0

      subroutine xcmaxr_1(a, mnflg)
      implicit none
c
      real,    intent(inout) :: a(:)
      integer, intent(in)    :: mnflg
c
c**********
c*
c  1) replace array a with its element-wise maximum over all tiles.
c
c  2) mnflg selects which nodes must return the minimum
c        = 0; all nodes
c        = n; node number n (mnproc=n)
c
c  3) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    a               real           in/out    target array
c    mnflg           integer        input     node return flag
c*
c**********
c
#if defined(MPI)
      include 'mpif.h'
      integer        mpierr,mpireq,mpistat
      common/xcmpii/ mpierr,mpireq(4),
     &               mpistat(mpi_status_size,4*iqr)
      save  /xcmpii/
#endif
c
      integer    nmax
      parameter (nmax=1024)
c
      real            b,c
      common/xcmaxr4/ b(nmax),c(nmax)
      save  /xcmaxr4/
c
      integer i,is0,isl,mn,mnp,n,nn
#if defined(TIMER)
c
      if     (nxc.eq.0) then
        call xctmr0(10)
        nxc = 10
      endif
#endif
c
c     stripmine a.
c
      n = size(a)
c
      do is0= 0,n-1,nmax
        isl = min(is0+nmax,n)
        nn = isl - is0
        do i= 1,nn
          b(i) = a(is0+i)
        enddo
c
#if defined(MPI)
        if     (mnflg.eq.0) then
          call mpi_allreduce(b,c,nn,MTYPER,mpi_max,
     &                       mpi_comm_hycom,mpierr)
        else
          call mpi_reduce(   b,c,nn,MTYPER,mpi_max,
     &                       idproc1(mnflg),
     &                       mpi_comm_hycom,mpierr)
        endif
#elif defined(SHMEM)
        BARRIER
        mnp = max(1,mnflg)
        if     (mnproc.eq.mnp) then
c         form global maximum on one processor
          do i= 1,nn
            c(i) = b(i)
          enddo
          do mn= 1,ijpr
            if     (mn.ne.mnp) then
              call SHMEM_GETR(b,b,nn, idproc1(mn))
              do i= 1,nn
                c(i) = max(c(i),b(i))
              enddo
            endif !.ne.mnp
          enddo
          BARRIER
        elseif (mnflg.eq.0) then
c         get global maximum from 1st processor
          BARRIER
          call SHMEM_GETR(c,c,nn, idproc1(mnp))
        endif
        ! no barrier needed here because using two buffers (b and c)
#endif
        if     (mnflg.eq.0 .or. mnflg.eq.mnproc) then
          do i= 1,nn
            a(is0+i) = c(i)
          enddo
        endif
      enddo  ! stripmine loop
#if defined(TIMER)
c
      if     (nxc.eq.10) then
        call xctmr1(10)
        nxc = 0
      endif
#endif
      return
      end subroutine xcmaxr_1

      subroutine xcmaxr_0o(a)
      implicit none
c
      real, intent(inout) :: a
c
c**********
c*
c  1) replace scalar a with its element-wise maximum over all tiles.
c
c  2) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    a               real           in/out    target variable
c*
c**********
c
      integer mnflg
      real    a1(1)
c
      mnflg = 0  !all nodes
      a1(1) = a
      call xcmaxr_1(a1, mnflg)
      a = a1(1)
      return
      end subroutine xcmaxr_0o

      subroutine xcmaxr_1o(a)
      implicit none
c
      real, intent(inout) :: a(:)
c
c**********
c*
c  1) replace array a with its element-wise maximum over all tiles.
c
c  2) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    a               real           in/out    target array
c*
c**********
c
      integer mnflg
c
      mnflg = 0  !all nodes
      call xcmaxr_1(a, mnflg)
      return
      end subroutine xcmaxr_1o

      subroutine xcminr_0(a, mnflg)
      implicit none
c
      real,    intent(inout) :: a
      integer, intent(in)    :: mnflg
c
c**********
c*
c  1) replace scalar a with its element-wise minimum over all tiles.
c
c  2) mnflg selects which nodes must return the minimum
c        = 0; all nodes
c        = n; node number n (mnproc=n)
c
c  3) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    a               real           in/out    target variable
c    mnflg           integer        input     node return flag
c*
c**********
c
      real a1(1)
c
      a1(1) = a
      call xcminr_1(a1, mnflg)
      a = a1(1)
      return
      end subroutine xcminr_0

      subroutine xcminr_1(a, mnflg)
      implicit none
c
      real,    intent(inout) :: a(:)
      integer, intent(in)    :: mnflg
c
c**********
c*
c  1) replace array a with its element-wise minimum over all tiles.
c
c  2) mnflg selects which nodes must return the minimum
c        = 0; all nodes
c        = n; node number n (mnproc=n)
c
c  2) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    a               real           in/out    target array
c    mnflg           integer        input     node return flag
c*
c**********
c
#if defined(MPI)
      include 'mpif.h'
      integer        mpierr,mpireq,mpistat
      common/xcmpii/ mpierr,mpireq(4),
     &               mpistat(mpi_status_size,4*iqr)
      save  /xcmpii/
#endif
c
      integer    nmax
      parameter (nmax=1024)
c
      real            b,c
      common/xcmaxr4/ b(nmax),c(nmax)
      save  /xcmaxr4/
c
      integer i,is0,isl,mn,mnp,n,nn
#if defined(TIMER)
c
      if     (nxc.eq.0) then
        call xctmr0(10)
        nxc = 10
      endif
#endif
c
c     stripmine a.
c
      n = size(a)
c
      do is0= 0,n-1,nmax
        isl = min(is0+nmax,n)
        nn = isl - is0
        do i= 1,nn
          b(i) = a(is0+i)
        enddo
c
#if defined(MPI)
        if     (mnflg.eq.0) then
          call mpi_allreduce(b,c,nn,MTYPER,mpi_min,
     &                       mpi_comm_hycom,mpierr)
        else
          call mpi_reduce(   b,c,nn,MTYPER,mpi_min,
     &                       idproc1(mnflg),
     &                       mpi_comm_hycom,mpierr)
        endif
#elif defined(SHMEM)
        BARRIER
        mnp = max(1,mnflg)
        if     (mnproc.eq.mnp) then
c         form global minimum on one processor
          do i= 1,nn
            c(i) = b(i)
          enddo
          do mn= 1,ijpr
            if     (mn.ne.mnp) then
              call SHMEM_GETR(b,b,nn, idproc1(mn))
              do i= 1,nn
                c(i) = min(c(i),b(i))
              enddo
            endif !.ne.mnp
          enddo
          BARRIER
        elseif (mnflg.eq.0) then
c         get global minimum from 1st processor
          BARRIER
          call SHMEM_GETR(c,c,nn, idproc1(mnp))
        endif
        ! no barrier needed here because using two buffers (b and c)
#endif
        if     (mnflg.eq.0 .or. mnflg.eq.mnproc) then
          do i= 1,nn
            a(is0+i) = c(i)
          enddo
        endif
      enddo
#if defined(TIMER)
c
      if     (nxc.eq.10) then
        call xctmr1(10)
        nxc = 0
      endif
#endif
      return
      end subroutine xcminr_1

      subroutine xcminr_0o(a)
      implicit none
c
      real, intent(inout) :: a
c
c**********
c*
c  1) replace scalar a with its element-wise minimum over all tiles.
c
c  2) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    a               real           in/out    target variable
c*
c**********
c
      integer mnflg
      real    a1(1)
c
      mnflg = 0  !all nodes
      a1(1) = a
      call xcminr_1(a1, mnflg)
      a = a1(1)
      return
      end subroutine xcminr_0o

      subroutine xcminr_1o(a)
      implicit none
c
      real, intent(inout) :: a(:)
c
c**********
c*
c  1) replace array a with its element-wise minimum over all tiles.
c
c  2) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    a               real           in/out    target array
c*
c**********
c
      integer mnflg
c
      mnflg = 0  !all nodes
      call xcminr_1(a, mnflg)
      return
      end subroutine xcminr_1o

#if defined(USE_ESMF)
      subroutine xcspmd(mpi_comm_vm)
      implicit none
c
      integer mpi_comm_vm
#else
      subroutine xcspmd
#if defined(USE_CCSM3)
      use ccsm3
      use ccsm3_io
#endif
      implicit none
#endif
c
c**********
c*
c  1) initialize data structures that identify the tiles.
c     version for 1-d uniform tiles only.
c
c  2) data structures (public):
c      ipr     - 1st 2-D node dimension (<=iqr)
c      jpr     - 2nd 2-D node dimension (<=jqr)
c      ijpr    -     1-D node dimension (ipr*jpr)
c      mproc   - 1st 2-D node index
c      nproc   - 2nd 2-D node index
c      mnproc  -     1-D node index
c      mp_1st  - 1st node in this row of 2-D nodes,        mpe_1(nproc)
c      i0      - 1st dimension offset for this tile, i0_pe(mproc,nproc)
c      ii      - 1st dimension extent for this tile, ii_pe(mproc,nproc)
c      j0      - 2nd dimension offset for this tile, j0_pe(mproc,nproc)
c      jj      - 2nd dimension extent for this tile, jj_pe(mproc,nproc)
c      nreg    -     region type
c      vland   -     fill value for land (standard value 0.0)
c
c  3) data structures (private):
c      idproc  -     2-D node addresses, with periodic wrap
c      idproc1 -     1-D node addresses, with periodic wrap
c      idhalo  -     left and right halo target nodes
c      i0_pe   - 1st dimension tile offsets
c      ii_pe   - 1st dimension tile extents (<=idm)
c      j0_pe   - 2nd dimension tile offsets
c      jj_pe   - 2nd dimension tile extents (<=jdm)
c      mpe_1   - 1st node in each row of 2-D nodes
c      mpe_e   - end node in each row of 2-D nodes
c      mpe_i   - mapping from 1st global dimension to 2-D nodes
c      npe_j   - mapping from 2nd global dimension to 2-D nodes
c      i1sum   - local index of 1st partial sum on each tile
c      iisum   - number of partial sums on each tile
c      m0_top  - tile offset:       top neighbors (0:jpr-1)
c      mm_top  - tile extent:       top neighbors (<=jpr)
c      i0_st   - halo offsets: send top neighbors
c      ii_st   - halo lengths: send top neighbors
c      i0_gt   - halo offsets:  get top neighbors
c      ii_gt   - halo lengths:  get top neighbors
c      m0_bot  - tile offset:       bot neighbors (0:jpr-1)
c      mm_bot  - tile extent:       bot neighbors (<=jpr)
c      i0_sb   - halo offsets: send bot neighbors
c      ii_sb   - halo lengths: send bot neighbors
c      i0_gb   - halo offsets:  get bot neighbors
c      ii_gb   - halo lengths:  get bot neighbors
c
c  4) all data structures are based on the processor number.
c*
c**********
c
      integer    mxsum
c
      integer i,idm_in,ixsum,idm_t,
     &        j,jdm_in,jdm_t,jtdm_in,
     &        l,m,mm,mn,mypei,n,npesi
c
#if defined(MPI)
      include 'mpif.h'
      integer        mpierr,mpireq,mpistat
      common/xcmpii/ mpierr,mpireq(4),
     &               mpistat(mpi_status_size,4*iqr)
      save  /xcmpii/
c
c     standard mpi (message passing) version.
c
       character cvarin*6
       
#if defined(USE_ESMF)
      call mpi_comm_dup(mpi_comm_vm,    mpi_comm_hycom, mpierr)
#elif defined(USE_CCSM3)
      call mpi_comm_dup(mpi_comm_ocn,   mpi_comm_hycom, mpierr)
#else
      call mpi_init(mpierr)
      call mpi_comm_dup(mpi_comm_world, mpi_comm_hycom, mpierr)
#endif
c
      call mpi_comm_rank(mpi_comm_hycom, mypei, mpierr)
      call mpi_comm_size(mpi_comm_hycom, npesi, mpierr)
c
      mnproc = mypei + 1  ! mnproc counts from 1
#if defined(DEBUG_ALL)
      lp = 6
      write(lp,'(a,i5)') 'mnproc =',mnproc
      call xcsync(flush_lp)
#endif
c
c     total array dimensions from regional.grid.b
c
      open(unit=11,file='regional.grid.b',form='formatted',
     &     status='old',action='read')
c
      read( 11,*) idm_t,cvarin
      if (cvarin.ne.'idm   ') then
        write(lp,*)
        write(lp,*) 'error in xcspmd - regional.grid.b input ',cvarin,
     &                        ' but should be idm   '
        write(lp,*)
        stop
      endif
      read( 11,*) jdm_t,cvarin
      if (cvarin.ne.'jdm   ') then
        write(lp,*)
        write(lp,*) 'error in xcspmd - regional.grid.b input ',cvarin,
     &                        ' but should be jdm   '
        write(lp,*)
        call flush(lp)
        stop
      endif
c
      close(unit=11)
c
      itdm=idm_t
      jtdm=jdm_t
      idm=itdm
      jdm=(jtdm+npesi-1)/npesi
c
c   Allocate various 'buffer' arrays needed by subroutines in mod_xc
c
      allocate(mpe_i(itdm,jqr))
      allocate(npe_j(jtdm))
c
c  arrays  for subroutines xcsum  &xcsumj
c
      mxsum=(idm+4*nbdy)/(2*nbdy+1)
      allocate(sum8t(mxsum*jdm))
      allocate(sum8j(jdm))
c
c  arrays  for subroutine xctilr
c      
      ilen= idm        *maxld*nbdy+64
      jlen=(jdm+2*nbdy)*maxld*nbdy+64
      allocate(ai(ilen,4))
      allocate(aj(jlen,4))
    	

c	
#elif defined(SHMEM)
      integer SHMEM_MYPE,SHMEM_NPES
c
c     shmem version.
c
      call start_pes(0)
c
      mypei = SHMEM_MYPE()
      npesi = SHMEM_NPES()
c
      mnproc = mypei + 1  ! mnproc counts from 1
      write(lp,*)'mnproc = ',mnproc
#else
      lp = 6
c
      write(lp,*)
      write(lp,*) '***** ERROR - UNDEFINED SPMD MACHINE TYPE *****'
      write(lp,*)
      call flush(lp)
      stop '(xcspmd)'
#endif
c
      vland  = 0.0
      vland4 = 0.0
c
      lp = 6
#if defined(T3E)
      open(unit=lp,delim='none')  ! forces open on stdout
#endif
c
c     calculate tile locations and sizes.
c     uniform 1-d distribution
c
      nreg    = 0      !assume closed
      ijpr    = npesi
      jpr     = npesi
      ipr     = 1
      idm_in  = itdm
      jdm_in  = (jtdm+jpr-1)/jpr
      jtdm_in = jpr*jdm_in  !>= jtdm
c
      n=1
        i0_pe(1,n) = 0
        ii_pe(1,n) = itdm
        j0_pe(1,n) = 0
        jj_pe(1,n) = jdm_in
      do n= 2,jpr
        i0_pe(1,n) = 0
        ii_pe(1,n) = itdm
        j0_pe(1,n) = j0_pe(1,n-1) + jj_pe(1,n-1)
        if     (n.le.jpr-(jtdm_in-jtdm)) then
          jj_pe(1,n) = jdm_in
        else
          jj_pe(1,n) = jdm_in-1
        endif
      enddo !j
      if     (j0_pe(1,jpr)+jj_pe(1,jpr).ne.jtdm) then
        call xcstop('xcspmd: bad jj_p calculation')
        stop '(xcspmd)'
      endif
c
#if defined(MPI)
c
c     mpi messages are sent and received by pe number (0:ijpr-1).
c
      null_tile = mpi_proc_null
c
      mn = 0
      do n= 1,jpr
        mpe_1(n) = 0
        do m= 1,ipr
          if     (ii_pe(m,n).eq.0) then
            idproc(m,n)   = null_tile
          else
            idproc1(mn+1) = mn
            idproc(m,n)   = mn
            mn = mn + 1
            if     (mnproc.eq.mn) then
              mproc = m
              nproc = n
            endif
            mpe_e(n) = m
            if     (mpe_1(n).eq.0) then
              mpe_1(n) = m
            endif
          endif
        enddo
      enddo
c
      mp_1st = mpe_1(nproc)  !1st node in this row (public)
c
c     mpi-2 i/o group (public), see mod_za.
c
      if     (mproc.eq.mp_1st) then
        i = 1
        call mpi_comm_split(mpi_comm_hycom, i,0,
     &                      group_1st_in_row, mpierr)
      else
        i = 0
        call mpi_comm_split(mpi_comm_hycom, i,0,
     &                      group_1st_in_row, mpierr)
        call mpi_comm_free( group_1st_in_row, mpierr)
      endif
#elif defined(SHMEM)
c
c     shmem messages are sent and received by pe number (0:ijpr-1).
c
      null_tile = -1
c
      mn = 0
      do n= 1,jpr
        mpe_1(n) = 0
        do m= 1,ipr
          if     (ii_pe(m,n).eq.0) then
            idproc(m,n)   = null_tile
          else
            idproc1(mn+1) = mn
            idproc(m,n)   = mn
            mn = mn + 1
            if     (mnproc.eq.mn) then
              mproc = m
              nproc = n
            endif
            mpe_e(n) = m
            if     (mpe_1(n).eq.0) then
              mpe_1(n) = m
            endif
          endif
        enddo
      enddo
c
      mp_1st = mpe_1(nproc)  !1st node in this row (public)
#endif
c
      if     (mn.ne.ijpr) then
        if     (mnproc.eq.1) then
          write(lp,'(a,i5)') 'input: ijpr =',ijpr
          write(lp,'(a,i5)') 'calc:  ijpr =',mn
          call flush(lp)
        endif
        call xcstop('xcspmd: wrong number of sea tiles')
        stop '(xcspmd)'
      endif
c
      if     (nreg.eq.0) then
c
c       longitudinal tile dimension is closed (not periodic)
c
        do n= 1,jpr
          idproc(    0,n) = null_tile
          idproc(ipr+1,n) = null_tile
        enddo
      else
c
c       longitudinal tile dimension is potentially periodic.
c
        do n= 1,jpr
          idproc(    0,n) = null_tile
          idproc(ipr+1,n) = null_tile
c
          i = maxval((i0_pe(1:ipr,n)+ii_pe(1:ipr,n)))
          if     (i0_pe(1,n).eq.0 .and. i.eq.itdm) then
            idproc(         0,n) = idproc(mpe_e(n),n)
            idproc(mpe_e(n)+1,n) = idproc(       1,n)
          endif
        enddo
      endif
c
c     latitudinal tile dimension is closed
c
      do m= 0,ipr+1
        idproc(m,    0) = null_tile
        idproc(m,jpr+1) = null_tile
      enddo
c
c     1-d tiling logic is easier if assumed periodic.
c
      idproc1(     0) = idproc1(ijpr)
      idproc1(ijpr+1) = idproc1(   1)
c
c     mapping from global i,j to mp,np.
c     ia,ja is on tile mpe_i(ia,npe_j(ja)),npe_j(ja), 
c     or on no tile if mpe_i(ia,npe_j(ja)) is 0 or -1.
c
      do n= 1,jpr
        mpe_i(1:itdm,n) = 0  ! default is an empty tile
        do m= 1,ipr  ! i-map potentially varies with n
          if     (ii_pe(m,n).gt.0) then
            do i= i0_pe(m,n)+1,i0_pe(m,n)+ii_pe(m,n)
              mpe_i(i,n) = m
            enddo
            if     (m.ne.ipr) then
              if     (ii_pe(m+1,n).gt.0) then
                do i= i0_pe(m,n)+ii_pe(m,n)+1,i0_pe(m+1,n)
                  mpe_i(i,n) = -1  ! gap between tiles
                enddo
              endif
            endif
          endif
        enddo
        m = 1  ! only one j-map
          do j= j0_pe(m,n)+1,j0_pe(m,n)+jj_pe(m,n)
            npe_j(j)   = n
          enddo
      enddo
c
c     do each partial sum on the tile that owns its center point.
c       i1sum - local index of 1st partial sum on each tile
c       iisum - number of partial sums on each tile
c     see xcsum for how i1sum and iisum are used.
c
      ixsum = 0
      do n= 1,jpr
        do m= 1,ipr
          if     (ii_pe(m,n).le.0) then
            i1sum(m,n) =  0
            iisum(m,n) =  0
          else
            idhalo(1) = idproc(m-1,n)
            idhalo(2) = idproc(m+1,n)
            if     (idhalo(1).ne.null_tile .and. m.ne.1) then
              if (i0_pe(m,n).ne.i0_pe(m-1,n)+ii_pe(m-1,n)) then
                idhalo(1) = null_tile
              endif
            endif
            if     (idhalo(2).ne.null_tile .and. m.ne.mpe_e(n)) then
              if     (i0_pe(m,n)+ii_pe(m,n).ne.i0_pe(m+1,n)) then
                idhalo(2) = null_tile
              endif
            endif
            i1sum(m,n) = -99
            iisum(m,n) =  0
            do i= 1+nbdy,itdm+nbdy,2*nbdy+1
              if     (i0_pe(m,n).lt.i .and.
     &                              i.le.i0_pe(m,n)+ii_pe(m,n)) then
                iisum(m,n) = iisum(m,n) + 1
                if     (iisum(m,n).eq.1) then
                  i1sum(m,n) = i - nbdy - i0_pe(m,n)
                endif
              elseif (idhalo(1).eq.null_tile .and.
     &                i.gt.i0_pe(m,n)-nbdy   .and.
     &                i.le.i0_pe(m,n)             ) then
                iisum(m,n) = iisum(m,n) + 1
                if     (iisum(m,n).eq.1) then
                  i1sum(m,n) = i - nbdy - i0_pe(m,n)
                endif
              elseif (idhalo(2).eq.null_tile          .and.
     &                i.gt.i0_pe(m,n)+ii_pe(m,n)      .and.
     &                i.le.i0_pe(m,n)+ii_pe(m,n)+nbdy      ) then
                iisum(m,n) = iisum(m,n) + 1
                if     (iisum(m,n).eq.1) then
                  i1sum(m,n) = i - nbdy - i0_pe(m,n)
                endif
              endif
            enddo
          endif
          ixsum = max( ixsum, iisum(m,n) )
        enddo !m
      enddo !n
c
c     local tile extents.
c
      i0 = i0_pe(mproc,nproc)
      ii = ii_pe(mproc,nproc)
      j0 = j0_pe(mproc,nproc)
      jj = jj_pe(mproc,nproc)
c
c     left and right halo targets
c
      idhalo(1) = idproc(mproc-1,nproc)
      idhalo(2) = idproc(mproc+1,nproc)
c
      if     (idhalo(1).ne.null_tile .and. mproc.ne.1) then
c
c       is the left tile touching this one?
c
        if (i0.ne.i0_pe(mproc-1,nproc)+ii_pe(mproc-1,nproc)) then
          idhalo(1) = null_tile
        endif
      endif
c
      if     (idhalo(2).ne.null_tile .and. mproc.ne.mpe_e(nproc)) then
c
c       is the right tile touching this one?
c
        if     (i0+ii.ne.i0_pe(mproc+1,nproc)) then
          idhalo(2) = null_tile
        endif
      endif
c
c     local halo exchange data structures
c
c     m0_top - tile offset:       top neighbors
c     mm_top - tile extent:       top neighbors (<=jpr)
c     i0_st  - halo offsets: send top neighbors
c     ii_st  - halo lengths: send top neighbors
c     i0_gt  - halo offsets:  get top neighbors
c     ii_gt  - halo lengths:  get top neighbors
c     m0_bot - tile offset:       bot neighbors
c     mm_bot - tile extent:       bot neighbors (<=jpr)
c     i0_sb  - halo offsets: send bot neighbors
c     ii_sb  - halo lengths: send bot neighbors
c     i0_gb  - halo offsets:  get bot neighbors
c     ii_gb  - halo lengths:  get bot neighbors
c
c     note that send is also receive, and is w.r.t. the local  tile.
c     similarly get  is also put,     and is w.r.t. the remote tile.
c
      if     (nproc.eq.jpr) then
c       no top nieghbor (closed boundary)
        m0_top = 0
        mm_top = 0
      else
        n = nproc + 1
        m0_top = 0
        mm_top = 0
        m      = 0
        do i= 1,ii
          if     (mpe_i(i0+i,n).ne.m) then
            if     (mm_top.eq.0) then
              m0_top = mpe_i(i0+i,n) - 1
            elseif (m.ne.-1) then
              ii_st(mm_top) = i-1 - i0_st(mm_top) 
              ii_gt(mm_top) = ii_st(mm_top)
            endif
            m = mpe_i(i0+i,n)
            if     (m.gt.0) then
              mm_top        = mm_top + 1
              i0_st(mm_top) = i-1
              i0_gt(mm_top) = i-1 + i0-i0_pe(m,n)
            elseif (m.eq.0) then
              mm_top        = mm_top + 1
              i0_st(mm_top) = i-1
              i0_gt(mm_top) = i0_gt(mm_top-1) + ii_gt(mm_top-1)
*           elseif (m.eq.-1) then  !do nothing
            endif
          endif
        enddo
        if     (mm_top.gt.0) then
          if     (m.gt.0) then
            ii_st(mm_top) = ii - i0_st(mm_top) 
            ii_gt(mm_top) = ii_st(mm_top)
          elseif (m.eq.0) then
            mm_top = mm_top-1
*         elseif (m.eq.-1) then  !do nothing
          endif
        endif
      endif  !nproc.eq.1:else
c
      if     (nproc.eq.1) then
c       no bottom nieghbor (closed boundary)
        m0_bot = 0
        mm_bot = 0
      else
        n = nproc - 1
        m0_bot = 0
        mm_bot = 0
        m      = 0
        do i= 1,ii
          if     (mpe_i(i0+i,n).ne.m) then
            if     (mm_bot.eq.0) then
              m0_bot = mpe_i(i0+i,n) - 1
            elseif (m.ne.-1) then
              ii_sb(mm_bot) = i-1 - i0_sb(mm_bot) 
              ii_gb(mm_bot) = ii_sb(mm_bot)
            endif
            m = mpe_i(i0+i,n)
            if     (m.gt.0) then
              mm_bot        = mm_bot + 1
              i0_sb(mm_bot) = i-1
              i0_gb(mm_bot) = i-1 + i0-i0_pe(m,n)
            elseif (m.eq.0) then
              mm_bot        = mm_bot + 1
              i0_sb(mm_bot) = i-1
              i0_gb(mm_bot) = i0_gb(mm_bot-1) + ii_gb(mm_bot-1)
*           elseif (m.eq.-1) then  !do nothing
            endif
          endif
        enddo
        if     (mm_bot.gt.0) then
          if     (m.gt.0) then
            ii_sb(mm_bot) = ii - i0_sb(mm_bot) 
            ii_gb(mm_bot) = ii_sb(mm_bot)
          elseif (m.eq.0) then
            mm_bot = mm_bot-1
*         elseif (m.eq.-1) then  !do nothing
          endif
        endif
      endif  !nproc.eq.1:else
c
c     printout the tile data structures.
c
      if     (mnproc.eq.1) then
        write(lp,'(/a)')
     &   'mnproc mproc nproc     i0   ii     j0   jj  i1sum iisum'
        mn = 0
        do n= 1,jpr
          do m= 1,ipr
            if     (ii_pe(m,n).ne.0) then
              mn= mn + 1
              write(lp,'(i6,2i6,i7,i5,i7,i5,i7,i6)')
     &           mn,m,n,
     &           i0_pe(m,n),ii_pe(m,n),
     &           j0_pe(m,n),jj_pe(m,n),
     &           i1sum(m,n),iisum(m,n)
            endif
          enddo !m
        enddo !n
        write(lp,*)
#if defined(DEBUG_ALL)
        do n= 1,jpr
          write(lp,*) 'mpe_1,mpe_e = ',mpe_1(n),mpe_e(n)
        enddo
        do n= 1,jpr
          write(lp,*) 'mpe_i = ',mpe_i(:,n)
        enddo
        write(lp,*)
        write(lp,*) 'npe_j = ',npe_j(:)
        write(lp,*)
#endif
      endif
      call xcsync(flush_lp)
c
#if defined(DEBUG_ALL)
      do n= 1,jpr
        do m= 1,ipr
          if     (mproc.eq.m .and. nproc.eq.n) then
            write(lp,'(a,2i3,i4,i3,16(i5,i4))')
     &         'm,n,_top,_st = ',
     &          m,n,m0_top,mm_top,(i0_st(l),ii_st(l), l= 1,mm_top)
#if defined(SHMEM)
            write(lp,'(a,2i3,i4,i3,16(i5,i4))')
     &         'm,n,_top,_gt = ',
     &          m,n,m0_top,mm_top,(i0_gt(l),ii_gt(l), l= 1,mm_top)
#endif
            if     (m.eq.ipr) then
              write(lp,*) !blank line
            endif
          endif
          call xcsync(flush_lp)
        enddo !m
        do m= 1,ipr
          if     (mproc.eq.m .and. nproc.eq.n) then
            write(lp,'(a,2i3,i4,i3,16(i5,i4))')
     &         'm,n,_bot,_sb = ',
     &          m,n,m0_bot,mm_bot,(i0_sb(l),ii_sb(l), l= 1,mm_bot)
#if defined(SHMEM)
            write(lp,'(a,2i3,i4,i3,16(i5,i4))')
     &         'm,n,_bot,_gb = ',
     &          m,n,m0_bot,mm_bot,(i0_gb(l),ii_gb(l), l= 1,mm_bot)
#endif
            if     (m.eq.ipr) then
              write(lp,*) !blank line
            endif
          endif
          call xcsync(flush_lp)
        enddo !m
        do m= 1,ipr
          if     (mproc.eq.m .and. nproc.eq.n) then
            write(lp,'(a,2i3,3i5)')
     &         'm,n,id,idhalo   = ',
     &          m,n,idproc(m,n),idhalo
            if     (m.eq.ipr) then
              write(lp,*) !blank line
            endif
          endif
          call xcsync(flush_lp)
        enddo !m
      enddo !n
#endif /* DEBUG_ALL */
#if defined(USE_ESMF)
c
c --- Currently ESMF only implemented for simplest block decomposition
c
      do n= 1,jpr
        countde2(n) = jj_pe(1,n)
      enddo !n
      do m= 1,ipr
        countde1(m) = ii_pe(m,1)
        if     (minval(ii_pe(m,1:jpr)).ne.
     &          maxval(ii_pe(m,1:jpr))    ) then
          if     (mnproc.eq.1) then
          write(lp,'(a,3i5)') 'm,min(ii),max(ii) = ',
     &                         m,minval(ii_pe(m,1:jpr)),
     &                           maxval(ii_pe(m,1:jpr))
          call flush(lp)
          endif
          call xcstop('xcspmd: bad ii for ESMF')
          stop '(xcspmd)'
        endif
      enddo !m
#endif
c
c     mxsum large enough?
c
      if     (ixsum.gt.mxsum) then
        if     (mnproc.eq.1) then
        write(lp,'(a,2i5)') 'mxsum,ixsum =',mxsum,ixsum
        call flush(lp)
        endif
        call xcstop('Error in xcspmd - mxsum too small')
        stop '(xcspmd)'
      endif
c
c     initialize timers.
c
      call xctmri
#if defined(TIMER)
      call xctmrn( 1,'xcaget')
      call xctmrn( 2,'xceget')
      call xctmrn( 3,'xclget')
      call xctmrn( 4,'xcaput')
      call xctmrn( 5,'xcXput')
      call xctmrn( 6,'xcsum ')
      call xctmrn( 9,'xcastr')
      call xctmrn(10,'xcmaxr')
      call xctmrn(12,'xctilr')
#endif
      return
      end subroutine xcspmd

      subroutine xcstop(cerror)
      implicit none
c
      character*(*), intent(in) :: cerror
c
c**********
c*
c  1) stop all processes.
c
c  2) all processes must call this routine.
c     use 'xchalt' for emergency stops.
c
c  3) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    cerror          char*(*)       input     error message
c*
c**********
c
#if defined(MPI)
      include 'mpif.h'
      integer        mpierr,mpireq,mpistat
      common/xcmpii/ mpierr,mpireq(4),
     &               mpistat(mpi_status_size,4*iqr)
      save  /xcmpii/
#endif
c
c     print active timers.
c
      call xctmrp
c
c     message passing version, set barrier and then stop everything.
c     use xcsync as the barrier so that stdout is flushed.
c     note that the system will hang unless all processes call xcstop.
c
      call xcsync(flush_lp)
      if     (mnproc.eq.1 .and. cerror.ne.' ') then
        write(lp,*) '**************************************************'
        write(lp,*) cerror
        write(lp,*) '**************************************************'
        write(lp,*)
      endif
      call xcsync(flush_lp)
c
#if defined(USE_ESMF) || defined(USE_CCSM3)
c
c --- we may not be running on all processes, so call mpi_abort
c
      if     (mnproc.eq.1) then
        call mpi_abort(mpi_comm_hycom,9)
      endif
      call xcsync(flush_lp)
#elif defined(MPI)
      write(lp,*) 'mpi_finalize called on processor ',mnproc
      call xcsync(flush_lp)
      call mpi_finalize(mpierr)
c
#endif
      stop '(xcstop)'
      end subroutine xcstop

      subroutine xcsum(sum, a,mask)
      implicit none
c
      real*8,  intent(out)   :: sum
      real,    intent(inout) :: a(   1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      integer, intent(in)    :: mask(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
c
c**********
c*
c  1) sum a 2-d array, where mask==1
c
c  2) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    sum             real*8         output    sum of a
c    a               real           input     source array
c    mask            integer        input     mask array
c
c  3) sum is bit for bit reproducable for the same halo size, nbdy.
c*
c**********
c
#if defined(MPI)
      include 'mpif.h'
      integer        mpierr,mpireq,mpistat
      common/xcmpii/ mpierr,mpireq(4),
     &               mpistat(mpi_status_size,4*iqr)
      save  /xcmpii/
#endif
c
      integer    mxsum
c
      real*8 , save  ::sum8s
c
      real*8     zero8
      parameter (zero8=0.0)
c
      real*8  sum8
      real    vsave
      integer i,i1,j,l,mp,np
C------------------------------------------------------------------
      mxsum=(idm+4*nbdy)/(2*nbdy+1)
#if defined(TIMER)
c
      if     (nxc.eq.0) then
        call xctmr0( 6)
        nxc = 6
      endif
#endif
c
c     halo update so that 2*nbdy+1 wide strips are on chip.
c
      vsave = vland
      vland = 0.0
      call xctilr(a,1,1, nbdy,0, halo_ps)
      vland = vsave
c
c     row sums in 2*nbdy+1 wide strips.
c
!$OMP PARALLEL DO PRIVATE(j,i1,i,l,sum8)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1,jj
        do l= 1,iisum(mproc,nproc)
          i1   = i1sum(mproc,nproc) + (l-1)*(2*nbdy+1)
          sum8 = zero8
          do i= max(i1,1-nbdy),min(i1+2*nbdy,ii+nbdy,itdm-i0)
            if     (mask(i,j).eq.1) then
              sum8 = sum8 + a(i,j)
            endif
          enddo
          sum8t(l + (j-1)*iisum(mproc,nproc)) = sum8
        enddo
      enddo
!$OMP END PARALLEL DO
c
c     complete row sums on first processor in each row.
c
#if defined(SHMEM)
      BARRIER
#endif
      if     (mproc.eq.mpe_1(nproc)) then
        do j=1,jj
          sum8j(j) = zero8
          do l= 1,iisum(mproc,nproc)
            sum8j(j) = sum8j(j) + sum8t(l + (j-1)*iisum(mproc,nproc))
          enddo
*         write(lp,'(a,i3,i5,f12.2)') 'xcsum: np,j,sum = ',
*    &                                1,j,sum8j(j)
        enddo
c
c       remote sums.
c
        do mp= mpe_1(nproc)+1,mpe_e(nproc)
          l = iisum(mp,nproc)*jj
          if     (l.gt.0) then
#if defined(MPI)
            call MPI_RECV(sum8t,l,MTYPED,
     &                    idproc(mp,nproc), 9900,
     &                    mpi_comm_hycom, mpistat, mpierr)
#elif defined(SHMEM)
            call SHMEM_GETD(sum8t,
     &                      sum8t,l,idproc(mp,nproc))
#endif
            do j=1,jj
              do l= 1,iisum(mp,nproc)
                sum8j(j) = sum8j(j) + sum8t(l + (j-1)*iisum(mp,nproc))
              enddo
*             write(lp,'(a,i3,i5,f12.2)') 'xcsum: np,j,sum = ',
*    &                                    mp,j,sum8j(j)
            enddo
          endif
        enddo
#if defined(MPI)
      else
        l = iisum(mproc,nproc)*jj
        if     (l.gt.0) then
          call MPI_SEND(sum8t,l,MTYPED,
     &                  idproc(mpe_1(nproc),nproc), 9900,
     &                  mpi_comm_hycom, mpierr)
        endif
#endif
      endif
c
c     sum of row sums, on first processor.
c
#if defined(SHMEM)
      BARRIER
#endif
      if     (mnproc.eq.1) then
        sum8 = zero8
        do j= 1,jj
          sum8 = sum8 + sum8j(j)
        enddo
*       write(lp,'(a,i5,f12.2)') 'xcsum: jj,sum = ',jj,sum8
c
        do np= 2,jpr
          mp = mpe_1(np)
#if defined(MPI)
          call MPI_RECV(sum8j,jj_pe(mp,np),MTYPED,
     &                  idproc(mp,np), 9901,
     &                  mpi_comm_hycom, mpistat, mpierr)
#elif defined(SHMEM)
          call SHMEM_GETD(sum8j,
     &                    sum8j,jj_pe(mp,np),idproc(mp,np))
#endif
          do j= 1,jj_pe(mp,np)
            sum8 = sum8 + sum8j(j)
          enddo
*         write(lp,'(a,i5,f12.2)') 'xcsum: jj,sum = ',
*    &                             jj_pe(mp,np),sum8
        enddo
        sum8s = sum8
#if defined(MPI)
      elseif (mproc.eq.mpe_1(nproc)) then
        call MPI_SEND(sum8j,jj,MTYPED,
     &                idproc1(1), 9901,
     &                mpi_comm_hycom, mpierr)
#endif
      endif
c
c     broadcast result to all processors.
c
#if defined(MPI)
      call mpi_bcast(sum8s,1,MTYPED,
     &               idproc1(1),mpi_comm_hycom,mpierr)
#elif defined(SHMEM)
      BARRIER
      if     (mnproc.ne.1) then
        call SHMEM_GETD(sum8s,
     &                  sum8s,1,idproc1(1))
      endif
#endif
c
      sum = sum8s
#if defined(TIMER)
c
      if     (nxc.eq. 6) then
        call xctmr1( 6)
        nxc = 0
      endif
#endif

      return
      end subroutine xcsum

      subroutine xcsumj(sumj, a,mask)
      implicit none
c
      real*8,  intent(out)   :: sumj(jtdm)
      real,    intent(inout) :: a(   1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      integer, intent(in)    :: mask(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
c
c**********
c*
c  1) rwo-sum of a 2-d array, where mask==1, on first processor only.
c
c  2) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    sumj            real*8         output    row-sum of a
c    a               real           input     source array
c    mask            integer        input     mask array
c
c  3) sum is bit for bit reproducable for the same halo size, nbdy.
c*
c**********
c
#if defined(MPI)
      include 'mpif.h'
      integer        mpierr,mpireq,mpistat
      common/xcmpii/ mpierr,mpireq(4),
     &               mpistat(mpi_status_size,4*iqr)
      save  /xcmpii/
#endif
c
      integer    mxsum
c
      real*8 , save             ::s um8s
c
      real*8     zero8
      parameter (zero8=0.0)
c
      real*8  sum8
      real    vsave
      integer i,i1,j,l,mp,np
      
      mxsum=(idm+4*nbdy)/(2*nbdy+1)

#if defined(TIMER)
c
      if     (nxc.eq.0) then
        call xctmr0( 6)
        nxc = 6
      endif
#endif
c
c     halo update so that 2*nbdy+1 wide strips are on chip.
c
      vsave = vland
      vland = 0.0
      call xctilr(a,1,1, nbdy,0, halo_ps)
      vland = vsave
c
c     row sums in 2*nbdy+1 wide strips.
c
!$OMP PARALLEL DO PRIVATE(j,i1,i,l,sum8)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1,jj
        do l= 1,iisum(mproc,nproc)
          i1   = i1sum(mproc,nproc) + (l-1)*(2*nbdy+1)
          sum8 = zero8
          do i= i1,min(i1+2*nbdy,ii+nbdy,itdm-i0)
            if     (mask(i,j).eq.1) then
              sum8 = sum8 + a(i,j)
            endif
          enddo
          sum8t(l + (j-1)*iisum(mproc,nproc)) = sum8
        enddo
      enddo
!$OMP END PARALLEL DO
c
c     complete row sums on first processor in each row.
c
#if defined(SHMEM)
      BARRIER
#endif
      if     (mproc.eq.mpe_1(nproc)) then
        do j=1,jj
          sum8j(j) = zero8
          do l= 1,iisum(mproc,nproc)
            sum8j(j) = sum8j(j) + sum8t(l + (j-1)*iisum(mproc,nproc))
          enddo
*         write(lp,'(a,i3,i5,f12.2)') 'xcsum: np,j,sum = ',
*    &                                1,j,sum8j(j)
        enddo
c
c       remote sums.
c
        do mp= mpe_1(nproc)+1,mpe_e(nproc)
          l = iisum(mp,nproc)*jj
          if     (l.gt.0) then
#if defined(MPI)
            call MPI_RECV(sum8t,l,MTYPED,
     &                    idproc(mp,nproc), 9900,
     &                    mpi_comm_hycom, mpistat, mpierr)
#elif defined(SHMEM)
            call SHMEM_GETD(sum8t,
     &                      sum8t,l,idproc(mp,nproc))
#endif
            do j=1,jj
              do l= 1,iisum(mp,nproc)
                sum8j(j) = sum8j(j) + sum8t(l + (j-1)*iisum(mp,nproc))
              enddo
*             write(lp,'(a,i3,i5,f12.2)') 'xcsum: np,j,sum = ',
*    &                                    mp,j,sum8j(j)
            enddo
          endif
        enddo
#if defined(MPI)
      else
        l = iisum(mproc,nproc)*jj
        if     (l.gt.0) then
          call MPI_SEND(sum8t,l,MTYPED,
     &                  idproc(mpe_1(nproc),nproc), 9900,
     &                  mpi_comm_hycom, mpierr)
        endif
#endif
      endif
c
c     send row sums to first processor.
c
#if defined(SHMEM)
      BARRIER
#endif
      if     (mnproc.eq.1) then
        do j= 1,jj
          sumj(j) = sum8j(j)
        enddo
c
        do np= 2,jpr
          mp = mpe_1(np)
#if defined(MPI)
          call MPI_RECV(sum8j,jj_pe(mp,np),MTYPED,
     &                  idproc(mp,np), 9901,
     &                  mpi_comm_hycom, mpistat, mpierr)
#elif defined(SHMEM)
          call SHMEM_GETD(sum8j,
     &                    sum8j,jj_pe(mp,np),idproc(mp,np))
#endif
          do j= 1,jj_pe(1,np)
            sumj(j+j0_pe(1,np)) = sum8j(j)
          enddo
        enddo
#if defined(MPI)
      elseif (mproc.eq.mpe_1(nproc)) then
        call MPI_SEND(sum8j,jj,MTYPED,
     &                idproc1(1), 9901,
     &                mpi_comm_hycom, mpierr)
#endif
      endif
#if defined(TIMER)
c
      if     (nxc.eq. 6) then
        call xctmr1( 6)
        nxc = 0
      endif
#endif

      return
      end subroutine xcsumj

      subroutine xcsync(lflush)
      implicit none
c
      logical, intent(in) :: lflush
c
c**********
c*
c  1) barrier, no processor exits until all arrive (and flush stdout).
c
c  2) some MPI implementations only flush stdout as a collective
c     operation, and hence the lflush=.true. option to flush stdout.
c
c  3) typically this is just a wrapper to the "BARRIER" macro.
c*
c**********
c
#if defined(MPI)
      include 'mpif.h'
      integer        mpierr,mpireq,mpistat
      common/xcmpii/ mpierr,mpireq(4),
     &               mpistat(mpi_status_size,4*iqr)
      save  /xcmpii/
#endif
c
      if     (lflush) then
#if defined(MPI) && defined(AIX) && ! (defined(AIX_NOFL) || defined(USE_ESMF) || defined(USE_CCSM3))
        call mp_flush(1)  ! flushes stdout, and implies a barrier
#else
        call flush(lp)
        BARRIER
#endif
      else
        BARRIER
      endif
      return
      end subroutine xcsync

#if defined(SHMEM) && defined(RINGB)
      subroutine xctbar(ipe1,ipe2)
      implicit none
c
      integer, intent(in) :: ipe1,ipe2
c
c**********
c*
c  1) sync with processors ipe1 and ipe2.
c
c  2) this is a global collective operation, and the calls on ipe1
c     and ipe2 must list this processor as one of the two targets.
c
c  3) this is used in place of a global barrier in halo operations,
c     but it only provides syncronization of one or two processors 
c     with the local processor.
c
c  4) ipe1 and/or ipe2 can be null_tile, to indicate no processor.
c*
c**********
c
      integer    cache_line,ilarge
      parameter (cache_line=32, ilarge=2**30)
c
      integer        ibp
      common/halobp/ ibp(cache_line,-1:ijpr-1)
      save  /halobp/
c
      integer i
c
      integer icount
      save    icount
      data    icount / -1 /
c
      icount = mod(icount+1,ilarge)
      if     (icount.eq.0) then
        call shmem_barrier_all()
        do i= -1,ijpr-1
          ibp(1,i) = -1
        enddo
        call shmem_barrier_all()
      endif
c
      ibp(1,-1) = icount
      if     (ipe1.ne.null_tile) then
        call shmem_integer_put(ibp(1,mnproc-1),icount,1,ipe1)
      endif
      if     (ipe2.ne.null_tile) then
        call shmem_integer_put(ibp(1,mnproc-1),icount,1,ipe2)
      endif
      call shmem_fence
c
      i = -1
      do while (i.lt.icount)
c       this assignment statement must not be optimized away.
cdir$   suppress
        i = min(ibp(1,ipe1),ibp(1,ipe2))
      enddo
      return
      end subroutine xctbar
#endif /* SHMEM && RINGB */

      recursive subroutine xctilr(a,l1,ld,mh,nh,itype)
      implicit none
c
      integer, intent(in)    :: l1,ld,mh,nh,itype
      real,    intent(inout) :: a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,ld)
c
c**********
c*
c  1) update the tile overlap halo of a real array.
c
c  2) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    a               real           in/out    target array
c    l1              integer        input     3rd dim. start index
c    ld              integer        input     3rd dimension of a
c    mh              integer        input     1st (EW) update halo size
c    nh              integer        input     2nd (NS) update halo size
c    itype           integer        input     grid and field type
c
c  3) itype selects both the grid and field type
c        itype= 1; p-grid, scalar field
c        itype= 2; q-grid, scalar field
c        itype= 3; u-grid, scalar field
c        itype= 4; v-grid, scalar field
c        itype=11; p-grid, vector field
c        itype=12; q-grid, vector field
c        itype=13; u-grid, vector field
c        itype=14; v-grid, vector field
c     it is ignored here because all types are the same unless
c      the grid includes the arctic ocean
c
c  4) the global variable vland is returned by halos over land.
c*
c**********
c
      integer i,j,k,l,lg0,ls0,lm,m,mhl,nhl
c
#if defined(MPI)
      include 'mpif.h'
      integer        mpierr,mpireq,mpistat
      common/xcmpii/ mpierr,mpireq(4),
     &               mpistat(mpi_status_size,4*iqr)
      save  /xcmpii/
c
c     persistent communication handles.
c
      integer mpireqa(4*iqr),mpireqb(4),ilold,jlold,nreqa
      save    mpireqa,mpireqb,ilold,jlold,nreqa
c
      data ilold,jlold / 0,0 /
#endif /* MPI */
c         1         2         3         4         5         6         7
c1234567890123456789012345678901234567890123456789012345678901234567890
c     WRITE(6,*)'XCTILR:2205,l1,ld,mh,nh,itype',l1,ld,mh,nh,itype
c
c --- split large requests into smaller pieces
c
c     ilold=0
c     jlold=0

      if     (ld-l1+1.gt.maxld) then
      WRITE(6,*)'XCTILR:2210,inside recursive loop'
        do k= l1,ld,maxld
          l = min(k+maxld-1,ld)
          call xctilr(a,k,l,mh,nh,itype)
        enddo
        return
c     else
c     WRITE(6,*)'XCTILR:2215 Not in recursive loop,kdm=',kdm
      endif
c
#if defined(TIMER)
c
      if     (nxc.eq.0) then
        call xctmr0(12)
        nxc = 12
      endif
#endif
c
      mhl = max(0,min(mh,nbdy))
      nhl = max(0,min(nh,nbdy))
c     WRITE(6,*)'XCTILR:2228 mhl,nhl=',mhl,nhl
c
      if     (nhl.gt.0) then
        if     (jpr.eq.1) then
          do k= l1,ld
            do j= 1,nhl
              do i= 1,ii
                a(i, 1-j,k) = vland
                a(i,jj+j,k) = vland
              enddo
            enddo
          enddo
        else
          l = 0
          do i= 1,ii  ! outer loop to simplify multiple neighbor case
            do k= l1,ld
              do j= 1,nhl
                l = l + 1
                ai(l,1) = a(i,jj+1-j,k)
                ai(l,2) = a(i,     j,k)
                ai(l,3) = vland
                ai(l,4) = vland
              enddo
            enddo
          enddo
c     WRITE(6,*)'XCTILR:2250'
c         write(lp,'(a,6i6)') 'xctilr - nhl,l1,ld,ii,l,mnproc = ',
c    &                                  nhl,l1,ld,ii,l,mnproc
c         call xcsync(flush_lp)
c
#if defined(MPI)
          if     (jlold.ne.l) then
            if     (jlold.ne.0) then
              do i= 1,nreqa
                call mpi_request_free(mpireqa(i), mpierr)
              enddo
            endif
            jlold = l
c
c           loop through all neigboring tiles.
c
            l = 0
            do m= 1,mm_top
              l   = l + 1
              ls0 = i0_st(m)*nhl*(ld-l1+1)
              lm  = ii_st(m)*nhl*(ld-l1+1)
              call mpi_send_init(
     &          ai(ls0+1,1),lm,MTYPER,idproc(m0_top+m,nproc+1), 9905,
     &          mpi_comm_hycom, mpireqa(l), mpierr)
            enddo
            do m= 1,mm_bot
              l   = l + 1
              ls0 = i0_sb(m)*nhl*(ld-l1+1)
              lm  = ii_sb(m)*nhl*(ld-l1+1)
              call mpi_send_init(
     &          ai(ls0+1,2),lm,MTYPER,idproc(m0_bot+m,nproc-1), 9906,
     &          mpi_comm_hycom, mpireqa(l), mpierr)
            enddo
            do m= 1,mm_top
              l   = l + 1
              ls0 = i0_st(m)*nhl*(ld-l1+1)
              lm  = ii_st(m)*nhl*(ld-l1+1)
              call mpi_recv_init(
     &          ai(ls0+1,4),lm,MTYPER,idproc(m0_top+m,nproc+1), 9906,
     &          mpi_comm_hycom, mpireqa(l), mpierr)
            enddo
            do m= 1,mm_bot
              l   = l + 1
              ls0 = i0_sb(m)*nhl*(ld-l1+1)
              lm  = ii_sb(m)*nhl*(ld-l1+1)
              call mpi_recv_init(
     &          ai(ls0+1,3),lm,MTYPER,idproc(m0_bot+m,nproc-1), 9905,
     &          mpi_comm_hycom, mpireqa(l), mpierr)
            enddo
            nreqa = l
          endif
          if     (nreqa.gt.0) then
            call mpi_startall(nreqa, mpireqa,          mpierr)
            call mpi_waitall( nreqa, mpireqa, mpistat, mpierr)
          endif
#elif defined(SHMEM)
          BARRIER
c
c         loop through all neigboring tiles.
c
          do m= 1,mm_top
            if     (idproc(m0_top+m,nproc+1).ne.null_tile) then
              lg0 = i0_gt(m)*nhl*(ld-l1+1)
              ls0 = i0_st(m)*nhl*(ld-l1+1)
              lm  = ii_st(m)*nhl*(ld-l1+1)
              call SHMEM_GETR(ai(ls0+1,4),
     &                        ai(lg0+1,2),lm, idproc(m0_top+m,nproc+1))
            endif
          enddo
c
          do m= 1,mm_bot
            if     (idproc(m0_bot+m,nproc-1).ne.null_tile) then
              lg0 = i0_gb(m)*nhl*(ld-l1+1)
              ls0 = i0_sb(m)*nhl*(ld-l1+1)
              lm  = ii_sb(m)*nhl*(ld-l1+1)
              call SHMEM_GETR(ai(ls0+1,3),
     &                        ai(lg0+1,1),lm, idproc(m0_bot+m,nproc-1))
            endif
          enddo
#endif  /* MPI:SHMEM */
c
          l = 0
          do i= 1,ii  ! outer loop to simplify multiple neighbor case
            do k= l1,ld
              do j= 1,nhl
                l = l + 1
                a(i, 1-j,k) = ai(l,3)
                a(i,jj+j,k) = ai(l,4)
              enddo
            enddo
          enddo
        endif  ! jpr.eq.1:else
      endif  ! nhl.gt.0
c
      if     (mhl.gt.0) then
        if     (ipr.eq.1) then
          if     (nreg.eq.0) then
            do k= l1,ld
              do j= 1-nhl,jj+nhl
                do i= 1,mhl
                  a( 1-i,j,k) = vland
                  a(ii+i,j,k) = vland
                enddo
              enddo
            enddo
          else
            do k= l1,ld
              do j= 1-nhl,jj+nhl
                do i= 1,mhl
                  a( 1-i,j,k) = a(ii+1-i,j,k)
                  a(ii+i,j,k) = a(     i,j,k)
                enddo
              enddo
            enddo
          endif
        else
          l = 0
          do k= l1,ld
            do j= 1-nhl,jj+nhl
              do i= 1,mhl
                l = l + 1
                aj(l,1) = a(ii+1-i,j,k)
                aj(l,2) = a(     i,j,k)
                aj(l,3) = vland
                aj(l,4) = vland
              enddo
            enddo
          enddo
*         write(lp,'(a,6i6)') 'xctilr - mhl,l1,ld,jj,l,mnproc = ',
*    &                                  mhl,l1,ld,jj,l,mnproc
*         call xcsync(flush_lp)
#if defined(MPISR)
          call mpi_sendrecv(
     &          aj(1,1),l,MTYPER,idhalo(2), 9907,
     &          aj(1,4),l,MTYPER,idhalo(2), 9908,
     &          mpi_comm_hycom, mpistat, mpierr)
          call mpi_sendrecv(
     &          aj(1,2),l,MTYPER,idhalo(1), 9908,
     &          aj(1,3),l,MTYPER,idhalo(1), 9907,
     &          mpi_comm_hycom, mpistat, mpierr)
#elif defined(MPI)
          if     (ilold.ne.l) then
            if     (ilold.ne.0) then
              call mpi_request_free(mpireqb(1), mpierr)
              call mpi_request_free(mpireqb(2), mpierr)
              call mpi_request_free(mpireqb(3), mpierr)
              call mpi_request_free(mpireqb(4), mpierr)
            endif
            ilold = l
            call mpi_send_init(
     &            aj(1,1),l,MTYPER,idhalo(2), 9907,
     &            mpi_comm_hycom, mpireqb(1), mpierr)
            call mpi_send_init(
     &            aj(1,2),l,MTYPER,idhalo(1), 9908,
     &            mpi_comm_hycom, mpireqb(2), mpierr)
            call mpi_recv_init(
     &            aj(1,3),l,MTYPER,idhalo(1), 9907,
     &            mpi_comm_hycom, mpireqb(3), mpierr)
            call mpi_recv_init(
     &            aj(1,4),l,MTYPER,idhalo(2), 9908,
     &            mpi_comm_hycom, mpireqb(4), mpierr)
          endif
          call mpi_startall(4, mpireqb,          mpierr)
          call mpi_waitall( 4, mpireqb, mpistat, mpierr)
#elif defined(SHMEM)
          BARRIER_MP
          if     (idhalo(1).ne.null_tile) then
            call SHMEM_GETR(aj(1,3),
     &                      aj(1,1),l,idhalo(1))
          endif
          if     (idhalo(2).ne.null_tile) then
            call SHMEM_GETR(aj(1,4),
     &                      aj(1,2),l,idhalo(2))
          endif
          BARRIER_MP
#endif
          l = 0
          do k= l1,ld
            do j= 1-nhl,jj+nhl
              do i= 1,mhl
                l = l + 1
                a( 1-i,j,k) = aj(l,3)
                a(ii+i,j,k) = aj(l,4)
              enddo
            enddo
          enddo
        endif  ! ipr.eq.1:else
      endif  ! mhl.gt.0
#if defined(TIMER)
c
      if     (nxc.eq.12) then
        call xctmr1(12)
        nxc = 0
      endif
#endif
      
      return
      end subroutine xctilr

      subroutine xctmri
      implicit none
c
c**********
c*
c  1) initialize timers.
c
c  2) timers  1:32 are for message passing routines,
c     timers 33:80 are for general hycom routines,
c     timers 81:96 are for user selected routines.
c     timer     97 is the total time.
c
c  3) call xctmri    to initialize timers (called in xcspmd),
c     call xctmr0(n) to start timer n,
c     call xctmr1(n) to stop  timer n and add event to timer sum,
c     call xctnrn(n,cname) to register a name for timer n,
c     call xctmrp to printout timer statistics (called by xcstop).
c
c  4) time every 50-th event above 1,000.
c*
c**********
c
      integer i
c
      real*8     zero8
      parameter (zero8=0.0)
c
      nxc = 0
      do i= 1,97
        cc(i) = '      '
        nc(i) = 0
        tc(i) = zero8
      enddo
c
      call xctmrn(97,'total ')
      call xctmr0(97)
      return
      end subroutine xctmri

      subroutine xctmr0(n)
      implicit none
c
      integer, intent(in) :: n
c
c**********
c*
c  1) start timer n.
c
c  2) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    n               integer        input     timer number
c
c  3) time every 50-th event above 1,000.
c*
c**********
c
      real*8 wtime
c
#if defined(DEBUG_TIMER_ALL)
      if     (              cc(n).ne.'      ') then
        write(lp,'(i5,2x,a,a)') mnproc,'call ',cc(n)
        call flush(lp)
      endif
#endif
#if defined(DEBUG_TIMER)
      if     (n.gt.32 .and. cc(n).ne.'      ') then
        if     (mnproc.eq.1) then
        write(lp,*) 'call ',cc(n)
        call flush(lp)
        endif
      endif
#endif
      if     (timer_on) then
        if     (mod(nc(n),50).eq.0 .or. nc(n).le.1000) then
          t0(n) = wtime()
        endif
      endif !timer_on
      return
      end subroutine xctmr0

      subroutine xctmr1(n)
      implicit none
c
      integer, intent(in) :: n
c
c**********
c*
c  1) add time since call to xctim0 to timer n.
c
c  2) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    n               integer        input     timer number
c
c  3) time every 50-th event above 1,000.
c*
c**********
c
      real*8  wtime
c
      if     (timer_on) then
        if     (nc(n).gt.1000) then
          if     (mod(nc(n),50).eq.0) then
            tc(n) = tc(n) + 50.0*(wtime() - t0(n))
          endif
        else
          tc(n) = tc(n) + (wtime() - t0(n))
        endif
        nc(n) = nc(n) + 1
      endif !timer_on
#if defined(DEBUG_TIMER_ALL)
      if     (              cc(n).ne.'      ') then
        write(lp,'(i5,2x,a,a)') mnproc,'exit ',cc(n)
        call flush(lp)
      endif
#endif
#if defined(DEBUG_TIMER)
      if     (n.gt.32 .and. cc(n).ne.'      ') then
        if     (mnproc.eq.1) then
        write(lp,*) 'exit ',cc(n)
        call flush(lp)
        endif
      endif
#endif
      return
      end subroutine xctmr1

      subroutine xctmrn(n,cname)
      implicit none
c
      character*6, intent(in) :: cname
      integer,     intent(in) :: n
c
c**********
c*
c  1) register name of timer n.
c
c  2) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    n               integer        input     timer number
c    cname           char*(6)       input     timer name
c*
c**********
c
      cc(n) = cname
      return
      end subroutine xctmrn

#if defined(TIMER_ALLOUT)
      subroutine xctmrp
      implicit none
c
c**********
c*
c  1) print all active timers, on all processors.
c
c  2) on exit all timers are reset to zero.
c*
c**********
c
#if defined(MPI)
      include 'mpif.h'
      integer        mpierr,mpireq,mpistat
      common/xcmpii/ mpierr,mpireq(4),
     &               mpistat(mpi_status_size,4*iqr)
      save  /xcmpii/
#endif
c
      integer i,mn
c
      real*8     zero8
      parameter (zero8=0.0)
c
c     get total time.
c
      call xctmr1(97)
c
      call xcsync(flush_lp)  !includes a barrier for shmem
      if     (mnproc.ne.1) then
#if defined(MPI)
        call MPI_SEND(tc,97,MTYPED,
     &                idproc1(1), 9949,
     &                mpi_comm_hycom, mpierr)
#endif
      else   !mnproc.eq.1
        do mn= 1,ijpr
          if     (mn.ne.1) then
#if defined(MPI)
            call MPI_RECV(tc,97,MTYPED,
     &                    idproc1(mn), 9949,
     &                    mpi_comm_hycom, mpistat, mpierr)
#elif defined(SHMEM)
            call SHMEM_GETD(tc,tc,97,idproc1(mn))
#endif
          endif
          write(lp,6000) mn,ijpr
          tcxc(1) = zero8
          do i= 1,32
            if     (nc(i).ne.0) then
              if     (cc(i).ne.'      ') then
                write(lp,6100) cc(i),nc(i),tc(i),tc(i)/nc(i)
              else
                write(lp,6150)    i, nc(i),tc(i),tc(i)/nc(i)
              endif
              if     (cc(i)(1:2).eq.'xc') then
                tcxc(1) = tcxc(1) + tc(i)  !communication overhead
              endif
            endif
          enddo !i
          write(lp,6100) 'xc****',1,tcxc(1),tcxc(1)
          do i= 33,97
            if     (nc(i).ne.0) then
              if     (cc(i).ne.'      ') then
                write(lp,6100) cc(i),nc(i),tc(i),tc(i)/nc(i)
              else
                write(lp,6150)    i, nc(i),tc(i),tc(i)/nc(i)
              endif
            endif
          enddo !i
        enddo !mn
        write(lp,6200)
      endif !mnproc.ne.1:else
      call xcsync(flush_lp)  !includes a barrier for shmem
c
c     reset timers to zero.
c
      do i= 1,97
        nc(i) = 0
        tc(i) = zero8
      enddo
      tcxc(1) = zero8
c
c     start a new total time measurement.
c
      call xctmr0(97)
      return
c
 6000 format(/ /
     &    3x,' timer statistics, processor',i5,' out of',i5 /
     &    3x,'-----------------------------------------------' /)
 6100 format(3x,a6,
     &   '   calls =',i9,
     &   '   time =',f11.5,
     &   '   time/call =',f14.8)
 6150 format(3x,'   #',i2,
     &   '   calls =',i9,
     &   '   time =',f11.5,
     &   '   time/call =',f14.8)
 6200 format(/ /)
      end subroutine xctmrp
#else
      subroutine xctmrp
      implicit none
c
c**********
c*
c  1) print all active timers.
c
c  2) on exit all timers are reset to zero.
c*
c**********
c
#if defined(MPI)
      include 'mpif.h'
      integer        mpierr,mpireq,mpistat
      common/xcmpii/ mpierr,mpireq(4),
     &               mpistat(mpi_status_size,4*iqr)
      save  /xcmpii/
#endif
c
      integer i,mn,mnloc
c
      real*8     zero8
      parameter (zero8=0.0)
c
c     get total time.
c
      call xctmr1(97)
c
c     report time on the processor with the least communication overhead
c
      tcxc(2) = mnproc
      tcxc(1) = zero8
      do i= 1,32
        if     (nc(i).ne.0 .and. cc(i)(1:2).eq.'xc') then
          tcxc(1) = tcxc(1) + tc(i)  !communication overhead
        endif
      enddo !i
c
      if     (ijpr.ne.1) then
#if defined(MPI)
        tcxl(1) = tcxc(1)
        tcxl(2) = tcxc(2)
        call mpi_allreduce(tcxl,tcxc,1,
     &                     mpi_2double_precision,mpi_minloc,
     &                     mpi_comm_hycom,mpierr)
        mnloc = tcxc(2)  !processor with the least comm. overhead
        if     (mnproc.eq.1) then
          if     (mnloc.ne.1) then
            call MPI_RECV(tc,97,MTYPED,
     &                    idproc1(mnloc), 9949,
     &                    mpi_comm_hycom, mpistat, mpierr)
          endif
        elseif (mnproc.eq.mnloc) then
          call MPI_SEND(tc,97,MTYPED,
     &                  idproc1(1), 9949,
     &                  mpi_comm_hycom, mpierr)
        endif
#elif defined(SHMEM)
        BARRIER
        if     (mnproc.eq.1) then
          mnloc = 1
          do mn= 2,ijpr
            call SHMEM_GETD(tcxc(2),tcxc(1),1,idproc1(mn))
            if     (tcxc(2).gt.tcxc(1)) then
              tcxc(1) = tcxc(2)
              mnloc   = mn
            endif
          enddo
          tcxc(2) = mnloc  !processor with the least comm. overhead
        endif
        if     (mnloc.ne.1) then
          call SHMEM_GETD(tc,tc,97,idproc1(mnloc))
        endif
        BARRIER
#endif
      endif
c
      call xcsync(flush_lp)
      if     (mnproc.eq.1) then
        write(lp,6000) mnloc,ijpr
        do i= 1,32
          if     (nc(i).ne.0) then
            if     (cc(i).ne.'      ') then
              write(lp,6100) cc(i),nc(i),tc(i),tc(i)/nc(i)
            else
              write(lp,6150)    i, nc(i),tc(i),tc(i)/nc(i)
            endif
          endif
        enddo !i
        write(lp,6100) 'xc****',1,tcxc(1),tcxc(1)
        do i= 33,97
          if     (nc(i).ne.0) then
            if     (cc(i).ne.'      ') then
              write(lp,6100) cc(i),nc(i),tc(i),tc(i)/nc(i)
            else
              write(lp,6150)    i, nc(i),tc(i),tc(i)/nc(i)
            endif
          endif
        enddo !i
        write(lp,6200)
      endif !mnproc.eq.1
      call xcsync(flush_lp)
c
c     reset timers to zero.
c
      do i= 1,97
        nc(i) = 0
        tc(i) = zero8
      enddo
      tcxc(1) = zero8
c
c     start a new total time measurement.
c
      call xctmr0(97)
      return
c
 6000 format(/ /
     &    3x,' timer statistics, processor',i5,' out of',i5 /
     &    3x,'-----------------------------------------------' /)
 6100 format(3x,a6,
     &   '   calls =',i9,
     &   '   time =',f11.5,
     &   '   time/call =',f14.8)
 6150 format(3x,'   #',i2,
     &   '   calls =',i9,
     &   '   time =',f11.5,
     &   '   time/call =',f14.8)
 6200 format(/ /)
      end subroutine xctmrp
#endif /* TIMER_ALLOUT:else */

      end module mod_xc
