      module mod_dimensions
c
c     module needed for CCSM3 integration
c
c-----------------------------------------------------------------------------
c --- START OF REGION AND TILING SPECIFIC PARAMETERS
c --- See: README.src.newregion for more details.
c
C==================================================================
C  Version for Dynamic Allocation Diagnostics
C   San Moore, QinetiQ,  July 2009
C------------------------------------------------------------------
c --- itdm  = total grid dimension in i direction
c --- jtdm  = total grid dimension in j direction
c --- kdm   =       grid dimension in k direction
c
      integer    itdm,jtdm,kdm  
c
c --- iqr   = maximum number of tiles in i direction
c --- jqr   = maximum number of tiles in j direction
      integer    iqr,jqr
      parameter (iqr=1,jqr=90)  ! multiple tiles (TYPE=ompi or mpi or shmem)
c
c --- idm   = maximum single tile grid dimension in i direction
c --- jdm   = maximum single tile grid dimension in j direction
c
      integer    idm,jdm
c
c ---   END OF REGION AND TILING SPECIFIC PARAMETERS
c-----------------------------------------------------------------------------
c
c --- halo size
c
      integer    nbdy
      parameter (nbdy=1)
c
c --- actual extent of this tile is (i0+1:i0+ii,j0+1:j0+jj,1:kk)
c
      integer      i0,j0,ii,jj
      common/dimi/ i0,j0,ii,jj
      save  /dimi/
c
c --- ijqr  = maximum total number of active tiles (= ipr*jpr)
c
      integer    ijqr
      parameter (ijqr=iqr*jqr)
c
c --- line printer unit (stdout)
c
      integer        lp
      common/linepr/ lp
      save  /linepr/
c-----------------------------------------------------------------------------
      end module mod_dimensions
