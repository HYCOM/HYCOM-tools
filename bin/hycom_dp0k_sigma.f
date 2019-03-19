      PROGRAM HYCOM_DP0K_SIGMA
      IMPLICIT NONE
C
C --- BUILD A TABLE OF DP0K VS SIGMA
C
      INTEGER KDM
      REAL    DP00,DP00X,DP00F,DP0K(9999),SIGMA(9999)
C
      INTEGER K
      REAL    DP0KF,Z
C
C --- INPUT, FROM STDIN
C
C --- 'kdm   ' = number of layers
C --- 'dp00  ' = z-level spacing minimum thickness (m)
C --- 'dp00x ' = z-level spacing maximum thickness (m)
C --- 'dp00f ' = z-level spacing stretching factor (1.0=const.z)
C --- 'sigma ' = target layer densities; k=1,kdm
C
C --- or, in place of 'dp00','dp00x','dp00f', specify:
C --- 'dp0k  ' = layer k z-level spacing minimum thickness (m); k=1,kdm
C
      CALL BLKINI(KDM,   'kdm   ')
      CALL BLKINR2(DP00,K,
     &                   'dp00  ','("# ",a6," =",f10.4," m")',
     &                   'dp0k  ','("# ",a6," =",f10.4," m")' )
      IF     (K.EQ.1) THEN
        CALL BLKINR(DP00X, 'dp00x ','("# ",a6," =",f10.4," m")')
        CALL BLKINR(DP00F, 'dp00f ','("# ",a6," =",f10.4,"  ")')
      ELSE
        DP0K(1) = DP00
        DP00    = -1.0  !signal that dp00 is not input
        DO K=2,KDM
          CALL BLKINR(DP0K(K), 'dp0k  ','("# ",a6," =",f10.4," m")')
        ENDDO !k
      ENDIF
      DO K=2,KDM
        CALL BLKINR(SIGMA(K), 'sigma ','("# ",a6," =",f10.4," sig")')
      ENDDO !k
C
      IF     (DP00.GT.0.0) THEN
        DP0KF=1.0/DP00F
        DO K=1,KDM
          DP0KF=DP0KF*DP00F
          DP0K(K)=MIN(DP00*DP0KF,DP00X)
        ENDDO !k
      ENDIF
C
C     DEPTH VS SIGMA
C
      WRITE(6,'(a)') '#'
      WRITE(6,'(a)') '#  Z-DEPTH     SIGMA'
      Z = 0.0
      DO K=1,KDM
        WRITE(6,'(F10.4,F10.3)') Z,SIGMA(K)
        Z = Z + DP0K(K)
        WRITE(6,'(F10.4,F10.3)') Z,SIGMA(K)
      ENDDO !k
      END

      subroutine blkinr(rvar,cvar,cfmt)
      implicit none
c
      real      rvar
      character cvar*6,cfmt*(*)
c
c     read in one real value
c
      character*6 cvarin
c
      read(5,*) rvar,cvarin
      write(6,cfmt) cvarin,rvar
c
      if     (cvar.ne.cvarin) then
        write(6,*) 
        write(6,*) 'error in blkinr - input ',cvarin,
     &                      ' but should be ',cvar
        write(6,*) 
        stop
      endif
      return
      end
      subroutine blkinr2(rvar,nvar,cvar1,cfmt1,cvar2,cfmt2)
      implicit none
c
      real      rvar
      integer   nvar
      character cvar1*6,cvar2*6,cfmt1*(*),cfmt2*(*)
c
c     read in one real value from stdin,
c     identified as either cvar1 (return nvar=1) or cvar2 (return nvar=2)
c
      character*6 cvarin
c
      read(5,*) rvar,cvarin
c
      if     (cvar1.eq.cvarin) then
        nvar = 1
        write(6,cfmt1) cvarin,rvar
      elseif (cvar2.eq.cvarin) then
        nvar = 2
        write(6,cfmt2) cvarin,rvar
      else
        write(6,*) 
        write(6,*) 'error in blkinr2 - input ',cvarin,
     +                      ' but should be ',cvar1,' or ',cvar2
        write(6,*) 
        stop
      endif
      return
      end
      subroutine blkini(ivar,cvar)
      implicit none
c
      integer     ivar
      character*6 cvar
c
c     read in one integer value
c
      character*6 cvarin
c
      read(5,*) ivar,cvarin
      write(6,6000) cvarin,ivar
c
      if     (cvar.ne.cvarin) then
        write(6,*) 
        write(6,*) 'error in blkini - input ',cvarin,
     &                      ' but should be ',cvar
        write(6,*) 
        stop
      endif
      return
 6000 format('# ',a6,' =',i6)
      end
