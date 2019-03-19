      PROGRAM HYCOM_FILL
      IMPLICIT NONE
C
C  hycom_fill - Usage:  hycom_fill fin.a fmask.a idm jdm fout.a [npass]
C
C                 Fills the fields in fin based on the first field from fmask.
C
C  fin.a and fmask.a are assumed to contain idm*jdm 32-bit IEEE real values
C   for each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C  if fmask.a is "NONE" or "NOARCTIC" then it is not read in and the
C   mask array has no data voids. if it is "NOARCTIC" then the domain
C   is not assumed to be a tripole grid.
C
C  the fields in fout.a will be copies of those in fin.a except that
C   they will have the same data void locations as in the first field
C   of fmask.a, with data being filled (extrapolated) to all locations
C   that are voids in fin.a but are not voids in fmask.a.
C  if npass is present the landfill will minimize the number of passes
C   (by filling voids in fmask and masking them at the end) and will
C   do no more than npass passes in total.
C
C  the 1st array dimension is assumed periodic, which will work for
C   closed domains with fmask(idm,:) as all data voids.
C  the 2nd array dimension is assumed closed, unless fmask(:,jdm)
C   contains non-voids when an Arctic dipole is assumed at j=jdm.
C
C  all fields in fout.a are assumed to have the same set of data voids.
C
C  note that any standard hycom bathymetry file has land masked, and is
C   therefore often a good choice for fmask.a.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  August 2002.
C
      REAL*4, ALLOCATABLE :: A(:,:),B(:,:)
      REAL*4              :: PAD(4096)
      INTEGER IOS
      INTEGER      IARGC
      INTEGER      NARG
      CHARACTER*240 CARG
C
      INTEGER      IDM,JDM,NPAD,NPASS
      CHARACTER*240 CFILE1,CFILE2,CFILEO
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.5) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CFILE2)
        CALL GETARG(3,CARG)
        READ(CARG,*) IDM
        CALL GETARG(4,CARG)
        READ(CARG,*) JDM
        CALL GETARG(5,CFILEO)
        NPASS = 0
      ELSEIF (NARG.EQ.6) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CFILE2)
        CALL GETARG(3,CARG)
        READ(CARG,*) IDM
        CALL GETARG(4,CARG)
        READ(CARG,*) JDM
        CALL GETARG(5,CFILEO)
        CALL GETARG(6,CARG)
        READ(CARG,*) NPASS
      ELSE
        WRITE(6,*)
     &    'Usage: hycom_fill fin.a fmask.a idm jdm fout.a [npass]'
        CALL EXIT(1)
      ENDIF
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( A(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_expr: could not allocate 1st ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( B(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_expr: could not allocate 2nd ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
C
      CALL FILL(A,B,IDM,JDM,PAD,NPAD,NPASS, CFILE1,CFILE2,CFILEO)
      CALL EXIT(0)
      END
      SUBROUTINE FILL(A,B,IDM,JDM,PAD,NPAD,NPASS, CFILE1,CFILE2,CFILEO)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE1,CFILE2,CFILEO
      INTEGER      IDM,JDM,NPAD,NPASS
      REAL*4       A(IDM,JDM),B(IDM,JDM),PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      CHARACTER*18 CASN
      LOGICAL      LARCTIC
      INTEGER      I,J,K,IOS,NRECL
      REAL*4       AMN,AMX
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
      INQUIRE( IOLENGTH=NRECL) A,PAD
#ifdef CRAY
#ifdef t3e
      IF     (MOD(NRECL,4096).EQ.0) THEN
        WRITE(CASN,8000) NRECL/4096
 8000   FORMAT('-F cachea:',I4.4,':1:0')
        IU8 = 11
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          write(6,*) 'Error: can''t asnunit 11'
          write(6,*) 'ios  = ',ios8
          write(6,*) 'casn = ',casn
          CALL EXIT(5)
        ENDIF
        IU8 = 12
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          write(6,*) 'Error: can''t asnunit 12'
          write(6,*) 'ios  = ',ios8
          write(6,*) 'casn = ',casn
          CALL EXIT(5)
        ENDIF
        IU8 = 21
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          write(6,*) 'Error: can''t asnunit 21'
          write(6,*) 'ios  = ',ios8
          write(6,*) 'casn = ',casn
          CALL EXIT(5)
        ENDIF
      ENDIF
#else
      CALL ASNUNIT(11,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 11'
        write(6,*) 'ios = ',ios
        CALL EXIT(5)
      ENDIF
      CALL ASNUNIT(12,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 12'
        write(6,*) 'ios = ',ios
        CALL EXIT(5)
      ENDIF
      CALL ASNUNIT(21,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 21'
        write(6,*) 'ios = ',ios
        CALL EXIT(5)
      ENDIF
#endif
#endif
      OPEN(UNIT=11, FILE=CFILE1, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILE1)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
      IF     (CFILE2.NE."NONE" .AND. CFILE2.NE."NOARCTIC") THEN
      OPEN(UNIT=12, FILE=CFILE2, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILE2)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
      ENDIF !.not.none
      OPEN(UNIT=21, FILE=CFILEO, FORM='UNFORMATTED', STATUS='NEW',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEO)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
c
      IF     (CFILE2.NE."NONE" .AND. CFILE2.NE."NOARCTIC") THEN
        READ(12,REC=1,IOSTAT=IOS) B
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(B,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'can''t read ',TRIM(CFILE2)
          CALL EXIT(4)
        ENDIF
      ELSE
        B(:,:) = 0.0  !no data voids
      ENDIF
C
      IF     (CFILE2.EQ."NOARCTIC") THEN
        LARCTIC = .FALSE.
      ELSE
        LARCTIC = MINVAL(B(:,JDM)).NE.SPVAL .AND.
     &            MAXVAL(B(:,JDM)).EQ.SPVAL
      ENDIF
      write(6,'(a,l1)') 'larctic = ',larctic
C
      DO 110 K= 1,9999
        READ(11,REC=K,IOSTAT=IOS) A
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          IF     (K.EQ.1) THEN
            WRITE(6,*) 'can''t read ',TRIM(CFILE1)
            CALL EXIT(4)
          ELSE
            GOTO 1110
          ENDIF
        ENDIF
C
        IF     (NPASS.EQ.0) THEN
          CALL LANDFILL(   A,B,IDM,JDM,        LARCTIC)
        ELSE
          CALL LANDFILL_NP(A,B,IDM,JDM, NPASS, LARCTIC)
        ENDIF
C
        AMN =  SPVAL
        AMX = -SPVAL
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (B(I,J).EQ.SPVAL) THEN
              A(I,J) = SPVAL
            ELSE
#ifdef sun
              IF     (IR_ISNAN(A(I,J)).NE.1) THEN
                AMN = MIN( AMN, A(I,J) )
                AMX = MAX( AMX, A(I,J) )
              ENDIF
#else
              AMN = MIN( AMN, A(I,J) )
              AMX = MAX( AMX, A(I,J) )
#endif
            ENDIF
          ENDDO
        ENDDO
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        WRITE(21,REC=K,IOSTAT=IOS) A
        WRITE(6,'(a,1p2g16.8)')
     &     'min, max = ',AMN,AMX
  110 CONTINUE
 1110 CONTINUE
      WRITE(6,*) 
      WRITE(6,*) K-1,' FIELDS PROCESSED'
      WRITE(6,*) 
C
      CLOSE(11)
      CLOSE(21)
C
      RETURN
      END
      subroutine landfill(a,amsk,m,n, larctic)
      implicit none
c
      logical larctic
      integer m,n
      real    a(m,n),amsk(m,n)
c
c --- extrapolate to fill voids in a that are not voids in amsk, using
c --- multiple pass's of a 9-point smoother based extrapolation scheme.
c --- do not extrapolate into points that are masked in amsk.
c
      real*4     spval
      parameter (spval=2.0**100)
c
      integer, allocatable :: mm(:,:,:),mmsumj(:)
c
      integer landa,landm,lande
      integer i,ii,ip0,ip1,ipass,j,jj,ki,kj,nup,nall,na
*     integer nrow
      real    sa,ss
c
      real s(-1:1,-1:1)
      data s / 1.0,2.0,1.0, 2.0,4.0,2.0, 1.0,2.0,1.0 /
c
      allocate( mm(0:m+1,0:n+1,0:1), mmsumj(n) )
c
      mm( : , : ,0) = 0
      do j= 1,n
        mmsumj(j) = 0
        do i= 1,m
          if     (   a(i,j).ne.spval) then
            mm(i,j,0) = 1
          elseif (amsk(i,j).ne.spval) then
            mmsumj(j) = mmsumj(j) + 1  !count of to be filled points
          endif
        enddo
      enddo

      if     (larctic) then
        na = n-1
      else
        na = n
      endif

      landa = 0
      landm = 0
      lande = 0
      do j= 1,na
        do i= 1,m
          if     (a(i,j).eq.spval) then
            landa = landa + 1
            if (amsk(i,j).ne.spval) then
              lande = lande + 1
            endif
          endif
          if     (amsk(i,j).eq.spval) then
            landm = landm + 1
          endif
        enddo
      enddo
      write(6,'(a,4i9)')
     &  '   pts,land = ',m*na,landa,landm,lande
c
      nall = 0
      do ipass= 1,99999
        ip0 = mod(ipass+1,2)
        ip1 = mod(ipass,  2)
        mm(:,:,ip1) = mm(:,:,ip0)
        nup  = 0
        do j= 1,na
          if     (mmsumj(j).ne.0) then
          do i= 1,m
            if     (mm(i,j,ip0).eq.0 .and. amsk(i,j).ne.spval) then
              sa = 0.0
              ss = 0.0
              do kj= -1,1
                jj = j+kj
                do ki= -1,1
                  ii = i+ki
                  if     (ii.eq.  0) then
                    ii = m  !periodic wrap
                  elseif (ii.eq.m+1) then
                    ii = 1  !periodic wrap
                  endif
                  if     (mm(ii,jj,ip0).eq.1) then
                    sa = sa + s(ki,kj)*a(ii,jj)
                    ss = ss + s(ki,kj)
                  endif
                enddo
              enddo
              if     (ss.ne.0.0) then
                a( i,j)     = sa/ss
                mm(i,j,ip1) = 1
                mmsumj(j)   = mmsumj(j) - 1
                nup         = nup + 1
*               if     (mod(nup,1000).eq.1) then
*                 write(6,'(a,2i5,f5.1,f10.3)') 
*    &              '   i,j,ss,a = ',i,j,ss,a(i,j)
*               endif
              endif
            endif
          enddo !i
          endif !mmsumj(j)
        enddo !j
        if     (larctic) then
          mmsumj(n) = mmsumj(n-1)
          do i= 1,m
            ii = m-mod(i-1,m)
             a(i,n)     =  a(ii,n-1)
            mm(i,n,ip1) = mm(ii,n-1,ip1)
          enddo !i
        endif
        nall = nall + nup
        write(6,'(a,i4,a,i6,a,a,i9,a)') 'landfill: pass',ipass,
     &                                  ' filled in',nup,' points',
     &                                  ' (total',nall,')'
*       nrow = 0
*       do j= 1,n
*         if     (mmsumj(j).eq.0) then
*           nrow = nrow + 1
*         endif
*       enddo !j
*       write(6,'(a,i4,a,i6,a,a,i9,a)') 'landfill: pass',ipass,
*    &                                  ' total of',nrow,' filled rows',
*    &                                  ' (left ',sum(mmsumj(1:na)),')'
        if     (nup.eq.0) then
          exit
        endif
      enddo  ! ipass
      write(6,*)
c
      deallocate( mm )
c
      return
      end
      subroutine landfill_np(a,amsk,m,n, npass, larctic)
      implicit none
c
      logical larctic
      integer m,n,npass
      real    a(m,n),amsk(m,n)
c
c --- extrapolate to fill voids in a that are not voids in amsk, using
c --- multiple pass's of a 9-point smoother based extrapolation scheme.
c --- extrapolate into points that are masked in amsk, but do no more
c --- than npass passes.
c
      real*4     spval
      parameter (spval=2.0**100)
c
      integer, allocatable :: mm(:,:,:),mmsumj(:)
c
      integer landa,landm,lande
      integer i,ii,ip0,ip1,ipass,j,jj,ki,kj,nup,nall,na
*     integer nrow
      real    sa,ss
c
      real s(-1:1,-1:1)
      data s / 1.0,2.0,1.0, 2.0,4.0,2.0, 1.0,2.0,1.0 /
c
      allocate( mm(0:m+1,0:n+1,0:1), mmsumj(n) )
c
      mm( : , : ,0) = 0
      do j= 1,n
        mmsumj(j) = 0
        do i= 1,m
          if     (   a(i,j).ne.spval) then
            mm(i,j,0) = 1
          else
            mmsumj(j) = mmsumj(j) + 1  !count of unfilled land points
          endif
        enddo
      enddo

      if     (larctic) then
        na = n-1
      else
        na = n
      endif

      landa = 0
      landm = 0
      lande = 0
      do j= 1,na
        do i= 1,m
          if     (a(i,j).eq.spval) then
            landa = landa + 1
            if (amsk(i,j).ne.spval) then
              lande = lande + 1
            endif
          endif
          if     (amsk(i,j).eq.spval) then
            landm = landm + 1
          endif
        enddo
      enddo
      write(6,'(a,4i9)')
     &  '   pts,land = ',m*na,landa,landm,lande
c
      nall = 0
      do ipass= 1,npass
        ip0 = mod(ipass+1,2)
        ip1 = mod(ipass,  2)
        mm(:,:,ip1) = mm(:,:,ip0)
        nup  = 0
        do j= 1,na
          if     (mmsumj(j).ne.0) then
          do i= 1,m
            if     (mm(i,j,ip0).eq.0) then
              sa = 0.0
              ss = 0.0
              do kj= -1,1
                jj = j+kj
                do ki= -1,1
                  ii = i+ki
                  if     (ii.eq.  0) then
                    ii = m  !periodic wrap
                  elseif (ii.eq.m+1) then
                    ii = 1  !periodic wrap
                  endif
                  if     (mm(ii,jj,ip0).eq.1) then
                    sa = sa + s(ki,kj)*a(ii,jj)
                    ss = ss + s(ki,kj)
                  endif
                enddo !ki
              enddo !kj
              if     (ss.ne.0.0) then
                a( i,j)     = sa/ss
                mm(i,j,ip1) = 1
                mmsumj(j)   = mmsumj(j) - 1
                if     (amsk(i,j).ne.spval) then
                  nup       = nup + 1
*                 if     (mod(nup,1000).eq.1) then
*                   write(6,'(a,2i5,f5.1,f10.3)') 
*    &                '   i,j,ss,a = ',i,j,ss,a(i,j)
*                 endif
                endif
              endif
            endif
          enddo !i
          endif !mmsumj(j)
        enddo !j
        if     (larctic) then
          mmsumj(n) = mmsumj(n-1)
          do i= 1,m
            ii = m-mod(i-1,m)
             a(i,n)     =  a(ii,n-1)
            mm(i,n,ip1) = mm(ii,n-1,ip1)
          enddo !i
        endif
        nall = nall + nup
        write(6,'(a,i4,a,i6,a,a,i9,a)') 'landfill: pass',ipass,
     &                                  ' filled in',nup,' points',
     &                                  ' (total',nall,')'
*       nrow = 0
*       do j= 1,n
*         if     (mmsumj(j).eq.0) then
*           nrow = nrow + 1
*         endif
*       enddo !j
*       write(6,'(a,i4,a,i6,a,a,i9,a)') 'landfill: pass',ipass,
*    &                                  ' total of',nrow,' filled rows',
*    &                                  ' (left ',sum(mmsumj(1:na)),')'
        if     (nup.eq.0) then
          exit
        endif
      enddo  ! ipass
      write(6,*)
c
      deallocate( mm )
c
      return
      end
