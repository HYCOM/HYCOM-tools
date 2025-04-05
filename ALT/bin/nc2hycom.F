      PROGRAM FNC2
      IMPLICIT NONE
C
C  nc2hycom - Usage:  nc2hycom   fnc.a vname idm jdm [spval] [i1 j1] fhycom.a
C                     nc2hycom3d fnc.a vname idm jdm kdm [spval] [i1 j1] fhycom.a
C
C  Outputs a HYCOM ".a" copy of a netCDF file.
C
C  The input array is (1:idm,1:jdm), output is (i1:idm,j1:jdm)
C
C  fnc.a is assumed to be a netcdf file containing variable vname with
C   axis either T,Y,X or Y,X or for 3D version T,Z,Y,X or Z,Y,X
C
C  fhycom.a will contain (idm-i1+1)*(jdm-i1+1) 32-bit IEEE real values
C   for each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and 2.0**100 indicating a data void.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  COAPS (FSU),  October 2022.
C
      REAL*4, ALLOCATABLE :: A(:,:)
      REAL*4              :: PAD(4096)
      INTEGER IOS
      INTEGER      IARGC
      INTEGER      NARG
      CHARACTER*240 CARG
C
      LOGICAL      LSPVAL,L3D
      REAL*4       SPVAL
      INTEGER      IDM,JDM,KDM,I1,J1,L,NPAD
      CHARACTER*240 CFILE1,CFILEO,VNAME
C
C     READ ARGUMENTS.
C
      CALL GETARG(0,CARG)
      L = LEN_TRIM(CARG)
*     WRITE(6,"(4a)") TRIM(CARG),'"',CARG(L-4:L),'"'
      L3D = CARG(L-1:L).EQ.'3d'
C
      NARG = IARGC()
C
      IF     (L3D) THEN
        IF     (NARG.EQ.9) THEN
          CALL GETARG(1,CFILE1)
          CALL GETARG(2,VNAME)
          CALL GETARG(3,CARG)
          READ(CARG,*) IDM
          CALL GETARG(4,CARG)
          READ(CARG,*) JDM
          CALL GETARG(5,CARG)
          READ(CARG,*) KDM
          LSPVAL = .TRUE.
          CALL GETARG(6,CARG)
          READ(CARG,*) SPVAL
          CALL GETARG(7,CARG)
          READ(CARG,*) I1
          CALL GETARG(8,CARG)
          READ(CARG,*) J1
          CALL GETARG(9,CFILEO)
        ELSEIF (NARG.EQ.8) THEN
          CALL GETARG(1,CFILE1)
          CALL GETARG(2,VNAME)
          CALL GETARG(3,CARG)
          READ(CARG,*) IDM
          CALL GETARG(4,CARG)
          READ(CARG,*) JDM
          CALL GETARG(5,CARG)
          READ(CARG,*) KDM
          LSPVAL = .FALSE.
          CALL GETARG(6,CARG)
          READ(CARG,*) I1
          CALL GETARG(7,CARG)
          READ(CARG,*) J1
          CALL GETARG(8,CFILEO)
        ELSEIF (NARG.EQ.7) THEN
          CALL GETARG(1,CFILE1)
          CALL GETARG(2,VNAME)
          CALL GETARG(3,CARG)
          READ(CARG,*) IDM
          CALL GETARG(4,CARG)
          READ(CARG,*) JDM
          CALL GETARG(5,CARG)
          READ(CARG,*) KDM
          LSPVAL = .TRUE.
          CALL GETARG(6,CARG)
          READ(CARG,*) SPVAL
          I1 = 1
          J1 = 1
          CALL GETARG(7,CFILEO)
        ELSEIF (NARG.EQ.6) THEN
          CALL GETARG(1,CFILE1)
          CALL GETARG(2,VNAME)
          CALL GETARG(3,CARG)
          READ(CARG,*) IDM
          CALL GETARG(4,CARG)
          READ(CARG,*) JDM
          CALL GETARG(5,CARG)
          READ(CARG,*) KDM
          LSPVAL = .FALSE.
          I1 = 1
          J1 = 1
          CALL GETARG(6,CFILEO)
        ELSE
          WRITE(6,'(2a)')
     &    'Usage: nc2hycom3d fnc.a vname idm jdm kdm ',
     &                       '[spval] [i1 j1] fhycom.a'
          CALL EXIT(1)
        ENDIF
      ELSE
        KDM = 0
        IF     (NARG.EQ.8) THEN
          CALL GETARG(1,CFILE1)
          CALL GETARG(2,VNAME)
          CALL GETARG(3,CARG)
          READ(CARG,*) IDM
          CALL GETARG(4,CARG)
          READ(CARG,*) JDM
          LSPVAL = .TRUE.
          CALL GETARG(5,CARG)
          READ(CARG,*) SPVAL
          CALL GETARG(6,CARG)
          READ(CARG,*) I1
          CALL GETARG(7,CARG)
          READ(CARG,*) J1
          CALL GETARG(8,CFILEO)
        ELSEIF (NARG.EQ.7) THEN
          CALL GETARG(1,CFILE1)
          CALL GETARG(2,VNAME)
          CALL GETARG(3,CARG)
          READ(CARG,*) IDM
          CALL GETARG(4,CARG)
          READ(CARG,*) JDM
          LSPVAL = .FALSE.
          CALL GETARG(5,CARG)
          READ(CARG,*) I1
          CALL GETARG(6,CARG)
          READ(CARG,*) J1
          CALL GETARG(7,CFILEO)
        ELSEIF (NARG.EQ.6) THEN
          CALL GETARG(1,CFILE1)
          CALL GETARG(2,VNAME)
          CALL GETARG(3,CARG)
          READ(CARG,*) IDM
          CALL GETARG(4,CARG)
          READ(CARG,*) JDM
          LSPVAL = .TRUE.
          CALL GETARG(5,CARG)
          READ(CARG,*) SPVAL
          I1 = 1
          J1 = 1
          CALL GETARG(6,CFILEO)
        ELSEIF (NARG.EQ.5) THEN
          CALL GETARG(1,CFILE1)
          CALL GETARG(2,VNAME)
          CALL GETARG(3,CARG)
          READ(CARG,*) IDM
          CALL GETARG(4,CARG)
          READ(CARG,*) JDM
          LSPVAL = .FALSE.
          I1 = 1
          J1 = 1
          CALL GETARG(5,CFILEO)
        ELSE
          WRITE(6,*)
     &    'Usage: nc2hycom fnc.a vname idm jdm [spval] [i1 j1] fhycom.a'
          CALL EXIT(1)
        ENDIF
      ENDIF !l3d:else
C
      NPAD = 4096 - MOD((IDM-I1+1)*(JDM-J1+1),4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( A(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in nc2hycom: could not allocate ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
C
      CALL NC(A,IDM,JDM,KDM,PAD,NPAD, I1,J1, LSPVAL,SPVAL,
     &        VNAME,CFILE1,CFILEO)
      CALL EXIT(0)
      END
      SUBROUTINE NC(A,IDM,JDM,KDM,PAD,NPAD, I1,J1,
     &              LSPVAL,SPVAL, VNAME,CFILE1,CFILEO)
      IMPLICIT NONE
C
      REAL*4     SPVALH
      PARAMETER (SPVALH=2.0**100)
C
      CHARACTER*240 VNAME,CFILE1,CFILEO
      LOGICAL       LSPVAL
      INTEGER       IDM,JDM,KDM,NPAD,I1,J1
      REAL*4        SPVAL
      REAL*4        A(IDM,JDM),PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C
      DOUBLE PRECISION TIME
      CHARACTER*18 CASN
      INTEGER      I,J,K,MR,N,IOS,NRECL
      INTEGER      NCX,NCY,NCZ,NCT
      REAL*4       AMN,AMX
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
      IF     (.NOT.LSPVAL) THEN
        SPVAL = SPVALH
      ENDIF
C
      IF     (NPAD.EQ.0) THEN
        INQUIRE( IOLENGTH=NRECL) A(I1:IDM,J1:JDM)
      ELSE
        INQUIRE( IOLENGTH=NRECL) A(I1:IDM,J1:JDM),PAD
        PAD(:) = SPVALH
      ENDIF
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
      ENDIF
#else
      CALL ASNUNIT(11,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 11'
        write(6,*) 'ios = ',ios
        CALL EXIT(5)
      ENDIF
#endif
#endif
      OPEN(UNIT=11, FILE=CFILEO, FORM='UNFORMATTED', STATUS='NEW',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEO)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
      if     (kdm.eq.0) then
        call rd_dimen2(ncx,ncy,nct, cfile1,vname)
        if     (ncx.ne.idm .or. ncy.ne.jdm) then
          write(6,*) 'Error: netcdf dimensions are wrong'
          write(6,*) 'ncx,ncy =',ncx,ncy
          write(6,*) 'idm,jdm =',idm,jdm
          CALL EXIT(3)
        endif
      else
        call rd_dimen3(ncx,ncy,ncz,nct, cfile1,vname)
        if     (ncx.ne.idm .or. ncy.ne.jdm .or. ncz.ne.kdm) then
          write(6,*) 'Error: netcdf dimensions are wrong'
          write(6,*) 'ncx,ncy,ncz =',ncx,ncy,ncz
          write(6,*) 'idm,jdm,kdm =',idm,jdm,kdm
          CALL EXIT(3)
        endif
      endif !kdm is 0 or not
C
      MR = 0
      DO N= MIN(NCT,1),NCT
        DO K= 1,MAX(KDM,1)
          MR = MR + 1
          IF     (KDM.EQ.0) THEN
            call rd_field(  ncx,ncy,  n, a,time, vname,cfile1)
          ELSE
            call rd_field3d(ncx,ncy,k,n, a,time, vname,cfile1)
          ENDIF
          AMN =  SPVALH
          AMX = -SPVALH
          DO 210 J= J1,JDM
            DO 212 I= I1,IDM
              IF     (A(I,J).GE.-HUGE(AMN) .AND.
     &                A(I,J).LE. HUGE(AMN)      ) THEN
                IF     (A(I,J).NE.SPVAL) THEN
                  AMN = MIN( AMN, A(I,J) )
                  AMX = MAX( AMX, A(I,J) )
                ELSEIF (LSPVAL) THEN
                  A(I,J) = SPVALH !hycom data void values
                ENDIF
              ELSE
                A(I,J) = SPVALH  !replace NaN or -Inf or +Inf with data void
              ENDIF
  212       CONTINUE
  210     CONTINUE
#ifdef ENDIAN_IO
          CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
          IF     (NPAD.EQ.0) THEN
            WRITE(11,REC=MR,IOSTAT=IOS) A(I1:IDM,J1:JDM)
          ELSE
            WRITE(11,REC=MR,IOSTAT=IOS) A(I1:IDM,J1:JDM),PAD
          ENDIF
          IF     (NCT.GE.1) THEN
            WRITE(6,"(a,' =',f11.3,1p2e16.7)")
     &        trim(vname),TIME,AMN,AMX
          ELSE
            WRITE(6,"(a,' =',1p2e16.7)")
     &        trim(vname),AMN,AMX
          ENDIF
        ENDDO !k
      ENDDO !n
*     WRITE(6,*) 
*     WRITE(6,*) MR,' FIELDS PROCESSED'
*     WRITE(6,*) 
C
      CLOSE(11)
C
      RETURN
      END

      subroutine rd_field(n,m,irec,
     &                    t,time,
     &                    name_t,flnm_t)
      use netcdf   ! NetCDF fortran 90 interface
      implicit none
c
      character*(*)    name_t,flnm_t
      integer          n,m,irec
      double precision time
      real             t(n,m)
c
c  subroutine to read a scalar 2-d field.
c
c  subroutine arguments:
c       n,m     = horizontal grid dimensions.
c       irec    = time record to input (0-> no time axis)
c       flnm_t  = filename of the netCDF input.
c       name_t  =     name of the netCDF field.
c
c       t       = field
c       time    = netcdf time axis time
c
      character*(NF90_MAX_NAME)     :: cD
      integer                       :: ncFID,ncVID,ncNID,ncSTATUS
      integer                       :: ncDIDs(nf90_max_var_dims)
c
*     write(6,*) "irec  = ",irec
c
        call nchek('nf90_open-T',
     &              nf90_open(trim(flnm_t), nf90_nowrite, ncFID))
        if     (irec.eq.0) then
          time = 0.0
        else
          call nchek('nf90_inq_varid-'//trim(name_t),
     &                nf90_inq_varid(ncFID,trim(name_t), ncVID))
          call nchek('nf90_inq_variable(dimids)',
     &                nf90_inquire_variable(ncFID,  ncVID,
     &                                           dimids=ncDIDs(:3)))
c
          call nchek('nf90_inquire_dimension-3',
     &              nf90_inquire_dimension(ncFID, ncDIDs(3), name=cD))
          call nchek('nf90_inq_varid-Time',
     &                nf90_inq_varid(ncFID,cD,  ncVID))
          call nchek('nf90_get_var-Time',
     &                nf90_get_var(  ncFID,         ncVID, time,
     &                                                (/ irec /) ))
        endif !irec
*       write(6,*) "name  = ",trim(name_t)
*       write(6,*) "time  = ",time
        call nchek('nf90_inq_varid-'//trim(name_t),
     &              nf90_inq_varid(ncFID,trim(name_t), ncVID))
        if     (irec.eq.0) then
          call nchek('nf90_get_var-'//trim(name_t),
     &                nf90_get_var(  ncFID,              ncVID, t(:,:),
     &                                                 (/ 1,1 /) ))
        else
          call nchek('nf90_get_var-'//trim(name_t),
     &                nf90_get_var(  ncFID,              ncVID, t(:,:),
     &                                                 (/ 1,1,irec /) ))
        endif !irec
        call nchek("nf90_close",
     &              nf90_close(ncFID))
c
      return
      end

      subroutine rd_field3d(n,m,k,irec,
     &                      t,time,
     &                      name_t,flnm_t)
      use netcdf   ! NetCDF fortran 90 interface
      implicit none
c
      character*(*)    name_t,flnm_t
      integer          n,m,k,irec
      double precision time
      real             t(n,m)
c
c  subroutine to read a scalar 2-d field.
c
c  subroutine arguments:
c       n,m     = horizontal grid dimensions.
c       k       = layer to input
c       irec    = time record to input (0-> no time axis)
c       flnm_t  = filename of the netCDF input.
c       name_t  =     name of the netCDF field.
c
c       t       = field
c       time    = netcdf time axis time
c
      character*(NF90_MAX_NAME)     :: cD
      integer                       :: ncFID,ncVID,ncNID,ncSTATUS
      integer                       :: ncDIDs(nf90_max_var_dims)
c
*     write(6,*) "irec  = ",irec
c
        call nchek('nf90_open-T',
     &              nf90_open(trim(flnm_t), nf90_nowrite, ncFID))
        if     (irec.eq.0) then
          time = 0.0
        else
          call nchek('nf90_inq_varid-'//trim(name_t),
     &                nf90_inq_varid(ncFID,trim(name_t), ncVID))
          call nchek('nf90_inq_variable(dimids)',
     &                nf90_inquire_variable(ncFID,  ncVID,
     &                                           dimids=ncDIDs(:4)))
c
          call nchek('nf90_inquire_dimension-4',
     &              nf90_inquire_dimension(ncFID, ncDIDs(4), name=cD))
          call nchek('nf90_inq_varid-Time',
     &                nf90_inq_varid(ncFID,cD,  ncVID))
          call nchek('nf90_get_var-Time',
     &                nf90_get_var(  ncFID,         ncVID, time,
     &                                                (/ irec /) ))
        endif !irec
*       write(6,*) "name  = ",trim(name_t)
*       write(6,*) "time  = ",time
        call nchek('nf90_inq_varid-'//trim(name_t),
     &              nf90_inq_varid(ncFID,trim(name_t), ncVID))
        if     (irec.eq.0) then
          call nchek('nf90_get_var-'//trim(name_t),
     &                nf90_get_var(  ncFID,            ncVID, t(:,:),
     &                                               (/ 1,1,k /) ))
        else
          call nchek('nf90_get_var-'//trim(name_t),
     &                nf90_get_var(  ncFID,            ncVID, t(:,:),
     &                                               (/ 1,1,k,irec /) ))
        endif !irec
        call nchek("nf90_close",
     &              nf90_close(ncFID))
c
      return
      end

      subroutine rd_dimen2(xto,yto,tto, cfile,cname)
      use netcdf   ! NetCDF fortran 90 interface
      implicit none
c
      integer       xto,yto,tto
      character*(*) cfile,cname
c
c  subroutine to read model 2-D dimensions
c       xto,yto= horizontal dimensions of entire grid
c       tto    = total number of time samples
c
      character*(NF90_MAX_NAME) :: cD
      integer  i,ncFID,ncDID,ncVID,ncNID,ncDIDs(nf90_max_var_dims)
      integer  nx(4),ny(4)
c
        call nchek('nf90_open-rd_dimen2',
     &              nf90_open(trim(cfile), nf90_nowrite, ncFID))
c
        call nchek('nf90_inq_varid',
     &              nf90_inq_varid(ncFID, trim(cname), ncVID))
c
        call nchek('nf90_inq_variable(ndims)',
     &              nf90_inquire_variable(ncFID, ncVID, ndims=ncNID))
c
        if     (ncNID.eq.4) then
          write(6,'(/ 3a /)')
     &  'error - variable ',trim(cname),' has 4 dimensions try 3d'
          stop
        elseif (ncNID.ne.2 .and. ncNID.ne.3) then
          write(6,'(/ 3a /)')
     &  'error - variable ',trim(cname),' does not have 2-3 dimensions'
          stop
        endif
c
        call nchek('nf90_inq_variable(dimids)',
     &              nf90_inquire_variable(ncFID,  ncVID,
     &                                         dimids=ncDIDs(:ncNID)))
c
        call nchek('nf90_inquire_dimension-1',
     &              nf90_inquire_dimension(ncFID, ncDIDs(1), len=xto))
c
        call nchek('nf90_inquire_dimension-2',
     &              nf90_inquire_dimension(ncFID, ncDIDs(2), len=yto))
c
        if     (ncNID.eq.2) then
          tto = 0
        else
          call nchek('nf90_inquire_dimension-3',
     &              nf90_inquire_dimension(ncFID, ncDIDs(3), len=tto))
        endif
c
        call nchek("nf90_close",
     &              nf90_close(ncFID))
      return 
      end

      subroutine rd_dimen3(xto,yto,zto,tto, cfile,cname)
      use netcdf   ! NetCDF fortran 90 interface
      implicit none
c
      integer       xto,yto,zto,tto
      character*(*) cfile,cname
c
c  subroutine to read model #-D dimensions
c       xto,yto,zto = horizontal dimensions of entire grid
c       tto         = total number of time samples
c
      character*(NF90_MAX_NAME) :: cD
      integer  i,ncFID,ncDID,ncVID,ncNID,ncDIDs(nf90_max_var_dims)
      integer  nx(4),ny(4)
c
        call nchek('nf90_open-rd_dimen3',
     &              nf90_open(trim(cfile), nf90_nowrite, ncFID))
c
        call nchek('nf90_inq_varid',
     &              nf90_inq_varid(ncFID, trim(cname), ncVID))
c
        call nchek('nf90_inq_variable(ndims)',
     &              nf90_inquire_variable(ncFID, ncVID, ndims=ncNID))
c
        if     (ncNID.eq.2) then
          write(6,'(/ 3a /)')
     &  'error - variable ',trim(cname),' has 2 dimensions try 2d'
          stop
        elseif (ncNID.ne.3 .and. ncNID.ne.4) then
          write(6,'(/ 3a /)')
     &  'error - variable ',trim(cname),' does not have 3-4 dimensions'
          stop
        endif
c
        call nchek('nf90_inq_variable(dimids)',
     &              nf90_inquire_variable(ncFID,  ncVID,
     &                                         dimids=ncDIDs(:ncNID)))
c
        call nchek('nf90_inquire_dimension-1',
     &              nf90_inquire_dimension(ncFID, ncDIDs(1), len=xto))
c
        call nchek('nf90_inquire_dimension-2',
     &              nf90_inquire_dimension(ncFID, ncDIDs(2), len=yto))
c
        call nchek('nf90_inquire_dimension-3',
     &              nf90_inquire_dimension(ncFID, ncDIDs(3), len=zto))
c
        if     (ncNID.eq.3) then
          tto = 0
        else
          call nchek('nf90_inquire_dimension-4',
     &              nf90_inquire_dimension(ncFID, ncDIDs(4), len=tto))
        endif
c
        call nchek("nf90_close",
     &              nf90_close(ncFID))
      return 
      end

      subroutine nchek(cnf90,status)
      use netcdf   ! NetCDF fortran 90 interface
      implicit none
c
      character*(*), intent(in) :: cnf90
      integer,       intent(in) :: status
c
c     subroutine to handle NetCDF errors
c
*     if     (.TRUE. ) then !debug
      if     (.FALSE.) then !nodebug
        write(6,'(a)') trim(cnf90)
      endif

      if (status /= nf90_noerr) then
        write(6,'(/a)')   'error from NetCDF library'
        write(6,'(a/)')   trim(cnf90)
        write(6,'(a/)')   trim(nf90_strerror(status))
        stop
      end if
      end subroutine nchek
