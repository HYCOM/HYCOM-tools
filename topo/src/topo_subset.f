      PROGRAM TOPO_SUBREG
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      integer    mxpe1,mxpe2
      parameter (mxpe1=128,mxpe2=32768)
c
      INTEGER*1     RAWPAL(0:255,3),PPMPAL(3,0:255)
      INTEGER       IH,JH,K,NC
      CHARACTER*5   CIH,CJH
      CHARACTER*1   CNWLIN
c
      logical   subreg(mxpe2)
      integer   iipx(mxpe1,mxpe1),ispx(mxpe1,mxpe1)
      integer   jjpx(mxpe1),jspx(mxpe1)  ! always separable
      integer   ispt(mxpe2),ilpt(mxpe2),jspt(mxpe2),jlpt(mxpe2)
      integer   idm_in,ibig,jdm_in,jbig,nmpe,npe,mpe
c
      real      hmaxa,hmaxb,hmina,hminb
      character preambl(5)*79,cline*80
c
      integer ksub
      integer i,j
      integer if,il,jf,jl,kreg,nreg,nsea
      logical ldebug
c
c --- This program reads in a standard HYCOM depth file, a patch
c --- distibution file, and a list of rectangular subregions,
c --- and writes out a list of the associated patches and the
c --- corresponding 24-bit PPM image file with the patch locations
c --- marked.
c
c --- stdin (unit 5) should have:
c      
c      pixel subsample factor
c      number of sub-regions, -ve for debugging printout
c      if il jf jl  (extent of sub-region, one sub-region per line)
c
      character*1, allocatable :: ib(:)
      integer,     allocatable :: ipsum(:),ip(:,:),iop(:,:)
      real,        allocatable :: depths(:,:)
c
      call xcspmd  !input idm,jdm
c
c --- ksub
      read(5,*) ksub
      write(6,'(a,i2)') 'ksub =',ksub
      call zhflsh(6)
c
      if     (ksub.gt.0) then
        allocate( ib(    19+3*((idm+ksub-1)/ksub)*((jdm+ksub-1)/ksub)) )
      endif
      allocate( ipsum( idm+jdm) )
      allocate( ip(    0:idm+1,0:jdm+1) )
      allocate( iop(     idm,    jdm) )
      allocate( depths(  idm,    jdm) )
c
c --- acquire basin depths from unit 51.
c
      call zhopnc(51, 'regional.depth.b', 'formatted', 'old', 0)
      read (51,'(a79)') preambl
      read (51,'(a)')   cline
      close(unit=51)
      write(6,'(/(a))') preambl,cline(1:len_trim(cline))
      write(6,*)
      call zhflsh(6)
c
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
c
      call zaiost
      call zaiopf('regional.depth.a','old', 51)
      call zaiord(depths,iop,.false., hmina,hmaxa, 51)
      call zaiocl(51)
c
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        write(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b topography files not consistent:',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        call zhflsh(6)
        stop
      endif
c
c --- subregions
c
      do j= 1,jdm
        do i= 1,idm
          iop(i,j) = 0
        enddo
      enddo
      read(5,*) nreg
      ldebug = nreg.lt.0  !negative nreg for debugging printout
      nreg = abs(nreg)
      nsea = 0
      do kreg= 1,nreg
        read(5,*) if,il,jf,jl
        write(6,'(a,4i5)') 'include subregion: ',if,il,jf,jl
        nsea = nsea + (il-if+1)*(jl-jf+1)
        do j= jf,jl
          do i= if,il
            if     (depths(i,j).lt.2.0**99) then
              iop(i,j) = 1
            endif
          enddo
        enddo
      enddo
      write(6,'(a,i11,f8.2)') 'number and % of subreg points =',
     &                        nsea,100.0*(real(nsea)/real(idm*jdm))
      call zhflsh(6)
c
c --- land mask
c
      do j= 1,jdm
        do i= 1,idm
          if     (depths(i,j).lt.2.0**99) then
            ip(    i,j) = 1
            if     (iop(i,j).eq.0) then
              depths(i,j) = -1.0
            endif
          else
            ip(    i,j) = 0
            depths(i,j) = 0.0
          endif
        enddo
      enddo
c
c     patch distibution file on unit 21 (patch.input).
c
        call zhopnc(21, 'patch.input', 'formatted', 'old', 0)
        read( 21,'(/7i6/)')   nmpe,npe,mpe,idm_in,jdm_in,ibig,jbig
        do j= 1,mpe
          read( 21,'(12x,8i6)') (ispx(i,j),i=1,npe)
          read( 21,'(12x,8i6)') (iipx(i,j),i=1,npe)
        enddo
        read( 21,*)
        read( 21,'(12x,8i6)') (jspx(j),j=1,mpe)
        read( 21,'(12x,8i6)') (jjpx(j),j=1,mpe)
        close(21)
c
        k    = 0
        kreg = 0
        do j= 1,mpe
          do i= 1,npe
            if     (iipx(i,j).ne.0) then
              k = k + 1
              ispt(k) = ispx(i,j)
              ilpt(k) = ispx(i,j) + iipx(i,j) - 1
              jspt(k) = jspx(  j)
              jlpt(k) = jspx(  j) + jjpx(  j) - 1
c
              if     (ilpt(k).le.idm) then
                nsea = sum(iop(ispt(k):ilpt(k),    jspt(k):jlpt(k)))
              else
                nsea = sum(iop(ispt(k):idm,        jspt(k):jlpt(k))) +
     &                 sum(iop(      1:ilpt(k)-idm,jspt(k):jlpt(k)))
              endif
              subreg(k) = nsea.gt.0
              if     (subreg(k)) then
                kreg = kreg + (ilpt(k)-ispt(k)+1)*(jlpt(k)-jspt(k)+1)
                if     (ldebug) then
                  write(6,'(a,5i5,i10)') 'pe,if,il,jf,jl,nsea = ',
     &                                k,ispt(k),ilpt(k),jspt(k),jlpt(k),
     &                                nsea
                endif !ldebug
              endif
            endif
          enddo
        enddo
        if     (k.ne.nmpe) then
          write(6,'(/ a,i5,a,i5 /)')
     &      'error - there are ',k,'sea-tiles, but nmpe = ',nmpe
          call zhflsh(6)
          stop
        endif
        write(6,'(a,i11,f8.2)') 'number and % of tiled  points =',
     &                          kreg,100.0*(real(kreg)/real(idm*jdm))
        write(6,*)
        call zhflsh(6)
C
      IF     (KSUB.GT.0) THEN
C
C       READ IN PSUDO-COLOR PALLETTE.
C
        OPEN(UNIT=41,FILE='xbathy.pal',STATUS='OLD',
     +       FORM='UNFORMATTED',ACCESS='DIRECT',RECL=768)
        READ( 41,REC=1) RAWPAL
        CLOSE(41)
        DO 110 K= 1,3
          DO 112 I= 0,255
            PPMPAL(K,I) = RAWPAL(I,K)
  112     CONTINUE
  110   CONTINUE
C
C       CONVERT TO IMAGE.
C
        IF     (MAX((IDM+KSUB-1)/KSUB,(JDM+KSUB-1)/KSUB).LE. 999) THEN
          NC = 15 + 3*((IDM+KSUB-1)/KSUB)*((JDM+KSUB-1)/KSUB)
        ELSEIF (MAX((IDM+KSUB-1)/KSUB,(JDM+KSUB-1)/KSUB).LE.9999) THEN
          NC = 17 + 3*((IDM+KSUB-1)/KSUB)*((JDM+KSUB-1)/KSUB)
        ELSE
          NC = 19 + 3*((IDM+KSUB-1)/KSUB)*((JDM+KSUB-1)/KSUB)
        ENDIF
C
        CALL TOPIMG(IB(NC-3*((IDM+KSUB-1)/KSUB)*
     &                      ((JDM+KSUB-1)/KSUB)+1),PPMPAL,
     &              DEPTHS,IDM,JDM,
     &              SUBREG,ISPT,ILPT,JSPT,JLPT,NMPE,KSUB)
C
C       OUTPUT AS PPM IMAGE FILE.
C
        IH = (IDM+KSUB-1)/KSUB
        JH = (JDM+KSUB-1)/KSUB
C
        CNWLIN = CHAR(10)
        WRITE(CIH,4000) IH
        WRITE(CJH,4000) JH
C
        IB(1) = 'P'
        IB(2) = '6'
        IB(3) = CNWLIN
        IF     (MAX(IH,JH).LE. 999) THEN
          IB( 4) = CIH(3:3)
          IB( 5) = CIH(4:4)

          IB( 6) = CIH(5:5)
          IB( 7) = ' '
          IB( 8) = CJH(3:3)
          IB( 9) = CJH(4:4)
          IB(10) = CJH(5:5)
          IB(11) = CNWLIN
          IB(12) = '2'
          IB(13) = '5'
          IB(14) = '5'
          IB(15) = CNWLIN
        ELSEIF (MAX(IH,JH).LE.9999) THEN
          IB( 4) = CIH(2:2)
          IB( 5) = CIH(3:3)
          IB( 6) = CIH(4:4)
          IB( 7) = CIH(5:5)
          IB( 8) = ' '
          IB( 9) = CJH(2:2)
          IB(10) = CJH(3:3)
          IB(11) = CJH(4:4)
          IB(12) = CJH(5:5)
          IB(13) = CNWLIN
          IB(14) = '2'
          IB(15) = '5'
          IB(16) = '5'
          IB(17) = CNWLIN
        ELSE
          IB( 4) = CIH(1:1)
          IB( 5) = CIH(2:2)
          IB( 6) = CIH(3:3)
          IB( 7) = CIH(4:4)
          IB( 8) = CIH(5:5)
          IB( 9) = ' '
          IB(10) = CJH(1:1)
          IB(11) = CJH(2:2)
          IB(12) = CJH(3:3)
          IB(13) = CJH(4:4)
          IB(14) = CJH(5:5)
          IB(15) = CNWLIN
          IB(16) = '2'
          IB(17) = '5'
          IB(18) = '5'
          IB(19) = CNWLIN
        ENDIF
C
        WRITE(6,*) 'WRITING IMAGE TO patch.ppm'
        WRITE(6,*)
        CALL ZHFLSH(6)
        OPEN(UNIT=31,FILE='patch.ppm',STATUS='NEW',
     +       FORM='UNFORMATTED',ACCESS='DIRECT',RECL=NC)
        WRITE(31,REC=1) IB(1:NC)
        CLOSE(31)
      ENDIF !ksub>0
C
      WRITE(6,*) 'WRITING TILES TO patch.subreg'
      WRITE(6,*)
      CALL ZHFLSH(6)
      OPEN(UNIT=32,FILE='patch.subreg',STATUS='NEW',
     +     FORM='FORMATTED')
      DO K= 1,NMPE
        IF     (SUBREG(K)) THEN
          WRITE(32,'(I5.5)') K
        ENDIF
      ENDDO
      CLOSE(32)
C
 4000 FORMAT( I5 )
C     END OF TOPO_SUBREG
      END
      SUBROUTINE TOPIMG(IB,IP,HD,IH,JH,
     +                  SUBREG,ISPT,ILPT,JSPT,JLPT,NMPE,KSUB)
C
      INTEGER       IH,JH,NMPE,KSUB
      LOGICAL       SUBREG(NMPE)
      INTEGER       ISPT(NMPE),ILPT(NMPE),JSPT(NMPE),JLPT(NMPE)
      REAL          HD(IH,JH)
      CHARACTER*(1) IB(3,(IH+KSUB-1)/KSUB,(JH+KSUB-1)/KSUB),IP(3,0:255)
C
C**********
C*
C 1)  CONVERT TOPOGRAPHY TO AN 8-BIT IMAGE.
C*
C**********
C
      LOGICAL LLAND
C
      WRITE(6,*) 
      WRITE(6,*) 'IP(:,47) = ',(ICHAR(IP(K,47)),K=1,3)
      WRITE(6,*) 'IP(:,48) = ',(ICHAR(IP(K,48)),K=1,3)
      WRITE(6,*) 'IP(:,49) = ',(ICHAR(IP(K,49)),K=1,3)
      WRITE(6,*) 'IP(:,50) = ',(ICHAR(IP(K,50)),K=1,3)
      WRITE(6,*) 
      CALL ZHFLSH(6)
      JB = (JH+KSUB-1)/KSUB
      N2 = 0
      N4 = 0
      N5 = 0
      D2 =  1.0/50.0
      DO J= 1,JH,KSUB
        DO I= 1,IH,KSUB
          IF     (HD(I,J).LE.  -1.0) THEN
            IB(1,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(1,49)
            IB(2,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(2,49)
            IB(3,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(3,49)
            N4 = N4 + 1
          ELSEIF (HD(I,J).LE.   0.0) THEN
            IB(1,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(1,50)
            IB(2,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(2,50)
            IB(3,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(3,50)
            N2 = N2 + 1
          ELSE
            II = MIN(248, 100 + NINT((HD(I,J)      )*D2) )
            IB(1,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(1,II)
            IB(2,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(2,II)
            IB(3,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(3,II)
            N5 = N5 + 1
          ENDIF
        ENDDO
      ENDDO
      WRITE(6,*)
      WRITE(6,*) 'TOPIMG STATISTICS:'
      WRITE(6,*) '           LAND POINTS = ',N2
      WRITE(6,*) '            SEA POINTS = ',N5+N4
      WRITE(6,*) '     NEEDED SEA POINTS = ',N5
      WRITE(6,*)
      CALL ZHFLSH(6)
C
C     MARK ALL LAND IN PATCHES
C
      DO N= 1,NMPE
        IF     (SUBREG(N)) THEN
          DO J= 1,JH,KSUB
            IF     (J.GE.JSPT(N) .AND. J.LE.JLPT(N)) THEN
              DO I= 1,IH,KSUB
                IF     (I.GE.ISPT(N) .AND. I.LE.ILPT(N) .AND.
     &                 HD(I,J).EQ.0.0) THEN
                  IB(1,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(1,48)
                  IB(2,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(2,48)
                  IB(3,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(3,48)
                ENDIF
              ENDDO
            ENDIF
          ENDDO
        ENDIF
      ENDDO
C
C     MARK PATCH BOUNDARIES.
C
      DO N= 1,NMPE
        IF     (SUBREG(N)) THEN
          I=ISPT(N)
          DO J= JSPT(N),JLPT(N)
            IB(1,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(1,47)
            IB(2,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(2,47)
            IB(3,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(3,47)
          ENDDO
          I=ILPT(N)
          DO J= JSPT(N),JLPT(N)
            IB(1,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(1,47)
            IB(2,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(2,47)
            IB(3,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(3,47)
          ENDDO
          I=ILPT(N)-KSUB+1
          DO J= JSPT(N),JLPT(N)
            IB(1,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(1,47)
            IB(2,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(2,47)
            IB(3,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(3,47)
          ENDDO
          J= JSPT(N)
          DO I= ISPT(N),ILPT(N)
            IB(1,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(1,47)
            IB(2,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(2,47)
            IB(3,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(3,47)
          ENDDO
          J= JLPT(N)
          DO I= ISPT(N),ILPT(N)
            IB(1,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(1,47)
            IB(2,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(2,47)
            IB(3,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(3,47)
          ENDDO
          J= JLPT(N)-KSUB+1
          DO I= ISPT(N),ILPT(N)
            IB(1,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(1,47)
            IB(2,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(2,47)
            IB(3,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(3,47)
          ENDDO
        ENDIF
      ENDDO
      RETURN
C     END OF TOPIMG.
      END
