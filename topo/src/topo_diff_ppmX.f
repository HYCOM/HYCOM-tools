      PROGRAM TOP_DIFF_PPMX
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      INTEGER*1     RAWPAL(0:255,3),PPMPAL(3,0:255)
      INTEGER       IH,JH,K,NC
      CHARACTER*5   CIH,CJH
      CHARACTER*1   CNWLIN
c
      real      hmaxa,hmaxb,hmina,hminb
      character preambl(5)*79,cline*80
c
      logical lexist
      integer i,j
      integer ksub
c
c --- This program reads in two standard HYCOM depth files, and 
c --- and writes out the a 24-bit PPM image file with the differences
c --- marked.
c
c --- the subsample factor, ksub, is read from stdin
c
      character*1, allocatable :: ib(:)
      integer,     allocatable :: iop(:,:)
      real,        allocatable :: depths_1(:,:),depths_2(:,:)
c
      read(5,*) ksub
      write(6,'(a,i3)') 'ksub =',ksub
      call zhflsh(6)
c
      call xcspmd  !input idm,jdm
      allocate( ib(    19+3*((idm+ksub-1)/ksub)*((jdm+ksub-1)/ksub)) )
      allocate( iop(     idm,    jdm) )
      allocate( depths_1(idm,    jdm) )
      allocate( depths_2(idm,    jdm) )
      ib(:) = char(0)
c
c --- acquire 1st basin depths from unit 51.
c
      call zhopen(51, 'formatted', 'old', 0)
      read (51,'(a79)') preambl
      read (51,'(a)')   cline
      close(unit=51)
      write(6,'(/(a))') preambl,cline(1:len_trim(cline))
c
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
c
      call zaiost
      call zaiopn('old', 51)
      call zaiord(depths_1,iop,.false., hmina,hmaxa, 51)
      call zaiocl(51)
c
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        write(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b 1st topography files not consistent:',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        call zhflsh(6)
        stop
      endif
c
      call zhopen(61, 'formatted', 'old', 0)
      read (61,'(a79)') preambl
      read (61,'(a)')   cline
      close(unit=61)
      write(6,'(/(a))') preambl,cline(1:len_trim(cline))
c
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
c
      call zaiost
      call zaiopn('old', 61)
      call zaiord(depths_2,iop,.false., hmina,hmaxa, 61)
      call zaiocl(61)
c
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        write(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b 2nd topography files not consistent:',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        call zhflsh(6)
        stop
      endif
c
c --- zero land
c
      do j= 1,jdm
        do i= 1,idm
          if     (depths_1(i,j).gt.2.0**99) then
            depths_1(i,j) = 0.0
          endif
          if     (depths_2(i,j).gt.2.0**99) then
            depths_2(i,j) = 0.0
          endif
        enddo
      enddo
C
C     READ IN PSUDO-COLOR PALLETTE.
C
      OPEN(UNIT=41,FILE='xbathy.pal',STATUS='OLD',
     +     FORM='UNFORMATTED',ACCESS='DIRECT',RECL=768)
      READ( 41,REC=1) RAWPAL
      CLOSE(41)
      DO 110 K= 1,3
        DO 112 I= 0,255
          PPMPAL(K,I) = RAWPAL(I,K)
  112   CONTINUE
  110 CONTINUE
C
C     CONVERT TO IMAGE.
C
      IF     (MAX((IDM+KSUB-1)/KSUB,(JDM+KSUB-1)/KSUB).LE. 999) THEN
        NC = 15 + 3*((IDM+KSUB-1)/KSUB)*((JDM+KSUB-1)/KSUB)
      ELSEIF (MAX((IDM+KSUB-1)/KSUB,(JDM+KSUB-1)/KSUB).LE.9999) THEN
        NC = 17 + 3*((IDM+KSUB-1)/KSUB)*((JDM+KSUB-1)/KSUB)
      ELSE
        NC = 19 + 3*((IDM+KSUB-1)/KSUB)*((JDM+KSUB-1)/KSUB)
      ENDIF
C
      WRITE(6,*) 
      WRITE(6,*) "NC = ",NC,NC-3*((IDM+KSUB-1)/KSUB)*
     &                           ((JDM+KSUB-1)/KSUB)+1
      CALL ZHFLSH(6)
C
      CALL TOPIMG(IB(NC-3*((IDM+KSUB-1)/KSUB)*
     &                    ((JDM+KSUB-1)/KSUB)+1),char(PPMPAL),
     &            DEPTHS_1,DEPTHS_2,IDM,JDM, KSUB)
C
C     OUTPUT AS PPM IMAGE FILE.
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
      WRITE(6,*) 'WRITING IMAGE TO fort.31'
      CALL ZHFLSH(6)
      OPEN(UNIT=31,FILE='fort.31',STATUS='NEW',
     +     FORM='UNFORMATTED',ACCESS='DIRECT',RECL=NC)
      WRITE(31,REC=1) IB(1:NC)
      CLOSE(31)
C
 4000 FORMAT( I5 )
C     END OF TOPPPM.
      END
      SUBROUTINE TOPIMG(IB,IP,HD1,HD2,IH,JH,KSUB)
C
      INTEGER       IH,JH,KSUB
      REAL          HD1(IH,JH),HD2(IH,JH)
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
      WRITE(6,*) 'IH,JH,KSUB = ',IH,JH,KSUB
      WRITE(6,*) 
      WRITE(6,*) 'IP(:,47) = ',(ICHAR(IP(K,47)),K=1,3)
      WRITE(6,*) 'IP(:,48) = ',(ICHAR(IP(K,48)),K=1,3)
      WRITE(6,*) 'IP(:,49) = ',(ICHAR(IP(K,49)),K=1,3)
      WRITE(6,*) 'IP(:,50) = ',(ICHAR(IP(K,50)),K=1,3)
      WRITE(6,*) 
      CALL ZHFLSH(6)
      JB = (JH+KSUB-1)/KSUB
      N0 = 0
      N1 = 0
      N2 = 0
      N5 = 0
      D2 =  1.0/50.0
      WRITE(6,*) 'JB,D2 = ',JB,D2
      WRITE(6,*) 'HD.1,1 = ',HD1(1,1),HD2(1,1)
      WRITE(6,*) 
      CALL ZHFLSH(6)
      DO J= 1,JH,KSUB
        DO I= 1,IH,KSUB
          IF     (HD1(I,J).LE.0.0 .AND. HD2(I,J).LE.0.0) THEN
            IB(1,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(1,48)
            IB(2,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(2,48)
            IB(3,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(3,48)
            N0 = N0 + 1
          ELSEIF (HD1(I,J).LE.0.0) THEN
            IB(1,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(1,49)
            IB(2,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(2,49)
            IB(3,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(3,49)
            N1 = N1 + 1
          ELSEIF (HD2(I,J).LE.0.0) THEN
            IB(1,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(1,50)
            IB(2,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(2,50)
            IB(3,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(3,50)
            N2 = N2 + 1
          ELSE
            II = MIN(248, 100 + NINT((HD1(I,J)      )*D2) )
            IB(1,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(1,II)
            IB(2,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(2,II)
            IB(3,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(3,II)
            N5 = N5 + 1
          ENDIF
        ENDDO
      ENDDO
      WRITE(6,*)
      WRITE(6,*) 'TOPIMG STATISTICS:'
      WRITE(6,*) '       LAND 1+2 POINTS = ',N0
      WRITE(6,*) '       LAND 1   POINTS = ',N1
      WRITE(6,*) '       LAND   2 POINTS = ',N2
      WRITE(6,*) '            SEA POINTS = ',N5
      WRITE(6,*)
      CALL ZHFLSH(6)
      RETURN
C     END OF TOPIMG.
      END
