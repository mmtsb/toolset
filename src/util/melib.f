c     MONSSTER energy function for use by external program
c     no simulation capabilities with this code
c     Michael Feig, Brooks group, TSRI, 2002
c
c     original MONSSTER code by:   
c     Andrzej Kolinski, Piotr Rotkiewicz, Jeffrey Skolnick
c     May 1999, May 2000
c     with modifications by Michael Feig, 1999, 2000, 2001
c     Exclusively for academic users.

      SUBROUTINE MSETUP(datdir,plenf,pseq,psec,
     $     paarep,pcompr,pesco,parlo,penone,pehbn,
     $     pars,pburek,penscal,pes3,ptemp)

      IMPLICIT REAL*8 (A-H)
      IMPLICIT INTEGER (I-Z)

      character*256 datdir
      INTEGER plenf
      REAL*8 paarep,pcompr,pesco,parlo,penone,pehbn
      REAL*8 pars,pburek,penscal,pes3,ptemp
      INTEGER pseq(*),psec(*)

      character*5 struct(5)
      character*3 aa(-1:20), NAME
      character*12 text

      PARAMETER(NDIM=500)
      PARAMETER(NBOX=150)

      common/icutt/IRESA,icutp(0:19,0:19),icuta(0:19,0:19),
     *     icutm(0:19,0:19)

      LOGICAL GOODC(646,646)

      DIMENSION vx(646),vy(646),vz(646)
      INTEGER vector(-5:5,-5:5,-5:5)
      DIMENSION eoinp(0:19,0:100)
      DIMENSION ECOVER(0:ndim, 0:24)

      COMMON/GEN/eoinp,ecover,atemp,GOODC,vx,vy,vz,vector

      COMMON/GPAR/aarep,compr,esco,arlo,enone,ehbn,ars,burek,
     $     enscal,es3,arest
      COMMON/BUR/IBUR(0:5),MINUSB(0:5),PLUSB(0:5),IBTARG(0:5)

      DIMENSION map(40)
      DIMENSION apabla(0:19,0:19), apablp(0:19,0:19)
      DIMENSION apablm(0:19,0:19)
      DIMENSION xt(ndim),yt(ndim),zt(ndim), nkb(0:19)

      COMMON /WORK/ xyz(NBOX,NBOX,NBOX)
      COMMON /CHAIN/ ICA(-2:NDIM), x(-2:NDIM), Y(-2:NDIM), Z(-2:NDIM)
      COMMON /SEQE/  SEQ(ndim), SEC(ndim), ENVIR(0:5,0:5,0:5,0:19)
      COMMON /pair/ apba(ndim,ndim),apbp(ndim,ndim)
      COMMON /pairm/apbm(ndim,ndim)
      COMMON /THREE/ ICONF(646,646)

      COMMON /LENGTHS/ LENF2,LENF1,lenf,LENGTH2(646)

      COMMON /SQUARES/ isqr(-646:646)
      COMMON /CENTER/  eone(ndim,0:100), ACRIT, compress, amid, mid

      common/randmc/arand,brand,crand,drand

      COMMON/U/ NCIJA(ndim),NCIJP(ndim),NCIJDA(1600),NCIJDM(1600),
     *     NCIJDP(1600) ,NOA(ndim),NOP(ndim), NCIJM(ndim), NOM(ndim)

      COMMON /SHORTE/ asr3(ndim,4),asr4(ndim,14),asr5(ndim,8),
     $     asr2(ndim,3),ESC,
     *     ibb2(0:199),ibb4(-100:103),ibb5(0:300),ibb3(0:199)

      DIMENSION asr14(0:19,0:19,14),
     *     asr15(0:19,0:19,8),
     *     asr13(0:19,0:19,8), asr12(0:19,0:19,8), abura(0:5)

      DIMENSION asr14H(0:19,0:19,14),
     *     asr15H(0:19,0:19,8),
     *     asr13H(0:19,0:19,8), asr12H(0:19,0:19,8)
      DIMENSION asr14E(0:19,0:19,14),
     *     asr15E(0:19,0:19,8),
     *     asr13E(0:19,0:19,8), asr12E(0:19,0:19,8)
      COMMON /HBN/ EHBOND, NHBN(ndim), prod(646,646)
      LOGICAL MAPR(ndim,ndim)
      COMMON /RES/ EREP,EREST,APHOBI,MRES(ndim),KRES(ndim,50)
      COMMON/RCN/MRESA(ndim),KRESA(ndim,4),iresta(ndim),jresta(ndim)
      DIMENSION ICONT(0:100)
      COMMON /NBN/ NBX(42),NBY(42),NBZ(42)
      LOGICAL PHOB, PHIL
      COMMON /KD/ PHOB(ndim), PHIL(ndim), nkbn(ndim)
      COMMON/CA/cgs(0:20)

      COMMON/CAC/CAR(646*646*3)
      COMMON/CAH/ICAH(646*646*3)
      COMMON/CNH/ICNH(646*646*3)
      COMMON /HB/ HBR(646*646*3)

      COMMON/kdh/   eonekd(0:19)

      COMMON /RENAB/NBPLUS,NCIJ(ndim),NCIJD(1600),INTO(1600)
      COMMON /BURAK/ EBUREK(ndim,0:24)

      DIMENSION STAMAP(ndim,ndim)
      DIMENSION KPBA(-1:19,-1:19)
      COMMON /SHELL1/ SHELL(680,24)
      COMMON /SHELL2/ SPHERE(NDIM,24), SPHEREM(NDIM,24)
      COMMON /SHELL3/ KPB(NDIM,NDIM), NVECTOR(-5:5,-5:5,-5:5)
      COMMON /CCCC/ ECORECT, EQ3(0:19,0:19,0:19)
      COMMON /RNN/ MRESH(ndim),KRESH(ndim,4)
      COMMON /FA/ faa(9:26,0:19,0:19), fra(9:26,0:19,0:19)
      COMMON /FM/ fam(9:26,0:19,0:19), frm(9:26,0:19,0:19)
      COMMON /FP/ fap(9:26,0:19,0:19), frp(9:26,0:19,0:19)
      COMMON /REP/erepp(9:26,0:19,0:19),erepa(9:26,0:19,0:19),
     *     erepm(9:26,0:19,0:19)
c
c     sphere -occupancy points covering (with multiplicity)
c     kpb    -points covering by a specific contact
c     shell  -how to read kpb closest points (ordered numbers
c     in respect to the connecting vector III (600
c     possibilities assumed - check)
c     ibux..- 210, and 200 vectors that define coordination
c     consider a different set (more spherical)
c
      INTEGER IBUX(48),IBUY(48),IBUZ(48)
      DIMENSION INDRI(30),INDII(30)
      DIMENSION BURPERCENT(-1:19,-1:19)
      DIMENSION arla(0:19,0:19),arlm(0:19,0:19),arlp(0:19,0:19)

      data nkb/3,3,4,5,6,5,7,3,7,5,5,7,4,5,6,6,6,7,7,8/
c
c     The set below detects all neighbors up to r=24 ? (included)
c
c
c     IBUX, IBUY, IBZ contain also some trash for long range moves
c
      DATA IBUX /2,2,2,2,-2,-2,-2,-2,
     *     0,0,1,-1,0,0,1,-1,1,-1,0,0,1,-1,0,0,
     *     1,-1,0,0,0,0,   2,-2,0,0,0,0,1,-1,1,-1,1,-1,1,-1,0,0,0,0/

      DATA IBUY /1,-1,0,0,1,-1,0,0,
     *     2,2,2,2,-2,-2,-2,-2,0,0,1,-1,0,0,1,-1,
     *     0,0,1,-1,0,0,   0,0,2,-2,0,0,-1,-1,1,1,0,0,0,0,1,-1,1,-1/

      DATA IBUZ /0,0,1,-1,0,0,1,-1,
     *     1,-1,0,0,1,-1,0,0,2,2,2,2,-2,-2,-2,-2,
     *     0,0,0,0,1,-1,   0,0,0,0,2,-2,0,0,0,0,-1,-1,1,1,-1,-1,1,1/

c
      data struct /' coil','helix',' turn',' beta','    ?'/
      data aa/ 'BCK','GLY','ALA','SER','CYS','VAL','THR','ILE',
     &                     'PRO','MET','ASP','ASN','LEU',
     &                     'LYS','GLU','GLN','ARG',
     &                     'HIS','PHE','TYR','TRP','CYX'/

C     ***********************************************************
C     PREPARATION OF THE CHAIN UNITS AND THEIR CORRELATION
C     (SIDE-CHAIN ONLY REPRESENTATION WITH A GUESS FOR  Ca's)
c
c     Assumed geometry:   8 < |v*v| < 31   and |vx| < 5
c     no co-linear segments, no too accute angles (r>9)
c     646 vectors used - suppresed set
c
C     ***********************************************************
C     PREPARATION OF THE CHAIN UNITS AND THEIR CORRELATION



      NWMAX=0
      DO ix=-5,5
         DO iy=-5,5
            DO iz=-5,5
               vector(ix,iy,iz)=0
               ir=ix*ix+iy*iy+iz*iz
               if(ir.gt.8) then
                  if(ir.lt.31) then
                     NWMAX=NWMAX+1
                     VX(NWMAX)=ix
                     VY(NWMAX)=iy
                     VZ(NWMAX)=iz
                     VECTOR(ix,iy,iz)=NWMAX
                  endif
               endif
            ENDDO
         ENDDO
      ENDDO

      ldatai=0
 212  ldatai=ldatai+1
      if (ichar(datdir(ldatai:ldatai)).ne.0) goto 212
      ldatai=ldatai-1

c
c     3-body **************************************************
c
      OPEN(UNIT=20,FILE=datdir(1:ldatai)//'/SCALE3B',STATUS='OLD')
      REWIND(20)
      do i=0,19
         do j=0,19
            do k=0,19
               EQ3(k,j,i)=0.0
            enddo
         enddo
      enddo

      esum=0.0
      Read(20,*) NUMBER
      do kk=1,NUMBER
         READ(20,1199) text, i,j,k, e
         esum=esum+e
         EQ3(i,j,k)=e
         EQ3(i,k,j)=e
         EQ3(j,i,k)=e
         EQ3(j,k,i)=e
         EQ3(k,i,j)=e
         EQ3(k,j,i)=e
      enddo
 1199 format(A12,I3,2i5,f9.3)
      CLOSE(20)

      esum=esum/dble(NUMBER)

      do i=0,19
         do j=0,19
            do k=0,19
               EQ3(k,j,i)=EQ3(k,j,i)-esum
            enddo
         enddo
      enddo
c
c     3-body **************************************************
c

      lenf=plenf

      OPEN(UNIT=8, FILE=datdir(1:ldatai)//'/PROFILE3',STATUS='OLD')
      OPEN(UNIT=14,FILE=datdir(1:ldatai)//'/S1234',STATUS='OLD')
      OPEN(UNIT=26,FILE=datdir(1:ldatai)//'/QUASI3',STATUS='OLD')
      OPEN(UNIT=27,FILE=datdir(1:ldatai)//'/ECOVERS_24',STATUS='OLD')
      OPEN(UNIT=28,FILE=datdir(1:ldatai)//'/BURIALS',STATUS='OLD')
c
c     ENVIROMENTAL PROFILE UPGRADE
c
      do i=0,19
         do im=0,5
            do ip=0,5
               do ia=0,5
                  envir(im,ip,ia,i)=0.0
               end do
            end do
         end do
      end do

      do i=0,19
         read(8,*)
c     (im is antiparallel,ip is parallel # of contacts)
         do im=0,4
            do ia=0,4
               read(8,*)(envir(ia,im,ip,i),ip=0,4)
               do ip=0,4
                  if(envir(ia,im,ip,i).gt.1.5) envir(ia,im,ip,i)=1.5
               enddo
            end do
            read(8,*)
         end do
         read(8,*)
      end do
c
      sxd=0
      syd=0
      szd=0
c
c
c     BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
c
      read(27,*)
      do i=0,19
         read(27,8765) NAME, (ECOVER(I,J),J=0,24)
         ECOVER(i,1)=ECOVER(i,0)
         do j=0,24
            if(ecover(i,j).gt.1.0)  ecover(i,j)=1.0
         enddo
      ENDDO
      CLOSE(27)

      read(28,*)
      read(28,*)
      read(28,*)
      do i=-1,19
         read(28,8766) NAME, (BURPERCENT(I,J),J=-1,19)
      ENDDO
      CLOSE(28)

c
c     Set_up the discretized scale  (NUMBER OF POINTS IN 24
C     POINT COORDINARTION SPHERE)
c
      do i=-1,19
         do j=-1,19
            KPBA(i,j)=INT(BURPERCENT(I,J)*0.24D0+0.4999D0)
            IF(kpba(i,j).eq.0) kpba(i,j)=1
c
c     this may appear inacurate, however is fine due to structure
c     of the data
c
         ENDDO
      ENDDO
C
C     COMPUTE THE "COVERING SETS" OF POINTS (210)
C
      III=0
      DO ix=-5,5
         DO iy=-5,5
            DO iz=-5,5
               ir=ix*ix+iy*iy+iz*iz
               if(ir.gt.8) then
                  if(ir.lt.31) then
c     may need larger distance for the backbone
                     III=III+1
                     NVECTOR(ix,iy,iz)=iii

C     generate the equivalence set (the closest vectors)
                     do i=1,24
                        ir=(ix-ibux(i))**2+(iy-ibuy(i))**2+
     $                       (iz-ibuz(i))**2
                        INDII(i)=i
                        INDRI(i)=ir
                     enddo
c     ordering according to proximity

                     do i=1,23
                        do j=i+1,24
                           jjr=INDRI(j)
                           iir=INDRI(i)
                           inj=indii(j)
                           ini=indii(i)
                           if(jjr.lt.iir) then
                              INDII(i)=inj
                              INDII(j)=ini
                              INDRI(i)=jjr
                              INDRI(j)=iir
                           endif
                        enddo
                     enddo
c     writting down the results
                     do i=1,24
                        SHELL(iii,i)=INDII(i)
                     enddo

                  endif
               endif
            ENDDO
         ENDDO
      ENDDO
c
c     BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
c

 8765 FORMAT(a3,f6.1,30f5.1)
 8766 FORMAT(a3,21f6.1)

c
      DO I=1,646
         ix=vx(i)
         iy=vy(i)
         iz=vz(i)
         length2(i)=ix*ix+iy*iy+iz*iz
         DO J=1,646
            jx=vx(j)
            jy=vy(j)
            jz=vz(j)
            ICONF(i,j)=(ix+jx)**2+(iy+jy)**2+(iz+jz)**2
            prod(i,j)= ix*jx+iy*jy+iz*jz
         ENDDO
      ENDDO

C     ..............................................................
C
c
c     GENERATE THE NEW BINNING PATTERNS FOR SHORT RANGE
c
c     3-bins of r12

      do i=0,199
         a=sqrt(dble(i)+0.00001D0)*1.45D0

c     >>> 4.6 instead 4.5 for a mobility reason
         if(a .lt. 5.0D0)                          ibb2(i)=1
         if(a .ge. 5.0D0 .and. a .lt. 6.5D0)       ibb2(i)=2
         if(a .ge. 6.5D0 )                         ibb2(i)=3
      end do

c
c     4- BINS OF R13
c
      do i=1,199
         a=sqrt(dble(i)+0.00001D0)*1.45D0

         if(a .lt. 6.5D0)                          ibb3(i)=1
         if(a .ge. 6.5D0 .and. a .lt. 8.0D0)       ibb3(i)=2
         if(a .ge. 8.0D0 .and. a .lt. 9.5D0)       ibb3(i)=3
         if(a .ge. 9.5D0)                          ibb3(i)=4
      end do
c
c     8- bins (chiral of r14)
c
      do i=-100,103
         a=sqrt(dble(iabs(i))+0.00001D0)*1.45D0
         if(i.lt.0) a=-a

         if(a.lt.-12.0D0)                       ibb4(i)=1
         if(a.ge.-12.0D0 .AND. a.lt.-10.0D0)    ibb4(i)=2
         if(a.ge.-10.0D0 .AND. a.lt.-8.0D0)     ibb4(i)=3
         if(a.ge.-8.0D0  .AND. a.lt.-6.0D0)     ibb4(i)=4
         if(a.ge.-6.0D0  .AND. a.lt.-4.0D0)     ibb4(i)=5
         if(a.ge.-4.0D0  .AND. a.lt.-2.0D0)     ibb4(i)=6
         if(a.ge.-2.0D0  .AND. a.lt. 0.0D0)     ibb4(i)=7
         if(a.ge. 0.0D0  .AND. a.lt. 2.0D0)     ibb4(i)=8
         if(a.ge. 2.0D0  .AND. a.lt. 4.0D0)     ibb4(i)=9
         if(a.ge. 4.0D0  .AND. a.lt. 6.0D0)     ibb4(i)=10
         if(a.ge. 6.0D0  .AND. a.lt. 8.0D0)     ibb4(i)=11
         if(a.ge. 8.0D0  .AND. a.lt.10.0D0)     ibb4(i)=12
         if(a.ge.10.0D0  .AND. a.lt.12.0D0)     ibb4(i)=13
         if(a.ge.12.0D0)                        ibb4(i)=14
      end do
c
c     7 BINS OF R15
c
      do i=1,299
         a=sqrt(dble(i)+0.00001D0)*1.45D0

         if(a.lt.5.5D0)                        ibb5(i)=1
         if(a.ge.5.5D0.and.a.lt.7.5D0)         ibb5(i)=2
         if(a.ge.7.5D0.and.a.lt.9.5D0)         ibb5(i)=3
         if(a.ge.9.5D0.and.a.lt.11.5D0)        ibb5(i)=4
         if(a.ge.11.5D0.and.a.lt.13.5D0)       ibb5(i)=5
         if(a.ge.13.5D0.and.a.lt.15.5D0)       ibb5(i)=6
         if(a.ge.15.5D0)                       ibb5(i)=7
      end do

      IICA=0
      DO I=1,646
         ix=vx(i)
         iy=vy(i)
         iz=vz(i)

         DO J=1,646

            If(ICONF(i,j).lt.9.OR.i.eq.j.or.ICONF(i,j).GT.68) then
c
c     limits extracted from the PDB distributions
c
               GOODC(i,j)=.FALSE.
            else
               GOODC(i,j)=.TRUE.

               jx=vx(j)
               jy=vy(j)
               jz=vz(j)
c
c     DESIGN OF "HYDROGEN BOND" GEOMETRY .............   LIST OF
c     EXPECTED POSITION OF H-BONDED NEIGHBORS !!!!!!!!!!!!!!!!!!!!!!
c
               kx=iy*jz-iz*jy
               ky=jx*iz-ix*jz
               kz=ix*jy-iy*jx
c
c     The lengt of H-bond set on ebout 5 Angstroms
c
c
               kkk=kx*kx+ky*ky+kz*kz
c
c     correction for 440 and 330 pairs of vectors (colinear)
c
               if(kkk.eq.0) then
                  GOODC(i,j)=.FALSE.
               else

                  ar=sqrt(dble(kkk))/3.0D0
                  cx=(dble(kx)/ar)
                  cy=(dble(ky)/ar)
                  cz=(dble(kz)/ar)

                  hbr(iica+1)=cx
                  hbr(iica+2)=cy
                  hbr(iica+3)=cz
c
c     GENERATE ALSO THE Ca POSITIONS LISTS  (more exact)
c
                  ax=dble(jx-ix)
                  ay=dble(jy-iy)
                  az=dble(jz-iz)

                  ab=sqrt(ax*ax+ay*ay+az*az)
                  fx=ax/ab
                  fy=ay/ab
                  fz=az/ab

                  car(iica+1)=fx
                  car(iica+2)=fy
                  car(iica+3)=fz

c     GENERATE POSITIONS (slightly expanded) OF Ca Hydrogens
c     Approximate (ca +/- 1.5A H positions on the other side
c     of the SG-SG-SG plane           (this may need a correction)
c

                  px=-ix
                  py=-iy
                  pz=-iz
                  ar=sqrt(dble(length2(i)))

                  bx=dble(px)/ar
                  by=dble(py)/ar
                  bz=dble(pz)/ar

                  ax=fx+cx/3.0D0-bx
                  ay=fy+cy/3.0D0-by
                  az=fz+cz/3.0D0-bz

                  icah(iica+1)=NINT(ax)
                  icah(iica+2)=NINT(ay)
                  icah(iica+3)=NINT(az)

                  ax=fx-cx/2.0D0
                  ay=fy-cy/2.0D0
                  az=fz-cz/2.0D0

                  icnh(iica+1)=NINT(ax)
                  icnh(iica+2)=NINT(ay)
                  icnh(iica+3)=NINT(az)

               endif
            endif
            IICA=IICA+3
         ENDDO
      ENDDO

      do i=0,19
         do j=0,5
            eoinp(i,j)=0.0
         enddo
         do j=0,2
            if(eonekd(i).gt.0.0) eoinp(i,j)=-eonekd(i)
         enddo
         eoinp(i,5)=2.0D0
         eoinp(i,6)=3.0D0
         do j=7,100
            eoinp(i,j)=(dble(j-6)*4.0D0)
         enddo
      enddo
c
c
C     INPUT  INPUT   INPUT   INPUT   INPUT  INPUT  INPUT
C     --------------------------------------------------
C     SET UP OF THE VECTOR REPRESENTATION OF THE CHAIN
C

      LENF1=LENF-1
      LENF2=LENF-2
      AL2=LENF2
      LENF3=LENF-3
      AL3=LENF3
      LENF4=LENF-4
      AL4=LENF4-0.00001D0
      LENF5=LENF-5
      AL5=LENF5
      AL7=lenf-7-0.00001D0
      AL16 =lenf-16
      LENHA=LENF/2
      LENFL=LENF+1

      do 121 i=1,lenf
         sec(i)=psec(i)
         seq(i)=pseq(i)
 121  continue

      do i=1,lenf
         ii=seq(i)
         do j=1,lenf
            stamap(j,i)=0
            jj=seq(j)
            KPB(i,j)=KPBA(ii,jj)
         enddo
c     backbone contribution
         KPB(i,1)=KPBA(ii,-1)
      enddo

c
c     READ THE SHORT RANGE INTERACTIONS (OVERALL)
c
c
      read(14,*)
      do i=0,19
         do j=0,19
            read(14,*)
            read(14,*) (asr12(i,j,k),k=1,3)
         enddo
      enddo

      read(14,*)
      do i=0,19
         do j=0,19
            read(14,*)
            read(14,*) (asr13(i,j,k),k=1,4)
         enddo
      enddo

      read(14,*)
      do i=0,19
         do j=0,19
            read(14,*)
            read(14,*) (asr14(i,j,k),k=1,14)
         enddo
      enddo

      read(14,*)
      do i=0,19
         do j=0,19
            read(14,*)
            read(14,*) (asr15(i,j,k),k=1,7)
         enddo
      enddo

      close(14)

      OPEN(UNIT=14,FILE=datdir(1:ldatai)//'/S1234H',STATUS='OLD')

c
c     READ THE SHORT RANGE INTERACTIONS (HELICAL)
c
c
      read(14,*)
      do i=0,19
         do j=0,19
            read(14,*)
            read(14,*) (asr12H(i,j,k),k=1,3)
         enddo
      enddo

      read(14,*)
      do i=0,19
         do j=0,19
            read(14,*)
            read(14,*) (asr13H(i,j,k),k=1,4)
         enddo
      enddo

      read(14,*)
      do i=0,19
         do j=0,19
            read(14,*)
            read(14,*) (asr14H(i,j,k),k=1,14)
         enddo
      enddo

      read(14,*)
      do i=0,19
         do j=0,19
            read(14,*)
            read(14,*) (asr15H(i,j,k),k=1,7)
         enddo
      enddo

      close(14)

      OPEN(UNIT=14,FILE=datdir(1:ldatai)//'/S1234E',STATUS='OLD')

c
c     READ THE SHORT RANGE INTERACTIONS (EXTENDED)
c
c
      read(14,*)
      do i=0,19
         do j=0,19
            read(14,*)
            read(14,*) (asr12E(i,j,k),k=1,3)
         enddo
      enddo

      read(14,*)
      do i=0,19
         do j=0,19
            read(14,*)
            read(14,*) (asr13E(i,j,k),k=1,4)
         enddo
      enddo

      read(14,*)
      do i=0,19
         do j=0,19
            read(14,*)
            read(14,*) (asr14E(i,j,k),k=1,14)
         enddo
      enddo

      read(14,*)
      do i=0,19
         do j=0,19
            read(14,*)
            read(14,*) (asr15E(i,j,k),k=1,7)
         enddo
      enddo

      close(14)

C
C     ***********************INPUT FILE******************************
C

      AAREP=paarep
      compr=pcompr
      ESCO=pesco
      ARLO=parlo
      ENONE=penone
      EHBN=pehbn
      ARS=pars
      BUREK=pburek
      ENSCAL=penscal
      ES3=pes3
      ATEMP=ptemp
c
c     Set a generalized shell model
c
      abura(0)= 2.42731D-02
      abura(1)= 0.1650727D0
      abura(2)= 0.3647117D0
      abura(3)= 0.342935D0
      abura(4)= 9.61324D-02
      abura(5)= 5.56900D-03

      DO i=0,5
         IBTARG(i)=INT(al2*abura(i)+0.5D0)
      enddo
C
      do i=1,lenf
         MRES(i)=0
         MRESA(i)=0
         MRESH(i)=0
         do j=1,lenf
            mapr(j,i)=.FALSE.
         enddo
      enddo

c
c     READING TERTARY INTERACTION SCALES WITH 2-body -
c     ORIENTATION DEPENDENT .............!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
      read(26,*)
      read(26,*)
      do i=0,19
         read(26,725) NAME, (apablp(i,j),j=0,19)
      enddo
      read(26,*)
      read(26,*)
      do i=0,19
         read(26,725) NAME, (apablm(i,j),j=0,19)
      enddo
      read(26,*)
      read(26,*)
      do i=0,19
         read(26,725) NAME, (apabla(i,j),j=0,19)
      enddo
      read(26,*)
      read(26,*)
      do i=0,19
         read(26,725) NAME, (arlp(i,j),j=0,19)
      enddo
      do i=0,19
         do j=0,19
            icutp(i,j)=INT(((arlp(i,j)+0.1D0)/1.45D0)**2+2.0D0)
            if(icutp(i,j).lt.17) icutp(i,j)=17
            if(icutp(i,j).gt.23) icutp(i,j)=23
         enddo
      enddo

      read(26,*)
      read(26,*)
      do i=0,19
         read(26,725) NAME, (arlm(i,j),j=0,19)
      enddo
      do i=0,19
         do j=0,19
            icutm(i,j)=INT(((arlm(i,j)+0.1D0)/1.45D0)**2+2.0D0)
            if(icutm(i,j).lt.17) icutm(i,j)=17
            if(icutm(i,j).gt.23) icutm(i,j)=23
         enddo
      enddo

      read(26,*)
      read(26,*)
      do i=0,19
         read(26,725) NAME, (arla(i,j),j=0,19)
      enddo
      do i=0,19
         do j=0,19
            icuta(i,j)=INT(((arla(i,j)+0.1D0)/1.45D0)**2+2.0D0)
            if(icuta(i,j).lt.17) icuta(i,j)=17
            if(icuta(i,j).gt.23) icuta(i,j)=23
         enddo
      enddo

 725  format(a3,1x,20f5.1)

      do i=0,19
         do j=0,19
            do ir=9,25
               erepa(ir,i,j)=0.0
               if(ir.lt.icuta(i,j)) then
                  faa(ir,i,j)=1.0
                  fra(ir,i,j)=1.0
               else
                  if(ir.lt.icuta(i,j)+3) faa(ir,i,j)=0.25D0
                  if(ir.eq.icuta(i,j)) faa(ir,i,j)=0.5D0
                  fra(ir,i,j)=0.0
               endif
            enddo
         enddo
      enddo

      do i=0,19
         do j=0,19
            do ir=9,25
               erepm(ir,i,j)=0.0
               if(ir.lt.icutm(i,j)) then
                  fam(ir,i,j)=1.0
                  frm(ir,i,j)=1.0
               else
                  if(ir.lt.icutm(i,j)+3) fam(ir,i,j)=0.25D0
                  if(ir.eq.icutm(i,j)) fam(ir,i,j)=0.5D0
                  frm(ir,i,j)=0.0
               endif
            enddo
         enddo
      enddo

      do i=0,19
         do j=0,19
            do ir=9,25
               erepp(ir,i,j)=0.0
               if(ir.lt.icutp(i,j)) then
                  fap(ir,i,j)=1.0
                  frp(ir,i,j)=1.0
               else
                  if(ir.lt.icutp(i,j)+3) fap(ir,i,j)=0.25D0
                  if(ir.eq.icutp(i,j)) fap(ir,i,j)=0.5D0
                  frp(ir,i,j)=0.0
               endif
            enddo
         enddo
      enddo

      do i=0,19
         do j=0,19
            do ir=9,13
               asr=sqrt(dble(ir))+1.125D0
               if(sqrt(dble(icuta(i,j))).gt.asr) erepa(ir,i,j)=1.0
               if(sqrt(dble(icutm(i,j))).gt.asr) erepm(ir,i,j)=1.0
               if(sqrt(dble(icutp(i,j))).gt.asr) erepp(ir,i,j)=1.0
            enddo
         enddo
      enddo

      Acrit= 2.2D0*exp(0.38D0*log(AL2))/1.45D0
c

      MAXIM=NBOX
      MID=MAXIM/2
      amid=dble(mid)

      close(26)

      do i=1,lenf
         ii=seq(i)
         nkbn(i)=nkb(ii)
         if(eonekd(ii).gt.0.0) then
            PHOB(i)=.TRUE.
            PHIL(i)=.FALSE.
         else
            PHOB(i)=.FALSE.
            PHIL(i)=.TRUE.
         endif
         do 128 j=1,lenf

            jj=seq(j)
            apba(i,j)=apabla(ii,jj)-0.25D0
            apbm(i,j)=apablm(ii,jj)-0.25D0
            apbp(i,j)=apablp(ii,jj)-0.25D0

            if(iabs(i-j).gt.4) then
               apba(i,j)=apba(i,j)-0.25D0
               apbm(i,j)=apbm(i,j)-0.25D0
               apbp(i,j)=apbp(i,j)-0.25D0

               if(i.gt.2.and.j.gt.2) then
                  if(eonekd(ii).gt.0. .AND. eonekd(jj).gt.0.) then
                     apba(i,j)=apba(i,j)-(eonekd(ii)+eonekd(jj))/32.0D0
                     apbm(i,j)=apbm(i,j)-(eonekd(ii)+eonekd(jj))/32.0D0
                     apbp(i,j)=apbp(i,j)-(eonekd(ii)+eonekd(jj))/32.0D0
                  endif
               endif
            endif

            if(iabs(i-j).lt.4) apba(i,j)=0.0
            if(iabs(i-j).lt.3) apbm(i,j)=0.0

            if(MAPR(i,j)) then
               apba(i,j)=min(-1.0D0,apba(i,j)-0.25D0)
               apbm(i,j)=min(-1.0D0,apbm(i,j)-0.25D0)
               apbp(i,j)=min(-1.0D0,apbp(i,j)-0.25D0)
            endif

            if(iabs(i-j).lt.2) then
               apba(i,j)=0.0
               apbm(i,j)=0.0
               apbp(i,j)=0.0
            endif

 128     continue
      enddo

c      atemp=1.0D0

c
c     3-body **************************************************
c
      do i=0,19
         do j=0,19
            do k=0,19
               EQ3(k,j,i)=ES3*EQ3(k,j,i)/ATEMP
            enddo
         enddo
      enddo
c
c     3-body **************************************************
c

c
c
c     ENVIROMENTAL PROFILE UPGRADE
c
      do i=0,19
         do im=0,4
            do ip=0,4
               do ia=0,4
                  envir(ia,im,ip,i)=ENSCAL*envir(ia,im,ip,i)/ATEMP
               end do
            end do
         end do
      end do

      do i=0,19
         do j=0,19
            do ir=9,13
               erepa(ir,i,j)=AAREP*erepa(ir,i,j)/atemp
               erepm(ir,i,j)=AAREP*erepm(ir,i,j)/atemp
               erepp(ir,i,j)=AAREP*erepp(ir,i,j)/atemp
            enddo
         enddo
      enddo

      do i=1,lenf
         do j=2,lenf1
            apba(j,i)=arlo*apba(j,i)/atemp
            apbm(j,i)=arlo*apbm(j,i)/atemp
            apbp(j,i)=arlo*apbp(j,i)/atemp
         enddo
      enddo

c
c
      do i=1,lenf-1
         do k=1,3
            asr2(i,k)=0.5D0*ars*asr12(seq(i),seq(i+1),k)/ATEMP
         enddo
         if(sec(i).eq.2.and.sec(i+1).eq.2) then
            do k=1,3
               asr2(i,k)=0.5D0*ars*asr12H(seq(i),seq(i+1),k)/ATEMP
            enddo
         endif
         if(sec(i).eq.4.and.sec(i+1).eq.4) then
            do k=1,3
               asr2(i,k)=0.5D0*ars*asr12E(seq(i),seq(i+1),k)/ATEMP
            enddo
         endif
      enddo

      do i=1,lenf-2
         do k=1,4
            asr3(i,k)=0.5D0*ars*asr13(seq(i),seq(i+2),k)/ATEMP
         enddo
         if(sec(i).eq.2.and.sec(i+1).eq.2.and.sec(i+2).eq.2) then
            do k=1,4
               asr3(i,k)=0.5D0*ars*asr13H(seq(i),seq(i+2),k)/ATEMP
            enddo
         endif
        if(sec(i).eq.4.and.sec(i+1).eq.4.and.sec(i+2).eq.4) then
           do k=1,4
              asr3(i,k)=0.5D0*ars*asr13E(seq(i),seq(i+2),k)/ATEMP
           enddo
        endif
      enddo

      do i=1,lenf-3
         do k=1,14
            asr4(i,k)=ars*(asr14(seq(i),seq(i+3),k))/ATEMP
         enddo
         if(sec(i+1).eq.2.and.sec(i+2).eq.2) then
            if(sec(i).eq.2.and.sec(i+3).eq.2) then
               do k=1,14
                  asr4(i,k)=ars*(asr14H(seq(i),seq(i+3),k))/ATEMP
               enddo
            endif
         endif
         if(sec(i+1).eq.4.and.sec(i+2).eq.4) then
            do k=1,14
               asr4(i,k)=ars*(asr14E(seq(i),seq(i+3),k))/ATEMP
            enddo
         endif
      enddo

      do i=1,lenf-4
         do k=1,7
            asr5(i,k)=ars*(asr15(seq(i),seq(i+4),k))/ATEMP
         enddo
         if(sec(i+1).eq.2.and.sec(i+2).eq.2.and.sec(i+3).eq.2) then
            if(sec(i).eq.2.and.sec(i+4).eq.2) then
               do k=1,7
                  asr5(i,k)=ars*(asr15H(seq(i),seq(i+4),k))/ATEMP
               enddo
            endif
         endif
         if(sec(i+1).eq.4.and.sec(i+2).eq.4.and.sec(i+3).eq.4) then
            do k=1,7
               asr5(i,k)=ars*(asr15E(seq(i),seq(i+4),k))/ATEMP
            enddo
         endif
      enddo

      DO i=1,NBOX
         DO j=1,NBOX
            DO k=1,NBOX
               xyz(k,j,i)=0
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END

c *************************************************************
c *************************************************************
c *************************************************************

      SUBROUTINE MENER(tcx,tcy,tcz,
     $     energ,esh,epr,ebur,ecor,ecent)

      IMPLICIT REAL*8 (A-H)
      IMPLICIT INTEGER (I-Z)

      INTEGER TCX(*),TCY(*),TCZ(*)
      REAL*8 energ,ecent,esh,epr,ebur,ecor

      PARAMETER(NDIM=500)
      PARAMETER(NBOX=150)

      LOGICAL GOODC(646,646)
      INTEGER vx(646),vy(646),vz(646)
      INTEGER vector(-5:5,-5:5,-5:5)
      REAL*8 eoinp(0:19,0:100)
      DIMENSION ECOVER(0:ndim, 0:24)
      COMMON/GEN/eoinp,ecover,atemp,GOODC,vx,vy,vz,vector
      COMMON/GPAR/aarep,compr,esco,arlo,enone,ehbn,ars,burek,
     $     enscal,es3,arest
      COMMON/BUR/IBUR(0:5),MINUSB(0:5),PLUSB(0:5),IBTARG(0:5)

      COMMON /WORK/ xyz(NBOX,NBOX,NBOX)
      COMMON /CHAIN/ ICA(-2:NDIM), x(-2:NDIM), Y(-2:NDIM), Z(-2:NDIM)
      COMMON /SEQE/  SEQ(ndim), SEC(ndim), ENVIR(0:5,0:5,0:5,0:19)
      COMMON /pair/ apba(ndim,ndim),apbp(ndim,ndim)
      COMMON /pairm/apbm(ndim,ndim)
      COMMON /THREE/ ICONF(646,646)

      COMMON /LENGTHS/ LENF2,LENF1,lenf,LENGTH2(646)

      COMMON /SQUARES/ isqr(-646:646)
      COMMON /CENTER/  eone(ndim,0:100), ACRIT, compress, amid, mid

      common/randmc/arand,brand,crand,drand

      COMMON/U/ NCIJA(ndim),NCIJP(ndim),NCIJDA(1600),NCIJDM(1600),
     *     NCIJDP(1600) ,NOA(ndim),NOP(ndim), NCIJM(ndim), NOM(ndim)

      COMMON /SHORTE/ asr3(ndim,4),asr4(ndim,14),asr5(ndim,8),
     $     asr2(ndim,3),ESC,
     *     ibb2(0:199),ibb4(-100:103),ibb5(0:300),ibb3(0:199)

      COMMON /HBN/ EHBOND, NHBN(ndim), prod(646,646)
      COMMON /RES/ EREP,EREST,APHOBI,MRES(ndim),KRES(ndim,50)
      COMMON/RCN/MRESA(ndim),KRESA(ndim,4),iresta(ndim),jresta(ndim)
      COMMON /NBN/ NBX(42),NBY(42),NBZ(42)
      LOGICAL PHOB, PHIL
      COMMON /KD/ PHOB(ndim), PHIL(ndim), nkbn(ndim)
      COMMON/CA/cgs(0:20)

      COMMON/CAC/CAR(646*646*3)
      COMMON/CAH/ICAH(646*646*3)
      COMMON/CNH/ICNH(646*646*3)
      COMMON /HB/ HBR(646*646*3)

      COMMON/kdh/   eonekd(0:19)

      COMMON /RENAB/NBPLUS,NCIJ(ndim),NCIJD(1600),INTO(1600)
      COMMON /BURAK/ EBUREK(ndim,0:24)

      COMMON /SHELL1/ SHELL(680,24)
      COMMON /SHELL2/ SPHERE(NDIM,24), SPHEREM(NDIM,24)
      COMMON /SHELL3/ KPB(NDIM,NDIM), NVECTOR(-5:5,-5:5,-5:5)
      COMMON /CCCC/ ECORECT, EQ3(0:19,0:19,0:19)
      COMMON /RNN/ MRESH(ndim),KRESH(ndim,4)
      COMMON /FA/ faa(9:26,0:19,0:19), fra(9:26,0:19,0:19)
      COMMON /FM/ fam(9:26,0:19,0:19), frm(9:26,0:19,0:19)
      COMMON /FP/ fap(9:26,0:19,0:19), frp(9:26,0:19,0:19)
      COMMON /REP/erepp(9:26,0:19,0:19),erepa(9:26,0:19,0:19),
     *     erepm(9:26,0:19,0:19)

c     sphere -occupancy points covering (with multiplicity)
c     kpb    -points covering by a specific contact
c     shell  -how to read kpb closest points (ordered numbers
c     in respect to the connecting vector III (600
c     possibilities assumed - check)
c     ibux..- 210, and 200 vectors that define coordination
c     consider a different set (more spherical)
c
c     The set below detects all neighbors up to r=24 ? (included)
c

      DO I=1,LENF
         X(I)=TCX(I)
         Y(I)=TCY(I)
         Z(I)=TCZ(I)
      ENDDO

      xshrf=0
      yshrf=0
      zshrf=0

      DO I=1,LENF1
         J=I+1
         WX=X(J)-X(I)
         WY=Y(J)-Y(I)
         WZ=Z(J)-Z(I)
         ICA(I)=VECTOR(WX,WY,WZ)
      ENDDO

      ICA(lenf+1)=0
      ICA(0)=0
      ICA(-1)=0
      ICA(-2)=0
c

      DO J=2,LENF-1
         I=J-1
         II=ICA(I)
         JJ=ICA(J)
         IF(.NOT.GOODC(II,JJ)) THEN
             WRITE(*,8011)I,J,vx(ii),vy(ii),vz(ii),
     $            vx(jj),vy(jj),vz(jj)
 8011        FORMAT(5X,
     $            'WARNING -WRONG INPUT CHAIN - VECTORS ',8I4)
         ENDIF
      ENDDO
      
      do i=1,lenf
         ii=seq(i)
         do j=0,100
            if(i.gt.1.and.i.lt.lenf) then
               eone(i,j)=enone*eoinp(ii,j)/atemp
            else
               eone(i,j)=0.0
            endif
         enddo
      enddo

c
c     BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
c

      do i=1,lenf
         ii=seq(i)
         do j=0, 24
c     j - number of contacts occupied points
            if(j.gt.1.and.j.lt.24) then
               eburek(i,j)=BUREK*(2.0D0*ECOVER(ii,j)+ECOVER(ii,j-1)
     *              +ECOVER(ii,j+1))/ATEMP/4.0D0
            else
               eburek(i,j)=BUREK*ECOVER(ii,j)/ATEMP
            endif
         enddo
      enddo
c
      ESC= ESCO/ATEMP
      EHBOND=ehbn/ATEMP
      EREST=AREST/ATEMP
      EREP=AAREP/ATEMP
      compress=compr/atemp

      do i=0,5
         ibur(i)=0
      enddo

      esh=ESHORT(1,lenf)
      
      SX=0
      SY=0
      SZ=0
      DO I=1,LENF
         SX=SX+X(I)
         SY=SY+Y(I)
         SZ=SZ+Z(I)
      ENDDO
      
c      endif
      
      SX=SX/LENF
      SY=SY/LENF
      SZ=SZ/LENF
      XSHIFT=MID-SX
      YSHIFT=MID-SY
      ZSHIFT=MID-SZ
      
      xshrf=xshrf+xshift
      yshrf=yshrf+yshift
      zshrf=zshrf+zshift
      
      DO I=1,LENF
         X(I)=X(I)+XSHIFT
         Y(I)=Y(I)+YSHIFT
         Z(I)=Z(I)+ZSHIFT
         
         do j=1,24
            sphere(i,j)=0
            spherem(i,j)=0
         enddo
         NHBN(i)=0
         NCIJ(i)=0
         NCIJA(i)=0
         NCIJP(i)=0
c     mod meikel
         NCIJM(i)=0
         NOA(i)=11
         NOP(i)=11
      ENDDO
      NBPLUS=0
         
      epr=EPAIRN(2,lenf1)
c
c     Correction for the initial 0-contact residues
c
      ebur=0.0D0
      do i=1,lenf
         if(ncij(i).eq.0)
     $        Ebur=Ebur+EBUREK(i,0)+ENVIR(0,0,0,seq(i))
      enddo
c     
      ecor=ecorect
      ECENT=ECENTN(1,lenf)
      ENERG=esh+epr+ebur+ecor+ecent

      return 
      END
c

c
c
c     **************************************************************

      BLOCK DATA foribm

      IMPLICIT REAL*8 (A-H)
      IMPLICIT INTEGER (I-Z)

      COMMON /NBN/ NBX(42),NBY(42),NBZ(42)
      COMMON/CA/cgs(0:20)
      COMMON/kdh/   eonekd(0:19)

      DATA NBX /3,-3,3,-3,0,0,0,0,3,-3,3,-3,
     *     3,3,3,3,-3,-3,-3,-3,2,-2,-2,2,2,-2,-2,2,2,-2,-2,2,
     *     2,-2,-2,2,                             0,0,0,0,4,-4/


      DATA NBY /3,-3,-3,3,3,-3,3,-3,0,0,0,0,
     *     2,-2,-2,2,2,-2,-2,2,3,3,3,3,-3,-3,-3,-3,2,-2,2,-2,
     *     2,-2,2,-2,                        0,0,4,-4,0,0/


      DATA NBZ /0,0,0,0,3,-3,-3,3,3,-3,-3,3,
     *     2,-2,2,-2,2,-2,2,-2,2,-2,2,-2,2,-2,2,-2,3,3,3,3,
     *     -3,-3,-3,-3,                        4,-4,0,0,0,0/


      data cgs/     0.0D0, 0.5D0, 0.75D0,  0.75D0, 1.0D0,
     &     1.0D0, 1.25D0,
     &     0.5D0, 1.5D0, 1.25D0, 1.25D0, 1.25D0,
     &     2.0D0, 1.5D0, 1.5D0, 2.0D0,
     &     1.5D0, 1.5D0, 1.75D0, 2.0D0, 0.75D0/


      data eonekd /-0.4D0, 1.8D0, -0.8D0, 2.5D0,
     &     4.2D0, -0.7D0, 4.5D0,
     &     -1.6D0, 1.9D0,  -3.5D0, -3.5D0, 3.8D0,
     &     -3.9D0, -3.5D0, -3.5D0, -4.5D0,
     &     -3.2D0, 2.8D0, -1.3D0, -0.9D0/

      end


C     **************************************************************

      SUBROUTINE SET(i, j, k, m)
      IMPLICIT REAL*8 (A-H)
      IMPLICIT INTEGER (I-Z)

      PARAMETER(NBOX=150)
      PARAMETER(NDIM=500)
      COMMON /CHAIN/ ICA(-2:NDIM), x(-2:NDIM), Y(-2:NDIM), Z(-2:NDIM)
      COMMON/WORK/ xyz(NBOX,NBOX,NBOX)

      XYZ(i,j,k)=m

      XYZ(i+1,j,k)=m
      XYZ(i-1,j,k)=m
      XYZ(i,j+1,k)=m
      XYZ(i,j-1,k)=m
      XYZ(i,j,k+1)=m
      XYZ(i,j,k-1)=m

      XYZ(i+1,j+1,k)=m
      XYZ(i+1,j-1,k)=m
      XYZ(i-1,j+1,k)=m
      XYZ(i-1,j-1,k)=m

      XYZ(i,j+1,k+1)=m
      XYZ(i,j+1,k-1)=m
      XYZ(i,j-1,k+1)=m
      XYZ(i,j-1,k-1)=m

      XYZ(i+1,j,k+1)=m
      XYZ(i+1,j,k-1)=m
      XYZ(i-1,j,k+1)=m
      XYZ(i-1,j,k-1)=m


      RETURN
      END

C     **************************************************************

      SUBROUTINE REM(i, j, k)
C
      IMPLICIT REAL*8 (A-H)
      IMPLICIT INTEGER (I-Z)

      PARAMETER(NBOX=150)
      PARAMETER(NDIM=500)
      COMMON /CHAIN/ ICA(-2:NDIM), x(-2:NDIM), Y(-2:NDIM), Z(-2:NDIM)
      COMMON/WORK/ xyz(NBOX,NBOX,NBOX)

      XYZ(i,j,k)=0

      XYZ(i+1,j,k)=0
      XYZ(i-1,j,k)=0
      XYZ(i,j+1,k)=0
      XYZ(i,j-1,k)=0
      XYZ(i,j,k+1)=0
      XYZ(i,j,k-1)=0

      XYZ(i+1,j+1,k)=0
      XYZ(i+1,j-1,k)=0
      XYZ(i-1,j+1,k)=0
      XYZ(i-1,j-1,k)=0

      XYZ(i,j+1,k+1)=0
      XYZ(i,j+1,k-1)=0
      XYZ(i,j-1,k+1)=0
      XYZ(i,j-1,k-1)=0

      XYZ(i+1,j,k+1)=0
      XYZ(i+1,j,k-1)=0
      XYZ(i-1,j,k+1)=0
      XYZ(i-1,j,k-1)=0

      RETURN
      END

C     **************************************************************

      FUNCTION ESHORT(iiii,jjjj)
      IMPLICIT REAL*8 (A-H)
      IMPLICIT INTEGER (I-Z)

      PARAMETER(NDIM=500)

      COMMON /SEQE/  SEQ(ndim), SEC(ndim), ENVIR(0:5,0:5,0:5,0:19)
      COMMON /THREE/ ICONF(646,646)
      COMMON/CA/cgs(0:20)

      COMMON/CAC/CAR(646*646*3)
      COMMON /HB/ HBR(646*646*3)

      COMMON /SHORTE/ asr3(ndim,4),asr4(ndim,14),asr5(ndim,8),
     $     asr2(ndim,3),ESC,
     *     ibb2(0:199),ibb4(-100:103),ibb5(0:300),ibb3(0:199)

      COMMON vx(646),vy(646),vz(646)
      COMMON /CHAIN/ ICA(-2:NDIM), x(-2:NDIM), Y(-2:NDIM), Z(-2:NDIM)
      COMMON /LENGTHS/ LENF2,LENF1,LENF,LENGTH2(646)

      COMMON /HBN/ EHBOND, NHBN(ndim), prod(646,646)
      COMMON /CENTER/  eone(ndim,0:100), ACRIT, compress, amid, mid
      COMMON /FR/ FRG(ndim,ndim)

      LOGICAL NOGLY

      ESHORT=0.0D0
      ESC1=ESC*2.0D0
      ESC2=ESC/2.0D0
      ESC3=ESC*1.5D0


c     TEST OF GENERIC STIFFNESS  (GENERAL approach)
c     R15 criterion assuming 1.45A spacing
c
c     SHORT RANGE SEQUENCE SPECIFFIC
c
      i1=iiii-1
      if(i1.lt.1) i1=1
      i2=jjjj
      if(i2.gt.lenf1) i2=lenf1

c     + add the backbone effect
c
c     Build a pseudo-backbone and evaluate
c
      IF(i1.gt.2) then
         fj=cgs(seq(i1-1))
         iica=ica(i1-2)*1938+ica(i1-1)*3-1941
         aax=(dble(x(i1-1))+CAR(iica+1)*fj)
         aay=(dble(y(i1-1))+CAR(iica+2)*fj)
         aaz=(dble(z(i1-1))+CAR(iica+3)*fj)
      endif

      if (i1.gt.1) THEN
         fj=cgs(seq(i1))
         iica=ica(i1-1)*1938+ica(i1)*3-1941
         asx=(dble(x(i1))+CAR(iica+1)*fj)
         asy=(dble(y(i1))+CAR(iica+2)*fj)
         asz=(dble(z(i1))+CAR(iica+3)*fj)
      ENDIF

      do i=i1,i2
         icai=ica(i)
         ican=ica(i-1)
         icaj=ica(i+1)

         j=i+1
         fj=cgs(seq(j))
         iica=ica(i)*1938+ica(i+1)*3-1941
         bsx=(dble(x(j))+CAR(iica+1)*fj)
         bsy=(dble(y(j))+CAR(iica+2)*fj)
         bsz=(dble(z(j))+CAR(iica+3)*fj)

         if (i.gt.1) THEN
            ar=(asx-bsx)**2+(asy-bsy)**2+(asz-bsz)**2
            ESHORT=ESHORT+ESC2*(sqrt(ar)-2.6D0)**2
            if(i.gt.2) then
               ar=sqrt((aax-bsx)**2+(aay-bsy)**2+(aaz-bsz)**2)
               if(ar.gt.5.5D0) ar=5.5D0+(ar-5.5D0)/3.0D0
               ESHORT=ESHORT+ESC2*abs((ar-3.6D0)*(4.5D0-ar))
            endif
         ENDIF

         aax=asx
         aay=asy
         aaz=asz
         asx=bsx
         asy=bsy
         asz=bsz

         ib2=ibb2(length2(icai))
         ESHORT=ESHORT+asr2(i,ib2)

         if (i.lt.lenf1 .and. i.gt.1) then
            WX1=VX(ican)
            WY1=VY(ican)
            WZ1=VZ(ican)
            WX2=VX(icai)
            WY2=VY(icai)
            WZ2=VZ(icai)
            WX3=VX(icaj)
            WY3=VY(icaj)
            WZ3=VZ(icaj)
            ir=(wx1+wx2+wx3)**2+(wy1+wy2+wy3)**2+(wz1+wz2+wz3)**2
            if(ir.gt.100) ir=100
            px=wy1*wz2-wy2*wz1
            py=wz1*wx2-wz2*wx1
            pz=wx1*wy2-wx2*wy1
            ihand=px*wx3+py*wy3+pz*wz3
            if(ihand.lt.0) ir=-ir

c
c     combined profile enters into the short range interactions
c
            if(i.gt.1) then
               ib3=ibb3(iconf(icai,icaj))
               ESHORT=ESHORT+asr3(i,ib3)

               if (i.gt.2) then
                  ib4=ibb4(ir)
                  ESHORT=ESHORT+asr4(i-1,ib4)
               endif
            endif
         endif
      enddo

      if(i2.lt.lenf1) then
         i=i2+1
         j=i+1
         fj=cgs(seq(j))
         iica=ica(i)*1938+ica(j)*3-1941
         bsx=(dble(x(j))+CAR(iica+1)*fj)
         bsy=(dble(y(j))+CAR(iica+2)*fj)
         bsz=(dble(z(j))+CAR(iica+3)*fj)
         ar=sqrt((aax-bsx)**2+(aay-bsy)**2+(aaz-bsz)**2)
         if(ar.gt.5.5D0) ar=5.5D0+(ar-5.5D0)/3.0D0
         ESHORT=ESHORT+ESC2*abs((ar-3.6D0)*(4.5D0-ar))
      endif

      i1=iiii-4
      if(i1.lt.2) i1=2
      i2=jjjj+1

      if(i2.gt.lenf1-4) i2=lenf1-4
      do i=i1,i2
         j=i+4
         icam4=ica(i-1)
         icam3=ica(i)
         icam2=ica(i+1)
         icam1=ica(i+2)
         icam=ica(i+3)
         icaj=ica(j)
         ix=x(j)-x(i)
         iy=y(j)-y(i)
         iz=z(j)-z(i)
     
         NOGLY=.FALSE.
         if(seq(i+1)*seq(i+2)*seq(i+3)*seq(j).gt.0) NOGLY=.TRUE.

         ax=dble(x(i+1)+x(i+2)+x(i+3)+x(j))/4.0D0
         ay=dble(y(i+1)+y(i+2)+y(i+3)+y(j))/4.0D0
         az=dble(z(i+1)+z(i+2)+z(i+3)+z(j))/4.0D0
         bx=ax-amid
         by=ay-amid
         bz=az-amid
         aract2=bx*bx+by*by+bz*bz
         ff=2.0D0-(aract2/acrit**2)
         if(ff.lt.0.5D0) ff=0.5D0
         if(ff.gt.1.0D0) ff=1.0D0

         jjj=(vx(icam3)-vx(icaj))**2+(vy(icam3)-vy(icaj))**2
     *        +(vz(icam3)-vz(icaj))**2-45

         if(jjj.lt.0) THEN
            if(jjj.lt.-30) then
               if(NOGLY) then
                  jjj=-45
               else
                  jjj=-30
               endif
            endif
            ESHORT=ESHORT+ff*ESC*dble(jjj)/90.0D0
         ENDIF

         iica1=ica(i)*1938+ica(i+1)*3-1941
         iica2=ica(i+2)*1938+ica(i+3)*3-1941

         ax=HBR(iica1+1)*HBR(iica2+1)
         ay=HBR(iica1+2)*HBR(iica2+2)
         az=HBR(iica1+3)*HBR(iica2+3)

         a=ax+ay+az

         if(a.gt.0.0D0) then
            if(a.gt.6.0D0) a=9.0D0
            ESHORT=ESHORT-ff*ESC*a/18.0D0
         endif

         iii=ix**2+iy**2+iz**2

         if(iii.lt.33) then
            if(sec(i+1).ne.4.AND.sec(i+2).ne.4.AND.sec(i+3).ne.4) then
               if(prod(icam3,icam) .gt.0) then
                  WX1=VX(icam3)
                  WY1=VY(icam3)
                  WZ1=VZ(icam3)
                  WX2=VX(icam2)
                  WY2=VY(icam2)
                  WZ2=VZ(icam2)
                  WX3=VX(icam1)
                  WY3=VY(icam1)
                  WZ3=VZ(icam1)
                  px=wy1*wz2-wy2*wz1
                  py=wz1*wx2-wz2*wx1
                  pz=wx1*wy2-wx2*wy1
                  ihand=px*wx3+py*wy3+pz*wz3
                  if(ihand.gt.0) then
                     ij=px*vx(icam)+py*vy(icam)+pz*vz(icam)
                     if(ij.gt.0) ESHORT=ESHORT-ESC-ff*ESC2
                  endif
               endif
            endif

         else

            if(iii.gt.60.and.iii.lt.125)         then
               if(sec(i+1).ne.2 .AND. sec(i+2).ne.2
     $              .AND. sec(i+3).ne.2) then
                  if(a.gt.6.0D0) ESHORT=ESHORT-ESC-ff*ESC2
               endif
            endif
         endif

         if(iii.gt.299) iii=299
         ib5=ibb5(iii)
         ESHORT=ESHORT+asr5(i,ib5)

      enddo

      I1=IIII-11
      if(i1.lt.1) i1=1
      i2=JJJJ
      if(i2.gt.lenf-12) i2=lenf-12

      do i=i1,i2
         j=i+4
         k=i+8
         wx1=x(j)-x(i)
         wy1=y(j)-y(i)
         wz1=z(j)-z(i)
         wx2=x(k)-x(j)
         wy2=y(k)-y(j)
         wz2=z(k)-z(j)
         kk=k+4
         wx3=x(kk)-x(k)
         wy3=y(kk)-y(k)
         wz3=z(kk)-z(k)
c
c     Prevents too clustered folding (too short secondary
c     structure fragments)
c
         if(wx1*wx2+wy1*wy2+wz1*wz2.lt.0) then
            if(wx1*wx3+wy1*wy3+wz1*wz3.gt.0) then
               if(wx2*wx3+wy2*wy3+wz2*wz3.lt.0) ESHORT=ESHORT+ESC1
            endif
         endif


         ir1=wx1*wx1+wy1*wy1+wz1*wz1
         ir2=wx2*wx2+wy2*wy2+wz2*wz2
c
c     Additional stifness of secondary structure elements
c

         if(ir1.lt.33.and. ir2.lt.33) then
            if(sec(j).ne.4) then
               wx3=wx2-wx1
               wy3=wy2-wy1
               wz3=wz2-wz1
               ir=wx3*wx3+wy3*wy3+wz3*wz3
               if(ir.lt.19)  ESHORT=ESHORT-ESC2
            endif
         else
            if(sec(j).ne.2) then
               if(ir1.gt.32.and.ir2.gt.32) then
                  if((ir1+ir2).gt.108) ESHORT=ESHORT-ESC2
               endif
            endif
         endif
      enddo

      ESHORT=ESHORT

      RETURN
      END

c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
      FUNCTION EPAIRN(ISTART,IEND)
      IMPLICIT REAL*8 (A-H)
      IMPLICIT INTEGER (I-Z)

      PARAMETER(NDIM=500)
      PARAMETER(NBOX=150)
      PARAMETER (IHCB=26)

      common/icutt/IRESA,icutp(0:19,0:19),icuta(0:19,0:19),
     *     icutm(0:19,0:19)
      COMMON /SEQE/  SEQ(ndim), SEC(ndim), ENVIR(0:5,0:5,0:5,0:19)
      COMMON /CHAIN/ ICA(-2:NDIM), x(-2:NDIM), Y(-2:NDIM), Z(-2:NDIM)
      COMMON /pair/ apba(ndim,ndim),apbp(ndim,ndim)
      COMMON /pairm/apbm(ndim,ndim)

      COMMON /LENGTHS/ LENF2,LENF1,LENF,LENGTH2(646)

      COMMON /SHORTE/ asr3(ndim,4),asr4(ndim,14),asr5(ndim,8),
     $     asr2(ndim,3),ESC,
     *     ibb2(0:199),ibb4(-100:103),ibb5(0:300),ibb3(0:199)

      COMMON/U/ NCIJA(ndim),NCIJP(ndim),NCIJDA(1600),NCIJDM(1600),
     *     NCIJDP(1600) ,NOA(ndim),NOP(ndim), NCIJM(ndim), NOM(ndim)

      COMMON/WORK/ xyz(NBOX,NBOX,NBOX)
      COMMON vx(646),vy(646),vz(646)
      COMMON /HBN/ EHBOND, NHBN(ndim), prod(646,646)
      COMMON /RES/ EREP,EREST,APHOBI,MRES(ndim),KRES(ndim,50)
      COMMON /NBN/ NBX(42),NBY(42),NBZ(42)
      DIMENSION ICONT(0:100)
      LOGICAL PHOB, PHIL
      COMMON /KD/ PHOB(ndim), PHIL(ndim), nkbn(ndim)
      COMMON/CA/cgs(0:20)

      COMMON/CAC/CAR(646*646*3)
      COMMON /HB/ HBR(646*646*3)

      COMMON /REHAB/NHPLUS,NHMIN,IGO(1600),ITO(1600)
      COMMON /RENAB/NBPLUS,NCIJ(ndim),NCIJD(1600),INTO(1600)
      COMMON /CENTER/  eone(ndim,0:100), ACRIT, compress, amid, mid
      COMMON /BURAK/ EBUREK(ndim,0:24)
      COMMON /SHELL1/ SHELL(680,24)
      COMMON /SHELL2/ SPHERE(NDIM,24), SPHEREM(NDIM,24)
      COMMON /SHELL3/ KPB(NDIM,NDIM), NVECTOR(-5:5,-5:5,-5:5)
      COMMON /CCCC/ ECORECT, EQ3(0:19,0:19,0:19)
      LOGICAL NOTUSED(1600)
      COMMON /REP/erepp(9:26,0:19,0:19),erepa(9:26,0:19,0:19),
     *     erepm(9:26,0:19,0:19)
c
      EPAIRN=0.
      NHPLUS=0
      EHBOND2=EHBOND/2.0D0

      DO K=ISTART,IEND
         icont(0)=0
         kx=x(k)
         ky=y(k)
         kz=z(k)

         CALL SET(kx,ky,kz,k)
c
c     Detect the contacts for the k-th residue
c

         NCONT=0
         do 1000 i=1,42
            ix=kx+NBX(i)
            iy=ky+NBY(i)
            iz=kz+NBZ(i)
            iii=XYZ(ix,iy,iz)
            if(iii.gt.0) then
               if(iabs(iii-k).gt.1) then
c     detect each contact only once
                  do kk=0,ncont
                     if(iii.eq.icont(kk)) go to 1000
                  enddo
                  NCONT=NCONT+1
                  ICONT(NCONT)=iii
               endif
            endif
 1000    continue

         if(NCONT.gt.0) THEN
            ican=ica(k-1)
            icak=ica(k)
            wx=vx(ican)-vx(icak)
            wy=vy(ican)-vy(icak)
            wz=vz(ican)-vz(icak)
            awx=sqrt(dble(wx*wx+wy*wy+wz*wz)+0.0001D0)

            fk=cgs(seq(k))
            iica=ica(k-1)*1938+ica(k)*3-1941
            bsx=(dble(kx)+CAR(iica+1)*fk)
            bsy=(dble(ky)+CAR(iica+2)*fk)
            bsz=(dble(kz)+CAR(iica+3)*fk)
c
c     !!!!        Hydrogen bonding of the k-th residue *********
c
            hx=HBR(iica+1)
            hy=HBR(iica+2)
            hz=HBR(iica+3)

            DO i=1,NCONT
               j=ICONT(i)
               jx=x(j)
               jy=y(j)
               jz=z(j)
               jjx=jx-kx
               jjy=jy-ky
               jjz=jz-kz
               ir=jjx*jjx + jjy*jjy + jjz*jjz

               if(ir.lt.IHCB) then

c     + add the backbone effect
c
c     Build a pseudo-backbone coordinates of j-th
                  fj=cgs(seq(j))
                  iica=ica(j-1)*1938+ica(j)*3-1941
                  asx=(dble(jx)+CAR(iica+1)*fj)
                  asy=(dble(jy)+CAR(iica+2)*fj)
                  asz=(dble(jz)+CAR(iica+3)*fj)
                  ajsx=asx-dble(kx)
                  ajsy=asy-dble(ky)
                  ajsz=asz-dble(kz)
                  aisr=ajsx*ajsx+ajsy*ajsy+ajsz*ajsz
                  if(aisr.lt.9.0D0) EPAIRN=EPAIRN+EREP

c     + add the backbone effect

                  ajsy=bsy-dble(jy)
                  ajsx=bsx-dble(jx)
                  ajsz=bsz-dble(jz)
                  aisr=ajsx*ajsx+ajsy*ajsy+ajsz*ajsz
                  if(aisr.lt.9.0D0) EPAIRN=EPAIRN+EREP

                  amx=asx-bsx
                  amy=asy-bsy
                  amz=asz-bsz
                  ar= amx*amx +amy*amy +amz*amz
                  if(ar.lt.8.0D0.AND.iabs(k-j).gt.3) EPAIRN=EPAIRN+EREP

                  ii=ica(j-1)
                  jj=ica(j)
c
c     Angular term for tertiary interactions
c
                  mx=vx(ii)-vx(jj)
                  my=vy(ii)-vy(jj)
                  mz=vz(ii)-vz(jj)
                  ip = wx*mx + wy*my + wz*mz
                  akx=sqrt(float(mx*mx+my*my+mz*mz)+0.0001)
                  akx=float(ip)/akx/awx
                  sj=seq(j)
                  sk=seq(k)
                  if(akx.gt.0.5) then
                     if(ir.lt.icutp(sj,sk)) then
                        NCIJ(k)=NCIJ(k)+1
                        NCIJ(j)=NCIJ(j)+1
                        NBPLUS=NBPLUS+1
                        INTO(NBPLUS)=j
                        NCIJD(NBPLUS)=1
                        NCIJDP(NBPLUS)=1
                        NCIJDA(NBPLUS)=0
                        NCIJDM(NBPLUS)=0
                        NCIJP(j)=NCIJP(j)+1
                        NBPLUS=NBPLUS+1
                        INTO(NBPLUS)=k
                        NCIJD(NBPLUS)=1
                        NCIJDP(NBPLUS)=1
                        NCIJDA(NBPLUS)=0
                        NCIJDM(NBPLUS)=0
                        NCIJP(k)=NCIJP(k)+1

                        kpoints=kpb(k,j)
c     kpb(k,j) coverage of k by j
                        nr=NVECTOR(jjx,jjy,jjz)
                        do kp=1,kpoints
                           kpp=shell(NR,kp)
                           sphere(k,kpp)=sphere(k,kpp)+1
                           spherem(k,kpp)=spherem(k,kpp)+1
                        enddo

                        jpoints=kpb(j,k)
c     kpb(j,k) coverage of j by k
                        nr=NVECTOR(-jjx,-jjy,-jjz)
                        do jp=1,jpoints
                           jpp=shell(NR,jp)
                           sphere(j,jpp)=sphere(j,jpp)+1
                           spherem(j,jpp)=spherem(j,jpp)+1
                        enddo

                        if(erepp(ir,sj,sk).gt.0.01) then
                           EPAIRN=EPAIRN+erepp(ir,sj,sk)
                        else
                           EPAIRN=EPAIRN+apbp(j,k)
                        endif
                     endif
                  else

                     if(akx.lt.-0.5) then
                        if(ir.lt.icuta(sj,sk))then
                           NCIJ(k)=NCIJ(k)+1
                           NCIJ(j)=NCIJ(j)+1
                           NBPLUS=NBPLUS+1
                           INTO(NBPLUS)=j
                           NCIJD(NBPLUS)=1
                           NCIJDA(NBPLUS)=1
                           NCIJDP(NBPLUS)=0
                           NCIJDM(NBPLUS)=0
                           NCIJA(j)=NCIJA(j)+1
                           NBPLUS=NBPLUS+1
                           INTO(NBPLUS)=k
                           NCIJD(NBPLUS)=1
                           NCIJDA(NBPLUS)=1
                           NCIJDP(NBPLUS)=0
                           NCIJDM(NBPLUS)=0
                           NCIJA(k)=NCIJA(k)+1

                           kpoints=kpb(k,j)
c     kpb(k,j) coverage of k by j
                           nr=NVECTOR(jjx,jjy,jjz)
                           do kp=1,kpoints
                              kpp=shell(NR,kp)
                              sphere(k,kpp)=sphere(k,kpp)+1
                              spherem(k,kpp)=spherem(k,kpp)+1
                           enddo

                           jpoints=kpb(j,k)
c     kpb(j,k) coverage of j by k
                           nr=NVECTOR(-jjx,-jjy,-jjz)
                           do jp=1,jpoints
                              jpp=shell(NR,jp)
                              sphere(j,jpp)=sphere(j,jpp)+1
                              spherem(j,jpp)=spherem(j,jpp)+1
                           enddo
                           if(erepa(ir,sj,sk).gt.0.01) then
                              EPAIRN=EPAIRN+erepa(ir,sj,sk)
                           else
                              EPAIRN=EPAIRN+apba(j,k)
                           endif
                        endif
                     else

                        if(ir.lt.icutm(sj,sk)) then
                           NCIJ(k)=NCIJ(k)+1
                           NCIJ(j)=NCIJ(j)+1
                           NBPLUS=NBPLUS+1
                           INTO(NBPLUS)=j
                           NCIJD(NBPLUS)=1
                           NCIJDP(NBPLUS)=0
                           NCIJDA(NBPLUS)=0
                           NCIJDM(NBPLUS)=1
                           NCIJM(j)=NCIJM(j)+1
                           NBPLUS=NBPLUS+1
                           INTO(NBPLUS)=k
                           NCIJD(NBPLUS)=1
                           NCIJDP(NBPLUS)=0
                           NCIJDA(NBPLUS)=0
                           NCIJDM(NBPLUS)=1
                           NCIJM(k)=NCIJM(k)+1

                           kpoints=kpb(k,j)
c     kpb(k,j) coverage of k by j
                           nr=NVECTOR(jjx,jjy,jjz)
                           do kp=1,kpoints
                              kpp=shell(NR,kp)
                              sphere(k,kpp)=sphere(k,kpp)+1
                              spherem(k,kpp)=spherem(k,kpp)+1
                           enddo

                           jpoints=kpb(j,k)
c     kpb(j,k) coverage of j by k
                           nr=NVECTOR(-jjx,-jjy,-jjz)
                           do jp=1,jpoints
                              jpp=shell(NR,jp)
                              sphere(j,jpp)=sphere(j,jpp)+1
                              spherem(j,jpp)=spherem(j,jpp)+1
                           enddo

                           if(erepm(ir,sj,sk).gt.0.01) then
                              EPAIRN=EPAIRN+erepm(ir,sj,sk)
                           else
                              EPAIRN=EPAIRN+apbm(j,k)
                           endif
                        endif
                     endif
                  endif

                  if(ir.lt.23) then
c
c     3-body **************************************************
c

                     IF(ir.lt.icutp(sk,sj)) then
                        IF(NCONT.GT.1) then
                           IF(i.gt.1) THEN
                              DO in=1,i-1
                                 i3=icont(in)
                                 IF(iabs(i3-j).gt.1) then
                                    si3=seq(i3)
                                    ikr=(kx-x(i3))**2+(ky-y(i3))**2+
     $                                   (kz-z(i3))**2
                                    if(ikr.lt.icutp(sk,si3)) then
                                       ijr=(jx-x(i3))**2+
     $                                      (jy-y(i3))**2+
     $                                      (jz-z(i3))**2
                                       if(ijr.lt.icutp(sj,si3)) then
                                          if(iabs(k-j)+iabs(i3-j)+
     $                                         iabs(i3-k).GT.7) then
                                             EPAIRN=EPAIRN+
     $                                            EQ3(sj,sk,si3)
                                          endif
                                       endif
                                    endif
                                 ENDIF
                              ENDDO
                           ENDIF
                        ENDIF
                     ENDIF
c
c     3-body **************************************************
c
                     if(ip.gt.0) then
                        if(ar.lt.20.0) then

c     parallel (compute "hydrogen bonds") !!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c     SECONDARY BIAS  H-bond SELECTIONS     !!!!!!!!!!!!!!!!!!!!!!!!!
c
                           if(iabs(j-k).lt.3) go to 987
                           if(sec(j)*sec(k).eq.8) Go to 987
                           if(iabs(j-k).gt.3) then
                              if(sec(j).eq.2.AND.sec(k).eq.2) go to 987
                              if(iabs(j-k).eq.4) Go to 987
                           endif

                           iica=ica(j-1)*1938+ica(j)*3-1941

                           hjx=HBR(iica+1)
                           hjy=HBR(iica+2)
                           hjz=HBR(iica+3)

                           coop= hjx*hx+hjy*hy+hjz*hz
                           if(abs(coop).gt.6.0) then

                              HPROP=0.0
                              iir=(x(k+1)-x(k-1))*(x(j+1)-x(j-1))+
     *                             (y(k+1)-y(k-1))*(y(j+1)-y(j-1))+
     *                             (z(k+1)-z(k-1))*(z(j+1)-z(j-1))
                              if(iir.gt.0) then
                                 jrp=(x(k+1)-x(j+1))**2+
     $                                (y(k+1)-y(j+1))**2+
     $                                (z(k+1)-z(j+1))**2
                                 if(jrp.lt.23)  then
                                    jrm=(x(k-1)-x(j-1))**2+
     $                                   (y(k-1)-y(j-1))**2+
     $                                   (z(k-1)-z(j-1))**2
                                    if(jrm.lt.23) HPROP=EHBOND
                                 endif
                              else
                                 jrp=(x(k+1)-x(j-1))**2+
     $                                (y(k+1)-y(j-1))**2+
     $                                (z(k+1)-z(j-1))**2
                                 if(jrp.lt.23)  then
                                    jrm=(x(k-1)-x(j+1))**2+
     $                                   (y(k-1)-y(j+1))**2+
     $                                   (z(k-1)-z(j+1))**2
                                    if(jrm.lt.23) HPROP=EHBOND
                                 endif
                              endif

                              anx=amx-hx
                              any=amy-hy
                              anz=amz-hz
                              air=anx*anx+any*any+anz*anz
                              if(air.lt.4.0) then
                                 NHBN(k)=NHBN(k)+1
                                 NHPLUS=NHPLUS+1
                                 ITO(NHPLUS)=k
                                 EPAIRN=EPAIRN+(EHBOND+HPROP)
                              endif

                              anx=amx+hjx
                              any=amy+hjy
                              anz=amz+hjz
                              air=anx*anx+any*any+anz*anz
                              if(air.lt.4.0) then
                                 NHBN(j)=NHBN(j)+1
                                 NHPLUS=NHPLUS+1
                                 ITO(NHPLUS)=j
                                 EPAIRN=EPAIRN+(EHBOND+HPROP)
                              endif

                              anx=amx+hx
                              any=amy+hy
                              anz=amz+hz
                              air=anx*anx+any*any+anz*anz
                              if(air.lt.4.0) then
                                 NHBN(k)=NHBN(k)+1
                                 NHPLUS=NHPLUS+1
                                 ITO(NHPLUS)=k
                                 EPAIRN=EPAIRN+(EHBOND+HPROP)
                              endif

                              anx=amx-hjx
                              any=amy-hjy
                              anz=amz-hjz
                              air=anx*anx+any*any+anz*anz
                              if(air.lt.4.0) then
                                 NHBN(j)=NHBN(j)+1
                                 NHPLUS=NHPLUS+1
                                 ITO(NHPLUS)=j
                                 EPAIRN=EPAIRN+(EHBOND+HPROP)
                              endif
                           endif

 987                       continue

                        endif
                     endif
                  endif

               ENDIF
            ENDDO
         ENDIF
      ENDDO


c
c     BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
c     BURIAL COMPUTATIONS - SURFACE COVERAGE
c     Differnece of energy (not state energy)
c
      ECORECT=0.
      IF(NBPLUS.GT.0) THEN
         do i=1,NBPLUS
            k=into(i)
            NOTUSED(k)=.TRUE.
         enddo

         DO i=1,NBPLUS
            k=into(i)

            npoints=0
            olp=0
            if(NOTUSED(k)) then

c
c     ENVIROMENTAL PROFILE UPGRADE
c
               ia=NCIJA(k)/2
               ip=NCIJP(k)/2
               im=NCIJM(k)/2
               is=seq(k)
               iaa=NOA(k)/2
               ipp=NOP(k)/2
               imm=NOM(k)/2
               EPAIRN=EPAIRN+ENVIR(ia,im,ip,is)-ENVIR(iaa,imm,ipp,is)

               do j=1,24
                  if(SPHERE(k,j).gt.0) npoints=npoints+1
                  if((SPHERE(k,j)-SPHEREM(k,j)).gt.0) olp=olp+1
               enddo
               ECORECT=ECORECT+EBUREK(k,olp)
               EPAIRN=EPAIRN+EBUREK(k,npoints)
            endif
            NOTUSED(k)=.FALSE.

         ENDDO
         EPAIRN=EPAIRN-ECORECT
      ENDIF


      DO K=ISTART,IEND
         kx=x(k)
         ky=y(k)
         kz=z(k)
         CALL REM(kx,ky,kz)
      ENDDO
c
c     BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
c

      RETURN
      END



      FUNCTION ECENTO(istart,iend)
      IMPLICIT REAL*8 (A-H)
      IMPLICIT INTEGER (i-z)
c     note that the contribution from this subroutine
c     is not included in total energy counter !!!!!!!

      PARAMETER(NDIM=500)
      PARAMETER(NBOX=150)
      PARAMETER (IHCB=26)

      common/icutt/IRESA,icutp(0:19,0:19),icuta(0:19,0:19),
     *     icutm(0:19,0:19)
      COMMON vx(646),vy(646),vz(646)
      COMMON /SEQE/  SEQ(ndim), SEC(ndim), ENVIR(0:5,0:5,0:5,0:19)

      COMMON/BUR/IBUR(0:5),MINUSB(0:5),PLUSB(0:5),IBTARG(0:5)
      COMMON /RES/ EREP,EREST,APHOBI,MRES(ndim),KRES(ndim,50)
      COMMON /RNN/ MRESH(ndim),KRESH(ndim,4)
      COMMON /CHAIN/ ICA(-2:NDIM), x(-2:NDIM), Y(-2:NDIM), Z(-2:NDIM)
      COMMON/RCN/MRESA(ndim),KRESA(ndim,4),iresta(ndim),jresta(ndim)
      COMMON /CENTER/  eone(ndim,0:100), ACRIT, compress, amid, mid

      COMMON /LENGTHS/ LENF2,LENF1,LENF,LENGTH2(646)

      COMMON/CA/cgs(0:20)
      COMMON/CAC/CAR(646*646*3)
      COMMON/CAH/ICAH(646*646*3)
      COMMON/CNH/ICNH(646*646*3)

      COMMON/kdh/  eonekd(0:19)

c
c     centrosymmetric one body potential
c

      ECENTO=0.
      do i=0,5
         ECENTO=ECENTO+compress*(MAX(2,IABS(IBUR(i)-IBTARG(i)))-2)
         MINUSB(i)=0
      enddo

      do i=istart, iend
         kx=x(i)
         ky=y(i)
         kz=z(i)

         ix=kx-mid
         iy=ky-mid
         iz=kz-mid

         aract2=dble(ix*ix+iy*iy+iz*iz)
         aract=sqrt(aract2)

         ff=aract/acrit
         ika=int(ff/0.333333D0)
         ECENTO=ECENTO+eone(i,ika)
         if(ika.gt.5) ika=5

         IBUR(ika)=IBUR(ika)-1
         MINUSB(ika)=MINUSB(ika)+1

         MM=MRES(i)
         if(MM.GT.0) THEN

            do k=1,MM
               j=kres(i,k)
               if(j.lt.i.OR.j.gt.iend) then

                  xj=x(j)
                  yj=y(j)
                  zj=z(j)
                  ir=(xj-kx)*(xj-kx)+(yj-ky)*(yj-ky)+(zj-kz)**2
                  if(ir.gt.IHCB)
     $                 ECENTO=ECENTO+EREST*sqrt(dble(ir-ihcb))
               endif
            enddo
         endif
      enddo
c
      do i=istart-1, iend+1
         MM=MRESH(i)
         if(MM.GT.0) THEN
            kx=x(i)
            ky=y(i)
            kz=z(i)
            iica=ica(i-1)*1938+ica(i)*3-1941
            ix=kx+ICNH(iica+1)
            iy=ky+ICNH(iica+2)
            iz=kz+ICNH(iica+3)
            do k=1,MM
               j=kresh(i,k)
               if(j.lt.i.OR.j.gt.iend+1) then
                  iica=ica(j-1)*1938+ica(j)*3-1941
                  xj=x(j)+ICNH(iica+1)
                  yj=y(j)+ICNH(iica+2)
                  zj=z(j)+ICNH(iica+3)
                  ir=(xj-ix)*(xj-ix)+(yj-iy)*(yj-iy)+(zj-iz)**2
                  if(ir.gt.16) ECENTO=ECENTO+EREST*sqrt(dble(ir)-3.0D0)
               endif
            enddo
         endif

         MM=MRESA(i)
         if(MM.GT.0) THEN
            kx=x(i)
            ky=y(i)
            kz=z(i)

            iica=ica(i-1)*1938+ica(i)*3-1941
            if(iresta(i).ne.i) then
c     means NH
               ix=kx+ICNH(iica+1)
               iy=ky+ICNH(iica+2)
               iz=kz+ICNH(iica+3)
            else
c     means CaH
               ix=kx+ICAH(iica+1)
               iy=ky+ICAH(iica+2)
               iz=kz+ICAH(iica+3)
            endif
            do k=1,MM
               j=kresa(i,k)
               if(j.lt.i.OR.j.gt.iend+1) then
                  iica=ica(j-1)*1938+ica(j)*3-1941
                  if(iresta(j).ne.j) then
c     means NH
                     xj=x(j)+ICNH(iica+1)
                     yj=y(j)+ICNH(iica+2)
                     zj=z(j)+ICNH(iica+3)
                  else
                     xj=x(j)+ICAH(iica+1)
                     yj=y(j)+ICAH(iica+2)
                     zj=z(j)+ICAH(iica+3)
                  endif
                  ir=(xj-ix)*(xj-ix)+(yj-iy)*(yj-iy)+(zj-iz)**2
                  if(ir.gt.16) ECENTO=ECENTO+EREST*sqrt(dble(ir)-3.0D0)
               endif
            enddo
         endif

      enddo



      RETURN
      END

c
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c

      FUNCTION ECENTN(istart,iend)
      IMPLICIT REAL*8 (A-H)
      IMPLICIT INTEGER (i-z)
c     note that the contribution from this subroutine
c     is not included in total energy counter !!!!!!!

      PARAMETER(NDIM=500)
      PARAMETER (IHCB=26)
      PARAMETER(NBOX=150)

      common/icutt/IRESA,icutp(0:19,0:19),icuta(0:19,0:19),
     *     icutm(0:19,0:19)
      COMMON vx(646),vy(646),vz(646)
      COMMON /SEQE/  SEQ(ndim), SEC(ndim), ENVIR(0:5,0:5,0:5,0:19)

      COMMON/BUR/IBUR(0:5),MINUSB(0:5),PLUSB(0:5),IBTARG(0:5)
      COMMON /RES/ EREP,EREST,APHOBI,MRES(ndim),KRES(ndim,50)
      COMMON /RNN/ MRESH(ndim),KRESH(ndim,4)
      COMMON/RCN/MRESA(ndim),KRESA(ndim,4),iresta(ndim),jresta(ndim)
      COMMON /CHAIN/ ICA(-2:NDIM), x(-2:NDIM), Y(-2:NDIM), Z(-2:NDIM)
      COMMON /CENTER/  eone(ndim,0:100), ACRIT, compress, amid, mid

      COMMON /LENGTHS/ LENF2,LENF1,LENF,LENGTH2(646)

      COMMON/CA/cgs(0:20)
      COMMON/CAC/CAR(646*646*3)
      COMMON/CAH/ICAH(646*646*3)
      COMMON/CNH/ICNH(646*646*3)
      COMMON/kdh/  eonekd(0:19)
c
c     centrosymmetric one body potential
c

      ECENTN=0.
      do i=0,5
         PLUSB(i)=0
      enddo


      do i=istart, iend
         kx=x(i)
         ky=y(i)
         kz=z(i)

         ix=kx-mid
         iy=ky-mid
         iz=kz-mid

         aract2=dble(ix*ix+iy*iy+iz*iz)
         aract=sqrt(aract2)

         ff=aract/acrit
         ika=int(ff/0.333333D0)

         ECENTN=ECENTN+eone(i,ika)
         if(ika.gt.5) ika=5

         IBUR(ika)=IBUR(ika)+1
         PLUSB(ika)=PLUSB(ika)+1


         MM=MRES(i)
         if(MM.GT.0) THEN

            do k=1,MM
               j=kres(i,k)
               if(j.lt.i.OR.j.gt.iend) then
                  xj=x(j)
                  yj=y(j)
                  zj=z(j)
                  ir=(xj-kx)*(xj-kx)+(yj-ky)*(yj-ky)+(zj-kz)**2
                  if(ir.gt.IHCB)
     $                 ECENTN=ECENTN+EREST*sqrt(dble(ir-ihcb))
               endif
            enddo
         endif
      enddo
c
c
      do i=istart-1, iend+1
         MM=MRESH(i)
         if(MM.GT.0) THEN
            kx=x(i)
            ky=y(i)
            kz=z(i)
            iica=ica(i-1)*1938+ica(i)*3-1941
            ix=kx+ICNH(iica+1)
            iy=ky+ICNH(iica+2)
            iz=kz+ICNH(iica+3)
            do k=1,MM
               j=kresh(i,k)
               if(j.lt.i.OR.j.gt.iend+1) then
                  iica=ica(j-1)*1938+ica(j)*3-1941

                  xj=x(j)+ICNH(iica+1)
                  yj=y(j)+ICNH(iica+2)
                  zj=z(j)+ICNH(iica+3)
                  ir=(xj-ix)*(xj-ix)+(yj-iy)*(yj-iy)+(zj-iz)**2
                  if(ir.gt.16) ECENTN=ECENTN+EREST*sqrt(dble(ir)-3.0D0)
               endif
            enddo
         endif

         MM=MRESA(i)
         if(MM.GT.0) THEN
            kx=x(i)
            ky=y(i)
            kz=z(i)
            iica=ica(i-1)*1938+ica(i)*3-1941
            if(iresta(i).ne.i) then
c     means NH
               ix=kx+ICNH(iica+1)
               iy=ky+ICNH(iica+2)
               iz=kz+ICNH(iica+3)
            else
c     means CaH
               ix=kx+ICAH(iica+1)
               iy=ky+ICAH(iica+2)
               iz=kz+ICAH(iica+3)
            endif
            do k=1,MM
               j=kresa(i,k)
               if(j.lt.i.OR.j.gt.iend+1) then
                  iica=ica(j-1)*1938+ica(j)*3-1941
                  if(iresta(j).ne.j) then
c     means NH
                     xj=x(j)+ICNH(iica+1)
                     yj=y(j)+ICNH(iica+2)
                     zj=z(j)+ICNH(iica+3)
                  else
c     means CaH
                     xj=x(j)+ICAH(iica+1)
                     yj=y(j)+ICAH(iica+2)
                     zj=z(j)+ICAH(iica+3)
                  endif
                  ir=(xj-ix)*(xj-ix)+(yj-iy)*(yj-iy)+(zj-iz)**2
                  if(ir.gt.16) ECENTN=ECENTN+EREST*sqrt(dble(ir)-3.0D0)
               endif
            enddo
         endif

      enddo

      do i=0,5
         ECENTN=ECENTN+compress*(MAX(2,IABS(IBUR(i)-IBTARG(i)))-2)
      enddo

      RETURN
      END

