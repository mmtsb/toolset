c     Fold.f: The source for Lattice Modeling Tool
c     (Andrzej Kolinski, May 1999, May 2000)
c     Exclusively for academic users.
c
c     modifications by Michael Feig, 1999, 2000
c     for fragment modeling
c
c     updated May 2000
c
c     positional restraints added 
c     August 2000 by Michael Feig
c
c     reorganization of arrays to improve SGI performance
c     September 2000, Michael Feig
c
c     I/O rewritten, mfold/mfragfold combined
c     October 2000, Michael Feig
c
c     mfold/mfoldrex combined
c     March 2001, Michael Feig

c     Applications:
C     (1) Ab initio folding
c     (2) Protein dynamics and thermodynamics
c     (3) Structure assembly from sparse experimental
c     data
c
c
c     program        requires the following input files:
c
c     INPUT  (a generic examle below)
c
c     6141 10  300  20 (Random, Ncycle, Icycle, Tcycle)
c     2.0   1.  4.0 1.0  (Temp1, Temp2, soft_core, generic_central)
c     1.0 1.75 0.125 -1. 0.5 (gen,pair,KD_centr,H-bond,short)
c     0.125  0.125 0.125  (burial 24-neib, multibody, 3-body)
c     0   0.0 (Nres, strength) Restraints - optional
c     0
c     0
c
c     SEQ,  CHAIN646,  + Force field datafiles (Andrzej Kolinski,
c     Piotr Rotkiewicz and Jeffrey Skolnick)
c     Other stuff as Kyte_Doolittle scale in data statements.
c
c     Adjust properly NDIM and NBOX values
c
c     contact A. Kolinski for updates: kolinski@chem.uw.edu.pl
c     phone: +48-22-8220211 ext. 320
c
c     See supplementary tools: chain.f, rmsCA145.f. Native.f
c     Seqgen.f and others

      IMPLICIT REAL*8 (A-H)
      IMPLICIT INTEGER (I-Z)

      character*5 struct(5)
      character*3 aa(-1:20), NAME
      character*12 text

      PARAMETER(NDIM=600)
      PARAMETER(NBOX=100)

      character*256 datdir

      common/icutt/IRESA,icutp(0:19,0:19),icuta(0:19,0:19),
     *     icutm(0:19,0:19)
      REAL*8 rrand

      LOGICAL LOOK,GOODC(646,646)

      COMMON/BUR/IBUR(0:5),MINUSB(0:5),PLUSB(0:5),IBTARG(0:5)

      COMMON /WORK/ xyz(NBOX,NBOX,NBOX)
      COMMON /CHAIN/ ICA(-2:NDIM), x(-2:NDIM), Y(-2:NDIM), Z(-2:NDIM)
      COMMON /SEQE/  SEQ(ndim), SEC(ndim), ENVIR(0:5,0:5,0:5,0:19)
      COMMON /pair/ apba(ndim,ndim),apbp(ndim,ndim)
      COMMON /pairm/apbm(ndim,ndim)
      COMMON /THREE/ ICONF(646,646)

      REAL*8 pres,rfac
      INTEGER xor,yor,zor,xshrf,yshrf,zshrf
      LOGICAL dorest
      COMMON /RESTR/ rfac(21,21,21),
     $     pres(NDIM),
     $     xor(NDIM),yor(NDIM),zor(NDIM),
     $     imx(NBOX*2), dorest,
     $     xshrf,yshrf,zshrf

      INTEGER FXSTART,FXEND,FXLEN,FXLEN5
      LOGICAL FXNT,FXCT
      COMMON /LENGTHS/ LENF2,LENF1,LENF,LENGTH2(646),
     $     FXSTART,FXEND,FXLEN,FXLEN5,FXNT,FXCT,nrest

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
      COMMON vx(646),vy(646),vz(646)
      DIMENSION vector(-5:5,-5:5,-5:5)
      DIMENSION map(40)
      DIMENSION apabla(0:19,0:19), apablp(0:19,0:19)
      DIMENSION eoinp(0:19,0:100), apablm(0:19,0:19)
      DIMENSION xt(ndim),yt(ndim),zt(ndim), nkb(0:19)
      COMMON /HBN/ EHBOND, NHBN(ndim), prod(646,646)
      LOGICAL MAPR(ndim,ndim)
      COMMON /RES/ EREP,EREST,APHOBI,MRES(ndim),KRES(ndim,50)
      COMMON/RCN/MRESA(ndim),KRESA(ndim,4),iresta(ndim),jresta(ndim)
      DIMENSION irest(ndim),jrest(ndim),ICONT(0:100),neib(ndim,ndim)
      DIMENSION iresth(ndim),jresth(ndim)
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
      DIMENSION ECOVER(0:ndim, 0:24)
      DIMENSION KPBA(-1:19,-1:19), STAMAP(ndim,ndim)
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

      COMMON/umbr/ rgsum,orgsum,nrgsum,tsrg,srgsum,
     $             rgyr,krg,rho,krho,scont,
     $             oscont,nscont,sscont
      REAL*8 rgsum,orgsum,nrgsum,tsrg,srgsum,
     $     rgyr,krg,rho,krho,scont,oscont,
     $     nscont,sscont

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

      CHARACTER*80 jobid
      CHARACTER*80 svname
      CHARACTER*80 svport
      CHARACTER*80 svid
      CHARACTER*80 svdir
      INTEGER mpflag,mpsend
      REAL*8 newtemp
      REAL*8 lstene,lstrg,lstrho
      LOGICAL rex
      LOGICAL start

c     0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9
c     data nkb/3,3,4,5,6,5,7,3,7,5,5,7,4,5,6,6,6,7,7,8/

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

      open(unit=78,err=8543,
     $     FILE='monsster.datadir',status='old')
      read(78,'(A)') datdir
      close(78)
      goto 8544

 8543 datdir='/usr/local/mmtsb/data'

 8544 do i=1,256
         if (ichar(datdir(i:i)).ne.32) then
            ldatai=i
         endif
      enddo

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

      open(UNIT=12, FILE='monsster.job',STATUS='OLD',ERR=1173)

      rex=.true.
      read(12,*) lenf
      read(12,*) mpflag
      read(12,'(A)') jobid
      read(12,'(A)') svname
      read(12,'(A)') svport
      read(12,'(A)') svid
      read(12,'(A)') svdir
      mpsend=0
      close(12)

      goto 1174

 1173 rex=.false.

      OPEN(UNIT=10,FILE='monsster.init.chain',STATUS='OLD')
      read(10,*) lenf
      CLOSE(10)

      OPEN(UNIT=9,FILE='monsster.tra',STATUS='UNKNOWN')

 1174 continue

      OPEN(UNIT=6, FILE='monsster.output',STATUS='UNKNOWN')
      OPEN(UNIT=5, FILE='monsster.input',STATUS='OLD')
      OPEN(UNIT=7, FILE='monsster.seq',STATUS='OLD')
      OPEN(UNIT=8, FILE=datdir(1:ldatai)//'/PROFILE3',STATUS='OLD')
      OPEN(UNIT=14,FILE=datdir(1:ldatai)//'/S1234',STATUS='OLD')
      OPEN(UNIT=26,FILE=datdir(1:ldatai)//'/QUASI3',STATUS='OLD')
      OPEN(UNIT=27,FILE=datdir(1:ldatai)//'/ECOVERS_24',STATUS='OLD')
      OPEN(UNIT=28,FILE=datdir(1:ldatai)//'/BURIALS',STATUS='OLD')
      OPEN(UNIT=18,ERR=6643,
     $     FILE='monsster.restraints',STATUS='OLD')
      DOREST=.TRUE.
      GOTO 6644
 6643 DOREST=.FALSE.
 6644 CONTINUE
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
      write(6,*)
      Write(6,*) ' number of various contact distances',iii
      write(6,*)

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

      IF (LENF.LT.8) THEN
         WRITE(6,*) 'input chain is too short'
         CLOSE(6)
         STOP
      ENDIF

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

      do 121 i=2,lenf1
         read(7,707) k, NAME,SEC(I)
c 707     format(i5,3x,a3,2i5)
 707     format(i5,3x,a3,i5)
         do j=0,20
            if(NAME.eq.aa(j)) then
               SEQ(i)=j
               go to 121
            endif
         enddo
 121  continue

      do i=2,lenf1
         ii=seq(i)
         do j=2,lenf1
            stamap(j,i)=0
            jj=seq(j)
            KPB(i,j)=KPBA(ii,jj)
         enddo
c     backbone contribution
         KPB(i,1)=KPBA(ii,-1)
      enddo

c     read positional restraints
      if (dorest) then
         do i=1,lenf
            pres(i)=0.0
         enddo
         read(18,717) npres
         do i=1,npres
            read(18,718) k, apres
            if (k.lt.1 .or. k.gt.lenf-2) then
               write(*,*) 'restraint index ',k,' out of range: ',
     $              1,lenf-2
            endif
            pres(k+1)=apres
            
         enddo
 717     format(i5)
 718     format(i5,f9.5)
         close(unit=18)
      endif

      do i1=-10,10
         ds1=dble(i1)*dble(i1)
         do i2=-10,10
            ds2=dble(i2)*dble(i2)
            do i3=-10,10
               ds3=dble(i3)*dble(i3)
               dd=sqrt(ds1+ds2+ds3)
               if (dd.lt.10.0D0) then
                  rfac(i1+11,i2+11,i3+11)=dd*dd/10.0D0
               else
                  rfac(i1+11,i2+11,i3+11)=10.0D0
               endif
            enddo
         enddo
      enddo



      do i=1,NBOX*2
         iv=i-NBOX
         if (iv.lt.-10) then
            imx(i)=1
         else if(iv.gt.10) then
            imx(i)=21
         else
            imx(i)=iv+11
         endif
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


      READ(5,*) RANDOM,NCYCLE,PHOT,TCYCLE

      READ(5,*) FXSTART,FXEND

      READ(5,*) RGYR, KRG, RHO, KRHO

      READ(5,*) ATEMP1,ATEMP2, AAREP,compr

      READ(5,*) ESCO, ARLO, ENONE, EHBN, ARS

      READ(5,*) BUREK, ENSCAL, ES3

      IF (FXSTART.EQ.1) THEN
         FXSTART=0
      ENDIF

      IF (FXEND.EQ.0 .OR. FXEND.EQ.LENF2) THEN
         FXEND=LENF1
      ENDIF

      IF (FXSTART.EQ.0) THEN
         FXNT=.TRUE.
      ELSE
         FXNT=.FALSE.
      ENDIF

      IF (FXEND.EQ.LENF1) THEN
         FXCT=.TRUE.
      ELSE
         FXCT=.FALSE.
      ENDIF

      IF (FXSTART.GT.1) THEN
         FXSTART=FXSTART-2
      ENDIF

      IF (FXEND.LT.LENF1) THEN
         FXEND=FXEND+2
      ENDIF

      FXLEN=FXEND-FXSTART+1
      FXAL4=FXLEN-4.0-0.00001D0

      FXLEN5=FXLEN-5
      FXAL7=FXLEN-7.0-0.00001D0
      FXAL16=FXLEN-16

      WRITE(6,*)  ' OUTPUT FROM Monster6.f for 646 vectors 1998'
      WRITE(6,*)  ' uses 42 shell-scanning vectors June 09, 1997'
      Write(6,*)  '   **  USES PREDICTED OR KNOWN SECONDARY  **'

C
      WRITE(6,8020) RANDOM,NCYCLE, PHOT, TCYCLE

 8020 FORMAT(1X,/,
     *     1X,' RANDOM SEED =',I6,'    NUMBER OF CYCLES',3I6,//)

      IF (REX) THEN
         WRITE(6,*) ' VARIABLE TEMPERATURE REPLICA EXCHANGE'
      ELSE
         WRITE(6,*) ' TEMPERATURE RANGE: ',ATEMP1,ATEMP2
      ENDIF

      WRITE(6,*)
      WRITE(6,*) ' THE SHORT RANGE PARAMETERS: ', ESCO, ars
      WRITE(6,*) ' THE LONG RANGE PARAMETER:   ', ARLO
      WRITE(6,*) ' ONE BODY CENTROSYMMETRIC:   ', ENONE
      WRITE(6,*) ' THREE BODY POTENTIAL:       ', ES3
      WRITE(6,*) ' HYDROGEN BONDING                  ', EHBN
      WRITE(6,*) ' Compressing force                  ', compr
      WRITE(6,*) ' NEW BURIAL                  ', BUREK
      WRITE(6,*) ' ENVIROMENTAL DESCRIPTOR     ', ENSCAL
      IF (KRG.GT.0.000001D0) THEN
       WRITE(6,*) ' Radius of gyration umbrella at    ',rgyr 
       WRITE(6,*) '          force constant           ',krg
      ENDIF
      IF (KRHO.GT.0000001D0) THEN
       WRITE(6,*) ' Fraction native contacts umbrella at ',rho
       WRITE(6,*) '          force constant             ',krho
      ENDIF

c
c     READS Side group - side group contacts (NOE or other)
c
      READ(5,*) NREST, AREST
      WRITE(6,*)
      WRITE(6,*) 'STRENGTH & # of SG-SG RESTRAINTS:  ',AREST,NREST
      WRITE(6,*)
      if(NREST.gt.0 ) THEN
         do k=1, NREST
            read(5,*) IREST(K),JREST(k)
         enddo
      ENDIF

c
c     READS Ca - HN contacts (NOE)
c     WRITE(6,*)' FOR BACKBONE RESTRAINTS USE UP TO ONE PER RES.'
c
      READ(5,*) NRESTA
      IF(NRESTA.gt.0 ) THEN
         do k=1, NRESTA
            read(5,*) IRESTA(K),JRESTA(k)
         enddo
      ENDIF

      WRITE(6,*)
      WRITE(6,*) '# of CaH - HN RESTRAINTS:  ',NRESTA
      WRITE(6,*)

c
c     READS NH - HN contacts (NOE)
c
      READ(5,*) NRESTH
      IF(NRESTH.gt.0) then
         do k=1, NRESTH
            read(5,*) IRESTH(K),JRESTH(k)
         enddo
      endif

      WRITE(6,*)
      WRITE(6,*) '# of NH - HN RESTRAINTS:  ',NRESTH
      WRITE(6,*)

      if (dorest) then
         write(6,*)
         write(6,*) 'Using harmonic restraints from monsster.restraints'
         write(6,*)
         write(6,*)
      endif

c
c     !!!!!!!!!!!!!!!!!!!!!!!!  Add others when required
c

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
         write(6,*) ' In shell # ',i,'   # of residues ',ibtarg(i)
      enddo

      WRITE(6,*)

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
      write(6,*)
      write(6,*)
      write(6,*) ' the parameters of the shell model '
      write(6,*) ' l.units  S= ', Acrit, '   One body = ',enone
      write(6,*)
c

      MAXIM=NBOX
      MID=MAXIM/2
      amid=dble(mid)
      WRITE(6,*)' Monte Carlo BOX and NDIM:', MAXIM, NDIM

c
c     RESTRAINTS - SIMPLE CHANGE OF THE POTENTIAL
c
      if(nrest.gt.0) then
         do k=1,NREST
            i=irest(k)+1
            j=jrest(k)+1
            ii=MRES(i)+1
            jj=MRES(j)+1
            MRES(i)=ii
            MRES(j)=jj
            KRES(i,ii)=j
            KRES(j,jj)=i
            irest(k)=i
            jrest(k)=j
            MAPR(i,j)=.TRUE.
            MAPR(j,I)=.TRUE.
         enddo
      endif

      if(nresta.gt.0) then
         do k=1,NRESTA
            i=iresta(k)+1
            j=jresta(k)+1
            ii=MRESA(i)+1
            jj=MRESA(j)+1
            MRESA(i)=ii
            MRESA(j)=jj
            KRESA(i,ii)=j
            KRESA(j,jj)=i
            iresta(k)=i
            jresta(k)=j
         enddo
      endif

      if(nresth.gt.0) then
         do k=1,NRESTH
            i=iresth(k)+1
            j=jresth(k)+1
            ii=MRESH(i)+1
            jj=MRESH(j)+1
            MRESH(i)=ii
            MRESH(j)=jj
            KRESH(i,ii)=j
            KRESH(j,jj)=i
            iresth(k)=i
            jresth(k)=j
         enddo
      endif

      close(26)

      do i=2,lenf1
         ii=seq(i)
         nkbn(i)=nkb(ii)
         if(eonekd(ii).gt.0.0) then
            PHOB(i)=.TRUE.
            PHIL(i)=.FALSE.
         else
            PHOB(i)=.FALSE.
            PHIL(i)=.TRUE.
         endif
         do 128 j=2,lenf1

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
            if(NREST.gt.0) then
               if(seq(j).eq.3.and.seq(i).eq.3) then
c     the non-crosslinked CYS have a weaker interactions
                  if(.NOT.MAPR(i,j)) then
                     apba(i,j)=0.
                     apbm(i,j)=0.
                     apbp(i,j)=0.
                  endif
               endif
            endif


 128     continue
      enddo
      
      IF (rex) THEN
         atemp=1.0
      ELSE
         IF (TCYCLE.GT.1) THEN
            dtem=(atemp2-atemp1)/dble(TCYCLE-1)
         ELSE
            dtem=0.0D0
         ENDIF
         atemp=atemp1
      ENDIF
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

      do i=2,lenf1
         do j=2,lenf1
            apba(j,i)=arlo*apba(j,i)/atemp
            apbm(j,i)=arlo*apbm(j,i)/atemp
            apbp(j,i)=arlo*apbp(j,i)/atemp
         enddo
      enddo

c
c
      do i=2,lenf-2
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

      do i=2,lenf-3
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

      do i=2,lenf-4
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

      do i=2,lenf-5
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

c     *******  original temperature scaling finished  ********
c
C     ******************INITIAL CONFORMATION******************
c
      arand=(dble(RANDOM)/0.90721181D0)
      arand=arand-aint(arand)
      brand=0.590211D0
      crand=0.646017D0
      drand=211.02317D0

      write(6,*)
      write(6,*) ' the first random number is: ', rrand(ise)
      write(6,*)

      lstene=0.0D0
      lstrg=0.0D0
      lstrho=0.0D0

      ITEMP=0

 8877 atempo=atemp

      itemp=itemp+1

      IF (rex) THEN
         atemp=newtemp(svname,svport,svid,jobid,
     $        svdir,lstene,lstrg,lstrho,mpsend)
         if (atemp.lt.-9999.0) goto 8879

         call getbias(argyr,akrg,arho,akrho)
         if (argyr.gt.0.0D0) then
            rgyr=argyr
            krg=0.5D0*akrg
         endif
         if (arho.gt.0.0D0) then
            rho=arho
            krho=0.5D0*akrho
         endif

         write(*,*) 'rg: ',rgyr,' krg: ',krg

         if (atemp.gt.0.0 .and. itemp.gt.1) goto 8878
         if (atemp.lt.0.0) then
            atemp=-atemp
         endif
      ELSE
         if (itemp.gt.tcycle) goto 8879
         atemp=atemp1+dble(itemp-1)*dtem
         if (itemp.gt.1) goto 8878
      ENDIF

      OPEN(UNIT=10,FILE='monsster.init.chain',STATUS='old')
      WRITE(6,*) 'reading init chain'
      read(10,*) k
      DO I=1,LENF
         READ(10,*) X(I),Y(I),Z(I)
         xor(i)=x(i)
         yor(i)=y(i)
         zor(i)=z(i)
      ENDDO
      close(10)

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

      ICA(lenf)=0
      ICA(0)=0
c
      write(6,*)
      iflip2=0
      iflip3=0
      qend=0

      DO i=1,NBOX
         DO j=1,NBOX
            DO k=1,NBOX
               xyz(k,j,i)=0
            ENDDO
         ENDDO
      ENDDO

      start=.true.

 8878 DO J=1,LENF
         I=J-1
         II=ICA(I)
         JJ=ICA(J)
         if(j.ne.1.and.j.ne.lenf) THEN
            IF(.NOT.GOODC(II,JJ)) THEN
               WRITE(6,8011)I,J,vx(ii),vy(ii),vz(ii),
     $              vx(jj),vy(jj),vz(jj)
 8011          FORMAT(5X,
     $              'WARNING -WRONG INPUT CHAIN - VECTORS ',8I4)
            ENDIF
         ENDIF
      ENDDO
      
      write(6,4901) itemp, atemp, rgyr, rho
 4901 format(1x,/,1x,'  step', i3,'   Temperature =',f8.4,f8.4,f8.4/)
      write(6,*)
     $     '    #    R2     S2            2_bond   ends  3_bond'
      write(6,*)
c
c
c     ****************** TEMPERATURE RESCALING *******************
c
c     3-body **************************************************
c
      do i=0,19
         do j=0,19
            do k=0,19
               EQ3(k,j,i)=ATEMPO*EQ3(k,j,i)/ATEMP
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
                  envir(ia,im,ip,i)=ATEMPO*envir(ia,im,ip,i)/ATEMP
               end do
            end do
         end do
      end do

         do i=2,lenf1
            ii=seq(i)
            do j=0,100
               if(i.gt.2.and.i.lt.lenf1) then
                  eone(i,j)=enone*eoinp(ii,j)/atemp
               else
                  eone(i,j)=0.0
               endif
            enddo
         enddo

c
c     BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
c

         do i=2,lenf1
            ii=seq(i)
            do j=0, 24
c     j - number of contacts occupied points
               if(j.gt.1.and.j.lt.24) then
                  eburek(i,j)=BUREK*(2.0D0*ECOVER(ii,j)+ECOVER(ii,j-1)
     *                 +ECOVER(ii,j+1))/ATEMP/4.0D0
               else
                  eburek(i,j)=BUREK*ECOVER(ii,j)/ATEMP
               endif
            enddo
         enddo

c

         do i=2,lenf-2
            do k=1,3
               asr2(i,k)=ATEMPO*asr2(i,k)/ATEMP
            enddo
         enddo

         do i=2,lenf-3
            do k=1,4
               asr3(i,k)=ATEMPO*asr3(i,k)/ATEMP
            enddo
         enddo

         do i=2,lenf-4
            do k=1,14
               asr4(i,k)=ATEMPO*asr4(i,k)/ATEMP
            enddo
         enddo

         do i=2,lenf-5
            do k=1,7
               asr5(i,k)=ATEMPO*asr5(i,k)/ATEMP
            enddo
         enddo

c
         do i=2,lenf1
            do j=2,lenf1
               apba(j,i)=ATEMPO*apba(j,i)/atemp
               apbm(j,i)=ATEMPO*apbm(j,i)/atemp
               apbp(j,i)=ATEMPO*apbp(j,i)/atemp
            enddo
         enddo

         ESC= ESCO/ATEMP
         EHBOND=ehbn/ATEMP
         EREST=AREST/ATEMP
         EREP=AAREP/ATEMP
         compress=compr/atemp

         krg=krg/atemp
         krho=krho/atemp

         do i=0,19
            do j=0,19
               do ir=9,13
                  erepa(ir,i,j)=ATEMPO*erepa(ir,i,j)/atemp
                  erepm(ir,i,j)=ATEMPO*erepm(ir,i,j)/atemp
                  erepp(ir,i,j)=ATEMPO*erepp(ir,i,j)/atemp
               enddo
            enddo
         enddo


C     *************************************************************
C     *                                                           *
C     *                 DYNAMICS OF THE CHAIN                     *
C     *                                                           *
C     *************************************************************
c
c
c
         abcd=dble(NCYCLE)

         qlong=0

         eaver=0.0D0
         eaver1=0.0D0         

         orgsum=0.0D0
         nrgsum=0.0D0
         srgsum=0.0D0
         rgsum=0.0D0
         oscont=0.0D0
         nscont=0.0D0
         sscont=0.0D0
         scont=0.0D0

         DO 7777 ICYCLE=1, NCYCLE
c
            DO 7700 IIIDU=1, PHOT

c
c     ******************   INITIAL ENERGY  ************************
c
               do i=0,5
                  ibur(i)=0
               enddo

               energ=ESHORT(2,lenf2)
c
c     .............................................................
c
               SX=0
               SY=0
               SZ=0
               DO I=1,LENF
                  SX=SX+X(I)
                  SY=SY+Y(I)
                  SZ=SZ+Z(I)
               ENDDO

               If(.not. start) then
                  DO I=2,LENF1
                     CALL REM(x(i),y(i),z(i))
                  ENDDO
               endif

               start=.false.

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
c mod meikel
                  NCIJM(i)=0
                  NOA(i)=11
                  NOP(i)=11
               ENDDO
               NBPLUS=0

               ENERG=ENERG+EPAIRN(2,lenf1)
c
c     Correction for the initial 0-contact residues
c
               do i=2,lenf1
                  if(ncij(i).eq.0)
     $                 ENERG=ENERG+EBUREK(i,0)+ENVIR(0,0,0,seq(i))
               enddo
c
               ENERG=ENERG+ECORECT
               orgsum=0.0D0
               oscont=0.0D0
               ECENT=ECENTN(2,lenf1)
               srgsum=nrgsum
               sscont=nscont

               CALL UNMAP

               do i=0,5
                  PLUSB(i)=0
                  MINUSB(i)=0
               enddo

               do i=1,lenf
                  if(NCIJA(i).eq.0) NOA(i)=0
                  if(NCIJP(i).eq.0) NOP(i)=0
               enddo
c
c
c     LONG DISTANCE MOVES (TWO BOND PERMUTATION BASED)
c     (  two-bond permutation alternatively executed )
c

               IF (FXAL16.LE.5) goto 99

               lik=0
 90            k1=FXSTART+INT(rrand(ise)*fxal16)+3
               lik=lik+1
               if(lik.gt.50) go to 99
               k2=k1+INT(rrand(ise)*18.0D0)+3
               if(k2.gt.fxlen-4) GO TO 90
               ii=ica(k1)
               jj=ica(k2)
               ix=vx(ii)
               iy=vy(ii)
               iz=vz(ii)
               k=int(rrand(ise)*47.999D0)+1
               ix=ix+ibux(k)
               iy=iy+ibuy(k)
               iz=iz+ibuz(k)
               ir=ix*ix+iy*iy+iz*iz
               if(ir.lt.9) go to 90
               if(ir.gt.30) go to 90
               jx=vx(jj)
               jy=vy(jj)
               jz=vz(jj)
               jx=jx-ibux(k)
               jy=jy-ibuy(k)
               jz=jz-ibuz(k)
               ir=jx*jx+jy*jy+jz*jz
               if(ir.lt.9) go to 90
               if(ir.gt.30) go to 90
               iiii=vector(ix,iy,iz)
               jjjj=vector(jx,jy,jz)
               if(MOD(lik,2).EQ.0) then
                  iii=jjjj
                  jjj=iiii
               else
                  iii=iiii
                  jjj=jjjj
               endif

               if(ii.eq.iii) go to 90

               if(.NOT.GOODC(ica(k1-1),iii)) go to 90
               if(.NOT.GOODC(ica(k2-1),jjj)) go to 90
               if(.NOT.GOODC(jjj,ica(k2+1))) go to 90
               if(.NOT.GOODC(iii,ica(k1+1))) go to 90

               k3=k2+1
               kp1=k1+1

               do k=kp1,k2
                  call rem(x(k),y(k),z(k))
               enddo

               ica(k1)=iii
               ica(k2)=jjj
               ix=vx(ii)-vx(iii)
               iy=vy(ii)-vy(iii)
               iz=vz(ii)-vz(iii)
               do k=kp1,k2
                  x(k)=x(k)-ix
                  y(k)=y(k)-iy
                  z(k)=z(k)-iz
               enddo

               do k=kp1,k2
                  if(look(x(k),y(k),z(k))) then
                     if(k.gt.kp1) then
                        do i=kp1,k-1
                           call rem(x(i),y(i),z(i))
                        enddo
                     endif
                     ica(k1)=ii
                     ica(k2)=jj
                     do kk=kp1,k2
                        x(kk)=x(kk)+ix
                        y(kk)=y(kk)+iy
                        z(kk)=z(kk)+iz
                     enddo
                     do kk=kp1,k2
                        call set(x(kk),y(kk),z(kk),kk)
                     enddo
                     go to 99
                  endif
                  call set(x(k),y(k),z(k),k)
               enddo
c     return the old conformation


               do k=kp1,k2
                  call rem(x(k),y(k),z(k))
               enddo
               ica(k1)=ii
               ica(k2)=jj
               do kk=kp1,k2
                  x(kk)=x(kk)+ix
                  y(kk)=y(kk)+iy
                  z(kk)=z(kk)+iz
               enddo

               do k=kp1,k2
                  call set(x(k),y(k),z(k),k)
               enddo

               ECO=ECENTO(kp1,k2)
               EOLD=ESHORT(k1,k3)+EPAIRO(k1,k3)


               ica(k1)=iii
               ica(k2)=jjj
               do k=kp1,k2
                  x(k)=x(k)-ix
                  y(k)=y(k)-iy
                  z(k)=z(k)-iz
               enddo

               ECN=ECENTN(kp1,k2)
               ENEW=ESHORT(k1,k3)+EPAIRN(k1,k3)

               de=ENEW-EOLD


               IF(exp(-DE-ECN+ECO). LT. rrand(ise)) then

c     return the old conformation

                  do k=kp1,k2
                     call rem(x(k),y(k),z(k))
                  enddo
                  ica(k1)=ii
                  ica(k2)=jj
                  do kk=kp1,k2
                     x(kk)=x(kk)+ix
                     y(kk)=y(kk)+iy
                     z(kk)=z(kk)+iz
                  enddo
                  do k=kp1,k2
                     call set(x(k),y(k),z(k),k)
                  enddo

                  do i=0,5
                     IBUR(i)=IBUR(i)+MINUSB(i)-PLUSB(i)
                  enddo
                  CALL REMAP
                  GO TO 99
               endif

               qlong=qlong+1
               energ=energ+de
               ecent=ecent+ECN-ECO
               srgsum=nrgsum
               sscont=nscont
               CALL UNMAP


c
c     .........KINKS AND THREE BOND MOVES................

 99            DO 777        IDUMI=1,FXLEN5

c     TWO BOND KINKS

                  j=FXSTART+INT(rrand(ise)*FXAL4)+3
                  i=j-1
                  ii=ica(i)
                  jj=ica(j)
                  ix=vx(ii)
                  iy=vy(ii)
                  iz=vz(ii)
                  k=J+1

c
c     REDEFINED: generate a move and make a biass against an unlike
c     distances
c
 3307             iix=ix+INT(rrand(ise)*3.0D0)-1
                  if(iix.gt.5) iix=5
                  if(iix.lt.-5)iix=-5
                  iiy=iy+INT(rrand(ise)*3.0D0)-1
                  if(iiy.gt.5) iiy=5
                  if(iiy.lt.-5) iiy=-5
                  iiz=iz+INT(rrand(ise)*3.0D0)-1
                  if(iiz.gt.5) iiz=5
                  if(iiz.lt.-5) iiz=-5
                  ir=(ix-iix)**2+(iy-iiy)**2+(iz-iiz)**2
                  if(ir.eq.0) go to 7774
                  irr=iix*iix+iiy*iiy+iiz*iiz
                  if(irr.lt.9) go to 3307
                  if(irr.gt.30) go to 3307
                  jjx=x(k)-x(i)-iix
                  if(iabs(jjx).gt.5) go to 3307
                  jjy=y(k)-y(i)-iiy
                  if(iabs(jjy).gt.5) go to 3307
                  jjz=z(k)-z(i)-iiz
                  if(iabs(jjz).gt.5) go to 3307
                  irr=jjx*jjx+jjy*jjy+jjz*jjz
                  if(irr.lt.9) go to 3307
                  if(irr.gt.30) go to 3307

                  nv1=vector(iix,iiy,iiz)
                  nv2=vector(jjx,jjy,jjz)
                  if(.NOT.GOODC(nv1,nv2)) go to 3307
                  ica1=ica(i-1)
                  if(.NOT.GOODC(ica1,nv1)) go to 3307
                  ica3=ica(k)
                  if(.NOT.GOODC(nv2,ica3)) go to 3307

                  call rem(x(j),y(j),z(j))
                  px=X(i)+vx(nv1)
                  py=Y(i)+vy(nv1)
                  pz=Z(i)+vz(nv1)

                  IF(LOOK(px,py,pz)) then
                     call set(x(j),y(j),z(j),j)
                     go to 7774
                  endif
                  call set(x(j),y(j),z(j),j)

                  ECO=ECENTO(j,j)
                  EOLD=ESHORT(i,k)+EPAIRO(i,k)

                  ICA(I)=nv1
                  ICA(j)=nv2
                  kx=x(j)
                  ky=y(j)
                  kz=z(j)
                  x(j)=px
                  y(j)=py
                  z(j)=pz

                  ECN=ECENTN(j,j)
                  ENEW=ESHORT(i,k)+EPAIRN(i,k)

                  de=ENEW-EOLD

                  IF(exp(-DE-ECN+ECO). LT. rrand(ise)) then
                     call rem(px,py,pz)
                     ICA(i)=ii
                     ICA(j)=jj
                     x(j)=kx
                     y(j)=ky
                     z(j)=kz
                     call set(kx,ky,kz,j)
                     do i=0,5
                        IBUR(i)=IBUR(i)+MINUSB(i)-PLUSB(i)
                     enddo
                     CALL REMAP
                  else
c

                     iflip2=iflip2+1
                     energ=energ+de
                     ecent=ecent+ECN-ECO
                     srgsum=nrgsum
                     sscont=nscont
                     CALL UNMAP

                  endif
c
c     three-bond kinks - REDEFINED


 7774             i=FXSTART+INT(rrand(ise)*FXAL7)+3
                  K=I+3
                  II=I-1
                  ICAI=ICA(I)
                  I1=I+1
                  I2=I+2

                  aaaaa=rrand(ise)
                  if(aaaaa.lt.0.25D0) THEN
                     ICA1=ICA(i2)
                     ICA2=ICA(i1)
                     ICA3=ICAI
                  else
                     ix=vx(icai)
                     iy=vy(icai)
                     iz=vz(icai)

 7884                px=INT(rrand(ise)*3.0D0)-1
                     py=INT(rrand(ise)*3.0D0)-1
                     pz=INT(rrand(ise)*3.0D0)-1
                     iix=ix+px
                     if(iix.gt.5) iix=5
                     if(iix.lt.-5)iix=-5
                     iiy=iy+py
                     if(iiy.gt.5) iiy=5
                     if(iiy.lt.-5) iiy=-5
                     iiz=iz+pz
                     if(iiz.gt.5) iiz=5
                     if(iiz.lt.-5) iiz=-5
                     irr=iix*iix+iiy*iiy+iiz*iiz
                     if(irr.lt.9) go to 7884
                     if(irr.gt.30) go to 7884
                     ir=(ix-iix)**2+(iy-iiy)**2+(iz-iiz)**2
                     if(ir.eq.0) go to 7774

                     ican=ica(i2)
                     jjx=vx(ican)-(iix-ix)
                     jjy=vy(ican)-(iiy-iy)
                     jjz=vz(ican)-(iiz-iz)
                     irr=jjx*jjx+jjy*jjy+jjz*jjz
                     if(irr.lt.9) go to 777
                     if(irr.gt.30) go to 777
                     ica2=ica(i1)
                     ica1=vector(iix,iiy,iiz)
                     ica3=vector(jjx,jjy,jjz)
                  endif

                  IF(.NOT.GOODC(ica1,ica2)) GO TO 777
                  IF(.NOT.GOODC(ica2,ica3)) GO TO 777
                  ICA0=ICA(II)
                  if(.NOT.GOODC(ica0,ica1)) go to 777
                  ICA4=ICA(K)
                  if(.NOT.GOODC(ica3,ica4)) go to 777


                  px=x(i1)
                  py=y(i1)
                  pz=z(i1)
                  kx=x(i2)
                  ky=y(i2)
                  kz=z(i2)

                  call rem(px,py,pz)
                  call rem(kx,ky,kz)

                  X1=X(I)+VX(ICA1)
                  Y1=Y(I)+VY(ICA1)
                  Z1=Z(I)+VZ(ICA1)

                  if(look(x1,y1,z1)) then
                     call set(px,py,pz,i1)
                     call set(kx,ky,kz,i2)
                     go to 777
                  endif
                  call set(x1,y1,z1,i1)

                  X2=X1+VX(ICA2)
                  Y2=Y1+VY(ICA2)
                  Z2=Z1+VZ(ICA2)

                  if(look(x2,y2,z2)) then
                     call rem(x1,y1,z1)
                     call set(px,py,pz,i1)
                     call set(kx,ky,kz,i2)
                     go to 777
                  endif

                  call rem(x1,y1,z1)

                  call set(px,py,pz,i1)
                  call set(kx,ky,kz,i2)

                  ECO=ECENTO(i1,i2)
                  EOLD=ESHORT(i,k)+EPAIRO(i,k)

                  iii=ica(i)
                  jjj=ica(i1)
                  kkk=ica(i2)
                  ICA(I)=ICA1
                  ICA(I1)=ICA2
                  ICA(I2)=ICA3

                  X(I1)=X1
                  Y(I1)=Y1
                  Z(I1)=Z1
                  X(I2)=X2
                  Y(I2)=Y2
                  Z(I2)=Z2

                  ECN=ECENTN(i1,i2)
                  ENEW=ESHORT(i,k)+EPAIRN(i,k)

                  de=ENEW-EOLD
                  IF(exp(-DE-ECN+ECO). LT. rrand(ise)) THEN
                     call rem (x1,y1,z1)
                     call rem (x2,y2,z2)
                     ICA(I)=III
                     ICA(I1)=JJJ
                     ICA(I2)=KKK
                     X(I1)=px
                     Y(I1)=py
                     Z(I1)=pz
                     X(I2)=kx
                     Y(I2)=ky
                     Z(I2)=kz
                     call set(px,py,pz,i1)
                     call set(kx,ky,kz,i2)
                     do i=0,5
                        IBUR(i)=IBUR(i)+MINUSB(i)-PLUSB(i)
                     enddo

                     CALL REMAP
                     go to 7771
                  endif

                  iflip3=iflip3+1
                  energ=energ+de
                  ecent=ecent+ECN-ECO
                  srgsum=nrgsum
                  sscont=nscont
                  CALL UNMAP
                  go to 777
c

 7771             lik=0
 190              k1=FXSTART+INT(rrand(ise)*FXAL7)+3
                  lik=lik+1
                  if(lik.gt.10) go to 777
                  k2=k1+3
                  if(k2.gt.fxlen-4) GO TO 190
                  ii=ica(k1)
                  jj=ica(k2)
                  ix=vx(ii)
                  iy=vy(ii)
                  iz=vz(ii)
                  k=int(rrand(ise)*47.999D0)+1
                  ix=ix+ibux(k)
                  iy=iy+ibuy(k)
                  iz=iz+ibuz(k)
                  ir=ix*ix+iy*iy+iz*iz
                  if(ir.lt.9) go to 190
                  if(ir.gt.30) go to 190
                  jx=vx(jj)
                  jy=vy(jj)
                  jz=vz(jj)
                  jx=jx-ibux(k)
                  jy=jy-ibuy(k)
                  jz=jz-ibuz(k)
                  ir=jx*jx+jy*jy+jz*jz
                  if(ir.lt.9) go to 190
                  if(ir.gt.30) go to 190
                  iiii=vector(ix,iy,iz)
                  jjjj=vector(jx,jy,jz)
                  if(MOD(lik,2).EQ.0) then
                     iii=jjjj
                     jjj=iiii
                  else
                     iii=iiii
                     jjj=jjjj
                  endif

                  if(ii.eq.iii) go to 190

                  if(.NOT.GOODC(ica(k1-1),iii)) go to 190
                  if(.NOT.GOODC(ica(k2-1),jjj)) go to 190
                  if(.NOT.GOODC(jjj,ica(k2+1))) go to 190
                  if(.NOT.GOODC(iii,ica(k1+1))) go to 190

                  k3=k2+1
                  kp1=k1+1

                  do k=kp1,k2
                     call rem(x(k),y(k),z(k))
                  enddo

                  ica(k1)=iii
                  ica(k2)=jjj
                  ix=vx(ii)-vx(iii)
                  iy=vy(ii)-vy(iii)
                  iz=vz(ii)-vz(iii)
                  do k=kp1,k2
                     x(k)=x(k)-ix
                     y(k)=y(k)-iy
                     z(k)=z(k)-iz
                  enddo

                  do k=kp1,k2
                     if(look(x(k),y(k),z(k))) then
                        if(k.gt.kp1) then
                           do i=kp1,k-1
                              call rem(x(i),y(i),z(i))
                           enddo
                        endif
                        ica(k1)=ii
                        ica(k2)=jj
                        do kk=kp1,k2
                           x(kk)=x(kk)+ix
                           y(kk)=y(kk)+iy
                           z(kk)=z(kk)+iz
                        enddo
                        do kk=kp1,k2
                           call set(x(kk),y(kk),z(kk),kk)
                        enddo
                        go to 777
                     endif
                     call set(x(k),y(k),z(k),k)
                  enddo
c     return the old conformation


                  do k=kp1,k2
                     call rem(x(k),y(k),z(k))
                  enddo
                  ica(k1)=ii
                  ica(k2)=jj
                  do kk=kp1,k2
                     x(kk)=x(kk)+ix
                     y(kk)=y(kk)+iy
                     z(kk)=z(kk)+iz
                  enddo

                  do k=kp1,k2
                     call set(x(k),y(k),z(k),k)
                  enddo

                  ECO=ECENTO(kp1,k2)
                  EOLD=ESHORT(k1,k3)+EPAIRO(k1,k3)


                  ica(k1)=iii
                  ica(k2)=jjj
                  do k=kp1,k2
                     x(k)=x(k)-ix
                     y(k)=y(k)-iy
                     z(k)=z(k)-iz
                  enddo

                  ECN=ECENTN(kp1,k2)
                  ENEW=ESHORT(k1,k3)+EPAIRN(k1,k3)

                  de=ENEW-EOLD

                  IF(exp(-DE-ECN+ECO). LT. rrand(ise)) then

c     return the old conformation

                     do k=kp1,k2
                        call rem(x(k),y(k),z(k))
                     enddo
                     ica(k1)=ii
                     ica(k2)=jj
                     do kk=kp1,k2
                        x(kk)=x(kk)+ix
                        y(kk)=y(kk)+iy
                        z(kk)=z(kk)+iz
                     enddo
                     do k=kp1,k2
                        call set(x(k),y(k),z(k),k)
                     enddo

                     do i=0,5
                        IBUR(i)=IBUR(i)+MINUSB(i)-PLUSB(i)
                     enddo
                     CALL REMAP
                     GO TO 777
                  endif
                  iflip3=iflip3+1
                  energ=energ+de
                  ecent=ecent+ECN-ECO
                  srgsum=nrgsum
                  sscont=nscont
                  CALL UNMAP

 777           CONTINUE


C
C     END FLIPS (TWO BONDS REARANGEMENTS)
C
C     N-TERMINUS (TAIL)
C
               IF (.NOT. FXNT) GOTO 79

               JV3=ICA(3)
 60            NV2=INT(rrand(ise)*645.99D0 )+1
               NV1=INT(rrand(ise)*645.99D0 )+1

               if(.NOT.GOODC(nv1,nv2)) go to 60
               if(.NOT.GOODC(nv2,jv3)) go to 60


               kx=x(1)
               ky=y(1)
               kz=z(1)
               px=x(2)
               py=y(2)
               pz=z(2)
               call rem(px,py,pz)
               WX2=VX(NV2)
               WY2=VY(NV2)
               WZ2=VZ(NV2)
               X2=X(3)-WX2
               Y2=Y(3)-WY2
               Z2=Z(3)-WZ2
               ica1=ica(1)
               ica2=ica(2)

               if(look(x2,y2,z2)) then
                  call set(px,py,pz,2)
                  go to 79
               endif

               call set(px,py,pz,2)

               ECO=ECENTO(2,2)
               EOLD=ESHORT(2,3)+EPAIRO(2,3)

               WX1=VX(NV1)
               WY1=VY(NV1)
               WZ1=VZ(NV1)
               X1=X2-WX1
               Y1=Y2-WY1
               Z1=Z2-WZ1
               ICA(1)=nv1
               ICA(2)=nv2

               X(2)=X2
               Y(2)=Y2
               Z(2)=Z2
               X(1)=X1
               Y(1)=Y1
               Z(1)=Z1

               ECN=ECENTN(2,2)
               ENEW=ESHORT(2,3)+EPAIRN(2,3)

               de=ENEW-EOLD
               IF(exp(-DE-ECN+ECO). LT. rrand(ise)) then
                  call rem(x2,y2,z2)
                  ICA(1)=ica1
                  ICA(2)=ica2
                  X(2)=px
                  Y(2)=py
                  Z(2)=pz
                  X(1)=kx
                  Y(1)=ky
                  Z(1)=kz
                  call set(px,py,pz,2)
                  do i=0,5
                     IBUR(i)=IBUR(i)+MINUSB(i)-PLUSB(i)
                  enddo

                  CALL REMAP
                  GO TO 79
               endif

               qend=qend+1
               energ=energ+de
               ecent=ecent+ECN-ECO
               srgsum=nrgsum
               sscont=nscont
               CALL UNMAP


c
C
C     C-TERMINUS (HEAD)
C

 79            IF (.NOT. FXCT) GOTO 7700

               JV3=ICA(LENF3)
 80            NV2=INT(rrand(ise)*645.99D0 )+1
               NV1=INT(rrand(ise)*645.99D0 )+1
               if(.NOT.GOODC(nv2,nv1)) go to 80
               if(.NOT.GOODC(jv3,nv2)) go to 80

               kx=x(lenf)
               ky=y(lenf)
               kz=z(lenf)
               px=x(lenf1)
               py=y(lenf1)
               pz=z(lenf1)
               call rem(px,py,pz)
               WX2=VX(NV2)
               WY2=VY(NV2)
               WZ2=VZ(NV2)
               X2=X(LENF2)+WX2
               Y2=Y(LENF2)+WY2
               Z2=Z(LENF2)+WZ2

               ica1=ica(lenf1)
               ica2=ica(lenf2)

               if(look(x2,y2,z2)) then
                  call set(px,py,pz,lenf1)
                  go to 7700
               endif

               call set(px,py,pz,lenf1)

               ECO=ECENTO(lenf1,lenf1)
               EOLD=ESHORT(lenf2,lenf1)+EPAIRO(lenf2,lenf1)

               WX1=VX(NV1)
               WY1=VY(NV1)
               WZ1=VZ(NV1)
               X1=X2+WX1
               Y1=Y2+WY1
               Z1=Z2+WZ1
               ICA(lenf1)=nv1
               ICA(lenf2)=nv2
               X(LENF1)=X2
               Y(LENF1)=Y2
               Z(LENF1)=Z2
               X(LENF)=X1
               Y(LENF)=Y1
               Z(LENF)=Z1


               ECN=ECENTN(lenf1,lenf1)
               ENEW=ESHORT(lenf2,lenf1)+EPAIRN(lenf2,lenf1)

               de=ENEW-EOLD
               IF(exp(-DE-ECN+ECO). LT. rrand(ise)) then
                  call rem(x2,y2,z2)
                  ICA(lenf1)=ica1
                  ICA(lenf2)=ica2
                  X(LENF1)=px
                  Y(LENF1)=py
                  Z(LENF1)=pz
                  X(LENF)=kx
                  Y(LENF)=ky
                  Z(LENF)=kz
                  do i=0,5
                     IBUR(i)=IBUR(i)+MINUSB(i)-PLUSB(i)
                  enddo

                  CALL REMAP
                  call set(px,py,pz,lenf1)
                  GO TO 7700
               endif

               qend=qend+1
               energ=energ+de
               ecent=ecent+ECN-ECO
               srgsum=nrgsum
               sscont=nscont
               CALL UNMAP

c
c     ................ALL MOVES COMPLETED........................
C

 7700       CONTINUE

            pppp=0
            qqqq=0

            do k=2,lenf1
               do j=1,24
                  if(SPHERE(k,j).lt.0) pppp=k
                  if(SPHEREM(k,j).gt.0) qqqq=k
               enddo
            enddo
            if(pppp.gt.0.or.qqqq.gt.0) then
               do k=2,lenf1
                  write(6,1195) k,aa(seq(k)),(SPHERE(k,j),  j=1,24)
                  write(6,1195) k,aa(seq(k)),(SPHEREM(k,j),  j=1,24)
               enddo
 1195          format(I3,1x,A3,1x,30i3)
               write(6,*) 'error in sphere or spherem.'
               write(6,*) 'pppp=',pppp,' qqqq=',qqqq
               close(6)
               STOP 'error'
            endif



c
c     Additional test - remove when the program is stable
c

            DO J=1,LENF
               I=J-1
               II=ICA(I)
               JJ=ICA(J)
               if(j.ne.1.and.j.ne.lenf) THEN
                  IF(.NOT.GOODC(II,JJ)) THEN
                     WRITE(6,8111)I,J,vx(ii),vy(ii),vz(ii),
     $                    vx(jj),vy(jj),vz(jj)
 8111                FORMAT(5X,
     $                    'WARNING -WRONG INPUT CHAIN - VECTORS ',8I4)
                     STOP 'error'
                  ENDIF
               ENDIF
            ENDDO
c
            aaa=ICYCLE+(ITEMP-1)*NCYCLE

            if (.not. rex) then
               write(9,814)   ICYCLE,aaa,xshrf,yshrf,zshrf,atemp

 814           format(i5,f6.0,3i5,f8.4)

c            do i=2,lenf1
c               kk=seq(i)
c               iica=ica(i-1)*1938+ica(i)*3-1941
c               ix=NINT(car(iica+1)*cgs(kk))
c               iy=NINT(car(iica+2)*cgs(kk))
c               iz=NINT(car(iica+3)*cgs(kk))
c               xt(i)=x(i)+ix
c               yt(i)=y(i)+iy
c               zt(i)=z(i)+iz
c            enddo

c            xt(1)=X(1)
c            yt(1)=y(1)
c            zt(1)=z(1)
c            zt(lenf)=z(lenf)
c            xt(lenf)=x(lenf)
c            yt(lenf)=y(lenf)

            
               if (dorest) then
                  do i=1,lenf
                     xt(i)=x(i)-xshrf
                     yt(i)=y(i)-yshrf
                     zt(i)=z(i)-zshrf
                  end do
                  
                  write(9,713)    (xt(i),yt(i),zt(i),i=1,lenf)
               else
                  write(9,713)    (x(i),y(i),z(i),i=1,lenf)
               endif

 713           format(12I5)

            endif

            R2=(X(LENF)-X(1))**2+(Y(LENF)-Y(1))**2+(Z(LENF)-Z(1))**2
c            asumr2=asumr2+r2
            AS2=0.
C     CENTRE OF GRAVITY COORDINATES
            SX=0
            SY=0
            SZ=0
            do i=1,lenf
               sx=sx+x(i)
               sy=sy+y(i)
               sz=sz+z(i)
            enddo
            ASX=dble(sx)/lenf
            ASY=dble(sy)/lenf
            ASZ=dble(sz)/lenf

            DO I=1,LENF
               BX=(ASX-X(I))**2
               BY=(ASY-Y(I))**2
               BZ=(ASZ-Z(I))**2
               AS2=AS2+BX+BY+BZ
            ENDDO

            AS2=AS2/LENF
c            asums2=asums2+as2
            aaa=ICYCLE+(ITEMP-1)*NCYCLE
            aflip=iflip3/(al5*aaa*PHOT)
            de=iflip2/(al5*aaa*PHOT)
            aend=qend/(2.0D0*aaa*PHOT)
            IF (NREST.GT.0) THEN
               tsrg=sscont/DBLE(NREST)
            ELSE 
               tsrg=0.0D0
            ENDIF
            WRITE(6,8009) ICYCLE,R2,AS2,qlong,de,aend,aflip,energ,
     $           ECENT,sqrt(srgsum/DBLE(lenf2))*1.45D0,tsrg
            eaver=eaver+energ
            eaver1=eaver1+energ +ECENT

            lstene=energ*atemp
            lstrg=sqrt(srgsum/DBLE(lenf2))*1.45D0
            lstrho=tsrg

 8009       format(I6,I6,F8.1,2x,i6,3f8.4,1x,2f8.2,F8.3,F8.3)
C
 7777    CONTINUE



         DO J=2,LENF1
            I=J-1
            II=ICA(I)
            JJ=ICA(J)
            IF(.NOT.GOODC(II,JJ)) THEN
               WRITE(6,801) I,J
 801           FORMAT(5X,'WRONG OUTPUT CHAIN - VECTORS ',2I4)

            ENDIF
         ENDDO

c
c
c     ****************** end of anealing **************************


         if (rex .or. itemp.eq.tcycle) then
            OPEN(UNIT=11,FILE='monsster.final.chain',STATUS='UNKNOWN')
            WRITE(11,8000) LENF
            
            if (dorest) then
               DO I=1,LENF
                  WRITE(11,8000) X(I)-xshrf,Y(I)-yshrf,Z(I)-zshrf
               ENDDO
            else
               DO I=1,LENF
                  WRITE(11,8000) X(I),Y(I),Z(I)
               ENDDO
            endif
            close(11)
            IF (mpflag .eq. 1) THEN
               mpsend=1
            ENDIF
 8000       format(1x,3i5)
         endif

      goto 8877

 8879 continue
      
      lik=0
      do i=1,nbox
         do k=1,nbox
            do j=1,nbox
               if(xyz(j,k,i).ne.0) lik=lik+1
            enddo
         enddo
      enddo

      kk=19*lenf2

      write(6,*)
      write(6,*) 'TEST OF EXCLUDED VOLUME', kk,LIK
      write(6,*)
c
c
c     ******************   FINAL   ENERGY  ************************
c
c

      do i=0,5
         ibur(i)=0
      enddo

      orgsum=0.0D0
      oscont=0.0D0
      ecent=ECENTN(2,lenf1)
      srgsum=nrgsum
      sscont=nscont

      NNNN=0
      NBBB=0
      do i=2,lenf1
         call rem(x(i),y(i),z(i))
         NNNN=NNNN+NHBN(i)
         NBBB=NBBB+NCIJ(I)
         NHBN(i)=0
         NCIJ(i)=0
         NCIJA(i)=0
         NCIJP(i)=0
c mf
         NCIJM(i)=0
         NOA(i)=11
         NOP(i)=11
         do j=1,24
            sphere(i,j)=0
         enddo
      enddo

      NBPLUS=0

      e1=ESHORT(2,lenf2)
      e2=EPAIRN(2,lenf1)

c
c     Correction for the initial 0-contact residues
c
      do i=2,lenf1
         if(ncij(i).eq.0) then 
            e2=e2+EBUREK(i,0)+ENVIR(0,0,0,seq(i))
         endif
      enddo

      e2=e2+ECORECT
c
      mmmm=0
      mbbb=0
      do i=2,lenf1
         mmmm=mmmm+nhbn(i)
         mbbb=mbbb+ncij(i)
      enddo

      energy=e1+e2
      write(6,*)
      write(6,*)' short range', e1
      write(6,*)' long range ', e2,'     contacts',nbbb,mbbb
      write(6,8910) energy,ecent, NNNN, MMMM
 8910 format(1x,/,1x,' FINAL ENERGY =', 2F10.2,'   H-bonds',2i5,/)

c
c
C     DETAILED ANALYSIS OF THE CHAIN STRUCTURE
C
      WRITE(6,8004)
 8004 FORMAT(1X,/,8X,'VECTOR1',/)

      DO I=2,LENF1
         II=ICA(I-1)
         JJ=ICA(I)
         KK=ICA(I+1)

         ICJ=ICONF(II,JJ)

         WX1=VX(II)
         WY1=VY(II)
         WZ1=VZ(II)
         WX2=VX(JJ)
         WY2=VY(JJ)
         WZ2=VZ(JJ)
         sss=1

         if(I.lt.lenf1) then
            WX3=VX(KK)
            WY3=VY(KK)
            WZ3=VZ(KK)

            R3=(wx1+wx2+wx3)**2+(wy1+wy2+wy3)**2+(wz1+wz2+wz3)**2

            px=wy1*wz2-wy2*wz1
            py=wz1*wx2-wz2*wx1
            pz=wx1*wy2-wx2*wy1
            ihand=px*wx3+py*wy3+pz*wz3
            if(r3.lt.64) sss=3

            i2=i-2
            j2=i+2
            r4=(x(i2)-x(j2))**2+(y(i2)-y(j2))**2+(z(i2)-z(j2))**2
            if(i2.lt.1) r4=0

            if(r4.gt.60 .AND. iconf(ii,kk).gt.36) sss=4

            if(R3.gt.0.and.R4.lt.30) then
               if(ihand.gt.0) sss=2
            endif
         endif
         if(j2.gt.lenf) r4=0
         kx=x(i)
         ky=y(i)
         kz=z(i)

         ix=kx-mid
         iy=ky-mid
         iz=kz-mid

         aract2=dble(ix*ix+iy*iy+iz*iz)
         aract=sqrt(aract2)

         ff=aract/acrit
         ika=int(ff/0.333333)
         if(ihand.lt.0) R3=-R3
         WRITE(6,8003)I,WX1,WY1,WZ1,ICJ,R3, R4,aa(seq(i)),
     *        sss,struct(sec(i)), ika, eone(i,ika), NHBN(i), ncij(i)

      ENDDO
C
 8003 FORMAT(1X,I4,1X,3I3,1X,2I4,i5,2x,a3,i4,1x,a5,
     *     2x,i4,f5.1,I4,i4)

      write(6,*)
      write(6,*) ' BBBBBBB      burial status      BBBBBBBB '
      write(6,*)
      do i=2, lenf1
         kk=0
         pp=0
         do j=1,24
            if(sphere(i,j).gt.0) kk=kk+1
         enddo
         write(6,8) i-1, aa(seq(i)), kk, (sphere(i,j),j=1,24)
      enddo
 8    format(I3,1x,A3,i3,1x,30i2)
      write(6,*)

      do i=fxstart+1,fxend+1
         do k=1,24
            if(sphere(i,k).ne.spherem(i,k)) then
               write(6,9) i-1, k, aa(seq(i)), (spherem(i,j),j=1,24)
               write(6,*) 'sphere ne spherem'
               close(6)
               STOP 'error'
            endif
         enddo
      enddo
 9    format(I3,I3,1x,A3,4x,30i2)

c
c     ===================== CONTACT MAP ========================
c
      do i=2,lenf1
         call set(x(i),y(i),z(i),i )
         do j=2,lenf1
            neib(j,i)=0
         enddo
      enddo


      write(6,6003)
 6003 format(1x,/,1x,' *********   CONTACT MAP    **********',/)

      ll=0
      do i=2,lenf1
         NCONT=0
         kx=x(i)
         ky=y(i)
         kz=z(i)
         do 1009  k=1,42
            ix=kx+NBX(k)
            iy=ky+NBY(k)
            iz=kz+NBZ(k)
            iii=XYZ(ix,iy,iz)
            if(iii.gt.0) then
               if(iabs(iii-i).gt.1) then
c     detect each contact only once
                  do kk=0,ncont
                     if(iii.eq.icont(kk)) go to 1009
                  enddo
                  NCONT=NCONT+1
                  ICONT(NCONT)=iii
               endif
            endif
 1009    continue

         kk=0
         if(NCONT.gt.0) THEN
            DO kkk=1,NCONT
               j=ICONT(kkk)
               neib(i,j)=1
               neib(j,i)=1
               jx=x(j)
               jy=y(j)
               jz=z(j)
               ir=(jx-kx)*(jx-kx)+(jy-ky)*(jy-ky)+(jz-kz)*(jz-kz)
               if(ir.lt.icutp(seq(i),seq(j))) then
                  kk=kk+1
                  map(kk)=j
               endif
            ENDDO

            ll=ll+ncont
            write(6,6004)
     *           i, aa(seq(i)),(map(j),aa(seq(map(j))),j=1,kk)
         endif

 6004    format(1x,i3,1x,a3,' >-< ',20(i3,1x,a3))
      enddo
      WRITE(6,*)
      WRITE(6,*) ' TOTAL NUMBER OF CONTACTS |i-j|>1 ...',ll/2
      WRITE(6,*)

      WRITE(6,5027)
 5027 FORMAT(1X,//,1X,38(2H .),/)

      if(nrest.gt.0) then
         write(6,*)
         write(6,*) ' SIDE GROUP - SIDE GROUP RESTRAINTS '
         write(6,*)

         do k=1, NREST
            ii=irest(k)
            jj=jrest(k)
            kk=0
            if(neib(ii,jj).gt.0) kk=1
            iii=(x(ii)-x(jj))**2+(y(ii)-y(jj))**2+(z(ii)-z(jj))**2
            Write(6,4911)k,ii-1,aa(SEQ(ii)),jj-1,aa(SEQ(jj)),iii,kk
         enddo
 4911    format(i10,i4,1x,a3,'  with', i4,1x,a3, I9,i6)
      endif

      WRITE(6,*)

      if(nresta.gt.0) then
         write(6,*)
         write(6,*) ' CaH   -  HN-            RESTRAINTS '
         write(6,*)

         do k=1, NRESTA
            i=iresta(k)
            kx=x(i)
            ky=y(i)
            kz=z(i)

            iica=ica(i-1)*1938+ica(i)*3-1941
c     means CaH
            ix=kx+ICAH(iica+1)
            iy=ky+ICAH(iica+2)
            iz=kz+ICAH(iica+3)

            j=jresta(k)
            iica=ica(j-1)*1938+ica(j)*3-1941
c     means NH
            xj=x(j)+ICNH(iica+1)
            yj=y(j)+ICNH(iica+2)
            zj=z(j)+ICNH(iica+3)

            iii=(ix-xj)**2+(iy-yj)**2+(iz-zj)**2
            Write(6,4911)k,i-1,aa(SEQ(i-1)),j-1,aa(SEQ(j-1)),iii
         enddo
      endif

      WRITE(6,*)
      if(nresth.gt.0) then
         write(6,*)
         write(6,*) ' -NH   -  HN-            RESTRAINTS '
         write(6,*)

         do k=1, NRESTH
            i=iresth(k)
            kx=x(i)
            ky=y(i)
            kz=z(i)
c     means NH
            iica=ica(i-1)*1938+ica(i)*3-1941
            ix=kx+ICNH(iica+1)
            iy=ky+ICNH(iica+2)
            iz=kz+ICNH(iica+3)
            j=jresth(k)
            iica=ica(j-1)*1938+ica(j)*3-1941
c     means NH
            xj=x(j)+ICNH(iica+1)
            yj=y(j)+ICNH(iica+2)
            zj=z(j)+ICNH(iica+3)

            iii=(ix-xj)**2+(iy-yj)**2+(iz-zj)**2
            Write(6,4911)k,i-1,aa(SEQ(i-1)),j-1,aa(SEQ(j-1)),iii
         enddo
      endif

      WRITE(6,*)

      aaa=dble(TCYCLE)*dble(NCYCLE)
      eaver=eaver/aaa
      eaver1=eaver1/aaa
      write(6,*)'  E-internal =',eaver
      write(6,*)'  E-total    =',eaver1

      close(6)
c
c     PRINT various MAPS
c

      IF (.NOT. rex) then
         CLOSE(9)
      ENDIF

      STOP ''
      END
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


C     ******************************


      REAL*8 FUNCTION rrand(ise)
      IMPLICIT REAL*8 (A-H)
      common/randmc/a,b,c,d
      g=(a+b+c)*d
      g=g-aint(g)
      a=b
      b=c
      c=g
      rrand=g
      return
      end

C     **************************************************************


      FUNCTION LOOK(i, j, k)
      IMPLICIT REAL*8 (A-H)
      IMPLICIT INTEGER (I-Z)

      PARAMETER(NBOX=100)
      PARAMETER(NDIM=600)
      COMMON /CHAIN/ ICA(-2:NDIM), x(-2:NDIM), Y(-2:NDIM), Z(-2:NDIM)
      LOGICAL LOOK
      COMMON/WORK/ xyz(NBOX,NBOX,NBOX)

      LOOK=.TRUE.

      if(XYZ(i,j,k).gt.0)     RETURN

      if(XYZ(i+1,j,k).gt.0)   RETURN
      if(XYZ(i-1,j,k).gt.0)   RETURN
      if(XYZ(i,j+1,k).gt.0)   RETURN
      if(XYZ(i,j-1,k).gt.0)   RETURN
      if(XYZ(i,j,k+1).gt.0)   RETURN
      if(XYZ(i,j,k-1).gt.0)   RETURN

      if(XYZ(i+1,j+1,k).gt.0) RETURN
      if(XYZ(i+1,j-1,k).gt.0) RETURN
      if(XYZ(i-1,j+1,k).gt.0) RETURN
      if(XYZ(i-1,j-1,k).gt.0) RETURN

      if(XYZ(i,j+1,k+1).gt.0) RETURN
      if(XYZ(i,j+1,k-1).gt.0) RETURN
      if(XYZ(i,j-1,k+1).gt.0) RETURN
      if(XYZ(i,j-1,k-1).gt.0) RETURN

      if(XYZ(i+1,j,k+1).gt.0) RETURN
      if(XYZ(i+1,j,k-1).gt.0) RETURN
      if(XYZ(i-1,j,k+1).gt.0) RETURN
      if(XYZ(i-1,j,k-1).gt.0) RETURN
c
      LOOK=.FALSE.

      RETURN
      END


C     **************************************************************

      SUBROUTINE SET(i, j, k, m)
      IMPLICIT REAL*8 (A-H)
      IMPLICIT INTEGER (I-Z)

      PARAMETER(NBOX=100)
      PARAMETER(NDIM=600)
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

      PARAMETER(NBOX=100)
      PARAMETER(NDIM=600)
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

      PARAMETER(NDIM=600)

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
      INTEGER FXSTART,FXEND,FXLEN,FXLEN5
      LOGICAL FXNT,FXCT
      COMMON /LENGTHS/ LENF2,LENF1,LENF,LENGTH2(646),
     $     FXSTART,FXEND,FXLEN,FXLEN5,FXNT,FXCT,nrest

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
      if(i1.lt.2) i1=2
      i2=jjjj
      if(i2.gt.lenf2) i2=lenf2


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

      fj=cgs(seq(i1))
      iica=ica(i1-1)*1938+ica(i1)*3-1941
      asx=(dble(x(i1))+CAR(iica+1)*fj)
      asy=(dble(y(i1))+CAR(iica+2)*fj)
      asz=(dble(z(i1))+CAR(iica+3)*fj)

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

         ar=(asx-bsx)**2+(asy-bsy)**2+(asz-bsz)**2
         ESHORT=ESHORT+ESC2*(sqrt(ar)-2.6D0)**2
         if(i.gt.2) then
            ar=sqrt((aax-bsx)**2+(aay-bsy)**2+(aaz-bsz)**2)
            if(ar.gt.5.5D0) ar=5.5D0+(ar-5.5D0)/3.0D0
            ESHORT=ESHORT+ESC2*abs((ar-3.6D0)*(4.5D0-ar))
         endif

         aax=asx
         aay=asy
         aaz=asz
         asx=bsx
         asy=bsy
         asz=bsz

         ib2=ibb2(length2(icai))
         ESHORT=ESHORT+asr2(i,ib2)

         if (i.lt.lenf2) then

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


      if(i2.lt.lenf2) then
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

      if(i2.gt.lenf2-3) i2=lenf2-3
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
      if(i1.lt.2) i1=2
      i2=JJJJ
      if(i2.gt.lenf1-12) i2=lenf1-12

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

      RETURN
      END

c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c

      FUNCTION EPAIRO(ISTART,IEND)
      IMPLICIT REAL*8 (A-H)
      IMPLICIT INTEGER (I-Z)

      PARAMETER(NDIM=600)
      PARAMETER(NBOX=100)
      PARAMETER (IHCB=26)

      common/icutt/IRESA,icutp(0:19,0:19),icuta(0:19,0:19),
     *     icutm(0:19,0:19)
      COMMON /SEQE/  SEQ(ndim), SEC(ndim), ENVIR(0:5,0:5,0:5,0:19)
      COMMON /CHAIN/ ICA(-2:NDIM), x(-2:NDIM), Y(-2:NDIM), Z(-2:NDIM)
      COMMON /pair/ apba(ndim,ndim),apbp(ndim,ndim)
      COMMON /pairm/apbm(ndim,ndim)

      INTEGER FXSTART,FXEND,FXLEN,FXLEN5
      LOGICAL FXNT,FXCT
      COMMON /LENGTHS/ LENF2,LENF1,LENF,LENGTH2(646),
     $     FXSTART,FXEND,FXLEN,FXLEN5,FXNT,FXCT,nrest

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
      COMMON /REP/erepp(9:26,0:19,0:19),erepa(9:26,0:19,0:19),
     *     erepm(9:26,0:19,0:19)
      COMMON /CCCC/ ECORECT, EQ3(0:19,0:19,0:19)
c
      EPAIRO=0.
      NHMIN=0
      NBPLUS=0
      EHBOND2=EHBOND/2.0D0

      DO K=ISTART, IEND
         icont(0)=0
         kx=x(k)
         ky=y(k)
         kz=z(k)
         CALL REM(kx,ky,kz)
c
c     Detect the contacts for the k-th residue
c
         NCONT=0
         do 1000  i=1,42
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
c     !!!!        Hydrogen bonding of the k-th residue ***********
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

               if(ir.lt.IHCB)        then

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
                  if(aisr.lt.9.0D0) EPAIRO=EPAIRO+EREP

c     + add the backbone effect

                  ajsy=bsy-dble(jy)
                  ajsx=bsx-dble(jx)
                  ajsz=bsz-dble(jz)
                  aisr=ajsx*ajsx+ajsy*ajsy+ajsz*ajsz
                  if(aisr.lt.9.0D0) EPAIRO=EPAIRO+EREP

                  amx=asx-bsx
                  amy=asy-bsy
                  amz=asz-bsz
                  ar= amx*amx +amy*amy +amz*amz
                  IF(ar.lt.8.0D0.AND.iabs(k-j).gt.3) EPAIRO=EPAIRO+EREP

                  ii=ica(j-1)
                  jj=ica(j)

c
c     Angular term for tertiary interactions
c
                  mx=vx(ii)-vx(jj)
                  my=vy(ii)-vy(jj)
                  mz=vz(ii)-vz(jj)
                  akx=sqrt(float(mx*mx+my*my+mz*mz)+0.0001)
                  ip = wx*mx + wy*my + wz*mz
                  akx=float(ip)/akx/awx
                  sj=seq(j)
                  sk=seq(k)
                  if(akx.gt.0.5) then

                     if(ir.lt.icutp(sj,sk)) then
                        NCIJ(k)=NCIJ(k)-1
                        NCIJ(j)=NCIJ(j)-1
                        NBPLUS=NBPLUS+1
                        INTO(NBPLUS)=j
                        NCIJD(NBPLUS)=-1
                        NCIJDP(NBPLUS)=-1
                        NCIJDA(NBPLUS)=0
                        NCIJDM(NBPLUS)=0
                        NCIJP(j)=NCIJP(j)-1
                        NBPLUS=NBPLUS+1
                        INTO(NBPLUS)=k
                        NCIJD(NBPLUS)=-1
                        NCIJDP(NBPLUS)=-1
                        NCIJDA(NBPLUS)=0
                        NCIJDM(NBPLUS)=0
                        NCIJP(k)=NCIJP(k)-1

                        kpoints=kpb(k,j)
c     kpb(k,j) coverage of k by j
                        nr=NVECTOR(jjx,jjy,jjz)
                        do kp=1,kpoints
                           kpp=shell(NR,kp)
                           sphere(k,kpp)=sphere(k,kpp)-1
                           spherem(k,kpp)=spherem(k,kpp)-1
                        enddo

                        jpoints=kpb(j,k)
c     kpb(j,k) coverage of j by k
                        nr=NVECTOR(-jjx,-jjy,-jjz)
                        do jp=1,jpoints
                           jpp=shell(NR,jp)
                           sphere(j,jpp)=sphere(j,jpp)-1
                           spherem(j,jpp)=spherem(j,jpp)-1
                        enddo

                        if(erepp(ir,sj,sk).gt.0.01) then
                           EPAIRO=EPAIRO+erepp(ir,sj,sk)
                        else
                           EPAIRO=EPAIRO+apbp(j,k)
                        endif
                     endif
                  else

                     if(akx.lt.-0.5) then

                        if(ir.lt.icuta(sj,sk)) then
                           NCIJ(k)=NCIJ(k)-1
                           NCIJ(j)=NCIJ(j)-1
                           NBPLUS=NBPLUS+1
                           INTO(NBPLUS)=j
                           NCIJD(NBPLUS)=-1
                           NCIJDA(NBPLUS)=-1
                           NCIJDP(NBPLUS)=0
                           NCIJDM(NBPLUS)=0
                           NCIJA(j)=NCIJA(j)-1
                           NBPLUS=NBPLUS+1
                           INTO(NBPLUS)=k
                           NCIJD(NBPLUS)=-1
                           NCIJDA(NBPLUS)=-1
                           NCIJDP(NBPLUS)=0
                           NCIJDM(NBPLUS)=0
                           NCIJA(k)=NCIJA(k)-1

                           kpoints=kpb(k,j)
c     kpb(k,j) coverage of k by j
                           nr=NVECTOR(jjx,jjy,jjz)
                           do kp=1,kpoints
                              kpp=shell(NR,kp)
                              sphere(k,kpp)=sphere(k,kpp)-1
                              spherem(k,kpp)=spherem(k,kpp)-1
                           enddo

                           jpoints=kpb(j,k)
c     kpb(j,k) coverage of j by k
                           nr=NVECTOR(-jjx,-jjy,-jjz)
                           do jp=1,jpoints
                              jpp=shell(NR,jp)
                              sphere(j,jpp)=sphere(j,jpp)-1
                              spherem(j,jpp)=spherem(j,jpp)-1
                           enddo
c
                           if(erepa(ir,sj,sk).gt.0.01) then
                              EPAIRO=EPAIRO+erepa(ir,sj,sk)
                           else
                              EPAIRO=EPAIRO+apba(j,k)
                           endif
                        endif
                     else

                        if(ir.lt.icutm(sj,sk)) then
                           NCIJ(k)=NCIJ(k)-1
                           NCIJ(j)=NCIJ(j)-1
                           NBPLUS=NBPLUS+1
                           INTO(NBPLUS)=j
                           NCIJD(NBPLUS)=-1
                           NCIJDP(NBPLUS)=0
                           NCIJDA(NBPLUS)=0
                           NCIJDM(NBPLUS)=-1
                           NCIJM(j)=NCIJM(j)-1
                           NBPLUS=NBPLUS+1
                           INTO(NBPLUS)=k
                           NCIJD(NBPLUS)=-1
                           NCIJDP(NBPLUS)=0
                           NCIJDA(NBPLUS)=0
                           NCIJDM(NBPLUS)=-1
                           NCIJM(k)=NCIJM(k)-1

                           kpoints=kpb(k,j)
c     kpb(k,j) coverage of k by j
                           nr=NVECTOR(jjx,jjy,jjz)
                           do kp=1,kpoints
                              kpp=shell(NR,kp)
                              sphere(k,kpp)=sphere(k,kpp)-1
                              spherem(k,kpp)=spherem(k,kpp)-1
                           enddo

                           jpoints=kpb(j,k)
c     kpb(j,k) coverage of j by k
                           nr=NVECTOR(-jjx,-jjy,-jjz)
                           do jp=1,jpoints
                              jpp=shell(NR,jp)
                              sphere(j,jpp)=sphere(j,jpp)-1
                              spherem(j,jpp)=spherem(j,jpp)-1
                           enddo

                           if(erepm(ir,sj,sk).gt.0.01) then
                              EPAIRO=EPAIRO+erepm(ir,sj,sk)
                           else
                              EPAIRO=EPAIRO+apbm(j,k)
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
     $                                   (jy-y(i3))**2+(jz-z(i3))**2
                                       if(ijr.lt.icutp(sj,si3)) then
                                          if(iabs(k-j)+iabs(i3-j)+
     $                                         iabs(i3-k).GT.7) then
                                             EPAIRO=EPAIRO+
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

c     parallel (compute "hydrogen bonds") !!!!!!!!!!!!!!!!!!!!

c
c     SECONDARY BIAS  H-bond SELECTIONS     !!!!!!!!!!!!!!!!!!
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
                                 EPAIRO=EPAIRO+(EHBOND+HPROP)
                                 NHBN(k)=NHBN(k)-1
                                 NHMIN=NHMIN+1
                                 IGO(NHMIN)=k
                              endif

                              anx=amx+hjx
                              any=amy+hjy
                              anz=amz+hjz
                              air=anx*anx+any*any+anz*anz
                              if(air.lt.4.0) then
                                 EPAIRO=EPAIRO+(EHBOND+HPROP)
                                 NHBN(j)=NHBN(j)-1
                                 NHMIN=NHMIN+1
                                 IGO(NHMIN)=j
                              endif

                              anx=amx+hx
                              any=amy+hy
                              anz=amz+hz
                              air=anx*anx+any*any+anz*anz
                              if(air.lt.4.0) then
                                 EPAIRO=EPAIRO+(EHBOND+HPROP)
                                 NHBN(k)=NHBN(k)-1
                                 NHMIN=NHMIN+1
                                 IGO(NHMIN)=k
                              endif

                              anx=amx-hjx
                              any=amy-hjy
                              anz=amz-hjz
                              air=anx*anx+any*any+anz*anz
                              if(air.lt.4.0) then
                                 EPAIRO=EPAIRO+(EHBOND+HPROP)
                                 NHBN(j)=NHBN(j)-1
                                 NHMIN=NHMIN+1
                                 IGO(NHMIN)=j
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

      RETURN
      END

c
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
      FUNCTION EPAIRN(ISTART,IEND)
      IMPLICIT REAL*8 (A-H)
      IMPLICIT INTEGER (I-Z)

      PARAMETER(NDIM=600)
      PARAMETER(NBOX=100)
      PARAMETER (IHCB=26)

      common/icutt/IRESA,icutp(0:19,0:19),icuta(0:19,0:19),
     *     icutm(0:19,0:19)
      COMMON /SEQE/  SEQ(ndim), SEC(ndim), ENVIR(0:5,0:5,0:5,0:19)
      COMMON /CHAIN/ ICA(-2:NDIM), x(-2:NDIM), Y(-2:NDIM), Z(-2:NDIM)
      COMMON /pair/ apba(ndim,ndim),apbp(ndim,ndim)
      COMMON /pairm/apbm(ndim,ndim)

      INTEGER FXSTART,FXEND,FXLEN,FXLEN5
      LOGICAL FXNT,FXCT
      COMMON /LENGTHS/ LENF2,LENF1,LENF,LENGTH2(646),
     $     FXSTART,FXEND,FXLEN,FXLEN5,FXNT,FXCT,nrest

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

      DO K=ISTART, IEND
         icont(0)=0
         kx=x(k)
         ky=y(k)
         kz=z(k)

         CALL SET(kx,ky,kz,k)
c
c     Detect the contacts for the k-th residue
c

         NCONT=0
         do 1000  i=1,42
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

      PARAMETER(NDIM=600)
      PARAMETER(NBOX=100)
      PARAMETER (IHCB=26)

      common/icutt/IRESA,icutp(0:19,0:19),icuta(0:19,0:19),
     *     icutm(0:19,0:19)
      COMMON vx(646),vy(646),vz(646)
      COMMON /SEQE/  SEQ(ndim), SEC(ndim), ENVIR(0:5,0:5,0:5,0:19)

      REAL*8 pres,rfac
      INTEGER xor,yor,zor,xshrf,yshrf,zshrf
      LOGICAL dorest
      COMMON /RESTR/ rfac(21,21,21),
     $     pres(NDIM),
     $     xor(NDIM),yor(NDIM),zor(NDIM),
     $     imx(NBOX*2), dorest,
     $     xshrf,yshrf,zshrf

      COMMON/BUR/IBUR(0:5),MINUSB(0:5),PLUSB(0:5),IBTARG(0:5)
      COMMON /RES/ EREP,EREST,APHOBI,MRES(ndim),KRES(ndim,50)
      COMMON /RNN/ MRESH(ndim),KRESH(ndim,4)
      COMMON /CHAIN/ ICA(-2:NDIM), x(-2:NDIM), Y(-2:NDIM), Z(-2:NDIM)
      COMMON/RCN/MRESA(ndim),KRESA(ndim,4),iresta(ndim),jresta(ndim)
      COMMON /CENTER/  eone(ndim,0:100), ACRIT, compress, amid, mid

      INTEGER FXSTART,FXEND,FXLEN,FXLEN5
      LOGICAL FXNT,FXCT
      COMMON /LENGTHS/ LENF2,LENF1,LENF,LENGTH2(646),
     $     FXSTART,FXEND,FXLEN,FXLEN5,FXNT,FXCT,nrest

      COMMON/CA/cgs(0:20)
      COMMON/CAC/CAR(646*646*3)
      COMMON/CAH/ICAH(646*646*3)
      COMMON/CNH/ICNH(646*646*3)

      COMMON/kdh/  eonekd(0:19)

      COMMON/umbr/ rgsum,orgsum,nrgsum,tsrg,srgsum,
     $             rgyr,krg,rho,krho,scont,
     $             oscont,nscont,sscont
      REAL*8 rgsum,orgsum,nrgsum,tsrg,srgsum,
     $     rgyr,krg,rho,krho,scont,oscont,
     $     nscont,sscont

c
c     centrosymmetric one body potential
c

      ECENTO=0.
      do i=0,5
         ECENTO=ECENTO+compress*(MAX(2,IABS(IBUR(i)-IBTARG(i)))-2)
         MINUSB(i)=0
      enddo

c positional restraint, meikel

      if (dorest) then
         aeval=0.0D0
         neval=0
         do i=istart, iend
            if (pres(i).gt.0.000001) then
               kx=x(i)-xor(i)-xshrf+NBOX
               ky=y(i)-yor(i)-yshrf+NBOX
               kz=z(i)-zor(i)-zshrf+NBOX

               eval=rfac(imx(kx),imx(ky),imx(kz))*pres(i)

               aeval=aeval+eval
               if (eval.gt.0.00001D0) then
                  neval=neval+1
               endif
               ECENTO=ECENTO+eval
            endif
         enddo
      endif

c mf, rgyr umbrella
      rgsum=0.0D0
      scont=0.0D0

      do i=istart, iend
         kx=x(i)
         ky=y(i)
         kz=z(i)

         ix=kx-mid
         iy=ky-mid
         iz=kz-mid

         aract2=dble(ix*ix+iy*iy+iz*iz)
         aract=sqrt(aract2)

         rgsum=rgsum+aract2

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

                  IF (KRHO.GT.0.00001) THEN
                     scont=scont+1.0D0/(1.0D0+
     $                    exp(20.0D0*(sqrt(dble(ir))-5.1)))
                  ENDIF
                  
               endif
            enddo
         endif
      enddo

      orgsum=srgsum-rgsum
      tsrg=sqrt(srgsum/DBLE(lenf2))*1.45D0-rgyr
      ECENTO=ECENTO+krg*(tsrg*tsrg)

      IF (nrest.gt.0 .and. krho.gt.0.000001) then
         oscont=sscont-scont
         tsrg=sscont/dble(nrest)-rho
         ECENTO=ECENTO+krho*(tsrg*tsrg)
      ENDIF
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

      PARAMETER(NDIM=600)
      PARAMETER (IHCB=26)
      PARAMETER(NBOX=100)

      common/icutt/IRESA,icutp(0:19,0:19),icuta(0:19,0:19),
     *     icutm(0:19,0:19)
      COMMON vx(646),vy(646),vz(646)
      COMMON /SEQE/  SEQ(ndim), SEC(ndim), ENVIR(0:5,0:5,0:5,0:19)

      REAL*8 pres,rfac
      INTEGER xor,yor,zor,xshrf,yshrf,zshrf
      LOGICAL dorest
      COMMON /RESTR/ rfac(21,21,21),
     $     pres(NDIM),
     $     xor(NDIM),yor(NDIM),zor(NDIM),
     $     imx(NBOX*2), dorest,
     $     xshrf,yshrf,zshrf

      COMMON/BUR/IBUR(0:5),MINUSB(0:5),PLUSB(0:5),IBTARG(0:5)
      COMMON /RES/ EREP,EREST,APHOBI,MRES(ndim),KRES(ndim,50)
      COMMON /RNN/ MRESH(ndim),KRESH(ndim,4)
      COMMON/RCN/MRESA(ndim),KRESA(ndim,4),iresta(ndim),jresta(ndim)
      COMMON /CHAIN/ ICA(-2:NDIM), x(-2:NDIM), Y(-2:NDIM), Z(-2:NDIM)
      COMMON /CENTER/  eone(ndim,0:100), ACRIT, compress, amid, mid

      INTEGER FXSTART,FXEND,FXLEN,FXLEN5
      LOGICAL FXNT,FXCT
      COMMON /LENGTHS/ LENF2,LENF1,LENF,LENGTH2(646),
     $     FXSTART,FXEND,FXLEN,FXLEN5,FXNT,FXCT,nrest

      COMMON/CA/cgs(0:20)
      COMMON/CAC/CAR(646*646*3)
      COMMON/CAH/ICAH(646*646*3)
      COMMON/CNH/ICNH(646*646*3)
      COMMON/kdh/  eonekd(0:19)

      COMMON/umbr/ rgsum,orgsum,nrgsum,tsrg,srgsum,
     $             rgyr,krg,rho,krho,scont,
     $             oscont,nscont,sscont
      REAL*8 rgsum,orgsum,nrgsum,tsrg,srgsum,
     $     rgyr,krg,rho,krho,scont,oscont,
     $     nscont,sscont

c
c     centrosymmetric one body potential
c

      ECENTN=0.
      do i=0,5
         PLUSB(i)=0
      enddo

c     positional restraint, meikel
      if (dorest) then
         aeval=0.0D0
         neval=0
         do i=istart, iend
            if (pres(i).gt.0.000001) then
               kx=x(i)-xor(i)-xshrf+NBOX
               ky=y(i)-yor(i)-yshrf+NBOX
               kz=z(i)-zor(i)-zshrf+NBOX
               
               eval=rfac(imx(kx),imx(ky),imx(kz))*pres(i)
               aeval=aeval+eval
               if (eval.gt.0.00001D0) then
                  neval=neval+1
               endif
               ECENTN=ECENTN+eval
            endif
         enddo
      endif


c mf, rgyr umbrella
      nrgsum=orgsum
      nscont=oscont

      do i=istart, iend
         kx=x(i)
         ky=y(i)
         kz=z(i)

         ix=kx-mid
         iy=ky-mid
         iz=kz-mid

         aract2=dble(ix*ix+iy*iy+iz*iz)

         nrgsum=nrgsum+aract2

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
                  IF (KRHO.GT.0.00001) THEN
                     nscont=nscont+1.0D0/(1.0D0+
     $                    exp(20.0D0*(sqrt(dble(ir))-5.1)))
                  ENDIF
               endif
            enddo
         endif
      enddo

      tsrg=sqrt(nrgsum/dble(lenf2))*1.45D0-rgyr
      ECENTN=ECENTN+krg*(tsrg*tsrg)

      IF (NREST.GT.0 .AND. KRHO.GT.0.000001) THEN
         tsrg=nscont/dble(nrest)-rho
         ECENTN=ECENTN+krho*(tsrg*tsrg)
      ENDIF
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

c
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE REMAP
      IMPLICIT REAL*8 (A-H)
      IMPLICIT INTEGER (I-Z)

      PARAMETER(NDIM=600)

      COMMON /REHAB/NHPLUS,NHMIN,IGO(1600),ITO(1600)
      COMMON /HBN/ EHBOND, NHBN(ndim), prod(646,646)
      COMMON /RENAB/NBPLUS,NCIJ(ndim),NCIJD(1600),INTO(1600)
      COMMON /SHELL2/ SPHERE(NDIM,24), SPHEREM(NDIM,24)

      COMMON/U/ NCIJA(ndim),NCIJP(ndim),NCIJDA(1600),NCIJDM(1600),
     *     NCIJDP(1600) ,NOA(ndim),NOP(ndim), NCIJM(ndim), NOM(ndim)

      if(NHPLUS.GT.0) then
         do i=1,NHPLUS
            k=ITO(i)
            NHBN(k)=NHBN(k)-1
         enddo
      endif

      if(NHMIN.GT.0) then
         do i=1,NHMIN
            k=IGO(i)
            NHBN(k)=NHBN(k)+1
         enddo
      endif


      if(NBPLUS.GT.0) then
         do i=1,NBPLUS
            k=INTO(i)
            NCIJ(k)=NCIJ(k)-NCIJD(i)
            NCIJA(k)=NCIJA(k)-NCIJDA(i)
            NCIJP(k)=NCIJP(k)-NCIJDP(i)

            do j=1,24
               SPHERE(k,j)=SPHERE(k,j)-SPHEREM(k,j)
               IF (SPHERE(k,j).LT.0) THEN
                  WRITE(6,*) 'k=',k,' j=',j,' sphere=',sphere(k,j)
               ENDIF
               SPHEREM(k,j)=0
            enddo
         enddo
      endif

      RETURN
      END


      SUBROUTINE UNMAP
      IMPLICIT REAL*8 (A-H)
      IMPLICIT INTEGER (I-Z)

      PARAMETER(NDIM=600)

      COMMON /REHAB/NHPLUS,NHMIN,IGO(1600),ITO(1600)
      COMMON /HBN/ EHBOND, NHBN(ndim), prod(646,646)
      COMMON /RENAB/NBPLUS,NCIJ(ndim),NCIJD(1600),INTO(1600)
      COMMON /SHELL2/ SPHERE(NDIM,24), SPHEREM(NDIM,24)

      COMMON/U/ NCIJA(ndim),NCIJP(ndim),NCIJDA(1600),NCIJDM(1600),
     *     NCIJDP(1600) ,NOA(ndim),NOP(ndim), NCIJM(ndim), NOM(ndim)

      if(NBPLUS.GT.0) then
         DO i=1,NBPLUS
            k=into(i)
            NOP(k)=NCIJP(k)
            NOA(k)=NCIJA(k)

            do j=1,24
               SPHEREM(k,j)=0
            enddo
         ENDDO
      endif

      RETURN
      END
