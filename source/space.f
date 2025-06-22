      subroutine space(numc,numf,maxn,ncel,ipsx,ipsy,ipsz,iloc,iocc, &
                       nnab,innb,isphcub,rad,wid,hei,siz,rcut,x,y,z,rdi)

!--------------------------------------------------------------------
!- new subroutine space: generate space location of cells on a fixed
!- lattice of randomly packed balls produced by the Bridson algorithm.
!- basic space-limiting box can be a sphere or a parallelepiped
!
!- CAREFUL: only the first NUMF out of NUMC are used because of
!- low-symmetry box cannot be completely filled by Voronoi polyh
!
!- each lattice site has a (x,y,z) coordinate in (-N/2,+N/2)
!- sites are labelled from 1 to NUMF
!  array IOCC numbers sites from 1 to NUMC and tells which cell
!  occupies a given site at a given time
!  array ILOC numbers cells from 1 to NCEL and tells which site
!  a cell occupies (NCEL can change over time but <= NUMF)
!  CAREFUL: arrays ILOC and IOCC need to be constantly updated
!  always at the same time, otherwise information is lost

!- in the initial configuration the first NCEL<=NUMF cells occupy the
!- sites from 1 to NCEL
!--------------------------------------------------------------------

      implicit none
      integer*4  numc,numf,maxn,ncel,k,isphcub
      integer*4  ipsx(numc),ipsy(numc),ipsz(numc),iloc(numc)
      integer*4  iocc(numc),innb(maxn,numc),nnab(numc)
      real*8     x(numc),y(numc),z(numc),rdi(numc)
      real*8     wid,hei,siz,rad,rcut

!- create simulation box: random-packed distribution with Bridson algorithm
!- sites are ordered w/ respect to closeness to (0,0,0) both for cube and sphere
!- neighbor list is built by 3D Voronoi polygon contact

      call brid3d(numc,numf,maxn,innb,nnab,x,y,z,rdi,wid,hei,siz,isphcub,rad,rcut)
!     call fcc3d (numc,numf,maxn,innb,nnab,x,y,z,rdi,wid,hei,siz,isphcub,rad,rcut)

!- now fill sites in spiral progression starting from center
!- inc=1 is central site (x,y,z)=(0,0,z0)
!- in the initial configuration IOCC=ILOC necessarily

      write(*,*) ncel,' STARTING CELLS PLACED AT THE CENTRE '
      do k=1,ncel
      iocc(k)=k
      iloc(k)=k
      end do
  
      return
      end

!--------------------------------------------------------------------------
      subroutine sqube(idt,numc,numf,maxn,ncel,ipsx,ipsy,ipsz,iloc,iocc, &
                       nnab,innb,isphcub,rad,wid,hei,siz,rcut,x,y,z,rdi)

!
!  generate a square/cubic lattice of points in 2D/3D
!
!- each lattice site has a (x,y,/z) coordinate in (-N/2,+N/2)
!- sites are labelled from 1 to NUMF
!  array IOCC numbers sites from 1 to NUMC and tells which cell
!  occupies a given site at a given time
!  array ILOC numbers cells from 1 to NCEL and tells which site
!  a cell occupies (NCEL can change over time but <= NUMF)
!  CAREFUL: arrays ILOC and IOCC need to be constantly updated
!  always at the same time, otherwise information is lost

!- in the initial configuration the first NCEL<=NUMF cells occupy the
!- sites from 1 to NCEL (assumed confluent)

      implicit none
      integer*4  numc,numf,maxn,nadv,ncel,isphcub,nw,nh,ns,nm,ii,nc
      integer*4  idt,i,j,k,l,k1,k2,m,n,mn,m0,nbb
      integer*4  ipsx(numc),ipsy(numc),ipsz(numc)
      integer*4  iloc(numc),iocc(numc),innb(maxn,numc),nnab(numc)
      real*8     wid,hei,siz,rad,rad2,rcut,rcut2,r2cut,pigr,duepi,sqrt2
      real*8     x(numc),y(numc),z(numc)
      real*8     w,dd,dx,dy,dz,xxx,yyy,zzz,rdi(numc)
      character  simu(2)*6
      parameter (nm=20)
      data simu/'SQUARE','3-CUBE'/

!   wid,hei,siz=box dimensions along x,y,z
!   rad=radius of cell ; isphcub=0 no cut, =1 sphere, =2 cube
!   ncel=number of active cells
      pigr=4.d0*datan(1.d0)
      duepi=2.d0*pigr
      sqrt2=dsqrt(2.d0)
      rad2=2.d0*rad
      nw=nint(wid/rad2)
      if (nw.eq.2*(nw/2)) nw=nw+1
      nh=nint(hei/rad2)
      if (nh.eq.2*(nh/2)) nh=nh+1
      ns=nint(siz/rad2)
      if (ns.eq.2*(ns/2)) ns=ns+1
      numf=nw*nh
      if (idt.eq.3) numf=numf*ns
      if (NCEL .GT. NUMF) then
          write(6,*) 'INPUT ERROR: NCEL .GT. NUMBER OF SITES'
          stop 77
      end if
      
!- start creating the square/cube lattice
      nc=0
      if (idt.eq.2) then
        do i=1,nw
        do j=1,nh
        nc=nc+1
        x(nc)=i*rad2 
        y(nc)=j*rad2 
        end do;end do      
        print *, nw/2,nh/2,ns/2
        write(6,*) x(nw*(nh/2)+nw/2+1),y(nw*(nh/2)+nw/2+1)
      elseif (idt.eq.3) then
        do i=1,nw
        do j=1,nh
        do k=1,ns
        nc=nc+1
        x(nc)=i*rad2 
        y(nc)=j*rad2 
        z(nc)=k*rad2 
        end do;end do;end do      
        print *, nw/2,nh/2,ns/2
      end if
      if (NC .ne. NUMF) print *, 'error nc'

!   CLOSING BLOCK FROM SPACE SUBR NEEDED HERE
      write(*,*) ncel,' STARTING CELLS PLACED AT THE CENTRE '
      do k=1,ncel
      iocc(k)=k
      iloc(k)=k
      end do

      return
      end



