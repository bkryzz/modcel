      subroutine brid3d(numc,numf,maxn,innb,nnab,x,y,z,rdi,wid,hei,siz,isphcub,rad,rcut)
!
!  generates a random densely packed distribution of lattice
!  points according to the proximity algorithm of Bridson
!
      implicit none
      integer*4  numc,numf,maxn,nadv,isphcub,nw,nh,ns,nm,ii,nclist
      integer*4  i,j,k,l,k1,k2,m,n,mn,m0,nbb
      integer*4  iloc(numc),innb(maxn,numc),nnab(numc)
      real*8     wid,hei,siz,rad,rad2,rcut,rcut2,r2cut,pigr,duepi,sqrt2
      real*8     widh,heih,sizh,dis,dx,dy,dz,fi,r2,rcyl,rr,w,xi,xxx,yyy,zzz
      real*8     x(numc),y(numc),z(numc),rdi(numc)
      character  simu(3)*6
      data simu/'SPHERE','3-CUBE','CYLNDR'/

!   wid,hei,siz=box dimensions along x,y,z
!   rad=radius of cell ; isphcub=0 no cut, =1 sphere, =2 cube, 3=cylindric
!   ncel=number of active cells
      pigr=4.d0*datan(1.d0)
      duepi=2.d0*pigr
      nw=nint(wid/rad)
      nh=nint(hei/rad)
      ns=nint(siz/rad)
      wid=nw*rad
      hei=nh*rad
      siz=ns*rad
      r2cut=0.385d0*(wid*hei*siz)**(2.d0/3.d0)
      rad2=rad*rad
      widh=0.5d0*wid
      heih=0.5d0*hei
      sizh=0.5d0*siz
      write(*,*) ' BUILDING LINKED-CELL LIST FOR NEIGHBOR SEARCH'
      nclist=nw*nh*ns
      do ii=1,nclist
      
      end do

      write(*,*) ' GENERATING RANDOM PACKED DISTRIBUTION W BRIDSON'
      write(*,'(1x,a,i3,2x,a)') ' SIMULATION BOX =',isphcub,simu(isphcub)
      write(16,*) numc,2,0,0
      write(16,*)
!- place the first cell
      x(1)=0.d0
      y(1)=0.d0
      z(1)=0.5d0*rad
      rdi(1)=x(1)*x(1)+y(1)*y(1)+z(1)*z(1)

!- start picking cells according to Bridson algorithm
      do k=2,numc
      m0=0
  99  continue
      rr=rad*(1.+rand())
      fi=duepi*rand()
      xi=pigr*(-1.d0+2.d0*rand())
      ii=rand()*k
      if (ii.eq.0) ii=1
      x(k)=x(ii)+rr*dcos(fi)*dsin(xi)
      y(k)=y(ii)+rr*dsin(fi)*dsin(xi)
      z(k)=z(ii)+rr*dcos(xi)
      r2=x(k)*x(k)+y(k)*y(k)+z(k)*z(k)
      if (isphcub.eq.1) then      ! apply sphere cutoff
         if (r2.gt.r2cut) go to 99
      else if (isphcub.eq.2) then ! apply cube cutoff
         if (x(k).ge.widh.or.x(k).le.-widh) go to 99
         if (y(k).ge.heih.or.y(k).le.-heih) go to 99
         if (z(k).ge.sizh.or.z(k).le.-sizh) go to 99
      else if (isphcub.eq.3) then ! apply cylindric cutoff
         rcyl=x(k)**2+y(k)**2
      end if
      m0=m0+1
        do l=1,k-1
        dx=x(k)-x(l)
        dy=y(k)-y(l)
        dz=z(k)-z(l)
        dis=dx*dx+dy*dy+dz*dz
        m0=m0+1
        if (dis.lt.rad2) go to 99
        end do
      rdi(k)=r2
      nadv=mod(k,100)
      if (nadv.eq.0) write(*,'(a,i7)',advance='no') '...',k
!     if (nadv.eq.0) write(*,'(2i7)') k,m0
      end do

!- now sort cells according to distance from the center
      write(*,'(/,a)') '  NOW SORTING CELL SITES CLOSER TO BOX CENTER'
      do i=1,numc
!     print *, i,rdi(i),x(i),y(i),z(i)
      end do
      do j=2,numc       ! Pick out each element in turn.
      w=rdi(j)
      xxx=x(j)
      yyy=y(j)
      zzz=z(j)
      do i=j-1,1,-1     ! Look for the place to insert it.
      if (rdi(i).le.w) go to 10
          rdi(i+1)=rdi(i)
          x(i+1)=x(i)
          y(i+1)=y(i)
          z(i+1)=z(i)
      end do
      i=0
  10  rdi(i+1)=w        ! Insert it.
      x(i+1)=xxx
      y(i+1)=yyy
      z(i+1)=zzz
      enddo
      do i=1,numc
!     print *, i,rdi(i),x(i),y(i),z(i)
      write(16,*) 'C  ',x(i),y(i),z(i)
      end do

!- now call voronoi routines
      numf=0.6d0*numc
      call voron3(32,numc,numf,maxn,innb,nnab,x,y,z,wid,rcut)
!     do nn=1,numf
!     print *, nn,nnab(nn)
!     print *, (innb(mm,nn),mm=1,nnab(nn))
!     end do
 
      return
      end
!--------------------------------------------------------------------------
      subroutine fcc3d(numc,numf,maxn,innb,nnab,x,y,z,rdi,wid,hei,siz,isphcub,rad,rcut)
!
!  generates a FCC close-packed distribution of lattice points
!
      implicit none
      integer*4  numc,numf,maxn,nadv,isphcub,nw,nh,ns,nm,ii,nclist
      integer*4  i,j,k,l,k1,k2,m,n,mn,m0,nbb
      integer*4  iloc(numc),innb(maxn,numc),nnab(numc)
      real*8     wid,hei,siz,rad,rad2,rcut,rcut2,r2cut,pigr,duepi,sqrt2
      real*8     x(numc),y(numc),z(numc)
      real*8     w,dd,dx,dy,dz,xxx,yyy,zzz,rdi(numc)
      character  simu(2)*6
      parameter (nm=20)
      data simu/'SPHERE','3-CUBE'/

!   wid,hei,siz=box dimensions along x,y,z
!   rad=radius of cell ; isphcub=0 no cut, =1 sphere, =2 cube
!   ncel=number of active cells
      pigr=4.d0*datan(1.d0)
      duepi=2.d0*pigr
      sqrt2=dsqrt(2.d0)
      nw=nint(wid*sqrt2/rad)
      nh=nint(hei*sqrt2/rad)
      ns=nint(siz*sqrt2/rad)
      nclist=nw*nh*ns
      if (nh.gt.nw) then;nw=nh;wid=hei;end if
      if (ns.gt.nw) then;nw=ns;wid=siz;end if
      rcut=2.*rad
      rcut2=rcut/sqrt2
      r2cut=0.5d0*wid*wid
      rad2=rad*rad
      write(*,*) ' BUILDING LINKED-CELL LIST FOR NEIGHBOR SEARCH'
      do ii=1,nclist
      end do

      write(*,*) ' GENERATING FCC CLOSE-PACKED DISTRIBUTION '
      write(*,'(1x,a,i3,2x,a)') ' SIMULATION BOX =',isphcub,simu(isphcub)
      write(16,*) numc,2,0,0
      write(16,*)
      write(*,*) ' INTEGER CELL SITES ',nw,'**3'

!- generate cells on a close-packed lattice
      n=0
      do k=-nw,nw
      do l=-nw,nw
      do m=-nw,nw
      xxx=(0.0d0+dfloat(k))*rcut
      yyy=(0.0d0+dfloat(l))*rcut
      zzz=(0.0d0+dfloat(m))*rcut
      dd=xxx*xxx+yyy*yyy+zzz*zzz
      if (dd.le.r2cut) then
         n=n+1
         x(n)=xxx;y(n)=yyy;z(n)=zzz
         rdi(n)=dd
      end if 
      xxx=(0.5d0+dfloat(k))*rcut
      yyy=(0.0d0+dfloat(l))*rcut
      zzz=(0.5d0+dfloat(m))*rcut
      dd=xxx*xxx+yyy*yyy+zzz*zzz
      if (dd.le.r2cut) then
         n=n+1
         x(n)=xxx;y(n)=yyy;z(n)=zzz
         rdi(n)=dd
      end if 
      xxx=(0.5d0+dfloat(k))*rcut
      yyy=(0.5d0+dfloat(l))*rcut
      zzz=(0.0d0+dfloat(m))*rcut
      dd=xxx*xxx+yyy*yyy+zzz*zzz
      if (dd.le.r2cut) then
         n=n+1
         x(n)=xxx;y(n)=yyy;z(n)=zzz
         rdi(n)=dd
      end if 
      xxx=(0.0d0+dfloat(k))*rcut
      yyy=(0.5d0+dfloat(l))*rcut
      zzz=(0.5d0+dfloat(m))*rcut
      dd=xxx*xxx+yyy*yyy+zzz*zzz
      if (dd.le.r2cut) then
         n=n+1
         x(n)=xxx;y(n)=yyy;z(n)=zzz
         rdi(n)=dd
      end if 
      end do;end do;end do
      numf=n
!- now sort cells according to distance from the center
      write(*,'(/,a)') '  NOW SORTING CELL SITES CLOSER TO BOX CENTER'

      do j=2,numf       ! Pick out each element in turn.
      w=rdi(j)
      xxx=x(j)
      yyy=y(j)
      zzz=z(j)
      do i=j-1,1,-1     ! Look for the place to insert it.
      if (rdi(i).le.w) go to 10
          rdi(i+1)=rdi(i)
          x(i+1)=x(i)
          y(i+1)=y(i)
          z(i+1)=z(i)
      end do
      i=0
  10  rdi(i+1)=w        ! Insert it.
      x(i+1)=xxx
      y(i+1)=yyy
      z(i+1)=zzz
      enddo
      do i=1,numf
!     print *, i,rdi(i),x(i),y(i),z(i)
!     write(16,*) 'C  ',x(i),y(i),z(i)
      end do
 
      write(*,*) ' CELL CONFIGURATION ',numf
      write(*,*) ' BOX LENGTH =',wid
      write(*,*) ' NEIGHBOUR CUTOFF =',rcut,rcut2,rad

      rcut2=rcut2*rcut2+1.d-3
      do k1=1,numf
      nbb=0
      do k2=1,numf
      if (k2.ne.k1) then
      dx=x(k1)-x(k2)
      dy=y(k1)-y(k2)
      dz=z(k1)-z(k2)
      dd=dx*dx+dy*dy+dz*dz
      if (dd.le.rcut2) then
        nbb=nbb+1
        nnab(k1)=nbb
        innb(nbb,k1)=k2
!       print *,k1,k2,dsqrt(dd)
      end if 
      end if
      end do
      if(nnab(k1).eq.12) then 
      write(16,*) 'C  ',x(k1),y(k1),z(k1)
      else
      write(16,*) 'B  ',x(k1),y(k1),z(k1)
      end if
      end do

      return
      end

