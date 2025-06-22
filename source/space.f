!  FILE CONTAINS ROUTINES FOR SETTING UP SQUARE LATTICE AND 
!  DEFINING NEIGHBOR TABLES 
!  DIAMOND SYMMETRY ONLY AVAILABLE IN CURRENT VERSION
!
!                           .
!                           .
!                          NB4
!                ..... NB3 IC NB1 ....
!                          NB2
!                           .
!                           .
!
!  ARRAYS IPSX,IPSY RANGING FROM -N/2 to +N/2 DESCRIBE DISCRETE
!  SPATIAL COORDINATES
!  ARRAY IOCC NUMBERS SITES FROM 1 to N^2 AND TELLS WHICH CELL IC
!  OCCUPIES A GIVEN LATTICE SITE AT A GIVEN TIME
!  ARRAY ILOC NUMBERS CELLS FROM 1 to NCEL AND TELLS WHICH LATTICE
!  SITE A CELL OCCUPIES (NCEL CAN CHANGE OVER TIME BUT <= N^2)
!  ATTENTION: ARRAYS ILOC AND IOCC NEED TO BE CONSTANTLY UPDATED
!  ALWAYS AT THE SAME TIME, OTHERWISE INFORMATION IS LOST
!
!                          *-----*
!                          |*---*|
!                          ||C-*||     sites are filled with cells
!                          ||  |||     according to spiral arrang
!                          |*--*||     from center
!                          *----*|
!                                *
!------------------------------------------------------------------
      subroutine space(numc,ncel,ipsx,ipsy,ipsz,iloc,iocc,innb)
      implicit real*8 (a-h,o-z)
!- new subroutine space: assign space location to cells on a fixed
!- lattice. each lattice site has a (x,y,z) coordinate in (-N/2,+N/2)
!- sites are labelled from 1 to (N+1)**2 (incl. 0)
      dimension ipsx(numc),ipsy(numc),ipsz(numc),iloc(numc),innb(numc,6)
      dimension iocc(numc)
      n=int(dsqrt(dfloat(numc)))
      if (n*n.ne.numc) stop 888
      nm=n/2
!     il=-n2/2+n/2
      il=0
      do k=0,0
        do j=-nm,nm 
          do i=-nm,nm
           il=il+1
           ipsx(il)=i
           ipsy(il)=j
           ipsz(il)=k
          end do
        end do
      end do

!- create complete nabor tables
      call nnbtab(numc,ipsx,ipsy,ipsz,innb)

!- now fill sites in RH spiral progression starting from center
      inc=(n*n+1)/2    !   is central site (ix,iy,iz)=(0,0,0)
       m=1
      iloc(m)=inc
      iocc(inc)=m
      do k=1,numc
      if (2*(k/2).eq.k) then   ! select odd and even branch of the spiral
       do k1=1,k
       m=m+1                   ! increase cell counter
       iloc(m)=innb(inc,1)     ! location of cell m
       inc=iloc(m)             
       iocc(inc)=m             ! lattice site inc is occupied by cell m
       if (m.eq.ncel) go to 99
       end do
       do k1=1,k
       m=m+1
       iloc(m)=innb(inc,4)
       inc=iloc(m)
       iocc(inc)=m  
       if (m.eq.ncel) go to 99
       end do
      else
       do k1=1,k
       m=m+1
       iloc(m)=innb(inc,2)
       inc=iloc(m)
       iocc(inc)=m  
       if (m.eq.ncel) go to 99
       end do
       do k1=1,k
       m=m+1
       iloc(m)=innb(inc,3)
       inc=iloc(m)
       iocc(inc)=m  
       if (m.eq.ncel) go to 99
       end do
      end if
      end do

  99  continue
      return
      end
!------------------------------------------------------------------
      subroutine nnbtab(numc,ipsx,ipsy,ipsz,innb)
      implicit real*8 (a-h,o-z)
!- new subroutine nnbtab: build nearest neighbors tables
      dimension ipsx(numc),ipsy(numc),ipsz(numc),innb(numc,6)
      n=int(dsqrt(dfloat(numc)))
      do k=1,numc
      innb(k,1)=k+1 
         if (abs(ipsx(k)-ipsx(k+1)).gt.1)  innb(k,1)=0
      innb(k,2)=k-1  
         if (abs(ipsx(k)-ipsx(k-1)).gt.1)  innb(k,2)=0
      innb(k,3)=k+n 
         if (abs(ipsy(k)-ipsy(k+n)).gt.1)  innb(k,3)=0
      innb(k,4)=k-n
         if (abs(ipsy(k)-ipsy(k-n)).gt.1)  innb(k,4)=0
!     print *, ipsx(k),ipsy(k),k,(innb(k,l),l=1,6)
      end do 
      return 
      end
