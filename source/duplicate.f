!  FILE CONTAINS ROUTINE FOR CELL DUPLICATION 
!  AT EACH DUPLICATION A NEW CELL TAKES THE SITE OF THE MOTHER
!  MOTHER IS DISPLACED TO 1...4 ACCORDING TO INDEX RANDOM VAR.
!  ALL NEIGHBOR CELLS ARE SHIFTED IN THE SAME DIRECTION
!------------------------------------------------------------------
      subroutine double(duprb,ntime_cy,itime,ncmax,numc,iphse,
     &                  issb,idsb,iloc,iocc,ityp,innb,idupl,indupl)
!- test for cell duplication
      implicit none
      integer*4 :: ntime_cy, itime, ncmax, numc, indupl, iss
      integer*4 :: k0, k1, k2, k3, k4, k5, kc, lc, index, isafe
      integer*4 :: issb(ncmax),idsb(ncmax),iloc(ncmax),ityp(ncmax)
      integer*4 :: iocc(ncmax),innb(ncmax,6),idupl(ncmax),iphse(ncmax)
      real*8 :: dupl_time, prob, cc, duprb
      integer*4, allocatable, save :: ndupl(:)
      
      if (.not.allocated(ndupl)) then
          allocate (ndupl(ncmax)) ; ndupl=0 ; endif
      
      do lc=1,numc
      iss=iloc(lc)
c     print *, lc,iocc(iss)
      end do

      indupl=0
      do k0=1,numc
      if (iphse(k0).eq.0) go to 99    !  cell is quiescent G0 or dead
      kc=k0
      dupl_time=dfloat(itime-idupl(kc))/(duprb*ntime_cy)
      isafe=iphse(kc)/4
      prob=dfloat(isafe)*(1.d0-dexp(-dupl_time))
c     write(66,*) itime,kc,ndupl(kc)
      cc=rand()
      if (cc.lt.prob) then
          indupl=indupl+1
          ndupl(kc)=ndupl(kc)+1
          idupl(kc)=itime
          iphse(kc)=1           !  set cell cycle to G1
          numc=numc+1
          if (numc.gt.ncmax) stop 990
          issb(numc)=issb(kc)
          idsb(numc)=idsb(kc)
          ityp(numc)=ityp(kc)
          iloc(numc)=iloc(kc)
          iocc(iloc(kc))=numc
          idupl(numc)=itime
          iphse(numc)=1

  77      continue
          index=1+4*rand()
          k1=iloc(kc)        ! lattice site of cell kc
          k2=innb(k1,index)  ! nabor of lattice site k1
          k3=iocc(k2)        ! cell occupying nabor site
          iloc(kc)=k2        ! new site of cell kc
          iocc(k2)=kc        ! kc now occupies site k2
          kc=k3              ! set k3 as the new moved cell
          if (kc.ne.0) go to 77
c         call smappa(itime,ncmax,iloc,iocc)  ! plot only for small lattices!!
      
       end if          !  close if-loop on duplication prob for cell kc

  99  continue
      end do           !  close do-loop on cell kc

      return
      end
