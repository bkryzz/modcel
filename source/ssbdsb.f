!  FILE CONTAINS SUBROUTINES TO TEST AGAINST 
!  PRODUCTION OF SSB AND DSB DEFECTS
!  REPARATION OF SSB AND DSB DEFECTS                  .FC 04/2014
!------------------------------------------------------------------
      subroutine ssb_rad(nst,ncel,ps,issb)
!- check for possible production of ssb in cells
      implicit real*8(a-h,o-z)
      dimension ps(0:nst,0:nst),issb(ncel)
      do icel=1,ncel
      isl0=0
      cc=rand()
      do j=0,nst
      if (cc.le.ps(issb(icel),j)) isl0=j
      end do
      if(isl0.gt.issb(icel)) issb(icel)=isl0
      end do
      return
      end
!------------------------------------------------------------------
      subroutine dsb_rad(ndt,ncel,pd,idsb)
!- check for possible production of dsb in cells
      implicit real*8(a-h,o-z)
      dimension pd(0:ndt,0:ndt),idsb(ncel)
      do icel=1,ncel
      idl0=0
      cc=rand()
      do j=0,ndt
      if (cc.le.pd(idsb(icel),j)) idl0=j
      end do
      if(idl0.gt.idsb(icel)) idsb(icel)=idl0
      end do
      return
      end
!------------------------------------------------------------------
      subroutine ssb_rep(nst,ncel,rps,issb)
!- check for possible repair of ssb in damaged cells
      implicit real*8(a-h,o-z)
      dimension rps(0:nst,0:nst),issb(ncel)
      do icel=1,ncel
      ind=issb(icel)
      isl0=0
      cc=rand()
      do j=ind,0,-1
      if (cc.le.rps(ind,j)) isl0=j
      end do
      if(isl0.lt.issb(icel)) issb(icel)=isl0
      end do
      return
      end
!------------------------------------------------------------------
      subroutine dsb_rep(ndt,ncel,rpd,idsb)
!- check for possible repair of dsb in damaged cells
      implicit real*8(a-h,o-z)
      dimension rpd(0:ndt,0:ndt),idsb(ncel)
      do icel=1,ncel
      ind=idsb(icel)
      idl0=0
      cc=rand()
      do j=ind,0,-1
      if (cc.le.rpd(ind,j)) idl0=j
      end do
      if(idl0.lt.idsb(icel)) idsb(icel)=idl0
      end do
      return
      end
