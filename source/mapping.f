!  FILE CONTAINS ROUTINES TO PRODUCE READABLE MAPS OF THE CELL
!  STRUCTURE
!  ROUTINE MAP PRODUCES A PDB FORMAT WITH CELLS LOCATED AT THEIR
!  LATTICE SITE COORDINATE IPSX,IPSY AND ATTACH A CODE GIVING FOR
!  EACH CELL THE NUMBER OF SSB/DSB IN COLUMNS 'BETA' AND 'OCC' OF
!  THE PDB FORMAT (THESE CAN BE PLOTTED, ALSO DYNAMICALLY WITH VMD)
!  ROUTINE SMAPPA PRODUCES A WIDE X-Y FORMAT ASCII TABLE IN WHICH
!  CELL NUMBER IS WRITTEN AT EACH OCCUPIED SITE (0 ELSEWHERE). USEFUL
!  ONLY FOR DEBUGGING PURPOSES.
!------------------------------------------------------------------
      subroutine map(ncmax,numc,issb,idsb,iloc,ipsx,ipsy,ipsz)
      implicit none
!- map cell damage in pdb format (OCC is SSB, BETA is DSB)
      integer*4 numc,ncmax,n,j
      integer*4 issb(ncmax),idsb(ncmax),iloc(ncmax)
      integer*4 ipsx(ncmax),ipsy(ncmax),ipsz(ncmax)

      write(20,100) 'REMARK'
      write(21,*) numc
      write(21,*) 'tit'
      do n=1,numc
      j=iloc(n)
      write(20,100) 'ATOM  ',n,'  C1  A   A   1',
     &     dfloat(ipsx(j)),dfloat(ipsy(j)),0.d0,       !  z=0
     &     dfloat(issb(j)),dfloat(idsb(j))/5000.
      write(21,*) 'C     ',
     &     dfloat(ipsx(n)),dfloat(ipsy(n)),0.d0        !  z=0
      end do
      write(20,100) 'END   '
  100 format(a6,i5,a15,4x,3f8.3,2f6.2) 
      return
      end
!------------------------------------------------------------------
      subroutine smappa(itime,ncmax,iloc,iocc)
!- temporary routine to produce a numbered map of cell growth
      implicit none
      integer*4 :: ntime_cy, itime, ncmax, numc, indupl, iss
      integer*4 :: i, j, k, l, n, nm
      integer*4 :: iocc(ncmax),iloc(ncmax)

      n=int(dsqrt(dfloat(ncmax)))
      if (n*n.ne.ncmax) stop 888
      nm=n/2
      
      write(30,99) itime
      do j=-nm,nm
      k=(nm+j)*n+nm+1
      write(30,100) (iocc(k+i),i=-nm,nm)
 100  format(50i4)
  99  format(200('*'),'time=',i6)
      end do

      return
      end
