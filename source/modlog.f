      PROGRAM MODLOG
!  MAIN PROGRAM OF THE CELL SIMULATION CODE
!                                                         .FC 04/2014
!-----------------------------------------------------------------------------------
      implicit none
      integer*4, parameter ::  nst=1000, ndt=5
      integer*4, parameter ::  ntime=10000,nsize=999,ncmax=nsize*nsize
!- careful to set always ncmax=(nsize)^2 with nsize=integer+1
      integer*4 :: ncel, ntime_ir, ntime_cy, itconv, kc, kt, kj, lc
      integer*4 :: ikt, jkt, izt, jzt, index, indupl, iss, itime, jtime
      integer*4 :: ncrit, ndead
      integer, allocatable :: issb(:),idsb(:),ipsx(:),ipsy(:),ipsz(:)
      integer, allocatable :: idupl(:),ityp(:),iloc(:),iocc(:),idens(:)
      integer, allocatable :: iphse(:),idead(:),innb(:,:)
      real*8 :: dose, ssb, dsb, rssb, rdsb, duprb, dupl_time, cc, prob 
      real*8 :: dcoff, dffc, death, dcsi
      real*8, allocatable :: ps(:,:),pd(:,:),rps(:,:),rpd(:,:) 
      real*8, allocatable :: dens(:),strl(:),binom(:,:)
            allocate ( strl(0:nst), binom(0:nst,0:nst) )
            allocate ( dens(1:ncmax),idead(1:ncmax) )
            allocate ( issb(1:ncmax),idsb(1:ncmax) )
            allocate ( ipsx(1:ncmax),ipsy(1:ncmax),ipsz(1:ncmax) )
            allocate ( idupl(1:ncmax),ityp(1:ncmax),iloc(1:ncmax) )
            allocate ( iocc(1:ncmax),idens(1:ncmax),iphse(1:ncmax) )
            allocate ( innb(1:ncmax,1:6) )          ! near neighb table (+/-x,+/-y,+/-z)
            allocate ( ps(0:nst,0:nst), pd(0:ndt,0:ndt) ) 
            allocate ( rps(0:nst,0:nst), rpd(0:ndt,0:ndt) )
      ps = 0.d0 ; pd = 0.d0 ; rps = 0.d0 ; rpd = 0.d0
      issb = 0 ; idsb = 0 ; ipsx = 0 ; ipsy = 0 ; ipsz = 0 
      idupl = 0 ; ityp = 1 ; iphse = 1 ; iocc = 0 ; iloc = 0
      ndead = 0 ; dens = 0.d0
 
      if (ndt.gt.nst) stop 999
!- input data (either from line or file)
      ncel=500                      ! initial # of cells
      dose=2.d0                     ! dose (Gy)
      ntime_ir=120                  ! irr. time (seconds)
      ntime_cy=1000                 ! cell time (minutes)
      itconv=60                     ! conv. sec to min
      ssb=dose/dfloat(ntime_ir)     ! ssb prob = 1000/Gy
      dsb=0.2d0*ssb                 ! dsb prob
      rssb=0.27d-2                  ! ssb repair prob. 
      rdsb=0.02d0*rssb              ! dsb repair prob.
      duprb=100.                    ! cell duplication probability
      dffc=1.                       ! diffusion coefficient for tracer
      ncrit=5                       ! critical DSB threshold for death
      OPEN(unit=5,file='input',status='old')
      read(5,*) ncel, dose, ntime_ir, ntime_cy, itconv,
     &          ssb,dsb,rssb,rdsb,duprb,dffc,ncrit

!-  compute factorial w Stirling approx
      call facto_s(nst,strl)
      
!-  compute binomial coefficient w Stirling approx
      call binom_s(nst,strl,binom)

!-  build probability matrix for ssb
      call prob_s(nst,ps,binom,ssb)

!-  build probability matrix for dsb
      call prob_d(nst,ndt,pd,binom,dsb)

!-  build probability matrix for ssb repair
      call prob_rs(nst,rps,binom,rssb)

!-  build probability matrix for dsb repair
      call prob_rd(nst,ndt,rpd,binom,rdsb)

!-  assign spatial location to cells
      call space(ncmax,ncel,ipsx,ipsy,ipsz,iloc,iocc,innb)

!-------------------IRRADIATION PHASE-----------------------
!- start irradiation time iteration
!- time evolution of cell pop. on unit 15

      do itime=1,ntime_ir

!- try make SSBs in every cell w probability ssb
      call ssb_rad(nst,ncel,ps,issb)
 
!- try make DSBs in every cell w probability dsb
      call dsb_rad(ndt,ncel,pd,idsb)

!- count cells above ssb and dsb threshold
!- count average ssb/dsb per cell
      ikt=0 ; jkt=0 ; izt=0 ; jzt=0
      do kt=1,ncel
      ikt=ikt+issb(kt)/nst
      jkt=jkt+issb(kt)
      izt=izt+idsb(kt)/ndt
      jzt=jzt+idsb(kt)
      end do
      write(15,*) itime,ikt,izt,dfloat(jkt)/dfloat(ncel),
     $                          dfloat(jzt)/dfloat(ncel) 

      end do      !    close loop on irradiation time

!- plot map of DNA breaks
!     call map(ncmax,ncel,issb,idsb,iloc,ipsx,ipsy,ipsz)
!- plot map of tracer diffusion
      call map(ncmax,ncel,issb,idens,iloc,ipsx,ipsy,ipsz)

!-------------------CELL CYCLE PHASE-----------------------
!- start cell cycle time iteration
      DO ITIME=1,NTIME_CY

!- cell cycle phase
      do kc=1,ncel
      if (idsb(kc).gt.0) then     !  cell arrested if z1>0
          iphse(kc)=0
          if (idead(kc).eq.0) then
          dcoff=(dfloat(idsb(kc))/dfloat(ncrit))**2
          death=dmin1(1.d0,dcoff)
          dcsi=rand()
          if (dcsi.le.death) then !  check for cell death
              idead(kc)=1
              ndead=ndead+1
          end if ; end if
      end if
      if (iphse(kc).eq.0) go to 66
      jtime=itime-idupl(kc)
      if (iphse(kc).eq.1)
     &    idupl(kc)=idupl(kc)-(rand()/10.d0)*11  !  this anticipates slightly the G1 phase
      lc=jtime-1440*(jtime/1440)
      if (lc.gt.600)  iphse(kc)=2    !  phase S
      if (lc.gt.1000) iphse(kc)=3    !  phase G2
      if (lc.gt.1300) iphse(kc)=4    !  phase M
  66  continue
c     write (16,*) itime,kc,lc,iphse(kc),idupl(kc)
      end do
!- try repair SSBs in damaged cells w probability rssb
      if (idead(kc).eq.0) call ssb_rep(nst,ncel,rps,issb)

!- try repair DSBs in damaged cells w probability rdsb
      if (idead(kc).eq.0) call dsb_rep(ndt,ncel,rpd,idsb)

!- count cells above ssb and dsb threshold
      ikt=0 ; jkt=0 ; izt=0 ; jzt=0
      do kt=1,ncel
      ikt=ikt+issb(kt)/nst
      jkt=jkt+issb(kt)
      izt=izt+idsb(kt)/ndt
      jzt=jzt+idsb(kt)
      end do
      write(15,*) ntime_ir+itconv*itime,
     $            ikt,izt,dfloat(jkt)/dfloat(ncel),
     $            dfloat(jzt)/dfloat(ncel) 

!- test for tracer diffusion
      call fick(dens,dffc,itime,ncmax,ncel,idead,iloc,iocc,
     &                                            idens,innb)

!- test for cell duplication
      call double(duprb,ntime_cy,itime,ncmax,ncel,iphse,
     &            issb,idsb,iloc,iocc,ityp,innb,idupl,indupl)
 
!     write (6,101) itime, ncel, indupl
      if(mod(itime,100).eq.0) write (6,*) itime, ncel, indupl, ndead
  101 format('time=',i6,' population=',i6,' dupl/time=',i6)
!- print pdb map of ssb and dsb every 100 time steps
!     if (mod(itime,100).eq.0) 
!    &    call map(ncmax,ncel,issb,idsb,iloc,ipsx,ipsy,ipsz)
!- print map of tracer diffusion every 100 time steps
      if (mod(itime,100).eq.0)
     &    call map(ncmax,ncel,issb,idens,iloc,ipsx,ipsy,ipsz)

      END DO      !    loop on cell cycle time

      stop
      end
!-----------------------------------------------------------------
      subroutine fick(dens,dffc,itime,ncmax,numc,idead,
     &   iloc,iocc,idens,innb)
!- simulate diffusion of a tracer among cells

      implicit none
      integer*4 ::  numc, ncmax, itime, i, j, k, kc, kt, kj, lc
      integer*4 :: iloc(ncmax),ityp(ncmax),iocc(ncmax),innb(ncmax,6)
      integer*4 :: idens(ncmax),idead(ncmax)
      real*8 :: dc, dtot, dffc, zdif, dens(ncmax)
      real*8, allocatable :: dens_t(:), dens_s(:)
            allocate ( dens_t(1:ncmax),dens_s(1:ncmax) )
            dens_t = 0.d0 ; dens_s = 0.d0 ; idens = 0
      do i=1,numc
      kc=iloc(i)            !  kc is site occupied by cell i
      if (idead(kc).eq.0) then
        do lc=1,4           !  loop over 4 nabors
          zdif=1.d0
          kt=innb(kc,lc)    !  kt is site occupied by nabor lc of cell i
          j=iocc(kt)        !  j is nabor lc of cell i
          if (j.gt.0) then
            if (idead(j).gt.0) zdif=0.d0   !  if receiving cell is dead D=0
            dc=(dens(i)-dens(j))*dffc*zdif
            dens_t(i)=dens_t(i)-dc 
            dens_t(j)=dens_t(j)+dc
          end if 
        end do 
      end if
      end do

      dtot=0.d0
      dens_s(130)=dfloat(itime)*5.
      do i=1,numc
      dens(i)=dens(i)+dens_t(i)
     &                         +dens_s(i)   !  source term
      j=iloc(i)
!     idens(i)=nint(dens(i)) ! careful: use idens(i) to plot cell property
      idens(j)=nint(dens(i)) !           or idens(j) to plot site property
      dtot=dtot+dens_t(i)
      end do
!     if (dabs(dtot).gt.1.d-10) stop 777
!     call smappa(ncmax,iloc,idens)   !  plot integer part of density
!     call smappa(ncmax,iloc,iocc)    !  plot map of cell numbers
      return
      end
