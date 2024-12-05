#if NDIM==3
subroutine star_formation(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  use cooling_module, ONLY: XH=>X
  use constants, only: Myr2sec, Gyr2sec, mH, pi, rhoc, twopi
  use random
  use mpi_mod
#ifdef RT
#ifdef NTRACEGROUPS
  use rt_parameters, only: minind,rt_tracemassbins
  use SED_module,     only: get_tracer_group_from_halo
#endif
#endif 
  implicit none
#ifndef WITHOUTMPI
  integer::info
  integer,parameter::tag=1120
#endif
  integer::ilevel
  !----------------------------------------------------------------------
  ! Description: This subroutine spawns star-particle of constant mass
  ! using a Poisson probability law if some gas condition are fulfilled.
  ! It modifies hydrodynamic variables according to mass conservation
  ! and assumes an isothermal transformation...
  ! On exit, the gas velocity and sound speed are unchanged.
  ! New star particles are synchronized with other collisionless particles.
  ! Array flag2 is used as temporary work space.
  ! Yann Rasera  10/2002-01/2003
  !----------------------------------------------------------------------
  ! local constants
  real(dp)::d0,mgas,mcell
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:twotondim,1:3)::xc
  ! other variables
  integer :: ncache,nnew,ivar,ngrid,icpu,index_star
  integer :: igrid,ix,iy,iz,ind,i,n,iskip,nx_loc
  integer :: ntot,ntot_all,nstar_corrected,ncell
  logical :: ok_free
  real(dp):: d,x,y,z,u,v,w,e
  real(dp):: mstar,dstar,tstar,nISM,nCOM,phi_t,theta,sigs
  real(dp):: T2,nH,T_poly,cs2
  real(dp):: ul,ur,divv,curlv
  real(dp):: birth_epoch,factG
  real(kind=8):: mlost_all,mtot_all
#ifndef WITHOUTMPI
  real(kind=8)::mlost,mtot
#endif
  real(kind=8)::PoissMean
  real(dp),dimension(1:3)::skip_loc
  real(dp):: dx,dx_loc,scale,vol_loc,dx_min,vol_min,d1,d2,d3,d4,d5,d6
  real(dp),dimension(1:nvector)::sfr_ff, alpha_fk=0.0
  real(dp):: alpha0,trgv,scrit,e_cts
  real(dp):: d_gmc,lamjt
  real(dp):: mstar_nsn,scale_msun
  integer ,dimension(1:ncpu,1:IRandNumSize)    :: allseed
  integer ,dimension(1:nvector),save           :: ind_grid,ind_cell,nstar
  integer ,dimension(1:nvector),save           :: ind_grid_new,ind_cell_new,ind_part
  integer ,dimension(1:nvector,0:twondim)      :: ind_nbor
  logical ,dimension(1:nvector),save           :: ok,ok_new=.true.
  integer ,dimension(1:ncpu)                   :: ntot_star_cpu,ntot_star_all
  logical::isConvergent
  integer::idir
  real(dp),dimension(0:twondim) ::darr,uarr,varr,warr
  real(dp)::px_div,py_div,pz_div,Jx,Jy,Jz,rho_local
  integer  ::idx3,idx3_new
!#ifdef RT
!  real(dp):: f_H2
!#endif
#ifdef SOLVERmhd
  real(dp)::bx1,bx2,by1,by2,bz1,bz2
#endif
#if NENER>0
  integer::irad
#endif
#ifdef NCHEM
  real(dp),dimension(1:nchem) :: chem1=0.0
  integer  :: iche
#endif
#ifdef NTRACEGROUPS
  integer::dumpid
#endif
#ifdef POP3
  ! for PopIII stars
  real(dp) :: r3pc
  logical  :: pop3_ok
  real(dp),dimension(1:1000),save              :: mass_pop3
  integer ,dimension(1:1000),save              :: idx3icell
#endif

  real(dp) :: zg,uavg,vavg,wavg,dtot,mstar1,alpha
  real(dp) :: vl,wl,vr,wr

  ! TODO: when f2008 is obligatory - remove this and replace erfc_pre_f08 below by
  ! the f2008 intrinsic erfc() function:
  real(dp) erfc_pre_f08
  
  if (numbtot(1,ilevel)==0) return
  if (.not. hydro) return
  if (ndim.ne.3) return

  if (ilevel.lt.nlevelmax_current) return ! Only form stars at levelmax

  if (verbose) write(*,*)' Entering star_formation'
  
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_msun = scale_l**3*scale_d/1.989d33

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim
  dx_min=(0.5D0**nlevelmax)*scale
  vol_min=dx_min**ndim

  ! ISM density threshold from H/cc to code units
  nISM = n_star
  if(cosmo)then
     nCOM = del_star*omega_b*rhoc*(h0/100)**2/aexp**3*XH/mH
     nISM = MAX(nCOM,nISM)
  endif
  d0   = nISM/scale_nH

  !------------------------------------------------------
  ! Set the star particle mass from the number of SN [TK]
  !------------------------------------------------------
  if(nsn2mass>0.and.fstar_min<0)then
     ! Mesh spacing
     mstar_nsn  = (nsn2mass*M_SNII)/eta_sn/scale_msun
     ! ISM density threshold from H/cc to code units    
     mstar      = n_star/(scale_nH*aexp**3)*vol_min
     fstar_min  = mstar_nsn/mstar
     if(myid==1) write(*,*) ">>>TKNOTE: Mstar,min=",mstar_nsn*scale_msun,fstar_min
  endif

  ! Initial star particle mass
  if(m_star < 0d0)then
     mstar=n_star/(scale_nH*aexp**3)*vol_min*fstar_min
  else
     mstar=m_star*mass_sph
  endif
  dstar=mstar/vol_loc

  factG = 1d0
  if(cosmo) factG = 3d0/4d0/twopi*omega_m*aexp

  ! Birth epoch as proper time
  if(use_proper_time)then
     birth_epoch=texp
  else
     birth_epoch=t
  endif
  d_gmc = MAX(nCOM,n_gmc)/scale_nH

  ! Cells center position relative to grid center position
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     xc(ind,1)=(dble(ix)-0.5D0)*dx
     xc(ind,2)=(dble(iy)-0.5D0)*dx
     xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  ! If necessary, initialize random number generator
  if(localseed(1)==-1)then
     call rans(ncpu,iseed,allseed)
     localseed=allseed(myid,1:IRandNumSize)
  end if

  !------------------------------------------------
  ! Convert hydro variables to primitive variables
  !------------------------------------------------
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        do i=1,ngrid
           d=uold(ind_cell(i),1)
           u=uold(ind_cell(i),2)/d
           v=uold(ind_cell(i),3)/d
           w=uold(ind_cell(i),4)/d
           e=uold(ind_cell(i),5)
#ifdef SOLVERmhd
           bx1=uold(ind_cell(i),6)
           by1=uold(ind_cell(i),7)
           bz1=uold(ind_cell(i),8)
           bx2=uold(ind_cell(i),nvar+1)
           by2=uold(ind_cell(i),nvar+2)
           bz2=uold(ind_cell(i),nvar+3)
           e=e-0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)
#endif
           e=e-0.5d0*d*(u**2+v**2+w**2)
#if NENER>0
           do irad=0,nener-1
              e=e-uold(ind_cell(i),inener+irad)
           end do
#endif
           uold(ind_cell(i),1)=d
           uold(ind_cell(i),2)=u
           uold(ind_cell(i),3)=v
           uold(ind_cell(i),4)=w
           uold(ind_cell(i),5)=e/d
        end do
        do ivar=imetal,nvar
           do i=1,ngrid
              d=uold(ind_cell(i),1)
              w=uold(ind_cell(i),ivar)/d
              uold(ind_cell(i),ivar)=w
           end do
        end do
     end do
  end do

! get values of uold for density and velocities in virtual boundaries
#ifndef WITHOUTMPI
  do ivar=1,4
     call make_virtual_fine_dp(uold(1,ivar),ilevel)
  end do
#endif

  !------------------------------------------------
  ! Compute number of new stars in each cell
  !------------------------------------------------
  ntot=0
  idx3=0  ! index for PopIII stars
  ! Loop over grids
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     ! Star formation criterion ---> logical array ok(i)
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        ! Flag leaf cells
        do i=1,ngrid
           ok(i)=son(ind_cell(i))==0
        end do

        ! Geometrical criterion
        if(ivar_refine>0)then
           do i=1,ngrid
              d=uold(ind_cell(i),ivar_refine)
              if(d<=var_cut_refine)ok(i)=.false.
           end do
        endif

!#ifdef RT
!        ! Check H2 fraction
!        if(H2crit_SF>1d-10 .and. isH2)then
!           do i=1,ngrid
!              f_H2 = (1d0 - uold(ind_cell(i),iIons+ixHII-1) - uold(ind_cell(i),iIons+ixHI-1))
!              if(f_H2<H2crit_SF) ok(i)=.false.
!           end do
!        endif
!#endif

        if (TRIM(star_maker)=='density') then
           do i=1,ngrid
              d=uold(ind_cell(i),1)
              ! Density criterion
              if(d<=d0)ok(i)=.false.
              ! Temperature criterion
              if(ok(i))then
                 T2=uold(ind_cell(i),5)*scale_T2*(gamma-1.0d0)
                 nH=max(uold(ind_cell(i),1),smallr)*scale_nH
                 T_poly=T2_star*(nH/nISM)**(g_star-1.0d0)
                 T2=T2-T_poly
                 if(T2>T2thres_SF)ok(i)=.false.
              endif
           end do          
            
        else if (TRIM(star_maker)=='denmax') then
           ! Density criterion
           do i=1,ngrid
              d=uold(ind_cell(i),1)
              if(d<=d0)then
                 ok(i)=.false.
              else
                 ! we need local density to determine if the cell we are interested 
                 ! is the local density maxima
                 ncell = 1 ! we just want the neighbors of that cell
                 call getnbor(ind_cell(i),ind_nbor,ncell,ilevel)
                 ! First calculate velocities for divergence
                 d1    = uold(ind_nbor(1,1),1)     ; d4    = uold(ind_nbor(1,4),1) 
                 d2    = uold(ind_nbor(1,2),1)     ; d5    = uold(ind_nbor(1,5),1) 
                 d3    = uold(ind_nbor(1,3),1)     ; d6    = uold(ind_nbor(1,6),1)
                 if(d<d1.or.d<d2.or.d<d3.or.d<d4.or.d<d5.or.d<d6)ok(i)=.false.
              endif
           end do           
        else if(TRIM(star_maker)=='hopkins')then
           ! Enforce self-gravitating criterion a la Hopkins et al 2013
           do i=1,ngrid
              ! if cell is a leaf cell
              if (ok(i)) then 
                 d = uold(ind_cell(i),1)
                 ! density threshold to avoid computing the criterion too often but otherwise not needed
                 if (d <= d_gmc) then
                    ok(i) = .false.
                 else
                    ! we need velocity divergence and curl estimates in the cell, so construct values of velocity field 
                    ! on the 6 faces of the cell using simple linear interpolation from neighbouring cell values and differentiate.
                    ! get neighbor cells if they exist, otherwise use straight injection from local cell
                    ncell = 1 ! we just want the neighbors of that cell
                    call getnbor(ind_cell(i),ind_nbor,ncell,ilevel)
                    ! First calculate velocities for divergence
                    d1    = uold(ind_nbor(1,1),1)     ; d4    = uold(ind_nbor(1,4),1) 
                    d2    = uold(ind_nbor(1,2),1)     ; d5    = uold(ind_nbor(1,5),1) 
                    d3    = uold(ind_nbor(1,3),1)     ; d6    = uold(ind_nbor(1,6),1)  
                    ul    = (d2*uold(ind_nbor(1,2),2) + d*uold(ind_cell(i),2))/(d2+d)
                    vl    = (d4*uold(ind_nbor(1,4),3) + d*uold(ind_cell(i),3))/(d4+d)
                    wl    = (d6*uold(ind_nbor(1,6),4) + d*uold(ind_cell(i),4))/(d6+d)
                    ur    = (d1*uold(ind_nbor(1,1),2) + d*uold(ind_cell(i),2))/(d1+d)
                    vr    = (d3*uold(ind_nbor(1,3),3) + d*uold(ind_cell(i),3))/(d3+d)
                    wr    = (d5*uold(ind_nbor(1,5),4) + d*uold(ind_cell(i),4))/(d5+d)
                    divv  = ( ((ur-ul)+(vr-vl)+(wr-wl)) / dx_loc )**2 
                    ! Then calculate velocities for curl
                    vl    = (d6*uold(ind_nbor(1,6),3) + d*uold(ind_cell(i),3))/(d6+d)
                    wl    = (d4*uold(ind_nbor(1,4),4) + d*uold(ind_cell(i),4))/(d4+d)
                    vr    = (d5*uold(ind_nbor(1,5),3) + d*uold(ind_cell(i),3))/(d5+d)
                    wr    = (d3*uold(ind_nbor(1,3),4) + d*uold(ind_cell(i),4))/(d3+d)
                    curlv = ((wr-wl)-(vr-vl))**2         ! first term
                    ul    = (d6*uold(ind_nbor(1,6),2) + d*uold(ind_cell(i),2))/(d6+d)
                    wl    = (d2*uold(ind_nbor(1,2),4) + d*uold(ind_cell(i),4))/(d2+d)
                    ur    = (d5*uold(ind_nbor(1,5),2) + d*uold(ind_cell(i),2))/(d5+d)
                    wr    = (d1*uold(ind_nbor(1,1),4) + d*uold(ind_cell(i),4))/(d1+d)
                    curlv = curlv + ((ur-ul)-(wr-wl))**2 ! second term
                    ul    = (d4*uold(ind_nbor(1,4),2) + d*uold(ind_cell(i),2))/(d4+d)
                    vl    = (d2*uold(ind_nbor(1,2),3) + d*uold(ind_cell(i),3))/(d2+d)
                    ur    = (d3*uold(ind_nbor(1,3),2) + d*uold(ind_cell(i),2))/(d3+d)
                    vr    = (d1*uold(ind_nbor(1,1),3) + d*uold(ind_cell(i),3))/(d1+d)
                    curlv = curlv + ((vr-vl)-(ur-ul))**2 ! third term
                    curlv = curlv / dx_loc**2 
                    ! Check if gas in cell is self-gravitating (alpha < 1)
                    alpha = 0.5d0*(divv+curlv)/(factG*d)
                    if (alpha > 1d0) ok(i) = .false. 
                 endif
              endif
           end do

        else if(TRIM(star_maker)=='federrath')then
           ! Enforce turbulence criterion + efficiency following Federrath & Klessen 2012
           do i=1,ngrid
              ! if cell is a leaf cell
              if (ok(i)) then 
                 d = uold(ind_cell(i),1)
                 ! set density threshold to avoid to compute tensor trace too often but otherwise not needed
                 if (d <= d_gmc) then
                    ok(i) = .false.
                 else
                    ! We need to estimate the norm of the gradient of the velocity field in the cell (tensor of 2nd rank)
                    ! i.e. || A ||^2 = trace( A A^T) where A = grad vec(v) is the tensor. 
                    ! So construct values of velocity field on the 6 faces of the cell using simple linear interpolation 
                    ! from neighbouring cell values and differentiate. 
                    ! Get neighbor cells if they exist, otherwise use straight injection from local cell
                    ncell = 1 ! we just want the neighbors of that cell
                    call getnbor(ind_cell(i),ind_nbor,ncell,ilevel)
                    d1    = uold(ind_nbor(1,1),1)     ; d4    = uold(ind_nbor(1,4),1) 
                    d2    = uold(ind_nbor(1,2),1)     ; d5    = uold(ind_nbor(1,5),1) 
                    d3    = uold(ind_nbor(1,3),1)     ; d6    = uold(ind_nbor(1,6),1)  
                    ul    = (d2*uold(ind_nbor(1,2),2) + d*uold(ind_cell(i),2))/(d2+d)
                    ur    = (d1*uold(ind_nbor(1,1),2) + d*uold(ind_cell(i),2))/(d1+d)
                    trgv  = (ur-ul)**2
                    ul    = (d4*uold(ind_nbor(1,4),3) + d*uold(ind_cell(i),3))/(d4+d)
                    ur    = (d3*uold(ind_nbor(1,3),3) + d*uold(ind_cell(i),3))/(d3+d)
                    trgv  = trgv + (ur-ul)**2
                    ul    = (d6*uold(ind_nbor(1,6),4) + d*uold(ind_cell(i),4))/(d6+d)
                    ur    = (d5*uold(ind_nbor(1,5),4) + d*uold(ind_cell(i),4))/(d5+d)
                    trgv  = trgv + (ur-ul)**2
                    ul    = (d6*uold(ind_nbor(1,6),3) + d*uold(ind_cell(i),3))/(d6+d)
                    ur    = (d5*uold(ind_nbor(1,5),3) + d*uold(ind_cell(i),3))/(d5+d)
                    trgv  = trgv + (ur-ul)**2
                    ul    = (d4*uold(ind_nbor(1,4),4) + d*uold(ind_cell(i),4))/(d4+d)
                    ur    = (d3*uold(ind_nbor(1,3),4) + d*uold(ind_cell(i),4))/(d3+d)
                    trgv  = trgv + (ur-ul)**2
                    ul    = (d6*uold(ind_nbor(1,6),2) + d*uold(ind_cell(i),2))/(d6+d)
                    ur    = (d5*uold(ind_nbor(1,5),2) + d*uold(ind_cell(i),2))/(d5+d)
                    trgv  = trgv + (ur-ul)**2
                    ul    = (d2*uold(ind_nbor(1,2),4) + d*uold(ind_cell(i),4))/(d2+d)
                    ur    = (d1*uold(ind_nbor(1,1),4) + d*uold(ind_cell(i),4))/(d1+d)
                    trgv  = trgv + (ur-ul)**2
                    ul    = (d4*uold(ind_nbor(1,4),2) + d*uold(ind_cell(i),2))/(d4+d)
                    ur    = (d3*uold(ind_nbor(1,3),2) + d*uold(ind_cell(i),2))/(d3+d)
                    trgv  = trgv + (ur-ul)**2
                    ul    = (d2*uold(ind_nbor(1,2),3) + d*uold(ind_cell(i),3))/(d2+d)
                    ur    = (d1*uold(ind_nbor(1,1),3) + d*uold(ind_cell(i),3))/(d1+d)
                    trgv  = trgv + (ur-ul)**2
                    trgv  = trgv 
                    ! now compute sound speed squared = cs^2 
                    cs2  = uold(ind_cell(i),5)*(gamma -1.0)*gamma ! Joki: not totally sure about gamma
#if NENER>0
                    ! Add non-thermal pressure to sound speed calculation
                    do irad = 0,nener-1
                          cs2 = cs2 + uold(ind_cell(i),inener+irad) / d * (gamma_rad(irad+1)-1.0) * gamma_rad(irad+1)
                    end do
#endif
                    ! prevent numerical crash due to negative temperature
                    cs2  = max(cs2,smallc**2)
                    ! Calculate "turbulent" Jeans length in cell units, lamjt 
                    ! (see e.g. Chandrasekhar 51, Bonazzola et al 87, Federrath & Klessen 2012 eq 36)
                    lamjt = (pi*trgv + sqrt(pi*pi*trgv*trgv + 36.0*pi*cs2*factG*d*dx_loc**2))/(6.0*factG*d*dx_loc**2)
                    if (lamjt > sf_lam) then ! Jeans length resolved: gas is stable 
                       ok(i) = .false. 
                    else ! Jeans length not resolved --> form stars to lower density and stabilise gas
                       ! corresponding virial parameter for homogeneous sphere <= 1.5 in turbulent dominated 
                       ! limit and <= 0.5 in pressure dominated limit (in good agreement with observations,
                       ! at least for massive (>= 10^5 M_sun) clouds see Kauffmann, Pillai, Goldsmith 2013, Fig.1)  
                       alpha0 = 5d0/(pi*factG*d)*(trgv + cs2)/dx_loc**2
                       ! compute star formation efficiency per free-fall time (Federrath & Klessen 2012 eq 41) 
                       ! e_cts is the unresolved (for us) proto-stellar feedback: i.e. the dense gas core-to-star efficiency
                       e_cts = 0.5  ! would be 1.0 without feedback (Federrath & Klessen 2012)
                       phi_t = 0.57 ; theta = 0.33 ! best fit values of the Padoan & Nordlund multi-scale sf model to GMC simulation data 
                       sigs  = log(1.0+0.16*trgv/cs2) ! factor 0.16 is b^2 where b=0.4 for a mixture of turbulence forcing modes
                       scrit = log(0.067/theta**2*alpha0*trgv/cs2) ! best fit from Padoan & Nordlund MS model again 
                       sfr_ff(i) = e_cts/2.0*phi_t*exp(3.0/8.0*sigs)*(2.0-erfc_pre_f08((sigs-scrit)/sqrt(2.0*sigs)))
                       alpha_fk(i)=alpha0
                    endif
                 endif
              endif
           end do

        else if(TRIM(star_maker)=='federrath3')then
           ! trgv for turbulence, no Jeans length criterion to form stars
           do i=1,ngrid
              ! if cell is a leaf cell
              if (ok(i)) then 
                 d = uold(ind_cell(i),1)
                 ! density threshold to avoid computing the criterion too often but otherwise not needed
                 if (d <= d_gmc) then
                    ok(i) = .false.
                 else

                    ncell = 1 ! we just want the neighbors of that cell
                    call getnbor(ind_cell(i),ind_nbor,ncell,ilevel)
                    d1    = uold(ind_nbor(1,1),1)     ; d4    = uold(ind_nbor(1,4),1) 
                    d2    = uold(ind_nbor(1,2),1)     ; d5    = uold(ind_nbor(1,5),1) 
                    d3    = uold(ind_nbor(1,3),1)     ; d6    = uold(ind_nbor(1,6),1)  
                    ul    = (d2*uold(ind_nbor(1,2),2) + d*uold(ind_cell(i),2))/(d2+d)
                    ur    = (d1*uold(ind_nbor(1,1),2) + d*uold(ind_cell(i),2))/(d1+d)
                    trgv  = (ur-ul)**2
                    ul    = (d4*uold(ind_nbor(1,4),3) + d*uold(ind_cell(i),3))/(d4+d)
                    ur    = (d3*uold(ind_nbor(1,3),3) + d*uold(ind_cell(i),3))/(d3+d)
                    trgv  = trgv + (ur-ul)**2
                    ul    = (d6*uold(ind_nbor(1,6),4) + d*uold(ind_cell(i),4))/(d6+d)
                    ur    = (d5*uold(ind_nbor(1,5),4) + d*uold(ind_cell(i),4))/(d5+d)
                    trgv  = trgv + (ur-ul)**2
                    ul    = (d6*uold(ind_nbor(1,6),3) + d*uold(ind_cell(i),3))/(d6+d)
                    ur    = (d5*uold(ind_nbor(1,5),3) + d*uold(ind_cell(i),3))/(d5+d)
                    trgv  = trgv + (ur-ul)**2
                    ul    = (d4*uold(ind_nbor(1,4),4) + d*uold(ind_cell(i),4))/(d4+d)
                    ur    = (d3*uold(ind_nbor(1,3),4) + d*uold(ind_cell(i),4))/(d3+d)
                    trgv  = trgv + (ur-ul)**2
                    ul    = (d6*uold(ind_nbor(1,6),2) + d*uold(ind_cell(i),2))/(d6+d)
                    ur    = (d5*uold(ind_nbor(1,5),2) + d*uold(ind_cell(i),2))/(d5+d)
                    trgv  = trgv + (ur-ul)**2
                    ul    = (d2*uold(ind_nbor(1,2),4) + d*uold(ind_cell(i),4))/(d2+d)
                    ur    = (d1*uold(ind_nbor(1,1),4) + d*uold(ind_cell(i),4))/(d1+d)
                    trgv  = trgv + (ur-ul)**2
                    ul    = (d4*uold(ind_nbor(1,4),2) + d*uold(ind_cell(i),2))/(d4+d)
                    ur    = (d3*uold(ind_nbor(1,3),2) + d*uold(ind_cell(i),2))/(d3+d)
                    trgv  = trgv + (ur-ul)**2
                    ul    = (d2*uold(ind_nbor(1,2),3) + d*uold(ind_cell(i),3))/(d2+d)
                    ur    = (d1*uold(ind_nbor(1,1),3) + d*uold(ind_cell(i),3))/(d1+d)
                    trgv  = trgv + (ur-ul)**2
                    trgv  = trgv

                    cs2  = uold(ind_cell(i),5)*(gamma -1.0)*gamma ! Joki: not totally sure about gamma
#if NENER>0
                    ! Add radiation pressure to sound speed calculation
                    do irad = 0,nener-1
                          cs2 = cs2 + uold(ind_cell(i),inener+irad) / d * (gamma_rad(irad+1)-1.0) * gamma_rad(irad+1)
                    end do
#endif
                    cs2  = max(cs2,smallc**2)
                    alpha0 = 5d0/(pi*factG*d)*(trgv + cs2)/dx_loc**2
                    e_cts = 0.5
                    phi_t = 0.57 ; theta = 0.33
                    sigs  = log(1.0+0.16*trgv/cs2)
                    scrit = log(0.067/theta**2*alpha0*trgv/cs2)
                    sfr_ff(i) = e_cts/2.0*phi_t*exp(3.0/8.0*sigs)*(2.0-erfc_pre_f08((sigs-scrit)/sqrt(2.0*sigs)))
                    alpha_fk(i)=alpha0
                 endif ! if d>d_gmc
              endif ! if leaf cell
           enddo

        else if(TRIM(star_maker)=='FK2')then
           ! Enforce turbulence criterion + efficiency following Federrath & Klessen 2012
           do i=1,ngrid
              ! if cell is a leaf cell and high density
              d = uold(ind_cell(i),1)
              ok(i) = ok(i).and.(d>=d_gmc)
              if (ok(i)) then 
                 ! We need to estimate the norm of the gradient of the velocity field in the cell (tensor of 2nd rank)
                 ! i.e. || A ||^2 = trace( A A^T) where A = grad vec(v) is the tensor. 
                 ! So construct values of velocity field on the 6 faces of the cell using simple linear interpolation 
                 ! from neighbouring cell values and differentiate. 
                 ! Get neighbor cells if they exist, otherwise use straight injection from local cell
                 ncell = 1 ! we just want the neighbors of that cell
                 call getnbor(ind_cell(i),ind_nbor,ncell,ilevel)

                 darr(0) = d
                 do idir=1,6
                    darr(idir) = uold(ind_nbor(1,idir),1)
                 end do

                 ! local density maxima
                 if(ok(i)) then
                   if(maxval(darr)>d) ok(i)=.false.
                 endif

                 ! converging flow check
                 if(ok(i))then
                    uarr(0) = uold(ind_cell(i),2)
                    varr(0) = uold(ind_cell(i),3)
                    warr(0) = uold(ind_cell(i),4)
                    do idir=1,6
                        uarr(idir) = uold(ind_nbor(1,idir),2) 
                        varr(idir) = uold(ind_nbor(1,idir),3) 
                        warr(idir) = uold(ind_nbor(1,idir),4) 
                    end do
                    divv  = (uarr(2)*darr(2)-uarr(1)*darr(1)) &
                        & + (varr(4)*darr(4)-varr(3)*darr(3)) & 
                        & + (warr(6)*darr(6)-warr(5)*darr(5))
                    if(divv>0) ok(i) = .false.  ! diverging flow
                 endif

                 if(ok(i)) then
                    !average velocity
                    dtot  = sum(darr)
                    uavg  = sum(darr*uarr)/dtot
                    vavg  = sum(darr*varr)/dtot
                    wavg  = sum(darr*warr)/dtot

                    ! subtract the mean velocity field
                    uarr(:) = uarr(:) - uavg
                    varr(:) = varr(:) - vavg
                    warr(:) = warr(:) - wavg

                    ! subtract the symmetric divergence field                    
                    ! ex)  (---->,<--): only subtract (-->,<--): result (-->,0) 
                    ! ex)  (<----,-->): only subtract (<--,-->): result (<--,0)
                    px_div = min( abs(darr(1)*uarr(1)),abs(darr(2)*uarr(2)))
                    py_div = min( abs(darr(3)*varr(3)),abs(darr(4)*varr(4)))
                    pz_div = min( abs(darr(5)*warr(5)),abs(darr(6)*warr(6)))

                    isConvergent = darr(2)*uarr(2) - darr(1)*uarr(1) < 0 
                    if (isConvergent) then
                       uarr(1) = uarr(1) - px_div/darr(1)
                       uarr(2) = uarr(2) + px_div/darr(2)
                    else ! comment out if you do not want to subtract outflows
                       uarr(1) = uarr(1) + px_div/darr(1)
                       uarr(2) = uarr(2) - px_div/darr(2)
                    endif 

                    isConvergent = darr(4)*varr(4) - darr(3)*varr(3) < 0
                    if (isConvergent) then
                       varr(3) = varr(3) - py_div/darr(3)
                       varr(4) = varr(4) + py_div/darr(4)
                    else ! comment out if you do not want to subtract outflows
                       varr(3) = varr(3) + py_div/darr(3)
                       varr(4) = varr(4) - py_div/darr(4)
                    endif

                    isConvergent = darr(6)*warr(6) - darr(5)*warr(5) < 0
                    if (isConvergent) then 
                       warr(5) = warr(5) - pz_div/darr(5)
                       warr(6) = warr(6) + pz_div/darr(6)
                    else ! comment out if you do not want to subtract outflows
                       warr(5) = warr(5) + pz_div/darr(5)
                       warr(6) = warr(6) - pz_div/darr(6)
                    endif

                    ! subtract the rotational velocity field (x-y) plane
                    ! ^y       <-        |4|        |-u|
                    ! |       |  |     |1| |2|   |-v|  |+v|
                    ! --->x    ->        |3|        |+u|
                    Jz  = - varr(1)*darr(1) + varr(2)*darr(2) &
                      &   + uarr(3)*darr(3) - uarr(4)*darr(4)
                    Jz  = Jz / 4.0

                    varr(1) = varr(1) + Jz/darr(1) 
                    varr(2) = varr(2) - Jz/darr(2) 
                    uarr(3) = uarr(3) - Jz/darr(3)
                    uarr(4) = uarr(4) + Jz/darr(4)

                    ! subtract the rotational velocity field (y-z) plane
                    ! ^z       <-        |6|        |-v|  
                    ! |       |  |     |3| |4|   |-w|  |+w|
                    ! --->y    ->        |5|        |+v|
                    Jx  = - warr(3)*darr(3) + warr(4)*darr(4) &
                      &   + varr(5)*darr(5) - varr(6)*darr(6)
                    Jx  = Jx / 4.0

                    warr(3) = warr(3) + Jx/darr(3) 
                    warr(4) = warr(4) - Jx/darr(4) 
                    varr(5) = varr(5) - Jx/darr(5)
                    varr(6) = varr(6) + Jx/darr(6)

                    ! subtract the rotational velocity field (x-z) plane
                    ! ^z       ->        |6|        |+u|  
                    ! |       |  |     |1| |2|   |+w|  |-w|
                    ! --->x    <-        |5|        |-u|
                    Jy  = + warr(1)*darr(1) - warr(2)*darr(2) &
                      &   - uarr(5)*darr(5) + uarr(6)*darr(6)
                    Jy  = Jy / 4.0

                    warr(1) = warr(1) - Jy/darr(1) 
                    warr(2) = warr(2) + Jy/darr(2) 
                    uarr(5) = uarr(5) + Jy/darr(5)
                    uarr(6) = uarr(6) - Jy/darr(6)

                    ! From this point, uarr,uarr,warr is just the turbulent velocity
                    trgv  = 0.0

                    !x-direc
                    ul    = (darr(2)*uarr(2) + d*uarr(0))/(darr(2)+d)
                    ur    = (darr(1)*uarr(1) + d*uarr(0))/(darr(1)+d)
                    trgv  = trgv + (ur-ul)**2
                    !y-direc
                    ul    = (darr(4)*varr(4) + d*varr(0))/(darr(4)+d)
                    ur    = (darr(3)*varr(3) + d*varr(0))/(darr(3)+d)
                    trgv  = trgv + (ur-ul)**2
                    !z-direc
                    ul    = (darr(6)*warr(6) + d*warr(0))/(darr(6)+d)
                    ur    = (darr(5)*warr(5) + d*warr(0))/(darr(5)+d)
                    trgv  = trgv + (ur-ul)**2
                    !z-direc; tangential component - y
                    ul    = (darr(6)*varr(6) + d*varr(0))/(darr(6)+d)
                    ur    = (darr(5)*varr(5) + d*varr(0))/(darr(5)+d)
                    trgv  = trgv + (ur-ul)**2
                    !y-direc; tangential component - z
                    ul    = (darr(4)*warr(4) + d*warr(0))/(darr(4)+d)
                    ur    = (darr(3)*warr(3) + d*warr(0))/(darr(3)+d)
                    trgv  = trgv + (ur-ul)**2
                    !z-direc; tangential component - x
                    ul    = (darr(6)*uarr(6) + d*uarr(0))/(darr(6)+d)
                    ur    = (darr(5)*uarr(5) + d*uarr(0))/(darr(5)+d)
                    trgv  = trgv + (ur-ul)**2
                    !x-direc; tangential component - z
                    ul    = (darr(2)*warr(2) + d*warr(0))/(darr(2)+d)
                    ur    = (darr(1)*warr(1) + d*warr(0))/(darr(1)+d)
                    trgv  = trgv + (ur-ul)**2
                    !y-direc; tangential component - x
                    ul    = (darr(4)*uarr(4) + d*uarr(0))/(darr(4)+d)
                    ur    = (darr(3)*uarr(3) + d*uarr(0))/(darr(3)+d)
                    trgv  = trgv + (ur-ul)**2
                    !x-direc; tangential component - y
                    ul    = (darr(2)*varr(2) + d*varr(0))/(darr(2)+d)
                    ur    = (darr(1)*varr(1) + d*varr(0))/(darr(1)+d)
                    trgv  = trgv + (ur-ul)**2
                    ! now compute sound speed squared = cs^2 
                    cs2  = uold(ind_cell(i),5)*(gamma -1.0)*gamma ! Joki: not totally sure about gamma
                    ! prevent numerical crash due to negative temperature
#if NENER>0
                    ! Add radiation pressure to sound speed calculation
                    do irad = 0,nener-1
                          cs2 = cs2 + uold(ind_cell(i),inener+irad) / d * (gamma_rad(irad+1)-1.0) * gamma_rad(irad+1)
                    end do
#endif
                    cs2  = max(cs2,smallc**2)
                    ! Calculate "turbulent" Jeans length in cell units, lamjt 
                    ! (see e.g. Chandrasekhar 51, Bonazzola et al 87, Federrath & Klessen 2012 eq 36)
                    !rho_local = d
                    rho_local = rho(ind_cell(i)) ! be careful with dark matter particles
                      
                    lamjt = (pi*trgv + sqrt(pi*pi*trgv*trgv + 36.0*pi*cs2*factG*rho_local*dx_loc**2))/(6.0*factG*rho_local*dx_loc**2)
                    if (lamjt > sf_lam) then ! Jeans length resolved: gas is stable 
                       ok(i) = .false. 
                    else ! Jeans length not resolved --> form stars to lower density and stabilise gas
                       ! corresponding virial parameter for homogeneous sphere <= 1.5 in turbulent dominated 
                       ! limit and <= 0.5 in pressure dominated limit (in good agreement with observations,
                       ! at least for massive (>= 10^5 M_sun) clouds see Kauffmann, Pillai, Goldsmith 2013, Fig.1)  
                       alpha0 = 5d0/(pi*factG*rho_local)*(trgv + cs2)/dx_loc**2
                       alpha_fk(i) =  alpha0
                       ! compute star formation efficiency per free-fall time (Federrath & Klessen 2012 eq 41) 
                       ! e_cts is the unresolved (for us) proto-stellar feedback: i.e. the dense gas core-to-star efficiency
                       e_cts = 0.5  ! would be 1.0 without feedback (Federrath & Klessen 2012)
                       phi_t = 0.57 ; theta = 0.33 ! best fit values of the Padoan & Nordlund multi-scale sf model to GMC simulation data 
                       sigs  = log(1.0+0.16*trgv/cs2) ! factor 0.16 is b^2 where b=0.4 for a mixture of turbulence forcing modes
                       scrit = log(0.067/theta**2*alpha0*trgv/cs2) ! best fit from Padoan & Nordlund MS model again 
                       sfr_ff(i) = e_cts/2.0*phi_t*exp(3.0/8.0*sigs)*(2.0-erfc_pre_f08((sigs-scrit)/sqrt(2.0*sigs)))
                    endif

                 endif
              endif
           end do

        else if(TRIM(star_maker)=='dencon')then
           ! Cen & Ostriker scheme (converging flow)
           do i=1,ngrid
              ! if cell is a leaf cell
              if (ok(i)) then 
                 d  = uold(ind_cell(i),1)
                 if (d <= d0) then  ! density threshold
                    ok(i) = .false.
                 endif

                 if(ok(i))then ! converging flow check
                    ncell = 1 ! we just want the neighbors of that cell
                    call getnbor(ind_cell(i),ind_nbor,ncell,ilevel)
                    ! First calculate velocities for divergence
                    d1    = uold(ind_nbor(1,1),1)     ; d4    = uold(ind_nbor(1,4),1) 
                    d2    = uold(ind_nbor(1,2),1)     ; d5    = uold(ind_nbor(1,5),1) 
                    d3    = uold(ind_nbor(1,3),1)     ; d6    = uold(ind_nbor(1,6),1) 
                    if(d1>d.or.d2>d.or.d3>d.or.d4>d.or.d5>d.or.d6>d) ok(i)=.false.
                    if(ok(i))then 
                       ul    = d1*uold(ind_nbor(1,1),2)
                       ur    = d2*uold(ind_nbor(1,2),2)
                       vl    = d3*uold(ind_nbor(1,3),3)
                       vr    = d4*uold(ind_nbor(1,4),3)
                       wl    = d5*uold(ind_nbor(1,5),4)
                       wr    = d6*uold(ind_nbor(1,6),4)
                       divv  = (ur-ul)+(vr-vl)+(wr-wl)
                       if(divv>0) ok(i) = .false.  ! diverging flow
                    endif
                 endif

              endif
           enddo !i=1,ngrid
        else
           if(myid.eq.1)then 
              write(*,*) ' star_maker does not match any of the pre-defined functions'
              write(*,*) star_maker 
           endif
           stop
        endif ! star maker


        ! Calculate number of new stars in each cell using Poisson statistics
        do i=1,ngrid
           nstar(i)=0
           if(ok(i))then
              ! Compute mean number of events
              d         = uold(ind_cell(i),1)
              mcell     = d*vol_loc
              tstar     = .5427*sqrt(1.0/(factG*max(d,smallr)))
              ! this is the free fall time of an homogeneous sphere 
              sfr_ff(i) = 1.0 
              if(TRIM(star_maker).ne.'hopkins' .and. TRIM(star_maker).ne.'federrath' &
                   .and. TRIM(star_maker).ne.'federrath3' .and. TRIM(star_maker).ne.'FK2') then
                 sfr_ff(i) = eps_star
              endif

#ifdef POP3
              ! for PopIII star formation
              if(metal) zg = uold(ind_cell(i),imetal)
              pop3_ok = pop3.and.(zg.lt.Zcrit_pop3)
              if(pop3_ok) then
                 if(pop3_mass>0)then
                    mstar1 = pop3_mass/scale_msun
                    PoissMean = dtnew(ilevel)*sfr_ff(i)/tstar*mcell/mstar1
                 else
                    call pop3_mass_random(mstar1)
                    mstar1 = mstar1/scale_msun
                    PoissMean = dtnew(ilevel)*sfr_ff(i)/tstar*mcell/(100./scale_msun) ! 100: characteristic scale, for random sampling
                 endif
              else
#endif
                 mstar1 = mstar
                 PoissMean = dtnew(ilevel)*sfr_ff(i)/tstar*mcell/mstar1
#ifdef POP3
              endif
#endif
              !! If catastrophic star formation (massive star cluster) wants to occur, we need to limit the 
              !! maximal mass of the star particle we want to create in a cell.
              !PoissMean = min(PoissMean,10.0)
              ! Compute Poisson realisation
              call poissdev(localseed,PoissMean,nstar(i))
              ! Compute depleted gas mass
              mgas      = nstar(i)*mstar1
#ifdef POP3
              ! Allow 1 popIII star formation per event
              if(pop3_ok.and.(nstar(i).ge.1)) then
                 if(mcell*scale_msun.ge.1d3.and.mgas<0.9*mcell) then ! for random sampling
                    nstar(i)=1
                    idx3 = idx3 + 1
                    if(idx3.gt.10000)then
                       write(*,*)' ERROR in star formation: increase PopIII reservoir'
                       call clean_stop
                    endif
                    mass_pop3(idx3) = mstar1
                    idx3icell(idx3) = ind_cell(i)
                 else
                    nstar(i)=0
                 endif
              endif 
#endif

              ! Security to prevent more than 90% of gas depletion
              if (mgas > 0.9d0*mcell) then
                 nstar_corrected = int(0.9d0*mcell/mstar1)
                 mstar_lost      = mstar_lost+(nstar(i)-nstar_corrected)*mstar1
                 nstar(i)        = nstar_corrected
              endif

#ifdef POP3
              ! check if there is any active Pop III stars within 3 pc radius
              r3pc = 3.*3.08d18/(scale_l*boxlen)
              r3pc = max(r3pc,3*0.5**nlevelmax)
              if(nstar(i)>0)then ! stupid way of doing this...but..
                 do kk=1,npartmax
                    if(idp(kk)<0)then
                       if(zp(kk).lt.Zcrit_pop3)then
                          x = (xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
                          y = (xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
                          z = (xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
                          if(  (abs(xp(kk,1)-x).le.r3pc).and.&
                          &    (abs(xp(kk,2)-y).le.r3pc).and.&
                          &    (abs(xp(kk,3)-z).le.r3pc)) nstar(i)=0
                       endif 
                    endif
                 end do
              endif
#endif

              ! Compute new stars local statistics
              mstar_tot=mstar_tot+nstar(i)*mstar1
              if (nstar(i)>0) then
                 ntot = ntot+1
                 if(log_sf) then
399                  format('SF z= ',f9.5,1x,' nH= ',f7.3,' N= ',I5,' eps_ff= ',f9.5,' alpha= ',f9.5)
                     write(*,399) 1./aexp-1,log10(d*scale_nH),nstar(i),sfr_ff(i),alpha_fk(i)
                 endif
              else
                 nstar(i)=0
              endif
           endif
        enddo
        ! Store nstar in array flag2
        do i=1,ngrid
           flag2(ind_cell(i))=nstar(i)
        end do
     end do
  end do

  !---------------------------------
  ! Check for free particle memory
  !---------------------------------
  ok_free = ((numbp_free-ntot)>=0)
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(numbp_free,numbp_free_tot,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  numbp_free_tot=numbp_free
#endif
  if(.not. ok_free)then
     write(*,*)'No more free memory for particles'
     write(*,*)'Increase npartmax'
#ifndef WITHOUTMPI
    call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
    stop
#endif
  end if

  !---------------------------------
  ! Compute global stars statistics
  !---------------------------------
#ifndef WITHOUTMPI
  mlost=mstar_lost; mtot=mstar_tot
  call MPI_ALLREDUCE(ntot,ntot_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(mtot,mtot_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(mlost,mlost_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  ntot_all=ntot
  mtot_all=mstar_tot
  mlost_all=mstar_lost
#endif
  ntot_star_cpu=0; ntot_star_all=0
  ntot_star_cpu(myid)=ntot
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(ntot_star_cpu,ntot_star_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  ntot_star_cpu(1)=ntot_star_all(1)
#endif
  do icpu=2,ncpu
     ntot_star_cpu(icpu)=ntot_star_cpu(icpu-1)+ntot_star_all(icpu)
  end do
  nstar_tot = nstar_tot+ntot_all
  if (myid==1) then
     if (ntot_all.gt.0) then
        !print *,'##################################################################################'
        write(*,'(" Level=",I6," New star=",I6," Tot=",I10," Mass=",1PE10.3," Lost=",0PF4.1,"%")')&
             & ilevel,ntot_all,nstar_tot,mtot_all,mlost_all/mtot_all*100.
        !print *,'##################################################################################'
     endif
  end if

  !------------------------------
  ! Create new star particles
  !------------------------------
  ! Starting identity number
  if(myid==1)then
     index_star=nstar_tot-ntot_all
  else
     index_star=nstar_tot-ntot_all+ntot_star_cpu(myid-1)
  end if

  ! Loop over grids
  idx3_new = 1  ! index for PopIII (ind_cell_new)
  idx3   = 1  ! index for PopIII (ind_cell)
  ncache = active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do

     ! Loop over cells
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do

        ! Flag cells with at least one new star
        do i=1,ngrid
           ok(i)=flag2(ind_cell(i))>0
        end do

        ! Gather new star arrays
        nnew=0
        do i=1,ngrid
           if (ok(i))then
              nnew=nnew+1
              ind_grid_new(nnew)=ind_grid(i)
              ind_cell_new(nnew)=ind_cell(i)
           end if
        end do

        ! Update linked list for stars
        call remove_free(ind_part,nnew)
        call add_list(ind_part,ind_grid_new,ok_new,nnew)

        ! Calculate new star particle and modify gas density
        do i=1,nnew
           index_star=index_star+1

           ! Get gas variables
           n=flag2(ind_cell_new(i))
           d=uold(ind_cell_new(i),1)
           u=uold(ind_cell_new(i),2)
           v=uold(ind_cell_new(i),3)
           w=uold(ind_cell_new(i),4)
           x=(xg(ind_grid_new(i),1)+xc(ind,1)-skip_loc(1))*scale
           y=(xg(ind_grid_new(i),2)+xc(ind,2)-skip_loc(2))*scale
           z=(xg(ind_grid_new(i),3)+xc(ind,3)-skip_loc(3))*scale
           if (metal) then
              zg = uold(ind_cell_new(i),imetal)
#ifdef NCHEM
              do iche=1,nchem
                 chem1(iche) = uold(ind_cell_new(i),ichem+iche-1)
              enddo
#endif
           endif
           mstar1 = mstar
#ifdef POP3
           if(pop3)then
              if(ind_cell_new(i).eq.idx3icell(idx3_new))then
                 mstar1 = mass_pop3(idx3_new)
                 idx3_new = idx3_new + 1
              endif
           endif
#endif

           ! Set star particle variables
           tp(ind_part(i)) = birth_epoch     ! Birth epoch
           mp(ind_part(i)) = n*mstar1        ! Mass
           if(use_initial_mass) then
              mp0(ind_part(i)) = mp(ind_part(i)) ! Initial Mass
           endif
           levelp(ind_part(i)) = ilevel      ! Level
           idp(ind_part(i))    = -index_star ! Star identity
           typep(ind_part(i))%family = FAM_STAR
           typep(ind_part(i))%tag = 0
           xp(ind_part(i),1)   = x
           xp(ind_part(i),2)   = y
           xp(ind_part(i),3)   = z
           vp(ind_part(i),1)   = u
           vp(ind_part(i),2)   = v
           vp(ind_part(i),3)   = w
           if (metal) then
              zp(ind_part(i)) = zg  ! Initial star metallicity
#ifdef NCHEM
              do iche=1,nchem
                 chp(ind_part(i),iche) = chem1(iche)  ! Initial chemical abudance
              enddo
#endif
           endif
#ifdef RT
#ifdef NTRACEGROUPS
     dumpid = 0
     if (rt_tracemassbins(1).gt.-1.d0) then
        call get_tracer_group_from_halo(xp(ind_part(i),1:ndim),dumpid,minind)
     endif
     ptracegroup(ind_part(i)) = dumpid
#endif
#endif

        end do
        ! End loop over new star particles

        ! Modify gas density according to mass depletion
        do i=1,ngrid
           if (flag2(ind_cell(i))>0) then
              n = flag2(ind_cell(i))
              d = uold(ind_cell(i),1)

              ! for PopIII stars
              mstar1 = mstar
#ifdef POP3
              if(pop3)then
                 if(ind_cell(i).eq.idx3icell(idx3))then
                    mstar1 = mass_pop3(idx3)
                    idx3 = idx3 + 1
                 endif
              endif
#endif
              dstar = mstar1/vol_loc

              uold(ind_cell(i),1) = d-n*dstar
           endif
        end do

     end do
     ! End loop over cells
  end do
  ! End loop over grids

  !---------------------------------------------------------
  ! Convert hydro variables back to conservative variables
  !---------------------------------------------------------
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        do i=1,ngrid
           d=uold(ind_cell(i),1)
           u=uold(ind_cell(i),2)
           v=uold(ind_cell(i),3)
           w=uold(ind_cell(i),4)
           e=uold(ind_cell(i),5)*d
#ifdef SOLVERmhd
           bx1=uold(ind_cell(i),6)
           by1=uold(ind_cell(i),7)
           bz1=uold(ind_cell(i),8)
           bx2=uold(ind_cell(i),nvar+1)
           by2=uold(ind_cell(i),nvar+2)
           bz2=uold(ind_cell(i),nvar+3)
           e=e+0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)
#endif
           e=e+0.5d0*d*(u**2+v**2+w**2)
#if NENER>0
           do irad=0,nener-1
              e=e+uold(ind_cell(i),inener+irad)
           end do
#endif
           uold(ind_cell(i),1)=d
           uold(ind_cell(i),2)=d*u
           uold(ind_cell(i),3)=d*v
           uold(ind_cell(i),4)=d*w
           uold(ind_cell(i),5)=e
        end do
        do ivar=imetal,nvar
           do i=1,ngrid
              d=uold(ind_cell(i),1)
              w=uold(ind_cell(i),ivar)
              uold(ind_cell(i),ivar)=d*w
           end do
        end do
     end do
  end do

end subroutine star_formation 
#endif
!################################################################
!################################################################
!################################################################
!################################################################
subroutine getnbor(ind_cell,ind_father,ncell,ilevel)
  use amr_commons
  implicit none
  integer::ncell,ilevel
  integer,dimension(1:nvector)::ind_cell
  integer,dimension(1:nvector,0:twondim)::ind_father
  !-----------------------------------------------------------------
  ! This subroutine determines the 2*ndim neighboring cells
  ! cells of the input cell (ind_cell).
  ! If for some reasons they don't exist, the routine returns
  ! the input cell.
  !-----------------------------------------------------------------
  integer::i,j,iok,ind
  integer,dimension(1:nvector),save::ind_grid_father,pos
  integer,dimension(1:nvector,0:twondim),save::igridn,igridn_ok
  integer,dimension(1:nvector,1:twondim),save::icelln_ok


  if(ilevel==1)then
     write(*,*) 'Warning: attempting to form stars on level 1 --> this is not allowed ...'
     return
  endif

  ! Get father cell
  do i=1,ncell
     ind_father(i,0)=ind_cell(i)
  end do

  ! Get father cell position in the grid
  do i=1,ncell
     pos(i)=(ind_father(i,0)-ncoarse-1)/ngridmax+1
  end do

  ! Get father grid
  do i=1,ncell
     ind_grid_father(i)=ind_father(i,0)-ncoarse-(pos(i)-1)*ngridmax
  end do

  ! Get neighboring father grids
  call getnborgrids(ind_grid_father,igridn,ncell)

  ! Loop over position
  do ind=1,twotondim

     ! Select father cells that sit at position ind
     do j=0,twondim
        iok=0
        do i=1,ncell
           if(pos(i)==ind)then
              iok=iok+1
              igridn_ok(iok,j)=igridn(i,j)
           end if
        end do
     end do

     ! Get neighboring cells for selected cells
     if(iok>0)call getnborcells(igridn_ok,ind,icelln_ok,iok)

     ! Update neighboring father cells for selected cells
     do j=1,twondim
        iok=0
        do i=1,ncell
           if(pos(i)==ind)then
              iok=iok+1
              if(icelln_ok(iok,j)>0)then
                 ind_father(i,j)=icelln_ok(iok,j)
              else
                 ind_father(i,j)=ind_cell(i)
              end if
           end if
        end do
     end do

  end do


end subroutine getnbor
!##############################################################
!##############################################################
!##############################################################
!##############################################################
function erfc_pre_f08(x)

! complementary error function
  use amr_commons, ONLY: dp
  implicit none
  real(dp) erfc_pre_f08
  real(dp) x, y
  real(kind=8) pv, ph
  real(kind=8) q0, q1, q2, q3, q4, q5, q6, q7
  real(kind=8) p0, p1, p2, p3, p4, p5, p6, p7
  parameter(pv= 1.26974899965115684d+01, ph= 6.10399733098688199d+00)
  parameter(p0= 2.96316885199227378d-01, p1= 1.81581125134637070d-01)
  parameter(p2= 6.81866451424939493d-02, p3= 1.56907543161966709d-02)
  parameter(p4= 2.21290116681517573d-03, p5= 1.91395813098742864d-04)
  parameter(p6= 9.71013284010551623d-06, p7= 1.66642447174307753d-07)
  parameter(q0= 6.12158644495538758d-02, q1= 5.50942780056002085d-01)
  parameter(q2= 1.53039662058770397d+00, q3= 2.99957952311300634d+00)
  parameter(q4= 4.95867777128246701d+00, q5= 7.41471251099335407d+00)
  parameter(q6= 1.04765104356545238d+01, q7= 1.48455557345597957d+01)

  y = x*x
  y = exp(-y)*x*(p7/(y+q7)+p6/(y+q6) + p5/(y+q5)+p4/(y+q4)+p3/(y+q3) &
       &       + p2/(y+q2)+p1/(y+q1)+p0/(y+q0))
  if (x < ph) y = y+2.0d0/(exp(pv*x)+1.0d0)
  erfc_pre_f08 = y

  return
  
end function erfc_pre_f08
!##############################################################
!##############################################################
!##############################################################
subroutine MaxwellianCDF(xx,ywant,fn,df)
   implicit none
   real(kind=8)::xx, fn, df, ywant
   real(kind=8)::pi=3.14159265358979323
   fn = derf(xx) - 2d0/dsqrt(pi)*xx*dexp(-xx*xx) - ywant
   df = dsqrt(8d0/pi)*xx*xx*dexp(-xx*xx)
end subroutine MaxwellianCDF
!##############################################################
!##############################################################
!##############################################################
subroutine newton_raphson_safe(x1,x2,ywant,xeps,rtsafe)
   IMPLICIT NONE
   REAL(KIND=8)::x1,x2,ywant,xeps,rtsafe
   !Using a combination of Newton-Raphson and bisection, 
   !and the root of a function bracketed between x1 and x2. 
   !The root, returned as the function value rtsafe, will be returned 
   !until its accuracy is known within xeps. funcd is a user-supplied 
   !subroutine which returns both the function value and the 
   !first derivative of the function.
   INTEGER::j,MAXIT=100
   REAL(KIND=8)::df,dx,dxold,f,fh,fl,temp,xh,xl

   call MaxwellianCDF(x1,ywant,fl,df)
   call MaxwellianCDF(x2,ywant,fh,df)
   if((fl.gt.0..and.fh.gt.0.).or.(fl.lt.0..and.fh.lt.0.)) &
      write(*,*) 'root must be bracketed in rtsafe'
   if(fl.eq.0.)then
      rtsafe=x1
      return
   else if(fh.eq.0.)then
      rtsafe=x2
      return
   else if(fl.lt.0.)then !Orient the search so that f(xl) < 0.
      xl=x1
      xh=x2
   else
      xh=x1
      xl=x2
   endif
   rtsafe=.5*(x1+x2) !Initialize the guess for root,
   dxold=dabs(x2-x1) !the \stepsize before last,"
   dx=dxold !and the last step.
   call MaxwellianCDF(rtsafe,ywant,f,df)
   do j=1,MAXIT !Loop over allowed iterations.
      if((((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f).gt.0)& !Bisect if Newton out of range,
         .or.(dabs(2d0*f).gt.dabs(dxold*df)) ) then !or not decreasing fast enough.
         dxold=dx
         dx=0.5d0*(xh-xl)
         rtsafe=xl+dx
         if(xl.eq.rtsafe)return !Change in root is negligible.
      else !Newton step acceptable. Take it.
         dxold=dx
         dx=f/df
         temp=rtsafe
         rtsafe=rtsafe-dx
         if(temp.eq.rtsafe)return
      endif
      if(dabs(dx).lt.xeps) return !Convergence criterion.
      call MaxwellianCDF(rtsafe,ywant,f,df) !The one new function evaluation per iteration.
      if(f.lt.0.) then !Maintain the bracket on the root.
         xl=rtsafe
      else
         xh=rtsafe
      endif
    enddo
    write(*,*) 'rtsafe exceeding maximum iterations'
end subroutine newton_raphson_safe
!##############################################################
!##############################################################
!##############################################################
subroutine pop3_analytic (m,proba)
   implicit none
   real(kind=8)::m,mcha,proba
   mcha = 100d0
   proba = m**(-1.3)*exp(-(mcha/m)**(1.6))
end subroutine pop3_analytic
!##############################################################
!##############################################################
!##############################################################
subroutine pop3_mass_random (mass)
   use amr_commons, ONLY:dp
   implicit none
   real(kind=dp)::mass
   real(kind=8)::fx_cap,m,logmin,logmax,fx,y
   logical::not_found
   integer,save::iseed=-399778
   real(kind=8),external::ran1
   m = 100. ! charateristic mass
   call pop3_analytic(m, fx_cap)
   fx_cap = fx_cap * 2

   logmin = 1.0
   logmax = 3.0

   not_found=.true.
   do while(not_found)
      m = ran1(iseed)
      m = m*(logmax-logmin)+logmin
      m = 10d0**m
      call pop3_analytic(m, fx)
      y = fx_cap*ran1(iseed)
      if (y.le.fx) not_found=.false.
   end do
   mass = m
end subroutine pop3_mass_random
!##############################################################
!##############################################################
!##############################################################
