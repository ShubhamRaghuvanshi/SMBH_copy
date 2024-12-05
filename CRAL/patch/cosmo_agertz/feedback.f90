!################################################################
!################################################################
!################################################################
!################################################################
#if NDIM==3
subroutine thermal_feedback(ilevel)
  use pm_commons
  use amr_commons
  use hydro_commons
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info2,dummy_io
#endif
  integer::ilevel
  !------------------------------------------------------------------------
  ! This routine computes the thermal energy, the kinetic energy and
  ! the metal mass dumped in the gas by stars (SNII, SNIa, winds).
  ! This routine is called every fine time step.
  !------------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part,ivar
  integer::ig,ip,npart1,npart2,icpu,ilun,idim
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part
  character(LEN=80)::filename,filedir,fileloc,filedirini
  character(LEN=5)::nchar,ncharcpu
  logical::file_exist
  integer,parameter::tag=1120

  if(sf_log_properties) then
     call title(ifout-1,nchar)
     if(IOGROUPSIZEREP>0) then
        call title(((myid-1)/IOGROUPSIZEREP)+1,ncharcpu)
        filedirini='output_'//TRIM(nchar)//'/'
        filedir='output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/'
     else
        filedir='output_'//TRIM(nchar)//'/'
     endif
     filename=TRIM(filedir)//'stars_'//TRIM(nchar)//'.out'
     ilun=myid+103
     call title(myid,nchar)
     fileloc=TRIM(filename)//TRIM(nchar)
     ! Wait for the token
#ifndef WITHOUTMPI
     if(IOGROUPSIZE>0) then
        if (mod(myid-1,IOGROUPSIZE)/=0) then
           call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
                & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
        end if
     endif
#endif

     inquire(file=fileloc,exist=file_exist)
     if(.not.file_exist) then
        open(ilun, file=fileloc, form='formatted')
        write(ilun,'(A24)',advance='no') '# event id  ilevel  mp  '
        do idim=1,ndim
           write(ilun,'(A2,I1,A2)',advance='no') 'xp',idim,'  '
        enddo
        do idim=1,ndim
           write(ilun,'(A2,I1,A2)',advance='no') 'vp',idim,'  '
        enddo
        do ivar=1,nvar
           if(ivar.ge.10) then
              write(ilun,'(A1,I2,A2)',advance='no') 'u',ivar,'  '
           else
              write(ilun,'(A1,I1,A2)',advance='no') 'u',ivar,'  '
           endif
        enddo
        write(ilun,'(A5)',advance='no') 'tag  '
        write(ilun,'(A1)') ' '
     else
        open(ilun, file=fileloc, status="old", position="append", action="write", form='formatted')
     endif
  endif


  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Gather star particles only.

  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0

        ! Count star particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if ( is_star(typep(ipart)) ) then
                 npart2=npart2+1
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif

        ! Gather star particles
        if(npart2>0)then
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              ! Select only star particles
              if ( is_star(typep(ipart)) ) then
                 if(ig==0)then
                    ig=1
                    ind_grid(ig)=igrid
                 end if
                 ip=ip+1
                 ind_part(ip)=ipart
                 ind_grid_part(ip)=ig
              endif
              if(ip==nvector)then
                 call feedbk(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
                 ip=0
                 ig=0
              end if
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if
        igrid=next(igrid)   ! Go to next grid
     end do
     ! End loop over grids
     if(ip>0)call feedbk(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
  end do
  ! End loop over cpus

  if(sf_log_properties) close(ilun)

111 format('   Entering thermal_feedback for level ',I2)

end subroutine thermal_feedback
#endif
!################################################################
!################################################################
!################################################################
!################################################################
#if NDIM==3
subroutine feedbk(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use random
  use constants, only: M_sun, Myr2sec, pc2cm
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !-----------------------------------------------------------------------
  ! This routine is called by subroutine feedback. Each stellar particle
  ! dumps mass, momentum and energy in the nearest grid cell using array
  ! unew.
  !-----------------------------------------------------------------------
  integer::i,j,idim,nx_loc,ilun,iii,ii,jj,kk
  real(kind=8)::RandNum
  real(dp)::dx_min,vol_min
  real(dp)::t0,ESN,mejecta,zloss,time_simu
  real(dp)::tauIR
  real(dp)::dx,dx_loc,scale,birth_time
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::L1,Cr,KappaIR,KappaIR_0
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle based arrays
  logical,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector),save::mloss,mzloss,ethermal,dteff
  real(dp),dimension(1:nvector)::mlossFe,mlossO
  real(dp),dimension(1:nvector)::ptot
  real(dp),dimension(1:nvector),save::vol_loc,dx_loc2
  real(dp),dimension(1:nvector,1:ndim),save::x
  integer ,dimension(1:nvector,1:ndim),save::id,igd,icd
  integer ,dimension(1:nvector),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc
#if NENER>0
  integer::irad
#endif
  real(dp)::MBmin,MBmax,Mremnant,SNIametals,numII,ENSN,SNIaFe,SNIaO
  real(dp)::mstarmin,mstarmax,NumSNIa
  real(dp)::masslossIa,ESNIa
  real(dp)::masslossW,metalsW
  real(dp)::theint,t1,t2
  real(dp)::Mwindmin,Mwindmax,mett
  real(dp)::vol_loc2
  real(dp):: NSNIa,MIMF,SNIIFe,SNIIO,SNIIej,NSNII,MIMFChabrier,IMFChabrier
  external NSNIa,MIMF,SNIIFe,SNIIO,SNIIej,NSNII,MIMFChabrier,IMFChabrier
  real(dp)::eta1,eta2,alpha1,alpha2,beta,Cr1,Cr2,mumax,eps_cl,tcl,Mclmin,Mclmax
  real(dp)::alpha,mtrans,Zscale,Zgas,tcut
  integer::iicell,iskip
  integer,dimension(1:nvector,1:2,1:2,1:2)::indcube2  !Oscar
  real(dp)::xcont,ycont,zcont,contr
  real(dp)::Lum,p10  
  real(dp)::tt,tekin,vkick,maxPrad,maxcl,maxTau,maxMet
  integer ,dimension(1:ncpu,1:IRandNumSize)::allseed
  real(dp),dimension(1:nvector)::mcl,Prad,agecl,Lumcl
  integer ,dimension(1:nvector)::ibinning
  integer::irad,icenter
  real(dp)::vmax,momx,momy,momz,pST,n0,imfboost
  integer::indpmax,iradmax
  integer,dimension(:),allocatable::ind
  integer ,dimension(1:nvector)::indrad
  real(dp)::momxt,momyt,momzt,momtot,momtott,pmiss,clustermass,clusterage,meanmass,meanmassM
  real(dp)::SNmin,SNmax,Tmaxfb,rsf,cellsize
  real(dp)::maxadv,vxold,vyold,vzold,vxnew,vynew,vznew,Emax
  integer::indd
  real(dp)::rhofb,scale_m,twind,numresidual
!  real(dp),dimension(1:nvector,1:5)::workarray

  if(sf_log_properties) ilun=myid+103
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m=1.98892d33/scale_d/scale_l/scale_l/scale_l !Msun into grams and then internal units

!  If necessary, initialize random number generator                                           
  if(localseed(1)==-1)then
     call rans(ncpu,iseed,allseed)
     localseed=allseed(myid,1:IRandNumSize)
  end if

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  dx_loc2(1:nvector)=dx*scale
  vol_loc(1:nvector)=dx_loc**ndim   !new approach
  vol_loc2=dx_loc**ndim
  dx_min=(0.5D0**nlevelmax)*scale
  vol_min=dx_min**ndim


! ---------------------------- SN feedback parameters ----------------------------
  ! Type Ia supernova parameters
  MBmin=3.0d0 !Observ. parameter of binaries going SNIa
  MBmax=16.0d0 !Observ. parameter of binaries going SNIa
  Mremnant=1.4d0   !Chandrasekhar mass, nothing remains
  SNIaFe=0.63d0
  SNIaO=0.13d0
  SNIametals=SNIaFe+SNIaO  
  !Mass in Msun of enriched material (Raiteri et al. 1996) 0.63 Msun Fe + 0.13Msun O (e.g. Theilemann 1986)

  ! Scale to internal units 
  SNIametals=SNIametals*scale_m
  SNIaFe=SNIaFe*scale_m
  SNIaO=SNIaO*scale_m
  Mremnant=Mremnant*scale_m

  ! AGB wind parameters
  Mwindmin=0.5d0
  Mwindmax=8.0d0
  twind=1.0d7 !'fast winds'

  ! Massive star lifetime from Myr to code units
  t0=t_sne*1d6*(365.*24.*3600.)/scale_t

!  SNenergy=1.d51 !erg/SN
  ! Type II supernova specific energy from cgs to code units
  ESN=E_SNIa/(10.*2d33)/scale_v/scale_v  !energy per 10 Msun in internal units 
  ENSN=E_SNIa/scale_d/scale_l/scale_l/scale_l/scale_v/scale_v   !energy in internal units
  p10=12.0d0*3000.0d5*scale_m/scale_v !12 Msun at 3000 km/s per SNII event, calibrated to get same \dot{p} as SB99 
  ESNIa=ENSN  !Oscar: assume same energy release as SNII

  SNmin=8.0d0
  SNmax=80.0d0   !Limits for SNII events

  ! Radiation pressure
  eta1=2.0d0 !See paper  
  eta2=eta_rap
  alpha2=0.0d0   
  alpha1=0.4d0  
  beta=1.7d0  !spectrum slope
  mumax=1.0d0  
  eps_cl=0.2d0 !Cluster formation efficiency, observed  
  tcl=6.0d6 !assumed clump lifetime in years 
  tcut=3.0d6 !3.0d6 !From constant Lumonisity to powerlaw                 
  mtrans=3.0d4*scale_m  
  Cr2=2.5d0*3.08568d18/scale_l !2.5 parsec in internal units
  Cr1=Cr2/(mtrans)**0.4 
  imfboost=0.3143d0/0.224468d0 !Kroupa to Chabrier
  KappaIR_0=5.0d0*scale_d*scale_l !g-1 cm2 into internal units. Depending on dust mix, it can be as high as 30 in IR
  Mclmin=100.0*scale_m  !100 Msol, internal units

  ! Radiation pressure, bolometric
  L1=imfboost*scale_t*1.9d-7/scale_v !3.0  !Specific Lbol/c: Lbol/1d6 Msun/c in CGS converted into internal. Units are cm/s^2

!-------------------------------------------------------------------------------------

  vmax=vmaxFB
  Tmaxfb=T2max
  maxadv=maxadvfb*1d5/scale_v !into internal units

!--- Oscar: move age caluclation here to only do book-keeping if necessary?

  ! Lower left corner of 3x3x3 grid-cube
  do idim=1,ndim
     do i=1,ng
        x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
     end do
  end do

  ! Gather 27 neighboring father cells (should be present anytime !)
  do i=1,ng
     ind_cell(i)=father(ind_grid(i))
  end do
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ng,ilevel)

  ! Rescale position at level ilevel
  do idim=1,ndim
     do j=1,np
        x(j,idim)=xp(ind_part(j),idim)/scale+skip_loc(idim)
        x(j,idim)=x(j,idim)-x0(ind_grid_part(j),idim)
        x(j,idim)=x(j,idim)/dx
        id(j,idim)=int(x(j,idim))   ! NGP at level ilevel
        igd(j,idim)=id(j,idim)/2    ! Compute parent grids
     end do
  end do
 
  ok(1:np)=.true.
  do j=1,np
     kg(j)=1+igd(j,1)+3*igd(j,2)+9*igd(j,3)
     igrid(j)=son(nbors_father_cells(ind_grid_part(j),kg(j)))
     ok(j)=ok(j).and.igrid(j)>0    ! Check if particles are entirely in level ilevel
  end do

  ! Compute parent cell position
  do idim=1,ndim
     do j=1,np
        if(ok(j))then
           icd(j,idim)=id(j,idim)-2*igd(j,idim)
        end if
     end do
  end do
  do j=1,np
     if(ok(j))then
        icell(j)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
     end if
  end do

  indpmax=1
  ! Compute parent cell adresses
  do j=1,np
     if(ok(j))then
        indp(j)=ncoarse+(icell(j)-1)*ngridmax+igrid(j)
      else
        indp(j)= nbors_father_cells(ind_grid_part(j),kg(j))  
        vol_loc(j)=vol_loc(j)*2**ndim 
        dx_loc2(j)=dx_loc2(j)*2.0    
     end if
    ! max address 
     indpmax=max(indpmax,indp(j)) 
  end do

  allocate(ind(1:indpmax))
  ind=0.0
  indrad=0.0

! Loop over cells 
  indd=0
  do kk=1,2
     do jj=1,2
        do ii=1,2
           indd=indd+1  !counter from 1,twotondim
           iskip=ncoarse+(indd-1)*ngridmax
           do j=1,np
              ind_cell(j)=iskip+igrid(j)        !Oscar: Oct indices. do feedback on these
              indcube2(j,ii,jj,kk)=ind_cell(j)  !Oscar change back
           end do
        enddo
     enddo
  enddo

!-------------------------------------------------------------

   ! Compute individual time steps
  do j=1,np
     dteff(j)=dtnew(levelp(ind_part(j)))   !Oscar: also new (c.f. old versions)
     if(use_proper_time)then
        dteff(j)=dteff(j)*aexp**2  !dteff is in conformal
     endif
  end do

  ! Reset ejected mass, metallicity, thermal energy
     mloss(:)=0.d0
     mzloss(:)=0.d0
     ethermal(:)=0.d0
     Prad(:)=0.0d0
     ptot(:)=0.0d0
     mlossFe(:)=0.0 !new
     mlossO(:)=0.0  !new

 if(cosmo) then
     ! Find neighboring expansion factors                                 
     i=1
     do while(aexp_frw(i)>aexp.and.i<n_frw)
        i=i+1
     end do
     ! Interploate time
     time_simu=t_frw(i)*(aexp-aexp_frw(i-1))/(aexp_frw(i)-aexp_frw(i-1))+ &
          & t_frw(i-1)*(aexp-aexp_frw(i))/(aexp_frw(i-1)-aexp_frw(i))
    endif

! --------- Prad binning---------
   maxPrad=0.0
   maxcl=0.0
   maxTau=0.0
   maxMet=0.0
   irad=0
   mcl=0.0
   ibinning=0
   agecl=0.0
   Lumcl=0.0
! -------------------------------
!   workarray=0.0

  ! Compute feedback
   do j=1,np
       mejecta=0.0
       numII=0.0
       zloss=0.0

       if(cosmo) then
          if(use_proper_time)then    !Correct as tp is in proper, and dteff is proper. Now get to yrs
             call getAgeGyr(tp(ind_part(j)), t1)          !  End-of-dt age [Gyrs]
             call getAgeGyr(tp(ind_part(j))-dteff(j), t2) !  End-of-dt age [Gyrs]
             t1=t1*1.0d9
             t2=t2*1.0d9
         else
           ! Compute star age in years,conformal time to years
           iii=1
           do while(tau_frw(iii)>tp(ind_part(j)).and.iii<n_frw)
              iii=iii+1
           end do
           
           t2=t_frw(iii)*(tp(ind_part(j))-tau_frw(iii-1))/(tau_frw(iii)-tau_frw(iii-1))+ &                       
                & t_frw(iii-1)*(tp(ind_part(j))-tau_frw(iii))/(tau_frw(iii-1)-tau_frw(iii))
           t2=(time_simu-t2)/(h0*1.d5/3.08d24)/(365.*24.*3600.)                        !Units of years; particle age 
           
           t1=t_frw(iii)*(tp(ind_part(j))+dteff(j)-tau_frw(iii-1))/(tau_frw(iii)-tau_frw(iii-1))+ &
                & t_frw(iii-1)*(tp(ind_part(j))+dteff(j)-tau_frw(iii))/(tau_frw(iii-1)-tau_frw(iii))
           t1=(time_simu-t1)/(h0*1.d5/3.08d24)/(365.*24.*3600.)                        !Units of years; particle age - dt  
           
           birth_time=tp(ind_part(j))
           t1=max(t1,0.0d0) 
           t2=max(t2,0.0d0)            
        endif
        else
           t2=t-tp(ind_part(j))           !Age at t
           t2=t2*scale_t/365./24./60./60. !For non-cosmo units, in years
           t1=t-tp(ind_part(j))-dteff(j)  !Age at t-dt, in years
           t1=t1*scale_t/365./24./60./60.
           birth_time=tp(ind_part(j))     !internal units
           t1=max(t1,0.0d0)
        endif
        
!--------------------------- Get dying stars ----------------------
        !mett=zp(ind_part(j))

!  ---------- metals is star and host cell
        ! imetal=Fe
        ! imetal+1=O
        mett=2.09d0*zp2(ind_part(j))+1.06d0*zp(ind_part(j))    !Asplund, to be used throughout for total Z assuming solar mixture
        Zgas=(2.09d0*unew(indp(j),imetal+1)+1.06d0*unew(indp(j),imetal))/ &
             & max(unew(indp(j),1),smallr)/0.02d0 !Average, solar mix in gas (Asplund)
        Zgas=max(Zgas,0.0)
     
        if (mett<0.0004) then ! Raiteri et al. 1996 limits
           mett=0.0004
        endif
        if(mett>0.05) then
           mett=0.05
        endif
        if(t1.gt.0.0)then
           call agemass(t1,mett,mstarmax)  !Get stellar masses that exit main sequence during dt
           call agemass(t2,mett,mstarmin)  
        else
           mstarmin=120.d0
           mstarmax=120.d0
        endif
        if(mstarmin.le.SNmax.and.mstarmax.ge.SNmax) then
           mstarmax=SNmax-1d-2
        endif
        if(mstarmin.le.SNmin.and.mstarmax.ge.SNmin) then
           mstarmin=SNmin+1d-2
        endif
        if(mstarmin.le.SNmin.and.mstarmax.ge.SNmax) then !if timestep is >40 Myr
           mstarmin=SNmin+1d-2
           mstarmax=SNmax-1d-2
        endif
        meanmass=(mstarmax+mstarmin)/2.0d0
!---------------------------
!--------------------------- Supernovae
!---------------------------
        if(supernovae) then
!---------------------------
!--------------------------- Type II events
!---------------------------
           if(mstarmax.le.SNmax.and.mstarmin.ge.SNmin) then 
              call SNIInum(mstarmin,mstarmax,theint)  !no numerical integral
              numII=mp0(ind_part(j))*theint/scale_m 
              call ranf(localseed,RandNum) 
              numresidual=numII-int(numII)  !int <1 --> 0
              numII=int(numII)
              if(RandNum<numresidual) then
                 numII=numII+1       ! Another star from residual
              endif

              if(numII>0.0) then

                 Zscale=(max(Zgas,0.01d0))**(-0.2)                    !---------- Thornton et al. 1998.
                 n0=max(unew(indp(j),1),smallr)*scale_nH              !---------- mH/cc   
                 pST=4.0d5*numII**(0.941)*n0**(-0.1176)*Zscale        !---------- km/sST momentum Blondin et al (1998). 
                 pST=pST*1.0d5*scale_m/scale_v                        !---------- Into g cm/s and then into code units
                 !--------------------------------------------------------------- Check cooling radius
                 rsf=30.*numII**0.29*(n0**(-0.43))*((Zgas+0.01)**(-0.18)) !------ Shell formation radius (Hopkins et al. 2013, Coiffi, et al.)
                 cellsize=dx_loc2(j)*scale_l/3.08d18                  !---------- local resolution element
                 
                 if(cellsize.lt.rsf/Nrcool)then                       !----------- Kim & Ostriker criterion. Momentum: 1/3, energy: 1/10. 
                    momST=.false.
                 else
                    momST=.true.
                 endif

                 rhofb=1.d-4 !hardcoded (..and very low)
                 if(n0<rhofb)momST=.false.                            !--- Optional momentum limiter in very diffuse gas

                 if(momST) then
                    ptot(j)=ptot(j)+pST                               !------- Post Sedov Taylor momentum
                 else
                    ptot(j)=ptot(j)+p10*numII                         !------- Conservative initial bblastwave momentum
                 endif

                 ! -------  Stellar mass loss and SNII massloading
                 mejecta=numII*(0.7682*meanmass**1.056)*scale_m       !------- Woosley Weaver 1995, Raiteri 1996 
                 mloss(j)=mloss(j)+mejecta/vol_loc(j)
                 ethermal(j)=ethermal(j)+numII*ENSN/vol_loc(j)        !------- ENSN per SNII event (cooled away in low res cases)

                 ! --- Metallicity (here only Fe and O, Raiteri 1996)
                 meanmassM=min(meanmass,40.0)                         !------ yields from very massive stars assumed to be same as
                 mlossFe(j)=mlossFe(j)+numII*0.375d0*exp(-17.94d0/meanmassM)*scale_m/vol_loc(j)  !------ Wosley & Heger (2007)
                 mlossO(j)=mlossO(j)+numII*27.66d0*exp(-51.81d0/meanmassM)*scale_m/vol_loc(j)    !------ Wosley & Heger (2007)
        
                 ! --- Reduce star particle mass
                 mp(ind_part(j))=mp(ind_part(j))-mejecta
                 if(mp(ind_part(j)).le.0.0)then
                    write(*,*) "Error, mp=0 from ejecta, type II"              !------ Oscar: at mp< few 100 Msun, stochastic SN can trigger this. To be solved...
                 endif
              endif
           endif
           
!---------------------------
!--------------------------- Type Ia
!---------------------------
         if(mstarmin.gt.MBmin/2.0.and.mstarmax.lt.MBmax/2.0) then  !----- Mass limits in solarmasses (see Raiteri et al. 1996)
         call SNIanum(mstarmax,mstarmin,theint)  
         NumSNIa=theint*mp0(ind_part(j))/scale_m

         ! Do random sampling of type Ia SNe
         call ranf(localseed,RandNum)
         numresidual=NumSNIa-int(NumSNIa)              !----- int <1 --> 0
         NumSNIa=int(NumSNIa)
         if(RandNum<numresidual) then
            NumSNIa=NumSNIa+1      
         endif

         if(NumSNIa>0.0) then         
            ! -- SNIA mass loss
            masslossIa=NumSNIa*Mremnant 
            mloss(j)=mloss(j)+masslossIa/vol_loc(j)

            ! ------ Momentum, same approch as for typeII
            if(numII<=0.)then                           !---------- if cell is already check for ST momentum in Type II event, we can save so time
               Zscale=(max(Zgas,0.01d0))**(-0.2)        !---------- Thornton et al. 1998. Zgas computed above assuming solar mixture
               n0=max(unew(indp(j),1),smallr)*scale_nH 
               pST=4.0d5*NumSNIa**(0.941)*n0**(-0.1176)*Zscale 
               pST=pST*1.0d5*scale_m/scale_v      
               rsf=30.*NumSNIa**0.29*(n0**(-0.43))*((Zgas+0.01)**(-0.18))
               cellsize=dx_loc2(j)*scale_l/3.08d18    
               if(cellsize.lt.rsf/Nrcool)then 
                  momST=.false.
               else
                  momST=.true.
               endif
            endif

            if(momST) then
               ptot(j)=ptot(j)+pST  !ST momentum
            else
               ptot(j)=ptot(j)+p10*NumSNIa  !Conservative initial momentum
            endif

            ! -- SNIA thermal energy
            ethermal(j)=ethermal(j)+NumSNIa*ESNIa/vol_loc(j)  
            
            if(metal)then
               mlossFe(j)=mlossFe(j)+NumSNIa*SNIaFe/vol_loc(j)
               mlossO(j)=mlossO(j)+NumSNIa*SNIaO/vol_loc(j)
            endif
         
            ! -- Reduce star particle mass
            mp(ind_part(j))=mp(ind_part(j))-masslossIa
            if(mp(ind_part(j)).le.0.0)then
               write(*,*) "Error, mp=0 from ejecta, type Ia"              !------ Oscar: at mp< few 100 Msun, stochastic SN can trigger this. To be solved...
            endif
         endif
      endif
   endif

!---------------------------
!--------------------------- Winds
!---------------------------      
   if(winds) then 
!---------------------------
!--------------------------- High mass stars, fast winds
!---------------------------
        if(t2.ge.0.0.and.t1.le.twind) then
           call fm_w(t1,t2,mett,theint)                                            !------- Get fraction os stellar mass lost in winds
           mloss(j)=mloss(j)+mp0(ind_part(j))*theint/vol_loc(j)  
           mlossFe(j)=mlossFe(j)+zp(ind_part(j))*mp0(ind_part(j))*theint/vol_loc(j)
           mlossO(j)=mlossO(j)+zp2(ind_part(j))*mp0(ind_part(j))*theint/vol_loc(j)  !No proper yields, just metals from stars (to be included in the future)

           ! -- Reduce star particle mass                                                                               
           mp(ind_part(j))=mp(ind_part(j))-mp0(ind_part(j))*theint
           if(mp(ind_part(j)).le.0.0)then
               write(*,*) "Error, mp=0 from ejecta, fast winds"              !------ Oscar: at mp< few 100 Msun, stochastic SN can trigger this. To be solved...
           endif
            
           call p_w(t1,t2,mett,theint)                                       !------ Get specific momentum
           ptot(j)=ptot(j)+mp0(ind_part(j))*theint/scale_v                   !------ scale by vol_loc later
           
           call E_w(t1,t2,mett,theint)                                       !------ Get energy
           ethermal(j)=ethermal(j)+mp0(ind_part(j))*theint/vol_loc(j)/scale_v/scale_v
        endif
!---------------------------
!--------------------------- Low mass stars (Kalirai et al. 2008)
!---------------------------
        if(mstarmin<=Mwindmax.and.mstarmax>=Mwindmin.and.mstarmin<=mstarmax) then     
           call AGBmassloss(mstarmin,mstarmax,theint)                        !------ Integration limits are from m_i to m_(i-1)
           masslossW=theint*mp0(ind_part(j))                     
           mloss(j)=mloss(j)+masslossW/vol_loc(j)
           if(metal)then
              metalsW=zp(ind_part(j))                                        !------ Same as host
              mlossFe(j)=mlossFe(j)+masslossW*metalsW/vol_loc(j)
              metalsW=zp2(ind_part(j))                                       !------ Same as host
              mlossO(j)=mlossO(j)+masslossW*metalsW/vol_loc(j)
           endif
           ! Reduce star particle mass
           mp(ind_part(j))=mp(ind_part(j))-masslossW
            if(mp(ind_part(j)).le.0.0)then
               write(*,*) "Error, mp=0 from ejecta, slow winds"              !------ Oscar: at mp< few 100 Msun, stochastic SN can trigger this. To be solved...
            endif
        endif
     endif
!---------------------------
!--------------------------- Radiation pressure, jot to be used when full RT is used
!---------------------------
     if(radpressure) then
        if(t2.gt.0.0.and.t1.lt.tcl) then                      ! -- Bin young star particles get "cluster mass" in cell. Use tcl~1-10Myr
           if(t2<=tcut) then                                  ! -- t2 is current age
              Lum=L1*mp0(ind_part(j))       
           else        
              Lum=(L1*(t2/tcut)**(-1.25d0))*mp0(ind_part(j)) 
           endif
           if(t2>40.0d6) then !No more
              Lum=0.0
           endif
           if(ind(indp(j)).gt.0) then                           ! -- Already a particle at location
              icenter=ind(indp(j)) 
              mcl(icenter)=mcl(icenter)+mp0(ind_part(j))
              agecl(icenter)=agecl(icenter)+t2*mp0(ind_part(j)) ! -- get average age of stars in cluster
              Lumcl(icenter)=Lumcl(icenter)+Lum*(t2-t1)*365.*24.*3600./scale_t 
           else
              irad=irad+1
              ind(indp(j))=irad !Associate entry (cell index) to Prad entry
              mcl(irad)=mcl(irad)+mp0(ind_part(j))
              agecl(irad)=agecl(irad)+t2*mp0(ind_part(j)) !to get average age of stars in cluster                               
              Lumcl(irad)=Lumcl(irad)+Lum*(t2-t1)*365.*24.*3600./scale_t
              indrad(irad)=j !We need book-keeping to get back to indcube.
           endif
        endif
     endif
  enddo

!---------------------------
!--------------------------- Radiation pressure magnitude (abbove step is finding mass in young stars)
!---------------------------
  iradmax=irad
  if(radpressure) then
     if(iradmax.gt.0) then
        do i=1,iradmax                       ! -- over cells now
           if(Lumcl(i).gt.0.0) then          ! -- only do for cells with actual young stars in them
              agecl(i)=agecl(i)/mcl(i)       ! -- Normalize cluster age
              iicell=indp(indrad(i))         ! -- get particle "j"
              
              if(metalscaling) then          ! -- for dust opacity
                 Zgas=(2.09d0*unew(iicell,imetal+1)+1.06d0*unew(iicell,imetal))/&
                      & max(unew(iicell,1),smallr)/0.02d0 !Average, solar mix in gas (Asplund)
                 Zgas=max(Zgas,0.01) 
              else
                 Zgas=1.0
              endif
              Mclmax=mumax*mcl(i)  
              KappaIR=KappaIR_0*Zgas !Scaled by Z/Z_sun to get dust-to-gas ratio dependency. Current SN ejecta is included
              if(mcl(i)<=mtrans) then 
                 Cr=Cr1        
                 alpha=alpha1     
              else        
                 Cr=Cr2    
                 alpha=alpha2     
              endif
              if(tau_IR.ge.0) then
                 tauIR=tau_IR !set in paramter file
              else
                 tauIR=KappaIR*((1.-eps_cl)*(2.-beta)/(2.*acos(-1.)*Cr**2)/(3.-beta-2.*alpha))
                 tauIR=tauIR*(mumax/eps_cl)**(1.-2.*alpha)*(1.-(Mclmin/Mclmax)**(3.-2.*alpha-beta))
                 tauIR=tauIR/(1.-(Mclmin/Mclmax)**(2.-beta))
                 tauIR=tauIR*mcl(i)**(1.-2.*alpha)   !Correct, as we multiply by mp below
              endif
              
              tauIR=min(tauIR,50.0d0)  !Prad limiter

              if(agecl(i).lt.tcl)then !If clump is still intact
                 Prad(i)=(eta1+eta2*tauIR)*Lumcl(i)  !Lumcl=L1*mcl, we don't need mcl here!! dteff is accounted for above!
              else
                 Prad(i)=(eta1+eta2*KappaIR*(unew(iicell,1))*dx_loc)*Lumcl(i)  !Lumcl=L1*mcl*dteff
              endif
              
              if(ysc_stats) then  ! Get YSC stats printed to screen (change to file)
                 if(agecl(i).lt.10.0d6)then 
                    clustermass=mcl(i)*scale_d*scale_l*scale_l*scale_l/2.0d33 !Young cluster mass
                    clusterage=agecl(i)/1.0d6  !In Myr
                    write(*,"(a,6e14.5,i9)") "YSO:",t,clusterage,clustermass,unew(iicell,1)*scale_nH,tauIR,Zgas,ilevel
                 endif
              endif
           endif
        enddo
     endif
  endif
  
  !----------- Inject feedback ----------------
  
  do j=1,np
     if(mloss(j)>0.or.ethermal(j)>0.or.ptot(j)>0.) then  ! -- only enter if star actually injects something
        if(ok(j)) then                                 ! -- for particles in oct
        do ii=1,2                                      ! -- Do feedback over 3x3 cube
           do jj=1,2
              do kk=1,2
                 iicell=indcube2(j,ii,jj,kk)
                 
                 if(iicell.gt.0) then
                    !-------------- Do momentum and mass-loss under assumption of no Etherm increase
                    tt=unew(iicell,ndim+2)  
                    tekin=0.0d0
                    do idim=1,ndim
                       tekin=tekin+0.5*unew(iicell,idim+1)**2/max(unew(iicell,1),smallr) ! -- Kinetic E
                    end do
                    tt=tt-tekin  !Etherm                                                                  
                    pmiss=0.0
                    
                    !------- Geomtrical factors and cell index ---------------------------------------
                    xcont=-1.0+2.0*(ii-1)
                    ycont=-1.0+2.0*(jj-1)
                    zcont=-1.0+2.0*(kk-1)
                    contr=8.0*(xcont**2+ycont**2+zcont**2)**0.5  ! -- Each cell get 1/8th of momentum
                    
                    !----- Return ejected mass -------------------------------------------------------

                    if(mloss(j)>0.) then
                       unew(iicell,1)=unew(iicell,1)+(mloss(j))/8.0                     ! -- Spread over 8 cells
                       unew(iicell,2)=unew(iicell,2)+(mloss(j))*vp(ind_part(j),1)/8.0   ! -- Ejected mass momentum due to stars motion
                       unew(iicell,3)=unew(iicell,3)+(mloss(j))*vp(ind_part(j),2)/8.0
                       unew(iicell,4)=unew(iicell,4)+(mloss(j))*vp(ind_part(j),3)/8.0
                    endif
                    
                    if(momentum) then
                       if(ptot(j)>0.) then
                       !------ Momentum from winds and SNe -----------------
                       vkick=scale_v*ptot(j)/8.d0/1.d5/max(unew(iicell,1),smallr)/vol_loc(j) !Velocity for exactly ptot/8/mas\
                       
                       vxold=unew(iicell,2)/max(unew(iicell,1),smallr) !extra safety
                       vyold=unew(iicell,3)/max(unew(iicell,1),smallr)
                       vzold=unew(iicell,4)/max(unew(iicell,1),smallr)
                       
                       if(vkick.gt.vmax) then  !Limit momentum in this way (for stabbility)
                          momx=xcont*8.*unew(iicell,1)*vmax*1.d5/scale_v/contr  !mom density, factor of 8 is to cancel contr. Vmax is the correct velocity
                          momy=ycont*8.*unew(iicell,1)*vmax*1.d5/scale_v/contr  !mom density                       
                          momz=zcont*8.*unew(iicell,1)*vmax*1.d5/scale_v/contr  !mom density                   
                          ! missing -------
                          momtot=(momx**2+momy**2+momz**2)**0.5
                          momxt=xcont*ptot(j)/contr/vol_loc(j)
                          momyt=ycont*ptot(j)/contr/vol_loc(j)
                          momzt=zcont*ptot(j)/contr/vol_loc(j)
                          momtott=(momxt**2+momyt**2+momzt**2)**0.5
                          pmiss=pmiss+(momtott-momtot)
                       else
                          momx=xcont*ptot(j)/contr/vol_loc(j)
                          momy=ycont*ptot(j)/contr/vol_loc(j)
                          momz=zcont*ptot(j)/contr/vol_loc(j)
                       endif
                       unew(iicell,2)=unew(iicell,2)+momx 
                       unew(iicell,3)=unew(iicell,3)+momy 
                       unew(iicell,4)=unew(iicell,4)+momz 
                    endif
                    ! ------------------------------
                    if(fbsafety) then                       
                       vxnew=unew(iicell,2)/max(unew(iicell,1),smallr)  !extra safety
                       vynew=unew(iicell,3)/max(unew(iicell,1),smallr)
                       vznew=unew(iicell,4)/max(unew(iicell,1),smallr)

                       if(abs(vxnew).gt.maxadv) then
                          unew(iicell,2)=sign(maxadv,vxnew)*unew(iicell,1) 
                       endif
                       if(abs(vynew).gt.maxadv) then
                          unew(iicell,3)=sign(maxadv,vynew)*unew(iicell,1)
                       endif
                       if(abs(vznew).gt.maxadv) then
                          unew(iicell,4)=sign(maxadv,vznew)*unew(iicell,1)
                       endif
                    endif
           
         ! ----- All momentum is now added, calculate new Ekin and update Etot. This is consistent with SNe mass loading
                  tekin=0.0d0
                  do idim=1,ndim
                     tekin=tekin+0.5*unew(iicell,idim+1)**2/max(unew(iicell,1),smallr)            ! -- Kin E 
                  enddo
                  unew(iicell,ndim+2)=tt+tekin    
               endif
                  !tt is old Etherm, which should not change due to momentum additions
               if(energy) then
                  if(ethermal(j)>0.) then
                     !------ Update thermal energy -----------------------------------------------------
                     if(iicell==indp(j))then  !"central cell"
                        
                        Emax=tekin+Tmaxfb*unew(iicell,1)/scale_T2/(gamma-1.0d0)
                        unew(iicell,ndim+2)=unew(iicell,ndim+2)+ethermal(j)
                                              
                 !       if(pnonthermal) then  !adding missed momentum to cell via pressure
                 !          unew(iicell,ndim+2)=unew(iicell,ndim+2)+pmiss/6.0/dx_loc2(j)/dteff(j)
                 !       endif
                        
                        unew(iicell,ndim+2)=min(unew(iicell,ndim+2),Emax)  !Oscar, safety
                        
                     endif
                  endif
               endif
                  ! --------- Finally Return metals -----
               if(mlossFe(j)>0.or.mlossO(j)>0) then
                  unew(iicell,imetal)=unew(iicell,imetal)+mlossFe(j)/8.0   !8 cell stencil, neighbours must have same volume (oct)
                  unew(iicell,imetal+1)=unew(iicell,imetal+1)+mlossO(j)/8.0   !8 cell stencil, neighbours must have same volume (oct)
               endif
            endif
            if(iicell.le.0) then
               write(*,*) "Error in oct, iicell=0",iicell
            endif
         enddo
      enddo
   enddo
else   ! if not OK ---------------------------------------  DRIFTERS
      
      iicell=indp(j)
      
      !-------------- Do momentum and mass-loss under assumption of no Etherm increase
      tt=unew(iicell,ndim+2)  !Tot E 
      tekin=0.0d0
      do idim=1,ndim
         !         tekin=tekin+0.5*unew(iicell,idim+1)**2/max(unew(iicell,1),smallr)  !Kin E                       
         tekin=tekin+0.5*unew(iicell,idim+1)**2/max(unew(iicell,1),smallr)  !Kin E
      end do
      tt=tt-tekin  !Etherm                                                                  
      pmiss=0.0
      !------- Geomtrical factors and cell index ---------------------------------------------
      xcont=-1.0+2.0*(ii-1) !rand this
      ycont=-1.0+2.0*(jj-1)
      zcont=-1.0+2.0*(kk-1)
      contr=(xcont**2+ycont**2+zcont**2)**0.5  !Each cell get 1/8th of momentum
      
      !----- Return ejected mass ------------------------------------------------------------- 
      if(mloss(j)>0.) then      
         unew(iicell,1)=unew(iicell,1)+(mloss(j))  !Spread over 1 cells
         unew(iicell,2)=unew(iicell,2)+(mloss(j))*vp(ind_part(j),1)  !-------- Ejected mass momentum due to stars motion
         unew(iicell,3)=unew(iicell,3)+(mloss(j))*vp(ind_part(j),2)
         unew(iicell,4)=unew(iicell,4)+(mloss(j))*vp(ind_part(j),3)
      endif
        !------ Momentum
        
        ! Future: add momentum as random kick
        
        ! ----- All momentum is now added, calculate new Ekin and update Etot. This is consistent with SNe mass loading
        tekin=0.0d0
        do idim=1,ndim
           !           tekin=tekin+0.5*unew(iicell,idim+1)**2/max(unew(iicell,1),smallr)  !Kin E                          
           tekin=tekin+0.5*unew(iicell,idim+1)**2/max(unew(iicell,1),smallr)  !Kin E                          \
           
        enddo
        unew(iicell,ndim+2)=tt+tekin    !tt is old Etherm, which should not change due to momentum additions
        
        if(energy) then
           !------ Update thermal energy -----------------------------------------------------
           
           Emax=tekin+Tmaxfb*unew(iicell,1)/scale_T2/(gamma-1.0d0)
           unew(iicell,ndim+2)=unew(iicell,ndim+2)+ethermal(j)
          
          !if(pnonthermal) then !momentum done in nonthermal way
          !   unew(iicell,ndim+2)=unew(iicell,ndim+2)+ptot(j)/6.0/dx_loc2(j)/dteff(j)
          !endif
          
          unew(iicell,ndim+2)=min(unew(iicell,ndim+2),Emax)  !Oscar, safety
       endif
         
       ! --------- Finally Return metals -----
       if(mlossFe(j)>0.or.mlossO(j)>0) then
          unew(iicell,imetal)=unew(iicell,imetal)+mlossFe(j)      !1 cell for drifters
          unew(iicell,imetal+1)=unew(iicell,imetal+1)+mlossO(j)   !1 cell for drifters
       endif
    endif
 endif
enddo
 
 !-------------------------- Radiation pressure ------------------
 if(radpressure) then !Put this with the rest later
    if(iradmax.gt.0) then 
       do i=1,iradmax                   
          if(Prad(i)>0) then
          !We now the number of entries in Prad array, not particles    
          icenter=indrad(i)                          ! -- Particle index   
          if(ok(icenter)) then                       ! -- Only do for particles on local oct
             do ii=1,2  !j is now picked
                do jj=1,2
                   do kk=1,2
                      iicell=indcube2(icenter,ii,jj,kk) 
                      
                      ! -- Kicks are applied assuming constant etherm
                      tt=unew(iicell,ndim+2)  !Tot E                                                                           
                      tekin=0.0d0
                      do idim=1,ndim
                         tekin=tekin+0.5*unew(iicell,idim+1)**2/max(unew(iicell,1),smallr)  !Kin E
                      end do
                      tt=tt-tekin  !Etherm     
                      
                      !------- Geomtrical factors and cell index ---------------------------------------------
                      xcont=-1.0+2.0*(ii-1)
                      ycont=-1.0+2.0*(jj-1)
                      zcont=-1.0+2.0*(kk-1)
                      contr=8.0*(xcont**2+ycont**2+zcont**2)**0.5  !Each cell get 1/8th of momentum
                      
                      !               vkick=scale_v*Prad(i)/8.d0/1.d5/max(unew(iicell,1),smallr)/vol_loc(icenter)
                      
                      vkick=scale_v*Prad(i)/8.d0/1.d5/max(unew(iicell,1),smallr)/vol_loc(icenter)
                      
                      if(vkick.gt.vmax) then  !Limit momentum
                         momx=xcont*8.*unew(iicell,1)*vmax*1.d5/scale_v/contr  !mom density, factor of 8 is to cancel contr
                         momy=ycont*8.*unew(iicell,1)*vmax*1.d5/scale_v/contr  !mom density                       
                         momz=zcont*8.*unew(iicell,1)*vmax*1.d5/scale_v/contr  !mom density                   
                      else
                         momx=xcont*Prad(i)/contr/vol_loc(icenter)
                         momy=ycont*Prad(i)/contr/vol_loc(icenter)
                         momz=zcont*Prad(i)/contr/vol_loc(icenter)
                      endif
                      unew(iicell,2)=unew(iicell,2)+momx 
                      unew(iicell,3)=unew(iicell,3)+momy 
                      unew(iicell,4)=unew(iicell,4)+momz 
                      
                      if(fbsafety) then
           
                         vxnew=unew(iicell,2)/max(unew(iicell,1),smallr)
                         vynew=unew(iicell,3)/max(unew(iicell,1),smallr)
                         vznew=unew(iicell,4)/max(unew(iicell,1),smallr)
                         
                         if(abs(vxnew).gt.maxadv) then
                            unew(iicell,2)=sign(maxadv,vxnew)*unew(iicell,1)
                         endif
                         if(abs(vynew).gt.maxadv) then
                            unew(iicell,3)=sign(maxadv,vynew)*unew(iicell,1)
                         endif
                         if(abs(vznew).gt.maxadv) then
                            unew(iicell,4)=sign(maxadv,vznew)*unew(iicell,1)
                         endif
                      endif
                      
                      ! ----- All momentum is now added, calculate new Ekin and update Etot
                      tekin=0.0d0
                      do idim=1,ndim
                         tekin=tekin+0.5*unew(iicell,idim+1)**2/max(unew(iicell,1),smallr)  !Kin E     
                      enddo
                      unew(iicell,ndim+2)=tt+tekin  
                      !--------------------------------- Oscar: add pnonthemal for missing momentum
                   enddo
                enddo
             enddo
          else
             iicell=indp(icenter)         
             ! -- Kicks are applied assuming constant etherm
             if(pnonthermal) then  !adding missed momentum to cell via pressure
                tt=unew(iicell,ndim+2)  !Tot E
                tekin=0.0d0
                do idim=1,ndim
                   tekin=tekin+0.5*unew(iicell,idim+1)**2/max(unew(iicell,1),smallr) !Kin E
                end do
                tt=tt-tekin  !Etherm + Erad
                
                ! nonthermal
                unew(iicell,ndim+2)=unew(iicell,ndim+2)+Prad(i)/6.0/dx_loc2(icenter)/dteff(icenter)
                
                ! -- limiter
                Emax=tekin+Tmaxfb*unew(iicell,1)/scale_T2/(gamma-1.0d0)
                unew(iicell,ndim+2)=min(unew(iicell,ndim+2),Emax)  !Oscar, safety 
                
             endif
          endif
       endif
    enddo
 endif
endif
 
deallocate(ind) 

! 777 continue !FBTEST
 
 ! Romain's
  ! Compute stellar mass loss and thermal feedback due to supernovae
!  if(f_w==0)then
!     do j=1,np
!        birth_time=tp(ind_part(j))
!        ! Make sure that we don't count feedback tiwce
!        if(birth_time.lt.(current_time-t0).and.birth_time.ge.(current_time-t0-dteff(j)))then !          
!           ! Stellar mass loss
!           mejecta=eta_sn*mp(ind_part(j))
!           mloss(j)=mloss(j)+mejecta/vol_loc(j)
!           ! Thermal energy
!           ethermal(j)=ethermal(j)+mejecta*ESN/vol_loc(j)
!           ! Metallicity
!           if(metal)then
!              zloss=yield+(1d0-yield)*zp(ind_part(j))
!              mzloss(j)=mzloss(j)+mejecta*zloss/vol_loc(j)
!           endif
!           ! Reduce star particle mass
!           mp(ind_part(j))=mp(ind_part(j))-mejecta
!           ! Boost SNII energy and depopulate accordingly
!           if(SN_BOOST>1d0)then
!              call ranf(localseed,RandNum)
!              if(RandNum<1d0/SN_BOOST)then
!                 mloss(j)=SN_BOOST*mloss(j)
!                 mzloss(j)=SN_BOOST*mzloss(j)
!                 ethermal(j)=SN_BOOST*ethermal(j)
!              else
!                 mloss(j)=0d0
!                 mzloss(j)=0d0
!                 ethermal(j)=0d0
!              endif
!           endif           
!        endif
!     end do
!  endif
  
end subroutine feedbk
#endif
!################################################################
!################################################################
!################################################################
!################################################################
subroutine kinetic_feedback
  use amr_commons
  use pm_commons
  use hydro_commons
  use constants, only:Myr2sec
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info
  integer,dimension(1:ncpu)::nSN_icpu_all
  real(dp),dimension(:),allocatable::mSN_all,sSN_all,ZSN_all
  real(dp),dimension(:,:),allocatable::xSN_all,vSN_all
#endif
  !----------------------------------------------------------------------
  ! This subroutine compute the kinetic feedback due to SNII and
  ! imolement this using exploding GMC particles.
  ! This routine is called only at coarse time step.
  ! Yohan Dubois
  !----------------------------------------------------------------------
  ! local constants
  integer::ip,icpu,igrid,jgrid,npart1,npart2,ipart,jpart,next_part
  integer::nSN,nSN_loc,nSN_tot,iSN,ilevel,ivar
  integer,dimension(1:ncpu)::nSN_icpu
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,t0
  real(dp)::current_time
  real(dp)::scale,dx_min,vol_min,mstar
  integer::nx_loc
  integer,dimension(:),allocatable::ind_part,ind_grid
  logical,dimension(:),allocatable::ok_free
  integer,dimension(:),allocatable::indSN
  real(dp),dimension(:),allocatable::mSN,sSN,ZSN,m_gas,vol_gas,ekBlast
  real(dp),dimension(:,:),allocatable::xSN,vSN,u_gas,dq

  if(.not. hydro)return
  if(ndim.ne.3)return

  if(verbose)write(*,*)'Entering make_sn'

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_min=(0.5D0**nlevelmax)*scale
  vol_min=dx_min**ndim

  ! Minimum star particle mass
  if(m_star < 0d0)then
     mstar=n_star/(scale_nH*aexp**3)*vol_min
  else
     mstar=m_star*mass_sph
  endif


  ! Lifetime of Giant Molecular Clouds from Myr to code units
  ! Massive star lifetime from Myr to code units
  if(use_proper_time)then
     t0=t_sne*Myr2sec/(scale_t/aexp**2)
     current_time=texp
  else
     t0=t_sne*Myr2sec/scale_t
     current_time=t
  endif

  !------------------------------------------------------
  ! Gather GMC particles eligible for disruption
  !------------------------------------------------------
  nSN_loc=0
  ! Loop over levels
  do icpu=1,ncpu
  ! Loop over cpus
     igrid=headl(icpu,levelmin)
     ! Loop over grids
     do jgrid=1,numbl(icpu,levelmin)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0
        ! Count old enough GMC particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if ( is_debris(typep(ipart)) .and. tp(ipart).lt.(current_time-t0) ) then
                 npart2=npart2+1
              endif
              ipart=next_part  ! Go to next particle
            end do
        endif
        nSN_loc=nSN_loc+npart2   ! Add SNe to the total
        igrid=next(igrid)   ! Go to next grid
     end do
  end do
  ! End loop over levels
  nSN_icpu=0
  nSN_icpu(myid)=nSN_loc
#ifndef WITHOUTMPI
  ! Give an array of number of SN on each cpu available to all cpus
  call MPI_ALLREDUCE(nSN_icpu,nSN_icpu_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  nSN_icpu=nSN_icpu_all
#endif

  nSN_tot=sum(nSN_icpu(1:ncpu))

  if (nSN_tot .eq. 0) return

  if(myid==1)then
     write(*,*)'-----------------------------------------------'
     write(*,*)'Number of GMC to explode=',nSN_tot
     write(*,*)'-----------------------------------------------'
  endif

  ! Allocate arrays for the position and the mass of the SN
  allocate(xSN(1:nSN_tot,1:3),vSN(1:nSN_tot,1:3))
  allocate(mSN(1:nSN_tot),sSN(1:nSN_tot),ZSN(1:nSN_tot))
  xSN=0; vSN=0; mSN=0; sSN=0; ZSN=0
  ! Allocate arrays for particles index and parent grid
  if(nSN_loc>0)then
     allocate(ind_part(1:nSN_loc),ind_grid(1:nSN_loc),ok_free(1:nSN_loc))
  endif

  !------------------------------------------------------
  ! Store position and mass of the GMC into the SN array
  !------------------------------------------------------
  if(myid==1)then
     iSN=0
  else
     iSN=sum(nSN_icpu(1:myid-1))
  endif
  ! Loop over levels
  ip=0
  do icpu=1,ncpu
     igrid=headl(icpu,levelmin)
     ! Loop over grids
     do jgrid=1,numbl(icpu,levelmin)
        npart1=numbp(igrid)  ! Number of particles in the grid
        ! Count old enough star particles that have not exploded
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if ( is_debris(typep(ipart)) .and. tp(ipart).lt.(current_time-t0) ) then
                 iSN=iSN+1
                 xSN(iSN,1)=xp(ipart,1)
                 xSN(iSN,2)=xp(ipart,2)
                 xSN(iSN,3)=xp(ipart,3)
                 vSN(iSN,1)=vp(ipart,1)
                 vSN(iSN,2)=vp(ipart,2)
                 vSN(iSN,3)=vp(ipart,3)
                 mSN(iSN)=mp(ipart)
                 sSN(iSN)=dble(-idp(ipart))*mstar
                 if(metal)ZSN(iSN)=zp(ipart)
                 ip=ip+1
                 ind_grid(ip)=igrid
                 ind_part(ip)=ipart
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif
        igrid=next(igrid)   ! Go to next grid
     end do
  end do
  ! End loop over levels

  ! Remove GMC particle
  if(nSN_loc>0)then
     ok_free=.true.
     call remove_list(ind_part,ind_grid,ok_free,nSN_loc)
     call add_free_cond(ind_part,ok_free,nSN_loc)
     deallocate(ind_part,ind_grid,ok_free)
  endif

#ifndef WITHOUTMPI
  allocate(xSN_all(1:nSN_tot,1:3),vSN_all(1:nSN_tot,1:3),mSN_all(1:nSN_tot),sSN_all(1:nSN_tot),ZSN_all(1:nSN_tot))
  call MPI_ALLREDUCE(xSN,xSN_all,nSN_tot*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(vSN,vSN_all,nSN_tot*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(mSN,mSN_all,nSN_tot  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(sSN,sSN_all,nSN_tot  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(ZSN,ZSN_all,nSN_tot  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  xSN=xSN_all
  vSN=vSN_all
  mSN=mSN_all
  sSN=sSN_all
  ZSN=ZSN_all
  deallocate(xSN_all,vSN_all,mSN_all,sSN_all,ZSN_all)
#endif

  nSN=nSN_tot
  allocate(m_gas(1:nSN),u_gas(1:nSN,1:3),vol_gas(1:nSN),dq(1:nSN,1:3),ekBlast(1:nSN))
  allocate(indSN(1:nSN))

  ! Compute the grid discretization effects
  call average_SN(xSN,vol_gas,dq,ekBlast,indSN,nSN)

  ! Modify hydro quantities to account for a Sedov blast wave
  call Sedov_blast(xSN,vSN,mSN,sSN,ZSN,indSN,vol_gas,dq,ekBlast,nSN)

  deallocate(xSN,vSN,mSN,sSN,ZSN,indSN,m_gas,u_gas,vol_gas,dq,ekBlast)

  ! Update hydro quantities for split cells
  do ilevel=nlevelmax,levelmin,-1
     call upload_fine(ilevel)
     do ivar=1,nvar
        call make_virtual_fine_dp(uold(1,ivar),ilevel)
     enddo
  enddo

end subroutine kinetic_feedback
!################################################################
!################################################################
!################################################################
!################################################################
subroutine average_SN(xSN,vol_gas,dq,ekBlast,ind_blast,nSN)
  use pm_commons
  use amr_commons
  use hydro_commons
  use constants, only: pc2cm
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
  !------------------------------------------------------------------------
  ! This routine average the hydro quantities inside the SN bubble
  !------------------------------------------------------------------------
  integer::ilevel,ncache,nSN,iSN,ind,ix,iy,iz,ngrid,iskip
  integer::i,nx_loc,igrid
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp)::x,y,z,dr_SN,u,v,w,u2,v2,w2,dr_cell
  real(dp)::scale,dx,dxx,dyy,dzz,dx_min,dx_loc,vol_loc,rmax2,rmax
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  integer ,dimension(1:nSN)::ind_blast
  real(dp),dimension(1:nSN)::vol_gas,ekBlast
  real(dp),dimension(1:nSN,1:3)::xSN,dq,u2Blast
#ifndef WITHOUTMPI
  real(dp),dimension(1:nSN)::vol_gas_all,ekBlast_all
  real(dp),dimension(1:nSN,1:3)::dq_all,u2Blast_all
#endif
  logical ,dimension(1:nvector),save::ok

  if(nSN==0)return
  if(verbose)write(*,*)'Entering average_SN'

  ! Mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Maximum radius of the ejecta
  rmax=MAX(2.0d0*dx_min*scale_l/aexp,rbubble*pc2cm)
  rmax=rmax/scale_l
  rmax2=rmax*rmax

  ! Initialize the averaged variables
  vol_gas=0; dq=0; u2Blast=0; ekBlast=0; ind_blast=-1

  ! Loop over levels
  do ilevel=levelmin,nlevelmax
     ! Computing local volume (important for averaging hydro quantities)
     dx=0.5D0**ilevel
     dx_loc=dx*scale
     vol_loc=dx_loc**ndim
     ! Cells center position relative to grid center position
     do ind=1,twotondim
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc(ind,1)=(dble(ix)-0.5D0)*dx
        xc(ind,2)=(dble(iy)-0.5D0)*dx
        xc(ind,3)=(dble(iz)-0.5D0)*dx
     end do

     ! Loop over grids
     ncache=active(ilevel)%ngrid
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

           ! Flag leaf cells
           do i=1,ngrid
              ok(i)=son(ind_cell(i))==0
           end do

           do i=1,ngrid
              if(ok(i))then
                 ! Get gas cell position
                 x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
                 y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
                 z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
                 do iSN=1,nSN
                    ! Check if the cell lies within the SN radius
                    dxx=x-xSN(iSN,1)
                    dyy=y-xSN(iSN,2)
                    dzz=z-xSN(iSN,3)
                    dr_SN=dxx**2+dyy**2+dzz**2
                    dr_cell=MAX(ABS(dxx),ABS(dyy),ABS(dzz))
                    if(dr_SN.lt.rmax2)then
                       vol_gas(iSN)=vol_gas(iSN)+vol_loc
                       ! Take account for grid effects on the conservation of the
                       ! normalized linear momentum
                       u=dxx/rmax
                       v=dyy/rmax
                       w=dzz/rmax
                       ! Add the local normalized linear momentum to the total linear
                       ! momentum of the blast wave (should be zero with no grid effect)
                       dq(iSN,1)=dq(iSN,1)+u*vol_loc
                       dq(iSN,2)=dq(iSN,2)+v*vol_loc
                       dq(iSN,3)=dq(iSN,3)+w*vol_loc
                       u2Blast(iSN,1)=u2Blast(iSN,1)+u*u*vol_loc
                       u2Blast(iSN,2)=u2Blast(iSN,2)+v*v*vol_loc
                       u2Blast(iSN,3)=u2Blast(iSN,3)+w*w*vol_loc
                    endif
                    if(dr_cell.le.dx_loc/2.0)then
                       ind_blast(iSN)=ind_cell(i)
                       ekBlast  (iSN)=vol_loc
                    endif
                 end do
              endif
           end do

        end do
        ! End loop over cells
     end do
     ! End loop over grids
  end do
  ! End loop over levels

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(vol_gas,vol_gas_all,nSN  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(dq     ,dq_all     ,nSN*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(u2Blast,u2Blast_all,nSN*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(ekBlast,ekBlast_all,nSN  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  vol_gas=vol_gas_all
  dq     =dq_all
  u2Blast=u2Blast_all
  ekBlast=ekBlast_all
#endif
  do iSN=1,nSN
     if(vol_gas(iSN)>0d0)then
        dq(iSN,1)=dq(iSN,1)/vol_gas(iSN)
        dq(iSN,2)=dq(iSN,2)/vol_gas(iSN)
        dq(iSN,3)=dq(iSN,3)/vol_gas(iSN)
        u2Blast(iSN,1)=u2Blast(iSN,1)/vol_gas(iSN)
        u2Blast(iSN,2)=u2Blast(iSN,2)/vol_gas(iSN)
        u2Blast(iSN,3)=u2Blast(iSN,3)/vol_gas(iSN)
        u2=u2Blast(iSN,1)-dq(iSN,1)**2
        v2=u2Blast(iSN,2)-dq(iSN,2)**2
        w2=u2Blast(iSN,3)-dq(iSN,3)**2
        ekBlast(iSN)=max(0.5d0*(u2+v2+w2),0.0d0)
     endif
  end do

  if(verbose)write(*,*)'Exiting average_SN'

end subroutine average_SN
!################################################################
!################################################################
!################################################################
!################################################################
subroutine Sedov_blast(xSN,vSN,mSN,sSN,ZSN,indSN,vol_gas,dq,ekBlast,nSN)
  use pm_commons
  use amr_commons
  use hydro_commons
  use constants, only: M_sun, pc2cm
  use mpi_mod
  implicit none
  !------------------------------------------------------------------------
  ! This routine merges SN using the FOF algorithm.
  !------------------------------------------------------------------------
  integer::ilevel,iSN,nSN,ind,ix,iy,iz,ngrid,iskip
  integer::i,nx_loc,igrid,ncache
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp)::x,y,z,dx,dxx,dyy,dzz,dr_SN,u,v,w,ESN,mstar,eta_sn2,msne_min,mstar_max
  real(dp)::scale,dx_min,dx_loc,vol_loc,rmax2,rmax,vol_min
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nSN)::mSN,sSN,ZSN,p_gas,d_gas,d_metal,vol_gas,uSedov,ekBlast
  real(dp),dimension(1:nSN,1:3)::xSN,vSN,dq
  integer ,dimension(1:nSN)::indSN
  logical ,dimension(1:nvector),save::ok

  if(nSN==0)return
  if(verbose)write(*,*)'Entering Sedov_blast'

  ! Mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax
  vol_min=dx_min**ndim

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Maximum radius of the ejecta
  rmax=MAX(2.0d0*dx_min*scale_l/aexp,rbubble*pc2cm)
  rmax=rmax/scale_l
  rmax2=rmax*rmax

  ! Minimum star particle mass
  if(m_star < 0d0)then
     mstar=n_star/(scale_nH*aexp**3)*vol_min
  else
     mstar=m_star*mass_sph
  endif
  msne_min=mass_sne_min*M_sun/(scale_d*scale_l**3)
  mstar_max=mass_star_max*M_sun/(scale_d*scale_l**3)
  ! Supernova specific energy from cgs to code units
  ESN=(1d51/(10d0*M_sun))/scale_v**2

  do iSN=1,nSN
     eta_sn2    = eta_sn
     if(sf_imf)then
        if(mSN(iSN).le.mstar_max)then
           if(mSN(iSN).ge.msne_min) eta_sn2 = eta_ssn
           if(mSN(iSN).lt.msne_min) eta_sn2 = 0
        endif
     endif
     if(vol_gas(iSN)>0d0)then
        d_gas(iSN)=mSN(iSN)/vol_gas(iSN)
        if(metal)d_metal(iSN)=ZSN(iSN)*mSN(iSN)/vol_gas(iSN)
        if(ekBlast(iSN)==0d0)then
           p_gas(iSN)=eta_sn2*sSN(iSN)*ESN/vol_gas(iSN)
           uSedov(iSN)=0d0
        else
           p_gas(iSN)=(1d0-f_ek)*eta_sn2*sSN(iSN)*ESN/vol_gas(iSN)
           uSedov(iSN)=sqrt(f_ek*eta_sn2*sSN(iSN)*ESN/mSN(iSN)/ekBlast(iSN))
        endif
     else
        d_gas(iSN)=mSN(iSN)/ekBlast(iSN)
        p_gas(iSN)=eta_sn2*sSN(iSN)*ESN/ekBlast(iSN)
        if(metal)d_metal(iSN)=ZSN(iSN)*mSN(iSN)/ekBlast(iSN)
     endif
  end do

  ! Loop over levels
  do ilevel=levelmin,nlevelmax
     ! Computing local volume (important for averaging hydro quantities)
     dx=0.5D0**ilevel
     dx_loc=dx*scale
     vol_loc=dx_loc**ndim
     ! Cells center position relative to grid center position
     do ind=1,twotondim
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc(ind,1)=(dble(ix)-0.5D0)*dx
        xc(ind,2)=(dble(iy)-0.5D0)*dx
        xc(ind,3)=(dble(iz)-0.5D0)*dx
     end do

     ! Loop over grids
     ncache=active(ilevel)%ngrid
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

           ! Flag leaf cells
           do i=1,ngrid
              ok(i)=son(ind_cell(i))==0
           end do

           do i=1,ngrid
              if(ok(i))then
                 ! Get gas cell position
                 x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
                 y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
                 z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
                 do iSN=1,nSN
                    ! Check if the cell lies within the SN radius
                    dxx=x-xSN(iSN,1)
                    dyy=y-xSN(iSN,2)
                    dzz=z-xSN(iSN,3)
                    dr_SN=dxx**2+dyy**2+dzz**2
                    if(dr_SN.lt.rmax2)then
                       ! Compute the mass density in the cell
                       uold(ind_cell(i),1)=uold(ind_cell(i),1)+d_gas(iSN)
                       ! Compute the metal density in the cell
                       if(metal)uold(ind_cell(i),imetal)=uold(ind_cell(i),imetal)+d_metal(iSN)
                       ! Velocity at a given dr_SN linearly interpolated between zero and uSedov
                       u=uSedov(iSN)*(dxx/rmax-dq(iSN,1))+vSN(iSN,1)
                       v=uSedov(iSN)*(dyy/rmax-dq(iSN,2))+vSN(iSN,2)
                       w=uSedov(iSN)*(dzz/rmax-dq(iSN,3))+vSN(iSN,3)
                       ! Add each momentum component of the blast wave to the gas
                       uold(ind_cell(i),2)=uold(ind_cell(i),2)+d_gas(iSN)*u
                       uold(ind_cell(i),3)=uold(ind_cell(i),3)+d_gas(iSN)*v
                       uold(ind_cell(i),4)=uold(ind_cell(i),4)+d_gas(iSN)*w
                       ! Finally update the total energy of the gas
                       uold(ind_cell(i),5)=uold(ind_cell(i),5)+0.5d0*d_gas(iSN)*(u*u+v*v+w*w)+p_gas(iSN)
                    endif
                 end do
              endif
           end do

        end do
        ! End loop over cells
     end do
     ! End loop over grids
  end do
  ! End loop over levels

  do iSN=1,nSN
     if(vol_gas(iSN)==0d0)then
        u=vSN(iSN,1)
        v=vSN(iSN,2)
        w=vSN(iSN,3)
        if(indSN(iSN)>0)then
           uold(indSN(iSN),1)=uold(indSN(iSN),1)+d_gas(iSN)
           uold(indSN(iSN),2)=uold(indSN(iSN),2)+d_gas(iSN)*u
           uold(indSN(iSN),3)=uold(indSN(iSN),3)+d_gas(iSN)*v
           uold(indSN(iSN),4)=uold(indSN(iSN),4)+d_gas(iSN)*w
           uold(indSN(iSN),5)=uold(indSN(iSN),5)+d_gas(iSN)*0.5d0*(u*u+v*v+w*w)+p_gas(iSN)
           if(metal)uold(indSN(iSN),imetal)=uold(indSN(iSN),imetal)+d_metal(iSN)
        endif
     endif
  end do

  if(verbose)write(*,*)'Exiting Sedov_blast'

end subroutine Sedov_blast
!###########################################################
!###########################################################
!###########################################################
!###########################################################
!---------------------------------------
subroutine SNIInum(m1,m2,NSNII)
!---------------------------------------
  implicit none
  REAL(kind=8),intent(out) :: NSNII
  REAL(kind=8), intent(in) :: m1,m2
  REAL(kind=8):: A,ind

 ! A=0.2244557d0  !K01
  A=0.31491d0 !Chabrier 2003
  ind=-2.3d0

  NSNII=(-A/1.3)*(m2**(-1.3)-m1**(-1.3))

END subroutine SNIInum
!###########################################################
!---------------------------------------
subroutine SNIanum(m1,m2,NSNIa)
!---------------------------------------
  implicit none
  REAL(kind=8),intent(out) :: NSNIa
  REAL(kind=8), intent(in) :: m1,m2
  REAL(kind=8):: A,Ap,N1,N2,SNIafrac

 ! A=0.2244557d0  !K01
  SNIafrac=0.16d0 !Bergh & McClure 1994, rate of SN per century in a MW type galaxy
  A=0.31491d0 !Chabrier 2003
  Ap=SNIafrac*A

  N1=(-Ap*m1**2/3.3)*((2.*m1)**(-3.3)-(m1+8.)**(-3.3))  ! (eq 14 in Agertz et al. 2013)
  N2=(-Ap*m2**2/3.3)*((2.*m2)**(-3.3)-(m2+8.)**(-3.3))

  NSNIa=(m2-m1)*(N1+N2)/2.  !Trapez. For relevant timesteps, error is epsilon

END subroutine SNIanum

!###########################################################
!---------------------------------------
subroutine AGBmassloss(m1,m2,AGB)
!---------------------------------------
  implicit none
  REAL(kind=8),intent(out) :: AGB
  REAL(kind=8), intent(in) :: m1,m2
  REAL(kind=8):: A,N1,N2

 ! A=0.2244557d0  !K01
  A=0.31491d0 !Chabrier 2003

  N1=(m1**(-0.3))*(0.3031/m1-2.97)  !Agertz et al. 2013
  N2=(m2**(-0.3))*(0.3031/m2-2.97)
  AGB=A*(N2-N1)

END subroutine AGBmassloss
!###########################################################
!------------------------------------------------
SUBROUTINE fm_w(t_1,t_2,smet,fmw)
!-----------------------------------------------
  implicit none
  real(kind=8),intent(in)::t_1,t_2,smet
  real(kind=8),intent(out)::fmw
  real(kind=8)::a,b,ts,metalscale,imfboost

  imfboost=0.3143d0/0.224468d0  !K01 to Chabrier
  
  !--- Fitting parameters --- 
  a=0.024357d0*imfboost
  b=0.000460697d0  
  ts=1.0d7
  metalscale=a*log(smet/b+1.d0)
  fmw=0.0d0

  if(t_2.le.ts) then
     fmw=metalscale*(t_2-t_1)/ts  !Linear mass loss, m(t2)-m(t1)
  endif

  if(t_1.lt.ts.and.t_2.gt.ts) then
     fmw=metalscale*(ts-t_1)/ts  !Linear mass loss, m(t2)-m(t1)
  endif

  if(t_1.ge.ts) then
     fmw=0.0
  endif

end SUBROUTINE fm_w

!###########################################################
!###########################################################
!###########################################################
!###########################################################
!------------------------------------------------
SUBROUTINE p_w(t_1,t_2,smet,momW)
!-----------------------------------------------
  implicit none
  real(kind=8),intent(in)::t_1,t_2,smet
  real(kind=8),intent(out)::momW
  real(kind=8)::a,b,c,ts,metalscale,imfboost

  imfboost=0.31430400d0/0.224468d0  !K01 to Chabrier
  !--- Fitting parameters --- 
  a=imfboost*1.8d46/1.0d6/2.d33 !Scale to per gram
  b=0.00961529d0
  c=0.363086d0
  ts=6.5d6 
  momW=0.0d0

  metalscale=a*(smet/b)**c

  if(t_2.le.ts) then
     momW=metalscale*(t_2-t_1)/ts  !Linear mass loss, m(t2)-m(t1)
  endif

  if(t_1.lt.ts.and.t_2.gt.ts) then
     momW=metalscale*(ts-t_1)/ts  !Linear mass loss, m(t2)-m(t1)
  endif

  if(t_1.ge.ts) then
     momW=0.0
  endif

end SUBROUTINE p_w

!###########################################################
!###########################################################
!###########################################################
!###########################################################
!------------------------------------------------
SUBROUTINE E_w(t_1,t_2,smet,EW)
!-----------------------------------------------
  implicit none
  real(kind=8),intent(in)::t_1,t_2,smet
  real(kind=8),intent(out)::EW
  real(kind=8)::a,b,c,ts,metalscale,imfboost

  imfboost=0.31430400d0/0.224468d0  !K01 to Chabrier

  !--- Fitting parameters --- 
  a=imfboost*1.9d54/1.0d6/2.d33 !Scale to per gram
  b=0.0101565d0
  c=0.41017d0
  ts=6.5d6 
  EW=0.0d0

  metalscale=a*(smet/b)**c

  if(t_2.le.ts) then
     EW=metalscale*(t_2-t_1)/ts  !Linear mass loss, m(t2)-m(t1)
  endif

  if(t_1.lt.ts.and.t_2.gt.ts) then
     EW=metalscale*(ts-t_1)/ts  !Linear mass loss, m(t2)-m(t1)
  endif

  if(t_1.ge.ts) then
     EW=0.0
  endif

end SUBROUTINE E_w

!###########################################################
!###########################################################
!###########################################################
!----------------------------------------------
SUBROUTINE agemass(time,met,mass)
!---------------------------------------
  implicit none
  real*8::a0,a1,a2,a,b,c,zzz
  real(kind=8),intent(in)::time, met
  real(kind=8),intent(out)::mass
  
  !      IMPLICIT REAL*8 (A-H,L-Z)
  !
  !     Following Raiteri et al.
  !
  !     Masses: 0.6-120.0 M_sun and Z: 7e-5 to 3e-2
  !
  if(met.lt.7.0d-5) then
     zzz=7.0d-5
  else
     zzz=met
  endif
  if(met.gt.3.0d-2) then
     zzz=3.0d-2
  else
     zzz=met
  endif
  
  a0=10.13+0.07547*log10(zzz)-0.008084*(log10(zzz))**2 
  a1=-4.424-0.7939*log10(zzz)-0.1187*(log10(zzz))**2
  a2=1.262+0.3385*log10(zzz)+0.05417*(log10(zzz))**2
  
  c=(-log10(time)+a0)
  b=a1
  a=a2
  
  if(b*b-4.*a*c.ge.0.0) then 
     mass=-b-sqrt(b*b-4.0*a*c)
     mass=mass/(2.0*a)
     mass=10.0**mass
  else
     mass=120.0
  endif
  
END SUBROUTINE agemass

!----------------------------------------------
SUBROUTINE trapzd2(func,a,b,s,n)
!---------------------------------------
  INTEGER,intent(in):: n
  REAL(kind=8), intent(out)::s
  real(kind=8), intent(in)::a,b
  real(kind=8) func
  EXTERNAL func
  INTEGER it,j
  REAL(kind=8) del,sum,tnm,x
  
  if (n.eq.1) then
     s=0.5*(b-a)*(func(a)+func(b))
  else
     it=2**(n-2)
     tnm=it
     del=(b-a)/tnm
     x=a+0.5*del
     sum=0.
     do j=1,it
        sum=sum+func(x)
        x=x+del
     enddo
     s=0.5*(s+(b-a)*sum/tnm)
  endif
  RETURN
END SUBROUTINE trapzd2

!-------------------------------------------- Simpson integration, s is returned as the integral
SUBROUTINE qsimp(func,a,b,s)
!---------------------------------------
  INTEGER JMAX
  REAL(kind=8), intent(in):: a,b
  REAL(kind=8), intent(out):: s
  REAL(kind=8):: EPS
  REAL(kind=8):: func
  EXTERNAL func
  PARAMETER (EPS=1.e-6, JMAX=20)
  INTEGER j
  REAL(kind=8) os,ost,st
  OST=-1.D30
  OS= -1.D30
  DO J=1,JMAX           
     CALL trapzd2(func,a,b,st,j)
     s=(4.*st-ost)/3.
     if (j.gt.5) then
        if (abs(s-os).lt.EPS*abs(os).or.(s.eq.0..and.os.eq.0.)) return
     endif
     os=s
     ost=st
  enddo
  !        pause 'qsimp did not converge'
END SUBROUTINE qsimp



!--------------------------------------- 
FUNCTION SNIIO(mass)  
!---------------------------------------
  implicit none
  REAL(kind=8) :: A,ind,SNIIO,MO
  REAL(kind=8), intent(in) :: mass
  
!  A=0.22446846861563022  !K01
  A=0.31430400074671366d0  !Chabrier

  ind=-2.3d0
  MO=4.586d-4*mass**2.721
  SNIIO=MO*A*mass**ind !IMF

  end FUNCTION SNIIO
!--------------------------------------- 
FUNCTION SNIIFe(mass)  
!---------------------------------------
  implicit none
  REAL(kind=8) :: A,ind,SNIIFe,MFe
  REAL(kind=8), intent(in) :: mass
  
!   A=0.22446846861563022   !K01
  A=0.31430400074671366d0  !Chabrier 

 ind=-2.3d0
  MFe=2.802d-4*mass**1.864
  SNIIFe=MFe*A*mass**ind !IMF

end FUNCTION SNIIFe
!--------------------------------------- 
FUNCTION SNIIej(mass)  
!---------------------------------------
  implicit none
  REAL(kind=8) :: A,ind,SNIIej,Mej
  REAL(kind=8), intent(in) :: mass
  
!  A=0.22446846861563022 !Kroupa
  A=0.31430400074671366d0  !Chabrier

  ind=-2.3d0
!  Mej=0.7682*mass**1.056 !original Raiteri
  Mej=0.5*mass**1.056 !Modified for Kroupa and calibrated on SB99
  SNIIej=Mej*A*mass**ind !IMF
end FUNCTION SNIIej
!--------------------------------------- 
FUNCTION NSNII(mass)  
!---------------------------------------
  implicit none
  REAL(kind=8) :: A,ind,NSNII
  REAL(kind=8), intent(in) :: mass
  
!  A=0.22446846861563022  !K01
  A=0.31430400074671366d0  !Chabrier
  ind=-2.3d0
  NSNII=A*mass**ind !IMF
   
END FUNCTION NSNII
!--------------------------------------- 
  FUNCTION NSNIa(msec)  
!---------------------------------------
    implicit none
    REAL(kind=8) :: Minf,Msup,ind,SNIafrac,A,NSNIa
    REAL(kind=8), intent(in) :: msec  

  SNIafrac=0.16d0 !Bergh & McClure 1994, rate of SN per century in a MW type galaxy
!  A=0.22446846861563022   !K01
  A=0.31430400074671366d0  !Chabrier

  Minf=max(2.0*msec,3.0d0)
  Msup=msec+8.0d0
  ind=-3.3d0  !Kroupa has -2.3 for M>Msun. Subtract 2 due 1/Mb**2 and the add 1 for prim. func.
  NSNIa=Msup**ind-Minf**ind !Assuming Kroupa IMF
  NSNIa=SNIafrac*A*(NSNIa/(ind))*msec**2  
  
    RETURN
  END FUNCTION NSNIa

!--------------------------------------- 
  FUNCTION MIMF(mass)  
!---------------------------------------
    implicit none
    REAL(kind=8) :: A,ind,mimf
    REAL(kind=8),intent(in)::mass    
    A=0.22446846861563022  !Normalized to Mmax=100 Msun
    if(mass>=0.5) then 
       ind=-2.3
       A=A*1.0
    endif
    if(mass>=0.1.and.mass<0.5) then
       ind=-1.3
       A=A*2.0
    endif
    MIMF=A*mass*mass**ind
  RETURN

  END FUNCTION MIMF

!--------------------------------------- 
FUNCTION MIMFChabrier(mass)  
!---------------------------------------
  implicit none
  REAL(kind=8) :: A,ind,MIMFChabrier
  REAL(kind=8),intent(in)::mass
  A=0.31430400074671366d0 
  if(mass>=1.0) then 
     ind=-2.3
     A=A*1.0
     MIMFChabrier=A*mass*mass**ind
  endif
  if(mass<1.0) then
     A=A*2.2620
     MIMFChabrier=A*mass*0.86d0*exp(-(log10(mass)-log10(0.22))**2/(2.0*0.57**2))
  endif
  RETURN
END FUNCTION MIMFCHABRIER

!--------------------------------------- 
FUNCTION IMFChabrier(mass)  
!---------------------------------------
  implicit none
  REAL(kind=8) :: A,ind,IMFChabrier
  REAL(kind=8),intent(in)::mass
  A=0.31430400074671366d0 
  if(mass>=1.0) then 
     ind=-2.3
     A=A*1.0
     IMFChabrier=A*mass**ind
  endif
  if(mass<1.0) then
     A=A*2.2620
     IMFChabrier=A*0.86d0*exp(-(log10(mass)-log10(0.22))**2/(2.0*0.57**2))
  endif
  RETURN
END FUNCTION IMFCHABRIER
