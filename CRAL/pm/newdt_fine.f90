subroutine newdt_fine(ilevel)
  use pm_commons
  use amr_commons
  use hydro_commons
  use poisson_commons, ONLY: gravity_type
#ifdef RT
  use rt_parameters, ONLY: rt_advect, rt_nsubcycle
#endif
  use constants, ONLY: pi
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
  integer::ilevel
  !-----------------------------------------------------------
  ! This routine compute the time step using 4 constraints:
  ! 1- a Courant-type condition using particle velocity
  ! 2- the gravity free-fall time
  ! 3- 10% maximum variation for aexp
  ! 4- maximum step time for ATON
  ! This routine also computes the particle kinetic energy.
  !-----------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart
  integer::npart1,ip,i
  !sinktest
#ifdef SINKTEST
  integer::isink,idim
  real(dp)::dx,dx_loc,nx_loc,scale
  real(dp)::vsink2,vsink_max
#endif  
  !sinktest 
  !tracer 
#ifdef MC_tracer
  logical :: ok
#endif 
  !tracer 
  integer,dimension(1:nvector),save::ind_part
  real(kind=8)::dt_loc,dt_all,ekin_loc,ekin_all
  real(dp)::tff,fourpi,threepi2,aexp_next,dt_levelmin,dt_want
#ifdef ATON
  real(dp)::aton_time_step,dt_aton
#endif
#ifdef RT
  real(dp)::dt_rt
#endif

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  threepi2=3.0d0*pi**2

  ! Save old time step
  dtold(ilevel)=dtnew(ilevel)

  ! Maximum time step
  dtnew(ilevel)=boxlen/smallc
  if(poisson.and.gravity_type<=0)then
     fourpi=4.0d0*pi
     if(cosmo)fourpi=1.5d0*omega_m*aexp

  !sinktest     
#ifdef SINKTEST     
  if (sink)then
    tff=sqrt(threepi2/8/fourpi/(rho_max(ilevel)+rho_sink_tff(ilevel)))
  else
    tff=sqrt(threepi2/8/fourpi/(rho_max(ilevel)+smallr))
  end if
#else      
  tff=sqrt(threepi2/8/fourpi/(rho_max(ilevel)+smallr))
#endif 
  !sinktest 


     dtnew(ilevel)=MIN(dtnew(ilevel),courant_factor*tff)
  end if
  if(cosmo)then
     dtnew(ilevel)=MIN(dtnew(ilevel),0.1d0/hexp)
  end if

  !sinktest 
#ifdef SINKTEST
  ! Check sink velocity
  if(sink)then
    dx=0.5d0**ilevel
    nx_loc=dble(icoarse_max-icoarse_min+1)
    scale=boxlen/nx_loc
    dx_loc=dx*scale
    vsink_max=0d0
    do isink=1,nsink
      if(.not. new_born(isink))then
        vsink2=0d0
        do idim=1,ndim
          vsink2=vsink2+vsink(isink,idim)**2
        end do
        vsink_max=MAX(vsink_max,sqrt(vsink2))
      endif
    end do
    if(vsink_max.GT.0d0)then
      dtnew(ilevel)=MIN(dtnew(ilevel),courant_factor*dx_loc/vsink_max)
    endif
  endif
#endif 
  !sinktest 

#ifdef ATON
  ! Maximum time step for ATON
  if(aton)then
     dt_aton = aton_time_step()
     if(dt_aton>0d0)then
        dtnew(ilevel)=MIN(dtnew(ilevel),dt_aton)
     end if
  end if
#endif

#ifdef RT
  ! Maximum time step for radiative transfer
  if(rt_advect)then
     call get_rt_courant_coarse(dt_rt,ilevel)
     dtnew(ilevel) = 0.99999 * &
          MIN(dtnew(ilevel), dt_rt/2.0**(ilevel-levelmin) * rt_nsubcycle)
     if(static) RETURN
  endif
#endif

  if(pic) then

    dt_all=dtnew(ilevel); dt_loc=dt_all
    ekin_all=0; ekin_loc=0

    ! Compute maximum time step on active region
    if(numbl(myid,ilevel)>0)then
      ! Loop over grids
      ip=0
      igrid=headl(myid,ilevel)
      do jgrid=1,numbl(myid,ilevel)
        npart1=numbp(igrid)   ! Number of particles in the grid
        if(npart1>0)then
          ! Loop over particles
          ipart=headp(igrid)
          do jpart=1,npart1
            !tracer 
#ifdef MC_tracer
            if(is_not_tracer(typep(ipart))) then 
#endif 
            !tracer 
            ip=ip+1
            ind_part(ip)=ipart
            if(ip==nvector)then
              call newdt2(ind_part,dt_loc,ekin_loc,ip,ilevel)
              ip=0
            end if
            !tracer 
#ifdef MC_tracer
            endif
#endif 
            !tracer 
            ipart=nextp(ipart)    ! Go to next particle
          end do
          ! End loop over particles
        end if
        igrid=next(igrid)   ! Go to next grid
      end do
      ! End loop over grids
      if(ip>0)call newdt2(ind_part,dt_loc,ekin_loc,ip,ilevel)
    end if

     ! Minimize time step over all cpus
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(dt_loc,dt_all,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
          & MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(ekin_loc,ekin_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
          & MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
     dt_all=dt_loc
     ekin_all=ekin_loc
#endif
     ekin_tot=ekin_tot+ekin_all
     dtnew(ilevel)=MIN(dtnew(ilevel),dt_all)

  end if

  if(hydro)call courant_fine(ilevel)
 
  ! added by Taysun Kimm to output exactly when we want
  if(ilevel.eq.levelmin)then
     dt_levelmin = dtnew(ilevel)
     if(cosmo)then
        ! Find neighboring times
        i=1 
        do while(tau_frw(i)>t+dt_levelmin.and.i<n_frw)
           i=i+1
        end do
        ! Interpolate expansion factor
        aexp_next = aexp_frw(i  )*(t+dt_levelmin-tau_frw(i-1))/(tau_frw(i  )-tau_frw(i-1))+ &
                  & aexp_frw(i-1)*(t+dt_levelmin-tau_frw(i  ))/(tau_frw(i-1)-tau_frw(i  ))  
        if(aexp_next.gt.aout(iout))then
           i=1 
           do while(aexp_frw(i)>aout(iout).and.i<n_frw)
              i=i+1
           end do
           ! Interpolate expansion factor
           dt_want = tau_frw(i  )*(aout(iout)-aexp_frw(i-1))/(aexp_frw(i  )-aexp_frw(i-1))+ &
                   & tau_frw(i-1)*(aout(iout)-aexp_frw(i  ))/(aexp_frw(i-1)-aexp_frw(i  ))  
           dt_want = (dt_want - t)/nsubcycle(ilevel)
           dtnew(ilevel)=min(dtnew(ilevel),max(0.1*dtnew(ilevel),dt_want))
        endif
     else
        dt_want = (tout(iout)-t)/nsubcycle(ilevel)
        dtnew(ilevel)=min(dtnew(ilevel),max(0.001*dtnew(ilevel),dt_want))
     endif
  endif

  !sinktest 
#ifdef SINKTEST
  !Write sink file fine
  if(sink) then 
    if(output_finestep_sink)then
      !if(sink .and. nsink>0)then
      if(nsink>0)then
        if(ilevel==nlevelmax)then
          call write_sink_fine
        else if(numbtot(1,ilevel+1)==0)then
          call write_sink_fine
        endif
      endif
    endif
  endif
#endif 
  !sinktest   

111 format('   Entering newdt_fine for level ',I2)

end subroutine newdt_fine
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine newdt2(ind_part,dt_loc,ekin_loc,nn,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
  real(kind=8)::dt_loc,ekin_loc
  integer::nn,ilevel
  integer,dimension(1:nvector)::ind_part

  integer::i,idim,nx_loc
  real(dp)::dx,dx_loc,scale,dtpart
  real(dp),dimension(1:nvector),save::v2,mmm
  real(dp),dimension(1:nvector,1:ndim)::vvv

  ! Compute time step
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  v2(1:nn)=0.0D0
  do idim=1,ndim
     do i=1,nn
        vvv(i, idim) = vp(ind_part(i), idim)
        v2(i)=max(v2(i),vvv(i, idim)**2)
        ! v2(i)=v2(i)+vp(ind_part(i),idim)**2
     end do
  end do
  do i=1,nn
     if(v2(i)>0.0D0)then
        dtpart=courant_factor*dx_loc/sqrt(v2(i))
        dt_loc=MIN(dt_loc,dtpart)
     end if
  end do

  ! Fetch mass
  do i = 1, nn
     mmm(i) = mp(ind_part(i))
  end do

  ! Compute kinetic energy
  do idim=1,ndim
     do i=1,nn
        ekin_loc=ekin_loc+0.5D0*mmm(i)*vvv(i, idim)**2
     end do
  end do

end subroutine newdt2




