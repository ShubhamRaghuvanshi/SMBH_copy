#ifdef SINKTEST
subroutine backup_sink(filename, filename_desc)
  use amr_commons
  use pm_commons
  use dump_utils, only : dump_header_info, generic_dump, dim_keys
#ifdef RT
  use rt_parameters,only: rt_AGN
#endif
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer,parameter::tag=1135
  integer::dummy_io,info2
#endif

  character(LEN=*), intent(in)::filename, filename_desc

  integer::idim,i,ilevel, unit_out, unit_info
  character(LEN=80)::fileloc
  character(LEN=5)::nchar
  character(len=100) :: field_name
  real(dp),allocatable,dimension(:)::xdp
  integer,allocatable,dimension(:)::ii

  logical :: dump_info
  integer :: ivar

  if(.not. sink) return

  if(verbose) write(*,*)'Entering backup_sink'

  call title(myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)

  ! Set ivar to 1 for first variable
  ivar = 1

  ! Wait for the token
#ifndef WITHOUTMPI
  if(IOGROUPSIZE>0) then
     if (mod(myid-1,IOGROUPSIZE)/=0) then
        call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
             & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
     end if
  endif
#endif

  fileloc = TRIM(filename) // TRIM(nchar)
  if (myid == 1) then
     open(newunit=unit_info, file=trim(filename_desc), form='formatted')
     call dump_header_info(unit_info)
     dump_info = .true.
  else
     dump_info = .false.
  end if

  open(newunit=unit_out,file=TRIM(fileloc),form='unformatted')
  rewind(unit_out)

  write(unit_out)nsink
  write(unit_out)nindsink
  if(nsink>0)then
     allocate(ii(1:nsink))
     ! Write identity sink
     do i=1,nsink
        ii(i)=idsink(i)
     end do
     call generic_dump("identity", ivar, ii, unit_out, dump_info, unit_info)
     deallocate(ii)
     allocate(xdp(1:nsink))
     ! Write mass
     do i=1,nsink
        xdp(i)=msink(i)
     end do
     call generic_dump("mass", ivar, xdp, unit_out, dump_info, unit_info)
     ! Write position
     do idim=1,ndim
        do i=1,nsink
           xdp(i)=xsink(i,idim)
        end do
        call generic_dump("position_"//dim_keys(idim), ivar, xdp, unit_out, dump_info, unit_info)
     enddo
     ! Write velocity
     do idim=1,ndim
        do i=1,nsink
           xdp(i)=vsink(i,idim)
        end do
        call generic_dump("velocity_"//dim_keys(idim), ivar, xdp, unit_out, dump_info, unit_info)
     enddo
     ! Write time
     do i=1,nsink
        xdp(i)=tsink(i)
     end do
     call generic_dump("birth_time", ivar, xdp, unit_out, dump_info, unit_info)
     ! Write real accretion
     do i=1,nsink
        xdp(i)=dMsmbhoverdt_out(i)
     enddo
     call generic_dump("dMsmbh", ivar, xdp, unit_out, dump_info, unit_info)
     ! Write Bondi accretion
     do i=1,nsink
        xdp(i)=dMBH_coarse(i)
     end do
     call generic_dump("dMBH_coarse", ivar, xdp, unit_out, dump_info, unit_info)
     ! Write Eddington accretion
     do i=1,nsink
        xdp(i)=dMEd_coarse(i)
     end do
     call generic_dump("dMEd_coarse", ivar, xdp, unit_out, dump_info, unit_info)
     ! Write Esave
     do i=1,nsink
        xdp(i)=Esave(i)
     end do
     call generic_dump("Esave", ivar, xdp, unit_out, dump_info, unit_info)
     ! Write gas spin axis
     do idim=1,ndim
        do i=1,nsink
           xdp(i)=jsink(i,idim)
        end do
        call generic_dump("jsink_" // dim_keys(idim), ivar, xdp, unit_out, dump_info, unit_info)
     enddo
     ! Write BH spin axis
     do idim=1,ndim
        do i=1,nsink
           xdp(i)=bhspin(i,idim)
        end do
        call generic_dump("spin_" // dim_keys(idim), ivar, xdp, unit_out, dump_info, unit_info)
     enddo
     ! Write BH spin amplitude (signed)
     do i=1,nsink
        xdp(i)=spinmag(i)
     end do
     call generic_dump("spin_magnitude", ivar, xdp, unit_out, dump_info, unit_info)
     ! Write BH efficiency
     do i=1,nsink
        xdp(i)=eps_sink(i)
     end do
     call generic_dump("eps_sink", ivar, xdp, unit_out, dump_info, unit_info)
     ! Write sink_stat
     do idim=1,ndim*2+1
        do ilevel=levelmin,nlevelmax
           do i=1,nsink
              xdp(i)=sink_stat(i,ilevel,idim)
           end do
           write(field_name, "('sink_stat_', i0.2, '_', i0.2)") ilevel, idim
           call generic_dump(field_name, ivar, xdp, unit_out, dump_info, unit_info)
        enddo
     enddo
     ! AGNRT
#ifdef RT
     if(rt_AGN) then
        ! Write AGN radiation to be released
        do i=1,nsink
           xdp(i)=LAGN_coarse(i)
        end do
        call generic_dump("LAGN_coarse", ivar, xdp, unit_out, dump_info, unit_info)
     endif
#endif
     !/AGNRT

     deallocate(xdp)
  endif
  close(unit_out)
  if (myid == 1) close(unit_info)

  ! Send the token
#ifndef WITHOUTMPI
  if(IOGROUPSIZE>0) then
     if(mod(myid,IOGROUPSIZE)/=0 .and.(myid.lt.ncpu))then
        dummy_io=1
        call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag, &
             & MPI_COMM_WORLD,info2)
     end if
  endif
#endif

end subroutine backup_sink



subroutine output_sink(filename)
  use amr_commons
  use pm_commons
  implicit none
  character(LEN=80)::filename

  integer::isink
  integer::nx_loc,ilun
  real(dp)::scale,dx_min
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m
  character(LEN=80)::fileloc

  if(verbose)write(*,*)'Entering output_sink'

  ilun=myid+10

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m=scale_d*scale_l**3d0
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  !dx_min=scale*0.5D0**(nlevelmax-nlevelsheld)/aexp
  dx_min=scale*0.5D0**(nlevelmax)/aexp
  if(verbose)write(*,*)'Entering output_sink'

  ilun=2*ncpu+myid+10

  fileloc=TRIM(filename)
  open(unit=ilun,file=TRIM(fileloc),form='formatted',status='replace')
  !======================
  ! Write sink properties
  !======================
  write(ilun,*)'Number of sink = ',nsink

  write(ilun,'(" ================================================================================================================================== ")')
  write(ilun,'("        Id       Mass(Msol)             x                y                z               vx               vy               vz      ")')
  write(ilun,'(" ================================================================================================================================== ")')

  do isink=1,nsink
     write(ilun,'(I10,7(2X,E23.15))')idsink(isink),msink(isink)*scale_m/2d33,xsink(isink,1:ndim),vsink(isink,1:ndim)
  end do
  write(ilun,'(" ================================================================================================================================== ")')
  close(ilun)

end subroutine output_sink

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_sink_csv(filename)
  use amr_commons
  use pm_commons
  implicit none
  character(LEN=80)::filename,fileloc

  integer::ilun,isink

  if(verbose)write(*,*)'Entering output_sink_csv'

  ilun=2*ncpu+myid+10

  fileloc=TRIM(filename)
  open(unit=ilun,file=TRIM(fileloc),form='formatted',status='replace', recl=500)
  !======================
  ! Write sink properties
  !======================
  do isink=1,nsink
     write(ilun,'(I10,9(A1,ES20.10))')idsink(isink),',',msink(isink),&
          ',',xsink(isink,1),',',xsink(isink,2),',',xsink(isink,3),&
          ',',vsink(isink,1),',',vsink(isink,2),',',vsink(isink,3),&
          ',',t-tsink(isink),',',dMBHoverdt(isink)
  end do

  close(ilun)

end subroutine output_sink_csv
!################################################################
!################################################################
!################################################################
 subroutine create_sink_fine_file(idsink_loc)
  use pm_commons
  use amr_commons
  implicit none

  !===========================================
  !This routine creates the SINKFILE filled at
  !each fine timestep using write_sink_fine
  !===========================================
  integer::idsink_loc
  character(LEN=80)::filename,filedir,filecmd,fileinfo
  integer::ilun,info,i
  character(LEN=5)::nchar

  filedir='SINK'
  if(sink.and.myid.eq.1)then
     filecmd='mkdir -p '//TRIM(filedir)
#ifdef NOSYSTEM
     call PXFMKDIR(TRIM(filedir),LEN(TRIM(filedir)),O'755',info)
#else
     call system(filecmd)
#endif
     ilun=13*ncpu+10+idsink_loc
     write(*,*)"CREATING SINK FILE",ilun
     call title(idsink_loc,nchar)
     write(*,*)filedir,idsink_loc,ilun,nchar
     filename='SINK/sink_'//TRIM(nchar)//'.txt'
     open(unit=ilun,file=TRIM(filename),form='formatted',status='unknown',position='append')
     write(ilun,*)'#1      nstep_fine'
     write(ilun,*)'#2      nstep_coarse'
     write(ilun,*)'#3      aexp'
     write(ilun,*)'#4      time [yr]'
     write(ilun,*)'#5      BHmass [Msun]'
     write(ilun,*)'#6      levelsink 1 '
     write(ilun,*)'#7      levelmax_current 1'
     write(ilun,*)'#8      position 3 [pc]'
     write(ilun,*)'#9     vsink 3 [km/s]'
     write(ilun,*)'#10     dotMbondi [Msun/yr]'
     write(ilun,*)'#11     dotMedd [Msun/yr]'
     write(ilun,*)'#12     rg_scale 1 [pc]  (MIN(rbondi,racc))'
     write(ilun,*)'#13     spinmag  (black hole spin magnitude)'
     write(ilun,*)'#14     dMsmbh 1 [Msun] (mass for feedback built up'
     write(ilun,*)'#15     Esave 1 [erg] (feedback energy leftover from last step)'
     write(ilun,*)'#16     bhspin 3 [?]   (bh spin vector)'
     write(ilun,*)'#17     jsink 3     (gas angular momentum vector)'
     write(ilun,*)'#18     Efeed 1 [erg] (feedback energy available for feedback)'
     write(ilun,*)'#19     dpacc 3 [code units] ( accreted momentum)'
     write(ilun,*)'#20     fsink 3 [code units] (force on sink)'
  endif

end subroutine create_sink_fine_file
!################################################################
!################################################################
!################################################################
subroutine write_sink_fine
  use pm_commons
  use amr_commons
  use hydro_commons
  implicit none

!==============================================
! This routine prints sink information on each
! fine timestep to follow black hole evolution
! in more detail
!=============================================
  integer::isink,ilun
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m,scale_E,scale_a
  real(dp)::unit_amu,unit_pc,unit_msun,unit_dotM

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m=scale_d*scale_l**3
  scale_E=scale_m*scale_v**2
  scale_a=scale_l/scale_t**2
  unit_amu=1.660538921e-24
  unit_pc=3.08567758096d18   !Possibly 3.08d18?
  unit_msun=1.98841586d33
  unit_dotM=(scale_m/scale_t)/unit_msun*3600*24*365

  do isink=1,nsink
     if(myid.eq.1)then
        ilun=13*ncpu+10+idsink(isink)
        write(ilun,42)nstep,nstep_coarse,aexp,t*scale_t/(3600*24*365.25),msink(isink)*scale_m/unit_msun,levelsink(isink),nlevelmax_current,xsink(isink,:)*scale_l/unit_pc,vsink(isink,:)*scale_v/1d5,dMBHoverdt(isink)*unit_dotM,dMEdoverdt(isink)*unit_dotM,rg_scale(isink)*scale_l/unit_pc,spinmag(isink),dMsmbh(isink)*scale_m/unit_msun,Esave(isink)*scale_E,bhspin(isink,:),jsink(isink,:),Efeed(isink)*scale_E,acc_mom_all(isink,:),fsink(isink,:)*scale_a
     endif
  enddo
42 format(2(i8,1x),3(e23.15,1x),2(i8,1x),36(e23.15,1x))

  Efeed=0d0   !Reset to zero as only for writing to file purposes
  drag_level=0

end subroutine write_sink_fine
!################################################################
!################################################################
!################################################################

subroutine create_sink_dir
  character(LEN=80)::filedir,filecmd
  logical :: directory_exists
  integer :: ierr

  filedir='sinkdata/'
  INQUIRE(file=TRIM(filedir), exist=directory_exists)
  if(.not. directory_exists) then 
    filecmd='mkdir -p '//TRIM(filedir)
    ierr=1
    call EXECUTE_COMMAND_LINE(filecmd,exitstat=ierr)
    if(ierr.ne.0 .and. ierr.ne.127)then
      write(*,*) 'Error - Could not create ',TRIM(filedir),' error code=',ierr
      call clean_stop
    endif   
  endif
end subroutine create_sink_dir  

subroutine write_sink_props
  use amr_commons
  use pm_commons
  implicit none
  character(LEN=80)::filename,filedir,fileloc
  character(LEN=5) :: sinknum_char
  integer::isink,sinkid
  logical :: directory_exists
  real(dp)::msink_total
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m
  
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m=scale_d*scale_l**3d0

  !if(verbose)then
    if(myid==1) write(*,*)'Entering write_sink_props', nsink
  !endif

  if(nsink > 0 ) then
    msink_total=0.0
    do isink=1,nsink
      msink_total = msink_total + msink(isink)
    end do 

    do isink=1,nsink
      sinkid = idsink(isink)
      call title(sinkid,sinknum_char)
      filedir='sinkdata/'
      filename=TRIM(filedir)//'sink_'//TRIM(sinknum_char)//'.csv' 
      fileloc=TRIM(filename)

      OPEN(unit=346, file=TRIM(fileloc), form='formatted', status="unknown", position="append", action="write") 
      write(346,'(I10,13(A1,ES20.10))')nsink,',',msink(isink),&
          ',',xsink(isink,1),',',xsink(isink,2),',',xsink(isink,3),&
          ',',vsink(isink,1),',',vsink(isink,2),',',vsink(isink,3),&
          ',',t-tsink(isink),',',dMsmbhoverdt_out(isink), ',' ,dMBHoverdt(isink) ',', msink_total, ',', mstar_total_all &
          , ',', aexp 
    enddo   
    close(346)
  endif 

  if(msink_total < 0) then 
    write(*,*) "Negative mass sink", nsink
    do isink=1,nsink
      write(*,*) "msink :", isink, idsink(isink), msink(isink)
    end do 
    call clean_stop
  endif  

end subroutine write_sink_props

#endif
