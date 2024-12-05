#############################################################################
# If you have problems with this makefile, contact Romain.Teyssier@gmail.com
#############################################################################
# Compilation time parameters

# Do we want a debug build? 1=Yes, 0=No
DEBUG=0
# Compiler flavor: GNU or INTEL
COMPILER=INTEL
# Size of vector cache
NVECTOR=256
# Number of dimensions
NDIM=3
# Float precision size
NPRE=8
# hydro/mhd solver
SOLVER=hydro
# Patch
PATCH=
# Use RT? 1=Yes, 0=No
RT=1
# Use MPI? 1=Yes, 0=No
MPI=1
MPIF90 = mpiifort
# Root name of executable
EXEC = ramses_rt_smbh
# Use Grackle cooling? 1=Yes, 0=No
GRACKLE=0
# Number of metal species
NMETALS=0
# ATON flages: Uncomment to enable ATON
ATON_FLAGS = #-DATON
# Number of ions for RT
NIONS=3
# Number of photon groups for RT
NGROUPS=0
# Number of cosmic ray groups
NCR=0
# Number of additional energies
NENER = 0+$(NCR)
# Number of passive scalars
NPSCAL = 2
# Activate SINK particles 
SINK=1
AGNRT=0
#Activate MC tracers 
MC_TRACER=0

# Compute NVAR
NVAR = 2+$(NDIM)#+3 for metallicity and hydroscalar for virial sf+IR Trap
ifeq ($(SOLVER),mhd)
   NVAR = 8
endif
NVAR := $(NVAR)+$(NENER)+$(NMETALS)+$(NPSCAL)
ifeq ($(RT),1)
   NVAR := $(NVAR)+$(NIONS)
endif

# Set to one to use 'include "mpif.h"' instead of more recent "use mpi"
OLD_MPI_SUPPORT = 0
#############################################################################
GITBRANCH = $(shell git rev-parse --abbrev-ref HEAD)
GITHASH = $(shell git log --pretty=format:'%H' -n 1)
GITREMOTE = $(shell git config --get branch.$(GITBRANCH).remote)
GITREPO = $(shell git config --get remote.$(GITREMOTE).url)
BUILDDATE = $(shell date +"%D-%T")
DEFINES = -DNVECTOR=$(NVECTOR) -DNDIM=$(NDIM) -DNPRE=$(NPRE) -DNENER=$(NENER) \
          -DNVAR=$(NVAR) -DSOLVER$(SOLVER) -DNCR=$(NCR) $(ATON_FLAGS)
ifeq ($(GRACKLE),1)
   DEFINES += -Dgrackle
endif
ifeq ($(OLD_MPI_SUPPORT),1)
   DEFINES += -DMPI_OLD
endif
ifeq ($(RT),1)
   DEFINES += -DRT -DNIONS=$(NIONS) -DNGROUPS=$(NGROUPS)
endif

ifeq ($(SINK),1)
   DEFINES += -DSINKTEST
endif

ifeq ($(AGNRT),1)
   DEFINES += -DSINKAGN_RT
endif

ifeq ($(MC_TRACER),1)
   DEFINES += -DMC_tracer
endif


#############################################################################
# Fortran compiler options and directives

# GNU compiler (gfortran)
ifeq ($(COMPILER),GNU)
   FFLAGS = -x f95-cpp-input $(DEFINES)
   ifeq ($(MPI),1)
      F90 = $(MPIF90)
   else
      F90 = gfortran
      FFLAGS += -DWITHOUTMPI
   endif
   F90 += -frecord-marker=4 -fbacktrace -ffree-line-length-none -g -fimplicit-none
   ifeq ($(DEBUG),1)
      F90 += -O0 -fbounds-check -Wuninitialized -Wall
      FFLAGS += -ffpe-trap=zero,underflow,overflow,invalid -finit-real=nan
   else
      F90 += -O3
   endif
endif

# Intel compiler
ifeq ($(COMPILER),INTEL)
   FFLAGS = -cpp $(DEFINES)
   ifeq ($(MPI),1)
      F90 = $(MPIF90)
      FFLAGS += -DNOSYSTEM
   else
      F90 = ifort
      FFLAGS += -DWITHOUTMPI
   endif
   F90 += -fp-model source
   ifeq ($(DEBUG),1)
      F90 += -warn all -O0 -g -traceback -check bounds
      FFLAGS += -fpe0 -ftrapuv -init=zero -init=snan -init=arrays
   else
      F90 += -O0
   endif
endif

#############################################################################
MOD = mod
#############################################################################
# MPI librairies
LIBMPI =
#LIBMPI = -lfmpi -lmpi -lelan

# --- CUDA libraries, for Titane ---
LIBCUDA = -L/opt/cuda/lib  -lm -lcuda -lcudart

ifeq ($(GRACKLE),1)
   # Add include and library install path for grackle and hdf5 here
   LIBS_GRACKLE = -L$(HOME)/local/lib -lgrackle -lhdf5 -lz -lgfortran -ldl
   LIBS_OBJ     = -I$(HOME)/local/include -DCONFIG_BFLOAT_8 -DH5_USE_16_API -fPIC
endif
LIBS = $(LIBMPI) $(LIBS_GRACKLE)
#############################################################################
# Sources directories are searched in this exact order
#VPATH = $(shell [ -z $(PATCH) ] || find $(PATCH) -type d):../$(SOLVER):../aton:
VPATH = $(PATCH):../$(SOLVER):../aton:
ifeq ($(RT),1)
   VPATH += ../rt:
endif
ifeq ($(SOLVER),mhd)
   VPATH += ../conduction:
endif
ifneq ($(NCR),0)
   VPATH += ../cr:
endif
VPATH += ../hydro:../pm:../poisson:../amr:../io:../amr/geometry
#############################################################################
# All objects
MODOBJ = mpi_mod.o amr_parameters.o amr_commons.o random.o pm_parameters.o \
         pm_commons.o poisson_parameters.o dump_utils.o constants.o
ifeq ($(GRACKLE),1)
   MODOBJ += grackle_parameters.o
endif
MODOBJ += poisson_commons.o hydro_parameters.o hydro_commons.o \
          cooling_module.o bisection.o sparse_mat.o clfind_commons.o \
          gadgetreadfile.o write_makefile.o write_patch.o write_gitinfo.o
ifeq ($(RT),1)
   MODOBJ += rt_parameters.o rt_hydro_commons.o coolrates_module.o \
             rt_spectra.o rt_cooling_module.o rt_flux_module.o
endif
AMROBJ = read_params.o init_amr.o init_time.o init_refine.o adaptive_loop.o \
         amr_step.o update_time.o output_amr.o flag_utils.o \
         physical_boundaries.o virtual_boundaries.o refine_utils.o \
         nbors_utils.o hilbert.o load_balance.o title.o sort.o cooling_fine.o \
         units.o light_cone.o movie.o memory.o end.o

# Particle-Mesh objects
PMOBJ = init_part.o output_part.o rho_fine.o synchro_fine.o move_fine.o \
        newdt_fine.o particle_tree.o add_list.o remove_list.o star_formation.o \
        feedback.o clump_finder.o clump_merger.o 
ifeq ($(SINK),1)
   PMOBJ += utils.o geom_types.o geom_utils.o geom_distances.o geom_volumes.o \
						sink_particle.o flag_formation_sites.o init_sink.o output_sink.o 
endif
ifeq ($(MC_TRACER),1)
   PMOBJ += tracer_utils.o sink_particle_tracer.o
endif

# Poisson solver objects
POISSONOBJ = init_poisson.o phi_fine_cg.o interpol_phi.o force_fine.o \
             multigrid_coarse.o multigrid_fine_commons.o multigrid_fine_fine.o \
             multigrid_fine_coarse.o gravana.o boundary_potential.o rho_ana.o \
             output_poisson.o
# Hydro objects
HYDROOBJ = init_hydro.o init_flow_fine.o write_screen.o output_hydro.o \
           courant_fine.o godunov_fine.o uplmde.o umuscl.o interpol_hydro.o \
           godunov_utils.o condinit.o hydro_flag.o hydro_boundary.o boundana.o \
           read_hydro_params.o synchro_hydro_fine.o
# Conduction objects
CONDOBJ =  conduction.o cond_split.o cond_split_semi.o conduction_fine.o \
          coupling_fine.o
# Cosmic ray objects
CROBJ =  crdiffusion.o crdiff_split.o crdiffusion_fine.o
# RT objects
RTOBJ = rt_init_hydro.o rt_init_xion.o rt_init.o rt_init_flow_fine.o \
        rt_output_hydro.o rt_godunov_fine.o rt_interpol_hydro.o \
        rt_godunov_utils.o rt_condinit.o rt_hydro_flag.o rt_hydro_boundary.o \
        rt_boundana.o rt_units.o
# Patch objects
sinclude $(PATCH)/Makefile

# All objects
AMRLIB = $(AMROBJ) $(HYDROOBJ) $(PMOBJ) $(POISSONOBJ)
ifeq ($(RT),1)
   AMRLIB += $(RTOBJ)
endif
ifeq ($(SOLVER),mhd)
   AMRLIB += $(CONDOBJ)
endif
ifneq ($(NCR),0)
   AMRLIB += $(CROBJ)
endif
# ATON objects
ATON_MODOBJ = timing.o radiation_commons.o rad_step.o
ATON_OBJ = observe.o init_radiation.o rad_init.o rad_boundary.o rad_stars.o \
           rad_backup.o ../aton/atonlib/libaton.a
#############################################################################
ramses:	$(MODOBJ) $(AMRLIB) ramses.o
	$(F90) $(MODOBJ) $(AMRLIB) ramses.o -o $(EXEC)$(NDIM)d $(LIBS)
	rm write_makefile.f90
	rm write_patch.f90
ramses_aton: $(MODOBJ) $(ATON_MODOBJ) $(AMRLIB) $(ATON_OBJ) ramses.o
	$(F90) $(MODOBJ) $(ATON_MODOBJ) $(AMRLIB) $(ATON_OBJ) ramses.o \
		-o $(EXEC)$(NDIM)d $(LIBS) $(LIBCUDA)
	rm write_makefile.f90
	rm write_patch.f90
#############################################################################
write_gitinfo.o: FORCE
	$(F90) $(FFLAGS) -DPATCH='"$(PATCH)"' -DGITBRANCH='"$(GITBRANCH)"' \
		-DGITHASH='"$(GITHASH)"' -DGITREPO='"$(GITREPO)"' \
		-DBUILDDATE='"$(BUILDDATE)"' -c ../amr/write_gitinfo.f90 -o $@
write_makefile.o: FORCE
	../utils/scripts/cr_write_makefile.sh $(MAKEFILE_LIST)
	$(F90) $(FFLAGS) -c write_makefile.f90 -o $@
write_patch.o: FORCE
	../utils/scripts/cr_write_patch.sh $(PATCH)
	$(F90) $(FFLAGS) -c write_patch.f90 -o $@
%.o:%.F
	$(F90) $(FFLAGS) -c $^ -o $@ $(LIBS_OBJ)
%.o:%.f90
	$(F90) $(FFLAGS) -c $^ -o $@ $(LIBS_OBJ)
FORCE:
#############################################################################
clean:
	rm -f *.o *.$(MOD) *.i
#############################################################################
