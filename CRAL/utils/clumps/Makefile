# Makefile for RAMSES utils
F90=mpiifort
BINDIR=.
FLAGS=-ffree-line-length-none

MODOBJ=random.o io_ramses.o

# Make all targets

all: read_sink read_star amr2map
amr2map: $(MODOBJ) amr2map.o
	$(F90) $(FLAGS) amr2map.o -o $(BINDIR)/amr2map
read_sink: $(MODOBJ) read_sink.o
	$(F90) $(FLAGS) read_sink.o -o $(BINDIR)/read_sink
read_star: $(MODOBJ) read_star.o
	$(F90) $(FLAGS) read_star.o -o $(BINDIR)/read_star

# Make a specific object file
%.o: %.f90
	$(F90) $(FLAGS) -c $^
clean:
	rm *.o *.mod

# ./read_sink -inp /u/sraghu/sim/SHUB/b128/smbh/output_00015 -out sinkdata -xmi 0 -xma 1 -ymi 0 -yma 1 -zmi 0 -zma 1 	
