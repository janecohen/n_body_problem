
# Nbody Makefile

FC = gfortran

FFLAGS = -O0 -fno-second-underscore -Wall -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow

.SUFFIXES: .f90

# Default target
ALL: Test
	$(MAKE) PartialClean
	@echo "!!! Compilation OK !!!"

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

OBJ = initialize.o forces.o integrate.o

# Module dependencies
integrate.o: initialize.o forces.o
main.o: initialize.o forces.o integrate.o

# Compile executable
Test: $(OBJ) main.o
	$(FC) $(FFLAGS) main.o $(OBJ) -o main

# Remove intermediate files
PartialClean:
	rm -f *.o *.mod

# Remove everything
clean:
	rm -f *.o *.mod main

.PHONY: ALL Test PartialClean clean