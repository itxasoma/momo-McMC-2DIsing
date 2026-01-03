FC = gfortran
FFLAGS = -O2 -g -Wall -Wextra -fcheck=all -fbacktrace
TARGET = mdinamics

all: $(TARGET)

$(TARGET): parameters.o lattice.o ising_energy.o main.o r1279.o ran2.o binning.o
	$(FC) $(FFLAGS) -o $@ $^

binning.o: binning.f90
	$(FC) $(FFLAGS) -c $<

ran2.o: ran2.f
	$(FC) $(FFLAGS) -c $<

r1279.o: r1279.f90
	$(FC) $(FFLAGS) -c $<

parameters.o: parameters.f90
	$(FC) $(FFLAGS) -c $<

lattice.o: lattice.f90 parameters.o
	$(FC) $(FFLAGS) -c $<
	
ising_energy.o: ising_energy.f90 parameters.o
	$(FC) $(FFLAGS) -c $<

main.o: main.f90 parameters.o lattice.o ising_energy.o binning.o r1279.o ran2.o
	$(FC) $(FFLAGS) -c $<

run: $(TARGET)
	./$(TARGET) DATA.in

clean:
	rm -f *.o *.mod $(TARGET)
