PRGM ?= evolve.so

all: $(PRGM)

%.so :	%.f90 interfaces.mod
	f2py -c -m $* --opt="-O3" --f90flags="-Wall -fopenmp -cpp -ffree-form" -lgomp -DF2PY_REPORT_ON_ARRAY_COPY=100 $*.f90

interfaces.mod : interfaces.f90
	gfortran -c interfaces.f90

clean:
	rm -rf *.so *.mod *.o
