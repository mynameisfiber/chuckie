%.so :	%.f90
	f2py -c -m $* --opt="-O3" --f90flags="-Wall -fopenmp" -lgomp -DF2PY_REPORT_ON_ARRAY_COPY=100 $*.f90

clean:
	rm -rf *.so
