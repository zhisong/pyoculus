F90=gfortran
F90FLAGS=-ffree-line-length-none
OPT=-O3

SOURCES=pyvariables.f90 pybasefn.f90 pybfield.f90 pycoords.f90 pyPJH.f90

all:
	cd pyoculus/problems/SPECfortran && f2py -c $(SOURCES) -m fortran_module --f90flags=$(F90FLAGS) --opt=$(OPT) --f90exec=$(F90)

clean:
	cd pyoculus/problems/SPECfortran && rm -r fortran_module*

doxygen:
	doxygen doc/Doxyfile
