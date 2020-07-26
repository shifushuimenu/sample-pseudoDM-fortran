.SUFFIXES: 
.SUFFIXES: .f90 .o

include ./make.sys.gfortran

module = 

default: all

all: sample_pseudo_DM

sample_pseudo_DM: sample_pseudo_DM.f90  make.sys.gfortran
	${F90} ${LFLAGS} -o sample_pseudo_DM sample_pseudo_DM.f90


.o.f90:
	${F90} ${FFLAGS} -c $*.f90

clean:
	rm -f *.o
	rm -f *.mod

clean-dat:
	rm -f *.dat 
	rm -f *.out 

clean-all: clean clean-dat
