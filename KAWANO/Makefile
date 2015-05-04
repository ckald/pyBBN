FLAGS = -fdefault-real-8 -fno-align-commons

all: clean python kawano

python:
	f2py -c -m interpolation interpolation.f90

kawano:
	gfortran kawano_steriles.f interpolation.f90 $(FLAGS)

clean:
	rm -rf *.o *.so *.mod *.out
