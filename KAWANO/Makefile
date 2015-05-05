FLAGS = -freal-4-real-8 -fno-align-commons -Ofast

all: clean python kawano

python:
	f2py -c -m interpolation interpolation.f90 --f90flags="$(FLAGS)"

kawano:
	gfortran kawano_steriles.f interpolation.f90 $(FLAGS)

clean:
	rm -rf *.o *.so *.mod *.out
