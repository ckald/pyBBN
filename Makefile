<<<<<<< HEAD
COMPILER ?= c++
=======
COMPILER ?= gcc
>>>>>>> 259ef64... Add a compiler flag
PYTHON ?= python3

all: clean compile test

compile: compile3

compile3:
	cd interactions/four_particle/cpp \
	&& $(COMPILER) -fopenmp -fPIC -O3 -shared -std=c++11 -I pybind11/include \
	`$(PYTHON)-config --cflags --ldflags` -lgsl -lgslcblas -lm -Wfatal-errors \
	integral.cpp Ds.cpp -o integral.so

clean:
	cd interactions/four_particle/cpp && rm -rf *.o *.so *.c *.so.DSYM

test:
	$(PYTHON) -m "nose" tests/unit