COMPILER ?= c++
PYTHON ?= python3.5
all: integral.so

integral.so: integral.cpp Ds.cpp integral.h
	$(COMPILER) -fopenmp -fPIC -O3 -shared -std=c++11 -I pybind11/include \
	`$(PYTHON)-config --cflags --ldflags` -lgsl -lgslcblas -lm -Werror -Wfatal-errors \
	integral.cpp Ds.cpp -o integral.so

clean:
	rm -rf *.o *.so *.c *.so.DSYM
