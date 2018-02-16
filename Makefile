COMPILER ?= c++
PYTHON ?= python3.5

all: clean compile_3p compile_4p test

compile_3p:
	cd interactions/three_particle/cpp \
	&& $(COMPILER) -fopenmp -fPIC -O0 -g -shared -std=c++11 -I pybind11/include \
	`$(PYTHON)-config --cflags --ldflags` -lm -Wfatal-errors \
	integral.cpp -o integral.so

compile_4p:
	cd interactions/four_particle/cpp \
	&& $(COMPILER) -fopenmp -fPIC -O3 -shared -std=c++11 -I pybind11/include \
	`$(PYTHON)-config --cflags --ldflags` -lm -Wfatal-errors \
	integral.cpp Ds.cpp -o integral.so

clean:
	cd interactions/four_particle/cpp && rm -rf *.o *.so *.c *.so.DSYM ;\
    cd ../../three_particle/cpp && rm -rf *.o *.so *.c *.so.DSYM

test:
	$(PYTHON) -m "nose" tests/unit
