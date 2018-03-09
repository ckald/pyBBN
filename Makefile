COMPILER ?= c++
PYTHON ?= python3.5
all: clean compile_3p compile_4p test

compile_3p:
	cd interactions/three_particle/cpp \
	&& $(COMPILER) -fopenmp -fPIC -shared -o integral.so integral.cpp -O0 -g -std=c++11 -I pybind11/include \
	`$(PYTHON)-config --cflags --ldflags` -lgsl -lgslcblas -lm -Wfatal-errors \

compile_4p:
	cd interactions/four_particle/cpp \
	&& $(COMPILER) -fopenmp -fPIC -shared -o integral.so integral.cpp Ds.cpp -O3 -std=c++11 -I pybind11/include \
	`$(PYTHON)-config --cflags --ldflags` -lgsl -lgslcblas -lm -Wfatal-errors \


clean:
	cd interactions/four_particle/cpp && rm -rf *.o *.so *.c *.so.DSYM ;\
	cd ../../three_particle/cpp && rm -rf *.o *.so *.c *.so.DSYM

test:
	$(PYTHON) -m "nose" tests/unit
