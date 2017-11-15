PYTHON ?= python3

all: clean compile test

compile: compile3

compile3:
	cd interactions/four_particle/cpp \
	&& c++ -fopenmp -fPIC -O3 -shared -std=c++11 -I pybind11/include \
	`$(PYTHON)-config --cflags --ldflags` integral.cpp -o integral.so

clean:
	cd interactions/four_particle/cpp && rm -rf *.o *.so *.c *.so.DSYM

test:
	$(PYTHON) -m "nose" tests/unit
