all: clean compile test

compile: compile3


compile3:
	cd interactions/four_particle/cpp && c++ -fopenmp -O3 -shared -std=c++11 -fPIC -I pybind11/include `python3.5-config --cflags --ldflags` integral.cpp -o integral.so

clean:
	cd interactions/four_particle/cpp && rm -rf *.o *.so *.c *.so.DSYM

test:
	python3.5 -m "nose" tests/unit