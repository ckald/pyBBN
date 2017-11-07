all: clean compile test

compile: compile3

compile3:
	cd interactions/four_particle/cpp && c++ -fopenmp -O3 -shared -std=c++11 -I pybind11/include `python3-config --cflags --ldflags` integral.cpp -o integral.so

clean:
	cd interactions/four_particle/cpp && rm -rf *.o *.so *.c *.so.DSYM

test:
	pip3 install -U nose
	python3 -m "nose" tests/unit