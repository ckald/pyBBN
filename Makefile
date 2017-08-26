all: clean compile test

compile: compile3

compile3:
	cd interactions/four_particle/cpp && c++ -O3 -shared -std=c++11 -I pybind11/include `python3-config --cflags --ldflags` integral.cpp -o integral.so
	# cd interactions/four_particle/cpp && c++ -fopenmp -O3 -shared -std=c++11 -I pybind11/include `python3-config --cflags --ldflags` integral.cpp -o integral.so

clean:
	cd interactions/four_particle/cpp && rm -rf *.o *.so *.c *.so.DSYM

# compile: compile3

# compile2:
# 	cd interactions/four_particle && python2 setup.py build_ext --inplace
# compile3:
# 	cd interactions/four_particle && python3 setup.py build_ext --inplace


# clean:
# 	cd interactions/four_particle && rm -f *.o *.so *.c *.cpp

test:
	pip3 install -U nose
	python3 -m "nose" tests/unit