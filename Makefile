all: clean compile test

compile: compile3

compile2:
	cd interactions/four_particle && python2 setup.py build_ext --inplace
compile3:
	cd interactions/four_particle && python3 setup.py build_ext --inplace


clean:
	cd interactions/four_particle && rm -f *.o *.so *.c *.cpp

test:
	pip3 install -U nose
	python3 -m "nose" tests/unit