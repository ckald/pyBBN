all: clean compile test

compile: compile_interpolation compile_four_particle

compile_interpolation:
	cd particles/interpolation && python setup.py build_ext --inplace

compile_four_particle:
	cd interactions/four_particle && python setup.py build_ext --inplace

clean:
	cd particles/interpolation && rm -f *.o *.so *.c *.cpp
	cd interactions/four_particle && rm -f *.o *.so *.c *.cpp

test:
	pip install -U nose
	nosetests tests/unit