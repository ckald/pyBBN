all: clean compile test

compile:
	cd interactions/four_particle && python setup.py build_ext --inplace

clean:
	cd interactions/four_particle && rm -f *.o *.so *.c *.cpp

test:
	pip install -U nose
	nosetests tests/unit