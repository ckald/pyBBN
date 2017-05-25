all: clean compile test

compile:
	cd interactions/four_particle && python3 setup.py build_ext --inplace


clean:
	cd interactions/four_particle && rm -f *.o *.so *.c *.cpp

test:
	pip3 install -U nose
	python3 -m "nose" tests/unit