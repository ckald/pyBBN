COMPILER ?= c++
PYTHON ?= python3.5
all: compile_4p compile_3p test

compile_3p:
	cd interactions/three_particle/cpp && make all

compile_4p:
	cd interactions/four_particle/cpp && make all

clean:
	cd interactions/four_particle/cpp && make clean
	cd interactions/three_particle/cpp && make clean

test:
	$(PYTHON) -m "nose" tests/unit
