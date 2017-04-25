import os
import json
import numpy

from interactions.four_particle.integral import D1, D2, D3, Dpy
from common import LogSpacedGrid, UNITS, Params
# process ν_μ + e ⟶ ν_μ + e
m = (0,0.5,0,0.5)
K213 = 0.00000000000000000000104249
K114 = 0.000000000000000000000896662
K112 = 0.00000000000000000000121202
sides = (-1, -1, 1, 1)

cwd = os.path.split(__file__)[0]


GRID = LogSpacedGrid(MOMENTUM_SAMPLES=20)

def grid_iterator(dimension=4, grid=GRID):
    value = numpy.zeros(dimension, dtype=numpy.int)

    dim = dimension - 1
    while dim >= 0:
        # print value
        value[dim] += 1
        if value[dim] >= grid.MOMENTUM_SAMPLES:
            while dim >= 0 and value[dim] >= grid.MOMENTUM_SAMPLES:
                value[dim] = 0
                value[dim-1] += 1
                dim -= 1
            if dim > -1:
                dim = dimension - 1
        yield value, map(lambda i: grid.TEMPLATE[i], value)

print "Computing Dpy112..."
dpy112_table = numpy.ndarray(shape=(GRID.MOMENTUM_SAMPLES, GRID.MOMENTUM_SAMPLES,
                                GRID.MOMENTUM_SAMPLES, GRID.MOMENTUM_SAMPLES))
order112 = (0, 1, 2, 3)
for (i, j, k, l), momenta in grid_iterator():
    p= (i,j,k,l)
    E= (i, numpy.sqrt(j**2 + 0.5**2), k, numpy.sqrt(l**2 + 0.5**2))
    dpy112_table[i, j, k, l] = Dpy(p, E, m, K112, 0, order112, sides)

with open(os.path.join(cwd, 'Dpy112.json'), 'w') as f:
    f.write(json.dumps(dpy112_table.tolist()))
print "Done"

print "Computing Dpy114..."
dpy114_table = numpy.ndarray(shape=(GRID.MOMENTUM_SAMPLES, GRID.MOMENTUM_SAMPLES,
                                GRID.MOMENTUM_SAMPLES, GRID.MOMENTUM_SAMPLES))
order114 = (0, 3, 1, 2)
for (i, j, k, l), momenta in grid_iterator():
    p= (i,j,k,l)
    E= (i, numpy.sqrt(j**2 + 0.5**2), k, numpy.sqrt(l**2 + 0.5**2))
    dpy114_table[i, j, k, l] = Dpy(p, E, m, K114, 0, order114, sides)

with open(os.path.join(cwd, 'Dpy114.json'), 'w') as f:
    f.write(json.dumps(dpy114_table.tolist()))
print "Done"

print "Computing Dpy213..."
dpy213_table = numpy.ndarray(shape=(GRID.MOMENTUM_SAMPLES, GRID.MOMENTUM_SAMPLES,
                                GRID.MOMENTUM_SAMPLES, GRID.MOMENTUM_SAMPLES))
order213 = (1, 3, 0, 2)
for (i, j, k, l), momenta in grid_iterator():
    p= (i,j,k,l)
    E= (i, numpy.sqrt(j**2 + 0.5**2), k, numpy.sqrt(l**2 + 0.5**2))
    dpy213_table[i, j, k, l] = Dpy(p, E, m, K213, 0, order213, sides)

with open(os.path.join(cwd, 'Dpy213.json'), 'w') as f:
    f.write(json.dumps(dpy213_table.tolist()))
print "Done"

print "Computing Dpy..."
dpy_table = numpy.ndarray(shape=(GRID.MOMENTUM_SAMPLES, GRID.MOMENTUM_SAMPLES,
                                GRID.MOMENTUM_SAMPLES, GRID.MOMENTUM_SAMPLES))
order213 = (1, 3, 0, 2)
for (i, j, k, l), momenta in grid_iterator():
    dpy_table[i, j, k, l] = dpy112_table[i, j, k, l] + dpy114_table[i, j, k, l] + dpy213_table[i, j, k, l]

with open(os.path.join(cwd, 'Dpy.json'), 'w') as f:
    f.write(json.dumps(dpy_table.tolist()))
print "Done"