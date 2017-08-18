import os
import json
import numpy

# from interactions.four_particle.integral import D1, D2, D3
from interactions.four_particle.cpp.integral import D1, D2, D3

from common import LogSpacedGrid

cwd = os.path.split(__file__)[0]


GRID = LogSpacedGrid(MOMENTUM_SAMPLES=20)

def grid_iterator(dimension=4, grid=GRID):
    value = numpy.zeros(dimension, dtype=numpy.int)

    dim = dimension - 1
    while dim >= 0:
        # print(value)
        value[dim] += 1
        if value[dim] >= grid.MOMENTUM_SAMPLES:
            while dim >= 0 and value[dim] >= grid.MOMENTUM_SAMPLES:
                value[dim] = 0
                value[dim-1] += 1
                dim -= 1
            if dim > -1:
                dim = dimension - 1
        yield value, map(lambda i: grid.TEMPLATE[i], value)


grid_json = json.dumps(GRID.TEMPLATE.tolist())
with open(os.path.join(cwd, 'grid.json'), 'w') as f:
    f.write(grid_json)


print("Computing D1...")
d1_table = numpy.ndarray(shape=(GRID.MOMENTUM_SAMPLES, GRID.MOMENTUM_SAMPLES,
                                GRID.MOMENTUM_SAMPLES, GRID.MOMENTUM_SAMPLES))

for (i, j, k, l), momenta in grid_iterator():
    d1_table[i, j, k, l] = D1(*momenta)

with open(os.path.join(cwd, 'D1.json'), 'w') as f:
    f.write(json.dumps(d1_table.tolist()))
print("Done")

print("Computing D2...")
d2_table = numpy.ndarray(shape=(GRID.MOMENTUM_SAMPLES, GRID.MOMENTUM_SAMPLES,
                                GRID.MOMENTUM_SAMPLES, GRID.MOMENTUM_SAMPLES))

for (i, j, k, l), momenta in grid_iterator():
    d2_table[i, j, k, l] = D2(*momenta)

with open(os.path.join(cwd, 'D2.json'), 'w') as f:
    f.write(json.dumps(d2_table.tolist()))
print("Done")

print("Computing D3...")
d3_table = numpy.ndarray(shape=(GRID.MOMENTUM_SAMPLES, GRID.MOMENTUM_SAMPLES,
                                GRID.MOMENTUM_SAMPLES, GRID.MOMENTUM_SAMPLES))

for (i, j, k, l), momenta in grid_iterator():
    d3_table[i, j, k, l] = D3(*momenta)

with open(os.path.join(cwd, 'D3.json'), 'w') as f:
    f.write(json.dumps(d3_table.tolist()))
print("Done")
