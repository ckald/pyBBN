import os
import json
import numpy

from interactions import ds
from common import GRID

cwd = os.path.split(__file__)[0]


def grid_iterator(dimension=4):
    value = numpy.zeros(dimension)

    dim = dimension - 1
    while dim >= 0:
        # print value
        value[dim] += 1
        if value[dim] >= GRID.MOMENTUM_SAMPLES:
            while dim >= 0 and value[dim] >= GRID.MOMENTUM_SAMPLES:
                value[dim] = 0
                value[dim-1] += 1
                dim -= 1
            if dim > -1:
                dim = dimension - 1
        yield value, map(lambda i: GRID.TEMPLATE[i], value)


grid_json = json.dumps(GRID.TEMPLATE.tolist())
with open(os.path.join(cwd, 'grid.json'), 'w') as f:
    f.write(grid_json)


d1_table = numpy.ndarray(shape=(GRID.MOMENTUM_SAMPLES, GRID.MOMENTUM_SAMPLES,
                                GRID.MOMENTUM_SAMPLES, GRID.MOMENTUM_SAMPLES))

for (i, j, k, l), momenta in grid_iterator():
    d1_table[i, j, k, l] = ds.D1(*momenta)

with open(os.path.join(cwd, 'D1.json'), 'w') as f:
    f.write(json.dumps(d1_table.tolist()))


d2_table = numpy.ndarray(shape=(GRID.MOMENTUM_SAMPLES, GRID.MOMENTUM_SAMPLES,
                                GRID.MOMENTUM_SAMPLES, GRID.MOMENTUM_SAMPLES))

for (i, j, k, l), momenta in grid_iterator():
    d2_table[i, j, k, l] = ds.D2(*momenta)

with open(os.path.join(cwd, 'D2.json'), 'w') as f:
    f.write(json.dumps(d2_table.tolist()))

d3_table = numpy.ndarray(shape=(GRID.MOMENTUM_SAMPLES, GRID.MOMENTUM_SAMPLES,
                                GRID.MOMENTUM_SAMPLES, GRID.MOMENTUM_SAMPLES))

for (i, j, k, l), momenta in grid_iterator():
    d3_table[i, j, k, l] = ds.D3(*momenta)

with open(os.path.join(cwd, 'D3.json'), 'w') as f:
    f.write(json.dumps(d3_table.tolist()))
