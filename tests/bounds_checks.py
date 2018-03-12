import os
import numpy as np
import json

from common import UNITS, LogSpacedGrid

cwd = os.path.split(__file__)[0]

GRID = LogSpacedGrid(MOMENTUM_SAMPLES=30, MAX_MOMENTUM=20)

masses = np.linspace(0, 20, 4)

def grid_iterator(dimension=None, grid=GRID):
    value = np.zeros(dimension, dtype=np.int)

    dim = dimension - 1
    while dim >= 0:
        yield list(map(lambda i: grid.TEMPLATE[i], value))
        value[dim] += 1
        if value[dim] >= grid.MOMENTUM_SAMPLES:
            while dim >= 0 and value[dim] >= grid.MOMENTUM_SAMPLES:
                value[dim] = 0
                value[dim-1] += 1
                dim -= 1
            if dim > -1:
                dim = dimension - 1


def mass_iterator(dimension=None, mass=masses):
    value = np.zeros(dimension, dtype=np.int)

    dim = dimension - 1
    while dim >= 0:
        yield list(map(lambda i: mass[i], value))
        value[dim] += 1
        if value[dim] >= len(mass):
            while dim >= 0 and value[dim] >= len(mass):
                value[dim] = 0
                value[dim-1] += 1
                dim -= 1
            if dim > -1:
                dim = dimension - 1


def energy(m, p):
    return np.sqrt(m**2 + p**2)


def three_particle_filter(E, p, m):
    q1 = p[0]
    q2 = p[1]
    q3 = p[2]

    return (E[2] >= m[2] and (np.sign(q1+q2-q3) + np.sign(q1-q2+q3) - np.sign(q1-q2-q3) - np.sign(q1+q2+q3)))


def four_particle_filter(E, p, m):
    q1 = p[0]
    q2 = p[1]
    q3 = p[2]
    q4 = p[3]

    if q1 < q2:
        q1, q2 = q2, q1
    if q3 < q4:
        q3, q4 = q4, q3

    return (E[3] >= m[3] and q1 <= q2 + q3 + q4 and q3 <= q1 + q2 + q4)


def three_particle_energy_momentum(mom=None, mass=None):
    if mass[0] == 0:
        return False

    E2 = energy(mass[0], mom[0]) - energy(mass[1], mom[1])

    if E2 >= mass[2]:
        p2 = np.sqrt(E2**2 - mass[2]**2)
    else:
        return False

    Energy = np.array([energy(mass[0], mom[0]), energy(mass[1], mom[1]), E2])
    Momenta = np.append(mom, p2)

    return three_particle_filter(Energy, Momenta, mass)


def four_particle_energy_momentum(mom=None, mass=None, s1=1, s2=-1):
    E3 = energy(mass[0], mom[0]) + s1 * energy(mass[1], mom[1]) + s2 * energy(mass[2], mom[2])
    if E3 >= mass[3]:
        p3 = np.sqrt(E3**2 - mass[3]**2)
    if s1 == -1 and s2 == -1 and mass[0] == 0:
        return False
    else:
        return False

    Energy = np.array([energy(mass[0], mom[0]), energy(mass[1], mom[1]), energy(mass[2], mom[2]), E3])
    Momenta = np.append(mom, p3)

    return four_particle_filter(Energy, Momenta, mass)


def p2_min_dec(i=0, j=1, k=2, l=3, mass=None):
    return np.sqrt(
                (((mass[i] - mass[j])**2 - mass[k]**2 - mass[l]**2)**2
                - 4 * (mass[k] * mass[l])**2
                ) / (4 * (mass[i] - mass[j])**2)
            )


def four_particle_decay_check(mass=None, p0=None, p1=None, p2=None):
    upper_bound = np.sqrt((energy(mass[0], p0) - energy(mass[1], p1) - mass[3])**2 - mass[2]**2)

    if p0 == 0:
        lower_bound = np.maximum(p2_min_dec(0, 1, 2, 3, mass=mass)
                    - p1 * (p2_min_dec(0, 1, 2, 3, mass=mass) / p2_min_dec(0, 2, 1, 3, mass=mass)), 0)
    else:
        lower_bound = 0

    if p2 > upper_bound or p2 < lower_bound:
        return False
    return True


def four_particle_scattering_lower_bound(mass=None, p0=None, p1=None, p2=None):

    if p0 != 0:
        return 0.

    if mass[0] != 0:
        temp1 = mass[0] + energy(mass[1], p1)
        temp2 = mass[3]**2 + p1**2
        temp3 = ((temp2 - temp1**2 - mass[2]**2) * p1
                + np.sqrt(
                        temp1**2 * (temp1**4 + temp2**2
                        + (4 * p1**2 - 2 * temp2 + mass[2]**2) * mass[2]**2
                        - 2 * temp1**2 * (temp2 + mass[2]**2))
                    )
                ) / (2 * (temp1**2 - p1**2))

        return np.maximum(temp3, 0.)

    if mass[1] != 0:
        temp = (p1**2 + mass[1]**2) * (
                mass[1]**4 + (mass[2]**2 - mass[3]**2)**2
                - 2 * mass[1]**2 * (mass[2]**2 + mass[3]**2))

        if temp < 0:
            temp = 0

        temp1 = (-1 * p1 * (mass[1]**2 + mass[2]**2 - mass[3]**2)
                    + np.sqrt(temp)
                ) / (2 * mass[1]**2)

        temp2 = -temp1

        return np.maximum(temp1, temp2)

    return 0.


def four_particle_scattering_upper_bound(mass=None, p0=None, p1=None, p2=None):

    if p0 != 0:
        print(mass[0], mass[1], mass[2], mass[3], np.sqrt(
                (energy(mass[0], p0)
                 + energy(mass[1], p1)
                 - mass[3])**2
                - mass[2]**2
                ))
        return np.sqrt(
                (energy(mass[0], p0)
                 + energy(mass[1], p1)
                 - mass[3])**2
                - mass[2]**2
                )

    if mass[0] != 0:
        temp1 = mass[0] + energy(mass[1], p1)
        temp2 = mass[3]**2 + p1**2
        temp3 = ((-temp2 + temp1**2 + mass[2]**2) * p1
                + np.sqrt(
                        temp1**2 * (temp1**4 + temp2**2
                        + (4 * p1**2 - 2 * temp2 + mass[2]**2) * mass[2]**2
                        - 2 * temp1**2 * (temp2 + mass[2]**2))
                    )
                ) / (2 * (temp1**2 - p1**2))

        return np.maximum(temp3, 0.)

    if mass[1] != 0:
        temp = (p1**2 + mass[1]**2) * (
                mass[1]**4 + (mass[2]**2 - mass[3]**2)**2
                - 2 * mass[1]**2 * (mass[2]**2 + mass[3]**2)
                )

        if temp < 0:
            temp = 0

        return (p1 * (mass[1]**2 + mass[2]**2 - mass[3]**2)
                + np.sqrt(temp)
                ) / (2 * mass[1]**2)

    return p1


def four_particle_scattering_check(mass=None, p0=None, p1=None, p2=None):
    upper_bound = four_particle_scattering_upper_bound(mass, p0, p1, p2)
    lower_bound = four_particle_scattering_lower_bound(mass, p0, p1, p2)

    if p2 > upper_bound or p2 < lower_bound:
        return False
    return True


def three_particle_decay_bounds(mass=None, p0=None, sign1=None, sign2=None):
    temp1 = (p0**2 + mass[0]**2) * (mass[1]**4 + (mass[2]**2 - mass[0]**2)**2
        - 2 * mass[1]**2 * (mass[2]**2 + mass[0]**2)
        )

    if temp1 < 0:
        temp1 = 0

    temp2 = p0 * (mass[2]**2 - mass[1]**2 - mass[0]**2)
    return  (sign1 * temp2 + sign2 * np.sqrt(temp1)) / (2 * mass[0]**2)


def three_particle_decay_check(mass=None, p0=None, p1=None):
    if p0 == 0:
        return True

    temp = three_particle_decay_bounds(mass, p0, 1, 1)
    lower_bound = np.maximum(temp, -temp)
    upper_bound = three_particle_decay_bounds(mass, p0, -1, 1)

    if p1 > upper_bound or p1 < lower_bound:
        return False
    return True


print("Checking bounds for four-particle decay and scattering...")
for mass in mass_iterator(dimension=4):
    for i, j, k in grid_iterator(dimension=3):
        if four_particle_energy_momentum(mom=[k,i,j], mass=mass, s1=-1, s2=-1):
            if not four_particle_decay_check(mass[::-1], k, i, j):
                raise AssertionError("Sample point p0 = {p0:e}\t p1 = {p1:e}\t p2 = {p2:e}\t out of approximate bounds for decays".format(p0=k, p1=i, p2=j))

        if four_particle_energy_momentum(mom=[k,i,j], mass=mass, s1=1, s2=-1):
            if not four_particle_scattering_check(mass[::-1], k, i, j):
                raise AssertionError("Sample point\t p0 = {p0:e}\t p1 = {p1:e}\t p2 = {p2:e}\t out of approximate bounds for scatterings".format(p0=k, p1=i, p2=j))

print("Bounds for four-particle decay and scatterings are OK", "\n")


print("Checking bounds for three-particle decay...")
for mass in mass_iterator(dimension=3):
    for i, j in grid_iterator(dimension=2):
        if three_particle_energy_momentum(mom=[j, i], mass=mass[::-1]):
            if not three_particle_decay_check(mass[::-1], j, i):
                raise AssertionError("Sample point\t p0 = {p0:e}\t p1 = {p1:e}\t out of approximate bounds for decays".format(p0=j, p1=i))

print("Bounds for three-particle decay are OK")