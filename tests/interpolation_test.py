from interaction import Interaction
from particles import Particle, STATISTICS
from evolution import Universe
from common import CONST, UNITS, PARAMS, GRID


PARAMS.T_initial = 3 * UNITS.MeV
PARAMS.T_final = 0.075 * UNITS.MeV
PARAMS.dx = 1e-1 * UNITS.MeV
PARAMS.infer()


Particles = []
Interactions = []
photon = Particle(name='Photon',
                  statistics=STATISTICS.BOSON)
Particles.append(photon)

# neutron = Particle(name='Neutron',
#                    statistics=STATISTICS.FERMION,
#                    mass=0.939 * UNITS.GeV,
#                    decoupling_temperature=1.3 * UNITS.MeV)
# Particles.append(neutron)

# proton = Particle(name='Proton',
#                   statistics=STATISTICS.FERMION,
#                   mass=0.938 * UNITS.GeV,
#                   )
# Particles.append(proton)

neutrino = Particle(name='Neutrino',
                    statistics=STATISTICS.FERMION,
                    dof=4,
                    decoupling_temperature=3 * UNITS.MeV
                    )
Particles.append(neutrino)

# electron = Particle(name='Electron',
#                     mass=0.511 * UNITS.MeV,
#                     statistics=STATISTICS.FERMION,
#                     dof=4)
# Particles.append(electron)

neutrino_scattering = Interaction(
    in_particles=[neutrino, neutrino],
    out_particles=[neutrino, neutrino],
    decoupling_temperature=0 * UNITS.MeV,
    K1=128 * CONST.G_F**2,
    K2=0.,
    order=(0, 1, 2, 3),
    symmetry_factor=0.5
)
Interactions.append(neutrino_scattering)

import numpy
neutrino._distribution += numpy.vectorize(lambda x: 0.1 * numpy.exp(-(x-5)**2),
                                          otypes=[numpy.float_])(GRID.TEMPLATE)


# def check(p=[]):
#     return neutrino_scattering.F_B(in_p=p[:2], out_p=p[2:]) * (1 - neutrino.distribution(p[0])) \
#         - neutrino_scattering.F_A(in_p=p[:2], out_p=p[2:]) * neutrino.distribution(p[0])


# # print [neutrino.distribution(p) - neutrino._distribution_interpolation(p)
# #        for p in numpy.linspace(GRID.MIN_MOMENTUM, GRID.MAX_MOMENTUM,
# #                                num=200, endpoint=True)]

import matplotlib.pyplot as plt

x = numpy.linspace(GRID.MIN_MOMENTUM, GRID.MAX_MOMENTUM, num=10000, endpoint=True)
y = numpy.vectorize(neutrino.distribution)(x)
z = numpy.vectorize(neutrino.distribution_function)(x / PARAMS.aT)
plt.plot(x, y)
plt.plot(x, z)
plt.show()

# res = []
# for p0 in GRID.TEMPLATE:
#     for p1 in GRID.TEMPLATE:
#         for p2 in GRID.TEMPLATE:
#             p3 = p0 + p1 - p2
#             if p3 >= 0 and p3 <= GRID.MAX_MOMENTUM:
#                 # val = check(p=[p0, p1, p2, p3])
#                 # f0 = neutrino.distribution(p0)
#                 # f1 = neutrino.distribution(p1)
#                 # f2 = neutrino.distribution(p2)
#                 # f3 = neutrino.distribution(p3)
#                 # val = f2 * f3 * (1 - f0) * (1 - f1) - f0 * f1 * (1 - f2) * (1 - f3)
#                 res.append(val)

from ds import D1, D2, D3
from plotting import plot_integrand

plot_integrand(lambda (p1, p2): D1(p1, p2, 1, 1), 'D1')
plot_integrand(lambda (p1, p2): D2(p1, p2, 1, 1), 'D2')
plot_integrand(lambda (p1, p2): D3(p1, p2, 1, 1), 'D3')

raw_input("...")
