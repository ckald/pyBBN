from library import StandardModelParticles as SMP
from interaction import Interaction, M
from particles import Particle
from evolution import Universe
from common import CONST, UNITS, Params, GRID


params = Params(T_initial=2 * UNITS.MeV,
                T_final=0.075 * UNITS.MeV,
                dx=1e-1 * UNITS.MeV)

photon = Particle(params=params, **SMP.photon)
neutrino = Particle(params=params, **SMP.neutrino_e)

neutrino_scattering = Interaction(
    in_particles=[neutrino, neutrino],
    out_particles=[neutrino, neutrino],
    decoupling_temperature=0 * UNITS.MeV,
    Ms=[M(K1=64 * CONST.G_F**2, order=(0, 1, 2, 3))]
)

universe = Universe(params=params, plotting=False)
universe.particles += [photon, neutrino]
universe.interactions += [neutrino_scattering]


import numpy
addition = numpy.vectorize(lambda x: 0.1 * numpy.exp(-(x/UNITS.MeV-5)**2),
                           otypes=[numpy.float_])
# neutrino._distribution += addition(GRID.TEMPLATE)

# def check(p=[]):
#     return neutrino_scattering.F_B(in_p=p[:2], out_p=p[2:]) * (1 - neutrino.distribution(p[0])) \
#         - neutrino_scattering.F_A(in_p=p[:2], out_p=p[2:]) * neutrino.distribution(p[0])


# # print [neutrino.distribution(p) - neutrino._distribution_interpolation(p)
# #        for p in numpy.linspace(GRID.MIN_MOMENTUM, GRID.MAX_MOMENTUM,
# #                                num=200, endpoint=True)]

import matplotlib.pyplot as plt

x = numpy.linspace(GRID.MIN_MOMENTUM, GRID.MAX_MOMENTUM*2,
                   num=GRID.MOMENTUM_SAMPLES*10, endpoint=True)
y = numpy.vectorize(neutrino.distribution)(x)
z = numpy.vectorize(neutrino.equilibrium_distribution_function)(x / universe.params.aT)
w = z + addition(x)
plt.plot(x, y)
# plt.plot(x, z)
# plt.plot(x, w)
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

plot_integrand(lambda (p1, p2): D1(p1, p2, 1, 1), 'D1', 0)
plot_integrand(lambda (p1, p2): D2(p1, p2, 1, 1), 'D2', 0)
plot_integrand(lambda (p1, p2): D3(p1, p2, 1, 1), 'D3', 0)

raw_input("...")
