import matplotlib.pyplot as plt
import numpy
from particle import Particle
from interaction import Interaction
from common import Distributions, PARAMS, CONST, UNITS, STATISTICS
import numericalunits as nu

Particles = []
Interactions = []

PARAMS.T = 1 * nu.MeV

neutron = Particle(name='Neutron',
                   statistics=STATISTICS.FERMION,
                   mass=0.939 * nu.GeV,
                   decoupling_temperature=2 * nu.MeV)
Particles.append(neutron)

proton = Particle(name='Proton', statistics=STATISTICS.FERMION, mass=0.938 * nu.GeV)
Particles.append(proton)

neutrino = Particle(name='Neutrino', statistics=STATISTICS.FERMION, dof=4)
Particles.append(neutrino)

electron = Particle(name='Electron', mass=0.511 * nu.MeV, statistics=STATISTICS.FERMION, dof=4)
Particles.append(electron)

neutron_decay = Interaction(
    in_particles=[neutron],
    out_particles=[proton, neutrino, electron],
    decoupling_temperature=0 * nu.MeV,

    K1=64 * CONST.G_F**2,
    K2=0
)
Interactions.append(neutron_decay)

figure, plots = plt.subplots(2, 2)
figure.subplots_adjust(hspace=0.5)
plt.ion()

plots[0][0].set_title("Temperature")

plots[0][1].set_title("Scale factor")


def evolve(t_f=PARAMS.t_f, t_i=PARAMS.t_i, dt=PARAMS.dt):
    PARAMS.dt = dt
    time_steps = int((t_f - t_i) / dt)

    Ts = numpy.zeros(time_steps)
    Ts[0] = PARAMS.T
    As = numpy.zeros(time_steps)
    As[0] = PARAMS.a
    ts = numpy.zeros(time_steps)
    ts[0] = t_i
    rhos = numpy.zeros(time_steps)
    rhos[0] = 0
    for particle in Particles:
        rhos[0] += particle.energy_density()

    for t in xrange(time_steps - 1):
        ts[t+1] = ts[t] + dt
        rhos[t+1] = 0
        for particle in Particles:
            rhos[t+1] += particle.energy_density()

        Ts[t + 1] = Ts[t] * (
            1. - numpy.sqrt(8/3*numpy.pi*CONST.G * rhos[t+1]) * dt
        )
        PARAMS.T = Ts[t + 1]

        for particle in Particles:
            particle.update()

        for interaction in Interactions:
            interaction.calculate()

        As[t+1] = As[t] * (2. - Ts[t+1] / Ts[t])
        PARAMS.a = As[t + 1]

        if t % (1. * time_steps/100) == 0:
            print ts[t+1] / UNITS.s, "\t", Ts[t+1] / nu.MeV, "\t", As[t+1]

    return ts, Ts, As, rhos

ts, Ts, As, rhos = evolve(t_i=0 * UNITS.s, t_f=8 * UNITS.s, dt=0.08 * UNITS.s)

plots[0][0].plot(ts / UNITS.s, Ts / nu.MeV)
plots[0][0].set_yscale("log")

plots[0][1].plot(ts / UNITS.s, As)

plots[1][0].plot(ts / UNITS.s, Ts / Ts[0] * As / As[0])
plots[1][0].set_title("T * a")

plots[1][1].plot(ts / UNITS.s, rhos / nu.eV**4)
plots[1][1].set_title("Total energy density")
plots[1][1].set_yscale("log")
figure.show()

plt.figure(2)
plt.ion()
plt.plot(neutron.distribution_noneq)
plt.plot(neutron.distribution_eq)
plt.plot(neutron.distribution)
plt.show()

for particle in Particles:
    print particle

raw_input("...")
