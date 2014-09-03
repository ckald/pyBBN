import matplotlib.pyplot as plt
import numpy
from particle import Particle
from common import PARAMS, CONST, UNITS, STATISTICS, benchmark
import numericalunits as nu

PARAMS.T = 2. * 1e2 * nu.MeV

Particles = []
photon = Particle(name='Photon',
                  statistics=STATISTICS.BOSON)
Particles.append(photon)

neutron = Particle(name='Neutron',
                   statistics=STATISTICS.FERMION,
                   mass=0.939 * nu.GeV)
Particles.append(neutron)

proton = Particle(name='Proton',
                  statistics=STATISTICS.FERMION,
                  mass=0.938 * nu.GeV)
Particles.append(proton)

neutrino = Particle(name='Neutrino',
                    statistics=STATISTICS.FERMION,
                    dof=4,
                    # decoupling_temperature=2
                    )
Particles.append(neutrino)

electron = Particle(name='Electron',
                    mass=0.511 * nu.MeV,
                    statistics=STATISTICS.FERMION,
                    dof=4)
Particles.append(electron)


figure, plots = plt.subplots(2, 2)
figure.subplots_adjust(hspace=0.5, wspace=0.5)
plt.ion()

plots[0][0].set_title("Temperature")
plots[0][0].set_xlabel("time, s")
plots[0][0].set_ylabel("T, MeV")
plots[0][1].set_title("Scale factor")
plots[0][1].set_xlabel("time, s")
plots[0][1].set_ylabel("a, 1")


def evolve(t_f=PARAMS.t_f, dt=PARAMS.dt):
    time_steps = int((t_f - PARAMS.t) / dt)

    Ts = numpy.zeros(time_steps)
    Ts[0] = PARAMS.T
    As = numpy.zeros(time_steps)
    As[0] = PARAMS.a
    ts = numpy.zeros(time_steps)
    ts[0] = PARAMS.t
    rhos = numpy.zeros(time_steps)
    rhos[0] = 0
    for particle in Particles:
        rhos[0] += particle.energy_density()

    for t in xrange(time_steps - 1):
        ts[t+1] = ts[t] + dt

        rho = 0  # total energy density
        p = 0    # total pressure

        # equilibrium energy density derivative divided by  1/T * dT/dt
        energy_density_rate_eq = 0
        # nonequilibrium energy density derivative including collision integrals
        energy_density_rate_neq = 0

        for particle in Particles:
                rho += particle.energy_density()
                p += particle.pressure()

                if particle.in_equilibrium:
                    energy_density_rate_eq += particle.energy_density_rate()
                else:
                    energy_density_rate_neq += particle.energy_density_rate()

        H = numpy.sqrt(8./3.*numpy.pi * CONST.G * rho)

        A = (3. * H * (rho + p) + energy_density_rate_neq) / energy_density_rate_eq

        Ts[t+1] = Ts[t] * (1. - A * dt)
        PARAMS.T = Ts[t + 1]

        for particle in Particles:
            particle.update()
            if not particle.in_equilibrium:
                particle._distribution *= (1. - dt * 3. * H)

        As[t+1] = As[t] * (2. - Ts[t+1] / Ts[t])
        PARAMS.a = As[t + 1]

        rhos[t+1] = rho

        if t % (time_steps/100) == 0:
            print ts[t+1] / UNITS.s, "\t", Ts[t+1] / nu.MeV, "\t", As[t+1]

    PARAMS.t = t_f

    return ts, Ts, As, rhos

ts, Ts, As, rhos = evolve(t_f=1e-1 * UNITS.s, dt=1e-6 * UNITS.s)
ts1, Ts1, As1, rhos1 = evolve(t_f=10 * UNITS.s, dt=1e-4 * UNITS.s)
ts2, Ts2, As2, rhos2 = evolve(t_f=100 * UNITS.s, dt=1e-3 * UNITS.s)

ts = numpy.append(ts, ts1)
Ts = numpy.append(Ts, Ts1)
As = numpy.append(As, As1)
rhos = numpy.append(rhos, rhos1)

ts = numpy.append(ts, ts2)
Ts = numpy.append(Ts, Ts2)
As = numpy.append(As, As2)
rhos = numpy.append(rhos, rhos2)

plots[0][0].plot(ts / UNITS.s, Ts / nu.MeV)
plots[0][0].set_yscale("log")

plots[0][1].plot(ts / UNITS.s, As)

plots[1][0].plot(ts / UNITS.s, Ts / Ts[0] * As / As[0])
plots[1][0].set_title("T * a")
plots[1][0].set_xlabel("time ,s")
plots[1][0].set_ylabel("T * a, 1")

plots[1][1].plot(ts / UNITS.s, rhos / nu.eV**4)
plots[1][1].set_title("Total energy density")
plots[1][1].set_yscale("log")
plots[1][1].set_xlabel("time, s")
plots[1][1].set_ylabel("rho, eV**4")
figure.show()

for particle in Particles:
    print particle

raw_input("...")
