# -*- coding: utf-8 -*-

import os
import sys
import time
import numpy
import shutil
from datetime import timedelta

import environment
from common import CONST, UNITS, Params, utils
from common.integrators import adams_bashforth_correction, MAX_ADAMS_BASHFORTH_ORDER

import kawano


class Universe(object):

    """ ## Universe
        The master object that governs the calculation. """

    # System state is rendered to the evolution.txt file each `export_freq` steps
    export_freq = 100
    log_throttler = None
    clock_start = None

    particles = None
    interactions = None

    kawano = None
    kawano_log = None

    oscillations = None

    step_monitor = None

    data = utils.DynamicRecArray([
        ['aT', 'MeV', UNITS.MeV],
        ['T', 'MeV', UNITS.MeV],
        ['a', None, 1],
        ['x', 'MeV', UNITS.MeV],
        ['t', 's', UNITS.s],
        ['rho', 'MeV^4', UNITS.MeV**4],
        ['N_eff', None, 1],
        ['fraction', None, 1],
        ['S', 'MeV^3', UNITS.MeV**3]
    ])

    def __init__(self, folder=None, params=None, max_log_rate=2):
        """
        :param folder: Log file path (current `datetime` by default)
        """

        self.particles = []
        self.interactions = []

        self.clock_start = time.time()
        self.log_throttler = utils.Throttler(max_log_rate)

        self.params = params
        if not self.params:
            self.params = Params()

        self.folder = folder
        if self.folder:
            if os.path.exists(folder):
                shutil.rmtree(folder)
            self.init_log(folder=folder)

        self.fraction = 0

        self.step = 1

    def init_kawano(self, datafile='s4.dat', **kwargs):
        kawano.init_kawano(**kwargs)
        if self.folder:
            self.kawano_log = open(os.path.join(self.folder, datafile), 'w')
            self.kawano_log.write("\t".join([col[0] for col in kawano.heading]) + "\n")
        self.kawano = kawano
        self.kawano_data = utils.DynamicRecArray(self.kawano.heading)

    def oscillation_parameters(self):

        if self.oscillation_matter:
            MSW_12 = CONST.MSW_constant * self.oscillation_particles[0].grid.TEMPLATE**2 * self.params.T**4 / CONST.delta_m12_sq / self.params.a**2
            MSW_13 = CONST.MSW_constant * self.oscillation_particles[0].grid.TEMPLATE**2 * self.params.T**4 / CONST.delta_m13_sq / self.params.a**2
            if not environment.get("NORMAL_HIERARCHY_NEUTRINOS"):
                MSW_13 *= -1

        else:
            MSW_12 = 0.
            MSW_13 = 0.

        pattern = self.pattern_function(MSW_12, MSW_13, self.oscillation_matter)
        self.oscillations = (pattern, self.oscillation_particles)

    def init_oscillations(self, pattern_function, particles, matter_effects=True):
        self.pattern_function = pattern_function
        self.oscillation_particles = particles
        self.oscillation_matter = matter_effects
        self.oscillation_parameters()

    def evolve(self, T_final, export=True, init_time=True):
        """
        ## Main computing routine

        Modeling is carried in steps of scale factor, starting from the initial conditions defined\
        by a single parameter: the initial temperature.

        Initial temperature (e.g. 10 MeV) corresponds to a point in the history of the Universe \
        long before then BBN. Then most particle species are in the thermodynamical equilibrium.

        """
        T_initial = self.params.T

        print("\n\n" + "#"*32 + " Initial states " + "#"*32 + "\n")
        for particle in self.particles:
            print(particle)
        print("\n\n" + "#"*34 + " Log output " + "#"*34 + "\n")

        for interaction in self.interactions:
            if interaction.integrals:
                print(interaction)
        print("\n")

        # TODO: test if changing updating particles beforehand changes the computed time
        if init_time:
            self.params.init_time(self.total_energy_density())

        if self.params.rho is None:
#            self.update_particles()
            self.params.update(self.total_energy_density(), self.total_entropy())
        self.save_params()

        while self.params.T > T_final:
            try:
                self.log()
                self.make_step()
                self.save()
                self.step += 1
                if self.folder and self.step % self.export_freq == 0:
                    with open(os.path.join(self.folder, "evolution.txt"), "wb") as f:
                        self.data.savetxt(f)
                    if self.kawano:
                        with open(os.path.join(self.folder, "kawano.txt"), "wb") as f:
                            self.kawano_data.savetxt(f)
            except KeyboardInterrupt:
                print("\nKeyboard interrupt!")
                sys.exit(1)
                break

        if not (self.params.T > 0):
            print("\n(T < 0): suspect numerical instability")
            sys.exit(1)

        self.log()
        if export:
            self.export()

        return self.data

    def export(self):
        print("\n\n" + "#"*33 + " Final states " + "#"*33 + "\n")
        for particle in self.particles:
            print(particle)
        print("\n")

        if self.folder:
            if self.kawano:

                self.kawano_log.close()
                print(kawano.run(self.folder))

                with open(os.path.join(self.folder, "kawano.txt"), "wb") as f:
                    self.kawano_data.savetxt(f)

            print("Execution log saved to file {}".format(self.logfile))

            with open(os.path.join(self.folder, "evolution.txt"), "wb") as f:
                self.data.savetxt(f)

    def make_step(self):
        self.integrand(self.params.x, self.params.aT)

        if self.step_monitor:
            self.step_monitor(self)

        if environment.get('ADAMS_BASHFORTH_TEMPERATURE_CORRECTION'):
            fs = (list(self.data['fraction'][-MAX_ADAMS_BASHFORTH_ORDER:]) + [self.fraction])

            self.params.aT += adams_bashforth_correction(fs=fs, h=self.params.h)
        else:
            self.params.aT += self.fraction * self.params.h

        self.params.x += self.params.dx
        self.params.update(self.total_energy_density(), self.total_entropy())

        self.log_throttler.update()

    def add_particles(self, particles):
        for particle in particles:
            particle.set_params(self.params)

        self.particles += particles

        self.particles = utils.particle_orderer(self.particles)

    def update_particles(self):
        """ ### 1. Update particles state
            Update particle species distribution functions, check for regime switching,\
            update precalculated variables like energy density and pressure. """

        for particle in self.particles:
            particle.update()

    def init_interactions(self):
        """ ### 2. Initialize non-equilibrium interactions
            Non-equilibrium interactions of different particle species are treated by a\
            numerical integration of the Boltzmann equation for distribution functions in\
            the expanding space-time.

            Depending on the regime of the particle species involved and cosmological parameters, \
            each `Interaction` object populates `Particle.collision_integrals` array with \
            currently active `Integral` objects.
        """
        for interaction in self.interactions:
            interaction.initialize()

    def calculate_collisions(self):
        """ ### 3. Calculate collision integrals """

        particles = [particle for particle in self.particles if particle.collision_integrals]

        with utils.printoptions(precision=3, linewidth=100):
            for particle in particles:
                # with (utils.benchmark(lambda: "δf/f ({}) = {}".format(particle.symbol, particle.collision_integral / particle._distribution * self.params.h),
                #       self.log_throttler.output)):
                particle.collision_integral = particle.integrate_collisions()
                if not numpy.all(numpy.isfinite(particle.collision_integral)):
                    particle.collision_integral[numpy.isnan(particle.collision_integral)] = 0.

    def update_distributions(self):
        """ ### 4. Update particles distributions """

        if self.oscillations:
            self.oscillation_parameters()
            pattern, particles = self.oscillations

            if any(self.params.T < A.decoupling_temperature for A in particles):

                integrals = {A.flavour: A.collision_integral for A in particles}

                for A in particles:
                    A.collision_integral = sum(pattern[(A.flavour, B.flavour)] * integrals[B.flavour]
                                               for B in particles)

        for particle in self.particles:
            particle.update_distribution()

    def calculate_temperature_terms(self):
        """ ### 5. Calculate temperature equation terms """

        numerator = 0
        denominator = 0

        for particle in self.particles:
            numerator += particle.numerator()
            denominator += particle.denominator()

        return numerator, denominator

    def integrand(self, t, y):
        """ ## Temperature equation integrand

            Master equation for the temperature looks like

            \begin{equation}
                \frac{d (aT)}{dx} = \frac{\sum_i{N_i}}{\sum_i{D_i}}
            \end{equation}

            Where $N_i$ and $D_i$ represent contributions from different particle species.

            See definitions for different regimes:
              * [[Radiation|particles/RadiationParticle.py#master-equation-terms]]
              * [[Intermediate|particles/IntermediateParticle.py#master-equation-terms]]
              * [[Dust|particles/DustParticle.py#master-equation-terms]]
              * [[Non-equilibrium|particles/NonEqParticle.py#master-equation-terms]]
        """

        # 1\. Update particles states
        self.update_particles()
        # 2\. Initialize non-equilibrium interactions
        self.init_interactions()
        # 3\. Calculate collision integrals
        self.calculate_collisions()
        # 4\. Update particles distributions
        self.update_distributions()
        # 5\. Calculate temperature equation terms
        numerator, denominator = self.calculate_temperature_terms()

        if environment.get('LOGARITHMIC_TIMESTEP'):
            self.fraction = self.params.x * numerator / denominator
        else:
            self.fraction = numerator / denominator

        return self.fraction

    def save_params(self):
        self.data.append({
            'aT': self.params.aT,
            'T': self.params.T,
            'a': self.params.a,
            'x': self.params.x,
            'rho': self.params.rho,
            'N_eff': self.params.N_eff,
            't': self.params.t,
            'fraction': self.fraction,
            'S': self.params.S
        })

    def save(self):
        """ Save current Universe parameters into the data arrays or output files """
        self.save_params()

        if self.kawano and self.params.T <= self.kawano.T_kawano:

            #     t[s]         x    Tg[10^9K]   dTg/dt[10^9K/s] rho_tot[g cm^-3]     H[s^-1]
            # n nue->p e  p e->n nue  n->p e nue  p e nue->n  n e->p nue  p nue->n e

            rates = self.kawano.baryonic_rates(self.params.a)

            if environment.get('LOGARITHMIC_TIMESTEP'):
                dTdt = (self.fraction - self.params.aT) * self.params.H / self.params.a
            else:
                dTdt = (self.fraction - self.params.T / self.params.m) * self.params.H * self.params.m

            row = {
                self.kawano_data.columns[0]: self.params.t,
                self.kawano_data.columns[1]: self.params.x,
                self.kawano_data.columns[2]: self.params.T,
                self.kawano_data.columns[3]: dTdt,
                self.kawano_data.columns[4]: self.params.rho,
                self.kawano_data.columns[5]: self.params.H
            }

            row.update({self.kawano_data.columns[i]: rate
                        for i, rate in enumerate(rates, 6)})

            self.kawano_data.append(row)

            if self.log_throttler.output:
                print("KAWANO", self.kawano_data.row_repr(-1, names=True))
            self.kawano_log.write(self.kawano_data.row_repr(-1) + "\n")

    def init_log(self, folder=''):
        self.logfile = utils.ensure_path(os.path.join(self.folder, 'log.txt'))
        sys.stdout = utils.Logger(self.logfile)

    def log(self):
        """ Runtime log output """

        # Print parameters every now and then
        if self.log_throttler.output:
            print('[{clock}] #{step}\tt = {t:e} s\taT = {aT:e} MeV\tT = {T:e} MeV'
                  '\tδaT/aT = {daT:e}\tS = {S:e} MeV^3'
                  .format(clock=timedelta(seconds=int(time.time() - self.clock_start)),
                          step=self.step,
                          t=self.params.t / UNITS.s,
                          aT=self.params.aT / UNITS.MeV,
                          T=self.params.T / UNITS.MeV,
                          daT=self.fraction * self.params.h / self.params.aT,
                          S=self.params.S / UNITS.MeV**3))

    def total_entropy(self):
        return sum(particle.entropy for particle in self.particles) #* (self.params.a_ini/self.params.a)**3

    def total_energy_density(self):
        return sum(particle.energy_density for particle in self.particles)

    def QCD_transition(self, hadrons=None, quarkic_interactions=None, hadronic_interactions=None, secondary_interactions=None):
        self.data.truncate()

        dof_before_qcd = Params.entropic_dof_eq(self)

        self.particles = utils.particle_filter(
                            ['Up quark', 'Down quark', 'Charm quark',
                            'Strange quark', 'Top quark', 'Bottom quark', 'Gluon'],
                            self.particles
                            )
        self.add_particles(hadrons)

        dof_after_qcd = Params.entropic_dof_eq(self)

        self.params.aT = self.data['aT'][-1] * (dof_before_qcd/dof_after_qcd)**(1/3)
#        self.params.x = self.data['x'][-1] * (dof_before_qcd/dof_after_qcd)**(1/3)
#        self.params.a = self.data['a'][-1] * (dof_before_qcd/dof_after_qcd)**(1/3)
        self.params.update(self.total_energy_density(), self.total_entropy()) # T will be updated here...
#        self.params.T = CONST.lambda_QCD # That's why we change it here again
        self.update_particles()

        if quarkic_interactions and self.interactions:
            for inter in quarkic_interactions:
                self.interactions.remove(inter)
        if hadronic_interactions:
            self.interactions += (hadronic_interactions)
        if secondary_interactions:
            self.interactions += (secondary_interactions)
