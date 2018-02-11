# -*- coding: utf-8 -*-


import argparse
import os
from collections import defaultdict

from particles import Particle
from library.SM import particles as SMP, interactions as SMI
from library.NuMSM import particles as NuP, interactions as NuI
from evolution import Universe
from common import UNITS, Params, utils, LogSpacedGrid
from scipy.integrate import simps, cumtrapz
import numpy as np

parser = argparse.ArgumentParser(description='Run simulation for given mass and mixing angle')
parser.add_argument('--mass', default=33.9)
parser.add_argument('--theta', default=0.031)
parser.add_argument('--comment', default='')
args = parser.parse_args()

mass = float(args.mass) * UNITS.MeV
theta = float(args.theta)

folder = utils.ensure_dir(
    os.path.split(__file__)[0],
    "output",
    "mass={:e}_theta={:e}".format(mass / UNITS.MeV, theta)
    + args.comment
)

T_initial = 20. * UNITS.MeV
T_final = 0.0008 * UNITS.MeV
params = Params(T=T_initial,
                dy=0.003125)

universe = Universe(params=params, folder=folder)


from common import LinearSpacedGrid
linear_grid = LogSpacedGrid(MOMENTUM_SAMPLES=51, MAX_MOMENTUM=20*UNITS.MeV)
linear_grid_s = LinearSpacedGrid(MOMENTUM_SAMPLES=51, MAX_MOMENTUM=0.1*UNITS.MeV)


photon = Particle(**SMP.photon)

electron = Particle(**SMP.leptons.electron, **{'grid': linear_grid})

neutrino_e = Particle(**SMP.leptons.neutrino_e, **{'grid': linear_grid})
neutrino_mu = Particle(**SMP.leptons.neutrino_mu, **{'grid': linear_grid})
neutrino_tau = Particle(**SMP.leptons.neutrino_tau, **{'grid': linear_grid})

for neutrino in [neutrino_e, neutrino_mu, neutrino_tau]:
    neutrino.decoupling_temperature = 0 * UNITS.MeV


sterile = Particle(**NuP.dirac_sterile_neutrino(mass), **{'grid': linear_grid_s})
sterile.decoupling_temperature = T_initial

universe.add_particles([
    photon,

    electron,

    neutrino_e,
    neutrino_mu,
    neutrino_tau,

    sterile,
])

thetas = defaultdict(float, {
    'tau': theta,
})

interaction = NuI.sterile_leptons_interactions(
    thetas=thetas, sterile=sterile,
    neutrinos=[neutrino_e, neutrino_mu, neutrino_tau],
#    leptons=[]
    leptons=[electron]
)


def reaction_type(reaction):
    return sum(reactant.side for reactant in reaction)

#interaction.integrals = [integral for integral in interaction.integrals
#                            if reaction_type(integral.reaction) in [2]] 
for inter in interaction:
    inter.integrals = [integral for integral in inter.integrals
                             if reaction_type(integral.reaction) in [2]]


universe.interactions += (interaction)

"""
def step_monitor(universe):
    import numpy
    for particle in universe.particles:
        data = particle.data['params']
        if particle.mass > 0 and not particle.in_equilibrium:
            momenta = particle.grid.TEMPLATE 
            density = particle.density
            density_c = particle.density * particle.params.a**3 
            AA = simps(momenta**2 * particle.collision_integral, momenta)
            BB = simps(momenta**2 * particle._distribution / particle.energy(momenta), momenta) 
            integrand = (AA * particle.params.H)  / (particle.mass * BB)
            decay_rate = -integrand
            print(decay_rate)
            with open(os.path.join(folder, particle.name.replace(' ', '_') + ".decay_rate6.txt"), 'a') as f1:
                f1.write('{:e}'.format(particle.params.T / UNITS.MeV) + '\t' + '{:e}'.format(decay_rate / UNITS.MeV) + '\n')

"""


#"""
def step_monitor(universe):
    import numpy
    for particle in universe.particles:
        data = particle.data['params']
        if particle.mass > 0 and not particle.in_equilibrium:
            momenta = particle.grid.TEMPLATE 
            density = particle.density
            density_c = particle.density * particle.params.a**3 
            integrand = (particle.collision_integral * particle.params.H * particle.conformal_energy(particle.grid.TEMPLATE) / particle.mass / particle.params.a / particle._distribution)
            decay_rate = -integrand
            print(decay_rate / UNITS.MeV / 1.27031e-21) 
            with open(os.path.join(folder, particle.name.replace(' ', '_') + ".decay_rate6.txt"), 'a') as f1:
                f1.write('{:e}'.format(particle.params.T / UNITS.MeV) + '\t' + '{:e}'.format(particle.params.a) + '\t' + '\t'.join(['{:e}'.format(x) for x in decay_rate / UNITS.MeV]) + '\n')

#"""


###### NEGLECT THE REST ########

"""

def step_monitor(universe):

    for particle in universe.particles:
        data = particle.data['params']
        if particle.mass > 0 and not particle.in_equilibrium:
            momenta = particle.grid.TEMPLATE # MeV
            density = particle.density
            density_c = particle.density * particle.params.a**3 #/ UNITS.MeV**3 # MeV^3
            integrand = (particle.collision_integral * particle.params.H) # MeV
            collision = simps(momenta**2 * integrand, momenta)/(2*np.pi**2)
            decay_rate = -(collision / density / particle.params.a**3)
#            print("\n",collision / particle.params.a**3,"\n")
#            print("{}: Γ ={: .3e} MeV, τ ={: .3e} s, Y = n/S ={: .3e}".format(
#                particle.symbol,
#                decay_rate / UNITS.MeV,
#                1 / decay_rate / UNITS.s,
#                particle.density * universe.params.a**3 / universe.params.S
#            ))
         
            with open(os.path.join(folder, particle.name.replace(' ', '_') + ".decay_rate2.txt"), 'a') as f1:
                f1.write('{:e}'.format(particle.params.T / UNITS.MeV) + '\t' + '{:e}'.format(decay_rate / UNITS.MeV) + '\n')
#            with open(os.path.join(folder, particle.name.replace(' ', '_') + ".collision_integral.txt"), 'a') as f1:
#                if universe.step == 3:
#                    f1.write('0.000000e+00' + '\t' + '0.000000e+00' + '\t' + '\t'.join(['{:e}'.format(x) for x in particle.grid.TEMPLATE / UNITS.MeV]) + '\n')
#                f1.write('{:e}'.format(universe.params.T / UNITS.MeV) + '\t' + '\t'.join(['{:e}'.format(x) for x in particle._distribution]) + '\n')
#"""


"""
def step_monitor(universe):

    for particle in universe.particles:
        data = particle.data['params']
        if particle.mass > 0 and not particle.in_equilibrium and len(particle.data['params']) > 3:
            momenta = particle.grid.TEMPLATE # MeV
            density = particle.density
            density_c = particle.density * universe.params.a**3 #/ UNITS.MeV**3 # MeV^3
            integrand = (particle.collision_integral * universe.params.H) # MeV
            collision = simps(momenta**2 * integrand, momenta)/(2*np.pi**2)
            decay_rate = -(collision / (density_c * UNITS.MeV**3))
#                print("\n\n",collision,"\n\n")
#            print("{}: Γ ={: .3e} MeV, τ ={: .3e} s, Y = n/S ={: .3e}".format(
#                particle.symbol,
#                decay_rate / UNITS.MeV,
#                1 / decay_rate / UNITS.s,
#                particle.density * universe.params.a**3 / universe.params.S
#            ))

#            with open(os.path.join(folder, particle.name.replace(' ', '_') + ".decay_rate.txt"), 'a') as f1:
#                f1.write('{:e}'.format(decay_rate / UNITS.MeV) + '\n')
            with open(os.path.join(folder, particle.name.replace(' ', '_') + ".collision_integral2.txt"), 'a') as f1:
                if universe.step == 3:
                    f1.write('0.000000e+00' + '\t' + '0.000000e+00' + '\t' + '\t'.join(['{:e}'.format(x) for x in particle.grid.TEMPLATE / UNITS.MeV]) + '\n')
                f1.write('{:e}'.format(universe.params.T / UNITS.MeV) + '\t' + '{:e}'.format(universe.params.t / UNITS.MeV) + '\t' + '\t'.join(['{:e}'.format(x) for x in particle._distribution]) + '\n')
"""



"""
def step_monitor(universe):
    import numpy

    for particle in universe.particles:
        data = particle.data['params']
        if particle.mass > 0 and not particle.in_equilibrium and len(particle.data['params']) > 3:
            decay_rate = -(
                (numpy.log(data['density'][-1]) - numpy.log(data['density'][-2]))
                / ((data['t'][-1] - data['t'][-2]))
                + 3 * universe.params.H
            )

            print("{}: Γ ={: .3e} MeV, τ ={: .3e} s, Y = n/S ={: .3e}".format(
                particle.symbol,
                decay_rate / UNITS.MeV,
                1 / decay_rate / UNITS.s,
                particle.density * universe.params.a**3 / universe.params.S
            ))

            with open(os.path.join(folder, particle.name.replace(' ', '_') + ".decay_rate6.txt"), 'a') as f1:
                f1.write('{:e}'.format(particle.params.T / UNITS.MeV) + '\t' + '{:e}'.format(decay_rate / UNITS.MeV) + '\n')
"""
universe.step_monitor = step_monitor

universe.evolve(T_final)

