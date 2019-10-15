# -*- coding: utf-8 -*-

import os
import numpy as np
import environment
from common import CONST, UNITS, utils
from collections import Counter
from scipy.integrate import simps
from interactions.four_particle.cpp.integral import CollisionIntegralKind


def cm_momentum(mass_1, mass_2, mass_3):
    return np.sqrt(
            (mass_1**2 - mass_2**2 - mass_3**2)**2
            - (2 * mass_2 * mass_3)**2
            ) / (2 * mass_1)

def return_function(interaction, ps):
    if interaction.kind in [CollisionIntegralKind.Full, CollisionIntegralKind.Full_vacuum_decay] or hasattr(interaction.particle, 'fast_decay'):
        return np.zeros(len(ps)), np.zeros(len(ps))
    return np.zeros(len(ps))

def three_particle_bounds_creation(reaction=None):
    " Returns kinematically allowed lower and upper momentum bounds for particle A in reaction A + B <-- C "
    specie_1 = reaction[0].specie
    specie_2 = reaction[1].specie
    specie_3 = reaction[2].specie

    beta = specie_3.grid.MAX_MOMENTUM / np.sqrt(specie_3.conformal_mass**2 + specie_3.grid.MAX_MOMENTUM**2)
    gamma = 1 / np.sqrt(1 - beta**2)

    momentum_1_cm = cm_momentum(specie_3.conformal_mass, specie_2.conformal_mass, specie_1.conformal_mass)

    energy_1_cm = np.sqrt(specie_1.conformal_mass**2 + momentum_1_cm**2)
    energy_1_lab_max = gamma * (energy_1_cm + beta * momentum_1_cm)

    if specie_1.mass == 0.:
        energy_1_lab_min = gamma * (energy_1_cm - beta * momentum_1_cm)
    else:
        energy_1_lab_min = specie_1.conformal_mass

    momentum_1_lab_max = np.sqrt(energy_1_lab_max**2 - specie_1.conformal_mass**2)
    momentum_1_lab_min = np.sqrt(energy_1_lab_min**2 - specie_1.conformal_mass**2)
    return momentum_1_lab_min, momentum_1_lab_max

def three_particle_grid_bounds_creation(reaction=None):
    """ Returns grid elements (slice_2) for which three particle collision integral will be computed.
        Collision integral is evaluated to zero for other slices. """
    lower_bound, upper_bound = three_particle_bounds_creation(reaction)
    grid = reaction[0].specie.grid.TEMPLATE
    lower_element_grid = max(0, np.searchsorted(grid, lower_bound, side='right') - 1)
    upper_element_grid = min(len(grid) - 1, np.searchsorted(grid, upper_bound))

    slice_1 = [0] * len(grid[:lower_element_grid])
    slice_2 = grid[lower_element_grid:upper_element_grid + 1]
    slice_3 = [0] * (len(grid[upper_element_grid:]) - 1)

    if lower_element_grid == upper_element_grid:
        slice_2 = []

    return slice_1, np.array(slice_2), slice_3

def four_particle_bounds_creation(reaction=None):
    " Returns kinematically allowed lower and upper momentum bounds for particle A in reaction A + B + C <-- D "
    specie_1 = reaction[0].specie
    specie_2 = reaction[1].specie
    specie_3 = reaction[2].specie
    specie_4 = reaction[3].specie

    max_momentum = np.sqrt(
                    (   np.sqrt(specie_4.grid.MAX_MOMENTUM**2 + specie_4.conformal_mass**2)
                        - specie_2.conformal_mass - specie_3.conformal_mass
                    )**2
                    - specie_1.conformal_mass**2)

    return max(max_momentum, environment.get('MAX_MOMENTUM_MEV') * UNITS.MeV)

def four_particle_grid_cutoff_creation(reaction=None):
    """ Returns grid elements (slice_1) for which four particle collision integral will be computed.
        Collision integral is evaluated to zero for other slice. """
    upper_bound = four_particle_bounds_creation(reaction)
    grid = reaction[0].specie.grid.TEMPLATE
    upper_element_grid = min(len(grid) - 1, np.searchsorted(grid, upper_bound))

    slice_1 = grid[:upper_element_grid + 1]
    slice_2 = [0] * (len(grid[upper_element_grid:]) - 1)

    return np.array(slice_1), slice_2

def four_particle_grid_cutoff_scattering(particle):
    grid = particle.grid.TEMPLATE

    if particle.name == 'Sterile neutrino (Dirac)':
        max_mom = particle.grid.MAX_MOMENTUM #3 * particle.params.T
    elif particle.name != 'Sterile neutrino (Dirac)' and environment.get('HNL_ENERGY'):
        HNL_energy = float(environment.get('HNL_ENERGY'))
        max_mom = max(np.sqrt(HNL_energy**2 - particle.conformal_mass**2), environment.get('MAX_MOMENTUM_MEV') * UNITS.MeV)
    else:
        max_mom = environment.get('MAX_MOMENTUM_MEV') * UNITS.MeV

    upper_element_grid = min(len(grid) - 1, np.searchsorted(grid, max_mom))

    slice_1 = grid[:upper_element_grid + 1]
    slice_2 = [0] * (len(grid[upper_element_grid:]) - 1)

    return np.array(slice_1), slice_2

def four_particle_grid_cutoff_decay(particle):
    inv_f_min = 1e100
    aT = 1 * UNITS.MeV

    if particle.name == 'Sterile neutrino (Dirac)':
        max_mom = particle.grid.MAX_MOMENTUM
    elif particle.name != 'Sterile neutrino (Dirac)' and environment.get('HNL_ENERGY'):
        HNL_energy = float(environment.get('HNL_ENERGY'))
        max_mom = max(np.sqrt(HNL_energy**2 - particle.conformal_mass**2), environment.get('MAX_MOMENTUM_MEV') * UNITS.MeV)
    else:
        max_mom = particle.grid.MAX_MOMENTUM

    grid = particle.grid.TEMPLATE
    upper_element_grid = min(len(grid) - 1, np.searchsorted(grid, max_mom))

    slice_1 = grid[:upper_element_grid + 1]
    slice_2 = [0] * (len(grid[upper_element_grid:]) - 1)

    return np.array(slice_1), slice_2

def grid_cutoff_4p(interaction):
    # Cut off grid for creation reactions
    if utils.reaction_type(interaction).CREATION:
        slice_1, slice_2 = four_particle_grid_cutoff_creation(interaction.reaction)
        ps = slice_1 / interaction.particle.params.aT

    # Cut off grid for scattering reactions
    if utils.reaction_type(interaction).SCATTERING:
        slice_1, slice_2 = four_particle_grid_cutoff_scattering(interaction.particle)
        ps = slice_1 / interaction.particle.params.aT

    # Cut off grid for decay reactions
    if utils.reaction_type(interaction).DECAY:
        slice_1, slice_2 = four_particle_grid_cutoff_decay(interaction.particle)
        ps = slice_1 / interaction.particle.params.aT

    return ps, slice_1, slice_2

def grid_cutoff_3p(interaction, ps):
    slice_1 = []
    slice_3 = []
    if utils.reaction_type(interaction).CREATION:
        slice_1, slice_2, slice_3 = three_particle_grid_bounds_creation(interaction.reaction)
        if not slice_2.any():
            return False
        ps = slice_2 / interaction.particle.params.aT

    return ps, slice_1, slice_3

def decoupling_temperature_relativistic(mass, mixing_angle):
    muon_mass = 105.658 * UNITS.MeV

    if mass > CONST.lambda_QCD:
        g_star = 61.75
    if muon_mass < mass < CONST.lambda_QCD:
        g_star = 14.25
    else:
        g_star = 10.75

    theta_sq = mixing_angle**2

    return (1.66 * np.sqrt(g_star) / (CONST.M_p * CONST.G_F**2 * theta_sq))**(1/3)

def decoupling_temperature(mass, mixing_angle):
    T_dec_rel = decoupling_temperature_relativistic(mass, mixing_angle)
    if 1.5 * T_dec_rel > 1.5 * mass:
        os.environ['Relativistic_decoupling'] = 'True'
        return 1.5 * T_dec_rel
    elif mass <= 1.5 * T_dec_rel <= 1.5 * mass:
        return 1.5 * mass
    return mass

def grid_resolution_HNL(mass, mixing_angle):
    aT = 1 * UNITS.MeV
    delta_y_HNL = np.sqrt(3 - 8 / np.pi) * aT

    return delta_y_HNL

def grid_resolution_neutrino(mass, mixing_angle):
    muon_mass = 105.658 * UNITS.MeV
    pion_mass = 134.98 * UNITS.MeV

    if mass < pion_mass:
        if mass > muon_mass and mass - muon_mass < 10 * UNITS.MeV:
            return environment.get('FOUR_PARTICLE_GRID_RESOLUTION') * UNITS.MeV / 2.
        return environment.get('FOUR_PARTICLE_GRID_RESOLUTION') * UNITS.MeV

    meson_masses = np.array([134.98, 139.57, 547.86, 775.11, 782.65, 957.78, 1019.46]) * UNITS.MeV
    meson_mass = meson_masses[np.searchsorted(meson_masses, mass) - 1]
    delta_y_nu = 2 * cm_momentum(mass, 0, meson_mass) * grid_resolution_HNL(mass, mixing_angle) / mass

    return min(delta_y_nu * .75, environment.get('FOUR_PARTICLE_GRID_RESOLUTION') * UNITS.MeV / 3.)

def grid_resolution_meson(mass, mixing_angle):
    if mass < 134.98 * UNITS.MeV:
        return environment.get('FOUR_PARTICLE_GRID_RESOLUTION') * UNITS.MeV

    meson_masses = np.array([134.98, 139.57, 493.677, 497.611, 547.86, 775.11, 782.65, 957.78, 1019.46]) * UNITS.MeV
    meson_mass = meson_masses[np.searchsorted(meson_masses, mass) - 1]
    momentum_cm = cm_momentum(mass, 0, meson_mass)
    energy_cm = np.sqrt(momentum_cm**2 + meson_mass**2)

    delta_y_meson = np.sqrt((energy_cm + grid_resolution_HNL(mass, mixing_angle) * momentum_cm / mass)**2 - meson_mass**2) \
                    - np.sqrt((energy_cm - grid_resolution_HNL(mass, mixing_angle) * momentum_cm / mass)**2 - meson_mass**2)

    if not delta_y_meson > 0:
        delta_y_meson = environment.get('FOUR_PARTICLE_GRID_RESOLUTION') * UNITS.MeV / 3.

    return min(delta_y_meson / 3., environment.get('FOUR_PARTICLE_GRID_RESOLUTION') * UNITS.MeV)

def grid_resolution_lepton(mass, mixing_angle, std=False):
    muon_mass = 105.658 * UNITS.MeV
    pion_mass = 139.57 * UNITS.MeV

    if mass < pion_mass or std:
        if mass > muon_mass and mass - muon_mass < 10 * UNITS.MeV:
            return environment.get('FOUR_PARTICLE_GRID_RESOLUTION') * UNITS.MeV / 2.
        return environment.get('FOUR_PARTICLE_GRID_RESOLUTION') * UNITS.MeV

    momentum_cm = cm_momentum(pion_mass, 0, muon_mass)
    energy_cm = np.sqrt(momentum_cm**2 + muon_mass**2)

    meson_resolution = 0.67 * np.sqrt(pion_mass * UNITS.MeV * 0.5)
    delta_y_lepton = np.sqrt((energy_cm + meson_resolution * momentum_cm / pion_mass)**2 - muon_mass**2) \
                    - np.sqrt((energy_cm - meson_resolution * momentum_cm / pion_mass)**2 - muon_mass**2)

    return min(delta_y_lepton, environment.get('FOUR_PARTICLE_GRID_RESOLUTION') * UNITS.MeV)

def grid_params_HNL(mass, mixing_angle):
    aT = 1 * UNITS.MeV
    inv_f_min = 1e100
    a_ini = aT / decoupling_temperature(mass, mixing_angle)
    MAX_MOMENTUM = np.sqrt((aT * np.log(inv_f_min))**2 - (mass * a_ini)**2)

    MOMENTUM_SAMPLES = MAX_MOMENTUM / grid_resolution_HNL(mass, mixing_angle)

    return int(MOMENTUM_SAMPLES), MAX_MOMENTUM

def grid_params_neutrino(mass, mixing_angle):
    MAX_MOMENTUM = mass * environment.get('MAX_SCALE_FACTOR')
    MOMENTUM_SAMPLES = MAX_MOMENTUM / grid_resolution_neutrino(mass, mixing_angle)

    return int(MOMENTUM_SAMPLES), MAX_MOMENTUM

def grid_params_meson(mass, mixing_angle):
    pion_mass = 134.98 * UNITS.MeV

    if mass < pion_mass:
        MAX_MOMENTUM = environment.get('MAX_MOMENTUM_MEV') * UNITS.MeV
    elif mass - pion_mass < 0.35 * UNITS.MeV:
        MAX_MOMENTUM = 10 * environment.get('MAX_SCALE_FACTOR') * UNITS.MeV
    else:
        MAX_MOMENTUM = np.sqrt(mass**2 - pion_mass**2) * environment.get('MAX_SCALE_FACTOR')

    MOMENTUM_SAMPLES = MAX_MOMENTUM / grid_resolution_meson(mass, mixing_angle)

    return int(MOMENTUM_SAMPLES), MAX_MOMENTUM

def grid_params_lepton(mass, mixing_angle, std=False):
    muon_mass = 105.658 * UNITS.MeV

    if mass < muon_mass or std:
        MAX_MOMENTUM = mass * environment.get('MAX_SCALE_FACTOR')
    elif mass - muon_mass < 0.35 * UNITS.MeV:
        MAX_MOMENTUM = 10 * environment.get('MAX_SCALE_FACTOR') * UNITS.MeV
    else:
        MAX_MOMENTUM = np.sqrt(mass**2 - muon_mass**2) * environment.get('MAX_SCALE_FACTOR')

    MOMENTUM_SAMPLES = MAX_MOMENTUM / grid_resolution_lepton(mass, mixing_angle, std=std)

    return int(MOMENTUM_SAMPLES), MAX_MOMENTUM

def Neglect3pInteraction(interaction, ps):
    if interaction.reaction[0].specie.mass == 0 and interaction.reaction[1].side == 1:
        return True

    if utils.reaction_type(interaction).CREATION and interaction.reaction[-1].specie.decayed:
        if interaction.kind in [CollisionIntegralKind.Full, CollisionIntegralKind.Full_vacuum_decay] or hasattr(interaction.particle, 'fast_decay'):
            return True

    if utils.reaction_type(interaction).CREATION and hasattr(interaction.reaction[-1].specie, 'fast_decay') and interaction.reaction[-1].specie.num_creation == 0\
    or utils.reaction_type(interaction).DECAY and hasattr(interaction.reaction[0].specie, 'fast_decay') and interaction.reaction[0].specie.num_creation == 0:
        if interaction.kind in [CollisionIntegralKind.Full, CollisionIntegralKind.Full_vacuum_decay] or hasattr(interaction.particle, 'fast_decay'):
            return True

    return False

def Neglect4pInteraction(interaction, ps):
    # If particle has decayed, don't calculate creation integrals for decay products
    if utils.reaction_type(interaction).CREATION and interaction.reaction[-1].specie.decayed:
        return True

    # If particles have diluted to such extent that they can be neglected
    scat_thr = 1e-10 * (interaction.particle.params.a_ini / 10)**3
    if utils.reaction_type(interaction).SCATTERING and (interaction.reaction[0].specie.density / interaction.reaction[0].specie.data['params']['density'][0] < scat_thr\
    or interaction.reaction[1].specie.density / interaction.reaction[1].specie.data['params']['density'][0] < scat_thr)\
    and (interaction.reaction[2].specie.density / interaction.reaction[2].specie.data['params']['density'][0] < scat_thr or\
    interaction.reaction[3].specie.density / interaction.reaction[3].specie.data['params']['density'][0] < scat_thr):
        return True
    #TODO: Improve this (zero initial density etc)

    # If there are no muons/mesons created yet, skip creation and decay reactions
    if utils.reaction_type(interaction).CREATION and hasattr(interaction.reaction[-1].specie, 'fast_decay') and interaction.reaction[-1].specie.num_creation == 0\
    or utils.reaction_type(interaction).DECAY and hasattr(interaction.reaction[0].specie, 'fast_decay') and interaction.reaction[0].specie.num_creation == 0:
        return True

    # Decoupling of scattering reactions involving HNL
    if utils.reaction_type(interaction).SCATTERING and any(item.specie.name == 'Sterile neutrino (Dirac)' for item in interaction.reaction)\
    and (environment.get('Relativistic_decoupling') and interaction.particle.params.T < interaction.particle.params.m / interaction.particle.params.a_ini / 15. or interaction.particle.params.T < 1. * UNITS.MeV):
        return True

    # If temperature is higher than HNL mass, skip decay reaction to prevent incorrect computation of collision integral
    if utils.reaction_type(interaction).DECAY and interaction.particle.name == 'Sterile neutrino (Dirac)' and interaction.particle.params.T > interaction.particle.mass:
       return True

    return False

def CollisionMultiplier4p(interaction):
    if 'Sterile neutrino (Dirac)' in [item.specie.name for item in interaction.reaction] or \
    (utils.reaction_type(interaction).CREATION and not \
    (interaction.reaction[3].specie.majorana and interaction.particle.Q)):
        left = Counter(item.specie for item in interaction.reaction if item.side == -1)
        right = Counter(item.specie for item in interaction.reaction if item.side == 1)
        if left[interaction.reaction[0].specie] == 2 and right[interaction.reaction[0].specie] in [0, 1]:
            return 2.
        if left[interaction.reaction[0].specie] == 3 and right[interaction.reaction[0].specie] == 0:
            return 3.
    return 1.

def CollisionMultiplier3p(interaction):
    if 'Sterile neutrino (Dirac)' in [item.specie.name for item in interaction.reaction] or \
    utils.reaction_type(interaction).CREATION and not \
    (interaction.reaction[2].specie.majorana and interaction.particle.Q):
        left = Counter(item.specie for item in interaction.reaction if item.side == -1)
        right = Counter(item.specie for item in interaction.reaction if item.side == 1)
        if left[interaction.reaction[0].specie] == 2 and right[interaction.reaction[0].specie] == 0:
            return 2.
    return 1.

def store_energy(interaction):
    if interaction.particle.name == 'Sterile neutrino (Dirac)':
        os.environ['HNL_ENERGY'] = str(np.sqrt((interaction.particle.grid.MAX_MOMENTUM/10)**2 + interaction.particle.conformal_mass**2))

def interpolation_4p(interaction, ps, slice_1):
    interpolate = False
    grid = interaction.particle.grid
    if grid.MAX_MOMENTUM / (grid.MOMENTUM_SAMPLES - 1) < environment.get('FOUR_PARTICLE_GRID_RESOLUTION') * UNITS.MeV:
        steps = np.ceil(slice_1[-1] / (environment.get('FOUR_PARTICLE_GRID_RESOLUTION') * UNITS.MeV))
        interp_pos = np.searchsorted(ps, np.linspace(ps[0], ps[-1], steps))
        ps = ps[interp_pos]
        interpolate = True
    return ps, interpolate

def scaling(interaction, fullstack, constant):
    grid = interaction.particle.grid
    if utils.reaction_type(interaction).CREATION and hasattr(interaction.reaction[-1].specie, 'fast_decay'):
        sym = ''.join([item.specie.symbol for item in interaction.reaction[:-1]])
        for key in interaction.reaction[-1].specie.BR:
            if Counter(sym) == Counter(key):
                BR = interaction.reaction[-1].specie.BR[key]
        dof = interaction.particle.dof if interaction.particle.majorana else interaction.particle.dof / 2.
        created = simps(fullstack * constant * dof * grid.TEMPLATE**2, grid.TEMPLATE)
        if created == 0.:
            return False
        scaling = BR * interaction.reaction[-1].specie.num_creation / created
        fullstack *= scaling
    return fullstack

def has_decayed(particle, ps):
    dec_thr = 1e-10 * (particle.params.a_ini / 10)**3
    if len(particle.data['distribution']) > 1 and particle.mass != 0. and particle.data['params']['density'][0] != 0. \
    and particle.density / particle.data['params']['density'][0] < dec_thr and not hasattr(particle, 'fast_decay'):
        particle.decayed = True
        particle._distribution = np.zeros(len(ps))
        if particle.name == 'Sterile neutrino (Dirac)':
            os.environ['STERILE_DECAYED'] = 'True'
        return True
    return False

def Icoll_fast_decay(particle, ps):
    F1s = []
    Ffs = []
    Ffs_temp = []
    dofs = []
    for integral in particle.collision_integrals:
        F1, Ff = integral.integrate(ps, stepsize=particle.params.h)
        F1s.append(F1)
        Ffs_temp.append(Ff)
        dofs.append(particle.dof if particle.majorana else particle.dof / 2)

    distr_backg = particle.equilibrium_distribution(conf_mass=particle.mass * particle.aT / particle.decoupling_temperature)\
                * np.exp(-(particle.params.t - particle.t_decoupling) / particle.lifetime) if particle.thermal_dyn else np.zeros(len(ps))

    if hasattr(particle, 'thermalization'):
        distr_ini = distr_backg + sum(F1s) * particle.params.h
        distr_bef = distr_backg + simps(sum(F1s) * particle.grid.TEMPLATE**2, particle.grid.TEMPLATE) * particle.params.h \
                    / (2 * np.pi * particle.conformal_mass * particle.aT)**(3/2) \
                    * np.exp(-particle.grid.TEMPLATE**2 / (2 * particle.conformal_mass * particle.aT))

    else:
        distr_bef = distr_backg + sum(F1s) * particle.params.h

    particle.num_creation = simps(sum(F1s*np.array(dofs)[:,None]) * particle.grid.TEMPLATE**2, particle.grid.TEMPLATE)

    for index, Ff in enumerate(Ffs_temp):
        integral = particle.collision_integrals[index]
        if utils.reaction_type(integral.reaction).DECAY: # line not necessary
            sym = ''.join([item.specie.symbol for item in integral.reaction[1:]])
            for key in integral.reaction[0].specie.BR:
                if Counter(sym) == Counter(key):
                    BR = integral.reaction[0].specie.BR[key]
            dof = particle.dof if particle.majorana else particle.dof / 2
            decayed = simps(-1 * distr_bef * Ff * dof * particle.grid.TEMPLATE**2, particle.grid.TEMPLATE)
            if decayed == 0.:
                Ffs.append(np.zeros(len(Ff)))
            else:
                scaling = BR * particle.num_creation / decayed
                Ffs.append(scaling * Ff * distr_bef)
        else:
            Ffs.append(Ff * distr_bef)

    if hasattr(particle, 'thermalization'):
        coll_creation_th = sum(F1s)
        coll_thermalization_th = (distr_bef - distr_ini)/particle.params.h
        coll_decay_th = sum(Ffs)
        coll_integral_th = coll_creation_th + coll_thermalization_th + coll_decay_th
        particle._distribution = distr_bef
        return coll_integral_th

    coll_creation = sum(F1s)
    coll_decay = sum(Ffs)
    coll_integral = coll_creation + coll_decay
    particle._distribution = distr_bef
    return coll_integral
