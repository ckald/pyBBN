import argparse
import os.path as op
from collections import defaultdict
import numpy as np

import environment
from particles import Particle
from library.SM import particles as SMP, interactions as SMI
from library.NuMSM import particles as NuP, interactions as NuI
from evolution import Universe
from common import CONST, UNITS, Params, utils, kinematics, LinearSpacedGrid
from interactions.four_particle.cpp.integral import CollisionIntegralKind

parser = argparse.ArgumentParser(description='Run simulation for given mass and mixing and angle')
parser.add_argument('--mass', required=True)
parser.add_argument('--mixing', required=True, choices=['electron', 'muon', 'tau'])
parser.add_argument('--theta', required=True)
args = parser.parse_args()

mass = float(args.mass) * UNITS.MeV
mixing = args.mixing
theta = float(args.theta)

folder = utils.ensure_dir(
    op.split(__file__)[0],
    "output_{}_mixing".format(mixing),
    "mass={:e}_theta={:e}_flav={}".format(mass / UNITS.MeV, theta, mixing)
)

T_initial = kinematics.decoupling_temperature(mass, theta)
T_freeze = 0.1 * UNITS.MeV
T_final = 0.0008 * UNITS.MeV
params = Params(T=T_initial,
                dy=0.003125)

universe = Universe(params=params, folder=folder)

samples_HNL, max_mom_HNL = kinematics.grid_params_HNL(mass, theta)
samples_neutrino, max_mom_neutrino = kinematics.grid_params_neutrino(mass, theta)
samples_muon, max_mom_muon = kinematics.grid_params_lepton(mass, theta)
samples_electron, max_mom_electron = kinematics.grid_params_lepton(mass, theta, std=True)
samples_meson, max_mom_meson = kinematics.grid_params_meson(mass, theta)

grid_HNL = LinearSpacedGrid(MOMENTUM_SAMPLES=samples_HNL, MAX_MOMENTUM=max_mom_HNL)
grid_neutrino = LinearSpacedGrid(MOMENTUM_SAMPLES=samples_neutrino, MAX_MOMENTUM=max_mom_neutrino)
grid_muon = LinearSpacedGrid(MOMENTUM_SAMPLES=samples_muon, MAX_MOMENTUM=max_mom_muon)
grid_electron = LinearSpacedGrid(MOMENTUM_SAMPLES=samples_electron, MAX_MOMENTUM=max_mom_electron)
grid_meson = LinearSpacedGrid(MOMENTUM_SAMPLES=samples_meson, MAX_MOMENTUM=max_mom_meson)

photon = Particle(**SMP.photon)

electron = Particle(**SMP.leptons.electron, grid=grid_electron)
muon = Particle(**SMP.leptons.muon, grid=grid_muon)
tau = Particle(**SMP.leptons.tau, grid=grid_muon)

neutrino_decoupling_temperature = 5. * UNITS.MeV
neutrino_e = Particle(**SMP.leptons.neutrino_e, grid=grid_neutrino)
neutrino_mu = Particle(**SMP.leptons.neutrino_mu, grid=grid_neutrino)
neutrino_tau = Particle(**SMP.leptons.neutrino_tau, grid=grid_neutrino)
neutrino_e.decoupling_temperature = neutrino_decoupling_temperature
neutrino_mu.decoupling_temperature = neutrino_decoupling_temperature
neutrino_tau.decoupling_temperature = neutrino_decoupling_temperature

sterile = Particle(**NuP.dirac_sterile_neutrino(mass, theta), grid=grid_HNL)
sterile.decoupling_temperature = T_initial

up = Particle(**SMP.quarks.up)
down = Particle(**SMP.quarks.down)
strange = Particle(**SMP.quarks.strange)
charm = Particle(**SMP.quarks.charm)
bottom = Particle(**SMP.quarks.bottom)
gluon = Particle(**SMP.gluon)

charged_pion = Particle(**SMP.hadrons.charged_pion, grid=grid_meson)
neutral_pion = Particle(**SMP.hadrons.neutral_pion, grid=grid_meson)

heavy_meson_decoupling_temperature = 7. * UNITS.MeV
eta = Particle(**SMP.hadrons.eta, grid=grid_meson, thermal_dyn=False, decoupling_temperature=heavy_meson_decoupling_temperature)
neutral_rho = Particle(**SMP.hadrons.neutral_rho, grid=grid_meson, thermal_dyn=False, decoupling_temperature=heavy_meson_decoupling_temperature)
charged_rho = Particle(**SMP.hadrons.charged_rho, grid=grid_meson, thermal_dyn=False, decoupling_temperature=heavy_meson_decoupling_temperature)
eta_prime = Particle(**SMP.hadrons.eta_prime, grid=grid_meson, thermal_dyn=False, decoupling_temperature=heavy_meson_decoupling_temperature)
omega = Particle(**SMP.hadrons.omega, grid=grid_meson, thermal_dyn=False, decoupling_temperature=heavy_meson_decoupling_temperature)

universe.add_particles([
    photon,

    electron,

    neutrino_e,
    neutrino_mu,
    neutrino_tau,

    sterile,
])

if T_initial >= 10 * UNITS.MeV or mass > muon.mass:
    universe.add_particles([muon])

if T_initial >= 150 * UNITS.MeV:
    universe.add_particles([tau])

if T_initial > CONST.lambda_QCD:
    universe.add_particles([
        up, down,
        strange, charm,
        bottom,
        gluon
    ])

decaying_mesons = []
if (T_initial >= 10 * UNITS.MeV or mass > neutral_pion.mass):
    decaying_mesons += [neutral_pion, charged_pion]

for meson in [eta, charged_rho, neutral_rho, eta_prime, omega]:
    if mass > meson.mass:
        decaying_mesons.append(meson)

if T_initial <= CONST.lambda_QCD:
    universe.add_particles(decaying_mesons)

thetas = defaultdict(float, {
    mixing: theta
})

SM_interactions = SMI.neutrino_interactions(
    leptons=[electron],
    neutrinos=[
        neutrino_e,
        neutrino_mu,
        neutrino_tau
    ],
)

SM_interactions = utils.interaction_filter(['Muon'], SM_interactions)

if T_initial >= 150 * UNITS.MeV:
    populated_charged_leptons = [electron, muon, tau]
else:
    populated_charged_leptons = [electron, muon]

HNL_leptons_interactions = NuI.sterile_leptons_interactions(
    thetas=thetas, sterile=sterile,
    neutrinos=[neutrino_e, neutrino_mu, neutrino_tau],
    leptons=populated_charged_leptons
)

sterile_hadrons_interactions = NuI.sterile_hadrons_interactions(
    thetas=thetas, sterile=sterile,
    neutrinos=[neutrino_e, neutrino_mu, neutrino_tau],
    leptons=populated_charged_leptons,
    mesons=decaying_mesons
)

secondary_interactions_leptons = NuI.interactions_decay_products(
    interactions_primary=[HNL_leptons_interactions],
    neutrinos=[neutrino_e, neutrino_mu, neutrino_tau],
    leptons=[electron, muon],
    mesons=[],
    photon=[],
    kind=CollisionIntegralKind.F_decay
)

secondary_interactions_mesons = NuI.interactions_decay_products(
    interactions_primary=[sterile_hadrons_interactions],
    muon_dec=True,
    neutrinos=[neutrino_e, neutrino_mu, neutrino_tau],
    leptons=[electron, muon],
    mesons=decaying_mesons,
    photon=[photon],
    kind=CollisionIntegralKind.F_decay
)


sterile_quark_interactions = NuI.sterile_quark_interactions(
    thetas=thetas, sterile=sterile,
    neutrinos=[neutrino_e, neutrino_mu, neutrino_tau],
    leptons=[electron, muon, tau],
    quarks=[up, down, strange, charm, bottom]
)

sterile_quark_interactions = utils.interaction_filter(
    ['Up quark', 'Down quark', 'Charm quark',
     'Strange quark', 'Top quark', 'Bottom quark', 'Gluon'],
    sterile_quark_interactions
)

universe.interactions += (SM_interactions + HNL_leptons_interactions)

if mass >= muon.mass or T_initial >= 10 * UNITS.MeV:
    universe.interactions += secondary_interactions_leptons

if T_initial > CONST.lambda_QCD:
    universe.interactions += sterile_quark_interactions
else:
    if (T_initial >= 10 * UNITS.MeV or mass > neutral_pion.mass):
        universe.interactions += sterile_hadrons_interactions

    if mass > neutral_pion.mass:
        universe.interactions += secondary_interactions_mesons

universe.init_oscillations(pattern_function=SMP.leptons.oscillations_map,
                           particles=(neutrino_e, neutrino_mu, neutrino_tau),
                           matter_effects=True)

universe.init_kawano(electron=electron, neutrino=neutrino_e)


for particle in [neutrino_e, neutrino_mu, neutrino_tau, sterile]:
    with open(op.join(folder, particle.name.replace(' ', '_') + ".grid.npy"), 'wb') as f:
        np.save(f, particle.grid.TEMPLATE / UNITS.MeV, allow_pickle=False)


def step_monitor(universe):
    if universe.step % 10 == 0:
        for particle in [neutrino_e, neutrino_mu, neutrino_tau, sterile]:
            with open(op.join(folder, particle.name.replace(' ', '_') + ".rho.txt"), 'a') as f:
                f.write('{:e}\t{:e}\t{:e}\t{:e}\t{:e}\n'.format(
                    universe.params.a, universe.params.T / UNITS.MeV,
                    universe.params.aT / UNITS.MeV,
                    particle.energy_density / UNITS.MeV**4,
                    particle.density / UNITS.MeV**3
                ))

    if universe.step % 100 == 0:
        for particle in [neutrino_e, neutrino_mu, neutrino_tau, sterile]:
            with open(op.join(folder, particle.name.replace(' ', '_') + ".f.npy"), 'wb') as f:
                np.save(f, particle.data['distribution'].data, allow_pickle=False)


universe.step_monitor = step_monitor

print('Mass:', mass / UNITS.MeV, 'theta^2:', theta**2)
print('T_initial:', T_initial / UNITS.MeV)
print('T_freeze:', T_freeze / UNITS.MeV)
print('T_final:', T_final / UNITS.MeV)
print('HNL\n', 'Samples: ', samples_HNL, 'MAX_MOM:', max_mom_HNL / UNITS.MeV, 'T_dec:', sterile.decoupling_temperature / UNITS.MeV)
print('neutrino\n', 'Samples: ', samples_neutrino, 'MAX_MOM:', max_mom_neutrino / UNITS.MeV, 'T_dec:', neutrino_e.decoupling_temperature / UNITS.MeV)
print('electron\n', 'Samples: ', samples_electron, 'MAX_MOM:', max_mom_electron / UNITS.MeV, 'T_dec:', electron.decoupling_temperature / UNITS.MeV)
print('muon\n', 'Samples: ', samples_muon, 'MAX_MOM:', max_mom_muon / UNITS.MeV, 'T_dec:', muon.decoupling_temperature / UNITS.MeV)
print('pion\n', 'Samples: ', samples_meson, 'MAX_MOM:', max_mom_meson / UNITS.MeV, 'T_dec:', neutral_pion.decoupling_temperature / UNITS.MeV)


if T_initial > CONST.lambda_QCD:
    universe.evolve(CONST.lambda_QCD, export=False)
    universe.QCD_transition(
        hadrons=[
            charged_pion, neutral_pion, eta,
            neutral_rho, charged_rho,
            eta_prime, omega
        ],
        quarkic_interactions=sterile_quark_interactions,
        hadronic_interactions=sterile_hadrons_interactions,
        secondary_interactions=secondary_interactions_mesons
    )

universe.evolve(T_freeze, export=False)
universe.interactions = tuple()
universe.evolve(T_final)

