# -*- coding: utf-8 -*-
"""
# Neutrino Minimal Standard Model particles and interactions
"""

from __future__ import division

import numpy as np
import itertools
from collections import Counter

from common import UNITS, CONST, utils, statistics as STATISTICS
from interactions import CrossGeneratingInteraction
from interactions.three_particle import ThreeParticleM, ThreeParticleIntegral
from interactions.four_particle import FourParticleIntegral
from interactions.four_particle.cpp.integral import CollisionIntegralKind
from library.SM import WeakM, particles as SM_particles, interactions as SMI


class SterileM(WeakM):

    """ ## Sterile interactions matrix element
        Sterile processes usually include a common factor of $32 θ^2 G_F^2$ """

    def __init__(self, theta=1., *args, **kwargs):
        super(SterileM, self).__init__(*args, **kwargs)

        self.theta = theta
        self.K1 *= self.theta**2
        self.K2 *= self.theta**2

    def __str__(self):
        """ String-like representation of the matrix element """
        ret = ""
        if self.K1:
            ret += "K1=32 θ^2 G_F^2 {: .3e} ".format(self.K1 / self.const / self.theta**2)
        if self.K2:
            ret += "K2=32 θ^2 G_F^2 {: .3e} ".format(self.K2 / self.const / self.theta**2)
        return ret + self.order_format() + " θ ={: .3e}".format(self.theta)


class particles(object):
    @staticmethod
    def sterile_neutrino(mass=33.9 * UNITS.MeV, mixing_angle=1.):
        return {
            'name': 'Sterile neutrino (Majorana)',
            'symbol': 'N',
            'statistics': STATISTICS.FERMION,
            'mass': mass,
            'mixing_angle': mixing_angle,
            'dof': 2,
            'decoupling_temperature': 1. * UNITS.GeV,
            'majorana': True,
            'Q': 0
        }

    @staticmethod
    def dirac_sterile_neutrino(mass=33.9 * UNITS.MeV, mixing_angle=1.):
        return {
            'name': 'Sterile neutrino (Dirac)',
            'symbol': 'N',
            'statistics': STATISTICS.FERMION,
            'mass': mass,
            'mixing_angle': mixing_angle,
            'dof': 4,
            'decoupling_temperature': 1. * UNITS.GeV,
            'majorana': False,
            'Q': 0
        }


class interactions(object):

    @staticmethod
    def sterile_active_scattering(theta=1., sterile=None, active_a=None,
                                  active_b=None, kind=None):
        """ \begin{align}
                N_S + \nu_{\beta} &\to \nu_{\alpha} + \nu_{\beta}
            \end{align}
        """

        # SOMETHING IS CLEARLY NOT RIGHT
        K1 = 2. if active_a == active_b else 1

        return [CrossGeneratingInteraction(
            name="Sterile-active neutrino scattering",
            particles=particles,
            antiparticles=((False, False), (False, False)),
            Ms=Ms,
            integral_type=FourParticleIntegral,
            kind=kind
        ) for particles, Ms in zip(
            [((sterile, active_b), (active_a, active_b)),
             ((active_a, active_b), (sterile, active_b))],
            [(SterileM(theta=theta, K1=K1, order=(0, 1, 2, 3)), ),
             (SterileM(theta=theta, K1=K1, order=(2, 3, 0, 1)), )]
        )]

    @staticmethod
    def sterile_active_to_leptons_NC(theta=1., g_L=CONST.g_R+0.5, sterile=None,
                                     active=None, lepton=None, kind=None):
        """ \begin{align}
                N_S + \overline{\nu_{\alpha}} \to l^+ + l^-
            \end{align}
        """

        return [CrossGeneratingInteraction(
            name="Neutral current sterile-lepton interactions",
            particles=((sterile, active), (lepton, lepton)),
            antiparticles=antiparticles,
            Ms=(
                SterileM(theta=theta, K1=4. * g_L**2, order=(0, 2, 1, 3)),
                SterileM(theta=theta, K1=4. * CONST.g_R**2, order=(0, 3, 1, 2)),
                SterileM(theta=theta, K2=4. * CONST.g_R * g_L, order=(2, 3, 0, 1))
            ),
            integral_type=FourParticleIntegral,
            kind=kind
        ) for antiparticles in
            [((False, True), (True, False)),
             ((True, False), (False, True))]
        ]

    @staticmethod
    def sterile_active_to_leptons_CC(theta=1., sterile=None,
                                     active=None, lepton1=None, lepton2=None, kind=None):
        """ \begin{align}
                N_S + \overline{\nu_{\alpha}} \to l1^+ + l2^-
            \end{align}
        """

        return [CrossGeneratingInteraction(
            name="Charged current sterile-lepton interactions",
            particles=particles,
            antiparticles=((False, True), (True, False)),
            Ms=(
                SterileM(theta=theta, K1=4., order=(0, 2, 1, 3)),
            ),
            integral_type=FourParticleIntegral,
            kind=kind
        ) for particles in [
            ((sterile, active), (lepton2, lepton1)),
            ((active, sterile), (lepton1, lepton2))
        ]]

    @classmethod
    def sterile_leptons_interactions(cls, thetas=None, sterile=None,
                                     neutrinos=None, leptons=None, kind=None):

        g_R = CONST.g_R
        inters = []

        # Neutrinos scatterings
        for neutrino_a, neutrino_b in itertools.product(neutrinos, neutrinos):
            if thetas[neutrino_a.flavour]:
                inter = cls.sterile_active_scattering(
                    theta=thetas[neutrino_a.flavour],
                    sterile=sterile,
                    active_a=neutrino_a,
                    active_b=neutrino_b,
                    kind=kind
                )
                inter[1].integrals = [
                    integral for integral in inter[1].integrals
                    if integral.reaction[0].specie.name not in [
                        integral2.reaction[0].specie.name for integral2 in inter[0].integrals
                    ]
                ]

                inters.extend(inter)

        # Interactions of neutrinos and leptons
        for lepton in leptons:
            for neutrino in neutrinos:
                if thetas[neutrino.flavour]:
                    g_L = g_R + 0.5 if lepton.flavour == neutrino.flavour else g_R - 0.5
                    # Neutral current
                    inter = cls.sterile_active_to_leptons_NC(
                        theta=thetas[neutrino.flavour],
                        g_L=g_L,
                        sterile=sterile,
                        active=neutrino,
                        lepton=lepton,
                        kind=kind
                    )
                    if lepton.name in ['Muon', 'Tau']:
                        # for interac in inter:
                        #     interac.integrals = [integral for integral in interac.integrals if not (integral.reaction[0].specie.name != sterile.name and\
                        #                         Counter([item.specie.name for item in integral.reaction if item.side == -1])[lepton.name] == 1)]
                        inter[1].integrals = [integral for integral in inter[1].integrals if integral.reaction[0].specie.name != lepton.name]
                    inters.extend(inter)

                if len(leptons) > 1:
                    other_lepton = [par for par in leptons if par != lepton]
                    for other_lep in other_lepton:
                        if thetas[lepton.flavour] and neutrino.flavour == other_lep.flavour:
                            # Charged current
                            inter = cls.sterile_active_to_leptons_CC(
                                theta=thetas[lepton.flavour],
                                sterile=sterile,
                                active=neutrino,
                                lepton1=lepton,
                                lepton2=other_lep,
                                kind=kind
                            )
                            for interac in inter:
                                interac.integrals = [integral for integral in interac.integrals if
                                                     sum([item.side for item in integral.reaction if item.specie.name in [other_lep.name, sterile.name]]) == 0]
                            inters.extend(inter)

        for inter in inters:
            for integral in inter.integrals:
                if hasattr(integral.reaction[0].specie, 'fast_decay'):
                    integral.kind = CollisionIntegralKind.F_creation

        return inters

    # ## Quark interactions

    @staticmethod
    def sterile_quark_interactions(thetas=None, sterile=None,
                                   neutrinos=None, leptons=None, quarks=None, kind=None):

        inters = []

        for quark in quarks:
            for active in neutrinos:
                if thetas[active.flavour]:
                    inters.append(interactions.sterile_quark_neutral(
                        theta=thetas[active.flavour], sterile=sterile,
                        active=active, quark=quark, kind=kind
                    ))

        for up, down in itertools.product([q for q in quarks if q.Q == 2./3.],
                                          [q for q in quarks if q.Q == -1./3.]):
            for lepton in leptons:
                if thetas[lepton.flavour]:
                    inters.append(interactions.sterile_quark_charged(
                        theta=thetas[lepton.flavour], sterile=sterile,
                        lepton=lepton, up=up, down=down, kind=kind
                    ))

        return inters

    @staticmethod
    def sterile_quark_neutral(theta=1., sterile=None, active=None, quark=None, kind=None):
        """ \begin{align}
                N_S + \overline{\nu} &\to q + \overline{q}
            \end{align}
        """

        if quark.Q == 2./3.:
            g_R = CONST.g_R
        if quark.Q == -1./3.:
            g_R = CONST.g_R / 2.

        x = (3 - 4 * g_R)

        return CrossGeneratingInteraction(
            name="Sterile-quarks neutral channel interaction",
            particles=((sterile, active), (quark, quark)),
            antiparticles=((False, True), (False, True)),
            Ms=(
                SterileM(theta=theta, K1=16./9. * g_R**2, order=(0, 2, 1, 3)),
                SterileM(theta=theta, K1=1./9. * x**2, order=(0, 3, 1, 2)),
                SterileM(theta=theta, K2=-4./9. * g_R * x, order=(2, 3, 0, 1)),
            ),
            integral_type=FourParticleIntegral,
            kind=kind
        )

    @staticmethod
    def sterile_quark_charged(theta=1., sterile=None, lepton=None, up=None, down=None, kind=None):
        """ \begin{align}
                N_S + l^+ &\to u + \overline{d}
            \end{align}
        """

        CKM = SM_particles.quarks.CKM[(up.family, down.family)]

        return CrossGeneratingInteraction(
            name="Sterile-quarks charged channel interaction",
            particles=((sterile, lepton), (up, down)),
            antiparticles=((False, True), (False, True)),
            Ms=(
                SterileM(theta=theta, K1=4. * CKM**2, order=(0, 3, 1, 2)),
            ),
            integral_type=FourParticleIntegral,
            kind=kind
        )

    @classmethod
    def sterile_hadrons_interactions(cls, thetas=None, sterile=None, neutrinos=None,
                                     leptons=None, mesons=None, kind=None):

        inters = []

        for neutrino in neutrinos:
            if thetas[neutrino.flavour]:
                for particle in mesons:
                    if particle.Q == 0:
                        if particle.type == "scalar":
                            inters += cls.neutral_scalar_meson(theta=thetas[neutrino.flavour],
                                                               sterile=sterile,
                                                               active=neutrino,
                                                               meson=particle,
                                                               kind=kind)

                        if particle.type == "vector":
                            inters += cls.neutral_vector_meson(theta=thetas[neutrino.flavour],
                                                               sterile=sterile,
                                                               active=neutrino,
                                                               meson=particle,
                                                               kind=kind)

        for lepton in leptons:
            if thetas[lepton.flavour]:
                for particle in mesons:
                    if particle.Q == -1:
                        if particle.type == "scalar":
                            inters += cls.charged_scalar_meson(theta=thetas[lepton.flavour],
                                                               sterile=sterile,
                                                               lepton=lepton,
                                                               meson=particle,
                                                               kind=kind)

                        if particle.type == "vector":
                            inters += cls.charged_vector_meson(theta=thetas[neutrino.flavour],
                                                               sterile=sterile,
                                                               lepton=lepton,
                                                               meson=particle,
                                                               kind=kind)

        for inter in inters:
            for integral in inter.integrals:
                if hasattr(integral.reaction[0].specie, 'fast_decay'):
                    if utils.reaction_type(integral).CREATION:
                        integral.kind = CollisionIntegralKind.F_creation
                    else:
                        integral.kind = CollisionIntegralKind.F_decay

        return inters

    @staticmethod
    def neutral_scalar_meson(theta=1., sterile=None, active=None, meson=None, kind=None):
        """ \begin{align}
                N \to \nu + meson
            \end{align}
        """
        return [CrossGeneratingInteraction(
            name="Sterile neutrino decay to neutral scalar meson and neutrino",
            particles=((sterile, ), (active, meson)),
            antiparticles=antiparticles,
            Ms=(ThreeParticleM(K=(
                (CONST.G_F * theta * meson.decay_constant)**2
                * sterile.mass**4 * np.abs(1 - (meson.mass / sterile.mass)**2))
            ), ),
            integral_type=ThreeParticleIntegral,
            kind=kind
        ) for antiparticles in [
            ((False, ), (False, False)),
            ((True, ), (True, False))
        ]]

    @staticmethod
    def charged_scalar_meson(theta=1., sterile=None, lepton=None, meson=None, kind=None):
        """ \begin{align}
                N \to l + meson
            \end{align}
        """

        CKM = SM_particles.quarks.CKM[(1, 1)]

        return [CrossGeneratingInteraction(
            name="Sterile neutrino decay to charged scalar meson and lepton",
            particles=((sterile, ), (lepton, meson)),
            antiparticles=antiparticles,
            Ms=(ThreeParticleM(
                K=2 * (CONST.G_F * theta * meson.decay_constant * CKM)**2 * sterile.mass**4 * np.abs(
                    (1 - (lepton.mass / sterile.mass)**2)**2
                    - (meson.mass / sterile.mass)**2 * (1 + (lepton.mass / sterile.mass)**2)
                )
            ),),
            integral_type=ThreeParticleIntegral,
            kind=kind
        ) for antiparticles in [
            ((False, ), (False, True)),
            ((True, ), (True, False))
        ]]

    @staticmethod
    def neutral_vector_meson(theta=1., sterile=None, active=None, meson=None, kind=None):
        """ \begin{align}
                N \to \nu + meson
            \end{align}
        """

        if meson.name == 'Neutral rho':
            kappa = 1. - 2. * CONST.sin_theta_w_2
        if meson.name == 'Omega':
            kappa = 4. * CONST.sin_theta_w_2 / 3.
        if meson.name == 'Phi':
            kappa = 4. * CONST.sin_theta_w_2 / 3. - 1.

        return [CrossGeneratingInteraction(
            name="Sterile neutrino decay to neutral vector meson and neutrino",
            particles=((sterile, ), (active, meson)),
            antiparticles=antiparticles,
            Ms=(ThreeParticleM(K=(
                (CONST.G_F * theta * meson.decay_constant * kappa)**2
                * sterile.mass**4 * (1 + 2 * (meson.mass / sterile.mass)**2)
                * (1 - (meson.mass / sterile.mass)**2))),),
            integral_type=ThreeParticleIntegral,
            kind=kind
        ) for antiparticles in [
            ((False, ), (False, False)),
            ((True, ), (True, False))
        ]]

    @staticmethod
    def charged_vector_meson(theta=1., sterile=None, lepton=None, meson=None, kind=None):
        """ \begin{align}
                N \to l + meson
            \end{align}
        """

        CKM = SM_particles.quarks.CKM[(1, 1)]
        # Only charged rho for now
        return [CrossGeneratingInteraction(
            name="Sterile neutrino decay to charged vector meson and lepton",
            particles=((sterile, ), (lepton, meson)),
            antiparticles=antiparticles,
            Ms=(ThreeParticleM(
                K=2 * (CONST.G_F * theta * meson.decay_constant * CKM)**2 * sterile.mass**4 * (
                    (1 - (lepton.mass / sterile.mass)**2)**2 + (meson.mass / sterile.mass)**2
                    * (1 + (lepton.mass / sterile.mass)**2) - 2 * (meson.mass / sterile.mass)**4
                )
            ),),
            integral_type=ThreeParticleIntegral,
            kind=kind
        ) for antiparticles in [
            ((False, ), (False, True)),
            ((True, ), (True, False))
        ]]

    @classmethod
    def interactions_decay_products(cls, interactions_primary=None, interactions_SM=None, muon_dec=False,
                                    neutrinos=None, leptons=None, mesons=None, photon=None, kind=None):

        if kind in [CollisionIntegralKind.F_1_vacuum_decay, CollisionIntegralKind.F_f_vacuum_decay]:
            kind_creation = CollisionIntegralKind.F_1_vacuum_decay
            kind_decay = CollisionIntegralKind.F_f_vacuum_decay
        else:
            kind_creation = CollisionIntegralKind.F_creation
            kind_decay = CollisionIntegralKind.F_decay

        interactions = []
        primary_mesons = []
        species = Counter()

        for main in interactions_primary:
            for inter in main:
                for integral in inter.integrals:
                    species.update(Counter(item.specie.name for item in integral.reaction if item.side == 1))

        if species['Muon'] and not muon_dec or species['Tau']:
            interactions += SMI.lepton_interactions(
                leptons=leptons,
                neutrinos=neutrinos,
                # SM_inters=already_there(interactions_SM, ['Muon', 'Tau']),
                kind=kind
            )

        if mesons:
            for meson in mesons:
                if species[meson.name]:
                    primary_mesons.append(meson)

        if muon_dec:
            species['Muon'] += 1

        if primary_mesons:
            interactions += SMI.meson_interactions(
                primary_mesons=primary_mesons,
                mesons=mesons,
                leptons=leptons,
                neutrinos=neutrinos,
                photon=photon,
                muon_tau=[species['Muon'], species['Tau']],
                kind=kind
            )

        for inter in interactions:
            inter.integrals = [
                integral for integral in inter.integrals
                if utils.reaction_type(integral).CREATION or utils.reaction_type(integral).DECAY
            ]
            for integral in inter.integrals:
                if utils.reaction_type(integral).CREATION:
                    integral.kind = kind_creation
                else:
                    integral.kind = kind_decay

        return interactions

