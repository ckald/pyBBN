# -*- coding: utf-8 -*-
"""
# Neutrino Minimal Standard Model particles and interactions
"""

from __future__ import division

import itertools

from common import UNITS, CONST
from particles import STATISTICS
from interactions import Interaction
from interactions.three_particle import ThreeParticleM, ThreeParticleIntegral
from interactions.four_particle import FourParticleIntegral
from library.SM import WeakM, particles as SM_particles


class particles(object):
    @staticmethod
    def sterile_neutrino(mass=33.9 * UNITS.MeV):
        return {
            'name': 'Sterile neutrino',
            'symbol': 'N',
            'statistics': STATISTICS.FERMION,
            'mass': mass,
            'dof': 2,
            'decoupling_temperature': 1. * UNITS.GeV,
            'majorana': True
        }


class interactions(object):

    @staticmethod
    def sterile_active_scattering(theta=1., sterile=None, active_a=None, active_b=None):
        """ \begin{align}
                N_S + \nu_{\beta} &\to \nu_{\alpha} + \nu_{\beta}
            \end{align}
        """

        # SOMETHING IS CLEARLY NOT RIGHT
        K1 = 2. if active_a == active_b else 1

        return Interaction(
            name="Sterile-active neutrino scattering",
            particles=((sterile, active_b), (active_a, active_b)),
            antiparticles=((False, True), (False, True)),
            Ms=(WeakM(K1=theta**2 * K1, order=(0, 1, 2, 3)), ),
            integral=FourParticleIntegral
        )

    @staticmethod
    def sterile_active_to_leptons(theta=1., g_L=CONST.g_R+0.5,
                                  sterile=None, active=None, lepton=None):
        """ \begin{align}
                N_S + \overline{\nu_{\alpha}} \to l^+ + l^-
            \end{align}
        """

        return Interaction(
            name="Sterile decay into neutrino and leptons pair",
            particles=((sterile, active), (lepton, lepton)),
            antiparticles=((False, True), (True, False)),
            Ms=(
                WeakM(K1=4. * g_L**2 * theta**2, order=(0, 3, 1, 2)),
                WeakM(K1=4. * CONST.g_R**2 * theta**2, order=(0, 2, 1, 3)),
                WeakM(K2=4. * CONST.g_R * g_L * theta**2, order=(2, 3, 0, 1))
            ),
            integral=FourParticleIntegral
        )

    @classmethod
    def sterile_leptons_interactions(cls, thetas=None, sterile=None,
                                     neutrinos=None, leptons=None):

        g_R = CONST.g_R
        inters = []

        # Neutrinos scatterings
        for neutrino_a, neutrino_b in itertools.combinations_with_replacement(neutrinos, 2):
            if thetas[neutrino_a.flavour]:
                inters.append(cls.sterile_active_scattering(
                    theta=thetas[neutrino_a.flavour],
                    sterile=sterile,
                    active_a=neutrino_a,
                    active_b=neutrino_b
                ))

        # Interactions of neutrinos and leptons
        for lepton in leptons:
            for neutrino in neutrinos:
                if thetas[neutrino.flavour]:
                    g_L = g_R + 0.5 if lepton.flavour == neutrino.flavour else g_R - 0.5
                    inters.append(cls.sterile_active_to_leptons(
                        theta=thetas[neutrino.flavour],
                        g_L=g_L,
                        sterile=sterile,
                        active=neutrino,
                        lepton=lepton
                    ))

        return inters

    # ## Quark interactions

    @staticmethod
    def sterile_quark_interactions(thetas=None, sterile=None,
                                   neutrinos=None, leptons=None, quarks=None):

        inters = []

        for quark in quarks:
            for active in neutrinos:
                if thetas[active.flavour]:
                    inters.append(interactions.sterile_quark_neutral(
                        theta=thetas[active.flavour], sterile=sterile,
                        active=active, quark=quark
                    ))

        for up, down in itertools.product([q for q in quarks if q.Q == 2./3.],
                                          [q for q in quarks if q.Q == -1./3.]):
            for lepton in leptons:
                if thetas[lepton.flavour]:
                    inters.append(interactions.sterile_quark_charged(
                         theta=thetas[lepton.flavour], sterile=sterile,
                         lepton=lepton, up=up, down=down
                     ))

        return inters

    @staticmethod
    def sterile_quark_neutral(theta=1., sterile=None, active=None, quark=None):
        """ \begin{align}
                N_S + \overline{\nu} &\to q + \overline{q}
            \end{align}
        """

        g_R = CONST.g_R

        if quark.Q == 2./3.:
            x = (3 - 4 * g_R)
        if quark.Q == -1./3.:
            x = (3 - 2 * g_R)

        return Interaction(
            name="Sterile-quarks neutral channel interaction",
            particles=((sterile, active), (quark, quark)),
            antiparticles=((False, True), (False, True)),
            Ms=(
                WeakM(K1=32./9. * g_R * theta**2, order=(0, 3, 1, 2)),
                WeakM(K1=2./9. * x**2 * theta**2, order=(0, 2, 1, 3)),
                WeakM(K2=16./9. * g_R * x * theta**2, order=(0, 1, 2, 3)),
            ),
            integral=FourParticleIntegral
        )

    @staticmethod
    def sterile_quark_charged(theta=1., sterile=None, lepton=None, up=None, down=None):
        """ \begin{align}
                N_S + l^+ &\to u + \overline{d}
            \end{align}
        """

        CKM = SM_particles.quarks.CKM[(up.family, down.family)]

        return Interaction(
            name="Sterile-quarks charged channel interaction",
            particles=((sterile, lepton), (up, down)),
            antiparticles=((False, True), (False, True)),
            Ms=(
                WeakM(K1=.5 * CKM**2 * theta**2, order=(0, 2, 1, 3)),
                WeakM(K1=.5 * CKM**2 * theta**2, order=(0, 3, 1, 2)),
            ),
            integral=FourParticleIntegral
        )

    @classmethod
    def sterile_hadrons_interactions(cls, thetas=None, sterile=None, neutrinos=None,
                                     leptons=None, hadrons=None):

        inters = []

        for neutrino in neutrinos:
            if thetas[neutrino.flavour]:
                for meson in hadrons:
                    if meson.Q == 0:
                        inters.append(cls.sterile_pion_neutral(
                            theta=thetas[neutrino.flavour],
                            sterile=sterile,
                            active=neutrino,
                            pion=meson))

        for lepton in leptons:
            if thetas[lepton.flavour]:
                for meson in hadrons:
                    if meson.Q == -1:
                        inters.append(cls.sterile_pion_charged(
                            theta=thetas[lepton.flavour],
                            sterile=sterile,
                            lepton=lepton,
                            pion=meson))

        return inters

    @staticmethod
    def sterile_pion_neutral(theta=1., sterile=None, active=None, pion=None):
        """ \begin{align}
                N \to \nu + \pi^0
            \end{align}
        """

        return Interaction(
            name="Sterile neutrino decay to neutral pion and neutrino",
            particles=((sterile, ), (active, pion)),
            antiparticles=((False, ), (False, False)),
            Ms=(ThreeParticleM(K=(CONST.G_F * theta * CONST.f_pi)**2
                               * sterile.mass**2 * (sterile.mass**2 - pion.mass**2)), ),
            integral=ThreeParticleIntegral
        )

    @staticmethod
    def sterile_pion_charged(theta=1., sterile=None, lepton=None, pion=None):
        """ \begin{align}
                N \to l + \pi^+
            \end{align}
        """

        CKM = SM_particles.quarks.CKM[(1, 1)]

        return Interaction(
            name="Sterile neutrino decay to charged pion and lepton",
            particles=((sterile, ), (lepton, pion)),
            antiparticles=((False, ), (False, True)),
            Ms=(ThreeParticleM(
                K=(CONST.G_F * theta * CONST.f_pi * CKM)**2
                * (
                   (sterile.mass**2 - lepton.mass**2)**2
                   - pion.mass**2 * (sterile.mass**2 + lepton.mass**2)
                )
            ),),
            integral=ThreeParticleIntegral
        )
