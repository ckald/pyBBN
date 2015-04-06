# -*- coding: utf-8 -*-
"""
# Neutrino Minimal Standard Model particles and interactions
"""

from __future__ import division

import itertools

from common import UNITS, CONST
from particles import STATISTICS
from interactions import Interaction
from library.SM import WeakM, particles as SM_particles


class SterileWeakM(WeakM):

    """ ## Weak interactions matrix element
        Weak processes usually include a common factor of $32 G_F^2$ """

    def __init__(self, *args, **kwargs):
        super(SterileWeakM, self).__init__(*args, **kwargs)

        self.K1 *= kwargs['theta']**2
        self.K2 *= kwargs['theta']**2


class particles(object):
    @staticmethod
    def sterile_neutrino(mass=33.9 * UNITS.MeV):
        return {
            'name': 'Sterile neutrino',
            'symbol': 'N',
            'statistics': STATISTICS.FERMION,
            'mass': mass,
            'dof': 2,
            'decoupling_temperature': 50 * UNITS.MeV,
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
            Ms=(SterileWeakM(K1=K1, theta=theta, order=(0, 1, 2, 3)), )
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
                SterileWeakM(K1=4. * g_L**2, theta=theta, order=(0, 3, 1, 2)),
                SterileWeakM(K1=4. * CONST.g_R**2, theta=theta, order=(0, 2, 1, 3)),
                SterileWeakM(K2=4. * CONST.g_R * g_L, theta=theta, order=(2, 3, 0, 1))
            )
        )

    @classmethod
    def sterile_leptons_interactions(cls, thetas=None, sterile=None,
                                     neutrinos=None, leptons=None):

        g_R = CONST.g_R
        inters = []

        # Neutrinos scatterings
        for neutrino_a, neutrino_b in itertools.combinations_with_replacement(neutrinos, 2):
            if thetas[neutrino_b.flavour]:
                inters.append(cls.sterile_active_scattering(
                    theta=thetas[neutrino_b.flavour],
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
            print (up, down)
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
                SterileWeakM(K1=32./9. * g_R, theta=theta, order=(0, 3, 1, 2)),
                SterileWeakM(K1=2./9. * x**2, theta=theta, order=(0, 2, 1, 3)),
                SterileWeakM(K2=16./9. * g_R * x, theta=theta, order=(0, 1, 2, 3)),
            )
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
                SterileWeakM(K1=1./2. * CKM**2, theta=theta, order=(0, 2, 1, 3)),
                SterileWeakM(K1=1./2. * CKM**2, theta=theta, order=(0, 3, 1, 2)),
            )
        )
