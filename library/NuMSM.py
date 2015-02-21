# -*- coding: utf-8 -*-
"""
== Neutrino Minimal Standard Model particles and interactions ==
"""

from common import UNITS, CONST
from particles import STATISTICS
from interactions import Interaction
from library.SM import WeakM


class SterileWeakM(WeakM):

    """ == Weak interactions matrix element ==
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
            'decoupling_temperature': 50 * UNITS.MeV
        }


class interactions(object):

    @staticmethod
    def sterile_decay_one_flavour(theta=1., sterile=None, active=None):
        """ \begin{align}
                N_S \to \nu_{\alpha} + \nu_{\alpha} + \overline{\nu_{\alpha}}
            \end{align}
        """
        return Interaction(
            name="Sterile decay into neutrinos (1 flavour)",
            in_particles=[sterile],
            out_particles=[active, active, active],
            Ms=[SterileWeakM(K1=2., theta=theta, order=(0, 3, 1, 2))]
        )

    @staticmethod
    def sterile_decay_two_flavours(theta=1., sterile=None, active_a=None, active_b=None):
        """ \begin{align}
                N_S \to \nu_{\alpha} + \nu_{\beta} + \overline{\nu_{\beta}}
            \end{align}
        """
        return Interaction(
            name="Sterile decay into neutrinos (2 flavours)",
            in_particles=[sterile],
            out_particles=[active_a, active_b, active_b],
            Ms=[SterileWeakM(K1=1., theta=theta, order=(0, 3, 1, 2))]
        )

    @staticmethod
    def sterile_decay_to_leptons(theta=1., g_L=CONST.g_R+0.5,
                                 sterile=None, active=None, lepton=None):
        """ \begin{align}
                N_S \to \nu_{\alpha} + l^+ + l^-
            \end{align}
        """
        return Interaction(
            name="Sterile decay into neutrino and leptons pair",
            in_particles=[sterile],
            out_particles=[active, lepton, lepton],
            Ms=[
                SterileWeakM(K1=4. * g_L**2, theta=theta, order=(0, 2, 1, 3)),
                SterileWeakM(K1=4. * CONST.g_R**2, theta=theta, order=(0, 3, 1, 2)),
                SterileWeakM(K2=4. * CONST.g_R * g_L, theta=theta, order=(0, 1, 2, 3))
            ]
        )

    @staticmethod
    def sterile_active_mixing(theta=1., sterile=None, active=None):
        """ \begin{align}
                N_S + \nu_{\alpha} &\to \nu_{\alpha} + \nu_{\alpha}
                \\\\ N_S + \overline{\nu_{\alpha}} &\to \nu_{\alpha} + \overline{\nu_{\alpha}}
            \end{align}
        """
        return Interaction(
            name="Sterile-active neutrino mixing",
            in_particles=[sterile, active],
            out_particles=[active, active],
            Ms=[
                SterileWeakM(K1=2., theta=theta, order=(0, 1, 2, 3)),
                SterileWeakM(K1=4., theta=theta, order=(0, 3, 1, 2)),
            ]
        )

###### SOMETHING FISHY IS GOING ON HERE

    @staticmethod
    def sterile_flavour_change(theta=1., sterile=None, active_a=None, active_b=None):
        """ \begin{align}
                N_S + \nu_{\alpha} &\to \nu_{\beta} + \overline{\nu_{\beta}}
                \\\\ N_S + \nu_{\alpha} &\to \nu_{\alpha} + \nu_{\alpha}
            \end{align}
        """

        K1 = 1. if active_a != active_b else 4.

        return Interaction(
            name="Sterile-active neutrino mixing",
            in_particles=[sterile, active_a],
            out_particles=[active_b, active_b],
            Ms=[
                SterileWeakM(K1=K1, theta=theta, order=(0, 3, 1, 2)),
            ]
        )

    @staticmethod
    def sterile_scattering(theta=1., sterile=None, active_a=None, active_b=None):
        """ \begin{align}
                N_S + \nu_{\beta} &\to \nu_{\alpha} + \nu_{\beta}
                \\\\ N_S + \overline{\nu_{\beta}} &\to \nu_{\alpha} + \overline{\nu_{\beta}}
            \end{align}
        """

        return Interaction(
            name="Sterile-active neutrino mixing",
            in_particles=[sterile, active_a],
            out_particles=[active_b, active_b],
            Ms=[
                SterileWeakM(K1=1., theta=theta, order=(0, 1, 2, 3)),
                SterileWeakM(K1=1., theta=theta, order=(0, 3, 1, 2)),
            ]
        )
