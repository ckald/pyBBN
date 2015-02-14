# -*- coding: utf-8 -*-
"""
== Neutrino Minimal Standard Model particles and interactions ==
"""

from common import UNITS
from particles import STATISTICS
from interactions import Interaction
from interactions.four_particle import WeakM


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
    def sterile_active_mixing(sterile=None, active=None):
        """ \begin{align}
                N_S + \nu_{\alpha} &\to \nu_{\alpha} + \nu_{\alpha}
                \\\\ N_S + \overline{\nu_{\alpha}} &\to \nu_{\alpha} + \overline{\nu_{\alpha}}
            \end{align}
        """
        return Interaction(
            name="Sterile-active neutrino mixing",
            in_particles=[sterile, active],
            out_particles=[active, active]
        )
